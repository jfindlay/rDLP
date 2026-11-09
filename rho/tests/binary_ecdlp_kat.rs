//! Known-answer tests for the binary-curve Pollard-rho ECDLP solver.
//!
//! # What is tested
//!
//! 1. **End-to-end `k·G = Q` recovery** — the solver returns a [`DlogSolution`] whose `k`
//!    satisfies `k·G = Q`, and whose evidence identity `a₁·G + b₁·Q = a₂·G + b₂·Q (mod n)`
//!    rechecks from the fields alone.
//! 2. **Walk-state invariant** — `W = a·G + b·Q` holds after construction and after a sequence of
//!    steps (the walk-state correctness contract).
//! 3. **Input checks and failure classes** — invalid instances are rejected before any walk, a
//!    degenerate collision consumes an attempt, and an exhausted budget is reported as such.
//! 4. **τ-orbit rho** — `solve_koblitz` recovers a valid `k`, the τ eigenvalues satisfy
//!    `τ^j(P) = λ_j·P`, and a corrupted eigenvalue ends in `Inconsistent`, never a wrong `Ok`.
//!
//! # Toy curve
//!
//! All tests use the curve `y²+xy = x³+x²+1` over GF(2^4) with `x⁴+x+1` (poly = 0x13), group
//! order n = 4:
//! - G  = (1, 6)
//! - 2G = (0, 1)
//! - 3G = (1, 7) = −G
//! - 4G = ∞
//!
//! The group order n = 4 is composite (2²), so the input check does not require a prime order and
//! a collision is degenerate when `gcd(b₂−b₁, 4) > 1`; such attempts are retried.
//!
//! On this curve every orbit of τ together with negation has a single x-coordinate, so the
//! canonical-representative tie-break always keeps the `k = 0` member; of the eigenvalue table only
//! `λ₀` is consulted by the walk, which is why the corrupted-eigenvalue test corrupts `λ₀`.

use crypto_bigint::Uint;
use rho::binary_curve::{BinaryAffinePoint, BinaryCurve};
use rho::binary_ecdlp::koblitz::{
    compute_tau_eigenvalues, solve_koblitz, solve_koblitz_with_eigenvalues, tau_pow,
};
use rho::binary_ecdlp::{
    check_binary_order_input, solve, solve_brent, verify_binary_log, BinaryAddendTable,
    BinaryWalkState,
};
use rho::ecdlp::{DlogError, DlogSolution, InputDefect, RhoEvidence};
use shared_gf2m::{F2m, F2mNaive};

type F = F2mNaive<1>;
type Pt = BinaryAffinePoint<F>;

/// Group order of the toy curve.
const N: u64 = 4;
/// Field degree of GF(2^4).
const M: usize = 4;

// ── Curve and field parameters ────────────────────────────────────────────────

/// GF(2^4) irreducible: x⁴+x+1 = 0x13.
fn poly4() -> Uint<1> {
    Uint::<1>::from(0x13u64)
}

/// Toy binary curve: y²+xy = x³+x²+1 over GF(2^4) with x⁴+x+1.
///
/// Group order 4: G=(1,6), 2G=(0,1), 3G=(1,7), 4G=∞.
fn toy_curve() -> BinaryCurve {
    BinaryCurve {
        poly: poly4(),
        a: Uint::<1>::ONE,
        b: Uint::<1>::ONE,
        n: Uint::<1>::from(N),
        gx: Uint::<1>::ONE,
        gy: Uint::<1>::from(6u64),
    }
}

/// `k·P` on the toy curve.
fn mul(curve: &BinaryCurve, p: &Pt, k: u64) -> Pt {
    curve.scalar_mul(p, &Uint::<1>::from(k))
}

/// `(G, k·G)` on the toy curve.
fn instance(curve: &BinaryCurve, k_target: u64) -> (Pt, Pt) {
    let g = curve.generator::<F>();
    let q = mul(curve, &g, k_target);
    (g, q)
}

/// A point with `x = 1` and `y = 0`, which does not satisfy the curve equation.
fn off_curve_point(curve: &BinaryCurve) -> Pt {
    let p = BinaryAffinePoint::new(F::from_u64(1, &poly4()), F::from_u64(0, &poly4()));
    assert!(!curve.is_on_curve(&p), "fixture point must be off the curve");
    p
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Assert that `sol` is a verified solution and that its evidence rechecks from the fields.
///
/// Checks `k < n` and `k·G = Q`, then for a collision the identity
/// `a₁·G + b₁·Q = (−1)^negated · (a₂·G + b₂·Q)` and that `k` solves the congruence the collision
/// implies, `(b₂ − b₁)·k ≡ a₁ − a₂ (mod n)`.
fn assert_rechecked(curve: &BinaryCurve, g: &Pt, q: &Pt, sol: &DlogSolution<RhoEvidence>) {
    assert!(sol.k < N, "k = {} is not reduced mod n", sol.k);
    assert_eq!(mul(curve, g, sol.k), *q, "k = {} gives k·G ≠ Q", sol.k);
    assert!(verify_binary_log(curve, g, q, sol.k));

    match &sol.evidence {
        RhoEvidence::Trivial => {
            assert!(q.is_infinity() || q == g, "Trivial evidence for a non-trivial target");
        }
        RhoEvidence::Collision { steps, a1, b1, a2, b2, negated, .. } => {
            assert!(*steps >= 1, "a collision needs at least one walk step");
            for scalar in [a1, b1, a2, b2] {
                assert!(*scalar < N, "tracked scalar {scalar} is not reduced mod n");
            }
            let first = curve.add(&mul(curve, g, *a1), &mul(curve, q, *b1));
            let second = curve.add(&mul(curve, g, *a2), &mul(curve, q, *b2));
            let second = if *negated { curve.negate(&second) } else { second };
            assert_eq!(first, second, "evidence identity a₁G+b₁Q = ±(a₂G+b₂Q) fails");

            let n = i128::from(N);
            let lhs = (i128::from(*b2) - i128::from(*b1)) * i128::from(sol.k);
            let rhs = i128::from(*a1) - i128::from(*a2);
            assert_eq!((lhs - rhs).rem_euclid(n), 0, "k does not solve the collision congruence");
        }
    }
}

/// Run `solve` on `Q = k_target·G` and recheck the result.
fn check_solver_on_binary_curve(k_target: u64) -> DlogSolution<RhoEvidence> {
    let curve = toy_curve();
    let (g, q) = instance(&curve, k_target);
    let sol = solve(&curve, &g, &q, N)
        .unwrap_or_else(|e| panic!("solve failed for k_target={k_target}: {e:?}"));
    assert_rechecked(&curve, &g, &q, &sol);
    sol
}

/// The `InputDefect` inside an `InvalidInput` error, panicking on any other result.
fn input_defect<T: std::fmt::Debug>(r: Result<T, DlogError>) -> InputDefect {
    match r {
        Err(DlogError::InvalidInput(d)) => d,
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

// ── End-to-end ECDLP KATs ─────────────────────────────────────────────────────

/// DLP k=1 on the toy binary curve: Q = G, answered without a walk.
#[test]
fn binary_dlog_k1() {
    let sol = check_solver_on_binary_curve(1);
    assert_eq!(sol.k, 1);
    assert_eq!(sol.evidence, RhoEvidence::Trivial);
}

/// DLP k=2 on the toy binary curve: Q = 2G = (0, 1), recovered from a walk collision.
#[test]
fn binary_dlog_k2() {
    let sol = check_solver_on_binary_curve(2);
    assert!(matches!(sol.evidence, RhoEvidence::Collision { negated: false, .. }));
}

/// DLP k=3 on the toy binary curve: Q = 3G = (1, 7) = −G.
#[test]
fn binary_dlog_k3() {
    let sol = check_solver_on_binary_curve(3);
    assert!(matches!(sol.evidence, RhoEvidence::Collision { negated: false, .. }));
}

/// DLP k=0 (Q = ∞): the identity case is `Trivial` with `k = 0`, and `k·G = ∞` is verified.
#[test]
fn binary_dlog_k0_identity() {
    let curve = toy_curve();
    let g = curve.generator::<F>();
    let q: Pt = BinaryAffinePoint::Infinity;

    let sol = solve(&curve, &g, &q, N).expect("solve failed for Q=∞");
    assert_eq!(sol.k, 0);
    assert_eq!(sol.evidence, RhoEvidence::Trivial);
    assert!(mul(&curve, &g, sol.k).is_infinity());
}

// ── Walk-state invariant KATs ─────────────────────────────────────────────────

/// Walk-state invariant: `W = a·G + b·Q` holds after construction.
#[test]
fn walk_invariant_after_construction() {
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);

    let mut rng = ChaCha20Rng::seed_from_u64(0xABCD_1234);
    let walk = BinaryWalkState::<F>::new_random(&curve, &g, &q, N, &mut rng);

    let reconstructed = curve.add(&mul(&curve, &g, walk.a), &mul(&curve, &q, walk.b));
    assert_eq!(
        walk.point, reconstructed,
        "walk invariant broken at construction: point ≠ a·G + b·Q"
    );
}

/// Walk-state invariant: `W = a·G + b·Q` holds after multiple steps.
#[test]
fn walk_invariant_across_steps() {
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    let curve = toy_curve();
    let (g, q) = instance(&curve, 3);

    let mut rng = ChaCha20Rng::seed_from_u64(0xDEAD_BEEF);
    let table = BinaryAddendTable::new(&curve, &g, &q, N, &mut rng);
    let mut walk = BinaryWalkState::<F>::new_random(&curve, &g, &q, N, &mut rng);

    for step in 0..10 {
        let reconstructed = curve.add(&mul(&curve, &g, walk.a), &mul(&curve, &q, walk.b));
        assert_eq!(
            walk.point, reconstructed,
            "walk invariant broken at step {step}: point ≠ a·G + b·Q"
        );
        walk.step(&curve, &table, N);
    }
}

// ── solve_brent KATs ──────────────────────────────────────────────────────────

/// solve_brent with a fixed seed is deterministic, and its evidence rechecks.
#[test]
fn solve_brent_deterministic() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);

    let first = solve_brent::<F>(&curve, &g, &q, N, 42, 50).expect("first call failed");
    let second = solve_brent::<F>(&curve, &g, &q, N, 42, 50).expect("second call failed");

    assert_eq!(first, second, "solve_brent: different results for same seed");
    assert_rechecked(&curve, &g, &q, &first);
}

/// solve_brent succeeds across several seeds and rechecks each one.
#[test]
fn solve_brent_multiple_seeds() {
    let curve = toy_curve();
    for k_target in [2u64, 3] {
        let (g, q) = instance(&curve, k_target);
        for seed in [0u64, 1, 42, 0xDEAD_BEEF, 0xCAFE_BABE] {
            let sol = solve_brent::<F>(&curve, &g, &q, N, seed, 50)
                .unwrap_or_else(|e| panic!("solve_brent failed for seed={seed}: {e:?}"));
            assert_rechecked(&curve, &g, &q, &sol);
        }
    }
}

/// A degenerate collision (`gcd(b₂−b₁, 4) > 1`) consumes one attempt and the next attempt runs
/// with seed `seed + 1`.
///
/// Seed 0 is degenerate on `Q = 2G`: one attempt exhausts the budget, and with more attempts the
/// accepted collision is the one attempt 1 finds, which a direct run at seed 1 reproduces.
#[test]
fn degenerate_collision_consumes_an_attempt() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);

    let one = solve_brent::<F>(&curve, &g, &q, N, 0, 1);
    assert!(
        matches!(one, Err(DlogError::BudgetExhausted { attempts: 1, steps }) if steps > 0),
        "{one:?}"
    );

    let retried = solve_brent::<F>(&curve, &g, &q, N, 0, 50).expect("retries must succeed");
    let RhoEvidence::Collision { attempt, a1, b1, a2, b2, .. } = retried.evidence else {
        panic!("expected a collision");
    };
    assert_eq!(attempt, 1, "seed 0 is degenerate, so attempt 1 is the first usable one");

    let direct = solve_brent::<F>(&curve, &g, &q, N, 1, 1).expect("seed 1 succeeds first try");
    let RhoEvidence::Collision { attempt: 0, a1: da1, b1: db1, a2: da2, b2: db2, .. } =
        direct.evidence
    else {
        panic!("expected a first-attempt collision");
    };
    assert_eq!((a1, b1, a2, b2), (da1, db1, da2, db2), "attempt i must use seed + i");
}

/// With no attempts allowed the budget is exhausted before any step.
#[test]
fn zero_retries_exhaust_the_budget_immediately() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);

    assert_eq!(
        solve_brent::<F>(&curve, &g, &q, N, 0, 0),
        Err(DlogError::BudgetExhausted { attempts: 0, steps: 0 })
    );
    assert_eq!(
        solve_koblitz::<F>(&curve, &g, &q, N, M, 0, 0),
        Err(DlogError::BudgetExhausted { attempts: 0, steps: 0 })
    );
}

// ── Input checks ──────────────────────────────────────────────────────────────

/// The input check rejects each defective instance with its own defect, in order, before a walk.
#[test]
fn input_check_rejects_defective_instances() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);
    let inf: Pt = BinaryAffinePoint::Infinity;
    let bad = off_curve_point(&curve);

    assert_eq!(check_binary_order_input(&curve, &g, &q, N), Ok(()));

    let r = check_binary_order_input(&curve, &inf, &q, N);
    assert_eq!(input_defect(r), InputDefect::GeneratorAtInfinity);
    let r = check_binary_order_input(&curve, &bad, &q, N);
    assert_eq!(input_defect(r), InputDefect::GeneratorNotOnCurve);
    let r = check_binary_order_input(&curve, &g, &bad, N);
    assert_eq!(input_defect(r), InputDefect::TargetNotOnCurve);
    let r = check_binary_order_input(&curve, &g, &q, 0);
    assert_eq!(input_defect(r), InputDefect::ZeroOrder);
    // None of these annihilates G (order 4); primality of n plays no part in the check.
    for n in [1u64, 3, 5] {
        let r = check_binary_order_input(&curve, &g, &q, n);
        assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder, "n = {n}");
    }
}

/// A composite multiple of the order is accepted: only `n·G = ∞` and `n·Q = ∞` are claimed.
#[test]
fn input_check_accepts_composite_orders() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 3);
    for n in [4u64, 8, 12] {
        assert_eq!(check_binary_order_input(&curve, &g, &q, n), Ok(()), "n = {n}");
    }
}

/// With `G'` of order 2 (`G' = 2G`) and a target `Q = G` of order 4, `2·Q ≠ ∞`.
#[test]
fn input_check_rejects_target_outside_subgroup() {
    let curve = toy_curve();
    let full = curve.generator::<F>();
    let g2 = mul(&curve, &full, 2);
    assert!(mul(&curve, &g2, 2).is_infinity() && !g2.is_infinity());

    let r = check_binary_order_input(&curve, &g2, &full, 2);
    assert_eq!(input_defect(r), InputDefect::TargetOutsideSubgroup);
    // The order-2 subgroup itself is a valid instance.
    assert_eq!(check_binary_order_input(&curve, &g2, &g2, 2), Ok(()));
}

/// Every solver rejects a defective instance with `InvalidInput` instead of walking or panicking.
#[test]
fn solvers_reject_invalid_input() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);
    let inf: Pt = BinaryAffinePoint::Infinity;
    let bad = off_curve_point(&curve);

    assert_eq!(
        input_defect(solve(&curve, &g, &q, 0)),
        InputDefect::ZeroOrder,
        "n = 0 must not panic"
    );
    assert_eq!(input_defect(solve(&curve, &inf, &q, N)), InputDefect::GeneratorAtInfinity);
    assert_eq!(
        input_defect(solve_brent::<F>(&curve, &g, &bad, N, 0, 5)),
        InputDefect::TargetNotOnCurve
    );
    assert_eq!(input_defect(solve(&curve, &g, &q, 3)), InputDefect::WrongGeneratorOrder);

    assert_eq!(
        input_defect(solve_koblitz::<F>(&curve, &g, &q, 0, M, 0, 5)),
        InputDefect::ZeroOrder
    );
    assert_eq!(
        input_defect(solve_koblitz::<F>(&curve, &inf, &q, N, M, 0, 5)),
        InputDefect::GeneratorAtInfinity
    );
    assert_eq!(
        input_defect(solve_koblitz::<F>(&curve, &g, &bad, N, M, 0, 5)),
        InputDefect::TargetNotOnCurve
    );
    assert_eq!(
        input_defect(solve_koblitz::<F>(&curve, &g, &q, 3, M, 0, 5)),
        InputDefect::WrongGeneratorOrder
    );
}

// ── τ-orbit rho KATs ──────────────────────────────────────────────────────────
//
// Toy curve parameters: GF(2^4) with x⁴+x+1, group order n = 4, field degree m = 4.
//
// The τ-orbit walk folds τ-powers and signs into the tracked scalars before it compares points,
// so its evidence is an ordinary collision identity with `negated: false`.

/// Run `solve_koblitz` on `Q = k_target·G` and recheck the result.
fn check_koblitz_on_binary_curve(k_target: u64, seed: u64) -> DlogSolution<RhoEvidence> {
    let curve = toy_curve();
    let (g, q) = instance(&curve, k_target);
    let sol = solve_koblitz::<F>(&curve, &g, &q, N, M, seed, 50)
        .unwrap_or_else(|e| panic!("solve_koblitz failed for k_target={k_target}: {e:?}"));
    assert_rechecked(&curve, &g, &q, &sol);
    sol
}

/// τ-orbit rho on Q = G is answered without a walk.
#[test]
fn koblitz_dlog_k1() {
    let sol = check_koblitz_on_binary_curve(1, 0);
    assert_eq!(sol.k, 1);
    assert_eq!(sol.evidence, RhoEvidence::Trivial);
}

/// τ-orbit rho recovers k for Q = 2G; the evidence carries no sign.
#[test]
fn koblitz_dlog_k2() {
    let sol = check_koblitz_on_binary_curve(2, 0);
    assert!(matches!(sol.evidence, RhoEvidence::Collision { negated: false, .. }));
}

/// τ-orbit rho recovers k for Q = 3G = −G; the sign is folded into the scalars.
#[test]
fn koblitz_dlog_k3() {
    let sol = check_koblitz_on_binary_curve(3, 0);
    assert!(matches!(sol.evidence, RhoEvidence::Collision { negated: false, .. }));
}

/// τ-orbit rho handles Q = ∞ (k = 0).
#[test]
fn koblitz_dlog_k0_identity() {
    let curve = toy_curve();
    let g = curve.generator::<F>();
    let q: Pt = BinaryAffinePoint::Infinity;

    let sol = solve_koblitz::<F>(&curve, &g, &q, N, M, 0, 50).expect("solve_koblitz failed");
    assert_eq!(sol.k, 0);
    assert_eq!(sol.evidence, RhoEvidence::Trivial);
    assert!(mul(&curve, &g, sol.k).is_infinity());
}

/// The τ-orbit variant and the baseline both return a verified `k` for the same targets.
#[test]
fn koblitz_agrees_with_baseline() {
    let curve = toy_curve();
    for k_target in [1u64, 2, 3] {
        let (g, q) = instance(&curve, k_target);

        let baseline = solve(&curve, &g, &q, N)
            .unwrap_or_else(|e| panic!("baseline failed for k_target={k_target}: {e:?}"));
        let koblitz = solve_koblitz::<F>(&curve, &g, &q, N, M, 0, 50)
            .unwrap_or_else(|e| panic!("koblitz failed for k_target={k_target}: {e:?}"));

        assert_rechecked(&curve, &g, &q, &baseline);
        assert_rechecked(&curve, &g, &q, &koblitz);
        // k < n = ord(G), so a verified k is exactly the planted one.
        assert_eq!(baseline.k, k_target);
        assert_eq!(koblitz.k, k_target);
    }
}

/// τ-orbit rho is deterministic: same seed gives the same solution.
#[test]
fn koblitz_deterministic() {
    let first = check_koblitz_on_binary_curve(2, 42);
    let second = check_koblitz_on_binary_curve(2, 42);
    assert_eq!(first, second, "solve_koblitz: different results for same seed");
}

/// The eigenvalues from `compute_tau_eigenvalues` satisfy `τ^j(P) = λ_j·P` for `P` in `⟨G⟩`.
///
/// On the toy curve `τ(G) = 3G = −G`, so the table is `[1, 3, 1, 3]`.
#[test]
fn tau_eigenvalues_satisfy_tau_p_equals_lambda_p() {
    let curve = toy_curve();
    let g = curve.generator::<F>();
    let lambdas = compute_tau_eigenvalues(&curve, &g, M, N).expect("toy eigenvalues");
    assert_eq!(lambdas, [1, 3, 1, 3]);

    for k_target in 0..N {
        let p = mul(&curve, &g, k_target);
        for (j, lambda) in lambdas.iter().enumerate() {
            assert_eq!(
                tau_pow(&p, j, &curve.poly),
                mul(&curve, &p, *lambda),
                "τ^{j}(P) ≠ λ_{j}·P for P = {k_target}·G"
            );
        }
    }
}

/// A corrupted eigenvalue gives `Inconsistent`, never a wrong `Ok`.
///
/// `λ₀` (the only entry the toy walk consults) is replaced by 3, so folded scalars no longer
/// describe the canonical point. The same seed with the true table succeeds, and every seed ends
/// either in a verified `Ok` or in `Inconsistent` whose candidate really is wrong.
#[test]
fn corrupted_eigenvalue_is_inconsistent_not_wrong() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 3);
    let good = compute_tau_eigenvalues(&curve, &g, M, N).expect("toy eigenvalues");
    let mut corrupt = good.clone();
    corrupt[0] = 3;

    // Seed 3 hits a non-degenerate collision whose scalars are wrong under the corrupted table.
    let seed = 3;
    let sol = solve_koblitz_with_eigenvalues::<F>(&curve, &g, &q, N, M, &good, seed, 50)
        .expect("true eigenvalues solve");
    assert_rechecked(&curve, &g, &q, &sol);

    match solve_koblitz_with_eigenvalues::<F>(&curve, &g, &q, N, M, &corrupt, seed, 50) {
        Err(DlogError::Inconsistent { k }) => {
            assert_ne!(mul(&curve, &g, k), q, "Inconsistent must carry a candidate that fails");
        }
        other => panic!("expected Inconsistent for a corrupted eigenvalue, got {other:?}"),
    }

    let mut inconsistent = 0;
    for seed in 0..40 {
        match solve_koblitz_with_eigenvalues::<F>(&curve, &g, &q, N, M, &corrupt, seed, 50) {
            Ok(sol) => assert_eq!(mul(&curve, &g, sol.k), q, "a returned k must verify"),
            Err(DlogError::Inconsistent { .. }) => inconsistent += 1,
            Err(other) => panic!("unexpected failure class {other:?} at seed {seed}"),
        }
    }
    assert!(inconsistent > 0, "no seed surfaced the corrupted eigenvalue");
}

/// An eigenvalue table shorter than the field degree is refused instead of indexed out of range.
#[test]
fn short_eigenvalue_table_is_unsupported() {
    let curve = toy_curve();
    let (g, q) = instance(&curve, 2);
    let r = solve_koblitz_with_eigenvalues::<F>(&curve, &g, &q, N, M, &[1, 3], 0, 5);
    assert!(matches!(r, Err(DlogError::Unsupported(_))), "{r:?}");
}
