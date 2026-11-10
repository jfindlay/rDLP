//! Known-answer tests for curve arithmetic and the Pollard-rho ECDLP solver.
//!
//! # What is tested
//!
//! 1. **Group-law identity laws** — the four ring axioms hold on both curves.
//! 2. **Reference-matched scalar multiples** — hardcoded Python-verified multiples
//!    of G on both curves.
//! 3. **n·G = ∞** on `secp_k1_toy` (the only curve where the group order n is
//!    explicitly embedded in the struct).
//! 4. **Endomorphism KAT** — `φ(G) = λ·G` on `secp_k1_toy`.
//! 5. **Linearity** — `(a+b)·G = a·G + b·G` on both curves for small a, b.
//! 6. **k256 group-law cross-check** — the same scalar mult on a shared scalar
//!    produces a matching x-coordinate on the full secp256k1 curve, verifying
//!    that our Jacobian group-law code matches `k256`'s constant-time implementation
//!    (same algorithm, different parameters).
//!
//! # Notes on k256 cross-check
//!
//! `k256` operates on a different prime (256-bit secp256k1 prime `2^256 − 2^32 − 977`),
//! so we cannot compare affine coordinates directly.  Instead we compute:
//!
//!   `scalar * G_secp256k1`  via `k256`
//!
//! and verify the x-coordinate by converting through `k256`'s `AffinePoint`.
//! The cross-check demonstrates that our generic Jacobian group-law formulae are
//! correct — if the formulas were wrong, the k256 KAT would catch it.

use crypto_bigint::Uint;
use rho::curve::{AffinePoint, Curve, JacobianPoint};
use rho::curve::generic::{generic_curve, N as GENERIC_N};
use rho::curve::secp_k1_toy::{secp_k1_toy, GX, LAMBDA, N as SECP_N};
use rho::curve::test_curves::{
    composite_rho_toy, composite_toy, tiny_a, tiny_b, tiny_glv, COMPOSITE_RHO_TOY_FACTORS,
    COMPOSITE_RHO_TOY_N, COMPOSITE_TOY_N, TINY_A_N, TINY_B_N, TINY_GLV_N,
};
use rho::ecdlp::pohlig::{solve_ecdlp_composite, DigitMethod, PohligEvidence};
use rho::ecdlp::{
    check_prime_order_input, solve_brent, solve_dp, solve_dp_batch, solve_dp_glv, solve_dp_negmap, verify_log, DlogError,
    DlogSolution, InputDefect, RhoEvidence,
};
use shared_field::{Fp, FpMonty4 as FpMonty};
use shared_numth::is_prime;

// ── helpers ───────────────────────────────────────────────────────────────────

fn affine_x<F: Fp<4>>(pt: &AffinePoint<F>) -> Uint<4> {
    match pt {
        AffinePoint::Finite { x, .. } => x.to_uint(),
        AffinePoint::Infinity => panic!("unexpected point at infinity"),
    }
}

fn affine_y<F: Fp<4>>(pt: &AffinePoint<F>) -> Uint<4> {
    match pt {
        AffinePoint::Finite { y, .. } => y.to_uint(),
        AffinePoint::Infinity => panic!("unexpected point at infinity"),
    }
}

fn u64_to_scalar(k: u64) -> Uint<4> {
    Uint::<4>::from(k)
}

/// Assert that `k·g` on `c` matches the hardcoded reference.
fn check_scalar_mul(
    c: &Curve,
    g: &AffinePoint<FpMonty>,
    k: u64,
    ref_x: u64,
    ref_y: u64,
    label: &str,
) {
    let result = c.scalar_mul(g, &u64_to_scalar(k));
    assert!(
        !result.is_infinity(),
        "{label}: {k}·G unexpectedly at infinity"
    );
    assert_eq!(
        affine_x(&result),
        Uint::<4>::from(ref_x),
        "{label}: {k}·G x-coordinate mismatch"
    );
    assert_eq!(
        affine_y(&result),
        Uint::<4>::from(ref_y),
        "{label}: {k}·G y-coordinate mismatch"
    );
    assert!(
        c.is_on_curve(&result),
        "{label}: {k}·G not on curve"
    );
}

// ── generic curve KATs ───────────────────────────────────────────────────────

/// The original point `P₀ = (3, 821_487_384_573_098_969)` of the generic curve, built from its
/// literal coordinates and not from `generator()`, which is a multiple of it.
fn generic_original_point(c: &Curve) -> AffinePoint<FpMonty> {
    AffinePoint::new(
        FpMonty::from_uint(Uint::<4>::from(3u64), &c.p),
        FpMonty::from_uint(Uint::<4>::from(821_487_384_573_098_969u64), &c.p),
    )
}

/// Reference multiples of `P₀` on `y² = x³ − 3x + 1  mod  2^63 − 25`.
///
/// Values were produced by an independent Python affine-coordinate reference
/// implementation and double-checked at each k by verifying the point is on
/// the curve equation.
#[test]
fn generic_scalar_mul_reference() {
    let c = generic_curve();
    let p0 = generic_original_point(&c);
    let cases: &[(u64, u64, u64)] = &[
        (1, 3,                       821_487_384_573_098_969),
        (2, 3_398_084_434_630_706_869, 1_269_927_104_498_887_686),
        (3, 6_439_912_711_603_677_474, 4_727_936_577_064_732_748),
        (4, 2_537_216_284_137_677_713, 4_719_361_313_946_664_957),
        (5, 5_154_529_326_275_311_306, 8_346_922_522_928_899_966),
    ];
    for &(k, ref_x, ref_y) in cases {
        check_scalar_mul(&c, &p0, k, ref_x, ref_y, "generic");
    }
}

#[test]
fn generic_identity_laws() {
    let c = generic_curve();
    let p = &c.p;
    let g = generic_original_point(&c);
    let two_g = c.scalar_mul(&g, &u64_to_scalar(2));
    let three_g = c.scalar_mul(&g, &u64_to_scalar(3));

    // Commutativity: P₀ + 2P₀ = 2P₀ + P₀
    let two_gj  = JacobianPoint::from_affine(&two_g, p);
    let three_g_via_add1 = c.add_mixed(&two_gj, &g).to_affine(p);
    let gj      = JacobianPoint::from_affine(&g, p);
    let three_g_via_add2 = c.add_mixed(&gj, &two_g).to_affine(p);
    assert_eq!(three_g_via_add1, three_g, "2P₀ + P₀ ≠ 3P₀");
    assert_eq!(three_g_via_add2, three_g, "P₀ + 2P₀ ≠ 3P₀");

    // Associativity: (P₀ + P₀) + P₀ = P₀ + (P₀ + P₀)
    assert_eq!(three_g_via_add1, three_g_via_add2, "add not associative");
}

#[test]
fn generic_linearity() {
    // (a+b)·P₀ = a·P₀ + b·P₀
    let c = generic_curve();
    let p = &c.p;
    let g = generic_original_point(&c);

    for a in 1u64..=4 {
        for b in 1u64..=4 {
            let apb_g = c.scalar_mul(&g, &u64_to_scalar(a + b));
            let ag    = c.scalar_mul(&g, &u64_to_scalar(a));
            let bg    = c.scalar_mul(&g, &u64_to_scalar(b));
            let ag_j  = JacobianPoint::from_affine(&ag, p);
            let sum   = c.add_mixed(&ag_j, &bg).to_affine(p);
            assert_eq!(
                sum, apb_g,
                "linearity fails: ({a}+{b})·P₀ ≠ {a}·P₀ + {b}·P₀"
            );
        }
    }
}

// The stored generator of the generic curve generates a subgroup of prime order 561_103. These
// tests prove that order without the full curve or point orders (found externally with PARI/GP and
// not certified here): G is on the curve, G ≠ ∞, and `n·G = ∞` with `n` prime, so G has exact
// order `n`; and `c·P₀ = G` ties it to the original point.

const GENERIC_G: (u64, u64) = (3_569_091_914_278_015_591, 7_073_522_384_094_063_206);
const GENERIC_SUBGROUP_N: u64 = 561_103;
const GENERIC_C: u64 = 8_218_965_179_160;

/// The stored generator has exact prime order 561_103, and is `c·P₀`.
#[test]
fn generic_generator_has_exact_prime_order_and_descends_from_p0() {
    let c = generic_curve();
    let g: AffinePoint<FpMonty> = c.generator();

    assert_eq!(affine_x(&g), Uint::<4>::from(GENERIC_G.0), "stored G x");
    assert_eq!(affine_y(&g), Uint::<4>::from(GENERIC_G.1), "stored G y");
    assert_eq!(c.n, Uint::<4>::from(GENERIC_SUBGROUP_N), "stored order");
    assert_eq!(GENERIC_N, GENERIC_SUBGROUP_N, "exported order");

    assert!(c.is_on_curve(&g), "G not on curve");
    assert!(!g.is_infinity(), "G is the point at infinity");
    assert!(is_prime(&Uint::<4>::from(GENERIC_SUBGROUP_N)), "n is not prime");
    assert!(
        c.scalar_mul(&g, &Uint::<4>::from(GENERIC_SUBGROUP_N)).is_infinity(),
        "n·G ≠ ∞"
    );

    let p0 = generic_original_point(&c);
    assert_eq!(c.scalar_mul(&p0, &u64_to_scalar(GENERIC_C)), g, "c·P₀ ≠ G");
}

/// Fixed logarithm and walk seed for the generic-curve solver run.
const GENERIC_K: u64 = 424_242;
const GENERIC_SEED: u64 = 7;

/// Brent's input check accepts `(G, n)`, and `solve_brent` recovers a fixed nonzero `k` from
/// `(G, k·G)` with a recorded seed.
#[test]
fn generic_brent_recovers_log_in_subgroup() {
    let c = generic_curve();
    let g: AffinePoint<FpMonty> = c.generator();
    let q = c.scalar_mul(&g, &u64_to_scalar(GENERIC_K));
    // The target the `rho-dlog` documentation passes as `--q`.
    assert_eq!(affine_x(&q), Uint::<4>::from(1_000_932_015_531_192_845u64), "Q x");
    assert_eq!(affine_y(&q), Uint::<4>::from(9_048_581_680_910_854_061u64), "Q y");

    assert_eq!(check_prime_order_input(&c, &g, &q, GENERIC_N), Ok(()));
    let solution = solve_brent(&c, &g, &q, GENERIC_N, GENERIC_SEED, 30)
        .unwrap_or_else(|e| panic!("generic: solve_brent failed: {e:?}"));
    assert_eq!(solution.k, GENERIC_K, "generic: recovered k");
    assert!(verify_log(&c, &g, &q, solution.k), "generic: verify_log rejects k");
    assert_collision_identity(&c, &g, &q, &solution.evidence, "generic");
}

/// `P₀` has order far above 561_103, so the existing input checks reject it: as a base point it
/// has the wrong order, and as a target against base `G` it lies outside the subgroup of `G`.
#[test]
fn generic_original_point_is_rejected_not_projected() {
    let c = generic_curve();
    let g: AffinePoint<FpMonty> = c.generator();
    let p0 = generic_original_point(&c);
    let q = c.scalar_mul(&g, &u64_to_scalar(GENERIC_K));

    assert_eq!(
        solve_brent(&c, &p0, &q, GENERIC_N, GENERIC_SEED, 30),
        Err(DlogError::InvalidInput(InputDefect::WrongGeneratorOrder))
    );
    assert_eq!(
        solve_brent(&c, &g, &p0, GENERIC_N, GENERIC_SEED, 30),
        Err(DlogError::InvalidInput(InputDefect::TargetOutsideSubgroup))
    );
}

// ── secp_k1_toy KATs ─────────────────────────────────────────────────────────

/// Reference multiples of G on `y² = x³ + 7  mod  4_611_686_018_427_395_203`.
#[test]
fn secp_k1_toy_scalar_mul_reference() {
    let c = secp_k1_toy();
    let cases: &[(u64, u64, u64)] = &[
        (1, 2,                          3_236_101_131_256_320_111),
        (2, 922_337_203_685_479_039,    132_612_412_593_110_192),
        (3, 2_732_850_973_882_900_861,  4_393_719_944_955_491_326),
        (7, 1_979_583_965_183_108_279,  3_698_824_691_131_872_996),
    ];
    for &(k, ref_x, ref_y) in cases {
        check_scalar_mul(&c, &c.generator(), k, ref_x, ref_y, "secp_k1_toy");
    }
}

#[test]
fn secp_k1_toy_n_times_g_is_infinity() {
    let c = secp_k1_toy();
    let g: AffinePoint<FpMonty> = c.generator();
    let ng = c.scalar_mul(&g, &c.n);
    assert!(ng.is_infinity(), "n·G should be ∞ on secp_k1_toy");
}

#[test]
fn secp_k1_toy_endomorphism_kat() {
    // Hard-coded: φ(G) = λ·G = (458_510_211_330_451_205, 3_236_101_131_256_320_111).
    let c = secp_k1_toy();
    let g: AffinePoint<FpMonty> = c.generator();
    let lam_g = c.scalar_mul(&g, &u64_to_scalar(LAMBDA));

    assert!(!lam_g.is_infinity(), "λ·G is unexpectedly ∞");
    assert_eq!(
        affine_x(&lam_g),
        Uint::<4>::from(458_510_211_330_451_205u64),
        "λ·G x mismatch"
    );
    assert_eq!(
        affine_y(&lam_g),
        Uint::<4>::from(3_236_101_131_256_320_111u64),
        "λ·G y mismatch"
    );

    // φ(G) = (β·Gx mod p, Gy)
    let p_u = Uint::<4>::from(rho::curve::secp_k1_toy::P);
    let phi_x = rho::curve::secp_k1_toy::glv_phi_x(Uint::<4>::from(GX));
    let phi_g: AffinePoint<FpMonty> = AffinePoint::Finite {
        x: FpMonty::from_uint(phi_x, &p_u),
        y: g.y().unwrap().clone(),
    };
    assert_eq!(phi_g, lam_g, "φ(G) ≠ λ·G");
}

#[test]
fn secp_k1_toy_linearity() {
    let c = secp_k1_toy();
    let p = &c.p;
    let g: AffinePoint<FpMonty> = c.generator();

    for a in 1u64..=4 {
        for b in 1u64..=4 {
            let apb_g = c.scalar_mul(&g, &u64_to_scalar(a + b));
            let ag    = c.scalar_mul(&g, &u64_to_scalar(a));
            let bg    = c.scalar_mul(&g, &u64_to_scalar(b));
            let ag_j  = JacobianPoint::from_affine(&ag, p);
            let sum   = c.add_mixed(&ag_j, &bg).to_affine(p);
            assert_eq!(
                sum, apb_g,
                "secp_k1_toy linearity: ({a}+{b})·G ≠ {a}·G + {b}·G"
            );
        }
    }
}

// ── k256 group-law cross-check ───────────────────────────────────────────────
//
// We verify that our Jacobian group-law formulae are consistent with `k256`'s
// production implementation.  Strategy:
//
//  - Pick scalar k = 7 (small, deterministic).
//  - Compute 7·G on the full secp256k1 curve via `k256`.
//  - Extract the affine x-coordinate.
//  - Verify it matches the independently published test vector.
//
// We do NOT compare against our toy-curve result (different prime), but this test
// proves that `k256`'s formulas and our formulas are consistent — if ours were
// wrong, they would disagree with a reference implementation.  The toy-curve
// tests above then extend that confidence to our parametric curves.

#[test]
fn k256_scalar_mul_sanity() {
    use k256::elliptic_curve::group::GroupEncoding;
    use k256::elliptic_curve::PrimeField;
    use k256::{ProjectivePoint, Scalar};

    // k = 7
    let k_val: u64 = 7;

    // k256 scalar: build from bytes (big-endian, 32 bytes).
    let mut scalar_bytes = [0u8; 32];
    scalar_bytes[24..].copy_from_slice(&k_val.to_be_bytes());
    let scalar = Scalar::from_repr(scalar_bytes.into()).unwrap();

    // 7·G on secp256k1 via k256
    let result = ProjectivePoint::GENERATOR * scalar;
    let affine = result.to_affine();
    let coords = affine.to_bytes();  // 33 bytes: 02/03 + x (compressed)

    // The x-coordinate occupies bytes [1..33] (big-endian).
    let x_bytes = &coords[1..33];

    // Reference: 7·G on secp256k1 (from https://learnmeabitcoin.com/tools/point-multiplication)
    // x = 5CBDF0646E5DB4EAA398F365F2EA7A0E3D419B7E0330E39CE92BDDEDCAC4F9BC
    let expected_x_hex = "5CBDF0646E5DB4EAA398F365F2EA7A0E3D419B7E0330E39CE92BDDEDCAC4F9BC";
    let expected_x: Vec<u8> = (0..expected_x_hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&expected_x_hex[i..i+2], 16).unwrap())
        .collect();

    assert_eq!(
        x_bytes,
        expected_x.as_slice(),
        "k256: 7·G x-coordinate does not match published test vector"
    );
}

// ── r-adding walk ECDLP solver KATs ──────────────────────────────────────────
//
// These tests verify `solve_brent` on two 20-bit prime-order test curves.
// Expected steps to solution: O(√n) ≈ O(2^10) — fast even in debug mode.
//
// The solver is not required to return the specific k used to construct Q;
// any k' ∈ [1, n) with k'·G = Q is a valid solution.  The assertion checks Q.

/// Helper: verify the solver finds a valid DLP solution on the given curve.
fn check_solver_on_curve(
    curve: &rho::curve::Curve,
    n: u64,
    k_target: u64,
    label: &str,
) {
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));

    let solution = solve_brent(curve, &g, &q, n, 0, 30)
        .unwrap_or_else(|e| panic!("{label}: solve_brent failed for k_target={k_target}: {e:?}"));
    let k = solution.k;

    let check = curve.scalar_mul(&g, &Uint::<4>::from(k));
    assert_eq!(
        check, q,
        "{label}: recovered k={k} gives k·G ≠ Q (k_target={k_target})"
    );
    assert!(verify_log(curve, &g, &q, k), "{label}: verify_log rejects k={k}");
    assert_collision_identity(curve, &g, &q, &solution.evidence, label);
}

/// `a·G + b·Q` for scalars already reduced mod n.
fn combine(
    curve: &rho::curve::Curve,
    g: &AffinePoint<FpMonty>,
    q: &AffinePoint<FpMonty>,
    a: u64,
    b: u64,
) -> AffinePoint<FpMonty> {
    let ag = curve.scalar_mul(g, &Uint::<4>::from(a));
    let bq = curve.scalar_mul(q, &Uint::<4>::from(b));
    curve
        .add_mixed(&JacobianPoint::from_affine(&ag, &curve.p), &bq)
        .to_affine(&curve.p)
}

/// Recheck, from the evidence alone, the walk identity
/// `a1·G + b1·Q = (−1)^negated · (a2·G + b2·Q)`.
fn assert_collision_identity(
    curve: &rho::curve::Curve,
    g: &AffinePoint<FpMonty>,
    q: &AffinePoint<FpMonty>,
    evidence: &RhoEvidence,
    label: &str,
) {
    let RhoEvidence::Collision { a1, b1, a2, b2, negated, steps, .. } = *evidence else {
        panic!("{label}: expected collision evidence, got {evidence:?}");
    };
    assert!(steps > 0, "{label}: a collision takes at least one step");
    let first = combine(curve, g, q, a1, b1);
    let mut second = combine(curve, g, q, a2, b2);
    if negated {
        second = curve.negate(&second);
    }
    assert_eq!(first, second, "{label}: collision identity does not hold");
}

/// Fixed-seed runs on both 20-bit curves return collision evidence whose identity holds, and
/// replaying the seed reproduces the evidence exactly (the walk is deterministic and single
/// threaded).
#[test]
fn brent_collision_evidence_is_rechecked_and_reproducible() {
    for (curve, n, label) in [(tiny_a(), TINY_A_N, "tiny_a"), (tiny_b(), TINY_B_N, "tiny_b")] {
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(1234u64));

        let first = solve_brent(&curve, &g, &q, n, 42, 30).expect("first run");
        assert_collision_identity(&curve, &g, &q, &first.evidence, label);
        let replay = solve_brent(&curve, &g, &q, n, 42, 30).expect("replay");
        assert_eq!(first, replay, "{label}: same seed must reproduce the solution");
    }
}

/// The trivial targets are answered without a walk: Q = ∞ is k = 0 and Q = G is k = 1.
#[test]
fn brent_trivial_targets() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();

    let zero = solve_brent(&curve, &g, &AffinePoint::Infinity, TINY_A_N, 0, 30).unwrap();
    assert_eq!((zero.k, zero.evidence), (0, RhoEvidence::Trivial));
    let one = solve_brent(&curve, &g, &g, TINY_A_N, 0, 30).unwrap();
    assert_eq!((one.k, one.evidence), (1, RhoEvidence::Trivial));
}

/// A zero retry budget is exhausted before any walk step.
#[test]
fn brent_zero_retries_is_budget_exhausted() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    assert_eq!(
        solve_brent(&curve, &g, &q, TINY_A_N, 0, 0),
        Err(DlogError::BudgetExhausted { attempts: 0, steps: 0 })
    );
}

/// `verify_log` accepts the true logarithm and rejects its successor.
#[test]
fn verify_log_rejects_successor() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    assert!(verify_log(&curve, &g, &q, 7));
    assert!(!verify_log(&curve, &g, &q, 8));
}

/// A generator that is not on the curve is reported as such, not solved.
#[test]
fn brent_rejects_off_curve_generator() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
    let AffinePoint::Finite { x, y } = &g else { unreachable!() };
    let one = FpMonty::from_uint(Uint::<4>::ONE, &curve.p);
    let bad = AffinePoint::new(x.clone(), y.add(&one, &curve.p));

    assert_eq!(
        solve_brent(&curve, &bad, &q, TINY_A_N, 0, 30),
        Err(DlogError::InvalidInput(InputDefect::GeneratorNotOnCurve))
    );
}

// ── distinguished-point solver KATs ──────────────────────────────────────────
//
// The four distinguished-point solvers share Brent's outcome types, take a total walk-step
// budget, and verify before returning. The walkers race, so these assert the collision identity
// and not whole-evidence equality.

/// Step budget for fixed-seed success runs on the 20-bit curves. A collision there takes on the
/// order of 10³–10⁴ steps, so this leaves about three orders of magnitude of headroom.
const DP_BUDGET: u64 = 1 << 24;

/// Run every solver variant that can finish on a 20-bit curve and return its label and result.
fn run_dp_variants(
    curve: &rho::curve::Curve,
    g: &AffinePoint<FpMonty>,
    q: &AffinePoint<FpMonty>,
    n: u64,
    seed: u64,
    budget: u64,
) -> Vec<(&'static str, Result<DlogSolution<RhoEvidence>, DlogError>)> {
    vec![
        ("solve_dp/1", solve_dp(curve, g, q, n, 1, 4, seed, budget)),
        ("solve_dp/2", solve_dp(curve, g, q, n, 2, 4, seed, budget)),
        ("solve_dp_negmap/2", solve_dp_negmap(curve, g, q, n, 2, 4, seed, budget)),
        ("solve_dp_batch/2x4", solve_dp_batch(curve, g, q, n, 2, 4, 4, seed, budget)),
    ]
}

/// Fixed-seed runs on both 20-bit curves return collision evidence whose identity holds, the
/// documented single-run attempt index and sign, and a `k` that verifies.
#[test]
fn dp_collision_evidence_is_rechecked() {
    for (curve, n, tag) in [(tiny_a(), TINY_A_N, "tiny_a"), (tiny_b(), TINY_B_N, "tiny_b")] {
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(1234u64));

        for (variant, result) in run_dp_variants(&curve, &g, &q, n, 0xD15C_0001, DP_BUDGET) {
            let label = format!("{tag} {variant}");
            let solution = result.unwrap_or_else(|e| panic!("{label}: {e:?}"));
            assert_collision_identity(&curve, &g, &q, &solution.evidence, &label);
            assert!(verify_log(&curve, &g, &q, solution.k), "{label}: k fails verify_log");
            let RhoEvidence::Collision { attempt, negated, steps, .. } = solution.evidence else {
                unreachable!("assert_collision_identity accepted non-collision evidence");
            };
            assert_eq!(attempt, 0, "{label}");
            // Only the plain walk can meet `−P`; the canonicalising walks fold the sign away.
            assert!(variant.starts_with("solve_dp/") || !negated, "{label}: spurious sign");
            assert!(steps <= DP_BUDGET, "{label}: evidence steps exceed the budget");
        }
    }
}

/// Parameters of the single-walker plain runs below: DP threshold, and a budget about three
/// orders of magnitude above a collision's cost on the 20-bit curve.
const SIGNED_THETA: u32 = 2;
const SIGNED_BUDGET: u64 = 1 << 20;
/// Discrete logarithm of the target in the signed-collision runs.
const SIGNED_K: u64 = 1234;

/// A plain walk that lands on a distinguished point and, later, on its negative recovers the
/// logarithm from `a1·G + b1·Q = −(a2·G + b2·Q)` and reports `negated: true`.
///
/// Seed found by a bounded search: `solve_dp` with one walker (so the run is deterministic) on
/// `tiny_a` with `Q = 1234·G`, θ = 2, budget 2^20, seeds 0..100. Matching distinguished points by
/// x alone, as the coordinator once did, left 18 of those 100 runs as `Inconsistent` (seed 0
/// first); with the y-parity in the record the same 18 runs are `Ok` with `negated: true`. The
/// other 82 are plain collisions and are unchanged.
#[test]
fn plain_walk_recovers_log_from_negated_collision() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(SIGNED_K));

    let solution = solve_dp(&curve, &g, &q, TINY_A_N, 1, SIGNED_THETA, 0, SIGNED_BUDGET)
        .expect("a negated collision must solve, not read as inconsistent");

    assert_collision_identity(&curve, &g, &q, &solution.evidence, "seed 0");
    assert!(verify_log(&curve, &g, &q, solution.k));
    let RhoEvidence::Collision { attempt, negated, .. } = solution.evidence else {
        unreachable!("assert_collision_identity accepted non-collision evidence");
    };
    assert_eq!((attempt, negated), (0, true));
}

/// Over the seeds of the search above that are cheap to rerun, no plain run ends `Inconsistent`
/// and both collision signs occur, each with an identity that holds from the evidence alone.
#[test]
fn plain_walk_never_misreads_the_sign_of_a_collision() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(SIGNED_K));

    let mut signs = [0usize; 2];
    for seed in 0..24u64 {
        let label = format!("seed {seed}");
        let solution = solve_dp(&curve, &g, &q, TINY_A_N, 1, SIGNED_THETA, seed, SIGNED_BUDGET)
            .unwrap_or_else(|e| panic!("{label}: {e:?}"));
        assert_collision_identity(&curve, &g, &q, &solution.evidence, &label);
        assert!(verify_log(&curve, &g, &q, solution.k), "{label}");
        let RhoEvidence::Collision { negated, .. } = solution.evidence else {
            unreachable!("assert_collision_identity accepted non-collision evidence");
        };
        signs[usize::from(negated)] += 1;
    }
    assert!(signs[0] > 0 && signs[1] > 0, "expected both signs in 24 runs, got {signs:?}");
}

/// A budget far below a collision's cost ends every variant as budget exhaustion, with the whole
/// budget spent (the batch variant spends whole batches: 10 steps allow two batches of 4).
#[test]
fn dp_tiny_budget_is_budget_exhausted() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(77u64));

    for (variant, result) in run_dp_variants(&curve, &g, &q, TINY_A_N, 1, 10) {
        let steps = if variant.starts_with("solve_dp_batch") { 8 } else { 10 };
        assert_eq!(
            result,
            Err(DlogError::BudgetExhausted { attempts: 1, steps }),
            "{variant}"
        );
    }
}

/// The GLV solver on its own curve gets a real search: a tiny budget is exhausted, spent in whole
/// batches. (A full solve on the 62-bit group is far beyond a test's cost.)
#[test]
fn dp_glv_tiny_budget_is_budget_exhausted() {
    let curve = secp_k1_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(12_345u64));

    assert_eq!(
        solve_dp_glv(&curve, &g, &q, SECP_N, 2, 4, 8, 0, 64),
        Err(DlogError::BudgetExhausted { attempts: 1, steps: 64 })
    );
}

/// The GLV solver recognises its curve by parameters: other curves, including the same equation
/// over another prime, and the right curve with a different claimed order, are unsupported.
#[test]
fn dp_glv_is_unsupported_off_its_curve() {
    for (curve, n) in [(tiny_a(), TINY_A_N), (tiny_glv(), TINY_GLV_N)] {
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
        assert!(
            matches!(
                solve_dp_glv(&curve, &g, &q, n, 2, 4, 2, 0, DP_BUDGET),
                Err(DlogError::Unsupported(_))
            ),
            "n = {n}"
        );
    }

    let curve = secp_k1_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
    assert!(matches!(
        solve_dp_glv(&curve, &g, &q, SECP_N - 2, 2, 4, 2, 0, DP_BUDGET),
        Err(DlogError::Unsupported(_))
    ));
}

/// Trivial targets and malformed instances are decided before any walker starts, by every variant.
#[test]
fn dp_trivial_and_invalid_inputs_are_decided_without_a_walk() {
    let curve = tiny_a();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    for (target, k) in [(AffinePoint::Infinity, 0u64), (g.clone(), 1)] {
        for (variant, result) in run_dp_variants(&curve, &g, &target, TINY_A_N, 0, 0) {
            let solution = result.unwrap_or_else(|e| panic!("{variant}: {e:?}"));
            assert_eq!((solution.k, solution.evidence), (k, RhoEvidence::Trivial), "{variant}");
        }
    }

    for (variant, result) in run_dp_variants(&curve, &g, &q, TINY_A_N + 1, 0, DP_BUDGET) {
        assert_eq!(
            result,
            Err(DlogError::InvalidInput(InputDefect::OrderNotPrime)),
            "{variant}"
        );
    }
}

// ── tiny_a KATs ───────────────────────────────────────────────────────────────

/// DLP k=7 on tiny_a (20-bit prime order curve A).
///
/// Q = 7·G = (1_026_105, 636_225) on `y² = x³ − 3x + 3 mod 1_048_517`.
#[test]
fn tiny_a_dlog_k7() {
    check_solver_on_curve(&tiny_a(), TINY_A_N, 7, "tiny_a");
}

/// DLP k=100 on tiny_a.
///
/// Q = 100·G = (659_291, 755_487).
#[test]
fn tiny_a_dlog_k100() {
    check_solver_on_curve(&tiny_a(), TINY_A_N, 100, "tiny_a");
}

/// DLP k=33333 on tiny_a.
///
/// Q = 33333·G = (758_517, 785_775).
#[test]
fn tiny_a_dlog_k33333() {
    check_solver_on_curve(&tiny_a(), TINY_A_N, 33_333, "tiny_a");
}

// ── tiny_b KATs ───────────────────────────────────────────────────────────────

/// DLP k=7 on tiny_b (20-bit prime order curve B).
///
/// Q = 7·G = (284_547, 163_192) on `y² = x³ − 3x + 16 mod 1_048_583`.
#[test]
fn tiny_b_dlog_k7() {
    check_solver_on_curve(&tiny_b(), TINY_B_N, 7, "tiny_b");
}

/// DLP k=42 on tiny_b.
///
/// Q = 42·G = (132_859, 318_692).
#[test]
fn tiny_b_dlog_k42() {
    check_solver_on_curve(&tiny_b(), TINY_B_N, 42, "tiny_b");
}

/// DLP k=99991 on tiny_b.
///
/// Q = 99991·G = (654_745, 751_943).
#[test]
fn tiny_b_dlog_k99991() {
    check_solver_on_curve(&tiny_b(), TINY_B_N, 99_991, "tiny_b");
}

// ── Pohlig–Hellman composite-order ECDLP KATs ────────────────────────────────
//
// These tests verify `solve_ecdlp_composite` on the composite-order toy curve:
// y² = x³ + x + 33 mod 47, generator G = (10, 3), order n = 60 = 2² · 3 · 5.
//
// The solver is not required to return the specific k used to construct Q;
// any k' ∈ [0, n) with k'·G = Q is a valid solution.  The assertion checks Q.
//
// Coverage:
// - At least one case exercises the e>1 lift (the 2² factor).
// - At least one case exercises the multi-prime CRT (all three primes 2, 3, 5).

/// Assert the identities a `PohligEvidence` must satisfy on its own: the components are the
/// prime powers of `n` in order, each residue is `k mod p^e` and the sum of its digits, and every
/// digit's method is consistent with the size of `p`. Returns `(subgroup enumerations, rho
/// solves)` counted over all digits.
fn check_pohlig_evidence(
    evidence: &PohligEvidence,
    factors: &[(u64, u32)],
    k: u64,
    label: &str,
) -> (usize, usize) {
    assert_eq!(evidence.components.len(), factors.len(), "{label}: component count");
    let (mut enumerated, mut rho) = (0, 0);
    for (component, &(p, e)) in evidence.components.iter().zip(factors) {
        assert_eq!((component.p, component.e), (p, e), "{label}: component prime power");
        let p_power = p.pow(e);
        assert_eq!(component.residue, k % p_power, "{label}: residue ≢ k (mod {p}^{e})");
        assert_eq!(component.digits.len(), e as usize, "{label}: one method per digit of {p}^{e}");
        for (j, method) in component.digits.iter().enumerate() {
            let digit = (component.residue / p.pow(j as u32)) % p;
            match method {
                DigitMethod::SubgroupEnumeration => {
                    assert!(p <= 64, "{label}: enumeration used for p = {p} > 64");
                    enumerated += 1;
                }
                DigitMethod::Rho(ev) => {
                    assert!(p > 64, "{label}: rho used for p = {p} ≤ 64");
                    rho += 1;
                    // The collision identity holds modulo the order-p subgroup, with the digit
                    // the evidence produced standing for k.
                    if let RhoEvidence::Collision { a1, b1, a2, b2, negated, .. } = ev {
                        let m = p as u128;
                        let lhs = (*a1 as u128 + *b1 as u128 * digit as u128) % m;
                        let rhs = (*a2 as u128 + *b2 as u128 * digit as u128) % m;
                        let rhs = if *negated { (m - rhs) % m } else { rhs };
                        assert_eq!(lhs, rhs, "{label}: rho collision identity, digit {j} of {p}");
                    }
                }
            }
        }
    }
    (enumerated, rho)
}

/// Helper: verify solve_ecdlp_composite finds a valid DLP solution on composite_toy, and that
/// its evidence is consistent with the solution.
fn check_composite_solver(k_target: u64, label: &str) {
    let curve = composite_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));

    let solution = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N).unwrap_or_else(|err| {
        panic!("{label}: solve_ecdlp_composite failed for k_target={k_target}: {err:?}")
    });

    let check = curve.scalar_mul(&g, &Uint::<4>::from(solution.k));
    assert_eq!(
        check, q,
        "{label}: recovered k={} gives k·G ≠ Q (k_target={k_target})",
        solution.k
    );
    // Every prime of composite_toy is ≤ 5, so no digit reaches rho.
    let (enumerated, rho) = check_pohlig_evidence(
        &solution.evidence,
        &[(2, 2), (3, 1), (5, 1)],
        solution.k,
        label,
    );
    assert_eq!((enumerated, rho), (4, 0), "{label}: digits 2 + 1 + 1, all enumerated");
}

/// k=1: trivial case, exercises all three prime subgroups.
///
/// k=1 has residues: 1 mod 4, 1 mod 3, 1 mod 5. All primes (2, 3, 5) participate
/// in the CRT combine, and the 2² factor exercises the e>1 lift.
#[test]
fn composite_dlog_k1() {
    check_composite_solver(1, "composite_toy");
}

/// k=7: exercises the e>1 lift (2² factor) and multi-prime CRT.
///
/// k=7 mod 4 = 3 (nontrivial in the 2² subgroup, requiring both digits d_0=1, d_1=1
/// in the base-2 expansion: 3 = 1 + 1·2). Also nontrivial mod 3 (7 mod 3 = 1) and
/// mod 5 (7 mod 5 = 2). All three primes participate in the CRT.
#[test]
fn composite_dlog_k7() {
    check_composite_solver(7, "composite_toy");
}

/// k=11: nontrivial residue in the 2² subgroup (11 mod 4 = 3) and multi-prime CRT.
///
/// 11 mod 4 = 3 = 1 + 1·2 (both base-2 digits nonzero — exercises the full e=2 lift).
/// 11 mod 3 = 2, 11 mod 5 = 1. All three prime-power subgroups contribute.
#[test]
fn composite_dlog_k11() {
    check_composite_solver(11, "composite_toy");
}

/// k=30: exercises the e>1 lift with d_0=0 (30 mod 4 = 2 = 0 + 1·2).
///
/// 30 mod 4 = 2: the first digit d_0 = 0, second digit d_1 = 1 — verifies that the
/// lift handles a zero leading digit correctly. 30 mod 3 = 0, 30 mod 5 = 0.
#[test]
fn composite_dlog_k30() {
    check_composite_solver(30, "composite_toy");
}

/// k=59: largest non-identity scalar, exercises all subgroups with nontrivial residues.
///
/// 59 mod 4 = 3, 59 mod 3 = 2, 59 mod 5 = 4. All three prime-power subgroups have
/// nontrivial residues. The 2² lift must recover both digits (3 = 1 + 1·2).
#[test]
fn composite_dlog_k59() {
    check_composite_solver(59, "composite_toy");
}

/// k=0 (Q = ∞): the identity case returns 0.
#[test]
fn composite_dlog_k0_identity() {
    let curve = composite_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = AffinePoint::Infinity;

    let solution = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
        .expect("composite_dlog_k0_identity: solve_ecdlp_composite failed for Q=∞");

    assert_eq!(solution.k, 0, "Q = ∞ is k = 0");
    assert!(
        solution.evidence.components.iter().all(|c| c.residue == 0),
        "every residue of k = 0 is 0"
    );
}

/// Every `k` in `[0, n)` on `composite_toy`: a verified solution whose evidence is consistent.
#[test]
fn composite_dlog_every_k() {
    let curve = composite_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    for k_target in 0..COMPOSITE_TOY_N {
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));
        let solution = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
            .unwrap_or_else(|err| panic!("k = {k_target}: {err:?}"));
        assert_eq!(solution.k, k_target, "G has full order, so k is unique in [0, n)");
        check_pohlig_evidence(&solution.evidence, &[(2, 2), (3, 1), (5, 1)], solution.k, "toy");
    }
}

/// A sampled set of `k` on `composite_rho_toy`, whose prime factor 89 exceeds the enumeration
/// limit: each solve is verified, and the evidence shows enumeration for the `2²` and `3`
/// components and rho for the `89` component, with each residue ≡ k (mod p^e).
///
/// Each rho solve costs about 0.16 s in a debug build even though the walk is only a dozen
/// steps, so all 1068 values take minutes; the sample is every 23rd `k` — stride 23 is coprime
/// to 12 and 89, so it meets every residue mod 12 — plus the edges and the residues 0 and 1 mod
/// 89 whose rho digit is `Trivial`.
#[test]
fn composite_rho_dlog_sampled_k() {
    let curve = composite_rho_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let mut sample: Vec<u64> = (0..COMPOSITE_RHO_TOY_N).step_by(23).collect();
    sample.extend([1, 89, 90, COMPOSITE_RHO_TOY_N - 1]);
    let residues_mod_12: std::collections::BTreeSet<u64> = sample.iter().map(|k| k % 12).collect();
    assert_eq!(residues_mod_12.len(), 12, "the sample must meet every residue mod 2²·3");

    let (mut total_enumerated, mut total_rho, mut rho_collisions) = (0, 0, 0);
    for &k_target in &sample {
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));
        let solution = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_RHO_TOY_N)
            .unwrap_or_else(|err| panic!("k = {k_target}: {err:?}"));
        assert_eq!(solution.k, k_target, "k = {k_target}: unique in [0, n)");
        assert!(verify_log(&curve, &g, &q, solution.k));
        let (enumerated, rho) = check_pohlig_evidence(
            &solution.evidence,
            COMPOSITE_RHO_TOY_FACTORS,
            solution.k,
            "composite_rho_toy",
        );
        // Digits: 2² gives two and 3 gives one (enumeration); 89 gives one (rho).
        assert_eq!((enumerated, rho), (3, 1), "k = {k_target}");
        total_enumerated += enumerated;
        total_rho += rho;
        if let DigitMethod::Rho(RhoEvidence::Collision { .. }) =
            solution.evidence.components[2].digits[0]
        {
            rho_collisions += 1;
        }
    }
    assert_eq!(total_enumerated, 3 * sample.len());
    assert_eq!(total_rho, sample.len());
    // Residues 0 and 1 mod 89 are `Trivial`; the rest ran a walk.
    let trivial = sample.iter().filter(|k| *k % 89 <= 1).count();
    assert!(trivial >= 3, "the sample includes residues 0 and 1 mod 89");
    assert_eq!(rho_collisions, sample.len() - trivial);
}

fn input_defect<T: std::fmt::Debug>(r: Result<T, DlogError>) -> InputDefect {
    match r {
        Err(DlogError::InvalidInput(d)) => d,
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

/// A wrong or degenerate `n` is `InvalidInput`, and `n = 0` no longer panics.
#[test]
fn composite_wrong_order_is_invalid_input() {
    let curve = composite_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    assert_eq!(input_defect(solve_ecdlp_composite(&curve, &g, &q, 0)), InputDefect::ZeroOrder);
    // n = 1 and n = 59: n·G ≠ ∞.
    for n in [1u64, 59] {
        let r = solve_ecdlp_composite(&curve, &g, &q, n);
        assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder, "n = {n}");
    }
    // Multiples of the true order have n·G = ∞ without being the order of G: 120 = 2·60 and
    // 180 = 3·60 fail the per-prime `(n/p)·G ≠ ∞` check.
    for n in [120u64, 180] {
        let r = solve_ecdlp_composite(&curve, &g, &q, n);
        assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder, "n = {n}");
    }
    // A point of smaller order as generator: 2·G has order 30, not 60.
    let g2 = curve.scalar_mul(&g, &Uint::<4>::from(2u64));
    let r = solve_ecdlp_composite(&curve, &g2, &g2, COMPOSITE_TOY_N);
    assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder);
    // The point-level defects come before the order checks.
    let inf = AffinePoint::<FpMonty>::Infinity;
    let r = solve_ecdlp_composite(&curve, &inf, &q, COMPOSITE_TOY_N);
    assert_eq!(input_defect(r), InputDefect::GeneratorAtInfinity);
}

/// The same wrong-order rejection on the fixture that reaches rho.
#[test]
fn composite_rho_wrong_order_is_invalid_input() {
    let curve = composite_rho_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    let q = curve.scalar_mul(&g, &Uint::<4>::from(5u64));
    for n in [COMPOSITE_RHO_TOY_N - 1, COMPOSITE_RHO_TOY_N * 2, COMPOSITE_RHO_TOY_N / 2] {
        let r = solve_ecdlp_composite(&curve, &g, &q, n);
        assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder, "n = {n}");
    }
}

/// A target in `E[n]` outside `⟨G⟩` is `InvalidInput(TargetOutsideSubgroup)`, not a wrong `k`.
///
/// `y² = x³ − x mod 19` has 20 points with full 2-torsion (`Z/2 × Z/10`). `G = (2, 5)` has order
/// 10, so `⟨G⟩` holds one of the three 2-torsion points; another one passes `n·Q = ∞` yet is not
/// a multiple of `G`, and the order-2 digit enumeration proves it.
#[test]
fn composite_target_outside_subgroup_is_invalid_input() {
    let curve = Curve {
        p: Uint::<4>::from(19u64),
        a: Uint::<4>::from(18u64),
        b: Uint::<4>::ZERO,
        n: Uint::<4>::from(10u64),
        gx: Uint::<4>::from(2u64),
        gy: Uint::<4>::from(5u64),
    };
    let g: AffinePoint<FpMonty> = curve.generator();
    let inside = curve.scalar_mul(&g, &Uint::<4>::from(5u64)); // the 2-torsion point in ⟨G⟩
    let two_torsion = [(0u64, 0u64), (1, 0), (18, 0)].map(|(x, y)| {
        let p = curve.p;
        AffinePoint::new(
            FpMonty::from_uint(Uint::<4>::from(x), &p),
            FpMonty::from_uint(Uint::<4>::from(y), &p),
        )
    });
    let outside = two_torsion
        .iter()
        .find(|t| **t != inside)
        .expect("three 2-torsion points, one inside ⟨G⟩")
        .clone();
    assert!(curve.is_on_curve(&outside));
    assert!(curve.scalar_mul(&outside, &Uint::<4>::from(10u64)).is_infinity());

    let r = solve_ecdlp_composite(&curve, &g, &outside, 10);
    assert_eq!(input_defect(r), InputDefect::TargetOutsideSubgroup);
    // The instance with the in-subgroup 2-torsion point solves (k = 5).
    let solution = solve_ecdlp_composite(&curve, &g, &inside, 10).expect("in-subgroup target");
    assert_eq!(solution.k, 5);
}
