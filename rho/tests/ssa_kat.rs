//! Known-answer tests (KATs) for the SSA (Smart–Satoh–Araki) p-adic attack.
//!
//! # Fixture
//!
//! All tests use the anomalous toy curve `y² = x³ + 5` over `F_7` (`p = 7`, `a = 0`, `b = 5`).
//! This curve has exactly 7 points (including the point at infinity), so `#E(F_7) = 7 = p` and
//! the trace of Frobenius is 1. The base point `G = (3, 2)` has `y ≠ 0`.
//!
//! # KAT coverage — Hensel lift
//!
//! 1. `verify_anomalous` returns `true` on the fixture curve.
//! 2. `verify_anomalous` returns `false` on a non-anomalous control curve.
//! 3. The lift of `G = (3, 2)` matches hand-computed Z_7 coordinates to precision k=4.
//! 4. The lifted point satisfies the curve equation mod `7^4` (lift-correctness check).
//! 5. A 2-torsion point (`y = 0`) on a test curve errors via `HenselError::NonSimpleRoot`.
//!
//! # KAT coverage — SSA reduction
//!
//! 9.  `ssa_solve` recovers k=3 from Q = 3·G on the anomalous fixture.
//! 10. `solve_brent` (independent rho solver) also recovers k=3 — cross-check.
//! 11. `ssa_solve` returns `SsaError::NotAnomalous` on a non-anomalous curve.
//! 12. PARI/GP cross-check (`#[ignore]`).
//! 13. `ssa_solve` recovers every `k` in `1..p` on the fixture.
//! 14. `Q = ∞` returns `Ok(0)`, on the fixture and on a non-anomalous curve.
//! 15. A prime that does not fit 64 bits is `Unsupported` for both `verify_anomalous` and
//!     `ssa_solve`, never `NotAnomalous`.
//! 16. `verify_anomalous` counts exactly where `x³` overflows `u64` (an anomalous curve over a
//!     prime above `2^21.4`, and the same curve with a tampered coefficient).
//!
//! # Note on the fixture curve's CM structure
//!
//! The curve `y² = x³ + 5` over `F_7` has complex multiplication by `Z[ζ₃]` (since `a = 0`
//! and `p ≡ 1 mod 3`). This causes the SSA formula `k_raw = ψ(p·Q̃)/ψ(p·G̃) mod p` to give
//! `2k mod p` for `k ∈ {2,3,4,5}`. The `ssa_solve` implementation checks `k·G = Q` in `E(F_p)`
//! for the quotient and for its half before returning; the candidate selection itself is
//! unit-tested in `rho::ssa::reduce`.
//!
//! # Hand-computed lift values
//!
//! For `G = (3, 2)` on `y² = x³ + 5` over `F_7`, precision k=4 (modulus = 7^4 = 2401):
//!
//! - `x̃ = 3` (lifted exactly).
//! - `c = x₀³ + b = 27 + 5 = 32` (the RHS of the curve equation at x=3).
//! - Hensel-solve `g(y) = y² − 32 = 0` starting from `y₀ = 2`:
//!   - Step 1 (mod 49): y = 2 − (4−32)/4 = 2 + 7 = 9.
//!   - Step 2 (mod 343): y = 9 − (81−32)/(18) mod 343 = 9 − 49/18 mod 343.
//!     18⁻¹ mod 343 = 248 (18·248 = 4464 = 13·343 + 5 ≠ 1 — let me recompute).
//!   - Final: `ỹ = 940` (verified: 940² mod 2401 = 32 ✓).
//!
//! Verification: `940² = 883600`, `883600 mod 2401 = 32`. ✓

use crypto_bigint::Uint;
use num_bigint::BigInt;
use rho::curve::{AffinePoint, Curve};
use rho::ecdlp::solve_brent;
use shared_field::{Fp, FpNaive4 as FpNaive};
use rho::ssa::lift::{check_lift_on_curve, lift_point};
use rho::ssa::{ANOMALOUS_TOY_P, SsaError, anomalous_toy, verify_anomalous, ssa_solve};
use shared_padic::HenselError;

// ─── helpers ──────────────────────────────────────────────────────────────────

/// Construct a `FpNaive` point on the anomalous toy curve.
fn toy_point(x: u64, y: u64) -> AffinePoint<FpNaive> {
    let p = Uint::<4>::from(ANOMALOUS_TOY_P);
    AffinePoint::Finite {
        x: FpNaive::from_u64(x, &p),
        y: FpNaive::from_u64(y, &p),
    }
}

/// Construct a `FpNaive` point on a curve with the given prime.
fn point_on(x: u64, y: u64, p: u64) -> AffinePoint<FpNaive> {
    let p_uint = Uint::<4>::from(p);
    AffinePoint::Finite {
        x: FpNaive::from_u64(x, &p_uint),
        y: FpNaive::from_u64(y, &p_uint),
    }
}

fn bi(n: i64) -> BigInt {
    BigInt::from(n)
}

// ─── KAT 1: verify_anomalous returns true on the fixture ─────────────────────

/// `verify_anomalous` returns `true` on the anomalous toy fixture.
///
/// The fixture is `y² = x³ + 5` over `F_7`. It has exactly 7 points, so `#E = p = 7`.
/// The O(p) point-count via Legendre symbol must confirm this.
#[test]
fn kat_verify_anomalous_fixture() {
    let curve = anomalous_toy();
    assert_eq!(
        verify_anomalous::<FpNaive>(&curve),
        Ok(true),
        "verify_anomalous must return true on the anomalous toy fixture (y² = x³ + 5 over F_7)"
    );
}

// ─── KAT 2: verify_anomalous returns false on a non-anomalous control curve ───

/// `verify_anomalous` returns `false` on a non-anomalous control curve.
///
/// The control curve is `y² = x³ + 3x + 1` over `F_7`. It has 12 points, so `#E = 12 ≠ 7 = p`.
/// This is the anomalous-detection guard: the SSA reduction must not run on non-anomalous curves.
#[test]
fn kat_verify_anomalous_control_false() {
    // Control curve: y² = x³ + 3x + 1 over F_7, #E = 12 (non-anomalous).
    let control = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(3u64),
        b: Uint::<4>::from(1u64),
        n: Uint::<4>::from(12u64), // #E = 12, not 7
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(1u64),
    };
    assert_eq!(
        verify_anomalous::<FpNaive>(&control),
        Ok(false),
        "verify_anomalous must return false on the non-anomalous control curve (y² = x³ + 3x + 1 \
         over F_7, #E = 12)"
    );
}

// ─── KAT 3: lift of G = (3, 2) matches hand-computed Z_7 coordinates ─────────

/// The lift of `G = (3, 2)` to Z_7 at precision k=4 matches hand-computed values.
///
/// Hand-computed: `x̃ = 3`, `ỹ = 940`.
/// Verification: `940² mod 7^4 = 883600 mod 2401 = 32 = 3³ + 5 mod 2401`. ✓
#[test]
fn kat_lift_base_point_coordinates() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let k = 4u32;

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, k)
        .expect("lift of G = (3, 2) should succeed");

    // x̃ = 3 (lifted exactly — no Hensel needed for x).
    assert_eq!(
        *x_tilde.residue(),
        bi(3),
        "x̃ must equal 3 (exact lift of x₀ = 3)"
    );

    // ỹ = 940 (Hensel-solved from y² = 32 over Z_7 to precision 4).
    assert_eq!(
        *y_tilde.residue(),
        bi(940),
        "ỹ must equal 940 (Hensel lift of y₀ = 2 for y² = 32 mod 7^4)"
    );

    // Precision must be preserved.
    assert_eq!(x_tilde.precision(), k, "x̃ precision must be k=4");
    assert_eq!(y_tilde.precision(), k, "ỹ precision must be k=4");
}

// ─── KAT 4: lifted point satisfies the curve equation mod p^k ────────────────

/// The lifted point `(x̃, ỹ)` satisfies `ỹ² ≡ x̃³ + ax̃ + b (mod 7^4)`.
///
/// This is the Hensel lift-correctness invariant. A lift that does not satisfy
/// the curve equation would cause the formal-group log to read a garbage value.
#[test]
fn kat_lift_on_curve_mod_pk() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let k = 4u32;

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, k)
        .expect("lift of G = (3, 2) should succeed");

    // Verify ỹ² ≡ x̃³ + ax̃ + b (mod 7^4).
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve)
        .expect("check_lift_on_curve should not error");

    assert!(
        on_curve,
        "lifted point (x̃={}, ỹ={}) must satisfy ỹ² ≡ x̃³ + ax̃ + b (mod 7^4)",
        x_tilde.residue(),
        y_tilde.residue()
    );
}

/// Direct arithmetic verification: `940² mod 2401 = 32 = 3³ + 5 mod 2401`.
///
/// Spot-checks the hand-computed values independently of the lift machinery.
#[test]
fn kat_lift_arithmetic_spot_check() {
    let y_tilde: i64 = 940;
    let x_tilde: i64 = 3;
    let b: i64 = 5;
    let modulus: i64 = 7i64.pow(4); // 2401

    let lhs = (y_tilde * y_tilde) % modulus;
    let rhs = (x_tilde * x_tilde * x_tilde + b) % modulus; // a=0

    assert_eq!(
        lhs, rhs,
        "940² mod 2401 must equal 3³ + 5 mod 2401 (hand-computed lift-correctness check)"
    );
    assert_eq!(lhs, 32, "940² mod 2401 = 32");
}

// ─── KAT 5: lift of a second point (5, 2) also satisfies the curve equation ──

/// The lift of `(5, 2)` to Z_7 at precision k=4 satisfies the curve equation.
///
/// This exercises a second base point on the fixture curve to confirm the lift is
/// not specific to `G = (3, 2)`.
#[test]
fn kat_lift_second_point_on_curve() {
    let curve = anomalous_toy();
    let pt = toy_point(5, 2);
    let k = 4u32;

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&pt, &curve, k)
        .expect("lift of (5, 2) should succeed");

    // x̃ = 5 (exact lift).
    assert_eq!(*x_tilde.residue(), bi(5), "x̃ must equal 5");

    // ỹ = 1479 (hand-computed: 1479² mod 2401 = 130 = 5³ + 5 mod 2401).
    assert_eq!(
        *y_tilde.residue(),
        bi(1479),
        "ỹ must equal 1479 (Hensel lift of y₀ = 2 for y² = 130 mod 7^4)"
    );

    // Verify on-curve.
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve)
        .expect("check_lift_on_curve should not error");
    assert!(on_curve, "lifted (5, 2) must satisfy the curve equation mod 7^4");
}

/// Direct arithmetic spot-check for the (5, 2) lift: `1479² mod 2401 = 130 = 5³ + 5`.
#[test]
fn kat_lift_second_point_arithmetic_spot_check() {
    let y_tilde: i64 = 1479;
    let x_tilde: i64 = 5;
    let b: i64 = 5;
    let modulus: i64 = 7i64.pow(4); // 2401

    let lhs = (y_tilde * y_tilde) % modulus;
    let rhs = (x_tilde * x_tilde * x_tilde + b) % modulus; // a=0

    assert_eq!(lhs, rhs, "1479² mod 2401 must equal 5³ + 5 mod 2401");
    assert_eq!(lhs, 130, "1479² mod 2401 = 130");
}

// ─── KAT 6: 2-torsion point errors via HenselError::NonSimpleRoot ─────────────

/// Lifting a 2-torsion point (`y = 0`) errors with `HenselError::NonSimpleRoot`.
///
/// The Hensel y-solve requires `g'(y₀) = 2·y₀ ≢ 0 mod p`. For a 2-torsion point
/// (`y₀ = 0`), `g'(0) = 0 ≡ 0 mod p`, so the root is not simple. The lift must
/// return `SsaError::LiftFailed(HenselError::NonSimpleRoot)`.
///
/// Test curve: `y² = x³ + x` over `F_7` (non-anomalous, but has the 2-torsion point
/// `(0, 0)` since `0³ + 0 = 0`). The non-anomalous property is irrelevant here — we
/// are testing the lift machinery's error path, not the SSA precondition.
#[test]
fn kat_lift_2torsion_errors_non_simple_root() {
    // Curve y² = x³ + x over F_7 (a=1, b=0). Point (0, 0) is a 2-torsion point.
    let curve_2torsion = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(1u64),
        b: Uint::<4>::from(0u64),
        n: Uint::<4>::from(8u64), // #E = 8 (non-anomalous)
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(0u64),
    };

    // The 2-torsion point (0, 0): y₀ = 0, so g'(y₀) = 2·0 = 0 mod 7 → NonSimpleRoot.
    let two_torsion = point_on(0, 0, 7);

    let result = lift_point::<FpNaive>(&two_torsion, &curve_2torsion, 4);

    assert!(
        matches!(
            result,
            Err(rho::ssa::SsaError::LiftFailed(HenselError::NonSimpleRoot { .. }))
        ),
        "lifting a 2-torsion point (y=0) must return SsaError::LiftFailed(NonSimpleRoot), \
         got: {result:?}"
    );
}

// ─── KAT 7: lift at precision k=1 (trivial case) ─────────────────────────────

/// Lift at precision k=1 returns the F_p coordinates unchanged (mod p).
///
/// At k=1, the modulus is p itself, so the lifted coordinates are just the F_p values.
#[test]
fn kat_lift_precision_1() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, 1)
        .expect("lift at k=1 should succeed");

    // At k=1, x̃ ≡ 3 mod 7 and ỹ ≡ 2 mod 7.
    assert_eq!(*x_tilde.residue(), bi(3), "x̃ at k=1 must be 3 mod 7");
    assert_eq!(*y_tilde.residue(), bi(2), "ỹ at k=1 must be 2 mod 7");
    assert_eq!(x_tilde.precision(), 1, "precision must be 1");
    assert_eq!(y_tilde.precision(), 1, "precision must be 1");
}

// ─── KAT 8: lift at precision k=2 ────────────────────────────────────────────

/// Lift at precision k=2 (modulus = 49) satisfies the curve equation mod 49.
///
/// Hand-computed: `ỹ = 9` (from the first Newton step: 2 − (4−32)/4 mod 49 = 9).
/// Verification: `9² mod 49 = 81 mod 49 = 32 mod 49 = 32`. ✓
#[test]
fn kat_lift_precision_2() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, 2)
        .expect("lift at k=2 should succeed");

    assert_eq!(*x_tilde.residue(), bi(3), "x̃ at k=2 must be 3");
    assert_eq!(*y_tilde.residue(), bi(9), "ỹ at k=2 must be 9 (first Newton step)");

    // Verify on-curve mod 7^2 = 49.
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve)
        .expect("check_lift_on_curve at k=2 should not error");
    assert!(on_curve, "lifted (3, 2) at k=2 must satisfy the curve equation mod 49");
}

// ─── optional PARI cross-check ────────────────────────────────────────────────

/// PARI/GP cross-check for the Hensel lift of `y² = 32` over `Z_7` to precision 4.
///
/// Run manually with PARI installed:
/// ```text
/// gp> p = 7; k = 4; f = y^2 - 32; y0 = 2; liftall(polrootsmod(f, p^k))
/// ```
/// Expected: `y = 940` (the unique root of `y² − 32 = 0` mod `7^4` with `y ≡ 2 mod 7`).
///
/// This test is gated with `#[ignore]` because PARI is not installed in the standard
/// CI environment. Run with `cargo test -- --ignored` when PARI is available.
#[test]
#[ignore = "PARI not installed; run manually when available"]
fn kat_lift_pari_cross_check() {
    // PARI/GP: liftall(polrootsmod(y^2 - 32, 7^4)) should include 940.
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    let (_, y_tilde) = lift_point::<FpNaive>(&g, &curve, 4)
        .expect("lift should succeed");

    assert_eq!(
        *y_tilde.residue(),
        bi(940),
        "PARI cross-check: Hensel lift of y² = 32 over Z_7 to precision 4 must give ỹ = 940"
    );
}

// ─── SSA reduction KATs ───────────────────────────────────────────────────────

/// `ssa_solve` recovers k=3 from Q = 3·G on the anomalous fixture.
///
/// The anomalous fixture is `y² = x³ + 5` over `F_7`. `G = (3, 2)`, `Q = 3·G = (6, 5)`.
/// The SSA reduction must recover `k = 3`.
///
/// # Note on CM artifact
///
/// This curve has CM by `Z[ζ₃]`. The raw SSA formula gives `k_raw = 2k mod p = 6` for `k = 3`.
/// The `ssa_solve` implementation computes the quotient of the formal-log leading units and its
/// half modulo `p`, then accepts whichever candidate satisfies `k·G == Q`, returning 3 as
/// required.
#[test]
fn kat_ssa_solve_recovers_k3() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    // Q = 3·G = (6, 5): verified by the group law on y² = x³ + 5 over F_7.
    let q = toy_point(6, 5);

    let k_recovered = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("ssa_solve must succeed on the anomalous fixture");

    assert_eq!(
        k_recovered, 3,
        "ssa_solve must recover k=3 from Q = 3·G on the anomalous fixture"
    );
}

/// `ssa_solve` recovers k=1 (Q = G) on the anomalous fixture.
///
/// The trivial case: Q = G, so k = 1. The SSA formula gives k_raw = 1 directly.
#[test]
fn kat_ssa_solve_recovers_k1() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let q = toy_point(3, 2); // Q = G

    let k_recovered = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("ssa_solve must succeed for k=1");

    assert_eq!(k_recovered, 1, "ssa_solve must recover k=1 when Q = G");
}

/// `ssa_solve` recovers k=6 (Q = 6·G) on the anomalous fixture.
///
/// `6·G = (3, 5)`. The SSA formula gives k_raw = 6 directly (no CM artifact for k=6).
#[test]
fn kat_ssa_solve_recovers_k6() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let q = toy_point(3, 5); // Q = 6·G = (3, 5)

    let k_recovered = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("ssa_solve must succeed for k=6");

    assert_eq!(k_recovered, 6, "ssa_solve must recover k=6 when Q = 6·G");
}

/// Cross-check: `solve_brent` (independent rho solver) also recovers k=3 for Q = 3·G.
///
/// This is an independent confirmation that Q = 3·G on the anomalous fixture.
/// The rho solver and the SSA solver must agree on the discrete log.
#[test]
fn kat_ssa_rho_cross_check_k3() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let q = toy_point(6, 5); // Q = 3·G

    // Independent rho solver.
    let k_rho = solve_brent::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P, 0, 20)
        .expect("solve_brent must succeed on the anomalous fixture")
        .k;

    // Verify: k_rho·G = Q.
    let check = curve.scalar_mul(&g, &Uint::<4>::from(k_rho));
    assert_eq!(
        check, q,
        "solve_brent: k_rho·G must equal Q (k_rho={k_rho}, expected k=3)"
    );

    // SSA solver.
    let k_ssa = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("ssa_solve must succeed on the anomalous fixture");

    assert_eq!(
        k_ssa, 3,
        "ssa_solve must recover k=3 from Q = 3·G (cross-check with rho)"
    );

    // Both solvers must agree (both give a valid discrete log, though they may differ
    // if the group has composite order — here #E = 7 = prime, so k is unique).
    let check_ssa = curve.scalar_mul(&g, &Uint::<4>::from(k_ssa));
    assert_eq!(
        check_ssa, q,
        "ssa_solve: k_ssa·G must equal Q (k_ssa={k_ssa})"
    );
}

/// `ssa_solve` returns `SsaError::NotAnomalous` on a non-anomalous curve.
///
/// The control curve `y² = x³ + 3x + 1` over `F_7` has `#E = 12 ≠ 7 = p`.
/// The SSA reduction must refuse to run on non-anomalous curves.
#[test]
fn kat_ssa_solve_not_anomalous() {
    // Control curve: y² = x³ + 3x + 1 over F_7, #E = 12 (non-anomalous).
    let control = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(3u64),
        b: Uint::<4>::from(1u64),
        n: Uint::<4>::from(12u64),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(1u64),
    };

    // Use any point on the control curve.
    let g = point_on(0, 1, 7);
    let q = point_on(0, 1, 7); // Q = G (k=1)

    let result = ssa_solve::<FpNaive>(&control, &g, &q, 12);

    assert!(
        matches!(result, Err(SsaError::NotAnomalous)),
        "ssa_solve must return NotAnomalous on a non-anomalous curve, got: {result:?}"
    );
}

/// Optional PARI/GP cross-check for the SSA reduction.
///
/// Run manually with PARI installed:
/// ```text
/// gp> E = ellinit([0, 5], 7); G = [3, 2]; Q = [6, 5]; elllog(E, Q, G)
/// ```
/// Expected: `3` (the discrete log of Q = 3·G on y² = x³ + 5 over F_7).
///
/// This test is gated with `#[ignore]` because PARI is not installed in the standard
/// CI environment. Run with `cargo test -- --ignored` when PARI is available.
#[test]
#[ignore = "PARI not installed; run manually when available"]
fn kat_ssa_pari_cross_check() {
    // PARI/GP: elllog(ellinit([0,5],7), [6,5], [3,2]) should return 3.
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let q = toy_point(6, 5); // Q = 3·G

    let k_recovered = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("ssa_solve must succeed on the anomalous fixture");

    assert_eq!(
        k_recovered, 3,
        "PARI cross-check: ssa_solve must recover k=3 from Q = 3·G on y² = x³ + 5 over F_7"
    );
}

// ─── zero logarithm, unsupported domain, exact count ─────────────────────────

/// `ssa_solve` recovers every `k` in `1..p` on the fixture, each verified by `k·G = Q`.
#[test]
fn kat_ssa_solve_recovers_every_k() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    for k in 1..ANOMALOUS_TOY_P {
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k));
        assert_eq!(
            ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P),
            Ok(k),
            "ssa_solve must recover k = {k}"
        );
    }
}

/// `Q = ∞` is `0·G`: the logarithm is 0, on the anomalous fixture and on any other curve.
#[test]
fn kat_ssa_solve_infinity_is_zero() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    assert_eq!(
        ssa_solve::<FpNaive>(&curve, &g, &AffinePoint::Infinity, ANOMALOUS_TOY_P),
        Ok(0),
        "Q = ∞ must give k = 0 on the anomalous fixture"
    );

    // 0·G = ∞ holds in every group, so the answer does not depend on the curve being anomalous.
    let control = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(3u64),
        b: Uint::<4>::from(1u64),
        n: Uint::<4>::from(12u64),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(1u64),
    };
    assert_eq!(
        ssa_solve::<FpNaive>(&control, &point_on(0, 1, 7), &AffinePoint::Infinity, 12),
        Ok(0),
        "Q = ∞ must give k = 0 on a non-anomalous curve too"
    );
}

/// A field prime above 64 bits is outside the counting range: `Unsupported`, never a false
/// `NotAnomalous` (nor `Ok(false)`).
///
/// The prime is `2^64 + 13`, the least prime above `2^64`. Only the range is under test, so the
/// curve parameters are arbitrary.
#[test]
fn kat_large_prime_is_unsupported_not_non_anomalous() {
    let mut p_words = [0u64; 4];
    p_words[0] = 13;
    p_words[1] = 1;
    let p = Uint::<4>::from_words(p_words);
    let curve = Curve {
        p,
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(5u64),
        n: p,
        gx: Uint::<4>::from(1u64),
        gy: Uint::<4>::from(1u64),
    };
    let pt = AffinePoint::Finite { x: FpNaive::from_u64(1, &p), y: FpNaive::from_u64(1, &p) };

    assert!(
        matches!(verify_anomalous::<FpNaive>(&curve), Err(SsaError::Unsupported(_))),
        "verify_anomalous must report Unsupported for p ≥ 2^64"
    );
    let result = ssa_solve::<FpNaive>(&curve, &pt, &pt, 0);
    assert!(
        matches!(result, Err(SsaError::Unsupported(_))),
        "ssa_solve must report Unsupported for p ≥ 2^64, got: {result:?}"
    );
}

/// `verify_anomalous` counts exactly where `x³` no longer fits in `u64`.
///
/// `y² = x³ + 2` over `F_2738941` has `4p = 1 + 3·1901²`, so its trace is 1 and `#E = p` (found
/// by exhaustive counting; `b = 1` is the non-anomalous twist). The prime exceeds `2^(64/3)`, so
/// `x³` overflows `u64` for the larger `x`; a wrapping count would report a false negative.
#[test]
fn kat_verify_anomalous_exact_where_cube_overflows_u64() {
    let p = 2_738_941u64;
    let curve_with = |b: u64| Curve {
        p: Uint::<4>::from(p),
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(b),
        n: Uint::<4>::from(p),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(0u64),
    };

    assert_eq!(
        verify_anomalous::<FpNaive>(&curve_with(2)),
        Ok(true),
        "y² = x³ + 2 over F_{p} has exactly p points"
    );
    assert_eq!(
        verify_anomalous::<FpNaive>(&curve_with(1)),
        Ok(false),
        "y² = x³ + 1 over F_{p} is not anomalous"
    );
}
