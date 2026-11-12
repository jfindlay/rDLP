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
//! - `verify_anomalous` returns `true` on the fixture curve.
//! - `verify_anomalous` returns `false` on a non-anomalous control curve.
//! - The lift of `G = (3, 2)` matches hand-computed Z_7 coordinates to precision k=4, and a second
//!   point `(5, 2)` also lifts onto the curve.
//! - The lifted point satisfies the curve equation mod `7^4` (lift-correctness check).
//! - A 2-torsion point (`y = 0`) on a test curve errors via `HenselError::NonSimpleRoot`; the
//!   identity errors as `LiftError::PointAtInfinity`; `LiftError` carries its cause as the
//!   `source`.
//! - Lifts at precision k=1 and k=2 agree with the hand-computed values.
//! - A test (`#[ignore]`) asserting the lift value recorded from PARI/GP; it does not call `gp`.
//!
//! # KAT coverage — SSA reduction
//!
//! - `ssa_solve` recovers k=3 (and k=1, k=6) from Q = k·G on the anomalous fixture.
//! - `solve_brent` (independent rho solver) also recovers k=3 — cross-check.
//! - `ssa_solve` returns `SsaError::NotAnomalous` on a non-anomalous curve.
//! - A test (`#[ignore]`) asserting the reduction's answer recorded from PARI/GP; it does not call
//!   `gp`.
//! - `ssa_solve` recovers every `k` in `1..p` on the fixture.
//! - `Q = ∞` returns `Ok(0)`, on the fixture and on a non-anomalous curve.
//! - A prime that does not fit 64 bits is `Unsupported` for both `verify_anomalous` and
//!   `ssa_solve`, never `NotAnomalous`.
//! - `verify_anomalous` counts exactly where `x³` overflows `u64` (an anomalous curve over a
//!   prime above `2^21.4`, and the same curve with a tampered coefficient).
//! - Every `SsaInputDefect` is reachable through `ssa_solve`, in the order the checks run:
//!   composite modulus, singular curve, `G = ∞`, `G` and `Q` off the curve, `n ≠ p`.
//!   `verify_anomalous` shares the first two. No `SsaError` variant has a `source`.
//! - On the lift of `anomalous_toy` with its stored `a = 0`, `p·G̃` lies in `E₂` (the formal
//!   parameter has valuation 2), so the quotient of logarithms carries no information.
//! - On each anomalous fixture (`anomalous_toy` and a curve with `a ≠ 0` found by a recorded,
//!   bounded search) and for every point, at most one class of shifts `s` mod `p` (`a ↦ a +
//!   p·s`) puts `p·G̃` in `E₂`; `ssa_solve` recovers every `k` on each.
//!
//! The invariant-failure and inconsistent-quotient outcomes are reached through unit seams in
//! `rho::ssa::formal_log` and `rho::ssa::reduce`, which hold the private helpers they exercise.
//!
//! # Note on the fixture's lift
//!
//! The curve `y² = x³ + 5` over `F_7` has `a = 0`, so its stored lift to Z_7 is degenerate:
//! `p·G̃ ∈ E₂` there, and the quotient of leading units is a ratio of second-order terms. The
//! reduction therefore lifts the curve as `a ↦ a + p·s` with `s ≢ 0 (mod p)`, on which `p·G̃`
//! has valuation exactly 1 and the raw quotient is `k`. The `rho::ssa::reduce` module
//! documentation gives the derivation and the source.
//!
//! # Hand-computed lift values
//!
//! For `G = (3, 2)` on `y² = x³ + 5` over `F_7`, precision k=4 (modulus = 7^4 = 2401):
//!
//! - `x̃ = 3` (lifted exactly).
//! - `c = x₀³ + b = 27 + 5 = 32` (the RHS of the curve equation at x=3).
//! - Hensel-solve `g(y) = y² − 32 = 0` starting from `y₀ = 2`:
//!   - Step 1 (mod 49): y = 2 − (4−32)/4 = 2 + 7 = 9.
//!   - Later steps lift `y` to modulus 343 and then 2401.
//!   - Final: `ỹ = 940` (verified: 940² mod 2401 = 32 ✓).
//!
//! Verification: `940² = 883600`, `883600 mod 2401 = 32`. ✓

use crypto_bigint::Uint;
use num_bigint::BigInt;
use rho::curve::{AffinePoint, Curve};
use rho::ecdlp::solve_brent;
use shared_field::{Fp, FpNaive4 as FpNaive};
use rho::ssa::lift::{LiftError, check_lift_on_curve, lift_point, lifted_a};
use rho::ssa::reduce::kernel_filtration_level;
use rho::ssa::{
    ANOMALOUS_TOY_P, SsaError, SsaInputDefect, anomalous_toy, ssa_solve, verify_anomalous,
};
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

/// The coefficient `a` of the verbatim lift of `curve` (shift 0): `a` as stored.
fn stored_a(curve: &Curve) -> BigInt {
    lifted_a(curve, 0)
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

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, &stored_a(&curve), k)
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

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, &stored_a(&curve), k)
        .expect("lift of G = (3, 2) should succeed");

    // Verify ỹ² ≡ x̃³ + ax̃ + b (mod 7^4).
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve, &stored_a(&curve))
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

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&pt, &curve, &stored_a(&curve), k)
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
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve, &stored_a(&curve))
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
/// return `LiftError::Hensel(HenselError::NonSimpleRoot)`.
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

    let a = stored_a(&curve_2torsion);
    let result = lift_point::<FpNaive>(&two_torsion, &curve_2torsion, &a, 4);

    assert!(
        matches!(result, Err(LiftError::Hensel(HenselError::NonSimpleRoot { .. }))),
        "lifting a 2-torsion point (y=0) must return LiftError::Hensel(NonSimpleRoot), \
         got: {result:?}"
    );
}

/// The identity has no affine coordinates: lifting it is `LiftError::PointAtInfinity`.
#[test]
fn kat_lift_infinity_errors() {
    let curve = anomalous_toy();
    let result = lift_point::<FpNaive>(&AffinePoint::Infinity, &curve, &stored_a(&curve), 4);
    assert_eq!(result.err(), Some(LiftError::PointAtInfinity));
}

/// `LiftError` wraps its cause as the `source`; the identity has none.
#[test]
fn kat_lift_error_source_chain() {
    use std::error::Error;

    let curve_2torsion = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(1u64),
        b: Uint::<4>::from(0u64),
        n: Uint::<4>::from(8u64),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(0u64),
    };
    let a = stored_a(&curve_2torsion);
    let hensel = lift_point::<FpNaive>(&point_on(0, 0, 7), &curve_2torsion, &a, 4)
        .expect_err("2-torsion point cannot be lifted");
    assert!(hensel.source().is_some(), "LiftError::Hensel must expose its cause");
    assert!(LiftError::PointAtInfinity.source().is_none());
    assert!(!hensel.to_string().is_empty());
}

// ─── KAT 7: lift at precision k=1 (trivial case) ─────────────────────────────

/// Lift at precision k=1 returns the F_p coordinates unchanged (mod p).
///
/// At k=1, the modulus is p itself, so the lifted coordinates are just the F_p values.
#[test]
fn kat_lift_precision_1() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, &stored_a(&curve), 1)
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

    let (x_tilde, y_tilde) = lift_point::<FpNaive>(&g, &curve, &stored_a(&curve), 2)
        .expect("lift at k=2 should succeed");

    assert_eq!(*x_tilde.residue(), bi(3), "x̃ at k=2 must be 3");
    assert_eq!(*y_tilde.residue(), bi(9), "ỹ at k=2 must be 9 (first Newton step)");

    // Verify on-curve mod 7^2 = 49.
    let on_curve = check_lift_on_curve(&x_tilde, &y_tilde, &curve, &stored_a(&curve))
        .expect("check_lift_on_curve at k=2 should not error");
    assert!(on_curve, "lifted (3, 2) at k=2 must satisfy the curve equation mod 49");
}

// ─── value recorded from PARI/GP ──────────────────────────────────────────────

/// Asserts a value recorded from PARI/GP for the Hensel lift of `y² = 32` over `Z_7`, precision 4.
///
/// The test does not call `gp`; the value was obtained once with:
/// ```text
/// gp> p = 7; k = 4; f = y^2 - 32; y0 = 2; liftall(polrootsmod(f, p^k))
/// ```
/// Expected: `y = 940` (the unique root of `y² − 32 = 0` mod `7^4` with `y ≡ 2 mod 7`).
///
/// It is `#[ignore]`d and needs no external tool when run with `cargo test -- --ignored`.
#[test]
#[ignore = "asserts a value recorded from PARI/GP; does not call gp"]
fn kat_lift_pari_cross_check() {
    // PARI/GP: liftall(polrootsmod(y^2 - 32, 7^4)) should include 940.
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    let (_, y_tilde) = lift_point::<FpNaive>(&g, &curve, &stored_a(&curve), 4)
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
/// # Note on the lift
///
/// `ssa_solve` lifts the curve with `a ↦ a + p·s` (not the stored `a = 0`, on which `p·G̃ ∈ E₂`),
/// takes the quotient of the log leading units, which is `k = 3` itself, and checks `k·G == Q`
/// before returning it.
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
/// `6·G = (3, 5)`. Since k = 6 ≡ −1 (mod 7), `Q̃ = −G̃` and `k_raw = 6` directly.
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

/// Asserts a value recorded from PARI/GP for the SSA reduction.
///
/// The test does not call `gp`; the value was obtained once with:
/// ```text
/// gp> E = ellinit([0, 5], 7); G = [3, 2]; Q = [6, 5]; elllog(E, Q, G)
/// ```
/// Expected: `3` (the discrete log of Q = 3·G on y² = x³ + 5 over F_7).
///
/// It is `#[ignore]`d and needs no external tool when run with `cargo test -- --ignored`.
#[test]
#[ignore = "asserts a value recorded from PARI/GP; does not call gp"]
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

// ─── input validation ─────────────────────────────────────────────────────────

/// A curve over `F_p` with the given coefficients and stated order.
fn curve_over(p: u64, a: u64, b: u64, n: u64) -> Curve {
    Curve {
        p: Uint::<4>::from(p),
        a: Uint::<4>::from(a),
        b: Uint::<4>::from(b),
        n: Uint::<4>::from(n),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(0u64),
    }
}

fn invalid(defect: SsaInputDefect) -> Result<u64, SsaError> {
    Err(SsaError::InvalidInput(defect))
}

/// A composite field modulus is `FieldNotPrime`, from both entry points and ahead of the
/// singularity check.
#[test]
fn kat_composite_modulus_is_field_not_prime() {
    let g = point_on(3, 2, 9);
    // 9 = 3²: y² = x³ + 5 is nonsingular mod 9 here, so only primality is wrong.
    let curve = curve_over(9, 0, 5, 9);
    assert_eq!(
        verify_anomalous::<FpNaive>(&curve),
        Err(SsaError::InvalidInput(SsaInputDefect::FieldNotPrime))
    );
    assert_eq!(ssa_solve::<FpNaive>(&curve, &g, &g, 9), invalid(SsaInputDefect::FieldNotPrime));

    // A composite modulus that is also singular still reports the modulus first.
    let singular = curve_over(9, 0, 0, 9);
    assert_eq!(
        ssa_solve::<FpNaive>(&singular, &g, &g, 9),
        invalid(SsaInputDefect::FieldNotPrime)
    );

    // A large composite below 2^64 is caught by the primality test, not by the count.
    let big = u64::MAX; // 2^64 − 1 = 3 · 5 · 17 · 257 · 641 · 65537 · 6700417
    assert_eq!(
        verify_anomalous::<FpNaive>(&curve_over(big, 0, 5, big)),
        Err(SsaError::InvalidInput(SsaInputDefect::FieldNotPrime))
    );
}

/// A singular cubic is `SingularCurve`, from both entry points, never an anomalous count.
///
/// `y² = x³ + x + 5 = (x − 3)²(x + 6)` over `F_7` has a split node at `x = 3` (the tangent slopes
/// satisfy `s² = 3·3 ≡ 2`, a square), so it has exactly `p` points: a count alone would call it
/// anomalous. The other cases are singular by `4a³ + 27b² ≡ 0`.
#[test]
fn kat_singular_curve_is_rejected_before_the_count() {
    let p = 7u64;
    let g = point_on(3, 2, p);

    // The split-node cubic really has p points (brute-force count, independent of the solver).
    let points = 1 + (0..p)
        .map(|x| {
            let rhs = (x * x * x + x + 5) % p;
            (0..p).filter(|y| y * y % p == rhs).count() as u64
        })
        .sum::<u64>();
    assert_eq!(points, p, "the split-node cubic counts exactly p");

    for (a, b) in [(1u64, 5u64), (0, 0), (4, 2)] {
        let curve = curve_over(p, a, b, p);
        assert_eq!(
            verify_anomalous::<FpNaive>(&curve),
            Err(SsaError::InvalidInput(SsaInputDefect::SingularCurve)),
            "a = {a}, b = {b}"
        );
        assert_eq!(
            ssa_solve::<FpNaive>(&curve, &g, &g, p),
            invalid(SsaInputDefect::SingularCurve),
            "a = {a}, b = {b}"
        );
    }
}

/// `G = ∞` is `GeneratorAtInfinity`, even when `Q = ∞` too.
#[test]
fn kat_generator_at_infinity() {
    let curve = anomalous_toy();
    let inf = AffinePoint::<FpNaive>::Infinity;
    for q in [toy_point(6, 5), AffinePoint::Infinity] {
        assert_eq!(
            ssa_solve::<FpNaive>(&curve, &inf, &q, ANOMALOUS_TOY_P),
            invalid(SsaInputDefect::GeneratorAtInfinity)
        );
    }
}

/// `G` off the curve is `GeneratorNotOnCurve`, ahead of any check on `Q`.
#[test]
fn kat_generator_off_curve() {
    let curve = anomalous_toy();
    // (3, 3): 9 ≢ 3³ + 5 = 4 mod 7.
    let g = toy_point(3, 3);
    for q in [toy_point(6, 5), toy_point(6, 6), AffinePoint::Infinity] {
        assert_eq!(
            ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P),
            invalid(SsaInputDefect::GeneratorNotOnCurve)
        );
    }
}

/// `Q` off the curve is `TargetNotOnCurve`, ahead of the count (so on a non-anomalous curve too).
#[test]
fn kat_target_off_curve() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    // (6, 6): 36 ≡ 1 ≢ 6³ + 5 = 221 ≡ 4 mod 7.
    assert_eq!(
        ssa_solve::<FpNaive>(&curve, &g, &toy_point(6, 6), ANOMALOUS_TOY_P),
        invalid(SsaInputDefect::TargetNotOnCurve)
    );

    let control = curve_over(7, 3, 1, 12);
    assert_eq!(
        ssa_solve::<FpNaive>(&control, &point_on(0, 1, 7), &point_on(0, 2, 7), 12),
        invalid(SsaInputDefect::TargetNotOnCurve)
    );
}

/// On an anomalous curve any `n ≠ p` is `OrderNotFieldPrime`; `Q = ∞` and a non-anomalous curve
/// are decided before `n` is looked at.
#[test]
fn kat_order_must_equal_field_prime() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);
    let q = toy_point(6, 5);
    for n in [0, 1, 5, 12, 14] {
        assert_eq!(
            ssa_solve::<FpNaive>(&curve, &g, &q, n),
            invalid(SsaInputDefect::OrderNotFieldPrime),
            "n = {n}"
        );
    }
    assert_eq!(ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P), Ok(3));

    // Q = ∞ is answered before n is checked.
    assert_eq!(ssa_solve::<FpNaive>(&curve, &g, &AffinePoint::Infinity, 5), Ok(0));
    // A non-anomalous curve is a proven negative whatever n says.
    let control = curve_over(7, 3, 1, 12);
    assert_eq!(
        ssa_solve::<FpNaive>(&control, &point_on(0, 1, 7), &point_on(0, 1, 7), 7),
        Err(SsaError::NotAnomalous)
    );
}

/// No `SsaError` variant wraps another error, so none has a source; each displays something.
#[test]
fn kat_ssa_error_has_no_source() {
    use std::error::Error;

    let all = [
        SsaError::InvalidInput(SsaInputDefect::FieldNotPrime),
        SsaError::InvalidInput(SsaInputDefect::SingularCurve),
        SsaError::InvalidInput(SsaInputDefect::GeneratorAtInfinity),
        SsaError::InvalidInput(SsaInputDefect::GeneratorNotOnCurve),
        SsaError::InvalidInput(SsaInputDefect::TargetNotOnCurve),
        SsaError::InvalidInput(SsaInputDefect::OrderNotFieldPrime),
        SsaError::NotAnomalous,
        SsaError::Unsupported("x"),
        SsaError::Inconsistent { k: 3 },
        SsaError::InvariantFailure("x"),
    ];
    for e in all {
        assert!(e.source().is_none(), "{e:?}");
        assert!(!e.to_string().is_empty(), "{e:?}");
    }
}

// ─── the lift of the curve: where p·G̃ lands in the filtration ─────────────────

/// The valuation of the formal parameter of `p·G̃` on the lift of the toy curve with its stored
/// `a = 0`: it lies in `E₂` (valuation 2), so the quotient of logarithms there is a ratio of
/// second-order terms and not `k`. Every shift `s ≢ 0 (mod 7)` puts it in `E₁ \ E₂`.
#[test]
fn kat_stored_lift_of_toy_curve_is_degenerate() {
    let curve = anomalous_toy();
    let g = toy_point(3, 2);

    assert_eq!(kernel_filtration_level::<FpNaive>(&curve, &g, &stored_a(&curve)), Ok(2));
    for shift in 1..ANOMALOUS_TOY_P {
        assert_eq!(
            kernel_filtration_level::<FpNaive>(&curve, &g, &lifted_a(&curve, shift)),
            Ok(1),
            "shift s = {shift}"
        );
    }
}

/// The shifts `s` in `0..p` for which `p·P̃ ∈ E₂` on the lift `a ↦ a + p·s` (shift 0 is the
/// stored coefficient).
fn degenerate_shifts(curve: &Curve, point: &AffinePoint<FpNaive>, p: u64) -> Vec<u64> {
    (0..p)
        .filter(|&s| {
            let level = kernel_filtration_level::<FpNaive>(curve, point, &lifted_a(curve, s))
                .expect("the lift of a point of an anomalous curve reaches the kernel");
            level >= 2
        })
        .collect()
}

/// Every finite point of the curve over `F_p`.
fn finite_points(p: u64, a: u64, b: u64) -> Vec<AffinePoint<FpNaive>> {
    let mut points = Vec::new();
    for x in 0..p {
        let rhs = (x * x % p * x + a * x + b) % p;
        points.extend((0..p).filter(|y| y * y % p == rhs).map(|y| point_on(x, y, p)));
    }
    points
}

/// At most one class of shifts `s` mod `p` puts `p·P̃` in `E₂`, for every point of each fixture,
/// and it is the same class for every point (the degeneracy belongs to the lifted curve).
///
/// The classes are `s = 0` (the stored coefficient) and `s = 1, …, p − 1`. Two degenerate
/// classes would make one retry with a second shift insufficient; none is observed. On both
/// fixtures the one degenerate class is the stored coefficient, `s = 0`.
#[test]
fn kat_at_most_one_shift_class_is_degenerate() {
    let second = second_fixture();
    let fixtures = [
        (anomalous_toy(), ANOMALOUS_TOY_P, 0u64, 5u64),
        (second.curve(), second.p, second.a, second.b),
    ];

    for (curve, p, a, b) in fixtures {
        let points = finite_points(p, a, b);
        assert_eq!(points.len() as u64, p - 1, "p − 1 finite points on an anomalous curve");
        let classes: Vec<Vec<u64>> =
            points.iter().map(|pt| degenerate_shifts(&curve, pt, p)).collect();

        for (pt, shifts) in points.iter().zip(&classes) {
            assert!(shifts.len() <= 1, "p = {p}, point {pt:?}: degenerate shifts {shifts:?}");
        }
        assert!(
            classes.windows(2).all(|w| w[0] == w[1]),
            "p = {p}: the degenerate class differs between points: {classes:?}"
        );
        assert_eq!(classes[0], vec![0], "p = {p}: the stored coefficient is the degenerate lift");
    }
}

// ─── a second anomalous fixture, with a ≠ 0 ───────────────────────────────────

/// An anomalous curve found by search, with its base point.
struct SearchedFixture {
    p: u64,
    a: u64,
    b: u64,
    gx: u64,
    gy: u64,
}

impl SearchedFixture {
    fn curve(&self) -> Curve {
        Curve {
            p: Uint::<4>::from(self.p),
            a: Uint::<4>::from(self.a),
            b: Uint::<4>::from(self.b),
            n: Uint::<4>::from(self.p),
            gx: Uint::<4>::from(self.gx),
            gy: Uint::<4>::from(self.gy),
        }
    }

    fn generator(&self) -> AffinePoint<FpNaive> {
        point_on(self.gx, self.gy, self.p)
    }
}

/// The number of points of `y² = x³ + ax + b` over `F_p`, by direct enumeration.
fn naive_point_count(p: u64, a: u64, b: u64) -> u64 {
    1 + (0..p)
        .map(|x| {
            let rhs = (x * x % p * x + a * x + b) % p;
            (0..p).filter(|y| y * y % p == rhs).count() as u64
        })
        .sum::<u64>()
}

/// Search for an anomalous curve with `a ≠ 0` by enumeration.
///
/// The search runs over primes `11 ≤ p ≤ 31` in increasing order and, for each, over `a` in
/// `1..p` and `b` in `1..p` in lexicographic order, keeping the first nonsingular curve whose
/// naive point count equals `p`. The base point is the finite point with the least `x`, then the
/// least `y`. The search is bounded by `p ≤ 31` and stops at the first hit.
fn search_second_fixture() -> Option<SearchedFixture> {
    const LOWEST_PRIME: u64 = 11;
    const HIGHEST_PRIME: u64 = 31;
    let is_prime = |n: u64| n > 1 && (2..n).all(|d| n % d != 0);

    for p in (LOWEST_PRIME..=HIGHEST_PRIME).filter(|&n| is_prime(n)) {
        for a in 1..p {
            for b in 1..p {
                let disc = (4 * (a * a % p * a % p) + 27 * (b * b % p)) % p;
                if disc == 0 || naive_point_count(p, a, b) != p {
                    continue;
                }
                let (gx, gy) = (0..p)
                    .flat_map(|x| (1..p).map(move |y| (x, y)))
                    .find(|&(x, y)| y * y % p == (x * x % p * x + a * x + b) % p)
                    .expect("an anomalous curve has finite points");
                return Some(SearchedFixture { p, a, b, gx, gy });
            }
        }
    }
    None
}

/// The second anomalous fixture: `y² = x³ + x + 5` over `F_11`, `G = (0, 4)`.
///
/// Recorded search: [`search_second_fixture`] over primes `11 ≤ p ≤ 31`, `a` in `1..p`, `b` in
/// `1..p`, first hit `p = 11, a = 1, b = 5`. The search is re-run by
/// `kat_second_fixture_matches_its_recorded_search`, so the constants below are not trusted.
fn second_fixture() -> SearchedFixture {
    SearchedFixture { p: 11, a: 1, b: 5, gx: 0, gy: 4 }
}

/// The recorded search finds the pinned curve, and `verify_anomalous` re-proves that it is
/// anomalous (an exact count, independent of the enumeration above) while it is not singular.
#[test]
fn kat_second_fixture_matches_its_recorded_search() {
    let found = search_second_fixture().expect("the search bound contains an anomalous curve");
    let pinned = second_fixture();
    assert_eq!(
        (found.p, found.a, found.b, found.gx, found.gy),
        (pinned.p, pinned.a, pinned.b, pinned.gx, pinned.gy)
    );
    assert_ne!(pinned.a % pinned.p, 0, "the fixture must have a ≢ 0");
    assert_eq!(verify_anomalous::<FpNaive>(&pinned.curve()), Ok(true));
}

/// `ssa_solve` recovers every `k` in `1..p` on the second fixture, each verified by `k·G = Q`.
#[test]
fn kat_second_fixture_recovers_every_k() {
    let fixture = second_fixture();
    let (curve, g) = (fixture.curve(), fixture.generator());

    for k in 1..fixture.p {
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k));
        assert_eq!(ssa_solve::<FpNaive>(&curve, &g, &q, fixture.p), Ok(k), "k = {k}");
    }
}
