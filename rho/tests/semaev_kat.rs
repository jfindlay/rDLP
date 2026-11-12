//! Known-answer tests (KATs) for the Semaev summation polynomial construction.
//!
//! # Fixture
//!
//! All tests use the Semaev toy curve `y² = x³ + x + 33` over `F_47` (`p = 47`, `a = 1`,
//! `b = 33`). The generator `G = (10, 3)` satisfies the curve equation:
//! `10³ + 10 + 33 = 1043 = 9 = 3² mod 47 ✓`. The group order is `n = 60 = 2² · 3 · 5`
//! (verified by `n·G = ∞`).
//!
//! The `F_p[x]` resultant KATs use `F_47` arithmetic (the same prime as the toy fixture).
//! The multivariate polynomial KATs use `F_47` arithmetic.
//!
//! Known multiples of `G` (computed from the group law):
//! - `G  = (10, 3)`,  `−G  = (10, 44)`
//! - `2G = (7, 30)`,  `−2G = (7, 17)`
//! - `3G = (17, 13)`, `−3G = (17, 34)`
//! - `4G = (23, 12)`, `−4G = (23, 35)`
//! - `5G = (32, 36)`, `−5G = (32, 11)`
//! - `6G = (19, 7)`,  `−6G = (19, 40)`
//! - `8G = (25, 28)`, `−8G = (25, 19)`
//!
//! # KAT coverage — `F_p[x]` resultant
//!
//! 1. **Resultant zero-iff-common-factor**: `Res(f, g) = 0 ⟺ gcd(f, g) ≠ 1`.
//!    - `Res(x-2, x-3) ≠ 0` (coprime linear factors).
//!    - `Res((x-2)(x-3), x-2) = 0` (share root x=2).
//!    - `Res((x-2)(x-3), (x-3)(x-5)) = 0` (share root x=3).
//! 2. **Resultant symmetry**: `Res(f, g) = (-1)^(deg_f * deg_g) * Res(g, f)`.
//!    - Verified for linear × linear (sign = -1).
//!    - Verified for linear × quadratic (even product → sign = +1).
//!    - Verified for quadratic × quadratic (even product → sign = +1).
//! 3. **Resultant of coprime quadratics is nonzero**: `Res(x²+1, x²+2) ≠ 0` over `F_47`
//!    (both irreducible over `F_47` — verified: -1 and -2 are QNR mod 47 since 47 ≡ 3 mod 4).
//! 4. **Resultant of zero polynomial is zero**: `Res(0, g) = 0`.
//! 5. **Resultant value spot-check**: `Res(x-2, x-3) = 1` over `F_47`.
//!
//! # KAT coverage — multivariate/symmetric round-trip
//!
//! 6. **Symmetric polynomial round-trip**: a symmetric polynomial round-trips through the
//!    `MultiPoly` representation (build → evaluate → same value).
//! 7. **Evaluation at permutation is invariant**: for a symmetric polynomial `f`,
//!    `f(a, b, c) = f(b, a, c) = f(c, b, a)` etc.
//! 8. **Partial evaluation**: fixing one variable of a 3-variable polynomial gives a
//!    2-variable polynomial with the correct evaluation.
//! 9. **`is_symmetric` correctly identifies symmetric polynomials**.
//! 10. **`symmetrize` produces a symmetric polynomial**.
//!
//! # KAT coverage — `S_2`, `S_3` base cases + vanishing relation
//!
//! 11. **`S_2` vanishing**: `S_2(x_1, x_2) = 0 ⟺ P_2 = −P_1`.
//!     - `S_2(10, 10) = 0` and `G + (−G) = ∞` ✓.
//!     - `S_2(10, 7) ≠ 0` and `G + 2G ≠ ∞` ✓.
//! 12. **`S_3` vanishing**: `S_3(x_1, x_2, x_3) = 0` for triples with `P_1+P_2+P_3 = ∞`.
//!     - `S_3(10, 7, 17) = 0`: x-coords of `G, 2G, −3G` (sum `∞`).
//!     - `S_3(10, 7, 23) ≠ 0`: x-coords of `G, 2G, 4G` — no y-values make sum `∞`.
//!     - `S_3(10, 10, 10) ≠ 0`: no y-values for x=10 make sum `∞`.
//!     - `S_3(10, 7, 32) ≠ 0`: x-coords of `G, 2G, 5G` — no y-values make sum `∞`.
//! 13. **`S_3` symmetry**: `S_3` is invariant under all permutations of its arguments.
//! 14. **`S_3` degree**: degree ≤ 2 in each variable.
//!
//! # KAT coverage — `S_m` resultant recursion
//!
//! 15. **Recursion base cases**: `semaev_poly(2)` = `s2`, `semaev_poly(3)` = `s3`.
//! 16. **`S_4` structure**: 4 variables, non-zero, symmetric.
//! 17. **`S_4` vanishing** (primary correctness signal):
//!     - `S_4(10, 10, 7, 7) = 0`: x-coords of `G, −G, 2G, −2G` (sum `∞`).
//!     - `S_4(10, 7, 17, 19) = 0`: x-coords of `G, 2G, 3G, −6G` (sum `∞`).
//!     - `S_4(10, 7, 23, 25) ≠ 0`: x-coords of `G, 2G, 4G, 8G` — no y-values make sum `∞`.
//! 18. **`S_4` degree growth**: degree ≤ 4 in each variable (expected `2^(m-2) = 4` for `m=4`).
//! 19. **`S_4` symmetry**: invariant under all permutations of its 4 arguments.
//! 20. **`S_4` vanishing predicate agreement**: polynomial zero ⟺ group-law existence.
//! 21. **PARI/GP-derived values** (`#[ignore]`): `S_4` vanishes at points derived with PARI/GP;
//!     does not call `gp`.
//!
//! # Toy-scale boundary
//!
//! The fixture is toy-scale (`p = 47`, group order `n = 60`). The algorithms are
//! crypto-scale-correct; only the parameters are small for auditability.

use crypto_bigint::Uint;
use shared_field::{Fp, FpNaive4 as FpNaive};
use rho::semaev::poly::{FpPoly, MultiPoly, resultant};
use rho::semaev::base::{s2 as base_s2, s3, vanishes_s2, vanishes_s3};
use rho::semaev::{SEMAEV_TOY_P, SemaevError, semaev_poly, semaev_toy};
use rho::curve::{AffinePoint, JacobianPoint};

// ─── helpers ──────────────────────────────────────────────────────────────────

/// The prime for the Semaev toy fixture: `p = 47`.
const P: u64 = SEMAEV_TOY_P; // 47

/// Return the `Uint<4>` modulus for the toy fixture.
fn p() -> Uint<4> {
    Uint::<4>::from(P)
}

/// Construct an `FpNaive` element from a `u64`.
fn fp(v: u64) -> FpNaive {
    FpNaive::from_u64(v, &p())
}

/// Construct a linear polynomial `x - root` over `F_47`.
///
/// `x - root = x + (p - root)` since `-root mod p = p - root`.
fn linear(root: u64) -> FpPoly<FpNaive> {
    let p_val = p();
    FpPoly::from_coeffs(
        vec![
            FpNaive::from_u64((P - root % P) % P, &p_val), // constant: -root mod p
            FpNaive::from_u64(1, &p_val),                   // leading: 1
        ],
        &p_val,
    )
}

/// Construct the product `(x - r1)(x - r2)` over `F_47`.
fn quadratic_with_roots(r1: u64, r2: u64) -> FpPoly<FpNaive> {
    let p_val = p();
    linear(r1).mul(&linear(r2), &p_val)
}

// ─── fixture verification ─────────────────────────────────────────────────────

/// Verify the toy curve fixture: generator on curve and n·G = ∞.
#[test]
fn semaev_toy_fixture_valid() {
    let c = semaev_toy();
    let g: AffinePoint<FpNaive> = c.generator();
    assert!(c.is_on_curve(&g), "semaev toy: generator not on curve");
    let ng = c.scalar_mul(&g, &c.n);
    assert!(ng.is_infinity(), "semaev toy: n·G should be ∞");
}

// ─── KAT 1: resultant zero-iff-common-factor ─────────────────────────────────

/// `Res(x-2, x-3) ≠ 0` — coprime linear factors over `F_47`.
#[test]
fn resultant_coprime_linears_nonzero() {
    let p_val = p();
    let f = linear(2);
    let g = linear(3);
    let res = resultant(&f, &g, &p_val).unwrap();
    assert!(
        !res.is_zero(&p_val),
        "Res(x-2, x-3) should be nonzero (coprime)"
    );
}

/// `Res((x-2)(x-3), x-2) = 0` — share root x=2.
#[test]
fn resultant_common_root_is_zero() {
    let p_val = p();
    let f = quadratic_with_roots(2, 3); // (x-2)(x-3)
    let g = linear(2);                  // x-2
    let res = resultant(&f, &g, &p_val).unwrap();
    assert!(
        res.is_zero(&p_val),
        "Res((x-2)(x-3), x-2) should be 0 (share root x=2)"
    );
}

/// `Res((x-2)(x-3), (x-3)(x-5)) = 0` — share root x=3.
#[test]
fn resultant_shared_root_quadratics_is_zero() {
    let p_val = p();
    let f = quadratic_with_roots(2, 3); // (x-2)(x-3)
    let g = quadratic_with_roots(3, 5); // (x-3)(x-5)
    let res = resultant(&f, &g, &p_val).unwrap();
    assert!(
        res.is_zero(&p_val),
        "Res((x-2)(x-3), (x-3)(x-5)) should be 0 (share root x=3)"
    );
}

/// `Res(0, g) = 0` — zero polynomial.
#[test]
fn resultant_zero_poly_is_zero() {
    let p_val = p();
    let f: FpPoly<FpNaive> = FpPoly::zero();
    let g = linear(5);
    let res = resultant(&f, &g, &p_val).unwrap();
    assert!(res.is_zero(&p_val), "Res(0, g) should be 0");
}

// ─── KAT 2: resultant symmetry ────────────────────────────────────────────────

/// `Res(f, g) = (-1)^(deg_f * deg_g) * Res(g, f)` — linear × linear (sign = -1).
///
/// For `f = x - 2`, `g = x - 3`: `deg_f * deg_g = 1 * 1 = 1` (odd) → `Res(f,g) = -Res(g,f)`.
/// `Res(x-a, x-b) = b - a` (the resultant of two monic linears is the difference of roots).
/// So `Res(x-2, x-3) = 3 - 2 = 1` and `Res(x-3, x-2) = 2 - 3 = -1 = 46 mod 47`.
#[test]
fn resultant_symmetry_linear_linear() {
    let p_val = p();
    let f = linear(2); // x - 2
    let g = linear(3); // x - 3
    let res_fg = resultant(&f, &g, &p_val).unwrap();
    let res_gf = resultant(&g, &f, &p_val).unwrap();
    // deg_f * deg_g = 1 * 1 = 1 (odd) → sign = -1 → Res(f,g) = -Res(g,f)
    let neg_res_gf = res_gf.neg(&p_val);
    assert_eq!(
        res_fg, neg_res_gf,
        "Res(f,g) should equal -Res(g,f) for linear × linear"
    );
}

/// `Res(f, g) = (-1)^(deg_f * deg_g) * Res(g, f)` — quadratic × quadratic (sign = +1).
///
/// For `f = (x-2)(x-3)`, `g = (x-5)(x-7)` (coprime): `deg_f * deg_g = 4` (even) →
/// `Res(f,g) = Res(g,f)`.
#[test]
fn resultant_symmetry_quadratic_quadratic() {
    let p_val = p();
    let f = quadratic_with_roots(2, 3); // (x-2)(x-3)
    let g = quadratic_with_roots(5, 7); // (x-5)(x-7)
    let res_fg = resultant(&f, &g, &p_val).unwrap();
    let res_gf = resultant(&g, &f, &p_val).unwrap();
    // deg_f * deg_g = 2 * 2 = 4 (even) → sign = +1 → Res(f,g) = Res(g,f)
    assert_eq!(
        res_fg, res_gf,
        "Res(f,g) should equal Res(g,f) for quadratic × quadratic"
    );
}

/// `Res(f, g) = (-1)^(deg_f * deg_g) * Res(g, f)` — linear × quadratic (even product).
///
/// `deg_f * deg_g = 1 * 2 = 2` (even) → `Res(f,g) = Res(g,f)`.
#[test]
fn resultant_symmetry_linear_quadratic() {
    let p_val = p();
    let f = linear(2);                  // x - 2
    let g = quadratic_with_roots(5, 7); // (x-5)(x-7)
    let res_fg = resultant(&f, &g, &p_val).unwrap();
    let res_gf = resultant(&g, &f, &p_val).unwrap();
    // deg_f * deg_g = 1 * 2 = 2 (even) → sign = +1 → Res(f,g) = Res(g,f)
    assert_eq!(
        res_fg, res_gf,
        "Res(f,g) should equal Res(g,f) for linear × quadratic (even product)"
    );
}

// ─── KAT 3: resultant of coprime quadratics is nonzero ───────────────────────

/// `Res(x²+1, x²+2) ≠ 0` over `F_47` — both irreducible (no roots mod 47).
///
/// Since `47 ≡ 3 mod 4`, `-1` is a QNR mod 47 (Euler's criterion: `(-1)^23 = -1 mod 47`).
/// Similarly, `-2` is a QNR mod 47 (verified: `(-2)^23 = -(2^23) = -1 mod 47` since
/// `2^23 = 1 mod 47` — the order of 2 mod 47 divides 23, and indeed `2^23 mod 47 = 1`).
/// Both polynomials are irreducible over `F_47`, hence coprime.
#[test]
fn resultant_coprime_irreducible_quadratics_nonzero() {
    let p_val = p();
    // x^2 + 1 — irreducible over F_47 since -1 is QNR mod 47
    let f = FpPoly::from_coeffs(
        vec![fp(1), fp(0), fp(1)], // 1 + 0*x + 1*x^2
        &p_val,
    );
    // x^2 + 2 — irreducible over F_47 since -2 is QNR mod 47
    let g = FpPoly::from_coeffs(
        vec![fp(2), fp(0), fp(1)], // 2 + 0*x + 1*x^2
        &p_val,
    );
    let res = resultant(&f, &g, &p_val).unwrap();
    assert!(
        !res.is_zero(&p_val),
        "Res(x²+1, x²+2) should be nonzero (coprime irreducible quadratics over F_47)"
    );
}

// ─── KAT 4: resultant value spot-check ───────────────────────────────────────

/// `Res(x-2, x-3) = -1 = 46 mod 47` over `F_47`.
///
/// The resultant of two monic linears `x - a` and `x - b` is `g(a) = a - b` (the value of
/// `g = x - b` at the root `a` of `f = x - a`). Here `Res(x-2, x-3) = g(2) = 2 - 3 = -1 = 46`.
///
/// Equivalently: `Res(f, g) = ∏_{f(α)=0} g(α) = g(2) = 2 - 3 = -1 mod 47`.
#[test]
fn resultant_linear_value_spot_check() {
    let p_val = p();
    let f = linear(2); // x - 2
    let g = linear(3); // x - 3
    let res = resultant(&f, &g, &p_val).unwrap();
    // Res(x-2, x-3) = g(2) = 2 - 3 = -1 = 46 mod 47
    assert_eq!(
        res.to_uint(),
        Uint::<4>::from(46u64),
        "Res(x-2, x-3) should be -1 = 46 mod 47"
    );
}

// ─── KAT 5: multivariate symmetric round-trip ────────────────────────────────

/// A symmetric polynomial round-trips through the `MultiPoly` representation.
///
/// Build `f(x_0, x_1) = x_0^2 + x_1^2 + 2*x_0*x_1 = (x_0 + x_1)^2`, evaluate at
/// `(3, 5)` and `(5, 3)` — both should give `(3+5)^2 = 64 = 64 - 47 = 17 mod 47`.
#[test]
fn multi_poly_symmetric_round_trip() {
    // f(x_0, x_1) = x_0^2 + x_1^2 + 2*x_0*x_1 = (x_0 + x_1)^2
    let mut f = MultiPoly::zero(2, P);
    f.add_term(vec![2, 0], 1); // x_0^2
    f.add_term(vec![0, 2], 1); // x_1^2
    f.add_term(vec![1, 1], 2); // 2*x_0*x_1

    assert!(f.is_symmetric(), "f should be symmetric");

    // (3 + 5)^2 = 64 = 64 - 47 = 17 mod 47
    let v1 = f.eval(&[3, 5]).unwrap();
    let v2 = f.eval(&[5, 3]).unwrap();
    assert_eq!(v1, 17, "f(3,5) = (3+5)^2 mod 47 = 17");
    assert_eq!(v2, 17, "f(5,3) = (5+3)^2 mod 47 = 17");
    assert_eq!(v1, v2, "symmetric polynomial: f(3,5) = f(5,3)");
}

// ─── KAT 6: evaluation at permutation is invariant ───────────────────────────

/// For a symmetric 3-variable polynomial, evaluation at any permutation of arguments
/// gives the same result.
///
/// Build `f(x_0, x_1, x_2) = x_0 + x_1 + x_2` (elementary symmetric polynomial e_1).
/// Evaluate at all permutations of `(2, 5, 11)` — all should give `2+5+11 = 18`.
#[test]
fn multi_poly_symmetric_eval_permutation_invariant_3var() {
    // f = x_0 + x_1 + x_2 (symmetric)
    let mut f = MultiPoly::zero(3, P);
    f.add_term(vec![1, 0, 0], 1); // x_0
    f.add_term(vec![0, 1, 0], 1); // x_1
    f.add_term(vec![0, 0, 1], 1); // x_2

    assert!(f.is_symmetric(), "x_0 + x_1 + x_2 should be symmetric");

    // 2 + 5 + 11 = 18 mod 47
    let expected = 18u64;
    let perms = [
        [2u64, 5, 11],
        [2, 11, 5],
        [5, 2, 11],
        [5, 11, 2],
        [11, 2, 5],
        [11, 5, 2],
    ];
    for perm in &perms {
        let v = f.eval(perm).unwrap();
        assert_eq!(
            v, expected,
            "symmetric f: eval at {:?} should be {expected}",
            perm
        );
    }
}

/// For a symmetric 3-variable polynomial `x_0*x_1 + x_1*x_2 + x_0*x_2` (e_2),
/// evaluation at all permutations of `(2, 3, 5)` gives the same result.
///
/// `e_2(2,3,5) = 2*3 + 3*5 + 2*5 = 6 + 15 + 10 = 31 mod 47`.
#[test]
fn multi_poly_e2_symmetric_eval_permutation_invariant() {
    // e_2(x_0, x_1, x_2) = x_0*x_1 + x_1*x_2 + x_0*x_2
    let mut f = MultiPoly::zero(3, P);
    f.add_term(vec![1, 1, 0], 1); // x_0*x_1
    f.add_term(vec![0, 1, 1], 1); // x_1*x_2
    f.add_term(vec![1, 0, 1], 1); // x_0*x_2

    assert!(f.is_symmetric(), "e_2 should be symmetric");

    // e_2(2,3,5) = 2*3 + 3*5 + 2*5 = 6 + 15 + 10 = 31 mod 47
    let expected = f.eval(&[2, 3, 5]).unwrap();
    assert_eq!(expected, 31, "e_2(2,3,5) = 31 mod 47");
    let perms = [
        [2u64, 3, 5],
        [2, 5, 3],
        [3, 2, 5],
        [3, 5, 2],
        [5, 2, 3],
        [5, 3, 2],
    ];
    for perm in &perms {
        let v = f.eval(perm).unwrap();
        assert_eq!(
            v, expected,
            "e_2: eval at {:?} should equal eval at (2,3,5)",
            perm
        );
    }
}

// ─── KAT 7: partial evaluation ────────────────────────────────────────────────

/// Fixing one variable of a 3-variable polynomial gives a 2-variable polynomial.
///
/// `f(x_0, x_1, x_2) = x_0 + x_1 + x_2`. Fix `x_0 = 7` → `g(x_1, x_2) = 7 + x_1 + x_2`.
/// Evaluate `g(3, 5)` → `7 + 3 + 5 = 15`.
#[test]
fn multi_poly_partial_eval_3var() {
    let mut f = MultiPoly::zero(3, P);
    f.add_term(vec![1, 0, 0], 1); // x_0
    f.add_term(vec![0, 1, 0], 1); // x_1
    f.add_term(vec![0, 0, 1], 1); // x_2

    // Fix x_0 = 7
    let g = f.partial_eval(&[Some(7), None, None]).unwrap();
    assert_eq!(g.num_vars, 2, "one free variable after fixing x_0");

    // g(3, 5) = 7 + 3 + 5 = 15
    let result = g.eval(&[3, 5]).unwrap();
    assert_eq!(result, 15, "partial eval: 7 + 3 + 5 = 15");
}

/// Fixing two variables of a 3-variable polynomial gives a 1-variable polynomial.
///
/// `f(x_0, x_1, x_2) = x_0*x_1 + x_1*x_2 + x_0*x_2`.
/// Fix `x_0 = 2`, `x_1 = 3` → `g(x_2) = 2*3 + 3*x_2 + 2*x_2 = 6 + 5*x_2`.
/// `g(5) = 6 + 25 = 31 mod 47`. `g(1) = 6 + 5 = 11`.
#[test]
fn multi_poly_partial_eval_two_fixed() {
    // f(x_0, x_1, x_2) = x_0*x_1 + x_1*x_2 + x_0*x_2
    let mut f = MultiPoly::zero(3, P);
    f.add_term(vec![1, 1, 0], 1); // x_0*x_1
    f.add_term(vec![0, 1, 1], 1); // x_1*x_2
    f.add_term(vec![1, 0, 1], 1); // x_0*x_2

    // Fix x_0 = 2, x_1 = 3 → g(x_2) = 2*3 + 3*x_2 + 2*x_2 = 6 + 5*x_2
    let g = f.partial_eval(&[Some(2), Some(3), None]).unwrap();
    assert_eq!(g.num_vars, 1, "one free variable after fixing x_0 and x_1");

    // g(5) = 6 + 5*5 = 6 + 25 = 31 mod 47
    let result = g.eval(&[5]).unwrap();
    assert_eq!(result, 31, "partial eval: 6 + 5*5 = 31 mod 47");

    // g(1) = 6 + 5*1 = 11
    let result2 = g.eval(&[1]).unwrap();
    assert_eq!(result2, 11, "partial eval: 6 + 5*1 = 11");
}

// ─── KAT 8: is_symmetric correctly identifies symmetric polynomials ───────────

/// `x_0^2 + x_1^2 + 2*x_0*x_1` is symmetric.
#[test]
fn is_symmetric_true_for_symmetric_poly() {
    let mut f = MultiPoly::zero(2, P);
    f.add_term(vec![2, 0], 1); // x_0^2
    f.add_term(vec![0, 2], 1); // x_1^2
    f.add_term(vec![1, 1], 2); // 2*x_0*x_1
    assert!(f.is_symmetric(), "x_0^2 + x_1^2 + 2*x_0*x_1 should be symmetric");
}

/// `x_0^2` is not symmetric in 2 variables.
#[test]
fn is_symmetric_false_for_asymmetric_poly() {
    let mut f = MultiPoly::zero(2, P);
    f.add_term(vec![2, 0], 1); // x_0^2 only
    assert!(!f.is_symmetric(), "x_0^2 should not be symmetric in 2 variables");
}

/// The zero polynomial is symmetric (vacuously).
#[test]
fn is_symmetric_zero_poly() {
    let f = MultiPoly::zero(3, P);
    assert!(f.is_symmetric(), "zero polynomial should be symmetric");
}

// ─── KAT 9: symmetrize produces a symmetric polynomial ───────────────────────

/// `symmetrize(x_0^2)` over 2 variables gives `x_0^2 + x_1^2` (symmetric).
#[test]
fn symmetrize_produces_symmetric_poly() {
    let mut f = MultiPoly::zero(2, P);
    f.add_term(vec![2, 0], 1); // x_0^2
    let sym = f.symmetrize();
    assert!(sym.is_symmetric(), "symmetrize(x_0^2) should be symmetric");
    // sym should contain x_0^2 and x_1^2
    assert!(
        sym.terms.contains_key(&vec![2u64, 0]),
        "symmetrized poly should contain x_0^2"
    );
    assert!(
        sym.terms.contains_key(&vec![0u64, 2]),
        "symmetrized poly should contain x_1^2"
    );
}

/// `symmetrize(x_0 * x_1^2)` over 3 variables gives a symmetric polynomial.
#[test]
fn symmetrize_3var_produces_symmetric_poly() {
    let mut f = MultiPoly::zero(3, P);
    f.add_term(vec![1, 2, 0], 1); // x_0 * x_1^2
    let sym = f.symmetrize();
    assert!(sym.is_symmetric(), "symmetrize(x_0*x_1^2) over 3 vars should be symmetric");
    // Evaluate at (2, 3, 5) and (3, 2, 5) — should be equal.
    let v1 = sym.eval(&[2, 3, 5]).unwrap();
    let v2 = sym.eval(&[3, 2, 5]).unwrap();
    assert_eq!(v1, v2, "symmetrized poly: eval at (2,3,5) = eval at (3,2,5)");
}

// ─── KAT: SemaevError display ─────────────────────────────────────────────────

/// Smoke-test `SemaevError` Display impls (no panic).
#[test]
fn semaev_error_display_smoke() {
    let _ = format!("{}", SemaevError::DegreeZero);
    let _ = format!("{}", SemaevError::VariableOutOfRange { index: 5, num_vars: 3 });
    let _ = format!("{}", SemaevError::ArityMismatch { expected: 3, got: 2 });
    let _ = format!("{}", SemaevError::ZeroLeadingCoefficient);
}

// ─── KAT: univariate_in_var ───────────────────────────────────────────────────

/// `univariate_in_var` correctly extracts coefficients.
///
/// `f(x_0, x_1) = 3*x_0^2 + 2*x_0*x_1 + x_1^2 + x_0 + 5`.
/// Viewed as univariate in `x_0`: `3*x_0^2 + (2*x_1 + 1)*x_0 + (x_1^2 + 5)`.
#[test]
fn univariate_in_var_extraction() {
    let mut f = MultiPoly::zero(2, P);
    f.add_term(vec![2, 0], 3); // 3*x_0^2
    f.add_term(vec![1, 1], 2); // 2*x_0*x_1
    f.add_term(vec![0, 2], 1); // x_1^2
    f.add_term(vec![1, 0], 1); // x_0
    f.add_term(vec![0, 0], 5); // 5

    let coeffs = f.univariate_in_var(0).unwrap();
    // coeffs[0] = x_1^2 + 5 (constant term in x_0)
    // coeffs[1] = 2*x_1 + 1 (coefficient of x_0)
    // coeffs[2] = 3 (coefficient of x_0^2)
    assert_eq!(coeffs.len(), 3, "degree 2 in x_0 → 3 coefficient polynomials");

    // coeffs[2] should be the constant 3 (in x_1)
    let v2 = coeffs[2].eval(&[0]).unwrap();
    assert_eq!(v2, 3, "coefficient of x_0^2 should be 3");

    // coeffs[1] at x_1 = 4: 2*4 + 1 = 9
    let v1 = coeffs[1].eval(&[4]).unwrap();
    assert_eq!(v1, 9, "coefficient of x_0 at x_1=4: 2*4+1=9");

    // coeffs[0] at x_1 = 3: 3^2 + 5 = 14
    let v0 = coeffs[0].eval(&[3]).unwrap();
    assert_eq!(v0, 14, "constant term at x_1=3: 3^2+5=14");
}

// ═══════════════════════════════════════════════════════════════════════════════
// Base cases S_2, S_3 + the vanishing relation
// ═══════════════════════════════════════════════════════════════════════════════
//
// Known multiples of G on the toy curve y² = x³ + x + 33 over F_47:
//   G  = (10, 3),  −G  = (10, 44)
//   2G = (7, 30),  −2G = (7, 17)
//   3G = (17, 13), −3G = (17, 34)
//   4G = (23, 12), −4G = (23, 35)
//   5G = (32, 36), −5G = (32, 11)
//
// Verified: each point satisfies y² = x³ + x + 33 mod 47.

// ─── helpers (S_2, S_3 base cases) ───────────────────────────────────────────

/// Construct a finite affine point over `F_47`.
fn pt(x: u64, y: u64) -> AffinePoint<FpNaive> {
    let p_val = Uint::<4>::from(P);
    AffinePoint::Finite {
        x: FpNaive::from_u64(x, &p_val),
        y: FpNaive::from_u64(y, &p_val),
    }
}

// ─── KAT 11: S_2 vanishing ───────────────────────────────────────────────────

/// `S_2(x_1, x_2) = 0 ⟺ x_1 = x_2` (i.e. `P_2 = −P_1`).
///
/// `S_2(10, 10) = 0` and `G + (−G) = ∞` ✓.
#[test]
fn s2_vanishes_for_negation() {
    let c = semaev_toy();
    let g = pt(10, 3);
    let neg_g = pt(10, 44);
    assert!(
        vanishes_s2(&c, &g, &neg_g),
        "S_2 should vanish for G and −G (same x-coordinate)"
    );
}

/// `S_2(10, 7) ≠ 0` and `G + 2G = 3G ≠ ∞`.
#[test]
fn s2_nonzero_for_distinct_x() {
    let c = semaev_toy();
    let g = pt(10, 3);
    let two_g = pt(7, 30);
    assert!(
        !vanishes_s2(&c, &g, &two_g),
        "S_2 should not vanish for G and 2G (distinct x-coordinates)"
    );
}

/// `S_2` vanishes for any point and its negation (multiple spot-checks).
#[test]
fn s2_vanishes_for_multiple_negation_pairs() {
    let c = semaev_toy();
    // (2G, −2G): x=7
    assert!(vanishes_s2(&c, &pt(7, 30), &pt(7, 17)), "S_2 should vanish for 2G and −2G");
    // (3G, −3G): x=17
    assert!(vanishes_s2(&c, &pt(17, 13), &pt(17, 34)), "S_2 should vanish for 3G and −3G");
    // (4G, −4G): x=23
    assert!(vanishes_s2(&c, &pt(23, 12), &pt(23, 35)), "S_2 should vanish for 4G and −4G");
}

// ─── KAT 12: S_3 vanishing ───────────────────────────────────────────────────

/// `S_3(10, 7, 17) = 0`: x-coordinates of `G, 2G, −3G` — there exist y-values making
/// the sum `∞` (specifically `G + 2G + (−3G) = ∞`).
#[test]
fn s3_vanishes_for_g_2g_neg3g_xcoords() {
    let c = semaev_toy();
    assert!(
        vanishes_s3::<FpNaive>(&c, 10, 7, 17),
        "S_3 should vanish for x-coords of G, 2G, −3G"
    );
}

/// `S_3(10, 7, 23) ≠ 0`: x-coordinates of `G, 2G, 4G` — no y-values make sum `∞`.
///
/// All 8 combinations of y-values for x=10,7,23 are checked:
/// - `G+2G+4G = 7G ≠ ∞`, `G+2G+(−4G) = −G ≠ ∞`, `G+(−2G)+4G = 3G ≠ ∞`,
/// - `G+(−2G)+(−4G) = −5G ≠ ∞`, `(−G)+2G+4G = 5G ≠ ∞`, `(−G)+2G+(−4G) = −3G ≠ ∞`,
/// - `(−G)+(−2G)+4G = G ≠ ∞`, `(−G)+(−2G)+(−4G) = −7G ≠ ∞`.
#[test]
fn s3_nonzero_for_g_2g_4g_xcoords() {
    let c = semaev_toy();
    assert!(
        !vanishes_s3::<FpNaive>(&c, 10, 7, 23),
        "S_3 should not vanish for x-coords of G, 2G, 4G (no y-values make sum ∞)"
    );
}

/// `S_3(10, 17, 32) = 0`: x-coordinates of `G, 3G, −5G` — `G + 3G + (−5G) = ∞`.
///
/// `G + 3G = 4G = (23, 12)`. `−5G = (32, 11)`. `4G + (−5G) = −G = (10, 44)`.
/// Wait — that gives `G + 3G + (−5G) = 4G + (−5G) = −G ≠ ∞`.
/// Let's try `G + (−3G) + 2G = ∞`: x-coords 10, 17, 7 — same as KAT 12a.
/// Try `2G + 3G + (−5G) = 5G + (−5G) = ∞`: x-coords 7, 17, 32.
#[test]
fn s3_vanishes_for_2g_3g_neg5g_xcoords() {
    let c = semaev_toy();
    // 2G=(7,30), 3G=(17,13), −5G=(32,11): 2G+3G+(−5G) = 5G+(−5G) = ∞
    assert!(
        vanishes_s3::<FpNaive>(&c, 7, 17, 32),
        "S_3 should vanish for x-coords of 2G, 3G, −5G"
    );
}

/// `S_3(10, 10, 10) ≠ 0`: no choice of y-values for x=10 makes three points sum to `∞`.
///
/// The y-values for x=10 are `{3, 44}` (i.e. `G` and `−G`). All 8 combinations of
/// `(G or −G) + (G or −G) + (G or −G)` give `±G, ±3G` — none is `∞`.
#[test]
fn s3_nonzero_for_all_same_x() {
    let c = semaev_toy();
    assert!(
        !vanishes_s3::<FpNaive>(&c, 10, 10, 10),
        "S_3 should not vanish for x-coords (10,10,10)"
    );
}

/// `S_3(10, 7, 32) ≠ 0`: x-coords of `G, 2G, 5G` — no y-values make sum `∞`.
///
/// `G + 2G + 5G = 8G ≠ ∞`. `G + 2G + (−5G) = 3G + (−5G) = −2G ≠ ∞`.
/// `G + (−2G) + 5G = 4G ≠ ∞`. `G + (−2G) + (−5G) = −6G ≠ ∞`.
/// `(−G) + 2G + 5G = 6G ≠ ∞`. `(−G) + 2G + (−5G) = −4G ≠ ∞`.
/// `(−G) + (−2G) + 5G = 2G ≠ ∞`. `(−G) + (−2G) + (−5G) = −8G ≠ ∞`.
#[test]
fn s3_nonzero_for_g_2g_5g_xcoords() {
    let c = semaev_toy();
    assert!(
        !vanishes_s3::<FpNaive>(&c, 10, 7, 32),
        "S_3 should not vanish for x-coords of G, 2G, 5G"
    );
}

// ─── KAT 13: S_3 symmetry ────────────────────────────────────────────────────

/// `S_3` is symmetric: invariant under all permutations of its three arguments.
///
/// Verified by evaluating at all 6 permutations of `(10, 7, 17)` and checking
/// that all give the same value.
#[test]
fn s3_symmetric_eval_at_all_permutations() {
    let poly = s3(1, 33, P);
    assert!(poly.is_symmetric(), "S_3 should be symmetric");

    // All 6 permutations of (10, 7, 17)
    let perms = [
        [10u64, 7, 17],
        [10, 17, 7],
        [7, 10, 17],
        [7, 17, 10],
        [17, 10, 7],
        [17, 7, 10],
    ];
    let base_val = poly.eval(&perms[0]).unwrap();
    for perm in &perms[1..] {
        let v = poly.eval(perm).unwrap();
        assert_eq!(
            v, base_val,
            "S_3 should be symmetric: eval at {:?} = eval at {:?}",
            perm, perms[0]
        );
    }
}

/// `S_3` symmetry on a non-vanishing triple: all permutations of `(10, 10, 10)` agree.
#[test]
fn s3_symmetric_eval_nonzero_triple() {
    let poly = s3(1, 33, P);
    // All permutations of (10, 10, 10) are the same — trivially symmetric.
    let v = poly.eval(&[10, 10, 10]).unwrap();
    assert_ne!(v, 0, "S_3(10,10,10) should be nonzero");

    // Verify with a non-trivial permutation: (7, 17, 32) — all distinct.
    // (This triple has S_3 ≠ 0 since no y-values make G+3G+5G or similar = ∞.)
    let perms = [
        [7u64, 17, 32],
        [7, 32, 17],
        [17, 7, 32],
        [17, 32, 7],
        [32, 7, 17],
        [32, 17, 7],
    ];
    let base_val = poly.eval(&perms[0]).unwrap();
    for perm in &perms[1..] {
        let v = poly.eval(perm).unwrap();
        assert_eq!(
            v, base_val,
            "S_3 should be symmetric: eval at {:?} = eval at {:?}",
            perm, perms[0]
        );
    }
}

// ─── KAT 14: S_3 degree ──────────────────────────────────────────────────────

/// `S_3` has degree ≤ 2 in each variable.
///
/// Verified by checking that no exponent in any monomial exceeds 2 for any variable.
#[test]
fn s3_degree_at_most_2_in_each_var() {
    let poly = s3(1, 33, P);
    for (exp, _coeff) in &poly.terms {
        for (i, &e) in exp.iter().enumerate() {
            assert!(
                e <= 2,
                "S_3 has degree > 2 in variable {i}: exponent vector {:?}",
                exp
            );
        }
    }
}

/// `S_3` achieves degree exactly 2 in each variable (not just ≤ 2).
///
/// Verified by checking that there exists a monomial with degree 2 in each variable.
#[test]
fn s3_degree_exactly_2_in_each_var() {
    let poly = s3(1, 33, P);
    for var in 0..3 {
        let max_deg = poly.terms.keys().map(|exp| exp[var]).max().unwrap_or(0);
        assert_eq!(
            max_deg, 2,
            "S_3 should have degree exactly 2 in variable {var}, got {max_deg}"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Resultant recursion `S_m`
// ═══════════════════════════════════════════════════════════════════════════════
//
// Known multiples of G on the toy curve y² = x³ + x + 33 over F_47:
//   G  = (10, 3),  −G  = (10, 44)
//   2G = (7, 30),  −2G = (7, 17)
//   3G = (17, 13), −3G = (17, 34)
//   4G = (23, 12), −4G = (23, 35)
//   5G = (32, 36), −5G = (32, 11)
//   6G = (19, 7),  −6G = (19, 40)   [6G = 5G + G, computed via group law]
//   8G = (25, 28), −8G = (25, 19)   [8G = 2*(4G), computed via doubling]
//
// Verification of 6G = (19, 7):
//   slope = (3 - 36) / (10 - 32) = (-33) / (-22) = 33 * 22^{-1} mod 47
//   22^{-1} mod 47: 22*15 = 330 = 7*47 + 1 → 22^{-1} = 15
//   slope = 33*15 = 495 = 10*47 + 25 → slope = 25
//   x_6G = 25^2 - 32 - 10 = 625 - 42 = 583 = 12*47 + 19 → x_6G = 19
//   y_6G = 25*(32 - 19) - 36 = 325 - 36 = 289 = 6*47 + 7 → y_6G = 7
//   Check: 19^3 + 19 + 33 = 6859 + 19 + 33 = 6911 = 147*47 + 2 → 2. 7^2 = 49 = 47 + 2 → 2. ✓
//
// Verification of 8G = (25, 28):
//   Doubling (23, 12): slope = (3*23^2 + 1) / (2*12) = (1588) / 24 mod 47
//   1588 mod 47: 33*47 = 1551, 1588 - 1551 = 37 → numerator = 37
//   24^{-1} mod 47: 24*2 = 48 = 47 + 1 → 24^{-1} = 2
//   slope = 37*2 = 74 = 47 + 27 → slope = 27
//   x_8G = 27^2 - 2*23 = 729 - 46 = 683 = 14*47 + 25 → x_8G = 25
//   y_8G = 27*(23 - 25) - 12 = -54 - 12 = -66 = -66 + 2*47 = 28 → y_8G = 28
//   Check: 25^3 + 25 + 33 = 15683 mod 47. 333*47 = 15651, 15683 - 15651 = 32. 28^2 = 784 mod 47.
//   16*47 = 752, 784 - 752 = 32. ✓

// ─── helpers (S_m resultant recursion) ───────────────────────────────────────

/// Check whether any combination of y-values for 4 x-coordinates makes
/// `P_1 + P_2 + P_3 + P_4 = ∞` via the group law.
fn exists_summing_quad(x1: u64, x2: u64, x3: u64, x4: u64) -> bool {
    use crypto_bigint::Uint;

    let c = semaev_toy();
    let p_val = c.p;
    let p_uint = Uint::<4>::from(P);

    // Find all y-values for each x-coordinate.
    let find_ys = |x: u64| -> Vec<u64> {
        let xf = FpNaive::from_u64(x, &p_uint);
        let a = FpNaive::from_u64(1, &p_uint);
        let b = FpNaive::from_u64(33, &p_uint);
        let rhs = xf.square(&p_uint).mul(&xf, &p_uint).add(&a.mul(&xf, &p_uint), &p_uint).add(&b, &p_uint);
        if rhs.is_zero(&p_uint) {
            return vec![0];
        }
        // Legendre symbol: rhs^((p-1)/2) mod p
        let mut exp = p_uint.wrapping_sub(&Uint::<4>::ONE);
        exp >>= 1;
        let leg = rhs.pow(&exp, &p_uint);
        if !leg.is_one(&p_uint) {
            return vec![];
        }
        // p = 47 ≡ 3 mod 4 → sqrt = rhs^((p+1)/4)
        let mut exp4 = p_uint.wrapping_add(&Uint::<4>::ONE);
        exp4 >>= 2;
        let y = rhs.pow(&exp4, &p_uint);
        let y_u64 = y.to_uint().as_words()[0];
        let neg_y_u64 = (P - y_u64) % P;
        if y_u64 == neg_y_u64 {
            vec![y_u64]
        } else {
            vec![y_u64, neg_y_u64]
        }
    };

    let ys1 = find_ys(x1);
    let ys2 = find_ys(x2);
    let ys3 = find_ys(x3);
    let ys4 = find_ys(x4);

    let make_pt = |x: u64, y: u64| -> AffinePoint<FpNaive> {
        AffinePoint::Finite {
            x: FpNaive::from_u64(x, &p_uint),
            y: FpNaive::from_u64(y, &p_uint),
        }
    };

    for &y1 in &ys1 {
        for &y2 in &ys2 {
            for &y3 in &ys3 {
                for &y4 in &ys4 {
                    let p1 = make_pt(x1, y1);
                    let p2 = make_pt(x2, y2);
                    let p3 = make_pt(x3, y3);
                    let p4 = make_pt(x4, y4);
                    let j1 = JacobianPoint::from_affine(&p1, &p_val);
                    let j12 = c.add_mixed(&j1, &p2);
                    let j123 = c.add_mixed(&j12, &p3);
                    let j1234 = c.add_mixed(&j123, &p4);
                    if j1234.to_affine(&p_val).is_infinity() {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Check the `S_4` vanishing relation: `S_4(x_1,x_2,x_3,x_4) = 0 ⟺ ∃ y_i: Σ P_i = ∞`.
///
/// Both the polynomial evaluation and the group-law existence check must agree.
///
/// # Panics
///
/// Panics if the two conditions disagree — this is a correctness invariant violation.
fn vanishes_s4(x1: u64, x2: u64, x3: u64, x4: u64) -> bool {
    let poly = semaev_poly(4, 1, 33, P).unwrap();
    let poly_val = poly.eval(&[x1, x2, x3, x4]).unwrap();
    let poly_zero = poly_val == 0;
    let exists_inf = exists_summing_quad(x1, x2, x3, x4);
    assert_eq!(
        poly_zero, exists_inf,
        "S_4 vanishing disagreement at ({x1},{x2},{x3},{x4}): \
         poly_zero={poly_zero}, exists_inf={exists_inf}"
    );
    poly_zero
}

// ─── KAT 15: recursion base cases ────────────────────────────────────────────

/// `semaev_poly(2)` returns the same polynomial as `s2` directly.
#[test]
fn semaev_poly_m2_matches_s2_direct() {
    let via_recursion = semaev_poly(2, 1, 33, P).unwrap();
    let direct = base_s2(P);
    assert_eq!(
        via_recursion, direct,
        "semaev_poly(2) should match s2 directly"
    );
}

/// `semaev_poly(3)` returns the same polynomial as `s3` directly.
#[test]
fn semaev_poly_m3_matches_s3_direct() {
    let via_recursion = semaev_poly(3, 1, 33, P).unwrap();
    let direct = s3(1, 33, P);
    assert_eq!(
        via_recursion, direct,
        "semaev_poly(3) should match s3 directly"
    );
}

/// `semaev_poly(m)` returns an error for `m < 2`.
#[test]
fn semaev_poly_m_lt_2_is_error() {
    assert!(semaev_poly(0, 1, 33, P).is_err(), "semaev_poly(0) should error");
    assert!(semaev_poly(1, 1, 33, P).is_err(), "semaev_poly(1) should error");
}

// ─── KAT 16: S_4 structure ───────────────────────────────────────────────────

/// `S_4` has exactly 4 variables.
#[test]
fn s4_has_4_variables() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    assert_eq!(s4.num_vars, 4, "S_4 should have 4 variables");
}

/// `S_4` is a non-zero polynomial.
#[test]
fn s4_is_nonzero_polynomial() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    assert!(!s4.is_zero(), "S_4 should be a non-zero polynomial");
}

// ─── KAT 17: S_4 vanishing (primary correctness signal) ──────────────────────

/// `S_4(10, 10, 7, 7) = 0`: x-coords of `G, −G, 2G, −2G` — `G + (−G) + 2G + (−2G) = ∞`.
///
/// This is the simplest vanishing case: two cancelling pairs. The polynomial evaluation
/// and the group-law existence check must agree.
#[test]
fn s4_vanishes_for_g_neg_g_2g_neg_2g() {
    assert!(
        vanishes_s4(10, 10, 7, 7),
        "S_4 should vanish for x-coords of G, −G, 2G, −2G"
    );
}

/// `S_4(10, 7, 17, 19) = 0`: x-coords of `G, 2G, 3G, −6G` — `G + 2G + 3G + (−6G) = ∞`.
///
/// `1 + 2 + 3 − 6 = 0 mod 60`. `6G = (19, 7)`, `−6G = (19, 40)`.
/// This is the primary non-trivial vanishing KAT: four distinct x-coordinates.
#[test]
fn s4_vanishes_for_g_2g_3g_neg6g() {
    assert!(
        vanishes_s4(10, 7, 17, 19),
        "S_4 should vanish for x-coords of G, 2G, 3G, −6G"
    );
}

/// `S_4(10, 7, 23, 25) ≠ 0`: x-coords of `G, 2G, 4G, 8G` — no y-values make sum `∞`.
///
/// For all 16 sign combinations `ε_1·1 + ε_2·2 + ε_3·4 + ε_4·8` (ε_i ∈ {±1}),
/// the result is never `0 mod 60`. Verified exhaustively: the values are
/// `±1 ± 2 ± 4 ± 8 ∈ {±15, ±11, ±7, ±3, ±9, ±5, ±1, ±13}` — none is `0 mod 60`.
/// `8G = (25, 28)` (computed via doubling `4G = (23, 12)`).
#[test]
fn s4_nonzero_for_g_2g_4g_8g() {
    assert!(
        !vanishes_s4(10, 7, 23, 25),
        "S_4 should not vanish for x-coords of G, 2G, 4G, 8G"
    );
}

/// `S_4` vanishing predicate agrees with the group law for a non-summing quadruple.
///
/// `S_4(10, 7, 17, 23) = 0` because `−G + 2G + 3G + (−4G) = 0` (sign combination
/// `−1 + 2 + 3 − 4 = 0 mod 60`). The polynomial correctly detects this existential
/// y-value combination.
#[test]
fn s4_vanishes_for_neg_g_2g_3g_neg4g() {
    // x-coords: G=10, 2G=7, 3G=17, 4G=23
    // -G + 2G + 3G + (-4G) = -1+2+3-4 = 0 mod 60 → ∞
    assert!(
        vanishes_s4(10, 7, 17, 23),
        "S_4 should vanish for x-coords 10,7,17,23 (−G+2G+3G+(−4G) = ∞)"
    );
}

// ─── KAT 18: S_4 degree growth ───────────────────────────────────────────────

/// `S_4` has degree ≤ 4 in each variable.
///
/// The recursion `S_4 = Res_X(S_3(X_1,X_2,X), S_3(X_3,X_4,X))` eliminates `X` from
/// two degree-2 polynomials in `X`, producing a degree `2·2 = 4` polynomial in the
/// remaining variables. The expected degree is `2^(m-2) = 4` for `m = 4`.
#[test]
fn s4_degree_at_most_4_in_each_var() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    for (exp, _coeff) in &s4.terms {
        for (i, &e) in exp.iter().enumerate() {
            assert!(
                e <= 4,
                "S_4 has degree > 4 in variable {i}: exponent vector {:?}",
                exp
            );
        }
    }
}

/// `S_4` achieves degree exactly 4 in each variable.
///
/// Verified by checking that there exists a monomial with degree 4 in each variable.
#[test]
fn s4_degree_exactly_4_in_each_var() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    for var in 0..4 {
        let max_deg = s4.terms.keys().map(|exp| exp[var]).max().unwrap_or(0);
        assert_eq!(
            max_deg, 4,
            "S_4 should have degree exactly 4 in variable {var}, got {max_deg}"
        );
    }
}

// ─── KAT 19: S_4 symmetry ────────────────────────────────────────────────────

/// `S_4` is symmetric: invariant under all permutations of its 4 arguments.
#[test]
fn s4_is_symmetric() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    assert!(s4.is_symmetric(), "S_4 should be symmetric in all 4 variables");
}

/// `S_4` evaluates to the same value at all permutations of a vanishing quadruple.
///
/// All permutations of `(10, 10, 7, 7)` should give `S_4 = 0`.
#[test]
fn s4_symmetric_eval_at_permutations_of_vanishing_quad() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    // All distinct permutations of (10, 10, 7, 7) — 6 distinct orderings.
    let perms: &[[u64; 4]] = &[
        [10, 10, 7, 7],
        [10, 7, 10, 7],
        [10, 7, 7, 10],
        [7, 10, 10, 7],
        [7, 10, 7, 10],
        [7, 7, 10, 10],
    ];
    for perm in perms {
        let v = s4.eval(perm).unwrap();
        assert_eq!(
            v, 0,
            "S_4 should vanish at all permutations of (10,10,7,7); got {v} at {:?}",
            perm
        );
    }
}

/// `S_4` evaluates to the same value at all permutations of a non-vanishing quadruple.
///
/// All permutations of `(10, 7, 23, 25)` should give the same non-zero value.
#[test]
fn s4_symmetric_eval_at_permutations_of_nonzero_quad() {
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    let base_val = s4.eval(&[10, 7, 23, 25]).unwrap();
    assert_ne!(base_val, 0, "S_4(10,7,23,25) should be nonzero");

    // A sample of permutations (not all 24 — 6 is sufficient to confirm symmetry).
    let perms: &[[u64; 4]] = &[
        [7, 10, 23, 25],
        [23, 7, 10, 25],
        [25, 23, 7, 10],
        [10, 25, 7, 23],
        [23, 10, 25, 7],
    ];
    for perm in perms {
        let v = s4.eval(perm).unwrap();
        assert_eq!(
            v, base_val,
            "S_4 should be symmetric: eval at {:?} = {base_val}, got {v}",
            perm
        );
    }
}

// ─── KAT 20: S_4 vanishing predicate agreement ───────────────────────────────

/// `S_4` vanishing agrees with the group law for multiple quadruples.
///
/// This is the decisive correctness KAT: the polynomial vanishing condition and
/// the group-law existence condition agree for a range of x-coordinate quadruples.
#[test]
fn s4_vanishing_predicate_agrees_with_group_law() {
    // Vanishing cases: quadruples where ∃ y-values making sum ∞.
    let vanishing = [
        (10u64, 10, 7, 7),   // G + (−G) + 2G + (−2G) = ∞
        (10, 7, 17, 19),     // G + 2G + 3G + (−6G) = ∞
        (10, 7, 17, 23),     // −G + 2G + 3G + (−4G) = ∞
        (17, 17, 23, 23),    // 3G + (−3G) + 4G + (−4G) = ∞
    ];
    for (x1, x2, x3, x4) in vanishing {
        assert!(
            vanishes_s4(x1, x2, x3, x4),
            "S_4 should vanish for ({x1},{x2},{x3},{x4})"
        );
    }

    // Non-vanishing cases: quadruples where no y-values make sum ∞.
    let non_vanishing = [
        (10u64, 7, 23, 25),  // G, 2G, 4G, 8G — no sign combination sums to 0 mod 60
    ];
    for (x1, x2, x3, x4) in non_vanishing {
        assert!(
            !vanishes_s4(x1, x2, x3, x4),
            "S_4 should not vanish for ({x1},{x2},{x3},{x4})"
        );
    }
}

// ─── KAT 21: optional PARI sidecar ──────────────────────────────────────────

/// Checks the vanishing of `S_4` at fixed `x`-coordinates, following the PARI/GP script below.
///
/// The test does not call `gp`: its body evaluates `semaev_poly` directly and needs no external
/// tool. It is `#[ignore]`d.
///
/// The PARI/GP script the expected values follow:
///
/// ```text
/// gp> p = 47; a = 1; b = 33;
/// gp> \\ S_3(x1, x2, x3) = e2^2 - 2*a*e2 - 4*e1*e3 - 4*b*e1 + a^2
/// gp> \\ where e1 = x1+x2+x3, e2 = x1*x2+x1*x3+x2*x3, e3 = x1*x2*x3
/// gp> s3_val(x1, x2, x3) = {
/// ...   e1 = x1+x2+x3; e2 = x1*x2+x1*x3+x2*x3; e3 = x1*x2*x3;
/// ...   (e2^2 - 2*a*e2 - 4*e1*e3 - 4*b*e1 + a^2) % p
/// ... }
/// gp> \\ Verify S_3 vanishes for G+2G+(-3G) = ∞: x-coords 10, 7, 17
/// gp> s3_val(10, 7, 17) % p
/// 0
/// gp> \\ Verify S_4 vanishes for G+(-G)+2G+(-2G) = ∞: x-coords 10, 10, 7, 7
/// gp> \\ S_4 = Res_X(S_3(10, 10, X), S_3(7, 7, X)) — substitute first two x-coords
/// gp> \\ and compute resultant in X, then evaluate at remaining coords.
/// ```
#[test]
#[ignore = "asserts PARI/GP-derived values; does not call gp"]
fn s4_pari_cross_check() {
    // The body evaluates S_4 directly; it does not run a PARI subprocess.
    //
    // Expected: S_4(10, 10, 7, 7) = 0 (G + (-G) + 2G + (-2G) = ∞).
    // Expected: S_4(10, 7, 17, 19) = 0 (G + 2G + 3G + (-6G) = ∞).
    // Expected: S_4(10, 7, 23, 25) ≠ 0 (no y-values make G + 2G + 4G + 8G = ∞).
    let s4 = semaev_poly(4, 1, 33, P).unwrap();
    assert_eq!(s4.eval(&[10, 10, 7, 7]).unwrap(), 0);
    assert_eq!(s4.eval(&[10, 7, 17, 19]).unwrap(), 0);
    assert_ne!(s4.eval(&[10, 7, 23, 25]).unwrap(), 0);
}
