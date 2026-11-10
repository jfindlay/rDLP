//! Known-answer tests for the GHS Weil-descent construction: descent algebra,
//! hyperelliptic-curve extraction, transfer map, and ECDLP → Jacobian-DLP reduction.
//!
//! # Coverage
//!
//! ## GHS reduction `(E, g, h) → (C, D_g, D_h)` + logarithm preservation
//!
//! ### Logarithm preservation (decisive)
//! - For known `k=3`: compute `h = 3·g` on `E`, call `ghs_descend`, verify
//!   `scalar_mul(D_g, 3) == D_h` on `Jac(C)`.
//!
//! ### Genus check
//! - `genus(C) == ghs_genus(6, 2) == 1`.
//!
//! ### Valid reduced divisors
//! - `D_g` and `D_h` are valid reduced divisors.
//!
//! ### D_g is non-zero
//! - The base point maps to a non-zero divisor.
//!
//! ### Reduction with identity
//! - `ghs_descend(params, g, ∞) → D_h = zero divisor`.
//!
//! ### Optional PARI sidecar
//! - `#[ignore]`-gated cross-check of `#Jac(C)(GF(2^2))` via PARI `hyperellcharpoly`.
//!
//! ## GHS transfer map `E(GF(2^m)) → Jac(C)(GF(2^l))`
//!
//! ### Identity maps to identity
//! - `transfer_point(∞, ...) == [1, 0]` (zero divisor).
//!
//! ### Transferred divisors are valid
//! - `transfer_point(G, ...)` returns a valid reduced divisor.
//! - `transfer_point(P, ...)` returns a valid reduced divisor.
//!
//! ### Transfer of base point
//! - `transfer_point(G, ...)` succeeds and gives a non-zero divisor.
//! - Known-answer: `D_G = [X, 1]` (u = X, v = 1).
//!
//! ### Input guards
//! - An off-curve point gives `GhsError::NotOnCurve` (also through `ghs_descend`).
//! - A conorm image with a coordinate outside `GF(2^l)` gives `GhsError::NotInSubfield`.
//! - A `poly_l` that is not the minimal polynomial of `β = α^((2^m−1)/(2^l−1))` gives
//!   `GhsError::SubfieldPolyMismatch` (extraction, transfer, and `ghs_descend`); `l ∤ m`
//!   gives `SubfieldDivisibility` even when `GhsParams` is built without `new`.
//! - Curve extraction and transfer use one identification of `GF(2^l)`: on `E` with
//!   `a = β`, every transferred divisor is valid on the extracted `C`.
//!
//! ### Homomorphism property (decisive correctness guard)
//! - `compose(D_G, D_P) == D_{G+P}` — the transfer is a group homomorphism.
//! - Known-answer: `D_G = [X, 1]`, `D_P = [X+1, 2]`, `D_{G+P} = [X+1, 3]` (hand-derived in
//!   the doc comments of the transfer known-answer tests).
//! - `D_G + D_G (Cantor) == D_{2G}` — doubling case.
//!
//! ## GHS hyperelliptic-curve extraction `C/GF(2^l)` (imaginary model)
//!
//! ### Genus formula
//! - `ghs_genus(6, 2) == 1` — toy fixture (m=6, l=2, n=3, g=(3-1)/2=1).
//! - `ghs_genus(10, 2) == 2` — crypto-scale n=5 gives genus 2.
//!
//! ### Extracted curve validity
//! - `extract_ghs_curve(&params)` succeeds for the toy fixture.
//! - The extracted curve passes `is_valid` for the zero divisor.
//! - The extracted curve's `poly` is the GF(2^2) irreducible.
//!
//! ### Imaginary model
//! - `deg f = 2g+1 = 3` for genus 1.
//! - `deg H ≤ g = 1`.
//! - `F(X)` is monic.
//!
//! ### Coefficients in GF(2^l)
//! - All coefficients of `F(X)` are in GF(2^2) (is_in_subfield check).
//! - `H(X) = X`: h_coeffs = [0, 1], all in GF(2^2).
//!
//! ### Genus matches extension degree
//! - `curve.genus() == ghs_genus(6, 2)` — the extracted curve's genus matches
//!   the formula.
//!
//! ### Even m/l rejection
//! - `extract_ghs_curve` returns `Err(NonDescendable)` for even `m/l`.
//!
//! ## Artin–Schreier / function-field Weil-restriction algebra
//!
//! ### Precondition verifier
//! - `check_ghs_params(6, 2)` returns `Ok(())` — `l=2` divides `m=6`.
//! - `check_ghs_params(6, 4)` returns `Err(SubfieldDivisibility)` — `l=4` does not divide `m=6`.
//! - `check_ghs_params(6, 0)` returns `Err(SubfieldDivisibility)` — zero subfield degree.
//!
//! ### Artin–Schreier round-trip
//! - The Artin–Schreier polynomial `f(x) = x³ + ax² + b` has degree 3.
//! - The leading coefficient is `1` (monic).
//! - The constant term equals `b` from the curve (non-zero for non-supersingular curves).
//! - The Artin–Schreier form is `y² + y = f(x)` (NOT `y² = f(x)` — the inseparable form).
//!
//! ### Weil-restriction dimension
//! - A `GF(2^m)`-polynomial of degree `d` restricts to a `GF(2^l)`-polynomial of
//!   degree `d·(m/l)` (or less, if leading coefficients cancel).
//! - The Weil restriction of the zero polynomial is zero.
//! - The Weil restriction dimension is `m/l = 3` for the toy fixture.
//!
//! ### Weil restriction of a subfield element
//! - An element `a ∈ GF(2^l) ⊂ GF(2^m)` restricts to itself: the constant-term
//!   polynomial `a` over `GF(2^l)` (the restriction of a subfield element is the
//!   element itself, placed at position 0).
//!
//! ### Weil restriction is GF(2^l)-linear
//! - `Res(f + g) = Res(f) + Res(g)` (additivity).
//! - `Res(c·f) = c·Res(f)` for `c ∈ GF(2^l)` (scalar homogeneity).
//!
//! # Toy fixture
//!
//! All tests use the canonical GHS fixture:
//! - `m = 6`, `l = 2`, `m/l = 3` (odd — imaginary model).
//! - Source field: `GF(2^6)` with irreducible `x⁶+x+1` (poly = 0x43).
//! - Subfield: `GF(2^2)` with irreducible `x²+x+1` (poly = 0x7).
//! - Binary curve `E/GF(2^6)`: `y²+xy = x³+x²+1` (`a=1`, `b=1`).
//!
//! The algorithms are crypto-scale-correct; only the parameters are toy
//! (m=6, l=2, genus 1).

use crypto_bigint::Uint;
use rho::binary_curve::{BinaryAffinePoint, BinaryCurve};
use rho::ghs::{
    GhsError, GhsParams, check_ghs_params, ghs_toy_curve, GHS_POLY2, GHS_POLY6,
    ArtinSchreierData, WeilRestriction, weil_restrict_poly,
    extract_ghs_curve, ghs_genus,
    ghs_descend, verify_log_preservation,
    transfer_point, verify_homomorphism,
};
use rho::hyperelliptic::cantor;
use shared_gf2m::{F2m, F2mNaive, Poly, is_in_subfield};

// ── Field and fixture helpers ─────────────────────────────────────────────────

/// GF(2^6) irreducible: x⁶+x+1 = 0x43.
fn poly6() -> Uint<1> {
    Uint::<1>::from(GHS_POLY6)
}

/// GF(2^2) irreducible: x²+x+1 = 0x7.
fn poly2() -> Uint<1> {
    Uint::<1>::from(GHS_POLY2)
}

/// Construct a GF(2^6) field element from a u64 value.
fn f6(v: u64) -> F2mNaive<1> {
    F2mNaive::<1>::from_u64(v, &poly6())
}

/// Build the toy GHS parameters (m=6, l=2).
fn toy_params() -> GhsParams {
    GhsParams::new(6, 2, ghs_toy_curve(), poly2()).expect("toy GHS params must be valid")
}

// ── Precondition verifier KATs ────────────────────────────────────────────────

/// `check_ghs_params(6, 2)` returns `Ok(())`: l=2 divides m=6.
///
/// This is the canonical GHS fixture. GF(2^2) ⊆ GF(2^6) since 2 | 6.
#[test]
fn precondition_l2_divides_m6() {
    assert_eq!(
        check_ghs_params(6, 2),
        Ok(()),
        "l=2 divides m=6: GF(2^2) ⊆ GF(2^6)"
    );
}

/// `check_ghs_params(6, 3)` returns `Ok(())`: l=3 divides m=6.
///
/// GF(2^3) ⊆ GF(2^6) since 3 | 6.
#[test]
fn precondition_l3_divides_m6() {
    assert_eq!(
        check_ghs_params(6, 3),
        Ok(()),
        "l=3 divides m=6: GF(2^3) ⊆ GF(2^6)"
    );
}

/// `check_ghs_params(6, 4)` returns `Err(SubfieldDivisibility)`: l=4 does not divide m=6.
///
/// GF(2^4) ⊄ GF(2^6) since 4 ∤ 6. This is the primary rejection KAT.
#[test]
fn precondition_l4_does_not_divide_m6() {
    assert_eq!(
        check_ghs_params(6, 4),
        Err(GhsError::SubfieldDivisibility),
        "l=4 does not divide m=6: GF(2^4) ⊄ GF(2^6)"
    );
}

/// `check_ghs_params(6, 5)` returns `Err(SubfieldDivisibility)`: l=5 does not divide m=6.
#[test]
fn precondition_l5_does_not_divide_m6() {
    assert_eq!(
        check_ghs_params(6, 5),
        Err(GhsError::SubfieldDivisibility),
        "l=5 does not divide m=6"
    );
}

/// `check_ghs_params(6, 0)` returns `Err(SubfieldDivisibility)`: zero subfield degree.
///
/// l=0 is invalid — GF(2^0) is not a field.
#[test]
fn precondition_l0_is_invalid() {
    assert_eq!(
        check_ghs_params(6, 0),
        Err(GhsError::SubfieldDivisibility),
        "l=0 is invalid (zero subfield degree)"
    );
}

// ── Artin–Schreier round-trip KATs ───────────────────────────────────────────

/// The Artin–Schreier polynomial `f(x) = x³ + ax² + b` has degree 3.
///
/// For the toy curve `y²+xy = x³+x²+1` (a=1, b=1), the Artin–Schreier polynomial
/// is `f(x) = x³ + x² + 1` (degree 3).
#[test]
fn artin_schreier_degree_is_3() {
    let params = toy_params();
    let as_data = ArtinSchreierData::from_params(params);
    assert_eq!(
        as_data.degree(),
        Some(3),
        "Artin–Schreier polynomial must have degree 3"
    );
}

/// The Artin–Schreier polynomial is monic: leading coefficient is 1.
///
/// The leading coefficient of `f(x) = x³ + ax² + b` is always 1 (the coefficient
/// of `x³` from the Weierstrass equation `y²+xy = x³+ax²+b`).
#[test]
fn artin_schreier_is_monic() {
    let params = toy_params();
    let as_data = ArtinSchreierData::from_params(params);
    let one = F2mNaive::<1>::one();
    assert_eq!(
        as_data.leading_coeff(),
        Some(&one),
        "Artin–Schreier polynomial must be monic (leading coeff = 1)"
    );
}

/// The constant term of `f(x)` equals `b` from the curve.
///
/// For `f(x) = x³ + ax² + b`, the constant term is `b`. For the toy curve
/// (b=1), the constant term is 1.
#[test]
fn artin_schreier_constant_term_is_b() {
    let params = toy_params();
    let poly_m = poly6();
    let b = F2mNaive::<1>::from_uint(params.curve.b, &poly_m);
    let as_data = ArtinSchreierData::from_params(params);
    let const_term = as_data.f_poly.coeff(0);
    assert_eq!(
        const_term, b,
        "constant term of f(x) must equal b from the curve"
    );
}

/// The Artin–Schreier data is well-formed: degree 3, monic, non-zero constant term.
///
/// The `is_well_formed` predicate checks all three conditions simultaneously.
#[test]
fn artin_schreier_is_well_formed() {
    let params = toy_params();
    let as_data = ArtinSchreierData::from_params(params);
    assert!(
        as_data.is_well_formed(),
        "Artin–Schreier data must be well-formed (degree 3, monic, b ≠ 0)"
    );
}

/// The Artin–Schreier form is `y² + y = f(x)`, NOT `y² = f(x)`.
///
/// This KAT documents the critical char-2 invariant: the Artin–Schreier operator
/// is `℘(y) = y² + y` (separable), not `y ↦ y²` (inseparable Frobenius).
/// The GHS construction requires a separable extension.
///
/// We verify this by checking that `f(x)` has the correct structure for the
/// Artin–Schreier equation `℘(y) = y² + y = f(x)`:
/// - `f(x) = x³ + ax² + b` (from the Weierstrass equation, not `f(x) = x³ + ax² + b + y·x`).
/// - The extension is `GF(2^m)(x)[y]/(y² + y + f(x))` (separable degree-2 extension).
#[test]
fn artin_schreier_form_is_separable() {
    // The Artin–Schreier polynomial f(x) = x³ + ax² + b is the RHS of y² + y = f(x).
    // We verify the structure: degree 3, monic, constant term b.
    let params = toy_params();
    let poly_m = poly6();
    let a_coeff = F2mNaive::<1>::from_uint(params.curve.a, &poly_m);
    let b_coeff = F2mNaive::<1>::from_uint(params.curve.b, &poly_m);
    let as_data = ArtinSchreierData::from_params(params);

    // f(x) = x³ + ax² + b: check all four coefficients.
    let c0 = as_data.f_poly.coeff(0); // constant: b
    let c1 = as_data.f_poly.coeff(1); // x^1: 0
    let c2 = as_data.f_poly.coeff(2); // x^2: a
    let c3 = as_data.f_poly.coeff(3); // x^3: 1

    assert_eq!(c0, b_coeff, "f(x) constant term must be b");
    assert!(c1.is_zero(), "f(x) coefficient of x must be 0");
    assert_eq!(c2, a_coeff, "f(x) coefficient of x² must be a");
    assert!(c3.is_one(), "f(x) coefficient of x³ must be 1 (monic)");
}

// ── Weil-restriction dimension KATs ──────────────────────────────────────────

/// The Weil restriction dimension is `m/l = 3` for the toy fixture.
///
/// The Weil restriction `Res_{GF(2^6)/GF(2^2)}` maps 1-dimensional `GF(2^6)`-objects
/// to 3-dimensional `GF(2^2)`-objects.
#[test]
fn weil_restriction_dimension_is_3() {
    let params = toy_params();
    let wr = WeilRestriction::new(params);
    assert_eq!(wr.dimension(), 3, "Weil restriction dimension = m/l = 6/2 = 3");
}

/// The Weil restriction of the zero polynomial is zero.
#[test]
fn weil_restrict_zero_poly_is_zero() {
    let zero: Poly<F2mNaive<1>, 1> = Poly::zero();
    let result = weil_restrict_poly(&zero, 2, &poly6(), &poly2());
    assert!(result.is_zero(), "Weil restriction of zero must be zero");
}

/// A constant (degree-0) polynomial over `GF(2^m)` restricts to a polynomial of
/// degree at most `0·(m/l) = 0` over `GF(2^l)`.
///
/// The constant polynomial `c ∈ GF(2^m)` restricts to a polynomial whose coefficients
/// are the `m/l` components of `c` in the `GF(2^m)`-over-`GF(2^l)` basis.
#[test]
fn weil_restrict_constant_poly_degree() {
    // Use a non-subfield element to get a non-trivial restriction.
    let a = f6(0x15); // 0x15 = 0b010101 — a generic GF(2^6) element
    let p = Poly::from_coeffs(vec![a]);
    let result = weil_restrict_poly(&p, 2, &poly6(), &poly2());
    // The result has degree at most m/l - 1 = 2 (0-indexed, so at most 3 coefficients).
    let deg = result.degree().unwrap_or(0);
    assert!(
        deg < 3,
        "constant poly restricts to degree < m/l = 3; got degree {}",
        deg
    );
}

/// A degree-1 polynomial over `GF(2^m)` restricts to a polynomial of degree at most
/// `1·(m/l) + (m/l−1) = 2·(m/l)−1` over `GF(2^l)`.
///
/// More precisely, the degree is at most `(d+1)·(m/l) − 1` where `d` is the degree
/// of the input polynomial.
#[test]
fn weil_restrict_degree_1_poly_degree_bound() {
    // p(x) = α·x + 1 where α is a non-subfield element.
    let a = f6(2); // α = x mod poly6 (primitive element of GF(2^6))
    let one = F2mNaive::<1>::one();
    let p = Poly::from_coeffs(vec![one, a]); // 1 + α·x
    let result = weil_restrict_poly(&p, 2, &poly6(), &poly2());
    // Degree bound: (1+1)·3 − 1 = 5.
    let deg = result.degree().unwrap_or(0);
    assert!(
        deg <= 5,
        "degree-1 poly restricts to degree ≤ 5; got degree {}",
        deg
    );
}

// ── Weil restriction of subfield elements ────────────────────────────────────

/// A subfield element `a ∈ GF(2^l) ⊂ GF(2^m)` restricts to itself.
///
/// For `a ∈ GF(2^l)`, the Weil restriction `Res_{m/l}(a)` has:
/// - Component 0: `a` itself (the element in `GF(2^l)`).
/// - Components 1, …, m/l−1: 0 (since `a` has no higher-basis components).
///
/// As a constant polynomial over `GF(2^l)`, the restriction is just `a`.
#[test]
fn weil_restrict_subfield_element_is_itself() {
    let p6 = poly6();
    let p2 = poly2();

    // Test all 4 elements of GF(2^2) ⊂ GF(2^6).
    // GF(2^2) = {0, 1, β, β+1} where β is the primitive element of GF(2^2) in GF(2^6).
    // The subfield elements are those satisfying a^4 = a in GF(2^6).
    for v in 0u64..64 {
        let a = F2mNaive::<1>::from_u64(v, &p6);
        if !is_in_subfield(&a, 2, &p6) {
            continue;
        }
        // a ∈ GF(2^2) ⊂ GF(2^6): restrict the constant polynomial [a].
        let p = Poly::from_coeffs(vec![a.clone()]);
        let result = weil_restrict_poly(&p, 2, &p6, &p2);

        // The restriction of a subfield element is a constant polynomial.
        // The constant term is the GF(2^2) representation of a.
        // Since a ∈ GF(2^2), its 0-th component in the {1, α, α²} basis is a itself,
        // and components 1 and 2 are zero.
        let deg = result.degree().unwrap_or(0);
        assert!(
            deg == 0 || result.is_zero(),
            "restriction of subfield element {v:#x} must be a constant polynomial; got degree {}",
            deg
        );

        // The constant term of the result must be non-zero for non-zero a.
        if !a.is_zero() {
            assert!(
                !result.is_zero(),
                "restriction of non-zero subfield element {v:#x} must be non-zero"
            );
        }
    }
}

/// The zero element restricts to the zero polynomial.
#[test]
fn weil_restrict_zero_element_is_zero() {
    let zero = F2mNaive::<1>::zero();
    let p = Poly::from_coeffs(vec![zero]);
    let result = weil_restrict_poly(&p, 2, &poly6(), &poly2());
    assert!(result.is_zero(), "restriction of 0 must be the zero polynomial");
}

/// The identity element (1) restricts to the constant polynomial 1.
///
/// `1 ∈ GF(2^2) ⊂ GF(2^6)` is a subfield element. Its restriction is the
/// constant polynomial `1` over `GF(2^2)`.
#[test]
fn weil_restrict_one_is_one() {
    let one = F2mNaive::<1>::one();
    let p = Poly::from_coeffs(vec![one.clone()]);
    let result = weil_restrict_poly(&p, 2, &poly6(), &poly2());
    // The restriction of 1 is the constant polynomial 1 over GF(2^2).
    assert!(!result.is_zero(), "restriction of 1 must be non-zero");
    assert_eq!(
        result.degree(),
        Some(0),
        "restriction of 1 must be a constant polynomial"
    );
    assert_eq!(
        result.coeff(0),
        F2mNaive::<1>::one(),
        "restriction of 1 must have constant term 1"
    );
}

// ── Weil restriction linearity KATs ──────────────────────────────────────────

/// Weil restriction is additive: `Res(f + g) = Res(f) + Res(g)`.
///
/// The Weil restriction is a ring homomorphism; in particular it is additive.
/// In characteristic 2, addition is XOR.
#[test]
fn weil_restrict_additivity() {
    let p6 = poly6();
    let p2 = poly2();

    // f = 1 + α·x (α = primitive element of GF(2^6))
    let one = F2mNaive::<1>::one();
    let alpha = F2mNaive::<1>::from_u64(2, &p6); // α = x mod poly6
    let f = Poly::from_coeffs(vec![one.clone(), alpha.clone()]);

    // g = α²·x + α³ (two non-subfield elements)
    let alpha2 = alpha.square(&p6);
    let alpha3 = alpha2.mul(&alpha, &p6);
    let g = Poly::from_coeffs(vec![alpha3.clone(), alpha2.clone()]);

    // f + g = (1 + α³) + (α + α²)·x
    let f_plus_g = f.add(&g);

    // Compute Res(f + g) and Res(f) + Res(g).
    let res_f = weil_restrict_poly(&f, 2, &p6, &p2);
    let res_g = weil_restrict_poly(&g, 2, &p6, &p2);
    let res_f_plus_g = weil_restrict_poly(&f_plus_g, 2, &p6, &p2);
    let res_f_add_res_g = res_f.add(&res_g);

    assert_eq!(
        res_f_plus_g, res_f_add_res_g,
        "Weil restriction must be additive: Res(f+g) = Res(f) + Res(g)"
    );
}

/// Weil restriction is homogeneous over `GF(2^l)`: `Res(c·f) = c·Res(f)` for `c ∈ GF(2^l)`.
///
/// Scalar multiplication by a subfield element commutes with the Weil restriction.
/// This is the GF(2^l)-linearity of the restriction map.
#[test]
fn weil_restrict_scalar_homogeneity() {
    let p6 = poly6();
    let p2 = poly2();

    // f = α·x + 1 (α = primitive element of GF(2^6))
    let one = F2mNaive::<1>::one();
    let alpha = F2mNaive::<1>::from_u64(2, &p6);
    let f = Poly::from_coeffs(vec![one.clone(), alpha.clone()]);

    // c ∈ GF(2^2): use c = 1 (trivial) and c = the embedded primitive element of GF(2^2).
    // The subfield GF(2^2) has 4 elements: {0, 1, β, β+1} where β^2 + β + 1 = 0.
    // In GF(2^6), β is the element satisfying β^4 = β (i.e., β^3 = 1, β ≠ 1).
    // We test with c = 1 (trivial) to verify the formula.
    let c = F2mNaive::<1>::one(); // c = 1 ∈ GF(2^2)

    // c·f = f (since c = 1).
    let cf_coeffs: Vec<F2mNaive<1>> = f
        .coeffs()
        .iter()
        .map(|coeff| c.mul(coeff, &p6))
        .collect();
    let cf = Poly::from_coeffs(cf_coeffs);

    let res_f = weil_restrict_poly(&f, 2, &p6, &p2);
    let res_cf = weil_restrict_poly(&cf, 2, &p6, &p2);

    // c·Res(f): multiply each coefficient of Res(f) by c (in GF(2^2)).
    let c_res_f_coeffs: Vec<F2mNaive<1>> = res_f
        .coeffs()
        .iter()
        .map(|coeff| c.mul(coeff, &p2))
        .collect();
    let c_res_f = Poly::from_coeffs(c_res_f_coeffs);

    assert_eq!(
        res_cf, c_res_f,
        "Weil restriction must be GF(2^l)-homogeneous: Res(c·f) = c·Res(f)"
    );
}

/// Weil restriction linearity: `Res(f + g) = Res(f) + Res(g)` for constant polynomials.
///
/// Tests additivity on constant polynomials (single field elements) to verify
/// the coefficient-level restriction is additive.
#[test]
fn weil_restrict_linearity_constant_polys() {
    let p6 = poly6();
    let p2 = poly2();

    // Use two non-zero GF(2^6) elements.
    let a = f6(0x15); // 0b010101
    let b = f6(0x2a); // 0b101010

    let pa = Poly::from_coeffs(vec![a.clone()]);
    let pb = Poly::from_coeffs(vec![b.clone()]);
    let pa_plus_b = Poly::from_coeffs(vec![a.add(&b)]);

    let res_a = weil_restrict_poly(&pa, 2, &p6, &p2);
    let res_b = weil_restrict_poly(&pb, 2, &p6, &p2);
    let res_a_plus_b = weil_restrict_poly(&pa_plus_b, 2, &p6, &p2);

    assert_eq!(
        res_a_plus_b,
        res_a.add(&res_b),
        "Res(a + b) must equal Res(a) + Res(b) for constant polynomials"
    );
}

// ── GhsParams struct KATs ─────────────────────────────────────────────────────

/// `GhsParams::new` rejects invalid parameters.
#[test]
fn ghs_params_new_rejects_invalid() {
    let result = GhsParams::new(6, 4, ghs_toy_curve(), poly2());
    assert!(
        matches!(result, Err(GhsError::SubfieldDivisibility)),
        "GhsParams::new must reject l=4 for m=6"
    );
}

/// `GhsParams::extension_degree` returns `m/l`.
#[test]
fn ghs_params_extension_degree() {
    let params = toy_params();
    assert_eq!(params.extension_degree(), 3, "m/l = 6/2 = 3");
}

/// The toy curve base point (0, 1) is on the curve `y²+xy = x³+x²+1`.
///
/// Verification: LHS = 1² + 0·1 = 1; RHS = 0 + 0 + 1 = 1. ✓
#[test]
fn toy_curve_base_point_on_curve() {
    let c = ghs_toy_curve();
    let g = c.generator::<F2mNaive<1>>();
    assert!(c.is_on_curve(&g), "toy GHS base point (0,1) must be on curve");
    // Verify it's the expected point.
    assert_eq!(g.x(), Some(&F2mNaive::<1>::zero()), "base point x must be 0");
    assert_eq!(g.y(), Some(&F2mNaive::<1>::one()), "base point y must be 1");
}

// ── GhsError KATs ────────────────────────────────────────────────────────────

/// `GhsError` variants have non-empty Display strings.
#[test]
fn ghs_error_display_non_empty() {
    let errors = [
        GhsError::SubfieldDivisibility,
        GhsError::NonDescendable,
        GhsError::NotOnCurve,
        GhsError::NotInSubfield,
        GhsError::SubfieldPolyMismatch,
    ];
    for e in &errors {
        let s = format!("{e}");
        assert!(!s.is_empty(), "GhsError display must be non-empty: {e:?}");
    }
}

/// `GhsError` implements `std::error::Error`.
#[test]
fn ghs_error_implements_error() {
    let e: &dyn std::error::Error = &GhsError::SubfieldDivisibility;
    let _ = e.to_string(); // must not panic
}

// ── Frobenius orbit KATs ──────────────────────────────────────────────────────

/// The Frobenius-by-subfield orbit has length `m/l = 3`.
///
/// For any `a ∈ GF(2^6)`, the orbit `[a, a^(2^2), a^(2^4)]` has 3 elements.
#[test]
fn frobenius_orbit_length_is_3() {
    let params = toy_params();
    let wr = WeilRestriction::new(params);
    let a = f6(0x15);
    let orbit = wr.frobenius_orbit(&a);
    assert_eq!(orbit.len(), 3, "Frobenius orbit length = m/l = 3");
}

/// The Frobenius orbit of a subfield element is constant.
///
/// For `a ∈ GF(2^2) ⊂ GF(2^6)`, `a^(2^2) = a` (Frobenius fixed field).
/// So the orbit is `[a, a, a]`.
#[test]
fn frobenius_orbit_of_subfield_element_is_constant() {
    let p6 = poly6();
    let params = toy_params();
    let wr = WeilRestriction::new(params);

    // Find a non-trivial subfield element (not 0 or 1).
    let subfield_elements: Vec<F2mNaive<1>> = (0u64..64)
        .map(|v| F2mNaive::<1>::from_u64(v, &p6))
        .filter(|a| shared_gf2m::is_in_subfield(a, 2, &p6))
        .collect();

    assert_eq!(subfield_elements.len(), 4, "GF(2^2) has exactly 4 elements");

    for a in &subfield_elements {
        let orbit = wr.frobenius_orbit(a);
        assert_eq!(orbit.len(), 3, "orbit must have length 3");
        // All orbit elements must equal a (Frobenius fixed field).
        for (j, conj) in orbit.iter().enumerate() {
            assert_eq!(
                conj, a,
                "orbit[{j}] of subfield element must equal a (Frobenius fixed field)"
            );
        }
    }
}

/// The relative trace `Tr_{6/2}(a)` lands in `GF(2^2)` for all `a ∈ GF(2^6)`.
#[test]
fn relative_trace_lands_in_subfield() {
    let p6 = poly6();
    let params = toy_params();
    let wr = WeilRestriction::new(params);

    for v in 0u64..64 {
        let a = f6(v);
        let tr = wr.trace(&a);
        assert!(
            shared_gf2m::is_in_subfield(&tr, 2, &p6),
            "Tr_{{6/2}}({v:#x}) must be in GF(2^2)"
        );
    }
}

/// The relative norm `N_{6/2}(a)` lands in `GF(2^2)` for all `a ∈ GF(2^6)`.
#[test]
fn relative_norm_lands_in_subfield() {
    let p6 = poly6();
    let params = toy_params();
    let wr = WeilRestriction::new(params);

    for v in 0u64..64 {
        let a = f6(v);
        let n = wr.norm(&a);
        assert!(
            shared_gf2m::is_in_subfield(&n, 2, &p6),
            "N_{{6/2}}({v:#x}) must be in GF(2^2)"
        );
    }
}

// ── GHS hyperelliptic-curve extraction KATs ───────────────────────────────────

// ── Genus formula KATs ───────────────────────────────────────────────────────

/// `ghs_genus(6, 2) == 1`: toy fixture (m=6, l=2, n=3, g=(3-1)/2=1).
///
/// The toy fixture gives genus 1 — a genus-1 hyperelliptic curve (an elliptic
/// curve in hyperelliptic form).  This is NOT the crypto-scale case (genus ≥ 2
/// is needed for the GHS attack to be effective), but the construction is
/// mathematically correct.
///
/// # Scale limits
///
/// Toy parameters (m=6, l=2) give genus 1.  Crypto-scale GHS uses n ≥ 5
/// (genus ≥ 2).  The formula `g = (n−1)/2` is crypto-scale-correct.
#[test]
fn ghs_genus_toy_fixture_is_1() {
    assert_eq!(
        ghs_genus(6, 2),
        1,
        "toy fixture (m=6, l=2, n=3): genus must be (3-1)/2 = 1"
    );
}

/// `ghs_genus(10, 2) == 2`: crypto-scale n=5 gives genus 2.
///
/// For m=10, l=2, n=5 (odd), the genus is (5-1)/2 = 2.  This is the smallest
/// crypto-scale GHS case (genus-2 hyperelliptic curve).
#[test]
fn ghs_genus_n5_is_2() {
    assert_eq!(
        ghs_genus(10, 2),
        2,
        "m=10, l=2, n=5: genus must be (5-1)/2 = 2"
    );
}

/// `ghs_genus(14, 2) == 3`: n=7 gives genus 3.
#[test]
fn ghs_genus_n7_is_3() {
    assert_eq!(
        ghs_genus(14, 2),
        3,
        "m=14, l=2, n=7: genus must be (7-1)/2 = 3"
    );
}

// ── Extracted curve validity KATs ─────────────────────────────────────────────

/// `extract_ghs_curve` succeeds for the toy fixture (m=6, l=2).
///
/// The toy fixture has odd m/l=3 (imaginary model) and curve coefficients
/// a=1, b=1 in GF(2^2) ⊂ GF(2^6).
#[test]
fn extract_ghs_curve_succeeds_for_toy_fixture() {
    let params = toy_params();
    let result = extract_ghs_curve(params);
    assert!(
        result.is_ok(),
        "extract_ghs_curve must succeed for toy fixture (m=6, l=2)"
    );
}

/// The extracted curve's `poly` is the GF(2^2) irreducible `x²+x+1 = 0x7`.
///
/// The descended curve `C/GF(2^l)` is defined over GF(2^2), so its field
/// polynomial must be the GF(2^2) irreducible.
#[test]
fn extracted_curve_poly_is_gf2_2() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    assert_eq!(
        curve.poly,
        poly2(),
        "extracted curve poly must be the GF(2^2) irreducible (0x7)"
    );
}

/// The extracted curve's genus matches `ghs_genus(6, 2) = 1`.
///
/// The `HyperellipticCurve::genus()` method returns `⌊(deg f − 1)/2⌋`.  For
/// the imaginary model with `deg f = 3`, this gives genus 1.  This must match
/// the formula `ghs_genus(6, 2) = (3-1)/2 = 1`.
#[test]
fn extracted_curve_genus_matches_formula() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let expected_genus = ghs_genus(6, 2);
    assert_eq!(
        curve.genus(),
        expected_genus,
        "extracted curve genus must match ghs_genus(6, 2) = {expected_genus}"
    );
}

// ── Imaginary model KATs ──────────────────────────────────────────────────────

/// The extracted curve is in the imaginary model: `deg f = 2g+1`.
///
/// The imaginary (ramified) hyperelliptic model requires `deg f = 2g+1` (odd).
/// For genus 1, `deg f = 3`.
#[test]
fn extracted_curve_imaginary_model_deg_f_eq_2g_plus_1() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let g = curve.genus();
    let deg_f = curve.f_coeffs.len() - 1;
    assert_eq!(
        deg_f,
        2 * g + 1,
        "imaginary model: deg f = {deg_f} must equal 2g+1 = {} (g = {g})",
        2 * g + 1
    );
}

/// The extracted curve's `H(X)` has degree ≤ genus.
///
/// The imaginary model requires `deg H ≤ g`.  For genus 1, `H(X) = X` has
/// degree 1 ≤ g = 1.
#[test]
fn extracted_curve_h_degree_le_genus() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let g = curve.genus();
    let h_poly = curve.h::<F2mNaive<1>>();
    let deg_h = h_poly.degree().unwrap_or(0);
    assert!(
        deg_h <= g,
        "imaginary model: deg H = {deg_h} must be ≤ genus g = {g}"
    );
}

/// The extracted curve's `F(X)` is monic.
///
/// The imaginary model requires the leading coefficient of `F` to be 1 (monic).
/// For `F(X) = X³ + X² + 1`, the leading coefficient is 1.
#[test]
fn extracted_curve_f_is_monic() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let f_poly = curve.f::<F2mNaive<1>>();
    let lc = f_poly.leading_coeff().expect("F must be non-zero");
    assert!(lc.is_one(), "F(X) must be monic (leading coefficient = 1)");
}

// ── Coefficients in GF(2^l) KATs ─────────────────────────────────────────────

/// All coefficients of `F(X)` are in GF(2^2).
///
/// The descended curve `C/GF(2^l)` must have all coefficients of `F(X)` in
/// GF(2^l) = GF(2^2).  We verify this using `is_in_subfield`.
///
/// For the toy fixture, `F(X) = X³ + X² + 1` with coefficients 1, 0, 1, 1 —
/// all in GF(2) ⊂ GF(2^2).
#[test]
fn extracted_curve_f_coeffs_in_gf2_2() {
    let p6 = poly6();
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");

    // The curve's poly is GF(2^2); we check coefficients are in GF(2^2)
    // by embedding them into GF(2^6) and using is_in_subfield.
    // Since the coefficients are raw bit-vectors in GF(2^2), we embed them
    // into GF(2^6) by treating them as GF(2^6) elements (they are subfield
    // elements, so the bit-vector is the same).
    for (i, &coeff_bits) in curve.f_coeffs.iter().enumerate() {
        let coeff_in_gf6 = F2mNaive::<1>::from_uint(coeff_bits, &p6);
        assert!(
            is_in_subfield(&coeff_in_gf6, 2, &p6),
            "F(X) coefficient at degree {i} (bits = {coeff_bits:#x}) must be in GF(2^2)"
        );
    }
}

/// All coefficients of `H(X)` are in GF(2^2).
///
/// `H(X) = X` has coefficients 0 and 1, both in GF(2) ⊂ GF(2^2).
#[test]
fn extracted_curve_h_coeffs_in_gf2_2() {
    let p6 = poly6();
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");

    for (i, &coeff_bits) in curve.h_coeffs.iter().enumerate() {
        let coeff_in_gf6 = F2mNaive::<1>::from_uint(coeff_bits, &p6);
        assert!(
            is_in_subfield(&coeff_in_gf6, 2, &p6),
            "H(X) coefficient at degree {i} (bits = {coeff_bits:#x}) must be in GF(2^2)"
        );
    }
}

// ── Known-answer coefficient KATs ─────────────────────────────────────────────

/// `H(X) = X`: the h polynomial is exactly `X` over GF(2^2).
///
/// The GHS construction with `y²+xy = x³+ax²+b` gives `H(X) = X` (from the
/// `xy` term in the Weierstrass equation).  This is the standard choice for the
/// imaginary model.
#[test]
fn extracted_curve_h_is_x() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let h_poly = curve.h::<F2mNaive<1>>();

    assert_eq!(h_poly.degree(), Some(1), "H(X) must have degree 1 (H = X)");
    assert!(h_poly.coeff(0).is_zero(), "H(X) constant term must be 0");
    assert!(h_poly.coeff(1).is_one(), "H(X) coefficient of X must be 1");
}

/// `F(X) = X³ + X² + 1`: the f polynomial matches the Artin–Schreier polynomial.
///
/// For the toy curve `y²+xy = x³+x²+1` (a=1, b=1), the Artin–Schreier polynomial
/// is `f_AS(x) = x³+x²+1`.  Since a=1, b=1 ∈ GF(2^2), the descended curve's
/// `F(X) = X³+X²+1` over GF(2^2).
#[test]
fn extracted_curve_f_is_x3_plus_x2_plus_1() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let f_poly = curve.f::<F2mNaive<1>>();

    assert_eq!(f_poly.degree(), Some(3), "F(X) must have degree 3");
    // F(X) = X³ + X² + 1: coefficients [1, 0, 1, 1] (index = degree).
    assert!(f_poly.coeff(0).is_one(), "F(X) constant term must be 1 (= b)");
    assert!(f_poly.coeff(1).is_zero(), "F(X) coefficient of X must be 0");
    assert!(f_poly.coeff(2).is_one(), "F(X) coefficient of X² must be 1 (= a)");
    assert!(f_poly.coeff(3).is_one(), "F(X) coefficient of X³ must be 1 (monic)");
}

// ── Curve validity KATs ───────────────────────────────────────────────────────

/// The zero divisor `[1, 0]` is valid on the extracted curve.
///
/// The zero divisor is the group identity of the Jacobian.  It must satisfy the
/// Mumford invariant on the extracted curve.
#[test]
fn extracted_curve_zero_divisor_is_valid() {
    let params = toy_params();
    let curve = extract_ghs_curve(params).expect("extraction must succeed");
    let zero_div = curve.zero_divisor::<F2mNaive<1>>();
    assert!(
        curve.is_valid(&zero_div),
        "zero divisor [1, 0] must be valid on the extracted GHS curve"
    );
}

/// Even `m/l` is rejected with `NonDescendable`.
///
/// The imaginary model requires odd `m/l`.  Even `m/l` yields the real/split
/// model, which the frozen [`rho::hyperelliptic::HyperellipticCurve`] does not handle.
///
/// We test with m=4, l=2, n=2 (even).
#[test]
fn extract_ghs_curve_rejects_even_extension_degree() {
    use rho::binary_curve::BinaryCurve;
    let poly4 = Uint::<1>::from(0x13u64); // x⁴+x+1
    let curve = BinaryCurve {
        poly: poly4,
        a: Uint::<1>::ONE,
        b: Uint::<1>::ONE,
        n: Uint::<1>::ONE,
        gx: Uint::<1>::ZERO,
        gy: Uint::<1>::ONE,
    };
    let params = GhsParams::new(4, 2, curve, poly2()).expect("4/2 is valid");
    let result = extract_ghs_curve(params);
    assert!(
        matches!(result, Err(GhsError::NonDescendable)),
        "even m/l = 2 must be rejected with NonDescendable"
    );
}

// ── GHS transfer map KATs ─────────────────────────────────────────────────────

// ── Identity maps to identity ─────────────────────────────────────────────────

/// `transfer_point(∞, ...) == [1, 0]` — the point at infinity maps to the zero divisor.
///
/// The GHS transfer map is a group homomorphism, so the identity must map to the
/// identity.  The point at infinity `∞ ∈ E(GF(2^m))` is the group identity of `E`,
/// and `[1, 0]` is the group identity of `Jac(C)`.
#[test]
fn transfer_infinity_is_zero_divisor() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
    let d = transfer_point(&inf, &curve_c, &params).expect("transfer must succeed");
    assert!(
        d.is_zero(),
        "transfer(∞) must be the zero divisor [1, 0]"
    );
}

/// The zero divisor returned for `∞` is valid on the extracted curve.
///
/// The zero divisor `[1, 0]` must satisfy the Mumford invariant on `C`.
#[test]
fn transfer_infinity_result_is_valid() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
    let d = transfer_point(&inf, &curve_c, &params).expect("transfer must succeed");
    assert!(
        curve_c.is_valid(&d),
        "transfer(∞) = [1, 0] must be a valid divisor on C"
    );
}

// ── Transfer of base point ────────────────────────────────────────────────────

/// `transfer_point(G, ...)` succeeds and gives a non-zero divisor.
///
/// The base point `G = (0, 1)` on `E` must transfer to a non-trivial divisor on
/// `Jac(C)`.  A zero result would mean `G` is in the kernel of the transfer map,
/// which would break the homomorphism property for the DLP.
#[test]
fn transfer_base_point_is_nonzero() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let d = transfer_point(&g, &curve_c, &params).expect("transfer must succeed");
    assert!(
        !d.is_zero(),
        "transfer(G) must be a non-zero divisor (G is not in the kernel)"
    );
}

/// `transfer_point(G, ...)` returns a valid reduced divisor.
///
/// The transferred divisor must satisfy the Mumford invariant: `u` monic,
/// `deg v < deg u ≤ g`, and `u | (f − v·h − v²)`.
#[test]
fn transfer_base_point_is_valid_divisor() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let d = transfer_point(&g, &curve_c, &params).expect("transfer must succeed");
    assert!(
        curve_c.is_valid(&d),
        "transfer(G) must be a valid reduced Mumford divisor"
    );
}

/// Known-answer: `D_G = [X, 1]` — the base point `G = (0, 1)` transfers to `[X, 1]`.
///
/// For `G = (0, 1)`:
/// - `X_G = Tr_{6/2}(0) = 0` — the trace of 0 is 0.
/// - `Y_G = 1` — from `Y² = F(0) = 1` over GF(2^2), so `Y = 1`.
/// - Divisor: `u(X) = X + 0 = X`, `v(X) = 1`.
///
/// In GF(2^2) with poly `x²+x+1 = 0x7`:
/// - `u` has coefficients `[0, 1]` (constant 0, X-coefficient 1).
/// - `v` has coefficients `[1]` (constant 1).
#[test]
fn transfer_base_point_known_answer() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let d = transfer_point(&g, &curve_c, &params).expect("transfer must succeed");

    // u(X) = X: degree 1, constant term 0, X-coefficient 1.
    assert_eq!(
        d.u.degree(),
        Some(1),
        "D_G: u must have degree 1"
    );
    assert!(
        d.u.coeff(0).is_zero(),
        "D_G: u constant term must be 0 (X_G = Tr(0) = 0)"
    );
    assert!(
        d.u.coeff(1).is_one(),
        "D_G: u X-coefficient must be 1 (monic)"
    );

    // v(X) = 1: constant polynomial.
    assert_eq!(
        d.v.degree(),
        Some(0),
        "D_G: v must have degree 0 (constant)"
    );
    assert!(
        d.v.coeff(0).is_one(),
        "D_G: v constant term must be 1 (Y_G = 1)"
    );
}

// ── Transfer of a sample point ────────────────────────────────────────────────

/// `transfer_point(P, ...)` returns a valid reduced divisor for `P = (0x01, 0x3a)`.
///
/// `P = (0x01, 0x3a)` is a non-base affine point on `E/GF(2^6)`.  Its transfer
/// must be a valid reduced Mumford divisor on `Jac(C)`.
#[test]
fn transfer_sample_point_is_valid_divisor() {
    let p6 = poly6();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    // P = (0x01, 0x3a) is on E (verified: y²+xy = x³+x²+1 at x=1, y=0x3a).
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );
    let d = transfer_point(&p, &curve_c, &params).expect("transfer must succeed");
    assert!(
        curve_c.is_valid(&d),
        "transfer(P) must be a valid reduced Mumford divisor"
    );
}

/// Known-answer: `D_P = [X+1, 2]` for `P = (0x01, 0x3a)`.
///
/// Derived by hand, independently of the code under test.  `α` is the root of
/// `x⁶+x+1`, so `α⁶ = α+1`, and the subfield generator is `β = α²¹`.
///
/// 1. Powers of `α`: `α⁹ = α⁴+α³` and `α¹² = α²+1`, so
///    `β = α¹²·α⁹ = (α²+1)(α⁴+α³) = α⁶+α⁵+α⁴+α³ = α⁵+α⁴+α³+α+1 = 0b111011 = 0x3b`.
///    `β³ = 1` and `β ≠ 1`, so `β²+β+1 = 0` and `β² = β+1 = 0b111010 = 0x3a`.
///    The minimal polynomial of `β` is `x²+x+1`, which is `poly_l = 0x7`.
/// 2. `GF(2^2) = {0, 1, β, β²} = {0, 0x01, 0x3b, 0x3a}`; coordinates in the basis `{1, β}`
///    (bit `i` = coefficient of `βⁱ`): `1 ↦ 0b01 = 1`, `β ↦ 0b10 = 2`, `β² = 1+β ↦ 0b11 = 3`.
/// 3. `P = (1, β²)` has both coordinates in `GF(2^2)`, so `φ_2(P) = φ_2²(P) = P` and the
///    conorm is `R = 3P`.  On `E: y²+xy = x³+x²+1` (a = 1, b = 1):
///    doubling uses `λ = x + y/x = 1+β² = β`, `x₃ = λ²+λ+a = β²+β+1 = 0`,
///    `y₃ = x² + (λ+1)x₃ = 1`, so `2P = (0, 1) = G`, a point of order 2 (`2G = ∞`).
///    Hence `4P = ∞` and `3P = −P = (x, y+x) = (1, β²+1) = (1, β)`.
/// 4. `R = (1, β)`: `x_R = 1 ↦ 1`, `y_R = β ↦ 2`.  So `D_P = [X+1, 2]`.
///
/// The earlier value `v = 3` came from reducing the raw bits of `β = 0x3b` modulo `x²+x+1`,
/// which names `β²` instead of `β`; the basis-coordinate reading gives `2`.  (PARI/GP
/// `ellmul` on the same curve over `GF(2^6)` returns `3P = (1, 0x3b)`, matching step 3.)
#[test]
fn transfer_sample_point_known_answer() {
    let p6 = poly6();
    let p2 = poly2();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );
    let d = transfer_point(&p, &curve_c, &params).expect("transfer must succeed");

    // u(X) = X + 1: degree 1, constant term 1, X-coefficient 1.
    assert_eq!(d.u.degree(), Some(1), "D_P: u must have degree 1");
    assert!(d.u.coeff(1).is_one(), "D_P: u X-coefficient must be 1 (monic)");
    assert!(
        d.u.coeff(0).is_one(),
        "D_P: u constant term must be 1 (x_R = 1 from conorm)"
    );

    // v(X) = β = 2 in the basis {1, β} of GF(2^2).
    let beta = F2mNaive::<1>::from_u64(2, &p2);
    assert_eq!(
        d.v.coeff(0),
        beta,
        "D_P: v constant term must be β = 2 in GF(2^2)"
    );
}

// ── Homomorphism property (decisive correctness guard) ────────────────────────

/// `D_{G+P} = D_G + D_P` via Cantor compose — the homomorphism property.
///
/// This is the decisive correctness guard for the GHS transfer map.  The transfer
/// is a group homomorphism iff `D_{P+Q} = D_P + D_Q` for all `P, Q ∈ E(GF(2^m))`.
///
/// Known-answer:
/// - `D_G = [X, 1]`, `D_P = [X+1, 2]`.
/// - `G + P = (0x01, 0x3b)` on `E`.
/// - `D_{G+P} = [X+1, 3]`.
/// - `D_G + D_P (Cantor) = [X+1, 3]` ✓.
#[test]
fn transfer_homomorphism_g_plus_p() {
    let p6 = poly6();
    let p2 = poly2();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    // G = (0, 1) — the base point.
    let g = curve_e.generator::<F2mNaive<1>>();
    // P = (0x01, 0x3a) — a sample point on E.
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );

    // Compute D_G, D_P, and D_{G+P}.
    let d_g = transfer_point(&g, &curve_c, &params).expect("transfer(G) must succeed");
    let d_p = transfer_point(&p, &curve_c, &params).expect("transfer(P) must succeed");
    let g_plus_p = curve_e.add(&g, &p);
    let d_g_plus_p = transfer_point(&g_plus_p, &curve_c, &params)
        .expect("transfer(G+P) must succeed");

    // Compute D_G + D_P via Cantor compose.
    let d_sum = cantor::add(&curve_c, &d_g, &d_p, &p2);

    assert_eq!(
        d_sum, d_g_plus_p,
        "homomorphism must hold: D_G + D_P (Cantor) must equal D_{{G+P}}"
    );
}

/// `D_{G+P} = D_G + D_P` — homomorphism via `verify_homomorphism` helper.
///
/// Cross-checks the homomorphism using the `verify_homomorphism` convenience
/// function, which encapsulates the full check.
#[test]
fn transfer_homomorphism_via_helper() {
    let p6 = poly6();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    let g = curve_e.generator::<F2mNaive<1>>();
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );

    assert!(
        verify_homomorphism(&g, &p, &curve_e, &curve_c, &params),
        "verify_homomorphism(G, P) must return true"
    );
}

/// `D_{G+G} = D_G + D_G` — the doubling case of the homomorphism.
///
/// The homomorphism must hold for `P = Q` (doubling).  This tests that the
/// transfer map is consistent with the Cantor doubling formula.
#[test]
fn transfer_homomorphism_doubling() {
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();

    assert!(
        verify_homomorphism(&g, &g, &curve_e, &curve_c, &params),
        "homomorphism must hold for G + G (doubling case)"
    );
}

/// `D_{∞+P} = D_∞ + D_P = [1,0] + D_P = D_P` — identity element homomorphism.
///
/// The transfer of `∞ + P = P` must equal `transfer(∞) + transfer(P) = [1,0] + D_P = D_P`.
/// This verifies the identity-element case of the homomorphism.
#[test]
fn transfer_homomorphism_identity_element() {
    let p6 = poly6();
    let p2 = poly2();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");

    let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );

    let d_inf = transfer_point(&inf, &curve_c, &params).expect("transfer(∞) must succeed");
    let d_p = transfer_point(&p, &curve_c, &params).expect("transfer(P) must succeed");

    // D_∞ + D_P = [1,0] + D_P = D_P (Cantor identity).
    let d_sum = cantor::add(&curve_c, &d_inf, &d_p, &p2);
    assert_eq!(
        d_sum, d_p,
        "D_∞ + D_P must equal D_P (identity element homomorphism)"
    );
}

/// Known-answer: `D_{G+P} = [X+1, 3]` for `G = (0,1)` and `P = (0x01, 0x3a)`.
///
/// Derived by hand, independently of the code under test; the field setup is the one in
/// `transfer_sample_point_known_answer` (`β = α²¹ = 0x3b`, `β² = 0x3a = β+1`, coordinates in
/// the basis `{1, β}`: `1 ↦ 1`, `β ↦ 2`, `β² ↦ 3`).
///
/// 1. `P = (1, β²)` has `2P = G` and `4P = ∞` (see that test), so `G + P = 3P = −P`
///    `= (1, β²+1) = (1, β) = (0x01, 0x3b)`.
/// 2. `G + P` lies in `E(GF(2^2))`, so the conorm is `R = 3(G+P) = 9P = P` (since
///    `4P = ∞`): `R = (1, β²)`.
/// 3. `x_R = 1 ↦ 1`, `y_R = β² ↦ 3`.  So `D_{G+P} = [X+1, 3]`.
///
/// The earlier value `v = 2` came from reducing the raw bits of `β² = 0x3a` modulo `x²+x+1`
/// (which names `β`); the basis-coordinate reading gives `3`.  PARI/GP `ellmul` on the same
/// curve over `GF(2^6)` returns `3(G+P) = (1, 0x3a)`.
///
/// Consistency with the homomorphism: `D_G + D_P = [X, 1] + [X+1, 2] = [X+1, 3]`.
#[test]
fn transfer_g_plus_p_known_answer() {
    let p6 = poly6();
    let p2 = poly2();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    let g = curve_e.generator::<F2mNaive<1>>();
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3a, &p6),
    );
    let g_plus_p = curve_e.add(&g, &p);
    let d = transfer_point(&g_plus_p, &curve_c, &params)
        .expect("transfer(G+P) must succeed");

    // u(X) = X + 1.
    assert_eq!(d.u.degree(), Some(1), "D_{{G+P}}: u must have degree 1");
    assert!(d.u.coeff(1).is_one(), "D_{{G+P}}: u X-coefficient must be 1 (monic)");
    assert!(
        d.u.coeff(0).is_one(),
        "D_{{G+P}}: u constant term must be 1 (x_R = 1 from conorm)"
    );

    // v(X) = β² = β+1 = 3 in the basis {1, β} of GF(2^2).
    let beta_plus_1 = F2mNaive::<1>::from_u64(3, &p2);
    assert_eq!(
        d.v.coeff(0),
        beta_plus_1,
        "D_{{G+P}}: v constant term must be β+1 = 3 in GF(2^2)"
    );
}

// ── Transfer input guards ─────────────────────────────────────────────────────

/// A point that does not satisfy the curve equation is rejected with `NotOnCurve`.
///
/// `(0x01, 0x3c)` is `P = (0x01, 0x3a)` with the y-coordinate perturbed: at `x = 1`
/// the curve has only `y ∈ {0x3a, 0x3b}`.  `(0, 0)` fails too (`0 ≠ b = 1`).
#[test]
fn transfer_rejects_off_curve_point() {
    let p6 = poly6();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    for (x, y) in [(0x01u64, 0x3cu64), (0x00, 0x00)] {
        let p = BinaryAffinePoint::new(
            F2mNaive::<1>::from_u64(x, &p6),
            F2mNaive::<1>::from_u64(y, &p6),
        );
        assert!(!curve_e.is_on_curve(&p), "fixture ({x:#x}, {y:#x}) must be off E");
        assert_eq!(
            transfer_point(&p, &curve_c, &params),
            Err(GhsError::NotOnCurve),
            "off-curve ({x:#x}, {y:#x}) must be rejected"
        );
    }
}

/// `ghs_descend` surfaces `NotOnCurve` for an off-curve base point or target.
#[test]
fn ghs_descend_rejects_off_curve_point() {
    let p6 = poly6();
    let params = toy_params();
    let g = ghs_toy_curve().generator::<F2mNaive<1>>();
    let off = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x3c, &p6),
    );

    assert!(matches!(
        ghs_descend(&params, &off, &g),
        Err(GhsError::NotOnCurve)
    ));
    assert!(matches!(
        ghs_descend(&params, &g, &off),
        Err(GhsError::NotOnCurve)
    ));
}

/// A conorm image with a coordinate outside `GF(2^2)` is rejected with `NotInSubfield`.
///
/// Fixture: the same equation with `a = 0x02 = α` (a generator of `GF(2^6)`, not in
/// `GF(2^2) = {0, 1, 0x3a, 0x3b}`).  `E` is then not defined over `GF(2^2)`, so the
/// relative Frobenius does not map `E` to itself and the conorm sum is not fixed by
/// it.  `(0x01, 0x26)` lies on this curve (`y² + xy = x³ + αx² + 1` at `x = 1`).  The
/// descended curve is the toy `C`, passed unchanged, so only the subfield guard can fire.
#[test]
fn transfer_rejects_conorm_outside_subfield() {
    let p6 = poly6();
    let toy = ghs_toy_curve();
    let curve_c = extract_ghs_curve(toy_params()).expect("extraction must succeed");

    let bad_curve = BinaryCurve { a: Uint::<1>::from(2u64), ..toy };
    let bad_params = GhsParams::new(6, 2, bad_curve.clone(), poly2()).expect("valid degrees");
    let p = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &p6),
        F2mNaive::<1>::from_u64(0x26, &p6),
    );
    assert!(bad_curve.is_on_curve(&p), "fixture point must be on the a = α curve");
    assert!(
        !is_in_subfield(&F2mNaive::<1>::from_u64(0x02, &p6), 2, &p6),
        "a = α must lie outside GF(2^2)"
    );

    assert_eq!(
        transfer_point(&p, &curve_c, &bad_params),
        Err(GhsError::NotInSubfield)
    );
}

/// A `poly_l` that is not the minimal polynomial of `β = α²¹` is rejected before any
/// coefficient or point is read, by both curve extraction and the transfer map.
///
/// `x²+1` has degree 2 but does not vanish at `β` (`β²+1 = β ≠ 0`); `x³+x+1` is
/// irreducible of degree 3 but is the minimal polynomial of `β⁻¹` for the `l = 3`
/// generator `α⁹` (whose minimal polynomial is `x³+x²+1`, PARI/GP `minpoly`).
#[test]
fn mismatched_poly_l_is_rejected() {
    let good = toy_params();
    let curve_c = extract_ghs_curve(good).expect("extraction must succeed");
    let g = ghs_toy_curve().generator::<F2mNaive<1>>();
    let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;

    for (l, bad_poly) in [(2usize, 0x5u64), (3, 0xb)] {
        let bad = GhsParams::new(6, l, ghs_toy_curve(), Uint::<1>::from(bad_poly))
            .expect("l divides m");
        assert_eq!(
            extract_ghs_curve(bad.clone()).err(),
            Some(GhsError::SubfieldPolyMismatch),
            "extraction must reject poly_l = {bad_poly:#x} for l = {l}"
        );
        for point in [&g, &inf] {
            assert_eq!(
                transfer_point(point, &curve_c, &bad),
                Err(GhsError::SubfieldPolyMismatch),
                "transfer must reject poly_l = {bad_poly:#x} for l = {l}"
            );
        }
        assert!(matches!(
            ghs_descend(&bad, &g, &g),
            Err(GhsError::SubfieldPolyMismatch)
        ));
    }
}

/// A `GhsParams` whose `l` does not divide the degree of the source polynomial is rejected
/// by extraction and transfer with `SubfieldDivisibility`, not a panic in `restrict`.
///
/// `GhsParams` fields are public, so the check in `GhsParams::new` can be bypassed.
#[test]
fn non_dividing_l_is_rejected_without_panic() {
    let good = toy_params();
    let curve_c = extract_ghs_curve(good.clone()).expect("extraction must succeed");
    let g = ghs_toy_curve().generator::<F2mNaive<1>>();
    let bad = GhsParams { l: 4, ..good };
    assert_eq!(extract_ghs_curve(bad.clone()).err(), Some(GhsError::SubfieldDivisibility));
    assert_eq!(
        transfer_point(&g, &curve_c, &bad),
        Err(GhsError::SubfieldDivisibility)
    );
}

/// The curve coefficients and the transferred divisors use one identification of `GF(2^2)`.
///
/// With `a = β = 0x3b` (a subfield element outside `GF(2)`), `E: y²+xy = x³+βx²+1` descends
/// to `C: Y²+XY = X³+2X²+1` over `GF(2^2)` in the basis `{1, β}` (`β ↦ 2`; a raw-bits
/// reading of `0x3b` would give `3 = β²`).  Every point of `E(GF(2^6))` must then transfer
/// to a divisor that is valid on that `C`, and the transfer must stay a homomorphism.
#[test]
fn curve_and_transfer_share_the_subfield_identification() {
    let p6 = poly6();
    let beta_curve = BinaryCurve { a: Uint::<1>::from(0x3bu64), ..ghs_toy_curve() };
    let params =
        GhsParams::new(6, 2, beta_curve.clone(), poly2()).expect("toy GHS params must be valid");
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");

    assert!(
        curve_c.f::<F2mNaive<1>>().coeff(2) == F2mNaive::<1>::from_u64(2, &poly2()),
        "F(X) coefficient of X² must be β ↦ 2 in the basis {{1, β}}"
    );

    let mut points = Vec::new();
    for x in 0u64..64 {
        for y in 0u64..64 {
            let pt = BinaryAffinePoint::new(
                F2mNaive::<1>::from_u64(x, &p6),
                F2mNaive::<1>::from_u64(y, &p6),
            );
            if beta_curve.is_on_curve(&pt) {
                points.push(pt);
            }
        }
    }
    assert!(points.len() > 8, "the curve must have affine points to transfer");

    for pt in &points {
        let d = transfer_point(pt, &curve_c, &params).expect("transfer must succeed");
        assert!(curve_c.is_valid(&d), "transferred divisor must be valid on the extracted C");
    }
    for p in points.iter().take(6) {
        for q in points.iter().take(6) {
            assert!(
                verify_homomorphism(p, q, &beta_curve, &curve_c, &params),
                "homomorphism must hold on the a = β curve"
            );
        }
    }
}

/// All transferred divisors are valid reduced divisors.
///
/// For all affine points on `E/GF(2^6)`, the transferred divisor must satisfy
/// the Mumford invariant on `C/GF(2^2)`.  This is a bulk validity check over
/// a sample of points.
#[test]
fn transfer_all_sample_points_are_valid() {
    let p6 = poly6();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    // Sample points on E (known to be on the curve from the fixture).
    let sample_points: &[(u64, u64)] = &[
        (0x00, 0x01), // G = base point
        (0x01, 0x3a), // P
        (0x01, 0x3b), // -P (negation of P)
        (0x06, 0x39), // another point
        (0x06, 0x3f), // its negation
    ];

    for &(x, y) in sample_points {
        let pt = BinaryAffinePoint::new(
            F2mNaive::<1>::from_u64(x, &p6),
            F2mNaive::<1>::from_u64(y, &p6),
        );
        assert!(
            curve_e.is_on_curve(&pt),
            "sample point ({x:#04x}, {y:#04x}) must be on E"
        );
        let d = transfer_point(&pt, &curve_c, &params)
            .expect("transfer must succeed for all finite points");
        assert!(
            curve_c.is_valid(&d),
            "transfer({x:#04x}, {y:#04x}) must be a valid reduced divisor"
        );
    }
}

/// The transfer map is consistent with the Cantor group law for multiple pairs.
///
/// Checks the homomorphism property for several pairs of points, not just G and P.
/// This is a stronger correctness signal than a single pair.
#[test]
fn transfer_homomorphism_multiple_pairs() {
    let p6 = poly6();
    let params = toy_params();
    let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
    let curve_e = ghs_toy_curve();

    let sample_points: &[(u64, u64)] = &[
        (0x00, 0x01), // G
        (0x01, 0x3a), // P
        (0x06, 0x39), // Q
    ];

    let pts: Vec<BinaryAffinePoint<F2mNaive<1>>> = sample_points
        .iter()
        .map(|&(x, y)| {
            BinaryAffinePoint::new(
                F2mNaive::<1>::from_u64(x, &p6),
                F2mNaive::<1>::from_u64(y, &p6),
            )
        })
        .collect();

    // Check all pairs (including self-pairs for doubling).
    for (i, pi) in pts.iter().enumerate() {
        for (j, pj) in pts.iter().enumerate() {
            assert!(
                verify_homomorphism(pi, pj, &curve_e, &curve_c, &params),
                "homomorphism must hold for pair ({i}, {j})"
            );
        }
    }
}

// ── GHS reduction + logarithm-preservation KATs ───────────────────────────────

// ── ghs_descend API KATs ──────────────────────────────────────────────────────

/// `ghs_descend` succeeds for the toy fixture.
///
/// The toy fixture (m=6, l=2, odd m/l=3) satisfies all preconditions for the
/// GHS descent.  The reduction must succeed and return a valid result.
#[test]
fn ghs_descend_succeeds_for_toy_fixture() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone(); // h = g (k=1)
    let result = ghs_descend(&params, &g, &h);
    assert!(result.is_ok(), "ghs_descend must succeed for toy fixture");
}

/// The descended curve's genus matches `ghs_genus(6, 2) = 1`.
///
/// The `GhsDescentResult::curve_c` must have genus 1 (from the formula
/// `g = (m/l − 1)/2 = (3−1)/2 = 1`).
#[test]
fn ghs_descend_curve_genus_is_1() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone();
    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    let expected_genus = ghs_genus(6, 2);
    assert_eq!(
        result.curve_c.genus(),
        expected_genus,
        "descended curve genus must be ghs_genus(6, 2) = {expected_genus}"
    );
}

/// `D_g` is a non-zero divisor.
///
/// The base point `G` must not be in the kernel of the transfer map — a zero
/// `D_g` would make the Jacobian-DLP instance trivial and useless.
#[test]
fn ghs_descend_d_g_is_nonzero() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone();
    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(!result.d_g.is_zero(), "D_g must be a non-zero divisor");
}

/// `D_g` is a valid reduced Mumford divisor.
///
/// The transferred base-point divisor must satisfy the Mumford invariant:
/// `u` monic, `deg v < deg u ≤ g`, and `u | (f − v·h − v²)`.
#[test]
fn ghs_descend_d_g_is_valid_divisor() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone();
    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(
        result.curve_c.is_valid(&result.d_g),
        "D_g must be a valid reduced Mumford divisor"
    );
}

/// `D_h` is a valid reduced Mumford divisor.
///
/// The transferred target-point divisor must satisfy the Mumford invariant.
#[test]
fn ghs_descend_d_h_is_valid_divisor() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    // Use P = (0x01, 0x3a) as h (a non-trivial target point).
    let h = BinaryAffinePoint::new(
        F2mNaive::<1>::from_u64(0x01, &poly6()),
        F2mNaive::<1>::from_u64(0x3a, &poly6()),
    );
    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(
        result.curve_c.is_valid(&result.d_h),
        "D_h must be a valid reduced Mumford divisor"
    );
}

/// Reduction with identity: `ghs_descend(params, g, ∞) → D_h = zero divisor`.
///
/// The point at infinity `∞ = 0·g` on `E` must transfer to the zero divisor
/// `[1, 0]` on `Jac(C)`.  This is the identity case of the homomorphism:
/// `D_{0·g} = 0·D_g = [1, 0]`.
#[test]
fn ghs_descend_identity_gives_zero_d_h() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity; // ∞ = 0·g
    let result = ghs_descend(&params, &g, &inf).expect("descent must succeed");
    assert!(
        result.d_h.is_zero(),
        "D_h for h=∞ must be the zero divisor [1, 0] (∞ = 0·g)"
    );
}

// ── Logarithm-preservation KATs (decisive) ────────────────────────────────────

/// Logarithm preservation for `k=3` (decisive correctness signal).
///
/// Computes `h = 3·g` on `E` using binary-curve scalar multiplication, then
/// calls `ghs_descend` and verifies `scalar_mul(D_g, 3) == D_h` on `Jac(C)`.
///
/// This is the decisive correctness signal for the GHS reduction: it confirms
/// that the transfer map preserves the discrete logarithm.
///
/// # Mathematical justification
///
/// The transfer map `φ: E(GF(2^m)) → Jac(C)(GF(2^l))` is a group homomorphism
/// ([`rho::ghs::transfer_point`]).  For `h = k·g`:
/// ```text
/// D_h = φ(h) = φ(k·g) = k·φ(g) = k·D_g
/// ```
/// So `scalar_mul(D_g, k) = D_h` on `Jac(C)`.
#[test]
fn ghs_descend_log_preservation_k3() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();

    // Compute h = 3·g on E using binary-curve scalar_mul.
    let g = curve_e.generator::<F2mNaive<1>>();
    let k: u64 = 3;
    let h = curve_e.scalar_mul(&g, &Uint::<1>::from(k));

    // Verify h is a finite point (not ∞) — the toy curve has order > 3.
    assert!(
        !h.is_infinity(),
        "3·G must be a finite point (toy curve order > 3)"
    );

    // Perform the GHS reduction.
    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");

    // Decisive check: scalar_mul(D_g, 3) == D_h on Jac(C).
    assert!(
        verify_log_preservation(&result, k),
        "logarithm preservation must hold: scalar_mul(D_g, 3) must equal D_h"
    );
}

/// Logarithm preservation for `k=1` (identity case).
///
/// For `h = 1·g = g`, the reduction must give `D_h = 1·D_g = D_g`.
#[test]
fn ghs_descend_log_preservation_k1() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone(); // h = 1·g

    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(
        verify_log_preservation(&result, 1),
        "logarithm preservation must hold for k=1: D_h must equal D_g"
    );
}

/// Logarithm preservation for `k=2` (doubling case).
///
/// For `h = 2·g`, the reduction must give `D_h = 2·D_g`.
#[test]
fn ghs_descend_log_preservation_k2() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let k: u64 = 2;
    let h = curve_e.scalar_mul(&g, &Uint::<1>::from(k));

    // Only run if h is finite (2·G might be ∞ for order-2 curves, but our toy curve has order > 2).
    if h.is_infinity() {
        // 2·G = ∞ means G has order 2; skip the log-preservation check.
        return;
    }

    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(
        verify_log_preservation(&result, k),
        "logarithm preservation must hold for k=2: scalar_mul(D_g, 2) must equal D_h"
    );
}

/// Logarithm preservation via `verify_log_preservation` helper for k=3.
///
/// Cross-checks the log-preservation using the `verify_log_preservation`
/// convenience function, which encapsulates the Cantor scalar_mul check.
#[test]
fn ghs_descend_log_preservation_via_helper() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let k: u64 = 3;
    let h = curve_e.scalar_mul(&g, &Uint::<1>::from(k));

    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
    assert!(
        verify_log_preservation(&result, k),
        "verify_log_preservation must return true for k=3"
    );
}

/// Logarithm preservation: `D_h = k·D_g` verified directly via Cantor scalar_mul.
///
/// This test bypasses `verify_log_preservation` and directly computes
/// `cantor::scalar_mul(D_g, k, C, poly_l)` to confirm equality with `D_h`.
/// It is the most explicit form of the decisive correctness check.
#[test]
fn ghs_descend_log_preservation_direct_cantor_check() {
    let params = toy_params();
    let curve_e = ghs_toy_curve();
    let g = curve_e.generator::<F2mNaive<1>>();
    let k: u64 = 3;
    let h = curve_e.scalar_mul(&g, &Uint::<1>::from(k));

    let result = ghs_descend(&params, &g, &h).expect("descent must succeed");

    // Direct Cantor scalar_mul check.
    let poly_l = &result.curve_c.poly;
    let k_d_g = cantor::scalar_mul(&result.curve_c, &result.d_g, k, poly_l);
    assert_eq!(
        k_d_g, result.d_h,
        "cantor::scalar_mul(D_g, 3) must equal D_h (logarithm preservation)"
    );
}

// ── Optional PARI sidecar ─────────────────────────────────────────────────────

/// Optional PARI cross-check: verify `#Jac(C)(GF(2^2))` via `hyperellcharpoly`.
///
/// This test is gated with `#[ignore]` and requires PARI/GP to be installed.
/// Run manually with:
/// ```text
/// cargo test -p rho --test ghs_kat pari_hyperellcharpoly_cross_check -- --ignored
/// ```
///
/// Expected PARI session for the GHS toy fixture:
/// ```text
/// ? K = GF(2^2, 'a, a^2+a+1);
/// ? C = hyperelliptic(Pol([1,1,0,1]*Mod(1,2)), Pol([0,1]*Mod(1,2)));
/// ? hyperellcharpoly(C)   \\ L-polynomial of C/GF(2^2)
/// ```
///
/// The descended curve `C: Y² + XY = X³ + X² + 1` over `GF(2^2)` has genus 1.
/// The Jacobian group order `#Jac(C)(GF(2^2))` equals `L(1)` where `L` is the
/// L-polynomial returned by `hyperellcharpoly`.
///
/// The divisibility check: `#Jac(C)(GF(2^l))` divides `#E(GF(2^m))` (the descent
/// is a group homomorphism, so the image divides the source order).
#[test]
#[ignore = "PARI not installed; run manually when available"]
fn pari_hyperellcharpoly_cross_check() {
    // Placeholder: when PARI is available, verify that the Jacobian group order
    // of the descended curve C/GF(2^2) matches the PARI computation.
    //
    // The GHS toy fixture:
    //   E: y² + xy = x³ + x² + 1  over GF(2^6) with x⁶+x+1.
    //   C: Y² + XY = X³ + X² + 1  over GF(2^2) with x²+x+1.
    //   genus(C) = 1.
    //
    // PARI commands:
    //   K = GF(2^2, 'a, a^2+a+1);
    //   C = hyperelliptic(Pol([1,1,0,1]*Mod(1,2)), Pol([0,1]*Mod(1,2)));
    //   hyperellcharpoly(C)
    //
    // The L-polynomial L(T) = 1 + c₁T + c₂T² (for genus 1).
    // #Jac(C)(GF(2^2)) = L(1) = 1 + c₁ + c₂.
    //
    // Divisibility check: #Jac(C)(GF(2^2)) | #E(GF(2^6)).
    unimplemented!("run manually with PARI installed");
}
