//! GHS hyperelliptic-curve extraction `C/GF(2^l)` from the descent algebra.
//!
//! This module implements the curve-extraction step of the GHS Weil-descent
//! construction: given the Artin–Schreier data and Weil-restriction parameters
//! from [`crate::ghs::descent`], it produces the descended hyperelliptic curve
//! `C/GF(2^l)` as a frozen [`HyperellipticCurve`].
//!
//! # Mathematical background
//!
//! The GHS construction descends the binary elliptic curve `E/GF(2^m)` to a
//! hyperelliptic curve `C/GF(2^l)` via the Weil restriction
//! `Res_{GF(2^m)/GF(2^l)}(E)`.  The descended curve is in the **imaginary
//! (ramified) model**:
//!
//! ```text
//! C: Y² + H(X)·Y = F(X)   over GF(2^l)
//! ```
//!
//! where:
//! - `deg F = 2g+1` (odd — imaginary model).
//! - `deg H ≤ g`.
//! - The genus `g = (n−1)/2` with `n = m/l` (the extension degree), valid for
//!   odd `n`.
//!
//! ## Curve extraction for the toy fixture
//!
//! For `E: y²+xy = x³+ax²+b` over `GF(2^m)` with `a, b ∈ GF(2^l)`, the
//! Artin–Schreier polynomial `f_AS(x) = x³+ax²+b` already has all coefficients
//! in `GF(2^l)`.  The descended curve is then:
//!
//! ```text
//! F(X) = X³ + aX² + b   (= f_AS viewed over GF(2^l))
//! H(X) = X              (from the xy term in the Weierstrass equation)
//! ```
//!
//! This is the correct imaginary-model curve: `deg F = 3 = 2·1+1` (genus 1 for
//! `n = m/l = 3`).  The construction is valid because the curve coefficients
//! `a, b ∈ GF(2^l)` ensure the Artin–Schreier polynomial descends without
//! further Weil-restriction work.
//!
//! ## Genus formula
//!
//! For the GHS construction with odd extension degree `n = m/l`:
//!
//! ```text
//! g = (n − 1) / 2
//! ```
//!
//! | n | g | Note |
//! |---|---|------|
//! | 3 | 1 | toy fixture (m=6, l=2); genus-1 = elliptic in hyperelliptic form |
//! | 5 | 2 | crypto-scale genus-2 |
//! | 7 | 3 | crypto-scale genus-3 |
//!
//! # Principle-4 annotation
//!
//! The toy fixture (m=6, l=2, n=3) gives genus 1 — a genus-1 hyperelliptic
//! curve (an elliptic curve in Weierstrass-hyperelliptic form).  This is NOT the
//! crypto-scale case (genus ≥ 2 is needed for the GHS attack to be effective).
//! The construction is mathematically correct; only the parameters are toy for
//! auditability.  Crypto-scale GHS uses n ≥ 5 (genus ≥ 2).
//!
//! # Contracts
//!
//! **Consumes:** [`GhsParams`], [`ArtinSchreierData`] from [`crate::ghs::descent`].
//! **Provides:** [`extract_ghs_curve`], [`ghs_genus`].

use crypto_bigint::Uint;
use shared_gf2m::{F2m, is_in_subfield, restrict};

use crate::ghs::{ArtinSchreierData, GhsError, GhsParams, check_subfield_identification};
use crate::hyperelliptic::HyperellipticCurve;

// ─── genus formula ────────────────────────────────────────────────────────────

/// Compute the genus of the GHS descended curve from the extension degree.
///
/// For the GHS construction with odd extension degree `n = m/l`, the descended
/// curve `C/GF(2^l)` has genus:
///
/// ```text
/// g = (n − 1) / 2
/// ```
///
/// This formula holds for the imaginary (ramified) hyperelliptic model, which
/// requires odd `n`.  Even `n` yields the real/split model (not handled by the
/// frozen [`crate::hyperelliptic::HyperellipticCurve`]).
///
/// # Arguments
///
/// - `m` — the extension degree of the source field `GF(2^m)`.
/// - `l` — the subfield degree `GF(2^l)`, with `l | m` and odd `m/l`.
///
/// # Panics
///
/// Panics if `l` is zero or does not divide `m`.
///
/// # Principle-4 annotation
///
/// For the toy fixture (m=6, l=2, n=3), this returns 1 — a genus-1 curve.
/// Crypto-scale GHS uses n ≥ 5 (genus ≥ 2).
pub fn ghs_genus(m: usize, l: usize) -> usize {
    assert!(l > 0 && m % l == 0, "ghs_genus: l must divide m (l={l}, m={m})");
    let n = m / l;
    (n - 1) / 2
}

// ─── curve extraction ─────────────────────────────────────────────────────────

/// Extract the GHS hyperelliptic curve `C/GF(2^l)` from the descent algebra.
///
/// Given the GHS parameters (source field `GF(2^m)`, subfield `GF(2^l)`, binary
/// elliptic curve `E/GF(2^m)`), constructs the descended hyperelliptic curve
/// `C: Y²+H(X)Y = F(X)` over `GF(2^l)` as a frozen [`HyperellipticCurve<1>`].
///
/// # Construction
///
/// 1. Derives the Artin–Schreier polynomial `f_AS(x) = x³+ax²+b` from `E`.
/// 2. Verifies that `f_AS` has all coefficients in `GF(2^l)` (required for the
///    imaginary-model descent without additional Weil-restriction work).
/// 3. Sets `F(X) = f_AS` over `GF(2^l)`, each coefficient replaced by its coordinates
///    in the subfield basis `{1, β, …, β^(l−1)}`.
/// 4. Sets `H(X) = X` (from the `xy` term in the Weierstrass equation).
/// 5. Returns `HyperellipticCurve` with `poly = poly_l`, `h_coeffs`, `f_coeffs`.
///
/// The resulting curve has genus `g = (m/l − 1)/2` (for odd `m/l`).
///
/// The `GF(2^l)` coefficients are the coordinates of the `GF(2^m)` coefficients in the
/// subfield basis `{1, β, …, β^(l−1)}` (`shared_gf2m::restrict`), the same identification
/// the transfer map uses, so the curve and the transferred divisors share one field.
///
/// # Errors
///
/// - [`GhsError::SubfieldDivisibility`] or [`GhsError::SubfieldPolyMismatch`] if `l` does
///   not divide the degree of the source-field polynomial, or `poly_l` is not the minimal
///   polynomial of `β = α^((2^m−1)/(2^l−1))`.
/// - [`GhsError::NonDescendable`] if `m/l` is even (the imaginary model requires
///   odd `m/l`; even `m/l` yields the real/split model).
/// - [`GhsError::NonDescendable`] if the Artin–Schreier polynomial is not
///   well-formed (degree ≠ 3, not monic, or zero constant term).
/// - [`GhsError::NonDescendable`] if any coefficient of `f_AS` is not in
///   `GF(2^l)` (the descent requires subfield coefficients for this construction).
///
/// # Principle-4 annotation
///
/// For the toy fixture (m=6, l=2, n=3), the extracted curve has genus 1 — a
/// genus-1 hyperelliptic curve (an elliptic curve in hyperelliptic form).  This
/// is NOT the crypto-scale case (genus ≥ 2).  The construction is correct; only
/// the parameters are toy.
pub fn extract_ghs_curve(params: GhsParams) -> Result<HyperellipticCurve<1>, GhsError> {
    // The subfield identification must be sound before any coefficient is read.
    check_subfield_identification(params.l, params.poly_m(), &params.poly_l)?;

    let n = params.extension_degree(); // m/l

    // The imaginary model requires odd m/l.
    if n % 2 == 0 {
        return Err(GhsError::NonDescendable);
    }

    let poly_m = *params.poly_m();
    let poly_l = params.poly_l;

    // Build the Artin–Schreier data: f_AS(x) = x³ + ax² + b.
    let as_data = ArtinSchreierData::from_params(params);

    // Verify the Artin–Schreier polynomial is well-formed.
    if !as_data.is_well_formed() {
        return Err(GhsError::NonDescendable);
    }

    // Verify all coefficients of f_AS are in GF(2^l).
    // This is required for the imaginary-model descent: the curve coefficients
    // a and b must be in the subfield GF(2^l) for the Artin–Schreier polynomial
    // to descend directly (without additional Weil-restriction work).
    let l = as_data.params.l;
    let deg = as_data.f_poly.degree().unwrap_or(0);
    for i in 0..=deg {
        let coeff = as_data.f_poly.coeff(i);
        if !is_in_subfield(&coeff, l, &poly_m) {
            return Err(GhsError::NonDescendable);
        }
    }

    // Extract the f_coeffs over GF(2^l).
    // A GF(2^m) bit-vector is not a GF(2^l) bit-vector: the subfield is GF(2)(β) with
    // β = α^((2^m−1)/(2^l−1)), and a GF(2^l) element's bits are its coordinates in the
    // basis {1, β, …, β^(l−1)}.  `restrict` computes those coordinates, the same
    // identification the transfer map applies to the conorm image.
    //
    // f_AS(x) = x³ + ax² + b  →  F(X) = X³ + aX² + b over GF(2^l)
    // coeffs: [b, 0, a, 1]  (index = degree)
    let f_coeffs: Vec<Uint<1>> = (0..=deg)
        .map(|i| {
            restrict(&as_data.f_poly.coeff(i), l, &poly_m)
                .map(|c| c.to_uint())
                .ok_or(GhsError::NonDescendable)
        })
        .collect::<Result<_, _>>()?;

    // H(X) = X: from the xy term in the Weierstrass equation y²+xy = x³+ax²+b.
    // In the Artin–Schreier form, the xy term contributes H(X) = X to the
    // descended curve's h polynomial.
    // h_coeffs = [0, 1]  (H(X) = 0·X⁰ + 1·X¹ = X)
    let h_coeffs: Vec<Uint<1>> = vec![Uint::<1>::ZERO, Uint::<1>::ONE];

    Ok(HyperellipticCurve::new(poly_l, h_coeffs, f_coeffs))
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_gf2m::F2mNaive;

    use crate::ghs::{GHS_POLY2, ghs_toy_curve};

    fn toy_params() -> GhsParams {
        GhsParams::new(6, 2, ghs_toy_curve(), Uint::<1>::from(GHS_POLY2))
            .expect("toy GHS params must be valid")
    }

    #[test]
    fn ghs_genus_toy_fixture() {
        // m=6, l=2, n=3 → genus = (3-1)/2 = 1.
        assert_eq!(ghs_genus(6, 2), 1, "toy fixture genus must be 1");
    }

    #[test]
    fn ghs_genus_n5() {
        // m=10, l=2, n=5 → genus = (5-1)/2 = 2 (crypto-scale).
        assert_eq!(ghs_genus(10, 2), 2, "n=5 genus must be 2");
    }

    #[test]
    fn ghs_genus_n7() {
        // m=14, l=2, n=7 → genus = (7-1)/2 = 3.
        assert_eq!(ghs_genus(14, 2), 3, "n=7 genus must be 3");
    }

    #[test]
    fn extract_ghs_curve_ok() {
        let params = toy_params();
        let result = extract_ghs_curve(params);
        assert!(result.is_ok(), "extract_ghs_curve must succeed for toy fixture");
    }

    #[test]
    fn extracted_curve_genus_is_1() {
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        assert_eq!(curve.genus(), 1, "extracted curve genus must be 1 (n=3, g=(3-1)/2=1)");
    }

    #[test]
    fn extracted_curve_imaginary_model() {
        // Imaginary model: deg f = 2g+1.
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let g = curve.genus();
        let deg_f = curve.f_coeffs.len() - 1;
        assert_eq!(deg_f, 2 * g + 1, "imaginary model: deg f must equal 2g+1");
    }

    #[test]
    fn extracted_curve_h_degree_le_genus() {
        // deg H ≤ g.
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let g = curve.genus();
        // h_coeffs has length deg_h + 1 (or 0 for zero polynomial).
        // H(X) = X has degree 1, and g = 1, so deg H = 1 ≤ g = 1.
        let h_poly = curve.h::<F2mNaive<1>>();
        let deg_h = h_poly.degree().unwrap_or(0);
        assert!(deg_h <= g, "deg H = {deg_h} must be ≤ genus g = {g}");
    }

    #[test]
    fn extracted_curve_f_is_monic() {
        // F(X) = X³ + X² + 1 is monic.
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let f_poly = curve.f::<F2mNaive<1>>();
        let lc = f_poly.leading_coeff().expect("f must be non-zero");
        assert!(lc.is_one(), "F(X) must be monic (leading coefficient = 1)");
    }

    #[test]
    fn extracted_curve_poly_is_subfield() {
        // The curve's poly must be the GF(2^2) irreducible.
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        assert_eq!(
            curve.poly,
            Uint::<1>::from(GHS_POLY2),
            "extracted curve poly must be the GF(2^2) irreducible"
        );
    }

    #[test]
    fn extract_ghs_curve_even_n_fails() {
        // Even m/l (n=2) must fail with NonDescendable.
        // Use m=4, l=2, n=2 (even — real/split model).
        use crate::binary_curve::BinaryCurve;
        let poly4 = Uint::<1>::from(0x13u64); // x⁴+x+1
        let poly2 = Uint::<1>::from(GHS_POLY2);
        let curve = BinaryCurve {
            poly: poly4,
            a: Uint::<1>::ONE,
            b: Uint::<1>::ONE,
            n: Uint::<1>::ONE,
            gx: Uint::<1>::ZERO,
            gy: Uint::<1>::ONE,
        };
        let params = GhsParams::new(4, 2, curve, poly2).expect("4/2 is valid");
        let result = extract_ghs_curve(params);
        assert!(
            matches!(result, Err(GhsError::NonDescendable)),
            "even m/l must fail with NonDescendable"
        );
    }

    #[test]
    fn extracted_curve_zero_divisor_is_valid() {
        // The zero divisor [1, 0] must be valid on the extracted curve.
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let zero_div = curve.zero_divisor::<F2mNaive<1>>();
        assert!(
            curve.is_valid(&zero_div),
            "zero divisor [1, 0] must be valid on the extracted curve"
        );
    }

    #[test]
    fn extracted_curve_h_is_x() {
        // H(X) = X: h_coeffs = [0, 1].
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let poly_l = Uint::<1>::from(GHS_POLY2);
        let h_poly = curve.h::<F2mNaive<1>>();
        assert_eq!(h_poly.degree(), Some(1), "H(X) must have degree 1");
        assert!(
            h_poly.coeff(0).is_zero(),
            "H(X) constant term must be 0 (H = X)"
        );
        assert!(
            h_poly.coeff(1).is_one(),
            "H(X) coefficient of X must be 1 (H = X)"
        );
        // Verify the poly is GF(2^2).
        assert_eq!(curve.poly, poly_l, "curve poly must be GF(2^2) irreducible");
    }

    #[test]
    fn extracted_curve_f_coeffs_match_artin_schreier() {
        // F(X) = X³ + X² + 1 (= f_AS for a=1, b=1 over GF(2^2)).
        // f_coeffs = [1, 0, 1, 1]  (index = degree).
        let params = toy_params();
        let curve = extract_ghs_curve(params).expect("extraction must succeed");
        let f_poly = curve.f::<F2mNaive<1>>();
        assert_eq!(f_poly.degree(), Some(3), "F(X) must have degree 3");
        // Constant term: b = 1.
        assert!(f_poly.coeff(0).is_one(), "F(X) constant term must be 1 (= b)");
        // X coefficient: 0.
        assert!(f_poly.coeff(1).is_zero(), "F(X) coefficient of X must be 0");
        // X² coefficient: a = 1.
        assert!(f_poly.coeff(2).is_one(), "F(X) coefficient of X² must be 1 (= a)");
        // X³ coefficient: 1 (monic).
        assert!(f_poly.coeff(3).is_one(), "F(X) coefficient of X³ must be 1 (monic)");
    }
}
