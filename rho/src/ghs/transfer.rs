//! GHS transfer map: `E(GF(2^m)) → Jac(C)(GF(2^l))`.
//!
//! This module implements the core of the GHS Weil-descent attack: the transfer
//! map that carries a point `P ∈ E(GF(2^m))` to a reduced Mumford divisor
//! `D_P ∈ Jac(C)(GF(2^l))`.
//!
//! # Mathematical background
//!
//! The GHS construction builds `C` as the Weil restriction of `E`, so there is
//! a natural group homomorphism `φ: E(GF(2^m)) → Jac(C)(GF(2^l))`.  This is
//! the **transfer map** (also called the GHS conorm/corestriction map).
//!
//! ## The conorm map
//!
//! The transfer map is the **conorm** (corestriction) from `E(GF(2^m))` to
//! `Jac(C)(GF(2^l))`.  For a point `P ∈ E(GF(2^m))`, the conorm is:
//!
//! ```text
//! φ(P) = P + φ_l(P) + φ_l²(P) + … + φ_l^{n-1}(P)
//! ```
//!
//! where `φ_l : (x, y) ↦ (x^(2^l), y^(2^l))` is the relative Frobenius
//! (the `l`-th power of the absolute Frobenius), and `n = m/l`.  The sum is
//! taken in the group law of `E` (= `C` over `GF(2^m)`).
//!
//! The result `R = φ(P)` is a point on `C(GF(2^l))` (it is Galois-invariant
//! under the relative Frobenius, hence defined over `GF(2^l)`).  The Mumford
//! divisor is `[X + x_R, y_R]` (degree-1 divisor for genus 1).
//!
//! ## Why this is a group homomorphism
//!
//! The relative Frobenius `φ_l` is a group endomorphism of `E` (it is a ring
//! homomorphism on the function field).  The conorm is a sum of group
//! endomorphisms, hence a group homomorphism:
//!
//! ```text
//! φ(P + Q) = (P+Q) + φ_l(P+Q) + … = (P + φ_l(P) + …) + (Q + φ_l(Q) + …) = φ(P) + φ(Q)
//! ```
//!
//! ## Homomorphism property
//!
//! The transfer is a group homomorphism:
//! - `φ(P + Q) = φ(P) + φ(Q)` (Cantor compose on `Jac(C)`).
//! - `φ(∞) = [1, 0]` (the zero divisor — identity maps to identity).
//!
//! This is the property that makes the GHS descent preserve discrete logarithms:
//! if `h = k·g` on `E`, then `D_h = k·D_g` on `Jac(C)`.
//!
//! ## Toy fixture
//!
//! For the toy fixture (m=6, l=2, n=3, genus=1):
//! - `E: y²+xy = x³+x²+1` over `GF(2^6)`.
//! - `C: Y²+XY = X³+X²+1` over `GF(2^2)` (same equation, subfield).
//! - `φ_2(x, y) = (x^4, y^4)` (the 2nd Frobenius power).
//! - `φ(P) = P + φ_2(P) + φ_2²(P)` (sum of 3 conjugates).
//!
//! # Contracts
//!
//! **Consumes:** [`extract_ghs_curve`] from [`crate::ghs::curve`];
//! [`cantor`] from [`crate::hyperelliptic::cantor`];
//! [`BinaryAffinePoint`] from [`crate::binary_curve`];
//! subfield arithmetic from `shared_gf2m`.
//! **Provides:** [`transfer_point`], [`verify_homomorphism`].

use crypto_bigint::Uint;
use shared_gf2m::{F2m, F2mNaive, is_in_subfield};

use crate::binary_curve::{BinaryAffinePoint, BinaryCurve};
use crate::ghs::{GhsError, GhsParams};
use crate::hyperelliptic::{HyperellipticCurve, MumfordDivisor};
use crate::hyperelliptic::cantor;
use shared_gf2m::Poly;

// ─── transfer_point ───────────────────────────────────────────────────────────

/// Transfer a point `P ∈ E(GF(2^m))` to a reduced Mumford divisor
/// `D_P ∈ Jac(C)(GF(2^l))`.
///
/// This is the GHS conorm (corestriction) map.  It is a group homomorphism:
/// `D_{P+Q} = D_P + D_Q` (Cantor compose on `Jac(C)`).
///
/// # Algorithm: conorm map
///
/// The transfer is the **conorm** from `E(GF(2^m))` to `Jac(C)(GF(2^l))`:
///
/// ```text
/// φ(P) = P + φ_l(P) + φ_l²(P) + … + φ_l^{n-1}(P)
/// ```
///
/// where `φ_l(x, y) = (x^(2^l), y^(2^l))` is the relative Frobenius and
/// `n = m/l`.  The sum is in the group law of `E`.  The result is a point
/// on `C(GF(2^l))` (Galois-invariant under `φ_l`).
///
/// # Arguments
///
/// - `p` — the source point on `E(GF(2^m))`.
/// - `curve_c` — the descended hyperelliptic curve `C/GF(2^l)` (from
///   [`extract_ghs_curve`]).
/// - `params` — the GHS parameters (source field `GF(2^m)`, subfield `GF(2^l)`,
///   curve `E`).
///
/// # Returns
///
/// A reduced Mumford divisor `[u(X), v(X)]` on `Jac(C)(GF(2^l))`.
/// - For `P = ∞`: returns the zero divisor `[1, 0]`.
/// - For a finite `P = (x_P, y_P)`: returns a degree-1 divisor `[X + x_R, y_R]`
///   where `R = P + φ_l(P) + … + φ_l^{n-1}(P)` is the conorm image in `C(GF(2^l))`.
///
/// # Errors
///
/// - [`GhsError::NotOnCurve`] if `p` does not satisfy the Weierstrass equation of
///   `params.curve`.  The conorm is a homomorphism only on points of `E`; an
///   off-curve point would be fed to the group law, whose output is meaningless.
/// - [`GhsError::NotInSubfield`] if a coordinate of the conorm image `R` is not
///   fixed by the `l`-th Frobenius power (`c^(2^l) ≠ c`).  For a point on a curve
///   defined over `GF(2^l)` the image is always fixed, so this signals that
///   `params.curve` has coefficients outside `GF(2^l)` and does not descend.
///
/// The point at infinity is not an error: it is the identity of `E`, and a
/// homomorphism sends the identity to the identity, so `∞ ↦ [1, 0]`.
///
/// # Homomorphism property
///
/// For any `P, Q ∈ E(GF(2^m))`:
/// ```text
/// cantor::add(curve_c, transfer(P), transfer(Q), poly_l) == transfer(curve_e.add(P, Q))
/// ```
///
/// # Principle-4 annotation
///
/// The toy fixture (m=6, l=2, n=3, genus=1) gives degree-1 divisors.  The
/// algorithm is correct for any GHS descent; crypto-scale GHS uses genus ≥ 2.
pub fn transfer_point(
    p: &BinaryAffinePoint<F2mNaive<1>>,
    curve_c: &HyperellipticCurve<1>,
    params: &GhsParams,
) -> Result<MumfordDivisor<F2mNaive<1>, 1>, GhsError> {
    // The point at infinity maps to the zero divisor [1, 0].
    if p.is_infinity() {
        return Ok(curve_c.zero_divisor::<F2mNaive<1>>());
    }

    // The conorm is only defined on points of E.
    if !params.curve.is_on_curve(p) {
        return Err(GhsError::NotOnCurve);
    }

    let poly_m = params.poly_m();
    let n = params.extension_degree(); // m/l

    // Compute the conorm: R = P + φ_l(P) + φ_l²(P) + … + φ_l^{n-1}(P).
    //
    // φ_l(x, y) = (x^(2^l), y^(2^l)) — the relative Frobenius.
    // The sum is in the group law of E (= C over GF(2^m)).
    // The result R is a point on C(GF(2^l)) (Galois-invariant).
    let r = conorm_sum(p, &params.curve, n, params.l, poly_m);

    // Extract the coordinates of R (which are in GF(2^l)).
    let (x_r, y_r) = match &r {
        BinaryAffinePoint::Infinity => {
            // The conorm sum is ∞: return the zero divisor.
            return Ok(curve_c.zero_divisor::<F2mNaive<1>>());
        }
        BinaryAffinePoint::Finite { x, y } => (x.clone(), y.clone()),
    };

    // Re-express x_R and y_R as GF(2^l) elements, rejecting a coordinate that is
    // not fixed by the l-th Frobenius power instead of truncating it into a
    // wrong subfield element.
    let x_r_l = restrict_to_subfield(&x_r, params)?;
    let y_r_l = restrict_to_subfield(&y_r, params)?;

    // Construct the degree-1 Mumford divisor [u(X), v(X)] where:
    //   u(X) = X + x_R  (monic degree-1; in char 2, X − x_R = X + x_R)
    //   v(X) = y_R      (constant)
    let u = Poly::from_coeffs(vec![x_r_l, F2mNaive::<1>::one()]);
    let v = Poly::from_coeffs(vec![y_r_l]);

    Ok(MumfordDivisor::new(u, v))
}

/// Re-express a `GF(2^m)` element of the subfield `GF(2^l)` in the `GF(2^l)`
/// representation (reduction polynomial `params.poly_l`).
///
/// Returns [`GhsError::NotInSubfield`] unless `c^(2^l) = c` (Galois-invariant under
/// the relative Frobenius).  The bit pattern of `c` is then reduced modulo
/// `poly_l`; on the toy fixture this is a field isomorphism of `GF(4)`.
fn restrict_to_subfield(
    c: &F2mNaive<1>,
    params: &GhsParams,
) -> Result<F2mNaive<1>, GhsError> {
    if !is_in_subfield(c, params.l, params.poly_m()) {
        return Err(GhsError::NotInSubfield);
    }
    Ok(F2mNaive::<1>::from_uint(c.to_uint(), &params.poly_l))
}

// ─── conorm_sum ───────────────────────────────────────────────────────────────

/// Compute the conorm sum `P + φ_l(P) + φ_l²(P) + … + φ_l^{n-1}(P)`.
///
/// The relative Frobenius `φ_l(x, y) = (x^(2^l), y^(2^l))` maps a point on `E`
/// to another point on `E` (since `E` is defined over `GF(2^l)`).  The conorm
/// sums `n = m/l` conjugates using the group law of `E`.
///
/// The result is a point on `C(GF(2^l))` (Galois-invariant under `φ_l`).
///
/// # Arguments
///
/// - `p` — the source point.
/// - `curve_e` — the binary elliptic curve `E`.
/// - `n` — the extension degree `m/l`.
/// - `l` — the subfield degree.
/// - `poly_m` — the irreducible polynomial for `GF(2^m)`.
fn conorm_sum(
    p: &BinaryAffinePoint<F2mNaive<1>>,
    curve_e: &BinaryCurve,
    n: usize,
    l: usize,
    poly_m: &Uint<1>,
) -> BinaryAffinePoint<F2mNaive<1>> {
    // Start with the identity.
    let mut acc = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
    // Current conjugate: starts at P, then φ_l(P), φ_l²(P), ...
    let mut cur = p.clone();

    for _ in 0..n {
        // Add the current conjugate to the accumulator.
        acc = curve_e.add(&acc, &cur);
        // Advance: apply φ_l (the l-th Frobenius power).
        cur = frobenius_l(&cur, l, poly_m);
    }

    acc
}

/// Apply the relative Frobenius `φ_l(x, y) = (x^(2^l), y^(2^l))` to a point.
///
/// The relative Frobenius is the `l`-th power of the absolute Frobenius
/// `(x, y) ↦ (x², y²)`.  It maps a point on `E/GF(2^l)` to another point on
/// `E/GF(2^l)` (since `E` is defined over `GF(2^l)`, the Frobenius fixes the
/// curve equation).
fn frobenius_l(
    p: &BinaryAffinePoint<F2mNaive<1>>,
    l: usize,
    poly_m: &Uint<1>,
) -> BinaryAffinePoint<F2mNaive<1>> {
    match p {
        BinaryAffinePoint::Infinity => BinaryAffinePoint::Infinity,
        BinaryAffinePoint::Finite { x, y } => {
            // Apply Frobenius l times: x ↦ x^(2^l), y ↦ y^(2^l).
            let mut xf = x.clone();
            let mut yf = y.clone();
            for _ in 0..l {
                xf = xf.frobenius(poly_m);
                yf = yf.frobenius(poly_m);
            }
            BinaryAffinePoint::Finite { x: xf, y: yf }
        }
    }
}

// ─── Homomorphism verification helper ─────────────────────────────────────────

/// Verify the homomorphism property: `D_{P+Q} = D_P + D_Q`.
///
/// Computes `transfer(P + Q)` and `cantor::add(transfer(P), transfer(Q))` and
/// checks they are equal.  This is the decisive correctness guard for the GHS
/// transfer map.
///
/// # Arguments
///
/// - `p`, `q` — two points on `E(GF(2^m))`.
/// - `curve_e` — the binary elliptic curve `E`.
/// - `curve_c` — the descended hyperelliptic curve `C/GF(2^l)`.
/// - `params` — the GHS parameters.
///
/// # Returns
///
/// `true` if the homomorphism holds for this pair of points.
pub fn verify_homomorphism(
    p: &BinaryAffinePoint<F2mNaive<1>>,
    q: &BinaryAffinePoint<F2mNaive<1>>,
    curve_e: &crate::binary_curve::BinaryCurve,
    curve_c: &HyperellipticCurve<1>,
    params: &GhsParams,
) -> bool {
    let poly_l = &params.poly_l;

    // Compute D_P and D_Q.
    let d_p = match transfer_point(p, curve_c, params) {
        Ok(d) => d,
        Err(_) => return false,
    };
    let d_q = match transfer_point(q, curve_c, params) {
        Ok(d) => d,
        Err(_) => return false,
    };

    // Compute D_P + D_Q via Cantor compose.
    let d_sum_cantor = cantor::add(curve_c, &d_p, &d_q, poly_l);

    // Compute P + Q on E, then transfer.
    let p_plus_q = curve_e.add(p, q);
    let d_p_plus_q = match transfer_point(&p_plus_q, curve_c, params) {
        Ok(d) => d,
        Err(_) => return false,
    };

    d_sum_cantor == d_p_plus_q
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_gf2m::F2mNaive;

    use crate::ghs::{GHS_POLY2, extract_ghs_curve, ghs_toy_curve};

    fn poly2() -> Uint<1> {
        Uint::<1>::from(GHS_POLY2)
    }

    fn toy_params() -> GhsParams {
        GhsParams::new(6, 2, ghs_toy_curve(), poly2()).expect("toy GHS params must be valid")
    }

    #[test]
    fn transfer_infinity_is_zero_divisor() {
        let params = toy_params();
        let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
        let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
        let d = transfer_point(&inf, &curve_c, &params).expect("transfer must succeed");
        assert!(d.is_zero(), "transfer(∞) must be the zero divisor [1, 0]");
    }

    #[test]
    fn transfer_base_point_is_valid() {
        let params = toy_params();
        let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let d = transfer_point(&g, &curve_c, &params).expect("transfer must succeed");
        assert!(
            curve_c.is_valid(&d),
            "transfer(G) must be a valid reduced divisor"
        );
    }

    #[test]
    fn transfer_base_point_is_nonzero() {
        let params = toy_params();
        let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let d = transfer_point(&g, &curve_c, &params).expect("transfer must succeed");
        assert!(!d.is_zero(), "transfer(G) must be a non-zero divisor");
    }

    #[test]
    fn transfer_homomorphism_base_point() {
        // D_{G+G} = D_G + D_G via Cantor compose.
        let params = toy_params();
        let curve_c = extract_ghs_curve(params.clone()).expect("extraction must succeed");
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        assert!(
            verify_homomorphism(&g, &g, &curve_e, &curve_c, &params),
            "homomorphism must hold for G + G"
        );
    }

}
