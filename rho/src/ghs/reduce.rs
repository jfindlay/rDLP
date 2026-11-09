//! GHS reduction: ECDLP on `E/GF(2^m)` → Jacobian-DLP on `Jac(C)/GF(2^l)`.
//!
//! This module implements the top-level GHS reduction step: given an ECDLP instance
//! `(E, g, h)` with `h = k·g` on `E`, produce a Jacobian-DLP instance `(C, D_g, D_h)`
//! such that `D_h = k·D_g` on `Jac(C)`.
//!
//! # The reduction
//!
//! The GHS reduction proceeds in three steps:
//! 1. **Verify preconditions**: `l | m` (subfield divisibility), odd `m/l` (imaginary model).
//! 2. **Extract the descended curve**: `C/GF(2^l)` from the GHS parameters.
//! 3. **Transfer the points**: `D_g = transfer(g, C, params)`, `D_h = transfer(h, C, params)`.
//!
//! The correctness property (logarithm preservation): if `h = k·g` on `E`, then
//! `D_h = k·D_g` on `Jac(C)`.  This follows from the transfer map being a group
//! homomorphism ([`transfer_point`] from [`crate::ghs::transfer`]).
//!
//! # Scope boundary
//!
//! This module produces `(C, D_g, D_h)` and verifies `D_h = k·D_g` for known `k`.
//! It does **not** solve for `k` via index calculus — that is [`crate::index_calculus`].
//!
//! # Dependencies
//!
//! **Consumes:**
//! - [`transfer_point`] from [`crate::ghs::transfer`].
//! - [`extract_ghs_curve`] from [`crate::ghs::curve`].
//! - [`cantor::scalar_mul`] from [`crate::hyperelliptic::cantor`].
//! - [`BinaryAffinePoint`] from [`crate::binary_curve`].
//!
//! **Provides:** [`GhsDescentResult`], [`ghs_descend`].

use shared_gf2m::F2mNaive;

use crate::binary_curve::BinaryAffinePoint;
use crate::ghs::{GhsError, GhsParams, extract_ghs_curve, transfer_point};
use crate::hyperelliptic::{HyperellipticCurve, MumfordDivisor};
use crate::hyperelliptic::cantor;

// ─── GhsDescentResult ─────────────────────────────────────────────────────────

/// The result of the GHS reduction: a Jacobian-DLP instance `(C, D_g, D_h)`.
///
/// Produced by [`ghs_descend`].  The logarithm-preservation invariant holds:
/// if `h = k·g` on `E`, then `d_h = k·d_g` on `Jac(C)`.
///
/// # Invariant
///
/// `d_g` and `d_h` are valid reduced Mumford divisors on `curve_c`.
/// The discrete logarithm is preserved: `log_{d_g} d_h = log_g h` on their
/// respective groups.
#[derive(Clone, Debug)]
pub struct GhsDescentResult {
    /// The descended hyperelliptic curve `C/GF(2^l)`.
    pub curve_c: HyperellipticCurve<1>,
    /// The image of the base point `g` under the GHS transfer map.
    ///
    /// `d_g = transfer(g, C, params)` — a reduced Mumford divisor on `Jac(C)`.
    pub d_g: MumfordDivisor<F2mNaive<1>, 1>,
    /// The image of the target point `h` under the GHS transfer map.
    ///
    /// `d_h = transfer(h, C, params)` — a reduced Mumford divisor on `Jac(C)`.
    pub d_h: MumfordDivisor<F2mNaive<1>, 1>,
}

// ─── ghs_descend ──────────────────────────────────────────────────────────────

/// Reduce an ECDLP instance `(E, g, h)` to a Jacobian-DLP instance `(C, D_g, D_h)`.
///
/// Given the GHS parameters (source field `GF(2^m)`, subfield `GF(2^l)`, curve `E`),
/// and an ECDLP instance `(g, h)` with `h = k·g` on `E`, produces a Jacobian-DLP
/// instance `(C, D_g, D_h)` such that `D_h = k·D_g` on `Jac(C)`.
///
/// # Algorithm
///
/// 1. Extract the descended curve `C/GF(2^l)` from the GHS parameters.
/// 2. Transfer `g → D_g = transfer_point(g, C, params)`.
/// 3. Transfer `h → D_h = transfer_point(h, C, params)`.
/// 4. Return `(C, D_g, D_h)`.
///
/// # Logarithm-preservation invariant
///
/// If `h = k·g` on `E`, then `D_h = k·D_g` on `Jac(C)`.  This follows from the
    /// transfer map being a group homomorphism ([`transfer_point`]):
/// ```text
/// D_h = transfer(h) = transfer(k·g) = k·transfer(g) = k·D_g
/// ```
///
/// # Arguments
///
/// - `params` — the GHS parameters (source field `GF(2^m)`, subfield `GF(2^l)`, curve `E`).
/// - `g` — the base point on `E(GF(2^m))`.
/// - `h` — the target point on `E(GF(2^m))` with `h = k·g` for some unknown `k`.
///
/// # Returns
///
/// A [`GhsDescentResult`] containing the descended curve `C` and the transferred
/// divisors `D_g`, `D_h`.
///
/// # Errors
///
/// - [`GhsError::NonDescendable`] if the curve does not admit a GHS descent (e.g.,
///   even `m/l` for the imaginary model).
/// - [`GhsError::NotOnCurve`] if `g` or `h` is not on `E`.
/// - [`GhsError::NotInSubfield`] if the conorm image of `g` or `h` has a coordinate
///   outside `GF(2^l)` (see [`transfer_point`]).
///
/// # Principle-4 annotation
///
/// The toy fixture (m=6, l=2, n=3, genus=1) gives degree-1 divisors.  The reduction
/// is crypto-scale-correct; only the parameters are toy for auditability.
pub fn ghs_descend(
    params: &GhsParams,
    g: &BinaryAffinePoint<F2mNaive<1>>,
    h: &BinaryAffinePoint<F2mNaive<1>>,
) -> Result<GhsDescentResult, GhsError> {
    // Step 1: extract the descended curve C/GF(2^l).
    let curve_c = extract_ghs_curve(params.clone())?;

    // Step 2: transfer g → D_g.
    let d_g = transfer_point(g, &curve_c, params)?;

    // Step 3: transfer h → D_h.
    let d_h = transfer_point(h, &curve_c, params)?;

    Ok(GhsDescentResult { curve_c, d_g, d_h })
}

// ─── logarithm-preservation verifier ─────────────────────────────────────────

/// Verify the logarithm-preservation invariant: `scalar_mul(D_g, k) == D_h`.
///
/// For a known scalar `k`, checks that `k·D_g = D_h` on `Jac(C)` via Cantor
/// scalar multiplication.  This is the decisive correctness signal for the GHS
/// reduction: it confirms that the transfer map preserves the discrete logarithm.
///
/// # Arguments
///
/// - `result` — the [`GhsDescentResult`] from [`ghs_descend`].
/// - `k` — the known scalar (the discrete logarithm `h = k·g` on `E`).
///
/// # Returns
///
/// `true` if `k·D_g = D_h` on `Jac(C)`.
pub fn verify_log_preservation(result: &GhsDescentResult, k: u64) -> bool {
    let poly_l = &result.curve_c.poly;
    let k_d_g = cantor::scalar_mul(&result.curve_c, &result.d_g, k, poly_l);
    k_d_g == result.d_h
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_gf2m::F2mNaive;

    use crate::binary_curve::BinaryCurve;
    use crate::ghs::{GHS_POLY2, ghs_toy_curve};

    fn poly2() -> Uint<1> {
        Uint::<1>::from(GHS_POLY2)
    }

    fn toy_params() -> GhsParams {
        GhsParams::new(6, 2, ghs_toy_curve(), poly2()).expect("toy GHS params must be valid")
    }

    #[test]
    fn ghs_descend_succeeds_for_toy_fixture() {
        let params = toy_params();
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let h = curve_e.generator::<F2mNaive<1>>(); // h = g (k=1)
        let result = ghs_descend(&params, &g, &h);
        assert!(result.is_ok(), "ghs_descend must succeed for toy fixture");
    }

    #[test]
    fn ghs_descend_d_g_is_nonzero() {
        let params = toy_params();
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let h = curve_e.generator::<F2mNaive<1>>();
        let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
        assert!(!result.d_g.is_zero(), "D_g must be a non-zero divisor");
    }

    #[test]
    fn ghs_descend_d_g_is_valid() {
        let params = toy_params();
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let h = curve_e.generator::<F2mNaive<1>>();
        let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
        assert!(
            result.curve_c.is_valid(&result.d_g),
            "D_g must be a valid reduced divisor"
        );
    }

    #[test]
    fn ghs_descend_infinity_gives_zero_d_h() {
        let params = toy_params();
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let inf = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
        let result = ghs_descend(&params, &g, &inf).expect("descent must succeed");
        assert!(result.d_h.is_zero(), "D_h for h=∞ must be the zero divisor");
    }

    #[test]
    fn ghs_descend_log_preservation_k1() {
        // h = 1·g → D_h = 1·D_g.
        let params = toy_params();
        let curve_e = ghs_toy_curve();
        let g = curve_e.generator::<F2mNaive<1>>();
        let h = g.clone(); // h = g (k=1)
        let result = ghs_descend(&params, &g, &h).expect("descent must succeed");
        assert!(
            verify_log_preservation(&result, 1),
            "log preservation must hold for k=1"
        );
    }

    #[test]
    fn ghs_descend_even_n_fails() {
        // Even m/l (n=2) must fail with NonDescendable.
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
        let g = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
        let h = BinaryAffinePoint::<F2mNaive<1>>::Infinity;
        let result = ghs_descend(&params, &g, &h);
        assert!(
            matches!(result, Err(GhsError::NonDescendable)),
            "even m/l must fail with NonDescendable"
        );
    }
}
