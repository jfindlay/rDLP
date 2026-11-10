//! GHS (Gaudry–Hess–Smart) Weil-descent attack on binary elliptic curves.
//!
//! This module implements the algebraic machinery for the GHS construction, which
//! transfers an ECDLP on a binary elliptic curve `E/GF(2^m)` to a DLP on the
//! Jacobian of a hyperelliptic curve `C/GF(2^l)` over a subfield.
//!
//! # Structure
//!
//! - [`mod`] (this file) — `GhsError` enum, the toy binary-curve fixture
//!   (m=6, l=2), and `check_ghs_params` (the `l | m` precondition verifier).
//! - [`descent`] — the Artin–Schreier extension and Weil restriction of scalars
//!   (`ArtinSchreierData`, `WeilRestriction`, `weil_restrict_poly`).
//! - [`curve`] — the GHS hyperelliptic-curve extraction `C/GF(2^l)` from the
//!   descent algebra (`extract_ghs_curve`, `ghs_genus`).
//! - [`reduce`] — the GHS reduction `(E, g, h) → (C, D_g, D_h)` and the
//!   logarithm-preservation verifier (`ghs_descend`, `GhsDescentResult`).
//!
//! # GHS construction overview
//!
//! The GHS attack proceeds in four stages:
//! 1. **Descent algebra** (this module + `descent`): build the Artin–Schreier
//!    extension `y² + y = f(x)` of the function field of `E/GF(2^m)`, then apply
//!    the Weil restriction `Res_{GF(2^m)/GF(2^l)}` to lower the field from
//!    `GF(2^m)` to `GF(2^l)` and raise the dimension from 1 to `m/l`.
//! 2. **Curve extraction** (`curve`): extract the hyperelliptic curve `C/GF(2^l)`
//!    from the descent algebra via [`extract_ghs_curve`].  The genus is
//!    `g = (m/l − 1)/2` for odd `m/l` (imaginary model).
//! 3. **Transfer map** ([`transfer`]): carry a point on `E` to a divisor on `Jac(C)`.
//! 4. **Reduction** (`reduce`): combine steps 2–3 into the top-level reduction
//!    `(E, g, h) → (C, D_g, D_h)` with logarithm-preservation guarantee.
//!
//! # Toy fixture
//!
//! The canonical GHS fixture: m=6, l=2, m/l=3 (odd — imaginary model).
//! - Source field: GF(2^6) with irreducible `x⁶+x+1` (poly = 0x43).
//! - Subfield: GF(2^2) with irreducible `x²+x+1` (poly = 0x7).
//! - Binary curve `E/GF(2^6)`: `y²+xy = x³+x²+1` (a=1, b=1 over GF(2^6)).
//!
//! Odd `m/l` keeps the descent in the imaginary/ramified hyperelliptic model
//! that [`crate::hyperelliptic::HyperellipticCurve`] handles.
//!
//! # Scale limits
//!
//! The fixture is toy-scale (m=6, l=2). The construction does not depend on it: the subfield
//! is identified through `shared_gf2m::restrict` and checked by [`check_subfield_identification`]
//! for any parameters, and the field arithmetic is single-limb, so `m ≤ 63`. Only genus ≥ 2
//! (`m/l ≥ 5`) would make the transfer useful against a real ECDLP, and no such fixture is tested.

pub mod curve;
pub mod descent;
pub mod reduce;
pub mod transfer;

pub use curve::{extract_ghs_curve, ghs_genus};
pub use descent::{ArtinSchreierData, GhsParams, WeilRestriction, weil_restrict_poly};
pub use reduce::{GhsDescentResult, ghs_descend, verify_log_preservation};
pub use transfer::{transfer_point, verify_homomorphism};

use crypto_bigint::Uint;
use shared_gf2m::{BitLen, F2m, F2mNaive};

use crate::binary_curve::BinaryCurve;

// ─── error type ──────────────────────────────────────────────────────────────

/// Errors from the GHS descent construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GhsError {
    /// `l` does not divide `m` — `GF(2^l) ⊆ GF(2^m)` requires `l | m`.
    ///
    /// The Weil restriction `Res_{GF(2^m)/GF(2^l)}` is only defined when
    /// `GF(2^l)` is a subfield of `GF(2^m)`, which holds iff `l | m`.
    SubfieldDivisibility,
    /// The curve does not admit a GHS descent.
    ///
    /// For example, the extension degree `m/l` is even for the imaginary model
    /// (even `m/l` yields the real/split model, which
    /// [`crate::hyperelliptic::HyperellipticCurve`] does not handle).
    NonDescendable,
    /// The point does not satisfy the Weierstrass equation of `E`.
    ///
    /// The transfer map is a homomorphism only on `E(GF(2^m))`.
    NotOnCurve,
    /// A coordinate of the conorm image is not in the subfield `GF(2^l)`.
    ///
    /// The image of a point of `E` is fixed by the `l`-th Frobenius power only when
    /// `E` is defined over `GF(2^l)`; otherwise the image cannot be read as a point
    /// of `C(GF(2^l))`.
    NotInSubfield,
    /// `poly_l` is not the minimal polynomial of `β = α^((2^m−1)/(2^l−1))`.
    ///
    /// `GF(2^l)` sits inside `GF(2^m)` as `GF(2)(β)`, and a `GF(2^l)` element
    /// `Σ cᵢ·xⁱ` (reduced modulo `poly_l`) names the element `Σ cᵢ·βⁱ` of `GF(2^m)`
    /// only when `poly_l(β) = 0` and `deg poly_l = l`.  Any other `poly_l` describes
    /// a different identification of the subfield, so no coordinate can be trusted.
    SubfieldPolyMismatch,
}

impl std::fmt::Display for GhsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GhsError::SubfieldDivisibility => {
                write!(f, "l does not divide m: GF(2^l) ⊆ GF(2^m) requires l | m")
            }
            GhsError::NonDescendable => {
                write!(f, "curve does not admit a GHS descent (e.g. even m/l for imaginary model)")
            }
            GhsError::NotOnCurve => {
                write!(f, "point does not satisfy the curve equation of E")
            }
            GhsError::NotInSubfield => {
                write!(f, "conorm image has a coordinate outside the subfield GF(2^l)")
            }
            GhsError::SubfieldPolyMismatch => {
                write!(f, "poly_l is not the minimal polynomial of the GF(2^l) generator β")
            }
        }
    }
}

impl std::error::Error for GhsError {}

// ─── precondition verifier ────────────────────────────────────────────────────

/// Check the GHS preconditions: `l | m` (subfield divisibility).
///
/// Returns `Ok(())` if `l` divides `m`, or `Err(GhsError::SubfieldDivisibility)`
/// otherwise.  This is the entry-point guard for all GHS descent operations.
///
/// # Arguments
///
/// - `m` — the extension degree of the source field `GF(2^m)`.
/// - `l` — the subfield degree `GF(2^l)`, with `l | m` required.
pub fn check_ghs_params(m: usize, l: usize) -> Result<(), GhsError> {
    if l == 0 || m % l != 0 {
        return Err(GhsError::SubfieldDivisibility);
    }
    Ok(())
}

// `shared_gf2m::restrict` is defined for m ≤ 63 (the exponent (2^m − 1)/(2^l − 1) must fit
// in a u64).  Every `GF(2^m)` here is a single-limb `Uint<1>` polynomial of at most 64 bits,
// so m ≤ 63 holds by construction; this keeps the guarantee from silently lapsing if the
// limb width ever grows.
const _: () = assert!(Uint::<1>::BITS <= 64, "restrict needs m <= 63");

/// Check that `GF(2^l) ⊂ GF(2^m)` can be identified through `shared_gf2m::restrict`
/// using `poly_l` as the coordinate field.
///
/// The subfield is identified by the generator `β = α^((2^m−1)/(2^l−1))`, where `α` is
/// the root of `poly_m`; a `GF(2^l)` element with bit-vector `c` (modulo `poly_l`) stands
/// for `Σ cᵢ·βⁱ`.  That reading is a field isomorphism exactly when `poly_l` is the
/// minimal polynomial of `β`, i.e. `deg poly_l = l` and `poly_l(β) = 0` (`β` has degree
/// `l` over `GF(2)`, so a monic degree-`l` polynomial vanishing at `β` equals its minimal
/// polynomial).
///
/// # Arguments
///
/// - `l` — the subfield degree.
/// - `poly_m` — the irreducible polynomial of `GF(2^m)`; `m` is its degree.
/// - `poly_l` — the claimed irreducible polynomial of `GF(2^l)`.
///
/// # Errors
///
/// - [`GhsError::SubfieldDivisibility`] if `l = 0`, `m = 0` or `l ∤ m`.
/// - [`GhsError::SubfieldPolyMismatch`] if `poly_l` is not the minimal polynomial of `β`.
pub fn check_subfield_identification(
    l: usize,
    poly_m: &Uint<1>,
    poly_l: &Uint<1>,
) -> Result<(), GhsError> {
    let m = poly_m.bit_len().saturating_sub(1);
    if l == 0 || m == 0 || m % l != 0 {
        return Err(GhsError::SubfieldDivisibility);
    }
    if poly_l.bit_len() != l + 1 {
        return Err(GhsError::SubfieldPolyMismatch);
    }

    // β = α^((2^m − 1)/(2^l − 1)); both factors fit in u64 for m ≤ 63.
    let alpha = F2mNaive::<1>::from_u64(2, poly_m);
    let exponent = ((1u64 << m) - 1) / ((1u64 << l) - 1);
    let beta = alpha.pow(&Uint::<1>::from(exponent), poly_m);

    // Horner evaluation of poly_l at β.
    let mut acc = F2mNaive::<1>::zero();
    for i in (0..=l).rev() {
        acc = acc.mul(&beta, poly_m);
        if poly_l.bit_at(i) {
            acc = acc.add(&F2mNaive::<1>::one());
        }
    }
    if !acc.is_zero() {
        return Err(GhsError::SubfieldPolyMismatch);
    }
    Ok(())
}

// ─── toy fixture ─────────────────────────────────────────────────────────────

/// The GF(2^6) irreducible polynomial: `x⁶+x+1` = 0x43.
///
/// This is the canonical GHS fixture field. The polynomial `x⁶+x+1` is
/// irreducible over GF(2) (verified: no roots in GF(2), no degree-2 or
/// degree-3 factor).
pub const GHS_POLY6: u64 = 0x43; // x^6 + x + 1

/// The GF(2^2) irreducible polynomial: `x²+x+1` = 0x7.
///
/// The subfield for the GHS fixture. `x²+x+1` is the unique irreducible
/// polynomial of degree 2 over GF(2).
pub const GHS_POLY2: u64 = 0x7; // x^2 + x + 1

/// Return the toy GHS binary curve fixture: `y²+xy = x³+x²+1` over GF(2^6).
///
/// Parameters:
/// - `poly = 0x43` (GF(2^6) with `x⁶+x+1`)
/// - `a = 1`, `b = 1` (curve coefficients)
/// - `n = 1` (placeholder group order — not used in the descent algebra)
/// - `gx = 0`, `gy = 1` (placeholder base point)
///
/// The curve equation `y²+xy = x³+x²+1` over GF(2^6) is a non-supersingular
/// binary elliptic curve (b = 1 ≠ 0). The composite extension degree m=6 with
/// subfield l=2 (m/l=3, odd) makes this curve amenable to the GHS descent in
/// the imaginary hyperelliptic model.
///
/// # Placeholder fields
///
/// The group order `n` and base point `(gx, gy)` are placeholders — the descent
/// algebra does not require them. The transfer map ([`transfer_point`]) will supply
/// a concrete point. The curve parameters `a` and `b` are the load-bearing data.
pub fn ghs_toy_curve() -> BinaryCurve {
    BinaryCurve {
        poly: Uint::<1>::from(GHS_POLY6),
        a: Uint::<1>::ONE,
        b: Uint::<1>::ONE,
        // Placeholder: group order not needed for the descent algebra.
        n: Uint::<1>::ONE,
        // Placeholder base point: (0, 1) satisfies y²+xy = x³+x²+1 at x=0:
        //   LHS: 1²+0·1 = 1; RHS: 0+0+1 = 1 ✓
        gx: Uint::<1>::ZERO,
        gy: Uint::<1>::ONE,
    }
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_gf2m::F2mNaive;

    #[test]
    fn check_ghs_params_ok() {
        // l=2 divides m=6: valid GHS parameters.
        assert_eq!(check_ghs_params(6, 2), Ok(()));
        // l=3 divides m=6: also valid.
        assert_eq!(check_ghs_params(6, 3), Ok(()));
        // l=1 divides m=6: trivial subfield.
        assert_eq!(check_ghs_params(6, 1), Ok(()));
        // l=6 divides m=6: trivial (same field).
        assert_eq!(check_ghs_params(6, 6), Ok(()));
    }

    #[test]
    fn check_ghs_params_err_subdivisibility() {
        // l=4 does not divide m=6: invalid.
        assert_eq!(check_ghs_params(6, 4), Err(GhsError::SubfieldDivisibility));
        // l=5 does not divide m=6: invalid.
        assert_eq!(check_ghs_params(6, 5), Err(GhsError::SubfieldDivisibility));
        // l=0: invalid (zero subfield degree).
        assert_eq!(check_ghs_params(6, 0), Err(GhsError::SubfieldDivisibility));
    }

    #[test]
    fn toy_curve_base_point_on_curve() {
        // The placeholder base point (0, 1) must satisfy y²+xy = x³+x²+1.
        let c = ghs_toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        assert!(c.is_on_curve(&g), "toy GHS base point must be on curve");
    }

    #[test]
    fn ghs_error_display() {
        // Smoke-test Display impls (no panic).
        let _ = format!("{}", GhsError::SubfieldDivisibility);
        let _ = format!("{}", GhsError::NonDescendable);
        let _ = format!("{}", GhsError::NotOnCurve);
        let _ = format!("{}", GhsError::NotInSubfield);
        let _ = format!("{}", GhsError::SubfieldPolyMismatch);
    }

    #[test]
    fn subfield_identification_accepts_minimal_polynomials() {
        let p6 = Uint::<1>::from(GHS_POLY6);
        // β = α^21 has minimal polynomial x²+x+1; β = α^9 has x³+x²+1.
        assert_eq!(
            check_subfield_identification(2, &p6, &Uint::<1>::from(GHS_POLY2)),
            Ok(())
        );
        assert_eq!(check_subfield_identification(3, &p6, &Uint::<1>::from(0xdu64)), Ok(()));
    }

    #[test]
    fn subfield_identification_rejects_other_polynomials() {
        let p6 = Uint::<1>::from(GHS_POLY6);
        // x³+x+1 is irreducible but is the reciprocal of the minimal polynomial of
        // β = α^9, so it is the minimal polynomial of β⁻¹, not of β.
        for bad in [0xbu64, 0x5, 0x3, 0x1f] {
            assert_eq!(
                check_subfield_identification(3, &p6, &Uint::<1>::from(bad)),
                Err(GhsError::SubfieldPolyMismatch),
                "poly_l = {bad:#x} must be rejected for l = 3"
            );
        }
    }

    #[test]
    fn subfield_identification_rejects_bad_degrees() {
        let p6 = Uint::<1>::from(GHS_POLY6);
        let p2 = Uint::<1>::from(GHS_POLY2);
        assert_eq!(
            check_subfield_identification(4, &p6, &p2),
            Err(GhsError::SubfieldDivisibility)
        );
        assert_eq!(
            check_subfield_identification(0, &p6, &p2),
            Err(GhsError::SubfieldDivisibility)
        );
    }
}
