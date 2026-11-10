//! Miller's algorithm for computing `f_{ℓ,P}(Q)`.
//!
//! Miller's algorithm accumulates the rational function `f_{ℓ,P}` on the curve
//! associated to the divisor `ℓ·(P) − ℓ·(∞)`, evaluated at a point `Q`.  The
//! result is an element of `F_{p^k}*` and is the core building block for both
//! the Weil pairing ([`crate::pairing::weil`]) and the Tate pairing ([`crate::pairing::tate`]).
//!
//! # Algorithm
//!
//! The double-and-add accumulation over the bits of `ℓ`:
//!
//! ```text
//! T = P,  f = 1
//! for bit in bits(ℓ) from second-highest to lowest:
//!     f = f² · line(T, T, Q) / vertical(2T, Q)   -- doubling step
//!     T = 2T
//!     if bit == 1:
//!         f = f · line(T, P, Q) / vertical(T+P, Q)  -- addition step
//!         T = T + P
//! return f
//! ```
//!
//! # Line and vertical functions
//!
//! `line(T1, T2, Q)` is the line through `T1` and `T2` evaluated at `Q`:
//! - **Tangent** (`T1 = T2`): slope `λ = (3·x_T² + a) / (2·y_T)`,
//!   value `y_Q − y_T − λ·(x_Q − x_T)`.
//! - **Secant** (`T1 ≠ T2`): slope `λ = (y_{T2} − y_{T1}) / (x_{T2} − x_{T1})`,
//!   value `y_Q − y_{T1} − λ·(x_Q − x_{T1})`.
//! - **Vertical** (`T1 = −T2`): returns the vertical function `x_Q − x_T`.
//!
//! `vertical(T, Q)` = `x_Q − x_T` (the vertical line through `T` evaluated at `Q`).
//!
//! # Infinity handling
//!
//! When the intermediate point `T+T` or `T+P` is the point at infinity, the
//! corresponding vertical function evaluates to 1 (the divisor is at infinity
//! and contributes nothing to the rational function).
//!
//! # All arithmetic in `F_{p^k}`
//!
//! Both `T` and `Q` are [`PairingPoint`]s over [`FpExt`].  The curve
//! coefficients `a` and `b` are lifted from `F_p` to `F_{p^k}` via
//! [`FpExt::from_base`].

use crypto_bigint::Uint;
use shared_field::Fp;

use crate::curve::Curve;
use crate::pairing::ecext::PairingPoint;
use crate::pairing::fpext::{FpExt, IrreducibleModulus};

// ── Public entry point ────────────────────────────────────────────────────────

/// Compute `f_{ℓ,P}(Q)` via Miller's algorithm.
///
/// Returns the value of the rational function `f_{ℓ,P}` associated to the
/// divisor `ℓ·(P) − ℓ·(∞)`, evaluated at the point `Q`.  The result is an
/// element of `F_{p^k}*`.
///
/// # Parameters
///
/// - `curve` — the short-Weierstrass curve `y² = x³ + ax + b`; `a` and `b`
///   are lifted from `F_p` to `F_{p^k}` internally.
/// - `modulus` — the irreducible polynomial defining `F_{p^k}`.
/// - `p_point` — the base point `P ∈ E[ℓ]`.
/// - `q_point` — the evaluation point `Q ∈ E[ℓ]`, linearly independent of `P`.
/// - `ell` — the torsion order; must be ≥ 2.
///
/// # Panics
///
/// Panics if `ell < 2`.
pub fn miller_loop<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    p_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    ell: u64,
) -> FpExt<F> {
    assert!(ell >= 2, "miller_loop: ell must be >= 2");

    let p = &curve.p;
    let k = extension_degree(p_point, q_point, modulus);

    // Lift curve coefficient a into F_{p^k}.
    let a_ext = FpExt::from_base(F::from_uint(curve.a, p), k, p);

    // T = P, f = 1.
    let mut t = p_point.clone();
    let mut f = FpExt::one(k, p);

    // Find the position of the highest set bit in ell.
    // The highest bit is always 1 and initialises T=P, f=1; we skip it.
    let bit_len = 64 - ell.leading_zeros(); // number of bits in ell
    // Iterate from the second-highest bit down to bit 0.
    for i in (0..bit_len - 1).rev() {
        // ── Doubling step ─────────────────────────────────────────────────────
        // f = f² · line(T, T, Q) / vertical(2T, Q)
        f = f.square(modulus, p);

        let line_val = line_tangent(&t, q_point, &a_ext, modulus, p);
        f = f.mul(&line_val, modulus, p);

        let t_doubled = t.double(&a_ext, modulus, p);
        let vert_val = vertical_at(&t_doubled, q_point, k, modulus, p);
        // Divide by vertical: multiply by its inverse.
        f = f.mul(&vert_val.inv(modulus, p), modulus, p);

        t = t_doubled;

        // ── Addition step (if current bit is 1) ───────────────────────────────
        if (ell >> i) & 1 == 1 {
            // f = f · line(T, P, Q) / vertical(T+P, Q)
            let line_val = line_secant(&t, p_point, q_point, &a_ext, modulus, p);
            f = f.mul(&line_val, modulus, p);

            let t_plus_p = t.add(p_point, &a_ext, modulus, p);
            let vert_val = vertical_at(&t_plus_p, q_point, k, modulus, p);
            f = f.mul(&vert_val.inv(modulus, p), modulus, p);

            t = t_plus_p;
        }
    }

    f
}

// ── Line and vertical helpers ─────────────────────────────────────────────────

/// Evaluate the tangent line at `T` (i.e., the line through `T` and `T`) at `Q`.
///
/// Uses the doubling formula slope `λ = (3·x_T² + a) / (2·y_T)`.
/// Returns `y_Q − y_T − λ·(x_Q − x_T)`.
///
/// If `T` is the point at infinity or `y_T = 0` (so `2T = ∞`), returns 1
/// (the vertical at infinity contributes nothing).
fn line_tangent<F: Fp<4>>(
    t: &PairingPoint<F>,
    q: &PairingPoint<F>,
    a: &FpExt<F>,
    modulus: &IrreducibleModulus<F>,
    p: &Uint<4>,
) -> FpExt<F> {
    let k = a.degree();
    match (t, q) {
        (PairingPoint::Infinity, _) | (_, PairingPoint::Infinity) => FpExt::one(k, p),
        (PairingPoint::Finite { x: xt, y: yt }, PairingPoint::Finite { x: xq, y: yq }) => {
            // If y_T = 0, the tangent is vertical: 2T = ∞.  Return 1.
            if yt.is_zero(p) {
                return FpExt::one(k, p);
            }
            // λ = (3·x_T² + a) / (2·y_T)
            let three = FpExt::from_base(F::from_u64(3, p), k, p);
            let two = FpExt::from_base(F::from_u64(2, p), k, p);
            let xt2 = xt.square(modulus, p);
            let num = three.mul(&xt2, modulus, p).add(a, p);
            let den = two.mul(yt, modulus, p);
            let lam = num.mul(&den.inv(modulus, p), modulus, p);

            // y_Q − y_T − λ·(x_Q − x_T)
            let xdiff = xq.sub(xt, p);
            let lam_xdiff = lam.mul(&xdiff, modulus, p);
            yq.sub(yt, p).sub(&lam_xdiff, p)
        }
    }
}

/// Evaluate the secant line through `T1` and `T2` at `Q`.
///
/// Uses the chord formula slope `λ = (y_{T2} − y_{T1}) / (x_{T2} − x_{T1})`.
/// Returns `y_Q − y_{T1} − λ·(x_Q − x_{T1})`.
///
/// Falls back to [`line_tangent`] when `T1 = T2`, and returns the vertical
/// function `x_Q − x_{T1}` when `T1 = −T2` (the line is vertical).
///
/// Returns 1 if either point is the point at infinity.
fn line_secant<F: Fp<4>>(
    t1: &PairingPoint<F>,
    t2: &PairingPoint<F>,
    q: &PairingPoint<F>,
    a: &FpExt<F>,
    modulus: &IrreducibleModulus<F>,
    p: &Uint<4>,
) -> FpExt<F> {
    let k = a.degree();
    match (t1, t2, q) {
        (PairingPoint::Infinity, _, _) | (_, PairingPoint::Infinity, _) => FpExt::one(k, p),
        (_, _, PairingPoint::Infinity) => FpExt::one(k, p),
        (
            PairingPoint::Finite { x: x1, y: y1 },
            PairingPoint::Finite { x: x2, y: y2 },
            PairingPoint::Finite { x: xq, y: yq },
        ) => {
            let xdiff = x2.sub(x1, p);
            if xdiff.is_zero(p) {
                let ysum = y1.add(y2, p);
                if ysum.is_zero(p) {
                    // T1 = −T2: the line is vertical; return x_Q − x_{T1}.
                    return xq.sub(x1, p);
                }
                // T1 = T2: use the tangent.
                return line_tangent(t1, q, a, modulus, p);
            }
            // λ = (y_{T2} − y_{T1}) / (x_{T2} − x_{T1})
            let ydiff = y2.sub(y1, p);
            let lam = ydiff.mul(&xdiff.inv(modulus, p), modulus, p);

            // y_Q − y_{T1} − λ·(x_Q − x_{T1})
            let xq_diff = xq.sub(x1, p);
            let lam_xdiff = lam.mul(&xq_diff, modulus, p);
            yq.sub(y1, p).sub(&lam_xdiff, p)
        }
    }
}

/// Evaluate the vertical line through `T` at `Q`: `x_Q − x_T`.
///
/// Returns 1 if `T` is the point at infinity (the divisor is at infinity and
/// contributes nothing to the rational function).
///
/// If `x_Q = x_T`, then `Q = ±T` lies on the support of the vertical's divisor and the Miller
/// value at `Q` is undefined (the denominator is 0). The function returns 1 there only so the
/// caller's inversion does not divide by zero; that 1 is not a pairing value. This happens when
/// the evaluation point `Q` is a multiple of the Miller base point `P`. For `Q = P` the first
/// tangent-line value `y_Q − y_T − λ·(x_Q − x_T)` is already 0, so `f_{ℓ,P}(Q)` comes out as 0 and
/// `mov_reduce` reports `MovError::DegeneratePairing`
/// (`mov_reduce_degenerate_r_is_degenerate_pairing`). Other `Q ∈ ⟨P⟩` are not exercised by any
/// test, and nothing here guarantees a 0 result for them.
fn vertical_at<F: Fp<4>>(
    t: &PairingPoint<F>,
    q: &PairingPoint<F>,
    k: usize,
    modulus: &IrreducibleModulus<F>,
    p: &Uint<4>,
) -> FpExt<F> {
    match (t, q) {
        (PairingPoint::Infinity, _) => FpExt::one(k, p),
        (_, PairingPoint::Infinity) => FpExt::one(k, p),
        (PairingPoint::Finite { x: xt, .. }, PairingPoint::Finite { x: xq, .. }) => {
            let diff = xq.sub(xt, p);
            if diff.is_zero(p) {
                FpExt::one(k, p)
            } else {
                // The vertical function is x_Q − x_T; the caller divides by it (multiplies by
                // the inverse), so it is returned uninverted.
                let _ = modulus; // modulus not needed for subtraction
                diff
            }
        }
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Determine the extension degree `k` from the points and modulus.
///
/// Uses the degree of the modulus polynomial as the canonical source of `k`.
fn extension_degree<F: Fp<4>>(
    _p: &PairingPoint<F>,
    _q: &PairingPoint<F>,
    modulus: &IrreducibleModulus<F>,
) -> usize {
    modulus.degree()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_field::FpNaive;

    use crate::pairing::test_curves::pairing_toy;

    /// `miller_loop(P, Q, ℓ)` returns a non-trivial element of `F_{p^k}*`.
    ///
    /// A trivial (= 1) result would indicate a degenerate pairing.
    #[test]
    fn miller_loop_non_trivial() {
        let (curve, modulus, ell, p_point, q_point) = pairing_toy();
        let p = curve.p;
        let result = miller_loop(&curve, &modulus, &p_point, &q_point, ell);
        assert!(
            !result.is_one(&p),
            "miller_loop(P, Q, ℓ) should be non-trivial (≠ 1)"
        );
    }

    /// `miller_loop(P, Q, ℓ)` and `miller_loop(Q, P, ℓ)` are both non-zero.
    #[test]
    fn miller_loop_both_directions_non_zero() {
        let (curve, modulus, ell, p_point, q_point) = pairing_toy();
        let p = curve.p;
        let f_pq = miller_loop::<FpNaive<4>>(&curve, &modulus, &p_point, &q_point, ell);
        let f_qp = miller_loop::<FpNaive<4>>(&curve, &modulus, &q_point, &p_point, ell);
        assert!(!f_pq.is_zero(&p), "f_{{ℓ,P}}(Q) should be non-zero");
        assert!(!f_qp.is_zero(&p), "f_{{ℓ,Q}}(P) should be non-zero");
    }
}
