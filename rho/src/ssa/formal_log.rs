//! The p-adic logarithm `log(1 + t)` at the formal parameter `t = −X/Y`, for the SSA reduction.
//!
//! Given a kernel-of-reduction point in standard projective coordinates `(X, Y, Z)` over Z/p^k,
//! this module evaluates the p-adic logarithm from `shared-padic` at `1 + t`, where
//! `t = −X/Y` is the formal group parameter in standard projective coordinates.
//!
//! This is the logarithm of the *multiplicative* formal group evaluated at `t`; it is not the
//! elliptic curve's own formal logarithm `ψ`. The two agree with `t` to leading order
//! (`log(1 + t) = t − t²/2 + …` and `ψ(t) = t + O(t²)`), and the quotient of leading units that
//! the reduction takes uses only that leading term. That agreement is inferred from the series,
//! not proved here; the reduction's verification step makes the answer independent of it.
//!
//! # Mathematical background
//!
//! For a point in the kernel of reduction on a short Weierstrass curve, the formal group
//! parameter is `t = −x/y` where `(x, y)` are affine coordinates. In standard projective
//! coordinates `(X:Y:Z)` with affine `(x, y) = (X/Z, Y/Z)`, this becomes `t = −X/Y`
//! (the Z factors cancel).
//!
//! For a kernel-of-reduction point, `v_p(t) ≥ 1`, satisfying the convergence condition of the
//! p-adic logarithm.
//!
//! # Projective coordinate representation
//!
//! After the multiply-by-p step in `reduce.rs`, the result is in standard projective coordinates
//! `(X, Y, Z)` where `v_p(Z) ≥ 1` (the point reduces to infinity mod p). The formal group
//! parameter `t = −X/Y` has `v_p(t) ≥ 1` because `v_p(X) ≥ v_p(Y) + 1` for kernel points.
//!
//! # Scale limits
//!
//! Toy precision only. The series is evaluated at the precision of the input `Zp` elements, and
//! the multiply-by-p step that produces those inputs is O(p), so the whole reduction is
//! practical only for toy primes.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use shared_padic::{PadicLogError, Zp, ZpError, padic_log};

use crate::ssa::SsaError;

// ─── log(1 + t) at t = −X/Y ───────────────────────────────────────────────────

/// Compute `log(1 + t)` at `t = −X/Y` from standard projective coordinates.
///
/// Given a kernel-of-reduction point `(X, Y, Z)` in standard projective coordinates over Z/p^k
/// (where affine `(x, y) = (X/Z, Y/Z)`), computes the p-adic logarithm of `1 + t` where
/// `t = −X/Y`.
///
/// The formal group parameter `t = −X/Y` has `v_p(t) ≥ 1` for kernel-of-reduction points,
/// satisfying the convergence condition of the p-adic logarithm.
///
/// # Arguments
///
/// - `x_proj` — the X projective coordinate (in Z/p^k).
/// - `y_proj` — the Y projective coordinate (in Z/p^k).
/// - `k` — the p-adic precision.
///
/// # Returns
///
/// `log(1 + t)` as a `Zp` element at precision `k − v_p(t)`.
///
/// # Errors
///
/// [`SsaError::InvariantFailure`], with a message naming the breach, if:
///
/// - `X` or `Y` is zero at working precision (its valuation is at least the precision; this is
///   never read as valuation 0);
/// - `v_p(t) < 1` (the point is not in the kernel of reduction, so the series does not converge);
/// - `t` has no precision left after its valuation is removed;
/// - the p-adic layer rejects `t`, `1 + t` or the logarithm, or reports the convergence
///   condition violated after the kernel check passed.
///
/// On a lift where `p·P̃` has valuation exactly 1 these cannot occur for a point of an anomalous
/// curve. On the stored `a = 0` lift `p·P̃` may lie in `E₂`, and the zero-residue and
/// precision breaches can then occur on valid input.
///
/// # Scale limits
///
/// SCALE: toy precision only. The series is truncated at the input precision `k`.
pub fn elliptic_log_proj(x_proj: &Zp, y_proj: &Zp) -> Result<Zp, SsaError> {
    let p = x_proj.prime().clone();
    let k = x_proj.precision().min(y_proj.precision());

    let x_res = x_proj.residue().clone();
    let y_res = y_proj.residue().clone();

    // A zero residue means the valuation is at least the precision; `padic_valuation_bigint`
    // returns 0 for zero, which would silently misread it, so it is a breach.
    if x_res.is_zero() {
        return Err(SsaError::InvariantFailure("projective X vanishes at working precision"));
    }
    if y_res.is_zero() {
        return Err(SsaError::InvariantFailure("projective Y vanishes at working precision"));
    }

    // Compute t = −X/Y in Z_p.
    //
    // X and Y are integers mod p^k. For a kernel-of-reduction point (after multiply-by-p),
    // v_p(X) ≥ v_p(Y) + 1, so t = −X/Y has v_p(t) ≥ 1.
    //
    // Factor out p-parts: X = p^vX * X', Y = p^vY * Y' (units X', Y').
    // t = −p^(vX−vY) * X' / Y'.
    //
    // Precision of t: k_t = k − (vX − vY).
    let vx = padic_valuation_bigint(&x_res, &p);
    let vy = padic_valuation_bigint(&y_res, &p);

    // vt = vX − vY. For kernel-of-reduction points, vt ≥ 1.
    let vt = vx as i64 - vy as i64;
    if vt < 1 {
        // The point is not in the kernel of reduction. #E = p puts every p·P̃ in the kernel of
        // reduction, so on a proven anomalous curve this is an invariant failure.
        return Err(SsaError::InvariantFailure(
            "formal parameter has valuation below 1: point outside the kernel of reduction",
        ));
    }
    let vt = vt as u32;

    // Precision of t: k_t = k − vt.
    let k_t = k.saturating_sub(vt);
    if k_t < 1 {
        return Err(SsaError::InvariantFailure(
            "formal parameter has no precision left after removing its valuation",
        ));
    }
    let pk_t = pow_bigint(&p, k_t);

    // Factor out p-parts.
    let p_pow_vx = pow_bigint(&p, vx);
    let p_pow_vy = pow_bigint(&p, vy);
    let x_unit = &x_res / &p_pow_vx; // X' (unit part of X)
    let y_unit = &y_res / &p_pow_vy; // Y' (unit part of Y)

    // t = −p^vt * X' / Y' mod p^k_t.
    let p_pow_vt = pow_bigint(&p, vt);
    let y_unit_inv = mod_inverse_bigint(&y_unit, &pk_t);
    let t_residue = (-(&p_pow_vt * &x_unit * &y_unit_inv)).mod_floor(&pk_t);

    // Construct t as a Zp element at precision k_t.
    let zp_breach = |why: &'static str| move |_: ZpError| SsaError::InvariantFailure(why);
    let t = Zp::new(&t_residue, &p, k_t)
        .map_err(zp_breach("p-adic layer rejected the formal parameter t"))?;

    // Construct z = 1 + t for the padic_log call.
    // v_p(z − 1) = v_p(t) = vt ≥ 1. ✓ (convergence condition satisfied)
    let one = Zp::new(&BigInt::one(), &p, k_t)
        .map_err(zp_breach("p-adic layer rejected the unit 1 at the precision of t"))?;
    let z = one.add(&t).map_err(zp_breach("p-adic layer rejected 1 + t"))?;

    // Apply the p-adic logarithm: log(1 + t) = log(z).
    let log_z = padic_log(&z).map_err(|e| match e {
        // vt ≥ 1 was checked above, so the convergence condition holds.
        PadicLogError::ConvergenceViolation { .. } => SsaError::InvariantFailure(
            "p-adic logarithm convergence violated after the kernel-of-reduction check",
        ),
        PadicLogError::Zp(_) => {
            SsaError::InvariantFailure("p-adic layer failed inside the logarithm of 1 + t")
        }
    })?;

    Ok(log_z)
}

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Compute the p-adic valuation of a `BigInt` `x` (number of times p divides x).
///
/// Returns 0 for x = 0 (callers handle the zero case separately).
pub(crate) fn padic_valuation_bigint(x: &BigInt, p: &BigInt) -> u32 {
    if x.is_zero() {
        return 0;
    }
    let mut x = x.abs();
    let mut e: u32 = 0;
    while (&x % p).is_zero() {
        x /= p;
        e += 1;
    }
    e
}

/// Compute `p^k` as a `BigInt`.
pub(crate) fn pow_bigint(p: &BigInt, k: u32) -> BigInt {
    num_traits::pow(p.clone(), k as usize)
}

/// Compute the modular inverse of `a` mod `m` using the extended Euclidean algorithm.
///
/// Precondition: `gcd(a, m) = 1`. Panics if not satisfied.
pub(crate) fn mod_inverse_bigint(a: &BigInt, m: &BigInt) -> BigInt {
    let (mut old_r, mut r) = (a.clone(), m.clone());
    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());

    while !r.is_zero() {
        let q = &old_r / &r;
        let new_r = &old_r - &q * &r;
        let new_s = &old_s - &q * &s;
        old_r = r;
        r = new_r;
        old_s = s;
        s = new_s;
    }

    assert!(
        old_r.is_one(),
        "mod_inverse_bigint called on non-unit: gcd({a}, {m}) = {old_r}"
    );

    old_s.mod_floor(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zp(residue: i64, p: i64, k: u32) -> Zp {
        Zp::new(&BigInt::from(residue), &BigInt::from(p), k).expect("valid Zp")
    }

    /// A point outside the kernel of reduction (X and Y both units) is an invariant failure, not
    /// a p-adic parameter error.
    #[test]
    fn elliptic_log_proj_rejects_point_outside_kernel() {
        let result = elliptic_log_proj(&zp(3, 7, 8), &zp(2, 7, 8));
        assert!(matches!(result, Err(SsaError::InvariantFailure(_))), "got {result:?}");
    }

    /// A kernel point (v(X) = v(Y) + 1) yields a logarithm of valuation at least 1.
    #[test]
    fn elliptic_log_proj_accepts_kernel_point() {
        let log = elliptic_log_proj(&zp(7, 7, 8), &zp(1, 7, 8)).expect("kernel point");
        assert_ne!(*log.residue(), BigInt::zero());
    }

    /// A zero residue is a breach, never read as valuation 0.
    #[test]
    fn elliptic_log_proj_rejects_zero_residues() {
        for (x, y) in [(0, 1), (7, 0)] {
            let result = elliptic_log_proj(&zp(x, 7, 8), &zp(y, 7, 8));
            assert!(matches!(result, Err(SsaError::InvariantFailure(_))), "got {result:?}");
        }
    }

    /// A valuation of t that consumes the working precision (the lower of the two input
    /// precisions) is a breach, distinct from the kernel check.
    #[test]
    fn elliptic_log_proj_rejects_exhausted_precision() {
        // v(X) = 5 at precision 8, Y a unit at precision 2: v(t) = 5 exceeds k = 2.
        let result = elliptic_log_proj(&zp(16807, 7, 8), &zp(1, 7, 2));
        assert!(matches!(result, Err(SsaError::InvariantFailure(_))), "got {result:?}");
    }
}
