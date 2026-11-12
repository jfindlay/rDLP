//! Norm computation and the signed-BigInt → Uint<4> bridge for NFS sieving.
//!
//! NFS sieving requires computing two values for each candidate pair ``(a, b)``:
//!
//! - **Rational norm** ``N_rat(a, b) = a − b·m``: the value of the rational-side polynomial
//!   ``g(x) = x − m`` evaluated at ``a/b``, cleared of denominators. This is the norm from
//!   ℚ to ℚ of the element ``a/b − m`` in the rational number field.
//!
//! - **Algebraic value** ``F(a, b) = b^d · f(a/b)`` (the "homogeneous norm value", named
//!   ``N_alg`` in this module): the homogeneous form of ``f`` evaluated at ``(a, b)``. Computed
//!   directly as ``Σ f.coeffs[i] · a^i · b^{d−i}`` to avoid rational arithmetic. With ``c`` the
//!   leading coefficient of ``f`` and ``α`` a root, ``F(a, b) = c · N_{K/ℚ}(a − b·α)``, so it is
//!   the field norm ``N`` of ``a − b·α`` only when ``f`` is monic; in general it is ``c`` times
//!   that norm. "Norm" alone names ``N``.
//!
//! Both are signed ``BigInt`` values. The smoothness predicate ``trial_smooth`` operates
//! on ``Uint<4>`` (unsigned 256-bit integers). The **norm bridge** converts a signed value to
//! its absolute value as ``Uint<4>``, with a range check: toy-scale norms fit 256 bits per the
//! C1 resolution in ROADMAP, but the bridge rejects out-of-range norms rather than silently
//! truncating.
//!
//! # Science↔engineering note (principle 4)
//!
//! At toy scale, these values are small (tens of bits). At cryptographic scale (RSA-768+), the
//! algebraic value ``F(a, b) = b^d · f(a/b)`` can be hundreds of bits. The ``Uint<4>`` (256-bit)
//! bridge is sized for toy-scale values; a production implementation would use a wider type or a
//! different smoothness predicate. The ``NormBridgeError::Overflow`` variant documents this
//! boundary.

use crypto_bigint::Uint;
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use shared_numfield::IntPoly;

// ─── Rational norm ────────────────────────────────────────────────────────────

/// Compute the rational norm ``N_rat(a, b) = a − b·m``.
///
/// This is the value on the rational (g-side) of NFS: ``g(x) = x − m``, so
/// ``g(a/b) · b = a − b·m``. The result is signed (can be negative when ``a < b·m``).
///
/// :param a: The a-coordinate of the sieve pair.
/// :param b: The b-coordinate of the sieve pair (positive).
/// :param m: The base-m value from the polynomial pair.
/// :returns: The rational norm ``a − b·m``.
pub fn rational_norm(a: &BigInt, b: &BigInt, m: &BigInt) -> BigInt {
    a - b * m
}

// ─── Algebraic norm ───────────────────────────────────────────────────────────

/// Compute the algebraic value ``F(a, b) = b^d · f(a/b)`` (named ``N_alg`` here).
///
/// This is the homogeneous form of ``f`` evaluated at ``(a, b)``; for ``f`` with leading
/// coefficient ``c`` and root ``α`` it equals ``c · N_{K/ℚ}(a − b·α)``:
///
/// ```text
/// N_alg(a, b) = Σ_{i=0}^{d} f.coeffs[i] · a^i · b^{d−i}
/// ```
///
/// Computed directly from coefficients to avoid rational arithmetic. The result is signed
/// (can be negative depending on ``f`` and ``(a, b)``).
///
/// # Mathematical note
///
/// ``F(a, b) = c · N_{K/ℚ}(a − b·α)``, where ``c`` is the leading coefficient of ``f``. It
/// equals ``Res(f, a − b·x)`` up to sign (the resultant is taken with ``f`` of degree ``d``
/// and ``a − b·x`` of degree 1, so ``Res(f, a − b·x) = ± F(a, b)``). Only for monic ``f``
/// (``c = 1``) is it the norm of the ideal ``(a − b·α)`` in ℤ[α]; for non-monic ``f`` the
/// ideal norm is ``F(a, b)/c`` and ``α`` is not an algebraic integer. Existing tests:
/// ``algebraic_norm_cubic_hand_computed``, ``algebraic_norm_with_b_greater_than_1`` and
/// ``algebraic_norm_matches_eval_when_b_is_one`` below (monic ``f``), and, for a non-monic
/// ``f = 3x² + x − 1``, ``kat_good_and_exceptional_relations_accepted`` in
/// ``tests/checked_relations_kat.rs``.
///
/// # Why not f.eval(a/b)?
///
/// ``f.eval(a/b)`` requires rational arithmetic (``a/b`` is generally not an integer). The
/// homogeneous form avoids this: each term ``f.coeffs[i] · a^i · b^{d−i}`` is an integer
/// product. This is the standard NFS implementation approach.
///
/// :param a: The a-coordinate of the sieve pair.
/// :param b: The b-coordinate of the sieve pair (positive).
/// :param f: The algebraic polynomial (from PolyPair).
/// :returns: The algebraic value ``b^d · f(a/b) = c · N_{K/ℚ}(a − b·α)``.
pub fn algebraic_norm(a: &BigInt, b: &BigInt, f: &IntPoly) -> BigInt {
    let d = match f.degree() {
        None => return BigInt::zero(), // zero polynomial
        Some(d) => d,
    };

    // Compute Σ f.coeffs[i] · a^i · b^{d−i} using Horner's method adapted for homogeneous form.
    // Direct summation: accumulate a^i and b^{d-i} powers.
    // We compute a_powers[i] = a^i and b_powers[j] = b^j, then sum.
    // For efficiency at toy scale, compute directly.
    let mut result = BigInt::zero();
    let mut a_pow = BigInt::from(1i64); // a^i, starting at a^0 = 1
    // b^{d-i}: start at b^d and divide by b each step — but b may be 0.
    // Instead, precompute b^{d-i} as b^d / b^i = b^{d-i}.
    // Simpler: compute b_pow_d_minus_i = b^{d-i} directly.
    // We iterate i from 0 to d; b^{d-i} decreases from b^d to b^0.
    // Precompute b^d, then divide by b each step (but b may be 0 or negative).
    // Safest: compute b^{d-i} = b^(d-i) for each i.

    for i in 0..=d {
        let coeff = f.coeffs.get(i).cloned().unwrap_or_else(BigInt::zero);
        if !coeff.is_zero() {
            // b^{d-i}
            let b_pow = pow_bigint(b, d - i);
            result += &coeff * &a_pow * &b_pow;
        }
        // Update a^i → a^{i+1}
        if i < d {
            a_pow *= a;
        }
    }

    result
}

// ─── Norm bridge ──────────────────────────────────────────────────────────────

/// Error type for norm-to-Uint conversion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NormBridgeError {
    /// The absolute value of the norm exceeds 256 bits (Uint<4> capacity).
    ///
    /// :param bits_required: The number of bits required to represent |norm|.
    Overflow { bits_required: usize },
}

impl std::fmt::Display for NormBridgeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Overflow { bits_required } => {
                write!(
                    f,
                    "norm absolute value requires {bits_required} bits, exceeding Uint<4> capacity \
                     (256 bits); shrink the sieve region or use a wider Uint<L>"
                )
            }
        }
    }
}

impl std::error::Error for NormBridgeError {}

/// Convert a signed BigInt norm to ``Uint<4>`` for ``trial_smooth``.
///
/// Returns the absolute value as ``Uint<4>``, or an error if ``|norm| ≥ 2^256``.
/// The sign is tracked separately in the ``Relation`` (for the −1 column in the linear algebra step).
///
/// # Design note (C1 resolution)
///
/// Toy-scale norms fit 256 bits per ROADMAP. If a chosen toy N/region overflows, shrink
/// the region or widen via the documented ``Uint<L>`` path. The bridge rejects out-of-range
/// norms rather than silently truncating (silent truncation would produce incorrect smoothness
/// witnesses).
///
/// :param norm: The signed norm value.
/// :returns: ``Ok(|norm|)`` as ``Uint<4>``, or ``Err(NormBridgeError::Overflow)`` if too large.
pub fn norm_to_uint(norm: &BigInt) -> Result<Uint<4>, NormBridgeError> {
    let abs_norm = norm.abs();

    // Check if the absolute value fits in 256 bits.
    // BigInt::bits() returns the number of bits in the magnitude (excluding sign).
    let bits = abs_norm.bits() as usize;
    if bits > 256 {
        return Err(NormBridgeError::Overflow { bits_required: bits });
    }

    // Convert BigInt magnitude to Uint<4>.
    // Extract the magnitude as a little-endian byte array and load into Uint<4>.
    let (_, bytes_be) = abs_norm.to_bytes_be();
    // Uint<4> is 32 bytes (256 bits). Pad to 32 bytes big-endian.
    if bytes_be.len() > 32 {
        // This shouldn't happen given the bits check above, but be defensive.
        return Err(NormBridgeError::Overflow { bits_required: bytes_be.len() * 8 });
    }

    let mut buf = [0u8; 32];
    let offset = 32 - bytes_be.len();
    buf[offset..].copy_from_slice(&bytes_be);

    Ok(Uint::<4>::from_be_slice(&buf))
}

/// Extract the sign of a norm: true if negative, false if non-negative.
///
/// :param norm: The signed norm value.
/// :returns: ``true`` if ``norm < 0``, ``false`` if ``norm ≥ 0``.
pub fn norm_sign(norm: &BigInt) -> bool {
    norm.is_negative()
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Compute ``base^exp`` for ``BigInt`` using repeated squaring.
fn pow_bigint(base: &BigInt, exp: usize) -> BigInt {
    if exp == 0 {
        return BigInt::from(1i64);
    }
    let mut result = BigInt::from(1i64);
    let mut b = base.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result *= &b;
        }
        b = &b * &b;
        e >>= 1;
    }
    result
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    /// f(x) = x³ − x − 1 (the classic NFS toy polynomial).
    fn f_cubic() -> IntPoly {
        IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)])
    }

    // ── rational_norm ────────────────────────────────────────────────────────

    #[test]
    fn rational_norm_positive() {
        // a=5, b=1, m=3: N_rat = 5 - 1*3 = 2
        assert_eq!(rational_norm(&bi(5), &bi(1), &bi(3)), bi(2));
    }

    #[test]
    fn rational_norm_negative() {
        // a=1, b=2, m=3: N_rat = 1 - 2*3 = -5
        assert_eq!(rational_norm(&bi(1), &bi(2), &bi(3)), bi(-5));
    }

    #[test]
    fn rational_norm_zero() {
        // a=6, b=2, m=3: N_rat = 6 - 2*3 = 0
        assert_eq!(rational_norm(&bi(6), &bi(2), &bi(3)), bi(0));
    }

    // ── algebraic_norm ───────────────────────────────────────────────────────

    #[test]
    fn algebraic_norm_cubic_hand_computed() {
        // f(x) = x³ − x − 1, a=2, b=1, d=3.
        // N_alg(2, 1) = 1·2³·1⁰ + 0·2²·1¹ + (-1)·2¹·1² + (-1)·2⁰·1³
        //             = 8 + 0 - 2 - 1 = 5
        // Equivalently: f(2) = 8 - 2 - 1 = 5 (since b=1, b^d = 1).
        let f = f_cubic();
        let norm = algebraic_norm(&bi(2), &bi(1), &f);
        assert_eq!(norm, bi(5));
    }

    #[test]
    fn algebraic_norm_with_b_greater_than_1() {
        // f(x) = x³ − x − 1, a=1, b=2, d=3.
        // N_alg(1, 2) = 1·1³·2⁰ + 0·1²·2¹ + (-1)·1¹·2² + (-1)·1⁰·2³
        //             = 1 + 0 - 4 - 8 = -11
        let f = f_cubic();
        let norm = algebraic_norm(&bi(1), &bi(2), &f);
        assert_eq!(norm, bi(-11));
    }

    #[test]
    fn algebraic_norm_matches_eval_when_b_is_one() {
        // When b=1, N_alg(a, 1) = f(a).
        let f = f_cubic();
        for a in -5i64..=5 {
            let norm = algebraic_norm(&bi(a), &bi(1), &f);
            let eval = f.eval(&bi(a));
            assert_eq!(norm, eval, "N_alg(a, 1) should equal f(a) for a={a}");
        }
    }

    // ── norm_to_uint ─────────────────────────────────────────────────────────

    #[test]
    fn norm_to_uint_positive() {
        let norm = bi(42);
        let u = norm_to_uint(&norm).expect("42 should fit in Uint<4>");
        assert_eq!(u, Uint::<4>::from(42u64));
    }

    #[test]
    fn norm_to_uint_negative_takes_abs() {
        let norm = bi(-42);
        let u = norm_to_uint(&norm).expect("-42 should fit in Uint<4> (abs = 42)");
        assert_eq!(u, Uint::<4>::from(42u64));
    }

    #[test]
    fn norm_to_uint_zero() {
        let norm = bi(0);
        let u = norm_to_uint(&norm).expect("0 should fit in Uint<4>");
        assert_eq!(u, Uint::<4>::ZERO);
    }

    #[test]
    fn norm_to_uint_overflow() {
        // Construct a number with 257 bits: 2^257.
        let big = BigInt::from(1i64) << 257;
        let result = norm_to_uint(&big);
        assert!(
            matches!(result, Err(NormBridgeError::Overflow { bits_required }) if bits_required > 256),
            "2^257 should overflow Uint<4>"
        );
    }

    #[test]
    fn norm_to_uint_exactly_256_bits_fits() {
        // 2^256 - 1 is the maximum Uint<4> value (256 bits, all ones).
        // 2^256 itself would require 257 bits and should overflow.
        let max_256 = (BigInt::from(1i64) << 256) - BigInt::from(1i64);
        let result = norm_to_uint(&max_256);
        assert!(result.is_ok(), "2^256 - 1 should fit in Uint<4>");
    }

    // ── norm_sign ────────────────────────────────────────────────────────────

    #[test]
    fn norm_sign_negative() {
        assert!(norm_sign(&bi(-1)));
        assert!(norm_sign(&bi(-100)));
    }

    #[test]
    fn norm_sign_non_negative() {
        assert!(!norm_sign(&bi(0)));
        assert!(!norm_sign(&bi(1)));
        assert!(!norm_sign(&bi(100)));
    }
}
