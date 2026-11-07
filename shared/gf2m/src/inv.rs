//! GF(2^m) multiplicative inversion: extended-Euclidean baseline and Itoh–Tsujii.
//!
//! Two independent algorithms are provided so they can be cross-checked against
//! each other in the KATs (a bug in `square`/`frobenius` shows up in Itoh–Tsujii
//! but not in extended-Euclidean, and vice versa for a bug in polynomial division).
//!
//! # Extended Euclidean over GF(2)[x]
//!
//! [`ext_euclid_inv`] computes `a⁻¹` via the extended polynomial GCD of `a(x)`
//! with the irreducible `poly(x)`.  Since `poly` is irreducible and `a ≠ 0`,
//! `gcd(a, poly) = 1`, and the Bézout cofactor `s` satisfies `a·s ≡ 1 (mod poly)`.
//!
//! Algorithm (all arithmetic in GF(2)[x], i.e. XOR for addition):
//! - Maintain `(r0, r1) = (poly, a)` and `(s0, s1) = (0, 1)`.
//! - At each step: divide `r0` by `r1` to get quotient `q` and remainder `r`.
//! - Update: `(r0, r1) ← (r1, r)` and `(s0, s1) ← (s1, s0 XOR q·s1)`.
//! - When `r1 = 1`, `s1` is the inverse.
//!
//! # Itoh–Tsujii via the Frobenius tower
//!
//! [`itoh_tsujii_inv`] computes `a⁻¹ = a^(2^m − 2)` using the Frobenius tower.
//! The key identity is `2^m − 2 = 2·(2^(m−1) − 1)`.  We build up the exponent
//! `2^(m−1) − 1` via the binary representation of `m−1`, using the fact that
//! `a^(2^k)` is computed by applying `frobenius` (squaring) `k` times.
//!
//! Concretely, for exponent `e = m − 1` written in binary `e = Σ eᵢ·2^i`:
//! - Maintain an accumulator `acc` and a running power `frob_k = a^(2^k)`.
//! - For each set bit `i` of `e`: multiply `acc` by `a^(2^(2^i − 1))` (built
//!   from the Frobenius tower).
//! - Finally square once: `a^(2^m − 2) = (a^(2^(m−1) − 1))^2`.
//!
//! Simpler equivalent used here: compute `a^(2^m − 2)` directly via
//! square-and-multiply on the exponent `2^m − 2`, but use `frobenius` (iterated
//! squaring) to compute each `a^(2^k)` step.  This is pedagogically the
//! Frobenius-tower approach: the exponent is decomposed into a sum of powers of 2,
//! each power-of-2 exponentiation is a chain of Frobenius maps, and the results
//! are combined by multiplication.

use crate::BitLen;
use crypto_bigint::Uint;

use crate::naive::F2mNaive;
use crate::F2m;

// ── Extended Euclidean inversion ──────────────────────────────────────────────

/// Compute `a⁻¹` in GF(2^m) via the extended Euclidean algorithm over GF(2)[x].
///
/// This is the auditable baseline: correct, O(m²), easy to verify by inspection.
///
/// # Panics
///
/// Panics if `a` is zero (zero has no multiplicative inverse).
pub fn ext_euclid_inv<const L: usize>(a: &F2mNaive<L>, poly: &Uint<L>) -> F2mNaive<L>
where
    F2mNaive<L>: F2m<L>,
{
    assert!(!a.is_zero(), "ext_euclid_inv: zero has no multiplicative inverse");

    // Work in Uint<L>; all arithmetic is XOR (GF(2)[x] addition).
    // r0 = poly, r1 = a; s0 = 0, s1 = 1.
    let mut r0: Uint<L> = *poly;
    let mut r1: Uint<L> = a.to_uint();
    let mut s0: Uint<L> = Uint::<L>::ZERO;
    let mut s1: Uint<L> = Uint::<L>::ONE;

    // Invariant: r0 = poly·(something) + a·s0  (mod 2)
    //            r1 = poly·(something) + a·s1  (mod 2)
    // When r1 = 1, s1 is the inverse of a mod poly.
    loop {
        // Check termination at the top: when r1 = 1, s1 is the inverse.
        if r1 == Uint::<L>::ONE {
            return F2mNaive::<L>::from_uint(s1, poly);
        }

        // Safety: r1 should never reach zero if a is non-zero and poly is irreducible.
        debug_assert!(
            r1 != Uint::<L>::ZERO,
            "ext_euclid_inv: gcd != 1 — input is zero or poly is reducible"
        );

        // Polynomial division of r0 by r1 over GF(2): compute quotient q and remainder r.
        let (q, r) = poly_divmod(r0, r1);

        // Update: (r0, r1) ← (r1, r)
        r0 = r1;
        r1 = r;

        // Update: (s0, s1) ← (s1, s0 XOR q·s1)
        // q * s1 in GF(2)[x] (no field reduction — this is polynomial arithmetic).
        let qs1 = poly_mul_reduce_wide(q, s1);
        let new_s1 = s0 ^ qs1;
        s0 = s1;
        s1 = new_s1;
    }
}

/// Polynomial division over GF(2)[x]: divide `a` by `b`, returning `(quotient, remainder)`.
///
/// Both `a` and `b` are coefficient bit-vectors.  All arithmetic is XOR.
/// Requires `b ≠ 0`.
fn poly_divmod<const L: usize>(mut a: Uint<L>, b: Uint<L>) -> (Uint<L>, Uint<L>) {
    debug_assert!(b != Uint::<L>::ZERO, "poly_divmod: divisor is zero");

    let deg_b = b.bit_len().saturating_sub(1); // degree of b
    let mut q = Uint::<L>::ZERO;

    loop {
        let bits_a = a.bit_len();
        if bits_a == 0 {
            break; // a = 0, remainder is 0
        }
        let deg_a = bits_a - 1;
        if deg_a < deg_b {
            break; // degree of remainder < degree of divisor: done
        }
        let shift = deg_a - deg_b;
        // q += x^shift (set bit `shift` of quotient)
        q ^= Uint::<L>::ONE.shl_bits(shift);
        // a -= b * x^shift = a XOR (b << shift)
        a ^= b.shl_bits(shift);
    }

    (q, a)
}

/// Multiply two polynomials over GF(2)[x] without reduction.
///
/// Both inputs fit in `Uint<L>`.  The product may have degree up to
/// `deg(a) + deg(b)`, which also fits in `Uint<L>` for the Bézout cofactor
/// computation (the cofactor degree is bounded by `deg(poly) - 1 < L*64`).
///
/// This is a carryless multiply (XOR-based), NOT integer multiplication.
fn poly_mul_reduce_wide<const L: usize>(a: Uint<L>, b: Uint<L>) -> Uint<L> {
    // Schoolbook carryless multiply: for each set bit i of b, XOR a<<i into acc.
    // The result fits in Uint<L> because the Bézout cofactor stays below degree m.
    let mut acc = Uint::<L>::ZERO;
    let b_bits = b.bit_len();
    for i in 0..b_bits {
        if b.bit_at(i) {
            acc ^= a.shl_bits(i);
        }
    }
    acc
}

// ── Itoh–Tsujii inversion ─────────────────────────────────────────────────────

/// Compute `a⁻¹` in GF(2^m) via the Itoh–Tsujii algorithm.
///
/// Uses the identity `a⁻¹ = a^(2^m − 2)` and the Frobenius tower to compute
/// the exponentiation with O(log m) multiplications.
///
/// The exponent `2^m − 2 = 2·(2^(m−1) − 1)` is decomposed as follows:
/// - Write `e = m − 1` in binary.
/// - Build `a^(2^e − 1)` using the addition chain for `e` via the Frobenius tower.
/// - Square once to get `a^(2^m − 2) = a^(2·(2^(m−1) − 1))`.
///
/// The addition chain for `2^e − 1` uses the identity:
/// `2^(i+j) − 1 = (2^i − 1)·2^j + (2^j − 1)`, so
/// `a^(2^(i+j) − 1) = a^(2^(i+j) − 1)` can be built from
/// `a^(2^i − 1)` and `a^(2^j − 1)` with one multiplication and `j` Frobenius maps.
///
/// For simplicity and correctness, we use the square-and-multiply approach on
/// the exponent `2^m − 2` directly, but compute each `a^(2^k)` via iterated
/// Frobenius maps (not via `pow`).  This is the Frobenius-tower decomposition.
///
/// # Panics
///
/// Panics if `a` is zero.
pub fn itoh_tsujii_inv<const L: usize>(a: &F2mNaive<L>, poly: &Uint<L>) -> F2mNaive<L>
where
    F2mNaive<L>: F2m<L>,
{
    assert!(!a.is_zero(), "itoh_tsujii_inv: zero has no multiplicative inverse");

    let m = poly.bit_len() - 1; // degree of the irreducible = field degree

    // Special case: GF(2^1) — only element is 1, its own inverse.
    if m == 1 {
        return a.clone();
    }

    // We compute a^(2^m − 2) using the Frobenius tower.
    //
    // Strategy: build a^(2^(m-1) - 1) using the binary method on the exponent (m-1),
    // then square once.
    //
    // The binary method for a^(2^e - 1) where e = m-1:
    //   Write e in binary: e = Σ eᵢ·2^i (from MSB to LSB, skipping leading 1).
    //   Maintain `acc = a^(2^k - 1)` where k grows from 1 to e.
    //   For each bit of e (after the leading 1):
    //     - Double: acc ← acc^(2^k) · acc = a^(2^(2k) - 1)
    //       (apply k Frobenius maps to acc, then multiply by acc)
    //     - If bit is 1: acc ← acc · a^(2^(2k)) = a^(2^(2k+1) - 1)
    //       (apply one more Frobenius map to acc, then multiply by a)
    //
    // This gives a^(2^(m-1) - 1) in O(log m) multiplications and O(m) Frobenius maps.

    let e = m - 1; // exponent for the tower: we want a^(2^e - 1)

    // Find the bit length of e (position of the highest set bit).
    let e_bits = usize::BITS as usize - e.leading_zeros() as usize;

    // Start with acc = a^(2^1 - 1) = a^1 = a, k = 1.
    let mut acc = a.clone();
    let mut k: usize = 1; // acc = a^(2^k - 1)

    // Process bits of e from the second-highest bit down to bit 0.
    for bit_pos in (0..e_bits - 1).rev() {
        // Double step: acc ← acc^(2^k) · acc = a^(2^(2k) - 1)
        // Apply k Frobenius maps to acc to get acc^(2^k).
        let mut frob_acc = acc.clone();
        for _ in 0..k {
            frob_acc = frob_acc.frobenius(poly);
        }
        acc = frob_acc.mul(&acc, poly);
        k *= 2; // now acc = a^(2^k - 1) with new k = 2 * old k

        // If bit bit_pos of e is set, do an extra step:
        // acc ← acc · a^(2^k) = a^(2^(k+1) - 1)
        if (e >> bit_pos) & 1 == 1 {
            // Compute a^(2^k) by applying k Frobenius maps to a.
            let mut frob_a = a.clone();
            for _ in 0..k {
                frob_a = frob_a.frobenius(poly);
            }
            acc = acc.mul(&frob_a, poly);
            k += 1; // now acc = a^(2^k - 1) with new k = old k + 1
        }
    }

    // Now acc = a^(2^(m-1) - 1).  Square once to get a^(2^m - 2).
    acc.frobenius(poly)
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::F2m;

    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64) // x⁴ + x + 1
    }

    fn poly8() -> Uint<1> {
        Uint::<1>::from(0x11bu64) // x⁸ + x⁴ + x³ + x + 1 (AES)
    }

    #[test]
    fn ext_euclid_one_is_own_inverse() {
        let p = poly4();
        let one = F2mNaive::<1>::one();
        assert_eq!(ext_euclid_inv(&one, &p), one);
    }

    #[test]
    fn itoh_tsujii_one_is_own_inverse() {
        let p = poly4();
        let one = F2mNaive::<1>::one();
        assert_eq!(itoh_tsujii_inv(&one, &p), one);
    }

    #[test]
    fn ext_euclid_aes_known() {
        // 0x53 and 0xca are inverses in AES GF(2^8).
        let p = poly8();
        let a = F2mNaive::<1>::from_u64(0x53, &p);
        let b = F2mNaive::<1>::from_u64(0xca, &p);
        assert_eq!(ext_euclid_inv(&a, &p), b);
        assert_eq!(ext_euclid_inv(&b, &p), a);
    }

    #[test]
    fn itoh_tsujii_aes_known() {
        let p = poly8();
        let a = F2mNaive::<1>::from_u64(0x53, &p);
        let b = F2mNaive::<1>::from_u64(0xca, &p);
        assert_eq!(itoh_tsujii_inv(&a, &p), b);
        assert_eq!(itoh_tsujii_inv(&b, &p), a);
    }

    #[test]
    fn both_algorithms_agree_gf4() {
        let p = poly4();
        for v in 1u64..16 {
            let a = F2mNaive::<1>::from_u64(v, &p);
            let ee = ext_euclid_inv(&a, &p);
            let it = itoh_tsujii_inv(&a, &p);
            assert_eq!(ee, it, "algorithms disagree for v={v:#x} in GF(2^4)");
        }
    }
}
