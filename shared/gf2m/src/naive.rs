//! Schoolbook GF(2^m) arithmetic over a polynomial basis, generic over the limb count.
//!
//! `F2mNaive<L>` stores a canonical coefficient bit-vector `c ∈ [0, 2^m)` as a
//! `Uint<L>` and performs all operations with straightforward shift-and-XOR logic.
//! It is the pedagogical baseline: correct, slow, easy to audit.
//!
//! # Storage contract
//!
//! **`Uint<L>` is a polynomial coefficient bit-vector, NOT an integer.**
//! Bit `i` of `c` is the coefficient of `x^i`.  Only XOR, shift, and bit
//! operations are meaningful.  Integer-arithmetic methods on `Uint<L>`
//! (`wrapping_add`, `rem`, `mul_wide`, etc.) are meaningless on GF(2)
//! coefficient vectors and are NEVER used here.  In particular, `mul_wide`
//! performs integer multiplication with carry — the carryless multiply uses
//! the comb algorithm (shift-and-XOR) instead.
//!
//! # Algorithms
//!
//! - **Add/Sub/Neg:** XOR (sub == add, neg == identity in char 2).
//! - **Carryless comb multiply:** for each set bit `i` of `rhs.c`, XOR
//!   `lhs.c` left-shifted by `i` into a `Uint<2L>` accumulator.  The
//!   intermediate has degree < 2m.
//! - **Modular reduction:** while the accumulator's degree ≥ m, XOR `poly`
//!   (widened to `Uint<2L>`) shifted to align its leading term with the
//!   accumulator's leading term.  Repeats until degree < m.
//! - **Frobenius/square:** insert a zero bit between each coefficient bit of
//!   `self.c` (bit-spread), producing a degree-<2m intermediate, then reduce.
//!
//! # Macro-generated impls
//!
//! The `mul` method requires widening to `Uint<{2*L}>`, which is not
//! expressible in stable Rust for arbitrary `L`.  The macro `impl_f2m_naive!`
//! generates a complete `impl F2m<$L> for F2mNaive<$L>` for each concrete
//! `($L, $DL)` pair where `$DL = 2 * $L`.  Pairs: (1,2), (2,4), (4,8), (8,16).

use crate::BitLen;
use crypto_bigint::Uint;

use super::F2m;

/// Schoolbook GF(2^m) element: canonical coefficient bit-vector stored as `Uint<L>`.
///
/// All arithmetic is performed in canonical form (no internal encoding).
/// This makes the implementation easy to audit at the cost of performance.
///
/// The irreducible polynomial is NOT stored on the struct — it is passed
/// per-call as `poly: &Uint<L>`, mirroring `FpNaive`'s per-call `p: &Uint<L>`.
/// Operations mixing elements reduced under different `poly` values are
/// meaningless; the caller is responsible for consistency.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct F2mNaive<const L: usize> {
    /// Canonical coefficient bit-vector: bit `i` is the coefficient of `x^i`.
    c: Uint<L>,
}

impl<const L: usize> F2mNaive<L> {
    /// Wrap a pre-reduced coefficient bit-vector (degree must be < degree of `poly`).
    #[inline]
    fn wrap(c: Uint<L>) -> Self {
        Self { c }
    }
}

/// Return the bit-length of a `Uint<L>` (position of the highest set bit + 1).
///
/// Returns 0 for the zero value.  This is the degree of the polynomial + 1.
/// Delegates to [`BitLen::bit_len`], since `Uint::bits()` returns `u32`.
#[inline]
pub(crate) fn uint_bits<const L: usize>(v: &Uint<L>) -> usize {
    v.bit_len()
}

/// Implement the complete `F2m<L>` trait for a specific limb count.
///
/// The `mul` and `square` methods require widening to `Uint<{2*L}>`, which is
/// not expressible in stable Rust for arbitrary `L`.  This macro generates a
/// complete `impl F2m<$L> for F2mNaive<$L>` for each concrete `($L, $DL)`
/// pair where `$DL = 2 * $L`.
macro_rules! impl_f2m_naive {
    ($L:literal, $DL:literal) => {
        impl F2m<$L> for F2mNaive<$L> {
            fn zero() -> Self {
                Self::wrap(Uint::<$L>::ZERO)
            }

            fn one() -> Self {
                Self::wrap(Uint::<$L>::ONE)
            }

            fn from_u64(v: u64, poly: &Uint<$L>) -> Self {
                let u = Uint::<$L>::from(v);
                Self::from_uint(u, poly)
            }

            fn from_uint(v: Uint<$L>, poly: &Uint<$L>) -> Self {
                // Reduce v mod poly if its degree >= degree of poly.
                // Degree of poly = poly.bit_len() - 1.  Degree of v = v.bit_len() - 1.
                // We need degree(v) < degree(poly), i.e. v.bit_len() <= poly.bit_len() - 1.
                let m = poly.bit_len(); // degree of poly + 1; degree of poly = m - 1
                if v.bit_len() < m {
                    // Already in range [0, 2^(m-1)); no reduction needed.
                    return Self::wrap(v);
                }
                // Reduce: widen to 2L and use the reduction routine.
                let wide = wide_from_low::<$L, $DL>(v);
                let poly_wide = wide_from_low::<$L, $DL>(*poly);
                let reduced = poly_reduce::<$DL>(wide, poly_wide, m - 1);
                let lo = low_half::<$L, $DL>(reduced);
                Self::wrap(lo)
            }

            fn to_uint(&self) -> Uint<$L> {
                self.c
            }

            fn add(&self, rhs: &Self) -> Self {
                // GF(2^m) addition = XOR of coefficient bit-vectors.
                // No reduction needed: XOR of two degree-<m polynomials has degree < m.
                Self::wrap(self.c ^ rhs.c)
            }

            /// Subtraction in GF(2^m) equals addition (XOR) — char-2 invariant.
            ///
            /// **Do NOT implement as `p − self`** (that is `Fp`'s convention).
            /// In characteristic 2, `a − b = a + b = a XOR b`.
            fn sub(&self, rhs: &Self) -> Self {
                // sub == add in characteristic 2.
                self.add(rhs)
            }

            /// Negation in GF(2^m) is the identity — char-2 invariant.
            ///
            /// **Do NOT implement as `p − self`** (that is `Fp`'s convention).
            /// In characteristic 2, `−a = a`.
            fn neg(&self) -> Self {
                // neg == identity in characteristic 2.
                self.clone()
            }

            fn mul(&self, rhs: &Self, poly: &Uint<$L>) -> Self {
                // Carryless comb multiply:
                // For each set bit i of rhs.c, XOR self.c left-shifted by i into
                // a Uint<2L> accumulator.  The intermediate has degree < 2m.
                //
                // IMPORTANT: this is carryless (XOR-based), NOT integer multiplication.
                // Do NOT use mul_wide here — that performs integer mul-with-carry.
                let m = poly.bit_len() - 1; // degree of the irreducible
                let rhs_bits = rhs.c.bit_len();
                let mut acc = Uint::<$DL>::ZERO;

                for i in 0..rhs_bits {
                    if rhs.c.bit_at(i) {
                        // XOR self.c shifted left by i into the accumulator.
                        let shifted = wide_shl::<$L, $DL>(self.c, i);
                        acc ^= shifted;
                    }
                }

                // Reduce the degree-<2m intermediate mod poly.
                let poly_wide = wide_from_low::<$L, $DL>(*poly);
                let reduced = poly_reduce::<$DL>(acc, poly_wide, m);
                let lo = low_half::<$L, $DL>(reduced);
                Self::wrap(lo)
            }

            fn square(&self, poly: &Uint<$L>) -> Self {
                // Frobenius bit-spread: insert a zero bit between each coefficient bit.
                // If self.c = Σ a_i x^i, then self.c² = Σ a_i x^(2i) in GF(2)[x]
                // (cross terms vanish because 2*a_i*a_j = 0 in char 2).
                // This is implemented by spreading each bit of self.c to every other
                // position in a Uint<2L>.
                let m = poly.bit_len() - 1;
                let self_bits = self.c.bit_len();
                let mut spread = Uint::<$DL>::ZERO;

                for i in 0..self_bits {
                    if self.c.bit_at(i) {
                        // Bit i of self maps to bit 2*i of the spread.
                        let bit_val = Uint::<$DL>::ONE;
                        let shifted = bit_val.shl_bits(2 * i);
                        spread ^= shifted;
                    }
                }

                // Reduce the degree-<2m intermediate mod poly.
                let poly_wide = wide_from_low::<$L, $DL>(*poly);
                let reduced = poly_reduce::<$DL>(spread, poly_wide, m);
                let lo = low_half::<$L, $DL>(reduced);
                Self::wrap(lo)
            }

            fn frobenius(&self, poly: &Uint<$L>) -> Self {
                // Frobenius endomorphism a → a² is exactly squaring in char 2.
                self.square(poly)
            }

            fn pow(&self, exp: &Uint<$L>, poly: &Uint<$L>) -> Self {
                // Square-and-multiply (right-to-left binary method).
                let mut result = Self::one();
                let mut base = self.clone();
                let mut e = *exp;
                while e != Uint::<$L>::ZERO {
                    if e.bit_at(0) {
                        result = result.mul(&base, poly);
                    }
                    base = base.square(poly);
                    e >>= 1;
                }
                result
            }

            fn trace(&self, poly: &Uint<$L>) -> Self {
                // Absolute trace: Tr(a) = Σ_{i=0}^{m-1} a^(2^i) in GF(2^m).
                //
                // Computed by iterating Frobenius (squaring) m-1 times and XOR-summing.
                // The result is always 0 or 1 (an element of GF(2) ⊂ GF(2^m)).
                let m = poly.bit_len() - 1; // degree of the irreducible
                let mut acc = self.clone();
                let mut cur = self.clone();
                for _ in 1..m {
                    cur = cur.frobenius(poly);
                    acc = acc.add(&cur);
                }
                acc
            }

            fn solve_quadratic(c: &Self, poly: &Uint<$L>) -> Self {
                // Solve x² + x = c in GF(2^m) (Artin–Schreier equation).
                //
                // Solvability precondition: trace(c) = 0.  The caller is responsible
                // for ensuring this; if trace(c) ≠ 0, the result is meaningless.
                //
                // For odd m: the half-trace H(c) = Σ_{i=0}^{(m-1)/2} c^(2^(2i)) is
                // a solution.  This is the standard closed form (HMV §2.3.6).
                //
                // For even m: no closed-form half-trace exists.  We use brute-force
                // search over the field.  This is correct and auditable at toy scale
                // (m ≤ 8 in the KATs); production use should target odd-m fields.
                let m = poly.bit_len() - 1; // degree of the irreducible
                if m % 2 == 1 {
                    // Odd m: half-trace formula.
                    let mut acc = c.clone(); // c^(2^0)
                    let mut cur = c.clone();
                    for _ in 0..(m - 1) / 2 {
                        cur = cur.frobenius(poly).frobenius(poly); // advance by 2
                        acc = acc.add(&cur);
                    }
                    acc
                } else {
                    // Even m: brute-force search (toy fields only).
                    // Find x such that x² + x = c.
                    let field_size = 1u64 << m;
                    for v in 0..field_size {
                        let x = Self::from_u64(v, poly);
                        let lhs = x.square(poly).add(&x);
                        if lhs == *c {
                            return x;
                        }
                    }
                    // No solution exists (trace(c) ≠ 0); return zero as sentinel.
                    Self::zero()
                }
            }

            fn inv(&self, poly: &Uint<$L>) -> Self {
                // Delegate to the extended-Euclidean baseline in inv.rs.
                // Panics if self is zero (no inverse exists).
                crate::inv::ext_euclid_inv(self, poly)
            }

            fn div(&self, rhs: &Self, poly: &Uint<$L>) -> Self {
                // a / b = a * b⁻¹.  Panics if rhs is zero.
                self.mul(&rhs.inv(poly), poly)
            }
        }
    };
}

// ── Helper functions for double-width arithmetic ──────────────────────────────

/// Widen a `Uint<L>` to `Uint<DL>` by zero-extending into the low half.
#[inline]
pub(crate) fn wide_from_low<const L: usize, const DL: usize>(v: Uint<L>) -> Uint<DL> {
    // Copy the L words of v into the low L words of a DL-word Uint.
    let src = v.as_words();
    let mut dst = [0u64; DL];
    dst[..L].copy_from_slice(src);
    Uint::<DL>::from_words(dst)
}

/// Extract the low `L` words of a `Uint<DL>`.
#[inline]
pub(crate) fn low_half<const L: usize, const DL: usize>(v: Uint<DL>) -> Uint<L> {
    let words = v.as_words();
    let mut lo = [0u64; L];
    lo.copy_from_slice(&words[..L]);
    Uint::<L>::from_words(lo)
}

/// Left-shift a `Uint<L>` by `shift` bits, producing a `Uint<DL>`.
///
/// Used in the carryless comb multiply to shift `self.c` into the double-width
/// accumulator without losing high bits.
#[inline]
pub(crate) fn wide_shl<const L: usize, const DL: usize>(v: Uint<L>, shift: usize) -> Uint<DL> {
    let wide = wide_from_low::<L, DL>(v);
    wide.shl_bits(shift)
}

/// Reduce a `Uint<DL>` polynomial modulo a `Uint<DL>` irreducible.
///
/// `poly_wide` is the irreducible polynomial widened to `DL` limbs.
/// `m` is the degree of the irreducible (so `poly_wide.bit_len() == m + 1`).
///
/// Algorithm: while the accumulator's degree ≥ m, XOR `poly_wide` shifted
/// to align its leading term (at bit `m`) with the accumulator's leading term
/// (at bit `deg`).  Each XOR reduces the degree by at least 1.  Terminates
/// when degree < m.
#[inline]
pub(crate) fn poly_reduce<const DL: usize>(mut acc: Uint<DL>, poly_wide: Uint<DL>, m: usize) -> Uint<DL> {
    loop {
        let deg = uint_bits(&acc);
        if deg == 0 || deg <= m {
            break;
        }
        // Shift poly_wide so its leading bit (at position m) aligns with
        // the accumulator's leading bit (at position deg - 1).
        let shift = deg - 1 - m;
        let shifted_poly = poly_wide.shl_bits(shift);
        acc ^= shifted_poly;
    }
    acc
}

// ── Generate impls for each supported limb count ──────────────────────────────

impl_f2m_naive!(1, 2);
impl_f2m_naive!(2, 4);
impl_f2m_naive!(4, 8);
impl_f2m_naive!(8, 16);

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// GF(2^4) irreducible: x⁴ + x + 1 = 0b10011 = 0x13.
    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64)
    }

    /// GF(2^8) AES irreducible: x⁸ + x⁴ + x³ + x + 1 = 0x11b.
    fn poly8() -> Uint<1> {
        Uint::<1>::from(0x11bu64)
    }

    #[test]
    fn add_is_xor() {
        let p = poly4();
        let a = F2mNaive::<1>::from_u64(0b0101, &p);
        let b = F2mNaive::<1>::from_u64(0b0011, &p);
        // 0101 XOR 0011 = 0110
        assert_eq!(a.add(&b).to_uint(), Uint::<1>::from(0b0110u64));
    }

    #[test]
    fn sub_equals_add() {
        // char-2 invariant: sub == add
        let p = poly4();
        let a = F2mNaive::<1>::from_u64(0b1010, &p);
        let b = F2mNaive::<1>::from_u64(0b0110, &p);
        assert_eq!(a.sub(&b), a.add(&b));
    }

    #[test]
    fn neg_is_identity() {
        // char-2 invariant: neg == identity
        let p = poly4();
        let a = F2mNaive::<1>::from_u64(0b1010, &p);
        assert_eq!(a.neg(), a);
    }

    #[test]
    fn add_self_is_zero() {
        // a + a = 0 in char 2
        let p = poly4();
        let a = F2mNaive::<1>::from_u64(0b1010, &p);
        assert_eq!(a.add(&a), F2mNaive::<1>::zero());
    }

    #[test]
    fn mul_by_one_is_identity() {
        let p = poly4();
        let a = F2mNaive::<1>::from_u64(0b0111, &p);
        let one = F2mNaive::<1>::one();
        assert_eq!(a.mul(&one, &p), a);
    }

    #[test]
    fn mul_gf4_known() {
        // In GF(2^4) with x⁴+x+1:
        // x * x = x² (no reduction needed, degree 2 < 4)
        let p = poly4();
        let x = F2mNaive::<1>::from_u64(0b0010, &p); // x
        let x2 = x.mul(&x, &p);
        assert_eq!(x2.to_uint(), Uint::<1>::from(0b0100u64)); // x²

        // x² * x² = x⁴ ≡ x + 1 (mod x⁴+x+1)
        let x4 = x2.mul(&x2, &p);
        assert_eq!(x4.to_uint(), Uint::<1>::from(0b0011u64)); // x + 1
    }

    #[test]
    fn square_equals_mul_self() {
        let p = poly8();
        for v in [0u64, 1, 2, 5, 0x53, 0xab, 0xff] {
            let a = F2mNaive::<1>::from_u64(v, &p);
            assert_eq!(
                a.square(&p),
                a.mul(&a, &p),
                "square != mul(self,self) for v={v:#x}"
            );
        }
    }

    #[test]
    fn frobenius_equals_square() {
        let p = poly8();
        let a = F2mNaive::<1>::from_u64(0x53, &p);
        assert_eq!(a.frobenius(&p), a.square(&p));
    }

    #[test]
    fn from_uint_reduces() {
        // A value with degree >= m should be reduced.
        let p = poly4(); // degree 4
        // 0b10000 has degree 4 = m, so it should be reduced mod x⁴+x+1.
        // x⁴ ≡ x + 1 (mod x⁴+x+1), so 0b10000 → 0b0011.
        let a = F2mNaive::<1>::from_uint(Uint::<1>::from(0b10000u64), &p);
        assert_eq!(a.to_uint(), Uint::<1>::from(0b0011u64));
    }

    #[test]
    fn aes_mul_known() {
        // AES MixColumns uses GF(2^8) with x⁸+x⁴+x³+x+1.
        // Known: 0x53 * 0xca = 0x01 (they are inverses in AES field).
        let p = poly8();
        let a = F2mNaive::<1>::from_u64(0x53, &p);
        let b = F2mNaive::<1>::from_u64(0xca, &p);
        let prod = a.mul(&b, &p);
        assert_eq!(prod.to_uint(), Uint::<1>::ONE, "0x53 * 0xca should be 1 in AES GF(2^8)");
    }
}
