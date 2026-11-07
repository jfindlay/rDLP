//! Optimized GF(2^m) arithmetic: single-level Karatsuba carryless multiplier.
//!
//! `F2mOpt<L>` is the optimized analogue of [`F2mNaive`]: same canonical
//! coefficient-bit-vector storage (`{ c: Uint<L> }`), same trait surface, but
//! with a Karatsuba-split `mul` instead of the schoolbook comb.
//!
//! # Algorithm: single-level Karatsuba split
//!
//! For a degree-`m` field element (`m = poly.bit_len() − 1`), split each operand
//! at the midpoint `h = m / 2`:
//!
//! ```text
//! a = a_lo + x^h · a_hi
//! b = b_lo + x^h · b_hi
//! ```
//!
//! The product is:
//!
//! ```text
//! a·b = a_lo·b_lo  +  x^h · mid  +  x^(2h) · a_hi·b_hi
//! ```
//!
//! where `mid = (a_lo + a_hi)·(b_lo + b_hi) XOR a_lo·b_lo XOR a_hi·b_hi`
//! (in characteristic 2, subtraction = addition = XOR).
//!
//! The three half-products (`lo·lo`, `hi·hi`, `(lo+hi)·(lo+hi)`) are computed
//! with the schoolbook comb at half-width, then assembled into a degree-<2m
//! intermediate in `Uint<DL>` and reduced mod `poly` via the same `poly_reduce`
//! routine used by `F2mNaive`.
//!
//! # Equivalence
//!
//! `F2mOpt` produces byte-identical results to `F2mNaive` on every input.
//! The naive↔optimized agreement KAT in `tests/gf2m_kat.rs` verifies this
//! exhaustively for GF(2^4) and GF(2^8).
//!
//! # pclmulqdq note
//!
//! A crypto-scale implementation would use the x86 carryless-multiply
//! intrinsic `core::arch::x86_64::_mm_clmulepi64_si128`, gated behind
//! `#[cfg(target_feature = "pclmulqdq")]` on an `unsafe fn`.  That path is
//! **omitted here** to preserve the crate's `unsafe_code = "forbid"` invariant:
//! at toy scale (m ≤ 8 in the KATs, ≤ 64 in the type) the intrinsic buys
//! near-zero observable demonstration value while the cost is a real erosion of
//! the zero-unsafe guarantee every `shared/` crate upholds.  The software
//! Karatsuba already exhibits the subquadratic-split algorithm; `pclmulqdq`
//! would exhibit only an x86 micro-architectural detail.  Adding the gated path
//! later (behind `#[cfg(target_feature = "pclmulqdq")]` + a localized
//! `#[allow(unsafe_code)]`) is a purely additive, reversible change.
//!
//! # López–Dahab comb note
//!
//! The `F2mNaive` schoolbook comb **is** the comb method (López–Dahab without windowing).
//! A windowed-comb variant would be a within-baseline micro-optimization at toy scale —
//! algorithmically the same speed-up story as Karatsuba but without an additional equivalence
//! contract.  It is deliberately not separately implemented here.

use crate::BitLen;
use crypto_bigint::Uint;

use super::F2m;
use crate::naive::{low_half, poly_reduce, uint_bits, wide_from_low, wide_shl};

/// Optimized GF(2^m) element: same canonical storage as `F2mNaive`, Karatsuba `mul`.
///
/// The coefficient bit-vector `c` is identical in representation to `F2mNaive::c`.
/// Only the `mul` implementation differs — all other operations are identical to
/// `F2mNaive`.  This means `F2mOpt::to_uint()` and `F2mNaive::to_uint()` return
/// byte-identical values for the same field element, enabling the naive↔optimized
/// agreement KAT.
///
/// **Storage contract:** `Uint<L>` is a polynomial coefficient bit-vector, NOT an
/// integer.  Only XOR, shift, and bit operations are meaningful.  See the crate-level
/// documentation for the full warning.
///
/// The irreducible polynomial is NOT stored on the struct — it is passed per-call
/// as `poly: &Uint<L>`, mirroring `F2mNaive` and `FpNaive`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct F2mOpt<const L: usize> {
    /// Canonical coefficient bit-vector: bit `i` is the coefficient of `x^i`.
    c: Uint<L>,
}

impl<const L: usize> F2mOpt<L> {
    /// Wrap a pre-reduced coefficient bit-vector (degree must be < degree of `poly`).
    #[inline]
    fn wrap(c: Uint<L>) -> Self {
        Self { c }
    }
}

/// Schoolbook carryless comb multiply of two `Uint<L>` values, producing `Uint<L>`.
///
/// Both inputs are assumed to have degree < `half_m` (i.e., they fit in `half_m`
/// bits), so the product has degree < `2 * half_m ≤ m < L * 64` and fits in
/// `Uint<L>` without widening.
///
/// This is the half-width multiply used by the Karatsuba recombination step.
/// It is the same schoolbook comb as `F2mNaive::mul`, but operating entirely
/// within `Uint<L>` (no `Uint<DL>` widening needed because the product fits).
#[inline]
fn comb_mul_narrow<const L: usize>(a: Uint<L>, b: Uint<L>) -> Uint<L> {
    let b_bits = uint_bits(&b);
    let mut acc = Uint::<L>::ZERO;
    for i in 0..b_bits {
        if b.bit_at(i) {
            acc ^= a.shl_bits(i);
        }
    }
    acc
}

/// Implement the complete `F2m<L>` trait for a specific limb count.
///
/// Mirrors `impl_f2m_naive!` exactly, but with a Karatsuba `mul` body.
/// All other methods are identical to `F2mNaive`'s implementations.
///
/// The `square`/`frobenius` bit-spread is kept identical to naive — it is
/// already the optimal bit-spread; Karatsuba does not improve squaring.
macro_rules! impl_f2m_opt {
    ($L:literal, $DL:literal) => {
        impl F2m<$L> for F2mOpt<$L> {
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
                let m = poly.bit_len();
                if v.bit_len() < m {
                    return Self::wrap(v);
                }
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
                Self::wrap(self.c ^ rhs.c)
            }

            /// Subtraction in GF(2^m) equals addition (XOR) — char-2 invariant.
            fn sub(&self, rhs: &Self) -> Self {
                self.add(rhs)
            }

            /// Negation in GF(2^m) is the identity — char-2 invariant.
            fn neg(&self) -> Self {
                self.clone()
            }

            fn mul(&self, rhs: &Self, poly: &Uint<$L>) -> Self {
                // Single-level Karatsuba carryless multiply.
                //
                // Let m = degree(poly), h = m / 2 (the split point).
                // Split: a = a_lo + x^h * a_hi, b = b_lo + x^h * b_hi.
                //
                // Karatsuba recombination (char 2, so sub == add == XOR):
                //   a*b = lo_lo + x^h * mid + x^(2h) * hi_hi
                // where:
                //   lo_lo = a_lo * b_lo
                //   hi_hi = a_hi * b_hi
                //   mid   = (a_lo + a_hi) * (b_lo + b_hi) XOR lo_lo XOR hi_hi
                //
                // The three half-products fit in Uint<L> (degree < m).
                // Assembly into the full degree-<2m product uses Uint<DL>.
                let m = poly.bit_len() - 1; // degree of the irreducible
                let h = m / 2; // split point

                // Build the half-mask: bits 0..h-1 set.
                // For h < 64 (always true when L=1, m ≤ 64), this is (1 << h) - 1.
                // We compute it as a Uint<L> by shifting ONE left by h then subtracting 1.
                // In crypto-bigint, wrapping_sub on Uint is integer subtraction — safe here
                // because (1 << h) is never zero for h > 0, and for h == 0 the mask is 0.
                let mask: Uint<$L> = if h == 0 {
                    Uint::<$L>::ZERO
                } else {
                    Uint::<$L>::ONE.shl_bits(h).wrapping_sub(&Uint::<$L>::ONE)
                };

                // Split operands.
                let a_lo = self.c & mask;
                let a_hi = self.c.shr_bits(h);
                let b_lo = rhs.c & mask;
                let b_hi = rhs.c.shr_bits(h);

                // Three half-products (all fit in Uint<L>: degree < 2*h ≤ m).
                let lo_lo = comb_mul_narrow::<$L>(a_lo, b_lo);
                let hi_hi = comb_mul_narrow::<$L>(a_hi, b_hi);
                let mid_in = comb_mul_narrow::<$L>(a_lo ^ a_hi, b_lo ^ b_hi);
                let mid = mid_in ^ lo_lo ^ hi_hi;

                // Assemble the degree-<2m product in Uint<DL>:
                //   result = lo_lo + x^h * mid + x^(2h) * hi_hi
                // In char 2, + is XOR.
                let lo_lo_wide = wide_from_low::<$L, $DL>(lo_lo);
                let mid_wide: Uint<$DL> = wide_shl::<$L, $DL>(mid, h);
                // 2h may equal m (when m is even) or m-1 (when m is odd).
                // Either way 2h ≤ m < 2m, so the shift fits in Uint<DL>.
                let hi_hi_wide: Uint<$DL> = wide_shl::<$L, $DL>(hi_hi, 2 * h);

                let acc = lo_lo_wide ^ mid_wide ^ hi_hi_wide;

                // Reduce mod poly.
                let poly_wide = wide_from_low::<$L, $DL>(*poly);
                let reduced = poly_reduce::<$DL>(acc, poly_wide, m);
                let lo = low_half::<$L, $DL>(reduced);
                Self::wrap(lo)
            }

            fn square(&self, poly: &Uint<$L>) -> Self {
                // Frobenius bit-spread: identical to F2mNaive::square.
                // Karatsuba does not improve squaring — the bit-spread is already optimal.
                let m = poly.bit_len() - 1;
                let self_bits = self.c.bit_len();
                let mut spread = Uint::<$DL>::ZERO;

                for i in 0..self_bits {
                    if self.c.bit_at(i) {
                        let bit_val = Uint::<$DL>::ONE;
                        let shifted = bit_val.shl_bits(2 * i);
                        spread ^= shifted;
                    }
                }

                let poly_wide = wide_from_low::<$L, $DL>(*poly);
                let reduced = poly_reduce::<$DL>(spread, poly_wide, m);
                let lo = low_half::<$L, $DL>(reduced);
                Self::wrap(lo)
            }

            fn frobenius(&self, poly: &Uint<$L>) -> Self {
                self.square(poly)
            }

            fn pow(&self, exp: &Uint<$L>, poly: &Uint<$L>) -> Self {
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
                // Identical algorithm to F2mNaive::trace — iterate Frobenius and XOR-sum.
                let m = poly.bit_len() - 1;
                let mut acc = self.clone();
                let mut cur = self.clone();
                for _ in 1..m {
                    cur = cur.frobenius(poly);
                    acc = acc.add(&cur);
                }
                acc
            }

            fn solve_quadratic(c: &Self, poly: &Uint<$L>) -> Self {
                // Solve x² + x = c in GF(2^m).  Identical algorithm to F2mNaive.
                // For odd m: half-trace formula.  For even m: brute-force search.
                let m = poly.bit_len() - 1;
                if m % 2 == 1 {
                    let mut acc = c.clone();
                    let mut cur = c.clone();
                    for _ in 0..(m - 1) / 2 {
                        cur = cur.frobenius(poly).frobenius(poly);
                        acc = acc.add(&cur);
                    }
                    acc
                } else {
                    let field_size = 1u64 << m;
                    for v in 0..field_size {
                        let x = Self::from_u64(v, poly);
                        let lhs = x.square(poly).add(&x);
                        if lhs == *c {
                            return x;
                        }
                    }
                    Self::zero()
                }
            }

            fn inv(&self, poly: &Uint<$L>) -> Self {
                // Delegate to the extended-Euclidean baseline in inv.rs.
                // Panics if self is zero (no inverse exists).
                //
                // inv.rs operates on F2mNaive, so we convert via to_uint/from_uint.
                // The canonical representation is identical, so this is zero-cost.
                let naive = crate::naive::F2mNaive::<$L>::from_uint(self.c, poly);
                let inv_naive = crate::inv::ext_euclid_inv(&naive, poly);
                Self::wrap(inv_naive.to_uint())
            }

            fn div(&self, rhs: &Self, poly: &Uint<$L>) -> Self {
                self.mul(&rhs.inv(poly), poly)
            }
        }
    };
}

// ── Generate impls for each supported limb count ──────────────────────────────

impl_f2m_opt!(1, 2);
impl_f2m_opt!(2, 4);
impl_f2m_opt!(4, 8);
impl_f2m_opt!(8, 16);

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::F2m;

    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64)
    }

    fn poly8() -> Uint<1> {
        Uint::<1>::from(0x11bu64)
    }

    #[test]
    fn add_is_xor() {
        let p = poly4();
        let a = F2mOpt::<1>::from_u64(0b0101, &p);
        let b = F2mOpt::<1>::from_u64(0b0011, &p);
        assert_eq!(a.add(&b).to_uint(), Uint::<1>::from(0b0110u64));
    }

    #[test]
    fn mul_by_one_is_identity() {
        let p = poly4();
        let one = F2mOpt::<1>::one();
        for v in 0u64..16 {
            let a = F2mOpt::<1>::from_u64(v, &p);
            assert_eq!(a.mul(&one, &p), a, "a * 1 != a for a={v:#x}");
        }
    }

    #[test]
    fn aes_mul_known() {
        // AES GF(2^8): 0x53 * 0xca = 0x01.
        let p = poly8();
        let a = F2mOpt::<1>::from_u64(0x53, &p);
        let b = F2mOpt::<1>::from_u64(0xca, &p);
        let prod = a.mul(&b, &p);
        assert_eq!(prod.to_uint(), Uint::<1>::ONE, "0x53 * 0xca should be 1 in AES GF(2^8)");
    }

    #[test]
    fn gf4_mul_x2_x2() {
        // x² * x² = x⁴ ≡ x + 1 (mod x⁴+x+1).
        let p = poly4();
        let x2 = F2mOpt::<1>::from_u64(0b0100, &p);
        let prod = x2.mul(&x2, &p);
        assert_eq!(prod.to_uint(), Uint::<1>::from(0b0011u64));
    }

    #[test]
    fn square_equals_mul_self() {
        let p = poly8();
        for v in [0u64, 1, 2, 5, 0x53, 0xab, 0xff] {
            let a = F2mOpt::<1>::from_u64(v, &p);
            assert_eq!(
                a.square(&p),
                a.mul(&a, &p),
                "square != mul(self,self) for v={v:#x}"
            );
        }
    }
}
