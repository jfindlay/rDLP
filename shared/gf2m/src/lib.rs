//! GF(2^m) polynomial-basis arithmetic, generic over the limb count.
//!
//! This crate exposes the [`F2m`] trait and one concrete implementation:
//! [`F2mNaive`] (schoolbook carryless-multiply + modular reduction).
//!
//! # Design: const-generic-on-trait approach
//!
//! The trait is parameterised as `F2m<const L: usize>` where `L` is the
//! number of 64-bit limbs in the underlying `Uint<L>`.  This mirrors `Fp<L>`
//! exactly and avoids the nightly-only `generic_const_exprs` feature that
//! would be required by an explicit degree const (see `shared/field/src/lib.rs`
//! lines 7–18 for the same design note).  The degree `m` is a *runtime*
//! property recovered from the irreducible as `m = poly.bit_len() − 1`.
//!
//! # Storage: `Uint<L>` as a polynomial coefficient bit-vector
//!
//! **CRITICAL:** `Uint<L>` is used here as a *polynomial coefficient
//! bit-vector*, NOT as an integer.  Bit `i` of the stored value represents the
//! coefficient of `x^i` in the polynomial.  Only XOR, shift, and bit
//! operations are meaningful on these values.  The integer-arithmetic methods
//! on `Uint<L>` (`wrapping_add`, `rem`, `mul_wide`, etc.) are meaningless on
//! GF(2) coefficient vectors and MUST NOT be used for field arithmetic.  In
//! particular, `mul_wide` performs integer multiplication with carry — the
//! carryless multiply in `mul`/`square` uses the comb algorithm (shift-and-XOR)
//! instead.
//!
//! # Irreducible polynomial
//!
//! Every operation that depends on the field's identity takes `poly: &Uint<L>`
//! as a parameter, mirroring `Fp`'s per-call `p: &Uint<L>`.  The irreducible
//! is the bit-vector of the degree-`m` reduction polynomial: bit `i` set ⟺
//! coefficient of `x^i` is 1.  Examples:
//! - GF(2^4) with `x⁴+x+1`:           `poly = 0b1_0011`  (= 0x13)
//! - GF(2^8) AES with `x⁸+x⁴+x³+x+1`: `poly = 0b1_0001_1011` (= 0x11b)
//!
//! Operations mixing elements reduced under *different* irreducibles are
//! meaningless — the caller is responsible for consistency (analogous to the
//! C-MovBridge modulus-consistency guard in `shared/field`).
//!
//! # Characteristic-2 invariants
//!
//! - **`sub == add`**: subtraction in GF(2^m) is XOR, identical to addition.
//!   Do NOT port `Fp`'s `sub` (which computes `p − self`).
//! - **`neg == identity`**: negation in characteristic 2 is the identity map
//!   (`−a = a`).  Do NOT port `Fp`'s `neg` (which computes `p − self`).
//!
//! # Toy field sizes
//!
//! The KATs use GF(2^4) and GF(2^8).  The algorithms are not toy — they are
//! the same algorithms used for crypto-scale `m` (163, 233, 571) — but the
//! parameters are small for auditability (principle-4 boundary).

pub mod convert;
pub mod inv;
pub mod naive;
pub mod normal;
pub mod opt;
pub mod poly;
pub mod subfield;

pub use convert::{find_normal_element, frobenius_orbit, is_normal_element, normal_to_poly, poly_to_normal};
pub use inv::{ext_euclid_inv, itoh_tsujii_inv};
pub use naive::F2mNaive;
pub use normal::F2mNormal;
pub use opt::F2mOpt;
pub use poly::Poly;
pub use subfield::{
    embed, frobenius_subfield_orbit, is_in_subfield, relative_norm, relative_trace, restrict,
    subfield_basis,
};

use crypto_bigint::Uint;

/// `usize`-indexed bit access and shifts over a `Uint`.
///
/// `crypto_bigint::Uint` indexes and measures bits in `u32`, while the degree arithmetic in this
/// crate (`m = poly.bit_len() − 1`) works in `usize`. The shifts are the panic-free
/// `unbounded_*` variants: a shift of `BITS` or more yields zero rather than panicking.
pub trait BitLen: Sized {
    /// Number of bits needed to represent the value (zero for zero).
    fn bit_len(&self) -> usize;

    /// Bit `i` (little-endian from bit 0); `false` when `i` is out of range.
    fn bit_at(&self, i: usize) -> bool;

    /// `self << n`, zero if `n` exceeds the precision. Variable time in `n` only.
    #[must_use]
    fn shl_bits(&self, n: usize) -> Self;

    /// `self >> n`, zero if `n` exceeds the precision. Variable time in `n` only.
    #[must_use]
    fn shr_bits(&self, n: usize) -> Self;
}

impl<const L: usize> BitLen for Uint<L> {
    #[inline]
    fn bit_len(&self) -> usize {
        self.bits() as usize
    }

    #[inline]
    fn bit_at(&self, i: usize) -> bool {
        u32::try_from(i).is_ok_and(|i| self.bit_vartime(i))
    }

    #[inline]
    fn shl_bits(&self, n: usize) -> Self {
        u32::try_from(n).map_or(Self::ZERO, |n| self.unbounded_shl_vartime(n))
    }

    #[inline]
    fn shr_bits(&self, n: usize) -> Self {
        u32::try_from(n).map_or(Self::ZERO, |n| self.unbounded_shr_vartime(n))
    }
}

/// GF(2^m) field arithmetic, generic over the limb count `L`.
///
/// All values are implicitly reduced mod the irreducible polynomial `poly`.
/// The irreducible is passed per-call (not stored on the struct), mirroring
/// `Fp`'s per-call `p: &Uint<L>`.
///
/// # Const-generic design
///
/// Parameterised as `F2m<const L: usize>` where `L` = ⌈m/64⌉ (the number of
/// 64-bit limbs).  The degree `m` is recovered at runtime as `poly.bit_len() − 1`.
///
/// # Storage contract
///
/// The `Uint<L>` backing is a *polynomial coefficient bit-vector*, not an
/// integer.  Only XOR, shift, and bit operations are meaningful.  See the
/// crate-level documentation for the full warning.
///
/// # Characteristic-2 invariants
///
/// - `sub(a, b) == add(a, b)` — subtraction is XOR in char 2.
/// - `neg(a) == a` — negation is the identity in char 2.
///
/// Both invariants are KAT-checked in `shared/gf2m/tests/gf2m_kat.rs`.
pub trait F2m<const L: usize>: Clone + PartialEq + Eq + std::fmt::Debug + Send + Sync + 'static {
    // ── Constructors / canonical form ─────────────────────────────────────────

    /// Additive identity: the zero polynomial.
    fn zero() -> Self;

    /// Multiplicative identity: the polynomial `1`.
    fn one() -> Self;

    /// Construct from a small `u64` value, reducing mod `poly` if needed.
    fn from_u64(v: u64, poly: &Uint<L>) -> Self;

    /// Construct from an arbitrary `Uint<L>` coefficient bit-vector, reducing
    /// mod `poly` if the degree of `v` is ≥ degree of `poly`.
    fn from_uint(v: Uint<L>, poly: &Uint<L>) -> Self;

    /// Return the canonical coefficient bit-vector in `[0, 2^m)`.
    fn to_uint(&self) -> Uint<L>;

    // ── Char-2 arithmetic (implemented in `F2mNaive`) ────────────────────────

    /// Addition: `a + b` in GF(2^m) = XOR of coefficient bit-vectors.
    ///
    /// No `poly` needed — XOR is always in-field.
    #[must_use]
    fn add(&self, rhs: &Self) -> Self;

    /// Subtraction: `a − b` in GF(2^m).
    ///
    /// **In characteristic 2, subtraction equals addition (XOR).**
    /// This method exists for API symmetry with `Fp`; it MUST equal `add`.
    /// Do NOT implement as `p − self` (that is the `Fp` convention, wrong here).
    ///
    /// Invariant: `sub(a, b) == add(a, b)` for all `a`, `b`.
    #[must_use]
    fn sub(&self, rhs: &Self) -> Self;

    /// Negation: `−a` in GF(2^m).
    ///
    /// **In characteristic 2, negation is the identity (`−a = a`).**
    /// This method exists for API symmetry with `Fp`; it MUST return `self`.
    /// Do NOT implement as `p − self` (that is the `Fp` convention, wrong here).
    ///
    /// Invariant: `neg(a) == a` for all `a`.
    #[must_use]
    fn neg(&self) -> Self;

    /// Multiplication: carryless comb multiply followed by reduction mod `poly`.
    ///
    /// The intermediate product has degree < 2m and is reduced by XOR-shifting
    /// `poly` to align with each high-degree term.  See `naive.rs` for the
    /// auditable baseline.
    #[must_use]
    fn mul(&self, rhs: &Self, poly: &Uint<L>) -> Self;

    /// Squaring: Frobenius map `a → a²` via bit-spread then reduction mod `poly`.
    ///
    /// Equivalent to `mul(a, a, poly)` but implemented via the Frobenius
    /// bit-spread (insert a zero bit between each coefficient bit, then reduce).
    /// This is the characteristic-2 Frobenius endomorphism.
    #[must_use]
    fn square(&self, poly: &Uint<L>) -> Self;

    /// Frobenius endomorphism: `a → a²` (first-class name for binary-curve and GHS consumers).
    ///
    /// Delegates to `square`.  Declared as a distinct method because binary-curve and GHS
    /// code reaches for it by this name and may iterate it as the Frobenius map.
    #[must_use]
    fn frobenius(&self, poly: &Uint<L>) -> Self;

    /// Exponentiation: `a^exp` via square-and-multiply.
    ///
    /// `exp` is treated as a non-negative integer exponent (not a field element).
    #[must_use]
    fn pow(&self, exp: &Uint<L>, poly: &Uint<L>) -> Self;

    // ── Curve-facing over-specification (declared-and-stubbed; filled in binary-curve layer) ──

    /// Absolute trace: `Tr(a) = Σ_{i=0}^{m-1} a^(2^i)` in GF(2^m).
    ///
    /// Returns 0 or 1 (as a field element).  Used by binary-curve point operations.
    ///
    /// # Panics
    ///
    /// Not yet implemented — filled in the binary-curve layer.
    #[must_use]
    fn trace(&self, poly: &Uint<L>) -> Self;

    /// Half-trace / solve quadratic: find `x` such that `x² + x = c` in GF(2^m).
    ///
    /// Solves the Artin–Schreier equation `x² + x = c`.  Used by binary-curve point operations
    /// (point decompression over GF(2^m)).
    ///
    /// # Panics
    ///
    /// Not yet implemented — filled in the binary-curve layer.
    #[must_use]
    fn solve_quadratic(c: &Self, poly: &Uint<L>) -> Self;

    // ── Deferred (declared-and-stubbed) ──────────────────────────────────────

    /// Multiplicative inverse: `a⁻¹` such that `a · a⁻¹ = 1`.
    ///
    /// # Panics
    ///
    /// Panics if `self` is zero (no inverse exists).
    /// Not yet implemented — deferred to a future extension.
    #[must_use]
    fn inv(&self, poly: &Uint<L>) -> Self;

    /// Division: `a / b = a · b⁻¹`.
    ///
    /// # Panics
    ///
    /// Panics if `rhs` is zero.
    /// Not yet implemented — deferred to a future extension.
    #[must_use]
    fn div(&self, rhs: &Self, poly: &Uint<L>) -> Self;

    // ── Default helpers ───────────────────────────────────────────────────────

    /// Return `true` if this element is the zero polynomial.
    fn is_zero(&self) -> bool {
        self.to_uint() == Uint::<L>::ZERO
    }

    /// Return `true` if this element is the polynomial `1`.
    fn is_one(&self) -> bool {
        self.to_uint() == Uint::<L>::ONE
    }
}
