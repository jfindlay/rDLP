//! Prime-field arithmetic over GF(p), generic over the limb count.
//!
//! This crate exposes the [`Fp`] trait and two concrete implementations:
//! [`FpNaive`] (schoolbook arithmetic) and [`FpMonty`] (Montgomery form via
//! ``crypto-bigint``'s ``FixedMontyForm``).
//!
//! # Design: const-generic-on-trait approach
//!
//! The trait is parameterised as ``Fp<const L: usize>`` where ``L`` is the
//! number of 64-bit limbs in the underlying ``Uint<L>``.  This is stable Rust
//! (edition 2024) and avoids the nightly-only ``generic_const_exprs`` feature
//! that would be required by the alternative approach of an associated constant
//! ``const LIMBS: usize`` with a ``type Uint = Uint<{Self::LIMBS}>`` associated
//! type.
//!
//! Callers that want the 256-bit (4-limb) field used throughout ``rho`` can
//! write ``F: Fp<4>`` or use the type aliases ``FpNaive4`` / ``FpMonty4``
//! exported from this crate.
//!
//! # Square roots and the Legendre symbol
//!
//! The [`Fp`] trait provides default implementations of [`Fp::legendre`] (via
//! Euler's criterion) and [`Fp::sqrt`] (via Tonelli–Shanks).  Both assume ``p``
//! is an odd prime; behaviour is unspecified for composite moduli.
//!
//! ## Tonelli–Shanks algorithm
//!
//! For `p ≡ 3 (mod 4)` the shortcut `a^((p+1)/4)` is used.  For `p ≡ 1 (mod 4)`
//! the full Tonelli–Shanks loop is used.  A quadratic non-residue (QNR) is found
//! by trial starting at `n = 2, 3, 5, …`; this converges within ~10 trials for
//! all project primes.

pub mod monty;
pub mod naive;

pub use monty::FpMonty;
pub use naive::FpNaive;

/// Type alias: schoolbook 256-bit field element (4 × 64-bit limbs).
pub type FpNaive4 = FpNaive<4>;

/// Type alias: Montgomery-form 256-bit field element (4 × 64-bit limbs).
pub type FpMonty4 = FpMonty<4>;

use crypto_bigint::Uint;

/// Prime-field arithmetic over GF(p), generic over the limb count ``L``.
///
/// All values are implicitly reduced mod p.  Implementations are allowed to use
/// internal representations (e.g., Montgomery form) as long as [`to_uint`]
/// returns the canonical residue in ``[0, p)``.
///
/// # Const-generic design
///
/// The trait is parameterised as ``Fp<const L: usize>`` where ``L`` is the
/// number of 64-bit limbs in ``Uint<L>``.  Choose ``L`` so that
/// ``L * 64 >= bit-width(p)``.  For example, ``L = 4`` covers 256-bit primes
/// (secp256k1, P-256, etc.).
///
/// # Square roots and the Legendre symbol
///
/// Default implementations of [`legendre`] and [`sqrt`] are provided.  Both
/// assume ``p`` is an odd prime; behaviour is unspecified for composite moduli.
/// Implementations may override these defaults for performance, but the
/// mathematical contract must be preserved.
pub trait Fp<const L: usize>: Clone + PartialEq + Eq + std::fmt::Debug + Send + Sync + 'static {
    /// Additive identity: 0 mod p.
    fn zero(p: &Uint<L>) -> Self;

    /// Multiplicative identity: 1 mod p.
    fn one(p: &Uint<L>) -> Self;

    /// Construct from a small ``u64`` value, reducing mod p.
    fn from_u64(v: u64, p: &Uint<L>) -> Self;

    /// Construct from an arbitrary ``Uint<L>``, reducing mod p.
    fn from_uint(v: Uint<L>, p: &Uint<L>) -> Self;

    /// Return the canonical residue in ``[0, p)``.
    fn to_uint(&self) -> Uint<L>;

    /// Modular addition: ``self + rhs mod p``.
    fn add(&self, rhs: &Self, p: &Uint<L>) -> Self;

    /// Modular subtraction: ``self - rhs mod p``.
    fn sub(&self, rhs: &Self, p: &Uint<L>) -> Self;

    /// Modular negation: ``-self mod p``.
    fn neg(&self, p: &Uint<L>) -> Self;

    /// Modular multiplication: ``self * rhs mod p``.
    fn mul(&self, rhs: &Self, p: &Uint<L>) -> Self;

    /// Modular squaring: ``self^2 mod p``.
    ///
    /// May be faster than ``mul(self, self)`` for implementations that exploit
    /// the squaring structure (e.g., Karatsuba squaring).
    fn square(&self, p: &Uint<L>) -> Self;

    /// Modular exponentiation: ``self^exp mod p``.
    fn pow(&self, exp: &Uint<L>, p: &Uint<L>) -> Self;

    /// Modular inverse via Fermat's little theorem: ``self^(p-2) mod p``.
    ///
    /// # Panics
    ///
    /// Panics if ``self`` is zero (no inverse exists).
    fn inv(&self, p: &Uint<L>) -> Self;

    /// Return ``true`` if this element is zero.
    ///
    /// Default implementation compares ``to_uint()`` to ``Uint::<L>::ZERO``.
    /// Implementations may override for efficiency.
    fn is_zero(&self, _p: &Uint<L>) -> bool {
        self.to_uint() == Uint::<L>::ZERO
    }

    /// Return ``true`` if this element is the multiplicative identity (1 mod p).
    ///
    /// Default implementation compares ``to_uint()`` to ``Uint::<L>::ONE``.
    fn is_one(&self, _p: &Uint<L>) -> bool {
        self.to_uint() == Uint::<L>::ONE
    }

    /// Modular doubling: ``2 * self mod p``.
    ///
    /// Default implementation delegates to ``add(self, self, p)``.
    /// Implementations may override with a faster shift-based path.
    fn double(&self, p: &Uint<L>) -> Self {
        self.add(self, p)
    }

    /// Legendre symbol ``(self / p)`` via Euler's criterion.
    ///
    /// Returns ``0`` if ``self ≡ 0 (mod p)``, ``1`` if ``self`` is a quadratic
    /// residue mod ``p``, or ``-1`` if ``self`` is a quadratic non-residue.
    ///
    /// Assumes ``p`` is an odd prime.  Behaviour is unspecified for composite
    /// moduli.
    fn legendre(&self, p: &Uint<L>) -> i8 {
        if self.is_zero(p) {
            return 0;
        }
        // Euler's criterion: a^((p-1)/2) mod p is 1 (QR) or p-1 (QNR).
        let exp = p.wrapping_sub(&Uint::<L>::ONE) >> 1;
        let r = self.pow(&exp, p).to_uint();
        if r == Uint::<L>::ONE {
            1
        } else {
            // r == p - 1 for a QNR (since p is prime and r^2 ≡ 1 mod p).
            -1
        }
    }

    /// Modular square root via Tonelli–Shanks.
    ///
    /// Returns ``Some(r)`` where ``r^2 ≡ self (mod p)`` if ``self`` is a
    /// quadratic residue mod ``p``, or ``None`` if it is a non-residue.
    /// Returns ``Some(0)`` for ``self ≡ 0``.
    ///
    /// Uses the ``p ≡ 3 (mod 4)`` shortcut ``a^((p+1)/4)`` when applicable;
    /// falls back to the full Tonelli–Shanks loop for ``p ≡ 1 (mod 4)``.
    ///
    /// Assumes ``p`` is an odd prime.  Behaviour is unspecified for composite
    /// moduli.
    fn sqrt(&self, p: &Uint<L>) -> Option<Self> {
        if self.is_zero(p) {
            return Some(Self::zero(p));
        }
        if self.legendre(p) != 1 {
            return None;
        }

        // p mod 4: check the two low bits of p.
        let p_mod4 = p.as_words()[0] & 3;

        if p_mod4 == 3 {
            // Shortcut: p ≡ 3 (mod 4) → sqrt = a^((p+1)/4).
            let exp = p.wrapping_add(&Uint::<L>::ONE) >> 2;
            return Some(self.pow(&exp, p));
        }

        // General Tonelli–Shanks for p ≡ 1 (mod 4).
        //
        // Factor p-1 = Q * 2^S with Q odd.
        let pm1 = p.wrapping_sub(&Uint::<L>::ONE);
        let s = pm1.trailing_zeros() as usize;
        let q = pm1 >> s; // Q = (p-1) / 2^S, odd

        // Find a quadratic non-residue z by trial: try 2, 3, 5, ...
        // For all project primes this converges within ~10 trials.
        let z = {
            let mut candidate = 2u64;
            loop {
                let c = Self::from_u64(candidate, p);
                if c.legendre(p) == -1 {
                    break c;
                }
                candidate += 1;
            }
        };

        // Initialise the Tonelli–Shanks state.
        let mut m: usize = s;
        let mut c = z.pow(&q, p);
        let mut t = self.pow(&q, p);
        // r = a^((Q+1)/2): Q is odd so (Q+1)/2 is exact.
        let mut r = self.pow(&q.wrapping_add(&Uint::<L>::ONE).shr_vartime(1), p);

        let one = Self::one(p);

        loop {
            if t == one {
                return Some(r);
            }

            // Find the least i (1 ≤ i < m) such that t^(2^i) ≡ 1.
            let mut i: usize = 1;
            let mut tmp = t.square(p);
            while tmp != one {
                tmp = tmp.square(p);
                i += 1;
            }

            // b = c^(2^(m-i-1))
            let mut b = c.clone();
            for _ in 0..(m - i - 1) {
                b = b.square(p);
            }

            m = i;
            c = b.square(p);
            t = t.mul(&c, p);
            r = r.mul(&b, p);
        }
    }
}
