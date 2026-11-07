//! Montgomery-form field arithmetic over GF(p), generic over the limb count.
//!
//! ``FpMonty<L>`` wraps ``crypto_bigint::modular::FixedMontyForm<L>``
//! to get the Montgomery multiplication speedup.  The internal representation
//! is ``aR mod p`` where ``R = 2^(64·L)``.  Conversion to/from canonical form
//! is transparent via [`Fp::from_uint`] / [`Fp::to_uint`].
//!
//! ``FixedMontyForm<L>`` and ``FixedMontyParams<L>`` are const-generic over the
//! limb count, so the generic-over-``L`` reshape is straightforward.

use crypto_bigint::{
    modular::{FixedMontyForm, FixedMontyParams},
    Odd, Uint,
};

use super::Fp;

/// Montgomery-form field element backed by ``crypto-bigint``'s ``FixedMontyForm``.
///
/// The ``FixedMontyForm`` carries its ``FixedMontyParams`` internally, which encodes
/// the modulus and the Montgomery constant.  We cache the params as a ``Copy``
/// field so we can reconstruct elements after operations that consume ``self``.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FpMonty<const L: usize> {
    /// Inner ``FixedMontyForm`` (Montgomery form ``v * R mod p``).
    inner: FixedMontyForm<L>,
    /// Cached params — re-used on every operation.
    params: FixedMontyParams<L>,
}

impl<const L: usize> FpMonty<L> {
    /// Build ``FixedMontyParams`` from a modulus ``p``.
    ///
    /// Panics if ``p`` is even: Montgomery reduction requires an odd modulus.
    #[inline]
    fn params_from(p: &Uint<L>) -> FixedMontyParams<L> {
        FixedMontyParams::<L>::new(Odd::new(*p).expect("Montgomery modulus must be odd"))
    }

    /// Wrap a ``FixedMontyForm`` together with its params.
    #[inline]
    fn wrap(inner: FixedMontyForm<L>, params: FixedMontyParams<L>) -> Self {
        Self { inner, params }
    }
}

impl<const L: usize> Fp<L> for FpMonty<L> {
    fn zero(p: &Uint<L>) -> Self {
        let params = Self::params_from(p);
        let inner = FixedMontyForm::zero(&params);
        Self::wrap(inner, params)
    }

    fn one(p: &Uint<L>) -> Self {
        let params = Self::params_from(p);
        let inner = FixedMontyForm::one(&params);
        Self::wrap(inner, params)
    }

    fn from_u64(v: u64, p: &Uint<L>) -> Self {
        let params = Self::params_from(p);
        let inner = FixedMontyForm::new(&Uint::<L>::from(v), &params);
        Self::wrap(inner, params)
    }

    fn from_uint(v: Uint<L>, p: &Uint<L>) -> Self {
        let params = Self::params_from(p);
        let inner = FixedMontyForm::new(&v, &params);
        Self::wrap(inner, params)
    }

    fn to_uint(&self) -> Uint<L> {
        self.inner.retrieve()
    }

    fn add(&self, rhs: &Self, _p: &Uint<L>) -> Self {
        let result = self.inner.add(&rhs.inner);
        Self::wrap(result, self.params)
    }

    fn sub(&self, rhs: &Self, _p: &Uint<L>) -> Self {
        let result = self.inner.sub(&rhs.inner);
        Self::wrap(result, self.params)
    }

    fn neg(&self, _p: &Uint<L>) -> Self {
        let result = self.inner.neg();
        Self::wrap(result, self.params)
    }

    fn mul(&self, rhs: &Self, _p: &Uint<L>) -> Self {
        let result = self.inner.mul(&rhs.inner);
        Self::wrap(result, self.params)
    }

    fn square(&self, _p: &Uint<L>) -> Self {
        let result = self.inner.square();
        Self::wrap(result, self.params)
    }

    fn pow(&self, exp: &Uint<L>, _p: &Uint<L>) -> Self {
        let result = self.inner.pow(exp);
        Self::wrap(result, self.params)
    }

    fn inv(&self, p: &Uint<L>) -> Self {
        assert!(!self.is_zero(p), "attempted inversion of zero");
        let exp = p.wrapping_sub(&Uint::<L>::from(2u64));
        self.pow(&exp, p)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn p13() -> Uint<4> {
        Uint::<4>::from(13u64)
    }

    fn secp_p() -> Uint<4> {
        Uint::<4>::from_be_hex(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
        )
    }

    #[test]
    fn add_small() {
        let p = p13();
        let a = FpMonty::<4>::from_u64(7, &p);
        let b = FpMonty::<4>::from_u64(9, &p);
        assert_eq!(a.add(&b, &p).to_uint(), Uint::<4>::from(3u64));
    }

    #[test]
    fn sub_wrap() {
        let p = p13();
        let a = FpMonty::<4>::from_u64(3, &p);
        let b = FpMonty::<4>::from_u64(7, &p);
        assert_eq!(a.sub(&b, &p).to_uint(), Uint::<4>::from(9u64));
    }

    #[test]
    fn mul_small() {
        let p = p13();
        let a = FpMonty::<4>::from_u64(5, &p);
        let b = FpMonty::<4>::from_u64(6, &p);
        assert_eq!(a.mul(&b, &p).to_uint(), Uint::<4>::from(4u64));
    }

    #[test]
    fn inv_small() {
        let p = p13();
        let a = FpMonty::<4>::from_u64(5, &p);
        let ai = a.inv(&p);
        let one = a.mul(&ai, &p);
        assert_eq!(one, FpMonty::<4>::one(&p));
    }

    #[test]
    fn neg_zero_is_zero() {
        let p = p13();
        let z = FpMonty::<4>::zero(&p);
        assert_eq!(z.neg(&p), FpMonty::<4>::zero(&p));
    }

    #[test]
    fn add_neg_is_zero() {
        let p = secp_p();
        let a = FpMonty::<4>::from_u64(0xDEAD_BEEF, &p);
        let neg_a = a.neg(&p);
        let sum = a.add(&neg_a, &p);
        assert_eq!(sum, FpMonty::<4>::zero(&p));
    }

    #[test]
    fn mul_secp_prime() {
        let p = secp_p();
        let two = FpMonty::<4>::from_u64(2, &p);
        let half = (p.wrapping_add(&Uint::<4>::ONE)) >> 1;
        let half_fp = FpMonty::<4>::from_uint(half, &p);
        let one = two.mul(&half_fp, &p);
        assert_eq!(one, FpMonty::<4>::one(&p));
    }

    #[test]
    fn is_zero_and_is_one() {
        let p = p13();
        let z = FpMonty::<4>::zero(&p);
        let o = FpMonty::<4>::one(&p);
        assert!(z.is_zero(&p));
        assert!(!o.is_zero(&p));
        assert!(o.is_one(&p));
        assert!(!z.is_one(&p));
    }

    #[test]
    fn double_matches_add_self() {
        let p = p13();
        let a = FpMonty::<4>::from_u64(5, &p);
        assert_eq!(a.double(&p), a.add(&a, &p));
    }

    #[test]
    fn matches_naive() {
        // Cross-check FpMonty against FpNaive for a set of operations on the secp256k1 prime.
        use crate::naive::FpNaive;

        let p = secp_p();
        let vals: &[u64] = &[0, 1, 2, 12345, 0xDEAD_BEEF, 0xFFFF_FFFF_FFFF];

        for &a_u in vals {
            for &b_u in vals {
                let an = FpNaive::<4>::from_u64(a_u, &p);
                let bn = FpNaive::<4>::from_u64(b_u, &p);
                let am = FpMonty::<4>::from_u64(a_u, &p);
                let bm = FpMonty::<4>::from_u64(b_u, &p);

                assert_eq!(
                    an.add(&bn, &p).to_uint(),
                    am.add(&bm, &p).to_uint(),
                    "add mismatch a={a_u} b={b_u}"
                );
                assert_eq!(
                    an.sub(&bn, &p).to_uint(),
                    am.sub(&bm, &p).to_uint(),
                    "sub mismatch a={a_u} b={b_u}"
                );
                assert_eq!(
                    an.mul(&bn, &p).to_uint(),
                    am.mul(&bm, &p).to_uint(),
                    "mul mismatch a={a_u} b={b_u}"
                );
                if b_u != 0 {
                    assert_eq!(
                        an.mul(&bn.inv(&p), &p).to_uint(),
                        am.mul(&bm.inv(&p), &p).to_uint(),
                        "inv-mul mismatch a={a_u} b={b_u}"
                    );
                }
            }
        }
    }

    /// Verify that FpMonty<1> works for a small 1-limb prime.
    #[test]
    fn monty_1limb_basic() {
        let p = Uint::<1>::from(13u64);
        let a = FpMonty::<1>::from_u64(7, &p);
        let b = FpMonty::<1>::from_u64(9, &p);
        // 7 + 9 = 16 ≡ 3 (mod 13)
        assert_eq!(a.add(&b, &p).to_uint(), Uint::<1>::from(3u64));
        // 5 * 6 = 30 ≡ 4 (mod 13)
        let c = FpMonty::<1>::from_u64(5, &p);
        let d = FpMonty::<1>::from_u64(6, &p);
        assert_eq!(c.mul(&d, &p).to_uint(), Uint::<1>::from(4u64));
    }
}
