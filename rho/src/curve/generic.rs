//! Generic Weierstrass curve over GF(p) with a small prime-order subgroup, for runnable ECDLP
//! demonstrations.
//!
//! The curve is `y² = x³ − 3x + 1` over the 63-bit prime
//! `p = 2^63 − 25 = 9_223_372_036_854_775_783`. The stored generator generates a subgroup of
//! prime order [`N`] = 561_103, a 20-bit number (about 2¹⁹). That makes this a field of 63 bits
//! carrying a subgroup of about 2¹⁹ elements: a rho collision there takes about a thousand steps
//! (from `√n`, not measured), so a solver run finishes in a fraction of a second. It is a
//! runnable fixture, not a large-group performance baseline.
//!
//! # Parameters
//!
//! - **p** = `9_223_372_036_854_775_783`  (63-bit prime)
//! - **a** = `p − 3`  (i.e. −3 mod p, same as NIST-style curves)
//! - **b** = `1`
//! - **G** = [`GX`], [`GY`] = `(3_569_091_914_278_015_591, 7_073_522_384_094_063_206)`
//! - **n** = [`N`] = `561_103`  (prime order of G)
//!
//! # Derivation of the generator
//!
//! The curve was found with a point `P₀ = (3, 821_487_384_573_098_969)` ([`ORIGINAL_GX`],
//! [`ORIGINAL_GY`]) whose order is large. G is the multiple `c·P₀` for
//! `c = 8_218_965_179_160` ([`PROJECTION_MULTIPLIER`]), chosen so that G has the largest prime
//! order that divides the order of `P₀`.
//!
//! PARI/GP found `#E = 9_223_372_037_844_426_960` and `ord P₀ = 4_611_686_018_922_213_480`; that
//! 561_103 is the largest prime factor of the latter is also a PARI/GP result. None of these full
//! orders or the largest-factor claim is certified in this tree. What the tests prove is
//! independent of them: G is on the curve, G ≠ ∞, `n·G = ∞` with `n` prime (so G has exact order
//! `n`), and `c·P₀ = G`.
//!
//! The group order is therefore the order of the subgroup, not of the curve. A base point or
//! target passed to a solver must lie in the subgroup `⟨G⟩`; the solvers' input check rejects a
//! point such as `P₀` (outside the subgroup) rather than projecting it.
//!
//! # Cross-check note
//!
//! `k256` / `p256` operate on different primes, so we cannot cross-check
//! this curve against those crates directly.  The KATs instead cross-check
//! the group law against an independent Python/affine-coordinate reference
//! implementation (recorded values for multiples of `P₀` are hard-coded in the test module).

use crypto_bigint::Uint;
use crate::curve::Curve;

/// Field characteristic `p = 2^63 − 25`.
pub const P: u64 = 9_223_372_036_854_775_783;
/// Prime order of the stored generator: `561_103`.
pub const N: u64 = 561_103;
/// Generator x-coordinate.
pub const GX: u64 = 3_569_091_914_278_015_591;
/// Generator y-coordinate.
pub const GY: u64 = 7_073_522_384_094_063_206;
/// x-coordinate of the original point `P₀`, of large order; the generator is a multiple of it.
pub const ORIGINAL_GX: u64 = 3;
/// y-coordinate of the original point `P₀`.
pub const ORIGINAL_GY: u64 = 821_487_384_573_098_969;
/// Scalar `c` with `c·P₀ =` the generator.
pub const PROJECTION_MULTIPLIER: u64 = 8_218_965_179_160;

/// Return the generic 63-bit Weierstrass test curve with its order-[`N`] generator.
///
/// `y² = x³ − 3x + 1 mod p`,  `p = 2^63 − 25`; `n` is the prime order [`N`] of the generator,
/// not the order of the curve.
pub fn generic_curve() -> Curve {
    Curve {
        p: Uint::<4>::from(P),
        // a = p - 3 ≡ -3 (mod p)
        a: Uint::<4>::from(P - 3),
        b: Uint::<4>::ONE,
        n: Uint::<4>::from(N),
        gx: Uint::<4>::from(GX),
        gy: Uint::<4>::from(GY),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::{Fp, FpMonty4 as FpMonty};

    use crate::curve::{AffinePoint, JacobianPoint};

    // Reference points computed by independent Python affine arithmetic:
    //   y² = x³ − 3x + 1  mod  9_223_372_036_854_775_783
    //
    //   2P₀ = (3_398_084_434_630_706_869, 1_269_927_104_498_887_686)
    //   3P₀ = (6_439_912_711_603_677_474, 4_727_936_577_064_732_748)
    //   4P₀ = (2_537_216_284_137_677_713, 4_719_361_313_946_664_957)
    //   5P₀ = (5_154_529_326_275_311_306, 8_346_922_522_928_899_966)

    fn ref_points() -> [(u64, u64); 5] {
        [
            (3, 821_487_384_573_098_969),                      // 1P₀
            (3_398_084_434_630_706_869, 1_269_927_104_498_887_686), // 2P₀
            (6_439_912_711_603_677_474, 4_727_936_577_064_732_748), // 3P₀
            (2_537_216_284_137_677_713, 4_719_361_313_946_664_957), // 4P₀
            (5_154_529_326_275_311_306, 8_346_922_522_928_899_966), // 5P₀
        ]
    }

    /// The original point `P₀`, built from its own constants (not from `generator()`).
    fn original_point(c: &Curve) -> AffinePoint<FpMonty> {
        AffinePoint::new(
            FpMonty::from_uint(Uint::<4>::from(ORIGINAL_GX), &c.p),
            FpMonty::from_uint(Uint::<4>::from(ORIGINAL_GY), &c.p),
        )
    }

    #[test]
    fn original_point_on_curve() {
        let c = generic_curve();
        assert!(c.is_on_curve(&original_point(&c)), "P₀ not on curve");
    }

    #[test]
    fn generator_on_curve() {
        let c = generic_curve();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "generator not on curve");
    }

    #[test]
    fn scalar_mul_small_matches_reference() {
        let c = generic_curve();
        let g = original_point(&c);
        let refs = ref_points();

        for (k, &(ref_x, ref_y)) in refs.iter().enumerate() {
            let k_val = (k + 1) as u64;
            let result = c.scalar_mul(&g, &Uint::<4>::from(k_val));
            match &result {
                AffinePoint::Infinity => panic!("{k_val}·P₀ should not be infinity"),
                AffinePoint::Finite { x, y } => {
                    assert_eq!(
                        x.to_uint(),
                        Uint::<4>::from(ref_x),
                        "{k_val}·P₀ x-coordinate mismatch"
                    );
                    assert_eq!(
                        y.to_uint(),
                        Uint::<4>::from(ref_y),
                        "{k_val}·P₀ y-coordinate mismatch"
                    );
                }
            }
            assert!(c.is_on_curve(&result), "{k_val}·P₀ not on curve");
        }
    }

    #[test]
    fn add_then_scalar_mul_consistent() {
        // 2P₀ via double must equal 2P₀ via scalar_mul.
        let c = generic_curve();
        let p = &c.p;
        let g = original_point(&c);
        let gj = JacobianPoint::from_affine(&g, p);

        let two_g_double = c.double_jacobian(&gj).to_affine(p);
        let two_g_scalar = c.scalar_mul(&g, &Uint::<4>::from(2u64));
        assert_eq!(two_g_double, two_g_scalar, "double ≠ 2·P₀ via scalar_mul");
    }

    #[test]
    fn negation_consistency() {
        // P₀ + (−P₀) = ∞
        let c = generic_curve();
        let p = &c.p;
        let g = original_point(&c);
        let neg_g = c.negate(&g);
        let gj = JacobianPoint::from_affine(&g, p);
        let neg_gj = JacobianPoint::from_affine(&neg_g, p);
        let sum = c.add_jacobian(&gj, &neg_gj).to_affine(p);
        assert!(sum.is_infinity(), "P₀ + (−P₀) should be ∞");
    }
}
