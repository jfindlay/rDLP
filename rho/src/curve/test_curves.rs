//! Small curves for unit tests and ECDLP KATs.
//!
//! Parameters were computed offline via point-counting (brute-force enumeration of
//! `y² = x³ + ax + b mod p` for each `x`) and verified by checking that `n·G = ∞`.
//!
//! # Curve A — `tiny_a`
//!
//! `y² = x³ − 3x + 3  mod  1_048_517`
//!
//! | Symbol | Value |
//! |--------|-------|
//! | p | `1_048_517` (20-bit prime) |
//! | a | `1_048_514` (= −3 mod p) |
//! | b | `3` |
//! | n | `1_048_051` (20-bit prime) |
//! | G | `(1, 1)` |
//!
//! # Curve B — `tiny_b`
//!
//! `y² = x³ − 3x + 16  mod  1_048_583`
//!
//! | Symbol | Value |
//! |--------|-------|
//! | p | `1_048_583` (20-bit prime) |
//! | a | `1_048_580` (= −3 mod p) |
//! | b | `16` |
//! | n | `1_048_387` (20-bit prime) |
//! | G | `(0, 4)` |
//!
//! # Composite-order toy curve — `composite_toy`
//!
//! `y² = x³ + x + 33  mod  47`
//!
//! A short-Weierstrass curve over GF(47) whose generator has **full composite**
//! order `n = 60 = 2² · 3 · 5`.  The repeated prime `2²` exercises the `e > 1`
//! lift in Pohlig–Hellman.  The three distinct prime factors exercise
//! the CRT combine.
//!
//! Parameters computed offline by brute-force point-counting and verified by:
//! - `n·G = ∞`
//! - `[n/2]·G ≠ ∞`, `[n/3]·G ≠ ∞`, `[n/5]·G ≠ ∞` (full-order checks)
//!
//! | Symbol | Value |
//! |--------|-------|
//! | p | `47` |
//! | a | `1` |
//! | b | `33` |
//! | n | `60 = 2² · 3 · 5` |
//! | G | `(10, 3)` |
//!
//! # Composite-order curve with a large prime factor — `composite_rho_toy`
//!
//! `y² = x³ + 2x + 3  mod  1009`
//!
//! Generator of full order `n = 1068 = 2² · 3 · 89`. The factor `89` exceeds the
//! subgroup-enumeration limit (64), so Pohlig–Hellman solves that component with rho.
//!
//! | Symbol | Value |
//! |--------|-------|
//! | p | `1009` |
//! | a | `2` |
//! | b | `3` |
//! | n | `1068 = 2² · 3 · 89` |
//! | G | `(10, 425)` |

use crypto_bigint::Uint;
use crate::curve::Curve;

// ── Curve A ───────────────────────────────────────────────────────────────────

/// Return the 20-bit test curve A.
///
/// `y² = x³ − 3x + 3 mod 1_048_517`.  Prime group order `n = 1_048_051`.
pub fn tiny_a() -> Curve {
    let p = Uint::<4>::from(1_048_517u64);
    Curve {
        p,
        a: Uint::<4>::from(1_048_514u64), // -3 mod p
        b: Uint::<4>::from(3u64),
        n: Uint::<4>::from(1_048_051u64),
        gx: Uint::<4>::from(1u64),
        gy: Uint::<4>::from(1u64),
    }
}

/// Group order of `tiny_a`.
pub const TINY_A_N: u64 = 1_048_051;

// ── Curve B ───────────────────────────────────────────────────────────────────

/// Return the 20-bit test curve B.
///
/// `y² = x³ − 3x + 16 mod 1_048_583`.  Prime group order `n = 1_048_387`.
pub fn tiny_b() -> Curve {
    let p = Uint::<4>::from(1_048_583u64);
    Curve {
        p,
        a: Uint::<4>::from(1_048_580u64), // -3 mod p
        b: Uint::<4>::from(16u64),
        n: Uint::<4>::from(1_048_387u64),
        gx: Uint::<4>::from(0u64),
        gy: Uint::<4>::from(4u64),
    }
}

/// Group order of `tiny_b`.
pub const TINY_B_N: u64 = 1_048_387;

// ── Tiny GLV curve ────────────────────────────────────────────────────────────

/// Return the tiny GLV-capable test curve.
///
/// `y² = x³ + 7 mod 1051`.  Prime group order `n = 1093`.
///
/// This is a secp256k1-style curve (`a = 0`, `b = 7`) over a small prime field
/// that admits the GLV endomorphism `φ(x, y) = (β·x mod p, y)` with
/// `φ(P) = λ·P` for all curve points P.
///
/// # Parameters
///
/// | Symbol | Value | Notes |
/// |--------|-------|-------|
/// | p | `1051` | prime, `p ≡ 3 (mod 4)`, `p ≡ 1 (mod 3)` |
/// | b | `7` | secp256k1-style |
/// | n | `1093` | prime group order, `n ≡ 1 (mod 3)` |
/// | G | `(3, 666)` | base point |
/// | β | `870` | cube root of unity mod p |
/// | λ | `151` | GLV scalar (`λ² + λ + 1 ≡ 0 mod n`) |
pub fn tiny_glv() -> Curve {
    Curve {
        p: Uint::<4>::from(1051u64),
        a: Uint::<4>::ZERO,
        b: Uint::<4>::from(7u64),
        n: Uint::<4>::from(1093u64),
        gx: Uint::<4>::from(3u64),
        gy: Uint::<4>::from(666u64),
    }
}

/// Group order of `tiny_glv`.
pub const TINY_GLV_N: u64 = 1093;
/// Cube root of unity mod p for `tiny_glv`.
pub const TINY_GLV_BETA: u64 = 870;
/// GLV eigenvalue for `tiny_glv`.
pub const TINY_GLV_LAMBDA: u64 = 151;

// ── Composite-order toy curve ─────────────────────────────────────────────────

/// Return the composite-order toy curve.
///
/// `y² = x³ + x + 33 mod 47`.  Composite group order `n = 60 = 2² · 3 · 5`.
/// The generator `G = (10, 3)` has full order `n` (not a subgroup generator).
///
/// Used to exercise the Pohlig–Hellman prime-power lift (`e = 2` for the
/// factor `2²`) and the CRT combine across three distinct primes.
pub fn composite_toy() -> Curve {
    Curve {
        p:  Uint::<4>::from(47u64),
        a:  Uint::<4>::from(1u64),
        b:  Uint::<4>::from(33u64),
        n:  Uint::<4>::from(60u64),
        gx: Uint::<4>::from(10u64),
        gy: Uint::<4>::from(3u64),
    }
}

/// Group order of `composite_toy`: `n = 60 = 2² · 3 · 5`.
pub const COMPOSITE_TOY_N: u64 = 60;

/// Prime-power factorization of `COMPOSITE_TOY_N`.
///
/// Sorted `(prime, exponent)` pairs.  `∏ pᵢ^{eᵢ} = 60`.
/// The `2²` entry exercises the `e > 1` lift in Pohlig–Hellman.
pub const COMPOSITE_TOY_FACTORS: &[(u64, u32)] = &[(2, 2), (3, 1), (5, 1)];

// ── Composite-order curve with a rho-sized prime factor ───────────────────────

/// Return the composite-order curve whose order has a prime factor above 64.
///
/// `y² = x³ + 2x + 3 mod 1009`.  Group order `n = 1068 = 2² · 3 · 89`, generator
/// `G = (10, 425)` of full order `n`.  The prime `89` exceeds the subgroup-enumeration limit of
/// 64, so Pohlig–Hellman solves that component with rho while the `2²` and `3` components use
/// enumeration. Every prime power is small enough for `factor_order`'s trial division.
///
/// Parameters computed offline by point-counting; verified by `n·G = ∞` and
/// `[n/2]·G, [n/3]·G, [n/89]·G ≠ ∞`.
pub fn composite_rho_toy() -> Curve {
    Curve {
        p:  Uint::<4>::from(1009u64),
        a:  Uint::<4>::from(2u64),
        b:  Uint::<4>::from(3u64),
        n:  Uint::<4>::from(1068u64),
        gx: Uint::<4>::from(10u64),
        gy: Uint::<4>::from(425u64),
    }
}

/// Group order of `composite_rho_toy`: `n = 1068 = 2² · 3 · 89`.
pub const COMPOSITE_RHO_TOY_N: u64 = 1068;

/// Prime-power factorization of `COMPOSITE_RHO_TOY_N`.
pub const COMPOSITE_RHO_TOY_FACTORS: &[(u64, u32)] = &[(2, 2), (3, 1), (89, 1)];

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::{Fp, FpMonty4 as FpMonty};

    use crate::curve::AffinePoint;

    #[test]
    fn tiny_a_generator_on_curve() {
        let c = tiny_a();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "tiny_a: generator not on curve");
    }

    #[test]
    fn tiny_a_n_times_g_is_infinity() {
        let c = tiny_a();
        let g: AffinePoint<FpMonty> = c.generator();
        let ng = c.scalar_mul(&g, &c.n);
        assert!(ng.is_infinity(), "tiny_a: n·G should be ∞");
    }

    /// Known DLP values computed by the Python reference implementation.
    ///
    /// | k | Q.x | Q.y |
    /// |---|-----|-----|
    /// | 7 | 1_026_105 | 636_225 |
    /// | 100 | 659_291 | 755_487 |
    /// | 1000 | 925_418 | 411_028 |
    #[test]
    fn tiny_a_scalar_mul_reference() {
        let c = tiny_a();
        let g: AffinePoint<FpMonty> = c.generator();
        let cases: &[(u64, u64, u64)] = &[
            (7,    1_026_105, 636_225),
            (100,  659_291,   755_487),
            (1000, 925_418,   411_028),
        ];
        for &(k, ref_x, ref_y) in cases {
            let result = c.scalar_mul(&g, &Uint::<4>::from(k));
            match &result {
                AffinePoint::Infinity => panic!("tiny_a: {k}·G is ∞"),
                AffinePoint::Finite { x, y } => {
                    assert_eq!(x.to_uint(), Uint::<4>::from(ref_x), "tiny_a: {k}·G x mismatch");
                    assert_eq!(y.to_uint(), Uint::<4>::from(ref_y), "tiny_a: {k}·G y mismatch");
                }
            }
        }
    }

    #[test]
    fn tiny_b_generator_on_curve() {
        let c = tiny_b();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "tiny_b: generator not on curve");
    }

    #[test]
    fn tiny_b_n_times_g_is_infinity() {
        let c = tiny_b();
        let g: AffinePoint<FpMonty> = c.generator();
        let ng = c.scalar_mul(&g, &c.n);
        assert!(ng.is_infinity(), "tiny_b: n·G should be ∞");
    }

    /// Known DLP values computed by the Python reference implementation.
    ///
    /// | k | Q.x | Q.y |
    /// |---|-----|-----|
    /// | 7 | 284_547 | 163_192 |
    /// | 42 | 132_859 | 318_692 |
    /// | 99991 | 654_745 | 751_943 |
    #[test]
    fn tiny_b_scalar_mul_reference() {
        let c = tiny_b();
        let g: AffinePoint<FpMonty> = c.generator();
        let cases: &[(u64, u64, u64)] = &[
            (7,     284_547, 163_192),
            (42,    132_859, 318_692),
            (99991, 654_745, 751_943),
        ];
        for &(k, ref_x, ref_y) in cases {
            let result = c.scalar_mul(&g, &Uint::<4>::from(k));
            match &result {
                AffinePoint::Infinity => panic!("tiny_b: {k}·G is ∞"),
                AffinePoint::Finite { x, y } => {
                    assert_eq!(x.to_uint(), Uint::<4>::from(ref_x), "tiny_b: {k}·G x mismatch");
                    assert_eq!(y.to_uint(), Uint::<4>::from(ref_y), "tiny_b: {k}·G y mismatch");
                }
            }
        }
    }

    /// tiny_glv generator is on the curve.
    #[test]
    fn tiny_glv_generator_on_curve() {
        let c = tiny_glv();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "tiny_glv: generator not on curve");
    }

    /// n·G = ∞ for tiny_glv.
    #[test]
    fn tiny_glv_n_times_g_is_infinity() {
        let c = tiny_glv();
        let g: AffinePoint<FpMonty> = c.generator();
        let ng = c.scalar_mul(&g, &c.n);
        assert!(ng.is_infinity(), "tiny_glv: n·G should be ∞");
    }

    /// φ(G) = λ·G for tiny_glv (verifies the GLV constants are correct).
    #[test]
    fn tiny_glv_endomorphism_matches_lambda() {
        use crate::ecdlp::glv::glv_phi;
        let c = tiny_glv();
        let g: AffinePoint<FpMonty> = c.generator();
        let phi_g = glv_phi(&g, &c.p, TINY_GLV_BETA);
        let lam_g = c.scalar_mul(&g, &Uint::<4>::from(TINY_GLV_LAMBDA));
        assert_eq!(phi_g, lam_g, "tiny_glv: φ(G) ≠ λ·G");
    }

    // ── composite_toy tests ───────────────────────────────────────────────────

    /// composite_toy generator is on the curve.
    #[test]
    fn composite_toy_generator_on_curve() {
        let c = composite_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "composite_toy: generator not on curve");
    }

    /// n·G = ∞ for composite_toy.
    #[test]
    fn composite_toy_n_times_g_is_infinity() {
        let c = composite_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        let ng = c.scalar_mul(&g, &c.n);
        assert!(ng.is_infinity(), "composite_toy: n·G should be ∞");
    }

    /// Full-order check: [n/pᵢ]·G ≠ ∞ for every prime pᵢ | n.
    ///
    /// Verifies that G has order exactly n (not a proper divisor).
    /// Primes dividing n = 60 = 2² · 3 · 5 are {2, 3, 5}.
    #[test]
    fn composite_toy_generator_full_order() {
        let c = composite_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        let n = COMPOSITE_TOY_N;
        // [n/2]·G = [30]·G must not be ∞
        let g30 = c.scalar_mul(&g, &Uint::<4>::from(n / 2));
        assert!(!g30.is_infinity(), "composite_toy: [n/2]·G should not be ∞ (full-order check)");
        // [n/3]·G = [20]·G must not be ∞
        let g20 = c.scalar_mul(&g, &Uint::<4>::from(n / 3));
        assert!(!g20.is_infinity(), "composite_toy: [n/3]·G should not be ∞ (full-order check)");
        // [n/5]·G = [12]·G must not be ∞
        let g12 = c.scalar_mul(&g, &Uint::<4>::from(n / 5));
        assert!(!g12.is_infinity(), "composite_toy: [n/5]·G should not be ∞ (full-order check)");
    }

    /// COMPOSITE_TOY_FACTORS product equals COMPOSITE_TOY_N.
    #[test]
    fn composite_toy_factors_product() {
        let product: u64 = COMPOSITE_TOY_FACTORS.iter().map(|&(p, e)| p.pow(e)).product();
        assert_eq!(
            product, COMPOSITE_TOY_N,
            "∏ pᵢ^eᵢ = {product}, expected {COMPOSITE_TOY_N}"
        );
    }

    /// Known scalar multiples of G on composite_toy, computed offline.
    ///
    /// | k  | Q.x | Q.y |
    /// |----|-----|-----|
    /// | 1  | 10  | 3   |
    /// | 7  | 24  | 43  |
    /// | 12 | 16  | 44  |
    /// | 20 | 8   | 41  |
    /// | 30 | 38  | 0   |
    #[test]
    fn composite_toy_scalar_mul_reference() {
        let c = composite_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        let cases: &[(u64, u64, u64)] = &[
            (1,  10, 3),
            (7,  24, 43),
            (12, 16, 44),
            (20, 8,  41),
        ];
        for &(k, ref_x, ref_y) in cases {
            let result = c.scalar_mul(&g, &Uint::<4>::from(k));
            match &result {
                AffinePoint::Infinity => panic!("composite_toy: {k}·G is ∞"),
                AffinePoint::Finite { x, y } => {
                    assert_eq!(
                        x.to_uint(), Uint::<4>::from(ref_x),
                        "composite_toy: {k}·G x mismatch"
                    );
                    assert_eq!(
                        y.to_uint(), Uint::<4>::from(ref_y),
                        "composite_toy: {k}·G y mismatch"
                    );
                }
            }
        }
    }

    // ── composite_rho_toy tests ───────────────────────────────────────────────

    #[test]
    fn composite_rho_toy_generator_on_curve() {
        let c = composite_rho_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "composite_rho_toy: generator not on curve");
    }

    /// `G` has order exactly `n`: `n·G = ∞` and `[n/p]·G ≠ ∞` for every prime `p | n`.
    #[test]
    fn composite_rho_toy_generator_full_order() {
        let c = composite_rho_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.scalar_mul(&g, &c.n).is_infinity(), "n·G should be ∞");
        for &(p, _) in COMPOSITE_RHO_TOY_FACTORS {
            let cofactor = Uint::<4>::from(COMPOSITE_RHO_TOY_N / p);
            assert!(!c.scalar_mul(&g, &cofactor).is_infinity(), "[n/{p}]·G should not be ∞");
        }
    }

    /// The recorded factorization multiplies back to `n`, and one prime exceeds the
    /// enumeration limit so the rho branch is reachable.
    #[test]
    fn composite_rho_toy_factors_product_and_rho_prime() {
        let product: u64 = COMPOSITE_RHO_TOY_FACTORS.iter().map(|&(p, e)| p.pow(e)).product();
        assert_eq!(product, COMPOSITE_RHO_TOY_N);
        assert!(COMPOSITE_RHO_TOY_FACTORS.iter().any(|&(p, _)| p > 64));
    }

    /// [30]·G = (38, 0) — a point of order 2 (y = 0), confirming the 2² factor.
    #[test]
    fn composite_toy_order_2_point() {
        let c = composite_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        // [30]·G should be (38, 0), which has order 2 (y=0 means 2P = ∞)
        let g30 = c.scalar_mul(&g, &Uint::<4>::from(30u64));
        match &g30 {
            AffinePoint::Infinity => panic!("composite_toy: [30]·G should not be ∞"),
            AffinePoint::Finite { x, y } => {
                assert_eq!(x.to_uint(), Uint::<4>::from(38u64), "[30]·G x mismatch");
                assert_eq!(y.to_uint(), Uint::<4>::from(0u64),  "[30]·G y should be 0 (order-2 point)");
            }
        }
        // Doubling it gives ∞
        let g60 = c.scalar_mul(&g, &Uint::<4>::from(60u64));
        assert!(g60.is_infinity(), "composite_toy: [60]·G should be ∞");
    }
}
