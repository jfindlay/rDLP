//! Downsized GLV-friendly secp256k1-style curve.
//!
//! This is a 63-bit prime-field curve with `a = 0` and the CM endomorphism
//! `φ(x, y) = (β·x mod p, y)`, which satisfies `φ(P) = λ·P` for all curve
//! points P.  The structure mirrors secp256k1 (`y² = x³ + 7`, j-invariant 0,
//! CM by ℤ[ω₃]) but over a much smaller prime — suitable for pedagogical and
//! benchmark experiments.
//!
//! # Parameters
//!
//! | Symbol | Value | Notes |
//! |--------|-------|-------|
//! | p | `4_611_686_018_427_395_203` | 63-bit prime, `p ≡ 1 (mod 3)` |
//! | a | `0` | secp256k1-style |
//! | b | `7` | secp256k1-style |
//! | n | `4_611_686_022_420_787_627` | prime group order, `n ≡ 1 (mod 3)` |
//! | Gx | `2` | base point x-coordinate |
//! | Gy | `3_236_101_131_256_320_111` | base point y-coordinate |
//! | β | `2_535_098_114_878_923_204` | cube root of unity mod p (`β³ ≡ 1 mod p`) |
//! | λ | `441_215_077_713_529_363` | GLV scalar (`λ² + λ + 1 ≡ 0 mod n`) |
//!
//! # GLV endomorphism
//!
//! For any point P = (x, y), `φ(P) = (β·x mod p, y)` and `φ(P) = λ·P`.
//! This lets any scalar k be split as `k = k₁ + k₂·λ (mod n)` with
//! `k₁, k₂ ≈ √n`, halving the effective scalar-multiplication cost (GLV endomorphism optimization).
//!
//! # Parameter derivation
//!
//! The curve was constructed via the CM method:
//!
//! 1. Find prime `p ≡ 1 (mod 3)` via `4p = t² + 3v²` (Cornacchia's algorithm).
//! 2. Set `n = p + 1 − t`; verify `n` is prime and `n ≡ 1 (mod 3)`.
//! 3. Find `β` as an element of order 3 in `(ℤ/pℤ)*` via `β = g^((p−1)/3)`.
//! 4. Find `λ` as a root of `x² + x + 1 ≡ 0 (mod n)` via the quadratic formula.
//! 5. Select the `(β, λ)` pair where `φ(G) = λ·G`.
//! 6. Verify: `Gy² ≡ Gx³ + 7 (mod p)`, `β³ ≡ 1 (mod p)`,
//!    `λ² + λ + 1 ≡ 0 (mod n)`, `φ(G) = λ·G`.

use crypto_bigint::Uint;
use crate::curve::Curve;

/// `p = 4_611_686_018_427_395_203 = 0x4000_0000_0000_1C83`.
pub const P: u64 = 4_611_686_018_427_395_203;
/// Prime group order `n = 4_611_686_022_420_787_627`.
pub const N: u64 = 4_611_686_022_420_787_627;
/// Generator x-coordinate.
pub const GX: u64 = 2;
/// Generator y-coordinate.
pub const GY: u64 = 3_236_101_131_256_320_111;
/// Cube root of unity mod p: `β = 2_535_098_114_878_923_204`.
///
/// Satisfies β³ ≡ 1 (mod p) and β ≠ 1.  The GLV endomorphism is
/// `φ(x, y) = (β·x mod p, y)`.
pub const BETA: u64 = 2_535_098_114_878_923_204;
/// GLV scalar eigenvalue: `λ = 441_215_077_713_529_363`.
///
/// Satisfies λ² + λ + 1 ≡ 0 (mod n) and `φ(P) = λ·P` for all P.
pub const LAMBDA: u64 = 441_215_077_713_529_363;

/// Return the secp256k1-toy GLV curve.
///
/// `y² = x³ + 7 mod p`,  `p = 4_611_686_018_427_395_203`.
pub fn secp_k1_toy() -> Curve {
    Curve {
        p:  Uint::<4>::from(P),
        a:  Uint::<4>::ZERO,
        b:  Uint::<4>::from(7u64),
        n:  Uint::<4>::from(N),
        gx: Uint::<4>::from(GX),
        gy: Uint::<4>::from(GY),
    }
}

/// Apply the GLV endomorphism to an x-coordinate: `x ↦ β·x mod p`.
///
/// Used in the GLV endomorphism optimization to evaluate `φ(P) = (glv_phi_x(P.x), P.y)` cheaply
/// (one multiplication, no inversion).
///
/// The `p` argument must equal [`P`]; it is accepted explicitly so the
/// caller can use whatever field type it already has in scope.
pub fn glv_phi_x(x_uint: Uint<4>) -> Uint<4> {
    // Compute β·x mod p using widening multiplication.
    let beta_u = Uint::<4>::from(BETA);
    let p_u    = Uint::<4>::from(P);
    let (lo, hi) = beta_u.widening_mul(&x_uint);
    let wide = crypto_bigint::Uint::<8>::from((lo, hi));
    let p8   = crypto_bigint::Uint::<8>::from((p_u, Uint::<4>::ZERO));
    let nz   = crypto_bigint::NonZero::new(p8).expect("p is non-zero");
    let rem8 = wide.rem(&nz);
    let words = rem8.as_words();
    Uint::<4>::from_words([words[0], words[1], words[2], words[3]])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::{Fp, FpMonty4 as FpMonty};

    use crate::curve::AffinePoint;

    // Reference points for `y² = x³ + 7 mod p`, computed by independent
    // Python affine-coordinate reference implementation:
    //
    //   1*G = (2,                    3_236_101_131_256_320_111)
    //   2*G = (922_337_203_685_479_039, 132_612_412_593_110_192)
    //   3*G = (2_732_850_973_882_900_861, 4_393_719_944_955_491_326)
    //   7*G = (1_979_583_965_183_108_279, 3_698_824_691_131_872_996)
    //  λ*G = (458_510_211_330_451_205, 3_236_101_131_256_320_111)

    #[test]
    fn generator_on_curve() {
        let c = secp_k1_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        assert!(c.is_on_curve(&g), "generator not on curve");
    }

    #[test]
    fn n_times_generator_is_infinity() {
        let c = secp_k1_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        let n = Uint::<4>::from(N);
        let ng = c.scalar_mul(&g, &n);
        assert!(ng.is_infinity(), "n·G should be the point at infinity");
    }

    #[test]
    fn scalar_mul_small_matches_reference() {
        let c = secp_k1_toy();
        let g: AffinePoint<FpMonty> = c.generator();

        let cases: &[(u64, u64, u64)] = &[
            (1, 2, 3_236_101_131_256_320_111),
            (2, 922_337_203_685_479_039, 132_612_412_593_110_192),
            (3, 2_732_850_973_882_900_861, 4_393_719_944_955_491_326),
            (7, 1_979_583_965_183_108_279, 3_698_824_691_131_872_996),
        ];

        for &(k, ref_x, ref_y) in cases {
            let result = c.scalar_mul(&g, &Uint::<4>::from(k));
            match &result {
                AffinePoint::Infinity => panic!("{k}·G should not be infinity"),
                AffinePoint::Finite { x, y } => {
                    assert_eq!(
                        x.to_uint(),
                        Uint::<4>::from(ref_x),
                        "{k}·G x mismatch"
                    );
                    assert_eq!(
                        y.to_uint(),
                        Uint::<4>::from(ref_y),
                        "{k}·G y mismatch"
                    );
                }
            }
        }
    }

    #[test]
    fn endomorphism_matches_lambda_scalar_mul() {
        // phi(G) = (beta * Gx mod p, Gy) must equal lambda * G.
        let c = secp_k1_toy();
        let g: AffinePoint<FpMonty> = c.generator();
        let p = &c.p;

        // Compute lambda * G via scalar multiplication.
        let lam_g = c.scalar_mul(&g, &Uint::<4>::from(LAMBDA));

        // Compute phi(G) = (beta * Gx mod p, Gy) directly.
        let phi_x = glv_phi_x(Uint::<4>::from(GX));
        let phi_g: AffinePoint<FpMonty> = AffinePoint::Finite {
            x: FpMonty::from_uint(phi_x, p),
            y: g.y().unwrap().clone(),
        };

        assert_eq!(lam_g, phi_g, "φ(G) ≠ λ·G — endomorphism eigenvalue mismatch");
    }

    #[test]
    fn endomorphism_is_on_curve() {
        let c = secp_k1_toy();
        let lam_g = c.scalar_mul::<FpMonty>(
            &c.generator(),
            &Uint::<4>::from(LAMBDA),
        );
        assert!(c.is_on_curve(&lam_g), "λ·G not on curve");
    }

    #[test]
    fn beta_cube_is_one() {
        // β³ ≡ 1 (mod p) — sanity-check the constant.
        let beta = Uint::<4>::from(BETA);
        let p    = Uint::<4>::from(P);
        let three = Uint::<4>::from(3u64);
        use shared_field::{Fp, FpMonty4 as FpMonty};
        let beta_fp = FpMonty::from_uint(beta, &p);
        let beta3   = beta_fp.pow(&three, &p);
        assert_eq!(beta3, FpMonty::one(&p), "β³ ≢ 1 (mod p)");
    }

    #[test]
    fn lambda_satisfies_minimal_polynomial() {
        // λ² + λ + 1 ≡ 0 (mod n) — sanity-check the constant.
        // Compute lam^2 + lam + 1 mod n using u128 to avoid overflow.
        let l  = LAMBDA as u128;
        let nl = N as u128;
        let val = ((l * l % nl) + l + 1) % nl;
        assert_eq!(val, 0, "λ² + λ + 1 ≢ 0 (mod n)");
    }
}
