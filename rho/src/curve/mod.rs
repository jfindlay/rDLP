//! Elliptic curve group law over GF(p).
//!
//! This module defines:
//! - [`Curve`] — a short Weierstrass curve y² = x³ + ax + b over GF(p).
//! - [`AffinePoint`] — a point in affine coordinates (x, y) or the point at infinity.
//! - [`JacobianPoint`] — a point in Jacobian projective coordinates (X:Y:Z) where
//!   the affine representative is (X/Z², Y/Z³). Efficient for doubling chains.
//!
//! # Design
//!
//! Curve arithmetic is generic over any prime field `F` implementing [`crate::field::Fp`].
//! This lets the same group-law code run with `FpNaive` (pedagogical baseline) or `FpMonty`
//! (benchmarked fast path), directly demonstrating the field-level speedup flowing up through
//! curve arithmetic.
//!
//! Two concrete curves are defined in sub-modules:
//! - [`generic`] — a Weierstrass curve over a 63-bit prime field whose stored generator has prime
//!   order 561_103 (a 20-bit subgroup, about 2¹⁹). It is a runnable demonstration fixture, not a
//!   large-group performance baseline; base and target points must lie in that subgroup.
//! - [`secp_k1_toy`] — a downsized GLV-friendly secp256k1-style curve (63-bit prime order) with
//!   an explicit order-3 endomorphism (λ, β) for the GLV endomorphism optimization.

pub mod generic;
pub mod secp_k1_toy;
pub mod test_curves;

use crypto_bigint::Uint;
use shared_field::Fp;

// ── Point types ──────────────────────────────────────────────────────────────

/// A point on a short Weierstrass curve in affine coordinates.
///
/// The point at infinity is represented by the `Infinity` variant. All finite
/// points satisfy the curve equation y² ≡ x³ + ax + b (mod p).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AffinePoint<F> {
    /// The point at infinity (group identity).
    Infinity,
    /// A finite affine point (x, y).
    Finite {
        /// x-coordinate (reduced mod p).
        x: F,
        /// y-coordinate (reduced mod p).
        y: F,
    },
}

impl<F: Fp<4>> AffinePoint<F> {
    /// Construct a finite affine point. Does not check the curve equation.
    #[inline]
    pub fn new(x: F, y: F) -> Self {
        AffinePoint::Finite { x, y }
    }

    /// Return `true` if this is the point at infinity.
    #[inline]
    pub fn is_infinity(&self) -> bool {
        matches!(self, AffinePoint::Infinity)
    }

    /// Negate a point: (x, y) → (x, -y); ∞ → ∞.
    pub fn negate(&self, p: &Uint<4>) -> Self {
        match self {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => AffinePoint::Finite {
                x: x.clone(),
                y: y.neg(p),
            },
        }
    }

    /// Return the x-coordinate, or `None` for the point at infinity.
    pub fn x(&self) -> Option<&F> {
        match self {
            AffinePoint::Infinity => None,
            AffinePoint::Finite { x, .. } => Some(x),
        }
    }

    /// Return the y-coordinate, or `None` for the point at infinity.
    pub fn y(&self) -> Option<&F> {
        match self {
            AffinePoint::Infinity => None,
            AffinePoint::Finite { y, .. } => Some(y),
        }
    }
}

/// A point in Jacobian projective coordinates (X : Y : Z).
///
/// The corresponding affine point is (X/Z², Y/Z³) when Z ≠ 0. Z = 0 represents
/// the point at infinity. Jacobian coordinates avoid field inversions during
/// point addition and doubling — one inversion is deferred until the final
/// `to_affine` conversion.
///
/// Formulae follow Hankerson-Menezes-Vanstone §3.2.2 (complete for distinct inputs).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JacobianPoint<F> {
    /// X-coordinate numerator.
    pub x: F,
    /// Y-coordinate numerator.
    pub y: F,
    /// Denominator factor Z (Z=0 ↔ point at infinity).
    pub z: F,
}

impl<F: Fp<4>> JacobianPoint<F> {
    /// Construct the point at infinity.
    pub fn infinity(p: &Uint<4>) -> Self {
        JacobianPoint {
            x: F::one(p),
            y: F::one(p),
            z: F::zero(p),
        }
    }

    /// Return `true` if this is the point at infinity (Z = 0).
    pub fn is_infinity(&self, p: &Uint<4>) -> bool {
        self.z.is_zero(p)
    }

    /// Convert an affine point to Jacobian (Z = 1 for finite points).
    pub fn from_affine(pt: &AffinePoint<F>, p: &Uint<4>) -> Self {
        match pt {
            AffinePoint::Infinity => Self::infinity(p),
            AffinePoint::Finite { x, y } => JacobianPoint {
                x: x.clone(),
                y: y.clone(),
                z: F::one(p),
            },
        }
    }

    /// Convert this Jacobian point to affine, performing one field inversion.
    ///
    /// Returns `AffinePoint::Infinity` if Z = 0.
    pub fn to_affine(&self, p: &Uint<4>) -> AffinePoint<F> {
        if self.is_infinity(p) {
            return AffinePoint::Infinity;
        }
        let z_inv = self.z.inv(p);
        let z_inv2 = z_inv.square(p);
        let z_inv3 = z_inv2.mul(&z_inv, p);
        let x = self.x.mul(&z_inv2, p);
        let y = self.y.mul(&z_inv3, p);
        AffinePoint::Finite { x, y }
    }
}

// ── Curve ────────────────────────────────────────────────────────────────────

/// A short Weierstrass curve y² = x³ + ax + b over GF(p).
///
/// Parameters are stored as `Uint<4>` field constants and converted to `F` on demand.
/// The curve does not own `PhantomData<F>` — operations are generic methods that take
/// the field type as a type parameter on each call.
#[derive(Clone, Debug)]
pub struct Curve {
    /// Field characteristic (prime).
    pub p: Uint<4>,
    /// Curve coefficient a.
    pub a: Uint<4>,
    /// Curve coefficient b.
    pub b: Uint<4>,
    /// Order of the base point (prime group order n).
    pub n: Uint<4>,
    /// Base point x-coordinate (generator G).
    pub gx: Uint<4>,
    /// Base point y-coordinate (generator G).
    pub gy: Uint<4>,
}

impl Curve {
    /// Return the base point G as an affine point over field `F`.
    pub fn generator<F: Fp<4>>(&self) -> AffinePoint<F> {
        AffinePoint::Finite {
            x: F::from_uint(self.gx, &self.p),
            y: F::from_uint(self.gy, &self.p),
        }
    }

    /// Check whether `pt` satisfies the curve equation y² = x³ + ax + b mod p.
    pub fn is_on_curve<F: Fp<4>>(&self, pt: &AffinePoint<F>) -> bool {
        match pt {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => {
                let p = &self.p;
                let a = F::from_uint(self.a, p);
                let b = F::from_uint(self.b, p);
                let lhs = y.square(p);
                let rhs = x.square(p).mul(x, p)          // x³
                    .add(&a.mul(x, p), p)                 // + ax
                    .add(&b, p);                          // + b
                lhs == rhs
            }
        }
    }

    /// Negate a point: (x, y) → (x, −y); ∞ → ∞.
    pub fn negate<F: Fp<4>>(&self, pt: &AffinePoint<F>) -> AffinePoint<F> {
        pt.negate(&self.p)
    }

    // ── Jacobian group law ────────────────────────────────────────────────

    /// Point doubling in Jacobian coordinates: 2P.
    ///
    /// Uses the standard Jacobian doubling formulae (4M + 4S + 6add):
    ///   S = 4·X·Y²
    ///   M = 3·X² + a·Z⁴
    ///   X' = M² − 2S
    ///   Y' = M·(S − X') − 8·Y⁴
    ///   Z' = 2·Y·Z
    pub fn double_jacobian<F: Fp<4>>(&self, pt: &JacobianPoint<F>) -> JacobianPoint<F> {
        let p = &self.p;
        if pt.is_infinity(p) {
            return JacobianPoint::infinity(p);
        }

        let a = F::from_uint(self.a, p);

        // S = 4·X·Y²
        let y2 = pt.y.square(p);
        let s = F::from_u64(4, p).mul(&pt.x.mul(&y2, p), p);

        // M = 3·X² + a·Z⁴
        let x2 = pt.x.square(p);
        let z2 = pt.z.square(p);
        let z4 = z2.square(p);
        let m = F::from_u64(3, p).mul(&x2, p).add(&a.mul(&z4, p), p);

        // X' = M² − 2S
        let m2 = m.square(p);
        let two_s = s.add(&s, p);
        let x3 = m2.sub(&two_s, p);

        // Y' = M·(S − X') − 8·Y⁴
        let y4 = y2.square(p);
        let eight_y4 = F::from_u64(8, p).mul(&y4, p);
        let y3 = m.mul(&s.sub(&x3, p), p).sub(&eight_y4, p);

        // Z' = 2·Y·Z
        let z3 = F::from_u64(2, p).mul(&pt.y.mul(&pt.z, p), p);

        JacobianPoint { x: x3, y: y3, z: z3 }
    }

    /// Point addition in Jacobian + Jacobian coordinates: P + Q.
    ///
    /// Uses the Jacobian mixed-addition formulae (16M + 4S). Falls back to
    /// doubling when P = Q, and returns the correct identity when P = −Q.
    ///
    /// Reference: HMV §3.2.2, Algorithm 3.22.
    pub fn add_jacobian<F: Fp<4>>(
        &self,
        p1: &JacobianPoint<F>,
        p2: &JacobianPoint<F>,
    ) -> JacobianPoint<F> {
        let p = &self.p;

        if p1.is_infinity(p) {
            return p2.clone();
        }
        if p2.is_infinity(p) {
            return p1.clone();
        }

        // U1 = X1·Z2², U2 = X2·Z1²
        let z1_2 = p1.z.square(p);
        let z2_2 = p2.z.square(p);
        let u1 = p1.x.mul(&z2_2, p);
        let u2 = p2.x.mul(&z1_2, p);

        // S1 = Y1·Z2³, S2 = Y2·Z1³
        let z1_3 = z1_2.mul(&p1.z, p);
        let z2_3 = z2_2.mul(&p2.z, p);
        let s1 = p1.y.mul(&z2_3, p);
        let s2 = p2.y.mul(&z1_3, p);

        let h = u2.sub(&u1, p);  // H = U2 − U1
        let r = s2.sub(&s1, p);  // R = S2 − S1

        if h.is_zero(p) {
            if r.is_zero(p) {
                // P1 = P2: use doubling
                return self.double_jacobian(p1);
            } else {
                // P1 = −P2: return infinity
                return JacobianPoint::infinity(p);
            }
        }

        // X3 = R² − H³ − 2·U1·H²
        let h2 = h.square(p);
        let h3 = h2.mul(&h, p);
        let u1h2 = u1.mul(&h2, p);
        let two_u1h2 = u1h2.add(&u1h2, p);
        let x3 = r.square(p).sub(&h3, p).sub(&two_u1h2, p);

        // Y3 = R·(U1·H² − X3) − S1·H³
        let y3 = r.mul(&u1h2.sub(&x3, p), p).sub(&s1.mul(&h3, p), p);

        // Z3 = H·Z1·Z2
        let z3 = h.mul(&p1.z.mul(&p2.z, p), p);

        JacobianPoint { x: x3, y: y3, z: z3 }
    }

    /// Mixed Jacobian + affine addition: P_jac + Q_aff.
    ///
    /// When the second operand has Z = 1, several multiplications drop out.
    /// Uses 8M + 3S (vs 16M + 4S for full Jacobian addition). This is the hot
    /// path in scalar multiplication and r-adding walks.
    ///
    /// Reference: HMV §3.2.2, Algorithm 3.23.
    pub fn add_mixed<F: Fp<4>>(
        &self,
        p1: &JacobianPoint<F>,
        p2: &AffinePoint<F>,
    ) -> JacobianPoint<F> {
        let p = &self.p;

        if p1.is_infinity(p) {
            return JacobianPoint::from_affine(p2, p);
        }
        if p2.is_infinity() {
            return p1.clone();
        }
        let (x2, y2) = match p2 {
            AffinePoint::Finite { x, y } => (x, y),
            AffinePoint::Infinity => unreachable!(),
        };

        // Z1² and Z1³
        let z1_2 = p1.z.square(p);
        let z1_3 = z1_2.mul(&p1.z, p);

        // U2 = X2·Z1²,  S2 = Y2·Z1³
        let u2 = x2.mul(&z1_2, p);
        let s2 = y2.mul(&z1_3, p);

        let h = u2.sub(&p1.x, p);
        let r = s2.sub(&p1.y, p);

        if h.is_zero(p) {
            if r.is_zero(p) {
                return self.double_jacobian(p1);
            } else {
                return JacobianPoint::infinity(p);
            }
        }

        let h2 = h.square(p);
        let h3 = h2.mul(&h, p);
        let x1h2 = p1.x.mul(&h2, p);
        let two_x1h2 = x1h2.add(&x1h2, p);

        let x3 = r.square(p).sub(&h3, p).sub(&two_x1h2, p);
        let y3 = r.mul(&x1h2.sub(&x3, p), p).sub(&p1.y.mul(&h3, p), p);
        let z3 = h.mul(&p1.z, p);

        JacobianPoint { x: x3, y: y3, z: z3 }
    }

    /// Scalar multiplication: k·P using double-and-add (left-to-right binary method).
    ///
    /// `scalar` is a `Uint<4>`; the high bits are consumed first. Uses mixed
    /// Jacobian+affine addition for the hot path (each add uses the affine input P).
    ///
    /// Returns the result in affine coordinates (one final inversion).
    pub fn scalar_mul<F: Fp<4>>(&self, pt: &AffinePoint<F>, scalar: &Uint<4>) -> AffinePoint<F> {
        let p = &self.p;
        let mut result = JacobianPoint::infinity(p);
        let bits = 256u32;

        for i in (0..bits).rev() {
            result = self.double_jacobian(&result);
            if scalar.bit(i).into() {
                result = self.add_mixed(&result, pt);
            }
        }

        result.to_affine(p)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_field::FpMonty4 as FpMonty;

    /// A tiny 5-point curve: y² = x³ + x + 1 over GF(7).
    ///
    /// Points: ∞, (2,2), (2,5), (4,0), (0,1) — wait, let's be precise.
    ///
    /// Over GF(7), y² = x³ + x + 1:
    ///   x=0: y²=1  → y=1,6
    ///   x=1: y²=3  → no solution
    ///   x=2: y²=11=4 → y=2,5
    ///   x=3: y²=31=3 → no solution
    ///   x=4: y²=69=6 → no solution
    ///   x=5: y²=131=5 → no solution
    ///   x=6: y²=229=5 → no solution
    /// So finite points: (0,1),(0,6),(2,2),(2,5) — a group of order 5.
    fn tiny_curve() -> Curve {
        let seven = Uint::<4>::from(7u64);
        // a=1, b=1, p=7
        // Generator (0,1), order=5 (prime)
        Curve {
            p: seven,
            a: Uint::<4>::ONE,
            b: Uint::<4>::ONE,
            n: Uint::<4>::from(5u64),
            gx: Uint::<4>::ZERO,
            gy: Uint::<4>::ONE,
        }
    }

    #[test]
    fn generator_is_on_curve() {
        let c = tiny_curve();
        let g = c.generator::<FpMonty>();
        assert!(c.is_on_curve(&g), "generator not on curve");
    }

    #[test]
    fn double_generator_is_on_curve() {
        let c = tiny_curve();
        let p = &c.p;
        let g = c.generator::<FpMonty>();
        let gj = JacobianPoint::from_affine(&g, p);
        let two_g = c.double_jacobian(&gj).to_affine(p);
        assert!(c.is_on_curve(&two_g), "2G not on curve");
    }

    #[test]
    fn add_point_to_itself_matches_double() {
        let c = tiny_curve();
        let p = &c.p;
        let g = c.generator::<FpMonty>();
        let gj = JacobianPoint::from_affine(&g, p);
        let via_double = c.double_jacobian(&gj).to_affine(p);
        let via_add = c.add_jacobian(&gj, &gj).to_affine(p);
        assert_eq!(via_double, via_add, "double ≠ add(P,P)");
    }

    #[test]
    fn n_times_generator_is_infinity() {
        let c = tiny_curve();
        let g = c.generator::<FpMonty>();
        let n = c.n;
        let ng = c.scalar_mul(&g, &n);
        assert!(ng.is_infinity(), "n·G should be the point at infinity");
    }

    #[test]
    fn scalar_mul_one_is_identity() {
        let c = tiny_curve();
        let g = c.generator::<FpMonty>();
        let result = c.scalar_mul(&g, &Uint::<4>::ONE);
        assert_eq!(result, g, "1·G should equal G");
    }

    #[test]
    fn scalar_mul_two_matches_double() {
        let c = tiny_curve();
        let p = &c.p;
        let g = c.generator::<FpMonty>();
        let gj = JacobianPoint::from_affine(&g, p);
        let two_g_direct = c.double_jacobian(&gj).to_affine(p);
        let two_g_scalar = c.scalar_mul(&g, &Uint::<4>::from(2u64));
        assert_eq!(two_g_direct, two_g_scalar, "2·G via double ≠ via scalar_mul");
    }

    #[test]
    fn negate_then_add_is_infinity() {
        let c = tiny_curve();
        let p = &c.p;
        let g = c.generator::<FpMonty>();
        let neg_g = c.negate(&g);
        let gj = JacobianPoint::from_affine(&g, p);
        let neg_gj = JacobianPoint::from_affine(&neg_g, p);
        let sum = c.add_jacobian(&gj, &neg_gj).to_affine(p);
        assert!(sum.is_infinity(), "G + (−G) should be infinity");
    }

    #[test]
    fn add_mixed_matches_add_jacobian() {
        let c = tiny_curve();
        let p = &c.p;
        let g = c.generator::<FpMonty>();
        let two_g = c.scalar_mul(&g, &Uint::<4>::from(2u64));
        let two_gj = JacobianPoint::from_affine(&two_g, p);
        let gj = JacobianPoint::from_affine(&g, p);
        let via_jj = c.add_jacobian(&two_gj, &gj).to_affine(p);
        let via_mixed = c.add_mixed(&two_gj, &g).to_affine(p);
        assert_eq!(via_jj, via_mixed, "add_jacobian ≠ add_mixed for 2G+G");
    }
}
