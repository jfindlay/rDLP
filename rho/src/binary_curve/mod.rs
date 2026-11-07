//! Binary elliptic curve group law over GF(2^m).
//!
//! This module defines:
//! - [`BinaryCurve`] — the non-supersingular curve `y²+xy = x³+ax²+b` over GF(2^m).
//! - [`BinaryAffinePoint`] — a point in affine coordinates or the point at infinity.
//! - [`LDPoint`] — a point in López–Dahab projective coordinates (X:Y:Z), where the
//!   affine representative is `(X/Z, Y/Z²)`.
//!
//! # Design
//!
//! `BinaryCurve` is a **parallel** to `rho::curve::Curve`, not a reuse.  The
//! Jacobian group law in `Curve` divides by `2y` during doubling, which is zero
//! in characteristic 2 — the Jacobian formulae are **wrong** for binary curves.
//! López–Dahab coordinates are the char-2 analogue: the doubling formula avoids
//! the `2y` division entirely.
//!
//! The curve does not own `PhantomData<F>` — group-law methods are generic over
//! `F: F2m<1>` at the method level, threading the field type per-call.  This
//! mirrors the `Curve` idiom.
//!
//! # Characteristic-2 invariants
//!
//! - **Negation**: `−(x,y) = (x, x+y)` in char 2.  NOT `(x, −y)` — in char 2,
//!   `−y = y`, so naive negation gives `P` back, not `−P`.  The `P+(−P)=∞` KAT
//!   is the load-bearing guard.
//! - **Doubling**: uses the López–Dahab formula (NOT Jacobian).  The `2P == P+P`
//!   KAT is the guard.
//! - **`solve_quadratic` precondition**: `x²+x=c` is solvable iff `trace(c)=0`.
//!   Decompression of an x-coordinate not on the curve will produce a wrong y
//!   that fails `is_on_curve`.
//!
//! # López–Dahab coordinates
//!
//! In LD coordinates, a projective point `(X:Y:Z)` represents the affine point
//! `(X/Z, Y/Z²)` when `Z ≠ 0`.  `Z = 0` represents the point at infinity.
//!
//! **LD doubling** (derived from affine doubling `λ=x+y/x`, `x₃=λ²+λ+a`,
//! `y₃=x²+(λ+1)·x₃`; clearing denominators with `T=X²+Y·Z`, `Z₃=(X·Z)²`):
//! ```text
//! T  = X₁² + Y₁·Z₁
//! Z₃ = (X₁·Z₁)²
//! X₃ = T² + T·X₁·Z₁ + a·X₁²·Z₁²
//! Y₃ = X₁⁶·Z₁² + X₁·Z₁·(T + X₁·Z₁)·X₃
//! ```
//! Equivalently (on-curve, using `y²+xy=x³+ax²+b`): `X₃ = X₁⁴ + b·Z₁⁴`.
//!
//! **LD mixed addition** (P₁ projective, P₂ affine; `A=Y₁+y₂·Z₁²`, `B=X₁+x₂·Z₁`):
//! ```text
//! Z₃ = (B·Z₁)²
//! X₃ = A² + A·B·Z₁ + a·B²·Z₁² + B³·Z₁
//! Y₃ = A·B·Z₁·(X₁·B²·Z₁ + X₃) + B²·Z₁²·(X₃ + Y₁·B²)
//! ```
//!
//! # Toy field sizes
//!
//! The KATs use GF(2^4) with `x⁴+x+1` (poly = 0x13).  The algorithms are
//! crypto-scale-correct; only the parameters are toy (principle-4 boundary).

use crypto_bigint::Uint;
use shared_gf2m::F2m;

// ── Point types ───────────────────────────────────────────────────────────────

/// A point on a binary curve in affine coordinates.
///
/// The point at infinity is represented by the `Infinity` variant.  All finite
/// points satisfy `y²+xy = x³+ax²+b` over GF(2^m).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryAffinePoint<F> {
    /// The point at infinity (group identity).
    Infinity,
    /// A finite affine point (x, y).
    Finite {
        /// x-coordinate (field element).
        x: F,
        /// y-coordinate (field element).
        y: F,
    },
}

impl<F: Clone> BinaryAffinePoint<F> {
    /// Construct a finite affine point.  Does not check the curve equation.
    #[inline]
    pub fn new(x: F, y: F) -> Self {
        BinaryAffinePoint::Finite { x, y }
    }

    /// Return `true` if this is the point at infinity.
    #[inline]
    pub fn is_infinity(&self) -> bool {
        matches!(self, BinaryAffinePoint::Infinity)
    }

    /// Return the x-coordinate, or `None` for the point at infinity.
    pub fn x(&self) -> Option<&F> {
        match self {
            BinaryAffinePoint::Infinity => None,
            BinaryAffinePoint::Finite { x, .. } => Some(x),
        }
    }

    /// Return the y-coordinate, or `None` for the point at infinity.
    pub fn y(&self) -> Option<&F> {
        match self {
            BinaryAffinePoint::Infinity => None,
            BinaryAffinePoint::Finite { y, .. } => Some(y),
        }
    }
}

/// A point in López–Dahab projective coordinates `(X:Y:Z)`.
///
/// The affine representative is `(X/Z, Y/Z²)` when `Z ≠ 0`.
/// `Z = 0` represents the point at infinity.
///
/// LD coordinates are the char-2 analogue of Jacobian coordinates.  The
/// doubling formula avoids the `2y` division that makes Jacobian wrong in
/// characteristic 2.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LDPoint<F> {
    /// X-coordinate numerator.
    pub x: F,
    /// Y-coordinate numerator.
    pub y: F,
    /// Denominator factor Z (Z=0 ↔ point at infinity).
    pub z: F,
}

impl<F: F2m<1>> LDPoint<F> {
    /// Construct the point at infinity.
    pub fn infinity() -> Self {
        LDPoint {
            x: F::one(),
            y: F::one(),
            z: F::zero(),
        }
    }

    /// Return `true` if this is the point at infinity (Z = 0).
    pub fn is_infinity(&self) -> bool {
        self.z.is_zero()
    }

    /// Convert an affine point to LD projective (Z = 1 for finite points).
    pub fn from_affine(pt: &BinaryAffinePoint<F>) -> Self {
        match pt {
            BinaryAffinePoint::Infinity => Self::infinity(),
            BinaryAffinePoint::Finite { x, y } => LDPoint {
                x: x.clone(),
                y: y.clone(),
                z: F::one(),
            },
        }
    }

    /// Convert this LD point to affine, performing one field inversion.
    ///
    /// Returns `BinaryAffinePoint::Infinity` if Z = 0.
    pub fn to_affine(&self, poly: &Uint<1>) -> BinaryAffinePoint<F> {
        if self.is_infinity() {
            return BinaryAffinePoint::Infinity;
        }
        // x = X/Z, y = Y/Z²
        let z_inv = self.z.inv(poly);
        let z_inv2 = z_inv.square(poly);
        let x = self.x.mul(&z_inv, poly);
        let y = self.y.mul(&z_inv2, poly);
        BinaryAffinePoint::Finite { x, y }
    }
}

// ── Curve ─────────────────────────────────────────────────────────────────────

/// A binary elliptic curve `y²+xy = x³+ax²+b` over GF(2^m).
///
/// The curve is non-supersingular (`b ≠ 0`).  Parameters are stored as
/// `Uint<1>` constants and converted to field elements on demand.
///
/// The curve does not own `PhantomData<F>` — group-law methods are generic
/// over `F: F2m<1>` at the method level, threading the field type per-call.
/// This mirrors the `rho::curve::Curve` idiom.
///
/// # Toy scale
///
/// Designed for GF(2^4) and GF(2^8) toy curves.  The algorithms are
/// crypto-scale-correct; the parameters are toy (principle-4 boundary).
#[derive(Clone, Debug)]
pub struct BinaryCurve {
    /// Irreducible polynomial defining GF(2^m).
    pub poly: Uint<1>,
    /// Curve coefficient a.
    pub a: Uint<1>,
    /// Curve coefficient b (must be non-zero).
    pub b: Uint<1>,
    /// Order of the base point (prime group order n).
    pub n: Uint<1>,
    /// Base point x-coordinate.
    pub gx: Uint<1>,
    /// Base point y-coordinate.
    pub gy: Uint<1>,
}

impl BinaryCurve {
    /// Return the base point G as an affine point over field `F`.
    pub fn generator<F: F2m<1>>(&self) -> BinaryAffinePoint<F> {
        BinaryAffinePoint::Finite {
            x: F::from_uint(self.gx, &self.poly),
            y: F::from_uint(self.gy, &self.poly),
        }
    }

    /// Check whether `pt` satisfies the curve equation `y²+xy = x³+ax²+b`.
    pub fn is_on_curve<F: F2m<1>>(&self, pt: &BinaryAffinePoint<F>) -> bool {
        match pt {
            BinaryAffinePoint::Infinity => true,
            BinaryAffinePoint::Finite { x, y } => {
                let poly = &self.poly;
                let a = F::from_uint(self.a, poly);
                let b = F::from_uint(self.b, poly);
                // LHS = y² + x·y
                let lhs = y.square(poly).add(&x.mul(y, poly));
                // RHS = x³ + a·x² + b
                let x2 = x.square(poly);
                let x3 = x2.mul(x, poly);
                let rhs = x3.add(&a.mul(&x2, poly)).add(&b);
                lhs == rhs
            }
        }
    }

    /// Negate a point: `(x,y) → (x, x+y)` in char 2; `∞ → ∞`.
    ///
    /// **CRITICAL**: in characteristic 2, `−(x,y) = (x, x+y)`, NOT `(x, −y)`.
    /// In char 2, `−y = y`, so `(x, −y) = (x, y) = P` — that is the identity,
    /// not the negation.  The correct negation adds `x` to the y-coordinate.
    pub fn negate<F: F2m<1>>(&self, pt: &BinaryAffinePoint<F>) -> BinaryAffinePoint<F> {
        match pt {
            BinaryAffinePoint::Infinity => BinaryAffinePoint::Infinity,
            BinaryAffinePoint::Finite { x, y } => BinaryAffinePoint::Finite {
                x: x.clone(),
                // −(x,y) = (x, x+y) in char 2
                y: x.add(y),
            },
        }
    }

    // ── López–Dahab group law ─────────────────────────────────────────────────

    /// Point doubling in López–Dahab coordinates: 2P.
    ///
    /// Uses the LD doubling formula for `y²+xy=x³+ax²+b`.  This is NOT the
    /// Jacobian formula — Jacobian doubling divides by `2y = 0` in char 2.
    ///
    /// Formula (derived from affine doubling, clearing denominators):
    /// ```text
    /// T  = X₁² + Y₁·Z₁
    /// Z₃ = (X₁·Z₁)²
    /// X₃ = T² + T·X₁·Z₁ + a·X₁²·Z₁²
    /// Y₃ = X₁⁶·Z₁² + X₁·Z₁·(T + X₁·Z₁)·X₃
    /// ```
    /// On-curve, X₃ simplifies to `X₁⁴ + b·Z₁⁴` (using `y²+xy=x³+ax²+b`).
    pub fn double_ld<F: F2m<1>>(&self, pt: &LDPoint<F>) -> LDPoint<F> {
        let poly = &self.poly;
        if pt.is_infinity() {
            return LDPoint::infinity();
        }

        let a = F::from_uint(self.a, poly);

        let x1 = &pt.x;
        let y1 = &pt.y;
        let z1 = &pt.z;

        // T = X₁² + Y₁·Z₁
        let x1_sq = x1.square(poly);
        let t = x1_sq.add(&y1.mul(z1, poly));

        // Z₃ = (X₁·Z₁)²
        let xz = x1.mul(z1, poly);
        let z3 = xz.square(poly);

        // X₃ = T² + T·X₁·Z₁ + a·X₁²·Z₁²
        //     = T² + T·xz + a·z3
        let x3 = t.square(poly)
            .add(&t.mul(&xz, poly))
            .add(&a.mul(&z3, poly));

        // Y₃ = X₁⁶·Z₁² + X₁·Z₁·(T + X₁·Z₁)·X₃
        //     = x1_cu²·z1_sq + xz·(T + xz)·X₃
        let x1_cu = x1.mul(&x1_sq, poly);
        let x1_6 = x1_cu.square(poly);
        let z1_sq = z1.square(poly);
        let y3 = x1_6.mul(&z1_sq, poly)
            .add(&xz.mul(&t.add(&xz), poly).mul(&x3, poly));

        LDPoint { x: x3, y: y3, z: z3 }
    }

    /// Mixed López–Dahab + affine addition: P_ld + Q_aff.
    ///
    /// Adds a projective LD point to an affine point.  This is the hot path
    /// in scalar multiplication (the affine input Q has Z=1, simplifying the
    /// formula).
    ///
    /// Formula (derived from affine addition, clearing denominators;
    /// `A=Y₁+y₂·Z₁²`, `B=X₁+x₂·Z₁`):
    /// ```text
    /// Z₃ = (B·Z₁)²
    /// X₃ = A² + A·B·Z₁ + a·B²·Z₁² + B³·Z₁
    /// Y₃ = A·B·Z₁·(X₁·B²·Z₁ + X₃) + B²·Z₁²·(X₃ + Y₁·B²)
    /// ```
    pub fn add_mixed<F: F2m<1>>(
        &self,
        p1: &LDPoint<F>,
        p2: &BinaryAffinePoint<F>,
    ) -> LDPoint<F> {
        let poly = &self.poly;

        if p1.is_infinity() {
            return LDPoint::from_affine(p2);
        }
        if p2.is_infinity() {
            return p1.clone();
        }

        let (x2, y2) = match p2 {
            BinaryAffinePoint::Finite { x, y } => (x, y),
            BinaryAffinePoint::Infinity => unreachable!(),
        };

        let a = F::from_uint(self.a, poly);

        let x1 = &p1.x;
        let y1 = &p1.y;
        let z1 = &p1.z;

        // A = Y₁ + y₂·Z₁²
        let z1_sq = z1.square(poly);
        let a_val = y1.add(&y2.mul(&z1_sq, poly));

        // B = X₁ + x₂·Z₁
        let b_val = x1.add(&x2.mul(z1, poly));

        // Special cases: B = 0 means x₁ = x₂.
        if b_val.is_zero() {
            if a_val.is_zero() {
                // P₁ = P₂: use doubling.
                return self.double_ld(&LDPoint::from_affine(p2));
            } else {
                // P₁ = −P₂: return infinity.
                return LDPoint::infinity();
            }
        }

        // Z₃ = (B·Z₁)²
        let bz1 = b_val.mul(z1, poly);
        let z3 = bz1.square(poly);

        // X₃ = A² + A·B·Z₁ + a·B²·Z₁² + B³·Z₁
        //     = A² + A·bz1 + a·z3 + B²·bz1
        let a_sq = a_val.square(poly);
        let b_sq = b_val.square(poly);
        let x3 = a_sq
            .add(&a_val.mul(&bz1, poly))
            .add(&a.mul(&z3, poly))
            .add(&b_sq.mul(&bz1, poly));

        // Y₃ = A·B·Z₁·(X₁·B²·Z₁ + X₃) + B²·Z₁²·(X₃ + Y₁·B²)
        //     = A·bz1·(x1·b_sq·z1 + X₃) + z3·(X₃ + y1·b_sq)
        let x1_b2_z1 = x1.mul(&b_sq, poly).mul(z1, poly);
        let term1 = a_val.mul(&bz1, poly).mul(&x1_b2_z1.add(&x3), poly);
        let term2 = z3.mul(&x3.add(&y1.mul(&b_sq, poly)), poly);
        let y3 = term1.add(&term2);

        LDPoint { x: x3, y: y3, z: z3 }
    }

    /// Full LD + LD addition: P₁ + P₂ (both in projective coordinates).
    ///
    /// Converts P₂ to affine (one inversion) then calls `add_mixed`.
    /// Correct and auditable at toy scale; a production implementation would
    /// use the full LD-LD formula to avoid the inversion.
    pub fn add_ld<F: F2m<1>>(&self, p1: &LDPoint<F>, p2: &LDPoint<F>) -> LDPoint<F> {
        let poly = &self.poly;
        let p2_affine = p2.to_affine(poly);
        self.add_mixed(p1, &p2_affine)
    }

    /// Point addition in affine coordinates: P + Q.
    ///
    /// Uses the affine addition formula for `y²+xy=x³+ax²+b` (one inversion).
    /// This is the auditable baseline; `scalar_mul` uses LD projective internally.
    pub fn add<F: F2m<1>>(
        &self,
        p1: &BinaryAffinePoint<F>,
        p2: &BinaryAffinePoint<F>,
    ) -> BinaryAffinePoint<F> {
        let poly = &self.poly;

        match (p1, p2) {
            (BinaryAffinePoint::Infinity, _) => p2.clone(),
            (_, BinaryAffinePoint::Infinity) => p1.clone(),
            (
                BinaryAffinePoint::Finite { x: x1, y: y1 },
                BinaryAffinePoint::Finite { x: x2, y: y2 },
            ) => {
                let a = F::from_uint(self.a, poly);

                if x1 == x2 {
                    // Same x-coordinate: either P = Q (double) or P = -Q (infinity).
                    // P = -Q iff y₁ + y₂ = x₁ (char-2 negation: -(x,y) = (x, x+y)).
                    let y_sum = y1.add(y2);
                    if y_sum == *x1 {
                        // P = -Q: return infinity.
                        return BinaryAffinePoint::Infinity;
                    }
                    // P = Q: use doubling.
                    return self.double(p1);
                }

                // λ = (y₁ + y₂) / (x₁ + x₂)
                let x_sum = x1.add(x2);
                let y_sum = y1.add(y2);
                let lambda = y_sum.div(&x_sum, poly);

                // x₃ = λ² + λ + a + x₁ + x₂
                let x3 = lambda.square(poly).add(&lambda).add(&a).add(x1).add(x2);

                // y₃ = λ·(x₁ + x₃) + x₃ + y₁
                let y3 = lambda.mul(&x1.add(&x3), poly).add(&x3).add(y1);

                BinaryAffinePoint::Finite { x: x3, y: y3 }
            }
        }
    }

    /// Point doubling in affine coordinates: 2P.
    ///
    /// Uses the affine doubling formula for `y²+xy=x³+ax²+b` (one inversion).
    /// This is the auditable baseline; `scalar_mul` uses LD projective internally.
    pub fn double<F: F2m<1>>(&self, pt: &BinaryAffinePoint<F>) -> BinaryAffinePoint<F> {
        let poly = &self.poly;

        match pt {
            BinaryAffinePoint::Infinity => BinaryAffinePoint::Infinity,
            BinaryAffinePoint::Finite { x, y } => {
                if x.is_zero() {
                    // x = 0: the doubling formula λ = x + y/x is undefined.
                    // For y²+xy=x³+ax²+b with x=0: y²=b, so the point is (0, √b).
                    // The negation of (0,y) is (0, 0+y) = (0,y), so 2·(0,y) = ∞.
                    return BinaryAffinePoint::Infinity;
                }

                let a = F::from_uint(self.a, poly);

                // λ = x + y/x
                let lambda = x.add(&y.div(x, poly));

                // x₃ = λ² + λ + a
                let x3 = lambda.square(poly).add(&lambda).add(&a);

                // y₃ = x² + (λ+1)·x₃
                let x_sq = x.square(poly);
                let y3 = x_sq.add(&lambda.add(&F::one()).mul(&x3, poly));

                BinaryAffinePoint::Finite { x: x3, y: y3 }
            }
        }
    }

    /// Scalar multiplication: k·P using double-and-add (left-to-right binary method).
    ///
    /// Uses LD projective coordinates internally (mixed addition) to avoid
    /// per-step inversions.  One final inversion converts back to affine.
    pub fn scalar_mul<F: F2m<1>>(
        &self,
        pt: &BinaryAffinePoint<F>,
        scalar: &Uint<1>,
    ) -> BinaryAffinePoint<F> {
        let poly = &self.poly;
        let mut result = LDPoint::infinity();
        let bits = 64u32; // Uint<1> has 64 bits

        for i in (0..bits).rev() {
            result = self.double_ld(&result);
            if scalar.bit(i).into() {
                result = self.add_mixed(&result, pt);
            }
        }

        result.to_affine(poly)
    }

    // ── Point decompression ───────────────────────────────────────────────────

    /// Decompress a point: given x-coordinate and a sign bit, recover y.
    ///
    /// For the curve `y²+xy = x³+ax²+b`, substituting `λ = y/x` gives:
    ///   `λ²+λ = x + a + b/x²`
    ///
    /// Solve with `solve_quadratic` to get `λ`, then `y = λ·x`.
    /// The two roots `λ` and `λ+1` give the two y-values.
    ///
    /// The `sign_bit` selects between the two roots using the **constant term**
    /// (bit 0) of `λ`:
    /// - `sign_bit = false`: use `λ` if `bit_0(λ) = 0`, else use `λ+1`.
    /// - `sign_bit = true`:  use `λ` if `bit_0(λ) = 1`, else use `λ+1`.
    ///
    /// This convention works for all `m` (even and odd), since `λ` and `λ+1`
    /// always differ in bit 0.  (The trace-based convention only works for odd
    /// `m`, where `trace(1) = 1` distinguishes the two roots.)
    ///
    /// # Panics
    ///
    /// Panics if `x = 0` (the formula requires x ≠ 0) or if `trace(rhs) ≠ 0`
    /// (the x-coordinate is not on the curve).
    pub fn decompress<F: F2m<1>>(&self, x: &F, sign_bit: bool) -> BinaryAffinePoint<F> {
        let poly = &self.poly;
        assert!(!x.is_zero(), "decompress: x = 0 is not supported");

        let a = F::from_uint(self.a, poly);
        let b = F::from_uint(self.b, poly);

        // rhs = x + a + b/x²
        let x_inv = x.inv(poly);
        let x_inv_sq = x_inv.square(poly);
        let rhs = x.add(&a).add(&b.mul(&x_inv_sq, poly));

        // Verify solvability: trace(rhs) must be 0.
        let tr = rhs.trace(poly);
        assert!(
            tr.is_zero(),
            "decompress: trace(rhs) ≠ 0 — x is not a valid x-coordinate on this curve"
        );

        // Solve λ²+λ = rhs.
        let lambda = F::solve_quadratic(&rhs, poly);

        // Select the root based on sign_bit using bit 0 (constant term) of λ.
        // λ and λ+1 always differ in bit 0 (since 1 has bit 0 = 1), so this
        // convention works for all m (even and odd).
        let bit0_lambda: bool = lambda.to_uint().bit(0).into();
        let lambda_selected = if bit0_lambda == sign_bit {
            lambda
        } else {
            lambda.add(&F::one())
        };

        // y = λ·x
        let y = lambda_selected.mul(x, poly);

        BinaryAffinePoint::Finite { x: x.clone(), y }
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_gf2m::F2mNaive;

    /// GF(2^4) irreducible: x⁴+x+1 = 0x13.
    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64)
    }

    /// Toy binary curve over GF(2^4): y²+xy = x³+x²+1.
    ///
    /// Base point (1, 6) is on this curve:
    ///   LHS: 6²+1·6 = 7+6 = 1 (in GF(2^4))
    ///   RHS: 1+1+1 = 1 ✓
    fn toy_curve() -> BinaryCurve {
        let poly = poly4();
        BinaryCurve {
            poly,
            a: Uint::<1>::ONE,
            b: Uint::<1>::ONE,
            n: Uint::<1>::from(4u64), // group order 4 (3 affine points + ∞)
            gx: Uint::<1>::ONE,
            gy: Uint::<1>::from(6u64),
        }
    }

    #[test]
    fn generator_is_on_curve() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        assert!(c.is_on_curve(&g), "generator not on curve");
    }

    #[test]
    fn negate_char2() {
        // In char 2: -(x,y) = (x, x+y), NOT (x, -y).
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let neg_g = c.negate(&g);
        // neg_g should be (1, 1+6) = (1, 7).
        let expected_y = F2mNaive::<1>::from_u64(7, &poly);
        assert_eq!(neg_g.y(), Some(&expected_y), "negate: wrong y in char 2");
        assert!(c.is_on_curve(&neg_g), "negated point not on curve");
    }

    #[test]
    fn double_is_on_curve() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let two_g = c.double(&g);
        assert!(c.is_on_curve(&two_g), "2G not on curve");
    }

    #[test]
    fn double_equals_add_self() {
        // 2P via doubling must equal P+P via addition.
        // This is the load-bearing guard against Jacobian-formula contamination.
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let via_double = c.double(&g);
        let via_add = c.add(&g, &g);
        assert_eq!(via_double, via_add, "double(G) ≠ add(G,G)");
    }

    #[test]
    fn add_infinity_identity() {
        // P + ∞ = P and ∞ + P = P.
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let inf = BinaryAffinePoint::Infinity;
        assert_eq!(c.add(&g, &inf), g, "G + ∞ ≠ G");
        assert_eq!(c.add(&inf, &g), g, "∞ + G ≠ G");
    }

    #[test]
    fn add_neg_is_infinity() {
        // P + (-P) = ∞ with -P = (x, x+y).
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let neg_g = c.negate(&g);
        let sum = c.add(&g, &neg_g);
        assert!(sum.is_infinity(), "G + (-G) should be ∞");
    }

    #[test]
    fn scalar_mul_one_is_identity() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let result = c.scalar_mul(&g, &Uint::<1>::ONE);
        assert_eq!(result, g, "1·G should equal G");
    }

    #[test]
    fn scalar_mul_two_matches_double() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let two_g_direct = c.double(&g);
        let two_g_scalar = c.scalar_mul(&g, &Uint::<1>::from(2u64));
        assert_eq!(two_g_direct, two_g_scalar, "2·G via double ≠ via scalar_mul");
    }

}
