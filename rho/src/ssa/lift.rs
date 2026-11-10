//! Affine F_p point → Z_p lift via Hensel for the SSA anomalous-curve attack.
//!
//! Given an affine point `P = (x₀, y₀) ∈ E(F_p)`, a lift of the curve coefficient `a` and a target
//! precision `k`, this module lifts `P` to a pair `(x̃, ỹ) ∈ Z/p^k × Z/p^k` satisfying
//! `ỹ² ≡ x̃³ + ã·x̃ + b (mod p^k)`, where `ã` is the lifted coefficient. Every point that is lifted
//! or multiplied together must use the same `ã` (see [`lifted_a`]).
//!
//! # Lift strategy
//!
//! The x-coordinate is lifted exactly: `x̃ = Zp::new(x₀, p, k)`. The y-coordinate is
//! Hensel-solved from the curve equation: build `g(y) = y² − c` where `c = x̃³ + ã·x̃ + b`
//! (computed over `BigInt` mod `p^k`), then call `hensel_lift(g, y₀, p, k)` with `y₀` as
//! the simple root mod `p`.
//!
//! # Simple-root precondition
//!
//! The Hensel y-solve requires `g'(y₀) = 2·y₀ ≢ 0 mod p`, i.e. `y₀ ≠ 0` and `p ≠ 2`.
//! A 2-torsion base point (`y₀ = 0`) or `p = 2` is a degenerate case; `hensel_lift` returns
//! `HenselError::NonSimpleRoot` in that case, which propagates as `LiftError::Hensel`.
//!
//! # Lift-correctness invariant
//!
//! The returned `(x̃, ỹ)` satisfies `ỹ² ≡ x̃³ + ã·x̃ + b (mod p^k)`.
//!
//! # Scale limit
//!
//! The `Uint<4>` coordinates are read through their lowest limb (`as_words()[0]`), which is
//! exact only for `p < 2^64`; `ssa_solve` reports a larger `p` as `SsaError::Unsupported`
//! before any lift. Supporting a larger `p` would need the full `Uint<4>→BigInt` path.

use num_bigint::BigInt;
use shared_numfield::poly::IntPoly;
use shared_padic::{HenselError, Zp, ZpError, hensel_lift};

use shared_field::Fp;

use crate::curve::{AffinePoint, Curve};
use crate::ssa::uint4_to_bigint;

// ─── error type ──────────────────────────────────────────────────────────────

/// Errors from lifting an affine point to Z_p.
///
/// These belong to the lift helpers, not to the SSA solver: after the solver has validated its
/// input none of them can occur (a finite point on an odd-prime-order curve has `y₀ ≠ 0`), and
/// the solver maps any that does to an invariant failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LiftError {
    /// The point is the identity, which has no affine coordinates to lift.
    PointAtInfinity,
    /// The Hensel y-solve failed (for example a non-simple root at a 2-torsion point `y₀ = 0`).
    Hensel(HenselError),
    /// A Z/p^k arithmetic error from the p-adic layer.
    Padic(ZpError),
}

impl std::fmt::Display for LiftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LiftError::PointAtInfinity => write!(f, "the point at infinity cannot be lifted"),
            LiftError::Hensel(e) => write!(f, "Hensel lift failed: {e}"),
            LiftError::Padic(e) => write!(f, "p-adic arithmetic error: {e}"),
        }
    }
}

impl std::error::Error for LiftError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            LiftError::PointAtInfinity => None,
            LiftError::Hensel(e) => Some(e),
            LiftError::Padic(e) => Some(e),
        }
    }
}

impl From<HenselError> for LiftError {
    fn from(e: HenselError) -> Self {
        LiftError::Hensel(e)
    }
}

impl From<ZpError> for LiftError {
    fn from(e: ZpError) -> Self {
        LiftError::Padic(e)
    }
}

// ─── point lift ──────────────────────────────────────────────────────────────

/// The coefficient `a + p·shift` of the lift of `curve` that shifts `a` by a multiple of `p`.
///
/// The lifted curve `y² = x³ + (a + p·shift)x + b` reduces to `curve` modulo `p` for every
/// `shift`. Whether `p·G̃` lands in the kernel of reduction's first or a deeper filtration step
/// depends on this choice (see [`crate::ssa::reduce::kernel_filtration_level`]), so the point
/// lifts and the multiplication by `p` must all use the one lifted coefficient. `a` is taken as
/// stored, not reduced.
pub fn lifted_a(curve: &Curve, shift: u64) -> BigInt {
    uint4_to_bigint(&curve.a) + uint4_to_bigint(&curve.p) * BigInt::from(shift)
}

/// Lift an affine F_p point to Z_p coordinates at precision `k` on a lift of the curve.
///
/// Given `point = (x₀, y₀) ∈ E(F_p)` and the lifted curve `y² = x³ + a_lift·x + b`, returns
/// `(x̃, ỹ)` in `Z/p^k × Z/p^k` satisfying `ỹ² ≡ x̃³ + a_lift·x̃ + b (mod p^k)`.
///
/// # Arguments
///
/// - `point` — an affine point on `curve` over `F_p` (must be finite, not the identity).
/// - `curve` — the short Weierstrass curve `y² = x³ + ax + b` over `F_p`; supplies `p` and `b`.
/// - `a_lift` — the coefficient `a` of the lifted curve, congruent to `curve.a` modulo `p` (see
///   [`lifted_a`]; `curve.a` itself gives the verbatim lift).
/// - `k` — the target p-adic precision.
///
/// # Returns
///
/// `(x̃, ỹ)` as a `(Zp, Zp)` pair at precision `k`.
///
/// # Errors
///
/// - [`LiftError::PointAtInfinity`] if the point is the identity.
/// - [`LiftError::Hensel`] for a 2-torsion point (`y₀ = 0`) or `p = 2` (non-simple root — Hensel
///   lift not applicable).
/// - [`LiftError::Padic`] for Z/p^k arithmetic errors.
///
/// # Scale limit
///
/// The `as_words()[0]` extraction of the `Uint<4>` coordinates is exact only when `p < 2^64`
/// (see the module docs); a larger `p` would need the full `Uint<4>→BigInt` conversion.
pub fn lift_point<F: Fp<4>>(
    point: &AffinePoint<F>,
    curve: &Curve,
    a_lift: &BigInt,
    k: u32,
) -> Result<(Zp, Zp), LiftError> {
    // Extract (x₀, y₀) from the affine point. The identity cannot be lifted.
    let (x0_uint, y0_uint) = match point {
        AffinePoint::Infinity => return Err(LiftError::PointAtInfinity),
        AffinePoint::Finite { x, y } => (x.to_uint(), y.to_uint()),
    };

    // Single-limb conversion: exact only for p < 2^64 (see the module docs).
    let x0 = uint4_to_bigint(&x0_uint);
    let y0 = uint4_to_bigint(&y0_uint);
    let p_big = uint4_to_bigint(&curve.p);
    let b_big = uint4_to_bigint(&curve.b);

    // Step 1: lift x exactly — x̃ = Zp::new(x₀, p, k).
    // The x-coordinate needs no Hensel iteration: it lifts to itself.
    let x_tilde = Zp::new(&x0, &p_big, k)?;

    // Step 2: compute c = x₀³ + a·x₀ + b over BigInt (the RHS of the curve equation).
    // This is the constant term for the Hensel polynomial g(y) = y² − c.
    // We compute it over BigInt (not mod p^k) so that hensel_lift can reduce internally.
    let c = &x0 * &x0 * &x0 + a_lift * &x0 + &b_big;

    // Step 3: build g(y) = y² − c as an IntPoly.
    // Coefficients: [−c, 0, 1] (constant, linear, quadratic).
    // g(y) = 1·y² + 0·y + (−c)
    let g = IntPoly::from_coeffs(vec![-c.clone(), BigInt::from(0i64), BigInt::from(1i64)]);

    // Step 4: Hensel-solve g(y) = 0 for ỹ, starting from y₀ as the root mod p.
    // g'(y) = 2y; g'(y₀) = 2·y₀. The simple-root condition requires y₀ ≠ 0 and p ≠ 2.
    // If y₀ = 0 (2-torsion point), hensel_lift returns HenselError::NonSimpleRoot.
    let y_tilde_big = hensel_lift(&g, &y0, &p_big, k)?;
    let y_tilde = Zp::new(&y_tilde_big, &p_big, k)?;

    Ok((x_tilde, y_tilde))
}

// ─── lift-correctness check ───────────────────────────────────────────────────

/// Verify the lift-correctness invariant: `ỹ² ≡ x̃³ + a_lift·x̃ + b (mod p^k)`.
///
/// Returns `true` iff the lifted point satisfies the equation of the lifted curve (coefficient
/// `a_lift`, as passed to [`lift_point`]) mod `p^k`. This is the Hensel lift-correctness defense.
#[allow(dead_code)] // consumed by the SSA reduction and KATs
pub fn check_lift_on_curve(
    x_tilde: &Zp,
    y_tilde: &Zp,
    curve: &Curve,
    a_lift: &BigInt,
) -> Result<bool, ZpError> {
    let p_big = uint4_to_bigint(&curve.p);
    let b_big = uint4_to_bigint(&curve.b);
    let k = x_tilde.precision().min(y_tilde.precision());

    // lhs = ỹ²
    let y_sq = y_tilde.mul(y_tilde)?;

    // rhs = x̃³ + a·x̃ + b
    let x_sq = x_tilde.mul(x_tilde)?;
    let x_cu = x_sq.mul(x_tilde)?;
    let a_zp = Zp::new(a_lift, &p_big, k)?;
    let b_zp = Zp::new(&b_big, &p_big, k)?;
    let ax = a_zp.mul(x_tilde)?;
    let rhs = x_cu.add(&ax)?.add(&b_zp)?;

    Ok(y_sq.residue() == rhs.residue())
}
