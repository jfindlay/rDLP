//! SSA reduction: multiply-by-p kernel step, F_p division, and `ssa_solve` entry point.
//!
//! This module implements the core of the Smart–Satoh–Araki reduction:
//!
//! 1. Verify the curve is anomalous.
//! 2. Lift the curve to Z_p by shifting `a` by a multiple of `p`, and lift G and Q to it via
//!    Hensel ([`crate::ssa::lift::lift_point`]).
//! 3. Multiply both lifts by p using projective-coordinate point addition over Z/p^k.
//! 4. Apply log(1 + t) at t = −X/Y to both results.
//! 5. Divide in F_p to recover k, then check `k·G = Q`.
//!
//! # Projective coordinates for multiply-by-p
//!
//! The p-th addition step in the multiply-by-p loop has a denominator divisible by p^2
//! (the affine addition formula breaks down). Standard projective coordinates `(X:Y:Z)`
//! with affine `(x, y) = (X/Z, Y/Z)` handle this correctly: the projective addition formula
//! works mod p^k even when Z ≡ 0 mod p (the point reduces to infinity mod p).
//!
//! # Formal group parameter from projective coordinates
//!
//! In standard projective coordinates, the formal group parameter is `t = −x/y = −X/Y`
//! (the Z factors cancel). After multiply-by-p, `v_p(Z) ≥ 1` and `v_p(t) ≥ 1`.
//!
//! # Which lift of the curve
//!
//! Let `E_n` be the points of the lifted curve with `v_p(t) ≥ n`. For a lift on which
//! `p·G̃ ∈ E_1 \ E_2` (`v_p(t) = 1`) the map `P ↦ ψ(p·P̃)/p mod p` is a homomorphism
//! `E(F_p) → F_p` independent of the lifted point (changing `P̃` by a point `κ` of `E_1` changes
//! `p·P̃` by `p·κ ∈ E_2`), nonzero on `G`, so `k = ψ(p·Q̃)/ψ(p·G̃) mod p` exactly. Whether
//! `p·G̃` lies in `E_2` depends on the lift of the *curve*: on the lifts for which it does (Smart's
//! canonical lift, which occurs with probability about `1/p` among lifts; Smart 1999, p. 195) the
//! quotient is a ratio of second-order terms with no meaning. The stored coefficients are not used
//! verbatim: both point lifts and the multiplication by `p` use the curve with `a` replaced by
//! `a + p·s`, `s ≢ 0 (mod p)`, which reduces to the same curve modulo `p`.
//!
//! If `p·G̃ ∈ E_2` on the first shift, exactly one retry is made with a second shift
//! `s₂ ≢ s₁ (mod p)`, `s₂ ≢ 0`; `p·G̃ ∈ E_2` on both is an invariant failure. The retry is
//! a second lift of the curve, not a second attempt with fresh randomness: the shifts are fixed
//! ([`FIRST_SHIFT`], [`SECOND_SHIFT`]) and the outcome does not depend on a random source.
//!
//! # Verification step
//!
//! The quotient `k_raw = ψ(p·Q̃) / ψ(p·G̃) mod p` is the answer, and `k_raw·G = Q` is still checked
//! in E(F_p) before anything is returned. A failed check cannot arise from valid input on a lift
//! with `p·G̃ ∈ E_1 \ E_2`, so it is reported as [`SsaError::Inconsistent`], an invariant failure.
//!
//! # Scale limits
//!
//! Toy scale only. The multiply-by-p loop is O(p), though double-and-add would make it
//! O(log p); the verification is O(log p) group operations.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use shared_padic::{Zp, ZpError};

use shared_field::Fp;

use crate::curve::{AffinePoint, Curve};
use crate::ssa::formal_log::{
    elliptic_log_proj, formal_parameter_valuation, padic_valuation_bigint, pow_bigint,
};
use crate::ssa::lift::{lift_point, lifted_a};
use crate::ssa::{SsaError, SsaInputDefect, check_curve, uint4_to_bigint, verify_anomalous};

// ─── parameters of the reduction ─────────────────────────────────────────────

/// The p-adic working precision `k`: all arithmetic is modulo `p^k`.
///
/// The formal parameter of `p·P̃` has valuation 1 on a usable lift, so its leading unit is
/// readable only once `t` is known modulo `p²`; `formal_log::elliptic_log_proj` keeps
/// `k − v_p(t)` digits of `t`, so `k ≥ 3` is the least precision that can work (Smart 1999 notes
/// that arithmetic modulo `p²` of the group suffices; the extra digit here is the one the
/// logarithm's parameter loses to its own valuation).
///
/// Measured, not derived: with `k = 3` (and `4`, `5`, `6`, `8`) every `k` in `1..p` was recovered
/// on each of the 481 anomalous curves `y² = x³ + ax + b` over `F_p`, `5 ≤ p ≤ 59`, `b ≠ 0`, with
/// the same degenerate shift class for every precision (exploratory scan, not retained as a
/// test), so precision `3` already suffices on that range. The value `4` keeps one digit of margin
/// for the projective multiplication by `p`, whose coordinates are not normalised and could lose
/// precision on inputs the scan did not reach.
///
/// SCALE: toy-scale only — crypto-scale would need `k` fixed independently of `p`, with the
/// multiplication by `p` done by double-and-add.
pub const WORKING_PRECISION: u32 = 4;

/// The shift `s₁` of the first lifted curve `a ↦ a + p·s₁`.
pub const FIRST_SHIFT: u64 = 1;

/// The shift `s₂` of the retry, used only if `p·G̃ ∈ E₂` on the first lift.
///
/// `s₂ ≢ s₁` and `s₂ ≢ 0 (mod p)` for every prime `p ≥ 5`, the smallest prime the reduction
/// accepts, so the retry is a different lift of the curve than the first.
pub const SECOND_SHIFT: u64 = 2;

// ─── SSA entry point ──────────────────────────────────────────────────────────

/// Solve the ECDLP `Q = k·G` on an anomalous curve via the Smart–Satoh–Araki reduction.
///
/// Implements the SSA reduction: lift the curve to Z_p with `a ↦ a + p·s`, lift G and Q to it,
/// multiply both by p to land in the kernel of reduction, apply log(1 + t) at t = −X/Y, and divide
/// in F_p.
///
/// # Arguments
///
/// - `curve` — the short Weierstrass curve `y² = x³ + ax + b` over F_p.
/// - `g` — the base point G ∈ E(F_p).
/// - `q` — the target point Q ∈ E(F_p) with `Q = k·G`.
/// - `n` — the group order (must equal p for anomalous curves).
///
/// # Returns
///
/// `k` such that `Q = k·G` in E(F_p), as a `u64`, checked by `k·G = Q` before returning.
/// `Q = ∞` returns `0` once the curve, `G` and `Q` have been validated, on any such curve.
///
/// # Errors
///
/// The first failing check is returned, in this order:
///
/// 1. [`SsaError::Unsupported`] if `p`, `a` or `b` does not fit in 64 bits or `p < 5`.
/// 2. [`SsaError::InvalidInput`] with [`SsaInputDefect::FieldNotPrime`] if `p` is composite.
/// 3. [`SsaInputDefect::SingularCurve`] if `4a³ + 27b² ≡ 0 (mod p)`.
/// 4. [`SsaInputDefect::GeneratorAtInfinity`] if `G = ∞`.
/// 5. [`SsaInputDefect::GeneratorNotOnCurve`] if `G` fails the curve equation.
/// 6. [`SsaInputDefect::TargetNotOnCurve`] if `Q` fails the curve equation.
/// 7. `Q = ∞` returns `Ok(0)`.
/// 8. [`SsaError::NotAnomalous`] if an exact count gives `#E(F_p) ≠ p`.
/// 9. [`SsaInputDefect::OrderNotFieldPrime`] if `n ≠ p`.
///
/// After these checks `E(F_p)` has prime order `p`, so every finite point on it generates the
/// group: `n·G = ∞`, `n·Q = ∞` and `Q ∈ ⟨G⟩` follow, and are not tested separately. Then the
/// reduction runs, and can fail only with:
///
/// - [`SsaError::Inconsistent`] if the quotient of the formal logarithms fails `k·G = Q`. On a
///   lift with `p·G̃ ∈ E₁ \ E₂` the quotient is exactly `k`, so this is a breach of the data.
/// - [`SsaError::InvariantFailure`] if a property that must hold on validated input is violated:
///   a lift or p-adic error, a point outside the kernel of reduction, a vanishing or non-unit
///   logarithm, or `p·G̃ ∈ E₂` on both shifted lifts (about `1/p` of the lifts have `p·G̃ ∈ E₂`,
///   and the two shifts are distinct lifts, so valid input is not expected to reach this).
///
/// # Scale limits
///
/// SCALE: toy-scale only. The multiply-by-p loop is O(p) BigInt work, and the working precision
/// [`WORKING_PRECISION`] is fixed for small toy primes.
pub fn ssa_solve<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Result<u64, SsaError> {
    // Validate the input in the documented order; the first failure is returned. The count is a
    // proven negative only on a nonsingular curve over a prime field, so the curve check precedes
    // it. A prime outside the counting range is `Unsupported`, not a negative.
    let (p, _, _) = check_curve(curve)?;
    if matches!(g, AffinePoint::Infinity) {
        return Err(SsaError::InvalidInput(SsaInputDefect::GeneratorAtInfinity));
    }
    if !curve.is_on_curve(g) {
        return Err(SsaError::InvalidInput(SsaInputDefect::GeneratorNotOnCurve));
    }
    if !curve.is_on_curve(q) {
        return Err(SsaError::InvalidInput(SsaInputDefect::TargetNotOnCurve));
    }
    // Q = ∞ is 0·G whatever the (valid) curve: the logarithm is 0, with no reduction to run.
    if matches!(q, AffinePoint::Infinity) {
        return Ok(0);
    }
    if !verify_anomalous::<F>(curve)? {
        return Err(SsaError::NotAnomalous);
    }
    if n != p {
        return Err(SsaError::InvalidInput(SsaInputDefect::OrderNotFieldPrime));
    }

    // One lift of the curve per attempt. If p·G̃ lands in E₂ the quotient of logarithms carries no
    // information; the lift is then the wrong one, and exactly one other lift is tried. Smart
    // (1999, p. 195) prescribes this: detect the failure, repeat with another lift of the curve.
    let k_raw = quotient_with_retry::<F>(curve, g, q, p, FIRST_SHIFT, SECOND_SHIFT)?;

    // The quotient is the answer, but `k·G = Q` is checked in E(F_p) before anything is returned,
    // so the result never rests on the logarithm alone.
    verify_candidate(k_raw, g, q, curve)
}

/// The raw quotient on the lift with shift `first`, or on the lift with shift `second` if
/// `p·G̃ ∈ E₂` on the first.
///
/// At most two lifts are made. `p·G̃ ∈ E₂` on both is an invariant failure.
fn quotient_with_retry<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    p: u64,
    first: u64,
    second: u64,
) -> Result<u64, SsaError> {
    match attempt::<F>(curve, g, q, p, first)? {
        Attempt::Quotient(k_raw) => Ok(k_raw),
        Attempt::PgInE2 => match attempt::<F>(curve, g, q, p, second)? {
            Attempt::Quotient(k_raw) => Ok(k_raw),
            Attempt::PgInE2 => Err(SsaError::InvariantFailure("p·G̃ in E₂ on two shifted lifts")),
        },
    }
}

/// The outcome of one lifted curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Attempt {
    /// `p·G̃ ∈ E₁ \ E₂`: the raw quotient `ψ(p·Q̃) / ψ(p·G̃) mod p`.
    Quotient(u64),
    /// `p·G̃ ∈ E₂`: this lift of the curve is the degenerate one.
    PgInE2,
}

/// The messages of the breaches that can occur lifting and multiplying one point.
#[derive(Clone, Copy)]
struct Breaches {
    /// The Hensel lift of the point failed.
    lift: &'static str,
    /// The p-adic layer rejected the coordinates of the multiple of the point.
    multiply: &'static str,
}

const G_BREACHES: Breaches = Breaches {
    lift: "lift of G failed on validated input",
    multiply: "p-adic layer rejected the coordinates of p·G̃",
};

const Q_BREACHES: Breaches = Breaches {
    lift: "lift of Q failed on validated input",
    multiply: "p-adic layer rejected the coordinates of p·Q̃",
};

/// Lift `point` to the curve with coefficient `a_lift` and multiply the lift by `p`.
///
/// Returns the projective `(X, Y)` of `p·P̃` at [`WORKING_PRECISION`].
///
/// On validated input the point is finite with `y₀ ≠ 0` (an odd prime order has no 2-torsion), so
/// the lift cannot fail; a failure is an invariant breach, mapped explicitly, never through
/// `From`.
fn lift_and_multiply_by_p<F: Fp<4>>(
    point: &AffinePoint<F>,
    curve: &Curve,
    a_lift: &BigInt,
    breaches: Breaches,
) -> Result<(Zp, Zp), SsaError> {
    let p_big = uint4_to_bigint(&curve.p);
    let p = crate::ssa::uint4_to_u64(&curve.p);
    let (x, y) = lift_point::<F>(point, curve, a_lift, WORKING_PRECISION)
        .map_err(|_| SsaError::InvariantFailure(breaches.lift))?;
    // The same `a_lift` is used for the multiplication by p as for the lift: the lifted point is
    // on that curve and on no other.
    multiply_by_p_proj(&x, &y, a_lift, &p_big, WORKING_PRECISION, p)
        .map_err(|_| SsaError::InvariantFailure(breaches.multiply))
}

/// Run the reduction on the lift of the curve with `a ↦ a + p·shift`.
fn attempt<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    p: u64,
    shift: u64,
) -> Result<Attempt, SsaError> {
    // The lifted curve is computed once and passed to every lift and multiplication below.
    let a_lift = lifted_a(curve, shift);
    let p_big = BigInt::from(p);

    let (xg, yg) = lift_and_multiply_by_p::<F>(g, curve, &a_lift, G_BREACHES)?;
    if formal_parameter_valuation(&xg, &yg)? >= 2 {
        return Ok(Attempt::PgInE2);
    }
    let (xq, yq) = lift_and_multiply_by_p::<F>(q, curve, &a_lift, Q_BREACHES)?;

    // t = −X/Y has v_p(t) ≥ 1 after multiply-by-p, so log(1 + t) converges. ψ(S) ≡ t(S) (mod p²)
    // for S in E₁, and so does log(1 + t); the leading units agree.
    let psi_g = elliptic_log_proj(&xg, &yg)?;
    let psi_q = elliptic_log_proj(&xq, &yq)?;

    // k ≡ ψ(p·Q̃) · ψ(p·G̃)⁻¹ (mod p): divide out the common p^v, then divide in F_p.
    Ok(Attempt::Quotient(fp_divide_logs(&psi_g, &psi_q, p, &p_big)?))
}

/// The step of the filtration of the kernel of reduction that `p·P̃` lies in, on a lift.
///
/// Lifts `point` to the curve with coefficient `a_lift` (see [`lifted_a`]) at
/// [`WORKING_PRECISION`], multiplies by `p`, and returns `v_p(t)` for the formal parameter
/// `t = −X/Y` of the result. The value 1 means `p·P̃ ∈ E₁ \ E₂`, so the reduction's quotient
/// is meaningful on that lift; a value of 2 or more means `p·P̃ ∈ E₂` and the lift is degenerate.
///
/// The curve and the point must be valid inputs of [`ssa_solve`] with `#E(F_p) = p`; this is a
/// diagnostic of a single lift, and checks none of that itself.
///
/// # Errors
///
/// [`SsaError::InvariantFailure`] if the lift or the multiplication fails, or the point is
/// outside the kernel of reduction.
pub fn kernel_filtration_level<F: Fp<4>>(
    curve: &Curve,
    point: &AffinePoint<F>,
    a_lift: &BigInt,
) -> Result<u32, SsaError> {
    let (x, y) = lift_and_multiply_by_p::<F>(point, curve, a_lift, G_BREACHES)?;
    formal_parameter_valuation(&x, &y)
}

// ─── multiply-by-p in projective coordinates ─────────────────────────────────

/// Multiply a lifted point by p using standard projective coordinates over Z/p^k.
///
/// Starting from the affine lift `(x̃, ỹ)` on the curve with coefficient `a`, computes
/// `p·(x̃, ỹ)` as a projective point `(X, Y)` (with Z factored out). Returns the projective X and
/// Y as `Zp` elements.
///
/// Uses standard projective coordinates `(X:Y:Z)` with affine `(x, y) = (X/Z, Y/Z)`.
/// The projective addition formula handles the degenerate p-th step where the affine
/// denominator is not a unit mod p^k.
///
/// # Scale limits
///
/// SCALE: toy-scale only. p additions is O(p) — correct only for small toy primes.
fn multiply_by_p_proj(
    x_tilde: &Zp,
    y_tilde: &Zp,
    a: &BigInt,
    p: &BigInt,
    k: u32,
    p_u64: u64,
) -> Result<(Zp, Zp), ZpError> {
    let pk = pow_bigint(p, k);

    // Starting point in standard projective: (X, Y, Z) = (x̃, ỹ, 1).
    let mut x_cur = x_tilde.residue().clone();
    let mut y_cur = y_tilde.residue().clone();
    let mut z_cur = BigInt::one();

    // The addend is the affine point (x̃, ỹ) = projective (x̃, ỹ, 1).
    let x_add = x_tilde.residue().clone();
    let y_add = y_tilde.residue().clone();
    let z_add = BigInt::one();

    // Perform p − 1 additions: result = p · (x̃, ỹ).
    for _ in 1..p_u64 {
        let (x_new, y_new, z_new) =
            proj_add_std(&x_cur, &y_cur, &z_cur, &x_add, &y_add, &z_add, a, &pk);
        x_cur = x_new;
        y_cur = y_new;
        z_cur = z_new;
    }

    // The result (x_cur, y_cur, z_cur) is p·(x̃, ỹ) in standard projective.
    // Return (X, Y) as Zp elements (Z is not needed for the formal group parameter t = −X/Y).
    let x_proj = Zp::new(&x_cur, p, k)?;
    let y_proj = Zp::new(&y_cur, p, k)?;

    Ok((x_proj, y_proj))
}

/// Standard projective point addition over Z/pk.
///
/// Computes `(X1:Y1:Z1) + (X2:Y2:Z2)` using the standard projective addition formulas.
/// Affine coordinates: `(x, y) = (X/Z, Y/Z)`.
///
/// Returns the result as `(X3, Y3, Z3)` in standard projective coordinates.
///
/// # Scale limits
///
/// SCALE: toy-scale only. BigInt arithmetic mod p^k is correct but not efficient.
fn proj_add_std(
    x1: &BigInt,
    y1: &BigInt,
    z1: &BigInt,
    x2: &BigInt,
    y2: &BigInt,
    z2: &BigInt,
    a: &BigInt,
    pk: &BigInt,
) -> (BigInt, BigInt, BigInt) {
    // Check if P1 = P2 (same point — use doubling formula).
    // In standard projective: P1 = P2 iff X1*Z2 = X2*Z1 and Y1*Z2 = Y2*Z1.
    let x1z2 = (x1 * z2).mod_floor(pk);
    let x2z1 = (x2 * z1).mod_floor(pk);
    let y1z2 = (y1 * z2).mod_floor(pk);
    let y2z1 = (y2 * z1).mod_floor(pk);

    if x1z2 == x2z1 {
        if y1z2 == y2z1 {
            // P1 = P2: use doubling.
            return proj_double_std(x1, y1, z1, a, pk);
        } else {
            // P1 = −P2: return infinity (0:1:0).
            return (BigInt::zero(), BigInt::one(), BigInt::zero());
        }
    }

    // General addition (P1 ≠ P2):
    // U1 = Y2·Z1, U2 = Y1·Z2
    // V1 = X2·Z1, V2 = X1·Z2
    // U = U1 − U2, V = V1 − V2
    // W = Z1·Z2
    // A = U²·W − V³ − 2·V²·V2
    // X3 = V·A
    // Y3 = U·(V²·V2 − A) − V³·U2
    // Z3 = V³·W
    let u1 = (y2 * z1).mod_floor(pk);
    let u2 = (y1 * z2).mod_floor(pk);
    let v1 = (x2 * z1).mod_floor(pk);
    let v2 = (x1 * z2).mod_floor(pk);
    let u = (&u1 - &u2).mod_floor(pk);
    let v = (&v1 - &v2).mod_floor(pk);
    let w = (z1 * z2).mod_floor(pk);

    let v2_sq = (&v * &v).mod_floor(pk);
    let v3 = (&v2_sq * &v).mod_floor(pk);
    let v2_v2 = (&v2_sq * &v2).mod_floor(pk);
    let u2_sq = (&u * &u).mod_floor(pk);
    let a_val = (&u2_sq * &w - &v3 - BigInt::from(2u64) * &v2_v2).mod_floor(pk);

    let x3 = (&v * &a_val).mod_floor(pk);
    let y3 = (&u * (&v2_v2 - &a_val) - &v3 * &u2).mod_floor(pk);
    let z3 = (&v3 * &w).mod_floor(pk);

    (x3, y3, z3)
}

/// Standard projective point doubling over Z/pk.
///
/// Computes `2·(X:Y:Z)` using the standard projective doubling formula.
fn proj_double_std(
    x: &BigInt,
    y: &BigInt,
    z: &BigInt,
    a: &BigInt,
    pk: &BigInt,
) -> (BigInt, BigInt, BigInt) {
    // W = a·Z² + 3·X²
    // S = Y·Z
    // B = X·Y·S
    // H = W² − 8·B
    // X3 = 2·H·S
    // Y3 = W·(4·B − H) − 8·Y²·S²
    // Z3 = 8·S³
    let z2 = (z * z).mod_floor(pk);
    let x2 = (x * x).mod_floor(pk);
    let w = (a * &z2 + BigInt::from(3u64) * &x2).mod_floor(pk);
    let s = (y * z).mod_floor(pk);
    let b = (x * y * &s).mod_floor(pk);
    let w2 = (&w * &w).mod_floor(pk);
    let h = (&w2 - BigInt::from(8u64) * &b).mod_floor(pk);

    let x3 = (BigInt::from(2u64) * &h * &s).mod_floor(pk);
    let y2 = (y * y).mod_floor(pk);
    let s2 = (&s * &s).mod_floor(pk);
    let y3 = (&w * (BigInt::from(4u64) * &b - &h) - BigInt::from(8u64) * &y2 * &s2)
        .mod_floor(pk);
    let z3 = (BigInt::from(8u64) * &s * &s2).mod_floor(pk);

    (x3, y3, z3)
}

// ─── F_p division of log values ──────────────────────────────────────────────

/// Divide two log(1 + t) values (at t = −X/Y) in F_p to recover k.
///
/// Both `psi_g` and `psi_q` are in `p^v·Z_p` for some `v ≥ 1`. Extracts the unit parts
/// (divides by `p^v`), then computes `(psi_q / p^v) · (psi_g / p^v)⁻¹ mod p`.
///
/// # Errors
///
/// [`SsaError::InvariantFailure`] if either log is zero at working precision (its valuation is at
/// least the precision, which is never read as valuation 0), or if the g-log is not a unit once
/// the common `p^v` is divided out (`v(ψ_G) > v(ψ_Q)`).
///
/// # Scale limits
///
/// SCALE: toy-scale only. The extraction of the unit part assumes small p.
fn fp_divide_logs(psi_g: &Zp, psi_q: &Zp, p_u64: u64, p: &BigInt) -> Result<u64, SsaError> {
    let g_res = psi_g.residue().clone();
    let q_res = psi_q.residue().clone();

    // `padic_valuation_bigint` returns 0 for zero, which would misread a vanishing log as a unit.
    if g_res.is_zero() {
        return Err(SsaError::InvariantFailure(
            "formal logarithm of p·G̃ vanishes at working precision",
        ));
    }
    if q_res.is_zero() {
        return Err(SsaError::InvariantFailure(
            "formal logarithm of p·Q̃ vanishes at working precision",
        ));
    }

    // Divide out the common power p^v, the smaller valuation.
    let vg = padic_valuation_bigint(&g_res, p);
    let vq = padic_valuation_bigint(&q_res, p);
    let v = vg.min(vq);

    let p_pow_v = pow_bigint(p, v);
    let g_unit = &g_res / &p_pow_v;
    let q_unit = &q_res / &p_pow_v;

    // Reduce mod p.
    let p_big = BigInt::from(p_u64);
    let g_unit_mod_p = g_unit.mod_floor(&p_big).to_u64_digits().1.first().copied().unwrap_or(0);
    let q_unit_mod_p = q_unit.mod_floor(&p_big).to_u64_digits().1.first().copied().unwrap_or(0);

    // The g-log must be a unit after the common power is removed; ψ_Q ≡ k·ψ_G forces
    // v(ψ_Q) ≥ v(ψ_G), so a non-unit means the reduction's data is contradictory.
    if g_unit_mod_p == 0 {
        return Err(SsaError::InvariantFailure(
            "formal logarithm of p·G̃ is not a unit after dividing out the common power of p",
        ));
    }

    let g_inv = mod_inv_u64(g_unit_mod_p, p_u64);
    Ok(mul_mod_u64(q_unit_mod_p, g_inv, p_u64))
}

// ─── candidate verification ──────────────────────────────────────────────────

/// Return `k_raw` if it satisfies `k·G = Q`.
///
/// One `scalar_mul` call: O(log p) group operations.
///
/// # Errors
///
/// [`SsaError::Inconsistent`] if `k_raw·G ≠ Q`. The quotient of the formal logarithms is the
/// answer on a lift with `p·G̃ ∈ E₁ \ E₂`, so a quotient that fails the check is a breach of the
/// data, surfaced and never retried.
fn verify_candidate<F: Fp<4>>(
    k_raw: u64,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    curve: &Curve,
) -> Result<u64, SsaError> {
    if &curve.scalar_mul(g, &crypto_bigint::Uint::<4>::from(k_raw)) == q {
        Ok(k_raw)
    } else {
        Err(SsaError::Inconsistent { k: k_raw })
    }
}

// ─── modular arithmetic helpers ───────────────────────────────────────────────

/// Compute `a · b mod p` without overflow for any `p < 2^64`.
fn mul_mod_u64(a: u64, b: u64, p: u64) -> u64 {
    (u128::from(a) * u128::from(b) % u128::from(p)) as u64
}

/// Compute the modular inverse of `a` mod `p` (prime) using Fermat's little theorem.
///
/// Requires `p` to be prime and `a ≠ 0 mod p`.
fn mod_inv_u64(a: u64, p: u64) -> u64 {
    // Fermat's little theorem: a^(p-2) mod p.
    let mut result: u64 = 1;
    let mut base: u64 = a % p;
    let mut exp: u64 = p - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % p as u128) as u64;
        }
        base = ((base as u128 * base as u128) % p as u128) as u64;
        exp >>= 1;
    }
    result
}


#[cfg(test)]
mod tests {
    use crypto_bigint::Uint;
    use shared_field::FpNaive4 as FpNaive;

    use super::*;
    use crate::ssa::{ANOMALOUS_TOY_P, anomalous_toy};

    fn toy_point(x: u64, y: u64) -> AffinePoint<FpNaive> {
        let p = Uint::<4>::from(ANOMALOUS_TOY_P);
        AffinePoint::Finite { x: FpNaive::from_u64(x, &p), y: FpNaive::from_u64(y, &p) }
    }

    fn times(curve: &Curve, g: &AffinePoint<FpNaive>, k: u64) -> AffinePoint<FpNaive> {
        curve.scalar_mul(g, &Uint::<4>::from(k))
    }

    /// A quotient satisfying `k·G = Q` is returned unchanged.
    #[test]
    fn verify_candidate_accepts_the_quotient_that_verifies() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);
        let q = times(&curve, &g, 3);

        assert_eq!(verify_candidate::<FpNaive>(3, &g, &q, &curve), Ok(3));
    }

    /// A quotient that fails `k·G = Q` is an invariant failure carrying that quotient: no other
    /// candidate (a half, a double) is substituted for it.
    #[test]
    fn verify_candidate_rejects_unverifiable_quotients() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);
        let q = times(&curve, &g, 3);

        // 6 = 2·3 and 5 = 3·2⁻¹ are the doubled and halved values; neither is accepted for k = 3.
        // 0 never verifies for Q ≠ ∞.
        for k_raw in [2, 0, 6, 5] {
            assert_eq!(
                verify_candidate::<FpNaive>(k_raw, &g, &q, &curve),
                Err(SsaError::Inconsistent { k: k_raw })
            );
        }
    }

    /// On the toy fixture the verbatim lift (`a = 0`, shift 0) has `p·G̃ ∈ E₂`, so the attempt
    /// reports the degenerate lift and no quotient.
    #[test]
    fn attempt_reports_degenerate_stored_lift() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);
        let q = times(&curve, &g, 3);

        assert_eq!(attempt::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P, 0), Ok(Attempt::PgInE2));
    }

    /// On a shifted lift the raw quotient of the formal logarithms is `k` itself, for every `k`.
    #[test]
    fn attempt_on_shifted_lift_returns_the_raw_quotient() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);

        for k in 1..ANOMALOUS_TOY_P {
            let q = times(&curve, &g, k);
            assert_eq!(
                attempt::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P, FIRST_SHIFT),
                Ok(Attempt::Quotient(k)),
                "k = {k}"
            );
        }
    }

    /// A degenerate first lift is followed by exactly one retry on the second shift.
    #[test]
    fn retry_runs_the_second_shift_after_a_degenerate_first() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);
        let q = times(&curve, &g, 3);

        // Shift 0 is the degenerate stored lift; shift 1 is not.
        assert_eq!(
            quotient_with_retry::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P, 0, 1),
            Ok(3)
        );
    }

    /// `p·G̃ ∈ E₂` on both lifts is an invariant failure, not a quotient.
    #[test]
    fn two_degenerate_lifts_are_an_invariant_failure() {
        let curve = anomalous_toy();
        let g = toy_point(3, 2);
        let q = times(&curve, &g, 3);

        // Shifts 0 and 7 give a = 0 and a = 49, the same curve modulo 7², so both are degenerate.
        assert_eq!(
            quotient_with_retry::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P, 0, ANOMALOUS_TOY_P),
            Err(SsaError::InvariantFailure("p·G̃ in E₂ on two shifted lifts"))
        );
    }

    /// `ssa_solve` recovers `k` when the first shifted lift is the degenerate one.
    ///
    /// Storing `a = 42 ≡ 0 (mod 7)` makes the first lift `a = 42 + 7 = 49 ≡ 0 (mod 7²)`, the
    /// degenerate class of the toy curve, so the second shift is used.
    #[test]
    fn ssa_solve_retries_when_the_first_shift_is_degenerate() {
        let mut curve = anomalous_toy();
        curve.a = Uint::<4>::from(42u64);
        let g = toy_point(3, 2);

        assert_eq!(
            attempt::<FpNaive>(&curve, &g, &g, ANOMALOUS_TOY_P, FIRST_SHIFT),
            Ok(Attempt::PgInE2),
            "the first shift of this stored coefficient must be degenerate"
        );
        for k in 1..ANOMALOUS_TOY_P {
            let q = times(&curve, &g, k);
            assert_eq!(
                ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P),
                Ok(k),
                "k = {k}"
            );
        }
    }

    fn zp(residue: i64) -> Zp {
        Zp::new(&BigInt::from(residue), &BigInt::from(ANOMALOUS_TOY_P), 8).expect("valid Zp")
    }

    /// A g-log of higher valuation than the q-log is not a unit once the common power is removed:
    /// an invariant failure.
    #[test]
    fn fp_divide_logs_rejects_non_unit_g_log() {
        // v(ψ_G) = 2 > v(ψ_Q) = 1.
        let result = fp_divide_logs(&zp(49), &zp(7), ANOMALOUS_TOY_P, &BigInt::from(7));
        assert!(matches!(result, Err(SsaError::InvariantFailure(_))), "got {result:?}");
    }

    /// A zero log, on either side, is a breach and is not read as valuation 0.
    #[test]
    fn fp_divide_logs_rejects_zero_logs() {
        for (g, q) in [(0, 7), (7, 0), (0, 0)] {
            let result = fp_divide_logs(&zp(g), &zp(q), ANOMALOUS_TOY_P, &BigInt::from(7));
            assert!(matches!(result, Err(SsaError::InvariantFailure(_))), "got {result:?}");
        }
    }

    /// The common power of p is divided out before the F_p quotient is taken.
    #[test]
    fn fp_divide_logs_divides_out_common_power() {
        // (3·7²) / (2·7²) = 3 · 2⁻¹ = 3 · 4 = 12 ≡ 5 (mod 7), and the same with a unit pair.
        let result = fp_divide_logs(&zp(2 * 49), &zp(3 * 49), ANOMALOUS_TOY_P, &BigInt::from(7));
        assert_eq!(result, Ok(5));
        let result = fp_divide_logs(&zp(2), &zp(3), ANOMALOUS_TOY_P, &BigInt::from(7));
        assert_eq!(result, Ok(5));
    }

    /// Multiplication modulo a 64-bit modulus does not overflow.
    #[test]
    fn mul_mod_u64_handles_full_width_operands() {
        let p = u64::MAX;
        assert_eq!(mul_mod_u64(p - 1, p - 1, p), 1);
        assert_eq!(mul_mod_u64(p - 1, 2, p), p - 2);
    }
}
