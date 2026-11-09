//! Smart–Satoh–Araki p-adic attack on anomalous elliptic curves.
//!
//! This module implements the SSA reduction, which solves the ECDLP in polynomial time on
//! anomalous curves (curves with trace of Frobenius = 1, i.e. `#E(F_p) = p`). It consumes
//! the p-adic substrate from `shared-padic` (p-adic arithmetic, Hensel lifting, and the
//! p-adic logarithm).
//!
//! # Structure
//!
//! - [`mod`] (this file) — `SsaError` enum, the hardcoded anomalous toy curve fixture, and
//!   `verify_anomalous` (exact O(p) point-count via Legendre symbol, any prime below 2^64).
//! - [`lift`] — affine F_p point → Z_p lift via Hensel for the SSA reduction.
//!
//! # Anomalous curve fixture
//!
//! The fixture is `y² = x³ + 5` over `F_7` (`p = 7`, `a = 0`, `b = 5`). This curve has
//! exactly 7 points (including the point at infinity), so `#E(F_7) = 7 = p` and the trace
//! of Frobenius is 1. The base point `G = (3, 2)` has `y ≠ 0`, ensuring the Hensel y-solve
//! has a simple root (`f'(y₀) = 2·y₀ = 4 ≢ 0 mod 7`).
//!
//! # Scale limits
//!
//! The fixture is hand-picked (not discovered via Schoof–SEA). The O(p) point-count verify
//! is a fixture check, not a general point-counting algorithm. Toy precision only. Primes or
//! coefficients that do not fit one 64-bit limb are reported as [`SsaError::Unsupported`], never
//! as a non-anomalous curve.

pub mod formal_log;
pub mod lift;
pub mod reduce;

pub use reduce::ssa_solve;

use crypto_bigint::Uint;
use num_bigint::BigInt;
use shared_padic::HenselError;
use shared_padic::ZpError;

use shared_field::Fp;

use crate::curve::Curve;

// ─── error type ──────────────────────────────────────────────────────────────

/// Errors from the SSA reduction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SsaError {
    /// The curve is not anomalous: an exact point count gave `#E(F_p) ≠ p`.
    ///
    /// The SSA reduction exists only for anomalous curves. This is a proven negative, reported
    /// only after [`verify_anomalous`] counted the points exactly; a curve whose size is outside
    /// the counting range is [`SsaError::Unsupported`], never `NotAnomalous`.
    NotAnomalous,
    /// The curve was proven anomalous and the points lifted, but neither candidate for the
    /// logarithm satisfied `k·G = Q`.
    ///
    /// This is not a statement about the curve, and it is surfaced, never retried. On a lift
    /// where `p·P̃` has valuation exactly 1 the quotient of formal logarithms determines `k`, so
    /// the failure indicates an invariant failure (internally contradictory data). On the
    /// canonical lift (for example an `a = 0` curve with complex multiplication, like the
    /// fixture) `p·P̃` lies deeper in the filtration and the quotient is a lift-dependent
    /// second-order artefact; the method's derivation does not hold there, so this error can
    /// occur on valid input. See the `reduce` module documentation for the candidates.
    Inconsistent {
        /// The quotient of formal logarithms modulo `p` that no candidate could be derived from.
        k_raw: u64,
    },
    /// The input lies outside the implementation's domain: a prime or curve coefficient that
    /// does not fit one 64-bit limb, or a prime below 5 (short Weierstrass form needs `p > 3`).
    Unsupported(&'static str),
    /// Hensel lifting failed (e.g. non-simple root — 2-torsion base point with y = 0).
    LiftFailed(HenselError),
    /// A Z/p^k arithmetic error from the p-adic layer.
    Padic(ZpError),
}

impl std::fmt::Display for SsaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SsaError::NotAnomalous => write!(f, "curve is not anomalous: #E(F_p) ≠ p"),
            SsaError::Inconsistent { k_raw } => write!(
                f,
                "anomalous curve but no candidate for k verified (log quotient {k_raw} mod p)"
            ),
            SsaError::Unsupported(why) => write!(f, "outside the supported domain: {why}"),
            SsaError::LiftFailed(e) => write!(f, "Hensel lift failed: {e}"),
            SsaError::Padic(e) => write!(f, "p-adic arithmetic error: {e}"),
        }
    }
}

impl std::error::Error for SsaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SsaError::LiftFailed(e) => Some(e),
            SsaError::Padic(e) => Some(e),
            SsaError::NotAnomalous | SsaError::Inconsistent { .. } | SsaError::Unsupported(_) => {
                None
            }
        }
    }
}

impl From<HenselError> for SsaError {
    fn from(e: HenselError) -> Self {
        SsaError::LiftFailed(e)
    }
}

impl From<ZpError> for SsaError {
    fn from(e: ZpError) -> Self {
        SsaError::Padic(e)
    }
}

// ─── anomalous toy curve fixture ─────────────────────────────────────────────

/// The prime for the anomalous toy fixture: p = 7.
pub const ANOMALOUS_TOY_P: u64 = 7;

/// Return the hardcoded anomalous toy curve fixture.
///
/// The curve is `y² = x³ + 5` over `F_7` (`p = 7`, `a = 0`, `b = 5`). It has exactly
/// `#E(F_7) = 7 = p` points (trace of Frobenius = 1). The base point `G = (3, 2)` has
/// `y ≠ 0`, so the Hensel y-solve has a simple root.
///
/// # Fixture limits
///
/// This fixture is hand-picked, not discovered via Schoof–SEA. The O(p) `verify_anomalous`
/// helper confirms the order; it is a fixture check, not a general point-counting algorithm.
pub fn anomalous_toy() -> Curve {
    Curve {
        p: Uint::<4>::from(ANOMALOUS_TOY_P),
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(5u64),
        // n = p = 7 (the group order equals the field characteristic for anomalous curves).
        n: Uint::<4>::from(ANOMALOUS_TOY_P),
        // Base point G = (3, 2): 3³ + 5 = 32 ≡ 4 mod 7, 2² = 4 ✓
        gx: Uint::<4>::from(3u64),
        gy: Uint::<4>::from(2u64),
    }
}

// ─── anomalous verification ───────────────────────────────────────────────────

/// Check whether `curve` is anomalous: `#E(F_p) = p` (trace of Frobenius = 1).
///
/// Counts points exactly via the Legendre symbol:
/// `#E = 1 + Σ_{x=0}^{p-1} (1 + (x³ + ax + b | p))`. Every intermediate product is reduced
/// modulo `p` in `u128` arithmetic, so the count is exact for every prime `p < 2^64`; the total
/// is held in `u128` as well, since `#E ≤ p + 1 + 2√p` can exceed `u64`. The cost is O(p) field
/// operations — this is a counting check, not Schoof–SEA, and is practical only for small `p`.
///
/// `p` is assumed prime, as for every [`Curve`]; the Legendre symbol is not meaningful otherwise.
///
/// # Errors
///
/// [`SsaError::Unsupported`] if `p`, `a` or `b` does not fit in 64 bits or `p < 5`. A curve
/// outside that range is never reported as non-anomalous.
pub fn verify_anomalous<F: Fp<4>>(curve: &Curve) -> Result<bool, SsaError> {
    let (p, a, b) = curve_params_u64(curve)?;

    // #E = 1 (infinity) + Σ_x (1 + Legendre(x³ + ax + b | p)); the summand is the number of
    // y with y² = x³ + ax + b, i.e. 0, 1 or 2.
    let p128 = u128::from(p);
    let mut count: u128 = 1;
    for x in 0..p128 {
        let x2 = x * x % p128;
        let x3 = x2 * x % p128;
        let rhs = (x3 + u128::from(a) * x % p128 + u128::from(b)) % p128;
        let leg = legendre_symbol(rhs as u64, p);
        count += (1 + leg) as u128;
    }

    Ok(count == p128)
}

/// Read `(p, a, b)` from `curve` as `u64`, with `a` and `b` reduced modulo `p`.
///
/// # Errors
///
/// [`SsaError::Unsupported`] if `p`, `a` or `b` has a non-zero limb above the lowest, or if
/// `p < 5` (outside short Weierstrass form, and the p-adic reduction needs an odd prime above 3).
pub(crate) fn curve_params_u64(curve: &Curve) -> Result<(u64, u64, u64), SsaError> {
    let single_limb = |u: &Uint<4>| {
        let w = u.as_words();
        (w[1] == 0 && w[2] == 0 && w[3] == 0).then_some(w[0])
    };
    let p = single_limb(&curve.p)
        .ok_or(SsaError::Unsupported("field prime does not fit in 64 bits"))?;
    if p < 5 {
        return Err(SsaError::Unsupported("field prime must exceed 3 for short Weierstrass form"));
    }
    let a = single_limb(&curve.a)
        .ok_or(SsaError::Unsupported("curve coefficient a does not fit in 64 bits"))?;
    let b = single_limb(&curve.b)
        .ok_or(SsaError::Unsupported("curve coefficient b does not fit in 64 bits"))?;
    Ok((p, a % p, b % p))
}

/// Compute the Legendre symbol (n | p) for odd prime p.
///
/// Returns 0 if p | n, 1 if n is a quadratic residue mod p, -1 otherwise.
/// Uses Euler's criterion: (n | p) = n^((p-1)/2) mod p.
fn legendre_symbol(n: u64, p: u64) -> i64 {
    if n % p == 0 {
        return 0;
    }
    let result = pow_mod(n, (p - 1) / 2, p);
    if result == 1 { 1 } else { -1 }
}

/// Modular exponentiation: base^exp mod modulus (u64, iterative square-and-multiply).
fn pow_mod(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result: u64 = 1;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exp >>= 1;
    }
    result
}

// ─── BigInt helpers (used by lift.rs) ────────────────────────────────────────

/// Extract the toy-scale u64 value from a `Uint<4>` field element.
///
/// # Panics (debug)
///
/// Panics in debug mode if the upper limbs are non-zero (p ≥ 2^64).
///
/// # Scale limits
///
/// SCALE: toy-scale only — crypto-scale p would need full Uint<4>→BigInt conversion.
pub(crate) fn uint4_to_u64(u: &Uint<4>) -> u64 {
    let words = u.as_words();
    debug_assert!(
        words[1] == 0 && words[2] == 0 && words[3] == 0,
        "uint4_to_u64: Uint<4> has non-zero upper limbs — p >= 2^64 not supported at toy scale \
         (single-limb limit)"
    );
    words[0]
}

/// Convert a `Uint<4>` field element to a `BigInt`.
///
/// # Scale limits
///
/// SCALE: toy-scale only — crypto-scale p would need the full limb vector.
pub(crate) fn uint4_to_bigint(u: &Uint<4>) -> BigInt {
    BigInt::from(uint4_to_u64(u))
}
