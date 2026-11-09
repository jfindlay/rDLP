//! Verified outcome of a prime-order rho discrete-log solve.
//!
//! A solver returns [`DlogSolution`] only after checking `k·G == Q` with [`verify_log`], so a
//! success is independent of the walk that found it. [`DlogError`] keeps the failure classes
//! apart: invalid input, unsupported domain, resource exhaustion (which proves nothing about the
//! answer), and invariant failure (internally contradictory data, surfaced and never retried).
//! [`RhoEvidence`] records which path produced the answer and enough of the collision for a test to
//! recheck the walk's own identity from the fields alone.

use crypto_bigint::Uint;
use shared_field::Fp;
use shared_numth::is_prime;

use crate::curve::{AffinePoint, Curve};

/// A discrete logarithm `k < n` with `k·G == Q`, checked by the solver before returning.
///
/// `#[non_exhaustive]` keeps the fields readable everywhere while only `rho` can construct a
/// solution, so an unverified answer cannot be built by a caller.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct DlogSolution<E> {
    /// The logarithm: `k·G == Q`.
    pub k: u64,
    /// Which path produced `k` and the intermediate data behind it.
    pub evidence: E,
}

impl<E> DlogSolution<E> {
    /// Build a solution; crate-internal so every construction follows a verification.
    pub(crate) fn new(k: u64, evidence: E) -> Self {
        Self { k, evidence }
    }
}

/// Why a discrete-log solve did not return a verified logarithm.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DlogError {
    /// The instance does not satisfy the solver's preconditions.
    InvalidInput(InputDefect),
    /// The instance is outside the implementation's domain.
    Unsupported(&'static str),
    /// The attempt or step budget ran out. This is a resource outcome and proves nothing.
    BudgetExhausted {
        /// Attempts consumed.
        attempts: usize,
        /// Walk steps taken across all attempts.
        steps: u64,
    },
    /// A non-degenerate collision produced a `k` that fails `k·G == Q`.
    ///
    /// This is an invariant failure: it is surfaced, never retried.
    Inconsistent {
        /// The candidate that failed verification.
        k: u64,
    },
    /// Execution broke an invariant that carries no candidate (a walker panic, a lost channel).
    InvariantFailure(&'static str),
}

/// A precondition of the discrete-log instance that does not hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputDefect {
    /// The generator is the point at infinity.
    GeneratorAtInfinity,
    /// The generator does not satisfy the curve equation.
    GeneratorNotOnCurve,
    /// The target does not satisfy the curve equation.
    TargetNotOnCurve,
    /// The group order is zero. Reserved for composite-order checks; the prime-order check
    /// reports [`InputDefect::OrderNotPrime`] for `n ∈ {0, 1}`.
    ZeroOrder,
    /// The claimed group order is not prime.
    OrderNotPrime,
    /// `n·G` is not the point at infinity.
    WrongGeneratorOrder,
    /// `n·Q` is not the point at infinity, so `Q` cannot lie in `⟨G⟩`.
    TargetOutsideSubgroup,
}

/// Provenance of a rho solution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RhoEvidence {
    /// `Q = ∞` (`k = 0`) or `Q = G` (`k = 1`); no walk ran.
    Trivial,
    /// A walk collision that determined `k`.
    ///
    /// The identity `a1·G + b1·Q = (−1)^negated · (a2·G + b2·Q)  (mod n)` holds for the fields
    /// alone, so a test can recheck it without the solver. Endomorphism walks fold the orbit map
    /// into the tracked scalars before comparing points, so no endomorphism factor appears.
    Collision {
        /// Zero-based retry index; `0` for single-run solvers.
        attempt: usize,
        /// Walk steps counted by the call when the collision was accepted.
        steps: u64,
        /// Scalar of `G` at the first colliding state.
        a1: u64,
        /// Scalar of `Q` at the first colliding state.
        b1: u64,
        /// Scalar of `G` at the second colliding state.
        a2: u64,
        /// Scalar of `Q` at the second colliding state.
        b2: u64,
        /// Whether the second state equals the negation of the first.
        negated: bool,
    },
}

/// Check `k·G == Q` exactly.
///
/// Independent of the method that produced `k`; every solver calls it before constructing a
/// [`DlogSolution`].
#[must_use]
pub fn verify_log<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>, q: &AffinePoint<F>, k: u64) -> bool {
    curve.scalar_mul(g, &Uint::<4>::from(k)) == *q
}

/// Validate a prime-order discrete-log instance.
///
/// Checks, in order: `G ≠ ∞`; `G` and `Q` on the curve; `n` prime (deterministic Miller–Rabin over
/// all of `u64`); `n·G = ∞`; `n·Q = ∞`. It costs two scalar multiplications and one primality
/// test.
///
/// `n·Q = ∞` does not prove `Q ∈ ⟨G⟩` when `n²` divides the curve order: `Q` may then lie in
/// `E[n]` outside `⟨G⟩`. That residual case passes this check and surfaces from a solver as
/// [`DlogError::Inconsistent`].
///
/// # Errors
///
/// [`DlogError::InvalidInput`] with the first [`InputDefect`] found, in the order above.
/// [`InputDefect::ZeroOrder`] is never returned here.
pub fn check_prime_order_input<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Result<(), DlogError> {
    let defect = |d: InputDefect| Err(DlogError::InvalidInput(d));
    if g.is_infinity() {
        return defect(InputDefect::GeneratorAtInfinity);
    }
    if !curve.is_on_curve(g) {
        return defect(InputDefect::GeneratorNotOnCurve);
    }
    if !curve.is_on_curve(q) {
        return defect(InputDefect::TargetNotOnCurve);
    }
    let order = Uint::<4>::from(n);
    if !is_prime(&order) {
        return defect(InputDefect::OrderNotPrime);
    }
    if !curve.scalar_mul(g, &order).is_infinity() {
        return defect(InputDefect::WrongGeneratorOrder);
    }
    if !curve.scalar_mul(q, &order).is_infinity() {
        return defect(InputDefect::TargetOutsideSubgroup);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use shared_field::FpMonty4 as FpMonty;

    use crate::curve::test_curves::{composite_toy, tiny_a, TINY_A_N};

    fn input_defect(r: Result<(), DlogError>) -> InputDefect {
        match r {
            Err(DlogError::InvalidInput(d)) => d,
            other => panic!("expected InvalidInput, got {other:?}"),
        }
    }

    /// `(g, q = 5·g)` on tiny_a, a valid prime-order instance.
    fn tiny_a_instance() -> (Curve, AffinePoint<FpMonty>, AffinePoint<FpMonty>) {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(5u64));
        (curve, g, q)
    }

    /// The same x with `y + 1`: never on the curve for the tiny fixtures.
    fn bump_y(curve: &Curve, pt: &AffinePoint<FpMonty>) -> AffinePoint<FpMonty> {
        let AffinePoint::Finite { x, y } = pt else {
            panic!("bump_y: point at infinity");
        };
        let one = FpMonty::from_uint(Uint::<4>::ONE, &curve.p);
        AffinePoint::new(x.clone(), y.add(&one, &curve.p))
    }

    #[test]
    fn verify_log_accepts_k_and_rejects_neighbours() {
        let (curve, g, q) = tiny_a_instance();
        assert!(verify_log(&curve, &g, &q, 5));
        assert!(!verify_log(&curve, &g, &q, 6), "k + 1 must be rejected");
        assert!(!verify_log(&curve, &g, &q, 4), "k − 1 must be rejected");
    }

    #[test]
    fn valid_instance_passes_input_check() {
        let (curve, g, q) = tiny_a_instance();
        assert_eq!(check_prime_order_input(&curve, &g, &q, TINY_A_N), Ok(()));
    }

    #[test]
    fn input_check_rejects_infinite_generator() {
        let (curve, _, q) = tiny_a_instance();
        let inf = AffinePoint::<FpMonty>::Infinity;
        let r = check_prime_order_input(&curve, &inf, &q, TINY_A_N);
        assert_eq!(input_defect(r), InputDefect::GeneratorAtInfinity);
    }

    #[test]
    fn input_check_rejects_off_curve_generator() {
        let (curve, g, q) = tiny_a_instance();
        let bad = bump_y(&curve, &g);
        let r = check_prime_order_input(&curve, &bad, &q, TINY_A_N);
        assert_eq!(input_defect(r), InputDefect::GeneratorNotOnCurve);
    }

    #[test]
    fn input_check_rejects_off_curve_target() {
        let (curve, g, q) = tiny_a_instance();
        let bad = bump_y(&curve, &q);
        let r = check_prime_order_input(&curve, &g, &bad, TINY_A_N);
        assert_eq!(input_defect(r), InputDefect::TargetNotOnCurve);
    }

    #[test]
    fn input_check_rejects_composite_and_degenerate_orders() {
        let (curve, g, q) = tiny_a_instance();
        for n in [0u64, 1, 4, TINY_A_N + 1] {
            let r = check_prime_order_input(&curve, &g, &q, n);
            assert_eq!(input_defect(r), InputDefect::OrderNotPrime, "n = {n}");
        }
    }

    #[test]
    fn input_check_rejects_wrong_prime_order() {
        let (curve, g, q) = tiny_a_instance();
        // 13 is prime but is not the order of G on tiny_a.
        let r = check_prime_order_input(&curve, &g, &q, 13);
        assert_eq!(input_defect(r), InputDefect::WrongGeneratorOrder);
    }

    /// On the order-60 curve, `12·G` has order 5 while `G` itself has order 60 and so lies outside
    /// the order-5 subgroup.
    #[test]
    fn input_check_rejects_target_outside_subgroup() {
        let curve = composite_toy();
        let full: AffinePoint<FpMonty> = curve.generator();
        let g5 = curve.scalar_mul(&full, &Uint::<4>::from(12u64));
        assert!(curve.scalar_mul(&g5, &Uint::<4>::from(5u64)).is_infinity());
        assert!(!g5.is_infinity());

        let r = check_prime_order_input(&curve, &g5, &full, 5);
        assert_eq!(input_defect(r), InputDefect::TargetOutsideSubgroup);
        // The order-5 subgroup itself is a valid prime-order instance.
        let inside = curve.scalar_mul(&g5, &Uint::<4>::from(3u64));
        assert_eq!(check_prime_order_input(&curve, &g5, &inside, 5), Ok(()));
    }
}
