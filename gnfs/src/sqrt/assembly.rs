//! Congruence assembly and the checked factor driver for GNFS.
//!
//! A checked dependency has an exact rational root `x` with `x² = P` and, when `γ` is a square, a
//! verified algebraic root `β` whose modular image `y = β(c·m) mod N` satisfies `y² ≡ P (mod N)`.
//! Both squares equal `P` modulo `N`, so `x² ≡ y² (mod N)` and `N | (x − y)(x + y)`. This module
//! turns that congruence into a divisor of `N`:
//!
//! - [`factor_from_congruence`] checks `N > 1` and the congruence itself, reports `x ≡ ±y` as
//!   trivial, and otherwise takes `gcd(x − y, N)`, which is then a proper divisor.
//! - [`attempt_factor`] runs one checked dependency through the rational root, the verified
//!   algebraic root, and the congruence, and reports the outcome as a typed [`AttemptOutcome`].
//! - [`factor`] retries a bounded list of candidates in order and returns a [`FactorReport`].
//!
//! # What a report means
//!
//! A [`VerifiedFactor`] exists only when every step succeeded, so a returned divisor rests on a
//! verified algebraic root and an exact rational root over the same product `P`. A divisor found
//! as a by-product of the setup (a leading coefficient or a root denominator sharing a factor with
//! `N`) is an [`IncidentalFactor`] and is reported separately. A dependency whose congruence is
//! `x ≡ ±y (mod N)` is a [`AttemptOutcome::TrivialCongruence`], which is an expected outcome and
//! not an error. A divisor of `N` that appears only in that case, such as `gcd(x + y, N) = 5` for
//! `N = 35` and `x = y = 5`, is deliberately not reported (see [`factor_from_congruence`]).
//!
//! The algebraic target is certified a nonsquare only by the certificate that
//! [`algebraic_sqrt`] returns. An invalid candidate, a violated invariant, or an exhausted resource
//! budget is never recorded as a nonsquare, and an invalid candidate aborts the search instead of
//! being retried.

use std::fmt;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;
use shared_bigint::gcd;

use crate::linalg::{KernelSearchStatus, KernelVector};
use crate::sqrt::algebraic::algebraic_sqrt;
use crate::sqrt::context::{
    AlgebraicSqrtOutcome, FailureClass, IncidentalFactor, InvariantViolation, NonSquareCertificate,
    ResourceExhaustion, SqrtError, SqrtLimits, UnsupportedDomain, VerifiedRoot,
};
use crate::sqrt::dependency::{CheckedDependency, DependencyContext, DependencyError};
use crate::sqrt::rational::{RationalRoot, RationalSqrtOutcome, rational_sqrt};

// ─── factor_from_congruence ───────────────────────────────────────────────────

/// The result of extracting a divisor from a valid congruence of squares.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CongruenceOutcome {
    /// A proper divisor `1 < divisor < N` and its cofactor, with `divisor · cofactor = N`.
    Proper {
        /// The divisor `gcd(|x − y|, N)`.
        divisor: BigInt,
        /// `N / divisor`.
        cofactor: BigInt,
    },
    /// `x ≡ y` or `x ≡ −y (mod N)`: the congruence gives no divisor.
    Trivial,
}

/// Why a congruence of squares could not be assembled.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CongruenceError {
    /// The modulus is not greater than one.
    ModulusNotGreaterThanOne {
        /// The offending modulus.
        n: BigInt,
    },
    /// `x² ≢ y² (mod N)`: the inputs are not a congruence of squares.
    NotCongruent,
}

impl fmt::Display for CongruenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ModulusNotGreaterThanOne { n } => write!(f, "modulus must exceed 1, got {n}"),
            Self::NotCongruent => write!(f, "x^2 and y^2 are not congruent modulo N"),
        }
    }
}

impl std::error::Error for CongruenceError {}

/// Extract a divisor of `n` from `x² ≡ y² (mod n)`.
///
/// Requires `n > 1` and the congruence. If `x ≡ y` or `x ≡ −y (mod n)` the congruence is
/// [`CongruenceOutcome::Trivial`]. Otherwise `n` divides `(x − y)(x + y)` but neither factor, so
/// `g = gcd(|x − y|, n)` is neither `1` (that would force `n | x + y`) nor `n` (that would force
/// `n | x − y`), and `g` is returned as the divisor. A divisor that arises only because `x ≡ ±y`
/// and `x` shares a factor with `n` (or because `n` is even and `x + y = 2x`) is deliberately not
/// reported: it says nothing about the congruence.
///
/// # Errors
///
/// [`CongruenceError::ModulusNotGreaterThanOne`] if `n <= 1`, and
/// [`CongruenceError::NotCongruent`] if `n` does not divide `x² − y²`.
pub fn factor_from_congruence(
    x: &BigInt,
    y: &BigInt,
    n: &BigInt,
) -> Result<CongruenceOutcome, CongruenceError> {
    if *n <= BigInt::one() {
        return Err(CongruenceError::ModulusNotGreaterThanOne { n: n.clone() });
    }
    if !(x * x - y * y).is_multiple_of(n) {
        return Err(CongruenceError::NotCongruent);
    }
    if (x - y).is_multiple_of(n) || (x + y).is_multiple_of(n) {
        return Ok(CongruenceOutcome::Trivial);
    }
    let divisor = gcd(&(x - y), n);
    let cofactor = n / &divisor;
    Ok(CongruenceOutcome::Proper { divisor, cofactor })
}

// ─── verified factor ──────────────────────────────────────────────────────────

/// A proper divisor of `N` backed by an exact rational root and a verified algebraic root.
///
/// The fields are private and the only constructor is this module's attempt routine, so a value
/// of this type records that `x² = P` exactly, `y² ≡ P (mod N)` for a verified `β`, both over one
/// product `P`, and that `divisor · cofactor = N` with `1 < divisor < N`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedFactor {
    divisor: BigInt,
    cofactor: BigInt,
    rational: RationalRoot,
    algebraic: VerifiedRoot,
    relations: Vec<usize>,
}

impl VerifiedFactor {
    /// The proper divisor of `N`.
    #[must_use]
    pub fn divisor(&self) -> &BigInt {
        &self.divisor
    }

    /// `N / divisor`.
    #[must_use]
    pub fn cofactor(&self) -> &BigInt {
        &self.cofactor
    }

    /// The exact rational root `x₀` and its image.
    #[must_use]
    pub fn rational_root(&self) -> &RationalRoot {
        &self.rational
    }

    /// The verified algebraic root and its modular image `y`.
    #[must_use]
    pub fn algebraic_root(&self) -> &VerifiedRoot {
        &self.algebraic
    }

    /// The lifted original relation indices of the dependency, ascending.
    #[must_use]
    pub fn relations(&self) -> &[usize] {
        &self.relations
    }
}

// ─── attempt ─────────────────────────────────────────────────────────────────

/// The typed result of running one checked dependency to a congruence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttemptOutcome {
    /// A verified proper factor.
    Factor(VerifiedFactor),
    /// The congruence holds but is trivial: `x ≡ ±y (mod N)`.
    ///
    /// A divisor that arises only in this case (for example `gcd(x + y, N)` when `x ≡ y` and `x`
    /// shares a factor with `N`) is deliberately not reported, because it says nothing about the
    /// congruence; see [`factor_from_congruence`].
    TrivialCongruence {
        /// The rational image `x`.
        x: BigInt,
        /// The algebraic image `y`.
        y: BigInt,
    },
    /// The rational product is not a positive exact square. Unreachable for a validated
    /// dependency, and typed so that it is never conflated with an algebraic nonsquare.
    RationalNonSquare {
        /// The rational product.
        product: BigInt,
    },
    /// The algebraic product is certified not to be a square. This is the certificate returned by
    /// the reconstruction, passed through unchanged.
    AlgebraicNonSquare(NonSquareCertificate),
    /// The domain or the map is unsupported for this dependency.
    Unsupported(UnsupportedDomain),
    /// A proper divisor of `N` appeared as a by-product of the setup, not of a congruence.
    Incidental(IncidentalFactor),
    /// A resource budget ran out, in chronological order; never evidence of a nonsquare.
    ResourceExhausted(Vec<ResourceExhaustion>),
}

impl AttemptOutcome {
    /// The failure class, or `None` for a verified factor.
    #[must_use]
    pub fn class(&self) -> Option<FailureClass> {
        match self {
            Self::Factor(_) => None,
            Self::TrivialCongruence { .. } => Some(FailureClass::TrivialCongruence),
            Self::RationalNonSquare { .. } => Some(FailureClass::RationalNonSquare),
            Self::AlgebraicNonSquare(_) => Some(FailureClass::AlgebraicNonSquare),
            Self::Unsupported(_) => Some(FailureClass::Unsupported),
            Self::Incidental(_) => Some(FailureClass::IncidentalFactor),
            Self::ResourceExhausted(_) => Some(FailureClass::ResourceExhausted),
        }
    }
}

/// The outcome of one candidate, with its position and lifted relations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttemptRecord {
    /// The candidate's position in the supplied slice.
    pub candidate: usize,
    /// The lifted original relation indices of the candidate's dependency, ascending.
    pub relations: Vec<usize>,
    /// What happened.
    pub outcome: AttemptOutcome,
}

/// Errors of the factor driver, as opposed to its mathematical outcomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FactorError {
    /// A candidate did not validate as a dependency of the context. Computed-kernel candidates are
    /// kernel vectors by construction, so this is an upstream defect and is never retried.
    InvalidCandidate {
        /// The candidate's position in the supplied slice.
        candidate: usize,
        /// Why validation failed.
        error: DependencyError,
    },
    /// The algebraic square root reported invalid input or a violated invariant.
    Sqrt(SqrtError),
    /// An identity that the verified roots guarantee failed.
    Invariant(InvariantViolation),
}

impl fmt::Display for FactorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCandidate { candidate, error } => {
                write!(f, "candidate {candidate} is not a valid dependency: {error}")
            }
            Self::Sqrt(error) => write!(f, "{error}"),
            Self::Invariant(violation) => write!(f, "{violation}"),
        }
    }
}

impl std::error::Error for FactorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidCandidate { error, .. } => Some(error),
            Self::Sqrt(error) => Some(error),
            Self::Invariant(violation) => Some(violation),
        }
    }
}

/// Run one checked dependency through the rational root, the verified algebraic root, and the
/// congruence.
///
/// The steps, each stopping at its first non-success:
///
/// 1. [`rational_sqrt`]; a non-square product is [`AttemptOutcome::RationalNonSquare`].
/// 2. [`algebraic_sqrt`]; every outcome but a verified root passes through unchanged into the
///    matching variant, so a nonsquare is reported only with the certificate the reconstruction
///    returned.
/// 3. The verified root must be over the same product and context as the rational root.
/// 4. [`factor_from_congruence`] on the two modular images: a trivial result is
///    [`AttemptOutcome::TrivialCongruence`] and a proper one mints the [`VerifiedFactor`].
///
/// # Errors
///
/// [`FactorError::Sqrt`] if the reconstruction reports an error, and [`FactorError::Invariant`]
/// if the two roots disagree on `P` or the context, or the images are not a congruence of squares
/// modulo `N`. Each is a defect, never a nonsquare.
pub fn attempt_factor(
    dependency: &CheckedDependency,
    limits: SqrtLimits,
) -> Result<AttemptOutcome, FactorError> {
    let rational = match rational_sqrt(dependency) {
        RationalSqrtOutcome::Root(root) => root,
        RationalSqrtOutcome::NonSquare { product } => {
            return Ok(AttemptOutcome::RationalNonSquare { product });
        }
    };
    let algebraic = algebraic_sqrt(dependency, limits).map_err(FactorError::Sqrt)?;
    conclude(dependency, rational, algebraic)
}

/// Steps 3 and 4 of [`attempt_factor`], given both square roots.
fn conclude(
    dependency: &CheckedDependency,
    rational: RationalRoot,
    algebraic: AlgebraicSqrtOutcome,
) -> Result<AttemptOutcome, FactorError> {
    let root = match algebraic {
        AlgebraicSqrtOutcome::Verified(root) => root,
        AlgebraicSqrtOutcome::NotSquare(certificate) => {
            return Ok(AttemptOutcome::AlgebraicNonSquare(certificate));
        }
        AlgebraicSqrtOutcome::Unsupported(domain) => return Ok(AttemptOutcome::Unsupported(domain)),
        AlgebraicSqrtOutcome::Incidental(factor) => return Ok(AttemptOutcome::Incidental(factor)),
        AlgebraicSqrtOutcome::ResourceExhausted(spent) => {
            return Ok(AttemptOutcome::ResourceExhausted(spent));
        }
    };
    if root.rational_product() != rational.product() {
        return Err(invariant("rational and algebraic roots use different products"));
    }
    if root.context() != dependency.context() {
        return Err(invariant("verified root belongs to a different context"));
    }
    let x = rational.modular_image().clone();
    let y = root.modular_image().clone();
    match factor_from_congruence(&x, &y, dependency.context().modulus()) {
        Err(_) => Err(invariant("verified roots do not give a congruence of squares")),
        Ok(CongruenceOutcome::Trivial) => Ok(AttemptOutcome::TrivialCongruence { x, y }),
        Ok(CongruenceOutcome::Proper { divisor, cofactor }) => {
            Ok(AttemptOutcome::Factor(VerifiedFactor {
                divisor,
                cofactor,
                rational,
                algebraic: root,
                relations: dependency.relations().to_vec(),
            }))
        }
    }
}

fn invariant(what: &'static str) -> FactorError {
    FactorError::Invariant(InvariantViolation { what })
}

// ─── factor report ───────────────────────────────────────────────────────────

/// The overall result of trying a list of candidates.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FactorOutcome {
    /// The first verified proper factor.
    Factor(VerifiedFactor),
    /// No attempt factored, and the first incidental divisor seen is reported. Pipeline evidence
    /// never rests on this variant.
    IncidentalOnly(IncidentalFactor),
    /// No attempt produced a divisor. A complete negative result for the supplied candidates only
    /// when [`FactorReport::resource_limited`] is false; it says nothing about the kernel.
    Exhausted,
}

/// The outcome of [`factor`], with every attempt and the status of the candidate search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FactorReport {
    outcome: FactorOutcome,
    attempts: Vec<AttemptRecord>,
    search: KernelSearchStatus,
}

impl FactorReport {
    /// The overall outcome.
    #[must_use]
    pub fn outcome(&self) -> &FactorOutcome {
        &self.outcome
    }

    /// Every attempt made, in candidate order, including the successful one.
    #[must_use]
    pub fn attempts(&self) -> &[AttemptRecord] {
        &self.attempts
    }

    /// The status of the search that produced the candidates, as supplied by the caller.
    #[must_use]
    pub fn search(&self) -> &KernelSearchStatus {
        &self.search
    }

    /// True iff the candidate search was cut short or any attempt ran out of resources.
    ///
    /// When this is false and the outcome is [`FactorOutcome::Exhausted`], every supplied candidate
    /// was fully examined.
    #[must_use]
    pub fn resource_limited(&self) -> bool {
        !self.search.is_complete()
            || self
                .attempts
                .iter()
                .any(|a| matches!(a.outcome, AttemptOutcome::ResourceExhausted(_)))
    }
}

/// Try `candidates` in order against `context` and report what happened.
///
/// Each candidate is validated as a dependency of the context, then attempted with
/// [`attempt_factor`]. The loop stops at the first verified factor. An incidental divisor is
/// recorded and the loop continues, so [`FactorOutcome::IncidentalOnly`] is returned only when no
/// attempt factored. `search` is the status of the search that produced `candidates`; it is
/// carried into the report unchanged.
///
/// # Errors
///
/// [`FactorError::InvalidCandidate`] at the first candidate that does not validate: computed-kernel
/// candidates are kernel vectors by construction, so this is an upstream defect and the search is
/// aborted rather than retried. [`FactorError::Sqrt`] and [`FactorError::Invariant`] as for
/// [`attempt_factor`].
pub fn factor(
    context: &DependencyContext,
    candidates: &[KernelVector],
    search: KernelSearchStatus,
    limits: SqrtLimits,
) -> Result<FactorReport, FactorError> {
    let mut attempts = Vec::new();
    let mut incidental: Option<IncidentalFactor> = None;
    for (position, candidate) in candidates.iter().enumerate() {
        let dependency = context
            .validate_kernel_vector(candidate)
            .map_err(|error| FactorError::InvalidCandidate { candidate: position, error })?;
        let outcome = attempt_factor(&dependency, limits)?;
        let record = AttemptRecord {
            candidate: position,
            relations: dependency.relations().to_vec(),
            outcome,
        };
        match &record.outcome {
            AttemptOutcome::Factor(found) => {
                let found = found.clone();
                attempts.push(record);
                return Ok(FactorReport { outcome: FactorOutcome::Factor(found), attempts, search });
            }
            AttemptOutcome::Incidental(found) => {
                incidental.get_or_insert_with(|| found.clone());
            }
            _ => {}
        }
        attempts.push(record);
    }
    let outcome = incidental.map_or(FactorOutcome::Exhausted, FactorOutcome::IncidentalOnly);
    Ok(FactorReport { outcome, attempts, search })
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use shared_numfield::IntPoly;

    use super::*;
    use crate::filter::build_checked_matrix;
    use crate::linalg::RootCharacters;
    use crate::sieve::{
        CheckedRelations, ExponentVector, FactorBase, Relation, RelationChecker,
    };
    use crate::sqrt::context::SqrtContext;
    use crate::sqrt::dependency::DependencyContext;
    use crate::sqrt::rational::root_of_product;
    use crate::sqrt::reconstruct_root_from_factors;

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn limits() -> SqrtLimits {
        use crate::sqrt::context::{
            CrtLimits, DirectLimits, PrecisionLimits, PrimeSearchLimits, RecombinationLimits,
        };
        use shared_numfield::{NonresidueSearchLimits, RootSearchLimits};
        SqrtLimits {
            precision: PrecisionLimits { max_modulus_bits: 4096 },
            direct: DirectLimits {
                prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 100 },
                nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
            },
            crt: CrtLimits {
                prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 100 },
                max_primes: 1,
                root_search: RootSearchLimits { max_evaluations: 100_000 },
                nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
            },
            recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 20 },
        }
    }

    /// x² − 2 at m = 23 modulo 527 with the single relation (3, −2): `a − 23b = 49 = 7²`,
    /// `γ = 3 + 2α = (1 + α)²`, `β = ±(1 + α)`, `y = ±24` and `24² = 576 = 527 + 49`. Its
    /// factor base `FB(7, 7)` gives the rational primes [2, 3, 5, 7], so the entry for 7 is
    /// index 3 with exponent 2, and the algebraic norm `9 − 8 = 1` has no entry.
    fn seam_dependency() -> CheckedDependency {
        let f = IntPoly::from_coeffs(vec![big(-2), big(0), big(1)]);
        let ctx = SqrtContext::new(&f, &big(23), &big(527)).unwrap();
        let fb = FactorBase::new(&f, 7, 7);
        let checker = RelationChecker::new(&ctx, &fb).unwrap();
        let relation = Relation {
            a: big(3),
            b: big(-2),
            rational_exponents: ExponentVector { entries: vec![(3, 2)] },
            algebraic_exponents: ExponentVector { entries: vec![] },
            rational_sign: false,
        };
        let checked = CheckedRelations::new(&checker, &[relation]).unwrap();
        let (matrix, layout) = build_checked_matrix(&checked, 0).unwrap();
        let characters = RootCharacters::new(checked.checker(), &[]).unwrap();
        let context = DependencyContext::new(checked, matrix, layout, characters).unwrap();
        context.validate(&[0]).unwrap()
    }

    fn verified_root(dependency: &CheckedDependency) -> VerifiedRoot {
        match algebraic_sqrt(dependency, limits()).unwrap() {
            AlgebraicSqrtOutcome::Verified(root) => root,
            other => panic!("expected a verified root, got {other:?}"),
        }
    }

    fn rational_root(product: i64, context: &SqrtContext) -> RationalRoot {
        match root_of_product(&big(product), context) {
            RationalSqrtOutcome::Root(root) => root,
            other => panic!("expected a root, got {other:?}"),
        }
    }

    #[test]
    fn roots_over_different_products_are_an_invariant_violation() {
        // The verified root has P = 49; the rational root is of 64 = 8².
        let dependency = seam_dependency();
        let algebraic = AlgebraicSqrtOutcome::Verified(verified_root(&dependency));
        let rational = rational_root(64, dependency.context());
        assert_eq!(
            conclude(&dependency, rational, algebraic),
            Err(invariant("rational and algebraic roots use different products"))
        );
    }

    #[test]
    fn root_from_another_context_is_an_invariant_violation() {
        // x² − 2 at m = 3 modulo 7 with the factors (9, 0): P = 9, β = ±3. The rational root has
        // the same product, so only the context check can reject the pair.
        let f = IntPoly::from_coeffs(vec![big(-2), big(0), big(1)]);
        let other = SqrtContext::new(&f, &big(3), &big(7)).unwrap();
        let factors =
            crate::sqrt::RelationFactors::new(&other, &[(big(9), big(0))]).unwrap();
        let AlgebraicSqrtOutcome::Verified(foreign) =
            reconstruct_root_from_factors(factors, limits()).unwrap()
        else {
            panic!("9 is a square");
        };
        assert_eq!(foreign.rational_product(), &big(9));
        let dependency = seam_dependency();
        let rational = rational_root(9, &other);
        assert_eq!(
            conclude(&dependency, rational, AlgebraicSqrtOutcome::Verified(foreign)),
            Err(invariant("verified root belongs to a different context"))
        );
    }

    #[test]
    fn inconsistent_rational_image_is_an_invariant_violation() {
        // Same product and context, but a rational image of 8 for the product 49: 8² = 64 is not
        // congruent to 24² = 576 modulo 527, which the real roots rule out.
        let dependency = seam_dependency();
        let algebraic = AlgebraicSqrtOutcome::Verified(verified_root(&dependency));
        let forged = RationalRoot::forged(big(49), big(8), big(8));
        assert_eq!(
            conclude(&dependency, forged, algebraic),
            Err(invariant("verified roots do not give a congruence of squares"))
        );
    }

    #[test]
    fn non_roots_pass_through_conclude_unchanged() {
        use crate::sqrt::context::{IncidentalSource, ResourceKind};
        let dependency = seam_dependency();
        let rational = rational_root(49, dependency.context());
        let spent = vec![ResourceExhaustion { kind: ResourceKind::Precision, consumed: 6, limit: 1 }];
        assert_eq!(
            conclude(
                &dependency,
                rational.clone(),
                AlgebraicSqrtOutcome::ResourceExhausted(spent.clone())
            ),
            Ok(AttemptOutcome::ResourceExhausted(spent))
        );
        assert_eq!(
            conclude(
                &dependency,
                rational.clone(),
                AlgebraicSqrtOutcome::Unsupported(UnsupportedDomain::MapDenominatorNotUnit)
            ),
            Ok(AttemptOutcome::Unsupported(UnsupportedDomain::MapDenominatorNotUnit))
        );
        let incidental =
            IncidentalFactor::new(&big(527), big(17), IncidentalSource::MapDenominator).unwrap();
        assert_eq!(
            conclude(
                &dependency,
                rational,
                AlgebraicSqrtOutcome::Incidental(incidental.clone())
            ),
            Ok(AttemptOutcome::Incidental(incidental))
        );
    }

    #[test]
    fn congruence_is_checked_before_any_gcd() {
        // 35 does not divide 2² − 3² = −5, so the inputs are not a congruence even though
        // gcd(2 − 3, 35) = 1 would otherwise read as trivial.
        assert_eq!(
            factor_from_congruence(&big(2), &big(3), &big(35)),
            Err(CongruenceError::NotCongruent)
        );
        for n in [-5, 0, 1] {
            assert_eq!(
                factor_from_congruence(&big(6), &big(1), &big(n)),
                Err(CongruenceError::ModulusNotGreaterThanOne { n: big(n) })
            );
        }
    }

    #[test]
    fn divisor_is_the_difference_gcd_and_trivial_cases_report_no_divisor() {
        // 6² = 36 ≡ 1 ≡ 1² (mod 35). gcd(5, 35) = 5 is the divisor; the sum gives gcd(7, 35) = 7,
        // which is the cofactor.
        assert_eq!(
            factor_from_congruence(&big(6), &big(1), &big(35)),
            Ok(CongruenceOutcome::Proper { divisor: big(5), cofactor: big(7) })
        );
        // x = y = 5: the sum 10 shares the factor 5 with 35, but x ≡ y, so there is no divisor.
        // Likewise x ≡ −y: 6 + 29 = 35.
        for (x, y) in [(5, 5), (6, 29), (-6, 6), (0, 0)] {
            assert_eq!(
                factor_from_congruence(&big(x), &big(y), &big(35)),
                Ok(CongruenceOutcome::Trivial),
                "({x}, {y})"
            );
        }
    }

    #[test]
    fn trivial_congruence_with_a_shared_factor_reports_no_divisor() {
        // x = y = 5 modulo 35: gcd(x + y, 35) = gcd(10, 35) = 5 is a proper divisor of 35, but it
        // arises only because x ≡ y and x shares a factor with 35, so it is deliberately not
        // reported and the congruence is classified as trivial.
        let (x, y, n) = (big(5), big(5), big(35));
        assert_eq!(gcd(&(&x + &y), &n), big(5));
        assert_eq!(factor_from_congruence(&x, &y, &n), Ok(CongruenceOutcome::Trivial));
    }

    #[test]
    fn display_covers_every_error() {
        let violation = InvariantViolation { what: "w" };
        let invalid = FactorError::InvalidCandidate {
            candidate: 2,
            error: DependencyError::Invariant(violation),
        };
        assert!(invalid.to_string().starts_with("candidate 2 is not a valid dependency"));
        let sqrt = FactorError::Sqrt(SqrtError::Invariant(violation));
        assert!(sqrt.to_string().contains("w"));
        assert!(FactorError::Invariant(violation).to_string().contains("w"));
        for error in [&invalid, &sqrt, &FactorError::Invariant(violation)] {
            assert!(std::error::Error::source(error).is_some());
        }
        assert!(CongruenceError::NotCongruent.to_string().contains("congruent"));
        assert!(
            CongruenceError::ModulusNotGreaterThanOne { n: big(1) }.to_string().contains("got 1")
        );
    }
}
