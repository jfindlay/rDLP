//! Square root stage for GNFS: rational and algebraic square roots, and final assembly.
//!
//! This module is the entry point for the `gnfs::sqrt` sub-crate. It provides:
//!
//! - [`rational`] — rational square root: given a `KernelVector`, recover the relation index
//!   set S via `expand_provenance`, form the product ∏(a − bm) over ℤ, extract its integer
//!   square root X via `isqrt`, and reduce X mod N. Entry point: [`rational::rational_sqrt`].
//! - [`algebraic`] — algebraic square root via Couveignes' CRT algorithm in ℤ[α].
//!   Entry point: [`algebraic::algebraic_sqrt`].
//! - [`assembly`] — final assembly: combine X and Y, compute gcd(X − Y, N) to extract a
//!   non-trivial factor. Entry points: [`assembly::factor_from_congruence`], [`assembly::factor`].
//! - [`context`] — checked polynomial/map context ([`SqrtContext`]), the relation-factor input
//!   boundary ([`RelationFactors`]), non-forgeable root witnesses ([`VerifiedRoot`]), failure
//!   classes, and resource limits for the algebraic square root.
//! - [`dependency`] — checked dependencies: a candidate set of filtered-matrix rows validated
//!   against the original relations ([`DependencyContext`]), lifted by symmetric difference of
//!   provenance, with the exact rational product `P` and algebraic product `γ`
//!   ([`CheckedDependency`]).
//! - [`padic`] — direct reconstruction of the algebraic square root by inverse-square-root Newton
//!   lifting at one inert prime. Entry point: [`padic::direct_root`].
//! - [`crt`] — signed local contributions of a square root at completely split primes, lifted to
//!   a common prime power with Lagrange interpolation and Chinese-remainder weights, and the exact
//!   meet-in-the-middle sign search that turns them into a root. Entry points:
//!   [`crt::split_contributions`], [`crt::recombine_signs`], [`crt::split_prime_root`].
//! - [`reconstruct_root`] — the method-neutral entry point: direct Newton lifting first, split-prime
//!   lifting only when the search for an inert prime is exhausted. Both methods return through the
//!   same exact verifier.
//!
//! # Background
//!
//! The square root stage is the terminal stage of the GNFS factoring pipeline. A `KernelVector`
//! from the linear algebra step encodes a subset S of relations whose combined rational
//! and algebraic norms are each perfect squares. This module recovers the original (a, b) pairs
//! through the filtering provenance map and computes:
//!
//! 1. **Rational square root X**: X² ≡ ∏_{i ∈ S}(a_i − b_i·m) (mod N).
//! 2. **Algebraic square root Y**: Y = |Norm(β)| mod N where β² = ∏_{i ∈ S}(a_i − b_i·α).
//! 3. **Assembly**: gcd(X − Y, N) is a non-trivial factor of N when X ≢ ±Y (mod N).
//!
//! # Linear algebra → square root seam
//!
//! `KernelVector::expand_provenance` is the linear algebra → square root seam, over-specified
//! at the linear algebra substrate precisely for this consumer. The rational square root step
//! is its first real client.

pub mod algebraic;
pub mod assembly;
pub mod context;
pub mod crt;
pub mod dependency;
pub mod padic;
pub mod rational;

pub use algebraic::algebraic_sqrt;
pub use assembly::{factor, factor_from_congruence};
pub use context::normalize::{
    CandidateRejection, CandidateVerdict, NormalizedTarget, TargetOutcome,
};
pub use context::{
    AlgebraicSqrtOutcome, CrtLimits, DirectLimits, FailureClass, IncidentalFactor,
    IncidentalSource, InputRejection, InvalidInput, InvariantViolation, NonSquareCertificate,
    NonSquareEvidence, PrecisionLimits, PrimeSearchLimits, RecombinationLimits,
    ReconstructionMethod, RelationFactors, ResourceExhaustion, ResourceKind, SqrtContext,
    SqrtError, SqrtLimits, UnsupportedDomain, VerifiedRoot,
};
pub use crt::{
    Contribution, ContributionSet, ContributionSign, LocalPrime, SplitOutcome, recombine_signs,
    split_contributions, split_prime_root,
};
pub use dependency::{
    CheckedDependency, DependencyContext, DependencyContextError, DependencyError,
    DependencyRejection,
};
pub use padic::direct_root;
pub use rational::rational_sqrt;

/// Reconstruct a verified algebraic square root with bounded primary-to-fallback dispatch.
///
/// The primary method is [`direct_root`]. The fallback, [`split_prime_root`], runs only when the
/// primary reports [`ResourceKind::InertPrimeSearch`] or [`ResourceKind::MachinePrimeRange`]
/// exhaustion, which says the inert-prime search stopped, not that no inert prime exists. Every
/// other primary outcome is returned unchanged: a verified root, a nonsquare certificate, an
/// incidental factor, an unsupported map, or an exhausted nonresidue or precision budget. If the
/// fallback is also exhausted, the two resource lists are concatenated in the order the methods
/// ran. Both methods accept a root only through [`NormalizedTarget::verify_candidate`], so a
/// [`AlgebraicSqrtOutcome::Verified`] result is the same kind of witness whichever produced it.
///
/// `limits.precision` is checked against the target's minimum modulus `2B + 1` before any prime
/// search, so a target built under a looser budget is still refused here.
///
/// # Errors
///
/// [`SqrtError::Invariant`] if either method detects a violated internal identity. These indicate
/// a defect, never a nonsquare.
pub fn reconstruct_root(
    target: &NormalizedTarget,
    limits: SqrtLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    let bits = target.minimum_modulus().bits();
    if bits > limits.precision.max_modulus_bits {
        return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![
            ResourceExhaustion {
                kind: ResourceKind::Precision,
                consumed: bits,
                limit: limits.precision.max_modulus_bits,
            },
        ]));
    }
    let primary = match direct_root(target, limits.direct)? {
        AlgebraicSqrtOutcome::ResourceExhausted(spent) if inert_search_ended(&spent) => spent,
        other => return Ok(other),
    };
    Ok(
        match split_prime_root(target, limits.crt, limits.recombination)? {
            AlgebraicSqrtOutcome::ResourceExhausted(fallback) => {
                AlgebraicSqrtOutcome::ResourceExhausted(
                    primary.into_iter().chain(fallback).collect(),
                )
            }
            other => other,
        },
    )
}

/// Build the normalized target for `factors` under `limits.precision`, then [`reconstruct_root`].
///
/// A precision budget too small for the target's bound is reported as
/// [`ResourceKind::Precision`] exhaustion before any prime search.
///
/// # Errors
///
/// [`SqrtError`] as for [`NormalizedTarget::new`] and [`reconstruct_root`].
pub fn reconstruct_root_from_factors(
    factors: RelationFactors,
    limits: SqrtLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    match NormalizedTarget::new(factors, limits.precision)? {
        TargetOutcome::Ready(target) => reconstruct_root(&target, limits),
        TargetOutcome::ResourceExhausted(exhausted) => {
            Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![exhausted]))
        }
    }
}

/// True for a non-empty exhaustion list made only of inert-prime-search outcomes.
fn inert_search_ended(spent: &[ResourceExhaustion]) -> bool {
    !spent.is_empty()
        && spent.iter().all(|e| {
            matches!(
                e.kind,
                ResourceKind::InertPrimeSearch | ResourceKind::MachinePrimeRange
            )
        })
}
