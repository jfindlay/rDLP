//! Square root stage for GNFS: checked dependencies, verified rational and algebraic square
//! roots, and congruence assembly into a divisor of `N`.
//!
//! This module is the entry point for the `gnfs::sqrt` sub-crate. It provides:
//!
//! - [`rational`] — the rational square root of a [`CheckedDependency`]: the exact integer
//!   product `P = ∏(a − b·m)` is required to be positive and a perfect square, and its root
//!   `x₀` is reduced modulo `N`. Entry point: [`rational::rational_sqrt`].
//! - [`algebraic`] — the verified algebraic square root of a [`CheckedDependency`]: exactly
//!   [`reconstruct_root_from_factors`] on its relation factors. Entry point:
//!   [`algebraic::algebraic_sqrt`].
//! - [`assembly`] — congruence assembly and the checked driver: `x² ≡ y² (mod N)` is checked
//!   before any gcd, one dependency is run to a typed [`AttemptOutcome`], and a bounded candidate
//!   list is retried into a [`FactorReport`]. Entry points: [`assembly::factor_from_congruence`],
//!   [`assembly::attempt_factor`], [`assembly::factor`].
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
//! from the linear algebra step selects filtered-matrix rows, and [`DependencyContext`] lifts the
//! selection to a set S of original relations by symmetric difference of provenance. The parity
//! and character columns make the product of S a *candidate* for a square: they are necessary
//! filters, so exact verification decides. For the lifted relations the stage computes:
//!
//! 1. **Rational square root x**: `x₀² = P = ∏_{i ∈ S}(a_i − b_i·m)` exactly over `ℤ`, with
//!    `P > 0` and `P` a perfect square checked here, and `x = x₀ mod N`.
//! 2. **Algebraic square root y**: a root `β` with `β² = γ = ∏_{i ∈ S}(a_i − b_i·α)` in
//!    `ℚ[α]/(f)`, reconstructed by direct Newton lifting at an inert prime or by split-prime
//!    lifting with exact sign recombination, and *verified* exactly. The image is
//!    `y = β(c·m) mod N`, the normalized root evaluated at `θ = c·m` with denominators checked
//!    modulo `N`. It satisfies `y² ≡ P (mod N)`; it is not `Norm(β)`.
//! 3. **Assembly**: `x² ≡ y² (mod N)` is checked, and `gcd(x − y, N)` is a proper factor of `N`
//!    when `x ≢ ±y (mod N)`.
//!
//! Every entry point that produces `y` (including [`algebraic::algebraic_sqrt`],
//! [`attempt_factor`], and [`factor`]) goes through verified reconstruction; no unverified root
//! and no real-embedding sign heuristic is used. A product that is not a square is reported by a
//! typed outcome, and only a nonsquare certificate from the reconstruction says that `γ` is not a
//! square.
//!
//! # Linear algebra → square root seam
//!
//! Candidates are untrusted at this seam. [`DependencyContext::validate_kernel_vector`] re-derives
//! kernel membership and lifted parity from the original relations and never trusts a
//! candidate's own claims. Linear combinations of computed kernel vectors are tried in a
//! bounded, deterministic order by [`factor`].

pub mod algebraic;
pub mod assembly;
pub mod context;
pub mod crt;
pub mod dependency;
pub mod padic;
pub mod rational;

pub use algebraic::algebraic_sqrt;
pub use assembly::{
    AttemptOutcome, AttemptRecord, CongruenceError, CongruenceOutcome, FactorError, FactorOutcome,
    FactorReport, VerifiedFactor, attempt_factor, factor, factor_from_congruence,
};
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
pub use rational::{RationalRoot, RationalSqrtOutcome, rational_sqrt};

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
