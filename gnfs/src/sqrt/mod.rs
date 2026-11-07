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
//! - [`padic`] — direct reconstruction of the algebraic square root by inverse-square-root Newton
//!   lifting at one inert prime. Entry point: [`padic::direct_root`].
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
pub use padic::direct_root;
pub use rational::rational_sqrt;
