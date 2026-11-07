//! Algebraic square root of a checked dependency, with exact verification.
//!
//! A checked dependency carries the relation pairs `(a, b)` whose algebraic factors `a − b·α`
//! multiply to `γ`. The root is recovered by [`crate::sqrt::reconstruct_root_from_factors`]: direct
//! inverse-square-root Newton lifting at an inert prime, with split-prime lifting and exact sign
//! recombination as the fallback. Whichever method runs, a root is accepted only through
//! [`crate::sqrt::NormalizedTarget::verify_candidate`], which checks `β² = γ` over `Q` and the
//! modular image `y² ≡ P (mod N)`. The result is a [`VerifiedRoot`] or a typed non-success.

use crate::sqrt::context::{AlgebraicSqrtOutcome, SqrtError, SqrtLimits};
use crate::sqrt::dependency::CheckedDependency;
use crate::sqrt::reconstruct_root_from_factors;

/// The verified algebraic square root of the dependency's product `γ`, or a typed non-success.
///
/// This is exactly [`reconstruct_root_from_factors`] applied to the dependency's relation
/// factors. Only a [`AlgebraicSqrtOutcome::NotSquare`] certificate says that `γ` is not a square;
/// resource exhaustion and unsupported maps never do.
///
/// # Errors
///
/// [`SqrtError`] if the factors are invalid or an internal identity fails. Neither is a nonsquare.
///
/// [`VerifiedRoot`]: crate::sqrt::VerifiedRoot
pub fn algebraic_sqrt(
    dependency: &CheckedDependency,
    limits: SqrtLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    reconstruct_root_from_factors(dependency.factors().clone(), limits)
}
