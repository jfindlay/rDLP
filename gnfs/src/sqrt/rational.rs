//! Rational square root of a checked dependency's rational product.
//!
//! A checked dependency fixes the exact integer product `P = ∏(a − b·m)` of its relations. The
//! rational side of the congruence of squares is an integer `x₀` with `x₀² = P`, taken modulo `N`.
//! This module checks both preconditions itself, `P > 0` and `P` an exact integer square, rather
//! than trusting the sign and factor-base columns that selected the dependency: those columns make
//! the condition expected, not certain, and a failure is reported as a typed outcome, never a
//! panic.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Signed;
use shared_bigint::isqrt;

use crate::sqrt::context::SqrtContext;
use crate::sqrt::dependency::CheckedDependency;

/// An exact integer square root `x₀` of the rational product `P`, with its image modulo `N`.
///
/// Invariants: `P > 0`, `x₀ > 0`, `x₀² = P` exactly, and the image is `x₀ mod N` in `[0, N)`. The
/// fields are private; the only constructor is [`rational_sqrt`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalRoot {
    product: BigInt,
    root: BigInt,
    image: BigInt,
}

impl RationalRoot {
    /// The rational product `P = ∏(a − b·m)` over `Z`.
    #[must_use]
    pub fn product(&self) -> &BigInt {
        &self.product
    }

    /// The exact positive root `x₀` with `x₀² = P`.
    #[must_use]
    pub fn root(&self) -> &BigInt {
        &self.root
    }

    /// The image `x = x₀ mod N`, in `[0, N)`.
    #[must_use]
    pub fn modular_image(&self) -> &BigInt {
        &self.image
    }
}

/// The result of the rational square root.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RationalSqrtOutcome {
    /// `P` is a positive exact square.
    Root(RationalRoot),
    /// `P` is not a positive exact square: it is negative, zero, or not a perfect square.
    NonSquare {
        /// The rational product that failed.
        product: BigInt,
    },
}

/// The exact square root of the rational product of a checked dependency.
///
/// The product `P` is `dependency.rational_product()`. A root is returned only when `P > 0` and
/// `P = x₀²` for an integer `x₀`; the image is reduced modulo the dependency's `N`. For a
/// validated dependency the sign and rational-prime columns make `NonSquare` unreachable, so it
/// is a typed report of an inconsistency in the input and never a certificate about the algebraic
/// side.
pub fn rational_sqrt(dependency: &CheckedDependency) -> RationalSqrtOutcome {
    root_of_product(dependency.rational_product(), dependency.context())
}

/// The root of `product` reduced modulo the context's `N`, or the offending product.
pub(crate) fn root_of_product(product: &BigInt, context: &SqrtContext) -> RationalSqrtOutcome {
    let root = if product.is_positive() { isqrt(product) } else { None };
    match root {
        Some(root) => RationalSqrtOutcome::Root(RationalRoot {
            product: product.clone(),
            image: root.mod_floor(context.modulus()),
            root,
        }),
        None => RationalSqrtOutcome::NonSquare { product: product.clone() },
    }
}

#[cfg(test)]
impl RationalRoot {
    /// A root with arbitrary parts, only to test that consumers re-check what they are given.
    pub(crate) fn forged(product: BigInt, root: BigInt, image: BigInt) -> Self {
        Self { product, root, image }
    }
}

#[cfg(test)]
mod tests {
    use shared_numfield::IntPoly;

    use super::*;

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    /// x² − 2 at m = 23 modulo 527: f(23) = 529 − 2 = 527.
    fn context() -> SqrtContext {
        let f = IntPoly::from_coeffs(vec![big(-2), big(0), big(1)]);
        SqrtContext::new(&f, &big(23), &big(527)).unwrap()
    }

    #[test]
    fn positive_square_gives_root_and_reduced_image() {
        // 49 = 7² and 7 < 527. 360000 = 600² and 600 = 527 + 73, so the image of x₀ = 600 is 73.
        let ctx = context();
        let RationalSqrtOutcome::Root(root) = root_of_product(&big(49), &ctx) else {
            panic!("49 is a square");
        };
        assert_eq!((root.product(), root.root(), root.modular_image()), (&big(49), &big(7), &big(7)));
        let RationalSqrtOutcome::Root(root) = root_of_product(&big(360_000), &ctx) else {
            panic!("360000 is a square");
        };
        assert_eq!(root.root(), &big(600));
        assert_eq!(root.modular_image(), &big(73));
    }

    #[test]
    fn nonpositive_and_nonsquare_products_are_reported() {
        let ctx = context();
        for product in [big(50), big(2), big(0), big(-49), big(-1)] {
            assert_eq!(
                root_of_product(&product, &ctx),
                RationalSqrtOutcome::NonSquare { product: product.clone() },
                "{product}"
            );
        }
    }
}
