//! Number-field substrate: ℚ(α), element arithmetic, norm via resultant, and ideal arithmetic.
//!
//! This crate provides:
//! - [`IntPoly`] — polynomials over ℤ with `BigInt` coefficients.
//! - [`RatPoly`] — polynomials over ℚ with `BigRational` coefficients.
//! - [`NumberField`] — a number field K = ℚ(α) defined by a monic irreducible f ∈ ℤ[x].
//! - [`NumberFieldElement`] — an element of K, represented as a reduced polynomial in α.
//! - [`Ideal`] — a fractional ideal of ℤ[α] in two-element primary representation `(p, α − r)`.
//! - [`dedekind_factor`] — Dedekind factorisation of the ideal (p) in ℤ[α].
//! - [`dedekind_factor_extended`] — Dedekind factorisation with bad-prime handling.
//! - [`DedekindResult`] — result type for `dedekind_factor_extended`.
//! - [`discriminant`] — discriminant of a monic polynomial f ∈ ℤ[x].
//! - [`is_bad_prime`] — test whether p | disc(f).
//! - [`FpPoly`] — polynomials over a proved-prime field F_p, with gcd, `pow_mod`, an exact
//!   irreducibility test, and budgeted root enumeration (see [`fp_poly`]).
//!
//! # Design
//!
//! Coefficient arithmetic delegates to `num_bigint::BigInt` and `num_rational::BigRational`
//! (correctness-oracle dependencies, analogous to CADO-NFS for the sieve). The number-field
//! abstraction and all algorithms above it are first-party.
//!
//! # Invariants
//!
//! - Polynomial trailing zeros are always trimmed.
//! - `NumberFieldElement::poly` always has degree strictly less than `field.degree()`.
//!   Multiplication eagerly reduces mod f to maintain this invariant.
//! - `Ideal::p` is always positive.

pub mod dedekind;
pub mod element;
pub mod fp_poly;
pub mod ideal;
pub mod poly;
pub mod resultant;

pub use dedekind::{
    dedekind_factor, dedekind_factor_extended, discriminant, is_bad_prime, DedekindResult,
};
pub use element::{NumberField, NumberFieldElement};
pub use fp_poly::{
    ExtendedGcd, FpPoly, FpPolyError, PrimeModulus, RootEnumeration, RootSearchLimits,
    is_prime_u64,
};
pub use ideal::Ideal;
pub use poly::{IntPoly, RatPoly};
pub use resultant::{resultant, subresultant_gcd};
