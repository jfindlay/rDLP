//! Number-field substrate: ℚ(α), element arithmetic, norm via resultant, and ideal arithmetic.
//!
//! This crate provides:
//! - [`IntPoly`] — polynomials over ℤ with `BigInt` coefficients.
//! - [`RatPoly`] — polynomials over ℚ with `BigRational` coefficients.
//! - [`NumberField`] — a number field K = ℚ(α) defined by a monic irreducible f ∈ ℤ[x].
//! - [`NumberFieldElement`] — an element of K, represented as a reduced polynomial in α.
//! - [`Ideal`] — a fractional ideal of ℤ[α] in two-element primary representation `(p, α − r)`.
//! - [`dedekind_factor`] — Dedekind factorisation of the ideal (p) in ℤ[α].
//! - [`dedekind_factor_extended`] — Dedekind factorisation reporting bad primes and whether p
//!   divides the index.
//! - [`DedekindResult`] — result type for `dedekind_factor_extended`.
//! - [`dedekind_index_test`] — the Dedekind criterion: does p divide [ℤ_K : ℤ[α]]? Returns an
//!   [`IndexTest`] or a [`DedekindError`].
//! - [`discriminant`] — discriminant of a monic polynomial f ∈ ℤ[x].
//! - [`is_bad_prime`] — test whether p | disc(f).
//! - [`FpPoly`] — polynomials over a proved-prime field F_p, with gcd, `pow_mod`, an exact
//!   irreducibility test, and budgeted root enumeration (see [`fp_poly`]).
//! - [`ResidueField`] — the field `F_p[x]/(f)` for irreducible `f`, with inverses, a squareness
//!   test, and Tonelli–Shanks square roots under a bounded nonresidue search (see
//!   [`residue_field`]).
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
pub mod residue_field;
pub mod resultant;

pub use dedekind::{
    dedekind_factor, dedekind_factor_extended, dedekind_index_test, discriminant, is_bad_prime,
    DedekindError, DedekindResult, IndexTest,
};
pub use element::{NumberField, NumberFieldElement};
pub use fp_poly::{
    ExtendedGcd, FpPoly, FpPolyError, PrimeModulus, RootEnumeration, RootSearchLimits,
    is_prime_u64,
};
pub use ideal::Ideal;
pub use poly::{IntPoly, RatPoly};
pub use residue_field::{
    NonresidueSearch, NonresidueSearchLimits, ResidueField, ResidueFieldError, SqrtOutcome,
};
pub use resultant::{resultant, subresultant_gcd};
