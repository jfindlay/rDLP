//! Checked relations: validated records, prime classification, and the context that binds them.
//!
//! The legacy collector in [`crate::sieve`] produces [`Relation`] records whose exponent vectors
//! are trusted. This module re-derives every identity a downstream consumer relies on, and it
//! mints the validated value only on success:
//!
//! - [`RelationChecker`] owns a [`SqrtContext`] and a validated snapshot of the factor base.
//! - [`CheckedRelation`] is a relation whose exponents are exactly the `p`-adic valuations of
//!   its exact norms, whose affine ideal labels satisfy `a ≡ b·r (mod p)`, and whose sign is the
//!   sign of `a − b·m`. Only [`RelationChecker::check`] constructs one.
//! - [`CheckedRelations`] is an ordered, duplicate-free list of checked relations, so that
//!   position, input index, and matrix provenance coincide.
//! - [`classify_prime`] decides which primes carry an algebraic parity constraint.
//!
//! # Which algebraic parities are constraints
//!
//! For a prime `p ∤ c·disc(F)` the order `Z[θ]` is maximal at `p`, so each simple root `r` of
//! `f mod p` names a degree-1 prime ideal `(p, θ − c·r)`. For a coprime pair, `p | f_hom(a, b)`
//! forces `p ∤ b`; then `r ≡ a/b` is the only root and `v_p(f_hom(a, b))` is the valuation of
//! `a − bα` at that ideal. Its parity is a genuine square constraint. At `p | c·disc(F)` the
//! recorded exponent is only a norm-valuation attribution, and for non-monic `f` the norm of
//! `a − bα` is `f_hom(a, b)/c`, so its parity is not a constraint at all. Such primes are
//! [`PrimeClass::Exceptional`]; their columns are kept in the matrix layout and never set. The
//! rational side is unaffected: every rational prime and the sign remain constraints.
//!
//! Omitting exceptional parities weakens the candidate filter. It does not weaken the exact
//! square-root verification that follows, which is authoritative.
//!
//! # Inherited limitation
//!
//! The primitive affine collector records only degree-1 affine ideals `(p, α − r)`. A coprime pair
//! with `p | c`, `p | b` and `p | f_hom(a, b)` has no affine root `r` with `a ≡ b·r`; the collector
//! never returns it, and a hand-built record claiming one is rejected here. Complete projective
//! and exceptional-ideal factorization is not provided.

use std::collections::HashMap;
use std::fmt;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use shared_numfield::is_prime_u64;

use super::{algebraic_norm, rational_norm, FactorBase, Relation};
use crate::sqrt::SqrtContext;

// ─── prime classification ────────────────────────────────────────────────────

/// Whether a prime contributes honest algebraic parity constraints.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrimeClass {
    /// `p ∤ c·disc(F)`: every ideal above `p` is a degree-1 prime with a valid valuation.
    Good,
    /// `p | c·disc(F)`: the whole prime is conservatively excluded from algebraic parity.
    Exceptional,
}

/// Classify `p` against the leading coefficient `c` and discriminant `disc(F)` of `ctx`.
///
/// The result is [`PrimeClass::Exceptional`] exactly when `p | c·disc(F)`, tested with exact
/// integer remainders, so for a prime `p` it is the same as `p | c` or `p | disc(F)`. Values below
/// 2 are not primes and are classified as exceptional, so that no constraint is claimed for them.
/// Every ideal above an exceptional prime is exceptional, even one with a simple root: the rule
/// is conservative and easy to justify, and its price is that some genuine constraints are
/// dropped.
#[must_use]
pub fn classify_prime(ctx: &SqrtContext, p: u64) -> PrimeClass {
    if p < 2 {
        return PrimeClass::Exceptional;
    }
    let product = ctx.leading_coeff() * ctx.discriminant();
    if (product % BigInt::from(p)).is_zero() { PrimeClass::Exceptional } else { PrimeClass::Good }
}

// ─── rejection types ─────────────────────────────────────────────────────────

/// Which exponent list of a relation an error refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExponentSide {
    /// The rational factor-base primes.
    Rational,
    /// The algebraic factor-base ideals.
    Algebraic,
}

impl fmt::Display for ExponentSide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Rational => f.write_str("rational"),
            Self::Algebraic => f.write_str("algebraic"),
        }
    }
}

/// A factor base that [`RelationChecker::new`] refuses. `position` is a slice position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FactorBaseRejection {
    /// The rational primes are not strictly ascending at `position`.
    RationalNotStrictlyAscending {
        /// Position of the first entry not greater than its predecessor.
        position: usize,
    },
    /// The rational entry at `position` is not prime.
    RationalNotPrime {
        /// Position of the entry.
        position: usize,
        /// The offending value.
        value: u64,
    },
    /// The ideals are not strictly ascending by `(p, r)` at `position`.
    IdealsNotStrictlyAscending {
        /// Position of the first ideal not greater than its predecessor.
        position: usize,
    },
    /// The ideal at `position` carries an `index` different from `position`.
    IdealIndexMismatch {
        /// Position of the ideal.
        position: usize,
        /// The stored index.
        index: usize,
    },
    /// The ideal at `position` lies above a non-prime `p`.
    IdealNotPrime {
        /// Position of the ideal.
        position: usize,
        /// The offending value.
        p: u64,
    },
    /// The ideal's root `r` is not below `p`.
    IdealRootOutOfRange {
        /// Position of the ideal.
        position: usize,
        /// The prime.
        p: u64,
        /// The root.
        r: u64,
    },
    /// `f(r) ≢ 0 (mod p)`, with `f` taken from the context.
    IdealNotRoot {
        /// Position of the ideal.
        position: usize,
        /// The prime.
        p: u64,
        /// The claimed root.
        r: u64,
    },
}

/// A relation that [`RelationChecker::check`] refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RelationRejection {
    /// `gcd(|a|, |b|) ≠ 1`.
    NotCoprime {
        /// The actual gcd.
        gcd: BigInt,
    },
    /// An exponent list is not strictly ascending in index, or has exponent 0.
    MalformedExponents {
        /// The list.
        side: ExponentSide,
        /// Position of the first offending entry.
        position: usize,
    },
    /// An index lies outside the factor base.
    IndexOutOfRange {
        /// The list.
        side: ExponentSide,
        /// The offending index.
        index: usize,
        /// The size of the corresponding factor base.
        size: usize,
    },
    /// `a − b·m = 0`.
    ZeroRationalNorm,
    /// The stored sign differs from the sign of `a − b·m`.
    SignMismatch {
        /// `a − b·m < 0`.
        expected: bool,
        /// The stored flag.
        stored: bool,
    },
    /// `Π p^e` over the rational entries is not `|a − b·m|`.
    RationalNormMismatch {
        /// `|a − b·m|`.
        norm: BigInt,
    },
    /// `f_hom(a, b) = 0`.
    ZeroAlgebraicNorm,
    /// An ideal `(p, r)` is listed although `a ≢ b·r (mod p)`.
    AffineLabelMismatch {
        /// Position of the entry in the algebraic list.
        position: usize,
        /// The prime.
        p: u64,
        /// The root.
        r: u64,
    },
    /// `Π p^e` over the algebraic entries is not `|f_hom(a, b)|`.
    AlgebraicNormMismatch {
        /// `|f_hom(a, b)|`.
        norm: BigInt,
    },
}

/// A relation list that [`CheckedRelations::new`] refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckedRelationsError {
    /// Relation `index` fails [`RelationChecker::check`].
    Relation {
        /// Position in the input slice.
        index: usize,
        /// The reason.
        rejection: RelationRejection,
    },
    /// Relations `first < second` have the identical pair `(a, b)`.
    DuplicatePair {
        /// Position of the earlier occurrence.
        first: usize,
        /// Position of the repeated occurrence.
        second: usize,
    },
}

impl fmt::Display for FactorBaseRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RationalNotStrictlyAscending { position } => {
                write!(f, "rational primes are not strictly ascending at position {position}")
            }
            Self::RationalNotPrime { position, value } => {
                write!(f, "rational entry {value} at position {position} is not prime")
            }
            Self::IdealsNotStrictlyAscending { position } => {
                write!(f, "ideals are not strictly ascending by (p, r) at position {position}")
            }
            Self::IdealIndexMismatch { position, index } => {
                write!(f, "ideal at position {position} carries index {index}")
            }
            Self::IdealNotPrime { position, p } => {
                write!(f, "ideal at position {position} lies above non-prime {p}")
            }
            Self::IdealRootOutOfRange { position, p, r } => {
                write!(f, "ideal ({p}, {r}) at position {position} has root not below p")
            }
            Self::IdealNotRoot { position, p, r } => {
                write!(f, "ideal ({p}, {r}) at position {position}: f({r}) is not 0 mod {p}")
            }
        }
    }
}

impl fmt::Display for RelationRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotCoprime { gcd } => write!(f, "gcd(a, b) = {gcd}, not 1"),
            Self::MalformedExponents { side, position } => {
                write!(f, "{side} exponents malformed at position {position}")
            }
            Self::IndexOutOfRange { side, index, size } => {
                write!(f, "{side} index {index} is outside a factor base of size {size}")
            }
            Self::ZeroRationalNorm => f.write_str("rational norm a - b*m is zero"),
            Self::SignMismatch { expected, stored } => {
                write!(f, "rational sign is negative={expected}, but negative={stored} is stored")
            }
            Self::RationalNormMismatch { norm } => {
                write!(f, "rational exponents do not reconstruct |a - b*m| = {norm}")
            }
            Self::ZeroAlgebraicNorm => f.write_str("homogeneous norm f_hom(a, b) is zero"),
            Self::AffineLabelMismatch { position, p, r } => {
                write!(f, "algebraic entry {position}: ideal ({p}, {r}) needs a = b*{r} mod {p}")
            }
            Self::AlgebraicNormMismatch { norm } => {
                write!(f, "algebraic exponents do not reconstruct |f_hom(a, b)| = {norm}")
            }
        }
    }
}

impl fmt::Display for CheckedRelationsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Relation { index, rejection } => write!(f, "relation {index}: {rejection}"),
            Self::DuplicatePair { first, second } => {
                write!(f, "relations {first} and {second} have the same pair (a, b)")
            }
        }
    }
}

impl std::error::Error for FactorBaseRejection {}
impl std::error::Error for RelationRejection {}
impl std::error::Error for CheckedRelationsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Relation { rejection, .. } => Some(rejection),
            Self::DuplicatePair { .. } => None,
        }
    }
}

// ─── the checker ─────────────────────────────────────────────────────────────

/// A [`SqrtContext`] together with a validated snapshot of a factor base.
///
/// The snapshot fixes the column order: [`rational_primes`](Self::rational_primes) and
/// [`ideals`](Self::ideals) are listed in factor-base index order. Equality of checkers is the
/// binding used by relation sets and matrix layouts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RelationChecker {
    context: SqrtContext,
    rational_primes: Vec<u64>,
    ideals: Vec<(u64, u64)>,
    classes: Vec<PrimeClass>,
}

impl RelationChecker {
    /// Validate `fb` against `ctx` and snapshot it.
    ///
    /// Completeness of the factor base is not required: a missing ideal is caught by exact norm
    /// reconstruction. `FactorBase::obstruction_count` and `AlgebraicPrime::is_bad_prime` are not
    /// read; prime classification uses [`classify_prime`].
    ///
    /// # Errors
    ///
    /// Each check runs over the whole list before the next begins, and the first failure is
    /// reported, in the order of the [`FactorBaseRejection`] variants: rational primes strictly
    /// ascending, then prime; ideals strictly ascending by `(p, r)`, then `index == position`,
    /// `p` prime, `r < p`, and `f(r) ≡ 0 (mod p)`.
    pub fn new(ctx: &SqrtContext, fb: &FactorBase) -> Result<Self, FactorBaseRejection> {
        let rational = &fb.rational_primes;
        let ideals = &fb.algebraic_ideals;
        if let Some(position) = (1..rational.len()).find(|&i| rational[i] <= rational[i - 1]) {
            return Err(FactorBaseRejection::RationalNotStrictlyAscending { position });
        }
        if let Some((position, &value)) =
            rational.iter().enumerate().find(|&(_, &p)| !is_prime_u64(p))
        {
            return Err(FactorBaseRejection::RationalNotPrime { position, value });
        }
        let key = |i: usize| (ideals[i].p, ideals[i].r);
        if let Some(position) = (1..ideals.len()).find(|&i| key(i) <= key(i - 1)) {
            return Err(FactorBaseRejection::IdealsNotStrictlyAscending { position });
        }
        if let Some((position, ideal)) = ideals.iter().enumerate().find(|(i, a)| a.index != *i) {
            return Err(FactorBaseRejection::IdealIndexMismatch { position, index: ideal.index });
        }
        if let Some((position, ideal)) =
            ideals.iter().enumerate().find(|(_, a)| !is_prime_u64(a.p))
        {
            return Err(FactorBaseRejection::IdealNotPrime { position, p: ideal.p });
        }
        if let Some((position, ideal)) = ideals.iter().enumerate().find(|(_, a)| a.r >= a.p) {
            return Err(FactorBaseRejection::IdealRootOutOfRange {
                position,
                p: ideal.p,
                r: ideal.r,
            });
        }
        let f = ctx.algebraic_poly();
        let not_root = ideals
            .iter()
            .enumerate()
            .find(|(_, a)| !f.eval_mod(&BigInt::from(a.r), &BigInt::from(a.p)).is_zero());
        if let Some((position, ideal)) = not_root {
            return Err(FactorBaseRejection::IdealNotRoot {
                position,
                p: ideal.p,
                r: ideal.r,
            });
        }
        Ok(Self {
            context: ctx.clone(),
            rational_primes: rational.clone(),
            ideals: ideals.iter().map(|a| (a.p, a.r)).collect(),
            classes: ideals.iter().map(|a| classify_prime(ctx, a.p)).collect(),
        })
    }

    /// The owned context: `N`, the exact `m`, `f`, `c`, and `disc(F)`.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        &self.context
    }

    /// The rational primes in column order.
    #[must_use]
    pub fn rational_primes(&self) -> &[u64] {
        &self.rational_primes
    }

    /// The algebraic ideals `(p, r)` in column order.
    #[must_use]
    pub fn ideals(&self) -> &[(u64, u64)] {
        &self.ideals
    }

    /// The class of the prime below ideal `index`, or `None` when `index` is out of range.
    #[must_use]
    pub fn ideal_class(&self, index: usize) -> Option<PrimeClass> {
        self.classes.get(index).copied()
    }

    /// Validate `relation` and mint the checked value.
    ///
    /// Any sign of `b` is accepted. Reconstruction divides `|norm|` by each `p` exactly `e`
    /// times and requires remainder 1; it never forms `p^e` from the untrusted exponent, so the
    /// work is bounded by the bit length of the norm.
    ///
    /// # Errors
    ///
    /// The first failing check is reported, in this order:
    ///
    /// 1. [`RelationRejection::NotCoprime`].
    /// 2. For the rational list and then the algebraic list:
    ///    [`RelationRejection::MalformedExponents`] (indices strictly ascending, exponents at
    ///    least 1), then [`RelationRejection::IndexOutOfRange`].
    /// 3. [`RelationRejection::ZeroRationalNorm`], [`RelationRejection::SignMismatch`],
    ///    [`RelationRejection::RationalNormMismatch`].
    /// 4. [`RelationRejection::ZeroAlgebraicNorm`], [`RelationRejection::AffineLabelMismatch`]
    ///    for every entry whether good or exceptional, then
    ///    [`RelationRejection::AlgebraicNormMismatch`].
    pub fn check(&self, relation: &Relation) -> Result<CheckedRelation, RelationRejection> {
        let (a, b) = (&relation.a, &relation.b);
        let g = a.gcd(b);
        if !g.is_one() {
            return Err(RelationRejection::NotCoprime { gcd: g });
        }
        let rational = &relation.rational_exponents.entries;
        let algebraic = &relation.algebraic_exponents.entries;
        check_list(ExponentSide::Rational, rational, self.rational_primes.len())?;
        check_list(ExponentSide::Algebraic, algebraic, self.ideals.len())?;

        let rat_norm = rational_norm(a, b, self.context.root());
        if rat_norm.is_zero() {
            return Err(RelationRejection::ZeroRationalNorm);
        }
        let expected = rat_norm.is_negative();
        if expected != relation.rational_sign {
            return Err(RelationRejection::SignMismatch {
                expected,
                stored: relation.rational_sign,
            });
        }
        let magnitude = rat_norm.abs();
        let rational_factors = rational.iter().map(|&(i, e)| (self.rational_primes[i], e));
        if !divides_out(magnitude.clone(), rational_factors) {
            return Err(RelationRejection::RationalNormMismatch { norm: magnitude });
        }

        let alg_norm = algebraic_norm(a, b, self.context.algebraic_poly());
        if alg_norm.is_zero() {
            return Err(RelationRejection::ZeroAlgebraicNorm);
        }
        for (position, &(i, _)) in algebraic.iter().enumerate() {
            let (p, r) = self.ideals[i];
            if !((a - b * BigInt::from(r)) % BigInt::from(p)).is_zero() {
                return Err(RelationRejection::AffineLabelMismatch { position, p, r });
            }
        }
        let magnitude = alg_norm.abs();
        let algebraic_factors = algebraic.iter().map(|&(i, e)| (self.ideals[i].0, e));
        if !divides_out(magnitude.clone(), algebraic_factors) {
            return Err(RelationRejection::AlgebraicNormMismatch { norm: magnitude });
        }

        Ok(CheckedRelation {
            a: a.clone(),
            b: b.clone(),
            rational_norm: rat_norm,
            algebraic_norm: alg_norm,
            rational_exponents: rational.clone(),
            algebraic_exponents: algebraic.clone(),
        })
    }
}

/// Entries strictly ascending with exponent at least 1, then every index below `size`.
fn check_list(
    side: ExponentSide,
    entries: &[(usize, u32)],
    size: usize,
) -> Result<(), RelationRejection> {
    let malformed = entries
        .iter()
        .enumerate()
        .find(|&(i, &(index, e))| e == 0 || (i > 0 && index <= entries[i - 1].0));
    if let Some((position, _)) = malformed {
        return Err(RelationRejection::MalformedExponents { side, position });
    }
    if let Some(&(index, _)) = entries.iter().find(|&&(index, _)| index >= size) {
        return Err(RelationRejection::IndexOutOfRange { side, index, size });
    }
    Ok(())
}

/// True iff dividing `remaining` by each `(p, e)` exactly `e` times leaves exactly 1.
///
/// Each successful division at least halves a positive value, so the loop stops after at most
/// `log2(remaining)` iterations whatever `e` is.
fn divides_out(mut remaining: BigInt, factors: impl Iterator<Item = (u64, u32)>) -> bool {
    for (p, e) in factors {
        let p = BigInt::from(p);
        for _ in 0..e {
            let (quotient, rest) = remaining.div_rem(&p);
            if !rest.is_zero() {
                return false;
            }
            remaining = quotient;
        }
    }
    remaining.is_one()
}

// ─── checked relation ────────────────────────────────────────────────────────

/// A relation validated by [`RelationChecker::check`].
///
/// Guarantees: `gcd(a, b) = 1`; `a − b·m` and `f_hom(a, b)` are nonzero; each exponent list is
/// strictly ascending with exponents at least 1 and in range; the rational exponents are the exact
/// valuations of `a − b·m`; every algebraic entry `(p, r)` satisfies `a ≡ b·r (mod p)`; and the
/// algebraic exponents are the exact valuations of `f_hom(a, b)`. No public constructor exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckedRelation {
    a: BigInt,
    b: BigInt,
    rational_norm: BigInt,
    algebraic_norm: BigInt,
    rational_exponents: Vec<(usize, u32)>,
    algebraic_exponents: Vec<(usize, u32)>,
}

impl CheckedRelation {
    /// The coordinate `a`.
    #[must_use]
    pub fn a(&self) -> &BigInt {
        &self.a
    }

    /// The coordinate `b`, of any sign.
    #[must_use]
    pub fn b(&self) -> &BigInt {
        &self.b
    }

    /// The signed, nonzero rational norm `a − b·m`.
    #[must_use]
    pub fn rational_norm(&self) -> &BigInt {
        &self.rational_norm
    }

    /// The signed, nonzero homogeneous norm `f_hom(a, b) = Σ f_i a^i b^(d−i)`.
    #[must_use]
    pub fn algebraic_norm(&self) -> &BigInt {
        &self.algebraic_norm
    }

    /// True iff the rational norm is negative.
    #[must_use]
    pub fn rational_sign(&self) -> bool {
        self.rational_norm.is_negative()
    }

    /// The `(rational prime index, exponent)` entries, strictly ascending, exponents at least 1.
    #[must_use]
    pub fn rational_exponents(&self) -> &[(usize, u32)] {
        &self.rational_exponents
    }

    /// The `(ideal index, exponent)` entries, strictly ascending, exponents at least 1.
    #[must_use]
    pub fn algebraic_exponents(&self) -> &[(usize, u32)] {
        &self.algebraic_exponents
    }
}

// ─── checked relation set ────────────────────────────────────────────────────

/// An ordered, duplicate-free list of checked relations bound to one [`RelationChecker`].
///
/// Position `i` is the input index and the provenance index of any matrix row built from it. Every
/// nonempty subset of the pairs is accepted by `RelationFactors::new` for the bound context, since
/// each pair has `a − b·m ≠ 0` and `f_hom ≠ 0`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckedRelations {
    checker: RelationChecker,
    relations: Vec<CheckedRelation>,
}

impl CheckedRelations {
    /// Check `relations` in order against a clone of `checker`; the list may be empty.
    ///
    /// Pairs `(a, b)` and `(−a, −b)` are distinct records and both accepted. Duplicate pairs are
    /// rejected, never removed or renumbered: the caller must deduplicate beforehand.
    ///
    /// # Errors
    ///
    /// [`CheckedRelationsError::Relation`] for the first relation that fails
    /// [`RelationChecker::check`]. After relation `i` passes, an exact match of `(a, b)` with an
    /// earlier pair gives [`CheckedRelationsError::DuplicatePair`] with `second = i`.
    pub fn new(
        checker: &RelationChecker,
        relations: &[Relation],
    ) -> Result<Self, CheckedRelationsError> {
        let mut seen: HashMap<(BigInt, BigInt), usize> = HashMap::with_capacity(relations.len());
        let mut checked = Vec::with_capacity(relations.len());
        for (index, relation) in relations.iter().enumerate() {
            let item = checker
                .check(relation)
                .map_err(|rejection| CheckedRelationsError::Relation { index, rejection })?;
            if let Some(&first) = seen.get(&(item.a.clone(), item.b.clone())) {
                return Err(CheckedRelationsError::DuplicatePair { first, second: index });
            }
            seen.insert((item.a.clone(), item.b.clone()), index);
            checked.push(item);
        }
        Ok(Self { checker: checker.clone(), relations: checked })
    }

    /// The bound checker.
    #[must_use]
    pub fn checker(&self) -> &RelationChecker {
        &self.checker
    }

    /// The bound context.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        self.checker.context()
    }

    /// The checked relations; position equals input index equals provenance index.
    #[must_use]
    pub fn relations(&self) -> &[CheckedRelation] {
        &self.relations
    }

    /// The number of relations.
    #[must_use]
    pub fn len(&self) -> usize {
        self.relations.len()
    }

    /// True iff there are no relations.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.relations.is_empty()
    }
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn pairs(list: &[(u64, u32)]) -> impl Iterator<Item = (u64, u32)> + '_ {
        list.iter().copied()
    }

    #[test]
    fn divides_out_requires_exact_valuations() {
        // 360 = 2³·3²·5.
        let norm = BigInt::from(360);
        assert!(divides_out(norm.clone(), pairs(&[(2, 3), (3, 2), (5, 1)])));
        // Under-claiming a valuation leaves a cofactor.
        assert!(!divides_out(norm.clone(), pairs(&[(2, 2), (3, 2), (5, 1)])));
        // Over-claiming one fails to divide.
        assert!(!divides_out(norm, pairs(&[(2, 4), (3, 2), (5, 1)])));
    }

    #[test]
    fn divides_out_terminates_for_a_huge_exponent() {
        // 2^u32::MAX is never formed: the first failed division ends the loop.
        assert!(!divides_out(BigInt::from(8), pairs(&[(2, u32::MAX)])));
        assert!(divides_out(BigInt::one(), std::iter::empty()));
    }

    #[test]
    fn check_list_orders_malformed_before_range() {
        // Out-of-order and out-of-range together: the shape error wins.
        let bad = check_list(ExponentSide::Rational, &[(5, 1), (2, 1)], 3);
        assert_eq!(
            bad,
            Err(RelationRejection::MalformedExponents {
                side: ExponentSide::Rational,
                position: 1
            })
        );
        let range = check_list(ExponentSide::Algebraic, &[(0, 1), (3, 2)], 3);
        assert_eq!(
            range,
            Err(RelationRejection::IndexOutOfRange {
                side: ExponentSide::Algebraic,
                index: 3,
                size: 3
            })
        );
        assert_eq!(check_list(ExponentSide::Rational, &[], 0), Ok(()));
    }

    #[test]
    fn error_messages_name_the_failing_item() {
        let err = CheckedRelationsError::DuplicatePair { first: 0, second: 2 };
        assert!(err.to_string().contains("0 and 2"));
        let err = CheckedRelationsError::Relation {
            index: 3,
            rejection: RelationRejection::ZeroRationalNorm,
        };
        assert!(err.to_string().starts_with("relation 3"));
        assert!(std::error::Error::source(&err).is_some());
        assert!(std::error::Error::source(&CheckedRelationsError::DuplicatePair {
            first: 0,
            second: 1
        })
        .is_none());
    }
}
