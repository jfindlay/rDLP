//! Quadratic-character columns for the GF(2) matrix.
//!
//! # Individual root characters
//!
//! A [`RootCharacter`] is a pair `(q, s)` of an odd prime `q` and a root `s` of `f` modulo `q`.
//! It sends a relation `(a, b)` to the Legendre symbol of `a − b·s` modulo `q`, written as a GF(2)
//! bit: 0 for a residue, 1 for a non-residue. One prime with several roots yields several
//! characters, and the character count of a matrix counts pairs, not primes.
//!
//! The pair is valid when `q` is odd, `q ∤ c·disc(F)` (checked by
//! [`classify_prime`](crate::sieve::classify_prime)), and `f(s) ≡ 0 (mod q)`. Then `f mod q` is
//! squarefree of full degree, `s` is a simple root, and `α ↦ s` is a ring homomorphism from the
//! localisation of `Z[α]` at `q` to `F_q`. Because `q ∤ disc(F)` the order is maximal at `q`, so
//! if `Π(a − bα)` is a square `β²` in the field, `β` is integral at `q` and `Π(a − b·s) = φ(β)²`
//! is a square modulo `q`. The bits of a dependency's original relations must therefore sum to
//! zero at every valid pair.
//!
//! A pair is usable for a relation set only if no `a − b·s` vanishes modulo `q`. A zero has no
//! Legendre bit, so a vanishing pair is rejected and reported, never read as a 0 bit.
//!
//! # What the characters prove
//!
//! Characters are necessary filters. They are not square certificates and do not choose the sign
//! of a square root: a candidate that passes every sampled character can still fail to be a
//! square, for example a norm-one unit, whose norm character is trivial at every prime.
//! Exact reconstruction and verification of the root remain authoritative.
//!
//! # Provenance
//!
//! A merged row is the XOR of the original rows its provenance names, in every column including
//! character columns. The bits are therefore always computed from the original checked relations
//! and XORed over the complete provenance. They are never recomputed from merged norms.
//! [`RootCharacters::fill_columns`] fills the character columns of a freshly built or an already
//! reduced matrix by that rule, and [`RootCharacters::verify_columns`] checks the identity.
//!
//! # Column layout
//!
//! The columns are those of [`MatrixLayout`]: the sign column at `obstruction_col_start`, and
//! character slot `k` at [`MatrixLayout::character_column`]`(k)`.
//!
//! # Legacy norm characters
//!
//! [`populate_qc_columns`] and [`select_qc_primes`] keep their original behaviour for the
//! algebraic callers that still use them. They form Legendre symbols of homogeneous *norms* at
//! completely split primes. Once the norm exponents have even sum these retain at most sign
//! information, and they do not guarantee that a square root exists.
//!
//! # Legendre symbol
//!
//! For a prime `q` and integer `a`, the Legendre symbol `(a/q)` is:
//! - 0 if `q | a`
//! - 1 if `a` is a quadratic residue mod `q`
//! - −1 if `a` is a quadratic non-residue mod `q`
//!
//! Computed via Euler's criterion: `(a/q) ≡ a^((q-1)/2) (mod q)`.
//! The parity is 1 iff the symbol is −1 (i.e., `a^((q-1)/2) ≡ q-1 (mod q)`).

use std::fmt;

use num_bigint::BigInt;
use num_traits::Zero;

use crate::filter::{MatrixLayout, SparseMatrix};
use crate::polyselect::PolyPair;
use crate::sieve::{
    classify_prime, CheckedRelation, CheckedRelations, FactorBase, PrimeClass, Relation,
    RelationChecker,
};
use crate::sqrt::PrimeSearchLimits;
use shared_numfield::{
    is_prime_u64, FpPoly, FpPolyError, IntPoly, PrimeModulus, RootEnumeration, RootSearchLimits,
};

// ─── individual root characters ──────────────────────────────────────────────

/// The value of a character on a relation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CharacterValue {
    /// `a − b·s` is a nonzero square modulo `q`; bit 0.
    Square,
    /// `a − b·s` is a nonsquare modulo `q`; bit 1.
    NonSquare,
    /// `a − b·s ≡ 0 (mod q)`: no Legendre bit exists.
    Zero,
}

impl CharacterValue {
    /// The GF(2) bit, or `None` for [`CharacterValue::Zero`].
    #[must_use]
    pub fn bit(self) -> Option<bool> {
        match self {
            Self::Square => Some(false),
            Self::NonSquare => Some(true),
            Self::Zero => None,
        }
    }
}

/// A pair `(q, s)`: an odd prime and a root of `f` modulo `q`.
///
/// Only [`RootCharacters`] mints one, after checking `q ∤ c·disc(F)` and `f(s) ≡ 0 (mod q)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RootCharacter {
    modulus: PrimeModulus,
    root: u64,
}

impl RootCharacter {
    /// The odd prime `q`.
    #[must_use]
    pub fn prime(&self) -> u64 {
        self.modulus.get()
    }

    /// The root `s` in `[0, q)`.
    #[must_use]
    pub fn root(&self) -> u64 {
        self.root
    }

    /// The Legendre symbol of `a − b·s` modulo `q`, by Euler's criterion.
    #[must_use]
    pub fn evaluate(&self, a: &BigInt, b: &BigInt) -> CharacterValue {
        let m = self.modulus;
        let value = m.sub(m.reduce_bigint(a), m.mul(m.reduce_bigint(b), self.root));
        if value == 0 {
            CharacterValue::Zero
        } else if m.pow(value, (m.get() - 1) / 2) == 1 {
            CharacterValue::Square
        } else {
            CharacterValue::NonSquare
        }
    }
}

/// Why a pair `(q, s)` is not a valid character.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairRejection {
    /// `q` is even, so the Legendre symbol does not exist.
    Even,
    /// `q` is not prime.
    NotPrime,
    /// `q | c·disc(F)`.
    Exceptional,
    /// `s ≥ q`.
    RootOutOfRange,
    /// `f(s) ≢ 0 (mod q)`.
    NotRoot,
    /// The pair repeats the one at position `first`.
    Duplicate {
        /// Position of the earlier occurrence.
        first: usize,
    },
}

/// The limit that stopped a bounded character search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CharacterResource {
    /// `max_primes_examined` primes were examined.
    PrimeBudget {
        /// The configured limit.
        limit: u64,
    },
    /// The machine-prime range ended.
    PrimeRange,
}

/// A pair that vanishes on a retained relation and was therefore rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RejectedPair {
    /// The prime `q`.
    pub prime: u64,
    /// The root `s`.
    pub root: u64,
    /// The first relation, by position, with `a − b·s ≡ 0 (mod q)`.
    pub relation: usize,
}

/// What a bounded search saw, so that a shortage is never silent.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SelectionReport {
    /// Primes examined; each counts against the prime budget.
    pub primes_examined: u64,
    /// Examined primes dividing `c·disc(F)`, ascending; they contribute no pair.
    pub exceptional_primes: Vec<u64>,
    /// Primes whose root enumeration was truncated, ascending. Roots beyond the evaluation budget
    /// were not offered, which is a resource outcome and not evidence that they do not exist.
    pub truncated_primes: Vec<u64>,
    /// Roots rejected because they vanish on a retained relation, in search order.
    pub rejected_pairs: Vec<RejectedPair>,
}

/// An error of character selection, filling, or verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RootCharacterError {
    /// The pair at `position` is not a valid character.
    Pair {
        /// Position in the input list.
        position: usize,
        /// The prime.
        prime: u64,
        /// The root.
        root: u64,
        /// The reason.
        reason: PairRejection,
    },
    /// A bounded search ended before finding `requested` pairs. A resource outcome only.
    Exhausted {
        /// The limit reached.
        resource: CharacterResource,
        /// The requested number of pairs.
        requested: usize,
        /// The number of pairs found.
        found: usize,
        /// Everything the search saw.
        report: Box<SelectionReport>,
    },
    /// The characters and the relation set are bound to different checkers.
    RelationsCheckerMismatch,
    /// The characters and the layout are bound to different checkers.
    LayoutCheckerMismatch,
    /// The layout has `layout` character slots but there are `characters` characters.
    LayoutCharacterCount {
        /// The layout's slot count.
        layout: usize,
        /// The number of characters.
        characters: usize,
    },
    /// The matrix width, obstruction block, or weight table differs from the layout.
    MatrixShapeMismatch,
    /// A relation index is not below the relation count.
    RelationOutOfRange {
        /// The index.
        index: usize,
        /// The number of relations.
        len: usize,
    },
    /// A character vanishes on a relation, so it has no bit there.
    ZeroValue {
        /// The character slot.
        slot: usize,
        /// The prime.
        prime: u64,
        /// The root.
        root: u64,
        /// The relation position.
        relation: usize,
    },
    /// A row's provenance is not strictly ascending.
    MalformedProvenance {
        /// The row.
        row: usize,
    },
    /// A row's provenance names a relation that does not exist.
    ProvenanceOutOfRange {
        /// The row.
        row: usize,
        /// The offending index.
        index: usize,
        /// The number of relations.
        relations: usize,
    },
    /// A row already has a set bit in or beyond the character columns.
    CharacterColumnNotZero {
        /// The row.
        row: usize,
        /// The offending column.
        column: usize,
    },
    /// A row's character bit differs from the XOR of its provenance's original bits.
    CharacterMismatch {
        /// The row.
        row: usize,
        /// The character slot.
        slot: usize,
    },
    /// An arithmetic error from the finite-field substrate.
    Arithmetic(FpPolyError),
}

impl fmt::Display for PairRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Even => f.write_str("the prime is even"),
            Self::NotPrime => f.write_str("the modulus is not prime"),
            Self::Exceptional => f.write_str("the prime divides c*disc(F)"),
            Self::RootOutOfRange => f.write_str("the root is not below the prime"),
            Self::NotRoot => f.write_str("f does not vanish at the root"),
            Self::Duplicate { first } => write!(f, "the pair repeats position {first}"),
        }
    }
}

impl fmt::Display for CharacterResource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PrimeBudget { limit } => write!(f, "prime budget of {limit} reached"),
            Self::PrimeRange => f.write_str("machine-prime range ended"),
        }
    }
}

impl fmt::Display for RootCharacterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pair { position, prime, root, reason } => {
                write!(f, "pair {position} = ({prime}, {root}) is invalid: {reason}")
            }
            Self::Exhausted { resource, requested, found, .. } => {
                write!(f, "found {found} of {requested} pairs: {resource}")
            }
            Self::RelationsCheckerMismatch => {
                f.write_str("characters and relations are bound to different checkers")
            }
            Self::LayoutCheckerMismatch => {
                f.write_str("characters and layout are bound to different checkers")
            }
            Self::LayoutCharacterCount { layout, characters } => {
                write!(f, "layout has {layout} character slots but there are {characters}")
            }
            Self::MatrixShapeMismatch => f.write_str("matrix shape differs from the layout"),
            Self::RelationOutOfRange { index, len } => {
                write!(f, "relation {index} is out of range for {len} relations")
            }
            Self::ZeroValue { slot, prime, root, relation } => {
                write!(f, "character {slot} = ({prime}, {root}) vanishes on relation {relation}")
            }
            Self::MalformedProvenance { row } => {
                write!(f, "row {row} has provenance that is not strictly ascending")
            }
            Self::ProvenanceOutOfRange { row, index, relations } => {
                write!(f, "row {row} names relation {index} of {relations}")
            }
            Self::CharacterColumnNotZero { row, column } => {
                write!(f, "row {row} already has a bit at column {column}")
            }
            Self::CharacterMismatch { row, slot } => {
                write!(f, "row {row} character {slot} differs from its provenance XOR")
            }
            Self::Arithmetic(e) => write!(f, "finite-field arithmetic: {e}"),
        }
    }
}

impl std::error::Error for PairRejection {}
impl std::error::Error for CharacterResource {}
impl std::error::Error for RootCharacterError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Arithmetic(e) => Some(e),
            _ => None,
        }
    }
}

/// Budgets for [`select_root_characters`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RootCharacterLimits {
    /// Ascending prime search; the first candidate is raised to 3 because 2 has no character.
    pub prime_search: PrimeSearchLimits,
    /// Exhaustive root enumeration in each prime field.
    pub root_search: RootSearchLimits,
}

impl Default for RootCharacterLimits {
    /// Start at 3, examine at most 1000 primes, and evaluate at most 100 000 residues per prime.
    fn default() -> Self {
        Self {
            prime_search: PrimeSearchLimits { first_candidate: 3, max_primes_examined: 1_000 },
            root_search: RootSearchLimits { max_evaluations: 100_000 },
        }
    }
}

/// An ordered list of characters bound to the checker that defines `f`, `c`, and `disc(F)`.
///
/// Slot `k` is `characters()[k]`. The bound checker lets filling and verification refuse a
/// relation set or layout that belongs to a different context by equality, not by convention.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootCharacters {
    checker: RelationChecker,
    characters: Vec<RootCharacter>,
}

/// The result of a complete bounded search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootCharacterSelection {
    characters: RootCharacters,
    report: SelectionReport,
}

impl RootCharacterSelection {
    /// The selected characters.
    #[must_use]
    pub fn characters(&self) -> &RootCharacters {
        &self.characters
    }

    /// What the search saw, including rejected pairs and truncated primes.
    #[must_use]
    pub fn report(&self) -> &SelectionReport {
        &self.report
    }

    /// Split into the characters and the report.
    #[must_use]
    pub fn into_parts(self) -> (RootCharacters, SelectionReport) {
        (self.characters, self.report)
    }
}

impl RootCharacters {
    /// Validate explicit pairs `(q, s)` against `checker`'s context.
    ///
    /// # Errors
    ///
    /// [`RootCharacterError::Pair`] for the first invalid pair, with the first failing reason in
    /// the order: [`PairRejection::Even`], [`PairRejection::NotPrime`],
    /// [`PairRejection::Exceptional`], [`PairRejection::RootOutOfRange`],
    /// [`PairRejection::NotRoot`], [`PairRejection::Duplicate`].
    pub fn new(
        checker: &RelationChecker,
        pairs: &[(u64, u64)],
    ) -> Result<Self, RootCharacterError> {
        let mut characters: Vec<RootCharacter> = Vec::with_capacity(pairs.len());
        for (position, &(prime, root)) in pairs.iter().enumerate() {
            let reject = |reason| RootCharacterError::Pair { position, prime, root, reason };
            if prime % 2 == 0 {
                return Err(reject(PairRejection::Even));
            }
            if !is_prime_u64(prime) {
                return Err(reject(PairRejection::NotPrime));
            }
            if classify_prime(checker.context(), prime) == PrimeClass::Exceptional {
                return Err(reject(PairRejection::Exceptional));
            }
            if root >= prime {
                return Err(reject(PairRejection::RootOutOfRange));
            }
            let modulus = PrimeModulus::new(prime).map_err(RootCharacterError::Arithmetic)?;
            let reduced = FpPoly::from_int_poly(checker.context().algebraic_poly(), modulus);
            if reduced.eval(root) != 0 {
                return Err(reject(PairRejection::NotRoot));
            }
            if let Some(first) = characters.iter().position(|c| c.prime() == prime && c.root == root)
            {
                return Err(reject(PairRejection::Duplicate { first }));
            }
            characters.push(RootCharacter { modulus, root });
        }
        Ok(Self { checker: checker.clone(), characters })
    }

    /// The bound checker.
    #[must_use]
    pub fn checker(&self) -> &RelationChecker {
        &self.checker
    }

    /// The characters in slot order.
    #[must_use]
    pub fn characters(&self) -> &[RootCharacter] {
        &self.characters
    }

    /// The number of characters, counting pairs and not primes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.characters.len()
    }

    /// True iff there are no characters.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.characters.is_empty()
    }

    /// The bits of one original relation, one per slot.
    ///
    /// # Errors
    ///
    /// [`RootCharacterError::RelationsCheckerMismatch`] when `relations` is bound to another
    /// checker; [`RootCharacterError::RelationOutOfRange`] when `index` is not a relation;
    /// [`RootCharacterError::ZeroValue`] when a character vanishes on the relation.
    pub fn relation_bits(
        &self,
        relations: &CheckedRelations,
        index: usize,
    ) -> Result<Vec<bool>, RootCharacterError> {
        self.check_relations(relations)?;
        let relation = relations
            .relations()
            .get(index)
            .ok_or(RootCharacterError::RelationOutOfRange { index, len: relations.len() })?;
        self.bits_of(relation, index)
    }

    /// The bits of every original relation: `table[i][slot]`.
    ///
    /// # Errors
    ///
    /// As [`RootCharacters::relation_bits`], for the first failing relation.
    pub fn bit_table(
        &self,
        relations: &CheckedRelations,
    ) -> Result<Vec<Vec<bool>>, RootCharacterError> {
        self.check_relations(relations)?;
        relations.relations().iter().enumerate().map(|(i, r)| self.bits_of(r, i)).collect()
    }

    /// Set the character columns of `matrix` from the original relations.
    ///
    /// The bit of row `r` at slot `k` is the XOR over `r.provenance` of the original relations'
    /// bits, so it equals what filling the freshly built matrix would have produced and then
    /// carrying through elimination by XOR. The matrix is unchanged on error. Character columns
    /// are the last columns, so setting them keeps each row's columns sorted.
    ///
    /// # Errors
    ///
    /// Binding: [`RootCharacterError::RelationsCheckerMismatch`],
    /// [`RootCharacterError::LayoutCheckerMismatch`],
    /// [`RootCharacterError::LayoutCharacterCount`], [`RootCharacterError::MatrixShapeMismatch`].
    /// Data: [`RootCharacterError::ZeroValue`] when any character vanishes on any original
    /// relation; [`RootCharacterError::MalformedProvenance`] and
    /// [`RootCharacterError::ProvenanceOutOfRange`] for a row; and
    /// [`RootCharacterError::CharacterColumnNotZero`] when a row already has a bit in or beyond
    /// the character columns, so that filling twice is refused.
    pub fn fill_columns(
        &self,
        matrix: &mut SparseMatrix,
        layout: &MatrixLayout,
        relations: &CheckedRelations,
    ) -> Result<(), RootCharacterError> {
        let table = self.checked_table(matrix, layout, relations)?;
        let first = layout.sign_column() + 1;
        let mut additions: Vec<Vec<usize>> = Vec::with_capacity(matrix.rows.len());
        for (row_index, row) in matrix.rows.iter().enumerate() {
            if let Some(&column) = row.cols.iter().find(|&&c| c >= first) {
                return Err(RootCharacterError::CharacterColumnNotZero { row: row_index, column });
            }
            let bits = xor_provenance(&table, row_index, &row.provenance, self.len())?;
            additions.push(
                (0..self.len()).filter(|&slot| bits[slot]).map(|slot| first + slot).collect(),
            );
        }
        for (row, added) in matrix.rows.iter_mut().zip(additions) {
            for column in added {
                row.cols.push(column);
                matrix.col_weights[column] += 1;
            }
        }
        Ok(())
    }

    /// Check that every row's character columns equal the XOR of its provenance's original bits.
    ///
    /// This is the provenance identity that makes filling after filtering equivalent to filling
    /// before it.
    ///
    /// # Errors
    ///
    /// The binding and data errors of [`RootCharacters::fill_columns`] except
    /// `CharacterColumnNotZero`, and [`RootCharacterError::CharacterMismatch`] for the first row
    /// and slot that differ.
    pub fn verify_columns(
        &self,
        matrix: &SparseMatrix,
        layout: &MatrixLayout,
        relations: &CheckedRelations,
    ) -> Result<(), RootCharacterError> {
        let table = self.checked_table(matrix, layout, relations)?;
        let first = layout.sign_column() + 1;
        for (row_index, row) in matrix.rows.iter().enumerate() {
            let bits = xor_provenance(&table, row_index, &row.provenance, self.len())?;
            for (slot, &expected) in bits.iter().enumerate() {
                if row.cols.contains(&(first + slot)) != expected {
                    return Err(RootCharacterError::CharacterMismatch { row: row_index, slot });
                }
            }
        }
        Ok(())
    }

    fn check_relations(&self, relations: &CheckedRelations) -> Result<(), RootCharacterError> {
        if relations.checker() == &self.checker {
            Ok(())
        } else {
            Err(RootCharacterError::RelationsCheckerMismatch)
        }
    }

    /// Validate the binding of characters, relations, layout, and matrix, then tabulate bits.
    fn checked_table(
        &self,
        matrix: &SparseMatrix,
        layout: &MatrixLayout,
        relations: &CheckedRelations,
    ) -> Result<Vec<Vec<bool>>, RootCharacterError> {
        self.check_relations(relations)?;
        if layout.checker() != &self.checker {
            return Err(RootCharacterError::LayoutCheckerMismatch);
        }
        if layout.character_count() != self.len() {
            return Err(RootCharacterError::LayoutCharacterCount {
                layout: layout.character_count(),
                characters: self.len(),
            });
        }
        if !layout.shape_matches(matrix) {
            return Err(RootCharacterError::MatrixShapeMismatch);
        }
        self.bit_table(relations)
    }

    fn bits_of(
        &self,
        relation: &CheckedRelation,
        index: usize,
    ) -> Result<Vec<bool>, RootCharacterError> {
        self.characters
            .iter()
            .enumerate()
            .map(|(slot, c)| {
                c.evaluate(relation.a(), relation.b()).bit().ok_or(RootCharacterError::ZeroValue {
                    slot,
                    prime: c.prime(),
                    root: c.root,
                    relation: index,
                })
            })
            .collect()
    }
}

/// XOR of `table[i]` over the strictly ascending, in-range `provenance` of row `row`.
fn xor_provenance(
    table: &[Vec<bool>],
    row: usize,
    provenance: &[usize],
    width: usize,
) -> Result<Vec<bool>, RootCharacterError> {
    if provenance.windows(2).any(|w| w[0] >= w[1]) {
        return Err(RootCharacterError::MalformedProvenance { row });
    }
    let mut bits = vec![false; width];
    for &index in provenance {
        let original = table.get(index).ok_or(RootCharacterError::ProvenanceOutOfRange {
            row,
            index,
            relations: table.len(),
        })?;
        for (bit, &original_bit) in bits.iter_mut().zip(original) {
            *bit ^= original_bit;
        }
    }
    Ok(bits)
}

/// Select `count` valid pairs `(q, s)` by an ascending, bounded search over primes.
///
/// For each prime `q` from `limits.prime_search.first_candidate` (raised to 3), in order:
/// `q | c·disc(F)` is skipped and listed in the report; otherwise the roots of `f` modulo `q` are
/// enumerated ascending under `limits.root_search`, and each root `s` that makes
/// `a − b·s ≢ 0 (mod q)` for **every** relation of `relations` becomes the next character. A
/// vanishing root is reported in [`SelectionReport::rejected_pairs`]. A prime can supply several
/// pairs, and the search stops as soon as `count` are found, possibly in the middle of a prime.
/// A truncated enumeration still offers the roots it found, since each is a proven root; the prime
/// is listed in [`SelectionReport::truncated_primes`].
///
/// # Errors
///
/// [`RootCharacterError::Exhausted`] when `max_primes_examined` primes are used up
/// ([`CharacterResource::PrimeBudget`]) or the machine-prime range ends
/// ([`CharacterResource::PrimeRange`]) before `count` pairs are found. This is a resource outcome
/// and says nothing about whether more valid pairs exist. [`RootCharacterError::Arithmetic`] is
/// returned only if the finite-field substrate reports an impossible condition.
pub fn select_root_characters(
    relations: &CheckedRelations,
    count: usize,
    limits: RootCharacterLimits,
) -> Result<RootCharacterSelection, RootCharacterError> {
    let mut chosen: Vec<RootCharacter> = Vec::with_capacity(count);
    let mut report = SelectionReport::default();
    let mut cursor = Some(limits.prime_search.first_candidate.max(3));
    while chosen.len() < count {
        if report.primes_examined >= limits.prime_search.max_primes_examined {
            let resource =
                CharacterResource::PrimeBudget { limit: limits.prime_search.max_primes_examined };
            return Err(exhausted(resource, count, &chosen, report));
        }
        let Some(prime) = cursor.and_then(next_prime_from) else {
            return Err(exhausted(CharacterResource::PrimeRange, count, &chosen, report));
        };
        cursor = prime.checked_add(1);
        report.primes_examined += 1;
        scan_prime(relations, prime, limits.root_search, count, &mut chosen, &mut report)?;
    }
    Ok(RootCharacterSelection {
        characters: RootCharacters {
            checker: relations.checker().clone(),
            characters: chosen,
        },
        report,
    })
}

fn exhausted(
    resource: CharacterResource,
    requested: usize,
    chosen: &[RootCharacter],
    report: SelectionReport,
) -> RootCharacterError {
    RootCharacterError::Exhausted {
        resource,
        requested,
        found: chosen.len(),
        report: Box::new(report),
    }
}

/// The least prime `>= start`, or `None` when the machine range ends first.
fn next_prime_from(start: u64) -> Option<u64> {
    let mut candidate = start;
    loop {
        if is_prime_u64(candidate) {
            return Some(candidate);
        }
        candidate = candidate.checked_add(1)?;
    }
}

/// Offer the valid roots modulo odd prime `prime` to `chosen` until it holds `count` pairs.
fn scan_prime(
    relations: &CheckedRelations,
    prime: u64,
    root_search: RootSearchLimits,
    count: usize,
    chosen: &mut Vec<RootCharacter>,
    report: &mut SelectionReport,
) -> Result<(), RootCharacterError> {
    let context = relations.context();
    if classify_prime(context, prime) == PrimeClass::Exceptional {
        report.exceptional_primes.push(prime);
        return Ok(());
    }
    let modulus = PrimeModulus::new(prime).map_err(RootCharacterError::Arithmetic)?;
    let reduced = FpPoly::from_int_poly(context.algebraic_poly(), modulus);
    let roots = match reduced.roots(root_search).map_err(RootCharacterError::Arithmetic)? {
        RootEnumeration::Complete { roots, .. } => roots,
        RootEnumeration::Truncated { found, .. } => {
            report.truncated_primes.push(prime);
            found
        }
    };
    for root in roots {
        if chosen.len() == count {
            break;
        }
        let candidate = RootCharacter { modulus, root };
        let vanishing = relations
            .relations()
            .iter()
            .position(|r| candidate.evaluate(r.a(), r.b()) == CharacterValue::Zero);
        match vanishing {
            Some(relation) => report.rejected_pairs.push(RejectedPair { prime, root, relation }),
            None => chosen.push(candidate),
        }
    }
    Ok(())
}

// ─── legacy norm characters ──────────────────────────────────────────────────

// ─── DEFAULT_NUM_QC ───────────────────────────────────────────────────────────

/// Default number of quadratic-character columns (demonstration fidelity).
///
/// Principle-4 annotation: at NFS scale this is typically 20–50; at toy scale 10 suffices.
pub const DEFAULT_NUM_QC: usize = 10;

// ─── populate_qc_columns ─────────────────────────────────────────────────────

/// Populate legacy norm-character columns in a matrix.
///
/// These are Legendre symbols of algebraic *norms*, not individual root characters; see
/// [`RootCharacters::fill_columns`] for the latter. They do not guarantee that a square root
/// exists.
///
/// For each row `i`, for each QC prime `qc_primes[k]`, sets column
/// `matrix.obstruction_col_start + 1 + k` to the Legendre-symbol parity of the
/// algebraic norm of `relations[provenance[0]]` at `qc_primes[k]`.
///
/// Preconditions:
/// - `matrix` was built with `fb.obstruction_count = 1 + qc_primes.len()`.
/// - `qc_primes` are primes > `fb.b_alg` that split completely in K.
/// - Each row's provenance is non-empty (true for any matrix from `build_matrix`).
///
/// # Legendre-symbol computation
///
/// For a relation with algebraic norm `N_alg(a, b)` and a QC prime `q`, the QC parity is
/// `(N_alg / q) mod 2` where `(· / q)` is the Legendre symbol. If `N_alg ≡ 0 (mod q)`,
/// the symbol is 0 (even parity). The parity is 1 iff the Legendre symbol is −1.
///
/// For merged rows (provenance has multiple original relations), the QC parity is the
/// XOR (GF(2) sum) of the individual relations' parities — consistent with how merged
/// rows' factor-base columns are the XOR of the originals.
pub fn populate_qc_columns(
    matrix: &mut SparseMatrix,
    relations: &[Relation],
    fb: &FactorBase,
    poly: &PolyPair,
    qc_primes: &[u64],
) {
    let num_qc = qc_primes.len();
    debug_assert_eq!(
        fb.obstruction_count,
        1 + num_qc,
        "populate_qc_columns: fb.obstruction_count must be 1 + qc_primes.len()"
    );

    for row in matrix.rows.iter_mut() {
        // Compute the QC parity for each QC prime by XOR-ing over all provenance relations.
        for (k, &q) in qc_primes.iter().enumerate() {
            let col = matrix.obstruction_col_start + 1 + k;
            // XOR the parities of all provenance relations.
            let mut parity = false;
            for &rel_idx in &row.provenance {
                let rel = &relations[rel_idx];
                let norm = algebraic_norm_mod_q(&rel.a, &rel.b, &poly.f, q);
                parity ^= legendre_parity(norm, q);
            }
            // Update the column set: insert or remove the column.
            if parity {
                // Insert col into row.cols (maintaining sorted order).
                match row.cols.binary_search(&col) {
                    Ok(_) => {} // already present (shouldn't happen for a fresh matrix)
                    Err(pos) => {
                        row.cols.insert(pos, col);
                        matrix.col_weights[col] += 1;
                    }
                }
            }
            // If parity is false, the column is already absent (zero); nothing to do.
        }
    }
}

// ─── select_qc_primes ────────────────────────────────────────────────────────

/// Select `num_qc` auxiliary primes for legacy norm-character columns.
///
/// Not the root-character contract: [`select_root_characters`] chooses pairs `(q, s)`. Returns the first `num_qc` primes `q > b_alg` such that `f(x)` has `deg(f)` distinct
/// roots mod `q` (i.e., `q` splits completely in K = ℚ[x]/(f)).
pub fn select_qc_primes(f: &IntPoly, b_alg: u64, num_qc: usize) -> Vec<u64> {
    let deg = match f.degree() {
        None => return Vec::new(),
        Some(0) => return Vec::new(),
        Some(d) => d,
    };

    let mut result = Vec::with_capacity(num_qc);
    let mut candidate = b_alg + 1;

    while result.len() < num_qc {
        if is_prime(candidate) && splits_completely(f, candidate, deg) {
            result.push(candidate);
        }
        candidate += 1;
    }

    result
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Compute `N_alg(a, b) mod q` as a `u64` in `[0, q)`.
///
/// Uses the homogeneous form: `Σ f.coeffs[i] · a^i · b^{d-i} mod q`.
/// All arithmetic is done mod `q` to avoid large intermediate values.
fn algebraic_norm_mod_q(a: &BigInt, b: &BigInt, f: &IntPoly, q: u64) -> u64 {
    let d = match f.degree() {
        None => return 0,
        Some(d) => d,
    };

    let q_big = BigInt::from(q);
    let a_mod = mod_reduce_u64(a, q);
    let b_mod = mod_reduce_u64(b, q);

    // Compute Σ f.coeffs[i] · a^i · b^{d-i} mod q using direct summation.
    let mut result = 0u64;
    let mut a_pow = 1u64; // a^i mod q
    // Precompute b^{d-i} for i from 0 to d: b^d, b^{d-1}, ..., b^0.
    // We compute b_pow_d_minus_i = b^{d-i} mod q.
    // Start with b^d and divide by b each step — but division mod q requires modular inverse.
    // Instead, precompute all b powers: b_pows[j] = b^j mod q.
    let mut b_pows = vec![0u64; d + 1];
    b_pows[0] = 1;
    for j in 1..=d {
        b_pows[j] = mul_mod(b_pows[j - 1], b_mod, q);
    }

    for i in 0..=d {
        let coeff_big = f.coeffs.get(i).cloned().unwrap_or_else(BigInt::zero);
        if !coeff_big.is_zero() {
            let coeff_mod = mod_reduce_u64_bigint(&coeff_big, q);
            let b_pow = b_pows[d - i];
            let term = mul_mod(mul_mod(coeff_mod, a_pow, q), b_pow, q);
            result = add_mod(result, term, q);
        }
        if i < d {
            a_pow = mul_mod(a_pow, a_mod, q);
        }
    }

    // Suppress unused variable warning for q_big.
    let _ = q_big;
    result
}

/// Compute the Legendre-symbol parity: true iff `(norm / q) = -1`.
///
/// Uses Euler's criterion: `norm^((q-1)/2) mod q`.
/// - If result is 0: symbol is 0 (parity = false).
/// - If result is 1: symbol is 1 (parity = false).
/// - If result is q-1: symbol is -1 (parity = true).
fn legendre_parity(norm: u64, q: u64) -> bool {
    if norm == 0 {
        return false; // q | norm → symbol is 0
    }
    let exp = (q - 1) / 2;
    let result = pow_mod(norm, exp, q);
    result == q - 1
}

/// Check if `q` is prime using trial division.
fn is_prime(q: u64) -> bool {
    if q < 2 {
        return false;
    }
    if q == 2 {
        return true;
    }
    if q % 2 == 0 {
        return false;
    }
    let mut i = 3u64;
    while i * i <= q {
        if q % i == 0 {
            return false;
        }
        i += 2;
    }
    true
}

/// Check if `f` splits completely mod `q`: has exactly `deg(f)` distinct roots.
///
/// A prime `q` splits completely in K = ℚ[x]/(f) iff `f mod q` has `deg(f)` distinct
/// linear factors. We check this by counting distinct roots of `f` in `[0, q)`.
fn splits_completely(f: &IntPoly, q: u64, deg: usize) -> bool {
    let mut root_count = 0usize;
    for r in 0..q {
        let val = eval_poly_mod(f, r, q);
        if val == 0 {
            root_count += 1;
        }
    }
    root_count == deg
}

/// Evaluate `f(r) mod q` for `r` in `[0, q)`.
fn eval_poly_mod(f: &IntPoly, r: u64, q: u64) -> u64 {
    // Horner's method mod q.
    let mut result = 0u64;
    for c in f.coeffs.iter().rev() {
        result = mul_mod(result, r, q);
        let c_mod = mod_reduce_u64_bigint(c, q);
        result = add_mod(result, c_mod, q);
    }
    result
}

/// Reduce a `BigInt` to `[0, q)` as a `u64`.
fn mod_reduce_u64(a: &BigInt, q: u64) -> u64 {
    let q_big = BigInt::from(q);
    let r = a % &q_big;
    let r = if r < BigInt::zero() { r + &q_big } else { r };
    // Convert to u64: the result is in [0, q), which fits in u64.
    use num_traits::ToPrimitive;
    r.to_u64().unwrap_or(0)
}

/// Reduce a `BigInt` coefficient to `[0, q)` as a `u64`.
fn mod_reduce_u64_bigint(a: &BigInt, q: u64) -> u64 {
    mod_reduce_u64(a, q)
}

/// Multiply `a * b mod q` without overflow using u128.
fn mul_mod(a: u64, b: u64, q: u64) -> u64 {
    ((a as u128 * b as u128) % q as u128) as u64
}

/// Add `a + b mod q`.
fn add_mod(a: u64, b: u64, q: u64) -> u64 {
    let s = a + b;
    if s >= q { s - q } else { s }
}

/// Compute `base^exp mod q` using fast exponentiation.
fn pow_mod(mut base: u64, mut exp: u64, q: u64) -> u64 {
    if q == 1 {
        return 0;
    }
    let mut result = 1u64;
    base %= q;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, q);
        }
        base = mul_mod(base, base, q);
        exp >>= 1;
    }
    result
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    /// f(x) = x³ − x − 1.
    fn f_cubic() -> IntPoly {
        IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)])
    }

    #[test]
    fn legendre_parity_known_values() {
        // (2/7): 2^3 = 8 ≡ 1 (mod 7) → symbol = 1 → parity = false.
        assert!(!legendre_parity(2, 7));
        // (3/7): 3^3 = 27 ≡ 6 ≡ -1 (mod 7) → symbol = -1 → parity = true.
        assert!(legendre_parity(3, 7));
        // (0/7): symbol = 0 → parity = false.
        assert!(!legendre_parity(0, 7));
    }

    #[test]
    fn next_prime_from_returns_the_least_prime_at_or_above() {
        assert_eq!(next_prime_from(14), Some(17));
        assert_eq!(next_prime_from(17), Some(17));
        // 2^64 − 59 is the largest prime below 2^64.
        assert_eq!(next_prime_from(u64::MAX - 58), Some(u64::MAX - 58));
        assert_eq!(next_prime_from(u64::MAX - 57), None);
    }

    #[test]
    fn pow_mod_correctness() {
        assert_eq!(pow_mod(2, 10, 1000), 24); // 2^10 = 1024, 1024 mod 1000 = 24
        assert_eq!(pow_mod(3, 0, 7), 1);
        assert_eq!(pow_mod(0, 5, 7), 0);
    }

    #[test]
    fn is_prime_small() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(!is_prime(9));
        assert!(is_prime(97));
    }

    #[test]
    fn splits_completely_cubic() {
        // f(x) = x³ − x − 1. Check a few primes.
        let f = f_cubic();
        // q = 23: f splits completely mod 23 iff it has 3 distinct roots.
        // f(x) = x³ - x - 1 mod 23: check by brute force.
        let roots_23: Vec<u64> = (0..23).filter(|&r| eval_poly_mod(&f, r, 23) == 0).collect();
        let sc_23 = splits_completely(&f, 23, 3);
        assert_eq!(sc_23, roots_23.len() == 3);
    }

    #[test]
    fn select_qc_primes_returns_correct_count() {
        let f = f_cubic();
        let primes = select_qc_primes(&f, 13, 3);
        assert_eq!(primes.len(), 3, "should return exactly 3 primes");
        // All returned primes should be > 13.
        for &p in &primes {
            assert!(p > 13, "all QC primes should be > b_alg");
        }
        // All returned primes should split completely.
        for &p in &primes {
            assert!(splits_completely(&f, p, 3), "prime {p} should split completely");
        }
    }

    #[test]
    fn algebraic_norm_mod_q_matches_direct() {
        // N_alg(2, 1) for f(x) = x³ - x - 1 = 8 - 2 - 1 = 5.
        let f = f_cubic();
        let norm = algebraic_norm_mod_q(&bi(2), &bi(1), &f, 7);
        assert_eq!(norm, 5 % 7, "N_alg(2,1) mod 7 should be 5");
    }
}
