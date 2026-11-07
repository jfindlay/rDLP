//! Checked dependencies: a candidate set of filtered-matrix rows validated against the original
//! relations, with the exact rational and algebraic products of the relations it lifts to.
//!
//! A kernel vector of the filtered matrix is only a *candidate*. Whatever produced it, nothing in
//! the matrix proves that its rows stand for the relations it claims. [`DependencyContext`] binds
//! the matrix to the checked relations, matrix layout, and root characters it was built from, and
//! checks the identity that makes lifting sound: each matrix row equals the GF(2) sum of the
//! original rows its provenance names, in every column including the character columns. A
//! [`CheckedDependency`] is then minted only by [`DependencyContext::validate`], which re-derives
//! everything from the selected row indices:
//!
//! 1. the selection is non-empty, strictly ascending, and in range;
//! 2. the selected rows sum to zero in the context's own matrix;
//! 3. the provenance of the selected rows is lifted by symmetric difference, so a relation that
//!    two rows share cancels instead of being counted twice;
//! 4. the parity of the lifted relations is recomputed from the original rows, not read from the
//!    matrix.
//!
//! # Exact products
//!
//! For the `k` relations `(a_i, b_i)` of a dependency, the rational product is
//! `P = Π(a_i − b_i·m)` over `Z` with the exact `m`, and
//! `γ = Π(c·a_i − b_i·θ) / c^k` is the product `Π(a_i − b_i·α)` in the `θ` basis of
//! `A = Q[t]/(F)`. Neither is a certificate of squareness. The parity checks guarantee that every
//! rational factor-base column and the sign have even weight, so `P` is an exact integer square.
//! The good-ideal and character parities are necessary filters only, so `γ` need not be a square
//! and exact reconstruction stays authoritative. For a non-monic `f` the raw homogeneous norm
//! product `Π f_hom(a_i, b_i) = c^k·Norm(γ)` need not be a square even when `γ` is, so nothing here
//! forms or tests it.
//!
//! # Duplicates
//!
//! [`CheckedRelations`] rejects a repeated `(a, b)` pair. A caller that may hold duplicates removes
//! them from the relation list *before* building the checked relations and the matrix, keeping
//! first occurrences, so provenance indexes the list that was checked. Nothing renumbers
//! provenance after validation.

use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::pow::pow;
use shared_numfield::{IntPoly, RatPoly};

use crate::filter::{MatrixLayout, SparseMatrix};
use crate::linalg::{KernelVector, RootCharacterError, RootCharacters};
use crate::sieve::CheckedRelations;

use super::{
    FailureClass, InvariantViolation, NormalizedTarget, PrecisionLimits, RelationFactors,
    SqrtContext, SqrtError, TargetOutcome,
};

// ─── errors ──────────────────────────────────────────────────────────────────

/// Why a relation set, matrix, layout, and characters do not form a [`DependencyContext`].
///
/// Checks run in the order of the variants, and the first failure is reported. The per-row
/// variants are checked row by row in ascending row order, each row in variant order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DependencyContextError {
    /// The layout and the relations are bound to different checkers.
    LayoutRelationsMismatch,
    /// The root characters and the relations are bound to different checkers.
    CharactersRelationsMismatch,
    /// The layout has `layout` character slots but there are `characters` characters.
    CharacterCountMismatch {
        /// The layout's slot count.
        layout: usize,
        /// The number of characters.
        characters: usize,
    },
    /// The matrix width, obstruction block, or weight table differs from the layout.
    MatrixShapeMismatch,
    /// A character has no bit on some original relation (it vanishes there), or is otherwise
    /// unusable for the relations.
    Characters(RootCharacterError),
    /// A row's columns are not strictly ascending.
    MalformedColumns {
        /// The row.
        row: usize,
    },
    /// A row names a column at or beyond the matrix width.
    ColumnOutOfRange {
        /// The row.
        row: usize,
        /// The first offending column.
        column: usize,
    },
    /// A row has no provenance. A row derived by XOR elimination and row drops always has some.
    EmptyProvenance {
        /// The row.
        row: usize,
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
        /// The first offending relation index.
        index: usize,
        /// The number of relations.
        relations: usize,
    },
    /// A row differs from the XOR of the original rows its provenance names.
    LiftedParityMismatch {
        /// The row.
        row: usize,
        /// The smallest column at which they differ.
        column: usize,
    },
    /// A recorded column weight differs from the number of rows that set the column.
    ColumnWeightMismatch {
        /// The smallest such column.
        column: usize,
    },
}

impl fmt::Display for DependencyContextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LayoutRelationsMismatch => {
                f.write_str("matrix layout and relations are bound to different checkers")
            }
            Self::CharactersRelationsMismatch => {
                f.write_str("root characters and relations are bound to different checkers")
            }
            Self::CharacterCountMismatch { layout, characters } => {
                write!(f, "layout has {layout} character slots but there are {characters}")
            }
            Self::MatrixShapeMismatch => f.write_str("matrix shape differs from the layout"),
            Self::Characters(e) => write!(f, "root characters: {e}"),
            Self::MalformedColumns { row } => {
                write!(f, "row {row} has columns that are not strictly ascending")
            }
            Self::ColumnOutOfRange { row, column } => {
                write!(f, "row {row} names column {column} beyond the matrix width")
            }
            Self::EmptyProvenance { row } => write!(f, "row {row} has empty provenance"),
            Self::MalformedProvenance { row } => {
                write!(f, "row {row} has provenance that is not strictly ascending")
            }
            Self::ProvenanceOutOfRange { row, index, relations } => {
                write!(f, "row {row} names relation {index} of {relations}")
            }
            Self::LiftedParityMismatch { row, column } => {
                write!(f, "row {row} differs from the XOR of its original rows at column {column}")
            }
            Self::ColumnWeightMismatch { column } => {
                write!(f, "column {column} has a recorded weight that is not its true weight")
            }
        }
    }
}

impl std::error::Error for DependencyContextError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Characters(e) => Some(e),
            _ => None,
        }
    }
}

/// Why a selection of rows is not a dependency of the context's matrix.
///
/// Checks run in the order of the variants, and the first failure is reported. Each is a
/// malformed candidate, not a statement about the relations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DependencyRejection {
    /// The selection is empty.
    EmptySelection,
    /// The selection is not strictly ascending: `selection[position] <= selection[position - 1]`.
    MalformedSelection {
        /// The first position that is not larger than its predecessor.
        position: usize,
    },
    /// A selected row is not a row of the matrix.
    RowOutOfRange {
        /// The first offending position in the selection.
        position: usize,
        /// The row index at that position.
        row: usize,
        /// The number of rows.
        rows: usize,
    },
    /// The selected rows do not sum to zero in the matrix.
    NotInKernel {
        /// The smallest column with odd weight over the selected rows.
        column: usize,
    },
    /// The provenance of the selected rows cancels to the empty set of relations.
    EmptyLift,
}

impl fmt::Display for DependencyRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySelection => f.write_str("selection is empty"),
            Self::MalformedSelection { position } => {
                write!(f, "selection position {position} is not larger than its predecessor")
            }
            Self::RowOutOfRange { position, row, rows } => {
                write!(f, "selection position {position} names row {row} of {rows}")
            }
            Self::NotInKernel { column } => {
                write!(f, "selected rows have odd weight at column {column}")
            }
            Self::EmptyLift => f.write_str("selected provenance cancels to no relations"),
        }
    }
}

impl std::error::Error for DependencyRejection {}

/// The failure of [`DependencyContext::validate`]: a rejected candidate, or a broken invariant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DependencyError {
    /// The candidate is not a dependency of the matrix.
    Rejected(DependencyRejection),
    /// An identity that the context guarantees was observed to fail; a defect, not a rejection.
    Invariant(InvariantViolation),
}

impl DependencyError {
    /// The failure class: [`FailureClass::InvalidInput`] for a rejection and
    /// [`FailureClass::InvariantViolation`] for a broken invariant.
    #[must_use]
    pub fn class(&self) -> FailureClass {
        match self {
            Self::Rejected(_) => FailureClass::InvalidInput,
            Self::Invariant(_) => FailureClass::InvariantViolation,
        }
    }
}

impl fmt::Display for DependencyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Rejected(e) => write!(f, "candidate rejected: {e}"),
            Self::Invariant(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DependencyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Rejected(e) => Some(e),
            Self::Invariant(e) => Some(e),
        }
    }
}

fn invariant(what: &'static str) -> DependencyError {
    DependencyError::Invariant(InvariantViolation { what })
}

// ─── the context ─────────────────────────────────────────────────────────────

/// A filtered matrix bound to the checked relations, layout, and characters it came from.
///
/// Relations, layout, and characters belong together by checker equality and by the character
/// count. The matrix belongs to them by the per-row identity checked in
/// [`DependencyContext::new`]: every row equals the XOR of its provenance's original rows, where an
/// original row is the layout's columns for the relation plus one column per set character bit.
/// Relation sets that share a checker can produce equal layouts, so that identity, not the layout,
/// ties the matrix to *these* relations. A mismatch that slipped through could only yield a
/// candidate that exact verification rejects, never a false root.
#[derive(Debug, Clone)]
pub struct DependencyContext {
    relations: CheckedRelations,
    layout: MatrixLayout,
    characters: RootCharacters,
    matrix: SparseMatrix,
    /// Per original relation, its ascending columns including the character columns.
    original_columns: Vec<Vec<usize>>,
}

impl DependencyContext {
    /// Bind `matrix` to the relations, layout, and characters it was derived from.
    ///
    /// The original row of relation `i` is `layout.row_columns(relation i)` together with
    /// `layout.character_column(slot)` for every set bit of `characters.relation_bits(relations,
    /// i)`. Each matrix row's columns must equal the XOR of the original rows over its
    /// provenance. This subsumes the character check of [`RootCharacters::verify_columns`]. With
    /// characters it rejects unfilled character columns, and it forces the columns of omitted
    /// exceptional ideals to stay zero. Filling characters before or after filtering both pass.
    ///
    /// Empty provenance is rejected: elimination and row drops keep provenance vectors linearly
    /// independent, so a genuine row never has none. A duplicated genuine row passes here and can
    /// only produce [`DependencyRejection::EmptyLift`].
    ///
    /// # Errors
    ///
    /// [`DependencyContextError`], the first failure in this order: the checker of the layout
    /// ([`DependencyContextError::LayoutRelationsMismatch`]) and of the characters
    /// ([`DependencyContextError::CharactersRelationsMismatch`]) differ from the relations'; the
    /// character count differs from the layout's; the matrix shape differs from the layout;
    /// [`DependencyContextError::Characters`] when a character vanishes on a relation; then, per
    /// row in ascending order, malformed or out-of-range columns, empty, malformed or out-of-range
    /// provenance, and a lifted-parity mismatch; and finally a column-weight mismatch.
    pub fn new(
        relations: CheckedRelations,
        matrix: SparseMatrix,
        layout: MatrixLayout,
        characters: RootCharacters,
    ) -> Result<Self, DependencyContextError> {
        check_binding(&relations, &matrix, &layout, &characters)?;
        let bits =
            characters.bit_table(&relations).map_err(DependencyContextError::Characters)?;
        let original_columns = original_rows(&relations, &layout, &bits)?;
        check_rows(&matrix, &original_columns)?;
        check_column_weights(&matrix)?;
        Ok(Self { relations, layout, characters, matrix, original_columns })
    }

    /// The checked polynomial and map context shared by every relation.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        self.relations.context()
    }

    /// The original relations; provenance indexes `relations().relations()`.
    #[must_use]
    pub fn relations(&self) -> &CheckedRelations {
        &self.relations
    }

    /// The column layout the matrix is read through.
    #[must_use]
    pub fn layout(&self) -> &MatrixLayout {
        &self.layout
    }

    /// The root characters whose bits fill the character columns.
    #[must_use]
    pub fn characters(&self) -> &RootCharacters {
        &self.characters
    }

    /// The filtered matrix whose row indices candidates select.
    #[must_use]
    pub fn matrix(&self) -> &SparseMatrix {
        &self.matrix
    }

    /// The ascending columns of original relation `relation`, character columns included, or
    /// `None` when there is no such relation.
    #[must_use]
    pub fn original_columns(&self, relation: usize) -> Option<&[usize]> {
        self.original_columns.get(relation).map(Vec::as_slice)
    }

    /// Validate a selection of filtered-matrix rows and compute its exact products.
    ///
    /// The selection is untrusted: it is re-checked, never passed through
    /// [`KernelVector::verify`] or [`KernelVector::expand_provenance`] (the latter panics on an
    /// out-of-range row). Kernel membership is checked in this context's matrix. The lift is the
    /// symmetric difference of the selected rows' provenance, so a relation named by two selected
    /// rows cancels and nothing is a union. The lifted parity is then recomputed from the original
    /// rows. Returned relations are ascending and non-empty, `P` and `γ` are exact (see the module
    /// documentation), and neither is required to be a square.
    ///
    /// # Errors
    ///
    /// [`DependencyError::Rejected`] with the first failing [`DependencyRejection`]: empty,
    /// unsorted, or out-of-range selection, then a nonzero sum over the matrix, then an empty
    /// lift. [`DependencyError::Invariant`] if the lifted relations have odd parity in some column
    /// after the matrix sum was zero, or if a lifted relation set is not a valid factor set. The
    /// context's construction makes both impossible, so they report a defect and are never
    /// rejections.
    pub fn validate(&self, selection: &[usize]) -> Result<CheckedDependency, DependencyError> {
        self.check_selection(selection)?;
        if let Some(&column) = self.selected_sum(selection)?.first() {
            return Err(DependencyError::Rejected(DependencyRejection::NotInKernel { column }));
        }
        let lifted = self.lift(selection)?;
        if lifted.is_empty() {
            return Err(DependencyError::Rejected(DependencyRejection::EmptyLift));
        }
        let parity = xor_columns(&self.original_columns, &lifted)
            .map_err(|_| invariant("lifted relation is not an original relation"))?;
        if !parity.is_empty() {
            return Err(invariant("lifted original rows do not sum to zero"));
        }
        let pairs = self.pairs_of(&lifted)?;
        let factors = checked_factors(self.context(), &pairs)?;
        Ok(CheckedDependency {
            selection: selection.to_vec(),
            relations: lifted,
            rational_product: rational_product(self.context(), &pairs),
            gamma: algebraic_product(self.context(), &pairs),
            factors,
        })
    }

    /// Validate a kernel vector's row indices, exactly as [`DependencyContext::validate`].
    ///
    /// # Errors
    ///
    /// As [`DependencyContext::validate`]. A vector built by struct literal, with unsorted or
    /// out-of-range indices, is rejected and never panics.
    pub fn validate_kernel_vector(
        &self,
        vector: &KernelVector,
    ) -> Result<CheckedDependency, DependencyError> {
        self.validate(&vector.row_indices)
    }

    fn check_selection(&self, selection: &[usize]) -> Result<(), DependencyError> {
        let reject = |rejection| Err(DependencyError::Rejected(rejection));
        if selection.is_empty() {
            return reject(DependencyRejection::EmptySelection);
        }
        if let Some(position) = selection.windows(2).position(|w| w[0] >= w[1]) {
            return reject(DependencyRejection::MalformedSelection { position: position + 1 });
        }
        let rows = self.matrix.rows.len();
        match selection.iter().position(|&row| row >= rows) {
            Some(position) => {
                let row = selection[position];
                reject(DependencyRejection::RowOutOfRange { position, row, rows })
            }
            None => Ok(()),
        }
    }

    /// The XOR of the selected rows' columns in the owned matrix; empty exactly for a kernel
    /// selection. Ascending, so its first entry is the smallest column of odd weight.
    fn selected_sum(&self, selection: &[usize]) -> Result<Vec<usize>, DependencyError> {
        let mut sum = Vec::new();
        for &row in selection {
            let row = self
                .matrix
                .rows
                .get(row)
                .ok_or_else(|| invariant("selected row is not a row of the matrix"))?;
            sum = symmetric_difference(&sum, &row.cols);
        }
        Ok(sum)
    }

    /// The symmetric difference of the selected rows' provenance, ascending.
    fn lift(&self, selection: &[usize]) -> Result<Vec<usize>, DependencyError> {
        self.matrix
            .lift_dependency(selection)
            .ok_or_else(|| invariant("selected row is not a row of the matrix"))
    }

    fn pairs_of(&self, lifted: &[usize]) -> Result<Vec<(BigInt, BigInt)>, DependencyError> {
        lifted
            .iter()
            .map(|&i| {
                self.relations
                    .relations()
                    .get(i)
                    .map(|r| (r.a().clone(), r.b().clone()))
                    .ok_or_else(|| invariant("lifted relation is not an original relation"))
            })
            .collect()
    }
}

/// Checker equality, character count, and matrix shape.
fn check_binding(
    relations: &CheckedRelations,
    matrix: &SparseMatrix,
    layout: &MatrixLayout,
    characters: &RootCharacters,
) -> Result<(), DependencyContextError> {
    if layout.checker() != relations.checker() {
        return Err(DependencyContextError::LayoutRelationsMismatch);
    }
    if characters.checker() != relations.checker() {
        return Err(DependencyContextError::CharactersRelationsMismatch);
    }
    if layout.character_count() != characters.len() {
        return Err(DependencyContextError::CharacterCountMismatch {
            layout: layout.character_count(),
            characters: characters.len(),
        });
    }
    if !layout.shape_matches(matrix) {
        return Err(DependencyContextError::MatrixShapeMismatch);
    }
    Ok(())
}

/// The ascending columns of every original relation, character columns included.
///
/// `bits[i][slot]` is the character bit of relation `i`. Character columns lie beyond the sign
/// column, so appending them in slot order keeps each list ascending.
fn original_rows(
    relations: &CheckedRelations,
    layout: &MatrixLayout,
    bits: &[Vec<bool>],
) -> Result<Vec<Vec<usize>>, DependencyContextError> {
    relations
        .relations()
        .iter()
        .zip(bits)
        .map(|(relation, relation_bits)| {
            // `None` means an index beyond the layout's sizes, which is what a foreign checker
            // looks like; the checker equality above makes it unreachable.
            let mut columns = layout
                .row_columns(relation)
                .ok_or(DependencyContextError::LayoutRelationsMismatch)?;
            columns.extend(
                relation_bits
                    .iter()
                    .enumerate()
                    .filter(|&(_, &bit)| bit)
                    .filter_map(|(slot, _)| layout.character_column(slot)),
            );
            Ok(columns)
        })
        .collect()
}

/// Per-row structure and the lifted-parity identity, rows in ascending order.
fn check_rows(
    matrix: &SparseMatrix,
    original_columns: &[Vec<usize>],
) -> Result<(), DependencyContextError> {
    for (row, entry) in matrix.rows.iter().enumerate() {
        if entry.cols.windows(2).any(|w| w[0] >= w[1]) {
            return Err(DependencyContextError::MalformedColumns { row });
        }
        if let Some(&column) = entry.cols.iter().find(|&&c| c >= matrix.num_cols) {
            return Err(DependencyContextError::ColumnOutOfRange { row, column });
        }
        if entry.provenance.is_empty() {
            return Err(DependencyContextError::EmptyProvenance { row });
        }
        if entry.provenance.windows(2).any(|w| w[0] >= w[1]) {
            return Err(DependencyContextError::MalformedProvenance { row });
        }
        let expected = xor_columns(original_columns, &entry.provenance).map_err(|index| {
            DependencyContextError::ProvenanceOutOfRange {
                row,
                index,
                relations: original_columns.len(),
            }
        })?;
        if let Some(&column) = symmetric_difference(&expected, &entry.cols).first() {
            return Err(DependencyContextError::LiftedParityMismatch { row, column });
        }
    }
    Ok(())
}

/// The recorded weight of every column equals the number of rows that set it.
fn check_column_weights(matrix: &SparseMatrix) -> Result<(), DependencyContextError> {
    let mut weights = vec![0u64; matrix.num_cols];
    for row in &matrix.rows {
        for &column in &row.cols {
            if let Some(weight) = weights.get_mut(column) {
                *weight += 1;
            }
        }
    }
    match matrix.col_weights.iter().zip(&weights).position(|(&w, &t)| u64::from(w) != t) {
        Some(column) => Err(DependencyContextError::ColumnWeightMismatch { column }),
        None => Ok(()),
    }
}

/// `RelationFactors` for a lifted relation set, which the checked relations always allow.
fn checked_factors(
    ctx: &SqrtContext,
    pairs: &[(BigInt, BigInt)],
) -> Result<RelationFactors, DependencyError> {
    RelationFactors::new(ctx, pairs)
        .map_err(|_| invariant("lifted relations are not a valid factor set"))
}

/// `P = Π(a − b·m)` over `Z` with the exact `m`.
fn rational_product(ctx: &SqrtContext, pairs: &[(BigInt, BigInt)]) -> BigInt {
    pairs.iter().map(|(a, b)| a - b * ctx.root()).product()
}

/// `γ = Π(c·a − b·t) mod F / c^k` in the `θ` basis, reduced after every factor over `Z`.
fn algebraic_product(ctx: &SqrtContext, pairs: &[(BigInt, BigInt)]) -> RatPoly {
    let scaled = pairs.iter().fold(IntPoly::one(), |acc, (a, b)| {
        let factor = IntPoly::from_coeffs(vec![ctx.leading_coeff() * a, -b]);
        acc.mul(&factor).pseudo_div_rem(ctx.monic_poly()).1
    });
    let denominator = pow(ctx.leading_coeff().clone(), pairs.len());
    scaled.to_rat_poly().scale(&BigRational::new(BigInt::one(), denominator))
}

// ─── the checked dependency ──────────────────────────────────────────────────

/// A validated dependency: selected rows, the original relations they lift to, and exact products.
///
/// Minted only by [`DependencyContext::validate`]. It is a *candidate* for a congruence of squares,
/// not a certificate: `P` is an exact integer square by the parity checks, while `γ` may still be
/// a nonsquare, which exact reconstruction decides.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckedDependency {
    selection: Vec<usize>,
    relations: Vec<usize>,
    factors: RelationFactors,
    rational_product: BigInt,
    gamma: RatPoly,
}

impl CheckedDependency {
    /// The checked polynomial and map context.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        self.factors.context()
    }

    /// The filtered-matrix rows, as validated.
    #[must_use]
    pub fn selection(&self) -> &[usize] {
        &self.selection
    }

    /// The lifted original relation indices: ascending and non-empty.
    #[must_use]
    pub fn relations(&self) -> &[usize] {
        &self.relations
    }

    /// The `(a, b)` pairs in [`CheckedDependency::relations`] order, owning the context.
    #[must_use]
    pub fn factors(&self) -> &RelationFactors {
        &self.factors
    }

    /// The number of lifted relations, `k >= 1`.
    #[must_use]
    pub fn relation_count(&self) -> usize {
        self.relations.len()
    }

    /// `P = Π(a − b·m)` over `Z` with the exact `m`.
    #[must_use]
    pub fn rational_product(&self) -> &BigInt {
        &self.rational_product
    }

    /// `γ = Π(a − b·α)` as `Π(c·a − b·θ) / c^k` modulo `F`, in the `θ` basis.
    #[must_use]
    pub fn gamma(&self) -> &RatPoly {
        &self.gamma
    }

    /// The normalized target of this dependency's factors under `limits`.
    ///
    /// # Errors
    ///
    /// [`SqrtError`] as for [`NormalizedTarget::new`].
    pub fn target(&self, limits: PrecisionLimits) -> Result<TargetOutcome, SqrtError> {
        NormalizedTarget::new(self.factors.clone(), limits)
    }
}

// ─── sorted-set helpers ──────────────────────────────────────────────────────

/// Symmetric difference of two strictly ascending lists.
fn symmetric_difference(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// The XOR of `columns[i]` over `indices`; `Err(i)` for the first index with no entry.
fn xor_columns(columns: &[Vec<usize>], indices: &[usize]) -> Result<Vec<usize>, usize> {
    let mut sum = Vec::new();
    for &index in indices {
        let original = columns.get(index).ok_or(index)?;
        sum = symmetric_difference(&sum, original);
    }
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::{MatrixRow, build_checked_matrix};
    use crate::sieve::{ExponentVector, FactorBase, Relation, RelationChecker};

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn poly(coeffs: &[i64]) -> IntPoly {
        IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
    }

    fn relation(a: i64, b: i64, rational: &[(usize, u32)], algebraic: &[(usize, u32)]) -> Relation {
        Relation {
            a: big(a),
            b: big(b),
            rational_exponents: ExponentVector { entries: rational.to_vec() },
            algebraic_exponents: ExponentVector { entries: algebraic.to_vec() },
            rational_sign: a - 3 * b < 0,
        }
    }

    /// f = 3x² + x − 1, m = 3, N = 29, FB(7, 17); relations (4,1), (7,1), (2,1) of the
    /// checked-relation fixture, with columns `[6]`, `[7]`, `[8]` when K = 0.
    fn fixture() -> (CheckedRelations, SparseMatrix, MatrixLayout, RootCharacters) {
        let ctx = SqrtContext::new(&poly(&[-1, 1, 3]), &big(3), &big(29)).unwrap();
        let fb = FactorBase::new(&poly(&[-1, 1, 3]), 7, 17);
        let checker = RelationChecker::new(&ctx, &fb).unwrap();
        let relations = CheckedRelations::new(
            &checker,
            &[
                relation(4, 1, &[], &[(0, 1), (2, 1)]),
                relation(7, 1, &[(0, 2)], &[(0, 2), (3, 1)]),
                relation(2, 1, &[], &[(1, 1)]),
            ],
        )
        .unwrap();
        let (matrix, layout) = build_checked_matrix(&relations, 0).unwrap();
        let characters = RootCharacters::new(&checker, &[]).unwrap();
        (relations, matrix, layout, characters)
    }

    fn context_from(matrix_edit: impl FnOnce(&mut SparseMatrix)) -> DependencyContext {
        let (relations, mut matrix, layout, characters) = fixture();
        matrix_edit(&mut matrix);
        DependencyContext::new(relations, matrix, layout, characters).unwrap()
    }

    /// A hand-built matrix whose rows are XORs of original rows with overlapping provenance:
    /// `[0,1]`, `[1,2]`, `[0,2]` over the original rows `[6]`, `[7]`, `[8]`.
    fn cyclic_context() -> DependencyContext {
        context_from(|m| {
            m.rows = vec![
                MatrixRow { cols: vec![6, 7], provenance: vec![0, 1] },
                MatrixRow { cols: vec![7, 8], provenance: vec![1, 2] },
                MatrixRow { cols: vec![6, 8], provenance: vec![0, 2] },
            ];
            m.col_weights = vec![0, 0, 0, 0, 0, 0, 2, 2, 2];
        })
    }

    #[test]
    fn symmetric_difference_cancels_shared_entries() {
        assert_eq!(symmetric_difference(&[0, 2, 5], &[2, 3, 5, 9]), vec![0, 3, 9]);
        assert_eq!(symmetric_difference(&[], &[4]), vec![4]);
        assert!(symmetric_difference(&[1, 2], &[1, 2]).is_empty());
    }

    #[test]
    fn xor_columns_reports_the_first_missing_index() {
        let columns = vec![vec![1, 2], vec![2, 3]];
        assert_eq!(xor_columns(&columns, &[0, 1]), Ok(vec![1, 3]));
        assert_eq!(xor_columns(&columns, &[0, 7, 9]), Err(7));
    }

    #[test]
    fn hand_built_overlapping_provenance_cancels_to_an_empty_lift() {
        // Every relation index appears in exactly two of the three rows, so the lift is empty
        // even though the three rows sum to zero: the pairs would be a cycle, not a dependency.
        let context = cyclic_context();
        assert_eq!(
            context.validate(&[0, 1, 2]),
            Err(DependencyError::Rejected(DependencyRejection::EmptyLift))
        );
        // Rows 0 and 1 sum to [6, 7] ⊕ [7, 8] = [6, 8], whose smallest odd column is 6.
        assert_eq!(
            context.validate(&[0, 1]),
            Err(DependencyError::Rejected(DependencyRejection::NotInKernel { column: 6 }))
        );
    }

    #[test]
    fn corrupted_original_columns_are_an_invariant_error() {
        // A genuine matrix plus a copy of row 2 (columns `[8]`, provenance `[2]`).
        let (relations, mut matrix, layout, characters) = fixture();
        matrix.rows.push(MatrixRow { cols: vec![8], provenance: vec![2] });
        matrix.col_weights[8] += 1;
        let mut context = DependencyContext::new(relations, matrix, layout, characters).unwrap();
        // Rows 2 and 3 are both `[8]` with provenance `[2]`: the lift of `{2, 3}` is empty.
        assert_eq!(
            context.validate(&[2, 3]),
            Err(DependencyError::Rejected(DependencyRejection::EmptyLift))
        );
        // Point row 3 at relation 0 while its columns stay `[8]`: `{2, 3}` is still a kernel
        // selection, but its lift is `[0, 2]` and the original rows `[6] ⊕ [8]` do not cancel.
        context.matrix.rows[3].provenance = vec![0];
        assert_eq!(
            context.validate(&[2, 3]),
            Err(invariant("lifted original rows do not sum to zero"))
        );
    }

    #[test]
    fn a_lifted_relation_without_an_original_row_is_an_invariant_error() {
        let (relations, mut matrix, layout, characters) = fixture();
        matrix.rows.push(MatrixRow { cols: vec![8], provenance: vec![2] });
        matrix.col_weights[8] += 1;
        let mut context = DependencyContext::new(relations, matrix, layout, characters).unwrap();
        context.matrix.rows[3].provenance = vec![99];
        assert_eq!(
            context.validate(&[2, 3]),
            Err(invariant("lifted relation is not an original relation"))
        );
    }

    #[test]
    fn row_lookups_past_the_matrix_are_invariant_errors() {
        let context = context_from(|_| {});
        let missing = Err(invariant("selected row is not a row of the matrix"));
        assert_eq!(context.selected_sum(&[99]), missing);
        assert_eq!(context.lift(&[99]), missing);
    }

    #[test]
    fn lifted_pairs_and_factor_sets_are_checked_not_assumed() {
        let context = context_from(|_| {});
        assert_eq!(
            context.pairs_of(&[9]),
            Err(invariant("lifted relation is not an original relation"))
        );
        // (3, 1) has a − b·m = 0, which `RelationFactors` rejects.
        assert_eq!(
            checked_factors(context.context(), &[(big(3), big(1))]),
            Err(invariant("lifted relations are not a valid factor set"))
        );
    }

    #[test]
    fn errors_display_classify_and_chain() {
        let context_errors = [
            DependencyContextError::LayoutRelationsMismatch,
            DependencyContextError::CharactersRelationsMismatch,
            DependencyContextError::CharacterCountMismatch { layout: 2, characters: 1 },
            DependencyContextError::MatrixShapeMismatch,
            DependencyContextError::Characters(RootCharacterError::MatrixShapeMismatch),
            DependencyContextError::MalformedColumns { row: 3 },
            DependencyContextError::ColumnOutOfRange { row: 3, column: 40 },
            DependencyContextError::EmptyProvenance { row: 3 },
            DependencyContextError::MalformedProvenance { row: 3 },
            DependencyContextError::ProvenanceOutOfRange { row: 3, index: 8, relations: 6 },
            DependencyContextError::LiftedParityMismatch { row: 3, column: 5 },
            DependencyContextError::ColumnWeightMismatch { column: 5 },
        ];
        for error in &context_errors {
            let text = error.to_string();
            assert!(!text.is_empty(), "{error:?}");
            let chained = std::error::Error::source(error).is_some();
            assert_eq!(chained, matches!(error, DependencyContextError::Characters(_)));
        }
        assert!(context_errors[5].to_string().contains("row 3"));
        assert!(context_errors[6].to_string().contains("column 40"));

        let rejections = [
            DependencyRejection::EmptySelection,
            DependencyRejection::MalformedSelection { position: 1 },
            DependencyRejection::RowOutOfRange { position: 2, row: 9, rows: 6 },
            DependencyRejection::NotInKernel { column: 6 },
            DependencyRejection::EmptyLift,
        ];
        for rejection in &rejections {
            assert!(!rejection.to_string().is_empty(), "{rejection:?}");
            let error = DependencyError::Rejected(rejection.clone());
            assert_eq!(error.class(), FailureClass::InvalidInput);
            assert!(error.to_string().contains(&rejection.to_string()));
            assert!(std::error::Error::source(&error).is_some());
        }
        let broken = invariant("an identity");
        assert_eq!(broken.class(), FailureClass::InvariantViolation);
        assert!(broken.to_string().contains("an identity"));
        assert!(std::error::Error::source(&broken).is_some());
    }

    #[test]
    fn accessors_expose_the_bound_parts() {
        let context = context_from(|_| {});
        assert_eq!(context.relations().len(), 3);
        assert_eq!(context.layout().num_cols(), 9);
        assert!(context.characters().is_empty());
        assert_eq!(context.matrix().rows.len(), 3);
        assert_eq!(context.original_columns(1), Some(&[7][..]));
        assert_eq!(context.original_columns(3), None);
        assert_eq!(context.context().modulus(), &big(29));
    }
}
