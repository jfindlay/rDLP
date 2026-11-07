//! Filter substrate for GNFS: sparse GF(2) matrix construction, singleton removal,
//! clique/excess pruning, and column merging.
//!
//! This module is the entry point for the ``gnfs::filter`` sub-crate. It provides:
//!
//! - [`matrix`] — the sparse GF(2) matrix type (``SparseMatrix``, ``MatrixRow``) and
//!   the ``EXCESS_FLOOR`` constant.
//! - [`singleton`] — singleton removal to fixpoint (``remove_singletons``).
//! - [`merge`] — clique/excess pruning (``prune_cliques``) and column merging
//!   (``merge_columns``) over the singleton-removed matrix (clique pruning and merging step).
//! - [`build_matrix`] — constructs the initial ``SparseMatrix`` from a relation corpus.
//! - [`build_checked_matrix`] and [`MatrixLayout`] — the same layout built from validated
//!   relations, with explicit classification of every column.
//!
//! # Background
//!
//! The filtering step takes the relation corpus from sieving and reduces it to a
//! well-overdetermined sparse GF(2) matrix suitable for linear algebra. The two main
//! operations are:
//!
//! 1. **Singleton removal**: primes/ideals appearing in only one surviving relation cannot
//!    contribute to a GF(2) dependency; remove their rows and iterate to a fixpoint.
//!
//! 2. **Clique pruning and merging**: reduce matrix weight while preserving excess >=
//!    ``EXCESS_FLOOR`` (pruning), then eliminate low-weight columns by XOR-merging the
//!    rows that contain them (merging).
//!
//! # Column layout
//!
//! Total columns = ``FactorBase::matrix_width()`` = rational_size + algebraic_size + obstruction_count.
//!
//! - ``[0, rational_size)``: rational factor-base columns (GF(2) parities of rational exponents).
//! - ``[rational_size, rational_size + algebraic_size)``: algebraic columns (GF(2) parities of
//!   algebraic exponents).
//! - ``[rational_size + algebraic_size, matrix_width)``: obstruction columns. The sign bit
//!   (``rational_sign``) occupies the first obstruction column; quadratic-character columns
//!   follow and are filled by the linear algebra step (carried as zeros here).
//!
//! [`build_matrix`] treats every factor-base ideal column as a parity constraint.
//! [`build_checked_matrix`] uses the identical column order but leaves the columns of
//! exceptional ideals permanently zero; see [`MatrixLayout`].
//!
//! # Filtering contract
//!
//! The types and functions in this module implement the sparse GF(2) matrix substrate.
//! The clique pruning and merging step, the linear algebra step, and the square root step
//! consume this interface directly.

pub mod matrix;
pub mod merge;
pub mod singleton;

pub use matrix::{MatrixRow, SparseMatrix, EXCESS_FLOOR};
pub use merge::{merge_columns, prune_cliques};
pub use singleton::remove_singletons;

use std::fmt;

use crate::sieve::{
    CheckedRelation, CheckedRelations, FactorBase, PrimeClass, Relation, RelationChecker,
};

// ─── build_matrix ─────────────────────────────────────────────────────────────

/// Build the initial sparse GF(2) matrix from a relation corpus and factor base.
///
/// Each relation in ``relations`` contributes one row. Column layout:
///
/// - Rational columns ``[0, fb.rational_size())``: GF(2) parities of rational exponents.
/// - Algebraic columns ``[fb.rational_size(), fb.rational_size() + fb.algebraic_size())``:
///   GF(2) parities of algebraic exponents.
/// - Obstruction columns ``[fb.rational_size() + fb.algebraic_size(), fb.matrix_width())``:
///   sign bit at ``fb.rational_size() + fb.algebraic_size()`` (from ``relation.rational_sign``);
///   remaining obstruction columns (quadratic characters) set to 0 — the linear algebra step
///   fills them.
///
/// Note: ``Relation::rational_row_gf2`` places the sign at local index 0 of its return
/// value; ``build_matrix`` re-maps it to the global obstruction column index.
///
/// Provenance for row ``i`` is ``[i]`` (the original relation index).
///
/// This builder trusts its input and treats every factor-base ideal column as a parity
/// constraint, including ideals above primes dividing ``c·disc(F)`` where the recorded exponent
/// is not an ideal valuation. [`build_checked_matrix`] validates the relations and omits those
/// constraints.
///
/// :param relations: The relation corpus from sieving.
/// :param fb: The factor base (for size and column layout).
/// :returns: The initial sparse GF(2) matrix.
pub fn build_matrix(relations: &[Relation], fb: &FactorBase) -> SparseMatrix {
    let num_cols = fb.matrix_width();
    let rat_size = fb.rational_size();
    let alg_size = fb.algebraic_size();
    let obstruction_col_start = rat_size + alg_size;
    let obstruction_count = fb.obstruction_count;

    let mut col_weights = vec![0u32; num_cols];
    let mut rows = Vec::with_capacity(relations.len());

    for (i, relation) in relations.iter().enumerate() {
        // rational_row_gf2 returns [sign_bit, rat_col_0, rat_col_1, ...].
        // Local index 0 = sign bit → global column obstruction_col_start.
        // Local index 1+k = rational column k → global column k.
        let rat_row = relation.rational_row_gf2(fb);

        // algebraic_row_gf2 returns [alg_col_0, ..., alg_col_{alg_size-1}, obstruction_zeros...].
        // Local index k (k < alg_size) → global column rat_size + k.
        // Trailing obstruction zeros are already zero; we don't need to read them.
        let alg_row = relation.algebraic_row_gf2(fb);

        let mut cols: Vec<usize> = Vec::new();

        // Rational columns: local indices 1..=rat_size → global 0..rat_size.
        for k in 0..rat_size {
            if rat_row[1 + k] {
                cols.push(k);
            }
        }

        // Algebraic columns: local indices 0..alg_size → global rat_size..rat_size+alg_size.
        for k in 0..alg_size {
            if alg_row[k] {
                cols.push(rat_size + k);
            }
        }

        // Sign bit (obstruction column 0): local rat_row[0] → global obstruction_col_start.
        if rat_row[0] {
            cols.push(obstruction_col_start);
        }

        // Remaining obstruction columns (quadratic chars): carried as zeros, nothing to push.
        // (obstruction_col_start + 1 .. num_cols are all zero for the filtering step.)

        // cols is already sorted because we pushed in ascending global-column order:
        // rational (0..rat_size) < algebraic (rat_size..obstruction_col_start) < obstruction.
        debug_assert!(cols.windows(2).all(|w| w[0] < w[1]), "cols must be sorted and unique");

        // Update column weights.
        for &c in &cols {
            col_weights[c] += 1;
        }

        rows.push(MatrixRow { cols, provenance: vec![i] });
    }

    SparseMatrix { rows, num_cols, obstruction_col_start, obstruction_count, col_weights }
}

// ─── checked matrix ──────────────────────────────────────────────────────────

/// The meaning of one column of a checked matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColumnClass {
    /// Parity of the exponent of a rational prime; columns ``[0, R)``.
    RationalPrime {
        /// The prime.
        prime: u64,
    },
    /// Parity of the exponent of an ideal above a good prime; columns ``[R, R + A)``.
    GoodIdeal {
        /// The prime ``p``.
        prime: u64,
        /// The root ``r`` of the ideal ``(p, r)``.
        root: u64,
    },
    /// An ideal above a prime dividing ``c·disc(F)``; same range, omitted: never set.
    ExceptionalIdeal {
        /// The prime ``p``.
        prime: u64,
        /// The root ``r`` of the ideal ``(p, r)``.
        root: u64,
    },
    /// The parity of the number of negative rational norms; column ``R + A``.
    RationalSign,
    /// A quadratic-character slot; column ``R + A + 1 + slot``. Zero when the matrix is built.
    RootCharacter {
        /// The slot, below the layout's character count.
        slot: usize,
    },
}

/// The column count ``R + A + 1 + K`` does not fit in ``usize``.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColumnCountOverflow;

impl fmt::Display for ColumnCountOverflow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("matrix column count overflows usize")
    }
}

impl std::error::Error for ColumnCountOverflow {}

/// The column layout of a checked matrix, bound to the checker that defines it.
///
/// With ``R`` rational primes, ``A`` ideals, and ``K`` character slots, the columns are
/// ``[0, R)`` rational primes, ``[R, R + A)`` ideals in factor-base index order, ``R + A`` the
/// rational sign, and ``R + A + 1 + slot`` the root characters. This is the order of
/// [`build_matrix`], so singleton removal, merging, and the sign-then-characters convention read
/// it unchanged. Columns of exceptional ideals are kept, and are zero in every row under any
/// XOR, so column indices coincide with factor-base indices. ``SparseMatrix::excess`` therefore
/// undercounts by ``omitted_columns().len()``, which only makes pruning more conservative.
///
/// The layout owns a clone of its checker; a matrix, layout, and relation set belong together
/// exactly when their checkers compare equal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixLayout {
    checker: RelationChecker,
    character_count: usize,
    num_cols: usize,
}

impl MatrixLayout {
    /// Lay out ``checker``'s factor base with ``character_count`` character slots.
    fn new(
        checker: &RelationChecker,
        character_count: usize,
    ) -> Result<Self, ColumnCountOverflow> {
        let num_cols = checker
            .rational_primes()
            .len()
            .checked_add(checker.ideals().len())
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_add(character_count))
            .ok_or(ColumnCountOverflow)?;
        Ok(Self { checker: checker.clone(), character_count, num_cols })
    }

    /// The bound checker.
    #[must_use]
    pub fn checker(&self) -> &RelationChecker {
        &self.checker
    }

    /// The number of rational-prime columns ``R``.
    #[must_use]
    pub fn rational_size(&self) -> usize {
        self.checker.rational_primes().len()
    }

    /// The number of ideal columns ``A``, including omitted ones.
    #[must_use]
    pub fn algebraic_size(&self) -> usize {
        self.checker.ideals().len()
    }

    /// The number of character slots ``K``.
    #[must_use]
    pub fn character_count(&self) -> usize {
        self.character_count
    }

    /// The rational-sign column ``R + A``, which is also the first obstruction column.
    #[must_use]
    pub fn sign_column(&self) -> usize {
        self.rational_size() + self.algebraic_size()
    }

    /// The first obstruction column ``R + A``.
    #[must_use]
    pub fn obstruction_col_start(&self) -> usize {
        self.sign_column()
    }

    /// The number of obstruction columns ``1 + K``.
    #[must_use]
    pub fn obstruction_count(&self) -> usize {
        1 + self.character_count
    }

    /// The total number of columns ``R + A + 1 + K``.
    #[must_use]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// The column of character ``slot``, or ``None`` when ``slot >= K``.
    #[must_use]
    pub fn character_column(&self, slot: usize) -> Option<usize> {
        (slot < self.character_count).then(|| self.sign_column() + 1 + slot)
    }

    /// The meaning of ``column``, or ``None`` when ``column >= num_cols``.
    #[must_use]
    pub fn class(&self, column: usize) -> Option<ColumnClass> {
        let rational = self.rational_size();
        let sign = self.sign_column();
        if column < rational {
            Some(ColumnClass::RationalPrime { prime: self.checker.rational_primes()[column] })
        } else if column < sign {
            let index = column - rational;
            let (prime, root) = self.checker.ideals()[index];
            match self.checker.ideal_class(index)? {
                PrimeClass::Good => Some(ColumnClass::GoodIdeal { prime, root }),
                PrimeClass::Exceptional => Some(ColumnClass::ExceptionalIdeal { prime, root }),
            }
        } else if column == sign {
            Some(ColumnClass::RationalSign)
        } else if column < self.num_cols {
            Some(ColumnClass::RootCharacter { slot: column - sign - 1 })
        } else {
            None
        }
    }

    /// The columns of exceptional ideals, ascending; they are zero in every row.
    #[must_use]
    pub fn omitted_columns(&self) -> Vec<usize> {
        let rational = self.rational_size();
        (0..self.algebraic_size())
            .filter(|&i| self.checker.ideal_class(i) == Some(PrimeClass::Exceptional))
            .map(|i| rational + i)
            .collect()
    }

    /// The ascending columns set by ``relation``, character slots not included.
    ///
    /// These are the rational indices with odd exponent, ``R + index`` for every *good* ideal with
    /// odd exponent, and the sign column when the rational norm is negative. Exceptional entries
    /// contribute nothing, and there is no algebraic-sign column. Returns ``None`` only when an
    /// index exceeds this layout's sizes, i.e. the relation comes from a different checker.
    #[must_use]
    pub fn row_columns(&self, relation: &CheckedRelation) -> Option<Vec<usize>> {
        let in_range = relation.rational_exponents().iter().all(|&(i, _)| i < self.rational_size())
            && relation.algebraic_exponents().iter().all(|&(i, _)| i < self.algebraic_size());
        in_range.then(|| self.odd_columns(relation))
    }

    /// True iff ``matrix`` has this layout's width, obstruction block, and weight table.
    ///
    /// Rows and provenance are not examined.
    #[must_use]
    pub fn shape_matches(&self, matrix: &SparseMatrix) -> bool {
        matrix.num_cols == self.num_cols
            && matrix.obstruction_col_start == self.obstruction_col_start()
            && matrix.obstruction_count == self.obstruction_count()
            && matrix.col_weights.len() == self.num_cols
    }

    /// Row columns of a relation whose indices are already known to be in range.
    fn odd_columns(&self, relation: &CheckedRelation) -> Vec<usize> {
        let rational = self.rational_size();
        let mut cols: Vec<usize> = relation
            .rational_exponents()
            .iter()
            .filter(|&&(_, e)| e % 2 == 1)
            .map(|&(i, _)| i)
            .collect();
        cols.extend(
            relation
                .algebraic_exponents()
                .iter()
                .filter(|&&(i, e)| {
                    e % 2 == 1 && self.checker.ideal_class(i) == Some(PrimeClass::Good)
                })
                .map(|&(i, _)| rational + i),
        );
        if relation.rational_sign() {
            cols.push(self.sign_column());
        }
        cols
    }
}

/// Build the matrix of validated relations, with the layout that interprets it.
///
/// Row ``i`` comes from ``relations.relations()[i]`` and has provenance ``[i]``; its columns are
/// those of [`MatrixLayout::row_columns`]. Character columns are zero, and an empty row is kept
/// because it is itself a dependency. ``col_weights`` is exact. The matrix has
/// ``R + A + 1 + character_count`` columns, with the obstruction block starting at ``R + A``.
///
/// # Errors
///
/// [`ColumnCountOverflow`] when ``R + A + 1 + character_count`` overflows ``usize``.
pub fn build_checked_matrix(
    relations: &CheckedRelations,
    character_count: usize,
) -> Result<(SparseMatrix, MatrixLayout), ColumnCountOverflow> {
    let layout = MatrixLayout::new(relations.checker(), character_count)?;
    let mut col_weights = vec![0u32; layout.num_cols];
    let mut rows = Vec::with_capacity(relations.len());
    for (i, relation) in relations.relations().iter().enumerate() {
        // The relations were checked by this very checker, so every index is in range.
        let cols = layout.odd_columns(relation);
        for &c in &cols {
            if let Some(weight) = col_weights.get_mut(c) {
                *weight += 1;
            }
        }
        rows.push(MatrixRow { cols, provenance: vec![i] });
    }
    let matrix = SparseMatrix {
        rows,
        num_cols: layout.num_cols,
        obstruction_col_start: layout.obstruction_col_start(),
        obstruction_count: layout.obstruction_count(),
        col_weights,
    };
    Ok((matrix, layout))
}
