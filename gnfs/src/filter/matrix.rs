//! Sparse GF(2) matrix with relation-provenance map for GNFS filtering.
//!
//! This module defines the core data structures consumed by singleton removal, clique pruning
//! and merging, and the linear algebra step. The representation is row-major (each row =
//! sorted ``Vec<usize>`` of set-column indices) with a maintained column-weight side table,
//! balancing the needs of row-XOR operations (merging) and column-weight queries (singleton
//! removal).
//!
//! # Column layout
//!
//! Total columns = rational_size + algebraic_size + obstruction_count, where
//! ``[0, rational_size)`` are rational factor-base columns, the next ``algebraic_size`` columns
//! are algebraic factor-base columns, and the remaining obstruction columns hold the sign bit
//! first, then quadratic-character columns.
//!
//! The width is ``FactorBase::matrix_width()`` only for the legacy ``build_matrix``. A matrix
//! from ``build_checked_matrix`` has ``MatrixLayout::num_cols()`` columns, whose obstruction
//! count is ``1 + K`` for ``K`` character slots rather than the factor base's own count; that
//! layout also keeps the columns of exceptional ideals as permanently zero columns. Always read
//! the width from ``SparseMatrix::num_cols``, never recompute it from a factor base.
//!
//! # Provenance
//!
//! Each row carries a sorted ``Vec<usize>`` of *original* relation indices (indices into the
//! relation list the matrix was built from). For a freshly built row, provenance =
//! ``[original_index]``. A row is always the GF(2) sum of the original rows its provenance
//! names, so XOR-combining two rows takes the *symmetric difference* of their provenance, exactly
//! as it does for their columns: an original index shared by both rows cancels. Overlap arises
//! once weight-three elimination has run, so a union would silently produce a wrong lift.
//! ``SparseMatrix::lift_dependency`` expands a set of dependent rows to the original relations.

// ─── EXCESS_FLOOR ─────────────────────────────────────────────────────────────

/// Minimum excess the clique pruning step must preserve.
///
/// excess = rows − (columns − obstruction_count); see [`SparseMatrix::excess`] for how omitted
/// columns of a checked matrix affect it. At toy scale any positive excess
/// suffices for a non-trivial nullspace; 20 is a conservative floor that keeps the
/// matrix well-overdetermined. Annotated as scale-dependent (principle-4): at
/// cryptographic scale the floor is typically set to ~200 or a fraction of the column
/// count. Singleton removal defines the constant; clique pruning enforces it.
pub const EXCESS_FLOOR: usize = 20;

// ─── MatrixRow ────────────────────────────────────────────────────────────────

/// One row of the sparse GF(2) matrix, with its relation-provenance record.
///
/// ``cols`` is a sorted ``Vec<usize>`` of column indices where this row has a 1 (GF(2)),
/// laid out as described in the module documentation: factor-base columns, then the sign
/// column, then quadratic-character columns (zero until filled in by the linear algebra step).
///
/// ``provenance`` is a sorted ``Vec<usize>`` of *original* relation indices with no repeats.
/// For a freshly built row, provenance = ``[original_index]``. The row's ``cols`` equal the
/// GF(2) sum of the original rows named by ``provenance``, and XOR-combining rows takes the
/// symmetric difference of both vectors. Provenance stores original indices, not pre-reduced
/// row sums, so the square root step can recover the actual ``(a, b)`` pairs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixRow {
    /// Sorted column indices where this row has a 1 in GF(2).
    pub cols: Vec<usize>,
    /// Sorted, deduplicated original relation indices this row derives from.
    pub provenance: Vec<usize>,
}

impl MatrixRow {
    /// XOR this row with ``other`` in GF(2).
    ///
    /// Both ``cols`` and ``provenance`` take the symmetric difference: a column set in both rows
    /// cancels, and an original relation index named by both rows cancels from the provenance.
    /// The result is therefore the GF(2) sum of the original rows its provenance names, even
    /// when the operands' provenance overlaps.
    ///
    /// :param other: The row to XOR with.
    /// :returns: A new row with the symmetric difference of columns and of provenance.
    pub fn xor_merge(&self, other: &MatrixRow) -> MatrixRow {
        MatrixRow {
            cols: symmetric_difference(&self.cols, &other.cols),
            provenance: symmetric_difference(&self.provenance, &other.provenance),
        }
    }
}

/// Symmetric difference of two strictly ascending index lists, in O(n + m).
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
                // Present in both: the GF(2) sum is zero here.
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

// ─── SparseMatrix ─────────────────────────────────────────────────────────────

/// Sparse GF(2) matrix over a factor base, with relation-provenance map.
///
/// The obstruction block starts at ``obstruction_col_start`` = rational_size + algebraic_size
/// and has ``obstruction_count`` columns. Singleton removal and elimination skip any column
/// >= ``obstruction_col_start``. See the module documentation for which constructor fixes
/// ``num_cols`` and ``obstruction_count``.
///
/// ``col_weights[c]`` = number of rows with a 1 in column ``c``. Maintained in sync with
/// ``rows`` by every mutating operation in this crate (``remove_row``, elimination).
///
/// ``excess()`` = rows.len() − (num_cols − obstruction_count). Singleton removal reports it;
/// clique pruning enforces >= ``EXCESS_FLOOR``.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    /// The rows of the matrix, each with its column set and provenance.
    pub rows: Vec<MatrixRow>,
    /// Total number of columns. ``FactorBase::matrix_width()`` for the legacy builder;
    /// ``MatrixLayout::num_cols()`` (``R + A + 1 + K``) for a checked matrix.
    pub num_cols: usize,
    /// Index of the first obstruction column = rational_size + algebraic_size.
    pub obstruction_col_start: usize,
    /// Number of obstruction columns: the factor base's ``obstruction_count`` for the legacy
    /// builder, ``1 + K`` for a checked matrix.
    pub obstruction_count: usize,
    /// col_weights[c] = number of rows with a 1 in column c.
    pub col_weights: Vec<u32>,
}

impl SparseMatrix {
    /// Current excess: rows.len() − (num_cols − obstruction_count).
    ///
    /// Positive excess means the matrix is overdetermined (more rows than non-obstruction
    /// columns), which is required for a non-trivial nullspace. Clique pruning enforces
    /// excess >= ``EXCESS_FLOOR`` during pruning.
    ///
    /// The count treats every non-obstruction column as active. A checked matrix keeps the
    /// columns of exceptional ideals as always-zero columns, so for it this value undercounts
    /// the true excess by the number of omitted columns, which only makes pruning conservative.
    pub fn excess(&self) -> isize {
        self.rows.len() as isize - (self.num_cols as isize - self.obstruction_count as isize)
    }

    /// Remove a row by index, updating ``col_weights``.
    ///
    /// Used by singleton removal and clique pruning.
    ///
    /// IMPORTANT: uses ordered removal (``Vec::remove``, not ``swap_remove``) to preserve
    /// row-index stability during singleton fixpoint iteration.
    ///
    /// :param row_idx: Index of the row to remove.
    /// :panics: If ``row_idx`` is out of bounds.
    pub fn remove_row(&mut self, row_idx: usize) {
        let row = self.rows.remove(row_idx);
        for &col in &row.cols {
            self.col_weights[col] -= 1;
        }
    }

    /// Expand a set of rows of this matrix to the original relations they stand for.
    ///
    /// The result is the GF(2) sum of the selected rows' provenance, as a sorted list of
    /// original relation indices. Repeated entries of ``selected`` cancel in pairs, since the
    /// selection is a GF(2) vector. When the selected rows sum to zero in this matrix, the
    /// returned relations sum to zero in the matrix the rows were derived from, provided every
    /// transformation between the two was a row XOR (elimination) or a row drop (pruning,
    /// singleton removal).
    ///
    /// :param selected: Row indices into ``self.rows``.
    /// :returns: Sorted original relation indices, or ``None`` if an index is out of range.
    pub fn lift_dependency(&self, selected: &[usize]) -> Option<Vec<usize>> {
        let mut lifted: Vec<usize> = Vec::new();
        for &row_idx in selected {
            lifted = symmetric_difference(&lifted, &self.rows.get(row_idx)?.provenance);
        }
        Some(lifted)
    }
}
