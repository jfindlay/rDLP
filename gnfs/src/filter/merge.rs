//! Clique/excess pruning and kernel-preserving column elimination for the sparse GF(2) matrix.
//!
//! This module implements the two post-singleton-removal reduction steps. They have different
//! guarantees and are reported separately:
//!
//! 1. **Clique/excess pruning** ([`prune_cliques`], [`prune_cliques_reported`]): greedily
//!    removes the heaviest rows (relations) while the matrix excess exceeds [`EXCESS_FLOOR`].
//!    Dropping a row discards every dependency that uses it, so pruning is a deliberate loss of
//!    candidates. [`PruneReport`] lists exactly the rows lost.
//!
//! 2. **Column elimination** ([`merge_columns`], [`eliminate_column`]): removes low-weight
//!    non-obstruction columns by replacing the rows that contain them with XOR combinations.
//!    It is lossless: no dependency is created or destroyed (see below).
//!
//! # Why elimination preserves the kernel
//!
//! A dependency is a set of rows summing to zero in every column, so a kernel vector `x`
//! selects an *even* number of the rows incident to any column. For a column `c`:
//!
//! - **Weight two**, incident rows `r0, r1`. Evenness forces `x0 = x1`, so the kernel vectors
//!   are exactly those that use the single row `r0 ⊕ r1`, with column `c` cancelled.
//! - **Weight three**, incident rows `r0, r1, r2`. Evenness allows `(x0, x1, x2)` in
//!   `{000, 110, 101, 011}`, which is `(y1 + y2, y1, y2)` for free `y1, y2`: exactly the
//!   combinations of the two rows `r0 ⊕ r1` and `r0 ⊕ r2`, with column `c` cancelled in both.
//!   Replacing the three rows by their *total* XOR instead would keep column `c` set and lose
//!   dependencies: three copies of the row `[1]` have a two-dimensional kernel, whereas the
//!   single row `[1]` has none.
//!
//! In both cases the map from the reduced kernel to the old kernel is a linear bijection, so
//! eliminated dependencies lift exactly. It is the identity on rows that were not touched, and
//! its composition over many steps is [`SparseMatrix::lift_dependency`] through the provenance.
//!
//! Because the weight-three replacement rows `r0 ⊕ r1` and `r0 ⊕ r2` share `r0`, provenance of
//! different rows overlaps; XOR-combining those rows later cancels the shared originals. Columns
//! and provenance therefore both use symmetric difference, which keeps every row equal to the
//! GF(2) sum of the original rows its provenance names. Obstruction columns and always-zero
//! columns are never eliminated, but are XORed along with every other column, so they remain
//! consistent with provenance.
//!
//! # Ordering
//!
//! In production NFS, merges are ordered by a weight-cost heuristic that prefers merges
//! minimising the increase in total matrix weight (Cavallar 1996). At toy scale any order
//! gives a tractable matrix, so the saving is invisible. [`merge_columns`] processes columns in
//! increasing column index, weight two first and then weight three, and takes the lightest
//! incident row as the weight-three pivot (lowest row index on ties), which keeps the new rows
//! sparse. Higher-weight columns are left for the linear algebra step.

use std::fmt;

use crate::filter::matrix::{MatrixRow, SparseMatrix, EXCESS_FLOOR};

// ─── prune_cliques ────────────────────────────────────────────────────────────

/// The rows discarded by clique pruning: the only lossy step of the filter.
///
/// Every dependency of the pre-pruning matrix that uses one of these rows is no longer a
/// candidate. Dependencies that use none of them are all preserved.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PruneReport {
    /// The discarded rows in removal order, each with its provenance.
    pub removed: Vec<MatrixRow>,
    /// The matrix excess before pruning.
    pub excess_before: isize,
    /// The matrix excess after pruning.
    pub excess_after: isize,
}

impl PruneReport {
    /// True iff pruning discarded no rows, so no candidate was lost.
    #[must_use]
    pub fn is_lossless(&self) -> bool {
        self.removed.is_empty()
    }
}

/// Clique/excess pruning: greedily remove heavy rows while excess > EXCESS_FLOOR.
///
/// Equivalent to [`prune_cliques_reported`] with the report discarded. Callers that need to know
/// which candidates were lost should use that function.
///
/// :param matrix: The singleton-removed matrix (consumed).
/// :returns: The pruned matrix with `excess() >= EXCESS_FLOOR` (or an empty or zero-weight one).
pub fn prune_cliques(matrix: SparseMatrix) -> SparseMatrix {
    prune_cliques_reported(matrix).0
}

/// Clique/excess pruning that reports the rows it discards.
///
/// The "graph view": relations are edges, primes (non-obstruction columns) are nodes. A clique
/// is a set of relations sharing a common prime. This greedy heuristic approximates clique
/// pruning by removing the highest-weight row (most columns set) as long as
/// `matrix.excess() > EXCESS_FLOOR as isize`. Ties in weight are broken by row index (lowest
/// index removed first) for determinism. It stops when `excess() == EXCESS_FLOOR`, when no rows
/// remain, or when every remaining row is empty.
///
/// Unlike elimination, pruning is not kernel-preserving: it keeps exactly the dependencies
/// supported on the surviving rows. The report names the discarded rows.
///
/// :param matrix: The singleton-removed matrix (consumed).
/// :returns: The pruned matrix, and the report of discarded rows.
pub fn prune_cliques_reported(mut matrix: SparseMatrix) -> (SparseMatrix, PruneReport) {
    let excess_before = matrix.excess();
    let mut removed = Vec::new();

    while matrix.excess() > EXCESS_FLOOR as isize {
        let max_weight = matrix.rows.iter().map(|r| r.cols.len()).max().unwrap_or(0);
        if max_weight == 0 {
            // Nothing but empty rows (or no rows): removing them would not reduce excess usefully.
            break;
        }
        let Some(row_idx) = matrix.rows.iter().position(|r| r.cols.len() == max_weight) else {
            break;
        };
        removed.push(matrix.rows[row_idx].clone());
        matrix.remove_row(row_idx);
    }

    let excess_after = matrix.excess();
    (matrix, PruneReport { removed, excess_before, excess_after })
}

// ─── elimination ──────────────────────────────────────────────────────────────

/// Why [`eliminate_column`] declined to change the matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EliminationRefusal {
    /// The column index is not below the length of `col_weights`.
    OutOfRange,
    /// The column lies in the obstruction block, which is never eliminated.
    Obstruction,
    /// The column weight is not two or three. Carries the weight found.
    UnsupportedWeight(u32),
    /// `col_weights` disagrees with the number of rows that actually contain the column.
    InconsistentWeights,
}

impl fmt::Display for EliminationRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange => f.write_str("column index is out of range"),
            Self::Obstruction => f.write_str("obstruction columns are never eliminated"),
            Self::UnsupportedWeight(w) => {
                write!(f, "column weight {w} is not 2 or 3, so it cannot be eliminated")
            }
            Self::InconsistentWeights => {
                f.write_str("column weight table disagrees with the rows of the matrix")
            }
        }
    }
}

impl std::error::Error for EliminationRefusal {}

/// Column elimination: eliminate low-weight non-obstruction columns by XOR-combining rows.
///
/// Performs two passes over the non-obstruction columns, each in increasing column index:
///
/// - **Weight two:** the two incident rows are replaced by their XOR.
/// - **Weight three:** the incident rows are replaced by `pivot ⊕ other` for each of the two
///   other rows, the pivot being the lightest (see [`eliminate_column`]).
///
/// A column is re-checked when its turn comes, because earlier steps in the pass change weights;
/// columns whose weight is no longer the pass's target are skipped, and columns that become
/// weight two or three only during the weight-three pass are left alone. The step is lossless:
/// every dependency of the input lifts to one of the output and vice versa, through provenance.
/// It does not enforce `EXCESS_FLOOR`; eliminating columns changes what excess means.
///
/// :param matrix: The pruned matrix (consumed).
/// :returns: The column-eliminated matrix.
pub fn merge_columns(mut matrix: SparseMatrix) -> SparseMatrix {
    merge_pass(&mut matrix, 2);
    merge_pass(&mut matrix, 3);
    matrix
}

/// Eliminate one non-obstruction column of weight two or three, keeping every dependency.
///
/// Rows containing the column are removed, and the replacement rows are appended in order:
///
/// - weight two, rows `r0, r1` (ascending index): the single row `r0 ⊕ r1`;
/// - weight three: with the lightest incident row as pivot `p` (fewest columns, lowest index on
///   ties) and the other two rows `a, b` in ascending index order, the rows `p ⊕ a` then
///   `p ⊕ b`.
///
/// Columns and provenance both take the symmetric difference, so overlapping provenance cancels
/// exactly. The column vanishes from every row, `col_weights` stays equal to the true column
/// weights, and the map from reduced to old dependencies is a bijection (see the module
/// documentation).
///
/// # Errors
///
/// [`EliminationRefusal`] when the column is out of range, in the obstruction block, not of
/// weight two or three, or when `col_weights` disagrees with the rows. The matrix is unchanged
/// on error.
///
/// # Panics
///
/// If a row names a column beyond `col_weights.len()`, which no constructor in this crate
/// produces.
pub fn eliminate_column(
    matrix: &mut SparseMatrix,
    col: usize,
) -> Result<(), EliminationRefusal> {
    let Some(&weight) = matrix.col_weights.get(col) else {
        return Err(EliminationRefusal::OutOfRange);
    };
    if col >= matrix.obstruction_col_start {
        return Err(EliminationRefusal::Obstruction);
    }
    if weight != 2 && weight != 3 {
        return Err(EliminationRefusal::UnsupportedWeight(weight));
    }

    // Ascending row indices of the rows containing `col`.
    let incident: Vec<usize> = matrix
        .rows
        .iter()
        .enumerate()
        .filter(|(_, r)| r.cols.binary_search(&col).is_ok())
        .map(|(i, _)| i)
        .collect();
    if incident.len() != weight as usize {
        return Err(EliminationRefusal::InconsistentWeights);
    }

    let replacements = match incident[..] {
        [i0, i1] => vec![matrix.rows[i0].xor_merge(&matrix.rows[i1])],
        [i0, i1, i2] => {
            // Lightest row is the pivot; strict `<` keeps the lowest index on ties.
            let mut pivot = i0;
            for i in [i1, i2] {
                if matrix.rows[i].cols.len() < matrix.rows[pivot].cols.len() {
                    pivot = i;
                }
            }
            incident
                .iter()
                .filter(|&&i| i != pivot)
                .map(|&i| matrix.rows[pivot].xor_merge(&matrix.rows[i]))
                .collect()
        }
        _ => return Err(EliminationRefusal::InconsistentWeights),
    };

    // Remove from the highest index down so lower indices stay valid; `remove_row` decrements
    // the weights of the departing rows' columns.
    for &i in incident.iter().rev() {
        matrix.remove_row(i);
    }
    for row in replacements {
        for &c in &row.cols {
            matrix.col_weights[c] += 1;
        }
        matrix.rows.push(row);
    }
    Ok(())
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Perform one elimination pass over the non-obstruction columns of the given weight.
///
/// Collects the columns currently at `target_weight`, then eliminates them in column-index
/// order, skipping any whose weight has since changed.
///
/// :param matrix: The matrix to reduce in place.
/// :param target_weight: The column weight to target (2 or 3).
fn merge_pass(matrix: &mut SparseMatrix, target_weight: u32) {
    let candidates: Vec<usize> = matrix
        .col_weights
        .iter()
        .enumerate()
        .take(matrix.obstruction_col_start)
        .filter(|&(_, &w)| w == target_weight)
        .map(|(c, _)| c)
        .collect();

    for col in candidates {
        if matrix.col_weights[col] == target_weight {
            // A refusal here means a malformed weight table; the column is left in place.
            let _ = eliminate_column(matrix, col);
        }
    }
}
