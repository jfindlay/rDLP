//! Nullspace vector representation: the linear algebra → square root seam.
//!
//! `KernelVector` represents a vector in the left nullspace of the filtered matrix as a
//! sorted, deduplicated list of row indices. The square root step expands this to original
//! relation indices via the provenance map.

use std::collections::BTreeMap;

use crate::filter::SparseMatrix;

// ─── KernelVector ─────────────────────────────────────────────────────────────

/// A vector in the left nullspace of the matrix: a subset of rows whose GF(2) sum is zero.
///
/// Representation: a sorted, deduplicated `Vec<usize>` of **filtered-matrix row indices**
/// (indices into `SparseMatrix::rows`). The square root step expands this to original
/// relation indices by collecting `matrix.rows[i].provenance` for each `i` in `row_indices`
/// and taking the symmetric difference (XOR union).
///
/// # Why row indices, not a bit-mask
///
/// - The square root step needs row indices to look up provenance; a bit-mask would require a scan.
/// - Kernel vectors are sparse (typically a small fraction of rows); a bit-mask wastes space.
/// - Solvers (Lanczos, Wiedemann) internally work with bit-packed block vectors, but they
///   convert to `KernelVector` on output — the conversion is O(rows) and happens once per
///   kernel vector, not in the inner loop.
///
/// # F_ℓ generalisation note
///
/// For F_ℓ, a kernel vector is still a subset of rows (those with nonzero coefficient in
/// the nullspace vector). The representation is identical; the F_ℓ extension may add a
/// `coefficients` field (`Vec<Scalar>`) for the non-GF(2) case, but the row-index spine
/// is stable.
///
/// # Invariants
///
/// - `row_indices` is sorted and deduplicated.
/// - Each index is < `matrix.rows.len()` (the filtered matrix, not the original relations).
/// - The GF(2) sum of the selected rows is the zero vector (verified by `verify`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KernelVector {
    /// Sorted, deduplicated row indices into the filtered matrix.
    pub row_indices: Vec<usize>,
}

impl KernelVector {
    /// Construct from a sorted, deduplicated list of row indices.
    ///
    /// # Panics
    ///
    /// Panics if `row_indices` is not sorted or contains duplicates.
    #[must_use]
    pub fn new(row_indices: Vec<usize>) -> Self {
        assert!(
            row_indices.windows(2).all(|w| w[0] < w[1]),
            "KernelVector::new: row_indices must be sorted and deduplicated"
        );
        Self { row_indices }
    }

    /// Construct from a bit-mask over rows (used by solvers internally).
    ///
    /// `mask[i]` is true iff row `i` is in the kernel vector.
    #[must_use]
    pub fn from_mask(mask: &[bool]) -> Self {
        let row_indices = mask.iter().enumerate().filter_map(|(i, &b)| if b { Some(i) } else { None }).collect();
        Self { row_indices }
    }

    /// Verify that this is a valid left-nullspace vector of the given matrix.
    ///
    /// Returns `true` iff the GF(2) sum of `matrix.rows[i].cols` for `i` in `row_indices`
    /// is the empty set (all columns cancel).
    #[must_use]
    pub fn verify(&self, matrix: &SparseMatrix) -> bool {
        // XOR all selected rows' column sets; the result should be empty.
        // Use a bit-set over columns for efficiency.
        let num_cols = matrix.num_cols;
        // Use a Vec<u64> as a bit-set.
        let words = (num_cols + 63) / 64;
        let mut bits = vec![0u64; words];
        for &i in &self.row_indices {
            if i >= matrix.rows.len() {
                return false;
            }
            for &c in &matrix.rows[i].cols {
                bits[c / 64] ^= 1u64 << (c % 64);
            }
        }
        bits.iter().all(|&w| w == 0)
    }

    /// Expand through the provenance map to original relation indices.
    ///
    /// Returns the symmetric difference (XOR union) of `matrix.rows[i].provenance` for
    /// each `i` in `row_indices`. This is the set of original relations whose product
    /// yields a congruence of squares.
    ///
    /// The result is sorted and deduplicated.
    #[must_use]
    pub fn expand_provenance(&self, matrix: &SparseMatrix) -> Vec<usize> {
        // Collect all provenance indices and compute the symmetric difference.
        // Since each provenance list is sorted and deduplicated, we XOR them pairwise.
        let mut result: Vec<usize> = Vec::new();
        for &i in &self.row_indices {
            let prov = &matrix.rows[i].provenance;
            // Symmetric difference of result and prov (both sorted).
            let mut merged = Vec::with_capacity(result.len() + prov.len());
            let mut ri = 0;
            let mut pi = 0;
            while ri < result.len() && pi < prov.len() {
                match result[ri].cmp(&prov[pi]) {
                    std::cmp::Ordering::Less => {
                        merged.push(result[ri]);
                        ri += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        merged.push(prov[pi]);
                        pi += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        // Symmetric difference: both present → cancel.
                        ri += 1;
                        pi += 1;
                    }
                }
            }
            merged.extend_from_slice(&result[ri..]);
            merged.extend_from_slice(&prov[pi..]);
            result = merged;
        }
        result
    }

    /// Number of selected rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.row_indices.len()
    }

    /// True if no rows are selected (the trivial kernel vector).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.row_indices.is_empty()
    }
}

// ─── Bounded combinations of kernel vectors ──────────────────────────────────

/// Budgets for [`search_kernel_combinations`], passed by value.
///
/// Both budgets are enforced before the step they would pay for, so the search never performs
/// more than `max_xor_operations` merges and never emits more than `max_candidates` vectors,
/// whatever the number of input vectors. The number of distinct combinations is `2^rank − 1`,
/// which is never enumerated unless the budgets allow it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelSearchLimits {
    /// Maximum XOR operations: merges of two sorted row-index lists, in both the reduction of
    /// the inputs to an independent set and the construction of combinations.
    pub max_xor_operations: u64,
    /// Maximum distinct nonzero candidates emitted.
    pub max_candidates: u64,
}

/// The budget that ended a search early.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelSearchResource {
    /// The XOR-operation budget was exhausted.
    XorOperations,
    /// The candidate-output budget was exhausted with further combinations remaining.
    Candidates,
}

/// The stage of the search that was running when a budget ended it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelSearchPhase {
    /// Reducing the inputs to an independent set. No candidate has been emitted, and the rank is
    /// not yet known.
    Reduction,
    /// Emitting combinations of the independent set.
    Enumeration,
}

/// Whether a search covered every distinct nonzero combination of its inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KernelSearchStatus {
    /// Every distinct nonzero XOR combination of the inputs was emitted: `2^rank − 1` of them.
    Complete,
    /// A budget ended the search; combinations remain unexamined. This is never a statement
    /// about the kernel or about the unexamined combinations.
    ResourceLimited {
        /// The budget that ran out.
        resource: KernelSearchResource,
        /// The stage that was running.
        phase: KernelSearchPhase,
        /// XOR operations performed so far.
        xor_operations: u64,
        /// Candidates emitted so far.
        candidates: u64,
    },
}

impl KernelSearchStatus {
    /// True iff the search is [`KernelSearchStatus::Complete`].
    #[must_use]
    pub fn is_complete(&self) -> bool {
        matches!(self, Self::Complete)
    }
}

/// The result of [`search_kernel_combinations`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KernelSearch {
    /// Distinct nonzero combinations, each with strictly ascending row indices, in emission
    /// order (see [`search_kernel_combinations`]).
    pub candidates: Vec<KernelVector>,
    /// Whether `candidates` is every combination.
    pub status: KernelSearchStatus,
    /// The GF(2) rank of the inputs when the reduction finished; while it was cut short
    /// (`phase == Reduction`), the number of independent inputs found so far.
    pub rank: usize,
}

/// Enumerate the distinct nonzero XOR combinations of `vectors` within `limits`.
///
/// Each input is read as a GF(2) vector over row indices. An index listed an even number of times
/// cancels, so the input contract is "any `Vec<usize>`": nothing is validated against a matrix,
/// index values are labels (an out-of-range index is carried through unchanged), and no input
/// makes the search panic. A vector that normalizes to the empty set is dropped.
///
/// # Method
///
/// 1. **Reduction.** The inputs are scanned in order and each one that is linearly independent
///    of those kept before it joins an independent set; a later duplicate, a zero vector, or a
///    vector that is the XOR of earlier ones is dropped. The kept vectors are the original
///    inputs, not echelon rows, so every emitted basis vector is one the caller supplied. Because
///    the kept vectors are independent, every nonempty subset of them XORs to a distinct nonzero
///    vector: deduplication and zero-cancellation need no seen-set, and `Complete` means exactly
///    `2^rank − 1` candidates.
/// 2. **Enumeration.** Subsets are produced by ascending size and, within a size, in
///    lexicographic order over positions in the independent set. All independent vectors come
///    first, then every pair, and so on. A small budget therefore covers the caller's own
///    vectors and their low-order sums before long sums, and an output of any length is a prefix
///    of the output under a larger budget. Subsets are generated lazily, one at a time; only
///    `rank` positions of state are held, never the `2^rank` subsets.
///
/// # Budgets
///
/// The XOR budget is charged one unit per merge, before the merge, in both stages. A subset of
/// size `s` costs `s − 1` merges, built from scratch, and a single vector costs none. The output
/// budget is checked before building the next candidate, so a budget exactly equal to the
/// number of combinations still reports `Complete`. When both are exhausted at the same step the
/// output budget is reported. There is no resume cursor: retry with larger limits, which
/// reproduces the earlier output as a prefix.
#[must_use]
pub fn search_kernel_combinations(
    vectors: &[KernelVector],
    limits: KernelSearchLimits,
) -> KernelSearch {
    let mut meter = Meter { limit: limits.max_xor_operations, used: 0 };
    let (basis, reduced) = reduce_to_independent(vectors, &mut meter);
    let rank = basis.len();
    if !reduced {
        let status =
            limited(KernelSearchResource::XorOperations, KernelSearchPhase::Reduction, &meter, 0);
        return KernelSearch { candidates: Vec::new(), status, rank };
    }
    let mut candidates: Vec<KernelVector> = Vec::new();
    let mut emitted: u64 = 0;
    for size in 1..=rank {
        let mut positions: Vec<usize> = (0..size).collect();
        loop {
            if emitted >= limits.max_candidates {
                let status = limited(
                    KernelSearchResource::Candidates,
                    KernelSearchPhase::Enumeration,
                    &meter,
                    emitted,
                );
                return KernelSearch { candidates, status, rank };
            }
            let Some(combination) = combine(&basis, &positions, &mut meter) else {
                let status = limited(
                    KernelSearchResource::XorOperations,
                    KernelSearchPhase::Enumeration,
                    &meter,
                    emitted,
                );
                return KernelSearch { candidates, status, rank };
            };
            candidates.push(KernelVector { row_indices: combination });
            emitted += 1;
            if !advance(&mut positions, rank) {
                break;
            }
        }
    }
    KernelSearch { candidates, status: KernelSearchStatus::Complete, rank }
}

/// Running count of XOR operations against a fixed budget.
struct Meter {
    limit: u64,
    used: u64,
}

impl Meter {
    /// Pay for one XOR; false (with nothing charged) when the budget is spent.
    fn charge(&mut self) -> bool {
        if self.used >= self.limit {
            return false;
        }
        self.used += 1;
        true
    }
}

fn limited(
    resource: KernelSearchResource,
    phase: KernelSearchPhase,
    meter: &Meter,
    candidates: u64,
) -> KernelSearchStatus {
    KernelSearchStatus::ResourceLimited {
        resource,
        phase,
        xor_operations: meter.used,
        candidates,
    }
}

/// Sort and cancel index pairs: the strictly ascending set of indices listed an odd number of
/// times.
fn normalize(row_indices: &[usize]) -> Vec<usize> {
    let mut sorted = row_indices.to_vec();
    sorted.sort_unstable();
    let mut odd: Vec<usize> = Vec::with_capacity(sorted.len());
    for index in sorted {
        if odd.last() == Some(&index) {
            odd.pop();
        } else {
            odd.push(index);
        }
    }
    odd
}

/// Symmetric difference of two strictly ascending lists.
fn xor_sorted(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut merged = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                merged.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&a[i..]);
    merged.extend_from_slice(&b[j..]);
    merged
}

/// The inputs that are independent of their predecessors, and whether the scan finished.
///
/// Independence is decided by an echelon table keyed by the largest index: a vector reduces
/// against the table until its largest index is unclaimed (independent) or it is empty
/// (dependent). Each reduction step removes the current largest index, so the loop ends.
fn reduce_to_independent(
    vectors: &[KernelVector],
    meter: &mut Meter,
) -> (Vec<Vec<usize>>, bool) {
    let mut echelon: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut basis: Vec<Vec<usize>> = Vec::new();
    for vector in vectors {
        let original = normalize(&vector.row_indices);
        let mut reduced = original.clone();
        while let Some(&top) = reduced.last() {
            let Some(pivot) = echelon.get(&top) else {
                echelon.insert(top, reduced);
                basis.push(original);
                break;
            };
            if !meter.charge() {
                return (basis, false);
            }
            reduced = xor_sorted(&reduced, pivot);
        }
    }
    (basis, true)
}

/// XOR of `basis[p]` over `positions`, or `None` when the budget runs out first.
///
/// `positions` is nonempty and in range by construction in the caller.
fn combine(basis: &[Vec<usize>], positions: &[usize], meter: &mut Meter) -> Option<Vec<usize>> {
    let mut sum = basis[positions[0]].clone();
    for &position in &positions[1..] {
        if !meter.charge() {
            return None;
        }
        sum = xor_sorted(&sum, &basis[position]);
    }
    Some(sum)
}

/// Step `positions` (strictly ascending, each below `count`) to the next subset of the same size
/// in lexicographic order; false when it was the last.
fn advance(positions: &mut [usize], count: usize) -> bool {
    let size = positions.len();
    let Some(slot) = (0..size).rev().find(|&s| positions[s] < count - size + s) else {
        return false;
    };
    positions[slot] += 1;
    for s in slot + 1..size {
        positions[s] = positions[s - 1] + 1;
    }
    true
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::{MatrixRow, SparseMatrix};

    fn make_matrix(rows: Vec<(Vec<usize>, Vec<usize>)>, num_cols: usize) -> SparseMatrix {
        let col_weights = vec![0u32; num_cols];
        SparseMatrix {
            rows: rows
                .into_iter()
                .map(|(cols, prov)| MatrixRow { cols, provenance: prov })
                .collect(),
            num_cols,
            obstruction_col_start: num_cols,
            obstruction_count: 0,
            col_weights,
        }
    }

    #[test]
    fn new_panics_on_unsorted() {
        let result = std::panic::catch_unwind(|| KernelVector::new(vec![2, 1]));
        assert!(result.is_err(), "new should panic on unsorted input");
    }

    #[test]
    fn new_panics_on_duplicate() {
        let result = std::panic::catch_unwind(|| KernelVector::new(vec![1, 1]));
        assert!(result.is_err(), "new should panic on duplicate input");
    }

    #[test]
    fn from_mask_roundtrip() {
        let mask = vec![true, false, true, false, true];
        let kv = KernelVector::from_mask(&mask);
        assert_eq!(kv.row_indices, vec![0, 2, 4]);
    }

    #[test]
    fn verify_valid_nullspace() {
        // Matrix: row 0 = {0, 1}, row 1 = {1, 2}, row 2 = {0, 2}.
        // XOR of all three: {0,1} XOR {1,2} XOR {0,2} = {} (all cancel).
        let m = make_matrix(
            vec![(vec![0, 1], vec![0]), (vec![1, 2], vec![1]), (vec![0, 2], vec![2])],
            3,
        );
        let kv = KernelVector::new(vec![0, 1, 2]);
        assert!(kv.verify(&m), "XOR of all three rows should be zero");
    }

    #[test]
    fn verify_invalid_nullspace() {
        // Matrix: row 0 = {0, 1}, row 1 = {1, 2}.
        // XOR = {0, 2} ≠ {}.
        let m = make_matrix(
            vec![(vec![0, 1], vec![0]), (vec![1, 2], vec![1])],
            3,
        );
        let kv = KernelVector::new(vec![0, 1]);
        assert!(!kv.verify(&m), "XOR of rows 0 and 1 should not be zero");
    }

    #[test]
    fn expand_provenance_simple() {
        // Row 0 provenance = [0], row 1 provenance = [1].
        // Symmetric difference = [0, 1].
        let m = make_matrix(
            vec![(vec![0, 1], vec![0]), (vec![1, 2], vec![1])],
            3,
        );
        let kv = KernelVector::new(vec![0, 1]);
        let expanded = kv.expand_provenance(&m);
        assert_eq!(expanded, vec![0, 1]);
    }

    #[test]
    fn expand_provenance_cancellation() {
        // Row 0 provenance = [0, 1], row 1 provenance = [1, 2].
        // Symmetric difference = [0, 2] (1 cancels).
        let m = make_matrix(
            vec![(vec![0], vec![0, 1]), (vec![0], vec![1, 2])],
            1,
        );
        let kv = KernelVector::new(vec![0, 1]);
        let expanded = kv.expand_provenance(&m);
        assert_eq!(expanded, vec![0, 2]);
    }
}
