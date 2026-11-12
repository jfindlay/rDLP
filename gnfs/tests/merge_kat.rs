//! Known-answer tests (KATs) for clique/excess pruning and column elimination.
//!
//! The first group below covers the staged pipeline; the second group (after the
//! "kernel-preserving elimination" banner) checks that elimination is lossless and that pruning
//! is reported separately. Both groups use hand-derived fixtures; every exhaustive comparison
//! enumerates GF(2) kernels by brute force, independently of the code under test.
//!
//! First group:
//!
//! - **KAT (a) — 2-way merge correctness**: a weight-2 column is eliminated; the merged
//!   row has the correct cols (symmetric difference) and provenance (union).
//!
//! - **KAT (b) — Determinism**: the full pipeline (build_matrix + remove_singletons +
//!   prune_cliques + merge_columns) is deterministic for a fixed corpus.
//!
//! - **KAT (c) — CADO-NFS oracle**: gated/ignored when CADO is absent.
//!
//! - **KAT (d) — End-to-end provenance**: for each row in the final matrix, the XOR of
//!   the original relations indexed by its provenance equals that row's column set.
//!
//! # Toy setup
//!
//! All KATs use ``f(x) = x³ − x − 1`` (coefficients: [−1, −1, 0, 1] least-significant
//! first) with ``B_rat = B_alg = 13``.
//!
//! Column layout for B_rat = B_alg = 13:
//! - Rational primes ≤ 13: [2, 3, 5, 7, 11, 13] → rat_size = 6.
//! - Algebraic ideals: (5,2) → index 0, (7,5) → index 1, (11,6) → index 2.
//! - alg_size = 3; obstruction_count = 1; matrix_width = 10.
//! - Col 0: p=2, Col 1: p=3, Col 2: p=5, Col 3: p=7, Col 4: p=11, Col 5: p=13.
//! - Col 6: ideal(5,2), Col 7: ideal(7,5), Col 8: ideal(11,6).
//! - Col 9: sign bit (obstruction column 0).

use gnfs::filter::merge::{
    eliminate_column, prune_cliques_reported, EliminationRefusal, PruneReport,
};
use gnfs::filter::{build_checked_matrix, MatrixRow, SparseMatrix};
use gnfs::sieve::{CheckedRelations, RelationChecker};
use gnfs::sqrt::SqrtContext;
use gnfs::{
    build_matrix, merge_columns, prune_cliques, remove_singletons, ExponentVector, FactorBase,
    Relation,
};
use num_bigint::BigInt;
use shared_numfield::IntPoly;

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn bi(n: i64) -> BigInt {
    BigInt::from(n)
}

/// f(x) = x³ − x − 1.
fn f_cubic() -> IntPoly {
    IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)])
}

/// Construct a Relation directly from field values (bypasses smoothness checks).
fn make_relation(
    a: i64,
    b: i64,
    rational_entries: Vec<(usize, u32)>,
    algebraic_entries: Vec<(usize, u32)>,
    rational_sign: bool,
) -> Relation {
    Relation {
        a: bi(a),
        b: bi(b),
        rational_exponents: ExponentVector { entries: rational_entries },
        algebraic_exponents: ExponentVector { entries: algebraic_entries },
        rational_sign,
    }
}

// ─── KAT (a): 2-way merge correctness ────────────────────────────────────────

/// KAT (a): A weight-2 column is eliminated by merging its two containing rows.
///
/// Relations (columns: p=7 → 3, p=11 → 4, p=13 → 5):
/// R0 = {3,4}, R1 = {3,5}, and three copies R2, R3, R4 = {4,5}.
///
/// Column weights by counting: col 3 → 2 (R0, R1), col 4 → 4 (R0, R2, R3, R4), col 5 → 4 (R1,
/// R2, R3, R4). Only column 3 is eliminable (weight two); the weight-4 columns are left alone.
/// R0 ⊕ R1 = {3,4} Δ {3,5} = {4,5} with provenance {0} ∪ {1} = {0,1}; rows go from 5 to 4.
#[test]
fn kat_a_two_way_merge_correctness() {
    let f = f_cubic();
    let fb = FactorBase::new(&f, 13, 13);

    let idx_p7 = fb.rational_index(7).expect("prime 7 in rational base");   // col 3
    let idx_p11 = fb.rational_index(11).expect("prime 11 in rational base"); // col 4
    let idx_p13 = fb.rational_index(13).expect("prime 13 in rational base"); // col 5

    // Verify column indices match our assumptions.
    assert_eq!(idx_p7, 3, "p=7 should be at col 3");
    assert_eq!(idx_p11, 4, "p=11 should be at col 4");
    assert_eq!(idx_p13, 5, "p=13 should be at col 5");

    // R0: {3, 4}, R1: {3, 5} — col 3 weight 2 (only these two rows).
    // R2, R3, R4: {4, 5} — col 4 weight 4, col 5 weight 4 (not weight-2 or weight-3).
    let relations = vec![
        // R0: primes {7, 11} → cols {3, 4}
        make_relation(7, 1, vec![(idx_p7, 1), (idx_p11, 1)], vec![], false),
        // R1: primes {7, 13} → cols {3, 5}
        make_relation(13, 1, vec![(idx_p7, 1), (idx_p13, 1)], vec![], false),
        // R2, R3, R4: primes {11, 13} → cols {4, 5} (three copies to give weight 4 after merge)
        make_relation(11, 1, vec![(idx_p11, 1), (idx_p13, 1)], vec![], false),
        make_relation(17, 1, vec![(idx_p11, 1), (idx_p13, 1)], vec![], false),
        make_relation(19, 1, vec![(idx_p11, 1), (idx_p13, 1)], vec![], false),
    ];

    let matrix = build_matrix(&relations, &fb);

    // Verify initial col_weights.
    assert_eq!(matrix.rows.len(), 5, "initial: 5 rows");
    assert_eq!(matrix.col_weights[idx_p7], 2, "col 3 (p=7) initial weight = 2");
    assert_eq!(matrix.col_weights[idx_p11], 4, "col 4 (p=11) initial weight = 4");
    assert_eq!(matrix.col_weights[idx_p13], 4, "col 5 (p=13) initial weight = 4");

    // remove_singletons is a no-op here (no weight-1 columns).
    let after_singletons = remove_singletons(matrix);
    assert_eq!(after_singletons.rows.len(), 5, "after singleton removal: still 5 rows");

    // prune_cliques: excess = 5 - (10 - 1) = -4 < EXCESS_FLOOR (20) → no pruning.
    let after_prune = prune_cliques(after_singletons);
    assert_eq!(after_prune.rows.len(), 5, "prune_cliques: no rows removed (excess < EXCESS_FLOOR)");

    // merge_columns: weight-2 pass processes col 3 (only weight-2 column).
    // R0 and R1 are merged: new_row.cols = {4, 5}, new_row.provenance = {0, 1}.
    // col 4 and col 5 have weight 4 (not 2 or 3) → not merged.
    // Weight-3 pass: no weight-3 columns → no-op.
    // Final: 4 rows (R2, R3, R4, merged_row). Row count 5→4.
    let after_merge = merge_columns(after_prune);

    // Row count decreased by 1 (two rows merged into one, net -1).
    assert_eq!(
        after_merge.rows.len(),
        4,
        "merge should reduce row count by 1: 5 → 4 (R0+R1 merged, R2/R3/R4 unchanged)"
    );

    // Column 3 (p=7) should now have weight 0 (eliminated by the 2-way merge).
    assert_eq!(
        after_merge.col_weights[idx_p7], 0,
        "col 3 (p=7) should have weight 0 after merge (eliminated)"
    );

    // Find the merged row: cols = {4, 5} (symmetric difference of {3,4} and {3,5}),
    // provenance = {0, 1} (union of R0=[0] and R1=[1]).
    let merged_row = after_merge.rows.iter().find(|r| r.provenance == vec![0usize, 1usize]);
    assert!(merged_row.is_some(), "merged row with provenance [0, 1] should exist");

    let merged_row = merged_row.unwrap();

    // Merged row cols = symmetric difference of {3,4} and {3,5} = {4,5}.
    assert_eq!(
        merged_row.cols,
        vec![idx_p11, idx_p13],
        "merged row cols should be {{4, 5}} (col 3 cancelled by XOR)"
    );
    assert!(!merged_row.cols.contains(&idx_p7), "merged row should not contain col 3 (p=7)");

    // Merged row provenance = union of [0] and [1] = [0, 1].
    assert_eq!(
        merged_row.provenance,
        vec![0usize, 1usize],
        "merged row provenance should be [0, 1] (union of R0 and R1)"
    );

    // Verify col_weights consistency: col_weights[c] should equal the number of rows containing c.
    for col in 0..after_merge.num_cols {
        let actual_weight =
            after_merge.rows.iter().filter(|r| r.cols.binary_search(&col).is_ok()).count();
        assert_eq!(
            after_merge.col_weights[col] as usize,
            actual_weight,
            "col_weights[{col}] should be consistent with actual row contents"
        );
    }
}

// ─── KAT (b): Determinism ────────────────────────────────────────────────────

/// KAT (b): The full pipeline is deterministic for a fixed corpus.
///
/// Calls build_matrix + remove_singletons + prune_cliques + merge_columns twice on the
/// same input and asserts identical matrix dimensions and col_weights.
#[test]
fn kat_b_determinism() {
    let f = f_cubic();
    let fb = FactorBase::new(&f, 13, 13);

    let idx_p2 = fb.rational_index(2).unwrap();
    let idx_p3 = fb.rational_index(3).unwrap();
    let idx_p5 = fb.rational_index(5).unwrap();
    let idx_p7 = fb.rational_index(7).unwrap();
    let idx_p11 = fb.rational_index(11).unwrap();
    let idx_p13 = fb.rational_index(13).unwrap();
    let idx_alg_5_2 = fb.algebraic_index(5, 2).unwrap();
    let idx_alg_7_5 = fb.algebraic_index(7, 5).unwrap();

    // A corpus with enough relations to survive singleton removal and have mergeable columns.
    let relations = vec![
        make_relation(1, 1, vec![(idx_p2, 1), (idx_p3, 1)], vec![], false),
        make_relation(3, 1, vec![(idx_p3, 1), (idx_p5, 1)], vec![], false),
        make_relation(5, 1, vec![(idx_p5, 1), (idx_p7, 1)], vec![], false),
        make_relation(7, 1, vec![(idx_p7, 1), (idx_p11, 1)], vec![(idx_alg_5_2, 1)], false),
        make_relation(11, 1, vec![(idx_p11, 1), (idx_p13, 1)], vec![(idx_alg_7_5, 1)], true),
        make_relation(13, 1, vec![(idx_p7, 1), (idx_p13, 1)], vec![(idx_alg_5_2, 1)], true),
        make_relation(17, 1, vec![(idx_p2, 1), (idx_p7, 1)], vec![], false),
        make_relation(19, 1, vec![(idx_p3, 1), (idx_p11, 1)], vec![], false),
    ];

    // First run.
    let m1 = build_matrix(&relations, &fb);
    let m1 = remove_singletons(m1);
    let m1 = prune_cliques(m1);
    let m1 = merge_columns(m1);

    // Second run (same input).
    let m2 = build_matrix(&relations, &fb);
    let m2 = remove_singletons(m2);
    let m2 = prune_cliques(m2);
    let m2 = merge_columns(m2);

    // Assert identical dimensions.
    assert_eq!(
        m1.rows.len(), m2.rows.len(),
        "determinism: row count must be identical across runs"
    );
    assert_eq!(
        m1.num_cols, m2.num_cols,
        "determinism: num_cols must be identical"
    );

    // Assert identical col_weights.
    assert_eq!(
        m1.col_weights, m2.col_weights,
        "determinism: col_weights must be identical across runs"
    );

    // Assert identical row contents (cols and provenance).
    for i in 0..m1.rows.len() {
        assert_eq!(
            m1.rows[i].cols, m2.rows[i].cols,
            "determinism: row {i} cols must be identical"
        );
        assert_eq!(
            m1.rows[i].provenance, m2.rows[i].provenance,
            "determinism: row {i} provenance must be identical"
        );
    }
}

// ─── KAT (c): CADO-NFS comparison (unimplemented placeholder) ────────────────

/// KAT (c): placeholder for a comparison against CADO-NFS; it panics if run.
///
/// The body is `unimplemented!()`: it invokes no external tool and checks nothing. A real test
/// would:
/// 1. Run the full GNFS filtering pipeline at matched parameters.
/// 2. Compare the filtered matrix dimensions (row count, column count, total weight)
///    against CADO's output within a tolerance (e.g., ±10% on row count).
///
/// It is ignored so that the default test run skips it.
#[test]
#[ignore = "unimplemented placeholder: the body is unimplemented!(), so this panics if run"]
fn kat_c_cado_nfs_oracle() {
    // Placeholder: when CADO-NFS is available, compare filtered matrix dimensions
    // against CADO's output at matched parameters (same polynomial, same factor base
    // bounds, same relation corpus). The filtered matrix row count, column count, and
    // total Hamming weight should be within ±10% of CADO's output.
    //
    // Steps (when implementing):
    // 1. Generate a relation corpus using the same polynomial and bounds as a CADO run.
    // 2. Run build_matrix + remove_singletons + prune_cliques + merge_columns.
    // 3. Parse CADO's filtered matrix output (e.g., from .mat file).
    // 4. Assert dimensions are within tolerance.
    unimplemented!("CADO-NFS oracle not yet implemented — requires CADO installation");
}

// ─── KAT (d): End-to-end provenance ──────────────────────────────────────────

/// KAT (d): End-to-end provenance correctness.
///
/// For each row in the final matrix, collect its provenance indices, look up the
/// original relations, reconstruct their GF(2) column sets, XOR them together, and
/// verify the result equals the final row's cols.
///
/// This confirms that the provenance map is correct end-to-end: expanding a row through
/// its provenance recovers the original relations whose GF(2) parities XOR to that row.
///
/// Setup: a corpus where singleton removal and merging both occur, so provenance sets
/// are non-trivial (some rows have provenance of length > 1 after merging).
#[test]
fn kat_d_end_to_end_provenance() {
    let f = f_cubic();
    let fb = FactorBase::new(&f, 13, 13);

    let idx_p7 = fb.rational_index(7).unwrap();   // col 3
    let idx_p11 = fb.rational_index(11).unwrap(); // col 4
    let idx_p13 = fb.rational_index(13).unwrap(); // col 5
    let idx_alg_5_2 = fb.algebraic_index(5, 2).unwrap(); // col 6
    let idx_alg_7_5 = fb.algebraic_index(7, 5).unwrap(); // col 7

    // Corpus: 5 relations, no singletons, with weight-2 columns for merging.
    // R0: cols {3, 4}         (p=7, p=11)
    // R1: cols {3, 5}         (p=7, p=13)
    // R2: cols {4, 5}         (p=11, p=13)
    // R3: cols {3, 6}         (p=7, ideal(5,2))
    // R4: cols {5, 6, 7}      (p=13, ideal(5,2), ideal(7,5))
    //
    // col_weights after build:
    //   col 3 (p=7): 3 (R0, R1, R3)
    //   col 4 (p=11): 2 (R0, R2)
    //   col 5 (p=13): 3 (R1, R2, R4)
    //   col 6 (ideal(5,2)): 2 (R3, R4)
    //   col 7 (ideal(7,5)): 1 (R4) ← singleton!
    //
    // After remove_singletons: R4 removed (col 7 is singleton → R4 removed).
    //   After R4 removed: col 5 weight = 2 (R1, R2), col 6 weight = 1 (R3) → R3 removed.
    //   After R3 removed: col 3 weight = 2 (R0, R1). Fixpoint.
    // Surviving: R0, R1, R2.
    //
    // After prune_cliques: excess = 3 - (10 - 1) = -6 < EXCESS_FLOOR → no pruning.
    //
    // After merge_columns: col 3 (weight 2, R0 and R1) → merge.
    //   merged: cols = {4, 5}, provenance = {0, 1}.
    //   col 4 (weight 2, R2 and merged) → merge.
    //   merged2: cols = {5}, provenance = {0, 1, 2}.
    //   col 5 (weight 1 now) → singleton in col_weights but merge_columns doesn't enforce floor.
    //   Actually after the weight-2 pass, col 5 has weight 1 in merged2 only.
    //   The weight-3 pass: no weight-3 columns remain.
    // Final matrix: 1 row with cols {5}, provenance {0, 1, 2}.

    let relations = vec![
        // R0: primes {7, 11} → cols {3, 4}
        make_relation(7, 1, vec![(idx_p7, 1), (idx_p11, 1)], vec![], false),
        // R1: primes {7, 13} → cols {3, 5}
        make_relation(13, 1, vec![(idx_p7, 1), (idx_p13, 1)], vec![], false),
        // R2: primes {11, 13} → cols {4, 5}
        make_relation(11, 1, vec![(idx_p11, 1), (idx_p13, 1)], vec![], false),
        // R3: primes {7}, ideal(5,2) → cols {3, 6}
        make_relation(5, 1, vec![(idx_p7, 1)], vec![(idx_alg_5_2, 1)], false),
        // R4: primes {13}, ideal(5,2), ideal(7,5) → cols {5, 6, 7}
        make_relation(3, 1, vec![(idx_p13, 1)], vec![(idx_alg_5_2, 1), (idx_alg_7_5, 1)], false),
    ];

    // Build the initial matrix (keep a clone for provenance verification).
    let initial_matrix = build_matrix(&relations, &fb);

    // Record the GF(2) column set for each original relation (from the initial matrix rows).
    // initial_matrix.rows[i].cols is the GF(2) column set for relation i.
    let original_cols: Vec<Vec<usize>> =
        initial_matrix.rows.iter().map(|r| r.cols.clone()).collect();

    // Run the full pipeline.
    let after_singletons = remove_singletons(initial_matrix);
    let after_prune = prune_cliques(after_singletons);
    let final_matrix = merge_columns(after_prune);

    // For each row in the final matrix, verify provenance.
    for (row_idx, row) in final_matrix.rows.iter().enumerate() {
        // XOR the GF(2) column sets of all original relations in the provenance.
        let mut expected_cols: Vec<usize> = Vec::new();
        for &orig_idx in &row.provenance {
            // XOR expected_cols with original_cols[orig_idx] (symmetric difference).
            let orig = &original_cols[orig_idx];
            let mut merged = Vec::new();
            let mut i = 0;
            let mut j = 0;
            while i < expected_cols.len() && j < orig.len() {
                match expected_cols[i].cmp(&orig[j]) {
                    std::cmp::Ordering::Less => {
                        merged.push(expected_cols[i]);
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        merged.push(orig[j]);
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        // XOR cancels.
                        i += 1;
                        j += 1;
                    }
                }
            }
            merged.extend_from_slice(&expected_cols[i..]);
            merged.extend_from_slice(&orig[j..]);
            expected_cols = merged;
        }

        assert_eq!(
            row.cols, expected_cols,
            "row {row_idx}: XOR of provenance relations should equal the row's cols"
        );
    }

    // Verify col_weights consistency in the final matrix.
    for col in 0..final_matrix.num_cols {
        let actual_weight =
            final_matrix.rows.iter().filter(|r| r.cols.binary_search(&col).is_ok()).count();
        assert_eq!(
            final_matrix.col_weights[col] as usize,
            actual_weight,
            "final matrix col_weights[{col}] should be consistent with row contents"
        );
    }
}

// ─── Additional: prune_cliques excess floor ───────────────────────────────────

/// Verify that prune_cliques respects EXCESS_FLOOR.
///
/// Build a matrix with excess well above EXCESS_FLOOR, run prune_cliques, and assert
/// that the resulting excess is exactly EXCESS_FLOOR (or as close as possible without
/// going below).
#[test]
fn kat_prune_cliques_respects_excess_floor() {
    use gnfs::EXCESS_FLOOR;

    let f = f_cubic();
    let fb = FactorBase::new(&f, 13, 13);

    // matrix_width = 10, obstruction_count = 1.
    // excess = rows - (10 - 1) = rows - 9.
    // To have excess > EXCESS_FLOOR (20), we need rows > 29.
    // Build 35 rows: each row has a unique pair of columns from {3, 4, 5} plus a
    // unique "heavy" column that appears only in that row (to give varying weights).
    // Actually, let's just build many rows with shared columns.

    let idx_p7 = fb.rational_index(7).unwrap();   // col 3
    let idx_p11 = fb.rational_index(11).unwrap(); // col 4
    let idx_p13 = fb.rational_index(13).unwrap(); // col 5

    // Build 35 relations, all sharing cols {3, 4, 5} (weight 3 each).
    // excess = 35 - 9 = 26 > 20 = EXCESS_FLOOR.
    let mut relations = Vec::new();
    for i in 0..35i64 {
        relations.push(make_relation(
            i + 100,
            1,
            vec![(idx_p7, 1), (idx_p11, 1), (idx_p13, 1)],
            vec![],
            false,
        ));
    }

    let matrix = build_matrix(&relations, &fb);
    assert_eq!(matrix.rows.len(), 35, "initial: 35 rows");
    assert_eq!(matrix.excess(), 35 - 9, "initial excess = 26");

    // remove_singletons: no singletons (all cols have weight 35).
    let after_singletons = remove_singletons(matrix);
    assert_eq!(after_singletons.rows.len(), 35, "no singletons removed");

    let after_prune = prune_cliques(after_singletons);

    // After pruning: excess should be exactly EXCESS_FLOOR (20).
    // rows needed = EXCESS_FLOOR + 9 = 29.
    let expected_rows = EXCESS_FLOOR + (10 - 1); // EXCESS_FLOOR + non-obstruction cols
    assert_eq!(
        after_prune.rows.len(),
        expected_rows,
        "after prune_cliques: rows should be EXCESS_FLOOR + non-obstruction cols = {expected_rows}"
    );
    assert_eq!(
        after_prune.excess(),
        EXCESS_FLOOR as isize,
        "after prune_cliques: excess should be exactly EXCESS_FLOOR"
    );
}

// ═══ Kernel-preserving elimination ═══════════════════════════════════════════
//
// Notation. A *dependency* of a matrix is a set of rows whose GF(2) sum is zero in every column,
// i.e. a kernel vector `x` with `x·M = 0`. Rows are written as column sets, e.g. `[0, 2]`; the
// provenance of row `i` of the starting matrix is `[i]`. A reduced row's provenance names the
// starting rows it is the XOR of, and a reduced dependency *lifts* to the XOR of the provenance of
// its rows. Elimination must be lossless: the lifted reduced kernel equals the starting kernel,
// with no vector lost and none counted twice.

// ─── Independent oracles ─────────────────────────────────────────────────────

/// Column bitmask of a row (bit `c` set iff column `c` is set); needs fewer than 64 columns.
fn row_mask(row: &MatrixRow) -> u64 {
    row.cols.iter().fold(0, |m, &c| m | 1u64 << c)
}

/// Bitmask of the starting rows named by a row's provenance.
fn prov_mask(row: &MatrixRow) -> u32 {
    row.provenance.iter().fold(0, |m, &i| m | 1u32 << i)
}

/// A matrix whose last column is the (only) obstruction column and whose row `i` has provenance
/// `[i]`. `col_weights` is counted here, not by the library.
fn matrix_from_masks(rows: &[u64], num_cols: usize) -> SparseMatrix {
    let mut col_weights = vec![0u32; num_cols];
    let rows: Vec<MatrixRow> = rows
        .iter()
        .enumerate()
        .map(|(i, &mask)| {
            let cols: Vec<usize> = (0..num_cols).filter(|&c| mask >> c & 1 == 1).collect();
            for &c in &cols {
                col_weights[c] += 1;
            }
            MatrixRow { cols, provenance: vec![i] }
        })
        .collect();
    SparseMatrix {
        rows,
        num_cols,
        obstruction_col_start: num_cols - 1,
        obstruction_count: 1,
        col_weights,
    }
}

/// Every kernel vector of the rows, as a sorted list of row-selection bitmasks, by exhaustive
/// Gray-code enumeration of all `2^n` selections (`n <= 24`).
fn kernel_by_enumeration(rows: &[u64]) -> Vec<u32> {
    assert!(rows.len() <= 24, "enumeration is exponential in the row count");
    let mut kernel = vec![0u32];
    let (mut selection, mut sum) = (0u32, 0u64);
    for step in 1u32..(1u32 << rows.len()) {
        // Gray code: step k flips bit `trailing_zeros(k)`, visiting every selection once.
        let bit = step.trailing_zeros();
        selection ^= 1 << bit;
        sum ^= rows[bit as usize];
        if sum == 0 {
            kernel.push(selection);
        }
    }
    kernel.sort_unstable();
    kernel
}

/// Lift a reduced selection to the starting rows by XORing provenance masks.
fn lift_mask(reduced: &SparseMatrix, selection: u32) -> u32 {
    reduced
        .rows
        .iter()
        .enumerate()
        .filter(|(i, _)| selection >> i & 1 == 1)
        .fold(0, |acc, (_, row)| acc ^ prov_mask(row))
}

/// `col_weights` equals the true column weights, and every row is a sorted, in-range set.
fn assert_weights_true(m: &SparseMatrix) {
    assert_eq!(m.col_weights.len(), m.num_cols, "weight table has one entry per column");
    for c in 0..m.num_cols {
        let actual = m.rows.iter().filter(|r| r.cols.contains(&c)).count();
        assert_eq!(m.col_weights[c] as usize, actual, "col_weights[{c}] must be the true weight");
    }
    for row in &m.rows {
        assert!(row.cols.windows(2).all(|w| w[0] < w[1]), "cols strictly ascending: {row:?}");
        assert!(row.cols.iter().all(|&c| c < m.num_cols), "cols in range: {row:?}");
        assert!(row.provenance.windows(2).all(|w| w[0] < w[1]), "provenance sorted: {row:?}");
    }
}

/// The reduced matrix is a lossless image of the starting rows `original`.
///
/// 1. Each reduced row is the XOR of the starting rows its provenance names.
/// 2. The reduced kernel and the starting kernel have the same size.
/// 3. The lifted reduced kernel, as a set, equals the starting kernel. With (2), the lift is a
///    bijection: no dependency is lost and none is produced twice.
fn assert_lossless(original: &[u64], reduced: &SparseMatrix) {
    for row in &reduced.rows {
        let expected =
            row.provenance.iter().fold(0u64, |acc, &i| acc ^ original[i]);
        assert_eq!(row_mask(row), expected, "row {row:?} must equal the XOR of its provenance");
    }
    let reduced_masks: Vec<u64> = reduced.rows.iter().map(row_mask).collect();
    let old_kernel = kernel_by_enumeration(original);
    let new_kernel = kernel_by_enumeration(&reduced_masks);
    assert_eq!(new_kernel.len(), old_kernel.len(), "kernel sizes must agree (lift is injective)");
    let mut lifted: Vec<u32> = new_kernel.iter().map(|&y| lift_mask(reduced, y)).collect();
    lifted.sort_unstable();
    assert_eq!(lifted, old_kernel, "lifted reduced kernel must equal the starting kernel");
}

/// Eliminate columns one step at a time to a fixpoint, checking weights after every step.
///
/// Returns `(weight-3 steps taken, steps whose incident rows had overlapping provenance)`, so
/// callers can confirm that the interesting cases were reached.
fn eliminate_stepwise(m: &mut SparseMatrix) -> (usize, usize) {
    let (mut weight_three, mut overlapping_steps) = (0, 0);
    loop {
        let mut progressed = false;
        for col in 0..m.obstruction_col_start {
            let (weight, overlapping) = {
                let incident: Vec<&MatrixRow> =
                    m.rows.iter().filter(|r| r.cols.contains(&col)).collect();
                let overlapping = incident.iter().enumerate().any(|(i, a)| {
                    incident[i + 1..]
                        .iter()
                        .any(|b| a.provenance.iter().any(|p| b.provenance.contains(p)))
                });
                (incident.len(), overlapping)
            };
            if eliminate_column(m, col).is_ok() {
                progressed = true;
                weight_three += usize::from(weight == 3);
                overlapping_steps += usize::from(overlapping);
            }
            assert_weights_true(m);
        }
        if !progressed {
            return (weight_three, overlapping_steps);
        }
    }
}

/// Run both elimination routes on one matrix and check each is lossless.
fn check_matrix(rows: &[u64], num_cols: usize) -> (usize, usize) {
    let start = matrix_from_masks(rows, num_cols);

    let mut stepped = start.clone();
    let stats = eliminate_stepwise(&mut stepped);
    assert_lossless(rows, &stepped);

    let merged = merge_columns(start);
    assert_weights_true(&merged);
    assert_lossless(rows, &merged);
    stats
}

/// Deterministic xorshift64 generator for the sampled sweep.
struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

// ─── Exhaustive and sampled kernel comparisons ───────────────────────────────

/// Every matrix of two shapes: 4 rows × 4 columns and 5 rows × 3 columns (last column is the
/// obstruction block; the others are eliminable). That is 2^16 + 2^15 matrices, including zero
/// rows, repeated rows, zero columns, and every weight-2 and weight-3 configuration. For each,
/// the lifted reduced kernel must equal the brute-force kernel of the starting matrix, through
/// both the stepwise driver and `merge_columns`. The sweep must actually reach weight-3
/// eliminations whose operands have overlapping provenance, or it proves nothing.
#[test]
fn kat_exhaustive_small_matrices_lift_exactly() {
    let (mut weight_three, mut overlapping) = (0, 0);
    for &(num_rows, num_cols) in &[(4usize, 4usize), (5, 3)] {
        let row_bits = num_cols;
        for code in 0u32..(1u32 << (num_rows * row_bits)) {
            let rows: Vec<u64> = (0..num_rows)
                .map(|i| u64::from(code >> (i * row_bits)) & ((1u64 << row_bits) - 1))
                .collect();
            let (w3, ov) = check_matrix(&rows, num_cols);
            weight_three += w3;
            overlapping += ov;
        }
    }
    assert!(weight_three > 1000, "sweep must exercise weight-3 elimination ({weight_three})");
    assert!(overlapping > 100, "sweep must exercise overlapping provenance ({overlapping})");
}

/// Sampled larger matrices (8 × 7 and 10 × 8, density about one third) from a fixed seed. These
/// reach chains of eliminations in which provenance overlaps repeatedly.
#[test]
fn kat_sampled_larger_matrices_lift_exactly() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let mut overlapping = 0;
    for &(num_rows, num_cols, samples) in &[(8usize, 7usize, 3000usize), (10, 8, 1500)] {
        for _ in 0..samples {
            let rows: Vec<u64> = (0..num_rows)
                .map(|_| (0..num_cols).fold(0u64, |m, c| m | u64::from(rng.next() % 3 == 0) << c))
                .collect();
            overlapping += check_matrix(&rows, num_cols).1;
        }
    }
    assert!(overlapping > 1000, "samples must reach overlapping provenance ({overlapping})");
}

// ─── Hand-derived counterexample and provenance cases ────────────────────────

/// Three identical rows `[1]`, provenance `[0]`, `[1]`, `[2]`, in a 3-column matrix whose last
/// column is the obstruction block.
///
/// Hand derivation. A set of identical rows sums to zero iff its size is even, so the kernel is
/// `{∅, {0,1}, {0,2}, {1,2}}` (dimension 2). Column 1 has weight three. Pivot-relative
/// elimination (pivot row 0; all rows tie in weight) gives `r0 ⊕ r1 = []` with provenance
/// `[0,1]` and `r0 ⊕ r2 = []` with provenance `[0,2]`: two empty rows, whose kernel is all four
/// selections, lifting to `∅, {0,1}, {0,2}`, and `{0,1} ⊕ {0,2} = {1,2}`: exactly the original.
/// Replacing the three rows by their total XOR instead gives the single row `[1]`, which has no
/// non-empty dependency, so two of three dependencies (all non-empty ones) would be lost.
#[test]
fn kat_three_identical_rows_keep_their_kernel() {
    let rows = [0b010u64; 3];

    // The premise: the old behaviour (one row [1]) has a trivial kernel; the original has four.
    assert_eq!(kernel_by_enumeration(&rows), vec![0b000, 0b011, 0b101, 0b110]);
    assert_eq!(kernel_by_enumeration(&[0b010]), vec![0]);

    let merged = merge_columns(matrix_from_masks(&rows, 3));
    assert_eq!(
        merged.rows,
        vec![
            MatrixRow { cols: vec![], provenance: vec![0, 1] },
            MatrixRow { cols: vec![], provenance: vec![0, 2] },
        ]
    );
    assert_eq!(merged.col_weights, vec![0, 0, 0]);
    assert_lossless(&rows, &merged);
    // Both empty rows are themselves dependencies, and so is their sum.
    assert_eq!(merged.lift_dependency(&[0, 1]), Some(vec![1, 2]));
}

/// Weight-3 elimination whose incident rows have disjoint provenance.
///
/// Matrix: 4 columns (column 0 eliminable, columns 1–3 obstruction), rows
/// `R0=[0,1] R1=[0,2] R2=[0,3] R3=[1,2]`, provenance `[0] [1] [2] [3]`. Column 0 has weight three
/// in R0, R1, R2 (provenance pairwise disjoint). R0, R1, R2 all have two columns, so the pivot is
/// R0 (lowest index). Replacements: `R0⊕R1 = [1,2]` prov `[0,1]`, `R0⊕R2 = [1,3]` prov `[0,2]`;
/// R3 is untouched. Expected order: `[R3, R0⊕R1, R0⊕R2]`. By hand, `R0 ⊕ R1 ⊕ R3 = 0`, the only
/// dependency, so the reduced dependency `R3 + (R0⊕R1)` lifts to `{0,1,3}`.
#[test]
fn kat_weight_three_elimination_with_disjoint_provenance() {
    let rows = [0b0011u64, 0b0101, 0b1001, 0b0110];
    assert_eq!(kernel_by_enumeration(&rows), vec![0b0000, 0b1011]);

    let mut m = matrix_from_masks(&rows, 4);
    m.obstruction_col_start = 1;
    m.obstruction_count = 3;
    assert_eq!(m.col_weights[0], 3);
    assert_eq!(eliminate_column(&mut m, 0), Ok(()));

    assert_eq!(
        m.rows,
        vec![
            MatrixRow { cols: vec![1, 2], provenance: vec![3] },
            MatrixRow { cols: vec![1, 2], provenance: vec![0, 1] },
            MatrixRow { cols: vec![1, 3], provenance: vec![0, 2] },
        ]
    );
    assert_eq!(m.col_weights, vec![0, 3, 2, 1]);
    assert_weights_true(&m);
    assert_lossless(&rows, &m);
    assert_eq!(m.lift_dependency(&[0, 1]), Some(vec![0, 1, 3]));
}

/// Two successive weight-3 eliminations produce overlapping provenance that must cancel.
///
/// Matrix: 6 columns (columns 0–1 eliminable, 2–5 obstruction), rows
/// `R0=[0] R1=[0,1] R2=[0,1,2] R3=[1,2,3,4] R4=[2,3,4]`, provenance `[0]..[4]`.
///
/// Column 0 has weight 3 (R0, R1, R2); column 1 has weight 3 (R1, R2, R3).
///
/// - Eliminate column 0. Pivot is the lightest incident row, R0. `A = R0⊕R1 = [1]` prov `[0,1]`,
///   `B = R0⊕R2 = [1,2]` prov `[0,2]`. Rows become `[R3, R4, A, B]`. Column 1 is still in R3, A,
///   B: weight 3, and A and B share original row 0.
/// - Eliminate column 1. Incident rows R3 (4 columns), A (1), B (2): pivot is A. Replacements in
///   ascending order of the others: `A⊕R3 = [2,3,4]` prov `{0,1}⊕{3} = [0,1,3]`, and
///   `A⊕B = [2]` with provenance `{0,1}⊕{0,2} = [1,2]` — original row 0 cancels. A union would
///   claim `[0,1,2]`; indeed `R1⊕R2 = [0,1]⊕[0,1,2] = [2]` involves only rows 1 and 2.
///
/// Final rows `[R4, A⊕R3, A⊕B]`. By hand the sole dependency of the start is `R0+R1+R3+R4`
/// (`[0]⊕[0,1]⊕[1,2,3,4]⊕[2,3,4] = []`), and in the reduced matrix `R4 ⊕ (A⊕R3) = 0` lifts to
/// `{4} ⊕ {0,1,3} = {0,1,3,4}`.
#[test]
fn kat_overlapping_provenance_cancels() {
    let rows = [0b000001u64, 0b000011, 0b000111, 0b011110, 0b011100];
    assert_eq!(kernel_by_enumeration(&rows), vec![0b00000, 0b11011]);

    let mut m = matrix_from_masks(&rows, 6);
    m.obstruction_col_start = 2;
    m.obstruction_count = 4;

    assert_eq!(eliminate_column(&mut m, 0), Ok(()));
    assert_weights_true(&m);
    assert_eq!(m.col_weights[0], 0);
    assert_eq!(
        m.rows[2..],
        [
            MatrixRow { cols: vec![1], provenance: vec![0, 1] },
            MatrixRow { cols: vec![1, 2], provenance: vec![0, 2] },
        ]
    );
    assert_eq!(m.col_weights[1], 3);

    assert_eq!(eliminate_column(&mut m, 1), Ok(()));
    assert_weights_true(&m);
    assert_eq!(
        m.rows,
        vec![
            MatrixRow { cols: vec![2, 3, 4], provenance: vec![4] },
            MatrixRow { cols: vec![2, 3, 4], provenance: vec![0, 1, 3] },
            MatrixRow { cols: vec![2], provenance: vec![1, 2] },
        ]
    );
    assert_lossless(&rows, &m);
    assert_eq!(m.lift_dependency(&[0, 1]), Some(vec![0, 1, 3, 4]));

    // The route through `merge_columns` takes the same two steps and ends at the same matrix.
    let mut start = matrix_from_masks(&rows, 6);
    start.obstruction_col_start = 2;
    start.obstruction_count = 4;
    assert_eq!(merge_columns(start).rows, m.rows);
}

/// Weight-2 elimination of rows whose provenance overlaps, and the row/lift primitives.
///
/// Rows `P = [0,3]` prov `[0,1]` and `Q = [0,4]` prov `[0,2]` are the only rows with column 0
/// (weight two). Their XOR is `[3,4]` with provenance `{0,1}⊕{0,2} = [1,2]`: column 0 and
/// original row 0 both cancel. (This is the shape produced by weight-three elimination.)
#[test]
fn kat_weight_two_elimination_with_overlapping_provenance() {
    let p = MatrixRow { cols: vec![0, 3], provenance: vec![0, 1] };
    let q = MatrixRow { cols: vec![0, 4], provenance: vec![0, 2] };
    assert_eq!(
        p.xor_merge(&q),
        MatrixRow { cols: vec![3, 4], provenance: vec![1, 2] },
        "symmetric difference on both vectors"
    );
    assert_eq!(p.xor_merge(&p), MatrixRow { cols: vec![], provenance: vec![] });

    let mut m = SparseMatrix {
        rows: vec![p, q],
        num_cols: 5,
        obstruction_col_start: 1,
        obstruction_count: 4,
        col_weights: vec![2, 0, 0, 1, 1],
    };
    assert_eq!(eliminate_column(&mut m, 0), Ok(()));
    assert_eq!(m.rows, vec![MatrixRow { cols: vec![3, 4], provenance: vec![1, 2] }]);
    assert_eq!(m.col_weights, vec![0, 0, 0, 1, 1]);
    assert_weights_true(&m);
}

/// `lift_dependency` is the GF(2) sum of provenance; duplicates cancel; bad indices are `None`.
#[test]
fn kat_lift_dependency_is_symmetric_difference() {
    let m = SparseMatrix {
        rows: vec![
            MatrixRow { cols: vec![], provenance: vec![0, 1, 3] },
            MatrixRow { cols: vec![], provenance: vec![1, 2] },
            MatrixRow { cols: vec![], provenance: vec![3] },
        ],
        num_cols: 1,
        obstruction_col_start: 1,
        obstruction_count: 0,
        col_weights: vec![0],
    };
    // {0,1,3} ⊕ {1,2} = {0,2,3};  ... ⊕ {3} = {0,2}.
    assert_eq!(m.lift_dependency(&[0, 1]), Some(vec![0, 2, 3]));
    assert_eq!(m.lift_dependency(&[0, 1, 2]), Some(vec![0, 2]));
    assert_eq!(m.lift_dependency(&[1, 1]), Some(vec![]), "a repeated row cancels");
    assert_eq!(m.lift_dependency(&[]), Some(vec![]));
    assert_eq!(m.lift_dependency(&[0, 3]), None, "row 3 does not exist");
}

/// `eliminate_column` refuses what it must refuse and leaves the matrix untouched when it does.
///
/// Matrix: 4 columns, column 3 the obstruction. Rows `[0]`, `[0,1]`, `[0,1,2]`, `[0,1,2,3]` give
/// column weights `0:4, 1:3, 2:2, 3:1` by counting. Column 0 (weight 4) is unsupported, column 3
/// is an obstruction column, and column 4 does not exist. A second matrix with rows `[2]`, `[3]`
/// has weights `0:0, 1:0, 2:1, 3:1`, giving the weight-0 and weight-1 refusals. A weight table
/// that disagrees with the rows is reported as inconsistent.
#[test]
fn kat_eliminate_column_refusals() {
    let rows = [0b0001u64, 0b0011, 0b0111, 0b1111];
    let m = matrix_from_masks(&rows, 4);
    assert_eq!(m.col_weights, vec![4, 3, 2, 1]);

    let refuse = |m: &SparseMatrix, col: usize| {
        let mut copy = m.clone();
        let err = eliminate_column(&mut copy, col).unwrap_err();
        assert_eq!(copy.rows, m.rows, "refusal leaves rows unchanged");
        assert_eq!(copy.col_weights, m.col_weights, "refusal leaves weights unchanged");
        err
    };
    assert_eq!(refuse(&m, 0), EliminationRefusal::UnsupportedWeight(4));
    assert_eq!(refuse(&m, 3), EliminationRefusal::Obstruction);
    assert_eq!(refuse(&m, 4), EliminationRefusal::OutOfRange);
    assert_eq!(refuse(&m, usize::MAX), EliminationRefusal::OutOfRange);

    let zero_and_one = matrix_from_masks(&[0b0100, 0b1000], 4);
    assert_eq!(zero_and_one.col_weights, vec![0, 0, 1, 1]);
    assert_eq!(refuse(&zero_and_one, 0), EliminationRefusal::UnsupportedWeight(0));
    assert_eq!(refuse(&zero_and_one, 2), EliminationRefusal::UnsupportedWeight(1));

    // Weight table says 2 for column 1, but three rows contain it.
    let mut tampered = m.clone();
    tampered.col_weights[1] = 2;
    assert_eq!(refuse(&tampered, 1), EliminationRefusal::InconsistentWeights);
    // And the table says 3 where two rows contain it.
    let mut tampered = m.clone();
    tampered.col_weights[2] = 3;
    assert_eq!(refuse(&tampered, 2), EliminationRefusal::InconsistentWeights);

    for refusal in [
        EliminationRefusal::OutOfRange,
        EliminationRefusal::Obstruction,
        EliminationRefusal::UnsupportedWeight(7),
        EliminationRefusal::InconsistentWeights,
    ] {
        let message = refusal.to_string();
        assert!(!message.is_empty());
        let as_error: &dyn std::error::Error = &refusal;
        assert_eq!(as_error.to_string(), message);
    }
    assert!(EliminationRefusal::UnsupportedWeight(7).to_string().contains('7'));
}

// ─── Checked-matrix layout ───────────────────────────────────────────────────

/// A checked matrix with omitted (always-zero) exceptional columns, an empty row, and two
/// character columns; elimination must be lossless and must never touch the zero columns.
///
/// Fixture: `f = 3x² + x − 1`, `m = 3`, `N = 29` (`f(3) = 27 + 3 − 1 = 29`), `c = 3`,
/// `F = t² + t − 3`, `disc(F) = 13`. `FactorBase::new(f, 7, 17)`: rational primes `[2,3,5,7]`
/// (columns 0–3); ideals `(3,1)#0 (13,2)#1 (17,4)#2 (17,7)#3` (columns 4–7), of which `(3,1)`
/// and `(13,2)` are exceptional (3 | c, 13 | disc) and so columns 4, 5 are omitted; the sign is
/// column 8; with `K = 2` characters the columns are 9, 10 and `num_cols = 11`.
///
/// Relations (hand-derived norms `a − 3b` and `f_hom(a,b) = 3a² + ab − b²`):
///
/// | (a,b)    | a−3b | f_hom | exponent data                 | checked row  |
/// |----------|------|-------|-------------------------------|--------------|
/// | (4,1)    | 1    | 51    | alg (3,1)^1, (17,4)^1         | `[6]`        |
/// | (7,1)    | 4    | 153   | rat 2^2; alg (3,1)^2, (17,7)  | `[7]`        |
/// | (2,1)    | −1   | 13    | alg (13,2)^1                  | `[8]`        |
/// | (−3,2)   | −9   | 17    | rat 3^2; alg (17,7)           | `[7,8]`      |
/// | (−9,2)   | −15  | 221   | rat 3,5; alg (13,2),(17,4)    | `[1,2,6,8]`  |
/// | (5,2)    | −1   | 81    | alg (3,1)^4                   | `[8]`        |
/// | (−2,−1)  | 1    | 13    | alg (13,2)^1                  | `[]` (empty) |
///
/// The last relation has `a − 3b = −2 + 3 = 1` (no rational prime, positive) and its only odd
/// exponent is at the exceptional ideal `(13,2)` (`−2 ≡ −1·2 mod 13`), so its checked row is
/// empty: a dependency by itself. Provenance is the relation index.
///
/// Elimination by hand. Column 6 has weight two (rows 0 and 4): `[6] ⊕ [1,2,6,8] = [1,2,8]` with
/// provenance `[0,4]`. Column 7 has weight two (rows 1 and 3): `[7] ⊕ [7,8] = [8]` with
/// provenance `[1,3]`. Columns 1, 2 have weight one, column 8 is the obstruction. Result rows,
/// in order: `[8]{2}`, `[8]{5}`, `[]{6}`, `[1,2,8]{0,4}`, `[8]{1,3}`. The starting kernel has
/// dimension `7 − rank`, where the rows span `e6, e7, e8, e1+e2`: rank 4, so 8 dependencies.
#[test]
fn kat_checked_matrix_elimination_keeps_omitted_columns_zero() {
    let f = IntPoly::from_coeffs(vec![bi(-1), bi(1), bi(3)]);
    let ctx = SqrtContext::new(&f, &bi(3), &bi(29)).expect("valid context");
    let fb = FactorBase::new(&f, 7, 17);
    let checker = RelationChecker::new(&ctx, &fb).expect("valid factor base");
    let relations = vec![
        make_relation(4, 1, vec![], vec![(0, 1), (2, 1)], false),
        make_relation(7, 1, vec![(0, 2)], vec![(0, 2), (3, 1)], false),
        make_relation(2, 1, vec![], vec![(1, 1)], true),
        make_relation(-3, 2, vec![(1, 2)], vec![(3, 1)], true),
        make_relation(-9, 2, vec![(1, 1), (2, 1)], vec![(1, 1), (2, 1)], true),
        make_relation(5, 2, vec![], vec![(0, 4)], true),
        make_relation(-2, -1, vec![], vec![(1, 1)], false),
    ];
    let checked = CheckedRelations::new(&checker, &relations).expect("relations validate");
    let (matrix, layout) = build_checked_matrix(&checked, 2).expect("column count fits");

    let expected_rows: Vec<Vec<usize>> =
        vec![vec![6], vec![7], vec![8], vec![7, 8], vec![1, 2, 6, 8], vec![8], vec![]];
    let built: Vec<Vec<usize>> = matrix.rows.iter().map(|r| r.cols.clone()).collect();
    assert_eq!(built, expected_rows);
    assert_eq!((matrix.num_cols, layout.omitted_columns()), (11, vec![4, 5]));
    assert_weights_true(&matrix);

    let original: Vec<u64> = matrix.rows.iter().map(row_mask).collect();
    let reduced = merge_columns(matrix.clone());

    assert_eq!(
        reduced.rows,
        vec![
            MatrixRow { cols: vec![8], provenance: vec![2] },
            MatrixRow { cols: vec![8], provenance: vec![5] },
            MatrixRow { cols: vec![], provenance: vec![6] },
            MatrixRow { cols: vec![1, 2, 8], provenance: vec![0, 4] },
            MatrixRow { cols: vec![8], provenance: vec![1, 3] },
        ]
    );
    assert_eq!(reduced.col_weights, vec![0, 1, 1, 0, 0, 0, 0, 0, 4, 0, 0]);
    assert_weights_true(&reduced);
    assert_lossless(&original, &reduced);
    assert_eq!(kernel_by_enumeration(&original).len(), 8);

    // Layout is untouched: same shape, omitted and character columns still zero.
    assert!(layout.shape_matches(&reduced));
    for col in layout.omitted_columns().into_iter().chain([9, 10]) {
        assert_eq!(reduced.col_weights[col], 0, "column {col} must remain zero");
    }
    // The empty row survives elimination as its own dependency, lifting to relation 6.
    assert_eq!(reduced.lift_dependency(&[2]), Some(vec![6]));
    // The stepwise driver agrees on losslessness for this matrix too.
    let mut stepped = matrix;
    eliminate_stepwise(&mut stepped);
    assert_lossless(&original, &stepped);
}

// ─── Pruning is a separately reported loss ───────────────────────────────────

/// Pruning drops a row and so loses dependencies; the report names the loss exactly.
///
/// Matrix: 2 columns (column 0 ordinary, column 1 the obstruction), 22 rows: `R0 = [0,1]`,
/// `R1 = [1]`, `R2..R21 = [0]`; provenance `[i]`. Excess `= 22 − (2 − 1) = 21 > EXCESS_FLOOR = 20`.
/// The unique heaviest row is R0 (2 columns), so one row is removed and the excess becomes 20.
///
/// Hand count of the loss. Column 1 forces `x0 = x1`; column 0 forces `x0 + x2 + … + x21 = 0`.
/// The rows span `e0` and `e1`, rank 2, so the starting kernel has `2^(22−2) = 2^20` vectors.
/// After pruning, 21 rows of rank 2 leave `2^19`. The `2^19` vectors with `x0 = x1 = 1` are
/// exactly the candidates lost; the survivors are those with `x0 = 0`, i.e. the restricted
/// kernel on the surviving rows.
#[test]
fn kat_pruning_reports_the_candidates_it_loses() {
    use gnfs::EXCESS_FLOOR;

    let mut masks = vec![0b11u64, 0b10];
    masks.extend(std::iter::repeat(0b01u64).take(20));
    let mut start = matrix_from_masks(&masks, 2);
    start.obstruction_col_start = 1;
    start.obstruction_count = 1;
    assert_eq!(start.excess(), 21);

    let (pruned, report) = prune_cliques_reported(start.clone());
    assert_eq!(
        report,
        PruneReport {
            removed: vec![MatrixRow { cols: vec![0, 1], provenance: vec![0] }],
            excess_before: 21,
            excess_after: EXCESS_FLOOR as isize,
        }
    );
    assert!(!report.is_lossless());
    assert_eq!(pruned.rows.len(), 21);
    assert_weights_true(&pruned);
    // The unreported entry point returns the same matrix.
    let plain = prune_cliques(start);
    assert_eq!((plain.rows.clone(), plain.col_weights.clone()), (pruned.rows.clone(), pruned.col_weights.clone()));

    let full = kernel_by_enumeration(&masks);
    let surviving: Vec<u64> = pruned.rows.iter().map(row_mask).collect();
    let kept = kernel_by_enumeration(&surviving);
    assert_eq!((full.len(), kept.len()), (1 << 20, 1 << 19));

    // Pruned rows keep singleton provenance, so lifting is just re-indexing.
    let restricted: Vec<u32> = full.iter().copied().filter(|x| x & 1 == 0).collect();
    let mut lifted: Vec<u32> = kept.iter().map(|&y| lift_mask(&pruned, y)).collect();
    lifted.sort_unstable();
    assert_eq!(lifted, restricted, "pruning keeps exactly the kernel supported on survivors");
    let lost = full.iter().filter(|x| *x & 1 == 1).count();
    assert_eq!(lost, 1 << 19, "every vector using the removed row is lost");
}

/// With excess at or below the floor, pruning removes nothing and reports no loss.
///
/// Matrix: 3 rows, 3 columns with one obstruction column, so excess `= 3 − (3 − 1) = 1 ≤ 20`.
#[test]
fn kat_pruning_below_floor_is_lossless() {
    let start = matrix_from_masks(&[0b001, 0b011, 0b101], 3);
    let (pruned, report) = prune_cliques_reported(start.clone());
    assert!(report.is_lossless());
    assert_eq!((report.excess_before, report.excess_after), (1, 1));
    assert_eq!(pruned.rows, start.rows);
}

/// Elimination after pruning keeps exactly the kernel restricted to the surviving rows.
///
/// Starting rows (4 columns, column 3 the obstruction): `R0=[0,1] R1=[0,2] R2=[0,1,2] R3=[1,2]
/// R4=[0]`, provenance `0..4`. Row 4 is dropped by hand, standing in for a pruned row.
///
/// Hand derivation. `R0⊕R1 = [1,2] = R3`, so `{0,1,3}` is a dependency. `R2⊕R3 = [0] = R4`, so
/// `{2,3,4}` is one too, and their sum `{0,1,2,4}`. The rows span `e0+e1, e0+e2, e0+e1+e2, e0`
/// (rank 3 of 5), so the kernel has `2^2 = 4` vectors: `∅, {0,1,3}, {2,3,4}, {0,1,2,4}`. The
/// kernel restricted to surviving rows (`x4 = 0`) is `{∅, {0,1,3}}`: the two vectors using row 4
/// are gone with it. The four surviving rows have weight three in every column 0–2, so all three
/// are eliminated, and the lifted reduced kernel must be exactly the restricted kernel, never
/// the two lost vectors and never a spurious one.
#[test]
fn kat_elimination_after_dropping_rows_keeps_the_restricted_kernel() {
    let rows = [0b0011u64, 0b0101, 0b0111, 0b0110, 0b0001];
    assert_eq!(
        kernel_by_enumeration(&rows),
        vec![0b00000, 0b01011, 0b10111, 0b11100],
        "full kernel: ∅, {{0,1,3}}, {{0,1,2,4}}, {{2,3,4}}"
    );
    let restricted: Vec<u32> =
        kernel_by_enumeration(&rows).into_iter().filter(|x| x >> 4 & 1 == 0).collect();
    assert_eq!(restricted, vec![0b00000, 0b01011]);

    let mut m = matrix_from_masks(&rows, 4);
    m.remove_row(4);
    assert_eq!(m.col_weights, vec![3, 3, 3, 0]);

    let reduced = merge_columns(m);
    assert_weights_true(&reduced);
    let reduced_masks: Vec<u64> = reduced.rows.iter().map(row_mask).collect();
    let mut lifted: Vec<u32> = kernel_by_enumeration(&reduced_masks)
        .iter()
        .map(|&y| lift_mask(&reduced, y))
        .collect();
    lifted.sort_unstable();
    assert_eq!(lifted, restricted, "lifted kernel equals the kernel restricted to survivors");
    assert_lossless(&rows[..4], &reduced);
}
