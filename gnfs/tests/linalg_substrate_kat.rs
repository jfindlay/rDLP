//! Known-answer tests (KATs) for the GF(2) linear-algebra substrate.
//!
//! Four KATs:
//!
//! - **KAT 1 — Operator correctness**: for a small hand-built `SparseMatrix`, `A·V` and
//!   `Aᵀ·V` match the hand-computed GF(2) products for several block vectors `V`.
//!
//! - **KAT 2 — QC column construction**: widening `obstruction_count` to `1 + num_qc`
//!   yields a matrix of width `matrix_width()` with the sign column at
//!   `obstruction_col_start` and `num_qc` QC columns following; each row's QC parity
//!   matches the hand-computed Legendre-symbol parity for a toy relation set.
//!
//! - **KAT 3 — Round-trip with provenance**: a hand-built kernel vector expands through
//!   the filtering provenance map to the expected set of original relation indices (the
//!   linear algebra → square root seam exercised early).
//!
//! - **KAT 4 — Determinism**: the operator products and QC columns are deterministic for
//!   a fixed matrix.
//!
//! # Toy setup
//!
//! KATs 1–4 use `f(x) = x³ − x − 1` (coefficients: [−1, −1, 0, 1] least-significant
//! first) with `B_rat = B_alg = 13`. KATs 2 and 4 exercise the legacy norm-character helpers
//! (`populate_qc_columns`, `select_qc_primes`), whose behaviour is unchanged.
//!
//! # Individual root characters
//!
//! The `kat_root_character_*` tests exercise the characters `(q, s)` that send a relation
//! `(a, b)` to the Legendre symbol of `a − b·s` modulo `q`. They use `f = x² − 3` with `m = 7`,
//! `N = 23` (`f(7) = 46`) and `f = 3x² + x − 1` with `m = 3`, `N = 29` (`f(3) = 29`). Every
//! root list, bit, and selection below was recomputed with an independent brute-force script
//! that shares no code with the library.

use gnfs::{
    build_matrix, ExponentVector, FactorBase, Relation,
    linalg::{BlockVec, BLOCK_WIDTH, KernelVector, MatrixOperator, populate_qc_columns,
             select_qc_primes},
};
use gnfs::filter::{build_checked_matrix, merge_columns, MatrixLayout, MatrixRow, SparseMatrix};
use gnfs::linalg::{
    select_root_characters, CharacterResource, CharacterValue, PairRejection, RejectedPair,
    RootCharacterError, RootCharacterLimits, RootCharacters,
};
use gnfs::sieve::{CheckedRelations, PrimeClass, RelationChecker};
use gnfs::sqrt::{
    reconstruct_root_from_factors, AlgebraicSqrtOutcome, CrtLimits, DirectLimits,
    PrecisionLimits, PrimeSearchLimits, RecombinationLimits, RelationFactors, SqrtContext,
    SqrtLimits,
};
use num_bigint::BigInt;
use shared_numfield::{IntPoly, NonresidueSearchLimits, RootSearchLimits};

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

/// Build a small hand-crafted 3×4 sparse matrix (no factor base needed).
///
/// ```text
///     col: 0 1 2 3
/// row 0:   1 0 1 0
/// row 1:   0 1 0 1
/// row 2:   1 1 0 0
/// ```
///
/// Provenance: row 0 → [0], row 1 → [1], row 2 → [2].
fn hand_matrix_3x4() -> SparseMatrix {
    SparseMatrix {
        rows: vec![
            MatrixRow { cols: vec![0, 2], provenance: vec![0] },
            MatrixRow { cols: vec![1, 3], provenance: vec![1] },
            MatrixRow { cols: vec![0, 1], provenance: vec![2] },
        ],
        num_cols: 4,
        obstruction_col_start: 4,
        obstruction_count: 0,
        col_weights: vec![2, 2, 1, 1],
    }
}

/// Compute the GF(2) matrix-vector product A·v for a dense matrix `a` (rows × cols)
/// and a single dense vector `v` (length cols). Returns a dense vector of length rows.
fn dense_matvec_gf2(a: &[Vec<bool>], v: &[bool]) -> Vec<bool> {
    a.iter()
        .map(|row| {
            row.iter().zip(v.iter()).filter(|&(&r, &vi)| r && vi).count() % 2 == 1
        })
        .collect()
}

// ─── KAT 1: Operator correctness ─────────────────────────────────────────────

/// KAT 1: `A·V` and `Aᵀ·V` match hand-computed GF(2) products.
///
/// Matrix:
/// ```text
///     col: 0 1 2 3
/// row 0:   1 0 1 0
/// row 1:   0 1 0 1
/// row 2:   1 1 0 0
/// ```
///
/// We test several block vectors and verify each column of the output matches the
/// hand-computed GF(2) matrix-vector product.
#[test]
fn kat_1_operator_correctness() {
    let matrix = hand_matrix_3x4();
    let op = MatrixOperator::new(&matrix);

    assert_eq!(op.num_rows(), 3);
    assert_eq!(op.num_cols(), 4);

    // Dense representation of the matrix for reference computation.
    let a_dense: Vec<Vec<bool>> = vec![
        vec![true, false, true, false],  // row 0
        vec![false, true, false, true],  // row 1
        vec![true, true, false, false],  // row 2
    ];

    // Dense representation of Aᵀ (4×3).
    let at_dense: Vec<Vec<bool>> = vec![
        vec![true, false, true],   // row 0 of Aᵀ = col 0 of A
        vec![false, true, true],   // row 1 of Aᵀ = col 1 of A
        vec![true, false, false],  // row 2 of Aᵀ = col 2 of A
        vec![false, true, false],  // row 3 of Aᵀ = col 3 of A
    ];

    // Test several input vectors for apply (A·V).
    let test_vecs_apply: Vec<Vec<bool>> = vec![
        // Standard basis vectors.
        vec![true, false, false, false],
        vec![false, true, false, false],
        vec![false, false, true, false],
        vec![false, false, false, true],
        // Mixed vectors.
        vec![true, true, false, false],
        vec![true, false, true, false],
        vec![true, true, true, true],
        vec![false, false, false, false],
    ];

    for (vi, test_vec) in test_vecs_apply.iter().enumerate() {
        // Build a BlockVec with this vector in column 0.
        let mut bv = BlockVec::zeros(4);
        for (r, &bit) in test_vec.iter().enumerate() {
            bv.set(r, 0, bit);
        }

        let result = op.apply(&bv);
        let expected = dense_matvec_gf2(&a_dense, test_vec);

        for r in 0..3 {
            assert_eq!(
                result.get(r, 0),
                expected[r],
                "apply: test vector {vi}, row {r}: expected {}, got {}",
                expected[r],
                result.get(r, 0)
            );
        }
    }

    // Test several input vectors for apply_transpose (Aᵀ·V).
    let test_vecs_transpose: Vec<Vec<bool>> = vec![
        vec![true, false, false],
        vec![false, true, false],
        vec![false, false, true],
        vec![true, true, false],
        vec![true, true, true],
        vec![false, false, false],
    ];

    for (vi, test_vec) in test_vecs_transpose.iter().enumerate() {
        let mut bv = BlockVec::zeros(3);
        for (r, &bit) in test_vec.iter().enumerate() {
            bv.set(r, 0, bit);
        }

        let result = op.apply_transpose(&bv);
        let expected = dense_matvec_gf2(&at_dense, test_vec);

        for r in 0..4 {
            assert_eq!(
                result.get(r, 0),
                expected[r],
                "apply_transpose: test vector {vi}, row {r}: expected {}, got {}",
                expected[r],
                result.get(r, 0)
            );
        }
    }

    // Test with multiple columns in a single BlockVec.
    // Use all 4 standard basis vectors simultaneously (columns 0..3).
    let mut bv_multi = BlockVec::zeros(4);
    for j in 0..4 {
        bv_multi.set(j, j, true); // column j = e_j
    }
    let result_multi = op.apply(&bv_multi);
    // Column j of result should be A·e_j = column j of A.
    for j in 0..4 {
        for r in 0..3 {
            assert_eq!(
                result_multi.get(r, j),
                a_dense[r][j],
                "apply multi: col {j}, row {r}"
            );
        }
    }

    // Verify Aᵀ·(A·e_j) = Aᵀ·(col j of A).
    for j in 0..4 {
        let mut bv_col = BlockVec::zeros(3);
        for r in 0..3 {
            bv_col.set(r, 0, a_dense[r][j]);
        }
        let result_at = op.apply_transpose(&bv_col);
        let expected_at = dense_matvec_gf2(&at_dense, &a_dense.iter().map(|row| row[j]).collect::<Vec<_>>());
        for r in 0..4 {
            assert_eq!(
                result_at.get(r, 0),
                expected_at[r],
                "apply_transpose of col {j}: row {r}"
            );
        }
    }
}

// ─── KAT 2: QC column construction ───────────────────────────────────────────

/// KAT 2: QC column construction with hand-computed Legendre-symbol parities.
///
/// Setup: `f(x) = x³ − x − 1`, `B_rat = B_alg = 13`.
/// We widen `obstruction_count` to `1 + num_qc` and call `populate_qc_columns`.
/// For each row, we verify the QC parity matches the hand-computed Legendre symbol.
///
/// # Hand computation
///
/// For a relation with `(a, b)`, the algebraic norm is `N_alg(a, b) = b³·f(a/b)`.
/// For `b = 1`: `N_alg(a, 1) = f(a) = a³ − a − 1`.
///
/// We use `select_qc_primes` to find the first QC prime `q > 13` that splits completely
/// in K, then verify the QC parity for each relation by computing the Legendre symbol
/// `(N_alg(a,b) / q)` directly.
#[test]
fn kat_2_qc_column_construction() {
    use gnfs::polyselect::PolyPair;

    let f = f_cubic();
    let mut fb = FactorBase::new(&f, 13, 13);

    // Select 2 QC primes.
    let num_qc = 2usize;
    let qc_primes = select_qc_primes(&f, 13, num_qc);
    assert_eq!(qc_primes.len(), num_qc, "should have {num_qc} QC primes");
    for &q in &qc_primes {
        assert!(q > 13, "QC prime {q} should be > b_alg = 13");
    }

    // Widen obstruction_count to 1 + num_qc.
    fb.obstruction_count = 1 + num_qc;

    // Verify matrix_width() reflects the widened obstruction_count.
    let expected_width = fb.rational_size() + fb.algebraic_size() + 1 + num_qc;
    assert_eq!(fb.matrix_width(), expected_width, "matrix_width should include QC columns");

    // Build a toy PolyPair for norm computation.
    // N = 503 (prime), m = 7: f(7) = 343 - 7 - 1 = 335; 335 mod 503 ≠ 0.
    // Use m = 8: f(8) = 512 - 8 - 1 = 503 ≡ 0 (mod 503). ✓
    let n = bi(503);
    let m = bi(8);
    let g = IntPoly::from_coeffs(vec![-m.clone(), bi(1)]);
    let poly = PolyPair::new(f.clone(), g, m, n);

    // Build relations with known (a, b) pairs.
    // We use b=1 so N_alg(a, 1) = f(a) = a³ - a - 1.
    // Relations: (a=2, b=1), (a=3, b=1), (a=4, b=1).
    // N_alg(2, 1) = 8 - 2 - 1 = 5
    // N_alg(3, 1) = 27 - 3 - 1 = 23
    // N_alg(4, 1) = 64 - 4 - 1 = 59
    let idx_p5 = fb.rational_index(5).unwrap_or(2);
    let idx_p23 = fb.rational_index(23).unwrap_or(0); // 23 > 13, not in base
    let idx_p59 = fb.rational_index(59).unwrap_or(0); // 59 > 13, not in base

    // Use simple relations with rational primes in the factor base.
    // The exact exponent vectors don't matter for QC testing; we just need valid (a, b).
    let relations = vec![
        make_relation(2, 1, vec![(idx_p5, 1)], vec![], false),
        make_relation(3, 1, vec![(idx_p5, 1)], vec![], false),
        make_relation(4, 1, vec![(idx_p5, 1)], vec![], false),
    ];
    let _ = (idx_p23, idx_p59); // suppress unused warnings

    // Build the matrix with the widened obstruction_count.
    let mut matrix = build_matrix(&relations, &fb);

    // Verify the matrix width matches fb.matrix_width().
    assert_eq!(matrix.num_cols, fb.matrix_width(), "matrix num_cols should equal fb.matrix_width()");

    // Verify obstruction_col_start is correct.
    let expected_obs_start = fb.rational_size() + fb.algebraic_size();
    assert_eq!(
        matrix.obstruction_col_start, expected_obs_start,
        "obstruction_col_start should be rational_size + algebraic_size"
    );

    // Before populate_qc_columns: QC columns should all be zero.
    for (row_idx, row) in matrix.rows.iter().enumerate() {
        for k in 0..num_qc {
            let qc_col = matrix.obstruction_col_start + 1 + k;
            assert!(
                !row.cols.contains(&qc_col),
                "row {row_idx}: QC column {qc_col} should be zero before populate_qc_columns"
            );
        }
    }

    // Populate QC columns.
    populate_qc_columns(&mut matrix, &relations, &fb, &poly, &qc_primes);

    // Verify QC parities match hand-computed Legendre symbols.
    // For each relation (a, b), compute N_alg(a, b) mod q and check the Legendre symbol.
    let norms = vec![5u64, 23u64, 59u64]; // N_alg(2,1), N_alg(3,1), N_alg(4,1)

    for (row_idx, row) in matrix.rows.iter().enumerate() {
        let norm = norms[row_idx];
        for (k, &q) in qc_primes.iter().enumerate() {
            let qc_col = matrix.obstruction_col_start + 1 + k;
            let norm_mod_q = norm % q;
            // Legendre symbol via Euler's criterion.
            let exp = (q - 1) / 2;
            let symbol = pow_mod_test(norm_mod_q, exp, q);
            let expected_parity = symbol == q - 1;
            let actual_parity = row.cols.contains(&qc_col);
            assert_eq!(
                actual_parity, expected_parity,
                "row {row_idx}, QC prime {q}: expected parity {expected_parity}, got {actual_parity} \
                 (norm={norm}, norm mod q={norm_mod_q}, symbol={symbol})"
            );
        }
    }

    // Verify col_weights are consistent with row contents.
    for col in 0..matrix.num_cols {
        let actual_weight =
            matrix.rows.iter().filter(|r| r.cols.binary_search(&col).is_ok()).count();
        assert_eq!(
            matrix.col_weights[col] as usize,
            actual_weight,
            "col_weights[{col}] should be consistent with row contents after populate_qc_columns"
        );
    }
}

// ─── KAT 3: Round-trip with provenance ───────────────────────────────────────

/// KAT 3: A hand-built kernel vector expands through the provenance map to the expected
/// set of original relation indices (the linear algebra → square root seam).
///
/// Setup: a matrix with merged rows (provenance sets of size > 1), and a kernel vector
/// that selects a subset of rows. We verify that `expand_provenance` returns the correct
/// symmetric difference of the provenance sets.
///
/// # Provenance structure
///
/// - Row 0: provenance = [0, 1] (merged from original relations 0 and 1)
/// - Row 1: provenance = [2]    (original relation 2)
/// - Row 2: provenance = [1, 3] (merged from original relations 1 and 3)
///
/// Kernel vector: rows {0, 2}.
/// Expected expansion: sym_diff([0,1], [1,3]) = [0, 3] (1 cancels).
///
/// Kernel vector: rows {0, 1, 2}.
/// Expected expansion: sym_diff([0,1], [2], [1,3]) = [0, 2, 3] (1 cancels from 0 and 2).
#[test]
fn kat_3_round_trip_provenance() {
    // Build a matrix with hand-crafted provenance.
    // The column sets are chosen so that rows {0, 1, 2} form a valid nullspace vector:
    // row 0 cols = {0, 1}, row 1 cols = {1, 2}, row 2 cols = {0, 2}.
    // XOR: {0,1} XOR {1,2} XOR {0,2} = {} ✓
    let matrix = SparseMatrix {
        rows: vec![
            MatrixRow { cols: vec![0, 1], provenance: vec![0, 1] },
            MatrixRow { cols: vec![1, 2], provenance: vec![2] },
            MatrixRow { cols: vec![0, 2], provenance: vec![1, 3] },
        ],
        num_cols: 3,
        obstruction_col_start: 3,
        obstruction_count: 0,
        col_weights: vec![2, 2, 2],
    };

    // KAT 3a: kernel vector {0, 2}.
    let kv_02 = KernelVector::new(vec![0, 2]);
    // Not a nullspace vector (rows 0 and 2 don't cancel), but we test provenance expansion.
    let expanded_02 = kv_02.expand_provenance(&matrix);
    // sym_diff([0,1], [1,3]) = [0, 3].
    assert_eq!(expanded_02, vec![0, 3], "KAT 3a: expansion of rows {{0,2}} should be [0, 3]");

    // KAT 3b: kernel vector {0, 1, 2} — a valid nullspace vector.
    let kv_012 = KernelVector::new(vec![0, 1, 2]);
    assert!(kv_012.verify(&matrix), "KAT 3b: rows {{0,1,2}} should form a valid nullspace vector");
    let expanded_012 = kv_012.expand_provenance(&matrix);
    // sym_diff([0,1], [2], [1,3]):
    // Step 1: sym_diff([0,1], [2]) = [0, 1, 2].
    // Step 2: sym_diff([0,1,2], [1,3]) = [0, 2, 3] (1 cancels).
    assert_eq!(
        expanded_012, vec![0, 2, 3],
        "KAT 3b: expansion of rows {{0,1,2}} should be [0, 2, 3]"
    );

    // KAT 3c: single-row kernel vector {1}.
    let kv_1 = KernelVector::new(vec![1]);
    let expanded_1 = kv_1.expand_provenance(&matrix);
    assert_eq!(expanded_1, vec![2], "KAT 3c: expansion of row {{1}} should be [2]");

    // KAT 3d: empty kernel vector.
    let kv_empty = KernelVector::new(vec![]);
    let expanded_empty = kv_empty.expand_provenance(&matrix);
    assert_eq!(expanded_empty, vec![], "KAT 3d: expansion of empty kernel vector should be []");

    // KAT 3e: from_mask constructor.
    let mask = vec![true, true, true];
    let kv_from_mask = KernelVector::from_mask(&mask);
    assert_eq!(kv_from_mask.row_indices, vec![0, 1, 2], "KAT 3e: from_mask should produce [0,1,2]");
    assert!(kv_from_mask.verify(&matrix), "KAT 3e: from_mask kernel vector should be valid");

    // KAT 3f: verify that a non-nullspace vector returns false.
    let kv_invalid = KernelVector::new(vec![0, 1]);
    assert!(!kv_invalid.verify(&matrix), "KAT 3f: rows {{0,1}} should not be a nullspace vector");
}

// ─── KAT 4: Determinism ──────────────────────────────────────────────────────

/// KAT 4: Operator products and QC columns are deterministic for a fixed matrix.
///
/// Calls the operator and QC functions twice on the same input and asserts identical
/// results. This verifies there is no hidden state or randomness.
#[test]
fn kat_4_determinism() {
    use gnfs::polyselect::PolyPair;

    let f = f_cubic();
    let mut fb = FactorBase::new(&f, 13, 13);

    let num_qc = 2usize;
    let qc_primes = select_qc_primes(&f, 13, num_qc);
    fb.obstruction_count = 1 + num_qc;

    let n = bi(503);
    let m = bi(8);
    let g = IntPoly::from_coeffs(vec![-m.clone(), bi(1)]);
    let poly = PolyPair::new(f.clone(), g, m, n);

    let idx_p5 = fb.rational_index(5).unwrap_or(2);
    let relations = vec![
        make_relation(2, 1, vec![(idx_p5, 1)], vec![], false),
        make_relation(3, 1, vec![(idx_p5, 1)], vec![], false),
    ];

    // First run: build matrix, populate QC columns, apply operator.
    let mut matrix1 = build_matrix(&relations, &fb);
    populate_qc_columns(&mut matrix1, &relations, &fb, &poly, &qc_primes);
    let op1 = MatrixOperator::new(&matrix1);

    // Build a test block vector.
    let mut bv1 = BlockVec::zeros(op1.num_cols());
    for j in 0..BLOCK_WIDTH.min(op1.num_cols()) {
        bv1.set(j % op1.num_cols(), j % BLOCK_WIDTH, true);
    }
    let result1_apply = op1.apply(&bv1);
    let mut bv1t = BlockVec::zeros(op1.num_rows());
    for j in 0..BLOCK_WIDTH.min(op1.num_rows()) {
        bv1t.set(j % op1.num_rows(), j % BLOCK_WIDTH, true);
    }
    let result1_transpose = op1.apply_transpose(&bv1t);

    // Second run: identical inputs.
    let mut matrix2 = build_matrix(&relations, &fb);
    populate_qc_columns(&mut matrix2, &relations, &fb, &poly, &qc_primes);
    let op2 = MatrixOperator::new(&matrix2);

    let mut bv2 = BlockVec::zeros(op2.num_cols());
    for j in 0..BLOCK_WIDTH.min(op2.num_cols()) {
        bv2.set(j % op2.num_cols(), j % BLOCK_WIDTH, true);
    }
    let result2_apply = op2.apply(&bv2);
    let mut bv2t = BlockVec::zeros(op2.num_rows());
    for j in 0..BLOCK_WIDTH.min(op2.num_rows()) {
        bv2t.set(j % op2.num_rows(), j % BLOCK_WIDTH, true);
    }
    let result2_transpose = op2.apply_transpose(&bv2t);

    // Assert identical results.
    assert_eq!(
        result1_apply, result2_apply,
        "KAT 4: apply results must be identical across runs"
    );
    assert_eq!(
        result1_transpose, result2_transpose,
        "KAT 4: apply_transpose results must be identical across runs"
    );

    // Assert identical matrix contents (QC columns are deterministic).
    assert_eq!(matrix1.rows.len(), matrix2.rows.len(), "KAT 4: row counts must match");
    for i in 0..matrix1.rows.len() {
        assert_eq!(
            matrix1.rows[i].cols, matrix2.rows[i].cols,
            "KAT 4: row {i} cols must be identical across runs"
        );
    }
    assert_eq!(
        matrix1.col_weights, matrix2.col_weights,
        "KAT 4: col_weights must be identical across runs"
    );

    // Assert QC primes are deterministic.
    let qc_primes2 = select_qc_primes(&f, 13, num_qc);
    assert_eq!(qc_primes, qc_primes2, "KAT 4: select_qc_primes must be deterministic");
}

// ─── Individual root characters ──────────────────────────────────────────────

/// `x² − 3`, coefficients least-significant first.
const F_UNIT: [i64; 3] = [-3, 0, 1];
/// `3x² + x − 1`.
const F_NONMONIC: [i64; 3] = [-1, 1, 3];

/// Valuation of `|n|` at `p`.
fn valuation(n: i64, p: u64) -> u32 {
    let (mut n, p, mut e) = (n.abs(), p as i64, 0);
    while n != 0 && n % p == 0 {
        n /= p;
        e += 1;
    }
    e
}

/// `Σ f_i a^i b^(d−i)`.
fn homogeneous(f: &[i64], a: i64, b: i64) -> i64 {
    let d = f.len() - 1;
    f.iter().enumerate().map(|(i, &c)| c * a.pow(i as u32) * b.pow((d - i) as u32)).sum()
}

/// A relation whose exponents are the honest valuations: the rational norm `a − b·m` over the
/// rational primes, and `f_hom(a, b)` attributed to the ideal `(p, r)` with `a ≡ b·r (mod p)`.
/// The checker re-derives both, so a wrong table here is rejected rather than trusted.
fn relation_for(f: &[i64], m: i64, fb: &FactorBase, a: i64, b: i64) -> Relation {
    let rational = a - b * m;
    let rational_entries = fb
        .rational_primes
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| (valuation(rational, p) > 0).then(|| (i, valuation(rational, p))))
        .collect();
    let norm = homogeneous(f, a, b);
    let algebraic_entries = fb
        .algebraic_ideals
        .iter()
        .enumerate()
        .filter(|(_, ideal)| (a - b * ideal.r as i64).rem_euclid(ideal.p as i64) == 0)
        .filter_map(|(i, ideal)| (valuation(norm, ideal.p) > 0).then(|| (i, valuation(norm, ideal.p))))
        .collect();
    make_relation(a, b, rational_entries, algebraic_entries, rational < 0)
}

/// The checker for `f`, `m`, `N` with factor base bounds `(b_rat, b_alg)`, and the checked
/// relations for `pairs`.
fn checked(
    f: &[i64],
    m: i64,
    n: i64,
    bounds: (u64, u64),
    pairs: &[(i64, i64)],
) -> (RelationChecker, CheckedRelations) {
    let poly = IntPoly::from_coeffs(f.iter().map(|&c| bi(c)).collect());
    let fb = FactorBase::new(&poly, bounds.0, bounds.1);
    let ctx = SqrtContext::new(&poly, &bi(m), &bi(n)).expect("fixture context is valid");
    let checker = RelationChecker::new(&ctx, &fb).expect("fixture factor base is valid");
    let relations: Vec<Relation> =
        pairs.iter().map(|&(a, b)| relation_for(f, m, &fb, a, b)).collect();
    let set = CheckedRelations::new(&checker, &relations).expect("fixture relations are honest");
    (checker, set)
}

/// Brute-force Legendre bit of `v` modulo odd prime `q`: `Some(true)` for a nonresidue, `None`
/// for zero. Uses the list of squares, not Euler's criterion.
fn nonresidue(v: i64, q: u64) -> Option<bool> {
    let r = v.rem_euclid(q as i64) as u64;
    (r != 0).then(|| !(1..q).any(|x| x * x % q == r))
}

fn pair(prime: u64, root: u64) -> (u64, u64) {
    (prime, root)
}

fn matrix_with_rows(template: &SparseMatrix, rows: Vec<MatrixRow>) -> SparseMatrix {
    let mut col_weights = vec![0u32; template.num_cols];
    for row in &rows {
        for &c in &row.cols {
            col_weights[c] += 1;
        }
    }
    SparseMatrix {
        rows,
        num_cols: template.num_cols,
        obstruction_col_start: template.obstruction_col_start,
        obstruction_count: template.obstruction_count,
        col_weights,
    }
}

fn character_bits(matrix: &SparseMatrix, layout: &MatrixLayout, row: usize) -> Vec<bool> {
    (0..layout.character_count())
        .map(|slot| {
            let column = layout.character_column(slot).expect("slot is below the count");
            matrix.rows[row].cols.contains(&column)
        })
        .collect()
}

fn roomy_sqrt_limits() -> SqrtLimits {
    let primes = PrimeSearchLimits { first_candidate: 0, max_primes_examined: 100 };
    let nonresidue = NonresidueSearchLimits { max_candidates: 1000 };
    SqrtLimits {
        precision: PrecisionLimits { max_modulus_bits: 4096 },
        direct: DirectLimits { prime_search: primes, nonresidue_search: nonresidue },
        crt: CrtLimits {
            prime_search: primes,
            max_primes: 1,
            root_search: RootSearchLimits { max_evaluations: 100_000 },
            nonresidue_search: nonresidue,
        },
        recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 20 },
    }
}

/// The unit `u = 2 + α` and the square `s = (1 + 2α)² = 13 + 4α` on `x² − 3`, `m = 7`.
///
/// Hand arithmetic. `u` is the pair `(2, −1)`: `a − bm = 2 + 7 = 9 = 3²`, `f_hom = 4 − 3 = 1`,
/// so its matrix row is empty and it satisfies every parity and sign constraint. Its norm is 1 and
/// it is totally positive, yet it is not a square in `ℚ(√3)`: `2 + √3 = (1 + √3)²/2` and `2` is
/// not a square there. `s` is the pair `(13, −4)`: `(1 + 2α)² = 1 + 4α + 4·3 = 13 + 4α` (using
/// `α² = 3`), `a − bm = 13 + 28 = 41`, `f_hom = 169 − 48 = 121 = 11²` at the ideal `(11, 5)`
/// since `13 + 4·5 = 33`. Factor base bounds `(41, 11)` hold both.
fn unit_and_square() -> (RelationChecker, CheckedRelations) {
    checked(&F_UNIT, 7, 23, (41, 11), &[(2, -1), (13, -4)])
}

/// A norm-one nonsquare unit passes the characters at `q = 23` but not the one at `q = 11`.
///
/// Roots of `t² − 3`: mod 23 they are 7 and 16 (49 = 46 + 3, 256 = 253 + 3); mod 11, 5 and 6
/// (25 = 22 + 3, 36 = 33 + 3); mod 13, 4 and 9 (16 = 13 + 3, 81 = 78 + 3). For `u`, `2 + s` is 9
/// and 18 at `q = 23` (squares: 3² and 8² = 64 = 46 + 18), so both bits are 0 and the filter
/// cannot see the obstruction. At `q = 11` the values are 7 and 8, nonresidues (squares mod 11
/// are 1, 3, 4, 5, 9), and at `q = 13` they are 6 and 11, nonresidues (squares mod 13 are 1, 3,
/// 4, 9, 10, 12), so the bits are 1. The square `s` has bit 0 everywhere, as `φ(β)²` must be.
#[test]
fn kat_root_character_unit_passes_some_characters_but_is_not_a_square() {
    let (checker, relations) = unit_and_square();
    let u = &relations.relations()[0];
    // The old norm characters see nothing: the algebraic norm of `u` is 1, a square at every q.
    assert_eq!(*u.algebraic_norm(), bi(1));
    assert_eq!(u.rational_norm(), &bi(9));

    let at_23 = RootCharacters::new(&checker, &[pair(23, 7), pair(23, 16)]).unwrap();
    assert_eq!(at_23.relation_bits(&relations, 0).unwrap(), vec![false, false]);
    assert_eq!(at_23.relation_bits(&relations, 1).unwrap(), vec![false, false]);

    let split = RootCharacters::new(&checker, &[pair(11, 6), pair(13, 4), pair(13, 9)]).unwrap();
    assert_eq!(split.relation_bits(&relations, 0).unwrap(), vec![true, true, true]);
    assert_eq!(split.relation_bits(&relations, 1).unwrap(), vec![false, false, false]);

    // The matrix row of `u` is empty, so `u` alone is a kernel vector of the plain matrix. Filling
    // the q = 23 characters leaves it there; filling the others removes it from the kernel.
    for (characters, survives) in [(&at_23, true), (&split, false)] {
        let (mut matrix, layout) = build_checked_matrix(&relations, characters.len()).unwrap();
        assert!(matrix.rows[0].cols.is_empty());
        characters.fill_columns(&mut matrix, &layout, &relations).unwrap();
        assert_eq!(matrix.rows[0].cols.is_empty(), survives);
        // `s` has the 41 column set and no character bit, so characters never certify it either.
        assert_eq!(matrix.rows[1].cols.len(), 1);
    }

    // Exact reconstruction is authoritative: `u` is certified a nonsquare and `s` has a root.
    let verdict = |pairs: &[(BigInt, BigInt)]| {
        let factors = RelationFactors::new(checker.context(), pairs).unwrap();
        reconstruct_root_from_factors(factors, roomy_sqrt_limits()).unwrap()
    };
    assert!(matches!(verdict(&[(bi(2), bi(-1))]), AlgebraicSqrtOutcome::NotSquare(_)));
    assert!(matches!(verdict(&[(bi(13), bi(-4))]), AlgebraicSqrtOutcome::Verified(_)));
}

/// A prime with several roots yields several characters, and a vanishing pair is reported.
///
/// For `x² − 3` and the relations `u`, `s`: 3 divides `disc(F) = 12`; 5, 7, 17, 19 have 3 as a
/// nonresidue (squares mod 5: 1, 4; mod 7: 1, 2, 4; mod 17: 1, 2, 4, 8, 9, 13, 15, 16; mod 19: 1,
/// 4, 5, 6, 7, 9, 11, 16, 17), so they have no roots; 11 has roots 5, 6 but `13 + 4·5 = 33 ≡ 0`,
/// so `(11, 5)` vanishes on `s` (position 1) and only `(11, 6)` remains; 13 has roots 4 and 9,
/// 23 has 7 and 16.
#[test]
fn kat_root_character_selection_counts_pairs_not_primes() {
    let (_, relations) = unit_and_square();
    let report_of = |count| {
        let selection =
            select_root_characters(&relations, count, RootCharacterLimits::default()).unwrap();
        let pairs: Vec<(u64, u64)> =
            selection.characters().characters().iter().map(|c| (c.prime(), c.root())).collect();
        (pairs, selection.report().clone(), selection.characters().len())
    };

    let (pairs, report, len) = report_of(3);
    assert_eq!(pairs, vec![(11, 6), (13, 4), (13, 9)]);
    assert_eq!(len, 3, "three pairs from only two primes");
    assert_eq!(report.rejected_pairs, vec![RejectedPair { prime: 11, root: 5, relation: 1 }]);
    assert_eq!(report.exceptional_primes, vec![3]);
    assert!(report.truncated_primes.is_empty());
    assert_eq!(report.primes_examined, 5); // 3, 5, 7, 11, 13

    // The search stops mid-prime once the count is reached, and continues to 23 when asked to.
    let (pairs, report, _) = report_of(4);
    assert_eq!(pairs, vec![(11, 6), (13, 4), (13, 9), (23, 7)]);
    assert_eq!(report.primes_examined, 8); // ... 17, 19, 23
    assert_eq!(report_of(5).0.last(), Some(&(23, 16)));

    // Zero characters need no search at all.
    let (pairs, report, _) = report_of(0);
    assert!(pairs.is_empty());
    assert_eq!(report.primes_examined, 0);

    // A start below 3 behaves as 3; a later start skips earlier pairs.
    for first in [0, 1, 2] {
        let limits = RootCharacterLimits {
            prime_search: PrimeSearchLimits { first_candidate: first, max_primes_examined: 100 },
            ..RootCharacterLimits::default()
        };
        let selection = select_root_characters(&relations, 1, limits).unwrap();
        assert_eq!(selection.characters().characters()[0].prime(), 11);
    }
    let later = RootCharacterLimits {
        prime_search: PrimeSearchLimits { first_candidate: 12, max_primes_examined: 100 },
        ..RootCharacterLimits::default()
    };
    let selection = select_root_characters(&relations, 1, later).unwrap();
    assert_eq!(selection.characters().characters()[0].root(), 4);
}

/// Primes dividing `c·disc(F)` and the prime 2 are not valid, and a vanishing relation retains
/// a pair's rejection.
///
/// `f = 3x² + x − 1` has `c = 3` and `F = t² + t − 3`, `disc(F) = 1 + 12 = 13`. Mod 3, `f ≡ x − 1`
/// has the root 1; mod 13, `f` has the double root 2 (`3·4 + 2 − 1 = 13`); mod 17 the roots are 4
/// (48 + 4 − 1 = 51) and 7 (147 + 7 − 1 = 153 = 9·17); mod 5, 7, 11 there are none. With the
/// relations `(2, 1)`, `(5, 2)`, `(−9, 2)`: `(−9, 2)` has `f_hom = 221 = 13·17` and
/// `−9 − 2·4 = −17`, so `(17, 4)` vanishes on it (position 2), while `(17, 7)` gives
/// `2 − 7`, `5 − 14`, `−9 − 14` = `−5`, `−9`, `−23`, none zero mod 17.
#[test]
fn kat_root_character_pairs_exclude_even_and_exceptional_primes() {
    let pairs = [(2, 1), (5, 2), (-9, 2)];
    let (checker, relations) = checked(&F_NONMONIC, 3, 29, (7, 17), &pairs);
    assert_eq!(
        gnfs::sieve::classify_prime(checker.context(), 3),
        PrimeClass::Exceptional,
        "3 | c"
    );
    assert_eq!(
        gnfs::sieve::classify_prime(checker.context(), 13),
        PrimeClass::Exceptional,
        "13 | disc(F)"
    );

    let reason = |list: &[(u64, u64)]| match RootCharacters::new(&checker, list) {
        Err(RootCharacterError::Pair { position, reason, .. }) => (position, reason),
        other => panic!("expected a pair rejection, got {other:?}"),
    };
    assert_eq!(reason(&[(2, 1)]), (0, PairRejection::Even));
    assert_eq!(reason(&[(1, 0)]), (0, PairRejection::NotPrime));
    assert_eq!(reason(&[(9, 0)]), (0, PairRejection::NotPrime));
    assert_eq!(reason(&[(3, 1)]), (0, PairRejection::Exceptional));
    assert_eq!(reason(&[(13, 2)]), (0, PairRejection::Exceptional));
    assert_eq!(reason(&[(17, 17)]), (0, PairRejection::RootOutOfRange));
    assert_eq!(reason(&[(17, 5)]), (0, PairRejection::NotRoot));
    assert_eq!(
        reason(&[(17, 4), (17, 7), (17, 4)]),
        (2, PairRejection::Duplicate { first: 0 })
    );
    let both = RootCharacters::new(&checker, &[(17, 4), (17, 7)]).unwrap();
    assert_eq!(both.len(), 2);

    // The search skips 3 and 13 silently as to pairs but lists them, rejects (17, 4), and uses
    // (17, 7) then the two roots of f mod 23 (3·18² + 18 − 1 = 989 = 43·23, 3·20² + 20 − 1 =
    // 1219 = 53·23).
    let selection = select_root_characters(&relations, 3, RootCharacterLimits::default()).unwrap();
    let found: Vec<(u64, u64)> =
        selection.characters().characters().iter().map(|c| (c.prime(), c.root())).collect();
    assert_eq!(found, vec![(17, 7), (23, 18), (23, 20)]);
    assert_eq!(selection.report().exceptional_primes, vec![3, 13]);
    assert_eq!(
        selection.report().rejected_pairs,
        vec![RejectedPair { prime: 17, root: 4, relation: 2 }]
    );
}

/// A zero is never read as a 0 bit.
///
/// On `x² − 3`, `(11, 5)` is a valid pair (25 = 22 + 3) that vanishes on `s = (13, −4)`:
/// `13 − (−4)·5 = 33 = 3·11`. Evaluation reports `Zero`, with no bit, and filling refuses.
#[test]
fn kat_root_character_zero_value_is_never_a_bit() {
    let (checker, relations) = unit_and_square();
    let characters = RootCharacters::new(&checker, &[pair(11, 6), pair(11, 5)]).unwrap();
    let vanishing = &characters.characters()[1];
    let (a, b) = (relations.relations()[1].a(), relations.relations()[1].b());
    assert_eq!(vanishing.evaluate(a, b), CharacterValue::Zero);
    assert_eq!(CharacterValue::Zero.bit(), None);
    assert_eq!(CharacterValue::Square.bit(), Some(false));
    assert_eq!(CharacterValue::NonSquare.bit(), Some(true));

    let zero = RootCharacterError::ZeroValue { slot: 1, prime: 11, root: 5, relation: 1 };
    assert_eq!(characters.relation_bits(&relations, 1), Err(zero.clone()));
    assert_eq!(characters.bit_table(&relations), Err(zero.clone()));
    // The other relation has honest bits: 2 + 6 = 8 and 2 + 5 = 7 are both nonresidues mod 11.
    assert_eq!(characters.relation_bits(&relations, 0).unwrap(), vec![true, true]);
    assert_eq!(
        characters.relation_bits(&relations, 2),
        Err(RootCharacterError::RelationOutOfRange { index: 2, len: 2 })
    );

    let (mut matrix, layout) = build_checked_matrix(&relations, 2).unwrap();
    let before = matrix.clone();
    assert_eq!(characters.fill_columns(&mut matrix, &layout, &relations), Err(zero));
    assert_eq!(matrix.rows, before.rows);
    assert_eq!(matrix.col_weights, before.col_weights);
}

/// Exhausted budgets are resource outcomes, reported with what the search saw.
///
/// On `x² − 3` with `u`, `s`: the first three primes 3, 5, 7 yield no pair. With four primes the
/// pair `(11, 6)` is found and `(11, 5)` rejected. A root budget of 6 evaluates the residues
/// 0..5, so mod 11 only the root 5 and mod 13 only the root 4 are found, while the roots 6 and 9
/// lie beyond the budget. The largest prime below 2^64 is 2^64 − 59, so a search that starts at
/// 2^64 − 58 meets the end of the machine range at once.
#[test]
fn kat_root_character_bounded_selection_reports_resource_outcomes() {
    let (_, relations) = unit_and_square();
    let exhausted = |first, primes, evaluations, count| {
        let limits = RootCharacterLimits {
            prime_search: PrimeSearchLimits { first_candidate: first, max_primes_examined: primes },
            root_search: RootSearchLimits { max_evaluations: evaluations },
        };
        match select_root_characters(&relations, count, limits) {
            Err(RootCharacterError::Exhausted { resource, requested, found, report }) => {
                (resource, requested, found, *report)
            }
            other => panic!("expected exhaustion, got {other:?}"),
        }
    };

    let (resource, requested, found, report) = exhausted(3, 3, 100_000, 2);
    assert_eq!(resource, CharacterResource::PrimeBudget { limit: 3 });
    assert_eq!((requested, found, report.primes_examined), (2, 0, 3));
    assert_eq!(report.exceptional_primes, vec![3]);

    let (resource, _, found, report) = exhausted(3, 4, 100_000, 2);
    assert_eq!(resource, CharacterResource::PrimeBudget { limit: 4 });
    assert_eq!(found, 1);
    assert_eq!(report.rejected_pairs, vec![RejectedPair { prime: 11, root: 5, relation: 1 }]);

    let (resource, _, found, report) = exhausted(3, 5, 6, 2);
    assert_eq!(resource, CharacterResource::PrimeBudget { limit: 5 });
    assert_eq!(found, 1, "only (13, 4) is offered; (11, 5) vanishes and the other roots are cut");
    assert_eq!(report.truncated_primes, vec![11, 13]);

    let (resource, _, found, report) = exhausted(u64::MAX - 57, 10, 100_000, 1);
    assert_eq!(resource, CharacterResource::PrimeRange);
    assert_eq!((found, report.primes_examined), (0, 0));

    // Zero requested pairs succeed even with no budget at all.
    let none = RootCharacterLimits {
        prime_search: PrimeSearchLimits { first_candidate: 3, max_primes_examined: 0 },
        root_search: RootSearchLimits { max_evaluations: 0 },
    };
    assert!(select_root_characters(&relations, 0, none).unwrap().characters().is_empty());
}

/// The corpus for `x² − 3`, `m = 7`, `N = 23` with factor base bounds `(41, 11)`: every coprime
/// pair with `1 ≤ b ≤ 7`, `|a| ≤ 60`, `|a − 7b|` 41-smooth, and `|a² − 3b²|` `{2, 3, 11}`-smooth,
/// listed by `b` then `a`. Found and checked with an independent script.
const CORPUS: [(i64, i64); 43] = [
    (-27, 1), (-6, 1), (-5, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 1),
    (5, 1), (6, 1), (27, 1), (-3, 2), (-1, 2), (1, 2), (3, 2), (-7, 3), (-5, 3), (-4, 3),
    (4, 3), (5, 3), (7, 3), (-13, 4), (-9, 4), (-7, 4), (7, 4), (9, 4), (13, 4), (-14, 5),
    (-9, 5), (-3, 5), (3, 5), (8, 5), (9, 5), (14, 5), (-53, 7), (-13, 7), (-9, 7), (9, 7),
    (12, 7), (13, 7), (53, 7),
];

/// Characters of merged rows equal the XOR of the original rows' characters.
///
/// Fixture: the 43-pair corpus. The search over primes 3..37 rejects `(11, 5)` (it vanishes on
/// position 1, `−6 − 5 = −11`) and `(11, 6)` (position 0, `−27 − 6 = −33`), and returns the pairs
/// below. The expected bits are recomputed here from the list of squares and the original pairs
/// only, XORed over each row's provenance, never from a merged norm.
#[test]
fn kat_root_character_merged_rows_equal_xor_of_original_bits() {
    let (_, relations) = checked(&F_UNIT, 7, 23, (41, 11), &CORPUS);
    let selection = select_root_characters(&relations, 6, RootCharacterLimits::default()).unwrap();
    let pairs = [(13, 4), (13, 9), (23, 7), (23, 16), (37, 15), (37, 22)];
    let chosen: Vec<(u64, u64)> =
        selection.characters().characters().iter().map(|c| (c.prime(), c.root())).collect();
    assert_eq!(chosen, pairs);
    assert_eq!(
        selection.report().rejected_pairs,
        vec![
            RejectedPair { prime: 11, root: 5, relation: 1 },
            RejectedPair { prime: 11, root: 6, relation: 0 },
        ]
    );
    assert_eq!(selection.report().primes_examined, 11);
    let characters = selection.characters();
    let expected_bits = |provenance: &[usize]| -> Vec<bool> {
        pairs
            .iter()
            .map(|&(q, s)| {
                provenance.iter().fold(false, |acc, &i| {
                    let (a, b) = CORPUS[i];
                    acc ^ nonresidue(a - b * s as i64, q).expect("no valid pair vanishes")
                })
            })
            .collect()
    };

    let (plain, layout) = build_checked_matrix(&relations, characters.len()).unwrap();

    // Filling before filtering, then merging by elimination.
    let mut filled = plain.clone();
    characters.fill_columns(&mut filled, &layout, &relations).unwrap();
    for (i, row) in filled.rows.iter().enumerate() {
        assert_eq!(character_bits(&filled, &layout, i), expected_bits(&row.provenance));
    }
    let merged_before = merge_columns(filled);

    // Filling after filtering: elimination never touches character columns, and the bits follow
    // from the complete provenance alone.
    let mut merged_after = merge_columns(plain);
    assert!(merged_after.rows.iter().all(|r| character_bits_absent(r, &layout)));
    characters.fill_columns(&mut merged_after, &layout, &relations).unwrap();

    assert!(
        merged_after.rows.iter().any(|r| r.provenance.len() >= 2),
        "elimination merged at least one pair of rows"
    );
    // The pivot of a weight-three elimination is the lightest incident row, and filled rows are
    // heavier, so the two matrices need not be identical; each must satisfy the identity.
    for matrix in [&merged_before, &merged_after] {
        for (i, row) in matrix.rows.iter().enumerate() {
            assert_eq!(character_bits(matrix, &layout, i), expected_bits(&row.provenance));
        }
        characters.verify_columns(matrix, &layout, &relations).unwrap();
        for col in 0..matrix.num_cols {
            let weight = matrix.rows.iter().filter(|r| r.cols.contains(&col)).count();
            assert_eq!(matrix.col_weights[col] as usize, weight, "weight of column {col}");
        }
    }

    // Dropping an original index from a multi-relation provenance is detected when its bits are
    // not all zero, so the identity genuinely depends on the complete provenance.
    let table = characters.bit_table(&relations).unwrap();
    let mut detected = false;
    for (r, row) in merged_after.rows.iter().enumerate() {
        for &dropped in &row.provenance {
            let Some(slot) = table[dropped].iter().position(|&b| b) else { continue };
            let mut corrupt = merged_after.clone();
            corrupt.rows[r].provenance.retain(|&i| i != dropped);
            // The row's bits are the XOR with the dropped relation, so they differ from the
            // shortened provenance exactly where that relation has a set bit.
            assert_eq!(
                characters.verify_columns(&corrupt, &layout, &relations),
                Err(RootCharacterError::CharacterMismatch { row: r, slot })
            );
            detected = true;
        }
    }
    assert!(detected);
}

/// True iff `row` has no bit in the character columns of `layout`.
fn character_bits_absent(row: &MatrixRow, layout: &MatrixLayout) -> bool {
    row.cols.iter().all(|&c| c <= layout.sign_column())
}

/// Overlapping provenance cancels: `[0, 1] ⊕ [1, 2] = [0, 2]`, and a character column is the XOR
/// of the original relations 0 and 2 only. Relation 1 has nonzero bits `[0, 1, 1, 0, 0, 0]`, so a
/// union would give a different answer. Bits of 0, 1, 2 on the corpus pairs `(13,4), (13,9),
/// (23,7), (23,16), (37,15), (37,22)`: `[1,0,0,0,1,0]`, `[0,1,1,0,0,0]`, `[0,0,1,0,1,0]`.
#[test]
fn kat_root_character_overlapping_provenance_cancels() {
    let (_, relations) = checked(&F_UNIT, 7, 23, (41, 11), &CORPUS);
    let pairs = [pair(13, 4), pair(13, 9), pair(23, 7), pair(23, 16), pair(37, 15), pair(37, 22)];
    let characters = RootCharacters::new(relations.checker(), &pairs).unwrap();
    let (plain, layout) = build_checked_matrix(&relations, characters.len()).unwrap();
    let r01 = plain.rows[0].xor_merge(&plain.rows[1]);
    let r12 = plain.rows[1].xor_merge(&plain.rows[2]);
    let r02 = r01.xor_merge(&r12);
    assert_eq!(r02.provenance, vec![0, 2]);
    let mut matrix = matrix_with_rows(&plain, vec![r01, r12, r02]);
    characters.fill_columns(&mut matrix, &layout, &relations).unwrap();
    assert_eq!(character_bits(&matrix, &layout, 0), [true, true, true, false, true, false]);
    assert_eq!(character_bits(&matrix, &layout, 1), [false, true, false, false, true, false]);
    // [1,0,0,0,1,0] ⊕ [0,0,1,0,1,0], not the union's extra [0,1,1,0,0,0].
    assert_eq!(character_bits(&matrix, &layout, 2), [true, false, true, false, false, false]);
    // Row 2 is also the XOR of rows 0 and 1 in every character column, as a merged row must be.
    for slot in 0..6 {
        let bit = |row| character_bits(&matrix, &layout, row)[slot];
        assert_eq!(bit(2), bit(0) ^ bit(1));
    }
    characters.verify_columns(&matrix, &layout, &relations).unwrap();
}

/// Mismatched checkers, layouts, matrices, and malformed rows are refused.
///
/// Two factor bases, bounds `(41, 11)` and `(37, 11)`, give unequal checkers on `x² − 3`; the pair
/// `(2, −1)` has norms 9 and 1 and is representable in both. The pairs `(13, 4)` and `(13, 9)` are
/// valid roots (16 = 13 + 3, 81 = 78 + 3) and do not vanish on it (`2 + 4 = 6`, `2 + 9 = 11`).
#[test]
fn kat_root_character_layout_and_checker_mismatch_rejected() {
    let (checker_a, relations_a) = unit_and_square();
    let (_, relations_b) = checked(&F_UNIT, 7, 23, (37, 11), &[(2, -1)]);
    assert_ne!(checker_a, *relations_b.checker());
    let characters = RootCharacters::new(&checker_a, &[pair(13, 4), pair(13, 9)]).unwrap();
    let (matrix_a, layout_a) = build_checked_matrix(&relations_a, 2).unwrap();

    // Characters against relations of the other checker.
    assert_eq!(
        characters.bit_table(&relations_b),
        Err(RootCharacterError::RelationsCheckerMismatch)
    );
    assert_eq!(
        characters.relation_bits(&relations_b, 0),
        Err(RootCharacterError::RelationsCheckerMismatch)
    );
    // A layout from the other checker.
    let (mut matrix_b, layout_b) = build_checked_matrix(&relations_b, 2).unwrap();
    assert_eq!(
        characters.fill_columns(&mut matrix_b, &layout_b, &relations_a),
        Err(RootCharacterError::LayoutCheckerMismatch)
    );
    // A layout with a different number of slots.
    let (mut three, layout_three) = build_checked_matrix(&relations_a, 3).unwrap();
    assert_eq!(
        characters.fill_columns(&mut three, &layout_three, &relations_a),
        Err(RootCharacterError::LayoutCharacterCount { layout: 3, characters: 2 })
    );
    // A matrix that does not have the layout's shape.
    let mut skewed = matrix_a.clone();
    skewed.obstruction_count += 1;
    assert_eq!(
        characters.fill_columns(&mut skewed, &layout_a, &relations_a),
        Err(RootCharacterError::MatrixShapeMismatch)
    );
    assert_eq!(
        characters.verify_columns(&skewed, &layout_a, &relations_a),
        Err(RootCharacterError::MatrixShapeMismatch)
    );

    // Malformed rows are refused without changing the matrix.
    let mut unsorted = matrix_a.clone();
    unsorted.rows[0].provenance = vec![1, 0];
    let mut outside = matrix_a.clone();
    outside.rows[1].provenance = vec![0, 7];
    let first = layout_a.character_column(0).unwrap();
    let mut occupied = matrix_a.clone();
    occupied.rows[0].cols.push(first);
    for (mut matrix, expected) in [
        (unsorted, RootCharacterError::MalformedProvenance { row: 0 }),
        (outside, RootCharacterError::ProvenanceOutOfRange { row: 1, index: 7, relations: 2 }),
        (occupied, RootCharacterError::CharacterColumnNotZero { row: 0, column: first }),
    ] {
        let before = matrix.rows.clone();
        assert_eq!(characters.fill_columns(&mut matrix, &layout_a, &relations_a), Err(expected));
        assert_eq!(matrix.rows, before);
    }

    // A matrix whose character columns were filled wrongly fails verification.
    let mut wrong = matrix_a.clone();
    characters.fill_columns(&mut wrong, &layout_a, &relations_a).unwrap();
    let table = characters.bit_table(&relations_a).unwrap();
    assert_eq!(table, vec![vec![true, true], vec![false, false]]);
    wrong.rows[1].cols.push(first);
    assert_eq!(
        characters.verify_columns(&wrong, &layout_a, &relations_a),
        Err(RootCharacterError::CharacterMismatch { row: 1, slot: 0 })
    );
}

/// Error values describe the failing item and expose their source.
#[test]
fn kat_root_character_errors_are_described() {
    let (_, relations) = unit_and_square();
    let message = |e: RootCharacterError| e.to_string();
    assert!(message(RootCharacterError::Pair {
        position: 2,
        prime: 3,
        root: 0,
        reason: PairRejection::Exceptional
    })
    .contains("pair 2"));
    assert!(message(RootCharacterError::ZeroValue { slot: 1, prime: 11, root: 5, relation: 1 })
        .contains("vanishes on relation 1"));
    let exhaustion = select_root_characters(
        &relations,
        2,
        RootCharacterLimits {
            prime_search: PrimeSearchLimits { first_candidate: 3, max_primes_examined: 1 },
            ..RootCharacterLimits::default()
        },
    )
    .unwrap_err();
    assert!(exhaustion.to_string().contains("found 0 of 2"));
    assert!(std::error::Error::source(&exhaustion).is_none());
    let range = CharacterResource::PrimeRange.to_string();
    assert!(range.contains("range"));
    for reason in [
        PairRejection::Even,
        PairRejection::NotPrime,
        PairRejection::Exceptional,
        PairRejection::RootOutOfRange,
        PairRejection::NotRoot,
        PairRejection::Duplicate { first: 1 },
    ] {
        assert!(!reason.to_string().is_empty());
    }
    for error in [
        RootCharacterError::RelationsCheckerMismatch,
        RootCharacterError::LayoutCheckerMismatch,
        RootCharacterError::LayoutCharacterCount { layout: 1, characters: 2 },
        RootCharacterError::MatrixShapeMismatch,
        RootCharacterError::RelationOutOfRange { index: 3, len: 2 },
        RootCharacterError::MalformedProvenance { row: 0 },
        RootCharacterError::ProvenanceOutOfRange { row: 0, index: 9, relations: 2 },
        RootCharacterError::CharacterColumnNotZero { row: 0, column: 4 },
        RootCharacterError::CharacterMismatch { row: 0, slot: 1 },
        RootCharacterError::Arithmetic(shared_numfield::FpPolyError::ZeroPolynomial),
    ] {
        assert!(!error.to_string().is_empty());
    }
    let arithmetic = RootCharacterError::Arithmetic(shared_numfield::FpPolyError::ZeroPolynomial);
    assert!(std::error::Error::source(&arithmetic).is_some());
}

/// The legacy helper is unchanged: `x³ − x − 1` has three distinct roots mod 59 (4, 13, 42: sum
/// 59 ≡ 0 and product 2184 = 37·59 + 1 ≡ 1, as Vieta requires) and mod 101 (20, 89, 93), 167
/// (73, 127, 134), and 173 (97, 110, 139). A brute-force count over `[14, 100)` finds 59 as the
/// only prime that splits completely there, so the first three primes above 13 are 59, 101, 167
/// and above 100 they are 101, 167, 173 (checked with an independent script).
#[test]
fn kat_legacy_select_qc_primes_unchanged() {
    let f = f_cubic();
    assert_eq!(select_qc_primes(&f, 13, 3), vec![59, 101, 167]);
    assert_eq!(select_qc_primes(&f, 100, 3), vec![101, 167, 173]);
    assert_eq!(select_qc_primes(&f, 13, 0), Vec::<u64>::new());
    let constant = IntPoly::from_coeffs(vec![bi(5)]);
    assert!(select_qc_primes(&constant, 13, 3).is_empty());
}

// ─── Helper: modular exponentiation (for KAT 2 hand computation) ─────────────

/// Compute `base^exp mod q` using fast exponentiation.
fn pow_mod_test(mut base: u64, mut exp: u64, q: u64) -> u64 {
    if q == 1 {
        return 0;
    }
    let mut result = 1u64;
    base %= q;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % q as u128) as u64;
        }
        base = ((base as u128 * base as u128) % q as u128) as u64;
        exp >>= 1;
    }
    result
}
