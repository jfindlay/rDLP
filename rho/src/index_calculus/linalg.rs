//! Z/ℓZ linear algebra for the index-calculus ECDLP solver.
//!
//! This module provides the relation→`FlSparseMatrix` adapter and the left-kernel solve
//! that form the linear-algebra step of the index-calculus pipeline. It reuses the
//! `gnfs::dl::linalg` engine over F_ℓ with an index-calculus-specific adapter in place of
//! the NFS-bound `build_fl_matrix(DLMatrix)`.
//!
//! # Matrix orientation
//!
//! Let `E` be the `#relations × |FB|` exponent matrix, `E[i][j] = e_ij`. The adapter builds
//! its transpose `Eᵀ`: one sparse row per factor-base point, one column per relation. The
//! engine returns right-kernel vectors `u` of `Eᵀ`, which are exactly the relation
//! combinations `u` with `Σᵢ uᵢ·eᵢ ≡ 0 (mod ℓ)`, where `eᵢ` is the exponent row of
//! relation `i`. Such a `u` cancels every unknown factor-base logarithm, leaving a linear
//! condition on `x = log_G(Q)` alone (see [`crate::index_calculus::solve`]).
//!
//! # Solver choice
//!
//! `solve_ek_linalg` uses `block_wiedemann_fl`, which at this scale row-reduces the dense
//! matrix exactly (deterministic, no probabilistic failure modes over the small field).
//! That routine returns one kernel vector per free column, i.e. a *complete* basis of the
//! kernel, and the caller relies on completeness: a statement that no relation combination
//! determines `x` is only meaningful when the basis spans the whole kernel.
//!
//! # Scale
//!
//! The toy fixture is `F_47`, `n = 60`, `ℓ = 5`, `m = 2`, `|FB| = 6`. Over `E(F_p)` this
//! demonstrates the mechanism; the asymptotic advantage of index calculus needs the
//! extension-field setting `E(F_{p^n})` and is not observable here.

use gnfs::dl::linalg::{FlMatrixOperator, FlSparseMatrix, FlSparseRow, block_wiedemann_fl};
use shared_field::FpNaive4 as FpNaive;
use crate::index_calculus::strategy::{IndexCalcStrategy, Relation};
use crate::index_calculus::IndexCalcError;

// ─── build_ek_matrix ─────────────────────────────────────────────────────────

/// Build the transposed relation matrix for the index-calculus linear algebra step.
///
/// Row `j` is the sparse row over F_ℓ of factor-base point `j`: its entries are
/// `(relation index i, e_ij)` for every relation `i` whose exponent of point `j` is
/// non-zero. The matrix has `strategy.fb_size()` rows and `relations.len()` columns, so
/// its right kernel is the set of relation combinations `u` with `Σᵢ uᵢ·eᵢ ≡ 0`.
///
/// Entries within a row are sorted by relation index, as `FlSparseRow` requires: relations
/// are visited in order, and each `Relation.exponents` has no duplicate factor-base index.
pub fn build_ek_matrix(relations: &[Relation], strategy: &IndexCalcStrategy) -> FlSparseMatrix<FpNaive> {
    let mut rows: Vec<FlSparseRow<FpNaive>> =
        (0..strategy.fb_size()).map(|_| FlSparseRow { entries: Vec::new() }).collect();

    for (i, rel) in relations.iter().enumerate() {
        for (j, exp) in &rel.exponents {
            rows[*j].entries.push((i, exp.clone()));
        }
    }

    FlSparseMatrix { rows, num_cols: relations.len() }
}

// ─── solve_ek_linalg ─────────────────────────────────────────────────────────

/// Compute a basis of the left kernel of the relation exponent matrix.
///
/// Each returned vector `u ∈ F_ℓ^{#relations}` satisfies `Σᵢ uᵢ·eᵢ ≡ 0 (mod ℓ)`, where
/// `eᵢ` is the exponent row of `relations[i]`. The basis is complete: it has dimension
/// `#relations − rank(E)` and spans every such combination. It is empty exactly when the
/// exponent rows are linearly independent (in particular for an empty relation set); that
/// is a valid answer, not an error.
///
/// The kernel is the right kernel of the transposed matrix from `build_ek_matrix`,
/// computed by the `gnfs::dl::linalg` `block_wiedemann_fl` engine.
///
/// # Errors
///
/// Currently infallible; the `Result` keeps the signature stable for engines that can
/// fail.
pub fn solve_ek_linalg(
    relations: &[Relation],
    strategy: &IndexCalcStrategy,
) -> Result<Vec<Vec<FpNaive>>, IndexCalcError> {
    let matrix = build_ek_matrix(relations, strategy);
    let ell = &strategy.ell;

    let op = FlMatrixOperator::<FpNaive, 4>::new(&matrix);

    // The rng_seed is unused by the exact row-reduction implementation.
    let solutions = block_wiedemann_fl(&op, ell, 42);

    Ok(solutions.into_iter().map(|sol| sol.coefficients).collect())
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::Fp;
    use crate::index_calculus::strategy::IndexCalcStrategy;

    fn ell() -> Uint<4> {
        Uint::<4>::from(crate::index_calculus::strategy::TOY_ELL)
    }

    #[test]
    fn build_ek_matrix_dimensions() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        // Construct two minimal relations.
        let r0 = Relation::from_decomposition(1, 0, &[0, 1], &ell());
        let r1 = Relation::from_decomposition(2, 0, &[1, 2], &ell());
        let relations = vec![r0, r1];

        let matrix = build_ek_matrix(&relations, &strategy);

        assert_eq!(matrix.rows.len(), strategy.fb_size(), "one row per factor-base point");
        assert_eq!(matrix.num_cols, 2, "one column per relation");
    }

    #[test]
    fn build_ek_matrix_is_transpose_of_exponents() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        // Relation 0 hits points {0, 2}; relation 1 hits {2, 2} (exponent 2 on point 2).
        let r0 = Relation::from_decomposition(1, 0, &[0, 2], &ell());
        let r1 = Relation::from_decomposition(2, 0, &[2, 2], &ell());
        let matrix = build_ek_matrix(&[r0.clone(), r1.clone()], &strategy);

        let entries = |j: usize| -> Vec<(usize, FpNaive)> { matrix.rows[j].entries.clone() };
        assert_eq!(entries(0), vec![(0, r0.exponent(0, &ell()))]);
        assert_eq!(
            entries(2),
            vec![(0, r0.exponent(2, &ell())), (1, r1.exponent(2, &ell()))],
            "row of point 2 lists its exponent in each relation, in relation order"
        );
        assert!(entries(1).is_empty(), "a point in no relation has an empty row");
    }

    #[test]
    fn build_ek_matrix_empty_relations() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let matrix = build_ek_matrix(&[], &strategy);
        assert_eq!(matrix.rows.len(), strategy.fb_size(), "one (empty) row per point");
        assert_eq!(matrix.num_cols, 0, "no relations → no columns");
    }

    #[test]
    fn left_kernel_of_duplicate_rows() {
        // Two identical exponent rows are dependent: the left kernel is spanned by their
        // difference; a third independent row adds nothing.
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let ell = ell();
        let rels = vec![
            Relation::from_decomposition(1, 0, &[0, 1], &ell),
            Relation::from_decomposition(2, 0, &[0, 1], &ell),
            Relation::from_decomposition(3, 0, &[2, 3], &ell),
        ];
        let basis = solve_ek_linalg(&rels, &strategy).unwrap();
        assert_eq!(basis.len(), 1, "kernel dimension is #relations − rank = 3 − 2");
        let u = &basis[0];
        let minus_one = FpNaive::one(&ell).neg(&ell);
        let expected = vec![FpNaive::one(&ell), minus_one, FpNaive::zero(&ell)];
        assert!(
            *u == expected || *u == expected.iter().map(|c| c.neg(&ell)).collect::<Vec<_>>(),
            "kernel vector should be ±(1, −1, 0): {u:?}"
        );
    }

    #[test]
    fn left_kernel_empty_for_independent_rows() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let rels = vec![
            Relation::from_decomposition(1, 0, &[0, 1], &ell()),
            Relation::from_decomposition(2, 0, &[2, 3], &ell()),
        ];
        assert!(solve_ek_linalg(&rels, &strategy).unwrap().is_empty());
        assert!(solve_ek_linalg(&[], &strategy).unwrap().is_empty());
    }
}
