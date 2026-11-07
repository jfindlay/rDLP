//! Known-answer tests for checked dependencies: binding a filtered matrix to its original
//! relations, validating candidate row selections, and the exact rational and algebraic products.
//!
//! Notation: a relation `(a, b)` has rational norm `a − b·m` and homogeneous norm
//! `f_hom(a, b) = Σ f_i a^i b^(d−i)`. For leading coefficient `c`, the monic form is
//! `F(t) = c^(d−1) f(t/c)`, `θ = c·α`, and `γ = Π(c·a − b·θ) / c^k` for `k` relations.
//!
//! **Fixture T (non-monic).** `f = 3x² + x − 1`, `m = 3`, `N = 29` (`f(3) = 29`), `c = 3`,
//! `F = t² + t − 3` so `θ² = 3 − θ`, `disc(F) = 13`. `FactorBase::new(f, 7, 17)` gives rational
//! primes `[2, 3, 5, 7]` (columns 0–3) and ideals `(3,1) (13,2) (17,4) (17,7)` (columns 4–7); 3 and
//! 13 are exceptional (`3 | c`, `13 | disc F`), so columns 4 and 5 stay zero. The sign is column 8
//! and character slot `k` is column `9 + k`. Its six relations, in row order, are
//!
//! | row | (a, b) | a − 3b | f_hom | checked columns (K = 0) |
//! |-----|--------|--------|-------|------------------------|
//! | 0 | (4, 1)   | 1   | 51 = 3·17     | [6]          |
//! | 1 | (7, 1)   | 4   | 153 = 3²·17   | [7]          |
//! | 2 | (2, 1)   | −1  | 13            | [8]          |
//! | 3 | (−3, 2)  | −9  | 17            | [7, 8]       |
//! | 4 | (−9, 2)  | −15 | 221 = 13·17   | [1, 2, 6, 8] |
//! | 5 | (5, 2)   | −1  | 81 = 3⁴       | [8]          |
//!
//! Every product, norm, character bit, and column in this file was recomputed with an independent
//! exact-arithmetic script (Python `fractions` and `sympy`, no Rust involved) before being
//! asserted; the derivation is written beside each use.

use gnfs::filter::merge::eliminate_column;
use gnfs::filter::{
    MatrixLayout, SparseMatrix, build_checked_matrix, merge_columns, remove_singletons,
};
use gnfs::linalg::{
    KernelVector, MatrixOperator, RootCharacterError, RootCharacters, block_lanczos,
};
use gnfs::polyselect::PolyPair;
use gnfs::sieve::{
    CheckedRelations, CheckedRelationsError, ExponentVector, FactorBase, LineSieveConfig,
    Relation, RelationChecker, line_sieve,
};
use gnfs::sqrt::{
    CheckedDependency, DependencyContext, DependencyContextError, DependencyError,
    DependencyRejection, FailureClass, PrecisionLimits, SqrtContext, TargetOutcome,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use shared_numfield::{IntPoly, RatPoly};

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

/// The rational polynomial with coefficients `n_i / d_i`, lowest degree first.
fn rat(coeffs: &[(i64, i64)]) -> RatPoly {
    RatPoly::from_coeffs(coeffs.iter().map(|&(n, d)| BigRational::new(big(n), big(d))).collect())
}

fn relation(
    a: i64,
    b: i64,
    rational: &[(usize, u32)],
    algebraic: &[(usize, u32)],
    negative: bool,
) -> Relation {
    Relation {
        a: big(a),
        b: big(b),
        rational_exponents: ExponentVector { entries: rational.to_vec() },
        algebraic_exponents: ExponentVector { entries: algebraic.to_vec() },
        rational_sign: negative,
    }
}

const F_T: [i64; 3] = [-1, 1, 3];

fn context_t() -> SqrtContext {
    SqrtContext::new(&poly(&F_T), &big(3), &big(29)).expect("fixture T is a valid context")
}

fn checker_with(fb: &FactorBase) -> RelationChecker {
    RelationChecker::new(&context_t(), fb).expect("fixture T factor base is valid")
}

fn checker_t() -> RelationChecker {
    checker_with(&FactorBase::new(&poly(&F_T), 7, 17))
}

/// Fixture-T relations with their hand-derived exponent vectors (rational index: 2→0, 3→1, 5→2,
/// 7→3; ideal index: (3,1)→0, (13,2)→1, (17,4)→2, (17,7)→3). Labels: `a ≡ b·r (mod p)`.
///
/// - (4,1): 1; 51 = 3·17; labels (3,1) (4 ≡ 1) and (17,4).
/// - (7,1): 4 = 2²; 153 = 3²·17; labels (3,1) and (17,7).
/// - (2,1): −1; 13; label (13,2) (2 ≡ 2).
/// - (−3,2): −9 = −3²; 17; (17,7) because 2·7 = 14 ≡ −3.
/// - (−9,2): −15 = −3·5; 221 = 13·17; (13,2) because 2·2 = 4 ≡ −9, (17,4) because 8 ≡ −9.
/// - (5,2): −1; 81 = 3⁴; (3,1) because 5 ≡ 2·1.
fn t_relations() -> Vec<Relation> {
    vec![
        relation(4, 1, &[], &[(0, 1), (2, 1)], false),
        relation(7, 1, &[(0, 2)], &[(0, 2), (3, 1)], false),
        relation(2, 1, &[], &[(1, 1)], true),
        relation(-3, 2, &[(1, 2)], &[(3, 1)], true),
        relation(-9, 2, &[(1, 1), (2, 1)], &[(1, 1), (2, 1)], true),
        relation(5, 2, &[], &[(0, 4)], true),
    ]
}

fn checked_t() -> CheckedRelations {
    CheckedRelations::new(&checker_t(), &t_relations()).expect("fixture T relations are valid")
}

/// Everything a context is built from, with `pairs` as the root characters.
struct Parts {
    relations: CheckedRelations,
    matrix: SparseMatrix,
    layout: MatrixLayout,
    characters: RootCharacters,
}

impl Parts {
    fn build(&self) -> Result<DependencyContext, DependencyContextError> {
        DependencyContext::new(
            self.relations.clone(),
            self.matrix.clone(),
            self.layout.clone(),
            self.characters.clone(),
        )
    }
}

/// Fixture T with `pairs` as characters. The matrix is the unfilled output of
/// `build_checked_matrix` with `pairs.len()` character slots.
fn parts_t(pairs: &[(u64, u64)]) -> Parts {
    parts_for(checked_t(), pairs)
}

fn parts_for(relations: CheckedRelations, pairs: &[(u64, u64)]) -> Parts {
    let (matrix, layout) = build_checked_matrix(&relations, pairs.len()).unwrap();
    let characters = RootCharacters::new(relations.checker(), pairs).unwrap();
    Parts { relations, matrix, layout, characters }
}

/// Like [`parts_t`], with the character columns filled from the original relations.
fn filled_t(pairs: &[(u64, u64)]) -> Parts {
    let mut parts = parts_t(pairs);
    parts
        .characters
        .fill_columns(&mut parts.matrix, &parts.layout, &parts.relations)
        .expect("characters do not vanish on fixture T");
    parts
}

/// Two characters at q = 23, derived by the independent script: 23 ∤ 3·13, and
/// f(18) = 3·324 + 18 − 1 = 989 = 43·23, f(20) = 3·400 + 20 − 1 = 1219 = 53·23, so 18 and 20 are
/// roots mod 23. Legendre symbols of `a − b·s` (residues mod 23 are 1,2,3,4,6,8,9,12,13,16,18):
///
/// | (a, b)   | s = 18: a − 18b mod 23 | symbol | s = 20: a − 20b mod 23 | symbol |
/// |----------|------------------------|--------|------------------------|--------|
/// | (4, 1)   | −14 ≡ 9  | +1 (bit 0) | −16 ≡ 7  | −1 (bit 1) |
/// | (7, 1)   | −11 ≡ 12 | +1 (0)     | −13 ≡ 10 | −1 (1)     |
/// | (2, 1)   | −16 ≡ 7  | −1 (1)     | −18 ≡ 5  | −1 (1)     |
/// | (−3, 2)  | −39 ≡ 7  | −1 (1)     | −43 ≡ 3  | +1 (0)     |
/// | (−9, 2)  | −45 ≡ 1  | +1 (0)     | −49 ≡ 20 | −1 (1)     |
/// | (5, 2)   | −31 ≡ 15 | −1 (1)     | −35 ≡ 11 | −1 (1)     |
const CHARS_23: [(u64, u64); 2] = [(23, 18), (23, 20)];

/// The original columns of fixture T with `CHARS_23` as slots 0 and 1 (columns 9 and 10): the K = 0
/// columns plus 9 for a set slot-0 bit and 10 for a set slot-1 bit, from the table above.
const T_ORIGINAL_K2: [&[usize]; 6] =
    [&[6, 10], &[7, 10], &[8, 9, 10], &[7, 8, 9], &[1, 2, 6, 8, 10], &[8, 9, 10]];

const T_ORIGINAL_K0: [&[usize]; 6] = [&[6], &[7], &[8], &[7, 8], &[1, 2, 6, 8], &[8]];

fn dependency_context(parts: &Parts) -> DependencyContext {
    parts.build().expect("genuine parts form a context")
}

fn rejected(rejection: DependencyRejection) -> DependencyError {
    DependencyError::Rejected(rejection)
}

fn precision() -> PrecisionLimits {
    PrecisionLimits { max_modulus_bits: 4096 }
}

/// The normalized target's `γ` and `P`, computed by a second implementation.
fn assert_target_agrees(dependency: &CheckedDependency) {
    let TargetOutcome::Ready(target) = dependency.target(precision()).unwrap() else {
        panic!("a 4096-bit budget is ample for the fixtures");
    };
    assert_eq!(target.gamma(), dependency.gamma());
    assert_eq!(target.rational_product(), dependency.rational_product());
    assert_eq!(target.factors(), dependency.factors());
}

// ─── binding ─────────────────────────────────────────────────────────────────

#[test]
fn kat_context_binds_relations_layout_characters() {
    // K = 0: the original rows are the checked columns of the table in the module header.
    let parts = parts_t(&[]);
    let context = dependency_context(&parts);
    assert_eq!(context.context(), &context_t());
    assert_eq!(context.relations(), &parts.relations);
    assert_eq!(context.layout(), &parts.layout);
    assert_eq!(context.characters(), &parts.characters);
    assert_eq!(context.matrix().rows.len(), 6);
    for (i, expected) in T_ORIGINAL_K0.iter().enumerate() {
        assert_eq!(context.original_columns(i), Some(*expected), "relation {i}");
        let from_layout = parts.layout.row_columns(&parts.relations.relations()[i]).unwrap();
        assert_eq!(context.original_columns(i), Some(from_layout.as_slice()));
    }
    assert_eq!(context.original_columns(6), None);

    // K = 2: slot 0 is (23, 18), slot 1 is (23, 20). Each original row gains column 9 for a set
    // slot-0 bit and column 10 for a set slot-1 bit (the table at CHARS_23).
    let parts = filled_t(&CHARS_23);
    let context = dependency_context(&parts);
    assert_eq!(context.layout().character_count(), 2);
    assert_eq!(context.characters().len(), 2);
    for (i, expected) in T_ORIGINAL_K2.iter().enumerate() {
        assert_eq!(context.original_columns(i), Some(*expected), "relation {i}");
        // The matrix row of a freshly built relation is its original row.
        assert_eq!(context.matrix().rows[i].cols.as_slice(), *expected);
        // The character part equals the bits the characters report.
        let bits = parts.characters.relation_bits(&parts.relations, i).unwrap();
        let columns: Vec<usize> = (0..2).filter(|&s| bits[s]).map(|s| 9 + s).collect();
        let tail: Vec<usize> = expected.iter().copied().filter(|&c| c >= 9).collect();
        assert_eq!(columns, tail);
    }
}

#[test]
fn kat_context_accepts_characters_filled_before_or_after_filtering() {
    // Column 6 has weight 2 (rows 0 and 4); eliminating it XORs those two rows. Filling the
    // character columns after the elimination (by provenance XOR) and before it (carried through
    // by xor_merge) must both satisfy the provenance identity.
    let mut after = parts_t(&CHARS_23);
    eliminate_column(&mut after.matrix, 6).unwrap();
    after
        .characters
        .fill_columns(&mut after.matrix, &after.layout, &after.relations)
        .unwrap();

    let mut before = filled_t(&CHARS_23);
    eliminate_column(&mut before.matrix, 6).unwrap();

    for parts in [&after, &before] {
        let context = dependency_context(parts);
        // The merged row is original row 0 ⊕ original row 4, with provenance [0, 4]:
        // [6, 10] ⊕ [1, 2, 6, 8, 10] = [1, 2, 8].
        let merged = context.matrix().rows.last().unwrap();
        assert_eq!(merged.provenance, vec![0, 4]);
        assert_eq!(merged.cols, vec![1, 2, 8]);
    }
}

#[test]
fn kat_context_rejects_binding_mismatches() {
    let parts = parts_t(&[]);
    // A second factor base adds the rational prime 11, so its checker and every layout derived
    // from it differ from the first. The six relations are valid under both (all norms are
    // 7-smooth and the ideals have the same indices).
    let other_checker = checker_with(&FactorBase::new(&poly(&F_T), 11, 17));
    assert_ne!(other_checker, checker_t());
    let other_relations = CheckedRelations::new(&other_checker, &t_relations()).unwrap();
    let other = parts_for(other_relations, &[]);

    let mut foreign_layout = parts_t(&[]);
    foreign_layout.layout = other.layout.clone();
    assert_eq!(foreign_layout.build().unwrap_err(), DependencyContextError::LayoutRelationsMismatch);

    let mut foreign_characters = parts_t(&[]);
    foreign_characters.characters = other.characters.clone();
    assert_eq!(
        foreign_characters.build().unwrap_err(),
        DependencyContextError::CharactersRelationsMismatch
    );

    // Both foreign: the layout is checked first.
    let mut both = parts_t(&[]);
    both.layout = other.layout.clone();
    both.characters = other.characters.clone();
    assert_eq!(both.build().unwrap_err(), DependencyContextError::LayoutRelationsMismatch);

    // Character count: a layout with 2 slots against 0 characters, and the reverse.
    let mut wide_layout = parts_t(&CHARS_23);
    wide_layout.characters = parts.characters.clone();
    assert_eq!(
        wide_layout.build().unwrap_err(),
        DependencyContextError::CharacterCountMismatch { layout: 2, characters: 0 }
    );
    let mut wide_characters = parts_t(&[]);
    wide_characters.characters = parts_t(&CHARS_23).characters;
    assert_eq!(
        wide_characters.build().unwrap_err(),
        DependencyContextError::CharacterCountMismatch { layout: 0, characters: 2 }
    );

    // Shape: a truncated weight table, a wider matrix, and a matrix built for K = 2.
    let mut truncated = parts_t(&[]);
    truncated.matrix.col_weights.pop();
    assert_eq!(truncated.build().unwrap_err(), DependencyContextError::MatrixShapeMismatch);
    let mut wider = parts_t(&[]);
    wider.matrix.num_cols += 1;
    assert_eq!(wider.build().unwrap_err(), DependencyContextError::MatrixShapeMismatch);
    let mut wrong_k = parts_t(&[]);
    wrong_k.matrix = parts_t(&CHARS_23).matrix;
    assert_eq!(wrong_k.build().unwrap_err(), DependencyContextError::MatrixShapeMismatch);
}

// ─── per-row structure ───────────────────────────────────────────────────────

/// Apply `corrupt` to a genuine K = 0 matrix and return the construction error.
fn corrupted(corrupt: impl FnOnce(&mut SparseMatrix)) -> DependencyContextError {
    let mut parts = parts_t(&[]);
    corrupt(&mut parts.matrix);
    parts.build().expect_err("a corrupted matrix is never a context")
}

#[test]
fn kat_context_rejects_malformed_rows() {
    // Genuine rows: 0 [6], 1 [7], 2 [8], 3 [7,8], 4 [1,2,6,8], 5 [8]; provenance [i]; 9 columns.
    assert_eq!(
        corrupted(|m| m.rows[4].cols = vec![2, 1, 6, 8]),
        DependencyContextError::MalformedColumns { row: 4 },
        "unsorted columns"
    );
    assert_eq!(
        corrupted(|m| m.rows[4].cols = vec![1, 1, 2, 6, 8]),
        DependencyContextError::MalformedColumns { row: 4 },
        "duplicate columns"
    );
    assert_eq!(
        corrupted(|m| m.rows[0].cols = vec![6, 9]),
        DependencyContextError::ColumnOutOfRange { row: 0, column: 9 },
        "column 9 = num_cols"
    );
    assert_eq!(
        corrupted(|m| m.rows[2].provenance = vec![]),
        DependencyContextError::EmptyProvenance { row: 2 }
    );
    assert_eq!(
        corrupted(|m| m.rows[3].provenance = vec![3, 0]),
        DependencyContextError::MalformedProvenance { row: 3 },
        "unsorted provenance"
    );
    assert_eq!(
        corrupted(|m| m.rows[3].provenance = vec![3, 3]),
        DependencyContextError::MalformedProvenance { row: 3 },
        "duplicate provenance"
    );
    assert_eq!(
        corrupted(|m| m.rows[1].provenance = vec![1, 6]),
        DependencyContextError::ProvenanceOutOfRange { row: 1, index: 6, relations: 6 },
        "6 relations, so index 6 is out of range"
    );
    assert_eq!(
        corrupted(|m| m.col_weights[8] += 1),
        DependencyContextError::ColumnWeightMismatch { column: 8 }
    );

    // The first failing row is reported, whatever follows it.
    assert_eq!(
        corrupted(|m| {
            m.rows[5].provenance = vec![];
            m.rows[2].cols = vec![8, 3];
            m.rows[4].cols = vec![2, 1];
        }),
        DependencyContextError::MalformedColumns { row: 2 }
    );
    // Within one row the structural checks run in order: columns, provenance, then parity.
    assert_eq!(
        corrupted(|m| {
            m.rows[1].cols = vec![7, 7];
            m.rows[1].provenance = vec![];
        }),
        DependencyContextError::MalformedColumns { row: 1 }
    );
    assert_eq!(
        corrupted(|m| {
            m.rows[1].cols = vec![7, 12];
            m.rows[1].provenance = vec![];
        }),
        DependencyContextError::ColumnOutOfRange { row: 1, column: 12 }
    );
    // A row error precedes a weight error; with several wrong weights the smallest column wins.
    assert_eq!(
        corrupted(|m| {
            m.col_weights[8] += 1;
            m.rows[3].provenance = vec![];
        }),
        DependencyContextError::EmptyProvenance { row: 3 }
    );
    assert_eq!(
        corrupted(|m| {
            m.col_weights[8] += 1;
            m.col_weights[3] = 1;
        }),
        DependencyContextError::ColumnWeightMismatch { column: 3 }
    );
}

#[test]
fn kat_context_rejects_lifted_parity_mismatch() {
    // A flipped bit in row 0 (column 0 set) with its weight adjusted so the weight table is
    // honest: only the parity identity can catch it. Expected row 0 is [6]; the smallest
    // differing column is 0.
    assert_eq!(
        corrupted(|m| {
            m.rows[0].cols = vec![0, 6];
            m.col_weights[0] = 1;
        }),
        DependencyContextError::LiftedParityMismatch { row: 0, column: 0 }
    );
    // A bit in an omitted exceptional column (4 = ideal (3,1), 3 | c) must stay zero.
    assert_eq!(
        corrupted(|m| {
            m.rows[2].cols = vec![4, 8];
            m.col_weights[4] = 1;
        }),
        DependencyContextError::LiftedParityMismatch { row: 2, column: 4 }
    );
    // Provenance naming a different relation: row 1 would be original row 0 = [6], but it has
    // [7]. The smallest differing column of {6} Δ {7} is 6.
    assert_eq!(
        corrupted(|m| m.rows[1].provenance = vec![0]),
        DependencyContextError::LiftedParityMismatch { row: 1, column: 6 }
    );

    // Unfilled character columns at K = 2: original row 0 is [6, 10] but the fresh row is [6],
    // so the smallest differing column is 10 (slot-0 bit of (4,1) is 0 and slot-1 bit is 1).
    let unfilled = parts_t(&CHARS_23);
    assert_eq!(
        unfilled.build().unwrap_err(),
        DependencyContextError::LiftedParityMismatch { row: 0, column: 10 }
    );
    // A single wrong character bit in an otherwise filled matrix: drop column 9 from row 3.
    let mut flipped = filled_t(&CHARS_23);
    flipped.matrix.rows[3].cols = vec![7, 8];
    flipped.matrix.col_weights[9] -= 1;
    assert_eq!(
        flipped.build().unwrap_err(),
        DependencyContextError::LiftedParityMismatch { row: 3, column: 9 }
    );

    // A character that vanishes on a relation has no bit. (17, 4) is a valid pair for T
    // (f(4) = 51 = 3·17, 17 ∤ 39), but a − b·s = 4 − 1·4 = 0 on relation (4,1) = row 0.
    let vanishing = parts_t(&[(17, 4)]);
    assert_eq!(
        vanishing.build().unwrap_err(),
        DependencyContextError::Characters(RootCharacterError::ZeroValue {
            slot: 0,
            prime: 17,
            root: 4,
            relation: 0,
        })
    );
}

// ─── candidate selections ────────────────────────────────────────────────────

#[test]
fn kat_validate_rejects_malformed_selection() {
    let parts = parts_t(&[]);
    let context = dependency_context(&parts);
    let rows = 6;

    assert_eq!(context.validate(&[]), Err(rejected(DependencyRejection::EmptySelection)));
    assert_eq!(
        context.validate(&[2, 1]),
        Err(rejected(DependencyRejection::MalformedSelection { position: 1 }))
    );
    assert_eq!(
        context.validate(&[1, 1]),
        Err(rejected(DependencyRejection::MalformedSelection { position: 1 }))
    );
    assert_eq!(
        context.validate(&[0, 5, 5]),
        Err(rejected(DependencyRejection::MalformedSelection { position: 2 }))
    );
    assert_eq!(
        context.validate(&[0, 6]),
        Err(rejected(DependencyRejection::RowOutOfRange { position: 1, row: 6, rows }))
    );
    assert_eq!(
        context.validate(&[usize::MAX]),
        Err(rejected(DependencyRejection::RowOutOfRange { position: 0, row: usize::MAX, rows }))
    );
    // Order of checks: an unsorted list is reported before a range error in the same list.
    assert_eq!(
        context.validate(&[9, 2]),
        Err(rejected(DependencyRejection::MalformedSelection { position: 1 }))
    );

    // A struct literal that bypasses KernelVector::new's assertion: no panic, a typed rejection.
    let literal = KernelVector { row_indices: vec![5, 2] };
    assert_eq!(
        context.validate_kernel_vector(&literal),
        Err(rejected(DependencyRejection::MalformedSelection { position: 1 }))
    );
    let empty = KernelVector { row_indices: vec![] };
    assert_eq!(
        context.validate_kernel_vector(&empty),
        Err(rejected(DependencyRejection::EmptySelection))
    );
    // An out-of-range literal is a typed rejection, not an indexing panic.
    assert_eq!(
        context.validate_kernel_vector(&KernelVector { row_indices: vec![99] }),
        Err(rejected(DependencyRejection::RowOutOfRange { position: 0, row: 99, rows }))
    );

    // Classes.
    assert_eq!(
        context.validate(&[]).unwrap_err().class(),
        FailureClass::InvalidInput,
        "a rejection is invalid input"
    );

    // A duplicated genuine row: a copy of row 0 (columns [6], provenance [0]) is appended as
    // row 6. It passes construction, and selecting both copies is a kernel selection whose
    // provenance cancels to the empty set.
    let mut duplicated = parts_t(&[]);
    let copy = duplicated.matrix.rows[0].clone();
    duplicated.matrix.rows.push(copy);
    duplicated.matrix.col_weights[6] += 1;
    let context = dependency_context(&duplicated);
    assert_eq!(context.validate(&[0, 6]), Err(rejected(DependencyRejection::EmptyLift)));
    // Either copy alone is not in the kernel: column 6 has odd weight.
    assert_eq!(
        context.validate(&[6]),
        Err(rejected(DependencyRejection::NotInKernel { column: 6 }))
    );
}

#[test]
fn kat_validate_rejects_non_kernel_selection() {
    // Rows 0 and 1 are [6] and [7]: both columns have odd weight, and 6 is the smaller.
    let context = dependency_context(&parts_t(&[]));
    assert_eq!(
        context.validate(&[0, 1]),
        Err(rejected(DependencyRejection::NotInKernel { column: 6 }))
    );
    // [7] ⊕ [8] = [7, 8]: smallest odd column 7.
    assert_eq!(
        context.validate(&[1, 2]),
        Err(rejected(DependencyRejection::NotInKernel { column: 7 }))
    );
    // [6] ⊕ [1,2,6,8] = [1,2,8]: smallest odd column 1.
    assert_eq!(
        context.validate(&[0, 4]),
        Err(rejected(DependencyRejection::NotInKernel { column: 1 }))
    );
    // Membership is judged in the owned matrix, character columns included: rows 2 and 5 are both
    // [8] at K = 0 but [8,9,10] each at K = 2, so they still cancel; rows 2 and 3 do not, since
    // [8,9,10] ⊕ [7,8,9] = [7,10] has smallest odd column 7.
    let filled = dependency_context(&filled_t(&CHARS_23));
    assert!(filled.validate(&[2, 5]).is_ok());
    assert_eq!(
        filled.validate(&[2, 3]),
        Err(rejected(DependencyRejection::NotInKernel { column: 7 }))
    );
    // Rows 0 and 1 share their slot bits, so only the factor-base columns differ:
    // [6,10] ⊕ [7,10] = [6,7].
    assert_eq!(
        filled.validate(&[0, 1]),
        Err(rejected(DependencyRejection::NotInKernel { column: 6 }))
    );
}

// ─── overlap cancellation ────────────────────────────────────────────────────

/// Four genuine fixture-T relations whose matrix has a weight-three column, as `(a, b, rational,
/// algebraic, negative)`. In this order they are rows 0–3.
///
/// - (−9, 2): a − 3b = −15 = −3·5; f_hom = 243 − 18 − 4 = 221 = 13·17; labels (13,2), (17,4).
///   Columns [1, 2, 6, 8].
/// - (−3, 4): a − 3b = −15; f_hom = 27 − 12 − 16 = −1, so no algebraic entry. Columns [1, 2, 8].
/// - (−2, 1): a − 3b = −5; f_hom = 12 − 2 − 1 = 9 = 3², label (3,1) (−2 ≡ 1 mod 3). Columns
///   [2, 8].
/// - (4, 1): a − 3b = 1; f_hom = 51 = 3·17. Columns [6].
fn overlap_relations() -> Vec<Relation> {
    vec![
        relation(-9, 2, &[(1, 1), (2, 1)], &[(1, 1), (2, 1)], true),
        relation(-3, 4, &[(1, 1), (2, 1)], &[], true),
        relation(-2, 1, &[(2, 1)], &[(0, 2)], true),
        relation(4, 1, &[], &[(0, 1), (2, 1)], false),
    ]
}

#[test]
fn kat_overlap_cancellation_lifts_symmetric_difference() {
    let relations = CheckedRelations::new(&checker_t(), &overlap_relations()).unwrap();
    let mut parts = parts_for(relations, &[]);
    let rows: Vec<&[usize]> = parts.matrix.rows.iter().map(|r| r.cols.as_slice()).collect();
    assert_eq!(rows, [&[1, 2, 6, 8][..], &[1, 2, 8], &[2, 8], &[6]]);
    // Column 2 (prime 5) is set in rows 0, 1, 2: weight three. Column 1 has weight 2, column 6
    // weight 2 and the sign column 8 weight 3 (obstruction columns are never eliminated).
    assert_eq!(parts.matrix.col_weights[2], 3);

    // Eliminate column 2. The lightest incident row is row 2 (2 columns against 4 and 3), so the
    // pivot is p = row 2, and the replacements are p ⊕ row 0 then p ⊕ row 1:
    //   [2,8] ⊕ [1,2,6,8] = [1,6] with provenance [2] Δ [0] = [0,2];
    //   [2,8] ⊕ [1,2,8]   = [1]   with provenance [2] Δ [1] = [1,2].
    // Rows 0-2 are removed and row 3 moves to index 0, so the matrix is
    //   0: [6] prov [3];  1: [1,6] prov [0,2];  2: [1] prov [1,2].
    eliminate_column(&mut parts.matrix, 2).unwrap();
    let shape: Vec<(&[usize], &[usize])> =
        parts.matrix.rows.iter().map(|r| (r.cols.as_slice(), r.provenance.as_slice())).collect();
    assert_eq!(shape, [(&[6][..], &[3][..]), (&[1, 6], &[0, 2]), (&[1], &[1, 2])]);

    let context = dependency_context(&parts);
    // Rows 0, 1, 2 sum to [6] ⊕ [1,6] ⊕ [1] = [], a kernel selection that contains both
    // replacement rows. Relation 2 (the pivot) is named by both and cancels:
    //   [3] Δ [0,2] Δ [1,2] = [0,1,3]; a union would give [0,1,2,3].
    let dependency = context.validate(&[0, 1, 2]).unwrap();
    assert_eq!(dependency.selection(), [0, 1, 2]);
    assert_eq!(dependency.relations(), [0, 1, 3]);
    assert!(!dependency.relations().contains(&2));
    assert_eq!(dependency.relation_count(), 3);
    assert_eq!(
        dependency.factors().pairs(),
        [(big(-9), big(2)), (big(-3), big(4)), (big(4), big(1))]
    );

    // P over the symmetric difference: (−9 − 6)(−3 − 12)(4 − 3) = (−15)(−15)(1) = 225 = 15². The
    // union would add the pivot's a − 3b = −5 and give −1125, which is negative and not a square.
    assert_eq!(dependency.rational_product(), &big(225));
    assert_ne!(dependency.rational_product(), &big(-1125));
    // γ for the three relations, from the independent script (reduce Π(3a − bt) mod t²+t−3, then
    // divide by 3³): 950/9 + (1267/27)·θ.
    assert_eq!(dependency.gamma(), &rat(&[(950, 9), (1267, 27)]));
    assert_target_agrees(&dependency);

    // The lifted original rows do cancel: [1,2,6,8] ⊕ [1,2,8] ⊕ [6] = [].
    let sum = [0usize, 1, 3].iter().fold(Vec::<usize>::new(), |acc, &i| {
        let columns = context.original_columns(i).unwrap();
        let mut next: Vec<usize> = acc.iter().copied().filter(|c| !columns.contains(c)).collect();
        next.extend(columns.iter().copied().filter(|c| !acc.contains(c)));
        next.sort_unstable();
        next
    });
    assert!(sum.is_empty());

    // The selection that would lift to the union is not even a kernel: any proper subset of the
    // three rows fails.
    assert_eq!(
        context.validate(&[0, 1]),
        Err(rejected(DependencyRejection::NotInKernel { column: 1 }))
    );
}

// ─── exact products ──────────────────────────────────────────────────────────

#[test]
fn kat_exact_products_seam_monic() {
    // x² − 2, m = 23, N = 527 = 17·31 (f(23) = 529 − 2). FB(7, 7): rational [2,3,5,7]; ideals
    // (2,0) (7,3) (7,4). The pair (3, −2) has a − b·m = 3 + 46 = 49 = 7² and f_hom = 9 − 8 = 1,
    // so its rational entry is {7: 2} and its row is empty.
    let ctx = SqrtContext::new(&poly(&[-2, 0, 1]), &big(23), &big(527)).unwrap();
    let fb = FactorBase::new(&poly(&[-2, 0, 1]), 7, 7);
    let checker = RelationChecker::new(&ctx, &fb).unwrap();
    let relations =
        CheckedRelations::new(&checker, &[relation(3, -2, &[(3, 2)], &[], false)]).unwrap();
    let parts = parts_for(relations, &[]);
    let context = dependency_context(&parts);
    assert_eq!(context.original_columns(0), Some(&[][..]));

    let dependency = context.validate(&[0]).unwrap();
    assert_eq!(dependency.selection(), [0]);
    assert_eq!(dependency.relations(), [0]);
    assert_eq!(dependency.relation_count(), 1);
    assert_eq!(dependency.context(), &ctx);
    assert_eq!(dependency.rational_product(), &big(49));
    // c = 1, so γ = 3 − (−2)θ = 3 + 2θ with no reduction needed.
    assert_eq!(dependency.gamma(), &rat(&[(3, 1), (2, 1)]));
    assert_target_agrees(&dependency);
}

/// Norm of `x + yθ` for `θ² + θ − 3 = 0`: `x² − xy − 3y²` (sum of roots −1, product −3).
fn norm_t(x: &BigRational, y: &BigRational) -> BigRational {
    x * x - x * y - BigRational::from_integer(big(3)) * y * y
}

fn is_square(n: &BigInt) -> bool {
    n >= &big(0) && &(n.sqrt() * n.sqrt()) == n
}

#[test]
fn kat_exact_products_non_monic_without_norm_square() {
    let context = dependency_context(&parts_t(&[]));

    // Selection {2, 5}: relations (2,1) and (5,2), columns [8] ⊕ [8] = []. P = (2−3)(5−6) = 1.
    // γ = (3·2 − t)(3·5 − 2t) / 9 mod t²+t−3: (6 − t)(15 − 2t) = 90 − 27t + 2t² and 2t² = 6 − 2t,
    // so 96 − 29t; hence γ = (96 − 29θ)/9 = 32/3 − (29/9)θ.
    let dependency = context.validate(&[2, 5]).unwrap();
    assert_eq!(dependency.relations(), [2, 5]);
    assert_eq!(dependency.rational_product(), &big(1));
    assert_eq!(dependency.gamma(), &rat(&[(32, 3), (-29, 9)]));
    assert_target_agrees(&dependency);
    // Π f_hom = 13·81 = 1053 = 3⁴·13 is not a square, and Norm(γ) = (96² + 96·29 − 3·29²)/81 =
    // 9477/81 = 117 = 3²·13 is not a square either. Validation succeeds regardless: neither
    // is a precondition.
    let f_hom_product = big(13 * 81);
    assert!(!is_square(&f_hom_product));
    let coefficients = &dependency.gamma().coeffs;
    let norm = norm_t(&coefficients[0], &coefficients[1]);
    assert_eq!(norm, BigRational::from_integer(big(117)));
    assert!(!is_square(&big(117)));

    // Selection {1, 2, 3}: (7,1), (2,1), (−3,2), columns [7] ⊕ [8] ⊕ [7,8] = []. k = 3, h = 2.
    // P = (7−3)(2−3)(−3−6) = 4·(−1)·(−9) = 36. Π(3a − bt) mod t²+t−3 over (21−t)(6−t)(−9−2t):
    // (21−t)(6−t) = 126 − 27t + t² = 129 − 28t, and (129−28t)(−9−2t) = −1161 − 258t + 252t +
    // 56t² = −1161 − 6t + 56(3 − t) = −993 − 62t. So γ = (−993 − 62θ)/27.
    let dependency = context.validate(&[1, 2, 3]).unwrap();
    assert_eq!(dependency.relations(), [1, 2, 3]);
    assert_eq!(dependency.relation_count(), 3);
    assert_eq!(dependency.rational_product(), &big(36));
    assert!(is_square(dependency.rational_product()));
    assert_eq!(dependency.gamma(), &rat(&[(-331, 9), (-62, 27)]));
    assert_target_agrees(&dependency);
    // Norm(−993 − 62θ) = Π Norm(3a − bθ) = Π 3·f_hom = 27·(153·13·17) = 27·33813 = 912951, and
    // Π f_hom = 33813 = 3²·13·17² is not a square (13 has an odd exponent).
    let numerator_norm = norm_t(
        &BigRational::from_integer(big(-993)),
        &BigRational::from_integer(big(-62)),
    );
    assert_eq!(numerator_norm, BigRational::from_integer(big(912_951)));
    assert_eq!(big(27 * 33_813), big(912_951));
    assert_eq!(big(153 * 13 * 17), big(33_813));
    assert!(!is_square(&big(33_813)));
}

// ─── duplicates ──────────────────────────────────────────────────────────────

/// Keep the first occurrence of every `(a, b)`, in order. This happens before the checked relation
/// set and the matrix are built, so provenance indexes the list that was checked.
fn deduplicate(relations: Vec<Relation>) -> Vec<Relation> {
    let mut seen: Vec<(BigInt, BigInt)> = Vec::new();
    let mut kept = Vec::new();
    for r in relations {
        let key = (r.a.clone(), r.b.clone());
        if !seen.contains(&key) {
            seen.push(key);
            kept.push(r);
        }
    }
    kept
}

#[test]
fn kat_duplicates_reconciled_before_matrix() {
    let checker = checker_t();
    let r41 = || relation(4, 1, &[], &[(0, 1), (2, 1)], false);
    let r71 = || relation(7, 1, &[(0, 2)], &[(0, 2), (3, 1)], false);
    let r21 = || relation(2, 1, &[], &[(1, 1)], true);

    // The checked relation set refuses the repeated pair, naming both positions.
    let with_duplicate = vec![r41(), r71(), r41()];
    assert_eq!(
        CheckedRelations::new(&checker, &with_duplicate).unwrap_err(),
        CheckedRelationsError::DuplicatePair { first: 0, second: 2 }
    );
    let kept = deduplicate(with_duplicate);
    assert_eq!(kept.len(), 2);
    let relations = CheckedRelations::new(&checker, &kept).unwrap();
    let parts = parts_for(relations, &[]);
    let context = dependency_context(&parts);
    let provenance: Vec<&[usize]> =
        context.matrix().rows.iter().map(|r| r.provenance.as_slice()).collect();
    assert_eq!(provenance, [&[0][..], &[1]]);

    // A duplicate in the middle: after deduplication (2,1) is relation 2, and every provenance
    // index and lifted relation index refers to the deduplicated list, not the original positions.
    let kept = deduplicate(vec![r41(), r71(), r41(), r21()]);
    let relations = CheckedRelations::new(&checker, &kept).unwrap();
    let context = dependency_context(&parts_for(relations, &[]));
    assert_eq!(context.relations().relations()[2].a(), &big(2));
    assert_eq!(context.matrix().rows[2].provenance, vec![2]);
    assert_eq!(context.original_columns(2), Some(&[8][..]));
    assert_eq!(context.original_columns(3), None);

    // (a, b) and (−a, −b) are distinct records: (−4, −1) has a − 3b = −4 + 3 = −1, so a negative
    // rational sign, and f_hom(−4, −1) = 48 + 4 − 1 = 51 as for (4, 1). Labels: −4 ≡ −1 (mod 3)
    // and −4 ≡ −1·4 (mod 17). Its row is [6, 8].
    let negated = relation(-4, -1, &[], &[(0, 1), (2, 1)], true);
    let both = deduplicate(vec![r41(), r71(), negated]);
    assert_eq!(both.len(), 3, "the negated pair is not a duplicate");
    let relations = CheckedRelations::new(&checker, &both).expect("(a,b) and (−a,−b) accepted");
    let context = dependency_context(&parts_for(relations, &[]));
    assert_eq!(context.original_columns(0), Some(&[6][..]));
    assert_eq!(context.original_columns(2), Some(&[6, 8][..]));
    assert_eq!(context.matrix().rows[2].provenance, vec![2]);
}

// ─── computed kernel vectors ─────────────────────────────────────────────────

#[test]
fn kat_validate_computed_kernel_vectors() {
    // x³ − x − 1, m = 2, N = 5 (f(2) = 5), FB(30, 30), LineSieveConfig::new(10, 3). The pinned
    // relation set is the one the checked-relation tests validate; no parameter search happens.
    let f = poly(&[-1, -1, 0, 1]);
    let pair = PolyPair::new(f.clone(), poly(&[-2, 1]), big(2), big(5));
    pair.verify().expect("f(2) = 5 = N");
    let fb = FactorBase::new(&f, 30, 30);
    let collected = line_sieve(&pair, &fb, &LineSieveConfig::new(10, 3));
    let ctx = SqrtContext::from_poly_pair(&pair).unwrap();
    let checker = RelationChecker::new(&ctx, &fb).unwrap();
    let relations = CheckedRelations::new(&checker, &collected).unwrap();
    let relation_count = relations.len();

    let mut parts = parts_for(relations, &[]);
    // Singleton removal and column merging drop rows and XOR rows; neither changes the kernel of
    // the surviving rows, and both keep provenance as the XOR of original relation indices.
    parts.matrix = merge_columns(remove_singletons(parts.matrix.clone()));
    let context = dependency_context(&parts);
    assert!(context.matrix().rows.len() <= relation_count);

    let operator = MatrixOperator::new(context.matrix());
    let vectors = block_lanczos(&operator, 42);
    assert!(!vectors.is_empty(), "positive control: the pinned matrix has a nontrivial kernel");

    let mut validated = 0;
    for vector in &vectors {
        if !vector.verify(context.matrix()) {
            continue;
        }
        let dependency = context
            .validate_kernel_vector(vector)
            .unwrap_or_else(|e| panic!("kernel vector {:?} does not validate: {e}", vector));
        assert_eq!(dependency.selection(), vector.row_indices.as_slice());
        // Lifted parity: the original rows of the lifted relations sum to zero in every column,
        // recomputed here from the context's own original rows.
        let mut weights = vec![0u32; context.layout().num_cols()];
        for &i in dependency.relations() {
            for &c in context.original_columns(i).unwrap() {
                weights[c] += 1;
            }
        }
        assert!(weights.iter().all(|w| w % 2 == 0), "lifted parity for {:?}", vector.row_indices);
        // P is the product of the lifted relations' rational norms.
        let norms: BigInt = dependency
            .relations()
            .iter()
            .map(|&i| relations_norm(&context, i))
            .product();
        assert_eq!(dependency.rational_product(), &norms);
        validated += 1;
    }
    assert!(validated >= 1, "at least one computed vector was exercised");
    // The control is not vacuous: merging produced rows that stand for several relations.
    assert!(context.matrix().rows.iter().any(|r| r.provenance.len() > 1));
}

fn relations_norm(context: &DependencyContext, index: usize) -> BigInt {
    context.relations().relations()[index].rational_norm().clone()
}
