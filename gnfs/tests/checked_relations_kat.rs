//! Known-answer tests for validated relations, prime classification, and the checked matrix.
//!
//! Notation: a relation `(a, b)` has rational norm `a − b·m` and homogeneous norm
//! `f_hom(a, b) = Σ f_i a^i b^(d−i)`. For a leading coefficient `c`, the monic form is
//! `F(t) = c^(d−1) f(t/c)`, and a prime is *exceptional* when it divides `c·disc(F)`.
//!
//! Two fixtures carry most cases; every number is derived by hand next to its use.
//!
//! **Fixture T (non-monic).** `f = 3x² + x − 1`, `m = 3`, `N = 29` (`f(3) = 27 + 3 − 1 = 29`).
//! `c = 3`, `F = 3·(3(t/3)² + t/3 − 1) = t² + t − 3`, `disc(F) = 1 + 12 = 13`. `FactorBase::new(f,
//! 7, 17)` gives rational primes `[2, 3, 5, 7]` and ideals `(3,1)#0, (13,2)#1, (17,4)#2,
//! (17,7)#3`: mod 3, `f ≡ x − 1` (root 1); mod 13, `f(2) = 13` and `f'(2) = 6·2 + 1 = 13`, so 2
//! is a double root; mod 17, `f(4) = 51 = 3·17` and `f(7) = 153 = 9·17`. The primes 2, 5, 7,
//! 11 have no roots (checked by evaluating `f(0..p)`: values `−1, 3, 13, 29, 51, 79, 113, 153,
//! 199, 251, 309` contain no multiple of 5, 7, or 11, and are all odd).
//!
//! **Fixture C (monic cubic).** `f = x³ − x − 1`, `m = 2`, `N = 5` (`f(2) = 5`), `c = 1`,
//! `disc = −4·(−1)³ − 27·(−1)² = 4 − 27 = −23`. Mod 23, `x³ − x − 1 = (x − 3)(x − 10)²` because
//! `27 − 3 − 1 = 23` and `1000 − 10 − 1 = 989 = 43·23`.
//!
//! Column layout for `R` rational primes, `A` ideals, `K` character slots: `[0, R)` rational,
//! `[R, R + A)` ideals, `R + A` sign, then `K` characters.

use gnfs::filter::{
    build_checked_matrix, build_matrix, ColumnClass, ColumnCountOverflow, MatrixRow, SparseMatrix,
};
use gnfs::polyselect::PolyPair;
use gnfs::sieve::{
    classify_prime, line_sieve, AlgebraicPrime, CheckedRelations, CheckedRelationsError,
    ExponentSide, ExponentVector, FactorBase, FactorBaseRejection, LineSieveConfig, PrimeClass,
    Relation, RelationChecker, RelationRejection,
};
use gnfs::sqrt::{
    AlgebraicSqrtOutcome, CrtLimits, DirectLimits, PrecisionLimits, PrimeSearchLimits,
    RecombinationLimits, RelationFactors, SqrtContext, SqrtLimits,
};
use num_bigint::BigInt;
use shared_numfield::{IntPoly, NonresidueSearchLimits, RootSearchLimits};
use shared_numth::SmoothWitness;

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

fn context(f: &[i64], m: i64, n: i64) -> SqrtContext {
    SqrtContext::new(&poly(f), &big(m), &big(n)).expect("fixture is a valid context")
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
const F_C: [i64; 4] = [-1, -1, 0, 1];

fn context_t() -> SqrtContext {
    context(&F_T, 3, 29)
}

fn factor_base_t() -> FactorBase {
    FactorBase::new(&poly(&F_T), 7, 17)
}

fn checker_t() -> RelationChecker {
    RelationChecker::new(&context_t(), &factor_base_t()).expect("fixture T factor base is valid")
}

/// One fixture-T relation with the numbers derived by hand.
struct Case {
    a: i64,
    b: i64,
    rational_norm: i64,
    homogeneous_norm: i64,
    rational: &'static [(usize, u32)],
    algebraic: &'static [(usize, u32)],
    /// Columns of the checked row (K = 2): ideal columns are `4 + ideal index`, sign is 8.
    checked_row: &'static [usize],
    /// Columns of the legacy row (same indices, but exceptional ideals included).
    legacy_row: &'static [usize],
}

/// Fixture-T relations. Ideal indices: `(3,1)#0` (exceptional, `3 | c`), `(13,2)#1` (exceptional,
/// `13 | disc F`), `(17,4)#2` and `(17,7)#3` (good). Rational indices: 2→0, 3→1, 5→2, 7→3.
///
/// - `(4,1)`: `a − 3b = 1`; `f_hom = 48 + 4 − 1 = 51 = 3·17`. Labels: `4 ≡ 1·1 (3)`, `4 ≡ 4 (17)`.
///   Good ideal `(17,4)` has exponent 1 → column 6; `(3,1)` is exceptional, so no bit.
/// - `(7,1)`: `a − 3b = 4 = 2²`; `f_hom = 147 + 7 − 1 = 153 = 3²·17`. Labels `(3,1)`, `(17,7)`
///   (`7 ≡ 7`). Even exponents everywhere except `(17,7)` → column 7.
/// - `(2,1)`: `a − 3b = −1`; `f_hom = 12 + 2 − 1 = 13`. Label `(13,2)`. Sign column 8 only.
/// - `(−3,2)`: `a − 3b = −9 = −3²`; `f_hom = 27 − 6 − 4 = 17`. Label `−3 ≡ 2r (17)` gives
///   `r = 7` (`14 ≡ −3`). Columns 7 and 8.
/// - `(−9,2)`: `a − 3b = −15 = −3·5`; `f_hom = 243 − 18 − 4 = 221 = 13·17`. Labels: `−9 ≡ 2r (13)`
///   gives `r = 2` (`4 ≡ −9`); `−9 ≡ 2r (17)` gives `r = 4` (`8 ≡ −9`). Columns 1, 2 (3 and 5),
///   6 (good `(17,4)`), 8 (sign).
/// - `(5,2)`: `a − 3b = −1`; `f_hom = 75 + 10 − 4 = 81 = 3⁴`. Label `5 ≡ 2r (3)` gives `r = 1`.
///   Exponent 4 is even anyway. Column 8 only.
const CASES: [Case; 6] = [
    Case {
        a: 4,
        b: 1,
        rational_norm: 1,
        homogeneous_norm: 51,
        rational: &[],
        algebraic: &[(0, 1), (2, 1)],
        checked_row: &[6],
        legacy_row: &[4, 6],
    },
    Case {
        a: 7,
        b: 1,
        rational_norm: 4,
        homogeneous_norm: 153,
        rational: &[(0, 2)],
        algebraic: &[(0, 2), (3, 1)],
        checked_row: &[7],
        legacy_row: &[7],
    },
    Case {
        a: 2,
        b: 1,
        rational_norm: -1,
        homogeneous_norm: 13,
        rational: &[],
        algebraic: &[(1, 1)],
        checked_row: &[8],
        legacy_row: &[5, 8],
    },
    Case {
        a: -3,
        b: 2,
        rational_norm: -9,
        homogeneous_norm: 17,
        rational: &[(1, 2)],
        algebraic: &[(3, 1)],
        checked_row: &[7, 8],
        legacy_row: &[7, 8],
    },
    Case {
        a: -9,
        b: 2,
        rational_norm: -15,
        homogeneous_norm: 221,
        rational: &[(1, 1), (2, 1)],
        algebraic: &[(1, 1), (2, 1)],
        checked_row: &[1, 2, 6, 8],
        legacy_row: &[1, 2, 5, 6, 8],
    },
    Case {
        a: 5,
        b: 2,
        rational_norm: -1,
        homogeneous_norm: 81,
        rational: &[],
        algebraic: &[(0, 4)],
        checked_row: &[8],
        legacy_row: &[8],
    },
];

fn case_relation(case: &Case) -> Relation {
    relation(case.a, case.b, case.rational, case.algebraic, case.rational_norm < 0)
}

fn t_relations() -> Vec<Relation> {
    CASES.iter().map(case_relation).collect()
}

fn checked_t() -> CheckedRelations {
    CheckedRelations::new(&checker_t(), &t_relations()).expect("fixture T relations are valid")
}

fn checked_row(layout_rows: &SparseMatrix, i: usize) -> &[usize] {
    &layout_rows.rows[i].cols
}

// ─── classification ──────────────────────────────────────────────────────────

#[test]
fn kat_classify_primes_by_leading_coefficient_and_discriminant() {
    // T: c = 3, disc(F) = 13. 3 | c (and 3 ∤ 13) and 13 | disc are exceptional; 2, 5, 17 divide
    // neither 3 nor 13 (nor their product 39).
    let t = context_t();
    assert_eq!(t.leading_coeff(), &big(3));
    assert_eq!(t.discriminant(), &big(13));
    assert_eq!(classify_prime(&t, 3), PrimeClass::Exceptional);
    assert_eq!(classify_prime(&t, 13), PrimeClass::Exceptional);
    for good in [2, 5, 17, 29] {
        assert_eq!(classify_prime(&t, good), PrimeClass::Good, "prime {good}");
    }

    // C: disc = −23, c = 1. Both ideals above 23 (roots 3 and 10, one of them simple) lie over the
    // exceptional prime 23, so the checker classifies both as exceptional. 5 ∤ 23 is good.
    let c = context(&F_C, 2, 5);
    assert_eq!(c.discriminant(), &big(-23));
    assert_eq!(classify_prime(&c, 23), PrimeClass::Exceptional);
    assert_eq!(classify_prime(&c, 5), PrimeClass::Good);
    // 2 ∤ c·disc = −23: good for C, but exceptional for x² − 2 (disc 8, 2 | 8).
    assert_eq!(classify_prime(&c, 2), PrimeClass::Good);
    let seam = context(&[-2, 0, 1], 23, 527);
    assert_eq!(seam.discriminant(), &big(8));
    assert_eq!(classify_prime(&seam, 2), PrimeClass::Exceptional);

    // Non-primes below 2 never claim a constraint.
    assert_eq!(classify_prime(&t, 0), PrimeClass::Exceptional);
    assert_eq!(classify_prime(&t, 1), PrimeClass::Exceptional);

    let fb = FactorBase::new(&poly(&F_C), 30, 30);
    let checker = RelationChecker::new(&c, &fb).unwrap();
    let classes: Vec<_> = (0..checker.ideals().len()).map(|i| checker.ideal_class(i)).collect();
    // Ideals: (5,2) (7,5) (11,6) (17,5) (19,6) (23,3) (23,10).
    assert_eq!(checker.ideals().len(), 7);
    for (i, class) in classes.iter().enumerate() {
        let expected =
            if checker.ideals()[i].0 == 23 { PrimeClass::Exceptional } else { PrimeClass::Good };
        assert_eq!(*class, Some(expected), "ideal {i}");
    }
    assert_eq!(checker.ideal_class(7), None);
}

// ─── factor-base validation ──────────────────────────────────────────────────

fn ideal(position: usize, p: u64, r: u64) -> AlgebraicPrime {
    AlgebraicPrime { p, r, index: position, is_bad_prime: false }
}

fn base(rational: &[u64], ideals: Vec<AlgebraicPrime>) -> FactorBase {
    FactorBase {
        rational_primes: rational.to_vec(),
        algebraic_ideals: ideals,
        b_rat: 100,
        b_alg: 100,
        obstruction_count: 1,
    }
}

#[test]
fn kat_checker_rejects_malformed_factor_base() {
    let ctx = context_t();
    let reject = |fb: FactorBase| RelationChecker::new(&ctx, &fb).expect_err("must be rejected");

    // [2, 5, 3]: position 2 has 3 ≤ 5.
    assert_eq!(
        reject(base(&[2, 5, 3], vec![])),
        FactorBaseRejection::RationalNotStrictlyAscending { position: 2 }
    );
    // A repeated prime is not strictly ascending either.
    assert_eq!(
        reject(base(&[2, 3, 3], vec![])),
        FactorBaseRejection::RationalNotStrictlyAscending { position: 2 }
    );
    // 9 = 3² is composite; position 1.
    assert_eq!(
        reject(base(&[2, 9], vec![])),
        FactorBaseRejection::RationalNotPrime { position: 1, value: 9 }
    );
    // (17,7) before (17,4): same p, descending r.
    assert_eq!(
        reject(base(&[2], vec![ideal(0, 17, 7), ideal(1, 17, 4)])),
        FactorBaseRejection::IdealsNotStrictlyAscending { position: 1 }
    );
    // Index 5 at position 0.
    assert_eq!(
        reject(base(&[2], vec![AlgebraicPrime { p: 3, r: 1, index: 5, is_bad_prime: false }])),
        FactorBaseRejection::IdealIndexMismatch { position: 0, index: 5 }
    );
    // p = 15 = 3·5 is composite (and 15 is ordered, indexed, so only primality fails).
    assert_eq!(
        reject(base(&[2], vec![ideal(0, 15, 1)])),
        FactorBaseRejection::IdealNotPrime { position: 0, p: 15 }
    );
    // Root 4 is not below p = 3.
    assert_eq!(
        reject(base(&[2], vec![ideal(0, 3, 4)])),
        FactorBaseRejection::IdealRootOutOfRange { position: 0, p: 3, r: 4 }
    );
    // f(2) = 12 + 2 − 1 = 13 ≡ 1 (mod 3), so 2 is not a root mod 3.
    assert_eq!(
        reject(base(&[2], vec![ideal(0, 3, 2)])),
        FactorBaseRejection::IdealNotRoot { position: 0, p: 3, r: 2 }
    );

    // The checks follow the variant order: a rational defect is reported before an ideal defect.
    assert_eq!(
        reject(base(&[2, 9], vec![ideal(0, 3, 2)])),
        FactorBaseRejection::RationalNotPrime { position: 1, value: 9 }
    );

    // Completeness is not required: omitting (17,7) and the whole rational side is accepted.
    let partial = base(&[], vec![ideal(0, 3, 1), ideal(1, 13, 2), ideal(2, 17, 4)]);
    let checker = RelationChecker::new(&ctx, &partial).expect("incomplete base is accepted");
    assert_eq!(checker.ideals(), &[(3, 1), (13, 2), (17, 4)]);
    assert_eq!(checker.rational_primes(), &[] as &[u64]);

    // Every rejection renders a message naming its position.
    let err = FactorBaseRejection::IdealNotRoot { position: 4, p: 3, r: 2 };
    assert!(err.to_string().contains("position 4"));
    let _: &dyn std::error::Error = &err;
}

#[test]
fn kat_checker_snapshot_matches_the_factor_base() {
    let checker = checker_t();
    assert_eq!(checker.rational_primes(), &[2, 3, 5, 7]);
    assert_eq!(checker.ideals(), &[(3, 1), (13, 2), (17, 4), (17, 7)]);
    assert_eq!(checker.context(), &context_t());
    assert_eq!(
        (0..4).map(|i| checker.ideal_class(i).unwrap()).collect::<Vec<_>>(),
        vec![
            PrimeClass::Exceptional,
            PrimeClass::Exceptional,
            PrimeClass::Good,
            PrimeClass::Good
        ]
    );
}

// ─── accepted relations ──────────────────────────────────────────────────────

#[test]
fn kat_good_and_exceptional_relations_accepted() {
    let checker = checker_t();
    for case in &CASES {
        let checked = checker.check(&case_relation(case)).unwrap_or_else(|e| {
            panic!("({}, {}) must be accepted: {e}", case.a, case.b);
        });
        assert_eq!(checked.a(), &big(case.a));
        assert_eq!(checked.b(), &big(case.b));
        assert_eq!(checked.rational_norm(), &big(case.rational_norm));
        assert_eq!(checked.algebraic_norm(), &big(case.homogeneous_norm));
        assert_eq!(checked.rational_sign(), case.rational_norm < 0);
        assert_eq!(checked.rational_exponents(), case.rational);
        assert_eq!(checked.algebraic_exponents(), case.algebraic);
    }
}

#[test]
fn kat_accepted_relations_are_unit_factors_for_the_bound_context() {
    // f_hom ≠ 0 and a − b·m ≠ 0 for every relation, so every nonempty subset is accepted by the
    // reconstruction boundary of the same context. 2⁶ − 1 = 63 subsets.
    let relations = checked_t();
    let pairs: Vec<(BigInt, BigInt)> =
        relations.relations().iter().map(|r| (r.a().clone(), r.b().clone())).collect();
    for mask in 1u32..(1 << pairs.len()) {
        let subset: Vec<_> = (0..pairs.len())
            .filter(|i| mask & (1 << i) != 0)
            .map(|i| pairs[i].clone())
            .collect();
        let factors = RelationFactors::new(relations.context(), &subset)
            .unwrap_or_else(|e| panic!("subset {mask:#b} rejected: {e}"));
        assert_eq!(factors.relation_count(), subset.len());
    }
}

#[test]
fn kat_negative_b_seam_relation_accepted() {
    // x² − 2, m = 23, N = 527 = 17·31 (f(23) = 529 − 2 = 527). FB(7, 7): rational [2, 3, 5, 7];
    // ideals above 2 (root 0: f(0) = −2 ≡ 0) and 7 (roots 3, 4: 9 − 2 = 7, 16 − 2 = 14).
    let ctx = context(&[-2, 0, 1], 23, 527);
    let fb = FactorBase::new(&poly(&[-2, 0, 1]), 7, 7);
    let checker = RelationChecker::new(&ctx, &fb).unwrap();
    assert_eq!(checker.ideals(), &[(2, 0), (7, 3), (7, 4)]);

    // (3, −2): a − b·m = 3 + 46 = 49 = 7², f_hom = 9 − 8 = 1. gcd(3, 2) = 1. So the rational entry
    // is {7: 2} (index 3), there is no algebraic entry, and the sign is positive.
    let seam = relation(3, -2, &[(3, 2)], &[], false);
    let checked = checker.check(&seam).expect("b < 0 is accepted");
    assert_eq!(checked.rational_norm(), &big(49));
    assert_eq!(checked.algebraic_norm(), &big(1));
    assert_eq!(checked.rational_exponents(), &[(3, 2)]);
    assert!(checked.algebraic_exponents().is_empty());

    let relations = CheckedRelations::new(&checker, &[seam]).unwrap();
    let (matrix, layout) = build_checked_matrix(&relations, 0).unwrap();
    assert_eq!(checked_row(&matrix, 0), &[] as &[usize]);
    assert_eq!(layout.num_cols(), 4 + 3 + 1);
    assert!(matrix.col_weights.iter().all(|&w| w == 0));

    // The canonical-b form (−3, 2): a − b·m = −3 − 46 = −49, f_hom = 9 − 8 = 1. It is a different
    // record with the opposite sign (γ changes by −1), accepted and not a duplicate.
    let flipped = relation(-3, 2, &[(3, 2)], &[], true);
    let both = CheckedRelations::new(&checker, &[relation(3, -2, &[(3, 2)], &[], false), flipped])
        .expect("(3, −2) and (−3, 2) are distinct pairs");
    let (matrix, layout) = build_checked_matrix(&both, 0).unwrap();
    assert_eq!(checked_row(&matrix, 1), &[layout.sign_column()]);
}

// ─── rows ────────────────────────────────────────────────────────────────────

#[test]
fn kat_checked_rows_omit_exceptional_columns() {
    let relations = checked_t();
    let (matrix, layout) = build_checked_matrix(&relations, 2).unwrap();

    // R = 4, A = 4, K = 2: columns 0..4 rational, 4..8 ideals, 8 sign, 9 and 10 characters.
    assert_eq!(layout.num_cols(), 11);
    assert_eq!(matrix.num_cols, 11);
    assert_eq!(layout.sign_column(), 8);
    assert_eq!(matrix.obstruction_col_start, 8);
    assert_eq!(matrix.obstruction_count, 3);
    assert_eq!((layout.character_column(0), layout.character_column(1)), (Some(9), Some(10)));
    assert_eq!(matrix.rows.len(), 6);

    let mut expected_weights = vec![0u32; 11];
    for (i, case) in CASES.iter().enumerate() {
        assert_eq!(matrix.rows[i].cols, case.checked_row, "checked row ({}, {})", case.a, case.b);
        assert_eq!(matrix.rows[i].provenance, vec![i]);
        assert_eq!(layout.row_columns(&relations.relations()[i]).unwrap(), case.checked_row);
        for &c in case.checked_row {
            expected_weights[c] += 1;
        }
    }
    assert_eq!(matrix.col_weights, expected_weights);
    // Columns 4 (3,1) and 5 (13,2) are exceptional and never set.
    assert_eq!((expected_weights[4], expected_weights[5]), (0, 0));
    assert_eq!(layout.omitted_columns(), vec![4, 5]);
    // Character columns are zero at build.
    assert_eq!((expected_weights[9], expected_weights[10]), (0, 0));

    // The legacy builder differs exactly at the exceptional columns: (4,1) [4,6], (2,1) [5,8],
    // (−9,2) [1,2,5,6,8]; the other three rows agree. T's factor base has obstruction_count 1.
    let fb = factor_base_t();
    let legacy = build_matrix(&t_relations(), &fb);
    for (i, case) in CASES.iter().enumerate() {
        assert_eq!(legacy.rows[i].cols, case.legacy_row, "legacy row ({}, {})", case.a, case.b);
        let without_exceptional: Vec<usize> =
            legacy.rows[i].cols.iter().copied().filter(|c| *c != 4 && *c != 5).collect();
        assert_eq!(without_exceptional, case.checked_row);
    }
}

#[test]
fn kat_omitting_exceptional_parity_admits_a_candidate_exact_verification_rejects() {
    // (2,1) has row [8] and (5,2) has row [8], so {(2,1), (5,2)} XORs to zero in the checked
    // matrix. In the legacy matrix the rows are [5,8] and [8], which do not cancel.
    let relations = checked_t();
    let (matrix, _) = build_checked_matrix(&relations, 0).unwrap();
    assert_eq!(matrix.rows[2].cols, matrix.rows[5].cols);

    // The rational product is 1 = (−1)·(−1), a square. But Norm(a − bα) = f_hom(a, b)/c for
    // d = 2, so Norm(γ) = (13/3)·(81/3) = 1053/9 = 117 = 3²·13. A square element has a square
    // norm, and no integer squares to 117 (10² = 100 < 117 < 121 = 11²), so γ is not a square.
    let numerator = big(13) * big(81);
    assert_eq!(&numerator % big(9), big(0));
    let norm = numerator / big(9);
    assert_eq!(norm, big(117));
    assert!((0..=11).all(|k| big(k * k) != norm));

    // Exact verification agrees: the algebraic square root is certified not to exist.
    let pairs = vec![(big(2), big(1)), (big(5), big(2))];
    let factors = RelationFactors::new(relations.context(), &pairs).unwrap();
    let limits = SqrtLimits {
        precision: PrecisionLimits { max_modulus_bits: 4096 },
        direct: DirectLimits {
            prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 100 },
            nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
        },
        crt: CrtLimits {
            prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 200 },
            max_primes: 1,
            root_search: RootSearchLimits { max_evaluations: 100_000 },
            nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
        },
        recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 20 },
    };
    let outcome = gnfs::sqrt::reconstruct_root_from_factors(factors, limits).unwrap();
    assert!(matches!(outcome, AlgebraicSqrtOutcome::NotSquare(_)), "got {outcome:?}");
}

// ─── rejections ──────────────────────────────────────────────────────────────

fn rejection(checker: &RelationChecker, r: &Relation) -> RelationRejection {
    checker.check(r).expect_err("relation must be rejected")
}

#[test]
fn kat_relation_rejections_are_typed() {
    let t = checker_t();

    // NotCoprime: gcd(4, 2) = 2. Zero is coprime to nothing but ±1: gcd(0, 0) = 0.
    assert_eq!(
        rejection(&t, &relation(4, 2, &[], &[], false)),
        RelationRejection::NotCoprime { gcd: big(2) }
    );
    assert_eq!(
        rejection(&t, &relation(0, 0, &[], &[], false)),
        RelationRejection::NotCoprime { gcd: big(0) }
    );
    // Negative coordinates: gcd(−6, 4) = 2.
    assert_eq!(
        rejection(&t, &relation(-6, 4, &[], &[], true)),
        RelationRejection::NotCoprime { gcd: big(2) }
    );

    // ZeroRationalNorm: (3,1) has a − 3b = 0.
    assert_eq!(
        rejection(&t, &relation(3, 1, &[], &[], false)),
        RelationRejection::ZeroRationalNorm
    );

    // SignMismatch: (4,1) has norm 1 > 0 but claims negative; (2,1) has norm −1 but claims +.
    assert_eq!(
        rejection(&t, &relation(4, 1, &[], &[(0, 1), (2, 1)], true)),
        RelationRejection::SignMismatch { expected: false, stored: true }
    );
    assert_eq!(
        rejection(&t, &relation(2, 1, &[], &[(1, 1)], false)),
        RelationRejection::SignMismatch { expected: true, stored: false }
    );

    // MalformedExponents on the valid pair (−9,2): unsorted, duplicate, zero exponent.
    let rational = ExponentSide::Rational;
    let algebraic = ExponentSide::Algebraic;
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(2, 1), (1, 1)], &[(1, 1), (2, 1)], true)),
        RelationRejection::MalformedExponents { side: rational, position: 1 }
    );
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(1, 1), (2, 1)], &[(2, 1), (2, 1)], true)),
        RelationRejection::MalformedExponents { side: algebraic, position: 1 }
    );
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(1, 1), (2, 0)], &[(1, 1), (2, 1)], true)),
        RelationRejection::MalformedExponents { side: rational, position: 1 }
    );
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(1, 1), (2, 1)], &[(1, 0)], true)),
        RelationRejection::MalformedExponents { side: algebraic, position: 0 }
    );

    // IndexOutOfRange: sizes are R = 4 and A = 4.
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(1, 1), (4, 1)], &[(1, 1), (2, 1)], true)),
        RelationRejection::IndexOutOfRange { side: rational, index: 4, size: 4 }
    );
    assert_eq!(
        rejection(&t, &relation(-9, 2, &[(1, 1), (2, 1)], &[(1, 1), (9, 1)], true)),
        RelationRejection::IndexOutOfRange { side: algebraic, index: 9, size: 4 }
    );

    // RationalNormMismatch: (7,1) has |a − 3b| = 4 but claims 2¹ = 2; (−3,2) has 9 but claims 3³.
    assert_eq!(
        rejection(&t, &relation(7, 1, &[(0, 1)], &[(0, 2), (3, 1)], false)),
        RelationRejection::RationalNormMismatch { norm: big(4) }
    );
    assert_eq!(
        rejection(&t, &relation(-3, 2, &[(1, 3)], &[(3, 1)], true)),
        RelationRejection::RationalNormMismatch { norm: big(9) }
    );
    // An absurd exponent is refused without forming 2^u32::MAX.
    assert_eq!(
        rejection(&t, &relation(7, 1, &[(0, u32::MAX)], &[(0, 2), (3, 1)], false)),
        RelationRejection::RationalNormMismatch { norm: big(4) }
    );

    // AlgebraicNormMismatch: (4,1) has |f_hom| = 51 = 3·17 but lists only (3,1) (label 4 ≡ 1 (3)
    // holds), so the product 3 ≠ 51.
    assert_eq!(
        rejection(&t, &relation(4, 1, &[], &[(0, 1)], false)),
        RelationRejection::AlgebraicNormMismatch { norm: big(51) }
    );

    // AffineLabelMismatch: (17,7) is claimed on (4,1), but 4 − 7·1 = −3 is not divisible by 17.
    assert_eq!(
        rejection(&t, &relation(4, 1, &[], &[(0, 1), (3, 1)], false)),
        RelationRejection::AffineLabelMismatch { position: 1, p: 17, r: 7 }
    );
    // An exceptional entry's label is checked too: (3,1) on (2,1) needs 2 ≡ 1·1 (mod 3), false.
    assert_eq!(
        rejection(&t, &relation(2, 1, &[], &[(0, 1)], true)),
        RelationRejection::AffineLabelMismatch { position: 0, p: 3, r: 1 }
    );

    // ZeroAlgebraicNorm: x² − 1, m = 4, N = 15 (f(4) = 15). FB(5, 5) has rational [2, 3, 5]. The
    // pair (1,1): a − 4b = −3 = −3¹ (index 1), but f_hom = 1 − 1 = 0.
    let ctx = context(&[-1, 0, 1], 4, 15);
    let fb = FactorBase::new(&poly(&[-1, 0, 1]), 5, 5);
    let reducible = RelationChecker::new(&ctx, &fb).unwrap();
    assert_eq!(
        rejection(&reducible, &relation(1, 1, &[(1, 1)], &[], true)),
        RelationRejection::ZeroAlgebraicNorm
    );
}

#[test]
fn kat_projective_collection_limit_is_rejected_not_repaired() {
    // Under T, (1,3): a − 3b = 1 − 9 = −8 = −2³, f_hom = 3 + 3 − 9 = −3. Here 3 | c, 3 | b and
    // 3 | f_hom: the factor is projective at 3. No affine root r gives 1 ≡ 3r (mod 3).
    let t = checker_t();
    let rational = [(0usize, 3u32)];
    // Claiming the affine ideal (3,1) fails its label: 1 − 3·1 = −2.
    assert_eq!(
        rejection(&t, &relation(1, 3, &rational, &[(0, 1)], true)),
        RelationRejection::AffineLabelMismatch { position: 0, p: 3, r: 1 }
    );
    // Listing nothing leaves the factor 3 of |f_hom| = 3 unexplained.
    assert_eq!(
        rejection(&t, &relation(1, 3, &rational, &[], true)),
        RelationRejection::AlgebraicNormMismatch { norm: big(3) }
    );

    // The collector's constructor refuses it too: no r with r·0 ≡ 1 (mod 3).
    let witness = |factors: Vec<(u64, u32)>| SmoothWitness {
        factors,
        cofactor: crypto_bigint::Uint::<4>::ONE,
    };
    let built = Relation::new(
        big(1),
        big(3),
        &witness(vec![(2, 3)]),
        &witness(vec![(3, 1)]),
        true,
        &factor_base_t(),
    );
    assert!(built.is_none());
}

#[test]
fn kat_rejection_messages_and_error_sources() {
    let rejection = RelationRejection::AffineLabelMismatch { position: 2, p: 17, r: 7 };
    let text = rejection.to_string();
    assert!(text.contains("17") && text.contains('7'), "{text}");
    let wrapped = CheckedRelationsError::Relation { index: 3, rejection };
    assert!(wrapped.to_string().starts_with("relation 3"));
    assert!(std::error::Error::source(&wrapped).is_some());
    let duplicate = CheckedRelationsError::DuplicatePair { first: 0, second: 2 };
    assert!(std::error::Error::source(&duplicate).is_none());
    assert!(ColumnCountOverflow.to_string().contains("overflow"));
}

// ─── duplicate pairs ─────────────────────────────────────────────────────────

#[test]
fn kat_duplicate_pairs_rejected() {
    let checker = checker_t();
    let r = |i: usize| case_relation(&CASES[i]);

    // [(4,1), (7,1), (4,1)]: relation 2 repeats relation 0.
    let err = CheckedRelations::new(&checker, &[r(0), r(1), r(0)]).unwrap_err();
    assert_eq!(err, CheckedRelationsError::DuplicatePair { first: 0, second: 2 });

    // A bad relation before the duplicate wins; a bad relation is reported with its own index.
    let bad = relation(4, 2, &[], &[], false);
    let err = CheckedRelations::new(&checker, &[r(0), bad.clone(), r(0)]).unwrap_err();
    assert_eq!(
        err,
        CheckedRelationsError::Relation {
            index: 1,
            rejection: RelationRejection::NotCoprime { gcd: big(2) }
        }
    );
    // A repeated pair whose record is itself invalid (its algebraic list is empty, but
    // |f_hom| = 51) is reported as invalid, not as a duplicate: the check precedes the comparison.
    let repeat = relation(4, 1, &[], &[], false);
    let err = CheckedRelations::new(&checker, &[r(0), repeat]).unwrap_err();
    assert!(matches!(err, CheckedRelationsError::Relation { index: 1, .. }));

    // (a, b) and (−a, −b) are distinct: (−4, −1) has a − 3b = −4 + 3 = −1, f_hom = 48 + 4 − 1 = 51
    // (even degree), labels (3,1) (−4 ≡ −1 mod 3 since 2 ≡ 2) and (17,4) (−4 ≡ −4).
    let negated = relation(-4, -1, &[], &[(0, 1), (2, 1)], true);
    let both = CheckedRelations::new(&checker, &[r(0), negated]).expect("distinct pairs");
    assert_eq!(both.len(), 2);
    assert_eq!(both.relations()[1].rational_norm(), &big(-1));
    assert_eq!(both.relations()[1].b(), &big(-1));
    // Their rows differ exactly in the sign column.
    let (matrix, layout) = build_checked_matrix(&both, 0).unwrap();
    assert_eq!(matrix.rows[0].cols, vec![6]);
    assert_eq!(matrix.rows[1].cols, vec![6, layout.sign_column()]);

    // An empty list is accepted and yields a matrix of the right width with no rows.
    let empty = CheckedRelations::new(&checker, &[]).unwrap();
    assert!(empty.is_empty());
    let (matrix, layout) = build_checked_matrix(&empty, 1).unwrap();
    assert!(matrix.rows.is_empty());
    assert_eq!(matrix.num_cols, 4 + 4 + 1 + 1);
    assert!(layout.shape_matches(&matrix));
}

// ─── line sieve control ──────────────────────────────────────────────────────

#[test]
fn kat_line_sieve_relations_validate() {
    // Fixture C with FB(30, 30): rational primes 2..29 (R = 10), ideals (5,2) (7,5) (11,6) (17,5)
    // (19,6) (23,3) (23,10) (A = 7). The two ideals above 23 sit at columns 10 + 5 and 10 + 6.
    let f = poly(&F_C);
    let g = poly(&[-2, 1]);
    let pair = PolyPair::new(f.clone(), g, big(2), big(5));
    pair.verify().expect("x³ − x − 1 with m = 2 has f(2) = 5 = N");
    let fb = FactorBase::new(&f, 30, 30);
    let collected = line_sieve(&pair, &fb, &LineSieveConfig::new(10, 3));
    assert!(collected.len() >= 5, "positive control: the region contains smooth pairs");

    let ctx = SqrtContext::from_poly_pair(&pair).unwrap();
    let checker = RelationChecker::new(&ctx, &fb).unwrap();
    for r in &collected {
        checker.check(r).unwrap_or_else(|e| panic!("({}, {}) fails: {e}", r.a, r.b));
    }
    let relations =
        CheckedRelations::new(&checker, &collected).expect("collector gives unique pairs");
    assert_eq!(relations.len(), collected.len());

    // With K = 0 and the factor base's single obstruction column, the legacy width is
    // R + A + 1 = 18, which equals the checked width.
    let (matrix, layout) = build_checked_matrix(&relations, 0).unwrap();
    let legacy = build_matrix(&collected, &fb);
    assert_eq!(layout.num_cols(), 18);
    assert!(layout.shape_matches(&legacy));
    assert_eq!(layout.omitted_columns(), vec![15, 16]);
    let mut differing = 0;
    for (i, (checked_row, legacy_row)) in matrix.rows.iter().zip(&legacy.rows).enumerate() {
        let filtered: Vec<usize> =
            legacy_row.cols.iter().copied().filter(|c| *c != 15 && *c != 16).collect();
        assert_eq!(checked_row.cols, filtered, "row {i}");
        assert_eq!(checked_row.provenance, legacy_row.provenance);
        differing += usize::from(checked_row.cols != legacy_row.cols);
    }
    // The control is not vacuous: (3,1) has algebraic norm 23, so column 15 is set in the legacy
    // matrix and removed from the checked one.
    assert!(differing >= 1);
    assert_eq!(matrix.col_weights[15] + matrix.col_weights[16], 0);
    assert!(legacy.col_weights[15] + legacy.col_weights[16] > 0);
}

// ─── layout ──────────────────────────────────────────────────────────────────

#[test]
fn kat_layout_accessors_and_binding() {
    let relations = checked_t();
    let (matrix, layout) = build_checked_matrix(&relations, 2).unwrap();

    assert_eq!(layout.checker(), relations.checker());
    assert_eq!(relations.context(), relations.checker().context());
    assert_eq!(layout.checker().context(), &context_t());
    assert_eq!(
        (layout.rational_size(), layout.algebraic_size(), layout.character_count()),
        (4, 4, 2)
    );
    assert_eq!(layout.obstruction_col_start(), 8);
    assert_eq!(layout.obstruction_count(), 3);
    assert_eq!(layout.num_cols(), 11);

    // class() for every column kind, and None past the end.
    assert_eq!(layout.class(0), Some(ColumnClass::RationalPrime { prime: 2 }));
    assert_eq!(layout.class(3), Some(ColumnClass::RationalPrime { prime: 7 }));
    assert_eq!(layout.class(4), Some(ColumnClass::ExceptionalIdeal { prime: 3, root: 1 }));
    assert_eq!(layout.class(5), Some(ColumnClass::ExceptionalIdeal { prime: 13, root: 2 }));
    assert_eq!(layout.class(6), Some(ColumnClass::GoodIdeal { prime: 17, root: 4 }));
    assert_eq!(layout.class(7), Some(ColumnClass::GoodIdeal { prime: 17, root: 7 }));
    assert_eq!(layout.class(8), Some(ColumnClass::RationalSign));
    assert_eq!(layout.class(9), Some(ColumnClass::RootCharacter { slot: 0 }));
    assert_eq!(layout.class(10), Some(ColumnClass::RootCharacter { slot: 1 }));
    assert_eq!(layout.class(11), None);
    assert_eq!(layout.character_column(2), None);
    assert_eq!(layout.character_column(usize::MAX), None);

    // shape_matches: true for the built matrix, false for each perturbed field.
    assert!(layout.shape_matches(&matrix));
    let perturbed = |edit: &dyn Fn(&mut SparseMatrix)| {
        let mut m = matrix.clone();
        edit(&mut m);
        layout.shape_matches(&m)
    };
    assert!(!perturbed(&|m| m.num_cols += 1));
    assert!(!perturbed(&|m| m.obstruction_col_start += 1));
    assert!(!perturbed(&|m| m.obstruction_count = 1));
    assert!(!perturbed(&|m| m.col_weights.push(0)));
    // Rows and provenance are not examined.
    assert!(perturbed(&|m| m.rows.push(MatrixRow { cols: vec![0], provenance: vec![99] })));
    // A layout with K = 0 has the legacy shape of the same factor base, but not K = 2's.
    let (_, narrow) = build_checked_matrix(&relations, 0).unwrap();
    assert!(narrow.shape_matches(&build_matrix(&t_relations(), &factor_base_t())));
    assert!(!narrow.shape_matches(&matrix));

    // Bound by equality: a checker from a larger factor base differs from this one.
    let larger = RelationChecker::new(&context_t(), &FactorBase::new(&poly(&F_T), 11, 17)).unwrap();
    assert_ne!(&larger, relations.checker());
    // (1,4): a − 3b = −11 (rational prime 11 at index 4, outside T's R = 4); f_hom = 3 + 4 − 16
    // = −9 = −3² with label 1 ≡ 4r (mod 3), r = 1.
    let foreign = larger
        .check(&relation(1, 4, &[(4, 1)], &[(0, 2)], true))
        .expect("accepted by the larger checker");
    assert_eq!(layout.row_columns(&foreign), None);
    // Its own layout places the same relation: 11 is column 4, and the sign is column 5 + 4 = 9.
    let (_, wide) = build_checked_matrix(&CheckedRelations::new(&larger, &[]).unwrap(), 0).unwrap();
    assert_eq!(wide.row_columns(&foreign), Some(vec![4, wide.sign_column()]));
    assert_eq!(wide.sign_column(), 9);

    // Overflow of R + A + 1 + K is an error, with no allocation attempted.
    assert_eq!(build_checked_matrix(&relations, usize::MAX).unwrap_err(), ColumnCountOverflow);
    // R + A = 8, so K = usize::MAX − 8 gives exactly usize::MAX + 1.
    assert_eq!(build_checked_matrix(&relations, usize::MAX - 8).unwrap_err(), ColumnCountOverflow);
}

// ─── legacy path ─────────────────────────────────────────────────────────────

#[test]
fn kat_legacy_matrix_path_is_unchanged() {
    // The public record and matrix shapes are built by literal here, so any added or removed field
    // is a compile error. Legacy interpretation: every ideal column is a parity constraint.
    let f = poly(&F_T);
    let pair = PolyPair::new(f.clone(), poly(&[-3, 1]), big(3), big(29));
    pair.verify().expect("f(3) = 29");
    let fb = factor_base_t();
    let relations = t_relations();
    for r in &relations {
        r.verify(&pair, &fb).expect("legacy verification accepts all six");
    }

    let legacy = build_matrix(&relations, &fb);
    assert_eq!(legacy.num_cols, 9);
    assert_eq!(legacy.obstruction_col_start, 8);
    assert_eq!(legacy.obstruction_count, 1);
    let mut weights = vec![0u32; 9];
    for (i, case) in CASES.iter().enumerate() {
        assert_eq!(legacy.rows[i].cols, case.legacy_row);
        assert_eq!(legacy.rows[i].provenance, vec![i]);
        for &c in case.legacy_row {
            weights[c] += 1;
        }
    }
    assert_eq!(legacy.col_weights, weights);

    // Rebuilding from the public fields reproduces the record exactly.
    let rebuilt = SparseMatrix {
        rows: legacy.rows.clone(),
        num_cols: legacy.num_cols,
        obstruction_col_start: legacy.obstruction_col_start,
        obstruction_count: legacy.obstruction_count,
        col_weights: legacy.col_weights.clone(),
    };
    assert_eq!(rebuilt.rows, legacy.rows);
    assert_eq!(rebuilt.col_weights, legacy.col_weights);
}
