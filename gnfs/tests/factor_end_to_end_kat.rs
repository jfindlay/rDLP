//! Known-answer tests for congruence assembly and the checked factor driver.
//!
//! A dependency gives an exact rational root `x` with `x² = P` and a verified algebraic root `β`
//! whose image `y = β(c·m) mod N` has `y² ≡ P (mod N)`. Then `x² ≡ y² (mod N)`, and `x ≢ ±y`
//! gives the proper divisor `gcd(x − y, N)`. `attempt_factor` runs one dependency; `factor`
//! retries a bounded candidate list and reports every attempt.
//!
//! **Fixtures** (every relation set is built through the checked pipeline; the matrix and
//! provenance are genuine, and each candidate is a validated kernel selection):
//!
//! - **S** (seam): `x² − 2`, `m = 23`, `N = 527 = 17·31` (`f(23) = 527`), relation (3,−2): a
//!   nonconstant `b`. `a − 23b = 49 = 7²` so `x = 7`; `γ = 3 + 2α = (1 + α)²` so `β = ±(1 + α)`
//!   and `y = ±24`. `gcd(7 − 24, 527) = gcd(17, 527) = 17`, and `7 + 24 = 31`.
//! - **R** (eight relations, `x² − 2`, `m = 23`, `N = 527`, `FB(13, 13)`): rows 0–7 are (10,−1),
//!   (−3,2), (−2,−1), (−2,1), (4,3), (2,−1), (10,1), (1,0), with `a − 23b` equal to 33, −49, 21,
//!   −25, −65, 25, −13, 1. Its kernel is spanned by {1,3}, {5}, {7}. The `γ` are
//!   {1,3}: `10 + 7α` (not a square: odd valuation at `√2`); {5}: `2 + α` (same); {7}: `1`;
//!   {1,3,5}: `34 + 24α = (4 + 3α)²` with `P = (−49)(−25)(25) = 175²`, so `x = 175`, `y = ±73`
//!   and `gcd(175 − 73, 527) = gcd(102, 527) = 17`. Also {1,3,5,7} gives the same square.
//!   The exceptional prime above 2 carries no matrix column, which is why the nonsquare vectors
//!   {1,3} and {5} are kernel vectors at all.
//! - **NM** (non-monic): `f = 2x² − 2x − 5`, `m = 5`, `N = 35`, relations (9,−11), (−1,−1):
//!   `γ = 37/2 + 9α = (4 + α)²`, `P = 64·4 = 256`, `x = 16`, `y = ±9`, `gcd(16 − 9, 35) = 7`.
//! - **IC**: `x² − 12`, `m = 164`, `N = 286 = 2·11·13` (`f(164) = 26884 = 94·286`), `FB(100, 100)`.
//!   Rows 0–2 are (148,−7), (1,0), (961,−276), each a single-row kernel vector: row 0 has
//!   `P = 148 + 1148 = 1296 = 36²` and `γ = 148 + 7α = (1 + 7α/2)²`, a root with denominator 2,
//!   and `gcd(2, 286) = 2` is an incidental divisor; row 1 is the unit with `x = y = 1`; row 2 has
//!   `P = 961 + 164·276 = 46225 = 215²`, `γ = 961 + 276α = (23 + 6α)²`, `x = 215`, `y = ±149`,
//!   and `gcd(215 − 149, 286) = gcd(66, 286) = 22` (or 26 for `−149`).
//!
//! Hand values were recomputed with an independent exact-arithmetic script (Python `fractions`,
//! no Rust), and the arithmetic is written beside each use. Pinned limits keep
//! `crt.max_primes = 1` with ample precision headroom.

use gnfs::filter::build_checked_matrix;
use gnfs::linalg::{
    KernelSearchLimits, KernelSearchResource, KernelSearchStatus, KernelVector, MatrixOperator,
    RootCharacters, block_lanczos, search_kernel_combinations,
};
use gnfs::sieve::{CheckedRelations, ExponentVector, FactorBase, Relation, RelationChecker};
use gnfs::sqrt::{
    AlgebraicSqrtOutcome, AttemptOutcome, CheckedDependency, CongruenceError, CongruenceOutcome,
    CrtLimits, DependencyContext, DependencyError, DependencyRejection, DirectLimits, FactorError,
    FactorOutcome, FailureClass, IncidentalFactor, IncidentalSource, NonSquareEvidence,
    PrecisionLimits, PrimeSearchLimits, RecombinationLimits, ResourceKind, SqrtContext,
    SqrtLimits, UnsupportedDomain, VerifiedFactor, algebraic_sqrt, attempt_factor, factor,
    factor_from_congruence,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, Zero};
use shared_numfield::{IntPoly, NonresidueSearchLimits, RatPoly, RootSearchLimits};

// ─── helpers ─────────────────────────────────────────────────────────────────

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

fn rat(coeffs: &[(i64, i64)]) -> RatPoly {
    RatPoly::from_coeffs(coeffs.iter().map(|&(n, d)| BigRational::new(big(n), big(d))).collect())
}

fn negated(p: &RatPoly) -> RatPoly {
    RatPoly::from_coeffs(p.coeffs.iter().map(|c| -c.clone()).collect())
}

fn strip(n: &mut BigInt, p: u64) -> u32 {
    let p = BigInt::from(p);
    let mut e = 0;
    while !n.is_zero() && (&*n % &p).is_zero() {
        *n /= &p;
        e += 1;
    }
    e
}

/// The relation `(a, b)` with exponent vectors read off its two norms by trial division over the
/// factor base; the checker re-verifies every norm, so a wrong factorization cannot pass.
fn relation_of(fb: &FactorBase, f: &IntPoly, m: i64, a: i64, b: i64) -> Relation {
    let (a_big, b_big) = (big(a), big(b));
    let mut rational_norm = &a_big - &b_big * big(m);
    let rational_sign = rational_norm.is_negative();
    rational_norm = rational_norm.abs();
    let rational: Vec<(usize, u32)> = fb
        .rational_primes
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| {
            let e = strip(&mut rational_norm, p);
            (e > 0).then_some((i, e))
        })
        .collect();
    assert_eq!(rational_norm, big(1), "the rational norm of ({a},{b}) is smooth over the base");

    let d = f.coeffs.len() - 1;
    let mut algebraic_norm = BigInt::zero();
    for (i, c) in f.coeffs.iter().enumerate() {
        algebraic_norm += c * a_big.pow(i as u32) * b_big.pow((d - i) as u32);
    }
    algebraic_norm = algebraic_norm.abs();
    let algebraic: Vec<(usize, u32)> = fb
        .algebraic_ideals
        .iter()
        .filter_map(|ideal| {
            let label = (&a_big - &b_big * big(ideal.r as i64)) % BigInt::from(ideal.p);
            if !label.is_zero() {
                return None;
            }
            let e = strip(&mut algebraic_norm, ideal.p);
            (e > 0).then_some((ideal.index, e))
        })
        .collect();
    assert_eq!(algebraic_norm, big(1), "the algebraic norm of ({a},{b}) is smooth over the base");

    Relation {
        a: a_big,
        b: b_big,
        rational_exponents: ExponentVector { entries: rational },
        algebraic_exponents: ExponentVector { entries: algebraic },
        rational_sign,
    }
}

fn dependency_context(
    f: &[i64],
    m: i64,
    n: i64,
    bounds: (u64, u64),
    pairs: &[(i64, i64)],
) -> DependencyContext {
    let f = poly(f);
    let ctx = SqrtContext::new(&f, &big(m), &big(n)).expect("valid context");
    let fb = FactorBase::new(&f, bounds.0, bounds.1);
    let checker = RelationChecker::new(&ctx, &fb).expect("valid factor base");
    let relations: Vec<Relation> =
        pairs.iter().map(|&(a, b)| relation_of(&fb, &f, m, a, b)).collect();
    let checked = CheckedRelations::new(&checker, &relations).expect("valid relations");
    let (matrix, layout) = build_checked_matrix(&checked, 0).expect("column count fits");
    let characters = RootCharacters::new(checked.checker(), &[]).expect("no characters");
    DependencyContext::new(checked, matrix, layout, characters).expect("genuine parts")
}

fn context_s() -> DependencyContext {
    dependency_context(&[-2, 0, 1], 23, 527, (13, 13), &[(3, -2)])
}

fn context_r() -> DependencyContext {
    dependency_context(
        &[-2, 0, 1],
        23,
        527,
        (13, 13),
        &[(10, -1), (-3, 2), (-2, -1), (-2, 1), (4, 3), (2, -1), (10, 1), (1, 0)],
    )
}

fn context_nm() -> DependencyContext {
    dependency_context(&[-5, -2, 2], 5, 35, (13, 13), &[(9, -11), (-1, -1)])
}

fn context_ic() -> DependencyContext {
    dependency_context(&[-12, 0, 1], 164, 286, (100, 100), &[(148, -7), (1, 0), (961, -276)])
}

fn kv(rows: &[usize]) -> KernelVector {
    KernelVector { row_indices: rows.to_vec() }
}

fn kvs(vectors: &[&[usize]]) -> Vec<KernelVector> {
    vectors.iter().map(|rows| kv(rows)).collect()
}

fn ample() -> SqrtLimits {
    SqrtLimits {
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
    }
}

const SEARCH_AMPLE: KernelSearchLimits =
    KernelSearchLimits { max_xor_operations: 1 << 20, max_candidates: 1 << 20 };

/// The divisor and cofactor are one of the listed splits of `n`, in either order. The sign of the
/// algebraic root is not fixed, so `y` and `−y` select different gcds.
fn assert_split(found: &VerifiedFactor, n: i64, splits: &[[i64; 2]]) {
    let pair = (found.divisor().clone(), found.cofactor().clone());
    assert!(
        splits.iter().any(|s| pair == (big(s[0]), big(s[1])) || pair == (big(s[1]), big(s[0]))),
        "{pair:?} is not among the splits {splits:?} of {n}"
    );
    assert_eq!(found.divisor() * found.cofactor(), big(n));
}

/// The congruence behind a verified factor holds, with `x` exact and `y` verified.
fn assert_congruence(found: &VerifiedFactor, n: i64) {
    let x = found.rational_root().modular_image();
    let y = found.algebraic_root().modular_image();
    assert_eq!((x * x - y * y) % big(n), big(0), "x² ≡ y² (mod {n})");
    assert_eq!(found.rational_root().root() % big(n), x.clone());
    assert_eq!(found.rational_root().product(), found.algebraic_root().rational_product());
}

fn classes(report: &gnfs::FactorReport) -> Vec<Option<FailureClass>> {
    report.attempts().iter().map(|a| a.outcome.class()).collect()
}

fn found_factor(report: &gnfs::FactorReport) -> &VerifiedFactor {
    match report.outcome() {
        FactorOutcome::Factor(found) => found,
        other => panic!("expected a factor, got {other:?}"),
    }
}

// ─── congruence assembly ─────────────────────────────────────────────────────

#[test]
fn kat_congruence_extracts_divisor_from_difference() {
    // N = 35: 6² = 36 ≡ 1 = 1², and gcd(6 − 1, 35) = gcd(5, 35) = 5, cofactor 7.
    assert_eq!(
        factor_from_congruence(&big(6), &big(1), &big(35)),
        Ok(CongruenceOutcome::Proper { divisor: big(5), cofactor: big(7) })
    );
    // 29² = 841 = 24·35 + 1 ≡ 1²; 29 ≡ −1 (mod 5) but 29 ≡ 1 (mod 7), so 29 ≢ ±1 (mod 35) and
    // gcd(|29 − 1|, 35) = gcd(28, 35) = 7, cofactor 5.
    assert_eq!(
        factor_from_congruence(&big(29), &big(1), &big(35)),
        Ok(CongruenceOutcome::Proper { divisor: big(7), cofactor: big(5) })
    );
    // Signs and magnitudes are free: (−6)² = 36 ≡ 1, and gcd(|−6 − 1|, 35) = gcd(7, 35) = 7.
    assert_eq!(
        factor_from_congruence(&big(-6), &big(1), &big(35)),
        Ok(CongruenceOutcome::Proper { divisor: big(7), cofactor: big(5) })
    );
}

#[test]
fn kat_congruence_trivial_when_x_is_plus_or_minus_y() {
    let n = big(35);
    // x = y: 6 − 6 = 0, gcd = 35.
    assert_eq!(factor_from_congruence(&big(6), &big(6), &n), Ok(CongruenceOutcome::Trivial));
    // x ≡ −y: 6 + 29 = 35, so 6 ≡ −29, and 29² = 841 = 24·35 + 1 ≡ 6² (36 ≡ 1).
    assert_eq!(factor_from_congruence(&big(6), &big(29), &n), Ok(CongruenceOutcome::Trivial));
    // x ≡ y with a common factor of x and N: 5 ≡ 5 and gcd(5 + 5, 35) = 5, but that divisor says
    // nothing about the congruence, so it is trivial. Likewise for an even modulus: x = y = 1
    // mod 22 would otherwise report gcd(2, 22) = 2.
    assert_eq!(factor_from_congruence(&big(5), &big(5), &n), Ok(CongruenceOutcome::Trivial));
    assert_eq!(factor_from_congruence(&big(1), &big(1), &big(22)), Ok(CongruenceOutcome::Trivial));
}

#[test]
fn kat_congruence_rejects_bad_input_without_panic() {
    for n in [-35, -1, 0, 1] {
        assert_eq!(
            factor_from_congruence(&big(6), &big(1), &big(n)),
            Err(CongruenceError::ModulusNotGreaterThanOne { n: big(n) }),
            "N = {n}"
        );
    }
    // 2² = 4 and 3² = 9 differ by 5, which 35 does not divide.
    assert_eq!(
        factor_from_congruence(&big(2), &big(3), &big(35)),
        Err(CongruenceError::NotCongruent)
    );
    assert!(CongruenceError::NotCongruent.to_string().contains("congruent"));
    assert!(
        CongruenceError::ModulusNotGreaterThanOne { n: big(0) }.to_string().contains("got 0")
    );
}

// ─── one dependency, whole path ──────────────────────────────────────────────

#[test]
fn kat_seam_fixture_gives_a_verified_factor_through_the_whole_path() {
    let context = context_s();
    let dependency = context.validate(&[0]).unwrap();
    assert_eq!(dependency.rational_product(), &big(49));
    let found = match attempt_factor(&dependency, ample()).unwrap() {
        AttemptOutcome::Factor(found) => found,
        other => panic!("expected a factor, got {other:?}"),
    };
    // x = 7 (7² = 49) and y = ±24 (24² = 576 = 527 + 49): 7 − 24 = −17 and 7 + 24 = 31, so the
    // divisors are 17 and 31.
    assert_eq!(found.rational_root().root(), &big(7));
    assert_eq!(found.rational_root().modular_image(), &big(7));
    assert!([big(24), big(503)].contains(found.algebraic_root().modular_image()));
    assert_split(&found, 527, &[[17, 31]]);
    assert_congruence(&found, 527);
    assert_eq!(found.relations(), [0]);
    assert!(
        found.algebraic_root().beta() == &rat(&[(1, 1), (1, 1)])
            || found.algebraic_root().beta() == &negated(&rat(&[(1, 1), (1, 1)]))
    );

    // The same dependency through `factor`.
    let report = factor(&context, &kvs(&[&[0]]), KernelSearchStatus::Complete, ample()).unwrap();
    assert_eq!(found_factor(&report), &found);
    assert_eq!(report.attempts().len(), 1);
    assert_eq!(report.attempts()[0].candidate, 0);
    assert_eq!(report.attempts()[0].relations, vec![0]);
    assert!(!report.resource_limited());
}

#[test]
fn kat_non_monic_fixture_gives_a_verified_factor() {
    // c = 2: β = ±(4 + θ/2) with θ = 2α, x = 16, y = ±9, and the divisor is 7 (y = 9) or 5
    // (y = 26: gcd(16 − 26, 35) = gcd(10, 35) = 5).
    let context = context_nm();
    let dependency = context.validate(&[0, 1]).unwrap();
    let found = match attempt_factor(&dependency, ample()).unwrap() {
        AttemptOutcome::Factor(found) => found,
        other => panic!("expected a factor, got {other:?}"),
    };
    assert_eq!(found.rational_root().root(), &big(16));
    assert!([big(9), big(26)].contains(found.algebraic_root().modular_image()));
    assert_split(&found, 35, &[[5, 7]]);
    assert_congruence(&found, 35);
    assert_eq!(found.relations(), [0, 1]);
    let beta = rat(&[(4, 1), (1, 2)]);
    assert!(found.algebraic_root().beta() == &beta || found.algebraic_root().beta() == &negated(&beta));
}

// ─── retry over candidates ───────────────────────────────────────────────────

#[test]
fn kat_retry_passes_nonsquares_to_a_later_success() {
    let context = context_r();
    // {5} and {1,3} are kernel vectors whose algebraic product is not a square; {1,3,5} is.
    let candidates = kvs(&[&[5], &[1, 3], &[1, 3, 5]]);
    let report = factor(&context, &candidates, KernelSearchStatus::Complete, ample()).unwrap();
    let found = found_factor(&report);
    assert_split(found, 527, &[[17, 31]]);
    assert_congruence(found, 527);
    assert_eq!(found.relations(), [1, 3, 5]);
    assert_eq!(found.rational_root().root(), &big(175));

    assert_eq!(
        classes(&report),
        vec![
            Some(FailureClass::AlgebraicNonSquare),
            Some(FailureClass::AlgebraicNonSquare),
            None
        ]
    );
    let positions: Vec<usize> = report.attempts().iter().map(|a| a.candidate).collect();
    assert_eq!(positions, [0, 1, 2]);
    let relations: Vec<&[usize]> = report.attempts().iter().map(|a| a.relations.as_slice()).collect();
    assert_eq!(relations, [&[5][..], &[1, 3], &[1, 3, 5]]);
    for attempt in &report.attempts()[..2] {
        assert!(matches!(attempt.outcome, AttemptOutcome::AlgebraicNonSquare(_)));
    }
    assert_eq!(report.search(), &KernelSearchStatus::Complete);
    assert!(!report.resource_limited());
}

#[test]
fn kat_factor_stops_at_the_first_factor() {
    // {1,3,5} and {1,3,5,7} give the same square; only the first is attempted.
    let report = factor(
        &context_r(),
        &kvs(&[&[1, 3, 5], &[1, 3, 5, 7]]),
        KernelSearchStatus::Complete,
        ample(),
    )
    .unwrap();
    assert_eq!(report.attempts().len(), 1);
    assert_eq!(found_factor(&report).relations(), [1, 3, 5]);
}

#[test]
fn kat_factor_found_only_in_a_combination_of_computed_vectors() {
    let context = context_r();
    // Kernel basis {5}, {1,3}, {7}: none of these three factors ({5}, {1,3} are nonsquares, {7}
    // is the trivial unit), but their sum {1,3,5} does. The search emits the singles first, then
    // pairs in order: {1,3,5} = {5} ⊕ {1,3} is the fourth candidate.
    let basis = kvs(&[&[5], &[1, 3], &[7]]);
    for vector in &basis {
        assert!(vector.verify(context.matrix()));
    }
    let search = search_kernel_combinations(&basis, SEARCH_AMPLE);
    assert_eq!(search.status, KernelSearchStatus::Complete);
    assert_eq!(search.rank, 3);
    assert_eq!(search.candidates.len(), 7);
    assert_eq!(search.candidates[3], kv(&[1, 3, 5]));
    assert!(!basis.contains(&search.candidates[3]), "the winner is a genuine sum");

    let report =
        factor(&context, &search.candidates, search.status.clone(), ample()).unwrap();
    assert_eq!(
        classes(&report),
        vec![
            Some(FailureClass::AlgebraicNonSquare),
            Some(FailureClass::AlgebraicNonSquare),
            Some(FailureClass::TrivialCongruence),
            None,
        ]
    );
    let found = found_factor(&report);
    assert_eq!(found.relations(), [1, 3, 5]);
    assert_split(found, 527, &[[17, 31]]);
    assert_eq!(report.attempts()[3].candidate, 3);
    assert!(!report.resource_limited());
}

#[test]
fn kat_factor_from_lanczos_vectors_and_their_combinations() {
    // The kernel of R is exactly {1,3}, {5}, {7} and their sums (rank 3, seven nonzero vectors).
    // The computed vectors come from block Lanczos with a fixed seed; whichever basis it returns,
    // the complete combination set is the whole kernel, and exactly {1,3,5} and {1,3,5,7} are
    // factor-bearing.
    let context = context_r();
    let operator = MatrixOperator::new(context.matrix());
    let computed: Vec<KernelVector> =
        block_lanczos(&operator, 42).into_iter().filter(|v| v.verify(context.matrix())).collect();
    let search = search_kernel_combinations(&computed, SEARCH_AMPLE);
    assert_eq!(search.status, KernelSearchStatus::Complete);
    assert_eq!(search.rank, 3, "Lanczos spans the kernel of the pinned matrix");
    assert_eq!(search.candidates.len(), 7);

    let report = factor(&context, &search.candidates, search.status.clone(), ample()).unwrap();
    let found = found_factor(&report);
    assert!([vec![1, 3, 5], vec![1, 3, 5, 7]].contains(&found.relations().to_vec()));
    assert_split(found, 527, &[[17, 31]]);
    // Every attempt before the factor is a typed non-success.
    let last = report.attempts().len() - 1;
    for attempt in &report.attempts()[..last] {
        assert!(attempt.outcome.class().is_some());
    }
}

// ─── trivial congruences, incidental divisors ────────────────────────────────

#[test]
fn kat_trivial_congruence_is_recorded_and_exhausts() {
    // R row 7 is (1,0): P = 1, β = 1, x = y = 1, so the congruence holds trivially.
    let context = context_r();
    let report =
        factor(&context, &kvs(&[&[7]]), KernelSearchStatus::Complete, ample()).unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert!(!report.resource_limited(), "a complete negative result for the supplied candidates");
    // y = β(23) = ±1, that is 1 or 526 modulo 527; either sign is trivial.
    assert!(
        [big(1), big(526)].iter().any(|y| {
            report.attempts()[0].outcome == AttemptOutcome::TrivialCongruence { x: big(1), y: y.clone() }
        }),
        "{:?}",
        report.attempts()[0].outcome
    );
    assert_eq!(report.attempts()[0].outcome.class(), Some(FailureClass::TrivialCongruence));

    // Exhausting a candidate list of nonsquares and trivial vectors is the same.
    let report = factor(
        &context,
        &kvs(&[&[5], &[1, 3], &[7]]),
        KernelSearchStatus::Complete,
        ample(),
    )
    .unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert_eq!(report.attempts().len(), 3);
    assert!(!report.resource_limited());
}

#[test]
fn kat_incidental_divisor_has_its_own_provenance() {
    let context = context_ic();
    // Row 0: β = 1 + 7α/2 has the denominator 2, and gcd(2, 286) = 2 is a proper divisor that
    // comes from the map, not from a congruence of squares.
    let dependency = context.validate(&[0]).unwrap();
    assert_eq!(dependency.rational_product(), &big(1296));
    match attempt_factor(&dependency, ample()).unwrap() {
        AttemptOutcome::Incidental(incidental) => {
            assert_eq!((incidental.divisor(), incidental.cofactor()), (&big(2), &big(143)));
            assert_eq!(incidental.source(), IncidentalSource::MapDenominator);
        }
        other => panic!("expected an incidental divisor, got {other:?}"),
    }

    // Incidental then trivial: the report names the incidental divisor but is not a factor.
    let report = factor(
        &context,
        &kvs(&[&[0], &[1]]),
        KernelSearchStatus::Complete,
        ample(),
    )
    .unwrap();
    match report.outcome() {
        FactorOutcome::IncidentalOnly(incidental) => {
            assert_eq!(incidental.divisor(), &big(2));
            assert_eq!(incidental.source(), IncidentalSource::MapDenominator);
        }
        other => panic!("expected incidental-only, got {other:?}"),
    }
    assert_eq!(
        classes(&report),
        vec![Some(FailureClass::IncidentalFactor), Some(FailureClass::TrivialCongruence)]
    );
    assert!(!report.resource_limited());

    // Incidental, trivial, then a verified congruence factor: the incidental divisor is recorded
    // and the loop continues; the outcome is the verified factor, never the incidental one. Row 2
    // has x = 215, y = ±149: 22·13 = 286 (gcd(66, 286) = 22) or 26·11 (y = 137: gcd(78, 286) = 26).
    let report = factor(
        &context,
        &kvs(&[&[0], &[1], &[2]]),
        KernelSearchStatus::Complete,
        ample(),
    )
    .unwrap();
    let found = found_factor(&report);
    assert_eq!(found.relations(), [2]);
    assert_split(found, 286, &[[22, 13], [26, 11]]);
    assert_congruence(found, 286);
    assert_eq!(
        classes(&report),
        vec![Some(FailureClass::IncidentalFactor), Some(FailureClass::TrivialCongruence), None]
    );
    assert!(matches!(report.attempts()[0].outcome, AttemptOutcome::Incidental(_)));
}

// ─── resource limits ─────────────────────────────────────────────────────────

#[test]
fn kat_exhausted_distinguishes_complete_from_resource_limited() {
    let context = context_r();
    let basis = kvs(&[&[5], &[1, 3], &[7]]);

    // A complete search over supplied candidates with no factor: a complete negative result for
    // those candidates, not resource limited.
    let report = factor(&context, &basis, KernelSearchStatus::Complete, ample()).unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert!(!report.resource_limited());

    // A search cut short: the real enumeration refuses its fourth candidate (the first pair), so
    // only the three singles are supplied, and the report says the search was not complete.
    let limits = KernelSearchLimits { max_xor_operations: 100, max_candidates: 3 };
    let search = search_kernel_combinations(&basis, limits);
    assert_eq!(search.candidates, basis);
    assert!(matches!(
        search.status,
        KernelSearchStatus::ResourceLimited { resource: KernelSearchResource::Candidates, .. }
    ));
    let report = factor(&context, &search.candidates, search.status.clone(), ample()).unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert_eq!(report.search(), &search.status);
    assert!(report.resource_limited(), "the missing combinations include the factor");
    assert_eq!(report.attempts().len(), 3);
}

#[test]
fn kat_attempt_resource_exhaustion_is_never_a_nonsquare() {
    let context = context_r();
    // One bit of modulus precision cannot hold the bound for any target: the reconstruction
    // reports a precision exhaustion for a dependency that does factor.
    let mut tiny = ample();
    tiny.precision = PrecisionLimits { max_modulus_bits: 1 };
    let dependency = context.validate(&[1, 3, 5]).unwrap();
    match attempt_factor(&dependency, tiny).unwrap() {
        AttemptOutcome::ResourceExhausted(spent) => {
            assert_eq!(spent.len(), 1);
            assert_eq!(spent[0].kind, ResourceKind::Precision);
            assert_eq!(spent[0].limit, 1);
        }
        other => panic!("expected resource exhaustion, got {other:?}"),
    }

    // Through `factor` with a search that was complete: the attempt alone makes the report
    // resource limited, and the dependency that would have factored is not recorded as a
    // nonsquare.
    let report = factor(&context, &kvs(&[&[1, 3, 5]]), KernelSearchStatus::Complete, tiny).unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert!(report.resource_limited());
    assert_eq!(classes(&report), vec![Some(FailureClass::ResourceExhausted)]);
    assert!(!matches!(report.attempts()[0].outcome, AttemptOutcome::AlgebraicNonSquare(_)));

    // The same dependency under the ample limits factors, so the exhaustion was the limit.
    let report = factor(&context, &kvs(&[&[1, 3, 5]]), KernelSearchStatus::Complete, ample())
        .unwrap();
    assert!(matches!(report.outcome(), FactorOutcome::Factor(_)));
}

// ─── malformed candidates ────────────────────────────────────────────────────

#[test]
fn kat_invalid_candidate_aborts_instead_of_retrying() {
    let context = context_r();
    let aborted = |candidates: &[KernelVector]| {
        factor(&context, candidates, KernelSearchStatus::Complete, ample())
            .expect_err("an invalid candidate aborts")
    };

    // Position 1 names row 9 of an eight-row matrix: the loop aborts there, and the factor-bearing
    // candidate after it is never tried.
    let error = aborted(&kvs(&[&[5], &[9], &[1, 3, 5]]));
    assert_eq!(
        error,
        FactorError::InvalidCandidate {
            candidate: 1,
            error: DependencyError::Rejected(DependencyRejection::RowOutOfRange {
                position: 0,
                row: 9,
                rows: 8
            }),
        }
    );
    assert!(error.to_string().contains("candidate 1"));
    assert!(std::error::Error::source(&error).is_some());

    // Row 0 has columns [1, 4]: it is not a kernel vector, and the first odd column is 1.
    assert_eq!(
        aborted(&kvs(&[&[0]])),
        FactorError::InvalidCandidate {
            candidate: 0,
            error: DependencyError::Rejected(DependencyRejection::NotInKernel { column: 1 }),
        }
    );
    // Empty and unsorted selections.
    assert!(matches!(
        aborted(&kvs(&[&[7], &[]])),
        FactorError::InvalidCandidate {
            candidate: 1,
            error: DependencyError::Rejected(DependencyRejection::EmptySelection)
        }
    ));
    assert!(matches!(
        aborted(&kvs(&[&[3, 1]])),
        FactorError::InvalidCandidate {
            candidate: 0,
            error: DependencyError::Rejected(DependencyRejection::MalformedSelection {
                position: 1
            })
        }
    ));
}

#[test]
fn kat_no_candidates_is_an_empty_report() {
    let context = context_r();
    let report = factor(&context, &[], KernelSearchStatus::Complete, ample()).unwrap();
    assert_eq!(report.outcome(), &FactorOutcome::Exhausted);
    assert!(report.attempts().is_empty());
    assert!(!report.resource_limited());
}

// ─── typed outcomes ──────────────────────────────────────────────────────────

#[test]
fn kat_attempt_outcome_classes() {
    let context = context_r();
    let dependency = |rows: &[usize]| -> CheckedDependency { context.validate(rows).unwrap() };

    let nonsquare = attempt_factor(&dependency(&[5]), ample()).unwrap();
    assert_eq!(nonsquare.class(), Some(FailureClass::AlgebraicNonSquare));
    // The certificate is the one the reconstruction returns, passed through unchanged.
    let AttemptOutcome::AlgebraicNonSquare(certificate) = &nonsquare else { panic!("{nonsquare:?}") };
    let AlgebraicSqrtOutcome::NotSquare(from_reconstruction) =
        algebraic_sqrt(&dependency(&[5]), ample()).unwrap()
    else {
        panic!("the reconstruction certifies the nonsquare")
    };
    assert_eq!(certificate, &from_reconstruction);
    assert!(matches!(
        certificate.evidence(),
        NonSquareEvidence::InertLocalNonresidue { .. }
            | NonSquareEvidence::DirectBoundExhausted { .. }
    ));

    let factored = attempt_factor(&dependency(&[1, 3, 5]), ample()).unwrap();
    assert!(factored.class().is_none());
    let trivial = attempt_factor(&dependency(&[7]), ample()).unwrap();
    assert_eq!(trivial.class(), Some(FailureClass::TrivialCongruence));

    // Variants that no validated dependency produces here are still classified.
    let incidental = IncidentalFactor::new(&big(286), big(2), IncidentalSource::MapDenominator)
        .expect("2 is a proper divisor of 286");
    for (outcome, class) in [
        (
            AttemptOutcome::RationalNonSquare { product: big(-1) },
            FailureClass::RationalNonSquare,
        ),
        (
            AttemptOutcome::Unsupported(UnsupportedDomain::MapDenominatorNotUnit),
            FailureClass::Unsupported,
        ),
        (AttemptOutcome::Incidental(incidental), FailureClass::IncidentalFactor),
        (AttemptOutcome::ResourceExhausted(vec![]), FailureClass::ResourceExhausted),
    ] {
        assert_eq!(outcome.class(), Some(class));
    }
}

// ─── oracle ──────────────────────────────────────────────────────────────────

/// Factor a published 80–100-bit challenge using the full GNFS pipeline.
///
/// This KAT requires CADO-NFS or msieve as an external oracle to verify the result. It is
/// ignored by default; the deterministic fixtures above carry reproducibility without the oracle.
///
/// To run manually:
/// ```
/// cargo test --test factor_end_to_end_kat kat_c_oracle -- --ignored
/// ```
#[test]
#[ignore = "CADO-NFS not installed; run manually when available to verify end-to-end \
            factorisation of 80-100 bit challenge"]
fn kat_c_oracle_80_100_bit_challenge() {
    // Placeholder: factor a published 80–100-bit semiprime using the full GNFS pipeline.
    // The pipeline stages (polyselect → sieve → filter → linalg → sqrt) must all be
    // wired end-to-end. When CADO-NFS is available, verify the recovered factor against its
    // output.
    todo!(
        "Oracle KAT: wire the full GNFS pipeline (polyselect → sieve → filter → linalg → sqrt) \
         and factor an 80–100-bit semiprime. Verify against CADO-NFS / msieve."
    );
}
