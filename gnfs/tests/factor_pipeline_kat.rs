//! Sieve-to-factor demonstration: relations produced by `line_sieve` are validated, reduced,
//! solved for kernel vectors, lifted, square-rooted exactly, and turned into a proper factor.
//!
//! Every stage is the library's own: `line_sieve`, exact relation validation and deduplication,
//! the checked matrix, singleton removal / clique pruning / column merging, eight root
//! characters, `block_lanczos`, `search_kernel_combinations`, `DependencyContext`, and `factor`
//! (verified rational root, verified algebraic root, congruence, gcd). No relation, dependency,
//! kernel vector, or factor is supplied by hand, and no brute-force or incidental divisor counts.
//!
//! **Fixture.** `N = 45113 = 197·229`, `m = 31`, `f = x³ + 15x² + 29x + 8` (the textbook worked
//! example of the number field sieve; its relations and dependency are not used here). Then
//! `f(31) = 29791 + 14415 + 899 + 8 = 45113` exactly, and `f` is irreducible: modulo 3 it is
//! `x³ + 2x + 2`, which has no root (the values at 0, 1, 2 are 2, 5, 14, all nonzero mod 3), and
//! a cubic without a root is irreducible. An independent Python (sympy) check agrees, as does
//! the crate's own Rabin test in `kat_fixture_is_exact_and_irreducible`. The discriminant is
//! `44581 = 109·409`, so 109 is an exceptional prime and gets no parity or character column.
//! Sieve region: `a ∈ [−400, 400]`, `b ∈ [1, 41]` (32 841 lattice points) over factor bases
//! `B_rat = 30`, `B_alg = 90`, default threshold scale.
//!
//! **Pinned options** (fixed before this configuration was run and not changed afterwards; this
//! is the first of a fixed list of polynomial and sieve-size configurations that factored under
//! them):
//!
//! - Primary reconstruction: `DirectLimits { prime_search: { first_candidate: 0,
//!   max_primes_examined: 100 }, nonresidue_search: { max_candidates: 1000 } }`.
//! - Fallback reconstruction: `CrtLimits { prime_search: { first_candidate: 0,
//!   max_primes_examined: 200 }, max_primes: 1, root_search: { max_evaluations: 100_000 },
//!   nonresidue_search: { max_candidates: 1000 } }`.
//! - Precision: `max_modulus_bits = 65536`, a cap on `bits(2B + 1)` far above what these
//!   dependencies need; recombination `{ max_list_entries: 1 << 20, max_matches: 1 << 20 }`.
//! - Root characters: eight pairs from `select_root_characters` under `RootCharacterLimits {
//!   prime_search: { first_candidate: 3, max_primes_examined: 1000 }, root_search: {
//!   max_evaluations: 100_000 } }`, filled after the filter by XOR over complete provenance.
//! - Filter: `remove_singletons`, then `prune_cliques_reported` (the library's excess floor),
//!   then `merge_columns`.
//! - Solver: `block_lanczos` with the fixed seeds 42, 43, 44. Every returned vector is checked
//!   against the filtered matrix and the verified vectors are pooled.
//! - Candidates: `KernelSearchLimits { max_xor_operations: 1 << 20, max_candidates: 4096 }`: at
//!   most 4096 distinct nonzero combinations of the pooled vectors, with a merge budget far above
//!   what reduction and enumeration of that many combinations need.
//!
//! The first candidate factors, so the run attempts one dependency. The exact expected record
//! (60 relations, a 60 × 39 matrix filtered to 47 rows, rank 18, factors 197 and 229) is a
//! known answer for this deterministic fixture. It is a tested scope for one input, not a claim
//! about all inputs or about production-scale performance.

use std::collections::{BTreeMap, HashSet};

use gnfs::filter::merge::prune_cliques_reported;
use gnfs::filter::{build_checked_matrix, merge_columns, remove_singletons};
use gnfs::linalg::{
    KernelSearch, KernelSearchLimits, KernelSearchStatus, KernelVector, MatrixOperator,
    RootCharacterLimits, block_lanczos, search_kernel_combinations, select_root_characters,
};
use gnfs::polyselect::PolyPair;
use gnfs::sieve::{
    CheckedRelations, FactorBase, LineSieveConfig, Relation, RelationChecker, line_sieve,
};
use gnfs::sqrt::{
    AttemptOutcome, CrtLimits, DependencyContext, DirectLimits, FactorOutcome, FactorReport,
    PrecisionLimits, PrimeSearchLimits, RecombinationLimits, ReconstructionMethod, SqrtContext,
    SqrtLimits, factor,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use shared_numfield::{
    FpPoly, IntPoly, NonresidueSearchLimits, PrimeModulus, RatPoly, RootSearchLimits, is_prime_u64,
};

// ─── pinned options ──────────────────────────────────────────────────────────

const CHARACTER_COUNT: usize = 8;
const SOLVER_SEEDS: [u64; 3] = [42, 43, 44];
const SEARCH_LIMITS: KernelSearchLimits = KernelSearchLimits {
    max_xor_operations: 1 << 20,
    max_candidates: 4096,
};
const B_RAT: u64 = 30;
const B_ALG: u64 = 90;
const A_BOUND: u64 = 400;
const B_MAX: u64 = 41;

fn sqrt_limits() -> SqrtLimits {
    SqrtLimits {
        precision: PrecisionLimits {
            max_modulus_bits: 65_536,
        },
        direct: DirectLimits {
            prime_search: PrimeSearchLimits {
                first_candidate: 0,
                max_primes_examined: 100,
            },
            nonresidue_search: NonresidueSearchLimits {
                max_candidates: 1000,
            },
        },
        crt: CrtLimits {
            prime_search: PrimeSearchLimits {
                first_candidate: 0,
                max_primes_examined: 200,
            },
            max_primes: 1,
            root_search: RootSearchLimits {
                max_evaluations: 100_000,
            },
            nonresidue_search: NonresidueSearchLimits {
                max_candidates: 1000,
            },
        },
        recombination: RecombinationLimits {
            max_list_entries: 1 << 20,
            max_matches: 1 << 20,
        },
    }
}

fn character_limits() -> RootCharacterLimits {
    RootCharacterLimits {
        prime_search: PrimeSearchLimits {
            first_candidate: 3,
            max_primes_examined: 1000,
        },
        root_search: RootSearchLimits {
            max_evaluations: 100_000,
        },
    }
}

// ─── fixture and helpers ─────────────────────────────────────────────────────

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

fn fixture() -> PolyPair {
    PolyPair::new(poly(&[8, 29, 15, 1]), poly(&[-31, 1]), big(31), big(45113))
}

/// A prime `p` not dividing the leading coefficient at which `f` keeps its degree and is
/// irreducible (Rabin's test). Such a prime certifies that `f` is irreducible over `Q`: a
/// factorization over `Z` into positive degrees would reduce to one modulo `p`.
fn irreducibility_certificate(f: &IntPoly) -> Option<u64> {
    let degree = f.degree()?;
    let lead = f.leading_coeff()?;
    (3u64..1000).filter(|&p| is_prime_u64(p)).find(|&p| {
        let modulus = PrimeModulus::new(p).expect("odd prime modulus");
        let reduced = FpPoly::from_int_poly(f, modulus);
        (lead % BigInt::from(p)) != BigInt::zero()
            && reduced.degree() == Some(degree)
            && reduced.is_irreducible()
    })
}

// ─── the pipeline ────────────────────────────────────────────────────────────

/// What a run did, stage by stage: plain data for the assertions and for failure messages.
#[derive(Debug, Default)]
struct Record {
    lattice_points: u64,
    sieved: usize,
    duplicates: usize,
    validation_failures: Vec<String>,
    relations: usize,
    matrix_rows_initial: usize,
    matrix_cols: usize,
    omitted_columns: usize,
    rows_after_singletons: usize,
    pruned_rows_removed: usize,
    excess_before_pruning: isize,
    excess_after_pruning: isize,
    rows_after_pruning: usize,
    rows_after_merge: usize,
    characters: Vec<(u64, u64)>,
    solver_vectors_returned: Vec<usize>,
    solver_vectors_verified: usize,
    kernel_rank: usize,
    combinations: usize,
    search_status: String,
    attempts: usize,
    attempt_classes: BTreeMap<String, usize>,
}

/// The objects a run leaves for inspection.
struct Run {
    record: Record,
    context: DependencyContext,
    search: KernelSearch,
    report: FactorReport,
}

/// Run the pinned pipeline end to end. A library error at a stage that cannot fail on valid input
/// panics with `DEFECT`, so a regression is never mistaken for an ordinary negative result.
fn run_pipeline() -> Run {
    let pair = fixture();
    pair.verify()
        .expect("DEFECT: the polynomial pair is invalid");
    let ctx = SqrtContext::from_poly_pair(&pair).expect("DEFECT: the pair is not a valid context");
    let mut record = Record {
        lattice_points: (2 * A_BOUND + 1) * B_MAX,
        ..Record::default()
    };

    // Sieve.
    let fb = FactorBase::new(&pair.f, B_RAT, B_ALG);
    let sieved = line_sieve(&pair, &fb, &LineSieveConfig::new(A_BOUND, B_MAX));
    record.sieved = sieved.len();

    // Exact validation and deduplication, keeping first occurrences in order.
    let checker = RelationChecker::new(&ctx, &fb).expect("DEFECT: invalid factor base");
    let mut seen: HashSet<(BigInt, BigInt)> = HashSet::new();
    let mut unique: Vec<Relation> = Vec::new();
    for relation in sieved {
        if !seen.insert((relation.a.clone(), relation.b.clone())) {
            record.duplicates += 1;
            continue;
        }
        match checker.check(&relation) {
            Ok(_) => unique.push(relation),
            Err(rejection) => record
                .validation_failures
                .push(format!("({}, {}): {rejection:?}", relation.a, relation.b)),
        }
    }
    assert!(
        record.validation_failures.is_empty(),
        "DEFECT: the sieve produced relations the checker rejects: {record:#?}"
    );
    let checked = CheckedRelations::new(&checker, &unique).expect("DEFECT: deduplicated set");
    record.relations = checked.len();

    // Matrix, then singleton removal, clique pruning, and column merging.
    let (matrix, layout) =
        build_checked_matrix(&checked, CHARACTER_COUNT).expect("DEFECT: column count overflow");
    record.matrix_rows_initial = matrix.rows.len();
    record.matrix_cols = layout.num_cols();
    record.omitted_columns = layout.omitted_columns().len();
    let matrix = remove_singletons(matrix);
    record.rows_after_singletons = matrix.rows.len();
    let (matrix, pruning) = prune_cliques_reported(matrix);
    record.pruned_rows_removed = pruning.removed.len();
    record.excess_before_pruning = pruning.excess_before;
    record.excess_after_pruning = pruning.excess_after;
    record.rows_after_pruning = matrix.rows.len();
    let mut matrix = merge_columns(matrix);
    record.rows_after_merge = matrix.rows.len();

    // Root characters, filled by XOR over complete provenance.
    let selection = select_root_characters(&checked, CHARACTER_COUNT, character_limits())
        .expect("eight valid root characters exist within the pinned search");
    let (characters, _report) = selection.into_parts();
    record.characters = characters
        .characters()
        .iter()
        .map(|c| (c.prime(), c.root()))
        .collect();
    characters
        .fill_columns(&mut matrix, &layout, &checked)
        .expect("DEFECT: characters could not be filled");
    let context = DependencyContext::new(checked, matrix, layout, characters)
        .expect("DEFECT: the reduced matrix is not bound to its relations");

    // Kernel vectors from the existing solver, then bounded combinations.
    let operator = MatrixOperator::new(context.matrix());
    let mut computed: Vec<KernelVector> = Vec::new();
    for seed in SOLVER_SEEDS {
        let vectors = block_lanczos(&operator, seed);
        record.solver_vectors_returned.push(vectors.len());
        computed.extend(vectors.into_iter().filter(|v| v.verify(context.matrix())));
    }
    record.solver_vectors_verified = computed.len();
    let search = search_kernel_combinations(&computed, SEARCH_LIMITS);
    record.kernel_rank = search.rank;
    record.combinations = search.candidates.len();
    record.search_status = format!("{:?}", search.status);

    // Verified roots and gcd.
    let report = factor(
        &context,
        &search.candidates,
        search.status.clone(),
        sqrt_limits(),
    )
    .unwrap_or_else(|e| panic!("DEFECT: factor aborted: {e:?}\n{record:#?}"));
    record.attempts = report.attempts().len();
    for attempt in &report.attempts()[..] {
        let class = attempt
            .outcome
            .class()
            .map_or_else(|| "Factor".to_string(), |c| format!("{c:?}"));
        *record.attempt_classes.entry(class).or_default() += 1;
    }
    Run {
        record,
        context,
        search,
        report,
    }
}

// ─── independent re-derivation ──────────────────────────────────────────────

/// `β` re-expressed in the `α` basis: the `θ`-basis coefficient of `θ^i` times `c^i`.
fn in_alpha_basis(beta: &RatPoly, c: &BigInt) -> RatPoly {
    let mut power = BigInt::one();
    let mut coeffs = Vec::with_capacity(beta.coeffs.len());
    for coeff in &beta.coeffs {
        coeffs.push(coeff * BigRational::from_integer(power.clone()));
        power *= c;
    }
    RatPoly::from_coeffs(coeffs)
}

/// `Π(a − b·α) mod f`, from the raw pairs.
fn product_in_number_field(pairs: &[(BigInt, BigInt)], f: &IntPoly) -> RatPoly {
    let modulus = f.to_rat_poly();
    pairs.iter().fold(RatPoly::one(), |acc, (a, b)| {
        let factor = RatPoly::from_coeffs(vec![
            BigRational::from_integer(a.clone()),
            BigRational::from_integer(-b.clone()),
        ]);
        acc.mul(&factor).rem(&modulus)
    })
}

/// `β(c·m) mod N` for a `θ`-basis `β`, with every denominator a unit modulo `N`.
fn image_modulo_n(beta: &RatPoly, c: &BigInt, m: &BigInt, n: &BigInt) -> BigInt {
    let point = BigRational::from_integer(c * m);
    let mut value = BigRational::zero();
    let mut power = BigRational::one();
    for coeff in &beta.coeffs {
        value += coeff * &power;
        power *= &point;
    }
    let gcd = value.denom().gcd(n);
    assert!(
        gcd.is_one(),
        "the denominator of β(c·m) shares {gcd} with N"
    );
    let inverse = value.denom().extended_gcd(n).x.mod_floor(n);
    (value.numer() * inverse).mod_floor(n)
}

// ─── tests ───────────────────────────────────────────────────────────────────

#[test]
fn kat_fixture_is_exact_and_irreducible() {
    let pair = fixture();
    // f(31) = 31³ + 15·31² + 29·31 + 8 = 29791 + 14415 + 899 + 8 = 45113, exactly.
    assert_eq!(pair.f.eval(&pair.m), pair.n);
    assert_eq!(pair.n, big(197 * 229));
    pair.verify().expect("f(m) ≡ 0 (mod N) and g = x − m");
    assert_eq!(pair.degree, 3);
    assert!(pair.f.is_monic());
    // x³ + 2x + 2 over F_3 has no root and is irreducible, so f is irreducible over Q.
    assert_eq!(irreducibility_certificate(&pair.f), Some(3));
    // The checked context accepts the fixture and exposes the exact m, N and degree.
    let ctx = SqrtContext::from_poly_pair(&pair).expect("valid context");
    assert_eq!(ctx.degree(), 3);
    assert_eq!(ctx.root(), &big(31));
    assert_eq!(ctx.modulus(), &big(45113));
}

#[test]
fn kat_pipeline_factors_from_generated_relations() {
    let run = run_pipeline();
    let Run {
        record,
        context,
        search,
        report,
    } = &run;

    // The generated relations: 60 of 32 841 sieve points, all valid, none duplicated.
    assert_eq!(record.lattice_points, 32_841);
    assert_eq!(record.sieved, 60, "{record:#?}");
    assert_eq!(
        (record.duplicates, record.validation_failures.len()),
        (0, 0)
    );
    assert_eq!(record.relations, 60);

    // The filter: 60 rows over 39 columns. The exceptional prime 109 exceeds B_alg = 90, so it has
    // no ideal and no column is omitted. No singleton; ten rows are pruned down to the excess
    // floor of 20, and column merging leaves 47 rows.
    assert_eq!(
        (
            record.matrix_rows_initial,
            record.matrix_cols,
            record.omitted_columns
        ),
        (60, 39, 0)
    );
    assert_eq!(record.rows_after_singletons, 60);
    assert_eq!(
        (
            record.pruned_rows_removed,
            record.excess_before_pruning,
            record.excess_after_pruning
        ),
        (10, 30, 20)
    );
    assert_eq!(
        (record.rows_after_pruning, record.rows_after_merge),
        (50, 47)
    );

    // Eight characters. Every root at an odd prime below 97 vanishes on some relation (a ≡ b·s
    // mod q) and is rejected, and 109 is skipped as exceptional.
    assert_eq!(
        record.characters,
        [
            (97, 28),
            (101, 87),
            (103, 47),
            (107, 4),
            (107, 8),
            (107, 80),
            (113, 108),
            (139, 33)
        ]
    );
    assert_eq!(context.characters().len(), CHARACTER_COUNT);
    assert_eq!(context.matrix().rows.len(), 47);

    // The computed kernel: every vector the solver returned that is a kernel vector was pooled,
    // their span has rank 18, and the combination search stops at its 4096-candidate bound.
    assert_eq!(record.solver_vectors_returned, [12, 12, 20]);
    assert_eq!(record.solver_vectors_verified, 44);
    assert_eq!(record.kernel_rank, 18);
    assert_eq!(record.combinations, 4096);
    assert!(search.candidates.len() as u64 <= SEARCH_LIMITS.max_candidates);
    assert!(
        matches!(search.status, KernelSearchStatus::ResourceLimited { .. }),
        "{record:#?}"
    );
    let distinct: HashSet<&[usize]> = search
        .candidates
        .iter()
        .map(|v| v.row_indices.as_slice())
        .collect();
    assert_eq!(
        distinct.len(),
        search.candidates.len(),
        "candidates are distinct"
    );
    assert!(search.candidates.iter().all(|v| !v.row_indices.is_empty()));

    // The checked path: one attempt, the first candidate, and it is a verified factor.
    let FactorOutcome::Factor(found) = report.outcome() else {
        panic!(
            "expected a verified factor, got {:?}\n{record:#?}",
            report.outcome()
        );
    };
    assert_eq!(record.attempts, 1);
    assert_eq!(
        record.attempt_classes,
        BTreeMap::from([("Factor".to_string(), 1)])
    );
    let last = report.attempts().last().expect("a factor needs an attempt");
    assert_eq!(last.candidate, 0);
    assert!(matches!(&last.outcome, AttemptOutcome::Factor(f) if f == found));
    assert!(
        report.attempts()[..report.attempts().len() - 1]
            .iter()
            .all(|a| a.outcome.class().is_some())
    );

    // Provenance: the candidate is a verified kernel vector of the filtered matrix that lifts,
    // by symmetric difference of provenance, to exactly the factor's 22 relations.
    let candidate = &search.candidates[last.candidate];
    assert!(candidate.verify(context.matrix()));
    assert_eq!(last.relations, found.relations());
    assert_eq!(
        context
            .matrix()
            .lift_dependency(&candidate.row_indices)
            .as_deref(),
        Some(found.relations())
    );
    assert_eq!(
        found.relations(),
        [
            2, 3, 9, 12, 16, 17, 19, 20, 21, 23, 26, 30, 33, 38, 44, 45, 47, 51, 53, 55, 56, 58
        ]
    );
    assert!(context.validate_kernel_vector(candidate).is_ok());

    // Proper factor: 45113 = 197·229, in either order, checked by multiplication.
    let (divisor, cofactor) = (found.divisor(), found.cofactor());
    let mut split = [divisor.clone(), cofactor.clone()];
    split.sort();
    assert_eq!(split, [big(197), big(229)]);
    assert_eq!(divisor * cofactor, big(45113));
    assert!(divisor > &BigInt::one() && cofactor > &BigInt::one());

    // Reconstruction: the primary method at the inert prime 3 (f is irreducible mod 3).
    assert!(matches!(
        found.algebraic_root().method(),
        ReconstructionMethod::DirectNewton { prime: 3, .. }
    ));
}

#[test]
fn kat_pipeline_factor_is_rederived_from_the_raw_relations() {
    let Run {
        context,
        search,
        report,
        ..
    } = run_pipeline();
    let FactorOutcome::Factor(found) = report.outcome() else {
        panic!("expected a verified factor, got {:?}", report.outcome());
    };
    let sqrt = context.context();
    let (n, m, c) = (sqrt.modulus(), sqrt.root(), sqrt.leading_coeff());
    assert!(sqrt.degree() > 1);
    assert!(irreducibility_certificate(sqrt.algebraic_poly()).is_some());
    assert!(
        search.candidates.len() > 1,
        "the factor was found among a bounded candidate list"
    );

    // Relation validity and lifted parity, recomputed from the validated original rows. Every
    // participating relation has b ≥ 1 and a coprime (a, b).
    let relations = context.relations().relations();
    let mut weights = vec![0u32; context.layout().num_cols()];
    let mut pairs: Vec<(BigInt, BigInt)> = Vec::new();
    for &i in found.relations() {
        let relation = &relations[i];
        assert!(relation.a().gcd(relation.b()).is_one(), "coprime pair");
        assert!(relation.b() >= &BigInt::one(), "sieve pairs have b ≥ 1");
        assert_eq!(relation.rational_norm(), &(relation.a() - relation.b() * m));
        for &col in context.original_columns(i).expect("bound relation") {
            weights[col] += 1;
        }
        pairs.push((relation.a().clone(), relation.b().clone()));
    }
    assert_eq!(pairs.len(), 22);
    assert!(
        weights.iter().all(|w| w % 2 == 0),
        "lifted parity is even in every column"
    );

    // Exact rational square: P = Π(a − b·m) > 0 and x0² = P over Z; x ≡ x0 (mod N).
    let product: BigInt = pairs.iter().map(|(a, b)| a - b * m).product();
    assert!(product.is_positive(), "P is positive");
    let x0 = found.rational_root().root();
    assert_eq!(x0 * x0, product, "P is an exact integer square");
    assert_eq!(found.rational_root().product(), &product);
    let x = found.rational_root().modular_image();
    assert_eq!(x, &x0.mod_floor(n));

    // Verified algebraic root: β² = Π(a − b·α) exactly, re-derived here from the raw pairs, and
    // the image y = β(c·m) mod N is recomputed and squares to P.
    let root = found.algebraic_root();
    assert_eq!(root.rational_product(), &product);
    assert_eq!(root.relation_count(), pairs.len());
    let beta = in_alpha_basis(root.beta(), c);
    let modulus = sqrt.algebraic_poly().to_rat_poly();
    let target = product_in_number_field(&pairs, sqrt.algebraic_poly());
    assert_eq!(
        beta.mul(&beta).rem(&modulus),
        target,
        "β² = Π(a − b·α) in Q[x]/f"
    );
    let y = root.modular_image();
    assert_eq!(y, &image_modulo_n(root.beta(), c, m, n));
    assert!(
        ((y * y) - &product).mod_floor(n).is_zero(),
        "y² ≡ P (mod N)"
    );

    // The root genuinely depends on b: with every b set to zero the product is different, and β
    // is not a rational constant.
    let without_b: Vec<(BigInt, BigInt)> = pairs
        .iter()
        .map(|(a, _)| (a.clone(), BigInt::zero()))
        .collect();
    assert_ne!(
        product_in_number_field(&without_b, sqrt.algebraic_poly()),
        target
    );
    assert!(root.beta().degree() >= Some(1), "β has a non-constant part");

    // Congruence and proper factor, with x ≢ ±y.
    assert!(
        ((x * x) - (y * y)).mod_floor(n).is_zero(),
        "x² ≡ y² (mod N)"
    );
    assert!(
        !(x - y).mod_floor(n).is_zero() && !(x + y).mod_floor(n).is_zero(),
        "x ≢ ±y"
    );
    assert_eq!(
        *found.divisor(),
        (x - y).abs().gcd(n),
        "the divisor is gcd(|x − y|, N)"
    );
}
