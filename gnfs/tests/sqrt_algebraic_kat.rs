//! Known-answer tests for the verified algebraic square root of a checked dependency.
//!
//! `algebraic_sqrt` reconstructs `β` with `β² = γ = ∏(a − b·α)` and accepts it only through the
//! exact verifier. The modular image is `y = β(c·m) mod N`, the normalized root evaluated at
//! `θ = c·m`; it satisfies `y² ≡ P = ∏(a − b·m) (mod N)`. It is not `Norm(β)`.
//!
//! **Fixtures** (relation sets are built through the checked pipeline; see the rational root
//! tests for the construction):
//!
//! - **S** (seam): `x² − 2`, `m = 23`, `N = 527 = 17·31`, relation (3,−2). `γ = 3 + 2α =
//!   (1 + α)²`, `P = 3 + 46 = 49`, `β = ±(1 + α)`, `y = ±(1 + 23) = ±24`, and `24² = 576 =
//!   527 + 49`.
//! - **R** (eight relations, `x² − 2`, `m = 23`, `N = 527`, `FB(13, 13)`): (10,−1), (−3,2),
//!   (−2,−1), (−2,1), (4,3), (2,−1), (10,1), (1,0). Row 5 is (2,−1): `γ = 2 + α = α(1 + α)`,
//!   where `1 + α` is a unit and `α` generates the prime above 2, so `γ` has odd valuation 1 and
//!   is not a square. That ideal is exceptional (`2 | disc`), so no matrix column constrains it
//!   and row 5 is a kernel vector anyway.
//! - **NM** (non-monic): `f = 2x² − 2x − 5`, `m = 5`, `N = 35 = 5·7` (`f(5) = 35`), `F = t² − 2t
//!   − 10`, `θ = 2α`, `disc = 44`. Relations (9,−11), (−1,−1): `γ = (9 + 11α)(−1 + α)`; with
//!   `α² = α + 5/2` this is `37/2 + 9α = (4 + α)²`. So `β = ±(4 + θ/2)`, `P = 64·4 = 256`,
//!   `y = ±(4 + 5) = ±9`, and `9² = 81 ≡ 11 ≡ 256 (mod 35)`.
//! - **Cubic**: `x³ − 2`, `m = 7`, `N = 341 = 11·31` (`f(7) = 341`), pairs (3,1), (−1,−1):
//!   `γ = (3 − α)(−1 + α) = −3 + 4α − α² = (1 + α − α²)²` (expanding with `α³ = 2`).
//!   `y = ±(1 + 7 − 49) = ±(−41)`, i.e. `300` or `41` modulo 341.
//! - **Imaginary**: `x² + 1`, `m = 8`, `N = 65` (`f(8) = 65`), pair (−3,−4): `γ = −3 + 4α =
//!   (1 + 2α)²`, `y = ±17`, `P = −3 + 32 = 29`, and `17² = 289 = 4·65 + 29`.
//!
//! Hand values were recomputed with an independent exact-arithmetic script (Python `fractions`
//! and `sympy`, no Rust). Pinned limits keep `crt.max_primes = 1` and ample precision headroom.

use gnfs::filter::build_checked_matrix;
use gnfs::linalg::RootCharacters;
use gnfs::sieve::{CheckedRelations, ExponentVector, FactorBase, Relation, RelationChecker};
use gnfs::sqrt::{
    AlgebraicSqrtOutcome, CheckedDependency, CrtLimits, DependencyContext, DirectLimits,
    NonSquareEvidence, NormalizedTarget, PrecisionLimits, PrimeSearchLimits, RecombinationLimits,
    ReconstructionMethod, RelationFactors, ResourceKind, SqrtContext, SqrtLimits, TargetOutcome,
    VerifiedRoot, algebraic_sqrt, direct_root, reconstruct_root, reconstruct_root_from_factors,
    split_prime_root,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, Zero};
use shared_numfield::{IntPoly, NonresidueSearchLimits, RatPoly, RootSearchLimits};

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
/// factor base; the checker re-verifies every norm.
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

fn factors(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> RelationFactors {
    let ctx = SqrtContext::new(&poly(f), &big(m), &big(n)).expect("valid context");
    let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
    RelationFactors::new(&ctx, &pairs).expect("valid factors")
}

fn validated(context: &DependencyContext, rows: &[usize]) -> CheckedDependency {
    context.validate(rows).unwrap_or_else(|e| panic!("{rows:?} is not a dependency: {e}"))
}

fn prime_search(first: u64, examined: u64) -> PrimeSearchLimits {
    PrimeSearchLimits { first_candidate: first, max_primes_examined: examined }
}

fn nonresidue() -> NonresidueSearchLimits {
    NonresidueSearchLimits { max_candidates: 1000 }
}

fn limits_with(direct_first: u64, direct_examined: u64, crt_first: u64) -> SqrtLimits {
    SqrtLimits {
        precision: PrecisionLimits { max_modulus_bits: 4096 },
        direct: DirectLimits {
            prime_search: prime_search(direct_first, direct_examined),
            nonresidue_search: nonresidue(),
        },
        crt: CrtLimits {
            prime_search: prime_search(crt_first, 200),
            max_primes: 1,
            root_search: RootSearchLimits { max_evaluations: 100_000 },
            nonresidue_search: nonresidue(),
        },
        recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 20 },
    }
}

/// Direct reconstruction from the first prime; split primes, one at a time, as the fallback.
fn ample() -> SqrtLimits {
    limits_with(0, 100, 0)
}

fn verified(outcome: AlgebraicSqrtOutcome) -> VerifiedRoot {
    match outcome {
        AlgebraicSqrtOutcome::Verified(root) => root,
        other => panic!("expected a verified root, got {other:?}"),
    }
}

fn target_of(factors: RelationFactors) -> NormalizedTarget {
    match NormalizedTarget::new(factors, PrecisionLimits { max_modulus_bits: 4096 }).unwrap() {
        TargetOutcome::Ready(target) => target,
        TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
    }
}

/// `β` equals `expected` or `−expected`, and the image is one of the two signed values.
fn assert_root(root: &VerifiedRoot, expected_beta: &RatPoly, images: [i64; 2]) {
    assert!(
        root.beta() == expected_beta || root.beta() == &negated(expected_beta),
        "β = {:?} is neither ±{expected_beta:?}",
        root.beta()
    );
    let image = root.modular_image();
    assert!(
        *image == big(images[0]) || *image == big(images[1]),
        "image {image} is neither {} nor {}",
        images[0],
        images[1]
    );
}

#[test]
fn kat_algebraic_root_of_the_seam_dependency_is_verified() {
    let dependency = validated(&context_s(), &[0]);
    let root = verified(algebraic_sqrt(&dependency, ample()).unwrap());
    assert_eq!(root.rational_product(), &big(49));
    assert_eq!(root.relation_count(), 1);
    assert_eq!(root.context(), dependency.context());
    // β = ±(1 + θ) and y = ±24, that is 24 or 527 − 24 = 503.
    assert_root(&root, &rat(&[(1, 1), (1, 1)]), [24, 503]);
    // The first prime p ≥ 3 at which x² − 2 is irreducible is 3 (2 is a nonresidue mod 3).
    assert!(matches!(root.method(), ReconstructionMethod::DirectNewton { prime: 3, .. }));
}

#[test]
fn kat_modular_image_is_the_normalized_map_not_the_norm() {
    // γ = 3 + 2α, β = 1 + α. Norm(β) = (1 + √2)(1 − √2) = −1 and Norm(γ) = 9 − 8 = 1. A map by
    // Norm would give |Norm(β)| = 1 with 1² ≢ P = 49 (mod 527). The normalized map gives
    // y = β(23) = 24, and 24² = 576 ≡ 49 = P.
    let dependency = validated(&context_s(), &[0]);
    let root = verified(algebraic_sqrt(&dependency, ample()).unwrap());
    let n = dependency.context().modulus();
    let y = root.modular_image();
    assert_eq!((y * y) % n, big(49) % n);
    let norm_gamma = big(1);
    assert_ne!((y * y) % n, norm_gamma % n, "y² is P, not Norm(γ)");
    assert_ne!(y, &big(1), "y is not |Norm(β)|");
}

#[test]
fn kat_legacy_fixture_mod_seven_uses_the_same_map() {
    // x² − 2, m = 3, N = 7 (f(3) = 7), pair (9,0): γ = 9 = 3², P = 9, β = ±3, y = ±3 so 3 or 4
    // mod 7. The former map gave |Norm(3)| = 9 ≡ 2 (mod 7), and 2² = 4 ≢ P = 9 ≡ 2. Here
    // 3² = 9 ≡ 2 ≡ P.
    let outcome = reconstruct_root_from_factors(factors(&[-2, 0, 1], 3, 7, &[(9, 0)]), ample());
    let root = verified(outcome.unwrap());
    assert_eq!(root.rational_product(), &big(9));
    assert_root(&root, &rat(&[(3, 1)]), [3, 4]);
    let y = root.modular_image();
    assert_eq!((y * y) % big(7), big(2));
    assert_ne!((big(2) * big(2)) % big(7), big(9) % big(7));
}

#[test]
fn kat_totally_imaginary_field_needs_no_real_embedding() {
    // x² + 1 has no real root, so a real-embedding sign rule is undefined. γ = −3 + 4α =
    // (1 + 2α)², P = 29, y = ±17 so 17 or 48 mod 65, and 17² = 289 ≡ 29 (mod 65).
    let outcome = reconstruct_root_from_factors(factors(&[1, 0, 1], 8, 65, &[(-3, -4)]), ample());
    let root = verified(outcome.unwrap());
    assert_eq!(root.rational_product(), &big(29));
    assert_root(&root, &rat(&[(1, 1), (2, 1)]), [17, 48]);
}

#[test]
fn kat_nonsquare_is_reported_only_with_a_certificate() {
    // R row 5 is (2,−1): γ = 2 + α has odd valuation at √2, so it is not a square in Q(√2). The
    // rational product is 25 = 5², so the kernel vector is genuine on the rational side.
    let r = context_r();
    let dependency = validated(&r, &[5]);
    assert_eq!(dependency.rational_product(), &big(25));
    match algebraic_sqrt(&dependency, ample()).unwrap() {
        AlgebraicSqrtOutcome::NotSquare(certificate) => {
            assert!(matches!(
                certificate.evidence(),
                NonSquareEvidence::InertLocalNonresidue { .. }
                    | NonSquareEvidence::DirectBoundExhausted { .. }
            ));
        }
        other => panic!("expected a nonsquare certificate, got {other:?}"),
    }
    // Rows {1,3}: (−3,2)(−2,1) = (−3 − 2α)(−2 − α) = 6 + 3α + 4α + 2α² = 10 + 7α. Norm =
    // 100 − 98 = 2, again an odd valuation at √2.
    let dependency = validated(&r, &[1, 3]);
    assert!(matches!(
        algebraic_sqrt(&dependency, ample()).unwrap(),
        AlgebraicSqrtOutcome::NotSquare(_)
    ));
}

#[test]
fn kat_entry_point_is_reconstruction_from_the_dependency_factors() {
    let r = context_r();
    for rows in [&[7usize][..], &[1, 3, 5], &[5]] {
        let dependency = validated(&r, rows);
        let direct = algebraic_sqrt(&dependency, ample());
        let via_factors = reconstruct_root_from_factors(dependency.factors().clone(), ample());
        assert_eq!(direct, via_factors, "{rows:?}");
        // Deterministic: the same call gives the same outcome.
        assert_eq!(direct, algebraic_sqrt(&dependency, ample()));
    }
    // Rows {1,3,5}: (−3 − 2α)(−2 − α) = 10 + 7α, and row 5 contributes 2 − (−1)α = 2 + α, so
    // γ = (10 + 7α)(2 + α) = 20 + 10α + 14α + 14 = 34 + 24α = (4 + 3α)² (16 + 24α + 18).
    // P = 30625, y = 4 + 3·23 = 73, and 73² = 5329 ≡ 59 ≡ 30625 (mod 527).
    let dependency = validated(&r, &[1, 3, 5]);
    let root = verified(algebraic_sqrt(&dependency, ample()).unwrap());
    assert_root(&root, &rat(&[(4, 1), (3, 1)]), [73, 454]);
}

#[test]
fn kat_precision_exhaustion_is_resource_limited_never_a_nonsquare() {
    let dependency = validated(&context_s(), &[0]);
    let mut limits = ample();
    limits.precision = PrecisionLimits { max_modulus_bits: 1 };
    match algebraic_sqrt(&dependency, limits).unwrap() {
        AlgebraicSqrtOutcome::ResourceExhausted(spent) => {
            assert_eq!(spent.len(), 1);
            assert_eq!(spent[0].kind, ResourceKind::Precision);
            assert_eq!(spent[0].limit, 1);
            assert!(spent[0].consumed > 1);
        }
        other => panic!("expected resource exhaustion, got {other:?}"),
    }
}

#[test]
fn kat_non_monic_target_through_direct_and_split_primes() {
    // NM: c = 2, β = ±(4 + θ/2), P = 256, y = ±9 so 9 or 26 mod 35. 11 ≡ disc/4 is a nonresidue
    // mod 3, so F is irreducible mod 3: the direct method succeeds at p = 3.
    let context = context_nm();
    let dependency = validated(&context, &[0, 1]);
    assert_eq!(dependency.rational_product(), &big(256));
    let expected = rat(&[(4, 1), (1, 2)]);
    let direct = verified(algebraic_sqrt(&dependency, ample()).unwrap());
    assert_root(&direct, &expected, [9, 26]);
    assert!(matches!(direct.method(), ReconstructionMethod::DirectNewton { prime: 3, .. }));

    // Forcing the fallback: the direct search examines only p = 5, which divides N and is
    // skipped, so it is exhausted and the split-prime method runs. F splits completely mod p
    // exactly when 11 is a square mod p. Among the primes not dividing c·N·disc = 2·35·44, 11 is a
    // nonresidue mod 3, 13 (squares {1,3,4,9,10,12}) and 17 (squares {1,2,4,8,9,13,15,16}), and a
    // residue mod 19 (7² = 49 = 2·19 + 11), so 19 is the first split prime.
    let forced = limits_with(5, 1, 0);
    let split = verified(algebraic_sqrt(&dependency, forced).unwrap());
    assert_root(&split, &expected, [9, 26]);
    assert!(
        matches!(split.method(), ReconstructionMethod::SplitPrimeCrt { primes, .. } if primes == &[19]),
        "{:?}",
        split.method()
    );

    // The two methods agree on the root up to sign, and on the image up to the same sign.
    assert!(split.beta() == direct.beta() || split.beta() == &negated(direct.beta()));
    assert_eq!(split.beta() == direct.beta(), split.modular_image() == direct.modular_image());
}

#[test]
fn kat_cross_method_agreement_for_a_cubic() {
    // disc(x³ − 2) = −108 = −2²·3³, so p = 3 is skipped. x³ − 2 has the root 3 mod 5 (27 ≡ 2),
    // so it is not irreducible there, and it is irreducible mod 7 (the cubes mod 7 are
    // {0, 1, 6}, which omit 2): the direct method runs at p = 7. The first completely split
    // prime coprime to N = 11·31 is 43. Cubing is a bijection mod p ≡ 2 (mod 3), so such p give
    // one root; 13, 19, 37 are ≡ 1 (mod 3) with 2 a non-cube (2¹² ≡ 26 mod 37; cubes mod 13 are
    // {1,5,8,12}, mod 19 {1,7,8,11,12,18}); 31 divides N; and 2 is a cube mod 43 because
    // 2⁷ = 128 ≡ −1 so 2¹⁴ ≡ 1. β = ±(1 + α − α²) and y = ±41, i.e. 41 or 300 mod 341.
    let target = target_of(factors(&[-2, 0, 0, 1], 7, 341, &[(3, 1), (-1, -1)]));
    let limits = ample();
    let direct = verified(direct_root(&target, limits.direct).unwrap());
    let split = verified(split_prime_root(&target, limits.crt, limits.recombination).unwrap());
    let dispatched = verified(reconstruct_root(&target, limits).unwrap());
    let expected = rat(&[(1, 1), (1, 1), (-1, 1)]);
    for root in [&direct, &split, &dispatched] {
        assert_root(root, &expected, [41, 300]);
        assert_eq!(root.rational_product(), &big(-24));
        let y = root.modular_image();
        assert_eq!((y * y) % big(341), big(317), "−24 ≡ 317 and 41² = 1681 = 4·341 + 317");
    }
    assert!(matches!(direct.method(), ReconstructionMethod::DirectNewton { prime: 7, .. }));
    assert!(
        matches!(split.method(), ReconstructionMethod::SplitPrimeCrt { primes, .. } if primes == &[43]),
        "{:?}",
        split.method()
    );
    assert_eq!(dispatched.method(), direct.method(), "the dispatcher prefers the direct method");
}
