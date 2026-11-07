//! Known-answer tests for direct p-adic square-root reconstruction.
//!
//! Every fixture states the arithmetic that justifies it, derived by hand and not by running the
//! code under test. Throughout, `θ = c·α` is the root of the monic `F(t) = c^(d−1) f(t/c)`,
//! `γ = Π(a − b·α)`, `η = c^(2h−k) Π(c·a − b·θ)` with `h = ⌈k/2⌉`, `H = F'(θ)² η`, and a root of
//! the target is `z = c^h F'(θ) β` with `β² = γ`. The coefficient bound is
//! `R = 1 + max|F_i|`, `E = Σ|η_i| R^i`, `U = ⌈√E⌉`, `B_i = d U Σ_{ℓ>i} |F_ℓ| R^(ℓ−1−i)`, and the
//! lifting modulus is the least power `p^e` above `2B`. Only a global sign of the root is free, so
//! each assertion accepts `z` and `−z`. The private quotient ring is not reachable from here, so
//! the per-precision Newton invariant `H·u² ≡ 1` is pinned by a unit test beside the code.

use gnfs::sqrt::{
    AlgebraicSqrtOutcome, DirectLimits, FailureClass, NonSquareEvidence, NormalizedTarget,
    PrecisionLimits, PrimeSearchLimits, ReconstructionMethod, RelationFactors, ResourceExhaustion,
    ResourceKind, SqrtContext, TargetOutcome, VerifiedRoot, direct_root,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::Signed;
use shared_numfield::{IntPoly, NonresidueSearchLimits, RatPoly};

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Integer polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

/// Rational polynomial from `(numerator, denominator)` pairs, lowest degree first.
fn rat_poly(coeffs: &[(i64, i64)]) -> RatPoly {
    RatPoly::from_coeffs(
        coeffs
            .iter()
            .map(|&(n, d)| BigRational::new(big(n), big(d)))
            .collect(),
    )
}

fn context(f: &[i64], m: i64, n: i64) -> SqrtContext {
    SqrtContext::new(&poly(f), &big(m), &big(n)).expect("fixture is a valid context")
}

fn factors(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> RelationFactors {
    let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
    RelationFactors::new(&context(f, m, n), &pairs).expect("fixture relations are units")
}

fn target_within(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)], max_bits: u64) -> TargetOutcome {
    NormalizedTarget::new(
        factors(f, m, n, pairs),
        PrecisionLimits {
            max_modulus_bits: max_bits,
        },
    )
    .expect("a validated target builds")
}

fn target(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> NormalizedTarget {
    match target_within(f, m, n, pairs, 4096) {
        TargetOutcome::Ready(target) => target,
        TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
    }
}

/// Primes from `first` on, at most `primes` of them examined, and a nonresidue budget.
fn limits(first: u64, primes: u64, nonresidue: u64) -> DirectLimits {
    DirectLimits {
        prime_search: PrimeSearchLimits {
            first_candidate: first,
            max_primes_examined: primes,
        },
        nonresidue_search: NonresidueSearchLimits {
            max_candidates: nonresidue,
        },
    }
}

/// Starts at 3 with room for a hundred primes and a thousand nonresidue candidates.
fn generous() -> DirectLimits {
    limits(0, 100, 1000)
}

fn run(target: &NormalizedTarget, limits: DirectLimits) -> AlgebraicSqrtOutcome {
    direct_root(target, limits).expect("no invariant violation")
}

fn verified(outcome: AlgebraicSqrtOutcome) -> VerifiedRoot {
    match outcome {
        AlgebraicSqrtOutcome::Verified(root) => root,
        other => panic!("expected a verified root, got {other:?}"),
    }
}

fn exhausted(outcome: AlgebraicSqrtOutcome) -> Vec<ResourceExhaustion> {
    match outcome {
        AlgebraicSqrtOutcome::ResourceExhausted(list) => list,
        other => panic!("expected resource exhaustion, got {other:?}"),
    }
}

fn evidence(outcome: AlgebraicSqrtOutcome) -> NonSquareEvidence {
    match outcome {
        AlgebraicSqrtOutcome::NotSquare(certificate) => certificate.evidence().clone(),
        other => panic!("expected a nonsquare certificate, got {other:?}"),
    }
}

fn direct(prime: u64, exponent: u32) -> ReconstructionMethod {
    ReconstructionMethod::DirectNewton { prime, exponent }
}

/// `z` equals `expected` or its negation.
fn assert_root_up_to_sign(root: &VerifiedRoot, expected: &[i64]) {
    let z = root.integral_root();
    assert!(
        *z == poly(expected) || *z == poly(expected).neg(),
        "z = {z:?} is neither ±{expected:?}"
    );
}

fn assert_image_is_one_of(root: &VerifiedRoot, images: [i64; 2]) {
    let y = root.modular_image();
    assert!(
        *y == big(images[0]) || *y == big(images[1]),
        "y = {y} is neither of {images:?}"
    );
}

/// `candidate² ≡ γ (mod F)` over `Q`, and the reconstructed `β` is `±candidate`.
fn assert_agrees_with_square(root: &VerifiedRoot, target: &NormalizedTarget, candidate: &RatPoly) {
    let modulus = target.context().monic_poly().to_rat_poly();
    assert_eq!(&candidate.mul(candidate).rem(&modulus), target.gamma());
    assert!(
        root.beta() == candidate || *root.beta() == candidate.neg(),
        "β = {:?} is neither ±{candidate:?}",
        root.beta()
    );
}

/// `x² − 2`, `m = 23`, `N = 527 = 17·31`, relation `(3, −2)`.
fn monic_quadratic() -> NormalizedTarget {
    target(&[-2, 0, 1], 23, 527, &[(3, -2)])
}

/// `x² − 2`, `m = 3`, `N = 7`, relation `(3, 0)` twice: `H = 72`, which vanishes modulo 3.
fn vanishing_mod_three() -> NormalizedTarget {
    target(&[-2, 0, 1], 3, 7, &[(3, 0), (3, 0)])
}

/// `x² + x + 1`, `m = 2`, `N = 7`, relation `(1, 1)` twice: `F` is irreducible modulo 2 and 5 but
/// not modulo 3.
fn irreducible_mod_two() -> NormalizedTarget {
    target(&[1, 1, 1], 2, 7, &[(1, 1), (1, 1)])
}

#[test]
fn kat_nonconstant_root_monic_quadratic() {
    // f = x² − 2, m = 23, N = 527, relation (3, −2): γ = 3 + 2θ = (1 + θ)² since
    // (1 + θ)² = 1 + 2θ + θ² = 3 + 2θ. With F' = 2t: z = 2θ(1 + θ) = 2θ + 2θ² = 4 + 2θ, and
    // z² = 16 + 16θ + 4θ² = 24 + 16θ = 8(3 + 2θ) = H. The root is not constant.
    // Bound: R = 3, E = 3 + 2·3 = 9, U = 3, B = (18, 6), so the minimum modulus is 37. Modulo 3
    // the polynomial t² − 2 ≡ t² + 1 is irreducible (−1 is a nonresidue) and H ≡ θ ≠ 0, so p = 3
    // is the first usable prime. 3³ = 27 < 37 ≤ 81 = 3⁴ gives e = 4.
    // y = β(23) = 24, P = 3 + 46 = 49 and 24² = 576 = 527 + 49.
    let t = monic_quadratic();
    assert_eq!(t.target(), &poly(&[24, 16]));
    assert_eq!(t.minimum_modulus(), &big(37));
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert_root_up_to_sign(&root, &[4, 2]);
    assert_image_is_one_of(&root, [24, 527 - 24]);
    assert_eq!(root.rational_product(), &big(49));
    assert_eq!(root.relation_count(), 1);
    assert_eq!(root.context(), t.context());
}

#[test]
fn kat_non_monic_root() {
    // f = 2x² − 1, m = 2, N = 7, relation (3, −4): c = 2, F = t² − 2, γ = 3 + 2θ = (1 + θ)²,
    // H = 96 + 64θ and z = 2·F'·β = 4θ(1 + θ) = 8 + 4θ with z² = 64 + 64θ + 16θ² = 96 + 64θ.
    // Bounds B = (36, 12), minimum modulus 73. The prime 3 does not divide c·N = 14, F ≡ t² + 1 is
    // irreducible mod 3, and H ≡ θ ≠ 0. 3³ = 27 < 73 ≤ 81 gives e = 4. y = β(c·m = 4) = 5 and
    // 5² = 25 ≡ 4 ≡ 11 = P (mod 7).
    let t = target(&[-1, 0, 2], 2, 7, &[(3, -4)]);
    assert_eq!(t.minimum_modulus(), &big(73));
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert_root_up_to_sign(&root, &[8, 4]);
    assert_image_is_one_of(&root, [5, 2]);
    assert_agrees_with_square(&root, &t, &rat_poly(&[(1, 1), (1, 1)]));
}

#[test]
fn kat_even_degree_quartic_root() {
    // f = x⁴ + x − 1 with m = 2 and N = 17 = f(2). F = t⁴ + t − 1. Mod 3 it has no root
    // (values −1, −1, 1 at 0, 1, 2) and is not a product of two monic quadratics, so it is
    // irreducible and 3 ∤ 17 is the first inert prime. Relation (1, −1) twice: γ = (1 + α)²,
    // k = 2, h = 1, η = (1 + θ)², β = 1 + θ. With F' = 4t³ + 1 and θ⁴ = 1 − θ:
    // z = F'(1 + θ) = 4θ³ + 4θ⁴ + 1 + θ = 4θ³ + 4 − 4θ + 1 + θ = 5 − 3θ + 4θ³.
    // Bound: R = 2, E = 1 + 2·2 + 1·4 = 9, U = 3; with d = 4 the sums Σ|F_ℓ|R^(ℓ−1−i) are
    // (9, 4, 2, 1), so B = 12·(9, 4, 2, 1) = (108, 48, 24, 12). The minimum modulus is 217;
    // 3⁴ = 81 < 217 ≤ 243 = 3⁵ gives e = 5. y = β(2) = 3 and 3² = 9 = P.
    let t = target(&[-1, 1, 0, 0, 1], 2, 17, &[(1, -1), (1, -1)]);
    assert_eq!(t.minimum_modulus(), &big(217));
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 5));
    assert_root_up_to_sign(&root, &[5, -3, 0, 4]);
    assert_image_is_one_of(&root, [3, 17 - 3]);
    assert_eq!(root.rational_product(), &big(9));
    assert_agrees_with_square(&root, &t, &rat_poly(&[(1, 1), (1, 1)]));
}

#[test]
fn kat_totally_imaginary_root() {
    // f = x² + 1, m = 2, N = 5 = f(2): the field Q(i) has no real embedding, and no sign or
    // real-root convention is used. Relation (3, −4): γ = 3 + 4θ = (2 + θ)² because
    // (2 + θ)² = 4 + 4θ + θ² = 3 + 4θ. F' = 2t and F'² = 4t² ≡ −4, so H = −4(3 + 4θ) =
    // −12 − 16θ and z = 2θ(2 + θ) = 4θ + 2θ² = −2 + 4θ, with z² = 4 − 16θ + 16θ² = −12 − 16θ.
    // Bound: R = 2, E = 3 + 4·2 = 11, U = 4, B = (16, 8), minimum modulus 33. Mod 3, t² + 1 is
    // irreducible and H ≡ 2θ ≠ 0. 3³ = 27 < 33 ≤ 81 gives e = 4. y = β(2) = 4, and
    // 4² = 16 ≡ 1 ≡ 11 = P (mod 5).
    let t = target(&[1, 0, 1], 2, 5, &[(3, -4)]);
    assert!(
        t.context().discriminant().is_negative(),
        "disc(t² + 1) = −4 < 0: no real roots"
    );
    assert_eq!(t.minimum_modulus(), &big(33));
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert_root_up_to_sign(&root, &[-2, 4]);
    assert_image_is_one_of(&root, [4, 1]);
    assert_agrees_with_square(&root, &t, &rat_poly(&[(2, 1), (1, 1)]));
}

#[test]
fn kat_root_outside_the_power_basis_order() {
    // f = x² − 8, m = 6, N = 7 (f(6) = 28), relation (2, 0): γ = 2, β = θ/2 because
    // θ²/4 = 8/4 = 2, and β lies outside Z[θ]. F' = 2t, so z = F'β = θ² = 8 and H = 4t²·2 = 64.
    // Bound: R = 9, E = 2, U = 2, B = (36, 4), minimum modulus 73; mod 3, t² − 8 ≡ t² + 1 is
    // irreducible and H ≡ 1 ≠ 0. 3⁴ = 81 ≥ 73 > 27 gives e = 4. y = β(6) = 6/2 = 3 and
    // 3² = 9 ≡ 2 = P (mod 7).
    let t = target(&[-8, 0, 1], 6, 7, &[(2, 0)]);
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert_root_up_to_sign(&root, &[8]);
    assert_image_is_one_of(&root, [3, 4]);
    assert_agrees_with_square(&root, &t, &rat_poly(&[(0, 1), (1, 2)]));
}

#[test]
fn kat_precision_boundary_just_above_and_below_twice_the_bound() {
    // f = x − 5 (d = 1, F = t − 5, F' = 1), m = 5, N = 11. A relation (n, 0) has γ = η = H = n
    // and bound B = ⌈√n⌉ (R = 6, B_0 = d·U·|F_1| = U). With n = s² the roots are ±s exactly.
    // The modulus must exceed 2B: only 3^e above 2s works.
    //
    // n = 169 = 13²: B = 13, minimum 2B + 1 = 27 = 3³. The modulus equals 2B + 1 and the root has
    // |z| = 13 = B: the centered range [−13, 13] of Z/27 holds it, so e = 3 succeeds.
    let exact = target(&[-5, 1], 5, 11, &[(169, 0)]);
    assert_eq!(exact.bound(), &big(13));
    assert_eq!(exact.minimum_modulus(), &big(27));
    let root = verified(run(&exact, generous()));
    assert_eq!(root.method(), &direct(3, 3));
    assert_root_up_to_sign(&root, &[13]);

    // n = 196 = 14²: B = 14, minimum 29. Now 3³ = 27 = 2B − 1 is just below 2B: the centered
    // range of Z/27 is [−13, 13], and 14 ≡ −13 (mod 27) would reconstruct −13, whose square 169
    // is not 196. So e = 4 (81) is required and succeeds with z = ±14.
    let above = target(&[-5, 1], 5, 11, &[(196, 0)]);
    assert_eq!(above.bound(), &big(14));
    assert_eq!(above.minimum_modulus(), &big(29));
    assert_eq!(centered_residue(14, 27), big(-13));
    assert_ne!(centered_residue(14, 27).pow(2), big(196));
    let root = verified(run(&above, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert_root_up_to_sign(&root, &[14]);
}

/// The representative of `value` in `[−(M−1)/2, (M−1)/2]` for odd `modulus`.
fn centered_residue(value: i64, modulus: i64) -> BigInt {
    let r = big(value).mod_floor(&big(modulus));
    if &r * 2 > big(modulus) {
        r - big(modulus)
    } else {
        r
    }
}

#[test]
fn kat_local_nonsquare_certificate() {
    // f = x² − 2, m = 4, N = 7 (f(4) = 14), relation (1, 1): γ = 1 − θ with
    // Norm = 1 − 2 = −1, not a rational square, so γ is not a square in Q(√2). H = 8(1 − θ) =
    // 8 − 8θ ≡ 2 + θ (mod 3). In F_9 = F_3[t]/(t² − 2) an element x is a square iff
    // x^((9−1)/2) = Norm(x)^((3−1)/2) = Norm(x) is 1 (Euler); Norm(2 + θ) = 4 − 2 = 2 ≡ −1, a
    // nonresidue. So p = 3 gives a local obstruction, with no lifting.
    let t = target(&[-2, 0, 1], 4, 7, &[(1, 1)]);
    let outcome = run(&t, generous());
    assert_eq!(outcome.class(), Some(FailureClass::AlgebraicNonSquare));
    assert_eq!(
        evidence(outcome),
        NonSquareEvidence::InertLocalNonresidue { prime: 3 }
    );
}

#[test]
fn kat_global_nonsquare_certificate_after_full_bound() {
    // f = x² − 2, m = 4, N = 7, relation (3, 1): γ = 3 − θ with Norm = 9 − 2 = 7, not a rational
    // square, so γ is not a square in Q(√2). Yet H = 8(3 − θ) ≡ θ (mod 3) is a square in F_9,
    // since Norm(θ) = −2 ≡ 1: there is no local obstruction, and the nonsquare is found only
    // after a root modulo 3^e is lifted and its centered form fails at M > 2B.
    // Bound: R = 3, E = 3 + 3 = 6, U = ⌈√6⌉ = 3, B = (18, 6), minimum modulus 37, e = 4.
    let t = target(&[-2, 0, 1], 4, 7, &[(3, 1)]);
    assert_eq!(t.minimum_modulus(), &big(37));
    let outcome = run(&t, generous());
    assert_eq!(outcome.class(), Some(FailureClass::AlgebraicNonSquare));
    assert_eq!(
        evidence(outcome),
        NonSquareEvidence::DirectBoundExhausted {
            prime: 3,
            exponent: 4
        }
    );
}

#[test]
fn kat_prime_search_exhaustion_is_a_resource_outcome() {
    // F = t⁴ − 10t² + 1 (f = x⁴ − 10x² + 1, m = 2, N = 23 = −f(2)) is the minimal polynomial of
    // √2 + √3 and has Galois group V4: Frobenius has order at most 2, so F is reducible modulo
    // every prime. No inert prime exists, and a bounded search must not call that a nonsquare.
    let t = target(&[1, 0, -10, 0, 1], 2, 23, &[(1, 1)]);
    let outcome = run(&t, limits(0, 50, 1000));
    assert_eq!(outcome.class(), Some(FailureClass::ResourceExhausted));
    assert_eq!(
        exhausted(outcome),
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 50,
            limit: 50
        }]
    );

    // A zero budget examines nothing.
    let none = exhausted(run(&monic_quadratic(), limits(0, 0, 1000)));
    assert_eq!(
        none,
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 0,
            limit: 0
        }]
    );
}

#[test]
fn kat_machine_prime_range_end_is_a_resource_outcome() {
    // 2^64 − 59 is the largest prime below 2^64, so a search starting above it ends the range
    // before any prime is examined; the budget (5) was not spent.
    let t = monic_quadratic();
    let outcome = run(&t, limits(u64::MAX - 57, 5, 1000));
    assert_eq!(outcome.class(), Some(FailureClass::ResourceExhausted));
    assert_eq!(
        exhausted(outcome),
        vec![ResourceExhaustion {
            kind: ResourceKind::MachinePrimeRange,
            consumed: 0,
            limit: 5
        }]
    );
}

#[test]
fn kat_machine_prime_near_u64_max_reconstructs() {
    // p = 2^64 − 59 ≡ 1 (mod 4). With F = t − 5 (d = 1) it is inert, 169 = 13² is a residue,
    // and Tonelli–Shanks needs a nonresidue. The minimum modulus 27 is far below p, so e = 1,
    // and the centered residue of ±13 is ±13.
    let largest = u64::MAX - 58;
    let t = target(&[-5, 1], 5, 11, &[(169, 0)]);
    let root = verified(run(&t, limits(largest, 1, 1000)));
    assert_eq!(root.method(), &direct(largest, 1));
    assert_root_up_to_sign(&root, &[13]);
}

#[test]
fn kat_nonresidue_budget_exhaustion_is_a_resource_outcome() {
    // In F_9 the odd part of q − 1 = 8 is 1, and H ≡ θ ≠ 1, so Tonelli–Shanks needs a nonresidue.
    // A zero budget cannot supply one; that is a resource outcome and says nothing about H.
    let outcome = run(&monic_quadratic(), limits(0, 100, 0));
    assert_eq!(outcome.class(), Some(FailureClass::ResourceExhausted));
    assert_eq!(
        exhausted(outcome),
        vec![ResourceExhaustion {
            kind: ResourceKind::NonresidueSearch,
            consumed: 0,
            limit: 0
        }]
    );
}

#[test]
fn kat_precision_limit_exhaustion_is_a_resource_outcome() {
    // For the (3, −2) relation over x² − 2 the minimum modulus 37 has 6 bits. One bit fewer is
    // refused before any prime search; it is a resource outcome, never a nonsquare. At exactly 6
    // bits the direct method runs, and its modulus 3⁴ = 81 has 7 bits: at most bits(2B + 1) +
    // bits(p) = 6 + 2.
    match target_within(&[-2, 0, 1], 23, 527, &[(3, -2)], 5) {
        TargetOutcome::ResourceExhausted(e) => assert_eq!(
            e,
            ResourceExhaustion {
                kind: ResourceKind::Precision,
                consumed: 6,
                limit: 5
            }
        ),
        TargetOutcome::Ready(_) => panic!("5 bits cannot hold 37"),
    }
    let TargetOutcome::Ready(t) = target_within(&[-2, 0, 1], 23, 527, &[(3, -2)], 6) else {
        panic!("6 bits hold 37");
    };
    let root = verified(run(&t, generous()));
    assert_eq!(root.method(), &direct(3, 4));
    assert!(big(81).bits() <= t.minimum_modulus().bits() + big(3).bits());
}

#[test]
fn kat_two_is_never_examined() {
    // f = x² + x + 1, m = 2, N = 7, relation (1, 1) twice: γ = (1 − α)², P = 1, η = (1 − θ)² =
    // 1 − 2θ + θ² = −3θ (θ² = −θ − 1), F' = 2θ + 1 and F'² = 4θ² + 4θ + 1 = −3, so H = 9θ and
    // z = F'(1 − θ) = 3θ + 3 (check: 9θ² + 18θ + 9 = 9θ). Bound: R = 2, E = 0 + 3·2 = 6, U = 3,
    // B = (18, 6), minimum modulus 37.
    // t² + t + 1 is irreducible mod 2, but 2 is excluded. Mod 3 it is (t − 1)², reducible; mod 5
    // its discriminant −3 ≡ 2 is a nonresidue, so it is irreducible and H ≡ 4θ ≠ 0. 5² = 25 < 37
    // ≤ 125 gives e = 3. Starting the search at 0 or 2 both begin at 3: the examined primes are
    // 3 and 5, so a budget of 1 runs out while a budget of 2 succeeds. β = z/F' = 3(1 + θ)/
    // (1 + 2θ) = ±(1 − θ), y = β(2) = −1 ≡ 6 (mod 7), and y² = 1 = P.
    let t = irreducible_mod_two();
    assert_eq!(t.target(), &poly(&[0, 9]));
    for first in [0, 1, 2, 3] {
        let root = verified(run(&t, limits(first, 2, 1000)));
        assert_eq!(root.method(), &direct(5, 3));
        assert_root_up_to_sign(&root, &[3, 3]);
        assert_image_is_one_of(&root, [6, 1]);
        assert_agrees_with_square(&root, &t, &rat_poly(&[(1, 1), (-1, 1)]));
    }
    assert_eq!(
        exhausted(run(&t, limits(2, 1, 1000))),
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 1,
            limit: 1
        }]
    );
}

#[test]
fn kat_primes_dividing_the_modulus_or_leading_coefficient_are_excluded() {
    // For d ≥ 2 a prime dividing N or c leaves F with a root, so only d = 1 (where F = t − cm is
    // irreducible at every prime) isolates the divisibility exclusion.
    //
    // f = x − 5, m = 5, N = 21 = 3·7, relation (169, 0): H = 169 = 13², B = 13, minimum 27.
    // The prime 3 divides N and is skipped (it would otherwise succeed: 169 ≡ 1). Next 5:
    // 5² = 25 < 27 ≤ 125 gives e = 3. A budget of 1 is spent on 3 alone.
    let by_modulus = target(&[-5, 1], 5, 21, &[(169, 0)]);
    let root = verified(run(&by_modulus, limits(0, 2, 1000)));
    assert_eq!(root.method(), &direct(5, 3));
    assert_root_up_to_sign(&root, &[13]);
    assert_eq!(
        exhausted(run(&by_modulus, limits(0, 1, 1000))),
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 1,
            limit: 1
        }]
    );

    // f = 3x − 15, m = 5, N = 7 (c = 3): θ = 3α, F = t − 15, and c^(2h−k) = 3 makes
    // η = 3·(3·169) = 1521 = 39², so z = ±39 and B = 39, minimum 79. The prime 3 divides c; at 5,
    // H = 1521 ≡ 1 ≠ 0 and 5² = 25 < 79 ≤ 125 gives e = 3.
    let by_leading = target(&[-15, 3], 5, 7, &[(169, 0)]);
    assert_eq!(by_leading.target(), &poly(&[1521]));
    let root = verified(run(&by_leading, limits(0, 2, 1000)));
    assert_eq!(root.method(), &direct(5, 3));
    assert_root_up_to_sign(&root, &[39]);
    assert_eq!(
        exhausted(run(&by_leading, limits(0, 1, 1000))),
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 1,
            limit: 1
        }]
    );
}

#[test]
fn kat_zero_residue_target_prime_is_skipped() {
    // f = x² − 2, m = 3, N = 7, relation (3, 0) twice: γ = η = 9, F' = 2t, F'² = 8 and H = 72.
    // Mod 3, t² − 2 is irreducible but H ≡ 0, so 3 is skipped before any square root is taken
    // (the zero element has no inverse square root). Mod 5, 2 is a nonresidue so F is irreducible
    // and H ≡ 2 ≠ 0. z = F'·3 = 6θ, with z² = 36·2 = 72, sits exactly on its bound B_1 = 6
    // (E = 9, U = 3, R = 3, B = (18, 6), minimum 37). 5² = 25 < 37 ≤ 125 gives e = 3.
    let t = vanishing_mod_three();
    assert_eq!(t.target(), &poly(&[72]));
    assert_eq!(t.coefficient_bounds(), &[big(18), big(6)]);
    let root = verified(run(&t, limits(0, 2, 1000)));
    assert_eq!(root.method(), &direct(5, 3));
    assert_root_up_to_sign(&root, &[0, 6]);
    assert_eq!(
        exhausted(run(&t, limits(0, 1, 1000))),
        vec![ResourceExhaustion {
            kind: ResourceKind::InertPrimeSearch,
            consumed: 1,
            limit: 1
        }]
    );
    // Starting above 3 never sees the skipped prime.
    let direct_start = verified(run(&t, limits(4, 1, 1000)));
    assert_eq!(direct_start.method(), &direct(5, 3));
}

#[test]
fn kat_reconstruction_agrees_with_an_independently_built_square() {
    // The element β₀ = 1 + θ is squared here by the test's own rational arithmetic and reduced
    // modulo t² − 2 (3 + 2θ), then matched against the target built from the relation (3, −2),
    // whose γ is that square. The reconstructed β is ±β₀, y² ≡ P (mod 527), and the rational root
    // √49 = 7 with y = 24 splits N = 527: gcd(7 − 24, 527) = 17 and gcd(7 + 24, 527) = 31.
    let t = monic_quadratic();
    let beta0 = rat_poly(&[(1, 1), (1, 1)]);
    let root = verified(run(&t, generous()));
    assert_agrees_with_square(&root, &t, &beta0);
    assert_eq!(
        &beta0
            .mul(&beta0)
            .rem(&t.context().monic_poly().to_rat_poly()),
        t.gamma()
    );
    let y = root.modular_image();
    assert!((y * y - root.rational_product()).is_multiple_of(&big(527)));
    let divisors = (
        big(7).gcd(&big(527)),
        (big(7) - y).gcd(&big(527)),
        (big(7) + y).gcd(&big(527)),
    );
    assert_eq!(divisors.0, big(1));
    let found = [divisors.1, divisors.2];
    assert!(found.contains(&big(17)) && found.contains(&big(31)));
}
