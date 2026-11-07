//! Known-answer tests for the normalized target, its certified bounds, and the exact root
//! verifier.
//!
//! Every fixture states the arithmetic that justifies it, derived by hand and not by running the
//! code under test. Throughout, `θ = c·α` is the root of the monic `F(t) = c^(d−1) f(t/c)`,
//! `γ = Π(a − b·α)`, `η = c^(2h−k) Π(c·a − b·θ)` with `h = ⌈k/2⌉`, `H = F'(θ)² η`, and a root
//! of the target is `z = c^h F'(θ) β` with `β² = γ`. Candidate roots `z` are supplied directly:
//! the verifier is exercised with no reconstruction algorithm present.

use gnfs::sqrt::{
    CandidateRejection, CandidateVerdict, FailureClass, IncidentalSource, InputRejection,
    NormalizedTarget, PrecisionLimits, ReconstructionMethod, RelationFactors, ResourceExhaustion,
    ResourceKind, SqrtContext, TargetOutcome, UnsupportedDomain, VerifiedRoot,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use shared_numfield::{IntPoly, RatPoly};

const BIG_BUDGET: u64 = 4096;

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Integer polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

/// Rational polynomial from `(numerator, denominator)` pairs, lowest degree first.
fn rat_poly(coeffs: &[(i64, i64)]) -> RatPoly {
    RatPoly::from_coeffs(coeffs.iter().map(|&(n, d)| BigRational::new(big(n), big(d))).collect())
}

fn bigs(values: &[i64]) -> Vec<BigInt> {
    values.iter().map(|&v| big(v)).collect()
}

fn method() -> ReconstructionMethod {
    ReconstructionMethod::DirectNewton { prime: 3, exponent: 4 }
}

fn outcome(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)], max_bits: u64) -> TargetOutcome {
    let ctx = SqrtContext::new(&poly(f), &big(m), &big(n)).expect("fixture is a valid context");
    let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
    let factors = RelationFactors::new(&ctx, &pairs).expect("fixture relations are units");
    NormalizedTarget::new(factors, PrecisionLimits { max_modulus_bits: max_bits })
        .expect("a validated target builds")
}

fn target(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> NormalizedTarget {
    match outcome(f, m, n, pairs, BIG_BUDGET) {
        TargetOutcome::Ready(target) => target,
        TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
    }
}

fn verdict(target: &NormalizedTarget, z: &[i64]) -> CandidateVerdict {
    target.verify_candidate(&poly(z), method()).expect("no invariant violation")
}

fn verified(target: &NormalizedTarget, z: &[i64]) -> VerifiedRoot {
    match verdict(target, z) {
        CandidateVerdict::Verified(root) => root,
        other => panic!("expected a verified root for {z:?}, got {other:?}"),
    }
}

fn rejected(target: &NormalizedTarget, z: &[i64]) -> CandidateRejection {
    match verdict(target, z) {
        CandidateVerdict::Rejected(why) => why,
        other => panic!("expected a rejection for {z:?}, got {other:?}"),
    }
}

/// `2x² − 1`, `m = 2`, `N = 7`, relation `(3, −4)`: odd count, non-monic.
fn non_monic() -> NormalizedTarget {
    target(&[-1, 0, 2], 2, 7, &[(3, -4)])
}

#[test]
fn kat_bound_formula_hand_computed() {
    // Case 1: F = t² − 8 (f = x² − 8, m = 6, N = 7: f(6) = 28 = 4·7), relation (2, 0), so
    // η = 2. R = 1 + max(|−8|, |0|) = 9. E = |η_0|·R⁰ = 2, U = ⌈√2⌉ = 2 (not a square: the
    // ceiling, not the floor 1). With d = 2 and F_2 = 1, F_1 = 0:
    //   B_0 = d·U·(|F_1|·R⁰ + |F_2|·R¹) = 2·2·9 = 36,   B_1 = d·U·(|F_2|·R⁰) = 2·2·1 = 4.
    // So B = 36 and the least modulus above 2B = 72 is 73.
    let t = target(&[-8, 0, 1], 6, 7, &[(2, 0)]);
    assert_eq!(t.eta(), &poly(&[2]));
    assert_eq!(t.coefficient_bounds(), bigs(&[36, 4]).as_slice());
    assert_eq!(t.bound(), &big(36));
    assert_eq!(t.minimum_modulus(), &big(73));

    // Case 2: F = t² − 2 (m = 3, N = 7: f(3) = 7), relation (9, 0): η = 9, R = 1 + 2 = 3,
    // E = 9, U = 3 (exact square: no round-up). B_0 = 2·3·(0 + 1·3¹) = 18, B_1 = 2·3·1 = 6.
    let t = target(&[-2, 0, 1], 3, 7, &[(9, 0)]);
    assert_eq!(t.eta(), &poly(&[9]));
    assert_eq!(t.coefficient_bounds(), bigs(&[18, 6]).as_slice());
    assert_eq!(t.bound(), &big(18));
    assert_eq!(t.minimum_modulus(), &big(37));

    // Case 3: F = t³ − 2t + 3 (m = 2, N = 7: f(2) = 8 − 4 + 3 = 7), relation (3, 1): η = 3 − t.
    // F = [3, −2, 0, 1], so R = 1 + max(3, 2, 0) = 4. E = 3·4⁰ + 1·4¹ + 0·4² = 7, U = ⌈√7⌉ = 3
    // (the floor would give 2). With d = 3 the common factor is d·U = 9:
    //   i = 0: |F_1|·4⁰ + |F_2|·4¹ + |F_3|·4² = 2 + 0 + 16 = 18  → B_0 = 162
    //   i = 1: |F_2|·4⁰ + |F_3|·4¹            = 0 + 4      = 4   → B_1 = 36
    //   i = 2: |F_3|·4⁰                                   = 1   → B_2 = 9.
    let t = target(&[3, -2, 0, 1], 2, 7, &[(3, 1)]);
    assert_eq!(t.eta(), &poly(&[3, -1]));
    assert_eq!(t.coefficient_bounds(), bigs(&[162, 36, 9]).as_slice());
    assert_eq!(t.bound(), &big(162));
    assert_eq!(t.minimum_modulus(), &big(325));

    // Case 4: an even count with a larger E. F = t² − 2, relations (3, −2) twice: η =
    // (3 + 2t)² = 9 + 12t + 4t² = 17 + 12t, R = 3, E = 17 + 12·3 = 53, U = ⌈√53⌉ = 8
    // (7² = 49 < 53 ≤ 64 = 8²). B_0 = 2·8·3 = 48, B_1 = 2·8·1 = 16.
    let t = target(&[-2, 0, 1], 23, 527, &[(3, -2), (3, -2)]);
    assert_eq!(t.eta(), &poly(&[17, 12]));
    assert_eq!(t.coefficient_bounds(), bigs(&[48, 16]).as_slice());
    assert_eq!(t.bound(), &big(48));
}

#[test]
fn kat_target_identities_hold() {
    // f = 2x² − 1, m = 2, N = 7, relation (3, −4): c = 2, F = t² − 2, k = 1, h = 1.
    // c·a − b·θ = 6 + 4θ, so γ = (6 + 4θ)/2 = 3 + 2θ and η = c^(2h−k)(6 + 4θ) = 2(6 + 4θ)
    // = 12 + 8θ = c^(2h)·γ = 4(3 + 2θ). F' = 2t, and F'² = 4t² ≡ 8, so H = 8η = 96 + 64θ.
    // The inverse of F' is θ/4, since 2θ·(θ/4) = 2θ²/4 = 4/4 = 1.
    // P = a − b·m = 3 + 8 = 11.
    let t = non_monic();
    assert_eq!(t.relation_count(), 1);
    assert_eq!(t.half_count(), 1);
    assert_eq!(t.rational_product(), &big(11));
    assert_eq!(t.gamma(), &rat_poly(&[(3, 1), (2, 1)]));
    assert_eq!(t.eta(), &poly(&[12, 8]));
    assert_eq!(t.derivative(), &poly(&[0, 2]));
    assert_eq!(t.derivative_inverse(), &rat_poly(&[(0, 1), (1, 4)]));
    assert_eq!(t.target(), &poly(&[96, 64]));
    // Bounds: R = 3, E = 12 + 8·3 = 36, U = 6 (a perfect square), B_0 = 2·6·3 = 36, B_1 = 12.
    assert_eq!(t.coefficient_bounds(), bigs(&[36, 12]).as_slice());
    assert_eq!(t.minimum_modulus(), &big(73));

    // The inverse identity, recomputed here from the public accessors on a cubic as well.
    for t in [non_monic(), target(&[3, -2, 0, 1], 2, 7, &[(3, 1)])] {
        let modulus = t.context().monic_poly().to_rat_poly();
        let product = t.derivative().to_rat_poly().mul(t.derivative_inverse()).rem(&modulus);
        assert_eq!(product, RatPoly::one());
    }
}

#[test]
fn kat_odd_cardinality_non_monic_accepts() {
    // f = 2x² − 1, m = 2, N = 7, relation (3, −4); F = t² − 2, θ ↦ c·m = 4. γ = 3 + 2θ =
    // (1 + θ)² since (1 + θ)² = 1 + 2θ + θ² = 3 + 2θ. So β = 1 + θ, and z = c^h F' β =
    // 2·2θ·(1 + θ) = 4θ + 4θ² = 8 + 4θ. Check z² = 64 + 64θ + 16θ² = 96 + 64θ = H.
    // y = β(4) = 5 and 5² = 25 ≡ 4 ≡ 11 = P (mod 7). P = 11 is not an integer square, so this is
    // a test of the algebraic map only.
    let t = non_monic();
    let root = verified(&t, &[8, 4]);
    assert_eq!(root.integral_root(), &poly(&[8, 4]));
    assert_eq!(root.beta(), &rat_poly(&[(1, 1), (1, 1)]));
    assert_eq!(root.modular_image(), &big(5));
    assert_eq!(root.rational_product(), &big(11));
    assert_eq!(root.relation_count(), 1);
    assert_eq!(root.method(), &method());
    assert_eq!(root.context(), t.context());
    let square_gap = root.modular_image() * root.modular_image() - root.rational_product();
    assert!(square_gap.is_multiple_of(&big(7)));

    // The other sign: z = −8 − 4θ gives β = −1 − θ and y = −5 ≡ 2 (mod 7), and 2² = 4 ≡ 11.
    let negated = verified(&t, &[-8, -4]);
    assert_eq!(negated.beta(), &rat_poly(&[(-1, 1), (-1, 1)]));
    assert_eq!(negated.modular_image(), &big(2));
}

#[test]
fn kat_even_cardinality_and_assembly_seam() {
    // f = x² − 2, m = 23, N = 527 = 17·31, relation (3, −2): γ = 3 + 2θ = (1 + θ)², η = 3 + 2θ
    // (c = 1), H = F'²η = 8(3 + 2θ) = 24 + 16θ and z = 2θ(1 + θ) = 4 + 2θ with z² = 16 + 16θ + 8
    // = 24 + 16θ. y = β(23) = 24, P = 3 + 46 = 49, and 24² = 576 = 527 + 49.
    let t = target(&[-2, 0, 1], 23, 527, &[(3, -2)]);
    let root = verified(&t, &[4, 2]);
    assert_eq!(root.modular_image(), &big(24));
    assert_eq!(root.rational_product(), &big(49));
    // The rational root is √49 = 7: gcd(7 − 24, 527) = gcd(17, 527) = 17 and gcd(7 + 24, 527) =
    // gcd(31, 527) = 31, each a proper factor.
    assert_eq!(big(7 - 24).gcd(&big(527)), big(17));
    assert_eq!(big(7 + 24).gcd(&big(527)), big(31));

    // Even k = 2: the same relation twice. Now h = 1, c^(2h−k) = 1, γ = (3 + 2θ)² = 17 + 12θ,
    // β = 3 + 2θ, z = F'β = 2θ(3 + 2θ) = 8 + 6θ and z² = 64 + 96θ + 72 = 136 + 96θ =
    // 8(17 + 12θ) = H. y = 3 + 2·23 = 49, P = 49² = 2401, y² = 2401 exactly.
    let t = target(&[-2, 0, 1], 23, 527, &[(3, -2), (3, -2)]);
    assert_eq!(t.relation_count(), 2);
    assert_eq!(t.half_count(), 1);
    assert_eq!(t.target(), &poly(&[136, 96]));
    let root = verified(&t, &[8, 6]);
    assert_eq!(root.beta(), &rat_poly(&[(3, 1), (2, 1)]));
    assert_eq!(root.modular_image(), &big(49));
    assert_eq!(root.rational_product(), &big(2401));
    assert_eq!(root.relation_count(), 2);
}

#[test]
fn kat_negative_leading_coefficient_accepts() {
    // f = −2x² + 1, m = 2, N = 7: f(2) = −7. c = −2, F = t² − 2, θ ↦ c·m = −4 ≡ 3 (mod 7).
    // Relation (3, −4): α = θ/c = −θ/2, so a − bα = 3 + 4α = 3 − 2θ = γ = (1 − θ)² since
    // (1 − θ)² = 1 − 2θ + 2. η = c^(2h−k)(c·a − b·θ) = (−2)(−6 + 4θ) = 12 − 8θ, H = 8η =
    // 96 − 64θ. z = c^h F' β = (−2)(2θ)(1 − θ) = −4θ + 4θ² = 8 − 4θ; z² = 64 − 64θ + 32 = H.
    // y = β(3) = 1 − 3 = −2 ≡ 5, 5² = 25 ≡ 4 ≡ 11 = P (mod 7).
    let t = target(&[1, 0, -2], 2, 7, &[(3, -4)]);
    assert_eq!(t.gamma(), &rat_poly(&[(3, 1), (-2, 1)]));
    assert_eq!(t.eta(), &poly(&[12, -8]));
    assert_eq!(t.target(), &poly(&[96, -64]));
    let root = verified(&t, &[8, -4]);
    assert_eq!(root.beta(), &rat_poly(&[(1, 1), (-1, 1)]));
    assert_eq!(root.modular_image(), &big(5));
    assert_eq!(root.rational_product(), &big(11));
}

#[test]
fn kat_degree_one_target() {
    // f = x − 5, m = 5, N = 21: F = t − 5, F' = 1. Relation (9, 1): γ = 9 − 5 = 4 = P (α = 5).
    // η = 4 = H, R = 1 + 5 = 6, E = 4, U = 2, B_0 = 1·2·|F_1|·R⁰ = 2. z = 2 sits exactly on the
    // bound; y = 2 and 2² = 4 = P. z = −2 gives y = 19.
    let t = target(&[-5, 1], 5, 21, &[(9, 1)]);
    assert_eq!(t.coefficient_bounds(), bigs(&[2]).as_slice());
    assert_eq!(t.derivative_inverse(), &rat_poly(&[(1, 1)]));
    assert_eq!(verified(&t, &[2]).modular_image(), &big(2));
    assert_eq!(verified(&t, &[-2]).modular_image(), &big(19));
    assert_eq!(rejected(&t, &[3]), CandidateRejection::CoefficientOutOfBound { index: 0 });
    assert_eq!(rejected(&t, &[1]), CandidateRejection::NotSquareRootOfTarget);
    assert_eq!(rejected(&t, &[0, 1]), CandidateRejection::DegreeTooHigh);

    // Relation (7, 1): γ = 7 − 5 = 2 is not a rational square: no z in [−2, 2] squares to H = 2.
    let t = target(&[-5, 1], 5, 21, &[(7, 1)]);
    for z in [-2, -1, 0, 1, 2] {
        assert_eq!(rejected(&t, &[z]), CandidateRejection::NotSquareRootOfTarget);
    }
}

#[test]
fn kat_outside_order_root_accepts() {
    // F = t² − 8 (m = 6, N = 7), relation (2, 0): γ = 2, η = 2, and β = θ/2 satisfies
    // β² = θ²/4 = 8/4 = 2, but β ∉ Z[θ]. A search for roots in Z[θ] alone would miss it; the
    // integral root is z = F'β = 2θ·θ/2 = θ² = 8, with z² = 64 = F'²η = 4·8·2 = H.
    // y = β(6) = 6/2 = 3, and 3² = 9 ≡ 2 = P (mod 7).
    let t = target(&[-8, 0, 1], 6, 7, &[(2, 0)]);
    assert_eq!(t.target(), &poly(&[64]));
    let root = verified(&t, &[8]);
    assert_eq!(root.beta(), &rat_poly(&[(0, 1), (1, 2)]));
    assert_eq!(root.modular_image(), &big(3));
    assert_eq!(root.rational_product(), &big(2));
    assert_eq!(verified(&t, &[-8]).modular_image(), &big(4));
}

#[test]
fn kat_zero_divisor_factors_are_rejected_before_a_target_exists() {
    // x² − 1 at m = 4, N = 15: f(4) = 15. f_hom(a, b) = a² − b². (1, 1) and (1, −1) make it 0,
    // so 1 ∓ θ kills the root ±1 of F = t² − 1: a nonzero zero divisor, not a unit. The boundary
    // type refuses them, so no NormalizedTarget (and no prime search) can follow.
    let ctx = SqrtContext::new(&poly(&[-1, 0, 1]), &big(4), &big(15)).unwrap();
    for pair in [(1, 1), (1, -1)] {
        let pairs = [(big(2), big(0)), (big(pair.0), big(pair.1))];
        let rejection = RelationFactors::new(&ctx, &pairs).unwrap_err();
        let expected = UnsupportedDomain::ZeroDivisorFactor { index: 1 };
        assert_eq!(rejection, InputRejection::Unsupported(expected));
        assert_eq!(rejection.class(), FailureClass::Unsupported);
    }
}

#[test]
fn kat_wrong_norm_map_counterexample() {
    // f = x² − 2, m = 3, N = 7 (f(3) = 7), relation (9, 0): γ = 9, P = 9 ≡ 2 (mod 7), F' = 2t,
    // η = 9, H = 4t²·9 = 72. β = 3, z = F'β = 6θ with z² = 36·2 = 72. The image is y = β(3) = 3
    // and 3² = 9 ≡ 2 ≡ P. The legacy return, the norm Norm(β) = β·β̄ = 3·3 = 9 ≡ 2, squares to
    // 4 ≢ 2 ≡ P: the norm is not the map.
    let t = target(&[-2, 0, 1], 3, 7, &[(9, 0)]);
    assert_eq!(t.target(), &poly(&[72]));
    let root = verified(&t, &[0, 6]);
    assert_eq!(root.beta(), &rat_poly(&[(3, 1)]));
    assert_eq!(root.modular_image(), &big(3));
    let n = big(7);
    let p_mod_n = root.rational_product().mod_floor(&n);
    assert_eq!((root.modular_image() * root.modular_image()).mod_floor(&n), p_mod_n);
    let norm_of_beta = big(9);
    assert_ne!((&norm_of_beta * &norm_of_beta).mod_floor(&n), p_mod_n);

    // A candidate built from the wrong equation β² = Norm(γ) = 81 has β = 9 and z = 18θ, whose
    // coefficient 18 exceeds B_1 = 2·3·1 = 6. The candidate built from y = 2 ≡ Norm(β) has
    // β = 2 and z = 4θ, within bounds but with z² = 32 ≠ 72.
    assert_eq!(rejected(&t, &[0, 18]), CandidateRejection::CoefficientOutOfBound { index: 1 });
    assert_eq!(rejected(&t, &[0, 4]), CandidateRejection::NotSquareRootOfTarget);
}

#[test]
fn kat_noninvertible_denominators_have_typed_outcomes() {
    // F = t² − 8, relation (2, 0), z = 8, β = θ/2: the only denominator is 2.
    // N = 4 (m = 2, f(2) = −4 = −1·4): gcd(2, 4) = 2 is a proper divisor, so the result is an
    // incidental factor 4 = 2·2 found in the map denominator.
    let t = target(&[-8, 0, 1], 2, 4, &[(2, 0)]);
    match verdict(&t, &[8]) {
        CandidateVerdict::Incidental(found) => {
            assert_eq!(found.divisor(), &big(2));
            assert_eq!(found.cofactor(), &big(2));
            assert_eq!(found.source(), IncidentalSource::MapDenominator);
        }
        other => panic!("expected an incidental factor, got {other:?}"),
    }

    // N = 2 (m = 0, f(0) = −8): gcd(2, 2) = N, so the map is unsupported, not a factor.
    let t = target(&[-8, 0, 1], 0, 2, &[(2, 0)]);
    assert_eq!(
        verdict(&t, &[8]),
        CandidateVerdict::Unsupported(UnsupportedDomain::MapDenominatorNotUnit)
    );

    // The same root with N = 7 has a unit denominator: it is the denominator, not the root,
    // that blocked the two cases above.
    let t = target(&[-8, 0, 1], 6, 7, &[(2, 0)]);
    assert!(matches!(verdict(&t, &[8]), CandidateVerdict::Verified(_)));
}

#[test]
fn kat_candidate_rejections_are_classified_in_order() {
    // Target of the non-monic fixture: H = 96 + 64θ, bounds B = [36, 12], d = 2.
    let t = non_monic();
    // Degree is checked before bounds: z = 1000θ² has degree 2 = d.
    assert_eq!(rejected(&t, &[0, 0, 1000]), CandidateRejection::DegreeTooHigh);
    // Bounds come next, reporting the first offending index; the sign is irrelevant.
    assert_eq!(rejected(&t, &[37]), CandidateRejection::CoefficientOutOfBound { index: 0 });
    assert_eq!(rejected(&t, &[-37]), CandidateRejection::CoefficientOutOfBound { index: 0 });
    assert_eq!(rejected(&t, &[0, 13]), CandidateRejection::CoefficientOutOfBound { index: 1 });
    assert_eq!(rejected(&t, &[37, 13]), CandidateRejection::CoefficientOutOfBound { index: 0 });
    // Within bounds but not roots: zero; (8 + 5θ)² = 64 + 80θ + 50 = 114 + 80θ; and the sign
    // slip 8 − 4θ, whose square 64 − 64θ + 32 = 96 − 64θ differs from H = 96 + 64θ.
    for z in [&[][..], &[8, 5], &[8, -4]] {
        assert_eq!(rejected(&t, z), CandidateRejection::NotSquareRootOfTarget);
    }
    // The boundary itself is inclusive: |z_i| = B_i is within the bound and rejected only by the
    // square test. z = 36 + 12θ: z² = 1296 + 864θ + 288 = 1584 + 864θ ≠ H.
    assert_eq!(rejected(&t, &[36, 12]), CandidateRejection::NotSquareRootOfTarget);
}

#[test]
fn kat_modulus_precheck_precedes_any_search() {
    // The non-monic fixture has B = 36, so M > 72 forces M ≥ 73, which has 7 bits
    // (64 ≤ 73 < 128). The budget is on bits(2B + 1): 7 is enough, 6 is not.
    let spec = (&[-1, 0, 2][..], 2, 7, &[(3, -4)][..]);
    assert!(matches!(outcome(spec.0, spec.1, spec.2, spec.3, 7), TargetOutcome::Ready(_)));
    for limit in [6, 0] {
        assert_eq!(
            outcome(spec.0, spec.1, spec.2, spec.3, limit),
            TargetOutcome::ResourceExhausted(ResourceExhaustion {
                kind: ResourceKind::Precision,
                consumed: 7,
                limit,
            })
        );
    }

    // The cubic has B = 162, minimum modulus 325, with 9 bits (256 ≤ 325 < 512).
    let cubic = |bits| outcome(&[3, -2, 0, 1], 2, 7, &[(3, 1)], bits);
    assert!(matches!(cubic(9), TargetOutcome::Ready(_)));
    assert_eq!(
        cubic(8),
        TargetOutcome::ResourceExhausted(ResourceExhaustion {
            kind: ResourceKind::Precision,
            consumed: 9,
            limit: 8,
        })
    );
}
