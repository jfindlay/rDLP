//! Known-answer tests for the checked square-root context, the relation-factor boundary, and the
//! failure-class vocabulary.
//!
//! Every fixture is justified by arithmetic stated beside it. Throughout, `f(m) ≡ 0 (mod N)` is
//! verified by hand, `F(t) = c^(d-1) f(t/c)` has coefficient `f_k · c^(d-1-k)` below the leading
//! term, and a relation factor `a - b·α` is a unit of `Q[t]/(F)` exactly when
//! `f_hom(a, b) = Σ f_i a^i b^(d-i)` is nonzero.

use std::collections::HashSet;

use gnfs::polyselect::PolyPair;
use gnfs::sqrt::{
    CrtLimits, DirectLimits, FailureClass, IncidentalFactor, IncidentalSource, InputRejection,
    InvalidInput, InvariantViolation, PrecisionLimits, PrimeSearchLimits, RecombinationLimits,
    RelationFactors, ResourceExhaustion, ResourceKind, SqrtContext, SqrtError, SqrtLimits,
    UnsupportedDomain, AlgebraicSqrtOutcome,
};
use num_bigint::BigInt;
use shared_numfield::{IntPoly, NonresidueSearchLimits, RootSearchLimits};

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

fn pairs(list: &[(i64, i64)]) -> Vec<(BigInt, BigInt)> {
    list.iter().map(|&(a, b)| (big(a), big(b))).collect()
}

fn invalid(e: InvalidInput) -> InputRejection {
    InputRejection::Invalid(e)
}

fn context(f: &[i64], m: i64, n: i64) -> SqrtContext {
    SqrtContext::new(&poly(f), &big(m), &big(n)).expect("fixture is a valid context")
}

#[test]
fn kat_context_accepts_monic_and_non_monic() {
    // x² − 2, m = 23, N = 527 = 17·31: f(23) = 529 − 2 = 527 = 1·N. Monic, so F = f; θ = α maps
    // to 23. disc(t² − 2) = b² − 4ac = 0 + 8 = 8.
    let monic = context(&[-2, 0, 1], 23, 527);
    assert_eq!(monic.modulus(), &big(527));
    assert_eq!(monic.root(), &big(23));
    assert_eq!(monic.algebraic_poly(), &poly(&[-2, 0, 1]));
    assert_eq!(monic.leading_coeff(), &big(1));
    assert_eq!(monic.monic_poly(), &poly(&[-2, 0, 1]));
    assert_eq!(monic.degree(), 2);
    assert_eq!(monic.theta_image(), &big(23));
    assert_eq!(monic.discriminant(), &big(8));

    // 2x² − 1, m = 2, N = 7: f(2) = 8 − 1 = 7. c = 2, so F = t² + 0·t + (−1)·2 = t² − 2, and
    // θ = 2α maps to c·m = 4. F(4) = 14 ≡ 0 (mod 7), as the homomorphism requires.
    let non_monic = context(&[-1, 0, 2], 2, 7);
    assert_eq!(non_monic.algebraic_poly(), &poly(&[-1, 0, 2]));
    assert_eq!(non_monic.leading_coeff(), &big(2));
    assert_eq!(non_monic.monic_poly(), &poly(&[-2, 0, 1]));
    assert_eq!(non_monic.theta_image(), &big(4));
    assert_eq!(non_monic.discriminant(), &big(8));

    // −2x² + 1, m = 2, N = 7: f(2) = −8 + 1 = −7. The sign of f is kept (c = −2); F is still
    // t² − 2, and c·m = −4 reduces to 3 in [0, 7).
    let negative = context(&[1, 0, -2], 2, 7);
    assert_eq!(negative.leading_coeff(), &big(-2));
    assert_eq!(negative.monic_poly(), &poly(&[-2, 0, 1]));
    assert_eq!(negative.theta_image(), &big(3));

    // m is kept exact: 550 = 23 + 527 gives f(550) = 302500 − 2 = 302498 = 574·527, and the
    // image c·m mod N is still 23.
    let shifted = context(&[-2, 0, 1], 550, 527);
    assert_eq!(shifted.root(), &big(550));
    assert_eq!(shifted.theta_image(), &big(23));
}

#[test]
fn kat_context_rejects_modulus_le_one() {
    // Modulus checks come first and never reach a division: N = 1, 0, −5 are all rejected, even
    // with a zero polynomial that would otherwise fail a later check.
    for n in [1, 0, -5] {
        let expected = invalid(InvalidInput::ModulusNotGreaterThanOne { n: big(n) });
        assert_eq!(SqrtContext::new(&poly(&[-2, 0, 1]), &big(23), &big(n)), Err(expected.clone()));
        assert_eq!(SqrtContext::new(&IntPoly::zero(), &big(0), &big(n)), Err(expected));
    }
}

#[test]
fn kat_context_rejects_constant_and_zero_polynomial() {
    for f in [IntPoly::zero(), poly(&[5])] {
        assert_eq!(
            SqrtContext::new(&f, &big(3), &big(7)),
            Err(invalid(InvalidInput::ConstantPolynomial))
        );
    }
}

#[test]
fn kat_context_rejects_root_mismatch() {
    // x² − 2, m = 6, N = 35: f(6) = 36 − 2 = 34, and 34 mod 35 = 34 ≠ 0.
    assert_eq!(
        SqrtContext::new(&poly(&[-2, 0, 1]), &big(6), &big(35)),
        Err(invalid(InvalidInput::RootMismatch { f_of_m_mod_n: big(34) }))
    );
    // Negative values reduce to the least nonnegative residue: x² − 2, m = 1, N = 7 gives
    // f(1) = −1 ≡ 6.
    assert_eq!(
        SqrtContext::new(&poly(&[-2, 0, 1]), &big(1), &big(7)),
        Err(invalid(InvalidInput::RootMismatch { f_of_m_mod_n: big(6) }))
    );
}

#[test]
fn kat_from_poly_pair_rejects_degree_field_and_rational_poly() {
    let pair = |f: &[i64], g: &[i64], m: i64, n: i64, degree: usize| PolyPair {
        f: poly(f),
        g: poly(g),
        m: big(m),
        n: big(n),
        degree,
        skew: None,
        factor_base_bounds: None,
    };

    // The valid pair x² − 2, g = x − 23, m = 23, N = 527 yields the same context as `new`.
    let good = pair(&[-2, 0, 1], &[-23, 1], 23, 527, 2);
    assert_eq!(SqrtContext::from_poly_pair(&good), Ok(context(&[-2, 0, 1], 23, 527)));

    // The declared degree 3 disagrees with deg f = 2.
    let wrong_degree = pair(&[-2, 0, 1], &[-23, 1], 23, 527, 3);
    assert_eq!(
        SqrtContext::from_poly_pair(&wrong_degree),
        Err(invalid(InvalidInput::DegreeFieldMismatch { field_degree: 3, poly_degree: 2 }))
    );

    // g = x − 5 is not x − m for m = 23.
    let wrong_g = pair(&[-2, 0, 1], &[-5, 1], 23, 527, 2);
    assert_eq!(
        SqrtContext::from_poly_pair(&wrong_g),
        Err(invalid(InvalidInput::RationalPolyNotLinear {
            expected: poly(&[-23, 1]),
            actual: poly(&[-5, 1]),
        }))
    );

    // N = 0 is reported first, even though g and degree are also wrong, and nothing panics.
    let zero_n = pair(&[-2, 0, 1], &[-5, 1], 23, 0, 7);
    assert_eq!(
        SqrtContext::from_poly_pair(&zero_n),
        Err(invalid(InvalidInput::ModulusNotGreaterThanOne { n: big(0) }))
    );

    // A constant f is rejected before the degree field is compared.
    let constant = pair(&[5], &[-3, 1], 3, 7, 4);
    assert_eq!(
        SqrtContext::from_poly_pair(&constant),
        Err(invalid(InvalidInput::ConstantPolynomial))
    );

    // Later checks are those of `new`: x² − 2 with m = 6, N = 35 has f(m) = 34.
    let bad_root = pair(&[-2, 0, 1], &[-6, 1], 6, 35, 2);
    assert_eq!(
        SqrtContext::from_poly_pair(&bad_root),
        Err(invalid(InvalidInput::RootMismatch { f_of_m_mod_n: big(34) }))
    );
}

#[test]
fn kat_context_rejects_repeated_factor() {
    // (x − 3)² = x² − 6x + 9, m = 3, N = 7: f(3) = 0. F = f, disc = 0.
    assert_eq!(
        SqrtContext::new(&poly(&[9, -6, 1]), &big(3), &big(7)),
        Err(invalid(InvalidInput::RepeatedFactor))
    );
    // (2x − 1)² = 4x² − 4x + 1, m = 4, N = 7: f(4) = 64 − 16 + 1 = 49 = 7·7. c = 4, so
    // F = t² − 4t + 4 = (t − 2)², again with disc 0.
    assert_eq!(
        SqrtContext::new(&poly(&[1, -4, 4]), &big(4), &big(7)),
        Err(invalid(InvalidInput::RepeatedFactor))
    );
}

#[test]
fn kat_context_incidental_leading_coefficient_factor() {
    // 3x² + x − 4 = (x − 1)(3x + 4), m = 16, N = 15: f(16) = 768 + 16 − 4 = 780 = 52·15. F =
    // t² + t − 12 = (t − 3)(t + 4), disc = 1 + 48 = 49 ≠ 0. gcd(3, 15) = 3 is proper, with
    // 3·5 = 15.
    let found = |rejection: Result<SqrtContext, InputRejection>| match rejection {
        Err(InputRejection::Incidental(found)) => found,
        other => panic!("expected an incidental factor, got {other:?}"),
    };
    let factor = found(SqrtContext::new(&poly(&[-4, 1, 3]), &big(16), &big(15)));
    assert_eq!(factor.divisor(), &big(3));
    assert_eq!(factor.cofactor(), &big(5));
    assert_eq!(factor.divisor() * factor.cofactor(), big(15));
    assert_eq!(factor.source(), IncidentalSource::LeadingCoefficient);

    // The negated polynomial −3x² − x + 4 has the same roots; gcd(|−3|, 15) = 3 again.
    let negated = found(SqrtContext::new(&poly(&[4, -1, -3]), &big(16), &big(15)));
    assert_eq!(negated.divisor(), &big(3));
    assert_eq!(negated.cofactor(), &big(5));
}

#[test]
fn kat_context_invalid_precedes_incidental() {
    // Same polynomial as the incidental case but m = 2: f(2) = 12 + 2 − 4 = 10, not 0 mod 15.
    // gcd(3, 15) = 3 is proper, yet the malformed root is reported.
    assert_eq!(
        SqrtContext::new(&poly(&[-4, 1, 3]), &big(2), &big(15)),
        Err(invalid(InvalidInput::RootMismatch { f_of_m_mod_n: big(10) }))
    );
    // 3x², m = 0, N = 15: f(0) = 0. F = t² has a repeated factor, reported before gcd(3, 15).
    assert_eq!(
        SqrtContext::new(&poly(&[0, 0, 3]), &big(0), &big(15)),
        Err(invalid(InvalidInput::RepeatedFactor))
    );
}

#[test]
fn kat_context_leading_coefficient_not_unit_unsupported() {
    // 5x² + x, m = 5, N = 5: f(5) = 125 + 5 = 130 = 26·5. F = t² + t, disc = 1, and
    // gcd(5, 5) = 5 = N, so the unscaled map does not exist.
    assert_eq!(
        SqrtContext::new(&poly(&[0, 1, 5]), &big(5), &big(5)),
        Err(InputRejection::Unsupported(UnsupportedDomain::LeadingCoefficientNotUnit))
    );
}

#[test]
fn kat_context_reducible_squarefree_accepted() {
    // x² − 1 = (x − 1)(x + 1), m = 4, N = 15: f(4) = 15. F = f, disc = 4 ≠ 0. The quotient
    // algebra is a product of two copies of Q, and the context still accepts it.
    let ctx = context(&[-1, 0, 1], 4, 15);
    assert_eq!(ctx.discriminant(), &big(4));
    assert_eq!(ctx.theta_image(), &big(4));
    assert_eq!(ctx.degree(), 2);
}

#[test]
fn kat_relation_factors_boundary() {
    // x² − 1, m = 4, N = 15: f_hom(a, b) = a² − b².
    let ctx = context(&[-1, 0, 1], 4, 15);
    let reject = |list: &[(i64, i64)]| RelationFactors::new(&ctx, &pairs(list));

    assert_eq!(reject(&[]), Err(invalid(InvalidInput::EmptyRelationSet)));

    // (4, 1): a − b·m = 0. Reported at its index, and before the algebraic test, even for the
    // zero pair (0, 0) which is also the zero element.
    assert_eq!(reject(&[(4, 1)]), Err(invalid(InvalidInput::ZeroRationalFactor { index: 0 })));
    assert_eq!(
        reject(&[(2, 1), (0, 0)]),
        Err(invalid(InvalidInput::ZeroRationalFactor { index: 1 }))
    );

    // (1, 1) and (1, −1): a − b·m = −3 and 5, but f_hom = 1 − 1 = 0, a nonzero element
    // 1 − θ (resp. 1 + θ) that vanishes at the root θ = 1 (resp. −1).
    assert_eq!(
        reject(&[(1, 1)]),
        Err(InputRejection::Unsupported(UnsupportedDomain::ZeroDivisorFactor { index: 0 }))
    );
    assert_eq!(
        reject(&[(1, -1)]),
        Err(InputRejection::Unsupported(UnsupportedDomain::ZeroDivisorFactor { index: 0 }))
    );

    // The first failing index is reported: (2, 1) is fine (4 − 1 = 3), (1, 1) fails at index 1,
    // and (0, 0) at index 2 is never reached.
    assert_eq!(
        reject(&[(2, 1), (1, 1), (0, 0)]),
        Err(InputRejection::Unsupported(UnsupportedDomain::ZeroDivisorFactor { index: 1 }))
    );

    // Duplicates, negative b, and b = 0 are accepted: (3, −1) has f_hom = 9 − 1 = 8, (5, 0) has
    // f_hom = 25.
    let accepted = pairs(&[(2, 1), (2, 1), (3, -1), (5, 0)]);
    let factors = RelationFactors::new(&ctx, &accepted).expect("all factors are units");
    assert_eq!(factors.relation_count(), 4);
    assert_eq!(factors.pairs(), accepted.as_slice());
    assert_eq!(factors.context(), &ctx);

    // Degree 1: x − 5, m = 12, N = 7 (f(12) = 7), so A = Q and F = t − 5. The pair (5, 1) has
    // f_hom = f_0·b + f_1·a = −5 + 5 = 0, and a − b·m = −7 ≠ 0: the zero element of A. The pair
    // (6, 1) has f_hom = 1 and is a unit.
    let linear = context(&[-5, 1], 12, 7);
    assert_eq!(linear.monic_poly(), &poly(&[-5, 1]));
    assert_eq!(
        RelationFactors::new(&linear, &pairs(&[(6, 1), (5, 1)])),
        Err(invalid(InvalidInput::ZeroAlgebraicFactor { index: 1 }))
    );
    assert!(RelationFactors::new(&linear, &pairs(&[(6, 1)])).is_ok());
}

#[test]
fn kat_incidental_factor_is_self_certifying() {
    let n = big(15);
    let made = IncidentalFactor::new(&n, big(5), IncidentalSource::MapDenominator)
        .expect("5 is a proper divisor of 15");
    assert_eq!((made.divisor(), made.cofactor()), (&big(5), &big(3)));
    assert_eq!(made.source(), IncidentalSource::MapDenominator);
    // Trivial, non-dividing, and out-of-range divisors are refused.
    for bad in [-3, 0, 1, 15, 30, 4] {
        assert_eq!(IncidentalFactor::new(&n, big(bad), IncidentalSource::MapDenominator), None);
    }
    assert_eq!(IncidentalFactor::new(&big(1), big(1), IncidentalSource::MapDenominator), None);
}

#[test]
fn kat_failure_classes_distinct() {
    let factor = IncidentalFactor::new(&big(15), big(3), IncidentalSource::LeadingCoefficient)
        .expect("proper divisor");
    let exhaustion = ResourceExhaustion { kind: ResourceKind::Precision, consumed: 9, limit: 8 };

    let rejections = [
        (invalid(InvalidInput::RepeatedFactor), FailureClass::InvalidInput),
        (
            InputRejection::Unsupported(UnsupportedDomain::MapDenominatorNotUnit),
            FailureClass::Unsupported,
        ),
        (InputRejection::Incidental(factor.clone()), FailureClass::IncidentalFactor),
    ];
    for (rejection, class) in &rejections {
        assert_eq!(rejection.class(), *class);
    }

    let invariant = SqrtError::Invariant(InvariantViolation { what: "z² ≡ H" });
    let errors = [
        (SqrtError::Invalid(InvalidInput::EmptyRelationSet), FailureClass::InvalidInput),
        (invariant, FailureClass::InvariantViolation),
    ];
    for (error, class) in &errors {
        assert_eq!(error.class(), *class);
    }

    let outcomes = [
        (
            AlgebraicSqrtOutcome::Unsupported(UnsupportedDomain::LeadingCoefficientNotUnit),
            FailureClass::Unsupported,
        ),
        (AlgebraicSqrtOutcome::Incidental(factor), FailureClass::IncidentalFactor),
        (AlgebraicSqrtOutcome::ResourceExhausted(vec![exhaustion]), FailureClass::ResourceExhausted),
    ];
    for (outcome, class) in &outcomes {
        assert_eq!(outcome.class(), Some(*class));
    }

    // Input errors and outcomes share the class vocabulary, and each distinct condition has its
    // own class.
    let classes: HashSet<FailureClass> = rejections
        .iter()
        .map(|(_, c)| *c)
        .chain(errors.iter().map(|(_, c)| *c))
        .chain(outcomes.iter().map(|(_, c)| *c))
        .collect();
    assert_eq!(classes.len(), 5);
}

#[test]
fn errors_display_their_cause() {
    let root = invalid(InvalidInput::RootMismatch { f_of_m_mod_n: big(34) });
    assert!(root.to_string().contains("34"));
    let unit = InputRejection::Unsupported(UnsupportedDomain::ZeroDivisorFactor { index: 3 });
    assert!(unit.to_string().contains('3'));
    let factor = IncidentalFactor::new(&big(15), big(3), IncidentalSource::LeadingCoefficient)
        .expect("proper divisor");
    assert!(InputRejection::Incidental(factor).to_string().contains("15 = 3 * 5"));
    let error = SqrtError::Invariant(InvariantViolation { what: "beta squares to gamma" });
    assert!(error.to_string().contains("beta squares to gamma"));
    let boxed: Box<dyn std::error::Error> = Box::new(SqrtError::Invalid(InvalidInput::RepeatedFactor));
    assert!(boxed.to_string().contains("repeated factor"));
}

#[test]
fn resource_limits_are_plain_copyable_values() {
    let primes = PrimeSearchLimits { first_candidate: 3, max_primes_examined: 100 };
    let nonresidue = NonresidueSearchLimits { max_candidates: 64 };
    let limits = SqrtLimits {
        precision: PrecisionLimits { max_modulus_bits: 4096 },
        direct: DirectLimits { prime_search: primes, nonresidue_search: nonresidue },
        crt: CrtLimits {
            prime_search: primes,
            max_primes: 3,
            root_search: RootSearchLimits { max_evaluations: 1000 },
            nonresidue_search: nonresidue,
        },
        recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 10 },
    };
    let copy = limits;
    assert_eq!(copy, limits);
    assert_eq!(copy.crt.prime_search.max_primes_examined, 100);
}
