//! Known-answer tests for `dedekind_factor` and `dedekind_index_test`: the degree-one primes
//! above p in ℤ[α], their multiplicities, and the index criterion.
//!
//! Tests cover inert, split, ramified and partially split primes, root-search exhaustion, the
//! degree invariant of the report, discriminant computation, bad-prime detection, and the
//! Dedekind index criterion with hand-derived obstruction polynomials.
//!
//! A prime with no roots mod p has no degree-one primes: its report is empty and
//! `nonlinear_degree` is deg f. No ideal stands in for it.

use num_bigint::BigInt;
use shared_numfield::{
    dedekind_factor, dedekind_index_test, discriminant, is_bad_prime, DedekindError,
    DedekindReport, FpPoly, IntPoly, NumberField, PrimeModulus, RootSearchLimits,
};

/// A budget far above any root search over the small primes used here.
const AMPLE: RootSearchLimits = RootSearchLimits { max_evaluations: 1_000 };

/// The (r, multiplicity) pairs of a report, in the order reported.
fn roots_and_multiplicities(report: &DedekindReport<'_>) -> Vec<(BigInt, usize)> {
    report.degree_one.iter().map(|d| (d.ideal.r.clone(), d.multiplicity)).collect()
}

fn bi(n: i64) -> BigInt {
    BigInt::from(n)
}

/// Integer polynomial from least-significant-first coefficients.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| bi(c)).collect())
}

fn prime(p: u64) -> PrimeModulus {
    PrimeModulus::new(p).unwrap()
}

/// The monic polynomial over 𝔽_p with the given least-significant-first coefficients.
fn fp(p: u64, coeffs: &[u64]) -> FpPoly {
    FpPoly::from_coeffs(prime(p), coeffs)
}

/// f = x² − 2 (defines ℚ(√2))
fn field_sqrt2() -> NumberField {
    // coeffs: [-2, 0, 1] → -2 + 0·x + 1·x²
    NumberField::new(IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]))
}

/// f = x³ − x − 1 (discriminant −23, a totally real cubic)
fn field_cubic() -> NumberField {
    // coeffs: [-1, -1, 0, 1] → -1 − x + 0·x² + x³
    NumberField::new(IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)]))
}

// ─── Degree-one prime report ──────────────────────────────────────────────────

/// f = x² − 2, p = 3: f(0) = −2 ≡ 1, f(1) = −1 ≡ 2, f(2) = 2 (mod 3), so f̄ has no roots and is
/// irreducible. (3) is inert: there are no degree-one primes and the whole degree is nonlinear.
/// In particular no ideal such as (3, α) is listed: that ideal is the unit ideal here, since it
/// contains 3 and α² = 2, hence 3 − 2 = 1.
#[test]
fn kat1_inert_prime_has_no_degree_one_primes() {
    let k = field_sqrt2();
    let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();

    assert!(report.degree_one.is_empty(), "x²−2 has no roots mod 3");
    assert_eq!(report.nonlinear_degree, 2);
    assert!(!report.is_bad_prime);
    assert!(!report.index.divides_index());
}

/// f = x² − 2, p = 7: 3² = 9 ≡ 2 and 4² = 16 ≡ 2 (mod 7), so f̄ = (x − 3)(x − 4), two distinct
/// roots of multiplicity 1, listed in ascending order.
#[test]
fn kat2_split_prime() {
    let k = field_sqrt2();
    let report = dedekind_factor(&k, prime(7), AMPLE).unwrap();

    assert_eq!(roots_and_multiplicities(&report), vec![(bi(3), 1), (bi(4), 1)]);
    for d in &report.degree_one {
        assert_eq!(d.ideal.p, bi(7));
    }
    assert_eq!(report.nonlinear_degree, 0);
    assert!(!report.is_bad_prime);
    assert!(!report.index.divides_index());
}

/// f = x³ − x − 1, p = 5: f(0) = −1, f(1) = −1, f(2) = 5 ≡ 0, f(3) = 23 ≡ 3, f(4) = 59 ≡ 4 (mod 5),
/// so the only root is 2, and (x − 2) divides f̄ with quotient x² + 2x + 3. That quadratic has
/// discriminant 4 − 12 = −8 ≡ 2 (mod 5), a non-residue, so it is irreducible: one degree-one
/// prime and nonlinear degree 2. Since 5 ∤ disc = −23, the multiplicity is 1.
#[test]
fn kat3_cubic_partial_split() {
    let k = field_cubic();
    let report = dedekind_factor(&k, prime(5), AMPLE).unwrap();

    assert_eq!(roots_and_multiplicities(&report), vec![(bi(2), 1)]);
    assert_eq!(report.nonlinear_degree, 2);
    assert!(!report.is_bad_prime);
}

/// f = x² − 2, p = 7 (totally split): both degree-one primes have norm 7, so the product of
/// N(𝔭)^e over the report is p^deg f = 49, as the degree invariant requires.
#[test]
fn kat4_norm_product() {
    let k = field_sqrt2();
    let report = dedekind_factor(&k, prime(7), AMPLE).unwrap();

    assert_eq!(report.nonlinear_degree, 0);
    let norm_product: BigInt =
        report.degree_one.iter().map(|d| d.ideal.norm().pow(d.multiplicity as u32)).product();
    assert_eq!(norm_product, bi(7).pow(2u32));
}

// ─── KAT 5 — Discriminant of x² − 2 ─────────────────────────────────────────

/// disc(x² − 2) = 8.
///
/// Derivation: f = x² − 2, f' = 2x.
/// Res(f, f') = Res(x² − 2, 2x).
/// Sylvester matrix (3×3):
///   [ 1   0  −2 ]
///   [ 2   0   0 ]
///   [ 0   2   0 ]
/// det = 1·(0·0 − 0·2) − 0·(2·0 − 0·0) + (−2)·(2·2 − 0·0) = −8.
/// Sign factor: d=2, d(d−1)/2 = 1 (odd) → disc = −(−8) = 8.
#[test]
fn kat5_discriminant_quadratic() {
    let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);
    let disc = discriminant(&f);
    assert_eq!(disc, bi(8), "disc(x²−2) should be 8, got {disc}");
}

// ─── KAT 6 — Discriminant of x³ − x − 1 ─────────────────────────────────────

/// disc(x³ − x − 1) = −23.
///
/// This is a standard reference value: the splitting field of x³ − x − 1 has discriminant
/// −23, which is also the discriminant of the cubic number field ℚ(α) where α³ − α − 1 = 0.
/// The field has class number 1 and is the unique cubic field of discriminant −23.
#[test]
fn kat6_discriminant_cubic() {
    let f = IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)]);
    let disc = discriminant(&f);
    assert_eq!(disc, bi(-23), "disc(x³−x−1) should be −23, got {disc}");
}

// ─── KAT 7 — Bad-prime detection ─────────────────────────────────────────────

/// `is_bad_prime(x² − 2, 2)` is true (2 | disc = 8).
/// `is_bad_prime(x² − 2, 3)` is false (3 ∤ disc = 8).
///
/// The prime 2 is the unique bad prime for x² − 2: it divides the discriminant 8 = 2³.
/// All odd primes are good primes for x² − 2.
#[test]
fn kat7_is_bad_prime() {
    let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);

    assert!(
        is_bad_prime(&f, &bi(2)),
        "p=2 should be a bad prime for x²−2 (2 | disc=8)"
    );
    assert!(
        !is_bad_prime(&f, &bi(3)),
        "p=3 should not be a bad prime for x²−2 (3 ∤ disc=8)"
    );
    assert!(
        !is_bad_prime(&f, &bi(5)),
        "p=5 should not be a bad prime for x²−2 (5 ∤ disc=8)"
    );
    assert!(
        !is_bad_prime(&f, &bi(7)),
        "p=7 should not be a bad prime for x²−2 (7 ∤ disc=8)"
    );
}

// ─── Ramified and bad primes ──────────────────────────────────────────────────

/// f = x² − 2, p = 2: f̄ = x² (since −2 ≡ 0), a single root r = 0 of multiplicity 2.
///
/// disc = 8 so 2 is bad. The Dedekind criterion has F = (f − x·x)/2 = −1, T = gcd(1, x, x) = 1,
/// so ℤ[√2] is maximal at 2 and the listed ideal (2, α) = (√2) is a prime of the ring of
/// integers with ramification index 2 (2 is totally ramified).
#[test]
fn kat8_ramified_prime_with_maximal_order() {
    let k = field_sqrt2();
    let report = dedekind_factor(&k, prime(2), AMPLE).unwrap();

    assert!(report.is_bad_prime);
    assert!(!report.index.divides_index());
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(0), 2)]);
    assert_eq!(report.nonlinear_degree, 0);
}

/// f = x³ − 10, p = 3: f̄ = x³ − 1 = (x − 1)³, so one root r = 1 of multiplicity 3, and 3 | disc =
/// −2700 is bad. Since 10 ≡ 1 (mod 9) the index test is positive (T = x + 2): the ideal (3, α − 1)
/// is a prime of ℤ[α] of norm 3, and the multiplicity is the multiplicity in f̄ only. It is not a
/// ramification index, because ℤ[α] is not maximal at 3: (3) = 𝔭₁²𝔭₂ in the ring of integers.
#[test]
fn kat9_degree_one_report_when_the_index_is_divisible() {
    let k = NumberField::new(poly(&[-10, 0, 0, 1]));
    let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();

    assert!(report.is_bad_prime);
    assert!(report.index.divides_index());
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(1), 3)]);
    assert_eq!(report.nonlinear_degree, 0);
}

/// f = x² − 12, p = 2: f̄ = x², root 0 of multiplicity 2, and the index test is positive.
#[test]
fn kat_divisible_index_with_a_double_root() {
    let k = NumberField::new(poly(&[-12, 0, 1]));
    let report = dedekind_factor(&k, prime(2), AMPLE).unwrap();

    assert!(report.is_bad_prime);
    assert!(report.index.divides_index());
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(0), 2)]);
}

/// f = x² − 3, p = 3: f̄ = x², root 0 of multiplicity 2; ℤ[√3] is 3-maximal, so the index test is
/// negative even though 3 divides disc = 12.
#[test]
fn kat_bad_prime_with_a_maximal_order() {
    let k = NumberField::new(poly(&[-3, 0, 1]));
    let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();

    assert!(report.is_bad_prime);
    assert!(!report.index.divides_index());
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(0), 2)]);
}

/// Φ₅ = x⁴ + x³ + x² + x + 1, p = 5: f̄ = (x − 1)⁴, root 1 of multiplicity 4. ℤ[ζ₅] is maximal,
/// so this is a totally ramified prime of the ring of integers.
#[test]
fn kat_totally_ramified_cyclotomic_prime() {
    let k = NumberField::new(poly(&[1, 1, 1, 1, 1]));
    let report = dedekind_factor(&k, prime(5), AMPLE).unwrap();

    assert!(!report.index.divides_index());
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(1), 4)]);
    assert_eq!(report.nonlinear_degree, 0);
}

/// f = x³ − x − 1, p = 2: f̄ = x³ + x + 1 has f̄(0) = f̄(1) = 1, so no roots and, being a cubic,
/// it is irreducible: the whole degree is nonlinear.
#[test]
fn kat_irreducible_cubic_reduction_is_all_nonlinear() {
    let k = field_cubic();
    let report = dedekind_factor(&k, prime(2), AMPLE).unwrap();

    assert!(report.degree_one.is_empty());
    assert_eq!(report.nonlinear_degree, 3);
}

/// The defining polynomial need not be irreducible for the data about f̄ to be checked: a
/// reducible monic f gives a hand-derivable mix of multiplicities. f = x²(x − 1)(x² + 1), p = 3,
/// has f̄ = x²·(x − 1)·(x² + 1) with x² + 1 irreducible (values 1, 2, 2 at 0, 1, 2): root 0 with
/// multiplicity 2 (found in the second squarefree part) listed before root 1 with multiplicity 1
/// (found in the first), and nonlinear degree 2.
#[test]
fn kat_roots_ascend_across_multiplicity_classes() {
    let f = poly(&[0, 0, 1]).mul(&poly(&[-1, 1])).mul(&poly(&[1, 0, 1]));
    let k = NumberField::new(f);
    let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();

    assert_eq!(roots_and_multiplicities(&report), vec![(bi(0), 2), (bi(1), 1)]);
    assert_eq!(report.nonlinear_degree, 2);
}

/// A multiplicity divisible by p comes from the p-th-power branch of the decomposition.
/// f = (x³ − 1)(x + 1) = x⁴ + x³ − x − 1 (reducible, used only for its reduction), p = 3:
/// f̄ = (x − 1)³(x + 1), so root 1 has multiplicity 3 and root 2 (≡ −1) multiplicity 1.
#[test]
fn kat_multiplicity_equal_to_the_characteristic() {
    let f = poly(&[-1, 0, 0, 1]).mul(&poly(&[1, 1]));
    let k = NumberField::new(f);
    let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();

    assert_eq!(roots_and_multiplicities(&report), vec![(bi(1), 3), (bi(2), 1)]);
    assert_eq!(report.nonlinear_degree, 0);
}

// ─── Budget exhaustion ────────────────────────────────────────────────────────

/// f = x² − 2, p = 7 has roots 3 and 4. The scan finds 3 on its 4th evaluation and 4 on its 5th,
/// so a budget of 4 is exhausted with both roots counted but not all found, and 5 succeeds.
/// Exhaustion is an error carrying the exact root count, never a (wrong) report.
#[test]
fn kat_truncated_budget_is_an_error_not_a_claim() {
    let k = field_sqrt2();
    for max_evaluations in [0, 1, 2, 3, 4] {
        let limits = RootSearchLimits { max_evaluations };
        let outcome = dedekind_factor(&k, prime(7), limits);
        assert!(
            matches!(
                outcome,
                Err(DedekindError::RootSearchExhausted {
                    distinct_root_count: 2,
                    evaluations,
                }) if evaluations == max_evaluations
            ),
            "budget {max_evaluations} must be reported as exhausted"
        );
    }
    let report = dedekind_factor(&k, prime(7), RootSearchLimits { max_evaluations: 5 }).unwrap();
    assert_eq!(roots_and_multiplicities(&report), vec![(bi(3), 1), (bi(4), 1)]);
}

/// The root count is computed exactly before any scan, so a rootless reduction needs no budget:
/// an inert prime is reported as inert even under a zero budget.
#[test]
fn kat_rootless_reduction_needs_no_budget() {
    let k = field_sqrt2();
    let report = dedekind_factor(&k, prime(3), RootSearchLimits { max_evaluations: 0 }).unwrap();
    assert!(report.degree_one.is_empty());
    assert_eq!(report.nonlinear_degree, 2);
}

// ─── Degree invariant ─────────────────────────────────────────────────────────

/// For every fixture the multiplicities of the degree-one primes plus the nonlinear degree sum to
/// deg f, the roots ascend, and every ideal is (p, α − r) with r the root's residue in [0, p).
#[test]
fn kat_multiplicities_plus_nonlinear_degree_equal_the_degree() {
    let fixtures: [(&[i64], u64); 13] = [
        (&[-2, 0, 1], 2),
        (&[-2, 0, 1], 3),
        (&[-2, 0, 1], 7),
        (&[-12, 0, 1], 2),
        (&[-3, 0, 1], 3),
        (&[-1, -1, 0, 1], 2),
        (&[-1, -1, 0, 1], 5),
        (&[-10, 0, 0, 1], 3),
        (&[-2, 0, 0, 1], 3),
        (&[-2, 0, 0, 1], 2),
        (&[1, 1, 1, 1, 1], 5),
        (&[1, 1, 1], 2),
        (&[1, 1, 1], 3),
    ];
    for (coeffs, p) in fixtures {
        let k = NumberField::new(poly(coeffs));
        let report = dedekind_factor(&k, prime(p), AMPLE).unwrap();

        let total: usize = report.degree_one.iter().map(|d| d.multiplicity).sum();
        assert_eq!(total + report.nonlinear_degree, k.degree(), "f = {coeffs:?}, p = {p}");

        let rs: Vec<&BigInt> = report.degree_one.iter().map(|d| &d.ideal.r).collect();
        assert!(rs.windows(2).all(|w| w[0] < w[1]), "ascending roots, f = {coeffs:?}, p = {p}");
        for d in &report.degree_one {
            assert_eq!(d.ideal.p, bi(p as i64));
            assert!(d.ideal.r >= bi(0) && d.ideal.r < bi(p as i64));
            assert_eq!(k.f.eval(&d.ideal.r) % bi(p as i64), bi(0), "r is a root mod p");
        }
    }
}

// ─── Dedekind index criterion ─────────────────────────────────────────────────

/// f = x² − 12, p = 2: 2 divides the index. α = 2√3, so ℤ[α] = ℤ[2√3] has index 2 in ℤ[√3].
///
/// Derivation: f̄ = x², ḡ = h̄ = x. Lifts g = h = x give g·h = x², F = (f − x²)/2 = −6, F̄ = 0.
/// T = gcd(0, x, x) = x, of positive degree.
#[test]
fn index_divisible_x2_minus_12_at_2() {
    let t = dedekind_index_test(&poly(&[-12, 0, 1]), prime(2)).unwrap();
    assert!(t.divides_index());
    assert_eq!(t.obstruction(), &fp(2, &[0, 1]));
    assert_eq!(t.prime(), prime(2));
}

/// f = x² − 3, p = 3: 3 does not divide the index (ℤ[√3] is maximal, since 3 ≡ 3 mod 4).
///
/// Derivation: f̄ = x², ḡ = h̄ = x, g = h = x, F = (f − x²)/3 = −1, F̄ = 2 (a unit),
/// so T = gcd(2, x, x) = 1. Reducing f mod 3 before forming F would give F = 0 and a false T = x.
#[test]
fn index_not_divisible_x2_minus_3_at_3() {
    let t = dedekind_index_test(&poly(&[-3, 0, 1]), prime(3)).unwrap();
    assert!(!t.divides_index());
    assert_eq!(t.obstruction(), &FpPoly::one(prime(3)));
}

/// f = x² − 2, p = 2: a bad prime, but the index is not divisible (F = (f − x²)/2 = −1, T = 1).
#[test]
fn index_not_divisible_x2_minus_2_at_2() {
    let t = dedekind_index_test(&poly(&[-2, 0, 1]), prime(2)).unwrap();
    assert!(!t.divides_index());
    assert_eq!(t.obstruction(), &FpPoly::one(prime(2)));
}

/// f = x² + 3, p = 2: 2 divides the index; ℤ[√−3] ⊂ ℤ[(1+√−3)/2] with index 2.
///
/// Derivation: f̄ = x² + 1 = (x + 1)², ḡ = h̄ = x + 1. Lifts g = h = x + 1 give
/// g·h = x² + 2x + 1, F = (f − g·h)/2 = (2 − 2x)/2 = 1 − x, F̄ = x + 1. T = x + 1.
#[test]
fn index_divisible_x2_plus_3_at_2() {
    let t = dedekind_index_test(&poly(&[3, 0, 1]), prime(2)).unwrap();
    assert!(t.divides_index());
    assert_eq!(t.obstruction(), &fp(2, &[1, 1]));
}

/// f = x³ − 10, p = 3: 3 divides the index because 10 ≡ 1 (mod 9) (Dedekind's pure-cubic case).
///
/// Derivation: f̄ = x³ − 1 = (x − 1)³ = (x + 2)³, ḡ = x + 2, h̄ = (x + 2)². Lifts g = x + 2 and
/// h = x² + 4x + 4 give g·h = x³ + 6x² + 12x + 8, F = (f − g·h)/3 = −2x² − 4x − 6, so
/// F̄ = x² + 2x = x(x + 2). T = gcd(x(x + 2), x + 2, (x + 2)²) = x + 2.
/// Cross-check: disc(f) = −2700 and disc(ℚ(∛10)) = −300, so the index is √(2700/300) = 3.
#[test]
fn index_divisible_x3_minus_10_at_3() {
    let t = dedekind_index_test(&poly(&[-10, 0, 0, 1]), prime(3)).unwrap();
    assert!(t.divides_index());
    assert_eq!(t.obstruction(), &fp(3, &[2, 1]));
}

/// f = x³ − 2, p = 3: ℤ[∛2] is maximal (2 ≢ ±1 mod 9), so 3 does not divide the index.
///
/// Derivation: f̄ = x³ + 1 = (x + 1)³, ḡ = x + 1, h̄ = (x + 1)². Lifts g = x + 1,
/// h = x² + 2x + 1 give g·h = x³ + 3x² + 3x + 1, F = −x² − x − 1, F̄ = 2(x − 1)². Since
/// (x − 1)² at x = −1 is 4 ≡ 1 ≠ 0, F̄ is coprime to x + 1 and T = 1.
#[test]
fn index_not_divisible_x3_minus_2_at_3() {
    let t = dedekind_index_test(&poly(&[-2, 0, 0, 1]), prime(3)).unwrap();
    assert!(!t.divides_index());
    assert_eq!(t.obstruction(), &FpPoly::one(prime(3)));
}

/// f = x³ − 2, p = 2: Eisenstein at 2, so ℤ[∛2] is 2-maximal.
///
/// Derivation: f̄ = x³, ḡ = x, h̄ = x²; g = x, h = x², F = (f − x³)/2 = −1, T = 1.
#[test]
fn index_not_divisible_x3_minus_2_at_2() {
    let t = dedekind_index_test(&poly(&[-2, 0, 0, 1]), prime(2)).unwrap();
    assert!(!t.divides_index());
    assert_eq!(t.obstruction(), &FpPoly::one(prime(2)));
}

/// The fifth cyclotomic polynomial at p = 5: f̄ = (x − 1)⁴ is totally ramified, yet ℤ[ζ₅] is
/// maximal.
///
/// Derivation: ḡ = x − 1, h̄ = (x − 1)³. With lifts g = x − 1, h = (x − 1)³ the difference
/// f − g·h = 5x³ − 5x² + 5x, so F = x³ − x² + x = x(x² − x + 1). At x = 1 this is 1 ≠ 0, so
/// F̄ is coprime to x − 1 and T = 1.
#[test]
fn index_not_divisible_cyclotomic_5_at_5() {
    let t = dedekind_index_test(&poly(&[1, 1, 1, 1, 1]), prime(5)).unwrap();
    assert!(!t.divides_index());
}

/// A prime not dividing disc(f) has f̄ squarefree, so h̄ = 1 and T = 1.
#[test]
fn index_not_divisible_at_good_primes() {
    for (f, p) in [(poly(&[-2, 0, 1]), 7), (poly(&[-1, -1, 0, 1]), 5), (poly(&[-1, -1, 0, 1]), 2)] {
        assert!(!is_bad_prime(&f, &bi(p)));
        let t = dedekind_index_test(&f, prime(p as u64)).unwrap();
        assert!(!t.divides_index());
        assert_eq!(t.obstruction(), &FpPoly::one(prime(p as u64)));
    }
}

/// The report carries the same index test as `dedekind_index_test` returns for the same inputs.
#[test]
fn report_index_matches_the_standalone_index_test() {
    for (coeffs, p) in [(&[-12, 0, 1][..], 2), (&[-3, 0, 1], 3), (&[-10, 0, 0, 1], 3)] {
        let f = poly(coeffs);
        let k = NumberField::new(f.clone());
        let report = dedekind_factor(&k, prime(p), AMPLE).unwrap();
        assert_eq!(report.index, dedekind_index_test(&f, prime(p)).unwrap());
    }
}

/// A non-monic polynomial is rejected, as is a zero or constant one.
#[test]
fn index_test_rejects_invalid_polynomials() {
    let p = prime(3);
    assert_eq!(dedekind_index_test(&poly(&[1, 0, 2]), p), Err(DedekindError::NotMonic));
    assert_eq!(dedekind_index_test(&poly(&[1, 0, -1]), p), Err(DedekindError::NotMonic));
    assert_eq!(dedekind_index_test(&poly(&[1]), p), Err(DedekindError::ConstantPolynomial));
    assert_eq!(dedekind_index_test(&poly(&[7]), p), Err(DedekindError::ConstantPolynomial));
    assert_eq!(dedekind_index_test(&IntPoly::zero(), p), Err(DedekindError::ConstantPolynomial));
}
