//! Known-answer tests for `dedekind_factor`, `dedekind_factor_extended`, and
//! `dedekind_index_test`: Dedekind factorisation of (p) in ℤ[α] and the index criterion.
//!
//! Tests cover inert primes, split primes, cubic fields, norm-product verification,
//! discriminant computation, bad-prime detection, and the Dedekind index criterion with
//! hand-derived obstruction polynomials.
//!
//! # Inert-prime convention
//!
//! When f has no roots mod p (all irreducible factors of f mod p have degree > 1), the ideal (p)
//! is inert (or factors only into higher-degree prime ideals). `dedekind_factor` returns a single
//! sentinel ideal `Ideal { p, r: 0 }` in this case. The sentinel represents the convention that
//! (p, α − 0) = (p, α) generates the same ideal as (p) when f is irreducible mod p.

use num_bigint::BigInt;
use shared_numfield::{
    dedekind_factor, dedekind_factor_extended, dedekind_index_test, discriminant, is_bad_prime,
    DedekindError, FpPoly, IntPoly, NumberField, PrimeModulus,
};

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

// ─── KAT 1 — Inert prime ──────────────────────────────────────────────────────

/// f = x² − 2, p = 3: 2 is a quadratic non-residue mod 3 (Legendre(2,3) = −1),
/// so f is irreducible mod 3 and (3) is inert in ℤ[√2].
///
/// Verification: f(0) = −2 ≡ 1, f(1) = −1 ≡ 2, f(2) = 2 ≡ 2 (mod 3) — no roots.
///
/// `dedekind_factor` returns the sentinel ideal (3, α − 0) indicating (3) is inert.
#[test]
fn kat1_inert_prime() {
    let k = field_sqrt2();
    let ideals = dedekind_factor(&k, &bi(3));

    // Inert: single sentinel ideal
    assert_eq!(
        ideals.len(),
        1,
        "inert prime p=3 should return one sentinel ideal, got {} ideals",
        ideals.len()
    );
    assert_eq!(ideals[0].p, bi(3), "sentinel ideal should have p = 3");
    assert_eq!(ideals[0].r, bi(0), "sentinel ideal for inert prime should have r = 0");
}

// ─── KAT 2 — Split prime ──────────────────────────────────────────────────────

/// f = x² − 2, p = 7: 2 is a quadratic residue mod 7 (3² = 9 ≡ 2, 4² = 16 ≡ 2 mod 7),
/// so f splits as (x − 3)(x − 4) mod 7 and (7) splits into two prime ideals in ℤ[√2].
///
/// `dedekind_factor` returns two ideals: (7, α − 3) and (7, α − 4).
#[test]
fn kat2_split_prime() {
    let k = field_sqrt2();
    let ideals = dedekind_factor(&k, &bi(7));

    assert_eq!(
        ideals.len(),
        2,
        "split prime p=7 should return two ideals, got {} ideals",
        ideals.len()
    );

    // Both ideals should have p = 7
    for ideal in &ideals {
        assert_eq!(ideal.p, bi(7), "each ideal above 7 should have p = 7");
    }

    // The r values should be 3 and 4 (in some order)
    let mut rs: Vec<BigInt> = ideals.iter().map(|i| i.r.clone()).collect();
    rs.sort();
    assert_eq!(rs, vec![bi(3), bi(4)], "roots of x²−2 mod 7 should be 3 and 4");
}

// ─── KAT 3 — Cubic field, partial split ───────────────────────────────────────

/// f = x³ − x − 1, p = 5: f factors as (x − 2)(irreducible quadratic) mod 5.
///
/// Verification: f(2) = 8 − 2 − 1 = 5 ≡ 0 (mod 5). The other factor x² + 2x + 3 has
/// discriminant 4 − 12 = −8 ≡ 2 (mod 5), and Legendre(2, 5) = −1, so it is irreducible.
///
/// `dedekind_factor` returns at least one ideal with r = 2 (the linear factor).
/// The irreducible quadratic factor also contributes a prime ideal above 5, but its
/// two-element representation is not handled by this implementation (linear-factor scope only).
#[test]
fn kat3_cubic_partial_split() {
    let k = field_cubic();
    let ideals = dedekind_factor(&k, &bi(5));

    // At least one ideal with r = 2 (from the linear factor x − 2)
    let has_r2 = ideals.iter().any(|i| i.r == bi(2));
    assert!(
        has_r2,
        "dedekind_factor for f=x³−x−1, p=5 should include an ideal with r=2 (root of f mod 5)"
    );

    // All returned ideals should have p = 5
    for ideal in &ideals {
        assert_eq!(ideal.p, bi(5), "each ideal above 5 should have p = 5");
    }

    // Verify r=2 is indeed a root: f(2) = 8 − 2 − 1 = 5 ≡ 0 (mod 5)
    let f = &k.f;
    let val = f.eval(&bi(2));
    assert_eq!(&val % bi(5), bi(0), "f(2) should be ≡ 0 (mod 5)");
}

// ─── KAT 4 — Norm product ─────────────────────────────────────────────────────

/// For f = x² − 2, p = 7 (totally split), the product of norms of the returned ideals
/// equals p^deg(f) = 7² = 49.
///
/// Each prime ideal (7, α − r) has norm 7 (residue degree 1). The product is 7 · 7 = 49.
/// This verifies the fundamental identity: ∏ N(𝔭ᵢ)^eᵢ = p^[K:ℚ].
#[test]
fn kat4_norm_product() {
    let k = field_sqrt2();
    let ideals = dedekind_factor(&k, &bi(7));

    // For a totally split prime in a degree-2 field, we get exactly 2 ideals.
    assert_eq!(ideals.len(), 2, "p=7 should split into 2 prime ideals in ℚ(√2)");

    // Product of norms = p^d = 7^2 = 49
    let norm_product: BigInt = ideals.iter().map(|i| i.norm()).product();
    let expected = bi(7).pow(2u32);
    assert_eq!(
        norm_product, expected,
        "product of norms should equal p^d = 7^2 = 49, got {norm_product}"
    );
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

// ─── KAT 8 — Bad prime extended: f = x² − 2, p = 2 ──────────────────────────

/// For f = x² − 2, p = 2: `dedekind_factor_extended` returns `is_bad_prime = true`.
///
/// Analysis:
/// - disc(x² − 2) = 8, and 2 | 8, so p = 2 is a bad prime.
/// - f mod 2 = x² (since −2 ≡ 0 mod 2), f' mod 2 = 0 (since 2x ≡ 0 mod 2).
/// - Dedekind criterion: f̄ = x², so ḡ = x and h̄ = x; lifts g = h = x give
///   F = (f − g·h)/2 = (−2)/2 = −1, F̄ = 1, and T = gcd(1, x, x) = 1.
///   So index_divisible = false: ℤ[√2] IS the full ring of integers at p = 2.
/// - Root of f mod 2: f(0) = −2 ≡ 0 (mod 2), f(1) = −1 ≡ 1 (mod 2). Root: r = 0.
/// - One ideal: (2, α − 0) = (2, α), representing the unique prime above 2 in ℤ[√2].
///
/// The prime ideal (2, α) = (√2) satisfies (2, α)² = (2) in ℤ[√2], confirming that 2
/// is totally ramified in ℚ(√2).
#[test]
fn kat8_bad_prime_extended() {
    let k = field_sqrt2();
    let result = dedekind_factor_extended(&k, &bi(2));

    assert!(result.is_bad_prime, "p=2 should be flagged as a bad prime for x²−2");
    // ℤ[√2] is the maximal order, so the Dedekind criterion gives T=1 (index not divisible).
    assert!(
        !result.index_divisible,
        "index should not be divisible by p=2 for x²−2 (ℤ[√2] is the maximal order)"
    );

    // One ideal with r = 0 (the unique prime above 2 in ℤ[√2]).
    assert_eq!(
        result.ideals.len(),
        1,
        "p=2 should give one prime ideal above 2 in ℤ[√2], got {}",
        result.ideals.len()
    );
    assert_eq!(result.ideals[0].p, bi(2), "ideal should have p = 2");
    assert_eq!(result.ideals[0].r, bi(0), "ideal above 2 in ℤ[√2] should have r = 0");
}

// ─── KAT 9 — Good prime extended: f = x² − 2, p = 7 ─────────────────────────

/// For f = x² − 2, p = 7: `dedekind_factor_extended` returns `is_bad_prime = false`,
/// `index_divisible = false`, and two ideals with r ∈ {3, 4}.
///
/// This is the same split as KAT 2, now verified through the extended interface.
/// Since 7 ∤ disc(x²−2) = 8, p = 7 is a good prime and the Dedekind criterion is
/// not applied (index_divisible is always false for good primes).
#[test]
fn kat9_good_prime_extended() {
    let k = field_sqrt2();
    let result = dedekind_factor_extended(&k, &bi(7));

    assert!(!result.is_bad_prime, "p=7 should not be a bad prime for x²−2");
    assert!(
        !result.index_divisible,
        "index should not be divisible by p=7 (good prime)"
    );

    assert_eq!(
        result.ideals.len(),
        2,
        "p=7 should split into 2 prime ideals in ℚ(√2), got {}",
        result.ideals.len()
    );

    // Both ideals should have p = 7
    for ideal in &result.ideals {
        assert_eq!(ideal.p, bi(7), "each ideal above 7 should have p = 7");
    }

    // The r values should be 3 and 4 (in some order) — same as KAT 2.
    let mut rs: Vec<BigInt> = result.ideals.iter().map(|i| i.r.clone()).collect();
    rs.sort();
    assert_eq!(rs, vec![bi(3), bi(4)], "roots of x²−2 mod 7 should be 3 and 4");
}

// ─── Panic test ───────────────────────────────────────────────────────────────

/// `dedekind_factor` must panic when p ≤ 0.
#[test]
fn panics_on_nonpositive_p() {
    let k = field_sqrt2();

    let result_zero = std::panic::catch_unwind(|| {
        dedekind_factor(&k, &bi(0));
    });
    assert!(result_zero.is_err(), "dedekind_factor should panic when p = 0");

    let result_neg = std::panic::catch_unwind(|| {
        dedekind_factor(&k, &bi(-7));
    });
    assert!(result_neg.is_err(), "dedekind_factor should panic when p < 0");
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

/// The extended factorisation reports the criterion's answer rather than recomputing its own.
#[test]
fn extended_reports_index_criterion() {
    let k = NumberField::new(poly(&[-12, 0, 1]));
    let r = dedekind_factor_extended(&k, &bi(2));
    assert!(r.is_bad_prime);
    assert!(r.index_divisible, "2 divides the index of ℤ[2√3]");

    let k = NumberField::new(poly(&[-3, 0, 1]));
    let r = dedekind_factor_extended(&k, &bi(3));
    assert!(r.is_bad_prime);
    assert!(!r.index_divisible, "ℤ[√3] is 3-maximal");
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
