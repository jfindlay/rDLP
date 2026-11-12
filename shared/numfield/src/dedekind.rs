//! Dedekind factorisation of the ideal (p) in ℤ[α], with the Dedekind index criterion.
//!
//! Given a number field K = ℚ(α) with defining polynomial f ∈ ℤ[x] and a rational prime p,
//! this module decomposes the ideal (p) into prime ideals of ℤ[α] using Dedekind's theorem
//! and decides, by the Dedekind criterion, whether ℤ[α] is the full ring of integers at p.
//!
//! # Dedekind's theorem (concrete form)
//!
//! For a monic f ∈ ℤ[x] and a prime p with p ∤ disc(f):
//! - Factor f mod p. Each irreducible factor of f mod p corresponds to a prime ideal above p.
//! - For a linear factor (x − r), the prime ideal is (p, α − r) in two-element form.
//! - For an irreducible factor of degree > 1, the prime ideal still exists but has no
//!   representation of the form (p, α − r) with r an integer.
//!
//! # Bad primes and primes dividing the index
//!
//! Two distinct notions are in play. A prime p is **bad** for f when p | disc(f). A prime p
//! **divides the index** when p | [ℤ_K : ℤ[α]]. Since disc(f) = [ℤ_K : ℤ[α]]² · disc(K), every
//! prime dividing the index is bad, but not conversely: p = 2 is bad for x² − 2 (disc 8), yet
//! ℤ[√2] is maximal, so 2 does not divide the index.
//!
//! Dedekind's theorem reads the factorisation of (p) off f mod p exactly when p does not divide
//! the index; [`dedekind_index_test`] decides that.
//!
//! # The Dedekind criterion
//!
//! For monic f ∈ ℤ[x] and a prime p (Cohen, *A Course in Computational Algebraic Number Theory*,
//! §6.1.2, Dedekind's criterion):
//! 1. Factor f̄ = f mod p as ∏ tᵢ^eᵢ over 𝔽_p.
//! 2. Let ḡ = ∏ tᵢ (the radical of f̄) and h̄ = f̄ / ḡ.
//! 3. Lift ḡ and h̄ to monic integer polynomials g and h (any lifts).
//! 4. Form F = (f − g·h) / p from the original integer polynomial f. The numerator is divisible
//!    by p because g·h ≡ f (mod p).
//! 5. Let T = gcd(F̄, ḡ, h̄), monic in 𝔽_p[x].
//! 6. p divides [ℤ_K : ℤ[α]] if and only if deg T > 0.
//!
//! F must be formed from the integer polynomial f itself: the quotient by p carries information
//! that reducing f mod p first would discard. The result does not depend on the lifts: replacing
//! g by g + p·u changes F̄ by a multiple of h̄, and replacing h by h + p·v changes it by a
//! multiple of ḡ, so the gcd is unchanged.
//!
//! # Scope of this implementation
//!
//! [`dedekind_factor`] reports the **degree-one primes** above p: the ideals (p, α − r) for the
//! roots r of f̄ in 𝔽_p, each with its multiplicity as a root of f̄. Irreducible factors of
//! higher degree are not represented by an ideal; their combined contribution to deg f is
//! reported only as [`DedekindReport::nonlinear_degree`]. A prime with no roots mod p therefore
//! yields an empty list and `nonlinear_degree = deg f`: no ideal stands in for it.
//!
//! When p does not divide the index, the multiplicities are ramification indices and the listed
//! ideals are primes of ℤ_K. When p divides the index, ℤ[α] is not maximal at p, so the listed
//! ideals are primes of ℤ[α] of norm p and nothing is claimed about the factorisation of (p) in
//! ℤ_K; that requires the Round 2 / HNF-basis algorithm.
//!
//! The roots come from a budgeted search; running out of budget is an error, never a claim that
//! a root is absent.

use std::fmt;

use num_bigint::BigInt;
use num_traits::Zero;

use crate::element::NumberField;
use crate::fp_poly::{FpPoly, FpPolyError, PrimeModulus, RootEnumeration, RootSearchLimits};
use crate::ideal::Ideal;
use crate::poly::IntPoly;
use crate::resultant::resultant;

// ─── Public types ─────────────────────────────────────────────────────────────

/// The degree-one primes above p, the part of the factorisation of (p) that this crate can
/// represent, with the bookkeeping needed to know what is missing.
///
/// Invariant: the sum of the multiplicities of `degree_one` plus `nonlinear_degree` equals
/// deg f.
///
/// If p does not divide the index (`!index.divides_index()`), each listed ideal is a prime of ℤ_K of residue degree 1
/// whose ramification index is its multiplicity, and `nonlinear_degree` is Σ eᵢfᵢ over the
/// remaining primes (those with fᵢ > 1). If it does divide the index, each listed ideal is a prime
/// of ℤ[α] of norm p and nothing is claimed about the factorisation of (p) in ℤ_K.
pub struct DedekindReport<'a> {
    /// Degree-one primes (p, α − r), in ascending r.
    pub degree_one: Vec<DegreeOnePrime<'a>>,

    /// deg f minus the total multiplicity of the degree-one primes.
    ///
    /// Zero iff f̄ splits into linear factors over 𝔽_p. An inert prime has no degree-one primes
    /// and `nonlinear_degree = deg f`.
    pub nonlinear_degree: usize,

    /// True iff p | disc(f), i.e. p is a bad prime for f.
    ///
    /// Every prime dividing the index is bad, but a bad prime need not divide the index; `index`
    /// gives the finer test.
    pub is_bad_prime: bool,

    /// The Dedekind criterion for p, deciding which reading of `degree_one` applies.
    pub index: IndexTest,
}

/// A degree-one prime (p, α − r) together with the multiplicity of r as a root of f̄.
pub struct DegreeOnePrime<'a> {
    /// The ideal (p, α − r), with r the least nonnegative residue.
    pub ideal: Ideal<'a>,

    /// The multiplicity of r as a root of f mod p.
    pub multiplicity: usize,
}

/// Failure of [`dedekind_index_test`] or [`dedekind_factor`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DedekindError {
    /// The polynomial is not monic.
    NotMonic,
    /// The polynomial is zero or has degree 0; no number field is defined.
    ConstantPolynomial,
    /// A prime-field polynomial operation failed.
    Field(FpPolyError),
    /// A budgeted root enumeration ran out of resources before finishing. This reports
    /// exhaustion of the budget, not a mathematical statement about the polynomial.
    RootSearchExhausted {
        /// Distinct roots found before the budget ran out.
        distinct_root_count: usize,
        /// Evaluations performed before the budget ran out.
        evaluations: u64,
    },
}

impl fmt::Display for DedekindError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotMonic => f.write_str("defining polynomial is not monic"),
            Self::ConstantPolynomial => {
                f.write_str("defining polynomial is zero or constant")
            }
            Self::Field(e) => write!(f, "prime-field polynomial error: {e}"),
            Self::RootSearchExhausted { distinct_root_count, evaluations } => write!(
                f,
                "root search exhausted its budget after {evaluations} evaluations \
                 ({distinct_root_count} distinct roots found)"
            ),
        }
    }
}

impl std::error::Error for DedekindError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Field(e) => Some(e),
            _ => None,
        }
    }
}

impl From<FpPolyError> for DedekindError {
    fn from(e: FpPolyError) -> Self {
        Self::Field(e)
    }
}

/// Outcome of the Dedekind criterion for one polynomial and prime.
///
/// Records T = gcd(F̄, ḡ, h̄) in 𝔽_p[x]; see the module documentation for the construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexTest {
    prime: PrimeModulus,
    obstruction: FpPoly,
}

impl IndexTest {
    /// The prime p that was tested.
    #[must_use]
    pub fn prime(&self) -> PrimeModulus {
        self.prime
    }

    /// The obstruction T = gcd(F̄, ḡ, h̄), monic in 𝔽_p[x].
    #[must_use]
    pub fn obstruction(&self) -> &FpPoly {
        &self.obstruction
    }

    /// True iff p divides [ℤ_K : ℤ[α]], i.e. iff deg T > 0.
    #[must_use]
    pub fn divides_index(&self) -> bool {
        self.obstruction.degree().is_some_and(|d| d > 0)
    }
}

// ─── Public API ───────────────────────────────────────────────────────────────

/// Compute (−1)^(d(d−1)/2) · Res(f, f') for f ∈ ℤ[x] of degree d, where f' is the formal
/// derivative.
///
/// For monic f this is disc(f). For non-monic f no division by lc(f) is performed, so the value
/// is lc(f)·disc(f), where disc(f) = (−1)^(d(d−1)/2) · Res(f, f') / lc(f) is the usual
/// discriminant. It is a multiple of disc(f), and the primes dividing it are those dividing lc(f)
/// or disc(f).
///
/// Returns 0 for polynomials of degree < 1 (constant or zero).
pub fn discriminant(f: &IntPoly) -> BigInt {
    let d = match f.degree() {
        None | Some(0) => return BigInt::zero(),
        Some(d) => d,
    };

    let f_prime = formal_derivative(f);
    let res = resultant(f, &f_prime);

    // Sign factor: (−1)^(d(d−1)/2)
    // d(d−1)/2 is odd iff d ≡ 2 or 3 (mod 4).
    let sign_exp = d * (d - 1) / 2;
    if sign_exp % 2 == 0 { res } else { -res }
}

/// Test whether p divides [`discriminant`]`(f)`, i.e., whether p is a bad prime for f.
///
/// Returns true iff p | [`discriminant`]`(f)`: for monic f, iff p | disc(f); for non-monic f,
/// iff p | lc(f)·disc(f), so a prime dividing only the leading coefficient is also bad. Every
/// prime for which the result is false is good, and for monic f Dedekind's theorem applies
/// directly there. A bad prime may or may not divide the index [ℤ_K : ℤ[α]] (for monic f; the
/// index is defined only for integral α); use [`dedekind_index_test`] (or the `index` of a
/// /// [`DedekindReport`]) to decide that. A polynomial of degree < 1 or with a repeated factor has
/// [`discriminant`] 0, and every prime is then bad.
pub fn is_bad_prime(f: &IntPoly, p: &BigInt) -> bool {
    let disc = discriminant(f);
    if disc.is_zero() {
        // disc = 0 means f is not squarefree; every prime divides 0.
        return true;
    }
    (&disc % p).is_zero()
}

/// Decide whether p divides the index [ℤ_K : ℤ[α]] by the Dedekind criterion.
///
/// `f` must be monic of positive degree; irreducibility is the caller's responsibility, as for
/// [`NumberField`]. The result does not depend on the integer lifts used internally.
///
/// # Errors
///
/// - [`DedekindError::NotMonic`] if `f` is not monic.
/// - [`DedekindError::ConstantPolynomial`] if `f` is zero or constant.
/// - [`DedekindError::Field`] if a prime-field operation fails.
pub fn dedekind_index_test(f: &IntPoly, p: PrimeModulus) -> Result<IndexTest, DedekindError> {
    if f.degree().is_none_or(|d| d == 0) {
        return Err(DedekindError::ConstantPolynomial);
    }
    if !f.is_monic() {
        return Err(DedekindError::NotMonic);
    }

    // f is monic, so f̄ is monic of the same degree; ḡ = rad(f̄) and h̄ = f̄ / ḡ are monic too.
    let f_bar = FpPoly::from_int_poly(f, p);
    let g_bar = f_bar.radical()?;
    let (h_bar, _zero_remainder) = f_bar.div_rem(&g_bar)?;

    let obstruction =
        obstruction_from_lifts(f, p, &g_bar.to_int_poly(), &h_bar.to_int_poly())?;
    Ok(IndexTest { prime: p, obstruction })
}

/// Report the degree-one primes above p, with their multiplicities, and what they leave out.
///
/// The roots of f̄ = f mod p are found on its radical under `limits`; each root's multiplicity
/// comes from the square-free decomposition of f̄. See [`DedekindReport`] for what the listed
/// ideals are (primes of ℤ_K or only of ℤ[α]) and for the invariant tying the report to deg f.
///
/// # Errors
///
/// - [`DedekindError::RootSearchExhausted`] if `limits` ran out before every root was found. No
///   report is returned, since the missing roots would make the multiplicities and
///   `nonlinear_degree` wrong.
/// - [`DedekindError::Field`] if a prime-field operation fails.
pub fn dedekind_factor<'a>(
    field: &'a NumberField,
    p: PrimeModulus,
    limits: RootSearchLimits,
) -> Result<DedekindReport<'a>, DedekindError> {
    let f = &field.f;
    let index = dedekind_index_test(f, p)?;
    let is_bad = is_bad_prime(f, &BigInt::from(p.get()));

    let parts = FpPoly::from_int_poly(f, p).squarefree_decomposition()?;
    let radical = parts.iter().try_fold(FpPoly::one(p), |acc, (s, _)| acc.mul(s))?;
    let roots = match radical.roots(limits)? {
        RootEnumeration::Complete { roots, .. } => roots,
        RootEnumeration::Truncated { distinct_root_count, evaluations, .. } => {
            return Err(DedekindError::RootSearchExhausted { distinct_root_count, evaluations });
        }
    };

    // The parts are pairwise coprime and multiply to the radical, so each root lies in exactly
    // one part and that part's exponent is the root's multiplicity in f̄.
    let mut rooted: Vec<(u64, usize)> = parts
        .iter()
        .flat_map(|(s, e)| roots.iter().filter(|&&r| s.eval(r) == 0).map(move |&r| (r, *e)))
        .collect();
    rooted.sort_unstable();

    let big_p = BigInt::from(p.get());
    let degree_one: Vec<DegreeOnePrime<'a>> = rooted
        .into_iter()
        .map(|(r, multiplicity)| DegreeOnePrime {
            ideal: Ideal::new(field, big_p.clone(), BigInt::from(r)),
            multiplicity,
        })
        .collect();
    let linear_degree: usize = degree_one.iter().map(|d| d.multiplicity).sum();

    Ok(DedekindReport {
        degree_one,
        nonlinear_degree: field.degree() - linear_degree,
        is_bad_prime: is_bad,
        index,
    })
}

// ─── Dedekind criterion ───────────────────────────────────────────────────────

/// Compute T = gcd(F̄, ḡ, h̄) for explicit monic integer lifts `g` and `h` of ḡ and h̄.
///
/// Requires ḡ·h̄ = f̄ so that p divides every coefficient of f − g·h.
fn obstruction_from_lifts(
    f: &IntPoly,
    p: PrimeModulus,
    g: &IntPoly,
    h: &IntPoly,
) -> Result<FpPoly, DedekindError> {
    let big_p = BigInt::from(p.get());
    let numerator = f.sub(&g.mul(h));
    let quotient: Vec<BigInt> = numerator
        .coeffs
        .iter()
        .map(|c| {
            debug_assert!((c % &big_p).is_zero(), "g·h must reduce to f mod p");
            c / &big_p
        })
        .collect();

    let f_cap = FpPoly::from_int_poly(&IntPoly::from_coeffs(quotient), p);
    let g_bar = FpPoly::from_int_poly(g, p);
    let h_bar = FpPoly::from_int_poly(h, p);
    Ok(f_cap.gcd(&g_bar.gcd(&h_bar)?)?)
}

// ─── Polynomial helpers ────────────────────────────────────────────────────────

/// Compute the formal derivative of f: if f = ∑ aᵢ xⁱ, then f' = ∑ i·aᵢ xⁱ⁻¹.
fn formal_derivative(f: &IntPoly) -> IntPoly {
    if f.coeffs.len() <= 1 {
        return IntPoly::zero();
    }
    let coeffs: Vec<BigInt> = f.coeffs[1..]
        .iter()
        .enumerate()
        .map(|(i, c)| BigInt::from(i as i64 + 1) * c)
        .collect();
    IntPoly::from_coeffs(coeffs)
}

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::IntPoly;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    fn poly(coeffs: &[i64]) -> IntPoly {
        IntPoly::from_coeffs(coeffs.iter().map(|&c| bi(c)).collect())
    }

    fn prime(p: u64) -> PrimeModulus {
        PrimeModulus::new(p).unwrap()
    }

    /// f = x² − 2 (defines ℚ(√2))
    fn field_sqrt2() -> NumberField {
        NumberField::new(IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]))
    }

    const AMPLE: RootSearchLimits = RootSearchLimits { max_evaluations: 1_000 };

    #[test]
    fn inert_prime_has_no_degree_one_primes() {
        // p = 3, f = x² − 2: f(0) = −2 ≡ 1, f(1) = −1 ≡ 2, f(2) = 2 (mod 3), so no roots.
        let k = field_sqrt2();
        let report = dedekind_factor(&k, prime(3), AMPLE).unwrap();
        assert!(report.degree_one.is_empty());
        assert_eq!(report.nonlinear_degree, 2);
    }

    #[test]
    fn split_prime_lists_roots_ascending() {
        // p = 7, f = x² − 2: 3² = 9 ≡ 2 and 4² = 16 ≡ 2 mod 7, so roots are 3 and 4.
        let k = field_sqrt2();
        let report = dedekind_factor(&k, prime(7), AMPLE).unwrap();
        let rs: Vec<(BigInt, usize)> =
            report.degree_one.iter().map(|d| (d.ideal.r.clone(), d.multiplicity)).collect();
        assert_eq!(rs, vec![(bi(3), 1), (bi(4), 1)]);
        assert_eq!(report.nonlinear_degree, 0);
    }

    #[test]
    fn exhausted_budget_is_an_error() {
        // Two roots (3 and 4) exist mod 7; a budget of 2 never reaches them.
        let k = field_sqrt2();
        let limits = RootSearchLimits { max_evaluations: 2 };
        let err = dedekind_factor(&k, prime(7), limits).err().unwrap();
        assert_eq!(
            err,
            DedekindError::RootSearchExhausted { distinct_root_count: 2, evaluations: 2 }
        );
    }

    #[test]
    fn formal_derivative_quadratic() {
        // f = x² − 2 → f' = 2x
        let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);
        let fp = formal_derivative(&f);
        assert_eq!(fp.coeffs, vec![bi(0), bi(2)]);
    }

    #[test]
    fn formal_derivative_cubic() {
        // f = x³ − x − 1 → f' = 3x² − 1
        let f = IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)]);
        let fp = formal_derivative(&f);
        assert_eq!(fp.coeffs, vec![bi(-1), bi(0), bi(3)]);
    }

    #[test]
    fn discriminant_quadratic() {
        // disc(x² − 2) = 8
        let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);
        assert_eq!(discriminant(&f), bi(8));
    }

    #[test]
    fn discriminant_cubic() {
        // disc(x³ − x − 1) = −23
        let f = IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)]);
        assert_eq!(discriminant(&f), bi(-23));
    }

    #[test]
    fn is_bad_prime_sqrt2() {
        // disc(x² − 2) = 8; 2 | 8 (bad), 3 ∤ 8 (good)
        let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);
        assert!(is_bad_prime(&f, &bi(2)), "2 should be a bad prime for x²−2");
        assert!(!is_bad_prime(&f, &bi(3)), "3 should not be a bad prime for x²−2");
    }

    #[test]
    fn index_test_good_prime_has_trivial_obstruction() {
        // p = 7 ∤ disc(x² − 2) = 8: f̄ is squarefree, so h̄ = 1 and T = 1.
        let f = poly(&[-2, 0, 1]);
        let t = dedekind_index_test(&f, prime(7)).unwrap();
        assert!(!t.divides_index());
        assert_eq!(t.obstruction(), &FpPoly::one(prime(7)));
        assert_eq!(t.prime(), prime(7));
    }

    #[test]
    fn index_test_bad_prime_sqrt2_is_maximal() {
        // p = 2, f = x² − 2: f̄ = x², ḡ = h̄ = x, lifts g = h = x, F = (f − x²)/2 = −1.
        // T = gcd(1, x, x) = 1, so ℤ[√2] is 2-maximal even though 2 | disc.
        let t = dedekind_index_test(&poly(&[-2, 0, 1]), prime(2)).unwrap();
        assert!(!t.divides_index());
    }

    #[test]
    fn obstruction_is_independent_of_lifts() {
        // f = x³ − 10 at p = 3: f̄ = (x − 1)³, ḡ = x − 1, h̄ = (x − 1)².
        // Shift each lift's coefficients by multiples of p; T must stay x + 2 = x − 1.
        let f = poly(&[-10, 0, 0, 1]);
        let p = prime(3);
        let expected = FpPoly::from_coeffs(p, &[2, 1]);

        // (g, h) as coefficient lists, least significant first; all lift (x−1, (x−1)²) mod 3.
        let lifts: [(&[i64], &[i64]); 4] = [
            (&[-1, 1], &[1, -2, 1]),
            (&[2, 1], &[4, 4, 1]),
            (&[-4, 1], &[1, 1, 1]),
            (&[5, 1], &[-2, 7, 1]),
        ];
        for (g, h) in lifts {
            let t = obstruction_from_lifts(&f, p, &poly(g), &poly(h)).unwrap();
            assert_eq!(t, expected, "lifts g = {g:?}, h = {h:?}");
        }
    }

    #[test]
    fn report_carries_the_index_test() {
        // x² − 12 at 2: bad, and 2 divides the index (ℤ[2√3] ⊊ ℤ[√3]).
        let k = NumberField::new(poly(&[-12, 0, 1]));
        let r = dedekind_factor(&k, prime(2), AMPLE).unwrap();
        assert!(r.is_bad_prime);
        assert!(r.index.divides_index());
    }
}
