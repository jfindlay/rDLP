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
//! - For an irreducible factor of degree > 1, the prime ideal still exists but its two-element
//!   representation requires a second generator that is not simply α − r for an integer r.
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
//! This implementation handles the **linear-factor case** for both good and bad primes:
//! it finds all r ∈ {0, …, p−1} with f(r) ≡ 0 (mod p) and returns the corresponding prime
//! ideals (p, α − r). When p divides the index, the result is flagged and a caveat is
//! documented: ℤ[α] is not the full ring of integers at p, so those ideals are primes of ℤ[α]
//! only, and the factorisation of (p) in ℤ_K requires the full Round 2 / HNF-basis algorithm.
//!
//! # Inert primes
//!
//! If f has no roots mod p (all factors of f mod p have degree > 1), then (p) is inert (or
//! factors only into higher-degree prime ideals). In this case, `dedekind_factor` returns a
//! single sentinel ideal `(p, α − 0)` to indicate that (p) is prime (or has no linear-factor
//! decomposition). This sentinel is documented by the `is_inert` flag on the return type.
//!
//! **Convention**: an empty `Vec` would be ambiguous (no ideals vs. inert). Instead, the caller
//! can check whether the returned `Vec` has length 1 with `r == 0` and the field has no roots
//! mod p — but for simplicity, the function returns the sentinel `Ideal { p, r: 0 }` for inert
//! primes. Callers that only care about linear-factor ideals should filter by checking
//! `f.eval(&r) % p == 0` for each returned ideal's `r`.
//!
//! # Panics
//!
//! Panics if `p ≤ 0`.

use std::fmt;

use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::element::NumberField;
use crate::fp_poly::{FpPoly, FpPolyError, PrimeModulus};
use crate::ideal::Ideal;
use crate::poly::IntPoly;
use crate::resultant::resultant;

// ─── Public types ─────────────────────────────────────────────────────────────

/// Result of `dedekind_factor_extended`: prime ideals above p with bad-prime metadata.
pub struct DedekindResult<'a> {
    /// Prime ideals above p in two-element form (p, α − r).
    ///
    /// When p does not divide the index [ℤ_K : ℤ[α]] (`index_divisible = false`), this is the
    /// complete factorisation of (p) restricted to linear factors of f mod p.
    ///
    /// When p divides the index (`index_divisible = true`), this list covers only the linear
    /// factors of f mod p, and those ideals are primes of ℤ[α], not necessarily of ℤ_K. The
    /// prime ideals of ℤ_K above p are not described by it: they require the full Round 2 /
    /// HNF-basis algorithm (out of scope for this implementation).
    pub ideals: Vec<Ideal<'a>>,

    /// True if p | disc(f), i.e., p is a bad prime for f.
    ///
    /// Every prime dividing the index is bad, but a bad prime need not divide the index;
    /// `index_divisible` gives the finer test.
    pub is_bad_prime: bool,

    /// True if p | [ℤ_K : ℤ[α]], as decided by the Dedekind criterion.
    ///
    /// When true, ℤ[α] is not the full ring of integers at p and `ideals` does not describe
    /// the factorisation of (p) in ℤ_K.
    ///
    /// Always false when `is_bad_prime` is false (p ∤ disc(f) implies p ∤ [ℤ_K : ℤ[α]]).
    pub index_divisible: bool,
}

/// Failure of [`dedekind_index_test`].
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

/// Compute the discriminant of a monic polynomial f ∈ ℤ[x].
///
/// For monic f of degree d, disc(f) = (−1)^(d(d−1)/2) · Res(f, f') where f' is the
/// formal derivative. Since f is monic, lc(f) = 1 and no division by lc(f) is needed.
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

/// Test whether p divides disc(f), i.e., whether p is a bad prime for f.
///
/// Returns true iff p | disc(f). For primes not dividing disc(f), Dedekind's theorem applies
/// directly. A bad prime may or may not divide the index [ℤ_K : ℤ[α]]; use
/// [`dedekind_index_test`] (or `dedekind_factor_extended`) to decide that.
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

/// Factorise the ideal (p) in ℤ[α] using Dedekind's theorem.
///
/// Given a number field K = ℚ(α) with defining polynomial f ∈ ℤ[x] and a rational prime p
/// not dividing disc(f), returns the prime ideals above p that correspond to linear factors
/// of f mod p.
///
/// Each returned ideal is in two-element form `(p, α − r)` where r ∈ {0, …, p−1} satisfies
/// f(r) ≡ 0 (mod p).
///
/// # Inert primes
///
/// If f has no roots mod p (all factors of f mod p have degree > 1), the function returns a
/// single sentinel ideal `Ideal { p, r: 0 }`. This represents the convention that (p) is
/// inert (prime) in ℤ[α] with respect to linear factors. The sentinel ideal `(p, α − 0) = (p, α)`
/// generates the same ideal as (p) when f is irreducible mod p.
///
/// # Panics
///
/// Panics if `p ≤ 0`.
///
/// # Note
///
/// This function does NOT handle primes dividing disc(f) (bad primes). For those, use
/// `dedekind_factor_extended`, which also reports whether p divides the index.
pub fn dedekind_factor<'a>(field: &'a NumberField, p: &BigInt) -> Vec<Ideal<'a>> {
    assert!(p.is_positive(), "prime p must be positive, got {p}");

    let f = &field.f;
    let p_usize = p_to_usize(p);

    // Find all roots r ∈ {0, …, p−1} with f(r) ≡ 0 (mod p).
    let mut roots: Vec<BigInt> = Vec::new();
    for r in 0..p_usize {
        let r_big = BigInt::from(r);
        let val = f.eval(&r_big);
        // val mod p, reduced to [0, p)
        let rem = mod_reduce(&val, p);
        if rem.is_zero() {
            roots.push(r_big);
        }
    }

    if roots.is_empty() {
        // Inert prime: no linear factors. Return the sentinel ideal (p, α − 0).
        // The ideal (p, α) = (p, α − 0) generates the same ideal as (p) when f is
        // irreducible mod p, representing the inert prime above p.
        vec![Ideal::new(field, p.clone(), BigInt::zero())]
    } else {
        roots.into_iter().map(|r| Ideal::new(field, p.clone(), r)).collect()
    }
}

/// Factorise the ideal (p) in ℤ[α], reporting whether p is bad and whether it divides the index.
///
/// Returns the prime ideals from linear factors of f mod p, as `dedekind_factor` does, together
/// with `is_bad_prime` (p | disc(f)) and `index_divisible` (p | [ℤ_K : ℤ[α]], decided by
/// [`dedekind_index_test`]). When `index_divisible` is true, ℤ[α] is not maximal at p and the
/// returned ideals are primes of ℤ[α] only; the factorisation of (p) in ℤ_K needs the full
/// Round 2 algorithm.
///
/// # Panics
///
/// Panics if `p` is not a prime that fits in a `u64`.
pub fn dedekind_factor_extended<'a>(field: &'a NumberField, p: &BigInt) -> DedekindResult<'a> {
    assert!(p.is_positive(), "prime p must be positive, got {p}");

    let f = &field.f;
    let modulus = p
        .to_u64()
        .and_then(|v| PrimeModulus::new(v).ok())
        .unwrap_or_else(|| panic!("p must be a prime that fits in a u64, got {p}"));

    let bad = is_bad_prime(f, p);
    // `NumberField` guarantees f is monic of positive degree, so the test cannot fail.
    let index = dedekind_index_test(f, modulus)
        .expect("a NumberField's defining polynomial is monic of positive degree");

    // Find linear factors (roots) of f mod p — same logic as dedekind_factor.
    let p_usize = p_to_usize(p);
    let mut roots: Vec<BigInt> = Vec::new();
    for r in 0..p_usize {
        let r_big = BigInt::from(r);
        let val = f.eval(&r_big);
        let rem = mod_reduce(&val, p);
        if rem.is_zero() {
            roots.push(r_big);
        }
    }

    let ideals = if roots.is_empty() {
        vec![Ideal::new(field, p.clone(), BigInt::zero())]
    } else {
        roots.into_iter().map(|r| Ideal::new(field, p.clone(), r)).collect()
    };

    DedekindResult { ideals, is_bad_prime: bad, index_divisible: index.divides_index() }
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

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Reduce `a` into the canonical range [0, m) for m > 0.
fn mod_reduce(a: &BigInt, m: &BigInt) -> BigInt {
    let r = a % m;
    if r.is_negative() { r + m } else { r }
}

/// Convert a `BigInt` prime to `usize` for iteration.
///
/// Panics if p does not fit in usize (primes this large are not supported by the
/// brute-force root-finding loop).
fn p_to_usize(p: &BigInt) -> usize {
    p.to_usize().expect("prime p must fit in usize for root-finding iteration")
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

    #[test]
    fn inert_prime_returns_sentinel() {
        // p = 3, f = x² − 2: 2 is a QNR mod 3, so f is irreducible mod 3.
        // f(0) = -2 ≡ 1, f(1) = -1 ≡ 2, f(2) = 2 ≡ 2 mod 3 — no roots.
        let k = field_sqrt2();
        let ideals = dedekind_factor(&k, &bi(3));
        assert_eq!(ideals.len(), 1, "inert prime should return one sentinel ideal");
        assert_eq!(ideals[0].p, bi(3));
        assert_eq!(ideals[0].r, bi(0), "sentinel ideal has r = 0");
    }

    #[test]
    fn split_prime_returns_two_ideals() {
        // p = 7, f = x² − 2: 3² = 9 ≡ 2 and 4² = 16 ≡ 2 mod 7, so roots are 3 and 4.
        let k = field_sqrt2();
        let ideals = dedekind_factor(&k, &bi(7));
        assert_eq!(ideals.len(), 2, "split prime should return two ideals");
        let mut rs: Vec<i64> = ideals.iter().map(|i| i.r.to_string().parse().unwrap()).collect();
        rs.sort();
        assert_eq!(rs, vec![3i64, 4], "roots of x²−2 mod 7 should be 3 and 4");
    }

    #[test]
    fn panics_on_zero_p() {
        let k = field_sqrt2();
        let result = std::panic::catch_unwind(|| {
            dedekind_factor(&k, &bi(0));
        });
        assert!(result.is_err(), "dedekind_factor should panic when p = 0");
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
    fn extended_routes_through_index_test() {
        // x² − 12 at 2: bad, and 2 divides the index (ℤ[2√3] ⊊ ℤ[√3]).
        let k = NumberField::new(poly(&[-12, 0, 1]));
        let r = dedekind_factor_extended(&k, &bi(2));
        assert!(r.is_bad_prime);
        assert!(r.index_divisible);
    }

    #[test]
    fn extended_panics_on_composite_p() {
        let k = field_sqrt2();
        let result = std::panic::catch_unwind(|| {
            dedekind_factor_extended(&k, &bi(4));
        });
        assert!(result.is_err(), "a composite modulus is not a prime");
    }
}
