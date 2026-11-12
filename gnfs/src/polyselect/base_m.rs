//! Base-m polynomial selection for GNFS.
//!
//! The simplest NFS polynomial generator: given `n` and degree `d`, choose
//! `m` as the smallest integer with `m^{d+1} > n` (i.e., `m = floor(n^{1/(d+1)}) + 1`)
//! and write `n` in base `m`:
//!
//! ```text
//! n = a_0 + a_1·m + a_2·m² + ... + a_d·m^d
//! ```
//!
//! This yields the algebraic-side polynomial `f(x) = Σ a_i x^i` with `f(m) = n`, so
//! `f(m) ≡ 0 (mod n)` exactly. The rational-side polynomial is `g(x) = x − m`.
//!
//! # Non-monic leading coefficient
//!
//! The leading digit `a_d` satisfies `1 ≤ a_d < m`, so `f` is generally non-monic.
//! `PolyPair::monic_f()` performs the standard NFS homogenisation when monic form is required.
//! `PolyPair::number_field()` builds the number field on that monic form; it is used by the NFS
//! discrete-logarithm relation code and by tests, not by the factor base or the square-root
//! stage.
//!
//! # Degree heuristic
//!
//! The optimal degree is the one that minimises the total cost of sieving and linear algebra;
//! the algebraic and rational values have different sizes, and the degree sets the trade-off
//! between them. The classical asymptotic optimum is `d ≈ (3 ln N / ln ln N)^{1/3}`, clamped to
//! `[3, 6]` for toy-scale N. See [`optimal_degree`] for the implementation.
//!
//! # Science↔engineering note
//!
//! Base-m is the *starting point* for polynomial selection, not the end. Murphy-E scoring
//! and root sieving improve upon the base-m polynomial. At toy scale (N < 2^100), the
//! improvement from root sieving is small; at cryptographic scale (RSA-768+), it is essential.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use shared_numfield::IntPoly;

use super::{PolyGenerator, PolyPair};

// ─── optimal_degree ──────────────────────────────────────────────────────────

/// Compute the optimal degree `d` for NFS polynomial selection given `n`.
///
/// Uses the heuristic `d ≈ (3 ln N / ln ln N)^{1/3}`, clamped to `[3, 6]` for toy-scale N.
///
/// The heuristic is the asymptotic degree that minimises the total cost of the NFS, trading the
/// size of the algebraic value against that of the rational value. At cryptographic scale
/// (RSA-768+), this returns 5–6; at toy scale (60–100 bit), it returns 3–4.
///
/// The bit length of `n` is used to approximate `ln N ≈ bits · ln 2`.
///
/// :param n: The integer to factor. Must be positive.
/// :returns: The recommended polynomial degree, in `[3, 6]`.
pub fn optimal_degree(n: &BigInt) -> usize {
    // Approximate ln(N) from the bit length: ln(N) ≈ bits * ln(2).
    let bits = n.bits() as f64;
    let ln_n = bits * std::f64::consts::LN_2;

    if ln_n <= 0.0 {
        return 3;
    }

    let ln_ln_n = ln_n.ln();
    if ln_ln_n <= 0.0 {
        return 3;
    }

    // d ≈ (3 ln N / ln ln N)^{1/3}
    let d_float = (3.0 * ln_n / ln_ln_n).cbrt();
    let d = d_float.round() as usize;

    // Clamp to [3, 6] for toy-scale N.
    d.clamp(3, 6)
}

// ─── base-m expansion ────────────────────────────────────────────────────────

/// Compute the base `m` for base-m polynomial selection of degree `d`.
///
/// Returns the smallest integer `m` such that `m^{d+1} > n`, i.e.,
/// `m = floor(n^{1/(d+1)}) + 1`. This guarantees that `n` has at most `d+1` digits
/// in base `m` (indices 0 through `d`), so the base-m expansion `n = Σ a_i m^i` with
/// `0 ≤ a_i < m` terminates exactly at degree `d`.
///
/// Uses Newton's method in `BigInt` arithmetic to compute `floor(n^{1/(d+1)})`.
///
/// :param n: The integer to expand. Must be positive.
/// :param d: The polynomial degree; the root taken is `(d+1)`-th.
/// :returns: The smallest `m` with `m^{d+1} > n`.
fn base_m_for_degree(n: &BigInt, d: usize) -> BigInt {
    let root_degree = d + 1;

    if n.is_zero() {
        return BigInt::from(2u32); // degenerate: any m >= 2 works
    }

    // Initial estimate via f64: use the bit length to avoid overflow.
    let bits = n.bits() as f64;
    let approx_bits = bits / root_degree as f64;
    // Start with 2^(approx_bits+1) as the initial guess (slightly above the root).
    let init_bits = approx_bits.ceil() as u64 + 2;
    let mut x = BigInt::one() << init_bits;

    // Newton's method for the (d+1)-th root:
    // x_{k+1} = ((root_degree - 1) * x_k + n / x_k^{root_degree-1}) / root_degree
    let r = BigInt::from(root_degree as u64);
    let r_minus_1 = BigInt::from((root_degree - 1) as u64);

    loop {
        let x_pow = pow_bigint_local(&x, root_degree - 1);
        let x_next = (&r_minus_1 * &x + n / &x_pow) / &r;

        if x_next >= x {
            // Converged (or overshot); x is our floor estimate.
            break;
        }
        x = x_next;
    }

    // x is now floor(n^{1/root_degree}) or possibly one too large; correct downward.
    // We want the largest x such that x^root_degree <= n.
    loop {
        let x_pow = pow_bigint_local(&x, root_degree);
        if x_pow <= *n {
            break;
        }
        x -= BigInt::one();
    }
    // Also verify x+1 is too large (shouldn't be needed after Newton, but be safe).
    loop {
        let x1 = &x + BigInt::one();
        let x1_pow = pow_bigint_local(&x1, root_degree);
        if x1_pow > *n {
            break;
        }
        x = x1;
    }

    // x = floor(n^{1/(d+1)}). Return x + 1 = the smallest m with m^{d+1} > n.
    x + BigInt::one()
}

/// Compute `base^exp` for `BigInt` using repeated squaring (module-local copy).
fn pow_bigint_local(base: &BigInt, exp: usize) -> BigInt {
    if exp == 0 {
        return BigInt::one();
    }
    let mut result = BigInt::one();
    let mut b = base.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result *= &b;
        }
        b = &b * &b;
        e >>= 1;
    }
    result
}

/// Write `n` in base `m` to obtain the coefficients of the algebraic-side polynomial.
///
/// Returns the coefficient vector `[a_0, a_1, ..., a_d]` (least-significant first) such that
/// `n = a_0 + a_1·m + ... + a_d·m^d`. The vector always has exactly `d + 1` entries; if
/// `n < m^d`, the leading entries are zero.
///
/// :param n: The integer to expand. Must be positive.
/// :param m: The base. Must be ≥ 2.
/// :param d: The polynomial degree.
/// :returns: Coefficient vector of length `d + 1`.
fn base_m_digits(n: &BigInt, m: &BigInt, d: usize) -> Vec<BigInt> {
    let mut coeffs = Vec::with_capacity(d + 1);
    let mut remainder = n.clone();

    for _ in 0..=d {
        let digit = &remainder % m;
        coeffs.push(digit.clone());
        remainder = (remainder - digit) / m;
    }

    // After d+1 digits, remainder should be zero (m^{d+1} > n by construction of m).
    // If not, the caller chose m too small; we still return what we have.
    coeffs
}

// ─── select_base_m ───────────────────────────────────────────────────────────

/// Generate a polynomial pair via base-m expansion.
///
/// Given `n` and degree `d`, takes `m` as the least integer with `m^{d+1} > n`, i.e.
/// `m = floor(n^{1/(d+1)}) + 1`, and writes `n` in base `m`:
///
/// ```text
/// n = a_0 + a_1·m + a_2·m² + ... + a_d·m^d
/// ```
///
/// yielding `f(x) = Σ a_i x^i` with `f(m) = n`. The rational side is `g(x) = x − m`.
///
/// The resulting `f` is generally non-monic (`a_d < m`). This is the simplest polynomial
/// generator; Murphy-E scoring and root sieving improve upon it.
///
/// :param n: The integer to factor. Must be > 1.
/// :param degree: The polynomial degree `d`. Typically 3–6; use [`optimal_degree`] to choose.
/// :returns: A `PolyPair` satisfying `f(m) = n` and `g(m) = 0`.
pub fn select_base_m(n: &BigInt, degree: usize) -> PolyPair {
    let m = base_m_for_degree(n, degree);
    select_base_m_with_m(n, &m, degree)
}

/// Generate a polynomial pair via base-m expansion with a specified `m`.
///
/// Used for reproducibility testing (e.g., matching CADO-NFS published polynomials where
/// the exact `m` is known). The caller is responsible for ensuring `m` is appropriate for
/// the given `n` and `degree`.
///
/// :param n: The integer to factor.
/// :param m: The base for the expansion.
/// :param degree: The polynomial degree `d`.
/// :returns: A `PolyPair` satisfying `f(m) = n` and `g(m) = 0`.
pub fn select_base_m_with_m(n: &BigInt, m: &BigInt, degree: usize) -> PolyPair {
    let coeffs = base_m_digits(n, m, degree);
    let f = IntPoly::from_coeffs(coeffs.into_iter().map(|c| c).collect());
    // g(x) = x − m: coefficients [-m, 1] (least-significant first).
    let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
    PolyPair::new(f, g, m.clone(), n.clone())
}

// ─── BaseMGenerator ──────────────────────────────────────────────────────────

/// Base-m generator: produces a single polynomial pair via base-m expansion.
///
/// Implements [`PolyGenerator`] so that base-m fits into the common score-and-rank pipeline
/// alongside root sieve and Coppersmith generators.
///
/// The generator produces exactly one candidate per `(n, degree)` pair. Use
/// `.take(1)` or collect the single element.
pub struct BaseMGenerator {
    /// The integer to factor.
    pub n: BigInt,
    /// The polynomial degree.
    pub degree: usize,
}

impl PolyGenerator for BaseMGenerator {
    fn generate(&self) -> impl Iterator<Item = PolyPair> {
        std::iter::once(select_base_m(&self.n, self.degree))
    }
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn base_m_for_degree_basic() {
        // For n=8, d=2 (root_degree=3): floor(8^{1/3}) = 2, so m = 3.
        // 3^3 = 27 > 8. ✓
        assert_eq!(base_m_for_degree(&bi(8), 2), bi(3));
        // For n=27, d=2: floor(27^{1/3}) = 3, so m = 4.
        // 4^3 = 64 > 27. ✓
        assert_eq!(base_m_for_degree(&bi(27), 2), bi(4));
        // For n=26, d=2: floor(26^{1/3}) = 2, so m = 3.
        // 3^3 = 27 > 26. ✓
        assert_eq!(base_m_for_degree(&bi(26), 2), bi(3));
    }

    #[test]
    fn base_m_for_degree_guarantees_fit() {
        // For any n and d, m^{d+1} > n must hold.
        for n_val in [15i64, 100, 1000, 1022117, 999983 * 999979] {
            for d in [2usize, 3, 4, 5] {
                let n = BigInt::from(n_val);
                let m = base_m_for_degree(&n, d);
                let m_pow = pow_bigint_local(&m, d + 1);
                assert!(
                    m_pow > n,
                    "m^{{d+1}} should be > n: m={m}, d={d}, n={n}, m^{{d+1}}={m_pow}"
                );
            }
        }
    }

    #[test]
    fn base_m_digits_basic() {
        // 15 in base 3, degree 2: 15 = 0 + 2*3 + 1*9 → [0, 2, 1]
        let digits = base_m_digits(&bi(15), &bi(3), 2);
        assert_eq!(digits, vec![bi(0), bi(2), bi(1)]);
    }

    #[test]
    fn select_base_m_round_trip() {
        // For any n, m, d: f(m) should equal n.
        let n = bi(1009 * 1013);
        let pair = select_base_m(&n, 3);
        assert_eq!(pair.f.eval(&pair.m), n);
        assert_eq!(pair.g.eval(&pair.m), BigInt::zero());
        assert_eq!(pair.verify(), Ok(()));
    }

    #[test]
    fn optimal_degree_toy() {
        // 60-bit N: optimal_degree should return 3 or 4
        let n_60bit = BigInt::from(1u64 << 60);
        let d = optimal_degree(&n_60bit);
        assert!(d == 3 || d == 4, "expected 3 or 4, got {d}");
    }

    #[test]
    fn base_m_generator_produces_one() {
        let n = bi(1009 * 1013);
        let generator = BaseMGenerator { n: n.clone(), degree: 3 };
        let pairs: Vec<_> = generator.generate().collect();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].verify(), Ok(()));
    }
}
