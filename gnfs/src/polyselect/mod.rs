//! Polynomial selection for GNFS.
//!
//! This module provides the NFS polynomial pair type (`PolyPair`), the common generator
//! trait (`PolyGenerator`), and the entry surfaces for polynomial selection.
//!
//! # Background
//!
//! GNFS requires two polynomials `f` (algebraic side, degree d ≥ 2) and `g` (rational side,
//! degree 1) sharing a common root `m` modulo `n`:
//!
//! ```text
//! f(m) ≡ 0 (mod n)
//! g(m) = 0,  g(x) = x − m
//! ```
//!
//! The simplest construction is **base-m expansion**: write `n` in base `m` to obtain the
//! coefficients of `f`, so that `f(m) = n ≡ 0 (mod n)` exactly. See [`base_m`] for details.
//!
//! # Non-monic seam
//!
//! Base-m expansion produces a polynomial `f` whose leading coefficient `a_d` is generally
//! less than `m` (i.e., `f` is non-monic). The number-field substrate (`NumberField::new`)
//! requires a monic defining polynomial. `PolyPair` stores the original non-monic `f` and
//! exposes [`PolyPair::number_field`] / [`PolyPair::monic_f`] to perform the standard NFS
//! homogenisation `f(x) → a_d^{d−1} f(x/a_d)` when monic form is needed. The consumers of
//! `number_field` are the NFS discrete-logarithm relation code (which takes the field as a
//! parameter) and tests; the factor base searches roots of the original `f`, the sieve evaluates
//! `F(a, b)` from the original `f`, and the square-root stage builds its own monic form.
//!
//! # Module layout
//!
//! - [`mod.rs`](self) — `PolyPair`, `PolyPairError`, `PolyGenerator`, re-exports.
//! - [`base_m`] — `select_base_m`, `optimal_degree`, `BaseMGenerator`.
//! - [`root_sieve`] — `root_sieve`, `RootSieveConfig`, `RootSieveGenerator`.
//! - [`coppersmith`] — `coppersmith_polys`, `coppersmith_best`, `CoppersmithConfig`,
//!   `CoppersmithGenerator`.

pub mod base_m;
pub mod coppersmith;
pub mod murphy;
pub mod root_sieve;
pub mod roots;

pub use base_m::{select_base_m, select_base_m_with_m, optimal_degree, BaseMGenerator};
pub use coppersmith::{coppersmith_polys, coppersmith_best, CoppersmithConfig, CoppersmithGenerator};
pub use murphy::score;
pub use root_sieve::{root_sieve, rotate, RootSieveConfig, RootSieveGenerator};

use num_bigint::BigInt;
use num_traits::{One, Zero};
use shared_numfield::{IntPoly, NumberField};

// ─── PolyPairError ───────────────────────────────────────────────────────────

/// Error type for [`PolyPair::verify`].
///
/// Each variant describes a specific invariant violation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolyPairError {
    /// `f(m)` is not divisible by `n`.
    ///
    /// :param f_of_m: The actual value of `f(m)`.
    /// :param n: The modulus.
    RootCheckFailed { f_of_m: BigInt, n: BigInt },

    /// `g` is not the expected `x − m` form.
    ///
    /// :param expected_g: The polynomial `x − m` that `g` should equal.
    /// :param actual_g: The actual `g` stored in the pair.
    RationalPolyMismatch { expected_g: IntPoly, actual_g: IntPoly },

    /// The `degree` field does not match `f.degree()`.
    ///
    /// :param field_degree: The value stored in `PolyPair::degree`.
    /// :param poly_degree: The actual degree of `f` (or `None` if `f` is zero).
    DegreeMismatch { field_degree: usize, poly_degree: Option<usize> },

    /// `f` is zero or constant (degree < 1).
    InvalidAlgebraicPoly,
}

impl std::fmt::Display for PolyPairError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RootCheckFailed { f_of_m, n } => {
                write!(f, "root check failed: f(m) = {f_of_m}, which is not divisible by n = {n}")
            }
            Self::RationalPolyMismatch { expected_g, actual_g } => {
                write!(
                    f,
                    "rational polynomial mismatch: expected g = {:?}, got {:?}",
                    expected_g.coeffs, actual_g.coeffs
                )
            }
            Self::DegreeMismatch { field_degree, poly_degree } => {
                write!(
                    f,
                    "degree mismatch: field says {field_degree}, f.degree() = {poly_degree:?}"
                )
            }
            Self::InvalidAlgebraicPoly => {
                write!(f, "algebraic polynomial f is zero or constant (degree < 1)")
            }
        }
    }
}

impl std::error::Error for PolyPairError {}

// ─── PolyPair ─────────────────────────────────────────────────────────────────

/// NFS polynomial pair: algebraic-side `f`, rational-side `g = x − m`, shared root `m` mod `n`.
///
/// Invariants (checked by [`verify`](PolyPair::verify)):
///
/// - `f.eval(&m) % &n == 0` — `f` has `m` as a root mod `n`.
/// - `g = x − m` — the rational side is always the linear polynomial `x − m`.
/// - `f.degree() == Some(degree)` — the `degree` field matches the polynomial.
/// - `f` has degree ≥ 1.
///
/// The algebraic polynomial `f` is stored in its *original* form (generally non-monic for
/// base-m expansion). Use [`number_field`](PolyPair::number_field) or
/// [`monic_f`](PolyPair::monic_f) to obtain the monic form required by `NumberField::new`;
/// `number_field` serves the NFS discrete-logarithm relation code and tests only (the sieve and
/// factor base work from the original `f`).
///
/// The `skew` and `factor_base_bounds` fields are `None` at construction; no stage in the
/// library sets either today (they are reserved slots).
#[derive(Debug, Clone)]
pub struct PolyPair {
    /// Algebraic-side polynomial `f ∈ ℤ[x]`. Generally non-monic for base-m; stored as-is.
    pub f: IntPoly,
    /// Rational-side polynomial `g = x − m ∈ ℤ[x]`.
    pub g: IntPoly,
    /// The shared root: `f(m) ≡ 0 (mod n)` and `g(m) = 0`.
    pub m: BigInt,
    /// The integer to factor.
    pub n: BigInt,
    /// Polynomial degree (redundant with `f.degree()` but convenient for pattern matching).
    pub degree: usize,
    /// Skew parameter `s`: the ratio that balances algebraic and rational norm sizes.
    ///
    /// `None` at construction; [`score`] reads it when set (default 1.0). No stage sets it today.
    pub skew: Option<f64>,
    /// Factor-base bounds `(rational_bound, algebraic_bound)`.
    ///
    /// `None` at construction; no stage sets it today.
    pub factor_base_bounds: Option<(u64, u64)>,
}

impl PolyPair {
    /// Construct a new polynomial pair without verifying invariants.
    ///
    /// The `degree` field is inferred from `f.degree()`. Panics if `f` is the zero polynomial.
    /// Call [`verify`](Self::verify) after construction to check all invariants.
    ///
    /// :param f: Algebraic-side polynomial.
    /// :param g: Rational-side polynomial (should be `x − m`).
    /// :param m: Shared root.
    /// :param n: Integer to factor.
    /// :returns: A new `PolyPair` with `skew = None` and `factor_base_bounds = None`.
    pub fn new(f: IntPoly, g: IntPoly, m: BigInt, n: BigInt) -> Self {
        let degree = f.degree().expect("PolyPair::new: f must not be the zero polynomial");
        Self { f, g, m, n, degree, skew: None, factor_base_bounds: None }
    }

    /// Verify the polynomial-pair invariants.
    ///
    /// Checks:
    ///
    /// 1. `f` has degree ≥ 1 (not zero or constant).
    /// 2. `f.degree() == Some(self.degree)` — degree field is consistent.
    /// 3. `g = x − m` — the rational polynomial is the expected linear form.
    /// 4. `f(m) ≡ 0 (mod n)` — the algebraic polynomial has `m` as a root mod `n`.
    ///
    /// :returns: `Ok(())` if all invariants hold, `Err(PolyPairError)` otherwise.
    pub fn verify(&self) -> Result<(), PolyPairError> {
        // Check 1 & 2: degree consistency.
        match self.f.degree() {
            None => return Err(PolyPairError::InvalidAlgebraicPoly),
            Some(0) => return Err(PolyPairError::InvalidAlgebraicPoly),
            Some(d) if d != self.degree => {
                return Err(PolyPairError::DegreeMismatch {
                    field_degree: self.degree,
                    poly_degree: Some(d),
                });
            }
            Some(_) => {}
        }

        // Check 3: g = x − m.
        // g(x) = x − m has coeffs [-m, 1] (least-significant first).
        let expected_g = IntPoly::from_coeffs(vec![-self.m.clone(), BigInt::one()]);
        if self.g != expected_g {
            return Err(PolyPairError::RationalPolyMismatch {
                expected_g,
                actual_g: self.g.clone(),
            });
        }

        // Check 4: f(m) ≡ 0 (mod n).
        let f_of_m = self.f.eval(&self.m);
        if !f_of_m.is_zero() && (&f_of_m % &self.n) != BigInt::zero() {
            return Err(PolyPairError::RootCheckFailed { f_of_m, n: self.n.clone() });
        }

        Ok(())
    }

    /// Return the monic form of `f` via homogenisation, without constructing the full `NumberField`.
    ///
    /// For non-monic `f` with leading coefficient `a_d`, the standard NFS homogenisation is:
    ///
    /// ```text
    /// f_monic(x) = a_d^{d−1} · f(x / a_d)
    /// ```
    ///
    /// The coefficient of `x^k` in `f_monic` is `a_k · a_d^{d−1−k}`, which is always an integer.
    /// The resulting polynomial is monic (leading coefficient 1).
    ///
    /// If `f` is already monic, returns a clone of `f`.
    ///
    /// :returns: The monic form of `f`.
    pub fn monic_f(&self) -> IntPoly {
        let d = self.degree;
        let a_d = self.f.leading_coeff().expect("f must be non-zero").clone();

        if a_d.is_one() {
            return self.f.clone();
        }

        // Compute a_d^{d-1-k} for each k from 0 to d.
        // Coefficient of x^k in f_monic is a_k * a_d^{d-1-k}.
        let mut monic_coeffs = Vec::with_capacity(d + 1);
        for k in 0..=d {
            let a_k = self.f.coeffs.get(k).cloned().unwrap_or_else(BigInt::zero);
            let exp = (d - 1).saturating_sub(k); // d-1-k, but k can equal d giving exp=0 (monic)
            // When k == d: exp = d-1-(d) would underflow; saturating_sub gives 0, and a_d^0 = 1.
            // That gives a_d * 1 = a_d for the leading term, but we want 1 (monic).
            // Handle k == d separately: the leading coefficient of f_monic is always 1.
            let coeff = if k == d {
                BigInt::one()
            } else {
                let power = pow_bigint(&a_d, exp);
                a_k * power
            };
            monic_coeffs.push(coeff);
        }

        IntPoly::from_coeffs(monic_coeffs)
    }

    /// Construct the number field `K = ℚ(α)` where `α` is a root of the monic form of `f`.
    ///
    /// Performs the standard homogenisation `f(x) → a_d^{d−1} f(x/a_d)` to produce a monic
    /// polynomial for `NumberField::new`. The roots of `f_monic` are `a_d · α_i` where `α_i`
    /// are the roots of `f`.
    ///
    /// This is the seam between polynomial selection (which produces non-monic `f`) and
    /// number-field arithmetic (which requires monic `f`). The field returned is built on the
    /// order ℤ[a_d·α], not ℤ[α]. Its consumers are the NFS discrete-logarithm relation code,
    /// which takes the field as a parameter (see `dl::relation`), and tests. The factor base does
    /// its own root search on the original `f`; the algebraic value `F(a, b)` is computed from
    /// the original `f`, not from this field; and the square-root stage builds its own monic
    /// form and never calls this method.
    ///
    /// :returns: The number field `K = ℚ(α)` defined by the monic form of `f`.
    pub fn number_field(&self) -> NumberField {
        NumberField::new(self.monic_f())
    }
}

// ─── PolyGenerator ───────────────────────────────────────────────────────────

/// A polynomial generator produces candidate `PolyPair`s for scoring and ranking.
///
/// All generators (base-m, root sieve, Coppersmith) implement this trait, feeding a common
/// score-and-rank pipeline. The generator is responsible for producing candidates; the
/// Murphy-E scorer ranks them.
pub trait PolyGenerator {
    /// Generate polynomial-pair candidates.
    ///
    /// Returns an iterator of `PolyPair` values. The iterator may be finite (base-m produces
    /// exactly one candidate per `(n, d)` pair) or unbounded (root sieve searches a grid).
    /// Callers should use `.take(limit)` or score-based early termination.
    fn generate(&self) -> impl Iterator<Item = PolyPair>;
}

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Compute `base^exp` for `BigInt` using repeated squaring.
pub(crate) fn pow_bigint(base: &BigInt, exp: usize) -> BigInt {
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

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn poly_pair_verify_ok() {
        // f(x) = x^2 + 1, g(x) = x - 0, m = 0, n = 1 (trivial: f(0) = 1 ≡ 0 mod 1)
        // Use a real example: N = 15, m = 3, d = 2
        // 15 in base 3: 15 = 0 + 2*3 + 1*9, so f(x) = 2x + x^2, f(3) = 6 + 9 = 15
        let n = bi(15);
        let m = bi(3);
        let f = IntPoly::from_coeffs(vec![bi(0), bi(2), bi(1)]);
        let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
        let pair = PolyPair::new(f, g, m, n);
        assert_eq!(pair.verify(), Ok(()));
    }

    #[test]
    fn poly_pair_verify_root_check_fails() {
        // f(x) = x^2 + 1, g(x) = x - 2, m = 2, n = 7
        // f(2) = 5, 5 % 7 != 0
        let n = bi(7);
        let m = bi(2);
        let f = IntPoly::from_coeffs(vec![bi(1), bi(0), bi(1)]);
        let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
        let pair = PolyPair::new(f, g, m, n);
        assert!(matches!(pair.verify(), Err(PolyPairError::RootCheckFailed { .. })));
    }

    #[test]
    fn monic_f_already_monic() {
        // f(x) = x^2 + 3x + 1 (monic), monic_f should return a clone
        let n = bi(15);
        let m = bi(3);
        let f = IntPoly::from_coeffs(vec![bi(1), bi(3), bi(1)]);
        let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
        let pair = PolyPair::new(f.clone(), g, m, n);
        assert_eq!(pair.monic_f(), f);
    }

    #[test]
    fn monic_f_non_monic() {
        // f(x) = 2x^2 + 3x + 1, a_d = 2, d = 2
        // f_monic(x) = a_d^{d-1} f(x/a_d) = 2^1 f(x/2)
        //            = 2 * (2*(x/2)^2 + 3*(x/2) + 1)
        //            = 2 * (x^2/2 + 3x/2 + 1)
        //            = x^2 + 3x + 2
        // Coefficients: k=0: a_0 * a_d^{d-1-0} = 1 * 2^1 = 2
        //               k=1: a_1 * a_d^{d-1-1} = 3 * 2^0 = 3
        //               k=2: 1 (monic)
        let n = bi(100);
        let m = bi(5);
        let f = IntPoly::from_coeffs(vec![bi(1), bi(3), bi(2)]);
        let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
        let pair = PolyPair::new(f, g, m, n);
        let mf = pair.monic_f();
        assert_eq!(mf.coeffs, vec![bi(2), bi(3), bi(1)]);
        assert!(mf.is_monic());
    }
}
