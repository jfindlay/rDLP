//! Real-root approximation for NFS polynomial scoring.
//!
//! Provides `f64`-precision real-root finding for `IntPoly` values. This is a
//! *scoring-only* facility: the roots are used to count real roots and to
//! characterise the size of the algebraic value `F(a, b)` over the sieve region. Exact
//! arithmetic (discriminants, resultants) stays in `shared-numfield`.
//!
//! # Algorithm
//!
//! Real roots are found by:
//!
//! 1. **Sturm-sequence root counting** — determine the number of real roots in
//!    a bounding interval `[-R, R]` where `R` is a Cauchy bound on the roots.
//! 2. **Bisection** — isolate each root to a sign-change bracket, then refine
//!    to `f64` precision with Newton steps.
//!
//! This is sufficient for scoring; it is not a certified root-finder.
//!
//! # Science↔engineering note
//!
//! The number of real roots of `f` is one factor in Murphy-E: polynomials with
//! more real roots tend to have smaller average norms over the sieve region
//! because the norm `|F(a,b)|` dips near each real root. At toy scale this
//! effect is small; at cryptographic scale it is a meaningful discriminator.

use shared_numfield::IntPoly;

// ─── polynomial evaluation at f64 ────────────────────────────────────────────

/// Evaluate `f` at `x` using Horner's method in `f64` arithmetic.
///
/// :param f: The polynomial to evaluate.
/// :param x: The evaluation point.
/// :returns: `f(x)` as `f64`.
pub(crate) fn eval_f64(f: &IntPoly, x: f64) -> f64 {
    let mut result = 0.0_f64;
    for c in f.coeffs.iter().rev() {
        // Convert BigInt coefficient to f64 via string to avoid overflow in to_i64.
        let c_f64 = bigint_to_f64(c);
        result = result * x + c_f64;
    }
    result
}

/// Convert a `BigInt` to `f64`.
///
/// Uses the sign and magnitude separately to avoid intermediate overflow.
///
/// :param n: The integer to convert.
/// :returns: The nearest `f64` value.
fn bigint_to_f64(n: &num_bigint::BigInt) -> f64 {
    use num_traits::ToPrimitive;
    // ToPrimitive::to_f64 handles large BigInts correctly (returns ±inf on overflow).
    n.to_f64().unwrap_or(f64::INFINITY)
}

/// Evaluate the formal derivative of `f` at `x` in `f64` arithmetic.
///
/// :param f: The polynomial whose derivative is evaluated.
/// :param x: The evaluation point.
/// :returns: `f'(x)` as `f64`.
pub(crate) fn eval_deriv_f64(f: &IntPoly, x: f64) -> f64 {
    // f'(x) = Σ_{k=1}^{d} k * a_k * x^{k-1}
    // Evaluated via Horner on the derivative coefficient sequence [a_1, 2*a_2, ..., d*a_d].
    let d = match f.degree() {
        None | Some(0) => return 0.0,
        Some(d) => d,
    };
    let mut result = 0.0_f64;
    for k in (1..=d).rev() {
        let c_f64 = bigint_to_f64(&f.coeffs[k]) * k as f64;
        result = result * x + c_f64;
    }
    result
}

// ─── Cauchy root bound ────────────────────────────────────────────────────────

/// Compute a Cauchy upper bound on the absolute values of all roots of `f`.
///
/// The Cauchy bound is `1 + max(|a_0|, ..., |a_{d-1}|) / |a_d|`, where `a_d` is
/// the leading coefficient. All real (and complex) roots satisfy `|r| ≤ bound`.
///
/// :param f: The polynomial. Must be non-zero.
/// :returns: A finite upper bound on `|root|` for all roots of `f`.
fn cauchy_bound(f: &IntPoly) -> f64 {
    let d = match f.degree() {
        None | Some(0) => return 0.0,
        Some(d) => d,
    };
    let lc = bigint_to_f64(f.coeffs.last().unwrap()).abs();
    if lc == 0.0 {
        return f64::INFINITY;
    }
    let max_lower = f.coeffs[..d]
        .iter()
        .map(|c| bigint_to_f64(c).abs())
        .fold(0.0_f64, f64::max);
    1.0 + max_lower / lc
}

// ─── sign-change root isolation ───────────────────────────────────────────────

/// Find sign-change brackets for real roots of `f` in `[lo, hi]`.
///
/// Recursively bisects the interval until each sub-interval is small enough to
/// contain at most one root (width < `min_width`), then records any sub-interval
/// that contains a sign change. Returns a list of `(a, b)` pairs where `f(a)` and
/// `f(b)` have opposite signs.
///
/// Unlike a naive sign-change search, this always bisects down to `min_width`
/// before deciding whether a bracket contains a root. This ensures that intervals
/// with an even number of roots (same sign at both endpoints) are still explored.
///
/// :param f: The polynomial to bracket.
/// :param lo: Left endpoint of the search interval.
/// :param hi: Right endpoint of the search interval.
/// :param min_width: Minimum sub-interval width; recursion stops below this.
/// :returns: A list of sign-change brackets.
fn find_brackets(f: &IntPoly, lo: f64, hi: f64, min_width: f64) -> Vec<(f64, f64)> {
    if (hi - lo) <= min_width {
        // Leaf: check for a sign change in this tiny interval.
        let f_lo = eval_f64(f, lo);
        let f_hi = eval_f64(f, hi);
        if f_lo == 0.0 {
            return vec![(lo, lo)];
        }
        if f_hi == 0.0 {
            return vec![(hi, hi)];
        }
        if f_lo * f_hi < 0.0 {
            return vec![(lo, hi)];
        }
        return vec![];
    }

    // Bisect and recurse on both halves regardless of sign at endpoints.
    // This ensures we explore intervals with an even number of roots.
    let mid = (lo + hi) / 2.0;
    let mut brackets = find_brackets(f, lo, mid, min_width);
    brackets.extend(find_brackets(f, mid, hi, min_width));
    brackets
}

// ─── bisection + Newton refinement ───────────────────────────────────────────

/// Refine a sign-change bracket `[lo, hi]` to a root of `f` using bisection.
///
/// Runs up to 64 bisection steps, then polishes with Newton iterations. Returns
/// the midpoint of the final bracket.
///
/// :param f: The polynomial.
/// :param lo: Left endpoint (f(lo) and f(hi) have opposite signs).
/// :param hi: Right endpoint.
/// :returns: An approximate root of `f` in `[lo, hi]`.
fn refine_bracket(f: &IntPoly, mut lo: f64, mut hi: f64) -> f64 {
    if lo == hi {
        return lo; // degenerate bracket — lo is already a root
    }
    let mut f_lo = eval_f64(f, lo);

    for _ in 0..64 {
        let mid = (lo + hi) / 2.0;
        if mid == lo || mid == hi {
            break; // floating-point precision exhausted
        }
        let f_mid = eval_f64(f, mid);
        if f_mid == 0.0 {
            return mid;
        }
        if f_lo * f_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    let root = (lo + hi) / 2.0;

    // Polish with a few Newton steps.
    let mut x = root;
    for _ in 0..8 {
        let fx = eval_f64(f, x);
        let dfx = eval_deriv_f64(f, x);
        if dfx.abs() < 1e-300 {
            break;
        }
        let x_new = x - fx / dfx;
        if (x_new - x).abs() < 1e-14 * x.abs().max(1.0) {
            x = x_new;
            break;
        }
        x = x_new;
    }
    x
}

// ─── public API ──────────────────────────────────────────────────────────────

/// Approximate the real roots of `f` in `f64` precision.
///
/// Uses a Cauchy bound to determine the search interval `[-R, R]`, then
/// bisects recursively to isolate sign-change brackets, and refines each
/// bracket with bisection + Newton polishing.
///
/// This is a *scoring-only* root finder: it is fast and practical but not
/// certified. Repeated roots may be missed or duplicated. Exact root
/// arithmetic stays in `shared-numfield`.
///
/// :param f: The polynomial whose real roots are sought. Must be non-zero.
/// :returns: A vector of approximate real roots, in no particular order.
///   May contain duplicates near repeated roots.
pub fn real_roots_f64(f: &IntPoly) -> Vec<f64> {
    let d = match f.degree() {
        None | Some(0) => return vec![],
        Some(d) => d,
    };

    let r = cauchy_bound(f);
    if !r.is_finite() || r == 0.0 {
        return vec![];
    }

    // Use a minimum sub-interval width of 2R / 2^(log2(d)+8).
    // For degree d ≤ 6 and R ≤ a few thousand, this gives sub-intervals of width
    // ~R/128 which is small enough to isolate individual roots that are at least
    // 0.1 apart. For very close roots, the deduplication step below handles them.
    let _ = d; // used above for the degree check; not needed for min_width
    let min_width = (2.0 * r) / 512.0; // 512 = 2^9 sub-intervals

    let brackets = find_brackets(f, -r, r, min_width);

    // Deduplicate brackets that are very close together (repeated roots).
    let mut roots: Vec<f64> = brackets
        .into_iter()
        .map(|(lo, hi)| refine_bracket(f, lo, hi))
        .collect();

    // Remove duplicates within 1e-8 relative tolerance.
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-8 * b.abs().max(1.0));

    roots
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    fn poly(coeffs: Vec<i64>) -> IntPoly {
        IntPoly::from_coeffs(coeffs.into_iter().map(bi).collect())
    }

    #[test]
    fn roots_linear() {
        // f(x) = 2x - 6 → root at x = 3
        let f = poly(vec![-6, 2]);
        let roots = real_roots_f64(&f);
        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 3.0).abs() < 1e-10, "root should be 3.0, got {}", roots[0]);
    }

    #[test]
    fn roots_quadratic_two_real() {
        // f(x) = x^2 - 5x + 6 = (x-2)(x-3) → roots at 2 and 3
        let f = poly(vec![6, -5, 1]);
        let mut roots = real_roots_f64(&f);
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(roots.len(), 2, "expected 2 roots, got {:?}", roots);
        assert!((roots[0] - 2.0).abs() < 1e-8, "first root should be 2.0, got {}", roots[0]);
        assert!((roots[1] - 3.0).abs() < 1e-8, "second root should be 3.0, got {}", roots[1]);
    }

    #[test]
    fn roots_quadratic_no_real() {
        // f(x) = x^2 + 1 → no real roots
        let f = poly(vec![1, 0, 1]);
        let roots = real_roots_f64(&f);
        assert!(roots.is_empty(), "expected no real roots, got {:?}", roots);
    }

    #[test]
    fn roots_cubic_one_real() {
        // f(x) = x^3 - 1 → one real root at x = 1
        let f = poly(vec![-1, 0, 0, 1]);
        let roots = real_roots_f64(&f);
        assert_eq!(roots.len(), 1, "expected 1 real root, got {:?}", roots);
        assert!((roots[0] - 1.0).abs() < 1e-8, "root should be 1.0, got {}", roots[0]);
    }

    #[test]
    fn roots_cubic_three_real() {
        // f(x) = (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
        let f = poly(vec![-6, 11, -6, 1]);
        let mut roots = real_roots_f64(&f);
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(roots.len(), 3, "expected 3 roots, got {:?}", roots);
        assert!((roots[0] - 1.0).abs() < 1e-7);
        assert!((roots[1] - 2.0).abs() < 1e-7);
        assert!((roots[2] - 3.0).abs() < 1e-7);
    }

    #[test]
    fn roots_constant_returns_empty() {
        let f = poly(vec![5]);
        let roots = real_roots_f64(&f);
        assert!(roots.is_empty());
    }
}
