//! Murphy-E polynomial scoring for NFS polynomial pairs.
//!
//! Murphy-E is the standard heuristic score for NFS polynomial pairs. A higher
//! score predicts more smooth relations during sieving, and therefore a faster
//! factorisation. The score is defined as:
//!
//! ```text
//! E(f, g) ≈ (1/|S|) Σ_{(a,b) ∈ S} ρ(log|F(a,b)| / log B_f) · ρ(log|G(a,b)| / log B_g)
//! ```
//!
//! where:
//!
//! - `S` is a sample of `(a, b)` pairs in the sieve region.
//! - `F(a,b) = b^d · f(a/b)` is the algebraic value (homogeneous form of `f`); with `c` the
//!   leading coefficient and `α` a root of `f`, `F(a,b) = c · N(a − b·α)`. Evaluated by
//!   `alg_norm_f64`, which no test calls directly; `score` is exercised by `tests/murphy_kat.rs`
//!   (for example `kat1_ordering_better_scores_higher`).
//! - `G(a,b) = a − b·m` is the rational norm.
//! - `B_f`, `B_g` are the algebraic and rational factor-base bounds.
//! - `ρ` is the Dickman rho function (probability that a random integer near `x`
//!   is `x^{1/u}`-smooth).
//!
//! # Science↔engineering disconnect (principle-4 annotation)
//!
//! Murphy-E's *predictive* value — that higher E implies more relations — only
//! manifests at sieve scale (N ≳ 2^100, sieve region ≳ 10^6 pairs). At toy
//! scale (N < 2^60), the sieve region is tiny and the smoothness probabilities
//! are dominated by the factor-base bounds rather than the polynomial shape.
//! At toy scale, Murphy-E is a *ranking heuristic* whose payoff is
//! under-exposed: the ordering it induces is correct in expectation, but the
//! absolute values are not meaningful. Downstream consumers (root sieve, Coppersmith,
//! NFS-DL polynomial selection) should treat the score as an ordinal, not a cardinal.
//!
//! # References
//!
//! - Murphy, B. (1999). *Polynomial selection for the number field sieve integer
//!   factorisation algorithm*. PhD thesis, ANU.
//! - Bai, S., Bouvier, C., Kruppa, A., Zimmermann, P. (2014). *Better polynomials
//!   for GNFS*. Mathematics of Computation.

use num_traits::ToPrimitive;

use super::PolyPair;

// ─── Dickman ρ ────────────────────────────────────────────────────────────────

/// Approximate the Dickman rho function `ρ(u)`.
///
/// `ρ(u)` is the probability that a random integer near `x` is `x^{1/u}`-smooth
/// (has no prime factor larger than `x^{1/u}`). It satisfies:
///
/// - `ρ(u) = 1` for `u ≤ 1`.
/// - `ρ(u) = 1 − ln u` for `1 < u ≤ 2`.
/// - For `u > 2`: the delay-differential equation `u·ρ'(u) = −ρ(u−1)`.
///
/// This implementation uses:
///
/// - Exact formula for `u ≤ 2`.
/// - A piecewise series expansion for `2 < u ≤ 3`:
///   `ρ(u) = 1 − ln u + (ln u)²/2 − (ln u)³/6 + ... + ∫ correction`.
///   We use the known closed form for `2 < u ≤ 3`:
///   `ρ(u) = 1 − ln u + (1/2)(ln u − ln(u−1))² + Li₂(1 − 1/(u−1))` (approx).
///   In practice we use the recurrence integrated numerically via Simpson's rule.
/// - Numerical integration of the recurrence for `3 < u ≤ 25`.
/// - `ρ(u) = 0` for `u > 25` (essentially zero; `ρ(25) ≈ 10^{-30}`).
///
/// :param u: The smoothness-exponent ratio. Must be non-negative.
/// :returns: An approximation of `ρ(u)` in `[0, 1]`.
pub fn dickman_rho(u: f64) -> f64 {
    if u <= 0.0 {
        return 1.0;
    }
    if u <= 1.0 {
        return 1.0;
    }
    if u <= 2.0 {
        return 1.0 - u.ln();
    }
    if u > 25.0 {
        return 0.0;
    }

    // For u > 2, use the recurrence ρ(u) = (1/u) ∫_{u-1}^{u} ρ(t) dt integrated
    // numerically by building up from u=2 in steps.
    //
    // We use a precomputed table at integer points and interpolate linearly for
    // intermediate values. The table is built by integrating the recurrence:
    //   u·ρ'(u) + ρ(u-1) = 0  ⟹  ρ'(u) = -ρ(u-1)/u
    // using a fine grid (step h = 0.01) from u=2 to u=25.
    dickman_rho_numerical(u)
}

/// Compute `ρ(u)` for `u > 2` by numerically integrating the recurrence.
///
/// The recurrence `ρ'(u) = -ρ(u-1)/u` is integrated from `u=2` (where `ρ(2) = 1 - ln 2`)
/// forward using a 4th-order Runge-Kutta scheme on a uniform grid of step `h = 0.01`.
/// The result at the requested `u` is obtained by linear interpolation between the two
/// nearest grid points.
///
/// The full table (`buf`) is allocated as a `Vec` of 2301 entries (≈18 KB). This avoids
/// the correctness hazard of a rolling window, which would overwrite values still needed
/// for the one-unit lookback `ρ(t-1)`.
///
/// :param u: The evaluation point. Must satisfy `2 < u ≤ 25`.
/// :returns: An approximation of `ρ(u)`.
fn dickman_rho_numerical(u: f64) -> f64 {
    // Grid: t_k = 2.0 + k * H for k = 0, 1, ..., N_STEPS.
    // buf[k] = ρ(t_k).
    const H: f64 = 0.01;
    const N_STEPS: usize = 2300; // t goes from 2.0 to 25.0

    // ρ(t) for t ∈ [1, 2]: exact closed form.
    let rho_exact = |t: f64| -> f64 {
        if t <= 1.0 { 1.0 } else { 1.0 - t.ln() }
    };

    // Allocate the full table. 2301 f64 values = ~18 KB — negligible.
    let mut buf = vec![0.0_f64; N_STEPS + 1];
    buf[0] = rho_exact(2.0); // = 1 - ln 2

    // ρ(t-1) lookup: for t = 2.0 + k*H, t-1 = 1.0 + k*H.
    // For k ≤ 100: t-1 ∈ [1.0, 2.0] → exact formula.
    // For k > 100: t-1 = 2.0 + (k-100)*H → buf[k-100] (always already computed).
    let rho_prev = |k: usize, buf: &[f64]| -> f64 {
        if k <= 100 {
            rho_exact(1.0 + k as f64 * H)
        } else {
            buf[k - 100]
        }
    };

    // For RK4 intermediate points at t + H/2 and t + H, we need ρ(t-1+H/2) and ρ(t-1+H).
    // These fall between grid points; interpolate linearly between buf[k-100] and buf[k-99].
    let rho_prev_half = |k: usize, buf: &[f64]| -> f64 {
        // t + H/2 - 1 = 1.0 + k*H + H/2 - 1.0 = k*H + H/2
        // = (k + 0.5) * H above 1.0
        // For k < 100: in [1, 2] range → exact at t = 1.0 + (k+0.5)*H
        if k < 100 {
            rho_exact(1.0 + (k as f64 + 0.5) * H)
        } else {
            // t-1+H/2 = 2.0 + (k-100)*H + H/2 → midpoint between buf[k-100] and buf[k-99]
            let lo = buf[k - 100];
            let hi = if k - 99 <= N_STEPS { buf[k - 99] } else { 0.0 };
            (lo + hi) / 2.0
        }
    };

    let rho_prev_full = |k: usize, buf: &[f64]| -> f64 {
        // t + H - 1 = 1.0 + (k+1)*H
        if k + 1 <= 100 {
            rho_exact(1.0 + (k + 1) as f64 * H)
        } else {
            buf[k + 1 - 100]
        }
    };

    for k in 0..N_STEPS {
        let t = 2.0 + k as f64 * H;
        let rho_k = buf[k];

        // RK4 for ρ'(t) = -ρ(t-1)/t.
        let k1 = -rho_prev(k, &buf) / t;
        let k2 = -rho_prev_half(k, &buf) / (t + H / 2.0);
        let k3 = -rho_prev_half(k, &buf) / (t + H / 2.0); // same as k2 (autonomous in ρ)
        let k4 = -rho_prev_full(k, &buf) / (t + H);

        let rho_next = (rho_k + H / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)).max(0.0);
        buf[k + 1] = rho_next;
    }

    // Linear interpolation between the two nearest grid points.
    let k_lo = ((u - 2.0) / H) as usize;
    let k_lo = k_lo.min(N_STEPS - 1);
    let k_hi = k_lo + 1;
    let frac = (u - 2.0 - k_lo as f64 * H) / H;
    (buf[k_lo] + frac * (buf[k_hi] - buf[k_lo])).max(0.0)
}

// ─── homogeneous value evaluation ────────────────────────────────────────────

/// Evaluate the homogeneous algebraic value `F(a, b) = b^d · f(a/b)`, which equals
/// `c · N(a − b·α)` for leading coefficient `c` and root `α` (not the field norm itself unless
/// `f` is monic).
///
/// Computed directly as `Σ_k c_k · a^k · b^{d-k}` to avoid division. It has no direct test; it
/// is exercised only through [`score`] by `tests/murphy_kat.rs`
/// (`kat1_ordering_better_scores_higher`).
///
/// :param f: The algebraic polynomial (degree `d`).
/// :param a: The `a`-coordinate.
/// :param b: The `b`-coordinate. Must be non-zero.
/// :returns: `F(a, b)` as `f64`.
fn alg_norm_f64(f: &shared_numfield::IntPoly, a: f64, b: f64) -> f64 {
    // F(a,b) = Σ_{k=0}^{d} c_k * a^k * b^{d-k}
    let d = match f.degree() {
        None => return 0.0,
        Some(d) => d,
    };

    let mut result = 0.0_f64;
    let mut a_pow = 1.0_f64; // a^k
    let mut b_pow = {
        // b^d
        let mut p = 1.0_f64;
        for _ in 0..d {
            p *= b;
        }
        p
    };

    for k in 0..=d {
        let c_k = coeff_f64(f, k);
        result += c_k * a_pow * b_pow;
        a_pow *= a;
        if b != 0.0 {
            b_pow /= b;
        }
    }
    result
}

/// Extract coefficient `k` of `f` as `f64`.
fn coeff_f64(f: &shared_numfield::IntPoly, k: usize) -> f64 {
    use num_traits::ToPrimitive;
    f.coeffs.get(k).and_then(|c| c.to_f64()).unwrap_or(0.0)
}

/// Evaluate the rational norm `G(a, b) = a − b·m`.
///
/// :param a: The `a`-coordinate.
/// :param b: The `b`-coordinate.
/// :param m: The shared root (as `f64`).
/// :returns: `G(a, b) = a - b*m`.
fn rat_norm_f64(a: f64, b: f64, m: f64) -> f64 {
    a - b * m
}

// ─── Murphy-E score ───────────────────────────────────────────────────────────

/// Compute the Murphy-E score for an NFS polynomial pair.
///
/// Murphy-E is a heuristic score predicting the density of smooth relations
/// produced by the NFS sieve. It is defined as the average over a sample of
/// `(a, b)` pairs in the sieve region of the product of Dickman-ρ values for
/// the algebraic value `F` and the rational norm `G`:
///
/// ```text
/// E(f, g) ≈ (1/|S|) Σ_{(a,b) ∈ S} ρ(log|F(a,b)| / log B_f) · ρ(log|G(a,b)| / log B_g)
/// ```
///
/// **Science↔engineering disconnect:** Murphy-E's predictive value only
/// manifests at sieve scale (N ≳ 2^100). At toy scale it is a ranking
/// heuristic: the ordering is correct in expectation, but the absolute values
/// are not meaningful. Downstream consumers (root sieve, Coppersmith, NFS-DL polynomial
/// selection) should treat the score as an ordinal.
///
/// **Skew:** The optimal skew `s` balances algebraic and rational norm sizes.
/// This implementation uses the skew stored in `pair.skew` if present, or
/// defaults to `s = 1.0`. The skew is applied by sampling `(a, b·s)` instead
/// of `(a, b)`, which rescales the algebraic value `F` relative to the rational norm `G`.
/// Since `score` takes `&PolyPair` (not `&mut`), the computed skew is internal
/// and not written back to the pair.
///
/// :param pair: The polynomial pair to score.
/// :returns: The Murphy-E score as a non-negative `f64`. Higher is better.
pub fn score(pair: &PolyPair) -> f64 {
    // Sieve region parameters.
    // M = 1000 gives a 50×50 sample grid (a ∈ [-M, M], b ∈ [1, M]).
    const M: f64 = 1000.0;
    const GRID: usize = 50; // number of sample points per axis

    // Factor-base bounds. At toy scale, use fixed bounds; at cryptographic scale
    // these would be derived from N's size. 10^6 is a reasonable default.
    let log_bf = (1_000_000.0_f64).ln();
    let log_bg = (1_000_000.0_f64).ln();

    // Skew: use pair.skew if set, otherwise default to 1.0.
    let skew = pair.skew.unwrap_or(1.0).max(0.01);

    // m as f64 for rational norm evaluation.
    let m_f64 = pair.m.to_f64().unwrap_or(0.0);

    let a_step = 2.0 * M / GRID as f64;
    let b_step = (M - 1.0) / GRID as f64;

    let mut total = 0.0_f64;
    let mut count = 0usize;

    for i in 0..=GRID {
        let a = -M + i as f64 * a_step;
        for j in 1..=GRID {
            let b_raw = 1.0 + (j - 1) as f64 * b_step;
            let b = b_raw * skew; // apply skew to b

            // Algebraic value: F(a, b) = b^d * f(a/b) = c * N(a - b*alpha).
            // We evaluate using the homogeneous form to avoid division.
            let f_norm = alg_norm_f64(&pair.f, a, b).abs();
            // Rational norm: G(a, b) = a - b_raw * m (skew does not affect rational side).
            let g_norm = rat_norm_f64(a, b_raw, m_f64).abs();

            // Skip degenerate points.
            if f_norm < 1.0 || g_norm < 1.0 {
                count += 1;
                // A norm < 1 means the point is very close to a root — treat as
                // maximally smooth (ρ = 1) to avoid log(0).
                total += 1.0;
                continue;
            }

            let u_alg = f_norm.ln() / log_bf;
            let u_rat = g_norm.ln() / log_bg;

            let rho_alg = dickman_rho(u_alg);
            let rho_rat = dickman_rho(u_rat);

            total += rho_alg * rho_rat;
            count += 1;
        }
    }

    if count == 0 {
        return 0.0;
    }
    total / count as f64
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dickman ρ unit tests ──────────────────────────────────────────────────

    #[test]
    fn dickman_rho_at_zero() {
        assert_eq!(dickman_rho(0.0), 1.0);
    }

    #[test]
    fn dickman_rho_at_one() {
        assert_eq!(dickman_rho(1.0), 1.0);
    }

    #[test]
    fn dickman_rho_at_half() {
        // ρ(0.5) = 1.0 (u ≤ 1)
        assert_eq!(dickman_rho(0.5), 1.0);
    }

    #[test]
    fn dickman_rho_at_1_5() {
        // ρ(1.5) = 1 - ln(1.5) ≈ 0.5945
        let expected = 1.0 - 1.5_f64.ln();
        let got = dickman_rho(1.5);
        assert!(
            (got - expected).abs() < 1e-12,
            "ρ(1.5) expected {expected}, got {got}"
        );
    }

    #[test]
    fn dickman_rho_at_2() {
        // ρ(2.0) = 1 - ln(2) ≈ 0.3069
        let expected = 1.0 - 2.0_f64.ln();
        let got = dickman_rho(2.0);
        assert!(
            (got - expected).abs() < 1e-12,
            "ρ(2.0) expected {expected}, got {got}"
        );
    }

    #[test]
    fn dickman_rho_at_3() {
        // ρ(3) ≈ 0.04860838 (known value from tables)
        // We allow 1% relative tolerance for the numerical approximation.
        let got = dickman_rho(3.0);
        let expected = 0.04860838_f64;
        let rel_err = (got - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "ρ(3.0) expected ≈ {expected}, got {got} (rel err {rel_err:.4})"
        );
    }

    #[test]
    fn dickman_rho_at_10() {
        // ρ(10) is extremely small (< 1e-6).
        let got = dickman_rho(10.0);
        assert!(
            got < 1e-6,
            "ρ(10.0) should be < 1e-6, got {got}"
        );
    }

    #[test]
    fn dickman_rho_beyond_25() {
        assert_eq!(dickman_rho(26.0), 0.0);
        assert_eq!(dickman_rho(100.0), 0.0);
    }

    #[test]
    fn dickman_rho_is_decreasing() {
        // ρ should be monotonically non-increasing.
        let us = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0];
        let mut prev = dickman_rho(us[0]);
        for &u in &us[1..] {
            let cur = dickman_rho(u);
            assert!(
                cur <= prev + 1e-12,
                "ρ should be non-increasing: ρ({u}) = {cur} > ρ(prev) = {prev}"
            );
            prev = cur;
        }
    }
}
