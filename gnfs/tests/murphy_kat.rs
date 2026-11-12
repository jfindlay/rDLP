//! Known-answer tests (KATs) for Murphy-E polynomial scoring.
//!
//! Murphy-E is a heuristic float, so the KATs test *ordering and self-consistency*
//! rather than exact values:
//!
//! 1. **Ordering KAT:** a polynomial pair with smaller coefficients (better) scores
//!    higher than one with larger coefficients (worse) for the same `N`.
//! 2. **Monotonicity KAT:** scaling all coefficients of `f` by a large constant
//!    decreases the score (larger norms → less smooth → lower Murphy-E).
//! 3. **Positivity KAT:** `score(pair) > 0.0` for any valid `PolyPair`.
//! 4. **Dickman ρ unit tests:** spot-checks at known values.
//!
//! # Science↔engineering note
//!
//! Murphy-E's predictive value (higher E → more relations) only manifests at sieve
//! scale (N ≳ 2^100). At toy scale, these KATs verify the ordering property and
//! self-consistency of the implementation, not the absolute values.

use gnfs::{score, select_base_m, PolyPair};
use gnfs::polyselect::murphy::dickman_rho;
use num_bigint::BigInt;
use shared_numfield::IntPoly;

// ─── helpers ─────────────────────────────────────────────────────────────────

fn bi(n: i64) -> BigInt {
    BigInt::from(n)
}

// ─── KAT 1: Ordering ─────────────────────────────────────────────────────────

/// KAT 1: A polynomial pair with smaller coefficients scores higher than one with
/// larger coefficients for the same `N`.
///
/// We use `N = 1009 * 1013 = 1022117` (a toy 20-bit semiprime). The base-m pair
/// is the "better" candidate; we construct a "worse" pair by multiplying all
/// coefficients of `f` by 50 (which inflates the norms and reduces smoothness).
///
/// The base-m pair has `f(m) = N` exactly, so it is a valid pair. The scaled pair
/// has `f_scaled(m) = 50 * N`, which is divisible by `N`, so `verify()` holds.
#[test]
fn kat1_ordering_better_scores_higher() {
    let n = bi(1009 * 1013); // = 1022117
    let better = select_base_m(&n, 3);

    // Construct a "worse" pair: scale f by 50.
    // f_worse(m) = 50 * f(m) = 50 * N ≡ 0 (mod N). ✓
    let f_worse_coeffs: Vec<BigInt> =
        better.f.coeffs.iter().map(|c| c * bi(50)).collect();
    let f_worse = IntPoly::from_coeffs(f_worse_coeffs);
    let g_worse = IntPoly::from_coeffs(vec![-better.m.clone(), bi(1)]);
    let worse = PolyPair::new(f_worse, g_worse, better.m.clone(), n.clone());

    // Both pairs must be valid.
    better.verify().expect("better pair should verify");
    worse.verify().expect("worse pair should verify");

    let score_better = score(&better);
    let score_worse = score(&worse);

    assert!(
        score_better > score_worse,
        "better pair (smaller coefficients) should score higher: \
         score(better) = {score_better}, score(worse) = {score_worse}"
    );
}

// ─── KAT 2: Monotonicity ─────────────────────────────────────────────────────

/// KAT 2: Scaling all coefficients of `f` by a large constant decreases the score.
///
/// For the same `N` and `m`, multiplying `f` by `k` multiplies the algebraic value `F(a, b)`
/// by `k` at every `(a, b)`, increasing `u_alg = log|F| / log B_f` and thus decreasing `ρ(u_alg)`.
/// The rational norm is unaffected. Therefore `score` must decrease as `k` grows.
#[test]
fn kat2_monotonicity_scaling_decreases_score() {
    let n = bi(1009 * 1013);
    let base = select_base_m(&n, 3);

    // score(1x) > score(10x) > score(100x)
    let score_1x = score(&base);

    let f_10x = IntPoly::from_coeffs(base.f.coeffs.iter().map(|c| c * bi(10)).collect());
    let g_10x = IntPoly::from_coeffs(vec![-base.m.clone(), bi(1)]);
    let pair_10x = PolyPair::new(f_10x, g_10x, base.m.clone(), n.clone());
    let score_10x = score(&pair_10x);

    let f_100x = IntPoly::from_coeffs(base.f.coeffs.iter().map(|c| c * bi(100)).collect());
    let g_100x = IntPoly::from_coeffs(vec![-base.m.clone(), bi(1)]);
    let pair_100x = PolyPair::new(f_100x, g_100x, base.m.clone(), n.clone());
    let score_100x = score(&pair_100x);

    assert!(
        score_1x > score_10x,
        "score should decrease when f is scaled by 10: \
         score(1x) = {score_1x}, score(10x) = {score_10x}"
    );
    assert!(
        score_10x > score_100x,
        "score should decrease when f is scaled by 100: \
         score(10x) = {score_10x}, score(100x) = {score_100x}"
    );
}

// ─── KAT 3: Positivity ───────────────────────────────────────────────────────

/// KAT 3: `score(pair) > 0.0` for any valid `PolyPair`.
///
/// Murphy-E is a sum of non-negative terms (products of ρ values). It is zero
/// only if every sample point has a norm so large that ρ = 0, which cannot
/// happen for a polynomial of bounded degree over a bounded sieve region.
#[test]
fn kat3_positivity() {
    // Toy N.
    let n_toy = bi(1009 * 1013);
    let pair_toy = select_base_m(&n_toy, 3);
    assert!(
        score(&pair_toy) > 0.0,
        "score should be positive for toy N pair, got {}",
        score(&pair_toy)
    );

    // Slightly larger N.
    let n_med = bi(999983i64 * 999979i64);
    let pair_med = select_base_m(&n_med, 3);
    assert!(
        score(&pair_med) > 0.0,
        "score should be positive for medium N pair, got {}",
        score(&pair_med)
    );

    // A manually constructed pair with small coefficients.
    // N = 35 = 5 * 7, m = 5, f(x) = x + 2 (f(5) = 7, not 35 — use base-m instead).
    // Use select_base_m for correctness.
    let n_small = bi(35);
    let pair_small = select_base_m(&n_small, 2);
    pair_small.verify().expect("small pair should verify");
    assert!(
        score(&pair_small) > 0.0,
        "score should be positive for small N pair, got {}",
        score(&pair_small)
    );
}

// ─── KAT 4: Dickman ρ unit tests ─────────────────────────────────────────────

/// KAT 4a: ρ(0.5) = 1.0 (u ≤ 1 region).
#[test]
fn kat4a_dickman_rho_half() {
    assert_eq!(dickman_rho(0.5), 1.0, "ρ(0.5) should be exactly 1.0");
}

/// KAT 4b: ρ(1.0) = 1.0 (boundary of u ≤ 1 region).
#[test]
fn kat4b_dickman_rho_one() {
    assert_eq!(dickman_rho(1.0), 1.0, "ρ(1.0) should be exactly 1.0");
}

/// KAT 4c: ρ(1.5) ≈ 1 − ln(1.5) ≈ 0.5945.
#[test]
fn kat4c_dickman_rho_1_5() {
    let expected = 1.0 - 1.5_f64.ln(); // ≈ 0.5945
    let got = dickman_rho(1.5);
    assert!(
        (got - expected).abs() < 1e-12,
        "ρ(1.5) expected {expected:.6}, got {got:.6}"
    );
}

/// KAT 4d: ρ(2.0) = 1 − ln(2) ≈ 0.3069.
#[test]
fn kat4d_dickman_rho_two() {
    let expected = 1.0 - 2.0_f64.ln(); // ≈ 0.3069
    let got = dickman_rho(2.0);
    assert!(
        (got - expected).abs() < 1e-12,
        "ρ(2.0) expected {expected:.6}, got {got:.6}"
    );
}

/// KAT 4e: ρ(10.0) is very small (< 1e-6).
///
/// The exact value is approximately 2.77 × 10^{-10} (from tables), but we only
/// require it to be below 1e-6 to avoid over-constraining the approximation.
#[test]
fn kat4e_dickman_rho_ten() {
    let got = dickman_rho(10.0);
    assert!(
        got < 1e-6,
        "ρ(10.0) should be < 1e-6 (essentially zero), got {got:.2e}"
    );
    assert!(got >= 0.0, "ρ(10.0) should be non-negative, got {got}");
}

// ─── KAT 5: Identity invariance ──────────────────────────────────────────────

/// KAT 5: `score` is deterministic — calling it twice on the same pair gives the
/// same result (no hidden mutable state).
#[test]
fn kat5_score_is_deterministic() {
    let n = bi(1009 * 1013);
    let pair = select_base_m(&n, 3);
    let s1 = score(&pair);
    let s2 = score(&pair);
    assert_eq!(s1, s2, "score should be deterministic");
}

// ─── KAT 6: Ordering with a genuinely better polynomial ──────────────────────

/// KAT 6: A degree-3 polynomial with small, balanced coefficients scores higher
/// than one with large, unbalanced coefficients for the same `N` and `m`.
///
/// We construct two pairs for `N = 1022117`, `m = 101` (a valid base):
///
/// - `f_good(x) = x^3 + x^2 - x + 1` — small coefficients,
///   `f_good(101) = 1030503 + 10201 - 101 + 1 = 1040604`.
///   This does not equal N, so we use a pair where `f(m) = N` exactly.
///
/// Instead, we compare the base-m pair (which is the "natural" good polynomial)
/// against a pair where `f` has been perturbed to have much larger coefficients
/// while still satisfying `f(m) ≡ 0 (mod N)`.
///
/// Perturbation: `f_perturbed = f_base + N * x^0` (add N to the constant term).
/// Then `f_perturbed(m) = f_base(m) + N = N + N = 2N ≡ 0 (mod N)`. ✓
/// The constant term is now `N` times larger, inflating the algebraic values `F(a, b)`.
#[test]
fn kat6_perturbed_polynomial_scores_lower() {
    let n = bi(1009 * 1013); // = 1022117
    let base = select_base_m(&n, 3);

    // Perturb: add N to the constant coefficient.
    let mut perturbed_coeffs = base.f.coeffs.clone();
    perturbed_coeffs[0] += n.clone();
    let f_perturbed = IntPoly::from_coeffs(perturbed_coeffs);
    let g_perturbed = IntPoly::from_coeffs(vec![-base.m.clone(), bi(1)]);
    let perturbed = PolyPair::new(f_perturbed, g_perturbed, base.m.clone(), n.clone());

    // Verify both pairs.
    base.verify().expect("base pair should verify");
    perturbed.verify().expect("perturbed pair should verify");

    let score_base = score(&base);
    let score_perturbed = score(&perturbed);

    assert!(
        score_base > score_perturbed,
        "base pair should score higher than perturbed pair: \
         score(base) = {score_base}, score(perturbed) = {score_perturbed}"
    );
}
