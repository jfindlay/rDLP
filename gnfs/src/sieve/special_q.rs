//! Special-q strategy for GNFS relation collection.
//!
//! The special-q strategy is an optimization layer over the line sieve. For each
//! **special prime** ``q`` in a chosen range ``[q_min, q_max]``, the sieve is restricted to
//! pairs ``(a, b)`` for which ``q | F(a, b)``, where ``F(a, b) = b^d·f(a/b) = c·N_{K/ℚ}(a − b·α)``
//! is the algebraic value (``N_alg`` in code). This restriction is enforced by the sieve
//! condition: ``q | F(a, b)`` iff ``a ≡ r_q·b (mod q)`` for some root ``r_q`` of ``f mod q``
//! (for ``q ∤ b``).
//!
//! # Algorithm
//!
//! For each special prime ``q`` in ``[q_min, q_max]``:
//!
//! 1. Find all roots ``r_q`` of ``f mod q`` (there may be 0, 1, 2, or ``deg(f)`` roots).
//! 2. For each root ``r_q``, run a restricted line sieve:
//!    - For each ``b`` in ``1..=B``, only consider ``a`` values with ``a ≡ r_q·b (mod q)``.
//!    - This means stepping through ``a`` in increments of ``q`` (starting from the first
//!      ``a ≥ −A`` with ``a ≡ r_q·b (mod q)``).
//!    - For each such ``a``, accumulate log contributions from the rational and algebraic
//!      factor bases (as in the line sieve), then trial-divide survivors.
//! 3. Every confirmed relation carries ``q`` in its algebraic exponent vector (since the
//!    sieve condition guarantees ``q | F(a, b)``).
//!
//! # Yield multiplier
//!
//! The key insight: by restricting to ``a ≡ r_q·b (mod q)``, every candidate already has
//! ``q`` as a known algebraic factor. The remaining cofactor ``F(a, b) / q`` is smaller
//! and therefore more likely to be smooth over the algebraic factor base. This increases the
//! probability that the remaining cofactor is smooth, yielding more relations per sieve area.
//!
//! # Toy scale versus cryptographic scale
//!
//! At toy scale (small ``A``, ``B``, ``q``), the yield advantage of the special-q strategy
//! over the plain line sieve is under-exposed. The yield multiplier is a scale phenomenon:
//! at cryptographic scale, the special-q strategy is the dominant sieving technique because
//! the algebraic value ``F(a, b)`` is large and the probability of smoothness is low
//! without the pre-guaranteed factor ``q``. At toy scale, the values are already small and
//! smooth, so the advantage is marginal. The KATs and ``gnfs/docs/PEDAGOGY.md`` state the
//! same caveat.
//!
//! # Relation to the line sieve
//!
//! The special-q sieve is a variant of the line sieve with an additional constraint on ``a``.
//! It reuses the same log-sieve + trial-divide pattern, the same ``FactorBase``, and the same
//! ``Relation`` type. The ``q``-restriction is implemented by stepping ``a`` in increments of
//! ``q`` rather than 1, and by pre-guaranteeing ``q`` in the algebraic exponent vector.
//!
//! The plain line sieve remains available for benchmarking: the baseline is not altered by the
//! optimization layer.
//!
//! # Entry surface
//!
//! - [`special_q_sieve`] — the main entry point.
//! - [`SpecialQConfig`] — sieve parameters ``(A, B, q_min, q_max, threshold_scale)``.
//! - [`SpecialQResult`] — per-``q`` sieve output (relations + metadata).

use num_bigint::BigInt;
use num_traits::Zero;
use shared_numth::trial_smooth;

use super::{
    factor_base::FactorBase,
    norms::{algebraic_norm, norm_sign, norm_to_uint, rational_norm},
    Relation,
};
use crate::polyselect::PolyPair;

// ─── SpecialQConfig ───────────────────────────────────────────────────────────

/// Configuration for the special-q sieve.
///
/// Controls the sieve region ``|a| ≤ A``, ``1 ≤ b ≤ B``, the special-q range
/// ``[q_min, q_max]``, and the threshold scale factor used to decide which sieve
/// values are worth trial-dividing.
///
/// # Special-q range
///
/// The special primes are drawn from the algebraic factor base: only primes ``q`` in the
/// algebraic factor base with ``q_min ≤ q ≤ q_max`` are used. This ensures that ``q``
/// is already in the factor base and its index is known for the exponent vector.
///
/// # Threshold
///
/// Same as ``LineSieveConfig``: the threshold is ``threshold_scale × log2(B_alg)``.
/// A value of ``0.8`` (default) is conservative.
///
/// # Toy scale versus cryptographic scale
///
/// At toy scale, the threshold barely filters anything — the asymptotic win of the sieve
/// is under-exposed. See the module-level documentation.
#[derive(Debug, Clone)]
pub struct SpecialQConfig {
    /// Half-width of the sieve region: ``a`` ranges over ``−A..=A``.
    pub a_bound: u64,
    /// Height of the sieve region: ``b`` ranges over ``1..=B``.
    pub b_bound: u64,
    /// Minimum special prime (inclusive). Must be ≥ 2.
    pub q_min: u64,
    /// Maximum special prime (inclusive).
    pub q_max: u64,
    /// Scale factor for the smoothness threshold (default: 0.8).
    ///
    /// The threshold is ``threshold_scale × log2(B_alg)``. Lower values accept more
    /// candidates for trial division (slower but no missed relations); higher values
    /// filter more aggressively (faster but may miss relations with small norms).
    pub threshold_scale: f64,
}

impl SpecialQConfig {
    /// Construct a ``SpecialQConfig`` with the given bounds and default threshold scale (0.8).
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :param q_min: Minimum special prime (inclusive).
    /// :param q_max: Maximum special prime (inclusive).
    /// :returns: A new ``SpecialQConfig``.
    pub fn new(a_bound: u64, b_bound: u64, q_min: u64, q_max: u64) -> Self {
        Self { a_bound, b_bound, q_min, q_max, threshold_scale: 0.8 }
    }

    /// Construct a ``SpecialQConfig`` with an explicit threshold scale.
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :param q_min: Minimum special prime (inclusive).
    /// :param q_max: Maximum special prime (inclusive).
    /// :param threshold_scale: Scale factor for the smoothness threshold (0.0–1.0).
    /// :returns: A new ``SpecialQConfig``.
    pub fn with_threshold(
        a_bound: u64,
        b_bound: u64,
        q_min: u64,
        q_max: u64,
        threshold_scale: f64,
    ) -> Self {
        Self { a_bound, b_bound, q_min, q_max, threshold_scale }
    }
}

// ─── SpecialQResult ───────────────────────────────────────────────────────────

/// Per-``q`` output from the special-q sieve.
///
/// Collects the relations found for a single special prime ``q`` and one of its roots
/// ``r_q``. The ``q`` and ``r_q`` fields identify the sieve restriction; every relation
/// in ``relations`` carries ``q`` in its algebraic exponent vector.
#[derive(Debug, Clone)]
pub struct SpecialQResult {
    /// The special prime ``q``.
    pub q: u64,
    /// The root ``r_q ∈ [0, q)`` of ``f mod q`` used for this sieve run.
    pub r_q: u64,
    /// Relations collected in this sieve run.
    ///
    /// Every relation satisfies ``Relation::verify()`` and carries ``q`` in its
    /// algebraic exponent vector.
    pub relations: Vec<Relation>,
    /// The sieve area covered: ``(2·A + 1) × B`` pairs per ``(q, r_q)`` run, but only
    /// ``⌈(2·A + 1) / q⌉ × B`` pairs are actually trial-divided (the ``q``-restricted subset).
    ///
    /// Used for yield comparison: ``relations.len() / restricted_area`` is the per-area yield.
    pub restricted_area: u64,
}

// ─── special_q_sieve ─────────────────────────────────────────────────────────

/// Run the special-q sieve over the ``(a, b)`` rectangle and return per-``q`` results.
///
/// For each special prime ``q`` in ``[config.q_min, config.q_max]`` that is in the algebraic
/// factor base, and for each root ``r_q`` of ``f mod q``, runs a restricted line sieve over
/// ``a ≡ r_q·b (mod q)``, ``b ∈ 1..=B``. Returns a ``Vec<SpecialQResult>`` with one entry
/// per ``(q, r_q)`` pair.
///
/// # Correctness
///
/// Every returned ``Relation`` satisfies ``Relation::verify(poly, fb)``:
///
/// - ``gcd(a, b) = 1`` (coprimality).
/// - Both norms are fully smooth over their factor bases (``cofactor = 1``).
/// - The exponent vectors reconstruct the norm magnitudes.
/// - ``q`` appears in the algebraic exponent vector (guaranteed by the sieve restriction).
///
/// # Determinism
///
/// The output is deterministic for fixed ``(poly, fb, config)``: the sieve is a pure function
/// of its inputs, with no randomness.
///
/// # Toy scale versus cryptographic scale
///
/// At toy scale, the yield advantage of the special-q strategy over the plain line sieve is
/// under-exposed. See the module-level documentation for details.
///
/// :param poly: The NFS polynomial pair (provides ``f``, ``m`` for norm computation).
/// :param fb: The two-sided factor base.
/// :param config: Sieve region, special-q range, and threshold parameters.
/// :returns: Per-``(q, r_q)`` sieve results.
pub fn special_q_sieve(
    poly: &PolyPair,
    fb: &FactorBase,
    config: &SpecialQConfig,
) -> Vec<SpecialQResult> {
    let a_bound = config.a_bound as i64;

    // Precompute log2 values for all rational and algebraic primes.
    let rat_logs: Vec<f32> =
        fb.rational_primes.iter().map(|&p| (p as f32).log2()).collect();
    let alg_logs: Vec<f32> =
        fb.algebraic_ideals.iter().map(|ap| (ap.p as f32).log2()).collect();

    // Compute the smoothness threshold (same as line sieve).
    let b_alg_f = fb.b_alg.max(2) as f64;
    let threshold = ((config.threshold_scale * b_alg_f.log2()) as f32).max(1.0);

    // Factor-base slices for trial_smooth.
    let rat_fb_slice: Vec<u64> = fb.rational_primes.clone();
    let alg_fb_slice: Vec<u64> = fb.algebraic_ideals.iter().map(|ap| ap.p).collect();

    let mut all_results: Vec<SpecialQResult> = Vec::new();

    // Iterate over special primes q in the algebraic factor base within [q_min, q_max].
    //
    // We collect the distinct special primes (with their roots) from the algebraic factor
    // base. Each (p, r) ideal in the algebraic base with q_min ≤ p ≤ q_max is a candidate
    // special-q run. We group by prime p so that each (q, r_q) pair gets its own run.
    let special_ideals: Vec<(u64, u64)> = fb
        .algebraic_ideals
        .iter()
        .filter(|ap| ap.p >= config.q_min && ap.p <= config.q_max)
        .map(|ap| (ap.p, ap.r))
        .collect();

    for (q, r_q) in special_ideals {
        let result = sieve_for_q(poly, fb, config, q, r_q, a_bound, &rat_logs, &alg_logs,
                                  threshold, &rat_fb_slice, &alg_fb_slice);
        all_results.push(result);
    }

    all_results
}

/// Run the restricted line sieve for a single ``(q, r_q)`` pair.
///
/// For each ``b`` in ``1..=B``, only considers ``a`` values with ``a ≡ r_q·b (mod q)``,
/// stepping through ``a`` in increments of ``q``. Accumulates log contributions from the
/// rational and algebraic factor bases, then trial-divides survivors.
///
/// :param poly: The NFS polynomial pair.
/// :param fb: The two-sided factor base.
/// :param config: Sieve configuration.
/// :param q: The special prime.
/// :param r_q: The root of ``f mod q`` used for this run.
/// :param a_bound: ``config.a_bound`` as ``i64`` (pre-cast for convenience).
/// :param rat_logs: Precomputed ``log2(p)`` for rational primes.
/// :param alg_logs: Precomputed ``log2(p)`` for algebraic ideals.
/// :param threshold: Smoothness threshold.
/// :param rat_fb_slice: Rational factor base as ``Vec<u64>`` for ``trial_smooth``.
/// :param alg_fb_slice: Algebraic factor base primes for ``trial_smooth``.
/// :returns: ``SpecialQResult`` for this ``(q, r_q)`` pair.
#[allow(clippy::too_many_arguments)]
fn sieve_for_q(
    poly: &PolyPair,
    fb: &FactorBase,
    config: &SpecialQConfig,
    q: u64,
    r_q: u64,
    a_bound: i64,
    rat_logs: &[f32],
    alg_logs: &[f32],
    threshold: f32,
    rat_fb_slice: &[u64],
    alg_fb_slice: &[u64],
) -> SpecialQResult {
    let sieve_len = (2 * config.a_bound + 1) as usize;
    let mut relations: Vec<Relation> = Vec::new();
    let mut restricted_area: u64 = 0;

    for b in 1u64..=config.b_bound {
        let b_big = BigInt::from(b);

        // ── Step 1: Initialise sieve array ────────────────────────────────────
        let mut sieve: Vec<f32> = vec![0.0f32; sieve_len];

        // ── Step 2: Rational side sieve ───────────────────────────────────────
        // For prime p, mark a ≡ b·m (mod p).
        for (pi, &p) in fb.rational_primes.iter().enumerate() {
            let log_p = rat_logs[pi];
            let bm_mod_p = mod_u64_bigint(&(&b_big * &poly.m), p);
            let start_a = first_a_in_range(bm_mod_p, p, -a_bound);
            let mut a = start_a;
            while a <= a_bound {
                let idx = (a + a_bound) as usize;
                sieve[idx] += log_p;
                a += p as i64;
            }
        }

        // ── Step 3: Algebraic side sieve ──────────────────────────────────────
        // For ideal (p, r), mark a ≡ r·b (mod p).
        for (ai, ap) in fb.algebraic_ideals.iter().enumerate() {
            let p = ap.p;
            let r = ap.r;
            let log_p = alg_logs[ai];
            let rb_mod_p = (r as u128 * b as u128 % p as u128) as u64;
            let start_a = first_a_in_range(rb_mod_p, p, -a_bound);
            let mut a = start_a;
            while a <= a_bound {
                let idx = (a + a_bound) as usize;
                sieve[idx] += log_p;
                a += p as i64;
            }
        }

        // ── Step 4: Collect candidates in the q-restricted set ────────────────
        //
        // The special-q restriction: only consider a ≡ r_q·b (mod q).
        // Compute the starting a for this b.
        let rb_q = (r_q as u128 * b as u128 % q as u128) as u64;
        let start_a_q = first_a_in_range(rb_q, q, -a_bound);

        // Count the restricted candidates for this b.
        let mut a_q = start_a_q;
        while a_q <= a_bound {
            restricted_area += 1;
            a_q += q as i64;
        }

        // Sieve over the q-restricted a values.
        let mut a_q = start_a_q;
        while a_q <= a_bound {
            let idx = (a_q + a_bound) as usize;

            // Apply threshold filter.
            if sieve[idx] < threshold {
                a_q += q as i64;
                continue;
            }

            let a_big = BigInt::from(a_q);

            // ── Step 5a: Check gcd(a, b) = 1 ─────────────────────────────────
            if !is_coprime_i64_u64(a_q, b) {
                a_q += q as i64;
                continue;
            }

            // ── Step 5b: Compute norms ────────────────────────────────────────
            let rat_norm = rational_norm(&a_big, &b_big, &poly.m);
            let alg_norm = algebraic_norm(&a_big, &b_big, &poly.f);

            if rat_norm.is_zero() || alg_norm.is_zero() {
                a_q += q as i64;
                continue;
            }

            let rat_uint = match norm_to_uint(&rat_norm) {
                Ok(u) => u,
                Err(_) => { a_q += q as i64; continue; }
            };
            let alg_uint = match norm_to_uint(&alg_norm) {
                Ok(u) => u,
                Err(_) => { a_q += q as i64; continue; }
            };

            // ── Step 5c: Trial-smooth both norms ──────────────────────────────
            let rat_witness = trial_smooth(&rat_uint, rat_fb_slice);
            if !rat_witness.is_smooth() {
                a_q += q as i64;
                continue;
            }
            let alg_witness = trial_smooth(&alg_uint, alg_fb_slice);
            if !alg_witness.is_smooth() {
                a_q += q as i64;
                continue;
            }

            // ── Step 5d: Construct the Relation ───────────────────────────────
            let rational_sign = norm_sign(&rat_norm);
            if let Some(rel) = Relation::new(
                a_big,
                b_big.clone(),
                &rat_witness,
                &alg_witness,
                rational_sign,
                fb,
            ) {
                relations.push(rel);
            }

            a_q += q as i64;
        }
    }

    SpecialQResult { q, r_q, relations, restricted_area }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Compute ``n mod p`` as a non-negative ``u64`` in ``[0, p)``, for a ``BigInt`` ``n``.
fn mod_u64_bigint(n: &BigInt, p: u64) -> u64 {
    let p_big = BigInt::from(p);
    let r = n % &p_big;
    let r = if r < BigInt::zero() { r + &p_big } else { r };
    use num_traits::ToPrimitive;
    r.to_u64().expect("residue mod p fits in u64 for p ≤ u64::MAX")
}

/// Find the first ``a ≥ lo`` with ``a ≡ residue (mod p)``.
fn first_a_in_range(residue: u64, p: u64, lo: i64) -> i64 {
    let residue_i = residue as i64;
    let p_i = p as i64;
    let diff = lo - residue_i;
    let k = if diff >= 0 {
        (diff + p_i - 1) / p_i
    } else {
        -((-diff) / p_i)
    };
    residue_i + k * p_i
}

/// Check if ``gcd(|a|, b) = 1`` for ``a: i64``, ``b: u64``.
fn is_coprime_i64_u64(a: i64, b: u64) -> bool {
    let a_abs = a.unsigned_abs();
    gcd_u64(a_abs, b) == 1
}

/// Compute ``gcd(a, b)`` for ``u64`` values using the Euclidean algorithm.
fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_a_in_range_basic() {
        // residue=3, p=5, lo=-10: first a ≥ -10 with a ≡ 3 (mod 5) is -7.
        let a = first_a_in_range(3, 5, -10);
        assert_eq!(a, -7);
        assert_eq!(((a % 5) + 5) % 5, 3);
    }

    #[test]
    fn first_a_in_range_lo_positive() {
        // residue=2, p=7, lo=5: first a ≥ 5 with a ≡ 2 (mod 7) is 9.
        let a = first_a_in_range(2, 7, 5);
        assert_eq!(a, 9);
    }

    #[test]
    fn gcd_u64_cases() {
        assert_eq!(gcd_u64(7, 13), 1);
        assert_eq!(gcd_u64(12, 8), 4);
    }

    #[test]
    fn is_coprime_cases() {
        assert!(is_coprime_i64_u64(5, 7));
        assert!(is_coprime_i64_u64(-5, 7));
        assert!(!is_coprime_i64_u64(4, 2));
    }

    #[test]
    fn mod_u64_bigint_negative() {
        let n = BigInt::from(-7i64);
        assert_eq!(mod_u64_bigint(&n, 5), 3); // -7 mod 5 = 3 (canonical)
    }
}
