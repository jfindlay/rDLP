//! Line sieve over the ``(a, b)`` rectangle for GNFS relation collection.
//!
//! For each ``b`` in ``1..=B``, sieves the range ``|a| ≤ A`` on both the rational and algebraic
//! sides, accumulates approximate ``log p`` contributions, and trial-divides survivors to confirm
//! full smoothness. Each confirmed coprime pair ``(a, b)`` with both values smooth becomes a
//! ``Relation``.
//!
//! # Algorithm
//!
//! For each ``b``:
//!
//! 1. Initialise a sieve array of size ``2·A + 1`` (indices for ``a ∈ −A..=A``) with zeros.
//! 2. **Rational side**: for each prime ``p`` in the rational factor base, find the starting
//!    ``a ≡ b·m (mod p)`` and mark all ``a`` in range with ``+= log2(p)``.
//! 3. **Algebraic side**: for each ideal ``(p, r)`` in the algebraic factor base, find the
//!    starting ``a ≡ r·b (mod p)`` and mark all ``a`` in range with ``+= log2(p)``.
//! 4. **Threshold**: collect candidates where the sieve value exceeds a smoothness threshold.
//! 5. For each candidate ``(a, b)``:
//!    a. Skip if ``gcd(a, b) ≠ 1``.
//!    b. Compute the rational norm and the algebraic value ``F(a, b)``; skip if either is zero
//!       or overflows ``Uint<4>``.
//!    c. Call ``trial_smooth`` on both; skip if either has ``cofactor > 1``.
//!    d. Construct a ``Relation`` from the witnesses.
//!
//! # Toy scale versus cryptographic scale
//!
//! At toy scale (small ``A``, ``B``, ``B_rat``, ``B_alg``), the log-sieve barely beats direct
//! trial division of every ``(a, b)`` pair: the sieve array is tiny, the factor bases are small,
//! and the asymptotic win of the sieve (avoiding trial division for most pairs) is under-exposed.
//! The engineering heart of NFS — the ``log p`` mark-then-confirm pattern — is present in code,
//! but its yield advantage over brute-force trial division only becomes visible at cryptographic
//! scale (``B ≈ 10^7``, ``A ≈ 10^7``, ``B_rat/B_alg ≈ 10^6``). ``gnfs/docs/PEDAGOGY.md`` states
//! the same caveat.
//!
//! # Entry surface
//!
//! - [`line_sieve`] — the main entry point.
//! - [`LineSieveConfig`] — sieve parameters ``(A, B, threshold_scale)``.

use num_bigint::BigInt;
use num_traits::Zero;
use shared_numth::trial_smooth;

use super::{
    factor_base::FactorBase,
    norms::{algebraic_norm, norm_sign, norm_to_uint, rational_norm},
    Relation,
};

// ─── LineSieveConfig ─────────────────────────────────────────────────────────

/// Configuration for the line sieve.
///
/// Controls the sieve region ``|a| ≤ A``, ``1 ≤ b ≤ B``, and the threshold scale factor
/// used to decide which sieve values are worth trial-dividing.
///
/// # Threshold
///
/// The threshold is computed as ``threshold_scale × log2(B_alg)``, where ``B_alg`` is the
/// algebraic smoothness bound. The sieve accumulates ``log2(p)`` for each prime ``p`` dividing
/// the sieved value; for a fully smooth value ``V``, the total equals ``log2(|V|)``. The
/// threshold accepts candidates where the algebraic sieve has accumulated at least
/// ``threshold_scale × log2(B_alg)`` worth of contributions.
///
/// A value of ``0.8`` (default) is conservative: it accepts candidates where the algebraic
/// value ``F(a, b)`` has at least one prime factor ≥ ``B_alg^0.8``. The exact trial-division
/// step is the correctness gate; the threshold is only a pre-filter.
///
/// At toy scale, the threshold barely filters anything — the asymptotic win of the sieve is
/// under-exposed (see the module-level toy-scale caveat).
#[derive(Debug, Clone)]
pub struct LineSieveConfig {
    /// Half-width of the sieve region: ``a`` ranges over ``−A..=A``.
    pub a_bound: u64,
    /// Height of the sieve region: ``b`` ranges over ``1..=B``.
    pub b_bound: u64,
    /// Scale factor for the smoothness threshold (default: 0.8).
    ///
    /// The threshold is ``threshold_scale × log2(B_alg)``. Lower values accept more
    /// candidates for trial division (slower but no missed relations); higher values
    /// filter more aggressively (faster but may miss relations with small norms).
    pub threshold_scale: f64,
}

impl LineSieveConfig {
    /// Construct a ``LineSieveConfig`` with the given bounds and default threshold scale (0.8).
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :returns: A new ``LineSieveConfig``.
    pub fn new(a_bound: u64, b_bound: u64) -> Self {
        Self { a_bound, b_bound, threshold_scale: 0.8 }
    }

    /// Construct a ``LineSieveConfig`` with an explicit threshold scale.
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :param threshold_scale: Scale factor for the smoothness threshold (0.0–1.0).
    /// :returns: A new ``LineSieveConfig``.
    pub fn with_threshold(a_bound: u64, b_bound: u64, threshold_scale: f64) -> Self {
        Self { a_bound, b_bound, threshold_scale }
    }
}

// ─── line_sieve ───────────────────────────────────────────────────────────────

/// Run the line sieve over the ``(a, b)`` rectangle and return all confirmed relations.
///
/// For each ``b`` in ``1..=config.b_bound``, sieves the rational and algebraic sides over
/// ``a ∈ −config.a_bound..=config.a_bound``, then trial-divides survivors to confirm full
/// smoothness on both sides. Returns a ``Vec<Relation>`` of all confirmed relations.
///
/// # Correctness
///
/// Every returned ``Relation`` satisfies ``Relation::verify(poly, fb)``:
///
/// - ``gcd(a, b) = 1`` (coprimality).
/// - Both norms are fully smooth over their factor bases (``cofactor = 1``).
/// - The exponent vectors reconstruct the norm magnitudes.
///
/// # Determinism
///
/// The output is deterministic for fixed ``(poly, fb, config)``: the sieve is a pure function
/// of its inputs, with no randomness.
///
/// # Toy scale versus cryptographic scale
///
/// At toy scale, the log-sieve barely beats direct trial division — the asymptotic win is
/// under-exposed. See the module-level documentation for details.
///
/// :param poly: The NFS polynomial pair (provides ``f``, ``m`` for norm computation).
/// :param fb: The two-sided factor base.
/// :param config: Sieve region and threshold parameters.
/// :returns: All confirmed relations in the sieve region.
pub fn line_sieve(
    poly: &crate::polyselect::PolyPair,
    fb: &FactorBase,
    config: &LineSieveConfig,
) -> Vec<Relation> {
    let a_bound = config.a_bound as i64;
    let sieve_len = (2 * config.a_bound + 1) as usize; // indices 0..sieve_len → a = -A..=A

    // Precompute log2 values for all rational and algebraic primes.
    // These are stored as f32 for speed; precision is sufficient for threshold comparison.
    let rat_logs: Vec<f32> =
        fb.rational_primes.iter().map(|&p| (p as f32).log2()).collect();
    let alg_logs: Vec<f32> =
        fb.algebraic_ideals.iter().map(|ap| (ap.p as f32).log2()).collect();

    // Compute the smoothness threshold.
    //
    // The sieve accumulates log2(p) for each prime p dividing the sieved value. For a fully
    // smooth value V, the total accumulated log equals log2(|V|). The threshold is set to
    // threshold_scale × log2(B_alg), where B_alg is the algebraic smoothness bound.
    //
    // Rationale: the algebraic value is typically the bottleneck for smoothness (F(a, b) =
    // b^d · f(a/b) = c · N(a − b·α) is larger than the rational norm a − b·m for most (a, b)).
    // Using log2(B_alg) as the base means we accept candidates where the algebraic sieve has
    // accumulated at least threshold_scale × log2(B_alg) worth of contributions — i.e.,
    // candidates where F(a, b) has at least one factor ≥ B_alg^threshold_scale. The value
    // sieved is F (computed by `algebraic_norm`, tested by `algebraic_norm_cubic_hand_computed`
    // in norms.rs); it differs from the field norm N(a − b·α) by the factor c.
    //
    // This threshold is intentionally conservative (low) to avoid missing smooth pairs where
    // one norm is small (e.g., N_rat = 1 contributes nothing to the sieve). The exact
    // trial-division step is the correctness gate; the threshold is only a pre-filter.
    //
    // Scale note: at toy scale, the threshold barely filters anything — the asymptotic
    // win of the sieve (avoiding trial division for most pairs) is under-exposed. At
    // cryptographic scale (B_alg ≈ 10^7), the threshold filters ~99% of candidates, making
    // the sieve dramatically faster than brute-force trial division.
    let b_alg_f = fb.b_alg.max(2) as f64;
    let threshold = ((config.threshold_scale * b_alg_f.log2()) as f32).max(1.0);

    // Precompute the rational factor base as a slice for trial_smooth.
    let rat_fb_slice: Vec<u64> = fb.rational_primes.clone();
    let alg_fb_slice: Vec<u64> = fb.algebraic_ideals.iter().map(|ap| ap.p).collect();

    let mut relations: Vec<Relation> = Vec::new();

    for b in 1u64..=config.b_bound {
        let b_big = BigInt::from(b);

        // ── Step 1: Initialise sieve array ────────────────────────────────────
        let mut sieve: Vec<f32> = vec![0.0f32; sieve_len];

        // ── Step 2: Rational side sieve ───────────────────────────────────────
        // For prime p, mark a ≡ b·m (mod p).
        // The rational norm N_rat(a, b) = a − b·m is divisible by p iff a ≡ b·m (mod p).
        for (pi, &p) in fb.rational_primes.iter().enumerate() {
            let log_p = rat_logs[pi];
            // Compute b·m mod p (as a non-negative residue in [0, p)).
            let bm_mod_p = mod_u64_bigint(&(&b_big * &poly.m), p);
            // Find the first a in [-A, A] with a ≡ bm_mod_p (mod p).
            // a = bm_mod_p + k*p for integer k; we want the smallest a ≥ -A.
            let start_a = first_a_in_range(bm_mod_p, p, -a_bound);
            // Mark all a in [-A, A] with a ≡ bm_mod_p (mod p).
            let mut a = start_a;
            while a <= a_bound {
                let idx = (a + a_bound) as usize;
                sieve[idx] += log_p;
                a += p as i64;
            }
        }

        // ── Step 3: Algebraic side sieve ──────────────────────────────────────
        // For ideal (p, r), mark a ≡ r·b (mod p).
        // N_alg(a, b) is divisible by the ideal (p, r) iff a ≡ r·b (mod p).
        for (ai, ap) in fb.algebraic_ideals.iter().enumerate() {
            let p = ap.p;
            let r = ap.r;
            let log_p = alg_logs[ai];
            // Compute r·b mod p.
            let rb_mod_p = (r as u128 * b as u128 % p as u128) as u64;
            // Find the first a in [-A, A] with a ≡ rb_mod_p (mod p).
            let start_a = first_a_in_range(rb_mod_p, p, -a_bound);
            let mut a = start_a;
            while a <= a_bound {
                let idx = (a + a_bound) as usize;
                sieve[idx] += log_p;
                a += p as i64;
            }
        }

        // ── Step 4: Collect candidates above threshold ────────────────────────
        for idx in 0..sieve_len {
            if sieve[idx] < threshold {
                continue;
            }
            let a = idx as i64 - a_bound;
            let a_big = BigInt::from(a);

            // ── Step 5a: Check gcd(a, b) = 1 ─────────────────────────────────
            if !is_coprime_i64_u64(a, b) {
                continue;
            }

            // ── Step 5b: Compute norms ────────────────────────────────────────
            let rat_norm = rational_norm(&a_big, &b_big, &poly.m);
            let alg_norm = algebraic_norm(&a_big, &b_big, &poly.f);

            // Skip if either norm is zero (degenerate).
            if rat_norm.is_zero() || alg_norm.is_zero() {
                continue;
            }

            // Convert to Uint<4>; skip if overflow (shouldn't happen at toy scale).
            let rat_uint = match norm_to_uint(&rat_norm) {
                Ok(u) => u,
                Err(_) => continue,
            };
            let alg_uint = match norm_to_uint(&alg_norm) {
                Ok(u) => u,
                Err(_) => continue,
            };

            // ── Step 5c: Trial-smooth both norms ──────────────────────────────
            let rat_witness = trial_smooth(&rat_uint, &rat_fb_slice);
            if !rat_witness.is_smooth() {
                continue;
            }
            let alg_witness = trial_smooth(&alg_uint, &alg_fb_slice);
            if !alg_witness.is_smooth() {
                continue;
            }

            // ── Step 5d: Construct the Relation ───────────────────────────────
            let rational_sign = norm_sign(&rat_norm);
            if let Some(rel) =
                Relation::new(a_big, b_big.clone(), &rat_witness, &alg_witness, rational_sign, fb)
            {
                relations.push(rel);
            }
        }
    }

    relations
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Compute ``n mod p`` as a non-negative ``u64`` in ``[0, p)``, for a ``BigInt`` ``n``.
///
/// Handles negative ``n`` correctly by reducing to the canonical representative.
fn mod_u64_bigint(n: &BigInt, p: u64) -> u64 {
    let p_big = BigInt::from(p);
    let r = n % &p_big;
    let r = if r < BigInt::zero() { r + &p_big } else { r };
    // r is now in [0, p); convert to u64.
    use num_traits::ToPrimitive;
    r.to_u64().expect("residue mod p fits in u64 for p ≤ u64::MAX")
}

/// Find the first ``a ≥ lo`` with ``a ≡ residue (mod p)``.
///
/// Returns the smallest integer ``a ≥ lo`` such that ``a mod p == residue``.
/// ``residue`` must be in ``[0, p)``.
fn first_a_in_range(residue: u64, p: u64, lo: i64) -> i64 {
    // We want the smallest a ≥ lo with a ≡ residue (mod p).
    // a = residue + k*p for some integer k.
    // k = ceil((lo - residue) / p).
    let residue_i = residue as i64;
    let p_i = p as i64;
    // Compute k = ceil((lo - residue_i) / p_i).
    let diff = lo - residue_i;
    let k = if diff >= 0 {
        (diff + p_i - 1) / p_i // ceiling division for positive diff
    } else {
        -((-diff) / p_i) // floor division for negative diff gives ceiling
    };
    residue_i + k * p_i
}

/// Check if ``gcd(|a|, b) = 1`` for ``a: i64``, ``b: u64``.
///
/// Uses the Euclidean algorithm on ``u64`` values.
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
        // residue=3, p=5, lo=-10: first a ≥ -10 with a ≡ 3 (mod 5).
        // a = 3 + k*5; k = ceil((-10 - 3)/5) = ceil(-13/5) = ceil(-2.6) = -2.
        // a = 3 + (-2)*5 = 3 - 10 = -7. Check: -7 mod 5 = -7 + 10 = 3. ✓
        let a = first_a_in_range(3, 5, -10);
        assert_eq!(a, -7);
        assert_eq!(((a % 5) + 5) % 5, 3);
    }

    #[test]
    fn first_a_in_range_lo_positive() {
        // residue=2, p=7, lo=5: first a ≥ 5 with a ≡ 2 (mod 7).
        // a = 2 + k*7; k = ceil((5-2)/7) = ceil(3/7) = 1.
        // a = 2 + 7 = 9. Check: 9 mod 7 = 2. ✓
        let a = first_a_in_range(2, 7, 5);
        assert_eq!(a, 9);
    }

    #[test]
    fn first_a_in_range_lo_equals_residue() {
        // residue=3, p=5, lo=3: first a ≥ 3 with a ≡ 3 (mod 5) is a=3 itself.
        let a = first_a_in_range(3, 5, 3);
        assert_eq!(a, 3);
    }

    #[test]
    fn gcd_u64_coprime() {
        assert_eq!(gcd_u64(7, 13), 1);
        assert_eq!(gcd_u64(1, 100), 1);
    }

    #[test]
    fn gcd_u64_non_coprime() {
        assert_eq!(gcd_u64(12, 8), 4);
        assert_eq!(gcd_u64(100, 25), 25);
    }

    #[test]
    fn is_coprime_i64_u64_cases() {
        assert!(is_coprime_i64_u64(5, 7));
        assert!(is_coprime_i64_u64(-5, 7));
        assert!(!is_coprime_i64_u64(4, 2));
        assert!(!is_coprime_i64_u64(-6, 3));
    }

    #[test]
    fn mod_u64_bigint_positive() {
        let n = BigInt::from(17i64);
        assert_eq!(mod_u64_bigint(&n, 5), 2); // 17 mod 5 = 2
    }

    #[test]
    fn mod_u64_bigint_negative() {
        let n = BigInt::from(-7i64);
        assert_eq!(mod_u64_bigint(&n, 5), 3); // -7 mod 5 = 3 (canonical)
    }
}
