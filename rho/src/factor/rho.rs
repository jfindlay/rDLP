//! Pollard rho factorization: Floyd, Brent, and batched-GCD variants.
//!
//! All variants share the pseudorandom map `f(x) = (x² + c) mod N` and detect
//! a non-trivial factor via `gcd(|xᵢ − xⱼ|, N)`.  Arithmetic is in `u128`, which
//! covers the supported range below 2^80 comfortably.  For inputs where
//! `N < 2^64`, products `x*x` fit directly in `u128`.  For 64–80-bit N the
//! `mulmod` helper uses `u128::carrying_mul` to form the 256-bit product and
//! reduces it mod N in two u128 steps.  N below `2^80` is the supported range
//! ([`MAX_SUPPORTED`]); [`factor`] reports larger N as [`FactorError::Unsupported`].
//!
//! ## Variants
//!
//! - [`floyd`] — baseline Floyd cycle detection.  Two pointers advance at speed
//!   1 and 2; GCD tested at every step.
//! - [`brent`] — Brent's improved cycle detection, ~24% fewer `f` evaluations
//!   than Floyd on average.
//! - [`brent_batched`] — Brent + Montgomery's batched-GCD trick: accumulate a
//!   product of `|xᵢ − y|` values over a batch of `batch_size` steps, take one
//!   GCD per batch.  Falls back to step-by-step on `gcd = N`.
//! - [`factor`] — public entry point: runs `brent_batched` in parallel over
//!   multiple `c` values using `rayon`, returning the first non-trivial factor
//!   found as a [`Factorization`], or a [`FactorError`] that keeps invalid input, primes,
//!   unsupported sizes, spent budgets, and a violated divisor invariant apart.
//!
//! ## Step bounds
//!
//! The three primitives take an explicit `max_steps`, counted in evaluations of `f`.  Their
//! `Option` result has exactly one failure meaning: *this attempt failed* (a degenerate
//! collision with `gcd = n`, or the bound was spent).  It says nothing about whether `n` is
//! prime; [`factor`] owns that classification.

use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use shared_numth::{is_prime, MILLER_RABIN_PROOF_BOUND};

use crypto_bigint::Uint;

/// Exclusive upper bound of the supported domain: `mulmod` is exact only for `n < 2^80`.
pub const MAX_SUPPORTED: u128 = 1 << 80;

// `shared_numth::is_prime` tests the first thirteen prime bases (2 through 41) at and above
// ψ₁₂, so its "prime" verdict is a proof below ψ₁₃ = `MILLER_RABIN_PROOF_BOUND` ≈ 2^81.5
// (Sorenson & Webster 2017, Theorem 1.1).  Every supported odd n is below that bound, which this
// assertion keeps true: a "prime" verdict from `is_prime` is a primality proof on the whole
// supported range.
const _: () = assert!(MAX_SUPPORTED <= MILLER_RABIN_PROOF_BOUND);

/// A verified non-trivial factorization: `1 < d < n` and `n % d == 0`.
///
/// `#[non_exhaustive]` keeps the fields readable everywhere while only this crate can construct
/// one, so an unchecked divisor cannot be built by a caller.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Factorization {
    /// A proper divisor of the factored number.
    pub d: u128,
    /// Which path produced `d`.
    pub evidence: FactorEvidence,
}

/// Provenance of a factor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FactorEvidence {
    /// The number is even and larger than 2, so `2` divides it; no rho walk ran.
    Even,
    /// Pollard rho with `f(x) = x² + c` found the divisor via `gcd(|xᵢ − xⱼ|, n)`.
    Rho {
        /// The polynomial constant of the successful walk.
        c: u64,
    },
}

/// Why [`factor`] did not return a factor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FactorError {
    /// `n ≤ 1` has no non-trivial factorization to ask for.
    InvalidInput,
    /// `n` is prime, certified by a primality test that is proven on `n`.
    Prime,
    /// `n` is outside the implementation's domain: at least [`MAX_SUPPORTED`].
    Unsupported,
    /// Every constant tried spent its step bound or hit a degenerate collision.  This is a
    /// resource outcome and proves nothing about `n`.
    BudgetExhausted {
        /// Number of constants `c` whose walks were run and failed.
        constants: u64,
    },
    /// A walk returned a candidate `d` that is neither `1` nor `n` yet does not divide `n`.
    /// A gcd with `n` always divides `n`, so this is a breach of the divisor invariant in the
    /// arithmetic, not a failed attempt; it is reported at once and never retried.
    Inconsistent {
        /// The candidate that failed `n mod d == 0`.
        d: u128,
    },
}

/// Step bound for one rho attempt on `n`: `1024 + 64·(⌊n^{1/4}⌋ + 1)` evaluations of `f`.
///
/// The smallest prime factor `p ≤ √n` is found in about `√p` walk steps, so Brent's roughly
/// `3·√p` evaluations stay well inside `64·n^{1/4}`; the constant term covers tiny `n`.
fn attempt_steps(n: u128) -> u64 {
    let quartic_root = n.isqrt().isqrt() + 1;
    // n < 2^80 ⟹ quartic_root ≤ 2^20 + 1, so the sum is far below u64::MAX.
    1024 + 64 * quartic_root as u64
}

/// Judge one walk's result: a proper divisor of `n`, a failed attempt, or an invariant breach.
///
/// `None` and the trivial candidates `1` and `n` are failed attempts (`Ok(None)`).  Any other
/// candidate is checked as `1 < d < n` and `n mod d == 0`, independently of how it was found;
/// one that fails is [`FactorError::Inconsistent`], never a failed attempt.
fn judge_candidate(n: u128, d: Option<u128>) -> Result<Option<u128>, FactorError> {
    match d {
        None => Ok(None),
        Some(d) if d == 1 || d == n => Ok(None),
        Some(d) if d > 1 && d < n && n % d == 0 => Ok(Some(d)),
        Some(d) => Err(FactorError::Inconsistent { d }),
    }
}

/// Run `walk(c)` over the constants `c = 1..=max_c` in parallel and return the first verified
/// divisor.
///
/// Constants with `c ≡ 0` or `c ≡ n − 2 (mod n)` are skipped.  A walk that returns no candidate
/// or a trivial one counts as a failed attempt; an inconsistent candidate ends the search with
/// [`FactorError::Inconsistent`].  `n` must be odd and composite.
fn search_constants<W>(n: u128, max_c: u64, walk: W) -> Result<Factorization, FactorError>
where
    W: Fn(u64) -> Option<u128> + Sync,
{
    let tried = AtomicU64::new(0);
    let found = (1u64..=max_c)
        .into_par_iter()
        .filter(|&c| {
            let c_mod = c as u128 % n;
            c_mod != 0 && c_mod != n - 2
        })
        .find_map_any(|c| match judge_candidate(n, walk(c)) {
            Ok(Some(d)) => Some(Ok(Factorization { d, evidence: FactorEvidence::Rho { c } })),
            Ok(None) => {
                tried.fetch_add(1, Ordering::Relaxed);
                None
            }
            Err(e) => Some(Err(e)),
        });
    found.unwrap_or_else(|| {
        Err(FactorError::BudgetExhausted { constants: tried.load(Ordering::Relaxed) })
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Compute `2^exp mod m` via iterative doubling-and-reduction.
///
/// At each iteration the accumulator stays below `m`, so doubling it gives a
/// value below `2m`.  For `m < 2^127` this fits in a `u128` and the reduction
/// is a single conditional subtract.  This avoids the need for any widening
/// multiplication during the exponentiation itself.
#[inline]
fn pow2_mod(exp: u32, m: u128) -> u128 {
    let mut r = 1u128;
    for _ in 0..exp {
        // r < m < 2^127; doubling gives r < 2^128, safe in u128.
        r <<= 1;
        if r >= m {
            r -= m;
        }
    }
    r
}

/// Modular multiplication `a * b mod m` for `a, b < m < 2^80`.
///
/// Uses `u128::carrying_mul` to form the exact 256-bit product `(lo, hi)`, then
/// reduces in two steps:
///
/// 1. If `hi == 0` the 128-bit `lo` is reduced directly.
/// 2. Otherwise: `a*b = hi * 2^128 + lo ≡ hi * r2_128 + lo (mod m)`, where
///    `r2_128 = 2^128 mod m`.  For `m < 2^80` the `hi` word satisfies
///    `hi < m^2 / 2^128 < 2^(160-128) = 2^32`, so `hi * r2_128 < 2^112` and
///    the final multiplication fits in a `u128`.
///
/// # Panics (debug)
///
/// Panics in debug builds if `a >= m || b >= m || m >= 2^127`.
#[inline]
fn mulmod(a: u128, b: u128, m: u128) -> u128 {
    debug_assert!(a < m && b < m && m < (1u128 << 127), "mulmod: inputs out of range");
    let (lo, hi) = a.carrying_mul(b, 0);
    if hi == 0 {
        return lo % m;
    }
    // r2_128 = 2^128 mod m.  For m < 2^80, r2_128 < m < 2^80.
    let r2_128 = pow2_mod(128, m);
    // hi < 2^32 (for m < 2^80), r2_128 < 2^80: product < 2^112, safe.
    (hi * r2_128 % m + lo % m) % m
}

/// Step function: `f(x) = (x² + c) mod n`.
#[inline]
fn f(x: u128, c: u128, n: u128) -> u128 {
    (mulmod(x, x, n) + c) % n
}

/// Euclidean GCD.
#[inline]
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        a %= b;
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// Absolute difference: `|a − b|`.
#[inline]
fn abs_diff(a: u128, b: u128) -> u128 {
    if a >= b { a - b } else { b - a }
}

// ---------------------------------------------------------------------------
// Floyd's cycle detection
// ---------------------------------------------------------------------------

/// Pollard rho with Floyd's cycle detection (baseline).
///
/// Two pointers — tortoise (speed 1) and hare (speed 2) — advance together;
/// at each step we test `gcd(|tortoise − hare|, n)`.  Terminates when a
/// non-trivial factor is found, the pointers collide with `gcd = n` (rare
/// degenerate failure), or `max_steps` evaluations of `f` are spent.
///
/// # Arguments
///
/// * `n` — composite number to factor (must be > 1, not prime, `n < 2^80`).
/// * `c` — polynomial constant; avoid 0 and `n − 2` (degenerate cycles).
/// * `x0` — starting point (typically 2).
/// * `max_steps` — bound on evaluations of `f` (each iteration costs three).
///
/// # Returns
///
/// A non-trivial factor of `n`, or `None` when this attempt failed: a
/// degenerate collision (retry with a different `c`) or the step bound was
/// spent.  `None` does not show that `n` is prime.
pub fn floyd(n: u128, c: u128, x0: u128, max_steps: u64) -> Option<u128> {
    let mut tortoise = x0;
    let mut hare = x0;
    let mut steps = 0u64;
    while steps < max_steps {
        steps += 3;
        tortoise = f(tortoise, c, n);
        hare = f(f(hare, c, n), c, n);
        let d = gcd(abs_diff(tortoise, hare), n);
        if d == n {
            return None;
        }
        if d != 1 {
            return Some(d);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Brent's cycle detection
// ---------------------------------------------------------------------------

/// Pollard rho with Brent's cycle detection (~24% fewer f-evaluations vs Floyd).
///
/// Brent's algorithm advances the hare freely for a power-of-2 window `r`, then
/// tests GCD against the tortoise (fixed at the start of the window).  This
/// eliminates the redundant second-pointer evaluation Floyd requires, and the
/// tortoise advances only O(log λ) times total where λ is the cycle length.
///
/// # Arguments
///
/// * `n`, `c`, `x0` — same as [`floyd`].
/// * `max_steps` — bound on evaluations of `f`, skip-ahead and comparison phases together.
///
/// # Returns
///
/// A non-trivial factor of `n`, or `None` when this attempt failed (degenerate collision
/// or the step bound was spent); `None` does not show that `n` is prime.
pub fn brent(n: u128, c: u128, x0: u128, max_steps: u64) -> Option<u128> {
    let mut y = x0;
    let mut r: u64 = 1;
    let mut steps = 0u64;
    loop {
        let x = y; // tortoise: fixed for this r-step window
        // Advance hare r steps (skip-ahead phase — not compared to tortoise yet).
        for _ in 0..r {
            if steps >= max_steps {
                return None;
            }
            steps += 1;
            y = f(y, c, n);
        }
        // Comparison phase: advance hare another r steps, testing against tortoise.
        let mut ys = y;
        for _ in 0..r {
            if steps >= max_steps {
                return None;
            }
            steps += 1;
            ys = f(ys, c, n);
            let d = gcd(abs_diff(x, ys), n);
            if d == n {
                return None;
            }
            if d != 1 {
                return Some(d);
            }
        }
        // The window doubles; overflow is unreachable below `max_steps` but fails the attempt.
        r = r.checked_mul(2)?;
    }
}

// ---------------------------------------------------------------------------
// Brent + Montgomery batched GCD
// ---------------------------------------------------------------------------

/// Pollard rho with Brent's cycle detection and Montgomery's batched-GCD trick.
///
/// Instead of testing `gcd(|x − y|, n)` at every step, accumulate the product
/// `Π |x − yᵢ|` over `batch_size` hare steps and take a single GCD per batch.
/// One GCD + `batch_size` multiplications dominates `batch_size` GCDs by a
/// factor of ~4–8× when `batch_size` is 64–128.
///
/// On batch failure (`gcd = n`), falls back to step-by-step [`brent_fallback`]
/// from the pre-batch snapshot to pinpoint the actual factor.
///
/// # Arguments
///
/// * `n`, `c`, `x0` — same as [`floyd`].
/// * `batch_size` — number of steps accumulated before each GCD (typical: 128); `0` is
///   treated as `1`.
/// * `max_steps` — bound on evaluations of `f` in the main walk.  The step-by-step re-walk of
///   one batch after a `gcd = n` batch is not charged, so the attempt can exceed the bound by at
///   most `batch_size` evaluations.
///
/// # Returns
///
/// A non-trivial factor of `n`, or `None` when this attempt failed (degenerate collision
/// or the step bound was spent); `None` does not show that `n` is prime.
pub fn brent_batched(
    n: u128,
    c: u128,
    x0: u128,
    batch_size: usize,
    max_steps: u64,
) -> Option<u128> {
    let batch_size = batch_size.max(1);
    let mut y = x0;
    let mut r: u64 = 1;
    let mut steps = 0u64;
    loop {
        let x = y; // tortoise fixed for this window
        for _ in 0..r {
            if steps >= max_steps {
                return None;
            }
            steps += 1;
            y = f(y, c, n);
        }
        let mut k: u64 = 0;
        let mut ys = y;
        while k < r {
            let snapshot = ys; // hare position before this batch
            let batch = batch_size.min((r - k) as usize);
            let mut product = 1u128;
            for _ in 0..batch {
                if steps >= max_steps {
                    return None;
                }
                steps += 1;
                ys = f(ys, c, n);
                let delta = abs_diff(x, ys);
                if delta == 0 {
                    // Zero delta ⟹ accumulated product stays zero; fall back to pinpoint.
                    return brent_fallback(n, c, snapshot, x, batch);
                }
                product = mulmod(product, delta, n);
            }
            let d = gcd(product, n);
            if d == n {
                return brent_fallback(n, c, snapshot, x, batch);
            }
            if d != 1 {
                return Some(d);
            }
            k += batch as u64;
        }
        // The window doubles; overflow is unreachable below `max_steps` but fails the attempt.
        r = r.checked_mul(2)?;
    }
}

/// Step-by-step fallback used by [`brent_batched`] when a batch GCD equals `n`.
///
/// Resumes from the hare snapshot `ys` one step at a time, testing GCD against
/// the tortoise `x`, to pinpoint the factor that the batch overshot.  The batch
/// contained a step whose `gcd(|x − yᵢ|, n)` is not 1, so at most `limit` (the
/// batch length) steps are needed; `None` is a degenerate collision.
fn brent_fallback(n: u128, c: u128, mut ys: u128, x: u128, limit: usize) -> Option<u128> {
    for _ in 0..limit {
        ys = f(ys, c, n);
        let d = gcd(abs_diff(x, ys), n);
        if d == n {
            return None;
        }
        if d != 1 {
            return Some(d);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Multi-c parallel entry point
// ---------------------------------------------------------------------------

/// Factor `n` using parallel Pollard rho (Brent + batched GCD) over multiple `c` values.
///
/// Classifies `n` first, in this order: `n ≤ 1` is invalid; `n ≥ 2^80` is outside the `mulmod`
/// domain; `2` is prime; an even `n > 2` yields `d = 2` with [`FactorEvidence::Even`] and no
/// walk; an odd `n` is tested for primality.  Composite odd `n` is then attacked with
/// rayon tasks over `c ∈ 1..=max_c`, each running [`brent_batched`] from `x0 = 2` under the
/// step bound `1024 + 64·(⌊n^{1/4}⌋ + 1)` evaluations of `f`.  Constants with `c ≡ 0` or
/// `c ≡ n − 2 (mod n)` are skipped: c = 0 gives a fixed point at 0; c = n−2 gives a trivial
/// cycle via the identity `x² − 2 = (x−1)² − 1`.  Perfect squares and prime powers need no
/// special case: the walk finds the prime `p` whenever the sequence collides mod `p` before it
/// does mod `n`, and the divisor is checked as `1 < d < n` with `n % d == 0` before it is
/// returned.
///
/// # Arguments
///
/// * `n` — number to factor.
/// * `max_c` — number of constants `c = 1, 2, …, max_c` to try (typical: logical CPU count).
/// * `batch_size` — batch size for GCD accumulation (typical: 128).
///
/// # Errors
///
/// * [`FactorError::InvalidInput`] — `n ≤ 1`.
/// * [`FactorError::Prime`] — `n` is prime (`n = 2`, or odd and accepted by `is_prime`, which
///   is a proof on the whole supported range because [`MAX_SUPPORTED`] does not exceed
///   `shared_numth::MILLER_RABIN_PROOF_BOUND`).
/// * [`FactorError::Unsupported`] — `n` is at least [`MAX_SUPPORTED`].
/// * [`FactorError::BudgetExhausted`] — every constant failed; `constants` counts the walks
///   that ran.  This proves nothing about `n`; `max_c = 0` runs none.  A walk whose gcd is `1`
///   or `n` is one such failed walk.
/// * [`FactorError::Inconsistent`] — a walk returned a candidate other than `1` or `n` that
///   does not divide `n`; reported at once, without trying further constants.
pub fn factor(n: u128, max_c: u64, batch_size: usize) -> Result<Factorization, FactorError> {
    if n <= 1 {
        return Err(FactorError::InvalidInput);
    }
    if n >= MAX_SUPPORTED {
        return Err(FactorError::Unsupported);
    }
    if n == 2 {
        return Err(FactorError::Prime);
    }
    // Fast path: an even n > 2 is divisible by 2, which is proper because n > 2.
    if n % 2 == 0 {
        return Ok(Factorization { d: 2, evidence: FactorEvidence::Even });
    }
    if is_prime(&Uint::<4>::from(n)) {
        // n < MAX_SUPPORTED ≤ MILLER_RABIN_PROOF_BOUND, so this verdict is a proof.
        return Err(FactorError::Prime);
    }

    let max_steps = attempt_steps(n);
    search_constants(n, max_c, |c| brent_batched(n, c as u128, 2, batch_size, max_steps))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Trivial candidates and missing candidates are failed attempts; proper divisors pass.
    #[test]
    fn judge_candidate_classifies_each_case() {
        let n = 143u128;
        assert_eq!(judge_candidate(n, None), Ok(None));
        assert_eq!(judge_candidate(n, Some(1)), Ok(None));
        assert_eq!(judge_candidate(n, Some(n)), Ok(None));
        assert_eq!(judge_candidate(n, Some(11)), Ok(Some(11)));
        assert_eq!(judge_candidate(n, Some(13)), Ok(Some(13)));
    }

    /// Any other candidate is an invariant breach carrying the candidate, including `0` and
    /// values above `n`.
    #[test]
    fn judge_candidate_flags_non_dividing_candidates() {
        let n = 143u128;
        for d in [0u128, 2, 7, 12, 142, 144, u128::MAX] {
            assert_eq!(judge_candidate(n, Some(d)), Err(FactorError::Inconsistent { d }), "d={d}");
        }
    }

    /// A non-dividing candidate from any walk ends the search at once with that candidate.
    #[test]
    fn search_surfaces_inconsistent_candidate_without_retrying() {
        let n = 143u128;
        assert_eq!(search_constants(n, 64, |_| Some(7)), Err(FactorError::Inconsistent { d: 7 }));
    }

    /// An inconsistent candidate from one constant wins over failed attempts from the others.
    #[test]
    fn search_inconsistent_beats_failed_attempts() {
        let n = 143u128;
        let result = search_constants(n, 8, |c| if c == 3 { Some(5) } else { None });
        assert_eq!(result, Err(FactorError::Inconsistent { d: 5 }));
    }

    /// A candidate equal to `n` (gcd = n) consumes an attempt like `None` does and does not stop
    /// the search.
    #[test]
    fn search_trivial_candidates_consume_attempts() {
        let n = 143u128;
        assert_eq!(
            search_constants(n, 5, |_| Some(n)),
            Err(FactorError::BudgetExhausted { constants: 5 })
        );
        assert_eq!(
            search_constants(n, 5, |c| if c % 2 == 0 { Some(1) } else { Some(n) }),
            Err(FactorError::BudgetExhausted { constants: 5 })
        );
        assert_eq!(
            search_constants(n, 5, |_| None),
            Err(FactorError::BudgetExhausted { constants: 5 })
        );
    }

    /// A proper divisor from a walk is returned with the constant that produced it.
    #[test]
    fn search_returns_verified_divisor() {
        let n = 143u128;
        let found = search_constants(n, 8, |c| if c == 4 { Some(11) } else { None })
            .expect("a proper divisor is found");
        assert_eq!(found.d, 11);
        assert_eq!(found.evidence, FactorEvidence::Rho { c: 4 });
    }

    /// Skipped constants `c ≡ 0` and `c ≡ n − 2` never run a walk and are not counted as
    /// attempts.
    #[test]
    fn search_skips_degenerate_constants() {
        // n = 9: c = 7 ≡ n − 2 and c = 9 ≡ 0 are skipped, leaving c ∈ {1..6, 8}.
        assert_eq!(
            search_constants(9, 9, |_| None),
            Err(FactorError::BudgetExhausted { constants: 7 })
        );
    }
}
