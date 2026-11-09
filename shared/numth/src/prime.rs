//! Miller–Rabin primality testing, proven for ``n < 3,317,044,064,679,887,385,961,981``.
//!
//! The core primitive is [`miller_rabin`], which runs the Miller–Rabin strong
//! pseudoprime test for a caller-supplied list of witness bases.  The public
//! entry point [`is_prime`] selects the minimal deterministic witness set for
//! the magnitude of ``n`` (following the published tables from Pomerance,
//! Selfridge, and Wagstaff, and later refinements by Jaeschke and Sorenson &
//! Webster).
//!
//! # Determinism
//!
//! The first twelve prime bases ``{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}`` are proven
//! below ``318,665,857,834,031,151,167,461`` (ψ₁₂, ≈ ``2^{78.1}``): that number is the least
//! composite passing all twelve.  Adding base 41 extends the proof to
//! ``3,317,044,064,679,887,385,961,981`` ([`MILLER_RABIN_PROOF_BOUND`], ψ₁₃, ≈ ``2^{81.5}``).
//! Above ψ₁₃ (which still fits in ``Uint<4>``) the same thirteen bases give only a
//! *probable-prime* verdict with no proven error bound; a fixed set of bases has composite
//! strong pseudoprimes, and ``is_prime`` can return ``true`` for one.  A ``false`` verdict is
//! always a proof of compositeness.
//!
//! # Modular arithmetic
//!
//! Modular exponentiation is performed via ``FpNaive::<4>::pow``, which uses
//! schoolbook square-and-multiply.  The modulus ``n`` is treated as an opaque
//! integer; no primality of the modulus is assumed.

use crypto_bigint::Uint;
use shared_field::{Fp, FpNaive};

// ── Deterministic witness-set bounds ──────────────────────────────────────────
//
// Each entry is (exclusive_upper_bound_as_u128, witnesses).  The bounds are
// taken from:
//   - Pomerance, Selfridge, Wagstaff (1980)
//   - Jaeschke (1993)
//   - Sorenson & Webster (2015) — the 3,474,749,660,383, 318,665,857,834,031,151,167,461
//     and 3,317,044,064,679,887,385,961,981 entries (constants recalled, not re-read; the
//     factorization of the twelve-base bound is checked in the tests)
//   - Wikipedia "Miller–Rabin primality test" (deterministic variants table)
//
// The last two bounds exceed u64::MAX, so bounds are stored as u128.  The thirteenth
// base (41) is what makes the final row a proof: the twelve-base bound is itself a
// strong pseudoprime to bases 2..37.

/// Exclusive bound below which the first twelve prime bases (2 through 37) prove primality.
///
/// ψ₁₂ ≈ 2^78.08 is the least composite that is a strong pseudoprime to all twelve bases.
const TWELVE_BASE_PROOF_BOUND: u128 = 318_665_857_834_031_151_167_461;

/// Exclusive bound below which a [`is_prime`] "prime" verdict is a proof.
///
/// ψ₁₃ ≈ 2^81.46 is the least composite that is a strong pseudoprime to the first thirteen
/// prime bases (2 through 41).  From this value up, a "prime" verdict is only "probable
/// prime", with no proven error bound.  A "composite" verdict is always a proof.
pub const MILLER_RABIN_PROOF_BOUND: u128 = 3_317_044_064_679_887_385_961_981;

/// Witness sets for deterministic Miller–Rabin, ordered by increasing bound.
///
/// Each entry is ``(exclusive_upper_bound, witnesses)``.
static DETERMINISTIC_WITNESSES: &[(u128, &[u64])] = &[
    (2_047, &[2]),
    (1_373_653, &[2, 3]),
    (3_215_031_751, &[2, 3, 5, 7]),
    (3_474_749_660_383, &[2, 3, 5, 7, 11, 13]),
    (341_550_071_728_321, &[2, 3, 5, 7, 11, 13, 17]),
    (3_825_123_056_546_413_051, &[2, 3, 5, 7, 11, 13, 17, 19, 23]),
    (TWELVE_BASE_PROOF_BOUND, &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]),
    (MILLER_RABIN_PROOF_BOUND, &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]),
];

/// The thirteen-base set, used for ``n`` at or above every proof bound.
static FULL_WITNESSES: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];

// ── Core algorithm ─────────────────────────────────────────────────────────────

/// Test whether ``n`` is probably prime using Miller–Rabin with the given witness bases.
///
/// Returns ``true`` if all witnesses say "probably prime", ``false`` if any
/// witness says "composite".  ``false`` is always a proof of compositeness;
/// ``true`` proves primality only when the witnesses form a published
/// deterministic set for ``n``'s magnitude (see [`is_prime`] and
/// [`MILLER_RABIN_PROOF_BOUND`]).
///
/// # Preconditions
///
/// - ``n`` must be odd and ``n >= 5``.  The caller is responsible for handling
///   ``n <= 4`` as special cases before calling this function.
/// - Each witness ``a`` must satisfy ``1 < a < n - 1``; witnesses that violate
///   this are skipped (they give trivially uninformative results).
///
/// # Algorithm
///
/// Write ``n - 1 = 2^s * d`` with ``d`` odd.  For each witness ``a``:
///
/// 1. Compute ``x = a^d mod n``.
/// 2. If ``x == 1`` or ``x == n - 1``, the witness passes (continue to next).
/// 3. Square ``x`` up to ``s - 1`` times.  If ``x == n - 1`` at any point,
///    the witness passes.
/// 4. If no squaring produced ``n - 1``, ``n`` is composite.
pub fn miller_rabin(n: &Uint<4>, witnesses: &[u64]) -> bool {
    // Factor out powers of 2 from n - 1.
    let n_minus_1 = n.wrapping_sub(&Uint::<4>::ONE);
    let s = n_minus_1.trailing_zeros();
    // d = (n - 1) / 2^s
    let d = n_minus_1 >> s;

    for &a in witnesses {
        // Skip witnesses that are >= n (uninformative for small n).
        let a_uint = Uint::<4>::from(a);
        if a_uint >= *n {
            continue;
        }
        // Also skip a == 0 or a == 1 (trivially uninformative).
        if a <= 1 {
            continue;
        }

        // x = a^d mod n via FpNaive schoolbook exponentiation.
        let fp_a = FpNaive::<4>::from_uint(a_uint, n);
        let mut x = fp_a.pow(&d, n);

        let one = FpNaive::<4>::one(n);
        let n_minus_1_fp = FpNaive::<4>::from_uint(n_minus_1, n);

        if x == one || x == n_minus_1_fp {
            continue; // This witness passes.
        }

        // Square up to s - 1 times looking for n - 1.
        let mut composite = true;
        for _ in 0..(s - 1) {
            x = x.square(n);
            if x == n_minus_1_fp {
                composite = false;
                break;
            }
        }

        if composite {
            return false; // Definite composite.
        }
    }

    true // All witnesses passed: probably prime.
}

/// Primality test for ``n`` using the published minimal witness sets.
///
/// Uses the deterministic witness set that is sufficient for all
/// ``n < 3,215,031,751`` (4 witnesses), scaling up to twelve bases below
/// ``318,665,857,834,031,151,167,461`` and thirteen (bases 2 through 41) below
/// [`MILLER_RABIN_PROOF_BOUND`] (``3,317,044,064,679,887,385,961,981`` ≈ ``2^{81.5}``).
/// Below that bound both verdicts are proofs.  At or above it the thirteen bases give a
/// "probable prime" verdict with no proven error bound: ``false`` is still a proof of
/// compositeness, but ``true`` is not a proof of primality.
///
/// # Special cases
///
/// - ``n <= 1``: returns ``false`` (not prime).
/// - ``n == 2`` or ``n == 3``: returns ``true``.
/// - Even ``n > 2``: returns ``false``.
/// - ``n < 9``: handled directly (3, 5, 7 are prime; 4, 6, 8 are not).
pub fn is_prime(n: &Uint<4>) -> bool {
    // Handle small and degenerate cases.
    if *n <= Uint::<4>::ONE {
        return false;
    }
    if *n == Uint::<4>::from(2u64) || *n == Uint::<4>::from(3u64) {
        return true;
    }
    // Even numbers > 2 are composite.  bit(0) returns CtChoice; convert to bool.
    // Low bit == 0 means even.
    if !bool::from(n.bit(0)) {
        return false;
    }
    // n == 5 or n == 7 are prime; n == 9 = 3*3 is composite.
    if *n == Uint::<4>::from(5u64) || *n == Uint::<4>::from(7u64) {
        return true;
    }
    if *n < Uint::<4>::from(9u64) {
        // Only remaining odd values < 9 are 3, 5, 7 (handled above) and 9.
        // Since we're here, n must be 9.
        return false;
    }

    // Select the minimal deterministic witness set for the magnitude of n.
    // We compare n to each bound using a u128 representation of n's low 128 bits.
    // For n that fits in 128 bits (i.e., high 128 bits are zero), this is exact.
    // For n > 2^128, the comparison falls through to the full witness set.
    let witnesses = select_witnesses(n);
    miller_rabin(n, witnesses)
}

/// Select the minimal witness set for ``n``.
///
/// Compares ``n`` against the published deterministic bounds.  Returns the
/// smallest witness set whose bound exceeds ``n``.  Falls back to the full
/// 13-witness set for ``n`` beyond all published bounds.
fn select_witnesses(n: &Uint<4>) -> &'static [u64] {
    // Extract n as a u128 for comparison against the bounds table.
    // Uint<4> is 256 bits; words() returns [lo64, ..., hi64] in little-endian order.
    let words = n.as_words();
    // words[0] is the least-significant 64-bit limb.
    let n_lo = words[0] as u128 | ((words[1] as u128) << 64);
    let n_hi = words[2] as u128 | ((words[3] as u128) << 64);

    // If the high 128 bits are non-zero, n > 2^128 > all table bounds.
    if n_hi != 0 {
        return FULL_WITNESSES;
    }

    // n fits in 128 bits; compare against each bound.
    for &(bound, witnesses) in DETERMINISTIC_WITNESSES {
        if n_lo < bound {
            return witnesses;
        }
    }

    FULL_WITNESSES
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn u(v: u64) -> Uint<4> {
        Uint::<4>::from(v)
    }

    // ── KAT: small primes and composites ──────────────────────────────────────

    /// Verify that ``is_prime`` correctly identifies small primes.
    #[test]
    fn small_primes_recognised() {
        let primes: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
        for &p in primes {
            assert!(is_prime(&u(p)), "expected {p} to be prime");
        }
    }

    /// Verify that ``is_prime`` correctly identifies small composites.
    #[test]
    fn small_composites_rejected() {
        let composites: &[u64] = &[0, 1, 4, 6, 8, 9, 10, 15, 25, 49, 77, 91, 121];
        for &c in composites {
            assert!(!is_prime(&u(c)), "expected {c} to be composite");
        }
    }

    // ── KAT: first 100 primes (2..=541) ──────────────────────────────────────

    /// Verify that all 100 primes up to 541 are recognised as prime.
    #[test]
    fn first_100_primes() {
        // The 100th prime is 541.
        let first_100: &[u64] = &[
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
            83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
            179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
            271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373,
            379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
            479, 487, 491, 499, 503, 509, 521, 523, 541,
        ];
        assert_eq!(first_100.len(), 100);
        for &p in first_100 {
            assert!(is_prime(&u(p)), "expected prime {p} to be recognised");
        }
    }

    // ── KAT: strong pseudoprimes for base 2 ──────────────────────────────────

    /// Verify that published strong pseudoprimes for base 2 are rejected by ``is_prime``.
    ///
    /// These numbers pass the Miller–Rabin test with witness 2 alone but are
    /// composite.  ``is_prime`` uses additional witnesses and must reject them.
    #[test]
    fn base2_pseudoprimes_rejected() {
        // Strong pseudoprimes to base 2 (OEIS A001262): 2047, 3277, 4033, 4681, 8321, ...
        let pseudoprimes: &[u64] = &[
            2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665, 80581,
            85489, 88357, 90751,
        ];
        for &p in pseudoprimes {
            assert!(!is_prime(&u(p)), "expected pseudoprime {p} to be rejected");
        }
    }

    // ── KAT: known 64-bit primes ──────────────────────────────────────────────

    /// Verify that known large primes are recognised.
    #[test]
    fn mersenne_primes_recognised() {
        // 2^31 - 1 = 2,147,483,647 (M31, Mersenne prime)
        assert!(is_prime(&u(2_147_483_647)));
        // 2^61 - 1 = 2,305,843,009,213,693,951 (M61, Mersenne prime)
        assert!(is_prime(&u(2_305_843_009_213_693_951)));
    }

    /// Verify that a known 64-bit composite near a Mersenne prime is rejected.
    #[test]
    fn near_mersenne_composite_rejected() {
        // 2^31 - 1 + 2 = 2,147,483,649 = 3 * 715,827,883 (composite)
        assert!(!is_prime(&u(2_147_483_649)));
    }

    // ── KAT: miller_rabin with explicit witnesses ─────────────────────────────

    /// Verify that ``miller_rabin`` with witness {2} accepts known primes.
    #[test]
    fn miller_rabin_base2_primes() {
        for &p in &[5u64, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] {
            assert!(
                miller_rabin(&u(p), &[2]),
                "miller_rabin base-2 should accept prime {p}"
            );
        }
    }

    /// Verify that ``miller_rabin`` with witness {2} rejects 2047 (strong pseudoprime
    /// to base 2 — this should PASS with base 2 alone, demonstrating why multiple
    /// witnesses are needed).
    #[test]
    fn miller_rabin_base2_accepts_2047() {
        // 2047 is a strong pseudoprime to base 2: miller_rabin with {2} returns true.
        assert!(
            miller_rabin(&u(2047), &[2]),
            "2047 is a strong pseudoprime to base 2; single-witness MR should accept it"
        );
        // But is_prime (multi-witness) correctly rejects it.
        assert!(!is_prime(&u(2047)), "is_prime must reject 2047");
    }

    /// Verify that ``miller_rabin`` with witnesses {2, 3} rejects 2047.
    #[test]
    fn miller_rabin_base23_rejects_2047() {
        assert!(
            !miller_rabin(&u(2047), &[2, 3]),
            "miller_rabin with bases {{2,3}} should reject 2047"
        );
    }

    // ── KAT: the twelve-base proof bound ──────────────────────────────────────

    /// ψ₁₂ is the least strong pseudoprime to the first twelve prime bases; it equals
    /// 399_165_290_221 · 798_330_580_441 and is caught only by base 41.
    #[test]
    fn psi_12_is_the_product_it_claims() {
        let p = 399_165_290_221u128;
        let q = 798_330_580_441u128;
        assert_eq!(p * q, TWELVE_BASE_PROOF_BOUND);
        assert!(is_prime(&Uint::<4>::from(p)));
        assert!(is_prime(&Uint::<4>::from(q)));
    }

    /// The first twelve bases accept ψ₁₂, so a twelve-base verdict is not a proof there.
    #[test]
    fn twelve_bases_accept_psi_12_but_base_41_rejects_it() {
        let psi = Uint::<4>::from(TWELVE_BASE_PROOF_BOUND);
        let twelve: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        assert!(miller_rabin(&psi, twelve));
        assert!(!miller_rabin(&psi, &[41]));
    }

    /// ``is_prime`` rejects ψ₁₂ because the thirteenth base is applied from there up.
    #[test]
    fn psi_12_rejected() {
        assert!(!is_prime(&Uint::<4>::from(TWELVE_BASE_PROOF_BOUND)));
    }

    /// A prime between the twelve-base and thirteen-base bounds is still accepted.
    #[test]
    fn prime_above_twelve_base_bound_accepted() {
        // 2^79 − 67 is prime.
        let n = (1u128 << 79) - 67;
        assert!(n > TWELVE_BASE_PROOF_BOUND && n < MILLER_RABIN_PROOF_BOUND);
        assert!(is_prime(&Uint::<4>::from(n)));
    }

    // ── Property: no false negatives in 2..1000 ───────────────────────────────

    /// Verify that ``is_prime`` has no false negatives for n in 2..1000.
    ///
    /// Cross-checks against a brute-force trial-division reference.
    #[test]
    fn no_false_negatives_up_to_1000() {
        for n in 2u64..1000 {
            let expected = is_prime_brute(n);
            let got = is_prime(&u(n));
            assert_eq!(
                got, expected,
                "is_prime({n}) = {got}, brute-force says {expected}"
            );
        }
    }

    /// Brute-force primality check via trial division (reference implementation).
    fn is_prime_brute(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        if n == 2 {
            return true;
        }
        if n % 2 == 0 {
            return false;
        }
        let mut i = 3u64;
        while i * i <= n {
            if n % i == 0 {
                return false;
            }
            i += 2;
        }
        true
    }
}
