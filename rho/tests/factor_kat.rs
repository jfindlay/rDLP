//! Known-answer tests for Pollard rho factorization.
//!
//! Each test factors a known semiprime (product of two primes) and verifies
//! that all variants (Floyd, Brent, brent_batched, and the parallel `factor`)
//! return a correct non-trivial factor.  The check is: `1 < d < n` and
//! `n % d == 0`.  `factor` additionally reports which path produced the divisor, and its
//! error variants (invalid input, prime, unsupported, exhausted budget, inconsistent candidate)
//! are pinned separately.
//!
//! The semiprimes span about 20–64 bits; the supported range reaches 2^80, and the only
//! composite above 64 bits that is factored is ψ₁₂ ≈ 2^78.1
//! (`factor_splits_psi_12_instead_of_certifying_it`); primes up to 2^80 are only certified.
//! For each semiprime the factors are listed but not asserted — rho may return
//! either factor; both are correct.

use rho::factor::{
    brent, brent_batched, factor, floyd, FactorError, FactorEvidence, Factorization,
    MAX_SUPPORTED,
};

/// Step bound generous enough that no known semiprime below exhausts it.
const STEPS: u64 = 1 << 28;

/// Assert that `d` is a non-trivial factor of `n`.
fn check(d: Option<u128>, n: u128, label: &str) {
    let d = d.unwrap_or_else(|| panic!("{label}: returned None for n={n}"));
    assert!(d > 1, "{label}: factor {d} ≤ 1 for n={n}");
    assert!(d < n, "{label}: factor {d} ≥ n={n}");
    assert_eq!(n % d, 0, "{label}: {d} does not divide {n}");
}

/// Unwrap a successful `factor`, asserting the divisor invariant.
fn factored(n: u128, max_c: u64) -> Factorization {
    let found = factor(n, max_c, 128)
        .unwrap_or_else(|e| panic!("factor({n}, {max_c}) failed with {e:?}"));
    check(Some(found.d), n, &format!("factor n={n}"));
    found
}

/// Known semiprimes: `(n, p, q)` where `n = p * q` and both `p`, `q` are prime.
/// All `n < 2^80` per the supported range of `mulmod`.
const SEMIPRIMES: &[(u128, u128, u128)] = &[
    // Small sanity cases (easy for any variant)
    (15, 3, 5),
    (77, 7, 11),
    (221, 13, 17),
    (3_599, 59, 61),
    // ~20 bit
    (1_000_003 * 1_000_033, 1_000_003, 1_000_033),
    // ~30 bit
    (32_452_843 * 32_452_867, 32_452_843, 32_452_867),
    // ~40 bit
    (549_755_813_881 * 2, 549_755_813_881, 2), // even; trivial for factor()
    // ~40 bit proper semiprime
    (1_000_000_007 * 1_000_000_009, 1_000_000_007, 1_000_000_009),
    // ~50 bit (two ~25-bit primes)
    (33_554_467 * 33_554_473, 33_554_467, 33_554_473),
    // ~60 bit (two ~30-bit primes)
    (1_073_741_827 * 1_073_741_789, 1_073_741_827, 1_073_741_789),
    // ~64 bit (two ~32-bit primes): x < n < 2^64 so x*x < 2^128 — still needs carrying_mul.
    (4_294_967_291 * 4_294_967_311, 4_294_967_291, 4_294_967_311),
];

#[test]
fn floyd_known_semiprimes() {
    for &(n, _, _) in SEMIPRIMES {
        // Skip the even-number entry: floyd may loop forever on it (rho works on
        // odd composites; the even fast-path lives in factor()).
        if n % 2 == 0 {
            continue;
        }
        let d = floyd(n, 1, 2, STEPS);
        check(d, n, &format!("floyd n={n}"));
    }
}

#[test]
fn brent_known_semiprimes() {
    for &(n, _, _) in SEMIPRIMES {
        if n % 2 == 0 {
            continue;
        }
        let d = brent(n, 1, 2, STEPS);
        check(d, n, &format!("brent n={n}"));
    }
}

#[test]
fn brent_batched_known_semiprimes() {
    for &(n, _, _) in SEMIPRIMES {
        if n % 2 == 0 {
            continue;
        }
        let d = brent_batched(n, 1, 2, 128, STEPS);
        check(d, n, &format!("brent_batched n={n}"));
    }
}

#[test]
fn primitives_fail_when_the_step_bound_is_spent() {
    let n: u128 = 1_073_741_827 * 1_073_741_789;
    // A zero bound allows no evaluation, so no factor can be reported.
    assert_eq!(floyd(n, 1, 2, 0), None);
    assert_eq!(brent(n, 1, 2, 0), None);
    assert_eq!(brent_batched(n, 1, 2, 128, 0), None);
    // A bound far below the ~√p steps rho needs also fails rather than looping on.
    assert_eq!(floyd(n, 1, 2, 100), None);
    assert_eq!(brent(n, 1, 2, 100), None);
    assert_eq!(brent_batched(n, 1, 2, 128, 100), None);
}

#[test]
fn brent_batched_zero_batch_size_terminates() {
    // A zero batch is treated as one step per GCD; it must not spin without progress.
    let n: u128 = 1_000_003 * 1_000_033;
    check(brent_batched(n, 1, 2, 0, STEPS), n, "brent_batched batch_size=0");
}

#[test]
fn factor_parallel_known_semiprimes() {
    for &(n, _, _) in SEMIPRIMES {
        factored(n, 8);
    }
}

#[test]
fn factor_even_composites_carry_even_evidence() {
    for n in [4u128, 6, 2 * 1_000_000_007, 549_755_813_881 * 2, (1 << 79) + 2] {
        let found = factored(n, 4);
        assert_eq!((found.d, found.evidence), (2, FactorEvidence::Even), "n={n}");
    }
    // The even path runs no walk, so it does not depend on the constant budget.
    assert_eq!(factored(4, 0).evidence, FactorEvidence::Even);
}

#[test]
fn factor_odd_composites_carry_rho_evidence_that_reproduces_the_divisor() {
    for &(n, _, _) in SEMIPRIMES {
        if n % 2 == 0 {
            continue;
        }
        let found = factored(n, 8);
        let FactorEvidence::Rho { c } = found.evidence else {
            panic!("n={n}: odd composite reported {:?}", found.evidence);
        };
        assert!((1..=8).contains(&c), "n={n}: constant {c} outside 1..=8");
        // The evidence alone reproduces the divisor: the same walk gives the same gcd.
        assert_eq!(brent_batched(n, u128::from(c), 2, 128, STEPS), Some(found.d), "n={n}");
    }
}

#[test]
fn factor_squares_and_prime_powers_give_proper_divisors() {
    let p = 1_000_003u128;
    let cases = [4u128, 8, 9, 25, 27, 49, 81, 121, 343, p * p, 3u128.pow(40), 7u128.pow(25)];
    for n in cases {
        let found = factored(n, 16);
        assert!(found.d > 1 && found.d < n && n % found.d == 0, "n={n}");
    }
}

#[test]
fn factor_small_odd_composites_all_succeed() {
    for n in (9u128..3_000).step_by(2) {
        if (3..n).step_by(2).take_while(|d| d * d <= n).any(|d| n % d == 0) {
            factored(n, 16);
        }
    }
}

#[test]
fn factor_rejects_zero_and_one() {
    for n in [0u128, 1] {
        assert_eq!(factor(n, 8, 128), Err(FactorError::InvalidInput), "n={n}");
    }
}

#[test]
fn factor_reports_primes_and_never_returns_two_as_a_factor_of_two() {
    let primes = [
        2u128,
        3,
        5,
        7,
        13,
        1_000_000_007,
        (1 << 61) - 1,
        18_446_744_073_709_551_557, // 2^64 − 59, largest prime below 2^64
    ];
    for n in primes {
        assert_eq!(factor(n, 8, 128), Err(FactorError::Prime), "n={n}");
    }
    // A prime is classified before any constant is tried, so the budget is irrelevant.
    assert_eq!(factor(1_000_000_007, 0, 128), Err(FactorError::Prime));
}

#[test]
fn factor_rejects_sizes_beyond_mulmod() {
    // 2^89 − 1 is a prime above the range: the size, not primality, decides the outcome.
    for n in [MAX_SUPPORTED, MAX_SUPPORTED + 1, (1u128 << 89) - 1, u128::MAX] {
        assert_eq!(factor(n, 8, 128), Err(FactorError::Unsupported), "n={n}");
    }
    assert_eq!(MAX_SUPPORTED, 1u128 << 80);
}

/// Primes between ψ₁₂ ≈ 2^78.1 and 2^80 are certified: the thirteen-base Miller–Rabin test that
/// `is_prime` runs there is a proof below ψ₁₃ ≈ 2^81.5, and every supported `n` is below that.
#[test]
fn factor_certifies_primes_up_to_the_supported_limit() {
    const PSI_12: u128 = 318_665_857_834_031_151_167_461;
    // 2^79 − 67 and 2^80 − 65 (the largest prime below 2^80).
    for n in [(1u128 << 79) - 67, MAX_SUPPORTED - 65] {
        assert!(PSI_12 < n && n < MAX_SUPPORTED, "n={n}");
        assert_eq!(factor(n, 8, 128), Err(FactorError::Prime), "n={n}");
    }
    // The classification precedes the walks, so a zero budget changes nothing.
    assert_eq!(factor((1u128 << 79) - 67, 0, 128), Err(FactorError::Prime));
}

/// ψ₁₂ = 399_165_290_221 · 798_330_580_441 is composite and a strong pseudoprime to the first
/// twelve prime bases; it must never be reported prime, and rho splits it.
#[test]
fn factor_splits_psi_12_instead_of_certifying_it() {
    const P: u128 = 399_165_290_221;
    const Q: u128 = 798_330_580_441;
    let n = P * Q;
    assert_eq!(n, 318_665_857_834_031_151_167_461);
    let found = factored(n, 8);
    assert!(found.d == P || found.d == Q, "d={}", found.d);
}

#[test]
fn factor_with_no_constants_is_budget_exhausted() {
    for n in [9u128, 15, 1_000_003 * 1_000_033] {
        assert_eq!(
            factor(n, 0, 128),
            Err(FactorError::BudgetExhausted { constants: 0 }),
            "n={n}"
        );
    }
}

#[test]
fn factor_counts_the_constants_whose_walks_failed() {
    // 11899 = 73 · 163: the walks for c = 1..=4 all collide modulo both primes at once
    // (gcd = n), so a budget of four constants is spent, while a larger budget succeeds.
    let n: u128 = 73 * 163;
    assert_eq!(factor(n, 4, 128), Err(FactorError::BudgetExhausted { constants: 4 }));
    factored(n, 16);
}

/// Compile-time check that an error enum stays `Copy`.
fn assert_copy<T: Copy>() {}

#[test]
fn factor_error_inconsistent_carries_the_candidate_and_stays_copy() {
    assert_copy::<FactorError>();
    let e = FactorError::Inconsistent { d: 7 };
    let copy = e;
    assert_eq!(e, copy);
    assert_ne!(e, FactorError::Inconsistent { d: 8 });
    assert_ne!(e, FactorError::BudgetExhausted { constants: 0 });
}

#[test]
fn factor_never_reports_inconsistent_on_real_inputs() {
    // The gcd with n always divides n, so real inputs never trip the divisor check: semiprimes,
    // prime powers, and small odd composites either factor or spend a budget, but never report
    // an invariant breach.
    let mut inputs: Vec<u128> = SEMIPRIMES.iter().map(|&(n, _, _)| n).collect();
    inputs.extend([9, 25, 49, 121, 27, 125, 343, 3u128.pow(20), 1_000_003u128.pow(2)]);
    inputs.extend((9u128..2_000).step_by(2));
    for n in inputs {
        for max_c in [1u64, 4, 16] {
            match factor(n, max_c, 128) {
                Err(FactorError::Inconsistent { d }) => {
                    panic!("factor({n}, {max_c}) reported inconsistent candidate {d}")
                }
                Ok(found) => check(Some(found.d), n, &format!("factor n={n}")),
                Err(_) => {}
            }
        }
    }
}

#[test]
fn brent_batched_various_batch_sizes() {
    // Verify that different batch sizes all produce correct results.
    let n: u128 = 1_073_741_827 * 1_073_741_789;
    for &bs in &[1usize, 16, 64, 128, 256] {
        let d = brent_batched(n, 1, 2, bs, STEPS);
        check(d, n, &format!("brent_batched batch_size={bs}"));
    }
}
