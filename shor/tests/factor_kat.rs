//! Known-answer tests for Shor's factoring algorithm and its order-finding step.
//!
//! # What is tested
//!
//! 1. **Order KATs** — the classical order `ord_a(N)` matches published values:
//!    - `ord₂(15) = 4` (2^4 = 16 ≡ 1 mod 15)
//!    - `ord₇(15) = 4` (7^4 = 2401 ≡ 1 mod 15)
//!    - `ord₂(21) = 6` (2^6 = 64 ≡ 1 mod 21)
//!    - `ord₂(35) = 12` (2^12 = 4096 ≡ 1 mod 35)
//!    - `ord₄(35) = 6` (4^6 = 4096 ≡ 1 mod 35)
//!    - `ord₂(91) = 12` (2^12 = 4096 ≡ 1 mod 91)
//!
//! 2. **Continued-fraction KATs** — a known measured phase `s/2^t` recovers the correct order
//!    `r` via the convergents (deterministic classical KAT, no quantum circuit needed).
//!
//! 3. **End-to-end factoring** — `factor(15)→{3,5}`, `factor(21)→{3,7}`, `factor(35)→{5,7}`
//!    with fixed seeds so the run is reproducible. Each asserts the `OrderFinding` evidence
//!    (base, seed, measured phase, order) and rechecks `a^r ≡ 1` and the gcd step from it.
//!
//! 4. **91 ceiling-stress KAT** — `factor(91)→{7,13}`. N=91 uses 14 qubits (exp_len=7,
//!    work=7), well within the ~25-qubit simulator ceiling. This is the principle-4
//!    ceiling-stress case: the algorithm is correct; the qubit budget fits.
//!
//! 5. **Classification** — primes and `n ≤ 1` are `NotComposite`, even `n` is `Even`, prime powers
//!    are `PrimePower`, and a circuit beyond the qubit ceiling is `Unsupported`; none of these
//!    runs a circuit.
//!
//! # Qubit budgets (ancilla-free permutation-synthesis implementation)
//!
//! | N  | n_bits(N) | exp_len (t) | work (n) | total qubits |
//! |----|-----------|-------------|----------|--------------|
//! | 15 | 4         | 4           | 4        | 8            |
//! | 21 | 5         | 5           | 5        | 10           |
//! | 35 | 6         | 6           | 6        | 12           |
//! | 91 | 7         | 7           | 7        | 14           |
//!
//! All within the ~25-qubit ceiling.
//!
//! # Published reference values
//!
//! Orders (from Shor 1994, Nielsen & Chuang §5.3):
//!   - ord₂(15) = 4:  2^1=2, 2^2=4, 2^3=8, 2^4=16≡1 (mod 15)
//!   - ord₇(15) = 4:  7^1=7, 7^2=49≡4, 7^3=28≡13, 7^4=91≡1 (mod 15)
//!   - ord₂(21) = 6:  2^6=64≡1 (mod 21)
//!   - ord₂(35) = 12: 2^12=4096≡1 (mod 35)
//!   - ord₄(35) = 6:  4^6=4096≡1 (mod 35)
//!   - ord₂(91) = 12: 2^12=4096≡1 (mod 91)
//!
//! Continued-fraction recovery (N=15, a=2, t=8):
//!   - Measured s=64 (= 1·256/4): phase 64/256 = 1/4, convergent 1/4, r=4 ✓
//!   - Measured s=128 (= 2·256/4): phase 128/256 = 1/2, convergent 1/2, r=2 (not order)
//!     → convergent 1/4 not reached; but 2^2=4≠1, so s=128 gives r=None
//!   - Measured s=32 (a=7): phase 32/256 = 1/8, denominator 8 satisfies 7^8≡1 but is twice the
//!     order; the recovered order is reduced to 4
//!   - Measured s=192 (= 3·256/4): phase 192/256 = 3/4, convergent 3/4, r=4 ✓
//!
//! Reference: Shor (1994); Nielsen & Chuang §5.3.

use shor::arith::{mod_pow, n_bits};
use shor::shor::{
    convergents, exp_len_for, factor, find_order, gcd, order_from_phase, ShorEvidence,
    ShorFactorError, MAX_QUBITS,
};

// ── helpers ───────────────────────────────────────────────────────────────────

/// Least `r > 0` with `a^r ≡ 1 mod n`, by exhaustive scan (`gcd(a, n) = 1`).
fn brute_force_order(a: u64, n: u64) -> u64 {
    (1..n).find(|&r| mod_pow(a, r, n) == 1).expect("a is a unit mod n")
}

/// Assert that `factor(n, seed)` returns `{p, q}` with `p * q = n` and the `OrderFinding`
/// evidence `(a, evidence seed, s, r)`, then recheck that evidence from its fields alone.
///
/// The returned pair is normalized to `(min, max)` for comparison.
fn assert_factor(n: u64, seed: u64, expected: (u64, u64), order_run: (u64, u64, u64, u64)) {
    let found = factor(n, seed)
        .unwrap_or_else(|e| panic!("factor({n}, seed={seed}) failed: {e:?}; want {expected:?}"));
    let got = if found.p <= found.q { (found.p, found.q) } else { (found.q, found.p) };
    assert_eq!(got, expected, "factor({n}, seed={seed})");
    assert_eq!(found.p * found.q, n, "factor({n}): {} * {} ≠ {n}", found.p, found.q);
    assert!(found.p > 1 && found.q > 1, "factor({n}): trivial factor in {got:?}");

    let ShorEvidence::OrderFinding { a, seed: run_seed, s, r } = found.evidence else {
        panic!("factor({n}, seed={seed}): want OrderFinding evidence, got {:?}", found.evidence);
    };
    assert_eq!((a, run_seed, s, r), order_run, "factor({n}, seed={seed}) evidence");

    // The recorded phase gives the recorded order, which is the true multiplicative order.
    assert_eq!(order_from_phase(s, exp_len_for(n), a, n), Some(r), "order from recorded phase");
    assert_eq!(mod_pow(a, r, n), 1, "{a}^{r} mod {n} ≠ 1");
    assert_eq!(r, brute_force_order(a, n), "r is not the least exponent for a={a}, n={n}");

    // The gcd step: r even, a^(r/2) ≢ -1, and p is one of gcd(a^(r/2) ∓ 1, n).
    assert_eq!(r % 2, 0, "order {r} is odd");
    let half = mod_pow(a, r / 2, n);
    assert_ne!(half, n - 1, "a^(r/2) ≡ -1 mod {n}");
    assert!(
        found.p == gcd(half - 1, n) || found.p == gcd(half + 1, n),
        "p={} is neither gcd({half}-1, {n}) nor gcd({half}+1, {n})",
        found.p
    );
}


// ── KAT group 1: order KATs (classical verification) ─────────────────────────
//
// These KATs verify that the classical order function is correct.
// The order is verified by checking a^r ≡ 1 mod N.

/// ord₂(15) = 4: 2^4 = 16 ≡ 1 (mod 15).
#[test]
fn order_classical_2_mod_15() {
    assert_eq!(mod_pow(2, 4, 15), 1, "2^4 mod 15 should be 1");
    assert_ne!(mod_pow(2, 1, 15), 1, "2^1 mod 15 should not be 1");
    assert_ne!(mod_pow(2, 2, 15), 1, "2^2 mod 15 should not be 1");
    assert_ne!(mod_pow(2, 3, 15), 1, "2^3 mod 15 should not be 1");
}

/// ord₇(15) = 4: 7^4 = 2401 ≡ 1 (mod 15).
#[test]
fn order_classical_7_mod_15() {
    // 7^1=7, 7^2=49≡4, 7^3=28≡13, 7^4=91≡1 (mod 15)
    assert_eq!(mod_pow(7, 4, 15), 1, "7^4 mod 15 should be 1");
    assert_ne!(mod_pow(7, 1, 15), 1, "7^1 mod 15 should not be 1");
    assert_ne!(mod_pow(7, 2, 15), 1, "7^2 mod 15 should not be 1");
    assert_ne!(mod_pow(7, 3, 15), 1, "7^3 mod 15 should not be 1");
}

/// ord₂(21) = 6: 2^6 = 64 ≡ 1 (mod 21).
#[test]
fn order_classical_2_mod_21() {
    assert_eq!(mod_pow(2, 6, 21), 1, "2^6 mod 21 should be 1");
    for k in 1..6 {
        assert_ne!(mod_pow(2, k, 21), 1, "2^{k} mod 21 should not be 1 (order is 6)");
    }
}

/// ord₂(35) = 12: 2^12 = 4096 ≡ 1 (mod 35).
#[test]
fn order_classical_2_mod_35() {
    assert_eq!(mod_pow(2, 12, 35), 1, "2^12 mod 35 should be 1");
    // Verify no smaller order: check divisors of 12 (1, 2, 3, 4, 6).
    for k in [1u64, 2, 3, 4, 6] {
        assert_ne!(mod_pow(2, k, 35), 1, "2^{k} mod 35 should not be 1 (order is 12)");
    }
}

/// ord₄(35) = 6: 4^6 = 4096 ≡ 1 (mod 35).
#[test]
fn order_classical_4_mod_35() {
    assert_eq!(mod_pow(4, 6, 35), 1, "4^6 mod 35 should be 1");
    for k in [1u64, 2, 3] {
        assert_ne!(mod_pow(4, k, 35), 1, "4^{k} mod 35 should not be 1 (order is 6)");
    }
}

/// ord₂(91) = 12: 2^12 = 4096 ≡ 1 (mod 91).
#[test]
fn order_classical_2_mod_91() {
    // 4096 mod 91: 91*45 = 4095, so 4096 mod 91 = 1. ✓
    assert_eq!(mod_pow(2, 12, 91), 1, "2^12 mod 91 should be 1");
    for k in [1u64, 2, 3, 4, 6] {
        assert_ne!(mod_pow(2, k, 91), 1, "2^{k} mod 91 should not be 1 (order is 12)");
    }
}

/// ord₃(91) = 6: 3^6 ≡ 1 (mod 91).
///
/// 3 is coprime to 91 (gcd(3, 91) = 1). 3^6 = 729 = 8×91 + 1 = 729. ✓
/// This is a useful base for factoring 91: ord₃(91) = 6 (even), 3^3 = 27 mod 91 = 27.
/// gcd(27-1, 91) = gcd(26, 91) = 13. ✓
#[test]
fn order_classical_3_mod_91() {
    // 3^6 = 729. 729 / 91 = 8 remainder 1. So 3^6 mod 91 = 1. ✓
    assert_eq!(mod_pow(3, 6, 91), 1, "3^6 mod 91 should be 1");
    // Verify no smaller order: check divisors of 6 (1, 2, 3).
    assert_ne!(mod_pow(3, 1, 91), 1, "3^1 mod 91 should not be 1");
    assert_ne!(mod_pow(3, 2, 91), 1, "3^2 mod 91 should not be 1");
    assert_ne!(mod_pow(3, 3, 91), 1, "3^3 mod 91 should not be 1");
    // Factor extraction: 3^3 = 27, gcd(27-1, 91) = gcd(26, 91) = 13. ✓
    assert_eq!(gcd(mod_pow(3, 3, 91) - 1, 91), 13, "gcd(3^3 - 1, 91) should be 13");
}

/// n_bits helper: correct bit counts for factoring targets.
#[test]
fn n_bits_for_factoring_targets() {
    // n_bits(N) = ⌈log₂(N)⌉ = bits needed to represent N-1
    assert_eq!(n_bits(15), 4, "n_bits(15) should be 4 (15 = 0b1111)");
    assert_eq!(n_bits(21), 5, "n_bits(21) should be 5 (21 = 0b10101)");
    assert_eq!(n_bits(35), 6, "n_bits(35) should be 6 (35 = 0b100011)");
    assert_eq!(n_bits(91), 7, "n_bits(91) should be 7 (91 = 0b1011011)");
}

/// exp_len_for: correct exponent register sizes.
///
/// Uses t = n_bits(N) (the ancilla-free permutation-synthesis qubit budget).
#[test]
fn exp_len_for_factoring_targets() {
    assert_eq!(exp_len_for(15), 4, "exp_len for N=15 should be 4 (n_bits(15))");
    assert_eq!(exp_len_for(21), 5, "exp_len for N=21 should be 5 (n_bits(21))");
    assert_eq!(exp_len_for(35), 6, "exp_len for N=35 should be 6 (n_bits(35))");
    assert_eq!(exp_len_for(91), 7, "exp_len for N=91 should be 7 (n_bits(91))");
}

/// gcd helper: basic correctness.
#[test]
fn gcd_basic() {
    assert_eq!(gcd(12, 8), 4);
    assert_eq!(gcd(15, 5), 5);
    assert_eq!(gcd(7, 13), 1);
    assert_eq!(gcd(0, 5), 5);
    assert_eq!(gcd(5, 0), 5);
    assert_eq!(gcd(6, 15), 3);
    // gcd(14, 91): 14 = 2×7, 91 = 7×13 → gcd = 7
    assert_eq!(gcd(14, 91), 7);
    // gcd(12, 91): 12 = 4×3, 91 = 7×13 → gcd = 1 (no common factors)
    assert_eq!(gcd(12, 91), 1);
    // gcd(26, 91): 26 = 2×13, 91 = 7×13 → gcd = 13
    assert_eq!(gcd(26, 91), 13);
}

// ── KAT group 2: continued-fraction KATs (classical, deterministic) ───────────
//
// Given a known measured phase s/2^t, verify that order_from_phase recovers the correct order r.
// These are purely classical KATs — no quantum circuit is run.

/// Continued-fraction convergents of 1/4 (s=64, t=8): denominator 4 is a convergent.
///
/// Phase 64/256 = 1/4. Convergents: 0/1, 1/4. Denominator 4 satisfies 2^4 ≡ 1 mod 15.
#[test]
fn cf_convergents_1_over_4() {
    // 64/256 = 1/4. Continued fraction: 0 + 1/4 → convergents: 0/1, 1/4.
    let convs = convergents(64, 256, 15);
    let denoms: Vec<u64> = convs.iter().map(|&(_, d)| d).collect();
    assert!(
        denoms.contains(&4),
        "convergents(64, 256, 15) should include denominator 4, got {denoms:?}"
    );
}

/// Continued-fraction convergents of 3/4 (s=192, t=8): denominator 4 is a convergent.
///
/// Phase 192/256 = 3/4. Convergents include 3/4. Denominator 4 satisfies 2^4 ≡ 1 mod 15.
#[test]
fn cf_convergents_3_over_4() {
    let convs = convergents(192, 256, 15);
    let denoms: Vec<u64> = convs.iter().map(|&(_, d)| d).collect();
    assert!(
        denoms.contains(&4),
        "convergents(192, 256, 15) should include denominator 4, got {denoms:?}"
    );
}

/// order_from_phase: s=64, t=8, a=2, N=15 → r=4.
///
/// Phase 64/256 = 1/4 ≈ 1/r with r=4. Convergent 1/4 has denominator 4.
/// 2^4 = 16 ≡ 1 mod 15. ✓
#[test]
fn cf_order_from_phase_n15_a2_s64() {
    // s=64, t=8: phase = 64/256 = 1/4 ≈ k/r with k=1, r=4.
    let r = order_from_phase(64, 8, 2, 15);
    assert_eq!(r, Some(4), "order_from_phase(64, 8, 2, 15) should be Some(4), got {r:?}");
}

/// order_from_phase: s=192, t=8, a=2, N=15 → r=4.
///
/// Phase 192/256 = 3/4 ≈ 3/r with r=4. Convergent 3/4 has denominator 4.
/// 2^4 = 16 ≡ 1 mod 15. ✓
#[test]
fn cf_order_from_phase_n15_a2_s192() {
    let r = order_from_phase(192, 8, 2, 15);
    assert_eq!(r, Some(4), "order_from_phase(192, 8, 2, 15) should be Some(4), got {r:?}");
}

/// order_from_phase: s=0, t=8, a=2, N=15 → None (uninformative measurement).
///
/// s=0 gives phase 0, which is uninformative (no period information).
#[test]
fn cf_order_from_phase_zero_phase() {
    let r = order_from_phase(0, 8, 2, 15);
    assert_eq!(r, None, "order_from_phase(0, ...) should be None (uninformative)");
}

/// order_from_phase: s=32, t=8, a=7, N=15 → r=4, not the convergent denominator 8.
///
/// Phase 32/256 = 1/8. Convergents: 0/1, 1/8. The only denominator with 7^d ≡ 1 mod 15 is
/// d = 8, a multiple of ord₇(15) = 4. It is reduced by dividing out the prime 2 while the
/// identity holds: 7^4 ≡ 1, but 7^2 = 49 ≡ 4, so the order is 4.
#[test]
fn cf_order_from_phase_n15_a7_s32() {
    assert_eq!(mod_pow(7, 8, 15), 1, "denominator 8 satisfies 7^8 ≡ 1 mod 15");
    let r = order_from_phase(32, 8, 7, 15);
    assert_eq!(r, Some(4), "order_from_phase(32, 8, 7, 15) should be Some(4), got {r:?}");
}

/// Whenever `order_from_phase` returns `Some(r)`, `r` is the multiplicative order: every phase
/// numerator, every base coprime to N, for the four register widths the factoring targets use.
#[test]
fn cf_order_from_phase_is_always_the_least_order() {
    let mut recovered = 0usize;
    for n in [15u64, 21, 35, 91] {
        let t = exp_len_for(n);
        for a in (2..n).filter(|&a| gcd(a, n) == 1) {
            let order = brute_force_order(a, n);
            for s in 0..(1u64 << t) {
                if let Some(r) = order_from_phase(s, t, a, n) {
                    assert_eq!(r, order, "order_from_phase({s}, {t}, {a}, {n}) = {r}");
                    recovered += 1;
                }
            }
        }
    }
    assert!(recovered > 0, "no phase recovered any order");
}

/// order_from_phase: s=64, t=8, a=7, N=15 → r=4.
///
/// Phase 64/256 = 1/4. Convergents: 0/1, 1/4. Denominator 4: 7^4 ≡ 1 mod 15. ✓
#[test]
fn cf_order_from_phase_n15_a7_s64() {
    let r = order_from_phase(64, 8, 7, 15);
    assert_eq!(r, Some(4), "order_from_phase(64, 8, 7, 15) should be Some(4), got {r:?}");
}

/// order_from_phase: s=5, t=5, a=2, N=21 → r=6.
///
/// Phase 5/32 ≈ 1/6 (since 32/6 ≈ 5.33). Convergents of 5/32 include denominator 6.
/// 2^6 = 64 ≡ 1 mod 21. ✓
/// Note: exp_len for N=21 is 5 (n_bits(21)), so 2^t = 32.
#[test]
fn cf_order_from_phase_n21_a2_s5() {
    // t=5 for N=21 (exp_len = n_bits(21) = 5), 2^5 = 32.
    // 32/6 ≈ 5.33, so s=5 ≈ 1*32/6.
    // Convergents of 5/32: 32=6×5+2, 5=2×2+1, 2=2×1. CF=[0;6,2,2].
    // Convergents: 0/1, 1/6, 2/13, 5/32. Denominator 6: 2^6=64 mod 21=1. ✓
    let r = order_from_phase(5, 5, 2, 21);
    assert_eq!(r, Some(6), "order_from_phase(5, 5, 2, 21) should be Some(6), got {r:?}");
}

/// order_from_phase: s=27, t=5, a=2, N=21 → r=6.
///
/// Phase 27/32 ≈ 5/6 (since 5*32/6 ≈ 26.67). Convergents of 27/32 include denominator 6.
/// 2^6 = 64 ≡ 1 mod 21. ✓
#[test]
fn cf_order_from_phase_n21_a2_s27() {
    // t=5, 2^t=32. s=27 ≈ 5*32/6. Phase 27/32 ≈ 5/6.
    // Convergents of 27/32: 32=1×27+5, 27=5×5+2, 5=2×2+1, 2=2×1.
    // CF=[0;1,5,2,2]. Convergents: 0/1, 1/1, 5/6, 11/13, 27/32. Denominator 6: ✓
    let r = order_from_phase(27, 5, 2, 21);
    assert_eq!(r, Some(6), "order_from_phase(27, 5, 2, 21) should be Some(6), got {r:?}");
}

// ── KAT group 3: end-to-end factoring (quantum circuit + classical post-processing) ──
//
// These KATs run the full Shor factoring algorithm end-to-end with fixed seeds.
// Seeds are chosen to land a successful (even-order, nontrivial-gcd) run on the first circuit
// run (base a=2, evidence seed = the given seed).

/// factor(15, seed=0) → {3, 5}.
///
/// N=15 = 3×5. Uses 8 qubits (exp_len=4, work=4). Seed 0 measures s=4 (phase 4/16 = 1/4) for
/// a=2, recovering ord₂(15)=4. Then 2^2 = 4 ≢ -1 and gcd(3, 15) = 3. Factor: 15/3 = 5. ✓
#[test]
fn factor_15_seed0() {
    assert_factor(15, 0, (3, 5), (2, 0, 4, 4));
}

/// factor(21, seed=0) → {3, 7}.
///
/// N=21 = 3×7. Uses 10 qubits (exp_len=5, work=5). Seed 0 measures s=27 (≈ 5×32/6) for a=2,
/// recovering ord₂(21)=6. Then 2^3 = 8 ≢ -1 and gcd(7, 21) = 7. Factor: 21/7 = 3. ✓
#[test]
fn factor_21_seed0() {
    assert_factor(21, 0, (3, 7), (2, 0, 27, 6));
}

/// factor(35, seed=7) → {5, 7}.
///
/// N=35 = 5×7. Uses 12 qubits (exp_len=6, work=6). Seed 7 measures s=59 (phase 59/64 ≈ 11/12)
/// for a=2, recovering ord₂(35)=12. Then 2^6 = 64 mod 35 = 29 ≠ 34 and gcd(28, 35) = 7.
/// Factor: 35/7 = 5. ✓ (Seed 8 does not recover an order on its first run.)
#[test]
fn factor_35_seed7() {
    assert_factor(35, 7, (5, 7), (2, 7, 59, 12));
}

// ── KAT group 4: 91 ceiling-stress KAT ───────────────────────────────────────
//
// N=91 = 7×13. Uses 14 qubits (exp_len=7, work=7), well within the ~25-qubit ceiling.
// This is the principle-4 ceiling-stress case: the algorithm is correct and the qubit
// budget fits. The ancilla-free permutation-synthesis implementation keeps the budget at
// exp_len+n=14.

/// Seed finder for N=91 (run once in release mode to find the right seed).
///
/// This test finds the first seed where run_order_finding_circuit(2, layout_91, seed)
/// gives a phase that recovers the order of 2 mod 91. Used to determine the seed for
/// the factor_91 KAT. Marked #[ignore] so it doesn't run in the normal test suite.
#[test]
#[ignore]
fn find_seed_for_91() {
    use shor::arith::ModExpLayout;
    use shor::shor::run_order_finding_circuit;
    let n = 91u64;
    let a = 2u64;
    let exp_len = exp_len_for(n);
    let layout = ModExpLayout::standard(n, exp_len);
    for seed in 0u64..100 {
        let s = run_order_finding_circuit(a, &layout, seed);
        let r = order_from_phase(s, exp_len, a, n);
        eprintln!("seed={seed}: s={s}, order={r:?}");
        if let Some(r_val) = r {
            assert_eq!(mod_pow(a, r_val, n), 1, "2^{r_val} mod 91 should be 1");
            eprintln!("FOUND: seed={seed}, s={s}, r={r_val}");
            return;
        }
    }
    panic!("No useful seed found in 0..100");
}

/// factor(91, seed=1) → {7, 13}.
///
/// N=91 = 7×13. Uses 14 qubits (exp_len=7, work=7). Seed 1.
/// Published factorization: 91 = 7 × 13.
/// Qubit budget: 14 qubits, within the ~25-qubit ceiling (principle 4).
///
/// Seed 1 is chosen because run_order_finding_circuit(2, layout_91, 1) measures s=53
/// (≈ 5×128/12), which recovers ord₂(91)=12 via continued fractions. Then:
/// 2^6 = 64 mod 91 = 64 ≠ 90, gcd(63, 91) = 7. Factor: 91/7 = 13. ✓
///
/// Cost: a 14-qubit circuit simulation, a few seconds in debug.
#[test]
fn factor_91_seed1() {
    assert_factor(91, 1, (7, 13), (2, 1, 53, 12));
}

/// Qubit budget for N=91 is within the ~25-qubit ceiling.
///
/// This test documents the principle-4 resource-scale annotation: N=91 requires 14 qubits
/// (exp_len=7, work=7), which fits within the ~25-qubit simulator ceiling. The algorithm is
/// correct at this scale; larger N would require more qubits than the simulator supports.
///
/// Qubit budget for N=91 (7-bit, t=7): 14 qubits total.
#[test]
fn qubit_budget_91_within_ceiling() {
    use shor::arith::ModExpLayout;
    let exp_len = exp_len_for(91);
    let layout = ModExpLayout::standard(91, exp_len);
    let total = layout.total_qubits();
    assert_eq!(exp_len, 7, "exp_len for N=91 should be 7");
    assert_eq!(layout.work_len, 7, "work register for N=91 should be 7 qubits");
    assert_eq!(total, 14, "total qubits for N=91 should be 14");
    assert!(
        total <= 25,
        "N=91 requires {total} qubits, within the ~25-qubit ceiling"
    );
}

/// Qubit budgets for all factoring targets are within the ~25-qubit ceiling.
///
/// Verifies the ancilla-free permutation-synthesis qubit budgets for all factoring targets.
#[test]
fn qubit_budgets_all_targets() {
    use shor::arith::ModExpLayout;
    // (N, expected_exp_len, expected_work_len, expected_total)
    // exp_len = n_bits(N), work_len = n_bits(N+1)
    let targets = [(15u64, 4, 4, 8), (21, 5, 5, 10), (35, 6, 6, 12), (91, 7, 7, 14)];
    for (n, exp_len, work_len, total) in targets {
        let layout = ModExpLayout::standard(n, exp_len_for(n));
        assert_eq!(layout.exp_len, exp_len, "N={n}: exp_len mismatch");
        assert_eq!(layout.work_len, work_len, "N={n}: work_len mismatch");
        assert_eq!(layout.total_qubits(), total, "N={n}: total qubits mismatch");
        assert!(
            layout.total_qubits() <= MAX_QUBITS,
            "N={n}: {total} qubits exceeds the {MAX_QUBITS}-qubit ceiling"
        );
    }
}

// ── Additional robustness KATs ────────────────────────────────────────────────

/// factor(N) is reproducible: same seed → same result, evidence included.
#[test]
fn factor_reproducible() {
    for (n, seed) in [(15u64, 0u64), (21, 0), (35, 7)] {
        let r1 = factor(n, seed);
        let r2 = factor(n, seed);
        assert_eq!(r1, r2, "factor({n}, seed={seed}) is not reproducible");
    }
}

/// Even N is reported as `Even` with p = 2, and no circuit runs.
#[test]
fn factor_even_short_circuit() {
    for (n, q) in [(14u64, 7u64), (22, 11), (4, 2), (4096, 2048)] {
        let found = factor(n, 0).unwrap_or_else(|e| panic!("factor({n}) failed: {e:?}"));
        assert_eq!((found.p, found.q), (2, q), "factor({n})");
        assert_eq!(found.evidence, ShorEvidence::Even, "factor({n})");
    }
    // Even beyond the simulator ceiling is still classified without a circuit.
    let found = factor(u64::MAX - 1, 0).expect("even n of any size");
    assert_eq!((found.p, found.evidence), (2, ShorEvidence::Even));
}

/// Prime N, and N ≤ 1, are `NotComposite`; the primes below the ceiling are rejected before any
/// circuit runs (a circuit for these would take seconds, so a pass here is fast by construction).
#[test]
fn factor_prime_is_not_composite() {
    for n in [0u64, 1, 2, 3, 5, 13, 97, 4093] {
        assert_eq!(factor(n, 0), Err(ShorFactorError::NotComposite), "factor({n})");
    }
}

/// An odd prime power is `PrimePower` with its base and exponent, without a circuit.
#[test]
fn factor_prime_power_short_circuit() {
    for (n, prime, exponent) in [(9u64, 3u64, 2u32), (27, 3, 3), (125, 5, 3), (3u64.pow(40), 3, 40)]
    {
        let found = factor(n, 0).unwrap_or_else(|e| panic!("factor({n}) failed: {e:?}"));
        assert_eq!((found.p, found.q), (prime, n / prime), "factor({n})");
        assert_eq!(found.evidence, ShorEvidence::PrimePower { prime, exponent }, "factor({n})");
        assert_eq!(prime.pow(exponent), n);
    }
}

/// A circuit beyond `MAX_QUBITS` is `Unsupported` and reports the qubits it would need.
///
/// 4097 = 17 × 241 needs 13 exponent + 13 work qubits = 26 > 25. Odd N that is not a prime power
/// is rejected before it is tested for primality, so the large prime 1000003 is `Unsupported` too.
#[test]
fn factor_beyond_the_simulator_is_unsupported() {
    assert_eq!(MAX_QUBITS, 25);
    assert_eq!(factor(4097, 0), Err(ShorFactorError::Unsupported { qubits: 26 }));
    assert_eq!(factor(1_000_003, 0), Err(ShorFactorError::Unsupported { qubits: 40 }));
    assert_eq!(factor(u64::MAX, 0), Err(ShorFactorError::Unsupported { qubits: 128 }));
}

/// find_order returns the least order: a^r ≡ 1 mod N and no smaller exponent does.
///
/// Each seed is one whose circuit run recovers an order (measured phase `s` in the comment).
#[test]
fn find_order_returns_valid_order() {
    // (a, N, seed, s): s is the phase numerator measured by the circuit with that seed.
    let cases = [
        (2u64, 15u64, 0u64, 4u64),
        (7, 15, 0, 4),
        (2, 21, 0, 27),
        (4, 35, 1, 11),
        (2, 35, 7, 59),
    ];
    for (a, n, seed, s) in cases {
        let r = find_order(a, n, seed)
            .unwrap_or_else(|| panic!("find_order({a}, {n}, seed={seed}) found no order"));
        assert_eq!(
            Some(r),
            order_from_phase(s, exp_len_for(n), a, n),
            "find_order({a}, {n}, seed={seed}) disagrees with its measured phase s={s}"
        );
        assert_eq!(mod_pow(a, r, n), 1, "find_order({a}, {n}): {a}^{r} mod {n} ≠ 1");
        assert_eq!(r, brute_force_order(a, n), "find_order({a}, {n}): {r} is not the least order");
    }
}

/// find_order is `None` for a measurement that carries no usable phase.
///
/// Seed 1 measures s=8 (phase 1/2) for a=2 mod 15: the only convergent denominator is 2, and
/// 2^2 = 4 ≢ 1.
#[test]
fn find_order_none_on_uninformative_measurement() {
    assert_eq!(find_order(2, 15, 1), None);
}


