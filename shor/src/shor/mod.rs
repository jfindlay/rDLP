//! Order-finding circuit, continued-fraction period extraction, and `factor(N)` driver.
//!
//! This module implements the order-finding circuit orchestration, the classical
//! continued-fraction post-processing, and the end-to-end `factor(N)` driver.
//!
//! # Algorithm outline
//!
//! 1. **Order-finding circuit** ([`find_order`]): prepare the exponent register in uniform
//!    superposition, apply the C-ModExp circuit (`|x⟩|1⟩ → |x⟩|aˣ mod N⟩`), apply iQFT to
//!    the exponent register, measure to get a phase numerator `s ≈ k·2^t/r`.
//!
//! 2. **Continued-fraction period extraction** ([`order_from_phase`]): expand `s/2^t` as a
//!    continued fraction, take convergents with denominator `≤ N`, and recover the order `r`
//!    (the least `r > 0` with `a^r ≡ 1 mod N`: the first convergent denominator `d` with
//!    `a^d ≡ 1` is reduced by dividing out prime factors while the identity holds).
//!
//! 3. **`factor(N)` driver** ([`factor`]): classify `N` classically (not composite, even, prime
//!    power, beyond the simulator), else pick `a` coprime to `N`, run order-finding, extract
//!    `r`; if `r` is even and `a^(r/2) ≢ -1 mod N`, return `gcd(a^(r/2) ± 1, N)` as a
//!    [`ShorFactorization`] with its [`ShorEvidence`]; else retry within a fixed budget.
//!
//! # Register layout (C-ModExp — consumed)
//!
//! The order-finding circuit uses the standard `ModExpLayout::standard(N, exp_len)` layout:
//! - Qubits `[0, exp_len)` — exponent register (t bits, put in superposition)
//! - Qubits `[exp_len, exp_len + n)` — work register (n bits, initialized to `|1⟩`)
//!
//! Total qubits: `exp_len + n`, where `N < 2^n` and `exp_len = n_bits(N)`.
//!
//! The textbook analysis takes `2^t ≥ N²` for continued-fraction recovery; `t = n_bits(N)` is
//! half that width, and its sufficiency here is empirical (see [`factor`]).
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
//! All within the 25-qubit simulator ceiling ([`MAX_QUBITS`]).
//!
//! # Public surface
//!
//! - [`find_order`] + [`order_from_phase`]: the order-finding and period-extraction interface.
//! - [`factor`]: the end-to-end factoring driver, returning a [`ShorFactorization`] with its
//!   [`ShorEvidence`] or a [`ShorFactorError`].
//!
//! # References
//!
//! - Shor, P.W. (1994). "Algorithms for quantum computation: discrete logarithms and factoring."
//!   FOCS 1994.
//! - Nielsen, M.A., Chuang, I.L. (2000). "Quantum Computation and Quantum Information."
//!   Cambridge University Press. §5.3.

use std::f64::consts::PI;

use crate::arith::{controlled_mod_exp, mod_pow, n_bits, ModExpLayout};
use crate::gates::{controlled_phase, h, swap};
use crate::measure::measure_all_seeded;
use crate::statevec::StateVec;

// ── classical helpers ─────────────────────────────────────────────────────────

/// Compute `gcd(a, b)` using the Euclidean algorithm.
///
/// Returns the greatest common divisor of `a` and `b`. Returns 0 if both inputs are 0.
#[must_use]
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Compute the exponent register size `t = n_bits(N)` for order-finding.
///
/// Uses `t = n_bits(N) = ⌈log₂(N)⌉` bits for the exponent register. This matches the
/// ancilla-free permutation-synthesis qubit budget:
/// - N=15: t=4, total=8 qubits
/// - N=21: t=5, total=10 qubits
/// - N=35: t=6, total=12 qubits
/// - N=91: t=7, total=14 qubits
///
/// Phase resolution: `1/2^t`. The textbook bound is `2^t ≥ N²`; this width is half of it, so
/// sufficiency is empirical: continued-fraction recovery of the order `r` works for the four
/// factoring targets, and no bound is proved.
///
/// # Panics
///
/// Panics if `n == 0`.
#[must_use]
pub fn exp_len_for(n: u64) -> usize {
    n_bits(n)
}

/// Compute the continued-fraction convergents of `p/q` up to denominator bound `max_denom`.
///
/// Returns a `Vec` of `(numerator, denominator)` pairs for each convergent, in order from
/// the first (coarsest) to the last (finest) convergent with denominator `≤ max_denom`.
///
/// The continued-fraction expansion of `p/q` is computed via the Euclidean algorithm.
/// Each convergent `h_k/k_k` satisfies `|p/q - h_k/k_k| ≤ 1/(k_k · k_{k+1})`.
///
/// # Panics
///
/// Panics if `q == 0`.
#[must_use]
pub fn convergents(p: u64, q: u64, max_denom: u64) -> Vec<(u64, u64)> {
    assert!(q > 0, "convergents: denominator q must be nonzero");
    let mut result = Vec::new();

    // Continued-fraction coefficients via the Euclidean algorithm.
    // p/q = a_0 + 1/(a_1 + 1/(a_2 + ...))
    // Convergents: h_{-1}=1, h_0=a_0, h_k = a_k*h_{k-1} + h_{k-2}
    //              k_{-1}=0, k_0=1,   k_k = a_k*k_{k-1} + k_{k-2}
    // Standard continued-fraction recurrence:
    //   h_k = a_k * h_{k-1} + h_{k-2},  with h_{-2}=0, h_{-1}=1
    //   k_k = a_k * k_{k-1} + k_{k-2},  with k_{-2}=1, k_{-1}=0
    let mut h_prev = 0u64; // h_{k-2}, starts at h_{-2}=0
    let mut h_curr = 1u64; // h_{k-1}, starts at h_{-1}=1
    let mut k_prev = 1u64; // k_{k-2}, starts at k_{-2}=1
    let mut k_curr = 0u64; // k_{k-1}, starts at k_{-1}=0

    let mut num = p;
    let mut den = q;

    loop {
        let a = num / den;
        let rem = num % den;

        // h_k = a * h_{k-1} + h_{k-2}
        // k_k = a * k_{k-1} + k_{k-2}
        let h_next = a.saturating_mul(h_curr).saturating_add(h_prev);
        let k_next = a.saturating_mul(k_curr).saturating_add(k_prev);

        if k_next > max_denom {
            break;
        }

        result.push((h_next, k_next));

        h_prev = h_curr;
        h_curr = h_next;
        k_prev = k_curr;
        k_curr = k_next;

        if rem == 0 {
            break;
        }
        num = den;
        den = rem;
    }

    result
}

/// Reduce `d`, a multiple of the multiplicative order of `a` mod `n`, to that order.
///
/// Divides out each prime factor `f` of `d` for as long as `a^(d/f) ≡ 1 mod n` still holds. The
/// order divides `d`, so no prime power of the order is lost by a failed division.
///
/// Requires `a^d ≡ 1 mod n`, `d ≥ 1` and `n < 2^32` (the products inside [`mod_pow`] must fit
/// in a `u64`).
fn reduce_to_order(a: u64, d: u64, n: u64) -> u64 {
    let mut order = d;
    let mut rest = d;
    let mut f = 2u64;
    while f * f <= rest {
        if rest % f == 0 {
            while rest % f == 0 {
                rest /= f;
            }
            while order % f == 0 && mod_pow(a, order / f, n) == 1 {
                order /= f;
            }
        }
        f += 1;
    }
    // A remaining cofactor > 1 is itself prime.
    if rest > 1 {
        while order % rest == 0 && mod_pow(a, order / rest, n) == 1 {
            order /= rest;
        }
    }
    order
}

/// Recover the multiplicative order `r` of `a` mod `N` from a measured phase numerator `s`.
///
/// The phase is `s / 2^t ≈ k/r` for some integer `k`. Expands `s / 2^t` as a continued
/// fraction and takes the first convergent denominator `d ≤ N` with `a^d ≡ 1 mod N`. A
/// denominator can be a multiple of the order (`k/r` need not be in lowest terms, and the
/// convergent of a coarse phase can be a multiple of `r`), so `d` is then reduced: each prime
/// factor of `d` is divided out for as long as `a^(d/f) ≡ 1 mod N` still holds. The result is
/// the least `r > 0` with `a^r ≡ 1 mod N`, never a multiple of it.
///
/// Returns `Some(r)` if some convergent denominator satisfies `a^d ≡ 1 mod N`, `None`
/// otherwise (including `s = 0`, whose phase 0 carries no period information).
///
/// # Arguments
///
/// - `s` — measured phase numerator (from the exponent register measurement)
/// - `t` — exponent register size (number of qubits)
/// - `a` — the base whose order we seek
/// - `n` — the modulus (below `2^32`, so that modular products fit in a `u64`)
#[must_use]
pub fn order_from_phase(s: u64, t: usize, a: u64, n: u64) -> Option<u64> {
    if s == 0 {
        return None; // s=0 gives phase 0, which is uninformative
    }
    let two_t = 1u64 << t; // 2^t
    convergents(s, two_t, n)
        .into_iter()
        .map(|(_num, denom)| denom)
        .find(|&denom| denom > 0 && mod_pow(a, denom, n) == 1)
        .map(|denom| reduce_to_order(a, denom, n))
}

// ── partial iQFT on a subset of qubits ───────────────────────────────────────

/// Apply the inverse QFT to a specified subset of qubits in-place.
///
/// Applies iQFT to the qubits listed in `qubits` (in the order given), treating them as a
/// sub-register. The other qubits in `sv` are not touched by the QFT gates, but they may
/// be entangled with the sub-register.
///
/// The convention matches [`crate::qft::iqft`]: the inverse controlled-phase + Hadamard
/// ladder (phases negated, ladder reversed), then an output bit-reversal over the sub-register.
///
/// # Arguments
///
/// - `sv` — state vector (modified in-place)
/// - `qubits` — qubit indices of the sub-register (little-endian, LSB first)
///
/// # Panics
///
/// Panics if `qubits` is empty or any qubit index is out of range.
fn iqft_on_qubits(sv: &mut StateVec, qubits: &[usize]) {
    let m = qubits.len();
    assert!(m > 0, "iqft_on_qubits: qubits must be nonempty");

    // Inverse of the H + controlled-phase ladder (reversed order, negated phases).
    // The ladder processes qubits[0], qubits[1], ..., qubits[m-1] in the forward QFT.
    // The inverse reverses this: process qubits[m-1] down to qubits[0].
    for j in (0..m).rev() {
        for k in ((j + 1)..m).rev() {
            let theta = -2.0 * PI / (1usize << (k - j + 1)) as f64;
            controlled_phase(sv, qubits[k], qubits[j], theta);
        }
        h(sv, qubits[j]);
    }

    // Output bit-reversal: swap qubits[i] with qubits[m-1-i].
    // This undoes the input bit-reversal of the forward QFT.
    for i in 0..(m / 2) {
        swap(sv, qubits[i], qubits[m - 1 - i]);
    }
}

// ── order-finding circuit ─────────────────────────────────────────────────────

/// Run the order-finding circuit for base `a` modulo `N`, returning the measured phase numerator.
///
/// Circuit steps:
/// 1. Allocate `exp_len + n` qubits; initialize work register to `|1⟩`.
/// 2. Apply `H` to every exponent qubit (uniform superposition).
/// 3. Apply `controlled_mod_exp(sv, a, &layout)` — maps `|x⟩|1⟩ → |x⟩|aˣ mod N⟩`.
/// 4. Apply iQFT to the exponent register.
/// 5. Measure with `measure_all_seeded(sv, seed)`.
/// 6. Extract the exponent register value: `s = basis_index & ((1 << exp_len) - 1)`.
///
/// Returns the measured phase numerator `s` (the lower `exp_len` bits of the measurement).
///
/// # Arguments
///
/// - `a` — the base (must satisfy `gcd(a, N) = 1`)
/// - `layout` — register layout (from `ModExpLayout::standard`)
/// - `seed` — RNG seed for reproducible measurement
///
/// # Panics
///
/// Panics if `gcd(a, N) != 1` or any qubit index is out of range.
#[must_use]
pub fn run_order_finding_circuit(a: u64, layout: &ModExpLayout, seed: u64) -> u64 {
    let exp_len = layout.exp_len;
    let total_qubits = layout.total_qubits();

    // Step 1: Allocate state vector; initialize work register to |1⟩.
    // Work register starts at qubit exp_len (little-endian LSB = qubit exp_len).
    // Basis index for |1⟩ in the work register = 1 << exp_len (bit exp_len set = LSB of work).
    let work_init_index = 1usize << layout.work_start;
    let mut sv = StateVec::basis(total_qubits, work_init_index);

    // Step 2: Apply H to every exponent qubit (uniform superposition over exponent register).
    for q in layout.exp_start..(layout.exp_start + exp_len) {
        h(&mut sv, q);
    }

    // Step 3: Apply controlled modular exponentiation: |x⟩|1⟩ → |x⟩|aˣ mod N⟩.
    controlled_mod_exp(&mut sv, a, layout);

    // Step 4: Apply iQFT to the exponent register.
    let exp_qubits: Vec<usize> = (layout.exp_start..(layout.exp_start + exp_len)).collect();
    iqft_on_qubits(&mut sv, &exp_qubits);

    // Step 5: Measure the full register.
    let outcome = measure_all_seeded(sv, seed);

    // Step 6: Extract the exponent register value (lower exp_len bits of basis_index).
    let exp_mask = (1usize << exp_len) - 1;
    let s = (outcome.basis_index >> layout.exp_start) & exp_mask;
    s as u64
}

/// Find the order of `a` modulo `N` using the quantum order-finding circuit.
///
/// Runs the order-finding circuit with the given seed, then applies continued-fraction
/// period extraction to recover the order `r` (the smallest `r > 0` such that `a^r ≡ 1 mod N`).
///
/// Returns `Some(r)` (the least such `r`, see [`order_from_phase`]) if the circuit measurement
/// yields a phase that recovers the order, `None` if the measurement was uninformative (e.g.,
/// `s = 0` or no convergent satisfies `a^d ≡ 1 mod N`).
///
/// # Arguments
///
/// - `a` — the base (must satisfy `gcd(a, N) = 1`)
/// - `n` — the modulus
/// - `seed` — RNG seed for the measurement
///
/// # Panics
///
/// Panics if `gcd(a, N) != 1`, `n == 0`, or `a == 0`.
#[must_use]
pub fn find_order(a: u64, n: u64, seed: u64) -> Option<u64> {
    assert!(n > 1, "find_order: modulus n must be > 1");
    assert!(a > 0, "find_order: a must be nonzero");
    assert!(gcd(a % n, n) == 1, "find_order: gcd(a, N) must be 1, got a={a}, N={n}");

    let exp_len = exp_len_for(n);
    let layout = ModExpLayout::standard(n, exp_len);

    let s = run_order_finding_circuit(a, &layout, seed);
    order_from_phase(s, exp_len, a, n)
}

// ── end-to-end factoring driver ───────────────────────────────────────────────

/// Check whether `n` is a prime power `p^k` for `k >= 2`.
///
/// Returns `Some((p, k))` if `n = p^k` for some prime `p` and `k >= 2`, else `None`.
/// Used as a classical short-circuit in the factoring driver.
fn prime_power_base(n: u64) -> Option<(u64, u32)> {
    if n < 4 {
        return None;
    }
    // Try each possible exponent k from 2 up to log2(n).
    let max_k = 63 - n.leading_zeros(); // floor(log2(n))
    for k in 2..=max_k {
        // Compute the k-th root of n (integer).
        let root = integer_kth_root(n, u64::from(k));
        if root >= 2 && root.pow(k) == n {
            // Verify root is prime (simple trial division at toy scale).
            if is_prime_trial(root) {
                return Some((root, k));
            }
        }
    }
    None
}

/// Compute the integer k-th root of `n` (floor(n^{1/k})).
fn integer_kth_root(n: u64, k: u64) -> u64 {
    if k == 1 {
        return n;
    }
    if n == 0 {
        return 0;
    }
    // Newton's method for integer k-th root.
    let mut x = (n as f64).powf(1.0 / k as f64) as u64 + 1;
    loop {
        // x_new = ((k-1)*x + n/x^{k-1}) / k
        let xk1 = x.saturating_pow((k - 1) as u32);
        if xk1 == 0 {
            break;
        }
        let x_new = ((k - 1) * x + n / xk1) / k;
        if x_new >= x {
            break;
        }
        x = x_new;
    }
    // Adjust down if overshoot.
    while x > 0 && x.saturating_pow(k as u32) > n {
        x -= 1;
    }
    x
}

/// Simple primality test by trial division (sufficient for toy-scale N).
fn is_prime_trial(n: u64) -> bool {
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
    while i <= n / i {
        if n % i == 0 {
            return false;
        }
        i += 2;
    }
    true
}

/// Largest circuit, in qubits, that [`factor`] will simulate.
///
/// The dense register holds `2^n` complex amplitudes; at 25 qubits that is about 512 MiB (see
/// [`crate::statevec`]). With `exp_len = n_bits(N)` and the matching work register, this admits
/// `N ≤ 4096`.
pub const MAX_QUBITS: usize = 25;

/// Number of distinct bases `a` (coprime to `N`, in increasing order) that [`factor`] tries.
pub const MAX_BASES: usize = 20;

/// Number of circuit runs, each with a fresh seed, that [`factor`] makes per base.
pub const SEEDS_PER_BASE: usize = 8;

/// A verified nontrivial factorization `p · q = n` of the integer given to [`factor`].
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive] // readable everywhere; constructible only inside this crate
pub struct ShorFactorization {
    /// A factor with `1 < p < n`.
    pub p: u64,
    /// The cofactor `n / p`, with `p · q == n`.
    pub q: u64,
    /// How the factorization was obtained.
    pub evidence: ShorEvidence,
}

/// The provenance of an accepted factorization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShorEvidence {
    /// `n` is even, so `p = 2`; no circuit ran.
    Even,
    /// `n = prime^exponent` with `exponent ≥ 2`, so `p = prime`; no circuit ran.
    PrimePower {
        /// The prime base.
        prime: u64,
        /// The exponent, at least 2.
        exponent: u32,
    },
    /// One circuit run measured `s` for base `a`, from which the order `r` was recovered.
    ///
    /// Recheckable from the fields alone: `a^r ≡ 1 mod n` with `r` the least such exponent
    /// (`r = order_from_phase(s, exp_len_for(n), a, n)`), `s` is the measurement of
    /// `run_order_finding_circuit(a, &ModExpLayout::standard(n, exp_len_for(n)), seed)`, `r` is
    /// even, `h = a^(r/2) mod n` is not `n − 1`, and `p` is `gcd(h − 1, n)` or `gcd(h + 1, n)`.
    OrderFinding {
        /// The base, coprime to `n`.
        a: u64,
        /// Seed of the circuit run that measured `s`.
        seed: u64,
        /// Measured phase numerator.
        s: u64,
        /// The multiplicative order of `a` modulo `n`.
        r: u64,
    },
}

/// Why [`factor`] returned no factorization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShorFactorError {
    /// `n` has no nontrivial factorization: `n ≤ 1` (neither prime nor composite) or `n` is prime.
    NotComposite,
    /// The order-finding circuit for `n` needs more than [`MAX_QUBITS`] qubits.
    Unsupported {
        /// Qubits the circuit would need: the exponent register plus the work register.
        qubits: usize,
    },
    /// Every circuit run measured a phase that gave no even, non-trivial order.
    ///
    /// A resource outcome: it says nothing about `n` beyond the attempts made.
    AttemptsExhausted {
        /// Number of circuit runs made.
        attempts: usize,
    },
}

/// Qubits the order-finding circuit for `n` needs: `n_bits(n)` exponent plus `n_bits(n + 1)` work.
///
/// Equals `ModExpLayout::standard(n, exp_len_for(n)).total_qubits()`, computed without
/// constructing the layout so that `n = u64::MAX` does not overflow `n + 1`.
fn circuit_qubits(n: u64) -> usize {
    n_bits(n) + n_bits(n.saturating_add(1))
}

/// Accept `p` as a factor of `n` only if it is a nontrivial divisor, checked by division.
fn certify(n: u64, p: u64, evidence: ShorEvidence) -> Option<ShorFactorization> {
    (p > 1 && p < n && n % p == 0).then(|| ShorFactorization { p, q: n / p, evidence })
}

/// The factor Shor's extraction takes from a base `a` of even order `r`, if it gives one.
///
/// With `h = a^(r/2) mod n`: none if `r` is odd or `h ≡ −1`; otherwise the first nontrivial of
/// `gcd(h − 1, n)` and `gcd(h + 1, n)`.
fn factor_from_order(a: u64, r: u64, n: u64) -> Option<u64> {
    if r % 2 != 0 {
        return None;
    }
    let half_power = mod_pow(a, r / 2, n);
    if half_power == n - 1 {
        return None; // trivial square root of unity: no factor
    }
    [gcd(half_power.saturating_sub(1), n), gcd(half_power + 1, n)]
        .into_iter()
        .find(|&g| g > 1 && g < n)
}

/// Factor `N` using Shor's algorithm, returning the factorization with its provenance.
///
/// The returned [`ShorFactorization`] satisfies `p · q == n` with `1 < p, q`, checked by
/// division before it is returned.
///
/// # Algorithm
///
/// 1. **Classification** (no circuit): `n ≤ 1` and primes are [`ShorFactorError::NotComposite`];
///    even `n` is [`ShorEvidence::Even`]; a prime power `p^k` is [`ShorEvidence::PrimePower`]
///    (cheap enough to recognise at any size, so it is reported even beyond the simulator).
/// 2. **Ceiling**: a circuit needing more than [`MAX_QUBITS`] qubits (`N > 4096`) is
///    [`ShorFactorError::Unsupported`]. This is checked before the primality test, so an
///    unsupported odd `n` that is not a prime power is not classified further.
/// 3. **Order finding**: for each of the first [`MAX_BASES`] bases `a` coprime to `N`, make up to
///    [`SEEDS_PER_BASE`] circuit runs with seeds `seed + 1000·(base index) + (run index)`. Each
///    measurement `s` yields the order `r` of `a` through [`order_from_phase`]. If `r` is even
///    and `a^(r/2) ≢ −1 mod N`, then `gcd(a^(r/2) ∓ 1, N)` is a nontrivial factor, and the first
///    one found is returned as [`ShorEvidence::OrderFinding`].
///
/// Bases are drawn coprime to `N`, so a base sharing a factor with `N` is never tried and there
/// is no common-factor path: the classical `gcd(a, N)` shortcut would turn the base scan into
/// trial division, which is not the algorithm.
///
/// # Register width
///
/// The exponent register has `t = n_bits(N)` qubits. The textbook analysis takes `2^t ≥ N²`
/// (twice as many qubits) to bound the failure probability of continued-fraction recovery;
/// `t = ⌈log₂ N⌉` is half that. Its sufficiency here is empirical: it recovers the order for the
/// toy targets 15, 21, 35 and 91, and no bound is proved.
///
/// # Errors
///
/// [`ShorFactorError::NotComposite`] for `n ≤ 1` or prime `n`;
/// [`ShorFactorError::Unsupported`] if the circuit would exceed [`MAX_QUBITS`];
/// [`ShorFactorError::AttemptsExhausted`] if no circuit run gave a factor.
pub fn factor(n: u64, seed: u64) -> Result<ShorFactorization, ShorFactorError> {
    factor_with_budget(n, seed, MAX_BASES, SEEDS_PER_BASE)
}

/// [`factor`] with an explicit budget of `max_bases` bases and `runs_per_base` runs per base.
fn factor_with_budget(
    n: u64,
    seed: u64,
    max_bases: usize,
    runs_per_base: usize,
) -> Result<ShorFactorization, ShorFactorError> {
    if n <= 2 {
        return Err(ShorFactorError::NotComposite); // 0 and 1 have no factorization; 2 is prime
    }

    if n % 2 == 0 {
        if let Some(found) = certify(n, 2, ShorEvidence::Even) {
            return Ok(found);
        }
    }

    if let Some((prime, exponent)) = prime_power_base(n) {
        if let Some(found) = certify(n, prime, ShorEvidence::PrimePower { prime, exponent }) {
            return Ok(found);
        }
    }

    let qubits = circuit_qubits(n);
    if qubits > MAX_QUBITS {
        return Err(ShorFactorError::Unsupported { qubits });
    }

    if is_prime_trial(n) {
        return Err(ShorFactorError::NotComposite);
    }

    let exp_len = exp_len_for(n);
    let layout = ModExpLayout::standard(n, exp_len);
    let bases = (2u64..n).filter(|&a| gcd(a, n) == 1).take(max_bases);

    let mut attempts = 0usize;
    for (base_idx, a) in bases.enumerate() {
        for run_idx in 0..runs_per_base {
            let trial_seed = seed
                .wrapping_add(base_idx as u64 * 1000)
                .wrapping_add(run_idx as u64);

            attempts += 1;
            let s = run_order_finding_circuit(a, &layout, trial_seed);
            let Some(r) = order_from_phase(s, exp_len, a, n) else {
                continue;
            };
            let Some(p) = factor_from_order(a, r, n) else {
                continue;
            };
            let evidence = ShorEvidence::OrderFinding { a, seed: trial_seed, s, r };
            if let Some(found) = certify(n, p, evidence) {
                return Ok(found);
            }
        }
    }

    Err(ShorFactorError::AttemptsExhausted { attempts })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Least `r > 0` with `a^r ≡ 1 mod n`, by exhaustive scan (`gcd(a, n) = 1`).
    fn brute_force_order(a: u64, n: u64) -> u64 {
        (1..n).find(|&r| mod_pow(a, r, n) == 1).expect("a is a unit mod n")
    }

    #[test]
    fn reduce_to_order_finds_least_exponent_from_any_multiple() {
        for n in 3u64..120 {
            for a in (2..n).filter(|&a| gcd(a, n) == 1) {
                let order = brute_force_order(a, n);
                for m in 1..=12 {
                    assert_eq!(
                        reduce_to_order(a, order * m, n),
                        order,
                        "a={a}, n={n}, d={}",
                        order * m
                    );
                }
            }
        }
    }

    #[test]
    fn circuit_qubits_matches_standard_layout() {
        for n in 2u64..5000 {
            let layout = ModExpLayout::standard(n, exp_len_for(n));
            assert_eq!(circuit_qubits(n), layout.total_qubits(), "n={n}");
        }
        // No overflow of n + 1 at the top of the range.
        assert_eq!(circuit_qubits(u64::MAX), 128);
    }

    #[test]
    fn factor_from_order_rejects_odd_order_and_minus_one() {
        // ord₂(15) = 4, 2² = 4 ≢ -1: gcd(3, 15) = 3.
        assert_eq!(factor_from_order(2, 4, 15), Some(3));
        // ord₄(15) = 2, 4¹ = 4 ≢ -1: gcd(3, 15) = 3.
        assert_eq!(factor_from_order(4, 2, 15), Some(3));
        // ord₁₄(15) = 2, 14 ≡ -1: trivial square root.
        assert_eq!(factor_from_order(14, 2, 15), None);
        // ord₂(7) = 3 is odd.
        assert_eq!(factor_from_order(2, 3, 7), None);
    }

    #[test]
    fn spent_budget_is_attempts_exhausted_with_the_runs_made() {
        // One run with seed 1 measures s = 8 for a = 2 mod 15 (phase 1/2), which gives no order.
        assert_eq!(
            factor_with_budget(15, 1, 1, 1),
            Err(ShorFactorError::AttemptsExhausted { attempts: 1 })
        );
        // Seed 2 measures s = 12 (phase 3/4): the second run succeeds, so the budget of two does.
        let found = factor_with_budget(15, 1, 1, 2).expect("second run recovers the order");
        assert!(matches!(
            found.evidence,
            ShorEvidence::OrderFinding { a: 2, seed: 2, s: 12, r: 4 }
        ));
        // A zero budget makes no run.
        assert_eq!(
            factor_with_budget(15, 1, 0, 8),
            Err(ShorFactorError::AttemptsExhausted { attempts: 0 })
        );
    }

    #[test]
    fn prime_power_base_reports_prime_and_exponent() {
        assert_eq!(prime_power_base(9), Some((3, 2)));
        assert_eq!(prime_power_base(81), Some((3, 4)));
        assert_eq!(prime_power_base(3u64.pow(40)), Some((3, 40)));
        assert_eq!(prime_power_base(15), None);
        assert_eq!(prime_power_base(u64::MAX), None);
    }
}
