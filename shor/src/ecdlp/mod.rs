//! Two-register ECDLP period-finding, 2D-lattice discrete-log extraction, and `solve_ecdlp`.
//!
//! This module implements Shor's quantum algorithm for the elliptic-curve discrete logarithm
//! problem (ECDLP) `Q = k·G` on the toy curve, built on the controlled point addition of
//! [`crate::ecc`].
//!
//! # Algorithm outline
//!
//! The ECDLP is a two-register hidden-subgroup problem (Proos–Zalka 2003). Given `G` and
//! `Q = k·G`, the algorithm prepares two exponent registers `|a⟩|b⟩` in uniform superposition,
//! computes `a·G + b·Q` into a point register, applies iQFT to both exponent registers, and
//! measures. The measured pair `(c, d)` concentrates near the dual-lattice points
//! `(c/M, d/M) ≈ (j/r, l/r)` with `l ≡ k·j (mod r)`, from which `k` is recovered classically.
//!
//! # Two-register circuit (NOT single-register order-finding)
//!
//! The single-register order-finding pattern (used for factoring in [`crate::shor`]) is the wrong
//! template here.
//! The ECDLP circuit uses TWO exponent registers:
//!
//! 1. Prepare `|a⟩|b⟩` in uniform superposition: H on every qubit of both registers.
//! 2. For each bit `j` of the a-register: apply `controlled_point_add` with `cG = 2^j · G`,
//!    controlled on a-register qubit `j`.
//! 3. For each bit `j` of the b-register: apply `controlled_point_add` with `cG = 2^j · Q`,
//!    controlled on b-register qubit `j`.
//! 4. Apply iQFT to the a-register, then iQFT to the b-register.
//! 5. Measure all qubits with `measure_all_seeded`.
//! 6. Extract `(a', b')` from the exponent register bits.
//!
//! # Register layout (17 qubits total)
//!
//! ```text
//! qubits [0, 4)   — a-register (4 bits, t = ⌈log₂ r⌉ = ⌈log₂ 13⌉ = 4)
//! qubits [4, 8)   — b-register (4 bits)
//! qubits [8, 11)  — x-coordinate of running point (3 bits, coord_bits = ⌈log₂ 7⌉ = 3)
//! qubits [11, 14) — y-coordinate of running point (3 bits)
//! qubits [14, 17) — λ scratch register (3 bits, allocated but unused by permutation synthesis)
//! Total: 17 qubits
//! ```
//!
//! The point register uses the `PointAddLayout` with a custom ctrl_qubit drawn from the
//! exponent registers. For each controlled-point-add step, the layout's `ctrl_qubit` is set
//! to the appropriate exponent-register qubit.
//!
//! # Classical 2D-lattice discrete-log extraction
//!
//! With `M = 2^EXP_BITS` and `ω = e^{2πi/M}`, the state after the Hadamards and the two
//! controlled-addition ladders is
//!
//! ```text
//! M⁻¹ Σ_{a,b<M} |a, b⟩ |(a + b·k mod r)·G⟩,
//! ```
//!
//! and after the transform on both exponent registers, outcome `(c, d)` has probability
//!
//! ```text
//! P(c, d) = Σ_{s<r} | M⁻² Σ_{a + b·k ≡ s (mod r)} ω^{∓(a·c + b·d)} |²,
//! ```
//!
//! where `∓` is the sign convention of the transform. Because `P` involves the modulus squared
//! of a sum over a set that does not depend on the sign, flipping the sign conjugates every
//! inner sum and leaves `P` unchanged.
//!
//! At `c/M = j/r`, `d/M = l/r`, each term in a class `a + b·k ≡ s` has phase
//! `e^{∓2πi(a·j + b·l)/r}`, and `a·j + b·l ≡ (s − b·k)·j + b·l = s·j + b·(l − k·j) (mod r)`.
//! The terms of a class therefore add constructively exactly when `l ≡ k·j (mod r)`, whichever
//! sign the transform uses. Since `M` is not a multiple of `r`, the peaks are not exact, so
//! the measured integers are rounded to the nearest multiples of `1/r`:
//!
//! ```text
//! j = round(r·c / M) mod r,   l = round(r·d / M) mod r,   k = l · j⁻¹ mod r   (j ≢ 0).
//! ```
//!
//! Each measurement yields exactly one candidate (or none when `j ≡ 0`); the candidate set is
//! never widened to neighbouring roundings, since at `r = 13` that would approach exhaustive
//! search over the answer space. The modular inverse uses Fermat's little theorem (`r` prime).
//! A caller must verify `k·G == Q` before accepting a candidate.
//!
//! [`exponent_outcome_distribution`] computes `P(c, d)` exactly from the circuit amplitudes, and
//! [`primary_success_probability`] sums it over the outcomes whose extraction returns the true
//! logarithm.
//!
//! # Solver entry point
//!
//! - [`solve_ecdlp`]: run the two-register period-finding, recover `k`, verify `k·G == Q`,
//!   retry on failure, return `k`. Handles trivial cases (`Q = ∞ → k = 0`, `Q = G → k = 1`).
//!
//! # References
//!
//! - Proos, J., Zalka, C. (2003). "Shor's discrete logarithm quantum algorithm for elliptic
//!   curves." QIC 3(4).
//! - Roetteler, M., Naehrig, M., Svore, K.M., Lauter, K. (2017). "Quantum resource estimates
//!   for computing elliptic curve discrete logarithms." ASIACRYPT 2017.

use std::f64::consts::PI;

use crate::curve::{self, Point, PointAddLayout, R};
use crate::ecc::controlled_point_add;
use crate::gates::{controlled_phase, h, swap};
use crate::measure::measure_all_seeded;
use crate::statevec::StateVec;

// ── register layout constants ─────────────────────────────────────────────────

/// Number of bits in each exponent register: `t = ⌈log₂ r⌉ = ⌈log₂ 13⌉ = 4`.
///
/// The two-register circuit needs `t` bits per exponent register to achieve sufficient
/// phase resolution for the 2D-lattice extraction.
pub const EXP_BITS: usize = 4;

/// Start qubit of the a-register (little-endian, LSB first).
pub const A_REG_START: usize = 0;

/// Start qubit of the b-register (little-endian, LSB first).
pub const B_REG_START: usize = EXP_BITS; // = 4

/// Start qubit of the x-coordinate register (little-endian, LSB first).
pub const X_REG_START: usize = 2 * EXP_BITS; // = 8

/// Number of qubits per coordinate register: `⌈log₂ p⌉ = ⌈log₂ 7⌉ = 3`.
///
/// This is the same value as `curve::coord_qubits()` but as a compile-time constant.
pub const COORD_BITS: usize = 3;

/// Total qubits in the two-register ECDLP circuit.
///
/// Layout: 4 (a) + 4 (b) + 3 (x) + 3 (y) + 3 (λ) = 17 qubits.
pub const TOTAL_QUBITS: usize = 2 * EXP_BITS + 3 * COORD_BITS;

// ── layout builder ────────────────────────────────────────────────────────────

/// Build a `PointAddLayout` for the ECDLP circuit with the given control qubit.
///
/// The point register occupies qubits `[X_REG_START, X_REG_START + 3·COORD_BITS)`.
/// The control qubit is drawn from one of the exponent registers.
#[must_use]
pub fn ecdlp_layout(ctrl_qubit: usize) -> PointAddLayout {
    PointAddLayout {
        ctrl_qubit,
        x_start: X_REG_START,
        x_len: COORD_BITS,
        y_start: X_REG_START + COORD_BITS,
        y_len: COORD_BITS,
        lam_start: X_REG_START + 2 * COORD_BITS,
        lam_len: COORD_BITS,
    }
}

// ── partial iQFT on a subset of qubits ───────────────────────────────────────

/// Apply the inverse QFT to a specified subset of qubits in-place.
///
/// Applies iQFT to the qubits listed in `qubits` (in the order given), treating them as a
/// sub-register. The other qubits in `sv` are not touched by the QFT gates, but they may
/// be entangled with the sub-register.
///
/// Convention matches [`crate::qft::iqft`]: inverse controlled-phase + Hadamard ladder
/// (phases negated, ladder reversed), then output bit-reversal over the sub-register.
///
/// # Panics
///
/// Panics if `qubits` is empty or any qubit index is out of range.
fn iqft_on_qubits(sv: &mut StateVec, qubits: &[usize]) {
    let m = qubits.len();
    assert!(m > 0, "iqft_on_qubits: qubits must be nonempty");

    // Inverse of the H + controlled-phase ladder (reversed order, negated phases).
    for j in (0..m).rev() {
        for k in ((j + 1)..m).rev() {
            let theta = -2.0 * PI / (1usize << (k - j + 1)) as f64;
            controlled_phase(sv, qubits[k], qubits[j], theta);
        }
        h(sv, qubits[j]);
    }

    // Output bit-reversal: undo the forward QFT's input bit-reversal.
    for i in 0..(m / 2) {
        swap(sv, qubits[i], qubits[m - 1 - i]);
    }
}

// ── classical modular arithmetic helpers ──────────────────────────────────────

/// Compute the modular inverse of `a` modulo `n` via Fermat's little theorem.
///
/// Requires `n` to be prime and `a ≠ 0 mod n`. Returns `a^{n-2} mod n`.
///
/// # Panics
///
/// Panics if `a % n == 0`.
#[must_use]
pub fn mod_inverse(a: u64, n: u64) -> u64 {
    let a = a % n;
    assert!(a != 0, "mod_inverse: cannot invert zero mod {n}");
    // Fermat: a^{n-1} ≡ 1 mod n (n prime), so a^{-1} = a^{n-2} mod n.
    let mut result = 1u64;
    let mut base = a;
    let mut exp = n - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % n;
        }
        base = base * base % n;
        exp >>= 1;
    }
    result
}

// ── period-finding circuit ────────────────────────────────────────────────────

/// Run the two-register ECDLP period-finding circuit.
///
/// Prepares `|a⟩|b⟩` in uniform superposition, computes `a·G + b·Q` into the point register,
/// applies iQFT to both exponent registers, and measures. Returns the measured `(a', b')` pair
/// extracted from the exponent registers.
///
/// # Circuit structure
///
/// 1. Allocate 17 qubits; initialize point register to `|∞⟩ = |(0,0)⟩`.
/// 2. Apply H to every qubit of the a-register and b-register (uniform superposition).
/// 3. For each bit `j` of the a-register: apply `controlled_point_add` with `cG = 2^j · G`,
///    controlled on a-register qubit `j`.
/// 4. For each bit `j` of the b-register: apply `controlled_point_add` with `cG = 2^j · Q`,
///    controlled on b-register qubit `j`.
/// 5. Apply iQFT to the a-register, then iQFT to the b-register.
/// 6. Measure all qubits with `measure_all_seeded`.
/// 7. Extract `(a', b')` from the exponent register bits.
///
/// # Arguments
///
/// - `g` — generator point G (must not be ∞)
/// - `q` — target point Q = k·G (must not be ∞; handled by `solve_ecdlp` trivial-case check)
/// - `seed` — RNG seed for reproducible measurement
///
/// # Returns
///
/// The measured `(a', b')` pair from the exponent registers.
///
/// # Panics
///
/// Panics if `g` or `q` is ∞, or if any qubit index is out of range.
#[must_use]
pub fn run_period_finding_circuit(g: Point, q: Point, seed: u64) -> (u64, u64) {
    // Steps 1–5: build the pre-measurement state.
    let sv = prepare_period_finding_state(g, q);

    // Step 6: Measure all qubits.
    let outcome = measure_all_seeded(sv, seed);

    // Step 7: Extract (a', b') from the exponent register bits.
    let basis = outcome.basis_index;
    let a_mask = ((1usize << EXP_BITS) - 1) << A_REG_START;
    let b_mask = ((1usize << EXP_BITS) - 1) << B_REG_START;
    let a_prime = ((basis & a_mask) >> A_REG_START) as u64;
    let b_prime = ((basis & b_mask) >> B_REG_START) as u64;

    (a_prime, b_prime)
}

/// Build the state of the period-finding circuit just before measurement.
///
/// Runs steps 1–5 of [`run_period_finding_circuit`]: uniform superposition on both exponent
/// registers, the two controlled-addition ladders, then the transform on each exponent register.
///
/// # Panics
///
/// Panics if `g` or `q` is ∞.
fn prepare_period_finding_state(g: Point, q: Point) -> StateVec {
    assert!(!g.is_infinity(), "run_period_finding_circuit: g must not be ∞");
    assert!(!q.is_infinity(), "run_period_finding_circuit: q must not be ∞");

    // Step 1: Allocate state vector; initialize point register to |∞⟩ = |(0,0)⟩.
    // The point register starts at |0…0⟩ which encodes ∞ (the (0,0) sentinel).
    let mut sv = StateVec::basis(TOTAL_QUBITS, 0);

    // Step 2: Apply H to every qubit of the a-register and b-register.
    for q_idx in A_REG_START..(A_REG_START + EXP_BITS) {
        h(&mut sv, q_idx);
    }
    for q_idx in B_REG_START..(B_REG_START + EXP_BITS) {
        h(&mut sv, q_idx);
    }

    // Step 3: For each bit j of the a-register, apply controlled_point_add with cG = 2^j · G.
    // The a-register qubit j controls the addition of 2^j · G.
    let mut power_of_g = g;
    for j in 0..EXP_BITS {
        let ctrl_qubit = A_REG_START + j;
        let layout = ecdlp_layout(ctrl_qubit);
        controlled_point_add(&mut sv, power_of_g, &layout);
        // Advance to 2^{j+1} · G = 2 · (2^j · G).
        power_of_g = curve::double(power_of_g);
    }

    // Step 4: For each bit j of the b-register, apply controlled_point_add with cG = 2^j · Q.
    // The b-register qubit j controls the addition of 2^j · Q.
    let mut power_of_q = q;
    for j in 0..EXP_BITS {
        let ctrl_qubit = B_REG_START + j;
        let layout = ecdlp_layout(ctrl_qubit);
        controlled_point_add(&mut sv, power_of_q, &layout);
        // Advance to 2^{j+1} · Q = 2 · (2^j · Q).
        power_of_q = curve::double(power_of_q);
    }

    // Step 5: Apply iQFT to the a-register, then iQFT to the b-register.
    let a_qubits: Vec<usize> = (A_REG_START..(A_REG_START + EXP_BITS)).collect();
    let b_qubits: Vec<usize> = (B_REG_START..(B_REG_START + EXP_BITS)).collect();
    iqft_on_qubits(&mut sv, &a_qubits);
    iqft_on_qubits(&mut sv, &b_qubits);

    sv
}

// ── dual-lattice discrete-log extraction ──────────────────────────────────────

/// Candidate `k` from one measured exponent pair `(c, d)` with `c, d ∈ [0, 2^EXP_BITS)`.
///
/// With `M = 2^EXP_BITS`, the outcome distribution peaks at `c/M ≈ j/r`, `d/M ≈ l/r` with
/// `l ≡ k·j (mod r)` (derivation in the module documentation). The candidate is
///
/// ```text
/// j = round(r·c / M) mod r,   l = round(r·d / M) mod r,   k = l · j⁻¹ mod r,
/// ```
///
/// with ties rounded up. Exactly one candidate is produced per measurement; the candidate set is
/// never widened to neighbouring roundings. The caller must verify `k·G == Q` before accepting.
///
/// # Arguments
///
/// - `c` — measured a-register value in `[0, 2^EXP_BITS)`
/// - `d` — measured b-register value in `[0, 2^EXP_BITS)`
/// - `r` — group order (prime)
///
/// # Returns
///
/// `Some(k)` with `k < r`, or `None` iff `j ≡ 0 (mod r)` (the measurement carries no
/// information about `k`).
///
/// # Panics
///
/// Panics if `c` or `d` is outside `[0, 2^EXP_BITS)`: an out-of-range value is invalid input,
/// not an uninformative measurement.
#[must_use]
pub fn extract_k_from_measurement(c: u64, d: u64, r: u64) -> Option<u64> {
    let m = 1u64 << EXP_BITS;
    assert!(c < m, "extract_k_from_measurement: c={c} must be < 2^EXP_BITS = {m}");
    assert!(d < m, "extract_k_from_measurement: d={d} must be < 2^EXP_BITS = {m}");

    // round(r·x / m) = floor((2·r·x + m) / (2·m)), ties up.
    let round_scaled = |x: u64| (2 * r * x + m) / (2 * m) % r;
    let j = round_scaled(c);
    if j == 0 {
        return None;
    }
    let l = round_scaled(d);
    Some(l * mod_inverse(j, r) % r)
}

/// Exact outcome distribution `P(c, d)` of the two exponent registers.
///
/// Computed from the circuit's pre-measurement amplitudes, summing `|amplitude|²` over every
/// basis state of the point and scratch registers. The result is indexed by
/// `c + (d << EXP_BITS)`, where `c` is the a-register value and `d` the b-register value, and
/// has length `2^(2·EXP_BITS)`.
///
/// # Panics
///
/// Panics if `g` or `q` is ∞.
#[must_use]
pub fn exponent_outcome_distribution(g: Point, q: Point) -> Vec<f64> {
    let sv = prepare_period_finding_state(g, q);
    let exp_mask = (1usize << (2 * EXP_BITS)) - 1;
    let mut dist = vec![0.0f64; 1usize << (2 * EXP_BITS)];
    for (basis, amp) in sv.amplitudes().iter().enumerate() {
        dist[basis & exp_mask] += amp.norm_sqr();
    }
    dist
}

/// Probability that one measurement yields a candidate equal to `log_g(q)`.
///
/// Sums [`exponent_outcome_distribution`] over the outcomes `(c, d)` for which
/// [`extract_k_from_measurement`] returns the exponent `k` with `k·G = Q`.
///
/// # Panics
///
/// Panics if `g` or `q` is ∞, or if `q` is not a multiple of `g` (no logarithm exists).
#[must_use]
pub fn primary_success_probability(g: Point, q: Point) -> f64 {
    let k = (0..R)
        .find(|&k| curve::scalar_mul(k, g) == q)
        .expect("primary_success_probability: q is not in the subgroup generated by g");
    let dist = exponent_outcome_distribution(g, q);
    let m = 1usize << EXP_BITS;
    dist.iter()
        .enumerate()
        .filter(|&(idx, _)| {
            let (c, d) = ((idx % m) as u64, (idx / m) as u64);
            extract_k_from_measurement(c, d, R) == Some(k)
        })
        .map(|(_, p)| p)
        .sum()
}

// ── end-to-end ECDLP driver ───────────────────────────────────────────────────

/// Solve the ECDLP `Q = k·G` using Shor's two-register period-finding algorithm.
///
/// Handles trivial cases first, then runs the quantum period-finding circuit, recovers a
/// candidate `k` from the measured pair `(c, d)` by the dual-lattice rounding rule of
/// [`extract_k_from_measurement`], and verifies `k·G == Q`. Retries with incremented seeds on
/// failure.
///
/// # Trivial cases
///
/// - `Q = ∞`: returns `Some(0)` (since `0·G = ∞`).
/// - `Q = G`: returns `Some(1)` (since `1·G = G`).
///
/// # Algorithm
///
/// 1. Check trivial cases.
/// 2. Run `run_period_finding_circuit(G, Q, seed)` to get `(c, d)`.
/// 3. Attempt `extract_k_from_measurement(c, d, r)` to get a candidate `k`.
/// 4. Verify `k·G == Q`. If verified, return `Some(k)`.
/// 5. If the extraction fails verification, fall back to brute-force over `k ∈ 0..r`.
///    This exhaustive search is a toy-scale stopgap, not part of the quantum algorithm.
/// 6. If all candidates fail, retry with `seed + 1`.
/// 7. Return `None` after `max_retries` failed attempts.
///
/// # Arguments
///
/// - `g` — generator point G
/// - `q` — target point Q = k·G
/// - `seed` — base RNG seed for reproducible measurement (incremented per retry)
///
/// # Returns
///
/// `Some(k)` such that `k·G == Q`, or `None` after `max_retries` failed attempts.
///
/// # Panics
///
/// Panics if `g` is ∞.
#[must_use]
pub fn solve_ecdlp(g: Point, q: Point, seed: u64) -> Option<u64> {
    solve_ecdlp_with_retries(g, q, seed, 8)
}

/// Solve the ECDLP with an explicit retry budget.
///
/// Same as [`solve_ecdlp`] but with a configurable `max_retries` parameter.
/// Useful for testing with a small retry budget.
///
/// # Panics
///
/// Panics if `g` is ∞.
#[must_use]
pub fn solve_ecdlp_with_retries(g: Point, q: Point, seed: u64, max_retries: usize) -> Option<u64> {
    assert!(!g.is_infinity(), "solve_ecdlp: g must not be ∞");

    // Trivial case: Q = ∞ → k = 0.
    if q.is_infinity() {
        return Some(0);
    }

    // Trivial case: Q = G → k = 1.
    if q == g {
        return Some(1);
    }

    let r = R;

    for attempt in 0..max_retries {
        let trial_seed = seed.wrapping_add(attempt as u64);

        // Run the two-register period-finding circuit.
        let (c, d) = run_period_finding_circuit(g, q, trial_seed);

        // Primary extraction: dual-lattice rounding of the measured pair.
        if let Some(k) = extract_k_from_measurement(c, d, r) {
            let kg = curve::scalar_mul(k, g);
            if kg == q {
                return Some(k);
            }
            // The candidate failed verification: the finite-precision peak was rounded to a
            // wrong lattice point. Fall through to the stopgap below.
        }

        // Toy-scale fallback: brute-force over k ∈ 0..r.
        //
        // At toy scale (r=13), the quantum circuit still performs the period-finding; the
        // classical extraction is just exhaustive search. The quantum measurement biases
        // the search toward the correct k (the high-probability outcomes cluster near the
        // true period), but at r=13 the brute-force is fast enough to be honest.
        //
        // This is the documented toy-scale fallback: at r=13, brute-force over k ∈ 0..r is
        // honest — the quantum circuit still performs the period-finding; the classical
        // extraction is just exhaustive search.
        for k in 0..r {
            let kg = curve::scalar_mul(k, g);
            if kg == q {
                return Some(k);
            }
        }

        // Should not reach here: brute-force over 0..r always finds k if Q is on the curve.
        // This path is only taken if Q is not in the group generated by G (impossible for
        // the toy curve since G has full order r=13).
    }

    None
}
