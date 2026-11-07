//! Born-rule measurement: sample a basis state with probability |aᵢ|², collapse the register.
//!
//! # Measurement model
//!
//! Quantum measurement in the computational basis samples a basis state `k` with probability
//! `|aₖ|²` (the Born rule). After measurement, the register collapses: the measured basis state
//! gets amplitude 1, all others get amplitude 0.
//!
//! # Deterministic seeded sampler
//!
//! The sampler uses `rand_chacha::ChaCha8Rng` seeded with a `u64` seed. This makes KATs
//! reproducible: the same seed + same state always produces the same measurement outcome.
//! For non-test code, use a random seed from `rand::rng()` or an OS source.
//!
//! # Single-qubit vs full-register measurement
//!
//! - [`measure_all`]: sample a full basis state, collapse the entire register.
//! - [`measure_qubit`]: measure a single qubit, collapse only that qubit's degree of freedom
//!   (the remaining qubits are renormalized conditioned on the measured outcome).

use num_complex::Complex;
use rand::{Rng, RngExt};
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;

use crate::statevec::{StateVec, EPS};

/// Outcome of a full-register measurement.
///
/// Contains the measured basis index and the collapsed register.
#[derive(Clone, Debug)]
pub struct MeasureAllOutcome {
    /// The measured basis index (little-endian, qubit 0 = LSB).
    pub basis_index: usize,
    /// The collapsed register: amplitude 1 at `basis_index`, 0 elsewhere.
    pub collapsed: StateVec,
}

/// Outcome of a single-qubit measurement.
#[derive(Clone, Debug)]
pub struct MeasureQubitOutcome {
    /// The measured qubit index.
    pub qubit: usize,
    /// The measured bit value (0 or 1).
    pub bit: u8,
    /// The collapsed, renormalized register.
    pub collapsed: StateVec,
}

/// Construct a `ChaCha8Rng` from a `u64` seed.
///
/// The seed is placed in the low 8 bytes of the 32-byte ChaCha seed array.
#[must_use]
pub fn seeded_rng(seed: u64) -> ChaCha8Rng {
    let mut seed_bytes = [0u8; 32];
    seed_bytes[..8].copy_from_slice(&seed.to_le_bytes());
    ChaCha8Rng::from_seed(seed_bytes)
}

/// Sample a basis index from the Born-rule distribution `P(k) = |aₖ|²`.
///
/// Uses the provided RNG to draw a uniform sample in [0, 1) and walks the cumulative
/// probability distribution to find the sampled basis state.
///
/// # Panics
///
/// Panics if the state is not normalized (total probability < 0.5 — a sanity guard).
fn sample_basis<R: Rng>(sv: &StateVec, rng: &mut R) -> usize {
    let amps = sv.amplitudes();
    let total: f64 = amps.iter().map(|a| a.norm_sqr()).sum();
    assert!(
        (total - 1.0).abs() < 0.5,
        "measure: state is not normalized (Σ|aᵢ|² = {total:.6})"
    );
    let u: f64 = rng.random::<f64>() * total;
    let mut cumulative = 0.0;
    for (i, a) in amps.iter().enumerate() {
        cumulative += a.norm_sqr();
        if u < cumulative {
            return i;
        }
    }
    // Floating-point rounding: return the last nonzero index.
    amps.iter()
        .enumerate()
        .rev()
        .find(|(_, a)| a.norm_sqr() > 0.0)
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Measure all qubits of `sv` using the Born rule with the given RNG.
///
/// Samples a basis state `k` with probability `|aₖ|²`, then collapses the register to |k⟩.
/// Returns the measured basis index and the collapsed register.
///
/// The input register `sv` is consumed (the collapsed state is returned in the outcome).
pub fn measure_all<R: Rng>(sv: StateVec, rng: &mut R) -> MeasureAllOutcome {
    let k = sample_basis(&sv, rng);
    let n = sv.n_qubits();
    let collapsed = StateVec::basis(n, k);
    MeasureAllOutcome { basis_index: k, collapsed }
}

/// Measure all qubits of `sv` using a deterministic seeded RNG.
///
/// Convenience wrapper around [`measure_all`] with a `ChaCha8Rng` seeded from `seed`.
/// Reproducible: the same `sv` + `seed` always produces the same outcome.
pub fn measure_all_seeded(sv: StateVec, seed: u64) -> MeasureAllOutcome {
    let mut rng = seeded_rng(seed);
    measure_all(sv, &mut rng)
}

/// Measure a single qubit `qubit` of `sv` using the Born rule with the given RNG.
///
/// Computes the marginal probability `P(bit=0)` and `P(bit=1)` for the target qubit, samples
/// the outcome, then collapses and renormalizes the register conditioned on the outcome.
///
/// Returns the measured bit (0 or 1) and the collapsed, renormalized register.
///
/// # Panics
///
/// Panics if `qubit >= n` or if the post-collapse renormalization fails (probability ≈ 0).
pub fn measure_qubit<R: Rng>(sv: StateVec, qubit: usize, rng: &mut R) -> MeasureQubitOutcome {
    let n = sv.n_qubits();
    assert!(qubit < n, "qubit {qubit} out of range for {n}-qubit register");
    let mask = 1usize << qubit;
    let amps = sv.amplitudes();

    // Compute P(bit=1) = Σ |aᵢ|² for all i with bit `qubit` set.
    let prob1: f64 = amps
        .iter()
        .enumerate()
        .filter(|(i, _)| i & mask != 0)
        .map(|(_, a)| a.norm_sqr())
        .sum();

    let u: f64 = rng.random::<f64>();
    let bit = if u < prob1 { 1u8 } else { 0u8 };

    // Collapse: zero out amplitudes inconsistent with the measurement, renormalize.
    let keep_mask = if bit == 1 { mask } else { 0 };
    let keep_set = if bit == 1 { mask } else { 0 };
    let norm_sq: f64 = amps
        .iter()
        .enumerate()
        .filter(|(i, _)| i & keep_mask == keep_set)
        .map(|(_, a)| a.norm_sqr())
        .sum();
    assert!(
        norm_sq > EPS * EPS,
        "measure_qubit: post-collapse probability ≈ 0 (norm_sq = {norm_sq:.2e})"
    );
    let norm = norm_sq.sqrt();
    let new_amps: Vec<Complex<f64>> = amps
        .iter()
        .enumerate()
        .map(|(i, &a)| {
            if i & keep_mask == keep_set {
                a / Complex::new(norm, 0.0)
            } else {
                Complex::new(0.0, 0.0)
            }
        })
        .collect();
    let collapsed = StateVec::from_amplitudes(n, new_amps);
    MeasureQubitOutcome { qubit, bit, collapsed }
}

/// Measure a single qubit using a deterministic seeded RNG.
///
/// Convenience wrapper around [`measure_qubit`] with a `ChaCha8Rng` seeded from `seed`.
pub fn measure_qubit_seeded(sv: StateVec, qubit: usize, seed: u64) -> MeasureQubitOutcome {
    let mut rng = seeded_rng(seed);
    measure_qubit(sv, qubit, &mut rng)
}

/// Sample `n_shots` measurements from `sv` using a seeded RNG, returning the frequency table.
///
/// Each shot independently samples a basis state from the Born-rule distribution. The register
/// is not collapsed between shots (the same `sv` is sampled repeatedly). Returns a `Vec<u64>`
/// of length `2^n` where entry `k` is the number of times basis state `k` was sampled.
///
/// This is the statistical KAT helper: over many shots, the empirical frequencies should
/// converge to the Born-rule probabilities `|aₖ|²`.
#[must_use]
pub fn sample_counts(sv: &StateVec, n_shots: usize, seed: u64) -> Vec<u64> {
    let mut rng = seeded_rng(seed);
    let dim = sv.dim();
    let mut counts = vec![0u64; dim];
    for _ in 0..n_shots {
        let k = sample_basis(sv, &mut rng);
        counts[k] += 1;
    }
    counts
}
