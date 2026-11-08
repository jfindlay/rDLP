//! Known-answer tests for Shor's ECDLP algorithm.
//!
//! # What is tested
//!
//! 1. **Order/group-law KATs** — classical verification that the toy curve has the expected
//!    structure: `r·G = ∞`, `k·G = Q` for the chosen discrete-log instances.
//!
//! 2. **Dual-lattice extraction KATs** — exact peak pairs `(c, d)` with `c/M = j/r`,
//!    `d/M = l/r`, `l ≡ k·j (mod r)` recover `k`; `j ≡ 0` yields `None`; out-of-range inputs
//!    are rejected. Deterministic classical KATs, no quantum circuit needed.
//!
//! 3. **Outcome distribution KATs** — the exact exponent-register distribution computed from
//!    circuit amplitudes equals the classical Fourier-sum formula, and the extraction rule has
//!    positive success probability.
//!
//! 4. **End-to-end ECDLP** — `solve_ecdlp(G, Q=k·G)` at fixed measurement seeds returns `k` with
//!    `Measurement` evidence: the reported pair `(c, d)` is in range and
//!    `extract_k_from_measurement(c, d, r) == Some(k)`. The trivial cases report `Trivial`, a
//!    generator at infinity is `InvalidGenerator`, and a budget spent on a measurement whose
//!    candidate is wrong is `AttemptsExhausted` with the exact attempt count.
//!
//! 5. **rho cross-check** (`#[ignore]`-gated for cost) — `rho::ecdlp::solve_brent` on the
//!    equivalent `rho::Curve` returns a `k'` satisfying `k'·G = Q`, and agrees with the Shor
//!    solver's `k` for the same instance.
//!
//! # Toy curve parameters
//!
//! ```text
//! p = 7,  a = 0,  b = 3,  y² = x³ + 3 mod 7
//! G = (1, 2),  r = 13 (prime group order)
//! ```
//!
//! Scalar-multiple table (from `shor::curve` module documentation):
//!
//! ```text
//! 0·G = ∞,      1·G = (1,2),  2·G = (6,3),  3·G = (2,2),
//! 4·G = (4,5),  5·G = (3,3),  6·G = (5,3),  7·G = (5,4),
//! 8·G = (3,4),  9·G = (4,2),  10·G = (2,5), 11·G = (6,4),
//! 12·G = (1,5), 13·G = ∞
//! ```
//!
//! # Qubit budget
//!
//! The two-register ECDLP circuit uses 17 qubits:
//! - 4 (a-register) + 4 (b-register) + 3 (x) + 3 (y) + 3 (λ scratch) = 17
//! - Well within the ~25-qubit simulator ceiling (principle 4).

use std::f64::consts::PI;
use std::sync::OnceLock;

use shor::curve::{self, G, Point, R};
use shor::ecdlp::{
    exponent_outcome_distribution, extract_k_from_measurement, primary_success_probability,
    run_period_finding_circuit, solve_ecdlp, solve_ecdlp_with_retries, EcdlpError, EcdlpEvidence,
    EcdlpSolution, DEFAULT_MAX_ATTEMPTS, EXP_BITS, TOTAL_QUBITS,
};

// ── KAT group 1: order/group-law KATs (classical verification) ────────────────
//
// These KATs verify that the toy curve has the expected structure.
// They are purely classical — no quantum circuit is run.

/// r·G = ∞: the group order is 13.
#[test]
fn order_r_times_g_is_infinity() {
    let result = curve::scalar_mul(R, G);
    assert!(result.is_infinity(), "r·G should be ∞ (group order r={R})");
}

/// (r+1)·G = G: scalar multiplication wraps correctly.
#[test]
fn order_r_plus_1_times_g_is_g() {
    let result = curve::scalar_mul(R + 1, G);
    assert_eq!(result, G, "(r+1)·G should equal G");
}

/// k·G = Q for the known scalar-multiple table.
///
/// Verifies the discrete-log instances used in the end-to-end KATs are well-formed.
#[test]
fn group_law_scalar_multiples() {
    // Published scalar-multiple table from the curve module documentation.
    let expected: &[(u64, Point)] = &[
        (0, Point::Infinity),
        (1, Point::Affine { x: 1, y: 2 }),
        (2, Point::Affine { x: 6, y: 3 }),
        (3, Point::Affine { x: 2, y: 2 }),
        (4, Point::Affine { x: 4, y: 5 }),
        (5, Point::Affine { x: 3, y: 3 }),
        (6, Point::Affine { x: 5, y: 3 }),
        (7, Point::Affine { x: 5, y: 4 }),
        (8, Point::Affine { x: 3, y: 4 }),
        (9, Point::Affine { x: 4, y: 2 }),
        (10, Point::Affine { x: 2, y: 5 }),
        (11, Point::Affine { x: 6, y: 4 }),
        (12, Point::Affine { x: 1, y: 5 }),
        (13, Point::Infinity),
    ];
    for &(k, expected_pt) in expected {
        let got = curve::scalar_mul(k, G);
        assert_eq!(
            got, expected_pt,
            "{k}·G = {got:?}, expected {expected_pt:?}"
        );
    }
}

/// k·G = Q for the specific instances used in the end-to-end KATs.
///
/// Verifies that the discrete-log instances are well-formed before running the quantum circuit.
#[test]
fn ecdlp_instances_well_formed() {
    // The end-to-end KATs use k = 3, 5, 7, 9.
    for k in [3u64, 5, 7, 9] {
        let q = curve::scalar_mul(k, G);
        assert!(!q.is_infinity(), "k={k}: k·G should not be ∞ (k < r={R})");
        assert!(curve::on_curve(q), "k={k}: k·G should be on the curve");
        // Verify the inverse: scalar_mul(k, G) = Q, and Q is not G (k ≠ 1).
        assert_ne!(q, G, "k={k}: k·G should not equal G (k ≠ 1)");
    }
}

/// Qubit budget: 17 qubits, within the ~25-qubit ceiling.
#[test]
fn qubit_budget_within_ceiling() {
    assert_eq!(EXP_BITS, 4, "EXP_BITS should be 4 (⌈log₂ 13⌉ = 4)");
    assert_eq!(TOTAL_QUBITS, 17, "TOTAL_QUBITS should be 17");
    assert!(
        TOTAL_QUBITS <= 25,
        "TOTAL_QUBITS={TOTAL_QUBITS} exceeds the ~25-qubit ceiling"
    );
}

// ── KAT group 2: dual-lattice extraction KATs (classical, deterministic) ──────
//
// Exact peak pairs (c, d) with c ≈ j·M/r, d ≈ l·M/r, l ≡ k·j (mod r) must recover k.
// These are purely classical KATs.

/// Number of values per exponent register, `M = 2^EXP_BITS`.
const M: u64 = 1 << EXP_BITS;

/// Nearest integer to `x·M/r`: the measured value closest to the dual-lattice point `x/r`.
fn peak(x: u64) -> u64 {
    (2 * x * M + R) / (2 * R)
}

/// Every peak value fits in the exponent register.
#[test]
fn extract_peaks_fit_in_register() {
    for x in 0..R {
        assert!(peak(x) < M, "peak({x}) = {} is out of range", peak(x));
    }
}

/// extract_k: for every k and every non-zero lattice coordinate j, the exact peak pair
/// `(peak(j), peak(k·j mod r))` recovers k.
#[test]
fn extract_exact_peak_pairs_recover_k() {
    for k in 0..R {
        for j in 1..R {
            let l = k * j % R;
            let got = extract_k_from_measurement(peak(j), peak(l), R);
            assert_eq!(got, Some(k), "k={k}, j={j}, l={l}: peak pair gave {got:?}");
        }
    }
}

/// extract_k: a known pair for k=3. With j=1, l=3: c=peak(1)=1, d=peak(3)=4.
#[test]
fn extract_k3_from_known_pair() {
    assert_eq!(extract_k_from_measurement(1, 4, R), Some(3));
}

/// extract_k: a known pair for k=7. With j=2, l=14 mod 13=1: c=peak(2)=2, d=peak(1)=1.
#[test]
fn extract_k7_from_known_pair() {
    assert_eq!(extract_k_from_measurement(2, 1, R), Some(7));
}

/// extract_k: j ≡ 0 is uninformative → None, whatever d is.
///
/// For M=16, r=13 only c=0 rounds to j ≡ 0 (c=15 rounds to 12).
#[test]
fn extract_j_zero_returns_none() {
    for d in 0..M {
        assert_eq!(extract_k_from_measurement(0, d, R), None, "c=0, d={d} should be None");
    }
    // Every other c rounds to a non-zero j, so yields a candidate below r.
    for c in 1..M {
        for d in 0..M {
            let k = extract_k_from_measurement(c, d, R)
                .unwrap_or_else(|| panic!("c={c}, d={d} should give a candidate"));
            assert!(k < R, "c={c}, d={d}: k={k} must be < r={R}");
        }
    }
}

/// extract_k: a-register value ≥ M is invalid input and is rejected.
#[test]
#[should_panic(expected = "must be < 2^EXP_BITS")]
fn extract_rejects_c_out_of_range() {
    let _ = extract_k_from_measurement(M, 0, R);
}

/// extract_k: b-register value ≥ M is invalid input and is rejected.
#[test]
#[should_panic(expected = "must be < 2^EXP_BITS")]
fn extract_rejects_d_out_of_range() {
    let _ = extract_k_from_measurement(1, M, R);
}

// ── KAT group 3: exact outcome distribution ───────────────────────────────────
//
// The circuit's exponent-register distribution is compared against the classical Fourier-sum
// formula P(c,d) = Σ_s |M⁻² Σ_{a+bk≡s (mod r)} ω^{-(ac+bd)}|², ω = e^{2πi/M}, computed here
// without any circuit code.

/// Classical Fourier-sum outcome distribution for logarithm `k`, indexed `c + (d << EXP_BITS)`.
fn fourier_sum_distribution(k: u64) -> Vec<f64> {
    let m = M as usize;
    let norm = 1.0 / (M * M) as f64;
    let mut dist = vec![0.0; m * m];
    for c in 0..m {
        for d in 0..m {
            // Per-class complex sums: (re, im) indexed by s = (a + b·k) mod r.
            let mut classes = vec![(0.0f64, 0.0f64); R as usize];
            for a in 0..M {
                for b in 0..M {
                    let s = ((a + b * k) % R) as usize;
                    let phase = -2.0 * PI * ((a as usize * c + b as usize * d) % m) as f64 / M as f64;
                    classes[s].0 += phase.cos();
                    classes[s].1 += phase.sin();
                }
            }
            dist[c + (d << EXP_BITS)] = classes
                .iter()
                .map(|&(re, im)| (re * norm).powi(2) + (im * norm).powi(2))
                .sum();
        }
    }
    dist
}

/// Circuit distribution for `k`, computed once per `k` and shared across tests.
fn cached_circuit_distribution(k: u64) -> &'static [f64] {
    static CACHE: [OnceLock<Vec<f64>>; R as usize] = [const { OnceLock::new() }; R as usize];
    CACHE[k as usize].get_or_init(|| {
        let q = curve::scalar_mul(k, G);
        exponent_outcome_distribution(G, q)
    })
}

/// Assert the circuit distribution for `k` equals the Fourier-sum formula and sums to 1.
fn assert_distribution_matches_formula(k: u64) {
    let circuit = cached_circuit_distribution(k);
    let formula = fourier_sum_distribution(k);
    assert_eq!(circuit.len(), 1 << (2 * EXP_BITS));
    let total: f64 = circuit.iter().sum();
    assert!((total - 1.0).abs() < 1e-9, "k={k}: circuit distribution sums to {total}");
    for (idx, (&got, &want)) in circuit.iter().zip(&formula).enumerate() {
        assert!(
            (got - want).abs() < 1e-9,
            "k={k}, c={}, d={}: circuit {got} vs formula {want}",
            idx % M as usize,
            idx / M as usize
        );
    }
}

/// Circuit outcome distribution equals the Fourier-sum formula for k=3.
#[test]
fn distribution_matches_fourier_formula_k3() {
    assert_distribution_matches_formula(3);
}

/// Circuit outcome distribution equals the Fourier-sum formula for k=7.
#[test]
fn distribution_matches_fourier_formula_k7() {
    assert_distribution_matches_formula(7);
}

/// Circuit outcome distribution equals the Fourier-sum formula for every k in 2..r.
#[test]
#[ignore = "cost: one full circuit simulation per k (~10 s each in debug)"]
fn distribution_matches_fourier_formula_all_k() {
    for k in 2..R {
        assert_distribution_matches_formula(k);
    }
}

/// One-candidate extraction succeeds with positive probability for k=3 and k=7.
#[test]
fn primary_success_probability_positive_k3_k7() {
    for k in [3u64, 7] {
        let q = curve::scalar_mul(k, G);
        let cached = cached_circuit_distribution(k);
        let p = primary_success_probability(G, q);
        assert!(p > 0.0, "k={k}: primary success probability is {p}");
        eprintln!("k={k}: primary success probability = {p}");
        // Cross-check against the cached distribution and the public extraction rule.
        let expected: f64 = cached
            .iter()
            .enumerate()
            .filter(|&(idx, _)| {
                let (c, d) = ((idx as u64) % M, (idx as u64) / M);
                extract_k_from_measurement(c, d, R) == Some(k)
            })
            .map(|(_, &p)| p)
            .sum();
        assert!((p - expected).abs() < 1e-12, "k={k}: {p} vs {expected}");
    }
}

/// Primary success probability is positive for every k in 2..r (values printed).
#[test]
#[ignore = "cost: one full circuit simulation per k (~10 s each in debug)"]
fn primary_success_probability_positive_all_k() {
    for k in 2..R {
        let q = curve::scalar_mul(k, G);
        let p = primary_success_probability(G, q);
        eprintln!("k={k}: primary success probability = {p}");
        assert!(p > 0.0, "k={k}: primary success probability is {p}");
    }
}

// ── KAT group 4: end-to-end ECDLP (quantum circuit + classical post-processing) ──
//
// These KATs run the full Shor ECDLP algorithm end-to-end at fixed seeds. A solution is
// `Measurement` evidence: the measured pair `(c, d)` must itself extract to `k`, so the test
// checks the mechanism's own identity and not only the self-certifying `k·G == Q`.

/// Assert `solution` is `Measurement` evidence for logarithm `k` of `q`, obtained from the
/// circuit run at `base_seed + attempt`, and return the measured pair.
fn assert_measurement_evidence(
    solution: &EcdlpSolution,
    k: u64,
    q: Point,
    base_seed: u64,
    label: &str,
) -> (u64, u64) {
    assert_eq!(solution.k, k, "{label}: wrong exponent");
    assert_eq!(curve::scalar_mul(solution.k, G), q, "{label}: k·G ≠ Q");
    let EcdlpEvidence::Measurement { seed, attempt, c, d } = solution.evidence else {
        panic!("{label}: expected Measurement evidence, got {:?}", solution.evidence);
    };
    assert_eq!(seed, base_seed + attempt as u64, "{label}: seed is not base seed + attempt");
    assert!(attempt < DEFAULT_MAX_ATTEMPTS, "{label}: attempt {attempt} beyond the default budget");
    assert!(c < M && d < M, "{label}: measured pair ({c}, {d}) out of range");
    assert_eq!(
        extract_k_from_measurement(c, d, R),
        Some(k),
        "{label}: reported pair ({c}, {d}) does not extract to k={k}"
    );
    (c, d)
}

/// Assert `solve_ecdlp(G, k·G, 0)` succeeds on its first measurement with the recorded pair.
///
/// Seed 0 was kept for every fixture after a recorded search over seeds 0..16 for
/// k ∈ {3, 5, 7, 9, 12}: seed 0 extracts the correct `k` on its first attempt in all five cases,
/// so each fixture costs one circuit run. The expected pair is the one seed 0 measures.
fn assert_solve_first_attempt(k: u64, expected_pair: (u64, u64)) {
    let q = curve::scalar_mul(k, G);
    let solution = solve_ecdlp(G, q, 0).unwrap_or_else(|e| panic!("k={k}: solve_ecdlp failed: {e:?}"));
    let pair = assert_measurement_evidence(&solution, k, q, 0, &format!("k={k}"));
    assert_eq!(pair, expected_pair, "k={k}: seed 0 measured a different pair");
    assert_eq!(
        solution.evidence,
        EcdlpEvidence::Measurement { seed: 0, attempt: 0, c: pair.0, d: pair.1 }
    );
}

/// Q=∞: k=0 with `Trivial` evidence, no circuit run.
#[test]
fn ecdlp_trivial_q_infinity() {
    assert_eq!(
        solve_ecdlp(G, Point::Infinity, 0),
        Ok(EcdlpSolution { k: 0, evidence: EcdlpEvidence::Trivial })
    );
}

/// Q=G: k=1 with `Trivial` evidence, no circuit run.
#[test]
fn ecdlp_trivial_q_equals_g() {
    assert_eq!(solve_ecdlp(G, G, 0), Ok(EcdlpSolution { k: 1, evidence: EcdlpEvidence::Trivial }));
}

/// G=∞ is an invalid generator, whatever Q is; no circuit run.
#[test]
fn ecdlp_invalid_generator() {
    for q in [Point::Infinity, G] {
        assert_eq!(solve_ecdlp(Point::Infinity, q, 0), Err(EcdlpError::InvalidGenerator));
    }
}

/// k=3, Q=3·G=(2,2), seed 0.
#[test]
fn ecdlp_k3_seed0() {
    assert_solve_first_attempt(3, (10, 13));
}

/// k=5, Q=5·G=(3,3), seed 0.
#[test]
fn ecdlp_k5_seed0() {
    assert_solve_first_attempt(5, (12, 13));
}

/// k=7, Q=7·G=(5,4), seed 0.
#[test]
fn ecdlp_k7_seed0() {
    assert_solve_first_attempt(7, (9, 12));
}

/// k=9, Q=9·G=(4,2), seed 0.
#[test]
fn ecdlp_k9_seed0() {
    assert_solve_first_attempt(9, (8, 13));
}

/// k=12, Q=12·G=(1,5), seed 0. k=12 = r−1 is the largest non-trivial scalar.
#[test]
fn ecdlp_k12_seed0() {
    assert_solve_first_attempt(12, (4, 12));
}

/// One attempt whose single candidate is wrong exhausts the budget: `AttemptsExhausted` with the
/// exact count, never a recovered `k`.
///
/// For k=3, seed 2 measures (2, 11), which extracts to 11 ≠ 3 (found by the same bounded search
/// over seeds 0..16 as above; seeds 2, 7, 9 and 11 fail for k=3). With a budget of one attempt
/// the failure is therefore certain, not probabilistic. The old exhaustive fallback would have
/// returned k=3 here.
#[test]
fn ecdlp_budget_exhausted_on_wrong_candidate() {
    let q = curve::scalar_mul(3, G);
    assert_eq!(
        solve_ecdlp_with_retries(G, q, 2, 1),
        Err(EcdlpError::AttemptsExhausted { attempts: 1 })
    );
}

/// A failed attempt is followed by a fresh-seed attempt: k=3 from seed 2 fails at seed 2
/// (candidate 11) and succeeds at seed 3 (pair (13, 9) extracts to 3), so the evidence reports
/// attempt 1 and seed 3.
#[test]
fn ecdlp_retry_uses_fresh_seed() {
    let q = curve::scalar_mul(3, G);
    let solution = solve_ecdlp_with_retries(G, q, 2, 2).expect("seed 3 extracts k=3");
    let (c, d) = assert_measurement_evidence(&solution, 3, q, 2, "retry k=3");
    assert_eq!(solution.evidence, EcdlpEvidence::Measurement { seed: 3, attempt: 1, c, d });
    assert_eq!((c, d), (13, 9));
}

/// A zero budget makes no measurement and reports exhaustion with zero attempts.
#[test]
fn ecdlp_zero_budget_is_exhausted() {
    let q = curve::scalar_mul(3, G);
    assert_eq!(
        solve_ecdlp_with_retries(G, q, 0, 0),
        Err(EcdlpError::AttemptsExhausted { attempts: 0 })
    );
}

/// Every non-trivial k (2..=12) is solved at seed 0 within the default budget, with
/// `Measurement` evidence that extracts to `k`.
///
/// Ignored for cost: at least one full circuit simulation per k (~10 s each in debug mode). The
/// fixtures above cover k = 3, 5, 7, 9, 12.
#[test]
#[ignore = "cost: at least one full circuit simulation per k (~10 s each in debug)"]
fn ecdlp_all_k_values_within_default_budget() {
    for k in 2..R {
        let q = curve::scalar_mul(k, G);
        let solution = solve_ecdlp(G, q, 0)
            .unwrap_or_else(|e| panic!("k={k}: solve_ecdlp failed: {e:?}"));
        assert_measurement_evidence(&solution, k, q, 0, &format!("k={k}"));
    }
}

/// solve_ecdlp is reproducible: same seed → same solution, evidence included.
///
/// Ignored for cost: two full solves per k (several circuit simulations, ~10 s each in debug).
/// Reproducibility comes from the deterministic seeded RNG in `measure_all_seeded`.
#[test]
#[ignore = "cost: two full solves per k (~10 s per circuit simulation in debug)"]
fn ecdlp_reproducible() {
    for k in [3u64, 7, 9] {
        let q = curve::scalar_mul(k, G);
        let r1 = solve_ecdlp(G, q, 0);
        let r2 = solve_ecdlp(G, q, 0);
        assert_eq!(r1, r2, "solve_ecdlp(G, {k}·G, seed=0) is not reproducible");
    }
}

/// run_period_finding_circuit: the measured pair is in range and has non-zero probability
/// under the exact outcome distribution.
///
/// Runs the circuit for k=3 with seed=0 and checks the sampled pair against the distribution.
#[test]
fn period_finding_sample_in_support_of_distribution() {
    let k = 3u64;
    let q = curve::scalar_mul(k, G);
    let (c, d) = run_period_finding_circuit(G, q, 0);
    assert!(c < M && d < M, "measured pair ({c}, {d}) out of range");
    let p = cached_circuit_distribution(k)[(c + (d << EXP_BITS)) as usize];
    assert!(p > 1e-12, "measured pair ({c}, {d}) has probability {p}");
}

// ── KAT group 5: rho cross-check (#[ignore]-gated for cost) ───────────────────
//
// These KATs run `rho::ecdlp::solve_brent` on the same toy instance and compare it with the
// Shor solver. They are not oracle tests: `rho` is an in-workspace dev-dependency, no external
// tool is needed. They are `#[ignore]`-gated as a group; the cost is in the Shor comparison, whose
// full solve runs circuit simulations (~10 s each in debug). The two rho-only tests take about
// 0.2 s and stay ignored only because they belong to this group.
//
// The rho::Curve struct has fields: p, a, b, n, gx, gy (all Uint<4>).
// The toy curve parameters: p=7, a=0, b=3, n=13, gx=1, gy=2.

/// rho cross-check: solve_brent recovers k'·G = Q for k=3 on the toy curve.
///
/// Builds a rho::Curve from the toy curve parameters and calls solve_brent.
/// The returned k' must satisfy k'·G = Q (relationship-preservation).
#[test]
#[ignore = "cross-check group kept off the default path; cheap on its own (~0.2 s), unlike the Shor comparison"]
fn rho_crosscheck_k3() {
    use crypto_bigint::Uint;
    use rho::curve::{AffinePoint, Curve};
    use rho::ecdlp::solve_brent;
    use shared_field::FpMonty4 as FpMonty;

    // Build the rho::Curve from the toy curve parameters.
    let toy_curve = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(3u64),
        n: Uint::<4>::from(13u64),
        gx: Uint::<4>::from(1u64),
        gy: Uint::<4>::from(2u64),
    };
    let g_rho: AffinePoint<FpMonty> = toy_curve.generator();

    let k_target = 3u64;
    let q_rho = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k_target));

    let k = solve_brent(&toy_curve, &g_rho, &q_rho, 13, 0, 20)
        .expect("rho::solve_brent failed for k=3 on toy curve");

    // Verify relationship-preservation: k'·G = Q.
    let check = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k));
    assert_eq!(check, q_rho, "rho cross-check: k'·G ≠ Q (k'={k}, expected k={k_target})");
}

/// rho cross-check: solve_brent recovers k'·G = Q for k=7 on the toy curve.
#[test]
#[ignore = "cross-check group kept off the default path; cheap on its own (~0.2 s), unlike the Shor comparison"]
fn rho_crosscheck_k7() {
    use crypto_bigint::Uint;
    use rho::curve::{AffinePoint, Curve};
    use rho::ecdlp::solve_brent;
    use shared_field::FpMonty4 as FpMonty;

    let toy_curve = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(3u64),
        n: Uint::<4>::from(13u64),
        gx: Uint::<4>::from(1u64),
        gy: Uint::<4>::from(2u64),
    };
    let g_rho: AffinePoint<FpMonty> = toy_curve.generator();

    let k_target = 7u64;
    let q_rho = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k_target));

    let k = solve_brent(&toy_curve, &g_rho, &q_rho, 13, 0, 20)
        .expect("rho::solve_brent failed for k=7 on toy curve");

    let check = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k));
    assert_eq!(check, q_rho, "rho cross-check: k'·G ≠ Q (k'={k}, expected k={k_target})");
}

/// rho cross-check: solve_brent and solve_ecdlp agree on k for k=5.
///
/// Both solvers must return a k' satisfying k'·G = Q = 5·G. The group has prime order 13, so the
/// logarithm in `0..13` is unique and the two values must be equal; the Shor result also carries
/// `Measurement` evidence that extracts to it.
#[test]
#[ignore = "cost: a full Shor solve (circuit simulations, ~10 s each in debug) plus a rho walk"]
fn rho_crosscheck_agrees_with_shor_k5() {
    use crypto_bigint::Uint;
    use rho::curve::{AffinePoint, Curve};
    use rho::ecdlp::solve_brent;
    use shared_field::FpMonty4 as FpMonty;

    let toy_curve = Curve {
        p: Uint::<4>::from(7u64),
        a: Uint::<4>::from(0u64),
        b: Uint::<4>::from(3u64),
        n: Uint::<4>::from(13u64),
        gx: Uint::<4>::from(1u64),
        gy: Uint::<4>::from(2u64),
    };
    let g_rho: AffinePoint<FpMonty> = toy_curve.generator();

    let k_target = 5u64;
    let q_rho = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k_target));

    // rho solver.
    let k_rho = solve_brent(&toy_curve, &g_rho, &q_rho, 13, 0, 20)
        .expect("rho::solve_brent failed for k=5");
    let check_rho = toy_curve.scalar_mul(&g_rho, &Uint::<4>::from(k_rho));
    assert_eq!(check_rho, q_rho, "rho: k'·G ≠ Q (k'={k_rho})");

    // Shor solver, with its measurement evidence.
    let q_shor = curve::scalar_mul(k_target, G);
    let solution = solve_ecdlp(G, q_shor, 0).expect("solve_ecdlp failed for k=5");
    assert_measurement_evidence(&solution, k_target, q_shor, 0, "shor k=5");

    // Prime group order: the logarithm in 0..r is unique, so both solvers return it.
    assert_eq!(solution.k, k_rho % R, "shor k={} vs rho k={k_rho}", solution.k);
}
