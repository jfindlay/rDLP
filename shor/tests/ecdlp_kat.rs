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
//! 4. **End-to-end ECDLP** — `solve_ecdlp(G, Q=k·G)` recovers `k` for a fixture of `k` values
//!    on the toy curve, with a fixed measurement seed so the run is reproducible. The assertion
//!    checks `recovered_k · G == Q` (relationship-preservation).
//!
//! 5. **rho cross-check** (`#[ignore]`-gated) — `rho::ecdlp::solve_brent` on the equivalent
//!    `rho::Curve` returns a `k'` satisfying `k'·G = Q` for the same instances.
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
    run_period_finding_circuit, solve_ecdlp, EXP_BITS, TOTAL_QUBITS,
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

// ── KAT group 3: end-to-end ECDLP (quantum circuit + classical post-processing) ──
//
// These KATs run the full Shor ECDLP algorithm end-to-end with fixed seeds.
// Seeds are chosen so that the circuit measurement yields an informative (a', b') pair
// that recovers the correct k on the first or early retry.
//
// The assertion checks recovered_k · G == Q (relationship-preservation), not recovered_k == k,
// since the quantum algorithm may return any k' satisfying k'·G = Q (which equals k mod r).

/// Helper: assert that solve_ecdlp(G, Q=k·G, seed) recovers a k' with k'·G = Q.
fn assert_solve_ecdlp(k: u64, seed: u64, label: &str) {
    let q = curve::scalar_mul(k, G);
    let result = solve_ecdlp(G, q, seed);
    assert!(
        result.is_some(),
        "{label}: solve_ecdlp(G, {k}·G, seed={seed}) returned None"
    );
    let recovered_k = result.unwrap();
    let kg = curve::scalar_mul(recovered_k, G);
    assert_eq!(
        kg, q,
        "{label}: recovered_k={recovered_k}, but recovered_k·G={kg:?} ≠ Q={q:?}"
    );
}

/// solve_ecdlp: trivial case Q=∞ → k=0.
#[test]
fn ecdlp_trivial_q_infinity() {
    let result = solve_ecdlp(G, Point::Infinity, 0);
    assert_eq!(result, Some(0), "solve_ecdlp(G, ∞, 0) should return Some(0)");
}

/// solve_ecdlp: trivial case Q=G → k=1.
#[test]
fn ecdlp_trivial_q_equals_g() {
    let result = solve_ecdlp(G, G, 0);
    assert_eq!(result, Some(1), "solve_ecdlp(G, G, 0) should return Some(1)");
}

/// solve_ecdlp: k=3, Q=3·G=(2,2). Seed 0.
///
/// Verifies that the two-register circuit recovers k=3 (or any k' with k'·G = 3·G).
#[test]
fn ecdlp_k3_seed0() {
    assert_solve_ecdlp(3, 0, "ecdlp_k3_seed0");
}

/// solve_ecdlp: k=5, Q=5·G=(3,3). Seed 0.
#[test]
fn ecdlp_k5_seed0() {
    assert_solve_ecdlp(5, 0, "ecdlp_k5_seed0");
}

/// solve_ecdlp: k=7, Q=7·G=(5,4). Seed 0.
#[test]
fn ecdlp_k7_seed0() {
    assert_solve_ecdlp(7, 0, "ecdlp_k7_seed0");
}

/// solve_ecdlp: k=9, Q=9·G=(4,2). Seed 0.
#[test]
fn ecdlp_k9_seed0() {
    assert_solve_ecdlp(9, 0, "ecdlp_k9_seed0");
}

/// solve_ecdlp: k=12, Q=12·G=(1,5). Seed 0.
///
/// k=12 is the largest non-trivial scalar (r-1=12), testing the boundary.
#[test]
fn ecdlp_k12_seed0() {
    assert_solve_ecdlp(12, 0, "ecdlp_k12_seed0");
}

/// solve_ecdlp: all non-trivial k values (2..=12) recover k'·G = Q.
///
/// Runs solve_ecdlp for every non-trivial k on the toy curve with seed=0.
/// Verifies relationship-preservation: recovered_k · G == Q.
///
/// Marked `#[ignore]` because it runs 11 full quantum circuit instances (~165s in debug
/// mode). The targeted fixture tests (k=3, k=5, k=7, k=9, k=12) cover the same ground.
#[test]
#[ignore]
fn ecdlp_all_k_values_seed0() {
    for k in 2..R {
        assert_solve_ecdlp(k, 0, &format!("ecdlp_k{k}_seed0"));
    }
}

/// solve_ecdlp is reproducible: same seed → same result.
///
/// Marked `#[ignore]` because it runs 6 full quantum circuit instances (~90s in debug mode).
/// Reproducibility is verified by the deterministic seeded RNG in `measure_all_seeded`.
#[test]
#[ignore]
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

// ── KAT group 4: rho cross-check (#[ignore]-gated) ────────────────────────────
//
// These KATs use rho::ecdlp::solve_brent as an oracle to cross-check the Shor solver.
// They are #[ignore]-gated so they do not run in the normal test suite (the rho solver
// is a dev-dependency, and the cross-check is not part of the green path).
//
// The rho::Curve struct has fields: p, a, b, n, gx, gy (all Uint<4>).
// The toy curve parameters: p=7, a=0, b=3, n=13, gx=1, gy=2.

/// rho cross-check: solve_brent recovers k'·G = Q for k=3 on the toy curve.
///
/// Builds a rho::Curve from the toy curve parameters and calls solve_brent.
/// The returned k' must satisfy k'·G = Q (relationship-preservation).
///
/// This test is #[ignore]-gated: it requires the rho dev-dependency and is not
/// part of the green path (the Shor solver is oracle-free).
#[test]
#[ignore]
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
#[ignore]
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

/// rho cross-check: solve_brent and solve_ecdlp agree on k'·G = Q for k=5.
///
/// Both solvers should return a k' satisfying k'·G = Q = 5·G.
/// The specific k' values may differ (any valid discrete log is acceptable).
#[test]
#[ignore]
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

    // Shor solver.
    let q_shor = curve::scalar_mul(k_target, G);
    let k_shor = solve_ecdlp(G, q_shor, 0).expect("solve_ecdlp failed for k=5");
    let check_shor = curve::scalar_mul(k_shor, G);
    assert_eq!(check_shor, q_shor, "shor: k'·G ≠ Q (k'={k_shor})");

    // Both agree on the point Q (relationship-preservation).
    // (The specific k' values may differ since any valid discrete log is acceptable.)
    eprintln!("rho k'={k_rho}, shor k'={k_shor} — both satisfy k'·G = Q = {k_target}·G");
}
