//! Lattice sieving for GNFS (demonstration fidelity).
//!
//! The lattice sieve is an optimization layer over the special-q strategy. For a
//! special prime ``q`` with root ``r_q`` (so ``q | F(a, b)`` iff ``a ≡ r_q·b (mod q)``, with
//! ``F`` the algebraic value, ``N_alg`` in code),
//! the sieve is restricted to the **lattice**
//!
//! ```text
//! L_q = { (a, b) ∈ ℤ² : a ≡ r_q·b (mod q) }
//! ```
//!
//! rather than the full ``(a, b)`` rectangle. The lattice sieve enumerates ``L_q`` via a
//! **reduced basis** (the short vectors from a 2D Gauss lattice reduction), which covers the
//! sieve region more efficiently than stepping through ``a`` in increments of ``q``.
//!
//! # Algorithm
//!
//! 1. **Lattice construction.** A basis for ``L_q`` is:
//!
//!    ```text
//!    v1 = (q, 0)
//!    v2 = (r_q, 1)
//!    ```
//!
//!    Check: ``r_q·1 ≡ r_q (mod q)`` ✓ and ``q·0 ≡ 0 (mod q)`` ✓. Every lattice point
//!    ``s·v1 + t·v2 = (s·q + t·r_q, t)`` satisfies ``a ≡ r_q·b (mod q)`` since
//!    ``s·q + t·r_q ≡ t·r_q ≡ r_q·t (mod q)``.
//!
//! 2. **Gauss lattice reduction (2D).** Apply Gauss's algorithm (the 2D analogue of LLL) to
//!    find a reduced basis ``(V1, V2)`` with shorter vectors:
//!
//!    ```text
//!    while |V1| > |V2|: swap V1, V2
//!    V1 = V1 - round(dot(V1, V2) / dot(V2, V2)) * V2
//!    ```
//!
//!    Repeat until convergence. The reduced basis has the property that ``|V1| ≤ |V2|`` and
//!    ``|dot(V1, V2)| ≤ |V2|² / 2`` (Gauss-reduced condition).
//!
//! 3. **Lattice enumeration.** Enumerate lattice points ``(a, b) = s·V1 + t·V2`` for integer
//!    ``(s, t)`` in a bounded region. The bounds ``S``, ``T`` are chosen so that the enumerated
//!    region covers ``|a| ≤ A`` and ``1 ≤ b ≤ B``.
//!
//! 4. **Sieve over lattice points.** For each enumerated ``(a, b)``:
//!
//!    - Check ``gcd(a, b) = 1``.
//!    - Check ``b ≥ 1`` (b must be positive).
//!    - Compute the rational norm and the algebraic value ``F(a, b) = c·N_{K/ℚ}(a − b·α)``.
//!    - Trial-divide both.
//!    - If both are fully smooth, construct a ``Relation``.
//!
//! # Toy scale versus cryptographic scale
//!
//! The lattice sieve's yield advantage over the line sieve comes from covering ``L_q`` more
//! efficiently: the reduced basis vectors are shorter than the original basis ``(v1, v2)``,
//! so the enumeration visits fewer lattice points outside the sieve region ``|a| ≤ A``,
//! ``1 ≤ b ≤ B``. At cryptographic scale, this efficiency gain is significant: the reduced
//! basis has vectors of length ``≈ √q``, so the enumeration covers ``≈ A·B / q`` lattice
//! points (the density of ``L_q`` in the rectangle), with minimal waste.
//!
//! **At toy scale, this advantage is not visible.** The lattice enumeration covers the same
//! ``(a, b)`` pairs as the special-q line sieve restricted to ``a ≡ r_q·b (mod q)`` — the
//! two algorithms are mathematically equivalent for the same ``(q, r_q)``. The efficiency
//! difference (reduced basis vs. stepping by ``q``) is a constant factor that is swamped by
//! the overhead of the reduction and enumeration at small ``q``. The KATs state the same
//! caveat.
//!
//! # Relation to special-q sieve
//!
//! The lattice sieve is a variant of the special-q sieve with a different enumeration
//! strategy. Both restrict to ``L_q``; the special-q sieve steps through ``a`` in increments
//! of ``q`` for each ``b``, while the lattice sieve enumerates via the reduced basis. The
//! ``Relation`` type, ``FactorBase``, and smoothness predicates are unchanged (the baseline
//! stays available for benchmarking).
//!
//! # Entry surface
//!
//! - [`lattice_sieve`] — the main entry point.
//! - [`LatticeSieveConfig`] — sieve parameters ``(A, B, q_min, q_max, threshold_scale)``.
//! - [`LatticeSieveResult`] — per-``(q, r_q)`` sieve output (relations + metadata).
//! - [`LatticeBasis`] — the reduced basis for ``L_q`` (exposed for KAT inspection).

use num_bigint::BigInt;
use num_traits::Zero;
use shared_numth::trial_smooth;

use super::{
    factor_base::FactorBase,
    norms::{algebraic_norm, norm_sign, norm_to_uint, rational_norm},
    Relation,
};
use crate::polyselect::PolyPair;

// ─── LatticeSieveConfig ───────────────────────────────────────────────────────

/// Configuration for the lattice sieve.
///
/// Controls the sieve region ``|a| ≤ A``, ``1 ≤ b ≤ B``, the special-q range
/// ``[q_min, q_max]``, and the threshold scale factor.
///
/// # Special-q range
///
/// The special primes are drawn from the algebraic factor base: only primes ``q`` in the
/// algebraic factor base with ``q_min ≤ q ≤ q_max`` are used. This ensures that ``q``
/// is already in the factor base and its index is known for the exponent vector.
///
/// # Threshold
///
/// Same as ``SpecialQConfig``: the threshold is ``threshold_scale × log2(B_alg)``.
/// A value of ``0.8`` (default) is conservative.
///
/// # Toy scale versus cryptographic scale
///
/// At toy scale, the threshold barely filters anything — the asymptotic win of the sieve
/// is under-exposed. See the module-level documentation.
#[derive(Debug, Clone)]
pub struct LatticeSieveConfig {
    /// Half-width of the sieve region: ``a`` ranges over ``−A..=A``.
    pub a_bound: u64,
    /// Height of the sieve region: ``b`` ranges over ``1..=B``.
    pub b_bound: u64,
    /// Minimum special prime (inclusive). Must be ≥ 2.
    pub q_min: u64,
    /// Maximum special prime (inclusive).
    pub q_max: u64,
    /// Scale factor for the smoothness threshold (default: 0.8).
    ///
    /// The threshold is ``threshold_scale × log2(B_alg)``. Lower values accept more
    /// candidates for trial division (slower but no missed relations); higher values
    /// filter more aggressively (faster but may miss relations with small norms).
    pub threshold_scale: f64,
}

impl LatticeSieveConfig {
    /// Construct a ``LatticeSieveConfig`` with the given bounds and default threshold scale (0.8).
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :param q_min: Minimum special prime (inclusive).
    /// :param q_max: Maximum special prime (inclusive).
    /// :returns: A new ``LatticeSieveConfig``.
    pub fn new(a_bound: u64, b_bound: u64, q_min: u64, q_max: u64) -> Self {
        Self { a_bound, b_bound, q_min, q_max, threshold_scale: 0.8 }
    }

    /// Construct a ``LatticeSieveConfig`` with an explicit threshold scale.
    ///
    /// :param a_bound: Half-width of the sieve region.
    /// :param b_bound: Height of the sieve region.
    /// :param q_min: Minimum special prime (inclusive).
    /// :param q_max: Maximum special prime (inclusive).
    /// :param threshold_scale: Scale factor for the smoothness threshold (0.0–1.0).
    /// :returns: A new ``LatticeSieveConfig``.
    pub fn with_threshold(
        a_bound: u64,
        b_bound: u64,
        q_min: u64,
        q_max: u64,
        threshold_scale: f64,
    ) -> Self {
        Self { a_bound, b_bound, q_min, q_max, threshold_scale }
    }
}

// ─── LatticeBasis ─────────────────────────────────────────────────────────────

/// A reduced basis for the special-q lattice ``L_q``.
///
/// The lattice ``L_q = { (a, b) : a ≡ r_q·b (mod q) }`` has a natural basis
/// ``v1 = (q, 0)``, ``v2 = (r_q, 1)``. After Gauss reduction, the basis vectors
/// ``V1``, ``V2`` are shorter and more orthogonal, enabling efficient enumeration.
///
/// # Gauss-reduced condition
///
/// After reduction: ``|V1| ≤ |V2|`` and ``|dot(V1, V2)| ≤ |V2|² / 2``.
///
/// # Exposed for KAT inspection
///
/// The ``LatticeBasis`` is returned in ``LatticeSieveResult`` so that KATs can verify
/// that the reduced basis vectors are shorter than the original basis and that all
/// enumerated ``(a, b)`` pairs lie in ``L_q``.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LatticeBasis {
    /// First reduced basis vector ``V1 = (v1a, v1b)``.
    pub v1: (i64, i64),
    /// Second reduced basis vector ``V2 = (v2a, v2b)``.
    pub v2: (i64, i64),
    /// The special prime ``q``.
    pub q: u64,
    /// The root ``r_q ∈ [0, q)`` of ``f mod q``.
    pub r_q: u64,
}

impl LatticeBasis {
    /// Construct the initial (unreduced) basis for ``L_q``.
    ///
    /// The initial basis is ``v1 = (q, 0)``, ``v2 = (r_q, 1)``.
    ///
    /// :param q: The special prime.
    /// :param r_q: The root of ``f mod q``.
    /// :returns: The initial (unreduced) basis.
    pub fn initial(q: u64, r_q: u64) -> Self {
        Self { v1: (q as i64, 0), v2: (r_q as i64, 1), q, r_q }
    }

    /// Apply Gauss lattice reduction to find a reduced basis.
    ///
    /// Gauss's algorithm (the 2D analogue of LLL):
    ///
    /// ```text
    /// while |V1| > |V2|: swap V1, V2
    /// V1 = V1 - round(dot(V1, V2) / dot(V2, V2)) * V2
    /// ```
    ///
    /// Terminates when ``|V1| ≤ |V2|`` and ``|dot(V1, V2)| ≤ |V2|² / 2``.
    ///
    /// # Correctness
    ///
    /// The algorithm preserves the lattice: each step replaces ``V1`` with
    /// ``V1 - k·V2`` for integer ``k``, which is a unimodular transformation.
    /// The determinant of the basis matrix is preserved (up to sign), so the
    /// lattice ``L_q`` is unchanged.
    ///
    /// :returns: A new ``LatticeBasis`` with reduced vectors.
    pub fn gauss_reduce(&self) -> Self {
        let mut v1 = self.v1;
        let mut v2 = self.v2;

        // Gauss reduction loop.
        loop {
            // Ensure |V1| ≤ |V2| (swap if needed).
            let norm1_sq = dot(v1, v1);
            let norm2_sq = dot(v2, v2);
            if norm1_sq > norm2_sq {
                std::mem::swap(&mut v1, &mut v2);
            }

            // Compute k = round(dot(V1, V2) / dot(V2, V2)).
            // Use integer arithmetic: k = round(d12 / d22) where d12 = dot(V1, V2),
            // d22 = dot(V2, V2). Round-to-nearest: k = (2*d12 + d22) / (2*d22) using
            // integer division (floor), adjusted for sign.
            let d12 = dot(v1, v2);
            let d22 = dot(v2, v2);

            if d22 == 0 {
                // Degenerate case: V2 is the zero vector. Should not happen for q ≥ 2.
                break;
            }

            // k = round(d12 / d22): round-to-nearest-integer.
            // For exact halves, round away from zero (standard rounding).
            let k = round_div(d12, d22);

            if k == 0 {
                // Already reduced: |dot(V1, V2)| ≤ |V2|² / 2.
                break;
            }

            // V1 ← V1 - k·V2.
            v1 = (v1.0 - k * v2.0, v1.1 - k * v2.1);
        }

        // Final swap to ensure |V1| ≤ |V2|.
        if dot(v1, v1) > dot(v2, v2) {
            std::mem::swap(&mut v1, &mut v2);
        }

        Self { v1, v2, q: self.q, r_q: self.r_q }
    }

    /// Verify that a lattice point ``(a, b) = s·V1 + t·V2`` lies in ``L_q``.
    ///
    /// Checks ``a ≡ r_q·b (mod q)``.
    ///
    /// :param a: The a-coordinate.
    /// :param b: The b-coordinate.
    /// :returns: ``true`` if ``(a, b) ∈ L_q``.
    pub fn in_lattice(&self, a: i64, b: i64) -> bool {
        let q = self.q as i64;
        let r_q = self.r_q as i64;
        let lhs = ((a % q) + q) % q;
        let rhs = ((r_q * b % q) + q) % q;
        lhs == rhs
    }
}

// ─── LatticeSieveResult ───────────────────────────────────────────────────────

/// Per-``(q, r_q)`` output from the lattice sieve.
///
/// Collects the relations found for a single special prime ``q`` and one of its roots
/// ``r_q``. The ``q`` and ``r_q`` fields identify the sieve restriction; every relation
/// in ``relations`` carries ``q`` in its algebraic exponent vector.
#[derive(Debug, Clone)]
pub struct LatticeSieveResult {
    /// The special prime ``q``.
    pub q: u64,
    /// The root ``r_q ∈ [0, q)`` of ``f mod q`` used for this sieve run.
    pub r_q: u64,
    /// The reduced lattice basis used for enumeration.
    pub basis: LatticeBasis,
    /// Relations collected in this sieve run.
    ///
    /// Every relation satisfies ``Relation::verify()`` and carries ``q`` in its
    /// algebraic exponent vector.
    pub relations: Vec<Relation>,
    /// The number of lattice points enumerated (including those filtered out).
    ///
    /// Used for yield comparison: ``relations.len() / enumerated_points`` is the
    /// per-point yield. At toy scale, this is comparable to the special-q sieve yield.
    pub enumerated_points: u64,
}

// ─── lattice_sieve ────────────────────────────────────────────────────────────

/// Run the lattice sieve over the ``(a, b)`` rectangle and return per-``(q, r_q)`` results.
///
/// For each special prime ``q`` in ``[config.q_min, config.q_max]`` that is in the algebraic
/// factor base, and for each root ``r_q`` of ``f mod q``, constructs the lattice
/// ``L_q = { (a, b) : a ≡ r_q·b (mod q) }``, reduces its basis via Gauss's algorithm, and
/// enumerates lattice points in the region ``|a| ≤ A``, ``1 ≤ b ≤ B``. Returns a
/// ``Vec<LatticeSieveResult>`` with one entry per ``(q, r_q)`` pair.
///
/// # Correctness
///
/// Every returned ``Relation`` satisfies ``Relation::verify(poly, fb)``:
///
/// - ``gcd(a, b) = 1`` (coprimality).
/// - Both norms are fully smooth over their factor bases (``cofactor = 1``).
/// - The exponent vectors reconstruct the norm magnitudes.
/// - ``q`` appears in the algebraic exponent vector (guaranteed by the lattice restriction).
///
/// # Determinism
///
/// The output is deterministic for fixed ``(poly, fb, config)``: the sieve is a pure function
/// of its inputs, with no randomness.
///
/// # Toy scale versus cryptographic scale
///
/// At toy scale, the lattice sieve produces the same ``(a, b)`` pairs as the special-q line
/// sieve for the same ``(q, r_q)``. The efficiency advantage of the reduced basis (fewer
/// wasted candidates outside the sieve region) is not visible at toy scale. See the
/// module-level documentation.
///
/// :param poly: The NFS polynomial pair (provides ``f``, ``m`` for norm computation).
/// :param fb: The two-sided factor base.
/// :param config: Sieve region, special-q range, and threshold parameters.
/// :returns: Per-``(q, r_q)`` lattice sieve results.
pub fn lattice_sieve(
    poly: &PolyPair,
    fb: &FactorBase,
    config: &LatticeSieveConfig,
) -> Vec<LatticeSieveResult> {
    // Precompute log2 values for all rational and algebraic primes.
    let rat_logs: Vec<f32> =
        fb.rational_primes.iter().map(|&p| (p as f32).log2()).collect();
    let alg_logs: Vec<f32> =
        fb.algebraic_ideals.iter().map(|ap| (ap.p as f32).log2()).collect();

    // Compute the smoothness threshold.
    let b_alg_f = fb.b_alg.max(2) as f64;
    let threshold = ((config.threshold_scale * b_alg_f.log2()) as f32).max(1.0);

    // Factor-base slices for trial_smooth.
    let rat_fb_slice: Vec<u64> = fb.rational_primes.clone();
    let alg_fb_slice: Vec<u64> = fb.algebraic_ideals.iter().map(|ap| ap.p).collect();

    let mut all_results: Vec<LatticeSieveResult> = Vec::new();

    // Iterate over special primes q in the algebraic factor base within [q_min, q_max].
    let special_ideals: Vec<(u64, u64)> = fb
        .algebraic_ideals
        .iter()
        .filter(|ap| ap.p >= config.q_min && ap.p <= config.q_max)
        .map(|ap| (ap.p, ap.r))
        .collect();

    for (q, r_q) in special_ideals {
        let result = sieve_lattice_for_q(
            poly, fb, config, q, r_q, &rat_logs, &alg_logs, threshold,
            &rat_fb_slice, &alg_fb_slice,
        );
        all_results.push(result);
    }

    all_results
}

/// Run the lattice sieve for a single ``(q, r_q)`` pair.
///
/// Constructs and reduces the lattice basis for ``L_q``, then enumerates lattice points
/// ``(a, b) = s·V1 + t·V2`` in the region ``|a| ≤ A``, ``1 ≤ b ≤ B``. For each point,
/// applies the log-sieve threshold and trial-divides survivors.
///
/// :param poly: The NFS polynomial pair.
/// :param fb: The two-sided factor base.
/// :param config: Sieve configuration.
/// :param q: The special prime.
/// :param r_q: The root of ``f mod q`` used for this run.
/// :param rat_logs: Precomputed ``log2(p)`` for rational primes.
/// :param alg_logs: Precomputed ``log2(p)`` for algebraic ideals.
/// :param threshold: Smoothness threshold.
/// :param rat_fb_slice: Rational factor base as ``Vec<u64>`` for ``trial_smooth``.
/// :param alg_fb_slice: Algebraic factor base primes for ``trial_smooth``.
/// :returns: ``LatticeSieveResult`` for this ``(q, r_q)`` pair.
#[allow(clippy::too_many_arguments)]
fn sieve_lattice_for_q(
    poly: &PolyPair,
    fb: &FactorBase,
    config: &LatticeSieveConfig,
    q: u64,
    r_q: u64,
    rat_logs: &[f32],
    alg_logs: &[f32],
    threshold: f32,
    rat_fb_slice: &[u64],
    alg_fb_slice: &[u64],
) -> LatticeSieveResult {
    let a_bound = config.a_bound as i64;
    let b_bound = config.b_bound as i64;

    // ── Step 1: Construct and reduce the lattice basis ────────────────────────
    let initial_basis = LatticeBasis::initial(q, r_q);
    let basis = initial_basis.gauss_reduce();
    let (v1a, v1b) = basis.v1;
    let (v2a, v2b) = basis.v2;

    // ── Step 2: Determine enumeration bounds ──────────────────────────────────
    //
    // We want to enumerate all (s, t) such that (a, b) = s·V1 + t·V2 satisfies
    // |a| ≤ A and 1 ≤ b ≤ B.
    //
    // The b-coordinate is: b = s·v1b + t·v2b.
    // The a-coordinate is: a = s·v1a + t·v2a.
    //
    // To cover the region, we bound s and t conservatively:
    //   |s·v1b + t·v2b| ≤ B  →  |s| ≤ (B + |t·v2b|) / |v1b|  (if v1b ≠ 0)
    //   |s·v1a + t·v2a| ≤ A  →  |s| ≤ (A + |t·v2a|) / |v1a|  (if v1a ≠ 0)
    //
    // A safe conservative bound: enumerate s, t in [-S, S] where S is chosen so
    // that any (a, b) in the target region is covered.
    //
    // The lattice has determinant q (since det([v1, v2]) = q·1 - 0·r_q = q for the
    // initial basis, and Gauss reduction preserves the determinant). The reduced basis
    // vectors have length ≈ √q. The number of lattice points in the region is ≈ A·B / q.
    //
    // Conservative bound: S = T = ceil(max(A, B) / min_basis_length) + 2.
    // For safety, use S = T = A + B + 2 (covers all cases at toy scale).
    let s_bound = a_bound + b_bound + 2;
    let t_bound = a_bound + b_bound + 2;

    // ── Step 3: Enumerate lattice points and sieve ────────────────────────────
    let mut relations: Vec<Relation> = Vec::new();
    let mut enumerated_points: u64 = 0;

    // For each (s, t) in [-S, S] × [-T, T], compute (a, b) = s·V1 + t·V2.
    for t in -t_bound..=t_bound {
        for s in -s_bound..=s_bound {
            let a = s * v1a + t * v2a;
            let b = s * v1b + t * v2b;

            // Filter: must be in the sieve region.
            if a < -a_bound || a > a_bound {
                continue;
            }
            if b < 1 || b > b_bound {
                continue;
            }

            enumerated_points += 1;

            let a_big = BigInt::from(a);
            let b_big = BigInt::from(b);

            // ── Step 3a: Check gcd(a, b) = 1 ─────────────────────────────────
            if !is_coprime_i64(a, b) {
                continue;
            }

            // ── Step 3b: Log-sieve threshold ──────────────────────────────────
            //
            // Compute the approximate log-smoothness score for this (a, b) by summing
            // log2(p) for each prime p in the factor bases that divides the respective norm.
            //
            // Rational side: p | N_rat(a, b) iff a ≡ b·m (mod p).
            // Algebraic side: p | N_alg(a, b) iff a ≡ r·b (mod p) for some root r.
            let mut score: f32 = 0.0;

            // Rational side.
            for (pi, &p) in fb.rational_primes.iter().enumerate() {
                let p_i = p as i64;
                // a ≡ b·m (mod p): check (a - b*m) mod p == 0.
                // Use i64 arithmetic for small p.
                let m_i64 = bigint_to_i64_mod(&poly.m, p);
                let bm = (b as i128 * m_i64 as i128).rem_euclid(p as i128) as i64;
                let a_mod = a.rem_euclid(p_i);
                if a_mod == bm {
                    score += rat_logs[pi];
                }
            }

            // Algebraic side.
            for (ai, ap) in fb.algebraic_ideals.iter().enumerate() {
                let p_i = ap.p as i64;
                let r_i = ap.r as i64;
                // a ≡ r·b (mod p).
                let rb = (r_i as i128 * b as i128).rem_euclid(ap.p as i128) as i64;
                let a_mod = a.rem_euclid(p_i);
                if a_mod == rb {
                    score += alg_logs[ai];
                }
            }

            if score < threshold {
                continue;
            }

            // ── Step 3c: Compute norms ────────────────────────────────────────
            let rat_norm = rational_norm(&a_big, &b_big, &poly.m);
            let alg_norm = algebraic_norm(&a_big, &b_big, &poly.f);

            if rat_norm.is_zero() || alg_norm.is_zero() {
                continue;
            }

            let rat_uint = match norm_to_uint(&rat_norm) {
                Ok(u) => u,
                Err(_) => continue,
            };
            let alg_uint = match norm_to_uint(&alg_norm) {
                Ok(u) => u,
                Err(_) => continue,
            };

            // ── Step 3d: Trial-smooth both norms ──────────────────────────────
            let rat_witness = trial_smooth(&rat_uint, rat_fb_slice);
            if !rat_witness.is_smooth() {
                continue;
            }
            let alg_witness = trial_smooth(&alg_uint, alg_fb_slice);
            if !alg_witness.is_smooth() {
                continue;
            }

            // ── Step 3e: Construct the Relation ───────────────────────────────
            let rational_sign = norm_sign(&rat_norm);
            if let Some(rel) = Relation::new(
                a_big,
                b_big,
                &rat_witness,
                &alg_witness,
                rational_sign,
                fb,
            ) {
                relations.push(rel);
            }
        }
    }

    LatticeSieveResult { q, r_q, basis, relations, enumerated_points }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Compute the integer dot product of two 2D vectors.
fn dot(u: (i64, i64), v: (i64, i64)) -> i64 {
    u.0 * v.0 + u.1 * v.1
}

/// Round ``numerator / denominator`` to the nearest integer (round half away from zero).
///
/// Used in Gauss reduction: ``k = round(dot(V1, V2) / dot(V2, V2))``.
fn round_div(numerator: i64, denominator: i64) -> i64 {
    // Standard rounding: add half the denominator (with sign) before dividing.
    // For positive denominator: k = (2*n + d) / (2*d) using floor division.
    // Handle sign carefully.
    if denominator == 0 {
        return 0;
    }
    // Use i128 to avoid overflow for large lattice vectors.
    let n = numerator as i128;
    let d = denominator as i128;
    // Round half away from zero: (2n + sign(n)*d) / (2d).
    // Equivalently: floor((n + d/2) / d) for d > 0, adjusted for sign.
    // Simplest correct formula: (2n + d) / (2d) using truncating division, then adjust.
    // Use: k = (n + d/2) / d for d > 0 (floor), negated for d < 0.
    let (n_adj, d_abs) = if d > 0 { (n, d) } else { (-n, -d) };
    // Round half away from zero for positive d_abs:
    // k = floor((n_adj + d_abs/2) / d_abs) for n_adj >= 0
    // k = -floor((-n_adj + d_abs/2) / d_abs) for n_adj < 0
    let k = if n_adj >= 0 {
        (n_adj + d_abs / 2) / d_abs
    } else {
        -((-n_adj + d_abs / 2) / d_abs)
    };
    // Adjust sign back if d was negative.
    let k = if denominator < 0 { -k } else { k };
    k as i64
}

/// Check if ``gcd(|a|, |b|) = 1`` for ``a: i64``, ``b: i64``.
fn is_coprime_i64(a: i64, b: i64) -> bool {
    let a_abs = a.unsigned_abs();
    let b_abs = b.unsigned_abs();
    gcd_u64(a_abs, b_abs) == 1
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

/// Reduce a ``BigInt`` modulo ``p`` and return the result as ``i64``.
///
/// Returns the canonical representative in ``[0, p)``.
fn bigint_to_i64_mod(n: &BigInt, p: u64) -> i64 {
    let p_big = BigInt::from(p);
    let r = n % &p_big;
    let r = if r < BigInt::zero() { r + &p_big } else { r };
    use num_traits::ToPrimitive;
    r.to_i64().expect("residue mod p fits in i64 for p ≤ u64::MAX")
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── LatticeBasis ─────────────────────────────────────────────────────────

    #[test]
    fn initial_basis_vectors_are_correct() {
        // For q=7, r_q=5: v1=(7,0), v2=(5,1).
        let basis = LatticeBasis::initial(7, 5);
        assert_eq!(basis.v1, (7, 0));
        assert_eq!(basis.v2, (5, 1));
    }

    #[test]
    fn initial_basis_points_in_lattice() {
        // v1 = (7, 0): 7 ≡ 5·0 (mod 7) → 7 ≡ 0 (mod 7) ✓
        // v2 = (5, 1): 5 ≡ 5·1 (mod 7) ✓
        let basis = LatticeBasis::initial(7, 5);
        assert!(basis.in_lattice(7, 0), "v1 = (7, 0) should be in L_q");
        assert!(basis.in_lattice(5, 1), "v2 = (5, 1) should be in L_q");
    }

    #[test]
    fn gauss_reduce_produces_shorter_vectors() {
        // For q=7, r_q=5: initial basis v1=(7,0), v2=(5,1).
        // After reduction, |V1|² + |V2|² ≤ |v1|² + |v2|² (Gauss reduction shortens vectors).
        let initial = LatticeBasis::initial(7, 5);
        let reduced = initial.gauss_reduce();

        let initial_sum_sq = dot(initial.v1, initial.v1) + dot(initial.v2, initial.v2);
        let reduced_sum_sq = dot(reduced.v1, reduced.v1) + dot(reduced.v2, reduced.v2);

        assert!(
            reduced_sum_sq <= initial_sum_sq,
            "Gauss reduction should not increase the sum of squared norms: \
             initial={initial_sum_sq}, reduced={reduced_sum_sq}"
        );
    }

    #[test]
    fn gauss_reduce_preserves_lattice_membership() {
        // The reduced basis vectors must still be in L_q.
        let initial = LatticeBasis::initial(7, 5);
        let reduced = initial.gauss_reduce();

        assert!(
            reduced.in_lattice(reduced.v1.0, reduced.v1.1),
            "reduced V1 = {:?} should be in L_q (q={}, r_q={})",
            reduced.v1, reduced.q, reduced.r_q
        );
        assert!(
            reduced.in_lattice(reduced.v2.0, reduced.v2.1),
            "reduced V2 = {:?} should be in L_q (q={}, r_q={})",
            reduced.v2, reduced.q, reduced.r_q
        );
    }

    #[test]
    fn gauss_reduce_satisfies_reduced_condition() {
        // After reduction: |V1| ≤ |V2|.
        let initial = LatticeBasis::initial(7, 5);
        let reduced = initial.gauss_reduce();

        let norm1_sq = dot(reduced.v1, reduced.v1);
        let norm2_sq = dot(reduced.v2, reduced.v2);

        assert!(
            norm1_sq <= norm2_sq,
            "Gauss-reduced condition: |V1| ≤ |V2|, but |V1|²={norm1_sq} > |V2|²={norm2_sq}"
        );
    }

    #[test]
    fn round_div_basic() {
        assert_eq!(round_div(3, 2), 2);   // 3/2 = 1.5 → rounds to 2
        assert_eq!(round_div(1, 2), 1);   // 1/2 = 0.5 → rounds to 1 (away from zero)
        assert_eq!(round_div(-1, 2), -1); // -1/2 = -0.5 → rounds to -1 (away from zero)
        assert_eq!(round_div(4, 3), 1);   // 4/3 ≈ 1.33 → rounds to 1
        assert_eq!(round_div(5, 3), 2);   // 5/3 ≈ 1.67 → rounds to 2
        assert_eq!(round_div(0, 5), 0);   // 0/5 = 0
    }

    #[test]
    fn in_lattice_checks_congruence() {
        let basis = LatticeBasis::initial(7, 5);
        // (a=5, b=1): 5 ≡ 5·1 (mod 7) ✓
        assert!(basis.in_lattice(5, 1));
        // (a=12, b=1): 12 mod 7 = 5 ≡ 5·1 (mod 7) ✓
        assert!(basis.in_lattice(12, 1));
        // (a=3, b=1): 3 ≢ 5 (mod 7) ✗
        assert!(!basis.in_lattice(3, 1));
        // (a=10, b=2): 10 mod 7 = 3, 5·2 mod 7 = 3 ✓
        assert!(basis.in_lattice(10, 2));
    }

    #[test]
    fn gcd_u64_cases() {
        assert_eq!(gcd_u64(7, 13), 1);
        assert_eq!(gcd_u64(12, 8), 4);
        assert_eq!(gcd_u64(0, 5), 5);
        assert_eq!(gcd_u64(5, 0), 5);
    }

    #[test]
    fn is_coprime_i64_cases() {
        assert!(is_coprime_i64(5, 7));
        assert!(is_coprime_i64(-5, 7));
        assert!(!is_coprime_i64(4, 2));
        assert!(!is_coprime_i64(-6, 4));
    }
}
