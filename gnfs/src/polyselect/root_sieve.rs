//! Kleinjung-style rotation and translation search for NFS polynomial selection (root sieve).
//!
//! Given a seed polynomial pair `(f, g)` from base-m selection, this module searches over
//! a grid of rotation parameters `(j, k)` to find a better-scoring pair. The rotation is:
//!
//! ```text
//! f'(x) = f(x) + (j·x + k)·g(x)
//! g'(x) = g(x)   (unchanged)
//! ```
//!
//! This preserves the shared root `m` mod `n` because:
//!
//! ```text
//! f'(m) = f(m) + (j·m + k)·g(m) = 0 + (j·m + k)·0 = 0  (mod n)
//! ```
//!
//! so `PolyPair::verify()` holds for every rotated pair. The search scores each candidate
//! with Murphy-E scoring and returns the best-scoring one.
//!
//! # Rotation arithmetic
//!
//! For `g(x) = x − m` (coefficients `[−m, 1]`):
//!
//! ```text
//! (j·x + k)·g(x) = (j·x + k)·(x − m)
//!                 = j·x² + (k − j·m)·x − k·m
//! ```
//!
//! This adds to `f`'s coefficients:
//!
//! - `f'[0] += −k·m`
//! - `f'[1] += k − j·m`
//! - `f'[2] += j`
//!
//! For degree `d ≥ 3` (always the case for toy N via base-m), adding a degree-2 term
//! does not change the degree of `f`.
//!
//! # Toy scale versus cryptographic scale
//!
//! At toy scale (N < 2^60), the Murphy-E improvement from rotation is small and may not
//! manifest as a strict improvement for every seed. The search is correct in expectation
//! (the ordering is right), but the absolute gain is under-exposed. At cryptographic scale
//! (RSA-768+), rotation search is essential and produces measurable improvements.
//!
//! # References
//!
//! - Kleinjung, T. (2006). *Polynomial selection*. CADO workshop on integer factorization.
//! - Bai, S., Bouvier, C., Kruppa, A., Zimmermann, P. (2014). *Better polynomials for GNFS*.
//!   Mathematics of Computation.

use num_bigint::BigInt;
use num_traits::One;
use shared_numfield::IntPoly;

use super::{murphy::score, PolyGenerator, PolyPair};

// ─── RootSieveConfig ─────────────────────────────────────────────────────────

/// Configuration for the Kleinjung rotation search.
///
/// The search covers the grid `j ∈ [−j_range, j_range]` × `k ∈ [−k_range, k_range]`,
/// yielding `(2·j_range + 1) × (2·k_range + 1)` candidate pairs in total.
///
/// Larger ranges find better polynomials at the cost of more scoring calls. For toy N,
/// the default `j_range = k_range = 10` (441 candidates) is sufficient. For cryptographic
/// N, ranges of 100–1000 are typical.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootSieveConfig {
    /// Half-width of the `j` (rotation coefficient) search range: `j ∈ [−j_range, j_range]`.
    pub j_range: i64,
    /// Half-width of the `k` (translation) search range: `k ∈ [−k_range, k_range]`.
    pub k_range: i64,
}

impl Default for RootSieveConfig {
    fn default() -> Self {
        Self { j_range: 10, k_range: 10 }
    }
}

// ─── rotation ────────────────────────────────────────────────────────────────

/// Apply the Kleinjung rotation `f' = f + (j·x + k)·g` to a polynomial pair.
///
/// Given `g(x) = x − m`, the rotation adds:
///
/// - `f'[0] += −k·m`
/// - `f'[1] += k − j·m`
/// - `f'[2] += j`
///
/// The result satisfies `f'(m) ≡ 0 (mod n)` because `g(m) = 0`. The degree of `f'`
/// equals the degree of `f` when `d ≥ 3` (the `j` contribution lands at index 2, below
/// the leading term).
///
/// :param seed: The seed polynomial pair. Must have `f.degree() ≥ 3`.
/// :param j: The rotation coefficient (multiplies `x·g`).
/// :param k: The translation coefficient (multiplies `g`).
/// :returns: A new `PolyPair` with the rotated `f'` and the same `g`, `m`, `n`.
pub fn rotate(seed: &PolyPair, j: i64, k: i64) -> PolyPair {
    let m = &seed.m;
    let n = &seed.n;

    // Clone f's coefficients, extending to at least degree 3 if needed.
    let mut coeffs: Vec<BigInt> = seed.f.coeffs.clone();
    // Ensure we have at least 3 coefficients (indices 0, 1, 2) to apply the rotation.
    while coeffs.len() < 3 {
        coeffs.push(BigInt::from(0i64));
    }

    let j_big = BigInt::from(j);
    let k_big = BigInt::from(k);

    // f'[0] += -k * m
    coeffs[0] -= &k_big * m;
    // f'[1] += k - j * m
    coeffs[1] += &k_big - &j_big * m;
    // f'[2] += j
    coeffs[2] += &j_big;

    // Trim trailing zeros to keep the degree consistent, but preserve at least d+1 coeffs
    // so the degree does not accidentally drop below d for the j=0, k=0 case.
    // For j ≠ 0 and d = 2, the degree would increase to 3 — but base-m always gives d ≥ 3,
    // so this path is unreachable in practice.
    let f_prime = IntPoly::from_coeffs(coeffs);

    // g is unchanged: g(x) = x - m.
    let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);

    PolyPair::new(f_prime, g, m.clone(), n.clone())
}

// ─── root_sieve ──────────────────────────────────────────────────────────────

/// Search for a better polynomial pair by Kleinjung rotation.
///
/// Searches over `(j, k) ∈ [−j_range, j_range] × [−k_range, k_range]`, scoring each
/// rotated pair with Murphy-E and returning the best-scoring one. The seed pair itself
/// is included in the search (at `j = 0, k = 0`), so the returned pair is always at
/// least as good as the seed.
///
/// The search is deterministic for a fixed seed and config: the grid is traversed in
/// row-major order `(j from −j_range to j_range, k from −k_range to k_range)`, and ties
/// are broken by keeping the first maximum found.
///
/// **Science↔engineering disconnect:** At toy scale (N < 2^60), the Murphy-E improvement
/// from rotation may be small. The returned pair is the best according to the Murphy-E
/// heuristic, which is an ordinal ranking at toy scale. See the module-level note.
///
/// :param seed: The seed polynomial pair (typically from base-m selection).
/// :param config: Search grid configuration.
/// :returns: The best-scoring `PolyPair` found in the grid.
pub fn root_sieve(seed: &PolyPair, config: &RootSieveConfig) -> PolyPair {
    let mut best_pair = rotate(seed, 0, 0); // j=0, k=0 is the seed itself
    let mut best_score = score(&best_pair);

    for j in -config.j_range..=config.j_range {
        for k in -config.k_range..=config.k_range {
            if j == 0 && k == 0 {
                continue; // already scored above
            }
            let candidate = rotate(seed, j, k);
            let s = score(&candidate);
            if s > best_score {
                best_score = s;
                best_pair = candidate;
            }
        }
    }

    best_pair
}

// ─── RootSieveGenerator ──────────────────────────────────────────────────────

/// Generator that yields all rotation candidates in the search grid.
///
/// Implements [`PolyGenerator`] so that the root sieve fits into the common score-and-rank
/// pipeline alongside base-m and Coppersmith generators. Unlike [`root_sieve`], which
/// returns only the best candidate, `RootSieveGenerator::generate` yields all
/// `(2·j_range + 1) × (2·k_range + 1)` candidates in row-major order.
///
/// Callers can use `.max_by(|a, b| score(a).partial_cmp(&score(b)).unwrap())` to recover
/// the best, or `.take(limit)` for early termination.
pub struct RootSieveGenerator {
    /// The seed polynomial pair to rotate.
    pub seed: PolyPair,
    /// Search grid configuration.
    pub config: RootSieveConfig,
}

impl PolyGenerator for RootSieveGenerator {
    /// Yield all rotation candidates in the search grid.
    ///
    /// Iterates over `j ∈ [−j_range, j_range]` (outer) and `k ∈ [−k_range, k_range]`
    /// (inner) in row-major order, yielding the rotated pair for each `(j, k)`.
    /// The `(j=0, k=0)` candidate (the seed itself) is included.
    fn generate(&self) -> impl Iterator<Item = PolyPair> {
        let seed = self.seed.clone();
        let j_range = self.config.j_range;
        let k_range = self.config.k_range;

        // Collect all (j, k) pairs into a Vec so we can return a concrete iterator.
        // The grid is small (≤ (2*20+1)^2 = 1681 for default config), so allocation is fine.
        let pairs: Vec<(i64, i64)> = (-j_range..=j_range)
            .flat_map(move |j| (-k_range..=k_range).map(move |k| (j, k)))
            .collect();

        pairs.into_iter().map(move |(j, k)| rotate(&seed, j, k))
    }
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polyselect::base_m::select_base_m;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn rotate_identity_is_seed() {
        // Rotating with j=0, k=0 should reproduce the seed (same coefficients).
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let rotated = rotate(&seed, 0, 0);
        assert_eq!(rotated.f.coeffs, seed.f.coeffs);
        assert_eq!(rotated.verify(), Ok(()));
    }

    #[test]
    fn rotate_preserves_root() {
        // Every rotation must satisfy PolyPair::verify().
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        for j in -3i64..=3 {
            for k in -3i64..=3 {
                let rotated = rotate(&seed, j, k);
                assert_eq!(
                    rotated.verify(),
                    Ok(()),
                    "rotate({j}, {k}) failed verify"
                );
            }
        }
    }

    #[test]
    fn rotate_coefficient_arithmetic() {
        // Verify the coefficient update formula directly.
        // For j=1, k=2, m=seed.m:
        //   f'[0] = f[0] - k*m = f[0] - 2*m
        //   f'[1] = f[1] + k - j*m = f[1] + 2 - m
        //   f'[2] = f[2] + j = f[2] + 1
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let m = &seed.m;
        let j = 1i64;
        let k = 2i64;
        let rotated = rotate(&seed, j, k);

        let expected_0 = &seed.f.coeffs[0] - bi(k) * m;
        let expected_1 = &seed.f.coeffs[1] + bi(k) - bi(j) * m;
        let expected_2 = &seed.f.coeffs[2] + bi(j);

        assert_eq!(rotated.f.coeffs[0], expected_0, "coeff[0] mismatch");
        assert_eq!(rotated.f.coeffs[1], expected_1, "coeff[1] mismatch");
        assert_eq!(rotated.f.coeffs[2], expected_2, "coeff[2] mismatch");
        // Higher coefficients are unchanged.
        assert_eq!(rotated.f.coeffs[3], seed.f.coeffs[3], "coeff[3] should be unchanged");
    }

    #[test]
    fn root_sieve_includes_seed() {
        // root_sieve with j_range=0, k_range=0 must return the seed itself.
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let config = RootSieveConfig { j_range: 0, k_range: 0 };
        let result = root_sieve(&seed, &config);
        assert_eq!(result.f.coeffs, seed.f.coeffs);
    }

    #[test]
    fn generator_count() {
        // RootSieveGenerator should yield exactly (2*j+1)*(2*k+1) candidates.
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let config = RootSieveConfig { j_range: 3, k_range: 4 };
        let generator = RootSieveGenerator { seed, config: config.clone() };
        let count = generator.generate().count();
        let expected = (2 * config.j_range + 1) as usize * (2 * config.k_range + 1) as usize;
        assert_eq!(count, expected, "generator should yield {expected} candidates");
    }
}
