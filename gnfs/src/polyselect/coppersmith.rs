//! Coppersmith multiple-polynomial method for NFS polynomial selection.
//!
//! Coppersmith's method generates several algebraic-side polynomials `f_1, f_2, ..., f_k`
//! that all share the same rational-side polynomial `g = x − m`. The key insight is that
//! for NFS, different algebraic polynomials can be used in different parts of the sieve
//! region, potentially improving the overall relation yield.
//!
//! # Mathematical construction
//!
//! Given a seed pair `(f_0, g)` with `f_0(m) ≡ 0 (mod N)` and `g(x) = x − m`, each
//! variant is produced by the Kleinjung rotation:
//!
//! ```text
//! f_i(x) = f_0(x) + (j_i · x + k_i) · g(x)
//! ```
//!
//! This preserves the root because `g(m) = 0`:
//!
//! ```text
//! f_i(m) = f_0(m) + (j_i · m + k_i) · g(m) = 0 + 0 = 0  (mod N)
//! ```
//!
//! so `PolyPair::verify()` holds for every variant. The rotation parameters `(j_i, k_i)`
//! are chosen systematically — see [`coppersmith_polys`] for the spiral pattern used here.
//!
//! # **Toy scale versus cryptographic scale**
//!
//! In production NFS (e.g., RSA-768), using different polynomials for different sieve
//! regions measurably improves the relation yield. At toy scale (60–100 bit N), the sieve
//! region is too small for this effect to manifest — the improvement is present in the
//! construction but invisible in the numbers. This is a **demonstration-fidelity**
//! implementation: the mathematical content is complete, but the engineering payoff
//! requires cryptographic-scale N to observe.
//!
//! Concretely: at toy scale the Murphy-E score improvement from multi-poly is typically
//! small (< 2×), and the absolute score values are not meaningful (see the `murphy` module
//! for the full toy-scale caveat). Downstream consumers should treat the score as an
//! ordinal ranking, not a cardinal measure.
//!
//! # References
//!
//! - Coppersmith, D. (1993). *Modifications to the number field sieve*. Journal of
//!   Cryptology, 6(3), 169–180.
//! - Bai, S., Bouvier, C., Kruppa, A., Zimmermann, P. (2014). *Better polynomials for
//!   GNFS*. Mathematics of Computation.

use super::{murphy::score, root_sieve::rotate, PolyGenerator, PolyPair};

// ─── CoppersmithConfig ───────────────────────────────────────────────────────

/// Configuration for Coppersmith multi-poly generation.
///
/// Controls how many polynomial variants are generated and the step size used
/// to produce the rotation parameters `(j_i, k_i)` for each variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoppersmithConfig {
    /// Number of polynomial variants to generate (including the seed at `i = 0`).
    pub num_polys: usize,
    /// The rotation step size for generating variants.
    ///
    /// Variant `i` uses rotation parameters derived from `i * step` — see
    /// [`coppersmith_polys`] for the exact spiral pattern.
    pub step: i64,
}

impl Default for CoppersmithConfig {
    fn default() -> Self {
        Self { num_polys: 5, step: 1 }
    }
}

// ─── rotation parameter schedule ─────────────────────────────────────────────

/// Compute the `(j, k)` rotation parameters for variant `i` with the given step size.
///
/// The schedule uses a symmetric spiral pattern that alternates between translation
/// (`k`-axis) and rotation (`j`-axis) moves, keeping the parameters small:
///
/// ```text
/// i = 0: (j=0,  k=0)          — the seed itself
/// i = 1: (j=0,  k=+step)      — translate by +step
/// i = 2: (j=0,  k=-step)      — translate by -step
/// i = 3: (j=+step, k=0)       — rotate by +step
/// i = 4: (j=-step, k=0)       — rotate by -step
/// i = 5: (j=0,  k=+2*step)    — translate by +2*step
/// i = 6: (j=0,  k=-2*step)    — translate by -2*step
/// i = 7: (j=+2*step, k=0)     — rotate by +2*step
/// i = 8: (j=-2*step, k=0)     — rotate by -2*step
/// ...
/// ```
///
/// For `i ≥ 1`, the group index `g = (i - 1) / 4` and the position within the group
/// `p = (i - 1) % 4` determine the parameters:
///
/// - `p = 0`: `(0, +(g+1)*step)`
/// - `p = 1`: `(0, -(g+1)*step)`
/// - `p = 2`: `(+(g+1)*step, 0)`
/// - `p = 3`: `(-(g+1)*step, 0)`
///
/// :param i: Variant index (0-based).
/// :param step: Step size.
/// :returns: `(j, k)` rotation parameters for variant `i`.
fn rotation_params(i: usize, step: i64) -> (i64, i64) {
    if i == 0 {
        return (0, 0);
    }
    let g = ((i - 1) / 4) as i64 + 1; // group index, 1-based
    let p = (i - 1) % 4;
    let magnitude = g * step;
    match p {
        0 => (0, magnitude),
        1 => (0, -magnitude),
        2 => (magnitude, 0),
        3 => (-magnitude, 0),
        _ => unreachable!(),
    }
}

// ─── coppersmith_polys ───────────────────────────────────────────────────────

/// Generate multiple polynomial pairs via Coppersmith's method.
///
/// Generates `config.num_polys` polynomial pairs by applying systematic rotations
/// to the seed pair. Each generated pair satisfies `PolyPair::verify()`.
///
/// The rotation schedule is a symmetric spiral (see [`rotation_params`]):
/// variant 0 is the seed itself, and subsequent variants alternate between
/// positive/negative translations and rotations at increasing magnitudes.
///
/// **Demonstration fidelity:** the mathematical construction is present, but the
/// yield improvement (using different polynomials for different sieve regions)
/// is under-exposed at toy scale — see the module-level annotation.
///
/// :param seed: The seed polynomial pair (typically from base-m or root sieve).
/// :param config: Generation configuration.
/// :returns: A `Vec` of exactly `config.num_polys` polynomial pairs, each satisfying
///   `PolyPair::verify()`. The first element is always the seed (rotation `(0, 0)`).
pub fn coppersmith_polys(seed: &PolyPair, config: &CoppersmithConfig) -> Vec<PolyPair> {
    (0..config.num_polys)
        .map(|i| {
            let (j, k) = rotation_params(i, config.step);
            rotate(seed, j, k)
        })
        .collect()
}

// ─── coppersmith_best ────────────────────────────────────────────────────────

/// Select the best polynomial from a Coppersmith multi-poly set.
///
/// Generates `config.num_polys` candidates via [`coppersmith_polys`], scores each
/// with Murphy-E, and returns the highest-scoring pair.
///
/// Because the seed itself is always included (as variant 0), the returned pair
/// is always at least as good as the seed under the Murphy-E heuristic.
///
/// :param seed: The seed polynomial pair.
/// :param config: Generation configuration.
/// :returns: The highest-scoring `PolyPair` from the generated set.
pub fn coppersmith_best(seed: &PolyPair, config: &CoppersmithConfig) -> PolyPair {
    coppersmith_polys(seed, config)
        .into_iter()
        .max_by(|a, b| score(a).partial_cmp(&score(b)).expect("Murphy-E score is finite"))
        .expect("coppersmith_polys returns at least one pair when num_polys >= 1")
}

// ─── CoppersmithGenerator ────────────────────────────────────────────────────

/// Generator that yields all Coppersmith multi-poly candidates.
///
/// Implements [`PolyGenerator`] so that the Coppersmith method fits into the common
/// score-and-rank pipeline alongside base-m and root-sieve generators. Unlike
/// [`coppersmith_best`], which returns only the best candidate, `generate` yields
/// all `config.num_polys` candidates in spiral order.
///
/// Callers can use `.max_by(|a, b| score(a).partial_cmp(&score(b)).unwrap())` to
/// recover the best, or `.take(limit)` for early termination.
pub struct CoppersmithGenerator {
    /// The seed polynomial pair to rotate.
    pub seed: PolyPair,
    /// Generation configuration.
    pub config: CoppersmithConfig,
}

impl PolyGenerator for CoppersmithGenerator {
    /// Yield all Coppersmith multi-poly candidates in spiral order.
    ///
    /// Yields exactly `config.num_polys` pairs. The first pair is always the seed
    /// (rotation `(0, 0)`); subsequent pairs use the spiral schedule from
    /// [`rotation_params`].
    fn generate(&self) -> impl Iterator<Item = PolyPair> {
        let seed = self.seed.clone();
        let num_polys = self.config.num_polys;
        let step = self.config.step;

        (0..num_polys).map(move |i| {
            let (j, k) = rotation_params(i, step);
            rotate(&seed, j, k)
        })
    }
}

// ─── unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polyselect::base_m::select_base_m;
    use num_bigint::BigInt;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn rotation_params_seed_is_identity() {
        // Variant 0 must always be the identity rotation.
        assert_eq!(rotation_params(0, 1), (0, 0));
        assert_eq!(rotation_params(0, 5), (0, 0));
    }

    #[test]
    fn rotation_params_spiral_pattern() {
        // Verify the spiral schedule for step=1.
        assert_eq!(rotation_params(1, 1), (0, 1));
        assert_eq!(rotation_params(2, 1), (0, -1));
        assert_eq!(rotation_params(3, 1), (1, 0));
        assert_eq!(rotation_params(4, 1), (-1, 0));
        assert_eq!(rotation_params(5, 1), (0, 2));
        assert_eq!(rotation_params(6, 1), (0, -2));
        assert_eq!(rotation_params(7, 1), (2, 0));
        assert_eq!(rotation_params(8, 1), (-2, 0));
    }

    #[test]
    fn coppersmith_polys_count() {
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let config = CoppersmithConfig { num_polys: 7, step: 1 };
        let polys = coppersmith_polys(&seed, &config);
        assert_eq!(polys.len(), 7);
    }

    #[test]
    fn coppersmith_polys_all_verify() {
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let config = CoppersmithConfig::default();
        for (i, poly) in coppersmith_polys(&seed, &config).iter().enumerate() {
            poly.verify().unwrap_or_else(|e| panic!("variant {i} failed verify: {e}"));
        }
    }

    #[test]
    fn coppersmith_best_ge_seed() {
        let n = bi(1009 * 1013);
        let seed = select_base_m(&n, 3);
        let seed_score = score(&seed);
        let config = CoppersmithConfig::default();
        let best = coppersmith_best(&seed, &config);
        let best_score = score(&best);
        assert!(
            best_score >= seed_score,
            "coppersmith_best should return score ≥ seed: seed={seed_score:.6e}, best={best_score:.6e}"
        );
    }
}
