//! Index-calculus strategy substrate: factor base, prime-order subgroup, relation type.
//!
//! This module defines the three core types for the index-calculus pipeline:
//! - [`FbPoint`] — a curve point with its stable factor-base index.
//! - [`IndexCalcStrategy`] — the factor base + prime-order subgroup + decomposition arity,
//!   bound to a curve.
//! - [`Relation`] — one index-calculus relation: a decomposition of `R = a·G + b·Q` over
//!   the factor base as a sparse exponent vector over `F_ℓ`.
//!
//! # Factor-base enumeration
//!
//! Points are enumerated by ascending x-coordinate: for `x = 0, 1, 2, …`, compute
//! `v = x³ + ax + b mod p`, check if `v` is a quadratic residue (Legendre symbol = 1),
//! lift the **canonical** y (the smaller of `±y` by `to_uint()` comparison), and record
//! the point. The first `fb_size` such points form the factor base.
//!
//! # Prime-order subgroup
//!
//! The toy fixture has group order `n = 60 = 2²·3·5`. The index-calculus linear algebra
//! requires the relation exponents to live in a field `Z/ℓℤ`; this requires ℓ to be prime
//! and to divide n. The ratified choice is `ℓ = 5` (the largest prime factor of n = 60).
//! The ℓ-order subgroup generator is `G_ℓ = (n/ℓ)·G = 12·G`.
//!
//! # Scale limits
//!
//! The toy fixture is `F_47`, `n = 60`, `ℓ = 5`, `m = 2`. The algorithms are
//! mechanism-correct; native enumeration of the factor base does not scale past toy parameters,
//! and the asymptotic index-calculus win (which needs `E(F_{p^n})`) is not observable over
//! `E(F_p)` — the mechanism is demonstrated, not the advantage.

use crypto_bigint::Uint;

use shared_field::{Fp, FpNaive4 as FpNaive};

use crate::curve::{AffinePoint, Curve};
use crate::index_calculus::IndexCalcError;
use crate::semaev::semaev_toy;

// ─── toy-fixture constants ────────────────────────────────────────────────────

/// The prime-order subgroup modulus for the toy fixture: ℓ = 5 (the largest prime factor
/// of n = 60 = 2²·3·5). Relation exponents and the linear algebra live over `F_ℓ = F_5`.
pub const TOY_ELL: u64 = 5;

/// Factor-base size for the toy fixture. Over-determinable at m = 2: the collection target
/// is ≥ FB_SIZE + 1 = 7 relations; the toy curve's ~30 affine x-coordinates with a QR
/// supply ample candidates.
pub const TOY_FB_SIZE: usize = 6;

/// Decomposition arity: a point decomposes as a sum of `m` factor-base points (via S_{m+1}).
/// m = 2 uses S_3 (the smallest non-trivial Semaev polynomial, already built and KAT-covered
/// in `rho::semaev`).
pub const TOY_M: usize = 2;

// ─── FbPoint ─────────────────────────────────────────────────────────────────

/// A factor-base point: a curve point plus its stable factor-base index.
///
/// The index is the point's position in the enumeration order (ascending x-coordinate,
/// canonical y). It is stable across the pipeline — relations reference factor-base points
/// by this index, and the linear-algebra adapter uses it as the column index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FbPoint {
    /// Index into the factor base (the relation exponent-vector column).
    pub index: usize,
    /// The curve point (a finite affine point on the frozen `Curve`).
    pub point: AffinePoint<FpNaive>,
}

// ─── IndexCalcStrategy ───────────────────────────────────────────────────────

/// The index-calculus strategy substrate: factor base + prime-order subgroup + arity.
///
/// Bound to a specific curve. The toy instance uses `semaev_toy()` with `ℓ = 5`,
/// `FB_SIZE = 6`, `m = 2`. The decomposition, collection, linear algebra, and recovery
/// steps all consume this as the fixed strategy.
#[derive(Clone, Debug)]
pub struct IndexCalcStrategy {
    /// The curve (the frozen `semaev_toy()` fixture for the toy instance).
    pub curve: Curve,
    /// The prime-order subgroup modulus ℓ | n (ℓ = 5 for the toy).
    pub ell: Uint<4>,
    /// The factor base: small-x-coordinate points, each with its stable index.
    pub factor_base: Vec<FbPoint>,
    /// The decomposition arity m (m = 2 for the toy; uses S_{m+1} = S_3).
    pub m: usize,
}

impl IndexCalcStrategy {
    /// Build the toy strategy: `semaev_toy()` curve, ℓ = 5, the first `TOY_FB_SIZE`
    /// small-x points, m = 2.
    ///
    /// Enumerates the factor base via QR-test + canonical-root lift over the frozen field.
    /// Returns `IndexCalcError::InvalidSubgroup` if ℓ does not divide n (guard — cannot
    /// happen for the hardcoded toy constants).
    pub fn toy() -> Result<Self, IndexCalcError> {
        let curve = semaev_toy();
        let ell = Uint::<4>::from(TOY_ELL);

        // Guard: ℓ must divide n (the matrix-over-a-field precondition).
        let n_words = curve.n.as_words();
        let n_u64 = n_words[0]; // toy n = 60 fits in u64
        if n_u64 % TOY_ELL != 0 {
            return Err(IndexCalcError::InvalidSubgroup { ell: TOY_ELL, n: n_u64 });
        }

        let factor_base = Self::enumerate_factor_base(&curve, TOY_FB_SIZE)?;

        Ok(IndexCalcStrategy { curve, ell, factor_base, m: TOY_M })
    }

    /// Enumerate factor-base points: ascending x ∈ [0, p), keep QR, lift canonical y,
    /// take the first `fb_size` such points.
    ///
    /// The canonical y is the smaller of the two square roots (by `to_uint()` comparison).
    /// Returns `IndexCalcError::FactorBaseTooSmall` if fewer than `fb_size` QR points
    /// exist in [0, p).
    pub fn enumerate_factor_base(
        curve: &Curve,
        fb_size: usize,
    ) -> Result<Vec<FbPoint>, IndexCalcError> {
        let p = &curve.p;
        let a = FpNaive::from_uint(curve.a, p);
        let b = FpNaive::from_uint(curve.b, p);

        let mut result = Vec::with_capacity(fb_size);
        let mut index = 0usize;

        // Toy-scale: p = 47 fits in u64; iterate x = 0, 1, 2, … until p.
        // SCALE: toy-scale only — crypto-scale p would need a Uint<4> iterator.
        let p_u64 = p.as_words()[0];
        debug_assert!(
            p.as_words()[1] == 0 && p.as_words()[2] == 0 && p.as_words()[3] == 0,
            "enumerate_factor_base: p >= 2^64 not supported (toy parameters only)"
        );

        for x_u64 in 0..p_u64 {
            let x = FpNaive::from_u64(x_u64, p);

            // v = x³ + ax + b mod p
            let v = x.square(p).mul(&x, p)   // x³
                .add(&a.mul(&x, p), p)        // + ax
                .add(&b, p);                  // + b

            // Check QR: Legendre symbol = 1.
            if v.legendre(p) != 1 {
                continue;
            }

            // Lift y: Tonelli–Shanks gives one root; take the canonical (smaller) one.
            let y = v.sqrt(p).expect("legendre = 1 implies sqrt exists");
            let y_neg = y.neg(p);

            // Canonical y: the smaller of y and -y by their Uint<4> representation.
            let canonical_y = if y.to_uint() <= y_neg.to_uint() { y } else { y_neg };

            result.push(FbPoint {
                index,
                point: AffinePoint::Finite { x, y: canonical_y },
            });
            index += 1;

            if result.len() == fb_size {
                return Ok(result);
            }
        }

        // Exhausted [0, p) without finding fb_size QR points.
        Err(IndexCalcError::FactorBaseTooSmall { requested: fb_size, found: result.len() })
    }

    /// The factor-base size (the relation exponent-vector dimension).
    pub fn fb_size(&self) -> usize {
        self.factor_base.len()
    }

    /// The ℓ-order subgroup generator `G_ℓ = (n/ℓ)·G`.
    ///
    /// For the toy fixture: `G_ℓ = 12·G` (since n/ℓ = 60/5 = 12). The subgroup
    /// `{0, G_ℓ, 2·G_ℓ, 3·G_ℓ, 4·G_ℓ}` has order ℓ = 5; `ℓ·G_ℓ = ∞`.
    pub fn subgroup_generator(&self) -> AffinePoint<FpNaive> {
        // n/ℓ: toy-scale division (both fit in u64).
        let n_u64 = self.curve.n.as_words()[0];
        let ell_u64 = self.ell.as_words()[0];
        let cofactor = n_u64 / ell_u64;
        let cofactor_uint = Uint::<4>::from(cofactor);

        let g: AffinePoint<FpNaive> = self.curve.generator();
        self.curve.scalar_mul(&g, &cofactor_uint)
    }
}

// ─── Relation ────────────────────────────────────────────────────────────────

/// One index-calculus relation: a decomposition of `R = a·G + b·Q` over the factor base.
///
/// The exponent vector is **sparse over `F_ℓ`** — `Vec<(fb_index, exp mod ℓ)>` — exactly
/// the row shape `build_ek_matrix` pushes into `FlSparseRow`, so the adapter is a
/// near-identity copy (no re-encoding). Provenance `(a, b)` is load-bearing for DLP
/// recovery; never dropped.
///
/// **Invariants:**
/// - `exponents` is sorted by factor-base index with no duplicate index.
/// - All exponents are non-zero (zero entries are dropped by `from_decomposition`).
/// - Exponents live in `F_ℓ` (the prime-order subgroup — the matrix-over-a-field
///   precondition).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Relation {
    /// Provenance: the multiple `R = a·G + b·Q` this relation decomposes.
    /// The DLP recovery step reads it to reconstruct `log_G(Q) mod ℓ`.
    pub a: u64,
    /// Provenance: the multiple `R = a·G + b·Q` this relation decomposes.
    pub b: u64,
    /// Sparse exponent vector over `F_ℓ`: `(factor-base index, exponent mod ℓ)` pairs,
    /// sorted by index with no duplicate index and no zero exponent.
    pub exponents: Vec<(usize, FpNaive)>,
}

impl Relation {
    /// Construct from a list of factor-base indices the decomposition `R = Σ P_{i_k}` hit.
    ///
    /// Accumulates repeated indices into `F_ℓ` exponents (count occurrences mod ℓ), sorts
    /// by index, drops zero entries. The `ell` parameter is the prime-order subgroup modulus
    /// (used as the field modulus for the exponents).
    pub fn from_decomposition(a: u64, b: u64, fb_indices: &[usize], ell: &Uint<4>) -> Self {
        // Accumulate: count occurrences of each factor-base index.
        // Use a BTreeMap to keep sorted order and merge duplicates.
        use std::collections::BTreeMap;
        let mut counts: BTreeMap<usize, u64> = BTreeMap::new();
        for &i in fb_indices {
            *counts.entry(i).or_insert(0) += 1;
        }

        // Convert counts to F_ℓ exponents, dropping zeros.
        let exponents: Vec<(usize, FpNaive)> = counts
            .into_iter()
            .filter_map(|(i, count)| {
                let exp = FpNaive::from_u64(count, ell);
                if exp.is_zero(ell) { None } else { Some((i, exp)) }
            })
            .collect();

        Relation { a, b, exponents }
    }

    /// The exponent of factor-base point `i` in this relation (zero `F_ℓ` element if absent).
    ///
    /// Performs a linear scan over the sparse `exponents` vector (small factor bases make
    /// this acceptable; the linear-algebra adapter uses the sorted structure directly).
    pub fn exponent(&self, i: usize, ell: &Uint<4>) -> FpNaive {
        self.exponents
            .iter()
            .find(|(idx, _)| *idx == i)
            .map(|(_, e)| e.clone())
            .unwrap_or_else(|| FpNaive::zero(ell))
    }
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;

    const ELL: u64 = TOY_ELL; // 5

    fn ell() -> Uint<4> {
        Uint::<4>::from(ELL)
    }

    fn fp_ell(v: u64) -> FpNaive {
        FpNaive::from_u64(v, &ell())
    }

    // ── IndexCalcStrategy ────────────────────────────────────────────────────

    #[test]
    fn toy_strategy_builds_and_has_correct_size() {
        let s = IndexCalcStrategy::toy().unwrap();
        assert_eq!(s.fb_size(), TOY_FB_SIZE);
        assert_eq!(s.m, TOY_M);
        assert_eq!(s.ell, Uint::<4>::from(TOY_ELL));
    }

    #[test]
    fn factor_base_points_are_on_curve() {
        let s = IndexCalcStrategy::toy().unwrap();
        for fb in &s.factor_base {
            assert!(
                s.curve.is_on_curve(&fb.point),
                "factor-base point {} not on curve: {:?}",
                fb.index,
                fb.point
            );
        }
    }

    #[test]
    fn factor_base_indices_are_sequential() {
        let s = IndexCalcStrategy::toy().unwrap();
        for (i, fb) in s.factor_base.iter().enumerate() {
            assert_eq!(fb.index, i, "factor-base index mismatch at position {i}");
        }
    }

    #[test]
    fn ell_divides_n() {
        let s = IndexCalcStrategy::toy().unwrap();
        let n_u64 = s.curve.n.as_words()[0];
        let ell_u64 = s.ell.as_words()[0];
        assert_eq!(n_u64 % ell_u64, 0, "ℓ = {ell_u64} does not divide n = {n_u64}");
    }

    #[test]
    fn subgroup_generator_has_order_ell() {
        let s = IndexCalcStrategy::toy().unwrap();
        let g_ell = s.subgroup_generator();
        // ℓ·G_ℓ should be the point at infinity.
        let result = s.curve.scalar_mul(&g_ell, &s.ell);
        assert!(
            result.is_infinity(),
            "ℓ·G_ℓ should be ∞ (order ℓ), got {:?}",
            result
        );
    }

    #[test]
    fn subgroup_generator_is_not_infinity() {
        let s = IndexCalcStrategy::toy().unwrap();
        let g_ell = s.subgroup_generator();
        assert!(!g_ell.is_infinity(), "G_ℓ = (n/ℓ)·G should not be ∞");
    }

    // ── Relation ─────────────────────────────────────────────────────────────

    #[test]
    fn relation_from_decomposition_no_repeats() {
        // Indices [0, 1]: each appears once → exponents [(0, 1), (1, 1)] in F_5.
        let r = Relation::from_decomposition(3, 7, &[0, 1], &ell());
        assert_eq!(r.a, 3);
        assert_eq!(r.b, 7);
        assert_eq!(r.exponents.len(), 2);
        assert_eq!(r.exponents[0], (0, fp_ell(1)));
        assert_eq!(r.exponents[1], (1, fp_ell(1)));
    }

    #[test]
    fn relation_from_decomposition_with_repeats() {
        // Indices [0, 0, 1]: index 0 appears twice → exponent 2 mod 5 = 2.
        let r = Relation::from_decomposition(1, 2, &[0, 0, 1], &ell());
        assert_eq!(r.exponents.len(), 2);
        assert_eq!(r.exponents[0], (0, fp_ell(2)));
        assert_eq!(r.exponents[1], (1, fp_ell(1)));
    }

    #[test]
    fn relation_from_decomposition_wrap_mod_ell() {
        // Index 0 appears 5 times → exponent 5 mod 5 = 0 → dropped.
        let r = Relation::from_decomposition(0, 0, &[0, 0, 0, 0, 0], &ell());
        assert!(r.exponents.is_empty(), "5 mod 5 = 0 should be dropped");
    }

    #[test]
    fn relation_exponent_lookup() {
        let r = Relation::from_decomposition(0, 0, &[2, 4], &ell());
        assert_eq!(r.exponent(2, &ell()), fp_ell(1));
        assert_eq!(r.exponent(4, &ell()), fp_ell(1));
        // Absent index returns zero.
        assert_eq!(r.exponent(0, &ell()), fp_ell(0));
        assert_eq!(r.exponent(99, &ell()), fp_ell(0));
    }

    #[test]
    fn relation_exponents_sorted_by_index() {
        // Provide indices out of order; result must be sorted.
        let r = Relation::from_decomposition(0, 0, &[3, 1, 2], &ell());
        let indices: Vec<usize> = r.exponents.iter().map(|(i, _)| *i).collect();
        let mut sorted = indices.clone();
        sorted.sort_unstable();
        assert_eq!(indices, sorted, "exponents must be sorted by index");
    }
}
