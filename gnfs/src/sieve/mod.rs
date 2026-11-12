//! Sieve substrate for GNFS: factor bases, norms, and the relation type.
//!
//! This module is the entry point for the ``gnfs::sieve`` sub-crate. It provides:
//!
//! - [`factor_base`] — two-sided factor base construction (rational + algebraic).
//! - [`norms`] — rational norm and algebraic value computation, and the bridge to ``Uint<4>``.
//! - [`ExponentVector`] — sparse exponent vector over a factor base (rational or algebraic side).
//! - [`Relation`] — a coprime pair ``(a, b)`` with the rational norm and the algebraic value both
//!   smooth over the factor bases.
//! - [`RelationError`] — error type for relation verification.
//! - [`checked`] — the checked relation context: [`RelationChecker`], [`CheckedRelation`],
//!   [`CheckedRelations`], and [`classify_prime`], which re-derive every norm, sign, and affine
//!   ideal label of a [`Relation`] before a matrix is built from it.
//!
//! # Background
//!
//! A **relation** in GNFS is a coprime pair ``(a, b)`` for which both the rational norm
//! ``N_rat(a, b) = a − b·m`` and the algebraic value ``N_alg(a, b) = F(a, b) = b^d · f(a/b)``
//! are smooth over their respective factor bases. For ``f`` with leading coefficient ``c`` and
//! root ``α``, ``F(a, b) = c · N_{K/ℚ}(a − b·α)``; "norm" alone names ``N``. The value is
//! computed by [`norms::algebraic_norm`] (tested by ``algebraic_norm_cubic_hand_computed`` in
//! ``norms.rs`` and, for a non-monic ``f``, ``kat_good_and_exceptional_relations_accepted`` in
//! ``tests/checked_relations_kat.rs``). Each relation yields one row of the relation matrix that
//! the filtering step and linear algebra step consume.
//!
//! # Shared by factoring and discrete logarithms
//!
//! The ``Relation`` type stores full integer exponents (``u32``), not pre-reduced GF(2) parities.
//! The linear algebra step reduces to parities via ``rational_row_gf2`` / ``algebraic_row_gf2``;
//! the NFS-DL relation collection reads the integer exponents directly for GF(ℓ) linear algebra.
//! This over-specification for NFS-DL is load-bearing for both the factoring and the
//! discrete-logarithm code: re-narrowing the relation type after either consumer uses it would be
//! a breaking change to a stable interface.

pub mod checked;
pub mod factor_base;
pub mod lattice;
pub mod line;
pub mod norms;
pub mod special_q;

pub use checked::{
    classify_prime, CheckedRelation, CheckedRelations, CheckedRelationsError, ExponentSide,
    FactorBaseRejection, PrimeClass, RelationChecker, RelationRejection,
};
pub use factor_base::{AlgebraicPrime, FactorBase};
pub use lattice::{lattice_sieve, LatticeBasis, LatticeSieveConfig, LatticeSieveResult};
pub use line::{line_sieve, LineSieveConfig};
pub use norms::{algebraic_norm, norm_sign, norm_to_uint, rational_norm, NormBridgeError};
pub use special_q::{special_q_sieve, SpecialQConfig, SpecialQResult};

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use shared_numth::SmoothWitness;

use crate::polyselect::PolyPair;

// ─── ExponentVector ───────────────────────────────────────────────────────────

/// Exponent vector for one side of a relation (rational or algebraic).
///
/// Stores full integer exponents, not GF(2) parities. The linear algebra step reduces to
/// parities for the nullspace computation; NFS-DL relation collection reads the integer
/// exponents directly for GF(ℓ) linear algebra.
///
/// # Representation
///
/// Sparse: only non-zero exponents are stored. The ``index`` field refers to the factor-base
/// index (column number in the matrix). Entries are sorted by index.
///
/// # Over-specification for NFS-DL
///
/// The exponent type is ``u32``, not ``u8`` or ``bool``. This accommodates:
///
/// - NFS-factoring (linear algebra step): exponents are small (typically 1–3), reduced mod 2.
/// - NFS-DL (relation collection): exponents are reduced mod ℓ where ℓ is the target group
///   order; ℓ can be large, but exponents before reduction are still small integers.
///
/// The ``u32`` type is the smallest that accommodates both without overflow risk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExponentVector {
    /// Sparse representation: (factor-base index, exponent) pairs.
    ///
    /// Sorted by index. Exponents are always > 0 (zeros are not stored).
    pub entries: Vec<(usize, u32)>,
}

impl ExponentVector {
    /// Construct an empty exponent vector.
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    /// Construct from a ``SmoothWitness`` and a factor-base index lookup.
    ///
    /// The ``index_fn`` maps each prime ``p`` to its factor-base index. Returns ``None`` if
    /// any prime in the witness is not in the factor base (i.e., the witness has a prime
    /// outside the factor base, which should not happen for a fully smooth norm).
    ///
    /// Note: this function does not check ``witness.cofactor``; the caller (``Relation::new``)
    /// is responsible for rejecting witnesses with ``cofactor > 1``.
    ///
    /// :param witness: The smoothness witness from ``trial_smooth``.
    /// :param index_fn: Maps a prime ``p`` to its factor-base index, or ``None`` if not present.
    /// :returns: ``Some(ExponentVector)`` if all primes map to indices, ``None`` otherwise.
    pub fn from_witness<F>(witness: &SmoothWitness, index_fn: F) -> Option<Self>
    where
        F: Fn(u64) -> Option<usize>,
    {
        let mut entries: Vec<(usize, u32)> = Vec::with_capacity(witness.factors.len());
        for &(p, e) in &witness.factors {
            let idx = index_fn(p)?;
            entries.push((idx, e));
        }
        // Sort by index for canonical representation.
        entries.sort_by_key(|&(idx, _)| idx);
        Some(Self { entries })
    }

    /// Get the exponent at a given factor-base index (0 if not present).
    ///
    /// :param index: The factor-base column index.
    /// :returns: The exponent at ``index``, or 0 if not present.
    pub fn get(&self, index: usize) -> u32 {
        // Binary search since entries are sorted by index.
        match self.entries.binary_search_by_key(&index, |&(idx, _)| idx) {
            Ok(pos) => self.entries[pos].1,
            Err(_) => 0,
        }
    }

    /// Iterate over (index, exponent) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (usize, u32)> + '_ {
        self.entries.iter().copied()
    }

    /// Number of non-zero entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True if all exponents are zero (empty vector).
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl Default for ExponentVector {
    fn default() -> Self {
        Self::new()
    }
}

// ─── RelationError ────────────────────────────────────────────────────────────

/// Error type for relation verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RelationError {
    /// gcd(a, b) ≠ 1.
    ///
    /// :param gcd: The actual gcd(a, b).
    NotCoprime { gcd: BigInt },

    /// Rational exponents do not reconstruct the norm.
    ///
    /// :param expected: The actual |N_rat(a, b)|.
    /// :param got: The product reconstructed from the exponent vector.
    RationalMismatch { expected: BigInt, got: BigInt },

    /// Algebraic exponents do not reconstruct the norm.
    ///
    /// :param expected: The actual |N_alg(a, b)|.
    /// :param got: The product reconstructed from the exponent vector.
    AlgebraicMismatch { expected: BigInt, got: BigInt },

    /// Sign does not match the actual norm sign.
    ///
    /// :param expected: The sign of the actual rational norm.
    /// :param actual: The sign stored in the relation.
    SignMismatch { expected: bool, actual: bool },

    /// A prime in the witness is not in the factor base.
    ///
    /// :param prime: The prime that was not found.
    /// :param side: "rational" or "algebraic".
    PrimeNotInBase { prime: u64, side: &'static str },
}

impl std::fmt::Display for RelationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotCoprime { gcd } => write!(f, "gcd(a, b) = {gcd} ≠ 1"),
            Self::RationalMismatch { expected, got } => {
                write!(f, "rational norm mismatch: expected |N_rat| = {expected}, got {got}")
            }
            Self::AlgebraicMismatch { expected, got } => {
                write!(f, "algebraic norm mismatch: expected |N_alg| = {expected}, got {got}")
            }
            Self::SignMismatch { expected, actual } => {
                write!(
                    f,
                    "sign mismatch: N_rat sign is {expected} (negative={expected}), \
                     stored rational_sign = {actual}"
                )
            }
            Self::PrimeNotInBase { prime, side } => {
                write!(f, "prime {prime} not in {side} factor base")
            }
        }
    }
}

impl std::error::Error for RelationError {}

// ─── Relation ─────────────────────────────────────────────────────────────────

/// A relation: a coprime pair (a, b) whose rational norm and algebraic value F(a, b) are both
/// smooth over the factor bases.
///
/// # Invariants (checked by ``verify``)
///
/// 1. gcd(a, b) = 1 (coprimality).
/// 2. The rational exponent vector reconstructs |N_rat(a, b)| over the rational FB.
/// 3. The algebraic exponent vector reconstructs |N_alg(a, b)| = |F(a, b)| over the algebraic FB.
/// 4. Both witnesses are fully smooth (cofactor = 1).
///
/// # Sign handling
///
/// The ``rational_sign`` field is true iff ``N_rat(a, b) = a − b·m < 0``. This is the
/// "−1 column" for the linear algebra step: the product of all selected relations' rational
/// norms must be a perfect square, which requires the sign product to be +1.
///
/// The algebraic value's sign is not stored separately and is not a matrix column. A field
/// need not have a real embedding (totally imaginary fields have none), so the algebraic
/// square root is reconstructed and verified exactly, with no sign convention taken from ℝ.
/// (The quadratic-character columns in the linear algebra step are necessary filters on
/// the algebraic side: they exclude candidates, but they do not guarantee that a square
/// root exists in K, so every candidate is still verified exactly.)
///
/// # Algebraic parity
///
/// The algebraic exponents record the valuation of the homogeneous value ``F(a, b) = c ·
/// N(a − b·α)`` at each degree-1 ideal ``(p, r)`` with ``a ≡ b·r (mod p)``. Their parity is a
/// square constraint only at primes ``p`` that divide neither the leading coefficient nor the
/// discriminant of the monic form; [`checked::RelationChecker`] classifies the other primes as
/// exceptional. The primitive affine collector never returns a pair with a projective factor
/// (``p`` dividing both the leading coefficient and ``b``).
///
/// # Reuse by NFS-DL
///
/// NFS-DL uses the same relation structure but interprets exponents mod ℓ instead of mod 2.
/// The integer exponents stored here support both interpretations without changing the relation
/// format. NFS-DL
/// relation collection may add Schirokauer-map columns; the ``ExponentVector`` structure
/// accommodates additional entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Relation {
    /// The a-coordinate of the relation (can be negative).
    pub a: BigInt,
    /// The b-coordinate of the relation (always positive: b ≥ 1).
    pub b: BigInt,
    /// Exponent vector over the rational factor base.
    pub rational_exponents: ExponentVector,
    /// Exponent vector over the algebraic factor base.
    pub algebraic_exponents: ExponentVector,
    /// True if the rational norm a − b·m is negative.
    pub rational_sign: bool,
}

impl Relation {
    /// Construct a relation from ``(a, b)`` and the smoothness witnesses.
    ///
    /// Returns ``None`` if either witness has ``cofactor > 1`` (not fully smooth) or if
    /// ``gcd(a, b) ≠ 1``.
    ///
    /// For the algebraic side, the mapping from primes to ideals uses the sieve condition:
    /// the ideal ``(p, α − r)`` divides ``(a − b·α)`` iff ``a ≡ r·b (mod p)``. For each
    /// prime ``p`` in the algebraic witness, we find the root ``r`` satisfying this condition
    /// and look up ``fb.algebraic_index(p, r)``.
    ///
    /// :param a: The a-coordinate.
    /// :param b: The b-coordinate (positive).
    /// :param rational_witness: Smoothness witness for |N_rat(a, b)|.
    /// :param algebraic_witness: Smoothness witness for |N_alg(a, b)|.
    /// :param rational_sign: True if N_rat(a, b) < 0.
    /// :param fb: The factor base (for index lookup).
    /// :returns: ``Some(Relation)`` if fully smooth and coprime, ``None`` otherwise.
    pub fn new(
        a: BigInt,
        b: BigInt,
        rational_witness: &SmoothWitness,
        algebraic_witness: &SmoothWitness,
        rational_sign: bool,
        fb: &FactorBase,
    ) -> Option<Self> {
        // Reject if either witness is not fully smooth.
        if !rational_witness.is_smooth() || !algebraic_witness.is_smooth() {
            return None;
        }

        // Reject if gcd(a, b) ≠ 1.
        let g = gcd_bigint(a.abs(), b.abs());
        if !g.is_one() {
            return None;
        }

        // Build rational exponent vector: map each prime p to fb.rational_index(p).
        let rational_exponents = ExponentVector::from_witness(rational_witness, |p| {
            fb.rational_index(p)
        })?;

        // Build algebraic exponent vector: for each prime p in the witness, find the root r
        // such that a ≡ r·b (mod p), then look up fb.algebraic_index(p, r).
        let algebraic_exponents = ExponentVector::from_witness(algebraic_witness, |p| {
            // Find r ∈ [0, p) with a ≡ r·b (mod p).
            // Equivalently: r ≡ a · b^{-1} (mod p), but we can just check each r in 0..p.
            // For small p (≤ B_alg, typically ≤ 10^6), this is fast.
            let p_big = BigInt::from(p);
            let a_mod = mod_reduce(&a, &p_big);
            let b_mod = mod_reduce(&b, &p_big);
            // Find r such that r * b_mod ≡ a_mod (mod p).
            let r = find_root_for_ideal(&a_mod, &b_mod, p)?;
            fb.algebraic_index(p, r)
        })?;

        Some(Relation { a, b, rational_exponents, algebraic_exponents, rational_sign })
    }

    /// Verify the relation's invariants against the polynomial pair and factor base.
    ///
    /// Checks:
    ///
    /// 1. gcd(a, b) = 1.
    /// 2. Rational exponents reconstruct |a − b·m| (via product of p^e).
    /// 3. Algebraic exponents reconstruct |N_alg(a, b)| (via product of p^e).
    /// 4. The sign matches: rational_sign == (a − b·m < 0).
    ///
    /// :param poly: The polynomial pair (provides f, m for norm computation).
    /// :param fb: The factor base (provides prime lookup for reconstruction).
    /// :returns: ``Ok(())`` if all checks pass, ``Err(RelationError)`` otherwise.
    pub fn verify(&self, poly: &PolyPair, fb: &FactorBase) -> Result<(), RelationError> {
        // Check 1: gcd(a, b) = 1.
        let g = gcd_bigint(self.a.abs(), self.b.abs());
        if !g.is_one() {
            return Err(RelationError::NotCoprime { gcd: g });
        }

        // Check 4: sign matches.
        let actual_rat_norm = norms::rational_norm(&self.a, &self.b, &poly.m);
        let actual_sign = norms::norm_sign(&actual_rat_norm);
        if actual_sign != self.rational_sign {
            return Err(RelationError::SignMismatch {
                expected: actual_sign,
                actual: self.rational_sign,
            });
        }

        // Check 2: rational exponents reconstruct |N_rat|.
        let expected_rat = actual_rat_norm.abs();
        let got_rat = reconstruct_norm_from_exponents(&self.rational_exponents, &fb.rational_primes);
        if expected_rat != got_rat {
            return Err(RelationError::RationalMismatch { expected: expected_rat, got: got_rat });
        }

        // Check 3: algebraic exponents reconstruct |N_alg|.
        let actual_alg_norm = norms::algebraic_norm(&self.a, &self.b, &poly.f);
        let expected_alg = actual_alg_norm.abs();
        let got_alg = reconstruct_norm_from_alg_exponents(&self.algebraic_exponents, &fb.algebraic_ideals);
        if expected_alg != got_alg {
            return Err(RelationError::AlgebraicMismatch { expected: expected_alg, got: got_alg });
        }

        Ok(())
    }

    /// Convert the rational exponents to a GF(2) row for the linear algebra step.
    ///
    /// Returns a bit vector where bit ``i`` is ``exponents[i] mod 2``.
    /// The sign column is prepended: bit 0 is 1 iff ``rational_sign`` is true.
    ///
    /// The length of the returned vector is ``fb.rational_size() + 1`` (sign column prepended).
    ///
    /// :param fb: The factor base (for size information).
    /// :returns: GF(2) row for the rational side, with sign column prepended.
    pub fn rational_row_gf2(&self, fb: &FactorBase) -> Vec<bool> {
        let mut row = vec![false; 1 + fb.rational_size()];
        // Bit 0: sign column.
        row[0] = self.rational_sign;
        // Bits 1..=rational_size: exponent parities.
        for (idx, exp) in self.rational_exponents.iter() {
            row[1 + idx] = (exp % 2) == 1;
        }
        row
    }

    /// Convert the algebraic exponents to a GF(2) row for the linear algebra step.
    ///
    /// Returns a bit vector where bit ``i`` is ``exponents[i] mod 2``.
    /// Obstruction columns (quadratic characters) are appended as zeros; the linear algebra
    /// step fills them in.
    ///
    /// The length of the returned vector is ``fb.algebraic_size() + fb.obstruction_count``.
    ///
    /// :param fb: The factor base (for size and obstruction count).
    /// :returns: GF(2) row for the algebraic side, with obstruction columns appended as zeros.
    pub fn algebraic_row_gf2(&self, fb: &FactorBase) -> Vec<bool> {
        let mut row = vec![false; fb.algebraic_size() + fb.obstruction_count];
        for (idx, exp) in self.algebraic_exponents.iter() {
            if idx < fb.algebraic_size() {
                row[idx] = (exp % 2) == 1;
            }
        }
        // Obstruction columns remain false (zeros); the linear algebra step fills them in.
        row
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Reduce ``a`` into the canonical range [0, m) for m > 0.
fn mod_reduce(a: &BigInt, m: &BigInt) -> BigInt {
    let r = a % m;
    if r < BigInt::zero() { r + m } else { r }
}

/// Compute gcd(a, b) using the Euclidean algorithm.
///
/// Both ``a`` and ``b`` should be non-negative (absolute values).
fn gcd_bigint(mut a: BigInt, mut b: BigInt) -> BigInt {
    while !b.is_zero() {
        let t = b.clone();
        b = a % &t;
        a = t;
    }
    a
}

/// Find the root r ∈ [0, p) such that r * b_mod ≡ a_mod (mod p).
///
/// Returns ``None`` if no such root exists in the factor base (which would indicate the prime
/// is not actually a factor of the algebraic value ``F(a, b)``, contradicting the smoothness
/// witness).
///
/// For small p (≤ B_alg, typically ≤ 10^6), linear search is fast.
fn find_root_for_ideal(a_mod: &BigInt, b_mod: &BigInt, p: u64) -> Option<u64> {
    let p_big = BigInt::from(p);
    for r in 0..p {
        let r_big = BigInt::from(r);
        let lhs = mod_reduce(&(r_big * b_mod), &p_big);
        if lhs == *a_mod {
            return Some(r);
        }
    }
    None
}

/// Reconstruct the norm magnitude from a rational exponent vector.
///
/// Computes ``Π p^e`` for each ``(index, e)`` in the exponent vector, where
/// ``p = rational_primes[index]``.
fn reconstruct_norm_from_exponents(ev: &ExponentVector, primes: &[u64]) -> BigInt {
    let mut product = BigInt::one();
    for (idx, exp) in ev.iter() {
        if idx < primes.len() {
            let p = BigInt::from(primes[idx]);
            let mut pe = BigInt::one();
            for _ in 0..exp {
                pe *= &p;
            }
            product *= pe;
        }
    }
    product
}

/// Reconstruct the norm magnitude from an algebraic exponent vector.
///
/// Computes ``Π p^e`` for each ``(index, e)`` in the exponent vector, where
/// ``p = algebraic_ideals[index].p``.
fn reconstruct_norm_from_alg_exponents(
    ev: &ExponentVector,
    ideals: &[AlgebraicPrime],
) -> BigInt {
    let mut product = BigInt::one();
    for (idx, exp) in ev.iter() {
        if idx < ideals.len() {
            let p = BigInt::from(ideals[idx].p);
            let mut pe = BigInt::one();
            for _ in 0..exp {
                pe *= &p;
            }
            product *= pe;
        }
    }
    product
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    // ── ExponentVector ───────────────────────────────────────────────────────

    #[test]
    fn exponent_vector_new_is_empty() {
        let ev = ExponentVector::new();
        assert!(ev.is_empty());
        assert_eq!(ev.len(), 0);
    }

    #[test]
    fn exponent_vector_get_missing_returns_zero() {
        let ev = ExponentVector::new();
        assert_eq!(ev.get(0), 0);
        assert_eq!(ev.get(42), 0);
    }

    #[test]
    fn exponent_vector_from_witness_maps_primes() {
        // Witness: [(2, 3), (5, 1)] — 2^3 * 5^1 = 40.
        let witness = SmoothWitness {
            factors: vec![(2u64, 3u32), (5u64, 1u32)],
            cofactor: crypto_bigint::Uint::<4>::ONE,
        };
        // index_fn: 2 → 0, 5 → 2, others → None.
        let ev = ExponentVector::from_witness(&witness, |p| match p {
            2 => Some(0),
            5 => Some(2),
            _ => None,
        });
        let ev = ev.expect("should succeed");
        assert_eq!(ev.get(0), 3); // 2^3
        assert_eq!(ev.get(1), 0); // not present
        assert_eq!(ev.get(2), 1); // 5^1
    }

    #[test]
    fn exponent_vector_from_witness_returns_none_for_missing_prime() {
        let witness = SmoothWitness {
            factors: vec![(11u64, 1u32)], // 11 not in factor base
            cofactor: crypto_bigint::Uint::<4>::ONE,
        };
        let ev = ExponentVector::from_witness(&witness, |p| match p {
            2 => Some(0),
            3 => Some(1),
            _ => None,
        });
        assert!(ev.is_none(), "prime 11 not in factor base should return None");
    }

    // ── gcd_bigint ───────────────────────────────────────────────────────────

    #[test]
    fn gcd_coprime() {
        assert_eq!(gcd_bigint(bi(7), bi(13)), bi(1));
        assert_eq!(gcd_bigint(bi(1), bi(100)), bi(1));
    }

    #[test]
    fn gcd_non_coprime() {
        assert_eq!(gcd_bigint(bi(12), bi(8)), bi(4));
        assert_eq!(gcd_bigint(bi(100), bi(25)), bi(25));
    }
}
