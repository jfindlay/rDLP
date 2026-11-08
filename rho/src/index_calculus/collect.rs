//! Relation collection for the index-calculus ECDLP solver.
//!
//! This module implements the relation-collection loop: generate the multiples
//! `R = a·G + b·Q`, attempt to decompose each over the factor base via
//! [`crate::index_calculus::decompose::decompose`], and record successful decompositions
//! as `Relation` rows. [`collect_relations`] stops once the system is over-determined
//! (≥ `strategy.fb_size() + 1` relations); [`RelationStream`] yields relations one at a
//! time so a caller can ask for more.
//!
//! # Algorithm
//!
//! For the pairs `(a, b)` in a fixed systematic order (anti-diagonals `a + b = 1, 2, …`,
//! see [`RelationStream`]):
//! 1. Compute `R = a·G + b·Q` using `Curve::scalar_mul` and `Curve::add_jacobian`.
//! 2. Skip if `R = ∞` (decomposition is undefined for the point at infinity).
//! 3. Call `decompose(R, strategy)`. On `Some(fb_points)`, extract the factor-base
//!    indices and construct a `Relation` via `Relation::from_decomposition`.
//!
//! The pairs `(a, b)` range over `a ∈ [1, n)`, `b ∈ [0, n)`, which is every multiple of
//! `G` and `Q` modulo the group order `n`: that is the *relation-enumeration bound*. It
//! bounds how many relations exist to be collected; it is not a search over candidate
//! values of the logarithm.
//!
//! # Scale
//!
//! The toy fixture has `n = 60`, `|FB| = 6`, `m = 2`, so the full enumeration is at most
//! `n²` pairs. A larger instance would sample `(a, b)` at random and use a stronger
//! decomposition oracle (Gröbner-basis solving of the Semaev system).

use crypto_bigint::Uint;

use shared_field::FpNaive4 as FpNaive;

use crate::curve::{AffinePoint, JacobianPoint};
use crate::index_calculus::decompose::decompose;
use crate::index_calculus::strategy::{IndexCalcStrategy, Relation};
use crate::index_calculus::IndexCalcError;

// ─── RelationStream ──────────────────────────────────────────────────────────

/// A deterministic, resumable stream of relations over the `(a, b)` enumeration.
///
/// The stream visits the pairs `(a, b)` with `a ∈ [1, n)` and `b ∈ [0, n)` along
/// anti-diagonals `a + b = 1, 2, …` (within a diagonal, `a` ascending), and yields one
/// [`Relation`] per pair whose multiple `a·G + b·Q` is finite and decomposes over the
/// factor base. Diagonals mix both coefficients early: relations that share one value of
/// `b` (or of `a`) can make `b` a combination of the exponent columns, which hides `x`,
/// whereas a diagonal supplies varied `(a, b)` from the start.
///
/// Any prefix of the stream is the same prefix of any longer collection, so a caller can
/// extend a relation set without recomputing or reordering it.
pub struct RelationStream<'a> {
    g: AffinePoint<FpNaive>,
    q: AffinePoint<FpNaive>,
    strategy: &'a IndexCalcStrategy,
    n: u64,
    /// The current anti-diagonal `a + b`; beyond `2n − 2` the enumeration is exhausted.
    diagonal: u64,
    /// The next `a` on the current diagonal (`b = diagonal − a`).
    a: u64,
}

impl<'a> RelationStream<'a> {
    /// Start a stream at `(a, b) = (1, 0)`.
    pub fn new(
        g: AffinePoint<FpNaive>,
        q: AffinePoint<FpNaive>,
        strategy: &'a IndexCalcStrategy,
    ) -> Self {
        // The group order n (toy scale: fits in u64).
        let n = strategy.curve.n.as_words()[0];
        RelationStream { g, q, strategy, n, diagonal: 1, a: 1 }
    }

    /// Yield the next relation, or `None` once every `(a, b)` pair has been visited.
    pub fn next_relation(&mut self) -> Option<Relation> {
        let curve = &self.strategy.curve;
        let p = &curve.p;

        // Pairs exist for diagonals 1 ..= 2n − 2.
        while self.diagonal <= 2 * self.n - 2 {
            let (a, b) = (self.a, self.diagonal - self.a);
            // Advance to the next pair: a ≤ min(diagonal, n − 1) keeps b ≥ 0 and a < n, and
            // a ≥ diagonal − (n − 1) keeps b < n.
            if self.a < self.diagonal.min(self.n - 1) {
                self.a += 1;
            } else {
                self.diagonal += 1;
                self.a = self.diagonal.saturating_sub(self.n - 1).max(1);
            }

            // R = a·G + b·Q. The sum is formed only when both terms are finite.
            let ag = curve.scalar_mul(&self.g, &Uint::<4>::from(a));
            let bq = curve.scalar_mul(&self.q, &Uint::<4>::from(b));
            let r = if bq.is_infinity() {
                ag
            } else if ag.is_infinity() {
                bq
            } else {
                let jag = JacobianPoint::from_affine(&ag, p);
                let jbq = JacobianPoint::from_affine(&bq, p);
                curve.add_jacobian(&jag, &jbq).to_affine(p)
            };

            if r.is_infinity() {
                continue;
            }

            if let Some(fb_points) = decompose(r, self.strategy) {
                let fb_indices: Vec<usize> = fb_points.iter().map(|fp| fp.index).collect();
                return Some(Relation::from_decomposition(a, b, &fb_indices, &self.strategy.ell));
            }
        }
        None
    }
}

// ─── collect_relations ───────────────────────────────────────────────────────

/// Collect an over-determined system of index-calculus relations.
///
/// Draws relations from a [`RelationStream`] until there are at least
/// `strategy.fb_size() + 1` of them.
///
/// Returns a `Vec<Relation>` with exactly `strategy.fb_size() + 1` entries.
///
/// # Errors
///
/// Returns `IndexCalcError::UnderdeterminedSystem` if the full `(a, b)` enumeration
/// yields fewer relations than that.
pub fn collect_relations(
    g: AffinePoint<FpNaive>,
    q: AffinePoint<FpNaive>,
    strategy: &IndexCalcStrategy,
) -> Result<Vec<Relation>, IndexCalcError> {
    let target = strategy.fb_size() + 1;
    let mut stream = RelationStream::new(g, q, strategy);
    let mut relations: Vec<Relation> = Vec::with_capacity(target);

    while relations.len() < target {
        match stream.next_relation() {
            Some(rel) => relations.push(rel),
            None => {
                return Err(IndexCalcError::UnderdeterminedSystem {
                    found: relations.len(),
                    needed: target,
                });
            }
        }
    }

    Ok(relations)
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index_calculus::strategy::IndexCalcStrategy;
    use crate::semaev::semaev_toy;

    #[test]
    fn collect_returns_enough_relations() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let curve = semaev_toy();
        let g: AffinePoint<FpNaive> = curve.generator();
        // Q = 7·G (arbitrary non-trivial multiple).
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        let relations = collect_relations(g, q, &strategy)
            .expect("collect_relations should succeed for the toy fixture");

        assert!(
            relations.len() >= strategy.fb_size() + 1,
            "expected at least {} relations, got {}",
            strategy.fb_size() + 1,
            relations.len()
        );
    }

    #[test]
    fn stream_extends_collection_prefix() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let curve = semaev_toy();
        let g: AffinePoint<FpNaive> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        let collected = collect_relations(g.clone(), q.clone(), &strategy).unwrap();
        let mut stream = RelationStream::new(g, q, &strategy);
        for (i, rel) in collected.iter().enumerate() {
            assert_eq!(stream.next_relation().as_ref(), Some(rel), "prefix differs at {i}");
        }
        assert!(stream.next_relation().is_some(), "stream must continue past the prefix");
    }

    #[test]
    fn stream_exhausts_the_enumeration_bound() {
        // Every (a, b) pair is visited exactly once, so the stream ends and stays ended.
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let curve = semaev_toy();
        let g: AffinePoint<FpNaive> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        let mut stream = RelationStream::new(g, q, &strategy);
        let n = strategy.curve.n.as_words()[0];
        let mut count = 0usize;
        let mut last_diagonal = 0u64;
        let mut seen = std::collections::BTreeSet::new();
        while let Some(rel) = stream.next_relation() {
            assert!((1..n).contains(&rel.a) && rel.b < n, "pair ({}, {}) out of range", rel.a, rel.b);
            assert!(seen.insert((rel.a, rel.b)), "pair ({}, {}) repeated", rel.a, rel.b);
            assert!(rel.a + rel.b >= last_diagonal, "diagonals must not decrease");
            last_diagonal = rel.a + rel.b;
            count += 1;
        }
        assert!(count > strategy.fb_size(), "full enumeration yields many relations");
        assert!(stream.next_relation().is_none(), "an exhausted stream stays exhausted");
    }
}
