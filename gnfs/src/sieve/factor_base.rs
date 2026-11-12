//! Two-sided factor base for NFS sieving.
//!
//! The factor base is the fundamental data structure for NFS sieving. It has two sides:
//!
//! - **Rational factor base**: the set of primes ``p ≤ B_rat``. A rational norm
//!   ``N_rat(a, b) = a − b·m`` is smooth over this base if all its prime factors are ≤ ``B_rat``.
//!
//! - **Algebraic factor base**: the set of degree-1 prime ideals ``(p, α − r)`` in ℤ[α], where
//!   ``p ≤ B_alg`` and ``f(r) ≡ 0 (mod p)``. Each prime ``p`` contributes one entry per root of
//!   ``f mod p``. The ideal ``(p, α − r)`` divides ``(a − b·α)`` iff ``a ≡ r·b (mod p)``; the
//!   algebraic value ``F(a, b) = c·N_{K/ℚ}(a − b·α)`` (``N_alg`` in code) is then divisible by
//!   ``p`` (for ``p ∤ c``, the ideal's exponent is part of ``F``'s ``p``-adic valuation).
//!
//! # Column indexing (for filtering and linear algebra)
//!
//! The exponent vectors in ``Relation`` use indices into these factor bases:
//!
//! - Rational side: ``rational_primes[i]`` is the prime at column ``i``.
//! - Algebraic side: ``algebraic_ideals[j].index == j`` gives the column.
//!
//! # Obstruction columns (over-specified for the linear algebra step)
//!
//! The linear algebra step over GF(2) requires additional columns beyond the factor-base
//! primes: the sign column (−1) and quadratic-character columns, which are necessary filters
//! on the algebraic side (they exclude candidates but do not guarantee that an algebraic
//! square root exists). The ``obstruction_count`` field reserves column indices for these;
//! the line sieve does not populate them, but the slot exists so the linear algebra step
//! can extend without changing the factor-base layout.
//!
//! # Bad primes
//!
//! A prime ``p`` is **bad** when ``p | c·disc(f)``, where ``c`` is the leading coefficient of
//! ``f``. This is exactly what [`shared_numfield::is_bad_prime`] computes on the (generally
//! non-monic) ``f`` passed to [`FactorBase::new`]: it tests ``p`` against
//! ``shared_numfield::discriminant(f)``, which for non-monic ``f`` is ``c·disc(f)``. A prime
//! dividing only ``c`` is therefore bad. At toy scale bad primes are prominent (e.g.,
//! disc(x³ − x − 1) = −23, so p = 23 is bad); at cryptographic scale they are marginal.
//!
//! The primes dividing the index [O_K : ℤ[α]] (for integral α) are a subset of the bad primes,
//! not the other way round: a bad prime may have ℤ[α] maximal (x² − 2 at p = 2, where
//! disc = 8 but ℤ[√2] is the full ring of integers). So the flag says that Dedekind's theorem
//! cannot be applied without further checking, not that ℤ[α] fails to be maximal at ``p``.
//!
//! This implementation includes bad primes in the algebraic factor base with the
//! ``is_bad_prime`` flag set, using direct root-finding (``f(r) ≡ 0 mod p``), which is correct
//! for linear factors even at bad primes. The flag is informational. The checked relation path
//! ([`super::checked`]) does not consult it: it classifies a prime by exact divisibility of
//! ``c·disc(F)``, where ``F`` is the monic form of ``f``. Since
//! ``disc(F) = c^((d−1)(d−2))·disc(f)`` for ``d = deg f``, both tests flag the same primes.
//!
//! # Projective limitation
//!
//! Only affine ideals ``(p, α − r)`` are represented. For ``p | c`` and ``p | b``, no affine
//! root ``r`` satisfies ``a ≡ b·r (mod p)``, so [`Relation::new`](super::Relation::new)
//! returns ``None`` and such relations are never collected.

use num_bigint::BigInt;
use num_traits::Zero;
use shared_numfield::{is_bad_prime, IntPoly};
use shared_numth::factor_base_up_to;

// ─── AlgebraicPrime ───────────────────────────────────────────────────────────

/// A degree-1 prime ideal (p, α − r) in the algebraic factor base.
///
/// Represents the prime ideal above ``p`` corresponding to the root ``r`` of ``f mod p``.
/// For good primes (p ∤ c·disc(f), ``c`` the leading coefficient), these are exactly the prime
/// ideals above ``p`` with residue degree 1. For bad primes (p | c·disc(f)), these are the
/// ideals from linear factors of ``f mod p``; higher-degree factors may contribute additional
/// ideals not captured here (see ``is_bad_prime`` flag).
///
/// The ideal ``(p, α − r)`` divides ``(a − b·α)`` iff ``a ≡ r·b (mod p)``, in which case ``p``
/// divides the algebraic value ``F(a, b) = c·N_{K/ℚ}(a − b·α)`` (``N_alg`` in code). This is the
/// sieve condition for the algebraic side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlgebraicPrime {
    /// The rational prime p.
    pub p: u64,
    /// The root r ∈ [0, p) with f(r) ≡ 0 (mod p).
    pub r: u64,
    /// Index of this ideal in the algebraic factor base (for column mapping).
    pub index: usize,
    /// True if p | c·disc(f) (bad prime), ``c`` being the leading coefficient of ``f``.
    ///
    /// A bad prime need not divide the index of ℤ[α]; the flag does not decide maximality.
    /// Bad primes are over-represented at toy scale relative to NFS scale, where their
    /// contribution is marginal. Root-finding via ``f(r) ≡ 0 mod p`` is correct
    /// for linear factors even at bad primes; higher-degree factors of ``f mod p`` may yield
    /// additional prime ideals not representable in two-element form without the full Round 2
    /// algorithm (out of scope).
    pub is_bad_prime: bool,
}

// ─── FactorBase ───────────────────────────────────────────────────────────────

/// Two-sided factor base for NFS sieving.
///
/// The rational factor base is the set of primes ``p ≤ B_rat``.
/// The algebraic factor base is the set of degree-1 prime ideals ``(p, r)`` with
/// ``p ≤ B_alg`` and ``f(r) ≡ 0 (mod p)``.
///
/// # Column indexing (for filtering and linear algebra)
///
/// The exponent vectors in ``Relation`` use indices into these factor bases:
///
/// - Rational side: ``rational_primes[i]`` is the prime at column ``i``.
/// - Algebraic side: ``algebraic_ideals[j].index == j`` gives the column.
///
/// # Obstruction columns (over-specified for the linear algebra step)
///
/// The linear algebra step over GF(2) requires additional columns beyond the factor-base
/// primes: the sign column (−1) and quadratic-character columns, which are necessary filters
/// on the algebraic side and do not guarantee that an algebraic square root exists. The
/// ``obstruction_count`` field reserves column indices for these;
/// the line sieve does not populate them, but the slot exists so the linear algebra step
/// can extend without changing the factor-base layout.
#[derive(Debug, Clone)]
pub struct FactorBase {
    /// Rational factor base: primes p ≤ B_rat, sorted ascending.
    pub rational_primes: Vec<u64>,
    /// Algebraic factor base: degree-1 prime ideals (p, r), sorted by (p, r).
    pub algebraic_ideals: Vec<AlgebraicPrime>,
    /// Rational smoothness bound B_rat.
    pub b_rat: u64,
    /// Algebraic smoothness bound B_alg.
    pub b_alg: u64,
    /// Number of obstruction columns reserved for the linear algebra step (sign + quadratic chars).
    ///
    /// Initialised to 1 (the sign column); the linear algebra step may increase it. The
    /// quadratic-character columns are necessary filters, not a guarantee of a square root.
    pub obstruction_count: usize,
}

impl FactorBase {
    /// Construct the two-sided factor base for polynomial ``f`` with bounds ``(B_rat, B_alg)``.
    ///
    /// Uses ``factor_base_up_to`` for the rational side. For the algebraic side, finds all
    /// roots of ``f mod p`` for each prime ``p ≤ B_alg`` by trial (checking ``f(r) ≡ 0 mod p``
    /// for each ``r ∈ 0..p``). Bad primes (p | c·disc(f), ``c`` the leading coefficient) are
    /// included with the ``is_bad_prime`` flag set.
    ///
    /// The ``obstruction_count`` is initialised to 1 (the sign/−1 column for the linear algebra
    /// step).
    ///
    /// :param f: The algebraic polynomial (from PolyPair).
    /// :param b_rat: Rational smoothness bound.
    /// :param b_alg: Algebraic smoothness bound.
    /// :returns: The two-sided factor base.
    pub fn new(f: &IntPoly, b_rat: u64, b_alg: u64) -> Self {
        // Rational side: all primes ≤ B_rat.
        let rational_primes = factor_base_up_to(b_rat);

        // Algebraic side: for each prime p ≤ B_alg, find all roots of f mod p.
        let alg_primes_list = factor_base_up_to(b_alg);
        let mut algebraic_ideals: Vec<AlgebraicPrime> = Vec::new();
        let mut idx = 0usize;

        for p in alg_primes_list {
            let p_big = BigInt::from(p);
            let bad = is_bad_prime(f, &p_big);

            // Find all roots r ∈ [0, p) with f(r) ≡ 0 (mod p).
            // This is correct for linear factors even at bad primes.
            for r in 0..p {
                let r_big = BigInt::from(r);
                let val = f.eval(&r_big);
                // Reduce mod p into [0, p).
                let rem = mod_reduce(&val, &p_big);
                if rem.is_zero() {
                    algebraic_ideals.push(AlgebraicPrime { p, r, index: idx, is_bad_prime: bad });
                    idx += 1;
                }
            }
        }

        FactorBase {
            rational_primes,
            algebraic_ideals,
            b_rat,
            b_alg,
            // 1 obstruction column: the sign/−1 column for the linear algebra step.
            obstruction_count: 1,
        }
    }

    /// Number of rational factor-base primes.
    pub fn rational_size(&self) -> usize {
        self.rational_primes.len()
    }

    /// Number of algebraic factor-base ideals.
    pub fn algebraic_size(&self) -> usize {
        self.algebraic_ideals.len()
    }

    /// Total matrix width for the linear algebra step: rational + algebraic + obstruction columns.
    pub fn matrix_width(&self) -> usize {
        self.rational_size() + self.algebraic_size() + self.obstruction_count
    }

    /// Look up the index of a rational prime in the factor base.
    ///
    /// Returns ``None`` if ``p`` is not in the rational factor base.
    pub fn rational_index(&self, p: u64) -> Option<usize> {
        self.rational_primes.binary_search(&p).ok()
    }

    /// Look up the index of an algebraic ideal (p, r) in the factor base.
    ///
    /// Returns ``None`` if ``(p, r)`` is not in the algebraic factor base.
    pub fn algebraic_index(&self, p: u64, r: u64) -> Option<usize> {
        // The ideals are sorted by (p, r), so binary search is possible.
        // We search by the (p, r) key and return the stored index.
        let pos = self.algebraic_ideals.binary_search_by_key(&(p, r), |ap| (ap.p, ap.r));
        pos.ok().map(|i| self.algebraic_ideals[i].index)
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Reduce ``a`` into the canonical range [0, m) for m > 0.
fn mod_reduce(a: &BigInt, m: &BigInt) -> BigInt {
    let r = a % m;
    if r < BigInt::zero() { r + m } else { r }
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt {
        BigInt::from(n)
    }

    /// f(x) = x³ − x − 1 (the classic NFS toy polynomial).
    fn f_cubic() -> IntPoly {
        IntPoly::from_coeffs(vec![bi(-1), bi(-1), bi(0), bi(1)])
    }

    #[test]
    fn rational_factor_base_primes_up_to_10() {
        let f = f_cubic();
        let fb = FactorBase::new(&f, 10, 2);
        assert_eq!(fb.rational_primes, vec![2, 3, 5, 7]);
    }

    #[test]
    fn algebraic_index_lookup() {
        let f = f_cubic();
        // f(x) = x³ − x − 1; roots mod 2: f(0) = -1 ≡ 1, f(1) = -1 ≡ 1 — no roots mod 2.
        // roots mod 3: f(0) = -1 ≡ 2, f(1) = -1 ≡ 2, f(2) = 5 ≡ 2 — no roots mod 3.
        let fb = FactorBase::new(&f, 5, 5);
        // For each ideal in the base, the lookup should return its index.
        for ap in &fb.algebraic_ideals {
            let idx = fb.algebraic_index(ap.p, ap.r);
            assert_eq!(idx, Some(ap.index), "algebraic_index lookup failed for ({}, {})", ap.p, ap.r);
        }
    }

    #[test]
    fn rational_index_lookup() {
        let f = f_cubic();
        let fb = FactorBase::new(&f, 10, 2);
        assert_eq!(fb.rational_index(2), Some(0));
        assert_eq!(fb.rational_index(3), Some(1));
        assert_eq!(fb.rational_index(5), Some(2));
        assert_eq!(fb.rational_index(7), Some(3));
        assert_eq!(fb.rational_index(11), None);
    }

    #[test]
    fn matrix_width_is_rational_plus_algebraic_plus_obstruction() {
        let f = f_cubic();
        let fb = FactorBase::new(&f, 10, 10);
        assert_eq!(fb.matrix_width(), fb.rational_size() + fb.algebraic_size() + 1);
    }
}
