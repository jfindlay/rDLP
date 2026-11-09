//! Coordinator: receives distinguished points from walker threads and detects collisions.
//!
//! The coordinator owns a hash table keyed by the full x-coordinate of the distinguished
//! point (all four words; the low word alone is not an identity).  When two records with the
//! same x but different `(a, b)` pairs arrive, the two walk trajectories passed through the
//! same point or through opposite points `P` and `−P` — a collision that lets us recover `k`.
//! The y-parity in each record tells the two cases apart.
//!
//! # Collision recovery
//!
//! Equal points: if `W = a₁·G + b₁·Q = a₂·G + b₂·Q`, then
//! `(b₂ − b₁)·Q = (a₁ − a₂)·G`, so `k = (a₁ − a₂) / (b₂ − b₁) mod n`
//! (provided `b₂ ≠ b₁ mod n`).
//!
//! Opposite points: if `a₁·G + b₁·Q = −(a₂·G + b₂·Q)`, then
//! `(b₁ + b₂)·Q = −(a₁ + a₂)·G`, so `k = −(a₁ + a₂) / (b₁ + b₂) mod n`
//! (provided `b₁ + b₂ ≠ 0 mod n`).
//!
//! A degenerate collision (`b₁ = b₂ mod n` for equal points, `b₁ + b₂ = 0 mod n` for opposite
//! points) cannot yield `k`; the coordinator discards it and waits for the next collision.
//!
//! The coordinator reports the colliding scalars together with the candidate `k`
//! they imply, as a [`DpCollision`]; it never certifies `k`.  The caller checks
//! `k·G = Q` before trusting it.

use std::collections::HashMap;

use crate::ecdlp::dp::DpRecord;

// ── Collision evidence ───────────────────────────────────────────────────────

/// Two distinguished points that coincide, with the logarithm they imply.
///
/// Invariant: `a1·G + b1·Q = (−1)^negated · (a2·G + b2·Q)`, with `b1 ≠ b2 (mod n)` when
/// `negated` is false and `b1 + b2 ≠ 0 (mod n)` when it is true.  `k` is the candidate
/// `(a1 − a2) / (b2 − b1) mod n` (equal points) or `−(a1 + a2) / (b1 + b2) mod n` (opposite
/// points); it is unverified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DpCollision {
    /// Candidate logarithm implied by the collision (unverified).
    pub k: u64,
    /// Whether the two points are negatives of each other rather than equal.
    pub negated: bool,
    /// Scalar of `G` in the stored (first) record.
    pub a1: u64,
    /// Scalar of `Q` in the stored (first) record.
    pub b1: u64,
    /// Scalar of `G` in the incoming (second) record.
    pub a2: u64,
    /// Scalar of `Q` in the incoming (second) record.
    pub b2: u64,
}

// ── Coordinator ───────────────────────────────────────────────────────────────

/// Coordinator for the van Oorschot–Wiener parallel rho algorithm.
///
/// Maintains a hash table of seen distinguished points and detects collisions
/// that allow the discrete logarithm to be recovered.
pub struct Coordinator {
    /// Distinguished-point table, keyed by the full x-coordinate (`None` at infinity).
    table: HashMap<Option<[u64; 4]>, DpRecord>,
    /// Prime group order.
    n: u64,
}

impl Coordinator {
    /// Create a new coordinator for a group of order `n`.
    pub fn new(n: u64) -> Self {
        Coordinator { table: HashMap::new(), n }
    }

    /// Insert a distinguished point and report a usable collision with a stored one.
    ///
    /// Returns `Some(collision)` if the incoming record collides with a stored one and
    /// the relevant `b` combination is invertible mod `n`, so a candidate `k` exists.
    ///
    /// A *collision* is two records with the same full x-coordinate but different `(a, b)`
    /// pairs: the same point, or opposite points when their y-parities differ (the collision
    /// is then `negated`).  Records with equal low words but different full x-coordinates are
    /// different points and never collide.  Records from the same walker (same `walk_id`) at
    /// the same point with the same `(a, b)` are a *fruitless cycle* — the walk looped back
    /// to a DP it already visited — and are discarded without updating the table.
    ///
    /// # Arguments
    ///
    /// * `dp` — the incoming distinguished-point record.
    pub fn insert(&mut self, dp: DpRecord) -> Option<DpCollision> {
        use std::collections::hash_map::Entry;

        match self.table.entry(dp.point.x_key()) {
            Entry::Vacant(e) => {
                // First time we see this x-coordinate; store and wait.
                e.insert(dp);
                None
            }
            Entry::Occupied(e) => {
                let stored = e.get();

                // Fruitless cycle: same walker revisiting its own DP.
                // This happens when a walk loops without hitting a new DP.
                // Discard without updating the table so the stored record
                // (from a potentially different walker) remains available.
                if stored.walk_id == dp.walk_id
                    && stored.point == dp.point
                    && stored.a == dp.a
                    && stored.b == dp.b
                {
                    return None;
                }

                // Genuine collision: two distinct (a, b) pairs at the same x.
                let (a1, b1, a2, b2) = (stored.a, stored.b, dp.a, dp.b);
                let negated = stored.point.is_opposite(&dp.point);
                let k = if negated {
                    recover_k_negated(a1, b1, a2, b2, self.n)?
                } else {
                    recover_k(a1, b1, a2, b2, self.n)?
                };
                Some(DpCollision { k, negated, a1, b1, a2, b2 })
            }
        }
    }

    /// Return the number of distinct distinguished points currently stored.
    #[inline]
    pub fn table_len(&self) -> usize {
        self.table.len()
    }
}

// ── k recovery ────────────────────────────────────────────────────────────────

/// Attempt to recover `k` from a collision between two walk states.
///
/// Given `a₁·G + b₁·Q = a₂·G + b₂·Q`, solves for `k = Q/G`:
/// `k = (a₁ − a₂) / (b₂ − b₁) mod n`.
///
/// Returns `None` if `b₂ = b₁ mod n` (degenerate; no information about `k`).
fn recover_k(a1: u64, b1: u64, a2: u64, b2: u64, n: u64) -> Option<u64> {
    let db = submod64(b2, b1, n); // b2 - b1 mod n
    if db == 0 {
        return None; // degenerate collision
    }
    let da = submod64(a1, a2, n); // a1 - a2 mod n
    let db_inv = inv_mod_prime(db, n);
    Some(mulmod64(da, db_inv, n))
}

/// Attempt to recover `k` from a collision between a point and the negative of another.
///
/// Given `a₁·G + b₁·Q = −(a₂·G + b₂·Q)`, solves for `k = Q/G`:
/// `k = −(a₁ + a₂) / (b₁ + b₂) mod n`.
///
/// Returns `None` if `b₁ + b₂ = 0 mod n` (degenerate; no information about `k`).
fn recover_k_negated(a1: u64, b1: u64, a2: u64, b2: u64, n: u64) -> Option<u64> {
    let db = addmod64(b1, b2, n); // b1 + b2 mod n
    if db == 0 {
        return None; // degenerate collision
    }
    let neg_da = submod64(0, addmod64(a1, a2, n), n); // −(a1 + a2) mod n
    let db_inv = inv_mod_prime(db, n);
    Some(mulmod64(neg_da, db_inv, n))
}

// ── Modular arithmetic helpers ────────────────────────────────────────────────

#[inline]
fn addmod64(a: u64, b: u64, n: u64) -> u64 {
    // Widen so that `a + b` cannot overflow u64 when n > 2^63.
    ((a as u128 + b as u128) % n as u128) as u64
}

#[inline]
fn mulmod64(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

#[inline]
fn submod64(a: u64, b: u64, n: u64) -> u64 {
    // `n - (b - a)` rather than `a + n - b`: the sum `a + n` overflows u64 when n > 2^63.
    if a >= b { a - b } else { n - (b - a) }
}

/// Modular inverse via Fermat's little theorem: `a^(n-2) mod n`.
///
/// Requires `n` prime and `a ≠ 0 mod n`.
fn inv_mod_prime(a: u64, n: u64) -> u64 {
    debug_assert!(a != 0, "inv_mod_prime: zero has no inverse");
    let mut result: u64 = 1;
    let mut base: u64 = a % n;
    let mut exp: u64 = n - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mulmod64(result, base, n);
        }
        base = mulmod64(base, base, n);
        exp >>= 1;
    }
    result
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecdlp::dp::DpPoint;

    /// Small prime for arithmetic tests.
    const N: u64 = 1_048_051; // TINY_A_N

    /// A finite point whose x-coordinate has low word `x_low` and zero upper words.
    fn at(x_low: u64, y_odd: bool) -> DpPoint {
        DpPoint::Finite { x: [x_low, 0, 0, 0], y_odd }
    }

    #[test]
    fn no_collision_on_first_insert() {
        let mut coord = Coordinator::new(N);
        let dp = DpRecord { point: at(0x100, false), a: 1, b: 2, walk_id: 0 };
        assert!(coord.insert(dp).is_none());
    }

    #[test]
    fn collision_recovers_k() {
        // Construct a synthetic collision: two DPs at the same x_low with
        // known (a, b) pairs.  Verify that the recovered k satisfies the
        // linear equation.
        let mut coord = Coordinator::new(N);

        // Choose a1, b1, k such that a1 + k*b1 = a2 + k*b2 (mod N).
        // Equivalently: k = (a1 - a2) / (b2 - b1) mod N.
        let k_expected: u64 = 42;
        let a1: u64 = 100;
        let b1: u64 = 7;
        // a2 = a1 - k*(b2-b1) mod N; pick b2 = 11.
        let b2: u64 = 11;
        let db = submod64(b2, b1, N); // b2 - b1 = 4
        let a2 = submod64(a1, mulmod64(k_expected, db, N), N);

        let dp1 = DpRecord { point: at(0xABCD, false), a: a1, b: b1, walk_id: 0 };
        let dp2 = DpRecord { point: at(0xABCD, false), a: a2, b: b2, walk_id: 1 };

        assert!(coord.insert(dp1).is_none(), "first insert should not collide");
        let collision = coord.insert(dp2).expect("second insert should collide");
        assert_eq!(collision.k, k_expected, "recovered k mismatch");
        assert!(!collision.negated, "equal points are not a negated collision");
        // The evidence carries the colliding scalars in insertion order.
        assert_eq!(
            (collision.a1, collision.b1, collision.a2, collision.b2),
            (a1, b1, a2, b2)
        );
    }

    /// The recovery arithmetic does not overflow for a modulus above 2^63.
    #[test]
    fn collision_recovers_k_for_modulus_above_two_pow_63() {
        // 2^64 − 59 is the largest prime below 2^64.
        let n = u64::MAX - 58;
        let mut coord = Coordinator::new(n);
        let k_expected: u64 = n - 5;
        let (a1, b1, b2) = (3u64, n - 2, 4u64);
        // a2 = a1 − k·(b2 − b1) mod n, so that a1 + k·b1 = a2 + k·b2.
        let db = submod64(b2, b1, n);
        let a2 = submod64(a1, mulmod64(k_expected, db, n), n);

        coord.insert(DpRecord { point: at(0x1, false), a: a1, b: b1, walk_id: 0 });
        let collision = coord
            .insert(DpRecord { point: at(0x1, false), a: a2, b: b2, walk_id: 1 })
            .expect("collision");
        assert_eq!(collision.k, k_expected);
    }

    /// `submod64` does not overflow when `a + n` exceeds `u64::MAX`.
    #[test]
    fn submod64_handles_moduli_above_two_pow_63() {
        let n = u64::MAX - 58;
        assert_eq!(submod64(n - 2, n - 1, n), n - 1);
        assert_eq!(submod64(0, n - 1, n), 1);
        assert_eq!(submod64(5, 3, n), 2);
        assert_eq!(submod64(3, 3, n), 0);
    }

    /// Records for `P` and `−P` (equal x, opposite y-parity) with hand-chosen scalars give the
    /// `k` of `a1 + b1·k ≡ −(a2 + b2·k)` and are flagged as a negated collision.
    #[test]
    fn opposite_points_recover_k_and_are_negated() {
        let mut coord = Coordinator::new(N);
        let k_expected: u64 = 42;
        let (a1, b1, b2) = (100u64, 7u64, 11u64);
        // a2 = −a1 − (b1 + b2)·k mod n, so that a1 + b1·k = −(a2 + b2·k).
        let a2 = submod64(
            submod64(0, a1, N),
            mulmod64(addmod64(b1, b2, N), k_expected, N),
            N,
        );

        let dp1 = DpRecord { point: at(0xABCD, false), a: a1, b: b1, walk_id: 0 };
        let dp2 = DpRecord { point: at(0xABCD, true), a: a2, b: b2, walk_id: 1 };

        assert!(coord.insert(dp1).is_none());
        let collision = coord.insert(dp2).expect("opposite points must collide");
        assert_eq!(
            collision,
            DpCollision { k: k_expected, negated: true, a1, b1, a2, b2 }
        );
    }

    /// The opposite-point arithmetic does not overflow for a modulus above 2^63.
    #[test]
    fn opposite_points_recover_k_for_modulus_above_two_pow_63() {
        let n = u64::MAX - 58;
        let mut coord = Coordinator::new(n);
        let k_expected: u64 = n - 5;
        let (a1, b1, b2) = (n - 3, n - 2, n - 7);
        let a2 = submod64(
            submod64(0, a1, n),
            mulmod64(addmod64(b1, b2, n), k_expected, n),
            n,
        );

        coord.insert(DpRecord { point: at(0x1, true), a: a1, b: b1, walk_id: 0 });
        let collision = coord
            .insert(DpRecord { point: at(0x1, false), a: a2, b: b2, walk_id: 1 })
            .expect("collision");
        assert_eq!((collision.k, collision.negated), (k_expected, true));
    }

    /// `b1 + b2 ≡ 0` between opposite points carries no information about `k`.
    #[test]
    fn degenerate_opposite_collision_is_not_reported() {
        let mut coord = Coordinator::new(N);
        coord.insert(DpRecord { point: at(0x400, false), a: 5, b: 3, walk_id: 0 });
        let dp2 = DpRecord { point: at(0x400, true), a: 9, b: N - 3, walk_id: 1 };
        assert!(coord.insert(dp2).is_none());
    }

    /// Points that share their low x word but differ in the upper words are different points,
    /// whatever their y-parities.
    #[test]
    fn equal_low_word_with_different_full_x_does_not_collide() {
        let mut coord = Coordinator::new(N);
        let point = |hi: u64, y_odd| DpPoint::Finite { x: [0xABCD, hi, 0, 0], y_odd };

        assert!(coord
            .insert(DpRecord { point: point(1, false), a: 1, b: 2, walk_id: 0 })
            .is_none());
        // Same low word, different upper word and opposite parity: neither collision kind.
        assert!(coord
            .insert(DpRecord { point: point(2, true), a: 7, b: 9, walk_id: 1 })
            .is_none());
        assert_eq!(coord.table_len(), 2, "both points are stored separately");
        // The full x still matches when it is equal.
        let dp = DpRecord { point: point(1, false), a: 3, b: 5, walk_id: 2 };
        assert!(coord.insert(dp).is_some());
    }

    /// Distinguished points at infinity match each other but never a finite point at x = 0.
    #[test]
    fn infinity_collides_only_with_infinity() {
        let mut coord = Coordinator::new(N);
        coord.insert(DpRecord { point: DpPoint::Infinity, a: 5, b: 3, walk_id: 0 });
        let finite = DpRecord { point: at(0, false), a: 9, b: 4, walk_id: 1 };
        assert!(coord.insert(finite).is_none());
        let inf = DpRecord { point: DpPoint::Infinity, a: 9, b: 4, walk_id: 2 };
        let collision = coord.insert(inf).expect("infinity matches infinity");
        assert!(!collision.negated);
    }

    #[test]
    fn degenerate_collision_is_not_reported() {
        // b1 = b2 mod N → no k can be recovered.
        let mut coord = Coordinator::new(N);
        let dp1 = DpRecord { point: at(0x200, false), a: 5, b: 3, walk_id: 0 };
        let dp2 = DpRecord { point: at(0x200, false), a: 9, b: 3, walk_id: 1 };
        coord.insert(dp1);
        assert!(coord.insert(dp2).is_none(), "degenerate collision should return None");
    }

    #[test]
    fn fruitless_cycle_ignored() {
        // Same walker, same (a, b) — fruitless cycle; table unchanged.
        let mut coord = Coordinator::new(N);
        let dp = DpRecord { point: at(0x300, false), a: 1, b: 2, walk_id: 0 };
        coord.insert(dp.clone());
        assert_eq!(coord.table_len(), 1);
        // Re-inserting the identical record should not collide.
        assert!(coord.insert(dp).is_none());
        assert_eq!(coord.table_len(), 1);
    }
}
