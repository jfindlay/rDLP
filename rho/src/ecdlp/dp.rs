//! Distinguished-point predicate and record type.
//!
//! A *distinguished point* (DP) is a walk state whose x-coordinate has at least
//! `theta` low-order zero bits.  The expected number of steps between DPs is
//! `2^theta`, so `theta` trades off memory (fewer DPs stored) against the
//! expected number of steps per walker before a DP is emitted.
//!
//! # van Oorschot–Wiener parallel rho
//!
//! Each walker thread emits a [`DpRecord`] whenever it lands on a distinguished
//! point.  The coordinator collects these records and detects collisions: two
//! records with the same full x-coordinate but different `(a, b)` pairs are the same point
//! (two walk trajectories merged) or opposite points, and either way `k` can be recovered.

use shared_field::Fp;

use crate::curve::AffinePoint;

// ── Distinguished-point predicate ────────────────────────────────────────────

/// Return `true` when the low `theta` bits of `x_low_word` are all zero.
///
/// A point is "distinguished" when its x-coordinate has at least `theta`
/// trailing zero bits.  The expected density of DPs in the walk is `2^{-theta}`.
///
/// # Arguments
///
/// * `x_low_word` — low 64-bit word of the point's x-coordinate.
/// * `theta` — number of required low-order zero bits (0 ≤ theta ≤ 63).
#[inline]
pub fn is_distinguished(x_low_word: u64, theta: u32) -> bool {
    if theta == 0 {
        // Every point is distinguished when theta = 0.
        return true;
    }
    // Mask the low `theta` bits; a DP has all of them zero.
    let mask = (1u64 << theta).wrapping_sub(1);
    (x_low_word & mask) == 0
}

// ── Distinguished-point record ────────────────────────────────────────────────

/// Where a distinguished point sits: enough to tell `P` from `−P`.
///
/// `x` is the full x-coordinate (all four words, not only the low word that decides whether the
/// point is distinguished); `y_odd` is the parity of the y-coordinate. For an odd field prime
/// and `y ≠ 0`, `y` and `p − y` have opposite parity, so two finite points with equal `x` are
/// equal when their parities agree and negatives of each other when they differ.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DpPoint {
    /// The point at infinity, its own negative.
    Infinity,
    /// A finite point `(x, y)`.
    Finite {
        /// Full x-coordinate, least-significant word first.
        x: [u64; 4],
        /// Whether the y-coordinate is odd.
        y_odd: bool,
    },
}

impl DpPoint {
    /// Describe an affine point.
    pub fn of<F: Fp<4>>(point: &AffinePoint<F>) -> Self {
        match point {
            AffinePoint::Infinity => DpPoint::Infinity,
            AffinePoint::Finite { x, y } => DpPoint::Finite {
                x: *x.to_uint().as_words(),
                y_odd: y.to_uint().as_words()[0] & 1 == 1,
            },
        }
    }

    /// The low word of the x-coordinate (`0` at infinity); the word the distinguished-point
    /// predicate looks at.
    pub fn x_low(&self) -> u64 {
        match self {
            DpPoint::Infinity => 0,
            DpPoint::Finite { x, .. } => x[0],
        }
    }

    /// The identity of the point up to sign: the full x-coordinate, `None` at infinity.
    pub fn x_key(&self) -> Option<[u64; 4]> {
        match self {
            DpPoint::Infinity => None,
            DpPoint::Finite { x, .. } => Some(*x),
        }
    }

    /// Whether `self` and `other` are negatives of each other (`P` and `−P`, `P ≠ −P`).
    ///
    /// Callers have already matched [`x_key`](Self::x_key); points at infinity are never opposite.
    pub fn is_opposite(&self, other: &Self) -> bool {
        match (self, other) {
            (DpPoint::Finite { y_odd: p, .. }, DpPoint::Finite { y_odd: q, .. }) => p != q,
            _ => false,
        }
    }
}

/// A distinguished point emitted by a walker thread.
///
/// Carries the walk's accumulated scalar coefficients `(a, b)` at the moment the DP was
/// detected and the point itself, with its full x-coordinate and y-parity, so the coordinator
/// can match points by x and tell equal points from opposite ones.
///
/// Invariant: the walk point satisfies `W = a·G + b·Q` at the time of emission.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DpRecord {
    /// The distinguished point `W`.
    pub point: DpPoint,
    /// Accumulated scalar coefficient for G: `W = a·G + b·Q`.
    pub a: u64,
    /// Accumulated scalar coefficient for Q: `W = a·G + b·Q`.
    pub b: u64,
    /// Index of the walker thread that produced this record.
    pub walk_id: usize,
}

impl DpRecord {
    /// Record the walk state `point = a·G + b·Q` of walker `walk_id`.
    pub fn new<F: Fp<4>>(point: &AffinePoint<F>, a: u64, b: u64, walk_id: usize) -> Self {
        DpRecord { point: DpPoint::of(point), a, b, walk_id }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn theta_zero_always_distinguished() {
        // Every x value is a DP when theta = 0.
        for x in [0u64, 1, 2, 3, 0xFFFF_FFFF_FFFF_FFFF] {
            assert!(is_distinguished(x, 0), "theta=0: x={x} should be distinguished");
        }
    }

    #[test]
    fn theta_one_even_only() {
        assert!(is_distinguished(0, 1));
        assert!(!is_distinguished(1, 1));
        assert!(is_distinguished(2, 1));
        assert!(!is_distinguished(3, 1));
    }

    #[test]
    fn theta_eight_low_byte_zero() {
        assert!(is_distinguished(0x0000_0000_0000_0100, 8));
        assert!(!is_distinguished(0x0000_0000_0000_0101, 8));
        assert!(!is_distinguished(0x0000_0000_0000_00FF, 8));
        assert!(is_distinguished(0x0000_0000_0000_0000, 8));
    }

    #[test]
    fn theta_63_only_zero() {
        assert!(is_distinguished(0, 63));
        assert!(!is_distinguished(1 << 62, 63));
    }
}
