//! Pollard rho ECDLP solver over binary curves (`F2m` + `BinaryCurve`).
//!
//! This module implements the baseline Pollard-rho discrete-logarithm solver
//! for the binary elliptic curve group law.  It is a **parallel** implementation
//! over [`BinaryCurve`] + [`F2m`], mirroring the structure of `rho::ecdlp` (which
//! is `Fp<4>` + `Curve`-bound) but consuming `BinaryCurve`'s group law instead.
//!
//! # Algorithm
//!
//! The r-adding walk (Teske 1998) maintains the invariant:
//!
//! ```text
//! W = a·G + b·Q
//! ```
//!
//! where `G` is the curve generator, `Q = k·G` is the target point, and
//! `(a, b)` are tracked scalars in `ℤ/nℤ`.  When two walk states collide at
//! the same point (`a₁G+b₁Q = a₂G+b₂Q`), the DLP is recovered as:
//!
//! ```text
//! k = (a₁ − a₂) / (b₂ − b₁) mod n
//! ```
//!
//! provided `gcd(b₂ − b₁, n) = 1`.  Degenerate collisions (where the gcd > 1)
//! are discarded and the walk is restarted.
//!
//! # Walk design
//!
//! The addend table has `R = 16` entries.  Each entry `R[i] = αᵢ·G + βᵢ·Q`
//! is precomputed from random scalars `(αᵢ, βᵢ) ∈ ℤ/nℤ`.  The partition
//! function maps the current point's x-coordinate to an index in `[0, R)`.
//!
//! Brent's cycle detection is applied on top of the walk.
//!
//! # Distinguished points
//!
//! The [`solve`] entry point uses Brent's cycle detection (single-threaded).
//! The distinguished-point predicate (low `theta` bits of x-coordinate zero)
//! is available for future parallel extensions.
//!
//! # Modular arithmetic
//!
//! The linear recovery uses extended Euclidean GCD for the modular inverse,
//! which handles both prime and composite group orders.  This is necessary
//! because toy binary curves may have composite group orders.
//!
//! # Walk-state invariant
//!
//! **The invariant `W = a·G + b·Q` is the load-bearing correctness contract.**
//! It must hold at every step.  A wrong addend table or group-law bug shows
//! up as a recovered `k` with `k·G ≠ Q`, which the end-to-end KAT catches.
//!
//! # Baseline rule
//!
//! This module is the **baseline** that the Koblitz τ-orbit variant in
//! [`koblitz`] reads.  It must never be altered by the Koblitz variant — the
//! Koblitz variant adds a new function, not a replacement.  The baseline stays
//! intact for benchmark comparison.

pub mod koblitz;

use crypto_bigint::Uint;
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;

use shared_gf2m::F2m;

use crate::binary_curve::{BinaryAffinePoint, BinaryCurve};

// ── Constants ─────────────────────────────────────────────────────────────────

/// Number of addends in the r-adding walk table.
///
/// Teske recommends r ≈ 20 for large groups.  We use 16 (a power of 2) so
/// the partition function can use a bitmask instead of a modulo.
pub const R: usize = 16;

// ── Modular arithmetic helpers ────────────────────────────────────────────────

/// Add two scalars modulo `n`.
#[inline]
fn add_mod_n(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 + b as u128) % n as u128) as u64
}

/// Subtract two scalars modulo `n`: `(a − b) mod n`, result in `[0, n)`.
#[inline]
fn sub_mod_n(a: u64, b: u64, n: u64) -> u64 {
    if a >= b { a - b } else { a + n - b }
}

/// Multiply two scalars modulo `n`.
#[inline]
fn mul_mod_n(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

/// Extended Euclidean algorithm: returns `(gcd(a, b), x)` such that `a·x ≡ gcd(a,b) (mod b)`.
///
/// Uses the iterative Bezout identity computation.  All inputs are non-negative.
fn extended_gcd_iter(a: u64, b: u64) -> (u64, i128) {
    let (mut old_r, mut r) = (a as i128, b as i128);
    let (mut old_s, mut s) = (1i128, 0i128);

    while r != 0 {
        let q = old_r / r;
        let tmp_r = old_r - q * r;
        old_r = r;
        r = tmp_r;
        let tmp_s = old_s - q * s;
        old_s = s;
        s = tmp_s;
    }
    // old_r = gcd, old_s = Bezout coefficient for a
    (old_r as u64, old_s)
}

/// Modular inverse of `a` modulo `n` via extended GCD.
///
/// Returns `Some(x)` such that `a·x ≡ 1 (mod n)` if `gcd(a, n) = 1`.
/// Returns `None` if `gcd(a, n) > 1` (no inverse exists).
fn inv_mod_n(a: u64, n: u64) -> Option<u64> {
    if a == 0 {
        return None;
    }
    let (g, x) = extended_gcd_iter(a % n, n);
    if g != 1 {
        return None;
    }
    // x may be negative; reduce to [0, n).
    let result = ((x % n as i128) + n as i128) as u64 % n;
    Some(result)
}

/// Sample a uniform random scalar from `[1, n)`.
fn random_nonzero_scalar<R: Rng>(n: u64, rng: &mut R) -> u64 {
    loop {
        let v = rng.next_u64() % n;
        if v != 0 {
            return v;
        }
    }
}

// ── Addend table ──────────────────────────────────────────────────────────────

/// One entry in the r-adding table.
///
/// Stores the affine point `R = α·G + β·Q` together with the scalar pair `(α, β)`.
/// When a walk step selects entry `i`, the current point is updated as
/// `W ← W + R[i]`, and the tracked scalars as `a ← a + α mod n`, `b ← b + β mod n`.
///
/// Invariant: `point = α·G + β·Q`.
#[derive(Clone, Debug)]
pub struct BinaryAddend<F> {
    /// Affine point `α·G + β·Q`.
    pub point: BinaryAffinePoint<F>,
    /// Scalar coefficient for G.
    pub alpha: u64,
    /// Scalar coefficient for Q.
    pub beta: u64,
}

/// Precomputed table of `R` random addend points for the binary r-adding walk.
///
/// Built once per DLP instance; shared (read-only) across all walk states.
#[derive(Clone, Debug)]
pub struct BinaryAddendTable<F> {
    /// The `R` addend entries.
    pub entries: Vec<BinaryAddend<F>>,
}

impl<F: F2m<1>> BinaryAddendTable<F> {
    /// Build the addend table for the binary r-adding walk.
    ///
    /// Randomly samples `R` scalar pairs `(αᵢ, βᵢ)` from `[1, n)` and computes
    /// `R[i] = αᵢ·G + βᵢ·Q` for each.
    ///
    /// # Arguments
    ///
    /// * `curve` — the binary curve definition.
    /// * `g` — base point G in affine form.
    /// * `q` — target point Q in affine form.
    /// * `n` — group order (prime or composite).
    /// * `rng` — a random number generator.
    pub fn new<R: Rng>(
        curve: &BinaryCurve,
        g: &BinaryAffinePoint<F>,
        q: &BinaryAffinePoint<F>,
        n: u64,
        rng: &mut R,
    ) -> Self {
        let mut entries = Vec::with_capacity(R);

        for _ in 0..R {
            let alpha = random_nonzero_scalar(n, rng);
            let beta  = random_nonzero_scalar(n, rng);

            // R[i] = alpha·G + beta·Q
            let alpha_g = curve.scalar_mul(g, &Uint::<1>::from(alpha));
            let beta_q  = curve.scalar_mul(q, &Uint::<1>::from(beta));
            let point   = curve.add(&alpha_g, &beta_q);

            entries.push(BinaryAddend { point, alpha, beta });
        }

        BinaryAddendTable { entries }
    }

    /// Partition function: map an affine point's x-coordinate to table index `[0, R)`.
    ///
    /// Uses the low bits of the x-coordinate's `Uint<1>` representation.
    /// The x-coordinate is a GF(2^m) element stored as a polynomial bit-vector;
    /// we use the low 64-bit word modulo R as the partition index.
    #[inline]
    pub fn partition(&self, pt: &BinaryAffinePoint<F>) -> usize {
        match pt {
            BinaryAffinePoint::Infinity => 0,
            BinaryAffinePoint::Finite { x, .. } => {
                let x_uint = x.to_uint();
                let low_word = x_uint.as_words()[0];
                (low_word % R as u64) as usize
            }
        }
    }
}

// ── Walk state ────────────────────────────────────────────────────────────────

/// State of one binary r-adding walk instance.
///
/// **Invariant**: `point = a·G + b·Q` at all times.
///
/// This invariant is the load-bearing correctness contract.  A violation shows
/// up as a recovered `k` with `k·G ≠ Q` in the end-to-end KAT.
#[derive(Clone, Debug)]
pub struct BinaryWalkState<F> {
    /// Current walk point in affine coordinates.
    pub point: BinaryAffinePoint<F>,
    /// Scalar coefficient for G: `point = a·G + b·Q`.
    pub a: u64,
    /// Scalar coefficient for Q: `point = a·G + b·Q`.
    pub b: u64,
}

impl<F: F2m<1>> BinaryWalkState<F> {
    /// Initialise a walk from random starting scalars `a₀, b₀`.
    ///
    /// Sets `point = a₀·G + b₀·Q` and records `(a₀, b₀)`.
    pub fn new_random<R: Rng>(
        curve: &BinaryCurve,
        g: &BinaryAffinePoint<F>,
        q: &BinaryAffinePoint<F>,
        n: u64,
        rng: &mut R,
    ) -> Self {
        let a0 = random_nonzero_scalar(n, rng);
        let b0 = random_nonzero_scalar(n, rng);

        let a0_g  = curve.scalar_mul(g, &Uint::<1>::from(a0));
        let b0_q  = curve.scalar_mul(q, &Uint::<1>::from(b0));
        let point = curve.add(&a0_g, &b0_q);

        BinaryWalkState { point, a: a0, b: b0 }
    }

    /// Advance the walk by one step.
    ///
    /// Selects the addend `R[i]` based on the current point's x-coordinate,
    /// updates the affine point via `BinaryCurve::add`, and updates `(a, b)` mod n.
    ///
    /// After this call, the invariant `point = a·G + b·Q` is maintained.
    pub fn step(
        &mut self,
        curve: &BinaryCurve,
        table: &BinaryAddendTable<F>,
        n: u64,
    ) {
        let idx = table.partition(&self.point);
        let addend = &table.entries[idx];

        // W ← W + R[i]  (affine addition)
        self.point = curve.add(&self.point, &addend.point);

        // Update scalar coefficients mod n.
        self.a = add_mod_n(self.a, addend.alpha, n);
        self.b = add_mod_n(self.b, addend.beta, n);
    }
}

// ── Collision recovery ────────────────────────────────────────────────────────

/// Attempt to recover `k` from a collision between two walk states.
///
/// Given `a₁·G + b₁·Q = a₂·G + b₂·Q`, solves for `k = Q/G`:
/// `k = (a₁ − a₂) / (b₂ − b₁) mod n`.
///
/// Returns `None` if `gcd(b₂ − b₁, n) > 1` (degenerate; no information about `k`).
fn recover_k(a1: u64, b1: u64, a2: u64, b2: u64, n: u64) -> Option<u64> {
    let db = sub_mod_n(b2, b1, n); // b2 - b1 mod n
    if db == 0 {
        return None; // degenerate collision: b1 = b2 mod n
    }
    let db_inv = inv_mod_n(db, n)?; // None if gcd(db, n) > 1
    let da = sub_mod_n(a1, a2, n); // a1 - a2 mod n
    Some(mul_mod_n(da, db_inv, n))
}

// ── Brent's cycle detection ───────────────────────────────────────────────────

/// Solve `Q = k·G` on a binary curve via Pollard rho with Brent's cycle detection.
///
/// Builds an r-adding walk table from `g` and `q`, then runs tortoise and hare
/// pointers according to Brent's algorithm.  When a collision is detected,
/// the DLP is extracted from the tracked `(a, b)` scalars.
///
/// The solver retries up to `max_retries` times on degenerate failures (i.e., when
/// `gcd(b_h − b_t, n) > 1` at collision time).
///
/// # Walk-state invariant
///
/// At every step, `walk.point = walk.a · G + walk.b · Q`.  This invariant is
/// the load-bearing correctness contract.  A violation shows up as a recovered
/// `k` with `k·G ≠ Q` in the end-to-end KAT.
///
/// # Arguments
///
/// * `curve` — the binary curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — group order (prime or composite; 64-bit).
/// * `seed` — RNG seed for reproducibility.
/// * `max_retries` — maximum number of fresh attempts (typical: 30).
///
/// # Returns
///
/// `Some(k)` such that `k·G = Q`, or `None` if all retries were degenerate.
pub fn solve_brent<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    q: &BinaryAffinePoint<F>,
    n: u64,
    seed: u64,
    max_retries: usize,
) -> Option<u64> {
    // Special case: Q = ∞ means k = 0.
    if q.is_infinity() {
        return Some(0);
    }

    for attempt in 0..max_retries {
        let mut rng = ChaCha20Rng::seed_from_u64(seed.wrapping_add(attempt as u64));

        let table = BinaryAddendTable::new(curve, g, q, n, &mut rng);

        // Initialise both pointers from the same starting state (standard Brent).
        let start = BinaryWalkState::<F>::new_random(curve, g, q, n, &mut rng);
        let mut tortoise = start.clone();
        let mut hare     = start;

        let mut r: u64 = 1;     // current window size (power of 2)
        let mut count: u64 = 0; // steps taken in the current window

        loop {
            // Advance hare one step.
            hare.step(curve, &table, n);
            count += 1;

            // Check for collision: same affine point.
            if hare.point == tortoise.point {
                let ta = tortoise.a;
                let tb = tortoise.b;
                let ha = hare.a;
                let hb = hare.b;
                if let Some(k) = recover_k(ta, tb, ha, hb, n) {
                    return Some(k);
                }
                break; // degenerate — retry outer attempt loop
            }

            // If hare has taken `r` steps from the last tortoise snapshot:
            // snap tortoise forward to hare's current position and double the window.
            if count == r {
                tortoise = hare.clone();
                count = 0;
                r <<= 1;

                // Safety cap: give up after 2^28 total evaluations on this attempt.
                if r > (1 << 28) {
                    break;
                }
            }
        }
    }

    None
}

/// Solve `Q = k·G` on a binary curve (primary entry point).
///
/// Wraps [`solve_brent`] with a default seed and retry count.
///
/// # Arguments
///
/// * `curve` — the binary curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — group order.
///
/// # Returns
///
/// `Some(k)` such that `k·G = Q`, or `None` if the solver failed.
pub fn solve<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    q: &BinaryAffinePoint<F>,
    n: u64,
) -> Option<u64> {
    solve_brent(curve, g, q, n, 0, 50)
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_gf2m::F2mNaive;

    /// GF(2^4) irreducible: x⁴+x+1 = 0x13.
    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64)
    }

    /// Toy binary curve: y²+xy = x³+x²+1 over GF(2^4) with x⁴+x+1.
    ///
    /// Group order 4: G=(1,6), 2G=(0,1), 3G=(1,7), 4G=∞.
    fn toy_curve() -> BinaryCurve {
        BinaryCurve {
            poly: poly4(),
            a: Uint::<1>::ONE,
            b: Uint::<1>::ONE,
            n: Uint::<1>::from(4u64),
            gx: Uint::<1>::ONE,
            gy: Uint::<1>::from(6u64),
        }
    }

    /// Verify the walk invariant `W = a·G + b·Q` holds after construction
    /// and after a sequence of steps.
    #[test]
    fn walk_invariant_holds() {
        let curve = toy_curve();
        let g = curve.generator::<F2mNaive<1>>();
        let n = 4u64;

        // Use k_target = 2: Q = 2·G = (0, 1).
        let k_target: u64 = 2;
        let q = curve.scalar_mul(&g, &Uint::<1>::from(k_target));

        let mut rng = ChaCha20Rng::seed_from_u64(0xdeadbeef);
        let table = BinaryAddendTable::new(&curve, &g, &q, n, &mut rng);
        let mut walk = BinaryWalkState::<F2mNaive<1>>::new_random(&curve, &g, &q, n, &mut rng);

        for _ in 0..8 {
            // Reconstruct `a·G + b·Q` from current scalars.
            let ag = curve.scalar_mul(&g, &Uint::<1>::from(walk.a));
            let bq = curve.scalar_mul(&q, &Uint::<1>::from(walk.b));
            let reconstructed = curve.add(&ag, &bq);
            assert_eq!(
                walk.point, reconstructed,
                "walk invariant broken: point ≠ a·G + b·Q"
            );
            walk.step(&curve, &table, n);
        }
    }

    /// Partition function returns indices in `[0, R)`.
    #[test]
    fn partition_in_range() {
        let curve = toy_curve();
        let g = curve.generator::<F2mNaive<1>>();
        let n = 4u64;
        let q = curve.scalar_mul(&g, &Uint::<1>::from(2u64));

        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let table = BinaryAddendTable::<F2mNaive<1>>::new(&curve, &g, &q, n, &mut rng);
        let mut walk = BinaryWalkState::<F2mNaive<1>>::new_random(&curve, &g, &q, n, &mut rng);

        for _ in 0..20 {
            let idx = table.partition(&walk.point);
            assert!(idx < R, "partition index {idx} out of range [0, {R})");
            walk.step(&curve, &table, n);
        }
    }

    /// `inv_mod_n` returns the correct inverse for coprime inputs.
    #[test]
    fn inv_mod_n_correct() {
        // 3 * 3 = 9 ≡ 1 (mod 4) — wait, 9 mod 4 = 1. ✓
        assert_eq!(inv_mod_n(3, 4), Some(3));
        // 1 * 1 = 1 (mod 4). ✓
        assert_eq!(inv_mod_n(1, 4), Some(1));
        // 2 has no inverse mod 4 (gcd(2,4) = 2).
        assert_eq!(inv_mod_n(2, 4), None);
        // 5 * 5 = 25 ≡ 1 (mod 12). ✓
        assert_eq!(inv_mod_n(5, 12), Some(5));
        // 7 * 7 = 49 ≡ 1 (mod 12). ✓
        assert_eq!(inv_mod_n(7, 12), Some(7));
    }

    /// `recover_k` returns `None` for degenerate collisions (b1 = b2 mod n).
    #[test]
    fn recover_k_degenerate() {
        // b1 = b2 mod n → degenerate.
        assert_eq!(recover_k(1, 3, 2, 3, 4), None);
        // b2 - b1 = 2 mod 4 → gcd(2, 4) = 2 → no inverse.
        assert_eq!(recover_k(1, 1, 3, 3, 4), None);
    }

    /// `recover_k` returns the correct k for a synthetic collision.
    #[test]
    fn recover_k_synthetic() {
        // k = 3, n = 4.
        // Choose a1=3, b1=1, a2=0, b2=2.
        // db = b2-b1 = 2-1 = 1 mod 4; da = a1-a2 = 3-0 = 3 mod 4.
        // k = da * inv(db) = 3 * 1 = 3 mod 4. ✓
        let k = recover_k(3, 1, 0, 2, 4);
        assert_eq!(k, Some(3));
    }
}
