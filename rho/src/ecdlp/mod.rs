//! ECDLP solver via Pollard rho.
//!
//! Optimization layers:
//! 4. r-adding walk (Teske, r≈20) + Brent's cycle detection.
//! 5. Distinguished points + parallel collision search (van Oorschot–Wiener).
//! 6. Negation map + fruitless-cycle escape (BKNS).
//! 7. Batched field inversion + affine coordinates.
//! 8. GLV endomorphism (order-3, secp-toy curve only).
//!
//! # r-adding walk solver
//!
//! [`solve_brent`] is the single-threaded entry point.  It runs Brent's cycle
//! detection on top of an r-adding walk, tracking `(a, b)` scalars throughout
//! so that when the tortoise and hare collide at the same point `W`, the DLP
//! `Q = k·G` can be recovered as `k = (a_t − a_h) / (b_h − b_t) mod n`.
//!
//! The solver retries with a fresh random walk table and starting point on
//! degenerate failures (`b_h = b_t mod n`), which occur with probability `~1/n`
//! per attempt. It validates the instance first and verifies every `k` it returns; the result is
//! a [`DlogSolution`] carrying [`RhoEvidence`], or a [`DlogError`] (see [`outcome`]).
//!
//! # Distinguished-point parallel solver
//!
//! [`solve_dp`] is the parallel entry point.  It spawns `num_walkers` threads,
//! each running an r-adding walk and emitting a [`dp::DpRecord`] whenever the
//! walk lands on a distinguished point (low `theta` bits of x-coordinate are
//! zero).  A coordinator on the calling thread collects DPs via a channel and
//! detects collisions that allow `k` to be recovered. The four distinguished-point solvers share
//! one total step budget (`max_steps`) that every walker claims from, so no input can make them
//! run unboundedly; each validates the instance, verifies the `k` a collision implies, and
//! returns the same [`DlogSolution`] / [`DlogError`] as [`solve_brent`]. A walker panic or a lost
//! channel is reported as [`DlogError::InvariantFailure`].

pub mod coordinator;
pub mod dp;
pub mod glv;
pub mod negmap;
pub mod outcome;
pub mod pohlig;
pub mod walk;

pub use outcome::{
    check_prime_order_input, verify_log, DlogError, DlogSolution, InputDefect, RhoEvidence,
};

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use crossbeam_channel::{bounded, Receiver, Sender};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use shared_field::Fp;

use crate::curve::{AffinePoint, Curve, JacobianPoint};
use coordinator::{Coordinator, DpCollision};
use dp::{DpPoint, DpRecord, is_distinguished};
use glv::glv_canonical;
use negmap::{FruitlessCycleDetector, canonical_rep, negate_scalars};
use walk::{AddendTable, AffineWalkState, BatchedWalker, WalkState};

// ── Modular arithmetic helpers ────────────────────────────────────────────────

/// Modular multiplication of two 64-bit values mod `n`.
#[inline]
fn mulmod64(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

/// Modular subtraction: `(a - b) mod n`, result in `[0, n)`.
#[inline]
fn submod64(a: u64, b: u64, n: u64) -> u64 {
    // `n - (b - a)` rather than `a + n - b`: the sum `a + n` overflows u64 when n > 2^63.
    if a >= b { a - b } else { n - (b - a) }
}

/// Modular inverse via Fermat's little theorem: `a^(n-2) mod n`.
///
/// Requires `n` to be prime and `a ≠ 0 mod n`.
///
/// # Panics
///
/// Panics if `a == 0`.
fn inv_mod_prime(a: u64, n: u64) -> u64 {
    assert!(a != 0, "inv_mod_prime: zero has no inverse");
    // Compute a^(n-2) mod n by iterative square-and-multiply.
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

// ── Collision extraction ─────────────────────────────────────────────────────

/// Attempt to extract `k` from a tortoise/hare collision.
///
/// When `W_tortoise = W_hare`, we have `a_t·G + b_t·Q = a_h·G + b_h·Q`, so
/// `(b_h − b_t)·Q = (a_t − a_h)·G`, giving `k = (a_t − a_h) / (b_h − b_t) mod n`.
///
/// Returns `None` if `b_h == b_t mod n` (degenerate — retry with new walk).
fn extract_k(a_t: u64, b_t: u64, a_h: u64, b_h: u64, n: u64) -> Option<u64> {
    let db = submod64(b_h, b_t, n); // b_h - b_t mod n
    if db == 0 {
        return None; // degenerate; caller should retry
    }
    let da = submod64(a_t, a_h, n); // a_t - a_h mod n
    let db_inv = inv_mod_prime(db, n);
    Some(mulmod64(da, db_inv, n))
}

// ── Brent's cycle detection ───────────────────────────────────────────────────

/// Largest Brent window an attempt may reach: an attempt that would double past this gives up,
/// after `2^29 − 1` walk steps.
const BRENT_MAX_WINDOW: u64 = 1 << 28;

/// How one Brent attempt ended.
enum BrentAttempt {
    /// The hare caught the tortoise; the scalars of both states at the collision.
    Collision { a_t: u64, b_t: u64, a_h: u64, b_h: u64 },
    /// The window grew past [`BRENT_MAX_WINDOW`] without a collision.
    StepCap,
}

/// Run one Brent attempt from the walk seeded by `rng_seed`.
///
/// Adds the hare steps taken to `steps`.
fn brent_attempt<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    rng_seed: u64,
    steps: &mut u64,
) -> BrentAttempt {
    // Fresh RNG per attempt so each retry explores a different trajectory.
    let mut rng = ChaCha20Rng::seed_from_u64(rng_seed);

    let table = AddendTable::new(curve, g, q, n, &mut rng);

    // Initialise both pointers from the same starting state (standard Brent).
    let start = WalkState::<F>::new_random(curve, g, q, n, &mut rng);
    let mut tortoise = start.clone();
    let mut hare = start;

    let mut r: u64 = 1; // current window size (power of 2)
    let mut count: u64 = 0; // steps taken in the current window

    let mut tortoise_affine = tortoise.to_affine(curve);

    loop {
        // Advance hare one step.
        hare.step(curve, &table, n);
        let hare_affine = hare.to_affine(curve);
        count += 1;
        *steps += 1;

        if hare_affine == tortoise_affine {
            return BrentAttempt::Collision {
                a_t: tortoise.a,
                b_t: tortoise.b,
                a_h: hare.a,
                b_h: hare.b,
            };
        }

        // If hare has taken `r` steps from the last tortoise snapshot:
        // snap tortoise forward to hare's current position and double the window.
        if count == r {
            tortoise = hare.clone();
            tortoise_affine = hare_affine;
            count = 0;
            r <<= 1;

            if r > BRENT_MAX_WINDOW {
                return BrentAttempt::StepCap;
            }
        }
    }
}

/// Return `k` as a solution after confirming `k·G == Q`.
fn verified<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    k: u64,
    evidence: RhoEvidence,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    if verify_log(curve, g, q, k) {
        Ok(DlogSolution::new(k, evidence))
    } else {
        Err(DlogError::Inconsistent { k })
    }
}

/// Solve `Q = k·G` via Pollard rho with Brent's cycle detection.
///
/// Validates the instance with [`check_prime_order_input`], answers `Q = ∞` (`k = 0`) and `Q = G`
/// (`k = 1`) as [`RhoEvidence::Trivial`] without a walk, then builds an r-adding walk table from
/// `g` and `q` and runs tortoise and hare pointers according to Brent's algorithm. When a
/// collision is detected, the DLP is extracted from the tracked `(a, b)` scalars and checked with
/// [`verify_log`] before it is returned.
///
/// Brent's algorithm for a pseudorandom sequence `x₀, x₁, x₂, …`:
///
/// - The tortoise is "frozen" at the start of each power-of-2 window.
/// - The hare advances one step at a time; at each step it is compared to the
///   frozen tortoise.
/// - After `r` comparisons without collision, the tortoise is snapped forward to
///   the hare's current position, and the window doubles: `r ← 2r`.
///
/// This terminates when the hare catches the tortoise on the cycle, using
/// ~24% fewer function evaluations than Floyd's algorithm on average.
///
/// A degenerate collision (`b_h ≡ b_t mod n`, no inverse) and an attempt that reaches the
/// per-attempt step cap (`2^29 − 1` steps) each consume one of the `max_retries` attempts; attempt
/// `i` uses seed `seed + i` (wrapping).
///
/// # Arguments
///
/// * `curve` — the curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — prime group order (64-bit).
/// * `seed` — RNG seed for reproducibility.
/// * `max_retries` — maximum number of fresh attempts (typical: 10).
///
/// # Errors
///
/// * [`DlogError::InvalidInput`] — the instance fails [`check_prime_order_input`].
/// * [`DlogError::BudgetExhausted`] — all `max_retries` attempts ended degenerate or at the step
///   cap (`max_retries = 0` gives `attempts: 0, steps: 0`). This proves nothing about `Q`.
/// * [`DlogError::Inconsistent`] — a non-degenerate collision gave a `k` with `k·G ≠ Q`. It is
///   never retried. Besides a walk bug, the one input that reaches it is a target in `E[n]` but
///   outside `⟨G⟩`, which the input check cannot exclude when `n²` divides the curve order.
pub fn solve_brent<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    seed: u64,
    max_retries: usize,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    check_prime_order_input(curve, g, q, n)?;

    if q.is_infinity() {
        return verified(curve, g, q, 0, RhoEvidence::Trivial);
    }
    if q == g {
        return verified(curve, g, q, 1, RhoEvidence::Trivial);
    }

    let mut steps: u64 = 0;
    for attempt in 0..max_retries {
        let rng_seed = seed.wrapping_add(attempt as u64);
        let BrentAttempt::Collision { a_t, b_t, a_h, b_h } =
            brent_attempt(curve, g, q, n, rng_seed, &mut steps)
        else {
            continue; // step cap reached
        };
        let Some(k) = extract_k(a_t, b_t, a_h, b_h, n) else {
            continue; // degenerate collision
        };
        let evidence = RhoEvidence::Collision {
            attempt,
            steps,
            a1: a_t,
            b1: b_t,
            a2: a_h,
            b2: b_h,
            negated: false,
        };
        return verified(curve, g, q, k, evidence);
    }

    Err(DlogError::BudgetExhausted { attempts: max_retries, steps })
}

// ── Distinguished-point search machinery ─────────────────────────────────────

/// Total walk-step budget shared by every walker of one search.
///
/// Walkers claim steps before taking them, so the steps taken across all walkers never exceed
/// `limit`, whatever the input.
struct StepBudget {
    limit: u64,
    used: AtomicU64,
}

impl StepBudget {
    fn new(limit: u64) -> Self {
        Self { limit, used: AtomicU64::new(0) }
    }

    /// Claim `count` steps; `false` (claiming nothing) when they would pass the limit.
    fn try_take(&self, count: u64) -> bool {
        self.used
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |used| {
                used.checked_add(count).filter(|total| *total <= self.limit)
            })
            .is_ok()
    }

    /// Steps claimed so far.
    fn used(&self) -> u64 {
        self.used.load(Ordering::Relaxed)
    }
}

/// Why a walker thread returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WalkerExit {
    /// The coordinator asked the walkers to stop.
    Stopped,
    /// The shared step budget has no room for another step.
    BudgetSpent,
    /// The coordinator's end of the channel vanished while it was still collecting.
    ChannelClosed,
}

/// A walker's connection to the search: the stop flag, the step budget and the DP channel.
struct WalkerLink<'a> {
    budget: &'a StepBudget,
    stop: &'a AtomicBool,
    tx: &'a Sender<DpRecord>,
}

impl WalkerLink<'_> {
    /// Claim `steps` before taking them; `Some(exit)` when the walker must return instead.
    fn claim(&self, steps: u64) -> Option<WalkerExit> {
        if self.stop.load(Ordering::Relaxed) {
            return Some(WalkerExit::Stopped);
        }
        if !self.budget.try_take(steps) {
            return Some(WalkerExit::BudgetSpent);
        }
        None
    }

    /// Send a distinguished point; `Some(exit)` when the coordinator is gone.
    fn emit(&self, dp: DpRecord) -> Option<WalkerExit> {
        self.tx.send(dp).err().map(|_| WalkerExit::ChannelClosed)
    }
}

/// Stops the other walkers if its thread unwinds, so one panic does not leave them spending the
/// whole budget before the failure is reported.
struct StopOnPanic<'a>(&'a AtomicBool);

impl Drop for StopOnPanic<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.store(true, Ordering::Relaxed);
        }
    }
}

/// Run `num_walkers` walker threads against a coordinator on the calling thread.
///
/// `walker(walk_id, link)` is the body of one walker. The search ends at the first usable
/// collision, or when every walker has returned because the `max_steps` budget is spent.
///
/// Returns the collision and the budget counter when it was accepted.
///
/// # Errors
///
/// * [`DlogError::BudgetExhausted`] — the walkers spent `max_steps` without a usable collision.
/// * [`DlogError::InvariantFailure`] — a walker panicked or lost its channel.
fn run_dp_search<W>(
    n: u64,
    num_walkers: usize,
    max_steps: u64,
    channel_capacity: usize,
    walker: W,
) -> Result<(DpCollision, u64), DlogError>
where
    W: Fn(usize, &WalkerLink<'_>) -> WalkerExit + Sync,
{
    let budget = StepBudget::new(max_steps);
    let stop = AtomicBool::new(false);
    let (tx, rx): (Sender<DpRecord>, Receiver<DpRecord>) = bounded(channel_capacity);

    let (hit, invariant) = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..num_walkers)
            .map(|walk_id| {
                let tx = tx.clone();
                let (budget, stop, walker) = (&budget, &stop, &walker);
                scope.spawn(move || {
                    let _stop_on_panic = StopOnPanic(stop);
                    walker(walk_id, &WalkerLink { budget, stop, tx: &tx })
                })
            })
            .collect();
        // Drop the original sender so the channel closes when every walker has returned.
        drop(tx);

        let mut coord = Coordinator::new(n);
        let mut hit = None;
        for dp in &rx {
            if let Some(collision) = coord.insert(dp) {
                hit = Some((collision, budget.used()));
                stop.store(true, Ordering::Relaxed);
                break;
            }
        }
        // Drain the channel so a walker mid-send is never left blocked.
        for _ in &rx {}

        let mut invariant = None;
        for handle in handles {
            match handle.join() {
                Ok(WalkerExit::Stopped | WalkerExit::BudgetSpent) => {}
                Ok(WalkerExit::ChannelClosed) => {
                    invariant.get_or_insert("a walker lost its channel to the coordinator");
                }
                Err(_) => {
                    invariant.get_or_insert("a walker thread panicked");
                }
            }
        }
        (hit, invariant)
    });

    if let Some(what) = invariant {
        return Err(DlogError::InvariantFailure(what));
    }
    hit.ok_or(DlogError::BudgetExhausted { attempts: 1, steps: budget.used() })
}

/// Validate the instance and answer the trivial targets (`Q = ∞`, `Q = G`) without a walk.
///
/// `Ok(None)` means a walk is needed.
fn dp_preflight<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Result<Option<DlogSolution<RhoEvidence>>, DlogError> {
    check_prime_order_input(curve, g, q, n)?;
    if q.is_infinity() {
        return verified(curve, g, q, 0, RhoEvidence::Trivial).map(Some);
    }
    if q == g {
        return verified(curve, g, q, 1, RhoEvidence::Trivial).map(Some);
    }
    Ok(None)
}

/// Reject a thread shape that cannot run a search.
fn check_walk_shape(num_walkers: usize, batch_size: Option<usize>) -> Result<(), DlogError> {
    if num_walkers == 0 {
        return Err(DlogError::Unsupported("at least one walker is required"));
    }
    if batch_size == Some(0) {
        return Err(DlogError::Unsupported("the batch size must be at least one"));
    }
    Ok(())
}

/// Turn a search hit into a solution, verifying the candidate `k` it implies.
///
/// Distinguished-point solvers report `attempt: 0`; `negated` is the collision's own sign (set
/// only by a plain walk meeting `−P`, since the canonicalising walks never emit both signs) and
/// `steps` is the shared counter when the collision was accepted.
fn dp_solution<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    (collision, steps): (DpCollision, u64),
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    let DpCollision { k, negated, a1, b1, a2, b2 } = collision;
    let evidence = RhoEvidence::Collision { attempt: 0, steps, a1, b1, a2, b2, negated };
    verified(curve, g, q, k, evidence)
}

// ── Distinguished-point parallel solver ──────────────────────────────────────

/// Solve `Q = k·G` via the van Oorschot–Wiener parallel rho algorithm.
///
/// Validates the instance with [`check_prime_order_input`] and answers `Q = ∞` and `Q = G` as
/// [`RhoEvidence::Trivial`] without a walk. Otherwise spawns `num_walkers` threads, each running
/// an independent r-adding walk. Walkers emit a [`DpRecord`] whenever the walk lands on a
/// distinguished point (low `theta` bits of the x-coordinate are zero). A coordinator on the
/// calling thread collects DPs and reports collisions; the solver verifies the `k` a collision
/// implies with [`verify_log`] before returning it.
///
/// Walker restart policy: if a walker takes more than `2^(theta + 10)` steps
/// without hitting a DP, it restarts from a fresh random `(a₀, b₀)`.  This
/// prevents walkers from being trapped in long fruitless cycles.
///
/// The walkers share one step budget, `max_steps`, claimed before each step: the call cannot run
/// longer than that, whatever the input.
///
/// The walk does not fold `P` and `−P` together, so two DPs can sit at opposite points. The
/// coordinator matches DPs by their full x-coordinate and reads the y-parity of each: equal
/// points give `a₁ + b₁·k = a₂ + b₂·k`, opposite points give `a₁ + b₁·k = −(a₂ + b₂·k)` and the
/// evidence is marked `negated: true`. Either way the recovered `k` is verified before it is
/// returned.
///
/// # Arguments
///
/// * `curve` — the curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — prime group order (64-bit).
/// * `num_walkers` — number of parallel walker threads (≥ 1).
/// * `theta` — DP threshold: a point is distinguished when its x-coordinate
///   has at least `theta` low-order zero bits.  Expected steps between DPs: `2^theta`.
/// * `seed` — RNG seed for reproducibility (each walker gets a derived seed).
/// * `max_steps` — total walk steps allowed across all walkers.
///
/// # Errors
///
/// * [`DlogError::InvalidInput`] — the instance fails [`check_prime_order_input`].
/// * [`DlogError::Unsupported`] — `num_walkers` is zero.
/// * [`DlogError::BudgetExhausted`] — `max_steps` were spent without a usable collision
///   (`attempts: 1`). This proves nothing about `Q`.
/// * [`DlogError::Inconsistent`] — a collision gave a `k` with `k·G ≠ Q`; never retried.
/// * [`DlogError::InvariantFailure`] — a walker panicked or lost its channel.
pub fn solve_dp<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    num_walkers: usize,
    theta: u32,
    seed: u64,
    max_steps: u64,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    if let Some(trivial) = dp_preflight(curve, g, q, n)? {
        return Ok(trivial);
    }
    check_walk_shape(num_walkers, None)?;

    // The addend table is built once and shared read-only across all walkers.
    let mut rng0 = ChaCha20Rng::seed_from_u64(seed);
    let table = AddendTable::<F>::new(curve, g, q, n, &mut rng0);

    // Capacity = 4 * num_walkers so walkers rarely block on a full channel while the
    // coordinator is processing a collision.
    let hit = run_dp_search(n, num_walkers, max_steps, 4 * num_walkers, |walk_id, link| {
        // Each walker gets a distinct seed so they explore different trajectories.
        let walker_seed = seed.wrapping_add(walk_id as u64 + 1);
        run_walker(curve, g, q, n, &table, theta, walk_id, walker_seed, link)
    })?;
    dp_solution(curve, g, q, hit)
}

// ── Walker thread body ────────────────────────────────────────────────────────

/// Body of a single walker thread.
///
/// Runs an r-adding walk, claiming each step from the shared budget and emitting a [`DpRecord`]
/// each time the walk lands on a distinguished point.  Restarts from a fresh random state
/// after `2^(theta + 10)` steps without a DP (dead-walk escape).
///
/// Returns when the coordinator stops the search, the budget is spent, or the DP channel closes.
fn run_walker<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    table: &AddendTable<F>,
    theta: u32,
    walk_id: usize,
    seed: u64,
    link: &WalkerLink<'_>,
) -> WalkerExit {
    // Dead-walk threshold: restart after this many steps without a DP.
    // 2^(theta + 10) gives ~1024 expected DPs per walk before restart.
    let max_steps_without_dp: u64 = 1u64 << (theta.min(53) + 10);

    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut walk = WalkState::<F>::new_random(curve, g, q, n, &mut rng);
    let mut steps_since_dp: u64 = 0;

    loop {
        if let Some(exit) = link.claim(1) {
            return exit;
        }

        walk.step(curve, table, n);
        steps_since_dp += 1;

        // Check for distinguished point.
        let pt = walk.to_affine(curve);
        let x_low = DpPoint::of(&pt).x_low();

        if is_distinguished(x_low, theta) {
            let dp = DpRecord::new(&pt, walk.a, walk.b, walk_id);
            if let Some(exit) = link.emit(dp) {
                return exit;
            }
            steps_since_dp = 0;
        }

        // Dead-walk escape: restart if we have gone too long without a DP.
        if steps_since_dp >= max_steps_without_dp {
            walk = WalkState::<F>::new_random(curve, g, q, n, &mut rng);
            steps_since_dp = 0;
        }
    }
}

// ── Negation-map parallel solver ──────────────────────────────────────────────

/// Solve `Q = k·G` via the van Oorschot–Wiener parallel rho algorithm with
/// the negation map (BKNS) applied.
///
/// Identical to [`solve_dp`] except each walker applies [`canonical_rep`] after
/// every step, adjusts `(a, b)` via [`negate_scalars`] when the negation fires,
/// and uses a [`FruitlessCycleDetector`] to escape period-2 cycles.
///
/// The negation map halves the effective group size, reducing the expected
/// number of steps to collision by a factor of ~√2.
///
/// # Walker count and the budget
///
/// Prefer at least 2 walkers.  With a single walker the coordinator can only find
/// a collision when the walk's rho trajectory has a non-empty "tail" (the walk
/// visits a DP before entering the cycle).  When the walk starts on the cycle
/// (probability ~1/2), no collision is found until the dead-walk restart, and a
/// run that never gets a tail ends as [`DlogError::BudgetExhausted`] when `max_steps` is
/// spent. Two walkers explore independent trajectories and can collide cross-walker even
/// when each individual walk is purely cyclic.
///
/// # Arguments
///
/// * `curve` — the curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — prime group order (64-bit).
/// * `num_walkers` — number of parallel walker threads (≥ 1; ≥ 2 recommended).
/// * `theta` — DP threshold: a point is distinguished when its x-coordinate
///   has at least `theta` low-order zero bits.
/// * `seed` — RNG seed for reproducibility.
/// * `max_steps` — total walk steps allowed across all walkers.
///
/// # Errors
///
/// As [`solve_dp`]. Negation-map walks compare canonical points, so every collision is between
/// identical points and a wrong `k` can only come from a walk bug.
pub fn solve_dp_negmap<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    num_walkers: usize,
    theta: u32,
    seed: u64,
    max_steps: u64,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    if let Some(trivial) = dp_preflight(curve, g, q, n)? {
        return Ok(trivial);
    }
    check_walk_shape(num_walkers, None)?;

    let mut rng0 = ChaCha20Rng::seed_from_u64(seed);
    let table = AddendTable::<F>::new(curve, g, q, n, &mut rng0);

    let hit = run_dp_search(n, num_walkers, max_steps, 4 * num_walkers, |walk_id, link| {
        let walker_seed = seed.wrapping_add(walk_id as u64 + 1);
        run_walker_negmap(curve, g, q, n, &table, theta, walk_id, walker_seed, link)
    })?;
    dp_solution(curve, g, q, hit)
}

// ── Negation-map walker thread body ──────────────────────────────────────────

/// Body of a single walker thread with the negation map applied.
///
/// After each walk step the current point is canonicalised via [`canonical_rep`].
/// When the negation fires, `(a, b)` are adjusted via [`negate_scalars`].
/// A [`FruitlessCycleDetector`] watches for period-2 cycles; when one is
/// detected the walk is perturbed by doubling the current point (BKNS escape):
///
/// ```text
/// W' = 2W,  a' = 2a mod n,  b' = 2b mod n
/// ```
///
/// then `canonical_rep` is applied to `W'` and the detector is reset.
///
/// Distinguished-point emission uses the canonical x-coordinate so that
/// collisions are between canonical representatives.
///
/// Steps are claimed from the shared budget one at a time, as in [`run_walker`].
fn run_walker_negmap<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    table: &AddendTable<F>,
    theta: u32,
    walk_id: usize,
    seed: u64,
    link: &WalkerLink<'_>,
) -> WalkerExit {
    let max_steps_without_dp: u64 = 1u64 << (theta.min(53) + 10);

    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut walk = WalkState::<F>::new_random(curve, g, q, n, &mut rng);
    let mut detector = FruitlessCycleDetector::new();
    let mut steps_since_dp: u64 = 0;

    // Apply the negation map to the initial point.
    {
        let pt = walk.to_affine(curve);
        let (canon, negated) = canonical_rep(&pt, &curve.p);
        if negated {
            let (na, nb) = negate_scalars(walk.a, walk.b, n);
            walk.a = na;
            walk.b = nb;
            walk.point_jac = JacobianPoint::from_affine(&canon, &curve.p);
        }
    }

    loop {
        if let Some(exit) = link.claim(1) {
            return exit;
        }

        walk.step(curve, table, n);
        steps_since_dp += 1;

        // Apply negation map to canonicalise the current point.
        let pt = walk.to_affine(curve);
        let (canon, negated) = canonical_rep(&pt, &curve.p);
        if negated {
            let (na, nb) = negate_scalars(walk.a, walk.b, n);
            walk.a = na;
            walk.b = nb;
            walk.point_jac = JacobianPoint::from_affine(&canon, &curve.p);
        }

        // Extract the canonical x-coordinate for DP detection and cycle tracking.
        let x_low = DpPoint::of(&canon).x_low();

        // Feed the cycle detector.
        detector.push(x_low);

        // BKNS escape: perturb by doubling when a fruitless 2-cycle is detected.
        // The escape does not emit a DP; the walk continues from the new point.
        let escaped = if detector.is_fruitless() {
            // W' = 2W; a' = 2a mod n; b' = 2b mod n.
            walk.point_jac = curve.double_jacobian(&walk.point_jac);
            walk.a = mulmod64(2, walk.a, n);
            walk.b = mulmod64(2, walk.b, n);

            // Re-canonicalise after the doubling.
            let pt2 = walk.to_affine(curve);
            let (canon2, negated2) = canonical_rep(&pt2, &curve.p);
            if negated2 {
                let (na, nb) = negate_scalars(walk.a, walk.b, n);
                walk.a = na;
                walk.b = nb;
                walk.point_jac = JacobianPoint::from_affine(&canon2, &curve.p);
            }

            detector.reset();
            true
        } else {
            false
        };

        // Distinguished-point check on the canonical x.
        // Skip DP emission immediately after a BKNS escape.
        if !escaped && is_distinguished(x_low, theta) {
            let dp = DpRecord::new(&canon, walk.a, walk.b, walk_id);
            if let Some(exit) = link.emit(dp) {
                return exit;
            }
            steps_since_dp = 0;
        }

        // Dead-walk escape: restart if too long without a DP.
        // Also fires if the BKNS escape has been running for too long without
        // the walk reaching a distinguished point.
        if steps_since_dp >= max_steps_without_dp {
            walk = WalkState::<F>::new_random(curve, g, q, n, &mut rng);
            // Re-canonicalise the fresh starting point.
            let pt_fresh = walk.to_affine(curve);
            let (canon_fresh, neg_fresh) = canonical_rep(&pt_fresh, &curve.p);
            if neg_fresh {
                let (na, nb) = negate_scalars(walk.a, walk.b, n);
                walk.a = na;
                walk.b = nb;
                walk.point_jac = JacobianPoint::from_affine(&canon_fresh, &curve.p);
            }
            detector.reset();
            steps_since_dp = 0;
        }
    }
}

// ── Batched-inversion parallel solver ────────────────────────────────────────

/// Solve `Q = k·G` via the van Oorschot–Wiener parallel rho algorithm with
/// batched field inversion.
///
/// Each thread owns a [`BatchedWalker`] with `batch_size` walk states.  After
/// each [`BatchedWalker::step_all`] call, each walk is checked for distinguished
/// points.  The negation map ([`canonical_rep`]) is applied after every step,
/// and a [`FruitlessCycleDetector`] is maintained per walk.
///
/// The batched-inversion optimization amortises the dominant per-step field inversion cost:
/// `batch_size` inversions are replaced by 1 inversion + 3(batch_size−1)
/// multiplications per `step_all` call.
///
/// A `step_all` call is `batch_size` walk steps and is claimed from the budget as a unit, so
/// the final `max_steps mod batch_size` steps of the budget go unused and a budget below
/// `batch_size` allows no walk at all.
///
/// # Arguments
///
/// * `curve` — the curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — prime group order (64-bit).
/// * `num_walkers` — number of parallel threads (≥ 1).
/// * `batch_size` — walks per thread (B ≥ 1; typical: 8–32).
/// * `theta` — DP threshold.
/// * `seed` — RNG seed for reproducibility.
/// * `max_steps` — total walk steps allowed across all walkers and batch members.
///
/// # Errors
///
/// As [`solve_dp`]; [`DlogError::Unsupported`] also covers `batch_size = 0`.
pub fn solve_dp_batch<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    num_walkers: usize,
    batch_size: usize,
    theta: u32,
    seed: u64,
    max_steps: u64,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    if let Some(trivial) = dp_preflight(curve, g, q, n)? {
        return Ok(trivial);
    }
    check_walk_shape(num_walkers, Some(batch_size))?;

    let mut rng0 = ChaCha20Rng::seed_from_u64(seed);
    let table = AddendTable::<F>::new(curve, g, q, n, &mut rng0);

    let capacity = 4 * num_walkers * batch_size;
    let hit = run_dp_search(n, num_walkers, max_steps, capacity, |walk_id_base, link| {
        let walker_seed = seed.wrapping_add(walk_id_base as u64 + 1);
        run_walker_batch(
            curve, g, q, n, &table, theta, walk_id_base, batch_size, walker_seed, link,
        )
    })?;
    dp_solution(curve, g, q, hit)
}

// ── Batched walker thread body ────────────────────────────────────────────────

/// Body of a single batched-walker thread.
///
/// Owns a [`BatchedWalker`] with `batch_size` walk states.  After each
/// [`BatchedWalker::step_all`], applies the negation map to each walk, checks
/// for distinguished points, and handles fruitless-cycle escape (BKNS doubling).
///
/// Walk IDs are `walk_id_base * batch_size + i` for walk index `i`, ensuring
/// globally unique IDs across threads.
///
/// Each `step_all` claims `batch_size` steps from the shared budget.
fn run_walker_batch<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    table: &AddendTable<F>,
    theta: u32,
    walk_id_base: usize,
    batch_size: usize,
    seed: u64,
    link: &WalkerLink<'_>,
) -> WalkerExit {
    let max_steps_without_dp: u64 = 1u64 << (theta.min(53) + 10);

    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    // Initialise B walk states.
    let walks: Vec<AffineWalkState<F>> = (0..batch_size)
        .map(|_| AffineWalkState::new_random(curve, g, q, n, &mut rng))
        .collect();

    let mut bw = BatchedWalker::new(walks, table.clone(), n);

    // Per-walk fruitless-cycle detectors and step counters.
    let mut detectors: Vec<negmap::FruitlessCycleDetector> =
        (0..batch_size).map(|_| negmap::FruitlessCycleDetector::new()).collect();
    let mut steps_since_dp: Vec<u64> = vec![0u64; batch_size];

    // Apply initial canonicalisation to each walk.
    for i in 0..batch_size {
        let (canon, negated) = canonical_rep(&bw.walks[i].point, &curve.p);
        if negated {
            let (na, nb) = negate_scalars(bw.walks[i].a, bw.walks[i].b, n);
            bw.walks[i].a = na;
            bw.walks[i].b = nb;
            bw.walks[i].point = canon;
        }
    }

    loop {
        if let Some(exit) = link.claim(batch_size as u64) {
            return exit;
        }

        // Advance all B walks by one step (batched inversion).
        bw.step_all(curve);

        // Post-step processing for each walk.
        for i in 0..batch_size {
            steps_since_dp[i] += 1;

            // Apply negation map.
            let (canon, negated) = canonical_rep(&bw.walks[i].point, &curve.p);
            if negated {
                let (na, nb) = negate_scalars(bw.walks[i].a, bw.walks[i].b, n);
                bw.walks[i].a = na;
                bw.walks[i].b = nb;
                bw.walks[i].point = canon.clone();
            }

            let x_low = DpPoint::of(&canon).x_low();

            detectors[i].push(x_low);

            // BKNS escape: perturb by doubling when a fruitless 2-cycle is detected.
            let escaped = if detectors[i].is_fruitless() {
                // W' = 2W; a' = 2a mod n; b' = 2b mod n.
                let pt_jac = JacobianPoint::from_affine(&bw.walks[i].point, &curve.p);
                let doubled = curve.double_jacobian(&pt_jac);
                bw.walks[i].a = mulmod64(2, bw.walks[i].a, n);
                bw.walks[i].b = mulmod64(2, bw.walks[i].b, n);

                let pt2 = doubled.to_affine(&curve.p);
                let (canon2, negated2) = canonical_rep(&pt2, &curve.p);
                if negated2 {
                    let (na, nb) = negate_scalars(bw.walks[i].a, bw.walks[i].b, n);
                    bw.walks[i].a = na;
                    bw.walks[i].b = nb;
                    bw.walks[i].point = canon2;
                } else {
                    bw.walks[i].point = pt2;
                }

                detectors[i].reset();
                true
            } else {
                false
            };

            // Distinguished-point check.
            if !escaped && is_distinguished(x_low, theta) {
                let walk_id = walk_id_base * batch_size + i;
                let dp = DpRecord::new(&canon, bw.walks[i].a, bw.walks[i].b, walk_id);
                if let Some(exit) = link.emit(dp) {
                    return exit;
                }
                steps_since_dp[i] = 0;
            }

            // Dead-walk escape: restart this walk if too long without a DP.
            if steps_since_dp[i] >= max_steps_without_dp {
                bw.walks[i] = AffineWalkState::new_random(curve, g, q, n, &mut rng);
                let (canon_fresh, neg_fresh) = canonical_rep(&bw.walks[i].point, &curve.p);
                if neg_fresh {
                    let (na, nb) = negate_scalars(bw.walks[i].a, bw.walks[i].b, n);
                    bw.walks[i].a = na;
                    bw.walks[i].b = nb;
                    bw.walks[i].point = canon_fresh;
                }
                detectors[i].reset();
                steps_since_dp[i] = 0;
            }
        }
    }
}

// ── GLV parallel solver ───────────────────────────────────────────────────────

/// Whether `curve` and `n` are the `secp_k1_toy` parameters that [`BETA`](secp_k1_toy::BETA) and
/// [`LAMBDA`](secp_k1_toy::LAMBDA) belong to.
///
/// Compares the field prime, the coefficients `a`, `b` and the group order, not any name. The
/// generator is free: φ acts as `λ` on every point of the prime-order group.
fn is_secp_k1_toy_group(curve: &Curve, n: u64) -> bool {
    use crate::curve::secp_k1_toy as toy;
    curve.p == crypto_bigint::Uint::<4>::from(toy::P)
        && curve.a == crypto_bigint::Uint::<4>::ZERO
        && curve.b == crypto_bigint::Uint::<4>::from(7u64)
        && n == toy::N
}

/// Solve `Q = k·G` via the van Oorschot–Wiener parallel rho algorithm with
/// the GLV endomorphism optimization.
///
/// **Restriction**: this solver uses the `secp_k1_toy` GLV constants (`BETA`, `LAMBDA`)
/// internally, so it applies only to that curve. The curve is recognised by its parameters
/// (`p`, `a`, `b` and the group order `n`); on any other curve the call returns
/// [`DlogError::Unsupported`] before anything else is checked.
///
/// Each walker applies [`glv_canonical`] after every step, collapsing the
/// 6-orbit `{W, φ(W), φ²(W), −W, −φ(W), −φ²(W)}` to a single canonical
/// representative.  This reduces the effective group size by a factor of 6, which heuristically
/// gives about a √6 speedup vs the plain walk and √3 vs the negation-map-only walk (expected,
/// not measured here).
///
/// The batched-inversion infrastructure is reused: each thread
/// owns a [`BatchedWalker`] with `batch_size` walk states, and a `step_all` claims `batch_size`
/// steps from the budget as a unit, as in [`solve_dp_batch`].
///
/// # Arguments
///
/// * `curve` — the curve definition (must be `secp_k1_toy`).
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — prime group order (must equal `secp_k1_toy::N`).
/// * `num_walkers` — number of parallel threads (≥ 1).
/// * `batch_size` — walks per thread (B ≥ 1; typical: 8–32).
/// * `theta` — DP threshold.
/// * `seed` — RNG seed for reproducibility.
/// * `max_steps` — total walk steps allowed across all walkers and batch members.
///
/// # Errors
///
/// As [`solve_dp_batch`], plus [`DlogError::Unsupported`] when `curve` and `n` are not the
/// `secp_k1_toy` group.
pub fn solve_dp_glv<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    num_walkers: usize,
    batch_size: usize,
    theta: u32,
    seed: u64,
    max_steps: u64,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    use crate::curve::secp_k1_toy::{BETA, LAMBDA};
    if !is_secp_k1_toy_group(curve, n) {
        return Err(DlogError::Unsupported(
            "the GLV constants belong to the secp_k1_toy curve and its group order",
        ));
    }
    solve_dp_glv_impl(
        curve, g, q, n, num_walkers, batch_size, theta, seed, max_steps, BETA, LAMBDA,
    )
}

/// Internal GLV solver that accepts explicit GLV constants.
///
/// Separated from [`solve_dp_glv`] so that tests can use small curves with
/// their own (beta, lambda) parameters without the secp_k1_toy restriction.
fn solve_dp_glv_impl<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    num_walkers: usize,
    batch_size: usize,
    theta: u32,
    seed: u64,
    max_steps: u64,
    beta: u64,
    lambda: u64,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    if let Some(trivial) = dp_preflight(curve, g, q, n)? {
        return Ok(trivial);
    }
    check_walk_shape(num_walkers, Some(batch_size))?;

    let mut rng0 = ChaCha20Rng::seed_from_u64(seed);
    let table = AddendTable::<F>::new(curve, g, q, n, &mut rng0);

    let capacity = 4 * num_walkers * batch_size;
    let hit = run_dp_search(n, num_walkers, max_steps, capacity, |walk_id_base, link| {
        let walker_seed = seed.wrapping_add(walk_id_base as u64 + 1);
        run_walker_glv(
            curve, g, q, n, &table, theta, walk_id_base, batch_size, walker_seed, link, beta,
            lambda,
        )
    })?;
    dp_solution(curve, g, q, hit)
}

// ── GLV walker thread body ────────────────────────────────────────────────────

/// Body of a single GLV-walker thread.
///
/// Owns a [`BatchedWalker`] with `batch_size` walk states.  After each
/// [`BatchedWalker::step_all`], applies [`glv_canonical`] to each walk,
/// adjusting `(a, b)` to track the canonical representative.
///
/// A [`FruitlessCycleDetector`] watches for period-2 cycles in the canonical
/// x-coordinate sequence.  The GLV orbit has size 6, so the canonical walk
/// can still exhibit period-2 fruitless cycles (when the walk oscillates
/// between two orbit members that map to the same canonical representative).
/// The BKNS doubling escape handles this case identically to the negmap walker.
///
/// Walk IDs are `walk_id_base * batch_size + i` for walk index `i`.
///
/// Each `step_all` claims `batch_size` steps from the shared budget.
#[allow(clippy::too_many_arguments)]
fn run_walker_glv<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    table: &AddendTable<F>,
    theta: u32,
    walk_id_base: usize,
    batch_size: usize,
    seed: u64,
    link: &WalkerLink<'_>,
    beta: u64,
    lambda: u64,
) -> WalkerExit {
    let max_steps_without_dp: u64 = 1u64 << (theta.min(53) + 10);
    let p_mod = &curve.p;

    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    // Initialise B walk states.
    let walks: Vec<AffineWalkState<F>> = (0..batch_size)
        .map(|_| AffineWalkState::new_random(curve, g, q, n, &mut rng))
        .collect();

    let mut bw = BatchedWalker::new(walks, table.clone(), n);

    // Per-walk fruitless-cycle detectors and step counters.
    let mut detectors: Vec<FruitlessCycleDetector> =
        (0..batch_size).map(|_| FruitlessCycleDetector::new()).collect();
    let mut steps_since_dp: Vec<u64> = vec![0u64; batch_size];

    // Apply initial GLV canonicalisation to each walk.
    for i in 0..batch_size {
        let (canon, new_a, new_b) = glv_canonical(
            &bw.walks[i].point, bw.walks[i].a, bw.walks[i].b, p_mod, n, beta, lambda,
        );
        bw.walks[i].point = canon;
        bw.walks[i].a = new_a;
        bw.walks[i].b = new_b;
    }

    loop {
        if let Some(exit) = link.claim(batch_size as u64) {
            return exit;
        }

        // Advance all B walks by one step (batched inversion).
        bw.step_all(curve);

        // Post-step processing for each walk.
        for i in 0..batch_size {
            steps_since_dp[i] += 1;

            // Apply GLV canonical map.
            let (canon, new_a, new_b) = glv_canonical(
                &bw.walks[i].point, bw.walks[i].a, bw.walks[i].b, p_mod, n, beta, lambda,
            );
            bw.walks[i].point = canon;
            bw.walks[i].a = new_a;
            bw.walks[i].b = new_b;

            let x_low = DpPoint::of(&bw.walks[i].point).x_low();

            detectors[i].push(x_low);

            // BKNS escape: perturb by doubling when a fruitless 2-cycle is detected.
            let escaped = if detectors[i].is_fruitless() {
                // W' = 2W; a' = 2a mod n; b' = 2b mod n.
                let pt_jac = JacobianPoint::from_affine(&bw.walks[i].point, p_mod);
                let doubled = curve.double_jacobian(&pt_jac);
                bw.walks[i].a = mulmod64(2, bw.walks[i].a, n);
                bw.walks[i].b = mulmod64(2, bw.walks[i].b, n);

                let pt2 = doubled.to_affine(p_mod);
                let (canon2, new_a2, new_b2) = glv_canonical(
                    &pt2, bw.walks[i].a, bw.walks[i].b, p_mod, n, beta, lambda,
                );
                bw.walks[i].point = canon2;
                bw.walks[i].a = new_a2;
                bw.walks[i].b = new_b2;

                detectors[i].reset();
                true
            } else {
                false
            };

            // Distinguished-point check on the canonical x.
            if !escaped && is_distinguished(x_low, theta) {
                let walk_id = walk_id_base * batch_size + i;
                // Not escaped, so the canonical point is still the walk's point.
                let dp = DpRecord::new(
                    &bw.walks[i].point, bw.walks[i].a, bw.walks[i].b, walk_id,
                );
                if let Some(exit) = link.emit(dp) {
                    return exit;
                }
                steps_since_dp[i] = 0;
            }

            // Dead-walk escape: restart this walk if too long without a DP.
            if steps_since_dp[i] >= max_steps_without_dp {
                bw.walks[i] = AffineWalkState::new_random(curve, g, q, n, &mut rng);
                let (canon_fresh, new_a_f, new_b_f) = glv_canonical(
                    &bw.walks[i].point, bw.walks[i].a, bw.walks[i].b, p_mod, n, beta, lambda,
                );
                bw.walks[i].point = canon_fresh;
                bw.walks[i].a = new_a_f;
                bw.walks[i].b = new_b_f;
                detectors[i].reset();
                steps_since_dp[i] = 0;
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::FpMonty4 as FpMonty;

    use crate::curve::test_curves::{tiny_a, TINY_A_N, tiny_glv, TINY_GLV_N, TINY_GLV_BETA, TINY_GLV_LAMBDA};
    use negmap::{canonical_rep, negate_scalars};

    // ── solve_brent tests (r-adding walk) ────────────────────────────────────

    /// Solve a known DLP on the 20-bit test curve A and verify the result.
    ///
    /// The solver does not need to find the specific k that was used to construct Q —
    /// any k' satisfying k'·G = Q is correct.  The assertion checks Q, not k.
    fn check_dlog(k_target: u64) {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));

        let solution = solve_brent(&curve, &g, &q, TINY_A_N, 0, 20)
            .unwrap_or_else(|e| panic!("solve_brent failed for k={k_target}: {e:?}"));
        let k = solution.k;

        let check = curve.scalar_mul(&g, &Uint::<4>::from(k));
        assert_eq!(check, q, "k·G ≠ Q for computed k={k} (expected k={k_target})");
        assert!(matches!(solution.evidence, RhoEvidence::Collision { .. }));
    }

    /// `Q = ∞` is the trivial logarithm 0; no walk runs.
    #[test]
    fn solve_brent_target_at_infinity_is_trivial_zero() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = AffinePoint::<FpMonty>::Infinity;

        let solution = solve_brent(&curve, &g, &q, TINY_A_N, 0, 20).expect("Q = ∞ must solve");
        assert_eq!((solution.k, solution.evidence), (0, RhoEvidence::Trivial));
    }

    /// `Q = G` is the trivial logarithm 1; no walk runs, even with no attempts allowed.
    #[test]
    fn solve_brent_target_equals_generator_is_trivial_one() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();

        let solution = solve_brent(&curve, &g, &g, TINY_A_N, 0, 0).expect("Q = G must solve");
        assert_eq!((solution.k, solution.evidence), (1, RhoEvidence::Trivial));
    }

    /// With no attempts allowed the budget is exhausted before any step is taken.
    #[test]
    fn solve_brent_zero_retries_exhausts_budget_immediately() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        assert_eq!(
            solve_brent(&curve, &g, &q, TINY_A_N, 0, 0),
            Err(DlogError::BudgetExhausted { attempts: 0, steps: 0 })
        );
    }

    /// An instance that fails validation is rejected before any walk, with its defect named.
    #[test]
    fn solve_brent_rejects_composite_order_before_walking() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        assert_eq!(
            solve_brent(&curve, &g, &q, TINY_A_N + 1, 0, 20),
            Err(DlogError::InvalidInput(InputDefect::OrderNotPrime))
        );
    }

    /// A candidate that fails `k·G == Q` is reported as `Inconsistent`, not returned.
    #[test]
    fn verified_reports_wrong_candidate_as_inconsistent() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

        assert!(verified(&curve, &g, &q, 7, RhoEvidence::Trivial).is_ok());
        assert_eq!(
            verified(&curve, &g, &q, 8, RhoEvidence::Trivial),
            Err(DlogError::Inconsistent { k: 8 })
        );
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

    #[test]
    fn solve_tiny_a_dlp_7() {
        check_dlog(7);
    }

    #[test]
    fn solve_tiny_a_dlp_42() {
        check_dlog(42);
    }

    #[test]
    fn solve_tiny_a_dlp_100() {
        check_dlog(100);
    }

    // ── solve_dp tests (distinguished-point search) ──────────────────────────

    /// Step budget for fixed-seed success runs on the 20-bit curves. A collision there takes on
    /// the order of 10³–10⁴ steps, so this leaves about three orders of magnitude of headroom
    /// and a failure to find one is a bug, not bad luck.
    const SUCCESS_BUDGET: u64 = 1 << 24;

    /// `a·G + b·Q` for scalars reduced mod n.
    fn combine(
        curve: &Curve,
        g: &AffinePoint<FpMonty>,
        q: &AffinePoint<FpMonty>,
        a: u64,
        b: u64,
    ) -> AffinePoint<FpMonty> {
        let ag = curve.scalar_mul(g, &Uint::<4>::from(a));
        let bq = curve.scalar_mul(q, &Uint::<4>::from(b));
        curve
            .add_mixed(&JacobianPoint::from_affine(&ag, &curve.p), &bq)
            .to_affine(&curve.p)
    }

    /// Recheck, from the evidence alone, that a distinguished-point solution is a collision
    /// `a1·G + b1·Q = (−1)^negated · (a2·G + b2·Q)` and that its `k` verifies. Only the plain
    /// walk (`may_be_negated`) can meet `−P`; the canonicalising walks never report `negated`.
    /// Returns the evidence step count.
    fn assert_dp_evidence(
        curve: &Curve,
        g: &AffinePoint<FpMonty>,
        q: &AffinePoint<FpMonty>,
        solution: &DlogSolution<RhoEvidence>,
        label: &str,
        may_be_negated: bool,
    ) -> u64 {
        let RhoEvidence::Collision { attempt, steps, a1, b1, a2, b2, negated } = solution.evidence
        else {
            panic!("{label}: expected collision evidence, got {:?}", solution.evidence);
        };
        assert_eq!(attempt, 0, "{label}: DP solvers report attempt 0");
        assert!(may_be_negated || !negated, "{label}: canonicalising walks never negate");
        assert!(steps > 0, "{label}: a collision takes at least one step");
        let first = combine(curve, g, q, a1, b1);
        let second = combine(curve, g, q, a2, b2);
        let second = if negated { curve.negate(&second) } else { second };
        assert_eq!(first, second, "{label}: collision identity does not hold");
        assert!(verify_log(curve, g, q, solution.k), "{label}: k={} fails verify_log", solution.k);
        steps
    }

    /// `(curve, G, Q = k·G)` on tiny_a.
    fn tiny_a_instance(k: u64) -> (Curve, AffinePoint<FpMonty>, AffinePoint<FpMonty>) {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k));
        (curve, g, q)
    }

    /// Solve a known 20-bit DLP on tiny_a with a single walker.
    ///
    /// Uses theta = 4 (1-in-16 DPs) so the coordinator sees enough DPs quickly
    /// on the small curve without excessive memory.
    #[test]
    fn solve_dp_tiny() {
        let (curve, g, q) = tiny_a_instance(77);

        let solution = solve_dp(&curve, &g, &q, TINY_A_N, 1, 4, 0xDEAD_BEEF, SUCCESS_BUDGET)
            .expect("solve_dp failed on tiny_a");

        let steps = assert_dp_evidence(&curve, &g, &q, &solution, "solve_dp_tiny", true);
        assert!(steps <= SUCCESS_BUDGET);
    }

    /// Solve a 20-bit DLP on tiny_a with 2 walkers.
    ///
    /// Uses theta = 4 (1-in-16 DPs) on the 20-bit test curve so the test
    /// completes quickly even in debug mode.  The two-walker setup exercises
    /// the parallel coordination path without the cost of a larger curve.
    #[test]
    fn solve_dp_parallel() {
        let (curve, g, q) = tiny_a_instance(500_000);

        let solution = solve_dp(&curve, &g, &q, TINY_A_N, 2, 4, 0xCAFE_BABE, SUCCESS_BUDGET)
            .expect("solve_dp failed on tiny_a with 2 walkers");

        assert_dp_evidence(&curve, &g, &q, &solution, "solve_dp_parallel", true);
    }

    /// The trivial targets are answered before any walk, even with no walkers and no budget.
    #[test]
    fn dp_solvers_answer_trivial_targets_without_a_walk() {
        let (curve, g, _) = tiny_a_instance(1);
        let inf = AffinePoint::<FpMonty>::Infinity;

        for (target, k) in [(&inf, 0u64), (&g, 1u64)] {
            let solutions = [
                solve_dp(&curve, &g, target, TINY_A_N, 0, 4, 0, 0),
                solve_dp_negmap(&curve, &g, target, TINY_A_N, 0, 4, 0, 0),
                solve_dp_batch(&curve, &g, target, TINY_A_N, 0, 0, 4, 0, 0),
            ];
            for solution in solutions {
                let solution = solution.expect("trivial target must solve");
                assert_eq!((solution.k, solution.evidence), (k, RhoEvidence::Trivial));
            }
        }
    }

    /// A malformed instance is rejected with its defect before any walker starts.
    #[test]
    fn dp_solvers_reject_invalid_input_before_walking() {
        let (curve, g, q) = tiny_a_instance(7);
        let composite = TINY_A_N + 1;
        let invalid = Err(DlogError::InvalidInput(InputDefect::OrderNotPrime));

        assert_eq!(solve_dp(&curve, &g, &q, composite, 1, 4, 0, SUCCESS_BUDGET), invalid);
        assert_eq!(solve_dp_negmap(&curve, &g, &q, composite, 2, 4, 0, SUCCESS_BUDGET), invalid);
        assert_eq!(
            solve_dp_batch(&curve, &g, &q, composite, 2, 4, 4, 0, SUCCESS_BUDGET),
            invalid
        );
    }

    /// A thread shape that cannot run a search is unsupported, not a panic or a hang.
    #[test]
    fn dp_solvers_reject_empty_thread_shapes() {
        let (curve, g, q) = tiny_a_instance(7);

        assert!(matches!(
            solve_dp(&curve, &g, &q, TINY_A_N, 0, 4, 0, SUCCESS_BUDGET),
            Err(DlogError::Unsupported(_))
        ));
        assert!(matches!(
            solve_dp_negmap(&curve, &g, &q, TINY_A_N, 0, 4, 0, SUCCESS_BUDGET),
            Err(DlogError::Unsupported(_))
        ));
        assert!(matches!(
            solve_dp_batch(&curve, &g, &q, TINY_A_N, 2, 0, 4, 0, SUCCESS_BUDGET),
            Err(DlogError::Unsupported(_))
        ));
    }

    /// A budget far below a collision's cost is spent exactly and reported as exhaustion, by
    /// every variant. One-step claims use the whole budget; batch claims use whole batches.
    #[test]
    fn dp_solvers_report_budget_exhaustion() {
        let (curve, g, q) = tiny_a_instance(77);
        let spent = |steps| Err(DlogError::BudgetExhausted { attempts: 1, steps });

        assert_eq!(solve_dp(&curve, &g, &q, TINY_A_N, 2, 4, 1, 10), spent(10));
        assert_eq!(solve_dp_negmap(&curve, &g, &q, TINY_A_N, 2, 4, 1, 10), spent(10));
        // Two claims of 4 fit in 10; a third does not.
        assert_eq!(solve_dp_batch(&curve, &g, &q, TINY_A_N, 2, 4, 4, 1, 10), spent(8));
        // A budget of zero allows no step at all.
        assert_eq!(solve_dp(&curve, &g, &q, TINY_A_N, 2, 4, 1, 0), spent(0));
    }

    /// A single negation-map walker may never see a collision; the budget, not a hang, ends it.
    #[test]
    fn single_negmap_walker_ends_by_budget_when_it_finds_nothing() {
        let (curve, g, q) = tiny_a_instance(77);
        // Theta above the walk length gives no DPs, so no collision can ever form.
        let r = solve_dp_negmap(&curve, &g, &q, TINY_A_N, 1, 40, 3, 5_000);
        assert_eq!(r, Err(DlogError::BudgetExhausted { attempts: 1, steps: 5_000 }));
    }

    // ── search machinery ─────────────────────────────────────────────────────

    /// The budget grants whole claims up to its limit and refuses the rest.
    #[test]
    fn step_budget_grants_claims_up_to_the_limit() {
        let budget = StepBudget::new(10);
        assert!(budget.try_take(4));
        assert!(budget.try_take(4));
        assert!(!budget.try_take(4), "a claim past the limit must be refused whole");
        assert_eq!(budget.used(), 8);
        assert!(budget.try_take(2));
        assert!(!budget.try_take(1));
        assert!(!budget.try_take(u64::MAX), "an overflowing claim must be refused");
        assert_eq!(budget.used(), 10);
    }

    /// A panicking walker is an invariant failure, and the others are stopped rather than left to
    /// spend the budget.
    #[test]
    fn walker_panic_is_an_invariant_failure() {
        let r = run_dp_search(1_048_051, 2, u64::MAX, 8, |walk_id, link| {
            if walk_id == 0 {
                panic!("walker failed on purpose");
            }
            loop {
                if let Some(exit) = link.claim(1) {
                    return exit;
                }
            }
        });
        assert_eq!(r, Err(DlogError::InvariantFailure("a walker thread panicked")));
    }

    /// A walker that loses its channel while the coordinator is collecting is an invariant failure.
    #[test]
    fn lost_channel_is_an_invariant_failure() {
        let r = run_dp_search(1_048_051, 2, 1_000, 8, |_, _| WalkerExit::ChannelClosed);
        assert_eq!(
            r,
            Err(DlogError::InvariantFailure("a walker lost its channel to the coordinator"))
        );
    }

    /// A real walker whose coordinator is gone returns `ChannelClosed` on its first DP.
    #[test]
    fn walker_reports_a_closed_channel() {
        let (curve, g, q) = tiny_a_instance(77);
        let mut rng = ChaCha20Rng::seed_from_u64(9);
        let table = AddendTable::<FpMonty>::new(&curve, &g, &q, TINY_A_N, &mut rng);
        let budget = StepBudget::new(1_000);
        let stop = AtomicBool::new(false);
        let (tx, rx) = bounded::<DpRecord>(1);
        drop(rx);
        let link = WalkerLink { budget: &budget, stop: &stop, tx: &tx };

        // Theta 0 makes every point distinguished, so the first step emits.
        let exit = run_walker(&curve, &g, &q, TINY_A_N, &table, 0, 0, 5, &link);
        assert_eq!(exit, WalkerExit::ChannelClosed);
    }

    // ── negation-map tests ────────────────────────────────────────────────────

    /// canonical_rep(P) and canonical_rep(−P) produce the same point.
    ///
    /// Verified for G and several multiples on tiny_a.
    #[test]
    fn negmap_canonical_rep() {
        let curve = tiny_a();
        let p = &curve.p;
        let g: AffinePoint<FpMonty> = curve.generator();

        for k in [1u64, 2, 7, 42, 100, 999] {
            let pt = curve.scalar_mul(&g, &Uint::<4>::from(k));
            let neg_pt = curve.negate(&pt);
            let (canon_pt,  _) = canonical_rep(&pt,     p);
            let (canon_neg, _) = canonical_rep(&neg_pt, p);
            assert_eq!(
                canon_pt, canon_neg,
                "canonical_rep(P) ≠ canonical_rep(−P) for k={k}"
            );
        }
    }

    /// negate_scalars gives (n−a, n−b) and the adjusted scalars still satisfy
    /// the walk invariant `a'·G + b'·Q = original_point`.
    #[test]
    fn negmap_scalar_adjust() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let k_target: u64 = 42;
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));
        let n = TINY_A_N;

        // Choose arbitrary (a, b) and compute W = a·G + b·Q.
        let a: u64 = 300;
        let b: u64 = 500;
        let ag = curve.scalar_mul(&g, &Uint::<4>::from(a));
        let bq = curve.scalar_mul(&q, &Uint::<4>::from(b));
        let w = curve.add_mixed(
            &JacobianPoint::from_affine(&ag, &curve.p),
            &bq,
        ).to_affine(&curve.p);

        // Negate scalars: (n−a)·G + (n−b)·Q should equal −W.
        let (na, nb) = negate_scalars(a, b, n);
        assert_eq!(na, n - a);
        assert_eq!(nb, n - b);

        let nag = curve.scalar_mul(&g, &Uint::<4>::from(na));
        let nbq = curve.scalar_mul(&q, &Uint::<4>::from(nb));
        let neg_w_reconstructed = curve.add_mixed(
            &JacobianPoint::from_affine(&nag, &curve.p),
            &nbq,
        ).to_affine(&curve.p);

        let neg_w_direct = curve.negate(&w);
        assert_eq!(
            neg_w_reconstructed, neg_w_direct,
            "negated scalars do not reconstruct −W"
        );
    }

    /// Verify that the negmap walk emits DPs within a bounded number of steps.
    ///
    /// This is a single-threaded simulation of the negmap walker.  It checks
    /// that the walk emits at least one DP within 10000 steps with theta=4.
    #[test]
    fn negmap_walk_emits_dps() {
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        use walk::AddendTable;

        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();
        let k_target: u64 = 77;
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));
        let n = TINY_A_N;
        let theta = 4u32;

        let mut rng = ChaCha20Rng::seed_from_u64(0xDEAD_BEEF);
        let table = AddendTable::<FpMonty>::new(&curve, &g, &q, n, &mut rng);
        let mut walk = walk::WalkState::<FpMonty>::new_random(&curve, &g, &q, n, &mut rng);
        let mut detector = negmap::FruitlessCycleDetector::new();

        // Apply initial canonicalization.
        {
            let pt = walk.to_affine(&curve);
            let (canon, negated) = canonical_rep(&pt, &curve.p);
            if negated {
                let (na, nb) = negate_scalars(walk.a, walk.b, n);
                walk.a = na;
                walk.b = nb;
                walk.point_jac = JacobianPoint::from_affine(&canon, &curve.p);
            }
        }

        let mut dp_count = 0usize;
        let max_steps = 10_000usize;

        for _ in 0..max_steps {
            walk.step(&curve, &table, n);

            let pt = walk.to_affine(&curve);
            let (canon, negated) = canonical_rep(&pt, &curve.p);
            if negated {
                let (na, nb) = negate_scalars(walk.a, walk.b, n);
                walk.a = na;
                walk.b = nb;
                walk.point_jac = JacobianPoint::from_affine(&canon, &curve.p);
            }

            let x_low = match &canon {
                AffinePoint::Infinity => 0u64,
                AffinePoint::Finite { x, .. } => x.to_uint().as_words()[0],
            };

            detector.push(x_low);

            if detector.is_fruitless() {
                // BKNS escape.
                walk.point_jac = curve.double_jacobian(&walk.point_jac);
                walk.a = mulmod64(2, walk.a, n);
                walk.b = mulmod64(2, walk.b, n);
                let pt2 = walk.to_affine(&curve);
                let (canon2, negated2) = canonical_rep(&pt2, &curve.p);
                if negated2 {
                    let (na, nb) = negate_scalars(walk.a, walk.b, n);
                    walk.a = na;
                    walk.b = nb;
                    walk.point_jac = JacobianPoint::from_affine(&canon2, &curve.p);
                }
                detector.reset();
            } else if is_distinguished(x_low, theta) {
                dp_count += 1;
            }
        }

        assert!(dp_count > 0, "negmap walk emitted no DPs in {max_steps} steps");
    }

    /// solve_dp_negmap solves a known 20-bit DLP on tiny_a.
    ///
    /// Uses 2 walkers so that cross-walker collisions are possible even when
    /// a single walker's rho trajectory has no tail DPs.
    #[test]
    fn solve_dp_negmap_tiny() {
        let (curve, g, q) = tiny_a_instance(77);

        let solution = solve_dp_negmap(&curve, &g, &q, TINY_A_N, 2, 4, 0xDEAD_BEEF, SUCCESS_BUDGET)
            .expect("solve_dp_negmap failed on tiny_a");

        assert_dp_evidence(&curve, &g, &q, &solution, "solve_dp_negmap_tiny", false);
    }

    /// solve_dp_batch solves a known 20-bit DLP on tiny_a.
    ///
    /// Uses 2 threads with batch_size=4 and theta=4.
    #[test]
    fn solve_dp_batch_tiny() {
        let (curve, g, q) = tiny_a_instance(77);

        let solution =
            solve_dp_batch(&curve, &g, &q, TINY_A_N, 2, 4, 4, 0xBA7C_5EED, SUCCESS_BUDGET)
                .expect("solve_dp_batch failed on tiny_a");

        assert_dp_evidence(&curve, &g, &q, &solution, "solve_dp_batch_tiny", false);
    }

    // ── GLV endomorphism tests ────────────────────────────────────────────────

    /// solve_dp_glv_impl solves a known DLP on the tiny GLV test curve.
    ///
    /// Uses the tiny GLV curve (n=1093) with 2 threads, batch_size=4, theta=2.
    /// Verifies that the evidence is a collision and that the returned k satisfies k·G = Q.
    ///
    /// The tiny GLV curve has the same structure as secp_k1_toy (y² = x³ + 7,
    /// a=0) but over a much smaller prime field, making the test fast.
    #[test]
    fn solve_dp_glv_tiny() {
        let curve = tiny_glv();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(77u64));

        let solution = solve_dp_glv_impl(
            &curve, &g, &q, TINY_GLV_N, 2, 4, 2, 0xBA7C_5EED, SUCCESS_BUDGET,
            TINY_GLV_BETA, TINY_GLV_LAMBDA,
        ).expect("solve_dp_glv failed on tiny_glv");

        assert_dp_evidence(&curve, &g, &q, &solution, "solve_dp_glv_tiny", false);
    }

    /// The GLV solver spends exactly the whole batches its budget allows and reports exhaustion.
    #[test]
    fn solve_dp_glv_impl_reports_budget_exhaustion() {
        let curve = tiny_glv();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(77u64));

        let r = solve_dp_glv_impl(
            &curve, &g, &q, TINY_GLV_N, 2, 4, 2, 0xBA7C_5EED, 10, TINY_GLV_BETA, TINY_GLV_LAMBDA,
        );
        assert_eq!(r, Err(DlogError::BudgetExhausted { attempts: 1, steps: 8 }));
    }

    /// `solve_dp_glv` is `Unsupported` off the `secp_k1_toy` group, checked by parameters: a curve
    /// with the same equation over another prime, and the right curve with another claimed order.
    #[test]
    fn solve_dp_glv_rejects_other_curves() {
        use crate::curve::secp_k1_toy::{secp_k1_toy, N as SECP_N};

        let curve = tiny_glv();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(77u64));
        assert!(matches!(
            solve_dp_glv(&curve, &g, &q, TINY_GLV_N, 2, 4, 2, 0, SUCCESS_BUDGET),
            Err(DlogError::Unsupported(_))
        ));

        let toy = secp_k1_toy();
        let tg: AffinePoint<FpMonty> = toy.generator();
        let tq = toy.scalar_mul(&tg, &Uint::<4>::from(5u64));
        assert!(matches!(
            solve_dp_glv(&toy, &tg, &tq, SECP_N - 2, 2, 4, 2, 0, SUCCESS_BUDGET),
            Err(DlogError::Unsupported(_))
        ));
        assert!(is_secp_k1_toy_group(&toy, SECP_N));
    }

    /// GLV's total walk steps stay within a no-regression bound of negmap's on the tiny GLV curve.
    ///
    /// Runs both solvers on the same DLP instance over many seeds and compares the step counts
    /// they report in their collision evidence. The test checks only that GLV's summed steps are
    /// at most 2× the negation-map-only walk's; the expected √3 reduction is a heuristic and is
    /// not asserted here.
    ///
    /// Uses the tiny GLV curve (n=1093) so the test completes quickly. The walkers race, so a
    /// single run's count varies; the sum over many seeds is compared with a generous margin.
    #[test]
    fn glv_fewer_steps_than_negmap() {
        let k_target: u64 = 500;
        let curve = tiny_glv();
        let g: AffinePoint<FpMonty> = curve.generator();
        let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));

        let theta = 2u32;
        let runs = 20u64;

        let mut negmap_steps = 0u64;
        let mut glv_steps = 0u64;

        for i in 0..runs {
            let seed = 0xDEAD_C0DE_u64.wrapping_add(i * 0x1234_5678);

            let nm = solve_dp_negmap(&curve, &g, &q, TINY_GLV_N, 2, theta, seed, SUCCESS_BUDGET)
                .expect("solve_dp_negmap failed in glv_fewer_steps_than_negmap");
            negmap_steps += assert_dp_evidence(&curve, &g, &q, &nm, "negmap", false);

            let glv = solve_dp_glv_impl(
                &curve, &g, &q, TINY_GLV_N, 2, 4, theta, seed, SUCCESS_BUDGET,
                TINY_GLV_BETA, TINY_GLV_LAMBDA,
            ).expect("solve_dp_glv failed in glv_fewer_steps_than_negmap");
            glv_steps += assert_dp_evidence(&curve, &g, &q, &glv, "glv", false);
        }

        assert!(
            glv_steps <= negmap_steps * 2,
            "GLV took {glv_steps} steps over {runs} runs against negmap's {negmap_steps}; \
             expected GLV ≤ 2× negmap"
        );
    }

    /// Stress test: solve_dp_negmap never returns a wrong k on a 20-bit DLP.
    ///
    /// Runs `solve_dp_negmap` across different seeds and k-targets,
    /// verifying every result.  The primary goal is correctness: the negation
    /// map must never produce a k that fails the `k·G = Q` check.
    ///
    /// Also runs `solve_dp` on the same instances to confirm both solvers agree
    /// on the answer (they may return different valid k values, but both must
    /// satisfy `k·G = Q`).
    #[test]
    fn solve_dp_negmap_stress_correctness() {
        let curve = tiny_a();
        let g: AffinePoint<FpMonty> = curve.generator();

        // Use several distinct (k_target, seed) pairs to exercise different
        // walk trajectories.
        let cases: &[(u64, u64)] = &[
            (77,      0xDEAD_0001),
            (314_159, 0xDEAD_0002),
            (500_000, 0xDEAD_0003),
            (999_999, 0xDEAD_0004),
            (42,      0xDEAD_0005),
        ];

        for &(k_target, seed) in cases {
            let q = curve.scalar_mul(&g, &Uint::<4>::from(k_target));

            // Verify solve_dp_negmap always returns a valid k.
            // Use 2 walkers: with a single walker the negmap walk can get stuck
            // in a cycle where the coordinator never finds a valid collision
            // (the walk starts on the cycle, so all DPs have the same scalars
            // on each visit).  Two walkers explore different trajectories and
            // can find cross-walker collisions.
            let negmap = solve_dp_negmap(&curve, &g, &q, TINY_A_N, 2, 4, seed, SUCCESS_BUDGET)
                .unwrap_or_else(|e| panic!("solve_dp_negmap k={k_target} seed={seed:#x}: {e:?}"));
            assert_dp_evidence(&curve, &g, &q, &negmap, "negmap stress", false);

            // Also verify solve_dp for the same instance.
            let plain = solve_dp(&curve, &g, &q, TINY_A_N, 1, 4, seed, SUCCESS_BUDGET)
                .unwrap_or_else(|e| panic!("solve_dp k={k_target} seed={seed:#x}: {e:?}"));
            assert_dp_evidence(&curve, &g, &q, &plain, "plain stress", true);
        }
    }
}
