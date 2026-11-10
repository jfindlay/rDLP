//! Integer factorization via Pollard rho.
//!
//! Optimization layers:
//! 1. Floyd's cycle detection — baseline.
//! 2. Brent's cycle detection — ~24% fewer f-evaluations.
//! 3. Montgomery batched GCD — accumulate products, single GCD per batch.
//! 4. Multi-c parallel restart — `rayon` over a set of `c` values.
//!
//! All layers share the pseudorandom function `f(x) = x² + c mod N`.
//! Public entry point: [`rho::factor`] (parallel, batched), which returns a verified
//! [`Factorization`] or a [`FactorError`] separating invalid input, primes, unsupported sizes,
//! spent budgets, and a candidate divisor that fails the divisor check.  The individual variants
//! take an explicit step bound and are also exported for benchmarking.

pub mod rho;

pub use rho::{
    brent, brent_batched, factor, floyd, FactorError, FactorEvidence, Factorization, MAX_SUPPORTED,
    PRIME_PROOF_BOUND,
};
