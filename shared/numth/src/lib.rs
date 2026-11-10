//! Number-theoretic primitives for the rDLP project.
//!
//! This crate provides three families of primitives consumed by NFS sieving
//! (`gnfs::sieve`), NFS-DL relation collection (`gnfs::dl`), and ECDLP index calculus
//! (`rho::index_calculus`):
//!
//! - **Primality testing** ([`prime`]): Miller–Rabin with deterministic witness
//!   sets, proven for ``n < 3,317,044,064,679,887,385,961,981`` (≈ ``2^{81.5}``,
//!   [`MILLER_RABIN_PROOF_BOUND`]).  Above that bound a ``true`` verdict is only
//!   "probable prime", with no proven error bound.
//!
//! - **B-smoothness detection** ([`smooth`]): trial-division factorisation over
//!   a prime factor base, returning a [`SmoothWitness`] that carries the
//!   complete factorisation and an unfactored cofactor.  Designed for all three
//!   downstream consumers (see the module-level note in [`smooth`]).
//!
//! - **Elliptic Curve Method** ([`ecm`]): Lenstra's ECM for integer
//!   factorization using Montgomery-form curves (Suyama parameterization).
//!   No other crate in the workspace calls it; it is exercised by its own unit tests.
//!   (Pohlig–Hellman in `rho::ecdlp` factors its u64 group order by trial division; ECM would
//!   be the natural step for larger orders, but none is wired.)
//!
//! # Re-exports
//!
//! The most commonly used items are re-exported at the crate root:
//!
//! ```rust
//! use shared_numth::{
//!     is_prime, miller_rabin, MILLER_RABIN_PROOF_BOUND,
//!     SmoothWitness, trial_smooth, factor_base_up_to,
//!     EcmError, EcmResult, ecm_one_curve, ecm_factor,
//! };
//! ```

pub mod ecm;
pub mod prime;
pub mod smooth;

pub use ecm::{ecm_factor, ecm_one_curve, EcmError, EcmResult};
pub use prime::{is_prime, miller_rabin, MILLER_RABIN_PROOF_BOUND};
pub use smooth::{factor_base_up_to, trial_smooth, SmoothWitness};
