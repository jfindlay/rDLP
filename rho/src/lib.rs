//! Pollard rho: integer factorization and ECDLP with all canonical optimizations,
//! plus eight algebraic ECDLP attack modules.
//!
//! # Structure
//!
//! ## Pollard-rho baseline
//!
//! - [`curve`] — elliptic curve group law, affine/Jacobian points, two concrete curves.
//! - [`factor`] — integer factorization rho (Floyd → Brent → batched GCD → multi-c).
//! - [`ecdlp`] — ECDLP rho (r-adding walk → DPs → negmap → batched inv → GLV).
//!
//! ## Binary-curve substrate (shared by GHS and Koblitz attacks)
//!
//! - [`binary_curve`] — binary curve `y²+xy=x³+ax²+b` group law in López–Dahab
//!   projective coordinates; point decompression via the half-trace.
//! - [`binary_ecdlp`] — Pollard-rho ECDLP solver over binary curves (`F2m`+`BinaryCurve`):
//!   r-adding walk, distinguished-point collision, linear recovery of `k`.
//! - [`hyperelliptic`] — hyperelliptic curve `y²+h(x)y=f(x)` over GF(2^m), genus,
//!   Mumford divisor `[u,v]` representation, validity predicate, divisor-from-points.
//!
//! ## Algebraic ECDLP attack modules
//!
//! - [`pairing`] — `F_{p^k}` extension-field arithmetic and bilinear pairings (Weil,
//!   Tate); MOV/Frey–Rück reduction via the reduced Tate pairing.
//! - [`ssa`] — Smart–Satoh–Araki p-adic attack on anomalous curves (`#E(F_p) = p`):
//!   Hensel lift, `log(1 + t)` formal-group logarithm, DLP recovery.
//! - [`ghs`] — GHS Weil-descent attack on binary elliptic curves: Artin–Schreier
//!   extension, Weil restriction of scalars, hyperelliptic curve extraction, transfer map.
//! - [`semaev`] — Semaev summation polynomials over a prime-field Weierstrass curve:
//!   `F_p[x]` univariate resultant, multivariate/symmetric-polynomial type `S_m`.
//! - [`index_calculus`] — Gaudry–Diem–Joux–Vitse index-calculus ECDLP solver over
//!   `E(F_p)`: factor base, prime-order subgroup, relation collection, Z/ℓZ linear
//!   algebra, DLP recovery.

pub mod binary_curve;
pub mod binary_ecdlp;
pub mod curve;
pub mod ecdlp;
pub mod factor;
pub mod ghs;
pub mod hyperelliptic;
pub mod index_calculus;
pub mod pairing;
pub mod semaev;
pub mod ssa;
