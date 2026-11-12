//! General Number Field Sieve (GNFS): polynomial selection, sieving, and factorization.
//!
//! This crate implements the GNFS pipeline for integer factorization and discrete logarithm
//! computation. It builds on the number-field substrate in `shared-numfield` (polynomials,
//! number fields, resultants, Dedekind factorisation) and the number-theory utilities in
//! `shared-numth` (primality, smoothness).
//!
//! # Structure
//!
//! - [`polyselect`] — polynomial selection: base-m expansion, Murphy-E scoring, root sieve,
//!   Coppersmith multi-poly. Entry point: [`polyselect::select_base_m`].
//! - [`sieve`] — sieving substrate: two-sided factor bases, the rational norm and the algebraic
//!   value `F(a, b) = c·N(a − b·α)`, the bridge to `Uint<4>`, and the `Relation` type (whose
//!   integer exponents serve both the GF(2) factoring path and the NFS-DL relation code).
//! - [`filter`] — filtering substrate: sparse GF(2) matrix construction, relation-provenance
//!   map, and singleton removal. Entry points: [`filter::build_matrix`],
//!   [`filter::remove_singletons`].
//! - [`sqrt`] — square root stage: checked dependencies; the exact rational square root of the
//!   rational product; the verified algebraic square root (direct Newton lifting at an inert prime,
//!   split-prime lifting with exact sign recombination); checked congruence assembly and the
//!   bounded-retry `factor` driver with typed reports. Entry points: [`sqrt::rational_sqrt`],
//!   [`sqrt::algebraic_sqrt`], [`sqrt::attempt_factor`], [`sqrt::factor`].
//!
//! # Pedagogical intent
//!
//! This is a reference library: correctness and clarity take precedence over performance.
//! Each module is annotated with the mathematical background and the science↔engineering
//! tradeoffs that appear at toy vs. cryptographic scale.

pub mod polyselect;
pub mod sieve;
pub mod filter;
pub mod linalg;
pub mod sqrt;
pub mod dl;

pub use polyselect::{
    select_base_m, select_base_m_with_m, optimal_degree,
    BaseMGenerator, PolyGenerator, PolyPair, PolyPairError,
    score,
    root_sieve, rotate, RootSieveConfig, RootSieveGenerator,
    coppersmith_polys, coppersmith_best, CoppersmithConfig, CoppersmithGenerator,
};

pub use sieve::{
    AlgebraicPrime, FactorBase,
    algebraic_norm, norm_sign, norm_to_uint, rational_norm, NormBridgeError,
    ExponentVector, Relation, RelationError,
    line_sieve, LineSieveConfig,
    special_q_sieve, SpecialQConfig, SpecialQResult,
    lattice_sieve, LatticeBasis, LatticeSieveConfig, LatticeSieveResult,
};

pub use filter::{
    MatrixRow, SparseMatrix, EXCESS_FLOOR,
    build_matrix, remove_singletons,
    prune_cliques, merge_columns,
};

pub use sqrt::{
    AttemptOutcome, AttemptRecord, CongruenceError, CongruenceOutcome, FactorError, FactorOutcome,
    FactorReport, RationalRoot, RationalSqrtOutcome, VerifiedFactor, algebraic_sqrt,
    attempt_factor, factor, factor_from_congruence, rational_sqrt,
};

pub use dl::{
    compute_schirokauer, DLRelation, PrimeIdeal, SchirokauerError,
    augment_relation, collect_dl_relations, DLMatrix,
};
