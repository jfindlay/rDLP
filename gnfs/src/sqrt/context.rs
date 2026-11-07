//! Checked polynomial context, relation-factor boundary, outcomes, and resource limits for the
//! algebraic square root.
//!
//! The algebraic square root is only meaningful for a polynomial pair that satisfies a handful of
//! algebraic preconditions, and it can fail in several mathematically distinct ways. This module
//! owns both halves of that interface:
//!
//! - [`SqrtContext`] validates a polynomial `f`, root `m` and modulus `N` once, and carries the
//!   derived monic polynomial `F(t) = c^(d-1) f(t/c)`, `θ = c·α`, and the image `c·m mod N`.
//! - [`RelationFactors`] validates the relation pairs `(a, b)` against a context, rejecting
//!   factors that are not units of `A = Q[t]/(F)` before any expensive reconstruction.
//! - [`VerifiedRoot`] is a witness that can only be produced by an exact verifier living in this
//!   module's subtree; [`AlgebraicSqrtOutcome`] and the error types separate invalid input,
//!   unsupported domains, incidental factors, nonsquares, and resource exhaustion.
//! - The `*Limits` structs bound every search; none has a `Default`, so callers state budgets.
//!
//! See `gnfs/docs/PEDAGOGY.md` for the derivation of the unit criterion and the normalization.

use std::fmt;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use shared_numfield::{discriminant, IntPoly, NonresidueSearchLimits, RatPoly, RootSearchLimits};

use crate::polyselect::PolyPair;

pub mod normalize;

// ─── classification vocabulary ───────────────────────────────────────────────

/// Coarse classification of every way the square-root stage can fail to produce a verified root.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FailureClass {
    /// The caller supplied malformed data.
    InvalidInput,
    /// An internal identity that must hold was observed to fail.
    InvariantViolation,
    /// The rational product is not a perfect square.
    RationalNonSquare,
    /// The algebraic product is certified not to be a square.
    AlgebraicNonSquare,
    /// The input is well formed but outside the supported domain or map.
    Unsupported,
    /// A proper factor of `N` was found incidentally.
    IncidentalFactor,
    /// A resource limit was reached; this is never evidence about squareness.
    ResourceExhausted,
    /// The congruence of squares produced only a trivial gcd.
    TrivialCongruence,
}

/// Malformed input: something the caller must correct.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InvalidInput {
    /// `N <= 1`.
    ModulusNotGreaterThanOne {
        /// The offending modulus.
        n: BigInt,
    },
    /// `f` is zero or of degree 0.
    ConstantPolynomial,
    /// A [`PolyPair`]'s `degree` field disagrees with `f`.
    DegreeFieldMismatch {
        /// The pair's declared degree.
        field_degree: usize,
        /// The degree of `f`.
        poly_degree: usize,
    },
    /// A [`PolyPair`]'s rational polynomial is not `x - m`.
    RationalPolyNotLinear {
        /// `x - m`.
        expected: IntPoly,
        /// The pair's `g`.
        actual: IntPoly,
    },
    /// `f(m)` is not divisible by `N`.
    RootMismatch {
        /// `f(m)` reduced to the least nonnegative residue; nonzero.
        f_of_m_mod_n: BigInt,
    },
    /// `disc(F) = 0`: the polynomial has a repeated factor.
    RepeatedFactor,
    /// No relation pairs were supplied.
    EmptyRelationSet,
    /// A relation factor `a - b·m` is zero.
    ZeroRationalFactor {
        /// Position of the first offending pair.
        index: usize,
    },
    /// A relation factor `a - b·θ/c` is the zero element of `A`.
    ZeroAlgebraicFactor {
        /// Position of the first offending pair.
        index: usize,
    },
}

/// Well-formed input that the supported map cannot handle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UnsupportedDomain {
    /// `gcd(c, N) = N`: the leading coefficient vanishes modulo `N`.
    LeadingCoefficientNotUnit,
    /// A nonzero relation factor is a zero divisor of `A`, so it is not a unit.
    ZeroDivisorFactor {
        /// Position of the first offending pair.
        index: usize,
    },
    /// A denominator of the root has `gcd = N` with the modulus.
    MapDenominatorNotUnit,
}

/// Where an incidental factor of `N` was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IncidentalSource {
    /// `gcd(c, N)` is a proper divisor.
    LeadingCoefficient,
    /// A denominator of the algebraic root shares a proper divisor with `N`.
    MapDenominator,
}

/// A proper factorization `N = divisor · cofactor` found as a by-product, with its provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IncidentalFactor {
    divisor: BigInt,
    cofactor: BigInt,
    source: IncidentalSource,
}

impl IncidentalFactor {
    /// Record `divisor` as a proper factor of `n`.
    ///
    /// Self-certifying: returns `None` unless `1 < divisor < n` and `divisor | n`.
    #[must_use]
    pub fn new(n: &BigInt, divisor: BigInt, source: IncidentalSource) -> Option<Self> {
        if divisor <= BigInt::one() || divisor >= *n || !(n % &divisor).is_zero() {
            return None;
        }
        let cofactor = n / &divisor;
        Some(Self { divisor, cofactor, source })
    }

    /// The proper divisor of `N`.
    #[must_use]
    pub fn divisor(&self) -> &BigInt {
        &self.divisor
    }

    /// `N / divisor`.
    #[must_use]
    pub fn cofactor(&self) -> &BigInt {
        &self.cofactor
    }

    /// Where the factor was found.
    #[must_use]
    pub fn source(&self) -> IncidentalSource {
        self.source
    }
}

/// Why validated input was not accepted into a context or relation set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputRejection {
    /// Malformed input.
    Invalid(InvalidInput),
    /// Outside the supported domain.
    Unsupported(UnsupportedDomain),
    /// A proper factor of `N` was found.
    Incidental(IncidentalFactor),
}

impl InputRejection {
    /// The failure class of this rejection.
    #[must_use]
    pub fn class(&self) -> FailureClass {
        match self {
            Self::Invalid(_) => FailureClass::InvalidInput,
            Self::Unsupported(_) => FailureClass::Unsupported,
            Self::Incidental(_) => FailureClass::IncidentalFactor,
        }
    }
}

/// An internal identity that must hold failed; this is a bug, never a nonsquare.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvariantViolation {
    /// Which identity failed.
    pub what: &'static str,
}

/// Errors, as opposed to mathematical outcomes: invalid input or a violated invariant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SqrtError {
    /// Malformed input.
    Invalid(InvalidInput),
    /// A violated internal invariant.
    Invariant(InvariantViolation),
}

impl SqrtError {
    /// The failure class of this error.
    #[must_use]
    pub fn class(&self) -> FailureClass {
        match self {
            Self::Invalid(_) => FailureClass::InvalidInput,
            Self::Invariant(_) => FailureClass::InvariantViolation,
        }
    }
}

impl fmt::Display for InvalidInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ModulusNotGreaterThanOne { n } => write!(f, "modulus must exceed 1, got {n}"),
            Self::ConstantPolynomial => write!(f, "algebraic polynomial is zero or constant"),
            Self::DegreeFieldMismatch { field_degree, poly_degree } => write!(
                f,
                "degree field is {field_degree} but the polynomial has degree {poly_degree}"
            ),
            Self::RationalPolyNotLinear { expected, actual } => {
                write!(f, "rational polynomial must be {expected:?}, got {actual:?}")
            }
            Self::RootMismatch { f_of_m_mod_n } => {
                write!(f, "f(m) is {f_of_m_mod_n} modulo N, not 0")
            }
            Self::RepeatedFactor => write!(f, "polynomial has a repeated factor (disc = 0)"),
            Self::EmptyRelationSet => write!(f, "relation set is empty"),
            Self::ZeroRationalFactor { index } => {
                write!(f, "relation {index} has a zero rational factor a - b*m")
            }
            Self::ZeroAlgebraicFactor { index } => {
                write!(f, "relation {index} has a zero algebraic factor")
            }
        }
    }
}

impl fmt::Display for UnsupportedDomain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LeadingCoefficientNotUnit => {
                write!(f, "leading coefficient is divisible by the modulus")
            }
            Self::ZeroDivisorFactor { index } => {
                write!(f, "relation {index} has a nonzero algebraic factor that is not a unit")
            }
            Self::MapDenominatorNotUnit => {
                write!(f, "a root denominator is divisible by the modulus")
            }
        }
    }
}

impl fmt::Display for IncidentalFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "incidental factorization {} = {} * {} ({:?})",
            &self.divisor * &self.cofactor,
            self.divisor,
            self.cofactor,
            self.source
        )
    }
}

impl fmt::Display for InputRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(e) => write!(f, "invalid input: {e}"),
            Self::Unsupported(e) => write!(f, "unsupported domain: {e}"),
            Self::Incidental(e) => write!(f, "{e}"),
        }
    }
}

impl fmt::Display for InvariantViolation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "internal invariant violated: {}", self.what)
    }
}

impl fmt::Display for SqrtError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(e) => write!(f, "invalid input: {e}"),
            Self::Invariant(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for InvalidInput {}
impl std::error::Error for UnsupportedDomain {}
impl std::error::Error for IncidentalFactor {}
impl std::error::Error for InputRejection {}
impl std::error::Error for InvariantViolation {}
impl std::error::Error for SqrtError {}

// ─── checked polynomial/map context ──────────────────────────────────────────

/// A validated polynomial `f`, root `m`, and modulus `N`, with derived normalization data.
///
/// Invariants: `N > 1`; `deg f = d >= 1`; `f(m) ≡ 0 (mod N)`; `c = lc(f)` is nonzero with
/// `gcd(c, N) = 1`; and `F(t) = c^(d-1) f(t/c)` is monic and squarefree. `f` is stored as given
/// (`c` may be negative) and `m` is not reduced modulo `N`. A squarefree `F` may be reducible, in
/// which case `A = Q[t]/(F)` is a product of number fields, not a field; the context makes no
/// irreducibility claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SqrtContext {
    n: BigInt,
    m: BigInt,
    f: IntPoly,
    c: BigInt,
    monic: IntPoly,
    degree: usize,
    theta_image: BigInt,
    disc: BigInt,
}

impl SqrtContext {
    /// Validate `f`, `m`, and `N` and derive the normalization data.
    ///
    /// Checks stop at the first failure, in order: `N > 1`; `deg f >= 1`; `f(m) ≡ 0 (mod N)`;
    /// `F` squarefree; `gcd(|c|, N)`. All invalid-input checks precede the domain checks, so a
    /// malformed input is never reported as an incidental factor.
    ///
    /// # Errors
    ///
    /// - [`InvalidInput::ModulusNotGreaterThanOne`], [`InvalidInput::ConstantPolynomial`],
    ///   [`InvalidInput::RootMismatch`], or [`InvalidInput::RepeatedFactor`] for malformed input.
    /// - [`InputRejection::Incidental`] with [`IncidentalSource::LeadingCoefficient`] when
    ///   `1 < gcd(c, N) < N`.
    /// - [`UnsupportedDomain::LeadingCoefficientNotUnit`] when `N | c`.
    pub fn new(f: &IntPoly, m: &BigInt, n: &BigInt) -> Result<Self, InputRejection> {
        if *n <= BigInt::one() {
            return Err(invalid(InvalidInput::ModulusNotGreaterThanOne { n: n.clone() }));
        }
        let degree = match f.degree() {
            Some(d) if d >= 1 => d,
            _ => return Err(invalid(InvalidInput::ConstantPolynomial)),
        };
        // eval_mod reduces with mod_floor at every step, so the result is the least nonnegative
        // residue even for negative `m` or coefficients.
        let residue = f.eval_mod(m, n);
        if !residue.is_zero() {
            return Err(invalid(InvalidInput::RootMismatch { f_of_m_mod_n: residue }));
        }
        let c = f.coeffs[degree].clone();
        let monic = monic_form(f, &c, degree);
        let disc = discriminant(&monic);
        if disc.is_zero() {
            return Err(invalid(InvalidInput::RepeatedFactor));
        }
        let g0 = c.gcd(n);
        if g0 == *n {
            return Err(InputRejection::Unsupported(UnsupportedDomain::LeadingCoefficientNotUnit));
        }
        if !g0.is_one() {
            // 1 < g0 < n by the two comparisons above, so the constructor cannot refuse.
            if let Some(found) = IncidentalFactor::new(n, g0, IncidentalSource::LeadingCoefficient)
            {
                return Err(InputRejection::Incidental(found));
            }
        }
        let theta_image = (&c * m).mod_floor(n);
        Ok(Self { n: n.clone(), m: m.clone(), f: f.clone(), c, monic, degree, theta_image, disc })
    }

    /// Validate a [`PolyPair`] by reading its public fields.
    ///
    /// Never calls `PolyPair::verify`, `new`, `number_field`, or `monic_f`: `verify` panics for
    /// `N = 0`, so `N > 1` is checked first. After `N > 1` and `deg f >= 1`, the pair's `degree`
    /// field and `g = x - m` are checked, then the remaining checks of [`SqrtContext::new`].
    ///
    /// # Errors
    ///
    /// As [`SqrtContext::new`], plus [`InvalidInput::DegreeFieldMismatch`] and
    /// [`InvalidInput::RationalPolyNotLinear`].
    pub fn from_poly_pair(pair: &PolyPair) -> Result<Self, InputRejection> {
        if pair.n <= BigInt::one() {
            return Err(invalid(InvalidInput::ModulusNotGreaterThanOne { n: pair.n.clone() }));
        }
        let poly_degree = match pair.f.degree() {
            Some(d) if d >= 1 => d,
            _ => return Err(invalid(InvalidInput::ConstantPolynomial)),
        };
        if pair.degree != poly_degree {
            return Err(invalid(InvalidInput::DegreeFieldMismatch {
                field_degree: pair.degree,
                poly_degree,
            }));
        }
        let expected = IntPoly::from_coeffs(vec![-pair.m.clone(), BigInt::one()]);
        if pair.g != expected {
            return Err(invalid(InvalidInput::RationalPolyNotLinear {
                expected,
                actual: pair.g.clone(),
            }));
        }
        Self::new(&pair.f, &pair.m, &pair.n)
    }

    /// The modulus `N > 1`.
    #[must_use]
    pub fn modulus(&self) -> &BigInt {
        &self.n
    }

    /// The root `m`, exactly as given (not reduced modulo `N`).
    #[must_use]
    pub fn root(&self) -> &BigInt {
        &self.m
    }

    /// The algebraic polynomial `f`, as given.
    #[must_use]
    pub fn algebraic_poly(&self) -> &IntPoly {
        &self.f
    }

    /// The leading coefficient `c` of `f`, nonzero and coprime to `N`.
    #[must_use]
    pub fn leading_coeff(&self) -> &BigInt {
        &self.c
    }

    /// The monic polynomial `F(t) = c^(d-1) f(t/c)`, the minimal polynomial of `θ = c·α`.
    #[must_use]
    pub fn monic_poly(&self) -> &IntPoly {
        &self.monic
    }

    /// The degree `d >= 1`.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The image `c·m mod N` of `θ`, in `[0, N)`.
    #[must_use]
    pub fn theta_image(&self) -> &BigInt {
        &self.theta_image
    }

    /// The nonzero discriminant of `F`.
    #[must_use]
    pub fn discriminant(&self) -> &BigInt {
        &self.disc
    }
}

fn invalid(e: InvalidInput) -> InputRejection {
    InputRejection::Invalid(e)
}

/// `F(t) = c^(d-1) f(t/c)`, with coefficient `f_k · c^(d-1-k)` for `k < d` and leading 1.
fn monic_form(f: &IntPoly, c: &BigInt, degree: usize) -> IntPoly {
    let powers = powers_of(c, degree);
    let mut coeffs: Vec<BigInt> =
        (0..degree).map(|k| &f.coeffs[k] * &powers[degree - 1 - k]).collect();
    coeffs.push(BigInt::one());
    IntPoly::from_coeffs(coeffs)
}

/// `[x^0, x^1, ..., x^top]`.
fn powers_of(x: &BigInt, top: usize) -> Vec<BigInt> {
    let mut powers = Vec::with_capacity(top + 1);
    powers.push(BigInt::one());
    for j in 0..top {
        powers.push(&powers[j] * x);
    }
    powers
}

// ─── relation-factor input boundary ──────────────────────────────────────────

/// Relation pairs `(a, b)` validated against a [`SqrtContext`], which this value owns.
///
/// Each pair stands for the rational factor `a - b·m` and the algebraic factor `a - b·α`, i.e.
/// `(c·a - b·θ)/c` in `A`. Every factor is guaranteed nonzero and, in `A`, a unit. The pairs form
/// a multiset: duplicates are allowed, and `b` may be zero or negative. Binding the context here
/// means a relation set cannot be paired with a different polynomial later.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RelationFactors {
    context: SqrtContext,
    pairs: Vec<(BigInt, BigInt)>,
}

impl RelationFactors {
    /// Validate `pairs` against `ctx`, before any prime search or reconstruction.
    ///
    /// A factor `c·a - b·θ` has norm `c^(d-1) · f_hom(a, b)` with `f_hom(a, b) = Σ f_i a^i
    /// b^(d-i)`, so it is a unit of `A` exactly when `f_hom(a, b) ≠ 0`.
    ///
    /// # Errors
    ///
    /// The first failing pair, in order, determines the error:
    /// [`InvalidInput::EmptyRelationSet`]; [`InvalidInput::ZeroRationalFactor`] when
    /// `a - b·m = 0`; [`InvalidInput::ZeroAlgebraicFactor`] when `f_hom = 0` and the factor is the
    /// zero element (`d = 1`, or `a = b = 0`); [`UnsupportedDomain::ZeroDivisorFactor`] when
    /// `f_hom = 0` and the factor is a nonzero zero divisor.
    pub fn new(ctx: &SqrtContext, pairs: &[(BigInt, BigInt)]) -> Result<Self, InputRejection> {
        if pairs.is_empty() {
            return Err(invalid(InvalidInput::EmptyRelationSet));
        }
        for (index, (a, b)) in pairs.iter().enumerate() {
            check_factor(ctx, index, a, b)?;
        }
        Ok(Self { context: ctx.clone(), pairs: pairs.to_vec() })
    }

    /// The owning context.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        &self.context
    }

    /// The validated `(a, b)` pairs, in input order.
    #[must_use]
    pub fn pairs(&self) -> &[(BigInt, BigInt)] {
        &self.pairs
    }

    /// The number of relations `k >= 1`.
    #[must_use]
    pub fn relation_count(&self) -> usize {
        self.pairs.len()
    }
}

fn check_factor(
    ctx: &SqrtContext,
    index: usize,
    a: &BigInt,
    b: &BigInt,
) -> Result<(), InputRejection> {
    if (a - b * &ctx.m).is_zero() {
        return Err(invalid(InvalidInput::ZeroRationalFactor { index }));
    }
    if homogeneous_value(&ctx.f, ctx.degree, a, b).is_zero() {
        // With d = 1, A = Q and every nonunit is the zero element; for d >= 2 the element
        // c·a - b·t of degree < d is zero only when a = b = 0.
        let zero_element = ctx.degree == 1 || (a.is_zero() && b.is_zero());
        return Err(if zero_element {
            invalid(InvalidInput::ZeroAlgebraicFactor { index })
        } else {
            InputRejection::Unsupported(UnsupportedDomain::ZeroDivisorFactor { index })
        });
    }
    Ok(())
}

/// `Σ f_i a^i b^(d-i)`.
fn homogeneous_value(f: &IntPoly, degree: usize, a: &BigInt, b: &BigInt) -> BigInt {
    let a_powers = powers_of(a, degree);
    let b_powers = powers_of(b, degree);
    f.coeffs.iter().enumerate().map(|(i, f_i)| f_i * &a_powers[i] * &b_powers[degree - i]).sum()
}

// ─── resource options ────────────────────────────────────────────────────────

/// Budget for an ascending search over machine primes.
///
/// The search examines primes `p >= max(first_candidate, 3)` in order, and every prime examined
/// counts against `max_primes_examined`, including primes later skipped for any reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimeSearchLimits {
    /// Lower bound of the search; values below 3 are treated as 3.
    pub first_candidate: u64,
    /// Maximum number of primes examined.
    pub max_primes_examined: u64,
}

/// Budget for the size of the lifting modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrecisionLimits {
    /// Maximum bit length of the reconstruction modulus.
    pub max_modulus_bits: u64,
}

/// Budgets for direct reconstruction at an inert prime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DirectLimits {
    /// Search for a prime at which `F` is irreducible.
    pub prime_search: PrimeSearchLimits,
    /// Search for a nonresidue in the residue field.
    pub nonresidue_search: NonresidueSearchLimits,
}

/// Budgets for split-prime reconstruction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CrtLimits {
    /// Search for completely split primes.
    pub prime_search: PrimeSearchLimits,
    /// Upper bound on the number of primes used, at least 1.
    pub max_primes: u32,
    /// Exhaustive root search in each prime field.
    pub root_search: RootSearchLimits,
    /// Search for a nonresidue modulo each split prime.
    pub nonresidue_search: NonresidueSearchLimits,
}

/// Budgets for sign recombination.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecombinationLimits {
    /// Maximum entries in each meet-in-the-middle list.
    pub max_list_entries: u64,
    /// Maximum candidate matches examined.
    pub max_matches: u64,
}

/// Every budget of the square-root stage, passed by value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SqrtLimits {
    /// Precision budget.
    pub precision: PrecisionLimits,
    /// Direct reconstruction budgets.
    pub direct: DirectLimits,
    /// Split-prime reconstruction budgets.
    pub crt: CrtLimits,
    /// Sign recombination budgets.
    pub recombination: RecombinationLimits,
}

/// The resource whose budget ran out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResourceKind {
    /// Prime search for an irreducible reduction.
    InertPrimeSearch,
    /// Prime search for completely split primes.
    SplitPrimeSearch,
    /// The machine-prime range ended.
    MachinePrimeRange,
    /// Root enumeration in a prime field.
    RootSearch,
    /// Nonresidue search.
    NonresidueSearch,
    /// Lifting-modulus size.
    Precision,
    /// Sign-recombination list size.
    RecombinationList,
    /// Sign-recombination match count.
    RecombinationMatches,
}

/// A reached limit. Never evidence about whether a square root exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResourceExhaustion {
    /// Which budget was exhausted.
    pub kind: ResourceKind,
    /// Amount consumed (or required, for precision).
    pub consumed: u64,
    /// The configured limit.
    pub limit: u64,
}

// ─── owned outcomes and witnesses ────────────────────────────────────────────

/// How a verified root was reconstructed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReconstructionMethod {
    /// Inverse-square-root Newton lifting at one inert prime.
    DirectNewton {
        /// The inert prime.
        prime: u64,
        /// The lifting exponent `e` with modulus `prime^e`.
        exponent: u32,
    },
    /// Split-prime lifting with Chinese-remainder recombination.
    SplitPrimeCrt {
        /// The completely split primes used.
        primes: Vec<u64>,
        /// The common lifting exponent.
        exponent: u32,
    },
}

/// The local or exhaustive reason a target is certified not to be a square.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NonSquareEvidence {
    /// The target is a nonresidue in the residue field at an inert prime.
    InertLocalNonresidue {
        /// The inert prime.
        prime: u64,
    },
    /// Both signs of the unique candidate at an inert prime exceed the proved bound.
    DirectBoundExhausted {
        /// The inert prime.
        prime: u64,
        /// The lifting exponent.
        exponent: u32,
    },
    /// The target is a nonresidue at a root of a completely split prime.
    SplitLocalNonresidue {
        /// The split prime.
        prime: u64,
        /// The root in `F_prime`.
        root: u64,
    },
    /// Every sign combination was matched and rejected within the proved bound.
    SignSearchExhausted {
        /// The split primes used.
        primes: Vec<u64>,
        /// The common lifting exponent.
        exponent: u32,
    },
}

/// A certificate that the algebraic target is not a square.
///
/// It is minted only inside [`crate::sqrt`]; the guarantee is crate-external, not enforced by the
/// compiler against sibling modules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NonSquareCertificate {
    evidence: NonSquareEvidence,
}

impl NonSquareCertificate {
    pub(super) fn new(evidence: NonSquareEvidence) -> Self {
        Self { evidence }
    }

    /// The evidence behind this certificate.
    #[must_use]
    pub fn evidence(&self) -> &NonSquareEvidence {
        &self.evidence
    }
}

/// Plain bundle of the fields of a [`VerifiedRoot`]; private so only `seal` consumes it.
struct RootParts {
    context: SqrtContext,
    relation_count: usize,
    rational_product: BigInt,
    integral_root: IntPoly,
    beta: RatPoly,
    modular_image: BigInt,
    method: ReconstructionMethod,
}

/// An algebraic square root that has passed the exact verification.
///
/// With `A = Q[t]/(F)`, `θ = c·α`, and `γ = Π(a - b·α)`, a `VerifiedRoot` records `β` with
/// `β² = γ`, the integral root `z` with `z² ≡ H (mod F)`, the rational product
/// `P = Π(a - b·m)`, and the image `y = β(c·m) mod N`, which satisfies `y² ≡ P (mod N)`.
///
/// # Non-forgeability
///
/// Fields are private and the only constructor, `seal`, is private to the `sqrt::context` module
/// subtree. Code outside that subtree, including other crates and tests, cannot build a
/// `VerifiedRoot`:
///
/// ```compile_fail,E0451
/// use gnfs::sqrt::VerifiedRoot;
/// // Private fields: a struct literal is rejected.
/// let _forged = VerifiedRoot { relation_count: 1 };
/// ```
///
/// ```compile_fail,E0624
/// use gnfs::sqrt::VerifiedRoot;
/// // The sealing function is private.
/// let _ = VerifiedRoot::seal;
/// ```
///
/// Cloning an existing value is not forgery. There is no `Default` and no deserialization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedRoot {
    context: SqrtContext,
    relation_count: usize,
    rational_product: BigInt,
    integral_root: IntPoly,
    beta: RatPoly,
    modular_image: BigInt,
    method: ReconstructionMethod,
}

impl VerifiedRoot {
    fn seal(parts: RootParts) -> Self {
        let RootParts {
            context,
            relation_count,
            rational_product,
            integral_root,
            beta,
            modular_image,
            method,
        } = parts;
        Self {
            context,
            relation_count,
            rational_product,
            integral_root,
            beta,
            modular_image,
            method,
        }
    }

    /// The context the root was computed in.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        &self.context
    }

    /// The number of relations `k` in the product.
    #[must_use]
    pub fn relation_count(&self) -> usize {
        self.relation_count
    }

    /// `P = Π(a - b·m)` over `Z`, with the exact `m`.
    #[must_use]
    pub fn rational_product(&self) -> &BigInt {
        &self.rational_product
    }

    /// The integral root `z`, of degree `< d`, with `z² ≡ H (mod F)`.
    #[must_use]
    pub fn integral_root(&self) -> &IntPoly {
        &self.integral_root
    }

    /// The root `β` in the `θ` basis, of degree `< d`, with `β² = γ`.
    #[must_use]
    pub fn beta(&self) -> &RatPoly {
        &self.beta
    }

    /// The image `y = β(c·m) mod N`, in `[0, N)`.
    #[must_use]
    pub fn modular_image(&self) -> &BigInt {
        &self.modular_image
    }

    /// How the root was reconstructed.
    #[must_use]
    pub fn method(&self) -> &ReconstructionMethod {
        &self.method
    }
}

/// The mathematical or resource result of an algebraic square-root attempt.
///
/// Only `NotSquare` certifies that no root exists; `ResourceExhausted` never does.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlgebraicSqrtOutcome {
    /// A verified root.
    Verified(VerifiedRoot),
    /// The target is certified not to be a square.
    NotSquare(NonSquareCertificate),
    /// The domain or map is unsupported.
    Unsupported(UnsupportedDomain),
    /// A proper factor of `N` appeared incidentally.
    Incidental(IncidentalFactor),
    /// One or more budgets ran out, in chronological order; never empty.
    ResourceExhausted(Vec<ResourceExhaustion>),
}

impl AlgebraicSqrtOutcome {
    /// The failure class, or `None` for a verified root.
    #[must_use]
    pub fn class(&self) -> Option<FailureClass> {
        match self {
            Self::Verified(_) => None,
            Self::NotSquare(_) => Some(FailureClass::AlgebraicNonSquare),
            Self::Unsupported(_) => Some(FailureClass::Unsupported),
            Self::Incidental(_) => Some(FailureClass::IncidentalFactor),
            Self::ResourceExhausted(_) => Some(FailureClass::ResourceExhausted),
        }
    }
}

#[cfg(test)]
mod tests {
    use num_rational::BigRational;

    use super::*;

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    /// x² − 2 at m = 23 modulo 527 = 17·31: f(23) = 529 − 2 = 527.
    fn quadratic_context() -> SqrtContext {
        let f = IntPoly::from_coeffs(vec![big(-2), big(0), big(1)]);
        SqrtContext::new(&f, &big(23), &big(527)).unwrap()
    }

    #[test]
    fn seal_preserves_parts_and_accessors_expose_them() {
        // One relation (a, b) = (3, -2): γ = 3 + 2α = (1 + α)² in Q(√2), P = 3 + 2·23 = 49 = 7²;
        // β = 1 + α, y = β(23) = 24 and 24² = 576 = 527 + 49.
        let ctx = quadratic_context();
        let beta = RatPoly::from_coeffs(vec![BigRational::from(big(1)), BigRational::from(big(1))]);
        let z = IntPoly::from_coeffs(vec![big(4), big(2)]);
        let method = ReconstructionMethod::DirectNewton { prime: 3, exponent: 4 };
        let root = VerifiedRoot::seal(RootParts {
            context: ctx.clone(),
            relation_count: 1,
            rational_product: big(49),
            integral_root: z.clone(),
            beta: beta.clone(),
            modular_image: big(24),
            method: method.clone(),
        });
        assert_eq!(root.context(), &ctx);
        assert_eq!(root.relation_count(), 1);
        assert_eq!(root.rational_product(), &big(49));
        assert_eq!(root.integral_root(), &z);
        assert_eq!(root.beta(), &beta);
        assert_eq!(root.modular_image(), &big(24));
        assert_eq!(root.method(), &method);
        assert_eq!(root.clone(), root);

        let outcome = AlgebraicSqrtOutcome::Verified(root);
        assert_eq!(outcome.class(), None);
    }

    #[test]
    fn certificate_minted_in_module_reports_its_evidence_and_class() {
        let evidence = NonSquareEvidence::InertLocalNonresidue { prime: 3 };
        let cert = NonSquareCertificate::new(evidence.clone());
        assert_eq!(cert.evidence(), &evidence);
        let outcome = AlgebraicSqrtOutcome::NotSquare(cert);
        assert_eq!(outcome.class(), Some(FailureClass::AlgebraicNonSquare));
    }

    #[test]
    fn monic_form_scales_coefficients_by_powers_of_leading_coefficient() {
        // 2x³ + 3x² − x + 5 with c = 2: F = t³ + 3t² + (−1·2)t + 5·4 = t³ + 3t² − 2t + 20.
        let f = IntPoly::from_coeffs(vec![big(5), big(-1), big(3), big(2)]);
        let expected = IntPoly::from_coeffs(vec![big(20), big(-2), big(3), big(1)]);
        assert_eq!(monic_form(&f, &big(2), 3), expected);
    }
}
