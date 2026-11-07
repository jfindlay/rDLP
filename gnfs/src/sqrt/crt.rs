//! Split-prime local contributions to an algebraic square root.
//!
//! For a [`NormalizedTarget`] with monic `F` of degree `d` and target `H`, choose odd primes `p_i`
//! at which `F` splits completely into distinct linear factors and `H` is nonzero at every root
//! `r_ij`. Lift each root to `R_ij` modulo `Q_i = p_i^e` and take a local square root `s_ij` with
//! `s_ij² ≡ H(R_ij) (mod Q_i)`. Any square root `z` of `H` modulo `F` satisfies
//! `z(R_ij) ≡ ±s_ij (mod Q_i)`, because `H(R_ij)` is a unit and `Z/p^e` has exactly the unit
//! square roots `±s`. With the Lagrange basis `L_ij` of the roots modulo `Q_i`, `F ≡ Π(t − R_ij)`
//! there, and the Chinese-remainder weights `W_i` for `M = Π Q_i`, every root is therefore
//!
//! ```text
//! z ≡ Σ_ij ε_ij · W_i · s_ij · L_ij   (mod M),   ε_ij ∈ {−1, +1}
//! ```
//!
//! and, centered modulo `M > 2B`, equals its true coefficients exactly.
//!
//! This module builds the signed contributions `C_ij = W_i · s_ij · L_ij mod M` and checks every
//! local identity it relies on. It never chooses the signs `ε_ij` and never produces a verified
//! root: a local square at every root says nothing about a global square. A prime where `H` is a
//! nonresidue at a root is a genuine obstruction (`z(r)² ≡ H(r) (mod p)` has no solution), so it
//! yields a nonsquare certificate; exhausted budgets are resource outcomes only. See
//! `gnfs/docs/PEDAGOGY.md` for the derivation.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::pow::pow;
use num_traits::{One, Zero};
use shared_numfield::{
    FpPoly, IntPoly, NonresidueSearchLimits, PrimeModulus, ResidueField, RootEnumeration,
    SqrtOutcome, is_prime_u64,
};
use shared_padic::hensel_lift;

use super::{
    CrtLimits, InvariantViolation, NonSquareCertificate, NonSquareEvidence, NormalizedTarget,
    PrimeSearchLimits, ReconstructionMethod, ResourceExhaustion, ResourceKind, SqrtError,
};

/// The sign `ε` applied to one local contribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ContributionSign {
    /// Add the contribution.
    Plus,
    /// Subtract the contribution.
    Minus,
}

/// A selected split prime with its prime power and Chinese-remainder weight.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalPrime {
    prime: u64,
    modulus: BigInt,
    weight: BigInt,
}

impl LocalPrime {
    /// The prime `p_i`.
    #[must_use]
    pub fn prime(&self) -> u64 {
        self.prime
    }

    /// The prime power `Q_i = p_i^e`.
    #[must_use]
    pub fn modulus(&self) -> &BigInt {
        &self.modulus
    }

    /// The weight `W_i` in `[0, M)`: `≡ 1 (mod Q_i)` and `≡ 0 (mod Q_k)` for `k ≠ i`.
    #[must_use]
    pub fn weight(&self) -> &BigInt {
        &self.weight
    }
}

/// One root of `F` modulo a split prime: its lifts, Lagrange polynomial, and signed coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contribution {
    prime_index: usize,
    prime: u64,
    root: u64,
    lifted_root: BigInt,
    local_root: BigInt,
    lagrange: IntPoly,
    coefficients: Vec<BigInt>,
}

impl Contribution {
    /// Index into [`ContributionSet::locals`] of this contribution's prime.
    #[must_use]
    pub fn prime_index(&self) -> usize {
        self.prime_index
    }

    /// The prime `p`.
    #[must_use]
    pub fn prime(&self) -> u64 {
        self.prime
    }

    /// The root `r` in `[0, p)` with `F(r) ≡ 0 (mod p)`.
    #[must_use]
    pub fn root(&self) -> u64 {
        self.root
    }

    /// The lift `R` in `[0, Q)` with `R ≡ r (mod p)` and `F(R) ≡ 0 (mod Q)`.
    #[must_use]
    pub fn lifted_root(&self) -> &BigInt {
        &self.lifted_root
    }

    /// The local square root `s` in `[0, Q)` of `H(R) (mod Q)`; its sign is unspecified.
    #[must_use]
    pub fn local_root(&self) -> &BigInt {
        &self.local_root
    }

    /// The Lagrange basis polynomial `L` modulo `Q`: coefficients in `[0, Q)`, degree exactly
    /// `d − 1`, equal to 1 at this root and 0 at every other root of the same prime.
    #[must_use]
    pub fn lagrange(&self) -> &IntPoly {
        &self.lagrange
    }

    /// The coefficients of `C = W·s·L mod M`: exactly `d` entries in `[0, M)`.
    #[must_use]
    pub fn coefficients(&self) -> &[BigInt] {
        &self.coefficients
    }
}

/// The primes, exponent, and signed contributions for one target.
///
/// The set owns a clone of its target, so a set can never be paired with a different one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContributionSet {
    target: NormalizedTarget,
    locals: Vec<LocalPrime>,
    exponent: u32,
    modulus: BigInt,
    contributions: Vec<Contribution>,
}

impl ContributionSet {
    /// The target these contributions approximate a root of.
    #[must_use]
    pub fn target(&self) -> &NormalizedTarget {
        &self.target
    }

    /// The selected primes in ascending order; at least one.
    #[must_use]
    pub fn locals(&self) -> &[LocalPrime] {
        &self.locals
    }

    /// The common lifting exponent `e`.
    #[must_use]
    pub fn exponent(&self) -> u32 {
        self.exponent
    }

    /// The odd modulus `M = Π Q_i`, at least `2B + 1`.
    #[must_use]
    pub fn modulus(&self) -> &BigInt {
        &self.modulus
    }

    /// The `ℓ·d` contributions, ordered by ascending prime and then ascending root `r`.
    #[must_use]
    pub fn contributions(&self) -> &[Contribution] {
        &self.contributions
    }

    /// The reconstruction method to record on a root assembled from this set.
    #[must_use]
    pub fn method(&self) -> ReconstructionMethod {
        ReconstructionMethod::SplitPrimeCrt {
            primes: self.locals.iter().map(LocalPrime::prime).collect(),
            exponent: self.exponent,
        }
    }

    /// The centered residue of `Σ ε_k·C_k` modulo `M`, coefficientwise, as a polynomial of degree
    /// below `d` suitable for [`NormalizedTarget::verify_candidate`].
    ///
    /// This selects no signs and verifies nothing: the result is a candidate, and only the
    /// verifier can accept it.
    ///
    /// # Errors
    ///
    /// [`SqrtError::Invariant`] if `signs` does not have exactly one entry per contribution.
    pub fn assemble(&self, signs: &[ContributionSign]) -> Result<IntPoly, SqrtError> {
        if signs.len() != self.contributions.len() {
            return Err(violation(
                "sign vector length differs from contribution count",
            ));
        }
        let mut acc = vec![BigInt::zero(); self.target.context().degree()];
        for (contribution, sign) in self.contributions.iter().zip(signs) {
            for (slot, c) in acc.iter_mut().zip(&contribution.coefficients) {
                match sign {
                    ContributionSign::Plus => *slot += c,
                    ContributionSign::Minus => *slot -= c,
                }
            }
        }
        let reduced = acc
            .iter()
            .map(|a| centered(&a.mod_floor(&self.modulus), &self.modulus));
        Ok(IntPoly::from_coeffs(reduced.collect()))
    }
}

/// The result of selecting split primes and building contributions.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitOutcome {
    /// Contributions for at least one prime; no sign has been chosen.
    Ready(ContributionSet),
    /// `H` is a nonresidue at a root of a split prime: [`NonSquareEvidence::SplitLocalNonresidue`].
    NotSquare(NonSquareCertificate),
    /// A budget ran out before any usable prime was found. Never evidence about squareness.
    ResourceExhausted(ResourceExhaustion),
}

/// Select completely split primes for `target` and build their signed contributions.
///
/// Primes `p ≥ max(first_candidate, 3)` are examined in ascending order and every prime examined
/// counts against `prime_search.max_primes_examined`. A prime is usable when it is odd, does not
/// divide `c·N`, `F` is squarefree modulo `p` with `d` distinct roots (decided from
/// `is_squarefree` and `distinct_root_count`, never from a list), the complete root list fits
/// `root_search`, and `H(r) ≢ 0` at every root; otherwise it is skipped. Roots are taken in
/// ascending order; the first `H(r)` with no square root modulo `p` returns
/// [`SplitOutcome::NotSquare`]. The search stops at `max_primes` usable primes; when the budget or
/// the machine-prime range ends first, the primes found so far are used (so `ℓ` may be smaller)
/// and no exhaustion is reported unless none was found.
///
/// The exponent `e` is the least value with `(Π p_i)^e ≥ 2B + 1`, so `M = Π Q_i` exceeds `2B`.
///
/// Outcomes other than `Ready` and `NotSquare` are [`SplitOutcome::ResourceExhausted`] with
/// [`ResourceKind::SplitPrimeSearch`] (also for `max_primes = 0`), [`ResourceKind::MachinePrimeRange`],
/// [`ResourceKind::RootSearch`], [`ResourceKind::NonresidueSearch`], or [`ResourceKind::Precision`]
/// (the exponent would not fit a `u32`).
///
/// # Errors
///
/// [`SqrtError::Invariant`] if a lifting identity, a Lagrange or Chinese-remainder identity, or a
/// substrate precondition fails. These indicate a defect, never a nonsquare.
pub fn split_contributions(
    target: &NormalizedTarget,
    limits: CrtLimits,
) -> Result<SplitOutcome, SqrtError> {
    if limits.max_primes == 0 {
        return Ok(SplitOutcome::ResourceExhausted(ResourceExhaustion {
            kind: ResourceKind::SplitPrimeSearch,
            consumed: 0,
            limit: 0,
        }));
    }
    let selected = match select_primes(target, limits)? {
        Selection::Primes(selected) => selected,
        Selection::NotSquare(certificate) => return Ok(SplitOutcome::NotSquare(certificate)),
        Selection::Exhausted(exhausted) => return Ok(SplitOutcome::ResourceExhausted(exhausted)),
    };

    let product: BigInt = selected.iter().map(|s| BigInt::from(s.prime)).product();
    let minimum = target.minimum_modulus();
    let Some(exponent) = least_exponent(&product, minimum, u32::MAX) else {
        return Ok(SplitOutcome::ResourceExhausted(ResourceExhaustion {
            kind: ResourceKind::Precision,
            consumed: minimum.bits(),
            limit: u64::from(u32::MAX),
        }));
    };

    let lifted = selected
        .iter()
        .map(|prime| lift_prime(target, prime, exponent))
        .collect::<Result<Vec<_>, _>>()?;
    let moduli: Vec<BigInt> = lifted.iter().map(|l| l.modulus.clone()).collect();
    let (modulus, weights) = crt_weights(&moduli)?;
    Ok(SplitOutcome::Ready(assemble_set(
        target, lifted, weights, exponent, modulus,
    )))
}

// ─── prime search ────────────────────────────────────────────────────────────

/// Ascending iterator over machine primes that counts every prime it hands out.
struct PrimeCursor {
    next: Option<u64>,
    examined: u64,
    budget: u64,
}

impl PrimeCursor {
    fn new(limits: PrimeSearchLimits) -> Self {
        Self {
            next: Some(limits.first_candidate.max(3)),
            examined: 0,
            budget: limits.max_primes_examined,
        }
    }

    /// The next prime, or the exhaustion that stopped the search.
    fn advance(&mut self) -> Result<u64, ResourceExhaustion> {
        if self.examined >= self.budget {
            return Err(self.exhaustion(ResourceKind::SplitPrimeSearch));
        }
        let mut candidate = self.next;
        while let Some(value) = candidate {
            if is_prime_u64(value) {
                self.examined += 1;
                self.next = value.checked_add(1);
                return Ok(value);
            }
            candidate = value.checked_add(1);
        }
        self.next = None;
        Err(self.exhaustion(ResourceKind::MachinePrimeRange))
    }

    fn exhaustion(&self, kind: ResourceKind) -> ResourceExhaustion {
        ResourceExhaustion {
            kind,
            consumed: self.examined,
            limit: self.budget,
        }
    }
}

/// A usable split prime: its roots of `F` and a residue square root of `H` at each.
struct SplitPrime {
    prime: u64,
    roots: Vec<u64>,
    residue_roots: Vec<u64>,
}

/// The classification of one examined prime.
enum Screen {
    Skip,
    Usable(SplitPrime),
    NotSquare(NonSquareCertificate),
    Exhausted(ResourceExhaustion),
}

/// What the prime search concluded.
enum Selection {
    Primes(Vec<SplitPrime>),
    NotSquare(NonSquareCertificate),
    Exhausted(ResourceExhaustion),
}

/// Examine primes in ascending order until `max_primes` are usable or the search ends.
fn select_primes(target: &NormalizedTarget, limits: CrtLimits) -> Result<Selection, SqrtError> {
    let wanted = usize::try_from(limits.max_primes).unwrap_or(usize::MAX);
    let mut cursor = PrimeCursor::new(limits.prime_search);
    let mut found: Vec<SplitPrime> = Vec::new();
    while found.len() < wanted {
        let prime = match cursor.advance() {
            Ok(prime) => prime,
            Err(exhausted) if found.is_empty() => return Ok(Selection::Exhausted(exhausted)),
            // At least one usable prime exists: proceed with fewer than `max_primes`.
            Err(_) => break,
        };
        match screen_prime(target, prime, limits)? {
            Screen::Skip => {}
            Screen::Usable(split) => found.push(split),
            Screen::NotSquare(certificate) => return Ok(Selection::NotSquare(certificate)),
            Screen::Exhausted(exhausted) => return Ok(Selection::Exhausted(exhausted)),
        }
    }
    Ok(Selection::Primes(found))
}

/// Classify the prime `p`: split test, root enumeration, zero test of `H`, then residue roots.
fn screen_prime(target: &NormalizedTarget, p: u64, limits: CrtLimits) -> Result<Screen, SqrtError> {
    if p.is_multiple_of(2) {
        return Ok(Screen::Skip);
    }
    let modulus = PrimeModulus::new(p).map_err(|_| violation("candidate prime is not prime"))?;
    let context = target.context();
    let p_big = BigInt::from(p);
    if context.leading_coeff().is_multiple_of(&p_big) || context.modulus().is_multiple_of(&p_big) {
        return Ok(Screen::Skip);
    }
    let degree = context.degree();
    let monic = FpPoly::from_int_poly(context.monic_poly(), modulus);
    if monic.degree() != Some(degree) {
        return Err(violation(
            "monic polynomial lost its degree modulo an odd prime",
        ));
    }
    let count = monic
        .distinct_root_count()
        .map_err(|_| violation("root count of nonzero F failed"))?;
    if !monic.is_squarefree() || count != degree {
        return Ok(Screen::Skip);
    }
    let roots = match monic
        .roots(limits.root_search)
        .map_err(|_| violation("root enumeration of nonzero F failed"))?
    {
        RootEnumeration::Complete { roots, .. } => roots,
        RootEnumeration::Truncated { evaluations, .. } => {
            return Ok(Screen::Exhausted(ResourceExhaustion {
                kind: ResourceKind::RootSearch,
                consumed: evaluations,
                limit: limits.root_search.max_evaluations,
            }));
        }
    };
    if roots.len() != degree {
        return Err(violation(
            "complete root list disagrees with the exact root count",
        ));
    }

    let residue = FpPoly::from_int_poly(target.target(), modulus);
    let values: Vec<u64> = roots.iter().map(|&r| residue.eval(r)).collect();
    if values.contains(&0) {
        return Ok(Screen::Skip);
    }
    residue_roots(modulus, &roots, &values, limits.nonresidue_search)
}

/// Square roots of the values `H(r)` in `F_p`, using the residue field `F_p[t]/(t)`.
fn residue_roots(
    modulus: PrimeModulus,
    roots: &[u64],
    values: &[u64],
    nonresidue: NonresidueSearchLimits,
) -> Result<Screen, SqrtError> {
    let linear = FpPoly::monomial(modulus, 1, 1);
    let field =
        ResidueField::new(&linear).map_err(|_| violation("linear modulus rejected as a field"))?;
    let mut residue_roots = Vec::with_capacity(values.len());
    for (&root, &value) in roots.iter().zip(values) {
        let element = FpPoly::constant(modulus, value);
        match field
            .sqrt(&element, nonresidue)
            .map_err(|_| violation("residue square root failed"))?
        {
            SqrtOutcome::Root(s) => {
                let s0 = s.coeffs().first().copied().unwrap_or(0);
                if modulus.mul(s0, s0) != value {
                    return Err(violation("residue square root does not square to H(r)"));
                }
                residue_roots.push(s0);
            }
            SqrtOutcome::NotSquare => {
                return Ok(Screen::NotSquare(NonSquareCertificate::new(
                    NonSquareEvidence::SplitLocalNonresidue {
                        prime: modulus.get(),
                        root,
                    },
                )));
            }
            SqrtOutcome::ResourceExhausted { candidates_tried } => {
                return Ok(Screen::Exhausted(ResourceExhaustion {
                    kind: ResourceKind::NonresidueSearch,
                    consumed: candidates_tried,
                    limit: nonresidue.max_candidates,
                }));
            }
        }
    }
    Ok(Screen::Usable(SplitPrime {
        prime: modulus.get(),
        roots: roots.to_vec(),
        residue_roots,
    }))
}

/// The least `e ≥ 1` with `product^e ≥ minimum`, or `None` if it would exceed `cap`.
fn least_exponent(product: &BigInt, minimum: &BigInt, cap: u32) -> Option<u32> {
    let mut power = product.clone();
    let mut exponent: u32 = 1;
    while power < *minimum {
        if exponent >= cap {
            return None;
        }
        power *= product;
        exponent += 1;
    }
    Some(exponent)
}

// ─── local lifting ───────────────────────────────────────────────────────────

/// The lifted data of one root of `F` modulo `Q = p^e`.
struct LiftedRoot {
    root: u64,
    lifted: BigInt,
    local: BigInt,
}

/// The lifted roots and Lagrange basis of one prime.
struct LiftedPrime {
    prime: u64,
    modulus: BigInt,
    roots: Vec<LiftedRoot>,
    lagrange: Vec<IntPoly>,
}

/// Lift every root of one prime to `Q = p^e`, check the local identities, and interpolate.
fn lift_prime(
    target: &NormalizedTarget,
    split: &SplitPrime,
    exponent: u32,
) -> Result<LiftedPrime, SqrtError> {
    let p = BigInt::from(split.prime);
    let modulus = pow(p.clone(), exponent as usize);
    let monic = target.context().monic_poly();
    let mut roots = Vec::with_capacity(split.roots.len());
    for (&root, &s0) in split.roots.iter().zip(&split.residue_roots) {
        let lifted = hensel_lift(monic, &BigInt::from(root), &p, exponent)
            .map_err(|_| violation("Hensel lift of a root of F failed"))?;
        check_lifted_root(monic, &lifted, root, &p, &modulus)?;
        let value = target.target().eval_mod(&lifted, &modulus);
        let square = IntPoly::from_coeffs(vec![-&value, BigInt::zero(), BigInt::one()]);
        let local = hensel_lift(&square, &BigInt::from(s0), &p, exponent)
            .map_err(|_| violation("Hensel lift of a square root of H(R) failed"))?;
        check_local_root(&local, &value, s0, &p, &modulus)?;
        roots.push(LiftedRoot {
            root,
            lifted,
            local,
        });
    }
    let points: Vec<BigInt> = roots.iter().map(|r| r.lifted.clone()).collect();
    let lagrange = lagrange_basis(&points, &modulus)?;
    check_lagrange(&lagrange, &points, &modulus)?;
    Ok(LiftedPrime {
        prime: split.prime,
        modulus,
        roots,
        lagrange,
    })
}

/// `R ∈ [0, Q)`, `R ≡ r (mod p)`, and `F(R) ≡ 0 (mod Q)`.
fn check_lifted_root(
    monic: &IntPoly,
    lifted: &BigInt,
    root: u64,
    p: &BigInt,
    modulus: &BigInt,
) -> Result<(), SqrtError> {
    if out_of_range(lifted, modulus) || (lifted - BigInt::from(root)).mod_floor(p) != BigInt::zero()
    {
        return Err(violation(
            "lifted root is out of range or does not reduce to its root",
        ));
    }
    if !monic.eval_mod(lifted, modulus).is_zero() {
        return Err(violation("lifted root does not satisfy F(R) ≡ 0"));
    }
    Ok(())
}

/// `s ∈ [0, Q)`, `s ≡ s0 (mod p)`, and `s² ≡ H(R) (mod Q)`.
fn check_local_root(
    local: &BigInt,
    value: &BigInt,
    s0: u64,
    p: &BigInt,
    modulus: &BigInt,
) -> Result<(), SqrtError> {
    if out_of_range(local, modulus) || (local - BigInt::from(s0)).mod_floor(p) != BigInt::zero() {
        return Err(violation(
            "local root is out of range or does not reduce to its seed",
        ));
    }
    if (local * local - value).mod_floor(modulus) != BigInt::zero() {
        return Err(violation("local root does not satisfy s² ≡ H(R)"));
    }
    Ok(())
}

/// True unless `value ∈ [0, modulus)`.
fn out_of_range(value: &BigInt, modulus: &BigInt) -> bool {
    *value < BigInt::zero() || value >= modulus
}

/// The inverse of `a` modulo `q`, in `[0, q)`; fails when `a` is not a unit.
fn mod_inverse(a: &BigInt, q: &BigInt) -> Result<BigInt, SqrtError> {
    let reduced = a.mod_floor(q);
    let gcd = reduced.extended_gcd(q);
    if !gcd.gcd.is_one() {
        return Err(violation("a value that must be a unit has no inverse"));
    }
    let inverse = gcd.x.mod_floor(q);
    if (&reduced * &inverse).mod_floor(q) != BigInt::one().mod_floor(q) {
        return Err(violation("modular inverse does not multiply to one"));
    }
    Ok(inverse)
}

/// Lagrange basis `L_j = Π_{k≠j}(t − R_k) · (Π_{k≠j}(R_j − R_k))⁻¹` modulo `q`.
fn lagrange_basis(points: &[BigInt], modulus: &BigInt) -> Result<Vec<IntPoly>, SqrtError> {
    let mut basis = Vec::with_capacity(points.len());
    for (j, point) in points.iter().enumerate() {
        let mut numerator = vec![BigInt::one()];
        let mut denominator = BigInt::one();
        for other in points
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != j)
            .map(|(_, o)| o)
        {
            numerator = times_linear(&numerator, other, modulus);
            denominator = (denominator * (point - other)).mod_floor(modulus);
        }
        let inverse = mod_inverse(&denominator, modulus)?;
        let scaled = numerator
            .iter()
            .map(|c| (c * &inverse).mod_floor(modulus))
            .collect();
        basis.push(IntPoly::from_coeffs(scaled));
    }
    Ok(basis)
}

/// `poly · (t − point)` modulo `modulus`, coefficients in `[0, modulus)`.
fn times_linear(poly: &[BigInt], point: &BigInt, modulus: &BigInt) -> Vec<BigInt> {
    let mut product = vec![BigInt::zero(); poly.len() + 1];
    for (i, c) in poly.iter().enumerate() {
        product[i + 1] += c;
        product[i] -= c * point;
    }
    product.iter().map(|c| c.mod_floor(modulus)).collect()
}

/// `L_j(R_k) ≡ δ_jk`, `Σ L_j ≡ 1`, and `deg L_j = d − 1`, all modulo `q`.
fn check_lagrange(basis: &[IntPoly], points: &[BigInt], modulus: &BigInt) -> Result<(), SqrtError> {
    let degree = points.len();
    let mut total = vec![BigInt::zero(); degree];
    for (j, poly) in basis.iter().enumerate() {
        if poly.coeffs.len() != degree {
            return Err(violation("Lagrange polynomial does not have degree d − 1"));
        }
        for (k, point) in points.iter().enumerate() {
            let expected = if j == k {
                BigInt::one()
            } else {
                BigInt::zero()
            };
            if poly.eval_mod(point, modulus) != expected.mod_floor(modulus) {
                return Err(violation(
                    "Lagrange polynomial is not 1 at its root and 0 at the others",
                ));
            }
        }
        for (slot, c) in total.iter_mut().zip(&poly.coeffs) {
            *slot += c;
        }
    }
    let one_then_zeros = total.iter().enumerate().all(|(i, c)| {
        let expected = if i == 0 {
            BigInt::one()
        } else {
            BigInt::zero()
        };
        (c - expected).mod_floor(modulus).is_zero()
    });
    if one_then_zeros {
        Ok(())
    } else {
        Err(violation("Lagrange polynomials do not sum to 1"))
    }
}

// ─── Chinese remainder weights ───────────────────────────────────────────────

/// `M = Π Q_i` and the weights `W_i = (M/Q_i)·((M/Q_i)⁻¹ mod Q_i) mod M`.
fn crt_weights(moduli: &[BigInt]) -> Result<(BigInt, Vec<BigInt>), SqrtError> {
    let modulus: BigInt = moduli.iter().product();
    let mut weights = Vec::with_capacity(moduli.len());
    for q in moduli {
        let cofactor = &modulus / q;
        let inverse = mod_inverse(&cofactor, q)?;
        weights.push((cofactor * inverse).mod_floor(&modulus));
    }
    check_weights(&weights, moduli, &modulus)?;
    Ok((modulus, weights))
}

/// `W_i ≡ δ_ik (mod Q_k)`, `W_i² ≡ W_i`, and `Σ W_i ≡ 1 (mod M)`.
fn check_weights(weights: &[BigInt], moduli: &[BigInt], modulus: &BigInt) -> Result<(), SqrtError> {
    for (i, weight) in weights.iter().enumerate() {
        for (k, q) in moduli.iter().enumerate() {
            let expected = if i == k {
                BigInt::one()
            } else {
                BigInt::zero()
            };
            if (weight - expected).mod_floor(q) != BigInt::zero() {
                return Err(violation(
                    "CRT weight is not an indicator of its prime power",
                ));
            }
        }
        if (weight * weight - weight).mod_floor(modulus) != BigInt::zero() {
            return Err(violation("CRT weight is not idempotent"));
        }
    }
    let sum: BigInt = weights.iter().sum();
    if (sum - BigInt::one()).mod_floor(modulus) != BigInt::zero() {
        return Err(violation("CRT weights do not sum to 1"));
    }
    Ok(())
}

// ─── assembly of the set ─────────────────────────────────────────────────────

/// Combine lifted primes and weights into `C_ij = W_i·s_ij·L_ij mod M`, padded to `d` entries.
fn assemble_set(
    target: &NormalizedTarget,
    lifted: Vec<LiftedPrime>,
    weights: Vec<BigInt>,
    exponent: u32,
    modulus: BigInt,
) -> ContributionSet {
    let degree = target.context().degree();
    let mut locals = Vec::with_capacity(lifted.len());
    let mut contributions = Vec::with_capacity(lifted.len() * degree);
    for (index, (prime, weight)) in lifted.into_iter().zip(weights).enumerate() {
        for (root, lagrange) in prime.roots.into_iter().zip(prime.lagrange) {
            let scale = (&weight * &root.local).mod_floor(&modulus);
            let coefficients = (0..degree)
                .map(|k| {
                    let l_k = lagrange.coeffs.get(k).cloned().unwrap_or_else(BigInt::zero);
                    (&scale * l_k).mod_floor(&modulus)
                })
                .collect();
            contributions.push(Contribution {
                prime_index: index,
                prime: prime.prime,
                root: root.root,
                lifted_root: root.lifted,
                local_root: root.local,
                lagrange,
                coefficients,
            });
        }
        locals.push(LocalPrime {
            prime: prime.prime,
            modulus: prime.modulus,
            weight,
        });
    }
    ContributionSet {
        target: target.clone(),
        locals,
        exponent,
        modulus,
        contributions,
    }
}

/// The residue of `c` in the symmetric range around zero; `modulus` is odd.
fn centered(c: &BigInt, modulus: &BigInt) -> BigInt {
    if c * 2 > *modulus {
        c - modulus
    } else {
        c.clone()
    }
}

fn violation(what: &'static str) -> SqrtError {
    SqrtError::Invariant(InvariantViolation { what })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn poly(coeffs: &[i64]) -> IntPoly {
        IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
    }

    fn is_invariant<T: std::fmt::Debug>(result: Result<T, SqrtError>) -> bool {
        matches!(result, Err(SqrtError::Invariant(_)))
    }

    #[test]
    fn prime_cursor_counts_only_primes_and_reports_split_search() {
        let mut cursor = PrimeCursor::new(PrimeSearchLimits {
            first_candidate: 0,
            max_primes_examined: 2,
        });
        assert_eq!(cursor.advance(), Ok(3));
        assert_eq!(cursor.advance(), Ok(5));
        let exhausted = cursor.advance().unwrap_err();
        assert_eq!(exhausted.kind, ResourceKind::SplitPrimeSearch);
        assert_eq!((exhausted.consumed, exhausted.limit), (2, 2));
    }

    #[test]
    fn prime_cursor_reports_the_end_of_the_machine_range() {
        // 2^64 − 59 is the largest prime below 2^64.
        let largest = u64::MAX - 58;
        let mut cursor = PrimeCursor::new(PrimeSearchLimits {
            first_candidate: largest,
            max_primes_examined: 5,
        });
        assert_eq!(cursor.advance(), Ok(largest));
        let end = cursor.advance().unwrap_err();
        assert_eq!(end.kind, ResourceKind::MachinePrimeRange);
        assert_eq!((end.consumed, end.limit), (1, 5));
        assert_eq!(
            cursor.advance().unwrap_err().kind,
            ResourceKind::MachinePrimeRange
        );
    }

    #[test]
    fn least_exponent_is_least_and_respects_its_cap() {
        // Powers of 7: 7, 49, 343. The minimum is 2B + 1, so equality is allowed.
        assert_eq!(least_exponent(&big(7), &big(7), u32::MAX), Some(1));
        assert_eq!(least_exponent(&big(7), &big(8), u32::MAX), Some(2));
        assert_eq!(least_exponent(&big(7), &big(343), u32::MAX), Some(3));
        assert_eq!(least_exponent(&big(7), &big(344), 3), None);
        assert_eq!(least_exponent(&big(7), &big(344), 4), Some(4));
    }

    #[test]
    fn centered_residue_is_symmetric_for_odd_modulus() {
        // M = 27: residues 0…13 stay, 14…26 map to −13…−1.
        assert_eq!(centered(&big(13), &big(27)), big(13));
        assert_eq!(centered(&big(14), &big(27)), big(-13));
        assert_eq!(centered(&big(26), &big(27)), big(-1));
    }

    #[test]
    fn modular_inverse_rejects_nonunits() {
        // 4·13 = 52 = 3·17 + 1, and −4 ≡ 13. 0 is not a unit, and 34 and 51 share the factor 17.
        assert_eq!(mod_inverse(&big(4), &big(17)).unwrap(), big(13));
        assert_eq!(mod_inverse(&big(-4), &big(17)).unwrap(), big(4));
        assert!(is_invariant(mod_inverse(&big(0), &big(17))));
        assert!(is_invariant(mod_inverse(&big(34), &big(51))));
    }

    #[test]
    fn identity_checks_reject_corrupted_data() {
        // F = t² − 2 modulo 49 has the root 10 over 3 (10² = 100 = 2·49 + 2); 39 = −10.
        let monic = poly(&[-2, 0, 1]);
        let p = big(7);
        let q = big(49);
        assert!(check_lifted_root(&monic, &big(10), 3, &p, &q).is_ok());
        assert!(is_invariant(check_lifted_root(&monic, &big(11), 3, &p, &q)));
        assert!(is_invariant(check_lifted_root(&monic, &big(59), 3, &p, &q)));
        assert!(is_invariant(check_lifted_root(&monic, &big(10), 4, &p, &q)));

        // 24² = 576 = 11·49 + 37, so 24 is a square root of 37 over the seed 3.
        assert!(check_local_root(&big(24), &big(37), 3, &p, &q).is_ok());
        assert!(is_invariant(check_local_root(
            &big(24),
            &big(38),
            3,
            &p,
            &q
        )));
        assert!(is_invariant(check_local_root(
            &big(24),
            &big(37),
            4,
            &p,
            &q
        )));
        assert!(is_invariant(check_local_root(
            &big(-25),
            &big(37),
            3,
            &p,
            &q
        )));

        // L_0 = 27t + 25 and L_1 = 22t + 25 interpolate the roots 10 and 39 modulo 49.
        let basis = [poly(&[25, 27]), poly(&[25, 22])];
        let points = [big(10), big(39)];
        assert!(check_lagrange(&basis, &points, &q).is_ok());
        assert!(is_invariant(check_lagrange(
            &[poly(&[25, 27]), poly(&[25, 21])],
            &points,
            &q
        )));
        assert!(is_invariant(check_lagrange(
            &[poly(&[1]), poly(&[0, 1])],
            &points,
            &q
        )));
        // L_0 alone is 1 at 10 and 0 at 39 but does not sum to 1 without L_1.
        assert!(is_invariant(check_lagrange(&basis[..1], &points, &q)));
    }

    #[test]
    fn weight_checks_reject_corrupted_data() {
        // Moduli 7 and 23: W = (92, 70) is 1 mod its own prime power and 0 mod the other.
        let moduli = [big(7), big(23)];
        let modulus = big(161);
        assert!(check_weights(&[big(92), big(70)], &moduli, &modulus).is_ok());
        assert!(is_invariant(check_weights(
            &[big(70), big(92)],
            &moduli,
            &modulus
        )));
        assert!(is_invariant(check_weights(
            &[big(92), big(0)],
            &moduli,
            &modulus
        )));
        let (m, w) = crt_weights(&moduli).unwrap();
        assert_eq!((m, w), (big(161), vec![big(92), big(70)]));
    }
}
