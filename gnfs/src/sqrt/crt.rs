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
//! [`split_contributions`] builds the signed contributions `C_ij = W_i · s_ij · L_ij mod M` and
//! checks every local identity it relies on. It never chooses the signs `ε_ij`: a local square at
//! every root says nothing about a global square. A prime where `H` is a nonresidue at a root is a
//! genuine obstruction (`z(r)² ≡ H(r) (mod p)` has no solution), so it yields a nonsquare
//! certificate; exhausted budgets are resource outcomes only.
//!
//! [`recombine_signs`] then finds the signs by meet-in-the-middle matching on the coefficient with
//! the smallest bound, and [`split_prime_root`] runs both steps. The first contribution is fixed
//! to `+`, which removes the `z ↔ −z` symmetry and nothing else. A sign vector is accepted only
//! through [`NormalizedTarget::verify_candidate`]; a search cut short by a budget is a resource
//! outcome, and only a complete search whose every match was rejected certifies a nonsquare. See
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
    AlgebraicSqrtOutcome, CandidateRejection, CandidateVerdict, CrtLimits, InvariantViolation,
    NonSquareCertificate, NonSquareEvidence, NormalizedTarget, PrimeSearchLimits,
    RecombinationLimits, ReconstructionMethod, ResourceExhaustion, ResourceKind, SqrtError,
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

/// Reconstruct a verified square root of `target` by split-prime lifting and sign recombination.
///
/// Runs [`split_contributions`] and, when it yields contributions, [`recombine_signs`]. A local
/// nonresidue becomes a nonsquare certificate; a budget that ends the prime search is a
/// single-entry [`AlgebraicSqrtOutcome::ResourceExhausted`].
///
/// # Errors
///
/// [`SqrtError::Invariant`] as for [`split_contributions`] and [`recombine_signs`].
pub fn split_prime_root(
    target: &NormalizedTarget,
    crt: CrtLimits,
    recombination: RecombinationLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    match split_contributions(target, crt)? {
        SplitOutcome::Ready(set) => recombine_signs(&set, recombination),
        SplitOutcome::NotSquare(certificate) => Ok(AlgebraicSqrtOutcome::NotSquare(certificate)),
        SplitOutcome::ResourceExhausted(exhausted) => {
            Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![exhausted]))
        }
    }
}

/// Find the signs `ε` for which `Σ ε_k·C_k` is a verified root, by meet-in-the-middle matching.
///
/// Let `n = ℓ·d` and `c_k = C_k[κ]`, where `κ` is the coefficient with the smallest bound `b = B_κ`
/// (ties go to the lowest index). The first contribution is fixed to `+`, which removes the
/// symmetry `z ↔ −z`; no other sign is normalised. The remaining indices split into a left half
/// `1..=⌊(n−1)/2⌋` and a right half. The left list holds `c_0 + Σ ±c_k` over the left signs and the
/// right list holds `Σ ±c_k` over the right signs, both modulo `M`, enumerated by sign mask (bit
/// `j` set means `−`). A pair is a match when `x + y` has centered residue in `[−b, b]`, that is,
/// `y` lies in the residue interval `[−x−b, −x+b] (mod M)`. An interval that wraps past `M` is
/// scanned as the two exact ranges `[lo, M−1]` and `[0, hi]`; since `2b + 1 ≤ M` they never
/// overlap. The right list is sorted by `(residue, mask)` and searched by binary search, so
/// matches are examined in ascending left mask, then interval order, then right residue.
///
/// Each match is assembled in full and passed to [`NormalizedTarget::verify_candidate`], which
/// checks every coefficient bound and `z² = H` exactly. A close match on `κ` alone proves nothing.
/// [`AlgebraicSqrtOutcome::Verified`], [`AlgebraicSqrtOutcome::Incidental`] and
/// [`AlgebraicSqrtOutcome::Unsupported`] are returned at once: a root exists. Rejected matches are
/// skipped.
///
/// Budgets: if the larger list would exceed `max_list_entries`, the result is
/// [`ResourceKind::RecombinationList`] with the required entry count (saturated at `u64::MAX`),
/// before any list is built. Before examining the match numbered `max_matches + 1`, the result is
/// [`ResourceKind::RecombinationMatches`]. Neither is evidence about squareness. Only a complete
/// search in which every match was rejected, including a search with no match at all, returns
/// [`NonSquareEvidence::SignSearchExhausted`]: some `±z` has `ε_0 = +` and `|z_κ| ≤ b`, so a root
/// would always have been matched.
///
/// # Errors
///
/// [`SqrtError::Invariant`] if the set is empty or malformed, `2b + 1 > M`, or the verifier
/// reports a candidate of degree `≥ d`. These indicate a defect, never a nonsquare.
pub fn recombine_signs(
    set: &ContributionSet,
    limits: RecombinationLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    let contributions = set.contributions.len();
    let bounds = set.target.coefficient_bounds();
    let Some(kappa) = (0..bounds.len()).min_by(|&a, &b| bounds[a].cmp(&bounds[b]).then(a.cmp(&b)))
    else {
        return Err(violation("coefficient bounds are empty"));
    };
    let bound = &bounds[kappa];
    if contributions == 0 || bound * 2 + 1 > set.modulus {
        return Err(violation(
            "contribution set is empty or its modulus does not exceed twice the bound",
        ));
    }
    let column = set
        .contributions
        .iter()
        .map(|c| c.coefficients.get(kappa))
        .collect::<Option<Vec<&BigInt>>>()
        .ok_or_else(|| violation("contribution has fewer coefficients than the degree"))?;

    let halves = Halves::of(contributions);
    let lists = match MeetInTheMiddle::build(&column, halves, bound, &set.modulus, limits) {
        Ok(lists) => lists,
        Err(exhausted) => return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![exhausted])),
    };

    let mut examined: u64 = 0;
    for (left_mask, right_mask) in lists.candidates() {
        if examined >= limits.max_matches {
            return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![
                ResourceExhaustion {
                    kind: ResourceKind::RecombinationMatches,
                    consumed: examined,
                    limit: limits.max_matches,
                },
            ]));
        }
        examined += 1;
        let signs = halves.signs(left_mask, right_mask);
        let candidate = set.assemble(&signs)?;
        match set.target.verify_candidate(&candidate, set.method())? {
            CandidateVerdict::Verified(root) => return Ok(AlgebraicSqrtOutcome::Verified(root)),
            CandidateVerdict::Incidental(factor) => {
                return Ok(AlgebraicSqrtOutcome::Incidental(factor));
            }
            CandidateVerdict::Unsupported(domain) => {
                return Ok(AlgebraicSqrtOutcome::Unsupported(domain));
            }
            CandidateVerdict::Rejected(CandidateRejection::DegreeTooHigh) => {
                return Err(violation("assembled candidate has degree at least d"));
            }
            CandidateVerdict::Rejected(_) => {}
        }
    }
    Ok(AlgebraicSqrtOutcome::NotSquare(NonSquareCertificate::new(
        NonSquareEvidence::SignSearchExhausted {
            primes: set.locals.iter().map(LocalPrime::prime).collect(),
            exponent: set.exponent,
        },
    )))
}

/// The split of the `n − 1` free signs into a left half and a right half.
///
/// Contribution `0` is fixed to `+`; contributions `1..=left` form the left half and the rest the
/// right half, so the right half is never the smaller one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Halves {
    total: usize,
    left: usize,
}

impl Halves {
    fn of(total: usize) -> Self {
        Self {
            total,
            left: total.saturating_sub(1) / 2,
        }
    }

    fn right(self) -> usize {
        self.total.saturating_sub(1) - self.left
    }

    /// The full sign vector: `+` first, then the left mask, then the right mask (bit set = `−`).
    fn signs(self, left_mask: u64, right_mask: u64) -> Vec<ContributionSign> {
        let sign = |mask: u64, bit: usize| {
            if mask >> bit & 1 == 1 {
                ContributionSign::Minus
            } else {
                ContributionSign::Plus
            }
        };
        let left = (0..self.left).map(|j| sign(left_mask, j));
        let right = (0..self.right()).map(|j| sign(right_mask, j));
        std::iter::once(ContributionSign::Plus)
            .chain(left)
            .chain(right)
            .collect()
    }
}

/// `2^len` as an entry count if it fits `usize`.
fn entry_count(len: usize) -> Option<usize> {
    u32::try_from(len)
        .ok()
        .and_then(|bits| 1usize.checked_shl(bits))
}

/// `2^len` as a `u64`, saturating at `u64::MAX` when it does not fit.
fn required_entries(len: usize) -> u64 {
    u32::try_from(len)
        .ok()
        .and_then(|bits| 1u64.checked_shl(bits))
        .unwrap_or(u64::MAX)
}

/// The meet-in-the-middle lists for one coefficient.
struct MeetInTheMiddle {
    /// `c_0 + Σ ±c_k` over the left signs, indexed by left mask.
    left: Vec<BigInt>,
    /// `(Σ ±c_k, mask)` over the right signs, ascending.
    right: Vec<(BigInt, u64)>,
    bound: BigInt,
    modulus: BigInt,
}

impl MeetInTheMiddle {
    /// Build both lists after checking their size against the budget, so nothing is allocated for
    /// a search that cannot fit.
    fn build(
        column: &[&BigInt],
        halves: Halves,
        bound: &BigInt,
        modulus: &BigInt,
        limits: RecombinationLimits,
    ) -> Result<Self, ResourceExhaustion> {
        let required = required_entries(halves.right());
        let exhausted = ResourceExhaustion {
            kind: ResourceKind::RecombinationList,
            consumed: required,
            limit: limits.max_list_entries,
        };
        if halves.right() >= 64 || required > limits.max_list_entries {
            return Err(exhausted);
        }
        let left_terms = &column[1..=halves.left];
        let right_terms = &column[1 + halves.left..];
        let left = signed_sums(column[0], left_terms, modulus).ok_or(exhausted)?;
        let sums = signed_sums(&BigInt::zero(), right_terms, modulus).ok_or(exhausted)?;
        let mut right = Vec::new();
        right.try_reserve_exact(sums.len()).map_err(|_| exhausted)?;
        right.extend(sums.into_iter().zip(0u64..));
        right.sort_unstable();
        Ok(Self {
            left,
            right,
            bound: bound.clone(),
            modulus: modulus.clone(),
        })
    }

    /// Every `(left mask, right mask)` whose sum is within the bound, in examination order.
    fn candidates(&self) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.left
            .iter()
            .zip(0u64..)
            .flat_map(move |(x, left_mask)| {
                residue_intervals(x, &self.bound, &self.modulus)
                    .into_iter()
                    .flat_map(move |(lo, hi)| {
                        let start = self.right.partition_point(|(r, _)| *r < lo);
                        let end = self.right.partition_point(|(r, _)| *r <= hi);
                        self.right[start..end]
                            .iter()
                            .map(move |&(_, right_mask)| (left_mask, right_mask))
                    })
            })
    }
}

/// `head + Σ ±term` modulo `modulus` for every sign mask, indexed by mask (bit `j` set means the
/// `j`-th term is subtracted). `None` if the `2^len` entries cannot be allocated.
fn signed_sums(head: &BigInt, terms: &[&BigInt], modulus: &BigInt) -> Option<Vec<BigInt>> {
    let mut sums = Vec::new();
    sums.try_reserve_exact(entry_count(terms.len())?).ok()?;
    sums.push(head.mod_floor(modulus));
    for term in terms {
        let minus: Vec<BigInt> = sums
            .iter()
            .map(|s| (s - *term).mod_floor(modulus))
            .collect();
        for s in &mut sums {
            *s = (&*s + *term).mod_floor(modulus);
        }
        sums.extend(minus);
    }
    Some(sums)
}

/// The residues `y` with `x + y` within `bound` of zero modulo `modulus`, as ascending disjoint
/// inclusive intervals in scan order: one interval, or two when the interval wraps past `modulus`.
///
/// Requires `2·bound + 1 ≤ modulus`, so the interval never covers a residue twice.
fn residue_intervals(x: &BigInt, bound: &BigInt, modulus: &BigInt) -> Vec<(BigInt, BigInt)> {
    let lo = (-x - bound).mod_floor(modulus);
    let hi = (bound - x).mod_floor(modulus);
    if lo <= hi {
        vec![(lo, hi)]
    } else {
        vec![(lo, modulus - BigInt::one()), (BigInt::zero(), hi)]
    }
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
    use shared_numfield::RootSearchLimits;

    use super::*;
    use crate::sqrt::{PrecisionLimits, RelationFactors, SqrtContext, TargetOutcome};

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn poly(coeffs: &[i64]) -> IntPoly {
        IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
    }

    fn is_invariant<T: std::fmt::Debug>(result: Result<T, SqrtError>) -> bool {
        matches!(result, Err(SqrtError::Invariant(_)))
    }

    /// Contributions for `f`, `m`, `N` and the relation pairs, at `ell` split primes from `first`.
    fn set_for(
        f: &[i64],
        m: i64,
        n: i64,
        pairs: &[(i64, i64)],
        first: u64,
        ell: u32,
    ) -> ContributionSet {
        let context = SqrtContext::new(&poly(f), &big(m), &big(n)).unwrap();
        let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
        let factors = RelationFactors::new(&context, &pairs).unwrap();
        let precision = PrecisionLimits {
            max_modulus_bits: 4096,
        };
        let TargetOutcome::Ready(target) = NormalizedTarget::new(factors, precision).unwrap()
        else {
            panic!("precision budget is generous");
        };
        let limits = CrtLimits {
            prime_search: PrimeSearchLimits {
                first_candidate: first,
                max_primes_examined: 200,
            },
            max_primes: ell,
            root_search: RootSearchLimits {
                max_evaluations: 100_000,
            },
            nonresidue_search: NonresidueSearchLimits {
                max_candidates: 1000,
            },
        };
        match split_contributions(&target, limits).unwrap() {
            SplitOutcome::Ready(set) => set,
            other => panic!("expected contributions, got {other:?}"),
        }
    }

    fn roomy() -> RecombinationLimits {
        RecombinationLimits {
            max_list_entries: 1 << 20,
            max_matches: 1 << 20,
        }
    }

    /// The coefficient index with the least bound, lowest index on ties, and that bound.
    fn narrowest(set: &ContributionSet) -> (usize, BigInt) {
        let bounds = set.target().coefficient_bounds();
        let kappa = (0..bounds.len())
            .min_by_key(|&i| (bounds[i].clone(), i))
            .unwrap();
        (kappa, bounds[kappa].clone())
    }

    /// Bit `k − 1` of the mask is set exactly when contribution `k` is subtracted.
    fn mask_of(signs: &[ContributionSign]) -> u64 {
        signs[1..]
            .iter()
            .enumerate()
            .filter(|&(_, s)| *s == ContributionSign::Minus)
            .fold(0, |mask, (j, _)| mask | 1 << j)
    }

    /// Every sign vector with the first sign `+` whose assembled coefficient `κ` is within its
    /// bound, found by assembling all `2^(n−1)` vectors.
    fn brute_force_matches(set: &ContributionSet) -> Vec<u64> {
        let (kappa, bound) = narrowest(set);
        let n = set.contributions().len();
        (0..1u64 << (n - 1))
            .filter(|&mask| {
                let signs: Vec<ContributionSign> = std::iter::once(ContributionSign::Plus)
                    .chain((0..n - 1).map(|j| {
                        if mask >> j & 1 == 1 {
                            ContributionSign::Minus
                        } else {
                            ContributionSign::Plus
                        }
                    }))
                    .collect();
                let z = set.assemble(&signs).unwrap();
                let value = z.coeffs.get(kappa).cloned().unwrap_or_else(BigInt::zero);
                -&bound <= value && value <= bound
            })
            .collect()
    }

    fn lists_for(set: &ContributionSet) -> (Halves, MeetInTheMiddle) {
        let (kappa, bound) = narrowest(set);
        let column: Vec<&BigInt> = set
            .contributions()
            .iter()
            .map(|c| &c.coefficients()[kappa])
            .collect();
        let halves = Halves::of(column.len());
        let lists =
            MeetInTheMiddle::build(&column, halves, &bound, set.modulus(), roomy()).unwrap();
        (halves, lists)
    }

    #[test]
    fn residue_intervals_split_a_wrapped_interval_into_two_exact_ranges() {
        // M = 11, b = 3: x = 0 needs y ∈ {−3, …, 3} = {8, 9, 10, 0, 1, 2, 3}, which wraps; x = 5
        // needs y ∈ {−8, …, −2} = {3, …, 9}, which does not.
        let (m, b) = (big(11), big(3));
        assert_eq!(
            residue_intervals(&big(0), &b, &m),
            vec![(big(8), big(10)), (big(0), big(3))]
        );
        assert_eq!(residue_intervals(&big(5), &b, &m), vec![(big(3), big(9))]);
        // The edges: the interval ends exactly at M − 1 (x = 8: y ∈ {0, …, 6}), starts exactly
        // at 0 (x = 3: y ∈ {5, …, 11 ≡ 0}), and with b = 0 is a single residue.
        assert_eq!(residue_intervals(&big(8), &b, &m), vec![(big(0), big(6))]);
        assert_eq!(
            residue_intervals(&big(3), &b, &m),
            vec![(big(5), big(10)), (big(0), big(0))]
        );
        assert_eq!(
            residue_intervals(&big(4), &big(0), &m),
            vec![(big(7), big(7))]
        );
        // M = 2b + 1 = 7: every residue is within the bound, exactly once, for every x.
        let (m, b) = (big(7), big(3));
        for x in 0..7 {
            let intervals = residue_intervals(&big(x), &b, &m);
            let mut covered: Vec<BigInt> = Vec::new();
            for (lo, hi) in &intervals {
                assert!(lo <= hi);
                let mut y = lo.clone();
                while y <= *hi {
                    covered.push(y.clone());
                    y += 1;
                }
            }
            covered.sort();
            assert_eq!(covered, (0..7).map(big).collect::<Vec<_>>(), "x = {x}");
        }
    }

    #[test]
    fn residue_intervals_cover_exactly_the_residues_within_the_bound() {
        // Brute force over all x and b for a small odd modulus.
        let m = big(13);
        for b in 0..=6 {
            for x in 0..13 {
                let expected: Vec<i64> = (0..13)
                    .filter(|y| {
                        let sum = centered(&big((x + y) % 13), &m);
                        big(-b) <= sum && sum <= big(b)
                    })
                    .collect();
                let mut got: Vec<i64> = Vec::new();
                for (lo, hi) in residue_intervals(&big(x), &big(b), &m) {
                    let (lo, hi) = (i64::try_from(lo).unwrap(), i64::try_from(hi).unwrap());
                    got.extend(lo..=hi);
                }
                got.sort_unstable();
                assert_eq!(got, expected, "x = {x}, b = {b}");
            }
        }
    }

    #[test]
    fn halves_put_the_larger_part_on_the_right_and_fix_the_first_sign() {
        let shapes: Vec<(usize, usize)> = [1, 2, 3, 4, 5, 8]
            .iter()
            .map(|&n| (Halves::of(n).left, Halves::of(n).right()))
            .collect();
        assert_eq!(shapes, vec![(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (3, 4)]);
        // n = 4: signs are [+, left bit 0, right bit 0, right bit 1].
        let signs = Halves::of(4).signs(0b1, 0b10);
        assert_eq!(
            signs,
            vec![
                ContributionSign::Plus,
                ContributionSign::Minus,
                ContributionSign::Plus,
                ContributionSign::Minus
            ]
        );
        assert_eq!(Halves::of(1).signs(0, 0), vec![ContributionSign::Plus]);
    }

    #[test]
    fn list_sizes_saturate_instead_of_overflowing() {
        assert_eq!(required_entries(0), 1);
        assert_eq!(required_entries(10), 1024);
        assert_eq!(required_entries(63), 1 << 63);
        assert_eq!(required_entries(64), u64::MAX);
        assert_eq!(required_entries(usize::MAX), u64::MAX);
        assert_eq!(entry_count(3), Some(8));
        assert_eq!(entry_count(usize::MAX), None);
        // A request that cannot fit is refused rather than allocated.
        assert_eq!(signed_sums(&big(0), &[&big(1); 70], &big(11)), None);
    }

    #[test]
    fn signed_sums_are_indexed_by_sign_mask() {
        // 5 ± 1 ± 2 ± 4 modulo 100, bit j set subtracting term j: mask 0b101 = 5 − 1 + 2 − 4 = 2,
        // mask 0b010 = 5 + 1 − 2 + 4 = 8, mask 0b111 = 5 − 7 = −2 ≡ 98.
        let terms = [big(1), big(2), big(4)];
        let refs: Vec<&BigInt> = terms.iter().collect();
        let sums = signed_sums(&big(5), &refs, &big(100)).unwrap();
        assert_eq!(sums.len(), 8);
        assert_eq!(sums[0], big(12));
        assert_eq!(sums[0b101], big(2));
        assert_eq!(sums[0b010], big(8));
        assert_eq!(sums[0b111], big(98));
        // No terms: just the head.
        assert_eq!(
            signed_sums(&big(-3), &[], &big(100)).unwrap(),
            vec![big(97)]
        );
    }

    #[test]
    fn meet_in_the_middle_finds_exactly_the_brute_force_matches() {
        let fixtures = [
            // 527: z = 4 + 2θ, one prime and two primes.
            set_for(&[-2, 0, 1], 23, 527, &[(3, -2)], 0, 1),
            set_for(&[-2, 0, 1], 23, 527, &[(3, -2)], 0, 2),
            // A negative integer is a local square at 7 and 23 and no global square.
            set_for(&[-2, 0, 1], 23, 527, &[(-5, 0)], 0, 1),
            set_for(&[-2, 0, 1], 23, 527, &[(-5, 0)], 0, 2),
            // V4 quartic at one prime and two, and the degree-one case with three primes.
            set_for(&[1, 0, -10, 0, 1], 4, 97, &[(1, -1), (1, -1)], 0, 1),
            set_for(&[1, 0, -10, 0, 1], 4, 97, &[(1, -1), (1, -1)], 0, 2),
            set_for(&[-5, 1], 5, 21, &[(169, 0)], 0, 3),
        ];
        for (index, set) in fixtures.iter().enumerate() {
            let (halves, lists) = lists_for(set);
            let mut found: Vec<u64> = lists
                .candidates()
                .map(|(l, r)| mask_of(&halves.signs(l, r)))
                .collect();
            found.sort_unstable();
            let mut distinct = found.clone();
            distinct.dedup();
            assert_eq!(found, distinct, "fixture {index}: a match was listed twice");
            assert_eq!(found, brute_force_matches(set), "fixture {index}");
        }
    }

    #[test]
    fn matches_inside_a_wrapped_interval_are_not_lost() {
        // A left sum within b of zero modulo M has a residue interval that wraps past M: it is
        // [lo, M − 1] followed by [0, hi], with lo > hi. Scanning only [lo, hi] when lo ≤ hi would
        // silently drop the matches of such a sum.
        let wrapped = [
            set_for(&[-2, 0, 1], 23, 527, &[(-5, 0)], 0, 2),
            set_for(&[1, 0, -10, 0, 1], 4, 97, &[(1, -1), (1, -1)], 0, 2),
        ];
        for set in &wrapped {
            let (_, lists) = lists_for(set);
            let wrapping = lists
                .left
                .iter()
                .filter(|x| residue_intervals(x, &lists.bound, &lists.modulus).len() == 2)
                .count();
            assert!(wrapping > 0, "the fixture must exercise a wrapped interval");
            let unwrapped_only = lists
                .left
                .iter()
                .map(|x| {
                    let lo = (-x - &lists.bound).mod_floor(&lists.modulus);
                    let hi = (&lists.bound - x).mod_floor(&lists.modulus);
                    if lo <= hi {
                        lists
                            .right
                            .iter()
                            .filter(|(r, _)| lo <= *r && *r <= hi)
                            .count()
                    } else {
                        0
                    }
                })
                .sum::<usize>();
            assert!(unwrapped_only < lists.candidates().count());
        }
    }

    #[test]
    fn recombination_reports_defective_sets_as_invariant_errors() {
        let set = set_for(&[-2, 0, 1], 23, 527, &[(3, -2)], 0, 2);
        assert!(recombine_signs(&set, roomy()).is_ok());

        let mut empty = set.clone();
        empty.contributions.clear();
        assert!(is_invariant(recombine_signs(&empty, roomy())));

        // A modulus at most 2b + 1 would let an interval cover a residue twice: b = 6, so 13.
        let mut small = set.clone();
        small.modulus = big(12);
        assert!(is_invariant(recombine_signs(&small, roomy())));

        let mut short = set.clone();
        short.contributions[2].coefficients.truncate(1);
        assert!(is_invariant(recombine_signs(&short, roomy())));
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
