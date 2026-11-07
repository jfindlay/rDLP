//! Direct p-adic reconstruction of an algebraic square root at one inert prime.
//!
//! For a [`NormalizedTarget`] with monic `F` of degree `d` and target `H`, a square root `z` of `H`
//! modulo `F` is sought. At an odd prime `p` coprime to `c·N` with `F` irreducible modulo `p` and
//! `H ≢ 0`, the quotient `(Z/p^e)[t]/(F)` is a local ring whose unit square roots come in exactly
//! one pair `±z`. The algorithm is:
//!
//! 1. take a square root `s` of `H` in the residue field `F_(p^d)` and set `u = s⁻¹`;
//! 2. lift `u` by `u ← u(3 − H·u²)/2`, which keeps `H·u² ≡ 1` and doubles the precision;
//! 3. center `z = H·u` modulo the least `p^e > 2B`, where `B` is the certified coefficient bound;
//! 4. hand `z` to [`NormalizedTarget::verify_candidate`], the only way to obtain a
//!    [`VerifiedRoot`](super::VerifiedRoot).
//!
//! A nonresidue in the residue field, or a centered `z` that the verifier rejects at `M > 2B`,
//! certifies that `H` is not a square. Exhausting a budget certifies nothing. The quotient-ring
//! arithmetic below is private to this module. See `gnfs/docs/PEDAGOGY.md` for the derivation.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::pow::pow;
use num_traits::{One, Zero};
use shared_numfield::{
    FpPoly, IntPoly, NonresidueSearchLimits, PrimeModulus, ResidueField, SqrtOutcome, is_prime_u64,
};

use super::{
    AlgebraicSqrtOutcome, CandidateRejection, CandidateVerdict, DirectLimits, InvariantViolation,
    NonSquareCertificate, NonSquareEvidence, NormalizedTarget, PrimeSearchLimits,
    ReconstructionMethod, ResourceExhaustion, ResourceKind, SqrtError,
};

/// Reconstruct a verified square root of `target` at the first usable inert prime.
///
/// Primes `p ≥ max(first_candidate, 3)` are examined in ascending order; every prime examined
/// counts against the budget. A prime is skipped when it divides `c·N`, when `F` is not irreducible
/// modulo `p`, or when `H ≡ 0 (mod p)`. At the first usable prime the root is lifted to the least
/// `p^e > 2B` and verified exactly.
///
/// Outcomes:
/// - [`AlgebraicSqrtOutcome::Verified`] with [`ReconstructionMethod::DirectNewton`];
/// - [`AlgebraicSqrtOutcome::NotSquare`] with [`NonSquareEvidence::InertLocalNonresidue`] when `H`
///   is a nonresidue in the residue field, or [`NonSquareEvidence::DirectBoundExhausted`] when the
///   centered candidate at `M > 2B` is rejected;
/// - [`AlgebraicSqrtOutcome::Incidental`] / [`AlgebraicSqrtOutcome::Unsupported`] when a root
///   exists but the map to `Z/N` fails;
/// - [`AlgebraicSqrtOutcome::ResourceExhausted`] with one entry: [`ResourceKind::InertPrimeSearch`]
///   (budget spent), [`ResourceKind::MachinePrimeRange`] (no prime left below `2^64`; `limit` is
///   the configured budget), [`ResourceKind::NonresidueSearch`], or [`ResourceKind::Precision`]
///   (the exponent would not fit a `u32`). None of these says anything about squareness.
///
/// The precision budget is enforced when the target is built: the lifting modulus has at most
/// `bits(2B + 1) + bits(p)` bits.
///
/// # Errors
///
/// [`SqrtError::Invariant`] if an identity that must hold fails (a lifting identity, a prime that
/// is not prime, or a reduced candidate of degree `≥ d`). These indicate a defect, never a
/// nonsquare.
pub fn direct_root(
    target: &NormalizedTarget,
    limits: DirectLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    let mut cursor = PrimeCursor::new(limits.prime_search);
    loop {
        let prime = match cursor.advance() {
            Ok(prime) => prime,
            Err(exhausted) => return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![exhausted])),
        };
        if let Screen::Inert { field, residue } = screen_prime(target, prime)? {
            return reconstruct_at(target, prime, &field, &residue, limits.nonresidue_search);
        }
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
            return Err(self.exhaustion(ResourceKind::InertPrimeSearch));
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

/// Why a prime was or was not usable for direct reconstruction.
enum Screen {
    /// `p` is even.
    Even,
    /// `p` divides the leading coefficient or the modulus `N`.
    DividesInput,
    /// `F` is not irreducible modulo `p`.
    NotInert,
    /// `H ≡ 0 (mod p)`, so `H` is not a unit.
    ZeroTarget,
    /// `F` is irreducible modulo `p` and `H` is a nonzero element of the residue field.
    Inert {
        field: ResidueField,
        residue: FpPoly,
    },
}

/// Classify the prime `p` for `target`. Rejects `H ≡ 0` before any square root is attempted.
fn screen_prime(target: &NormalizedTarget, p: u64) -> Result<Screen, SqrtError> {
    if p.is_multiple_of(2) {
        return Ok(Screen::Even);
    }
    let modulus = PrimeModulus::new(p).map_err(|_| violation("candidate prime is not prime"))?;
    let context = target.context();
    let p_big = BigInt::from(p);
    if context.leading_coeff().is_multiple_of(&p_big) || context.modulus().is_multiple_of(&p_big) {
        return Ok(Screen::DividesInput);
    }
    let monic = FpPoly::from_int_poly(context.monic_poly(), modulus);
    if monic.degree() != Some(context.degree()) {
        return Err(violation(
            "monic polynomial lost its degree modulo an odd prime",
        ));
    }
    if !monic.is_irreducible() {
        return Ok(Screen::NotInert);
    }
    let residue = FpPoly::from_int_poly(target.target(), modulus);
    if residue.is_zero() {
        return Ok(Screen::ZeroTarget);
    }
    let field = ResidueField::new(&monic).map_err(|_| violation("irreducible modulus rejected"))?;
    Ok(Screen::Inert { field, residue })
}

// ─── reconstruction at one prime ─────────────────────────────────────────────

/// Take the residue-field square root, lift it, center it, and verify it.
fn reconstruct_at(
    target: &NormalizedTarget,
    prime: u64,
    field: &ResidueField,
    residue: &FpPoly,
    nonresidue: NonresidueSearchLimits,
) -> Result<AlgebraicSqrtOutcome, SqrtError> {
    let root = match field
        .sqrt(residue, nonresidue)
        .map_err(|_| violation("residue-field square root failed"))?
    {
        SqrtOutcome::Root(root) => root,
        SqrtOutcome::NotSquare => {
            return Ok(AlgebraicSqrtOutcome::NotSquare(NonSquareCertificate::new(
                NonSquareEvidence::InertLocalNonresidue { prime },
            )));
        }
        SqrtOutcome::ResourceExhausted { candidates_tried } => {
            return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![
                ResourceExhaustion {
                    kind: ResourceKind::NonresidueSearch,
                    consumed: candidates_tried,
                    limit: nonresidue.max_candidates,
                },
            ]));
        }
    };
    let seed_poly = field
        .inv(&root)
        .map_err(|_| violation("square root of a unit is not invertible"))?;
    let degree = target.context().degree();
    let seed: Element = (0..degree)
        .map(|i| BigInt::from(seed_poly.coeffs().get(i).copied().unwrap_or(0)))
        .collect();

    let minimum = target.minimum_modulus();
    let Some(exponent) = lifting_exponent(prime, minimum) else {
        return Ok(AlgebraicSqrtOutcome::ResourceExhausted(vec![
            ResourceExhaustion {
                kind: ResourceKind::Precision,
                consumed: minimum.bits(),
                limit: u64::from(u32::MAX),
            },
        ]));
    };
    let monic = target.context().monic_poly();
    let stages = lift_inverse_sqrt(target.target(), monic, prime, seed, exponent)?;
    let last = match stages.last() {
        Some(stage) if stage.exponent == exponent => stage,
        _ => return Err(violation("lifting did not reach the requested precision")),
    };

    let modulus = pow(BigInt::from(prime), exponent as usize);
    let ring = PrimePowerRing::new(monic, modulus.clone());
    let product = ring.mul(&ring.element(target.target()), &last.inverse_root);
    let z = IntPoly::from_coeffs(product.iter().map(|c| centered(c, &modulus)).collect());

    let method = ReconstructionMethod::DirectNewton { prime, exponent };
    Ok(match target.verify_candidate(&z, method)? {
        CandidateVerdict::Verified(root) => AlgebraicSqrtOutcome::Verified(root),
        CandidateVerdict::Rejected(CandidateRejection::DegreeTooHigh) => {
            return Err(violation("reduced candidate has degree at least d"));
        }
        // At M > 2B the centered residue equals the unique ±root when one exists, so a rejection
        // here proves that H has no square root.
        CandidateVerdict::Rejected(_) => AlgebraicSqrtOutcome::NotSquare(
            NonSquareCertificate::new(NonSquareEvidence::DirectBoundExhausted { prime, exponent }),
        ),
        CandidateVerdict::Incidental(factor) => AlgebraicSqrtOutcome::Incidental(factor),
        CandidateVerdict::Unsupported(domain) => AlgebraicSqrtOutcome::Unsupported(domain),
    })
}

/// The least `e ≥ 1` with `p^e ≥ minimum`, or `None` if `e` would not fit a `u32`.
fn lifting_exponent(prime: u64, minimum: &BigInt) -> Option<u32> {
    let p = BigInt::from(prime);
    let mut power = p.clone();
    let mut exponent: u32 = 1;
    while power < *minimum {
        power *= &p;
        exponent = exponent.checked_add(1)?;
    }
    Some(exponent)
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

// ─── Newton lifting in (Z/p^e)[t]/(F) ────────────────────────────────────────

/// A residue of `(Z/p^e)[t]/(F)`: exactly `d` coefficients, each in `[0, p^e)`.
type Element = Vec<BigInt>;

/// An inverse square root known modulo `p^exponent`.
struct LiftStage {
    exponent: u32,
    inverse_root: Element,
}

/// Lift the inverse square root `seed` (valid modulo `p`) to precision `p^final_exponent`.
///
/// Each stage doubles the previous precision up to the ceiling-halving chain of `final_exponent`.
/// The returned stages run from exponent 1 to `final_exponent` and each satisfies `H·u² ≡ 1`
/// modulo its own `p^exponent`.
fn lift_inverse_sqrt(
    h: &IntPoly,
    monic: &IntPoly,
    prime: u64,
    seed: Element,
    final_exponent: u32,
) -> Result<Vec<LiftStage>, SqrtError> {
    let p = BigInt::from(prime);
    let chain = exponent_chain(final_exponent);
    let mut previous_modulus = p.clone();
    let mut current = seed;
    let mut stages = vec![LiftStage {
        exponent: 1,
        inverse_root: current.clone(),
    }];
    for &exponent in chain.iter().skip(1) {
        let modulus = pow(p.clone(), exponent as usize);
        let ring = PrimePowerRing::new(monic, modulus.clone());
        current = ring.newton_step(&ring.element(h), &current, &previous_modulus)?;
        stages.push(LiftStage {
            exponent,
            inverse_root: current.clone(),
        });
        previous_modulus = modulus;
    }
    let ring = PrimePowerRing::new(monic, previous_modulus.clone());
    if !ring.is_unit_square(&ring.element(h), &current, &previous_modulus) {
        return Err(violation(
            "lifted inverse square root does not satisfy H·u² ≡ 1",
        ));
    }
    Ok(stages)
}

/// Exponents `[1, …, e]` where each is the ceiling half of the next, so each step at most doubles.
fn exponent_chain(final_exponent: u32) -> Vec<u32> {
    let mut chain = vec![final_exponent];
    let mut exponent = final_exponent;
    while exponent > 1 {
        exponent = exponent.div_ceil(2);
        chain.push(exponent);
    }
    chain.reverse();
    chain
}

/// Arithmetic in `(Z/M)[t]/(F)` for an odd modulus `M = p^e` and monic `F` of degree `d`.
struct PrimePowerRing {
    modulus: BigInt,
    /// `F_0, …, F_(d−1)` reduced modulo `M`; the leading coefficient is 1.
    monic_tail: Vec<BigInt>,
}

impl PrimePowerRing {
    fn new(monic: &IntPoly, modulus: BigInt) -> Self {
        let degree = monic.coeffs.len().saturating_sub(1);
        let monic_tail = monic
            .coeffs
            .iter()
            .take(degree)
            .map(|c| c.mod_floor(&modulus))
            .collect();
        Self {
            modulus,
            monic_tail,
        }
    }

    fn degree(&self) -> usize {
        self.monic_tail.len()
    }

    /// The residue class of an integer polynomial of degree `< d`.
    fn element(&self, poly: &IntPoly) -> Element {
        (0..self.degree())
            .map(|i| {
                poly.coeffs
                    .get(i)
                    .map_or_else(BigInt::zero, |c| c.mod_floor(&self.modulus))
            })
            .collect()
    }

    /// `a·b` reduced modulo `F` and `M`.
    fn mul(&self, a: &Element, b: &Element) -> Element {
        let d = self.degree();
        let mut acc = vec![BigInt::zero(); 2 * d - 1];
        for (i, x) in a.iter().enumerate().filter(|(_, x)| !x.is_zero()) {
            for (j, y) in b.iter().enumerate() {
                acc[i + j] += x * y;
            }
        }
        // Eliminate t^k for k = 2d−2 … d using t^d = −(F_0 + … + F_(d−1) t^(d−1)).
        for k in (d..2 * d - 1).rev() {
            let lead = acc[k].mod_floor(&self.modulus);
            for (j, f_j) in self.monic_tail.iter().enumerate() {
                acc[k - d + j] -= &lead * f_j;
            }
        }
        acc.truncate(d);
        acc.iter().map(|c| c.mod_floor(&self.modulus)).collect()
    }

    /// `3 − a`.
    fn three_minus(&self, a: &Element) -> Element {
        a.iter()
            .enumerate()
            .map(|(i, c)| {
                let constant = if i == 0 {
                    BigInt::from(3)
                } else {
                    BigInt::zero()
                };
                (constant - c).mod_floor(&self.modulus)
            })
            .collect()
    }

    /// True when `a ≡ 1` modulo `known`, which divides `M`.
    fn is_one_modulo(a: &Element, known: &BigInt) -> bool {
        a.iter().enumerate().all(|(i, c)| {
            let expected = if i == 0 {
                BigInt::one()
            } else {
                BigInt::zero()
            };
            (c - expected).is_multiple_of(known)
        })
    }

    /// `H·u² ≡ 1` modulo `known`.
    fn is_unit_square(&self, h: &Element, u: &Element, known: &BigInt) -> bool {
        let square = self.mul(h, &self.mul(u, u));
        Self::is_one_modulo(&square, known)
    }

    /// One step `u ← u(3 − H·u²)/2` from a root valid modulo `previous` to one valid modulo `M`.
    ///
    /// With `H·u² = 1 + ε` the result satisfies `H·u'² = 1 − (3/4)ε² + ε³/4`, so the error `ε`
    /// squares: precision doubles. Dividing by 2 multiplies by `(M + 1)/2` because `M` is odd.
    fn newton_step(
        &self,
        h: &Element,
        u: &Element,
        previous: &BigInt,
    ) -> Result<Element, SqrtError> {
        let square = self.mul(h, &self.mul(u, u));
        if !Self::is_one_modulo(&square, previous) {
            return Err(violation(
                "inverse square root lost its precision before a lift",
            ));
        }
        let corrected = self.mul(u, &self.three_minus(&square));
        let half: BigInt = (&self.modulus + 1) / 2;
        Ok(corrected
            .iter()
            .map(|c| (c * &half).mod_floor(&self.modulus))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sqrt::{PrecisionLimits, RelationFactors, SqrtContext, TargetOutcome};

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn poly(coeffs: &[i64]) -> IntPoly {
        IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
    }

    fn target(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> NormalizedTarget {
        let ctx = SqrtContext::new(&poly(f), &big(m), &big(n)).unwrap();
        let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
        let factors = RelationFactors::new(&ctx, &pairs).unwrap();
        match NormalizedTarget::new(
            factors,
            PrecisionLimits {
                max_modulus_bits: 4096,
            },
        )
        .unwrap()
        {
            TargetOutcome::Ready(target) => target,
            TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
        }
    }

    fn kind(screen: &Screen) -> &'static str {
        match screen {
            Screen::Even => "even",
            Screen::DividesInput => "divides",
            Screen::NotInert => "reducible",
            Screen::ZeroTarget => "zero",
            Screen::Inert { .. } => "inert",
        }
    }

    /// `a·b mod F` over `Z`, reduced modulo `modulus`; independent of `PrimePowerRing`.
    fn product_mod(a: &IntPoly, b: &IntPoly, f: &IntPoly, modulus: &BigInt) -> IntPoly {
        let remainder = a.mul(b).pseudo_div_rem(f).1;
        IntPoly::from_coeffs(
            remainder
                .coeffs
                .iter()
                .map(|c| c.mod_floor(modulus))
                .collect(),
        )
    }

    #[test]
    fn newton_doubles_precision_keeping_unit_square_invariant() {
        // F = t² − 2 and H = 24 + 16θ (the target of the relation (3, −2)). Modulo 3, H ≡ θ.
        // In F_9 = F_3[t]/(t²−2), s = 1 + 2θ has s² = 1 + 4θ + 4·2 = 9 + 4θ ≡ θ, and
        // s⁻¹ = (1 − 2θ)/Norm(s) with Norm(s) = 1 − 8 ≡ 2, so u₀ = 2(1 − 2θ) ≡ 2 + 2θ.
        // Check: u₀² = 4 + 8θ + 8 ≡ 2θ and θ·2θ = 4 ≡ 1, so H·u₀² ≡ 1 (mod 3).
        let f = poly(&[-2, 0, 1]);
        let h = poly(&[24, 16]);
        let seed = vec![big(2), big(2)];
        let stages = lift_inverse_sqrt(&h, &f, 3, seed, 5).unwrap();
        let exponents: Vec<u32> = stages.iter().map(|s| s.exponent).collect();
        assert_eq!(exponents, vec![1, 2, 3, 5]);
        for stage in &stages {
            let modulus = pow(big(3), stage.exponent as usize);
            let u = IntPoly::from_coeffs(stage.inverse_root.clone());
            let u2 = product_mod(&u, &u, &f, &modulus);
            let unit = product_mod(&h, &u2, &f, &modulus);
            assert_eq!(unit, poly(&[1]), "H·u² ≢ 1 modulo 3^{}", stage.exponent);
        }
        // z = H·u is ±(4 + 2θ), the root of the relation, at every precision.
        let last = stages.last().unwrap();
        let modulus = pow(big(3), 5);
        let u = IntPoly::from_coeffs(last.inverse_root.clone());
        let z = product_mod(&h, &u, &f, &modulus);
        let centered_z =
            IntPoly::from_coeffs(z.coeffs.iter().map(|c| centered(c, &modulus)).collect());
        assert!(centered_z == poly(&[4, 2]) || centered_z == poly(&[-4, -2]));
    }

    #[test]
    fn exponent_chain_halves_with_ceiling() {
        assert_eq!(exponent_chain(1), vec![1]);
        assert_eq!(exponent_chain(2), vec![1, 2]);
        assert_eq!(exponent_chain(4), vec![1, 2, 4]);
        assert_eq!(exponent_chain(5), vec![1, 2, 3, 5]);
        assert_eq!(exponent_chain(9), vec![1, 2, 3, 5, 9]);
    }

    #[test]
    fn lifting_exponent_is_least_with_power_at_least_minimum() {
        // Powers of 3: 3, 9, 27, 81, 243. The minimum is 2B + 1, so equality is allowed.
        assert_eq!(lifting_exponent(3, &big(2)), Some(1));
        assert_eq!(lifting_exponent(3, &big(3)), Some(1));
        assert_eq!(lifting_exponent(3, &big(4)), Some(2));
        assert_eq!(lifting_exponent(3, &big(27)), Some(3));
        assert_eq!(lifting_exponent(3, &big(29)), Some(4));
        assert_eq!(lifting_exponent(3, &big(37)), Some(4));
        assert_eq!(lifting_exponent(3, &big(243)), Some(5));
    }

    #[test]
    fn lifting_modulus_is_less_than_one_prime_above_the_minimum() {
        // p^(e−1) < minimum ≤ p^e implies p^e < p·minimum, so bits(p^e) ≤ bits(minimum) + bits(p).
        for (prime, minimum) in [
            (3u64, 37i64),
            (5, 37),
            (7, 1000),
            (101, 2),
            (1009, 10_000_000),
        ] {
            let minimum = big(minimum);
            let e = lifting_exponent(prime, &minimum).unwrap();
            let power = pow(BigInt::from(prime), e as usize);
            assert!(power >= minimum);
            assert!(power.bits() <= minimum.bits() + BigInt::from(prime).bits());
        }
    }

    #[test]
    fn centered_residue_is_symmetric_for_odd_modulus() {
        // M = 27: residues 0…13 stay, 14…26 map to −13…−1.
        assert_eq!(centered(&big(13), &big(27)), big(13));
        assert_eq!(centered(&big(14), &big(27)), big(-13));
        assert_eq!(centered(&big(26), &big(27)), big(-1));
        assert_eq!(centered(&big(0), &big(27)), big(0));
    }

    #[test]
    fn prime_cursor_counts_only_primes_and_skips_two() {
        let mut cursor = PrimeCursor::new(PrimeSearchLimits {
            first_candidate: 0,
            max_primes_examined: 3,
        });
        assert_eq!(cursor.advance(), Ok(3));
        assert_eq!(cursor.advance(), Ok(5));
        assert_eq!(cursor.advance(), Ok(7));
        let exhausted = cursor.advance().unwrap_err();
        assert_eq!(exhausted.kind, ResourceKind::InertPrimeSearch);
        assert_eq!((exhausted.consumed, exhausted.limit), (3, 3));

        let mut cursor = PrimeCursor::new(PrimeSearchLimits {
            first_candidate: 8,
            max_primes_examined: 2,
        });
        assert_eq!(cursor.advance(), Ok(11));
        assert_eq!(cursor.advance(), Ok(13));
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
        // Once the range has ended it stays ended.
        assert_eq!(
            cursor.advance().unwrap_err().kind,
            ResourceKind::MachinePrimeRange
        );
    }

    #[test]
    fn screen_excludes_even_dividing_reducible_and_zero_primes() {
        // F = t² + t + 1 (f = x² + x + 1, m = 2, N = 7), relation (1, 1) twice: H = 9θ.
        let quadratic = target(&[1, 1, 1], 2, 7, &[(1, 1), (1, 1)]);
        assert_eq!(kind(&screen_prime(&quadratic, 2).unwrap()), "even");
        // t² + t + 1 = (t − 1)² mod 3 is reducible; mod 7 it splits (7 | N as well).
        assert_eq!(kind(&screen_prime(&quadratic, 3).unwrap()), "reducible");
        assert_eq!(kind(&screen_prime(&quadratic, 7).unwrap()), "divides");
        // Mod 5 the discriminant −3 ≡ 2 is a nonresidue, so F is irreducible, and 9θ ≡ 4θ ≠ 0.
        assert_eq!(kind(&screen_prime(&quadratic, 5).unwrap()), "inert");

        // F = t² − 2 with relation (3, 0) twice: H = 72, which vanishes mod 3 but not mod 5.
        let zero_mod_three = target(&[-2, 0, 1], 3, 7, &[(3, 0), (3, 0)]);
        assert_eq!(kind(&screen_prime(&zero_mod_three, 3).unwrap()), "zero");
        assert_eq!(kind(&screen_prime(&zero_mod_three, 5).unwrap()), "inert");

        // d = 1 keeps F irreducible at every prime, so only the divisibility guard excludes p.
        // f = x − 5 with N = 21 = 3·7 excludes 3 and 7; f = 3x − 15 (c = 3, N = 7) excludes 3.
        let linear = target(&[-5, 1], 5, 21, &[(169, 0)]);
        assert_eq!(kind(&screen_prime(&linear, 3).unwrap()), "divides");
        assert_eq!(kind(&screen_prime(&linear, 7).unwrap()), "divides");
        assert_eq!(kind(&screen_prime(&linear, 5).unwrap()), "inert");
        let scaled = target(&[-15, 3], 5, 7, &[(169, 0)]);
        assert_eq!(kind(&screen_prime(&scaled, 3).unwrap()), "divides");
    }
}
