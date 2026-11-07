//! Normalized algebraic target, certified coefficient bounds, and the exact root verifier.
//!
//! For relation pairs `(a, b)` over a [`SqrtContext`] with `c = lc(f)`, `θ = c·α`, and `F` monic,
//! write `k` for the relation count and `h = ⌈k/2⌉`. A [`NormalizedTarget`] holds
//!
//! - `γ = Π(a − b·α) = Π(c·a − b·θ) / c^k` in `A = Q[t]/(F)`;
//! - `η = c^(2h−k) Π(c·a − b·θ) = c^(2h) γ` and `H = F'(θ)² η`, both integral;
//! - the verified inverse of `F'` modulo `F`, and the certified coefficient bounds `B_i`.
//!
//! Every square root `β` of `γ` yields an integral `z = c^h F'(θ) β` with `z² ≡ H (mod F)` and
//! `|z_i| ≤ B_i`. [`NormalizedTarget::verify_candidate`] is the only way to obtain a
//! [`VerifiedRoot`]: it re-derives `β` from a proposed `z`, checks `β² = γ` over `Q` exactly,
//! checks that every denominator of `β` is a unit modulo `N`, and evaluates `y = β(c·m) mod N`,
//! which satisfies `y² ≡ P = Π(a − b·m) (mod N)`.
//!
//! This module is a child of `context` so that it can call the private `VerifiedRoot::seal`;
//! `verify_candidate` is its single call site. See `gnfs/docs/PEDAGOGY.md` for the derivation of
//! the bounds and the map.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::pow::pow;
use num_traits::{One, Signed, Zero};
use shared_numfield::{IntPoly, RatPoly};

use super::{
    IncidentalFactor, IncidentalSource, InvariantViolation, PrecisionLimits, ReconstructionMethod,
    RelationFactors, ResourceExhaustion, ResourceKind, RootParts, SqrtContext, SqrtError,
    UnsupportedDomain, VerifiedRoot, powers_of,
};

/// The target square `H` together with everything needed to bound and verify its roots.
///
/// Built once per relation set by [`NormalizedTarget::new`]; owns its [`RelationFactors`] (and so
/// its [`SqrtContext`]), so a target cannot be paired with a different polynomial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedTarget {
    factors: RelationFactors,
    half: usize,
    rational_product: BigInt,
    gamma: RatPoly,
    eta: IntPoly,
    target: IntPoly,
    derivative: IntPoly,
    derivative_inverse: RatPoly,
    bounds: Vec<BigInt>,
    bound: BigInt,
    minimum_modulus: BigInt,
}

/// Result of building a [`NormalizedTarget`] within a precision budget.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TargetOutcome {
    /// The target was built and its minimum modulus fits the budget.
    Ready(NormalizedTarget),
    /// The modulus needed to certify a root exceeds the budget; this is never evidence about
    /// squareness.
    ResourceExhausted(ResourceExhaustion),
}

/// Why a proposed `z` is not a square root of the target. Never a certificate of nonsquareness.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateRejection {
    /// `deg z >= d`: not a reduced representative.
    DegreeTooHigh,
    /// `|z_i|` exceeds its certified bound `B_i`, so `z` cannot be any root.
    CoefficientOutOfBound {
        /// The first coefficient that exceeds its bound.
        index: usize,
    },
    /// `z² ≢ H (mod F)` over `Z`.
    NotSquareRootOfTarget,
}

/// The result of checking a proposed root.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateVerdict {
    /// `z` is a root and its image modulo `N` was computed and checked.
    Verified(VerifiedRoot),
    /// `z` is not a root of the target.
    Rejected(CandidateRejection),
    /// A root exists but a denominator shares a proper factor with `N`.
    Incidental(IncidentalFactor),
    /// A root exists but a denominator is divisible by `N`, so the coefficientwise map fails.
    Unsupported(UnsupportedDomain),
}

impl NormalizedTarget {
    /// Build `γ`, `η`, `H`, the verified inverse of `F'`, and the coefficient bounds.
    ///
    /// The precision budget is checked before the inverse is computed and before any caller
    /// searches for primes: if `bits(2B + 1) > limits.max_modulus_bits` the outcome is
    /// [`TargetOutcome::ResourceExhausted`] with kind [`ResourceKind::Precision`].
    ///
    /// # Errors
    ///
    /// [`SqrtError::Invariant`] if `F'` is not invertible modulo `F` or the computed inverse fails
    /// its identity. Both are impossible for a validated context, whose `F` is squarefree.
    pub fn new(
        factors: RelationFactors,
        limits: PrecisionLimits,
    ) -> Result<TargetOutcome, SqrtError> {
        let ctx = factors.context();
        let monic = ctx.monic_poly();
        let count = factors.relation_count();
        let half = count.div_ceil(2);
        let leading = ctx.leading_coeff();

        let scaled_product = scaled_factor_product(ctx, factors.pairs());
        let eta = scaled_product.scale(&pow(leading.clone(), 2 * half - count));
        let gamma = scaled_product
            .to_rat_poly()
            .scale(&BigRational::new(BigInt::one(), pow(leading.clone(), count)));
        let derivative = monic.derivative();
        let target = reduce_monic(&derivative.mul(&derivative).mul(&eta), monic);
        let (bounds, bound) = coefficient_bounds(monic, &eta);
        let minimum_modulus: BigInt = &bound * 2 + 1;

        let required_bits = minimum_modulus.bits();
        if required_bits > limits.max_modulus_bits {
            return Ok(TargetOutcome::ResourceExhausted(ResourceExhaustion {
                kind: ResourceKind::Precision,
                consumed: required_bits,
                limit: limits.max_modulus_bits,
            }));
        }
        let derivative_inverse = verified_inverse(&derivative, monic)?;
        let rational_product: BigInt =
            factors.pairs().iter().map(|(a, b)| a - b * ctx.root()).product();
        Ok(TargetOutcome::Ready(Self {
            factors,
            half,
            rational_product,
            gamma,
            eta,
            target,
            derivative,
            derivative_inverse,
            bounds,
            bound,
            minimum_modulus,
        }))
    }

    /// Check a proposed root `z` of `H` and, if it is one, map it to a residue modulo `N`.
    ///
    /// Checks run in order: `deg z < d`; `|z_i| ≤ B_i`; `z² ≡ H (mod F)` over `Z` (these three
    /// are rejections). Then `β = z / (c^h F')` and `β² = γ` over `Q`; every reduced denominator
    /// of `β` is a unit modulo `N`; `y = β(c·m) mod N`; `y² ≡ P (mod N)`. If any denominator
    /// shares a proper factor with `N` the result is [`CandidateVerdict::Incidental`] (preferred
    /// over an unsupported denominator found earlier); if only whole-`N` denominators occur it is
    /// [`CandidateVerdict::Unsupported`]. Both mean that a root exists but the map fails.
    ///
    /// # Errors
    ///
    /// [`SqrtError::Invariant`] when `β² ≠ γ` or `y² ≢ P (mod N)` after `z² ≡ H` passed. Both
    /// follow from the earlier checks for an untampered target, so they indicate a defect and are
    /// never reported as a rejection.
    pub fn verify_candidate(
        &self,
        z: &IntPoly,
        method: ReconstructionMethod,
    ) -> Result<CandidateVerdict, SqrtError> {
        if let Some(rejection) = self.reject_candidate(z) {
            return Ok(CandidateVerdict::Rejected(rejection));
        }
        let beta = self.recover_beta(z)?;
        let image = match self.modular_image(&beta) {
            Ok(image) => image,
            Err(blocked) => return Ok(blocked),
        };
        let ctx = self.context();
        if !(&image * &image - &self.rational_product).mod_floor(ctx.modulus()).is_zero() {
            return Err(violation("modular image does not square to the rational product"));
        }
        Ok(CandidateVerdict::Verified(VerifiedRoot::seal(RootParts {
            context: ctx.clone(),
            relation_count: self.relation_count(),
            rational_product: self.rational_product.clone(),
            integral_root: z.clone(),
            beta,
            modular_image: image,
            method,
        })))
    }

    /// The owning relation set.
    #[must_use]
    pub fn factors(&self) -> &RelationFactors {
        &self.factors
    }

    /// The owning context.
    #[must_use]
    pub fn context(&self) -> &SqrtContext {
        self.factors.context()
    }

    /// The number of relations `k`.
    #[must_use]
    pub fn relation_count(&self) -> usize {
        self.factors.relation_count()
    }

    /// `h = ⌈k/2⌉`.
    #[must_use]
    pub fn half_count(&self) -> usize {
        self.half
    }

    /// `P = Π(a − b·m)` over `Z`, with the exact `m`.
    #[must_use]
    pub fn rational_product(&self) -> &BigInt {
        &self.rational_product
    }

    /// `γ = Π(a − b·α)` in the `θ` basis, reduced modulo `F`.
    #[must_use]
    pub fn gamma(&self) -> &RatPoly {
        &self.gamma
    }

    /// `η = c^(2h) γ`, integral and reduced modulo `F`.
    #[must_use]
    pub fn eta(&self) -> &IntPoly {
        &self.eta
    }

    /// `H = F'(θ)² η`, integral and reduced modulo `F`: the polynomial whose root is sought.
    #[must_use]
    pub fn target(&self) -> &IntPoly {
        &self.target
    }

    /// The derivative `F'`.
    #[must_use]
    pub fn derivative(&self) -> &IntPoly {
        &self.derivative
    }

    /// The inverse of `F'` modulo `F`, verified at construction by `F'·inverse ≡ 1`.
    #[must_use]
    pub fn derivative_inverse(&self) -> &RatPoly {
        &self.derivative_inverse
    }

    /// The certified bounds `B_i` with `|z_i| ≤ B_i` for every root `z`, one per coefficient.
    #[must_use]
    pub fn coefficient_bounds(&self) -> &[BigInt] {
        &self.bounds
    }

    /// `B = max_i B_i`.
    #[must_use]
    pub fn bound(&self) -> &BigInt {
        &self.bound
    }

    /// The least integer modulus `M` with `M > 2B`, namely `2B + 1`.
    #[must_use]
    pub fn minimum_modulus(&self) -> &BigInt {
        &self.minimum_modulus
    }

    /// The first failing shape, bound, or square check on `z`, in that order.
    fn reject_candidate(&self, z: &IntPoly) -> Option<CandidateRejection> {
        if z.degree().is_some_and(|degree| degree >= self.context().degree()) {
            return Some(CandidateRejection::DegreeTooHigh);
        }
        // Coefficients beyond z's length are zero and trivially within any bound.
        let excess = z.coeffs.iter().zip(&self.bounds).position(|(zi, bi)| zi.abs() > *bi);
        if let Some(index) = excess {
            return Some(CandidateRejection::CoefficientOutOfBound { index });
        }
        let square = reduce_monic(&z.mul(z), self.context().monic_poly());
        (square != self.target).then_some(CandidateRejection::NotSquareRootOfTarget)
    }

    /// `β = z / (c^h F')` reduced modulo `F`, checked against `β² = γ`.
    fn recover_beta(&self, z: &IntPoly) -> Result<RatPoly, SqrtError> {
        let ctx = self.context();
        let modulus = ctx.monic_poly().to_rat_poly();
        let c_to_h = pow(ctx.leading_coeff().clone(), self.half);
        let unscale = BigRational::new(BigInt::one(), c_to_h);
        let beta = z.to_rat_poly().mul(&self.derivative_inverse).rem(&modulus).scale(&unscale);
        if beta.mul(&beta).rem(&modulus) != self.gamma {
            return Err(violation("recovered root does not square to the target element"));
        }
        Ok(beta)
    }

    /// `y = β(c·m) mod N`, or the verdict explaining why the coefficientwise map is unavailable.
    fn modular_image(&self, beta: &RatPoly) -> Result<BigInt, CandidateVerdict> {
        let ctx = self.context();
        let n = ctx.modulus();
        let mut residues = Vec::with_capacity(beta.coeffs.len());
        let mut blocked = false;
        for coeff in &beta.coeffs {
            let inverse = coeff.denom().extended_gcd(n);
            if inverse.gcd.is_one() {
                residues.push((coeff.numer() * &inverse.x).mod_floor(n));
                continue;
            }
            // `new` refuses divisor = N, which is exactly the unsupported case.
            match IncidentalFactor::new(n, inverse.gcd, IncidentalSource::MapDenominator) {
                Some(found) => return Err(CandidateVerdict::Incidental(found)),
                None => blocked = true,
            }
        }
        if blocked {
            return Err(CandidateVerdict::Unsupported(UnsupportedDomain::MapDenominatorNotUnit));
        }
        let theta = ctx.theta_image();
        Ok(residues.iter().rev().fold(BigInt::zero(), |acc, r| (acc * theta + r).mod_floor(n)))
    }
}

fn violation(what: &'static str) -> SqrtError {
    SqrtError::Invariant(InvariantViolation { what })
}

/// Remainder of `p` modulo the monic `monic`; exact over `Z` because the divisor is monic.
fn reduce_monic(p: &IntPoly, monic: &IntPoly) -> IntPoly {
    p.pseudo_div_rem(monic).1
}

/// `Π(c·a − b·t) mod F` over `Z`, reduced after every factor.
fn scaled_factor_product(ctx: &SqrtContext, pairs: &[(BigInt, BigInt)]) -> IntPoly {
    pairs.iter().fold(IntPoly::one(), |acc, (a, b)| {
        let factor = IntPoly::from_coeffs(vec![ctx.leading_coeff() * a, -b]);
        reduce_monic(&acc.mul(&factor), ctx.monic_poly())
    })
}

/// The bounds `B_i` and `B = max B_i` for the integral roots of `eta`, in exact integers.
///
/// `R = 1 + max_{i<d}|F_i|` bounds every root of `F` (Cauchy). Then `|η(θ_j)| ≤ E =
/// Σ|η_i| R^i`, so `|w_j| ≤ U = ⌈√E⌉`, and `B_i = d·U·Σ_{ℓ>i} |F_ℓ| R^(ℓ−1−i)` bounds the trace
/// expression for `z_i`. The inner sum satisfies `S_(i−1) = |F_i| + R·S_i`.
fn coefficient_bounds(monic: &IntPoly, eta: &IntPoly) -> (Vec<BigInt>, BigInt) {
    let degree = monic.coeffs.len() - 1;
    let root_radius = BigInt::one()
        + monic.coeffs[..degree].iter().fold(BigInt::zero(), |acc, c| acc.max(c.abs()));
    let powers = powers_of(&root_radius, degree);
    let eta_weight: BigInt = eta.coeffs.iter().zip(&powers).map(|(e, r)| e.abs() * r).sum();
    let scale = BigInt::from(degree) * ceil_sqrt(&eta_weight);

    let mut bounds = vec![BigInt::zero(); degree];
    let mut inner = BigInt::zero();
    for i in (0..degree).rev() {
        inner = monic.coeffs[i + 1].abs() + &root_radius * inner;
        bounds[i] = &scale * &inner;
    }
    let bound = bounds.iter().fold(BigInt::zero(), |acc, b| acc.max(b.clone()));
    (bounds, bound)
}

/// `⌈√e⌉` for `e >= 0`, exactly: the floor root, plus one unless it squares to `e`.
fn ceil_sqrt(e: &BigInt) -> BigInt {
    let floor = e.sqrt();
    if &floor * &floor < *e { floor + 1 } else { floor }
}

/// The inverse of `derivative` modulo `monic` over `Q` by the extended Euclidean algorithm,
/// returned only after `derivative · inverse ≡ 1 (mod monic)` holds exactly.
fn verified_inverse(derivative: &IntPoly, monic: &IntPoly) -> Result<RatPoly, SqrtError> {
    let modulus = monic.to_rat_poly();
    let derivative = derivative.to_rat_poly();
    // Invariant: r_i ≡ s_i · derivative (mod modulus).
    let (mut r0, mut r1) = (modulus.clone(), derivative.clone());
    let (mut s0, mut s1) = (RatPoly::zero(), RatPoly::one());
    while r1.degree().is_some() {
        let (quotient, remainder) = r0.div_rem(&r1);
        let next = s0.sub(&quotient.mul(&s1));
        (r0, r1, s0, s1) = (r1, remainder, s1, next);
    }
    let unit = match (r0.degree(), r0.leading_coeff()) {
        (Some(0), Some(unit)) => unit.clone(),
        _ => return Err(violation("derivative is not invertible modulo the monic polynomial")),
    };
    let inverse = s0.rem(&modulus).scale(&unit.recip());
    if inverse.mul(&derivative).rem(&modulus) != RatPoly::one() {
        return Err(violation("derivative inverse fails its identity"));
    }
    Ok(inverse)
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

    /// `2x² − 1` at `m = 2`, `N = 7`, relation `(3, −4)`, whose root is `z = 8 + 4θ`: see the
    /// integration KATs for the arithmetic.
    fn non_monic_target() -> NormalizedTarget {
        let ctx = SqrtContext::new(&poly(&[-1, 0, 2]), &big(2), &big(7)).unwrap();
        let factors = RelationFactors::new(&ctx, &[(big(3), big(-4))]).unwrap();
        match NormalizedTarget::new(factors, PrecisionLimits { max_modulus_bits: 64 }).unwrap() {
            TargetOutcome::Ready(target) => target,
            TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
        }
    }

    fn method() -> ReconstructionMethod {
        ReconstructionMethod::DirectNewton { prime: 3, exponent: 4 }
    }

    #[test]
    fn ceil_sqrt_is_exact_on_squares_and_neighbours() {
        // ⌈√e⌉ for small e: 0→0, 1→1, 2→2, 3→2, 4→2, 5→3, 9→3, 10→4, 16→4, 17→5.
        let expected = [(0, 0), (1, 1), (2, 2), (3, 2), (4, 2), (5, 3), (9, 3), (10, 4), (16, 4)];
        for (e, root) in expected {
            assert_eq!(ceil_sqrt(&big(e)), big(root), "e = {e}");
        }
        assert_eq!(ceil_sqrt(&big(17)), big(5));
        // 10^40 = (10^20)²; one more must round up, one less must not drop below 10^20.
        let root = BigInt::from(10u32).pow(20);
        let square = &root * &root;
        assert_eq!(ceil_sqrt(&square), root);
        assert_eq!(ceil_sqrt(&(&square + 1)), &root + 1);
        assert_eq!(ceil_sqrt(&(&square - 1)), root);
    }

    #[test]
    fn corrupted_derivative_inverse_is_an_invariant_error() {
        // With the inverse replaced by 1, β = z/c^h = 4 + 2θ squares to 24 + 16θ, not γ = 3 + 2θ.
        let mut target = non_monic_target();
        target.derivative_inverse = RatPoly::one();
        let outcome = target.verify_candidate(&poly(&[8, 4]), method());
        assert!(matches!(outcome, Err(SqrtError::Invariant(_))), "{outcome:?}");
    }

    #[test]
    fn corrupted_rational_product_is_an_invariant_error() {
        // y = 5 gives y² = 25 ≡ 4 (mod 7); with P replaced by 12 ≡ 5 the final congruence fails.
        let mut target = non_monic_target();
        target.rational_product = big(12);
        let outcome = target.verify_candidate(&poly(&[8, 4]), method());
        assert!(matches!(outcome, Err(SqrtError::Invariant(_))), "{outcome:?}");
    }

    #[test]
    fn non_invertible_derivative_is_an_invariant_error() {
        // (t − 1)² has derivative 2(t − 1), which shares the factor t − 1 with it, so the gcd is
        // not a unit and no inverse exists. A validated context never reaches this case.
        let repeated = poly(&[1, -2, 1]);
        let result = verified_inverse(&repeated.derivative(), &repeated);
        assert!(matches!(result, Err(SqrtError::Invariant(_))), "{result:?}");
    }
}
