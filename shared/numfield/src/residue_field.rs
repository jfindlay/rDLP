//! Residue fields `F_p[x]/(f)` for an irreducible `f`, with inversion and square roots.
//!
//! A [`ResidueField`] is validated at construction: the modulus polynomial must be nonzero, of
//! positive degree, and irreducible over a proved-prime `F_p` (decided exactly by
//! [`FpPoly::is_irreducible`]). Elements are plain [`FpPoly`] values of degree below the field
//! degree; every operation first reduces its operands modulo `f`, so an unreduced polynomial is an
//! accepted spelling of its class and no element ever borrows a field.
//!
//! The field has `q = p^d` elements. Cardinalities and exponents (`q`, `(q − 1)/2`, the odd part of
//! `q − 1`) are [`BigUint`], so they never overflow a machine word.
//!
//! Square roots follow Tonelli–Shanks: an Euler-criterion test decides squareness exactly, and a
//! nonresidue is searched only when the algorithm actually needs one. The search enumerates
//! candidate elements in a fixed order under an explicit budget ([`NonresidueSearchLimits`]).
//! Running out of budget is the resource outcome [`SqrtOutcome::ResourceExhausted`]; it never
//! claims anything about the input. Characteristic 2 is handled directly, because the Frobenius map
//! makes every element a square with a unique root.

use std::fmt;

use num_bigint::BigUint;

use crate::fp_poly::{FpPoly, FpPolyError, PrimeModulus};
use crate::poly::IntPoly;

// ─── errors ──────────────────────────────────────────────────────────────────

/// Error returned by checked residue-field construction and arithmetic.
///
/// Each variant is invalid input or a violated internal invariant; none is a mathematical
/// certificate about a field element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueFieldError {
    /// An underlying prime-field error: a non-prime characteristic, or an element over a different
    /// characteristic than the field.
    Poly(FpPolyError),
    /// The modulus polynomial is zero.
    ZeroModulus,
    /// The modulus polynomial is a nonzero constant, which defines the zero ring.
    ConstantModulus,
    /// The modulus polynomial is not irreducible, so the quotient is not a field.
    ReducibleModulus,
    /// The zero element has no inverse.
    DivisionByZero,
    /// An internal consistency check failed. A proved-irreducible modulus makes this unreachable;
    /// it is reported instead of returning an unverified result.
    InvariantViolation,
}

impl From<FpPolyError> for ResidueFieldError {
    fn from(e: FpPolyError) -> Self {
        Self::Poly(e)
    }
}

impl fmt::Display for ResidueFieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Poly(e) => write!(f, "{e}"),
            Self::ZeroModulus => f.write_str("the modulus polynomial is zero"),
            Self::ConstantModulus => f.write_str("the modulus polynomial is a nonzero constant"),
            Self::ReducibleModulus => f.write_str("the modulus polynomial is not irreducible"),
            Self::DivisionByZero => f.write_str("the zero element has no inverse"),
            Self::InvariantViolation => f.write_str("residue-field invariant violated"),
        }
    }
}

// thiserror is not a numfield dependency, so the impl is by hand.
impl std::error::Error for ResidueFieldError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Poly(e) => Some(e),
            _ => None,
        }
    }
}

// ─── search budget and outcomes ──────────────────────────────────────────────

/// Candidate budget for the nonresidue search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonresidueSearchLimits {
    /// Maximum number of candidate elements tested; the search never exceeds it.
    pub max_candidates: u64,
}

/// Outcome of a bounded nonresidue search.
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub enum NonresidueSearch {
    /// A verified nonsquare of the field.
    Found {
        /// The first nonsquare in the fixed candidate order.
        element: FpPoly,
        /// Candidates tested, including the one found.
        candidates_tried: u64,
    },
    /// Characteristic 2: every element is a square, so there is nothing to find. A certificate.
    Absent,
    /// The budget ran out first. Resource exhaustion, never a proof about the field.
    BudgetExhausted {
        /// Candidates tested; equals the budget.
        candidates_tried: u64,
    },
}

/// Outcome of a square-root computation.
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub enum SqrtOutcome {
    /// A verified root `r` with `r² = a`. The other root is `−r` (equal to `r` for zero and in
    /// characteristic 2). Which of the two is returned is deterministic but unspecified.
    Root(FpPoly),
    /// `a` is a proved nonsquare: the Euler criterion gave `a^((q−1)/2) = −1`.
    NotSquare,
    /// Tonelli–Shanks needed a nonresidue and the budget ran out before one was found. Resource
    /// exhaustion, never evidence about `a`.
    ResourceExhausted {
        /// Candidates tested; equals the budget.
        candidates_tried: u64,
    },
}

// ─── ResidueField ────────────────────────────────────────────────────────────

/// The field `F_p[x]/(f)` for an irreducible `f`, with `q = p^d` elements.
///
/// Owned and immutable: elements are [`FpPoly`] values that do not borrow the field.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueField {
    /// The monic irreducible modulus polynomial.
    modulus: FpPoly,
    characteristic: PrimeModulus,
    degree: usize,
    /// `q = p^d`.
    order: BigUint,
}

impl ResidueField {
    /// Build `F_p[x]/(f)`, scaling `f` to be monic and proving it irreducible.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::ZeroModulus`] for the zero polynomial,
    /// [`ResidueFieldError::ConstantModulus`] for a nonzero constant, and
    /// [`ResidueFieldError::ReducibleModulus`] when `f` has a nontrivial factorisation (including
    /// repeated factors).
    pub fn new(f: &FpPoly) -> Result<Self, ResidueFieldError> {
        let degree = match f.degree() {
            None => return Err(ResidueFieldError::ZeroModulus),
            Some(0) => return Err(ResidueFieldError::ConstantModulus),
            Some(n) => n,
        };
        if !f.is_irreducible() {
            return Err(ResidueFieldError::ReducibleModulus);
        }
        let characteristic = f.modulus();
        Ok(Self {
            modulus: f.to_monic()?,
            characteristic,
            degree,
            order: num_traits::pow(characteristic.to_biguint(), degree),
        })
    }

    /// Build `F_p[x]/(f)` from an integer polynomial reduced modulo a candidate prime `p`.
    ///
    /// The field degree is that of the reduction, which is lower than `deg f` when `p` divides the
    /// leading coefficient.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::Poly`] carrying [`FpPolyError::NotPrime`] if `p` is not prime; otherwise
    /// as [`ResidueField::new`].
    pub fn from_int_poly(f: &IntPoly, p: u64) -> Result<Self, ResidueFieldError> {
        Self::new(&FpPoly::from_int_poly(f, PrimeModulus::new(p)?))
    }

    /// The characteristic `p`.
    #[must_use]
    pub fn characteristic(&self) -> PrimeModulus {
        self.characteristic
    }

    /// The extension degree `d = deg f`.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The monic irreducible modulus polynomial `f`.
    #[must_use]
    pub fn modulus_poly(&self) -> &FpPoly {
        &self.modulus
    }

    /// The number of elements `q = p^d`.
    #[must_use]
    pub fn order(&self) -> &BigUint {
        &self.order
    }

    /// The additive identity.
    #[must_use]
    pub fn zero(&self) -> FpPoly {
        FpPoly::zero(self.characteristic)
    }

    /// The multiplicative identity.
    #[must_use]
    pub fn one(&self) -> FpPoly {
        FpPoly::one(self.characteristic)
    }

    /// The canonical representative of `a`: its remainder modulo `f`, of degree below `d`.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::Poly`] carrying [`FpPolyError::ModulusMismatch`] if `a` is over a
    /// different characteristic.
    pub fn element(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        Ok(a.rem(&self.modulus)?)
    }

    /// `a + b`.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`], for either operand.
    pub fn add(&self, a: &FpPoly, b: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        Ok(self.element(a)?.add(&self.element(b)?)?)
    }

    /// `a − b`.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`], for either operand.
    pub fn sub(&self, a: &FpPoly, b: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        Ok(self.element(a)?.sub(&self.element(b)?)?)
    }

    /// `−a`.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`].
    pub fn neg(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        Ok(self.element(a)?.neg())
    }

    /// `a · b`.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`], for either operand.
    pub fn mul(&self, a: &FpPoly, b: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        Ok(a.mul_mod(b, &self.modulus)?)
    }

    /// `a^exp`, with `a^0 = 1`. The exponent is arbitrary precision.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`].
    pub fn pow(&self, a: &FpPoly, exp: &BigUint) -> Result<FpPoly, ResidueFieldError> {
        Ok(a.pow_mod(exp, &self.modulus)?)
    }

    /// The multiplicative inverse of `a`, by the extended Euclidean algorithm against `f`.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::DivisionByZero`] if `a` is zero in the field;
    /// [`ResidueFieldError::InvariantViolation`] if the Bézout identity does not yield a unit
    /// (impossible for an irreducible `f`); otherwise as [`ResidueField::element`].
    pub fn inv(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        let a = self.element(a)?;
        if a.is_zero() {
            return Err(ResidueFieldError::DivisionByZero);
        }
        let bezout = a.extended_gcd(&self.modulus)?;
        if bezout.gcd.degree() != Some(0) {
            return Err(ResidueFieldError::InvariantViolation);
        }
        self.element(&bezout.s)
    }

    /// Whether `a` is a square in the field. Zero is a square (of zero).
    ///
    /// In characteristic 2 every element is a square. Otherwise the Euler criterion decides it
    /// exactly: a nonzero `a` is a square iff `a^((q−1)/2) = 1`.
    ///
    /// # Errors
    ///
    /// As [`ResidueField::element`].
    pub fn is_square(&self, a: &FpPoly) -> Result<bool, ResidueFieldError> {
        let a = self.element(a)?;
        if a.is_zero() || self.is_even_characteristic() {
            return Ok(true);
        }
        Ok(self.euler(&a)? == self.one())
    }

    /// Search for a nonsquare, testing candidates in a fixed order within a budget.
    ///
    /// The search is deterministic, uses no randomness, and tests at most `limits.max_candidates`
    /// candidates. In characteristic 2 there is no nonsquare and the answer is
    /// [`NonresidueSearch::Absent`] without testing anything.
    ///
    /// Squareness in `F_q` equals squareness of the norm in `F_p`, and a scalar `c ∈ F_p*` is a
    /// square in `F_q` iff `d` is even or `c` is a square in `F_p`. So the candidate order depends
    /// on the parity of `d`, and each list is guaranteed to contain a nonsquare:
    ///
    /// - odd `d`: the constants `1, 2, …, p−1`; a nonresidue of `F_p` stays a nonresidue in `F_q`.
    /// - even `d`: the monic polynomials of degree `1 … d−1`, by degree and then by the base-`p`
    ///   value of the lower coefficients (constant term least significant). Constants are all
    ///   squares here, and a nonsquare `a = c·m` with monic `m` forces `m` to be a nonsquare too.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::InvariantViolation`] if the candidates are exhausted without a
    /// nonsquare in an odd-characteristic field (impossible: `(q−1)/2` nonsquares exist).
    pub fn find_nonresidue(
        &self,
        limits: NonresidueSearchLimits,
    ) -> Result<NonresidueSearch, ResidueFieldError> {
        if self.is_even_characteristic() {
            return Ok(NonresidueSearch::Absent);
        }
        for k in 1..=limits.max_candidates {
            let Some(candidate) = self.candidate(k) else {
                return Err(ResidueFieldError::InvariantViolation);
            };
            if self.euler(&candidate)? != self.one() {
                return Ok(NonresidueSearch::Found { element: candidate, candidates_tried: k });
            }
        }
        Ok(NonresidueSearch::BudgetExhausted { candidates_tried: limits.max_candidates })
    }

    /// A square root of `a` by Tonelli–Shanks, verified by squaring before it is returned.
    ///
    /// Zero has root zero. In characteristic 2 the root is `a^(q/2)`. Otherwise squareness is
    /// decided first by the Euler criterion ([`SqrtOutcome::NotSquare`] is a proof), and a
    /// nonresidue is searched under `limits` only if the algorithm needs one: when the odd part
    /// `t` of `q − 1` satisfies `a^t = 1` (always the case for `q ≡ 3 mod 4`) no search occurs, so
    /// even a zero budget succeeds there.
    ///
    /// # Errors
    ///
    /// [`ResidueFieldError::InvariantViolation`] if an internal consistency check fails
    /// (unreachable for a validated field); otherwise as [`ResidueField::element`].
    pub fn sqrt(
        &self,
        a: &FpPoly,
        limits: NonresidueSearchLimits,
    ) -> Result<SqrtOutcome, ResidueFieldError> {
        let a = self.element(a)?;
        if a.is_zero() {
            return Ok(SqrtOutcome::Root(a));
        }
        let outcome = if self.is_even_characteristic() {
            // x ↦ x² is a bijection (Frobenius), and (a^(q/2))² = a^q = a.
            SqrtOutcome::Root(self.pow(&a, &(&self.order >> 1u32))?)
        } else if self.euler(&a)? == self.one() {
            self.tonelli_shanks(&a, limits)?
        } else {
            return Ok(SqrtOutcome::NotSquare);
        };
        match outcome {
            SqrtOutcome::Root(r) if self.mul(&r, &r)? != a => {
                Err(ResidueFieldError::InvariantViolation)
            }
            other => Ok(other),
        }
    }

    fn is_even_characteristic(&self) -> bool {
        self.characteristic.get() == 2
    }

    /// `a^((q−1)/2)` for odd `q`: 1 for a nonzero square, `−1` for a nonsquare.
    fn euler(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError> {
        self.pow(a, &((&self.order - 1u32) >> 1u32))
    }

    /// The `k`-th candidate nonsquare (`k ≥ 1`), or `None` once the candidate list is spent.
    ///
    /// For odd `d` the list is the constants `1 … p−1`; for even `d` it is the monic polynomials of
    /// degree `1 … d−1`, by degree and then by the base-`p` value of the lower coefficients. See
    /// [`ResidueField::find_nonresidue`] for why each list contains a nonsquare.
    fn candidate(&self, k: u64) -> Option<FpPoly> {
        let p = self.characteristic.get();
        let m = self.characteristic;
        if !self.degree.is_multiple_of(2) {
            return (k < p).then(|| FpPoly::constant(m, k));
        }
        let mut idx = k - 1;
        for j in 1..self.degree {
            // A block of p^j candidates shares the leading term x^j; an overflowing block is
            // larger than any index.
            let block = u32::try_from(j).ok().and_then(|e| p.checked_pow(e));
            if block.is_none_or(|b| idx < b) {
                let mut coeffs = Vec::with_capacity(j + 1);
                for _ in 0..j {
                    coeffs.push(idx % p);
                    idx /= p;
                }
                coeffs.push(1);
                return Some(FpPoly::from_coeffs(m, &coeffs));
            }
            idx -= block.unwrap_or(0);
        }
        None
    }

    /// Tonelli–Shanks for a nonzero square `a` in an odd-order field.
    ///
    /// Write `q − 1 = 2^s·t` with `t` odd. The invariants are `r² = a·tt`, `tt` of order dividing
    /// `2^(m−1)`, and `c` of order exactly `2^m`; each pass strictly lowers `m`.
    fn tonelli_shanks(
        &self,
        a: &FpPoly,
        limits: NonresidueSearchLimits,
    ) -> Result<SqrtOutcome, ResidueFieldError> {
        let q1 = &self.order - 1u32;
        let s = q1.trailing_zeros().ok_or(ResidueFieldError::InvariantViolation)?;
        let t = &q1 >> s;
        let one = self.one();
        let mut r = self.pow(a, &((&t + 1u32) >> 1u32))?;
        let mut tt = self.pow(a, &t)?;
        let mut m = s;
        let mut c: Option<FpPoly> = None;
        while tt != one {
            let c_cur = match c.take() {
                Some(c) => c,
                None => match self.find_nonresidue(limits)? {
                    NonresidueSearch::Found { element, .. } => self.pow(&element, &t)?,
                    NonresidueSearch::BudgetExhausted { candidates_tried } => {
                        return Ok(SqrtOutcome::ResourceExhausted { candidates_tried });
                    }
                    NonresidueSearch::Absent => return Err(ResidueFieldError::InvariantViolation),
                },
            };
            // Least i with tt^(2^i) = 1; i < m because tt has order dividing 2^(m−1).
            let mut i = 0u64;
            let mut x = tt.clone();
            while x != one {
                x = self.mul(&x, &x)?;
                i += 1;
                if i >= m {
                    return Err(ResidueFieldError::InvariantViolation);
                }
            }
            let mut b = c_cur;
            for _ in 0..(m - i - 1) {
                b = self.mul(&b, &b)?;
            }
            r = self.mul(&r, &b)?;
            let b2 = self.mul(&b, &b)?;
            tt = self.mul(&tt, &b2)?;
            c = Some(b2);
            m = i;
        }
        Ok(SqrtOutcome::Root(r))
    }
}
