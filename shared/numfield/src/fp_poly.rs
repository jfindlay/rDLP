//! Checked polynomial arithmetic over a prime field `F_p`, for any prime p < 2^64.
//!
//! [`PrimeModulus`] is a prime proved at construction. [`FpPoly`] is a polynomial over `F_p` in
//! canonical form: coefficients least-significant first, each in `[0, p)`, trailing zeros trimmed,
//! zero polynomial = empty slice. The fields are private, so every value is canonical and equality
//! and hashing include the modulus.
//!
//! Scalar sums and products widen to `u128` and reduce immediately; polynomial products reduce
//! after every multiply-add. Operations on mixed moduli and zero divisors return [`FpPolyError`];
//! no input panics. Root enumeration is exhaustive and resource-accounted ([`RootEnumeration`]):
//! a truncated search never claims anything about the roots it did not reach.
//!
//! Primality is decided deterministically for all of `u64` by [`is_prime_u64`].

use std::collections::BTreeMap;
use std::fmt;

use num_bigint::{BigInt, BigUint};
use num_traits::ToPrimitive;

use crate::poly::IntPoly;

// ─── errors ──────────────────────────────────────────────────────────────────

/// Error returned by checked finite-field operations.
///
/// `NotPrime` is invalid input; `ModulusMismatch` is a caller invariant violation;
/// `DivisionByZero` and `ZeroPolynomial` are invalid arguments. None is a mathematical certificate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FpPolyError {
    /// The candidate modulus is not a prime.
    NotPrime {
        /// The rejected candidate.
        candidate: u64,
    },
    /// A binary operation received operands over different moduli.
    ModulusMismatch {
        /// Modulus of the left operand.
        left: u64,
        /// Modulus of the right operand.
        right: u64,
    },
    /// A zero scalar was inverted or a zero polynomial was used as divisor or modulus.
    DivisionByZero,
    /// The operation is undefined for the zero polynomial.
    ZeroPolynomial,
}

impl fmt::Display for FpPolyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPrime { candidate } => write!(f, "{candidate} is not prime"),
            Self::ModulusMismatch { left, right } => {
                write!(f, "operands are over different moduli ({left} and {right})")
            }
            Self::DivisionByZero => f.write_str("division by zero"),
            Self::ZeroPolynomial => f.write_str("operation undefined for the zero polynomial"),
        }
    }
}

// thiserror is not a numfield dependency, so the impl is by hand.
impl std::error::Error for FpPolyError {}

// ─── deterministic u64 primality ─────────────────────────────────────────────

/// The first twelve primes: trial-division set and strong-probable-prime bases.
const PRIME_BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Narrow a `u128` known to be below 2^64.
#[allow(clippy::cast_possible_truncation)] // callers reduce modulo a u64 value first
fn narrow(x: u128) -> u64 {
    x as u64
}

fn mul_mod_u64(a: u64, b: u64, n: u64) -> u64 {
    narrow(u128::from(a) * u128::from(b) % u128::from(n))
}

fn pow_mod_u64(mut base: u64, mut exp: u64, n: u64) -> u64 {
    let mut acc = 1 % n;
    base %= n;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod_u64(acc, base, n);
        }
        base = mul_mod_u64(base, base, n);
        exp >>= 1;
    }
    acc
}

/// Strong-probable-prime test of odd `n = d·2^s + 1` to base `a`, with `1 < a < n`.
fn is_strong_probable_prime(n: u64, a: u64, d: u64, s: u32) -> bool {
    let mut x = pow_mod_u64(a, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 1..s {
        x = mul_mod_u64(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

/// Decide whether `n` is prime. Deterministic and exact for every `u64`.
///
/// Returns false for `n < 2`. Trial division by the first twelve primes (2..=37) settles every `n`
/// divisible by one of them (`n` equal to one of them is prime). Otherwise `n > 37` is coprime to
/// all of them, and `n` is prime iff it is a strong probable prime to all twelve bases, using
/// `u128` products. This is a proof rather than probable primality because the least composite
/// strong pseudoprime to all of the first twelve prime bases is
/// ψ₁₂ = 318665857834031151167461 > 2^64 (Sorenson and Webster, "Strong pseudoprimes to twelve
/// prime bases", Math. Comp. 86 (2017)), so no composite below 2^64 passes. Fewer prime bases do
/// not suffice: ψ₉ = ψ₁₀ = ψ₁₁ = 3825123056546413051 < 2^64 (Jiang and Deng, Math. Comp. 83
/// (2014)). The same table lives in `shared/numth/src/prime.rs`; it is duplicated here because
/// depending on that crate would add a numfield dependency.
#[must_use]
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &PRIME_BASES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    PRIME_BASES.iter().all(|&a| is_strong_probable_prime(n, a, d, s))
}

// ─── PrimeModulus ────────────────────────────────────────────────────────────

/// A prime modulus p < 2^64, proved prime at construction by [`is_prime_u64`].
///
/// All scalar methods accept any `u64` (or integer) operand and return the least nonnegative
/// residue; intermediate products widen to `u128`, so results are exact for every prime p < 2^64.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PrimeModulus {
    p: u64,
}

impl PrimeModulus {
    /// Validate `p` as a prime modulus.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::NotPrime`] if `p` is 0, 1, or composite.
    pub fn new(p: u64) -> Result<Self, FpPolyError> {
        if is_prime_u64(p) { Ok(Self { p }) } else { Err(FpPolyError::NotPrime { candidate: p }) }
    }

    /// The modulus as a machine integer.
    #[must_use]
    pub fn get(self) -> u64 {
        self.p
    }

    /// The modulus as an arbitrary-precision integer.
    #[must_use]
    pub fn to_biguint(self) -> BigUint {
        BigUint::from(self.p)
    }

    /// Least nonnegative residue of `a`.
    #[must_use]
    pub fn reduce(self, a: u64) -> u64 {
        a % self.p
    }

    /// Least nonnegative residue of a signed integer.
    #[must_use]
    pub fn reduce_i64(self, a: i64) -> u64 {
        let r = i128::from(a).rem_euclid(i128::from(self.p));
        // rem_euclid gives 0 <= r < p < 2^64.
        narrow(r.unsigned_abs())
    }

    /// Least nonnegative residue of an arbitrary-precision integer.
    #[must_use]
    pub fn reduce_bigint(self, a: &BigInt) -> u64 {
        let r = (a.magnitude() % self.p).to_u64().unwrap_or(0); // remainder < p fits u64
        if a.sign() == num_bigint::Sign::Minus && r != 0 { self.p - r } else { r }
    }

    /// `a + b mod p`.
    #[must_use]
    pub fn add(self, a: u64, b: u64) -> u64 {
        narrow((u128::from(a) + u128::from(b)) % u128::from(self.p))
    }

    /// `a − b mod p`.
    #[must_use]
    pub fn sub(self, a: u64, b: u64) -> u64 {
        let p = u128::from(self.p);
        narrow((u128::from(a) % p + p - u128::from(b) % p) % p)
    }

    /// `−a mod p`.
    #[must_use]
    pub fn neg(self, a: u64) -> u64 {
        self.sub(0, a)
    }

    /// `a · b mod p`.
    #[must_use]
    pub fn mul(self, a: u64, b: u64) -> u64 {
        mul_mod_u64(a, b, self.p)
    }

    /// `a^e mod p`, with `a^0 = 1`.
    #[must_use]
    pub fn pow(self, a: u64, e: u64) -> u64 {
        pow_mod_u64(a, e, self.p)
    }

    /// Multiplicative inverse of `a` by Fermat's little theorem.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::DivisionByZero`] if `a ≡ 0 mod p`.
    pub fn inv(self, a: u64) -> Result<u64, FpPolyError> {
        if self.reduce(a) == 0 {
            return Err(FpPolyError::DivisionByZero);
        }
        Ok(self.pow(a, self.p - 2))
    }
}

// ─── result types ────────────────────────────────────────────────────────────

/// Bézout data: `s·a + t·b = gcd`, with `gcd` monic (or zero when `a = b = 0`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtendedGcd {
    /// The canonical monic greatest common divisor.
    pub gcd: FpPoly,
    /// Cofactor of the left operand.
    pub s: FpPoly,
    /// Cofactor of the right operand.
    pub t: FpPoly,
}

/// Evaluation budget for [`FpPoly::roots`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RootSearchLimits {
    /// Maximum number of candidate evaluations; the search never exceeds it.
    pub max_evaluations: u64,
}

/// Outcome of an exhaustive, budgeted root search in `F_p`.
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub enum RootEnumeration {
    /// Every root was found: distinct and ascending.
    Complete {
        /// All roots in `F_p`.
        roots: Vec<u64>,
        /// Candidates evaluated.
        evaluations: u64,
    },
    /// The budget ran out first. Resource exhaustion, never a proof of absence.
    Truncated {
        /// Exactly the roots in `[0, evaluations)`.
        found: Vec<u64>,
        /// The exact number of distinct roots in `F_p` (greater than `found.len()`).
        distinct_root_count: usize,
        /// Candidates evaluated; equals the budget.
        evaluations: u64,
    },
}

// ─── FpPoly ──────────────────────────────────────────────────────────────────

/// A polynomial over `F_p` in canonical form (see the module documentation).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FpPoly {
    modulus: PrimeModulus,
    coeffs: Vec<u64>,
}

fn trim(v: &mut Vec<u64>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

/// Prime divisors of `n` in ascending order, by trial division.
fn prime_divisors(mut n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut r = 2;
    while r * r <= n {
        if n.is_multiple_of(r) {
            out.push(r);
            while n.is_multiple_of(r) {
                n /= r;
            }
        }
        r += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

impl FpPoly {
    /// Wrap coefficients already reduced into `[0, p)`, trimming trailing zeros.
    fn from_reduced(modulus: PrimeModulus, mut coeffs: Vec<u64>) -> Self {
        trim(&mut coeffs);
        Self { modulus, coeffs }
    }

    fn same_modulus(&self, other: &Self) -> Result<(), FpPolyError> {
        if self.modulus == other.modulus {
            Ok(())
        } else {
            Err(FpPolyError::ModulusMismatch {
                left: self.modulus.get(),
                right: other.modulus.get(),
            })
        }
    }

    /// The zero polynomial.
    #[must_use]
    pub fn zero(m: PrimeModulus) -> Self {
        Self { modulus: m, coeffs: Vec::new() }
    }

    /// The constant polynomial 1.
    #[must_use]
    pub fn one(m: PrimeModulus) -> Self {
        Self::constant(m, 1)
    }

    /// The constant polynomial `c mod p`.
    #[must_use]
    pub fn constant(m: PrimeModulus, c: u64) -> Self {
        Self::from_reduced(m, vec![m.reduce(c)])
    }

    /// The monomial `c·x^degree`, with `c` reduced mod p.
    #[must_use]
    pub fn monomial(m: PrimeModulus, c: u64, degree: usize) -> Self {
        let mut coeffs = vec![0; degree + 1];
        coeffs[degree] = m.reduce(c);
        Self::from_reduced(m, coeffs)
    }

    /// Build from least-significant-first coefficients, reducing each mod p and trimming.
    #[must_use]
    pub fn from_coeffs(m: PrimeModulus, coeffs: &[u64]) -> Self {
        Self::from_reduced(m, coeffs.iter().map(|&c| m.reduce(c)).collect())
    }

    /// Reduce an integer polynomial mod p (least nonnegative residues).
    ///
    /// Reducing a non-monic polynomial may lower its degree; callers needing the degree preserved
    /// check [`FpPoly::degree`].
    #[must_use]
    pub fn from_int_poly(f: &IntPoly, m: PrimeModulus) -> Self {
        Self::from_reduced(m, f.coeffs.iter().map(|c| m.reduce_bigint(c)).collect())
    }

    /// The integer polynomial with representatives in `[0, p)`.
    #[must_use]
    pub fn to_int_poly(&self) -> IntPoly {
        IntPoly::from_coeffs(self.coeffs.iter().map(|&c| BigInt::from(c)).collect())
    }

    /// The coefficient field's modulus.
    #[must_use]
    pub fn modulus(&self) -> PrimeModulus {
        self.modulus
    }

    /// Canonical coefficients, least-significant first.
    #[must_use]
    pub fn coeffs(&self) -> &[u64] {
        &self.coeffs
    }

    /// Degree, or `None` for the zero polynomial.
    #[must_use]
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// Leading coefficient, or `None` for the zero polynomial.
    #[must_use]
    pub fn leading_coeff(&self) -> Option<u64> {
        self.coeffs.last().copied()
    }

    /// True iff this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// True iff the leading coefficient is 1.
    #[must_use]
    pub fn is_monic(&self) -> bool {
        self.leading_coeff() == Some(1)
    }

    /// Evaluate at `x mod p` by Horner's rule.
    #[must_use]
    pub fn eval(&self, x: u64) -> u64 {
        let m = self.modulus;
        let x = m.reduce(x);
        self.coeffs.iter().rev().fold(0, |acc, &c| m.add(m.mul(acc, x), c))
    }

    /// Formal derivative; coefficient `i·c_i` is taken mod p, so `d/dx x^p = 0`.
    #[must_use]
    pub fn derivative(&self) -> Self {
        let m = self.modulus;
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| m.mul(c, m.reduce(i as u64)))
            .collect();
        Self::from_reduced(m, coeffs)
    }

    /// Additive inverse.
    #[must_use]
    pub fn neg(&self) -> Self {
        let m = self.modulus;
        Self::from_reduced(m, self.coeffs.iter().map(|&c| m.neg(c)).collect())
    }

    /// Multiply by the scalar `c mod p`.
    #[must_use]
    pub fn scale(&self, c: u64) -> Self {
        let m = self.modulus;
        Self::from_reduced(m, self.coeffs.iter().map(|&a| m.mul(a, c)).collect())
    }

    /// Sum of two polynomials over the same modulus.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ.
    pub fn add(&self, rhs: &Self) -> Result<Self, FpPolyError> {
        self.same_modulus(rhs)?;
        let m = self.modulus;
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let at = |v: &[u64], i: usize| v.get(i).copied().unwrap_or(0);
        let coeffs = (0..len).map(|i| m.add(at(&self.coeffs, i), at(&rhs.coeffs, i))).collect();
        Ok(Self::from_reduced(m, coeffs))
    }

    /// Difference of two polynomials over the same modulus.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ.
    pub fn sub(&self, rhs: &Self) -> Result<Self, FpPolyError> {
        self.same_modulus(rhs)?;
        let m = self.modulus;
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let at = |v: &[u64], i: usize| v.get(i).copied().unwrap_or(0);
        let coeffs = (0..len).map(|i| m.sub(at(&self.coeffs, i), at(&rhs.coeffs, i))).collect();
        Ok(Self::from_reduced(m, coeffs))
    }

    /// Product of two polynomials over the same modulus (schoolbook; reduces after every step).
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ.
    pub fn mul(&self, rhs: &Self) -> Result<Self, FpPolyError> {
        self.same_modulus(rhs)?;
        let m = self.modulus;
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero(m));
        }
        let mut out = vec![0; self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in rhs.coeffs.iter().enumerate() {
                out[i + j] = m.add(out[i + j], m.mul(a, b));
            }
        }
        Ok(Self::from_reduced(m, out))
    }

    /// Quotient and remainder: `self = q·divisor + r` with `r = 0` or `deg r < deg divisor`.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ; [`FpPolyError::DivisionByZero`] if
    /// `divisor` is zero.
    pub fn div_rem(&self, divisor: &Self) -> Result<(Self, Self), FpPolyError> {
        self.same_modulus(divisor)?;
        let m = self.modulus;
        let lead_inv = m.inv(divisor.leading_coeff().ok_or(FpPolyError::DivisionByZero)?)?;
        let dlen = divisor.coeffs.len();
        let mut rem = self.coeffs.clone();
        let mut quo = vec![0; rem.len().saturating_sub(dlen - 1)];
        while rem.len() >= dlen {
            let shift = rem.len() - dlen;
            let factor = m.mul(rem[rem.len() - 1], lead_inv);
            quo[shift] = factor;
            for (j, &dc) in divisor.coeffs.iter().enumerate() {
                rem[shift + j] = m.sub(rem[shift + j], m.mul(factor, dc));
            }
            // The leading term cancelled exactly, so trimming strictly shortens `rem`.
            trim(&mut rem);
        }
        Ok((Self::from_reduced(m, quo), Self::from_reduced(m, rem)))
    }

    /// Remainder of division by `divisor`.
    ///
    /// # Errors
    ///
    /// As [`FpPoly::div_rem`].
    pub fn rem(&self, divisor: &Self) -> Result<Self, FpPolyError> {
        self.div_rem(divisor).map(|(_, r)| r)
    }

    /// Scale to leading coefficient 1.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ZeroPolynomial`] for the zero polynomial.
    pub fn to_monic(&self) -> Result<Self, FpPolyError> {
        let lead = self.leading_coeff().ok_or(FpPolyError::ZeroPolynomial)?;
        Ok(self.scale(self.modulus.inv(lead)?))
    }

    /// Monic greatest common divisor; `gcd(0, 0) = 0`.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ.
    pub fn gcd(&self, rhs: &Self) -> Result<Self, FpPolyError> {
        self.same_modulus(rhs)?;
        let (mut a, mut b) = (self.clone(), rhs.clone());
        while !b.is_zero() {
            let r = a.rem(&b)?;
            a = b;
            b = r;
        }
        if a.is_zero() { Ok(a) } else { a.to_monic() }
    }

    /// Canonical gcd with Bézout cofactors satisfying `s·self + t·rhs = gcd` exactly.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] if the moduli differ.
    pub fn extended_gcd(&self, rhs: &Self) -> Result<ExtendedGcd, FpPolyError> {
        self.same_modulus(rhs)?;
        let m = self.modulus;
        let (mut r0, mut r1) = (self.clone(), rhs.clone());
        let (mut s0, mut s1) = (Self::one(m), Self::zero(m));
        let (mut t0, mut t1) = (Self::zero(m), Self::one(m));
        while !r1.is_zero() {
            let (q, r2) = r0.div_rem(&r1)?;
            let s2 = s0.sub(&q.mul(&s1)?)?;
            let t2 = t0.sub(&q.mul(&t1)?)?;
            (r0, r1) = (r1, r2);
            (s0, s1) = (s1, s2);
            (t0, t1) = (t1, t2);
        }
        let Some(lead) = r0.leading_coeff() else {
            return Ok(ExtendedGcd { gcd: r0, s: Self::zero(m), t: Self::zero(m) });
        };
        let inv = m.inv(lead)?;
        Ok(ExtendedGcd { gcd: r0.scale(inv), s: s0.scale(inv), t: t0.scale(inv) })
    }

    /// `self · rhs` reduced modulo `modulus`.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] on mixed moduli; [`FpPolyError::DivisionByZero`] if
    /// `modulus` is zero.
    pub fn mul_mod(&self, rhs: &Self, modulus: &Self) -> Result<Self, FpPolyError> {
        self.mul(rhs)?.rem(modulus)
    }

    /// `self^exp` reduced modulo `modulus`, by left-to-right square-and-multiply.
    ///
    /// `exp` is arbitrary precision, so exponents such as `p^d` or `(p^d − 1)/2` do not overflow;
    /// `self` need not be reduced. `exp = 0` gives `1 rem modulus`, and a nonzero constant
    /// `modulus` gives zero.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ModulusMismatch`] on mixed moduli; [`FpPolyError::DivisionByZero`] if
    /// `modulus` is zero.
    pub fn pow_mod(&self, exp: &BigUint, modulus: &Self) -> Result<Self, FpPolyError> {
        self.pow_mod_bits((0..exp.bits()).rev().map(|i| exp.bit(i)), modulus)
    }

    /// [`FpPoly::pow_mod`] for a machine-word exponent (used for Frobenius powers).
    fn pow_mod_u64(&self, exp: u64, modulus: &Self) -> Result<Self, FpPolyError> {
        let bits = u64::BITS - exp.leading_zeros();
        self.pow_mod_bits((0..bits).rev().map(|i| (exp >> i) & 1 == 1), modulus)
    }

    fn pow_mod_bits(
        &self,
        bits: impl Iterator<Item = bool>,
        modulus: &Self,
    ) -> Result<Self, FpPolyError> {
        self.same_modulus(modulus)?;
        if modulus.is_zero() {
            return Err(FpPolyError::DivisionByZero);
        }
        let base = self.rem(modulus)?;
        let mut acc = Self::one(self.modulus).rem(modulus)?;
        for bit in bits {
            acc = acc.mul_mod(&acc, modulus)?;
            if bit {
                acc = acc.mul_mod(&base, modulus)?;
            }
        }
        Ok(acc)
    }

    /// True iff the polynomial is nonzero and has no repeated irreducible factor.
    ///
    /// Zero is not squarefree; a nonzero constant is. Otherwise `gcd(f, f') = 1`. A positive-degree
    /// p-th power has `f' = 0`, so the gcd is `f` itself and it is correctly rejected.
    #[must_use]
    pub fn is_squarefree(&self) -> bool {
        if self.is_zero() {
            return false;
        }
        matches!(self.gcd(&self.derivative()), Ok(g) if g.degree() == Some(0))
    }

    /// Square-free decomposition: `f = lc(f) · ∏ sᵢ^eᵢ`.
    ///
    /// Each `sᵢ` is monic, of positive degree, squarefree, and pairwise coprime to the others; the
    /// pairs are in strictly ascending `eᵢ`, one entry per distinct multiplicity (so `sᵢ` is the
    /// product of all irreducible factors of exact multiplicity `eᵢ`). A nonzero constant gives an
    /// empty list.
    ///
    /// Yun's algorithm peels off the parts whose multiplicity is not divisible by p. What remains
    /// of `gcd(f, f')` is then a p-th power `h(x)^p` (Frobenius fixes `F_p`, so `h` takes the
    /// coefficients of `f` at indices divisible by p) and is decomposed recursively, multiplying
    /// its multiplicities by p. A polynomial with `f' = 0` is handled the same way from the start.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ZeroPolynomial`] for the zero polynomial.
    pub fn squarefree_decomposition(&self) -> Result<Vec<(Self, usize)>, FpPolyError> {
        let mut parts = BTreeMap::new();
        self.to_monic()?.collect_squarefree_parts(1, &mut parts)?;
        Ok(parts.into_iter().map(|(e, s)| (s, e)).collect())
    }

    /// Monic product of the distinct irreducible factors (`one` for a nonzero constant).
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ZeroPolynomial`] for the zero polynomial.
    pub fn radical(&self) -> Result<Self, FpPolyError> {
        self.squarefree_decomposition()?
            .iter()
            .try_fold(Self::one(self.modulus), |acc, (s, _)| acc.mul(s))
    }

    /// Accumulate the squarefree parts of the monic `self`, scaling every multiplicity by `scale`.
    fn collect_squarefree_parts(
        &self,
        scale: usize,
        parts: &mut BTreeMap<usize, Self>,
    ) -> Result<(), FpPolyError> {
        if self.degree().unwrap_or(0) == 0 {
            return Ok(());
        }
        let m = self.modulus;
        let mut c = self.gcd(&self.derivative())?;
        let mut w = self.div_rem(&c)?.0;
        let mut i = 1;
        while w.degree() != Some(0) {
            let y = w.gcd(&c)?;
            let z = w.div_rem(&y)?.0;
            if z.degree() != Some(0) {
                let slot = parts.entry(i * scale).or_insert_with(|| Self::one(m));
                *slot = slot.mul(&z)?;
            }
            c = c.div_rem(&y)?.0;
            w = y;
            i += 1;
        }
        if c.degree() != Some(0) {
            // A p-th power. A modulus beyond `usize` cannot divide any exponent present, so the
            // saturated stride keeps only the constant term, which the degree check excludes.
            let stride = usize::try_from(m.get()).unwrap_or(usize::MAX);
            let root = c.coeffs.iter().copied().step_by(stride).collect();
            Self::from_reduced(m, root).collect_squarefree_parts(scale * stride, parts)?;
        }
        Ok(())
    }

    /// Rabin's exact irreducibility test.
    ///
    /// Degree < 1 is false (zero and units are not irreducible); degree 1 is true. For monic
    /// `g = to_monic(self)` of degree `n`, `g` is irreducible iff `x^(p^n) ≡ x mod g` and
    /// `gcd(x^(p^(n/r)) − x, g) = 1` for every prime `r | n`. The first condition also rejects
    /// repeated factors, since `x^(p^n) − x` is squarefree. The powers come from `n` successive
    /// p-th powers modulo `g`. A root search never substitutes for this test.
    #[must_use]
    pub fn is_irreducible(&self) -> bool {
        // Every operation inside is on one modulus with a nonzero polynomial modulus, so it
        // cannot fail; `false` is the conservative answer if that ever stopped holding.
        self.rabin().unwrap_or(false)
    }

    fn rabin(&self) -> Result<bool, FpPolyError> {
        let n = match self.degree() {
            None | Some(0) => return Ok(false),
            Some(1) => return Ok(true),
            Some(n) => n,
        };
        let m = self.modulus;
        let g = self.to_monic()?;
        let x = Self::monomial(m, 1, 1); // degree 1 < n, already reduced mod g
        let checkpoints: Vec<usize> = prime_divisors(n).into_iter().map(|r| n / r).collect();
        let mut frobenius_at = Vec::with_capacity(checkpoints.len());
        let mut h = x.clone();
        for k in 1..=n {
            h = h.pow_mod_u64(m.get(), &g)?;
            if checkpoints.contains(&k) {
                frobenius_at.push(h.clone());
            }
        }
        if h != x {
            return Ok(false);
        }
        for hk in frobenius_at {
            if hk.sub(&x)?.gcd(&g)?.degree() != Some(0) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Number of distinct roots in `F_p`: `deg gcd(f, (x^p mod f) − x)`.
    ///
    /// Costs `O(log p)` products of polynomials of degree below `deg f` (polynomial, unbudgeted).
    /// Complete splitting into distinct linear factors holds iff [`FpPoly::is_squarefree`] and the
    /// count equals the degree.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ZeroPolynomial`] for the zero polynomial (every element is a root).
    pub fn distinct_root_count(&self) -> Result<usize, FpPolyError> {
        if self.is_zero() {
            return Err(FpPolyError::ZeroPolynomial);
        }
        let x = Self::monomial(self.modulus, 1, 1);
        let xp = x.pow_mod_u64(self.modulus.get(), self)?;
        Ok(self.gcd(&xp.sub(&x)?)?.degree().unwrap_or(0))
    }

    /// Exhaustively enumerate the roots in `F_p` within an evaluation budget.
    ///
    /// The exact distinct-root count is computed first; zero roots returns an empty `Complete`
    /// with no evaluations. Otherwise candidates `0, 1, 2, …` are evaluated in ascending order,
    /// one evaluation each, until the count-th root is found (`Complete`) or the budget is spent
    /// (`Truncated`, whose `found` lists exactly the roots below `evaluations`). The budget is
    /// never exceeded and no randomness is used; retry with a larger budget to continue.
    ///
    /// # Errors
    ///
    /// [`FpPolyError::ZeroPolynomial`] for the zero polynomial.
    pub fn roots(&self, limits: RootSearchLimits) -> Result<RootEnumeration, FpPolyError> {
        let count = self.distinct_root_count()?;
        let mut found = Vec::with_capacity(count.min(self.coeffs.len()));
        let mut evaluations = 0u64;
        if count == 0 {
            return Ok(RootEnumeration::Complete { roots: found, evaluations });
        }
        for a in 0..self.modulus.get() {
            if evaluations == limits.max_evaluations {
                return Ok(RootEnumeration::Truncated {
                    found,
                    distinct_root_count: count,
                    evaluations,
                });
            }
            evaluations += 1;
            if self.eval(a) == 0 {
                found.push(a);
                if found.len() == count {
                    break;
                }
            }
        }
        // Reaching here with `found.len() < count` would mean all of `F_p` was scanned, in which
        // case `found` is complete regardless.
        Ok(RootEnumeration::Complete { roots: found, evaluations })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduce_bigint_handles_sign_and_exact_multiples() {
        let m = PrimeModulus::new(7).unwrap();
        assert_eq!(m.reduce_bigint(&BigInt::from(-14)), 0);
        assert_eq!(m.reduce_bigint(&BigInt::from(-15)), 6);
        assert_eq!(m.reduce_bigint(&BigInt::from(15)), 1);
    }

    #[test]
    fn prime_divisors_of_composite_degrees() {
        assert_eq!(prime_divisors(12), vec![2, 3]);
        assert_eq!(prime_divisors(7), vec![7]);
        assert_eq!(prime_divisors(8), vec![2]);
    }

    #[test]
    fn zero_is_trimmed_to_empty() {
        let m = PrimeModulus::new(5).unwrap();
        assert!(FpPoly::from_coeffs(m, &[5, 10, 0]).coeffs().is_empty());
        assert!(FpPoly::monomial(m, 5, 3).is_zero());
    }

    #[test]
    fn display_covers_every_variant() {
        let errors = [
            FpPolyError::NotPrime { candidate: 9 },
            FpPolyError::ModulusMismatch { left: 3, right: 5 },
            FpPolyError::DivisionByZero,
            FpPolyError::ZeroPolynomial,
        ];
        for e in errors {
            assert!(!e.to_string().is_empty());
        }
    }
}
