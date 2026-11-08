//! Known-answer tests for checked polynomial arithmetic over prime fields.
//!
//! Every fixture is justified independently of the code under test: by hand computation, by a
//! structural fact (a factorisation written out explicitly, a Frobenius identity, the order of the
//! unit group), or by exhaustive brute force over a small field.

use num_bigint::{BigInt, BigUint};
use shared_numfield::{
    FpPoly, FpPolyError, IntPoly, PrimeModulus, RootEnumeration, RootSearchLimits, is_prime_u64,
};

/// Largest prime below 2^64; 2^64 − 59 ≡ 1 (mod 4).
const P64: u64 = u64::MAX - 58;
/// Mersenne prime 2^61 − 1; ≡ 3 (mod 4).
const P61: u64 = (1 << 61) - 1;

fn pm(p: u64) -> PrimeModulus {
    PrimeModulus::new(p).unwrap()
}

/// Polynomial over F_p from signed least-significant-first coefficients.
fn fp(p: u64, coeffs: &[i64]) -> FpPoly {
    let m = pm(p);
    let reduced: Vec<u64> = coeffs.iter().map(|&c| m.reduce_i64(c)).collect();
    FpPoly::from_coeffs(m, &reduced)
}

fn big(n: u64) -> BigUint {
    BigUint::from(n)
}

/// `base^e mod modulus` by naive repeated multiplication: independent of `pow_mod`.
fn naive_pow_mod(base: &FpPoly, e: u64, modulus: &FpPoly) -> FpPoly {
    let mut acc = FpPoly::one(base.modulus()).rem(modulus).unwrap();
    for _ in 0..e {
        acc = acc.mul_mod(base, modulus).unwrap();
    }
    acc
}

// ─── primality and modulus validation ────────────────────────────────────────

#[test]
fn kat_prime_modulus_accepts_small_and_boundary_primes() {
    // 2, 3, 5, 7: smallest primes (2 is the only even one). 2^31 − 1 and 2^61 − 1 are Mersenne
    // primes. 2^64 − 59 is the largest prime below 2^64 (standard table of primes just below
    // powers of two: 2^64 − k is prime first at k = 59).
    for p in [2, 3, 5, 7, (1u64 << 31) - 1, P61, P64] {
        let m = PrimeModulus::new(p).unwrap();
        assert_eq!(m.get(), p);
        assert_eq!(m.to_biguint(), big(p));
        assert!(is_prime_u64(p));
    }
}

#[test]
fn kat_prime_modulus_rejects_units_and_composites() {
    // 0 and 1 are not prime; 4 = 2², 9 = 3², 561 = 3·11·17 (least Carmichael number);
    // 2^64 − 1 = 3·5·17·257·641·65537·6700417 (the prime factors of F_0 … F_5 combined, since
    // 2^64 − 1 = F_0·F_1·F_2·F_3·F_4·F_5 with F_5 = 641·6700417); verified by product below.
    for n in [0, 1, 4, 9, 561, u64::MAX] {
        assert_eq!(PrimeModulus::new(n), Err(FpPolyError::NotPrime { candidate: n }));
    }
    let product: u128 = [3u128, 5, 17, 257, 641, 65537, 6_700_417].iter().product();
    assert_eq!(product, u128::from(u64::MAX));
    assert_eq!(561, 3 * 11 * 17);
}

#[test]
fn kat_prime_modulus_rejects_strong_pseudoprimes() {
    // Each n is the least strong pseudoprime to the listed prime bases (or the bound for the
    // first nine/ten/eleven bases), so a Miller–Rabin test with too few bases accepts it. The
    // factorisations are explicit and checked by multiplication, independently of any
    // Miller–Rabin computation: a nontrivial divisor proves compositeness.
    let cases: [(u64, &[u64]); 8] = [
        (2047, &[23, 89]),                          // psi_1, base 2
        (1_373_653, &[829, 1657]),                  // psi_2, bases 2, 3
        (25_326_001, &[2251, 11251]),               // psi_3, bases 2, 3, 5
        (3_215_031_751, &[151, 751, 28351]),        // psi_4, bases 2, 3, 5, 7
        (2_152_302_898_747, &[6763, 10627, 29947]), // psi_5, bases 2..11
        (3_474_749_660_383, &[1303, 16927, 157_543]), // psi_6, bases 2..13
        (341_550_071_728_321, &[10_670_053, 32_010_157]), // psi_7, bases 2..17
        (3_825_123_056_546_413_051, &[149_491, 747_451, 34_233_211]), // psi_9 = psi_10 = psi_11
    ];
    for (n, factors) in cases {
        assert!(factors.len() >= 2, "{n}: need a nontrivial factorisation");
        let divisor = factors[0];
        assert!(divisor > 1 && divisor < n && n % divisor == 0, "{n}: {divisor} is not a divisor");
        let product: u128 = factors.iter().map(|&f| u128::from(f)).product();
        assert_eq!(product, u128::from(n), "{n} is not the product of {factors:?}");
        assert!(!is_prime_u64(n), "{n} accepted as prime");
        assert_eq!(PrimeModulus::new(n), Err(FpPolyError::NotPrime { candidate: n }));
    }
}

#[test]
fn kat_is_prime_u64_matches_trial_division_below_100000() {
    // Independent oracle: sieve of Eratosthenes.
    const LIMIT: usize = 100_000;
    let mut sieve = vec![true; LIMIT];
    sieve[0] = false;
    sieve[1] = false;
    let mut i = 2;
    while i * i < LIMIT {
        if sieve[i] {
            for multiple in (i * i..LIMIT).step_by(i) {
                sieve[multiple] = false;
            }
        }
        i += 1;
    }
    for (n, &expected) in sieve.iter().enumerate() {
        assert_eq!(is_prime_u64(n as u64), expected, "disagreement at {n}");
    }
}

// ─── scalar arithmetic ───────────────────────────────────────────────────────

#[test]
fn kat_scalar_arithmetic_near_u64_max() {
    let m = pm(P64);
    let top = P64 - 1; // ≡ −1
    // (−1) + (−1) = −2; the u64 sum would overflow.
    assert_eq!(m.add(top, top), P64 - 2);
    // (−1)² = 1; the u64 product would overflow.
    assert_eq!(m.mul(top, top), 1);
    assert_eq!(m.sub(0, 1), top);
    assert_eq!(m.sub(1, 0), 1);
    assert_eq!(m.neg(0), 0);
    assert_eq!(m.neg(1), top);
    // (−1) is its own inverse, and a · a⁻¹ = 1.
    assert_eq!(m.inv(top).unwrap(), top);
    for a in [2u64, 3, 12345, P64 / 2, top - 1] {
        assert_eq!(m.mul(a, m.inv(a).unwrap()), 1);
        // Fermat's little theorem: a^(p−1) = 1 and a^p = a.
        assert_eq!(m.pow(a, P64 - 1), 1);
        assert_eq!(m.pow(a, P64), a);
    }
    assert_eq!(m.pow(7, 0), 1);
    assert_eq!(m.inv(0), Err(FpPolyError::DivisionByZero));
    assert_eq!(m.inv(P64), Err(FpPolyError::DivisionByZero));
    // Operands need not be reduced: u64::MAX = P64 + 58.
    assert_eq!(m.reduce(u64::MAX), 58);
    assert_eq!(m.add(u64::MAX, u64::MAX), 116);
    assert_eq!(m.mul(u64::MAX, u64::MAX), 58 * 58);
    // −2^63 ≡ p − 2^63 = 2^63 − 59.
    assert_eq!(m.reduce_i64(i64::MIN), (1u64 << 63) - 59);
    assert_eq!(m.reduce_i64(i64::MAX), (1u64 << 63) - 1);
    // 2^64 ≡ 59, so 2^100 = 2^36 · 2^64 ≡ 59 · 2^36 < p, and −2^100 ≡ p − 59 · 2^36.
    let two_100 = BigInt::from(1u8) << 100;
    let residue = 59 * (1u64 << 36);
    assert_eq!(m.reduce_bigint(&two_100), residue);
    assert_eq!(m.reduce_bigint(&-two_100), P64 - residue);
    assert_eq!(m.reduce_bigint(&BigInt::from(P64)), 0);
    assert_eq!(m.reduce_bigint(&-BigInt::from(P64)), 0);
    // p = 2: inverse of 1 is 1 (exponent p − 2 = 0).
    assert_eq!(pm(2).inv(1).unwrap(), 1);
    assert_eq!(pm(2).inv(3).unwrap(), 1);
}

// ─── canonical representation ────────────────────────────────────────────────

#[test]
fn kat_canonical_representation() {
    let m5 = pm(5);
    // 7 ≡ 2; 5 and 10 ≡ 0, so the trailing zeros are trimmed away.
    let f = FpPoly::from_coeffs(m5, &[7, 5, 0, 10]);
    assert_eq!(f.coeffs(), &[2]);
    assert_eq!(f.degree(), Some(0));
    // Negative coefficients take their least nonnegative residue: −1 ≡ 4, −7 ≡ 3 (mod 5).
    let int = IntPoly::from_coeffs(vec![BigInt::from(-1), BigInt::from(0), BigInt::from(1)]);
    let g = FpPoly::from_int_poly(&int, m5);
    assert_eq!(g.coeffs(), &[4, 0, 1]);
    assert_eq!(fp(5, &[-7, 1]).coeffs(), &[3, 1]);
    // Huge integer coefficients: 10^30 ≡ 0 and −10^30 − 1 ≡ 4 (mod 5).
    let ten30: BigInt = BigInt::from(10).pow(30);
    let h = IntPoly::from_coeffs(vec![ten30.clone(), -ten30 - 1]);
    assert_eq!(FpPoly::from_int_poly(&h, m5).coeffs(), &[0, 4]);
    // Degree drop: 2x² + x + 1 mod 2 = x + 1.
    let two = IntPoly::from_coeffs(vec![BigInt::from(1), BigInt::from(1), BigInt::from(2)]);
    let dropped = FpPoly::from_int_poly(&two, pm(2));
    assert_eq!(dropped.coeffs(), &[1, 1]);
    assert_eq!(dropped.degree(), Some(1));
    assert_eq!(FpPoly::from_int_poly(&two, pm(3)).degree(), Some(2));
    // Round trip through the integer representation uses representatives in [0, p).
    let back = g.to_int_poly();
    assert_eq!(back.coeffs, vec![BigInt::from(4), BigInt::from(0), BigInt::from(1)]);
    assert_eq!(FpPoly::from_int_poly(&back, m5), g);
    // Zero is the empty slice; monic detection follows the canonical leading coefficient.
    assert!(FpPoly::zero(m5).coeffs().is_empty());
    assert_eq!(FpPoly::zero(m5).degree(), None);
    assert_eq!(FpPoly::zero(m5).leading_coeff(), None);
    assert!(!FpPoly::zero(m5).is_monic());
    assert!(FpPoly::monomial(m5, 1, 3).is_monic());
    assert!(!FpPoly::monomial(m5, 2, 3).is_monic());
    assert_eq!(FpPoly::monomial(m5, 7, 2).coeffs(), &[0, 0, 2]);
    assert_eq!(FpPoly::constant(m5, 12).coeffs(), &[2]);
    assert_eq!(FpPoly::constant(m5, 10), FpPoly::zero(m5));
    // Equality includes the modulus: the same coefficients over F_3 and F_5 differ.
    assert_ne!(fp(3, &[1, 1]), fp(5, &[1, 1]));
    assert_eq!(fp(3, &[1, 1]).modulus().get(), 3);
    assert_ne!(FpPoly::zero(pm(3)), FpPoly::zero(pm(5)));
}

#[test]
fn kat_modulus_mismatch_is_error() {
    fn mismatch<T>() -> Result<T, FpPolyError> {
        Err(FpPolyError::ModulusMismatch { left: 3, right: 5 })
    }
    let a = fp(3, &[1, 1]);
    let b = fp(5, &[1, 1]);
    assert_eq!(a.add(&b), mismatch());
    assert_eq!(a.sub(&b), mismatch());
    assert_eq!(a.mul(&b), mismatch());
    assert_eq!(a.div_rem(&b), mismatch());
    assert_eq!(a.rem(&b), mismatch());
    assert_eq!(a.gcd(&b), mismatch());
    assert_eq!(a.extended_gcd(&b), mismatch());
    assert_eq!(a.mul_mod(&a, &b), mismatch());
    assert_eq!(a.mul_mod(&b, &a), mismatch());
    assert_eq!(a.pow_mod(&big(2), &b), mismatch());
    // The mismatch is reported even when the right operand is zero.
    assert_eq!(a.div_rem(&FpPoly::zero(pm(5))), mismatch());
    assert_eq!(
        b.add(&a),
        Err(FpPolyError::ModulusMismatch { left: 5, right: 3 }),
        "left and right follow operand order"
    );
}

// ─── division, gcd, powers ───────────────────────────────────────────────────

#[test]
fn kat_div_rem_identity_and_degree() {
    // Synthetic division of x³ + 2x + 3 by x − 1 over F_7: quotient x² + x + 3 and remainder
    // f(1) = 6.
    let f = fp(7, &[3, 2, 0, 1]);
    let (q, r) = f.div_rem(&fp(7, &[-1, 1])).unwrap();
    assert_eq!(q, fp(7, &[3, 1, 1]));
    assert_eq!(r, fp(7, &[6]));
    assert_eq!(f.rem(&fp(7, &[-1, 1])).unwrap(), r);
    // Including non-monic and constant divisors, every quotient/remainder pair must reassemble
    // and have a smaller remainder degree.
    let dividends = [fp(7, &[1, 2, 3, 4, 5, 6, 1]), fp(7, &[0, 0, 0, 3]), fp(7, &[5]), fp(7, &[])];
    let divisors = [fp(7, &[3, 5]), fp(7, &[1, 0, 2]), fp(7, &[4]), fp(7, &[2, 1, 1, 3])];
    for a in &dividends {
        for b in &divisors {
            let (q, r) = a.div_rem(b).unwrap();
            assert_eq!(&q.mul(b).unwrap().add(&r).unwrap(), a);
            assert!(r.is_zero() || r.degree() < b.degree(), "remainder degree for {a:?} / {b:?}");
        }
    }
    // Lower degree divided by higher: quotient 0, remainder is the dividend.
    let small = fp(7, &[1, 1]);
    let (q, r) = small.div_rem(&fp(7, &[2, 0, 1])).unwrap();
    assert!(q.is_zero());
    assert_eq!(r, small);
    // Division by a nonzero constant is exact and scales.
    let (q, r) = fp(7, &[2, 4]).div_rem(&fp(7, &[2])).unwrap();
    assert_eq!(q, fp(7, &[1, 2]));
    assert!(r.is_zero());
    // Zero divisor.
    assert_eq!(f.div_rem(&FpPoly::zero(pm(7))), Err(FpPolyError::DivisionByZero));
    assert_eq!(f.rem(&FpPoly::zero(pm(7))), Err(FpPolyError::DivisionByZero));
}

#[test]
fn kat_arithmetic_basics() {
    // Over F_5: (x + 1)(x + 4) = x² + 5x + 4 = x² + 4. Hand expansion.
    let a = fp(5, &[1, 1]);
    let b = fp(5, &[4, 1]);
    assert_eq!(a.mul(&b).unwrap(), fp(5, &[4, 0, 1]));
    assert_eq!(a.add(&b).unwrap(), fp(5, &[0, 2]));
    assert_eq!(a.sub(&a).unwrap(), FpPoly::zero(pm(5)));
    // Leading coefficients cancel and the sum is trimmed: (x² + 1) + (4x² + 2) = 3.
    assert_eq!(fp(5, &[1, 0, 1]).add(&fp(5, &[2, 0, 4])).unwrap(), fp(5, &[3]));
    assert_eq!(a.neg(), fp(5, &[4, 4]));
    assert_eq!(a.scale(3), fp(5, &[3, 3]));
    assert!(a.scale(5).is_zero());
    assert!(a.mul(&FpPoly::zero(pm(5))).unwrap().is_zero());
    // Products near u64::MAX do not overflow: ((p − 1)x + (p − 1))² = x² + 2x + 1 = (x + 1)².
    let neg_one_x_plus_one = FpPoly::from_coeffs(pm(P64), &[P64 - 1, P64 - 1]);
    assert_eq!(
        neg_one_x_plus_one.mul(&neg_one_x_plus_one).unwrap(),
        FpPoly::from_coeffs(pm(P64), &[1, 2, 1])
    );
    // Horner evaluation: x² + 4 at 1 over F_5 is 0; at 2 is 8 ≡ 3; the argument is reduced.
    let f = fp(5, &[4, 0, 1]);
    assert_eq!(f.eval(1), 0);
    assert_eq!(f.eval(2), 3);
    assert_eq!(f.eval(7), 3);
    assert_eq!(FpPoly::zero(pm(5)).eval(3), 0);
    // (x + 1)^2 at P64 − 1 over F_P64 is 0.
    assert_eq!(FpPoly::from_coeffs(pm(P64), &[1, 1]).eval(P64 - 1), 0);
}

#[test]
fn kat_gcd_and_extended_gcd() {
    // (x − 1)(x − 2) and (x − 1)(x − 3) over F_7 share exactly x − 1 (2 ≠ 3 mod 7).
    let m = pm(7);
    let x1 = fp(7, &[-1, 1]);
    let a = x1.mul(&fp(7, &[-2, 1])).unwrap();
    let b = x1.mul(&fp(7, &[-3, 1])).unwrap();
    assert_eq!(a.gcd(&b).unwrap(), x1);
    // Non-monic multiples still give the monic gcd.
    assert_eq!(a.scale(3).gcd(&b.scale(5)).unwrap(), x1);
    let eg = a.extended_gcd(&b).unwrap();
    assert_eq!(eg.gcd, x1);
    let bezout = eg.s.mul(&a).unwrap().add(&eg.t.mul(&b).unwrap()).unwrap();
    assert_eq!(bezout, eg.gcd);
    // Coprime: x and x + 1 have gcd 1, and Bézout holds exactly.
    let p = fp(7, &[0, 1]);
    let q = fp(7, &[1, 1]);
    assert_eq!(p.gcd(&q).unwrap(), FpPoly::one(m));
    let eg = p.extended_gcd(&q).unwrap();
    assert_eq!(eg.gcd, FpPoly::one(m));
    assert_eq!(eg.s.mul(&p).unwrap().add(&eg.t.mul(&q).unwrap()).unwrap(), FpPoly::one(m));
    // Degenerate operands: gcd(0, 0) = 0; gcd(0, f) and gcd(f, 0) are the monic form of f.
    let zero = FpPoly::zero(m);
    assert!(zero.gcd(&zero).unwrap().is_zero());
    let eg = zero.extended_gcd(&zero).unwrap();
    assert!(eg.gcd.is_zero());
    let f = fp(7, &[3, 3]); // 3x + 3 = 3(x + 1)
    assert_eq!(zero.gcd(&f).unwrap(), fp(7, &[1, 1]));
    assert_eq!(f.gcd(&zero).unwrap(), fp(7, &[1, 1]));
    for (u, v) in [(&zero, &f), (&f, &zero)] {
        let eg = u.extended_gcd(v).unwrap();
        assert_eq!(eg.gcd, fp(7, &[1, 1]));
        assert_eq!(eg.s.mul(u).unwrap().add(&eg.t.mul(v).unwrap()).unwrap(), eg.gcd);
    }
    // Bézout over a large prime: gcd of (x + 1)(x + 2) and (x + 2)(x + 3).
    let big_m = pm(P64);
    let lin = |c: u64| FpPoly::from_coeffs(big_m, &[c, 1]);
    let a = lin(1).mul(&lin(2)).unwrap();
    let b = lin(2).mul(&lin(3)).unwrap();
    let eg = a.extended_gcd(&b).unwrap();
    assert_eq!(eg.gcd, lin(2));
    assert_eq!(eg.s.mul(&a).unwrap().add(&eg.t.mul(&b).unwrap()).unwrap(), eg.gcd);
}

#[test]
fn kat_to_monic() {
    assert_eq!(fp(7, &[3, 6]).to_monic().unwrap(), fp(7, &[4, 1])); // 6⁻¹ = 6; 3·6 = 18 ≡ 4
    assert_eq!(fp(7, &[3, 1]).to_monic().unwrap(), fp(7, &[3, 1]));
    assert_eq!(FpPoly::zero(pm(7)).to_monic(), Err(FpPolyError::ZeroPolynomial));
}

#[test]
fn kat_pow_mod_matches_repeated_multiplication_and_frobenius() {
    // f = x² + 1 is irreducible over F_3 (−1 is a nonresidue mod 3), so F_3[x]/f is the field
    // F_9, whose unit group has order 8, and Frobenius x ↦ x³ has order 2: x^9 = x.
    let f = fp(3, &[1, 0, 1]);
    let x = fp(3, &[0, 1]);
    assert_eq!(x.pow_mod(&big(9), &f).unwrap(), x, "x^(p^d) = x");
    assert_eq!(x.pow_mod(&big(3), &f).unwrap(), fp(3, &[0, 2]), "x^3 = −x since x² = −1");
    // pow_mod agrees with naive repeated multiplication.
    let base = fp(3, &[1, 1]);
    for e in 0..=12 {
        assert_eq!(base.pow_mod(&big(e), &f).unwrap(), naive_pow_mod(&base, e, &f), "e = {e}");
    }
    // An exponent above 2^64 reduces mod the unit-group order p^d − 1 = 8: 2^70 + 5 ≡ 5.
    let huge = (BigUint::from(1u8) << 70) + 5u8;
    assert_eq!(base.pow_mod(&huge, &f).unwrap(), naive_pow_mod(&base, 5, &f));
    // Frobenius in a larger field: x⁴ + x + 1 is irreducible over F_2, so x^(2^4) = x; the same
    // exponent as p^d with d = 4 and with a BigUint power.
    let g = fp(2, &[1, 1, 0, 0, 1]);
    let x2 = fp(2, &[0, 1]);
    assert_eq!(x2.pow_mod(&big(16), &g).unwrap(), x2);
    // x² + 1 is irreducible over F_(2^61−1) (≡ 3 mod 4), so x^(p²) = x with a 122-bit exponent.
    let p_to_d = BigUint::from(P61).pow(2);
    let h = FpPoly::from_coeffs(pm(P61), &[1, 0, 1]);
    let xh = FpPoly::monomial(pm(P61), 1, 1);
    assert_eq!(xh.pow_mod(&p_to_d, &h).unwrap(), xh, "Frobenius over a 61-bit prime");
    // Exponent 0 is 1 rem modulus; a nonzero constant modulus collapses everything to zero.
    assert_eq!(x.pow_mod(&BigUint::from(0u8), &f).unwrap(), FpPoly::one(pm(3)));
    assert!(x.pow_mod(&BigUint::from(0u8), &fp(3, &[2])).unwrap().is_zero());
    assert!(x.pow_mod(&big(5), &fp(3, &[2])).unwrap().is_zero());
    // Unreduced base: x^3 + 2 mod (x² + 1) over F_3 reduces as (x·x² + 2) = −x + 2.
    let unreduced = fp(3, &[2, 0, 0, 1]);
    assert_eq!(unreduced.pow_mod(&big(1), &f).unwrap(), fp(3, &[2, 2]));
    // Zero modulus is an error even for exponent 0.
    assert_eq!(
        x.pow_mod(&BigUint::from(0u8), &FpPoly::zero(pm(3))),
        Err(FpPolyError::DivisionByZero)
    );
    assert_eq!(x.mul_mod(&x, &FpPoly::zero(pm(3))), Err(FpPolyError::DivisionByZero));
}

#[test]
fn kat_derivative_characteristic_p() {
    // d/dx x^p = p·x^(p−1) = 0 in characteristic p.
    assert!(fp(5, &[0, 0, 0, 0, 0, 1]).derivative().is_zero());
    // x^5 + x over F_5: derivative 5x⁴ + 1 = 1.
    assert_eq!(fp(5, &[0, 1, 0, 0, 0, 1]).derivative(), fp(5, &[1]));
    // x^7 + 2x³ over F_7: 7x⁶ + 6x² = 6x².
    assert_eq!(fp(7, &[0, 0, 0, 2, 0, 0, 0, 1]).derivative(), fp(7, &[0, 0, 6]));
    // 3x² + x + 7 over F_5: 6x + 1 = x + 1.
    assert_eq!(fp(5, &[7, 1, 3]).derivative(), fp(5, &[1, 1]));
    // Constants and zero differentiate to zero.
    assert!(fp(5, &[4]).derivative().is_zero());
    assert!(FpPoly::zero(pm(5)).derivative().is_zero());
    // Over F_2 the derivative of x² + x + 1 is 1 (2x vanishes), so it is squarefree.
    assert_eq!(fp(2, &[1, 1, 1]).derivative(), fp(2, &[1]));
}

// ─── squarefreeness and irreducibility ───────────────────────────────────────

#[test]
fn kat_irreducible_fixtures() {
    // Each fixture is irreducible by an independent argument.
    let fixtures: [(u64, Vec<i64>, &str); 9] = [
        // −1 is a nonresidue mod 3 (squares are 0, 1).
        (3, vec![1, 0, 1], "x²+1 / F_3"),
        // 2 is a nonresidue mod 3 and mod 5 (squares mod 5 are 0, 1, 4).
        (3, vec![-2, 0, 1], "x²−2 / F_3"),
        (5, vec![-2, 0, 1], "x²−2 / F_5"),
        // A cubic without a root is irreducible: x³−x−1 takes the values −1, −1 at 0, 1 over F_2
        // and F_3 and 5 ≡ 2 at 2 over F_3, never 0.
        (2, vec![-1, -1, 0, 1], "x³−x−1 / F_2"),
        (3, vec![-1, -1, 0, 1], "x³−x−1 / F_3"),
        // The standard primitive polynomial for GF(16).
        (2, vec![1, 1, 0, 0, 1], "x⁴+x+1 / F_2"),
        // The standard primitive polynomial for GF(64).
        (2, vec![1, 1, 0, 0, 0, 0, 1], "x⁶+x+1 / F_2"),
        // 2^61 − 1 ≡ 3 (mod 4), so −1 is a nonresidue.
        (P61, vec![1, 0, 1], "x²+1 / F_(2^61−1)"),
        // Linear polynomials: degree one is irreducible.
        (7, vec![3, 1], "x+3 / F_7"),
    ];
    for (p, coeffs, label) in fixtures {
        let f = fp(p, &coeffs);
        assert!(f.is_irreducible(), "{label} must be irreducible");
        assert!(f.is_squarefree(), "{label} must be squarefree");
        // A nonzero scalar multiple is irreducible too (units do not change irreducibility);
        // over F_2 the only unit is 1, so the check is vacuous there.
        let unit = if p == 2 { 1 } else { 2 };
        assert!(f.scale(unit).is_irreducible(), "{label} scaled");
        // An irreducible polynomial of degree ≥ 2 has no roots.
        if f.degree() >= Some(2) {
            assert_eq!(f.distinct_root_count().unwrap(), 0, "{label} has no roots");
        }
    }
    // Non-monic linear polynomial.
    assert!(fp(5, &[1, 2]).is_irreducible());
    // Cross-check irreducibility with the Frobenius characterisation: an irreducible f of degree
    // d satisfies x^(p^d) = x mod f; (x²+1)(x²+x+2)'s Frobenius fails.
    let f = fp(3, &[1, 0, 1]);
    let x = fp(3, &[0, 1]);
    assert_eq!(x.pow_mod(&big(9), &f).unwrap(), x);
}

#[test]
fn kat_reducible_fixtures_without_roots() {
    // (x²+1)(x²+x+2) over F_3: both factors are irreducible (discriminants −4 ≡ 2 and
    // 1 − 8 ≡ 2, nonresidues mod 3), and the product expands to x⁴ + x³ + 3x² + x + 2 =
    // x⁴ + x³ + x + 2.
    let a = fp(3, &[1, 0, 1]);
    let b = fp(3, &[2, 1, 1]);
    assert!(a.is_irreducible() && b.is_irreducible());
    let product = a.mul(&b).unwrap();
    assert_eq!(product, fp(3, &[2, 1, 0, 1, 1]));
    // t⁴ − 10t² + 1 ≡ t⁴ + 1 = (t²+2)(t²+3) over F_5 (t⁴ + 5t² + 6); −2 ≡ 3 and −3 ≡ 2 are
    // nonresidues mod 5, so both quadratics are irreducible.
    let t4 = fp(5, &[1, 0, -10, 0, 1]);
    assert_eq!(t4, fp(5, &[2, 0, 1]).mul(&fp(5, &[3, 0, 1])).unwrap());
    assert!(fp(5, &[2, 0, 1]).is_irreducible() && fp(5, &[3, 0, 1]).is_irreducible());
    // (x³+x+1)(x³+x²+1) = x⁶ + x⁵ + x⁴ + x³ + x² + x + 1 = (x⁷ − 1)/(x − 1) over F_2.
    let c = fp(2, &[1, 1, 0, 1]);
    let d = fp(2, &[1, 0, 1, 1]);
    assert!(c.is_irreducible() && d.is_irreducible());
    let seventh = fp(2, &[1, 1, 1, 1, 1, 1, 1]);
    assert_eq!(c.mul(&d).unwrap(), seventh);
    // Squarefree and reducible, and with no roots (all factors have degree ≥ 2 or f(a) ≠ 0).
    for f in [product, t4, seventh] {
        assert!(f.is_squarefree());
        assert!(!f.is_irreducible(), "{f:?} is a product of two irreducibles");
        assert_eq!(f.distinct_root_count().unwrap(), 0);
    }
}

#[test]
fn kat_repeated_factor_fixtures() {
    // (x²+x+1)² = x⁴ + x² + 1 over F_2 (cross terms 2x³ and 2x vanish).
    let square = fp(2, &[1, 1, 1]).mul(&fp(2, &[1, 1, 1])).unwrap();
    assert_eq!(square, fp(2, &[1, 0, 1, 0, 1]));
    // t⁴ − 10t² + 1 ≡ t⁴ − t² + 1 = t⁴ + 2t² + 1 = (t² + 1)² over F_3.
    let t4 = fp(3, &[1, 0, -10, 0, 1]);
    assert_eq!(t4, fp(3, &[1, 0, 1]).mul(&fp(3, &[1, 0, 1])).unwrap());
    // x³ over F_3.
    let cube = fp(3, &[0, 0, 0, 1]);
    for f in [square, t4, cube] {
        assert!(!f.is_squarefree(), "{f:?} has a repeated factor");
        assert!(!f.is_irreducible(), "{f:?} is not irreducible");
    }
    // A p-th power has a zero derivative: (x + 1)^5 = x^5 + 1 over F_5.
    let fifth = fp(5, &[1, 0, 0, 0, 0, 1]);
    assert!(fifth.derivative().is_zero());
    assert!(!fifth.is_squarefree());
    assert!(!fifth.is_irreducible());
    assert_eq!(fifth.distinct_root_count().unwrap(), 1);
}

#[test]
fn kat_reducible_with_roots() {
    // x² + 1 = (x − 2)(x + 2) over F_5 (4 + 1 = 5).
    let f = fp(5, &[1, 0, 1]);
    assert!(!f.is_irreducible());
    assert_eq!(f.distinct_root_count().unwrap(), 2);
    // x³ − x − 1 = (x − 2)(x² + 2x + 3) over F_5: expanding gives x³ − x − 6 ≡ x³ − x − 1, and
    // the quadratic has discriminant 4 − 12 = −8 ≡ 2, a nonresidue, so the only root is 2.
    let g = fp(5, &[-1, -1, 0, 1]);
    assert_eq!(g, fp(5, &[-2, 1]).mul(&fp(5, &[3, 2, 1])).unwrap());
    assert!(!g.is_irreducible());
    assert!(g.is_squarefree());
    assert_eq!(g.distinct_root_count().unwrap(), 1);
    // 2^64 − 59 ≡ 1 (mod 4) makes −1 a quadratic residue, so x² + 1 has two roots.
    let h = fp(P64, &[1, 0, 1]);
    assert!(!h.is_irreducible());
    assert!(h.is_squarefree());
    assert_eq!(h.distinct_root_count().unwrap(), 2);
    // A product of distinct linear factors is squarefree, reducible, and fully split.
    let split = fp(7, &[-1, 1]).mul(&fp(7, &[-2, 1])).unwrap().mul(&fp(7, &[-3, 1])).unwrap();
    assert!(split.is_squarefree() && !split.is_irreducible());
    assert_eq!(split.distinct_root_count(), Ok(3));
}

#[test]
fn kat_degree_and_zero_handling() {
    let m = pm(5);
    let zero = FpPoly::zero(m);
    // Zero: every element is a root, so root-based queries are rejected.
    assert!(!zero.is_irreducible());
    assert!(!zero.is_squarefree());
    assert_eq!(zero.distinct_root_count(), Err(FpPolyError::ZeroPolynomial));
    let limits = RootSearchLimits { max_evaluations: 10 };
    assert_eq!(zero.roots(limits), Err(FpPolyError::ZeroPolynomial));
    assert_eq!(zero.to_monic(), Err(FpPolyError::ZeroPolynomial));
    // Nonzero constants are units: squarefree, not irreducible, no roots.
    let unit = fp(5, &[3]);
    assert!(!unit.is_irreducible());
    assert!(unit.is_squarefree());
    assert_eq!(unit.distinct_root_count(), Ok(0));
    assert_eq!(unit.degree(), Some(0));
    // A degree-2 integer polynomial can drop to a unit: 5x² + 3 mod 5.
    let int = IntPoly::from_coeffs(vec![BigInt::from(3), BigInt::from(0), BigInt::from(5)]);
    let dropped = FpPoly::from_int_poly(&int, m);
    assert_eq!(dropped.degree(), Some(0));
    assert!(!dropped.is_irreducible());
    // Linear polynomials are irreducible and have exactly one root.
    let lin = fp(5, &[2, 3]);
    assert!(lin.is_irreducible() && lin.is_squarefree());
    assert_eq!(lin.distinct_root_count(), Ok(1));
    // Exact degree: a product has the sum of the degrees (no zero divisors in F_p[x]).
    let a = fp(P64, &[1, 2, 3]);
    let b = fp(P64, &[5, 0, 0, 7]);
    assert_eq!(a.mul(&b).unwrap().degree(), Some(5));
}

// ─── root counting and enumeration ───────────────────────────────────────────

#[test]
fn kat_distinct_root_count() {
    // x⁵ − x vanishes on all of F_5 (Fermat).
    assert_eq!(fp(5, &[0, -1, 0, 0, 0, 1]).distinct_root_count(), Ok(5));
    // (x − 1)²(x − 2) has distinct roots {1, 2}.
    let rep = fp(5, &[-1, 1]).mul(&fp(5, &[-1, 1])).unwrap().mul(&fp(5, &[-2, 1])).unwrap();
    assert_eq!(rep.distinct_root_count(), Ok(2));
    // t⁴ + 1 ≡ (t²+2)(t²+3) over F_5, both irreducible: no roots.
    assert_eq!(fp(5, &[1, 0, 0, 0, 1]).distinct_root_count(), Ok(0));
    // x² + 1 over 2^64 − 59 ≡ 1 (mod 4): two roots, ± a square root of −1.
    assert_eq!(fp(P64, &[1, 0, 1]).distinct_root_count(), Ok(2));
    // Splitting into distinct linear factors is decided by squarefree + count == degree.
    let split = fp(5, &[0, -1, 0, 0, 0, 1]);
    assert!(split.is_squarefree() && split.distinct_root_count() == Ok(split.degree().unwrap()));
    assert!(!rep.is_squarefree());
}

#[test]
fn kat_root_enumeration_complete() {
    let ample = RootSearchLimits { max_evaluations: u64::MAX };
    // x² − 2 over F_7: 3² = 9 ≡ 2 and 4² = 16 ≡ 2; candidates 0..=4 are evaluated.
    assert_eq!(
        fp(7, &[-2, 0, 1]).roots(ample),
        Ok(RootEnumeration::Complete { roots: vec![3, 4], evaluations: 5 })
    );
    // Exactly the needed budget still completes.
    assert_eq!(
        fp(7, &[-2, 0, 1]).roots(RootSearchLimits { max_evaluations: 5 }),
        Ok(RootEnumeration::Complete { roots: vec![3, 4], evaluations: 5 })
    );
    // x⁵ − x over F_5 vanishes everywhere.
    assert_eq!(
        fp(5, &[0, -1, 0, 0, 0, 1]).roots(ample),
        Ok(RootEnumeration::Complete { roots: vec![0, 1, 2, 3, 4], evaluations: 5 })
    );
    // (x − 1)²(x − 2) over F_5: distinct roots 1, 2; found after evaluating 0, 1, 2.
    let rep = fp(5, &[-1, 1]).mul(&fp(5, &[-1, 1])).unwrap().mul(&fp(5, &[-2, 1])).unwrap();
    assert_eq!(
        rep.roots(ample),
        Ok(RootEnumeration::Complete { roots: vec![1, 2], evaluations: 3 })
    );
    // No roots: nothing is evaluated, even with a zero budget.
    for f in [fp(5, &[1, 0, 0, 0, 1]), fp(5, &[3]), fp(7, &[1, 0, 1])] {
        let expected = RootEnumeration::Complete { roots: vec![], evaluations: 0 };
        assert_eq!(f.roots(ample), Ok(expected.clone()));
        assert_eq!(f.roots(RootSearchLimits { max_evaluations: 0 }), Ok(expected));
    }
    // Exhaustive cross-check: against brute-force evaluation for every nonzero polynomial of
    // degree ≤ 3 over F_2, F_3, F_5 and F_7.
    for p in [2u64, 3, 5, 7] {
        for code in 1..p.pow(4) {
            let coeffs: Vec<u64> = (0..4).map(|i| code / p.pow(i) % p).collect();
            let f = FpPoly::from_coeffs(pm(p), &coeffs);
            if f.is_zero() {
                continue;
            }
            let brute: Vec<u64> = (0..p).filter(|&a| f.eval(a) == 0).collect();
            match f.roots(ample).unwrap() {
                RootEnumeration::Complete { roots, evaluations } => {
                    assert_eq!(roots, brute, "p = {p}, f = {coeffs:?}");
                    assert_eq!(evaluations, brute.last().map_or(0, |last| last + 1));
                }
                other @ RootEnumeration::Truncated { .. } => panic!("unexpected {other:?}"),
            }
            assert_eq!(f.distinct_root_count().unwrap(), brute.len());
        }
    }
}

#[test]
fn kat_root_enumeration_truncated() {
    // x² − 2 over F_7 has roots {3, 4}. Budget 4 evaluates 0..=3 and finds only 3.
    let f = fp(7, &[-2, 0, 1]);
    assert_eq!(
        f.roots(RootSearchLimits { max_evaluations: 4 }),
        Ok(RootEnumeration::Truncated { found: vec![3], distinct_root_count: 2, evaluations: 4 })
    );
    // Budget 0 evaluates nothing but still reports the exact count.
    assert_eq!(
        f.roots(RootSearchLimits { max_evaluations: 0 }),
        Ok(RootEnumeration::Truncated { found: vec![], distinct_root_count: 2, evaluations: 0 })
    );
    // Budget 3 stops before reaching 3.
    assert_eq!(
        f.roots(RootSearchLimits { max_evaluations: 3 }),
        Ok(RootEnumeration::Truncated { found: vec![], distinct_root_count: 2, evaluations: 3 })
    );
    // x² + 1 over 2^64 − 59 has two roots, but a budget of 1000 cannot reach them: the result is
    // resource exhaustion with the exact count, not a claim of absence. Every found value
    // (almost certainly none) is a true root.
    let h = fp(P64, &[1, 0, 1]);
    match h.roots(RootSearchLimits { max_evaluations: 1000 }).unwrap() {
        RootEnumeration::Truncated { found, distinct_root_count, evaluations } => {
            assert_eq!(distinct_root_count, 2);
            assert_eq!(evaluations, 1000);
            assert!(found.len() < distinct_root_count);
            assert!(found.iter().all(|&a| h.eval(a) == 0));
            let brute: Vec<u64> = (0..evaluations).filter(|&a| h.eval(a) == 0).collect();
            assert_eq!(found, brute, "found lists exactly the roots in [0, evaluations)");
        }
        other @ RootEnumeration::Complete { .. } => panic!("unexpected {other:?}"),
    }
}

// ─── square-free decomposition ───────────────────────────────────────────────

/// `lc · ∏ sᵢ^eᵢ` by repeated multiplication: independent of `pow_mod`.
fn recompose(lead: u64, parts: &[(FpPoly, usize)], m: PrimeModulus) -> FpPoly {
    let mut acc = FpPoly::constant(m, lead);
    for (s, e) in parts {
        for _ in 0..*e {
            acc = acc.mul(s).unwrap();
        }
    }
    acc
}

/// Assert the structural contract of a decomposition of `f` (nonzero, positive degree or not).
fn assert_decomposition_contract(f: &FpPoly, parts: &[(FpPoly, usize)]) {
    let lead = f.leading_coeff().unwrap();
    assert_eq!(&recompose(lead, parts, f.modulus()), f, "f = lc · ∏ sᵢ^eᵢ");
    for (s, _) in parts {
        assert!(s.is_monic() && s.degree().unwrap() > 0 && s.is_squarefree(), "{s:?}");
    }
    assert!(parts.windows(2).all(|w| w[0].1 < w[1].1), "multiplicities strictly ascending");
    for (i, (a, _)) in parts.iter().enumerate() {
        for (b, _) in &parts[i + 1..] {
            assert_eq!(a.gcd(b).unwrap().degree(), Some(0), "pairwise coprime");
        }
    }
}

#[test]
fn kat_squarefree_decomposition_small_characteristics() {
    // x² over F_2: f' = 0, so f = (x)² via the p-th root.
    assert_eq!(fp(2, &[0, 0, 1]).squarefree_decomposition(), Ok(vec![(fp(2, &[0, 1]), 2)]));
    // x³ + 1 over F_3 is (x + 1)³ (Frobenius: the cross terms 3x² and 3x vanish).
    assert_eq!(fp(3, &[1, 0, 0, 1]).squarefree_decomposition(), Ok(vec![(fp(3, &[1, 1]), 3)]));
    // x²(x+1)³(x+2)⁶ over F_3: multiplicities 2 (prime to 3), 3 and 6 (multiples of 3).
    let f = fp(3, &[0, 1]).mul(&fp(3, &[0, 1])).unwrap();
    let f = f.mul(&recompose(1, &[(fp(3, &[1, 1]), 3), (fp(3, &[2, 1]), 6)], pm(3))).unwrap();
    assert_eq!(
        f.squarefree_decomposition(),
        Ok(vec![(fp(3, &[0, 1]), 2), (fp(3, &[1, 1]), 3), (fp(3, &[2, 1]), 6)])
    );
    // x(x+1)⁵ over F_5: (x+1)⁵ = x⁵ + 1, so f = x⁶ + x with f' = x⁵ + 1 ≠ 0.
    let f = fp(5, &[0, 1, 0, 0, 0, 0, 1]);
    assert_eq!(
        f.squarefree_decomposition(),
        Ok(vec![(fp(5, &[0, 1]), 1), (fp(5, &[1, 1]), 5)])
    );
    // Same multiplicity shared by two different irreducibles merges into one entry:
    // (x+1)²(x+2)² over F_5 is ((x+1)(x+2))².
    let both = fp(5, &[1, 1]).mul(&fp(5, &[2, 1])).unwrap();
    let f = both.mul(&both).unwrap();
    assert_eq!(f.squarefree_decomposition(), Ok(vec![(both, 2)]));
}

#[test]
fn kat_squarefree_decomposition_edge_inputs() {
    // A squarefree polynomial is its own monic form with multiplicity 1; the leading
    // coefficient is discarded: 3(x² + 1) over F_7 (x² + 1 has no roots mod 7).
    let g = fp(7, &[1, 0, 1]);
    assert_eq!(g.scale(3).squarefree_decomposition(), Ok(vec![(g.clone(), 1)]));
    assert_eq!(g.radical(), Ok(g));
    // Nonzero constants have no factors; the radical is 1.
    let c = fp(7, &[4]);
    assert_eq!(c.squarefree_decomposition(), Ok(vec![]));
    assert_eq!(c.radical(), Ok(FpPoly::one(pm(7))));
    // The zero polynomial is rejected.
    let z = FpPoly::zero(pm(7));
    assert_eq!(z.squarefree_decomposition(), Err(FpPolyError::ZeroPolynomial));
    assert_eq!(z.radical(), Err(FpPolyError::ZeroPolynomial));
    // (x+1)² over p = 2^61 − 1: f' = 2x + 2 ≠ 0, so the p-th-root path is not taken.
    assert_eq!(
        fp(P61, &[1, 2, 1]).squarefree_decomposition(),
        Ok(vec![(fp(P61, &[1, 1]), 2)])
    );
    // Radical of x²(x+1)³(x+2)⁶ over F_3 is x(x+1)(x+2) = x³ + 2x (the product of all of F_3's
    // linear factors is x³ − x).
    let f = recompose(
        2,
        &[(fp(3, &[0, 1]), 2), (fp(3, &[1, 1]), 3), (fp(3, &[2, 1]), 6)],
        pm(3),
    );
    assert_eq!(f.radical(), Ok(fp(3, &[0, -1, 0, 1])));
}

#[test]
fn kat_squarefree_decomposition_exhaustive_small_fields() {
    // Every nonzero polynomial up to the stated degree: structural contract, plus the radical
    // being the squarefree monic product of the parts.
    for (p, max_len) in [(2u64, 9u32), (3, 7), (5, 5), (7, 4)] {
        for code in 1..p.pow(max_len) {
            let coeffs: Vec<u64> = (0..max_len).map(|i| code / p.pow(i) % p).collect();
            let f = FpPoly::from_coeffs(pm(p), &coeffs);
            let parts = f.squarefree_decomposition().unwrap();
            assert_decomposition_contract(&f, &parts);
            let rad = f.radical().unwrap();
            assert!(rad.is_monic() && rad.is_squarefree(), "p = {p}, f = {coeffs:?}");
            assert_eq!(f.rem(&rad).unwrap(), FpPoly::zero(pm(p)), "radical divides f");
        }
    }
}
