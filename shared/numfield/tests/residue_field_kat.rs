//! Known-answer tests for residue fields `F_p[x]/(f)`: inverses, squareness and square roots.
//!
//! Every fixture is justified independently of the code under test: by exhaustive enumeration of a
//! small field with a private schoolbook multiplier that shares no code with `FpPoly`, by
//! group-theoretic structure (the 2-adic valuation of `q − 1`, the norm criterion for squareness),
//! or by an explicit factorisation.

use std::collections::HashSet;

use num_bigint::{BigInt, BigUint};
use shared_numfield::{
    FpPoly, FpPolyError, IntPoly, NonresidueSearch, NonresidueSearchLimits, PrimeModulus,
    ResidueField, ResidueFieldError, SqrtOutcome,
};

/// Largest prime below 2^64; 2^64 − 59 ≡ 5 (mod 8), so 2 is a nonresidue and x² − 2 is irreducible.
const P64: u64 = u64::MAX - 58;
/// Mersenne prime 2^61 − 1; ≡ 7 (mod 8), so −1 is a nonresidue and x² + 1 is irreducible.
const P61: u64 = (1 << 61) - 1;

const GENEROUS: NonresidueSearchLimits = NonresidueSearchLimits { max_candidates: 100_000 };
const ZERO_BUDGET: NonresidueSearchLimits = NonresidueSearchLimits { max_candidates: 0 };

fn pm(p: u64) -> PrimeModulus {
    PrimeModulus::new(p).unwrap()
}

/// Polynomial over F_p from signed least-significant-first coefficients.
fn fp(p: u64, coeffs: &[i64]) -> FpPoly {
    let m = pm(p);
    let reduced: Vec<u64> = coeffs.iter().map(|&c| m.reduce_i64(c)).collect();
    FpPoly::from_coeffs(m, &reduced)
}

fn field(p: u64, f: &[i64]) -> ResidueField {
    ResidueField::new(&fp(p, f)).unwrap()
}

// ─── independent oracle: schoolbook arithmetic on raw coefficient vectors ───────────────────────

/// Product in `F_p[x]/(f)` of two length-`d` vectors, for monic `f` of degree `d`, using only
/// `u128` arithmetic. Shares no code with `FpPoly`.
fn naive_mul(p: u64, f: &[u64], a: &[u64], b: &[u64]) -> Vec<u64> {
    let d = f.len() - 1;
    let p128 = u128::from(p);
    let mut prod = vec![0u64; 2 * d - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            let t = (u128::from(prod[i + j]) + u128::from(x) * u128::from(y)) % p128;
            prod[i + j] = u64::try_from(t).unwrap();
        }
    }
    for deg in (d..prod.len()).rev() {
        let c = u128::from(prod[deg]);
        for (i, &fc) in f.iter().enumerate() {
            let sub = c * u128::from(fc) % p128;
            let cur = u128::from(prod[deg - d + i]);
            prod[deg - d + i] = u64::try_from((cur + p128 - sub) % p128).unwrap();
        }
    }
    prod.truncate(d);
    prod
}

/// Monic coefficients of `f` (length `d + 1`), reduced into `[0, p)`.
fn monic_vec(p: u64, f: &[i64]) -> Vec<u64> {
    let poly = fp(p, f).to_monic().unwrap();
    poly.coeffs().to_vec()
}

/// All `p^d` coefficient vectors of length `d`, as base-`p` digit expansions of `0 .. p^d`.
fn all_elements(p: u64, d: usize) -> Vec<Vec<u64>> {
    let q = p.pow(u32::try_from(d).unwrap());
    (0..q)
        .map(|mut n| {
            (0..d)
                .map(|_| {
                    let digit = n % p;
                    n /= p;
                    digit
                })
                .collect()
        })
        .collect()
}

fn to_poly(p: u64, v: &[u64]) -> FpPoly {
    FpPoly::from_coeffs(pm(p), v)
}

/// Canonical length-`d` coefficient vector of a field element.
fn padded(poly: &FpPoly, d: usize) -> Vec<u64> {
    let mut v = poly.coeffs().to_vec();
    v.resize(d, 0);
    v
}

/// Cross-check `is_square`, `sqrt`, `inv` and the nonresidue search against enumeration of all
/// `p^d` elements of `F_p[x]/(f)`.
fn exhaustive_field_check(p: u64, f: &[i64]) {
    let fld = field(p, f);
    let d = fld.degree();
    let fv = monic_vec(p, f);
    let elems = all_elements(p, d);
    assert_eq!(fld.order(), &BigUint::from(elems.len()));

    let one: Vec<u64> = (0..d).map(|i| u64::from(i == 0)).collect();
    let squares: HashSet<Vec<u64>> = elems.iter().map(|e| naive_mul(p, &fv, e, e)).collect();
    let nonsquares: Vec<&Vec<u64>> = elems.iter().filter(|e| !squares.contains(*e)).collect();
    // Odd q: exactly the (q − 1)/2 nonzero squares and as many nonsquares; q even: all squares.
    if p == 2 {
        assert!(nonsquares.is_empty());
    } else {
        assert_eq!(nonsquares.len(), (elems.len() - 1) / 2);
    }

    for e in &elems {
        let a = to_poly(p, e);
        let is_sq = squares.contains(e);
        assert_eq!(fld.is_square(&a).unwrap(), is_sq, "is_square({e:?}) mod {f:?}/F_{p}");
        match fld.sqrt(&a, GENEROUS).unwrap() {
            SqrtOutcome::Root(r) => {
                assert!(is_sq, "root reported for nonsquare {e:?}");
                let rv = padded(&r, d);
                assert_eq!(&naive_mul(p, &fv, &rv, &rv), e, "root of {e:?}");
            }
            SqrtOutcome::NotSquare => assert!(!is_sq, "square {e:?} reported as nonsquare"),
            SqrtOutcome::ResourceExhausted { .. } => panic!("generous budget exhausted on {e:?}"),
        }
        if e.iter().all(|&c| c == 0) {
            assert_eq!(fld.inv(&a), Err(ResidueFieldError::DivisionByZero));
        } else {
            let inv = padded(&fld.inv(&a).unwrap(), d);
            assert_eq!(naive_mul(p, &fv, e, &inv), one, "inverse of {e:?}");
        }
    }

    // The search returns the first nonsquare in its documented order; every nonsquare it
    // returns is one of the brute-force nonsquares.
    match fld.find_nonresidue(GENEROUS).unwrap() {
        NonresidueSearch::Found { element, .. } => {
            assert!(nonsquares.contains(&&padded(&element, d)));
        }
        NonresidueSearch::Absent => assert_eq!(p, 2),
        NonresidueSearch::BudgetExhausted { .. } => panic!("generous budget exhausted"),
    }
}

// ─── field construction and validation ───────────────────────────────────────────────────────────

#[test]
fn kat_valid_fields_report_characteristic_degree_and_order() {
    // x² + 1 over F_3: no roots (0, 1, 2 give 1, 2, 2), so irreducible; 3² = 9 elements.
    let f = field(3, &[1, 0, 1]);
    assert_eq!(f.characteristic().get(), 3);
    assert_eq!(f.degree(), 2);
    assert_eq!(f.order(), &BigUint::from(9u32));
    assert_eq!(f.modulus_poly(), &fp(3, &[1, 0, 1]));
    // x³ − 2 over F_7: the cubes mod 7 are {0, 1, 6}, so 2 is not a cube and the cubic has no
    // root; a reducible cubic must have a root, so it is irreducible. 7³ = 343.
    assert_eq!(field(7, &[-2, 0, 0, 1]).order(), &BigUint::from(343u32));
    // A linear modulus gives F_p itself.
    let linear = field(17, &[0, 1]);
    assert_eq!((linear.degree(), linear.order()), (1, &BigUint::from(17u32)));
    // Characteristic 2: x⁴ + x + 1 is the standard primitive polynomial for GF(16).
    assert_eq!(field(2, &[1, 1, 0, 0, 1]).order(), &BigUint::from(16u32));
}

#[test]
fn kat_non_monic_modulus_is_scaled_to_monic() {
    // 2x² + 2 = 2·(x² + 1) over F_3 (2 is a unit), so both define the same field.
    let scaled = field(3, &[2, 0, 2]);
    assert_eq!(scaled.modulus_poly(), &fp(3, &[1, 0, 1]));
    assert_eq!(scaled, field(3, &[1, 0, 1]));
}

#[test]
fn kat_cardinality_exceeds_machine_words() {
    // p = 2^64 − 59 and d = 2: q = p² ≈ 2^128 cannot fit u64; (q − 1)/2 is needed as an exponent.
    let f = field(P64, &[-2, 0, 1]);
    let p = BigUint::from(P64);
    assert_eq!(f.order(), &(&p * &p));
    assert!(f.order() > &BigUint::from(u64::MAX));
    // p = 2^61 − 1 and d = 3 gives q ≈ 2^183, using the first non-cube c as x³ − c.
    let m = pm(P61);
    let non_cube = (2..).find(|&c| m.pow(c, (P61 - 1) / 3) != 1).unwrap();
    let cubic = field(P61, &[-i64::try_from(non_cube).unwrap(), 0, 0, 1]);
    let p = BigUint::from(P61);
    assert_eq!(cubic.order(), &(&p * &p * &p));
}

#[test]
fn kat_invalid_moduli_are_rejected_with_the_right_error() {
    let m5 = pm(5);
    // Zero polynomial: every element would be 0 in the quotient.
    assert_eq!(ResidueField::new(&FpPoly::zero(m5)), Err(ResidueFieldError::ZeroModulus));
    // A nonzero constant is a unit, so the quotient is the zero ring.
    assert_eq!(ResidueField::new(&fp(5, &[3])), Err(ResidueFieldError::ConstantModulus));
    // Reducible: x² − 1 = (x − 1)(x + 1) has roots ±1 in F_5.
    assert_eq!(ResidueField::new(&fp(5, &[-1, 0, 1])), Err(ResidueFieldError::ReducibleModulus));
    // Reducible without roots: (x² + 1)(x² + x + 2) = x⁴ + x³ + x + 2 over F_3 (product written
    // out; both quadratics have no root mod 3 so the quartic has none either).
    assert_eq!(
        ResidueField::new(&fp(3, &[2, 1, 0, 1, 1])),
        Err(ResidueFieldError::ReducibleModulus)
    );
    // Repeated factor: x² over F_3 and (x² + x + 1)² over F_2 are not squarefree.
    assert_eq!(ResidueField::new(&fp(3, &[0, 0, 1])), Err(ResidueFieldError::ReducibleModulus));
    assert_eq!(
        ResidueField::new(&fp(2, &[1, 0, 1, 0, 1])),
        Err(ResidueFieldError::ReducibleModulus)
    );
}

#[test]
fn kat_integer_polynomial_entry_checks_the_prime_and_the_reduction() {
    let ints = |c: &[i64]| IntPoly::from_coeffs(c.iter().map(|&x| BigInt::from(x)).collect());
    let x2_plus_1 = ints(&[1, 0, 1]);
    // 0, 1, 15 = 3·5 and 561 (a Carmichael number) are not prime moduli.
    for n in [0, 1, 15, 561] {
        assert_eq!(
            ResidueField::from_int_poly(&x2_plus_1, n),
            Err(ResidueFieldError::Poly(FpPolyError::NotPrime { candidate: n }))
        );
    }
    assert_eq!(ResidueField::from_int_poly(&x2_plus_1, 3).unwrap(), field(3, &[1, 0, 1]));
    // x² + 1 = (x − 2)(x + 2) mod 5 is reducible.
    assert_eq!(
        ResidueField::from_int_poly(&x2_plus_1, 5),
        Err(ResidueFieldError::ReducibleModulus)
    );
    // 3x² + x + 1 mod 3 drops to x + 1: a degree-1 field, because the leading coefficient dies.
    let dropped = ResidueField::from_int_poly(&ints(&[1, 1, 3]), 3).unwrap();
    assert_eq!(dropped.degree(), 1);
    // A constant divisible by p reduces to the zero polynomial.
    assert_eq!(ResidueField::from_int_poly(&ints(&[3]), 3), Err(ResidueFieldError::ZeroModulus));
}

#[test]
fn kat_error_display_and_source() {
    let not_prime = ResidueFieldError::Poly(FpPolyError::NotPrime { candidate: 15 });
    assert_eq!(not_prime.to_string(), "15 is not prime");
    let source = std::error::Error::source(&not_prime).unwrap();
    assert_eq!(source.to_string(), "15 is not prime");
    for e in [
        ResidueFieldError::ZeroModulus,
        ResidueFieldError::ConstantModulus,
        ResidueFieldError::ReducibleModulus,
        ResidueFieldError::DivisionByZero,
        ResidueFieldError::InvariantViolation,
    ] {
        assert!(!e.to_string().is_empty());
        assert!(std::error::Error::source(&e).is_none());
    }
}

#[test]
fn kat_elements_over_another_characteristic_are_rejected() {
    let f = field(3, &[1, 0, 1]);
    let foreign = fp(5, &[1, 1]);
    let own = fp(3, &[1, 1]);
    let mismatch = |r: Result<_, ResidueFieldError>| {
        matches!(r, Err(ResidueFieldError::Poly(FpPolyError::ModulusMismatch { .. })))
    };
    assert!(mismatch(f.element(&foreign)));
    assert!(mismatch(f.add(&foreign, &own)));
    assert!(mismatch(f.add(&own, &foreign)));
    assert!(mismatch(f.sub(&own, &foreign)));
    assert!(mismatch(f.neg(&foreign)));
    assert!(mismatch(f.mul(&foreign, &own)));
    assert!(mismatch(f.pow(&foreign, &BigUint::from(2u32))));
    assert!(mismatch(f.inv(&foreign)));
    assert!(matches!(
        f.is_square(&foreign),
        Err(ResidueFieldError::Poly(FpPolyError::ModulusMismatch { left: 5, right: 3 }))
    ));
    assert!(matches!(
        f.sqrt(&foreign, GENEROUS),
        Err(ResidueFieldError::Poly(FpPolyError::ModulusMismatch { .. }))
    ));
}

// ─── arithmetic and inverses ─────────────────────────────────────────────────────────────────────

#[test]
fn kat_arithmetic_reduces_modulo_the_modulus() {
    let f = field(3, &[1, 0, 1]); // x² = −1
    let x = fp(3, &[0, 1]);
    // x · x = −1 = 2.
    assert_eq!(f.mul(&x, &x).unwrap(), fp(3, &[2]));
    // (1 + x)² = 1 + 2x + x² = 2x.
    let u = fp(3, &[1, 1]);
    assert_eq!(f.mul(&u, &u).unwrap(), fp(3, &[0, 2]));
    // Unreduced spellings denote the same element: x² + x ≡ x − 1.
    assert_eq!(f.element(&fp(3, &[0, 1, 1])).unwrap(), fp(3, &[-1, 1]));
    // add/sub/neg: (1 + x) + (2 + 2x) = 0; (1 + x) − (1 + x) = 0; −(1 + x) = 2 + 2x.
    assert_eq!(f.add(&u, &fp(3, &[2, 2])).unwrap(), f.zero());
    assert_eq!(f.sub(&u, &u).unwrap(), f.zero());
    assert_eq!(f.neg(&u).unwrap(), fp(3, &[2, 2]));
    // pow: the unit group has order 8, so u⁸ = 1 and u⁰ = 1.
    assert_eq!(f.pow(&u, &BigUint::from(8u32)).unwrap(), f.one());
    assert_eq!(f.pow(&u, &BigUint::from(0u32)).unwrap(), f.one());
    assert_ne!(f.zero(), f.one());
}

#[test]
fn kat_inverse_by_multiplication() {
    // F_5[x]/(x² − 2): (1 + x)(a + bx) = a + 2b + (a + b)x = 1 forces a = −b, and a + 2b = b = 1,
    // so the inverse of 1 + x is −1 + x = 4 + x. Checked by hand, then by multiplication.
    let f = field(5, &[-2, 0, 1]);
    let u = fp(5, &[1, 1]);
    let inv = f.inv(&u).unwrap();
    assert_eq!(inv, fp(5, &[4, 1]));
    assert_eq!(f.mul(&u, &inv).unwrap(), f.one());
    // The inverse of a unit that is only an unreduced spelling: f ≡ 0, f + 1 ≡ 1.
    assert_eq!(f.inv(&fp(5, &[-1, 0, 1])).unwrap(), f.one());
    // Over a large prime, inverses survive 128-bit intermediate products.
    let big = field(P64, &[-2, 0, 1]);
    let w = fp(P64, &[7, -3]);
    assert_eq!(big.mul(&w, &big.inv(&w).unwrap()).unwrap(), big.one());
}

#[test]
fn kat_zero_has_no_inverse_but_is_a_square_with_root_zero() {
    for (p, coeffs) in [(3u64, vec![1, 0, 1]), (7, vec![-2, 0, 0, 1]), (2, vec![1, 1, 0, 0, 1])] {
        let f = field(p, &coeffs);
        assert_eq!(f.inv(&f.zero()), Err(ResidueFieldError::DivisionByZero));
        // The modulus itself is an unreduced zero.
        assert_eq!(f.inv(f.modulus_poly()), Err(ResidueFieldError::DivisionByZero));
        assert!(f.is_square(&f.zero()).unwrap());
        // Zero needs no nonresidue, so even a zero budget finds its root.
        assert_eq!(f.sqrt(&f.zero(), ZERO_BUDGET).unwrap(), SqrtOutcome::Root(f.zero()));
        assert_eq!(f.sqrt(f.modulus_poly(), ZERO_BUDGET).unwrap(), SqrtOutcome::Root(f.zero()));
    }
}

// ─── exhaustive cross-checks against enumeration ─────────────────────────────────────────────────

#[test]
fn kat_exhaustive_even_degree_fields() {
    // Quadratic extensions where constants are squares, so the nonresidue search leaves F_p.
    // Moduli: x²+1/F_3 (q−1 = 8 = 2³, large 2-adic part), x²−2/F_5 (24 = 2³·3), x²+1/F_7
    // (48 = 2⁴·3), x²+1/F_11, x²−2/F_13 (168 = 2³·21), x²−3/F_17 (288 = 2⁵·9); each is
    // irreducible because the constant term is a nonresidue (−1 mod 3, 7, 11; 2 mod 5, 13; 3
    // mod 17, from the square tables {1,4,9,16,8,2,15,13}).
    for (p, f) in [
        (3, [1, 0, 1]),
        (5, [-2, 0, 1]),
        (7, [1, 0, 1]),
        (11, [1, 0, 1]),
        (13, [-2, 0, 1]),
        (17, [-3, 0, 1]),
    ] {
        exhaustive_field_check(p, &f);
    }
}

#[test]
fn kat_exhaustive_odd_degree_and_prime_fields() {
    // Odd degree: x³ − x − 1 over F_3 (q − 1 = 26 = 2·13) and x³ − 2 over F_7 (q − 1 = 342). The
    // prime fields F_17 = F_17[x]/(x) (16 = 2⁴), F_257 (256 = 2⁸) and F_97 (96 = 2⁵·3) give long
    // Tonelli–Shanks loops because the odd part of q − 1 is 1, 1 and 3.
    exhaustive_field_check(3, &[-1, -1, 0, 1]);
    exhaustive_field_check(7, &[-2, 0, 0, 1]);
    for p in [17, 97, 257] {
        exhaustive_field_check(p, &[0, 1]);
    }
}

#[test]
fn kat_exhaustive_characteristic_two() {
    // GF(2), GF(8) (x³+x+1) and GF(16) (x⁴+x+1): every element is a square because Frobenius is
    // a bijection of a finite field of characteristic 2, and the root is unique.
    exhaustive_field_check(2, &[0, 1]);
    exhaustive_field_check(2, &[1, 1, 0, 1]);
    exhaustive_field_check(2, &[1, 1, 0, 0, 1]);
}

// ─── square roots: structure and large cases ─────────────────────────────────────────────────────

#[test]
fn kat_roots_of_degree_two_squares_are_plus_or_minus_the_witness() {
    // Squares of known witnesses: the only roots of w² in a field are ±w.
    let cases: [(u64, &[i64], &[i64]); 3] = [
        (3, &[1, 0, 1], &[1, 1]),          // F_9, w = 1 + x
        (5, &[-2, 0, 1], &[3, 4]),         // F_25, w = 3 + 4x
        (P64, &[-2, 0, 1], &[1_000_003, 5]), // F_{p²} with p = 2^64 − 59
    ];
    for (p, modulus, witness) in cases {
        let f = field(p, modulus);
        let w = fp(p, witness);
        let a = f.mul(&w, &w).unwrap();
        assert!(f.is_square(&a).unwrap());
        let SqrtOutcome::Root(r) = f.sqrt(&a, GENEROUS).unwrap() else {
            panic!("square of a witness must have a root");
        };
        assert_eq!(f.mul(&r, &r).unwrap(), a);
        assert!(r == f.element(&w).unwrap() || r == f.neg(&w).unwrap());
    }
}

/// Whether `u + v·x` is a square in `F_p[x]/(x² − n)`: its norm `u² − n·v²` is a square in
/// `F_p`, because `a^((q−1)/2) = N(a)^((p−1)/2)`. Independent of the field code.
fn quadratic_element_is_square(p: u64, n: u64, u: u64, v: u64) -> bool {
    let m = pm(p);
    let norm = m.sub(m.mul(u, u), m.mul(n, m.mul(v, v)));
    assert_ne!(norm, 0, "irreducible modulus has no zero norm");
    m.pow(norm, (p - 1) / 2) == 1
}

#[test]
fn kat_deep_two_adic_valuation_over_a_machine_prime() {
    // p = 2^61 − 1, x² + 1: q − 1 = (p − 1)(p + 1) = 2·(2^60 − 1)·2^61, so the 2-adic valuation of
    // q − 1 is 62 and Tonelli–Shanks runs up to 61 reduction passes on 128-bit-scale elements.
    let f = field(P61, &[1, 0, 1]); // x² = −1, so n = p − 1
    let n = P61 - 1;
    let mut squares = 0;
    let mut nonsquares = 0;
    for (u, v) in [(1, 1), (2, 1), (3, 1), (5, 2), (7, 3), (11, 4), (123, 456), (P61 - 2, 17)] {
        let a = fp(P61, &[i64::try_from(u).unwrap(), i64::try_from(v).unwrap()]);
        let expect_square = quadratic_element_is_square(P61, n, u, v);
        assert_eq!(f.is_square(&a).unwrap(), expect_square, "norm criterion for {u}+{v}x");
        match f.sqrt(&a, GENEROUS).unwrap() {
            SqrtOutcome::Root(r) => {
                squares += 1;
                assert!(expect_square);
                assert_eq!(f.mul(&r, &r).unwrap(), a);
            }
            SqrtOutcome::NotSquare => {
                nonsquares += 1;
                assert!(!expect_square);
            }
            SqrtOutcome::ResourceExhausted { .. } => panic!("generous budget exhausted"),
        }
    }
    // The sample must exercise both outcomes, and squares of fixed witnesses must be roots.
    assert!(squares > 0 && nonsquares > 0, "{squares} squares, {nonsquares} nonsquares");
    let w = fp(P61, &[1_234_567, 7_654_321]);
    let w2 = f.mul(&w, &w).unwrap();
    let SqrtOutcome::Root(r) = f.sqrt(&w2, GENEROUS).unwrap() else { panic!("w² has a root") };
    assert!(r == w || r == f.neg(&w).unwrap());
}

#[test]
fn kat_odd_degree_over_a_machine_prime() {
    // p = 2^61 − 1, x³ − c with c a non-cube: q = p³ ≈ 2^183. In an odd-degree extension a scalar
    // is a square iff it is a square in F_p (c^((q−1)/2) = Legendre(c)^(1+p+p²), and 1+p+p² is
    // odd), so the first F_p nonresidue stays a nonsquare.
    let m = pm(P61);
    let non_cube = (2..).find(|&c| m.pow(c, (P61 - 1) / 3) != 1).unwrap();
    let f = field(P61, &[-i64::try_from(non_cube).unwrap(), 0, 0, 1]);
    let w = fp(P61, &[3, -5, 11]);
    let a = f.mul(&w, &w).unwrap();
    let SqrtOutcome::Root(r) = f.sqrt(&a, GENEROUS).unwrap() else { panic!("w² has a root") };
    assert_eq!(f.mul(&r, &r).unwrap(), a);
    // The least nonresidue of F_p is found among the constants and is a nonsquare of F_q.
    let least = (2..).find(|&c| m.pow(c, (P61 - 1) / 2) != 1).unwrap();
    assert!(!f.is_square(&FpPoly::constant(m, least)).unwrap());
    assert_eq!(f.sqrt(&FpPoly::constant(m, least), GENEROUS).unwrap(), SqrtOutcome::NotSquare);
    // Constants of F_p that are squares in F_p stay squares of F_q.
    assert!(f.is_square(&FpPoly::constant(m, 4)).unwrap());
}

#[test]
fn kat_large_prime_field_roots_and_nonsquares() {
    // F_p for p = 2^64 − 59 (q − 1 = 2²·odd, so s = 2) with modulus x − 5: elements are constants.
    let f = field(P64, &[-5, 1]);
    let m = pm(P64);
    for a in [4u64, 9, 1_000_000_007, P64 - 1] {
        let expect_square = m.pow(a, (P64 - 1) / 2) == 1;
        match f.sqrt(&FpPoly::constant(m, a), GENEROUS).unwrap() {
            SqrtOutcome::Root(r) => {
                assert!(expect_square);
                assert_eq!(f.mul(&r, &r).unwrap(), FpPoly::constant(m, a));
            }
            SqrtOutcome::NotSquare => assert!(!expect_square),
            SqrtOutcome::ResourceExhausted { .. } => panic!("generous budget exhausted"),
        }
    }
    // 2 is a nonresidue mod p ≡ 5 (mod 8).
    assert_eq!(f.sqrt(&FpPoly::constant(m, 2), GENEROUS).unwrap(), SqrtOutcome::NotSquare);
    // p ≡ 1 (mod 4), so −1 is a square (its roots are the order-4 elements).
    let SqrtOutcome::Root(i) = f.sqrt(&FpPoly::constant(m, P64 - 1), GENEROUS).unwrap() else {
        panic!("−1 is a square mod a prime ≡ 1 mod 4")
    };
    assert_eq!(f.mul(&i, &i).unwrap(), FpPoly::constant(m, P64 - 1));
}

#[test]
fn kat_unreduced_inputs_are_accepted() {
    // 2 + x² ≡ 1 (mod x² + 1) over F_3, so the unreduced spelling has the roots ±1.
    let f = field(3, &[1, 0, 1]);
    let SqrtOutcome::Root(r) = f.sqrt(&fp(3, &[2, 0, 1]), GENEROUS).unwrap() else {
        panic!("1 is a square")
    };
    assert!(r == f.one() || r == f.neg(&f.one()).unwrap());
    assert!(f.is_square(&fp(3, &[2, 0, 1])).unwrap());
}

// ─── characteristic 2 ────────────────────────────────────────────────────────────────────────────

#[test]
fn kat_characteristic_two_has_unique_roots_and_no_nonresidues() {
    // GF(16) with x⁴ + x + 1: the root of a is a^(q/2) = a⁸, unique since −r = r.
    let f = field(2, &[1, 1, 0, 0, 1]);
    assert_eq!(f.find_nonresidue(GENEROUS).unwrap(), NonresidueSearch::Absent);
    assert_eq!(f.find_nonresidue(ZERO_BUDGET).unwrap(), NonresidueSearch::Absent);
    let x = fp(2, &[0, 1]);
    // x² has root x; zero budget is irrelevant because no search is performed.
    let x2 = f.mul(&x, &x).unwrap();
    assert_eq!(f.sqrt(&x2, ZERO_BUDGET).unwrap(), SqrtOutcome::Root(x.clone()));
    assert_eq!(f.neg(&x).unwrap(), x);
    // Every element is a square: x⁴ = x + 1 gives x⁸ = x² + 1, and (x² + 1)² = x⁴ + 1 = x.
    assert_eq!(f.sqrt(&x, ZERO_BUDGET).unwrap(), SqrtOutcome::Root(fp(2, &[1, 0, 1])));
}

// ─── resource accounting ─────────────────────────────────────────────────────────────────────────

#[test]
fn kat_search_order_and_budget_boundary() {
    // F_17 (odd degree): the constants 1, 2, 3, … are tried in order; the squares mod 17 are
    // {1, 2, 4, 8, 9, 13, 15, 16}, so 3 is the first nonresidue, at candidate 3.
    let f = field(17, &[0, 1]);
    let found = |limit| f.find_nonresidue(NonresidueSearchLimits { max_candidates: limit });
    assert_eq!(found(0).unwrap(), NonresidueSearch::BudgetExhausted { candidates_tried: 0 });
    assert_eq!(found(2).unwrap(), NonresidueSearch::BudgetExhausted { candidates_tried: 2 });
    assert_eq!(
        found(3).unwrap(),
        NonresidueSearch::Found { element: fp(17, &[3]), candidates_tried: 3 }
    );
    assert_eq!(
        found(1_000).unwrap(),
        NonresidueSearch::Found { element: fp(17, &[3]), candidates_tried: 3 }
    );
    // F_7[x]/(x³ − 2) is an odd-degree field too: squares mod 7 are {1, 2, 4}, so 3 again.
    let g = field(7, &[-2, 0, 0, 1]);
    assert_eq!(
        g.find_nonresidue(GENEROUS).unwrap(),
        NonresidueSearch::Found { element: fp(7, &[3]), candidates_tried: 3 }
    );
}

#[test]
fn kat_even_degree_search_skips_constants() {
    // F_9 = F_3[x]/(x² + 1): the constants 1 and 2 = x² are squares, so the candidates are the
    // monic linear polynomials x + 0, x + 1, x + 2 in that order. The expected position of the
    // first nonsquare among them comes from enumerating all squares of F_9 by brute force.
    let f = field(3, &[1, 0, 1]);
    let fv = monic_vec(3, &[1, 0, 1]);
    let squares: HashSet<Vec<u64>> =
        all_elements(3, 2).iter().map(|e| naive_mul(3, &fv, e, e)).collect();
    let expected = (0..3u64)
        .map(|c| vec![c, 1])
        .position(|v| !squares.contains(&v))
        .expect("a monic nonsquare exists");
    let k = u64::try_from(expected).unwrap() + 1;
    let limits = |n| NonresidueSearchLimits { max_candidates: n };
    assert_eq!(
        f.find_nonresidue(limits(k)).unwrap(),
        NonresidueSearch::Found { element: to_poly(3, &[k - 1, 1]), candidates_tried: k }
    );
    assert_eq!(
        f.find_nonresidue(limits(k - 1)).unwrap(),
        NonresidueSearch::BudgetExhausted { candidates_tried: k - 1 }
    );
    // F_{p²} with p = 2^61 − 1: constants (p − 1 of them) are all squares, yet the search finds a
    // nonsquare among x + c within a few candidates, verified by the norm criterion.
    let big = field(P61, &[1, 0, 1]);
    let NonresidueSearch::Found { element, candidates_tried } = big.find_nonresidue(GENEROUS).unwrap()
    else {
        panic!("a nonsquare exists among x + c")
    };
    assert!(candidates_tried < 100);
    let (u, v) = (element.coeffs()[0], element.coeffs()[1]);
    assert_eq!(element.degree(), Some(1));
    assert!(!quadratic_element_is_square(P61, P61 - 1, u, v));
}

#[test]
fn kat_budget_exhaustion_is_distinct_from_nonsquare() {
    // F_5[x]/(x² − 2): q − 1 = 24 = 2³·3. A nonzero square with a^3 ≠ 1 forces a nonresidue
    // search. With budget 0 that is a resource outcome. Nonsquares and zero never need one.
    let f = field(5, &[-2, 0, 1]);
    let fv = monic_vec(5, &[-2, 0, 1]);
    let elems = all_elements(5, 2);
    let squares: HashSet<Vec<u64>> = elems.iter().map(|e| naive_mul(5, &fv, e, e)).collect();
    let mut exhausted = 0;
    for e in &elems {
        let a = to_poly(5, e);
        let outcome = f.sqrt(&a, ZERO_BUDGET).unwrap();
        let is_zero = e.iter().all(|&c| c == 0);
        if is_zero {
            assert_eq!(outcome, SqrtOutcome::Root(f.zero()));
        } else if !squares.contains(e) {
            // A proved nonsquare, independent of the budget.
            assert_eq!(outcome, SqrtOutcome::NotSquare);
        } else {
            match outcome {
                SqrtOutcome::Root(r) => assert_eq!(f.mul(&r, &r).unwrap(), a),
                SqrtOutcome::ResourceExhausted { candidates_tried } => {
                    exhausted += 1;
                    assert_eq!(candidates_tried, 0);
                    // A larger budget resolves exactly the exhausted cases.
                    assert!(matches!(f.sqrt(&a, GENEROUS).unwrap(), SqrtOutcome::Root(_)));
                }
                SqrtOutcome::NotSquare => panic!("square {e:?} reported as nonsquare"),
            }
        }
    }
    // The 12 nonzero squares form the subgroup of order 12. Those with a³ = 1 (the 3 cube roots
    // of unity, which lie in it) need no search, so the other 12 − 3 = 9 are exhausted.
    assert_eq!(exhausted, 9);
}

#[test]
fn kat_too_small_budget_reports_candidates_spent() {
    // F_17: a = 4 has a^1 ≠ 1 (odd part t = 1), so a nonresidue is needed and the first one is the
    // third candidate. Budgets 1 and 2 are exhausted with that many candidates tried; 3 succeeds.
    let f = field(17, &[0, 1]);
    let four = FpPoly::constant(pm(17), 4);
    for n in [0u64, 1, 2] {
        assert_eq!(
            f.sqrt(&four, NonresidueSearchLimits { max_candidates: n }).unwrap(),
            SqrtOutcome::ResourceExhausted { candidates_tried: n }
        );
    }
    let SqrtOutcome::Root(r) = f.sqrt(&four, NonresidueSearchLimits { max_candidates: 3 }).unwrap()
    else {
        panic!("budget 3 reaches the first nonresidue")
    };
    assert!(r == FpPoly::constant(pm(17), 2) || r == FpPoly::constant(pm(17), 15));
}

#[test]
fn kat_three_mod_four_needs_no_nonresidue() {
    // q = 3 mod 4 means q − 1 = 2·t with t odd, so for a square a, a^t = a^((q−1)/2) = 1 and
    // Tonelli–Shanks returns a^((t+1)/2) immediately: even a zero budget succeeds.
    // F_7 (p = 3 mod 4, d = 1) and F_3[x]/(x³ − x − 1) (q = 27 = 3 mod 4).
    for (p, coeffs) in [(7u64, vec![0, 1]), (3, vec![-1, -1, 0, 1])] {
        let f = field(p, &coeffs);
        for a in all_elements(p, f.degree()) {
            let a = to_poly(p, &a);
            match f.sqrt(&a, ZERO_BUDGET).unwrap() {
                SqrtOutcome::Root(r) => assert_eq!(f.mul(&r, &r).unwrap(), a),
                SqrtOutcome::NotSquare => assert!(!f.is_square(&a).unwrap()),
                SqrtOutcome::ResourceExhausted { .. } => panic!("no search is needed"),
            }
        }
    }
}
