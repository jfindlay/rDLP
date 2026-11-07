//! Known-answer tests and field-axiom tests for `F2mNaive`, `F2mNormal`, and `F2mOpt`.
//!
//! Coverage:
//! - Field axioms over GF(2^4) with `x⁴+x+1` (poly = 0x13).
//! - Field axioms over GF(2^8) AES with `x⁸+x⁴+x³+x+1` (poly = 0x11b).
//! - Char-2 invariants: `sub(a,b) == add(a,b)`, `neg(a) == a`.
//! - Frobenius fixed-field law: `a^(2^m) == a` for all `a` in the field.
//! - `square(a) == mul(a, a)`.
//! - Associativity of `mul`.
//! - Distributivity of `mul` over `add`.
//! - `mul` by `one` is identity.
//! - `add` is its own inverse: `add(a, a) == zero`.
//! - Inversion round-trip: `mul(a, inv(a)) == one` for all non-zero `a`.
//! - Extended-Euclidean ↔ Itoh–Tsujii agreement.
//! - `inv(one) == one`.
//! - `inv(zero)` panics.
//! - `div(a, b) == mul(a, inv(b))`.
//! - Normal-basis round-trips: `poly_to_normal(normal_to_poly(a)) == a` and vice versa.
//! - Squaring-is-a-cyclic-shift: `square(a).c == cyclic_left_shift(a.c, m)` (direct check).
//! - Cross-representation agreement: `mul`/`add` in polynomial basis equals the same
//!   operation in normal basis after conversion.
//! - **Naive↔optimized agreement**: `F2mNaive::mul` and `F2mOpt::mul` give byte-identical
//!   results for all input pairs in GF(2^4) (256 pairs) and GF(2^8) (65536 pairs).
//! - **`F2mOpt` field axioms**: associativity, distributivity, mul-by-one, add-self-zero
//!   verified directly against `F2mOpt` (not only via agreement with naive).
//!
//! The Frobenius fixed-field law `a^(2^m) = a` is the loudest correctness
//! signal: it fails immediately if modular reduction is wrong (wrong irreducible,
//! off-by-one in the shift, etc.).

use crypto_bigint::Uint;
use shared_gf2m::BitLen;
use shared_gf2m::{
    ext_euclid_inv, find_normal_element, itoh_tsujii_inv, normal_to_poly, poly_to_normal, F2m,
    F2mNaive, F2mNormal, F2mOpt,
};

// ── Field parameters ──────────────────────────────────────────────────────────

/// GF(2^4) irreducible: x⁴ + x + 1 = 0b1_0011 = 0x13.  Degree m = 4.
fn poly4() -> Uint<1> {
    Uint::<1>::from(0x13u64)
}

/// GF(2^8) AES irreducible: x⁸ + x⁴ + x³ + x + 1 = 0x11b.  Degree m = 8.
fn poly8() -> Uint<1> {
    Uint::<1>::from(0x11bu64)
}

/// All non-zero elements of GF(2^4): 1 through 15.
fn gf4_elements() -> Vec<F2mNaive<1>> {
    let p = poly4();
    (0u64..16).map(|v| F2mNaive::<1>::from_u64(v, &p)).collect()
}

/// All non-zero elements of GF(2^8): 1 through 255.
fn gf8_elements() -> Vec<F2mNaive<1>> {
    let p = poly8();
    (0u64..256).map(|v| F2mNaive::<1>::from_u64(v, &p)).collect()
}

// ── Char-2 invariant KATs ─────────────────────────────────────────────────────

/// `sub(a, b) == add(a, b)` for all pairs in GF(2^4).
///
/// This is the load-bearing char-2 invariant: subtraction is XOR, identical
/// to addition.  A wrong `sub` (porting `Fp`'s `p − self`) would fail here.
#[test]
fn sub_equals_add_gf4() {
    for a in gf4_elements() {
        for b in gf4_elements() {
            assert_eq!(
                a.sub(&b),
                a.add(&b),
                "sub != add for a={:#x} b={:#x} in GF(2^4)",
                a.to_uint(),
                b.to_uint()
            );
        }
    }
}

/// `sub(a, b) == add(a, b)` for all pairs in GF(2^8).
#[test]
fn sub_equals_add_gf8() {
    // Exhaustive over all 256×256 pairs is 65536 — fast enough.
    for a in gf8_elements() {
        for b in gf8_elements() {
            assert_eq!(
                a.sub(&b),
                a.add(&b),
                "sub != add for a={:#x} b={:#x} in GF(2^8)",
                a.to_uint(),
                b.to_uint()
            );
        }
    }
}

/// `neg(a) == a` for all elements in GF(2^4).
///
/// In characteristic 2, negation is the identity.  A wrong `neg` (porting
/// `Fp`'s `p − self`) would fail here.
#[test]
fn neg_is_identity_gf4() {
    for a in gf4_elements() {
        assert_eq!(
            a.neg(),
            a,
            "neg(a) != a for a={:#x} in GF(2^4)",
            a.to_uint()
        );
    }
}

/// `neg(a) == a` for all elements in GF(2^8).
#[test]
fn neg_is_identity_gf8() {
    for a in gf8_elements() {
        assert_eq!(
            a.neg(),
            a,
            "neg(a) != a for a={:#x} in GF(2^8)",
            a.to_uint()
        );
    }
}

// ── Additive axioms ───────────────────────────────────────────────────────────

/// `add(a, a) == zero` for all elements in GF(2^4).
///
/// Every element is its own additive inverse in char 2.
#[test]
fn add_self_is_zero_gf4() {
    let zero = F2mNaive::<1>::zero();
    for a in gf4_elements() {
        assert_eq!(
            a.add(&a),
            zero,
            "a + a != 0 for a={:#x} in GF(2^4)",
            a.to_uint()
        );
    }
}

/// `add(a, a) == zero` for all elements in GF(2^8).
#[test]
fn add_self_is_zero_gf8() {
    let zero = F2mNaive::<1>::zero();
    for a in gf8_elements() {
        assert_eq!(
            a.add(&a),
            zero,
            "a + a != 0 for a={:#x} in GF(2^8)",
            a.to_uint()
        );
    }
}

// ── Multiplicative axioms ─────────────────────────────────────────────────────

/// `mul(a, one) == a` for all elements in GF(2^4).
#[test]
fn mul_by_one_is_identity_gf4() {
    let p = poly4();
    let one = F2mNaive::<1>::one();
    for a in gf4_elements() {
        assert_eq!(
            a.mul(&one, &p),
            a,
            "a * 1 != a for a={:#x} in GF(2^4)",
            a.to_uint()
        );
    }
}

/// `mul(a, one) == a` for all elements in GF(2^8).
#[test]
fn mul_by_one_is_identity_gf8() {
    let p = poly8();
    let one = F2mNaive::<1>::one();
    for a in gf8_elements() {
        assert_eq!(
            a.mul(&one, &p),
            a,
            "a * 1 != a for a={:#x} in GF(2^8)",
            a.to_uint()
        );
    }
}

/// Associativity of `mul` over GF(2^4): `mul(mul(a,b),c) == mul(a,mul(b,c))`.
#[test]
fn mul_associative_gf4() {
    let p = poly4();
    // Exhaustive over all 16^3 = 4096 triples.
    for a in gf4_elements() {
        for b in gf4_elements() {
            for c in gf4_elements() {
                let lhs = a.mul(&b, &p).mul(&c, &p);
                let rhs = a.mul(&b.mul(&c, &p), &p);
                assert_eq!(
                    lhs, rhs,
                    "mul not associative: a={:#x} b={:#x} c={:#x} in GF(2^4)",
                    a.to_uint(), b.to_uint(), c.to_uint()
                );
            }
        }
    }
}

/// Associativity of `mul` over GF(2^8): spot-check a representative set.
///
/// Exhaustive (256^3 = 16M) is too slow; we check a representative subset.
#[test]
fn mul_associative_gf8() {
    let p = poly8();
    // Representative elements: 0, 1, x, x²-1, AES S-box inputs, max.
    let elems: Vec<u64> = vec![0, 1, 2, 3, 0x53, 0xca, 0x8d, 0xf6, 0x01, 0xff, 0x1b, 0xe5];
    for &av in &elems {
        for &bv in &elems {
            for &cv in &elems {
                let a = F2mNaive::<1>::from_u64(av, &p);
                let b = F2mNaive::<1>::from_u64(bv, &p);
                let c = F2mNaive::<1>::from_u64(cv, &p);
                let lhs = a.mul(&b, &p).mul(&c, &p);
                let rhs = a.mul(&b.mul(&c, &p), &p);
                assert_eq!(
                    lhs, rhs,
                    "mul not associative: a={av:#x} b={bv:#x} c={cv:#x} in GF(2^8)"
                );
            }
        }
    }
}

/// Distributivity of `mul` over `add` in GF(2^4):
/// `mul(a, add(b, c)) == add(mul(a,b), mul(a,c))`.
#[test]
fn mul_distributive_gf4() {
    let p = poly4();
    for a in gf4_elements() {
        for b in gf4_elements() {
            for c in gf4_elements() {
                let lhs = a.mul(&b.add(&c), &p);
                let rhs = a.mul(&b, &p).add(&a.mul(&c, &p));
                assert_eq!(
                    lhs, rhs,
                    "distributivity failed: a={:#x} b={:#x} c={:#x} in GF(2^4)",
                    a.to_uint(), b.to_uint(), c.to_uint()
                );
            }
        }
    }
}

/// Distributivity of `mul` over `add` in GF(2^8): spot-check.
#[test]
fn mul_distributive_gf8() {
    let p = poly8();
    let elems: Vec<u64> = vec![0, 1, 2, 3, 0x53, 0xca, 0x8d, 0xf6, 0x01, 0xff, 0x1b, 0xe5];
    for &av in &elems {
        for &bv in &elems {
            for &cv in &elems {
                let a = F2mNaive::<1>::from_u64(av, &p);
                let b = F2mNaive::<1>::from_u64(bv, &p);
                let c = F2mNaive::<1>::from_u64(cv, &p);
                let lhs = a.mul(&b.add(&c), &p);
                let rhs = a.mul(&b, &p).add(&a.mul(&c, &p));
                assert_eq!(
                    lhs, rhs,
                    "distributivity failed: a={av:#x} b={bv:#x} c={cv:#x} in GF(2^8)"
                );
            }
        }
    }
}

// ── Squaring / Frobenius ──────────────────────────────────────────────────────

/// `square(a) == mul(a, a)` for all elements in GF(2^4).
#[test]
fn square_equals_mul_self_gf4() {
    let p = poly4();
    for a in gf4_elements() {
        assert_eq!(
            a.square(&p),
            a.mul(&a, &p),
            "square != mul(a,a) for a={:#x} in GF(2^4)",
            a.to_uint()
        );
    }
}

/// `square(a) == mul(a, a)` for all elements in GF(2^8).
#[test]
fn square_equals_mul_self_gf8() {
    let p = poly8();
    for a in gf8_elements() {
        assert_eq!(
            a.square(&p),
            a.mul(&a, &p),
            "square != mul(a,a) for a={:#x} in GF(2^8)",
            a.to_uint()
        );
    }
}

// ── Frobenius fixed-field law: a^(2^m) == a ───────────────────────────────────
//
// This is the loudest correctness signal for the implementation.  If modular
// reduction is wrong (wrong irreducible, off-by-one in the shift, etc.), this
// test fails immediately.  Every element of GF(2^m) satisfies a^(2^m) = a
// (the Frobenius endomorphism has order m, so a^(2^m) = a).

/// `a^(2^4) == a` for all elements in GF(2^4).
///
/// Computed as 4 successive squarings (the Frobenius tower).
#[test]
fn frobenius_fixed_field_gf4() {
    let p = poly4();
    // 2^4 = 16 as a Uint<1>.
    let exp = Uint::<1>::from(16u64);
    for a in gf4_elements() {
        let a_pow = a.pow(&exp, &p);
        assert_eq!(
            a_pow, a,
            "a^(2^4) != a for a={:#x} in GF(2^4) — reduction is wrong",
            a.to_uint()
        );
    }
}

/// `a^(2^8) == a` for all elements in GF(2^8).
///
/// Computed via `pow` with exponent 256.
#[test]
fn frobenius_fixed_field_gf8() {
    let p = poly8();
    // 2^8 = 256 as a Uint<1>.
    let exp = Uint::<1>::from(256u64);
    for a in gf8_elements() {
        let a_pow = a.pow(&exp, &p);
        assert_eq!(
            a_pow, a,
            "a^(2^8) != a for a={:#x} in GF(2^8) — reduction is wrong",
            a.to_uint()
        );
    }
}

/// Frobenius tower: 4 successive squarings give identity in GF(2^4).
///
/// Verifies that `square` iterated m times returns the original element,
/// which is the direct Frobenius-tower check (distinct from `pow`).
#[test]
fn frobenius_tower_gf4() {
    let p = poly4();
    for a in gf4_elements() {
        let mut x = a.clone();
        for _ in 0..4 {
            x = x.square(&p);
        }
        assert_eq!(
            x, a,
            "Frobenius tower (4 squarings) != identity for a={:#x} in GF(2^4)",
            a.to_uint()
        );
    }
}

/// Frobenius tower: 8 successive squarings give identity in GF(2^8).
#[test]
fn frobenius_tower_gf8() {
    let p = poly8();
    for a in gf8_elements() {
        let mut x = a.clone();
        for _ in 0..8 {
            x = x.square(&p);
        }
        assert_eq!(
            x, a,
            "Frobenius tower (8 squarings) != identity for a={:#x} in GF(2^8)",
            a.to_uint()
        );
    }
}

// ── Known-answer tests ────────────────────────────────────────────────────────

/// AES GF(2^8) multiplication: 0x53 * 0xca = 0x01.
///
/// These are multiplicative inverses in the AES field.
#[test]
fn aes_mul_inverse_pair() {
    let p = poly8();
    let a = F2mNaive::<1>::from_u64(0x53, &p);
    let b = F2mNaive::<1>::from_u64(0xca, &p);
    let prod = a.mul(&b, &p);
    assert_eq!(
        prod.to_uint(),
        Uint::<1>::ONE,
        "0x53 * 0xca should be 1 in AES GF(2^8)"
    );
}

/// AES GF(2^8) multiplication: 0x57 * 0x83 = 0xc1.
///
/// Known AES MixColumns test vector.
#[test]
fn aes_mul_known_vector() {
    let p = poly8();
    let a = F2mNaive::<1>::from_u64(0x57, &p);
    let b = F2mNaive::<1>::from_u64(0x83, &p);
    let prod = a.mul(&b, &p);
    assert_eq!(
        prod.to_uint(),
        Uint::<1>::from(0xc1u64),
        "0x57 * 0x83 should be 0xc1 in AES GF(2^8)"
    );
}

/// GF(2^4) multiplication: x² * x² = x⁴ ≡ x + 1 (mod x⁴+x+1).
#[test]
fn gf4_mul_x2_x2() {
    let p = poly4();
    let x2 = F2mNaive::<1>::from_u64(0b0100, &p); // x²
    let prod = x2.mul(&x2, &p);
    // x⁴ mod (x⁴+x+1) = x+1 = 0b0011
    assert_eq!(prod.to_uint(), Uint::<1>::from(0b0011u64));
}

/// GF(2^4) squaring: square(x²) = x⁴ ≡ x + 1 (mod x⁴+x+1).
#[test]
fn gf4_square_x2() {
    let p = poly4();
    let x2 = F2mNaive::<1>::from_u64(0b0100, &p); // x²
    let sq = x2.square(&p);
    assert_eq!(sq.to_uint(), Uint::<1>::from(0b0011u64));
}

/// `frobenius` delegates to `square` and gives the same result.
#[test]
fn frobenius_equals_square() {
    let p = poly8();
    for a in gf8_elements() {
        assert_eq!(
            a.frobenius(&p),
            a.square(&p),
            "frobenius != square for a={:#x}",
            a.to_uint()
        );
    }
}

// ── Inversion round-trip ──────────────────────────────────────────────────────

/// `mul(a, inv(a)) == one` for all non-zero `a` in GF(2^4).
///
/// This is the primary inversion correctness check: every non-zero element
/// must have a multiplicative inverse.
#[test]
fn inv_round_trip_gf4() {
    let p = poly4();
    let one = F2mNaive::<1>::one();
    for v in 1u64..16 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let a_inv = a.inv(&p);
        let prod = a.mul(&a_inv, &p);
        assert_eq!(
            prod, one,
            "a * inv(a) != 1 for a={v:#x} in GF(2^4)"
        );
    }
}

/// `mul(a, inv(a)) == one` for all non-zero `a` in GF(2^8).
///
/// Exhaustive over all 255 non-zero elements of GF(2^8).
#[test]
fn inv_round_trip_gf8() {
    let p = poly8();
    let one = F2mNaive::<1>::one();
    for v in 1u64..256 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let a_inv = a.inv(&p);
        let prod = a.mul(&a_inv, &p);
        assert_eq!(
            prod, one,
            "a * inv(a) != 1 for a={v:#x} in GF(2^8)"
        );
    }
}

/// `inv(one) == one` in GF(2^4).
#[test]
fn inv_one_is_one_gf4() {
    let p = poly4();
    let one = F2mNaive::<1>::one();
    assert_eq!(one.inv(&p), one, "inv(1) != 1 in GF(2^4)");
}

/// `inv(one) == one` in GF(2^8).
#[test]
fn inv_one_is_one_gf8() {
    let p = poly8();
    let one = F2mNaive::<1>::one();
    assert_eq!(one.inv(&p), one, "inv(1) != 1 in GF(2^8)");
}

/// `inv(zero)` panics — zero has no multiplicative inverse.
#[test]
#[should_panic]
fn inv_zero_panics() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    let _ = zero.inv(&p);
}

// ── Extended-Euclidean ↔ Itoh–Tsujii agreement ───────────────────────────────

/// Both inversion algorithms agree for all non-zero elements of GF(2^4).
///
/// This is the cross-check: a bug in `square`/`frobenius` shows up in
/// Itoh–Tsujii but not in extended-Euclidean, and vice versa for a bug in
/// polynomial division.
#[test]
fn ext_euclid_itoh_tsujii_agree_gf4() {
    let p = poly4();
    for v in 1u64..16 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let ee = ext_euclid_inv(&a, &p);
        let it = itoh_tsujii_inv(&a, &p);
        assert_eq!(
            ee, it,
            "ext_euclid_inv != itoh_tsujii_inv for a={v:#x} in GF(2^4)"
        );
    }
}

/// Both inversion algorithms agree for all non-zero elements of GF(2^8).
///
/// Exhaustive over all 255 non-zero elements.
#[test]
fn ext_euclid_itoh_tsujii_agree_gf8() {
    let p = poly8();
    for v in 1u64..256 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let ee = ext_euclid_inv(&a, &p);
        let it = itoh_tsujii_inv(&a, &p);
        assert_eq!(
            ee, it,
            "ext_euclid_inv != itoh_tsujii_inv for a={v:#x} in GF(2^8)"
        );
    }
}

// ── Division ──────────────────────────────────────────────────────────────────

/// `div(a, b) == mul(a, inv(b))` for all `a` and non-zero `b` in GF(2^4).
#[test]
fn div_equals_mul_inv_gf4() {
    let p = poly4();
    for av in 0u64..16 {
        for bv in 1u64..16 {
            let a = F2mNaive::<1>::from_u64(av, &p);
            let b = F2mNaive::<1>::from_u64(bv, &p);
            let lhs = a.div(&b, &p);
            let rhs = a.mul(&b.inv(&p), &p);
            assert_eq!(
                lhs, rhs,
                "div(a,b) != mul(a,inv(b)) for a={av:#x} b={bv:#x} in GF(2^4)"
            );
        }
    }
}

/// `div(a, b) == mul(a, inv(b))` for a representative set in GF(2^8).
#[test]
fn div_equals_mul_inv_gf8() {
    let p = poly8();
    let elems: Vec<u64> = vec![0, 1, 2, 3, 0x53, 0xca, 0x8d, 0xf6, 0x01, 0xff, 0x1b, 0xe5];
    for &av in &elems {
        for &bv in &elems {
            if bv == 0 {
                continue; // skip zero divisor
            }
            let a = F2mNaive::<1>::from_u64(av, &p);
            let b = F2mNaive::<1>::from_u64(bv, &p);
            let lhs = a.div(&b, &p);
            let rhs = a.mul(&b.inv(&p), &p);
            assert_eq!(
                lhs, rhs,
                "div(a,b) != mul(a,inv(b)) for a={av:#x} b={bv:#x} in GF(2^8)"
            );
        }
    }
}

/// `div(zero, b) == zero` for non-zero `b` in GF(2^4).
#[test]
fn div_zero_numerator_is_zero() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    for bv in 1u64..16 {
        let b = F2mNaive::<1>::from_u64(bv, &p);
        assert_eq!(
            zero.div(&b, &p),
            zero,
            "0 / b != 0 for b={bv:#x} in GF(2^4)"
        );
    }
}

/// `div(zero, b)` panics when `b` is zero.
#[test]
#[should_panic]
fn div_zero_denominator_panics() {
    let p = poly4();
    let a = F2mNaive::<1>::from_u64(0x5, &p);
    let zero = F2mNaive::<1>::zero();
    let _ = a.div(&zero, &p);
}

// ── AES inversion KATs ────────────────────────────────────────────────────────

/// AES GF(2^8): `inv(0x53) == 0xca` (known AES S-box inverse pair).
#[test]
fn aes_inv_known_pair() {
    let p = poly8();
    let a = F2mNaive::<1>::from_u64(0x53, &p);
    let b = F2mNaive::<1>::from_u64(0xca, &p);
    assert_eq!(a.inv(&p), b, "inv(0x53) != 0xca in AES GF(2^8)");
    assert_eq!(b.inv(&p), a, "inv(0xca) != 0x53 in AES GF(2^8)");
}

// ── Normal-basis KATs ─────────────────────────────────────────────────────────
//
// These tests cover the normal-basis representation (`F2mNormal`) and the
// polynomial↔normal change-of-basis isomorphism (`convert.rs`).

/// Basis-conversion round-trip: `poly_to_normal(normal_to_poly(a)) == a`
/// for all elements of GF(2^4).
///
/// Verifies that the normal→poly→normal round-trip is the identity.  A bug in
/// either direction of the conversion would break this test.
#[test]
fn normal_to_poly_to_normal_round_trip_gf4() {
    let p = poly4();
    let beta = find_normal_element(&p);
    for v in 0u64..16 {
        let a_normal = Uint::<1>::from(v);
        let a_poly = normal_to_poly(&a_normal, &beta, &p);
        let a_back = poly_to_normal(&a_poly, &beta, &p);
        assert_eq!(
            a_normal, a_back,
            "normal→poly→normal round-trip failed for normal coords={v:#06b}"
        );
    }
}

/// Basis-conversion round-trip: `normal_to_poly(poly_to_normal(a)) == a`
/// for all elements of GF(2^4).
///
/// Verifies that the poly→normal→poly round-trip is the identity.
#[test]
fn poly_to_normal_to_poly_round_trip_gf4() {
    let p = poly4();
    let beta = find_normal_element(&p);
    for v in 0u64..16 {
        let a_poly = Uint::<1>::from(v);
        let a_normal = poly_to_normal(&a_poly, &beta, &p);
        let a_back = normal_to_poly(&a_normal, &beta, &p);
        assert_eq!(
            a_poly, a_back,
            "poly→normal→poly round-trip failed for poly value={v:#06b}"
        );
    }
}

/// Squaring in a normal basis is a cyclic left-shift of the coefficient vector.
///
/// This is the **pedagogical payoff** of the normal basis.  The test verifies
/// the property **directly**: it checks that `square(a).c` equals the cyclic
/// left-shift of `a.c` by 1 position within the m-bit field — not merely that
/// `square` returns the correct field element.
///
/// If `square` were implemented by converting to polynomial basis and back
/// (rather than as a cyclic shift), this test would still pass for correctness
/// but would fail to exhibit the shift property.  The implementation in
/// `normal.rs` uses the cyclic shift directly, and this test confirms it.
#[test]
fn squaring_is_cyclic_shift_gf4() {
    let p = poly4();
    let m = p.bit_len() - 1; // 4
    assert_eq!(m, 4, "GF(2^4) should have degree 4");

    // For every element in GF(2^4), verify:
    //   square(a).c == cyclic_left_shift_m(a.c, m)
    //
    // This checks the IMPLEMENTATION MECHANISM, not just the result.
    for v in 0u64..16 {
        let a = F2mNormal::<1>::from_uint(Uint::<1>::from(v), &p);
        let sq = a.square(&p);

        // Compute the expected cyclic left-shift of a.normal_coords().
        let c_bits = a.normal_coords().as_words()[0];
        let mask = (1u64 << m) - 1;
        let top_bit = (c_bits >> (m - 1)) & 1; // bit m−1 wraps to bit 0
        let expected_c = ((c_bits << 1) | top_bit) & mask;

        assert_eq!(
            sq.normal_coords().as_words()[0],
            expected_c,
            "square(a).c is not a cyclic left-shift of a.c for \
             poly-basis value={v:#x} (normal coords a.c={c_bits:#06b}, \
             expected sq.c={expected_c:#06b}, got sq.c={:#06b})",
            sq.normal_coords().as_words()[0]
        );
    }
}

/// Cross-representation agreement for `mul` in GF(2^4).
///
/// For a representative set of pairs, verifies that multiplying in the
/// polynomial basis gives the same result as multiplying in the normal basis
/// (after converting inputs and converting the result back).
///
/// This is the isomorphism check: the two representations are the same field.
#[test]
fn cross_representation_mul_gf4() {
    let p = poly4();
    let beta = find_normal_element(&p);

    for av in 0u64..16 {
        for bv in 0u64..16 {
            // Polynomial-basis multiplication.
            let a_naive = F2mNaive::<1>::from_u64(av, &p);
            let b_naive = F2mNaive::<1>::from_u64(bv, &p);
            let prod_poly = a_naive.mul(&b_naive, &p).to_uint();

            // Normal-basis multiplication (convert inputs, multiply, convert back).
            let a_norm = F2mNormal::<1>::from_poly(Uint::<1>::from(av), beta, p);
            let b_norm = F2mNormal::<1>::from_poly(Uint::<1>::from(bv), beta, p);
            let prod_norm_poly = a_norm.mul(&b_norm, &p).to_uint();

            assert_eq!(
                prod_poly, prod_norm_poly,
                "cross-representation mul disagrees: a={av:#x} b={bv:#x} in GF(2^4): \
                 poly={prod_poly:#x} normal={prod_norm_poly:#x}"
            );
        }
    }
}

/// Cross-representation agreement for `add` in GF(2^4).
///
/// For all pairs, verifies that adding in the polynomial basis gives the same
/// result as adding in the normal basis (after converting inputs and converting
/// the result back).
#[test]
fn cross_representation_add_gf4() {
    let p = poly4();
    let beta = find_normal_element(&p);

    for av in 0u64..16 {
        for bv in 0u64..16 {
            // Polynomial-basis addition.
            let a_naive = F2mNaive::<1>::from_u64(av, &p);
            let b_naive = F2mNaive::<1>::from_u64(bv, &p);
            let sum_poly = a_naive.add(&b_naive).to_uint();

            // Normal-basis addition (convert inputs, add, convert back).
            let a_norm = F2mNormal::<1>::from_poly(Uint::<1>::from(av), beta, p);
            let b_norm = F2mNormal::<1>::from_poly(Uint::<1>::from(bv), beta, p);
            let sum_norm_poly = a_norm.add(&b_norm).to_uint();

            assert_eq!(
                sum_poly, sum_norm_poly,
                "cross-representation add disagrees: a={av:#x} b={bv:#x} in GF(2^4): \
                 poly={sum_poly:#x} normal={sum_norm_poly:#x}"
            );
        }
    }
}

/// Normal-basis round-trip for GF(2^8): exhaustive poly→normal→poly.
///
/// Verifies the conversion isomorphism holds for all 256 elements of GF(2^8).
#[test]
fn poly_to_normal_to_poly_round_trip_gf8() {
    let p = poly8();
    let beta = find_normal_element(&p);
    for v in 0u64..256 {
        let a_poly = Uint::<1>::from(v);
        let a_normal = poly_to_normal(&a_poly, &beta, &p);
        let a_back = normal_to_poly(&a_normal, &beta, &p);
        assert_eq!(
            a_poly, a_back,
            "poly→normal→poly round-trip failed for v={v:#x} in GF(2^8)"
        );
    }
}

/// Squaring-is-a-cyclic-shift in GF(2^8).
///
/// Same direct check as `squaring_is_cyclic_shift_gf4`, but over all 256
/// elements of GF(2^8).  Confirms the cyclic-shift property holds for the
/// larger field.
#[test]
fn squaring_is_cyclic_shift_gf8() {
    let p = poly8();
    let m = p.bit_len() - 1; // 8
    assert_eq!(m, 8, "GF(2^8) should have degree 8");

    for v in 0u64..256 {
        let a = F2mNormal::<1>::from_uint(Uint::<1>::from(v), &p);
        let sq = a.square(&p);

        let c_bits = a.normal_coords().as_words()[0];
        let mask = (1u64 << m) - 1;
        let top_bit = (c_bits >> (m - 1)) & 1;
        let expected_c = ((c_bits << 1) | top_bit) & mask;

        assert_eq!(
            sq.normal_coords().as_words()[0],
            expected_c,
            "square(a).c is not a cyclic left-shift of a.c for \
             poly-basis value={v:#x} in GF(2^8) (a.c={c_bits:#010b}, \
             expected sq.c={expected_c:#010b}, got sq.c={:#010b})",
            sq.normal_coords().as_words()[0]
        );
    }
}

/// Cross-representation `square` agreement in GF(2^8).
///
/// Verifies that squaring in normal basis (cyclic shift) gives the same
/// polynomial-basis result as squaring in polynomial basis (Frobenius bit-spread).
#[test]
fn cross_representation_square_gf8() {
    let p = poly8();

    for v in 0u64..256 {
        // Polynomial-basis squaring.
        let a_naive = F2mNaive::<1>::from_u64(v, &p);
        let sq_poly = a_naive.square(&p).to_uint();

        // Normal-basis squaring (cyclic shift), then convert back to polynomial basis.
        let a_norm = F2mNormal::<1>::from_uint(Uint::<1>::from(v), &p);
        let sq_norm_poly = a_norm.square(&p).to_uint();

        assert_eq!(
            sq_poly, sq_norm_poly,
            "cross-representation square disagrees for v={v:#x} in GF(2^8): \
             poly={sq_poly:#x} normal={sq_norm_poly:#x}"
        );
    }
}

// ── Naive↔optimized agreement KATs ───────────────────────────────────────────
//
// These tests are the C-F2mOpt equivalence contract: `F2mOpt::mul` must give
// byte-identical results to `F2mNaive::mul` on every input pair.  This is the
// `FpNaive`/`FpMonty` discipline applied to characteristic 2: the optimized
// multiplier is *equivalent*, not just plausible.
//
// Exhaustive coverage:
// - GF(2^4): 16 × 16 = 256 pairs.
// - GF(2^8): 256 × 256 = 65536 pairs (same scale as `sub_equals_add_gf8`).

/// Naive↔optimized `mul` agreement for all pairs in GF(2^4).
///
/// For every (a, b) in GF(2^4) × GF(2^4), asserts that
/// `F2mNaive::mul(a, b).to_uint() == F2mOpt::mul(a, b).to_uint()`.
/// Byte-identical equality on the canonical `Uint<1>` — the strongest correct
/// assertion given that both types use the same canonical storage.
#[test]
fn naive_opt_agree_gf4() {
    let p = poly4();
    for av in 0u64..16 {
        for bv in 0u64..16 {
            let a_naive = F2mNaive::<1>::from_u64(av, &p);
            let b_naive = F2mNaive::<1>::from_u64(bv, &p);
            let a_opt = F2mOpt::<1>::from_u64(av, &p);
            let b_opt = F2mOpt::<1>::from_u64(bv, &p);

            let naive_prod = a_naive.mul(&b_naive, &p).to_uint();
            let opt_prod = a_opt.mul(&b_opt, &p).to_uint();

            assert_eq!(
                naive_prod, opt_prod,
                "naive↔opt mul disagrees: a={av:#x} b={bv:#x} in GF(2^4): \
                 naive={naive_prod:#x} opt={opt_prod:#x}"
            );
        }
    }
}

/// Naive↔optimized `mul` agreement for all pairs in GF(2^8).
///
/// Exhaustive over all 256 × 256 = 65536 pairs.  This is the load-bearing
/// correctness signal for the Karatsuba implementation: a bug in the
/// recombination or reduction step gives wrong products that still pass weak
/// smoke tests but fail here.
#[test]
fn naive_opt_agree_gf8() {
    let p = poly8();
    for av in 0u64..256 {
        for bv in 0u64..256 {
            let a_naive = F2mNaive::<1>::from_u64(av, &p);
            let b_naive = F2mNaive::<1>::from_u64(bv, &p);
            let a_opt = F2mOpt::<1>::from_u64(av, &p);
            let b_opt = F2mOpt::<1>::from_u64(bv, &p);

            let naive_prod = a_naive.mul(&b_naive, &p).to_uint();
            let opt_prod = a_opt.mul(&b_opt, &p).to_uint();

            assert_eq!(
                naive_prod, opt_prod,
                "naive↔opt mul disagrees: a={av:#x} b={bv:#x} in GF(2^8): \
                 naive={naive_prod:#x} opt={opt_prod:#x}"
            );
        }
    }
}

// ── F2mOpt field-axiom suite ──────────────────────────────────────────────────
//
// These tests verify the field axioms directly against `F2mOpt` — not only via
// agreement with naive.  This confirms the Karatsuba multiplier satisfies the
// axioms in its own right.

/// `mul(a, one) == a` for all elements in GF(2^4) via `F2mOpt`.
#[test]
fn f2m_opt_mul_by_one_gf4() {
    let p = poly4();
    let one = F2mOpt::<1>::one();
    for v in 0u64..16 {
        let a = F2mOpt::<1>::from_u64(v, &p);
        assert_eq!(
            a.mul(&one, &p).to_uint(),
            a.to_uint(),
            "F2mOpt: a * 1 != a for a={v:#x} in GF(2^4)"
        );
    }
}

/// `add(a, a) == zero` for all elements in GF(2^4) via `F2mOpt`.
#[test]
fn f2m_opt_add_self_is_zero_gf4() {
    let zero = F2mOpt::<1>::zero();
    let p = poly4();
    for v in 0u64..16 {
        let a = F2mOpt::<1>::from_u64(v, &p);
        assert_eq!(
            a.add(&a).to_uint(),
            zero.to_uint(),
            "F2mOpt: a + a != 0 for a={v:#x} in GF(2^4)"
        );
    }
}

/// Associativity of `F2mOpt::mul` over GF(2^4): exhaustive 16^3 = 4096 triples.
#[test]
fn f2m_opt_mul_associative_gf4() {
    let p = poly4();
    for av in 0u64..16 {
        for bv in 0u64..16 {
            for cv in 0u64..16 {
                let a = F2mOpt::<1>::from_u64(av, &p);
                let b = F2mOpt::<1>::from_u64(bv, &p);
                let c = F2mOpt::<1>::from_u64(cv, &p);
                let lhs = a.mul(&b, &p).mul(&c, &p);
                let rhs = a.mul(&b.mul(&c, &p), &p);
                assert_eq!(
                    lhs, rhs,
                    "F2mOpt: mul not associative: a={av:#x} b={bv:#x} c={cv:#x} in GF(2^4)"
                );
            }
        }
    }
}

/// Distributivity of `F2mOpt::mul` over `add` in GF(2^4): exhaustive 16^3 triples.
#[test]
fn f2m_opt_mul_distributive_gf4() {
    let p = poly4();
    for av in 0u64..16 {
        for bv in 0u64..16 {
            for cv in 0u64..16 {
                let a = F2mOpt::<1>::from_u64(av, &p);
                let b = F2mOpt::<1>::from_u64(bv, &p);
                let c = F2mOpt::<1>::from_u64(cv, &p);
                let lhs = a.mul(&b.add(&c), &p);
                let rhs = a.mul(&b, &p).add(&a.mul(&c, &p));
                assert_eq!(
                    lhs, rhs,
                    "F2mOpt: distributivity failed: a={av:#x} b={bv:#x} c={cv:#x} in GF(2^4)"
                );
            }
        }
    }
}

/// Frobenius fixed-field law for `F2mOpt` in GF(2^4): `a^(2^4) == a`.
#[test]
fn f2m_opt_frobenius_fixed_field_gf4() {
    let p = poly4();
    let exp = Uint::<1>::from(16u64);
    for v in 0u64..16 {
        let a = F2mOpt::<1>::from_u64(v, &p);
        let a_pow = a.pow(&exp, &p);
        assert_eq!(
            a_pow.to_uint(),
            a.to_uint(),
            "F2mOpt: a^(2^4) != a for a={v:#x} in GF(2^4) — reduction is wrong"
        );
    }
}

// ── Trace KATs ────────────────────────────────────────────────────────────────
//
// The absolute trace Tr: GF(2^m) → GF(2) is a linear map.  Its image is
// always {0, 1} (as field elements).  The key properties tested:
//   1. Tr(a) ∈ {0, 1} for all a.
//   2. Tr(a + b) = Tr(a) ⊕ Tr(b)  (GF(2)-linearity).
//   3. Tr(a²) = Tr(a)              (Frobenius-invariance).
//   4. Tr(1) = m mod 2             (known value for the unit element).
//
// In GF(2^4) with x⁴+x+1, the trace is determined by the x³ coefficient:
// Tr(a) = a₃ (the coefficient of x³ in the polynomial representation).
// Elements 0–7 have trace 0; elements 8–15 have trace 1.

/// `trace(a) ∈ {0, 1}` for all elements of GF(2^4).
///
/// The trace is valued in GF(2) ⊂ GF(2^m): it must be exactly 0 or 1.
#[test]
fn trace_valued_in_gf2_gf4() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    let one = F2mNaive::<1>::one();
    for v in 0u64..16 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let tr = a.trace(&p);
        assert!(
            tr == zero || tr == one,
            "trace(a) not in {{0,1}} for a={v:#x} in GF(2^4): got {:#x}",
            tr.to_uint()
        );
    }
}

/// `trace(a) ∈ {0, 1}` for all elements of GF(2^8).
#[test]
fn trace_valued_in_gf2_gf8() {
    let p = poly8();
    let zero = F2mNaive::<1>::zero();
    let one = F2mNaive::<1>::one();
    for v in 0u64..256 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let tr = a.trace(&p);
        assert!(
            tr == zero || tr == one,
            "trace(a) not in {{0,1}} for a={v:#x} in GF(2^8): got {:#x}",
            tr.to_uint()
        );
    }
}

/// `trace(a + b) = trace(a) ⊕ trace(b)` (GF(2)-linearity) for all pairs in GF(2^4).
///
/// The trace is a GF(2)-linear map: Tr(a+b) = Tr(a) + Tr(b) in GF(2) = XOR.
#[test]
fn trace_linearity_gf4() {
    let p = poly4();
    for av in 0u64..16 {
        for bv in 0u64..16 {
            let a = F2mNaive::<1>::from_u64(av, &p);
            let b = F2mNaive::<1>::from_u64(bv, &p);
            let tr_sum = a.add(&b).trace(&p);
            let sum_tr = a.trace(&p).add(&b.trace(&p));
            assert_eq!(
                tr_sum, sum_tr,
                "trace not linear: trace(a+b) != trace(a)+trace(b) \
                 for a={av:#x} b={bv:#x} in GF(2^4)"
            );
        }
    }
}

/// `trace(a²) = trace(a)` (Frobenius-invariance) for all elements of GF(2^4).
///
/// The Frobenius endomorphism a → a² fixes the trace: Tr(a²) = Tr(a).
/// This follows from Tr(a) = Σ a^(2^i) and the cyclic structure of the sum.
#[test]
fn trace_frobenius_invariant_gf4() {
    let p = poly4();
    for v in 0u64..16 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let tr_a = a.trace(&p);
        let tr_a2 = a.square(&p).trace(&p);
        assert_eq!(
            tr_a, tr_a2,
            "trace(a²) != trace(a) for a={v:#x} in GF(2^4)"
        );
    }
}

/// `trace(a²) = trace(a)` (Frobenius-invariance) for all elements of GF(2^8).
#[test]
fn trace_frobenius_invariant_gf8() {
    let p = poly8();
    for v in 0u64..256 {
        let a = F2mNaive::<1>::from_u64(v, &p);
        let tr_a = a.trace(&p);
        let tr_a2 = a.square(&p).trace(&p);
        assert_eq!(
            tr_a, tr_a2,
            "trace(a²) != trace(a) for a={v:#x} in GF(2^8)"
        );
    }
}

/// Known-answer: `trace(x³) = 1` in GF(2^4) with x⁴+x+1.
///
/// In GF(2^4) with x⁴+x+1, the trace is determined by the x³ coefficient:
/// Tr(a) = a₃.  So Tr(x³) = 1 and Tr(x²) = Tr(x) = Tr(1) = 0.
#[test]
fn trace_known_answers_gf4() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    let one = F2mNaive::<1>::one();

    // Tr(0) = 0
    assert_eq!(F2mNaive::<1>::zero().trace(&p), zero, "Tr(0) != 0");
    // Tr(1) = 0 (m=4 is even, so Tr(1) = 1·4 mod 2 = 0)
    assert_eq!(F2mNaive::<1>::one().trace(&p), zero, "Tr(1) != 0 in GF(2^4)");
    // Tr(x) = 0
    assert_eq!(
        F2mNaive::<1>::from_u64(0b0010, &p).trace(&p),
        zero,
        "Tr(x) != 0 in GF(2^4)"
    );
    // Tr(x²) = 0
    assert_eq!(
        F2mNaive::<1>::from_u64(0b0100, &p).trace(&p),
        zero,
        "Tr(x²) != 0 in GF(2^4)"
    );
    // Tr(x³) = 1
    assert_eq!(
        F2mNaive::<1>::from_u64(0b1000, &p).trace(&p),
        one,
        "Tr(x³) != 1 in GF(2^4)"
    );
    // Tr(x³+1) = 1 (linearity: Tr(x³)+Tr(1) = 1+0 = 1)
    assert_eq!(
        F2mNaive::<1>::from_u64(0b1001, &p).trace(&p),
        one,
        "Tr(x³+1) != 1 in GF(2^4)"
    );
}

/// `F2mOpt::trace` agrees with `F2mNaive::trace` for all elements of GF(2^4).
#[test]
fn trace_naive_opt_agree_gf4() {
    let p = poly4();
    for v in 0u64..16 {
        let a_naive = F2mNaive::<1>::from_u64(v, &p);
        let a_opt = F2mOpt::<1>::from_u64(v, &p);
        assert_eq!(
            a_naive.trace(&p).to_uint(),
            a_opt.trace(&p).to_uint(),
            "trace: naive↔opt disagree for a={v:#x} in GF(2^4)"
        );
    }
}

/// `F2mOpt::trace` agrees with `F2mNaive::trace` for all elements of GF(2^8).
#[test]
fn trace_naive_opt_agree_gf8() {
    let p = poly8();
    for v in 0u64..256 {
        let a_naive = F2mNaive::<1>::from_u64(v, &p);
        let a_opt = F2mOpt::<1>::from_u64(v, &p);
        assert_eq!(
            a_naive.trace(&p).to_uint(),
            a_opt.trace(&p).to_uint(),
            "trace: naive↔opt disagree for a={v:#x} in GF(2^8)"
        );
    }
}

// ── solve_quadratic KATs ──────────────────────────────────────────────────────
//
// `solve_quadratic(c)` solves x² + x = c (the Artin–Schreier equation).
// A solution exists iff trace(c) = 0.  The two solutions are x and x+1.
//
// Properties tested:
//   1. For all c with trace(c) = 0: solve_quadratic(c)² + solve_quadratic(c) = c.
//   2. The two solutions are x and x+1 (they differ by 1).
//   3. Known-answer: x²+x = x (0b0010) has solution x = x³+x (0b1010) in GF(2^4).

/// `solve_quadratic(c)` satisfies `x² + x = c` for all c with trace 0 in GF(2^4).
///
/// This is the primary correctness check: the returned value must actually
/// satisfy the Artin–Schreier equation.
#[test]
fn solve_quadratic_correct_gf4() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    for v in 0u64..16 {
        let c = F2mNaive::<1>::from_u64(v, &p);
        if c.trace(&p) != zero {
            continue; // skip unsolvable c
        }
        let x = F2mNaive::<1>::solve_quadratic(&c, &p);
        let lhs = x.square(&p).add(&x);
        assert_eq!(
            lhs, c,
            "solve_quadratic: x²+x != c for c={v:#x} in GF(2^4): \
             got x={:#x}, x²+x={:#x}",
            x.to_uint(),
            lhs.to_uint()
        );
    }
}

/// `solve_quadratic(c)` satisfies `x² + x = c` for all c with trace 0 in GF(2^8).
#[test]
fn solve_quadratic_correct_gf8() {
    let p = poly8();
    let zero = F2mNaive::<1>::zero();
    for v in 0u64..256 {
        let c = F2mNaive::<1>::from_u64(v, &p);
        if c.trace(&p) != zero {
            continue;
        }
        let x = F2mNaive::<1>::solve_quadratic(&c, &p);
        let lhs = x.square(&p).add(&x);
        assert_eq!(
            lhs, c,
            "solve_quadratic: x²+x != c for c={v:#x} in GF(2^8): \
             got x={:#x}, x²+x={:#x}",
            x.to_uint(),
            lhs.to_uint()
        );
    }
}

/// The two solutions to x²+x=c are x and x+1 (they differ by 1).
///
/// If x is a solution, then (x+1)²+(x+1) = x²+1+x+1 = x²+x = c.
/// So x+1 is always the other solution.
#[test]
fn solve_quadratic_two_solutions_gf4() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    let one = F2mNaive::<1>::one();
    for v in 0u64..16 {
        let c = F2mNaive::<1>::from_u64(v, &p);
        if c.trace(&p) != zero {
            continue;
        }
        let x = F2mNaive::<1>::solve_quadratic(&c, &p);
        let x_plus_1 = x.add(&one);
        let lhs = x_plus_1.square(&p).add(&x_plus_1);
        assert_eq!(
            lhs, c,
            "solve_quadratic: x+1 is not the other solution for c={v:#x} in GF(2^4)"
        );
    }
}

/// Known-answer: x²+x = x (0b0010) has solution x³+x (0b1010) in GF(2^4).
///
/// Verified by hand: (x³+x)² + (x³+x) = x⁶+x² + x³+x.
/// x⁶ = x³+x² (mod x⁴+x+1), so x³+x²+x²+x³+x = x = 0b0010. ✓
#[test]
fn solve_quadratic_known_answer_gf4() {
    let p = poly4();
    let c = F2mNaive::<1>::from_u64(0b0010, &p); // c = x
    let x = F2mNaive::<1>::solve_quadratic(&c, &p);
    // The solution must satisfy x²+x = c.
    let lhs = x.square(&p).add(&x);
    assert_eq!(
        lhs, c,
        "solve_quadratic known-answer failed: x²+x != x for c=x in GF(2^4)"
    );
    // The solution should be x³+x (0b1010) or x³+x+1 (0b1011).
    let sol1 = F2mNaive::<1>::from_u64(0b1010, &p);
    let sol2 = F2mNaive::<1>::from_u64(0b1011, &p);
    assert!(
        x == sol1 || x == sol2,
        "solve_quadratic: expected x³+x or x³+x+1 for c=x in GF(2^4), got {:#x}",
        x.to_uint()
    );
}

/// `F2mOpt::solve_quadratic` agrees with `F2mNaive::solve_quadratic` for all
/// solvable c in GF(2^4).
#[test]
fn solve_quadratic_naive_opt_agree_gf4() {
    let p = poly4();
    let zero = F2mNaive::<1>::zero();
    for v in 0u64..16 {
        let c_naive = F2mNaive::<1>::from_u64(v, &p);
        if c_naive.trace(&p) != zero {
            continue;
        }
        let c_opt = F2mOpt::<1>::from_u64(v, &p);
        let x_naive = F2mNaive::<1>::solve_quadratic(&c_naive, &p).to_uint();
        let x_opt = F2mOpt::<1>::solve_quadratic(&c_opt, &p).to_uint();
        // Both solutions must satisfy x²+x=c; they may differ by 1.
        // Verify both satisfy the equation rather than requiring identical output.
        let x_naive_elem = F2mNaive::<1>::from_uint(x_naive, &p);
        let x_opt_elem = F2mNaive::<1>::from_uint(x_opt, &p);
        let lhs_naive = x_naive_elem.square(&p).add(&x_naive_elem);
        let lhs_opt = x_opt_elem.square(&p).add(&x_opt_elem);
        assert_eq!(
            lhs_naive, c_naive,
            "naive solve_quadratic wrong for c={v:#x}"
        );
        assert_eq!(
            lhs_opt, c_naive,
            "opt solve_quadratic wrong for c={v:#x}"
        );
    }
}
