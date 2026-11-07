//! Known-answer tests for split-prime local contributions to an algebraic square root.
//!
//! Every fixture states the arithmetic that justifies it, derived by hand and not by running the
//! code under test. Throughout, `θ = c·α` is the root of the monic `F(t) = c^(d−1) f(t/c)`,
//! `η = c^(2h−k) Π(c·a − b·θ)` with `h = ⌈k/2⌉`, `H = F'(θ)² η`, and a root of the target is
//! `z = c^h F'(θ) β` with `β² = Π(a − b·α)`. A prime `p` is usable when `F` splits into distinct
//! linear factors modulo `p` and `H(r) ≢ 0` at each root `r`; its roots lift to `R ≡ r` modulo
//! `Q = p^e`, and `s² ≡ H(R) (mod Q)`. The Lagrange polynomials `L` of the lifted roots and the
//! Chinese-remainder weights `W` (`≡ 1` modulo its own `Q`, `≡ 0` modulo the others) give the
//! contributions `C = W·s·L mod M`. Any root `z` is `Σ ε·C` for some signs `ε`, because
//! `z(R) ≡ ±s`. The tests derive those signs from a known `z` and never ask the code under test
//! to choose them.

use gnfs::sqrt::{
    CandidateVerdict, ContributionSet, ContributionSign, CrtLimits, NonSquareEvidence,
    NormalizedTarget, PrecisionLimits, PrimeSearchLimits, ReconstructionMethod, RelationFactors,
    ResourceExhaustion, ResourceKind, SplitOutcome, SqrtContext, SqrtError, TargetOutcome,
    split_contributions,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::pow::pow;
use shared_numfield::{IntPoly, NonresidueSearchLimits, RootSearchLimits};

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Integer polynomial from coefficients listed lowest degree first.
fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

fn target_big(f: &[i64], m: i64, n: i64, pairs: &[(BigInt, BigInt)]) -> NormalizedTarget {
    let ctx = SqrtContext::new(&poly(f), &big(m), &big(n)).expect("fixture is a valid context");
    let factors = RelationFactors::new(&ctx, pairs).expect("fixture relations are units");
    match NormalizedTarget::new(
        factors,
        PrecisionLimits {
            max_modulus_bits: 4096,
        },
    )
    .expect("a validated target builds")
    {
        TargetOutcome::Ready(target) => target,
        TargetOutcome::ResourceExhausted(e) => panic!("unexpected exhaustion {e:?}"),
    }
}

fn target(f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> NormalizedTarget {
    let pairs: Vec<(BigInt, BigInt)> = pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
    target_big(f, m, n, &pairs)
}

/// Primes from `first` on, at most `examined` examined, `max_primes` used, and the two inner
/// budgets.
fn limits(first: u64, examined: u64, max_primes: u32, evaluations: u64, nonres: u64) -> CrtLimits {
    CrtLimits {
        prime_search: PrimeSearchLimits {
            first_candidate: first,
            max_primes_examined: examined,
        },
        max_primes,
        root_search: RootSearchLimits {
            max_evaluations: evaluations,
        },
        nonresidue_search: NonresidueSearchLimits {
            max_candidates: nonres,
        },
    }
}

/// Starts at 3 with room for 200 primes and generous inner budgets.
fn generous(max_primes: u32) -> CrtLimits {
    limits(0, 200, max_primes, 100_000, 1000)
}

fn run(target: &NormalizedTarget, limits: CrtLimits) -> SplitOutcome {
    split_contributions(target, limits).expect("no invariant violation")
}

fn ready(outcome: SplitOutcome) -> ContributionSet {
    match outcome {
        SplitOutcome::Ready(set) => set,
        other => panic!("expected contributions, got {other:?}"),
    }
}

fn exhausted(outcome: SplitOutcome) -> ResourceExhaustion {
    match outcome {
        SplitOutcome::ResourceExhausted(e) => e,
        other => panic!("expected resource exhaustion, got {other:?}"),
    }
}

fn evidence(outcome: SplitOutcome) -> NonSquareEvidence {
    match outcome {
        SplitOutcome::NotSquare(certificate) => certificate.evidence().clone(),
        other => panic!("expected a nonsquare certificate, got {other:?}"),
    }
}

fn primes_of(set: &ContributionSet) -> Vec<u64> {
    set.locals().iter().map(|l| l.prime()).collect()
}

/// `a·b mod (F, m)` computed with polynomial division, independent of the code under test.
fn product_mod(a: &IntPoly, b: &IntPoly, f: &IntPoly, m: &BigInt) -> IntPoly {
    let remainder = a.mul(b).pseudo_div_rem(f).1;
    IntPoly::from_coeffs(remainder.coeffs.iter().map(|c| c.mod_floor(m)).collect())
}

/// The sign vector with `Minus` exactly where the bit of `mask` is set.
fn signs_from_mask(n: usize, mask: u32) -> Vec<ContributionSign> {
    (0..n)
        .map(|k| {
            if mask >> k & 1 == 1 {
                ContributionSign::Minus
            } else {
                ContributionSign::Plus
            }
        })
        .collect()
}

/// Signs `ε` with `z(R) ≡ ε·s (mod Q)` at every contribution, from the lifted roots alone.
fn signs_of_known_root(set: &ContributionSet, z: &IntPoly) -> Vec<ContributionSign> {
    set.contributions()
        .iter()
        .map(|c| {
            let q = set.locals()[c.prime_index()].modulus();
            let value = z.eval_mod(c.lifted_root(), q);
            if value == *c.local_root() {
                ContributionSign::Plus
            } else if value == q - c.local_root() {
                ContributionSign::Minus
            } else {
                panic!(
                    "z(R) = {value} is neither ±s for the root {} of {}",
                    c.root(),
                    c.prime()
                )
            }
        })
        .collect()
}

fn is_root_up_to_sign(candidate: &IntPoly, z: &IntPoly) -> bool {
    candidate == z || *candidate == z.neg()
}

/// How many of the `2^n` sign vectors assemble to `±z`.
fn count_assemblies_of(set: &ContributionSet, z: &IntPoly) -> usize {
    let n = set.contributions().len();
    (0..1u32 << n)
        .filter(|&mask| {
            let candidate = set
                .assemble(&signs_from_mask(n, mask))
                .expect("length matches");
            is_root_up_to_sign(&candidate, z)
        })
        .count()
}

/// Every local identity of every contribution, recomputed here from public data only.
fn assert_local_identities(set: &ContributionSet) {
    let context = set.target().context();
    let monic = context.monic_poly();
    let h = set.target().target();
    let d = context.degree();
    let modulus = set.modulus();
    let zero = big(0);
    let one = big(1);

    // M = Π Q_i ≥ 2B + 1, odd, with Q_i = p_i^e; e is the least common exponent.
    let product_of_primes: BigInt = set.locals().iter().map(|l| big(l.prime() as i64)).product();
    let e = set.exponent() as usize;
    assert_eq!(*modulus, pow(product_of_primes.clone(), e));
    assert!(modulus >= set.target().minimum_modulus());
    assert!(e == 1 || pow(product_of_primes, e - 1) < *set.target().minimum_modulus());
    assert!(modulus.is_odd());
    assert_eq!(set.contributions().len(), set.locals().len() * d);

    // Chinese-remainder weights: W_i ≡ δ_ik mod Q_k, W_i² ≡ W_i, Σ W_i ≡ 1 mod M.
    for (i, wi) in set.locals().iter().enumerate() {
        assert_eq!(*wi.modulus(), pow(big(wi.prime() as i64), e));
        assert!(*wi.weight() >= zero && wi.weight() < modulus);
        for (k, wk) in set.locals().iter().enumerate() {
            let expected = if i == k { &one } else { &zero };
            assert_eq!(
                wi.weight().mod_floor(wk.modulus()),
                *expected,
                "W_{i} mod Q_{k}"
            );
        }
        assert_eq!(
            (wi.weight() * wi.weight() - wi.weight()).mod_floor(modulus),
            zero
        );
    }
    let weight_sum: BigInt = set.locals().iter().map(|l| l.weight().clone()).sum();
    assert_eq!((weight_sum - &one).mod_floor(modulus), zero);

    for c in set.contributions() {
        let local = &set.locals()[c.prime_index()];
        let (p, q) = (big(c.prime() as i64), local.modulus());
        assert_eq!(c.prime(), local.prime());
        // R ∈ [0, Q), R ≡ r mod p, F(R) ≡ 0 mod Q.
        assert!(*c.lifted_root() >= zero && c.lifted_root() < q);
        assert_eq!(c.lifted_root().mod_floor(&p), big(c.root() as i64));
        assert_eq!(monic.eval_mod(c.lifted_root(), q), zero, "F(R) mod Q");
        // s ∈ [0, Q), s² ≡ H(R) mod Q, and the other sign is equally a square root.
        let value = h.eval_mod(c.lifted_root(), q);
        assert!(*c.local_root() >= zero && c.local_root() < q);
        assert_eq!((c.local_root() * c.local_root()).mod_floor(q), value);
        assert_ne!(value.mod_floor(&p), zero, "H(R) must be a unit");
        let negated = q - c.local_root();
        assert_eq!((&negated * &negated).mod_floor(q), value);
        // Lagrange polynomial: degree d−1, coefficients in [0, Q), 1 at R and 0 at the other
        // roots of the same prime.
        assert_eq!(c.lagrange().coeffs.len(), d);
        assert!(c.lagrange().coeffs.iter().all(|x| *x >= zero && x < q));
        for other in set
            .contributions()
            .iter()
            .filter(|o| o.prime_index() == c.prime_index())
        {
            let at = c.lagrange().eval_mod(other.lifted_root(), q);
            let expected = if std::ptr::eq(other, c) { &one } else { &zero };
            assert_eq!(
                at,
                *expected,
                "L(R) at root {} of {}",
                other.root(),
                c.prime()
            );
        }
        // C = W·s·L mod M, d entries in [0, M).
        assert_eq!(c.coefficients().len(), d);
        for (k, coefficient) in c.coefficients().iter().enumerate() {
            let l_k = c
                .lagrange()
                .coeffs
                .get(k)
                .cloned()
                .unwrap_or_else(|| big(0));
            let expected = (local.weight() * c.local_root() * l_k).mod_floor(modulus);
            assert_eq!(*coefficient, expected, "C[{k}]");
        }
    }
    // The Lagrange polynomials of one prime sum to 1 modulo Q.
    for (index, local) in set.locals().iter().enumerate() {
        let mut total = vec![big(0); d];
        for c in set
            .contributions()
            .iter()
            .filter(|c| c.prime_index() == index)
        {
            for (slot, x) in total.iter_mut().zip(&c.lagrange().coeffs) {
                *slot += x;
            }
        }
        for (i, x) in total.iter().enumerate() {
            let expected = if i == 0 { 1 } else { 0 };
            assert_eq!(
                x.mod_floor(local.modulus()),
                big(expected),
                "Σ L, coefficient {i}"
            );
        }
    }
}

/// Every sign vector assembles to a square root of `H` modulo `(F, M)`.
fn assert_every_sign_vector_is_a_local_square(set: &ContributionSet) {
    let f = set.target().context().monic_poly();
    let m = set.modulus();
    let h = IntPoly::from_coeffs(
        set.target()
            .target()
            .coeffs
            .iter()
            .map(|c| c.mod_floor(m))
            .collect(),
    );
    let n = set.contributions().len();
    for mask in 0..1u32 << n {
        let z = set
            .assemble(&signs_from_mask(n, mask))
            .expect("length matches");
        assert!(
            z.coeffs.len() <= set.target().context().degree(),
            "degree below d"
        );
        assert_eq!(
            product_mod(&z, &z, f, m),
            h,
            "z² ≢ H for sign mask {mask:b}"
        );
    }
}

fn verified_root_up_to_sign(set: &ContributionSet, z: &IntPoly, candidate: &IntPoly) {
    match set
        .target()
        .verify_candidate(candidate, set.method())
        .expect("no invariant violation")
    {
        CandidateVerdict::Verified(root) => {
            assert!(is_root_up_to_sign(root.integral_root(), z));
            assert_eq!(*root.method(), set.method());
        }
        other => panic!("expected a verified root, got {other:?}"),
    }
}

/// `x² − 2`, `m = 23`, `N = 527 = 17·31`, relation `(3, −2)`: `3 + 2α = (1 + α)²`, `z = 4 + 2θ`.
fn target_527() -> NormalizedTarget {
    target(&[-2, 0, 1], 23, 527, &[(3, -2)])
}

/// `F = t⁴ − 10t² + 1`, `m = 4`, `N = 97` (`f(4) = 256 − 160 + 1 = 97`), relation `(1, −1)` twice:
/// `β = 1 + θ`, `z = F'β = (4θ³ − 20θ)(1 + θ) = −4 − 20θ + 20θ² + 4θ³` after `θ⁴ = 10θ² − 1`.
fn target_v4() -> NormalizedTarget {
    target(&[1, 0, -10, 0, 1], 4, 97, &[(1, -1), (1, -1)])
}

fn z_v4() -> IntPoly {
    poly(&[-4, -20, 20, 4])
}

#[test]
fn kat_split_selection_skips_unusable_primes() {
    // F = t² − 2, c = 1, relation (3, −2): H = F'²η = 4t²(3 + 2t) = 8(3 + 2t) = 24 + 16t, and
    // B = (18, 6), so the least modulus above 2B is 37. Modulo p, t² = 2 has a root iff 2 is a
    // square: 3 and 5 (squares {1}, {1, 4}) and 11 ({1, 3, 4, 5, 9}), 13 ({1, 3, 4, 9, 10, 12}),
    // 19 ({1, 4, 5, 6, 7, 9, 11, 16, 17}) have none. 7 has 3² = 9 ≡ 2, 17 has 6² = 36 ≡ 2 (but
    // 17 | 527), and 23 has 5² = 25 ≡ 2. Primes examined: 3, 5, 7, 11, 13, 17, 19, 23 (8).
    let t = target_527();
    assert_eq!(t.minimum_modulus(), &big(37));
    let h = t.target();
    assert_eq!(*h, poly(&[24, 16]));

    // One prime wanted: stops at 7 after three primes. Roots 3, 4; H(3) = 24 + 48 = 72 ≡ 2 and
    // H(4) = 24 + 64 = 88 ≡ 4 modulo 7, both squares (3² = 2, 2² = 4).
    let set = ready(run(&t, generous(1)));
    assert_eq!(primes_of(&set), vec![7]);
    let roots: Vec<u64> = set.contributions().iter().map(|c| c.root()).collect();
    assert_eq!(roots, vec![3, 4]);
    assert_eq!(h.eval_mod(&big(3), &big(7)), big(2));
    assert_eq!(h.eval_mod(&big(4), &big(7)), big(4));

    // Two primes wanted: 17 splits but divides N, so the second is 23 after 8 examined primes.
    // Roots 5 and 18 (5² = 25 ≡ 2, 18 = −5); H(5) = 104 ≡ 12 (9² = 81 = 3·23 + 12) and
    // H(18) = 312 ≡ 13 (6² = 36 = 23 + 13), both squares.
    let set = ready(run(&t, limits(0, 8, 2, 100_000, 1000)));
    assert_eq!(primes_of(&set), vec![7, 23]);
    let layout: Vec<(usize, u64, u64)> = set
        .contributions()
        .iter()
        .map(|c| (c.prime_index(), c.prime(), c.root()))
        .collect();
    assert_eq!(layout, vec![(0, 7, 3), (0, 7, 4), (1, 23, 5), (1, 23, 18)]);
    assert_eq!(h.eval_mod(&big(5), &big(23)), big(12));
    assert_eq!(h.eval_mod(&big(18), &big(23)), big(13));

    // The prime budget counts every examined prime: a budget of 6 or 7 ends before 23, so the
    // partial prime set [7] is used and no exhaustion is reported; a budget of 8 reaches 23.
    for examined in [6, 7] {
        let set = ready(run(&t, limits(0, examined, 2, 100_000, 1000)));
        assert_eq!(primes_of(&set), vec![7], "budget {examined}");
    }

    // 17 | N is skipped although F splits modulo 17; starting at 17 the single prime is 23.
    assert_eq!(limits_ready_primes(&t, 17, 1), vec![23]);
    // Starting after 7 the same holds: 11, 13 and 19 are inert.
    assert_eq!(limits_ready_primes(&t, 8, 1), vec![23]);
    // The set owns the target it was built from, and records the method for the verifier.
    assert_eq!(set.target(), &t);
    assert_eq!(
        set.method(),
        ReconstructionMethod::SplitPrimeCrt {
            primes: vec![7, 23],
            exponent: 1
        }
    );
}

/// The primes selected from `first` on with a generous budget and `max_primes` wanted.
fn limits_ready_primes(target: &NormalizedTarget, first: u64, max_primes: u32) -> Vec<u64> {
    primes_of(&ready(run(
        target,
        limits(first, 200, max_primes, 100_000, 1000),
    )))
}

#[test]
fn kat_lifted_roots_satisfy_polynomial_mod_prime_power() {
    let t = target_527();
    // ℓ = 1: 7 < 37 ≤ 7² = 49, so e = 2 and Q = M = 49. F(t) = t² − 2 lifts the roots 3, 4 of
    // F mod 7: 3 + 7k with (3 + 7k)² = 9 + 42k ≡ 2 mod 49 ⇒ 7 + 42k ≡ 0 ⇒ 1 + 6k ≡ 0 mod 7 ⇒ k = 1,
    // so R = 10 (10² = 100 = 2·49 + 2); the other root is 49 − 10 = 39 ≡ 4 mod 7.
    let set = ready(run(&t, generous(1)));
    assert_eq!(set.exponent(), 2);
    assert_eq!(set.modulus(), &big(49));
    assert_eq!(set.locals()[0].modulus(), &big(49));
    let lifted: Vec<(u64, BigInt)> = set
        .contributions()
        .iter()
        .map(|c| (c.root(), c.lifted_root().clone()))
        .collect();
    assert_eq!(lifted, vec![(3, big(10)), (4, big(39))]);
    assert_local_identities(&set);

    // ℓ = 2: 7·23 = 161 ≥ 37, so e = 1 and the lifts are the residues themselves.
    let set = ready(run(&t, generous(2)));
    assert_eq!(set.exponent(), 1);
    assert_eq!(set.modulus(), &big(161));
    let lifted: Vec<BigInt> = set
        .contributions()
        .iter()
        .map(|c| c.lifted_root().clone())
        .collect();
    assert_eq!(lifted, vec![big(3), big(4), big(5), big(18)]);
    assert_local_identities(&set);

    // V4: B = 69168, minimum modulus 138337. 23 is the first completely split prime, with roots
    // ±2, ±11 (2⁴ − 40 + 1 = −23, 11² = 121 ≡ 6 and 36 − 60 + 1 = −23). 23³ = 12167 < 138337 ≤
    // 23⁴ = 279841, so e = 4.
    let v4 = target_v4();
    assert_eq!(v4.minimum_modulus(), &big(138_337));
    let set = ready(run(&v4, generous(1)));
    assert_eq!(primes_of(&set), vec![23]);
    assert_eq!(set.exponent(), 4);
    assert_eq!(set.modulus(), &big(279_841));
    let roots: Vec<u64> = set.contributions().iter().map(|c| c.root()).collect();
    assert_eq!(roots, vec![2, 11, 12, 21]);
    assert_local_identities(&set);
}

#[test]
fn kat_local_roots_satisfy_target_identity() {
    let t = target_527();
    let set = ready(run(&t, generous(1)));
    // H(R) = 24 + 16R mod 49: H(10) = 184 = 3·49 + 37 and H(39) = 648 = 13·49 + 11. Square roots:
    // 24² = 576 = 11·49 + 37 gives ±24 = {24, 25}; 16² = 256 = 5·49 + 11 gives {16, 33}.
    let local: Vec<BigInt> = set
        .contributions()
        .iter()
        .map(|c| c.local_root().clone())
        .collect();
    assert!(
        local[0] == big(24) || local[0] == big(25),
        "s at 10 is {}",
        local[0]
    );
    assert!(
        local[1] == big(16) || local[1] == big(33),
        "s at 39 is {}",
        local[1]
    );
    for (c, value) in set.contributions().iter().zip([37, 11]) {
        assert_eq!(t.target().eval_mod(c.lifted_root(), &big(49)), big(value));
    }
    assert_local_identities(&set);

    // ℓ = 2, e = 1: H ≡ 12, 13 modulo 23 with 9² = 81 ≡ 12 and 6² = 36 ≡ 13; modulo 7 the values
    // 2 and 4 have roots {3, 4} and {2, 5}.
    let set = ready(run(&t, generous(2)));
    let local: Vec<BigInt> = set
        .contributions()
        .iter()
        .map(|c| c.local_root().clone())
        .collect();
    let allowed: [[i64; 2]; 4] = [[3, 4], [2, 5], [9, 14], [6, 17]];
    for (s, pair) in local.iter().zip(allowed) {
        assert!(
            *s == big(pair[0]) || *s == big(pair[1]),
            "s = {s}, expected one of {pair:?}"
        );
    }
    assert_local_identities(&set);
}

#[test]
fn kat_lagrange_basis_identities() {
    // R = 10, 39 modulo 49. L_0 = (t − 39)/(10 − 39) = (t − 39)/20 (−29 ≡ 20), with 20·27 = 540 =
    // 11·49 + 1 so 20⁻¹ = 27, giving L_0 = 27t − 1053 ≡ 27t + 25 (1053 = 21·49 + 24). L_1 =
    // (t − 10)/29 with 29·22 = 638 = 13·49 + 1, giving 22t − 220 ≡ 22t + 25 (220 = 4·49 + 24).
    // Checks: L_0(10) = 295 = 6·49 + 1, L_0(39) = 1078 = 22·49, and L_0 + L_1 = 49t + 50 ≡ 1.
    let set = ready(run(&target_527(), generous(1)));
    let basis: Vec<IntPoly> = set
        .contributions()
        .iter()
        .map(|c| c.lagrange().clone())
        .collect();
    assert_eq!(basis, vec![poly(&[25, 27]), poly(&[25, 22])]);
    assert_local_identities(&set);

    // V4: degree 3 basis over four roots at 23⁴ and at 23²·47², all identities by definition.
    assert_local_identities(&ready(run(&target_v4(), generous(1))));
    let set = ready(run(&target_v4(), generous(2)));
    assert_eq!(set.contributions().len(), 8);
    assert_local_identities(&set);
    assert!(
        set.contributions()
            .iter()
            .all(|c| c.lagrange().degree() == Some(3))
    );

    // d = 1: the empty products make L = 1 exactly.
    let degree_one = target(&[-5, 1], 5, 21, &[(169, 0)]);
    let set = ready(run(&degree_one, generous(3)));
    assert!(
        set.contributions()
            .iter()
            .all(|c| *c.lagrange() == poly(&[1]))
    );
    assert_local_identities(&set);
}

#[test]
fn kat_crt_idempotents() {
    // Q = 7 and 23, M = 161. M/7 = 23 ≡ 2 mod 7 with 2·4 = 8 ≡ 1, so W_7 = 23·4 = 92 (92 = 13·7 + 1,
    // 92 = 4·23). M/23 = 7 with 7·10 = 70 = 3·23 + 1, so W_23 = 70 (70 = 10·7, 70 ≡ 1 mod 23).
    // 92 + 70 = 162 ≡ 1 and 92² = 8464 = 52·161 + 92.
    let set = ready(run(&target_527(), generous(2)));
    let weights: Vec<BigInt> = set.locals().iter().map(|l| l.weight().clone()).collect();
    assert_eq!(weights, vec![big(92), big(70)]);
    assert_local_identities(&set);

    // ℓ = 1: M = Q, the cofactor is 1, and W = 1.
    let set = ready(run(&target_527(), generous(1)));
    assert_eq!(set.locals()[0].weight(), &big(1));

    // d = 1 with primes 5, 11, 17 (3 and 7 divide N = 21, 13 gives H = 169 ≡ 0): M = 935.
    // M/5 = 187 ≡ 2, 2·3 ≡ 1, W_5 = 561 (561 = 33·17 = 51·11, ≡ 1 mod 5). M/11 = 85 ≡ 8,
    // 8·7 = 56 ≡ 1, W_11 = 595 = 35·17 = 119·5, ≡ 1 mod 11. M/17 = 55 ≡ 4, 4·13 = 52 ≡ 1,
    // W_17 = 715 = 143·5 = 65·11, ≡ 1 mod 17. 561 + 595 + 715 = 1871 = 2·935 + 1.
    let degree_one = target(&[-5, 1], 5, 21, &[(169, 0)]);
    let set = ready(run(&degree_one, generous(3)));
    assert_eq!(set.modulus(), &big(935));
    let weights: Vec<BigInt> = set.locals().iter().map(|l| l.weight().clone()).collect();
    assert_eq!(weights, vec![big(561), big(595), big(715)]);
    assert_local_identities(&set);
}

#[test]
fn kat_every_sign_vector_is_a_local_square() {
    // z ≡ Σ ε·s·L has z(R) = ±s at every root, so z(R)² = H(R) modulo Q, and F ≡ Π(t − R) modulo
    // Q makes z² ≡ H modulo (F, Q); the Chinese remainder theorem lifts that to M. This holds for
    // every one of the 2^n sign vectors, not only the one that is a global root. Computed here by
    // polynomial division, not by the code under test.
    assert_every_sign_vector_is_a_local_square(&ready(run(&target_527(), generous(1))));
    assert_every_sign_vector_is_a_local_square(&ready(run(&target_527(), generous(2))));
    assert_every_sign_vector_is_a_local_square(&ready(run(&target_v4(), generous(1))));
    let degree_one = target(&[-5, 1], 5, 21, &[(169, 0)]);
    assert_every_sign_vector_is_a_local_square(&ready(run(&degree_one, generous(3))));
}

#[test]
fn kat_no_sign_is_selected() {
    // Both signs of every contribution are valid local data, and flipping all signs negates the
    // assembled candidate exactly (M is odd, so the centered residue is symmetric).
    let set = ready(run(&target_527(), generous(2)));
    let n = set.contributions().len();
    let plus = set.assemble(&vec![ContributionSign::Plus; n]).unwrap();
    let minus = set.assemble(&vec![ContributionSign::Minus; n]).unwrap();
    assert_eq!(minus, plus.neg());
    // Each contribution is C = W·s·L. Replacing s by Q − s changes C to −C modulo M, because
    // W·Q ≡ 0 modulo every Q_k: the stored sign of s is arbitrary, the pair {C, −C} is not.
    let m = set.modulus();
    for c in set.contributions() {
        let w = set.locals()[c.prime_index()].weight();
        let negated: Vec<BigInt> = c
            .lagrange()
            .coeffs
            .iter()
            .map(|l| {
                (w * (set.locals()[c.prime_index()].modulus() - c.local_root()) * l).mod_floor(m)
            })
            .collect();
        for (stored, other) in c.coefficients().iter().zip(&negated) {
            assert_eq!((stored + other).mod_floor(m), big(0));
        }
    }
}

#[test]
fn kat_known_root_sign_vector_assembles() {
    // z = 4 + 2θ squares to 16 + 16θ + 4θ² = 24 + 16θ = H modulo t² − 2, so z(R) is a square root
    // of H(R) at every lifted root; the signs come from z(R) ≡ ±s alone.
    let z = poly(&[4, 2]);
    let t = target_527();

    // ℓ = 1, e = 2: z(10) = 24 and z(39) = 82 ≡ 33 modulo 49. Hand sum with s = 24, 33: C_0 =
    // 24·(27t + 25) = 648t + 600 ≡ 11t + 12 (648 = 13·49 + 11, 600 = 12·49 + 12) and C_1 =
    // 33·(22t + 25) = 726t + 825 ≡ 40t + 41; C_0 + C_1 = 51t + 53 ≡ 2t + 4.
    let set = ready(run(&t, generous(1)));
    let z_at: Vec<BigInt> = [10, 39]
        .iter()
        .map(|&r| z.eval_mod(&big(r), &big(49)))
        .collect();
    assert_eq!(z_at, vec![big(24), big(33)]);
    let signs = signs_of_known_root(&set, &z);
    let assembled = set.assemble(&signs).unwrap();
    assert!(is_root_up_to_sign(&assembled, &z));
    assert_eq!(count_assemblies_of(&set, &z), 2);
    verified_root_up_to_sign(&set, &z, &assembled);
    let up_to_sign = |actual: &[BigInt], expected: [i64; 2]| {
        let plain = actual == [big(expected[0]), big(expected[1])];
        let negated = actual == [big(49 - expected[0]), big(49 - expected[1])];
        assert!(
            plain || negated,
            "{actual:?} is not ±{expected:?} modulo 49"
        );
    };
    up_to_sign(set.contributions()[0].coefficients(), [12, 11]);
    up_to_sign(set.contributions()[1].coefficients(), [41, 40]);

    // ℓ = 2, e = 1, M = 161: the same z, with 2 of the 16 sign vectors giving ±z.
    let set = ready(run(&t, generous(2)));
    let signs = signs_of_known_root(&set, &z);
    let assembled = set.assemble(&signs).unwrap();
    assert!(is_root_up_to_sign(&assembled, &z));
    assert_eq!(count_assemblies_of(&set, &z), 2);
    verified_root_up_to_sign(&set, &z, &assembled);

    // V4, no inert prime. z = −4 − 20θ + 20θ² + 4θ³ squares to H = −384 + 4224θ² + 768θ³ modulo
    // F: checked here directly. 2 of 16 sign vectors at 23⁴ and 2 of 256 at 23²·47².
    let v4 = target_v4();
    let z = z_v4();
    assert_eq!(
        z.mul(&z).pseudo_div_rem(v4.context().monic_poly()).1,
        *v4.target()
    );
    for max_primes in [1, 2] {
        let set = ready(run(&v4, generous(max_primes)));
        let signs = signs_of_known_root(&set, &z);
        let assembled = set.assemble(&signs).unwrap();
        assert!(is_root_up_to_sign(&assembled, &z));
        assert_eq!(count_assemblies_of(&set, &z), 2);
        match v4.verify_candidate(&assembled, set.method()).unwrap() {
            CandidateVerdict::Verified(root) => {
                // β = z/(c^h F') = 1 + θ, whose image at cm = 4 is 5 with 5² = 25 = P.
                assert_eq!(root.rational_product(), &big(25));
                let beta: Vec<BigRational> = root.beta().coeffs.clone();
                let one = BigRational::from_integer(big(1));
                let one_plus_theta = [one.clone(), one];
                assert!(beta == one_plus_theta || beta == one_plus_theta.map(|c| -c));
                assert_eq!(root.method(), &set.method());
            }
            other => panic!("expected a verified V4 root, got {other:?}"),
        }
    }

    // d = 1, F = t − 5, H = 169 = 13²: primes 5, 11, 17, z = 13, M = 935 ≥ 27. Verified at once.
    let degree_one = target(&[-5, 1], 5, 21, &[(169, 0)]);
    let z = poly(&[13]);
    let set = ready(run(&degree_one, generous(3)));
    assert_eq!(primes_of(&set), vec![5, 11, 17]);
    let assembled = set.assemble(&signs_of_known_root(&set, &z)).unwrap();
    assert!(is_root_up_to_sign(&assembled, &z));
    assert_eq!(count_assemblies_of(&set, &z), 2);
    verified_root_up_to_sign(&set, &z, &assembled);
}

#[test]
fn kat_arbitrary_precision_moduli_beyond_u64() {
    // a = 10²⁰, relation (a, 1) twice, F = t² − 2, m = 23, N = 527: γ = (a − α)², β = a − θ,
    // η = (a − θ)² = a² + 2 − 2aθ, and z = F'β = 2θ(a − θ) = −4 + 2aθ with z² = 16 − 16aθ +
    // 4a²θ² = 8a² + 16 − 16aθ = 8η = H. E = a² + 2 + 6a so U = a + 3 and B = (6U, 2U), which
    // contains |z_1| = 2a ≤ 2a + 6. Hence 2B + 1 = 12a + 37 ≈ 1.2·10²¹ > 2⁶⁴.
    let a = BigInt::from(100_000_000_000_000_000_000u128);
    let pairs = vec![(a.clone(), big(1)), (a.clone(), big(1))];
    let t = target_big(&[-2, 0, 1], 23, 527, &pairs);
    let minimum = &a * 12 + 37;
    assert_eq!(t.minimum_modulus(), &minimum);
    let z = IntPoly::from_coeffs(vec![big(-4), &a * 2]);
    let two_64 = BigInt::from(u64::MAX) + 1;

    // ℓ = 1 at 7 (a ≡ 2 mod 7, not a root 3 or 4): 7²⁴ ≈ 1.9·10²⁰ < 12a + 37 ≤ 7²⁵.
    let set = ready(run(&t, generous(1)));
    assert_eq!(primes_of(&set), vec![7]);
    assert_eq!(set.exponent(), 25);
    assert_eq!(set.modulus(), &pow(big(7), 25));
    assert_eq!(
        set.modulus(),
        &"1341068619663964900807".parse::<BigInt>().unwrap()
    );
    assert!(set.modulus() > &two_64);
    assert_local_identities(&set);
    let assembled = set.assemble(&signs_of_known_root(&set, &z)).unwrap();
    assert!(is_root_up_to_sign(&assembled, &z));
    assert_eq!(count_assemblies_of(&set, &z), 2);
    verified_root_up_to_sign(&set, &z, &assembled);

    // ℓ = 2 at 7 and 23 (a ≡ 3 mod 23, not a root 5 or 18): 161⁹ < 12a + 37 ≤ 161¹⁰ with
    // 161¹⁰ = 11701964070276793473601.
    let set = ready(run(&t, generous(2)));
    assert_eq!(primes_of(&set), vec![7, 23]);
    assert_eq!(set.exponent(), 10);
    assert_eq!(
        set.modulus(),
        &"11701964070276793473601".parse::<BigInt>().unwrap()
    );
    // Each Q_i fits a word (7¹⁰ = 282475249, 23¹⁰ = 41426511213649) but the product M does not.
    assert!(set.modulus() > &two_64);
    assert!(set.locals().iter().all(|l| l.modulus() < &two_64));
    assert_local_identities(&set);
    let assembled = set.assemble(&signs_of_known_root(&set, &z)).unwrap();
    assert!(is_root_up_to_sign(&assembled, &z));
    assert_eq!(count_assemblies_of(&set, &z), 2);
    verified_root_up_to_sign(&set, &z, &assembled);
}

#[test]
fn kat_local_nonresidue_certifies_not_square() {
    // F = t² − 3, m = 4, N = 13 (f(4) = 13), relation (2, −1): γ = 2 + θ, H = F'²η = 4t²(2 + t) =
    // 12(2 + t) = 24 + 12t. Modulo 3, t² is not squarefree; modulo 5 and 7, 3 is not a square
    // (squares {1, 4} and {1, 2, 4}). Modulo 11, 5² = 25 ≡ 3 gives roots 5 and 6, and
    // H(5) = 24 + 60 = 84 = 7·11 + 7: squares modulo 11 are {1, 3, 4, 5, 9}, so 7 is a nonresidue.
    let t = target(&[-3, 0, 1], 4, 13, &[(2, -1)]);
    assert_eq!(t.target(), &poly(&[24, 12]));
    assert_eq!(t.target().eval_mod(&big(5), &big(11)), big(7));
    assert_eq!(
        evidence(run(&t, generous(3))),
        NonSquareEvidence::SplitLocalNonresidue { prime: 11, root: 5 }
    );
    // The first nonresidue ends the search even when more primes are wanted.
    assert_eq!(
        evidence(run(&t, generous(1))),
        NonSquareEvidence::SplitLocalNonresidue { prime: 11, root: 5 }
    );

    // A reducible quotient is not an obstruction to the certificate: f = 3x² + x + 2 over Q(√−23)
    // has c = 3, F = t² + t + 6, m = 7, N = 13 (f(7) = 156 = 12·13), relation (2, 1): η = 3(6 − θ),
    // H = (2t + 1)²(18 − 3t) ≡ 69t − 414 modulo F. 3 | c skips 3; 13 | N splits F (4² + 4 + 6 =
    // 26, 8² + 8 + 6 = 78) with H(4) ≡ 5, a nonresidue, so only the divisibility rule prevents a
    // certificate at 13. 23 divides disc F = −23 and is skipped. 29 has the root 10 with
    // H(10) = 690 − 414 = 276 ≡ 15, and squares modulo 29 do not contain 15.
    let c3 = target(&[2, 1, 3], 7, 13, &[(2, 1)]);
    assert_eq!(c3.target(), &poly(&[-414, 69]));
    assert_eq!(c3.target().eval_mod(&big(4), &big(13)), big(5));
    assert_eq!(
        evidence(run(&c3, generous(3))),
        NonSquareEvidence::SplitLocalNonresidue {
            prime: 29,
            root: 10
        }
    );
}

#[test]
fn kat_local_false_positive_is_not_a_certificate() {
    // The same target with first_candidate 14: 17 (squares {1, 2, 4, 8, 9, 13, 15, 16}) and 19
    // have no square root of 3. Modulo 23, 7² = 49 ≡ 3, so the roots are 7 and 16, and
    // H(7) = 24 + 84 = 108 ≡ 16 = 4², H(16) = 24 + 192 = 216 ≡ 9 = 3²: both squares, so the
    // outcome is Ready. 2 + θ is a norm-one unit that is not a square in Q(√3), so no sign vector
    // is a root: B = (24, 6), 23 < 49 ≤ 23² makes M = 529, and every candidate is rejected.
    let t = target(&[-3, 0, 1], 4, 13, &[(2, -1)]);
    let set = ready(run(&t, limits(14, 200, 1, 100_000, 1000)));
    assert_eq!(primes_of(&set), vec![23]);
    assert_eq!(set.exponent(), 2);
    assert_eq!(set.modulus(), &big(529));
    let roots: Vec<u64> = set.contributions().iter().map(|c| c.root()).collect();
    assert_eq!(roots, vec![7, 16]);
    assert_eq!(t.target().eval_mod(&big(7), &big(23)), big(16));
    assert_eq!(t.target().eval_mod(&big(16), &big(23)), big(9));
    assert_local_identities(&set);
    assert_every_sign_vector_is_a_local_square(&set);
    for mask in 0..4 {
        let candidate = set.assemble(&signs_from_mask(2, mask)).unwrap();
        match t.verify_candidate(&candidate, set.method()).unwrap() {
            CandidateVerdict::Rejected(_) => {}
            other => panic!("mask {mask}: expected a rejection, got {other:?}"),
        }
    }
}

#[test]
fn kat_split_resource_outcomes() {
    let t = target_527();

    // Prime budget: 3 and 5 are inert, so after 1 or 2 examined primes none is usable; a spent
    // budget with no usable prime reports consumption equal to the limit.
    for examined in [0, 1, 2] {
        let e = exhausted(run(&t, limits(0, examined, 1, 100_000, 1000)));
        assert_eq!(
            e,
            ResourceExhaustion {
                kind: ResourceKind::SplitPrimeSearch,
                consumed: examined,
                limit: examined
            }
        );
    }
    // max_primes = 0 stops before any search, whatever the prime budget.
    let e = exhausted(run(&t, limits(0, 100, 0, 100_000, 1000)));
    assert_eq!(
        e,
        ResourceExhaustion {
            kind: ResourceKind::SplitPrimeSearch,
            consumed: 0,
            limit: 0
        }
    );

    // The machine-prime range: 2⁶⁴ − 59 is the largest prime below 2⁶⁴ and ≡ 5 mod 8, so 2 is a
    // nonsquare there and F = t² − 2 is inert; no prime follows.
    let end = u64::MAX - 58;
    let e = exhausted(run(&t, limits(end, 5, 1, 100_000, 1000)));
    assert_eq!(
        e,
        ResourceExhaustion {
            kind: ResourceKind::MachinePrimeRange,
            consumed: 1,
            limit: 5
        }
    );

    // Root search: at 7 the roots 3 and 4 need evaluations of 0 … 4, so 5; a budget of 4 finds only
    // 3 and 0 finds none. Both are truncated.
    let e = exhausted(run(&t, limits(0, 100, 1, 0, 1000)));
    assert_eq!(
        e,
        ResourceExhaustion {
            kind: ResourceKind::RootSearch,
            consumed: 0,
            limit: 0
        }
    );
    let e = exhausted(run(&t, limits(0, 100, 1, 4, 1000)));
    assert_eq!(
        e,
        ResourceExhaustion {
            kind: ResourceKind::RootSearch,
            consumed: 4,
            limit: 4
        }
    );
    assert!(matches!(
        run(&t, limits(0, 100, 1, 5, 1000)),
        SplitOutcome::Ready(_)
    ));

    // Nonresidue search at p ≡ 1 mod 4: F = t² − 2, m = 3, N = 7 (f(3) = 7), relation (3, −2), the
    // first usable prime is 17 (3, 5, 11, 13 are inert, 7 | N) with 6² = 36 ≡ 2 and roots 6, 11.
    // H = 24 + 16t: H(6) = 120 = 7·17 + 1 is 1, a square needing no search. H(11) = 200 = 11·17 +
    // 13 is 8² = 64 = 3·17 + 13, but 17 − 1 = 2⁴ has odd part 1 and 13¹ ≠ 1, so Tonelli–Shanks
    // needs a nonresidue and a budget of 0 cannot supply one.
    let sevens = target(&[-2, 0, 1], 3, 7, &[(3, -2)]);
    let e = exhausted(run(&sevens, limits(0, 100, 1, 100_000, 0)));
    assert_eq!(
        e,
        ResourceExhaustion {
            kind: ResourceKind::NonresidueSearch,
            consumed: 0,
            limit: 0
        }
    );
    let set = ready(run(&sevens, limits(0, 100, 1, 100_000, 1000)));
    assert_eq!(primes_of(&set), vec![17]);
    let roots: Vec<u64> = set.contributions().iter().map(|c| c.root()).collect();
    assert_eq!(roots, vec![6, 11]);
    // 17 < 37 ≤ 17² = 289 gives e = 2.
    assert_eq!(set.modulus(), &big(289));
    assert_local_identities(&set);

    // The precision budget is enforced when the target is built, before any prime search.
    let ctx = SqrtContext::new(&poly(&[-2, 0, 1]), &big(23), &big(527)).unwrap();
    let factors = RelationFactors::new(&ctx, &[(big(3), big(-2))]).unwrap();
    match NormalizedTarget::new(
        factors,
        PrecisionLimits {
            max_modulus_bits: 3,
        },
    )
    .unwrap()
    {
        TargetOutcome::ResourceExhausted(e) => assert_eq!(e.kind, ResourceKind::Precision),
        TargetOutcome::Ready(_) => panic!("2B + 1 = 37 has 6 bits"),
    }
}

#[test]
fn kat_skipped_primes_are_divisors_and_zero_target_values() {
    // Zero H at a root: F = t² − 2, m = 23, N = 527, relation (3, 1) twice, so η = (3 − t)² and
    // H = 8(3 − t)². Modulo 7 the root 3 gives H(3) = 0, so 7 is skipped even though it splits;
    // 17 divides N; 3, 5, 11, 13, 19 are inert. 23 has the roots 5, 18 with H(5) = 8·4 = 32 ≡ 9
    // and H(18) = 8·225 = 1800 ≡ 6 (1800 = 78·23 + 6), both nonzero (3 − 18 = −15 ≢ 0).
    let t = target(&[-2, 0, 1], 23, 527, &[(3, 1), (3, 1)]);
    assert_eq!(t.target().eval_mod(&big(3), &big(7)), big(0));
    let set = ready(run(&t, generous(1)));
    assert_eq!(primes_of(&set), vec![23]);
    let roots: Vec<u64> = set.contributions().iter().map(|c| c.root()).collect();
    assert_eq!(roots, vec![5, 18]);
    assert_eq!(t.target().eval_mod(&big(5), &big(23)), big(9));
    assert_eq!(t.target().eval_mod(&big(18), &big(23)), big(6));
    assert_local_identities(&set);
    // With two primes wanted 7 is still skipped and 31 divides N = 17·31; the second prime is 41,
    // where 17² = 289 = 7·41 + 2 gives the roots 17 and 24, and 3 − 17, 3 − 24 are nonzero.
    assert_eq!(limits_ready_primes(&t, 0, 2), vec![23, 41]);

    // p | N although F splits: 17 divides 527 and x² − 2 has 6² = 36 ≡ 2 there, so starting at
    // 17 the prime is 23 (tested above); and 3 and 7 divide N = 21 for d = 1, F = t − 5, while 13
    // is skipped because H = 169 ≡ 0 at its root 5: the usable primes are 5, 11, 17.
    let degree_one = target(&[-5, 1], 5, 21, &[(169, 0)]);
    let set = ready(run(&degree_one, generous(3)));
    assert_eq!(primes_of(&set), vec![5, 11, 17]);
    let n = set.contributions().len();
    assert_eq!(n, 3);
    assert!(
        set.contributions()
            .iter()
            .all(|c| c.root() == if c.prime() == 5 { 0 } else { 5 })
    );
    // A budget of 5 examines 3, 5, 7, 11, 13 and stops: ℓ = 2 with no exhaustion, M = 55 ≥ 27.
    let set = ready(run(&degree_one, limits(0, 5, 3, 100_000, 1000)));
    assert_eq!(primes_of(&set), vec![5, 11]);
    assert_eq!(set.modulus(), &big(55));
}

#[test]
fn kat_degree_one_every_prime_splits() {
    // F = t − 5 splits modulo every prime, with the single root 5, and the Lagrange polynomial is
    // the empty-product constant 1. So n = ℓ contributions of one coefficient each.
    let t = target(&[-5, 1], 5, 21, &[(169, 0)]);
    for ell in 1..=3u32 {
        let set = ready(run(&t, generous(ell)));
        assert_eq!(set.locals().len(), ell as usize);
        assert_eq!(set.contributions().len(), ell as usize);
        assert!(
            set.contributions()
                .iter()
                .all(|c| c.coefficients().len() == 1)
        );
        assert!(
            set.contributions()
                .iter()
                .all(|c| *c.lagrange() == poly(&[1]))
        );
        assert_local_identities(&set);
    }
    // B = 13, so the minimum modulus is 27: with ℓ = 1 at 5, 5² = 25 < 27 ≤ 5³ = 125 and e = 3.
    let set = ready(run(&t, generous(1)));
    assert_eq!(set.exponent(), 3);
    assert_eq!(set.modulus(), &big(125));
}

#[test]
fn kat_assemble_rejects_wrong_sign_count() {
    let set = ready(run(&target_527(), generous(2)));
    let n = set.contributions().len();
    assert_eq!(n, 4);
    for len in [0, 1, n - 1, n + 1] {
        let result = set.assemble(&vec![ContributionSign::Plus; len]);
        assert!(
            matches!(result, Err(SqrtError::Invariant(_))),
            "length {len}"
        );
    }
    assert!(set.assemble(&vec![ContributionSign::Plus; n]).is_ok());
}
