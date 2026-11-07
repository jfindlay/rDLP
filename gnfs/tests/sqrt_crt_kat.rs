//! Known-answer tests for split-prime local contributions to an algebraic square root, the exact
//! sign recombination that turns them into a verified root, and the method-neutral dispatcher.
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
//!
//! The recombination tests fix `ε_0 = +` to remove the `z ↔ −z` symmetry and leave every other
//! sign free. A *match* is such a sign vector whose assembled coefficient of least bound lies
//! within that bound; the tests enumerate all `2^(n−1)` vectors through `assemble` to count the
//! matches independently of the meet-in-the-middle search. The stored sign of each local root `s`
//! is arbitrary (`s ↔ −s` turns `C` into `−C`), which permutes the matches without changing their
//! number, so the tests pin match counts and never an order.

use gnfs::sqrt::{
    AlgebraicSqrtOutcome, CandidateRejection, CandidateVerdict, ContributionSet, ContributionSign,
    CrtLimits, DirectLimits, FailureClass, IncidentalSource, NonSquareEvidence, NormalizedTarget,
    PrecisionLimits, PrimeSearchLimits, RecombinationLimits, ReconstructionMethod, RelationFactors,
    ResourceExhaustion, ResourceKind, SplitOutcome, SqrtContext, SqrtError, SqrtLimits,
    TargetOutcome, UnsupportedDomain, VerifiedRoot, direct_root, recombine_signs, reconstruct_root,
    reconstruct_root_from_factors, split_contributions, split_prime_root,
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

// ─── sign recombination and dispatch ─────────────────────────────────────────

fn recombination(max_list_entries: u64, max_matches: u64) -> RecombinationLimits {
    RecombinationLimits {
        max_list_entries,
        max_matches,
    }
}

/// Room for every list and match count that the small fixtures can produce.
fn roomy() -> RecombinationLimits {
    recombination(1 << 20, 1 << 20)
}

fn direct_limits(first: u64, examined: u64, nonresidue: u64) -> DirectLimits {
    DirectLimits {
        prime_search: PrimeSearchLimits {
            first_candidate: first,
            max_primes_examined: examined,
        },
        nonresidue_search: NonresidueSearchLimits {
            max_candidates: nonresidue,
        },
    }
}

fn sqrt_limits(
    direct: DirectLimits,
    crt: CrtLimits,
    recombination: RecombinationLimits,
) -> SqrtLimits {
    SqrtLimits {
        precision: PrecisionLimits {
            max_modulus_bits: 4096,
        },
        direct,
        crt,
        recombination,
    }
}

/// Direct reconstruction from the first prime; split primes, one at a time, as the fallback.
fn generous_sqrt() -> SqrtLimits {
    sqrt_limits(direct_limits(0, 100, 1000), generous(1), roomy())
}

fn recombined(set: &ContributionSet, limits: RecombinationLimits) -> AlgebraicSqrtOutcome {
    recombine_signs(set, limits).expect("no invariant violation")
}

fn split_root(
    target: &NormalizedTarget,
    crt: CrtLimits,
    limits: RecombinationLimits,
) -> AlgebraicSqrtOutcome {
    split_prime_root(target, crt, limits).expect("no invariant violation")
}

fn dispatched(target: &NormalizedTarget, limits: SqrtLimits) -> AlgebraicSqrtOutcome {
    reconstruct_root(target, limits).expect("no invariant violation")
}

fn verified(outcome: AlgebraicSqrtOutcome) -> VerifiedRoot {
    match outcome {
        AlgebraicSqrtOutcome::Verified(root) => root,
        other => panic!("expected a verified root, got {other:?}"),
    }
}

fn resources(outcome: AlgebraicSqrtOutcome) -> Vec<ResourceExhaustion> {
    match outcome {
        AlgebraicSqrtOutcome::ResourceExhausted(list) => list,
        other => panic!("expected resource exhaustion, got {other:?}"),
    }
}

fn nonsquare(outcome: AlgebraicSqrtOutcome) -> NonSquareEvidence {
    match outcome {
        AlgebraicSqrtOutcome::NotSquare(certificate) => certificate.evidence().clone(),
        other => panic!("expected a nonsquare certificate, got {other:?}"),
    }
}

fn exhaustion(kind: ResourceKind, consumed: u64, limit: u64) -> ResourceExhaustion {
    ResourceExhaustion {
        kind,
        consumed,
        limit,
    }
}

fn split_method(primes: &[u64], exponent: u32) -> ReconstructionMethod {
    ReconstructionMethod::SplitPrimeCrt {
        primes: primes.to_vec(),
        exponent,
    }
}

fn assert_root_up_to_sign(root: &VerifiedRoot, expected: &[i64]) {
    let z = root.integral_root();
    assert!(
        is_root_up_to_sign(z, &poly(expected)),
        "z = {z:?} is neither ±{expected:?}"
    );
}

/// The modular image is `y` or `N − y`: `β` and `−β` give the same pair of gcds.
fn assert_image_up_to_sign(root: &VerifiedRoot, n: i64, y: i64) {
    let image = root.modular_image();
    assert!(
        *image == big(y) || *image == big(n - y),
        "image {image} is neither {y} nor {}",
        n - y
    );
}

/// The index of the least coefficient bound (lowest index on ties) and that bound.
fn narrowest(target: &NormalizedTarget) -> (usize, BigInt) {
    let bounds = target.coefficient_bounds();
    let kappa = (0..bounds.len())
        .min_by_key(|&i| (bounds[i].clone(), i))
        .expect("a target has at least one coefficient");
    (kappa, bounds[kappa].clone())
}

/// Every sign vector with `ε_0 = +` whose assembled coefficient of least bound lies within that
/// bound, found by assembling all `2^(n−1)` vectors with `assemble` alone.
fn matching_vectors(set: &ContributionSet) -> Vec<Vec<ContributionSign>> {
    let (kappa, bound) = narrowest(set.target());
    let n = set.contributions().len();
    (0..1u32 << (n - 1))
        .map(|mask| signs_from_mask(n, mask << 1))
        .filter(|signs| {
            let z = set.assemble(signs).expect("length matches");
            let value = z.coeffs.get(kappa).cloned().unwrap_or_else(|| big(0));
            -&bound <= value && value <= bound
        })
        .collect()
}

/// `signs` negated as a whole when its first sign is `−`, so that `ε_0 = +`.
fn with_first_sign_plus(signs: Vec<ContributionSign>) -> Vec<ContributionSign> {
    if signs[0] == ContributionSign::Plus {
        return signs;
    }
    signs
        .into_iter()
        .map(|s| match s {
            ContributionSign::Plus => ContributionSign::Minus,
            ContributionSign::Minus => ContributionSign::Plus,
        })
        .collect()
}

/// Whether the residue interval `[−x−b, −x+b] (mod M)` of this vector's left-half sum `x` wraps
/// past `M`. The left half is contributions `1..=⌊(n−1)/2⌋`, with contribution 0 added to `x`.
fn left_interval_wraps(set: &ContributionSet, signs: &[ContributionSign]) -> bool {
    let (kappa, bound) = narrowest(set.target());
    let modulus = set.modulus();
    let column = |k: usize| set.contributions()[k].coefficients()[kappa].clone();
    let mut x = column(0);
    for (k, sign) in signs
        .iter()
        .enumerate()
        .take((signs.len() - 1) / 2 + 1)
        .skip(1)
    {
        match sign {
            ContributionSign::Plus => x += column(k),
            ContributionSign::Minus => x -= column(k),
        }
    }
    let x = x.mod_floor(modulus);
    let lo = (-&x - &bound).mod_floor(modulus);
    let hi = (&bound - &x).mod_floor(modulus);
    lo > hi
}

/// `x² − 2`, `m = 23`, `N = 527`, relation `(−5, 0)`: `γ = −5`, which is negative and so not a
/// square in the real field `Q(√2)`. `η = −5`, `F' = 2t` and `H = 4t²·(−5) = −40` modulo
/// `t² − 2`. Modulo 7, `−40 ≡ 2 = 3²` and `t² = 2` has the roots 3, 4; modulo 23, `−40 ≡ 6 = 11²`
/// with roots 5, 18. So `H` is a nonzero square at every root, with no local obstruction. The
/// bound is `B = (18, 6)` as for `(3, −2)`, giving `M = 49` at one prime and `161` at two.
fn target_negative_five() -> NormalizedTarget {
    target(&[-2, 0, 1], 23, 527, &[(-5, 0)])
}

/// `x − 5`, `m = 5`, `N = 21`, relation `(169, 0)`: `F = t − 5`, `H = 169 = 13²`, `B = 13`.
fn target_degree_one() -> NormalizedTarget {
    target(&[-5, 1], 5, 21, &[(169, 0)])
}

/// `x² − 3`, `m = 4`, `N = 13`, relation `(2, −1)`, with split primes searched from 14: `γ = 2 + θ`
/// is a nonsquare unit of `Q(√3)` that is a square at the prime 23.
fn target_false_positive() -> NormalizedTarget {
    target(&[-3, 0, 1], 4, 13, &[(2, -1)])
}

fn false_positive_limits() -> CrtLimits {
    limits(14, 200, 1, 100_000, 1000)
}

#[test]
fn kat_single_prime_recombination_verifies_a_root() {
    // x² − 2, N = 527, relation (3, −2), z = 4 + 2θ, ℓ = 1 at p = 7: M = 49 and n = 2, so the left
    // half is empty and the right half has the one free sign. The two contributions are
    // C_0 = 11t + 12 and C_1 = 40t + 41 up to sign (see the lifted-root KAT), so the sign vectors
    // (+, ±) have coefficient 1 equal to 11 ± 40: 51 ≡ 2 and −29 ≡ 20 modulo 49. Only 2 is within
    // b = 6, and it is the coefficient of z.
    let set = ready(run(&target_527(), generous(1)));
    assert_eq!(matching_vectors(&set).len(), 1);
    let root = verified(recombined(&set, roomy()));
    assert_root_up_to_sign(&root, &[4, 2]);
    assert_eq!(*root.method(), split_method(&[7], 2));
    assert_eq!(root.rational_product(), &big(49));
    // β = ±(1 + θ), y = β(23) = 24 and 24² = 576 = 527 + 49.
    assert_image_up_to_sign(&root, 527, 24);

    // d = 1, ℓ = 1: F = t − 5, p = 5, e = 3, M = 125, n = 1. There is no free sign at all: the
    // candidate is C_0 = s with s² ≡ 169 (mod 125), that is s ∈ {13, 112}, centered ±13.
    let degree_one = target_degree_one();
    let set = ready(run(&degree_one, generous(1)));
    assert_eq!(set.contributions().len(), 1);
    assert_eq!(set.modulus(), &big(125));
    assert_eq!(matching_vectors(&set).len(), 1);
    // One list entry suffices when there is no free sign, and a single match.
    let root = verified(recombined(&set, recombination(1, 1)));
    assert_root_up_to_sign(&root, &[13]);
    assert_eq!(*root.method(), split_method(&[5], 3));
}

#[test]
fn kat_multiple_primes_recombine_without_cross_prime_sign_normalisation() {
    // x² − 2, N = 527, z = 4 + 2θ at primes 7 and 23: M = 161, n = 4. The true vector assigns one
    // sign to each of the four contributions, and the search may not assume, for instance, that
    // the two primes agree: only ε_0 is fixed. It is the only vector with ε_0 = + that matches.
    let z = poly(&[4, 2]);
    let set = ready(run(&target_527(), generous(2)));
    let known = with_first_sign_plus(signs_of_known_root(&set, &z));
    let matches = matching_vectors(&set);
    assert_eq!(matches, vec![known]);
    let root = verified(recombined(&set, roomy()));
    assert_root_up_to_sign(&root, &[4, 2]);
    assert_eq!(*root.method(), split_method(&[7, 23], 1));

    // d = 1 with the primes 5, 11, 17: M = 935, n = 3, z = ±13. The one match is the true vector.
    let set = ready(run(&target_degree_one(), generous(3)));
    let known = with_first_sign_plus(signs_of_known_root(&set, &poly(&[13])));
    assert_eq!(matching_vectors(&set), vec![known]);
    let root = verified(recombined(&set, roomy()));
    assert_root_up_to_sign(&root, &[13]);
    assert_eq!(*root.method(), split_method(&[5, 11, 17], 1));

    // V4 at 23 and 47: M = 23²·47² = 1168561 and n = 8 with 3 matches, one of which is the root.
    let set = ready(run(&target_v4(), generous(2)));
    let known = with_first_sign_plus(signs_of_known_root(&set, &z_v4()));
    assert!(matching_vectors(&set).contains(&known));
    let root = verified(recombined(&set, roomy()));
    assert_root_up_to_sign(&root, &[-4, -20, 20, 4]);
    assert_eq!(*root.method(), split_method(&[23, 47], 2));
}

#[test]
fn kat_v4_quartic_without_an_inert_prime_reconstructs_by_split_primes() {
    // F = t⁴ − 10t² + 1 has no inert prime, so the direct method can only exhaust its budget: it
    // examines 50 primes and every one is a reducible reduction.
    let v4 = target_v4();
    let direct = direct_root(&v4, direct_limits(0, 50, 1000)).expect("no invariant violation");
    assert_eq!(
        resources(direct),
        vec![exhaustion(ResourceKind::InertPrimeSearch, 50, 50)]
    );
    // The split-prime method needs no inert prime. z = −4 − 20θ + 20θ² + 4θ³, β = 1 + θ, whose
    // image at cm = 4 is 5, and 5² = 25 = P = (1 + 1·4)·(1 + 1·4) (relation (1, −1) twice, a − bm).
    for (ell, primes, exponent) in [(1, vec![23], 4), (2, vec![23, 47], 2)] {
        let root = verified(split_root(&v4, generous(ell), roomy()));
        assert_root_up_to_sign(&root, &[-4, -20, 20, 4]);
        assert_eq!(*root.method(), split_method(&primes, exponent));
        assert_eq!(root.rational_product(), &big(25));
        assert_image_up_to_sign(&root, 97, 5);
        let one = BigRational::from_integer(big(1));
        let beta = root.beta().coeffs.clone();
        let one_plus_theta = vec![one.clone(), one];
        assert!(
            beta == one_plus_theta || beta == one_plus_theta.iter().map(|c| -c).collect::<Vec<_>>()
        );
    }
}

#[test]
fn kat_close_matches_are_not_roots() {
    // A match is close on one coefficient only. For (−5, 0) over x² − 2 the candidates are local
    // squares at every root yet never a global square: γ = −5 < 0. At ℓ = 1 the sole match is
    // z = ±3 with z² = 9 ≡ −40 = H (mod 49) but 9 ≠ −40 over Z. At ℓ = 2 (M = 161) the two
    // matches are z = ±(80, 0), whose first coefficient exceeds B_0 = 18, and ±(−11, 0), with
    // 121 ≡ −40 (mod 161) but 121 ≠ −40.
    let t = target_negative_five();
    let rejection_of = |set: &ContributionSet, signs: &[ContributionSign]| {
        let candidate = set.assemble(signs).expect("length matches");
        match t
            .verify_candidate(&candidate, set.method())
            .expect("no invariant violation")
        {
            CandidateVerdict::Rejected(reason) => reason,
            other => panic!("expected a rejection, got {other:?}"),
        }
    };
    let set = ready(run(&t, generous(1)));
    let found = matching_vectors(&set);
    assert_eq!(found.len(), 1);
    assert_eq!(
        rejection_of(&set, &found[0]),
        CandidateRejection::NotSquareRootOfTarget
    );
    let set = ready(run(&t, generous(2)));
    let found = matching_vectors(&set);
    assert_eq!(found.len(), 2);
    let reasons: Vec<CandidateRejection> = found.iter().map(|s| rejection_of(&set, s)).collect();
    assert!(reasons.contains(&CandidateRejection::CoefficientOutOfBound { index: 0 }));
    assert!(reasons.contains(&CandidateRejection::NotSquareRootOfTarget));

    // In the V4 field a close match can also be rejected while a genuine root is present: of the
    // 2 matches at 23⁴ and the 3 at 23²·47², exactly one is the root, and the search still finds
    // it whatever order the rejected ones come in.
    let v4 = target_v4();
    for (ell, matches) in [(1, 2), (2, 3)] {
        let set = ready(run(&v4, generous(ell)));
        let found = matching_vectors(&set);
        assert_eq!(found.len(), matches);
        let roots = found
            .iter()
            .filter(|signs| {
                let candidate = set.assemble(signs).expect("length matches");
                matches!(
                    v4.verify_candidate(&candidate, set.method())
                        .expect("no invariant violation"),
                    CandidateVerdict::Verified(_)
                )
            })
            .count();
        assert_eq!(roots, 1);
        assert_root_up_to_sign(&verified(recombined(&set, roomy())), &[-4, -20, 20, 4]);
    }
}

#[test]
fn kat_complete_search_with_every_match_rejected_certifies_nonsquare() {
    // γ = −5 is not a square, and its matches are all rejected: the search is complete, so it is a
    // certificate. The primes and exponent recorded are those searched.
    let t = target_negative_five();
    for (ell, primes, exponent) in [(1, vec![7], 2), (2, vec![7, 23], 1)] {
        let set = ready(run(&t, generous(ell)));
        let expected = NonSquareEvidence::SignSearchExhausted { primes, exponent };
        assert_eq!(nonsquare(recombined(&set, roomy())), expected);
        assert_eq!(nonsquare(split_root(&t, generous(ell), roomy())), expected);
    }
}

#[test]
fn kat_zero_match_complete_search_certifies_nonsquare() {
    // x² − 3 with γ = 2 + θ, prime 23 (first_candidate 14): 2 + θ = (1 + θ)²/2, so it is a square
    // in Q(√3) only if 2 is, and √2 ∉ Q(√3). B = (24, 6) so b = 6 at κ = 1, and M = 529. The two
    // vectors (+, ±) assemble, centered modulo 529 and up to sign, to 61 − 156θ and 126 + 126θ,
    // whose θ-coefficients 156 and 126 both exceed 6: the match list is empty. A search that
    // examined nothing is nonetheless complete, even with no match budget at all and the least
    // list budget (the right half has one sign, so two entries).
    let t = target_false_positive();
    let set = ready(run(&t, false_positive_limits()));
    assert_eq!(set.modulus(), &big(529));
    assert!(matching_vectors(&set).is_empty());
    let expected = NonSquareEvidence::SignSearchExhausted {
        primes: vec![23],
        exponent: 2,
    };
    assert_eq!(nonsquare(recombined(&set, recombination(2, 0))), expected);
    let outcome = split_root(&t, false_positive_limits(), recombination(2, 0));
    assert_eq!(outcome.class(), Some(FailureClass::AlgebraicNonSquare));
    assert_eq!(nonsquare(outcome), expected);
}

#[test]
fn kat_list_budget_exhaustion_is_a_resource_outcome() {
    // V4 at two primes: n = 8, seven free signs, left half 3 and right half 4, so the larger list
    // has 2⁴ = 16 entries. 15 is refused and reports the 16 required; 16 suffices.
    let set = ready(run(&target_v4(), generous(2)));
    let outcome = recombined(&set, recombination(15, u64::MAX));
    assert_eq!(outcome.class(), Some(FailureClass::ResourceExhausted));
    assert_eq!(
        resources(outcome),
        vec![exhaustion(ResourceKind::RecombinationList, 16, 15)]
    );
    assert_root_up_to_sign(
        &verified(recombined(&set, recombination(16, 1 << 20))),
        &[-4, -20, 20, 4],
    );
    // The list size is checked before any match is considered.
    assert_eq!(
        resources(recombined(&set, recombination(0, 0))),
        vec![exhaustion(ResourceKind::RecombinationList, 16, 0)]
    );
    // V4 at one prime: n = 4, right half of 2 signs, 4 entries.
    let set = ready(run(&target_v4(), generous(1)));
    assert_eq!(
        resources(recombined(&set, recombination(3, 1 << 20))),
        vec![exhaustion(ResourceKind::RecombinationList, 4, 3)]
    );
    // n = 1 still needs its one entry, so a zero budget is refused.
    let set = ready(run(&target_degree_one(), generous(1)));
    assert_eq!(
        resources(recombined(&set, recombination(0, 1 << 20))),
        vec![exhaustion(ResourceKind::RecombinationList, 1, 0)]
    );
    // The false-positive fixture is not a certificate when its list does not fit.
    let set = ready(run(&target_false_positive(), false_positive_limits()));
    assert_eq!(
        resources(recombined(&set, recombination(1, 0))),
        vec![exhaustion(ResourceKind::RecombinationList, 2, 1)]
    );

    // A huge sign space is refused without allocating: d = 1 with 130 split primes has n = 130,
    // a right half of 65 signs, and 2⁶⁵ > u64::MAX entries, which saturates the count even when
    // the budget is u64::MAX. 3, 7 divide N = 21 and 13 gives H = 169 ≡ 0, so 130 usable primes
    // are found among the first 139 primes above 3.
    let crowd = limits(0, 400, 130, 100_000, 1000);
    let outcome = split_root(
        &target_degree_one(),
        crowd,
        recombination(u64::MAX, u64::MAX),
    );
    assert_eq!(
        resources(outcome),
        vec![exhaustion(
            ResourceKind::RecombinationList,
            u64::MAX,
            u64::MAX
        )]
    );
}

#[test]
fn kat_match_budget_exhaustion_is_a_resource_outcome() {
    // γ = −5 has 1 match at one prime and 2 at two. A budget equal to the count completes the
    // search (NotSquare); one fewer stops before examining the last match, which is a resource
    // outcome reporting the examined count, not a certificate.
    let t = target_negative_five();
    for (ell, matches) in [(1, 1u64), (2, 2)] {
        let set = ready(run(&t, generous(ell)));
        assert_eq!(matching_vectors(&set).len() as u64, matches);
        assert!(matches!(
            nonsquare(recombined(&set, recombination(1 << 20, matches))),
            NonSquareEvidence::SignSearchExhausted { .. }
        ));
        for limit in 0..matches {
            let outcome = recombined(&set, recombination(1 << 20, limit));
            assert_eq!(outcome.class(), Some(FailureClass::ResourceExhausted));
            assert_eq!(
                resources(outcome),
                vec![exhaustion(ResourceKind::RecombinationMatches, limit, limit)]
            );
        }
    }

    // V4 at two primes has 3 matches, one of them the root. With no budget the first match
    // already exceeds it; with room for all three the root is found. In between, the answer is
    // either the root or an exhaustion that consumed exactly its budget, and a larger budget never
    // loses a root that a smaller one found.
    let set = ready(run(&target_v4(), generous(2)));
    assert_eq!(
        resources(recombined(&set, recombination(1 << 20, 0))),
        vec![exhaustion(ResourceKind::RecombinationMatches, 0, 0)]
    );
    let mut found = false;
    for limit in 1..=3 {
        match recombined(&set, recombination(1 << 20, limit)) {
            AlgebraicSqrtOutcome::Verified(root) => {
                assert_root_up_to_sign(&root, &[-4, -20, 20, 4]);
                found = true;
            }
            AlgebraicSqrtOutcome::ResourceExhausted(list) => {
                assert!(!found, "a larger budget lost the root");
                assert_eq!(
                    list,
                    vec![exhaustion(ResourceKind::RecombinationMatches, limit, limit)]
                );
            }
            other => panic!("limit {limit}: {other:?}"),
        }
    }
    assert!(found, "three matches are room for the root");
}

#[test]
fn kat_wrapped_intervals_are_scanned_exactly() {
    // γ = −5 at 7 and 23: M = 161, b = 6, n = 4, so the left half is the one contribution 1. H = −40
    // is constant, so s is the same up to sign at both roots of a prime, and the t-coefficients of
    // the two Lagrange polynomials of one prime are opposite: for one choice of the sign the two
    // contributions of the prime 7 cancel at κ = 1. The left sum is then 0 and its interval
    // [−6, 6] wraps past M as [155, 160] together with [0, 6].
    let set = ready(run(&target_negative_five(), generous(2)));
    let found = matching_vectors(&set);
    assert_eq!(found.len(), 2);
    assert!(
        found.iter().any(|signs| left_interval_wraps(&set, signs)),
        "a match must lie in a wrapped interval for this test to mean anything"
    );
    // Both matches are found, no more and no fewer: a scan of only [lo, M − 1] or only [lo, hi]
    // would report fewer than 2 and complete the search early.
    assert!(matches!(
        nonsquare(recombined(&set, recombination(1 << 20, 2))),
        NonSquareEvidence::SignSearchExhausted { .. }
    ));
    assert_eq!(
        resources(recombined(&set, recombination(1 << 20, 1))),
        vec![exhaustion(ResourceKind::RecombinationMatches, 1, 1)]
    );
}

#[test]
fn kat_direct_and_split_prime_roots_agree_up_to_sign() {
    // Both methods return through the same verifier, so each result is a root of the same target;
    // since a root is determined up to sign at either kind of prime, they agree up to sign. For
    // x² − 2: direct at the inert prime 3 (e = 4) and split at 7 (e = 2). For x − 5: both at 5.
    for (t, n) in [(target_527(), 527), (target_degree_one(), 21)] {
        let primary =
            verified(direct_root(&t, direct_limits(0, 100, 1000)).expect("no invariant violation"));
        let fallback = verified(split_root(&t, generous(1), roomy()));
        assert!(matches!(
            primary.method(),
            ReconstructionMethod::DirectNewton { .. }
        ));
        assert!(matches!(
            fallback.method(),
            ReconstructionMethod::SplitPrimeCrt { .. }
        ));
        assert!(is_root_up_to_sign(
            fallback.integral_root(),
            primary.integral_root()
        ));
        assert!(
            fallback.beta() == primary.beta() || *fallback.beta() == primary.beta().neg(),
            "β differs beyond sign"
        );
        let y = primary.modular_image();
        let other = fallback.modular_image();
        assert!(other == y || *other == big(n) - y);
        assert_eq!(primary.rational_product(), fallback.rational_product());
    }
}

#[test]
fn kat_dispatch_prefers_the_direct_method() {
    // x² − 2 at N = 527: the first prime 3 is inert, and 3⁴ is the least power above 2B = 36.
    let root = verified(dispatched(&target_527(), generous_sqrt()));
    assert_eq!(
        *root.method(),
        ReconstructionMethod::DirectNewton {
            prime: 3,
            exponent: 4
        }
    );
    assert_root_up_to_sign(&root, &[4, 2]);
}

#[test]
fn kat_dispatch_falls_back_when_the_inert_prime_search_is_exhausted() {
    // Starting at 7 with one prime examined, the direct method sees only 7, at which t² − 2 splits.
    // That is a spent inert-prime budget, not a statement about inert primes, so the split-prime
    // method runs and finds z = ±(4 + 2θ) at 7.
    let limits = sqrt_limits(direct_limits(7, 1, 1000), generous(1), roomy());
    let t = target_527();
    let primary = direct_root(&t, limits.direct).expect("no invariant violation");
    assert_eq!(
        resources(primary),
        vec![exhaustion(ResourceKind::InertPrimeSearch, 1, 1)]
    );
    let root = verified(dispatched(&t, limits));
    assert_eq!(*root.method(), split_method(&[7], 2));
    assert_root_up_to_sign(&root, &[4, 2]);

    // The end of the machine-prime range is the other fallback trigger: no prime above 2⁶⁴ − 59,
    // so a search from just above it ends at once without spending its budget.
    let limits = sqrt_limits(direct_limits(u64::MAX - 57, 5, 1000), generous(1), roomy());
    let primary = direct_root(&t, limits.direct).expect("no invariant violation");
    assert_eq!(
        resources(primary),
        vec![exhaustion(ResourceKind::MachinePrimeRange, 0, 5)]
    );
    let root = verified(dispatched(&t, limits));
    assert_eq!(*root.method(), split_method(&[7], 2));

    // V4 has no inert prime at all, so with the generous direct budget the fallback is the only
    // way to a root, at either number of primes.
    for (ell, primes, exponent) in [(1, vec![23], 4), (2, vec![23, 47], 2)] {
        let limits = sqrt_limits(direct_limits(0, 50, 1000), generous(ell), roomy());
        let root = verified(dispatched(&target_v4(), limits));
        assert_eq!(*root.method(), split_method(&primes, exponent));
        assert_root_up_to_sign(&root, &[-4, -20, 20, 4]);
    }
}

#[test]
fn kat_dispatch_reports_primary_and_fallback_exhaustion_in_order() {
    let v4 = target_v4();
    let spent = exhaustion(ResourceKind::InertPrimeSearch, 50, 50);
    let direct = direct_limits(0, 50, 1000);

    // The fallback finds no split prime at all: no prime budget for it.
    let outcome = dispatched(&v4, sqrt_limits(direct, generous(0), roomy()));
    assert_eq!(
        resources(outcome),
        vec![spent, exhaustion(ResourceKind::SplitPrimeSearch, 0, 0)]
    );
    // The fallback's root search is truncated at 23, the first split prime.
    let outcome = dispatched(
        &v4,
        sqrt_limits(direct, limits(0, 200, 1, 0, 1000), roomy()),
    );
    assert_eq!(
        resources(outcome),
        vec![spent, exhaustion(ResourceKind::RootSearch, 0, 0)]
    );
    // The fallback's sign lists do not fit (V4 at one prime needs 4 entries), or its matches
    // run out (it has 2 and the first one already exceeds a budget of 0).
    let outcome = dispatched(
        &v4,
        sqrt_limits(direct, generous(1), recombination(3, 1 << 20)),
    );
    assert_eq!(
        resources(outcome),
        vec![spent, exhaustion(ResourceKind::RecombinationList, 4, 3)]
    );
    let outcome = dispatched(
        &v4,
        sqrt_limits(direct, generous(1), recombination(1 << 20, 0)),
    );
    assert_eq!(
        resources(outcome),
        vec![spent, exhaustion(ResourceKind::RecombinationMatches, 0, 0)]
    );
    // A machine-range end in the primary comes first as well.
    let range = direct_limits(u64::MAX - 57, 5, 1000);
    let outcome = dispatched(&target_527(), sqrt_limits(range, generous(0), roomy()));
    assert_eq!(
        resources(outcome),
        vec![
            exhaustion(ResourceKind::MachinePrimeRange, 0, 5),
            exhaustion(ResourceKind::SplitPrimeSearch, 0, 0)
        ]
    );
}

#[test]
fn kat_dispatch_never_falls_back_on_other_outcomes() {
    // The split-prime budget is zero, so a fallback would show up as a SplitPrimeSearch
    // exhaustion; every outcome below is the primary's own, returned unchanged.
    let limits = sqrt_limits(direct_limits(0, 100, 1000), generous(0), roomy());

    // A local nonresidue at the inert prime 3: x² − 2, m = 4, N = 7, relation (1, 1), where
    // γ = 1 − θ has norm −1 and H ≡ 2 + θ has norm 2 ≡ −1 in F_9, a nonresidue.
    let t = target(&[-2, 0, 1], 4, 7, &[(1, 1)]);
    assert_eq!(
        nonsquare(dispatched(&t, limits)),
        NonSquareEvidence::InertLocalNonresidue { prime: 3 }
    );
    // A global nonsquare found only at the full bound: relation (3, 1), γ = 3 − θ of norm 7,
    // H ≡ θ is a square in F_9; B = (18, 6), 3⁴ = 81 ≥ 37 gives e = 4.
    let t = target(&[-2, 0, 1], 4, 7, &[(3, 1)]);
    assert_eq!(
        nonsquare(dispatched(&t, limits)),
        NonSquareEvidence::DirectBoundExhausted {
            prime: 3,
            exponent: 4
        }
    );
    // A nonresidue budget of 0 cannot find the nonresidue F_9 needs for H ≡ θ at x² − 2, N = 527.
    let spent_nonresidue = sqrt_limits(direct_limits(0, 100, 0), generous(0), roomy());
    assert_eq!(
        resources(dispatched(&target_527(), spent_nonresidue)),
        vec![exhaustion(ResourceKind::NonresidueSearch, 0, 0)]
    );
    // A root exists but the map fails: x² − 8, relation (2, 0), z = 8 and β = θ/2, whose only
    // denominator 2 is not a unit modulo N. N = 4 (m = 2, f(2) = −4): gcd(2, 4) = 2 is a proper
    // divisor, an incidental factor 4 = 2·2. N = 2 (m = 0, f(0) = −8): gcd(2, 2) = N, unsupported.
    let incidental = target(&[-8, 0, 1], 2, 4, &[(2, 0)]);
    match dispatched(&incidental, limits) {
        AlgebraicSqrtOutcome::Incidental(found) => {
            assert_eq!(found.divisor(), &big(2));
            assert_eq!(found.cofactor(), &big(2));
            assert_eq!(found.source(), IncidentalSource::MapDenominator);
        }
        other => panic!("expected an incidental factor, got {other:?}"),
    }
    let unsupported = target(&[-8, 0, 1], 0, 2, &[(2, 0)]);
    assert_eq!(
        dispatched(&unsupported, limits),
        AlgebraicSqrtOutcome::Unsupported(UnsupportedDomain::MapDenominatorNotUnit)
    );
    // The split-prime method reaches the same verdicts through its own matching: t² − 8 splits at
    // 7 (8 ≡ 1, roots ±1), where H = 64 ≡ 1 is a nonzero square.
    assert!(matches!(
        split_root(&incidental, generous(1), roomy()),
        AlgebraicSqrtOutcome::Incidental(_)
    ));
    assert_eq!(
        split_root(&unsupported, generous(1), roomy()),
        AlgebraicSqrtOutcome::Unsupported(UnsupportedDomain::MapDenominatorNotUnit)
    );
}

#[test]
fn kat_dispatch_precision_budget_applies_before_any_search() {
    // The minimum modulus 2B + 1 = 37 has 6 bits. A dispatcher budget of 5 refuses a target that
    // was built under a looser one, and builds none from factors.
    let mut tight = generous_sqrt();
    tight.precision = PrecisionLimits {
        max_modulus_bits: 5,
    };
    assert_eq!(
        resources(dispatched(&target_527(), tight)),
        vec![exhaustion(ResourceKind::Precision, 6, 5)]
    );
    let context = SqrtContext::new(&poly(&[-2, 0, 1]), &big(23), &big(527)).unwrap();
    let factors = RelationFactors::new(&context, &[(big(3), big(-2))]).unwrap();
    assert_eq!(
        resources(reconstruct_root_from_factors(factors.clone(), tight).unwrap()),
        vec![exhaustion(ResourceKind::Precision, 6, 5)]
    );
    // At exactly 6 bits the target builds and the direct method runs.
    let mut exact = generous_sqrt();
    exact.precision = PrecisionLimits {
        max_modulus_bits: 6,
    };
    let root = verified(reconstruct_root_from_factors(factors, exact).unwrap());
    assert_eq!(
        *root.method(),
        ReconstructionMethod::DirectNewton {
            prime: 3,
            exponent: 4
        }
    );
    assert_root_up_to_sign(&root, &[4, 2]);
}
