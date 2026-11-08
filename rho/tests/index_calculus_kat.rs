//! Known-answer tests (KATs) for the index-calculus ECDLP solver over `E(F_p)`.
//!
//! # Fixture
//!
//! All tests use the `semaev_toy()` curve: `y² = x³ + x + 33` over `F_47` (`p = 47`,
//! `a = 1`, `b = 33`). The generator `G = (10, 3)` has group order `n = 60 = 2²·3·5`.
//! The prime-order subgroup uses `ℓ = 5` (the largest prime factor of n); the ℓ-order
//! subgroup generator is `G_ℓ = (n/ℓ)·G = 12·G`.
//!
//! # KAT coverage — factor base and relation types
//!
//! 1. **factor-base-on-curve**: every enumerated factor-base point satisfies `is_on_curve`
//!    via the frozen `Curve`. Guards the QR-test + canonical-root-lift enumeration.
//! 2. **subgroup-validity**: `ℓ | n` (i.e., `n % ℓ == 0`); the subgroup generator
//!    `G_ℓ = 12·G` satisfies `ℓ·G_ℓ = ∞` (the identity). Guards the prime-order-subgroup
//!    precondition for the `Z/ℓℤ` linear algebra.
//! 3. **relation-round-trip**: construct a `Relation` from a known decomposition
//!    (factor-base points 0 and 1), then verify that `Σ e_i·P_i` (using the frozen group
//!    law) equals the expected point sum. Guards the exponent-vector encoding and the
//!    `from_decomposition` / `exponent` API.
//!
//! # KAT coverage — point decomposition via Semaev polynomial
//!
//! 4. **decomposition-correctness**: `decompose(Q)` returns `Some(pts)` for a point `Q`
//!    that is a sum of two factor-base points; the returned decomposition sums to `Q` via
//!    the frozen group law.
//! 5. **decomposition-none-for-non-decomposable**: `decompose(Q)` returns `None` for a
//!    point `Q` that is not a sum of any two factor-base points (verified by exhaustive
//!    enumeration of all factor-base pairs).
//! 6. **decomposition-sum-check**: for any returned decomposition `[P_i, P_j]`, verify
//!    `P_i + P_j = Q` using the group law (the primary correctness signal).
//!
//! # KAT coverage — relation collection
//!
//! 7. **relation-validity**: for every relation returned by `collect_relations`, verify
//!    that `Σ e_i·P_i = a·G + b·Q` via the frozen group law. Guards the provenance
//!    `(a, b)` and the exponent-vector encoding.
//! 8. **over-determination**: the returned collection has at least `fb_size + 1` relations.
//!    Guards the loop termination condition.
//!
//! # KAT coverage — Z/ℓZ linear algebra
//!
//! 9. **adapter_fidelity**: `build_ek_matrix` produces the transposed layout: row `j`
//!    (one per factor-base point) lists `(relation index, exponent)` for exactly the
//!    relations that contain point `j`, in relation order. Guards the adapter from
//!    `Relation.exponents` to `FlSparseRow`.
//! 10. **left_kernel_annihilates_exponents**: every vector `u` in the basis returned by
//!    `solve_ek_linalg` satisfies `Σᵢ uᵢ·eᵢ ≡ 0 (mod ℓ)`, is non-zero, and the basis size
//!    equals `#relations − rank` (a complete basis). Guards the Z/ℓZ linear algebra step.
//!
//! # KAT coverage — DLP recovery and end-to-end pipeline
//!
//! 11. **ecdlp_agreement** (decisive): `index_calculus_dlp(G, Q)` returns `x = log_G(Q) mod ℓ`
//!    that agrees with the `rho::ecdlp::pohlig::solve_ecdlp_composite` oracle: `k % ℓ == x`.
//!    Also verifies `x·G_ℓ = Q_ℓ` (the ℓ-subgroup projection check).
//! 12. **recovery_soundness**: `x·G_ℓ = Q_ℓ` via `scalar_mul`, independently of the rho
//!    cross-check.
//! 13. **end_to_end**: the full pipeline solves a non-trivial toy ECDLP instance end-to-end
//!    (factor base → decomposition → collection → linalg → recovery → verification).
//! 14. **recovers_every_residue**: for every `k ∈ 0..60`, `Q = k·G` gives `x = k mod ℓ` with
//!    `uᵀb ≢ 0` and `Σᵢ uᵢ·(aᵢ + bᵢ·x) ≡ 0` checked explicitly from the returned provenance.
//! 15. **undetermined_when_no_combination_sees_b**: honest relations that all have
//!    `b ≡ 0 (mod ℓ)` give `Undetermined`.
//! 16. **corrupted_relation_is_inconsistent**: one corrupted relation gives `Inconsistent`.
//! 17. **toy_fixture_parameters**: the fixture constants that bound what is demonstrated.
//!
//! # Scale: index calculus over `E(F_p)` at toy scale
//!
//! The fixture is toy-scale (`p = 47`, `n = 60`). Over `E(F_p)` these tests demonstrate
//! the mechanism, not the asymptotic advantage: that advantage needs `E(F_{p^n})` and is not
//! observable at this scale.

use crypto_bigint::Uint;
use rho::curve::{AffinePoint, JacobianPoint};
use rho::ecdlp::pohlig::solve_ecdlp_composite;
use shared_field::{Fp, FpNaive4 as FpNaive};
use rho::index_calculus::solve::recover_log_from_relations;
use rho::index_calculus::{
    build_ek_matrix, collect_relations, decompose, index_calculus_dlp, solve_ek_linalg,
    IndexCalcError, IndexCalcStrategy, Relation, RelationStream, TOY_ELL, TOY_FB_SIZE,
};

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Return the `Uint<4>` prime-order subgroup modulus for the toy fixture.
fn ell() -> Uint<4> {
    Uint::<4>::from(TOY_ELL)
}

/// Add two affine points using the frozen `Curve` group law.
///
/// Converts to Jacobian, adds, converts back to affine.
fn add_points(
    curve: &rho::curve::Curve,
    p1: &AffinePoint<FpNaive>,
    p2: &AffinePoint<FpNaive>,
) -> AffinePoint<FpNaive> {
    let p = &curve.p;
    let j1 = JacobianPoint::from_affine(p1, p);
    let j2 = JacobianPoint::from_affine(p2, p);
    curve.add_jacobian(&j1, &j2).to_affine(p)
}

/// Scalar-multiply an affine point by a `u64` scalar.
fn scalar_mul_u64(
    curve: &rho::curve::Curve,
    pt: &AffinePoint<FpNaive>,
    k: u64,
) -> AffinePoint<FpNaive> {
    curve.scalar_mul(pt, &Uint::<4>::from(k))
}

// ─── KAT 1: factor-base-on-curve ─────────────────────────────────────────────

/// Every enumerated factor-base point satisfies `is_on_curve` via the frozen `Curve`.
///
/// Guards the QR-test + canonical-root-lift enumeration in `enumerate_factor_base`.
/// If any point is off-curve, the enumeration has a bug (wrong y, wrong x, or wrong
/// QR test).
#[test]
fn factor_base_on_curve() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    assert_eq!(
        strategy.fb_size(),
        TOY_FB_SIZE,
        "factor base should have exactly TOY_FB_SIZE = {TOY_FB_SIZE} points"
    );

    for fb in &strategy.factor_base {
        assert!(
            strategy.curve.is_on_curve(&fb.point),
            "factor-base point {} is not on the curve: {:?}",
            fb.index,
            fb.point
        );
        // Also verify it is a finite point (not the point at infinity).
        assert!(
            !fb.point.is_infinity(),
            "factor-base point {} is the point at infinity (unexpected)",
            fb.index
        );
    }
}

// ─── KAT 2: subgroup-validity ─────────────────────────────────────────────────

/// `ℓ | n` and `ℓ·G_ℓ = ∞` (the prime-order-subgroup precondition).
///
/// Guards the `Z/ℓℤ` linear-algebra precondition: the relation exponents live in a
/// field `F_ℓ` only if ℓ is prime and divides n. The KAT checks both the divisibility
/// condition and the order of the subgroup generator `G_ℓ = (n/ℓ)·G`.
#[test]
fn subgroup_validity() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    // ℓ | n.
    let n_u64 = curve.n.as_words()[0];
    let ell_u64 = strategy.ell.as_words()[0];
    assert_eq!(
        n_u64 % ell_u64,
        0,
        "ℓ = {ell_u64} must divide n = {n_u64}"
    );

    // G_ℓ = (n/ℓ)·G is not the identity.
    let g_ell = strategy.subgroup_generator();
    assert!(
        !g_ell.is_infinity(),
        "G_ℓ = (n/ℓ)·G should not be ∞ (it has order ℓ = {ell_u64})"
    );

    // ℓ·G_ℓ = ∞.
    let ell_times_g_ell = curve.scalar_mul(&g_ell, &strategy.ell);
    assert!(
        ell_times_g_ell.is_infinity(),
        "ℓ·G_ℓ should be ∞ (G_ℓ has order ℓ = {ell_u64}), got {:?}",
        ell_times_g_ell
    );

    // Sanity: (ℓ-1)·G_ℓ ≠ ∞ (G_ℓ has exact order ℓ, not a proper divisor).
    let ell_minus_1 = Uint::<4>::from(ell_u64 - 1);
    let ell_minus_1_times_g_ell = curve.scalar_mul(&g_ell, &ell_minus_1);
    assert!(
        !ell_minus_1_times_g_ell.is_infinity(),
        "(ℓ-1)·G_ℓ should not be ∞ (G_ℓ has exact order ℓ = {ell_u64})"
    );
}

// ─── KAT 3: Relation round-trip ──────────────────────────────────────────────

/// Relation round-trip: `Σ e_i·P_i` reconstructed from the relation equals the
/// expected point sum.
///
/// Constructs a `Relation` from a known decomposition (factor-base points 0 and 1,
/// each appearing once), then reconstructs the sum `e_0·P_0 + e_1·P_1` using the
/// frozen group law and verifies it equals `P_0 + P_1` (the direct sum).
///
/// This guards:
/// - `Relation::from_decomposition` correctly encodes the decomposition as a sparse
///   exponent vector.
/// - `Relation::exponent` correctly retrieves exponents (including zero for absent indices).
/// - The exponent vector faithfully represents the decomposition (the round-trip invariant
///   the whole pipeline relies on: a relation's exponent vector reconstructs its recorded
///   factor-base-point sum).
#[test]
fn relation_round_trip() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let ell = ell();

    // Take the first two factor-base points.
    let p0 = &strategy.factor_base[0].point;
    let p1 = &strategy.factor_base[1].point;

    // Compute the expected sum P_0 + P_1 directly via the frozen group law.
    let expected_sum = add_points(curve, p0, p1);

    // Construct a Relation encoding the decomposition [P_0, P_1] (each once).
    // Provenance (a=0, b=0) is placeholder — the round-trip KAT checks the exponent
    // vector, not the provenance.
    let relation = Relation::from_decomposition(0, 0, &[0, 1], &ell);

    // Verify the exponent vector: e_0 = 1, e_1 = 1, all others = 0.
    assert_eq!(
        relation.exponent(0, &ell),
        FpNaive::from_u64(1, &ell),
        "exponent of factor-base point 0 should be 1"
    );
    assert_eq!(
        relation.exponent(1, &ell),
        FpNaive::from_u64(1, &ell),
        "exponent of factor-base point 1 should be 1"
    );
    assert_eq!(
        relation.exponent(2, &ell),
        FpNaive::from_u64(0, &ell),
        "exponent of absent factor-base point 2 should be 0"
    );

    // Reconstruct the sum Σ e_i·P_i from the relation's exponent vector.
    // For each (i, e_i) in the relation, add e_i·P_i to the running sum.
    let p = &curve.p;
    let mut reconstructed = JacobianPoint::<FpNaive>::infinity(p);
    for (i, exp) in &relation.exponents {
        let fb_point = &strategy.factor_base[*i].point;
        // e_i is in F_ℓ; convert to u64 for scalar_mul.
        let exp_u64 = exp.to_uint().as_words()[0];
        if exp_u64 == 0 {
            continue;
        }
        let contribution = scalar_mul_u64(curve, fb_point, exp_u64);
        let contrib_jac = JacobianPoint::from_affine(&contribution, p);
        reconstructed = curve.add_jacobian(&reconstructed, &contrib_jac);
    }
    let reconstructed_affine = reconstructed.to_affine(p);

    assert_eq!(
        reconstructed_affine, expected_sum,
        "round-trip failed: Σ e_i·P_i = {:?}, expected P_0 + P_1 = {:?}",
        reconstructed_affine, expected_sum
    );
}

// ─── KAT 3b: Relation round-trip with repeated indices ───────────────────────

/// Relation round-trip with a repeated factor-base index.
///
/// Constructs a `Relation` from the decomposition [P_0, P_0] (factor-base point 0
/// appearing twice), then verifies that `Σ e_i·P_i = 2·P_0`.
///
/// Guards the accumulation logic in `from_decomposition` (repeated indices → exponent 2).
#[test]
fn relation_round_trip_repeated() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let ell = ell();

    let p0 = &strategy.factor_base[0].point;

    // Expected: 2·P_0.
    let expected_sum = scalar_mul_u64(curve, p0, 2);

    // Decomposition: [P_0, P_0] → exponent vector [(0, 2)].
    let relation = Relation::from_decomposition(0, 0, &[0, 0], &ell);

    assert_eq!(
        relation.exponent(0, &ell),
        FpNaive::from_u64(2, &ell),
        "exponent of factor-base point 0 (appearing twice) should be 2"
    );
    assert_eq!(relation.exponents.len(), 1, "only one non-zero exponent");

    // Reconstruct Σ e_i·P_i.
    let p = &curve.p;
    let mut reconstructed = JacobianPoint::<FpNaive>::infinity(p);
    for (i, exp) in &relation.exponents {
        let fb_point = &strategy.factor_base[*i].point;
        let exp_u64 = exp.to_uint().as_words()[0];
        if exp_u64 == 0 {
            continue;
        }
        let contribution = scalar_mul_u64(curve, fb_point, exp_u64);
        let contrib_jac = JacobianPoint::from_affine(&contribution, p);
        reconstructed = curve.add_jacobian(&reconstructed, &contrib_jac);
    }
    let reconstructed_affine = reconstructed.to_affine(p);

    assert_eq!(
        reconstructed_affine, expected_sum,
        "round-trip (repeated) failed: 2·P_0 via relation = {:?}, direct = {:?}",
        reconstructed_affine, expected_sum
    );
}

// ─── KAT 4: decomposition-correctness ───────────────────────────────────────

/// Decomposition correctness: `decompose(Q)` returns `Some(pts)` for a point `Q` that
/// is a sum of two factor-base points; the returned decomposition sums to `Q` via the
/// frozen group law.
///
/// Constructs `Q = P_0 + P_1` (the sum of the first two factor-base points) directly
/// via the frozen group law, then calls `decompose(Q)` and verifies:
/// - The result is `Some(...)` (the point is decomposable).
/// - The returned decomposition has exactly `m = 2` points.
/// - The sum of the returned points equals `Q` via the frozen group law.
///
/// Guards the Semaev root-finding step (the primary correctness signal for point decomposition).
#[test]
fn decomposition_correctness() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    // Construct Q = P_0 + P_1 directly via the frozen group law.
    // This Q is guaranteed to be decomposable (it is a sum of two factor-base points).
    let p0 = &strategy.factor_base[0].point;
    let p1 = &strategy.factor_base[1].point;
    let q = add_points(curve, p0, p1);

    assert!(
        !q.is_infinity(),
        "P_0 + P_1 should not be ∞ for the toy fixture (they are distinct non-negation points)"
    );

    // Call decompose(Q) — must return Some(...).
    let decomp = decompose(q.clone(), &strategy)
        .expect("decompose(P_0 + P_1) should return Some(...): Q is a sum of two FB points");

    // The decomposition must have exactly m = 2 points.
    assert_eq!(
        decomp.len(),
        strategy.m,
        "decomposition should have exactly m = {} points",
        strategy.m
    );

    // All returned points must be factor-base points (valid indices).
    for fb_pt in &decomp {
        assert!(
            fb_pt.index < strategy.fb_size(),
            "decomposition point index {} is out of factor-base range [0, {})",
            fb_pt.index,
            strategy.fb_size()
        );
        assert!(
            curve.is_on_curve(&fb_pt.point),
            "decomposition point {} is not on the curve: {:?}",
            fb_pt.index,
            fb_pt.point
        );
    }

    // The sum of the returned points must equal Q via the frozen group law.
    let sum = add_points(curve, &decomp[0].point, &decomp[1].point);
    assert_eq!(
        sum, q,
        "decompose returned a decomposition that does not sum to Q: \
         decomp[0] = {:?}, decomp[1] = {:?}, sum = {:?}, Q = {:?}",
        decomp[0].point, decomp[1].point, sum, q
    );
}

// ─── KAT 5: decomposition-none-for-non-decomposable ─────────────────────────

/// Decomposition returns `None` for a point that is not a sum of any two factor-base
/// points.
///
/// Finds a point `Q` (a multiple of `G`) that is not the sum of any two factor-base
/// points by exhaustive enumeration of all factor-base pairs. Then calls `decompose(Q)`
/// and verifies it returns `None`.
///
/// Guards the `None` path of `decompose` — the Semaev polynomial correctly identifies
/// non-decomposable points.
#[test]
fn decomposition_none_for_non_decomposable() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    // Build the set of all points that ARE sums of two factor-base points.
    // This is the set { P_i + P_j : 0 ≤ i, j < FB_SIZE, P_i + P_j ≠ ∞ }.
    let mut decomposable_points: Vec<AffinePoint<FpNaive>> = Vec::new();
    for fb_i in &strategy.factor_base {
        for fb_j in &strategy.factor_base {
            let sum = add_points(curve, &fb_i.point, &fb_j.point);
            if !sum.is_infinity() {
                decomposable_points.push(sum);
            }
        }
    }

    // Find a multiple of G that is NOT in the decomposable set.
    // Try k·G for k = 1, 2, ..., n-1 until we find a non-decomposable point.
    let g: AffinePoint<FpNaive> = curve.generator();
    let n_u64 = curve.n.as_words()[0]; // n = 60

    let mut non_decomposable_q: Option<AffinePoint<FpNaive>> = None;
    for k in 1..n_u64 {
        let q = scalar_mul_u64(curve, &g, k);
        if q.is_infinity() {
            continue;
        }
        if !decomposable_points.contains(&q) {
            non_decomposable_q = Some(q);
            break;
        }
    }

    let q = non_decomposable_q
        .expect("should find a non-decomposable multiple of G in the toy fixture");

    // decompose(Q) must return None.
    let result = decompose(q.clone(), &strategy);
    assert!(
        result.is_none(),
        "decompose(Q) should return None for a non-decomposable point Q = {:?}, \
         but returned {:?}",
        q,
        result
    );
}

// ─── KAT 6: decomposition-sum-check ─────────────────────────────────────────

/// For every returned decomposition `[P_i, P_j]`, verify `P_i + P_j = Q` using the
/// frozen group law.
///
/// Exhaustively tries all factor-base pairs as candidate Q values, calls `decompose(Q)`
/// for each, and verifies that every returned decomposition sums to Q. This is the
/// primary correctness signal: the Semaev polynomial correctly identifies decomposable
/// points and the returned decomposition is valid.
///
/// Guards the sum-to-Q invariant across all decomposable points in the factor-base
/// span (not just one example).
#[test]
fn decomposition_sum_check() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let mut checked = 0usize;

    // Try all sums of two factor-base points as candidate Q values.
    for fb_i in &strategy.factor_base {
        for fb_j in &strategy.factor_base {
            let q = add_points(curve, &fb_i.point, &fb_j.point);
            if q.is_infinity() {
                continue;
            }

            // Call decompose(Q).
            if let Some(decomp) = decompose(q.clone(), &strategy) {
                assert_eq!(
                    decomp.len(),
                    strategy.m,
                    "decomposition should have exactly m = {} points",
                    strategy.m
                );

                // Verify P_i + P_j = Q via the frozen group law.
                let sum = add_points(curve, &decomp[0].point, &decomp[1].point);
                assert_eq!(
                    sum, q,
                    "decomposition sum check failed: decomp[0] = {:?}, decomp[1] = {:?}, \
                     sum = {:?}, Q = {:?}",
                    decomp[0].point, decomp[1].point, sum, q
                );

                // Verify all returned points are on the curve.
                for fb_pt in &decomp {
                    assert!(
                        curve.is_on_curve(&fb_pt.point),
                        "decomposition point {} is not on the curve: {:?}",
                        fb_pt.index,
                        fb_pt.point
                    );
                }

                checked += 1;
            }
        }
    }

    assert!(
        checked > 0,
        "no decomposable Q found among all factor-base pairs — \
         unexpected for the toy fixture (factor base should span some decomposable points)"
    );
}

// ─── KAT 7: relation-validity ────────────────────────────────────────────────

/// For every relation returned by `collect_relations`, verify `Σ e_i·P_i = a·G + b·Q`.
///
/// Calls `collect_relations(G, Q, &strategy)` for the toy fixture. For each returned
/// `Relation`, reconstructs `Σ e_i·P_i` from the exponent vector using the frozen group
/// law, and verifies it equals `a·G + b·Q` (computed from the provenance fields).
///
/// Guards:
/// - The provenance `(a, b)` is faithfully recorded.
/// - The exponent vector correctly encodes the decomposition.
/// - The `collect_relations` loop only records valid decompositions.
#[test]
fn relation_validity() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let p = &curve.p;

    // Use G and Q = 7·G as the toy fixture points.
    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let relations = collect_relations(g.clone(), q.clone(), &strategy)
        .expect("collect_relations should succeed for the toy fixture");

    assert!(
        !relations.is_empty(),
        "collect_relations returned an empty relation set"
    );

    for (idx, relation) in relations.iter().enumerate() {
        // Compute a·G + b·Q from the provenance fields.
        let ag = scalar_mul_u64(curve, &g, relation.a);
        let bq = scalar_mul_u64(curve, &q, relation.b);
        let expected = if bq.is_infinity() {
            ag.clone()
        } else if ag.is_infinity() {
            bq.clone()
        } else {
            add_points(curve, &ag, &bq)
        };

        // Reconstruct Σ e_i·P_i from the relation's exponent vector.
        let mut reconstructed = JacobianPoint::<FpNaive>::infinity(p);
        for (i, exp) in &relation.exponents {
            let fb_point = &strategy.factor_base[*i].point;
            let exp_u64 = exp.to_uint().as_words()[0];
            if exp_u64 == 0 {
                continue;
            }
            let contribution = scalar_mul_u64(curve, fb_point, exp_u64);
            let contrib_jac = JacobianPoint::from_affine(&contribution, p);
            reconstructed = curve.add_jacobian(&reconstructed, &contrib_jac);
        }
        let reconstructed_affine = reconstructed.to_affine(p);

        assert_eq!(
            reconstructed_affine, expected,
            "relation {} (a={}, b={}): Σ e_i·P_i = {:?}, expected a·G + b·Q = {:?}",
            idx, relation.a, relation.b, reconstructed_affine, expected
        );
    }
}

// ─── KAT 8: over-determination ───────────────────────────────────────────────

/// The returned collection has at least `fb_size + 1` relations.
///
/// Calls `collect_relations(G, Q, &strategy)` for the toy fixture and verifies that
/// the returned `Vec<Relation>` has at least `strategy.fb_size() + 1` entries. This
/// is the over-determination condition required for the index-calculus linear algebra
/// the linear-algebra step to have a non-trivial kernel.
///
/// Guards the loop termination condition in `collect_relations`.
#[test]
fn over_determination() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let relations = collect_relations(g, q, &strategy)
        .expect("collect_relations should succeed for the toy fixture");

    let min_required = strategy.fb_size() + 1;
    assert!(
        relations.len() >= min_required,
        "over-determination check failed: got {} relations, need at least {} (fb_size + 1 = {} + 1)",
        relations.len(),
        min_required,
        strategy.fb_size()
    );
}

// ─── KAT 9: adapter-fidelity ─────────────────────────────────────────────────

/// Adapter fidelity for the transposed layout.
///
/// Calls `collect_relations` to get a real relation set, then `build_ek_matrix`. The matrix
/// has one row per factor-base point and one column per relation, and row `j` holds exactly
/// the `(relation index i, exponent)` pairs for which `relations[i]` contains point `j`,
/// sorted by `i`. Every non-zero relation exponent appears in exactly one matrix entry.
#[test]
fn adapter_fidelity() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let relations = collect_relations(g, q, &strategy)
        .expect("collect_relations should succeed for the toy fixture");

    let matrix = build_ek_matrix(&relations, &strategy);

    // Verify dimensions: factor-base points are rows, relations are columns.
    assert_eq!(
        matrix.rows.len(),
        strategy.fb_size(),
        "matrix should have one row per factor-base point"
    );
    assert_eq!(
        matrix.num_cols,
        relations.len(),
        "matrix should have one column per relation"
    );

    // Row j lists, in relation order, the (i, e_ij) of every relation i containing point j.
    let ell = ell();
    for (j, row) in matrix.rows.iter().enumerate() {
        let expected: Vec<(usize, FpNaive)> = relations
            .iter()
            .enumerate()
            .filter_map(|(i, rel)| {
                let e = rel.exponent(j, &ell);
                if e.is_zero(&ell) { None } else { Some((i, e)) }
            })
            .collect();
        assert_eq!(
            row.entries, expected,
            "row {j}: entries should be the exponents of point {j} across relations"
        );
    }

    // Nothing is lost or invented: entry count equals the total non-zero exponent count.
    let matrix_entries: usize = matrix.rows.iter().map(|r| r.entries.len()).sum();
    let exponent_entries: usize = relations.iter().map(|r| r.exponents.len()).sum();
    assert_eq!(matrix_entries, exponent_entries, "total entry count must be preserved");
}

// ─── KAT 10: left_kernel_annihilates_exponents ───────────────────────────────

/// Rank of an integer matrix over `F_ℓ` by plain row reduction (independent of the
/// production linear-algebra engine; ℓ is prime).
fn rank_mod_ell(mut rows: Vec<Vec<u64>>, ell: u64) -> usize {
    let inv = |a: u64| (1..ell).find(|x| a * x % ell == 1).expect("ℓ prime, a ≠ 0");
    let cols = rows.first().map_or(0, Vec::len);
    let mut rank = 0;
    for col in 0..cols {
        let Some(p) = (rank..rows.len()).find(|&r| rows[r][col] % ell != 0) else { continue };
        rows.swap(rank, p);
        let scale = inv(rows[rank][col] % ell);
        for c in 0..cols {
            rows[rank][c] = rows[rank][c] * scale % ell;
        }
        for r in 0..rows.len() {
            let f = rows[r][col] % ell;
            if r != rank && f != 0 {
                for c in 0..cols {
                    rows[r][c] = (rows[r][c] + ell * ell - f * rows[rank][c]) % ell;
                }
            }
        }
        rank += 1;
    }
    rank
}

/// Left-kernel correctness: every `u` returned by `solve_ek_linalg` satisfies
/// `Σᵢ uᵢ·eᵢ ≡ 0 (mod ℓ)`.
///
/// Calls `collect_relations` then `solve_ek_linalg`. For each basis vector `u` and each
/// factor-base point `j`, computes `Σᵢ uᵢ·e_ij` over F_ℓ and verifies it is zero.
///
/// Guards:
/// - `build_ek_matrix` encodes the relation system in the orientation the engine solves.
/// - Each returned vector has length `#relations`, is non-zero, and is a genuine
///   combination (not a spurious solution).
/// - The basis is complete: its size is `#relations − rank(E)`, with the rank computed by an
///   independent row reduction, and the vectors are linearly independent (so they span the
///   whole kernel).
#[test]
fn left_kernel_annihilates_exponents() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let ell = ell();
    let ell_u64 = TOY_ELL;

    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let relations = collect_relations(g, q, &strategy)
        .expect("collect_relations should succeed for the toy fixture");

    let basis = solve_ek_linalg(&relations, &strategy)
        .expect("solve_ek_linalg should succeed for the over-determined toy system");

    // Over-determined (#relations > |FB|) forces a non-trivial kernel.
    assert!(!basis.is_empty(), "an over-determined system has a non-trivial left kernel");

    for (k, u) in basis.iter().enumerate() {
        assert_eq!(u.len(), relations.len(), "basis vector {k} has one entry per relation");
        assert!(u.iter().any(|c| !c.is_zero(&ell)), "basis vector {k} must be non-zero");

        // For every factor-base point j: Σᵢ uᵢ·e_ij ≡ 0 (mod ℓ).
        for j in 0..strategy.fb_size() {
            let mut sum = FpNaive::zero(&ell);
            for (ui, rel) in u.iter().zip(&relations) {
                sum = sum.add(&ui.mul(&rel.exponent(j, &ell), &ell), &ell);
            }
            assert!(
                sum.is_zero(&ell),
                "basis vector {k}: Σ uᵢ·e_i{j} = {sum:?}, expected 0 mod ℓ"
            );
        }
    }

    // Completeness: dim(left kernel) = #relations − rank(E), vectors independent.
    let dense = |rows: &[Vec<FpNaive>]| -> Vec<Vec<u64>> {
        rows.iter()
            .map(|r| r.iter().map(|c| c.to_uint().as_words()[0]).collect())
            .collect()
    };
    let exponent_rows: Vec<Vec<u64>> = relations
        .iter()
        .map(|rel| {
            (0..strategy.fb_size())
                .map(|j| rel.exponent(j, &ell).to_uint().as_words()[0])
                .collect()
        })
        .collect();
    let rank_e = rank_mod_ell(exponent_rows, ell_u64);
    assert_eq!(
        basis.len(),
        relations.len() - rank_e,
        "kernel basis size must equal #relations − rank(E)"
    );
    assert_eq!(
        rank_mod_ell(dense(&basis), ell_u64),
        basis.len(),
        "basis vectors must be linearly independent"
    );
}

// ─── KAT 11: ecdlp_agreement (decisive cross-check) ─────────────────────────

/// Decisive cross-check: `index_calculus_dlp(G, Q)` agrees with `rho::ecdlp` on the
/// same toy instance.
///
/// The index-calculus recovered `x = log_G(Q) mod ℓ` must satisfy `k % ℓ == x` where
/// `k` is the full discrete log returned by the
/// `rho::ecdlp::pohlig::solve_ecdlp_composite` oracle. Also verifies the ℓ-subgroup
/// projection: `x·G_ℓ = Q_ℓ`.
///
/// # Scale
///
/// Over `E(F_p)` at toy scale this demonstrates the index-calculus mechanism, not its
/// asymptotic advantage over Pollard rho, which needs `E(F_{p^n})`.
///
/// # Cross-check design
///
/// `index_calculus_dlp` returns `x = log_G(Q) mod ℓ` (mod 5 for the toy).
/// `solve_ecdlp_composite` returns `k` such that `k·G = Q` (the full discrete log mod n).
/// The cross-check: `k % ℓ == x` — the index-calculus result is the correct residue of
/// the full discrete log modulo the prime-order subgroup modulus ℓ.
///
/// Additionally: `x·G_ℓ = Q_ℓ` where `G_ℓ = (n/ℓ)·G` and `Q_ℓ = (n/ℓ)·Q` (the
/// ℓ-order subgroup projections). This is the ℓ-subgroup projection check — the
/// index-calculus result is the correct discrete log in the ℓ-order subgroup.
#[test]
fn ecdlp_agreement() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    // Use Q = 7·G as the non-trivial target (not in the factor base).
    // log_G(Q) = 7; log_G(Q) mod ℓ = 7 mod 5 = 2.
    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    // ── Step 1: Index-calculus DLP recovery ──────────────────────────────────
    let x = index_calculus_dlp(g.clone(), q.clone(), &strategy)
        .expect("index_calculus_dlp should recover log_G(Q) mod ℓ for Q = 7·G")
        .x;

    // x must be in [0, ℓ).
    let ell_u64 = strategy.ell.as_words()[0];
    assert!(
        x < ell_u64,
        "recovered x = {x} is not in [0, ℓ = {ell_u64})"
    );

    // ── Step 2: rho::ecdlp cross-check ───────────────────────────────────────
    // solve_ecdlp_composite returns k such that k·G = Q (full discrete log mod n).
    let n_u64 = curve.n.as_words()[0];
    let k = solve_ecdlp_composite(curve, &g, &q, n_u64)
        .expect("solve_ecdlp_composite should succeed for Q = 7·G on the toy fixture");

    // Verify the rho oracle is correct: k·G = Q.
    let k_g = scalar_mul_u64(curve, &g, k);
    assert_eq!(
        k_g, q,
        "rho oracle sanity: k·G ≠ Q (k = {k})"
    );

    // Cross-check: k % ℓ == x (the index-calculus residue matches the full log mod ℓ).
    assert_eq!(
        k % ell_u64,
        x,
        "ecdlp_agreement FAILED: index-calculus x = {x}, rho k = {k}, k % ℓ = {} ≠ x",
        k % ell_u64
    );

    // ── Step 3: ℓ-subgroup projection check ──────────────────────────────────
    // G_ℓ = (n/ℓ)·G, Q_ℓ = (n/ℓ)·Q. Verify x·G_ℓ = Q_ℓ.
    let cofactor = n_u64 / ell_u64;
    let g_ell = scalar_mul_u64(curve, &g, cofactor);
    let q_ell = scalar_mul_u64(curve, &q, cofactor);

    // x·G_ℓ should equal Q_ℓ.
    let x_g_ell = scalar_mul_u64(curve, &g_ell, x);
    assert_eq!(
        x_g_ell, q_ell,
        "ℓ-subgroup projection check FAILED: x·G_ℓ = {:?}, Q_ℓ = {:?} (x = {x})",
        x_g_ell, q_ell
    );
}

// ─── KAT 12: recovery_soundness ──────────────────────────────────────────────

/// Recovery soundness: `x·G_ℓ = Q_ℓ` via the frozen `scalar_mul`.
///
/// Verifies the ℓ-subgroup projection check independently of the rho cross-check:
/// the index-calculus recovered `x = log_G(Q) mod ℓ` satisfies `x·G_ℓ = Q_ℓ` where
/// `G_ℓ = (n/ℓ)·G` and `Q_ℓ = (n/ℓ)·Q` are the ℓ-order subgroup projections.
///
/// This guards the DLP recovery step (the relation combination) independently of the rho
/// oracle. The check is: the recovered `x` is the correct discrete log in the ℓ-order
/// subgroup, verified by the group law.
#[test]
fn recovery_soundness() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let x = index_calculus_dlp(g.clone(), q.clone(), &strategy)
        .expect("index_calculus_dlp should recover log_G(Q) mod ℓ")
        .x;

    let n_u64 = curve.n.as_words()[0];
    let ell_u64 = strategy.ell.as_words()[0];
    let cofactor = n_u64 / ell_u64;

    // Compute G_ℓ = (n/ℓ)·G and Q_ℓ = (n/ℓ)·Q.
    let g_ell = scalar_mul_u64(curve, &g, cofactor);
    let q_ell = scalar_mul_u64(curve, &q, cofactor);

    // G_ℓ must not be ∞ (it has order ℓ = 5).
    assert!(
        !g_ell.is_infinity(),
        "G_ℓ = (n/ℓ)·G should not be ∞ (it has order ℓ = {ell_u64})"
    );

    // Q_ℓ must not be ∞ (since log_G(Q) mod ℓ ≠ 0 for Q = 7·G, 7 mod 5 = 2 ≠ 0).
    assert!(
        !q_ell.is_infinity(),
        "Q_ℓ = (n/ℓ)·Q should not be ∞ for Q = 7·G (7 mod 5 = 2 ≠ 0)"
    );

    // x·G_ℓ must equal Q_ℓ.
    let x_g_ell = scalar_mul_u64(curve, &g_ell, x);
    assert_eq!(
        x_g_ell, q_ell,
        "recovery_soundness FAILED: x·G_ℓ = {:?}, Q_ℓ = {:?} (x = {x})",
        x_g_ell, q_ell
    );
}

// ─── KAT 13: end_to_end ──────────────────────────────────────────────────────

/// End-to-end: the full index-calculus pipeline solves a non-trivial toy ECDLP.
///
/// Runs the complete pipeline (factor base → decomposition → collection → linalg →
/// recovery) on the toy fixture for Q = 7·G and verifies:
/// - The pipeline returns `Ok(solution)` (no error).
/// - `x` is in `[0, ℓ)` (a valid F_ℓ element).
/// - `x·G_ℓ = Q_ℓ` (the ℓ-subgroup projection check — the primary soundness signal).
#[test]
fn end_to_end() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    // Run the full pipeline.
    let result = index_calculus_dlp(g.clone(), q.clone(), &strategy);
    assert!(
        result.is_ok(),
        "end_to_end: index_calculus_dlp returned an error: {:?}",
        result
    );

    let x = result.unwrap().x;

    let ell_u64 = strategy.ell.as_words()[0];
    assert!(
        x < ell_u64,
        "end_to_end: recovered x = {x} is not in [0, ℓ = {ell_u64})"
    );

    // ℓ-subgroup projection check: x·G_ℓ = Q_ℓ.
    let n_u64 = curve.n.as_words()[0];
    let cofactor = n_u64 / ell_u64;
    let g_ell = scalar_mul_u64(curve, &g, cofactor);
    let q_ell = scalar_mul_u64(curve, &q, cofactor);
    let x_g_ell = scalar_mul_u64(curve, &g_ell, x);

    assert_eq!(
        x_g_ell, q_ell,
        "end_to_end: ℓ-subgroup check FAILED: x·G_ℓ = {:?}, Q_ℓ = {:?} (x = {x})",
        x_g_ell, q_ell
    );
}

// ─── KAT 14: recovers_every_residue ──────────────────────────────────────────

/// Evaluate `uᵀa`, `uᵀb` over F_ℓ from a solution's provenance.
fn combination_sums(
    solution: &rho::index_calculus::IndexCalcSolution,
    ell: &Uint<4>,
) -> (FpNaive, FpNaive) {
    let mut ua = FpNaive::zero(ell);
    let mut ub = FpNaive::zero(ell);
    for (u, rel) in solution.combination.iter().zip(&solution.relations) {
        ua = ua.add(&u.mul(&FpNaive::from_u64(rel.a, ell), ell), ell);
        ub = ub.add(&u.mul(&FpNaive::from_u64(rel.b, ell), ell), ell);
    }
    (ua, ub)
}

/// For every `k ∈ 0..60`, `Q = k·G` recovers `x = k mod ℓ` through a combination with
/// `uᵀb ≢ 0`, and the identity `Σᵢ uᵢ·(aᵢ + bᵢ·x) ≡ 0 (mod ℓ)` holds.
///
/// The identity is evaluated explicitly from the returned provenance (`combination` and
/// `relations`), and `u` is checked to annihilate the exponent rows, so the test rechecks
/// the mechanism rather than just the final answer. The cost is small (a toy group of
/// order 60), so it stays on the default path.
#[test]
fn recovers_every_residue() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let ell = ell();
    let g: AffinePoint<FpNaive> = curve.generator();
    let n = curve.n.as_words()[0];
    assert_eq!(n, 60);

    for k in 0..n {
        let q = scalar_mul_u64(curve, &g, k);
        let solution = index_calculus_dlp(g.clone(), q, &strategy)
            .unwrap_or_else(|e| panic!("k = {k}: index_calculus_dlp failed: {e}"));

        assert_eq!(solution.x, k % TOY_ELL, "k = {k}: x should be k mod ℓ");
        assert_eq!(solution.combination.len(), solution.relations.len(), "k = {k}");

        let (ua, ub) = combination_sums(&solution, &ell);
        assert!(!ub.is_zero(&ell), "k = {k}: the combination used must have uᵀb ≢ 0");

        // Σᵢ uᵢ·(aᵢ + bᵢ·x) = uᵀa + x·uᵀb ≡ 0, evaluated term by term.
        let x = FpNaive::from_u64(solution.x, &ell);
        let mut identity = FpNaive::zero(&ell);
        for (u, rel) in solution.combination.iter().zip(&solution.relations) {
            let a = FpNaive::from_u64(rel.a, &ell);
            let b = FpNaive::from_u64(rel.b, &ell);
            identity = identity.add(&u.mul(&a.add(&b.mul(&x, &ell), &ell), &ell), &ell);
        }
        assert!(identity.is_zero(&ell), "k = {k}: Σ uᵢ(aᵢ + bᵢx) should be 0 mod ℓ");
        assert_eq!(identity, ua.add(&x.mul(&ub, &ell), &ell), "k = {k}: identity expansion");

        // The combination cancels every factor-base logarithm: Σᵢ uᵢ·e_ij ≡ 0 for each j.
        for j in 0..strategy.fb_size() {
            let mut sum = FpNaive::zero(&ell);
            for (u, rel) in solution.combination.iter().zip(&solution.relations) {
                sum = sum.add(&u.mul(&rel.exponent(j, &ell), &ell), &ell);
            }
            assert!(sum.is_zero(&ell), "k = {k}: combination must annihilate column {j}");
        }
    }
}

// ─── KAT 15: undetermined_when_no_combination_sees_b ─────────────────────────

/// Honest relations that all have `b ≡ 0 (mod ℓ)` give `Undetermined`.
///
/// Every `u` has `uᵀb ≡ 0` for such a set, so no combination can reveal `x`. The relations
/// are real (drawn from the enumeration, kept only when `b ≡ 0 mod ℓ`) and there are more
/// than `|FB|` of them, so the left kernel is non-trivial and the outcome comes from the
/// `uᵀb` condition rather than from an empty kernel.
#[test]
fn undetermined_when_no_combination_sees_b() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    let mut stream = RelationStream::new(g.clone(), q.clone(), &strategy);
    let mut relations: Vec<Relation> = Vec::new();
    while relations.len() < TOY_FB_SIZE + 2 {
        let rel = stream.next_relation().expect("enough relations with b ≡ 0 mod ℓ");
        if rel.b % TOY_ELL == 0 {
            relations.push(rel);
        }
    }

    let basis = solve_ek_linalg(&relations, &strategy).unwrap();
    assert!(!basis.is_empty(), "more relations than factor-base points → non-trivial kernel");

    assert_eq!(
        recover_log_from_relations(&relations, &strategy, &g, &q),
        Err(IndexCalcError::Undetermined { relations: relations.len() })
    );
}

// ─── KAT 16: corrupted_relation_is_inconsistent ──────────────────────────────

/// One corrupted relation gives `Inconsistent`.
///
/// Starting from honest relations, picks a relation `i` that appears (`uᵢ ≠ 0`) in some
/// basis vector with `uᵀb ≢ 0`, and shifts its recorded `a` by one. For such a `u` the
/// recovered `x` moves by `−uᵢ/uᵀb ≠ 0`, so the result must disagree with another
/// combination or fail the `x·G_ℓ = Q_ℓ` check. It must never come back as a success or as
/// a different failure class.
#[test]
fn corrupted_relation_is_inconsistent() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;
    let ell = ell();
    let g: AffinePoint<FpNaive> = curve.generator();
    let q = scalar_mul_u64(curve, &g, 7);

    // Honest relations: a few batches so the kernel has several vectors.
    let mut stream = RelationStream::new(g.clone(), q.clone(), &strategy);
    let mut relations: Vec<Relation> = Vec::new();
    while relations.len() < 2 * (TOY_FB_SIZE + 1) {
        relations.push(stream.next_relation().expect("enough relations"));
    }
    let honest = recover_log_from_relations(&relations, &strategy, &g, &q)
        .expect("the honest relations recover x");
    assert_eq!(honest.x, 7 % TOY_ELL);

    // Find a relation that some x-revealing basis vector actually uses.
    let basis = solve_ek_linalg(&relations, &strategy).unwrap();
    let used = (0..relations.len())
        .find(|&i| {
            basis.iter().any(|u| {
                let ub = u.iter().zip(&relations).fold(FpNaive::zero(&ell), |acc, (c, r)| {
                    acc.add(&c.mul(&FpNaive::from_u64(r.b, &ell), &ell), &ell)
                });
                !ub.is_zero(&ell) && !u[i].is_zero(&ell)
            })
        })
        .expect("some relation is used by an x-revealing combination");

    let mut corrupted = relations.clone();
    corrupted[used].a += 1;
    assert_eq!(
        recover_log_from_relations(&corrupted, &strategy, &g, &q),
        Err(IndexCalcError::Inconsistent),
        "corrupting relation {used} must be detected"
    );
}

// ─── KAT 17: toy_fixture_parameters ──────────────────────────────────────────

/// The toy fixture constants that bound what this solver demonstrates.
///
/// Index calculus here recovers `log_G(Q) mod ℓ`, not mod `n`, on a fixture far too small
/// for the asymptotic comparison with Pollard rho (which needs `E(F_{p^n})`): the tests show
/// the mechanism, not the advantage.
///
/// - ℓ = 5 (prime-order subgroup modulus, the linear-algebra field), with ℓ ∥ n.
/// - n = 60 = 2²·3·5 (composite group order — the log is recovered mod ℓ, not mod n).
/// - m = 2 (decomposition arity, using S_3).
/// - p = 47.
#[test]
fn toy_fixture_parameters() {
    let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
    let curve = &strategy.curve;

    let ell_u64 = strategy.ell.as_words()[0];
    assert_eq!(ell_u64, 5, "toy ℓ should be 5 (the largest prime factor of n = 60)");

    let n_u64 = curve.n.as_words()[0];
    assert_eq!(n_u64, 60, "toy n should be 60 = 2²·3·5");
    assert_eq!(n_u64 % ell_u64, 0, "ℓ must divide n");
    assert_ne!((n_u64 / ell_u64) % ell_u64, 0, "ℓ divides n exactly once, so the ℓ-Sylow is cyclic");

    assert_eq!(strategy.m, 2, "toy m should be 2");

    let p_u64 = curve.p.as_words()[0];
    assert_eq!(p_u64, 47, "toy p should be 47 (the semaev_toy prime field)");

    // The collection target is fb_size + 1 relations.
    assert_eq!(
        strategy.fb_size(),
        TOY_FB_SIZE,
        "toy FB_SIZE should be {TOY_FB_SIZE}"
    );

    // The ℓ-order subgroup generator has order exactly ℓ.
    let g_ell = strategy.subgroup_generator();
    assert!(!g_ell.is_infinity(), "G_ℓ = (n/ℓ)·G should not be ∞");
    let ell_g_ell = curve.scalar_mul(&g_ell, &strategy.ell);
    assert!(
        ell_g_ell.is_infinity(),
        "ℓ·G_ℓ should be ∞ (G_ℓ has order ℓ = {ell_u64})"
    );
}
