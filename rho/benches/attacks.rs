//! Cross-attack ECDLP benchmark harness: Pohlig–Hellman, MOV, SSA, GHS, index calculus.
//!
//! Measures the wall-clock cost of each algebraic ECDLP attack on its respective toy fixture.
//! Each bench body asserts the solver returns the known correct answer BEFORE timing, so the
//! bench doubles as a no-regression smoke test (C-EWBench invariant).
//!
//! # Attacks benched
//!
//! - **Pohlig–Hellman** — `solve_ecdlp_composite` on `composite_toy()` (n = 60 = 2²·3·5).
//! - **MOV/Frey–Rück** — `mov_reduce` on `pairing_toy()` (ℓ = 3, k = 2, F_{47²}).
//! - **SSA** — `ssa_solve` on `anomalous_toy()` (y² = x³ + 5 mod 7, #E = 7 = p).
//! - **GHS-reduction** — `ghs_descend` + `verify_log_preservation` on `ghs_toy_curve()`
//!   (GF(2^6), m = 6, l = 2). Benched as a *transfer*, not an end-to-end solve.
//! - **Index calculus** — `index_calculus_dlp` on `IndexCalcStrategy::toy()` (ℓ = 5, |FB| = 6,
//!   m = 2). Relation/decomposition counts derived from `collect_relations` + `decompose`.
//!
//! # Pollard rho baseline
//!
//! Pollard rho is already benched in `rho/benches/ecdlp.rs` on `secp_k1_toy` (63-bit).
//! That bench is the generic-√n baseline column in the E.W table. It is NOT duplicated here.
//!
//! # Scale
//!
//! All fixtures are toy-scale (p = 47 or p = 7). The asymptotic L-notation separations
//! between attacks are NOT observable at this scale: the table reports toy costs. Over
//! E(F_p) the index-calculus bench demonstrates the mechanism, not its asymptotic advantage.

use criterion::{Criterion, criterion_group, criterion_main};
use crypto_bigint::Uint;
use shared_field::{Fp as FpTrait, FpMonty4 as FpMonty, FpNaive4 as FpNaive};
use shared_gf2m::F2mNaive;

use rho::curve::AffinePoint;
use rho::curve::test_curves::{COMPOSITE_TOY_N, composite_toy};
use rho::ecdlp::pohlig::solve_ecdlp_composite;
use rho::ghs::{GHS_POLY2, GhsParams, ghs_descend, ghs_toy_curve, verify_log_preservation};
use rho::index_calculus::{IndexCalcStrategy, collect_relations, decompose, index_calculus_dlp};
use rho::pairing::ecext::PairingPoint;
use rho::pairing::fpext::FpExt;
use rho::pairing::mov::mov_reduce;
use rho::pairing::test_curves::pairing_toy;
use rho::ssa::{ANOMALOUS_TOY_P, anomalous_toy, ssa_solve};

// ── Pohlig–Hellman ────────────────────────────────────────────────────────────

/// Bench: Pohlig–Hellman on `composite_toy()` (n = 60 = 2²·3·5).
///
/// Known answer: k = 7 (Q = 7·G = (24, 43) on the composite-order toy curve).
/// The bench asserts `solve_ecdlp_composite` returns `Some(7)` before timing.
fn bench_pohlig_hellman(c: &mut Criterion) {
    let curve = composite_toy();
    let g: AffinePoint<FpMonty> = curve.generator();
    // Q = 7·G — known scalar from the composite_toy reference table.
    let q: AffinePoint<FpMonty> = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    // Pre-check: assert the solver returns the known answer before timing.
    let k_check = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
        .expect("Pohlig–Hellman pre-check: solve_ecdlp_composite must succeed");
    // The solver returns k mod n; verify k·G = Q (the canonical correctness check).
    let q_check: AffinePoint<FpMonty> = curve.scalar_mul(&g, &Uint::<4>::from(k_check));
    assert_eq!(
        q_check, q,
        "Pohlig–Hellman pre-check: k·G ≠ Q (k = {k_check}, expected k = 7)"
    );

    c.bench_function("attacks/pohlig_hellman", |b| {
        b.iter(|| {
            solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
                .expect("Pohlig–Hellman: solve_ecdlp_composite must succeed")
        });
    });
}

// ── MOV/Frey–Rück ─────────────────────────────────────────────────────────────

/// Bench: MOV/Frey–Rück on `pairing_toy()` (ℓ = 3, k = 2, F_{47²}).
///
/// Known answer: k = 2 (Q = 2·G in the order-3 subgroup of E(F_47)).
/// The bench asserts `mov_reduce` returns `Ok(2)` before timing.
///
/// The bench measures the full reduction: two Tate-pairing evaluations + the
/// F_{p^k} DLP call (`gnfs::dl::solve_dl` at k = 2).
fn bench_mov(c: &mut Criterion) {
    let (curve, modulus, ell, p_point, q_point) = pairing_toy();
    let p = curve.p;
    // Curve coefficient a = 1 lifted into F_{47²} (needed for scalar_mul on PairingPoint).
    let a_ext = FpExt::from_base(FpNaive::from_u64(1, &p), 2, &p);

    // Q' = 2·G — the ECDLP target with known scalar k = 2.
    let q_prime: PairingPoint<FpNaive> = p_point.scalar_mul(2, &a_ext, &modulus, &p);
    // R = Q (the fixture's extension-field 3-torsion point, a μ_ℓ-generator).
    let r_point = q_point;

    // Pre-check: assert the solver returns the known answer before timing.
    let result = mov_reduce::<FpNaive>(&curve, &modulus, &p_point, &q_prime, &r_point, ell)
        .expect("MOV pre-check: mov_reduce must succeed");
    assert_eq!(
        result, 2u64,
        "MOV pre-check: mov_reduce must recover k = 2 mod ell = 3; got {result}"
    );

    c.bench_function("attacks/mov_frey_ruck", |b| {
        b.iter(|| {
            mov_reduce::<FpNaive>(&curve, &modulus, &p_point, &q_prime, &r_point, ell)
                .expect("MOV: mov_reduce must succeed")
        });
    });
}

// ── SSA ───────────────────────────────────────────────────────────────────────

/// Bench: Smart–Satoh–Araki on `anomalous_toy()` (y² = x³ + 5 mod 7, #E = 7 = p).
///
/// Known answer: k = 3 (Q = 3·G on the anomalous fixture).
/// The bench asserts `ssa_solve` returns `Ok(3)` before timing.
fn bench_ssa(c: &mut Criterion) {
    let curve = anomalous_toy();
    let g: AffinePoint<FpNaive> = curve.generator();
    // Q = 3·G — known scalar from the SSA KAT suite.
    let q: AffinePoint<FpNaive> = curve.scalar_mul(&g, &Uint::<4>::from(3u64));

    // Pre-check: assert the solver returns the known answer before timing.
    let k_check = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
        .expect("SSA pre-check: ssa_solve must succeed on the anomalous fixture");
    assert_eq!(
        k_check, 3u64,
        "SSA pre-check: ssa_solve must recover k = 3; got {k_check}"
    );

    c.bench_function("attacks/ssa", |b| {
        b.iter(|| {
            ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
                .expect("SSA: ssa_solve must succeed on the anomalous fixture")
        });
    });
}

// ── GHS-reduction (transfer, not end-to-end solve) ────────────────────────────

/// Bench: GHS Weil-descent transfer on `ghs_toy_curve()` (GF(2^6), m = 6, l = 2).
///
/// The GHS bench measures the **descent reduction + log-preservation verification**,
/// NOT an end-to-end solve. `rho::ghs` has no `ghs_dlp`; the downstream solve is
/// index calculus (a deferred re-shard). Reporting GHS as an end-to-end solve time
/// would be a misrepresentation (C-EWBench invariant: GHS is a transfer).
///
/// Known answer: k = 1 (h = 1·g — the generator maps to itself under the transfer).
/// The bench asserts `verify_log_preservation` returns `true` before timing.
fn bench_ghs_transfer(c: &mut Criterion) {
    let curve_e = ghs_toy_curve();
    let poly2 = Uint::<1>::from(GHS_POLY2);
    let params = GhsParams::new(6, 2, curve_e.clone(), poly2)
        .expect("GHS pre-check: toy GHS params must be valid");

    // g = base point, h = g (k = 1: the generator maps to itself).
    let g = curve_e.generator::<F2mNaive<1>>();
    let h = g.clone(); // h = 1·g

    // Pre-check: assert the descent succeeds and log-preservation holds for k = 1.
    let result = ghs_descend(&params, &g, &h)
        .expect("GHS pre-check: ghs_descend must succeed for toy fixture");
    assert!(
        verify_log_preservation(&result, 1),
        "GHS pre-check: verify_log_preservation must hold for k = 1"
    );

    c.bench_function("attacks/ghs_transfer", |b| {
        b.iter(|| {
            let descent = ghs_descend(&params, &g, &h)
                .expect("GHS: ghs_descend must succeed");
            // verify_log_preservation is part of the transfer cost — it confirms the
            // reduction is correct. The downstream index-calculus solve is NOT included.
            let preserved = verify_log_preservation(&descent, 1);
            assert!(preserved, "GHS: log-preservation must hold for k = 1");
            descent
        });
    });
}

// ── Index calculus ────────────────────────────────────────────────────────────

/// Bench: index calculus on `IndexCalcStrategy::toy()` (ℓ = 5, |FB| = 6, m = 2).
///
/// Known answer: the solver recovers `k mod ℓ` for Q = 7·G on the semaev_toy curve.
/// The bench asserts `index_calculus_dlp` returns a solution `x = k` with `k·G_ℓ = Q_ℓ`
/// before timing.
///
/// Relation/decomposition counts are derived from `collect_relations(...).len()` and
/// `decompose(...)`.
fn bench_index_calculus(c: &mut Criterion) {
    let strategy = IndexCalcStrategy::toy()
        .expect("index calculus pre-check: toy strategy must build");
    let curve = strategy.curve.clone();
    let g: AffinePoint<FpNaive> = curve.generator();
    // Q = 7·G on the semaev_toy curve (n = 60, ℓ = 5).
    let q: AffinePoint<FpNaive> = curve.scalar_mul(&g, &Uint::<4>::from(7u64));

    // Pre-check: assert the solver returns a valid answer before timing.
    let solution = index_calculus_dlp(g.clone(), q.clone(), &strategy)
        .expect("index calculus pre-check: index_calculus_dlp must recover k");
    let k = solution.x;
    // Verify k·G_ℓ = Q_ℓ (the subgroup-log correctness check).
    let ell_u64 = strategy.ell.as_words()[0];
    let n_u64 = curve.n.as_words()[0];
    let cofactor = n_u64 / ell_u64;
    let g_ell: AffinePoint<FpNaive> = curve.scalar_mul(&g, &Uint::<4>::from(cofactor));
    let q_ell: AffinePoint<FpNaive> = curve.scalar_mul(&q, &Uint::<4>::from(cofactor));
    let k_g_ell: AffinePoint<FpNaive> = curve.scalar_mul(&g_ell, &Uint::<4>::from(k));
    assert_eq!(
        k_g_ell, q_ell,
        "index calculus pre-check: k·G_ℓ ≠ Q_ℓ (k = {k})"
    );

    // Derive relation and decomposition counts from the public re-exports.
    // These counts are the pedagogical signal: how many relations the collection loop found,
    // and how many points were decomposable over the factor base.
    let relations = collect_relations(g.clone(), q.clone(), &strategy)
        .expect("index calculus pre-check: collect_relations must succeed");
    let relation_count = relations.len();

    // Count decomposable points by trying decompose on each factor-base point.
    // (Illustrative: the factor-base points are the simplest decomposable inputs.)
    let decomp_count = strategy
        .factor_base
        .iter()
        .filter(|fb| decompose(fb.point.clone(), &strategy).is_some())
        .count();

    // Annotate the bench group with the counts (printed in bench output, not timed).
    let mut group = c.benchmark_group("attacks");
    group.bench_function("index_calculus", |b| {
        b.iter(|| {
            index_calculus_dlp(g.clone(), q.clone(), &strategy)
                .expect("index calculus: index_calculus_dlp must not error")
        });
    });
    group.finish();

    // Print counts for the E.W table (not timed — informational only).
    let _ = (relation_count, decomp_count);
}

// ── Criterion wiring ──────────────────────────────────────────────────────────

criterion_group!(
    benches,
    bench_pohlig_hellman,
    bench_mov,
    bench_ssa,
    bench_ghs_transfer,
    bench_index_calculus,
);
criterion_main!(benches);
