//! Factorization benchmark: Floyd vs Brent vs Brent+batched-GCD.
//!
//! Measures total wall time to factor three semiprimes of increasing size
//! (30, 50, 64 bits) with each variant.  The pedagogical signal is the
//! wall-clock ordering: Floyd > Brent > brent_batched, with the batched-GCD
//! advantage growing with semiprime size (larger N ⟹ more rho steps per run
//! ⟹ more GCD savings per run).
//!
//! We do not include `factor` (parallel multi-c) in the microbenchmark because
//! Criterion runs each function many times in the same process and rayon thread
//! pools interact poorly with per-iteration timing; the parallel speedup is
//! instead demonstrated in the KAT timing output.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rho::factor::{brent, brent_batched, floyd};

/// Step bound far above what any input below needs, so the bound never ends a run early.
const STEPS: u64 = 1 << 28;

/// Three semiprimes spanning 30–64 bits.  Each is a product of two primes
/// of roughly equal size so that rho's expected steps are ~N^(1/4).
struct Input {
    label: &'static str,
    n: u128,
}

const INPUTS: &[Input] = &[
    Input { label: "30-bit", n: 32_452_843 * 32_452_867 },
    Input { label: "50-bit", n: 33_554_467 * 33_554_473 },
    Input { label: "64-bit", n: 4_294_967_291 * 4_294_967_311 },
];

fn bench_floyd(c: &mut Criterion) {
    let mut group = c.benchmark_group("floyd");
    for input in INPUTS {
        group.bench_with_input(BenchmarkId::new("n", input.label), &input.n, |b, &n| {
            b.iter(|| floyd(n, 1, 2, STEPS).expect("floyd failed"));
        });
    }
    group.finish();
}

fn bench_brent(c: &mut Criterion) {
    let mut group = c.benchmark_group("brent");
    for input in INPUTS {
        group.bench_with_input(BenchmarkId::new("n", input.label), &input.n, |b, &n| {
            b.iter(|| brent(n, 1, 2, STEPS).expect("brent failed"));
        });
    }
    group.finish();
}

fn bench_brent_batched(c: &mut Criterion) {
    let mut group = c.benchmark_group("brent_batched");
    for input in INPUTS {
        group.bench_with_input(BenchmarkId::new("n", input.label), &input.n, |b, &n| {
            b.iter(|| brent_batched(n, 1, 2, 128, STEPS).expect("brent_batched failed"));
        });
    }
    group.finish();
}

criterion_group!(factor_benches, bench_floyd, bench_brent, bench_brent_batched);
criterion_main!(factor_benches);
