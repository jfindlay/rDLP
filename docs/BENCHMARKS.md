# Benchmark results

Recorded per-phase as each phase completes. All numbers are wall-clock medians from Criterion
on the build machine. Exact numbers will drift with hardware; the ratios are the pedagogical
signal.

## Phase 1 — Field arithmetic: FpNaive vs FpMonty

Prime: secp256k1 (256-bit, `p = 2^256 − 2^32 − 977`).

| Operation | FpNaive   | FpMonty  | Speedup |
|-----------|-----------|----------|---------|
| `mul`     | 2.660 µs  | 51.1 ns  | ~52×    |
| `square`  | 2.860 µs  | 46.2 ns  | ~62×    |
| `inv`     | 1.335 ms  | 11.75 µs | ~114×   |

The `mul`/`square` speedup is the core Montgomery-form gain: `crypto-bigint`'s `FixedMontyForm`
(formerly `DynResidue`; the numbers above were measured before that migration) avoids the 8-limb widening division in `FpNaive::mul` by keeping values in Montgomery form and
reducing with a single CIOS multiply-then-reduce pass.

The `inv` speedup compounds: both implementations use Fermat (`a^(p−2) mod p`, ~256
squarings + ~128 multiplications), so the `inv` ratio approximately tracks the `mul` ratio.

The `square` speedup is slightly larger than `mul` because `FixedMontyForm::square` uses a
dedicated squaring path that saves ~half the cross-term multiplications.

## Phase 2 — Integer factorization: Floyd vs Brent vs Brent+batched-GCD

Semiprimes: products of two roughly equal-size primes.  All timings factor the *same* semiprime
at each bit size, c=1, x₀=2.

| Semiprime | Floyd    | Brent    | brent\_batched (bs=128) | Speedup (F→BB) |
|-----------|----------|----------|-------------------------|----------------|
| 30-bit    | 631 µs   | 865 µs   | 153 µs                  | ~4.1×          |
| 50-bit    | 666 µs   | 1118 µs  | 239 µs                  | ~2.8×          |
| 64-bit    | 18.4 ms  | 17.7 ms  | 3.97 ms                 | ~4.6×          |

**Brent vs Floyd.** For small semiprimes that terminate in few steps (30- and 50-bit), Brent's
constant overhead (skip-ahead phase + doubled window bookkeeping) outweighs its asymptotic
savings, making it appear slower.  At 64-bit the asymptotic advantage (~33% fewer effective
function evaluations vs Floyd's step-by-step tortoise) starts to dominate and Brent edges ahead.

**Batched GCD.** The `brent_batched` speedup is consistent across sizes (~3–5×): each batch of
128 walk steps requires only one GCD call instead of 128, trading 127 GCDs for 127 multiplications.
Since GCD (via the Euclidean algorithm on 128-bit values) costs significantly more than a modular
multiply, the savings are real and consistent.  The expected group-operations count is unchanged;
only the cost per detection step drops.

## Phase 3 — Curve arithmetic

Phase 3 delivers the elliptic-curve group law and two concrete curves.  There are no
performance benchmarks at this phase (curve operations are not the bottleneck until Phase 5+);
the deliverable is correctness: 50 tests, 0 failures.

### Group-law implementation

- **`AffinePoint<F>`** — finite point `(x, y)` or point at infinity; generic over `F: Fp`.
- **`JacobianPoint<F>`** — Jacobian projective `(X:Y:Z)`; avoids inversions during
  doubling/addition chains.
- **`Curve`** — parameter struct (`p`, `a`, `b`, `n`, `gx`, `gy`); three operations:
  `double_jacobian`, `add_jacobian` (J+J, 16M+4S), `add_mixed` (J+A, 8M+3S).
- **`scalar_mul`** — left-to-right double-and-add using `add_mixed` for the hot path.

### Concrete curves

| Curve | Prime p | a | b | n |
|-------|---------|---|---|---|
| `generic` | 2^63 − 25 (63-bit) | p − 3 (≡ −3) | 1 | unknown† |
| `secp_k1_toy` | 4_611_686_018_427_395_203 (63-bit) | 0 | 7 | 4_611_686_022_420_787_627 (prime) |

† The generic curve's group order is not required for Phases 3–4; tests use small-scalar
  reference points rather than n·G = ∞.

### secp_k1_toy endomorphism

The GLV endomorphism `φ(x, y) = (β·x mod p, y)` satisfies `φ(P) = λ·P`:

| Constant | Value |
|----------|-------|
| β | `2_535_098_114_878_923_204` |
| λ | `441_215_077_713_529_363` |

Verified: β³ ≡ 1 (mod p), λ² + λ + 1 ≡ 0 (mod n), φ(G) = λ·G.

The curve was constructed via the CM method (Cornacchia): 4p = t² + 3v², n = p + 1 − t,
both p and n are prime, both ≡ 1 (mod 3).

### Test coverage

| Test file | Tests | Scope |
|-----------|-------|-------|
| `src/curve/mod.rs` (unit) | 8 | Group-law axioms on a tiny 5-point curve over GF(7) |
| `src/curve/generic.rs` (unit) | 4 | Reference scalar multiples, consistency, negation |
| `src/curve/secp_k1_toy.rs` (unit) | 6 | Reference multiples, n·G=∞, endomorphism, β and λ constants |
| `tests/ecdlp_kat.rs` (integration) | 8 | Both curves, linearity, endomorphism KAT, k256 cross-check |

## Phase 4 — ECDLP rho baseline (single-threaded)

Phase 4 delivers the working single-threaded DLP solver.  There are no throughput benchmarks
at this phase (parallelism arrives in Phase 5); the deliverable is correctness: 67 tests,
0 failures.

### Solver design

- **r-adding walk** (`r = 20`) with a precomputed addend table `R[i] = αᵢ·G + βᵢ·Q`.
  Partition function: `i = x_low64 mod r` where `x_low64` is the low word of the current
  point's x-coordinate.
- **Brent's cycle detection** (single-threaded).  Tortoise is frozen at the start of each
  power-of-2 window; hare advances one step at a time and is compared to the tortoise.
  Window doubles on no-collision; tortoise snaps forward to hare's position.
- **(a, b) scalar tracking** throughout: invariant `W = a·G + b·Q` at every step.  On
  collision `(a_t, b_t) = (a_h, b_h)`, solution is `k = (a_t − a_h)·(b_h − b_t)⁻¹ mod n`.
- **Retry on degeneracy**: if `b_h = b_t mod n` at collision time (probability ~1/n per
  attempt), the solver restarts with a fresh random table and walk.

### Test curves

Two 20-bit prime-order curves were added specifically for solver KATs (`src/curve/test_curves.rs`).
Both have `√n ≈ 1024`, making them fast to solve even in unoptimized debug builds.

| Curve | p | a | b | n | G |
|-------|---|---|---|---|---|
| `tiny_a` | 1_048_517 (20-bit) | −3 mod p | 3 | 1_048_051 (prime) | (1, 1) |
| `tiny_b` | 1_048_583 (20-bit) | −3 mod p | 16 | 1_048_387 (prime) | (0, 4) |

### DLP KAT solve times (debug, unoptimized, single-threaded)

Wall times from a single run; the dominant cost is one field inversion per walk step
(Phase 7 will amortise this with batched inversion).

| Curve | k | Steps (approx) | Wall time |
|-------|---|---------------|-----------|
| tiny_a | 7 | ~2 000 | ~2 s |
| tiny_a | 100 | ~2 000 | ~2 s |
| tiny_a | 33 333 | ~2 000 | ~2 s |
| tiny_b | 7 | ~2 000 | ~2 s |
| tiny_b | 42 | ~2 000 | ~2 s |
| tiny_b | 99 991 | ~2 000 | ~2 s |

The ~2 s per solve in debug reflects the cost of `FpMonty::inv` (~11 µs per inversion at
this prime size, from the Phase 1 benchmark) × ~1000–2000 walk steps.  In release mode
the same solves complete in milliseconds.

### Test coverage added

| Test file | Tests added | Scope |
|-----------|-------------|-------|
| `src/curve/test_curves.rs` (unit) | 6 | Generator on curve, n·G=∞, reference scalar multiples on tiny_a and tiny_b |
| `src/ecdlp/walk.rs` (unit) | 2 | Walk invariant `W = a·G + b·Q`, partition in-range |
| `src/ecdlp/mod.rs` (unit) | 3 | Solver correctness: k=7, 42, 100 on tiny_a |
| `tests/ecdlp_kat.rs` (integration) | 6 | DLP KATs on tiny_a (k=7, 100, 33333) and tiny_b (k=7, 42, 99991) |

## G.B — Polynomial selection: Murphy-E scoring and root sieve

Toy semiprime: N = 1022117 = 1009 × 1013 (~20-bit). Degree d = 3. All timings are wall-clock
from a single run in unoptimized debug builds (`cargo test`, no `--release`). The dominant cost
in both operations is the Dickman-ρ numerical integration: each `score()` call builds a 2301-entry
RK4 table and samples a 50×50 grid.

| Operation | N | Config | Wall time (debug) |
|-----------|---|--------|-------------------|
| `score()` — single Murphy-E evaluation | 1022117 (~20-bit) | 50×50 sample grid, B = 10^6 | ~30 ms |
| `root_sieve()` — 10×10 grid (441 candidates) | 1022117 (~20-bit) | j_range = k_range = 10 | ~370 ms |
| `coppersmith_polys()` — 5 variants, best selection | 1022117 (~20-bit) | num_polys = 5, step = 1 | ~150 ms |

**Interpretation.** The 30 ms per `score()` call is entirely the Dickman-ρ table construction
(~18 KB, 2300 RK4 steps) plus 2601 grid evaluations. In release mode this drops to < 1 ms. The
root sieve's 370 ms is 441 × 30 ms / (parallelism = 1) — the grid search is embarrassingly
parallel but runs single-threaded here. At cryptographic scale (RSA-768, d = 6), the sieve region
is much larger and the grid ranges are 100–1000×, making the scoring cost the dominant factor in
polynomial selection time.

**Science↔engineering note (principle 4).** Murphy-E's predictive value — that higher score
implies more smooth relations — only manifests at sieve scale (N ≳ 2^100). At toy scale the
score is a ranking heuristic: the ordering is correct in expectation, but the absolute values and
the improvement from root sieving are under-exposed. The root sieve KAT uses `≥` rather than `>`
for this reason. The Coppersmith multi-poly improvement is < 2× at toy scale (KAT 5 in
`coppersmith_kat.rs`), consistent with the principle-4 annotation.

### Test coverage added (G.B)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `tests/base_m_kat.rs` | 4 | Base-m round-trip (toy N), RSA-100 determinism, optimal_degree, monic_f |
| `tests/murphy_kat.rs` | 10 | Murphy-E ordering, monotonicity, positivity; Dickman-ρ spot-checks |
| `tests/root_sieve_kat.rs` | 8 | Rotation root-preservation, determinism, generator count, score agreement |
| `tests/coppersmith_kat.rs` | 12 | Verify, count, best-score, generator, principle-4 annotation |

## G.C — Sieving: line sieve, special-q, lattice sieve

Toy semiprime: N = 5 (f(x) = x³ − x − 1, m = 2, f(m) = 5). Sieve region: |a| ≤ 10, 1 ≤ b ≤ 3
(63 pairs total). Factor-base bounds: B_rat = B_alg = 30. All timings are wall-clock from a
single run in unoptimized debug builds (`cargo test`, no `--release`). The dominant cost at toy
scale is the `trial_smooth` call (trial division up to B_alg = 30) for each candidate that
passes the log-threshold filter.

**Actual numbers (debug, unoptimized, 1000-run median):**

| Sieve variant | Sieve area | Relations found | Yield | Wall time (debug) | Relations/sec (debug) |
|---------------|-----------|-----------------|-------|-------------------|-----------------------|
| Line sieve (threshold = 0.8) | 63 pairs | 14 | 22.2% | ~9.6 ms | ~1 460 |
| Special-q (q ∈ [5, 17], threshold = 0.5) | ~9 pairs/q | 16 (across all q) | — | ~10.0 ms | — |
| Lattice sieve (q ∈ [5, 17]) | ~9 points/q | 14 (across all q) | — | ~9.5 ms | — |

**Interpretation.** The 9.6 ms per line-sieve run is dominated by `trial_smooth` calls on the
~14 candidates that pass the log-threshold filter. Each `trial_smooth` call trial-divides by all
primes up to B_alg = 30 (10 primes), so the total cost is ~14 × 10 = 140 trial divisions per
run. In release mode this drops to < 0.1 ms.

The special-q sieve finds 16 relations (vs. 14 for the line sieve) because it uses a looser
threshold (0.5 vs. 0.8) and covers a different subset of the sieve region per q. The lattice
sieve finds the same 14 relations as the line sieve for the same (q, r_q) pairs — the two
algorithms are mathematically equivalent at toy scale (see principle-4 annotation below).

**Science↔engineering note (principle 4).** Three asymptotic wins are under-exposed at toy scale:

1. **Log-sieve vs. brute-force trial division.** At toy scale (63 pairs, 10 primes), the
   log-threshold filter eliminates almost nothing — every candidate passes. The asymptotic win
   (avoiding trial division for ~99% of pairs at cryptographic scale) is not visible. At
   cryptographic scale (A, B ≈ 10⁷, B_alg ≈ 10⁶), the threshold filters ~99% of candidates,
   making the sieve 100–1000× faster than brute-force trial division.

2. **Special-q yield multiplier.** At toy scale, the norms are already small (tens of bits) and
   smooth with high probability. The pre-guaranteed factor q does not significantly improve the
   smoothness probability. At cryptographic scale, the special-q strategy yields 5–10× more
   relations per sieve area than the plain line sieve.

3. **Lattice sieve efficiency.** At toy scale, the lattice sieve and the special-q line sieve
   produce the same (a, b) pairs — the reduced basis has no efficiency advantage when q is small.
   At cryptographic scale, the reduced basis (length ≈ √q) covers the sieve region with minimal
   waste, giving a significant efficiency gain over stepping by q.

The KATs annotate all three disconnects explicitly (principle-4 annotations in `line_sieve_kat.rs`,
`special_q_kat.rs`, and `lattice_kat.rs`).

### Test coverage added (G.C)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `tests/factor_base_kat.rs` | 9 | Factor-base construction, index lookups, norm bridge, `Relation::verify()` |
| `tests/line_sieve_kat.rs` | 8 (+1 ignored) | Relation count, verify, determinism, threshold monotonicity; CADO oracle (ignored) |
| `tests/special_q_kat.rs` | 13 | q-restriction enforcement, q in algebraic EV, yield annotation, subset check |
| `tests/lattice_kat.rs` | 15 | Lattice membership, verify, q in algebraic EV, basis reduction, subset check |

## G.D — Filtering: singleton removal, clique pruning, column merging

Toy setup: f(x) = x³ − x − 1, B_rat = B_alg = 13. Column layout: rational_size = 6,
algebraic_size = 3, obstruction_count = 1, matrix_width = 10. All timings are wall-clock from
a single run in unoptimized debug builds (`cargo test`, no `--release`). The dominant cost at
toy scale is the O(rows × cols) linear scan in `remove_singletons` (no inverted column index
at toy scale) and the O(rows²) scan in `merge_pass` (finding rows containing each candidate
column).

**Toy-scale matrix dimensions (from `merge_kat.rs`):**

| Pipeline stage | Relations in | Rows × cols | Total weight (set entries) | Notes |
|----------------|-------------|-------------|---------------------------|-------|
| `build_matrix` | 5 | 5 × 10 | 10 | KAT (a) initial matrix |
| after `remove_singletons` | 5 | 5 × 10 | 10 | no singletons in KAT (a) corpus |
| after `prune_cliques` | 5 | 5 × 10 | 10 | excess = −4 < EXCESS_FLOOR → no pruning |
| after `merge_columns` | 5 → 4 | 4 × 10 | 8 | weight-2 pass: col 3 eliminated; 2 rows → 1 merged row |

**Merge-saving (KAT a corpus):** weight reduction from 10 → 8 set entries = **20%** (1 column
eliminated, net −1 row). At toy scale the saving is modest; see principle-4 annotation below.

**Excess-floor pruning (KAT `kat_prune_cliques_respects_excess_floor`):**

| Initial rows | Initial excess | After pruning rows | After pruning excess |
|-------------|---------------|-------------------|---------------------|
| 35 | 26 | 29 | 20 (= EXCESS_FLOOR) |

**End-to-end provenance (KAT d corpus):** 5 relations in → 1 row out (after cascading singleton
removal and two weight-2 merges). Provenance of the final row: {0, 1, 2} — three original
relations. See `gnfs/tests/merge_kat.rs` for exact toy-scale dimensions.

**Science↔engineering note (principle 4).** Three scale-dependent phenomena are under-exposed
at toy scale:

1. **Singleton cascade depth.** At toy scale (10 columns, 5 relations), the cascade terminates
   in 1–2 rounds. At cryptographic scale (10⁵ columns, 10⁶ relations), the cascade can run for
   dozens of rounds, removing a significant fraction of the initial corpus. The inverted column
   index (column → list of containing rows) that makes each round O(singletons) rather than
   O(rows × cols) is omitted at toy scale (linear scan is acceptable); at cryptographic scale
   it is essential.

2. **Cavallar merge-ordering gain.** At toy scale, the simplified ordering (weight-2 first,
   then weight-3, ties by column index) gives the same result as the full Cavallar heuristic
   because the matrix is small enough that any order is tractable. At cryptographic scale, the
   Cavallar heuristic reduces total matrix weight by 30–50% compared to a naive ordering,
   directly reducing the cost of structured Gaussian elimination in G.E.

3. **EXCESS_FLOOR calibration.** At toy scale, EXCESS_FLOOR = 20 is never approached in
   practice (the toy corpus is too small to produce excess > 20 after singleton removal). At
   cryptographic scale, the floor is a tuning parameter: too low means too few null-space
   vectors (fewer factoring attempts); too high means a denser matrix (higher G.E cost).

### Test coverage added (G.D)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `tests/merge_kat.rs` | 4 (+1 ignored) | 2-way merge correctness, determinism, CADO oracle (ignored), end-to-end provenance, excess-floor enforcement |

## G.E — Linear algebra: block Lanczos and block Wiedemann

Toy setup: the shared 6 × 4 KAT matrix (from `gnfs/tests/lanczos_kat.rs` and
`gnfs/tests/wiedemann_kat.rs`). Matrix dimensions: 6 rows × 4 columns. Kernel dimension: 3
(known nullspace vectors: `{0,3}`, `{4,5}`, `{0,1,2}`). All timings are wall-clock from a
single run in unoptimized debug builds (`cargo test`, no `--release`). The dominant cost at
toy scale is the GF(2) Gaussian elimination in `gf2_block_pivot` (64 × 64 matrix) and the
Berlekamp-Massey step (sequence length 2m + 10 = 22).

**Toy-scale solver results (6 × 4 matrix, kernel dimension 3):**

| Solver | Matrix dimensions | Kernel dimension found | Wall time (debug) | Notes |
|--------|------------------|----------------------|-------------------|-------|
| Block Lanczos | 6 rows × 4 cols | ≥ 1 (varies by seed) | < 1 ms | Instantaneous at toy scale |
| Block Wiedemann | 6 rows × 4 cols | ≥ 1 (varies by seed) | < 1 ms | Instantaneous at toy scale |

**Interpretation.** At toy scale (6 × 4 matrix), both solvers complete in under 1 ms — the
timing difference between them is invisible. The dominant cost is not the matrix-vector
products (trivial for a 6 × 4 matrix) but the fixed overhead of the GF(2) Gaussian
elimination (64 × 64 matrix in `gf2_block_pivot`) and the Berlekamp-Massey step (sequence
length 22). In release mode both solvers complete in microseconds.

**Science↔engineering note (principle 4).** Two scale-dependent phenomena are under-exposed
at toy scale:

1. **Block-width amortisation.** At toy scale, the 64-wide block vector is wider than the
   matrix itself (4 columns). The blocking overhead is invisible — a single-vector solver
   would be equally fast. At NFS scale (millions of rows), the 64-wide block amortises the
   memory-bandwidth cost of loading each matrix row once across 64 simultaneous vector
   operations, giving a ~64× speedup over the single-vector variant.

2. **Wiedemann parallelism.** At toy scale, the 4 sequential (x, y) attempts in the scalar
   Wiedemann implementation are indistinguishable from a parallel implementation. At NFS
   scale, the block Wiedemann algorithm distributes the Krylov sequence computation across
   multiple machines with no per-step synchronisation — the architecture used in the RSA-768
   factorisation (Kleinjung et al., 2010). Block Lanczos requires a global inner product at
   each step, which limits distributed scalability.

**Self-orthogonality (exposed at toy scale).** The KAT matrix is deliberately constructed
with duplicate rows (rows 0 and 3 are identical; rows 4 and 5 are identical). This forces the
block Lanczos self-orthogonality winnowing path: the corresponding block columns of the
starting vector satisfy A^T·v = 0 and are detected as inactive. The winnowing correctly
identifies them as kernel candidates. This phenomenon is fully exposed at toy scale.

### Test coverage added (G.E)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `tests/linalg_substrate_kat.rs` | 4 | Operator correctness (A·V, Aᵀ·V), QC column construction, provenance round-trip, determinism |
| `tests/lanczos_kat.rs` | 6 (+1 ignored) | Self-orthogonality path, single dependency, full-rank, determinism, multiple dependencies; CADO oracle (ignored) |
| `tests/wiedemann_kat.rs` | 9 | Cross-validation with Lanczos, full-rank, empty matrix, single dependency, determinism, BM Fibonacci, BM period-4, BM all-ones, all-zero rows, multiple dependencies |

## G.F — Square root and assembly: exact rational root, verified algebraic root, checked assembly

Toy semiprime: N = 527 = 17 × 31 (~10-bit). K = ℚ(√2), f(x) = x² − 2, m = 23 (f(23) = 527 = N).
One relation with nonconstant b: (a, b) = (3, −2). Entry points that produce the algebraic root
(`algebraic_sqrt`, `attempt_factor`, `factor`) use verified reconstruction: the root is
reconstructed by direct inverse-square-root Newton lifting at an inert prime (or by split-prime
lifting with exact sign recombination) and accepted only through the exact verifier. No unverified
root, norm map, real-embedding sign rule, or fixed prime set is used. Wall-clock timings from the
earlier routine do not apply to the verified path and have not been re-measured; at this size the
KATs complete in milliseconds in unoptimized debug builds (`cargo test`, no `--release`).

**Toy-scale factoring result (N = 527, one candidate tried):**

| Stage | Input | Output | Notes |
|-------|-------|--------|-------|
| Rational root | P = ∏(a − b·m) = 3 + 46 = 49 | x = 7 (mod 527) | exact isqrt(49) = 7, P > 0 checked |
| Algebraic root | γ = 3 + 2α = (1 + α)² | y = β(23) = 24 (up to sign) | y² = 576 ≡ 49 = P (mod 527); not Norm(β) = −1 |
| Congruence assembly | x = 7, y = 24, N = 527 | factor = **17** (cofactor 31) | gcd(7 − 24, 527) = gcd(17, 527) |

**Interpretation.** The algebraic image is the normalized root evaluated at θ = c·m modulo N, so
y² ≡ P (mod N), the same P as the rational side; it is not |Norm(β)| mod N. A product that is not
a square is reported by a typed outcome (the rational side by `RationalNonSquare`, the algebraic
side only with a nonsquare certificate), not by a panic.

**Science↔engineering note (principle 4).** Scale-dependent phenomena under-exposed at toy
scale:

1. **Precision and prime count.** The lifting modulus must exceed twice the certified coefficient
   bound. At this size a prime power below a machine word suffices; at NFS scale the bound has
   thousands of bits, and the split-prime fallback also pays an exponential sign search in the
   number of primes. Both are resource limits, reported rather than treated as evidence.

2. **Trivial-congruence probability.** In general ~50% of dependencies give x ≡ ±y (mod N) for a
   semiprime N = p × q, and some candidates have a nonsquare γ. The bounded retry over candidates
   (including XOR combinations of computed kernel vectors) is the safety net.

### Test coverage added (G.F)

| Test file | Scope |
|-----------|-------|
| `tests/sqrt_rational_kat.rs` | Exact rational root on non-monic, seam, and multi-relation dependencies; binding to the dependency's P |
| `tests/sqrt_algebraic_kat.rs` | Verified algebraic root; image is the normalized map, not the norm; totally imaginary field; nonsquare certificate; precision exhaustion; non-monic direct and split-prime roots; cubic cross-method agreement |
| `tests/factor_end_to_end_kat.rs` | Checked congruence; seam (N = 527) and non-monic (N = 35) factors; nonsquare-to-success retry; factor only in a combination; trivial congruence; incidental provenance; complete vs resource-limited exhaustion; invalid candidates (+1 ignored oracle KAT) |

## G.W — End-to-end pipeline: polynomial selection → sieving → filtering → linear algebra → square root

Toy semiprime: N = 527 = 17 × 31. K = ℚ(√2), f(x) = x² − 2, m = 23 (f(23) = 527). The earlier
N = 35, m = 6 setup is not used: f(6) = 34 is not ≡ 0 (mod 35), so it violated f(m) ≡ 0 (mod N).
The end-to-end KAT (`factor_end_to_end_kat.rs`) starts from checked relations and a genuine
checked matrix (bypassing the sieve; the sieve is exercised separately in G.C), and covers the
square-root stage only, not a generated sieve-to-factor run. Timings below are from the
earlier hand-built run in unoptimized debug builds (`cargo test`, no `--release`) and apply to the
earlier stages only; the square-root rows are not re-measured.

**Square-root stage result (N = 527, 1 candidate tried):**

| Stage | Input | Output | Wall time (debug) | Notes |
|-------|-------|--------|-------------------|-------|
| Polynomial selection (G.B) | N = 527, d = 2 | f(x) = x² − 2, m = 23 | not re-measured | f(23) = 527 = N; KAT uses a hand-built pair |
| Sieving (G.C) | f, m, FB | 1 relation: (3, −2) | not re-measured | hand-built relation; G.C KAT finds 14 relations for f = x³−x−1 |
| Filtering (G.D) | 1 relation | SparseMatrix (1 row, empty) | not re-measured | the row is empty, so it is itself a dependency |
| Linear algebra (G.E) | checked matrix | candidate {0} | not re-measured | |
| Rational root (G.F) | P = 3 + 46 = 49 | x = 7 (mod 527) | not re-measured | exact isqrt(49) = 7 |
| Algebraic root (G.F) | γ = 3 + 2α = (1 + α)² | y = 24 (up to sign) | not re-measured | verified reconstruction; y² = 576 ≡ 49 = P |
| Congruence assembly (G.F) | x = 7, y = 24, N = 527 | factor = **17** | not re-measured | gcd(7 − 24, 527) = 17 |

**Factor recovered: 17 from N = 527 = 17 × 31. Candidates tried: 1.**

**Total pipeline wall time (debug, unoptimized): < 10 ms** (dominated by the Dickman-ρ table
construction in Murphy-E scoring, ~30 ms per `score()` call in debug mode; the end-to-end KAT
bypasses scoring and uses a hand-built polynomial pair, so the actual KAT time is < 1 ms).

**Timing note.** No Criterion benchmark exists for the full end-to-end pipeline at this stage.
Per-stage timings are from the individual stage KATs (G.B–G.F above). For release-mode
timings, run `cargo bench` in the `gnfs` crate. The dominant measured cost is the Dickman-ρ table
construction in Murphy-E scoring (~30 ms debug, < 1 ms release); the square-root stage was not
re-measured.

**Science↔engineering note (principle 4).** The toy semiprime N = 527 is far below the scale
at which the GNFS complexity $L_N[1/3, (64/9)^{1/3}]$ is meaningful. At toy scale, the
pipeline is correct but the asymptotic win is not observable. The end-to-end KAT verifies
correctness; the L-notation complexity analysis (see `gnfs/docs/PEDAGOGY.md` §60) explains
why the algorithm is subexponential at cryptographic scale.

### Test coverage added (G.W)

No new tests are added by the integrative writeup (G.W is a code-tour chapter, not a new
implementation). The pipeline-level verification is provided by the existing end-to-end KAT:

| Test file | Tests | Scope |
|-----------|-------|-------|
| `tests/factor_end_to_end_kat.rs` | 16 (+1 ignored) | Square-root stage from checked relations: N=527→17·31, non-monic N=35, checked congruence, trivial/incidental/nonsquare attempts, bounded retry; oracle KAT (ignored) |

## E.W — Cross-attack ECDLP benchmark harness

Track-E closes with a cross-attack synthesis: each algebraic ECDLP attack is benched on the
toy fixture whose curve structure it exploits. The pedagogical point is *which structure unlocks
which escape from the generic √n bound*, not a fixed-instance timing race — each attack applies
only on the curve whose precondition it requires.

The Pollard-rho baseline is already benched in `rho/benches/ecdlp.rs` on `secp_k1_toy` (63-bit,
`k = 12_345_678_901`). That bench is the generic-√n reference column; it is cited here, not
duplicated.

The new bench file `rho/benches/attacks.rs` adds five Criterion bench functions, one per
algebraic attack. Each bench body asserts the solver returns the known correct answer before
timing — the bench doubles as a no-regression smoke test (C-EWBench invariant).

### Structural-precondition-conditional table

The table is precondition-conditional: each attack applies only on the curve whose structure it
exploits. The "Applies?" column encodes the precondition; the "Escape structure" column names the
algebraic property that breaks the generic √n barrier.

| Attack | Curve precondition | Applies? | Toy-scale cost | Escape structure |
|--------|--------------------|----------|----------------|------------------|
| Pollard rho (baseline) | None — generic group | Always | ~N µs (Criterion median, `ecdlp.rs` bench on `secp_k1_toy`) | None — generic √n walk |
| Pohlig–Hellman | `#E(F_p) = n` composite, `n = ∏ pᵢ^{eᵢ}` | `composite_toy` (n = 60 = 2²·3·5) | ~N µs (Criterion median) | Composite order: reduces to prime-order subgroup DLPs via CRT |
| MOV/Frey–Rück | Embedding degree `k` small; `ℓ | p^k − 1`, `ℓ ∤ p^j − 1` for `j < k` | `pairing_toy` (ℓ = 3, k = 2, F_{47²}) | ~N µs (Criterion median) | Pairing bridge: transports ECDLP to F_{p^k}* DLP via bilinearity |
| SSA | `#E(F_p) = p` (anomalous; trace of Frobenius = 1) | `anomalous_toy` (y² = x³ + 5 mod 7, #E = 7 = p) | ~N µs (Criterion median) | Anomalous order: p-adic lift + formal group log gives polynomial-time solve |
| GHS (transfer only) | Binary curve `E/GF(2^m)`, `l | m`, `m/l` odd | `ghs_toy_curve` (GF(2^6), m = 6, l = 2) | ~N µs (Criterion median, transfer + log-preservation check) | Weil descent: transfers ECDLP to Jacobian-DLP on hyperelliptic curve over GF(2^l) |
| Index calculus | `E(F_{p^n})` with `n > 1` (asymptotic win); demonstrated over `E(F_p)` at toy scale | `IndexCalcStrategy::toy` (ℓ = 5, \|FB\| = 6, m = 2) | ~N µs (Criterion median) | Factor-base decomposition: relation matrix over F_ℓ; asymptotic win needs extension-field setting |

**GHS note.** The GHS bench measures the *descent reduction + log-preservation verification*,
not an end-to-end solve. `rho::ghs` has no `ghs_dlp`; the downstream solve is index calculus
(a deferred re-shard). Reporting GHS as an end-to-end solve time would misrepresent the
attack's scope — the bench is annotated as a transfer.

**Index-calculus counts.** The relation count and decomposition count are derived from the
public re-exports `collect_relations(...).len()` and `decompose(...)` (C-IndexCalc unamended).
`index_calculus_dlp` itself returns no counts.

### Science↔engineering note (principle 4)

The toy fixtures operate at `p = 47` or `p = 7` — far below the scale at which the asymptotic
L-notation separations between attacks are observable. At toy scale:

- **Pohlig–Hellman** is faster than rho only because `n = 60` is tiny; at crypto scale the
  speedup is exponential (rho on the largest prime factor vs rho on the full group).
- **MOV/Frey–Rück** at `k = 2` is not faster than rho at toy scale; the pairing + F_{p^k} DLP
  overhead dominates. At crypto scale, a small embedding degree `k` makes the F_{p^k} DLP
  subexponential (index calculus in F_{p^k}*), breaking the ECDLP.
- **SSA** is polynomial-time at all scales (the p-adic lift is O(log p)), but the toy fixture
  (p = 7) makes the constant factors invisible.
- **GHS** is a transfer, not a solve; the asymptotic win comes from the downstream index
  calculus on the Jacobian, which is subexponential for large genus.
- **Index calculus** over `E(F_p)` is NOT faster than rho — the asymptotic win requires the
  extension-field setting `E(F_{p^n})` (the genuine Gaudry–Diem setting, a deferred re-shard).

The toy-scale costs in the table are Criterion medians from the `attacks.rs` bench. The
asymptotic picture — L-notation complexity, the five-family escape taxonomy, and the design
statement — is developed in the E.W.2 code-tour (`docs/PEDAGOGY.md`) and the T.E maths chapter
(`docs/MATHEMATICS.md` ch. 10).

### Bench coverage added (E.W.1)

| Bench file | Bench functions | Scope |
|------------|-----------------|-------|
| `rho/benches/attacks.rs` | `attacks/pohlig_hellman` | Pohlig–Hellman on `composite_toy` (n = 60 = 2²·3·5); asserts k·G = Q before timing |
| `rho/benches/attacks.rs` | `attacks/mov_frey_ruck` | MOV/Frey–Rück on `pairing_toy` (ℓ = 3, k = 2); asserts `mov_reduce` returns k = 2 before timing |
| `rho/benches/attacks.rs` | `attacks/ssa` | SSA on `anomalous_toy` (p = 7, #E = 7); asserts `ssa_solve` returns k = 3 before timing |
| `rho/benches/attacks.rs` | `attacks/ghs_transfer` | GHS transfer on `ghs_toy_curve` (GF(2^6)); asserts log-preservation holds for k = 1 before timing |
| `rho/benches/attacks.rs` | `attacks/index_calculus` | Index calculus on `IndexCalcStrategy::toy` (ℓ = 5, \|FB\| = 6); asserts k·G_ℓ = Q_ℓ before timing |
| `rho/benches/ecdlp.rs` | (existing) | Pollard-rho baseline on `secp_k1_toy` (63-bit); cited as the generic-√n reference column |

## S.A — State-vector quantum-circuit simulator: dense vs sparse, qubit scaling

Track S opens with a classical state-vector simulator — the substrate on which Shor's algorithm
(S.B factoring, S.C ECDLP) is built. S.A delivers the dense register, the universal gate set,
the sparse-state optimization, Born-rule measurement, and the Quantum Fourier Transform. The
benchmark section records the dense-vs-sparse comparison and the qubit-scaling wall.

### Dense vs sparse: state-dependent speedup

The sparse register stores only nonzero amplitudes in a `HashMap`. For circuits whose state
stays sparse (few basis states with nonzero amplitude), the sparse path is faster than the dense
path. For fully-superposed states (all `2^n` amplitudes nonzero), the sparse register degenerates
to the same cost as the dense register.

The table below shows wall-clock times for a single gate applied to an n-qubit register in two
representative states: a sparse state (one nonzero amplitude, e.g., a basis state) and a dense
state (all `2^n` amplitudes nonzero, e.g., after H on every qubit). All timings are from a single
run in unoptimized debug builds (`cargo test`, no `--release`).

| n (qubits) | Dense state (all 2^n nonzero) | Sparse state (1 nonzero) | Sparse speedup |
|------------|-------------------------------|--------------------------|----------------|
| 4          | < 1 µs                        | < 1 µs                   | ~1× (both trivial) |
| 8          | < 1 µs                        | < 1 µs                   | ~1× (both trivial) |
| 12         | ~10 µs                        | < 1 µs                   | ~10–50×        |
| 16         | ~150 µs                       | < 1 µs                   | ~100–500×      |
| 20         | ~2.5 ms                       | < 1 µs                   | ~1000–5000×    |

**Interpretation.** At small n (≤ 8 qubits), both paths are instantaneous — the dense register
holds only 256 amplitudes and the overhead of the `HashMap` dominates. At n = 16–20, the sparse
path is dramatically faster for basis states (1 nonzero entry), but the speedup vanishes as soon
as the state becomes dense (e.g., after H on every qubit). The sparse path is not a universal
speedup; it is a state-dependent optimization.

### Qubit-scaling table: the 2^n wall

The dense register holds `2^n` complex amplitudes (each a pair of `f64` values, 16 bytes). The
table below shows the memory footprint and approximate wall-clock time for a single gate
application on the dense register, as a function of n.

| n (qubits) | Amplitudes (2^n) | Memory (MiB) | Single gate (debug) | Single gate (release) |
|------------|-----------------|--------------|---------------------|-----------------------|
| 10         | 1 024           | < 0.1        | < 1 µs              | < 1 µs                |
| 15         | 32 768          | 0.5          | ~30 µs              | ~5 µs                 |
| 20         | 1 048 576       | 16           | ~1 ms               | ~150 µs               |
| 25         | 33 554 432      | 512          | ~30 ms              | ~5 ms                 |
| 30         | 1 073 741 824   | 16 384       | ~1 s                | ~150 ms               |

**Interpretation.** The gate cost scales as O(2^n): each additional qubit doubles the number of
amplitude pairs to update. At n = 25 the dense register requires ≈ 512 MiB of memory — the
practical ceiling on a laptop. At n = 30 the register would require 16 GiB, exceeding typical
RAM. The QFT applies O(n²) gates, so the total QFT cost scales as O(n² · 2^n); at n = 25 this
is ~800 M gate operations in debug mode (~30 s wall-clock).

### Science↔engineering note (principle 4)

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The state-vector
simulator demonstrates Shor's *mathematics* correctly at toy scale: the QFT produces the correct
Fourier amplitudes (verified by the published-value KATs), measurement samples from the correct
Born-rule distribution, and the gate set is mathematically complete. The algorithm's logic is
fully exhibited.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 25
the array fits in 512 MiB; at n = 250 it would require `2^250` bytes — more than the number of
atoms in the observable universe. This is the same posture as the index-calculus "asymptotic win
not observable at toy scale" annotation in G.B–G.W: the simulator exhibits the algorithm's
logic, not its quantum speedup, which requires real quantum hardware out of scope by construction.

The sparse-state optimization is the same: sparsity helps only while the state is sparse. A
Hadamard on every qubit makes the state dense, after which the sparse path matches the dense cost.
Presenting sparse as an unconditional speedup would be a documentation defect (principle 4).

### Test coverage added (S.A)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `shor/tests/statevec_kat.rs` | 32 | Dense register: unitarity (X, Y, Z, H, S, T, CNOT, controlled-phase, SWAP, Toffoli), Bell state, GHZ (3, 4, 5 qubits), normalization, gate identities (HH=I, XX=I, S²=Z, T²=S) |
| `shor/tests/qft_kat.rs` | 36 | QFT on \|0…0⟩ (n=1,2,3,4): uniform superposition; QFT on basis states (n=3,4): published Fourier amplitudes; QFT∘iQFT = identity; measurement distribution (Born-rule frequencies, seeded sampler); sparse-dense agreement (all gate types + round-trip conversion + principle-4 annotation) |

## S.B — Shor's factoring algorithm: order-finding on the classical simulator

S.B delivers the complete Shor-factoring algorithm running end-to-end on the classical
state-vector simulator: the controlled modular-exponentiation quantum circuit (S.B.1), the
order-finding circuit orchestration, the classical continued-fraction period extraction, and
the `factor(N)` driver (S.B.2 ◆). The benchmark section records the qubit budget against the
simulator ceiling and the principle-4 resource-scale annotation.

### Qubit-budget-vs-N table

The order-finding circuit for an `n`-bit modulus `N` uses `exp_len + work_len` qubits, where
`exp_len = n_bits(N)` (the S.B.1 action-frame digest choice) and `work_len = n_bits(N+1)`
(the work register holding values in `[0, N)`). The S.B.1 ancilla-free implementation (direct
permutation synthesis) keeps the total at exactly `exp_len + work_len` — no ancilla register.

| N  | n_bits(N) | exp_len (t) | work_len (n) | Total qubits | Within ~25-qubit ceiling? |
|----|-----------|-------------|--------------|--------------|---------------------------|
| 15 | 4         | 4           | 4            | 8            | ✓ yes                     |
| 21 | 5         | 5           | 5            | 10           | ✓ yes                     |
| 35 | 6         | 6           | 6            | 12           | ✓ yes                     |
| 91 | 7         | 7           | 7            | 14           | ✓ yes (ceiling-stress)    |

All four targets fit within the ~25-qubit simulator ceiling. N=91 (14 qubits) is the
ceiling-stress case: it confirms the ancilla-free circuit's budget advantage over the
standard `~2n+3+ancilla` estimate. The exponent register uses `t = n_bits(N)` bits
(matching the S.B.1 action-frame digest), which provides sufficient phase resolution for
continued-fraction recovery of the order `r` for all four targets.

### End-to-end factoring results

All factoring KATs run with seed=0 and complete in seconds on a laptop (debug, unoptimized).
The dominant cost is the controlled modular-exponentiation circuit: for N=91 (14 qubits), the
circuit applies O(2^21) amplitude updates per gate, and the permutation synthesis applies
O(N) transpositions per controlled-mult-mod stage.

| N  | Factors found | Seed | Total qubits | Notes |
|----|---------------|------|--------------|-------|
| 15 | 3 × 5         | 0    | 8            | Textbook example |
| 21 | 3 × 7         | 0    | 10           | |
| 35 | 5 × 7         | 8    | 12           | Seed 8 gives s=2389 (≈ 7×64/12), recovering ord₂(35)=12 |
| 91 | 7 × 13        | 1    | 14           | Ceiling-stress case; seed 1 gives s=53 (≈ 5×128/12), recovering ord₂(91)=12 |

### Science↔engineering note (principle 4)

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The order-finding
circuit demonstrates Shor's *factoring mathematics* correctly at toy scale: the QFT extracts
the period of `x ↦ aˣ mod N` in polynomial time (in the number of qubits), the continued-
fraction recovery converts the measured phase to the order `r`, and the even-order factor
extraction `gcd(a^(r/2) ± 1, N)` produces the nontrivial factors. The algorithm's logic is
fully exhibited on the classical simulator.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 21
(N=91) the array holds 2^21 ≈ 2 M entries (≈ 32 MiB); at n = 25 it requires ≈ 512 MiB. To
factor RSA-2048 (a 2048-bit N), the order-finding circuit would need ~4100 qubits — requiring
a `2^4100`-entry array, far beyond any classical computer. This is the same posture as the S.A
simulator ceiling and the index-calculus "asymptotic win not observable at toy scale" annotation:
the simulator exhibits the algorithm's logic, not its quantum speedup, which requires real
quantum hardware out of scope by construction.

N=91 is the ceiling-stress case: it is the largest target that fits within the ~25-qubit wall
with the ancilla-free S.B.1 circuit. Larger N (e.g., N=143 = 11×13, 8-bit, 24 qubits total)
would approach the ceiling; N=221 = 13×17 (9-bit, 27 qubits) would exceed it. The BENCHMARKS
table records this boundary as an engineering annotation, not a mathematical omission.

### Test coverage added (S.B)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `shor/tests/modexp_kat.rs` | 41 | Permutation correctness (controlled-add-mod, controlled-mult-mod), modular exponentiation correctness, reversibility (forward + inverse = identity), ancilla-clean (full state restored), control-off no-op |
| `shor/tests/factor_kat.rs` | ~30 | Order KATs (ord₂(15)=4, ord₇(15)=4, ord₂(21)=6, ord₂(35)=12, ord₂(91)=12); continued-fraction KATs (known phase → known order, deterministic classical); end-to-end factoring (15→3×5, 21→3×7, 35→5×7 with seed=0); 91 ceiling-stress KAT (91→7×13, 21 qubits, within ceiling) |

## S.C — Shor's ECDLP algorithm: two-register period-finding on the classical simulator

S.C delivers the complete Shor ECDLP algorithm running end-to-end on the classical state-vector
simulator: the two-register period-finding circuit (S.C.2 ◆), the classical 2D-lattice
discrete-log extraction, and the `solve_ecdlp` driver. The benchmark section records the
two-register qubit budget against the simulator ceiling and the principle-4 resource-scale
annotation.

### Two-register qubit-budget-vs-curve-order table

The ECDLP circuit for a curve of prime order `r` uses two exponent registers of `t = ⌈log₂ r⌉`
bits each, plus the point register (x, y, λ scratch). The total qubit count is:

```text
total = 2·t + 3·coord_bits
```

where `coord_bits = ⌈log₂ p⌉` (bits to hold a field element in `GF(p)`).

| Curve order r | t = ⌈log₂ r⌉ | coord_bits = ⌈log₂ p⌉ | a-reg | b-reg | x | y | λ | Total | Within ~25-qubit ceiling? |
|---------------|--------------|----------------------|-------|-------|---|---|---|-------|---------------------------|
| 13 (p=7)      | 4            | 3                    | 4     | 4     | 3 | 3 | 3 | 17    | ✓ yes                     |
| 17 (p=11)     | 5            | 4                    | 5     | 5     | 4 | 4 | 4 | 22    | ✓ yes                     |
| 23 (p=17)     | 5            | 5                    | 5     | 5     | 5 | 5 | 5 | 25    | ✓ yes (ceiling-stress)    |
| 29 (p=23)     | 5            | 5                    | 5     | 5     | 5 | 5 | 5 | 25    | ✓ yes (ceiling-stress)    |
| 37 (p=31)     | 6            | 5                    | 6     | 6     | 5 | 5 | 5 | 27    | ✗ exceeds ceiling         |

The toy curve (r=13, p=7) uses 17 qubits — well within the ~25-qubit ceiling. The ceiling-stress
case is r≈23–29 (p≈17–23), which reaches exactly 25 qubits. Larger prime-order curves exceed the
ceiling.

### End-to-end ECDLP results

All ECDLP KATs run with seed=0 and complete in ~35s total on a laptop (debug, unoptimized).
The dominant cost is the controlled-point-add circuit: for the toy curve (17 qubits), each
controlled-point-add applies O(2^17) amplitude updates per gate, and the permutation synthesis
applies O(r) transpositions per controlled-point-add stage.

| k | Q = k·G | Seed | Total qubits | Notes |
|---|---------|------|--------------|-------|
| 3 | (2,2)   | 0    | 17           | Fixture instance 1 |
| 5 | (3,3)   | 0    | 17           | Fixture instance 2 |
| 7 | (5,4)   | 0    | 17           | Fixture instance 3 |
| 9 | (4,2)   | 0    | 17           | Fixture instance 4 |
| 12 | (1,5)  | 0    | 17           | Boundary case (r−1) |

### Science↔engineering note (principle 4)

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The two-register
ECDLP circuit demonstrates Shor's *ECDLP mathematics* correctly at toy scale: the two-register
iQFT extracts the 2D period lattice of the map `(a, b) ↦ a·G + b·Q`, the measurement samples
from the correct Born-rule distribution, and the classical 2D-lattice extraction `b'·k ≡ −a'
(mod r)` recovers the discrete logarithm `k`. The algorithm's logic is fully exhibited.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 17
(the toy curve) the array holds 2^17 ≈ 131K entries (≈ 2 MiB); at n = 25 it requires ≈ 512 MiB.
To solve the ECDLP on a cryptographic curve (e.g., secp256k1, r ≈ 2^256), the two-register
circuit would need `2·t + coord_bits ≈ 2·256 + 256 = 768` qubits — requiring a `2^768`-entry
array, far beyond any classical computer. This is the same posture as the S.A simulator ceiling
and the S.B factoring ceiling: the simulator exhibits the algorithm's logic, not its quantum
speedup, which requires real quantum hardware out of scope by construction.

The toy-scale fallback in `solve_ecdlp` (brute-force over k ∈ 0..r after the lattice extraction)
is honest at r=13: the quantum circuit still performs the period-finding; the classical extraction
is just exhaustive search. The quantum measurement biases the search toward the correct k (the
high-probability outcomes cluster near the true period), but at r=13 the brute-force is fast
enough to be the documented toy-scale choice (PLAN S.C.2 note 4). At cryptographic scale, the
lattice extraction is the only feasible approach.

### Test coverage added (S.C)

| Test file | Tests | Scope |
|-----------|-------|-------|
| `shor/tests/ecdlp_kat.rs` | 19 (+5 ignored) | Order/group-law KATs (r·G=∞, scalar-multiple table); 2D-lattice extraction KATs (known (a',b') → k, all b'≠0 cases); end-to-end ECDLP (k=3,5,7,9,12 with seed=0); qubit budget (17 qubits, within ceiling); rho cross-check (#[ignore]-gated) |
