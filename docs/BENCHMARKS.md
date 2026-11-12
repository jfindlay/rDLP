# Benchmark results

Results are grouped by component. All numbers are wall-clock medians from Criterion
on the build machine. Exact numbers will drift with hardware; the ratios are the pedagogical
signal.

## Field arithmetic: FpNaive vs FpMonty

Prime: secp256k1 (256-bit, `p = 2^256 − 2^32 − 977`).

| Operation | FpNaive   | FpMonty  | Speedup |
|-----------|-----------|----------|---------|
| `mul`     | 2.660 µs  | 51.1 ns  | ~52×    |
| `square`  | 2.860 µs  | 46.2 ns  | ~62×    |
| `inv`     | 1.335 ms  | 11.75 µs | ~114×   |

The `mul`/`square` speedup is the core Montgomery-form gain: `crypto-bigint`'s `FixedMontyForm`
(formerly `DynResidue`; the numbers above were measured before that migration) avoids the 8-limb
widening division in `FpNaive::mul` by keeping values in Montgomery form and reducing with a single
CIOS multiply-then-reduce pass.

The `inv` speedup compounds: both implementations use Fermat (`a^(p−2) mod p`, ~256
squarings + ~128 multiplications), so the `inv` ratio approximately tracks the `mul` ratio.

The `square` speedup is slightly larger than `mul` because `FixedMontyForm::square` uses a
dedicated squaring path that saves ~half the cross-term multiplications.

## Integer factorization: Floyd vs Brent vs Brent+batched-GCD

Semiprimes: products of two roughly equal-size primes.  All timings factor the *same* semiprime
at each bit size, c=1, x₀=2.

| Semiprime | Floyd    | Brent    | brent\_batched (bs=128) | Speedup (F→BB) |
|-----------|----------|----------|-------------------------|----------------|
| 30-bit    | 631 µs   | 865 µs   | 153 µs                  | ~4.1×          |
| 50-bit    | 666 µs   | 1118 µs  | 239 µs                  | ~2.8×          |
| 64-bit    | 18.4 ms  | 17.7 ms  | 3.97 ms                 | ~4.6×          |

Each bench passes an explicit step bound of 2²⁸ evaluations, far above what these inputs need
(`rho/benches/factor.rs`), and unwraps the result, so a run that spent the bound would abort the
bench rather than be timed. The bench calls the primitives (`floyd`, `brent`, `brent_batched`), not
`factor`, whose outcomes are typed (`docs/PEDAGOGY.md` §5).

**Brent vs Floyd.** For small semiprimes that terminate in few steps (30- and 50-bit), Brent's
constant overhead (skip-ahead phase + doubled window bookkeeping) outweighs its asymptotic
savings, making it appear slower.  At 64-bit the asymptotic advantage (~33% fewer effective
function evaluations vs Floyd's step-by-step tortoise) starts to dominate and Brent edges ahead.

**Batched GCD.** The `brent_batched` speedup is consistent across sizes (~3–5×): each batch of
128 walk steps requires only one GCD call instead of 128, trading 127 GCDs for 127 multiplications.
Since GCD (via the Euclidean algorithm on 128-bit values) costs significantly more than a modular
multiply, the savings are real and consistent.  The expected group-operations count is unchanged;
only the cost per detection step drops.

## Curve arithmetic

The curve-arithmetic component delivers the elliptic-curve group law and two concrete curves.  There
are no performance benchmarks for it (curve operations are not the bottleneck until the parallel
distinguished-point solver); the deliverable is correctness, covered by the unit and integration tests listed under "Test
coverage" below.

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
| `generic` | 2^63 − 25 (63-bit) | p − 3 (≡ −3) | 1 | 561_103 (prime; order of the stored generator's subgroup)† |
| `secp_k1_toy` | 4_611_686_018_427_395_203 (63-bit) | 0 | 7 | 4_611_686_022_420_787_627 (prime) |
| `prime_order_25` | 16_777_259 (25-bit) | p − 3 (≡ −3) | 72 | 16_777_729 (prime; the whole group)‡ |

† `generic` is a 63-bit field carrying a subgroup of about 2¹⁹ elements, not a 63-bit group. Its
  stored generator `G = (3_569_091_914_278_015_591, 7_073_522_384_094_063_206)` is the multiple
  `8_218_965_179_160 · P₀` of the original point `P₀ = (3, 821_487_384_573_098_969)`, and the tests
  prove in-tree that `G` is on the curve, `G ≠ ∞`, `561_103·G = ∞` with 561_103 prime (so `G` has
  exact order 561_103), and that the multiple relation holds. PARI/GP supplied the rest as
  discovery provenance, not certified here: `#E = 9_223_372_037_844_426_960`,
  `ord P₀ = 4_611_686_018_922_213_480`, and that 561_103 is the largest prime factor of the
  latter. The group order is the order of the subgroup, so a solver base point and target must lie
  in `⟨G⟩`; `P₀` is rejected, not projected. A rho collision there takes about a thousand steps
  (arithmetic from `√n`), so it is a small runnable subgroup, not a representative large-group or
  parallel-performance fixture, and its `a = p − 3`, `b = 1` form has no GLV endomorphism.

‡ `prime_order_25` is the fixture the `ecdlp` bench uses (see the cross-attack section). Its
  numbers were found with PARI/GP 2.17.3 (`ellcard`) by a bounded search: `a = −3`, `p` the first
  prime at or above 2²⁴, `b = 1, 2, …, 200` in order skipping the singular `b = 2`, stopping at the
  first `b` whose `#E` is prime and `≠ p` (that was `b = 72`, the 71st curve tried); `G` is the
  point with the smallest `x`. The tests prove in-tree that `G` is on the curve, `n` is prime,
  `n·G = ∞`, and, via the Hasse bound, that `n` is the order of the whole group.

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

| Test file | Scope |
|-----------|-------|
| `src/curve/mod.rs` (unit) | Group-law axioms on a tiny 5-point curve over GF(7) |
| `src/curve/generic.rs` (unit) | Reference scalar multiples, consistency, negation |
| `src/curve/secp_k1_toy.rs` (unit) | Reference multiples, n·G=∞, endomorphism, β and λ constants |
| `tests/ecdlp_kat.rs` (integration) | Both curves, linearity, endomorphism KAT, k256 cross-check |

## ECDLP rho baseline (single-threaded)

The single-threaded baseline delivers the working DLP solver.  There are no throughput benchmarks
of it on the KAT curves; the deliverable is correctness. The `ecdlp` bench (see the cross-attack
section) runs it on the 25-bit prime-order `prime_order_25` group.

### Solver design

- **r-adding walk** (`r = 20`) with a precomputed addend table `R[i] = αᵢ·G + βᵢ·Q`.
  Partition function: `i = x_low64 mod r` where `x_low64` is the low word of the current
  point's x-coordinate.
- **Brent's cycle detection** (single-threaded).  Tortoise is frozen at the start of each
  power-of-2 window; hare advances one step at a time and is compared to the tortoise.
  Window doubles on no-collision; tortoise snaps forward to hare's position.
- **(a, b) scalar tracking** throughout: invariant `W = a·G + b·Q` at every step.  On
  collision `(a_t, b_t) = (a_h, b_h)`, solution is `k = (a_t − a_h)·(b_h − b_t)⁻¹ mod n`.
- **Attempts, not retries without bound**: if `b_h = b_t mod n` at collision time (probability
  ~1/n per attempt), or an attempt reaches its step cap (2²⁹ − 1 steps), that attempt is spent and
  the solver restarts with a fresh random table and walk, up to `max_retries` attempts. Running out
  is `BudgetExhausted`, which proves nothing about `Q`.
- **Verified result**: the solver returns `k` only after checking `k·G == Q`, with the collision
  scalars as evidence (`brent_collision_evidence_is_rechecked_and_reproducible`,
  `brent_zero_retries_is_budget_exhausted` in `rho/tests/ecdlp_kat.rs`).

### Test curves

Two 20-bit prime-order curves exist specifically for solver KATs (`src/curve/test_curves.rs`).
Both have `√n ≈ 1024`, making them fast to solve even in unoptimized debug builds.

| Curve | p | a | b | n | G |
|-------|---|---|---|---|---|
| `tiny_a` | 1_048_517 (20-bit) | −3 mod p | 3 | 1_048_051 (prime) | (1, 1) |
| `tiny_b` | 1_048_583 (20-bit) | −3 mod p | 16 | 1_048_387 (prime) | (0, 4) |

### DLP KAT solve times (debug, unoptimized, single-threaded)

Wall times from a single run; the dominant cost is one field inversion per walk step
(batched inversion amortises this).

| Curve | k | Steps (approx) | Wall time |
|-------|---|---------------|-----------|
| tiny_a | 7 | ~2 000 | ~2 s |
| tiny_a | 100 | ~2 000 | ~2 s |
| tiny_a | 33 333 | ~2 000 | ~2 s |
| tiny_b | 7 | ~2 000 | ~2 s |
| tiny_b | 42 | ~2 000 | ~2 s |
| tiny_b | 99 991 | ~2 000 | ~2 s |

The ~2 s per solve in debug reflects the cost of `FpMonty::inv` (~11 µs per inversion at
this prime size, from the field-arithmetic benchmark) × ~1000–2000 walk steps.  In release mode
the same solves complete in milliseconds.

### Test coverage added

| Test file | Scope |
|-----------|-------|
| `src/curve/test_curves.rs` (unit) | Generator on curve, n·G=∞, reference scalar multiples on tiny_a and tiny_b |
| `src/ecdlp/walk.rs` (unit) | Walk invariant `W = a·G + b·Q`, partition in-range |
| `src/ecdlp/mod.rs` (unit) | Solver correctness: k=7, 42, 100 on tiny_a |
| `tests/ecdlp_kat.rs` (integration) | DLP KATs on tiny_a (k=7, 100, 33333) and tiny_b (k=7, 42, 99991) |

## Polynomial selection: Murphy-E scoring and root sieve

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

**Science↔engineering note.** Murphy-E's predictive value — that higher score
implies more smooth relations — only manifests at sieve scale (N ≳ 2^100). At toy scale the
score is a ranking heuristic: the ordering is correct in expectation, but the absolute values and
the improvement from root sieving are under-exposed. The root sieve KAT uses `≥` rather than `>`
for this reason. The Coppersmith multi-poly improvement is < 2× at toy scale (KAT 5 in
`coppersmith_kat.rs`), consistent with the scale annotation.

### Test coverage added (polynomial selection)

| Test file | Scope |
|-----------|-------|
| `tests/base_m_kat.rs` | Base-m round-trip (toy N), RSA-100 determinism, optimal_degree, monic_f |
| `tests/murphy_kat.rs` | Murphy-E ordering, monotonicity, positivity; Dickman-ρ spot-checks |
| `tests/root_sieve_kat.rs` | Rotation root-preservation, determinism, generator count, score agreement |
| `tests/coppersmith_kat.rs` | Verify, count, best-score, generator, scale annotation |

## Sieving: line sieve, special-q, lattice sieve

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
algorithms are mathematically equivalent at toy scale (see scale annotation below).

**Science↔engineering note.** Three asymptotic wins are under-exposed at toy scale:

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

The KATs annotate all three disconnects explicitly (scale annotations in `line_sieve_kat.rs`,
`special_q_kat.rs`, and `lattice_kat.rs`).

### Test coverage added (sieving)

| Test file | Scope |
|-----------|-------|
| `tests/factor_base_kat.rs` | Factor-base construction, index lookups, norm bridge, `Relation::verify()` |
| `tests/line_sieve_kat.rs` | Relation count, verify, determinism, threshold monotonicity; an `#[ignore]`d CADO-NFS placeholder whose assertion compares the sieve with itself (no oracle is run) |
| `tests/special_q_kat.rs` | q-restriction enforcement, q in algebraic EV, yield annotation, subset check |
| `tests/lattice_kat.rs` | Lattice membership, verify, q in algebraic EV, basis reduction, subset check |

## Filtering: singleton removal, clique pruning, column merging

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
eliminated, net −1 row). At toy scale the saving is modest; see scale annotation below.

**Excess-floor pruning (KAT `kat_prune_cliques_respects_excess_floor`):**

| Initial rows | Initial excess | After pruning rows | After pruning excess |
|-------------|---------------|-------------------|---------------------|
| 35 | 26 | 29 | 20 (= EXCESS_FLOOR) |

**End-to-end provenance (KAT d corpus):** 5 relations in → 1 row out (after cascading singleton
removal and two weight-2 merges). Provenance of the final row: {0, 1, 2} — three original
relations. See `gnfs/tests/merge_kat.rs` for exact toy-scale dimensions.

**Science↔engineering note.** Three scale-dependent phenomena are under-exposed
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
   directly reducing the cost of structured Gaussian elimination in the linear-algebra stage.

3. **EXCESS_FLOOR calibration.** At toy scale, EXCESS_FLOOR = 20 is never approached in
   practice (the toy corpus is too small to produce excess > 20 after singleton removal). At
   cryptographic scale, the floor is a tuning parameter: too low means too few null-space
   vectors (fewer factoring attempts); too high means a denser matrix (higher linear-algebra cost).

### Test coverage added (filtering)

| Test file | Scope |
|-----------|-------|
| `tests/merge_kat.rs` | 2-way merge correctness, determinism, end-to-end provenance, excess-floor enforcement; an `#[ignore]`d CADO-NFS placeholder whose body is `unimplemented!()` (no oracle is run) |

## Linear algebra: block Lanczos and block Wiedemann

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

**Science↔engineering note.** Two scale-dependent phenomena are under-exposed
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

### Test coverage added (linear algebra)

| Test file | Scope |
|-----------|-------|
| `tests/linalg_substrate_kat.rs` | Operator correctness (A·V, Aᵀ·V), root-character columns on a checked matrix, provenance round-trip, determinism, and the individual root-character tests |
| `tests/lanczos_kat.rs` | Self-orthogonality path, single dependency, full-rank, determinism, multiple dependencies; an `#[ignore]`d test named for CADO-NFS that checks kernel vectors and provenance on a hand-built matrix and never invokes CADO-NFS |
| `tests/wiedemann_kat.rs` | Cross-validation with Lanczos, full-rank, empty matrix, single dependency, determinism, BM Fibonacci, BM period-4, BM all-ones, all-zero rows, multiple dependencies |

## Square root and assembly: exact rational root, verified algebraic root, checked assembly

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
| Algebraic root | γ = 3 + 2α = (1 + α)² | y = β(c·m) = β(23) = 24 (up to sign; c = 1) | y² = 576 ≡ 49 = P (mod 527); not Norm(β) = −1 |
| Congruence assembly | x = 7, y = 24, N = 527 | factor = **17** (cofactor 31) | gcd(7 − 24, 527) = gcd(17, 527) |

**Interpretation.** The algebraic image is the normalized root evaluated at θ = c·m modulo N, so
y² ≡ P (mod N), the same P as the rational side; it is not |Norm(β)| mod N. A product that is not
a square is reported by a typed outcome (the rational side by `RationalNonSquare`, the algebraic
side only with a nonsquare certificate), not by a panic.

**Science↔engineering note.** Scale-dependent phenomena under-exposed at toy
scale:

1. **Precision and prime count.** The lifting modulus must exceed twice the certified coefficient
   bound. At this size a prime power below a machine word suffices; at NFS scale the bound has
   thousands of bits, and the split-prime fallback also pays an exponential sign search in the
   number of primes. Both are resource limits, reported rather than treated as evidence.

2. **Trivial-congruence probability.** In general ~50% of dependencies give x ≡ ±y (mod N) for a
   semiprime N = p × q, and some candidates have a nonsquare γ. The bounded retry over candidates
   (including XOR combinations of computed kernel vectors) is the safety net.

### Test coverage added (square root and assembly)

| Test file | Scope |
|-----------|-------|
| `tests/sqrt_rational_kat.rs` | Exact rational root on non-monic, seam, and multi-relation dependencies; binding to the dependency's P |
| `tests/sqrt_algebraic_kat.rs` | Verified algebraic root; image is the normalized map, not the norm; totally imaginary field; nonsquare certificate; precision exhaustion; non-monic direct and split-prime roots; cubic cross-method agreement |
| `tests/factor_end_to_end_kat.rs` | Checked congruence; seam (N = 527) and non-monic (N = 35) factors; nonsquare-to-success retry; factor only in a combination; trivial congruence; incidental provenance; complete vs resource-limited exhaustion; invalid candidates; `#[ignore]`d PARI/GP cross-checks (`kat_pari_*`, which need `gp`) and one `#[ignore]`d placeholder (`todo!()`) for an 80–100-bit challenge that is not an oracle check |

## End-to-end pipeline: polynomial selection → sieving → filtering → linear algebra → square root

Toy semiprime: N = 527 = 17 × 31. K = ℚ(√2), f(x) = x² − 2, m = 23 (f(23) = 527). The earlier
N = 35, m = 6 setup is not used: f(6) = 34 is not ≡ 0 (mod 35), so it violated f(m) ≡ 0 (mod N).
The end-to-end KAT (`factor_end_to_end_kat.rs`) starts from checked relations and a genuine
checked matrix (bypassing the sieve; the sieve is exercised separately by the sieving KATs), and
covers the square-root stage only. The run from generated relations to a verified factor is
`factor_pipeline_kat.rs`, on a different fixture (N = 45113 = 197 · 229); no wall-clock figure is
recorded for it here. The table below describes the N = 527 square-root stage, and none of its
rows has been re-measured.

**Square-root stage result (N = 527, 1 candidate tried):**

| Stage | Input | Output | Wall time (debug) | Notes |
|-------|-------|--------|-------------------|-------|
| Polynomial selection | N = 527, d = 2 | f(x) = x² − 2, m = 23 | not re-measured | f(23) = 527 = N; KAT uses a hand-built pair |
| Sieving | f, m, FB | 1 relation: (3, −2) | not re-measured | hand-built relation; the sieving KAT finds 14 relations for f = x³−x−1 |
| Filtering | 1 relation | SparseMatrix (1 row, empty) | not re-measured | the row is empty, so it is itself a dependency |
| Linear algebra | checked matrix | candidate {0} | not re-measured | |
| Rational root | P = 3 + 46 = 49 | x = 7 (mod 527) | not re-measured | exact isqrt(49) = 7 |
| Algebraic root | γ = 3 + 2α = (1 + α)² | y = 24 (up to sign) | not re-measured | verified reconstruction; y² = 576 ≡ 49 = P |
| Congruence assembly | x = 7, y = 24, N = 527 | factor = **17** | not re-measured | gcd(7 − 24, 527) = 17 |

**Factor recovered: 17 from N = 527 = 17 × 31. Candidates tried: 1.**

**Timing note.** No end-to-end wall time is recorded: the table rows are not re-measured, and no
Criterion benchmark exists for the `gnfs` crate. The only per-stage figures in this file are those
of the earlier sections (polynomial selection through linear algebra) and were measured on
their own fixtures. The end-to-end KAT bypasses Murphy-E scoring, whose Dickman-ρ table
construction is the one measured cost in polynomial selection (~30 ms per `score()` call in debug
builds), and uses a hand-built polynomial pair.

**Science↔engineering note.** The toy semiprime N = 527 is far below the scale
at which the GNFS complexity $L_N[1/3, (64/9)^{1/3}]$ is meaningful. At toy scale, the
pipeline is correct but the asymptotic win is not observable. The end-to-end KAT verifies
correctness; the L-notation complexity analysis (see `gnfs/docs/PEDAGOGY.md` §60) explains
why the algorithm is subexponential at cryptographic scale.

### Test coverage added (end-to-end pipeline)

This section defines no tests of its own (the end-to-end writeup is a code-tour chapter, not an
implementation). The pipeline-level verification is provided by the end-to-end KAT:

| Test file | Scope |
|-----------|-------|
| `tests/factor_end_to_end_kat.rs` | Square-root stage from checked relations: N=527→17·31, non-monic N=35, checked congruence, trivial/incidental/nonsquare attempts, bounded retry; `#[ignore]`d PARI/GP cross-checks (`kat_pari_*`) and an `#[ignore]`d `todo!()` placeholder that is not an oracle check |

## Cross-attack ECDLP benchmark harness

The algebraic ECDLP attacks are benched together: each attack runs on the toy fixture whose curve
structure it exploits. The pedagogical point is *which structure unlocks
which escape from the generic √n bound*, not a fixed-instance timing race — each attack applies
only on the curve whose precondition it requires.

The Pollard-rho baseline is benched in `rho/benches/ecdlp.rs`, which covers Brent, `solve_dp`
(1, 2, 4 walkers), `solve_dp_negmap` and `solve_dp_batch`. Setup:

- **Group.** `prime_order_25`: `y² = x³ − 3x + 72` over the 25-bit prime field
  `p = 16_777_259`, whole group of prime order `n = 16_777_729` (25 bits), generator
  `G = (1, 6_299_712)`. Provenance is in the "Concrete curves" footnote. A plain rho collision
  takes about 1.25·√n ≈ 5·10³ steps (arithmetic from the group order, not a measurement).
- **Target.** `k = 12_345_678_901 mod n`, `Q = k·G`.
- **Tuning.** `theta = 4` (one point in 16 distinguished); distinguished-point step budget
  2²² ≈ 4.2·10⁶ (a safety cap, about 800 times the expected plain-walk cost); Brent 20 attempts
  of at most 2²⁹ − 1 steps; walkers 1, 2, 4 for `solve_dp`, 2 for the negation-map and
  batched-inversion variants; batch size 16; Criterion sample size 10; fixed seed.
- **Check before timing.** Every bench solves its instance once, outside the timed loop, and
  asserts the solve succeeded and returned `k`, so a solver that cannot finish fails the bench
  rather than timing an error.
- **Not benched.** `solve_dp_glv` supports only `secp_k1_toy` (group order about 2⁶², on the order
  of 10⁹ rho steps per solve) and the fixture has no GLV endomorphism, so there is no GLV bench.

The group is small, so the bench says nothing about performance at cryptographic sizes, and each
comparison holds only for this one group. The group size was chosen so that every variant finishes
in about two seconds or less per solve in a release build; a sizing run on a 31-bit prime-order
group showed `solve_dp_negmap` not finishing within 60–120 s (see the bench's module comment for
the sizing observations). The sizing runs were wall-clock checks that each solve finishes, not
benchmark results. They also showed `solve_dp_negmap` taking several times the step count of
`solve_dp` on the 25-bit group (about 3·10⁴–7·10⁴ against 2·10³–10⁴), well above the expected
cost, so the √2 saving expected from the negation map is not observed here and nothing in this
section claims a negation-map or batched-inversion speedup. The cause is not investigated.

An earlier version of the bench ran on the 63-bit `secp_k1_toy` group (about 2.7·10⁹ expected
steps per solve), where the solves could not finish in practice; no timing from that setup exists.
No Criterion median is recorded for the current setup either: the `~N µs` entries below stand for
figures not yet measured.

The bench file `rho/benches/attacks.rs` has five Criterion bench functions, one per
algebraic attack. Each bench body asserts the solver returns the known correct answer before
timing — the bench doubles as a no-regression smoke test (every benched attack still solves, or
for GHS transfers, its toy instance).

### Structural-precondition-conditional table

The table is precondition-conditional: each attack applies only on the curve whose structure it
exploits. The "Applies?" column encodes the precondition; the "Escape structure" column names the
algebraic property that breaks the generic √n barrier.

| Attack | Curve precondition | Applies? | Toy-scale cost | Escape structure |
|--------|--------------------|----------|----------------|------------------|
| Pollard rho (baseline) | None — generic group | Always | ~N µs (Criterion median, `ecdlp.rs` bench on `prime_order_25`; not yet measured) | None — generic √n walk |
| Pohlig–Hellman | `#E(F_p) = n` composite, `n = ∏ pᵢ^{eᵢ}` | `composite_toy` (n = 60 = 2²·3·5) | ~N µs (Criterion median) | Composite order: reduces to prime-order subgroup DLPs via CRT |
| MOV/Frey–Rück | Embedding degree `k` small; `ℓ | p^k − 1`, `ℓ ∤ p^j − 1` for `j < k` | `pairing_toy` (ℓ = 3, k = 2, F_{47²}) | ~N µs (Criterion median) | Pairing bridge: transports ECDLP to F_{p^k}* DLP via bilinearity |
| SSA | `#E(F_p) = p` (anomalous; trace of Frobenius = 1) | `anomalous_toy` (y² = x³ + 5 mod 7, #E = 7 = p) | ~N µs (Criterion median) | Anomalous order: p-adic lift + formal group log gives polynomial-time solve |
| GHS (transfer only) | Binary curve `E/GF(2^m)`, `l | m`, `m/l` odd | `ghs_toy_curve` (GF(2^6), m = 6, l = 2) | ~N µs (Criterion median, transfer + log-preservation check) | Weil descent: transfers ECDLP to Jacobian-DLP on hyperelliptic curve over GF(2^l) |
| Index calculus | `E(F_{p^n})` with `n > 1` (asymptotic win); demonstrated over `E(F_p)` at toy scale | `IndexCalcStrategy::toy` (ℓ = 5, \|FB\| = 6, m = 2) | ~N µs (Criterion median) | Factor-base decomposition: relation matrix over F_ℓ; asymptotic win needs extension-field setting |

**GHS note.** The GHS bench measures the *descent reduction + log-preservation verification*,
not an end-to-end solve. `rho::ghs` has no `ghs_dlp`; the downstream solve is index calculus on
the Jacobian, which is not implemented. Reporting GHS as an end-to-end solve time would
misrepresent the attack's scope — the bench is annotated as a transfer.

**Index-calculus counts.** `index_calculus_dlp` returns an `IndexCalcSolution` carrying `x`
(the verified `log_G(Q) mod ℓ`), the relation combination `u` that produced it, and the relations
that combination combines, so the relation count used is `solution.relations.len()`. The bench
also prints the number of relations `collect_relations(...)` returns and the number of
factor-base points for which `decompose(...)` (which returns `Option<Vec<FbPoint>>`) finds a
decomposition; these are informational counts printed by the bench, not asserted values. The
bench pre-check asserts `x·G_ℓ == Q_ℓ` for `Q = 7·G`; the recovery mechanism itself (left kernel,
every residue `k mod ℓ`, failure classes) is tested in `rho/tests/index_calculus_kat.rs`
(`left_kernel_annihilates_exponents`, `recovers_every_residue`,
`undetermined_when_no_combination_sees_b`, `corrupted_relation_is_inconsistent`).

### Science↔engineering note

The toy fixtures operate at `p = 47` or `p = 7` — far below the scale at which the asymptotic
L-notation separations between attacks are observable. At toy scale:

- **Pohlig–Hellman** is faster than rho only because `n = 60` is tiny; at crypto scale the
  speedup is exponential (rho on the largest prime factor vs rho on the full group).
- **MOV/Frey–Rück** at `k = 2` is not faster than rho at toy scale; the pairing + F_{p^k} DLP
  overhead dominates. At crypto scale, a small embedding degree `k` makes the F_{p^k} DLP
  subexponential (index calculus in F_{p^k}*), breaking the ECDLP.
- **SSA** is polynomial in log p as an algorithm, but this implementation is not: it multiplies by
  p as p − 1 repeated additions and confirms anomalousness by an O(p) point count. At the toy
  fixture (p = 7) the constant factors are invisible either way.
- **GHS** is a transfer, not a solve; the asymptotic win comes from the downstream index
  calculus on the Jacobian, which is subexponential for large genus.
- **Index calculus** over `E(F_p)` is NOT faster than rho — the asymptotic win requires the
  extension-field setting `E(F_{p^n})` (the genuine Gaudry–Diem setting), which this code does
  not implement.

The `~N µs` cells in the table are placeholders: no timing is recorded in this document. The
benches report Criterion medians when run (`cargo bench --bench attacks`); this document makes no
speed claim for any attack. The asymptotic picture — L-notation complexity and the five-family
escape taxonomy — is developed in the attack code-tour (`docs/PEDAGOGY.md`) and the algebraic
attacks chapter of the textbook (`docs/MATHEMATICS.md` ch. 10).

### Bench coverage

| Bench file | Bench functions | Scope |
|------------|-----------------|-------|
| `rho/benches/attacks.rs` | `attacks/pohlig_hellman` | Pohlig–Hellman on `composite_toy` (n = 60 = 2²·3·5); asserts k·G = Q before timing |
| `rho/benches/attacks.rs` | `attacks/mov_frey_ruck` | MOV/Frey–Rück on `pairing_toy` (ℓ = 3, k = 2); asserts the verified solution's `k_mod_ell` is 2 (the log modulo ℓ = 3) before timing |
| `rho/benches/attacks.rs` | `attacks/ssa` | SSA on `anomalous_toy` (p = 7, #E = 7); asserts `ssa_solve` returns `Ok(3)` before timing |
| `rho/benches/attacks.rs` | `attacks/ghs_transfer` | GHS transfer on `ghs_toy_curve` (GF(2^6)); asserts log-preservation holds for k = 1 before timing |
| `rho/benches/attacks.rs` | `attacks/index_calculus` | Index calculus on `IndexCalcStrategy::toy` (ℓ = 5, \|FB\| = 6); asserts x·G_ℓ = Q_ℓ before timing |
| `rho/benches/ecdlp.rs` | `ecdlp/brent`, `ecdlp/solve_dp`, `ecdlp/negmap_vs_dp`, `ecdlp/batch_vs_negmap` | Pollard-rho baseline on `prime_order_25` (25-bit field, prime group order 16_777_729); cited as the generic-√n reference column; asserts each solve succeeds before timing; no GLV bench (see above) |

## Shor simulator — State-vector quantum-circuit simulator: dense vs sparse, qubit scaling

The Shor chapter opens with a classical state-vector simulator — the substrate on which Shor's
algorithm (factoring and ECDLP) is built. The simulator provides the dense register, the
universal gate set, the sparse-state optimization, Born-rule measurement, and the Quantum Fourier
Transform. This section records the dense-vs-sparse comparison and the qubit-scaling wall.

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

### Science↔engineering note

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The state-vector
simulator demonstrates Shor's *mathematics* correctly at toy scale: the QFT produces the correct
Fourier amplitudes (verified by the published-value KATs), measurement samples from the correct
Born-rule distribution, and the gate set is mathematically complete. The algorithm's logic is
fully exhibited.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 25
the array fits in 512 MiB; at n = 250 it would require `2^250` bytes — more than the number of
atoms in the observable universe. This is the same posture as the other scale annotations in
this document (for example, that index calculus over `E(F_p)` shows no asymptotic win at toy
scale): the simulator exhibits the algorithm's logic, not its quantum speedup, which requires
real quantum hardware out of scope by construction.

The sparse-state optimization is the same: sparsity helps only while the state is sparse. A
Hadamard on every qubit makes the state dense, after which the sparse path matches the dense cost.
Presenting sparse as an unconditional speedup would be a documentation defect.

### Test coverage (Shor simulator)

| Test file | Scope |
|-----------|-------|
| `shor/tests/statevec_kat.rs` | Dense register: unitarity (X, Y, Z, H, S, T, CNOT, controlled-phase, SWAP, Toffoli), Bell state, GHZ (3, 4, 5 qubits), normalization, gate identities (HH=I, XX=I, S²=Z, T²=S) |
| `shor/tests/qft_kat.rs` | QFT on \|0…0⟩ (n=1,2,3,4): uniform superposition; QFT on basis states (n=3,4): published Fourier amplitudes; QFT∘iQFT = identity; measurement distribution (Born-rule frequencies, seeded sampler); sparse-dense agreement (all gate types + round-trip conversion + the dense-state scale annotation) |

## Shor factoring — order-finding on the classical simulator

This section covers the complete Shor-factoring algorithm running end-to-end on the classical
state-vector simulator: the controlled modular-exponentiation quantum circuit, the order-finding
circuit orchestration, the classical continued-fraction period extraction, and the `factor(N)`
driver. The section records the qubit budget against the simulator ceiling and the
resource-scale annotation.

### Qubit-budget-vs-N table

The order-finding circuit for a modulus `N` uses `exp_len + work_len` qubits, where
`exp_len = n_bits(N) = ⌈log₂ N⌉` (the exponent-register width chosen for this circuit) and
`work_len = n_bits(N+1)` (the work register holding values in `[0, N)`). The ancilla-free
implementation (direct permutation synthesis) keeps the total at exactly `exp_len + work_len` —
no ancilla register.

| N  | n_bits(N) | exp_len (t) | work_len (n) | Total qubits | Within the 25-qubit ceiling? |
|----|-----------|-------------|--------------|--------------|------------------------------|
| 15 | 4         | 4           | 4            | 8            | yes                          |
| 21 | 5         | 5           | 5            | 10           | yes                          |
| 35 | 6         | 6           | 6            | 12           | yes                          |
| 91 | 7         | 7           | 7            | 14           | yes                          |

(`qubit_budgets_all_targets` in `shor/tests/factor_kat.rs` checks every row.) The driver's ceiling
is the constant `MAX_QUBITS = 25`, which admits `N ≤ 4096`; a larger `N` that is not even and not a
prime power is reported as `Unsupported` before any circuit is built
(`factor_beyond_the_simulator_is_unsupported`: 4097 needs 13 + 13 = 26 qubits). The four targets
use 8 to 14 qubits, well inside it. The exponent register uses `t = ⌈log₂ N⌉` bits, half of the
textbook `2^t ≥ N²` width. Whether that gives enough phase resolution for continued-fraction
recovery is not proved: it is an empirical observation that the four targets each recover their
order at a recorded seed, nothing more.

### End-to-end factoring results

The factoring KATs run at the seeds below and complete in seconds on a laptop (debug,
unoptimized). The dominant cost is the controlled modular-exponentiation circuit: for N=91
(14 qubits) the state vector holds 2^14 amplitudes, and the permutation synthesis applies
O(N) transpositions per controlled-mult-mod stage.

| N  | Factors found | Seed | Total qubits | Evidence `(a, seed, s, r)` | Notes |
|----|---------------|------|--------------|----------------------------|-------|
| 15 | 3 × 5         | 0    | 8            | (2, 0, 4, 4)               | Textbook example |
| 21 | 3 × 7         | 0    | 10           | (2, 0, 27, 6)              | |
| 35 | 5 × 7         | 7    | 12           | (2, 7, 59, 12)             | s = 59 ≈ 11·64/12 recovers ord₂(35) = 12; seed 8 does not recover an order on its first run |
| 91 | 7 × 13        | 1    | 14           | (2, 1, 53, 12)             | Largest fixture; s = 53 ≈ 5·128/12 recovers ord₂(91) = 12 |

Each row is asserted by its own test (`factor_15_seed0`, `factor_21_seed0`, `factor_35_seed7`,
`factor_91_seed1` in `shor/tests/factor_kat.rs`), which rechecks the evidence from its fields:
`r` is the least order of `a`, is even, and the factor is one of `gcd(a^(r/2) ∓ 1, N)`. The seeds
are ones that succeed on the first circuit run, so the table records fixtures, not a success rate.

### Science↔engineering note

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The order-finding
circuit demonstrates Shor's *factoring mathematics* correctly at toy scale: the QFT extracts
the period of `x ↦ aˣ mod N` in polynomial time (in the number of qubits), the continued-
fraction recovery converts the measured phase to the order `r`, and the even-order factor
extraction `gcd(a^(r/2) ± 1, N)` produces the nontrivial factors. The algorithm's logic is
fully exhibited on the classical simulator.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 14
(N=91) the array holds 2^14 = 16 384 entries; at n = 25 (the ceiling) it requires ≈ 512 MiB. To
factor RSA-2048 (a 2048-bit N), this layout (`exp_len + work_len ≈ 2·2048`) would need about
4100 qubits — arithmetic on the layout, not a published estimate — requiring a `2^4100`-entry
array, far beyond any classical computer. This is the same posture as the simulator ceiling and
the index-calculus "asymptotic win not observable at toy scale" annotation: the simulator
exhibits the algorithm's logic, not its quantum speedup, which requires real quantum hardware
out of scope by construction.

N=91 is the largest fixture, not the largest admitted modulus: with this layout the 25-qubit
ceiling (`MAX_QUBITS`) admits `N ≤ 4096`, and a modulus needing more is reported as `Unsupported`
(4097 needs 13 + 13 = 26 qubits). Only the four fixtures are run end to end. The boundary is an
engineering annotation, not a mathematical omission.

### Test coverage (Shor factoring)

| Test file | Scope |
|-----------|-------|
| `shor/tests/modexp_kat.rs` | Permutation correctness (controlled-add-mod, controlled-mult-mod), modular exponentiation correctness, reversibility (forward + inverse = identity), ancilla-clean (full state restored), control-off no-op |
| `shor/tests/factor_kat.rs` | Order KATs (ord₂(15)=4, ord₇(15)=4, ord₂(21)=6, ord₂(35)=12, ord₂(91)=12); continued-fraction KATs (known phase → known order, the least order for every phase and coprime base, deterministic classical); end-to-end factoring (15→3×5 and 21→3×7 at seed 0, 35→5×7 at seed 7, 91→7×13 at seed 1, 14 qubits) with the evidence rechecked; the `Even`, `PrimePower`, `NotComposite` and `Unsupported` outcomes |

## Shor ECDLP — two-register period-finding on the classical simulator

This section covers Shor's ECDLP algorithm running end-to-end on the classical state-vector
simulator: the two-register period-finding circuit, the classical dual-lattice discrete-log
extraction, and the `solve_ecdlp` driver. The section records the two-register qubit budget
against the simulator ceiling and the resource-scale annotation. It records no success
probability or timing of its own; the quantities that matter are named below with the test or
function that computes them.

### Two-register qubit-budget-vs-curve-order table

The ECDLP circuit for a curve of prime order `r` uses two exponent registers of `t = ⌈log₂ r⌉`
bits each, plus the point register (x, y, λ scratch). The total qubit count is:

```text
total = 2·t + 3·coord_bits
```

where `coord_bits = ⌈log₂ p⌉` (bits to hold a field element in `GF(p)`).

Only the first row is a circuit this project builds; its 17 qubits are asserted by
`qubit_budget_within_ceiling` (`EXP_BITS == 4`, `TOTAL_QUBITS == 17`). The other rows are
arithmetic on the formula above for hypothetical prime orders and fields; the project does not
construct curves of those orders, and no test checks those rows.

| Curve order r | t = ⌈log₂ r⌉ | coord_bits = ⌈log₂ p⌉ | a-reg | b-reg | x | y | λ | Total | Within ~25-qubit ceiling? |
|---------------|--------------|----------------------|-------|-------|---|---|---|-------|---------------------------|
| 13 (p=7)      | 4            | 3                    | 4     | 4     | 3 | 3 | 3 | 17    | ✓ yes                     |
| 17 (p=11)     | 5            | 4                    | 5     | 5     | 4 | 4 | 4 | 22    | ✓ yes                     |
| 23 (p=17)     | 5            | 5                    | 5     | 5     | 5 | 5 | 5 | 25    | ✓ yes (ceiling-stress)    |
| 29 (p=23)     | 5            | 5                    | 5     | 5     | 5 | 5 | 5 | 25    | ✓ yes (ceiling-stress)    |
| 37 (p=31)     | 6            | 5                    | 6     | 6     | 5 | 5 | 5 | 27    | ✗ exceeds ceiling         |

The toy curve (r=13, p=7) uses 17 qubits — well within the ~25-qubit ceiling. By the formula, the
ceiling-stress case is r≈23–29 (p≈17–23), which reaches exactly 25 qubits. Larger prime-order
curves exceed the ceiling.

### End-to-end ECDLP results

The fixture instances below are solved by `solve_ecdlp(G, k·G, 0)`; each is asserted by its own
test, which checks the reported `k`, that the first attempt succeeds, the measured pair `(c, d)`
shown, and that `extract_k_from_measurement(c, d, r) == Some(k)`. No wall-clock time is recorded
here: the cost is dominated by simulating the 17-qubit circuit once per attempt (the test file's
ignore reasons give roughly 10 s per simulation in a debug build).

| k  | Q = k·G | Seed | Measured (c, d) | Total qubits | Test | Notes |
|----|---------|------|-----------------|--------------|------|-------|
| 3  | (2,2)   | 0    | (10, 13)        | 17           | `ecdlp_k3_seed0`  | Fixture instance 1 |
| 5  | (3,3)   | 0    | (12, 13)        | 17           | `ecdlp_k5_seed0`  | Fixture instance 2 |
| 7  | (5,4)   | 0    | (9, 12)         | 17           | `ecdlp_k7_seed0`  | Fixture instance 3 |
| 9  | (4,2)   | 0    | (8, 13)         | 17           | `ecdlp_k9_seed0`  | Fixture instance 4 |
| 12 | (1,5)   | 0    | (4, 12)         | 17           | `ecdlp_k12_seed0` | Boundary case (r−1) |

Seed 0 was kept for these fixtures because it succeeds on the first attempt in all five cases;
the table is a record of fixtures, not a success rate.

**Success probability (estimate).** One measurement yields the true `k` with probability
`primary_success_probability(G, Q)`, which is computed exactly from the circuit amplitudes. No
test asserts its numeric value (`primary_success_probability_positive_k3_k7` and
`primary_success_probability_positive_all_k` assert only that it is positive). A single run of
the ignored test printed values of about 0.64–0.66 for `k = 2..12`, so the budget of
`DEFAULT_MAX_ATTEMPTS = 8` measurements is exhausted with probability roughly `0.36⁸ ≈ 3·10⁻⁴`
(treating attempts as independent). Both figures are estimates from that printout. The test
`ecdlp_all_k_values_within_default_budget` (ignored for cost) asserts that every `k` in `2..=12`
is solved at seed 0 within the default budget.

### Science↔engineering note

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The two-register
ECDLP circuit demonstrates Shor's *ECDLP mathematics* correctly at toy scale: the two-register
iQFT concentrates the exponent pair near the dual lattice `l ≡ k·j (mod r)` of the hidden
subgroup of `(a, b) ↦ a·G + b·Q` (the circuit's exact distribution equals the Fourier-sum formula,
`distribution_matches_fourier_formula_k3` and `_k7`), and the classical extraction
`k = l·j⁻¹ (mod r)` from the nearest multiples `j/r`, `l/r` of the measured `(c, d)/M` recovers
the discrete logarithm. The algorithm's logic is exhibited at this scale; nothing here claims
the single-rounding extraction rule works for general `r`.

The ceiling is purely engineering: the `2^n`-amplitude array is the resource wall. At n = 17
(the toy curve) the array holds 2^17 ≈ 131K entries (≈ 2 MiB); at n = 25 it requires ≈ 512 MiB.
For a cryptographic curve (e.g., secp256k1, r ≈ 2^256), this project's layout
`2·t + 3·coord_bits` gives `2·256 + 3·256 = 1280` qubits — arithmetic on the layout, not a
published resource estimate — requiring a `2^1280`-entry array, far beyond any classical
computer. This is the same posture as the simulator ceiling and the factoring ceiling: the
simulator exhibits the algorithm's logic, not its quantum speedup, which requires real quantum
hardware out of scope by construction.

`solve_ecdlp` has no classical fallback. It accepts a candidate only after `k·G == Q`, and a
failed measurement is repeated with a fresh seed (the same randomized mechanism, which is the
algorithm itself). Exhausting the attempt budget is reported as `AttemptsExhausted`, a resource
outcome rather than a statement about the instance; no exponent is ever enumerated, and no
success reported by the solver comes from a path other than a measured pair.

### Test coverage (Shor ECDLP)

| Test file | Scope |
|-----------|-------|
| `shor/tests/ecdlp_kat.rs` | Order/group-law KATs (r·G=∞, scalar-multiple table); dual-lattice extraction KATs (exact peak pairs → k for every k and j≠0, `None` iff j≡0, out-of-range inputs rejected); exact outcome distribution equals the Fourier-sum formula (k=3, 7; all k ignored for cost); positive one-candidate success probability (k=3, 7; all k ignored); end-to-end ECDLP (k=3,5,7,9,12 with seed=0 and `Measurement` evidence; trivial cases; invalid generator; exhausted and zero budgets; retry at a fresh seed); qubit budget (17 qubits, within ceiling); in-workspace `rho` cross-checks (`#[ignore]`-gated: the cheap ones and one running a full Shor solve) |
