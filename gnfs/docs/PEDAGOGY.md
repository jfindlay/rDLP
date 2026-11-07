# Polynomial Selection: A Code-Tour Chapter

This chapter explains the `gnfs/src/polyselect` module — the polynomial-selection stage of the
General Number Field Sieve. It is organised by the implementation, not by mathematical abstraction:
each section introduces the mathematics, then shows how it is realised in code. A reader who has
read the `shared/numfield` PEDAGOGY chapter can read this one without consulting the source.

The chapter covers the base-m construction, degree selection, Murphy-E scoring, the root sieve,
Coppersmith multi-poly, and the public contracts that downstream stages consume.

---

## 1. Introduction: Why Polynomial Selection Matters

The General Number Field Sieve factors a large integer N by finding many pairs (a, b) ∈ ℤ² such
that both the rational norm G(a, b) = a − bm and the algebraic norm F(a, b) = b^d · f(a/b) are
B-smooth — divisible only by primes up to a bound B. When enough such pairs are collected, linear
algebra over 𝔽₂ produces a congruence of squares mod N, and hence a non-trivial factor.

The key insight is that the smoothness probability of F(a, b) depends heavily on the polynomial f.
A polynomial with small coefficients, many real roots, and many roots modulo small primes produces
smaller norms on average, which are more likely to be smooth. Polynomial selection is the stage
that finds such an f before the sieve begins.

The output of polynomial selection is a **polynomial pair** (f, g): the algebraic-side polynomial
f ∈ ℤ[x] of degree d ≥ 2, and the rational-side polynomial g(x) = x − m, where m is a shared
root mod N:

```text
f(m) ≡ 0 (mod N)
g(m) = 0,  g(x) = x − m
```

This shared-root invariant is the load-bearing contract that the sieve consumes. The sieve does not
care how f was found — only that f(m) ≡ 0 (mod N) and that f has small norms over the sieve region.

---

## 2. The `PolyPair` Contract (C-PolyPair)

The central type in this module is `PolyPair`, defined in `gnfs/src/polyselect/mod.rs`. It carries
the polynomial pair and all metadata that downstream stages need:

```rust
pub struct PolyPair {
    /// Algebraic-side polynomial f ∈ ℤ[x]. Generally non-monic for base-m; stored as-is.
    pub f: IntPoly,
    /// Rational-side polynomial g = x − m ∈ ℤ[x].
    pub g: IntPoly,
    /// The shared root: f(m) ≡ 0 (mod n) and g(m) = 0.
    pub m: BigInt,
    /// The integer to factor.
    pub n: BigInt,
    /// Polynomial degree (redundant with f.degree() but convenient for pattern matching).
    pub degree: usize,
    /// Skew parameter s: balances algebraic and rational norm sizes.
    /// None until Murphy-E scoring (G.B.2) computes it.
    pub skew: Option<f64>,
    /// Factor-base bounds (rational_bound, algebraic_bound).
    /// None until sieving (G.C) sets them.
    pub factor_base_bounds: Option<(u64, u64)>,
}
```

The invariants are checked by `PolyPair::verify()`:

1. `f` has degree ≥ 1 (not zero or constant).
2. `f.degree() == Some(self.degree)` — the `degree` field is consistent with the polynomial.
3. `g = x − m` — the rational polynomial is always the linear form `x − m`.
4. `f(m) ≡ 0 (mod n)` — the algebraic polynomial has m as a root mod n.

```rust
let pair = select_base_m(&n, 3);
assert_eq!(pair.verify(), Ok(()));
```

The `skew` and `factor_base_bounds` fields are `None` at construction and populated by later
pipeline stages. This is the substrate-over-specifies pattern: adding them now is cheaper than
amending the contract after G.C consumes it.

### The non-monic seam

Base-m expansion (§3) produces a polynomial f whose leading coefficient a_d satisfies 1 ≤ a_d < m,
so f is generally non-monic. The `shared-numfield` crate's `NumberField::new` requires a monic
defining polynomial. `PolyPair` resolves this seam with two methods:

```rust
// Return the monic form of f via the standard NFS homogenisation:
// f_monic(x) = a_d^{d−1} · f(x / a_d)
pub fn monic_f(&self) -> IntPoly;

// Construct the number field K = ℚ(α) where α is a root of the monic form of f.
pub fn number_field(&self) -> NumberField;
```

The homogenisation `f(x) → a_d^{d−1} f(x/a_d)` is monic and has roots a_d · αᵢ where αᵢ are the
roots of f. The coefficient of x^k in f_monic is a_k · a_d^{d−1−k}, which is always an integer.

The design decision is to store the original non-monic f and expose the monic form on demand:

- **Sieving uses the original f.** Norm computation in G.C uses the original polynomial
  coefficients; the monic transformation would change the norms.
- **The transformation is internal to `number_field()`.** Consumers that need `NumberField`
  element arithmetic call `poly_pair.number_field()`; consumers that need the original f
  (sieving, scoring) access `poly_pair.f` directly.

---

## 3. Base-m Construction

### The algorithm

The simplest NFS polynomial generator is **base-m expansion**, implemented in
`gnfs/src/polyselect/base_m.rs`. Given N and degree d:

1. Choose m = ⌊N^{1/(d+1)}⌋ + 1 — the smallest integer with m^{d+1} > N.
2. Write N in base m: N = a_0 + a_1·m + a_2·m² + … + a_d·m^d.
3. Define f(x) = Σ aᵢ xⁱ.

The key invariant is immediate: f(m) = N ≡ 0 (mod N). The rational side is g(x) = x − m.

```rust
pub fn select_base_m(n: &BigInt, degree: usize) -> PolyPair {
    let m = base_m_for_degree(n, degree);
    select_base_m_with_m(n, &m, degree)
}
```

The `base_m_for_degree` function computes m using Newton's method in `BigInt` arithmetic to find
⌊N^{1/(d+1)}⌋, then adds 1. The `base_m_digits` function performs the base-m expansion by
repeated division:

```rust
fn base_m_digits(n: &BigInt, m: &BigInt, d: usize) -> Vec<BigInt> {
    let mut coeffs = Vec::with_capacity(d + 1);
    let mut remainder = n.clone();
    for _ in 0..=d {
        let digit = &remainder % m;
        coeffs.push(digit.clone());
        remainder = (remainder - digit) / m;
    }
    coeffs
}
```

The resulting coefficient vector has exactly d + 1 entries (indices 0 through d), with the
constant term first (least-significant coefficient first, matching `IntPoly`'s storage convention).

### Why f(m) = N exactly

The base-m construction gives f(m) = N exactly, not just f(m) ≡ 0 (mod N). This is stronger than
required — the sieve only needs the congruence — but it is the natural consequence of writing N in
base m. The KAT confirms this:

```rust
let pair = select_base_m(&n, 3);
assert_eq!(pair.f.eval(&pair.m), n);  // exact equality, not just mod n
```

### The non-monic leading coefficient

The leading digit a_d satisfies 1 ≤ a_d < m (since N < m^{d+1} by construction of m). So f is
non-monic in general. For example, for N = 1022117 = 1009 × 1013 and d = 3, the base-m expansion
gives m = 101 and f(x) = c₀ + c₁x + c₂x² + c₃x³ with c₃ < 101. The `monic_f()` method
(§2) handles the conversion when needed.

---

## 4. Degree Selection

### The L-notation balance

The optimal degree d for NFS polynomial selection is determined by balancing the smoothness bounds
on the algebraic and rational sides. The classical heuristic derives from the L-notation complexity
of NFS:

> The NFS running time is L_N[1/3, c] = exp((c + o(1)) · (ln N)^{1/3} · (ln ln N)^{2/3})

Balancing the algebraic and rational sides gives the condition:

> d ≈ (3 ln N / ln ln N)^{1/3}

This is implemented in `optimal_degree`:

```rust
pub fn optimal_degree(n: &BigInt) -> usize {
    let bits = n.bits() as f64;
    let ln_n = bits * std::f64::consts::LN_2;
    let ln_ln_n = ln_n.ln();
    let d_float = (3.0 * ln_n / ln_ln_n).cbrt();
    let d = d_float.round() as usize;
    d.clamp(3, 6)
}
```

The bit length approximates ln N ≈ bits · ln 2. The result is clamped to [3, 6] for toy-scale N.

### Why degree 5 for RSA-100

For RSA-100 (330-bit, 100 decimal digits), the heuristic gives d ≈ 5. The KAT confirms this:

```rust
let n_rsa100 = BigInt::from_str(RSA_100_STR).unwrap();
assert_eq!(optimal_degree(&n_rsa100), 5);
```

For toy N (60–100 bit), the heuristic gives d = 3 or 4. The clamping to [3, 6] ensures that
even for very small N (where the formula might suggest d = 2), the degree is at least 3 — the
minimum that makes the NFS algebraic side non-trivial.

### The degree-smoothness tradeoff

Higher degree means:
- Larger algebraic norms (the coefficients of f grow as N^{1/d}, so larger d means smaller
  coefficients, but the norm F(a, b) = b^d · f(a/b) grows with d for fixed sieve region).
- More roots mod small primes (a degree-d polynomial can have up to d roots mod p, giving more
  factor-base elements).

The L-notation balance is the point where these effects cancel. At toy scale, the balance is
approximate and the degree choice matters less than at cryptographic scale.

---

## 5. Murphy-E Scoring

### What Murphy-E predicts

Murphy-E is the standard heuristic score for NFS polynomial pairs, implemented in
`gnfs/src/polyselect/murphy.rs`. It predicts the density of smooth relations produced by the
sieve. A higher score means more smooth pairs (a, b) in the sieve region, which means fewer
sieve pairs needed to collect enough relations for the linear algebra step.

The score is defined as the average over a sample of (a, b) pairs in the sieve region of the
product of Dickman-ρ values for the algebraic and rational norms:

```text
E(f, g) ≈ (1/|S|) Σ_{(a,b) ∈ S} ρ(log|F(a,b)| / log B_f) · ρ(log|G(a,b)| / log B_g)
```

where:
- F(a, b) = b^d · f(a/b) is the algebraic norm (homogeneous form of f).
- G(a, b) = a − bm is the rational norm.
- B_f, B_g are the algebraic and rational factor-base bounds.
- ρ is the Dickman rho function.

The implementation samples a 50×50 grid over the sieve region [-M, M] × [1, M] with M = 1000:

```rust
pub fn score(pair: &PolyPair) -> f64 {
    const M: f64 = 1000.0;
    const GRID: usize = 50;
    let log_bf = (1_000_000.0_f64).ln();
    let log_bg = (1_000_000.0_f64).ln();
    // ... sample grid, compute rho products, average ...
}
```

### The Dickman ρ function

The Dickman rho function ρ(u) is the probability that a random integer near x is x^{1/u}-smooth
(has no prime factor larger than x^{1/u}). It satisfies:

- ρ(u) = 1 for u ≤ 1.
- ρ(u) = 1 − ln u for 1 < u ≤ 2.
- For u > 2: the delay-differential equation u·ρ'(u) = −ρ(u−1).

The implementation uses the exact formula for u ≤ 2 and 4th-order Runge-Kutta numerical
integration of the recurrence for 2 < u ≤ 25, with ρ(u) = 0 for u > 25:

```rust
pub fn dickman_rho(u: f64) -> f64 {
    if u <= 1.0 { return 1.0; }
    if u <= 2.0 { return 1.0 - u.ln(); }
    if u > 25.0 { return 0.0; }
    dickman_rho_numerical(u)
}
```

The numerical integrator builds a table of 2301 entries (step h = 0.01) from u = 2 to u = 25,
using the one-unit lookback ρ(t−1) at each step. The full table is allocated as a `Vec` to avoid
the correctness hazard of a rolling window that would overwrite values still needed for the
lookback.

Key known values (verified by KATs):
- ρ(1.5) = 1 − ln(1.5) ≈ 0.5945
- ρ(2.0) = 1 − ln(2) ≈ 0.3069
- ρ(3.0) ≈ 0.04861 (from tables)
- ρ(10.0) ≈ 2.77 × 10⁻¹⁰ (essentially zero)

### The skew parameter

The optimal skew s balances algebraic and rational norm sizes. The score function applies the
skew by sampling (a, b·s) instead of (a, b), which rescales the algebraic norm relative to the
rational norm. The skew stored in `pair.skew` is used if present; otherwise it defaults to 1.0.

### Science↔engineering disconnect at toy scale (principle-4 annotation)

**Murphy-E's predictive value — that higher E implies more relations — only manifests at sieve
scale (N ≳ 2^100, sieve region ≳ 10^6 pairs).** At toy scale (N < 2^60), the sieve region is
tiny and the smoothness probabilities are dominated by the factor-base bounds rather than the
polynomial shape. At toy scale, Murphy-E is a *ranking heuristic* whose payoff is under-exposed:
the ordering it induces is correct in expectation, but the absolute values are not meaningful.

This is why the KATs test ordering and self-consistency rather than absolute values:

```rust
// KAT 1: smaller coefficients → higher score (ordering is correct)
assert!(score_better > score_worse);

// KAT 2: scaling f by 10 decreases the score (monotonicity)
assert!(score_1x > score_10x);

// KAT 3: score is always positive (positivity)
assert!(score(&pair) > 0.0);
```

Downstream consumers (root sieve, Coppersmith, D.A NFS-DL poly selection) should treat the score
as an ordinal ranking, not a cardinal measure.

---

## 6. Root Sieve

### The Kleinjung rotation

The root sieve, implemented in `gnfs/src/polyselect/root_sieve.rs`, searches for a better
polynomial by applying **Kleinjung rotation** to a seed pair. The rotation is:

```text
f'(x) = f(x) + (j·x + k)·g(x)
g'(x) = g(x)   (unchanged)
```

This preserves the shared root m mod N because g(m) = 0:

```text
f'(m) = f(m) + (j·m + k)·g(m) = 0 + (j·m + k)·0 = 0  (mod N)
```

So `PolyPair::verify()` holds for every rotated pair. The rotation is implemented in `rotate`:

```rust
pub fn rotate(seed: &PolyPair, j: i64, k: i64) -> PolyPair {
    let m = &seed.m;
    let mut coeffs: Vec<BigInt> = seed.f.coeffs.clone();
    while coeffs.len() < 3 { coeffs.push(BigInt::from(0i64)); }

    let j_big = BigInt::from(j);
    let k_big = BigInt::from(k);

    // (j·x + k)·(x − m) = j·x² + (k − j·m)·x − k·m
    coeffs[0] -= &k_big * m;       // f'[0] += −k·m
    coeffs[1] += &k_big - &j_big * m;  // f'[1] += k − j·m
    coeffs[2] += &j_big;           // f'[2] += j

    let f_prime = IntPoly::from_coeffs(coeffs);
    let g = IntPoly::from_coeffs(vec![-m.clone(), BigInt::one()]);
    PolyPair::new(f_prime, g, m.clone(), n.clone())
}
```

For g(x) = x − m, the rotation adds a degree-2 perturbation to f. Since base-m always gives
d ≥ 3, the degree of f is unchanged by the rotation (the j contribution lands at index 2, below
the leading term at index d ≥ 3).

### The search strategy

`root_sieve` searches over a grid of (j, k) values, scoring each candidate with Murphy-E and
returning the best:

```rust
pub fn root_sieve(seed: &PolyPair, config: &RootSieveConfig) -> PolyPair {
    let mut best_pair = rotate(seed, 0, 0);  // seed itself
    let mut best_score = score(&best_pair);

    for j in -config.j_range..=config.j_range {
        for k in -config.k_range..=config.k_range {
            if j == 0 && k == 0 { continue; }
            let candidate = rotate(seed, j, k);
            let s = score(&candidate);
            if s > best_score {
                best_score = s;
                best_pair = candidate;
            }
        }
    }
    best_pair
}
```

The default `RootSieveConfig` uses `j_range = k_range = 10`, giving 441 candidates. The seed
itself (at j = 0, k = 0) is always included, so the returned pair is always at least as good as
the seed. The search is deterministic: the grid is traversed in row-major order and ties are
broken by keeping the first maximum found.

### The `RootSieveGenerator` trait implementation

`RootSieveGenerator` implements the `PolyGenerator` trait, yielding all candidates in the grid
rather than just the best:

```rust
pub struct RootSieveGenerator {
    pub seed: PolyPair,
    pub config: RootSieveConfig,
}

impl PolyGenerator for RootSieveGenerator {
    fn generate(&self) -> impl Iterator<Item = PolyPair> {
        // Yields all (2·j_range + 1) × (2·k_range + 1) candidates in row-major order.
        // ...
    }
}
```

This fits into the common score-and-rank pipeline alongside `BaseMGenerator` and
`CoppersmithGenerator`. Callers can use `.max_by(...)` to recover the best, or `.take(limit)` for
early termination.

### Why rotation preserves the root: the algebraic argument

The root-preservation property is not just a numerical coincidence — it is an algebraic identity.
For any polynomial h(x) ∈ ℤ[x]:

```text
f'(m) = f(m) + h(m)·g(m) = 0 + h(m)·0 = 0  (mod N)
```

The Kleinjung rotation uses h(x) = j·x + k, but any h would work. The rotation is a change of
basis in the space of polynomials with root m mod N: the set {f + h·g : h ∈ ℤ[x]} is the coset
of g in ℤ[x] that contains f. The root sieve searches this coset for the element with the best
Murphy-E score.

### Science↔engineering note

At toy scale (N < 2^60), the Murphy-E improvement from rotation is small and may not manifest as
a strict improvement for every seed. The KAT uses `≥` rather than `>`:

```rust
// KAT 1: root_sieve returns score ≥ seed score (not necessarily strictly better)
assert!(best_score >= seed_score);
```

At cryptographic scale (RSA-768+), rotation search is essential and produces measurable
improvements in relation yield.

---

## 7. Coppersmith Multi-Poly

### The mathematical construction

Coppersmith's multiple-polynomial method, implemented in `gnfs/src/polyselect/coppersmith.rs`,
generates several algebraic-side polynomials f₁, f₂, …, f_k that all share the same rational-side
polynomial g = x − m. Each variant is produced by a Kleinjung rotation:

```text
f_i(x) = f_0(x) + (j_i · x + k_i) · g(x)
```

The rotation parameters (j_i, k_i) are chosen by a symmetric spiral schedule:

```text
i = 0: (j=0,  k=0)          — the seed itself
i = 1: (j=0,  k=+step)      — translate by +step
i = 2: (j=0,  k=-step)      — translate by -step
i = 3: (j=+step, k=0)       — rotate by +step
i = 4: (j=-step, k=0)       — rotate by -step
i = 5: (j=0,  k=+2*step)    — translate by +2*step
...
```

The spiral keeps the parameters small and symmetric, exploring the neighbourhood of the seed in
a balanced way. The implementation:

```rust
pub fn coppersmith_polys(seed: &PolyPair, config: &CoppersmithConfig) -> Vec<PolyPair> {
    (0..config.num_polys)
        .map(|i| {
            let (j, k) = rotation_params(i, config.step);
            rotate(seed, j, k)
        })
        .collect()
}

pub fn coppersmith_best(seed: &PolyPair, config: &CoppersmithConfig) -> PolyPair {
    coppersmith_polys(seed, config)
        .into_iter()
        .max_by(|a, b| score(a).partial_cmp(&score(b)).expect("Murphy-E score is finite"))
        .expect("coppersmith_polys returns at least one pair when num_polys >= 1")
}
```

### Demonstration-fidelity role (principle-4 annotation)

**In production NFS (e.g., RSA-768), using different polynomials for different sieve regions
measurably improves the relation yield.** The idea is that different parts of the sieve region
may be better served by different polynomials — a polynomial that is good on average may be poor
in some sub-region, and a different rotation may be better there.

**At toy scale (60–100 bit N), the sieve region is too small for this effect to manifest.** The
Murphy-E improvement from multi-poly is typically small (< 2×), and the absolute score values are
not meaningful. This is a **demonstration-fidelity** implementation: the mathematical content is
complete, but the engineering payoff requires cryptographic-scale N to observe.

The KAT documents this explicitly:

```rust
// KAT 5: principle-4 annotation — improvement is < 2× at toy scale
let improvement = best_score / seed_score;
assert!(improvement >= 1.0);  // always at least as good as seed
assert!(improvement < 2.0);   // but not dramatically better at toy scale
```

This is not a failure condition — it is the expected behaviour at toy scale, and the assertion
documents the science↔engineering disconnect.

### The `CoppersmithGenerator` trait implementation

`CoppersmithGenerator` implements `PolyGenerator`, yielding all variants in spiral order:

```rust
pub struct CoppersmithGenerator {
    pub seed: PolyPair,
    pub config: CoppersmithConfig,
}
```

The default config generates 5 variants with step = 1. Callers can increase `num_polys` and
`step` to explore a larger neighbourhood.

---

## 8. The `PolyGenerator` Trait: A Common Pipeline

All three generators — base-m, root sieve, and Coppersmith — implement the `PolyGenerator` trait:

```rust
pub trait PolyGenerator {
    fn generate(&self) -> impl Iterator<Item = PolyPair>;
}
```

This trait is the common interface for the score-and-rank pipeline. The generator produces
candidates; the scorer (`murphy::score`) ranks them. The three generators differ in how many
candidates they produce:

| Generator | Candidates | Use case |
|-----------|-----------|----------|
| `BaseMGenerator` | 1 (exactly) | Starting point; always valid |
| `RootSieveGenerator` | (2·j+1)×(2·k+1) | Grid search for better polynomial |
| `CoppersmithGenerator` | `num_polys` | Multi-poly set for sieve region partitioning |

A typical pipeline:

```rust
// 1. Start with base-m.
let seed = select_base_m(&n, optimal_degree(&n));

// 2. Improve with root sieve.
let config = RootSieveConfig::default();
let best = root_sieve(&seed, &config);

// 3. Optionally generate a multi-poly set.
let copper_config = CoppersmithConfig::default();
let polys = coppersmith_polys(&best, &copper_config);
```

---

## 9. C-PolyPair and C-Score: Public Contracts

### C-PolyPair

The `PolyPair` type is the **C-PolyPair contract** — the interface between polynomial selection
and all downstream stages. It is frozen at G.B.1 and consumed by:

- **G.B.2 (Murphy-E scoring):** reads `pair.f`, `pair.m`, `pair.skew`; writes `pair.skew`.
- **G.B.3 (root sieve):** reads `pair.f`, `pair.g`, `pair.m`, `pair.n`; produces new `PolyPair`.
- **G.B.4 (Coppersmith):** same as root sieve.
- **G.C (sieving):** reads `pair.f`, `pair.m`, `pair.n`, `pair.degree`; writes
  `pair.factor_base_bounds`. The sieve evaluates the algebraic norm F(a, b) = b^d · f(a/b) for
  each sieve pair (a, b) using the original non-monic f.
- **D.A (NFS-DL poly selection):** adapts `PolyPair` for the discrete logarithm setting. The
  polynomial-pair structure is identical; the difference is in the linear algebra step.

The `verify()` method is the contract enforcement mechanism: any code that produces a `PolyPair`
should call `verify()` before passing it downstream. The KATs enforce this:

```rust
// Every rotation must satisfy verify().
for j in -3i64..=3 {
    for k in -3i64..=3 {
        let rotated = rotate(&seed, j, k);
        assert_eq!(rotated.verify(), Ok(()));
    }
}
```

### C-Score

The `score` function is the **C-Score contract** — the Murphy-E scoring interface consumed by:

- **G.B.3 (root sieve):** scores each rotation candidate to find the best.
- **G.B.4 (Coppersmith):** scores each variant to select `coppersmith_best`.
- **D.A (NFS-DL):** adapts the scoring for the DL polynomial selection problem.
- **E.K (factor-base balancing):** uses the score to balance algebraic and rational factor-base
  sizes. The Murphy-E score's sensitivity to the skew parameter is the signal that E.K uses to
  choose the optimal factor-base bounds.

The C-Score interface is:

```rust
pub fn score(pair: &PolyPair) -> f64;
```

It takes `&PolyPair` (not `&mut`), so it is safe to call from multiple contexts without
mutation. The score is deterministic for a fixed pair (no hidden mutable state).

### Cross-track C3 reuse (D.A and E.K)

The C3 cross-track reuse is the observation that the same polynomial-selection infrastructure
serves multiple downstream tracks:

**D.A (NFS-DL poly selection):** The discrete logarithm variant of NFS uses the same polynomial
pair structure. The base-m construction, Murphy-E scoring, and root sieve are all reused directly.
The only difference is the target: instead of factoring N, the goal is to compute discrete
logarithms in a finite field, which requires a polynomial pair adapted to the field's structure.

**E.K (factor-base balancing):** The factor-base balancing stage uses the Murphy-E score's
sensitivity to the skew parameter to choose the optimal split between algebraic and rational
factor-base sizes. The `pair.skew` field (set by G.B.2 scoring) is the signal that E.K reads.
The `pair.factor_base_bounds` field (set by G.C sieving) is the output that E.K produces.

Both reuses are enabled by the `PolyPair` struct's over-specification: the `skew` and
`factor_base_bounds` fields are `None` at construction and populated by later stages, so the
same type serves all three tracks without modification.

---

## 10. KAT Summary

The following table lists the key known-answer tests across the module, with the mathematical
fact each one verifies. All tests are in `gnfs/tests/`.

| Test | File | Mathematical fact verified |
|------|------|---------------------------|
| `kat1_toy_base_m_round_trip` | `base_m_kat.rs` | f(m) = N exactly for toy N; g(m) = 0; verify() holds |
| `kat2_rsa100_base_m_deterministic` | `base_m_kat.rs` | Same (N, d) → same f; f(m) = N for RSA-100 |
| `kat3_optimal_degree` | `base_m_kat.rs` | d = 3–4 for toy N; d = 5 for RSA-100 (330-bit) |
| `monic_f_and_number_field_from_base_m` | `base_m_kat.rs` | monic_f() is monic; number_field() does not panic |
| `kat1_ordering_better_scores_higher` | `murphy_kat.rs` | Smaller coefficients → higher Murphy-E score |
| `kat2_monotonicity_scaling_decreases_score` | `murphy_kat.rs` | Scaling f by 10 decreases score (monotonicity) |
| `kat3_positivity` | `murphy_kat.rs` | score(pair) > 0 for any valid PolyPair |
| `kat4a–4e_dickman_rho_*` | `murphy_kat.rs` | ρ(0.5) = 1; ρ(1.5) ≈ 0.5945; ρ(2) ≈ 0.3069; ρ(10) < 1e-6 |
| `kat5_score_is_deterministic` | `murphy_kat.rs` | score() is deterministic (no hidden mutable state) |
| `kat1_improvement_score_ge_seed` | `root_sieve_kat.rs` | root_sieve returns score ≥ seed score |
| `kat2_determinism` | `root_sieve_kat.rs` | root_sieve is deterministic for fixed seed and config |
| `kat3_verify` | `root_sieve_kat.rs` | root_sieve result satisfies PolyPair::verify() |
| `kat3b_all_candidates_verify` | `root_sieve_kat.rs` | All 121 candidates in 5×5 grid satisfy verify() |
| `kat4_generator_candidate_count` | `root_sieve_kat.rs` | Generator yields exactly (2j+1)×(2k+1) candidates |
| `kat5_generator_best_matches_root_sieve` | `root_sieve_kat.rs` | Generator best score agrees with root_sieve |
| `kat1_all_polys_verify` | `coppersmith_kat.rs` | All coppersmith_polys variants satisfy verify() |
| `kat2_exact_count` | `coppersmith_kat.rs` | coppersmith_polys returns exactly num_polys pairs |
| `kat2b_single_poly_is_seed` | `coppersmith_kat.rs` | num_polys=1 returns the seed (identity rotation) |
| `kat3_best_score_ge_seed` | `coppersmith_kat.rs` | coppersmith_best returns score ≥ seed score |
| `kat3c_best_matches_max_of_polys` | `coppersmith_kat.rs` | coppersmith_best agrees with max of coppersmith_polys |
| `kat5_principle4_improvement_small_at_toy_scale` | `coppersmith_kat.rs` | Multi-poly improvement < 2× at toy scale (principle-4) |

---

## 11. What Comes Next

Polynomial selection is the first stage of the GNFS pipeline. The output — a `PolyPair` with a
good Murphy-E score — is the input to the sieve.

**G.C (sieving)** consumes `C-PolyPair` directly. For each sieve pair (a, b) in the sieve region,
the sieve evaluates:
- The algebraic norm F(a, b) = b^d · f(a/b) using the original non-monic f.
- The rational norm G(a, b) = a − bm.

Both norms are trial-divided over their respective factor bases. A pair (a, b) is **smooth** if
both norms factor completely over their factor bases. The factor-base bounds (stored in
`pair.factor_base_bounds`) are set by G.C based on the polynomial's Murphy-E score and the
target number of relations.

The `number_field()` method on `PolyPair` is the seam between polynomial selection and the
algebraic side of the sieve: it constructs the number field K = ℚ(α) where α is a root of the
monic form of f, which the sieve uses for ideal-theoretic computations (Dedekind factorisation,
factor-base construction).

The polynomial pair is also the input to the linear algebra step (G.D): the smooth pairs collected
by the sieve form a matrix over 𝔽₂, and the null space of this matrix gives the congruence of
squares that factors N.

---

## Further Reading

1. **Murphy, B. (1999).** *Polynomial selection for the number field sieve integer factorisation
   algorithm.* PhD thesis, Australian National University. The original source for Murphy-E scoring,
   the Dickman-ρ integral, and the skew parameter.

2. **Kleinjung, T. (2006).** *Polynomial selection.* CADO workshop on integer factorization. The
   rotation `f' = f + (j·x + k)·g` and the root-sieve search strategy.

3. **Bai, S., Bouvier, C., Kruppa, A., Zimmermann, P. (2014).** *Better polynomials for GNFS.*
   Mathematics of Computation. The state of the art in polynomial selection, including the
   alpha-value (root property) term in Murphy-E and the size-optimisation step.

4. **Coppersmith, D. (1993).** *Modifications to the number field sieve.* Journal of Cryptology,
   6(3), 169–180. The original multi-polynomial method.

5. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective.*
   2nd ed. Springer. Chapter 6 covers the Number Field Sieve, including polynomial selection and
   the role of the factor base.

6. **Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) (1993).** *The Development of the Number Field
   Sieve.* Springer LNM 1554. The original papers on GNFS, including the algebraic side and the
   role of the defining polynomial.

---

# Sieving: A Code-Tour Chapter

This chapter explains the `gnfs/src/sieve` module — the sieving stage of the General Number Field
Sieve. It is organised by the implementation, not by mathematical abstraction: each section
introduces the mathematics, then shows how it is realised in code. A reader who has read the
polynomial selection chapter (§1–§11 above) can read this one without consulting the source.

The chapter covers the relation as the unit of NFS data, the two-sided factor base, the norm
bridge, the line sieve, the special-q strategy, lattice sieving, the exponent vector and relation
matrix, and the downstream reuse in G.D, G.E, and D.A.

---

## 12. The Relation as the Unit of NFS Data

### What a relation is

A **relation** in GNFS is a coprime integer pair (a, b) with b ≥ 1 for which both the rational
norm and the algebraic norm are smooth over their respective factor bases:

```text
N_rat(a, b) = a − b·m          (rational side: the value of g(x) = x − m at a/b, cleared)
N_alg(a, b) = b^d · f(a/b)     (algebraic side: the homogeneous form of f at (a, b))
```

Both norms must be smooth — divisible only by primes up to their respective bounds B_rat and B_alg.
A pair that is smooth on only one side is not a relation; it contributes nothing to the factoring
computation.

The `Relation` type in `gnfs/src/sieve/mod.rs` carries the pair and both smoothness witnesses:

```rust
pub struct Relation {
    /// The a-coordinate (can be negative).
    pub a: BigInt,
    /// The b-coordinate (always positive: b ≥ 1).
    pub b: BigInt,
    /// Exponent vector over the rational factor base.
    pub rational_exponents: ExponentVector,
    /// Exponent vector over the algebraic factor base.
    pub algebraic_exponents: ExponentVector,
    /// True if the rational norm a − b·m is negative.
    pub rational_sign: bool,
}
```

The `verify()` predicate checks all four invariants:

1. gcd(a, b) = 1 (coprimality — non-coprime pairs are degenerate).
2. The rational exponent vector reconstructs |N_rat(a, b)| (product of p^e matches the norm).
3. The algebraic exponent vector reconstructs |N_alg(a, b)| (same check on the algebraic side).
4. The sign matches: `rational_sign == (a − b·m < 0)`.

```rust
rel.verify(&poly, &fb).expect("relation must be valid");
```

### Why smoothness on both sides is required

The factoring computation (G.F) requires finding a subset S of relations such that the product
of all rational norms is a perfect square in ℤ and the product of all algebraic norms is a
perfect square in the number field K = ℚ(α). This is the congruence of squares that yields a
non-trivial factor of N.

For the product to be a perfect square, every prime must appear to an even total exponent across
all selected relations. This is a linear algebra problem over GF(2): find a subset S such that
the sum of all exponent vectors (mod 2) is the zero vector. The relation matrix has one row per
relation and one column per prime (plus sign and quadratic-character columns). The null space of
this matrix over GF(2) gives the subsets S.

A relation that is smooth on only one side has no algebraic (or rational) exponent vector, so it
cannot participate in the square-root computation. Both sides must be smooth for the relation to
contribute a complete row to the matrix.

---

## 13. The Two-Sided Factor Base

### The rational factor base

The rational factor base is the set of primes p ≤ B_rat. A rational norm N_rat(a, b) = a − b·m
is smooth over this base if all its prime factors are ≤ B_rat.

### The algebraic factor base: degree-1 prime ideals

The algebraic factor base is not simply the set of primes p ≤ B_alg. It is the set of
**degree-1 prime ideals** (p, α − r) in ℤ[α], where α is a root of f and r ∈ [0, p) satisfies
f(r) ≡ 0 (mod p). Each prime p contributes as many ideals as f has roots mod p (up to deg(f)).

The `AlgebraicPrime` type represents one such ideal:

```rust
pub struct AlgebraicPrime {
    /// The rational prime p.
    pub p: u64,
    /// The root r ∈ [0, p) with f(r) ≡ 0 (mod p).
    pub r: u64,
    /// Index of this ideal in the algebraic factor base (for column mapping).
    pub index: usize,
    /// True if p | disc(f) (bad prime).
    pub is_bad_prime: bool,
}
```

The key sieve condition: the ideal (p, α − r) divides N_alg(a, b) if and only if a ≡ r·b (mod p).
This is the algebraic analogue of the rational condition a ≡ b·m (mod p) for the rational side.

### Why degree-1 ideals

At toy scale, the degree-1 ideal (p, r) representation suffices because:

1. The algebraic norm N_alg(a, b) = b^d · f(a/b) factors in ℤ[α] as a product of prime ideals.
2. For a prime p with f(r) ≡ 0 (mod p), the ideal (p, α − r) divides the principal ideal (a − b·α)
   in ℤ[α], and hence divides N_alg(a, b) = Norm(a − b·α).
3. The sieve condition a ≡ r·b (mod p) is the efficient way to enumerate which (a, b) pairs have
   the ideal (p, r) dividing their algebraic norm.

### The `FactorBase` type

The `FactorBase` type in `gnfs/src/sieve/factor_base.rs` holds both sides:

```rust
pub struct FactorBase {
    /// Rational factor base: primes p ≤ B_rat, sorted ascending.
    pub rational_primes: Vec<u64>,
    /// Algebraic factor base: degree-1 prime ideals (p, r), sorted by (p, r).
    pub algebraic_ideals: Vec<AlgebraicPrime>,
    /// Rational smoothness bound B_rat.
    pub b_rat: u64,
    /// Algebraic smoothness bound B_alg.
    pub b_alg: u64,
    /// Number of obstruction columns reserved for G.E (sign + quadratic chars).
    pub obstruction_count: usize,
}
```

Construction:

```rust
let fb = FactorBase::new(&poly.f, b_rat, b_alg);
```

`FactorBase::new` uses `factor_base_up_to` (from `shared-numth`) for the rational side, and
brute-force root-finding (checking f(r) ≡ 0 mod p for each r ∈ 0..p) for the algebraic side.

### Bad primes and the principle-4 annotation

A prime p is **bad** if p | disc(f). At bad primes, Dedekind's theorem does not apply directly:
ℤ[α] may not be the full ring of integers at p, and the ideal factorisation above p may be more
complex than the root structure of f mod p suggests.

At toy scale, bad primes are prominent. For f(x) = x³ − x − 1, the discriminant is −23, so p = 23
is a bad prime. The implementation includes bad primes in the algebraic factor base with the
`is_bad_prime: true` flag, using direct root-finding (correct for linear factors even at bad
primes). The flag documents the principle-4 over-exposure:

**Principle-4 annotation.** At cryptographic scale (RSA-768+), bad primes are marginal: the
discriminant has few prime factors, and those primes are a negligible fraction of the algebraic
factor base. At toy scale, bad primes can be a significant fraction of the factor base and are
unavoidable for hand-picked polynomials. The `is_bad_prime` flag makes this visible without
changing the algorithm.

### Column indexing for G.D and G.E

The `AlgebraicPrime::index` field and the lookup methods `rational_index(p)` and
`algebraic_index(p, r)` provide the mapping from factor-base elements to matrix columns. G.D
(filtering) and G.E (linear algebra) use these indices to build and manipulate the relation
matrix without re-scanning the factor base.

The `matrix_width()` method returns the total number of columns:

```rust
pub fn matrix_width(&self) -> usize {
    self.rational_size() + self.algebraic_size() + self.obstruction_count
}
```

The `obstruction_count` is initialised to 1 (the sign/−1 column). G.E will add quadratic-character
columns; the slot exists so the matrix-width calculation is stable across stages.

---

## 14. The Norm Bridge

### Why norms are signed BigInt but trial_smooth takes Uint<4>

The rational norm N_rat(a, b) = a − b·m can be negative (when a < b·m). The algebraic norm
N_alg(a, b) = b^d · f(a/b) can also be negative depending on f and (a, b). Both are computed as
signed `BigInt` values.

The smoothness predicate `trial_smooth` (from `shared-numth`) operates on `Uint<4>` — unsigned
256-bit integers. The **norm bridge** in `gnfs/src/sieve/norms.rs` converts a signed norm to its
absolute value as `Uint<4>`:

```rust
pub fn norm_to_uint(norm: &BigInt) -> Result<Uint<4>, NormBridgeError> {
    let abs_norm = norm.abs();
    let bits = abs_norm.bits() as usize;
    if bits > 256 {
        return Err(NormBridgeError::Overflow { bits_required: bits });
    }
    // ... convert to Uint<4> via big-endian byte array ...
}
```

The bridge rejects out-of-range norms with `NormBridgeError::Overflow` rather than silently
truncating. Silent truncation would produce incorrect smoothness witnesses — a correctness hazard.

### The sign column in the relation matrix

The sign of the rational norm is tracked separately in the `rational_sign: bool` field of
`Relation`. This is the "−1 column" for G.E's linear algebra: the product of all selected
relations' rational norms must be a perfect square, which requires the sign product to be +1
(an even number of negative norms). The sign column encodes this constraint.

The algebraic norm's sign is not stored separately because the algebraic square root computation
(G.F) handles sign via the real embedding of K, not via a matrix column.

### The algebraic norm as a resultant

The algebraic norm N_alg(a, b) = b^d · f(a/b) is computed as the homogeneous form:

```text
N_alg(a, b) = Σ_{i=0}^{d} f.coeffs[i] · a^i · b^{d−i}
```

This equals Res(f, a − b·x) up to sign and leading-coefficient factors. The resultant
relationship is the algebraic hook: N_alg(a, b) is the norm of the ideal (a − b·α) in ℤ[α],
which equals the resultant of f and a − b·x. The implementation computes the homogeneous form
directly (avoiding rational arithmetic) rather than via the resultant, which would be slower.

---

## 15. The Line Sieve

### The log-p mark-then-confirm pattern

The line sieve in `gnfs/src/sieve/line.rs` is the baseline relation-collection algorithm. For
each b in 1..=B, it sieves the range a ∈ −A..=A on both sides simultaneously:

**Step 1.** Initialise a sieve array of size 2A + 1 with zeros.

**Step 2 (rational side).** For each prime p in the rational factor base, find the starting
a ≡ b·m (mod p) and mark all a in range with += log₂(p). The rational norm N_rat(a, b) = a − b·m
is divisible by p iff a ≡ b·m (mod p).

**Step 3 (algebraic side).** For each ideal (p, r) in the algebraic factor base, find the
starting a ≡ r·b (mod p) and mark all a in range with += log₂(p). The algebraic norm N_alg(a, b)
is divisible by the ideal (p, r) iff a ≡ r·b (mod p).

**Step 4 (threshold filter).** Collect candidates where the sieve value exceeds a smoothness
threshold (default: 0.8 × log₂(B_alg)). The threshold accepts candidates where the algebraic
norm has accumulated enough log contributions to plausibly be smooth.

**Step 5 (confirm).** For each candidate (a, b): check gcd(a, b) = 1; compute both norms; call
`trial_smooth` on both; if both are fully smooth, construct a `Relation`.

```rust
pub fn line_sieve(
    poly: &PolyPair,
    fb: &FactorBase,
    config: &LineSieveConfig,
) -> Vec<Relation> {
    // ...
    for b in 1u64..=config.b_bound {
        let mut sieve: Vec<f32> = vec![0.0f32; sieve_len];
        // Step 2: rational side sieve
        for (pi, &p) in fb.rational_primes.iter().enumerate() {
            let bm_mod_p = mod_u64_bigint(&(&b_big * &poly.m), p);
            let start_a = first_a_in_range(bm_mod_p, p, -a_bound);
            let mut a = start_a;
            while a <= a_bound { sieve[(a + a_bound) as usize] += rat_logs[pi]; a += p as i64; }
        }
        // Step 3: algebraic side sieve
        for (ai, ap) in fb.algebraic_ideals.iter().enumerate() {
            let rb_mod_p = (ap.r as u128 * b as u128 % ap.p as u128) as u64;
            let start_a = first_a_in_range(rb_mod_p, ap.p, -a_bound);
            let mut a = start_a;
            while a <= a_bound { sieve[(a + a_bound) as usize] += alg_logs[ai]; a += ap.p as i64; }
        }
        // Steps 4–5: threshold filter and confirm
        // ...
    }
}
```

### Why this is the engineering heart of NFS

The log-p mark-then-confirm pattern is the key engineering insight of NFS sieving. The sieve
array accumulates approximate log contributions cheaply (one addition per prime per sieve
position). The threshold filter eliminates most candidates before the expensive `trial_smooth`
call. At cryptographic scale (B ≈ 10⁷, A ≈ 10⁷, B_rat/B_alg ≈ 10⁶), the threshold filters
~99% of candidates, making the sieve dramatically faster than brute-force trial division of
every (a, b) pair.

### Principle-4 annotation: asymptotic win under-exposed at toy scale

**At toy scale (A = 10, B = 3, B_rat = B_alg = 30), the log-sieve barely beats direct trial
division.** The sieve array has only 21 entries; the factor bases have ~10 primes each; the
threshold filters almost nothing (every candidate passes). The asymptotic win — avoiding trial
division for most pairs — is not visible because there are so few pairs to begin with.

The KAT uses a conservative lower bound (≥ 5 relations) rather than an exact count, because the
threshold heuristic may miss some smooth pairs at toy scale. The exact count is pinned by the
determinism KAT (KAT b3 in `line_sieve_kat.rs`).

At cryptographic scale, the sieve array has ~10⁷ entries, the factor bases have ~10⁵ elements,
and the threshold filters ~99% of candidates. The log-sieve is then 100–1000× faster than
brute-force trial division.

---

## 16. The Special-q Strategy

### The yield multiplier

The special-q strategy in `gnfs/src/sieve/special_q.rs` is an optimization layer over the line
sieve. For each **special prime** q in a chosen range [q_min, q_max], the sieve is restricted to
pairs (a, b) with q | N_alg(a, b). The restriction is enforced by the sieve condition:

```text
q | N_alg(a, b)  iff  a ≡ r_q·b (mod q)
```

for some root r_q of f mod q. Every candidate in the restricted set already has q as a known
algebraic factor. The remaining cofactor N_alg(a, b) / q is smaller and therefore more likely
to be smooth over the algebraic factor base.

The `SpecialQResult` type records the per-q output:

```rust
pub struct SpecialQResult {
    /// The special prime q.
    pub q: u64,
    /// The root r_q ∈ [0, q) of f mod q used for this sieve run.
    pub r_q: u64,
    /// Relations collected in this sieve run.
    pub relations: Vec<Relation>,
    /// The sieve area covered: ⌈(2A+1)/q⌉ × B pairs.
    pub restricted_area: u64,
}
```

Every relation in `relations` carries q in its algebraic exponent vector — this is the structural
invariant of the special-q strategy, checked by the KAT.

### The restricted sieve loop

The implementation runs the same log-sieve as the line sieve, but only visits a values satisfying
a ≡ r_q·b (mod q) for each b. This is implemented by stepping a in increments of q:

```rust
let rb_q = (r_q as u128 * b as u128 % q as u128) as u64;
let start_a_q = first_a_in_range(rb_q, q, -a_bound);
let mut a_q = start_a_q;
while a_q <= a_bound {
    // apply threshold filter and trial-divide
    a_q += q as i64;
}
```

The log-sieve array is still computed over the full range (for all primes in the factor bases),
but only the q-restricted positions are trial-divided. This is the correct approach: the sieve
contributions from other primes are still needed to filter candidates.

### Principle-4 annotation: yield advantage under-exposed at toy scale

**At toy scale, the yield advantage of the special-q strategy over the plain line sieve is
under-exposed.** The yield multiplier is a scale phenomenon:

- At cryptographic scale (B_alg ≈ 10⁶, A ≈ 10⁷), the algebraic norm N_alg(a, b) is large
  (hundreds of bits) and the probability of smoothness is low. The pre-guaranteed factor q
  significantly reduces the cofactor, making smoothness much more likely. The special-q strategy
  yields 5–10× more relations per sieve area than the plain line sieve.

- At toy scale (B_alg = 30, A = 10, B = 3), the norms are already small (tens of bits) and
  smooth with high probability. The pre-guaranteed factor q does not significantly improve the
  smoothness probability. The yield advantage is not observable.

The KAT (`kat_b_yield_comparison_principle4_annotated` in `special_q_kat.rs`) checks the
structural property (q in the algebraic exponent vector) and annotates the yield comparison as
under-exposed at toy scale, rather than asserting a yield improvement that would fail.

---

## 17. Lattice Sieving (Demonstration Fidelity)

### The lattice L_q and its natural basis

For a special prime q with root r_q, the set of (a, b) pairs satisfying the sieve restriction is
the **lattice**:

```text
L_q = { (a, b) ∈ ℤ² : a ≡ r_q·b (mod q) }
```

A natural basis for L_q is:

```text
v1 = (q, 0)    (check: q ≡ r_q·0 (mod q) ✓)
v2 = (r_q, 1)  (check: r_q ≡ r_q·1 (mod q) ✓)
```

Every lattice point s·v1 + t·v2 = (s·q + t·r_q, t) satisfies a ≡ r_q·b (mod q) since
s·q + t·r_q ≡ t·r_q ≡ r_q·t (mod q).

### Gauss 2D lattice reduction

The natural basis (v1, v2) is not efficient for enumeration: v1 = (q, 0) has length q, which
is large for large special primes. The **Gauss lattice reduction** algorithm finds a shorter,
more orthogonal basis (V1, V2):

```text
while |V1| > |V2|: swap V1, V2
V1 = V1 - round(dot(V1, V2) / dot(V2, V2)) * V2
```

Repeat until convergence. The reduced basis satisfies |V1| ≤ |V2| and
|dot(V1, V2)| ≤ |V2|² / 2 (the Gauss-reduced condition). The algorithm preserves the lattice
because each step is a unimodular transformation (determinant ±1).

The `LatticeBasis` type in `gnfs/src/sieve/lattice.rs` implements this:

```rust
pub struct LatticeBasis {
    pub v1: (i64, i64),
    pub v2: (i64, i64),
    pub q: u64,
    pub r_q: u64,
}

impl LatticeBasis {
    pub fn initial(q: u64, r_q: u64) -> Self { /* v1=(q,0), v2=(r_q,1) */ }
    pub fn gauss_reduce(&self) -> Self { /* Gauss reduction loop */ }
    pub fn in_lattice(&self, a: i64, b: i64) -> bool { /* a ≡ r_q·b (mod q) */ }
}
```

### Enumeration via the reduced basis

The lattice sieve enumerates (a, b) = s·V1 + t·V2 for integer (s, t) in a bounded region,
rather than stepping through a in increments of q for each b. The reduced basis vectors are
shorter than the original basis, so the enumeration visits fewer lattice points outside the
sieve region |a| ≤ A, 1 ≤ b ≤ B.

For each enumerated (a, b), the sieve applies the same log-threshold filter and trial-division
as the line sieve. The `LatticeSieveResult` records the per-(q, r_q) output:

```rust
pub struct LatticeSieveResult {
    pub q: u64,
    pub r_q: u64,
    pub basis: LatticeBasis,  // exposed for KAT inspection
    pub relations: Vec<Relation>,
    pub enumerated_points: u64,
}
```

### Principle-4 annotation: yield-per-area improvement is a scale phenomenon

**At toy scale, the lattice sieve produces the same (a, b) pairs as the special-q line sieve
for the same (q, r_q).** The two algorithms are mathematically equivalent: both enumerate L_q
in the region |a| ≤ A, 1 ≤ b ≤ B. The efficiency difference — the reduced basis visits fewer
wasted candidates outside the sieve region — is a constant factor that is swamped by the
overhead of the reduction and enumeration at small q.

At cryptographic scale, the reduced basis has vectors of length ≈ √q, so the enumeration covers
≈ A·B / q lattice points with minimal waste. For q ≈ 10⁶ and A, B ≈ 10⁷, this is a significant
efficiency gain over stepping by q.

The KAT (`kat_b_lattice_sieve_subset_of_special_q_sieve` in `lattice_kat.rs`) verifies that the
lattice sieve produces a subset of the special-q sieve relations for the same (q, r_q), and
annotates the yield-per-area comparison as under-exposed at toy scale.

---

## 18. The Exponent Vector and the Relation Matrix

### The `ExponentVector` type

The `ExponentVector` type in `gnfs/src/sieve/mod.rs` stores the sparse exponent representation
for one side of a relation:

```rust
pub struct ExponentVector {
    /// Sparse: (factor-base index, exponent) pairs, sorted by index.
    /// Exponents are always > 0 (zeros are not stored).
    pub entries: Vec<(usize, u32)>,
}
```

The exponent type is `u32`, not `u8` or `bool`. This accommodates:

- **NFS-factoring (G.E):** exponents are small (typically 1–3), reduced mod 2 for the GF(2)
  nullspace computation.
- **NFS-DL (D.A):** exponents are reduced mod ℓ where ℓ is the target group order; ℓ can be
  large, but exponents before reduction are still small integers.

The `u32` type is the smallest that accommodates both without overflow risk. This is the
over-specification for D.A: storing integer exponents rather than pre-reduced GF(2) parities
means D.A can read the exponents directly without resharding the relation format.

### The GF(2) row methods

The `Relation` type provides two methods for G.E:

```rust
/// GF(2) row for the rational side: sign column prepended, then exponent parities.
pub fn rational_row_gf2(&self, fb: &FactorBase) -> Vec<bool> {
    let mut row = vec![false; 1 + fb.rational_size()];
    row[0] = self.rational_sign;  // sign/−1 column
    for (idx, exp) in self.rational_exponents.iter() {
        row[1 + idx] = (exp % 2) == 1;
    }
    row
}

/// GF(2) row for the algebraic side: exponent parities, obstruction columns appended as zeros.
pub fn algebraic_row_gf2(&self, fb: &FactorBase) -> Vec<bool> {
    let mut row = vec![false; fb.algebraic_size() + fb.obstruction_count];
    for (idx, exp) in self.algebraic_exponents.iter() {
        if idx < fb.algebraic_size() { row[idx] = (exp % 2) == 1; }
    }
    // Obstruction columns remain false; G.E fills them in.
    row
}
```

The sign column is prepended to the rational row: bit 0 is 1 iff `rational_sign` is true. The
obstruction columns (quadratic characters) are appended to the algebraic row as zeros; G.E fills
them in when it constructs the full matrix.

### Why integer exponents, not GF(2) parities

The design decision to store full integer exponents (not pre-reduced GF(2) parities) is the
load-bearing cross-track call in the C-Relation contract. Three consumers need different views:

- **G.D (filtering):** reads the exponent vectors to detect duplicate or linearly dependent
  relations. Integer exponents are needed to identify exact duplicates.
- **G.E (linear algebra):** reduces exponents mod 2 via `rational_row_gf2` / `algebraic_row_gf2`
  to build the GF(2) matrix. The `u32` exponents are reduced on demand.
- **D.A (NFS-DL):** reads the integer exponents directly for GF(ℓ) linear algebra, where ℓ is
  the target group order. Pre-reducing to GF(2) would destroy the information D.A needs.

Re-narrowing C-Relation after G.E or D.A consumes it would be a destructive reshard. The
over-specification (integer exponents) is the correct design for a cross-track contract.

---

## 19. Downstream Reuse: G.D, G.E, and D.A

### G.D (filtering)

The relation corpus from the sieve feeds G.D (filtering), which removes:

- **Duplicate relations:** the same (a, b) pair found by multiple sieve runs (e.g., the line
  sieve and the special-q sieve may both find the same relation).
- **Linearly dependent relations:** relations whose exponent vectors are linearly dependent over
  GF(2), which would not contribute to the null space.

G.D reads the `ExponentVector` entries and the `FactorBase` column indices. The `AlgebraicPrime::index`
field and the `rational_index` / `algebraic_index` lookup methods provide the column mapping.

### G.E (linear algebra)

G.E (linear algebra over GF(2)) takes the filtered relation corpus and finds a subset S such that
the sum of all exponent vectors (mod 2) is the zero vector. This is the null space computation:

1. Build the relation matrix: one row per relation, one column per factor-base prime (plus sign
   and quadratic-character columns). Each row is `rational_row_gf2` concatenated with
   `algebraic_row_gf2`.
2. Find the null space of this matrix over GF(2) (e.g., via structured Gaussian elimination).
3. Each null-space vector identifies a subset S of relations whose exponent product is a perfect
   square on both sides.

The `matrix_width()` method on `FactorBase` gives the total column count:
rational_size + algebraic_size + obstruction_count.

### D.A (NFS-DL)

The NFS discrete logarithm algorithm (D.A) uses the same relation structure as NFS factoring,
but interprets exponents mod ℓ (the target group order) instead of mod 2. The integer exponents
stored in `ExponentVector` support this interpretation directly.

D.A may add **Schirokauer map** columns to the algebraic exponent vector. The sparse
`(index, exponent)` representation accommodates additional entries without resharding the
`Relation` type. This is the cross-track over-specification: the `Relation` type is designed so
D.A can read it without modification.

---

## 20. KAT Summary (G.C)

The following table lists the key known-answer tests across the sieve module, with the
mathematical fact each one verifies. All tests are in `gnfs/tests/`.

| Test | File | Mathematical fact verified |
|------|------|---------------------------|
| `kat1_algebraic_factor_base_matches_brute_force` | `factor_base_kat.rs` | Algebraic FB = brute-force root enumeration mod each p ≤ B_alg |
| `kat1b_index_lookups_consistent` | `factor_base_kat.rs` | `rational_index` / `algebraic_index` round-trip correctly |
| `kat2_norm_reconstruction_and_relation_verify` | `factor_base_kat.rs` | `Relation::verify()` holds for hand-constructed smooth relation |
| `kat2b_verify_fails_on_perturbed_exponent` | `factor_base_kat.rs` | `verify()` fails when an exponent is perturbed |
| `kat2c_relation_new_rejects_non_coprime` | `factor_base_kat.rs` | `Relation::new` returns None for non-coprime (a, b) |
| `kat2d_relation_new_rejects_partial_smoothness` | `factor_base_kat.rs` | `Relation::new` returns None when cofactor > 1 |
| `kat3_norm_bridge_range` | `factor_base_kat.rs` | Toy-scale norm fits Uint<4>; 2^257 overflows with error |
| `kat3b_norm_sign` | `factor_base_kat.rs` | `norm_sign` returns true iff norm < 0 |
| `kat3c_norm_bridge_round_trip_with_trial_smooth` | `factor_base_kat.rs` | `norm_to_uint` + `trial_smooth` round-trip for known smooth norm |
| `kat_a_sieve_produces_relations_all_verify` | `line_sieve_kat.rs` | Sieve finds ≥ 5 relations; all satisfy `verify()` |
| `kat_a2_all_relations_are_coprime_and_smooth` | `line_sieve_kat.rs` | All returned relations are coprime and fully smooth |
| `kat_b_deterministic_relation_count` | `line_sieve_kat.rs` | Same parameters → same relations (determinism) |
| `kat_b2_lower_threshold_finds_at_least_as_many_relations` | `line_sieve_kat.rs` | Lower threshold ≥ higher threshold relation count |
| `kat_b3_exact_count_is_stable` | `line_sieve_kat.rs` | Exact count is stable across runs (pins the implementation) |
| `kat_a_relations_verify_and_carry_q` | `special_q_kat.rs` | All special-q relations verify and carry q in algebraic EV |
| `kat_a2_sieve_restriction_enforced` | `special_q_kat.rs` | All (a, b) satisfy a ≡ r_q·b (mod q) |
| `kat_a3_spot_check_known_relation_q7` | `special_q_kat.rs` | (a=5, b=1) found for q=7, r_q=5; carries q=7 in algebraic EV |
| `kat_b_yield_comparison_principle4_annotated` | `special_q_kat.rs` | Yield advantage annotated as under-exposed at toy scale |
| `kat_b2_special_q_relations_subset_of_line_sieve` | `special_q_kat.rs` | Special-q relations ⊆ line sieve relations for same q |
| `kat_a_all_relations_lie_in_lattice` | `lattice_kat.rs` | All lattice-sieve (a, b) satisfy a ≡ r_q·b (mod q) |
| `kat_a2_all_relations_verify` | `lattice_kat.rs` | All lattice-sieve relations satisfy `verify()` |
| `kat_a3_relations_carry_q_in_algebraic_exponents` | `lattice_kat.rs` | All lattice-sieve relations carry q in algebraic EV |
| `kat_a4_reduced_basis_q7_r5_is_correct` | `lattice_kat.rs` | Gauss-reduced basis for q=7, r_q=5 is correct |
| `kat_b_lattice_sieve_subset_of_special_q_sieve` | `lattice_kat.rs` | Lattice-sieve relations ⊆ special-q sieve relations |
| `kat_c_yield_comparison_principle4_annotated` | `lattice_kat.rs` | Yield-per-area improvement annotated as under-exposed |

---

## 21. What Comes Next

Sieving is the third stage of the GNFS pipeline. The output — a corpus of `Relation` objects,
each with both exponent vectors and the sign bit — is the input to the filtering and linear
algebra stages.

**G.D (filtering)** takes the raw relation corpus and removes duplicates and linearly dependent
relations. The `ExponentVector` entries and the `FactorBase` column indices are the interface
between G.C and G.D.

**G.E (linear algebra over GF(2))** takes the filtered corpus and finds a subset S of relations
whose exponent product is a perfect square on both sides. The `rational_row_gf2` and
`algebraic_row_gf2` methods on `Relation` produce the GF(2) matrix rows. The `matrix_width()`
method on `FactorBase` gives the column count.

**G.F (square root)** takes the subset S from G.E and computes the congruence of squares
mod N that yields the non-trivial factor. This stage uses the number-field arithmetic from
`shared-numfield` and the algebraic structure of K = ℚ(α).

**D.A (NFS-DL)** adapts the same relation structure for discrete logarithm computation. The
integer exponents in `ExponentVector` are read directly (not reduced mod 2) and used in GF(ℓ)
linear algebra. The `Relation` type is designed to support this adaptation without resharding.

---

## Further Reading (G.C)

1. **Pomerance, C. (1996).** *A tale of two sieves.* Notices of the AMS, 43(12), 1473–1485.
   An accessible introduction to the quadratic sieve and the number field sieve, with the
   sieving step explained at a high level.

2. **Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) (1993).** *The Development of the Number
   Field Sieve.* Springer LNM 1554. The original papers on GNFS, including the sieving step
   and the algebraic side.

3. **Briggs, M. E. (1998).** *An introduction to the general number field sieve.* Master's
   thesis, Virginia Tech. A self-contained exposition of the GNFS pipeline including sieving,
   with worked examples at toy scale.

4. **Kleinjung, T., Aoki, K., Franke, J., et al. (2010).** *Factorization of a 768-bit RSA
   modulus.* CRYPTO 2010. The RSA-768 factorisation, which used the special-q lattice sieve
   as the primary relation-collection strategy.

5. **Schirokauer, O. (1993).** *Discrete logarithms and local units.* Philosophical Transactions
   of the Royal Society A, 345(1676), 409–423. The Schirokauer maps used in D.A (NFS-DL) to
   handle the units in the number field — the D.A extension of the C-Relation contract.

---

# Filtering: A Code-Tour Chapter

This chapter explains the `gnfs/src/filter` module — the filtering stage of the General Number
Field Sieve. It is organised by the implementation, not by mathematical abstraction: each section
introduces the mathematics, then shows how it is realised in code. A reader who has read the
sieving chapter (§12–§21 above) can read this one without consulting the source.

The chapter covers the transition from a relation corpus to a sparse GF(2) matrix, the graph view
of that matrix, singleton removal, excess and clique pruning, column merging, and the provenance
map that threads through to G.F.

---

## 22. The Relation Corpus as a Sparse Matrix

### From `Vec<Relation>` to a GF(2) matrix

The sieve (G.C) produces a `Vec<Relation>` — a corpus of smooth pairs (a, b). Each relation
carries two exponent vectors: one over the rational factor base, one over the algebraic factor
base. The filtering stage (G.D) converts this corpus into a sparse matrix over GF(2).

The conversion is straightforward:

- **One row per relation.** Relation i contributes row i.
- **One column per factor-base element.** Each rational prime p and each algebraic ideal (p, r)
  contributes one column.
- **The entry is the GF(2) parity of the exponent.** If prime p divides the rational norm to an
  odd exponent, the entry is 1; if even (or absent), it is 0.

The column layout is fixed by `FactorBase::matrix_width()`:

$$
\text{columns} = \underbrace{\text{rational\_size}}_{\text{rational primes}} + \underbrace{\text{algebraic\_size}}_{\text{algebraic ideals}} + \underbrace{\text{obstruction\_count}}_{\text{sign + quad. chars.}}
$$

Concretely, for the toy setup in `merge_kat.rs` (f(x) = x³ − x − 1, B_rat = B_alg = 13):

- Rational primes ≤ 13: {2, 3, 5, 7, 11, 13} → columns 0–5 (`rational_size = 6`).
- Algebraic ideals: (5, 2), (7, 5), (11, 6) → columns 6–8 (`algebraic_size = 3`).
- Obstruction column: sign bit → column 9 (`obstruction_count = 1`).
- Total: `matrix_width = 10`.

The `build_matrix` function in `gnfs/src/filter/mod.rs` performs this construction:

```rust
pub fn build_matrix(relations: &[Relation], fb: &FactorBase) -> SparseMatrix {
    // ...
    for (i, relation) in relations.iter().enumerate() {
        let rat_row = relation.rational_row_gf2(fb);  // [sign_bit, rat_col_0, ...]
        let alg_row = relation.algebraic_row_gf2(fb); // [alg_col_0, ..., zeros...]

        let mut cols: Vec<usize> = Vec::new();
        // Rational columns: local index 1+k → global column k.
        for k in 0..rat_size {
            if rat_row[1 + k] { cols.push(k); }
        }
        // Algebraic columns: local index k → global column rat_size + k.
        for k in 0..alg_size {
            if alg_row[k] { cols.push(rat_size + k); }
        }
        // Sign bit → global obstruction column.
        if rat_row[0] { cols.push(obstruction_col_start); }

        rows.push(MatrixRow { cols, provenance: vec![i] });
    }
    // ...
}
```

The `rational_row_gf2` and `algebraic_row_gf2` methods on `Relation` (defined in G.C) produce
the GF(2) parities. The sign bit is placed at the first obstruction column, not at column 0 of
the rational block — this is the re-mapping that `build_matrix` performs.

### Why GF(2): the parity condition

The goal of the linear algebra step (G.E) is to find a subset S of relations such that the
product of all rational norms is a perfect square in ℤ and the product of all algebraic norms
is a perfect square in the number field K = ℚ(α). A product of integers is a perfect square
iff every prime appears to an even total exponent. This is a **parity condition**: we need the
sum of all exponent vectors (mod 2) to be the zero vector. Working over GF(2) encodes exactly
this condition.

The sign column encodes the sign parity: the product of all selected rational norms must be
positive (an even number of negative norms). The quadratic-character columns (filled by G.E)
encode additional algebraic-side parity conditions needed to guarantee that the algebraic square
root is well-defined.

---

## 23. The Graph View

### Relations as edges, primes as nodes

The sparse GF(2) matrix has a natural graph interpretation that illuminates the filtering
algorithms:

- **Nodes** are the non-obstruction columns — the rational primes and algebraic ideals.
- **Edges** are the rows — the relations. A relation "connects" the primes it contains: if
  relation R has a 1 in columns c₁ and c₂, then R is an edge between nodes c₁ and c₂.

A **GF(2) dependency** is a subset of relations whose column-sum (mod 2) is the zero vector.
In graph terms, this is a **cycle**: a set of edges such that every node is incident to an even
number of edges in the set. The null space of the matrix over GF(2) is exactly the set of all
such cycles.

### Why singletons cannot appear in a dependency

A **singleton column** is a column of Hamming weight 1 — a prime or ideal that appears in
exactly one surviving relation. In graph terms, it is a node of degree 1 (a leaf).

A leaf cannot be part of any cycle. If a node has degree 1, the unique edge incident to it must
appear in any cycle that includes that node. But then the node is incident to exactly one edge
in the cycle — an odd count — which violates the cycle condition. Therefore, no cycle can
include a leaf node, and the relation (edge) containing the singleton column can never be part
of a GF(2) dependency.

This is the mathematical justification for singleton removal: the row containing a singleton
column contributes nothing to the null space and can be safely discarded.

---

## 24. Singleton Removal

### The fixpoint algorithm

Singleton removal is implemented in `gnfs/src/filter/singleton.rs` as a fixpoint iteration:

```rust
pub fn remove_singletons(mut matrix: SparseMatrix) -> SparseMatrix {
    loop {
        let mut rows_to_remove: Vec<usize> = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for col in 0..matrix.obstruction_col_start {
            if matrix.col_weights[col] == 1 {
                // Find the unique row containing this column.
                for (row_idx, row) in matrix.rows.iter().enumerate() {
                    if row.cols.binary_search(&col).is_ok() {
                        if seen.insert(row_idx) {
                            rows_to_remove.push(row_idx);
                        }
                        break;
                    }
                }
            }
        }

        if rows_to_remove.is_empty() { break; }  // fixpoint reached

        // Remove in descending order to preserve index validity.
        rows_to_remove.sort_unstable();
        rows_to_remove.dedup();
        for &row_idx in rows_to_remove.iter().rev() {
            matrix.remove_row(row_idx);
        }
    }
    matrix
}
```

The loop collects all rows that contain a weight-1 non-obstruction column, removes them, and
repeats. It terminates when no weight-1 non-obstruction column remains.

### Why a single pass is insufficient: cascading singletons

Removing a row decrements the weight of every column it contains. A column that had weight 2
becomes weight 1 after its partner row is removed — a new singleton. A single pass over the
column list would miss these newly created singletons. The fixpoint loop handles cascading
singletons correctly: after each batch of removals, the column weights are updated and the
scan restarts.

The KAT in `merge_kat.rs` (KAT d) demonstrates cascading removal: relation R4 contains the
only occurrence of ideal (7, 5) (a singleton), so R4 is removed. After R4's removal, ideal
(5, 2) drops to weight 1 (only R3 remains), so R3 is removed. After R3's removal, prime 7
drops to weight 2 (R0 and R1 remain) — the cascade stops.

### Obstruction columns are exempt

The obstruction columns (sign bit and quadratic-character columns) are structural: they encode
constraints that G.E needs, not factor-base primes. A weight-1 sign column does not mean that
the relation containing it is useless — it means that exactly one relation has a negative
rational norm. Treating it as a singleton would incorrectly discard that relation. The loop
scans only `col in 0..matrix.obstruction_col_start`.

### Provenance is unchanged by singleton removal

Singleton removal drops rows, never merges them. Each surviving row's provenance remains its
original singleton set `[i]` — the index of the original relation. This is the design decision
documented in the `remove_singletons` docstring:

> Provenance is preserved unchanged: singleton removal drops rows, never merges them, so each
> surviving row's provenance is its original singleton set.

The provenance map is the thread that G.F uses to recover original (a, b) pairs from a null-space
vector. Singleton removal does not disturb it.

---

## 25. Excess and the Clique/Pruning Step

### Defining excess

The **excess** of the matrix is:

$$
\text{excess} = \text{rows} - (\text{columns} - \text{obstruction\_count})
$$

implemented as:

```rust
pub fn excess(&self) -> isize {
    self.rows.len() as isize - (self.num_cols as isize - self.obstruction_count as isize)
}
```

The denominator `columns − obstruction_count` is the number of non-obstruction columns — the
"active" dimension of the matrix. For the null space to be non-trivial (i.e., for there to
exist at least one non-zero null-space vector), the matrix must have more rows than active
columns: excess must be positive.

The constant `EXCESS_FLOOR = 20` is the minimum excess that the pruning step must preserve:

```rust
/// Minimum excess G.D.2 pruning must preserve.
///
/// excess = rows − (columns − obstruction_count). At toy scale any positive excess
/// suffices for a non-trivial nullspace; 20 is a conservative floor that keeps the
/// matrix well-overdetermined.
pub const EXCESS_FLOOR: usize = 20;
```

**Principle-4 annotation (scale-dependent floor).** At toy scale, any positive excess suffices
for a non-trivial null space — the floor of 20 is conservative and is never approached in
practice because the toy corpus is small. At cryptographic scale (RSA-768+), the excess floor
is typically set to ~200 or a fraction of the column count, balancing the need for many
independent null-space vectors (each gives one congruence-of-squares candidate) against the
cost of carrying extra rows through the linear algebra step. The floor is the knob that trades
matrix density against null-space dimension; at toy scale this trade-off is invisible.

### The clique view and greedy pruning

In the graph view, a **clique** is a set of relations (edges) that all share a common prime
(node). A clique is dense: every pair of relations in the clique shares at least one column.
Dense cliques increase the total matrix weight (total number of set entries) without
proportionally increasing the null-space dimension.

The `prune_cliques` function in `gnfs/src/filter/merge.rs` approximates clique pruning by a
greedy heuristic: repeatedly remove the heaviest row (the relation with the most set columns)
as long as the excess exceeds `EXCESS_FLOOR`:

```rust
pub fn prune_cliques(mut matrix: SparseMatrix) -> SparseMatrix {
    loop {
        if matrix.excess() <= EXCESS_FLOOR as isize { break; }
        if matrix.rows.is_empty() { break; }

        // Find the row with the maximum weight (number of set columns).
        let max_weight = matrix.rows.iter().map(|r| r.cols.len()).max().unwrap_or(0);
        if max_weight == 0 { break; }

        let row_idx = matrix.rows.iter().enumerate()
            .find(|(_, r)| r.cols.len() == max_weight)
            .map(|(i, _)| i).unwrap();

        matrix.remove_row(row_idx);
    }
    matrix
}
```

The heaviest rows are the most "connected" relations — those that share columns with many other
relations. Removing them reduces the total matrix weight while preserving the excess floor.
Ties are broken by row index (lowest index first) for determinism.

The KAT `kat_prune_cliques_respects_excess_floor` in `merge_kat.rs` verifies that after
pruning, `excess() == EXCESS_FLOOR` exactly (when the initial excess is above the floor).

---

## 26. Merging

### Eliminating a column by XOR

A column of weight 2 appears in exactly two rows. XOR-merging those two rows eliminates the
column: the shared column cancels in the symmetric difference (GF(2) addition), and the merged
row contains only the columns that appeared in exactly one of the two source rows.

The `MatrixRow::xor_merge` method computes the symmetric difference of column sets and the
union of provenance sets in a single O(n + m) merge pass:

```rust
pub fn xor_merge(&self, other: &MatrixRow) -> MatrixRow {
    // Symmetric difference of two sorted vecs — O(n + m).
    let mut cols = Vec::new();
    let mut i = 0; let mut j = 0;
    while i < self.cols.len() && j < other.cols.len() {
        match self.cols[i].cmp(&other.cols[j]) {
            Ordering::Less    => { cols.push(self.cols[i]); i += 1; }
            Ordering::Greater => { cols.push(other.cols[j]); j += 1; }
            Ordering::Equal   => { /* XOR = 0, skip both */ i += 1; j += 1; }
        }
    }
    cols.extend_from_slice(&self.cols[i..]);
    cols.extend_from_slice(&other.cols[j..]);
    // ... union of provenance sets (same pattern) ...
    MatrixRow { cols, provenance }
}
```

### 2-way merge: two rows sharing a column

For a weight-2 column c with rows R₀ and R₁:

1. Compute `merged = R₀.xor_merge(R₁)`.
2. Remove R₀ and R₁ from the matrix.
3. Append `merged`.
4. Update `col_weights`: decrement for all columns in R₀ and R₁, increment for all columns
   in `merged`. The shared column c is decremented twice and not incremented → weight drops
   to 0 (eliminated).

The net effect on row count: two rows become one (net −1). The net effect on column count:
column c is eliminated (weight 0). The excess formula changes: both the numerator and the
effective denominator decrease, so the excess may increase or decrease depending on the
specific column weights.

### k-way merge for k = 3

For a weight-3 column c with rows R₀, R₁, R₂, the merge is:

```text
tmp = R₀ ⊕ R₁
merged = tmp ⊕ R₂
```

implemented in `merge_three_rows`. The three source rows are removed and `merged` is appended.
The shared column c is decremented three times and not incremented → weight drops to 0.

The `merge_columns` function performs two passes:

```rust
pub fn merge_columns(mut matrix: SparseMatrix) -> SparseMatrix {
    matrix = merge_pass(&mut matrix, 2);  // weight-2 pass
    matrix = merge_pass(&mut matrix, 3);  // weight-3 pass (after re-scan)
    matrix
}
```

The weight-2 pass is performed first; after it completes, some weight-3 columns may have
become weight-2 or weight-1 (because earlier merges in the weight-2 pass changed their
weights). The weight-3 pass re-scans `col_weights` and processes the remaining weight-3
columns.

### The Cavallar weight-cost heuristic

In production NFS, merges are ordered by a **weight-cost heuristic** (Cavallar 1996): prefer
merges that minimise the increase in total matrix weight. The total weight is the sum of all
`col_weights` — equivalently, the total number of set entries in the matrix. A merge of a
weight-k column eliminates that column (saving k entries) but the merged row may be heavier
than the sum of the source rows minus the shared column (adding entries). The Cavallar
heuristic picks the merge with the smallest net weight increase.

This implementation uses a simplified Cavallar ordering: process columns in order of increasing
weight (weight-2 first, then weight-3), breaking ties by column index. This is correct in
principle — lower-weight columns are cheaper to merge — but the full heuristic would also
consider the weight of the merged row.

**Principle-4 annotation (Cavallar weight-cost ordering).** At toy scale, the matrix is small
enough that any merge order gives a tractable matrix for G.E. The weight-saving the Cavallar
heuristic buys is under-exposed: the difference between the simplified ordering and the optimal
ordering is a few percent of total weight, which is invisible when the matrix has tens of rows.
At cryptographic scale (RSA-768+), the matrix has millions of rows and hundreds of thousands
of columns; the Cavallar heuristic can reduce the total weight by 30–50% compared to a naive
ordering, which directly reduces the cost of the structured Gaussian elimination in G.E. The
module docstring in `merge.rs` annotates this disconnect explicitly.

---

## 27. The Provenance Map as the Thread to G.F

### What provenance tracks

Every `MatrixRow` carries a `provenance: Vec<usize>` — a sorted, deduplicated list of original
relation indices (indices into the `Vec<Relation>` passed to `build_matrix`):

```rust
pub struct MatrixRow {
    /// Sorted column indices where this row has a 1 in GF(2).
    pub cols: Vec<usize>,
    /// Sorted, deduplicated original relation indices this row derives from.
    pub provenance: Vec<usize>,
}
```

For a freshly built row, `provenance = [i]` (the original relation index). Under merging,
`provenance` is combined by sorted union: if R₀ has provenance `[0, 2]` and R₁ has provenance
`[1, 3]`, then `R₀.xor_merge(R₁)` has provenance `[0, 1, 2, 3]`.

Singleton removal drops rows without merging, so provenance is unchanged. Clique pruning also
drops rows without merging. Only `merge_columns` grows provenance sets.

### The G.F expansion

G.E (linear algebra) finds a null-space vector — a subset of rows in the final matrix whose
GF(2) column-sum is the zero vector. G.F (square root) needs to expand this subset back to
the original (a, b) pairs to compute the congruence of squares.

The expansion is:

1. G.E selects a subset of rows in the final matrix (the null-space vector).
2. For each selected row, collect its `provenance` set.
3. Take the union of all collected provenance sets — this gives the set of original relation
   indices.
4. Look up the original `Relation` objects by index to recover the (a, b) pairs.

The KAT `kat_d_end_to_end_provenance` in `merge_kat.rs` verifies this property: for each row
in the final matrix, the XOR of the GF(2) column sets of all original relations indexed by its
provenance equals that row's column set.

### Why provenance stores original indices, not pre-reduced row sums

The provenance map stores original relation indices, not the GF(2) column set of the merged
row. This is the over-specification documented in the C-Matrix contract:

> Over-specify the provenance map now: store original relation indices, not pre-reduced row
> sums, so G.F can recover the actual (a, b) pairs.

The reason is that G.F needs the actual (a, b) pairs to compute the algebraic square root —
the product of (a − b·α) in the number field K = ℚ(α). The GF(2) column set of a merged row
is the XOR of the original exponent parities; it does not contain the original (a, b) values.
Storing the GF(2) column set instead of the original indices would make G.F impossible without
re-deriving the (a, b) pairs from the exponent vectors, which would require keeping the full
`Vec<Relation>` in scope anyway. The provenance map is the clean interface: G.F receives the
final matrix and the original `Vec<Relation>`, and uses the provenance map to bridge them.

---

## 28. What G.E Inherits

### The `SparseMatrix` contract (C-Matrix)

The output of the full filtering pipeline — `build_matrix` → `remove_singletons` →
`prune_cliques` → `merge_columns` — is a `SparseMatrix` that G.E consumes directly:

```rust
pub struct SparseMatrix {
    /// The rows of the matrix, each with its column set and provenance.
    pub rows: Vec<MatrixRow>,
    /// Total number of columns = FactorBase::matrix_width().
    pub num_cols: usize,
    /// Index of the first obstruction column = rational_size + algebraic_size.
    pub obstruction_col_start: usize,
    /// Number of obstruction columns = FactorBase::obstruction_count.
    pub obstruction_count: usize,
    /// col_weights[c] = number of rows with a 1 in column c.
    pub col_weights: Vec<u32>,
}
```

The key properties G.E inherits:

- **Dimensions.** `rows.len()` rows, `num_cols` columns. The excess `rows.len() − (num_cols −
  obstruction_count)` is ≥ `EXCESS_FLOOR` (guaranteed by `prune_cliques`), ensuring a
  non-trivial null space.

- **Obstruction columns.** The sign column (at `obstruction_col_start`) is populated by
  `build_matrix` from `relation.rational_sign`. The quadratic-character columns (at
  `obstruction_col_start + 1` through `num_cols − 1`) are carried as zeros by G.D — G.E fills
  them in before running the null-space computation.

- **Row-major store + `col_weights` side table.** The row-major representation (each row is a
  sorted `Vec<usize>` of set column indices) supports the row-XOR operations that structured
  Gaussian elimination performs. The `col_weights` side table supports column-weight queries
  without scanning all rows — G.E uses it to choose pivot columns.

- **Provenance map.** Each row carries its provenance set. G.E does not modify provenance; it
  passes the final matrix (with provenance intact) to G.F.

### The pipeline in full

The complete G.D pipeline, as exercised by the KATs in `merge_kat.rs`:

```rust
// 1. Build the initial matrix from the relation corpus.
let matrix = build_matrix(&relations, &fb);

// 2. Remove singleton columns to fixpoint.
let matrix = remove_singletons(matrix);

// 3. Prune heavy rows while excess > EXCESS_FLOOR.
let matrix = prune_cliques(matrix);

// 4. Eliminate weight-2 and weight-3 columns by XOR-merging.
let matrix = merge_columns(matrix);

// matrix is now ready for G.E (linear algebra over GF(2)).
```

Each step is a pure function (consumes and returns `SparseMatrix`), making the pipeline
composable and testable in isolation.

---

## 29. KAT Summary (G.D)

The following table lists the key known-answer tests across the filter module, with the
mathematical fact each one verifies. All tests are in `gnfs/tests/merge_kat.rs`.

| Test | Mathematical fact verified |
|------|---------------------------|
| `kat_a_two_way_merge_correctness` | Weight-2 column eliminated; merged row has correct cols (symmetric difference) and provenance (union) |
| `kat_b_determinism` | Full pipeline (build → singletons → prune → merge) is deterministic for a fixed corpus |
| `kat_c_cado_nfs_oracle` | (Ignored) CADO-NFS oracle for filtered matrix dimensions — gated when CADO absent |
| `kat_d_end_to_end_provenance` | For each final row, XOR of provenance relations' GF(2) column sets equals the row's cols |
| `kat_prune_cliques_respects_excess_floor` | After pruning, `excess() == EXCESS_FLOOR` exactly when initial excess > floor |

---

## 30. What Comes Next

Filtering is the fourth stage of the GNFS pipeline. The output — a `SparseMatrix` with
provenance map — is the input to the linear algebra and square-root stages.

**G.E (linear algebra over GF(2))** takes the filtered matrix and finds its null space. It
fills in the quadratic-character obstruction columns (currently zero in the G.D output), then
runs structured Gaussian elimination to find a basis for the null space. Each null-space vector
identifies a subset of rows whose GF(2) column-sum is zero — a candidate congruence of squares.

**G.F (square root)** takes a null-space vector from G.E and expands it through the provenance
map to recover the original (a, b) pairs. It then computes the product of all rational norms
(a perfect square in ℤ) and the product of all algebraic norms (a perfect square in K = ℚ(α)),
and extracts the square roots to form the congruence x² ≡ y² (mod N). A non-trivial GCD of
x − y and N yields a factor.

**D.A (NFS-DL)** adapts the same filtering infrastructure for discrete logarithm computation.
The `SparseMatrix` and provenance map are reused directly; the difference is in the linear
algebra step (GF(ℓ) instead of GF(2)) and the square-root step (Schirokauer maps instead of
the algebraic square root).

---

## Further Reading (G.D)

1. **Cavallar, S. (2002).** *Strategies in filtering in the number field sieve.* In: Bosma, W.
   (ed.) Algorithmic Number Theory (ANTS-IV), LNCS 1838. Springer. The original source for the
   weight-cost merge heuristic and the clique/excess balance.

2. **Buhler, J., Lenstra, H. W. Jr., and Pomerance, C. (1993).** *Factoring integers with the
   number field sieve.* In: Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) The Development of
   the Number Field Sieve, LNM 1554. Springer. The original description of the filtering step
   and the relation matrix.

3. **Franke, J., and Kleinjung, T. (2006).** *Continued fractions and lattice sieving.* CADO
   workshop on integer factorization. Describes the filtering pipeline used in CADO-NFS,
   including the singleton removal and merge strategies.

4. **Pomerance, C. (1996).** *A tale of two sieves.* Notices of the AMS, 43(12), 1473–1485.
   An accessible introduction to the relation matrix and the linear algebra step, with the
   filtering step explained at a high level.

5. **Wiedemann, D. (1986).** *Solving sparse linear equations over finite fields.* IEEE
   Transactions on Information Theory, 32(1), 54–62. The Wiedemann algorithm used in G.E for
   the null-space computation — the downstream consumer of the filtered `SparseMatrix`.

---

# Linear Algebra: A Code-Tour Chapter

This chapter explains the `gnfs/src/linalg` module — the linear algebra stage of the General
Number Field Sieve. It is organised by the implementation, not by mathematical abstraction: each
section introduces the mathematics, then shows how it is realised in code. A reader who has read
the filtering chapter (§21–§30 above) can read this one without consulting the source.

The chapter covers the filtered matrix as a GF(2) linear system, the nullspace as a congruence
of squares, the role of sign and quadratic-character columns, block Lanczos, block Wiedemann,
the kernel-vector → provenance → original-relations thread, and the C-LinAlg contract.

---

## 31. The Filtered Matrix as a GF(2) Linear System

### From filtering to linear algebra

The output of G.D filtering is a `SparseMatrix` — a sparse GF(2) matrix A with m rows and n
columns. Each row corresponds to a filtered relation (or a merged aggregate of relations), and
each column corresponds to a prime or obstruction entry in the factor base. The entry A[i][j]
is 1 if the j-th prime appears to an odd total exponent in the i-th relation's norm product,
and 0 otherwise.

The goal of G.E is to find vectors in the **left nullspace** of A: subsets of rows whose GF(2)
column-sum is the zero vector. Formally, we seek non-zero vectors v ∈ GF(2)^m such that:

$$
v^T A = 0 \pmod{2}
$$

Equivalently, we seek vectors in the nullspace of A^T (the right nullspace of the transpose).

### Why the left nullspace

A left nullspace vector v is a 0/1 indicator over rows. The rows it selects are a subset S of
filtered relations. The GF(2) column-sum being zero means: for every prime p in the factor
base, the total exponent of p across all relations in S is even. This is exactly the condition
for the product of the corresponding norms to be a perfect square.

### The `SparseMatrix` contract (C-Matrix)

G.E inherits the `SparseMatrix` from G.D unchanged:

```rust
pub struct SparseMatrix {
    /// The rows of the matrix, each with its column set and provenance.
    pub rows: Vec<MatrixRow>,
    /// Total number of columns = FactorBase::matrix_width().
    pub num_cols: usize,
    /// Index of the first obstruction column = rational_size + algebraic_size.
    pub obstruction_col_start: usize,
    /// Number of obstruction columns = FactorBase::obstruction_count.
    pub obstruction_count: usize,
    /// col_weights[c] = number of rows with a 1 in column c.
    pub col_weights: Vec<u32>,
}
```

The column layout is:

```text
[  rational primes  |  algebraic ideals  |  sign  |  QC columns  ]
 0 .. rational_size   rational_size ..      obs_start  obs_start+1 ..
                       obs_start            (col 0)    (cols 1..num_qc)
```

G.D populates the rational and algebraic columns and the sign column. G.E fills in the
quadratic-character (QC) columns before running the nullspace solver.

### Toy-scale matrix dimensions

The KAT matrices in `gnfs/tests/lanczos_kat.rs` and `gnfs/tests/wiedemann_kat.rs` use a
hand-built 6 × 4 matrix:

```text
    col: 0 1 2 3
row 0:   1 1 0 0   ← provenance [0]
row 1:   0 1 1 0   ← provenance [1]
row 2:   1 0 1 0   ← provenance [2]
row 3:   1 1 0 0   ← provenance [3]  (duplicate of row 0)
row 4:   0 0 0 1   ← provenance [4]
row 5:   0 0 0 1   ← provenance [5]  (duplicate of row 4)
```

Known left nullspace vectors: `{0, 3}` (rows 0 and 3 are identical), `{4, 5}` (rows 4 and 5
are identical), and `{0, 1, 2}` (the XOR of rows 0, 1, 2 is zero). The kernel dimension is 3.

---

## 32. The Nullspace as a Congruence of Squares

### From nullspace vector to congruence

A left nullspace vector v — a subset S of filtered-matrix rows — has a precise algebraic
meaning. Let the rows in S correspond to relations with rational norms N_rat(aᵢ, bᵢ) and
algebraic norms N_alg(aᵢ, bᵢ). The GF(2) column-sum being zero means:

$$
\prod_{i \in S} N_{\mathrm{rat}}(a_i, b_i) = X^2 \quad \text{for some } X \in \mathbb{Z}
$$

$$
\prod_{i \in S} N_{\mathrm{alg}}(a_i, b_i) = Y^2 \quad \text{for some } Y \in K = \mathbb{Q}(\alpha)
$$

The first equation holds because the rational exponent parities all cancel (every prime appears
to an even total exponent). The second holds for the same reason on the algebraic side — with
the additional guarantee provided by the sign and QC columns (§33).

### The congruence of squares

The rational norm N_rat(a, b) = a − b·m is a linear form in (a, b). The product of the
rational norms over S is:

$$
\prod_{i \in S} (a_i - b_i \cdot m) = X^2 \pmod{N}
$$

The algebraic norm N_alg(a, b) = b^d · f(a/b) is the norm of the ideal (a − b·α) in ℤ[α].
The product of the algebraic norms over S is the norm of the product ideal:

$$
\prod_{i \in S} (a_i - b_i \cdot \alpha) = \beta^2 \quad \text{in } K = \mathbb{Q}(\alpha)
$$

Taking the norm of both sides: Norm(β²) = Norm(β)² = Y² where Y = Norm(β) ∈ ℤ. Setting
x = X mod N and y = Y mod N gives x² ≡ y² (mod N). If x ≢ ±y (mod N), then gcd(x − y, N)
is a non-trivial factor of N.

### Why the matrix encodes the right thing

Each row of the matrix is the GF(2) exponent vector of a relation: the parity of each prime's
exponent in the norm factorisation. The GF(2) column-sum of a subset S being zero is exactly
the condition that every prime appears to an even total exponent — i.e., the product of norms
is a perfect square. The linear algebra step is the bridge between the sieve output (individual
smooth pairs) and the factoring step (a congruence of squares).

---

## 33. Why the Sign and Quadratic-Character Columns Are Needed

### The sign column

The rational norm N_rat(a, b) = a − b·m can be negative (when a < b·m). A perfect square must
be positive, so the product of the rational norms over S must be positive. The sign column
(at `obstruction_col_start`) records the parity of the sign: it is 1 if the rational norm is
negative, 0 if positive. For the product to be a positive perfect square, the total number of
negative norms in S must be even — i.e., the sign column must sum to zero over S.

The `Relation` type records this directly:

```rust
pub struct Relation {
    // ...
    /// True if the rational norm a − b·m is negative.
    pub rational_sign: bool,
}
```

`build_matrix` populates the sign column from `relation.rational_sign`. A nullspace vector
that satisfies the sign column constraint has an even number of negative rational norms, so
their product is positive.

### The quadratic-character columns

The algebraic side requires more care. Even if the product of algebraic norms is a perfect
square integer, the algebraic square root β (satisfying Norm(β) = Y) might not exist in K.
The issue is that the product of the ideals (aᵢ − bᵢ·α) might be a square in the ideal group
of ℤ[α] without being a square of a principal ideal — i.e., the product might not be a
principal ideal generated by a square element.

The quadratic-character (QC) columns resolve this. For each auxiliary prime q that splits
completely in K = ℚ[x]/(f), the QC column records the Legendre-symbol parity of the algebraic
norm at q:

$$
\text{QC}_{q}(a, b) = \left(\frac{N_{\mathrm{alg}}(a, b)}{q}\right) \bmod 2
$$

where (·/q) is the Legendre symbol. A nullspace vector that satisfies all QC column constraints
guarantees that the product of the algebraic norms is a square in K, not merely a square
integer. Without QC columns, the algebraic square root might not exist in K even if the norm
is a square.

### The `populate_qc_columns` function

G.E fills in the QC columns before running the nullspace solver:

```rust
pub fn populate_qc_columns(
    matrix: &mut SparseMatrix,
    relations: &[Relation],
    fb: &FactorBase,
    poly: &PolyPair,
    qc_primes: &[u64],
) {
    for row in matrix.rows.iter_mut() {
        for (k, &q) in qc_primes.iter().enumerate() {
            let col = matrix.obstruction_col_start + 1 + k;
            let mut parity = false;
            for &rel_idx in &row.provenance {
                let rel = &relations[rel_idx];
                let norm = algebraic_norm_mod_q(&rel.a, &rel.b, &poly.f, q);
                parity ^= legendre_parity(norm, q);
            }
            if parity {
                // Insert col into row.cols (maintaining sorted order).
                match row.cols.binary_search(&col) {
                    Ok(_) => {}
                    Err(pos) => { row.cols.insert(pos, col); }
                }
            }
        }
    }
}
```

The QC primes are selected by `select_qc_primes`: the first `num_qc` primes q > B_alg such
that f splits completely mod q (i.e., f has deg(f) distinct roots mod q). Splitting completely
ensures the Legendre symbol is well-defined and non-degenerate.

The Legendre symbol is computed via Euler's criterion: (a/q) ≡ a^((q−1)/2) (mod q). The
parity is 1 iff the symbol is −1 (a is a quadratic non-residue mod q):

```rust
fn legendre_parity(norm: u64, q: u64) -> bool {
    if norm == 0 { return false; }
    let exp = (q - 1) / 2;
    let result = pow_mod(norm, exp, q);
    result == q - 1
}
```

**Principle-4 annotation.** At toy scale, `DEFAULT_NUM_QC = 10` QC columns suffice. At NFS
scale, 20–50 QC columns are typical. The number of QC columns is a scale knob: more columns
give a stronger guarantee that the algebraic square root exists, at the cost of a wider matrix.

---

## 34. The C-LinAlg Substrate

### The three frozen types

The C-LinAlg contract (frozen at G.E.1) defines three types that all downstream consumers
(G.E.2, G.E.3, G.F, D.B) use:

**`BlockVec`** — a block of `BLOCK_WIDTH = 64` GF(2) vectors, each of length `num_rows`,
packed into a `Vec<u64>`. The layout is "row of words": `data[i]` is a `u64` whose bit j is
the j-th vector's value at row i. This makes iterating over rows contiguous in memory and
keeps the 64 vectors interleaved bit-by-bit within each word.

```rust
pub struct BlockVec {
    /// Packed row data: `data[row]` bit `j` = vector `j`'s value at `row`.
    pub data: Vec<u64>,
    /// Number of rows (the vector dimension).
    pub num_rows: usize,
}
```

**`MatrixOperator`** — a view of a `SparseMatrix` as a linear operator over GF(2). Provides
`apply` (A·V) and `apply_transpose` (Aᵀ·V) for block vectors. Both solvers consume this
interface exclusively; they never read `SparseMatrix` fields directly.

```rust
pub struct MatrixOperator<'a> {
    matrix: &'a SparseMatrix,
}
```

**`KernelVector`** — a vector in the left nullspace of the matrix, represented as a sorted,
deduplicated `Vec<usize>` of filtered-matrix row indices.

```rust
pub struct KernelVector {
    /// Sorted, deduplicated row indices into the filtered matrix.
    pub row_indices: Vec<usize>,
}
```

### The `MatrixOperator` apply operations

The `apply` operation (A·V) computes the GF(2) matrix-vector product for all 64 block vectors
simultaneously. For each row i of the matrix, the result is the XOR of `v.data[c]` for each
column c in that row:

```rust
pub fn apply(&self, v: &BlockVec) -> BlockVec {
    let mut result = BlockVec::zeros(self.num_rows());
    for (i, row) in self.matrix.rows.iter().enumerate() {
        let mut word = 0u64;
        for &c in &row.cols {
            word ^= v.data[c];
        }
        result.data[i] = word;
    }
    result
}
```

The key insight: `v.data[c]` is a `u64` whose bit j is vector j's value at column c. XOR-ing
all `v.data[c]` for c in the row gives a `u64` whose bit j is the GF(2) dot product of the
row with vector j. All 64 dot products are computed in a single XOR loop.

The `apply_transpose` operation (Aᵀ·V) scatters contributions on-the-fly: for each row i,
for each column c in that row, XOR `v.data[i]` into `result.data[c]`:

```rust
pub fn apply_transpose(&self, v: &BlockVec) -> BlockVec {
    let mut result = BlockVec::zeros(self.num_cols());
    for (i, row) in self.matrix.rows.iter().enumerate() {
        let vi = v.data[i];
        if vi == 0 { continue; }
        for &c in &row.cols {
            result.data[c] ^= vi;
        }
    }
    result
}
```

### The `BLOCK_WIDTH` scale knob

```rust
/// Block width: 64 vectors packed into machine words.
///
/// Principle-4 annotation: at toy scale a single word suffices and the blocking overhead
/// is invisible; at NFS scale the word-wide block is the inner loop's cache-friendly unit.
/// The width is the scale knob — D.B may widen to 128 or parameterise over block width.
pub const BLOCK_WIDTH: usize = 64;
```

At toy scale (6 × 4 matrix), a single 64-bit word is more than enough to represent the entire
block. The blocking overhead is invisible. At NFS scale (millions of rows), the 64-wide block
is the cache-friendly inner loop unit: each `apply` call processes 64 vectors in a single pass
over the matrix, amortising the memory-bandwidth cost of loading each row once.

---

## 35. Block Lanczos: Krylov Subspaces over GF(2)

### The symmetric reduction

Both solvers work with the symmetric matrix B = A · Aᵀ (m × m), whose left nullspace equals
the left nullspace of A. This is because:

$$
A^T v = 0 \implies A A^T v = 0
$$

and conversely, if B·v = 0 then A·Aᵀ·v = 0, so Aᵀ·v is in the nullspace of A, and if A has
full column rank then Aᵀ·v = 0. Working with B avoids the need to handle the non-square
structure of A directly.

The B·v product is computed as two operator applications:

```rust
// w_cur = B * v_cur = A * (A^T * v_cur).
let at_v = op.apply_transpose(&v_cur);
let w_cur = op.apply(&at_v);
```

### The Krylov subspace

The block Lanczos algorithm builds a Krylov basis for B. Starting from a random block vector
V₀ ∈ GF(2)^{m × 64}, the Krylov sequence is:

$$
V_0,\; B V_0,\; B^2 V_0,\; \ldots
$$

The algorithm constructs an orthogonal basis for the Krylov subspace span{V₀, BV₀, B²V₀, …}
using a three-term recurrence. At each step, the recurrence advances the block vector using
the A-orthogonality condition.

### The A-orthogonality condition and the three-term recurrence

The block Lanczos recurrence is:

$$
V_{k+1} = B V_k - V_k \alpha_k - V_{k-1} \beta_k
$$

where the coefficient matrices α_k and β_k are chosen so that V_{k+1} is A-orthogonal to all
previous block vectors. The A-orthogonality condition is:

$$
V_i^T B V_j = 0 \quad \text{for } i \neq j
$$

The coefficients are:

$$
\alpha_k = S_k^{-1} (V_k^T B^2 V_k), \quad \beta_k = S_{k-1}^{-1} (V_{k-1}^T B V_k)
$$

where S_k = V_k^T B V_k is the BLOCK_WIDTH × BLOCK_WIDTH inner product matrix. The inverse
S_k^{-1} is computed by GF(2) Gaussian elimination.

In the implementation:

```rust
// alpha_raw = v_cur^T * B * w_cur (BLOCK_WIDTH × BLOCK_WIDTH).
let alpha_raw = v_cur.inner_product_matrix(&bw_cur);
// alpha = s_inv * alpha_raw (restricted to active columns).
let alpha = gf2_matmul_block(s_inv, alpha_raw, active_mask);

// beta_raw = v_prev^T * w_cur (BLOCK_WIDTH × BLOCK_WIDTH).
let beta_raw = v_prev.inner_product_matrix(&w_cur);
// beta = s_prev_inv * beta_raw (restricted to previously active columns).
let beta = gf2_matmul_block(s_prev_inv, beta_raw, active_mask_prev);

// v_next = w_cur - v_cur * alpha - v_prev * beta (over GF(2), - = +).
let mut v_next = w_cur.clone();
block_vec_sub_matmul(&mut v_next, &v_cur, alpha);
block_vec_sub_matmul(&mut v_next, &v_prev, beta);
```

### The GF(2) self-orthogonality problem and winnowing

Over ℝ, the Lanczos recurrence is guaranteed to produce orthogonal vectors because the inner
product is positive-definite. Over GF(2), this guarantee fails: a **nonzero** vector v can
satisfy v^T B v = 0. This is the self-orthogonality problem.

Self-orthogonality occurs exactly when B·v = 0, i.e., when v is already in the nullspace of B
(and hence in the left nullspace of A). When the inner product matrix S_k = V_k^T B V_k has
a zero column, the corresponding block column of V_k is self-orthogonal — it is already a
nullspace vector.

The **winnowing** step handles this: at each iteration, the algorithm identifies which block
columns are "active" (those where S_k has a non-zero pivot) and which are "inactive" (those
where S_k has a zero column). Inactive columns are checked for nullspace membership:

```rust
// Find active columns via GF(2) Gaussian elimination on s.
let (active_mask, s_inv) = gf2_block_pivot(s);

// Collect kernel candidates from inactive columns.
let inactive_mask = !active_mask;
let mut bits = inactive_mask;
while bits != 0 {
    let j = bits.trailing_zeros() as usize;
    bits &= bits - 1;

    let col = v_cur.column(j);
    if col.iter().any(|&b| b) {
        let col_bv = BlockVec::from_columns(&[col.clone()]);
        let at_col = op.apply_transpose(&col_bv);
        let is_zero = at_col.data.iter().all(|&x| x == 0);
        if is_zero {
            let kv = KernelVector::from_mask(&col);
            if !kv.is_empty() {
                results.push(kv);
            }
        }
    }
}
```

The `gf2_block_pivot` function performs GF(2) Gaussian elimination on the 64 × 64 inner
product matrix S, returning the pivot column mask (active columns) and the inverse of the
pivot submatrix (used in the recurrence):

```rust
fn gf2_block_pivot(s: [u64; BLOCK_WIDTH]) -> (u64, [u64; BLOCK_WIDTH]) {
    let mut s_work = s;
    let mut inv = [0u64; BLOCK_WIDTH];
    for i in 0..BLOCK_WIDTH { inv[i] = 1u64 << i; } // identity
    let mut pivot_mask = 0u64;
    for col in 0..BLOCK_WIDTH {
        // Find pivot row for this column.
        let mut found = BLOCK_WIDTH;
        for row in col..BLOCK_WIDTH {
            if (s_work[row] >> col) & 1 == 1 { found = row; break; }
        }
        if found == BLOCK_WIDTH { continue; } // no pivot — inactive column
        s_work.swap(col, found);
        inv.swap(col, found);
        pivot_mask |= 1u64 << col;
        // Full reduced row echelon form.
        for row in 0..BLOCK_WIDTH {
            if row != col && (s_work[row] >> col) & 1 == 1 {
                s_work[row] ^= s_work[col];
                inv[row] ^= inv[col];
            }
        }
    }
    // ...
    (pivot_mask, s_inv)
}
```

**Principle-4 annotation (self-orthogonality).** Self-orthogonality is a phenomenon that IS
exposed at toy scale. The KAT matrix in `lanczos_kat.rs` is deliberately constructed with
duplicate rows (rows 0 and 3 are identical; rows 4 and 5 are identical). Duplicate rows
produce self-orthogonal block columns in the starting vector, exercising the winnowing path.
The KAT comment explains this explicitly:

> To force this path in a KAT, we use a matrix with duplicate rows. If rows i and j are
> identical, then e_i XOR e_j is in the left nullspace of A. When the random starting block
> vector has a component along e_i XOR e_j, that component is self-orthogonal under B and
> triggers the winnowing.

**Principle-4 annotation (block width).** The block width `BLOCK_WIDTH = 64` is the scale
knob. At toy scale (6 × 4 matrix), the entire matrix fits in a handful of words and the
blocking overhead is invisible — a single-vector Lanczos would be equally fast. At NFS scale
(millions of rows), the 64-wide block amortises the cost of loading each matrix row once
across 64 simultaneous vector operations, giving a ~64× speedup over the single-vector variant.

---

## 36. Block Wiedemann: Krylov Sequences and Berlekamp-Massey

### The Wiedemann approach

Block Wiedemann takes a different route to the same destination. Instead of building an
orthogonal Krylov basis (Lanczos), it computes a scalar Krylov sequence and uses
Berlekamp-Massey to find the minimal polynomial of B restricted to a random direction.

For a random pair (x, y) ∈ GF(2)^m × GF(2)^m, the scalar Krylov sequence is:

$$
s_i = x^T B^i y, \quad i = 0, 1, 2, \ldots, 2m + 10
$$

The Berlekamp-Massey algorithm finds the shortest LFSR that generates this sequence — the
minimal polynomial f(z) of B in the direction (x, y). Then the kernel vector is extracted by
evaluating f(B)·y:

$$
w = f(B) \cdot y = \sum_{k=0}^{d} f_k \cdot B^k \cdot y
$$

If A^T·w = 0 and w ≠ 0, then w is a left nullspace vector of A.

### The Krylov sequence computation

```rust
fn krylov_sequence(
    op: &MatrixOperator<'_>,
    m: usize,
    x: &[bool],
    y: &[bool],
    seq_len: usize,
) -> Vec<bool> {
    let y_bv = bool_vec_to_blockvec(y, m);
    let mut v = y_bv;
    let mut sequence = Vec::with_capacity(seq_len);

    for _ in 0..seq_len {
        // s_i = x^T * v = XOR of v[j] for j where x[j] = true.
        let s = inner_product_scalar(x, &v);
        sequence.push(s);

        // v = B * v = apply(apply_transpose(v)).
        let at_v = op.apply_transpose(&v);
        v = op.apply(&at_v);
    }
    sequence
}
```

Each step applies B = A·Aᵀ to the current vector, then takes the scalar inner product with x.
The sequence length 2m + 10 is sufficient for Berlekamp-Massey to find the minimal polynomial
of degree at most m.

### The Berlekamp-Massey algorithm

Berlekamp-Massey finds the shortest LFSR (linear feedback shift register) that generates a
given GF(2) sequence. The output is the minimal polynomial f(z) = 1 + f₁z + … + f_d z^d such
that:

$$
\sum_{k=0}^{d} f_k \cdot s_{n-k} = 0 \quad \text{for all } n \geq d
$$

The implementation follows the standard GF(2) Berlekamp-Massey algorithm:

```rust
pub fn berlekamp_massey(s: &[bool]) -> Vec<bool> {
    let n = s.len();
    let mut c = vec![true]; // current LFSR polynomial C = 1
    let mut b = vec![true]; // previous polynomial B = 1
    let mut l: usize = 0;   // current LFSR length
    let mut m: usize = 1;   // steps since last length change

    for n_idx in 0..n {
        // Compute discrepancy d = s[n_idx] XOR sum_{i=1}^{L} C[i] * s[n_idx - i].
        let mut d = s[n_idx];
        for i in 1..=l {
            if n_idx >= i && c.len() > i && c[i] { d ^= s[n_idx - i]; }
        }

        if !d {
            m += 1; // no update needed
        } else if 2 * l <= n_idx {
            // Length must increase.
            let t = c.clone();
            xor_shifted(&mut c, &b, m); // C = C XOR z^m * B
            l = n_idx + 1 - l;
            b = t;
            m = 1;
        } else {
            xor_shifted(&mut c, &b, m); // C = C XOR z^m * B
            m += 1;
        }
    }
    c
}
```

The KAT in `wiedemann_kat.rs` verifies this on the Fibonacci sequence mod 2 (0, 1, 1, 0, 1,
1, …), which has minimal polynomial f(z) = 1 + z + z² (degree 2):

```rust
let s: Vec<bool> = vec![false, true, true, false, true, true, false, true, true, false];
let f = berlekamp_massey(&s);
assert_eq!(f.len(), 3); // degree 2
assert!(f[0] && f[1] && f[2]); // f(z) = 1 + z + z^2
```

### Kernel extraction via Horner's method

Once the minimal polynomial f(z) is known, the kernel vector is extracted by evaluating
f(B)·y via Horner's method:

$$
f(B) \cdot y = f_d \cdot y + B \cdot (f_{d-1} \cdot y + B \cdot (\ldots + B \cdot (f_1 \cdot y + B \cdot (f_0 \cdot y)) \ldots))
$$

```rust
fn eval_poly_on_krylov(op: &MatrixOperator<'_>, m: usize, y: &[bool], f: &[bool]) -> BlockVec {
    let d = f.len() - 1;
    let y_bv = bool_vec_to_blockvec(y, m);
    let mut w = if f[d] { y_bv.clone() } else { BlockVec::zeros(m) };

    for k in (0..d).rev() {
        // w = B * w.
        let at_w = op.apply_transpose(&w);
        w = op.apply(&at_w);
        // w = w XOR (f_k * y).
        if f[k] { w.xor_assign(&y_bv); }
    }
    w
}
```

The result w = f(B)·y satisfies B·w = 0 (since f is the minimal polynomial of B in the
direction (x, y)), and hence A^T·w = 0 — w is in the left nullspace of A.

### The parallelism payoff

The key architectural difference between Wiedemann and Lanczos is where the parallelism lives.
In block Lanczos, each iteration requires a global inner product (S_k = V_k^T B V_k) that
synchronises all 64 block columns. In block Wiedemann, the Krylov sequence computation for
each (x, y) pair is independent: the sequence {x^T B^i y} can be computed in parallel across
multiple (x, y) pairs with no global synchronisation per step.

At production NFS scale, the block Wiedemann algorithm uses `BLOCK_WIDTH` (x, y) pairs
simultaneously, distributing the Krylov sequence computation across multiple machines. Each
machine computes a subset of the sequences independently; the Berlekamp-Massey step is then
run on the combined sequences. This is the architecture used in the RSA-768 factorisation
(Kleinjung et al., 2010).

**Principle-4 annotation (Wiedemann parallelism).** At toy scale, the parallelism payoff is
invisible. The implementation uses a single (x, y) pair (the scalar Wiedemann variant) at
demonstration fidelity; the block variant would use `BLOCK_WIDTH` pairs simultaneously. The
module docstring annotates this explicitly:

> Block Wiedemann's payoff is distributed/parallel: the Krylov sequence {x^T A^i y} can be
> computed in parallel across multiple (x, y) pairs, with no global synchronisation per step
> (unlike block Lanczos, which requires a global inner product at each step). At toy scale,
> this parallelism is invisible — Lanczos is simpler and just as fast.

The implementation runs 4 independent (x, y) attempts to improve the probability of finding
a kernel vector:

```rust
let num_attempts = 4;
for attempt in 0..num_attempts {
    // ...
    let x = random_gf2_vec(m, rng_state);
    let y = random_gf2_vec(m, rng_state);
    if let Some(kv) = wiedemann_attempt(op, m, &x, &y) {
        if !results.iter().any(|r| r.row_indices == kv.row_indices) {
            results.push(kv);
        }
    }
}
```

---

## 37. The Kernel-Vector → Provenance → Original-Relations Thread

### What G.F consumes

G.F (the square-root step) needs the original (a, b) pairs — not the filtered-matrix row
indices. The bridge is the provenance map: each `MatrixRow` carries a sorted, deduplicated
list of original relation indices that merged into it during G.D filtering.

The `KernelVector` type provides the `expand_provenance` method that performs this expansion:

```rust
pub fn expand_provenance(&self, matrix: &SparseMatrix) -> Vec<usize> {
    let mut result: Vec<usize> = Vec::new();
    for &i in &self.row_indices {
        let prov = &matrix.rows[i].provenance;
        // Symmetric difference of result and prov (both sorted).
        let mut merged = Vec::with_capacity(result.len() + prov.len());
        let mut ri = 0;
        let mut pi = 0;
        while ri < result.len() && pi < prov.len() {
            match result[ri].cmp(&prov[pi]) {
                std::cmp::Ordering::Less    => { merged.push(result[ri]); ri += 1; }
                std::cmp::Ordering::Greater => { merged.push(prov[pi]); pi += 1; }
                std::cmp::Ordering::Equal   => { ri += 1; pi += 1; } // cancel
            }
        }
        merged.extend_from_slice(&result[ri..]);
        merged.extend_from_slice(&prov[pi..]);
        result = merged;
    }
    result
}
```

The symmetric difference (XOR union) is the correct operation: if a relation index appears in
two provenance sets, it means the relation was merged into two different filtered rows, and
those two rows' contributions cancel in the GF(2) sum. The symmetric difference removes
cancelled indices, leaving only the original relations that contribute a net odd number of
times.

### The full thread

The complete kernel-vector → original-relations thread:

1. **G.E** runs `block_lanczos` or `block_wiedemann` on the filtered matrix, returning a
   `Vec<KernelVector>`. Each `KernelVector` is a sorted list of filtered-matrix row indices.

2. **G.F** calls `kv.expand_provenance(&matrix)` to get the set of original relation indices.
   This is the symmetric difference of the provenance sets of all selected rows.

3. **G.F** looks up the original `Relation` objects by index to recover the (a, b) pairs.

4. **G.F** computes the product of all rational norms (a perfect square in ℤ) and the product
   of all algebraic norms (a perfect square in K = ℚ(α)), and extracts the square roots to
   form the congruence x² ≡ y² (mod N).

The KAT `kat_3_round_trip_provenance` in `linalg_substrate_kat.rs` verifies this thread end-
to-end before G.F exists. For a matrix with hand-crafted provenance:

```rust
// Row 0: provenance = [0, 1], Row 1: provenance = [2], Row 2: provenance = [1, 3].
// Kernel vector {0, 2}: sym_diff([0,1], [1,3]) = [0, 3] (1 cancels).
let kv_02 = KernelVector::new(vec![0, 2]);
let expanded = kv_02.expand_provenance(&matrix);
assert_eq!(expanded, vec![0, 3]);

// Kernel vector {0, 1, 2} (valid nullspace vector):
// sym_diff([0,1], [2], [1,3]) = [0, 2, 3] (1 cancels from rows 0 and 2).
let kv_012 = KernelVector::new(vec![0, 1, 2]);
assert!(kv_012.verify(&matrix));
let expanded = kv_012.expand_provenance(&matrix);
assert_eq!(expanded, vec![0, 2, 3]);
```

### Why row indices, not a bit-mask

The `KernelVector` stores row indices rather than a bit-mask for three reasons:

1. G.F needs row indices to look up provenance; a bit-mask would require a scan.
2. Kernel vectors are sparse (typically a small fraction of rows); a bit-mask wastes space.
3. Solvers internally work with bit-packed block vectors, but convert to `KernelVector` on
   output — the conversion is O(rows) and happens once per kernel vector, not in the inner
   loop.

---

## 38. The C-LinAlg Contract and D.B Generalisation

### The frozen seam

The C-LinAlg contract (frozen at G.E.1) defines the interface that all downstream consumers
use. The three types — `BlockVec`, `MatrixOperator`, and `KernelVector` — are the frozen seam.
G.E.2 (block Lanczos), G.E.3 (block Wiedemann), G.F (square root), and D.B (NFS-DL linear
algebra over F_ℓ) all consume this interface directly.

The module docstring in `gnfs/src/linalg/mod.rs` states the contract explicitly:

```rust
//! # C-LinAlg contract
//!
//! The types and functions in this module implement the C-LinAlg contract frozen at G.E.1.
//! G.E.2 (block Lanczos), G.E.3 (Wiedemann), G.F (square root), and D.B (GF(ℓ) extension)
//! consume this interface directly.
```

### The D.B generalisation

D.B (NFS-DL linear algebra) generalises the G.E linear algebra from GF(2) to F_ℓ (the field
with ℓ elements, where ℓ is a prime). The generalisation touches three points:

**`BlockVec` generalisation.** The GF(2)-specific packing (64 scalars per `u64`) is the scale
knob. For F_ℓ, the natural generalisation is `data: Vec<[Scalar; BLOCK_WIDTH]>` where `Scalar`
is the field element type. The `BlockVec` docstring annotates this:

> For F_ℓ (ℓ > 2), the natural generalisation is `data: Vec<[Scalar; BLOCK_WIDTH]>` where
> `Scalar` is the field element type. The GF(2) specialisation packs 64 scalars into one
> `u64`. D.B may introduce a `BlockVec<S>` generic or a parallel `BlockVecFl` type; the
> *interface* (inner products, apply, apply_transpose) is the stable seam.

**`MatrixOperator` generalisation.** The operator interface is already abstract: `apply` and
`apply_transpose` take and return `BlockVec`. For F_ℓ, the same interface shape applies with
scalar multiplication replacing the GF(2) XOR. The `MatrixOperator` docstring annotates:

> For F_ℓ, the operator needs the same shape but with scalar multiplication. The natural
> generalisation is a trait `LinearOperator<V>` with `apply(&self, v: &V) -> V` and
> `apply_transpose(&self, v: &V) -> V`. G.E implements the concrete GF(2) version; D.B
> may introduce the trait and have `MatrixOperator` implement it.

**`KernelVector` generalisation.** For F_ℓ, a kernel vector is still a subset of rows (those
with nonzero coefficient in the nullspace vector). The row-index spine is stable; D.B may add
a `coefficients` field (`Vec<Scalar>`) for the non-GF(2) case.

### The interface shape is stable

The key design decision is that the GF(2)-specific packing is isolated in `BlockVec::data`
(the `Vec<u64>` representation) and the inner-loop operations (`xor_assign`,
`inner_product_matrix`). The interface shape — `apply`, `apply_transpose`, `KernelVector`,
`expand_provenance` — is stable across the GF(2) → F_ℓ generalisation. D.B can reuse the
`MatrixOperator` and `KernelVector` types directly, or introduce a thin generic wrapper, without
changing the solver structure.

---

## 39. Principle-4 Annotations Summary (G.E)

The following table collects all principle-4 annotations from the G.E module — phenomena that
are scale-dependent and either under-exposed or over-exposed at toy scale.

| Phenomenon | Toy-scale behaviour | NFS-scale behaviour | Scale knob |
|------------|--------------------|--------------------|------------|
| **Block width** | Single 64-bit word suffices; blocking overhead invisible | 64-wide block amortises memory bandwidth; ~64× speedup over single-vector | `BLOCK_WIDTH = 64` |
| **Self-orthogonality** | Exposed by duplicate rows in KAT matrices; winnowing path exercised | Same algorithm; more frequent at scale due to larger nullspace | — (algorithm feature) |
| **Wiedemann parallelism** | Invisible; 4 sequential attempts, no distributed computation | Krylov sequences computed in parallel across machines; no per-step sync | `num_attempts` / block width |
| **QC column count** | `DEFAULT_NUM_QC = 10` suffices | 20–50 QC columns typical | `DEFAULT_NUM_QC` |
| **Excess floor** | `EXCESS_FLOOR = 20` never approached in practice | Floor is a tuning parameter: too low → few null-space vectors; too high → denser matrix | `EXCESS_FLOOR` |

The self-orthogonality phenomenon is the one that IS fully exposed at toy scale: the KAT
matrices are deliberately constructed to trigger the winnowing path, and the algorithm handles
it correctly. The block-width and Wiedemann-parallelism phenomena are the ones most
under-exposed: at toy scale, both solvers are instantaneous and the performance difference
between them is invisible.

---

## 40. KAT Summary (G.E)

The following table lists the key known-answer tests across the linalg module, with the
mathematical fact each one verifies.

| Test file | Test | Mathematical fact verified |
|-----------|------|---------------------------|
| `linalg_substrate_kat.rs` | `kat_1_operator_correctness` | `A·V` and `Aᵀ·V` match hand-computed GF(2) products for a 3×4 matrix |
| `linalg_substrate_kat.rs` | `kat_2_qc_column_construction` | QC column parities match hand-computed Legendre symbols for toy relations |
| `linalg_substrate_kat.rs` | `kat_3_round_trip_provenance` | `expand_provenance` returns the correct symmetric difference of provenance sets |
| `linalg_substrate_kat.rs` | `kat_4_determinism` | Operator products and QC columns are deterministic for a fixed matrix |
| `lanczos_kat.rs` | `kat_a_correctness_with_self_orthogonality` | Block Lanczos finds valid kernel vectors for a 6×4 matrix with duplicate rows (self-orthogonality path exercised) |
| `lanczos_kat.rs` | `kat_a2_single_dependency` | Block Lanczos finds the single dependency `{0,3}` for a 4×3 matrix |
| `lanczos_kat.rs` | `kat_a3_full_rank_no_nullspace` | Block Lanczos returns only valid vectors for a full-rank 3×3 identity matrix |
| `lanczos_kat.rs` | `kat_b_determinism` | Block Lanczos is deterministic for a fixed matrix and seed |
| `lanczos_kat.rs` | `kat_c_cado_oracle_n35` | (Ignored) CADO oracle for N=35 — gated when CADO absent |
| `lanczos_kat.rs` | `kat_multiple_dependencies` | Block Lanczos finds at least one of three known dependencies in an 8×6 matrix |
| `wiedemann_kat.rs` | `kat_a_cross_validation_with_lanczos` | Wiedemann and Lanczos cross-validate: both find valid kernel vectors for the shared 6×4 matrix |
| `wiedemann_kat.rs` | `kat_a4_single_dependency` | Wiedemann finds the single dependency `{0,3}` for a 4×3 matrix |
| `wiedemann_kat.rs` | `kat_b_determinism` | Wiedemann is deterministic for a fixed matrix and seed |
| `wiedemann_kat.rs` | `kat_c_bm_fibonacci_degree` | Berlekamp-Massey returns degree-2 polynomial `1 + z + z²` for the Fibonacci mod 2 sequence |
| `wiedemann_kat.rs` | `kat_c2_bm_period4_sequence` | Berlekamp-Massey returns degree-4 polynomial `1 + z⁴` for a period-4 sequence |

---

## 41. What Comes Next

Linear algebra is the fifth stage of the GNFS pipeline. The output — a `Vec<KernelVector>`,
each with its provenance expansion thread back to the original relations — is the input to the
square-root stage.

**G.F (square root)** takes a `KernelVector` from G.E and expands it through the provenance
map to recover the original (a, b) pairs. It then computes the product of all rational norms
(a perfect square in ℤ) and the product of all algebraic norms (a perfect square in K = ℚ(α)),
and extracts the square roots to form the congruence x² ≡ y² (mod N). A non-trivial GCD of
x − y and N yields a factor.

**D.B (NFS-DL linear algebra)** generalises the G.E linear algebra from GF(2) to F_ℓ. The
`BlockVec`, `MatrixOperator`, and `KernelVector` types are the frozen seam; D.B may introduce
a generic wrapper or a parallel `BlockVecFl` type for the F_ℓ case. The solver structure
(block Lanczos or block Wiedemann) is unchanged; only the field arithmetic differs.

---

## Further Reading (G.E)

1. **Montgomery, P. L. (1995).** *A block Lanczos algorithm for finding dependencies over
   GF(2).* In: Guillou, L. C., and Quisquater, J.-J. (eds.) Advances in Cryptology —
   EUROCRYPT '95, LNCS 921. Springer. The original source for the block Lanczos algorithm
   used in G.E.2, including the self-orthogonality winnowing and the three-term recurrence.

2. **Coppersmith, D. (1994).** *Solving homogeneous linear equations over GF(2) via block
   Wiedemann algorithm.* Mathematics of Computation, 62(205), 333–350. The original source
   for the block Wiedemann algorithm used in G.E.3, including the Berlekamp-Massey step and
   the parallelism analysis.

3. **Wiedemann, D. (1986).** *Solving sparse linear equations over finite fields.* IEEE
   Transactions on Information Theory, 32(1), 54–62. The scalar Wiedemann algorithm that
   block Wiedemann generalises.

4. **Kleinjung, T., Aoki, K., Franke, J., Lenstra, A. K., Thomé, E., Bos, J. W., Gaudry, P.,
   Kruppa, A., Montgomery, P. L., Osvik, D. A., te Riele, H., Timofeev, A., and Zimmermann,
   P. (2010).** *Factorization of a 768-bit RSA modulus.* In: Rabin, T. (ed.) Advances in
   Cryptology — CRYPTO 2010, LNCS 6223. Springer. The RSA-768 factorisation, which used
   block Wiedemann for the linear algebra step across a distributed cluster.

5. **Berlekamp, E. R. (1968).** *Algebraic Coding Theory.* McGraw-Hill. The original source
   for the Berlekamp-Massey algorithm used in the Wiedemann step.

6. **Massey, J. L. (1969).** *Shift-register synthesis and BCH decoding.* IEEE Transactions
    on Information Theory, 15(1), 122–127. The Massey formulation of the LFSR synthesis
    problem, which is the GF(2) Berlekamp-Massey algorithm used in `wiedemann.rs`.

---

# Square Root and Assembly: A Code-Tour Chapter

This chapter explains the `gnfs/src/sqrt` module — the square-root and assembly stage of the
General Number Field Sieve. It is organised by the implementation, not by mathematical
abstraction: each section introduces the mathematics, then shows how it is realised in code.
A reader who has read the linear algebra chapter (§31–§41 above) can read this one without
consulting the source.

The chapter covers the kernel vector as a congruence of squares, the rational square root, the
algebraic square root via Couveignes' CRT algorithm, the final GCD assembly, and the whole
factoring arc that this stage closes.

---

## 42. The Kernel Vector as a Congruence of Squares

### What a `KernelVector` represents

The linear algebra step (G.E, §31–§41) produces a `KernelVector` — a sorted list of
filtered-matrix row indices whose GF(2) column-sum is the zero vector. Each row in the
filtered matrix corresponds to a relation (or a merged aggregate of relations), and each
column corresponds to a prime or obstruction entry in the factor base.

The GF(2) column-sum being zero has a precise algebraic meaning. Let S be the set of original
relation indices recovered from the kernel vector via `expand_provenance` (§37). For every
prime p in the factor base, the total exponent of p across all relations in S is even. This
means:

$$
\prod_{i \in S} (a_i - b_i \cdot m) = X^2 \quad \text{for some } X \in \mathbb{Z}
$$

$$
\prod_{i \in S} (a_i - b_i \cdot \alpha) = \beta^2 \quad \text{for some } \beta \in K = \mathbb{Q}(\alpha)
$$

The first equation holds because the rational exponent parities all cancel (every prime appears
to an even total exponent). The second holds for the same reason on the algebraic side — with
the additional guarantee provided by the sign and quadratic-character (QC) columns.

### Why the QC columns guarantee β exists in K

The QC columns (populated by `populate_qc_columns` in G.E, §33) are the load-bearing
guarantee that the algebraic square root β exists *in K*, not merely as an ideal square. For
each auxiliary prime q that splits completely in K, the QC column records the Legendre-symbol
parity of the algebraic norm at q:

$$
\text{QC}_{q}(a, b) = \left(\frac{N_{\mathrm{alg}}(a, b)}{q}\right) \bmod 2
$$

A kernel vector that satisfies all QC column constraints guarantees that the product of the
algebraic norms is a square in K (a principal ideal generated by a square element), not merely
a square integer. Without QC columns, the algebraic square root might not exist in K even if
the norm is a square.

The QC primes are selected by `select_qc_primes` in `gnfs/src/linalg/qc.rs`: the first
`DEFAULT_NUM_QC = 10` primes q > B_alg such that f splits completely mod q. The same prime
selection machinery is reused by Couveignes' algorithm (§44).

### The `expand_provenance` seam

The `KernelVector::expand_provenance` method (§37) bridges the filtered-matrix row indices to
the original relation indices. It computes the symmetric difference of the provenance sets of
all selected rows:

```rust
pub fn expand_provenance(&self, matrix: &SparseMatrix) -> Vec<usize> {
    // For each selected row, take the symmetric difference of its provenance set
    // with the accumulated result. Indices that appear an even number of times cancel.
    // ...
}
```

This seam was over-specified at G.E.1 precisely for this consumer. G.F is its first real
client: the recovered index set S is the input to both the rational and algebraic square root
computations.

---

## 43. The Rational Square Root

### The algorithm

The rational square root is implemented in `gnfs/src/sqrt/rational.rs`. Given a
`KernelVector`, it computes X such that $X^2 \equiv \prod_{i \in S}(a_i - b_i \cdot m) \pmod{N}$.

The algorithm has five steps:

1. **Recover S.** Call `kv.expand_provenance(matrix)` to get the relation index set S.
2. **Form the product.** Compute $P = \prod_{i \in S}(a_i - b_i \cdot m)$ over ℤ.
3. **Assert P > 0.** The sign column in G.E guarantees an even count of negative rational
   norms, so P must be positive. A negative P is an upstream kernel/sign-column bug.
4. **Extract the integer square root.** Call `isqrt(&P)` from `shared_bigint`. If `None`,
   panic: P is not a perfect square, indicating an upstream bug.
5. **Reduce mod N.** Return X = isqrt(P) mod N.

```rust
pub fn rational_sqrt(
    kv: &KernelVector,
    matrix: &SparseMatrix,
    relations: &[Relation],
    poly: &PolyPair,
) -> BigInt {
    let s = kv.expand_provenance(matrix);

    let m = &poly.m;
    let mut product = BigInt::from(1i64);
    for &i in &s {
        let rel = &relations[i];
        let factor = &rel.a - &rel.b * m;
        product *= factor;
    }

    if product.is_negative() {
        panic!("rational_sqrt: rational norm product is negative — upstream kernel/sign-column bug");
    }

    let x = isqrt(&product).unwrap_or_else(|| {
        panic!("rational_sqrt: rational norm product is not a perfect square — upstream bug")
    });

    let x_mod_n = x % &poly.n;
    if x_mod_n.is_negative() { x_mod_n + &poly.n } else { x_mod_n }
}
```

### Why the product is a perfect square

The matrix nullspace + sign column guarantee it. The GF(2) column-sum being zero means every
prime appears to an even total exponent across all relations in S. The sign column (at
`obstruction_col_start`) ensures the count of negative rational norms is even, so the product
is positive. Together, these two conditions guarantee that P is a positive perfect square.

The `isqrt` function in `shared_bigint` computes the exact integer square root: `isqrt(n²) =
Some(n)`, `isqrt(k) = None` if k is not a perfect square. A `None` return from `isqrt` is
therefore a detected upstream bug, not a normal path.

### The sign invariant

The `rational_sign: bool` field on each `Relation` records whether $a - b \cdot m < 0$. The
G.E linear algebra step selects only subsets S where the count of negative rational norms is
even, so the product is guaranteed positive. The `rational_sqrt` function checks this
defensively and panics loudly if violated — surfacing the upstream bug rather than papering
over it with an incorrect result.

---

## 44. Couveignes' Algorithm: The Algebraic Square Root

### The problem

The algebraic square root is the harder of the two computations. We need $\beta \in K =
\mathbb{Q}(\alpha)$ such that:

$$
\beta^2 = \prod_{i \in S}(a_i - b_i \cdot \alpha) = \gamma \in K
$$

Then $Y = |\mathrm{Norm}(\beta)| \bmod N$ is the algebraic square root. The QC columns
guarantee $\beta$ exists in K; the challenge is to compute it.

Couveignes' algorithm solves this by Chinese remaindering in $\mathbb{Z}[\alpha]$: reduce
$\gamma$ modulo many split prime ideals, take the square root in each residue field
$\mathbb{F}_p$ (Tonelli–Shanks), and CRT-lift the per-prime results to recover $\beta$.

The algorithm is implemented in `gnfs/src/sqrt/algebraic.rs`.

### Step 1: Form γ

The first step forms $\gamma = \prod_{i \in S}(a_i - b_i \cdot \alpha) \in K$ using
`NumberFieldElement::mul`:

```rust
fn form_gamma<'a>(nf: &'a NumberField, relations: &[Relation], s: &[usize])
    -> NumberFieldElement<'a>
{
    let mut gamma = nf.from_int(BigInt::from(1i64));
    for &i in s {
        let rel = &relations[i];
        // Factor: a_i − b_i·α = a_i·1 + (−b_i)·α.
        let factor_poly = RatPoly::from_coeffs(vec![
            BigRational::from(rel.a.clone()),
            BigRational::from(-rel.b.clone()),
        ]);
        let factor = NumberFieldElement { field: nf, poly: factor_poly };
        gamma = gamma.mul(&factor);
    }
    gamma
}
```

Each factor $(a_i - b_i \cdot \alpha)$ is a degree-1 element of K. The product $\gamma$ is
an element of $\mathbb{Z}[\alpha]$ of degree $< d$ (reduced modulo the defining polynomial f).

### Step 2: Select CRT primes

Couveignes' algorithm requires primes p that **split completely** in K — primes for which f
has d distinct roots mod p. The same `select_qc_primes` machinery from G.E is reused:

```rust
fn select_couveignes_primes(f: &IntPoly, b_alg: u64, num_primes: usize) -> Vec<u64> {
    select_qc_primes(f, b_alg, num_primes)
}
```

**Principle-4 annotation (the prime-count scale knob).** The constant
`DEFAULT_COUVEIGNES_PRIMES = 10` suffices at toy scale. The coefficient bit-length of $\beta$
is bounded by $|S| \cdot \max(\log|a_i|, \log|b_i|) \cdot d$, which at toy scale (|S| ~ 2,
coefficients ~ 4 bits, d ~ 2) is tiny — a handful of primes is ample margin. At NFS scale
(|S| ~ 10⁵, coefficients ~ 100 bits, d ~ 5), the prime count grows as
$O(\text{coefficient\_bits} / 64)$. The algorithm is identical at all scales; only the count
changes.

### Step 3: Per-prime square root and Lagrange interpolation

For each split prime p with roots $r_1, \ldots, r_d$ of f mod p:

**3a. Reduce γ to 𝔽_p.** For each root $r_j$, call `gamma.reduce_mod_ideal(&p_big, &r_j_big)`
to get $\gamma_j \in \mathbb{F}_p$. The `reduce_mod_ideal` method (added to
`NumberFieldElement` at G.F.1) evaluates the element's polynomial at $\alpha \equiv r_j
\pmod{p}$, clearing denominators mod p.

**3b. Compute the per-root square root.** Call `tonelli_shanks(gamma_j, p)` to get
$\beta_j \in \mathbb{F}_p$ with $\beta_j^2 \equiv \gamma_j \pmod{p}$. The standalone
`tonelli_shanks` function implements the full Tonelli–Shanks algorithm:

```rust
pub(crate) fn tonelli_shanks(n: u64, p: u64) -> Option<u64> {
    if n == 0 { return Some(0); }
    // Euler's criterion: n^((p-1)/2) mod p must be 1.
    let legendre = pow_mod(n, (p - 1) / 2, p);
    if legendre != 1 { return None; }
    // Special case: p ≡ 3 (mod 4) → sqrt = n^((p+1)/4) mod p.
    if p % 4 == 3 { return Some(pow_mod(n, (p + 1) / 4, p)); }
    // General Tonelli–Shanks: factor p-1 = Q * 2^S with Q odd, ...
    // ...
}
```

If `tonelli_shanks` returns `None`, $\gamma$ is not a quadratic residue mod p — an upstream
kernel bug (the QC columns should guarantee $\gamma$ is a square in K, hence a QR mod every
split prime).

**3c. Lagrange interpolation mod p.** Given the d pairs $(r_j, \beta_j)$, reconstruct the
unique polynomial $\beta \bmod p$ of degree $< d$ such that $\beta(r_j) = \beta_j$ for all j.
This is the unique polynomial in $\mathbb{F}_p[x]/(f \bmod p)$ that squares to $\gamma \bmod p$:

$$
\beta(x) = \sum_{j=0}^{d-1} \beta_j \cdot L_j(x), \quad L_j(x) = \prod_{k \neq j} \frac{x - r_k}{r_j - r_k}
$$

The implementation in `lagrange_interp_mod_p` builds the result polynomial coefficient by
coefficient, accumulating the contribution of each basis polynomial $L_j$.

**Sign consistency.** The per-root square roots have a $\pm 1$ ambiguity: Tonelli–Shanks
returns one of $\{\beta_j, p - \beta_j\}$. The implementation resolves this incrementally:
after fixing $\beta_0$, for each subsequent root $r_j$, it evaluates the partial Lagrange
polynomial at $r_j$ and chooses the square root that matches. After interpolation, the result
is verified unconditionally:

```rust
// Verify unconditionally: the result polynomial must square to γ at each root.
for (j, (&r, &gamma_j)) in roots.iter().zip(gamma_vals.iter()).enumerate() {
    let beta_at_r = eval_poly_u64(&result, r, p);
    let beta_sq = mul_mod(beta_at_r, beta_at_r, p);
    assert_eq!(
        beta_sq, gamma_j,
        "per_prime_beta: sign inconsistency — β(r_{j})² ≠ γ(r_{j}) (mod p={p}); \
         upstream kernel bug"
    );
}
```

This unconditional check (not a debug-only assert) is the G.F.3 review juncture's finding:
a wrong sign at this step is a silent failure that would propagate into the CRT lift and
produce a wrong $\beta$.

### Step 4: CRT lift via Garner's algorithm

After processing all primes, we have for each coefficient index k (0 ≤ k < d) a set of
residues $\{c_k \bmod p_1, c_k \bmod p_2, \ldots\}$. Garner's algorithm lifts these to a
single $c_k \in \mathbb{Z}$:

$$
x = a_1 + p_1(a_2 + p_2(a_3 + \cdots))
$$

where the mixed-radix coefficients $a_i$ are computed iteratively. The implementation in
`garner_crt` returns the unique $x$ with $0 \le x < \prod p_i$; a centering step then maps
$x$ to the symmetric range $(-\prod p_i / 2, \prod p_i / 2]$:

```rust
let prod: BigInt = primes.iter().map(|&p| BigInt::from(p)).product();
let half_prod = &prod / BigInt::from(2u64);
let ck = if ck_raw > half_prod { ck_raw - &prod } else { ck_raw };
```

The centered representation gives $\beta$'s coefficients as signed integers, which is the
natural form for the norm computation.

### Step 5: Sign resolution via the real embedding

After the CRT lift, $\beta$ is determined up to a global sign: both $\beta$ and $-\beta$
satisfy $\beta^2 = \gamma$. The wrong sign gives $Y = \mathrm{Norm}(-\beta) = (-1)^d \cdot
\mathrm{Norm}(\beta)$, which for odd d flips the sign of Y. Even if $|Y|$ is correct, the
wrong sign can yield a trivial GCD at the assembly step.

**The resolution.** Use the **real embedding** of K. Since f is the GNFS algebraic polynomial,
it has at least one real root $\theta$ (the polynomial arises from base-m expansion of N, which
is positive, so f has a real root near $m^{1/d}$). Evaluate $\beta$ at $\theta$ to get a real
number $\beta(\theta)$. Choose the sign of $\beta$ such that $\beta(\theta) > 0$:

```rust
let m_f64 = poly.m.to_f64().unwrap_or(1.0);
let x0 = m_f64.powf(1.0 / d as f64);
let sign_positive = match find_real_root(&f, x0, 1e-10, 200) {
    Some(theta) => {
        let val = eval_at_real(&coeffs, theta);
        val >= 0.0  // β(θ) > 0 means we have the positive square root
    }
    None => true,  // fallback: keep the CRT result as-is
};
if !sign_positive {
    for c in &mut coeffs { *c = -c.clone(); }
}
```

The real root $\theta$ is found by Newton's method starting from $m^{1/d}$ (a good initial
guess for base-m polynomials). At toy scale, f64 precision suffices.

**Why this works.** The real embedding $\alpha \mapsto \theta$ is a ring homomorphism $K \to
\mathbb{R}$. For $\gamma = \beta^2$, we have $\gamma(\theta) = \beta(\theta)^2 \ge 0$. By
choosing $\beta(\theta) > 0$, we select the "positive" square root under the real embedding,
consistent with the rational side's sign convention (X > 0 from `isqrt`).

**Principle-4 annotation (embedding-sign resolution).** The embedding-sign resolution is
**not** a scale artifact — it is a correctness obligation present even at toy scale. The
G.F.3 review juncture identified this as the silent-failure locus: a wrong-sign $\beta$
returns a plausible Y that produces a trivial GCD at assembly, not a red test. The G.F.4
retry loop (§46) is the safety net, but the sign resolution is the primary defence.

### Step 6: Compute Y = |Norm(β)| mod N

Once $\beta$'s coefficients are determined and sign-resolved, compute the norm:

```rust
let beta_poly = RatPoly::from_coeffs(
    coeffs.iter().map(|c| BigRational::from(c.clone())).collect(),
);
let beta = NumberFieldElement { field: &nf, poly: beta_poly };
let norm_rat = beta.norm();
let norm_abs = norm_rat.numer().clone().abs();
norm_abs % n
```

The norm of an algebraic integer is a rational integer; the denominator of `norm_rat` should
be 1. The absolute value is taken because the norm may be negative (for odd-degree fields);
the assembly step uses $|Y|$ mod N.

---

## 45. The CRT Congruences: A Summary

The Couveignes algorithm is most clearly understood through the CRT congruences it satisfies.
For each split prime p with roots $r_1, \ldots, r_d$ of f mod p:

$$
\beta(r_j) \equiv \sqrt{\gamma(r_j)} \pmod{p}, \quad j = 1, \ldots, d
$$

These d congruences (one per root) determine $\beta \bmod p$ as a polynomial of degree $< d$
via Lagrange interpolation. Across all CRT primes $p_1, \ldots, p_k$:

$$
\beta \equiv \beta_1 \pmod{p_1}, \quad \beta \equiv \beta_2 \pmod{p_2}, \quad \ldots, \quad \beta \equiv \beta_k \pmod{p_k}
$$

where each $\beta_i$ is a polynomial in $\mathbb{F}_{p_i}[x]/(f \bmod p_i)$. The CRT lift
recovers $\beta \in \mathbb{Z}[\alpha]$ coefficient by coefficient: for each coefficient index
k, the residues $\{c_k \bmod p_i\}$ are lifted to $c_k \in \mathbb{Z}$ via Garner's algorithm.

The key identity that makes this work is:

$$
\beta^2 = \prod_{i \in S}(a_i - b_i \cdot \alpha)
$$

which holds in K by the QC-column guarantee. The CRT lift recovers $\beta$ (up to global sign)
from its residues mod the split primes.

---

## 46. The Final GCD Assembly

### The congruence of squares identity

The assembly step is implemented in `gnfs/src/sqrt/assembly.rs`. Given X (from `rational_sqrt`)
and Y (from `algebraic_sqrt`), the GNFS congruence-of-squares identity guarantees:

$$
X^2 \equiv Y^2 \pmod{N}
$$

This means $N \mid (X - Y)(X + Y)$. If $X \not\equiv \pm Y \pmod{N}$, then neither factor
$(X - Y)$ nor $(X + Y)$ is divisible by N, so $\gcd(X - Y, N)$ is a non-trivial factor of N.

### The `factor_from_congruence` function

```rust
pub fn factor_from_congruence(x: &BigInt, y: &BigInt, n: &BigInt) -> Option<BigInt> {
    // Primary check: gcd(X − Y, N).
    let diff = x - y;
    let g = gcd(&diff.abs(), n);
    if is_nontrivial(&g, n) { return Some(g); }

    // Fallback: gcd(X + Y, N). Catches the X ≡ −Y (mod N) case.
    let sum = x + y;
    let g2 = gcd(&sum.abs(), n);
    if is_nontrivial(&g2, n) { return Some(g2); }

    None
}
```

The `gcd` function is from `shared_bigint` (exported at G.F.1). The `is_nontrivial` predicate
checks $1 < g < N$: $g = 1$ means the GCD is trivial (the congruence gave no information);
$g = N$ means $N \mid (X - Y)$, i.e., $X \equiv Y \pmod{N}$ (also trivial).

### Why the trivial GCD is expected

The trivial-GCD outcome ($X \equiv \pm Y \pmod{N}$) is *expected* for some kernel vectors —
it is not a bug. The probability that a random congruence of squares yields a trivial GCD is
roughly 50% for a semiprime $N = p \cdot q$ (because $X \equiv \pm Y \pmod{p}$ and
$X \equiv \pm Y \pmod{q}$ are independent events, each with probability 1/2). The driver
loops over multiple kernel vectors and returns the first non-trivial factor.

### The `factor` driver

The top-level `factor` function chains the whole pipeline:

```rust
pub fn factor(
    poly: &PolyPair,
    matrix: &SparseMatrix,
    relations: &[Relation],
    kernel_vectors: &[KernelVector],
) -> Option<BigInt> {
    for kv in kernel_vectors {
        let x = rational_sqrt(kv, matrix, relations, poly);
        let y = algebraic_sqrt(kv, matrix, relations, poly);
        if let Some(f) = factor_from_congruence(&x, &y, &poly.n) {
            return Some(f);
        }
    }
    None
}
```

For each kernel vector: compute X, compute Y, attempt the GCD. Return the first non-trivial
factor, or `None` if all kernel vectors yield trivial GCDs.

---

## 47. The Whole Factoring Arc

This chapter closes the GNFS pipeline. The five stages, end to end:

1. **Polynomial selection (G.B, §1–§11).** Choose a polynomial pair $(f, g)$ with $f(m) \equiv
   0 \pmod{N}$ and good Murphy-E score. Output: `PolyPair`.

2. **Sieving (G.C, §12–§21).** For each $(a, b)$ in the sieve region, check whether both the
   rational norm $a - b \cdot m$ and the algebraic norm $b^d \cdot f(a/b)$ are smooth over
   their respective factor bases. Output: `Vec<Relation>`.

3. **Filtering (G.D, §22–§30).** Remove duplicate and linearly dependent relations; merge
   weight-2 and weight-3 columns; build the sparse GF(2) matrix with provenance map. Output:
   `SparseMatrix`.

4. **Linear algebra (G.E, §31–§41).** Find the left nullspace of the filtered matrix over
   GF(2) via block Lanczos or block Wiedemann. Fill in the QC columns to guarantee the
   algebraic square root exists in K. Output: `Vec<KernelVector>`.

5. **Square root and assembly (G.F, §42–§47, this chapter).** For each kernel vector: recover
   the original $(a, b)$ pairs via `expand_provenance`; compute the rational square root X and
   the algebraic square root Y via Couveignes' CRT algorithm; form $\gcd(X - Y, N)$. Output:
   a non-trivial factor of N.

The pipeline is visible end-to-end in the `factor` function in `gnfs/src/sqrt/assembly.rs`,
which chains `rational_sqrt` and `algebraic_sqrt` and calls `factor_from_congruence`. The
end-to-end KAT in `gnfs/tests/factor_end_to_end_kat.rs` verifies this for N = 35 = 5 × 7:

- **Setup.** $f(x) = x^2 - 2$ (K = ℚ(√2)), $m = 6$, $N = 35$.
- **Relations.** $(a=4, b=0)$ and $(a=9, b=0)$: rational factors 4 and 9.
- **Rational product.** $4 \times 9 = 36 = 6^2$. $X = 6$. $X \bmod 35 = 6$.
- **Algebraic product.** $\gamma = 4 \times 9 = 36 = 6^2$ in K (rational). $\beta = 6$.
  $\mathrm{Norm}(6) = 6^2 = 36$ (degree-2 field). $Y = 36 \bmod 35 = 1$.
- **Assembly.** $\gcd(X - Y, N) = \gcd(5, 35) = 5$. Non-trivial factor: **5**. ✓

---

## 48. Code-Tour Passages

### `gnfs/src/sqrt/mod.rs` — module entry point

The `sqrt` module re-exports the three public entry points:

```rust
pub use algebraic::algebraic_sqrt;
pub use assembly::{factor, factor_from_congruence};
pub use rational::rational_sqrt;
```

The module docstring explains the C-LinAlg seam: `KernelVector::expand_provenance` was
over-specified at G.E.1 precisely for this consumer.

### `gnfs/src/sqrt/rational.rs` — rational square root

Key functions:
- `rational_sqrt(kv, matrix, relations, poly) -> BigInt` — the public entry point. Calls
  `expand_provenance`, forms the product, calls `isqrt`, reduces mod N.
- The sign check (`if product.is_negative() { panic!(...) }`) is the upstream-bug detector.
- The `isqrt` call is from `shared_bigint::isqrt` (added at G.F.1).

### `gnfs/src/sqrt/algebraic.rs` — Couveignes' algorithm

Key functions:
- `algebraic_sqrt(kv, matrix, relations, poly) -> BigInt` — the public entry point.
- `form_gamma(nf, relations, s) -> NumberFieldElement` — forms $\gamma = \prod(a_i - b_i\alpha)$.
- `select_couveignes_primes(f, b_alg, num_primes) -> Vec<u64>` — delegates to `select_qc_primes`.
- `per_prime_beta(gamma, f, p, d) -> Vec<u64>` — per-prime square root + Lagrange interpolation.
  Contains the unconditional sign-consistency check (the G.F.3 review juncture's finding).
- `tonelli_shanks(n, p) -> Option<u64>` — standalone Tonelli–Shanks (deviation from the
  C-AlgSqrt design's `FpNaive::<4>::sqrt` path; same algorithm, simpler interface at toy scale).
- `lagrange_interp_mod_p(points, p, d) -> Vec<u64>` — Lagrange interpolation mod p.
- `garner_crt(residues) -> BigInt` — CRT lift via Garner's algorithm.
- `find_real_root(f, x0, tol, max_iter) -> Option<f64>` — Newton's method for sign resolution.
- `eval_at_real(coeffs, theta) -> f64` — evaluate $\beta$ at the real root $\theta$.

### `gnfs/src/sqrt/assembly.rs` — final assembly

Key functions:
- `factor_from_congruence(x, y, n) -> Option<BigInt>` — computes $\gcd(X - Y, N)$ and
  $\gcd(X + Y, N)$; returns the first non-trivial factor.
- `factor(poly, matrix, relations, kernel_vectors) -> Option<BigInt>` — the top-level driver.
  Iterates over kernel vectors, calls both square roots, calls `factor_from_congruence`.
- `is_nontrivial(g, n) -> bool` — checks $1 < g < N$.

---

## 49. Principle-4 Annotations Summary (G.F)

The following table collects all principle-4 annotations from the G.F module — phenomena that
are scale-dependent and either under-exposed or over-exposed at toy scale.

| Phenomenon | Toy-scale behaviour | NFS-scale behaviour | Scale knob |
|------------|--------------------|--------------------|------------|
| **CRT prime count** | 10 primes suffice; coefficient bit-length is tiny | $O(\text{coeff\_bits} / 64)$ primes needed | `DEFAULT_COUVEIGNES_PRIMES` |
| **Embedding-sign resolution** | Correctness obligation at all scales; not scale-dependent | Same algorithm; same obligation | — (correctness, not scale) |
| **Rational product size** | Product fits in a few hundred bits | Product can be thousands of bits; `isqrt` must handle large integers | — (BigInt arithmetic) |
| **Trivial-GCD probability** | ~50% per kernel vector for a semiprime | Same probability; more kernel vectors available from larger nullspace | Number of kernel vectors |

The embedding-sign resolution is the one phenomenon that is **not** a scale artifact: it is a
correctness obligation present even at toy scale, and the primary target of the G.F.3 review
juncture. The CRT prime count is the primary scale knob: the algorithm is identical at all
scales; only the count changes.

---

## 50. KAT Summary (G.F)

The following table lists the key known-answer tests across the sqrt module, with the
mathematical fact each one verifies.

| Test file | Test | Mathematical fact verified |
|-----------|------|---------------------------|
| `sqrt_rational_kat.rs` | `kat_a_x_squared_congruent_to_product_mod_n` | X² ≡ ∏(a_i − b_i·m) (mod N) for a hand-built relation set |
| `sqrt_rational_kat.rs` | `kat_b_correct_index_set_from_expand_provenance` | The index set S matches a hand-traced `expand_provenance` |
| `sqrt_rational_kat.rs` | `kat_c_non_square_product_panics` | A non-square product triggers the upstream-bug panic |
| `sqrt_rational_kat.rs` | `single_relation_perfect_square` | Single relation with a perfect-square rational norm |
| `sqrt_rational_kat.rs` | `m_dependent_factors_correct` | Product uses the correct m from `PolyPair` |
| `sqrt_algebraic_kat.rs` | `kat_a_known_square_norm_matches` | For a known-square γ, Couveignes recovers β with Norm(β) matching the hand-computed Y |
| `sqrt_algebraic_kat.rs` | `kat_a_y_squared_congruent_to_norm_gamma` | Y² ≡ Norm(γ) (mod N) for a hand-built K and γ |
| `sqrt_algebraic_kat.rs` | `kat_b_x_squared_congruent_to_y_squared_mod_n` | X² ≡ Y² (mod N) for a toy N (congruence holds before the full driver) |
| `sqrt_algebraic_kat.rs` | `kat_c_determinism` | Couveignes is deterministic for a fixed γ and prime set |
| `sqrt_algebraic_kat.rs` | `bonus_degree2_two_relations` | Two-relation degree-2 case: β recovered correctly |
| `factor_end_to_end_kat.rs` | `kat_a_factor_from_congruence_known_values` | gcd(6 − 1, 35) = 5; non-trivial factor of N = 35 |
| `factor_end_to_end_kat.rs` | `kat_a_factor_driver_recovers_known_factor` | Full driver recovers 5 or 7 from N = 35 |
| `factor_end_to_end_kat.rs` | `kat_b_trivial_gcd_x_equals_y` | X = Y gives trivial GCD (None) |
| `factor_end_to_end_kat.rs` | `kat_b_retry_loop_skips_trivial_and_finds_factor` | Retry loop skips KV1 (trivial) and returns 5 from KV2 |
| `factor_end_to_end_kat.rs` | `kat_b_all_trivial_returns_none` | All trivial kernel vectors → `factor` returns None |

---

## The Checked Square-Root Context

The square-root stage is correct only under algebraic preconditions that the legacy path never
checks. `gnfs::sqrt::context` makes those preconditions a type: a value of `SqrtContext` is a proof
that they held when it was built, so later code can rely on them without re-checking.

### What a context certifies

Given `f` of degree `d ≥ 1`, an integer `m`, and `N`, `SqrtContext::new` accepts exactly when

1. `N > 1`;
2. `f(m) ≡ 0 (mod N)`;
3. the leading coefficient `c` is nonzero and `gcd(c, N) = 1`;
4. the monic polynomial `F(t) = c^(d−1) f(t/c)` is squarefree, i.e. `disc(F) ≠ 0`.

`F` has coefficient `f_k · c^(d−1−k)` below the leading term and is monic, so `θ = c·α` is an
algebraic integer. Because `F(c·m) = c^(d−1) f(m) ≡ 0 (mod N)`, sending `θ ↦ c·m` is a ring
homomorphism from `Z[θ]` to `Z/N`, which is how an algebraic root is carried back to a residue
modulo `N`. The context stores `f` as given (so `c` can be negative) and `m` unreduced, because the
rational product `Π(a − b·m)` is taken over the integers.

Squarefree `F` does not make `A = Q[t]/(F)` a field. For `x² − 1` at `m = 4`, `N = 15`, `A ≅ Q × Q`.
The context accepts that case and makes no irreducibility claim; a demonstration that needs a field
must certify irreducibility separately.

### Failure is classified, not collapsed

The checks run in a fixed order and stop at the first failure: `N > 1` is first, so nothing is
reduced modulo zero; malformed input (constant `f`, wrong root, repeated factor) precedes the domain
checks, so a bad input is never reported as luck. A proper `gcd(c, N)` is an *incidental factor of
`N`*, reported with its provenance and never treated as a reconstruction result. `N | c` is an
unsupported domain: the unscaled map does not exist. `SqrtError` carries invalid input and violated
internal invariants, which are never retried as nonsquares; `AlgebraicSqrtOutcome` carries
mathematical and resource results. Only `NotSquare` certifies that no root exists; resource
exhaustion is never evidence.

### Relation factors must be units

A factor `a − b·α` enters the product `γ = Π(a − b·α)`. Its scaled form `c·a − b·θ` has norm

    Π_j (c·a − b·θ_j) = b^d F(c·a/b) = c^(d−1) · f_hom(a, b),   f_hom(a, b) = Σ f_i a^i b^(d−i),

so the factor is a unit of `A` exactly when `f_hom(a, b) ≠ 0`. A zero `f_hom` means the factor
vanishes at some root of `F`: it is a zero divisor (or the zero element, always so when `d = 1`),
and no prime would ever make that component invertible, so searching for one would not terminate.
`RelationFactors::new` therefore checks every pair against the context before any prime search. In
`Q[t]/(t² − 1)` the pair `(1, 1)` gives `f_hom = 1 − 1 = 0`: the factor `1 − θ` kills the root
`θ = 1`. For `f = x − 5` the pair `(5, 1)` gives the zero element of `A = Q`.

### Why `VerifiedRoot` cannot be forged

`VerifiedRoot` has private fields and one private constructor, `seal`. Rust private items are
visible only in their module and its descendants, so every module outside the `sqrt::context`
subtree, and every other crate, is excluded by the compiler from creating one. The module that
performs the exact check `β² = γ` over `Q` and `y² ≡ Π(a − b·m) (mod N)` is meant to be a child of
`context`, which is why it can call `seal` while sibling modules such as the reconstruction code
cannot. A sibling module would need `pub(super)` visibility on the constructor, which would also be
visible to every other module under `sqrt`. Two `compile_fail` doctests on `VerifiedRoot` pin this
property. `NonSquareCertificate` is weaker: its constructor is `pub(super)`, so only the compiler
protects it from outside `gnfs::sqrt`.

### Resource options

`PrimeSearchLimits`, `PrecisionLimits`, `DirectLimits`, `CrtLimits`, `RecombinationLimits`, and
`SqrtLimits` are plain `Copy` structs with no `Default`: a caller must state its budgets. A prime
search examines primes `p ≥ max(first_candidate, 3)` in ascending order and every prime examined
counts, including primes skipped later. Reaching a limit yields a `ResourceExhaustion` naming the
resource, the amount consumed, and the limit.

### KAT summary

| Test file | Test | Mathematical fact verified |
|-----------|------|---------------------------|
| `sqrt_context_kat.rs` | `kat_context_accepts_monic_and_non_monic` | `x²−2` (`N=527`) and `2x²−1` (`N=7`) give `F=t²−2`, `disc 8`, `θ ↦ c·m mod N` |
| `sqrt_context_kat.rs` | `kat_context_rejects_modulus_le_one` | `N ∈ {1, 0, −5}` is rejected first, without division |
| `sqrt_context_kat.rs` | `kat_context_rejects_constant_and_zero_polynomial` | Degree `< 1` is invalid |
| `sqrt_context_kat.rs` | `kat_context_rejects_root_mismatch` | `x²−2` at `m=6`, `N=35` has `f(m)=34` |
| `sqrt_context_kat.rs` | `kat_from_poly_pair_rejects_degree_field_and_rational_poly` | A literal pair is checked without panicking, including `N=0` |
| `sqrt_context_kat.rs` | `kat_context_rejects_repeated_factor` | `(x−3)²` and `(2x−1)²` have `disc 0` |
| `sqrt_context_kat.rs` | `kat_context_incidental_leading_coefficient_factor` | A proper `gcd(c, N)` is an incidental factor with provenance |
| `sqrt_context_kat.rs` | `kat_context_invalid_precedes_incidental` | Malformed input is never reported as an incidental factor |
| `sqrt_context_kat.rs` | `kat_context_leading_coefficient_not_unit_unsupported` | `N \| c` is an unsupported domain |
| `sqrt_context_kat.rs` | `kat_context_reducible_squarefree_accepted` | `x²−1` is accepted; `A` is not a field |
| `sqrt_context_kat.rs` | `kat_relation_factors_boundary` | Empty set, zero factors, zero divisors, `d=1`, and first-failing index |
| `sqrt_context_kat.rs` | `kat_failure_classes_distinct` | Each condition maps to its own `FailureClass` |

---

## The Normalized Target and Exact Verification

A reconstruction algorithm proposes an integer polynomial `z`; it never produces a trusted root.
`gnfs::sqrt::context::normalize` owns the two things every algorithm shares: the *target* that `z`
must square to, with a certified bound on its coefficients, and the *verifier* that is the only way
to obtain a `VerifiedRoot`. It is a child of `context` so that it can call the private `seal`.

### From relations to an integral target

With `A = Q[t]/(F)`, `θ = c·α`, `k` relations and `h = ⌈k/2⌉`, the algebraic product is
`γ = Π(a − b·α) = Π(c·a − b·θ) / c^k`. It need not be integral, but

    η = c^(2h−k) · Π(c·a − b·θ) = c^(2h) · γ      (in Z[θ])
    H = F'(θ)² · η                                (in Z[θ])

are. If `β² = γ` then `w = c^h β` satisfies `w² = η` and is integral at every conjugate, and
`z = F'(θ)·w` satisfies `z² = H` exactly modulo `F`. Multiplying by `F'` is what puts `z` in the power
basis `Z[θ]` even when `Z[θ]` is not the maximal order: `z_i = Tr(w·h_i(θ))` with
`h_i(t) = Σ_{ℓ>i} F_ℓ t^(ℓ−1−i)`, and each `w·h_i(θ)` is an algebraic integer, so its trace is an
integer. The unscaled root is recovered over `Q` as `β = z / (c^h F'(θ))`.

The root need not lie in `Z[θ]`. For `F = t² − 8` and `γ = 2`, the root is `β = θ/2`, with `β² =
θ²/4 = 2`, but the integral root is `z = F'β = θ² = 8`; a search restricted to `w ∈ Z[θ]` would
wrongly reject this square.

### Certified coefficient bound

All arithmetic is exact integer arithmetic:

    R = 1 + max_{i<d} |F_i|           (Cauchy: every root of F has |θ_j| ≤ R)
    E = Σ_{i<d} |η_i| R^i             (so |η(θ_j)| ≤ E and |w_j| ≤ √E)
    U = ⌈√E⌉                          (ceiling of an exact integer square root)
    B_i = d · U · Σ_{ℓ>i} |F_ℓ| R^(ℓ−1−i),    B = max_i B_i

`|z_i| ≤ B_i` for every possible root, so a reconstruction modulus `M` with `M > 2B` determines `z`
from its centered residues; the least such integer is `2B + 1`. `⌈√E⌉` is computed from the floor
root, adding one unless it squares to `E`. The square root of a nonsquare is therefore neither
truncated nor approximated, and floating point is never used. The inner sum is evaluated downward as
`S_(i−1) = |F_i| + R·S_i`. For `F = t³ − 2t + 3` with `η = 3 − t`, `R = 4`, `E = 7`, `U = 3`
(the floor would give 2), and `B = (162, 36, 9)`. The bound is deliberately conservative, trading
extra precision for a proof that needs no numerical embeddings.

`NormalizedTarget::new` checks `bits(2B + 1)` against `PrecisionLimits` before returning, so a
caller learns of an unaffordable precision before it searches for a single prime. It also verifies
the inverse of `F'` modulo `F` (computed over `Q` by the extended Euclidean algorithm) rather than
trusting a generic field inverse; `F` squarefree makes it exist.

### What `verify_candidate` checks, and what it concludes

The checks run in a fixed order:

1. `deg z < d`; 2. `|z_i| ≤ B_i`; 3. `z² ≡ H (mod F)` over `Z`. A failure here is a
   `CandidateRejection`. It says that this `z` is not a root. It never says that no root exists.
4. `β = z/(c^h F')` and `β² = γ` over `Q`. This follows from step 3 for an untampered target, so a
   failure is an `Err(Invariant)`, a defect, never a rejection.
5. Every reduced denominator of `β` is a unit modulo `N`. A denominator sharing a proper factor with
   `N` yields an `Incidental` factor (source: map denominator), preferred over an unsupported one
   found earlier. A denominator divisible by `N` yields `Unsupported(MapDenominatorNotUnit)`. A root
   exists in both cases, only the coefficientwise map fails, and neither is a nonsquare
   certificate. For `F = t² − 8` and `β = θ/2` the only denominator is 2: modulo `N = 4` that gives the
   factor `2·2`, modulo `N = 2` it is unsupported, and modulo `N = 7` it is a unit.
6. `y = β(c·m) mod N`. Writing `β = Σ b_i θ^i`, each `b_i` is reduced modulo `N` and the polynomial is
   evaluated at `c·m`; this is a ring homomorphism on the localization because `F(c·m) = c^(d−1) f(m)
   ≡ 0 (mod N)`.
7. `y² ≡ P = Π(a − b·m) (mod N)`, again an invariant for the same reason.

### The map is not the norm

The image `y` squares to the rational product `P`, not to `Norm(γ)`. For `f = x² − 2`, `m = 3`,
`N = 7` and the factor `9`: `β = 3`, `y = β(3) = 3` and `3² = 9 ≡ 2 ≡ P`. The norm `Norm(β) = 9 ≡ 2`
squares to `4 ≢ 2`. A candidate built to satisfy the norm equation, such as `β = 2` (`z = 4θ`), has
`z² = 32 ≠ 72 = H` and is rejected.

### Non-monic fixtures

For `f = 2x² − 1`, `m = 2`, `N = 7` and the factor `(3, −4)`: `θ = 2α`, `F = t² − 2`,
`γ = 3 + 2θ = (1 + θ)²`, `η = 12 + 8θ`, and `z = 8 + 4θ` with `z² = 96 + 64θ = H`. The image is
`y = β(4) = 5` with `5² ≡ 4 ≡ 11 = P`. A negative leading coefficient is handled the same way, with
`c^h` carrying the sign.

### KAT summary

| Test file | Test | Mathematical fact verified |
|-----------|------|---------------------------|
| `sqrt_normalize_kat.rs` | `kat_bound_formula_hand_computed` | `R`, `E`, `U = ⌈√E⌉`, `B_i`, `B` for four hand-computed targets, including nonsquare `E` |
| `sqrt_normalize_kat.rs` | `kat_target_identities_hold` | `γ`, `η = c^(2h)γ`, `H = F'²η`, and `F'·F'⁻¹ ≡ 1` for a non-monic input |
| `sqrt_normalize_kat.rs` | `kat_odd_cardinality_non_monic_accepts` | `2x²−1` accepts `z = 8 + 4θ`, `y = 5`, both signs |
| `sqrt_normalize_kat.rs` | `kat_even_cardinality_and_assembly_seam` | Even `k`, and `gcd(7 ∓ 24, 527)` gives 17 and 31 |
| `sqrt_normalize_kat.rs` | `kat_negative_leading_coefficient_accepts` | `c = −2` carries the sign through `c^h` and the map |
| `sqrt_normalize_kat.rs` | `kat_degree_one_target` | `d = 1`, `F' = 1`, a root exactly on its bound |
| `sqrt_normalize_kat.rs` | `kat_outside_order_root_accepts` | `β = θ/2` is outside `Z[θ]` but `z = 8` verifies |
| `sqrt_normalize_kat.rs` | `kat_zero_divisor_factors_are_rejected_before_a_target_exists` | A zero divisor is stopped at the relation boundary |
| `sqrt_normalize_kat.rs` | `kat_wrong_norm_map_counterexample` | `y² ≡ P`, while `Norm(β)² ≢ P`; norm-equation candidates are rejected |
| `sqrt_normalize_kat.rs` | `kat_noninvertible_denominators_have_typed_outcomes` | Proper gcd is incidental; `N \| denominator` is unsupported |
| `sqrt_normalize_kat.rs` | `kat_candidate_rejections_are_classified_in_order` | Degree, bound, and square checks in order; inclusive bounds |
| `sqrt_normalize_kat.rs` | `kat_modulus_precheck_precedes_any_search` | `bits(2B + 1)` against the precision budget |
| `normalize.rs` unit tests | `corrupted_*`, `non_invertible_derivative_*`, `ceil_sqrt_*` | Invariant failures are errors, not rejections; exact ceiling root |

---

## Direct p-Adic Reconstruction

`gnfs::sqrt::padic::direct_root` finds the algebraic square root at a single inert prime. It takes a
`NormalizedTarget` and `DirectLimits` and returns an `AlgebraicSqrtOutcome`. It proposes an integer
polynomial `z` and never trusts it: success exists only through `verify_candidate`.

### Why an inert prime makes the root unique up to sign

Take an odd prime `p` that does not divide `c·N`, for which `F` is irreducible modulo `p`, and with
`H ≢ 0 (mod p)`. Then `(Z/p^e)[t]/(F)` is a local ring with residue field `F_(p^d)`, and `H` is a
unit. A unit square root is unique up to sign, because a second root `z'` would give
`(z/z')² = 1` with `z/z' ≡ ±1` modulo the maximal ideal, and Hensel's lemma lifts only those two
signs. So the only choice the algorithm faces is a global `±`, and both signs pass the exact
verifier identically. There is no real-sign convention and no fixed number of primes, which is why
a totally imaginary field such as `Q(i)` needs no special handling.

The prime is found by ascending search from `max(first_candidate, 3)`. Each prime is classified: it
is skipped when it is even, divides `c` or `N`, leaves `F` reducible (the Rabin test also rejects a
repeated factor), or makes `H ≡ 0`. Every prime examined counts against `max_primes_examined`,
including skipped ones. `H ≡ 0` is rejected before any square root is taken, because zero has no
inverse square root. Running out of budget, or of machine primes below `2^64`, is a resource
outcome: it does not prove that an inert prime is absent. For `F = t⁴ − 10t² + 1`, whose Galois group
is `V4`, no inert prime exists at all, and the search ends only by budget.

For `d ≥ 2` a prime dividing `N` or `c` already leaves `F` with a root, so the divisibility guard
changes the outcome only for `d = 1`, where `F = t − c·m` is irreducible at every prime.

### Seed and Newton step

In the residue field a Tonelli–Shanks square root `s` of `H` gives the seed `u₀ = s⁻¹`, with
`H·u₀² ≡ 1 (mod p)`. The lift is

    u ← u·(3 − H·u²)/2

and if `H·u² = 1 + ε` the new value satisfies `H·u'² = 1 − (3/4)ε² + ε³/4`. The error is squared, so
each step doubles the precision. Division by 2 multiplies by `(M + 1)/2`, which is valid because
`M = p^e` is odd. The precision runs through the ceiling-halving chain of the final exponent
(`5 → 3 → 2 → 1`), so a step never more than doubles. Each step checks `H·u² ≡ 1` at the precision it
starts from, and the last result is checked at the final precision. A failed check is an
`Err(Invariant)`, never a nonsquare.

The ring arithmetic is private to this module. An element is `d` coefficients in `[0, p^e)`;
multiplication is the schoolbook product followed by elimination of `t^(2d−2), …, t^d` through
`t^d = −(F_0 + … + F_(d−1)t^(d−1))`. Moduli are `BigInt`, never limited to a machine word. The
existing scalar `Zp` type is not this ring: it is not a quotient by `F`.

### Centering and the bound

`z = H·u` satisfies `z² = H²u² = H` modulo `(F, p^e)`. Choose the least `e` with `p^e > 2B`, where
`B` is the certified coefficient bound of the target. Every true root has `|z_i| ≤ B ≤ (M − 1)/2`,
so its centered residues are its coefficients exactly. Since `p^(e−1) ≤ 2B`, we have
`p^e < p·(2B + 1)`, so the modulus has at most `bits(2B + 1) + bits(p)` bits. The precision budget
is enforced on `2B + 1` when the target is built, so `direct_root` takes no precision argument. The boundary is sharp in both directions. For `F = t − 5` and `H = 169 = 13²`
(`B = 13`) the least modulus `3³ = 27 = 2B + 1` works, with the root `13` exactly on the bound. For
`H = 196 = 14²` (`B = 14`) `3³ = 27 < 2B + 1`: `14 ≡ −13 (mod 27)` would reconstruct `−13`, whose
square is `169`, so `3⁴ = 81` is required.

### What is certified

- **Local obstruction.** `H` is a nonresidue in `F_(p^d)`. By the Euler criterion this holds
  exactly when `H^((q−1)/2) ≠ 1`. Over `F = t² − 2` at `p = 3`, `8(1 − θ) ≡ 2 + θ` has norm
  `2 ≡ −1`: a nonresidue, with no lifting at all. The certificate is `InertLocalNonresidue`.
- **Full-bound failure.** `H` is a square modulo `p^e` but the centered `z` is rejected at
  `M > 2B`. Since a real root would equal `±z`, none exists. For `3 − θ` over `t² − 2` the norm is
  `7`, not a rational square, yet `H ≡ θ` is a square in `F_9`; only the lift and the exact check
  expose it. The certificate is `DirectBoundExhausted`.
- **Not certified.** Exhausted prime, machine-prime, or nonresidue budgets are resource outcomes.
  An incidental factor or an unsupported map denominator means that a root exists, so they are not
  nonsquare results either.

### Fixtures

- `x² − 2`, `N = 527`, relation `(3, −2)`: `z = 4 + 2θ`, `p = 3`, `e = 4` (`B = 18`, `37 ≤ 81`).
- Non-monic `2x² − 1`, `N = 7`, relation `(3, −4)`: `z = 8 + 4θ`.
- `x⁴ + x − 1`, `N = 17`, relation `(1, −1)` twice: `z = 5 − 3θ + 4θ³` at `p = 3`, `e = 5`.
- `x² + 1`, `N = 5`, relation `(3, −4)`: no real embedding, `z = −2 + 4θ`.
- `x² − 8`: the root `β = θ/2` is outside `Z[θ]`, yet `z = 8`.
- `x² + x + 1` is irreducible modulo 2 but 2 is excluded, reducible modulo 3 (`(t − 1)²`), and
  inert at 5.

### KAT summary

| Test file | Test | Mathematical fact verified |
|-----------|------|---------------------------|
| `sqrt_padic_kat.rs` | `kat_nonconstant_root_monic_quadratic` | `z = 4 + 2θ`, `p = 3`, `e = 4`, `y = 24` for `x² − 2` |
| `sqrt_padic_kat.rs` | `kat_non_monic_root` | `c = 2`: `z = 8 + 4θ`, `β = 1 + θ` |
| `sqrt_padic_kat.rs` | `kat_even_degree_quartic_root` | `d = 4`, `p = 3`, `e = 5`, `z = 5 − 3θ + 4θ³` |
| `sqrt_padic_kat.rs` | `kat_totally_imaginary_root` | `Q(i)` needs no real-sign convention |
| `sqrt_padic_kat.rs` | `kat_root_outside_the_power_basis_order` | `β = θ/2 ∉ Z[θ]` still gives `z = 8` |
| `sqrt_padic_kat.rs` | `kat_precision_boundary_just_above_and_below_twice_the_bound` | `M = 2B + 1` suffices, `M = 2B − 1` would reconstruct the wrong sign |
| `sqrt_padic_kat.rs` | `kat_local_nonsquare_certificate` | Nonresidue in `F_9`: `InertLocalNonresidue` |
| `sqrt_padic_kat.rs` | `kat_global_nonsquare_certificate_after_full_bound` | Local square, global nonsquare: `DirectBoundExhausted` |
| `sqrt_padic_kat.rs` | `kat_prime_search_exhaustion_is_a_resource_outcome` | `V4` field has no inert prime; budget end is not a certificate |
| `sqrt_padic_kat.rs` | `kat_machine_prime_range_end_is_a_resource_outcome` | No prime above `2^64 − 59` |
| `sqrt_padic_kat.rs` | `kat_machine_prime_near_u64_max_reconstructs` | The largest `u64` prime works as an inert prime |
| `sqrt_padic_kat.rs` | `kat_nonresidue_budget_exhaustion_is_a_resource_outcome` | A zero nonresidue budget is a resource outcome |
| `sqrt_padic_kat.rs` | `kat_precision_limit_exhaustion_is_a_resource_outcome` | `bits(2B + 1)` over budget is refused before prime search |
| `sqrt_padic_kat.rs` | `kat_two_is_never_examined` | `p = 2` is excluded and never counted |
| `sqrt_padic_kat.rs` | `kat_primes_dividing_the_modulus_or_leading_coefficient_are_excluded` | `p \| N` and `p \| c` are skipped (`d = 1`) |
| `sqrt_padic_kat.rs` | `kat_zero_residue_target_prime_is_skipped` | `H ≡ 0 (mod 3)` skips to `p = 5` |
| `sqrt_padic_kat.rs` | `kat_reconstruction_agrees_with_an_independently_built_square` | `β = ±(1 + θ)` and `gcd(7 ∓ 24, 527) = 17, 31` |
| `padic.rs` unit tests | `newton_doubles_precision_keeping_unit_square_invariant` | `H·u² ≡ 1` modulo `3, 9, 27, 243` and the root `±(4 + 2θ)` |
| `padic.rs` unit tests | `exponent_chain_*`, `lifting_exponent_*`, `centered_*`, `prime_cursor_*`, `screen_*` | Precision chain, least exponent, symmetric range, prime accounting, and prime classification |

---

## 51. What Comes Next

The square-root and assembly stage is the final stage of the GNFS factoring pipeline. The
output — a non-trivial factor of N — closes the GNFS arc.

**G.W (integrative writeup)** will articulate the whole cross-phase pipeline in a single
narrative, connecting polynomial selection through sieving, filtering, linear algebra, and
square root. The L-notation complexity analysis and the Couveignes correctness proof (if
treated as a payoff proof) land in G.W, not here.

**D.B (NFS-DL linear algebra)** generalises the G.E linear algebra from GF(2) to $\mathbb{F}_\ell$.
The `reduce_mod_ideal` method added at G.F.1 is the natural re-consumer: D.B uses the same
residue-field reduction to compute Schirokauer maps in $\mathbb{F}_\ell$. The C-AlgSqrt
contract (frozen at G.F.3) documents the D.B-friendly generalisation paths.

---

## Further Reading (G.F)

1. **Couveignes, J.-M. (1993).** *Computing a square root for the number field sieve.*
   Unpublished manuscript, cited in Lenstra et al. (1993). The original description of the
   CRT-based algebraic square root algorithm used in G.F.3.

2. **Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) (1993).** *The Development of the Number
   Field Sieve.* Springer LNM 1554. The original papers on GNFS, including the square-root
   step and the role of the quadratic-character columns.

3. **Cohen, H. (1993).** *A Course in Computational Algebraic Number Theory.* Springer GTM
   138. Chapter 4 covers the number field sieve, including the algebraic square root and the
   Chinese Remainder Theorem in number fields.

4. **Pomerance, C. (1996).** *A tale of two sieves.* Notices of the AMS, 43(12), 1473–1485.
   An accessible introduction to the GNFS pipeline, including the square-root step and the
   congruence-of-squares identity.

5. **Buhler, J., Lenstra, H. W. Jr., and Pomerance, C. (1993).** *Factoring integers with the
   number field sieve.* In: Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) The Development of
   the Number Field Sieve, LNM 1554. Springer. The original description of the full GNFS
   pipeline, including the square-root stage and the role of the QC columns.

6. **Kleinjung, T., Aoki, K., Franke, J., et al. (2010).** *Factorization of a 768-bit RSA
   modulus.* CRYPTO 2010. The RSA-768 factorisation, which used the full GNFS pipeline
   including the algebraic square root via Couveignes' algorithm.

---

# The GNFS Pipeline End-to-End: An Integrative Chapter

This chapter is the integrative code-tour for the complete GNFS factoring pipeline. The five
preceding chapters (§1–§51) each explained one stage in isolation. This chapter traces the
whole data-flow from N to a factor in a single narrative, names the cross-phase contracts that
connect the stages, verifies the design statement against what was actually implemented, and
derives the L-notation complexity of GNFS as a payoff.

A reader who has read §1–§51 can read this chapter as a synthesis. A reader who has not can
use it as a map: each section names the stage, states its contract, and points to the detailed
chapter.

---

## 52. The Pipeline at a Glance

The GNFS pipeline has five stages. Each stage consumes a contract from the previous stage and
produces a contract for the next. The contracts are the frozen interfaces that allow the stages
to be developed, tested, and reasoned about independently.

```text
N (integer to factor)
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  G.B — Polynomial Selection                                         │
│  Input:  N, degree d                                                │
│  Output: PolyPair (f, g, m, n, skew, factor_base_bounds)           │
│  Contracts: C-PolyPair (frozen G.B.1), C-Score (frozen G.B.2)      │
└─────────────────────────────────────────────────────────────────────┘
    │  C-PolyPair
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  G.C — Sieving                                                      │
│  Input:  PolyPair, sieve region [−A, A] × [1, B]                   │
│  Output: Vec<Relation>                                              │
│  Contracts: C-Relation (frozen G.C.1), C-FactorBase (frozen G.C.1) │
└─────────────────────────────────────────────────────────────────────┘
    │  C-Relation, C-FactorBase
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  G.D — Filtering                                                    │
│  Input:  Vec<Relation>, FactorBase                                  │
│  Output: SparseMatrix (with provenance map)                         │
│  Contract: C-Matrix (frozen G.D.1)                                  │
└─────────────────────────────────────────────────────────────────────┘
    │  C-Matrix
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  G.E — Linear Algebra                                               │
│  Input:  SparseMatrix, Vec<Relation>, PolyPair                      │
│  Output: Vec<KernelVector>                                          │
│  Contract: C-LinAlg (frozen G.E.1)                                  │
└─────────────────────────────────────────────────────────────────────┘
    │  C-LinAlg (KernelVector + expand_provenance seam)
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  G.F — Square Root and Assembly                                     │
│  Input:  KernelVector, SparseMatrix, Vec<Relation>, PolyPair        │
│  Output: Option<BigInt> (a non-trivial factor of N)                 │
│  Contract: C-AlgSqrt (frozen G.F.3)                                 │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
factor of N
```

The pipeline is realised end-to-end in `gnfs/src/sqrt/assembly.rs` (`factor` function), which
chains `rational_sqrt` and `algebraic_sqrt` over a loop of kernel vectors. The KAT in
`gnfs/tests/factor_end_to_end_kat.rs` verifies the full chain for N = 35 = 5 × 7.

---

## 53. Stage 1: Polynomial Selection and the C-PolyPair Contract

### What polynomial selection does

Polynomial selection (G.B, §1–§11) finds a polynomial pair (f, g) such that:

- f ∈ ℤ[x] has degree d ≥ 2, small coefficients, and many roots mod small primes.
- g(x) = x − m is the rational-side polynomial.
- The **shared-root invariant** holds: f(m) ≡ 0 (mod N).

The shared-root invariant is the load-bearing mathematical fact that makes the whole pipeline
work. It means that for any sieve pair (a, b), the rational norm a − bm and the algebraic norm
b^d · f(a/b) are connected: both are norms of the same ideal (a − bα) in ℤ[α], where α is a
root of f. This connection is what allows the linear algebra step to produce a congruence of
squares mod N.

### The C-PolyPair contract (frozen G.B.1, commit 2f43f99)

The `PolyPair` struct is the **C-PolyPair contract** — the interface between polynomial
selection and all downstream stages:

```rust
pub struct PolyPair {
    pub f: IntPoly,                          // algebraic-side polynomial
    pub g: IntPoly,                          // rational-side polynomial g = x − m
    pub m: BigInt,                           // shared root: f(m) ≡ 0 (mod n)
    pub n: BigInt,                           // the integer to factor
    pub degree: usize,                       // degree of f
    pub skew: Option<f64>,                   // Murphy-E skew (set by G.B.2)
    pub factor_base_bounds: Option<(u64, u64)>, // (B_rat, B_alg) (set by G.C)
}
```

The `verify()` method enforces the invariants: f has degree ≥ 1, `degree` is consistent with
f, g = x − m, and f(m) ≡ 0 (mod n). Every generator (base-m, root sieve, Coppersmith) calls
`verify()` before returning a pair.

The `skew` and `factor_base_bounds` fields are `None` at construction and populated by later
stages. This is the over-specification pattern: adding them now is cheaper than amending the
contract after G.C consumes it.

### The C-Score contract (frozen G.B.2, commit 00aa32d)

The `score` function is the **C-Score contract** — the Murphy-E scoring interface:

```rust
pub fn score(pair: &PolyPair) -> f64;
```

Murphy-E predicts the density of smooth relations the sieve will produce. A higher score means
more smooth pairs in the sieve region. The score is the average over a sample of (a, b) pairs
of the product of Dickman-ρ values for the algebraic and rational norms:

```text
E(f, g) ≈ (1/|S|) Σ_{(a,b) ∈ S} ρ(log|F(a,b)| / log B_f) · ρ(log|G(a,b)| / log B_g)
```

The C-Score contract is consumed by the root sieve (G.B.3), Coppersmith multi-poly (G.B.4),
and the D.A NFS-DL polynomial selection. It is also the signal that E.K (factor-base balancing)
uses to choose the optimal factor-base bounds.

**Principle-4 annotation.** Murphy-E's predictive value only manifests at sieve scale
(N ≳ 2^100). At toy scale, the score is a ranking heuristic: the ordering is correct in
expectation, but the absolute values are not meaningful.

### The non-monic seam

Base-m expansion produces a polynomial f whose leading coefficient satisfies 1 ≤ a_d < m, so
f is generally non-monic. The `shared-numfield` crate's `NumberField::new` requires a monic
defining polynomial. `PolyPair` resolves this seam with two methods:

```rust
pub fn monic_f(&self) -> IntPoly;       // f_monic(x) = a_d^{d−1} · f(x / a_d)
pub fn number_field(&self) -> NumberField; // K = ℚ(α) where α is a root of monic_f
```

The sieve uses the original non-monic f for norm computation; the square-root stage uses
`number_field()` for algebraic arithmetic. The seam is internal to `PolyPair`.

---

## 54. Stage 2: Sieving and the C-Relation / C-FactorBase Contracts

### What sieving does

Sieving (G.C, §12–§21) collects smooth pairs (a, b) — pairs for which both the rational norm
a − bm and the algebraic norm b^d · f(a/b) factor completely over their respective factor
bases. Each smooth pair is a **relation**: a row in the eventual GF(2) matrix.

The sieve uses the log-p mark-then-confirm pattern: for each b in 1..=B, initialise a sieve
array of size 2A + 1, mark each position with the log of each factor-base prime that divides
the corresponding norm, and confirm candidates where the accumulated log exceeds a smoothness
threshold. The threshold filter eliminates most candidates before the expensive trial-division
confirmation step.

### The C-FactorBase contract (frozen G.C.1, commit c1dc0b6)

The `FactorBase` struct is the **C-FactorBase contract** — the two-sided factor base:

```rust
pub struct FactorBase {
    pub rational_primes: Vec<u64>,       // primes p ≤ B_rat
    pub algebraic_ideals: Vec<AlgebraicPrime>, // degree-1 prime ideals (p, r)
    pub b_rat: u64,                      // rational smoothness bound
    pub b_alg: u64,                      // algebraic smoothness bound
    pub obstruction_count: usize,        // sign + QC columns (initially 1)
}
```

The algebraic factor base consists of **degree-1 prime ideals** (p, α − r) where r is a root
of f mod p. The sieve condition for the algebraic side is: the ideal (p, r) divides the
algebraic norm N_alg(a, b) if and only if a ≡ r·b (mod p). This is the algebraic analogue of
the rational condition a ≡ b·m (mod p).

The `obstruction_count` field is initialised to 1 (the sign/−1 column). G.E adds quadratic-
character columns; the slot exists so the matrix-width calculation is stable across stages.

**Why degree-1 ideals.** The algebraic norm N_alg(a, b) = b^d · f(a/b) factors in ℤ[α] as a
product of prime ideals. For a prime p with f(r) ≡ 0 (mod p), the ideal (p, α − r) divides
the principal ideal (a − b·α) in ℤ[α], and hence divides N_alg(a, b) = Norm(a − b·α). The
sieve condition a ≡ r·b (mod p) is the efficient way to enumerate which (a, b) pairs have the
ideal (p, r) dividing their algebraic norm.

### The C-Relation contract (frozen G.C.1, commit c1dc0b6)

The `Relation` struct is the **C-Relation contract** — the unit of NFS data:

```rust
pub struct Relation {
    pub a: BigInt,                           // a-coordinate (can be negative)
    pub b: BigInt,                           // b-coordinate (always ≥ 1)
    pub rational_exponents: ExponentVector,  // exponents over rational factor base
    pub algebraic_exponents: ExponentVector, // exponents over algebraic factor base
    pub rational_sign: bool,                 // true if a − b·m < 0
}
```

The `ExponentVector` stores sparse (factor-base index, exponent) pairs with `u32` exponents —
not pre-reduced GF(2) parities. This is the load-bearing cross-track decision: G.D needs
integer exponents to detect exact duplicates; G.E reduces them mod 2 on demand; D.A reads them
mod ℓ for GF(ℓ) linear algebra. Pre-reducing to GF(2) would destroy the information D.A needs.

The `rational_sign` field is the "−1 column" for G.E's linear algebra: the product of all
selected relations' rational norms must be positive (an even number of negative norms). The
sign column encodes this constraint.

The `verify()` method checks four invariants: gcd(a, b) = 1; the rational exponent vector
reconstructs |N_rat(a, b)|; the algebraic exponent vector reconstructs |N_alg(a, b)|; and the
sign matches.

### The three sieve variants

The sieve module provides three relation-collection strategies:

| Variant | Description | Principle-4 annotation |
|---------|-------------|------------------------|
| `line_sieve` | For each b, sieve a ∈ [−A, A] | Asymptotic win under-exposed at toy scale |
| `special_q_sieve` | Restrict to pairs with q | N_alg(a, b); yield multiplier | Yield advantage under-exposed at toy scale |
| `lattice_sieve` | Enumerate the lattice L_q via Gauss-reduced basis | Efficiency gain under-exposed at toy scale |

All three produce `Vec<Relation>` satisfying the C-Relation contract. The choice of variant
affects the yield and efficiency but not the downstream interface.

---

## 55. Stage 3: Filtering and the C-Matrix Contract

### What filtering does

Filtering (G.D, §22–§30) converts the raw relation corpus into a sparse GF(2) matrix suitable
for the linear algebra step. The pipeline is:

```rust
let matrix = build_matrix(&relations, &fb);     // Vec<Relation> → SparseMatrix
let matrix = remove_singletons(matrix);          // fixpoint singleton removal
let matrix = prune_cliques(matrix);              // greedy excess-floor pruning
let matrix = merge_columns(matrix);              // weight-2 and weight-3 XOR merges
```

Each step is a pure function (consumes and returns `SparseMatrix`), making the pipeline
composable and testable in isolation.

### The C-Matrix contract (frozen G.D.1, commit a0e854b)

The `SparseMatrix` struct is the **C-Matrix contract** — the filtered sparse GF(2) matrix with
provenance map:

```rust
pub struct SparseMatrix {
    pub rows: Vec<MatrixRow>,            // rows with column sets and provenance
    pub num_cols: usize,                 // total columns = FactorBase::matrix_width()
    pub obstruction_col_start: usize,    // index of first obstruction column
    pub obstruction_count: usize,        // number of obstruction columns
    pub col_weights: Vec<u32>,           // col_weights[c] = number of rows with 1 in col c
}
```

Each `MatrixRow` carries:
- `cols: Vec<usize>` — sorted column indices where this row has a 1 in GF(2).
- `provenance: Vec<usize>` — sorted, deduplicated original relation indices.

The provenance map is the thread to G.F: it records which original (a, b) pairs contributed
to each filtered-matrix row. G.E does not modify provenance; it passes the final matrix (with
provenance intact) to G.F.

### Why the provenance map is over-specified

The provenance map stores original relation indices, not the GF(2) column set of the merged
row. G.F needs the actual (a, b) pairs to compute the algebraic square root — the product of
(a − b·α) in the number field K = ℚ(α). The GF(2) column set of a merged row is the XOR of
the original exponent parities; it does not contain the original (a, b) values. Storing the
GF(2) column set instead of the original indices would make G.F impossible without re-deriving
the (a, b) pairs from the exponent vectors, which would require keeping the full `Vec<Relation>`
in scope anyway.

### The graph view and singleton removal

The sparse GF(2) matrix has a natural graph interpretation: nodes are the non-obstruction
columns (primes and ideals), and edges are the rows (relations). A GF(2) dependency is a cycle
in this graph. A singleton column (weight 1) is a leaf node — it cannot be part of any cycle.
Singleton removal is therefore a correctness step, not just an optimisation: rows containing
singleton columns contribute nothing to the null space and can be safely discarded.

Singleton removal is a fixpoint iteration: removing a row decrements the weight of every column
it contains, potentially creating new singletons. The loop terminates when no weight-1
non-obstruction column remains.

---

## 56. Stage 4: Linear Algebra and the C-LinAlg Contract

### What linear algebra does

Linear algebra (G.E, §31–§41) finds the left nullspace of the filtered matrix over GF(2). A
left nullspace vector v is a 0/1 indicator over rows: the rows it selects form a subset S of
filtered relations such that, for every prime p in the factor base, the total exponent of p
across all relations in S is even. This is the condition for the product of the corresponding
norms to be a perfect square.

Before running the nullspace solver, G.E fills in the quadratic-character (QC) columns. These
columns are the load-bearing guarantee that the algebraic square root β exists *in K*, not
merely as an ideal square. For each auxiliary prime q that splits completely in K, the QC
column records the Legendre-symbol parity of the algebraic norm at q. A kernel vector that
satisfies all QC column constraints guarantees that the product of the algebraic norms is a
square in K.

### The C-LinAlg contract (frozen G.E.1, commit 416f6db)

The C-LinAlg contract defines three types:

**`BlockVec`** — a block of `BLOCK_WIDTH = 64` GF(2) vectors, packed into a `Vec<u64>`. The
layout is "row of words": `data[i]` is a `u64` whose bit j is the j-th vector's value at row
i. This makes iterating over rows contiguous in memory and keeps the 64 vectors interleaved
bit-by-bit within each word.

**`MatrixOperator`** — a view of a `SparseMatrix` as a linear operator over GF(2). Provides
`apply` (A·V) and `apply_transpose` (Aᵀ·V) for block vectors. Both solvers consume this
interface exclusively; they never read `SparseMatrix` fields directly.

**`KernelVector`** — a vector in the left nullspace of the matrix, represented as a sorted,
deduplicated `Vec<usize>` of filtered-matrix row indices.

The `expand_provenance` method on `KernelVector` is the seam to G.F:

```rust
pub fn expand_provenance(&self, matrix: &SparseMatrix) -> Vec<usize> {
    // Symmetric difference of the provenance sets of all selected rows.
    // Indices that appear an even number of times cancel.
    // ...
}
```

This seam was over-specified at G.E.1 precisely for G.F: it bridges the filtered-matrix row
indices to the original relation indices that G.F needs.

### The two solvers

G.E provides two nullspace solvers:

**Block Lanczos** (`block_lanczos`): builds a Krylov basis for B = A·Aᵀ using a three-term
recurrence. The self-orthogonality winnowing step handles the GF(2)-specific problem that a
nonzero vector v can satisfy v^T B v = 0 (when v is already in the nullspace). The KAT matrix
is deliberately constructed with duplicate rows to exercise this path.

**Block Wiedemann** (`block_wiedemann`): computes a scalar Krylov sequence x^T B^i y and uses
Berlekamp-Massey to find the minimal polynomial of B in the direction (x, y). The kernel vector
is extracted by evaluating f(B)·y via Horner's method.

The key architectural difference: block Lanczos requires a global inner product at each step
(synchronisation barrier); block Wiedemann's Krylov sequences are independent across (x, y)
pairs (no per-step synchronisation). At production NFS scale, block Wiedemann distributes
across multiple machines — the architecture used in the RSA-768 factorisation.

---

## 57. Stage 5: Square Root, Assembly, and the C-AlgSqrt Contract

### What the square-root stage does

The square-root stage (G.F, §42–§51) takes a `KernelVector` from G.E and produces a
non-trivial factor of N. The stage has three sub-steps:

1. **Rational square root.** Expand the kernel vector via `expand_provenance` to recover the
   original relation index set S. Compute the product P = ∏_{i ∈ S} (a_i − b_i·m) over ℤ.
   The GF(2) column-sum being zero guarantees P is a positive perfect square. Compute X =
   isqrt(P) mod N.

2. **Algebraic square root (Couveignes' algorithm).** Compute γ = ∏_{i ∈ S} (a_i − b_i·α)
   in K = ℚ(α). The QC columns guarantee γ is a square in K. Use Chinese remaindering over
   split primes to recover β ∈ K with β² = γ. Resolve the global sign via the real embedding.
   Compute Y = |Norm(β)| mod N.

3. **GCD assembly.** Compute gcd(X − Y, N). If non-trivial, return the factor. Otherwise try
   gcd(X + Y, N). If all kernel vectors yield trivial GCDs, return None.

### The C-AlgSqrt contract (frozen G.F.3, commits c80a855 + ec69a1f)

The C-AlgSqrt contract defines the algebraic square root interface:

```rust
pub fn algebraic_sqrt(
    kv: &KernelVector,
    matrix: &SparseMatrix,
    relations: &[Relation],
    poly: &PolyPair,
) -> BigInt;
```

The contract guarantees: given a valid kernel vector (satisfying the QC column constraints),
`algebraic_sqrt` returns Y such that Y² ≡ Norm(γ) (mod N), where γ = ∏_{i ∈ S}(a_i − b_i·α).

The Couveignes algorithm is the load-bearing implementation:

1. **Form γ.** Compute γ = ∏_{i ∈ S}(a_i − b_i·α) using `NumberFieldElement::mul`.
2. **Select CRT primes.** Find primes p that split completely in K (f has d distinct roots
   mod p). The same `select_qc_primes` machinery from G.E is reused.
3. **Per-prime square root.** For each split prime p with roots r_1, …, r_d of f mod p:
   reduce γ to 𝔽_p at each root, compute the square root via Tonelli–Shanks, and reconstruct
   the unique polynomial β mod p of degree < d via Lagrange interpolation.
4. **CRT lift.** Use Garner's algorithm to lift the per-prime residues to β ∈ ℤ[α].
5. **Sign resolution.** Evaluate β at the real root θ of f (Newton's method from m^{1/d}).
   Choose the sign such that β(θ) > 0.
6. **Compute Y.** Return |Norm(β)| mod N.

The unconditional sign-consistency check in `per_prime_beta` (the G.F.3 review juncture's
finding) is the primary defence against silent failures: a wrong sign at the Lagrange
interpolation step would propagate into the CRT lift and produce a wrong β, yielding a trivial
GCD at assembly rather than a red test.

### The top-level driver

The `factor` function in `gnfs/src/sqrt/assembly.rs` chains the whole pipeline:

```rust
pub fn factor(
    poly: &PolyPair,
    matrix: &SparseMatrix,
    relations: &[Relation],
    kernel_vectors: &[KernelVector],
) -> Option<BigInt> {
    for kv in kernel_vectors {
        let x = rational_sqrt(kv, matrix, relations, poly);
        let y = algebraic_sqrt(kv, matrix, relations, poly);
        if let Some(f) = factor_from_congruence(&x, &y, &poly.n) {
            return Some(f);
        }
    }
    None
}
```

The retry loop is the safety net for trivial GCDs: roughly 50% of kernel vectors yield trivial
GCDs for a semiprime N = p × q (because X ≡ ±Y mod p and X ≡ ±Y mod q are independent events,
each with probability 1/2). The linear algebra step produces many kernel vectors (the nullspace
dimension grows with the excess), so the retry loop terminates quickly.

---

## 58. The Cross-Phase Contracts: A Unified View

The seven cross-phase contracts are the frozen interfaces that allow the GNFS pipeline to be
developed, tested, and reasoned about in stages. This section names them in one place for the
first time.

| Contract | Frozen at | What it defines | Consumed by |
|----------|-----------|-----------------|-------------|
| **C-PolyPair** | G.B.1 (2f43f99) | `PolyPair`: polynomial pair + number-field constructor | G.B.2–4, G.C, G.F, D.A |
| **C-Score** | G.B.2 (00aa32d) | `score(pair) -> f64`: Murphy-E scoring | G.B.3–4, D.A, E.K |
| **C-Relation** | G.C.1 (c1dc0b6) | `Relation`: (a, b) + two exponent vectors + sign bit | G.D, G.E, G.F, D.A |
| **C-FactorBase** | G.C.1 (c1dc0b6) | `FactorBase`: two-sided factor base + sign/QC columns | G.D, G.E, G.F, D.A |
| **C-Matrix** | G.D.1 (a0e854b) | `SparseMatrix`: filtered sparse GF(2) matrix + provenance map | G.E, G.F, D.B |
| **C-LinAlg** | G.E.1 (416f6db) | `BlockVec`, `MatrixOperator`, `KernelVector` + `expand_provenance` | G.E.2–3, G.F, D.B |
| **C-AlgSqrt** | G.F.3 (c80a855 + ec69a1f) | `algebraic_sqrt`: Couveignes CRT algebraic square root | G.F.4, D.B |

### The data-flow in one sentence per contract

**C-PolyPair** carries the polynomial pair from selection to sieving: the sieve evaluates
F(a, b) = b^d · f(a/b) using the original non-monic f, and the square-root stage constructs
K = ℚ(α) using `number_field()`.

**C-Score** carries the Murphy-E score from scoring to the root sieve and Coppersmith: higher
score means more smooth relations, so the search maximises score.

**C-Relation** carries each smooth pair from sieving to filtering: the exponent vectors are
the rows of the GF(2) matrix; the sign bit is the "−1 column"; the integer exponents (not
pre-reduced GF(2) parities) support D.A's GF(ℓ) linear algebra.

**C-FactorBase** carries the column layout from sieving to filtering and linear algebra: the
`matrix_width()` method gives the total column count; the `rational_index` / `algebraic_index`
lookup methods provide the column mapping without re-scanning the factor base.

**C-Matrix** carries the filtered matrix from filtering to linear algebra and square root: the
provenance map is the thread that G.F uses to recover original (a, b) pairs from a null-space
vector.

**C-LinAlg** carries the nullspace substrate from the substrate layer to both solvers and to
G.F: `expand_provenance` is the seam that bridges filtered-matrix row indices to original
relation indices.

**C-AlgSqrt** carries the algebraic square root from G.F.3 to the assembly step: the contract
guarantees Y² ≡ Norm(γ) (mod N) for any valid kernel vector.

### The over-specification pattern

Every contract is over-specified relative to the immediate consumer's needs. The pattern is
consistent across all seven contracts:

- **C-PolyPair** includes `skew` and `factor_base_bounds` fields that are `None` at
  construction and populated by later stages. Adding them at construction is cheaper than
  amending the contract after G.C consumes it.
- **C-Relation** stores integer exponents rather than GF(2) parities. G.E only needs parities,
  but D.A needs the full integers. Pre-reducing would destroy the information D.A needs.
- **C-Matrix** stores original relation indices in the provenance map rather than the GF(2)
  column set of the merged row. G.E only needs the column set, but G.F needs the original
  (a, b) pairs.
- **C-LinAlg** stores row indices in `KernelVector` rather than a bit-mask. G.F needs row
  indices to look up provenance; a bit-mask would require a scan.

The over-specification rule is: when a contract is consumed by multiple downstream stages with
different needs, over-specify for the most demanding consumer. The cost of over-specification
is a slightly larger data structure; the cost of under-specification is a destructive reshard
of a frozen contract.

---

## 59. Design-Statement Verification

The ROADMAP names G.W as "the moment where the design statement is verified against the actual
implementation." The design statement has three scoping principles. This section walks each
principle against what G.A–G.F actually shipped.

### Principle 1: Algorithmic content complete

**Statement.** The full GNFS algorithm is implemented, not a stub. Every stage — polynomial
selection, sieving, filtering, linear algebra, square root, and assembly — is present and
correct.

**Verification.**

- **G.A (number-field substrate).** `shared/numfield` implements polynomial arithmetic
  (`IntPoly`, `RatPoly`), number fields and elements (`NumberField`, `NumberFieldElement`),
  norm via resultant, trace via companion matrix, ideal representation (two-element primary
  form), Dedekind factorisation, and the Dedekind criterion. The `reduce_mod_ideal` method
  (added at G.F.1) is the seam that Couveignes' algorithm consumes. **Complete.**

- **G.B (polynomial selection).** Base-m expansion, degree selection, Murphy-E scoring
  (Dickman-ρ via 4th-order Runge-Kutta), root sieve (Kleinjung rotation), and Coppersmith
  multi-poly are all implemented and KAT-verified. **Complete.**

- **G.C (sieving).** Two-sided factor base (rational primes + algebraic degree-1 ideals),
  norm bridge (signed BigInt → Uint<4>), line sieve (log-p mark-then-confirm), special-q
  strategy, and lattice sieve (Gauss 2D reduction + enumeration) are all implemented and
  KAT-verified. **Complete.**

- **G.D (filtering).** `build_matrix`, `remove_singletons` (fixpoint), `prune_cliques`
  (greedy excess-floor), and `merge_columns` (weight-2 and weight-3 XOR merges) are all
  implemented and KAT-verified. The provenance map is maintained throughout. **Complete.**

- **G.E (linear algebra).** QC column construction (`populate_qc_columns`), block Lanczos
  (three-term recurrence + self-orthogonality winnowing), block Wiedemann (Krylov sequence +
  Berlekamp-Massey + Horner extraction), and `expand_provenance` are all implemented and
  KAT-verified. **Complete.**

- **G.F (square root and assembly).** Rational square root (`rational_sqrt`), Couveignes CRT
  algebraic square root (`algebraic_sqrt`: form γ, split-prime CRT, Garner lift, real-embedding
  sign resolution), and GCD assembly (`factor`, `factor_from_congruence`) are all implemented
  and KAT-verified. The end-to-end KAT recovers factor 5 from N = 35. **Complete.**

**Verdict: Principle 1 — pass.** The full GNFS algorithm is implemented end-to-end. No stage
is a stub.

### Principle 3: No new engineering optimisations

**Statement.** CADO-NFS / msieve remain dev-only oracles, unchanged. No FFT sieve, no GPU
acceleration, no large-prime variation, no distributed computation.

**Verification.**

- **CADO-NFS / msieve.** Both remain dev-only oracles, gated behind `#[ignore]` KATs. No
  production dependency was added. The oracle KATs (`kat_c_cado_nfs_oracle` in `merge_kat.rs`,
  `kat_c_cado_oracle_n35` in `lanczos_kat.rs`) are ignored in the standard test run.

- **Sieve optimisations.** The line sieve uses the log-p mark-then-confirm pattern (the
  standard NFS sieve). No FFT-based sieve, no bucket sieve, no large-prime variation. The
  special-q and lattice sieves are present as demonstration-fidelity implementations, not
  production optimisations.

- **Linear algebra.** Block Lanczos and block Wiedemann are implemented at demonstration
  fidelity: `BLOCK_WIDTH = 64` (the standard word width), no distributed computation, no
  GPU acceleration.

- **Polynomial selection.** Base-m, root sieve, and Coppersmith multi-poly are present. No
  size-optimisation step (Bai et al. 2014), no alpha-value term in Murphy-E.

**Verdict: Principle 3 — pass.** No new engineering optimisations were added. The
implementation is at demonstration fidelity throughout.

### Principle 4: Scale-only at demonstration fidelity

**Statement.** The implementation is correct but not production-scale. Scale-dependent
phenomena are annotated explicitly; the annotations document the science↔engineering gap
rather than papering over it.

**Verification.** The principle-4 annotations are present throughout the codebase:

| Stage | Phenomenon | Annotation location |
|-------|-----------|---------------------|
| G.B | Murphy-E predictive value under-exposed at toy scale | §5, §7 of this document |
| G.C | Log-sieve asymptotic win under-exposed | §15 |
| G.C | Special-q yield multiplier under-exposed | §16 |
| G.C | Lattice sieve efficiency under-exposed | §17 |
| G.D | Singleton cascade depth under-exposed | §25 |
| G.D | Cavallar merge-ordering gain under-exposed | §26 |
| G.D | EXCESS_FLOOR calibration under-exposed | §25 |
| G.E | Block-width amortisation under-exposed | §34, §39 |
| G.E | Wiedemann parallelism under-exposed | §36, §39 |
| G.F | CRT prime count under-exposed | §44, §49 |
| G.F | Trivial-GCD probability correctly exposed | §46, §49 |

The embedding-sign resolution (G.F.3) is the one phenomenon that is **not** a scale artifact:
it is a correctness obligation present even at toy scale. The G.F.3 review juncture identified
this as the silent-failure locus.

**One discovery: bad primes over-exposed at toy scale.** The `is_bad_prime` flag in
`AlgebraicPrime` (§13) documents a phenomenon that is *over-exposed* at toy scale: bad primes
(primes dividing disc(f)) are marginal at cryptographic scale but unavoidable for hand-picked
polynomials at toy scale. This is a principle-4 annotation in the opposite direction — the
phenomenon is more prominent at toy scale than at real scale. The annotation is present in the
code and in §13 of this document.

**Verdict: Principle 4 — pass.** The implementation is correct at demonstration fidelity.
Scale-dependent phenomena are annotated explicitly throughout. No phenomenon is silently
omitted or misrepresented.

### Design-statement verification summary

**Design-statement verified: pass on all three principles.**

- Principle 1 (algorithmic content complete): **pass** — full GNFS pipeline implemented
  end-to-end, no stubs.
- Principle 3 (no new engineering optimisations): **pass** — CADO-NFS / msieve remain
  dev-only oracles; no FFT sieve, GPU, or distributed computation.
- Principle 4 (scale-only at demonstration fidelity): **pass** — all scale-dependent
  phenomena annotated explicitly; one over-exposed phenomenon (bad primes) noted.

No divergence requiring a corrective follow-on session was found. The pipeline is complete
and KAT-green.

---

## 60. L-Notation Complexity Analysis of GNFS

This section derives the heuristic complexity $L_N[1/3, (64/9)^{1/3}]$ of GNFS. The derivation
is the payoff: it shows *why* the smoothness-probability / sieving-cost tradeoff optimises to
exponent 1/3, and where the constant $(64/9)^{1/3}$ comes from. For the full payoff proof at
textbook depth, see `docs/MATHEMATICS.md` §GNFS (the T.G chapter, to be appended).

### The L-notation

Recall the L-notation from `docs/MATHEMATICS.md` §Prerequisites:

$$L_N[\alpha, c] = \exp\!\left(c \cdot (\log N)^\alpha \cdot (\log \log N)^{1-\alpha}\right)$$

This interpolates between polynomial ($\alpha = 0$) and fully exponential ($\alpha = 1$). The
intermediate regime $0 < \alpha < 1$ is *subexponential*. The GNFS achieves $\alpha = 1/3$,
which is the best known exponent for the general integer factorisation problem.

Write $\log N = n$ for brevity (natural logarithm). Then $L_N[\alpha, c] = e^{c \cdot n^\alpha
\cdot (\log n)^{1-\alpha}}$.

### The three costs

The GNFS complexity is determined by three costs:

1. **Sieving cost.** The sieve must collect enough smooth relations to fill the matrix. The
   number of relations needed is approximately equal to the number of factor-base elements,
   which is $\pi(B) \approx B / \log B$ for a smoothness bound $B$. The sieve examines pairs
   (a, b) in a region of size $\sim B^2$ (the sieve region must be large enough to find enough
   smooth pairs). The sieving cost is therefore:

   $$C_{\text{sieve}} \approx B^2 / \log B \approx B^2$$

   (ignoring logarithmic factors, which are absorbed into the $o(1)$ in the L-notation).

2. **Smoothness probability.** A random integer near $N$ is $B$-smooth with probability
   approximately $\rho(u)$ where $u = \log N / \log B$ (Canfield–Erdős–Pomerance). For the
   algebraic norm $N_{\text{alg}}(a, b) \approx N^{1/d}$ (the norm of a degree-d polynomial
   with coefficients $\sim N^{1/d}$ evaluated at a sieve pair), the smoothness probability is
   $\rho(u')$ where $u' = \log(N^{1/d}) / \log B = n / (d \log B)$. The number of sieve pairs
   needed to find one smooth relation is $\sim \rho(u')^{-2}$ (both sides must be smooth,
   approximately independently). The total sieve cost is:

   $$C_{\text{sieve}} \approx \rho(u')^{-2} \cdot (\text{cost per pair})$$

   where the cost per pair is $O(\log B)$ (the log-sieve step).

3. **Linear algebra cost.** The matrix has $\sim B / \log B$ rows and columns. Structured
   Gaussian elimination (or block Lanczos/Wiedemann) costs $\sim (B / \log B)^2$ operations.
   Ignoring logarithmic factors:

   $$C_{\text{linalg}} \approx B^2$$

### The optimal smoothness bound

The total cost is dominated by the maximum of $C_{\text{sieve}}$ and $C_{\text{linalg}}$. Both
are $\sim B^2$ (up to logarithmic factors), so the total cost is $\sim B^2$ regardless of the
smoothness probability — as long as the sieve can find enough smooth pairs.

The constraint is that the sieve must find enough smooth pairs. The number of smooth pairs in
the sieve region of size $\sim B^2$ is $\sim B^2 \cdot \rho(u')^2$. For this to be at least
$B / \log B$ (the number of relations needed), we need:

$$B^2 \cdot \rho(u')^2 \gtrsim B / \log B \implies B \cdot \rho(u')^2 \gtrsim 1 / \log B$$

The smoothness probability $\rho(u')$ for $u' = n / (d \log B)$ is, for large $u'$:

$$\rho(u') \approx u'^{-u'} = \exp(-u' \log u')$$

Setting $B = L_N[1/3, c] = e^{c \cdot n^{1/3} \cdot (\log n)^{2/3}}$, we get $\log B = c \cdot
n^{1/3} \cdot (\log n)^{2/3}$ and:

$$u' = \frac{n}{d \log B} = \frac{n}{d \cdot c \cdot n^{1/3} \cdot (\log n)^{2/3}}
     = \frac{n^{2/3}}{d \cdot c \cdot (\log n)^{2/3}}
     = \frac{1}{dc} \cdot \left(\frac{n}{\log n}\right)^{2/3}$$

For large $N$, $u' \to \infty$, and $\rho(u') \approx e^{-u' \log u'}$. The smoothness
probability is:

$$\rho(u')^2 \approx e^{-2 u' \log u'}$$

The total cost is:

$$C_{\text{total}} \approx B^2 \cdot \rho(u')^{-2} = e^{2c \cdot n^{1/3} (\log n)^{2/3}}
\cdot e^{2 u' \log u'}$$

To minimise this, we need to choose $c$ and $d$ to balance the two exponentials. The key
observation is that both $B^2 = e^{2c \cdot n^{1/3} (\log n)^{2/3}}$ and $\rho(u')^{-2}$
are of the form $L_N[1/3, \cdot]$. The total cost is minimised when both are equal — i.e.,
when the sieving cost and the smoothness-probability cost are balanced.

### The exponent-1/3 derivation

The exponent $\alpha = 1/3$ arises from the balance condition. Let $B = L_N[\alpha, c]$ for
some $\alpha$ to be determined. The sieving cost is $B^2 = L_N[\alpha, 2c]$. The smoothness
probability for norms of size $\sim N^{1/d}$ is:

$$\rho(u') \approx L_N\!\left[-\alpha, -\frac{1}{dc}\right]$$

(using the Corollary in `docs/MATHEMATICS.md` §Prerequisites with $x = N^{1/d}$ and $y = B$).
The number of sieve pairs needed is $\rho(u')^{-2} = L_N[\alpha, 2/(dc)]$. The total cost is:

$$C_{\text{total}} \approx L_N\!\left[\alpha,\, 2c + \frac{2}{dc}\right]$$

The linear algebra cost is also $B^2 = L_N[\alpha, 2c]$. The total cost is dominated by the
maximum of the sieving and linear algebra costs:

$$C_{\text{total}} \approx L_N\!\left[\alpha,\, \max\!\left(2c,\, 2c + \frac{2}{dc}\right)\right]
= L_N\!\left[\alpha,\, 2c + \frac{2}{dc}\right]$$

To minimise over $c$ (for fixed $\alpha$ and $d$), differentiate with respect to $c$:

$$\frac{d}{dc}\!\left(2c + \frac{2}{dc}\right) = 2 - \frac{2}{dc^2} = 0 \implies c^2 = \frac{1}{d}
\implies c = \frac{1}{\sqrt{d}}$$

At the optimal $c = 1/\sqrt{d}$, the total cost is:

$$C_{\text{total}} \approx L_N\!\left[\alpha,\, \frac{2}{\sqrt{d}} + \frac{2\sqrt{d}}{d}\right]
= L_N\!\left[\alpha,\, \frac{4}{\sqrt{d}}\right]$$

Now minimise over $d$: the cost $4/\sqrt{d}$ decreases as $d$ increases, but the norm size
$N^{1/d}$ also decreases, which changes the smoothness probability. The correct balance
requires the norm size to be $L_N[\alpha, \cdot]$, which means:

$$N^{1/d} = L_N[\alpha, \cdot] \implies \frac{1}{d} \log N = O(n^\alpha (\log n)^{1-\alpha})
\implies d = O\!\left(\frac{n^{1-\alpha}}{(\log n)^{1-\alpha}}\right)$$

For the smoothness probability to be $L_N[\alpha, \cdot]$, we need $u' = \log(N^{1/d}) / \log B$
to be $O(n^{1-2\alpha} / (\log n)^{1-2\alpha})$. The balance condition $u' = O(1)$ (so that
$\rho(u')$ is a constant, not subexponentially small) requires $1 - 2\alpha = 0$, i.e.,
$\alpha = 1/2$. But this gives $C_{\text{total}} = L_N[1/2, \cdot]$, which is the quadratic
sieve complexity — not the GNFS complexity.

The GNFS improvement comes from the **two-sided sieve**: both the rational norm and the
algebraic norm must be smooth. The algebraic norm has size $\sim N^{1/d}$ (much smaller than
$N$ for large $d$), so the smoothness probability is much higher. The optimal degree $d$ is
chosen to make the algebraic norm size $L_N[1/3, \cdot]$, which requires:

$$N^{1/d} = L_N[1/3, \cdot] \implies \frac{n}{d} = O(n^{1/3} (\log n)^{2/3})
\implies d = O\!\left(\frac{n^{2/3}}{(\log n)^{2/3}}\right)$$

This is the optimal degree formula $d \sim (3n / \log n)^{1/3}$ (as implemented in
`optimal_degree`). With this choice, the algebraic norm is $L_N[1/3, \cdot]$-smooth with
probability $L_N[-1/3, \cdot]$, and the total cost is:

$$C_{\text{total}} \approx L_N\!\left[\frac{1}{3},\, 2c + \frac{2}{dc}\right]$$

The optimal $c$ for fixed $d$ is $c = 1/\sqrt{d}$, giving cost $L_N[1/3, 4/\sqrt{d}]$. With
$d \sim (3n/\log n)^{1/3}$, the cost is:

$$C_{\text{total}} \approx L_N\!\left[\frac{1}{3},\, \frac{4}{\sqrt{(3n/\log n)^{1/3}}}\right]$$

In the L-notation, the $n$ and $\log n$ factors in $d$ contribute to the constant $c$ in
$L_N[1/3, c]$. The precise computation (see `docs/MATHEMATICS.md` §GNFS for the full
derivation) gives:

$$C_{\text{total}} = L_N\!\left[\frac{1}{3},\, \left(\frac{64}{9}\right)^{1/3}\right]$$

### Why the exponent is 1/3, not 1/2

The quadratic sieve achieves $L_N[1/2, 1]$ by sieving for smooth values of a single quadratic
polynomial $f(x) = x^2 - N$. The GNFS achieves $L_N[1/3, (64/9)^{1/3}]$ by sieving for smooth
values of *two* polynomials — the rational norm and the algebraic norm — where the algebraic
norm has size $\sim N^{1/d}$ for a carefully chosen degree $d$.

The key improvement is that the algebraic norm is much smaller than $N$ (by a factor of
$N^{1 - 1/d}$), so it is smooth with much higher probability. The optimal degree $d \sim
(3n/\log n)^{1/3}$ balances the smoothness probability against the cost of the linear algebra
step. The exponent $1/3$ is the result of this balance: it is the unique $\alpha$ such that
the sieving cost $L_N[\alpha, \cdot]$ and the smoothness-probability cost $L_N[\alpha, \cdot]$
are both $L_N[1/3, \cdot]$ at the optimal degree.

### The constant $(64/9)^{1/3}$

The constant $(64/9)^{1/3} \approx 1.923$ arises from the precise computation of the optimal
$c$ and $d$. The key steps are:

1. The optimal degree is $d = (3n / \log n)^{1/3}$ (from the balance condition).
2. The optimal smoothness bound is $B = L_N[1/3, c]$ with $c = (8/9)^{1/3}$ (from the
   per-side optimisation).
3. The total cost is $L_N[1/3, 2c + 2/(dc)]$ evaluated at the optimal $c$ and $d$, which
   gives $L_N[1/3, (64/9)^{1/3}]$.

The derivation is heuristic: it assumes that the smoothness probabilities for the rational and
algebraic norms are approximately independent, and that the Dickman-ρ approximation is accurate
in the relevant range. These assumptions are believed to be correct in the asymptotic regime
but are not proved. The constant $(64/9)^{1/3}$ is therefore a heuristic constant, not a
theorem.

For the complete payoff proof — the full derivation with all steps, the precise statement of
the heuristic assumptions, and the comparison with the quadratic sieve — see
`docs/MATHEMATICS.md` §GNFS (the T.G chapter, to be appended). That chapter is the maths-first
sibling to this code-tour: it develops the same derivation for a reader who does not already
know the implementation.

### The L-notation hierarchy in context

The GNFS complexity $L_N[1/3, (64/9)^{1/3}]$ sits in the L-notation hierarchy:

| Algorithm | Problem | Complexity |
|-----------|---------|-----------|
| Baby-step/giant-step | ECDLP | $L_n[1, 1/2] = \sqrt{n}$ (fully exponential) |
| Quadratic sieve | Factoring | $L_N[1/2, 1]$ (subexponential) |
| **GNFS** | **Factoring** | $L_N[1/3, (64/9)^{1/3}]$ **(subexponential, best known)** |
| NFS-DL | DLP in $\mathbb{F}_p^*$ | $L_p[1/3, (64/9)^{1/3}]$ (same exponent) |
| BGJT | DLP in $\mathbb{F}_{p^n}$, small $p$ | $L[0, c]$ (quasi-polynomial) |
| Shor | Factoring + DLP | $O((\log N)^3)$ (polynomial, quantum) |

The exponent $1/3$ has been stable since the early 1990s. Three decades of work have improved
the constant and the engineering — not the exponent. Within the NFS paradigm, the exponent has
converged. Whether a classical algorithm can break the $L[1/3]$ barrier for the general
factoring problem is an open question.

---

## 61. Cross-References

### Per-stage code-tour chapters (this document)

- **G.A (number-field substrate):** `shared/numfield/docs/PEDAGOGY.md` — polynomial arithmetic,
  number fields, norm via resultant, ideal representation, Dedekind factorisation.
- **G.B (polynomial selection):** §1–§11 of this document — base-m construction, Murphy-E
  scoring, root sieve, Coppersmith multi-poly, C-PolyPair and C-Score contracts.
- **G.C (sieving):** §12–§21 of this document — the relation as the unit of NFS data, the
  two-sided factor base, the norm bridge, the line sieve, special-q, lattice sieving,
  C-Relation and C-FactorBase contracts.
- **G.D (filtering):** §22–§30 of this document — the relation corpus as a sparse matrix, the
  graph view, singleton removal, excess and clique pruning, column merging, C-Matrix contract.
- **G.E (linear algebra):** §31–§41 of this document — the filtered matrix as a GF(2) linear
  system, sign and QC columns, block Lanczos, block Wiedemann, the kernel-vector →
  provenance → original-relations thread, C-LinAlg contract.
- **G.F (square root and assembly):** §42–§51 of this document — the kernel vector as a
  congruence of squares, rational square root, Couveignes' algorithm, GCD assembly,
  C-AlgSqrt contract.

### Mathematical textbook

- **`docs/MATHEMATICS.md` §GNFS** (T.G chapter, to be appended) — the maths-first sibling to
  this integrative chapter. Develops the GNFS algorithm mathematically from the
  difference-of-squares idea through the number-field bridge, factor-base smoothness, linear-
  algebra dependency, and square-root recovery. Carries the full L-notation subexponentiality
  derivation as the payoff proof. Cross-references this chapter for the code realisation.

- **`docs/MATHEMATICS.md` §Prerequisites** — the L-notation definition, the Canfield–Erdős–
  Pomerance smooth-number density theorem, and the Dickman-ρ function. The prerequisites for
  the L-notation derivation in §60 above.

- **`docs/MATHEMATICS.md` §On Scale** — the three axes of scale (resource/operational,
  mathematical-dimension, structural), the three couplings, and the honest science↔engineering
  gap. The natural-philosophy context for the principle-4 annotations throughout this document.

---

## 62. KAT Summary (G.W — Integrative)

The integrative chapter does not add new KATs (it is a code-tour, not a new implementation).
The KATs that verify the cross-phase pipeline end-to-end are:

| Test | File | What it verifies |
|------|------|-----------------|
| `kat_a_factor_driver_recovers_known_factor` | `factor_end_to_end_kat.rs` | Full pipeline: N = 35 → factor 5 or 7 |
| `kat_a_factor_from_congruence_known_values` | `factor_end_to_end_kat.rs` | gcd(6 − 1, 35) = 5; congruence of squares identity |
| `kat_b_retry_loop_skips_trivial_and_finds_factor` | `factor_end_to_end_kat.rs` | Retry loop: KV1 trivial, KV2 yields factor 5 |
| `kat_d_end_to_end_provenance` | `merge_kat.rs` | Provenance thread: XOR of original relations' GF(2) column sets equals the merged row |
| `kat_3_round_trip_provenance` | `linalg_substrate_kat.rs` | `expand_provenance` returns the correct symmetric difference |
| `kat_a_cross_validation_with_lanczos` | `wiedemann_kat.rs` | Block Lanczos and block Wiedemann cross-validate on the shared 6×4 matrix |

The end-to-end KAT (`kat_a_factor_driver_recovers_known_factor`) is the pipeline-level
verification: it exercises the full chain from a hand-built `PolyPair` and `Vec<Relation>`
through the square-root stage to a recovered factor. The provenance KATs verify the
cross-stage thread that connects G.D, G.E, and G.F.

---

## Further Reading (G.W — Integrative)

1. **Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) (1993).** *The Development of the Number
   Field Sieve.* Springer LNM 1554. The original papers on GNFS, including the full pipeline
   from polynomial selection through the square-root step.

2. **Buhler, J., Lenstra, H. W. Jr., and Pomerance, C. (1993).** *Factoring integers with the
   number field sieve.* In: Lenstra and Lenstra (eds.), LNM 1554. The original description of
   the full GNFS pipeline, including the L-notation complexity analysis.

3. **Pomerance, C. (1996).** *A tale of two sieves.* Notices of the AMS, 43(12), 1473–1485.
   An accessible introduction to the quadratic sieve and the number field sieve, with the
   L-notation comparison.

4. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective.*
   2nd ed. Springer. Chapter 6 covers the full GNFS pipeline, including the L-notation
   complexity analysis and the algebraic square root.

5. **Kleinjung, T., Aoki, K., Franke, J., et al. (2010).** *Factorization of a 768-bit RSA
   modulus.* CRYPTO 2010. The RSA-768 factorisation — the largest GNFS factorisation to date
   — which used the full pipeline described in this chapter.

6. **Granville, A. (2008).** *Smooth numbers: computational number theory and beyond.* MSRI
   Publications 44. The smooth-number density theorem (Canfield–Erdős–Pomerance) and the
   Dickman-ρ function — the mathematical engine of the L-notation derivation in §60.

7. **`docs/MATHEMATICS.md` §GNFS** (T.G chapter, to be appended). The maths-first sibling to
   this chapter: the full L-notation subexponentiality derivation as a payoff proof, the
   Couveignes-correctness sketch, and the GNFS algorithm developed mathematically from first
   principles.

---

# The NFS-DL Pipeline End-to-End: An Integrative Chapter

This chapter is the integrative code-tour for the complete NFS-DL (discrete-logarithm via the
Number Field Sieve) pipeline. The four preceding Track-D sessions (D.A–D.C) each implemented
one stage in isolation. This chapter traces the whole data-flow from $(g, h, p, \ell)$ to
$\log_g(h) \bmod \ell$ in a single narrative, names the cross-phase contracts that connect the
stages, compares the NFS-DL pipeline explicitly with the NFS-factoring pipeline it parallels,
and annotates the scale-dependent phenomena at demonstration fidelity.

A reader who has read the G.W integrative chapter (§52–§62) can read this chapter as the
one-phase-over analogue: every section has a factoring counterpart, and the comparison is
explicit. A reader who has not can use this chapter as a map: each section names the stage,
states its contract, and points to the source.

For the mathematics — the Schirokauer-map algebra, the F_ℓ linear-algebra theory, and the
L-notation complexity of NFS-DL — see `docs/MATHEMATICS.md` ch. 9 (the T.D chapter, to be
appended). This chapter is code-first: it assumes the reader knows the mathematics and shows
how the code realises it.

---

## 63. The NFS-DL Pipeline at a Glance

The NFS-DL pipeline has four stages. Each stage consumes a contract from the previous stage
and produces a contract for the next. The contracts are the frozen interfaces that allow the
stages to be developed, tested, and reasoned about independently.

```text
(g, h, p, ℓ)  — generator, target, prime modulus, subgroup order
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  D.A — DL Relation Collection                                       │
│  Input:  PolyPair, FactorBase, sieve config, ℓ, Schirokauer ideals  │
│  Output: DLMatrix (DLRelation rows: exponent vecs + Schirokauer cols)│
│  Contracts: C-DLRelation (frozen D.A.1), C-Schirokauer (frozen D.A.1)│
└─────────────────────────────────────────────────────────────────────┘
    │  C-DLRelation (DLMatrix)
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  D.B — F_ℓ Linear Algebra                                           │
│  Input:  DLMatrix, ℓ                                                │
│  Output: VirtualLogTable (log_g(p_i) mod ℓ for each factor-base p_i)│
│  Contract: C-LinAlgFl (frozen D.B.1, extended D.B.2)               │
└─────────────────────────────────────────────────────────────────────┘
    │  C-LinAlgFl (VirtualLogTable)
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│  D.C — Individual-Logarithm Descent                                 │
│  Input:  g, h, p, ℓ, VirtualLogTable, PolyPair, FactorBase         │
│  Output: log_g(h) mod ℓ                                             │
│  Contracts: C-Descent (frozen D.C.1), C2 (frozen D.C.3)            │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
log_g(h) mod ℓ
```

The pipeline is realised end-to-end in `solve_dl_full` (in `gnfs/src/dl/descent/solve.rs`),
which chains initialization-smoothing, `run_descent`, and `assemble_log` over a single call.
The frozen C2 interface `solve_dl` is the cross-track entry point that E.C (the MOV bridge)
will call. The end-to-end KAT in `gnfs/tests/dl_individual_log_kat.rs` verifies the full chain
for $p = 11$, $g = 2$, $\ell = 5$.

---

## 64. Stage 1: DL Relation Collection and the C-DLRelation / C-Schirokauer Contracts

### What DL relation collection does

DL relation collection (D.A, frozen D.A.1) takes smooth factoring relations from the sieve and
augments each with **Schirokauer virtual-log columns** to produce the DL relation matrix that
D.B's F_ℓ linear algebra will consume.

The key difference from NFS-factoring relation collection is the augmentation step. In
NFS-factoring, a smooth relation (a, b) contributes a row of GF(2) exponent parities to the
matrix. In NFS-DL, the same smooth relation contributes a row of **integer exponents mod ℓ**
(the rational and algebraic factor-base columns) plus **Schirokauer columns** (the virtual-log
correction terms that resolve the unit-group obstruction). The Schirokauer columns are the
algebraic heart of the DL setting — see §68 for the explicit comparison.

### The C-DLRelation contract (frozen D.A.1, commit f2dbf0a)

The `DLRelation` struct is the **C-DLRelation contract** — the unit of NFS-DL data:

```rust
pub struct DLRelation {
    /// The factoring relation (u32 exponent vectors, rational + algebraic sides).
    /// Reused directly from the NFS-factoring pipeline.
    pub relation: Relation,

    /// Schirokauer virtual-log columns: one BigInt per prime ideal in the Schirokauer set.
    /// Each entry is λ_i(β) ∈ ℤ/ℓ, the ℓ-adic virtual-log coordinate for the i-th ideal.
    pub schirokauer_cols: Vec<BigInt>,
}
```

The `DLRelation` wraps the factoring `Relation` rather than re-typing it. This is the
over-specification pattern from C-Relation: the factoring relation stores integer exponents
(not pre-reduced GF(2) parities), so D.A can read them directly mod ℓ without resharding.

The `DLMatrix` assembles a collection of `DLRelation` rows with metadata:

```rust
pub struct DLMatrix {
    pub relations: Vec<DLRelation>,
    pub rational_size: usize,    // rational exponent columns
    pub algebraic_size: usize,   // algebraic exponent columns
    pub schirokauer_rank: usize, // Schirokauer columns (virtual-log dimension r)
}
```

The column layout is:

$$
\text{columns} = \underbrace{\text{rational\_size}}_{\text{rational primes}} +
\underbrace{\text{algebraic\_size}}_{\text{algebraic ideals}} +
\underbrace{\text{schirokauer\_rank}}_{\text{Schirokauer correction}}
$$

The `collect_dl_relations` function is the D.A entry point:

```rust
pub fn collect_dl_relations<'a>(
    poly: &PolyPair,
    fb: &FactorBase,
    sieve_config: &LineSieveConfig,
    nf: &'a NumberField,
    ell: &BigInt,
    schirokauer_ideals: &[PrimeIdeal<'a>],
) -> Vec<DLRelation>
```

It runs `line_sieve` to collect smooth factoring relations, then augments each with
`augment_relation` (which constructs $\beta = a + b\alpha$ and evaluates the Schirokauer map).
Relations for which the Schirokauer map is undefined (degenerate elements) are silently skipped,
matching NFS-DL practice.

### The C-Schirokauer contract (frozen D.A.1, commit f2dbf0a)

The Schirokauer map is the DL-specific algebraic ingredient with no factoring analogue. The
public interface is:

```rust
pub fn schirokauer<'a>(
    elt: &NumberFieldElement<'a>,
    ell: &BigInt,
    ideals: &[PrimeIdeal<'a>],
) -> Result<Vec<BigInt>, SchirokauerError>
```

Re-exported as `compute_schirokauer` at the `dl` module level. The `PrimeIdeal<'a>` type is
`Ideal<'a>` from `shared-numfield`, re-exported under the DL module's preferred name.

For each prime ideal $\varphi_i = (p_i, \alpha - r_i)$ with $p_i \equiv 1 \pmod{\ell}$, the
map computes the ℓ-adic virtual-log coordinate $\lambda_i(\beta) \in \mathbb{Z}/\ell$:

1. $\varepsilon_i = (p_i - 1)/\ell$ (an integer, since $\ell \mid p_i - 1$).
2. Compute $\beta^{\varepsilon_i}$ in $K = \mathbb{Q}[\alpha]/(f(\alpha))$ (exact arithmetic).
3. Reduce $\beta^{\varepsilon_i} \bmod \ell^2$ (coefficient-wise).
4. $\delta = \beta^{\varepsilon_i} - 1$ (subtract 1 from the constant term).
5. Divide each coefficient of $\delta$ by $\ell$ (exact; valid since $\beta^{\varepsilon_i}
   \equiv 1 \pmod{\ell}$ for $\beta$ coprime to $\ell$).
6. Evaluate at $\alpha \equiv r_i \pmod{\ell}$ to get $\lambda_i(\beta) \in \mathbb{Z}/\ell$.

The map is a group homomorphism: $\lambda(\beta\gamma) = \lambda(\beta) + \lambda(\gamma)
\bmod \ell$. This is the defining algebraic property verified by the KAT in
`gnfs/tests/dl_relation_kat.rs`.

The `SchirokauerError` type covers the three failure modes: `RamifiedPrime` (p ≢ 1 mod ℓ),
`ElementDivisibleByEll` (ℓ divides a coefficient of β), and `NotDivisibleByEll` (β^ε − 1 is
not divisible by ℓ, indicating a non-integer-coefficient element or arithmetic error).

---

## 65. Stage 2: F_ℓ Linear Algebra and the C-LinAlgFl Contract

### What F_ℓ linear algebra does

F_ℓ linear algebra (D.B, frozen D.B.1, extended D.B.2) takes the `DLMatrix` from D.A and
solves the linear system over the prime field $\mathbb{F}_\ell$ to recover the **virtual-log
table**: the discrete logarithms $\log_g(p_i) \bmod \ell$ for each factor-base prime $p_i$.

The key difference from NFS-factoring linear algebra is the field. In NFS-factoring (G.E), the
linear system is over GF(2): the exponent parities must sum to zero, and the null-space vector
identifies a congruence of squares. In NFS-DL, the linear system is over $\mathbb{F}_\ell$: the
exponents (reduced mod ℓ) plus the Schirokauer correction columns must sum to zero, and the
solution vector gives the virtual logs. See §68 for the explicit comparison.

### The C-LinAlgFl contract (frozen D.B.1, extended D.B.2, commit 652cfa6)

The C-LinAlgFl contract defines the F_ℓ analogues of the GF(2) types in `gnfs::linalg`. The
module is **parallel** to `gnfs::linalg` — no shared trait is introduced; the GF(2) types
remain frozen and untouched.

**`FlBlockVec<F, L>`** — a block of `FL_BLOCK_WIDTH = 32` F_ℓ vectors, each of length
`num_rows`, stored as `data: Vec<[F; FL_BLOCK_WIDTH]>`. The layout is row-major with vectors
interleaved: `data[row][j]` is the j-th vector's value at `row`. This is the F_ℓ analogue of
GF(2)'s bit-packed `BlockVec` — the layout is identical, but each scalar is a field element
rather than a bit.

```rust
pub struct FlBlockVec<F: Fp<L>, const L: usize> {
    pub data: Vec<[F; FL_BLOCK_WIDTH]>,
    pub num_rows: usize,
}
```

The block width `FL_BLOCK_WIDTH = 32` is smaller than GF(2)'s `BLOCK_WIDTH = 64` because field
elements are larger than bits. The principle-4 annotation applies: at toy scale the blocking
overhead is invisible; at NFS scale this is the cache-friendly inner-loop unit.

**`FlSparseMatrix<F>`** — a sparse F_ℓ matrix in CSR-like format, built from `DLMatrix` by
`build_fl_matrix`. Rational and algebraic exponents (u32) are reduced mod ℓ via `from_u64`;
Schirokauer columns (BigInt, already in ℤ/ℓ) are converted via `bigint_to_fp`.

**`FlMatrixOperator<'a, F, L>`** — a view of `FlSparseMatrix<F>` as a linear operator over
F_ℓ. Provides `apply` (A·V) and `apply_transpose` (Aᵀ·V) for block vectors. Mirrors the GF(2)
`MatrixOperator` interface.

**`FlSolution<F>`** — the solver return type: a coefficient vector (the solution to A·x = 0
over F_ℓ) plus an `is_kernel` flag.

**`VirtualLogTable<F>`** — the virtual-log table extracted from an F_ℓ solver solution:

```rust
pub struct VirtualLogTable<F> {
    /// Virtual logs of rational factor-base primes: rational_logs[i] = log_g(p_i) mod ℓ.
    pub rational_logs: Vec<F>,
    /// Virtual logs of algebraic factor-base ideals: algebraic_logs[i] = log_g(φ_i) mod ℓ.
    pub algebraic_logs: Vec<F>,
}
```

The `recover_virtual_logs` function extracts the virtual-log table from a solver solution:

```rust
pub fn recover_virtual_logs<F: Clone>(
    solution: &FlSolution<F>,
    num_rational: usize,
    num_algebraic: usize,
) -> VirtualLogTable<F>
```

It slices the solution's coefficient vector: columns `0..num_rational` are the rational virtual
logs; columns `num_rational..num_rational+num_algebraic` are the algebraic virtual logs. The
Schirokauer correction columns (the remaining entries) are not stored in the table — they are
correction terms, not logs of factor-base elements.

### The two F_ℓ solvers

D.B provides two solvers, mirroring G.E's block Lanczos and block Wiedemann:

**`block_lanczos_fl`**: builds a Krylov basis for $B = A \cdot A^T$ over F_ℓ using a three-term
recurrence. The self-orthogonality winnowing step handles the F_ℓ-specific problem that a
nonzero vector v can satisfy $v^T B v = 0$ (when v is already in the kernel).

**`block_wiedemann_fl`**: computes a scalar Krylov sequence $x^T B^i y$ over F_ℓ and uses
Berlekamp-Massey to find the minimal polynomial of B in the direction (x, y). The kernel vector
is extracted by evaluating $f(B) \cdot y$ via Horner's method.

The end-to-end DL KAT in `gnfs/tests/dl_end_to_end_kat.rs` verifies both solvers recover the
correct virtual logs for the toy setup ($p = 11$, $g = 2$, $\ell = 5$, factor base $\{2, 3\}$).

---

## 66. Stage 3: Individual-Logarithm Descent and the C-Descent / C2 Contracts

### What individual-logarithm descent does

Individual-logarithm descent (D.C, frozen D.C.1, finalised D.C.3) takes the virtual-log table
from D.B and computes $\log_g(h) \bmod \ell$ for a specific target $h$. This is the stage with
**no NFS-factoring analogue** — it has no counterpart in the G.A–G.F pipeline.

The descent has three sub-steps:

1. **Initialization-smoothing.** Find an exponent $e$ such that $g^e \cdot h \bmod p$ is smooth
   over primes up to a medium bound $B'$. This reduces the problem to computing the log of a
   smooth number.

2. **Special-q descent.** For each medium prime $q$ in the smooth factorisation of $g^e \cdot h$
   that is above the factor-base bound, find a relation that rewrites $\log q$ as a combination
   of logs of smaller primes. Repeat until all primes are factor-base leaves with known virtual
   logs.

3. **Log assembly.** Walk the descent tree from leaves (known virtual logs) up to the root,
   accumulating $\log q = \sum \log(\text{child}) \bmod \ell$ at each node. Then
   $\log_g(h) = (\log_g(g^e \cdot h) - e) \bmod \ell$.

### The C-Descent contract (frozen D.C.1, commit 1651993)

The C-Descent contract defines the descent-tree data structures:

**`DescentTarget`** — a prime or ideal being descended:

```rust
pub enum DescentTarget {
    /// A rational prime p.
    Rational(u64),
    /// An algebraic prime ideal (p, r) with f(r) ≡ 0 (mod p).
    Algebraic { p: u64, r: u64 },
}
```

The `prime()` method returns `p` in both cases — the value that must strictly decrease at each
descent step (the termination invariant).

**`DescentNode<F>`** — a node in the descent tree:

```rust
pub struct DescentNode<F> {
    pub target: DescentTarget,
    pub rewriting_relation: Option<DLRelation>,
    pub children: Vec<DescentNode<F>>,
    pub known_log: Option<F>,
}
```

Leaf nodes have `known_log: Some(F)` (looked up from `VirtualLogTable`) and empty `children`.
Interior nodes have `rewriting_relation: Some(rel)` (the special-q descent relation) and
non-empty `children` (the smaller primes from the relation). The invariant
`target.prime() > child.target.prime()` for all children is the termination guarantee.

**`DescentFrontier<F>`** — a max-heap of pending descent targets, ordered by prime descending:

```rust
pub struct DescentFrontier<F> { /* max-heap of HeapEntry<F> */ }
```

The largest-first ordering is the termination invariant: each descent step pops the largest
prime $q$, finds a relation rewriting $q$ as smaller primes, and pushes those smaller primes
back. Since each step strictly reduces the largest prime, the descent terminates when all
frontier elements are factor-base leaves.

The three key functions:

```rust
pub fn init_descent_frontier<F: Clone>(
    g: &BigInt, h: &BigInt, p: &BigInt,
    medium_bound: u64, max_attempts: u64,
) -> Result<(BigInt, DescentFrontier<F>), InitSmoothingError>

pub fn descend_node<F: Clone>(
    target: &DescentTarget,
    poly: &PolyPair,
    fb: &FactorBase,
    vtable: &VirtualLogTable<F>,
    sieve_cfg: &DescentSieveConfig,
) -> Result<DescentNode<F>, DescentStepError>

pub fn run_descent<F: Clone>(
    frontier: DescentFrontier<F>,
    poly: &PolyPair,
    fb: &FactorBase,
    vtable: &VirtualLogTable<F>,
    sieve_cfg: &DescentSieveConfig,
) -> Result<Vec<DescentNode<F>>, SolveDlError>
```

`init_descent_frontier` iterates $e$ from 0 to `max_attempts − 1`, computing
$g^e \cdot h \bmod p$ at each step and trial-dividing by all primes up to `medium_bound`. On
success it returns the smoothing exponent $e$ and the initial frontier (one `DescentNode` per
prime factor, with multiplicity).

`descend_node` finds a rewriting relation for a single medium prime $q$ via the special-q sieve.
It builds a temporary `FactorBase` with $b_\text{alg} \geq q$ (so that $q$ is included as an
algebraic ideal), calls `special_q_sieve` with $q_\text{min} = q_\text{max} = q$, and selects
the first relation that strictly reduces the largest prime. If no such relation is found within
the sieve bounds, it returns `DescentStepError::NoStrictReduction`.

`run_descent` drives the full descent loop: pop the largest frontier prime, call `descend_node`,
push children back, repeat until all primes are factor-base leaves. Returns the flat list of
completed nodes in largest-prime-first order (the order they were popped from the frontier).

### The C2 contract (frozen D.C.3, commit 9d07c51)

The C2 contract is the cross-track interface that E.C (the MOV bridge) will consume:

```rust
pub fn solve_dl(
    g: &BigInt,
    h: &BigInt,
    p: &BigInt,
    k: usize,
    ell: &BigInt,
) -> Result<BigInt, SolveDlError>
```

The `SolveDlError` taxonomy is **frozen at D.C.3**:

```rust
pub enum SolveDlError {
    /// Extension field F_{p^k} (k > 1) is not yet supported.
    Unsupported { k: usize },
    /// Initialization-smoothing failed.
    InitSmoothingFailed { attempts: u64 },
    /// Descent failed: a medium prime could not be rewritten.
    DescentFailed { stuck_prime: u64 },
}
```

**k = 1 (prime field F_p):** The k = 1 path runs initialization-smoothing. The frozen C2
`solve_dl` does not take a `SolveDlContext` (factor base, virtual-log table, sieve config) —
those are not part of the cross-track interface. For the full pipeline, use `solve_dl_full`
with a `SolveDlContext`.

**k > 1 (extension field F_{p^k}):** Returns `SolveDlError::Unsupported { k }` immediately.
The F_{p^k} NFS-DL extension is deferred to an E.C-prep session. This is an
**engineering-scale boundary, not a mathematical one**: the mathematics of NFS-DL over
extension fields is well-understood (see §70 for the principle-4 annotation), but the
implementation is out of D.C scope.

The context-bearing implementation `solve_dl_full` runs the full pipeline:

```rust
pub fn solve_dl_full<F: Clone>(
    g: &BigInt, h: &BigInt, p: &BigInt, k: usize, ell: &BigInt,
    ctx: &SolveDlContext<'_, F>,
) -> Result<BigInt, SolveDlError>
```

where `SolveDlContext<F>` bundles `poly`, `fb`, `vtable`, `sieve_cfg`, and a `to_bigint`
closure for converting virtual logs to `BigInt` for assembly arithmetic.

---

## 67. The Cross-Phase Contracts: A Unified View

The five cross-phase contracts are the frozen interfaces that allow the NFS-DL pipeline to be
developed, tested, and reasoned about in stages.

| Contract | Frozen at | What it defines | Consumed by |
|----------|-----------|-----------------|-------------|
| **C-DLRelation** | D.A.1 (f2dbf0a) | `DLRelation`: factoring `Relation` + Schirokauer cols; `DLMatrix` | D.B |
| **C-Schirokauer** | D.A.1 (f2dbf0a) | `schirokauer` / `compute_schirokauer`, `PrimeIdeal`, `SchirokauerError` | D.A.2, D.B |
| **C-LinAlgFl** | D.B.1 (652cfa6) | `FlBlockVec`, `FlSparseMatrix`, `FlMatrixOperator`, `FlSolution`, `VirtualLogTable`, `recover_virtual_logs` | D.B.2, D.C |
| **C-Descent** | D.C.1 (1651993) | `DescentNode`, `DescentFrontier`, `init_descent_frontier`, `descend_node`, `run_descent` | D.C.2, D.C.3 |
| **C2** | D.C.3 (9d07c51) | `solve_dl(g, h, p, k, ell) -> Result<BigInt, SolveDlError>` | E.C (MOV bridge) |

### The data-flow in one sentence per contract

**C-DLRelation** carries each smooth pair from relation collection to F_ℓ linear algebra: the
integer exponents (not pre-reduced mod ℓ) support D.B's `build_fl_matrix` directly; the
Schirokauer columns are the correction terms that make the DL linear system solvable.

**C-Schirokauer** carries the Schirokauer map from the algebraic substrate to the relation
augmentation step: `compute_schirokauer` evaluates the ℓ-adic virtual-log coordinates for each
prime ideal, and `augment_relation` calls it for each smooth relation.

**C-LinAlgFl** carries the F_ℓ block-solver substrate from the linear-algebra layer to the
descent layer: `VirtualLogTable` is the output of D.B and the input to D.C's descent; the
`recover_virtual_logs` function extracts it from the solver solution.

**C-Descent** carries the descent-tree substrate from the substrate layer to the recursion and
assembly layers: `DescentFrontier` is the work queue; `DescentNode` is the tree node; `run_descent`
drives the full loop; `assemble_log` walks the completed tree.

**C2** carries the NFS-DL solver interface from Track D to Track E: `solve_dl` is the
parameter-minimal frozen API; `SolveDlError` is the frozen error taxonomy.

### The over-specification pattern

The over-specification pattern from the factoring pipeline carries through to NFS-DL:

- **C-DLRelation** wraps the factoring `Relation` rather than re-typing it. The integer
  exponents (not pre-reduced GF(2) parities) are the over-specification that D.B needs.
- **C-LinAlgFl** stores the full coefficient vector in `FlSolution`, not just the virtual-log
  slice. `recover_virtual_logs` extracts the relevant portion; the Schirokauer correction
  columns are available if needed.
- **C-Descent** stores the full `DLRelation` in `DescentNode::rewriting_relation`, not just the
  prime factors. The assembly step uses only the children's logs, but the relation is available
  for debugging and KAT inspection.

---

## 68. The Explicit NFS-Factoring Comparison

This section is the ROADMAP's named requirement for D.W.1: an explicit, structural comparison
of what changes when the target is $\log_g(h)$ rather than a factor of $N$. The comparison is
organised by pipeline stage.

### What stays the same

The NFS-DL pipeline reuses the NFS-factoring infrastructure at every stage:

- **Polynomial selection.** The same `PolyPair` (base-m construction, Murphy-E scoring, root
  sieve, Coppersmith multi-poly) is used without modification. The shared-root invariant
  $f(m) \equiv 0 \pmod{N}$ is the same; the polynomial pair serves both factoring and DL.
- **Sieving.** The same `line_sieve` / `special_q_sieve` / `lattice_sieve` infrastructure
  collects smooth pairs (a, b). The `Relation` type (integer exponent vectors, sign bit) is
  reused directly — the DL-ready over-specification of C-Relation (§54) pays off here.
- **Factor base.** The same `FactorBase` (rational primes + algebraic degree-1 ideals) is used.
  The column layout (rational | algebraic | obstruction) is the same; the DL pipeline adds
  Schirokauer columns after the algebraic columns.

### What changes: the four structural differences

**1. Two number fields and the shared rational side.**

In NFS-factoring, the pipeline uses one algebraic number field $K = \mathbb{Q}(\alpha)$ where
$\alpha$ is a root of the algebraic polynomial $f$. The rational side is the linear polynomial
$g(x) = x - m$.

In NFS-DL, the same two-polynomial structure is used, but the interpretation is different: both
the rational side and the algebraic side contribute to the discrete-log linear system. The
rational norm $a - bm$ and the algebraic norm $b^d \cdot f(a/b)$ are both smooth, and both
contribute exponent columns to the DL matrix. The "two number fields" of NFS-DL are the rational
field $\mathbb{Q}$ (via the rational norm) and the algebraic field $K = \mathbb{Q}(\alpha)$
(via the algebraic norm) — the same two fields as in factoring, but now both contribute to a
single linear system over $\mathbb{F}_\ell$ rather than to a congruence of squares.

**2. Schirokauer maps replacing quadratic-character columns.**

In NFS-factoring (G.E), the obstruction columns are the **sign column** (one column, encoding
the parity of the rational norm's sign) and the **quadratic-character (QC) columns** (several
columns, encoding Legendre-symbol parities). These columns guarantee that the algebraic square
root $\beta$ exists in $K$, not merely as an ideal square. The obstruction is the
**class-group obstruction**: the product of the algebraic norms may be a square in the ideal
group without being a square of a principal ideal.

In NFS-DL, the obstruction is different: the **unit-group obstruction**. The units of $\mathcal{O}_K$
(the ring of integers of $K$) are not accounted for by the exponent vectors alone. The
**Schirokauer map** $\lambda: K^* \to (\mathbb{Z}/\ell)^r$ is the correction that resolves this
obstruction. For each prime ideal $\varphi_i = (p_i, \alpha - r_i)$ with $p_i \equiv 1
\pmod{\ell}$, the map computes the ℓ-adic virtual-log coordinate $\lambda_i(\beta) \in
\mathbb{Z}/\ell$ — the "fractional part" of the discrete log that the exponent vectors miss.

The Schirokauer columns replace the QC columns in the DL matrix. Both serve the same structural
role (resolving an obstruction to the linear system's solvability), but the obstruction is
different (unit group vs class group) and the correction is different (ℓ-adic log vs Legendre
symbol).

**3. F_ℓ (prime-field) linear algebra replacing GF(2).**

In NFS-factoring (G.E), the linear system is over GF(2): find a subset of relations whose
exponent parities sum to zero. The null-space vector identifies a congruence of squares.

In NFS-DL (D.B), the linear system is over $\mathbb{F}_\ell$: find a vector $x$ such that
$A \cdot x = 0$ over $\mathbb{F}_\ell$, where $A$ is the DL relation matrix (exponents reduced
mod ℓ, Schirokauer columns already in $\mathbb{Z}/\ell$). The solution vector gives the virtual
logs: $x_i = \log_g(p_i) \bmod \ell$ for each factor-base prime $p_i$.

The structural difference is the field: GF(2) vs $\mathbb{F}_\ell$. The algorithmic structure
(block Lanczos / block Wiedemann, Krylov subspaces, self-orthogonality winnowing) is the same.
The implementation difference is the scalar type: GF(2) uses bit-packed `u64` words; F_ℓ uses
`[F; FL_BLOCK_WIDTH]` arrays of field elements. The block width changes from 64 (GF(2)) to 32
(F_ℓ) because field elements are larger than bits.

**4. Virtual logs replacing the congruence of squares.**

In NFS-factoring, the output of the linear algebra step is a `KernelVector` — a subset of
relations whose exponent product is a perfect square on both sides. The square-root stage
(G.F) extracts the congruence $X^2 \equiv Y^2 \pmod{N}$ and computes $\gcd(X - Y, N)$.

In NFS-DL, the output of the linear algebra step is a `VirtualLogTable` — the discrete
logarithms $\log_g(p_i) \bmod \ell$ for each factor-base prime $p_i$. There is no square-root
stage: the virtual logs are the answer to the factor-base sub-problem. The individual-logarithm
descent (D.C) then uses the virtual-log table to compute $\log_g(h)$ for a specific target $h$.

### The individual-logarithm special-q descent: the stage with no factoring analogue

The individual-logarithm descent (D.C) is the stage that has **no counterpart** in the
NFS-factoring pipeline. In NFS-factoring, once the linear algebra step produces a kernel vector,
the square-root stage directly computes the factor — there is no "individual target" computation.

In NFS-DL, the linear algebra step gives the virtual logs of the factor-base elements, but not
the log of the specific target $h$. The descent bridges this gap: it finds a path from $h$ down
to the factor base by repeatedly rewriting medium primes as combinations of smaller primes, until
all primes are factor-base leaves with known virtual logs.

The descent is a **tree search** guided by the special-q sieve: for each medium prime $q$, the
sieve finds a smooth relation that contains $q$ as a special prime, rewriting $\log q$ as a
linear combination of logs of smaller primes. The termination invariant (each step strictly
reduces the largest prime) guarantees the tree is finite.

The side-by-side comparison:

| Stage | NFS-factoring | NFS-DL |
|-------|--------------|--------|
| Relation collection | Smooth (a, b) pairs → exponent vectors | Same, plus Schirokauer columns |
| Obstruction columns | Sign + QC columns (class-group correction) | Schirokauer columns (unit-group correction) |
| Linear algebra | GF(2) null space → kernel vector | F_ℓ null space → virtual-log table |
| Output of linalg | Congruence of squares (subset of relations) | Virtual logs of factor-base elements |
| Final step | Square root + GCD assembly | Individual-log descent + log assembly |
| No-factoring-analogue | — | **Special-q descent** |

---

## 69. Design-Statement Note and Principle-4 Annotations

The full design-statement verification for the NFS-DL arc (principles 1/3/4) is D.W.2's
responsibility. This section carries the toy-scale annotations that the code-tour documents.

### Descent-tree breadth and medium-prime tuning

**Principle-4 annotation (descent-tree breadth).** At toy scale ($p = 11$, factor base
$\{2, 3\}$, $\ell = 5$), the descent tree has depth 0 or 1: the initialization-smoothing step
finds $g^e \cdot h$ that is already smooth over the factor base, so no medium primes need to be
descended. The `run_descent` function returns immediately with all frontier nodes as leaves.

At NFS scale ($p \approx 2^{512}$, factor base up to $B \approx 10^6$, medium-prime bound
$B' \approx 10^{12}$), the descent tree has depth $O(\log \log p)$ — typically 3–5 levels. Each
level requires a special-q sieve run to find a rewriting relation. The descent cost is
asymptotically subdominant to the linear algebra cost (see `docs/MATHEMATICS.md` ch. 9 for the
L-notation analysis), but it is the most algorithmically novel stage.

The `DescentSieveConfig` is the scale knob:

```rust
pub struct DescentSieveConfig {
    pub a_bound: u64,    // half-width of the sieve region
    pub b_bound: u64,    // height of the sieve region
    pub threshold_scale: f64, // smoothness threshold (0.5 for descent)
}
```

At toy scale, `a_bound = 10`, `b_bound = 5` suffice. At NFS scale, the bounds are calibrated to
the medium-prime bound and the descent depth.

**Principle-4 annotation (medium-prime tuning).** The medium-prime bound $B'$ (the smoothness
bound for initialization-smoothing) is a critical parameter at NFS scale: too small and the
smoothing fails; too large and the descent tree is too deep. At toy scale, the frozen C2
`solve_dl` uses a hardcoded `medium_bound = 100`; `solve_dl_full` uses the factor-base bound
`ctx.fb.b_alg`. At NFS scale, $B'$ is calibrated to the factor-base bound and the expected
descent depth.

### The C2 F_{p^k} k > 1 `Unsupported` debt

**Engineering-scale boundary (not mathematical).** The `solve_dl` function returns
`SolveDlError::Unsupported { k }` immediately for $k > 1$. This is an **engineering-scale
boundary, not a mathematical one**:

- The mathematics of NFS-DL over extension fields $\mathbb{F}_{p^k}$ is well-understood. The
  algorithm generalises: the number field $K$ is replaced by a degree-$k$ extension; the
  Schirokauer map is adapted to the extension; the descent works over the extended field.
- The implementation gap is that the `PolyPair` / `NumberField` / `FactorBase` infrastructure
  is built for the prime-field case ($k = 1$). Extending to $k > 1$ requires adapting the
  polynomial pair to the extension field structure — a non-trivial engineering task, but not a
  new mathematical idea.

The debt is recorded here as a principle-4 annotation: the code is correct at demonstration
fidelity for $k = 1$; the $k > 1$ path is a known gap, not a silent omission. The E.C-prep
session (Track E) owns the resolution.

---

## 70. KAT Summary (D.W — Integrative)

The integrative chapter does not add new KATs (it is a code-tour, not a new implementation).
The KATs that verify the NFS-DL pipeline end-to-end are:

| Test | File | What it verifies |
|------|------|-----------------|
| `kat_descent_frontier_ordering` | `dl_descent_kat.rs` | `DescentFrontier` pops targets in largest-prime-first order (termination invariant) |
| `kat_init_smoothing_finds_exponent` | `dl_descent_kat.rs` | `init_descent_frontier` finds smooth exponent for toy input (p=101, g=2, h=50) |
| `kat_solve_dl_unsupported_k2` | `dl_descent_kat.rs` | `solve_dl` returns `Unsupported { k: 2 }` for k > 1 |
| `kat_single_node_descent` | `dl_descent_kat.rs` | One medium prime q=17 descends to smaller primes via special-q sieve |
| `kat_multi_level_descent` | `dl_descent_kat.rs` | Two-level descent: frontier {17, 7} descends to factor-base leaves |
| `kat_undescendable_returns_error` | `dl_descent_kat.rs` | Undescendable prime surfaces `DescentFailed` (not a loop) |
| `kat_assembly_correctness` | `dl_descent_kat.rs` | Hand-built descent tree assembles to correct `log_g(h) mod ell` |
| `kat_end_to_end_shape` | `dl_descent_kat.rs` | `solve_dl_full` returns `Ok(2)` for p=11, g=2, h=4, ell=5 |
| `kat_a_fl_arithmetic` | `dl_linalg_kat.rs` | `FlBlockVec` inner-product matrix matches hand-computed F_ℓ values (ℓ=7) |
| `kat_b_matrix_build` | `dl_linalg_kat.rs` | `build_fl_matrix` produces correct F_ℓ values (exponents mod ℓ, Schirokauer cols) |
| `kat_c_block_lanczos_fl` | `dl_linalg_kat.rs` | `block_lanczos_fl` finds kernel vector for known 4×3 F_ℓ system (ℓ=7) |
| `kat_d_block_wiedemann_fl` | `dl_linalg_kat.rs` | `block_wiedemann_fl` cross-validates with Lanczos on the same system |
| `kat_e_virtual_log_recovery` | `dl_linalg_kat.rs` | `recover_virtual_logs` extracts correct rational/algebraic log entries |
| `kat_end_to_end_dl_linalg` | `dl_end_to_end_kat.rs` | Full D.B pipeline: DLMatrix → F_ℓ solve → virtual logs [1, 3] mod 5 for p=11, g=2 |
| `kat_individual_log_h4` | `dl_individual_log_kat.rs` | `solve_dl_full` recovers log_2(4) = 2 mod 5 for p=11, g=2, ell=5 |
| `kat_individual_log_h8` | `dl_individual_log_kat.rs` | `solve_dl_full` recovers log_2(8) = 3 mod 5 |
| `kat_individual_log_cross_check` | `dl_individual_log_kat.rs` | g^result ≡ h (mod p) for all recovered logs (cross-check) |

The end-to-end KATs (`dl_end_to_end_kat.rs`, `dl_individual_log_kat.rs`) are the
pipeline-level verification: they exercise the full chain from a hand-built `DLMatrix` and
`VirtualLogTable` through the descent and assembly stages to a recovered discrete log. The
descent KATs (`dl_descent_kat.rs`) verify the termination invariant, the smoothing step, and
the assembly arithmetic.

---

## 71. Cross-References

### Per-stage code-tour chapters (this document)

- **G.B (polynomial selection):** §1–§11 of this document — base-m construction, Murphy-E
  scoring, root sieve, Coppersmith multi-poly, C-PolyPair and C-Score contracts. The same
  polynomial pair is reused by D.A without modification.
- **G.C (sieving):** §12–§21 of this document — the relation as the unit of NFS data, the
  two-sided factor base, the norm bridge, the line sieve, special-q, lattice sieving,
  C-Relation and C-FactorBase contracts. D.A reuses `line_sieve` and `special_q_sieve`
  directly; the C-Relation over-specification (integer exponents) is the seam D.A consumes.
- **G.E (linear algebra):** §31–§41 of this document — the filtered matrix as a GF(2) linear
  system, block Lanczos, block Wiedemann, C-LinAlg contract. D.B is the F_ℓ analogue: same
  algorithmic structure, different field.
- **G.W (integrative):** §52–§62 of this document — the NFS-factoring pipeline end-to-end.
  The present chapter (§63–§71) is the one-phase-over analogue for NFS-DL.

### Mathematical textbook

- **`docs/MATHEMATICS.md` ch. 9** (T.D chapter) — the maths-first sibling to this integrative
  chapter. Develops the NFS-DL algorithm mathematically: the DLP and the number-field bridge,
  the factor-base / virtual-logarithm construction, Schirokauer maps (proof-sketch depth), F_ℓ
  linear algebra, and the individual-logarithm special-q descent. Carries the L-notation
  complexity of NFS-DL as the designated payoff proof (§9.7). Cross-references this chapter for
  the code realisation.

- **`docs/MATHEMATICS.md` §GNFS** (T.G chapter) — the NFS-factoring maths chapter. The
  L-notation derivation in T.D (§9.7) is a delta on T.G's §7 derivation: same exponent 1/3,
  same constant; the DL delta is descent-cost-subdominance + F_ℓ linear-algebra shape.

- **`docs/MATHEMATICS.md` §On Scale** — the three axes of scale and the honest
  science↔engineering gap. The natural-philosophy context for the principle-4 annotations
  throughout this document.

### Contract definitions

- **C-DLRelation / C-Schirokauer:** `gnfs/src/dl/relation.rs`, `gnfs/src/dl/schirokauer.rs`
  (frozen D.A.1, commit f2dbf0a).
- **C-LinAlgFl:** `gnfs/src/dl/linalg/blockvec_fl.rs` (frozen D.B.1, commit 652cfa6).
- **C-Descent / C2:** `gnfs/src/dl/descent/` (node.rs, recurse.rs, solve.rs) (frozen D.C.1,
  commit 1651993; C2 finalised D.C.3, commit 9d07c51).

---

## Further Reading (D.W — Integrative)

1. **Schirokauer, O. (1993).** *Discrete logarithms and local units.* Philosophical Transactions
   of the Royal Society A, 345(1676), 409–423. The original source for the Schirokauer map —
   the ℓ-adic virtual-log correction that resolves the unit-group obstruction in NFS-DL.

2. **Gordon, D. M. (1993).** *Discrete logarithms in GF(p) using the number field sieve.*
   SIAM Journal on Discrete Mathematics, 6(1), 124–138. The original NFS-DL algorithm for
   prime fields, including the individual-logarithm descent.

3. **Adleman, L. M. (1994).** *The function field sieve.* In: Adleman, L. M., and Huang, M.-D.
   (eds.) Algorithmic Number Theory (ANTS-I), LNCS 877. Springer. The function-field analogue
   of NFS-DL; context for the descent stage.

4. **Joux, A., Lercier, R., Smart, N., and Vercauteren, F. (2006).** *The number field sieve
   in the medium prime case.* In: Dwork, C. (ed.) Advances in Cryptology — CRYPTO 2006, LNCS
   4117. Springer. The medium-prime NFS-DL, including the special-q descent and the
   medium-prime bound calibration.

5. **Barbulescu, R., Gaudry, P., Joux, A., and Thomé, E. (2014).** *A heuristic quasi-polynomial
   algorithm for discrete logarithm in finite fields of small characteristic.* In: Nguyen, P. Q.,
   and Oswald, E. (eds.) Advances in Cryptology — EUROCRYPT 2014, LNCS 8441. Springer. The
   BGJT algorithm — the quasi-polynomial DLP algorithm for small-characteristic fields that
   broke the NFS-DL paradigm for that case.

6. **Kleinjung, T., Aoki, K., Franke, J., et al. (2010).** *Factorization of a 768-bit RSA
   modulus.* CRYPTO 2010. The RSA-768 factorisation — the NFS-factoring counterpart to the
   NFS-DL pipeline described in this chapter.

7. **`docs/MATHEMATICS.md` ch. 9** (T.D chapter). The maths-first sibling to this chapter: the
   full L-notation NFS-DL complexity derivation as a payoff proof (§9.7), the Schirokauer-map
   algebra at proof-sketch depth, and the NFS-DL algorithm developed mathematically from first
   principles.
