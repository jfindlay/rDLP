# rDLP: A Survey of Discrete-Logarithm and Integer-Factorisation Algorithms

rDLP is a complete pedagogical treatment of discrete-logarithm solutions and integer-factorisation
strategies — classical and quantum — organised by a single through-line: **structure-based escape
from search**. Every algorithm in the library finds exploitable algebraic structure that escapes the
generic $\sqrt{n}$ or $L$-notation search bound. The project surveys that escape hierarchy from the
birthday-paradox baseline (Pollard rho, $L[1, 1/2]$) through the number-field sieve
($L[1/3, (64/9)^{1/3}]$) to Shor's quantum period-finding ($L[0]$, polynomial time), with the
algebraic ECDLP attacks filling the middle ground.

The goal is a clean, fast, *readable* Rust implementation where each algorithm demonstrates the
mathematics of its escape at toy scale. Not a record-attempt platform. Not constant-time. Not a
cryptographic library — these are *attacks on* cryptographic problems, run on toy parameters.

---

## The five tracks

### Track ρ — Pollard rho: the birthday-paradox baseline

**Crate:** `rho`  
**What structure it exploits:** the birthday paradox. A pseudorandom walk on a cyclic group of
order $n$ collides in $O(\sqrt{n})$ steps — not because of any special curve structure, but because
of the pigeonhole principle alone. This is the generic lower bound: every algorithm in the library
is measured against it. Track ρ builds the full optimisation stack (Brent's cycle detection,
distinguished points, parallel collision search, negation map, batched inversion, GLV endomorphism)
and establishes the $L[1, 1/2]$ baseline for both integer factorisation and ECDLP.

**Code tour:** [`docs/PEDAGOGY.md`](docs/PEDAGOGY.md) §1–7 (rho phases) and §8–18 (Track E
algebraic attacks, which live in the same `rho` crate).

---

### Track G — GNFS: integer factorisation via the number field sieve

**Crate:** `gnfs` (shared with Track D)  
**What structure it exploits:** *smoothness in two rings simultaneously*. The General Number Field
Sieve maps the factoring problem into a number field $K = \mathbb{Q}(\alpha)$ and collects pairs
$(a, b)$ whose norms are $B$-smooth on both the rational and algebraic sides. Linear algebra over
$\mathbb{F}_2$ then produces a congruence of squares mod $N$. The smoothness density in the number
field — controlled by the choice of polynomial $f$ — is the structural lever that drives the
complexity from $L[1/2]$ (quadratic sieve) down to $L[1/3, (64/9)^{1/3}]$.

**Code tour:** [`gnfs/docs/PEDAGOGY.md`](gnfs/docs/PEDAGOGY.md) §1–71 (polynomial selection,
sieving, filtering, linear algebra, square root).

**Tested end to end:** `gnfs/tests/factor_pipeline_kat.rs` factors $N = 45113 = 197 \cdot 229$ from
relations that `line_sieve` generates at run time (degree 3, $m = 31$, $f = x^3 + 15x^2 + 29x + 8$,
certified irreducible). The relations are validated exactly, filtered, given eight root characters,
and solved with block Lanczos. A verified rational root and a verified algebraic root of the first
computed dependency then give a congruence of squares and a proper divisor. That is one
deterministic fixture at toy scale, run through the checked path (`build_checked_matrix`, root
characters, `DependencyContext`, `factor`); it is not a claim that arbitrary inputs factor or a
statement about production-scale performance. Other cases (a non-monic polynomial, a single
relation with $b \ne 0$, retry past nonsquare candidates) are covered by the square-root and
assembly KATs on hand-built relation sets.

---

### Track D — NFS-DL: discrete logarithm via the number field sieve

**Crate:** `gnfs` (shared with Track G; the `gnfs::dl` module)  
**What structure it exploits:** the same smoothness structure as GNFS, adapted to the
discrete-logarithm problem. The key additions are Schirokauer maps (to handle the unit-group
obstruction in the algebraic side) and individual-logarithm special-$q$ descent (to recover a
specific logarithm from the factor-base logarithms). The complexity is the same $L[1/3]$ class as
GNFS.

**Code tour:** [`gnfs/docs/PEDAGOGY.md`](gnfs/docs/PEDAGOGY.md) (NFS-DL sections, §§9.x in the
mathematical textbook).

---

### Track E — Algebraic ECDLP attacks: structure-specific escapes

**Crates:** `rho` (attack implementations) + `shared/padic` + `shared/gf2m`  
**What structure it exploits:** each attack in Track E finds a *different* algebraic property of
the curve that breaks the generic $\sqrt{n}$ bound:

- **Pohlig–Hellman** — composite group order: factors the ECDLP over prime-order subgroups via CRT.
- **MOV/Frey–Rück** — small embedding degree: a bilinear pairing transports the ECDLP to
  $\mathbb{F}_{p^k}^*$, where index calculus applies.
- **Smart–Satoh–Araki (SSA)** — anomalous curve ($\#E(\mathbb{F}_p) = p$): the $p$-adic formal
  group logarithm gives a polynomial-time solve.
- **GHS/Weil descent** — binary-curve subfield tower: Weil restriction transfers the ECDLP to a
  hyperelliptic Jacobian DLP over $\mathbb{F}_{2^l}$ (a transfer, not an end-to-end solve).
- **Semaev / index calculus** — factor-base decomposability: Semaev summation polynomials build a
  relation matrix whose solution recovers the discrete log.

**Code tour:** [`docs/PEDAGOGY.md`](docs/PEDAGOGY.md) §8–18 (the E.W integrative chapter).

---

### Track S — Shor's algorithm: quantum period-finding

**Crate:** `shor`  
**What structure it exploits:** *quantum periodicity*. The quantum Fourier transform extracts the
period of the modular-exponentiation function in polynomial time — a feat classically impossible
under standard assumptions. This collapses the factoring and discrete-logarithm problems to
polynomial time on a quantum computer, placing them in $L[0]$. Track S implements a classical
state-vector simulator and runs Shor's algorithm (both factoring and ECDLP variants) on toy
parameters, demonstrating the quantum escape without requiring quantum hardware.

**Code tour:** [`shor/docs/PEDAGOGY.md`](shor/docs/PEDAGOGY.md) (§§1–6, three sub-tracks: S.A
state-vector simulator, S.B Shor factoring, S.C Shor ECDLP).

---

## The shared substrate

Six shared crates underpin all five tracks:

| Crate | Role |
|-------|------|
| `shared/field` | Prime-field arithmetic: `Fp` trait, `FpNaive` (schoolbook), `FpMonty` (Montgomery form) |
| `shared/bigint` | Multi-precision integer primitives |
| `shared/numth` | Number-theoretic substrate: Miller–Rabin primality, $B$-smoothness detection, ECM, Tonelli–Shanks |
| `shared/numfield` | Algebraic number-field substrate: polynomial arithmetic, norms, resultants, ideal factorisation (Dedekind) |
| `shared/padic` | $p$-adic arithmetic: Hensel lifting, $p$-adic logarithm — used by the SSA attack in Track E |
| `shared/gf2m` | Binary-field $\mathbb{F}_{2^m}$ arithmetic — used by the GHS attack in Track E |

**Code tours:** [`shared/numth/docs/PEDAGOGY.md`](shared/numth/docs/PEDAGOGY.md) (the α-substrate:
primality, smoothness, ECM) and [`shared/numfield/docs/PEDAGOGY.md`](shared/numfield/docs/PEDAGOGY.md)
(number-field substrate: polynomial arithmetic, norms, Dedekind factorisation).

---

## Workspace crate map

The workspace (`Cargo.toml`) has nine members:

```
gnfs/           Track G (GNFS) + Track D (NFS-DL)
rho/            Track ρ (Pollard rho) + Track E (algebraic ECDLP attacks)
shor/           Track S (Shor's algorithm)
shared/field/   Prime-field arithmetic substrate
shared/bigint/  Multi-precision integer substrate
shared/numth/   Number-theoretic substrate (primality, smoothness, ECM)
shared/numfield/ Number-field substrate (polynomial arithmetic, ideals)
shared/padic/   p-adic arithmetic substrate
shared/gf2m/    Binary-field F_{2^m} arithmetic substrate
```

---

## Deep dives

Two complementary artifacts cover the full library:

- **[`docs/PEDAGOGY.md`](docs/PEDAGOGY.md)** — the master code-first tour. Starts with the Pollard
  rho phase-by-phase implementation (§1–7), then the Track E algebraic attacks (§8–18). Z.1.2 will
  extend this into the library-wide synthesis threading all five tracks.
- **[`docs/MATHEMATICS.md`](docs/MATHEMATICS.md)** — the maths-first textbook. A 12-chapter
  survey (documentation register: undergraduate maths background, survey with proof-sketch depth,
  Markdown + MathJax) covering the escape-from-search through-line (§3), prerequisites (§4), Pollard
  rho (§6), the α-substrate (§7), GNFS (§8), NFS-DL (§9), algebraic ECDLP attacks (§10),
  Shor + post-quantum context (§11), and the modularity speculation (§12).
- **[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md)** — per-phase benchmark results: Phases 1–4 (rho),
  G.B–G.W (GNFS pipeline), the cross-attack ECDLP harness, and the Shor simulator, factoring,
  and ECDLP sections.

---

## Build and test

```sh
# Requires stable Rust (rust-toolchain.toml in rho/ pins to stable; no nightly needed)
cargo build --workspace
cargo test --workspace        # all KATs across all nine crates
cargo build --release --workspace
```

## Benchmarks

```sh
cargo bench                          # all bench targets
cargo bench --bench field            # FpNaive vs FpMonty mul/square/inv
cargo bench --bench factor           # Floyd vs Brent vs Brent+batched-GCD
cargo bench --bench ecdlp            # rho baseline → DPs → negmap → batched-inv → GLV
cargo bench --bench attacks          # algebraic ECDLP attacks (E.W harness)
```

The bench profile uses `opt-level = 3, lto = "thin"`.

## CLIs

**Factor an integer** (Track ρ):

```sh
cargo run --release --bin rho-factor -- <N>
cargo run --release --bin rho-factor -- <N> --threads 4 --batch-size 128
```

**Solve an ECDLP** (`Q = k·P`, returns `k`) (Track ρ):

```sh
# Single-threaded Brent baseline
cargo run --release --bin rho-dlog -- --curve secp-toy --q <x,y>

# Parallel distinguished-point solver
cargo run --release --bin rho-dlog -- --curve secp-toy --q <x,y> --walkers 4 --theta 8

# Add negation map
cargo run --release --bin rho-dlog -- --curve secp-toy --q <x,y> --walkers 4 --negmap

# Add batched inversion (B=16 walks per thread)
cargo run --release --bin rho-dlog -- --curve secp-toy --q <x,y> --walkers 4 --batch-size 16

# Add GLV endomorphism (secp-toy only)
cargo run --release --bin rho-dlog -- --curve secp-toy --q <x,y> --walkers 4 --batch-size 16 --glv
```

Curves: `secp-toy` (63-bit GLV-friendly secp256k1-style) or `generic` (63-bit Weierstrass).
Points are given as `x,y` in decimal with no spaces.

---

## Track ρ: rho crate structure and optimisation inventory

The `rho` crate hosts both Track ρ (Pollard rho) and Track E (algebraic ECDLP attacks).

### Crate layout

```
rho/
  src/
    lib.rs
    field/
      mod.rs            # Fp trait
      naive.rs          # FpNaive: schoolbook mod-p on Uint<4>
      monty.rs          # FpMonty: Montgomery form via DynResidue
    curve/
      mod.rs            # AffinePoint, JacobianPoint, Curve, group law
      generic.rs        # generic Weierstrass y²=x³+ax+b over GF(p)
      secp_k1_toy.rs    # downsized GLV-friendly secp256k1-style curve
      test_curves.rs    # small curves used in unit tests
    factor/
      mod.rs
      rho.rs            # Floyd + Brent + Montgomery batched GCD + multi-c
      cli.rs            # rho-factor binary
    ecdlp/
      mod.rs            # solve_brent, solve_dp, solve_dp_negmap, solve_dp_batch, solve_dp_glv
      walk.rs           # r-adding walk, AffineWalkState, BatchedWalker
      dp.rs             # distinguished-point predicate + DpRecord
      coordinator.rs    # DP hash table, collision detection
      negmap.rs         # canonical representative, FruitlessCycleDetector, BKNS escape
      glv.rs            # GLV endomorphism φ, glv_canonical (6-orbit)
      cli.rs            # rho-dlog binary
    pairing/            # MOV/Frey–Rück pairing reduction (Track E)
    ssa/                # Smart–Satoh–Araki anomalous-curve attack (Track E)
    ghs/                # GHS/Weil descent binary-curve transfer (Track E)
    semaev/             # Semaev summation polynomials (Track E)
    index_calculus/     # Index-calculus DLP solver (Track E)
    pohlig/             # Pohlig–Hellman composite-order reduction (Track E)
    binary_curve/       # Binary-curve group law for GHS (Track E)
    binary_ecdlp/       # Binary-curve ECDLP primitives (Track E)
    hyperelliptic/      # Hyperelliptic Jacobian arithmetic (Track E / GHS)
    util/
      mod.rs
      batch_inv.rs      # Montgomery's batched inversion trick
  benches/
    field.rs            # FpNaive vs FpMonty
    factor.rs           # Floyd vs Brent vs Brent+batched-GCD
    ecdlp.rs            # baseline → DPs → negmap → batched-inv → GLV
    attacks.rs          # algebraic ECDLP attack harness (E.W)
  tests/
    factor_kat.rs       # known semiprimes 15–2^64, all four factor variants
    ecdlp_kat.rs        # group-law KATs, k256 cross-check, ECDLP solver KATs
```

### Optimisation inventory (Track ρ)

**Factoring rho:**

| Optimisation                    | Notes                                              |
|---------------------------------|----------------------------------------------------|
| Floyd's cycle detection         | Baseline; kept for benchmark comparison only       |
| Brent's cycle detection         | ~24% fewer group ops than Floyd                    |
| Polynomial choice `x² + c`      | Reject c ∈ {0, −2}                                 |
| Multi-c parallel restart        | `rayon` over c-values                              |
| Montgomery batched GCD          | Accumulate ∏(xᵢ−xⱼ) mod N over batch ~128, one GCD |
| Pollard–Brent backup on `gcd=N` | Failure-mode recovery                              |

**ECDLP rho:**

| Optimisation                        | Notes                                                 |
|-------------------------------------|-------------------------------------------------------|
| r-adding walk (Teske)               | r = 20; partition by low bits of x-coord              |
| Brent's cycle detection             | Single-threaded baseline; DPs supersede it            |
| Distinguished points (vOW)          | Predicate: low θ bits of x-coord = 0                  |
| Parallel collision search           | N walkers + 1 coordinator; speedup scales as N        |
| DP hash table + match-on-insert     | `HashMap<x_low, DpRecord>`                            |
| Negation map (±P equivalence)       | Canonical rep = point with smaller y (≤ p/2)          |
| Fruitless-cycle detection + escape  | 16-entry sliding window; BKNS deterministic escape    |
| Batched field inversion             | Montgomery's trick: B inversions → 1 inv + 3(B−1) mul |
| Affine coordinates                  | Viable once inversions are amortised via batching     |
| GLV endomorphism (order-3)          | 6-orbit equivalence on secp-toy; ~√3 further speedup  |

---

## Design statement and honest scope

rDLP demonstrates the *mathematics* of each structure-based escape from the search bound at toy
scale. It is not a cryptographically-relevant implementation: the parameters are small enough to
run in seconds on a laptop, the code is not constant-time, and no attack here threatens real-world
key sizes. The project's honest terminus is this: every escape family is surveyed, implemented, and
benchmarked at demonstration fidelity — the mathematics is real, the scale is not. That gap is the
subject of [`docs/MATHEMATICS.md`](docs/MATHEMATICS.md) §5 ("On Scale: A Natural-Philosophy
Interlude"), and it is the honest front-page note for the project as a whole.
