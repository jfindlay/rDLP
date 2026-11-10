# rDLP — Library-Wide Code-Tour: Structure-Based Escape from Search

> **Maths-first sibling.** For the mathematical foundations of every algorithm surveyed here —
> the birthday-paradox collision argument, the $L$-notation hierarchy, the number-field bridge,
> the pairing reduction, the quantum period-finding argument, and the full taxonomy of
> structure-based escapes — see `docs/MATHEMATICS.md`. That file is the maths-first textbook;
> this file is the code-first tour. Neither is a prerequisite for the other.

This document is the library-wide code-tour for the rDLP workspace. It surveys all five algorithm
families of the project — Pollard rho, GNFS, NFS-DL, the algebraic ECDLP attacks, and Shor's
algorithm — as chapters in one story: the story of how each family finds a different kind
of exploitable structure that escapes the generic search bound. The code is the example; this
document is the lesson.

**How to read this document.** The umbrella synthesis (§0) names the through-line once and maps
all five families. The per-family chapters (§§1–7 for Pollard rho, §§8–18 for the algebraic ECDLP
attacks) are the depth for those two families; the remaining three are threaded here as synthesis
chapters that cite their per-crate tours for the implementation detail. The per-crate tours are:

| Family | Code-tour | Sections here |
|--------|-----------|---------------|
| Pollard rho + ECDLP | *this file* §§1–7 | Full depth |
| α — Arithmetic substrate | `shared/numth/docs/PEDAGOGY.md` | §0.3 (synthesis) |
| NF — Number-field substrate | `shared/numfield/docs/PEDAGOGY.md` | §0.3 (synthesis) |
| GNFS + NFS-DL | `gnfs/docs/PEDAGOGY.md` §§1–71 | §0.4 (synthesis) |
| Algebraic ECDLP attacks | *this file* §§8–18 | Full depth |
| Shor's algorithm | `shor/docs/PEDAGOGY.md` | §0.5 (synthesis) |

---

## §0. The Umbrella: Structure-Based Escape from Search

### §0.1 The through-line

Every algorithm in this library is a story about finding exploitable structure — a group
homomorphism, a smoothness phenomenon, an endomorphism, a pairing, a quantum period — that
escapes the generic $\sqrt{n}$ or $L$-notation search bound. This through-line is named once
here and threaded through every chapter below.

The generic bound is the starting point. For a group of order $n$, exhaustive search costs $O(n)$
operations. The birthday paradox gives the first non-trivial improvement: by sampling $O(\sqrt{n})$
random elements, one expects a collision. Pollard's rho algorithm (1978) meets this bound with
$O(\sqrt{n})$ group operations — and for a *generic* group, this is the best any classical
algorithm can do. The $\sqrt{n}$ barrier is the **generic search bound**.

The algorithms in this library escape that barrier by finding structure. The structures fall into
five families (the full taxonomy is in `docs/MATHEMATICS.md` §"Escape from Search: The
Through-Line"):

1. **Group homomorphisms / smoothness** — index calculus, NFS, NFS-DL.
2. **Endomorphisms** — GLV, Koblitz (the rho GLV optimisation).
3. **Pairings** — MOV/Frey–Rück (§10).
4. **Number-field structure** — GNFS, NFS-DL.
5. **Quantum period-finding** — Shor's algorithm.

Each family in this library demonstrates one or more of these escapes in working Rust code, at toy
scale, with known-answer tests. The code is correct; the scale is not cryptographic (see
`docs/MATHEMATICS.md` §"On Scale" for the honest science↔engineering gap).

### §0.2 The comparative L-notation table

The complexity of the best known algorithms for the main problems forms a hierarchy in
$L$-notation:

$$L_N[\alpha, c] = \exp\!\left(c \cdot (\log N)^\alpha \cdot (\log \log N)^{1-\alpha}\right)$$

The table below gives the leading complexity for each algorithm, read from the per-algorithm
derivations in `docs/MATHEMATICS.md` (§6 for ρ, §7 for GNFS, §9.7 for NFS-DL, §10.6 for algebraic
ECDLP, §11.3.4 and §11.4.4 for Shor). **Do not recall these constants from memory — read the cited
sections.**

| Algorithm | Problem | L-notation complexity | Structure exploited | Notes |
|-----------|---------|----------------------|---------------------|-------|
| Pollard rho | ECDLP (generic curve) | $L_n[1,\, 1/2] = \Theta(\sqrt{n})$ | None — generic $\sqrt{n}$ walk | Fully exponential in $\log n$; the baseline all other families escape |
| ECM (Lenstra) | Integer factoring (sub-step) | $L_p[1/2,\, 1]$ (largest prime factor $p$) | Group-order smoothness | Substrate primitive; no other solver in the workspace calls it |
| GNFS | Integer factoring | $L_N[1/3,\, (64/9)^{1/3}]$ | Number-field bridge + smoothness | Exponent $1/3$ from sieve/linear-algebra balance; constant from full analysis [LLMP93] |
| NFS-DL | DLP in $\mathbb{F}_p^*$ | $L_p[1/3,\, (64/9)^{1/3}]$ | Same number-field bridge | Same exponent and constant as GNFS; descent is asymptotically subdominant (§9.7 delta) |
| MOV/Frey–Rück | ECDLP (small embedding degree) | $L_{p^k}[1/3,\, (64/9)^{1/3}]$ (via NFS-DL) | Pairing → $\mathbb{F}_{p^k}^*$ DLP | Cross-family bridge: ECDLP → NFS-DL; precondition $\ell \mid p^k - 1$, $k$ small |
| Smart–Satoh–Araki | ECDLP (anomalous curve) | $L_p[0] = O(\log p)$ (polynomial time) | $p$-adic formal group logarithm | Sharpest escape among the algebraic attacks; precondition $\#E(\mathbb{F}_p) = p$ |
| Index calculus | ECDLP ($E/\mathbb{F}_{p^n}$, $n > 1$) | $L_{p^n}[1/2,\, c]$ (subexponential) | Factor-base decomposability | Asymptotic win requires extension-field setting; toy fixture is $n = 1$ |
| Pohlig–Hellman | ECDLP (composite order) | $L_{p_{\max}}[1,\, 1/2]$ (largest prime factor) | CRT group-order homomorphism | Reduces to prime-order subgroup DLPs; speedup requires $p_{\max} \ll n$ |
| Shor (factoring) | Integer factoring | $O((\log N)^3) = L_N[0,\, 3]$ | Quantum period-finding | Polynomial in $\log N$; dissolves the $L$-notation framework entirely |
| Shor (ECDLP) | ECDLP | $O((\log r)^3) = L_r[0,\, 3]$ (gate count, from the literature [PZ03]) | Quantum period-finding (2-register HSP) | Same polynomial collapse; this project's layout gives 2·256 + 3·256 = 1280 qubits for a 256-bit curve (arithmetic on the layout, not a published estimate; see the ECDLP qubit-layout discussion in docs/MATHEMATICS.md §11.4.4) |

**The quantitative spine of the library.** Reading the table top-to-bottom is the story of the
project: from the generic $\sqrt{n}$ baseline (Pollard rho, $\alpha = 1$) through the subexponential
classical attacks (GNFS/NFS-DL at $\alpha = 1/3$; index calculus at $\alpha = 1/2$; MOV bridging
ECDLP to the $\alpha = 1/3$ regime) through the polynomial-time anomalous-curve attack (SSA,
$\alpha = 0$ classically) to the quantum dissolution of the $L$-notation framework entirely (Shor,
$\alpha = 0$ polynomially). Each step down the $\alpha$ hierarchy is achieved by finding a new kind
of exploitable structure.

**Scale honesty.** The $L$-notation separations in this table are NOT observable at the toy
scale of the project's fixtures. The table reports theoretical complexity classes; the toy-scale
timings are in `docs/BENCHMARKS.md`. The separation between $L[1, 1/2]$ (rho), $L[1/3]$ (GNFS/
NFS-DL/MOV), $L[1/2]$ (index calculus), and $L[0]$ (SSA, Shor) is a statement about asymptotic
behaviour at cryptographic scale — not a ranking of toy-scale timings.

### §0.3 The shared substrate: α and NF chapters

Before the five families diverge, two shared substrate crates provide the arithmetic and algebraic
infrastructure that every family sits on.

**The α-substrate** (`shared/numth/docs/PEDAGOGY.md`). Three crates — `shared-field`,
`shared-bigint`, and `shared-numth` — provide prime-field arithmetic, number-theoretic primitives,
and cross-cutting helpers. The key algorithms: Miller–Rabin primality testing (needed to build
factor bases and to classify an input as prime before factoring it), B-smoothness detection (the
engine of every sieve),
Lenstra's ECM (factoring by group-order smoothness — a direct application of the group-homomorphism
structure; a standalone primitive that no other solver in the workspace calls), and Tonelli–Shanks
square roots (needed for curve-point recovery and the GNFS square-root step). The α-substrate does
not itself escape the generic bound — it provides the tools that other algorithms use to do so. For
the full code-tour, see `shared/numth/docs/PEDAGOGY.md`; for the mathematical development, see
`docs/MATHEMATICS.md` §"The α-Substrate: Primality, Smoothness, and ECM".

**The number-field substrate** (`shared/numfield/docs/PEDAGOGY.md`). The `shared-numfield` crate
provides the algebraic arithmetic of a number field $K = \mathbb{Q}(\alpha)$: polynomial
arithmetic (`IntPoly`, `RatPoly`), number-field elements and the norm map, the resultant and
subresultant GCD, ideal representation, and Dedekind factorisation (including bad primes and the
Dedekind criterion). This substrate is the algebraic engine of GNFS and NFS-DL: both
algorithms exploit the norm map $N_{K/\mathbb{Q}}$ to connect smoothness in $K$ to
smoothness in $\mathbb{Z}$. For the full code-tour, see `shared/numfield/docs/PEDAGOGY.md`; for
the mathematical development, see `docs/MATHEMATICS.md` §GNFS §3 "The Number-Field Bridge".

### §0.4 GNFS and NFS-DL

**The code-tour for GNFS and NFS-DL is `gnfs/docs/PEDAGOGY.md` §§1–71.** This section is the
synthesis chapter; the per-crate tour is the depth.

Both algorithms live in the `gnfs` crate and exploit the same structure: the number-field bridge.
The integers $\mathbb{Z}$ embed into a number field $K = \mathbb{Q}(\alpha)$, and the norm map
$N_{K/\mathbb{Q}}$ connects factorisation in $K$ to factorisation in $\mathbb{Z}$. By sieving for
pairs $(a, b)$ whose norms on both sides are $B$-smooth, collecting enough such *relations*, and
solving a linear system, both algorithms extract the answer — a factor (GNFS) or a discrete
logarithm (NFS-DL).

**GNFS** (`gnfs/docs/PEDAGOGY.md` §§1–62). The General Number Field Sieve factors a
large integer $N$ in five stages: polynomial selection (§§1–11), sieving (§§12–21), filtering
(§§22–30), linear algebra over GF(2) (§§31–41), and square-root recovery (§§42–51). The
integrative chapter (§§52–62) surveys the full pipeline end-to-end. The complexity is
$L_N[1/3, (64/9)^{1/3}]$ — the derivation is the designated payoff proof in
`docs/MATHEMATICS.md` §7.

**NFS-DL** (`gnfs/docs/PEDAGOGY.md` §§63–71). The NFS-DL algorithm adapts GNFS to the
discrete-logarithm problem in $\mathbb{F}_p^*$. The three DL-specific additions are: Schirokauer
maps (resolving the unit-group obstruction that replaces GNFS's class-group obstruction),
$\mathbb{F}_\ell$ linear algebra (replacing GF(2)), and individual-logarithm descent (a
special-$q$ sieve that rewrites each target as a combination of factor-base elements). The
complexity is $L_p[1/3, (64/9)^{1/3}]$ — the same exponent and constant as GNFS, because the
descent is asymptotically subdominant (the delta proof is in `docs/MATHEMATICS.md` §9.7).

**The shared substrate threading GNFS and NFS-DL.** Both algorithms consume the α-substrate
(smoothness detection) and the number-field substrate (polynomial
arithmetic, norm map, ideal factorisation, Dedekind criterion). The `PolyPair` type is the
load-bearing interface: it carries the polynomial pair $(f, g)$ with the shared-root invariant
$f(m) \equiv 0 \pmod{N}$, and is consumed unchanged by the NFS-DL relation collection.

### §0.5 Shor's Algorithm

**The code-tour for Shor's algorithm is `shor/docs/PEDAGOGY.md`.** This section is the synthesis
chapter; the per-crate tour is the depth.

Shor's algorithm is the fifth and final escape from search: quantum period-finding. Where
the classical families find algebraic structure that *reduces* the search cost, Shor's algorithm
finds a quantum period that **dissolves** the $L$-notation bound to polynomial time. The complexity
drops from $L_N[1/3, (64/9)^{1/3}]$ (GNFS) to $O((\log N)^3)$ (Shor factoring) and from
$L_n[1, 1/2]$ (Pollard rho) to $O((\log r)^3)$ (Shor ECDLP): from subexponential or fully
exponential to polynomial in the bit-length.

Shor's algorithm is organised into three parts (`shor/docs/PEDAGOGY.md` §§1–8):

- **The simulator** — the state-vector quantum-circuit simulator (`statevec`, `gates`, `sparse`,
  `measure`, `qft`). The substrate on which Shor's algorithm runs: a classical program that
  faithfully simulates the evolution of a quantum register by maintaining the full
  $2^n$-amplitude array. Its gate set is complete before the algorithms are built on it: neither
  factoring nor ECDLP adds a gate.

- **Factoring** — Shor's factoring algorithm (`arith`, `shor`). The complete order-finding circuit
  (reversible modular exponentiation via permutation synthesis, iQFT, continued-fraction period
  extraction) and the `factor(N)` driver. Factors 15 → {3,5}, 21 → {3,7}, 35 → {5,7}, 91 →
  {7,13} using 8–14 qubits.

- **ECDLP** — Shor's ECDLP algorithm (`curve`, `ecc`, `ecdlp`). The two-register hidden-subgroup
  construction: the $a$-register and $b$-register each hold $t$ qubits; the work register holds
  the running point on $E$. The controlled point-addition circuit (permutation synthesis, ancilla-
  free) evaluates $f(a, b) = a \cdot G + b \cdot Q$ in superposition; the two-dimensional iQFT
  concentrates amplitude near the dual lattice of the hidden subgroup, the pairs $(j/r, l/r)$ with
  $l \equiv k j \pmod{r}$; the extraction rounds the measured pair $(c, d)$ to the nearest
  multiples of $1/r$ and returns the single candidate $k = l \cdot j^{-1} \bmod r$, which the
  solver accepts only after checking $k \cdot G = Q$ and otherwise measures again. Solves the
  ECDLP on a toy curve ($r = 13$) using 17 qubits.

For the mathematical development — the QFT phase-estimation argument, the order-finding →
factoring reduction, the two-register hidden-subgroup proof, and the post-quantum migration
landscape — see `docs/MATHEMATICS.md` ch. 11.

### §0.6 Cross-family connections

Three connections thread the five families into one story.

**The MOV bridge (algebraic ECDLP → NFS-DL).** The MOV/Frey–Rück reduction (§10) is the cross-family
bridge of the library. Given an elliptic curve with small embedding degree $k$
($\ell \mid p^k - 1$), the Weil/Tate pairing
$e: E[\ell] \times E[\ell] \to \mu_\ell \subset \mathbb{F}_{p^k}^*$ transports the ECDLP to a DLP in
$\mathbb{F}_{p^k}^*$, where NFS-DL applies. In code: `mov_reduce` evaluates the reduced Tate pairing
twice, encodes both outputs via the bridge encoder (`fpext_to_bigint`), calls
`gnfs::dl::solve_dl` (the `gnfs` crate's discrete-log entry point), and returns the result only
after checking $k \cdot G = Q$ on the curve. The MOV bridge is the one place
in the library where two families are directly composed at the code level.

**The shared substrate threading GNFS, NFS-DL, and the algebraic ECDLP attacks.** The α-substrate
(smoothness detection) and the number-field substrate (polynomial arithmetic, norm map, ideal
factorisation) are consumed by all three classical-algebraic families. The index-calculus attack
(§§13–14) uses the same smoothness-detection infrastructure as GNFS and NFS-DL; the Semaev summation
polynomial is the ECDLP analogue of the NFS factor base. The substrate is specified once, for all
consumers — the substrate-over-specifies discipline that runs through the whole project.

**The classical → quantum arc (Pollard rho/GNFS/NFS-DL/algebraic ECDLP → Shor).** The five families
form a historical and conceptual arc. Pollard rho (1978) established the $\sqrt{n}$ generic bound
and the first optimisations (distinguished points, negation map, GLV). GNFS (1993) broke the
$L[1/2]$ quadratic-sieve barrier with the number-field bridge, reaching $L[1/3]$. NFS-DL extended
the same structure to discrete logarithms. The algebraic ECDLP attacks showed that the $\sqrt{n}$
ECDLP bound breaks whenever the curve has exploitable structure — pairing, anomalous order, binary
field tower, or factor-base decomposability. Shor's algorithm (1994) dissolved the $L$-notation
framework entirely in the quantum model. The arc is: generic search → subexponential classical →
polynomial quantum. The quantum chapter is the honest coda: the post-quantum migration landscape
(NIST PQC, the SIDH/Castryck–Decru break) shows that structure is a double-edged sword — the same
structural richness that enables escapes also enables attacks on post-quantum candidates.

**The honest terminus.** The whole library demonstrates the *mathematics* of each escape at toy
scale, never at cryptographically relevant scale. The toy fixtures (secp_k1_toy at 63 bits,
GNFS on small semiprimes, NFS-DL at $p = 11$, Shor factoring $N \leq 91$) are correct
demonstrations of the algorithms; they are not attacks. The $L$-notation separations are
asymptotic statements; the toy-scale timings are in `docs/BENCHMARKS.md`. This scale honesty
runs through every chapter.

---

# Pollard Rho: Optimisation by Optimisation

> **Maths-first treatment.** For the mathematical foundations — the birthday-paradox collision
> argument, Floyd's cycle detection, the group-homomorphism structure that makes rho work on
> elliptic curves, orbit-collapsing maps, and the $L$-notation bound — see
> `docs/MATHEMATICS.md §Pollard Rho for ECDLP`. This chapter is the code-tour sibling: it assumes
> the reader knows the mathematics and focuses on the implementation, one optimisation at a time, in
> Rust.

This chapter explains the optimization sequence implemented in the `rho` crate. It is aimed at a
reader who knows group theory and has seen the discrete logarithm problem before, but is new to the
specific algorithmic techniques used here. The code is the example; this chapter is the lesson.

---

## 1. Introduction

### The discrete logarithm problem

Let G be a cyclic group of prime order n, and let G be a generator. Given a target element Q, the
*discrete logarithm problem* (DLP) asks for the integer k ∈ [0, n) such that Q = k·G. On an
elliptic curve over a finite field, this is the *elliptic curve discrete logarithm problem* (ECDLP).

The hardness of the ECDLP is the security foundation of elliptic-curve cryptography. The best known
general-purpose algorithm for solving it runs in O(√n) group operations — not O(n), not O(log n).
That square-root barrier is why a 256-bit elliptic curve provides roughly 128-bit security: an
attacker needs ~2^128 operations, not ~2^256.

### Why Pollard rho gives O(√n) expected time

The key insight is the birthday paradox. If you sample ~√n random elements from a group of size n,
you expect a collision — two samples that land on the same element. Pollard's rho algorithm
(Pollard 1978) exploits this by generating a pseudorandom walk on the group and detecting when the
walk revisits a point. Because the walk is deterministic given its starting state, a revisit implies
a cycle, and the cycle structure encodes the discrete logarithm.

More precisely: the walk maintains the invariant W = a·G + b·Q for tracked scalars a, b ∈ ℤ/nℤ.
When two walk states W₁ = a₁·G + b₁·Q and W₂ = a₂·G + b₂·Q collide at the same group element,
we have (a₁ − a₂)·G = (b₂ − b₁)·Q = (b₂ − b₁)·k·G, so k = (a₁ − a₂)/(b₂ − b₁) mod n (when
b₂ ≠ b₁ mod n). The expected number of steps to collision is O(√n) by the birthday bound.

### What this codebase is

This crate implements two Pollard rho variants — integer factorization and ECDLP — with all
canonical optimizations, built up one optimization at a time. Each stage adds one optimization,
delivers a runnable artifact, and produces a measurable speedup. The goal is a clean, readable CPU
implementation that a researcher can extend or benchmark. It is not a constant-time cryptographic
library; these are *attacks on* cryptographic problems, run on toy parameters.

---

## 2. Optimisation-by-Optimisation Narrative

### Crate skeleton and `FpNaive`

**What it adds:** the crate scaffolding and a schoolbook modular arithmetic implementation.

The `Fp` trait defines the field operations needed throughout: `add`, `sub`, `mul`, `square`,
`inv`, `pow`, `from_u64`, `to_uint`, `zero`, `one`. `FpNaive` implements this trait using
schoolbook arithmetic on `crypto_bigint::Uint<4>` — straightforward reduction after every
operation. This is correct but slow: every multiplication requires a full modular reduction.

**Code reference:** `src/field/naive.rs` (the `FpNaive` implementation).

**Speedup:** baseline. No optimization yet; this stage establishes correctness.

---

### `FpMonty` and the Montgomery-form speedup

**What it adds:** Montgomery-form field arithmetic, the first measurable optimization.

The key insight is that modular reduction is expensive because it requires division. Montgomery's
representation (Montgomery 1985) avoids division by working in a transformed domain: instead of
storing x, store x·R mod p, where R = 2^256 (or another power of 2 chosen for the word size). In
this domain, multiplication of two Montgomery-form elements x·R and y·R gives x·y·R² mod p, which
can be reduced to (x·y·R) mod p using only shifts and additions — no division. The cost is a
one-time conversion into and out of Montgomery form per computation session.

`FpMonty` wraps `crypto_bigint::modular::MontyForm`, which handles the representation internally.
Both `FpNaive` and `FpMonty` implement the same `Fp` trait and pass the same known-answer tests.
The benchmark at `benches/field.rs` shows the speedup directly.

**Code reference:** `src/field/monty.rs` (the `FpMonty` implementation).

**Speedup:** ~2–4× for field multiplication over `FpNaive`, depending on the prime size. This
speedup propagates to every later stage because all curve arithmetic is built on top of it.

---

### Integer factorization rho

**What it adds:** Pollard rho for integer factorization, with Floyd, Brent, and batched-GCD
variants.

The factoring variant uses the pseudorandom map f(x) = (x² + c) mod N. Floyd's cycle detection
runs two pointers (tortoise at speed 1, hare at speed 2) and tests gcd(|x_tortoise − x_hare|, N)
at every step. Brent's improvement freezes the tortoise at the start of each power-of-2 window and
advances only the hare, reducing the number of f-evaluations by ~24% on average.

The batched-GCD trick (Montgomery 1987) amortises the GCD cost further. Instead of computing a GCD
at every step, accumulate the product ∏(xᵢ − y) over a batch of B steps, then take a single GCD
of that product against N. A non-trivial factor is almost always preserved in the product. The only
failure mode is gcd = N (the product is divisible by all prime factors simultaneously), which is
handled by falling back to step-by-step GCD on the offending batch.

**Code reference:** `src/factor/rho.rs` — `floyd`, `brent`, `brent_batched`, and `factor`.

**Speedup:** Brent gives ~24% fewer f-evaluations than Floyd. Batched GCD (batch size ~100) reduces
wall time by another 2–5× by replacing expensive GCD calls with cheap multiplications.

---

### Curve arithmetic

**What it adds:** elliptic curve group law in affine and Jacobian projective coordinates.

An elliptic curve over GF(p) in short Weierstrass form is y² = x³ + ax + b. The group law (point
addition and doubling) is defined by explicit formulas. Affine coordinates are conceptually simple
but require a field inversion per addition — the dominant cost. Jacobian projective coordinates
represent a point (X : Y : Z) with affine coordinates (X/Z², Y/Z³), allowing addition and doubling
with only multiplications and squarings (no inversions), at the cost of a single inversion when
converting back to affine.

Two curves are defined: a generic Weierstrass curve over a ~64-bit prime (used for the
single-threaded baseline through batched inversion), and `secp_k1_toy` — a 63-bit prime-field curve
with a = 0 and an explicit order-3 endomorphism (used for the GLV endomorphism). The secp_k1_toy
curve mirrors the structure of secp256k1 (y² = x³ + 7, j-invariant
0) but over a much smaller prime, making it suitable for pedagogical experiments.

**Code reference:** `src/curve/mod.rs` (group law), `src/curve/secp_k1_toy.rs` (GLV curve
parameters and derivation).

**Speedup:** baseline for ECDLP. The choice of Jacobian vs affine coordinates becomes important with
batched inversion.

---

### ECDLP rho baseline: r-adding walk and Brent's cycle detection

**What it adds:** a working single-threaded ECDLP solver.

The naive pseudorandom walk on an elliptic curve — "add a random point at each step" — has poor
mixing properties. Teske (1998) showed that an *r-adding walk* mixes much better: precompute a
table of r random addends R[0], …, R[r−1], where R[i] = αᵢ·G + βᵢ·Q with random scalars αᵢ, βᵢ.
At each step, select the addend by i = x mod r (where x is the low word of the current point's
x-coordinate), then set W ← W + R[i], a ← a + αᵢ mod n, b ← b + βᵢ mod n. With r ≈ 20, the
walk behaves like a random function on the group, which is the assumption the birthday-bound
analysis requires.

Brent's cycle detection is applied on top of this walk. The tortoise is frozen at the start of each
power-of-2 window; the hare advances one step at a time and is compared to the frozen tortoise at
each step. After r comparisons without collision, the tortoise snaps forward to the hare's position
and the window doubles. This terminates when the hare catches the tortoise on the cycle, using ~24%
fewer walk steps than Floyd's algorithm.

When a collision is detected at walk states (a_t, b_t) and (a_h, b_h), the DLP is recovered as:

```
k = (a_t − a_h) / (b_h − b_t)  mod n
```

The degenerate case b_h = b_t mod n (probability ~1/n) discards the attempt; the next attempt
starts from a fresh random walk table and starting point (attempt `i` uses seed `seed + i`).

**What `solve_brent` returns.** The signature is

```rust
pub fn solve_brent<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>, q: &AffinePoint<F>, n: u64,
    seed: u64, max_retries: usize) -> Result<DlogSolution<RhoEvidence>, DlogError>;
```

The quotient above is only a candidate. The solver validates the instance first, then returns `k`
only after checking `k·G == Q` (`verify_log`), so a success does not depend on the walk:

| Outcome | Meaning | Exercised by |
|---------|---------|--------------|
| `Ok(DlogSolution { k, evidence: Collision { attempt, steps, a1, b1, a2, b2, negated } })` | `k·G == Q`; the evidence identity `a1·G + b1·Q = (−1)^negated · (a2·G + b2·Q)` holds from the fields alone, and a fixed seed reproduces the solution | `brent_collision_evidence_is_rechecked_and_reproducible` (`rho/tests/ecdlp_kat.rs`) |
| `Ok(.. evidence: Trivial)` | `Q = ∞` (`k = 0`) or `Q = G` (`k = 1`); no walk ran | `brent_trivial_targets` |
| `Err(InvalidInput(defect))` | `G = ∞`, `G` or `Q` off the curve, `n` not prime, `n·G ≠ ∞`, or `n·Q ≠ ∞`, found before any walk, in that order | `brent_rejects_off_curve_generator`; `solve_brent_rejects_composite_order_before_walking` (unit, `src/ecdlp/mod.rs`) |
| `Err(BudgetExhausted { attempts, steps })` | all `max_retries` attempts ended degenerate or at the per-attempt step cap; proves nothing about `Q`; `max_retries = 0` gives `{ 0, 0 }` | `brent_zero_retries_is_budget_exhausted` |
| `Err(Inconsistent { k })` | a non-degenerate collision gave a `k` with `k·G ≠ Q`; an invariant failure, never retried | `verified_reports_wrong_candidate_as_inconsistent` (unit; checks the verifier, no fixture reaches it through the walk) |

Each attempt also stops once Brent's window would double past 2²⁸, after 2²⁹ − 1 steps
(`BRENT_MAX_WINDOW` in `src/ecdlp/mod.rs`); that consumes an attempt rather than failing at once.
No test reaches this cap. The input check cannot exclude every bad target: when n² divides the
curve order, a `Q` in `E[n]` outside ⟨G⟩ passes it and surfaces as `Inconsistent` (documented on
`check_prime_order_input`; no test constructs such an instance).

**Code reference:** `src/ecdlp/walk.rs::AddendTable`, `src/ecdlp/walk.rs::WalkState::step`,
`src/ecdlp/mod.rs::solve_brent`, `src/ecdlp/outcome.rs` (`DlogSolution`, `DlogError`,
`RhoEvidence`, `verify_log`, `check_prime_order_input`).

**Speedup:** O(√n) expected walk steps — this is the theoretical baseline. All subsequent
optimizations reduce the constant factor.

---

### Distinguished points and parallel collision search

**What it adds:** the van Oorschot–Wiener (vOW) parallel rho architecture.

Brent's cycle detection is inherently single-threaded: the tortoise and hare must stay in sync.
Distinguished points (van Oorschot and Wiener 1999) break this dependency. A point is
*distinguished* when its x-coordinate has at least θ low-order zero bits — a property that can be
checked independently by any walker. The expected number of steps between distinguished points is
2^θ.

The architecture is a client/server design: N walker threads each run an independent r-adding walk
and emit a `DpRecord` (carrying the point as a `DpPoint` — the full four-word x-coordinate and the
parity of y — and the accumulated scalars a, b) whenever the walk lands on a distinguished point.
A coordinator on the calling thread maintains a hash table keyed on the full x-coordinate (the low
word that decides whether a point is distinguished is not an identity). When two DPs with the same
x-coordinate but different (a, b) pairs arrive, the coordinator solves the linear system for k and
signals all walkers to stop.

Equal x means equal points *or* opposite points, since P and −P share x. For odd p and y ≠ 0, y
and p − y have opposite parity, so the parity bit tells the cases apart:

```
equal points     (parities agree):  a1 + b1·k ≡   a2 + b2·k   k = (a1 − a2) / (b2 − b1)
opposite points  (parities differ): a1 + b1·k ≡ −(a2 + b2·k)  k = −(a1 + a2) / (b1 + b2)
```

The second case is reported as `negated: true` in the evidence, and a pair with b1 + b2 ≡ 0 (or
b1 ≡ b2 in the first) is degenerate and skipped. Only the plain walk can meet −P: the negation-map,
batched and GLV walks canonicalise every point before emitting it, so their collisions are always
between identical points (`dp_collision_evidence_is_rechecked` asserts `negated` is false for
them). Matching by x alone and reading every match as equality was a wrong-answer path: on `tiny_a`
with Q = 1234·G, one walker, θ = 2 and a budget of 2²⁰, 18 of the first 100 seeds ended in a
candidate that failed `k·G == Q`; with the parity bit those runs succeed with `negated: true`
(`plain_walk_recovers_log_from_negated_collision`, which fixes seed 0, and
`plain_walk_never_misreads_the_sign_of_a_collision`, both in `rho/tests/ecdlp_kat.rs`).

**Budget and outcomes.** All four distinguished-point solvers
(`solve_dp`, `solve_dp_negmap`, `solve_dp_batch`, `solve_dp_glv`) return
`Result<DlogSolution<RhoEvidence>, DlogError>` and take a final `max_steps` argument, the total
walk steps across all walkers. Walkers claim steps from one shared counter before taking them, so
the budget is never exceeded (`step_budget_grants_claims_up_to_the_limit`, unit,
`src/ecdlp/mod.rs`); the batched and GLV walkers claim a whole batch at a time, so up to
`batch_size − 1` steps of the budget can go unused. Each solver validates the instance and answers
the trivial targets without starting a walker, and verifies the `k` a collision implies, exactly as
`solve_brent` does:

| Outcome | Meaning | Exercised by |
|---------|---------|--------------|
| `Ok(.. Collision { attempt: 0, steps, .. })` | verified `k`; `attempt` is always 0 for a single run, `steps` is the shared counter when the collision was accepted | `dp_collision_evidence_is_rechecked` |
| `Err(BudgetExhausted { attempts: 1, steps })` | `max_steps` spent without a usable collision; proves nothing about `Q`; the batch variant spends whole batches (a budget of 10 with batch 4 spends 8) | `dp_tiny_budget_is_budget_exhausted`, `dp_glv_tiny_budget_is_budget_exhausted` (`rho/tests/ecdlp_kat.rs`); `dp_solvers_report_budget_exhaustion`, `single_negmap_walker_ends_by_budget_when_it_finds_nothing` (unit, `src/ecdlp/mod.rs`) |
| `Err(InvalidInput(..))` | checked before any walker starts | `dp_trivial_and_invalid_inputs_are_decided_without_a_walk` |
| `Err(Unsupported(..))` | no walkers, a zero batch size, or (`solve_dp_glv`) a curve other than the `secp_k1_toy` group | `dp_solvers_reject_empty_thread_shapes` (unit); `dp_glv_is_unsupported_off_its_curve` |
| `Err(InvariantFailure(..))` | a walker thread panicked or lost its channel | `walker_panic_is_an_invariant_failure`, `lost_channel_is_an_invariant_failure` (unit) |
| `Err(Inconsistent { k })` | a collision gave a `k` with `k·G ≠ Q`; never retried | `verified_reports_wrong_candidate_as_inconsistent` (unit; checks the verifier) |

A search that never finds a collision is bounded by the budget, not by luck: a single walker may
never see one, and `single_negmap_walker_ends_by_budget_when_it_finds_nothing` pins the bounded
outcome (θ above the walk length, so no distinguished point is ever emitted). `solve_dp_glv`
recognises its curve by the parameters p, a, b and the group order, not by name, and returns
`Unsupported` before checking anything else on any other curve, because it applies the
`secp_k1_toy` constants β and λ.

The key insight is that two walks that ever visit the same point will, from that point on, follow
identical trajectories and will therefore hit the same distinguished points. So a collision in the
DP table implies a collision in the underlying walks, which implies the DLP can be recovered.

The θ parameter trades memory against time: larger θ means fewer DPs stored (less memory) but more
steps per walker between DPs (more time before a collision is detected). The textbook choice is
θ such that 2^θ ≈ √n / (N · target-walks), giving ~1000 DPs at solution time.

**Code reference:** `src/ecdlp/dp.rs::is_distinguished`, `src/ecdlp/dp.rs::DpPoint`,
`src/ecdlp/dp.rs::DpRecord`, `src/ecdlp/coordinator.rs`, `src/ecdlp/mod.rs::solve_dp`.

**Speedup:** N× (linear in thread count). This is the second pedagogical moment: the speedup from
parallelism is linear, not √N. Each additional walker independently explores the group, and the
birthday bound applies to the combined pool of DPs. Doubling the number of walkers halves the
expected time to collision.

---

### Negation map and fruitless-cycle escape

**What it adds:** a ~√2 reduction in expected walk length by collapsing {P, −P} pairs.

On an elliptic curve, P and −P share the same x-coordinate (negation flips only the y-coordinate).
The negation map exploits this: instead of treating P and −P as distinct walk states, collapse them
to a single canonical representative. The canonical element is the one with the numerically smaller
y-coordinate (equivalently, y ≤ p/2). After each walk step, apply the canonical map: if the
current point is not canonical, replace it with its negation and flip the signs of the accumulated
scalars (a ← n − a, b ← n − b, since −W = (n−a)·G + (n−b)·Q).

This halves the effective group size: the walk now explores only the ~n/2 canonical representatives
rather than all n group elements. By the birthday bound, halving the group size reduces the
expected number of steps to collision by √2.

The complication is *fruitless cycles*. When the walk oscillates between a point W and its negation
−W, the canonical map maps both to the same representative, so the canonical x-coordinate repeats
every two steps. The walk is trapped in a length-2 cycle and will never hit a distinguished point.
The BKNS escape (Bos, Kleinjung, Niederhagen, Schwabe) handles this deterministically: detect the
period-2 pattern in a sliding window of recent canonical x-coordinates, then perturb by doubling
the current point (W ← 2W, a ← 2a mod n, b ← 2b mod n) and re-canonicalising. The doubling
breaks the cycle because 2W is generically not in the same {W, −W} pair.

**Code reference:** `src/ecdlp/negmap.rs::canonical_rep`, `src/ecdlp/negmap.rs::negate_scalars`,
`src/ecdlp/negmap.rs::FruitlessCycleDetector`, `src/ecdlp/mod.rs::solve_dp_negmap`.

**Speedup:** ~√2 ≈ 1.41× reduction in expected walk steps vs distinguished points alone. The
fruitless-cycle overhead is small in practice (escapes are rare for a well-mixing walk), so the net
speedup is close to the theoretical √2.

---

### Batched field inversion and affine coordinates

**What it adds:** 2–4× reduction in per-step cost by amortising field inversions across a batch of
walks.

In the Jacobian walk (the single-threaded baseline), each step requires one field inversion to
convert the current point to affine coordinates for the partition function. Field inversion is the
most expensive field operation — typically 50–100× the cost of a multiplication. This cost dominates
the per-step budget.

The key insight is that if B walks are advanced in lock-step, their B inversion denominators can be
computed together using Montgomery's batched inversion trick: compute the prefix products
p[0] = d[0], p[1] = d[0]·d[1], …, p[B−1] = d[0]·…·d[B−1]; invert the last prefix product (one
inversion); then recover each individual inverse in a backward pass using 3(B−1) multiplications.
Total cost: 1 inversion + 3(B−1) multiplications, versus B inversions individually.

This makes affine coordinates viable throughout the walk: instead of maintaining Jacobian
coordinates and converting to affine only for partitioning, keep the current point in affine form
at all times. The affine addition formula is:

```
λ = (y₂ − y₁) / (x₂ − x₁)
x₃ = λ² − x₁ − x₂
y₃ = λ·(x₁ − x₃) − y₁
```

The denominator (x₂ − x₁) is the quantity that requires inversion. Collecting B such denominators
and batch-inverting them reduces the inversion count from B to 1 per round of B steps.

**Code reference:** `src/util/batch_inv.rs::batch_invert`, `src/ecdlp/walk.rs::BatchedWalker`,
`src/ecdlp/walk.rs::BatchedWalker::step_all`, `src/ecdlp/mod.rs::solve_dp_batch`.

**Speedup:** 2–4× reduction in wall time per step, depending on batch size B and the relative cost
of inversion vs multiplication on the target hardware. Typical batch sizes are B = 8–32; beyond
B ≈ 32 the marginal gain diminishes as the inversion cost becomes negligible relative to the
3(B−1) multiplications.

---

### GLV endomorphism and the 6-orbit canonical map

**What it adds:** a further ~√3 reduction in expected walk length by collapsing 6-element orbits.

The secp_k1_toy curve has an order-3 endomorphism φ: for any point P = (x, y), φ(P) = (β·x mod p,
y), where β is a cube root of unity mod p (β³ ≡ 1 mod p, β ≠ 1). This endomorphism satisfies
φ(P) = λ·P for all curve points P, where λ is a root of x² + x + 1 ≡ 0 mod n (the minimal
polynomial of a primitive cube root of unity in ℤ/nℤ). Because φ has order 3, the orbit of any
generic point P under the group generated by φ and negation is:

```
{P,  φ(P),  φ²(P),  −P,  −φ(P),  −φ²(P)}
```

This orbit has size 6. Collapsing it to a single canonical representative (the orbit member with
the smallest x-coordinate) reduces the effective group size by a factor of 6, dropping the expected
walk length by √6 vs the plain walk, or √3 vs the negation-map-only walk.

The scalar bookkeeping is more involved than for the negation map alone. If the current walk state
satisfies W = a·G + b·Q, then:

| Orbit member | Point   | Adjusted a    | Adjusted b    |
|--------------|---------|---------------|---------------|
| W            | W       | a             | b             |
| φ(W)         | φ(W)    | λa mod n      | λb mod n      |
| φ²(W)        | φ²(W)   | λ²a mod n     | λ²b mod n     |
| −W           | −W      | n − a         | n − b         |
| −φ(W)        | −φ(W)   | n − λa mod n  | n − λb mod n  |
| −φ²(W)       | −φ²(W)  | n − λ²a mod n | n − λ²b mod n |

After each walk step, all six orbit members are computed (the three x-coordinates require two
multiplications by β; the y-coordinates are shared with their negations), the minimum-x member is
selected, and the scalars are adjusted accordingly. The batched-inversion infrastructure from
the batched-inversion stage is reused: each thread still owns a `BatchedWalker` with B walk states.

Fruitless cycles can still occur in the GLV walk (the canonical walk can oscillate between two
orbit members that share the same canonical representative). The BKNS doubling escape handles this
identically to the negation-map walk.

**Code reference:** `src/ecdlp/glv.rs::glv_phi`, `src/ecdlp/glv.rs::glv_canonical`,
`src/curve/secp_k1_toy.rs` (curve constants β, λ), `src/ecdlp/mod.rs::solve_dp_glv`.

**Speedup:** ~√3 ≈ 1.73× reduction in expected walk steps vs the negation map alone, or ~√6 ≈
2.45× vs distinguished points (plain DP walk). The GLV endomorphism is essentially free to evaluate
— one multiplication by β per orbit member — so the theoretical speedup should be nearly fully
realised in practice. These are expected ratios, not tested ones: the only comparison in the test
suite is `glv_fewer_steps_than_negmap` (unit, `src/ecdlp/mod.rs`, on the 1093-element `tiny_glv`
group), which asserts only that the GLV step total over 20 seeds is at most twice the negation-map
total. A ratio of about 0.62–0.77 was observed when that test was written; the observation is not
asserted, and the test calls the private `solve_dp_glv_impl` because the public `solve_dp_glv` is
`Unsupported` off the `secp_k1_toy` group.

---

### Binary-field and Koblitz variants

**What it adds:** the same rho walk over binary curves y² + xy = x³ + ax² + b on GF(2^m)
(`rho::binary_ecdlp::{solve_brent, solve}`), and a τ-orbit variant for Koblitz curves
(`rho::binary_ecdlp::koblitz::solve_koblitz`).

The binary solvers use the same r-adding walk (16 addends) and Brent detection with the invariant
W = a·G + b·Q, over `BinaryCurve` and `BinaryAffinePoint<F: F2m<1>>` with a `u64` group order, and
return the same `Result<DlogSolution<RhoEvidence>, DlogError>`. The differences are the ones that
matter for the outcomes:

- **Composite orders are allowed.** `check_binary_order_input` checks G ≠ ∞, both points on the
  curve, `n ≠ 0`, `n·G = ∞` and `n·Q = ∞`, and never requires `n` prime (it has no `OrderNotPrime`
  case): the toy curve has order 4. The modular inverse is by extended gcd, so a collision with
  gcd(b2 − b1, n) > 1 is degenerate; it discards the attempt and consumes one of the `max_retries`
  attempts, exactly as for the prime-field solver. With gcd = 1 and `n·G = n·Q = ∞` the unique
  solution satisfies `k·G = Q`, so a failed verification there is a genuine fault, not bad luck (an
  argument from the module documentation; no test reaches that branch). Tests, all in
  `rho/tests/binary_ecdlp_kat.rs`: `input_check_rejects_defective_instances`,
  `input_check_accepts_composite_orders`, `input_check_rejects_target_outside_subgroup`,
  `solvers_reject_invalid_input`, `degenerate_collision_consumes_an_attempt`,
  `zero_retries_exhaust_the_budget_immediately`.
- **The τ-orbit walk folds the orbit map into the scalars.** The Frobenius τ(x, y) = (x², y²) is an
  automorphism of the group. `solve_koblitz` replaces each walk state by the canonical member of its
  τ-orbit and negation, and `tau_canonical` adjusts (a, b) by the eigenvalue λ_j with
  τ^j(G) = λ_j·G (and by the sign) so that the invariant still holds. A collision is therefore
  between identical canonical points: the evidence has `negated: false`, and there is no separate
  τ-power field. A wrong eigenvalue corrupts the scalars and surfaces as `Inconsistent`, never as a
  wrong `Ok` (`corrupted_eigenvalue_is_inconsistent_not_wrong`).
- **Eigenvalues come from enumerating the group.** `compute_tau_eigenvalues` gives [1, 3, 1, 3] on
  the toy curve, where τ(G) = −G (`tau_eigenvalues_satisfy_tau_p_equals_lambda_p`). It enumerates
  all `n` multiples of G, so it reports `Unsupported`, not a panic, above `TAU_EIGENVALUE_MAX_ORDER`
  (2²⁰) and when some τ^j(G) lies outside ⟨G⟩ (`tau_eigenvalues_refuse_huge_orders` and
  `tau_eigenvalues_outside_subgroup_are_unsupported`, unit tests in
  `rho/src/binary_ecdlp/koblitz.rs`); a table shorter than the field degree is `Unsupported` too
  (`short_eigenvalue_table_is_unsupported`).

**What the tests do not show.** On the toy curve every orbit of τ together with negation has a
single x-coordinate, so the canonical tie-break always keeps the k = 0 member and the walk consults
only λ₀; the corrupted-eigenvalue test therefore corrupts λ₀. τ-powers above 0 are not exercised
end to end, and the expected √(2m) saving of the orbit collapse is not measured. A target outside
⟨G⟩ with `n·Q = ∞` is not tested either (there is no non-cyclic binary fixture): the module
documentation says the baseline walk then ends as budget exhaustion and the τ fold can end as
`Inconsistent`, and nothing here claims more.

---

## 3. Cumulative Speedup Table

The table below summarises the expected speedup of each optimization relative to the previous one
and relative to the single-threaded Brent baseline. Speedups marked with √ are reductions in
expected walk length (and therefore in expected time, assuming constant per-step cost); the
batched-inversion speedup is a reduction in per-step cost rather than walk length.

| Optimization                     | Speedup vs prev            | Cumulative vs baseline |
|----------------------------------|----------------------------|------------------------|
| r-adding walk + Brent (baseline) | —                          | 1×                     |
| Distinguished points + N walkers | N× (linear parallelism)    | N×                     |
| Negation map + BKNS escape       | ~√2 ≈ 1.41×                | ~N√2                   |
| Batched inversion + affine walks | ~2–4× (per-step cost)      | ~N√2 · B-factor        |
| GLV 6-orbit canonical map        | ~√3 vs negation map alone  | ~N√6 · B-factor        |

Notes:

- The distinguished-points speedup is N (linear in thread count) because each additional walker
  independently contributes DPs to the collision pool. This is not a constant-factor improvement —
  it scales with hardware.
- The negation map and the GLV endomorphism give √2 and √3 reductions in expected walk *length*.
  These multiply together because they reduce the effective group size independently: the negation
  map halves it (factor of 2), the GLV endomorphism reduces it by 3 (factor of 3), for a combined
  factor of 6 and a combined walk-length reduction of √6.
- Batched inversion reduces per-step *cost* rather than walk length. Its speedup factor (the
  "B-factor") depends on the ratio of inversion cost to multiplication cost on the target hardware,
  and on the batch size B. It composes multiplicatively with the walk-length reductions from the
  negation map and the GLV endomorphism.

---

## 4. How Speedups Compose

The speedups from distinguished points, the negation map, batched inversion, and the GLV
endomorphism multiply rather than add because they act on independent
axes of the total computation cost.

The total expected work is (expected walk steps) × (cost per step). Distinguished points reduce
expected walk steps by a factor of N (more walkers means the birthday collision happens sooner in
wall time). The negation map and the GLV endomorphism each reduce expected walk steps by reducing
the effective group size: the negation map by a factor of 2, the GLV endomorphism by a factor of 3
(GLV orbit). Because these reductions apply to the same quantity — expected walk steps — they
multiply: the combined reduction is √(2 · 3) = √6. Batched inversion reduces the cost per step by
amortising inversions; this multiplies with the walk-step reductions because it is a separate
factor.

The practical ceiling is set by the birthday bound itself: no matter how many optimizations are
applied, the expected number of *collisions* needed is O(1), and each collision requires O(√(n/m))
walk steps where m is the effective group size after all orbit-collapsing maps. For the GLV curve
with negation map, m = n/6, so the expected steps per walker is O(√(n/6)). With N walkers, the
expected wall-time steps is O(√(n/6) / N). Beyond this, further speedup requires either a larger
orbit-collapsing map (not available for this curve family) or more parallelism.

---

## 5. Integer Factorization Rho

The factoring side of the codebase implements three variants of Pollard rho for
integers. The pseudorandom map is f(x) = (x² + c) mod N; a non-trivial factor is detected when
gcd(|xᵢ − xⱼ|, N) ∈ (1, N).

Floyd's algorithm (the baseline) runs two pointers at speeds 1 and 2 and tests a GCD at every
step. Brent's algorithm freezes the slow pointer at the start of each power-of-2 window, reducing
f-evaluations by ~24%. The batched-GCD variant (Montgomery 1987) accumulates the product of B
differences before taking a single GCD, replacing B expensive GCD calls with B−1 cheap
multiplications and one GCD. Multi-c parallelism runs independent instances with different c values
via `rayon::par_iter`, with the first successful thread winning.

The factoring and ECDLP sides share only the high-level structure (pseudorandom walk, cycle
detection, birthday-bound analysis). The specific optimizations diverge: Brent's cycle detection
is used in both, but distinguished points are specific to ECDLP (factoring uses Brent's algorithm
directly), and the negation map and GLV endomorphism have no factoring analogue.

**What `factor` returns.** The entry point is

```rust
pub fn factor(n: u128, max_c: u64, batch_size: usize) -> Result<Factorization, FactorError>;
// Factorization { d, evidence } with 1 < d < n and n % d == 0, checked before it is returned
// FactorEvidence::{ Even, Rho { c } };  FactorError::{ InvalidInput, Prime, Unsupported,
//                                                     BudgetExhausted { constants } }
```

`factor` classifies `n` first, in this order: `n ≤ 1` is `InvalidInput`
(`factor_rejects_zero_and_one`); `n ≥ 2^80`, beyond the exact range of `mulmod`, is `Unsupported`
(`factor_rejects_sizes_beyond_mulmod`); `n = 2` is `Prime`; an even `n > 2` yields `d = 2` with
`Even` evidence and no walk (`factor_even_composites_carry_even_evidence`); an odd `n` is tested
for primality and a prime is reported as `Prime`, never as a factor of itself
(`factor_reports_primes_and_never_returns_two_as_a_factor_of_two`). The primality verdict is a
proof on the whole supported range: `shared_numth::is_prime` is proven below ψ₁₃ ≈ 2^81.5
(`MILLER_RABIN_PROOF_BOUND`), 2^80 < ψ₁₃, and a compile-time assertion keeps `MAX_SUPPORTED` at
most that bound, so even 2^79 − 67 and 2^80 − 65 are `Prime`
(`factor_certifies_primes_up_to_the_supported_limit`). The composite ψ₁₂ ≈ 2^78.1 passes the first
twelve Miller–Rabin bases but not the thirteenth; it is split by rho, never certified
(`factor_splits_psi_12_instead_of_certifying_it`). Only an odd composite reaches the walks.

Each walk is bounded by `1024 + 64·(⌊n^{1/4}⌋ + 1)` evaluations of f (read from the code,
`attempt_steps`; no test pins the formula). `max_c` is the number of constants c = 1, …, `max_c`
tried; c ≡ 0 and c ≡ n − 2 (mod n) are skipped and not counted (by code reading). A spent budget is
`BudgetExhausted { constants }`, where `constants` counts the walks that ran and failed, and proves
nothing about `n`: it is not a primality verdict (`factor_with_no_constants_is_budget_exhausted`;
`factor_counts_the_constants_whose_walks_failed`, where n = 73·163 exhausts four constants and a
budget of sixteen succeeds). The primitives `floyd`, `brent` and `brent_batched` take an explicit
`max_steps` and return an `Option` with one meaning, *this attempt failed*
(`primitives_fail_when_the_step_bound_is_spent`); a batch size of zero is treated as one
(`brent_batched_zero_batch_size_terminates`). Squares and prime powers need no special case: the
returned divisor is checked, and `factor_squares_and_prime_powers_give_proper_divisors` covers them.
All of these tests are in `rho/tests/factor_kat.rs`. The `rho-factor` CLI reports each class
separately (README, "CLIs").

**Code reference:** `src/factor/rho.rs` — `floyd`, `brent`, `brent_batched`, `factor`.

---

## 6. Running the Experiments Yourself

### Prerequisites

```
cargo build --release
```

### Known-answer tests

```
cargo test
```

This runs all KATs: field arithmetic (both `FpNaive` and `FpMonty` against reference values),
curve arithmetic (scalar multiplication cross-checked against reference points), factoring (known
semiprimes), and ECDLP (known discrete logarithms on the toy curves). All tests must pass before
trusting benchmark results.

### ECDLP benchmark suite

```
cargo bench --bench ecdlp
```

This benches each solver variant: `solve_brent` (single-threaded baseline), `solve_dp`
(distinguished points), `solve_dp_negmap` (negation map), `solve_dp_batch` (batched inversion), and
`solve_dp_glv` (GLV endomorphism), all on `secp_k1_toy` with `Q = 12_345_678_901·G`. That group has
a 63-bit prime order, and a collision in it takes about 1.25·√n ≈ 2.7·10⁹ walk steps (arithmetic
from n, not a measurement), so these benches are not practically runnable: the cumulative speedup
is demonstrated at test scale by the 20-bit KATs, not by this suite. Each bench body unwraps its
solver, so a budget-exhausted or otherwise failed solve aborts the bench rather than being
timed. The distinguished-point benches pass a step budget of 2⁴⁰ ≈ 1.1·10¹², a safety cap a few
hundred times the expected cost; the Brent bench passes 20 attempts, each capped at 2²⁹ − 1 ≈
5.4·10⁸ steps, below the expected cost of one collision, so by that arithmetic it would most
likely end as budget exhaustion and abort. These figures are arithmetic on the code and the group
order; the benches were not run to produce them.

### CLI: stepping through the optimization layers

The `rho-dlog` binary exposes the solver variants via command-line flags. Every invocation needs
the target `--q <x,y>`. The following sequence demonstrates the progression from the single-threaded
baseline to the fully-optimized solver:

```sh
# Single-threaded Brent baseline (--walkers 0, the default)
rho-dlog --curve secp-toy --q <x,y>

# 4 parallel walkers with distinguished points
rho-dlog --curve secp-toy --q <x,y> --walkers 4

# 4 walkers with batched inversion (batch size 16)
rho-dlog --curve secp-toy --q <x,y> --walkers 4 --batch-size 16

# 4 walkers, batch size 16, GLV endomorphism (secp-toy only)
rho-dlog --curve secp-toy --q <x,y> --walkers 4 --batch-size 16 --glv
```

`--walkers N` with `N ≥ 1` selects the distinguished-point solvers; `--walkers 0` (the default)
selects `solve_brent`, so `--walkers 1` is one distinguished-point walker, not the Brent baseline.
The `--theta` flag controls the distinguished-point threshold (default: 10). Increasing θ reduces
memory usage at the cost of more steps per walker; decreasing it increases memory usage but
detects collisions sooner. `--max-steps S` caps the total walk steps of the distinguished-point
solvers (default: 64·√n plus a start-up allowance of 64·walkers·batch·2^θ); reaching it prints a
budget-exhausted message, which says nothing about whether a logarithm exists.

```sh
# Explicit theta: 1-in-256 DPs (theta=8)
rho-dlog --curve secp-toy --q <x,y> --walkers 4 --batch-size 16 --glv --theta 8
```

On success the CLI prints `k` and a verification line; the solver has already checked
`k·G == Q`, so the CLI does not recheck it. Each `DlogError` class gets its own message and exits
with status 1: invalid input (naming the defect), unsupported domain, budget exhaustion,
inconsistent result, and solver invariant failure (`describe_error` in `src/ecdlp/cli.rs`). The
`generic` curve's built-in order, 9 223 372 037 218 517 353, is composite
(7 · 17 · 77 507 328 043 853 087; the product can be checked by hand), so `--curve generic`
reports `invalid input: group order is not prime` for the Brent and the distinguished-point
solvers alike: they all check the input first, and
`solve_brent_rejects_composite_order_before_walking` (unit, `src/ecdlp/mod.rs`) tests that outcome
for a composite `n`. No test runs the CLI itself, and it does not project onto a prime-order
subgroup.

---

## 7. Further Reading

1. **Pollard, J. M. (1978).** "Monte Carlo methods for index computation (mod p)." *Mathematics of
   Computation*, 32(143), 918–924. The original rho algorithm for discrete logarithms.

2. **Teske, E. (1998).** "Speeding up Pollard's rho method for computing discrete logarithms."
   *Algorithmic Number Theory Symposium (ANTS-III)*, LNCS 1423, 541–554. Introduces the r-adding
   walk and proves it mixes better than the naive random walk.

3. **van Oorschot, P. C., and Wiener, M. J. (1999).** "Parallel collision search with
   cryptanalytic applications." *Journal of Cryptology*, 12(1), 1–28. The definitive reference for
   distinguished points and the parallel rho architecture.

4. **Bos, J. W., Kleinjung, T., Niederhagen, R., and Schwabe, P. (2012).** "ECC2K-130 on Cell
   CPUs." *Progress in Cryptology — AFRICACRYPT 2010*, LNCS 6055, 225–242. Introduces the
   deterministic fruitless-cycle escape (the BKNS doubling perturbation) used in the negation-map
   walk.

5. **Gallant, R. P., Lambert, R. J., and Vanstone, S. A. (2001).** "Faster point multiplication on
   elliptic curves with efficient endomorphisms." *Advances in Cryptology — CRYPTO 2001*, LNCS
   2139, 190–200. The GLV decomposition; the endomorphism used in the GLV walk is the order-3
   special case.

6. **Montgomery, P. L. (1987).** "Speeding the Pollard and elliptic curve methods of
   factorization." *Mathematics of Computation*, 48(177), 243–264. Introduces both the batched-GCD
   trick (factoring rho) and the batched-inversion trick (batched-inversion walk).

7. **Bernstein, D. J., and Lange, T. (2007).** "Faster addition and doubling on elliptic curves."
   *Advances in Cryptology — ASIACRYPT 2007*, LNCS 4833, 29–50. A comprehensive reference for
   explicit addition formulas in various coordinate systems, including the affine formulas used in
   the batched-inversion walk.

8. **Brent, R. P. (1980).** "An improved Monte Carlo factorization algorithm." *BIT*, 20(2),
   176–184. Brent's cycle detection algorithm, used in both the factoring and ECDLP
   solvers.

---

# Algebraic ECDLP Attacks: An Integrative Chapter

This chapter is the integrative code-tour for the complete algebraic ECDLP attack survey. It is
the depth chapter for the algebraic attacks in the umbrella synthesis (§0 above). Each attack or
primitive is implemented in isolation in the `rho` crate. This chapter surveys the whole attack
landscape in a single narrative, names the interfaces that connect the attacks, states what each
attack implements and claims, and annotates the scale-dependent phenomena at demonstration
fidelity.

The algebraic attacks form an *attack survey*, not a linear pipeline. The attacks do not compose
sequentially —
each exploits a different curve structure and applies only on the curve whose precondition it
requires. The "at a glance" view is therefore a **taxonomy table** (which structure unlocks which
escape), not a data-flow diagram.

For the mathematics — the pairing reduction, the p-adic logarithm, the Weil descent, and the
per-attack L-notation complexity — see `docs/MATHEMATICS.md` ch. 10. This
chapter is code-first: it assumes the reader knows the mathematics and shows how the code
realises it.

---

## 8. Algebraic ECDLP Attacks: At a Glance

The through-line of the algebraic attacks is the §"Escape from Search" framing from
`docs/MATHEMATICS.md`: every attack finds a curve structure that escapes the generic $\sqrt{n}$
bound. The table below names the structure, the module that realises it, and the toy fixture.

| Attack | Curve structure exploited | Module surface | Toy fixture |
|--------|--------------------------|----------------|-------------|
| Pollard rho (baseline) | None — generic group | `rho::ecdlp` (`solve_brent`, `solve_dp`, `solve_dp_negmap`, `solve_dp_batch`, `solve_dp_glv`) | `tiny_a`, `tiny_b` (20-bit prime order) in the KATs; `secp_k1_toy` (63-bit, prime order) in the benches and CLI |
| Pohlig–Hellman | Composite group order $n = \prod p_i^{e_i}$ | `rho::ecdlp::pohlig` (`solve_ecdlp_composite`) | `composite_toy()` ($n = 60 = 2^2 \cdot 3 \cdot 5$); `composite_rho_toy()` ($n = 1068 = 2^2 \cdot 3 \cdot 89$) |
| MOV/Frey–Rück | Small embedding degree $k$; $\ell \mid p^k - 1$ | `rho::pairing::mov` (`mov_reduce`) | `pairing_toy()` ($\ell = 3$, $k = 2$, $\mathbb{F}_{47^2}$) |
| Smart–Satoh–Araki (SSA) | Anomalous curve: $\#E(\mathbb{F}_p) = p$ (trace = 1) | `rho::ssa` (`ssa_solve`) | `anomalous_toy()` ($y^2 = x^3 + 5 \bmod 7$, $\#E = 7 = p$) |
| GHS/Weil descent (**transfer**) | Binary curve $E/\mathbb{F}_{2^m}$ with subfield tower | `rho::ghs` (`ghs_descend`, `verify_log_preservation`, `transfer_point`) | `ghs_toy_curve()` ($\mathbb{F}_{2^6}$, $m = 6$, $l = 2$) |
| Semaev / index calculus | Semaev-decomposable points over a factor base | `rho::semaev` + `rho::index_calculus` (`semaev_poly`, `index_calculus_dlp`, `collect_relations`, `decompose`) | `IndexCalcStrategy::toy()` ($\ell = 5$, $|FB| = 6$, $m = 2$) |

**The through-line in one sentence per row.** Pollard rho is the generic $\sqrt{n}$ baseline —
no structure exploited. Pohlig–Hellman finds a *group-order homomorphism*: composite order
factors the ECDLP over prime-order subgroups via CRT. MOV finds a *pairing*: the bilinear map
$e: E[n] \times E[n] \to \mu_n$ transports the ECDLP to $\mathbb{F}_{p^k}^*$, where index
calculus applies. SSA finds an *anomalous endomorphism*: the p-adic formal group logarithm
gives a polynomial-time solve. GHS finds a *field tower*: Weil restriction transfers the ECDLP
to a hyperelliptic Jacobian DLP over $\mathbb{F}_{2^l}$ (a transfer, not an end-to-end solve).
Index calculus finds *factor-base decomposability*: Semaev summation polynomials build a
relation matrix whose left-kernel combination over $\mathbb{F}_\ell$ recovers the discrete log
modulo a prime $\ell \mid n$.

---

## 9. Pohlig–Hellman: CRT to Prime-Order Subgroups

### What it exploits

The group order $n = \#E(\mathbb{F}_p)$ is composite: $n = \prod p_i^{e_i}$. The ECDLP in a
group of composite order reduces to independent DLPs in each prime-power subgroup via the
Chinese Remainder Theorem. A prime-power subgroup is solved digit by digit in base $p_i$, each
digit being a DLP in a group of prime order $p_i$. For $p_i \leq 64$ that DLP is solved by
enumerating the subgroup, which is the digit step of Pohlig and Hellman (1978, §IV); for larger
$p_i$ it is solved by Pollard rho (`solve_brent`), which requires prime order. Pohlig and Hellman's
own large-prime step is a time–memory search, not rho, so the paper is cited only for the
small-subgroup step. The CRT reassembles the answer.

The escape: instead of running rho on the full group of order $n$, solve each prime-order digit
problem on its own. The cost drops from $O(\sqrt{n})$ to $O(\sum e_i (\log n + \sqrt{p_i}))$ when
every digit is a rho walk — a dramatic speedup when $n$ has small prime factors.

### The module surface

```rust
// rho::ecdlp::pohlig

/// Decompose a u64 group order into prime-power factors.
pub fn factor_order(n: u64) -> Vec<(u64, u32)>;

/// Project (G, Q) to the order-p^e subgroup via [n / p^e]-scalar-multiplication.
pub fn project_to_subgroup<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    p_power: u64,
) -> (AffinePoint<F>, AffinePoint<F>);

/// Composite-order ECDLP via prime-power lift + CRT combine, verified k·G == Q.
pub fn solve_ecdlp_composite<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Result<DlogSolution<PohligEvidence>, DlogError>;

pub struct PohligEvidence { pub components: Vec<PrimePowerComponent> }
pub struct PrimePowerComponent { pub p: u64, pub e: u32, pub residue: u64,
                                 pub digits: Vec<DigitMethod> }
pub enum DigitMethod { SubgroupEnumeration, Rho(RhoEvidence) }
```

`solve_ecdlp_composite` is the entry point. It checks the instance, calls `factor_order` to
decompose $n$, projects to each prime-power subgroup via `project_to_subgroup`, recovers
$k \bmod p^e$ digit by digit (each digit by `SubgroupEnumeration` when $p \leq 64$, otherwise by
`solve_brent` on the order-$p$ subgroup, with its evidence carried unchanged as `Rho`), reassembles
via CRT, and checks `k·G == Q` before returning. Every digit solve receives the prime order $p$,
never $p^e$.

| Outcome | Meaning | Exercised by |
|---------|---------|--------------|
| `Ok(DlogSolution { k, evidence })` | verified `k`; each component's `residue` is `k mod p^e`, the sum of its digits, with one `DigitMethod` per digit, `SubgroupEnumeration` only for `p ≤ 64` and `Rho` only above | `composite_dlog_every_k` (all `k` in `0..60`), `composite_rho_dlog_sampled_k` (via `check_pohlig_evidence`, `rho/tests/ecdlp_kat.rs`) |
| `Err(InvalidInput(ZeroOrder))`, `(WrongGeneratorOrder)` | `n = 0` (no longer a panic); `n·G ≠ ∞`, or `(n/p)·G = ∞` for some prime `p ∣ n` (a multiple of the true order, or a `G` of smaller order) | `composite_wrong_order_is_invalid_input`, `composite_rho_wrong_order_is_invalid_input` |
| `Err(InvalidInput(TargetOutsideSubgroup))` | `n·Q ≠ ∞` (by code reading; no composite test), or an order-`p` enumeration finds no match, which proves `Q ∉ ⟨G⟩` | `composite_target_outside_subgroup_is_invalid_input` (the enumeration case) |
| `Err(BudgetExhausted { .. })`, `Err(Inconsistent { .. })` from a rho digit | passed through unchanged, so exhaustion is never reported as anything else | not exercised: no fixture's rho digit exhausts or fails |
| `Err(Inconsistent { k })` after the CRT step | every digit solved but `k·G ≠ Q`; judged unreachable on a verified instance | not exercised |

`composite_toy()` has only primes up to 5, so no digit reaches rho there; the rho path is exercised
on `composite_rho_toy()` (below). The `DigitMethod` name follows the source: the order-`p` step in
Pohlig and Hellman (1978, §IV) is a function with `p` possible values computed by table lookup or
by trial and error, and the code takes the trial-and-error form, so it stores no table.

### The toy KAT/fixture

The fixture is `composite_toy()`: $y^2 = x^3 + x + 33 \bmod 47$, $n = 60 = 2^2 \cdot 3 \cdot 5$.
The generator $G = (10, 3)$. The known scalar is $k = 7$: $Q = 7 \cdot G$.

```rust
let curve = composite_toy();
let g = curve.generator::<FpMonty>();
let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
let solution = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
    .expect("Pohlig–Hellman must succeed on composite_toy");
assert_eq!(curve.scalar_mul(&g, &Uint::<4>::from(solution.k)), q);
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_pohlig_hellman`) asserts this before
timing. The Pohlig–Hellman KATs live in `rho/tests/ecdlp_kat.rs` (`composite_dlog_*`,
`composite_rho_*`, `composite_wrong_order_is_invalid_input`,
`composite_target_outside_subgroup_is_invalid_input`).

The second fixture, `composite_rho_toy()`, is $y^2 = x^3 + 2x + 3 \bmod 1009$ with
$G = (10, 425)$ and $n = 1068 = 2^2 \cdot 3 \cdot 89$. The factor $89$ exceeds the enumeration
limit, so its digit is solved by rho. Each rho solve costs about 0.16 s in a debug build, so
`composite_rho_dlog_sampled_k` samples every 23rd $k$ plus the edges and the residues 0 and 1
modulo 89 (whose rho digit is `Trivial`), not all 1068 values.

### Cross-reference

For the mathematical development — the CRT reduction, the prime-power lift, and the L-notation
cost — see `docs/MATHEMATICS.md` §10.1.

---

## 10. MOV/Frey–Rück: The Pairing Reduction

### What it exploits

The curve has small embedding degree $k$: $\ell \mid p^k - 1$ but $\ell \nmid p^j - 1$ for
$j < k$. The Weil/Tate pairing $e: E[\ell] \times E[\ell] \to \mu_\ell \subset \mathbb{F}_{p^k}^*$
is bilinear and non-degenerate. Given $Q = k \cdot G$, compute $e(G, R)$ and $e(Q, R)$ for a
$\mu_\ell$-generator $R$. By bilinearity, $e(Q, R) = e(k \cdot G, R) = e(G, R)^k$. Writing
$g_0 = e(G, R)$ and $h_0 = e(Q, R)$, the ECDLP scalar $k$ is exactly $\log_{g_0}(h_0)$ in
$\mathbb{F}_{p^k}^*$ — where index calculus applies.

This is the cross-crate bridge: the MOV reduction transports the ECDLP to the NFS-DL setting
(the `gnfs` crate's discrete-log solver), which the `gnfs::dl::solve_dl` entry point solves.

### The module surface

```rust
// rho::pairing::mov

/// MOV/Frey–Rück reduction: ECDLP → F_{p^k}* DLP via the Tate pairing.
///
/// Given (curve, G, Q, R, ell) where Q = k·G in the order-ℓ subgroup and R is a
/// μ_ℓ-generator with e(G, R) ≠ 1, returns k mod ℓ with the pairing values behind it,
/// after checking k·G == Q directly in E[ℓ].
pub fn mov_reduce<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    g_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    r_point: &PairingPoint<F>,
    ell: u64,
) -> Result<MovSolution<F>, MovError>;

pub struct MovSolution<F> { pub k_mod_ell: u64, pub ell: u64,
                            pub g_pair: FpExt<F>, pub h_pair: FpExt<F> }  // #[non_exhaustive]
pub enum MovError { InvalidInput(MovInputDefect), DegeneratePairing, Unsupported(&'static str),
                    Dl(gnfs::dl::SolveDlError), Inconsistent { k_mod_ell: u64 } }
```

`mov_reduce` validates the instance, evaluates the reduced Tate pairing twice (via `reduced_tate`
in `rho::pairing::tate`), encodes both outputs through the bridge encoder (`fpext_to_bigint`), and
calls `gnfs::dl::solve_dl` to recover $k \bmod \ell$. It returns the result only after checking
`k·G == Q` on the curve, independently of the pairing and the log engine, and `g_pair^k == h_pair`
holds in the returned solution.

| Outcome | Meaning | Exercised by |
|---------|---------|--------------|
| `Ok(MovSolution { k_mod_ell, ell, g_pair, h_pair })` | verified `k mod ℓ`; `g_pair^k = h_pair`, `g_pair ≠ 1`, `g_pair^ℓ = 1` | `mov_reduce_recovers_k1`, `mov_reduce_recovers_k2`, `mov_reduce_recovers_k0_for_target_at_infinity` (`rho/tests/mov_kat.rs`) |
| `Err(InvalidInput(defect))` | `ℓ` not prime, `ℓ ∤ p^k − 1`, `G = ∞`, or `G`, `Q`, `R` off the curve or outside `E[ℓ]` | `mov_reduce_rejects_composite_ell`, `mov_reduce_rejects_ell_not_dividing_field_order`, `mov_reduce_rejects_infinite_generator`, `mov_reduce_rejects_off_curve_points`, `mov_reduce_rejects_points_outside_e_ell` |
| `Err(DegeneratePairing)` | `R = ∞`, `e(G, R) = 1`, or a Miller evaluation of `0` (as for `R = G`): `R` is not a μ_ℓ-generator for `G` | `mov_reduce_degenerate_r_is_degenerate_pairing`, `mov_reduce_infinite_r_is_degenerate_pairing` |
| `Err(Unsupported(..))` | `p ≥ 2^64` (a non-zero word above the lowest), reported before any arithmetic on `p` | `mov_reduce_rejects_field_above_one_limb` |
| `Err(Dl(..))` | the `solve_dl` error, passed through unchanged | not exercised |
| `Err(Inconsistent { k_mod_ell })` | the engine's value fails `k·G == Q`; an invariant failure | `log_matches_accepts_only_the_true_log` (unit, `src/pairing/mov.rs`; checks the verifier, no fixture reaches it) |

The single-limb limit `p < 2^64` of the bridge is enforced by `mov_reduce` as an error
(`Unsupported`), in release builds too, so no later conversion truncates. The final exponent
`(p^k − 1)/ℓ` is computed exactly (it can exceed one limb even for a one-limb `p`);
`reduced_tate_with_exponent_above_one_limb_lands_in_mu_ell` (`src/pairing/mov.rs`) checks it on
`y² = x³ + 1` over `F_68719476767`.

The pairing infrastructure:
- `rho::pairing::miller` — Miller's algorithm (the efficient pairing evaluation).
- `rho::pairing::tate` — the reduced Tate pairing (Miller loop + final exponentiation).
- `rho::pairing::weil` — the Weil pairing (for reference; the Tate pairing is used in `mov_reduce`).
- `rho::pairing::fpext` — $\mathbb{F}_{p^k}$ arithmetic (the target field).
- `rho::pairing::ecext` — elliptic curve arithmetic over $\mathbb{F}_{p^k}$ (for the extension-field
  torsion points).

### The toy KAT/fixture

The fixture is `pairing_toy()`: the base curve $y^2 = x^3 + x + 33 \bmod 47$, $\ell = 3$,
$k = 2$, $\mathbb{F}_{47^2}$. The known scalar is $k = 2$: $Q' = 2 \cdot G$ in the order-3
subgroup.

```rust
let (curve, modulus, ell, p_point, q_point) = pairing_toy();
let q_prime = p_point.scalar_mul(2, &a_ext, &modulus, &p);
let solution = mov_reduce::<FpNaive4<4>>(&curve, &modulus, &p_point, &q_prime, &q_point, ell)
    .expect("MOV must succeed on pairing_toy");
assert_eq!(solution.k_mod_ell, 2u64); // the log modulo ℓ = 3, verified k·G == Q'
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_mov`) asserts this before timing. The
`2` here is `k_mod_ell`, the logarithm modulo $\ell = 3$.

### Cross-reference

For the full proof of the MOV reduction — the bilinearity, non-degeneracy, and the reduction
argument — see `docs/MATHEMATICS.md` §10.5 (the designated payoff proof).

---

## 11. Smart–Satoh–Araki: The Anomalous-Curve Polynomial-Time Attack

### What it exploits

The curve is *anomalous*: $\#E(\mathbb{F}_p) = p$ (trace of Frobenius = 1). The p-adic formal
group logarithm turns the ECDLP into a division in $\mathbb{F}_p$. The algorithm is polynomial in
$\log p$ — the sharpest escape among the algebraic attacks — although this implementation is not
(see below).

The key steps: (1) lift the curve to $\mathbb{Z}_p$ by shifting $a \mapsto a + p\,s$ with
$s \not\equiv 0 \pmod p$ (the same curve modulo $p$), and lift the affine points
$G, Q \in E(\mathbb{F}_p)$ to points $\tilde{G}, \tilde{Q}$ on it via Hensel's lemma; (2) multiply
each lift by $p$ on that lifted curve: because $\#E(\mathbb{F}_p) = p$, $p \cdot \tilde{G}$ and
$p \cdot \tilde{Q}$ reduce to $\mathcal{O}$, so they lie in the kernel of reduction
$E_1 \cong \hat{E}(p\mathbb{Z}_p)$ (whether they lie in the deeper $E_2$ depends on the lift of the
curve — the logarithm would then carry no information); (3) apply the formal group logarithm
$\psi$ to both, which lands in $p\mathbb{Z}_p$; (4) recover
$k \equiv \psi(p \cdot \tilde{Q}) / \psi(p \cdot \tilde{G}) \pmod{p}$ from the quotient of the
leading terms. Step (4) is justified when $p \cdot \tilde{G} \in E_1 \setminus E_2$ for the lift of
the curve used; the derivation, why the requirement is on the lift of the curve, and the source for
the failing (canonical) lift are in `docs/MATHEMATICS.md` §10.2.

### The module surface

```rust
// rho::ssa

/// Solve the ECDLP on an anomalous curve; the result is checked k·G == Q.
///
/// Requires #E(F_p) = p (anomalous), established by an exact point count, and n == p.
/// Returns k mod p such that Q = k·G; Q = ∞ gives 0.
pub fn ssa_solve<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Result<u64, SsaError>;

pub enum SsaError {
    InvalidInput(SsaInputDefect), NotAnomalous, Unsupported(&'static str),
    Inconsistent { k: u64 }, InvariantFailure(&'static str),
}
pub enum SsaInputDefect { FieldNotPrime, SingularCurve, GeneratorAtInfinity,
                          GeneratorNotOnCurve, TargetNotOnCurve, OrderNotFieldPrime }
```

The internal structure:
- `rho::ssa::lift` — Hensel lift: $F_p$ point → $\mathbb{Z}_p$ point on a lifted curve
  (`lift_point(point, curve, a_lift, k)`; $x$ exactly, $y$ by Hensel's lemma), the coefficient of
  the lifted curve (`lifted_a(curve, shift)` $= a + p \cdot$ `shift`), and the lift check
  (`check_lift_on_curve(x, y, curve, a_lift)`). The helpers report `LiftError`, which is not an
  `SsaError`.
- `rho::ssa::formal_log` — the logarithm applied to $p \cdot \tilde{P}$ (`elliptic_log_proj`).
  Despite the name it evaluates $\log(1 + t)$ at $t = -X/Y$, the logarithm of the multiplicative
  formal group, not the series $\psi$ of $E$'s formal group; the two agree modulo $p^2$ on $E_1$.
  `formal_parameter_valuation` reads the filtration step $v_p(t)$.
- `rho::ssa::reduce` — the full SSA reduction assembling the steps (`ssa_solve`): one lifted curve
  per attempt, passed to both the point lifts and the multiplication by $p$
  (`multiply_by_p_proj`, $p - 1$ repeated projective additions), logarithm, division in
  $\mathbb{F}_p$ (`fp_divide_logs`), then verification (`verify_candidate`).
  `kernel_filtration_level(curve, point, a_lift)` exposes the filtration step of $p \cdot \tilde{P}$
  on one lift, for tests.

The p-adic substrate (arithmetic, Hensel lifting, and the p-adic logarithm) is provided by
`shared-padic`.

**Cost.** The multiply-by-$p$ loop and the exact anomalous check (`verify_anomalous`, a Legendre
sum over all $x$) are both $O(p)$, so this implementation is toy-scale; it does not realise the
polynomial-time bound, which needs double-and-add for the multiplication by $p$ (the loop here is
$O(p)$, double-and-add is $O(\log p)$) and a point-counting method better than summation. The
precision is fixed at $p^4$ (`WORKING_PRECISION`).

**Outcomes.** The quotient is the answer; `ssa_solve` still returns a `k` only after checking
`k·G == Q`.

| Outcome | Meaning | Exercised by |
|---------|---------|--------------|
| `Ok(k)` | verified `k` with `k·G == Q`; `Q = ∞` gives `0` immediately, on any curve | `kat_ssa_solve_recovers_every_k` (every `k` in `1..7`), `kat_second_fixture_recovers_every_k` (every `k` in `1..11`, `a ≠ 0`), `kat_ssa_solve_infinity_is_zero` (`rho/tests/ssa_kat.rs`) |
| `Err(NotAnomalous)` | an exact count gives `#E ≠ p`; a proven negative | `kat_ssa_solve_not_anomalous` |
| `Err(Unsupported(..))` | `p`, `a` or `b` does not fit one 64-bit limb, or `p < 5`; never reported as non-anomalous | `kat_large_prime_is_unsupported_not_non_anomalous` |
| `Err(InvalidInput(..))` | composite modulus, singular curve, `G = ∞`, `G` or `Q` off the curve, `n ≠ p` | `kat_composite_modulus_is_field_not_prime`, `kat_singular_curve_is_rejected_before_the_count`, `kat_generator_at_infinity`, `kat_generator_off_curve`, `kat_target_off_curve`, `kat_order_must_equal_field_prime` |
| `Err(Inconsistent { k })` | the quotient of logarithms, the only candidate, fails `k·G == Q`; an invariant failure, never retried | `verify_candidate_rejects_unverifiable_quotients` (unit, `src/ssa/reduce.rs`); no fixture reaches it through `ssa_solve` |
| `Err(InvariantFailure(..))` | a lift, p-adic or logarithm breach on validated input, or `p·G̃ ∈ E₂` on both shifted lifts | `two_degenerate_lifts_are_an_invariant_failure` and the `formal_log` / `fp_divide_logs` unit tests |

**A degenerate lift is detected and replaced, not compensated for.** On the toy fixture the stored
lift ($a = 0$) has $p \cdot \tilde{G} \in E_2$, with parameter valuation $2$
(`kat_stored_lift_of_toy_curve_is_degenerate`), and the quotient of leading units there is not
$\lambda$. `ssa_solve` lifts the curve with $a \mapsto a + p\,s$ instead; on every shift
$s \not\equiv 0$ the raw quotient is exactly $k$
(`attempt_on_shifted_lift_returns_the_raw_quotient`), and there is no halved candidate. If the
first shift is the degenerate one, one second shift is tried
(`ssa_solve_retries_when_the_first_shift_is_degenerate`). That at most one class of shifts is
degenerate is measured on the fixtures, for every point, by
`kat_at_most_one_shift_class_is_degenerate`; see `docs/MATHEMATICS.md` §10.2 for the evidence and
its limits.

### The toy KAT/fixture

The fixture is `anomalous_toy()`: $y^2 = x^3 + 5 \bmod 7$, $p = 7$, $\#E = 7 = p$. The
generator $G = (3, 2)$. The known scalar is $k = 3$: $Q = 3 \cdot G$. A second fixture with
$a \neq 0$, $y^2 = x^3 + x + 5$ over $\mathbb{F}_{11}$ with $G = (0, 4)$, comes from a recorded
bounded enumeration in `rho/tests/ssa_kat.rs` (`kat_second_fixture_matches_its_recorded_search`),
and `verify_anomalous` re-proves it anomalous.

```rust
let curve = anomalous_toy();
let g = curve.generator::<FpNaive>();
let q = curve.scalar_mul(&g, &Uint::<4>::from(3u64));
let k = ssa_solve::<FpNaive>(&curve, &g, &q, ANOMALOUS_TOY_P)
    .expect("SSA must succeed on anomalous_toy");
assert_eq!(k, 3u64);
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_ssa`) asserts this before timing.

### Cross-reference

For the mathematical development — the formal group, the lift-independent homomorphism that makes
the quotient meaningful, and the complexity of the algorithm — see `docs/MATHEMATICS.md` §10.2.

---

## 12. GHS/Weil Descent: The Binary-Curve Transfer

### What it exploits

The curve is defined over a binary field $E/\mathbb{F}_{2^m}$ with a subfield tower
$\mathbb{F}_{2^l} \subset \mathbb{F}_{2^m}$ (with $l \mid m$). The Weil restriction
$\mathrm{Res}_{\mathbb{F}_{2^m}/\mathbb{F}_{2^l}}$ transfers the ECDLP on $E$ to a DLP on the
Jacobian of a hyperelliptic curve $C/\mathbb{F}_{2^l}$, where index calculus applies.

**GHS is a transfer, not an end-to-end solve.** `rho::ghs` exposes the descent reduction and
log-preservation verification, but has no `ghs_dlp` solver. The downstream solve is index
calculus on the Jacobian, which is not implemented. Representing GHS as an end-to-end solve time
would misrepresent the attack's scope.

### The module surface

```rust
// rho::ghs

/// Top-level GHS descent: (E, g, h) → (C, D_g, D_h) with log-preservation guarantee.
pub fn ghs_descend(
    params: &GhsParams,
    g: &BinaryAffinePoint<F2mNaive<1>>,
    h: &BinaryAffinePoint<F2mNaive<1>>,
) -> Result<GhsDescentResult, GhsError>;

/// Verify that the descent preserves the discrete logarithm.
/// Returns true iff k·D_g = D_h on Jac(C) for the known scalar k.
pub fn verify_log_preservation(result: &GhsDescentResult, k: u64) -> bool;

/// Transfer a single point on E to a divisor on Jac(C).
pub fn transfer_point(
    p: &BinaryAffinePoint<F2mNaive<1>>,
    curve_c: &HyperellipticCurve<1>,
    params: &GhsParams,
) -> Result<MumfordDivisor<F2mNaive<1>, 1>, GhsError>;

pub enum GhsError {
    SubfieldDivisibility, NonDescendable, NotOnCurve, NotInSubfield, SubfieldPolyMismatch,
}
```

`transfer_point` guards its input. A point that does not satisfy the curve equation of $E$ is
`NotOnCurve` (`transfer_rejects_off_curve_point`; `ghs_descend_rejects_off_curve_point` for the
base point or the target), because the conorm map is a homomorphism only on points of $E$. The
point at infinity is not an error: it maps to the identity divisor `[1, 0]`
(`transfer_infinity_is_zero_divisor`). A conorm image with a coordinate not fixed by the relative
Frobenius $c \mapsto c^{2^l}$ is `NotInSubfield`: the image is a point of $C(\mathbb{F}_{2^l})$ only
when $E$ is defined over $\mathbb{F}_{2^l}$ (`transfer_rejects_conorm_outside_subfield`, with
coefficient $a = \alpha \notin \mathbb{F}_{2^2}$). All of these tests are in
`rho/tests/ghs_kat.rs`.

**Identifying the subfield.** $\mathbb{F}_{2^l}$ sits inside $\mathbb{F}_{2^m}$ as
$\mathbb{F}_2(\beta)$ with $\beta = \alpha^{(2^m-1)/(2^l-1)}$, and the bits of an
$\mathbb{F}_{2^l}$ element are its coordinates in the basis $\{1, \beta, \dots, \beta^{l-1}\}$.
Both `extract_ghs_curve` (the coefficients of $F$) and `transfer_point` (the conorm image) read
elements through `shared_gf2m::restrict`, so the curve and the divisors live in the same field.
Reducing the raw bit pattern modulo `poly_l` would not do: it is not a field homomorphism
($x^2 + x + 1$ does not divide $x^6 + x + 1$) and on $\mathbb{F}_4 \subset \mathbb{F}_{2^6}$ it
exchanges the conjugates $\beta$ and $\beta^2$. The reading is an isomorphism only when `poly_l` is
the minimal polynomial of $\beta$, so both functions first reject, before touching any coefficient
or point, an $l$ that does not divide $m$ (`SubfieldDivisibility`) and a `poly_l` that is not that
minimal polynomial (`SubfieldPolyMismatch`; `mismatched_poly_l_is_rejected`). The fixture's two
transfer known answers, $D_P = [X+1, 2]$ and $D_{G+P} = [X+1, 3]$, are derived by hand in the test
doc comments and cross-checked with PARI/GP; `docs/MATHEMATICS.md` §10.3 gives the worked
derivation. The homomorphism and log-preservation tests
(`transfer_homomorphism_*`, `ghs_descend_log_preservation_*`) hold under either reading because the
fixture curve has coefficients in $\mathbb{F}_2$;
`curve_and_transfer_share_the_subfield_identification` is the test that exercises the
identification, on a curve with $a = \beta$.

The internal structure:
- `rho::ghs::descent` — Artin–Schreier extension and Weil restriction of scalars
  (`ArtinSchreierData`, `WeilRestriction`, `weil_restrict_poly`).
- `rho::ghs::curve` — hyperelliptic curve extraction from the descent algebra
  (`extract_ghs_curve`, `ghs_genus`).
- `rho::ghs::reduce` — the GHS reduction assembling descent + transfer (`ghs_descend`,
  `GhsDescentResult`).
- `rho::ghs::transfer` — the point-to-divisor transfer map (`transfer_point`,
  `verify_homomorphism`).

### The toy KAT/fixture

The fixture is `ghs_toy_curve()`: $E/\mathbb{F}_{2^6}$ with $m = 6$, $l = 2$, $m/l = 3$ (odd —
imaginary hyperelliptic model). The source field uses irreducible $x^6 + x + 1$ (poly = 0x43);
the subfield uses $x^2 + x + 1$ (poly = 0x7). The binary curve is $y^2 + xy = x^3 + x^2 + 1$.

```rust
let curve_e = ghs_toy_curve();
let params = GhsParams::new(6, 2, curve_e.clone(), Uint::<1>::from(GHS_POLY2))
    .expect("GHS params must be valid");
let g = curve_e.generator::<F2mNaive<1>>();
let h = g.clone(); // h = 1·g
let result = ghs_descend(&params, &g, &h)
    .expect("GHS descent must succeed");
assert!(verify_log_preservation(&result, 1));
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_ghs_transfer`) asserts this before
timing. The bench measures the *descent reduction + log-preservation verification*, annotated
as a transfer.

### Cross-reference

For the mathematical development — the Weil restriction, the hyperelliptic Jacobian, and the
conditional L-notation — see `docs/MATHEMATICS.md` §10.3.

---

## 13. Semaev Summation Polynomials: The Index-Calculus Primitive

### What it exploits

The Semaev summation polynomial $S_m(X_1, \ldots, X_m)$ is a symmetric multivariate polynomial
over $\mathbb{F}_p$ that vanishes on $(x_1, \ldots, x_m)$ precisely when there exist $y_i$ such
that $P_i = (x_i, y_i) \in E(\mathbb{F}_p)$ with $P_1 + \cdots + P_m = \mathcal{O}$. This is
the combinatorial primitive at the heart of the Gaudry–Diem–Joux–Vitse index calculus.

The factor base $\mathcal{F}$ is a set of points on $E$ with small x-coordinates. A point $P$
is *decomposable* over $\mathcal{F}$ if $P = \sum_{i=1}^m F_i$ for some $F_i \in \mathcal{F}$.
The Semaev polynomial detects this: $S_m(x_{F_1}, \ldots, x_{F_{m-1}}, x_P) = 0$ iff $P$
decomposes over $\mathcal{F}$ with the given $F_i$.

### The module surface

```rust
// rho::semaev

/// Compute the m-th Semaev summation polynomial S_m via the resultant ladder.
///
/// S_3 is the base case (computed from the curve equation).
/// S_m = Res_X(S_{m-1}(X_1, ..., X_{m-2}, X), S_3(X_{m-1}, X_m, X)) for m > 3.
pub fn semaev_poly(m: usize, curve: &Curve<FpNaive>) -> Result<MultiPoly, SemaevError>;
```

The `MultiPoly` type is a symmetric multivariate polynomial over $\mathbb{F}_p$, supporting
evaluation, resultant computation, and variable substitution. The recursion ladder builds
$S_m$ from $S_3$ via iterated resultants.

### Cross-reference

The Semaev primitive is consumed by the index-calculus solver (§14 below). For the mathematical
development — the summation polynomial definition and the relation-collection loop — see
`docs/MATHEMATICS.md` §10.4.

---

## 14. Index Calculus: Semaev Decomposition over the Factor Base

### What it exploits

The index-calculus attack collects *relations*: points $a_i \cdot G + b_i \cdot Q$ that decompose
over the factor base $\mathcal{F}$ (detected via the Semaev polynomial). Each decomposition
$a_i G + b_i Q = \sum_j e_{ij} F_j$ gives a linear equation over $\mathbb{Z}/\ell\mathbb{Z}$,

$$a_i + b_i \cdot x \equiv \sum_j e_{ij} \cdot l_j \pmod{\ell},$$

in the unknown $x = \log_G(Q) \bmod \ell$ and the unknown factor-base logarithms $l_j$. The
equation holds modulo $\ell$ because $\ell \,\|\, n$: the $\ell$-Sylow subgroup of
$\langle G\rangle$ is cyclic, generated by $G_\ell = (n/\ell) \cdot G$, so every point's
$\ell$-component has a logarithm modulo $\ell$. This is an invariant of the fixture (asserted by
`toy_fixture_parameters`), not something a test proves independently.

The code does not solve for the $l_j$. It computes the *left kernel* of the exponent matrix — the
vectors $u$ with $\sum_i u_i e_{ij} \equiv 0$ for every $j$ — and any $u$ with
$u^{\top} b \not\equiv 0$ gives $x \equiv -u^{\top} a \,/\, u^{\top} b$, because the weighted sum of
the relations has no $l_j$ left in it. The value of $x$ comes only from a relation combination; no
candidate logarithm is ever tried.

**Scale boundary.** Over $E(\mathbb{F}_p)$, index calculus is NOT faster than Pollard rho — the
asymptotic win requires the extension-field setting $E(\mathbb{F}_{p^n})$ with $n > 1$ (the genuine
Gaudry–Diem setting), which this code does not implement. The toy fixture demonstrates the
mechanism; the asymptotic win is not observable at $p = 47$.

### The module surface

```rust
// rho::index_calculus

/// Collect relations and recover log_G(Q) mod ℓ from a left-kernel combination.
pub fn index_calculus_dlp(
    g: AffinePoint<FpNaive>,
    q: AffinePoint<FpNaive>,
    strategy: &IndexCalcStrategy,
) -> Result<IndexCalcSolution, IndexCalcError>;

/// The verified answer and the evidence that produced it.
pub struct IndexCalcSolution {
    pub x: u64,                    // log_G(Q) mod ℓ, verified x·G_ℓ == Q_ℓ
    pub combination: Vec<FpNaive>, // the kernel vector u used; uᵀb ≢ 0
    pub relations: Vec<Relation>,  // the relation system u combines
}

/// Left-kernel basis of the relation exponent matrix over F_ℓ: each u satisfies Σᵢ uᵢ·eᵢ ≡ 0.
pub fn solve_ek_linalg(
    relations: &[Relation],
    strategy: &IndexCalcStrategy,
) -> Result<Vec<Vec<FpNaive>>, IndexCalcError>;

/// Collect relations over the (a, b) enumeration until the system is over-determined.
pub fn collect_relations(
    g: AffinePoint<FpNaive>,
    q: AffinePoint<FpNaive>,
    strategy: &IndexCalcStrategy,
) -> Result<Vec<Relation>, IndexCalcError>;

/// Decompose a point as a sum of m factor-base points (None if no such decomposition exists).
pub fn decompose(
    q: AffinePoint<FpNaive>,
    strategy: &IndexCalcStrategy,
) -> Option<Vec<FbPoint>>;
```

`decompose` returns `Option<Vec<FbPoint>>`: `Some` holds the $m$ factor-base points (with
repetition allowed; $m = 2$ is the only implemented arity) whose sum is `q` under the group law,
and `None` means no decomposition exists over the current factor base. It returns the points
themselves, not their indices; the relation records the indices.

The `IndexCalcStrategy` bundles the curve, factor base, subgroup modulus $\ell$, and Semaev
polynomial degree $m$. The `toy()` constructor builds the canonical toy fixture.

The internal structure:
- `rho::index_calculus::strategy` — `FbPoint`, `IndexCalcStrategy`, `Relation`.
- `rho::index_calculus::collect` — relation collection (`collect_relations`, and
  `RelationStream`, which yields relations one at a time).
- `rho::index_calculus::decompose` — point decomposition via Semaev polynomial (`decompose`).
- `rho::index_calculus::linalg` — $\mathbb{F}_\ell$ linear algebra (`build_ek_matrix`,
  `solve_ek_linalg`). At this scale the kernel is computed by exact dense row reduction inside
  `block_wiedemann_fl`, not by a Krylov (Wiedemann or Lanczos) iteration; the returned basis is
  complete.
- `rho::index_calculus::solve` — full pipeline assembler (`index_calculus_dlp`).

**How the solver draws relations.** `index_calculus_dlp` starts from `fb_size + 1` relations and
draws a further `fb_size + 1` per round while no kernel vector has $u^{\top} b \not\equiv 0$. The
relations come from a deterministic `RelationStream` that visits the pairs $(a, b)$ with
$a \in [1, n)$, $b \in [0, n)$ along the anti-diagonals $a + b = 1, 2, \ldots$. That enumeration
bounds how many relations exist to be collected; it is not a search over candidate logarithms.

**Outcomes.** `index_calculus_dlp` returns `Ok` only after checking that every basis vector with
$u^{\top} b \not\equiv 0$ gives the same $x$ and that $x \cdot G_\ell = Q_\ell$ by the group law.
Its errors keep distinct classes apart:

- `Undetermined { relations }` — the relation enumeration is exhausted and no kernel vector has
  $u^{\top} b \not\equiv 0$. It is a statement about the relation set, not a negative about the
  logarithm.
- `Inconsistent` — an invariant failure: combinations disagree on $x$, verification fails, or a
  kernel vector has $u^{\top} b \equiv 0$ but $u^{\top} a \not\equiv 0$. It is surfaced, never
  retried or masked.
- `UnderdeterminedSystem`, `InvalidSubgroup`, `FactorBaseTooSmall`, `Semaev` — input and
  setup errors propagated from the earlier steps.

**Counts.** The relation count is `solution.relations.len()`. `collect_relations(...).len()` and
`decompose(...)` remain available for a caller that wants counts without running the solver.

### The toy KAT/fixture

The fixture is `IndexCalcStrategy::toy()`: $y^2 = x^3 + x + 33 \bmod 47$, $n = 60$, $\ell = 5$,
$|FB| = 6$, $m = 2$ (asserted by `toy_fixture_parameters`). The known scalar is $k = 7$:
$Q = 7 \cdot G$ (with $Q_\ell = (n/\ell) \cdot Q$ in the order-$\ell$ subgroup).

```rust
let strategy = IndexCalcStrategy::toy().expect("toy strategy must build");
let curve = strategy.curve.clone();
let g = curve.generator::<FpNaive>();
let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
let solution = index_calculus_dlp(g.clone(), q.clone(), &strategy)
    .expect("index_calculus_dlp must recover x");
let x = solution.x;
// Verify x·G_ℓ = Q_ℓ (subgroup-log correctness).
let cofactor = n / ell;
let g_ell = curve.scalar_mul(&g, &Uint::<4>::from(cofactor));
let q_ell = curve.scalar_mul(&q, &Uint::<4>::from(cofactor));
assert_eq!(curve.scalar_mul(&g_ell, &Uint::<4>::from(x)), q_ell);
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_index_calculus`) asserts this before
timing. The tests in `rho/tests/index_calculus_kat.rs` that pin the mechanism are:

- `left_kernel_annihilates_exponents` — every returned basis vector satisfies
  $\sum_i u_i e_i \equiv 0$, is nonzero, and the basis size equals #relations − rank (complete).
- `adapter_fidelity` — the matrix handed to the linear algebra is the transposed exponent matrix.
- `recovers_every_residue` — for every $k \in 0..60$, the solver returns $x = k \bmod 5$ with
  $u^{\top} b \not\equiv 0$, and $\sum_i u_i (a_i + b_i x) \equiv 0$ is checked explicitly from the
  returned combination and relations.
- `undetermined_when_no_combination_sees_b` — honest relations that all have
  $b \equiv 0 \pmod{\ell}$ give `Undetermined`.
- `corrupted_relation_is_inconsistent` — one corrupted relation gives `Inconsistent`.
- `toy_fixture_parameters` — the fixture constants above.

`Undetermined` is not reachable for an honest $Q$ on this fixture; the tests exercise it only on
constructed relation sets. Over $E(\mathbb{F}_p)$ these tests demonstrate the mechanism, not the
asymptotic advantage.

### Cross-reference

For the mathematical development — the Semaev polynomial, the relation-collection loop, the
left-kernel recovery, and the L-notation — see `docs/MATHEMATICS.md` §10.4.

---

## 15. The Attack Interfaces

Each algebraic attack exposes a small, separate interface, which lets the attacks be developed,
tested, and reasoned about independently. This section names them in one place.

| Interface | What it exposes | Used by |
|-----------|-----------------|---------|
| Pollard rho (`rho::ecdlp`) | `solve_brent`, `solve_dp`, `solve_dp_negmap`, `solve_dp_batch`, `solve_dp_glv` — the generic-$\sqrt{n}$ rho solvers, each returning `Result<DlogSolution<RhoEvidence>, DlogError>` | the benchmark baseline column; Pohlig–Hellman (rho for digits above 64) |
| Pohlig–Hellman (`rho::ecdlp::pohlig`) | `solve_ecdlp_composite` (returning `DlogSolution<PohligEvidence>`), `factor_order`, `project_to_subgroup` — composite-order ECDLP via CRT | the Pohlig–Hellman bench and code-tour |
| MOV (`rho::pairing::mov`) | `mov_reduce` — MOV/Frey–Rück reduction to $\mathbb{F}_{p^k}^*$ DLP, returning `Result<MovSolution, MovError>` | the MOV bench and code-tour; the payoff proof in `docs/MATHEMATICS.md` §10.5 |
| SSA (`rho::ssa`) | `ssa_solve` — SSA anomalous-curve attack, returning `Result<u64, SsaError>` | the SSA bench and code-tour |
| GHS (`rho::ghs`) | `ghs_descend`, `verify_log_preservation`, `transfer_point` (returning `Result<_, GhsError>`) — binary-curve transfer (NOT a solver) | the GHS bench and code-tour, annotated as a transfer |
| Semaev (`rho::semaev`) | `semaev_poly` + `MultiPoly` operations — the index-calculus primitive | index calculus; the Semaev code-tour |
| Index calculus (`rho::index_calculus`) | `index_calculus_dlp` returning `IndexCalcSolution`, `collect_relations`, `solve_ek_linalg`, `decompose` — the relation-collection and left-kernel pipeline | the index-calculus bench and code-tour |

The benches add only Criterion timing and prose; no solver surface is touched by them.

---

## 16. What Each Attack Implements and Claims

This section states, attack by attack, what the code implements, what exercises it, and what it
deliberately does not claim.

### Algorithmic content

Every attack — Pollard rho, Pohlig–Hellman, MOV, SSA, GHS, Semaev, and index calculus — is
implemented head-on, not stubbed, and exercised by KATs on a toy fixture.

- **Pollard rho.** `solve_brent`, `solve_dp`, `solve_dp_negmap`, `solve_dp_batch`,
  `solve_dp_glv` — all five solver variants are exercised by KATs, on the 20-bit `tiny_a` and
  `tiny_b` curves (and, for GLV, the 1093-element `tiny_glv` group through the private
  `solve_dp_glv_impl`). On `secp_k1_toy` itself only the budget-exhaustion and unsupported-domain
  outcomes of the GLV solver are tested (`dp_glv_tiny_budget_is_budget_exhausted`,
  `dp_glv_is_unsupported_off_its_curve`); a full solve of that 63-bit group is not practically
  runnable. Every solver returns a verified logarithm with its collision evidence or a typed error
  that keeps invalid input, unsupported domain, budget exhaustion and invariant failure apart.

- **Pohlig–Hellman.** `solve_ecdlp_composite` with `factor_order` and
  `project_to_subgroup` — composite-order ECDLP via CRT, exercised on `composite_toy()`
  ($n = 60$, every digit by subgroup enumeration) and `composite_rho_toy()` ($n = 1068$, one digit
  by rho). The result is verified and its evidence records the method of each digit.

- **MOV/Frey–Rück.** `mov_reduce` with the full pairing infrastructure (Miller's
  algorithm, reduced Tate pairing, $\mathbb{F}_{p^k}$ arithmetic, the bridge to
  `gnfs::dl::solve_dl`) — exercised on `pairing_toy()` ($\ell = 3$, $k = 2$). The result is
  verified `k·G == Q`, and a degenerate pairing is reported as such.

- **Smart–Satoh–Araki.** `ssa_solve` with Hensel lift and an evaluation of $\log(1 + t)$ on
  $p \cdot \tilde{P}$ — anomalous-curve attack, exercised on `anomalous_toy()` ($p = 7$,
  $\#E = 7$). The algorithm is polynomial in $\log p$; this implementation is $O(p)$. Every
  returned value is verified, and the solver's second candidate (the quotient times $2^{-1}$) is an
  empirical compensation for that fixture (§11).

- **GHS/Weil descent.** `ghs_descend`, `verify_log_preservation`, `transfer_point` —
  binary-curve transfer to hyperelliptic Jacobian, exercised on `ghs_toy_curve()`
  ($\mathbb{F}_{2^6}$). It is a transfer only (no `ghs_dlp`); the downstream solve on the
  Jacobian is not implemented. The transfer rejects off-curve points, curves that do not
  descend, and a subfield polynomial that does not match the subfield generator; it identifies
  $\mathbb{F}_{2^l} \subset \mathbb{F}_{2^m}$ through `shared_gf2m::restrict` (§12).

- **Semaev summation polynomials.** `semaev_poly` via the resultant ladder — the
  index-calculus primitive, exercised by `rho/tests/semaev_kat.rs`.

- **Index calculus.** `index_calculus_dlp` with `collect_relations`, `decompose`,
  `build_ek_matrix`, `solve_ek_linalg` — the relation-collection and left-kernel mechanism of
  the Gaudry–Diem–Joux–Vitse algorithm, exercised on `IndexCalcStrategy::toy()`. It recovers
  $x = \log_G(Q) \bmod \ell$ from a relation combination and verifies $x \cdot G_\ell = Q_\ell$;
  the tests that pin this are listed in §14. It demonstrates the mechanism over $E(\mathbb{F}_p)$,
  not the asymptotic advantage, and `Undetermined` is a statement about a relation set, not a
  negative about the logarithm.

### No engineering optimisation

Oracle tests stay opt-in; there is no production solver acceleration, no GPU, and no distributed
computation.

- **Oracle tests.** The oracle KATs that need an external tool (PARI/GP, msolve) are gated behind
  `#[ignore = "PARI not installed; run manually when available"]` (or the msolve analogue).
  None are on the default test path. No production dependency was added.

- **Pairing computation.** Miller's algorithm is implemented at demonstration fidelity: the
  standard Miller loop with the final exponentiation. No optimised pairing (ate pairing, optimal
  ate) was introduced.

- **Index calculus.** The relation-collection loop uses the Semaev polynomial at degree $m = 2$
  (the toy setting). No Gröbner-basis acceleration, no Weil descent optimisation, no
  large-prime variation.

- **GHS descent.** The Weil restriction and Artin–Schreier extension are implemented at
  demonstration fidelity. No optimised descent (no Magma/PARI oracle for the hyperelliptic
  Jacobian DLP).

All oracle cross-checks remain `#[ignore]`-gated.

### Scale: demonstration fidelity

The implementation is correct but not production-scale. Scale-dependent phenomena are annotated
explicitly, and the annotations document the science↔engineering gap:

| Attack | Phenomenon | Annotation |
|--------|-----------|------------|
| Pohlig–Hellman | Speedup over rho only exponential at crypto scale; at toy scale ($n = 60$) the gain is invisible | §10.1 of the textbook; the cross-attack harness in BENCHMARKS.md |
| MOV/Frey–Rück | At $k = 2$, $p = 47$: pairing + $\mathbb{F}_{p^2}$ DLP overhead dominates; asymptotic win requires crypto-scale $p$ | §10.5 of the textbook; the cross-attack harness in BENCHMARKS.md |
| SSA | Polynomial-time at all scales; toy fixture ($p = 7$) makes constant factors invisible | §10.2 of the textbook |
| GHS | Transfer only; asymptotic win comes from downstream index calculus on the Jacobian, which is not implemented | §10.3 of the textbook; the cross-attack harness in BENCHMARKS.md |
| Index calculus | Over $E(\mathbb{F}_p)$: NOT faster than rho; asymptotic win requires $E(\mathbb{F}_{p^n})$, $n > 1$ | §10.4 of the textbook; the cross-attack harness in BENCHMARKS.md |
| L-notation separations | $L[1, 1/2]$ (rho) vs $L[0]$ (SSA) vs $L[1/3]$ (MOV via NFS-DL) are NOT observable at $p = 47$/$p = 7$ | §10.6 of the textbook |

---

## 17. KAT Summary

This chapter adds no new KATs (it is a code-tour, not a new implementation). The KATs that
exercise the attacks are the existing per-attack KAT suites plus the `attacks.rs` bench
pre-check asserts (no-regression smoke tests).

### Per-attack KATs

| Test file | Tests | Scope |
|-----------|-------|-------|
| `rho/tests/ecdlp_kat.rs` | rho and Pohlig–Hellman KATs | Brent and distinguished-point solvers on `tiny_a`/`tiny_b` (verified `k`, rechecked collision evidence, signed collisions, budget exhaustion, GLV `Unsupported` off its curve); `solve_ecdlp_composite` on `composite_toy()` ($n = 60$) and `composite_rho_toy()` ($n = 1068$); CRT correctness; prime-power lift; wrong-order and outside-subgroup inputs |
| `rho/tests/binary_ecdlp_kat.rs` | binary-field and Koblitz rho KATs | verified solutions with rechecked evidence on the order-4 toy curve; composite-order input checks; degenerate collisions; τ eigenvalues; a corrupted eigenvalue ending in `Inconsistent` |
| `rho/tests/factor_kat.rs` | rho factoring KATs | known semiprimes; `Even`/`Rho` evidence; `Prime`, `InvalidInput`, `Unsupported` and `BudgetExhausted` outcomes |
| `rho/tests/mov_kat.rs` | MOV/Frey–Rück KATs | `mov_reduce` on `pairing_toy()` ($\ell = 3$, $k = 2$); pairing bilinearity; recovery of $k \bmod \ell$ for $k = 0, 1, 2$ with the pairing identity rechecked; input defects and degenerate pairings |
| `rho/tests/ssa_kat.rs` | SSA KATs | `ssa_solve` on `anomalous_toy()` ($p = 7$); Hensel lift; every $k \in 1..7$ recovered and verified; `NotAnomalous`, `Unsupported`, and $Q = \infty$ |
| `rho/tests/ghs_kat.rs` | GHS KATs | `ghs_descend` + `verify_log_preservation` on `ghs_toy_curve()` ($\mathbb{F}_{2^6}$); log-preservation; `NotOnCurve`, `NotInSubfield` and `SubfieldPolyMismatch` guards; hand-derived transfer known answers |
| `rho/tests/semaev_kat.rs` | Semaev KATs | `semaev_poly` correctness; $S_3$ base case; resultant ladder |
| `rho/tests/index_calculus_kat.rs` | Index-calculus KATs | factor base and relations; Semaev decomposition; the transposed matrix and its complete left kernel (`left_kernel_annihilates_exponents`); recovery of every residue $k \bmod 5$ for $k \in 0..60$ from a relation combination (`recovers_every_residue`); the `Undetermined` and `Inconsistent` outcomes; the fixture constants (`toy_fixture_parameters`) |

### Bench pre-check asserts (no-regression smoke tests)

Each bench in `rho/benches/attacks.rs` asserts the solver returns the known correct answer
before timing. These asserts are a no-regression invariant: every benched attack still solves
(or, for GHS, transfers) its toy instance.

| Bench function | Pre-check assert | Known answer |
|----------------|-----------------|--------------|
| `attacks/pohlig_hellman` | `k·G = Q` for $k$ returned by `solve_ecdlp_composite` | $k = 7$ on `composite_toy()` |
| `attacks/mov_frey_ruck` | the verified solution's `k_mod_ell` is $2$ | $k \bmod \ell = 2$ on `pairing_toy()` |
| `attacks/ssa` | `ssa_solve` returns `Ok(3)` | $k = 3$ on `anomalous_toy()` |
| `attacks/ghs_transfer` | `verify_log_preservation` returns `true` for $k = 1$ | $k = 1$ on `ghs_toy_curve()` |
| `attacks/index_calculus` | $x \cdot G_\ell = Q_\ell$ for the $x$ in the `IndexCalcSolution` returned by `index_calculus_dlp` | $x \cdot G_\ell = Q_\ell$ on `IndexCalcStrategy::toy()` |

---

## 18. Cross-References

### Mathematical textbook

- **`docs/MATHEMATICS.md` ch. 10** (algebraic ECDLP attacks) — the maths-first sibling to this
  integrative chapter. Develops the five algebraic ECDLP attacks mathematically:
  Pohlig–Hellman (§10.1), Smart–Satoh–Araki (§10.2), GHS/Weil descent (§10.3), index calculus
  (§10.4), and the **MOV/Frey–Rück reduction (§10.5 — the designated payoff proof)**. Carries the
  per-attack L-notation comparison (§10.6). Cross-references this chapter for the code realisation.

- **`docs/MATHEMATICS.md` §Pollard Rho for ECDLP** — the rho baseline this chapter's
  bound-breaking extends. The $\sqrt{n}$ generic bound established there is the baseline all
  algebraic attacks escape.

- **`docs/MATHEMATICS.md` §"Escape from Search: The Through-Line"** — the five-family structure
  taxonomy and the L-notation hierarchy table. The scaffold the algebraic attacks realise.

### Benchmark data

- **`docs/BENCHMARKS.md`, the cross-attack ECDLP harness** — the structural-precondition-conditional
  benchmark table (Attack / Curve-precondition / Applies? / Toy-scale cost / Escape structure).
  The empirical substrate this code-tour stands on.

### Interface definitions

- **Pollard rho:** `rho/src/ecdlp/mod.rs`.
- **Pohlig–Hellman:** `rho/src/ecdlp/pohlig.rs`.
- **MOV:** `rho/src/pairing/mov.rs`.
- **SSA:** `rho/src/ssa/mod.rs`.
- **GHS:** `rho/src/ghs/mod.rs`.
- **Semaev:** `rho/src/semaev/mod.rs`.
- **Index calculus:** `rho/src/index_calculus/mod.rs`.

---

## Further Reading

1. **Menezes, A. J., Okamoto, T., and Vanstone, S. A. (1993).** "Reducing elliptic curve
   logarithms to logarithms in a finite field." *IEEE Transactions on Information Theory*,
   39(5), 1639–1646. The original MOV reduction paper.

2. **Frey, G., and Rück, H.-G. (1994).** "A remark concerning m-divisibility and the discrete
   logarithm in the divisor class group of curves." *Mathematics of Computation*, 62(206),
   865–874. The Frey–Rück variant of the pairing reduction.

3. **Smart, N. P. (1999).** "The discrete logarithm problem on elliptic curves of trace one."
   *Journal of Cryptology*, 12(3), 193–196. The SSA polynomial-time attack.

4. **Satoh, T., and Araki, K. (1998).** "Fermat quotients and the polynomial time discrete log
   algorithm for anomalous elliptic curves." *IEICE Transactions on Fundamentals*, E81-A(6),
   1228–1233. The Satoh–Araki variant of the anomalous-curve attack.

5. **Gaudry, P., Hess, F., and Smart, N. P. (2002).** "Constructive and destructive facets of
   Weil descent on elliptic curves." *Journal of Cryptology*, 15(1), 19–46. The GHS Weil-descent
   attack.

6. **Semaev, I. (2004).** "Summation polynomials and the discrete logarithm problem on elliptic
   curves." Cryptology ePrint Archive, Report 2004/031. The summation polynomial primitive.

7. **Gaudry, P. (2009).** "Index calculus for abelian varieties of small dimension and the
   elliptic curve discrete logarithm problem." *Journal of Symbolic Computation*, 44(12),
   1690–1702. The Gaudry index-calculus algorithm.

8. **`docs/MATHEMATICS.md` ch. 10** (algebraic ECDLP attacks). The maths-first sibling to this
   chapter: the full MOV payoff proof (§10.5), the per-attack L-notation comparison (§10.6), and the
   algebraic ECDLP attacks developed mathematically from first principles.
