# rDLP — Library-Wide Code-Tour: Structure-Based Escape from Search

> **Maths-first sibling.** For the mathematical foundations of every algorithm surveyed here —
> the birthday-paradox collision argument, the $L$-notation hierarchy, the number-field bridge,
> the pairing reduction, the quantum period-finding argument, and the full taxonomy of
> structure-based escapes — see `docs/MATHEMATICS.md`. That file is the maths-first textbook;
> this file is the code-first tour. Neither is a prerequisite for the other.

This document is the library-wide code-tour for the rDLP workspace. It surveys all five tracks
of the project — ρ (Pollard rho), G (GNFS), D (NFS-DL), E (algebraic ECDLP attacks), and S
(Shor's algorithm) — as chapters in one story: the story of how each track finds a different kind
of exploitable structure that escapes the generic search bound. The code is the example; this
document is the lesson.

**How to read this document.** The umbrella synthesis (§0) names the through-line once and maps
all five tracks. The per-track chapters (§§1–7 for ρ, §§8–18 for E) are the depth for those two
tracks; the remaining three tracks are threaded here as synthesis chapters that cite their
per-crate tours for the implementation detail. The per-crate tours are:

| Track | Code-tour | Sections here |
|-------|-----------|---------------|
| ρ — Pollard rho + ECDLP | *this file* §§1–7 | Full depth |
| α — Arithmetic substrate | `shared/numth/docs/PEDAGOGY.md` | §0.3 (synthesis) |
| NF — Number-field substrate | `shared/numfield/docs/PEDAGOGY.md` | §0.3 (synthesis) |
| G+D — GNFS + NFS-DL | `gnfs/docs/PEDAGOGY.md` §§1–71 | §0.4 (synthesis) |
| E — Algebraic ECDLP attacks | *this file* §§8–18 | Full depth |
| S — Shor's algorithm | `shor/docs/PEDAGOGY.md` | §0.5 (synthesis) |

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

1. **Group homomorphisms / smoothness** — index calculus, NFS, NFS-DL (Tracks G, D).
2. **Endomorphisms** — GLV, Koblitz (Track ρ, Phase 8).
3. **Pairings** — MOV/Frey–Rück (Track E, §10).
4. **Number-field structure** — GNFS, NFS-DL (Tracks G, D).
5. **Quantum period-finding** — Shor's algorithm (Track S).

Each track in this library demonstrates one or more of these escapes in working Rust code, at toy
scale, with known-answer tests. The code is correct; the scale is not cryptographic (see
`docs/MATHEMATICS.md` §"On Scale" for the honest science↔engineering gap).

### §0.2 The comparative L-notation table

The complexity of the best known algorithms for the main problems forms a hierarchy in
$L$-notation:

$$L_N[\alpha, c] = \exp\!\left(c \cdot (\log N)^\alpha \cdot (\log \log N)^{1-\alpha}\right)$$

The table below gives the leading complexity for each track, read from the per-track derivations
in `docs/MATHEMATICS.md` (§6 for ρ, §7 for GNFS, §9.7 for NFS-DL, §10.6 for algebraic ECDLP,
§11.3.4 and §11.4.4 for Shor). **Do not recall these constants from memory — read the cited
sections.**

| Track | Algorithm | Problem | L-notation complexity | Structure exploited | Notes |
|-------|-----------|---------|----------------------|---------------------|-------|
| ρ | Pollard rho | ECDLP (generic curve) | $L_n[1,\, 1/2] = \Theta(\sqrt{n})$ | None — generic $\sqrt{n}$ walk | Fully exponential in $\log n$; the baseline all other tracks escape |
| α | ECM (Lenstra) | Integer factoring (sub-step) | $L_p[1/2,\, 1]$ (largest prime factor $p$) | Group-order smoothness | Substrate; used inside NFS large-prime variations |
| G | GNFS | Integer factoring | $L_N[1/3,\, (64/9)^{1/3}]$ | Number-field bridge + smoothness | Exponent $1/3$ from sieve/linear-algebra balance; constant from full analysis [LLMP93] |
| D | NFS-DL | DLP in $\mathbb{F}_p^*$ | $L_p[1/3,\, (64/9)^{1/3}]$ | Same number-field bridge | Same exponent and constant as GNFS; descent is asymptotically subdominant (§9.7 delta) |
| E | MOV/Frey–Rück | ECDLP (small embedding degree) | $L_{p^k}[1/3,\, (64/9)^{1/3}]$ (via NFS-DL) | Pairing → $\mathbb{F}_{p^k}^*$ DLP | Cross-track bridge: ECDLP → NFS-DL; precondition $\ell \mid p^k - 1$, $k$ small |
| E | Smart–Satoh–Araki | ECDLP (anomalous curve) | $L_p[0] = O(\log p)$ (polynomial time) | $p$-adic formal group logarithm | Sharpest escape among the algebraic attacks; precondition $\#E(\mathbb{F}_p) = p$ |
| E | Index calculus | ECDLP ($E/\mathbb{F}_{p^n}$, $n > 1$) | $L_{p^n}[1/2,\, c]$ (subexponential) | Factor-base decomposability | Asymptotic win requires extension-field setting; toy fixture is $n = 1$ |
| E | Pohlig–Hellman | ECDLP (composite order) | $L_{p_{\max}}[1,\, 1/2]$ (largest prime factor) | CRT group-order homomorphism | Reduces to prime-order subgroup DLPs; speedup requires $p_{\max} \ll n$ |
| S | Shor (factoring) | Integer factoring | $O((\log N)^3) = L_N[0,\, 3]$ | Quantum period-finding | Polynomial in $\log N$; dissolves the $L$-notation framework entirely |
| S | Shor (ECDLP) | ECDLP | $O((\log r)^3) = L_r[0,\, 3]$ | Quantum period-finding (2-register HSP) | Same polynomial collapse; ~768 qubits for secp256k1 |

**The quantitative spine of the library.** Reading the table top-to-bottom is the story of the
project: from the generic $\sqrt{n}$ baseline (ρ, $\alpha = 1$) through the subexponential
classical attacks (G/D at $\alpha = 1/3$; index calculus at $\alpha = 1/2$; MOV bridging ECDLP
to the $\alpha = 1/3$ regime) through the polynomial-time anomalous-curve attack (SSA, $\alpha = 0$
classically) to the quantum dissolution of the $L$-notation framework entirely (Shor, $\alpha = 0$
polynomially). Each step down the $\alpha$ hierarchy is achieved by finding a new kind of
exploitable structure.

**Scale honesty.** The $L$-notation separations in this table are NOT observable at the toy
scale of the project's fixtures. The table reports theoretical complexity classes; the toy-scale
timings are in `docs/BENCHMARKS.md`. The separation between $L[1, 1/2]$ (rho), $L[1/3]$ (GNFS/
NFS-DL/MOV), $L[1/2]$ (index calculus), and $L[0]$ (SSA, Shor) is a statement about asymptotic
behaviour at cryptographic scale — not a ranking of toy-scale timings.

### §0.3 The shared substrate: α and NF chapters

Before the five tracks diverge, two shared substrate crates provide the arithmetic and algebraic
infrastructure that every track sits on.

**The α-substrate** (`shared/numth/docs/PEDAGOGY.md`). Three crates — `shared-field`,
`shared-bigint`, and `shared-numth` — provide prime-field arithmetic, number-theoretic primitives,
and cross-cutting helpers. The key algorithms: Miller–Rabin primality testing (needed to build
factor bases and certify found factors), B-smoothness detection (the engine of every sieve),
Lenstra's ECM (factoring by group-order smoothness — a direct application of the group-homomorphism
structure, used inside NFS large-prime variations and Pohlig–Hellman), and Tonelli–Shanks square
roots (needed for curve-point recovery and the GNFS square-root step). The α-substrate does not
itself escape the generic bound — it provides the tools that other algorithms use to do so. For
the full code-tour, see `shared/numth/docs/PEDAGOGY.md`; for the mathematical development, see
`docs/MATHEMATICS.md` §"The α-Substrate: Primality, Smoothness, and ECM".

**The number-field substrate** (`shared/numfield/docs/PEDAGOGY.md`). The `shared-numfield` crate
provides the algebraic arithmetic of a number field $K = \mathbb{Q}(\alpha)$: polynomial
arithmetic (`IntPoly`, `RatPoly`), number-field elements and the norm map, the resultant and
subresultant GCD, ideal representation, and Dedekind factorisation (including bad primes and the
Dedekind criterion). This substrate is the algebraic engine of Tracks G and D: the GNFS and
NFS-DL algorithms both exploit the norm map $N_{K/\mathbb{Q}}$ to connect smoothness in $K$ to
smoothness in $\mathbb{Z}$. For the full code-tour, see `shared/numfield/docs/PEDAGOGY.md`; for
the mathematical development, see `docs/MATHEMATICS.md` §GNFS §3 "The Number-Field Bridge".

### §0.4 Tracks G and D: GNFS and NFS-DL

**The code-tour for Tracks G and D is `gnfs/docs/PEDAGOGY.md` §§1–71.** This section is the
synthesis chapter; the per-crate tour is the depth.

Both tracks live in the `gnfs` crate and exploit the same structure: the number-field bridge. The
integers $\mathbb{Z}$ embed into a number field $K = \mathbb{Q}(\alpha)$, and the norm map
$N_{K/\mathbb{Q}}$ connects factorisation in $K$ to factorisation in $\mathbb{Z}$. By sieving for
pairs $(a, b)$ whose norms on both sides are $B$-smooth, collecting enough such *relations*, and
solving a linear system, both algorithms extract the answer — a factor (Track G) or a discrete
logarithm (Track D).

**Track G — GNFS** (`gnfs/docs/PEDAGOGY.md` §§1–62). The General Number Field Sieve factors a
large integer $N$ in five stages: polynomial selection (§§1–11), sieving (§§12–21), filtering
(§§22–30), linear algebra over GF(2) (§§31–41), and square-root recovery (§§42–51). The
integrative chapter (§§52–62) surveys the full pipeline end-to-end. The complexity is
$L_N[1/3, (64/9)^{1/3}]$ — the derivation is the designated payoff proof in
`docs/MATHEMATICS.md` §7.

**Track D — NFS-DL** (`gnfs/docs/PEDAGOGY.md` §§63–71). The NFS-DL algorithm adapts GNFS to the
discrete-logarithm problem in $\mathbb{F}_p^*$. The three DL-specific additions are: Schirokauer
maps (resolving the unit-group obstruction that replaces GNFS's class-group obstruction),
$\mathbb{F}_\ell$ linear algebra (replacing GF(2)), and individual-logarithm descent (a
special-$q$ sieve that rewrites each target as a combination of factor-base elements). The
complexity is $L_p[1/3, (64/9)^{1/3}]$ — the same exponent and constant as GNFS, because the
descent is asymptotically subdominant (the delta proof is in `docs/MATHEMATICS.md` §9.7).

**The shared substrate threading G and D.** Both tracks consume the α-substrate (smoothness
detection, ECM for large-prime variations) and the number-field substrate (polynomial arithmetic,
norm map, ideal factorisation, Dedekind criterion). The `PolyPair` type is the load-bearing
interface: it carries the polynomial pair $(f, g)$ with the shared-root invariant
$f(m) \equiv 0 \pmod{N}$, and is consumed unchanged by the NFS-DL relation collection.

### §0.5 Track S: Shor's Algorithm

**The code-tour for Track S is `shor/docs/PEDAGOGY.md`.** This section is the synthesis chapter;
the per-crate tour is the depth.

Track S is the fifth and final escape from search: quantum period-finding. Where the classical
tracks find algebraic structure that *reduces* the search cost, Shor's algorithm finds a quantum
period that **dissolves** the $L$-notation bound to polynomial time. The complexity drops from
$L_N[1/3, (64/9)^{1/3}]$ (GNFS) to $O((\log N)^3)$ (Shor factoring) and from $L_n[1, 1/2]$
(Pollard rho) to $O((\log r)^3)$ (Shor ECDLP): from subexponential or fully exponential to
polynomial in the bit-length.

Track S is organised into three parts (`shor/docs/PEDAGOGY.md` §§1–8):

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

### §0.6 Cross-track connections

Three connections thread the five tracks into one story.

**The MOV bridge (E → D).** The MOV/Frey–Rück reduction (Track E, §10) is the cross-track bridge
of the library. Given an elliptic curve with small embedding degree $k$ ($\ell \mid p^k - 1$), the
Weil/Tate pairing $e: E[\ell] \times E[\ell] \to \mu_\ell \subset \mathbb{F}_{p^k}^*$ transports
the ECDLP to a DLP in $\mathbb{F}_{p^k}^*$, where NFS-DL (Track D) applies. In code: `mov_reduce`
evaluates the reduced Tate pairing twice, encodes both outputs via the bridge encoder
(`fpext_to_bigint`), and calls `gnfs::dl::solve_dl` (the `gnfs` crate's discrete-log entry
point). The MOV bridge is the one place in the library where two tracks are directly composed at the code
level.

**The shared substrate threading G, D, and E.** The α-substrate (smoothness detection, ECM) and
the number-field substrate (polynomial arithmetic, norm map, ideal factorisation) are consumed by
all three classical-algebraic tracks. Track E's index-calculus attack (§§13–14) uses the same
smoothness-detection infrastructure as Tracks G and D; the Semaev summation polynomial is the
ECDLP analogue of the NFS factor base. The substrate is specified once, for all consumers —
the substrate-over-specifies discipline that runs through the whole project.

**The classical → quantum arc (ρ/G/D/E → S).** The five tracks form a historical and conceptual
arc. Pollard rho (1978) established the $\sqrt{n}$ generic bound and the first optimisations
(distinguished points, negation map, GLV). GNFS (1993) broke the $L[1/2]$ quadratic-sieve barrier
with the number-field bridge, reaching $L[1/3]$. NFS-DL extended the same structure to discrete
logarithms. The algebraic ECDLP attacks (Track E) showed that the $\sqrt{n}$ ECDLP bound breaks
whenever the curve has exploitable structure — pairing, anomalous order, binary field tower, or
factor-base decomposability. Shor's algorithm (1994) dissolved the $L$-notation framework
entirely in the quantum model. The arc is: generic search → subexponential classical → polynomial
quantum. Track S is the honest coda: the post-quantum migration landscape (NIST PQC, the
SIDH/Castryck–Decru break) shows that structure is a double-edged sword — the same structural
richness that enables escapes also enables attacks on post-quantum candidates.

**The honest terminus.** The whole library demonstrates the *mathematics* of each escape at toy
scale, never at cryptographically relevant scale. The toy fixtures (secp_k1_toy at 63 bits,
GNFS on small semiprimes, NFS-DL at $p = 11$, Shor factoring $N \leq 91$) are correct
demonstrations of the algorithms; they are not attacks. The $L$-notation separations are
asymptotic statements; the toy-scale timings are in `docs/BENCHMARKS.md`. This scale honesty
runs through every chapter.

---

# Chapter ρ — Pollard Rho: Phase-by-Phase

> **Maths-first treatment.** For the mathematical foundations — the birthday-paradox collision
> argument, Floyd's cycle detection, the group-homomorphism structure that makes rho work on
> elliptic curves, orbit-collapsing maps, and the $L$-notation bound — see
> `docs/MATHEMATICS.md §Pollard Rho for ECDLP`. This chapter is the code-tour sibling: it assumes
> the reader knows the mathematics and focuses on the phase-by-phase implementation in Rust.

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
canonical optimizations, built up phase by phase. Each phase adds one optimization, delivers a
runnable artifact, and produces a measurable speedup. The goal is a clean, readable CPU
implementation that a researcher can extend or benchmark. It is not a constant-time cryptographic
library; these are *attacks on* cryptographic problems, run on toy parameters.

---

## 2. Phase-by-Phase Narrative

### Phase 0 — Crate skeleton and `FpNaive`

**What it adds:** the crate scaffolding and a schoolbook modular arithmetic implementation.

The `Fp` trait defines the field operations needed throughout: `add`, `sub`, `mul`, `square`,
`inv`, `pow`, `from_u64`, `to_uint`, `zero`, `one`. `FpNaive` implements this trait using
schoolbook arithmetic on `crypto_bigint::Uint<4>` — straightforward reduction after every
operation. This is correct but slow: every multiplication requires a full modular reduction.

**Code reference:** `src/field/naive.rs` (the `FpNaive` implementation).

**Speedup:** baseline. No optimization yet; this phase establishes correctness.

---

### Phase 1 — `FpMonty` and the Montgomery-form speedup

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
speedup propagates to every subsequent phase because all curve arithmetic is built on top of it.

---

### Phase 2 — Integer factorization rho

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

### Phase 3 — Curve arithmetic

**What it adds:** elliptic curve group law in affine and Jacobian projective coordinates.

An elliptic curve over GF(p) in short Weierstrass form is y² = x³ + ax + b. The group law (point
addition and doubling) is defined by explicit formulas. Affine coordinates are conceptually simple
but require a field inversion per addition — the dominant cost. Jacobian projective coordinates
represent a point (X : Y : Z) with affine coordinates (X/Z², Y/Z³), allowing addition and doubling
with only multiplications and squarings (no inversions), at the cost of a single inversion when
converting back to affine.

Two curves are defined: a generic Weierstrass curve over a ~64-bit prime (used for Phases 4–7), and
`secp_k1_toy` — a 63-bit prime-field curve with a = 0 and an explicit order-3 endomorphism (used
for Phase 8). The secp_k1_toy curve mirrors the structure of secp256k1 (y² = x³ + 7, j-invariant
0) but over a much smaller prime, making it suitable for pedagogical experiments.

**Code reference:** `src/curve/mod.rs` (group law), `src/curve/secp_k1_toy.rs` (GLV curve
parameters and derivation).

**Speedup:** baseline for ECDLP. The choice of Jacobian vs affine coordinates becomes important in
Phase 7.

---

### Phase 4 — ECDLP rho baseline: r-adding walk and Brent's cycle detection

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

The degenerate case b_h = b_t mod n (probability ~1/n) triggers a retry with a fresh random walk
table and starting point.

**Code reference:** `src/ecdlp/walk.rs::AddendTable`, `src/ecdlp/walk.rs::WalkState::step`,
`src/ecdlp/mod.rs::solve_brent`.

**Speedup:** O(√n) expected walk steps — this is the theoretical baseline. All subsequent phases
reduce the constant factor.

---

### Phase 5 — Distinguished points and parallel collision search

**What it adds:** the van Oorschot–Wiener (vOW) parallel rho architecture.

Brent's cycle detection is inherently single-threaded: the tortoise and hare must stay in sync.
Distinguished points (van Oorschot and Wiener 1999) break this dependency. A point is
*distinguished* when its x-coordinate has at least θ low-order zero bits — a property that can be
checked independently by any walker. The expected number of steps between distinguished points is
2^θ.

The architecture is a client/server design: N walker threads each run an independent r-adding walk
and emit a `DpRecord` (carrying the current x-coordinate and the accumulated scalars a, b) whenever
the walk lands on a distinguished point. A single coordinator thread maintains a hash table keyed
on the x-coordinate. When two DPs with the same x-coordinate but different (a, b) pairs arrive,
the coordinator solves the linear system for k and signals all walkers to stop.

The key insight is that two walks that ever visit the same point will, from that point on, follow
identical trajectories and will therefore hit the same distinguished points. So a collision in the
DP table implies a collision in the underlying walks, which implies the DLP can be recovered.

The θ parameter trades memory against time: larger θ means fewer DPs stored (less memory) but more
steps per walker between DPs (more time before a collision is detected). The textbook choice is
θ such that 2^θ ≈ √n / (N · target-walks), giving ~1000 DPs at solution time.

**Code reference:** `src/ecdlp/dp.rs::is_distinguished`, `src/ecdlp/dp.rs::DpRecord`,
`src/ecdlp/coordinator.rs`, `src/ecdlp/mod.rs::solve_dp`.

**Speedup:** N× (linear in thread count). This is the second pedagogical moment: the speedup from
parallelism is linear, not √N. Each additional walker independently explores the group, and the
birthday bound applies to the combined pool of DPs. Doubling the number of walkers halves the
expected time to collision.

---

### Phase 6 — Negation map and fruitless-cycle escape

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

**Speedup:** ~√2 ≈ 1.41× reduction in expected walk steps vs Phase 5. The fruitless-cycle
overhead is small in practice (escapes are rare for a well-mixing walk), so the net speedup is
close to the theoretical √2.

---

### Phase 7 — Batched field inversion and affine coordinates

**What it adds:** 2–4× reduction in per-step cost by amortising field inversions across a batch of
walks.

In the Jacobian walk (Phase 4), each step requires one field inversion to convert the current point
to affine coordinates for the partition function. Field inversion is the most expensive field
operation — typically 50–100× the cost of a multiplication. This cost dominates the per-step
budget.

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

### Phase 8 — GLV endomorphism and the 6-orbit canonical map

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
Phase 7 is reused: each thread still owns a `BatchedWalker` with B walk states.

Fruitless cycles can still occur in the GLV walk (the canonical walk can oscillate between two
orbit members that share the same canonical representative). The BKNS doubling escape handles this
identically to Phase 6.

**Code reference:** `src/ecdlp/glv.rs::glv_phi`, `src/ecdlp/glv.rs::glv_canonical`,
`src/curve/secp_k1_toy.rs` (curve constants β, λ), `src/ecdlp/mod.rs::solve_dp_glv`.

**Speedup:** ~√3 ≈ 1.73× reduction in expected walk steps vs Phase 6 (negmap alone), or ~√6 ≈
2.45× vs Phase 5 (plain DP walk). The GLV endomorphism is essentially free to evaluate — one
multiplication by β per orbit member — so the theoretical speedup is nearly fully realised in
practice.

---

## 3. Cumulative Speedup Table

The table below summarises the expected speedup of each phase relative to the previous phase and
relative to the Phase 4 single-threaded baseline. Speedups marked with √ are reductions in
expected walk length (and therefore in expected time, assuming constant per-step cost); the
batched-inversion speedup is a reduction in per-step cost rather than walk length.

| Phase | Optimization                     | Speedup vs prev        | Cumulative vs Phase 4  |
|-------|----------------------------------|------------------------|------------------------|
| 4     | r-adding walk + Brent (baseline) | —                      | 1×                     |
| 5     | Distinguished points + N walkers | N× (linear parallelism)| N×                     |
| 6     | Negation map + BKNS escape       | ~√2 ≈ 1.41×            | ~N√2                   |
| 7     | Batched inversion + affine walks | ~2–4× (per-step cost)  | ~N√2 · B-factor        |
| 8     | GLV 6-orbit canonical map        | ~√3 vs Phase 6 alone   | ~N√6 · B-factor        |

Notes:

- The Phase 5 speedup is N (linear in thread count) because each additional walker independently
  contributes DPs to the collision pool. This is not a constant-factor improvement — it scales with
  hardware.
- Phases 6 and 8 give √2 and √3 reductions in expected walk *length*. These multiply together
  because they reduce the effective group size independently: Phase 6 halves it (factor of 2),
  Phase 8 reduces it by 3 (factor of 3), for a combined factor of 6 and a combined walk-length
  reduction of √6.
- Phase 7 reduces per-step *cost* rather than walk length. Its speedup factor (the "B-factor")
  depends on the ratio of inversion cost to multiplication cost on the target hardware, and on the
  batch size B. It composes multiplicatively with the walk-length reductions from Phases 6 and 8.

---

## 4. How Speedups Compose

The speedups from Phases 5, 6, 7, and 8 multiply rather than add because they act on independent
axes of the total computation cost.

The total expected work is (expected walk steps) × (cost per step). Phase 5 reduces expected walk
steps by a factor of N (more walkers means the birthday collision happens sooner in wall time).
Phases 6 and 8 each reduce expected walk steps by reducing the effective group size: Phase 6 by a
factor of 2 (negation map), Phase 8 by a factor of 3 (GLV orbit). Because these reductions apply
to the same quantity — expected walk steps — they multiply: the combined reduction is √(2 · 3) =
√6. Phase 7 reduces the cost per step by amortising inversions; this multiplies with the walk-step
reductions because it is a separate factor.

The practical ceiling is set by the birthday bound itself: no matter how many optimizations are
applied, the expected number of *collisions* needed is O(1), and each collision requires O(√(n/m))
walk steps where m is the effective group size after all orbit-collapsing maps. For the GLV curve
with negation map, m = n/6, so the expected steps per walker is O(√(n/6)). With N walkers, the
expected wall-time steps is O(√(n/6) / N). Beyond this, further speedup requires either a larger
orbit-collapsing map (not available for this curve family) or more parallelism.

---

## 5. Integer Factorization Rho (Phases 0–2)

The factoring side of the codebase (Phases 0–2) implements three variants of Pollard rho for
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

This runs the full ECDLP benchmark suite, measuring walk throughput (steps/second) and time-to-
solution for each solver variant: `solve_brent` (Phase 4), `solve_dp` (Phase 5), `solve_dp_negmap`
(Phase 6), `solve_dp_batch` (Phase 7), and `solve_dp_glv` (Phase 8). The benchmark output shows
the cumulative speedup as each optimization is added.

### CLI: stepping through the optimization layers

The `rho-dlog` binary exposes the solver variants via command-line flags. The following sequence
demonstrates the progression from the single-threaded baseline to the fully-optimized solver:

```sh
# Phase 4: single-threaded Brent baseline
rho-dlog --curve secp-toy --walkers 1

# Phase 5: 4 parallel walkers with distinguished points
rho-dlog --curve secp-toy --walkers 4

# Phase 7: 4 walkers with batched inversion (batch size 16)
rho-dlog --curve secp-toy --walkers 4 --batch-size 16

# Phase 8: 4 walkers, batch size 16, GLV endomorphism (secp-toy only)
rho-dlog --curve secp-toy --walkers 4 --batch-size 16 --glv
```

The `--theta` flag controls the distinguished-point threshold (default: auto-tuned to give ~1000
DPs at solution time). Increasing θ reduces memory usage at the cost of more steps per walker;
decreasing it increases memory usage but detects collisions sooner.

```sh
# Explicit theta: 1-in-256 DPs (theta=8)
rho-dlog --curve secp-toy --walkers 4 --batch-size 16 --glv --theta 8
```

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
   deterministic fruitless-cycle escape (the BKNS doubling perturbation) used in Phase 6.

5. **Gallant, R. P., Lambert, R. J., and Vanstone, S. A. (2001).** "Faster point multiplication on
   elliptic curves with efficient endomorphisms." *Advances in Cryptology — CRYPTO 2001*, LNCS
   2139, 190–200. The GLV decomposition; the endomorphism used in Phase 8 is the order-3 special
   case.

6. **Montgomery, P. L. (1987).** "Speeding the Pollard and elliptic curve methods of
   factorization." *Mathematics of Computation*, 48(177), 243–264. Introduces both the batched-GCD
   trick (Phase 2) and the batched-inversion trick (Phase 7).

7. **Bernstein, D. J., and Lange, T. (2007).** "Faster addition and doubling on elliptic curves."
   *Advances in Cryptology — ASIACRYPT 2007*, LNCS 4833, 29–50. A comprehensive reference for
   explicit addition formulas in various coordinate systems, including the affine formulas used in
   Phase 7.

8. **Brent, R. P. (1980).** "An improved Monte Carlo factorization algorithm." *BIT*, 20(2),
   176–184. Brent's cycle detection algorithm, used in both the factoring (Phase 2) and ECDLP
   (Phase 4) solvers.

---

# Chapter E — Algebraic ECDLP Attacks: An Integrative Chapter

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
| Pollard rho (baseline) | None — generic group | `rho::ecdlp` (`solve_brent`, `solve_dp`, `solve_dp_negmap`, `solve_dp_batch`, `solve_dp_glv`) | `secp_k1_toy` (63-bit, prime order) |
| Pohlig–Hellman | Composite group order $n = \prod p_i^{e_i}$ | `rho::ecdlp::pohlig` (`solve_ecdlp_composite`) | `composite_toy()` ($n = 60 = 2^2 \cdot 3 \cdot 5$) |
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
relation matrix whose solution recovers the discrete log.

---

## 9. Pohlig–Hellman: CRT to Prime-Order Subgroups

### What it exploits

The group order $n = \#E(\mathbb{F}_p)$ is composite: $n = \prod p_i^{e_i}$. The ECDLP in a
group of composite order reduces to independent DLPs in each prime-power subgroup via the
Chinese Remainder Theorem. Each prime-order subgroup DLP is solved by the rho solvers
(which require prime order). The CRT reassembles the answer.

The escape: instead of running rho on the full group of order $n$, run rho on the largest
prime-power subgroup of order $p_{\max}^{e_{\max}}$. The cost drops from $O(\sqrt{n})$ to
$O(\sum e_i (\log n + \sqrt{p_i}))$ — a dramatic speedup when $n$ has small prime factors.

### The module surface

```rust
// rho::ecdlp::pohlig

/// Decompose a u64 group order into prime-power factors.
pub fn factor_order(n: u64) -> Vec<(u64, u32)>;

/// Project (G, Q) to the order-p^e subgroup via [n / p^e]-scalar-multiplication.
pub fn project_to_subgroup<F: Fp>(
    curve: &Curve<F>,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    p_power: u64,
) -> (AffinePoint<F>, AffinePoint<F>);

/// Composite-order ECDLP via prime-power lift + CRT combine.
pub fn solve_ecdlp_composite<F: Fp>(
    curve: &Curve<F>,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Option<u64>;
```

`solve_ecdlp_composite` is the entry point. It calls `factor_order` to decompose $n$, projects
to each prime-power subgroup via `project_to_subgroup`, solves each subgroup DLP by calling the
frozen `solve_brent` (which requires prime order), and reassembles via CRT.

### The toy KAT/fixture

The fixture is `composite_toy()`: $y^2 = x^3 + x + 33 \bmod 47$, $n = 60 = 2^2 \cdot 3 \cdot 5$.
The generator $G = (10, 3)$. The known scalar is $k = 7$: $Q = 7 \cdot G$.

```rust
let curve = composite_toy();
let g = curve.generator::<FpMonty>();
let q = curve.scalar_mul(&g, &Uint::<4>::from(7u64));
let k = solve_ecdlp_composite(&curve, &g, &q, COMPOSITE_TOY_N)
    .expect("Pohlig–Hellman must succeed on composite_toy");
assert_eq!(curve.scalar_mul(&g, &Uint::<4>::from(k)), q);
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_pohlig_hellman`) asserts this before
timing.

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
/// μ_ℓ-generator with e(G, R) ≠ 1, returns k mod ℓ.
pub fn mov_reduce<F: Fp<4>>(
    curve: &Curve<FpNaive4<4>>,
    modulus: &IrreducibleModulus<F>,
    p_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    r_point: &PairingPoint<F>,
    ell: u64,
) -> Result<u64, MovError>;
```

`mov_reduce` evaluates the reduced Tate pairing twice (via `reduced_tate` in `rho::pairing::tate`),
encodes both outputs through the bridge encoder (`fpext_to_bigint`), and calls `gnfs::dl::solve_dl`
to recover $k \bmod \ell$.

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
let result = mov_reduce::<FpNaive4<4>>(&curve, &modulus, &p_point, &q_prime, &q_point, ell)
    .expect("MOV must succeed on pairing_toy");
assert_eq!(result, 2u64);
```

The bench pre-check in `rho/benches/attacks.rs` (`bench_mov`) asserts this before timing.

### Cross-reference

For the full proof of the MOV reduction — the bilinearity, non-degeneracy, and the reduction
argument — see `docs/MATHEMATICS.md` §10.5 (the designated payoff proof).

---

## 11. Smart–Satoh–Araki: The Anomalous-Curve Polynomial-Time Attack

### What it exploits

The curve is *anomalous*: $\#E(\mathbb{F}_p) = p$ (trace of Frobenius = 1). The p-adic formal
group logarithm lifts the ECDLP to $\mathbb{Z}_p$ via Hensel's lemma, where it is trivially
linear. The escape is polynomial time — the sharpest escape among the algebraic attacks.

The key steps: (1) lift the affine point $G \in E(\mathbb{F}_p)$ to a point $\tilde{G}$ on the
formal group $\hat{E}(\mathbb{Z}_p)$ via Hensel's lemma; (2) apply the formal group logarithm
$\log_{\hat{E}}$ to get $\log_{\hat{E}}(\tilde{G}) \in p\mathbb{Z}_p$; (3) the ECDLP scalar
$k$ is recovered as $\log_{\hat{E}}(\tilde{Q}) / \log_{\hat{E}}(\tilde{G}) \bmod p$.

### The module surface

```rust
// rho::ssa

/// Solve the ECDLP on an anomalous curve via the p-adic formal group logarithm.
///
/// Requires #E(F_p) = p (anomalous). Returns k mod p such that Q = k·G.
pub fn ssa_solve<F: Fp>(
    curve: &Curve<F>,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    p: u64,
) -> Result<u64, SsaError>;
```

The internal structure:
- `rho::ssa::lift` — Hensel lift: $F_p$ point → $\mathbb{Z}_p$ point (`hensel_lift_point`).
- `rho::ssa::formal_log` — formal group logarithm: $\hat{E}(\mathbb{Z}_p) \to p\mathbb{Z}_p$
  (`formal_group_log`).
- `rho::ssa::reduce` — the full SSA reduction assembling the three steps (`ssa_solve`).

The p-adic substrate (arithmetic, Hensel lifting, and the p-adic logarithm) is provided by
`shared-padic`.

### The toy KAT/fixture

The fixture is `anomalous_toy()`: $y^2 = x^3 + 5 \bmod 7$, $p = 7$, $\#E = 7 = p$. The
generator $G = (3, 2)$. The known scalar is $k = 3$: $Q = 3 \cdot G$.

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

For the mathematical development — the formal group logarithm, the Hensel lift, and the
polynomial-time complexity — see `docs/MATHEMATICS.md` §10.2.

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
    g: &BinaryPoint<F2mNaive<1>>,
    h: &BinaryPoint<F2mNaive<1>>,
) -> Result<GhsDescentResult, GhsError>;

/// Verify that the descent preserves the discrete logarithm.
/// Returns true iff log_{D_g}(D_h) = k (the known scalar).
pub fn verify_log_preservation(result: &GhsDescentResult, k: u64) -> bool;

/// Transfer a single point on E to a divisor on Jac(C).
pub fn transfer_point(
    params: &GhsParams,
    point: &BinaryPoint<F2mNaive<1>>,
) -> Result<HyperellipticDivisor, GhsError>;
```

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
$\langle G\rangle$ is cyclic, generated by $G_\ell = (n/\ell) \cdot G$, so every point's $\ell$-component has a
logarithm modulo $\ell$. This is an invariant of the fixture (asserted by `toy_fixture_parameters`),
not something a test proves independently.

The code does not solve for the $l_j$. It computes the *left kernel* of the exponent matrix — the
vectors $u$ with $\sum_i u_i e_{ij} \equiv 0$ for every $j$ — and any $u$ with
$u^{\top} b \not\equiv 0$ gives $x \equiv -u^{\top} a \,/\, u^{\top} b$, because the weighted sum of
the relations has no $l_j$ left in it. The value of $x$ comes only from a relation combination; no candidate logarithm is ever
tried.

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
asks for `fb_size` more per round while no kernel vector has $u^{\top} b \not\equiv 0$. The
relations come from a deterministic `RelationStream` that visits the pairs $(a, b)$ with $a \in [1, n)$,
$b \in [0, n)$ along the anti-diagonals $a + b = 1, 2, \ldots$. That enumeration bounds how many
relations exist to be collected; it is not a search over candidate logarithms.

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
| Pollard rho (`rho::ecdlp`) | `solve_brent`, `solve_dp`, `solve_dp_negmap`, `solve_dp_batch`, `solve_dp_glv` — the generic-$\sqrt{n}$ rho solvers | the benchmark baseline column; Pohlig–Hellman (rho as sub-solver) |
| Pohlig–Hellman (`rho::ecdlp::pohlig`) | `solve_ecdlp_composite`, `factor_order`, `project_to_subgroup` — composite-order ECDLP via CRT | the Pohlig–Hellman bench and code-tour |
| MOV (`rho::pairing::mov`) | `mov_reduce` — MOV/Frey–Rück reduction to $\mathbb{F}_{p^k}^*$ DLP | the MOV bench and code-tour; the payoff proof in `docs/MATHEMATICS.md` §10.5 |
| SSA (`rho::ssa`) | `ssa_solve` — SSA polynomial-time anomalous-curve attack | the SSA bench and code-tour |
| GHS (`rho::ghs`) | `ghs_descend`, `verify_log_preservation`, `transfer_point` — binary-curve transfer (NOT a solver) | the GHS bench and code-tour, annotated as a transfer |
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
  `solve_dp_glv` — all five solver variants are exercised by KATs on `secp_k1_toy`.

- **Pohlig–Hellman.** `solve_ecdlp_composite` with `factor_order` and
  `project_to_subgroup` — composite-order ECDLP via CRT, exercised on `composite_toy()`
  ($n = 60$).

- **MOV/Frey–Rück.** `mov_reduce` with the full pairing infrastructure (Miller's
  algorithm, reduced Tate pairing, $\mathbb{F}_{p^k}$ arithmetic, the bridge to
  `gnfs::dl::solve_dl`) — exercised on `pairing_toy()` ($\ell = 3$, $k = 2$).

- **Smart–Satoh–Araki.** `ssa_solve` with Hensel lift and formal group logarithm —
  polynomial-time anomalous-curve attack, exercised on `anomalous_toy()`
  ($p = 7$, $\#E = 7$).

- **GHS/Weil descent.** `ghs_descend`, `verify_log_preservation`, `transfer_point` —
  binary-curve transfer to hyperelliptic Jacobian, exercised on `ghs_toy_curve()`
  ($\mathbb{F}_{2^6}$). It is a transfer only (no `ghs_dlp`); the downstream solve on the
  Jacobian is not implemented.

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
| `rho/tests/pohlig_kat.rs` | Pohlig–Hellman KATs | `solve_ecdlp_composite` on `composite_toy()` ($n = 60$); CRT correctness; prime-power lift |
| `rho/tests/mov_kat.rs` | MOV/Frey–Rück KATs | `mov_reduce` on `pairing_toy()` ($\ell = 3$, $k = 2$); pairing bilinearity; $k = 2$ recovery |
| `rho/tests/ssa_kat.rs` | SSA KATs | `ssa_solve` on `anomalous_toy()` ($p = 7$); Hensel lift; formal group log; $k = 3$ recovery |
| `rho/tests/ghs_kat.rs` | GHS KATs | `ghs_descend` + `verify_log_preservation` on `ghs_toy_curve()` ($\mathbb{F}_{2^6}$); log-preservation |
| `rho/tests/semaev_kat.rs` | Semaev KATs | `semaev_poly` correctness; $S_3$ base case; resultant ladder |
| `rho/tests/index_calculus_kat.rs` | Index-calculus KATs | factor base and relations; Semaev decomposition; the transposed matrix and its complete left kernel (`left_kernel_annihilates_exponents`); recovery of every residue $k \bmod 5$ for $k \in 0..60$ from a relation combination (`recovers_every_residue`); the `Undetermined` and `Inconsistent` outcomes; the fixture constants (`toy_fixture_parameters`) |

### Bench pre-check asserts (no-regression smoke tests)

Each bench in `rho/benches/attacks.rs` asserts the solver returns the known correct answer
before timing. These asserts are a no-regression invariant: every benched attack still solves
(or, for GHS, transfers) its toy instance.

| Bench function | Pre-check assert | Known answer |
|----------------|-----------------|--------------|
| `attacks/pohlig_hellman` | `k·G = Q` for $k$ returned by `solve_ecdlp_composite` | $k = 7$ on `composite_toy()` |
| `attacks/mov_frey_ruck` | `mov_reduce` returns $k = 2$ | $k = 2$ on `pairing_toy()` |
| `attacks/ssa` | `ssa_solve` returns $k = 3$ | $k = 3$ on `anomalous_toy()` |
| `attacks/ghs_transfer` | `verify_log_preservation` returns `true` for $k = 1$ | $k = 1$ on `ghs_toy_curve()` |
| `attacks/index_calculus` | $x \cdot G_\ell = Q_\ell$ for the $x$ in the `IndexCalcSolution` returned by `index_calculus_dlp` | $x \cdot G_\ell = Q_\ell$ on `IndexCalcStrategy::toy()` |

---

## 18. Cross-References

### Mathematical textbook

- **`docs/MATHEMATICS.md` ch. 10** (algebraic ECDLP attacks) — the maths-first sibling to this
  integrative chapter. Develops the five algebraic ECDLP attacks mathematically:
  Pohlig–Hellman (§10.1), Smart–Satoh–Araki (§10.2), GHS/Weil descent (§10.3), index calculus (§10.4), and the
  **MOV/Frey–Rück reduction (§10.5 — the designated payoff proof)**. Carries the per-attack
  L-notation comparison (§10.6). Cross-references this chapter for the code realisation.

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
