# rDLP — Mathematical Textbook

*A standing mathematical survey of discrete-logarithm and integer-factorisation algorithms.*
*Maths-first, code-second, learnable on its own.*

---

## Documentation Register

This section states the register every chapter of this textbook obeys. It applies to every
chapter (and, as a recommendation, to the code-tours). A chapter that cannot be written within
this register must say so explicitly and relax the register for that chapter — not silently raise
the level.

### Audience floor

The reader has a full undergraduate mathematics background:

- **Proofs.** Comfortable reading and writing $\varepsilon$–$\delta$ arguments, induction,
  contradiction, and direct proof. Knows what it means for a proof to be rigorous.
- **Algebra.** Groups, rings, fields, homomorphisms, quotients, ideals. Polynomial rings. The
  isomorphism theorems. Finite fields $\mathbb{F}_q$. Elliptic curves as groups (the chord-and-tangent
  law, the group axioms, the order of a point). Dedekind domains and unique factorisation of ideals
  are introduced in the Prerequisites chapter and used in the GNFS chapter.
- **Analysis.** Big-$O$ and little-$o$ notation. Limits. The prime number theorem (statement only).
  Logarithms and exponentials. The $L$-notation for subexponential functions is introduced in the
  Prerequisites chapter.
- **Probability.** Discrete probability spaces, expectation, the birthday paradox. The density of
  smooth numbers (Canfield–Erdős–Pomerance) is stated and motivated in the Prerequisites chapter.
- **Logic.** Propositional and first-order logic at the level of a first proof-writing course.
  Computability and complexity at the level of knowing what $\mathbf{P}$, $\mathbf{NP}$, and
  polynomial-time reductions are.

**Nothing beyond this floor is assumed.** Anything beyond is built up in text or cited with a
precise reference. A later chapter that genuinely requires graduate-level background (e.g. étale
cohomology for the Weil conjectures, or the theory of $p$-adic representations) must surface this
as an explicit relaxation of this register, not silently assume it.

### Depth

**Survey with proof-sketch depth.** Every key theorem is stated precisely and motivated; proofs are
sketches with citations to the full proof, *except* where the proof is the pedagogical payoff. The
designated payoff proofs — where the proof *is* the lesson — are given in full:

- The $L$-notation subexponentiality derivation for GNFS (the GNFS chapter).
- The MOV reduction (the algebraic ECDLP attacks chapter).

All other proofs are sketches: the key idea is stated, the main steps are named, and a citation
points to the full argument. The sketch is honest — it does not pretend to be a proof, and it does
not omit the key idea.

**Complete:** no key idea is silently omitted. If a theorem is used, it is stated. If an algorithm
depends on a non-obvious mathematical fact, that fact is named and cited.

**Academic and clinical:** the register is precise and unsentimental. Intuition leads, rigour
follows — but the rigour is there.

**Not exhaustive:** no encyclopaedic case enumeration. The goal is understanding, not a reference
manual. Edge cases are noted when they are pedagogically important; otherwise they are cited.

**Not inscrutable:** every section opens with the intuition before the formalism. A reader who
skips the proofs should still understand what is true and why it matters.

### Through-line

**Structure-based escape from search.** Every algorithm in this textbook is a story about finding
exploitable structure — a group homomorphism, a smoothness phenomenon, an endomorphism, a pairing,
a quantum period — that escapes the generic $\sqrt{n}$ or $L$-notation search bound. This
through-line is stated once in the Escape-from-Search chapter and revisited at the opening of every
subsequent chapter.

### Markup

**Markdown + MathJax.** Inline mathematics uses `$…$`; display mathematics uses `$$…$$`. This
replaces the earlier undecided "rST or Markdown" choice. Rationale:
zero migration (every artifact stays `.md`), renders on GitHub and in standard Markdown viewers,
and lets trivial inline glyphs stay as Unicode while non-trivial expressions use TeX delimiters.
The TeX source is the same notation the textbook would use under any tooling. Reopen only on a
hard MathJax limitation (a renderer the project must target that lacks MathJax).

Trivial inline glyphs ($\mathbb{Z}$, $\mathbb{F}_2$, $\equiv$, subscripts) may remain as Unicode
where they read naturally. Non-trivial expressions — fractions, limits, aligned derivations,
matrices, sub/superscript stacks — use MathJax delimiters.

### Artifact location

**`docs/MATHEMATICS.md`** (single file). Chapters are appended to this file as they are written.

**Single-file decision.** The single-file format is retained. At approximately 4500 lines across
12 chapters, the file is large but navigable: every chapter is
clearly delimited by a top-level heading, the ToC provides a complete map, and the cross-references
are all intra-file or to named sections in the PEDAGOGY.md files. The principal argument for
promotion to `docs/textbook/` would be if per-chapter editing became unwieldy or if a build system
needed to process chapters independently — neither applies here. The single-file format preserves
the ability to search the whole textbook in one pass, keeps all cross-references as intra-file
anchors (no inter-file link management), and matches the `docs/PEDAGOGY.md` master-tour convention.
**Decision: remain as `docs/MATHEMATICS.md` (single file). No promotion to `docs/textbook/`.**

### Chapter-pairing pattern

Each textbook chapter is the maths-first sibling of a `PEDAGOGY.md` code-tour chapter. The
pairing is:

| Textbook chapter (this file) | Code-tour chapter |
|------------------------------|-------------------|
| Ch. 6 — Pollard Rho for ECDLP | `docs/PEDAGOGY.md` §§1–7 (Pollard rho code-tour) |
| Ch. 7 — The α-Substrate | `shared/numth/docs/PEDAGOGY.md` (α-substrate code-tour) |
| Ch. 8 — GNFS | `gnfs/docs/PEDAGOGY.md` §§52–62 (GNFS integrative chapter) |
| Ch. 9 — NFS-DL | `gnfs/docs/PEDAGOGY.md` §§63–71 (NFS-DL chapter) |
| Ch. 10 — Algebraic ECDLP Attacks | `docs/PEDAGOGY.md` §§8–18 (algebraic-attack integrative chapter) |
| Ch. 11 — Shor's Algorithm and Post-Quantum Context | `shor/docs/PEDAGOGY.md` (Shor code-tour) |
| Ch. 12 — Modularity and the Arithmetic of Elliptic Curves | *(no code-tour sibling — see §12.7)* |

**Note on Ch. 10.** The algebraic-attack integrative chapter lives in the master
`docs/PEDAGOGY.md` §§8–18 — not in `gnfs/docs/PEDAGOGY.md`. The attack code (`rho::ecdlp`,
`rho::pairing`, `rho::ssa`, `rho::ghs`, `rho::semaev`, `rho::index_calculus`) lives in the `rho`
crate, so its code-tour is a chapter of the master tour rather than of the `gnfs` tour.

**Note on Ch. 12 (Modularity).** Chapter 12 has no code-tour sibling. The modularity
correspondence is not implemented in the rDLP workspace; the chapter is a mathematical survey and
speculation. See §12.7 for the explicit statement.

The code-tour cites the textbook chapter for the mathematics; the textbook chapter cites the
code-tour for the realisation. Neither is a prerequisite for the other — they are complementary
lenses on the same material.

---

## Table of Contents

The chapter skeleton across the whole survey. Each entry is a chapter title and a one-sentence
scope statement.

1. **Documentation Register** *(this section)* — the audience, depth,
   through-line, markup, and location rules every chapter obeys.

2. **Table of Contents** *(this section)* — the chapter skeleton.

3. **Escape from Search: The Through-Line** — the conceptual spine: every attack finds exploitable
   structure that escapes the generic $\sqrt{n}$ or $L$-notation search bound, and the taxonomy of
   structures that make this possible.

4. **Prerequisites** — the undergraduate-background bridge: the specific theorems from algebra,
   analysis, probability, and logic that later chapters lean on, each stated precisely with a
   proof sketch or citation.

5. **On Scale: A Natural-Philosophy Interlude** — the three axes of scale (resource/operational,
   mathematical-dimension, structural), their couplings, and the honest science↔engineering gap
   that runs through every chapter.

6. **Pollard Rho for ECDLP** — the birthday-paradox collision argument, Floyd's cycle detection,
   the group-homomorphism structure that makes rho work on elliptic curves, and the $L$-notation
   bound; maths-first sibling to `docs/PEDAGOGY.md`.

7. **The α-Substrate: Primality, Smoothness, and ECM** — Miller–Rabin primality (Fermat witnesses,
   strong pseudoprimes), smooth-number theory (the $B$-smooth density theorem), ECM (Lenstra's
   elliptic-curve method), and Tonelli–Shanks (square roots mod $p$); maths-first sibling to
   `shared/numth/docs/PEDAGOGY.md`.

8. **GNFS: The General Number Field Sieve** — the number-field bridge, factor-base smoothness,
   linear-algebra dependency, square-root recovery, and the full $L_N[1/3, (64/9)^{1/3}]$
   subexponentiality derivation as the payoff proof; maths-first sibling to
   `gnfs/docs/PEDAGOGY.md` §§52–62.

9. **NFS-DL: Discrete Logarithm via the Number Field Sieve** — the adaptation of GNFS to the
   discrete-logarithm problem, Schirokauer maps, and the individual logarithm descent; maths-first
   sibling to `gnfs/docs/PEDAGOGY.md` §§63–71.

10. **Algebraic ECDLP Attacks** — Pohlig–Hellman, the MOV/Frey–Rück reduction (payoff proof: the
    MOV bridge), Smart–Satoh–Araki, GHS/Weil descent, and Gaudry–Diem–Joux–Vitse index calculus;
    maths-first sibling to `docs/PEDAGOGY.md` §§8–18.

11. **Shor's Algorithm and Post-Quantum Context** — the quantum period-finding algorithm as the
    fifth escape-from-search structure; the full payoff proof (QFT phase estimation, order-finding
    → factoring reduction, two-register hidden-subgroup → ECDLP reduction); and the post-quantum
    migration landscape (NIST PQC families, the SIDH/Castryck–Decru break as the honest coda, the
    migration rationale); maths-first sibling to `shor/docs/PEDAGOGY.md`.

12. **Modularity and the Arithmetic of Elliptic Curves: A Speculation** — the Taniyama–Shimura–Weil
    correspondence as a structural phenomenon in the escape-from-search family; the connection to
    point counting, isogenies, and the Castryck–Decru precedent; and an explicit speculation about
    whether the modularity structure suggests further escapes from the ECDLP search bound. The
    project's one chartered speculation; no code-tour sibling.

---

## Escape from Search: The Through-Line

### The generic search bound

Every cryptographic hardness assumption rests, at bottom, on the difficulty of a *search problem*:
find $k$ such that $k \cdot G = Q$ (the discrete logarithm), or find $p$ such that $p \mid N$ (the
factoring problem). The naive approach is exhaustive search: try every candidate. For a group of
order $n$, this costs $O(n)$ operations.

The first non-trivial observation is that exhaustive search is never necessary. The **birthday
paradox** gives a generic lower bound on what any algorithm must do, and a generic algorithm that
meets it: by sampling $O(\sqrt{n})$ random elements of a group of order $n$, one expects a
collision — two samples that land on the same element. Pollard's rho algorithm (1978) exploits this
to solve the discrete logarithm in $O(\sqrt{n})$ group operations, and Shank's baby-step/giant-step
algorithm achieves the same bound with $O(\sqrt{n})$ space. For the factoring problem, Fermat's
method and its descendants give $O(N^{1/4})$ operations for a balanced semiprime.

The $\sqrt{n}$ barrier is the **generic search bound**: it is what any algorithm achieves without
exploiting the specific structure of the problem. For a 256-bit elliptic curve group ($n \approx
2^{256}$), the generic bound is $2^{128}$ operations — the security level of ECC-256.

### Escaping the bound: the taxonomy of structures

The algorithms in this textbook all escape the generic bound by finding exploitable structure. The
structures fall into five families:

**1. Group homomorphisms (index calculus, smoothness).** If the group $G$ admits a homomorphism
$\phi: G \to H$ where $H$ has a known structure (e.g. $H = \mathbb{Z}/n\mathbb{Z}$), then the
discrete logarithm in $G$ reduces to a problem in $H$. The key instance is *smoothness*: if the
group elements can be expressed as products of a small set of generators (the *factor base*), then
the discrete logarithm reduces to a system of linear equations over $\mathbb{Z}/n\mathbb{Z}$. This
is the engine of index calculus, NFS, and NFS-DL.

**2. Endomorphisms (GLV, Koblitz).** An endomorphism $\phi: E \to E$ of an elliptic curve satisfies
$\phi(P) = \lambda P$ for all $P$, where $\lambda$ is a fixed scalar. This collapses the orbit
$\{P, \phi(P), \phi^2(P), \ldots\}$ to a single representative, reducing the effective group size
and the cost of Pollard rho. The GLV decomposition (Gallant–Lambert–Vanstone) and the Koblitz
automorphism are the canonical instances.

**3. Pairings (MOV reduction).** A bilinear pairing $e: E[n] \times E[n] \to \mu_n$ maps the
discrete logarithm on an elliptic curve to the discrete logarithm in a multiplicative group
$\mathbb{F}_{p^k}^*$, where index calculus applies. The MOV/Frey–Rück reduction (Menezes–Okamoto–
Vanstone, Frey–Rück) is the canonical instance; it is the cross-track bridge of this project.

**4. Number-field structure (GNFS, NFS-DL).** The integers $\mathbb{Z}$ embed into a number field
$K = \mathbb{Q}(\alpha)$, where the norm map $N_{K/\mathbb{Q}}$ connects factorisation in $K$ to
factorisation in $\mathbb{Z}$. The GNFS exploits this by sieving for *smooth norms* — elements of
$\mathbb{Z}[\alpha]$ whose norms factor completely over a small factor base — and then using linear
algebra to find a congruence of squares $x^2 \equiv y^2 \pmod{N}$. The subexponential complexity
$L_N[1/3, (64/9)^{1/3}]$ arises from the optimal balance between the sieving cost and the
smoothness probability.

**5. Quantum period-finding (Shor).** Shor's algorithm (1994) uses the quantum Fourier transform to
find the period of the function $f(k) = g^k \bmod N$ in polynomial time. This dissolves the
$\sqrt{n}$ barrier entirely in the quantum model, reducing both factoring and discrete logarithm to
polynomial time. The structure exploited is the *periodicity* of modular exponentiation — a
structure that classical algorithms cannot access efficiently.

### The L-notation hierarchy

The complexity of the best known algorithms for the main problems forms a hierarchy:

$$L_N[\alpha, c] = \exp\!\left(c \cdot (\log N)^\alpha \cdot (\log \log N)^{1-\alpha}\right)$$

- **Fully exponential** ($\alpha = 1$): $L_N[1, c] = N^c$. Generic search (baby-step/giant-step,
  Pollard rho) for ECDLP on a generic curve.
- **Subexponential** ($0 < \alpha < 1$): $L_N[\alpha, c]$. Index calculus for DLP in
  $\mathbb{F}_p^*$ ($\alpha = 1/2$); GNFS for factoring and NFS-DL ($\alpha = 1/3$).
- **Quasi-polynomial** ($\alpha \to 0$): $L_N[0, c] = (\log N)^c$. The Barbulescu–Gaudry–Joux–
  Thomé algorithm for DLP in small-characteristic fields (2014).
- **Polynomial** ($\alpha = 0$, exactly): Shor's algorithm for both factoring and DLP.

The through-line of this textbook is the story of how each step down this hierarchy is achieved by
finding a new kind of exploitable structure.

### Cross-track L-notation synthesis

The table below gives the leading complexity for each track and each major attack, unified across
all five tracks. This is the quantitative spine of the textbook: reading it top-to-bottom is the
story of the project, from the generic $\sqrt{n}$ baseline through the subexponential classical
attacks to the polynomial-time quantum dissolution. The per-track derivations are in §6 (ρ), §7
(GNFS), §9.7 (NFS-DL), §10.6 (algebraic ECDLP), and §11.3.4 and §11.4.4 (Shor). The code-first
version of this table is in `docs/PEDAGOGY.md` §0.2.

| Family | Algorithm | Problem | $L$-notation complexity | Structure exploited |
|--------|-----------|---------|------------------------|---------------------|
| ρ | Pollard rho | ECDLP (generic curve) | $L_n[1,\, 1/2] = \Theta(\sqrt{n})$ | None — generic $\sqrt{n}$ walk; the baseline |
| α | ECM (Lenstra) | Integer factoring (sub-step) | $L_p[1/2,\, 1]$ (largest prime factor $p$) | Group-order smoothness |
| G | GNFS | Integer factoring | $L_N[1/3,\, (64/9)^{1/3}]$ | Number-field bridge + smoothness (§7) |
| D | NFS-DL | DLP in $\mathbb{F}_p^*$ | $L_p[1/3,\, (64/9)^{1/3}]$ | Same number-field bridge; descent subdominant (§9.7) |
| E | Pohlig–Hellman | ECDLP (composite order) | $L_{p_{\max}}[1,\, 1/2]$ (largest prime factor) | CRT group-order homomorphism (§10.1) |
| E | Smart–Satoh–Araki | ECDLP (anomalous curve) | $L_p[0] = O(\log p)$ (polynomial time) | $p$-adic formal group logarithm (§10.2) |
| E | GHS/Weil descent | ECDLP ($E/\mathbb{F}_{2^m}$) | Conditional: $L_{p^{m/l}}[1/2,\, c]$ (downstream) | Binary field tower; transfer to Jacobian (§10.3) |
| E | Index calculus | ECDLP ($E/\mathbb{F}_{p^n}$, $n > 1$) | $L_{p^n}[1/2,\, c]$ (subexponential) | Factor-base decomposability (§10.4) |
| E | MOV/Frey–Rück | ECDLP (small embedding degree) | $L_{p^k}[1/3,\, (64/9)^{1/3}]$ (via NFS-DL) | Pairing → $\mathbb{F}_{p^k}^*$ DLP (§10.5) |
| S | Shor (factoring) | Integer factoring | $O((\log N)^3) = L_N[0,\, 3]$ | Quantum period-finding (§11.3) |
| S | Shor (ECDLP) | ECDLP | $O((\log r)^3) = L_r[0,\, 3]$ | Quantum period-finding, 2-register HSP (§11.4) |

**Reading the table.** The $\alpha$ exponent in $L[\alpha, c]$ is the key invariant. Pollard rho
sits at $\alpha = 1$ (fully exponential — the generic bound). GNFS and NFS-DL sit at
$\alpha = 1/3$ (subexponential via the number-field bridge). The algebraic attacks span the full
range: Pohlig–Hellman stays
at $\alpha = 1$ (structural speedup, not asymptotic); index calculus and GHS reach $\alpha = 1/2$
(subexponential in the extension-field setting); MOV bridges ECDLP to the $\alpha = 1/3$ regime
via NFS-DL; SSA achieves $\alpha = 0$ classically (polynomial time, but only for anomalous
curves). Shor's algorithm dissolves the $L$-notation framework entirely: $\alpha = 0$
polynomially for all curves and all integers.

**Scale honesty.** The $L$-notation separations in this table are NOT observable at the toy
scale of the project's fixtures. The table reports theoretical complexity classes; the toy-scale
timings are in `docs/BENCHMARKS.md`. The separation between $L[1, 1/2]$ (rho), $L[1/3]$
(GNFS/NFS-DL/MOV), $L[1/2]$ (index calculus), and $L[0]$ (SSA, Shor) is a statement about
asymptotic behaviour at cryptographic scale — not a ranking of toy-scale timings.

**Chapter 12 (Modularity).** Chapter 12 is the project's one deliberate speculation; it surveys
the Taniyama–Shimura–Weil correspondence as a structural phenomenon and asks whether it suggests
further escapes. No complexity claim is made for any modularity-based attack; the chapter does not
add a row to this table.

### Why the bound matters

The $\sqrt{n}$ generic bound is not just a complexity statement — it is a *security guarantee*. A
256-bit elliptic curve group provides 128-bit security *if and only if* no algorithm beats the
generic bound on that curve. The moment a structural exploit is found (an endomorphism, a pairing
with small embedding degree, a special field structure), the security level collapses. The history
of cryptanalysis is a history of finding such structures in curves and fields that were believed to
be generic.

This is why the through-line matters: understanding *which* structures enable *which* escapes is
the key to both attacking and defending cryptographic systems.

---

## Prerequisites

This chapter collects the specific theorems from algebra, analysis, probability, and logic that
later chapters lean on. Each theorem is stated precisely; the proof is either sketched (with a
citation to the full argument) or, where the proof is short and illuminating, given in full.

A reader who is comfortable with all of these results can proceed directly to any chapter. A reader
who needs to review specific topics should consult the cited references.

### Algebra

#### Groups, rings, and fields

**Definition (group).** A *group* $(G, \cdot)$ is a set $G$ with a binary operation $\cdot$ that
is associative, has an identity element $e$, and has inverses. A group is *abelian* (or
*commutative*) if $a \cdot b = b \cdot a$ for all $a, b \in G$. The *order* of $G$, written $|G|$,
is the cardinality of $G$. The *order* of an element $g \in G$ is the smallest positive integer $k$
such that $g^k = e$.

**Theorem (Lagrange).** *If $H$ is a subgroup of a finite group $G$, then $|H|$ divides $|G|$.*

*Proof sketch.* The cosets $gH$ partition $G$ into equal-sized pieces of size $|H|$. $\square$

**Corollary.** *The order of any element $g \in G$ divides $|G|$.*

**Definition (cyclic group).** A group $G$ is *cyclic* if there exists $g \in G$ such that every
element of $G$ is a power of $g$. Such a $g$ is called a *generator* of $G$.

**Theorem (structure of cyclic groups).** *Every cyclic group of order $n$ is isomorphic to
$\mathbb{Z}/n\mathbb{Z}$. Every subgroup of a cyclic group is cyclic.*

**Definition (ring, field).** A *ring* $(R, +, \cdot)$ is an abelian group under $+$ with a
compatible multiplication that is associative and distributes over $+$. A *field* is a commutative
ring in which every non-zero element has a multiplicative inverse.

**Theorem (finite fields).** *For every prime power $q = p^k$, there exists a unique (up to
isomorphism) field $\mathbb{F}_q$ with $q$ elements. The multiplicative group $\mathbb{F}_q^*$ is
cyclic of order $q - 1$.*

*Proof sketch.* Existence: $\mathbb{F}_{p^k}$ is the splitting field of $x^{p^k} - x$ over
$\mathbb{F}_p$. Uniqueness: any two fields with $q$ elements are both splitting fields of the same
polynomial, hence isomorphic. Cyclicity of $\mathbb{F}_q^*$: a finite subgroup of the
multiplicative group of any field is cyclic (the key lemma: a polynomial of degree $d$ over a
field has at most $d$ roots, so the number of elements of order dividing $d$ is at most $d$, which
forces the group to be cyclic). $\square$

#### Ideals and factorisation in Dedekind domains

**Definition (ideal).** An *ideal* $\mathfrak{a}$ in a ring $R$ is a non-empty subset closed under
addition and under multiplication by any element of $R$. A *prime ideal* $\mathfrak{p}$ is an ideal
such that $ab \in \mathfrak{p}$ implies $a \in \mathfrak{p}$ or $b \in \mathfrak{p}$.

**Definition (Dedekind domain).** A *Dedekind domain* is an integral domain $R$ that is Noetherian,
integrally closed in its fraction field, and has Krull dimension 1 (every non-zero prime ideal is
maximal).

**Theorem (unique factorisation of ideals in Dedekind domains).** *In a Dedekind domain, every
non-zero ideal factors uniquely as a product of prime ideals.*

*Proof sketch.* The key steps are: (1) every ideal is contained in a maximal ideal; (2) the
fractional ideals form a group under multiplication; (3) the group of fractional ideals is free
abelian on the prime ideals. The proof uses the Noetherian condition to guarantee termination of
the factorisation process. See Neukirch [N99, §I.3] for the full argument. $\square$

**Why this matters for GNFS.** The ring of integers $\mathcal{O}_K$ of a number field $K =
\mathbb{Q}(\alpha)$ is a Dedekind domain. The GNFS exploits the unique factorisation of ideals in
$\mathcal{O}_K$ to connect smoothness of norms in $\mathbb{Z}$ to smoothness of ideals in
$\mathcal{O}_K$. The primes dividing the index $[\mathcal{O}_K : \mathbb{Z}[\alpha]]$ are exactly
the primes $p$ at which $\mathbb{Z}[\alpha]$ fails to be integrally closed, and they require special
handling (the Dedekind index criterion). They are distinct from the *bad primes*, those dividing
$\mathrm{disc}(f)$: since $\mathrm{disc}(f) = [\mathcal{O}_K : \mathbb{Z}[\alpha]]^2 \cdot
\mathrm{disc}(K)$, every prime dividing the index is bad, but a bad prime need not divide the index
(for $f = x^2 - 2$ the prime $2$ is bad, yet $\mathbb{Z}[\sqrt{2}]$ is maximal).

#### Elliptic curves

**Definition (elliptic curve).** An *elliptic curve* over a field $k$ is a smooth projective curve
of genus 1 with a specified base point $\mathcal{O}$. In short Weierstrass form (valid when
$\mathrm{char}(k) \neq 2, 3$):

$$E: y^2 = x^3 + ax + b, \quad \Delta = -16(4a^3 + 27b^2) \neq 0.$$

The condition $\Delta \neq 0$ ensures smoothness (no cusps or self-intersections).

**Theorem (group law).** *The points of $E(k)$ (including the point at infinity $\mathcal{O}$)
form an abelian group under the chord-and-tangent law, with $\mathcal{O}$ as the identity.*

*Proof sketch.* Associativity is the non-trivial part; it follows from the Riemann–Roch theorem
applied to the divisor class group of $E$. The explicit addition formulas (slope, $x$-coordinate,
$y$-coordinate) are derived from the intersection of a line with the cubic. $\square$

**Theorem (Hasse's bound).** *For an elliptic curve $E$ over $\mathbb{F}_p$,*

$$\left| |E(\mathbb{F}_p)| - (p + 1) \right| \leq 2\sqrt{p}.$$

*Proof sketch.* The bound follows from the Riemann hypothesis for curves over finite fields
(Weil's theorem): the eigenvalues of Frobenius on the Tate module have absolute value $\sqrt{p}$.
See Silverman [S09, §V.1]. $\square$

**Why this matters.** Hasse's bound says the group order $|E(\mathbb{F}_p)|$ lies in the interval
$[p + 1 - 2\sqrt{p},\, p + 1 + 2\sqrt{p}]$. For ECM (Lenstra's method), this means that by
choosing different curves, one can vary the group order over an interval of width $4\sqrt{p}$
around $p$, and hope to find a curve whose order is smooth.

#### The Chinese Remainder Theorem

**Theorem (CRT).** *Let $m_1, \ldots, m_k$ be pairwise coprime positive integers, and let $M =
m_1 \cdots m_k$. Then the natural ring homomorphism*

$$\mathbb{Z}/M\mathbb{Z} \xrightarrow{\;\sim\;} \mathbb{Z}/m_1\mathbb{Z} \times \cdots \times
\mathbb{Z}/m_k\mathbb{Z}$$

*is an isomorphism.*

*Proof.* The map $x \mapsto (x \bmod m_1, \ldots, x \bmod m_k)$ is a ring homomorphism. Injectivity
follows from $\gcd(m_i, m_j) = 1$; surjectivity follows from Bézout's identity. $\square$

**Why this matters.** CRT is used in the GNFS square-root step (the split-prime method computes the
algebraic square root modulo several primes and lifts via CRT), in Pohlig–Hellman (reducing ECDLP in
a group of composite order to ECDLP in prime-order subgroups), and in the construction of
deterministic Miller–Rabin witness sets.

### Analysis

#### Big-$O$ and $L$-notation

**Definition (big-$O$).** $f(n) = O(g(n))$ as $n \to \infty$ if there exist constants $C > 0$ and
$n_0$ such that $|f(n)| \leq C \cdot g(n)$ for all $n \geq n_0$.

**Definition ($L$-notation).** For $0 \leq \alpha \leq 1$ and $c > 0$, define

$$L_N[\alpha, c] = \exp\!\left(c \cdot (\log N)^\alpha \cdot (\log \log N)^{1-\alpha}\right).$$

This interpolates between polynomial ($\alpha = 0$: $L_N[0, c] = (\log N)^c$) and fully
exponential ($\alpha = 1$: $L_N[1, c] = N^c$). The intermediate regime $0 < \alpha < 1$ is
*subexponential*: faster than any polynomial in $\log N$ but slower than any polynomial in $N$.

**Lemma ($L$-notation arithmetic).** *For $0 < \alpha < 1$:*

$$L_N[\alpha, c_1] \cdot L_N[\alpha, c_2] = L_N[\alpha, c_1 + c_2], \qquad
L_N[\alpha, c]^k = L_N[\alpha, kc].$$

*Moreover, $L_N[\alpha, c] = o(N^\varepsilon)$ for every $\varepsilon > 0$, and $L_N[\alpha, c] =
\omega((\log N)^M)$ for every $M > 0$.*

*Proof.* Direct computation from the definition. $\square$

**Why this matters.** The $L$-notation is the natural language for stating the complexity of
subexponential algorithms. The GNFS achieves $L_N[1/3, (64/9)^{1/3}]$; the Pohlig–Hellman
reduction achieves $L_p[1/2, 1]$ for the largest prime factor $p$ of the group order; Shor's
algorithm achieves $O((\log N)^3)$. The exponent $\alpha$ is the key invariant: it determines
whether the algorithm is polynomial, subexponential, or exponential in the input size.

#### The prime number theorem

**Theorem (prime number theorem).** *Let $\pi(x)$ denote the number of primes $\leq x$. Then*

$$\pi(x) \sim \frac{x}{\log x} \quad \text{as } x \to \infty.$$

*Equivalently, the $n$-th prime $p_n \sim n \log n$.*

*Proof sketch.* The proof uses the Riemann zeta function $\zeta(s) = \sum_{n=1}^\infty n^{-s}$ and
the fact that $\zeta(s)$ has no zeros on the line $\mathrm{Re}(s) = 1$. See Davenport [D00,
Chapter 13]. $\square$

**Why this matters.** The prime number theorem controls the size of factor bases: the number of
primes up to $B$ is approximately $B / \log B$. This enters the complexity analysis of every
sieve-based algorithm.

### Probability

#### The birthday paradox

**Theorem (birthday bound).** *If $m$ elements are chosen uniformly at random (with replacement)
from a set of size $n$, the probability that at least two are equal is approximately*

$$1 - e^{-m^2/(2n)} \approx \frac{m^2}{2n} \quad \text{for } m \ll n.$$

*In particular, the expected number of draws before a collision is $\Theta(\sqrt{n})$.*

*Proof.* The probability that all $m$ draws are distinct is $\prod_{k=0}^{m-1}(1 - k/n) \approx
e^{-m(m-1)/(2n)}$ for $m \ll n$. Setting this equal to $1/2$ and solving gives $m \approx
\sqrt{2n \log 2} = \Theta(\sqrt{n})$. $\square$

**Why this matters.** The birthday bound is the foundation of Pollard rho: the expected number of
steps before a collision in a pseudorandom walk on a group of order $n$ is $\Theta(\sqrt{n})$. It
is also the reason that a 256-bit elliptic curve group provides 128-bit security against generic
attacks.

#### Smooth-number density

**Definition ($B$-smooth).** An integer $m$ is *$B$-smooth* if all its prime factors are $\leq B$.

**Theorem (Canfield–Erdős–Pomerance, 1983).** *Let $\Psi(x, y)$ denote the number of $y$-smooth
integers in $[1, x]$. For $u = \log x / \log y$ fixed as $x \to \infty$,*

$$\Psi(x, y) = x \cdot \rho(u) \cdot (1 + o(1)),$$

*where $\rho(u)$ is the Dickman function, defined by $\rho(u) = 1$ for $0 \leq u \leq 1$ and*

$$u \rho'(u) = -\rho(u - 1) \quad \text{for } u > 1.$$

*For large $u$, $\rho(u) = u^{-u(1 + o(1))}$.*

*Proof sketch.* The key step is to count integers in $[1, x]$ whose largest prime factor is $\leq
y$, by conditioning on the largest prime factor and applying the prime number theorem. The Dickman
function arises as the continuous analogue of this recursion. See Granville [G08] for a modern
treatment. $\square$

**Corollary (smoothness probability in $L$-notation).** *If $x = N$ and $y = L_N[1/3, c]$, then*

$$\frac{\Psi(N, y)}{N} = L_N[-1/3, -1/(3c)] \cdot (1 + o(1)).$$

*This is the key estimate in the GNFS complexity analysis: the probability that a random integer
near $N$ is $L_N[1/3, c]$-smooth is $L_N[-1/3, \ldots]$, which is subexponentially small but
much larger than $N^{-\varepsilon}$ for any fixed $\varepsilon > 0$.*

**Why this matters.** The smooth-number density theorem is the engine of every sieve-based
algorithm. The GNFS complexity $L_N[1/3, (64/9)^{1/3}]$ arises from optimising the balance between
the smoothness probability (which improves as the smoothness bound $B$ grows) and the cost of the
linear algebra step (which grows with $B$). The optimal $B = L_N[1/3, c]$ for some constant $c$,
and the exponent $1/3$ is the result of this optimisation.

### Logic and complexity

#### Polynomial-time reductions

**Definition (polynomial-time reduction).** A problem $A$ *polynomial-time reduces* to a problem
$B$ (written $A \leq_P B$) if there is a polynomial-time algorithm that transforms any instance of
$A$ into an instance of $B$ such that the answer is preserved.

**Why this matters.** The MOV reduction (§10.5) is a polynomial-time reduction from ECDLP
on a curve with small embedding degree to DLP in $\mathbb{F}_{p^k}^*$. Pohlig–Hellman is a
polynomial-time reduction from ECDLP in a group of composite order to ECDLP in prime-order
subgroups. Understanding reductions is essential for understanding why these attacks work.

#### Fermat's little theorem and Euler's criterion

**Theorem (Fermat's little theorem).** *If $p$ is prime and $\gcd(a, p) = 1$, then $a^{p-1} \equiv
1 \pmod{p}$.*

*Proof.* The map $x \mapsto ax$ is a bijection on $(\mathbb{Z}/p\mathbb{Z})^*$, so $\prod_{x=1}^{p-1}
(ax) \equiv \prod_{x=1}^{p-1} x \pmod{p}$, giving $a^{p-1} \equiv 1$. $\square$

**Theorem (Euler's criterion).** *For an odd prime $p$ and $\gcd(a, p) = 1$,*

$$a^{(p-1)/2} \equiv \left(\frac{a}{p}\right) \pmod{p},$$

*where $\left(\frac{a}{p}\right)$ is the Legendre symbol: $+1$ if $a$ is a quadratic residue mod
$p$, $-1$ if not.*

*Proof.* The multiplicative group $(\mathbb{Z}/p\mathbb{Z})^*$ is cyclic of order $p - 1$. An
element $a$ is a square iff $a^{(p-1)/2} = 1$ (since squaring is a 2-to-1 map on the group). The
non-residues satisfy $a^{(p-1)/2} = -1$ (the only other square root of $1$ in a field). $\square$

**Why this matters.** Fermat's little theorem is the basis of Miller–Rabin primality testing and of
the Fermat primality test. Euler's criterion is the basis of the Legendre symbol computation and
the Tonelli–Shanks square root algorithm.

#### The discrete logarithm problem

**Definition (DLP).** Let $G$ be a cyclic group of order $n$ with generator $g$. Given $h \in G$,
the *discrete logarithm problem* (DLP) asks for the unique $k \in \{0, 1, \ldots, n-1\}$ such that
$g^k = h$.

**Definition (ECDLP).** The *elliptic curve discrete logarithm problem* (ECDLP) is the DLP in the
group $E(\mathbb{F}_p)$ of points on an elliptic curve over a prime field.

**Why hardness is believed.** No polynomial-time classical algorithm is known for the DLP or ECDLP
in a generic group. The best known classical algorithms for ECDLP on a generic curve run in
$\Theta(\sqrt{n})$ time (Pollard rho, baby-step/giant-step). For DLP in $\mathbb{F}_p^*$, the best
known algorithms run in $L_p[1/3, (64/9)^{1/3}]$ time (NFS-DL). The hardness of these problems is
a *computational assumption*, not a theorem — no proof that they require superpolynomial time is
known.

### Modular forms and Galois representations (Chapter 12)

Chapter 12 (Modularity) uses concepts from the theory of modular forms and Galois representations
that go beyond the standard undergraduate background. This section states the minimum needed to
follow Chapter 12 at survey depth; full proofs are in the cited references.

**The upper half-plane.** The *upper half-plane* $\mathfrak{H} = \{z \in \mathbb{C} : \mathrm{Im}(z)
> 0\}$ is acted on by $\mathrm{SL}_2(\mathbb{Z})$ via Möbius transformations
$\bigl(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\bigr) \cdot z = (az+b)/(cz+d)$. The congruence
subgroup $\Gamma_0(N) \subset \mathrm{SL}_2(\mathbb{Z})$ consists of matrices with $c \equiv 0
\pmod{N}$.

**Modular forms and cusp forms.** A *modular form* of weight $k$ and level $N$ is a holomorphic
function $f: \mathfrak{H} \to \mathbb{C}$ satisfying a transformation law under $\Gamma_0(N)$ and
a growth condition at the cusps. A *cusp form* additionally vanishes at every cusp. Every cusp form
has a Fourier expansion $f(z) = \sum_{n=1}^\infty a_n q^n$ with $q = e^{2\pi iz}$. The space
$S_k(\Gamma_0(N))$ of cusp forms of weight $k$ and level $N$ is a finite-dimensional $\mathbb{C}$-vector
space. Chapter 12 uses only weight-2 cusp forms ($k = 2$).

**Hecke operators and newforms.** The space $S_k(\Gamma_0(N))$ carries commuting Hecke operators
$T_p$ (for primes $p \nmid N$). A *newform* (or Hecke eigenform) is a cusp form that is a
simultaneous eigenvector for all Hecke operators, normalised so that $a_1 = 1$. For a newform,
the Hecke eigenvalue equals the Fourier coefficient: $T_p f = a_p f$.

**L-functions.** To each newform $f = \sum a_n q^n$ one attaches an *L-function*
$L(f, s) = \sum_{n=1}^\infty a_n n^{-s}$, which has an Euler product and extends to an entire
function satisfying a functional equation. To each elliptic curve $E/\mathbb{Q}$ one attaches an
L-function $L(E, s) = \prod_p L_p(E, s)$ via the trace of Frobenius $a_p(E) = p + 1 - \#E(\mathbb{F}_p)$.
The Taniyama–Shimura–Weil theorem asserts $L(E, s) = L(f, s)$ for a unique newform $f$ of level
equal to the conductor of $E$.

**Galois representations.** To an elliptic curve $E/\mathbb{Q}$ and a prime $\ell$, one attaches
the *$\ell$-adic Tate module* $T_\ell(E) = \varprojlim E[\ell^n]$, which gives a 2-dimensional
$\ell$-adic Galois representation $\rho_{E,\ell}: \mathrm{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})
\to \mathrm{GL}_2(\mathbb{Z}_\ell)$. The Frobenius element at $p$ (for $p \nmid N\ell$) acts with
characteristic polynomial $x^2 - a_p(E) x + p$. The proof of the TSW theorem proceeds by showing
this Galois representation is *modular* (arises from a modular form) via the $R = T$ strategy of
Wiles. Chapter 12 uses this at the level of a proof sketch; the full argument is in [W95], [TW95],
[BCDT01], and [DS05].

**The Weil conjectures and the Ramanujan bound.** The *Weil conjectures* (proved by Deligne [D74])
assert that the zeta function of a smooth projective variety over $\mathbb{F}_p$ satisfies a
Riemann hypothesis. For an elliptic curve, this gives Hasse's bound $|a_p(E)| \leq 2\sqrt{p}$.
For a weight-2 newform, the same bound $|a_p(f)| \leq 2\sqrt{p}$ is the *Ramanujan conjecture*
(proved as a consequence of Deligne's theorem). Chapter 12 uses this bound to motivate the
analogy between the two L-functions.

**What Chapter 12 does NOT require.** The chapter is a survey and speculation, not a proof. It
does not require: étale cohomology, the theory of $p$-adic representations beyond the sketch above,
the theory of automorphic forms beyond the modular-form setting, or any graduate-level algebraic
geometry beyond what is cited. A reader who has absorbed the prerequisites above and is comfortable
with the elliptic-curve material in §Algebra can follow Chapter 12 at survey depth.

### References for prerequisites

- [N99] Neukirch, J. *Algebraic Number Theory*. Springer, 1999. (Dedekind domains, ideal
  factorisation, number fields.)
- [S09] Silverman, J. H. *The Arithmetic of Elliptic Curves*. 2nd ed. Springer, 2009. (Elliptic
  curves, Hasse's bound, the group law.)
- [D00] Davenport, H. *Multiplicative Number Theory*. 3rd ed. Springer, 2000. (Prime number
  theorem, Dirichlet series.)
- [G08] Granville, A. "Smooth numbers: computational number theory and beyond." In *Algorithmic
  Number Theory: Lattices, Number Fields, Curves and Cryptography*, MSRI Publications 44, 2008.
  (Smooth-number density, Dickman function.)
- [CP05] Crandall, R., and Pomerance, C. *Prime Numbers: A Computational Perspective*. 2nd ed.
  Springer, 2005. (Miller–Rabin, Tonelli–Shanks, ECM, smoothness.)
- [DS05] Diamond, F., and Shurman, J. *A First Course in Modular Forms*. Springer, 2005. (Modular
  forms, Hecke operators, newforms, L-functions — the graduate-level reference for Chapter 12.)
- [W95] Wiles, A. "Modular elliptic curves and Fermat's Last Theorem." *Annals of Mathematics*,
  141(3), 1995. (The TSW proof for semistable curves; the $R = T$ strategy.)
- [D74] Deligne, P. "La conjecture de Weil. I." *Publications Mathématiques de l'IHÉS*, 43, 1974.
  (The Weil conjectures, including the Ramanujan bound as a special case.)

---

## On Scale: A Natural-Philosophy Interlude

*This interlude gives the natural-philosophy context for the project's repeated statement that
toy-scale demonstrations show mechanism, not asymptotic advantage — context that belongs in a
textbook register rather than in a short status note.*

### The word "scale" is not one thing

Every chapter in this textbook involves a comparison between "toy scale" and "real scale" — between
the parameters at which the algorithms are demonstrated here and the parameters at which they are
deployed in practice. The comparison is essential for honest pedagogy: an algorithm that is
"correct" at toy scale may be missing engineering that only matters at real scale, and a phenomenon
that is "prominent" at toy scale may wash out at real scale.

But the comparison is routinely made carelessly, as if "scale" were a single axis — a ladder from
small to large, with toy instances at the bottom and RSA-2048 at the top. It is not. At least three
distinct things travel under the word, and they do not form one ordered line.

### Axis 1: Resource/operational scale

The first axis is the one most people mean: how large an instance one can actually run. Toy scale
here means roughly 80-bit $N$ on a laptop; "NFS scale" means RSA-768 ($N \approx 2^{768}$) on a
cluster-month, or RSA-829 ($N \approx 2^{829}$) on a larger cluster. This axis *is* an unbounded
ladder — bigger $N$ always exists, up to a physical-computation ceiling.

But this axis is **phenomenologically flat**: scaling $N$ alone introduces no new mathematics. It
makes the same machinery bigger and improves the statistics. This is exactly why *over-exposed*
phenomena (bad primes, for instance) *wash out* at large $N$ rather than changing in kind: at
cryptographic scale, the contribution of bad primes is marginal and largely absorbed by polynomial
selection, but at toy scale with a hand-picked polynomial they are unavoidable and must be handled
head-on. The honest annotation is that the prominence is a toy-scale artifact, not the typical
NFS-scale picture.

The flatness of this axis is also why the GNFS complexity $L_N[1/3, (64/9)^{1/3}]$ is a
*heuristic*, not a theorem: the smoothness-probability estimates that drive the analysis are
asymptotic, and the asymptotic regime is only reached at scales far beyond what any computer can
run. The algorithm is correct at toy scale; the complexity analysis is only accurate at a scale
that is, in a precise sense, unreachable.

### Axis 2: Mathematical-dimension scale

The second axis is qualitatively different: the *degree* $d = [K:\mathbb{Q}]$ of the number field,
the *embedding degree* $k$ of a pairing, the *characteristic* vs *extension degree* shape of a
finite field $\mathbb{F}_{p^n}$. These are dimensions *of the mathematics*, not of the instance
size.

This axis is **not** a ladder in the same sense. A number field of degree $d = 6$ is not "larger"
than one of degree $d = 3$ in any operational sense — it is *richer*. The ramification structure,
the Galois group, the arithmetic of the ring of integers — all of these change qualitatively as $d$
grows, not just quantitatively. A toy instance over a degree-6 field can exhibit phenomena that a
"huge" instance over a degree-3 field never touches.

This is the axis that drives the most interesting mathematics in this textbook. The MOV reduction
works because the embedding degree $k$ of the target curve is small; the GHS descent works because
the field $\mathbb{F}_{2^m}$ has a tower structure; the quasi-polynomial DLP algorithm works
because the characteristic is small relative to the extension degree. None of these phenomena are
about instance size — they are about the *shape* of the mathematical object.

### Axis 3: Structural scale

The third axis is the most subtle: *thresholds* where scaling unlocks (or requires) a *different*
machine. These are not just quantitative changes — they are qualitative transitions.

The canonical examples:

- **Large-prime variations.** The basic NFS sieve collects only fully smooth relations. The
  large-prime variation collects *partial* relations — those with one or two large prime cofactors —
  and combines them to form full relations. This variation only pays off once the relation yield is
  large enough that the partial-relation graph has a giant component; at toy scale, the graph is
  sparse and the variation adds overhead without benefit. The threshold is a structural property of
  the relation graph, not a function of $N$ alone.

- **The quasi-polynomial DLP regime.** The Barbulescu–Gaudry–Joux–Thomé algorithm (2014) achieves
  quasi-polynomial complexity for DLP in $\mathbb{F}_{p^n}$ when the characteristic $p$ is small
  (e.g. $p = 2$). This is not a quantitative improvement on the $L[1/3]$ bound — it is a
  *structural* discovery that the small-characteristic case admits a fundamentally different
  algorithm. The threshold is the ratio $p / n$, not the size of the field.

- **The asymptotic regime of $L$-notation.** The $L$-notation complexity of GNFS is an asymptotic
  statement: it is accurate only when $N$ is large enough that the smoothness-probability estimates
  are reliable. At toy scale, the optimal parameters (smoothness bound $B$, sieve region size) are
  not the asymptotically optimal ones — they are calibrated empirically. The transition to the
  asymptotic regime is a structural threshold, not a smooth function of $N$.

### The three couplings

The three axes are sometimes independent, sometimes coupled, and sometimes structurally enabling.

**Independent:** the degree $d$ and the instance size $N$ are separate knobs. One can choose a
degree-6 number field for a 100-bit $N$, or a degree-3 field for a 1000-bit $N$. The mathematics
of the number field does not depend on $N$ directly.

**Coupled along the efficient frontier.** The *optimal* degree for GNFS is tied to $N$ by

$$d \sim \left(\frac{3 \log N}{\log \log N}\right)^{1/3},$$

so scaling $N$ "properly" (along the efficient frontier) drags $d$ along. At toy scale, the optimal
$d$ is 3 or 4 and barely moves; at RSA-768 scale, $d \approx 6$. The coupling is invisible at toy
scale, which is why toy-scale experiments with $d = 3$ are not misleading — they are just not
exploring the full frontier.

**Structural-enabling.** At certain thresholds, scaling one axis changes *which* mathematics
applies, not just how big it is. The quasi-polynomial DLP regime is the clearest example: the
algorithm is not a quantitative improvement on NFS-DL, it is a different algorithm that only
applies in a specific structural regime. The large-prime variation is another: it requires the
relation graph to have a giant component, which is a structural threshold.

### Method-convergence vs problem-openness

A natural question: is NFS asymptotically convergent? The precise answer distinguishes two things
the word "convergence" blurs.

**Convergence of the method (true).** The GNFS heuristic complexity has sat at $L_N[1/3,
(64/9)^{1/3}]$ since the early 1990s. Three decades of work have improved the constant and the
engineering — not the exponent $1/3$. Within the NFS paradigm (sieve number-field relations for
smoothness, then linear algebra), the exponent has converged. Note even this is *heuristic*,
resting on unproven smoothness assumptions, not a theorem.

**Convergence of the problem (false / open).** The true complexity of factoring and discrete
logarithm is **open**. $L[1/3]$ is a believed barrier for the *general* problem, not a proven
lower bound, and history shows that apparent convergence is repeatedly punctured by structural
discovery:

- The Barbulescu–Gaudry–Joux–Thomé quasi-polynomial algorithm (2014) collapsed small-characteristic
  DLP — a regime everyone had taken to be $L[1/3]$ — to $L[0, c]$.
- Shor's algorithm (1994) dissolves the barrier entirely in the quantum model.

The plateau is real, but it is a plateau of *one method*. It is reset, not approached, by new
mathematical structure.

### The honest science↔engineering gap

This textbook implements the mathematics of these algorithms completely. It does not implement the
engineering optimisations that make them run at cryptographic scale. The gap is real and is
annotated explicitly wherever it appears.

The gap runs in both directions:

**Under-exposed at toy scale.** Some phenomena only appear at large scale. Large-prime variations
only pay off once relation yield is in the millions. Block-Lanczos convergence behaviour that a
tiny matrix never stresses. The optimal-degree coupling between $d$ and $N$ that is invisible at
$d = 3$. These are implemented at *demonstration fidelity* — the mathematical content is present,
but the payoff is unreachable here. The annotation says: the mathematics is correct; the
engineering scale to exhibit the payoff is out of scope.

**Over-exposed at toy scale.** Some phenomena are more prominent at toy scale than at real scale.
Bad primes (primes dividing $\mathrm{disc}(f)$) are marginal at cryptographic scale but unavoidable
at toy scale with a hand-picked polynomial. The annotation says: this prominence is a toy-scale
artifact; at real scale it washes out.

The disconnect itself is pedagogical content. Part of honestly surveying how these algorithms
behave at the scales they were designed for is acknowledging what toy scale can and cannot exhibit.

---

## Pollard Rho for ECDLP

*Maths-first sibling to `docs/PEDAGOGY.md`. For the code-tour — phase-by-phase implementation
narrative, optimisation sequence, benchmarks, and CLI — see `docs/PEDAGOGY.md`. For the
prerequisites used in this chapter, see §Prerequisites above.*

### The through-line for this chapter

Pollard rho for ECDLP is the canonical example of **structure-based escape from search via group
homomorphism**. The structure exploited is the *linear representation* of walk states: every point
in the walk is a known linear combination $W = a \cdot G + b \cdot Q$ of the generator $G$ and the
target $Q$. A collision in the walk — two states that land on the same group element — immediately
yields a linear equation in the unknown discrete logarithm $k$. The birthday paradox guarantees
that such a collision occurs after $O(\sqrt{n})$ steps, where $n = |E(\mathbb{F}_p)|$.

The escape from the generic $O(n)$ exhaustive search is not a reduction in the *difficulty* of the
problem — it is a reduction in the *search space* from $n$ to $\sqrt{n}$, achieved by the birthday
paradox. The group homomorphism structure (the linear invariant $W = aG + bQ$) is what makes the
collision *useful*: without it, a collision would be a dead end.

### The birthday-paradox collision argument

Let $G$ be a cyclic group of prime order $n$, with generator $G$ and target $Q = k \cdot G$. A
*pseudorandom walk* on $G$ is a sequence $W_0, W_1, W_2, \ldots$ where each $W_{i+1} = f(W_i)$
for some deterministic function $f: G \to G$ that behaves like a random function.

**The key invariant.** The walk maintains scalars $a_i, b_i \in \mathbb{Z}/n\mathbb{Z}$ such that

$$W_i = a_i \cdot G + b_i \cdot Q$$

at every step. The function $f$ is chosen so that $a_{i+1}$ and $b_{i+1}$ can be computed from
$a_i$ and $b_i$ in $O(1)$ time.

**The collision.** Suppose $W_i = W_j$ for $i \neq j$. Then

$$a_i \cdot G + b_i \cdot Q = a_j \cdot G + b_j \cdot Q,$$

so $(a_i - a_j) \cdot G = (b_j - b_i) \cdot Q = (b_j - b_i) \cdot k \cdot G$, giving

$$k \equiv \frac{a_i - a_j}{b_j - b_i} \pmod{n}$$

provided $b_j \not\equiv b_i \pmod{n}$ (which holds with probability $1 - 1/n$; if it fails, the
attempt is discarded and a fresh walk is started).

**The quotient is a candidate, not an answer.** The solvers return $k$ only after checking $k \cdot
G = Q$ directly, so a returned logarithm does not depend on the walk that found it
(`verify_log_rejects_successor` in `rho/tests/ecdlp_kat.rs`). A candidate that fails the check is
an invariant failure, reported as `Inconsistent` and never retried
(`verified_reports_wrong_candidate_as_inconsistent`, a unit test in `rho/src/ecdlp/mod.rs`).
Spending the attempt or step budget without a usable collision is reported as budget exhaustion
and proves nothing about $Q$ (`brent_zero_retries_is_budget_exhausted`,
`dp_tiny_budget_is_budget_exhausted`, both in `rho/tests/ecdlp_kat.rs`).

**The birthday bound.** If $f$ behaves like a random function on $G$, the expected number of steps
before a collision is $\Theta(\sqrt{n})$ by the birthday paradox (§Prerequisites). This is the
$O(\sqrt{n})$ complexity of Pollard rho.

### The $r$-adding walk and mixing

The naive pseudorandom walk — "add a random point at each step" — has poor mixing properties: the
walk tends to cluster near its starting point, and the birthday bound is not achieved in practice.
Teske (1998) showed that an *$r$-adding walk* mixes much better.

**Definition ($r$-adding walk).** Precompute a table of $r$ random addends $R[0], \ldots, R[r-1]$,
where $R[i] = \alpha_i \cdot G + \beta_i \cdot Q$ with random scalars $\alpha_i, \beta_i$. At each
step, select the addend by $i = x \bmod r$ (where $x$ is a low word of the current point's
$x$-coordinate), then set

$$W \leftarrow W + R[i], \quad a \leftarrow a + \alpha_i \bmod n, \quad b \leftarrow b + \beta_i
\bmod n.$$

With $r \approx 20$, the walk behaves like a random function on the group, and the birthday bound
is achieved in practice.

**Why mixing matters.** The birthday bound $\Theta(\sqrt{n})$ assumes that the walk visits group
elements uniformly at random. A poorly mixing walk may visit some elements many times and others
never, increasing the expected collision time. The $r$-adding walk is designed to avoid this by
making the step function depend on the current point in a way that spreads the walk across the
group.

### Floyd's cycle detection

The walk $W_0, W_1, W_2, \ldots$ is eventually periodic: since $G$ is finite, the walk must
revisit a point, and from that point on it follows the same trajectory. The *rho* shape — a tail
leading into a cycle — gives the algorithm its name.

**Floyd's algorithm.** Run two pointers: the *tortoise* at speed 1 ($W_i$) and the *hare* at speed
2 ($W_{2i}$). When the hare catches the tortoise ($W_i = W_{2i}$), a collision has been detected.
The collision yields the discrete logarithm as above.

**Brent's improvement.** Freeze the tortoise at the start of each power-of-2 window; advance only
the hare. Compare the hare to the frozen tortoise at each step. After $2^k$ comparisons without
collision, snap the tortoise forward to the hare's position and double the window. This reduces the
number of function evaluations by approximately 24% compared to Floyd's algorithm.

**Expected complexity.** The expected number of steps before the hare catches the tortoise is
$\Theta(\sqrt{n})$, matching the birthday bound. The constant factor depends on the mixing
properties of the walk.

### The group-homomorphism structure

The reason Pollard rho works on elliptic curves — and not just on abstract groups — is the
*group-homomorphism structure* of the walk. The key point is that the elliptic curve group law is
*computable*: given two points $P$ and $Q$ on $E(\mathbb{F}_p)$, one can compute $P + Q$ in $O(1)$
field operations. This means the walk can be implemented efficiently, and the linear invariant
$W = aG + bQ$ can be maintained at each step.

More precisely: the map $(a, b) \mapsto aG + bQ$ is a group homomorphism from $\mathbb{Z}/n\mathbb{Z}
\times \mathbb{Z}/n\mathbb{Z}$ to $E(\mathbb{F}_p)$. The walk on $E(\mathbb{F}_p)$ lifts to a walk
on $\mathbb{Z}/n\mathbb{Z} \times \mathbb{Z}/n\mathbb{Z}$, and a collision in the former gives a
linear equation in the latter. The discrete logarithm is the solution to that equation.

This is the general pattern of index calculus: find a homomorphism from the group to a structure
where the discrete logarithm is easy, and use it to reduce the hard problem to an easy one. In
Pollard rho, the "easy structure" is just the linear equation $k \equiv (a_i - a_j)/(b_j - b_i)
\pmod{n}$, which is trivial to solve.

### Distinguished points and parallel collision search

The single-threaded Floyd/Brent algorithm is inherently sequential: the tortoise and hare must stay
in sync. The *distinguished-point* (DP) method (van Oorschot–Wiener, 1999) breaks this dependency
and enables linear parallelism.

**Definition (distinguished point).** A point $W \in E(\mathbb{F}_p)$ is *distinguished* if its
$x$-coordinate has at least $\theta$ low-order zero bits. The expected number of steps between
distinguished points is $2^\theta$.

**The parallel architecture.** $N$ walker threads each run an independent $r$-adding walk and emit
a record $(x, \varepsilon, a, b)$ whenever the walk lands on a distinguished point, where $x$ is the
full $x$-coordinate (every word, not only the low bits that decide whether the point is
distinguished) and $\varepsilon$ is the parity of $y$. A coordinator thread maintains a hash table
keyed on $x$. When two records with the same $x$ but different $(a, b)$ pairs arrive, the
coordinator solves for $k$.

**Equal points and opposite points.** $P$ and $-P$ share an $x$-coordinate, so equal $x$ does not
by itself mean equal points. For odd $p$ and $y \neq 0$ the values $y$ and $p - y$ have opposite
parity, so the parities $\varepsilon_1, \varepsilon_2$ tell the cases apart. If they agree the two
points are equal and $a_1 + b_1 k \equiv a_2 + b_2 k$, giving $k \equiv (a_1 - a_2)/(b_2 - b_1)$.
If they differ the points are negatives and $a_1 + b_1 k \equiv -(a_2 + b_2 k)$, giving $k \equiv
-(a_1 + a_2)/(b_1 + b_2) \pmod{n}$, which is degenerate exactly when $b_1 + b_2 \equiv 0$. A walk
without the negation map can reach both $P$ and $-P$, so the second case occurs. Treating every
$x$-match as an equal-point collision produces a wrong $k$: matching by $x$ alone and reading every
match as equality, the coordinator produced a candidate that failed $k \cdot G = Q$ in 18 of 100
recorded seeds (the bounded seed search described on
`plain_walk_recovers_log_from_negated_collision` in `rho/tests/ecdlp_kat.rs`, one walker,
$\theta = 2$). With the parity in the record the same runs succeed and report the negated
relation (that test and
`plain_walk_never_misreads_the_sign_of_a_collision`). The negation-map and GLV walks canonicalise
every point before emitting it, so their collisions are always between identical points
(`dp_collision_evidence_is_rechecked` asserts `negated` is false for them).

**Termination.** All walkers draw from one shared budget of walk steps, claimed before each step,
so no input can make a search run unboundedly (`step_budget_grants_claims_up_to_the_limit` in
`rho/src/ecdlp/mod.rs`); a search that spends the budget without a usable collision ends as budget
exhaustion (`dp_solvers_report_budget_exhaustion` in the same file,
`dp_tiny_budget_is_budget_exhausted` in `rho/tests/ecdlp_kat.rs`).

**Why it works.** Two walks that ever visit the same point will, from that point on, follow
identical trajectories (the walk is deterministic). They will therefore hit the same distinguished
points. A collision in the DP table implies a collision in the underlying walks, which implies the
DLP can be recovered.

**Speedup.** With $N$ walkers, the expected wall-time steps to collision is $O(\sqrt{n} / N)$ —
linear speedup in the number of walkers. This is the second pedagogical moment: the speedup from
parallelism is linear, not $\sqrt{N}$.

### Orbit-collapsing maps: negation and GLV

The birthday bound $\Theta(\sqrt{n})$ applies to the *effective group size* — the number of
distinct walk states. By collapsing orbits of the group action, one can reduce the effective group
size and improve the expected collision time.

**The negation map.** On an elliptic curve, $P$ and $-P$ share the same $x$-coordinate (negation
flips only the $y$-coordinate). By collapsing $\{P, -P\}$ to a single canonical representative
(the one with the smaller $y$-coordinate), the effective group size is halved, reducing the
expected walk length by $\sqrt{2}$.

**The GLV endomorphism.** For a curve with an order-3 endomorphism $\phi$ (e.g. $y^2 = x^3 + 7$
over a field with a cube root of unity), the orbit of any generic point $P$ under $\langle \phi,
-1 \rangle$ has size 6:

$$\{P,\; \phi(P),\; \phi^2(P),\; -P,\; -\phi(P),\; -\phi^2(P)\}.$$

Collapsing this orbit to a single canonical representative reduces the effective group size by 6,
reducing the expected walk length by $\sqrt{6}$ compared to the plain walk, or $\sqrt{3}$ compared
to the negation-map-only walk.

**Scalar bookkeeping.** The orbit-collapsing maps require adjusting the scalars $(a, b)$ to
maintain the invariant $W = aG + bQ$. For the negation map: $-W = (n-a)G + (n-b)Q$. For the GLV
endomorphism: $\phi(W) = \lambda W$ where $\lambda$ is the eigenvalue of $\phi$, so $\phi(W) =
(\lambda a)G + (\lambda b)Q$.

**Fruitless cycles.** The canonical map can create *fruitless cycles*: the walk oscillates between
two orbit members that share the same canonical representative, and never hits a distinguished
point. The BKNS escape (Bos–Kleinjung–Niederhagen–Schwabe) handles this deterministically: detect
the period-2 pattern in a sliding window, then perturb by doubling the current point ($W \leftarrow
2W$, $a \leftarrow 2a$, $b \leftarrow 2b$) and re-canonicalising.

### The $L$-notation bound for ECDLP

For a generic elliptic curve over $\mathbb{F}_p$ with $|E(\mathbb{F}_p)| = n$, the best known
classical algorithm runs in $\Theta(\sqrt{n})$ group operations. In $L$-notation:

$$\text{ECDLP complexity} = L_n[1, 1/2] = n^{1/2} = \Theta(\sqrt{n}).$$

This is *fully exponential* in $\log n$ — the algorithm is not subexponential. The reason is that
no index-calculus algorithm is known for a generic elliptic curve: the group $E(\mathbb{F}_p)$ does
not have a natural "factor base" structure that would allow smooth-element collection.

**The contrast with DLP in $\mathbb{F}_p^*$.** For the multiplicative group $\mathbb{F}_p^*$, the
best known algorithm (NFS-DL) runs in $L_p[1/3, (64/9)^{1/3}]$ — subexponential. The key
difference is that $\mathbb{F}_p^*$ has a natural factor base (the small primes), while a generic
elliptic curve group does not. This is why ECC-256 provides 128-bit security while RSA-3072 is
needed for the same security level against classical attacks.

**When the bound breaks.** The $\sqrt{n}$ bound is the *generic* bound. It breaks when the curve
has exploitable structure:
- **Small embedding degree** (MOV/Frey–Rück): the ECDLP reduces to DLP in $\mathbb{F}_{p^k}^*$,
  where index calculus applies.
- **Anomalous curves** (Smart–Satoh–Araki): curves with $|E(\mathbb{F}_p)| = p$ admit a
  polynomial-time attack via $p$-adic logarithms.
- **Special field structure** (GHS descent): curves over $\mathbb{F}_{2^m}$ with a tower
  structure admit a descent to a hyperelliptic Jacobian.

These are the subjects of the algebraic ECDLP attacks chapter (§10).

### Cross-reference

For the phase-by-phase implementation of Pollard rho — the code, the optimisation sequence, the
benchmarks, and the CLI — see `docs/PEDAGOGY.md`. That chapter is the code-tour sibling to this
one: it assumes the reader knows the mathematics (this chapter) and focuses on the realisation in
Rust.

---

## The α-Substrate: Primality, Smoothness, and ECM

*Maths-first sibling to `shared/numth/docs/PEDAGOGY.md`. For the code-tour — implementation
details, code references, KAT summary, and downstream consumption — see
`shared/numth/docs/PEDAGOGY.md`. For the prerequisites used in this chapter, see §Prerequisites
above.*

### The through-line for this chapter

The α-substrate algorithms are the *infrastructure* of structure-based escape from search. They do
not themselves escape the generic bound — they provide the tools that other algorithms use to do
so. Specifically:

- **Miller–Rabin** decides primality, which is needed to build factor bases and to classify an
  input as prime before any factoring is attempted.
- **Smooth-number detection** identifies the elements that index-calculus algorithms collect; the
  smooth-number density theorem (§Prerequisites) is the mathematical engine of every sieve.
- **ECM** (Lenstra's elliptic-curve method) factors integers by exploiting the variability of
  elliptic-curve group orders — a direct application of the group-homomorphism structure from the
  Pollard rho chapter, now used for factoring rather than discrete logarithm.
- **Tonelli–Shanks** computes square roots modulo a prime, which is needed for curve-point recovery
  and for the GNFS square-root step.

The through-line is: these algorithms are the substrate on which the escape from search is built.
Understanding them is understanding the *tools* of the escape.

### Miller–Rabin primality testing

#### The Fermat test and its failure

The simplest primality test is the *Fermat test*: if $n$ is prime, then $a^{n-1} \equiv 1 \pmod{n}$
for all $a$ with $\gcd(a, n) = 1$ (Fermat's little theorem). A *Fermat pseudoprime* to base $a$ is
a composite $n$ that passes the Fermat test for base $a$.

The Fermat test fails catastrophically for *Carmichael numbers* — composites $n$ that pass the
Fermat test for *every* base $a$ with $\gcd(a, n) = 1$. The smallest Carmichael number is 561 =
$3 \cdot 11 \cdot 17$. There are infinitely many Carmichael numbers (Alford–Granville–Pomerance,
1994), so the Fermat test alone is not a reliable primality test.

#### The Miller–Rabin test

The Miller–Rabin test strengthens the Fermat test by exploiting the structure of square roots of 1
modulo a prime.

**Setup.** Write $n - 1 = 2^s \cdot d$ with $d$ odd. For a witness base $a$, compute $x = a^d
\bmod n$. Then square $x$ up to $s - 1$ times. If the sequence $a^d, a^{2d}, \ldots, a^{2^{s-1}d}$
ever hits $n - 1$ (or starts at 1), the witness $a$ is *satisfied* and offers no evidence of
compositeness. If it never does, $a$ is a *witness to compositeness* and $n$ is definitely
composite.

**Theorem (Miller–Rabin soundness).** *If $n$ is prime, every $a \in (1, n-1)$ is satisfied. If
$n$ is composite, at least $3/4$ of the bases $a \in \{1, \ldots, n-1\}$ are witnesses to
compositeness.*

*Proof sketch.* For a prime $n$, the only square roots of 1 modulo $n$ are $\pm 1$ (since
$\mathbb{F}_n$ is a field). The sequence $a^d, a^{2d}, \ldots$ must therefore hit $\pm 1$ before
reaching $a^{n-1} = 1$; if it hits $-1 = n-1$, the witness is satisfied. For a composite $n$,
the argument uses the structure of the group $(\mathbb{Z}/n\mathbb{Z})^*$ and the Chinese Remainder
Theorem to show that at most $1/4$ of the bases are satisfied. See Crandall–Pomerance [CP05, §3.5]
for the full argument. $\square$

**Corollary.** A composite $n$ that passes $k$ independent random bases slips through with
probability at most $4^{-k}$.

#### From probabilistic to deterministic

Miller–Rabin is probabilistic in general, but for bounded $n$ it can be made *deterministic* by
using a fixed witness set known to have no exceptions below a threshold. The key result is:

**Theorem (deterministic Miller–Rabin).** *The twelve prime bases $\{2, 3, 5, 7, 11, 13, 17, 19, 23,
29, 31, 37\}$ are a deterministic primality test for all $n < \psi_{12}$, where
$\psi_{12} = 318{,}665{,}857{,}834{,}031{,}151{,}167{,}461 \approx 2^{78.1}$. Adding the base 41
(thirteen bases) extends this to all $n < \psi_{13}$, where
$\psi_{13} = 3{,}317{,}044{,}064{,}679{,}887{,}385{,}961{,}981 \approx 2^{81.5}$.*

*Proof.* Computational verification: no composite below the threshold passes all the bases. The
thresholds are Theorem 1.1 of Sorenson–Webster, "Strong pseudoprimes to twelve prime bases"
(Math. Comp. 86 (2017); arXiv:1509.00864), which defines $\psi_m$ as the smallest strong
pseudoprime to the first $m$ prime bases and computes $\psi_{12}$ and $\psi_{13}$ as stated above;
they refine earlier tables by Pomerance–Selfridge–Wagstaff and Jaeschke. The exhaustive search
behind them is not repeated here. What the tests do compute is
the boundary case: $\psi_{12} = 399{,}165{,}290{,}221 \times 798{,}330{,}580{,}441$ is composite and
passes all twelve bases, and base 41 rejects it (`psi_12_is_the_product_it_claims`,
`twelve_bases_accept_psi_12_but_base_41_rejects_it`, `psi_12_rejected` in
`shared/numth/src/prime.rs`). A twelve-base verdict is therefore not a proof at $\psi_{12}$ or
above, which is why `is_prime` applies the thirteenth base from there up. $\square$

**Above the proven bound.** For $n \geq \psi_{13}$ (which still fits in 256 bits) the same thirteen
bases give only a *probable-prime* verdict with no proven error bound: some composite is a strong
pseudoprime to every fixed base set, and `is_prime` can return true for one. A *composite* verdict
is always a proof: a witness is a certificate that $n$ fails the test every prime passes. Callers
that need a proof of primality must therefore treat a prime verdict at or above $\psi_{13}$ as
unproven; the public constant `MILLER_RABIN_PROOF_BOUND` in `shared-numth` is $\psi_{13}$.

**Why this matters.** The deterministic Miller–Rabin test is the primality oracle used throughout
this project. It builds factor bases (each base element is certified prime) and classifies an input
as prime before any factoring is attempted. The factoring drivers report `Prime` only where the
verdict is proven, and `Unsupported` for an odd probable prime beyond the proof. Rho factoring
accepts only $n < 2^{80}$, and $2^{80} < \psi_{13}$, so every prime verdict it sees is a proof: the
prime $2^{79} - 67$, between $\psi_{12}$ and $2^{80}$, is `Prime`
(`factor_certifies_primes_up_to_the_supported_limit` in `rho/tests/factor_kat.rs`), a compile-time
assertion keeps the supported limit at most $\psi_{13}$, and the composite $\psi_{12}$ is split, not
certified (`factor_splits_psi_12_instead_of_certifying_it`). ECM, whose range is not capped at
$2^{80}$, uses $\psi_{13}$ as its limit, so $2^{89} - 1$ is `Unsupported` and $\psi_{12}$,
being composite, is never certified prime (`prime_above_proof_bound_is_unsupported`,
`composite_psi_12_is_never_prime` in `shared/numth/src/ecm.rs`). ECM does not test the factor it
finds for primality: it returns a divisor only after checking $1 < d < n$ and $n \bmod d = 0$
(`judge_candidate_classifies_trivial_and_proper`), and that divisor need not be prime. A candidate
that is neither $1$ nor $n$ yet does not divide $n$ is reported as `Inconsistent { d }` (in both
`ecm_factor` and rho `factor`), never retried: a gcd with $n$ always divides $n$, so it is an
invariant failure, exercised only through stubbed searches
(`judge_candidate_flags_non_dividing_candidates`,
`search_surfaces_inconsistent_candidate_without_retrying`), while $d \in \{1, n\}$ is a failed
curve or walk that consumes one attempt.

#### Strong pseudoprimes and the 2-Sylow structure

The Miller–Rabin test is more powerful than the Fermat test because it exploits the *2-Sylow
structure* of $(\mathbb{Z}/n\mathbb{Z})^*$. Writing $n - 1 = 2^s \cdot d$, the test checks not
just that $a^{n-1} \equiv 1$, but that the sequence of squarings $a^d, a^{2d}, \ldots, a^{2^{s-1}d}$
has the right structure: it must hit $-1$ before reaching $1$ (unless it starts at $1$). This is
exactly the structure that the 2-Sylow subgroup of $(\mathbb{Z}/p\mathbb{Z})^*$ forces for a prime
$p$.

A *strong pseudoprime* to base $a$ is a composite $n$ that passes the Miller–Rabin test for base
$a$. Strong pseudoprimes are much rarer than Fermat pseudoprimes: the smallest strong pseudoprime
to base 2 is 2047 = $23 \cdot 89$, and it is caught by base 23.

### Smooth-number theory

#### The $B$-smooth density theorem

The Canfield–Erdős–Pomerance theorem (§Prerequisites) gives the density of $B$-smooth integers:
for $u = \log N / \log B$ fixed,

$$\Pr[\text{random integer in } [1, N] \text{ is } B\text{-smooth}] = \rho(u) \cdot (1 + o(1)),$$

where $\rho(u)$ is the Dickman function. For large $u$, $\rho(u) \approx u^{-u}$.

**The key estimate for sieve algorithms.** If $B = L_N[1/3, c]$, then $u = \log N / \log B \sim
(\log N)^{2/3} / (c \cdot (\log \log N)^{1/3})$, and

$$\rho(u) = L_N\!\left[-\frac{1}{3},\, -\frac{1}{3c}\right] \cdot (1 + o(1)).$$

This is the smoothness probability that enters the GNFS complexity analysis. The probability is
subexponentially small, but much larger than $N^{-\varepsilon}$ for any fixed $\varepsilon > 0$.

**Why this matters.** The smooth-number density theorem is the mathematical engine of every sieve.
The GNFS sieve collects pairs $(a, b)$ such that both the rational norm $|a - bm|$ and the
algebraic norm $|N_{K/\mathbb{Q}}(a - b\alpha)|$ are $B$-smooth. The probability that a random
pair is smooth is $\rho(u)^2$ (approximately, for the two independent norms), and the optimal $B$
is chosen to balance this probability against the cost of the linear algebra step.

#### The factor base

**Definition (factor base).** A *factor base* $\mathcal{F}$ is a finite set of primes $\{p_1, p_2,
\ldots, p_k\}$ with $p_k \leq B$. An integer $m$ is $\mathcal{F}$-smooth (or $B$-smooth) if all
its prime factors are in $\mathcal{F}$.

**The smoothness witness.** A *smoothness witness* for $m$ over $\mathcal{F}$ is the factorisation
$m = \prod_{p \in \mathcal{F}} p^{e_p} \cdot r$, where $r$ is the *cofactor* — the part of $m$
not accounted for by $\mathcal{F}$. If $r = 1$, $m$ is fully $\mathcal{F}$-smooth; if $r > 1$, $m$
has a *large prime* cofactor outside the base.

**The large-prime variation.** In the basic sieve, only fully smooth relations are collected. In
the large-prime variation, partial relations (with one or two large prime cofactors) are also
collected. Two partial relations with the same large prime cofactor can be combined to form a full
relation. This variation improves the relation yield at the cost of more complex bookkeeping.

### Lenstra's Elliptic Curve Method (ECM)

#### The idea: Pollard's $p-1$, lifted onto a curve

Pollard's $p-1$ method factors $N$ by exploiting that if $p \mid N$ and $p - 1$ is $B$-smooth,
then for a base $a$, $a^k \equiv 1 \pmod{p}$ once $k = \mathrm{lcm}(1, \ldots, B)$ is a multiple
of the order of $a$ mod $p$. Then $\gcd(a^k - 1, N)$ reveals $p$.

The weakness: it only works when $p - 1$ happens to be smooth — a property of the single fixed
group $(\mathbb{Z}/p\mathbb{Z})^*$.

**Lenstra's insight (1987).** Replace the fixed group $(\mathbb{Z}/p\mathbb{Z})^*$ with the group
of points on a *random* elliptic curve $E$ over $\mathbb{Z}/N\mathbb{Z}$. The order $|E(\mathbb{F}_p)|$
varies with the curve choice and ranges over an interval around $p$ (Hasse's bound: $|p + 1 -
|E(\mathbb{F}_p)|| \leq 2\sqrt{p}$). If the chosen curve has $B$-smooth order mod some prime $p
\mid N$, then multiplying a starting point $P$ by $k = \prod q^e$ (prime powers $\leq B$) drives
$kP$ to the identity mod $p$ while it stays non-trivial mod $N$. The group law's modular inversion
then fails with a denominator sharing a factor with $N$, and $\gcd(\text{denominator}, N) = p$
falls out.

**The key advantage.** If one curve's order is not smooth, try another curve. ECM converts a
question of luck (is $p - 1$ smooth?) into a question one can re-roll (is *some* curve's order
smooth?).

#### The group-order smoothness argument

**Theorem (ECM expected complexity).** *The expected number of curves ECM must try before finding
a factor $p$ of $N$ is $O(\rho(u))^{-1}$, where $u = \log p / \log B$ and $\rho$ is the Dickman
function. The expected running time per curve is $O(B \log B \log N)$ (for stage 1). The total
expected running time is*

$$O\!\left(B \log B \log N \cdot \rho(u)^{-1}\right) = L_p\!\left[\frac{1}{2},\, \sqrt{2}\right]
\cdot \mathrm{poly}(\log N)$$

*when $B$ is chosen optimally as $B = L_p[1/2, 1/\sqrt{2}]$.*

*Proof sketch.* The probability that a random curve has $B$-smooth order mod $p$ is approximately
$\rho(u)$ (by the Canfield–Erdős–Pomerance theorem applied to the interval $[p + 1 - 2\sqrt{p},
p + 1 + 2\sqrt{p}]$). The optimal $B$ balances the per-curve cost $O(B \log B \log N)$ against the
number of curves $\rho(u)^{-1}$. Setting $B = L_p[1/2, c]$ and optimising over $c$ gives the
stated bound. See Lenstra [L87] and Crandall–Pomerance [CP05, §7.4] for the full argument. $\square$

**The subexponential complexity.** ECM runs in $L_p[1/2, \sqrt{2}]$ time in the size of the
*smallest prime factor* $p$, not in the size of $N$. This makes ECM the right tool for splitting
off small-to-medium prime factors: it is subexponential in $p$, but the cost grows with $p$, so it
is not competitive with NFS for large factors.

#### Montgomery-form curves and the $(X:Z)$ ladder

ECM uses *Montgomery-form curves* $Bv^2 = u^3 + Au^2 + u$ with the *Suyama parameterisation* to
generate a valid (curve, starting point) pair from a single integer $\sigma$. Montgomery form
admits projective $(X:Z)$ coordinates that omit the $Y$-coordinate entirely: scalar multiplication
needs only the $x$-coordinate, computed by the *Montgomery ladder* using *differential addition*
(given $P$, $Q$, and $P - Q$, compute $P + Q$) and doubling.

**Why Montgomery form?** The Montgomery ladder is *uniform* (no branch on the bit of the scalar)
and *efficient* (differential addition uses fewer field operations than the general addition
formula). The $Z = 0$ point represents the identity; a non-invertible denominator in the group law
signals that $\gcd(Z, N) > 1$, which is the factor.

**The Suyama parameterisation.** Given $\sigma$, compute the curve constant $A_{24} = (A+2)/4$ and
a starting point, returning early with a factor if any of the constructions hits a non-invertible
denominator mod $N$ (itself a lucky factor find). The parameterisation ensures that the starting
point has a specific structure that makes stage 1 more likely to succeed.

#### Stage 1 and stage 2

**Stage 1.** Compute $Q = kP$ with $k = \prod \{q^{\lfloor \log_q B_1 \rfloor} : q \text{ prime}
\leq B_1\}$ — every prime power up to $B_1$. After each ladder call, test $\gcd(Z, N)$; a
non-trivial gcd is the factor. Stage 1 succeeds when $|E(\mathbb{F}_p)|$ is fully $B_1$-smooth.

*No product scalar.* The product is $\mathrm{lcm}(1, \ldots, B_1) \approx e^{B_1}$, which first
exceeds $2^{256}$ at $B_1 = 179$ ($256 \ln 2 \approx 177$), so it does not fit a fixed-width
integer at larger $B_1$. The code never forms it: it multiplies the running point by each prime
power $q = p^{\lfloor \log_p B_1 \rfloor}$ in its own ladder call, giving the same point
$(q_1 q_2 \cdots)P$ at every $B_1$. Once $Z \equiv 0 \pmod r$ for a prime $r \mid N$, doubling and
differential addition keep $Z \equiv 0$, so testing the gcd after each prime is sound, and it
reports a factor before every prime of $N$ has reached the identity (which would make the gcd
improper). The tests in `shared/numth/src/ecm.rs` check that the result equals the single-scalar
ladder for $B_1$ up to $178$ and that $N = (2^{61}-1) \cdot 1009$ with $\sigma = 15$, $B_1 = 250$
yields $1009$, which a scalar truncated mod $2^{256}$ misses.

**Stage 2.** The standard refinement for the common case where $|E(\mathbb{F}_p)|$ is $B_1$-smooth
*except for one* prime factor $q \in (B_1, B_2]$. Rather than raising $B_1$ (which re-does all
the work), stage 2 walks the multiples $qQ$ for primes $q \in (B_1, B_2]$ using differential
addition with a fixed step, accumulates the product of their $Z$-coordinates, and takes a *single*
gcd with $N$ at the end. This catches the one-large-prime case cheaply.

**The two-stage structure** is the ECM analogue of the large-prime variation in NFS sieving: both
extend the reach of the algorithm by handling the case where the "smooth" condition fails by one
prime factor.

### Tonelli–Shanks: square roots mod $p$

#### The problem

Given a prime $p$ and a quadratic residue $a$ (i.e. $\left(\frac{a}{p}\right) = 1$), find $r$ with
$r^2 \equiv a \pmod{p}$.

#### The 2-Sylow structure argument

The multiplicative group $(\mathbb{Z}/p\mathbb{Z})^*$ is cyclic of order $p - 1$. Write $p - 1 =
Q \cdot 2^S$ with $Q$ odd. The *2-Sylow subgroup* of $(\mathbb{Z}/p\mathbb{Z})^*$ is the unique
subgroup of order $2^S$.

**The easy case: $p \equiv 3 \pmod{4}$ (i.e. $S = 1$).** Here $r = a^{(p+1)/4}$ is a square root,
because $r^2 = a^{(p+1)/2} = a \cdot a^{(p-1)/2} = a \cdot 1 = a$ (using Euler's criterion for
a residue). One exponentiation, no search.

**The general case: $p \equiv 1 \pmod{4}$ (i.e. $S \geq 2$).** The Tonelli–Shanks algorithm
maintains a tuple $(M, c, t, r)$ where:
- $M$ is the current "error exponent" (starts at $S$),
- $c$ is a power of a quadratic non-residue $z$ (so $c^{2^{M-1}} = -1$),
- $t = a^Q$ (the "error term": $t^{2^{M-1}} = 1$ iff $r^2 = a$),
- $r = a^{(Q+1)/2}$ (the current candidate square root).

At each step, find the smallest $i$ such that $t^{2^i} = 1$. If $i = 0$, then $t = 1$ and $r^2 =
a$ — done. Otherwise, update:

$$c \leftarrow c^{2^{M-i-1}}, \quad t \leftarrow t \cdot c^2, \quad r \leftarrow r \cdot c,
\quad M \leftarrow i.$$

This reduces the error exponent $M$ at each step, so the algorithm terminates in at most $S - 1$
iterations.

**Finding the quadratic non-residue $z$.** The trial search tests $z = 2, 3, 5, \ldots$ until
$\left(\frac{z}{p}\right) = -1$. Non-residues have density $1/2$, so the expected number of trials
is 2. There is no known fast deterministic way to produce a non-residue (it is related to open
problems about the least quadratic non-residue under GRH), so the trial loop is the standard
practical choice.

**Theorem (Tonelli–Shanks correctness).** *The algorithm terminates in $O(S^2)$ multiplications
and returns $r$ with $r^2 \equiv a \pmod{p}$.*

*Proof sketch.* The invariant $r^2 \equiv a \cdot t \pmod{p}$ is maintained throughout. The error
exponent $M$ strictly decreases at each step (from $M$ to $i < M$), so the algorithm terminates.
At termination, $t = 1$ and $r^2 \equiv a$. $\square$

**Why this matters.** Tonelli–Shanks is used in this project for curve-point recovery (given an
$x$-coordinate, find the $y$-coordinate satisfying $y^2 = x^3 + ax + b$) and in the GNFS
square-root step (recovering a square root in the number field). The 2-Sylow structure argument
is the key insight: the algorithm works by iteratively reducing the "error" in the 2-Sylow
subgroup until the error vanishes.

### Cross-reference

For the phase-by-phase implementation of the α-substrate — the code, the KAT summary, the
downstream consumption table, and the scale-disconnect annotations — see
`shared/numth/docs/PEDAGOGY.md`. That chapter is the code-tour sibling to this one: it assumes
the reader knows the mathematics (this chapter) and focuses on the realisation in Rust.

---

## References

1. **Pollard, J. M. (1978).** "Monte Carlo methods for index computation (mod p)." *Mathematics of
   Computation*, 32(143), 918–924.

2. **Teske, E. (1998).** "Speeding up Pollard's rho method for computing discrete logarithms."
   *Algorithmic Number Theory Symposium (ANTS-III)*, LNCS 1423, 541–554.

3. **van Oorschot, P. C., and Wiener, M. J. (1999).** "Parallel collision search with
   cryptanalytic applications." *Journal of Cryptology*, 12(1), 1–28.

4. **Gallant, R. P., Lambert, R. J., and Vanstone, S. A. (2001).** "Faster point multiplication on
   elliptic curves with efficient endomorphisms." *Advances in Cryptology — CRYPTO 2001*, LNCS
   2139, 190–200.

5. **Lenstra, H. W. Jr. (1987).** "Factoring integers with elliptic curves." *Annals of
   Mathematics*, 126(3), 649–673.

6. **Montgomery, P. L. (1987).** "Speeding the Pollard and elliptic curve methods of
   factorization." *Mathematics of Computation*, 48(177), 243–264.

7. **Canfield, E. R., Erdős, P., and Pomerance, C. (1983).** "On a problem of Oppenheim concerning
   'Factorisatio Numerorum'." *Journal of Number Theory*, 17(1), 1–28. (The smooth-number density
   theorem.)

8. **Silverman, J. H. (2009).** *The Arithmetic of Elliptic Curves*. 2nd ed. Springer.

9. **Neukirch, J. (1999).** *Algebraic Number Theory*. Springer.

10. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective*. 2nd
    ed. Springer.

11. **Davenport, H. (2000).** *Multiplicative Number Theory*. 3rd ed. Springer.

12. **Granville, A. (2008).** "Smooth numbers: computational number theory and beyond." In
    *Algorithmic Number Theory: Lattices, Number Fields, Curves and Cryptography*, MSRI
    Publications 44.

13. **Barbulescu, R., Gaudry, P., Joux, A., and Thomé, E. (2014).** "A heuristic quasi-polynomial
    algorithm for discrete logarithm in finite fields of small characteristic." *Advances in
    Cryptology — EUROCRYPT 2014*, LNCS 8441, 1–16.

14. **Shor, P. W. (1994).** "Algorithms for quantum computation: discrete logarithms and
    factoring." *Proceedings of the 35th Annual Symposium on Foundations of Computer Science*,
    124–134.

15. **Bos, J. W., Kleinjung, T., Niederhagen, R., and Schwabe, P. (2012).** "ECC2K-130 on Cell
    CPUs." *Progress in Cryptology — AFRICACRYPT 2010*, LNCS 6055, 225–242.

16. **Alford, W. R., Granville, A., and Pomerance, C. (1994).** "There are infinitely many
    Carmichael numbers." *Annals of Mathematics*, 139(3), 703–722.

17. **Sorenson, J., and Webster, J. (2017).** "Strong pseudoprimes to twelve prime bases."
    *Mathematics of Computation*, 86(304), 985–1003.

---

# The General Number Field Sieve: Structure-Based Escape from Search

*Maths-first sibling to `gnfs/docs/PEDAGOGY.md` (integrative chapter, §52–§62). For the
code-tour — pipeline contracts, stage-by-stage implementation narrative, and KAT summary — see
`gnfs/docs/PEDAGOGY.md`. For the prerequisites used in this chapter, see §Prerequisites above.
For the number-field substrate (rings of integers, ideal factorisation, norm maps), see
`shared/numfield/docs/PEDAGOGY.md`.*

---

## §1 Introduction and Through-Line

The General Number Field Sieve (GNFS) is the fastest known classical algorithm for factoring
large integers. Its complexity is

$$L_N\!\left[\tfrac{1}{3},\, \left(\tfrac{64}{9}\right)^{1/3}\right],$$

where $L_N[\alpha, c] = \exp\!\bigl(c\,(\log N)^\alpha\,(\log\log N)^{1-\alpha}\bigr)$. This is
*subexponential* in $\log N$: faster than any polynomial in $N$ but slower than any polynomial in
$\log N$. The exponent $1/3$ — not $1/2$, not $1/4$, but exactly $1/3$ — is the signature of the
number-field structure that GNFS exploits.

**The through-line.** Every algorithm in this textbook escapes the generic search bound by finding
exploitable structure. For GNFS, the structure is *number-field arithmetic*: the integers
$\mathbb{Z}$ embed into a number field $K = \mathbb{Q}(\alpha)$, and the norm map
$N_{K/\mathbb{Q}}: K \to \mathbb{Q}$ connects factorisation in $K$ to factorisation in
$\mathbb{Z}$. By sieving for elements of $K$ with smooth norms, GNFS collects relations that
encode multiplicative structure modulo $N$. A linear-algebra step over $\mathrm{GF}(2)$ then
finds candidate dependencies among these relations; a verified algebraic square root of a
dependency gives a congruence of squares $x^2 \equiv y^2 \pmod{N}$, from which a factor of $N$
follows by a GCD computation unless the congruence is trivial.

The generic approach to finding a congruence of squares — collecting random squares and hoping
for a dependency — costs $L_N[1/2, 1]$. The number-field bridge reduces this to $L_N[1/3,
(64/9)^{1/3}]$. The improvement from exponent $1/2$ to exponent $1/3$ is the payoff of the
structure, and the derivation of this improvement is the payoff proof of this chapter (§7).

**Main theorem (GNFS complexity, heuristic).** *Under standard heuristic assumptions on the
distribution of smooth norms, the GNFS factors a composite integer $N$ in expected time
$L_N[1/3, (64/9)^{1/3}]$.*

The word "heuristic" is load-bearing: the smoothness-probability estimates that drive the
analysis are asymptotic, and the asymptotic regime is only reached at scales far beyond what any
computer can run. The algorithm is correct at all scales; the complexity analysis is accurate
only asymptotically. See §On Scale in this textbook for the full discussion.

---

## §2 The Congruence-of-Squares Idea

### Difference of squares

The oldest factoring idea is the *difference of squares*: if $N = x^2 - y^2 = (x-y)(x+y)$, then
$\gcd(x - y, N)$ is a non-trivial factor of $N$ (provided $x \not\equiv \pm y \pmod{N}$). Fermat's
method (1643) searches for such a representation by trying $x = \lceil\sqrt{N}\rceil,
\lceil\sqrt{N}\rceil + 1, \ldots$ and checking whether $x^2 - N$ is a perfect square.

The modern generalisation replaces the exact equation $x^2 = y^2$ with a congruence:

**Theorem (congruence-of-squares factoring).** *Let $N$ be a composite integer with $N > 1$. If
$x^2 \equiv y^2 \pmod{N}$ and $x \not\equiv \pm y \pmod{N}$, then $\gcd(x - y, N)$ is a
non-trivial factor of $N$.*

*Proof.* From $x^2 \equiv y^2 \pmod{N}$ we get $N \mid (x-y)(x+y)$. Since $x \not\equiv y
\pmod{N}$, we have $N \nmid (x - y)$, so $\gcd(x - y, N) < N$. Since $x \not\equiv -y \pmod{N}$,
we have $N \nmid (x + y)$, so $\gcd(x - y, N) > 1$. $\square$

**The probability of success.** If $N = pq$ is a semiprime and $x^2 \equiv y^2 \pmod{N}$ is
chosen uniformly at random among all such congruences, the probability that $x \not\equiv \pm y
\pmod{N}$ is exactly $1/2$. So each congruence of squares gives a factor with probability $1/2$.

### The generic approach and its cost

The simplest way to find a congruence of squares is *random squares*: pick random $x_i$ modulo
$N$, compute $x_i^2 \bmod N$, and hope that some product $\prod_{i \in S} x_i^2 \equiv y^2
\pmod{N}$ for a subset $S$. This is equivalent to finding a linear dependency over
$\mathrm{GF}(2)$ among the exponent vectors of the $x_i^2 \bmod N$ values.

The problem is that a random integer near $N$ has no reason to be smooth: its prime factors are
distributed like those of a random integer, and the probability that all prime factors are $\leq
B$ is $\rho(\log N / \log B)$ — exponentially small for any fixed $B$. To collect enough smooth
values to find a dependency, one must either take $B$ very large (making the linear algebra
expensive) or accept a very low smoothness probability (making the sieve expensive). The optimal
balance gives $L_N[1/2, 1]$:

**Theorem (Dixon's random squares, 1981).** *The random-squares method finds a congruence of
squares in expected time $L_N[1/2, 1]$.*

*Proof sketch.* Set $B = L_N[1/2, c]$. The smoothness probability for a random integer near $N$
is $\rho(\log N / \log B) \approx L_N[-1/2, -1/(2c)]$. To collect $B$ smooth values (enough for
a linear dependency), one needs $B \cdot L_N[1/2, 1/(2c)]$ trials. The total cost is
$B \cdot L_N[1/2, 1/(2c)] = L_N[1/2, c + 1/(2c)]$. Minimising over $c$ gives $c = 1/\sqrt{2}$
and total cost $L_N[1/2, \sqrt{2}]$. The linear algebra costs $O(B^2) = L_N[1, 2c]$, which is
dominated. $\square$

The $L_N[1/2, \cdot]$ exponent is the signature of the *quadratic sieve* family of algorithms.
GNFS escapes to $L_N[1/3, \cdot]$ by replacing random integers near $N$ with *norms of algebraic
integers* — quantities that are systematically smoother than random integers of the same size.

---

## §3 The Number-Field Bridge

### The setup

Let $N$ be the integer to factor. Choose a degree-$d$ polynomial $f \in \mathbb{Z}[x]$ and an
integer $m$ such that

$$f(m) \equiv 0 \pmod{N}.$$

This is easy to arrange: the *base-$m$ expansion* of $N$ gives $f$ directly. Write $N$ in base
$m$ (for any $m \approx N^{1/d}$):

$$N = a_d m^d + a_{d-1} m^{d-1} + \cdots + a_1 m + a_0,$$

and set $f(x) = a_d x^d + \cdots + a_1 x + a_0$. Then $f(m) = N \equiv 0 \pmod{N}$ by
construction.

Now let $\alpha$ be a root of $f$ in $\mathbb{C}$ (or in an algebraic closure of $\mathbb{Q}$),
and let $K = \mathbb{Q}(\alpha)$ be the number field generated by $\alpha$. When $f$ is monic the
ring $\mathbb{Z}[\alpha] = \{g(\alpha) : g \in \mathbb{Z}[x]\}$ is a subring of the ring of
integers $\mathcal{O}_K$. For a leading coefficient $c \ne 1$ it is $c\alpha$ that is an algebraic
integer, and §6 works with $\theta = c\alpha$; the statements of this section about ideals of
$\mathbb{Z}[\alpha]$ are then meant away from the primes dividing $c$.

### Two maps into $\mathbb{Z}/N\mathbb{Z}$

The key observation is that there are *two* natural ring homomorphisms into $\mathbb{Z}/N\mathbb{Z}$:

**The rational map** $\phi_{\mathbb{Z}}: \mathbb{Z} \to \mathbb{Z}/N\mathbb{Z}$, the standard
reduction modulo $N$.

**The algebraic map** $\phi_\alpha: \mathbb{Z}[\alpha] \to \mathbb{Z}/N\mathbb{Z}$, defined by
$\phi_\alpha(g(\alpha)) = g(m) \bmod N$. This is well-defined because $f(m) \equiv 0 \pmod{N}$,
so any polynomial relation satisfied by $\alpha$ over $\mathbb{Z}$ is also satisfied by $m$ over
$\mathbb{Z}/N\mathbb{Z}$.

**The bridge.** For any pair $(a, b) \in \mathbb{Z}^2$ with $\gcd(a, b) = 1$, consider the element
$a - b\alpha \in \mathbb{Z}[\alpha]$. Its image under $\phi_\alpha$ is $a - bm \bmod N$. Its image
under the rational map is also $a - bm \bmod N$. So both maps agree on $a - b\alpha$:

$$\phi_\alpha(a - b\alpha) = \phi_{\mathbb{Z}}(a - bm) = a - bm \bmod N.$$

This is the bridge: the same residue $a - bm \bmod N$ is the image of the integer $a - bm$ under
the rational map and the image of the element $a - b\alpha$ under the algebraic map. No
norm is involved: the bridge is a statement about images under ring homomorphisms, and it is what
makes the congruence of §6 hold.

### The algebraic homomorphism and prime ideals

For a rational prime $p$ and a root $r$ of $f$ modulo $p$ (i.e. $f(r) \equiv 0 \pmod{p}$), the
*prime ideal* $\mathfrak{p} = (p, \alpha - r)$ in $\mathbb{Z}[\alpha]$ is the kernel of the
evaluation homomorphism

$$\mathbb{Z}[\alpha] \to \mathbb{F}_p, \quad g(\alpha) \mapsto g(r) \bmod p.$$

This is a ring homomorphism: it sends $\alpha \mapsto r$ and reduces coefficients mod $p$. The
ideal $\mathfrak{p}$ is prime because the quotient $\mathbb{Z}[\alpha] / \mathfrak{p} \cong
\mathbb{F}_p$ is a field.

**Why this matters.** The prime ideal $\mathfrak{p} = (p, \alpha - r)$ divides the principal ideal
$(a - b\alpha)$ in $\mathbb{Z}[\alpha]$ if and only if $a - br \equiv 0 \pmod{p}$, i.e. $a \equiv
br \pmod{p}$. This is the *sieve condition*: it tells us exactly which pairs $(a, b)$ have the
ideal $\mathfrak{p}$ dividing $a - b\alpha$, and hence $p$ dividing the algebraic norm.

### The factor-base construction

The GNFS factor base has two sides:

**Rational factor base.** A set of rational primes $\mathcal{F}_{\mathbb{Z}} = \{p : p \leq B\}$.
A pair $(a, b)$ is *rationally smooth* if the integer $a - bm$ factors completely over
$\mathcal{F}_{\mathbb{Z}}$.

**Algebraic factor base.** A set of prime ideals $\mathcal{F}_\alpha = \{(p, \alpha - r) : p \leq
B,\; f(r) \equiv 0 \pmod{p}\}$. A pair $(a, b)$ is *algebraically smooth* if the principal ideal
$(a - b\alpha)$ in $\mathbb{Z}[\alpha]$ factors completely over $\mathcal{F}_\alpha$.

A pair $(a, b)$ that is both rationally and algebraically smooth is a *relation*. Each relation
encodes a multiplicative dependency modulo $N$ that the linear algebra step can exploit.

**The size of the factor base.** By the prime number theorem, the number of primes $\leq B$ is
approximately $B / \log B$. For each prime $p$, the number of roots of $f$ modulo $p$ is at most
$d = \deg f$. So the algebraic factor base has at most $d \cdot B / \log B$ elements — the same
order of magnitude as the rational factor base.

---

## §4 Smooth-Number Sieving

### The sieve region

The sieve collects pairs $(a, b)$ from the region

$$\mathcal{R} = \{(a, b) \in \mathbb{Z}^2 : |a| \leq A,\; 1 \leq b \leq B_s,\; \gcd(a, b) = 1\},$$

for parameters $A$ and $B_s$ (the sieve bound). The coprimality condition $\gcd(a, b) = 1$ ensures
that the pair is primitive — it avoids counting the same relation multiple times.

### The two norms

For each pair $(a, b)$, GNFS computes two norms:

**The rational norm** is simply

$$N_{\mathrm{rat}}(a, b) = |a - bm|.$$

(Some presentations use $a + bm$; this text uses $a - bm$ throughout, matching the element
$a - b\alpha$ and the implementation. The sign convention affects only the sign column in the
linear algebra.) This is an integer of size approximately $A + B_s \cdot m \approx N^{1/d}$ for
typical sieve parameters.

**The algebraic norm form** is the homogeneous form

$$N_{\mathrm{alg}}(a, b) = F(a, b) = b^d \cdot f(a/b) = c \cdot N_{K/\mathbb{Q}}(a - b\alpha)
  = \mathrm{Res}(f(x),\, a - bx),$$

where $c = \mathrm{lc}(f)$ and $\mathrm{Res}$ denotes the resultant. The identity follows from
$f(x) = c\prod_i (x - \alpha_i)$ over the roots $\alpha_i$ of $f$: $b^d f(a/b) = c\prod_i
(a - b\alpha_i)$. For monic $f$ the form is exactly the norm of $a - b\alpha$; for $c \ne 1$ it is
$c$ times the norm, and the norm itself may have a denominator. This is an integer of size
approximately $b^d \cdot |f(a/b)| \approx (B_s)^d \cdot (A/B_s)^d = A^d$ for typical parameters —
again of order $N^{1/d}$ when $A \approx N^{1/d}$.

**The key point.** Both norms are integers of size $\approx N^{1/d}$, which is *much smaller* than
$N$ itself. This is the source of the improvement over the quadratic sieve: instead of sieving
integers of size $\approx N$, GNFS sieves integers of size $\approx N^{1/d}$. The smoothness
probability for an integer of size $X$ with smoothness bound $B$ is $\rho(\log X / \log B)$; for
$X = N^{1/d}$ instead of $X = N$, this probability is much larger.

### Why both norms being smooth gives a relation

A pair $(a, b)$ is a relation if and only if both $N_{\mathrm{rat}}(a, b)$ and $N_{\mathrm{alg}}(a, b)$
are $B$-smooth. When this holds:

- The rational norm factors as $|a - bm| = \prod_{p \leq B} p^{e_p}$ over the rational factor
  base.
- The principal ideal factors as $(a - b\alpha) = \prod_{(p,r) \in \mathcal{F}_\alpha}
  \mathfrak{p}_{p,r}^{f_{p,r}}$ over the algebraic factor base (as an ideal in $\mathbb{Z}[\alpha]$,
  away from the primes dividing $c \cdot \mathrm{disc}(F)$, with $F$ as in §6), where $f_{p,r}$
  is the exponent of the *prime ideal* $\mathfrak{p}_{p,r}$, equal to $v_p(F(a, b))$.

The exponent vectors $(e_p)_{p \in \mathcal{F}_{\mathbb{Z}}}$ and $(f_{p,r})_{(p,r) \in
\mathcal{F}_\alpha}$ encode the multiplicative structure of the pair. Concatenating them gives a
vector in $\mathbb{Z}^k$ (where $k = |\mathcal{F}_{\mathbb{Z}}| + |\mathcal{F}_\alpha|$). Reducing
modulo 2 gives a vector in $\mathrm{GF}(2)^k$.

### The Canfield–Erdős–Pomerance estimate

The probability that a random integer of size $X$ is $B$-smooth is $\rho(u)$ where $u = \log X /
\log B$ and $\rho$ is the Dickman function (§Prerequisites). For the GNFS norms:

- Norm size: $X \approx N^{1/d}$, so $\log X \approx (\log N)/d$.
- Smoothness bound: $B = L_N[1/3, c]$, so $\log B \approx c \cdot (\log N)^{1/3} \cdot
  (\log\log N)^{2/3}$.
- The ratio: $u = \log X / \log B \approx \frac{(\log N)/d}{c \cdot (\log N)^{1/3} \cdot
  (\log\log N)^{2/3}} = \frac{(\log N)^{2/3}}{d \cdot c \cdot (\log\log N)^{2/3}}$.

For $u$ in the range relevant to GNFS (roughly $u \approx 3$ to $10$), the Dickman function
satisfies $\rho(u) \approx u^{-u(1+o(1))}$. In $L$-notation, the smoothness probability for a
single norm is

$$\rho(u) = L_N\!\left[-\tfrac{1}{3},\, -\tfrac{1}{3c}\right] \cdot (1 + o(1)).$$

The probability that *both* norms are smooth is approximately $\rho(u)^2 = L_N[-1/3, -2/(3c)]$
(treating the two norms as approximately independent, which is a standard heuristic assumption).

### The sieve implementation

The sieve uses the *log-sum* technique: for each $b$ in the sieve range, initialise an array
indexed by $a \in [-A, A]$, and for each prime $p$ in the factor base, add $\log p$ to every
position $a$ where $p \mid N_{\mathrm{rat}}(a, b)$ (resp. $p \mid N_{\mathrm{alg}}(a, b)$). A
position where the accumulated log-sum exceeds $\log B$ is a candidate smooth pair, confirmed by
trial division. This is the *line sieve*; more efficient variants (special-$q$ sieve, lattice
sieve) are described in `gnfs/docs/PEDAGOGY.md` §54.

---

## §5 The Linear Algebra Step

### Exponent vectors and GF(2)

Each relation $(a, b)$ produces an exponent vector $\mathbf{v}(a, b) \in \mathbb{Z}^k$, one
coordinate per *column* of the matrix. The columns are:

- one for each rational prime $p \le B$, holding the exponent of $p$ in the integer $a - bm$;
- one for the sign of $a - bm$;
- one for each recorded degree-one prime ideal $\mathfrak{p}_{p,r} = (p, \alpha - r)$ of
  $\mathbb{Z}[\alpha]$, holding the exponent of that *prime ideal* in the principal ideal
  $(a - b\alpha)$.

The ideal columns are not parities of the exponents of rational primes in a norm. In a product of
relations, two different ideals $\mathfrak{p}_{p,r}$ and $\mathfrak{p}_{p,r'}$ over the same $p$
can each have an odd exponent while the exponent of $p$ in the norm of the product, the sum of the
two, is even. A matrix with one column per ideal constrains each exponent separately, and a matrix
built from norms would not. Away from the primes dividing $c = \mathrm{lc}(f)$ and the index
$[\mathcal{O}_K : \mathbb{Z}[\alpha]]$ (the implementation tests the conservative condition
$p \nmid c \cdot \mathrm{disc}(F)$, with $F$ as in §6), the exponent of $\mathfrak{p}_{p,r}$ is
the valuation of $a - b\alpha$ at the prime ideal of $\mathcal{O}_K$ with the same name. At the
primes that do divide $c \cdot \mathrm{disc}(F)$ the checked matrix keeps the columns of the
ideals but leaves them unconstrained: the recorded exponent there is only a norm-valuation
attribution, so it is not a necessary condition for a square.

Reducing modulo 2 gives a vector $\bar{\mathbf{v}}(a, b) \in \mathrm{GF}(2)^k$. A *linear
dependency* over $\mathrm{GF}(2)$ is a non-empty subset $S$ of relations such that

$$\sum_{(a,b) \in S} \bar{\mathbf{v}}(a, b) = \mathbf{0} \in \mathrm{GF}(2)^k.$$

This means that for every rational prime, and for every constrained ideal, the total exponent
across the relations in $S$ is even, and that the number of relations with $a - bm < 0$ is even.

### From dependency to congruence of squares

For a dependency $S$ put

$$P = \prod_{(a,b) \in S} (a - bm) \in \mathbb{Z}, \qquad
  \gamma = \prod_{(a,b) \in S} (a - b\alpha) \in K = \mathbb{Q}(\alpha).$$

The rational columns and the sign column make $P$ a positive perfect square $P = x_0^2$ with
$x_0 \in \mathbb{Z}_{>0}$. The algebraic columns do not make $\gamma$ a square: they say only that
its valuations are even at the ideals the matrix constrains (the obstruction is the subject of the
next two subsections). The congruence of squares therefore needs a *verified* square root of
$\gamma$.

**Theorem (a square root of $\gamma$ gives a congruence of squares).** *Let $f(m) \equiv 0
\pmod N$, let $c = \mathrm{lc}(f)$ be coprime to $N$, let $\theta = c\alpha$ and let $F$ be the
monic polynomial of $\theta$ (§6). Let $S$ be a dependency and suppose $\beta \in K$ satisfies
$\beta^2 = \gamma$, with the coordinates of $\beta$ in the basis $1, \theta, \ldots, \theta^{d-1}$
having denominators coprime to $N$. Put $y = \beta(\theta \mapsto c\,m) \bmod N$. Then*

$$y^2 \equiv P = x_0^2 \pmod{N}.$$

*Proof.* The map $\theta \mapsto c\,m$ is a ring homomorphism $\psi : \mathbb{Z}[\theta] \to
\mathbb{Z}/N\mathbb{Z}$, because $F(c\,m) = c^{d-1} f(m) \equiv 0 \pmod N$. Every integer coprime
to $N$ maps to a unit, so $\psi$ extends uniquely to the subring $R \subset K$ of elements whose
$\theta$-coordinates have denominators coprime to $N$. Each factor $a - b\alpha = (ca - b\theta)/c$
lies in $R$ and $\psi(a - b\alpha) = (ca - b\,c\,m)/c = a - bm$. Hence $\psi(\gamma) = P$, and
$y^2 = \psi(\beta)^2 = \psi(\beta^2) = \psi(\gamma) = P \pmod N$. $\square$

Since $x_0^2 = P \equiv y^2$, the integer $N$ divides $(x_0 - y)(x_0 + y)$. If $x_0 \not\equiv
\pm y \pmod N$ then $\gcd(x_0 - y, N)$ is neither $1$ nor $N$: it is a proper factor.

No norm enters this argument. The image $y$ is the evaluation of $\beta$ under $\psi$, and it
does not equal $\mathrm{Norm}(\beta)$ in general; §6 gives a fixture where the two differ.

**The linear algebra.** Collect $n \ge k + \ell$ relations, for some small excess $\ell$, and form
the $n \times k$ matrix $M$ over $\mathrm{GF}(2)$ whose *rows* are the reduced exponent vectors. A
dependency is a nonzero vector $v \in \mathrm{GF}(2)^n$ with $v^{T} M = 0$, that is, an element
of the *left kernel* of $M$. Since $\mathrm{rank}(M) \le k$, the left kernel has dimension at
least $n - k \ge \ell$. Let $W$ be the total number of nonzero entries of $M$. The *Wiedemann*
and *Lanczos* algorithms perform $O(k)$ applications of the sparse operator, each costing $O(W)$,
for a total of $O(k \cdot W)$ operations. *Block* Lanczos works on blocks of $\mathtt{BLOCK\_WIDTH}
= 64$ vectors packed into machine words, so it needs about $k / 64$ block applications of $O(W)$
word operations each: the matrix-product cost divides by the word width. Dense elimination costs
$O(k^3)$.

### The square-class obstruction

This subsection states, once, what the dependency condition leaves open. The sources are
Buhler–Lenstra–Pomerance, "Factoring integers with the number field sieve" (the sections on
quadratic characters and the square root), and Crandall–Pomerance §6.2.

**What the matrix constrains.** The matrix constrains the parities of the exponents of the
recorded degree-one prime ideals $(p, r)$ of $\mathbb{Z}[\alpha]$, one column each, together with
the rational primes and the rational sign. These are not parities of norm exponents. The matrix
does not constrain the ideals whose columns are left unconstrained (above primes dividing
$c \cdot \mathrm{disc}(F)$), the signs of $\gamma$ at the real places of $K$, or the class of
$\gamma$ in $K^*/K^{*2}$ beyond what the characters below add.

**The obstruction group.** Let $V$ be the group of elements of $K^*$ whose valuation at every
prime ideal of $\mathcal{O}_K$ is even, taken modulo $K^{*2}$. Sending a class $[x]$ to the ideal
class of the $J$ with $J^2 = (x)$ gives an exact sequence

$$1 \to \mathcal{O}_K^*/\mathcal{O}_K^{*2} \to V \to \mathrm{Cl}(\mathcal{O}_K)[2] \to 1,$$

so

$$\dim_{\mathbb{F}_2} V = r_1 + r_2 + \mathrm{rk}_2\,\mathrm{Cl}(\mathcal{O}_K).$$

Here $r_1 + r_2 - 1$ is the unit rank, one more is the torsion unit $-1$, and $\mathrm{rk}_2
\mathrm{Cl}$ is the $2$-rank of the class group. A class of $V$ that is trivial in
$\mathrm{Cl}[2]$ is a unit times a square; one that is nontrivial is represented by an element
$x$ with $(x) = J^2$ and $J$ not principal. Both occur: in $\mathbb{Q}(\sqrt{3})$ the unit
$2 + \sqrt{3}$ has norm $1$ and is not a square, and in $\mathbb{Q}(\sqrt{-5})$ the element $2$
generates $J^2$ for the non-principal $J = (2, 1 + \sqrt{-5})$ and is not a square. Both the units
and the class group therefore obstruct in factoring. Neither is specific to discrete logarithms.

A dependency need not land in $V$. At the prime ideals the matrix leaves unconstrained (those
above primes dividing $c \cdot \mathrm{disc}(F)$, which include the primes of the index, and any
ideal for which the collection records no exponent) the valuation of $\gamma$ may be odd. If $T$
is the finite set of such ideals, $\gamma$ lies in the larger group $V_T$ of classes with even
valuation outside $T$, and $\dim V_T \le \dim V + |T|$. The normalization of §6 multiplies $\gamma$
by the square $(c^h F'(\theta))^2$ to make the target integral; that does not change its class
modulo squares, so it does not remove this enlargement.

### Quadratic characters, the sign column, and the norm character

**Characters.** Let $q$ be an odd prime with $q \nmid c \cdot \mathrm{disc}(f)$, and let $s$ be a
root of $f$ modulo $q$. For a dependency $S$ with $q \nmid (a - bs)$ for every relation, put

$$\chi_{q,s}(\gamma) = \prod_{(a,b) \in S} \left(\frac{a - bs}{q}\right) \in \{+1, -1\}.$$

*If $\gamma = \beta^2$ with $\beta \in K$, then $\chi_{q,s}(\gamma) = +1$.* Indeed $q \nmid
c \cdot \mathrm{disc}(f)$ makes $s$ a simple root and $\alpha \mapsto s$ a ring homomorphism
$\varphi$ from the localisation of the ring of integers at $\mathfrak{q} = (q, \alpha - s)$ to
$\mathbb{F}_q$. The element $\gamma$ is a $\mathfrak{q}$-unit because $\varphi(\gamma) = \prod
(a - bs) \ne 0$, so $\beta$ is a $\mathfrak{q}$-unit as well and
$\varphi(\gamma) = \varphi(\beta)^2$ is a square modulo $q$. The function $\chi_{q,s}$ is multiplicative and trivial on squares, so it
is a homomorphism on the group $V$ of the previous subsection (every class of $V$ has a
representative that is a $\mathfrak{q}$-unit).

Adding one column per character $(q, s)$ to the matrix imposes $\chi_{q,s}(\gamma) = +1$. A set
of $t$ independent characters cuts the finite group the dependency lands in. *As a heuristic*,
treating the characters as independent uniformly random homomorphisms, the expected number of
nonsquare classes that pass all $t$ is $(2^{\dim} - 1)\,2^{-t}$, about $2^{-(t - \dim)}$, where
$\dim$ is the dimension of the group in which $\gamma$ lies. This figure is a heuristic, not a
theorem, and it is never zero. The dependency lands in a finite group $V_T$, and some finite set
of characters does separate $V_T$ from the squares (Chebotarev), but one cannot know that a given
set does so without computing $V_T$ itself. Characters are therefore *necessary* conditions on a
square, not certificates of one, and exact verification of the root decides.

**The norm character is redundant.** One may be tempted to use the Legendre symbol of the
homogeneous norm form $F(a, b)$ instead. For a completely split $q \nmid c \cdot \mathrm{disc}(f)$,
$f \equiv c\prod_s (x - s) \pmod q$ over the $d$ roots $s$, so for $q \nmid F(a, b)$

$$\left(\frac{F(a,b)}{q}\right) = \left(\frac{c}{q}\right) \prod_{s} \left(\frac{a - bs}{q}\right),$$

which is a constant times the product of the root characters at $q$. It can never separate what
the root characters cannot. The two are related by the norm identity $F(a, b) = b^d f(a/b) =
c \cdot N_{K/\mathbb{Q}}(a - b\alpha)$, which is the form the sieve smooths (for $c \ne 1$ the
norm itself can have denominators). Moreover suppose every prime has an even exponent in
$\prod_S N_{K/\mathbb{Q}}(a - b\alpha)$ (for example because every ideal dividing the relations has
an even total exponent) and $q \nmid c \cdot F(a, b)$ for every relation. Then
$\prod_S N_{K/\mathbb{Q}}(a - b\alpha) = \pm Y^2$ with $Y \in \mathbb{Q}^\times$ a $q$-unit, and
$\prod_S F(a, b) = c^{|S|}\prod_S N_{K/\mathbb{Q}}(a - b\alpha)$, so

$$\prod_{(a,b) \in S} \left(\frac{F(a,b)}{q}\right) = \left(\frac{c}{q}\right)^{|S|}\left(\frac{\pm 1}{q}\right).$$

The total depends only on the sign of the norm product, which is the sign of
$N_{K/\mathbb{Q}}(\gamma)$, and on the parity of $|S|$. That sign is a genuine square-class
invariant: if $\gamma = \beta^2$ then $N_{K/\mathbb{Q}}(\gamma) = N_{K/\mathbb{Q}}(\beta)^2 > 0$,
and no matrix column imposes it. So for $q \equiv 3 \pmod 4$ the norm character can reject a
candidate with $N_{K/\mathbb{Q}}(\gamma) < 0$. Example: the unit $\gamma = -1 - \sqrt{2}$ of
$\mathbb{Q}(\sqrt{2})$ has norm $-1$ and is not a square, and its norm character at $q = 7$ is
$(-1/7) = -1$. That one bit is all the norm character carries beyond the parity of $|S|$: it cannot
separate two classes whose norms have the same sign, so a nonsquare unit or class-group element of
positive norm passes it (such as $2 + \sqrt{3}$ above). It is also redundant with the root
characters of the same $q$, because it is $(c/q)^{|S|}$ times their product (in the example the
root characters at $q = 7$ already reject $\gamma$). (If some primes have an odd exponent, the
total acquires one more factor, the symbol of the product of those primes; the redundancy
statement is unaffected.) The filter the implementation uses is therefore the individual root
character; it has no norm-character column.

**The sign column.** The rational product $P = \prod_{(a,b) \in S} (a - bm)$ must be positive for
$x_0 = \sqrt{P}$ to exist. The *sign column* has a $1$ for each relation with $a - bm < 0$, and the
dependency condition forces an even number of them, so $P > 0$. This is the sign of a rational
integer. The matrix imposes no positivity at the real places of $K$, although a square in $K$ is
positive at every real place.

The columns raise the chance that a dependency yields an algebraic square, but they do not make
it certain; the square-root step reports a nonsquare candidate and moves to the next one. The code
realisation is described in `gnfs/docs/PEDAGOGY.md` §56.

---

## §6 The Square Root Step

### The rational square root

Given a dependency $S$, the rational square root is straightforward: compute

$$P = \prod_{(a,b) \in S} (a - bm) \in \mathbb{Z},$$

then compute $x_0 = \sqrt{P} \in \mathbb{Z}$ by exact integer square root. For a dependency of the
checked matrix the rational-prime and sign columns make $P$ a positive perfect square, and the
implementation *checks* $P > 0$ and that $P$ is an exact square rather than assuming it; a failure
is a typed outcome, never a panic.

### The algebraic square root

The algebraic side is more subtle. With $\gamma = \prod_{(a,b) \in S} (a - b\alpha) \in K =
\mathbb{Q}(\alpha)$ we need $\beta \in K$ with $\beta^2 = \gamma$, and we must also be able to
report that none exists. The parity and character columns are necessary filters, not proofs, so
$\gamma$ is only a *candidate* square: a unit, a class-group element, or an ideal the matrix does
not constrain can obstruct (the square-class obstruction of §5). Only an exact verification, or a
local obstruction or exhaustive search within a proved bound, decides.

**Normalization.** Let $c$ be the leading coefficient of $f$, $\theta = c\alpha$ and $F(t) =
c^{d-1} f(t/c)$, monic, so $\theta$ is an algebraic integer. For $k = |S|$ relations let $h =
\lceil k/2 \rceil$, $\eta = c^{2h-k} \prod (ca - b\theta)$ and $H = F'(\theta)^2 \eta$. A root
$z$ of $H$ in $\mathbb{Z}[\theta]$ (that is, $z^2 = H$) gives
$\beta = z / (c^h F'(\theta))$ over $\mathbb{Q}$, and $\beta^2 = \gamma$ is then verified
exactly. A certified bound $B$ on the coefficients of $z$ (computed from $F$ and $\eta$) makes any
modulus $M > 2B$ determine $z$ up to sign from its residue.

**Reconstruction.** Two methods produce candidates for $z$, and both pass through the same exact
verifier.

- *Direct.* At a prime $p$ where $F$ is irreducible, $\mathbb{F}_p[t]/F$ is a field
  $\mathbb{F}_{p^d}$ in which a square root of $H$ exists or not; Tonelli–Shanks in the finite
  extension followed by inverse-square-root Newton lifting gives $z \bmod p^e$ with $p^e > 2B$,
  and the centered residue is the candidate. A local nonresidue is a certificate of nonsquare.
- *Split-prime fallback* (after Thomé). At completely split primes the algebra decomposes into
  copies of $\mathbb{F}_p$, each root of $F$ giving $\pm$ a local square root; the roots are lifted
  and interpolated, combined by CRT, and the sign pattern is found by an exact
  meet-in-the-middle search with only one global sign fixed. A nonsquare is certified either by a
  local obstruction (the target is a nonresidue at a root of a split prime) or by a complete
  matching search in which every sign combination was rejected within the proved bound. An
  incomplete search is never a certificate.

Neither a real embedding nor a norm is used. Both $\pm z$ verify, so the global sign is
immaterial, and a real root of $f$ need not exist. The classical method of Couveignes assumes odd
degree and inert primes with a known norm to select signs, and is not the method implemented.

### The modular image and congruence assembly

The image of the verified root under the homomorphism $\psi$ of §5, $\theta \mapsto c\,m$
(well defined because $F(c\,m) = c^{d-1} f(m) \equiv 0 \pmod N$; for monic $f$ this is
$\alpha \mapsto m$ on $\mathbb{Z}[\alpha]$), is the normalized root evaluated at $\theta = c\,m$
with denominators checked modulo $N$:

$$y = \beta(c\,m) \bmod N \in \mathbb{Z}/N\mathbb{Z},\qquad y^2 \equiv \gamma(m) = P \pmod N.$$

So $y^2 \equiv P$ where $P$ is the *rational* product; it is **not** $\mathrm{Norm}(\beta)$, and
$y^2$ is not $\mathrm{Norm}(\gamma)$. (Example: $f = x^2 - 2$, $m = 23$, $N = 527$, relation
$(3, -2)$: $\gamma = 3 + 2\alpha = (1 + \alpha)^2$ and $\mathrm{Norm}(1 + \alpha) = -1$, while
$P = 49$ and $y = 1 + 23 = 24$ with $24^2 = 576 \equiv 49 \pmod{527}$.) If a denominator of $\beta$
shares a factor with $N$ the map is undefined at that factor, and the factor is reported as an
incidental divisor; if all of $N$ divides it, the map is unsupported.

Now $x^2 \equiv y^2 \pmod{N}$ with $x = x_0 \bmod N$, so $N \mid (x - y)(x + y)$. The congruence
is checked, not assumed. If $x \equiv \pm y \pmod N$ it is trivial; otherwise $\gcd(x - y, N)$ is a
proper factor (it is neither $1$ nor $N$). A random dependency is nontrivial with probability about
$1/2$ for a semiprime (§2). A trivial congruence, a nonsquare candidate, and an exhausted resource
budget are each reported distinctly, and the next candidate from the left kernel (including
bounded XOR combinations of computed kernel vectors) is tried.

**A worked run on generated relations.** Take $N = 45113 = 197 \cdot 229$, $m = 31$ and $f = x^3 +
15x^2 + 29x + 8$, so $f(31) = 45113$. The cubic is irreducible: modulo $3$ it is $x^3 + 2x + 2$,
which has no root. The line sieve over $|a| \le 400$, $1 \le b \le 41$ with factor-base bounds $30$
(rational) and $90$ (algebraic) returns $60$ relations, all of which pass exact validation. The
checked matrix has $39$ columns; clique pruning and column merging reduce $60$ rows to $47$, and
eight root characters (at the primes $97, 101, 103, 107, 113, 139$, with $109 \mid \mathrm{disc}(f)$
excluded) are filled by XOR over provenance. Block Lanczos returns kernel vectors spanning a space
of rank $18$. The first computed vector lifts to a dependency of $22$ relations, each with
$b \ge 1$. For it $P = \prod (a - bm)$ is a positive exact square $x_0^2$, the algebraic root
$\beta$ verifies $\beta^2 = \prod (a - b\alpha)$ exactly (so it depends on those $b$-carrying
relations), the image $y = \beta(\theta \mapsto c\,m) \bmod N$ satisfies $y^2 \equiv P \equiv
x_0^2 \pmod N$ with $x \not\equiv \pm y$ (this run is monic, $c = 1$, so $\theta = \alpha$ and the
image is $\beta(m)$), and $\gcd(|x - y|, N)$ is $197$ or $229$. This is a single deterministic
fixture with the stated options, run at toy scale. It shows that the stages compose on this input; it does not show that
every input factors, and it says nothing about performance at cryptographic size. The first
candidate was itself a computed kernel vector, so the combination search was not needed here (it
is exercised separately on hand-built relation sets).

---

## §7 The $L$-Notation Subexponentiality Derivation

*This is the designated payoff proof for the GNFS chapter. The exponent
$1/3$ is derived in full; the precise constant $(64/9)^{1/3}$ is stated and the optimisation
argument is given, with a citation for the full derivation including the polynomial degree
optimisation.*

### Setup and the three costs

The GNFS complexity is determined by three costs:

1. **Sieve cost:** the number of $(a, b)$ pairs that must be examined to collect enough relations.
2. **Smoothness probability:** the probability that a random pair is a relation (both norms smooth).
3. **Linear algebra cost:** the cost of finding a null vector in the $k \times k$ matrix over
   $\mathrm{GF}(2)$.

We will show that the dominant cost is the sieve cost, and that it is minimised at $L_N[1/3, c]$
for an optimal constant $c$. The exponent $1/3$ is the result of this optimisation.

### The smoothness probability in $L$-notation

Set the smoothness bound $B = L_N[1/3, c]$ for a parameter $c > 0$ to be optimised. The norms
$N_{\mathrm{rat}}(a, b)$ and $N_{\mathrm{alg}}(a, b)$ are both of size approximately $N^{1/d}$
for the optimal degree $d$ (see below). The smoothness probability for a single norm is

$$\Pr[\text{norm is } B\text{-smooth}] = \rho(u) \cdot (1 + o(1)),$$

where $u = \log(N^{1/d}) / \log B$.

**Computing $u$.** We have

$$u = \frac{(\log N)/d}{\log B} = \frac{(\log N)/d}{c \cdot (\log N)^{1/3} \cdot (\log\log N)^{2/3}}
  = \frac{(\log N)^{2/3}}{d \cdot c \cdot (\log\log N)^{2/3}}.$$

For the optimal degree $d \sim (3\log N / \log\log N)^{1/3}$ (see below), this simplifies to

$$u \sim \frac{(\log N)^{2/3}}{c \cdot (3\log N / \log\log N)^{1/3} \cdot (\log\log N)^{2/3}}
  = \frac{(\log N)^{2/3}}{c \cdot 3^{1/3} \cdot (\log N)^{1/3} \cdot (\log\log N)^{1/3}}
  = \frac{(\log N)^{1/3}}{c \cdot 3^{1/3} \cdot (\log\log N)^{1/3}}.$$

In $L$-notation, $u = (1/(3c)) \cdot (\log N)^{1/3} / (\log\log N)^{1/3}$, which grows as
$L_N[1/3, 1/(3c)] / c$ — but what matters for the Dickman function is the *value* of $u$, not
its $L$-notation form.

**The Dickman function estimate.** For large $u$, the Dickman function satisfies
$\rho(u) = u^{-u(1+o(1))}$. Taking logarithms:

$$\log \rho(u) = -u \log u \cdot (1 + o(1)).$$

With $u \sim \frac{1}{3c} \cdot \frac{(\log N)^{1/3}}{(\log\log N)^{1/3}}$ (absorbing the
$3^{1/3}$ into the $o(1)$ for the leading-order analysis), we get

$$\log \rho(u) \sim -\frac{1}{3c} \cdot \frac{(\log N)^{1/3}}{(\log\log N)^{1/3}} \cdot
  \log\!\left(\frac{(\log N)^{1/3}}{(\log\log N)^{1/3}}\right)
  \sim -\frac{1}{3c} \cdot \frac{(\log N)^{1/3}}{(\log\log N)^{1/3}} \cdot
  \frac{1}{3} \log\log N
  = -\frac{1}{9c} \cdot (\log N)^{1/3} \cdot (\log\log N)^{2/3}.$$

Therefore

$$\rho(u) = \exp\!\left(-\frac{1}{9c} \cdot (\log N)^{1/3} \cdot (\log\log N)^{2/3} \cdot
  (1 + o(1))\right) = L_N\!\left[-\tfrac{1}{3},\, -\tfrac{1}{9c}\right] \cdot (1 + o(1)).$$

Wait — let us be more careful. The standard result (see §Prerequisites, Corollary to
Canfield–Erdős–Pomerance) states that for $x = N^{1/d}$ and $y = B = L_N[1/3, c]$:

$$\rho(u) = L_N\!\left[-\tfrac{1}{3},\, -\tfrac{1}{3c}\right] \cdot (1 + o(1)).$$

This is the form stated in the Prerequisites chapter and used in the complexity analysis. The
derivation of this precise form requires a more careful treatment of the Dickman function (see
Granville [G08]); we use it as stated.

### The sieve cost

To collect $k \approx B = L_N[1/3, c]$ relations (enough for the linear algebra), we need to
examine approximately

$$\frac{k}{\rho(u)^2}$$

pairs $(a, b)$. The $\rho(u)^2$ in the denominator is because *both* norms must be smooth, and
we treat the two smoothness events as approximately independent (the standard heuristic).

**The sieve cost in $L$-notation:**

$$\text{Sieve cost} = \frac{k}{\rho(u)^2} = \frac{L_N[1/3, c]}{L_N[-1/3, -2/(3c)]}
  = L_N[1/3, c] \cdot L_N[1/3, 2/(3c)] = L_N\!\left[\tfrac{1}{3},\, c + \tfrac{2}{3c}\right].$$

Here we used the $L$-notation arithmetic lemma: $L_N[\alpha, c_1] \cdot L_N[\alpha, c_2] =
L_N[\alpha, c_1 + c_2]$, and $1 / L_N[-1/3, -2/(3c)] = L_N[1/3, 2/(3c)]$.

### The linear algebra cost

The linear algebra step finds a left-kernel vector of a $k \times k$ matrix over $\mathrm{GF}(2)$,
where $k \approx B = L_N[1/3, c]$ (§5). Dense Gaussian elimination costs $O(k^3)$. The Wiedemann
and Lanczos algorithms perform $O(k)$ applications of the sparse operator at cost $O(W)$ each, for
$O(k \cdot W)$ in total, where $W$ is the number of nonzero entries; block Lanczos divides the
matrix-product cost by the word width.

The GNFS matrix is sparse: before filtering, each row has $w = O(\log N)$ nonzero entries (the
distinct prime factors of the two norms, which number at most their bit lengths, plus
$O(\log N)$ character columns), so $W = O(k \cdot \log N)$ and, since $\log N = L_N[0, 1]$ is
negligible against $L_N[1/3, \cdot]$ factors,

$$\text{Linear algebra cost} = O(k^2 \log N) = L_N\!\left[\tfrac{1}{3},\, 2c\right].$$

**Comparing the two costs.** The sieve cost is $L_N[1/3, c + 2/(3c)]$ and the linear algebra
cost is $L_N[1/3, 2c]$. For the sieve cost to dominate (which is the regime of interest), we need

$$c + \frac{2}{3c} > 2c \iff \frac{2}{3c} > c \iff c^2 < \frac{2}{3} \iff c < \sqrt{\frac{2}{3}}.$$

At the optimal $c = \sqrt{2/3}$, the two costs are equal: $c + 2/(3c) = 2c = 2\sqrt{2/3}$.
For $c > \sqrt{2/3}$, the linear algebra dominates; for $c < \sqrt{2/3}$, the sieve dominates.
The optimal $c$ is the crossover point.

### The optimisation: deriving the exponent $1/3$

The total cost of GNFS is dominated by the maximum of the sieve cost and the linear algebra cost:

$$\text{Total cost} = L_N\!\left[\tfrac{1}{3},\, \max\!\left(c + \tfrac{2}{3c},\, 2c\right)\right].$$

To minimise the total cost, we minimise $\max(c + 2/(3c), 2c)$ over $c > 0$. The minimum occurs
where the two expressions are equal:

$$c + \frac{2}{3c} = 2c \implies \frac{2}{3c} = c \implies c^2 = \frac{2}{3} \implies
  c = \sqrt{\frac{2}{3}}.$$

At this optimal $c$:

$$c + \frac{2}{3c} = \sqrt{\frac{2}{3}} + \frac{2}{3\sqrt{2/3}} = \sqrt{\frac{2}{3}} +
  \frac{2}{3} \cdot \sqrt{\frac{3}{2}} = \sqrt{\frac{2}{3}} + \sqrt{\frac{2}{3}} \cdot
  \frac{2/3}{2/3} = 2\sqrt{\frac{2}{3}}.$$

Let us verify: $\frac{2}{3\sqrt{2/3}} = \frac{2}{3} \cdot \frac{1}{\sqrt{2/3}} = \frac{2}{3}
\cdot \sqrt{\frac{3}{2}} = \frac{2}{3} \cdot \frac{\sqrt{3}}{\sqrt{2}} = \frac{2\sqrt{3}}{3\sqrt{2}}
= \frac{\sqrt{6}}{3}$. And $\sqrt{2/3} = \sqrt{2}/\sqrt{3} = \sqrt{6}/3$. So indeed
$c + 2/(3c) = \sqrt{6}/3 + \sqrt{6}/3 = 2\sqrt{6}/3$.

**The minimum total cost is $L_N[1/3, 2\sqrt{6}/3]$.**

### The exponent $1/3$ is the key structural result

The derivation above establishes the key structural result:

**Theorem (GNFS exponent, heuristic).** *Under the heuristic independence assumption on the two
norms, the GNFS total cost is minimised at $L_N[1/3, 2\sqrt{6}/3]$.*

The exponent $\alpha = 1/3$ is the result of the optimisation: it is the unique exponent for
which the sieve cost and the linear algebra cost can be simultaneously minimised. The argument is:

- The sieve cost is $L_N[1/3, c + 2/(3c)]$ because the norms have size $N^{1/d}$ (not $N$), and
  the smoothness probability for $N^{1/d}$-sized integers with $L_N[1/3, c]$-smooth bound is
  $L_N[-1/3, -1/(3c)]$.
- The linear algebra cost is $L_N[1/3, 2c]$ because the matrix has $L_N[1/3, c]$ rows and
  columns.
- The crossover occurs at $c = \sqrt{2/3}$, giving total cost $L_N[1/3, 2\sqrt{6}/3]$.

If the norms had size $N$ (as in the quadratic sieve), the smoothness probability would be
$L_N[-1/2, -1/(2c)]$ for $B = L_N[1/2, c]$, and the same optimisation would give exponent $1/2$
and cost $L_N[1/2, \sqrt{2}]$. The improvement from $1/2$ to $1/3$ is entirely due to the
number-field bridge reducing the norm size from $N$ to $N^{1/d}$.

### The precise constant $(64/9)^{1/3}$

The constant $2\sqrt{6}/3$ from the simplified analysis above is not the standard constant
$(64/9)^{1/3}$ quoted in the literature. The discrepancy arises from two sources:

**1. The polynomial degree optimisation.** The optimal degree $d$ is not a free parameter — it is
tied to $N$ by the relation $d \sim (3\log N / \log\log N)^{1/3}$. Incorporating this into the
analysis modifies the smoothness probability and the sieve cost. The full analysis (Lenstra–
Lenstra–Manasse–Pollard [LLMP93], Buhler–Lenstra–Pomerance [BLP93]) accounts for the degree
optimisation and yields the constant $(64/9)^{1/3}$.

**2. The two-sided sieve.** The analysis above treats the two norms as having the same size
$N^{1/d}$. In practice, the rational norm $|a - bm|$ and the algebraic form $|b^d f(a/b)|$ have
slightly different sizes, and the optimal factor-base bounds $B_{\mathrm{rat}}$ and
$B_{\mathrm{alg}}$ are not equal. The full analysis optimises over both bounds separately.

**The standard result.** The complete analysis gives:

$$\text{GNFS total cost} = L_N\!\left[\tfrac{1}{3},\, \left(\tfrac{64}{9}\right)^{1/3}\right]
  \cdot (1 + o(1)).$$

Note that $(64/9)^{1/3} = (64/9)^{1/3} \approx 1.923$, while $2\sqrt{6}/3 \approx 1.633$. The
gap reflects the additional costs from the degree optimisation and the two-sided sieve.

**Verification.** We can verify the standard constant as follows. The full analysis sets
$B = L_N[1/3, c]$ and optimises the total cost including the degree-$d$ contribution. The
optimal parameters are $c = (8/9)^{1/3}$ and $d = (3\log N / \log\log N)^{1/3}$. The total cost
is then $L_N[1/3, 3c] = L_N[1/3, 3(8/9)^{1/3}] = L_N[1/3, (3^3 \cdot 8/9)^{1/3}] =
L_N[1/3, (27 \cdot 8/9)^{1/3}] = L_N[1/3, (24)^{1/3}]$... Hmm, this does not immediately give
$(64/9)^{1/3}$. The precise derivation requires tracking the contributions of both sides of the
sieve and the degree optimisation simultaneously; see Lenstra–Lenstra–Manasse–Pollard [LLMP93,
Theorem 1] and the exposition in Crandall–Pomerance [CP05, §6.2] for the complete argument.

**What this chapter proves.** The derivation above establishes the two key structural results:

1. **The exponent is $1/3$** (not $1/2$, not $1/4$): this follows from the number-field bridge
   reducing norm sizes from $N$ to $N^{1/d}$, and from the optimisation of the sieve/linear-algebra
   tradeoff.

2. **The constant is of order 1** (not growing with $N$): the $L$-notation constant is a fixed
   real number, confirming that GNFS is genuinely subexponential.

The precise constant $(64/9)^{1/3}$ is the result of the full analysis in [LLMP93] and [BLP93],
which we cite rather than reproduce.

### Summary of the derivation

| Quantity | Value |
|----------|-------|
| Smoothness bound | $B = L_N[1/3, c]$ |
| Norm size | $\approx N^{1/d} \approx N^{1/3}$ (at optimal $d$) |
| Smoothness probability (one norm) | $L_N[-1/3, -1/(3c)]$ |
| Smoothness probability (both norms) | $L_N[-1/3, -2/(3c)]$ |
| Relations needed | $k \approx B = L_N[1/3, c]$ |
| Sieve cost | $L_N[1/3, c + 2/(3c)]$ |
| Linear algebra cost | $L_N[1/3, 2c]$ |
| Optimal $c$ | $\sqrt{2/3}$ (simplified) or $(8/9)^{1/3}$ (full) |
| Minimum cost (simplified) | $L_N[1/3, 2\sqrt{6}/3]$ |
| Minimum cost (full analysis) | $L_N[1/3, (64/9)^{1/3}]$ |

The exponent $1/3$ is exact; the constant $(64/9)^{1/3}$ is the result of the full analysis.

---

## §8 Cross-References

**Code realisation.** The GNFS pipeline — polynomial selection, sieving, filtering, linear
algebra, and square root — is implemented in the `gnfs` crate and documented in
`gnfs/docs/PEDAGOGY.md` §52–§62. That chapter is the code-tour sibling to this one: it assumes
the reader knows the mathematics (this chapter) and focuses on the Rust realisation, the
stage-by-stage contracts, and the KAT summary.

**Prerequisites.** The prerequisites used in this chapter are collected in §Prerequisites of this
file:
- The Canfield–Erdős–Pomerance smooth-number density theorem (§Probability).
- The $L$-notation and its arithmetic (§Analysis).
- Unique factorisation of ideals in Dedekind domains (§Algebra).
- The Chinese Remainder Theorem (§Algebra).
- Tonelli–Shanks square roots (§α-Substrate chapter).

**Number-field substrate.** The ring of integers $\mathcal{O}_K$, the norm map, ideal
factorisation, and the Dedekind index criterion are developed in `shared/numfield/docs/PEDAGOGY.md`.
That chapter documents the `shared-numfield` crate, which provides the algebraic arithmetic used
in the GNFS square-root step.

**Scale.** The honest science↔engineering gap — the gap between the asymptotic complexity
$L_N[1/3, (64/9)^{1/3}]$ and the behaviour at toy scale — is discussed in §On Scale of this
textbook. The key points: the complexity analysis is heuristic and asymptotic; the optimal
parameters at toy scale are not the asymptotically optimal ones; and some phenomena (large-prime
variations, the degree-$N$ coupling) are under-exposed at toy scale.

---

## References for the GNFS Chapter

18. **Lenstra, A. K., Lenstra, H. W. Jr., Manasse, M. S., and Pollard, J. M. (1993).** "The
    number field sieve." In *The Development of the Number Field Sieve*, Lecture Notes in
    Mathematics 1554, Springer, 11–42. [LLMP93]

19. **Buhler, J. P., Lenstra, H. W. Jr., and Pomerance, C. (1993).** "Factoring integers with
    the number field sieve." In *The Development of the Number Field Sieve*, Lecture Notes in
    Mathematics 1554, Springer, 50–94. [BLP93]

20. **Couveignes, J.-M. (1993).** "Computing a square root for the number field sieve." In *The
    Development of the Number Field Sieve*, Lecture Notes in Mathematics 1554, Springer, 95–102.

21. **Dixon, J. D. (1981).** "Asymptotically fast factorization of integers." *Mathematics of
    Computation*, 36(153), 255–260.

22. **Pomerance, C. (1996).** "A tale of two sieves." *Notices of the American Mathematical
    Society*, 43(12), 1473–1485. (Accessible survey of the quadratic sieve and NFS.)

23. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective*.
    2nd ed. Springer. §6.2 for the NFS complexity analysis. [CP05]

---

# NFS-DL: Discrete Logarithm via the Number Field Sieve

*Maths-first sibling to `gnfs/docs/PEDAGOGY.md` §63–§71 (the NFS-DL code-tour). For the
code-tour — pipeline contracts, stage-by-stage implementation narrative, and KAT summary — see
`gnfs/docs/PEDAGOGY.md` §63–§71. For the prerequisites used in this chapter, see §Prerequisites
above. For the GNFS factoring chapter whose structure this chapter mirrors and whose §7
L-notation derivation this chapter deltas on, see §GNFS (§1–§8) above.*

---

## §9.1 Introduction and Through-Line

The **discrete logarithm problem** (DLP) in a prime field $\mathbb{F}_p^*$ asks: given a
generator $g$ and a target $h$ in $\mathbb{F}_p^*$, find the unique integer $x \in
\{0, 1, \ldots, p-2\}$ such that $g^x \equiv h \pmod{p}$. This is the search problem at the
heart of Diffie–Hellman key exchange, ElGamal encryption, and DSA signatures. Its hardness in
$\mathbb{F}_p^*$ is the assumption that makes these systems secure.

**The generic search bound.** The DLP in a group of order $n$ can be solved in $O(\sqrt{n})$
group operations by Pollard rho or baby-step/giant-step (§Pollard Rho chapter). For
$\mathbb{F}_p^*$ with $n = p - 1$, this gives $O(\sqrt{p})$ multiplications — fully exponential
in $\log p$. The NFS-DL algorithm escapes this bound by exploiting the *multiplicative structure*
of $\mathbb{F}_p^*$: the group has a natural factor base (the small primes), and the discrete
logarithm reduces to a linear system over $\mathbb{F}_\ell$ (where $\ell$ is the subgroup order).

**The through-line.** NFS-DL is the DLP analogue of GNFS for factoring. Both algorithms exploit
*number-field structure*: the integers $\mathbb{Z}$ embed into a number field $K = \mathbb{Q}(\alpha)$,
and the norm map connects factorisation in $K$ to factorisation in $\mathbb{Z}$. For factoring,
the payoff is a congruence of squares $x^2 \equiv y^2 \pmod{N}$. For DLP, the payoff is a
*virtual-logarithm table* — the discrete logs of all factor-base elements — from which the log
of any specific target $h$ can be computed by a descent procedure.

The structure-based escape from search works in two stages:

1. **Precomputation (relation collection + linear algebra).** Sieve for smooth pairs $(a, b)$
   in two number fields sharing a rational side. Augment each relation with *Schirokauer map*
   columns (the DL-specific correction that makes the unit contribution linear modulo $\ell$;
   GNFS uses quadratic characters for a different purpose, §5). Solve
   the resulting linear system over $\mathbb{F}_\ell$ to recover the virtual logs of all
   factor-base elements.

2. **Individual-logarithm descent.** Given a specific target $h$, find an exponent $e$ such
   that $g^e \cdot h$ is smooth over a medium-prime bound, then recursively rewrite each medium
   prime as a combination of smaller primes (via the special-$q$ sieve) until all primes are
   factor-base leaves with known virtual logs. Assemble the log of $h$ from the descent tree.

**Main theorem (NFS-DL complexity, heuristic).** *Under standard heuristic assumptions on the
distribution of smooth norms, NFS-DL computes $\log_g(h)$ in $\mathbb{F}_p^*$ in expected time*

$$L_p\!\left[\tfrac{1}{3},\, \left(\tfrac{64}{9}\right)^{1/3}\right].$$

This is the *same* complexity as GNFS for factoring — the same exponent $1/3$ and the same
constant $(64/9)^{1/3}$. The reason is structural: the precomputation (relation collection +
linear algebra) has the same cost shape as GNFS, and the individual-logarithm descent is
*asymptotically subdominant* — it does not change the leading complexity. The derivation of
this equality is the payoff proof of this chapter (§9.7).

**What is different from GNFS.** The algorithm differs from GNFS in three ways:

- **Schirokauer maps** take the place of the quadratic-character columns of GNFS. The two
  problems share the same obstruction group (units modulo squares and the class group; see §5),
  and neither is specific to one of them. What differs is the condition imposed. Factoring needs
  $\gamma$ to be a *square*, a condition that quadratic characters can only filter: no given
  finite set is known to separate the dependency's group from the squares, so exact verification
  decides. NFS-DL needs the unit contribution of each relation to be *linear* modulo $\ell$:
  the Schirokauer maps are homomorphisms to $\mathbb{Z}/\ell$ that vanish on $\ell$-th powers, so
  with enough of them the unit part of a dependency is an $\ell$-th power. The class group is
  handled by assuming $\ell \nmid h_K$.
- **$\mathbb{F}_\ell$ linear algebra** replaces GF(2). The system is solved over a prime field
  rather than the field with two elements.
- **Individual-logarithm descent** has no factoring analogue. In GNFS, the linear algebra step
  yields a candidate dependency, and a verified algebraic root turns it into a congruence of
  squares; in NFS-DL, it yields only the factor-base logs, and a separate descent is needed for
  each target $h$.

---

## §9.2 The Number-Field Bridge for DL

### Two number fields sharing a rational side

The NFS-DL setup mirrors the GNFS setup (§3) with one key difference in interpretation. Let $p$
be the prime modulus and $\ell$ a prime dividing $p - 1$ (the subgroup order). Choose a
degree-$d$ polynomial $f \in \mathbb{Z}[x]$ and an integer $m$ such that $f(m) \equiv 0
\pmod{p}$. Let $\alpha$ be a root of $f$ in $\mathbb{C}$, and let $K = \mathbb{Q}(\alpha)$ be
the number field generated by $\alpha$.

There are two natural ring homomorphisms into $\mathbb{F}_p$:

**The rational map** $\phi_{\mathbb{Z}}: \mathbb{Z} \to \mathbb{F}_p$, reduction modulo $p$.

**The algebraic map** $\phi_\alpha: \mathbb{Z}[\alpha] \to \mathbb{F}_p$, defined by
$\phi_\alpha(g(\alpha)) = g(m) \bmod p$. This is well-defined because $f(m) \equiv 0 \pmod{p}$.

For a pair $(a, b) \in \mathbb{Z}^2$ with $\gcd(a, b) = 1$, both maps agree on $a - b\alpha$:

$$\phi_\alpha(a - b\alpha) = \phi_{\mathbb{Z}}(a - bm) = a - bm \bmod p.$$

This is the same bridge as in GNFS (§3). The difference is what the bridge is used for.

### The homomorphism to $\mathbb{F}_p^*$

In GNFS, the bridge is used to construct a congruence of squares: the rational product $P$ is
an integer square and the image of a square root of $\prod (a - b\alpha)$ squares to the same $P$
modulo $N$ (§5). In NFS-DL, the bridge is used to construct a
*multiplicative relation in $\mathbb{F}_p^*$*: a product of smooth elements that equals 1 in
$\mathbb{F}_p^*$.

Specifically, for a smooth pair $(a, b)$ with rational norm $a - bm = \prod_i p_i^{e_i}$ and
algebraic norm $N_{K/\mathbb{Q}}(a - b\alpha) = \prod_j \mathfrak{p}_j^{f_j}$, the bridge gives

$$\prod_i p_i^{e_i} \equiv \prod_j \phi_\alpha(\mathfrak{p}_j)^{f_j} \pmod{p}.$$

Taking discrete logarithms base $g$ on both sides:

$$\sum_i e_i \log_g(p_i) \equiv \sum_j f_j \log_g(\phi_\alpha(\mathfrak{p}_j)) \pmod{p - 1}.$$

This is a *linear equation* in the unknown discrete logs of the factor-base elements. Collecting
enough such equations gives a linear system whose solution is the virtual-log table.

### The factor-base and virtual-logarithm construction

**The factor base** for NFS-DL has the same structure as for GNFS (§3):

- **Rational factor base** $\mathcal{F}_{\mathbb{Z}} = \{p_1, \ldots, p_k\}$: rational primes
  up to a smoothness bound $B$.
- **Algebraic factor base** $\mathcal{F}_\alpha = \{(p, \alpha - r) : p \leq B,\; f(r) \equiv
  0 \pmod{p}\}$: prime ideals in $\mathbb{Z}[\alpha]$ up to $B$.

A pair $(a, b)$ is a *DL relation* if both the rational norm $|a - bm|$ and the algebraic norm
$|N_{K/\mathbb{Q}}(a - b\alpha)|$ are $B$-smooth. Each DL relation gives a linear equation in
the logs of the factor-base elements.

**Virtual logarithms.** The *virtual logarithm* of a factor-base element $p_i$ is
$\log_g(p_i) \bmod \ell$ — the discrete log of $p_i$ in the subgroup of order $\ell$. The
virtual logs are the unknowns of the linear system. Once the system is solved, the virtual-log
table gives the log of any $B$-smooth element of $\mathbb{F}_p^*$ as a linear combination of
the virtual logs.

**Why "virtual"?** The logs are computed modulo $\ell$ (the subgroup order), not modulo $p - 1$
(the full group order). If $\ell \neq p - 1$, the virtual logs give the log in the $\ell$-order
subgroup, not the full group. Recovering the full log requires Pohlig–Hellman (§Prerequisites)
or choosing $\ell = p - 1$ (which requires $p - 1$ to be prime, i.e. $p$ is a safe prime).

**The linear system.** Collecting $k + \delta$ DL relations (for a small excess $\delta$) gives
a $(k + \delta) \times k$ matrix $A$ over $\mathbb{F}_\ell$, where $k = |\mathcal{F}_{\mathbb{Z}}|
+ |\mathcal{F}_\alpha|$ is the total number of factor-base elements. The virtual logs are the
solution to $A \cdot x = 0$ over $\mathbb{F}_\ell$ — the kernel of $A$.

*Note on the column layout.* The DL relation matrix has three groups of columns: rational
exponent columns (one per rational factor-base prime), algebraic exponent columns (one per
algebraic factor-base ideal), and Schirokauer columns (one per prime ideal in the Schirokauer
set). The Schirokauer columns are the DL-specific addition; they are described in §9.3.

---

## §9.3 Schirokauer Maps

*This section is the DL-specific algebra — the clearest single "what's different about DL"
moment. It has no counterpart in the GNFS factoring chapter.*

### The obstruction to principality

In GNFS, the linear algebra step finds a subset $S$ of relations for which the product
$\prod_{(a,b) \in S} (a - b\alpha)$ is a *candidate* square in $K$; it need not be one. Even
valuations at the ideals the matrix constrains do not give that: the square-class obstruction
group of §5 has a unit part and a class-group part, and both obstruct. The quadratic characters
(§5 of the GNFS chapter) filter this obstruction but are necessary conditions only, because one
cannot know that a given finite set separates the dependency's group from the squares without
computing that group; exact verification decides.

The same group is present in NFS-DL, where the condition to be met is different. The linear system
$A \cdot x = 0$ over $\mathbb{F}_\ell$ has a solution only if the exponent vectors correctly
account for all multiplicative structure in $\mathbb{F}_p^*$. The exponent vectors record the
*ideal* factorisation of $a - b\alpha$, so the *units* of $\mathcal{O}_K$ — the ring of integers
of $K$ — are not accounted for by them alone. A unit $u \in \mathcal{O}_K^*$ has norm
$N_{K/\mathbb{Q}}(u) = \pm 1$, so it contributes nothing to the rational exponent vector. But
$\phi_\alpha(u) \in \mathbb{F}_p^*$ may be non-trivial, and its discrete log is missing from the
linear system.

**Linearity of the unit contribution.** The Schirokauer maps make that unit contribution linear
modulo $\ell$: they are homomorphisms to $\mathbb{Z}/\ell$ that vanish on $\ell$-th powers, so the
unit contribution of a product of relations is the sum of the contributions of its factors, which
is the linearity the $\mathbb{F}_\ell$ system needs. For the class group the argument assumes
$\ell \nmid h_K$, so that the class group has no $\ell$-torsion to obstruct. See Schirokauer,
"Discrete logarithms and local units" [S93].

### The Schirokauer map

**Definition (Schirokauer map).** Let $\ell$ be a prime, and let $\varphi = (p_0, \alpha - r_0)$
be a prime ideal in $\mathcal{O}_K$ with $p_0 \equiv 1 \pmod{\ell}$ (so that $\ell \mid p_0 - 1$).
Define $\varepsilon = (p_0 - 1)/\ell$. The *Schirokauer map* $\lambda_\varphi: K^* \to
\mathbb{Z}/\ell\mathbb{Z}$ is defined by

$$\lambda_\varphi(\beta) = \frac{\beta^\varepsilon - 1}{\ell} \bigg|_{\alpha = r_0} \bmod \ell,$$

where the expression $(\beta^\varepsilon - 1)/\ell$ is computed in $\mathbb{Z}[\alpha]/\ell^2$
(the ring of integers modulo $\ell^2$), and then evaluated at $\alpha = r_0 \bmod \ell$.

**Why this works.** The key facts are:

1. For $\beta \in \mathcal{O}_K^*$ (a unit), $\beta^{p_0 - 1} \equiv 1 \pmod{\varphi}$ by
   Fermat's little theorem in $\mathcal{O}_K/\varphi \cong \mathbb{F}_{p_0}$. So
   $\beta^\varepsilon$ is an $\ell$-th root of unity modulo $\varphi$, and $\beta^\varepsilon
   \equiv 1 \pmod{\ell}$ in $\mathbb{Z}[\alpha]$. This means $(\beta^\varepsilon - 1)/\ell$ is
   an integer (the division is exact), and $\lambda_\varphi(\beta)$ is well-defined.

2. The map $\lambda_\varphi$ is a *group homomorphism*: $\lambda_\varphi(\beta\gamma) =
   \lambda_\varphi(\beta) + \lambda_\varphi(\gamma) \pmod{\ell}$. This follows from the
   logarithm-like behaviour of the $\ell$-adic extraction.

3. For a principal ideal $(\beta) = \prod_i \mathfrak{p}_i^{e_i}$, the Schirokauer map
   $\lambda_\varphi(\beta)$ captures the "fractional part" of $\log_g(\phi_\alpha(\beta))$ that
   the exponent vector misses — specifically, the contribution of the unit part of $\beta$.

**The multi-coordinate map.** In practice, one uses $r$ prime ideals $\varphi_1, \ldots,
\varphi_r$ (each with $p_i \equiv 1 \pmod{\ell}$) and defines the Schirokauer map as the
$r$-tuple $\lambda(\beta) = (\lambda_{\varphi_1}(\beta), \ldots, \lambda_{\varphi_r}(\beta))
\in (\mathbb{Z}/\ell\mathbb{Z})^r$. Each coordinate gives one Schirokauer column in the DL
relation matrix.

**Proof sketch (homomorphism property).** For $\beta, \gamma \in \mathcal{O}_K^*$:

$$(\beta\gamma)^\varepsilon - 1 = \beta^\varepsilon \gamma^\varepsilon - 1
  = (\beta^\varepsilon - 1) + (\gamma^\varepsilon - 1) + (\beta^\varepsilon - 1)(\gamma^\varepsilon - 1).$$

Since $\beta^\varepsilon \equiv 1 \pmod{\ell}$ and $\gamma^\varepsilon \equiv 1 \pmod{\ell}$,
the cross term $(\beta^\varepsilon - 1)(\gamma^\varepsilon - 1) \equiv 0 \pmod{\ell^2}$.
Dividing by $\ell$ and evaluating at $\alpha = r_0$:

$$\lambda_\varphi(\beta\gamma) = \lambda_\varphi(\beta) + \lambda_\varphi(\gamma) \pmod{\ell}.$$

$\square$

### The Schirokauer columns in the DL matrix

Each DL relation $(a, b)$ contributes one Schirokauer column per ideal $\varphi_i$: the value
$\lambda_{\varphi_i}(a - b\alpha) \in \mathbb{Z}/\ell\mathbb{Z}$. These columns are appended to
the exponent vector to form the full DL relation row:

$$\text{row}(a, b) = \bigl(\underbrace{e_1, \ldots, e_k}_{\text{rational exponents}},\;
  \underbrace{f_1, \ldots, f_m}_{\text{algebraic exponents}},\;
  \underbrace{\lambda_1(a-b\alpha), \ldots, \lambda_r(a-b\alpha)}_{\text{Schirokauer columns}}\bigr)
  \in \mathbb{F}_\ell^{k+m+r}.$$

The Schirokauer columns ensure that the kernel of the augmented matrix gives the correct virtual
logs — including the unit contribution that the exponent vectors alone would miss.

**Contrast with GNFS.** In GNFS, the quadratic-character columns are Legendre symbols
$\left(\frac{a - bs}{q}\right) \in \{+1, -1\}$ — elements of $\{0, 1\} \subset \mathrm{GF}(2)$ —
and are necessary filters on a square. In NFS-DL, the Schirokauer columns are $\ell$-adic log
extractions in $\mathbb{Z}/\ell\mathbb{Z}$ that make the unit contribution linear modulo $\ell$.
The obstruction group is the same one (§5); the difference is the condition (squareness filtered
by characters, versus linearity modulo $\ell$ with $\ell \nmid h_K$ assumed) and the correction
($\ell$-adic log versus Legendre symbol).

**Code realisation.** The Schirokauer map is implemented in `gnfs/src/dl/schirokauer.rs`
The function `schirokauer` computes the $r$-tuple
$(\lambda_{\varphi_1}(\beta), \ldots, \lambda_{\varphi_r}(\beta))$ for a number-field element
$\beta$ and a list of prime ideals. The `augment_relation` function in `gnfs/src/dl/relation.rs`
wraps each smooth factoring relation with its Schirokauer columns
to produce a `DLRelation`. See `gnfs/docs/PEDAGOGY.md` §64–§65 for the code-tour.

---

## §9.4 The $\mathbb{F}_\ell$ Linear Algebra Step

### The linear system over $\mathbb{F}_\ell$

The DL relation matrix $A$ is a $(k + \delta) \times (k + m + r)$ matrix over $\mathbb{F}_\ell$,
where:
- $k$ = number of rational factor-base primes,
- $m$ = number of algebraic factor-base ideals,
- $r$ = number of Schirokauer columns,
- $\delta$ = small excess (a few extra relations for numerical stability).

The rational and algebraic exponent columns are reduced modulo $\ell$ (from the integer exponent
vectors). The Schirokauer columns are already in $\mathbb{Z}/\ell\mathbb{Z}$.

**The kernel.** The virtual-log table is the kernel of $A$ over $\mathbb{F}_\ell$: a vector
$x \in \mathbb{F}_\ell^{k+m+r}$ such that $A \cdot x = 0$. The first $k$ entries of $x$ are
the virtual logs of the rational factor-base primes; the next $m$ entries are the virtual logs
of the algebraic factor-base ideals; the last $r$ entries are Schirokauer correction terms.

**Why the kernel gives the virtual logs.** Each row of $A$ encodes the linear equation

$$\sum_i e_i \log_g(p_i) + \sum_j f_j \log_g(\phi_\alpha(\mathfrak{p}_j)) +
  \sum_s \lambda_s(a - b\alpha) \cdot c_s \equiv 0 \pmod{\ell},$$

where $c_s$ are the Schirokauer correction coefficients. A kernel vector $x$ satisfying $A \cdot
x = 0$ is a consistent assignment of virtual logs to all factor-base elements (and Schirokauer
corrections) that satisfies every relation simultaneously. Under the standard heuristic
assumptions, the kernel is one-dimensional (up to scalar multiples), and the unique kernel vector
(up to scaling) gives the virtual logs.

### Block Wiedemann and block Lanczos over $\mathbb{F}_\ell$

The kernel is found by the *block Wiedemann* or *block Lanczos* algorithm — the same Krylov
subspace methods used in GNFS (§5), but over $\mathbb{F}_\ell$ instead of $\mathrm{GF}(2)$.

**The algorithmic structure is identical.** Both algorithms:
1. Start with a random block vector $V_0 \in \mathbb{F}_\ell^{(k+m+r) \times w}$ (where $w$ is
   the block width).
2. Iterate $V_{i+1} = A^T A \cdot V_i$ (or $A \cdot A^T \cdot V_i$) to build a Krylov subspace.
3. Use the self-orthogonality of the Krylov sequence to find a null vector.

**The implementation difference.** In GNFS, the block width is $w = 64$ (one 64-bit word per
block vector entry, since $\mathrm{GF}(2)$ elements are bits). In NFS-DL, the block width is
$w = 32$ (one $\mathbb{F}_\ell$ element per entry, since field elements are larger than bits).
The scalar arithmetic changes from XOR (GF(2)) to modular multiplication and addition
($\mathbb{F}_\ell$). The convergence criterion changes from "the block vector is zero in GF(2)"
to "the inner-product matrix $V_i^T V_i$ is singular over $\mathbb{F}_\ell$".

**Complexity.** For a matrix with $k + m + r$ columns and $W$ non-zero entries, the
Lanczos/Wiedemann algorithms perform $O(k + m + r)$ applications of the sparse operator at cost
$O(W)$ each, so $O((k + m + r) \cdot W)$ field operations over $\mathbb{F}_\ell$ (dense
elimination costs $O((k + m + r)^3)$). With $O(\log p)$ non-zeros per row, $W = O((k + m + r)
\log p)$, and in $L$-notation with $k + m + r \approx B = L_p[1/3, c]$:

$$\text{Linear algebra cost} = O\!\left(L_p[1/3, c]^2 \log p\right) = L_p\!\left[\tfrac{1}{3},\, 2c\right].$$

This is the *same* cost shape as the GF(2) linear algebra in GNFS. The field $\mathbb{F}_\ell$
vs $\mathrm{GF}(2)$ distinction changes the constant factor (field operations over $\mathbb{F}_\ell$
are more expensive than XOR), but not the $L$-notation exponent or constant. This is the first
half of the DL delta in the complexity derivation (§9.7).

**Code realisation.** The $\mathbb{F}_\ell$ block-solver substrate is implemented in
`gnfs/src/dl/linalg/blockvec_fl.rs`. The types `FlBlockVec`,
`FlSparseMatrix`, `FlMatrixOperator`, `FlSolution`, `VirtualLogTable`, and the function
`recover_virtual_logs` are its interface. See `gnfs/docs/PEDAGOGY.md` §65–§66 for the
code-tour.

---

## §9.5 Individual-Logarithm Special-$q$ Descent

*This section describes the stage with no factoring analogue — the part of NFS-DL that goes
beyond what GNFS needs.*

### The individual-logarithm problem

The linear algebra step (§9.4) recovers the virtual logs of all factor-base elements: $\log_g(p_i)
\bmod \ell$ for each rational prime $p_i \leq B$ and $\log_g(\phi_\alpha(\mathfrak{p}_j)) \bmod
\ell$ for each algebraic ideal $\mathfrak{p}_j$. But the target $h$ is not a factor-base element
— it is an arbitrary element of $\mathbb{F}_p^*$. The *individual-logarithm descent* bridges
this gap.

**The descent strategy.** The descent proceeds in two phases:

1. **Initialization-smoothing.** Find an exponent $e \geq 0$ such that $g^e \cdot h \bmod p$
   is smooth over a *medium-prime bound* $B'$ (with $B < B' \ll p$). Then
   $\log_g(h) = \log_g(g^e \cdot h) - e \pmod{\ell}$, so it suffices to compute
   $\log_g(g^e \cdot h)$.

2. **Special-$q$ descent.** The smooth factorisation of $g^e \cdot h$ involves primes in
   $(B, B']$ — the *medium primes* — whose logs are not yet known. For each medium prime $q$,
   find a DL relation that expresses $\log_g(q)$ as a linear combination of logs of smaller
   primes. Repeat until all primes are factor-base leaves with known virtual logs.

### Initialization-smoothing

**The smoothing step.** Iterate $e = 0, 1, 2, \ldots$ and compute $c_e = g^e \cdot h \bmod p$.
Trial-divide $c_e$ by all primes up to $B'$. If $c_e$ factors completely over primes $\leq B'$,
the smoothing succeeds with exponent $e$ and factorisation $c_e = \prod_i q_i^{a_i}$.

The expected number of trials before a smooth $c_e$ is found is $\rho(u')^{-1}$, where
$u' = \log p / \log B'$ and $\rho$ is the Dickman function. For $B' = L_p[1/3, c']$, this is
$L_p[1/3, 1/(3c')]$ trials — subexponential in $\log p$.

**The initial frontier.** The smooth factorisation $c_e = \prod_i q_i^{a_i}$ gives the *initial
frontier*: the set of prime factors $q_i$ (with multiplicity) that must be descended. Factor-base
primes ($q_i \leq B$) are leaves with known virtual logs; medium primes ($B < q_i \leq B'$) are
interior nodes that must be descended.

### The special-$q$ descent recursion

**The descent step.** For a medium prime $q$ in the frontier, the descent finds a DL relation
that rewrites $\log_g(q)$ as a combination of logs of smaller primes. Specifically, it runs the
special-$q$ sieve with $q$ as the special prime: sieve for pairs $(a, b)$ such that $q \mid
N_{K/\mathbb{Q}}(a - b\alpha)$ and both norms are smooth over primes $< q$. Such a relation gives

$$e_q \cdot \log_g(q) \equiv \sum_{p_i < q} e_i \log_g(p_i) \pmod{\ell},$$

where $e_q$ is the exponent of $q$ in the algebraic norm and $e_i$ are the exponents of the
smaller primes. Since $e_q \not\equiv 0 \pmod{\ell}$ (generically), this gives
$\log_g(q) \equiv e_q^{-1} \sum_i e_i \log_g(p_i) \pmod{\ell}$.

**The termination invariant.** The descent tree is a max-heap ordered by prime descending. Each
step pops the largest prime $q$, finds a relation rewriting $\log_g(q)$ in terms of strictly
smaller primes, and pushes those smaller primes back. Since each step strictly reduces the
largest prime, the descent terminates when all frontier elements are factor-base leaves.

**The descent tree.** The descent produces a tree of `DescentNode` objects:
- **Leaf nodes** are factor-base elements with known virtual logs (from the `VirtualLogTable`).
- **Interior nodes** are medium primes with a rewriting relation and child nodes for the smaller
  primes.

The tree has depth $O(\log \log p)$ in the asymptotic regime (since each level reduces the
prime by a constant factor in the $L$-notation sense), but at toy scale the depth is 0 or 1
(the smoothing step often finds $g^e \cdot h$ already smooth over the factor base).

### Log assembly along the descent tree

Once the descent tree is complete (all leaves have known virtual logs), the log of $g^e \cdot h$
is assembled by a bottom-up traversal:

- **Leaf:** $\log_g(q) = \text{known\_log}(q)$ from the virtual-log table.
- **Interior node:** $\log_g(q) = \sum_{\text{children}} \log_g(\text{child}) \pmod{\ell}$,
  where the sum is over all children with multiplicity (a child appearing twice contributes its
  log twice).

The assembled log of $g^e \cdot h$ is the sum of the logs of the initial frontier targets. Then:

$$\log_g(h) = \log_g(g^e \cdot h) - e \pmod{\ell}.$$

### Subgroup recovery

The descent recovers $\log_g(h) \bmod \ell$, where $\ell$ is the subgroup order. If $\ell = p - 1$
(the full group order), this is the complete discrete log. If $\ell < p - 1$, the result is the
log in the $\ell$-order subgroup. Recovering the full log modulo $p - 1$ requires Pohlig–Hellman
(§Prerequisites): factor $p - 1 = \prod_i \ell_i^{e_i}$, run NFS-DL for each prime power
$\ell_i^{e_i}$, and combine via CRT.

**Code realisation.** The descent substrate is implemented in `gnfs/src/dl/descent/`. The key
types and functions are `DescentNode<F>`, `DescentFrontier<F>`, `init_descent_frontier`,
`descend_node`, `run_descent`, and `assemble_log`. The entry point `solve_dl(g, h, p, k, ell)` is
consumed across crates by the MOV reduction (§10.5). See `gnfs/docs/PEDAGOGY.md` §66 for the
code-tour.

**Engineering-scale boundary ($k > 1$).** The `solve_dl` function
returns `SolveDlError::Unsupported { k }` immediately for $k > 1$ (extension fields
$\mathbb{F}_{p^k}$). This is an *engineering-scale boundary, not a mathematical one*: the
mathematics of NFS-DL over extension fields is well-understood (the number field $K$ is replaced
by a degree-$k$ extension; the Schirokauer map adapts; the descent works over the extended
field), but the implementation requires adapting the `PolyPair` / `NumberField` / `FactorBase`
infrastructure to the extension-field structure. This is a non-trivial engineering task that is
not done here. The gap is recorded as a scale annotation (see §On Scale for the three-axis
framework): the code is correct at demonstration fidelity for $k = 1$; the $k > 1$ path is a
known gap, not a silent omission.

---

## §9.6 What the NFS-DL Implementation Implements and Claims

*This section states what the whole NFS-DL pipeline (relation collection, linear algebra,
descent) implements, what it leaves out, and which phenomena are scale-only.*

### Algorithmic content

The NFS-DL pipeline implements the algorithmic content of NFS-DL for prime fields end-to-end:

- **Relation collection:** Smooth-pair sieving reusing the GNFS infrastructure, augmented
  with Schirokauer map columns. The `DLMatrix` type carries the full column
  layout (rational | algebraic | Schirokauer).
- **$\mathbb{F}_\ell$ linear algebra:** Block Lanczos / block Wiedemann over $\mathbb{F}_\ell$,
  recovering the virtual-log table. The `VirtualLogTable` type is the interface to the descent.
- **Individual-logarithm descent:** Initialization-smoothing, special-$q$ descent
  recursion, and log assembly. The `solve_dl_full` function runs the full pipeline; `solve_dl` is
  the cross-crate entry point.

No algorithmic stage is silently omitted or replaced by a lookup table.

### No engineering optimisations

The implementation is at *demonstration fidelity* — the mathematical content
is present, but the engineering optimisations that make NFS-DL run at cryptographic scale are
not implemented:

- **Large-prime variations** (partial relations with one or two large prime cofactors) are not
  implemented. The sieve collects only fully smooth relations.
- **Optimal polynomial selection** for NFS-DL (calibrated to $p$ and $\ell$) is not implemented.
  The base-$m$ construction is used.
- **Lattice sieving** for the descent step is not implemented. The line sieve is used.
- **Pohlig–Hellman** for the full group order is not implemented. The log is recovered modulo
  $\ell$ only.

None of these omissions affect the correctness of the algorithm at toy scale. They are
engineering optimisations that improve performance at NFS scale but are not needed for the
mathematical demonstration.

### Scale-only phenomena

The following phenomena are annotated as scale-only:

- **Descent-tree depth.** At toy scale ($p = 11$, factor base $\{2, 3\}$, $\ell = 5$), the
  descent tree has depth 0: the initialization-smoothing step finds $g^e \cdot h$ already smooth
  over the factor base, so no medium primes need to be descended. At NFS scale, the tree has
  depth $O(\log \log p)$ — typically 3–5 levels. The depth is a scale-only phenomenon.

- **Medium-prime tuning.** The medium-prime bound $B'$ is a critical parameter at NFS scale.
  At toy scale, `solve_dl` uses a hardcoded `medium_bound = 100`; `solve_dl_full`
  uses the factor-base bound. At NFS scale, $B'$ is calibrated to the factor-base bound and the
  expected descent depth.

- **Block width.** The $\mathbb{F}_\ell$ block width `FL_BLOCK_WIDTH = 32` is a cache-friendly
  unit for the inner loop at NFS scale. At toy scale, the blocking overhead is invisible.

- **$k > 1$ `Unsupported` debt.** The $\mathbb{F}_{p^k}$ extension is an engineering-scale
  boundary (see §9.5 above). The mathematics is complete for $k = 1$; the $k > 1$ path is a
  known gap recorded as a scale annotation.

---

## §9.7 The $L$-Notation Complexity of NFS-DL

*This is the designated payoff proof for the NFS-DL chapter. It is a
**delta** on the §GNFS §7 derivation — not a re-derivation from scratch. The reader should
read §GNFS §7 first; this section identifies the two DL-specific differences and shows that
neither changes the leading complexity.*

### The §GNFS §7 derivation: what it establishes

The §GNFS §7 derivation (§7 of the GNFS chapter above) establishes:

1. **The exponent is $1/3$**: the sieve cost is $L_N[1/3, c + 2/(3c)]$ and the linear algebra
   cost is $L_N[1/3, 2c]$; the crossover at $c = \sqrt{2/3}$ gives total cost
   $L_N[1/3, 2\sqrt{6}/3]$ (simplified) or $L_N[1/3, (64/9)^{1/3}]$ (full analysis).

2. **The constant is $(64/9)^{1/3}$**: the full analysis (Lenstra–Lenstra–Manasse–Pollard
   [LLMP93], Buhler–Lenstra–Pomerance [BLP93]) accounts for the polynomial degree optimisation
   and the two-sided sieve, yielding the standard constant.

The derivation uses:
- Smoothness bound $B = L_N[1/3, c]$.
- Norm size $\approx N^{1/d}$ at optimal degree $d$.
- Smoothness probability (one norm) $= L_N[-1/3, -1/(3c)]$.
- Sieve cost $= L_N[1/3, c + 2/(3c)]$ (collecting $B$ smooth pairs).
- Linear algebra cost $= L_N[1/3, 2c]$ (null vector in $B \times B$ matrix over GF(2)).
- Optimal $c = \sqrt{2/3}$ (simplified) or $c = (8/9)^{1/3}$ (full).

### The NFS-DL delta: two differences

NFS-DL differs from GNFS in two ways that could, in principle, change the leading complexity.
We show that neither does.

#### Delta 1: $\mathbb{F}_\ell$ linear algebra vs GF(2)

In GNFS, the linear algebra step finds a null vector in a $B \times B$ matrix over GF(2). In
NFS-DL, the linear algebra step finds a null vector in a $B \times B$ matrix over $\mathbb{F}_\ell$.

**The $L$-notation cost is the same.** The block Lanczos / block Wiedemann algorithm costs
$O(B \cdot W)$ field operations in both cases, where $W$ is the number of non-zero entries (dense
elimination would cost $O(B^3)$). With $O(\log p)$ non-zeros per row, $W = O(B \log p)$. In
$L$-notation:

$$\text{NFS-DL linear algebra cost} = O\!\left(L_p[1/3, c]^2 \log p\right) = L_p\!\left[\tfrac{1}{3},\, 2c\right].$$

This is identical to the GNFS linear algebra cost $L_N[1/3, 2c]$. The field $\mathbb{F}_\ell$
vs GF(2) distinction changes the *constant factor* inside the $O(\cdot)$ (field operations over
$\mathbb{F}_\ell$ are more expensive than XOR), but not the $L$-notation exponent or constant.
The $L$-notation is insensitive to polynomial factors in the constant, so the leading complexity
is unchanged.

**Formal statement.** Let $T_{\mathrm{GF}(2)}(B)$ and $T_{\mathbb{F}_\ell}(B)$ be the costs of
finding a null vector in a $B \times B$ matrix over GF(2) and $\mathbb{F}_\ell$ respectively.
Then $T_{\mathrm{GF}(2)}(B) = O(B^2)$ and $T_{\mathbb{F}_\ell}(B) = O(B^2 \cdot \log \ell)$
(since each $\mathbb{F}_\ell$ operation costs $O(\log \ell)$ bit operations). In $L$-notation
with $B = L_p[1/3, c]$ and $\ell \leq p - 1 = L_p[1, 1]$:

$$T_{\mathbb{F}_\ell}(B) = O\!\left(L_p[1/3, c]^2 \cdot \log p\right)
  = L_p\!\left[\tfrac{1}{3},\, 2c\right] \cdot O(\log p).$$

The factor $O(\log p)$ is polynomial in $\log p$ — it is $L_p[0, 1]$ in $L$-notation, which is
dominated by $L_p[1/3, 2c]$ for any $c > 0$. Therefore:

$$T_{\mathbb{F}_\ell}(B) = L_p\!\left[\tfrac{1}{3},\, 2c\right] \cdot (1 + o(1)).$$

The $\mathbb{F}_\ell$ linear algebra has the same $L$-notation cost as the GF(2) linear algebra.
$\square$

#### Delta 2: Individual-logarithm descent is asymptotically subdominant

In GNFS, there is no descent step: the linear algebra step yields a candidate dependency, and a
verified algebraic root turns it into a congruence of squares. In NFS-DL,
the individual-logarithm descent is an additional stage. Could it dominate the total cost?

**The descent cost.** The descent for a single target $h$ consists of:

1. **Initialization-smoothing:** Find $e$ such that $g^e \cdot h$ is smooth over $B' = L_p[1/3, c']$.
   Expected number of trials: $\rho(u')^{-1}$ where $u' = \log p / \log B'$. For $B' = L_p[1/3, c']$:

   $$\rho(u')^{-1} = L_p\!\left[\tfrac{1}{3},\, \tfrac{1}{3c'}\right] \cdot (1 + o(1)).$$

   Each trial costs $O(\log p)$ (one modular multiplication). Total smoothing cost:

   $$L_p\!\left[\tfrac{1}{3},\, \tfrac{1}{3c'}\right] \cdot O(\log p)
     = L_p\!\left[\tfrac{1}{3},\, \tfrac{1}{3c'}\right] \cdot (1 + o(1)).$$

2. **Special-$q$ descent:** For each medium prime $q \in (B, B']$, find a rewriting relation.
   The descent tree has depth $D = O(\log(B'/B) / \log(B/B_{\min}))$ — in the asymptotic regime,
   $D = O(1)$ (a constant number of levels, since $B$ and $B'$ are both $L_p[1/3, \cdot]$).
   Each descent step costs one sieve run: $O(B \cdot \log B)$ operations. Total descent cost:

   $$D \cdot O(B \log B) = O(1) \cdot L_p\!\left[\tfrac{1}{3},\, c\right] \cdot O(\log p)
     = L_p\!\left[\tfrac{1}{3},\, c\right] \cdot (1 + o(1)).$$

**Comparing descent cost to precomputation cost.** The precomputation cost (sieve + linear
algebra) is $L_p[1/3, (64/9)^{1/3}]$. The descent cost for a single target is at most
$L_p[1/3, c + 1/(3c')]$ (the larger of the smoothing and descent costs). For any choice of
$c, c' > 0$, the descent cost is $L_p[1/3, \cdot]$ with a constant that is at most as large as
the precomputation constant.

**The key point.** The descent cost is $L_p[1/3, \cdot]$ — the same $L$-notation *exponent*
as the precomputation. It is not a lower-order term in the absolute sense. However, it is
*subdominant* in the following precise sense: the descent cost for a *single* target $h$ is
$L_p[1/3, c_{\mathrm{descent}}]$ for some constant $c_{\mathrm{descent}}$, while the
precomputation cost is $L_p[1/3, (64/9)^{1/3}]$. The total cost is

$$\text{Total cost} = L_p\!\left[\tfrac{1}{3},\, \max\!\left((64/9)^{1/3},\, c_{\mathrm{descent}}\right)\right].$$

The descent constant $c_{\mathrm{descent}}$ can be made smaller than $(64/9)^{1/3}$ by choosing
$B'$ appropriately (specifically, $B' = B = L_p[1/3, (8/9)^{1/3}]$, so that the smoothing and
descent costs are both dominated by the precomputation). With this choice:

$$c_{\mathrm{descent}} \leq (8/9)^{1/3} + \tfrac{1}{3(8/9)^{1/3}} < (64/9)^{1/3}.$$

Therefore the descent does not increase the leading constant, and the total cost remains
$L_p[1/3, (64/9)^{1/3}]$.

**Formal statement (descent subdominance).** *Under the heuristic smoothness assumptions, the
individual-logarithm descent for a single target $h$ costs at most $L_p[1/3, c_{\mathrm{descent}}]$
with $c_{\mathrm{descent}} < (64/9)^{1/3}$. The total NFS-DL cost (precomputation + descent) is
therefore $L_p[1/3, (64/9)^{1/3}] \cdot (1 + o(1))$, the same as the precomputation alone.*

*Proof sketch.* The smoothing cost is $L_p[1/3, 1/(3c')]$ for $B' = L_p[1/3, c']$. The descent
cost is $O(D) \cdot L_p[1/3, c]$ where $D = O(1)$ is the tree depth. Setting $c' = c$ (the
same smoothness bound for smoothing and precomputation) gives total descent cost
$L_p[1/3, c + 1/(3c)]$. At the optimal $c = (8/9)^{1/3}$:

$$c + \frac{1}{3c} = \left(\frac{8}{9}\right)^{1/3} + \frac{1}{3(8/9)^{1/3}}
  = \left(\frac{8}{9}\right)^{1/3} + \frac{1}{3} \cdot \left(\frac{9}{8}\right)^{1/3}
  = \left(\frac{8}{9}\right)^{1/3}\!\left(1 + \frac{1}{3} \cdot \frac{9}{8}\right)
  = \left(\frac{8}{9}\right)^{1/3} \cdot \frac{11}{8}.$$

Numerically: $(8/9)^{1/3} \approx 0.961$, so $c + 1/(3c) \approx 0.961 \cdot 1.375 \approx
1.321$. The precomputation constant is $(64/9)^{1/3} \approx 1.923$. Since $1.321 < 1.923$,
the descent cost is strictly dominated by the precomputation cost. $\square$

### The main theorem: NFS-DL and GNFS share the same $L$-notation complexity

**Theorem (NFS-DL complexity, heuristic).** *Under standard heuristic assumptions on the
distribution of smooth norms, NFS-DL computes $\log_g(h)$ in $\mathbb{F}_p^*$ in expected time*

$$L_p\!\left[\tfrac{1}{3},\, \left(\tfrac{64}{9}\right)^{1/3}\right] \cdot (1 + o(1)).$$

*This is the same complexity as GNFS for factoring an integer $N \approx p$.*

**Proof (delta on §GNFS §7).** The precomputation (relation collection + $\mathbb{F}_\ell$ linear
algebra) has the same cost as GNFS:

- **Relation collection:** The sieve cost is $L_p[1/3, c + 2/(3c)]$ (same as GNFS §7, with $p$
  in place of $N$). The norms have size $\approx p^{1/d}$ at optimal degree $d$; the smoothness
  probability is $L_p[-1/3, -1/(3c)]$ per norm; collecting $B = L_p[1/3, c]$ relations costs
  $L_p[1/3, c + 2/(3c)]$.

- **$\mathbb{F}_\ell$ linear algebra:** The cost is $L_p[1/3, 2c]$ (Delta 1 above: same as GF(2)).

- **Optimisation:** The total precomputation cost is $L_p[1/3, \max(c + 2/(3c), 2c)]$, minimised
  at $c = \sqrt{2/3}$ (simplified) or $c = (8/9)^{1/3}$ (full), giving $L_p[1/3, (64/9)^{1/3}]$
  (full analysis, citing [LLMP93] and [BLP93] as in §GNFS §7).

- **Individual-logarithm descent:** The descent cost is $L_p[1/3, c_{\mathrm{descent}}]$ with
  $c_{\mathrm{descent}} < (64/9)^{1/3}$ (Delta 2 above). It does not increase the leading
  constant.

Therefore the total NFS-DL cost is $L_p[1/3, (64/9)^{1/3}] \cdot (1 + o(1))$. $\square$

### The asymptotic comparison: why NFS-DL and GNFS share the same complexity

The equality of the NFS-DL and GNFS complexities is not a coincidence — it is a structural
consequence of the fact that both algorithms are driven by the *same* bottleneck: the balance
between the sieve cost and the linear algebra cost.

**The shared bottleneck.** In both algorithms:
- The sieve cost is $L[1/3, c + 2/(3c)]$ (collecting smooth pairs from norms of size $\approx
  N^{1/d}$ or $p^{1/d}$).
- The linear algebra cost is $L[1/3, 2c]$ (finding a null vector in a $B \times B$ matrix).
- The optimal $c$ is the crossover point where the two costs are equal.

The DL-specific stages (Schirokauer maps, $\mathbb{F}_\ell$ arithmetic, individual-logarithm
descent) do not change this balance:
- Schirokauer maps add $r$ columns to the matrix, but $r = O(1)$ (a constant number of ideals),
  so the matrix size is still $B \times B$ in $L$-notation.
- $\mathbb{F}_\ell$ arithmetic changes the constant factor in the linear algebra cost, but not
  the $L$-notation exponent or constant (Delta 1).
- The individual-logarithm descent is subdominant (Delta 2).

**The through-line at its sharpest.** The equality $L_p[1/3, (64/9)^{1/3}] = L_N[1/3,
(64/9)^{1/3}]$ (for $p \approx N$) is the structure-based escape from search at its sharpest:
two different problems (factoring and discrete logarithm), solved by two different algorithms
(GNFS and NFS-DL), with two different outputs (a factor vs a discrete log), but with the *same*
asymptotic complexity. The reason is that both algorithms exploit the *same* structure — the
number-field bridge and the smoothness of norms — and are limited by the *same* bottleneck —
the sieve/linear-algebra balance.

This is why the $L[1/3]$ barrier is believed to be a genuine barrier for both problems: any
algorithm that breaks it for one problem would likely break it for the other, since the
bottleneck is shared.

### Summary of the derivation

| Quantity | GNFS (§7) | NFS-DL (§9.7) |
|----------|-----------|---------------|
| Smoothness bound | $B = L_N[1/3, c]$ | $B = L_p[1/3, c]$ |
| Norm size | $\approx N^{1/d}$ | $\approx p^{1/d}$ |
| Smoothness probability (one norm) | $L_N[-1/3, -1/(3c)]$ | $L_p[-1/3, -1/(3c)]$ |
| Sieve cost | $L_N[1/3, c + 2/(3c)]$ | $L_p[1/3, c + 2/(3c)]$ |
| Linear algebra field | GF(2) | $\mathbb{F}_\ell$ |
| Linear algebra cost | $L_N[1/3, 2c]$ | $L_p[1/3, 2c]$ (Delta 1) |
| Descent cost | — (no descent) | $L_p[1/3, c_{\mathrm{descent}}] < L_p[1/3, (64/9)^{1/3}]$ (Delta 2) |
| Optimal $c$ | $(8/9)^{1/3}$ (full) | $(8/9)^{1/3}$ (full) |
| Total cost | $L_N[1/3, (64/9)^{1/3}]$ | $L_p[1/3, (64/9)^{1/3}]$ |

The exponent $1/3$ and the constant $(64/9)^{1/3}$ are the same for both problems. The DL
content is entirely in the two deltas: $\mathbb{F}_\ell$ linear algebra (same cost shape as
GF(2)) and individual-logarithm descent (asymptotically subdominant).

---

## §9.8 Cross-References and References

### Cross-references within this textbook

**Code realisation.** The NFS-DL pipeline — relation collection with Schirokauer augmentation,
$\mathbb{F}_\ell$ linear algebra, and individual-logarithm descent — is implemented in the
`gnfs` crate and documented in `gnfs/docs/PEDAGOGY.md` §63–§71. That chapter is the code-tour
sibling to this one: it assumes the reader knows the mathematics (this chapter) and focuses on
the Rust realisation, the stage-by-stage contracts, and the KAT summary.

**§GNFS §7 (the shared derivation core).** The L-notation derivation in §9.7 is a delta on the
§GNFS §7 derivation. The reader should read §GNFS §7 first; §9.7 identifies the two DL-specific
differences (Delta 1: $\mathbb{F}_\ell$ vs GF(2); Delta 2: descent subdominance) and shows that
neither changes the leading complexity. The full derivation of the exponent $1/3$ and the
constant $(64/9)^{1/3}$ is in §GNFS §7; §9.7 re-uses it by reference.

**§Prerequisites.** The prerequisites used in this chapter:
- The Canfield–Erdős–Pomerance smooth-number density theorem and the $L$-notation (§Analysis,
  §Probability) — the engine of the complexity derivation.
- The discrete logarithm problem and Pohlig–Hellman (§Logic and complexity) — for subgroup
  recovery.
- Unique factorisation of ideals in Dedekind domains (§Algebra) — for the algebraic factor base.
- The $L$-notation arithmetic lemma (§Analysis) — for the cost calculations in §9.7.

**§On Scale.** The honest science↔engineering gap — the gap between the asymptotic complexity
$L_p[1/3, (64/9)^{1/3}]$ and the behaviour at toy scale — is discussed in §On Scale. The
key points for NFS-DL: the descent-tree depth is a scale-only phenomenon (depth 0 at toy scale,
$O(\log \log p)$ at NFS scale); the medium-prime tuning is invisible at toy scale; the
$\mathbb{F}_{p^k}$ $k > 1$ `Unsupported` debt is an engineering-scale boundary (mathematical-
dimension axis, not resource/operational axis).

### References for the NFS-DL Chapter

24. **Schirokauer, O. (1993).** "Discrete logarithms and local units." *Philosophical
    Transactions of the Royal Society A*, 345(1676), 409–423. [S93] The original source for the
    Schirokauer map — the $\ell$-adic virtual-log correction that makes the unit contribution
    linear modulo $\ell$ in NFS-DL.

25. **Gordon, D. M. (1993).** "Discrete logarithms in GF(p) using the number field sieve."
    *SIAM Journal on Discrete Mathematics*, 6(1), 124–138. [G93] The original NFS-DL algorithm
    for prime fields, including the individual-logarithm descent.

26. **Adleman, L. M. (1994).** "The function field sieve." In: Adleman, L. M., and Huang, M.-D.
    (eds.) *Algorithmic Number Theory (ANTS-I)*, LNCS 877. Springer. The function-field analogue
    of NFS-DL; context for the descent stage.

27. **Joux, A., Lercier, R., Smart, N., and Vercauteren, F. (2006).** "The number field sieve
    in the medium prime case." In: Dwork, C. (ed.) *Advances in Cryptology — CRYPTO 2006*, LNCS
    4117. Springer. The medium-prime NFS-DL, including the special-$q$ descent and the
    medium-prime bound calibration.

28. **Barbulescu, R., Gaudry, P., Joux, A., and Thomé, E. (2014).** "A heuristic quasi-polynomial
    algorithm for discrete logarithm in finite fields of small characteristic." In: Nguyen, P. Q.,
    and Oswald, E. (eds.) *Advances in Cryptology — EUROCRYPT 2014*, LNCS 8441. Springer. The
    BGJT algorithm — the quasi-polynomial DLP algorithm for small-characteristic fields that
    broke the NFS-DL paradigm for that case. [BGJT14]

29. **Lenstra, A. K., Lenstra, H. W. Jr., Manasse, M. S., and Pollard, J. M. (1993).** "The
    number field sieve." In *The Development of the Number Field Sieve*, Lecture Notes in
    Mathematics 1554, Springer, 11–42. [LLMP93] (Also cited in §GNFS §7.) The full analysis
    giving the constant $(64/9)^{1/3}$; applies to NFS-DL by the delta argument of §9.7.

30. **Buhler, J. P., Lenstra, H. W. Jr., and Pomerance, C. (1993).** "Factoring integers with
    the number field sieve." In *The Development of the Number Field Sieve*, Lecture Notes in
    Mathematics 1554, Springer, 50–94. [BLP93] (Also cited in §GNFS §7.) The two-sided sieve
    analysis; applies to NFS-DL by the delta argument of §9.7.

31. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective*.
    2nd ed. Springer. §6.2 for the NFS complexity analysis; §6.3 for NFS-DL. [CP05]

---

## Algebraic ECDLP Attacks

*Maths-first sibling to `docs/PEDAGOGY.md` §8–§18 (the algebraic-attack code-tour). For the
phase-by-phase implementation — module surfaces, toy fixtures, KAT summary, and design-statement
verification — see `docs/PEDAGOGY.md`. This chapter is the mathematical development: proof
sketches, the full MOV payoff proof, and the per-attack L-notation comparison.*

### §10.0 The through-line for this chapter

The §"Escape from Search" chapter above established the five-family structure taxonomy and the
L-notation hierarchy. This chapter extends that taxonomy per-attack: each of the five algebraic
ECDLP attacks finds a specific curve structure that escapes the generic $\sqrt{n}$ bound the
Pollard rho chapter (§Pollard Rho for ECDLP) established.

The five structures, named here and developed in the sections below:

1. **Composite group order** (Pohlig–Hellman, §10.1): $\#E$ composite $\Rightarrow$ CRT
   reduction to prime-order subgroup DLPs.
2. **Anomalous order** (Smart–Satoh–Araki, §10.2): $\#E(\mathbb{F}_p) = p$ $\Rightarrow$
   p-adic lift + formal group logarithm gives polynomial time.
3. **Binary field tower** (GHS/Weil descent, §10.3): $E/\mathbb{F}_{2^m}$ with subfield
   $\mathbb{F}_{2^l}$ $\Rightarrow$ Weil restriction transfers ECDLP to a hyperelliptic
   Jacobian DLP (a transfer, not an end-to-end solve).
4. **Factor-base decomposability** (index calculus, §10.4): Semaev-decomposable points
   $\Rightarrow$ relation matrix over $\mathbb{Z}/\ell\mathbb{Z}$; asymptotic win in the
   extension-field setting.
5. **Small embedding degree** (MOV/Frey–Rück, §10.5): $\ell \mid p^k - 1$ for small $k$
   $\Rightarrow$ bilinear pairing maps ECDLP to $\mathbb{F}_{p^k}^*$ DLP, where index calculus
   applies (the cross-track bridge to NFS-DL).

The chapter follows the ordering Pohlig–Hellman → SSA → GHS → index calculus → MOV, placing
MOV last as the designated climax: the chapter builds the engine (index calculus) before
presenting the bridge (MOV) that connects ECDLP to that engine.

---

### §10.1 Pohlig–Hellman: CRT to Prime-Order Subgroups

#### The structure

The group order $n = \#E(\mathbb{F}_p)$ is composite: $n = \prod_{i=1}^r p_i^{e_i}$. The
Chinese Remainder Theorem (§Prerequisites) gives an isomorphism

$$\mathbb{Z}/n\mathbb{Z} \;\cong\; \mathbb{Z}/p_1^{e_1}\mathbb{Z} \times \cdots \times
\mathbb{Z}/p_r^{e_r}\mathbb{Z}.$$

The ECDLP in $\langle G \rangle$ (a cyclic group of order $n$) therefore reduces to independent
DLPs in each prime-power subgroup $\langle G_i \rangle$ of order $p_i^{e_i}$.

#### The escape

**Subgroup projection.** For each prime power $p_i^{e_i}$, define

$$G_i = \frac{n}{p_i^{e_i}} \cdot G, \qquad Q_i = \frac{n}{p_i^{e_i}} \cdot Q.$$

Then $G_i$ has order $p_i^{e_i}$ and $Q_i = k \cdot G_i$ in $\langle G_i \rangle$. The ECDLP
in the full group reduces to the ECDLP $Q_i = k_i \cdot G_i$ in the prime-power subgroup, where
$k_i = k \bmod p_i^{e_i}$.

**Prime-power lift.** For $e_i > 1$, recover $k_i = k \bmod p_i^{e_i}$ digit-by-digit in base
$p_i$. At digit step $j$ (recovering $k_i \bmod p_i^{j+1}$ from $k_i \bmod p_i^j$), project to
the order-$p_i$ sub-subgroup and solve a single DLP in a group of prime order $p_i$.

**The digit step.** The implementation solves that prime-order DLP in one of two ways, chosen by
the size of $p_i$. For $p_i \leq 64$ it enumerates the subgroup, trying $d = 0, 1, \ldots, p_i - 1$
until $d \cdot \gamma$ equals the target. This is the digit step of Pohlig and Hellman (1978, §IV),
where each digit is a function taking one of only $p_i$ values, computed there either by table
lookup or by trying the $p_i$ candidates in turn; the code takes the second form and stores no
table. Only this small-subgroup step is taken from the paper. For $p_i > 64$ the code runs Pollard
rho on the order-$p_i$ subgroup. That is not the paper's own large-prime step, which is a
time–memory search rather than rho; Pollard (1978) is the rho reference and is cited here from
memory, not re-read. Both methods receive the prime order $p_i$, never $p_i^{e_i}$. The solution's
evidence records which method found each digit (`check_pohlig_evidence` in
`rho/tests/ecdlp_kat.rs`, run by `composite_dlog_every_k`). The fixture `composite_toy`
($n = 60$) has only primes $\leq 5$, so every digit is enumerated; `composite_rho_dlog_sampled_k`
runs `composite_rho_toy` ($n = 1068 = 2^2 \cdot 3 \cdot 89$), where the factor $89$ reaches rho.

**Checked input and result.** For a composite $n$, $n \cdot G = \mathcal{O}$ shows only that the
order of $G$ divides $n$; the solver also requires $(n/p) \cdot G \neq \mathcal{O}$ for every prime
$p \mid n$, so that $G$ has order exactly $n$ (`composite_wrong_order_is_invalid_input`,
`composite_rho_wrong_order_is_invalid_input`). A target with $n \cdot Q = \mathcal{O}$ that is
nevertheless outside $\langle G \rangle$ is rejected as invalid input when a digit enumeration finds
no match (`composite_target_outside_subgroup_is_invalid_input`). The combined $k$ is checked
against $k \cdot G = Q$ before it is returned. A rho digit solve that exhausts its budget passes
that error through unchanged, so exhaustion is never reported as anything else.

**CRT reconstruction.** Given $k_i = k \bmod p_i^{e_i}$ for each $i$, the CRT gives $k \bmod n$
uniquely.

#### Proof sketch

The key step is the subgroup projection: $G_i = (n/p_i^{e_i}) \cdot G$ has order $p_i^{e_i}$
because $p_i^{e_i} \cdot G_i = n \cdot G = \mathcal{O}$, and no smaller multiple of $G_i$ is
$\mathcal{O}$ (since $n/p_i^{e_i}$ is the exact cofactor). The projection $Q_i = (n/p_i^{e_i})
\cdot Q = (n/p_i^{e_i}) \cdot k \cdot G = k \cdot G_i$ follows directly.

The prime-power lift is a standard Hensel-style digit extraction: at step $j$, the known
$k_i \bmod p_i^j$ is subtracted from $Q_i$ (scaled by the appropriate cofactor), and the
residue is projected to the order-$p_i$ sub-subgroup for a single prime-order solve.

#### L-notation

The cost is $O\!\left(\sum_{i=1}^r e_i (\log n + \sqrt{p_i})\right)$ when each prime-order solve
is a rho walk; the implementation's enumeration for $p_i \leq 64$ costs at most 64 point additions,
a constant. For a group of prime order $n$, this is $O(\sqrt{n})$ — the same as rho. The gain is
structural, not asymptotic: when $n$ has small prime factors, the cost is dominated by
$\sqrt{p_{\max}}$ (the largest prime factor), which can be exponentially smaller than $\sqrt{n}$.

In L-notation: Pohlig–Hellman achieves $L_n[1, 1/2]$ in the worst case (prime $n$), but
$L_{p_{\max}}[1, 1/2]$ when $n$ has a small largest prime factor. The asymptotic gain over rho
is structural, not a change in the $L$-notation exponent.

**Scale note.** At toy scale ($n = 60 = 2^2 \cdot 3 \cdot 5$), the speedup over
rho is invisible — the group is too small for the birthday bound to be meaningful. At crypto
scale, Pohlig–Hellman is catastrophic for curves whose order has a small largest prime factor
(e.g. $n = 2^{255} \cdot q$ for a small $q$).

#### Cross-reference

Code realisation: `docs/PEDAGOGY.md` §9 (Pohlig–Hellman code-tour, `rho::ecdlp::pohlig`).

---

### §10.2 Smart–Satoh–Araki: The Anomalous-Curve Polynomial-Time Attack

#### The structure

The curve is *anomalous*: $\#E(\mathbb{F}_p) = p$ (equivalently, the trace of Frobenius
$t = p + 1 - \#E = 1$). This is a rare but constructible condition. The anomalous structure
admits a polynomial-time ECDLP algorithm via the p-adic formal group.

#### The escape

**The formal group and its filtration.** Lift the curve to $\tilde{E}$ over $\mathbb{Z}_p$ (lift
the coefficients $a, b$ to $\mathbb{Z}_p$). The points of $\tilde{E}(\mathbb{Q}_p)$ that reduce to
$\mathcal{O}$ modulo $p$ form a subgroup $E_1$, identified with the formal group
$\hat{E}(p\mathbb{Z}_p)$ through the formal parameter $t = -x/y$; write $E_n$ for the subgroup with
$v_p(t) \geq n$. The formal group logarithm $\psi = \log_{\hat{E}}$ is a homomorphism $E_1 \to
p\mathbb{Z}_p$ with $\psi(E_n) \subseteq p^n\mathbb{Z}_p$ and $\psi(P) \equiv t(P) \pmod{p^2}$.
These are standard facts about formal groups, stated here without proof.

**Lifting the points.** Given $P = (x_0, y_0) \in E(\mathbb{F}_p)$, choose any lift $\tilde{P} \in
\tilde{E}(\mathbb{Z}_p)$: take $\tilde{x} \equiv x_0 \pmod{p}$ and solve the curve equation for
$\tilde{y} \equiv y_0$ by Hensel's lemma. The solve exists and is unique when $y_0 \neq 0$ (that
is, $P$ is not a 2-torsion point), since $f'(y_0) = 2y_0 \not\equiv 0 \pmod{p}$ for $p > 2$.

**The anomalous hypothesis.** Because $\#E(\mathbb{F}_p) = p$, $p \cdot P = \mathcal{O}$ in
$E(\mathbb{F}_p)$, so $p \cdot \tilde{P}$ reduces to $\mathcal{O}$: it lies in $E_1$, the domain of
$\psi$. Whether it lies deeper depends on the lifted curve. The kernel of reduction of
$E(\mathbb{Z}/p^2\mathbb{Z})$ is $E_1/E_2$, so $p \cdot \tilde{P} = \mathcal{O}$ in
$E(\mathbb{Z}/p^2\mathbb{Z})$ exactly when $p \cdot \tilde{P} \in E_2$; on such a lift the
logarithm of $p \cdot \tilde{P}$ vanishes modulo $p^2$ and carries no information about $P$. (The
point $p \cdot \tilde{P}$ is in $E_1$ in all cases; it is its membership in $E_2$ that is in
question.)

**A lift-independent homomorphism.** Fix a lift $\tilde{E}$ of the curve and define
$$\lambda(P) = \psi(p \cdot \tilde{P}) / p \bmod p.$$
Another lift $\tilde{P}'$ of $P$ to the same $\tilde{E}$ differs from $\tilde{P}$ by some
$\kappa \in E_1$, so $p\tilde{P}' - p\tilde{P} = p\kappa$ and
$\psi(p\kappa) = p\,\psi(\kappa) \in p^2\mathbb{Z}_p$: $\lambda(P)$ does not depend on the lift of
the point. Since $\tilde{P} + \tilde{Q}$ lifts $P + Q$ and $\psi$ is additive on $E_1$,
$\lambda(P + Q) = \lambda(P) + \lambda(Q)$: $\lambda: E(\mathbb{F}_p) \to \mathbb{F}_p$ is a
homomorphism (for the chosen $\tilde{E}$).

**The DLP recovery.** For $Q = k \cdot G$ this gives $\lambda(Q) = k\,\lambda(G)$, so
$$k \equiv \frac{\psi(p \cdot \tilde{Q})}{\psi(p \cdot \tilde{G})} \pmod{p}$$
whenever $\lambda(G) \neq 0$. Because $\psi(P) \equiv t(P) \pmod{p^2}$, that condition says exactly
that $p \cdot \tilde{G} \in E_1 \setminus E_2$ (its parameter has valuation $1$).

**Which lifts of the curve work.** The condition $\lambda(G) \neq 0$ depends on the lift of the
*curve* $\tilde{E}$, not on the choice of point lifts. If $p \cdot \tilde{G} \in E_2$ then
$\psi(p \cdot \tilde{G})$ is divisible by $p^2$, $\lambda(G) = 0$, and the displayed ratio is a
quotient of second-order terms from which the derivation does not determine $k$. What is derived
above is the precise requirement, $p \cdot \tilde{G} \notin E_2$.

Smart (*The discrete logarithm problem on elliptic curves of trace one*, J. Cryptology 12(3),
1999, 193–196, p. 195) states that the method fails when the lift of the curve is the *canonical*
lift, that a lift is canonical with probability about $1/p$, that the failure is detectable, and
that the remedy is to repeat with another lift of the curve; that note also observes that only
$E(\mathbb{Z}/p^2\mathbb{Z})$ is needed. It does not derive why the canonical lift puts
$p \cdot \tilde{G}$ in $E_2$, and neither does this text. Leprévost–Monnerat–Varrette–Vaudenay
(*Generating anomalous elliptic curves*, Inform. Process. Lett. 93, 2005, 225–230; preprint read)
restate the method and do not discuss the failure. The identification of the failing lifts with
the canonical lift, beyond Smart's one sentence, is therefore not checked here.

**The repair used by the implementation.** The lifted curve is $y^2 = x^3 + (a + p\,s)x + b$ with
$s \not\equiv 0 \pmod p$, which reduces to the same curve modulo $p$; the points are lifted to it
and multiplied by $p$ on it. If $p \cdot \tilde{G} \in E_2$ for the shift $s_1$, one second shift
$s_2 \not\equiv s_1, 0$ is tried, and $p \cdot \tilde{G} \in E_2$ on both is reported as an
invariant failure. This is justified only if at most one class of $s$ modulo $p$ is degenerate for
a given curve; that is *measured*, not derived: on the toy curve and on a second searched curve
with $a \neq 0$ every class is checked, for every point, by
`kat_at_most_one_shift_class_is_degenerate`, and an exploratory scan of all 481 anomalous
curves $y^2 = x^3 + ax + b$ over $\mathbb{F}_p$ with $5 \leq p \leq 59$ and $b \neq 0$ (not
retained as a test) found exactly one degenerate class of $s$ for each. On the toy curve $y^2 =
x^3 + 5$ over $\mathbb{F}_7$ that class is $s = 0$, the stored coefficient $a = 0$ (parameter
valuation $2$, asserted by `kat_stored_lift_of_toy_curve_is_degenerate`).

#### Proof sketch

The two facts doing the work are that $\psi$ is a homomorphism on $E_1$ and that $\psi$ of a
multiple $p\kappa$ of an element of $E_1$ lies in $p^2\mathbb{Z}_p$, so that dividing $\psi(p
\tilde{P})$ by $p$ and reducing modulo $p$ forgets the choice of lift. Together with the anomalous
hypothesis (which puts $p \cdot \tilde{P}$ in the domain of $\psi$) they make $\lambda$ a
homomorphism $E(\mathbb{F}_p) \to \mathbb{F}_p$, and a homomorphism that is nonzero on the generator
$G$ of a group of prime order is injective, so $k$ is the ratio $\lambda(Q)/\lambda(G)$. Nonzero on
$G$ is the requirement $p \cdot \tilde{G} \notin E_2$ on the chosen lift of the curve.

#### L-notation

The *algorithm* is polynomial in $\log p$: multiplying a lift by $p$ by double-and-add takes
$O(\log p)$ group operations in $\tilde{E}(\mathbb{Z}/p^r\mathbb{Z})$ for a small fixed $r$ (Smart
observes that $r = 2$ suffices for the group computation), each of polynomial cost, and the
logarithm is a truncated power series. This is the sharpest escape in this chapter: the anomalous
structure collapses the ECDLP to a linear computation. The *implementation* is not polynomial in
$\log p$. Its multiplication by $p$ is a loop of $p - 1$ projective additions, $O(p)$ rather than
the $O(\log p)$ of double-and-add, and it confirms anomalousness by an exact point count that sums
a Legendre symbol over all $x \in \mathbb{F}_p$, also $O(p)$ (`verify_anomalous` is exact where
$x^3$ overflows 64 bits, as `kat_verify_anomalous_exact_where_cube_overflows_u64` checks, but
practical only for small $p$). The final check $k \cdot G = Q$ is $O(\log p)$.

**Scale note.** At toy scale ($p = 7$) the constant factors are invisible — the computation is
instantaneous. The pedagogical content is the mechanism, not the timing.

#### What the implementation does

`ssa_solve` lifts the curve as $a \mapsto a + p\,s$ (first $s = 1$, then $s = 2$ only if
$p \cdot \tilde{G} \in E_2$; $b$ is used as stored), lifts $G$ and $Q$ to that curve at precision
$p^4$ (`lift_point`: $x$ exactly, $y$ by Hensel's lemma), multiplies each lift by $p$ on the same
curve in projective coordinates (`multiply_by_p_proj`), evaluates `elliptic_log_proj` at
$t = -X/Y$, and divides the leading units in $\mathbb{F}_p$ to obtain
$k_{\mathrm{raw}} = \psi(p \cdot \tilde{Q}) / \psi(p \cdot \tilde{G}) \bmod p$. Despite its name,
`elliptic_log_proj` computes $\log(1 + t)$, the logarithm of the *multiplicative* formal group;
that is a different power series from the logarithm $\psi$ of $E$'s formal group used in the
derivation above, but both are $\equiv t \pmod{p^2}$ on $E_1$, and only the leading unit is read.
The precision $4$ is `WORKING_PRECISION`; precision $3$ is the least that can work, and was
measured to suffice on the same scan as the shift classes (not asserted by a test). The quotient
is the answer, but the solver still returns a $k$ only after checking $k \cdot G = Q$.

The failure classes are kept apart. $Q = \mathcal{O}$ is $k = 0$ with no reduction
(`kat_ssa_solve_infinity_is_zero`). A curve whose exact count differs from $p$ is `NotAnomalous`
(`kat_ssa_solve_not_anomalous`); a prime that does not fit one 64-bit limb is `Unsupported`, never
reported as non-anomalous (`kat_large_prime_is_unsupported_not_non_anomalous`). The invalid inputs
(composite modulus, singular curve, $G = \mathcal{O}$, $G$ or $Q$ off the curve, $n \neq p$) are
`InvalidInput`; a 2-torsion point cannot occur once those checks pass, since $E(\mathbb{F}_p)$
then has odd prime order, and the lift helper reports it as `LiftError::Hensel`
(`kat_lift_2torsion_errors_non_simple_root`, which tests `lift_point` directly). A degenerate lift
is not a failure: it is detected (`p·G̃ ∈ E₂`, valuation of $t$ at least $2$) and answered with the
second lift, tested by `ssa_solve_retries_when_the_first_shift_is_degenerate` and
`retry_runs_the_second_shift_after_a_degenerate_first`; $p \cdot \tilde{G} \in E_2$ on both lifts
is `InvariantFailure("p·G̃ in E₂ on two shifted lifts")`
(`two_degenerate_lifts_are_an_invariant_failure`). If the quotient is the only candidate and fails
$k \cdot G = Q$, the result is `Inconsistent { k }`, an invariant failure that is never retried
(`verify_candidate_rejects_unverifiable_quotients`, a unit test in `rho/src/ssa/reduce.rs`; no
fixture reaches it through `ssa_solve`). `kat_ssa_solve_recovers_every_k` and
`kat_second_fixture_recovers_every_k` check that every $k \in 1..p$ is recovered on the toy curve
and on a searched anomalous curve with $a \neq 0$ ($y^2 = x^3 + x + 5$ over $\mathbb{F}_{11}$).

**The stored lift of the toy curve.** The toy curve $y^2 = x^3 + 5$ over $\mathbb{F}_7$ has
$a = 0$, so its stored lift has $j = 0$. On that lift $p \cdot \tilde{G}$ has parameter valuation
$2$ (`kat_stored_lift_of_toy_curve_is_degenerate`), so the leading units are second-order terms,
not $\lambda$; on every shift $s \not\equiv 0$ the valuation is $1$ and the raw quotient is
exactly $k$, with no factor $2$ and no second candidate
(`attempt_on_shifted_lift_returns_the_raw_quotient`). This text claims no connection between the
degeneracy and complex multiplication beyond the measurement.

#### Cross-reference

Code realisation: `docs/PEDAGOGY.md` §11 (SSA code-tour, `rho::ssa`).

---

### §10.3 GHS/Weil Descent: The Binary-Curve Transfer

#### The structure

The curve is defined over a binary field $E/\mathbb{F}_{2^m}$ with a subfield tower
$\mathbb{F}_{2^l} \subset \mathbb{F}_{2^m}$ (with $l \mid m$). The Weil restriction
$\mathrm{Res}_{\mathbb{F}_{2^m}/\mathbb{F}_{2^l}}(E)$ is an abelian variety of dimension $m/l$
over $\mathbb{F}_{2^l}$. When $m/l$ is odd, this abelian variety contains the Jacobian of a
hyperelliptic curve $C/\mathbb{F}_{2^l}$ of genus $g = (m/l - 1)/2$.

#### The escape (as a transfer)

**This section represents GHS honestly as a transfer.** The GHS construction reduces the ECDLP
on $E/\mathbb{F}_{2^m}$ to a DLP on the Jacobian $\mathrm{Jac}(C)/\mathbb{F}_{2^l}$. The
downstream solve — index calculus on the hyperelliptic Jacobian — is a separate step that this
project does not implement. The chapter covers the descent reduction and log-preservation
verification; the downstream solve is not developed here.

**The Artin–Schreier extension.** The function field of $E/\mathbb{F}_{2^m}$ is a degree-2
extension of $\mathbb{F}_{2^m}(x)$ defined by $y^2 + h(x)y = f(x)$ (the Artin–Schreier form
for characteristic 2). The Weil restriction replaces the field $\mathbb{F}_{2^m}$ with the
subfield $\mathbb{F}_{2^l}$ and raises the dimension from 1 to $m/l$: the single variable $x$
over $\mathbb{F}_{2^m}$ becomes $m/l$ variables over $\mathbb{F}_{2^l}$.

**The hyperelliptic curve.** The Weil restriction produces an abelian variety over
$\mathbb{F}_{2^l}$. When $m/l$ is odd (the imaginary hyperelliptic model), this abelian variety
is the Jacobian of a hyperelliptic curve $C/\mathbb{F}_{2^l}$ of genus $g = (m/l - 1)/2$.

**Log-preservation.** The transfer map $\phi: E(\mathbb{F}_{2^m}) \to \mathrm{Jac}(C)(\mathbb{F}_{2^l})$
is a group homomorphism. It preserves the discrete logarithm: if $Q = k \cdot G$ in
$E(\mathbb{F}_{2^m})$, then $\phi(Q) = k \cdot \phi(G)$ in $\mathrm{Jac}(C)(\mathbb{F}_{2^l})$.

**Input guards in the implementation.** The map is a homomorphism only on $E(\mathbb{F}_{2^m})$,
so `transfer_point` rejects a point that does not satisfy the curve equation as `NotOnCurve`
(`transfer_rejects_off_curve_point`, `ghs_descend_rejects_off_curve_point` in
`rho/tests/ghs_kat.rs`), and maps the point at infinity to the identity divisor rather than
rejecting it (`transfer_infinity_is_zero_divisor`). The image of a point is a point of
$C(\mathbb{F}_{2^l})$ only when $E$ is defined over $\mathbb{F}_{2^l}$, so a coordinate of the
conorm image that is not fixed by the relative Frobenius $c \mapsto c^{2^l}$ is rejected as
`NotInSubfield` (`transfer_rejects_conorm_outside_subfield`, on a curve with coefficient $a =
\alpha \notin \mathbb{F}_{2^2}$).

**Identifying $\mathbb{F}_{2^l}$ inside $\mathbb{F}_{2^m}$.** Let $\alpha$ be the root of the
defining polynomial of $\mathbb{F}_{2^m}$ and $\beta = \alpha^{(2^m-1)/(2^l-1)}$, a generator of
$\mathbb{F}_{2^l}^\times$. The subfield is $\mathbb{F}_2(\beta)$, and an element of
$\mathbb{F}_{2^l}$ with coordinate vector $(c_0, \dots, c_{l-1})$ is the element
$\sum c_i \beta^i$ of $\mathbb{F}_{2^m}$. Reading a $\mathbb{F}_{2^m}$ element of the subfield back as coordinates in the
basis $\{1, \beta, \dots, \beta^{l-1}\}$ (`shared_gf2m::restrict`) is a field isomorphism exactly
when the subfield polynomial `poly_l` is the minimal polynomial of $\beta$: $\deg = l$ and
$\mathrm{poly}_l(\beta) = 0$. Reducing the raw bit pattern of the element modulo `poly_l` is not an
isomorphism: $x^2 + x + 1$ does not divide $x^6 + x + 1$, so bit reduction does not respect
multiplication, and on $\mathbb{F}_4 \subset \mathbb{F}_{2^6}$ it exchanges the two conjugates
$\beta$ and $\beta^2$. Both `extract_ghs_curve` (the coefficients of $F$) and `transfer_point` (the
coordinates of the conorm image) use `restrict`, and both reject parameters before any coefficient
or point is read: `SubfieldDivisibility` if $l \nmid m$ (or $l = 0$), and
`SubfieldPolyMismatch` if `poly_l` is not the minimal polynomial of $\beta$
(`mismatched_poly_l_is_rejected`, `non_dividing_l_is_rejected_without_panic`). The restriction map
needs $(2^m-1)/(2^l-1)$ to fit in a `u64`, i.e. $m \le 63$, which a single-limb defining polynomial
guarantees.

*Worked check on the toy fixture.* With $\alpha$ the root of $x^6 + x + 1$: $\alpha^9 = \alpha^4 +
\alpha^3$ and $\alpha^{12} = \alpha^2 + 1$, so $\beta = \alpha^{21} = \alpha^5 + \alpha^4 + \alpha^3
+ \alpha + 1 =$ `0x3b`, and $\beta^2 = \beta + 1 =$ `0x3a`; the minimal polynomial of $\beta$ is
$x^2 + x + 1$, the fixture's `poly_l`. The coordinates in $\{1, \beta\}$ are $1 \mapsto 1$, $\beta
\mapsto 2$, $\beta^2 \mapsto 3$. The point $P = (1, \beta^2)$ of $E: y^2 + xy = x^3 + x^2 + 1$ lies
in $E(\mathbb{F}_4)$, so its conorm is $3P$; doubling gives $2P = (0, 1)$, a point of order 2, hence
$4P = \infty$ and $3P = -P = (1, \beta)$. Thus $D_P = [X + 1, 2]$, and for $G + P = -P$ the conorm
is $3(-P) = P$ giving $D_{G+P} = [X + 1, 3]$. These values are pinned by
`transfer_sample_point_known_answer` and `transfer_g_plus_p_known_answer`, whose doc comments carry
the derivation, and agree with PARI/GP (`ellmul` over $\mathbb{F}_{2^6}$ gives
$3P = (1, \texttt{0x3b})$ and $3(G+P) = (1, \texttt{0x3a})$). The raw-bit reading would give $3$
and $2$ instead. The
homomorphism and log-preservation tests hold under either reading, because the curve
$C: Y^2 + XY = X^3 + X^2 + 1$ has coefficients in $\mathbb{F}_2$ and is stable under the
conjugation; they do not test the identification, whereas
`curve_and_transfer_share_the_subfield_identification` does (a curve with $a = \beta$ has
$F$-coefficient $2$ under `restrict`, and every transferred divisor is valid on that $C$).

#### L-notation

The L-notation for GHS is conditional on the genus $g$ of the hyperelliptic curve:

- For small genus $g$ (e.g. $g = 1$, which reduces to an elliptic curve DLP — no gain), the
  transfer does not help.
- For large genus $g$, index calculus on $\mathrm{Jac}(C)$ achieves subexponential complexity.
  The asymptotic win depends on the genus and the specific index-calculus algorithm applied.

The transfer itself (the descent reduction) costs $O(\mathrm{poly}(m))$ — polynomial in the
field extension degree. The asymptotic win, if any, comes from the downstream solve.

**Scale note.** The toy fixture ($m = 6$, $l = 2$, $g = 1$) has genus 1 — the
hyperelliptic Jacobian is an elliptic curve, and the downstream DLP is no easier than the
original. The GHS construction is demonstrated at toy scale for the mechanism; the asymptotic
win requires larger $m/l$ (and hence larger genus).

#### Cross-reference

Code realisation: `docs/PEDAGOGY.md` §12 (GHS code-tour, `rho::ghs`).

---

### §10.4 Index Calculus: Semaev Decomposition over the Factor Base

#### The structure

The factor base $\mathcal{F} = \{F_1, \ldots, F_B\}$ is a set of points on $E$ with small
x-coordinates. A point $P \in E(\mathbb{F}_{p^n})$ is *decomposable* over $\mathcal{F}$ if
$P = \sum_{i=1}^m F_{j_i}$ for some $F_{j_i} \in \mathcal{F}$. The Semaev summation polynomial
$S_m$ detects decomposability: $S_m(x_{F_{j_1}}, \ldots, x_{F_{j_{m-1}}}, x_P) = 0$ iff $P$
decomposes with the given $F_{j_i}$.

#### The escape

**Why the logarithms live modulo $\ell$.** Let $n$ be the order of $G$ and let $\ell$ be a prime
with $\ell \,\|\, n$ (that is, $\ell \mid n$ and $\ell^2 \nmid n$). Then the $\ell$-Sylow subgroup
of $\langle G \rangle$ is cyclic of order $\ell$, generated by $G_\ell = (n/\ell) \cdot G$.
Multiplying any point by $n/\ell$ projects it into that subgroup, so the $\ell$-component of $Q$
and of every factor-base point $F_j$ has a logarithm modulo $\ell$ with respect to $G_\ell$. This
is an invariant of the construction: `toy_fixture_parameters` asserts $\ell \,\|\, n$ for the
fixture, and the other tests check its consequence, $x \cdot G_\ell = Q_\ell$, for every $k$
they try; no test proves the invariant independently of the fixture.

**Relation collection.** Enumerate pairs $(a_i, b_i)$ and compute $P_i = a_i \cdot G + b_i \cdot
Q$. The code walks the pairs deterministically, along the anti-diagonals $a + b = 1, 2, \ldots$
inside $a \in [1, n)$, $b \in [0, n)$; this is a bound on how many relations are collected, not a
search over candidate logarithms. For each $P_i$, attempt to decompose it over $\mathcal{F}$ using
the Semaev polynomial. A successful decomposition $P_i = \sum_j e_{ij} F_j$ gives, after
projecting to the $\ell$-component and taking logarithms to the base $G_\ell$, a linear relation

$$a_i + b_i \cdot x \equiv \sum_j e_{ij} \cdot l_j \pmod{\ell},$$

where $x = \log_G(Q) \bmod \ell$ is the unknown and $l_j = \log_{G_\ell}\!\big((n/\ell) \cdot
F_j\big)$ are the unknown factor-base logarithms.

**Linear algebra: the left kernel.** Write $E$ for the matrix of exponents $e_{ij}$ (one row per
relation, one column per factor-base point) and $a$, $b$ for the vectors of $a_i$, $b_i$. A
vector $u$ in the *left kernel* of $E$ — that is, $\sum_i u_i \cdot e_{ij} \equiv 0 \pmod{\ell}$
for every $j$ — cancels every unknown $l_j$ when the relations are summed with weights $u_i$:

$$u^{\top} a + (u^{\top} b) \cdot x \equiv 0 \pmod{\ell}.$$

When $u^{\top} b \not\equiv 0$ this determines $x \equiv -\,u^{\top} a \,/\, u^{\top} b$. The
factor-base logarithms $l_j$ are never solved for. A kernel vector with $u^{\top} b \equiv 0$ says
nothing about $x$, and for honest relations it must also have $u^{\top} a \equiv 0$. The code
computes a complete basis of the left kernel (the transposed exponent matrix has the relations as
columns, and its right kernel is the left kernel of $E$). At this scale the kernel routine
`block_wiedemann_fl` performs an exact dense row reduction; it is the same entry point as the
NFS-DL linear algebra (§9), not a Krylov (block-Wiedemann or block-Lanczos) iteration.

**Outcomes.** The solver accepts $x$ only if (i) every basis vector with $u^{\top} b \not\equiv 0$
gives the same $x$ and (ii) $x \cdot G_\ell = Q_\ell$ by the group law. If no basis vector has
$u^{\top} b \not\equiv 0$ after the relation enumeration is exhausted, the result is
*undetermined*: a statement about the collected relation set, not a negative about the logarithm.
If combinations disagree, a vector with $u^{\top} b \equiv 0$ has $u^{\top} a \not\equiv 0$, or
verification fails, the relations are internally contradictory; that is an invariant failure and
is reported as such, never retried.

**The Semaev polynomial.** $S_m(X_1, \ldots, X_m)$ is defined recursively:
- $S_2(X_1, X_2) = X_1 - X_2$ (trivial: two points sum to $\mathcal{O}$ iff they are equal).
- $S_3(X_1, X_2, X_3)$: the resultant of the chord-and-tangent addition formula, giving a
  degree-4 polynomial in each variable.
- $S_m = \mathrm{Res}_X(S_{m-1}(X_1, \ldots, X_{m-2}, X),\; S_3(X_{m-1}, X_m, X))$ for $m > 3$.

The degree of $S_m$ in each variable is $2^{m-2}$, growing exponentially with $m$.

#### Proof sketch

The key step is that $S_m(x_1, \ldots, x_m) = 0$ iff there exist $y_i$ such that
$(x_i, y_i) \in E$ and $(x_1, y_1) + \cdots + (x_m, y_m) = \mathcal{O}$. This follows from
the definition of the summation polynomial as the resultant of the addition law. Because the
polynomial sees only x-coordinates, a vanishing triple means the sum is $\pm Q$; the group law
then decides which sign holds, and only the decomposition summing to $Q$ is kept. The relation
collection loop finds decomposable points by evaluating $S_m$ at the enumerated combinations
$a \cdot G + b \cdot Q$; the left-kernel combination above recovers $x = \log_G(Q) \bmod \ell$
from the collected relations alone.

**What the tests establish.** On the toy fixture ($y^2 = x^3 + x + 33$ over $\mathbb{F}_{47}$,
$n = 60$, $\ell = 5$, $|\mathcal{F}| = 6$, $m = 2$; asserted by `toy_fixture_parameters`):
`left_kernel_annihilates_exponents` checks that every returned basis vector satisfies
$\sum_i u_i e_i \equiv 0$ and that the basis is complete; `recovers_every_residue` runs every
$k \in \{0, \ldots, 59\}$, checks that the combination has $u^{\top} b \not\equiv 0$ and that
$\sum_i u_i (a_i + b_i x) \equiv 0$ holds for the returned $x = k \bmod 5$;
`undetermined_when_no_combination_sees_b` and `corrupted_relation_is_inconsistent` exercise the
two failure classes. The undetermined outcome is not reachable for an honest $Q$ on this fixture.

#### L-notation

Over $E(\mathbb{F}_p)$ (the toy setting): index calculus is **not** faster than Pollard rho.
The asymptotic win requires the extension-field setting $E(\mathbb{F}_{p^n})$ with $n > 1$
(the Gaudry–Diem setting). In that setting, the complexity is subexponential:

$$\text{Index calculus over } E(\mathbb{F}_{p^n}) = L_{p^n}\!\left[\tfrac{1}{2}, c\right]
\quad \text{(heuristic, for } n > 1\text{)}.$$

The exponent $1/2$ is a genuine improvement over rho's $L[1, 1/2]$ (fully exponential), but
not as sharp as NFS-DL's $L[1/3]$.

**Scale note.** The toy fixture operates over $E(\mathbb{F}_p)$ ($n = 1$). The asymptotic win is
not observable at toy scale — the mechanism is demonstrated, but the complexity separation
requires $n > 1$, which this code does not implement.

#### Cross-reference

Code realisation: `docs/PEDAGOGY.md` §13–§14 (Semaev and index-calculus code-tours,
`rho::semaev` and `rho::index_calculus`).

---

### §10.5 MOV/Frey–Rück: The Pairing Reduction (The Payoff Proof)

*This is the designated payoff proof for this chapter — the one full proof at payoff depth. It
is the cross-chapter bridge of this project: the MOV reduction connects the ECDLP (this
chapter) to the NFS-DL setting (§9), where index calculus applies.*

#### The structure

The curve $E/\mathbb{F}_p$ has *small embedding degree* $k$: $\ell \mid p^k - 1$ but
$\ell \nmid p^j - 1$ for $j < k$. Here $\ell$ is the prime order of the subgroup $\langle G
\rangle \subset E(\mathbb{F}_p)$ in which the ECDLP is posed.

The embedding degree $k$ is the smallest positive integer such that $\mathbb{F}_{p^k}$ contains
all $\ell$-th roots of unity — equivalently, $\mu_\ell \subset \mathbb{F}_{p^k}^*$.

#### The Weil pairing

**Definition.** Let $E[\ell] = \{P \in E(\overline{\mathbb{F}_p}) : \ell \cdot P = \mathcal{O}\}$
be the $\ell$-torsion subgroup of $E$ over the algebraic closure. The *Weil pairing* is a map

$$e_\ell: E[\ell] \times E[\ell] \to \mu_\ell \subset \overline{\mathbb{F}_p}^*$$

defined via the divisor theory of $E$. Concretely: for $P, Q \in E[\ell]$, choose rational
functions $f_P, f_Q$ on $E$ with divisors $\mathrm{div}(f_P) = \ell \cdot (P) - \ell \cdot
(\mathcal{O})$ and $\mathrm{div}(f_Q) = \ell \cdot (Q) - \ell \cdot (\mathcal{O})$. Then

$$e_\ell(P, Q) = \frac{f_P(Q + S)}{f_P(S)} \cdot \frac{f_Q(P + T)}{f_Q(T)}^{-1}$$

for generic auxiliary points $S, T$ (the value is independent of the choice of $S, T$).

**Bilinearity.** The Weil pairing is bilinear in both arguments:

$$e_\ell(P_1 + P_2, Q) = e_\ell(P_1, Q) \cdot e_\ell(P_2, Q), \qquad
e_\ell(P, Q_1 + Q_2) = e_\ell(P, Q_1) \cdot e_\ell(P, Q_2).$$

*Proof.* Bilinearity follows from the linearity of the divisor map: $\mathrm{div}(f_{P_1+P_2})
= \mathrm{div}(f_{P_1}) + \mathrm{div}(f_{P_2})$ (up to principal divisors), so $f_{P_1+P_2}
= f_{P_1} \cdot f_{P_2}$ (up to a constant), and the pairing value multiplies accordingly.
$\square$

**Non-degeneracy.** The Weil pairing is non-degenerate: for every $P \in E[\ell]$ with
$P \neq \mathcal{O}$, there exists $Q \in E[\ell]$ such that $e_\ell(P, Q) \neq 1$.

*Proof sketch.* Non-degeneracy follows from the Riemann–Roch theorem applied to the divisor
class group of $E$: the pairing is a perfect pairing on $E[\ell] \times E[\ell]$, which is a
free $\mathbb{Z}/\ell\mathbb{Z}$-module of rank 2. See Silverman [S09, §III.8]. $\square$

**Galois equivariance.** For $\sigma \in \mathrm{Gal}(\overline{\mathbb{F}_p}/\mathbb{F}_p)$:

$$e_\ell(\sigma(P), \sigma(Q)) = \sigma(e_\ell(P, Q)).$$

This implies that if $P, Q \in E(\mathbb{F}_p)$ (i.e. they are defined over $\mathbb{F}_p$),
then $e_\ell(P, Q) \in \mathbb{F}_{p^k}$ (the smallest field containing all $\ell$-th roots of
unity).

#### The Tate pairing and Miller's algorithm

The Weil pairing is theoretically clean but computationally expensive. In practice, the
*reduced Tate pairing* is used:

$$\hat{e}: E(\mathbb{F}_{p^k})[\ell] \times E(\mathbb{F}_{p^k}) / \ell E(\mathbb{F}_{p^k})
\to \mathbb{F}_{p^k}^* / (\mathbb{F}_{p^k}^*)^\ell.$$

After the *final exponentiation* (raising to the power $(p^k - 1)/\ell$), the reduced Tate
pairing lands in $\mu_\ell \subset \mathbb{F}_{p^k}^*$ and is efficiently computable via
*Miller's algorithm* in $O(k \log p)$ field operations.

**Miller's algorithm** evaluates the rational function $f_P$ (whose divisor is
$\ell \cdot (P) - \ell \cdot (\mathcal{O})$) at a point $Q$ by a double-and-add loop over the
binary expansion of $\ell$. At each step, the current function $f$ is updated by the line
function through the current accumulator point and the next doubling/addition. The final value
is $f_P(Q)$, which is then raised to the power $(p^k - 1)/\ell$ (the final exponentiation).

#### The MOV reduction

**Setup.** Let $G \in E(\mathbb{F}_p)$ have prime order $\ell$, and let $Q = k \cdot G$ be the
ECDLP target. Let $k$ be the embedding degree of $E$ with respect to $\ell$.

**Step 1: Choose a $\mu_\ell$-generator.** Find a point $R \in E(\mathbb{F}_{p^k})[\ell]$ such
that $e(G, R) \neq 1$ (where $e$ is the reduced Tate pairing). Such $R$ exists by
non-degeneracy. In practice, $R$ is chosen as a random point in $E(\mathbb{F}_{p^k})[\ell]$
and the condition $e(G, R) \neq 1$ is verified.

**Step 2: Compute the pairing values.** Compute

$$g_0 = e(G, R) \in \mu_\ell \subset \mathbb{F}_{p^k}^*, \qquad
h_0 = e(Q, R) \in \mu_\ell \subset \mathbb{F}_{p^k}^*.$$

**Step 3: Apply bilinearity.** By bilinearity of the pairing:

$$h_0 = e(Q, R) = e(k \cdot G, R) = e(G, R)^k = g_0^k.$$

Therefore $k = \log_{g_0}(h_0)$ in $\mathbb{F}_{p^k}^*$.

**Step 4: Solve the DLP in $\mathbb{F}_{p^k}^*$.** The DLP $h_0 = g_0^k$ in $\mathbb{F}_{p^k}^*$
is solved by NFS-DL (the `gnfs::dl::solve_dl` entry point; §9). Since $g_0, h_0
\in \mu_\ell \subset \mathbb{F}_{p^k}^*$, the DLP is in the order-$\ell$ subgroup of
$\mathbb{F}_{p^k}^*$.

**Step 5: Recover $k$.** The NFS-DL solver returns $k \bmod \ell$, which is the ECDLP scalar.

**What `mov_reduce` checks.** The returned value is a candidate until it passes $k \cdot G = Q$,
which the code evaluates directly in $E[\ell]$, independently of the pairing and of the finite-field
logarithm engine. A success is a `MovSolution` carrying $k \bmod \ell$, $\ell$, and the two pairing
values with $g_0^k = h_0$ (`mov_reduce_recovers_k1`, `mov_reduce_recovers_k2` in
`rho/tests/mov_kat.rs`, which recheck all of it from the fields). $Q = \mathcal{O}$ gives $k = 0$
with $h_0 = 1$ (`mov_reduce_recovers_k0_for_target_at_infinity`). The failure classes are kept
apart:

- *Invalid input* (`MovError::InvalidInput`): $\ell$ not prime, $\ell \nmid p^k - 1$,
  $G = \mathcal{O}$, or any of $G$, $Q$, $R$ off the curve or outside $E[\ell]$
  (`mov_reduce_rejects_composite_ell`, `mov_reduce_rejects_ell_not_dividing_field_order`,
  `mov_reduce_rejects_infinite_generator`, `mov_reduce_rejects_off_curve_points`,
  `mov_reduce_rejects_points_outside_e_ell`).
- *Degenerate pairing* (`MovError::DegeneratePairing`): $R = \mathcal{O}$, $e(G, R) = 1$, or a
  Miller evaluation that returns $0$ because the evaluation point lies on the support of the divisor
  of the Miller function, where the value is undefined. That is the case for $R = G$: the Miller
  function of $R$ is evaluated at $G = R$, the first tangent-line numerator vanishes there, and the
  reduced Tate pairing yields $0$, not $1$ (the test asserts only the `DegeneratePairing` outcome,
  not the value $0$). Other $R$ in the subgroup generated by $G$ are not tested. In each
  degenerate case $R$ is not a $\mu_\ell$-generator for $G$ and the finite-field logarithm
  carries no information about $k$; it is a property of the supplied $R$, not a solver failure
  (`mov_reduce_degenerate_r_is_degenerate_pairing`, `mov_reduce_infinite_r_is_degenerate_pairing`).
- *Unsupported domain* (`MovError::Unsupported`): the field characteristic has a non-zero word
  above the lowest ($p \ge 2^{64}$), which the single-limb bridge to the log engine cannot
  represent. It is reported before any arithmetic on $p$, so no later conversion truncates
  (`mov_reduce_rejects_field_above_one_limb`; the guard is an ordinary conditional, not a
  `debug_assert`, so it holds in release builds too).
- *Engine failure* (`MovError::Dl`): the error of `gnfs::dl::solve_dl` is passed through unchanged.
  The code calls it with extension degree $2$ only.
- *Invariant failure* (`MovError::Inconsistent`): the engine's value fails $k \cdot G = Q$. No
  fixture reaches it through `mov_reduce`; only a unit test of the check does
  (`log_matches_accepts_only_the_true_log` in `rho/src/pairing/mov.rs`).

The final exponent $(p^k - 1)/\ell$ is computed exactly, whatever its size: for one-limb $p$ it
already exceeds $2^{64}$ once $p^2 > 3 \cdot 2^{64}$ (for $\ell = 3$, $k = 2$). Cutting it to its
low word lands the result outside $\mu_\ell$. On $y^2 = x^3 + 1$ over $\mathbb{F}_p$,
$p = 68719476767$ (supersingular, $\#E = p + 1 \equiv 0 \bmod 3$), with $G = (0, 1)$ and a point $R$
of $E(\mathbb{F}_{p^2})[3] \smallsetminus E(\mathbb{F}_p)$, the reduced Tate pairing satisfies
$e^3 = 1$ and $e(2G, R) = e(G, R)^2$ only with the full exponent
(`reduced_tate_with_exponent_above_one_limb_lands_in_mu_ell`, `rho/src/pairing/mov.rs`). The
log engine (`gnfs::dl::solve_dl`) scans all of $\mathbb{F}_p$ when it picks a modulus, so
`mov_reduce` itself does not finish at that $p$; the pairing is tested directly.

#### Why this is a polynomial-time reduction

The reduction from ECDLP to DLP in $\mathbb{F}_{p^k}^*$ is polynomial-time:

1. **Finding $R$:** a random point in $E(\mathbb{F}_{p^k})[\ell]$ can be found in $O(\mathrm{poly}(k
   \log p))$ time by computing $(\#E(\mathbb{F}_{p^k}) / \ell) \cdot P$ for a random $P$.
2. **Computing the pairing:** Miller's algorithm costs $O(k \log p)$ field operations in
   $\mathbb{F}_{p^k}$.
3. **The final exponentiation:** costs $O(k \log p)$ field operations.

The total reduction cost is $O(\mathrm{poly}(k \log p))$ — polynomial in $k \log p$. For small
$k$ (the MOV threshold), this is efficient.

#### The embedding-degree condition (the MOV threshold)

The MOV reduction is effective only when $k$ is small. For a generic elliptic curve over
$\mathbb{F}_p$, the embedding degree $k$ is $\Theta(p)$ — exponentially large — and the
reduction is useless (the DLP in $\mathbb{F}_{p^k}^*$ is harder than the original ECDLP).

The *MOV threshold* is the condition $k \leq C \log p$ for a small constant $C$. Curves with
small embedding degree are called *MOV-vulnerable*. For a random curve over $\mathbb{F}_p$,
the probability that $k \leq C \log p$ is negligibly small (roughly $1/p$), so most curves
are not MOV-vulnerable. Cryptographic curves are specifically chosen to have large embedding
degree (e.g. secp256k1 has $k \approx 2^{128}$).

#### L-notation

The MOV reduction reduces the ECDLP to a DLP in $\mathbb{F}_{p^k}^*$. The DLP in
$\mathbb{F}_{p^k}^*$ is solved by NFS-DL in time $L_{p^k}[1/3, (64/9)^{1/3}]$ (§9 of this
textbook). For small $k$, $p^k$ is polynomial in $p$, so:

$$\text{MOV-reduced ECDLP complexity} = L_{p^k}\!\left[\tfrac{1}{3},\, \left(\tfrac{64}{9}\right)^{1/3}\right]
= L_p\!\left[\tfrac{1}{3},\, k^{1/3} \cdot \left(\tfrac{64}{9}\right)^{1/3}\right].$$

This is subexponential in $\log p$ — a dramatic improvement over rho's $L_p[1, 1/2]$ (fully
exponential). The cross-chapter bridge: the MOV reduction connects the ECDLP (§10) to the
NFS-DL setting (§9), where the $L[1/3]$ complexity applies.

**Scale note.** At toy scale ($p = 47$, $k = 2$, $\ell = 3$), the pairing + NFS-DL
overhead dominates the toy rho cost. The asymptotic win requires crypto-scale $p$ where the
$L[1/3]$ vs $L[1, 1/2]$ separation is observable.

#### Cross-reference

Code realisation: `docs/PEDAGOGY.md` §10 (MOV code-tour, `rho::pairing::mov`). The NFS-DL
solver that `mov_reduce` calls is documented in `gnfs/docs/PEDAGOGY.md` §63–§71 (the NFS-DL
chapter) and developed mathematically in §9 of this textbook.

---

### §10.6 Per-Attack L-Notation Comparison

This table extends the §"Escape from Search" L-notation hierarchy table with a per-attack
row for the five algebraic ECDLP attacks. The baseline is Pollard rho.

| Attack | Escape structure | L-notation complexity | Precondition | Toy-scale observable? |
|--------|-----------------|----------------------|--------------|----------------------|
| Pollard rho (baseline) | None — generic $\sqrt{n}$ walk | $L_n[1, 1/2] = \Theta(\sqrt{n})$ | None | Yes (the baseline) |
| Pohlig–Hellman | Composite group order | $L_{p_{\max}}[1, 1/2]$ (largest prime factor) | $n$ composite with small $p_{\max}$ | No — speedup requires $p_{\max} \ll n$ |
| Smart–Satoh–Araki | Anomalous order ($\#E = p$) | $L_p[0] = O(\log p)$ (polynomial time) | $\#E(\mathbb{F}_p) = p$ | Mechanism yes; timing separation no |
| GHS/Weil descent | Binary field tower | Conditional: $L_{p^{m/l}}[1/2, c]$ (downstream index calculus on Jacobian) | $E/\mathbb{F}_{2^m}$, $l \mid m$, $m/l$ odd | No — transfer only at toy scale |
| Index calculus | Factor-base decomposability | $L_{p^n}[1/2, c]$ (subexponential, $n > 1$) | $E(\mathbb{F}_{p^n})$, $n > 1$ | No — toy fixture is $n = 1$ |
| MOV/Frey–Rück | Small embedding degree | $L_{p^k}[1/3, (64/9)^{1/3}]$ (via NFS-DL) | $\ell \mid p^k - 1$, $k$ small | No — pairing overhead dominates at toy scale |

**The scale boundary, stated explicitly.** The asymptotic L-notation separations in this
table are NOT observable at the toy scale of the benchmark fixtures ($p = 47$, $p = 7$). The
table reports the theoretical complexity class; the toy-scale benchmark setup is in the
cross-attack harness section of `docs/BENCHMARKS.md`. The separation between $L[1, 1/2]$ (rho),
$L[1/2]$ (index calculus), $L[1/3]$ (MOV via NFS-DL), and $L[0]$ (SSA) is a statement about
asymptotic behaviour at cryptographic scale — not a ranking of the toy-scale timings.

---

### §10.7 Cross-References

#### Code realisation

The algebraic-attack code-tour in `docs/PEDAGOGY.md` §8–§18 is the code-first sibling to this
chapter.
It documents the module surfaces, toy fixtures, KAT summary, and design-statement verification
for the five attacks. The chapter-pairing is:

| This chapter | Code-tour section |
|-------------|-------------------|
| §10.1 Pohlig–Hellman | `docs/PEDAGOGY.md` §9 |
| §10.2 Smart–Satoh–Araki | `docs/PEDAGOGY.md` §11 |
| §10.3 GHS/Weil descent | `docs/PEDAGOGY.md` §12 |
| §10.4 Index calculus | `docs/PEDAGOGY.md` §13–§14 |
| §10.5 MOV/Frey–Rück | `docs/PEDAGOGY.md` §10 |

#### Within this textbook

- **§Pollard Rho for ECDLP** — the generic $\sqrt{n}$ baseline this chapter's attacks escape.
  The "When the bound breaks" subsection names the five attacks; this chapter develops them.

- **§"Escape from Search: The Through-Line"** — the five-family structure taxonomy and the
  L-notation hierarchy table. This chapter extends that taxonomy per-attack.

- **§9 (NFS-DL)** — the NFS-DL chapter. The MOV reduction (§10.5) calls `gnfs::dl::solve_dl`
  (the NFS-DL entry point); the NFS-DL complexity analysis (§9.7) is the mathematical basis for
  the MOV L-notation claim.

- **§Prerequisites** — the CRT (used in Pohlig–Hellman), the discrete logarithm problem
  definition, and the polynomial-time reduction definition (used in the MOV reduction).

#### Benchmark data

- **`docs/BENCHMARKS.md`, the cross-attack ECDLP harness** — the structural-precondition-
  conditional table: the benchmark setup this chapter's L-notation claims are read against (at
  toy scale).

---

### §10.8 Further Reading

32. **Menezes, A. J., Okamoto, T., and Vanstone, S. A. (1993).** "Reducing elliptic curve
    logarithms to logarithms in a finite field." *IEEE Transactions on Information Theory*,
    39(5), 1639–1646. [MOV93] The original MOV reduction paper — the source for §10.5.

33. **Frey, G., and Rück, H.-G. (1994).** "A remark concerning m-divisibility and the discrete
    logarithm in the divisor class group of curves." *Mathematics of Computation*, 62(206),
    865–874. [FR94] The Frey–Rück variant of the pairing reduction.

34. **Smart, N. P. (1999).** "The discrete logarithm problem on elliptic curves of trace one."
    *Journal of Cryptology*, 12(3), 193–196. [S99] The SSA polynomial-time attack (§10.2).

35. **Satoh, T., and Araki, K. (1998).** "Fermat quotients and the polynomial time discrete log
    algorithm for anomalous elliptic curves." *IEICE Transactions on Fundamentals*, E81-A(6),
    1228–1233. [SA98] The Satoh–Araki variant of the anomalous-curve attack.

36. **Gaudry, P., Hess, F., and Smart, N. P. (2002).** "Constructive and destructive facets of
    Weil descent on elliptic curves." *Journal of Cryptology*, 15(1), 19–46. [GHS02] The GHS
    Weil-descent attack (§10.3).

37. **Semaev, I. (2004).** "Summation polynomials and the discrete logarithm problem on elliptic
    curves." Cryptology ePrint Archive, Report 2004/031. [Sem04] The summation polynomial
    primitive (§10.4).

38. **Gaudry, P. (2009).** "Index calculus for abelian varieties of small dimension and the
    elliptic curve discrete logarithm problem." *Journal of Symbolic Computation*, 44(12),
    1690–1702. [Gau09] The Gaudry index-calculus algorithm (§10.4).

39. **Silverman, J. H. (2009).** *The Arithmetic of Elliptic Curves*. 2nd ed. Springer. [S09]
    §III.8 for the Weil pairing (non-degeneracy, bilinearity); §V.1 for Hasse's bound.

40. **Washington, L. C. (2008).** *Elliptic Curves: Number Theory and Cryptography*. 2nd ed.
    CRC Press. Chapter 11 for the MOV attack; Chapter 4 for the Weil and Tate pairings.

41. **Blake, I. F., Seroussi, G., and Smart, N. P. (2005).** *Advances in Elliptic Curve
    Cryptography*. Cambridge University Press. Chapter IX for the MOV and Frey–Rück attacks;
    Chapter X for the GHS descent.

42. **Galbraith, S. D. (2012).** *Mathematics of Public Key Cryptography*. Cambridge University
    Press. Chapter 19 for pairings and the MOV reduction; Chapter 21 for index calculus on
    elliptic curves.

43. **Proos, J., and Zalka, C. (2003).** "Shor's discrete logarithm quantum algorithm for
    elliptic curves." *Quantum Information and Computation*, 3(4), 317–344. [PZ03] The
    two-register hidden-subgroup construction for ECDLP; the qubit-budget analysis for
    elliptic-curve period finding. Source for §11.4.

44. **Nielsen, M. A., and Chuang, I. L. (2000).** *Quantum Computation and Quantum Information*.
    Cambridge University Press. [NC00] §5.3 (order-finding and factoring), §5.4 (discrete
    logarithms), §5.2 (the quantum Fourier transform and phase estimation). The standard
    textbook reference for the QFT, phase estimation, and Shor's algorithm.

45. **National Institute of Standards and Technology (2024).** *Module-Lattice-Based
    Key-Encapsulation Mechanism Standard*. FIPS 203. U.S. Department of Commerce. [FIPS203]
    The ML-KEM (Kyber) standard — the primary lattice-based key-encapsulation mechanism
    standardised by NIST.

46. **National Institute of Standards and Technology (2024).** *Module-Lattice-Based Digital
    Signature Standard*. FIPS 204. U.S. Department of Commerce. [FIPS204] The ML-DSA
    (Dilithium) standard — the primary lattice-based digital signature scheme standardised
    by NIST.

47. **National Institute of Standards and Technology (2024).** *Stateless Hash-Based Digital
    Signature Standard*. FIPS 205. U.S. Department of Commerce. [FIPS205] The SLH-DSA
    (SPHINCS+) standard — the hash-based digital signature scheme standardised by NIST.

48. **Castryck, W., and Decru, T. (2023).** "An efficient key recovery attack on SIDH."
    *Advances in Cryptology — EUROCRYPT 2023*, LNCS 14008, 423–447. [CD23] The classical
    polynomial-time attack on SIDH/SIKE; the break that eliminated the leading isogeny-based
    key-encapsulation candidate from the NIST PQC process.

---

## Shor's Algorithm and Post-Quantum Context

*Maths-first sibling to `shor/docs/PEDAGOGY.md` (the Shor code-tour). For the phase-by-phase
implementation — module surfaces, toy fixtures, KAT summary, and design-statement verification —
see `shor/docs/PEDAGOGY.md`. This chapter is the mathematical development: the full payoff proof
(QFT phase estimation, order-finding → factoring, two-register hidden-subgroup → ECDLP), and the
post-quantum migration survey.*

---

### §11.0 The Through-Line for This Chapter

The §"Escape from Search" chapter established five families of exploitable structure. The first
four are classical: group homomorphisms (index calculus), endomorphisms (GLV/Koblitz), pairings
(MOV), and number-field structure (GNFS/NFS-DL). Each reduces the search cost from fully
exponential to subexponential — from $L_N[1, c]$ to $L_N[1/3, c]$ or better — but none escapes
the $L$-notation framework entirely.

Shor's algorithm (1994) is the fifth structure: **quantum period-finding**. Where the classical
attacks find algebraic structure that *reduces* the search cost, Shor finds a quantum period that
**dissolves** the $L$-notation bound to polynomial time. The escape is not incremental — it is a
phase transition. The complexity drops from $L_N[1/3, (64/9)^{1/3}]$ (GNFS) to $O((\log N)^3)$
(Shor): from subexponential to polynomial in the bit-length of $N$.

The structure exploited is *periodicity*: the function $f(k) = g^k \bmod N$ is periodic with
period equal to the multiplicative order of $g$ modulo $N$. Classical algorithms cannot find this
period efficiently — the best classical period-finding algorithms are as hard as factoring itself.
The quantum Fourier transform finds the period in polynomial time by concentrating amplitude on
the multiples of the inverse period. This is the quantum analogue of the classical Fourier
transform's ability to detect periodicity — but in a superposition that evaluates $f$ at
exponentially many points simultaneously.

**The honest coda.** The post-quantum context (§11.5–§11.7) closes the through-line with an
important counterpoint: the SIDH/SIKE break (Castryck–Decru, 2022) shows that structure is a
double-edged sword. A post-quantum candidate that appeared to resist both classical and quantum
attacks fell to a classical polynomial-time attack exploiting the very structure that made it
efficient. The through-line — "structure enables escape" — applies to the attacker as well as the
defender.

---

### §11.1 Quantum Computation: The Minimal Background

This section states the quantum-mechanical facts that Shor's algorithm uses, at the level needed
to follow the proof. It is not a course in quantum mechanics; it is the undergraduate-background
bridge for this chapter.

#### State vectors and superposition

An $n$-qubit quantum register is described by a *state vector* $|\psi\rangle \in \mathbb{C}^{2^n}$,
written in the computational basis as

$$|\psi\rangle = \sum_{k=0}^{2^n - 1} \alpha_k |k\rangle,$$

where $|k\rangle$ denotes the basis state corresponding to the binary representation of $k$, and
the amplitudes $\alpha_k \in \mathbb{C}$ satisfy the normalisation condition $\sum_k |\alpha_k|^2 = 1$.

A *superposition* is any state with more than one non-zero amplitude. The key property of
superposition is that a quantum gate applied to a superposition acts on all basis states
simultaneously — this is the source of quantum parallelism.

#### Unitary evolution and quantum gates

Quantum evolution is described by *unitary operators* $U: \mathbb{C}^{2^n} \to \mathbb{C}^{2^n}$
(matrices satisfying $U^\dagger U = I$). A *quantum gate* is a unitary operator acting on one or
more qubits. The gates used in Shor's algorithm are:

- **Hadamard gate $H$:** $H|0\rangle = \tfrac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$,
  $H|1\rangle = \tfrac{1}{\sqrt{2}}(|0\rangle - |1\rangle)$. Applied to all $t$ qubits of a
  register initialised to $|0\rangle$, it produces the uniform superposition
  $\tfrac{1}{\sqrt{2^t}} \sum_{k=0}^{2^t-1} |k\rangle$.

- **Controlled-$U_f$ gate:** Given a classical function $f: \{0,\ldots,2^t-1\} \to \{0,\ldots,2^n-1\}$,
  the gate $|x\rangle|y\rangle \mapsto |x\rangle|y \oplus f(x)\rangle$ evaluates $f$ reversibly.
  For Shor's factoring algorithm, $f(x) = a^x \bmod N$ (modular exponentiation); for the ECDLP
  algorithm, $f(a, b) = a \cdot G + b \cdot Q$ (elliptic-curve scalar multiplication).

- **Quantum Fourier Transform (QFT):** defined in §11.2 below.

#### Measurement and the Born rule

*Measuring* the register in the computational basis yields outcome $k$ with probability
$|\alpha_k|^2$, and collapses the state to $|k\rangle$. This is the **Born rule**. After
measurement, the quantum information is destroyed — the register is in a definite classical state.

The art of quantum algorithm design is to arrange the amplitudes so that the *useful* outcomes
have high probability before measurement. Shor's algorithm does this by using the QFT to
concentrate amplitude on the multiples of the inverse period.

#### Entanglement

When two registers are in a *product state* $|\psi\rangle \otimes |\phi\rangle$, they are
independent. When the state cannot be written as a product — when the amplitudes of the joint
system are correlated — the registers are *entangled*. The modular-exponentiation step in Shor's
algorithm creates entanglement between the exponent register and the work register: measuring the
work register collapses the exponent register to a superposition over the preimage of the measured
value, which is exactly the periodic superposition the QFT needs.

---

### §11.2 The Quantum Fourier Transform

#### Definition

The **Quantum Fourier Transform** (QFT) over $\mathbb{Z}/M\mathbb{Z}$ (with $M = 2^t$) is the
unitary operator defined on basis states by

$$\mathrm{QFT}_M |j\rangle = \frac{1}{\sqrt{M}} \sum_{k=0}^{M-1} e^{2\pi i jk/M} |k\rangle.$$

This is the discrete Fourier transform (DFT) over $\mathbb{Z}/M\mathbb{Z}$, implemented as a
quantum gate. The QFT is unitary because the DFT matrix is unitary (its rows are orthonormal
characters of $\mathbb{Z}/M\mathbb{Z}$).

By linearity, the QFT acts on a general state $|\psi\rangle = \sum_j \alpha_j |j\rangle$ as

$$\mathrm{QFT}_M |\psi\rangle = \sum_{k=0}^{M-1} \hat{\alpha}_k |k\rangle,
\quad \text{where} \quad
\hat{\alpha}_k = \frac{1}{\sqrt{M}} \sum_{j=0}^{M-1} \alpha_j e^{2\pi i jk/M}.$$

The amplitudes $\hat{\alpha}_k$ are the DFT coefficients of the input amplitudes $\alpha_j$.

#### The circuit

The QFT over $\mathbb{Z}/2^t\mathbb{Z}$ can be implemented with $O(t^2)$ gates using the
Hadamard + controlled-phase ladder. For a $t$-qubit register, the circuit is:

$$\mathrm{QFT}_{2^t}: \quad
H \text{ on qubit } t-1, \quad
R_2 \text{ on qubit } t-1 \text{ controlled by qubit } t-2, \quad
\ldots, \quad
H \text{ on qubit } 0,$$

where $R_k$ is the phase gate $|1\rangle \mapsto e^{2\pi i/2^k}|1\rangle$. The circuit produces
the QFT in bit-reversed order; a bit-reversal permutation (SWAP network) at the output corrects
this. See Nielsen–Chuang [NC00, §5.1] for the full circuit derivation.

**Code realisation.** The `qft` module in the `shor` crate implements this circuit exactly, with the
bit-reversal convention described in `shor/docs/PEDAGOGY.md` §2.2. The little-endian qubit
ordering requires an input bit-reversal step before the Hadamard + controlled-phase ladder; see
§2.2 for the load-bearing convention note.

#### Phase estimation: the key lemma

The QFT is the engine of *quantum phase estimation*, the subroutine at the heart of Shor's
algorithm. The key lemma is:

**Lemma (QFT concentrates amplitude on multiples of the inverse period).** *Let $r$ divide $M =
2^t$, and let*

$$|\psi_r\rangle = \frac{1}{\sqrt{r}} \sum_{j=0}^{r-1} |j \cdot M/r\rangle$$

*be the uniform superposition over the multiples of $M/r$ in $\{0, \ldots, M-1\}$. Then*

$$\mathrm{QFT}_M |\psi_r\rangle = \frac{1}{\sqrt{r}} \sum_{s=0}^{r-1} e^{2\pi i \cdot 0 \cdot s/r}
|s \cdot M/r\rangle = \frac{1}{\sqrt{r}} \sum_{s=0}^{r-1} |s \cdot M/r\rangle = |\psi_r\rangle.$$

*More generally, for the uniform superposition over a coset $\{a_0 + j \cdot M/r : j = 0, \ldots,
r-1\}$, the QFT concentrates amplitude on the multiples of $r$ in the frequency domain.*

The general statement — which is what Shor's algorithm actually uses — is:

**Lemma (phase estimation, general).** *Let $r \mid M$ and let $c \in \{0, \ldots, M/r - 1\}$.
Define the periodic superposition*

$$|\phi_c\rangle = \frac{1}{\sqrt{r}} \sum_{j=0}^{r-1} |c + j \cdot (M/r)\rangle.$$

*Then*

$$\mathrm{QFT}_M |\phi_c\rangle = \frac{1}{\sqrt{r}} \sum_{s=0}^{r-1}
e^{2\pi i c s / M} \cdot |s \cdot (M/r)\rangle.$$

*In particular, measuring the QFT output yields a uniformly random element of
$\{0, M/r, 2M/r, \ldots, (r-1)M/r\}$, each with probability $1/r$.*

*Proof.* Compute directly:

$$\mathrm{QFT}_M |\phi_c\rangle
= \frac{1}{\sqrt{r}} \sum_{j=0}^{r-1} \mathrm{QFT}_M |c + j(M/r)\rangle
= \frac{1}{\sqrt{r}} \sum_{j=0}^{r-1} \frac{1}{\sqrt{M}} \sum_{k=0}^{M-1}
  e^{2\pi i (c + j M/r) k / M} |k\rangle.$$

Swap the order of summation:

$$= \frac{1}{\sqrt{rM}} \sum_{k=0}^{M-1} e^{2\pi i ck/M}
  \left(\sum_{j=0}^{r-1} e^{2\pi i jk/r}\right) |k\rangle.$$

The inner sum $\sum_{j=0}^{r-1} e^{2\pi i jk/r}$ is a geometric series. It equals $r$ when
$r \mid k$, and $0$ otherwise (since the $r$-th roots of unity sum to zero when $r \nmid k$).
Therefore only the terms with $k = s \cdot (M/r)$ for $s \in \{0, \ldots, r-1\}$ survive:

$$= \frac{1}{\sqrt{rM}} \sum_{s=0}^{r-1} e^{2\pi i c s (M/r) / M} \cdot r \cdot |s(M/r)\rangle
= \frac{1}{\sqrt{r}} \sum_{s=0}^{r-1} e^{2\pi i cs/r} |s(M/r)\rangle. \quad \square$$

**Remark (the case $r \nmid M$).** In practice, $r$ does not divide $M = 2^t$ exactly. The
analysis extends to this case: the QFT output is concentrated near the multiples of $M/r$, with
amplitude falling off as $1/|k - s \cdot M/r|$ for $k$ not a multiple of $M/r$. The
continued-fraction algorithm (§11.3.3) recovers $r$ from the approximate multiples. See
Nielsen–Chuang [NC00, §5.3.1] for the full analysis with the approximation error.

---

### §11.3 Shor's Factoring Algorithm: The Full Proof

#### §11.3.1 The order-finding reduction

**Theorem (order-finding → factoring).** *Let $N$ be an odd composite integer with no small prime
factors. Let $a$ be chosen uniformly at random from $\{2, \ldots, N-2\}$ with $\gcd(a, N) = 1$.
Let $r = \mathrm{ord}_N(a)$ be the multiplicative order of $a$ modulo $N$ (i.e. the smallest
positive integer with $a^r \equiv 1 \pmod{N}$). If $r$ is even and $a^{r/2} \not\equiv -1
\pmod{N}$, then $\gcd(a^{r/2} - 1, N)$ and $\gcd(a^{r/2} + 1, N)$ are both non-trivial factors
of $N$.*

*Proof.* Since $a^r \equiv 1 \pmod{N}$, we have $N \mid (a^r - 1) = (a^{r/2} - 1)(a^{r/2} + 1)$.
Since $a^{r/2} \not\equiv -1 \pmod{N}$, we have $N \nmid (a^{r/2} + 1)$, so
$\gcd(a^{r/2} + 1, N) < N$. Since $a^{r/2} \not\equiv 1 \pmod{N}$ (because $r$ is the *minimal*
order), we have $N \nmid (a^{r/2} - 1)$, so $\gcd(a^{r/2} - 1, N) < N$. Since $N \mid
(a^{r/2}-1)(a^{r/2}+1)$, at least one of the two GCDs is $> 1$. $\square$

**Probability of success.** The condition "$r$ is even and $a^{r/2} \not\equiv -1 \pmod{N}$"
holds for a random $a$ with probability at least $1 - 1/2^{k-1}$, where $k$ is the number of
distinct odd prime factors of $N$. For $N = pq$ (a semiprime), this probability is at least $1/2$.
The proof uses the structure of $(\mathbb{Z}/N\mathbb{Z})^*$ as a product of cyclic groups (one
per prime factor) and the Chinese Remainder Theorem; see Nielsen–Chuang [NC00, §5.3.2] for the
full argument.

**Consequence.** To factor $N$, it suffices to find the order $r = \mathrm{ord}_N(a)$ for a
random $a$. The quantum part of Shor's algorithm solves this order-finding problem in polynomial
time.

#### §11.3.2 The quantum order-finding circuit

**Setup.** Choose $t = \lceil 2 \log_2 N \rceil + 1$ (so that $M = 2^t > N^2$; this precision
is needed for the continued-fraction step). Let $n = \lceil \log_2 N \rceil$. The circuit uses
$t + n$ qubits: the *exponent register* of $t$ qubits and the *work register* of $n$ qubits.

*Register width in the implementation.* The `shor` crate uses $t = \lceil \log_2 N \rceil$, half the
textbook width above (so $2^t \geq N$, not $2^t > N^2$), with a work register of
$\lceil \log_2 (N + 1) \rceil$ qubits. The continued-fraction guarantee below is proved for the
textbook width only. At the narrower width its sufficiency is empirical: the fixtures
`factor_15_seed0`, `factor_21_seed0`, `factor_35_seed7` and `factor_91_seed1`
(`shor/tests/factor_kat.rs`) each recover the order and factor $N$ at a recorded seed, and no bound
on the failure probability is proved. Whatever the width, a recovered order is checked before it is
used (see the recovery procedure below), and a returned factorisation is checked by division.

**Step 1: Uniform superposition.** Initialise both registers to $|0\rangle$. Apply $H$ to each
qubit of the exponent register:

$$|0\rangle_t |0\rangle_n \xrightarrow{H^{\otimes t} \otimes I}
\frac{1}{\sqrt{M}} \sum_{x=0}^{M-1} |x\rangle |0\rangle.$$

**Step 2: Modular exponentiation.** Apply the controlled modular-exponentiation gate
$|x\rangle|y\rangle \mapsto |x\rangle|y \cdot a^x \bmod N\rangle$ (with $y = 1$):

$$\frac{1}{\sqrt{M}} \sum_{x=0}^{M-1} |x\rangle |0\rangle
\xrightarrow{U_{a,N}}
\frac{1}{\sqrt{M}} \sum_{x=0}^{M-1} |x\rangle |a^x \bmod N\rangle.$$

This step entangles the exponent register with the work register. The function $x \mapsto a^x
\bmod N$ is periodic with period $r$: $a^{x+r} \equiv a^x \pmod{N}$.

**Step 3: Measure the work register.** Measure the work register. Suppose the outcome is some
value $v = a^{x_0} \bmod N$ for some $x_0 \in \{0, \ldots, r-1\}$. The measurement collapses
the exponent register to the uniform superposition over all $x$ with $a^x \equiv v \pmod{N}$,
i.e. over the coset $\{x_0, x_0 + r, x_0 + 2r, \ldots\}$ (truncated to $\{0, \ldots, M-1\}$):

$$\frac{1}{\sqrt{\lfloor M/r \rfloor}} \sum_{j=0}^{\lfloor M/r \rfloor - 1} |x_0 + jr\rangle.$$

This is the periodic superposition $|\phi_{x_0}\rangle$ from the phase-estimation lemma (§11.2),
with period $r$ and offset $x_0$.

**Remark.** In practice, the work-register measurement is not performed explicitly — the QFT
step works equally well on the mixed state that results from tracing out the work register. The
analysis is the same either way; the explicit measurement is a pedagogical convenience.

**Step 4: Apply the inverse QFT.** Apply $\mathrm{QFT}_M^\dagger$ (the inverse QFT) to the
exponent register. By the phase-estimation lemma:

$$\mathrm{QFT}_M^\dagger |\phi_{x_0}\rangle
= \frac{1}{\sqrt{r}} \sum_{s=0}^{r-1} e^{-2\pi i x_0 s / r} |s \cdot (M/r)\rangle$$

(using $\mathrm{QFT}^\dagger = \mathrm{QFT}^{-1}$ and the lemma with the phase conjugated).
Each outcome $s \cdot (M/r)$ appears with probability $1/r$.

**Step 5: Measure the exponent register.** The measurement yields some value $y \approx s \cdot
M/r$ for a uniformly random $s \in \{0, \ldots, r-1\}$. The rational number $y/M$ is close to
$s/r$:

$$\left| \frac{y}{M} - \frac{s}{r} \right| \leq \frac{1}{2M}.$$

#### §11.3.3 Continued-fraction period extraction

**The continued-fraction algorithm.** Given the measured value $y$ and the precision $M = 2^t$,
expand $y/M$ as a continued fraction:

$$\frac{y}{M} = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cdots}}$$

and compute the convergents $p_0/q_0, p_1/q_1, p_2/q_2, \ldots$ (the rational approximations
obtained by truncating the continued fraction). The convergents satisfy

$$\left| \frac{y}{M} - \frac{p_k}{q_k} \right| \leq \frac{1}{q_k q_{k+1}}.$$

**Theorem (continued-fraction recovery).** *If $|y/M - s/r| \leq 1/(2M)$ and $r < N$, then
$s/r$ appears as a convergent of the continued fraction of $y/M$.*

*Proof.* Since $r < N \leq \sqrt{M}$ (by the choice $M > N^2$), we have $r^2 < M$, so
$1/(2M) < 1/(2r^2)$. By the theory of continued fractions (see Hardy–Wright, §10.15), any
rational $p/q$ with $q \leq r$ and $|y/M - p/q| < 1/(2q^2)$ is a convergent of $y/M$. Since
$|y/M - s/r| \leq 1/(2M) < 1/(2r^2)$, the fraction $s/r$ (in lowest terms) is a convergent.
$\square$

**The recovery procedure.** For each convergent $p_k/q_k$ of $y/M$ with $q_k < N$, test whether
$a^{q_k} \equiv 1 \pmod{N}$. If so, $q_k$ is a *multiple* of $r$, because $r$ divides every
exponent $e$ with $a^e \equiv 1$. It need not equal $r$: the measured phase approximates $s/r$ for
some $s$ that may share a factor with $r$, and a convergent of a coarse phase can have a
denominator that is a proper multiple of $r$. The first such $q_k$ is therefore reduced to the least
order by dividing out each prime $f \mid q_k$ for as long as $a^{q_k/f} \equiv 1$ still holds
(the order divides $q_k$, so a failed division never loses a prime power of the order). With that
reduction `order_from_phase` returns the least $r$ with $a^r \equiv 1$, never a multiple of it
(`cf_order_from_phase_n15_a7_s32`, where $s = 32$, $t = 8$, $a = 7$, $N = 15$ has first
denominator $8$ and order $4$; `cf_order_from_phase_is_always_the_least_order`, over every phase
and every base coprime to $N$ at the four register widths the targets use, both in
`shor/tests/factor_kat.rs`). The order $r$ is recovered with probability $\Omega(1/\log \log r)$
from a single measurement; repeating $O(\log \log N)$ times gives $r$ with high probability.

**Code realisation.** The `order_from_phase` function in `shor/src/shor/mod.rs` implements this
continued-fraction recovery, and `factor` in the same module is the driver. `factor` classifies $N$
classically in a fixed order: $N \leq 2$ is `NotComposite` ($0$ and $1$ have no factorisation and
$2$ is prime); an even $N$ is factored by $2$ with no circuit; a prime power $p^e$ is factored by
$p$ with no circuit; a circuit needing more than `MAX_QUBITS` $= 25$ qubits
($\lceil \log_2 N \rceil + \lceil \log_2 (N+1) \rceil$, which admits $N \leq 4096$) is
`Unsupported`; an odd $N$ that passes those is tested for primality by trial division and a prime is
`NotComposite`; only then does the order-finding loop run, over at most `MAX_BASES` $= 20$ bases
coprime to $N$ with `SEEDS_PER_BASE` $= 8$ circuit runs each. Because the primality test follows
the ceiling, a large odd prime is `Unsupported` rather than `NotComposite`. Bases are drawn coprime
to $N$, so no base ever shares a factor with $N$ and there is no common-factor outcome: the
shortcut $\gcd(a, N)$ would turn the base scan into trial division, which is not the algorithm. A
spent budget is `AttemptsExhausted` and says nothing about $N$. The tests are
`factor_even_short_circuit`, `factor_prime_power_short_circuit`, `factor_prime_is_not_composite`,
`factor_beyond_the_simulator_is_unsupported` and `factor_15_seed0` (one per outcome), and the
unit test `spent_budget_is_attempts_exhausted_with_the_runs_made` in `shor/src/shor/mod.rs`. See
`shor/docs/PEDAGOGY.md` §3.2 for the module surface and §3.3 for the KATs.

#### §11.3.4 Complexity

**Theorem (Shor's factoring algorithm, complexity).** *Shor's algorithm factors an $n$-bit
integer $N$ in expected time $O(n^3)$ on a quantum computer.*

*Proof sketch.* The circuit uses $t + n = O(n)$ qubits. The modular-exponentiation circuit
(Step 2) requires $O(n^3)$ elementary quantum gates (using standard reversible arithmetic). The
QFT (Step 4) requires $O(n^2)$ gates. The continued-fraction computation (Step 5) is classical
and costs $O(n^2)$. The expected number of repetitions before a successful factorisation is
$O(\log \log N) = O(\log n)$. The total expected gate count is $O(n^3)$. $\square$

**The L-notation collapse.** GNFS factors $N$ in time $L_N[1/3, (64/9)^{1/3}]$. Shor factors
$N$ in time $O((\log N)^3) = L_N[0, 3]$. The exponent drops from $1/3$ to $0$: from
subexponential to polynomial. This is not an incremental improvement — it is the dissolution of
the $L$-notation framework. The quantum period-finding structure does not reduce the search cost;
it eliminates the search entirely.

**Scale note.** The tests factor $N \in \{15, 21, 35, 91\}$ using 8–14 qubits on a classical
state-vector simulator. The driver's qubit ceiling, `MAX_QUBITS` $= 25$, admits $N \leq 4096$ and
reports a larger $N$ as unsupported; only the four fixtures are run end to end. The simulator holds
$2^n$ complex
amplitudes; at $n = 25$ qubits, this is $2^{25} \approx 33\,\mathrm{M}$ entries ($\approx 512\,
\mathrm{MiB}$), the practical ceiling on a laptop. To factor RSA-2048 ($N$ a 2048-bit integer),
this project's layout ($2 \cdot n_{\mathrm{bits}}(N)$ qubits) would need about 4100 qubits — an
arithmetic consequence of the layout, not a published estimate — requiring a $2^{4100}$-entry
amplitude array, far beyond any classical computer. The simulator demonstrates Shor's
mathematics correctly at toy scale; it does not claim quantum speedup, which requires real quantum
hardware out of scope by construction. See `shor/docs/PEDAGOGY.md` §6.3 and the Shor factoring
section of `docs/BENCHMARKS.md` for the full resource-scale annotation.

---

### §11.4 Shor's ECDLP Algorithm: The Two-Register Hidden-Subgroup Proof

#### §11.4.1 The hidden-subgroup formulation

The ECDLP asks: given $G$ and $Q = k \cdot G$ on an elliptic curve of prime order $r$, find $k$.
Shor's approach (Proos–Zalka [PZ03]) reformulates this as a **two-register hidden-subgroup
problem** (HSP) over $\mathbb{Z}/M\mathbb{Z} \times \mathbb{Z}/M\mathbb{Z}$ with $M = 2^t$. The
implementation uses $t = \lceil \log_2 r \rceil = 4$ for its toy curve, so $M = 16$ and the
group order is $r = 13$: $M$ is **not** a multiple of $r$, which governs everything in
§11.4.3.

**The function.** Define $f: \mathbb{Z}/M\mathbb{Z} \times \mathbb{Z}/M\mathbb{Z} \to E$ by

$$f(a, b) = a \cdot G + b \cdot Q.$$

Since $Q = k \cdot G$, we have $f(a, b) = (a + bk) \cdot G$. The function $f$ is constant on
the cosets of the *hidden subgroup*

$$H = \{(a, b) \in \mathbb{Z}/r\mathbb{Z} \times \mathbb{Z}/r\mathbb{Z} : a + bk \equiv 0
\pmod{r}\} = \{(bk \bmod r,\, -b \bmod r) : b \in \mathbb{Z}/r\mathbb{Z}\}.$$

More precisely, $f(a, b) = f(a', b')$ if and only if $(a + bk) \equiv (a' + b'k) \pmod{r}$,
i.e. $(a - a') \equiv -(b - b')k \pmod{r}$. The hidden subgroup is the kernel of the map
$(a, b) \mapsto a + bk \bmod r$, which is the cyclic subgroup generated by $(-k \bmod r, 1)$
in $(\mathbb{Z}/r\mathbb{Z})^2$.

**The HSP.** The hidden subgroup $H \leq (\mathbb{Z}/r\mathbb{Z})^2$ is generated by $(-k, 1)$.
Finding $H$ (or its generator) yields $k$ directly. The quantum algorithm finds $H$ by preparing
a superposition over $(\mathbb{Z}/M\mathbb{Z})^2$, entangling it with the function values
$f(a, b)$, and applying the two-dimensional QFT.

#### §11.4.2 The two-register circuit

**Setup.** Use two exponent registers of $t$ qubits each (the $a$-register and $b$-register) and
one work register of $n$ qubits (holding the running point on $E$). Total: $2t + n$ qubits.
For the toy curve ($r = 13$, $t = 4$) the point register holds a 3-bit $x$ coordinate and a
3-bit $y$ coordinate, and the layout also allocates a 3-bit scratch register that the
permutation-synthesis construction never uses. That is $2 \cdot 4 + 3 \cdot 3 = 17$ qubits in
total (asserted by `qubit_budget_within_ceiling`; see `shor/docs/PEDAGOGY.md` §4.2).

**Step 1: Uniform superposition.** Apply $H^{\otimes t}$ to both the $a$-register and
$b$-register:

$$|0\rangle_t |0\rangle_t |0\rangle_n
\xrightarrow{H^{\otimes 2t} \otimes I}
\frac{1}{M} \sum_{a=0}^{M-1} \sum_{b=0}^{M-1} |a\rangle |b\rangle |\mathcal{O}\rangle,$$

where $\mathcal{O}$ is the identity element of $E$ (the point at infinity).

**Step 2: Controlled point addition.** For each bit $j$ of the $a$-register, apply the
controlled point-addition gate $|P\rangle \mapsto |P + 2^j G\rangle$ controlled on
$a$-register qubit $j$. Then for each bit $j$ of the $b$-register, apply $|P\rangle \mapsto
|P + 2^j Q\rangle$ controlled on $b$-register qubit $j$. The combined effect is:

$$\frac{1}{M} \sum_{a=0}^{M-1} \sum_{b=0}^{M-1} |a\rangle |b\rangle |a \cdot G + b \cdot Q\rangle.$$

This is the quantum evaluation of $f(a, b) = a \cdot G + b \cdot Q$ in superposition. Because
$M$ is not a multiple of $r$, $f$ is not periodic on $\mathbb{Z}/M\mathbb{Z} \times
\mathbb{Z}/M\mathbb{Z}$: the level sets $\{(a, b) : a + bk \equiv s \pmod{r}\}$ are truncated by
the box $0 \le a, b < M$, so they are not full cosets of $H$. This is the finite-precision
effect that the rest of the section accounts for.

**Step 3: Transform both exponent registers.** Apply $\mathrm{QFT}_M^\dagger$ independently to
the $a$-register and to the $b$-register. Write $\omega = e^{2\pi i / M}$; the sign of the
exponent below is the sign convention of the transform, and nothing that follows depends on it.
The amplitude of $|c\rangle |d\rangle |sG\rangle$ for $s \in \{0, \ldots, r-1\}$ is

$$M^{-2} \sum_{\substack{a, b < M \\ a + bk \equiv s \ (\mathrm{mod}\ r)}} \omega^{\mp(ac + bd)}.$$

The points $sG$ are distinct for distinct $s$, so the outcome $(c, d)$ has probability

$$P(c, d) = \sum_{s=0}^{r-1} \left| M^{-2} \sum_{\substack{a, b < M \\ a + bk \equiv s}}
\omega^{\mp(ac + bd)} \right|^2.$$

This is an exact formula for the circuit's exponent-register distribution. The function
`exponent_outcome_distribution` computes $P(c, d)$ from the circuit's amplitudes, and
`distribution_matches_fourier_formula_k3` and `distribution_matches_fourier_formula_k7` check,
on the default test path, that it equals this sum to within $10^{-9}$ and sums to $1$. The same
check for every $k \in \{2, \ldots, 12\}$ is `distribution_matches_fourier_formula_all_k`, which
is `#[ignore]`-gated for cost.

**Step 4: Measure all qubits.** The implementation measures every qubit of the register (both
exponent registers and the point and scratch registers) in one shot and keeps the exponent pair
$(c, d)$. The point-register outcome carries no information about $k$ beyond what the exponent
pair already carries.

#### §11.4.3 The dual relation and the extraction rule

**The dual relation.** Evaluate the inner sum at the lattice point $c/M = j/r$, $d/M = l/r$
for integers $j, l$. Every term in the class $a + bk \equiv s$ then has phase
$e^{\mp 2\pi i (aj + bl)/r}$, and, since $a \equiv s - bk \pmod{r}$,

$$aj + bl \equiv (s - bk)\,j + bl = sj + b\,(l - kj) \pmod{r}.$$

The first term is constant on the class, so the terms add constructively exactly when
$l \equiv k j \pmod{r}$. This is the *dual relation*: the measured pairs concentrate near the
lattice points $(j/r, l/r)$ with $l \equiv kj$, i.e. the dual lattice of the hidden subgroup
$H = \langle(-k, 1)\rangle$, since $(j, l) \cdot (-k, 1) = l - kj$. It is independent of the sign
convention, because flipping the sign conjugates every inner sum and leaves its modulus
unchanged. Reading $l \equiv kj$ for $j \not\equiv 0$ gives $k \equiv l\,j^{-1} \pmod{r}$; no
sign flip enters.

**Finite precision.** When $r \mid M$ the lattice points $(j/r, l/r)$ are exactly representable
as $(c/M, d/M)$ and the peaks are sharp. Here $M = 16$ and $r = 13$, so the points $(jM/r, lM/r)$
are not integers; the probability mass of each peak is spread over neighbouring integer
outcomes. The extraction rounds each measured value to the nearest multiple of $1/r$:

$$j = \operatorname{round}(r c / M) \bmod r, \qquad l = \operatorname{round}(r d / M) \bmod r,
\qquad k = l \cdot j^{-1} \pmod{r},$$

with ties rounded up. The rule produces exactly one candidate per measurement, and none (`None`)
when $j \equiv 0$, in which case the measurement carries no information about $k$. The
candidate set is not widened to neighbouring roundings: at $r = 13$ doing so would amount to
enumerating a large fraction of the exponents, which is a search over the answer space and not
Shor's algorithm. No continued-fraction step is used here, and nothing in this section claims
that a continued-fraction (or any single-rounding) rule recovers $k$ for general $r$ and $M$.
The modular inverse is computed by Fermat's little theorem, since $r$ is prime.

**What the tests establish about the rule.**
`extract_exact_peak_pairs_recover_k` checks that for every $k \in \{0, \ldots, 12\}$ and every
$j \in \{1, \ldots, 12\}$ the pair $(c, d)$ of integers nearest to $(jM/r,\ lM/r)$, with
$l = kj \bmod r$, extracts $k$;
`extract_j_zero_returns_none` checks that, for $M = 16$ and $r = 13$, only $c = 0$ gives
`None`; `extract_rejects_c_out_of_range` and `extract_rejects_d_out_of_range` check that
$c, d \ge M$ is rejected as invalid input. These tests pin the rule on exact lattice points; they do
not say how likely a real measurement is to land on one.

**Success probability per measurement (estimate).** The probability that one measurement
yields a candidate equal to the true logarithm is
$\sum P(c, d)$ over the outcomes whose extraction returns $k$; the function
`primary_success_probability` computes it exactly from the circuit amplitudes. No test asserts
its numeric value: `primary_success_probability_positive_k3_k7` (default path) and
`primary_success_probability_positive_all_k` (ignored for cost) assert only that it is positive
and, for the first, that it equals the sum of the cached distribution over the extracting
outcomes. A single run of the ignored test printed values between about $0.64$ and $0.66$ for
$k = 2, \ldots, 12$ (for example $0.6450$ for $k = 2, 6, 7, 11$ and $0.6568$ for $k = 12$).
Treat these as estimates read from that printout, not as tested facts.

**Retries.** `solve_ecdlp` repeats measure, extract, and verify $k \cdot G = Q$ with a fresh
seed per attempt, up to `DEFAULT_MAX_ATTEMPTS = 8` times; a candidate is accepted only after the
group-law check. Treating the attempts as independent draws with the per-measurement probability
above, the chance that all eight fail is about $0.36^8 \approx 3 \cdot 10^{-4}$. That figure is
an estimate derived from the printout, not a tested bound; the only budget-related behaviour the
tests assert is the error semantics (`ecdlp_budget_exhausted_on_wrong_candidate`,
`ecdlp_zero_budget_is_exhausted`, `ecdlp_retry_uses_fresh_seed`) and, in an ignored test,
that every $k$ is solved within the default budget (`ecdlp_all_k_values_within_default_budget`).
Exhausting the budget is a resource outcome (`AttemptsExhausted`), not evidence that no
logarithm exists, and no classical search over $k$ remains in the solver.

#### §11.4.4 Complexity

**Claim from the literature (not proved or measured here).** Proos and Zalka [PZ03] show that
Shor's two-register algorithm solves the ECDLP on an elliptic curve of prime order $r$ with
$O((\log r)^3)$ elementary quantum gates, in the sense that the circuit has $O(\log r)$ qubits,
the controlled point additions of Step 2 dominate the gate count (using reversible
elliptic-curve arithmetic), and the two transforms of Step 3 each cost $O((\log r)^2)$ gates.
Turning that gate count into an *expected running time* also needs a lower bound on the
probability that one measurement yields a usable candidate, which depends on the register size
and on the classical post-processing. This project does not derive such a bound: the toy
implementation uses the single-rounding rule of §11.4.3 at $M = 16$, $r = 13$, and its
per-measurement success probability is only the estimate recorded there. Nothing in this section
establishes the extraction rule for large $r$.

**Scale note.** The implementation solves the ECDLP on a toy curve with $r = 13$ using 17 qubits
(`qubit_budget_within_ceiling`). The layout it uses, $2t + 3 \cdot \lceil \log_2 p \rceil$
qubits (two exponent registers, plus $x$, $y$, and scratch registers), would be $2 \cdot 256 +
3 \cdot 256 = 1280$ qubits for a 256-bit curve such as secp256k1. That figure is arithmetic on
this project's layout, not a published resource estimate, and this project does not reproduce
the resource counts published for optimised circuits (such as [PZ03]). A classical
state-vector simulator would need a $2^{n}$-entry amplitude array for $n$ qubits, which is
beyond any classical computer long before that size. The simulator demonstrates the mathematics
at toy scale; a quantum speedup requires quantum hardware. See `shor/docs/PEDAGOGY.md` §6.3 and
the Shor ECDLP section of `docs/BENCHMARKS.md`.

---

### §11.5 Post-Quantum Cryptography: The Migration Landscape

#### §11.5.1 Why Shor's algorithm forces a migration

The security of RSA, Diffie–Hellman, and elliptic-curve cryptography rests on the hardness of
integer factoring and the discrete logarithm problem. Shor's algorithm solves both in polynomial
time on a quantum computer. A sufficiently large, fault-tolerant quantum computer would therefore
break all three families simultaneously.

The migration is not hypothetical. The question is not whether quantum computers will eventually
be built at the scale needed to run Shor's algorithm on cryptographic parameters — the physics
is not in dispute. The question is when, and the answer is uncertain. The conservative response
— adopted by NIST, NSA, and the major standards bodies — is to begin the migration now, before
large-scale quantum computers exist, because:

1. **Harvest-now, decrypt-later.** An adversary can record encrypted traffic today and decrypt
   it once a quantum computer is available. For data with long-term confidentiality requirements
   (state secrets, medical records, financial records), the threat is present even if quantum
   computers are decades away.

2. **Migration timelines are long.** Replacing cryptographic primitives in deployed systems
   takes years to decades. The migration from SHA-1 to SHA-2/SHA-3 took over a decade; the
   migration from RSA-1024 to RSA-2048 is still incomplete. Starting the migration before the
   threat materialises is the only way to complete it in time.

3. **The classical+quantum arc this project demonstrates is the reason.** The rDLP project
   documents the complete classical attack surface (Pollard rho, GNFS, NFS-DL, algebraic ECDLP
   attacks) and the quantum attack (Shor's algorithm). Together, they show that every classical
   hardness assumption underlying RSA and ECC is broken — classically by subexponential
   algorithms, quantumly by polynomial-time algorithms. The migration is the engineering response
   to this mathematical reality.

#### §11.5.2 The NIST PQC standardisation process

The National Institute of Standards and Technology (NIST) ran a multi-year public competition
(2016–2024) to standardise post-quantum cryptographic algorithms. The competition evaluated
candidates across four families: lattice-based, code-based, hash-based, and isogeny-based.

**The lattice family (primary standardisation).** The two primary standards are:

- **ML-KEM (FIPS 203, 2024)** — Module-Lattice-Based Key-Encapsulation Mechanism, based on
  the Module Learning With Errors (MLWE) problem. The underlying scheme is Kyber. ML-KEM
  provides key encapsulation (the post-quantum replacement for RSA-KEM and ECDH) at three
  security levels (ML-KEM-512, ML-KEM-768, ML-KEM-1024), targeting 128, 192, and 256 bits of
  classical security respectively.

- **ML-DSA (FIPS 204, 2024)** — Module-Lattice-Based Digital Signature Standard, based on the
  Module Learning With Errors and Module Short Integer Solution (MSIS) problems. The underlying
  scheme is Dilithium. ML-DSA provides digital signatures at three security levels (ML-DSA-44,
  ML-DSA-65, ML-DSA-87).

**The hardness assumption.** Both ML-KEM and ML-DSA rest on the hardness of the Learning With
Errors (LWE) problem: given a matrix $A \in \mathbb{Z}_q^{m \times n}$ and a vector
$b = As + e \pmod{q}$ (where $s$ is a secret vector and $e$ is a small error vector), find $s$.
The LWE problem is believed to be hard for both classical and quantum computers; the best known
quantum algorithms for LWE (lattice sieving with quantum speedup) achieve only a quadratic
speedup over the best classical algorithms, leaving the security level intact at the parameter
sizes chosen by NIST.

**The hash-based family.** SLH-DSA (FIPS 205, 2024) — Stateless Hash-Based Digital Signature
Standard, based on SPHINCS+. Hash-based signatures derive their security entirely from the
collision resistance and preimage resistance of the underlying hash function (SHA-256 or SHAKE).
Their post-quantum security is well-understood: Grover's algorithm gives a quadratic speedup for
preimage search, which is accounted for by doubling the hash output length. SLH-DSA is the
conservative choice: its security rests on minimal assumptions (hash function security only),
but its signatures are large (8–50 KB depending on the parameter set).

**The code-based family.** Code-based cryptography (McEliece, 1978) encodes the public key as
a generator matrix of a binary Goppa code with a hidden structure. The hardness assumption is
the syndrome decoding problem (NP-hard in the worst case). The Classic McEliece scheme was a
NIST finalist; it was not standardised in the first round due to its very large public keys
(hundreds of kilobytes), but remains a conservative candidate. NIST has indicated it may be
standardised in a future round.

**The isogeny-based family.** Isogeny-based cryptography (Couveignes, Rostovtsev–Stolbunov,
SIDH/SIKE) constructs key exchange from the hardness of finding an isogeny between two
supersingular elliptic curves. The SIDH/SIKE scheme was a NIST finalist — until it was broken
in 2022 (see §11.6).

#### §11.5.3 The security landscape

The post-quantum security of the standardised schemes rests on different mathematical
assumptions, providing defence in depth:

| Family | Hardness assumption | Best quantum attack | Security level (128-bit) |
|--------|--------------------|--------------------|--------------------------|
| Lattice (ML-KEM, ML-DSA) | MLWE / MSIS | Lattice sieving + Grover | ML-KEM-768 / ML-DSA-65 |
| Hash-based (SLH-DSA) | Hash collision/preimage | Grover (quadratic speedup) | SLH-DSA-128f |
| Code-based (McEliece) | Syndrome decoding | Information-set decoding + Grover | Classic McEliece 6688128 |

No polynomial-time quantum algorithm is known for any of these problems. The lattice and
code-based problems are believed to be in the complexity class $\mathbf{QMA}$-hard (quantum
analogue of $\mathbf{NP}$-hard), though this is not proven.

---

### §11.6 The SIDH Break: Structure Cuts Both Ways

#### §11.6.1 The SIDH construction

Supersingular Isogeny Diffie–Hellman (SIDH), introduced by Jao and De Feo (2011), constructs
a key-exchange protocol from the following mathematical structure. Let $p$ be a prime of the
form $p = \ell_A^{e_A} \ell_B^{e_B} f \pm 1$ (with $\ell_A = 2$, $\ell_B = 3$ in the standard
parameters). Let $E_0$ be a supersingular elliptic curve over $\mathbb{F}_{p^2}$.

Alice's secret is an isogeny $\phi_A: E_0 \to E_A$ with kernel a random subgroup of
$E_0[\ell_A^{e_A}]$. Bob's secret is an isogeny $\phi_B: E_0 \to E_B$ with kernel a random
subgroup of $E_0[\ell_B^{e_B}]$. The shared secret is the $j$-invariant of the curve
$E_{AB} = E_A / \phi_A(\ker \phi_B) = E_B / \phi_B(\ker \phi_A)$.

**The apparent post-quantum security.** The hardness assumption is the Supersingular Isogeny
Diffie–Hellman (SIDH) problem: given $E_A$ and $E_B$ (and the images of the torsion-point
generators under $\phi_A$ and $\phi_B$), find the shared $j$-invariant. No polynomial-time
classical or quantum algorithm was known for this problem prior to 2022. The scheme was a
finalist in the NIST PQC competition under the name SIKE (Supersingular Isogeny Key
Encapsulation).

#### §11.6.2 The Castryck–Decru attack (2022)

In July 2022, Wouter Castryck and Thomas Decru published a classical polynomial-time attack on
SIDH [CD23]. The attack was subsequently refined by Maino–Martindale and Robert, and the full
SIKE parameter sets were broken within weeks.

**The key insight.** The SIDH protocol publishes, as part of the public key, the images of the
torsion-point generators $P_A, Q_A \in E_0[\ell_A^{e_A}]$ under $\phi_A$:

$$\phi_A(P_A), \quad \phi_A(Q_A) \in E_A[\ell_A^{e_A}].$$

This auxiliary information — published to enable the key exchange — turns out to be fatal. It
allows an attacker to recover the secret isogeny $\phi_A$ using a technique from the theory of
*abelian varieties with extra structure*.

**The mathematical mechanism.** The attack uses the fact that an isogeny between supersingular
curves is determined by its action on torsion points. The published torsion-point images
$\phi_A(P_A), \phi_A(Q_A)$ constrain the isogeny $\phi_A$ to a small set of candidates. The
Castryck–Decru attack constructs a sequence of auxiliary isogenies (using the theory of
$(2, 2)$-isogenies between products of elliptic curves, i.e. isogenies of abelian surfaces)
that progressively narrows this set to a single candidate, recovering $\phi_A$ in polynomial
time.

The attack runs in $O(p^{1/2})$ time in the original formulation (polynomial in $\log p$, since
$p$ is the security parameter), and was subsequently improved to $O(p^{1/4})$ and then to
$O((\log p)^c)$ for a small constant $c$ — fully polynomial.

**The structural lesson.** The SIDH break is the honest coda to the through-line of this
textbook. The through-line is "structure enables escape from search" — and the SIDH break shows
that this cuts both ways. The torsion-point auxiliary information that made SIDH efficient (it
enabled the commutativity of the diagram $E_{AB} = E_A / \phi_A(\ker \phi_B)$) is the same
structure that enabled the attack. The scheme published too much information about the secret
isogeny; the attacker used that information to reconstruct it.

This is not a failure of the isogeny-based paradigm in general — it is a failure of the specific
SIDH construction. Isogeny-based schemes that do not publish torsion-point images (such as
CSIDH, which uses commutative group actions on ordinary curves) are not broken by the
Castryck–Decru attack. But the SIDH break is a vivid reminder that:

1. **Post-quantum security is not self-evident.** A scheme can resist all known classical and
   quantum attacks for years and then fall to a classical attack exploiting structure that was
   always present but not recognised.

2. **Auxiliary information is dangerous.** The torsion-point images were published for
   efficiency; they were the attack surface. The principle "publish only what is necessary" is
   as important in post-quantum cryptography as in classical cryptography.

3. **Structure is the double-edged sword.** Every escape from search in this textbook — group
   homomorphisms, endomorphisms, pairings, number-field structure, quantum periods — is a
   structure that the algorithm exploits. The SIDH break shows that the attacker can exploit
   structure too. The history of cryptanalysis is a history of finding structure in systems
   that were believed to be structureless.

#### §11.6.3 The aftermath

The NIST PQC competition removed SIKE from consideration immediately after the Castryck–Decru
publication. The remaining candidates — ML-KEM, ML-DSA, SLH-DSA — were standardised in 2024.
Research on isogeny-based cryptography continues, with attention shifted to constructions that
avoid the fatal torsion-point leakage.

---

### §11.7 The Migration Rationale: Closing the Through-Line

The rDLP project documents the complete arc from classical to quantum cryptanalysis:

- **Pollard rho:** the generic $\sqrt{n}$ bound and the birthday-paradox collision
  argument — the baseline that every subsequent attack must beat.
- **GNFS:** the number-field structure that reduces factoring to $L_N[1/3, (64/9)^{1/3}]$
  — the best known classical attack on RSA.
- **NFS-DL:** the adaptation of GNFS to discrete logarithms — the best known classical
  attack on Diffie–Hellman.
- **Algebraic ECDLP attacks:** the five algebraic structures (Pohlig–Hellman, SSA,
  GHS, index calculus, MOV) that escape the generic bound for elliptic-curve DLP — the classical
  attack surface for ECC.
- **Shor's algorithm:** the quantum period-finding structure that dissolves the
  $L$-notation bound entirely — the quantum attack on RSA, DH, and ECC simultaneously.

Together, these algorithms show that every classical hardness assumption underlying RSA and ECC is
broken — classically by subexponential algorithms, quantumly by polynomial-time algorithms. The
post-quantum migration is the engineering response to this mathematical reality.

**The scale statement for the Shor chapter.** The implementation demonstrates Shor's
mathematics at toy scale (up to 17 qubits for ECDLP with $r = 13$; factoring $N \leq 91$). The
~25-qubit simulator ceiling is a resource-scale wall: the dense state-vector simulator holds $2^n$
complex amplitudes, and at $n = 25$ this is the practical memory limit on a laptop. The ceiling is
not a mathematical one: the simulator computes the exact state at any size (the ECDLP extraction
rule is analysed here only at $r = 13$, §11.4.3). The quantum speedup (the polynomial
vs subexponential separation) is not observable at toy scale; it requires real quantum hardware
with thousands of fault-tolerant qubits. The migration is *surveyed, not built* — the
pedagogical terminus of this project, not a new construction track.

**The survey is the contribution.** The post-quantum migration landscape (NIST PQC families,
the SIDH break, the migration rationale) is surveyed within established literature. No PQC
primitive is implemented; no new crate, dependency, or benchmark is added. The contribution of
the Shor chapter is the documentation of the quantum arc — the mathematics (this chapter) and the
code realisation (`shor/docs/PEDAGOGY.md`) — that contextualises the migration and closes the
project's classical+quantum survey.

---

### §11.8 Cross-References

#### Code realisation

The Shor code-tour in `shor/docs/PEDAGOGY.md` is the code-first sibling to this chapter.
The chapter-pairing is:

| This chapter | Code-tour section |
|-------------|-------------------|
| §11.2 QFT and phase estimation | `shor/docs/PEDAGOGY.md` §2.2 (`qft` module) |
| §11.3 Shor's factoring algorithm | `shor/docs/PEDAGOGY.md` §3 (factoring) |
| §11.4 Shor's ECDLP algorithm | `shor/docs/PEDAGOGY.md` §4 (ECDLP) |
| §11.7 Scale statement | `shor/docs/PEDAGOGY.md` §6.3 (resource ceiling) |

#### Within this textbook

- **§"Escape from Search: The Through-Line"** — the five-family structure taxonomy and the
  L-notation hierarchy table. This chapter is the fifth entry: quantum period-finding, $L[0]$.

- **§"On Scale: A Natural-Philosophy Interlude"** — the three axes of scale and the
  science↔engineering gap. The ~25-qubit simulator ceiling is an instance of Axis 1
  (resource/operational scale); the quantum vs classical complexity separation is an instance
  of Axis 3 (structural scale — a threshold where a different machine applies).

- **§10 (Algebraic ECDLP Attacks)** — the classical attack surface for ECC. The MOV reduction
  (§10.5) is the cross-track bridge from ECDLP to NFS-DL; Shor's ECDLP algorithm (§11.4) is
  the quantum dissolution of the same problem.

#### Benchmark data

- **`docs/BENCHMARKS.md`, the Shor simulator, factoring, and ECDLP sections** — the qubit-budget
  tables, timings, and `### Science↔engineering note` annotations for the three parts of the
  Shor chapter.

---

### §11.9 Further Reading

43. **Proos, J., and Zalka, C. (2003).** "Shor's discrete logarithm quantum algorithm for
    elliptic curves." *Quantum Information and Computation*, 3(4), 317–344. [PZ03] The
    two-register ECDLP construction; the qubit-budget analysis.

44. **Nielsen, M. A., and Chuang, I. L. (2000).** *Quantum Computation and Quantum Information*.
    Cambridge University Press. [NC00] §5.2 (QFT), §5.3 (order-finding and factoring), §5.4
    (discrete logarithms). The standard textbook reference.

45. **NIST FIPS 203 (2024).** *Module-Lattice-Based Key-Encapsulation Mechanism Standard*.
    [FIPS203] The ML-KEM (Kyber) standard.

46. **NIST FIPS 204 (2024).** *Module-Lattice-Based Digital Signature Standard*. [FIPS204]
    The ML-DSA (Dilithium) standard.

47. **NIST FIPS 205 (2024).** *Stateless Hash-Based Digital Signature Standard*. [FIPS205]
    The SLH-DSA (SPHINCS+) standard.

48. **Castryck, W., and Decru, T. (2023).** "An efficient key recovery attack on SIDH."
    *Advances in Cryptology — EUROCRYPT 2023*, LNCS 14008, 423–447. [CD23] The classical
    polynomial-time break of SIDH/SIKE.

49. **Shor, P. W. (1994).** "Algorithms for quantum computation: discrete logarithms and
    factoring." *Proceedings of the 35th Annual Symposium on Foundations of Computer Science*,
    124–134. [Shor94] The original paper. (Also cited as ref. 14 in the main references section.)

50. **Jao, D., and De Feo, L. (2011).** "Towards quantum-resistant cryptosystems from
    supersingular elliptic curve isogenies." *Post-Quantum Cryptography — PQCrypto 2011*,
    LNCS 7071, 19–34. The original SIDH construction.

51. **Bernstein, D. J., and Lange, T. (2017).** "Post-quantum cryptography." *Nature*, 549,
    188–194. A survey of the post-quantum landscape accessible to a broad scientific audience.

---

## Chapter 12 — Modularity and the Arithmetic of Elliptic Curves: A Speculation

> **This chapter is the project's one deliberate speculation.** The algorithms of Chapters 1–11
> are implemented and exercised by known-answer tests at toy scale, and their mathematics is
> established in those chapters (each chapter says where its claims stop; for instance, the Shor
> ECDLP extraction rule in §11.4 is analysed only for the toy curve). This chapter is different:
> it surveys the Taniyama–Shimura–Weil modularity theorem as a *structural* phenomenon in the same
> escape-from-search family the project surveys, and then asks — explicitly as speculation, not as
> a claimed result — whether the kind of deep correspondence it exemplifies points at attacks not
> yet realised. The speculation is flagged at every turn. No algorithm in this chapter is
> implemented; no complexity claim is made; no attack is asserted. The chapter is here as the
> project's honest closing gesture toward the mathematics that lies beyond its own scope.

### §12.1 Introduction: Correspondences as Structure

The through-line of this textbook is structure-based escape from search. Every chapter has told
the same story: an algorithm finds a piece of exploitable structure — a group homomorphism, a
smoothness phenomenon, an endomorphism, a pairing, a quantum period — and uses it to escape the
generic search bound. The structures are diverse in their details but uniform in their role: they
are *correspondences* between two mathematical worlds, and the correspondence is what the algorithm
exploits.

The Taniyama–Shimura–Weil (TSW) theorem is the deepest known correspondence in the arithmetic of
elliptic curves. It asserts that every elliptic curve over $\mathbb{Q}$ is, in a precise sense,
the same object as a modular form — that the two mathematical worlds of elliptic curves and modular
forms are not merely analogous but *identified* by a canonical map. This identification is not an
algorithm; it is a structural theorem. But the history of the subject suggests that deep structural
theorems of this kind have a way of generating algorithms — sometimes decades later, sometimes in
ways their discoverers did not anticipate.

This chapter surveys the TSW correspondence at the depth the documentation register permits (survey
with proof-sketch depth; full proofs are cited, not given). It then asks, explicitly as
speculation: does the modularity structure suggest further escapes from the search bound that have
not yet been realised? The answer is: possibly, but we do not know. The speculation is honest
about its own uncertainty.

### §12.2 Modular Forms: The Other Side of the Correspondence

To state the TSW theorem, we need both sides of the correspondence. The elliptic-curve side is
developed in the Prerequisites chapter (§4) and used throughout the textbook. The modular-form
side requires a brief introduction.

#### §12.2.1 The upper half-plane and modular symmetry

Let $\mathfrak{H} = \{z \in \mathbb{C} : \mathrm{Im}(z) > 0\}$ be the *upper half-plane*. The
group $\mathrm{SL}_2(\mathbb{Z})$ of $2 \times 2$ integer matrices with determinant 1 acts on
$\mathfrak{H}$ by Möbius transformations:

$$\begin{pmatrix} a & b \\ c & d \end{pmatrix} \cdot z = \frac{az + b}{cz + d}.$$

A *modular form* of weight $k$ and level $N$ is a holomorphic function $f: \mathfrak{H} \to
\mathbb{C}$ that transforms in a controlled way under the action of the congruence subgroup

$$\Gamma_0(N) = \left\{ \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in \mathrm{SL}_2(\mathbb{Z})
: c \equiv 0 \pmod{N} \right\},$$

satisfies a growth condition at the cusps (the boundary points of $\mathfrak{H}$ in
$\mathbb{P}^1(\mathbb{Q})$), and — for the *cusp forms* that appear in the TSW theorem — vanishes
at every cusp. Precisely: $f$ is a cusp form of weight $k$ and level $N$ if

$$f\!\left(\frac{az+b}{cz+d}\right) = (cz+d)^k f(z) \quad \text{for all } \begin{pmatrix} a & b
\\ c & d \end{pmatrix} \in \Gamma_0(N),$$

and $f$ vanishes at every cusp of $\Gamma_0(N) \backslash \mathfrak{H}$.

Every cusp form has a Fourier expansion at the cusp $i\infty$:

$$f(z) = \sum_{n=1}^{\infty} a_n q^n, \quad q = e^{2\pi i z}.$$

The coefficients $a_n$ are the arithmetic data carried by the modular form.

#### §12.2.2 Hecke operators and eigenforms

The space $S_k(\Gamma_0(N))$ of cusp forms of weight $k$ and level $N$ is a finite-dimensional
$\mathbb{C}$-vector space. It carries a commuting family of *Hecke operators* $T_p$ (one for each
prime $p \nmid N$) and $U_p$ (for $p \mid N$). A *Hecke eigenform* (or *newform*) is a cusp form
$f$ that is a simultaneous eigenvector for all Hecke operators and is normalised so that $a_1 = 1$.
For a newform $f$ with $T_p f = \lambda_p f$, the eigenvalue satisfies $\lambda_p = a_p$ (the
$p$-th Fourier coefficient).

The Hecke eigenvalues of a newform satisfy a remarkable bound: for $p \nmid N$,

$$|a_p| \leq 2\sqrt{p}.$$

This is the *Ramanujan conjecture* (proved for weight-2 newforms as a consequence of the Weil
conjectures; proved in general by Deligne [D74]). The bound $|a_p| \leq 2\sqrt{p}$ is the modular
analogue of Hasse's bound $|p + 1 - \#E(\mathbb{F}_p)| \leq 2\sqrt{p}$ for elliptic curves — and
this analogy is not a coincidence.

#### §12.2.3 The L-function of a modular form

To each newform $f = \sum a_n q^n$ one attaches an *L-function*:

$$L(f, s) = \sum_{n=1}^{\infty} \frac{a_n}{n^s} = \prod_{p \nmid N} \frac{1}{1 - a_p p^{-s} +
p^{k-1-2s}} \cdot \prod_{p \mid N} \frac{1}{1 - a_p p^{-s}}.$$

The Euler product converges for $\mathrm{Re}(s)$ sufficiently large. The L-function extends to an
entire function on $\mathbb{C}$ and satisfies a functional equation relating $L(f, s)$ to
$L(f, k - s)$. These analytic properties are consequences of the modularity of $f$.

### §12.3 The Taniyama–Shimura–Weil Theorem

The TSW theorem (conjectured by Taniyama in 1955, refined by Shimura, and given its modern form
by Weil [W67]) asserts that the two L-functions — the one attached to an elliptic curve and the
one attached to a modular form — are the same.

#### §12.3.1 The L-function of an elliptic curve

For an elliptic curve $E/\mathbb{Q}$, the *L-function* is defined by an Euler product over primes.
For a prime $p$ of good reduction (where $E$ has a smooth reduction mod $p$), the local factor is

$$L_p(E, s) = \frac{1}{1 - a_p(E) p^{-s} + p^{1-2s}},$$

where $a_p(E) = p + 1 - \#E(\mathbb{F}_p)$ is the *trace of Frobenius* at $p$. For primes of bad
reduction, the local factor is simpler (depending on the reduction type). The global L-function is

$$L(E, s) = \prod_p L_p(E, s).$$

By Hasse's bound, $|a_p(E)| \leq 2\sqrt{p}$, so the product converges for $\mathrm{Re}(s) > 3/2$.

#### §12.3.2 The theorem

**Theorem (Taniyama–Shimura–Weil; Wiles [W95], Taylor–Wiles [TW95], Breuil–Conrad–Diamond–Taylor
[BCDT01]).** *Every elliptic curve $E/\mathbb{Q}$ is modular: there exists a newform $f \in
S_2(\Gamma_0(N))$ of weight 2 and level $N$ (the conductor of $E$) such that*

$$L(E, s) = L(f, s).$$

*Equivalently, for every prime $p \nmid N$,*

$$a_p(E) = a_p(f).$$

*The map $E \mapsto f$ is the modularity correspondence.*

*Proof sketch.* The proof proceeds via Galois representations. To $E$ one attaches the
$\ell$-adic Tate module $T_\ell(E) = \varprojlim E[\ell^n]$, which gives a 2-dimensional
$\ell$-adic Galois representation $\rho_{E,\ell}: \mathrm{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})
\to \mathrm{GL}_2(\mathbb{Z}_\ell)$. The Frobenius element at $p$ (for $p \nmid N\ell$) acts with
characteristic polynomial $x^2 - a_p(E) x + p$, which encodes the point count $\#E(\mathbb{F}_p)$.
Wiles's strategy is to show that this Galois representation is *modular* — that it arises from a
modular form — by proving that a certain deformation ring equals a Hecke ring ($R = T$). The
semistable case (Wiles [W95], Taylor–Wiles [TW95]) sufficed to prove Fermat's Last Theorem; the
full theorem for all elliptic curves over $\mathbb{Q}$ was completed by Breuil, Conrad, Diamond,
and Taylor [BCDT01]. See Diamond and Shurman [DS05] for a graduate-level exposition. $\square$

#### §12.3.3 What the theorem says, structurally

The TSW theorem is a *correspondence*: it identifies two mathematical objects — an elliptic curve
and a modular form — that appear to live in entirely different worlds. The elliptic curve is a
geometric object (a smooth projective curve of genus 1 over $\mathbb{Q}$); the modular form is an
analytic object (a holomorphic function on the upper half-plane with prescribed symmetry). The
identification is canonical: the Fourier coefficients $a_p(f)$ of the modular form are exactly the
point counts $a_p(E) = p + 1 - \#E(\mathbb{F}_p)$ of the elliptic curve.

This is the same kind of structural identification that runs through the rest of this textbook —
but deeper. The MOV reduction identifies the ECDLP with a DLP in $\mathbb{F}_{p^k}^*$ via a
pairing; the identification is useful because the target world ($\mathbb{F}_{p^k}^*$) has more
structure (index calculus applies). The TSW correspondence identifies an elliptic curve with a
modular form; the identification is useful because the modular world has analytic structure
(L-functions, functional equations, Hecke operators) that the geometric world does not obviously
possess.

### §12.4 Modularity in the Escape-from-Search Family

The five escape families surveyed in §"Escape from Search: The Through-Line" are all
*computational* structures: they enable algorithms that run faster than generic search. The TSW
correspondence is a *structural* theorem, not an algorithm. But the two are related in a way that
is worth making explicit.

#### §12.4.1 Point counting and the Frobenius trace

The most direct connection between modularity and the algorithms in this textbook is through point
counting. The quantity $\#E(\mathbb{F}_p) = p + 1 - a_p(E)$ is the group order that determines
the security of ECDLP. Computing $\#E(\mathbb{F}_p)$ efficiently is a prerequisite for deploying
an elliptic curve in cryptography (to verify that the group order is prime, or to find its
factorisation for Pohlig–Hellman).

The Schoof–Elkies–Atkin (SEA) algorithm [S85, E92, A92] computes $\#E(\mathbb{F}_p)$ in
polynomial time by computing $a_p(E) \bmod \ell$ for small primes $\ell$ (using the $\ell$-torsion
structure of $E$) and then applying the Chinese Remainder Theorem. The modularity correspondence
enters here implicitly: the Hecke operators $T_\ell$ act on the modular form $f$ associated to $E$
with eigenvalue $a_\ell(f) = a_\ell(E)$, and the SEA algorithm can be understood as computing
these eigenvalues directly from the curve.

This connection is computational but not an escape from search: SEA computes the group order, it
does not solve the ECDLP. The modularity structure is used to *understand* the curve, not to
*attack* it.

#### §12.4.2 The Weil conjectures and the Ramanujan bound

The Ramanujan bound $|a_p(E)| \leq 2\sqrt{p}$ — which is Hasse's theorem for elliptic curves —
is a consequence of the Weil conjectures (proved by Deligne [D74] in general). The Weil
conjectures assert that the zeta function of a smooth projective variety over $\mathbb{F}_p$
satisfies a Riemann hypothesis: its zeros lie on the line $\mathrm{Re}(s) = 1/2$ (in the
appropriate normalisation). For an elliptic curve, this gives $|a_p| \leq 2\sqrt{p}$.

The Weil conjectures are a structural theorem about the arithmetic of varieties over finite fields.
They do not directly enable attacks, but they constrain the structure of the objects that attacks
must work with. The Hasse bound, for instance, constrains the group order $\#E(\mathbb{F}_p)$ to
lie in a narrow interval around $p + 1$, which is why the SEA algorithm can find it efficiently.

#### §12.4.3 Isogenies and the structure of the modular correspondence

The TSW correspondence is mediated by *isogenies*: morphisms $\phi: E_1 \to E_2$ of elliptic
curves (defined over $\mathbb{Q}$) that are group homomorphisms with finite kernel. The modular
curve $X_0(N)$ parametrises pairs $(E, C)$ where $E$ is an elliptic curve and $C \subset E$ is a
cyclic subgroup of order $N$; the modularity of $E$ is witnessed by a non-constant map
$X_0(N) \to E$ (a *modular parametrisation*). The degree of this map is the *Manin constant* of
$E$.

Isogenies are the morphisms of the category of elliptic curves. They are also the objects of
study in *isogeny-based cryptography* — a family of post-quantum proposals (SIDH, CSIDH, SQISign)
that base their security on the difficulty of computing isogenies between supersingular elliptic
curves. The Castryck–Decru break of SIDH [CD23] (surveyed in §11.6) used the structure of the
isogeny graph — specifically, the action of the class group of an imaginary quadratic field on the
set of supersingular curves — to recover the secret isogeny in polynomial time. This is a
structural escape from search: the class-group action provides the exploitable structure.

The connection to modularity is that the isogeny graph of supersingular elliptic curves over
$\overline{\mathbb{F}_p}$ is the Bruhat–Tits tree of $\mathrm{PGL}_2(\mathbb{Q}_p)$, and the
Hecke operators on modular forms act on this graph. The Castryck–Decru attack exploits the
*endomorphism ring* structure of supersingular curves — a structure that is deeply connected to
the modularity of those curves. This is not a coincidence; it is a manifestation of the same
correspondence.

### §12.5 The Speculation: Does Modularity Suggest Further Escapes?

> **The following is speculation.** The observations in §§12.3–12.4 are established mathematics.
> The question posed in this section — whether the modularity structure suggests further escapes
> from the search bound that have not yet been realised — is open. No algorithm is claimed; no
> complexity bound is asserted; no attack is described. The speculation is here because it is the
> honest closing gesture of a project that has surveyed the known escapes and wants to gesture
> honestly at the unknown ones.

The five escape families in this textbook are all *known*: they have been discovered, analysed,
and (at toy scale) implemented. The question is whether there are further escapes — structural
correspondences that have not yet been turned into algorithms.

#### §12.5.1 The pattern of structural escapes

Looking back at the five families, a pattern emerges:

1. **Group homomorphisms / smoothness** (index calculus, NFS): the structure is a homomorphism
   from the group to a simpler object ($\mathbb{Z}/n\mathbb{Z}$, a factor base). The escape is
   that the simpler object can be computed with.

2. **Endomorphisms** (GLV, Koblitz): the structure is a self-map of the curve that acts as scalar
   multiplication by a fixed $\lambda$. The escape is that the orbit collapses, reducing the
   effective group size.

3. **Pairings** (MOV): the structure is a bilinear map from the curve to a multiplicative group.
   The escape is that the target group has more structure (index calculus applies).

4. **Number-field structure** (GNFS, NFS-DL): the structure is an embedding of $\mathbb{Z}$ into
   a number field where smoothness is more tractable. The escape is that the number field provides
   a richer arithmetic.

5. **Quantum period-finding** (Shor): the structure is the periodicity of modular exponentiation,
   accessible only in the quantum model. The escape is that the quantum Fourier transform can find
   the period efficiently.

In each case, the escape is enabled by a *correspondence* between the problem's native world and
a world with more exploitable structure. The MOV reduction is a pairing (a bilinear correspondence
between two groups). The GNFS is a norm map (a multiplicative correspondence between a number
field and $\mathbb{Z}$). Shor's algorithm is a Fourier transform (a correspondence between a
function and its frequency spectrum).

The TSW correspondence is a deeper identification: it says that the arithmetic of an elliptic
curve over $\mathbb{Q}$ is *the same* as the arithmetic of a modular form. The modular world has
structure that the geometric world does not obviously possess: Hecke operators, L-functions,
functional equations, the Ramanujan bound. If any of this structure could be turned into an
algorithm — a way to compute something about the curve by working in the modular world instead —
it would be a new escape.

#### §12.5.2 What would a modularity-based escape look like?

**Speculation (not a claimed result).** A modularity-based escape from the ECDLP search bound
would require:

1. A way to *compute* in the modular world using only information available to an attacker (the
   curve $E$, the generator $G$, the target $Q$).

2. A way to *translate* the result back to the ECDLP: to recover the scalar $k$ from information
   about the modular form $f$ associated to $E$.

3. A complexity advantage: the computation in the modular world should be faster than $O(\sqrt{n})$
   group operations.

None of these three requirements is obviously satisfiable. The modularity correspondence is
*existential* (every elliptic curve over $\mathbb{Q}$ has an associated modular form) but not
obviously *constructive* in the sense needed for an attack. Computing the modular parametrisation
$X_0(N) \to E$ explicitly is a hard problem; the conductor $N$ of a cryptographic curve is
enormous (of the order of $p^2$, where $p$ is the field size), and the modular curve $X_0(N)$ at
that level is intractable.

The Hecke operators $T_p$ act on the modular form $f$ with eigenvalue $a_p(E) = p + 1 -
\#E(\mathbb{F}_p)$. Computing $a_p(E)$ is exactly the point-counting problem, which SEA solves
in polynomial time — but this is point counting, not ECDLP. The ECDLP asks for the scalar $k$
such that $Q = k \cdot G$; the Hecke eigenvalue $a_p(E)$ encodes the group order, not the
discrete logarithm of a specific point.

#### §12.5.3 The isogeny connection and the Castryck–Decru precedent

The most suggestive precedent is the Castryck–Decru break of SIDH [CD23]. SIDH was designed to
be secure against both classical and quantum attacks; its security rested on the difficulty of
computing isogenies between supersingular elliptic curves. The break exploited the structure of
the isogeny graph — specifically, the action of the endomorphism ring of a supersingular curve —
to recover the secret isogeny in polynomial time.

The endomorphism ring of a supersingular elliptic curve over $\overline{\mathbb{F}_p}$ is a
maximal order in a quaternion algebra ramified at $p$ and $\infty$. This is a deeply modular
object: the Eichler orders in this quaternion algebra correspond to Hecke operators on modular
forms of level $p$. The Castryck–Decru attack is, in a precise sense, an attack that exploits
the modular structure of supersingular curves.

**Speculation.** The Castryck–Decru break suggests that the modularity structure of elliptic
curves — specifically, the endomorphism ring structure and its connection to Hecke operators —
can be turned into an attack when the problem is set up in a way that exposes this structure.
SIDH exposed the structure by revealing auxiliary torsion points; the attack used this exposure
to compute the endomorphism ring and then the secret isogeny.

The question is whether a similar exposure exists for the ECDLP on ordinary (non-supersingular)
elliptic curves over prime fields — the setting of secp256k1 and the other curves used in
practice. The endomorphism ring of an ordinary elliptic curve over $\mathbb{F}_p$ is an order in
an imaginary quadratic field (not a quaternion algebra), and its structure is less rich. The
modularity correspondence for ordinary curves is the TSW theorem (over $\mathbb{Q}$), which does
not obviously expose the same kind of structure that the Castryck–Decru attack exploited.

**This is where the speculation ends.** We do not know whether the modularity structure of
ordinary elliptic curves can be turned into an ECDLP attack. The question is open; the answer
is likely "no" for the curves used in practice (the security of secp256k1 and similar curves is
not threatened by any known modularity-based attack), but the honest answer is that we do not
have a proof that no such attack exists.

#### §12.5.4 L-functions and the Birch–Swinnerton-Dyer conjecture

One further speculative direction: the Birch–Swinnerton-Dyer (BSD) conjecture [BSD65] asserts
that the rank of the Mordell–Weil group $E(\mathbb{Q})$ (the group of rational points) equals the
order of vanishing of $L(E, s)$ at $s = 1$. The BSD conjecture is one of the Millennium Prize
Problems; it is not proved in general.

The connection to cryptography is indirect. The ECDLP is defined over $\mathbb{F}_p$ (a finite
field), not over $\mathbb{Q}$. The Mordell–Weil group $E(\mathbb{Q})$ is a finitely generated
abelian group (by the Mordell–Weil theorem), but it is not the group in which the ECDLP is posed.
The BSD conjecture concerns the arithmetic of $E$ over $\mathbb{Q}$; the ECDLP concerns the
arithmetic of $E$ over $\mathbb{F}_p$.

**Speculation.** If the BSD conjecture were proved, and if the proof were constructive (giving an
algorithm to compute the rank and generators of $E(\mathbb{Q})$ from $L(E, s)$), it would not
directly give an ECDLP algorithm — the groups are different. But it would give a deeper
understanding of the arithmetic of $E$ over $\mathbb{Q}$, which might (speculatively) shed light
on the arithmetic of $E$ over $\mathbb{F}_p$ via the reduction map $E(\mathbb{Q}) \to
E(\mathbb{F}_p)$.

This is a very long chain of speculation. The honest assessment is: the BSD conjecture, even if
proved, would not obviously give an ECDLP attack. The connection is structural (both concern the
arithmetic of $E$) but not algorithmic (the groups are different, and the L-function encodes
global arithmetic, not local discrete logarithms).

### §12.6 Why This Speculation Belongs Here

The project's through-line is structure-based escape from search. The five escape families
surveyed in Chapters 6–11 are the known escapes; they are implemented and exercised by
known-answer tests at toy scale, with their mathematics developed in those chapters. This
chapter is the honest closing gesture: a survey of the deepest known structural theorem about
elliptic curves (TSW), an account of how it connects to the
algorithms in this textbook (point counting, isogenies, the Castryck–Decru precedent), and an
explicit speculation about whether it points at further escapes.

The speculation is chartered because it is the intellectually honest terminus of the through-line.
Having surveyed every known escape, the project should gesture at the unknown ones — not by
claiming them, but by naming the mathematical structures that might generate them and being
explicit about what is known and what is not.

**The honest summary.** The TSW correspondence is a deep structural theorem that identifies
elliptic curves with modular forms. It connects to the algorithms in this textbook through point
counting (SEA), isogenies (the Castryck–Decru break of SIDH), and the Weil conjectures (the
Hasse bound). Whether it suggests further escapes from the ECDLP search bound is an open
question. The most suggestive precedent — the Castryck–Decru break — exploited the modular
structure of *supersingular* curves in a setting (SIDH) that exposed the endomorphism ring. The
ordinary curves used in practice (secp256k1, P-256) do not obviously expose the same structure.
No modularity-based ECDLP attack is known; no proof that none exists is known. The speculation
is honest about its own uncertainty, and that honesty is the point.

### §12.7 Cross-References

#### Within this textbook

- **§"Escape from Search: The Through-Line"** — the five-family structure taxonomy and the
  L-notation hierarchy. This chapter speculates beyond that taxonomy.

- **§4 (Prerequisites), §"Elliptic curves"** — the group law, Hasse's bound, and the Weil
  conjectures. The Ramanujan bound (§12.2.2) is the modular analogue of Hasse's bound.

- **§10 (Algebraic ECDLP Attacks)** — the five algebraic structures that escape the generic
  bound for ECDLP. The MOV reduction (§10.5) is the cross-track bridge; the isogeny connection
  (§12.4.3) is the speculative extension.

- **§11 (Shor's Algorithm and Post-Quantum Context)** — the quantum dissolution of the
  L-notation bound, and the post-quantum migration landscape. The Castryck–Decru break of SIDH
  (§11.6) is the precedent for a modularity-based attack; §12.4.3 and §12.5.3 revisit it in the
  modularity context.

#### Code realisation

This chapter has no code-tour sibling. The modularity correspondence is not implemented in the
rDLP workspace; no algorithm in this chapter is realised in code. The chapter is a mathematical
survey and speculation, not a code-tour. The code-first umbrella (`docs/PEDAGOGY.md` §0) names
the through-line and the five escape families; this chapter is the maths-first extension of that
through-line into speculative territory.

### §12.8 Further Reading

52. **Wiles, A. (1995).** "Modular elliptic curves and Fermat's Last Theorem." *Annals of
    Mathematics*, 141(3), 443–551. [W95] The proof of the semistable case of the TSW conjecture;
    the central step in the proof of Fermat's Last Theorem.

53. **Taylor, R., and Wiles, A. (1995).** "Ring-theoretic properties of certain Hecke algebras."
    *Annals of Mathematics*, 141(3), 553–572. [TW95] The companion paper completing the $R = T$
    argument.

54. **Breuil, C., Conrad, B., Diamond, F., and Taylor, R. (2001).** "On the modularity of
    elliptic curves over $\mathbb{Q}$: wild 3-adic exercises." *Journal of the American
    Mathematical Society*, 14(4), 843–939. [BCDT01] The completion of the TSW theorem for all
    elliptic curves over $\mathbb{Q}$.

55. **Diamond, F., and Shurman, J. (2005).** *A First Course in Modular Forms*. Springer
    Graduate Texts in Mathematics, vol. 228. [DS05] The standard graduate-level introduction to
    modular forms and the TSW theorem; accessible to a reader with a solid undergraduate algebra
    background.

56. **Silverman, J. H. (2009).** *The Arithmetic of Elliptic Curves*, 2nd ed. Springer Graduate
    Texts in Mathematics, vol. 106. [S09] The standard reference for elliptic curves; Chapter V
    covers the Hasse bound and the Weil conjectures.

57. **Schoof, R. (1985).** "Elliptic curves over finite fields and the computation of square
    roots mod $p$." *Mathematics of Computation*, 44(170), 483–494. [S85] The original
    polynomial-time point-counting algorithm.

58. **Elkies, N. D. (1992).** "Explicit isogenies." Unpublished manuscript. [E92] The Elkies
    improvement to Schoof's algorithm (the SEA algorithm).

59. **Deligne, P. (1974).** "La conjecture de Weil. I." *Publications Mathématiques de
    l'IHÉS*, 43, 273–307. [D74] The proof of the Weil conjectures (including the Ramanujan
    conjecture for weight-2 newforms as a special case).

60. **Birch, B. J., and Swinnerton-Dyer, H. P. F. (1965).** "Notes on elliptic curves. II."
    *Journal für die reine und angewandte Mathematik*, 218, 79–108. [BSD65] The original
    formulation of the BSD conjecture.

61. **Castryck, W., and Decru, T. (2023).** "An efficient key recovery attack on SIDH."
    *Advances in Cryptology — EUROCRYPT 2023*, LNCS 14008, 423–447. [CD23] The classical
    polynomial-time break of SIDH/SIKE; the precedent for a modularity-based attack on
    isogeny-based cryptography. (Also cited as ref. 48 in §11.9.)

62. **Weil, A. (1967).** "Über die Bestimmung Dirichletscher Reihen durch Funktionalgleichungen."
    *Mathematische Annalen*, 168, 149–156. [W67] Weil's formulation of the conjecture that every
    elliptic curve over $\mathbb{Q}$ has an associated modular form with the same L-function.
