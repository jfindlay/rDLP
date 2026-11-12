# rDLP — Exposition corrections (non-NFS-DL)

Execution graph for the ROADMAP scope "Mathematical exposition", excluding NFS-DL claims. The
previous (completed) solver-outcome follow-up plan is preserved in git at `76d26d3:docs/PLAN.md`.

## Purpose (design intent)

Make the mathematical exposition outside NFS-DL say only what the mathematics, the cited sources
and the code support. The GNFS complexity is derived, not patched: one smoothness estimate, one
norm-size model, one optimisation, ending at the published constant. Non-monic GNFS prose
distinguishes the homogeneous value F(a, b) = c·N(a − bα) from the norm and never applies an index
[O_K : ℤ[α]] to a non-integral α. "Bad prime" (p | lc(f)·disc(f), what the code flags) is kept
distinct from "index divisor". Linear algebra yields candidates; exact reconstruction decides.
Testing documentation describes the tests that exist. Heuristic speedups are labelled as heuristic.
Quantum-resource statements carry their source's evidence scope. Where a claim can be pinned by a
test, it is pinned.

Scope boundary: `docs/MATHEMATICS.md` (§§ Prerequisites, α-substrate smoothness, Pollard rho GLV,
GNFS §3–§7, §9.7's restatement of §7 only, §11 resources and bibliography), `docs/PEDAGOGY.md`,
`gnfs/docs/PEDAGOGY.md`, `shor/docs/PEDAGOGY.md`, `README.md`, `docs/development.md`,
`docs/BENCHMARKS.md` (counts only), rustdoc/comments in `gnfs/src/{sieve,polyselect,sqrt}`,
`shared/numfield/src/dedekind.rs`, `rho/src/ecdlp/glv.rs`, plus two test files. Out of scope:
NFS-DL capability/equation claims (`docs/MATHEMATICS.md` §9.1–§9.6, gnfs PEDAGOGY NFS-DL parts,
`gnfs/src/dl/**`) — they wait for the NFS-DL design; `docs/SURVEY.md` (retirement scope); README
navigation paths and track vocabulary (project-wide anneal scope); behaviour changes to any solver;
new benchmarks or measurements; new crate or dependency edges.

## Verify gate

- VERIFY_TYPES: `cargo check --workspace`
- VERIFY_TEST: `cargo test --workspace`
- Cadence (authority: `docs/development.md`, "Verification cadence"): VERIFY_TYPES and
  VERIFY_TEST once on the final inputs of each completed row (rustdoc edits count as potentially
  executable input); `cargo coverage` (80% line floor) at the ◆ boundary and at final delivery,
  deferred with that checkpoint named at interior rows; `cargo dep-policy` and `cargo audit` when
  dependency inputs change and on the first active run of each day. A row with a deferred check
  is complete under its per-unit gate, never a full-gate pass.
- Not gates: `cargo fmt-check`, `cargo lint` (unadopted; do not fix drive-by).
- Oracles: none required. Row 3 may cross-check its KAT facts with `gp` manually; never a gate.

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 @architect | Re-derive the GNFS L[1/3, (64/9)^{1/3}] complexity from one smoothness estimate | A | T1 (juncture T0) | C-NonMonicNorm | `docs/MATHEMATICS.md` (Prerequisites smooth-number density, α-substrate §B-smooth density, §4 CEP estimate, §7, L-notation synthesis table, §9.7 restatement of §7), `docs/PEDAGOGY.md` (§7 cites), `gnfs/docs/PEDAGOGY.md` (complexity passages) |
| 2 | State the non-monic homogeneous value as c·N(a − bα) and localise the index | B | T1 | C-NonMonicNorm | `docs/MATHEMATICS.md` (§3, Prerequisites Dedekind), `gnfs/docs/PEDAGOGY.md`, `gnfs/src/sieve/{norms,mod,line}.rs`, `gnfs/src/polyselect/murphy.rs` |
| 3 | Separate bad primes from index divisors and name `number_field()`'s consumers | B | T1 | C-NonMonicNorm | `gnfs/src/sieve/factor_base.rs`, `shared/numfield/src/dedekind.rs` (rustdoc), `gnfs/src/polyselect/{mod,base_m}.rs` (rustdoc), `gnfs/docs/PEDAGOGY.md`, `shared/numfield/tests/` (one KAT) |
| 4 | Pin the trivial-congruence policy and say linear algebra yields candidates | B | T1 | — | `gnfs/src/sqrt/assembly.rs`, `README.md` (GNFS intro), `gnfs/docs/PEDAGOGY.md` (~:4447) |
| 5 | Make the testing documentation describe the tests that exist | B | T1 | — | `docs/development.md`, `docs/BENCHMARKS.md`, `shor/docs/PEDAGOGY.md`, `shared/numth/src/ecm.rs` (one unit test), `rho/tests/ssa_kat.rs` (header) |
| 6 | Label the GLV, negation-map and Koblitz ratios as expected, not measured | B | T1 | — | `docs/MATHEMATICS.md` (Pollard rho orbit maps), `docs/PEDAGOGY.md` (rho optimisation tour), `README.md`, `rho/src/ecdlp/glv.rs` |
| 7 | Rework the quantum-resource exposition from primary sources | B | T1 | — | `docs/MATHEMATICS.md` (§11.3.4, §11.4.4, §11.7, Shor bibliography), `shor/docs/PEDAGOGY.md` (references, §6.3), `docs/PEDAGOGY.md` (Shor rows) |
| 8 ◆ | Correct surviving claims the exposition rows contradicted | I | T1 (◆ fork T0) | all | durable files touched by rows 1–7 |
| 9 | Anneal cleanup (skippable) | I | — | none | durable files touched by rows 1–8 |

Order: 1 → 2 → 3 (shared GNFS notation and the 6,000-line gnfs PEDAGOGY, edited serially); 4, 5, 6
and 7 are independent; row order is the run order. Row 1 is the inflection point: the juncture
fork derives the calculation (C-GnfsComplexity) before row 1 writes prose.

## Session detail

### 1 @architect — Re-derive the GNFS L[1/3, (64/9)^{1/3}] complexity from one smoothness estimate

Deliverable per C-GnfsComplexity (frozen by the juncture). Rewrite MATHEMATICS §7 (~:2000–2222) as
one derivation: norm sizes for both sides with degree d and coefficient bound from base-m
selection (carrying c per C-NonMonicNorm, or stating the monic simplification explicitly); u for
each side; one smoothness estimate ρ(u) = u^{−u(1+o(1))} and its L-form; sieve and linear-algebra
cost; optimisation over d and B yielding (64/9)^{1/3}. Constants are never absorbed into o(1)
unless they are o(1) in the exponent's leading term. Remove "Wait — let us be more careful",
"Hmm…", the 24^{1/3} attempt and the √(2/3) / 2√6/3 one-sided result unless the derivation keeps
it as an explicitly labelled simplification. Reconcile the Prerequisites (~:509–540), α-substrate
(~:1189–1211) and §4 (~:1688–1706) statements of the smoothness estimate to the same form. Fix the
summary table (~:2207–2222, "(8/9)^{1/3} (full)") and the L-notation synthesis table (~:262–305).
Keep the `§7` heading (≥15 consumers cite it); update §9.7's restatement of what §7 establishes
(~:2791–2810) only where the result's statement changed — no NFS-DL capability text. State the
heuristic status (CEP applied to polynomial values is a heuristic, not a theorem).

KAT substitute (prose row): every step names its reference (section/page of the source read) or is
labelled "derived here"; the juncture's step list is the checkable skeleton; each consumer cite of
§7 is re-read and still true.

Subtleties: the canonical references (Buhler–Lenstra–Pomerance 1993; Crandall–Pomerance §6.2;
Lenstra–Lenstra 1993) may not be readable in full. If not, the derivation must be self-contained
and checkable step by step, with the source cited only for the published final statement and
labelled "recalled".

### 2 — State the non-monic homogeneous value as c·N(a − bα) and localise the index

Deliverable per C-NonMonicNorm. Across gnfs PEDAGOGY (~:16, 216, 242, 582, 679, 921, 959, 3246,
4447, 4521, 4811, 5814), rustdoc in `sieve/norms.rs` (:9, 24, 50, 76; vague :12–13), `sieve/mod.rs`
:17, `sieve/line.rs` :157, `polyselect/murphy.rs` (:14, 178, 288): F(a, b) = b^d·f(a/b) is "the
algebraic value F" or "homogeneous norm value", equal to c·N(a − bα); "norm" alone names N. In
MATHEMATICS, [O_K : ℤ[α]] (~:387–392) is stated for monic f, with the non-monic replacement (ℤ[cα]
or localisation away from c) given once; the ℤ[α]-ideal uses at ~:1574–1679 get the convention of
the correct passages (~:1559–1564, 1653–1663). Identifiers unchanged. Do not touch §7 (row 1 owns
it) or NFS-DL chapters.

KAT substitute: each corrected rustdoc names the function computing the stated quantity and that
function's existing test; any sieve correctness argument that needs N rather than F is flagged as
a discovery, not reworded.

### 3 — Separate bad primes from index divisors and name `number_field()`'s consumers

Deliverable per C-BadPrime. `factor_base.rs` module comment (~:31–43, 61–82, 131): bad prime = p |
lc(f)·disc(f) as `is_bad_prime` computes it on non-monic f; index divisors are a subset; a bad
prime may be maximal (x² − 2 at 2); drop "p | c can be invisible". `dedekind.rs` rustdoc
(`discriminant` ~:192–197, `is_bad_prime`): state what it returns for non-monic f (lc(f)·disc(f)
up to sign) since `FactorBase::new` passes non-monic f. gnfs PEDAGOGY bad-prime passages
(~:879–894, 4025, 4967–4972) likewise. `number_field()` passages (gnfs PEDAGOGY ~:687–690,
4511–4512, 4812; `polyselect/mod.rs` :248–258; `base_m.rs` :17): its consumers are NFS-DL relation
code and tests; the factor base does its own root search; the square-root stage does not call it.

KATs (in `shared/numfield/tests/`): `is_bad_prime` on a non-monic f flags a prime dividing lc(f)
but not disc of the monic transform; `is_bad_prime(x² − 2, 2)` is true while `dedekind_factor`
reports ℤ[α] 2-maximal (use its existing report; if it exposes no maximality result, assert
the factorisation shape instead and say so).

### 4 — Pin the trivial-congruence policy and say linear algebra yields candidates

Deliverable per C-TrivialCongruence. `AttemptOutcome::TrivialCongruence` rustdoc (assembly.rs
:171) and the module doc (~:19–21) state: x ≡ ±y (mod N); a divisor that arises only from that
case is deliberately not reported (link the detailed explanation at ~:84–91). README :39–40 and
gnfs PEDAGOGY ~:4447–4449: linear algebra produces candidate dependencies; a congruence of squares
exists only after both square roots are reconstructed and verified.

KAT: `factor_from_congruence` (or the narrowest public/private seam that classifies attempts) on
N = 35, x = y = 5 returns the trivial-congruence outcome, not a factor, though gcd(x + y, N) = 5.

### 5 — Make the testing documentation describe the tests that exist

Deliverable per C-TestingDocs. `docs/development.md` "KAT-first coverage" (~:144–159): describe the
actual mix (KATs primary; property tests in `shared/field/tests/sqrt_legendre_kat.rs`; opt-in
oracles), drop "no proptest/quickcheck layer" and "one file per algorithm family" where untrue.
Remove hard-coded test counts from BENCHMARKS (~:59, 118, 179, 212–215, 268–271, 329, 379–381,
472, 662–663, 743–744, 842) and `shor/docs/PEDAGOGY.md` (~:584–589), keeping each file's coverage
description. `ecm.rs` `stage2_finds_factor_stage1_misses` (~:838–857): make the test assert what
its name says (stage 1 alone fails, stage 1 + 2 succeeds) through an existing private seam with a
bounded, recorded parameter search; if no seam can express "stage 1 alone", rename and re-doc the
test to what it checks. `ssa_kat.rs` header (~:9–40): number to match the body or drop numbering.

KAT: the strengthened ECM test (or its honest rename); prose claims cite the files.

### 6 — Label the GLV, negation-map and Koblitz ratios as expected, not measured

Deliverable per C-HeuristicRatio. MATHEMATICS ~:1006–1011, PEDAGOGY ~:522, 533, 599, 613, 710–722,
745, README :336 and `glv.rs` :10: each ratio is "expected (heuristic), not measured"; the
negation-map passages note the open step-count anomaly without quantifying it beyond BENCHMARKS'
record. No measurement, no bench change.

KAT substitute: the labels cite `glv_fewer_steps_than_negmap` for what is asserted (a 2×
no-regression bound) where they mention testing.

### 7 — Rework the quantum-resource exposition from primary sources

Deliverable per C-QuantumEvidence. Read full texts where obtainable: Babbush et al.
arXiv:2603.28846v2; Cain et al. arXiv:2603.28627v1; Schrottenloher arXiv:2606.02235v1; Proos–Zalka
2003 (QIC 3(4)); Roetteler et al. 2017 (ASIACRYPT). Rewrite MATHEMATICS §11.4.4 and §11.7 (and
§11.3.4 if it states resources) so every quantitative claim names its source and scope: logical vs
physical qubits; per-point-addition vs full-algorithm counts; architecture/runtime assumptions;
secp256k1 vs generic prime field; Babbush's withheld-circuit/zero-knowledge-proof evidence
boundary; Schrottenloher's circuits compared against Babbush, not independent estimates of the
same scope. Add bibliography entries; cite Shor's journal version (SIAM J. Comput. 26(5),
1484–1509, 1997) beside FOCS 1994 in MATHEMATICS refs 14/49 and `shor/docs/PEDAGOGY.md`. Mention
and link; do not expound new algorithms (scope boundary). Do not edit SURVEY.

KAT substitute: each number cites section/table of a source read; anything from an abstract only is
labelled so and not compared quantitatively. Record which versions were read.

### 8 ◆ — Correct surviving claims the exposition rows contradicted

Deliverable: (a) stale-claim sweep, including unchanged passages, for each corrected claim: the
one-sided √(2/3) constant and the 1/(3c) exponent; "norm" for F without c; index for non-integral
α; "p | c invisible", "bad prime ⇒ non-maximal"; `number_field()` consumers; "linear algebra
produces a congruence"; "both gcds trivial"; "no proptest"; test counts; unlabelled ratios;
unscoped quantum numbers. (b) Every consumer cite of §7 and of renamed/re-anchored sections
resolves. (c) Anneal the arc's durable files against the denylist. NFS-DL text is report-only.

KAT substitute: every changed behaviour claim cites the test that exercises it or is labelled.

### 9 — Anneal cleanup (skippable)

`/plan/anneal` over the arc's durable files: translate leaked plan coordinates and rewrap ledgered
line-length overruns. The driver may mark it `skipped` if the ◆ anneal and its fork ran clean, the
review had no scope limit, and no overruns were ledgered.

## Cross-session contracts

### C-NonMonicNorm

Prose-enforced. Defined by the existing MATHEMATICS §3–§4 passages (~:1559–1564, 1653–1663);
propagated at 2. Consumed by 1, 2, 3, 8. For f with leading coefficient c and root α:
F(a, b) = b^d·f(a/b) = c·N_{K/ℚ}(a − bα); "norm" names N; an index [O_K : ℤ[α]] is used only for
integral α, otherwise ℤ[cα] or a localisation away from c. Identifiers are not renamed.

### C-GnfsComplexity

Prose-enforced. Frozen at 1 by the inflection juncture. Consumed by 1, 8. The `§7` heading and its
position stay stable.

**Notation.** λ = log N, μ = log log N (natural logs). The definition of L_N[α, c] is kept; every
asymptotic result is written L_N[α, c + o(1)] (the o(1) inside the constant) and a probability as
1/L_N[α, γ + o(1)]. Never "· (1 + o(1))" outside an L, never L_N[−1/3, ·]. One sentence stating
this convention is added at the L-notation definition (Prerequisites).

**Sources.** READ: [Bai11] S. Bai, *Polynomial Selection for the Number Field Sieve*, PhD thesis,
ANU, 2011, §2.3, pp. 22–24 (CEP form and range as its Theorem 2.5; box |a| ≤ U, 0 < b ≤ U; bounds
|N(a − bα)| ≤ (d+1)n^{1/d}U^d, |a − bm| ≤ n^{1/d}U; random-integer assumption; factor base
B^{1+o(1)}; sieve U^{2+o(1)}, linear algebra B^{2+o(1)}, B ≈ U; d = (3^{1/3} + o(1))(λ/μ)^{1/3},
U = B = L(1/3, (8/9)^{1/3}), time L(1/3, (64/9)^{1/3}); it assumes f, g monic). Add it to the GNFS
reference list. RECALLED (cite only for the published final statement, no section number):
[BLP93] (the Leiden PDF is a scan without text, not read), [CP05] §6.2, [LLMP93]. CEP 1983 is
cited "as stated in [Bai11, Thm 2.5]". Pomerance's "A Tale of Two Sieves" was not read; do not
cite it as read. Every step below is "derived here" unless it names [Bai11].

**Step list (the checkable skeleton of §7).**

1. Polynomials and box ([Bai11] read; c carried here). Base-m: m = N^{(1+o(1))/d} (the code takes
   the least m with m^{d+1} > N; the difference from N^{1/d} is λ/d² = o(λ/d), o(1) in the
   leading term). Digits |f_i| < m; leading coefficient c = f_d ∈ [1, m); g(x) = x − m. Box
   |a| ≤ M, 0 < b ≤ M, gcd(a, b) = 1.
2. Sizes, both sides. |G(a, b)| = |a − bm| ≤ 2Mm, so log|G| ≤ λ/d + log M + O(1).
   |F(a, b)| = |b^d f(a/b)| ≤ (d+1)·m·M^d, so log|F| ≤ λ/d + d·log M + log(d+1). The quantity
   sieved and sized is F(a, b) = c·N_{K/ℚ}(a − bα) (C-NonMonicNorm), not N(a − bα); c enters
   only through max|f_i| < m, so no monic simplification is needed, and the prose says so.
3. Ansatz B = L_N[1/3, β], M = L_N[1/3, η], d = δ(λ/μ)^{1/3}. Then
   log|G| = (1/δ + o(1))λ^{2/3}μ^{1/3} and log|F| = (1/δ + δη + o(1))λ^{2/3}μ^{1/3}. The only
   absorptions are log M on the rational side and log(d+1), both O(λ^{1/3}μ^{2/3}), smaller than
   the leading term by (μ/λ)^{1/3}. δ, β, η, 3 and 3^{1/3} are never absorbed.
4. u per side: u_G = (1/(δβ) + o(1))(λ/μ)^{1/3}, u_F = ((1/δ + δη)/β + o(1))(λ/μ)^{1/3};
   log u = (1/3 + o(1))μ on both sides.
5. The one smoothness estimate ([Bai11] Thm 2.5, read): Ψ(x, x^{1/u}) = x·u^{−u(1+o(1))} as
   u → ∞, uniformly for x ≥ u^{u/(1−ε)}; here u^u = exp(O(λ^{1/3}μ^{2/3})), inside the range.
   Its L-form (derived here): if x = L_N[a, κ] and y = L_N[b, β] with 0 < b < a ≤ 1, then
   Ψ(x, y)/x = 1/L_N[a − b, (a − b)κ/β + o(1)]. GNFS instance: a = 2/3, b = 1/3 gives
   1/L_N[1/3, κ/(3β) + o(1)]. Instance x = N: 1/L_N[2/3, 2/(3β) + o(1)]. Dickman's ρ(u) is the
   fixed-u limit of Ψ(x, x^{1/u})/x; u^{−u(1+o(1))} is its large-u asymptotic; at practical u
   (≈ 3–10) ρ is evaluated numerically, never by u^{−u}.
6. Joint probability (heuristics H1, H2 below): Pr = 1/L_N[1/3, (2/δ + δη)/(3β) + o(1)]. This
   equals treating G·F as one random integer of log-size (2/δ + δη)λ^{2/3}μ^{1/3} ([Bai11]'s
   form), because log u is (1/3)μ on both sides; the prose says so.
7. Relations needed ([Bai11] read): rational primes ≤ B, degree-1 prime ideals (p, r) with
   p ≤ B, and O(λ) character columns: B^{1+o(1)} columns. Yield condition M²·Pr ≥ B^{1+o(1)},
   i.e. 2η − (2/δ + δη)/(3β) ≥ β.
8. Costs ([Bai11] read): sieving the box costs M^{2+o(1)} = L_N[1/3, 2η + o(1)] (each p ≤ B hits
   about M²/p cells; Σ_{p≤B} 1/p = O(μ)); sparse linear algebra (Wiedemann or block Lanczos) on a
   B^{1+o(1)}-square matrix with O(λ) entries per row costs B^{2+o(1)} = L_N[1/3, 2β + o(1)].
   Total L_N[1/3, 2·max(β, η) + o(1)].
9. Optimisation over B, M, d. Take M = B (η = β, sieve and linear algebra balanced, as [Bai11]
   "B ≈ U"). Step 7 becomes 3β² ≥ δβ + 2/δ; minimal β at equality. Differentiating
   3β² = δβ + 2/δ with dβ/dδ = 0 gives β = 2/δ²; substituting, 12/δ⁴ = 4/δ, so δ³ = 3.
   Hence δ = 3^{1/3}, β = 2·3^{−2/3} = (8/9)^{1/3} ≈ 0.961, and total cost
   L_N[1/3, 2β + o(1)] = L_N[1/3, (64/9)^{1/3} + o(1)], (64/9)^{1/3} ≈ 1.923.
10. Values at the optimum (summary tables use these): d = (3^{1/3} + o(1))(λ/μ)^{1/3} (the
    `optimal_degree` formula before its [3, 6] clamp); B = M = L_N[1/3, (8/9)^{1/3} + o(1)];
    |G| = N^{(1+o(1))/d} = L_N[2/3, 3^{−1/3} + o(1)]; M^d = N^{(2+o(1))/d}, so
    |F| = N^{(3+o(1))/d}; u_G = (3^{1/3}/2 + o(1))(λ/μ)^{1/3},
    u_F = ((3/2)·3^{1/3} + o(1))(λ/μ)^{1/3}; Pr = 1/L_N[1/3, (8/9)^{1/3} + o(1)] = B^{−1+o(1)}; sieve = linear algebra =
    L_N[1/3, (64/9)^{1/3} + o(1)]. (8/9)^{1/3} is the constant of B, not an alternative cost
    constant. The two sides are not the same size.
11. Why the exponent is 1/3: two balances force it. Degree balance λ/d ≍ d·log M gives
    d ≍ (λ/log B)^{1/2}; smoothness balance u·log u ≍ log B with u ≍ (λ/d)/log B then gives
    log B ≍ λ^{1/3}μ^{2/3}. With d fixed the values are N^{Θ(1)} and the exponent is 1/2. A
    quadratic-sieve contrast, if kept, says QS values are N^{1/2+o(1)} giving L_N[1/2, 1]
    (recalled, [CP05]); the "values of size N give L_N[1/2, √2]" case is labelled hypothetical
    (derived here: max(β + 1/(2β), 2β) is least at β = 1/√2).

**Heuristic status (carried in §7 and consistent with the main theorem at the GNFS chapter
head).** The result is a heuristic expected running time, not a theorem. CEP is a theorem about
all integers ≤ x. Unproved: (H1) the values G(a, b), F(a, b) over the box are B-smooth as often as
random integers of the same size; (H2) the two sides are independent; (H3) B^{1+o(1)} relations
give dependencies whose square roots yield a nontrivial congruence with probability bounded below
(the quadratic characters handle the unit and class-group obstructions).

**Removed, not kept as a labelled simplification:** "Wait — let us be more careful", "Hmm…",
the 3c / 24^{1/3} computation, the √(2/3) / 2√6/3 result (it rests on the wrong L-form and on
equal-size sides, so it is not a valid simplification), "absorbing the 3^{1/3} into o(1)", the
claim that degree optimisation and the two-sided sieve explain a gap from 2√6/3, "the sieve cost
dominates", "(8/9)^{1/3} (full)" as an optimal-c alternative, "both norms ≈ N^{1/d}" and
"N^{1/d} ≈ N^{1/3}", and "N^{1/d} = L_N[1/3, ·]".

**Reconciled statements (one estimate, step 5, everywhere):** Prerequisites smooth-number density
(~:509–541): the Dickman fixed-u statement and the CEP large-u statement with its range, then the
step 5 lemma replacing the Corollary. α-substrate B-smooth density (~:1189–1212): same estimate,
u for x = N, two-norm sentence under H1–H2. §4 (~:1650–1669 size sentences, ~:1688–1706 CEP
estimate): sizes per steps 2 and 10, joint probability per step 6. Size sentences only:
C-NonMonicNorm's c·N wording there is kept. §7: steps 1–11, the heuristic status, and a summary
table of step 10's values. L-notation synthesis table (~:262–305): GNFS row unchanged in value.
`docs/PEDAGOGY.md` GNFS row note: the exponent comes from the degree/smoothness balance (§7), the
constant is derived in §7. `gnfs/docs/PEDAGOGY.md` §60 (~:4995–5190): summary consistent with
steps 1–11, dropping 2/(dc), c = 1/√d, 4/√d, L_N[−α, ·] and the "u′ = O(1) ⇒ α = 1/2" argument.
§4 degree passages (~:174–221) are re-read and kept where consistent.

**§9.7 (restatement only).** Update ~:2791–2810. Every other §9.7 statement of a §7 quantity
(sieve/linear-algebra/smoothness bullets ~:2934–2943 and ~:2957–2961, and the comparison table's
norm-size, smoothness, sieve, linear-algebra and optimal-B rows in both columns, the p-column
being the same formula in p) counts as a restatement and gets the step 10 values. Descent
(Delta 2, c_descent ~:2880–2920), Delta 1, Schirokauer and every NFS-DL capability claim stay
untouched and are reported.

### C-BadPrime

Test + prose-enforced. To be frozen at 3. Consumed by 8. Bad prime = p | lc(f)·disc(f) (what
`is_bad_prime` computes); index divisors ⊂ bad primes; bad primes may be maximal;
`PolyPair::number_field()` serves NFS-DL relation code and tests only.

### C-TrivialCongruence

Test + prose-enforced. To be frozen at 4. Consumed by 8. `TrivialCongruence` means x ≡ ±y (mod N);
an incidental divisor in that case is deliberately not reported; N = 35, x = y = 5 pins it.

### C-TestingDocs

Prose + test-enforced. To be frozen at 5. Consumed by 8. Orientation docs carry no test counts;
development.md describes KATs, property tests and opt-in oracles as they exist; every test's doc
states what its body checks.

### C-HeuristicRatio

Prose-enforced. To be frozen at 6. Consumed by 8. Every GLV/negation-map/Koblitz ratio in durable
text is labelled expected (heuristic) and not measured.

### C-QuantumEvidence

Prose-enforced. To be frozen at 7. Consumed by 8. Every quantitative quantum-resource statement
names its source, the version read, and its scope (logical/physical, per-operation/full,
architecture); abstract-only claims are labelled and not compared.

## Progress ledger

Bindings: PLAN = `docs/PLAN.md`; VERIFY_TYPES = `cargo check --workspace`; VERIFY_TEST =
`cargo test --workspace`; boundary check = `cargo coverage`; freshness checks = `cargo
dep-policy`, `cargo audit` (separate commands). Cadence per "Verify gate".

Stash snapshot at run start: `6c65ca58890b36e241ed0d7923f4ae9a9af4641a` (single pre-existing entry).

Gate evidence:

- Freshness, 2026-11-11, tree `3002c52` clean (dependency inputs unchanged by the chain; toolchain
  rustc/cargo 1.95.0; cargo-deny 0.20.2; cargo-audit 0.22.2; advisory db `fd6f398`):
  - `cargo dep-policy` (workspace root): passed, exit 0, completed 18:29:49 −05:00, log
    `/tmp/opencode/dep-policy.log`, advisories/bans/licenses/sources ok.
  - `cargo audit` (workspace root): passed, exit 0, completed 18:29:50 −05:00, log
    `/tmp/opencode/audit.log`, db refreshed, 1458 advisories; one allowed warning
    (RUSTSEC-2026-0486, rustix 1.1.4 unsound, informational).
- Row 1: `cargo check`/`cargo test`: not run, `non-executable diff` (three Markdown files; grep of
  `*.rs`/`*.toml` finds no `include_str!`/doc-include of Markdown). `cargo coverage`: deferred to the
  ◆ boundary (row 8).
- Row 2 (rustdoc edits in `.rs`: potentially executable input), worker-run, driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, completed 2026-11-11 19:48:08
    −05:00, log `/tmp/opencode/s2-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, completed 20:03:12 −05:00, log
    `/tmp/opencode/s2-test.log`, 2032 passed, 0 failed.
  - Input identity: base `fe7132c` + unstaged diff, `git diff HEAD -- . ':!docs/PLAN.md' |
    sha256sum` = `f0e0ade8…472ffe2` (recomputed by the driver before commit: identical); no
    untracked files; rustc/cargo 1.95.0. Diff in `.rs` files verified comment-only.
  - `cargo coverage`: deferred to the ◆ boundary (row 8).
- Row 3 (rustdoc/comment edits in `.rs` plus two new tests), worker-run, driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, completed 2026-11-11 20:37:22
    −05:00, log `/tmp/opencode/s3-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, completed 20:50:11 −05:00, log
    `/tmp/opencode/s3-test.log`, 2034 passed, 0 failed (both new tests listed `ok`).
  - Input identity: base `dc84b74` + unstaged diff sha256 `4e32cd9a…6ea709` (recomputed by the
    driver before commit: identical); no untracked files; rustc/cargo 1.95.0. Non-test `.rs` diff
    verified comment-only.
  - `cargo coverage`: deferred to the ◆ boundary (row 8).
- Row 4 (rustdoc/comment edits in `.rs` plus one new unit test), worker-run, driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, completed 2026-11-11 20:52:15
    −05:00, log `/tmp/opencode/s4-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, completed 20:56:48 −05:00, log
    `/tmp/opencode/s4-test.log`, 2035 passed, 0 failed (new test listed `ok`).
  - Input identity: base `4356059` + unstaged diff sha256 `679d3a78…4e403` (recomputed by the
    driver before commit: identical); no untracked files; rustc/cargo 1.95.0. Non-test `.rs` diff
    verified comment-only.
  - `cargo coverage`: deferred to the ◆ boundary (row 8).
- Freshness, 2026-11-12 (first active run of the new day), tree `f75fe2f` clean (dependency inputs
  unchanged since `3002c52`: no manifest, lockfile, `deny.toml` or toolchain diff; cargo-deny 0.20.2;
  cargo-audit 0.22.2; advisory db `af7b3ec`, 1490 advisories):
  - `cargo dep-policy` (workspace root): passed, exit 0, completed 10:28:47 −05:00, log
    `/tmp/opencode/dep-policy-1112.log`.
  - `cargo audit` (workspace root): passed, exit 0, completed 10:28:48 −05:00, log
    `/tmp/opencode/audit-1112.log`, db refreshed; same single allowed warning as 2026-11-11.
- Row 5 (one ECM test and a test-header comment in `.rs`, Markdown docs), worker-run,
  driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, completed 2026-11-12 10:44:49
    −05:00, log `/tmp/opencode/s5-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, completed 10:49:35 −05:00, log
    `/tmp/opencode/s5-test.log`, 2035 passed, 0 failed, 36 ignored;
    `stage2_finds_factor_stage1_misses ... ok`.
  - Input identity: base `f75fe2f` + unstaged diff sha256 `138b6976…c33c0f` (recomputed by the
    driver before commit: identical); no untracked files; rustc/cargo 1.95.0. Non-test `.rs` diff
    verified to be test and comment lines only.
  - `cargo coverage`: deferred to the ◆ boundary (row 8).
- Row 6 (one module doc comment in `.rs`, Markdown docs), worker-run, driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, completed 2026-11-12 10:53:02
    −05:00, log `/tmp/opencode/s6-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, completed 11:00:00 −05:00, log
    `/tmp/opencode/s6-test.log`, 2035 passed, 0 failed; `glv_fewer_steps_than_negmap ... ok`.
  - Input identity: base `95b2bde` + unstaged diff sha256 `42b3c672…b45eb` (recomputed by the
    driver before commit: identical); no untracked files; rustc/cargo 1.95.0. `.rs` diff verified
    comment-only. No benchmark run.
  - `cargo coverage`: deferred to the ◆ boundary (row 8).
- Row 7: `cargo check`/`cargo test`: not run, `non-executable diff` (three Markdown files; no source
  includes Markdown, confirmed at row 1 and unchanged: no `.rs`/`.toml` in this diff).
  `cargo coverage`: deferred to the ◆ boundary (row 8).
- Row 8 (comment/rustdoc/test-doc edits in 14 `.rs` files plus Markdown), worker-run,
  driver-validated:
  - `cargo check --workspace` (workspace root): `reused`, exit 0, log `/tmp/opencode/s8-check.log`.
  - `cargo test --workspace` (workspace root): `reused`, exit 0, log `/tmp/opencode/s8-test.log`,
    2035 passed, 0 failed, 36 ignored (doc-tests included).
  - Input identity: base `267ced1` + unstaged diff sha256 `3f8e4514…b6972` (recomputed by the
    driver before commit: identical); no untracked files; rustc/cargo 1.95.0. `.rs` diff verified
    to contain only comment/doc lines.
  - `cargo coverage`: the ◆ boundary check, run at the boundary (see boundary note).

- ◆ boundary (rows 1–8), 2026-11-12:
  - `cargo coverage` (workspace root, `cargo llvm-cov --workspace --fail-under-lines 80`): first run
    on tree `937425c` passed (exit 0, 94.46% lines, log `/tmp/opencode/coverage-8.log`) but was
    invalidated by the boundary correction (it edited `.rs` comments and `#[ignore]` reason
    strings; classified as potentially executable). Rerun on tree `154755c` (clean): `passed`, exit
    0, completed 18:25:28 −05:00, 94.46% lines (95.17% regions, 95.62% functions; floor 80%), log
    `/tmp/opencode/coverage-8b.log`; rustc/cargo 1.95.0.
  - Boundary correction (one `@build` dispatch, commit `154755c`): `cargo check --workspace` exit 0
    (completed 18:09:03), `cargo test --workspace` exit 0 (completed 18:13:09), logs
    `/tmp/opencode/s8fix-check.log`, `s8fix-test.log`; `reused`; identity base `937425c` + diff
    sha256 `5385a456…e4d99` (recomputed by the driver before commit: identical). `.rs` diff verified
    to be comments, doc comments and `#[ignore = "…"]` reason strings only.
  - `cargo dep-policy`, `cargo audit`: valid for 2026-11-12 (advisory db `af7b3ec`; dependency inputs
    unchanged); `reused`.

Boundary note (rows 1–8):
  - Anneal gate: denylist grep over the 33 changed durable files (then 45 paths with the
    correction). Zero hits for sub-track / frozen contract / action-frame / ◆ / Cat patterns.
    Legitimate matches: bibliography keys `[S09]`, `[S93]`, `[S99]`, `[S85]` (pattern 1),
    "Appendix A.2/A.3/A.5/A.6" paper-appendix names (pattern 4). All other hits (~150 lines) are in
    `gnfs/docs/PEDAGOGY.md` §§63–71 and "Further Reading", NFS-DL text, report-only by the plan,
    unedited by this arc (genuine leaks there: "C-DLRelation (frozen D.A.1)", "D.B", "G.E", "Track-D",
    "principle-4", "the ROADMAP's named requirement"). No uncertain referent passed to the fork
    beyond those below.
  - Fork (`boundary-transform`): `still-on-intent`. §7 follows C-GnfsComplexity step by step (algebra
    re-checked: δ³ = 3, β = η = (8/9)^{1/3}, |F| = L_N[2/3, 3^{2/3}]); Dixon trial-division line
    correct (2c + 1/(2c) least at c = 1/2); ECM L_p[1/2, √2] and Pohlig–Hellman L_p[1, 1/2] clean;
    chapter-12 renumbering: every citation resolves (52–55 unique); `--help` change harmless.
  - Contradictions (C-TestingDocs: placeholder/self-checking ignored tests described as external
    oracles) and planning-vocabulary leaks ("reshard", "cross-track", two "(G.E)", numfield §10
    heading) were corrected in one bounded dispatch, `154755c`; also: Dixon cost-model labelled
    "derived here" plus a modelling-choice sentence; translated source quotations in gnfs PEDAGOGY
    now introduced as paraphrases. Driver re-grep after the correction: no `reshard`/`cross-track`
    hits in the non-NFS-DL changed files except MATHEMATICS "track" vocabulary (below).
  - Still standing (reported, not fixed): third "(G.E)" ~gnfs PEDAGOGY :5818 and "cross-track
    interface/entry point" ~:5335, 5652, 5679 (NFS-DL chapter); MATHEMATICS "track"/"cross-track"
    vocabulary (~:229, 263–275, 3258, 4832, 4987, 5460; track vocabulary is out of scope by the
    plan); `docs/development.md:74` "plan-ledger edits"; BENCHMARKS "Test coverage added" headings
    and "deliverable" (~:59, 124); `shared/numfield/docs/PEDAGOGY.md:190, 219` "contract" (probably
    API-invariant use); literal test identifiers containing "principle4"
    (`kat5_principle4_improvement_small_at_toy_scale`, `kat_b_yield_comparison_principle4_annotated`,
    `kat_c_yield_comparison_principle4_annotated`; renaming needs human sign-off); labels
    `principle-4`/"scoping principle 3" in `gnfs/src/linalg/{blockvec,operator,wiedemann}.rs`,
    `gnfs/src/filter/matrix.rs:39`, `gnfs/tests/{lattice_kat,coppersmith_kat,root_sieve_kat}.rs`
    (not enrolled); `unlabelled non-contract figures` (`line.rs` "~99%", Montgomery 2–4×, Brent 24%,
    batched GCD 2–5×, Cavallar 30–50%).
  - Scope limit (review incomplete, never clean): the fork read the diffs of all changed files over
    `3002c52..937425c` in full but not every surviving line of MATHEMATICS (~5500 lines) or gnfs
    PEDAGOGY §§1–62; unchanged prose there was checked by targeted claim and vocabulary greps only.
    [Bai11], [PZ03] and the 2026 sources were not re-fetched by the fork.

Pending checks: `cargo coverage` current for tree `154755c` (reusable at final delivery if no
executable input changes afterwards). `cargo dep-policy` and `cargo audit` fresh for 2026-11-12
only; rerun if the chain crosses a day or touches dependency inputs.

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Re-derive the GNFS L[1/3, (64/9)^{1/3}] complexity from one smoothness estimate | done | a0b12ad | C-GnfsComplexity (frozen by the inflection fork, `design-confident`). §9.7 restatements widened per the contract; allowed extras: none. Out-of-scope stale text reported to row 8 (see digest). |
| 2 | State the non-monic homogeneous value as c·N(a − bα) and localise the index | done | 6e67f92 | Propagated C-NonMonicNorm (no new freeze). Extras: none. |
| 3 | Separate bad primes from index divisors and name `number_field()`'s consumers | done | 9bfbd4e | C-BadPrime (tests `kat_is_bad_prime_non_monic_flags_a_prime_dividing_only_the_leading_coefficient`, `kat_bad_prime_that_is_maximal_has_ramified_factorisation` in `shared/numfield/tests/dedekind_kat.rs`). Extras: none. |
| 4 | Pin the trivial-congruence policy and say linear algebra yields candidates | done | 161c7b3 | C-TrivialCongruence (test `trivial_congruence_with_a_shared_factor_reports_no_divisor`, unit test in `gnfs/src/sqrt/assembly.rs`; no new public API). Extras: none. |
| 5 | Make the testing documentation describe the tests that exist | done | 24b0b88 | C-TestingDocs (strengthened ECM test `stage2_finds_factor_stage1_misses`: N = 1009·3541, B1 = 20, B2 = 200, sigma 7 found by a bounded 200-curve scan). Extras: none. |
| 6 | Label the GLV, negation-map and Koblitz ratios as expected, not measured | done | ac3e451 | C-HeuristicRatio. Extras: none. |
| 7 | Rework the quantum-resource exposition from primary sources | done | b329f83 | C-QuantumEvidence (versions read: Babbush arXiv:2603.28846v2 full; Cain 2603.28627v1 main text, appendices unread; Schrottenloher 2606.02235v1 partial; Proos–Zalka arXiv quant-ph/0301141v2 full; Roetteler arXiv:1706.06752v3 full). Extras: none. |
| 8 | Correct surviving claims the exposition rows contradicted ◆ | done | de204f8 (+ boundary correction 154755c) | Allowed extras (stale-claim sweep targets, comments/rustdoc/docs only): `docs/PEDAGOGY.md`, `shared/numfield/docs/PEDAGOGY.md`, `gnfs/src/{lib.rs,linalg/kernel.rs,sieve/{special_q,lattice}.rs,polyselect/{roots,root_sieve,coppersmith}.rs,sqrt/context.rs}`, `gnfs/tests/{special_q_kat,murphy_kat,lanczos_kat}.rs`, `rho/src/ecdlp/cli.rs`. Boundary complete: `still-on-intent` (see boundary note). |
| 9 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

### 1 — 2026-11-11
Discovery/flex: the inflection fork (`design-confident`, reading [Bai11] §2.3) found the old one-norm
L-form L_N[−1/3, −1/(3c)] wrong and replaced it by the two-sided form; (8/9)^{1/3} is the constant
of B, not an alternative cost constant; the √(2/3) / 2√6/3 result was removed, not kept as a
simplification. The contract reads every §9.7 statement of a §7 quantity as a restatement (wider
than ~:2791–2810); row 1 followed it. Row 1's worker named hypotheses (H1)–(H3) locally inside §7.
Affected: C-GnfsComplexity (frozen); C-NonMonicNorm (size sentences now correct; row 2 must keep
them and change only the c·N wording).
Deferred: yes — row 8 should re-examine: (a) the descent bound c_descent ≤ (8/9)^{1/3} +
1/(3(8/9)^{1/3}) (~MATHEMATICS :2895–2919) and NFS-DL text still on the old `ρ(u′)^{-1} = …(1+o(1))`
form (report-only, NFS-DL design); (b) MATHEMATICS §2 Dixon (~:1546–1558: c = 1/√2 against L_N[1/2,
1]; "L_N[1, 2c]" should be L_N[1/2, 2c]) and the ECM synthesis row (L_p[1/2, 1] vs ECM theorem
L_p[1/2, √2]); (c) `gnfs/src/polyselect/base_m.rs:22–24, 44–46` rustdoc ("balances the smoothness
bounds on both sides"; :188 `m = floor(n^{1/(d+1)})` vs the code's least m with m^{d+1} > n; row 3
touches this file's rustdoc); (d) leaked labels in `gnfs/docs/PEDAGOGY.md` ~:5133–5136, 5148, 5164,
5219, 5942–5943 ("T.G", "T.D", "C-PolyPair", "C-Score", "G.A–G.F", "G.W"); (e) the Prerequisites CEP
statement is now the large-u form with its range; confirm no consumer relied on the old Corollary;
(f) [Bai11] page numbers carried from the contract, not re-fetched.
Texture: references 18, 19, 22, 23 now carry recalled/not-read notes; [Bai11] added as "23a".

### 2 — 2026-11-11
Discovery/flex: no sieve or Murphy-score argument needs N rather than F (the sieve smooths F; the
valuation of F equals the ideal exponent only for p ∤ c, which `Relation`'s parity paragraph
already restricts to). Test gap, not a wording issue: private `alg_norm_f64` (murphy.rs) is
reached only through `score`, whose test `kat1_ordering_better_scores_higher` checks ordering, not
the value; the rustdoc says so. One outright wrong sentence fixed (gnfs PEDAGOGY §53: both sides
"norms of the same ideal (a − bα) in ℤ[α]").
Affected: none (C-NonMonicNorm propagated).
Deferred: yes — row 8 should sweep surviving "algebraic norm" used for F outside row 2's files:
`gnfs/src/sieve/{factor_base.rs:10,66; special_q.rs:34; lattice.rs:47}`, `polyselect/roots.rs:5`,
`sqrt/assembly.rs:488` (row 4 touches this file), `lib.rs:12`, `dl/ext/relation.rs:18,113,172`
(NFS-DL, report-only), gnfs tests (`special_q_kat`, `factor_base_kat`, `sqrt_*_kat`,
`factor_end_to_end_kat`, `checked_relations_kat:733`, `murphy_kat`, `dl_descent_kat`),
`shared/numfield/docs/PEDAGOGY.md:22,691`, MATHEMATICS :1906 ("homogeneous norm form"), :2514,
:2534, :2790 (NFS-DL §9, report-only); `sieve/mod.rs:195` is a runtime error string (leave).
Leaked plan labels in files row 2 edited, pre-existing: `norms.rs` module comment ("C1 resolution
in ROADMAP", "principle 4"), `line.rs` ("Principle-4 note"), gnfs PEDAGOGY ("G.B.x", "C-PolyPair",
"principle-4").
Texture: gnfs PEDAGOGY bad-prime passage (~:892–896) left for row 3.

### 3 — 2026-11-11
Discovery/flex: the entry's claims held. `discriminant` returns lc(f)·disc(f) up to sign (no division
by lc); the checked path's c·disc(F) test flags the same primes; `dedekind_factor` exposes a
maximality result (`report.index.divides_index()`), so the KAT uses it. `number_field()` callers are
tests only (`gnfs/tests/dl_relation_kat.rs`, `base_m_kat.rs`); NFS-DL relation code takes `nf` as a
parameter; `FactorBase::new` does its own root search; `sqrt/context.rs` builds its own monic form.
The old PEDAGOGY claim that the flag "cannot see p | c" for d = 2 was false and is corrected. Also
fixed in `base_m.rs` rustdoc: degree heuristic comes from minimising total cost (sides differ in
size), and m is the least integer with m^{d+1} > n (this closes item (c) from row 1's digest).
Affected: C-BadPrime (frozen).
Deferred: yes — row 8 should sweep: "principle-4" labels throughout gnfs PEDAGOGY (~:307, 504, 519,
678, 2869, 4991, …; plan vocabulary in durable text); `polyselect/mod.rs` module header "Non-monic
seam" (~:19–25) and `PolyPair` doc (~:118–120) do not name `number_field()`'s consumers;
`sieve/mod.rs` ~:243 bad-prime comment unread; `dedekind.rs` module header and
`DedekindReport::is_bad_prime` "p | disc(f)" is correct there (monic-only) — keep.
Texture: worker reflowed one over-width PEDAGOGY line and did not re-measure.

### 4 — 2026-11-11
Discovery/flex: the doc claim held: `factor_from_congruence` returns `Trivial` on
`(x − y) % n == 0 || (x + y) % n == 0` before any gcd; no code path reports a divisor there. The old
rustdoc claim "both gcds are trivial" was false for x = y = 5, N = 35 (gcd(x + y, N) = 5) and is
dropped. An existing test already covered (5, 5) but not the gcd; the new unit test asserts both
(private seam, as the plan's discovery clause anticipated, no public API).
Affected: C-TrivialCongruence (frozen).
Deferred: yes — row 8 should sweep: `docs/ROADMAP.md:253-255` (rolling doc; not edited by the chain)
still describes the old README/`assembly.rs` wording; `gnfs/docs/PEDAGOGY.md` ~:1892 ("compute the
congruence of squares") and ~:2028 ("candidate congruence of squares") unchecked in context; other
README lines over 100 characters (~:118, 120–122, 333) untouched (may be byte-count false positives).
Texture: calendar day rolled to 2026-11-12 during row 4; freshness checks rerun before row 5.

### 5 — 2026-11-12
Discovery/flex: the ECM test was strengthened, not renamed: the existing public `ecm_one_curve`
skips stage 2 when `b2 == b1`, so the same curve runs with stage 2 off and on (no new API). Two
`#[ignore]` placeholders are not real oracle checks: `gnfs/tests/lanczos_kat.rs:308` (CADO
placeholder that never invokes CADO-NFS) and `gnfs/tests/factor_end_to_end_kat.rs:858` (`todo!()`
body); development.md's "opt-in oracles are `#[ignore]`-gated" is true of them only formally, and
BENCHMARKS still lists a CADO oracle for lanczos and an ignored oracle for `factor_end_to_end_kat`
as coverage. BENCHMARKS' "21 qubits" discrepancy no longer exists (14 qubits, matching its table).
Affected: C-TestingDocs (frozen).
Deferred: yes — row 8 should sweep: the BENCHMARKS coverage lines for the two placeholder ignored
tests (claims vs reality; likely a ROADMAP/GNFS-assurance finding, not a prose fix here);
`stage2_extends_reach` doc ("stage 1 alone rarely succeeds") is not checked by its body (outside the
row's one-test scope; C-TestingDocs says every test's doc states what its body checks); outside
permitted files: `docs/NOTES.md:134` ("no proptest/quickcheck", rolling doc), `docs/SURVEY.md`
(:792 "no proptest, no quickcheck"; test counts :646–660, 730, 740, 750, 762–764, 854, 886, 905;
SURVEY is retirement-scope, report-only), `docs/ROADMAP.md:272` (stale "about 40 tests … (45
exist)"); `docs/development.md:70–74` over 100 columns (pre-existing).
Texture: Tests/Tests-added columns removed from ten BENCHMARKS tables; no measurements touched.

### 6 — 2026-11-12
Discovery/flex: `glv_fewer_steps_than_negmap` (`rho/src/ecdlp/mod.rs`) asserts only
`glv_steps <= 2 × negmap_steps` on the 1093-element `tiny_glv` group over 20 seeds, not √3 or √6;
no test asserts the Koblitz √(2m) saving (`koblitz_*` check logarithms; `binary_ecdlp_kat.rs`
asserts `steps >= 1`). The worker removed an unsourced PEDAGOGY claim ("ratio about 0.62–0.77
observed when that test was written") and the step-count figures it had copied from a bench
comment, pointing at BENCHMARKS instead. BENCHMARKS records the step-count anomaly for
`solve_dp_negmap` only, not for batch; the passages claim no more.
Affected: C-HeuristicRatio (frozen).
Deferred: yes — row 8 should sweep unlabelled ratios outside row 6's files:
`rho/src/ecdlp/cli.rs:96` ("Expected to reduce wall time by ~√2 vs plain DP", unlabelled);
`cli.rs:45` and `rho/src/binary_ecdlp/koblitz.rs:351` are already labelled. Two README table rows
(~:338 and the negation-map row) exceed 100 columns (cannot wrap a table row).
Texture: ratios in MATHEMATICS now sit only in "Orbit-collapsing maps: negation and GLV"; the old
~:1006–1011 passage states none.

### 7 — 2026-11-12
Discovery/flex: every source was obtainable, so no abstract-only fallback. Corrections: the old
§11.4.4 credited PZ03 with an O((log r)²) cost for "the two transforms" that PZ03 does not state
(it is the chapter's own O(t²) result); PZ03's gate/qubit counts are logical, whole-algorithm,
prime-field, and exclude exponent registers (the project's 1280/4100 layout counts include them, so
they are not comparable); "thousands of fault-tolerant qubits" removed as unscoped. Cain v1's
abstract says "a few days" for 26,000 qubits while its body and conclusion say 10 days; the doc
follows the body and notes the mismatch; the 9,739-qubit headline pairs with ≈1000 days at a 1 ms
cycle. Babbush's zero-knowledge proof certifies only the secp256k1 point-addition circuits; full
counts, generic-n claim, physical qubits and runtimes are composed/uncertified, and the doc says so.
Schrottenloher's Babbush rows are Babbush's stated figures (circuits not inspectable) and its
simulation is word-level; presented as a comparison with a stated target. No cross-source ratios.
Source-version caveats: PZ03 read as arXiv v2 (the QIC typeset version was not read); Roetteler as
arXiv v3 (not the ASIACRYPT proceedings); Cain appendices unread; Schrottenloher's arithmetic
sections not checked line by line.
Affected: C-QuantumEvidence (frozen). No new algorithm requirement.
Deferred: yes — row 8 should sweep: `docs/SURVEY.md` (~:2139–2153, 2189, 2465 "768 qubits"
attributed to PZ03, which states ≈1500 (1800) logical at n = 256; ~:2168–2173; ~:2180–2181 "millions
of physical qubits" for logical-qubit papers; SURVEY is report-only); `shor/docs/PEDAGOGY.md` §6.2
~:535 attributes "Brent–Kung adder trees" to Roetteler et al., not found in their paper (outside
the permitted sections of row 7); `docs/MATHEMATICS.md` ~:4003 second PZ03 entry (§10 list) lacks
the arXiv-version note; Chevignard/Litinski figures are labelled secondhand from Babbush §II.B.
Texture: bibliographic details of Shor SIAM (26(5), 1484–1509) and Roetteler (pp. 241–270) checked
against Crossref; §11.7.1 added (tables + evidence-boundary + sources-read note).

### 8 — 2026-11-12
Discovery/flex: row 8 ran as a sweep with allowed extras (list in the ledger row). Items for the
boundary fork to weigh:
- MATHEMATICS §2 Dixon: the proof gives L_N[1/2, √2] (not 1); theorem, §1 intro, ECM synthesis row
  (L_p[1/2, √2], smallest prime) and `docs/PEDAGOGY.md:75` corrected; Pohlig–Hellman
  L_p[1/2, 1] → L_p[1, 1/2] = √p. The worker added a trial-division variant (L_N[1/2, 2] at
  c = 1/2) from its own short derivation, **not sourced**.
- `PolyPair::skew` and `factor_base_bounds` are never set by any code (`score` reads `skew`); docs
  claiming "set by Murphy-E scoring / sieving" were false and a described "factor-base balancing"
  step does not exist; docs now say "not implemented". A code-state finding, no behaviour change.
- numfield PEDAGOGY §10 falsely said sieving evaluates `NumberFieldElement::norm()` (never called in
  `gnfs/src`); corrected. "x ≡ ±y mod p and mod q independent events" (gnfs PEDAGOGY §46, §57)
  rewritten as independent sign choices, labelled a modelling assumption.
- Ignored tests: `todo!()`/`unimplemented!()` placeholders (`factor_end_to_end_kat`
  `kat_c_oracle_80_100_bit_challenge`, `dl_individual_log_kat`, `dl_end_to_end_kat`, `merge_kat`,
  `ghs_kat`); `line_sieve_kat::kat_c_cado_nfs_oracle` compares the sieve with itself;
  `lanczos_kat::kat_c_cado_oracle_n35` never invokes CADO-NFS; `semaev_kat`, `ssa_kat` (2) and padic
  `log_kat` assert hard-coded PARI values without calling `gp`. Only three files start `gp`; no
  test starts CADO-NFS. BENCHMARKS (five rows) and `development.md` now say so; no test changed. This
  is a GNFS/oracle-assurance finding for the ROADMAP.
- `rho/src/ecdlp/cli.rs:96` doc comment is clap `--help` text: the help output changed (expected-
  not-measured label); no solver behaviour. MATHEMATICS reference lists: §11.9 holds 52–55 and the
  chapter-12 list was renumbered 52–62 → 56–66 so the shor PEDAGOGY cite resolves uniquely.
- gnfs PEDAGOGY: ~430 denylist hits (G.x/D.x/E.x labels, `C-*`, `principle-4`, `T.G`) translated in
  §§1–62, §58/§59 rewritten (contract table, principle verdicts); §63+ (NFS-DL) untouched.
Affected: none (all seven contracts upheld; no wrong `Ok`; no new dependency edge).
Deferred: yes — boundary fork should weigh the unsourced Dixon trial-division line and the
chapter-12 renumbering. Report-only (NFS-DL, rolling and retirement-scope files, not edited):
MATHEMATICS §9: "algebraic norm" for F ~:2529, 2549, 2805; old `ρ(u′)^{-1}` / `L_p[1/3,
1/(3c′)]·(1+o(1))` form ~:2787–2789, 3010–3026, 3057–3076, 3101 and the `c_descent` bound ~:3047,
`(1+o(1))` ~:2796; `gnfs/src/dl/**`: `ext/relation.rs:18,113,172` "algebraic norm", labels
`schirokauer.rs:35`, `mod.rs:116,118` (cites a "C-Relation contract" no longer in `sieve/mod.rs`),
`linalg/blockvec_fl.rs:33`, `linalg/wiedemann_fl.rs:24,61`, `ext/descent.rs:35,252`; tests
`dl_descent_kat.rs:35,416,490,549` ("the PLAN"), `schirokauer_kat.rs:11,162,164,231`; gnfs PEDAGOGY
§§63–71 and "Further Reading (D.W)" (~100 lines of leaked labels, ~:5270); ROADMAP :253–255, 270,
272; NOTES :134; SURVEY :792, test counts ~:646–905, "768 qubits" ~:2139–2153, 2189, 2465,
"millions of physical qubits" ~:2174–2181. Needs-human: `special_q_kat.rs`
`kat_b_yield_comparison_principle4_annotated` (identifier); `coppersmith_kat.rs`, `root_sieve_kat.rs`
carry `principle-4` strings (some in assertion messages; not enrolled, unedited); `line.rs` "~99%"
and the PEDAGOGY "Montgomery 2–4×, Brent 24%, batched GCD 2–5×", Cavallar 30–50% figures are
unlabelled but outside the heuristic-ratio contract; `docs/PEDAGOGY.md:68`, MATHEMATICS:268 mix
chapter numbers with GNFS-internal §7 (not changed). Over-width prose left (pre-existing, regex
`.{101,}`): gnfs PEDAGOGY 19 outside tables, numfield PEDAGOGY 5, shor PEDAGOGY 5, BENCHMARKS 1,
MATHEMATICS 1; tables unwrappable.
Texture: committer noted its own model is Haiku 5.5 (mechanical, not credited); trailers unchanged.

### ◆ boundary (rows 1–8) — 2026-11-12
Discovery/flex: fork verdict `still-on-intent`; §7 faithful to C-GnfsComplexity. The review found
incomplete propagation of C-TestingDocs (ignored placeholders and recorded-value tests described
as external oracles, with misleading `#[ignore]` reasons) and residual planning vocabulary
("reshard", "cross-track"); both corrected in `154755c` (a medium deviation: `#[ignore]` reason
strings, not only comments, were edited). Two NFS-DL sentences were edited for terminology in row 8
(gnfs PEDAGOGY ~:5459–5460, ~:5779–5782), no capability claim changed. gnfs PEDAGOGY quotes of
source comments are now paraphrases because the sources still carry `principle-4` labels.
Affected: none.
Deferred: yes — the review is incomplete (scope limit above); row 9 (anneal) should take the
unenrolled `principle-4` source comments (`linalg/{blockvec,operator,wiedemann}.rs`,
`filter/matrix.rs`, `lattice_kat`, `coppersmith_kat`, `root_sieve_kat`) only if they are in the arc's
files (they are not; report) and the line-length pass over the arc's files; GNFS/oracle assurance
(no test invokes `gp` outside three files, none invokes CADO-NFS; several ignored tests are
placeholders) and the DP step-count anomaly remain ROADMAP findings.
Texture: coverage 94.46% lines at the boundary on tree `154755c`.

## Discoveries & risks

Carried down from ROADMAP and the sharding survey; phrased for adjudication.

- **Row 1:** if no readable source and no self-contained derivation reaches (64/9)^{1/3} under the
  norm-size model the chapter uses: juncture `design-uncertain-HALT`. If the corrected result changes
  §9.7's restatement: internal-continue, restatement only.
- **Row 2:** if a sieve or Murphy-score argument actually requires N rather than F (or vice versa)
  and the code computes the other: additive-reshard (surface; behaviour is out of scope).
- **Row 3:** if the factor base discards maximal bad primes it could use: a coverage limit, not a
  defect — record for the GNFS assurance scope (internal-continue).
- **Row 4:** if no seam exposes the N = 35 classification without new public API: add a private
  unit test in `assembly.rs` (internal-continue); no new public API.
- **Row 5:** if the ECM parameter search fails inside its bound: rename and re-doc the test
  (internal-continue).
- **Row 7:** a full text not obtainable → label abstract-only and drop quantitative comparison
  (internal-continue). A paper implying a new algorithm requirement: surface, never implement.
- **Additive-reshard (surface, do not absorb):** NFS-DL capability/equation claims (MATHEMATICS
  §9.6, 2710–2731, 3007, 3031, 4733, 5092; gnfs PEDAGOGY 4880–4905, 5367–5370, 5701–5991;
  `gnfs/src/dl/**` rustdoc); SURVEY's 768-qubit and wrong-journal text; README navigation and track
  vocabulary; the DP step-count anomaly.
- **Destructive-HALT:** a corrected claim exposes a wrong `Ok` on an existing fixture; a row needs a
  new crate or dependency edge.

## Notes for executors

**Sizing.**
- (1) Ambient complexity is high in the docs (MATHEMATICS ≈ 5,200 lines, gnfs PEDAGOGY ≈ 6,100) but
  low in code.
- (2) The floor is row 1: one derivation must be internally consistent, so §7 and the three
  restated smoothness estimates are one unit.
- (3) Cost of a design error is highest in row 1 (a re-published wrong derivation), which drove its
  T0 inflection juncture.
- (4) Correctness-criticality is high for rows 1–4 (mathematical claims); rows 5–6 are moderate.
- (5) Inner-loop bandwidth for prose is low (no test catches a wrong sentence), which is why rows
  are split by claim family rather than by file.

Outlier check: rows 2 and 3 share files and notation; they are kept separate because 3 freezes a
test-enforced classification and 2 is notation-only (split lever: 3, 5). Rows 4 and 6 are small
but each closes one claim family with its own pin. Row 8 is the only cross-cutting row. No row uses
F/C: rows are prose-dominated, and the coverage floor is workspace-wide.

**Routing.** T1 implements every row. T0 runs the row-1 inflection juncture and the ◆ boundary
fork. Row 1's juncture writes C-GnfsComplexity before row 1 is dispatched.

**Register.** State the property or invariant, never the plan coordinate. Correct claims in place;
no new chapters; keep section numbers and anchors stable unless a row says otherwise.

**Anneal denylist** (grep enrolled durable files at ◆, excluding PLAN/ROADMAP/NOTES/SURVEY/AGENTS):

    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat(egory)?[ -]?[ABCI]\b
    \bC-[A-Z][A-Za-z0-9]+\b
    \b[A-H]\.[A-Z0-9]+\b|(?i)\bprinciples?[ -]?[0-9]+\b
    \bTrack [A-Z]\b|\bPhase [0-9]+\b|(?i)\b(phase|session|stage)[ -][α-ωΑ-Ω]\b|[α-ω]\.[0-9]+
    (PLAN|ROADMAP|NOTES|SURVEY)(\.md)?\b|plan-(run|shard|juncture)|(?i)\b(the|per the) plan('s|’s)?\b

KAT, algorithmic step numbers, group notation (S₃), bibliography keys ([S09], [LLMP93]), E₁/E₂,
section signs (§7) and quantum phase terminology are domain terms: classify collisions, do not
rename. `\bS[0-9]+\b` is known to collide with bibliography keys and Jacobian variables.

**First run:** after committing the planning docs, from `@plan-admin`:

    /plan/run docs/PLAN.md halt-at-boundaries
