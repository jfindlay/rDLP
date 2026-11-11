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

Prose-enforced. To be frozen at 1 by the inflection juncture, which writes here: the step list
(norm-size model, u, the one smoothness estimate, cost terms, optimisation), the reconciled
constants, the heuristic-status statement, and which reference each step follows (read or
recalled). Consumed by 1, 8. The `§7` heading and its position stay stable.

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

Gate evidence: *(none yet)*

Pending checks: `cargo dep-policy` and `cargo audit` due on the first active run (new calendar day
since the 2026-11-10 records).

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Re-derive the GNFS L[1/3, (64/9)^{1/3}] complexity from one smoothness estimate | pending | | |
| 2 | State the non-monic homogeneous value as c·N(a − bα) and localise the index | pending | | |
| 3 | Separate bad primes from index divisors and name `number_field()`'s consumers | pending | | |
| 4 | Pin the trivial-congruence policy and say linear algebra yields candidates | pending | | |
| 5 | Make the testing documentation describe the tests that exist | pending | | |
| 6 | Label the GLV, negation-map and Koblitz ratios as expected, not measured | pending | | |
| 7 | Rework the quantum-resource exposition from primary sources | pending | | |
| 8 | Correct surviving claims the exposition rows contradicted ◆ | pending | | |
| 9 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

*(none yet)*

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
