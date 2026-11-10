# rDLP — Solver-outcome follow-ups

Execution graph for the ROADMAP scope "Solver-outcome audit", follow-up items. The previous
(completed) solver-outcome conformance plan is preserved in git at `fc434fb:docs/PLAN.md`.

## Purpose (design intent)

Close the follow-ups the solver-outcome repairs found but were not allowed to change, under
`docs/development.md`'s "Solver outcomes and fallbacks": SSA must recover k from a lift on which
the formal-log quotient is exact (no empirical halving), and must type bad input and invariant
breaches distinctly; the GHS transfer must identify GF(2^l) inside GF(2^m) by a field
isomorphism that holds beyond the fixture; MOV must reject fields it cannot represent instead of
truncating them in release builds; ECM stage 1 must actually be B1-power-smooth at every B1;
factoring solvers must surface a broken divisor invariant instead of retrying past it; rho
factoring must certify every prime its domain admits; and the `rho-dlog` CLI and ECDLP benches
must run instances that can finish. Each repair corrects the exposition it directly falsifies.

Scope boundary: `rho/src/{ssa,ghs,pairing,factor}`, `rho/src/ecdlp/cli.rs`, `rho/src/curve`,
`rho/benches/ecdlp.rs`, `shared/numth/src/{ecm,prime}.rs`, their tests and exposition. Out of
scope: NFS-DL and MOV's k=2 log engine, Shor factoring's register width, removing ECM,
re-measuring published benchmark numbers (tooling scope), new algorithms, new crate or dependency
edges (dev-dependencies included), git-history housekeeping.

## Verify gate

- VERIFY_TYPES: `cargo check --workspace`
- VERIFY_TEST: `cargo test --workspace`
- Cadence (authority: `docs/development.md`, "Verification cadence"): VERIFY_TYPES and
  VERIFY_TEST once on the final inputs of each completed code row; `cargo coverage` (80% line
  floor) at the ◆ boundary and at final delivery, earlier after substantial test removal (rows 2
  and 7 delete tests and branches — the driver judges), deferred with that checkpoint named at
  interior rows; `cargo dep-policy` and `cargo audit` when dependency inputs change, on the first
  active run of each day, and before release. A row with a deferred check is complete under its
  per-unit gate, never a full-gate pass.
- Not gates: `cargo fmt-check`, `cargo lint` (unadopted; do not fix drive-by).
- Oracles (manual, only if `gp` is on `PATH`; never a gate): `kat_ssa_pari_cross_check`,
  `kat_lift_pari_cross_check` at rows 1–2; any `kat_pari_*` row 3 adds.
- Bench smoke (row 9 only, not a test substitute): `cargo bench -p rho --bench ecdlp -- --test`.

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 @architect | Type SSA input and invariant failures | A | T1 (juncture T0) | C-OutcomeShape | `rho/src/ssa/{mod.rs,reduce.rs,formal_log.rs,lift.rs}`, `rho/tests/ssa_kat.rs` |
| 2 | Lift SSA curves non-canonically and drop the halved candidate | B | T1 | C-SsaOutcome | `rho/src/ssa/{reduce.rs,lift.rs,mod.rs}`, `rho/tests/ssa_kat.rs`, `rho/benches/attacks.rs`, `docs/MATHEMATICS.md` §10.2, `docs/PEDAGOGY.md` §11 |
| 3 | Identify the GHS conorm subfield with `restrict` | B | T1 | C-OutcomeShape | `rho/src/ghs/{transfer.rs,curve.rs,mod.rs}`, `rho/tests/ghs_kat.rs`, `docs/MATHEMATICS.md` §10.3, `docs/PEDAGOGY.md` §12 |
| 4 | Guard MOV's single-limb field domain in release builds | B | T1 | C-OutcomeShape | `rho/src/pairing/{mov.rs,tate.rs}`, `rho/tests/mov_kat.rs`, `docs/MATHEMATICS.md` §10.5, `docs/PEDAGOGY.md` §10 |
| 5 | Apply ECM stage-1 prime powers one ladder call at a time | B | T1 | — | `shared/numth/src/ecm.rs`, `shared/numth/docs/PEDAGOGY.md` §5, `docs/MATHEMATICS.md` (ECM) |
| 6 | Report non-dividing factor candidates as `Inconsistent` | B | T1 | C-OutcomeShape, C-EcmStage1 | `rho/src/factor/{rho.rs,mod.rs}`, `rho/tests/factor_kat.rs`, `shared/numth/src/ecm.rs` |
| 7 | Certify rho-factoring primes below 2^80 with the 13-base bound | B | T1 | C-FactorInconsistent | `rho/src/factor/rho.rs`, `rho/tests/factor_kat.rs`, `shared/numth/src/prime.rs` (rustdoc), `docs/MATHEMATICS.md` (primality) |
| 8 | Project the generic rho-dlog curve onto its prime-order subgroup | B | T1 | C-DlogOutcome (prior arc) | `rho/src/curve/generic.rs`, `rho/src/ecdlp/cli.rs`, `rho/tests/ecdlp_kat.rs`, `README.md`, `docs/PEDAGOGY.md` §6 |
| 9 | Retarget the ECDLP benches to groups that finish | C | T1 | C-DpBudget (prior arc), C-GenericSubgroup | `rho/benches/ecdlp.rs`, `rho/src/curve/test_curves.rs`, `rho/tests/ecdlp_kat.rs`, `docs/BENCHMARKS.md` |
| 10 ◆ | Correct the follow-up exposition | I | T1 (◆ fork T0) | all | `docs/MATHEMATICS.md`, `docs/PEDAGOGY.md`, `shared/numth/docs/PEDAGOGY.md`, `docs/BENCHMARKS.md`, `README.md`, `rho/src/pairing/{mov.rs,miller.rs}` (rustdoc/comments) |
| 11 | Anneal cleanup (skippable) | I | — | none | durable files touched by rows 1–10 |

Order: 1→2 (SSA); 5→6→7 (`ecm.rs`, then `rho.rs`); 8→9 (bench may use the projected curve);
3 and 4 are independent; row order is the run order. Row 1 is the inflection point: the juncture
fork designs C-SsaOutcome against `rho/src/ssa/*` and the row-2 lift mechanism before row 1
implements it.

## Session detail

### 1 @architect — Type SSA input and invariant failures

Deliverable per C-SsaOutcome (frozen here by the juncture). `ssa_solve` validates input first
(G ≠ ∞, G and Q on the curve, `n` used: the anomalous order equals p, n·G = ∞, n·Q = ∞) and
returns `InvalidInput(SsaInputDefect)`; invariant breaches after the curve is proven anomalous —
`formal_log.rs`'s v_p(t) < 1 branches (~:92–103, :127–131), `fp_divide_logs`'s non-unit g-log —
become the dedicated invariant-failure class, not `Padic(..)`. `LiftFailed`/`Padic` survive only
where the juncture shows they are reachable on validated input. The halving stays until row 2;
do not change the lift here. `Display`/`source` updated.

KATs: each `SsaInputDefect` reachable (off-curve G and Q, G = ∞, wrong n); the invariant class is
reachable through a unit seam (like `resolve_candidate_*`); existing solve KATs unchanged.

Subtleties: the juncture also fixes the outcome of a second-lift failure (row 2): if two distinct
non-canonical shifts provably cannot both land in E₂, failure is the invariant class; otherwise a
typed exhaustion outcome. Pin it in the contract now so row 2 does not reopen the enum.

### 2 — Lift SSA curves non-canonically and drop the halved candidate

Deliverable per C-SsaLift. Compute one lifted curve a ↦ a + p·s (s ≢ 0 mod p) and pass it to both
`lift_point` and `multiply_by_p_proj` (today each reads `curve.a` itself — the single seam); k is
the raw quotient ψ(p·Q̃)/ψ(p·G̃) mod p on a lift with p·G̃ ∈ E₁ \ E₂; retry one second shift if the
first lands in E₂, per the contract. Delete `Candidate`, `resolve_candidate`'s halving and their
unit tests; re-check whether precision k = 8 (`reduce.rs:132`) is still needed and name the
constant. Read a primary source for the canonical-lift failure (Smart 1999 and/or Leprévost–
Monnerat–Varrette–Vaudenay 2005) and cite it, or label the claim unverified. Correct MATHEMATICS
§10.2 (log_Ê applied to G̃ ∉ E₁; "p·G̃ = O in E(Z/p²)"; the halving paragraph; O(log p) vs the
O(p) multiply loop) and PEDAGOGY §11 (names and signatures).

KATs: the stored a = 0 lift of `anomalous_toy` lands in E₂ (valuation 2) via a test seam; every k
on `anomalous_toy` recovered from the raw quotient on the shifted lift; a second anomalous fixture
with a ≠ 0 found by a bounded, recorded search (naive point count, small p; `verify_anomalous`
re-proves it in-tree), every k recovered; `Inconsistent` now invariant-only.

Subtleties: `lift_point` is `pub` and the KATs call it; a signature change touches `ssa_kat.rs`.
`bench_ssa` asserts `Ok(3)` before timing — keep it true.

### 3 — Identify the GHS conorm subfield with `restrict`

Deliverable per C-GhsSubfieldIdentification. Replace the raw-bits reinterpretation
(`restrict_to_subfield`, `transfer.rs:188–196`) with `shared_gf2m::restrict`, and apply the same
identification to `extract_ghs_curve`'s coefficients (`curve.rs:182–184`) so curve and transfer
agree; keep `descent.rs:554`'s basis-coordinate extraction consistent. Reject parameters where
`poly_l` is not the minimal polynomial of β = α^((2^m−1)/(2^l−1)), or `l ∤ m`, or m > 63
(`restrict` panics there) with a typed `GhsError` before any transfer.

KATs: re-derive `transfer_sample_point_known_answer` and `transfer_g_plus_p_known_answer`
independently of the code (hand basis expansion written in the test doc comment, or an opt-in
`kat_pari_*`); the expected 2↔3 swap is a hypothesis, not the oracle. Homomorphism and
log-preservation KATs unchanged; mismatched `poly_l` → the new error.

### 4 — Guard MOV's single-limb field domain in release builds

Deliverable per C-MovDomain. `MovError::Unsupported(&'static str)`; `mov_reduce` rejects p ≥ 2^64
(any non-zero high word of `curve.p`) before `check_mov_input`'s arithmetic; every single-limb
conversion on the path — `check_mov_input` (~:216), `tate.rs:132`, and the final exponent
`exp_u128 as u64` at `tate.rs:147` — is guarded upstream or made exact. The `debug_assert`s may
remain as internal invariants.

KATs: a curve with p ≥ 2^64 → `Unsupported` (asserted before any debug assertion fires); a k = 2
case whose final exponent (p²−1)/ℓ exceeds 2^64, if a fixture exists or is cheaply constructed,
verifies; existing MOV KATs unchanged.

Subtleties: `tate.rs:147` may be a live wrong-pairing path on supported k = 2 inputs (see
Discoveries); every `Ok` is verified, so it would show as `DegeneratePairing`/`Inconsistent`, not
a wrong k. Record what the fixtures show.

### 5 — Apply ECM stage-1 prime powers one ladder call at a time

Deliverable per C-EcmStage1: `run_stage1` applies q = p^⌊log_p B1⌋ for each prime p ≤ B1 with its
own `montgomery_ladder` call (it already takes `&Uint<4>`); `stage1_scalar` and its wrapping
product go; optional early gcd per step. Correct the stage-1 limit text in numth PEDAGOGY §5
(~:336–344) and MATHEMATICS (~:1286–1295), and the Montgomery-ladder comment beside it.

KATs: at B1 ≤ 177 the per-prime point equals the old single-product ladder (projectively); a
fixture that factors only when a prime power above the old overflow point is applied (bounded,
recorded search for (n, sigma, B1)); existing ECM KATs unchanged. No new dependency for a big
reference product.

### 6 — Report non-dividing factor candidates as `Inconsistent`

Deliverable per C-FactorInconsistent: `FactorError::Inconsistent { d: u128 }` and
`EcmError::Inconsistent { d: Uint<4> }` (both stay `Copy`); d ∈ {1, n} remains a failed attempt; any
other d failing 1 < d < n, n mod d == 0 returns `Inconsistent` at once, never retried.

KATs: each `Inconsistent` reachable through a private seam (the gcd is correct on real inputs);
d = n still consumes an attempt; existing outcomes unchanged.

### 7 — Certify rho-factoring primes below 2^80 with the 13-base bound

Deliverable per C-RhoPrimeBound: read Sorenson and Webster, "Strong pseudoprimes to twelve prime
bases" (Math. Comp. 2017; arXiv:1509.00864) and confirm ψ₁₂ and ψ₁₃; record the citation in
`prime.rs` rustdoc (replace "recalled, not re-read"). Then every odd n < `MAX_SUPPORTED` with
`is_prime` → `Prime`, justified by `shared_numth::MILLER_RABIN_PROOF_BOUND` with a compile-time
assertion `MAX_SUPPORTED ≤ MILLER_RABIN_PROOF_BOUND`; remove `PRIME_PROOF_BOUND` and the
`Unsupported`-prime branch; `Unsupported` only for n ≥ 2^80. Correct MATHEMATICS's primality
passage (~:1095) if it states the rho range.

KATs: 2^79 − 67 → `Prime` (was `Unsupported`); ψ₁₂ (composite) factors or is not `Prime`; n ≥ 2^80
→ `Unsupported`.

Subtleties: if the paper cannot be read, do not raise the bound (see Discoveries).

### 8 — Project the generic rho-dlog curve onto its prime-order subgroup

Deliverable per C-GenericSubgroup: the `generic` generator becomes G′ = 119·G with order
n′ = 77_507_328_043_853_087, stored with the curve like `secp_k1_toy`; `GENERIC_N` (composite) and
its "Computed by Sage" provenance go. README and PEDAGOGY §6 describe `generic` as a prime-order
subgroup and state which default runs can exhaust their budget (secp-toy Brent: per-attempt cap
2²⁹ − 1 vs ≈ 2.7·10⁹ expected steps, arithmetic only).

KATs: G′ ≠ ∞, n′·G′ = ∞, `is_prime(n′)` (together they prove ord G′ = n′ without the external
order); Brent's input check accepts (G′, n′).

### 9 — Retarget the ECDLP benches to groups that finish

Deliverable per C-BenchFixture: every bench in `rho/benches/ecdlp.rs` asserts its solve succeeds
once before timing (as `bench_ssa` does) on a group sized so one iteration finishes in seconds in
release; reuse an existing prime-order fixture if one fits, else add one (prime order proved
in-tree; generator and order found by a bounded, recorded search — an external tool may find,
never certify). Budgets per C-DpBudget. BENCHMARKS: replace the ecdlp-bench setup description;
mark old numbers as not re-measured (re-measurement belongs to the tooling scope).

KATs: the fixture's prime-order test; `cargo bench -p rho --bench ecdlp -- --test` completes
(smoke, recorded; not a test substitute).

### 10 ◆ — Correct the follow-up exposition

Deliverable: (a) MOV R = G: `DegeneratePairing`'s rustdoc, the `vertical_at` comment
(`miller.rs:235–249`) and MATHEMATICS §10.5 / PEDAGOGY §10 say the evaluation point lies on the
Miller divisor's support (the value is undefined; the tangent-line numerator vanishes), reported
as `DegeneratePairing`. (b) Sweep superseded claims across all durable surfaces these solvers
touch, including unchanged passages: the SSA halving and "`LiftFailed` means bad input"; GHS
"valid only on the fixture"; MOV `debug_assert`-only domain; ECM "B1 ≳ 177 wraps"; `generic`
composite order; rho factoring "probable prime in [ψ₁₂, 2^80) is `Unsupported`"; `pohlig.rs`'s
ECM mention. (c) Anneal the arc's durable files against the denylist.

KAT substitute (prose row): every changed behaviour claim cites the test that exercises it or is
labelled; flagged as this row's contract.

### 11 — Anneal cleanup (skippable)

`/plan-anneal` over the arc's durable files: translate leaked plan coordinates and rewrap
ledgered line-length overruns. The driver may mark it `skipped` if the ◆ anneal and its fork ran
clean and no overruns were ledgered.

## Cross-session contracts

### C-OutcomeShape

Prose + test-enforced. Carried from the prior arc unchanged. Consumed by 1–9. Every solver
follows `docs/development.md` rules 1–6: success only after an exact independent check; typed
errors separating invalid input, unsupported domain, resource exhaustion, invariant failure
(surfaced, never retried) and proven negatives. No crate edges to share types.

### C-SsaOutcome

Compiler + test-enforced. **To be frozen at 1 (inflection juncture designs it).** Consumed by 2,
10. Draft for the juncture:

```rust
pub enum SsaError {
    InvalidInput(SsaInputDefect),
    NotAnomalous,                    // exact count only
    Unsupported(&'static str),
    Inconsistent { k: u64 },         // valuation-1 lift, k fails k·G == Q (invariant class)
    InvariantFailure(&'static str),  // breached formal-log / divisor-log invariants
    // juncture: second-lift exhaustion variant, if not provably impossible
    // juncture: whether LiftFailed / Padic remain reachable after validation
}
pub enum SsaInputDefect { GeneratorAtInfinity, GeneratorNotOnCurve, TargetNotOnCurve,
                          OrderNotFieldPrime, WrongGeneratorOrder, TargetOutsideSubgroup }
```

### C-SsaLift

Test + prose-enforced. To be frozen at 2. Consumed by 10. `ssa_solve` lifts the curve once
(a ↦ a + p·s) and uses that one lifted curve for point lifting and multiplication by p; k is the
raw quotient on a lift with p·G̃ ∈ E₁ \ E₂; no halved candidate exists.

### C-GhsSubfieldIdentification

Compiler + test-enforced. To be frozen at 3. Consumed by 10. Transfer and curve extraction
identify GF(2^l) ⊂ GF(2^m) through `shared_gf2m::restrict`; parameters outside `restrict`'s domain
or with `poly_l` ≠ minpoly(β) are a typed `GhsError`; the transfer constants have provenance
independent of the code.

### C-MovDomain

Compiler + test-enforced. To be frozen at 4. Consumed by 10. `MovError` gains
`Unsupported(&'static str)`; p ≥ 2^64 is rejected before arithmetic; no silent single-limb
truncation remains on the `mov_reduce` path.

### C-EcmStage1

Test + prose-enforced. To be frozen at 5. Consumed by 6, 10. Stage 1 applies each prime power
≤ B1 separately; no product scalar is formed; correct at every B1.

### C-FactorInconsistent

Compiler + test-enforced. To be frozen at 6. Consumed by 7, 10. `FactorError::Inconsistent { d }`,
`EcmError::Inconsistent { d }`; trivial d is a failed attempt, any other failed divisor check is
surfaced at once.

### C-RhoPrimeBound

Compiler + test-enforced. To be frozen at 7. Consumed by 10. Odd probable primes below
`MAX_SUPPORTED` are `Prime`, justified by a compile-time `MAX_SUPPORTED ≤
MILLER_RABIN_PROOF_BOUND`; `PRIME_PROOF_BOUND` removed.

### C-GenericSubgroup

Test-enforced. To be frozen at 8. Consumed by 9, 10. The `generic` curve's generator has proven
prime order n′ = 77_507_328_043_853_087 (checked in-tree).

### C-BenchFixture

Test + prose-enforced. To be frozen at 9. Consumed by 10. Every ECDLP bench solves its instance
before timing; any fixture it adds has in-tree prime-order proof and recorded provenance.

## Progress ledger

Bindings: PLAN = `docs/PLAN.md`; VERIFY_TYPES = `cargo check --workspace`; VERIFY_TEST =
`cargo test --workspace`; boundary check = `cargo coverage`; freshness checks = `cargo
dep-policy`, `cargo audit` (separate commands, never combined). Cadence per "Verify gate".

Gate evidence (one line per check: check, outcome passed/reused/deferred, command and working
directory, exit status, time, log reference, input identity, environment; deferred → next
checkpoint): *(none yet)*

Pending checks: *(none)*

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Type SSA input and invariant failures | pending | | |
| 2 | Lift SSA curves non-canonically and drop the halved candidate | pending | | |
| 3 | Identify the GHS conorm subfield with `restrict` | pending | | |
| 4 | Guard MOV's single-limb field domain in release builds | pending | | |
| 5 | Apply ECM stage-1 prime powers one ladder call at a time | pending | | |
| 6 | Report non-dividing factor candidates as `Inconsistent` | pending | | |
| 7 | Certify rho-factoring primes below 2^80 with the 13-base bound | pending | | |
| 8 | Project the generic rho-dlog curve onto its prime-order subgroup | pending | | |
| 9 | Retarget the ECDLP benches to groups that finish | pending | | |
| 10 | Correct the follow-up exposition ◆ | pending | | |
| 11 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

*(none yet)*

## Discoveries & risks

Carried down from ROADMAP and the sharding survey; phrased for adjudication.

- **Inferred (survey), enrolled 4:** `tate.rs:147` truncates the final exponent (p^k−1)/ℓ to u64;
  for k = 2 it exceeds 2^64 once p ≳ 2^32·√ℓ, so supported inputs may get a wrong pairing value.
  Verification catches any wrong k. Internal-continue; record whether a fixture was affected.
- **Confirmed (survey), enrolled 3:** the two GHS transfer constants have no independent
  provenance; `extract_ghs_curve` uses the same raw-bits reinterpretation. If an independent
  derivation disagrees with both the raw and the `restrict` outputs: **destructive-HALT** (the
  transfer itself is suspect).
- **Inferred (ROADMAP, recalled literature), enrolled 2:** a non-canonical lift gives valuation 1.
  If shifted lifts still land in E₂ generically: **destructive-HALT** (the repair's premise is
  wrong). If the a ≠ 0 fixture search fails inside its bound: internal-continue, record, narrow
  the KAT.
- **Recalled (ROADMAP), enrolled 7:** Sorenson–Webster constants. If the paper cannot be read or
  disagrees: additive-reshard — keep ψ₁₂ behaviour, surface; never raise the bound on a recalled
  constant.
- **Confirmed (survey), enrolled 8:** `GENERIC_N` = 7·17·77_507_328_043_853_087 ("Computed by
  Sage"). If n′·G′ ≠ ∞: additive-reshard (surface; dropping `generic` is the user's fallback
  decision).
- **Inferred, enrolled 5:** Z ≡ 0 modulo a factor persists through later ladder calls, so the
  final gcd still detects it; if stage 2 depends on stage 1's point semantics, internal-continue.
- **Internal-continue:** call-site adaptation outside a row's files when the compiler forces it;
  bounded recorded fixture searches; test-seam additions for unreachable invariant variants.
- **Additive-reshard (surface, do not absorb):** NFS-DL `solve_dl` and MOV's k=2 engine; Shor
  register width; ECM removal; benchmark re-measurement; plan-label survivors outside the rows'
  files (project-wide anneal scope); `@git-editor` trailer fix and the stale stash.
- **Destructive-HALT:** an `Ok` returns a wrong answer on an existing fixture that the module
  cannot fix; a contract needs a new crate or dependency edge.
- **Preflight:** commit the working-tree cadence edits (`AGENTS.md`, `docs/development.md`) and
  the new PLAN/ROADMAP before `/plan-run`.

## Notes for executors

**Sizing.**
- (1) Ambient complexity is moderate: SSA and GHS are small modules with intricate mathematics.
- (2) The floor is row 2: the shared lift, `Candidate` removal, fixture search and §10.2
  correction are one unit — not split.
- (3) Cost of a design error is highest for the `SsaError` shape (consumed by row 2), which drove
  the 1|2 split at the enum freeze and the inflection juncture on row 1.
- (4) Correctness-criticality is high, as in the prior arc.
- (5) Inner-loop bandwidth is high: rho and numth tests are fast; no row touches Shor circuits.

Outlier check: rows 6 and 7 are small and share `rho.rs`; kept separate because 7 carries its own
literature read and behaviour flip (user decision). Rows 8 and 9 kept separate (different
surfaces; 9's fixture search is its own experiment). Row 10 is the only cross-module row and is
prose. No row uses F/C: rows are small and the coverage floor is workspace-wide.

**Routing.** T1 implements every row; T0 runs the row-1 inflection juncture and the ◆ boundary
fork.

**Register.** State the property or invariant, never the plan coordinate. Each row corrects the
exposition its change falsifies, in the sections named; no new chapters.

**Anneal denylist** (grep enrolled durable files at ◆, excluding PLAN/ROADMAP/NOTES/SURVEY/AGENTS):

    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat-?[ABCI]\b
    \bC-[A-Z][A-Za-z0-9]+\b
    \b[A-H]\.[A-Z0-9]+\b|[Pp]rinciple[- ]?[0-9]+
    \bTrack [A-Z]\b|\bPhase [0-9]+\b|\bPhase [α-ω]|[α-ω]\.[0-9]+
    (PLAN|ROADMAP|NOTES|SURVEY)(\.md)?\b|plan-(run|shard|juncture)|\b(the|per the) plan\b

KAT, algorithmic step numbers, group notation (S₃, S4), bibliography keys ([S09]), E₁/E₂
filtration indices and quantum phase terminology are domain terms: classify collisions, do not
rename.

**First run:** after committing the planning docs, from `@plan-admin`:

    /plan-run docs/PLAN.md
