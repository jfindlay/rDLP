# ALIGN — workspace development policy

## Purpose (design intent)

Reconcile rDLP's development formalities with rust-template without changing mathematical behavior.
Implement SURVEY findings F-D1-01 through F-D1-09 here; prepare a separate rust-template handoff for
F-D1-10 through F-D1-12 and the approved root-toolchain policy. Do not mix repositories in one run.

This supersedes the completed REFACTOR ledger, preserved in git history. Its layout/coherence
contracts remain binding. Completing this graph closes the rDLP portion, not all bidirectional ALIGN.

Use this graph as the session-economics step-3 pilot described in
/home/justin/Source/opencode-config/PLAN.md. Pilot evidence establishes compatibility, not savings.

## Verify gate

VERIFY_TYPES: `cargo check --workspace`
VERIFY_TEST: `cargo test --workspace`

Both must pass at every session boundary. Preserve existing mathematical KATs and ignored-oracle
behavior. Do not add `--include-ignored` or silently substitute benchmark execution for tests.

User-approved ALIGN-only acceptance rule: configuration/documentation rows use the explicit native
assertions below instead of /plan-run's requirement for a new or changed test file. The driver must
bind this interpretation before dispatch. No compiler, regression, scope, or commit gate is waived.

New gates become mandatory after their introducing row succeeds:
- B.3: `cargo deny check`.
- B.4: `cargo fmt --all -- --check`.
- B.5: `cargo clippy --workspace --all-targets -- -D warnings`.
- B.6: `cargo llvm-cov --workspace --fail-under-lines <approved-floor>`.

Prospective gate failures are baseline evidence before adoption, not successful gate results.
A missing required tool blocks adoption. Tool installation requires separate approval.

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| B.1 | Record workspace quality-gate readiness | A | T1 | C-Scope | docs/NOTES.md; AGENTS.md |
| B.2 | Declare inherited workspace licensing | A | T1 | C-Scope; C-Baseline | Manifests M; LICENSE |
| B.3a | Remediate advisory-blocking transitive pins | A | T1 | C-Baseline | Cargo.lock |
| B.3 | Add workspace dependency audit policy | A | T1 | C-Baseline; C-License; B.3a | deny.toml; .cargo/config.toml |
| B.4 | Normalize workspace formatting and toolchain scope | A | T1 | C-Baseline | rust-toolchain.toml; rho/rust-toolchain.toml; rustfmt.toml; .cargo/config.toml |
| B.5 | Enforce inherited workspace lint policy | A | T1 | C-Baseline; C-Toolchain | Manifests M; .cargo/config.toml |
| B.6 | Add the measured workspace coverage gate | A | T1 | C-Baseline; C-Toolchain | .cargo/config.toml; docs/NOTES.md |
| B.7 | ◆ Publish the contributor guide and template-seeding handoff | I | T1 | All contracts above | docs/development.md; AGENTS.md; docs/NOTES.md |

Manifests M means exactly these paths; expand them in dispatch and commit filesets:
`Cargo.toml`, `gnfs/Cargo.toml`, `rho/Cargo.toml`, `shor/Cargo.toml`,
`shared/field/Cargo.toml`, `shared/bigint/Cargo.toml`, `shared/numth/Cargo.toml`,
`shared/numfield/Cargo.toml`, `shared/padic/Cargo.toml`, `shared/gf2m/Cargo.toml`.

All implementation workers are T1. The configured T0 @plan-juncture handles discovery adjudication
and the final boundary review. No speculative @architect inflection marker is needed.

## Session detail

### B.1 — Verification readiness
Deliver an appended baseline in NOTES and a minimal AGENTS guide naming working workspace commands,
oracle policy, and register discipline. Do not copy unrelated template application conventions.
Inventory installed toolchains/components and developer tools; select an exact stable release for
the root pin. Measure prospective format/lint/coverage/audit readiness where tools are available.
Probe the intended workspace-wide lint policy, not merely the three existing opted-in manifests.
Record missing tools and prospective failures honestly; group remediation by real file/contract seams.
Acceptance/KAT: both VERIFY commands pass; Cargo metadata resolves the nine declared workspace
members; recorded commands/results correspond to this tree and selected toolchain.
Missing optional tools may remain explicit blockers for later rows; B.1 does not certify those gates.

### B.2 — Licensing
Add the template's GPL license text and GPL-3.0-or-later workspace package metadata, inherited by every
member. Preserve package identities, dependencies, editions, and layout.
Acceptance/KAT: `cargo metadata --no-deps --format-version 1` resolves GPL-3.0-or-later for every
workspace package; LICENSE matches the intended GPL text; VERIFY passes.
The broad manifest fileset is one atomic inheritance change, not unrelated per-crate cleanup.

### B.3a — Advisory remediation (inserted by discovery-adjudication ahead of B.3)
A live cargo-deny/cargo-audit run against the actual dependency graph found RUSTSEC-2026-0204 (an
invalid-pointer-dereference vulnerability in `crossbeam-epoch`, reached via `rayon`, a direct
non-dev dependency of `rho`) and RUSTSEC-2026-0190 (a lower-severity "unsound" finding in `anyhow`,
reachable only through a wasm-target, proc-macro-time path). Both patches are available as ordinary
semver-compatible bumps. Bump `crossbeam-epoch` to >=0.9.20 and `anyhow` to >=1.0.103 via
`cargo update -p <crate>`, confined to `Cargo.lock`. No manifest edit, no `[advisories] ignore`
entry, no source/registry widening, no `rayon` demotion or removal, and no edit to the pre-existing
`rho/Cargo.lock` residue (that stays out of scope; a later toolchain-scope row owns it).
Acceptance/KAT: `cargo update -p crossbeam-epoch` and `-p anyhow` land only in `Cargo.lock`; VERIFY
passes; `cargo metadata` still resolves every license in the dependency-policy allow-list unchanged.
This row exists solely to unblock B.3's own stated acceptance criteria; it freezes no contract.

### B.3 — Dependency audit
Adapt the template's advisory/license/source policy to the installed cargo-deny schema and actual
dependency graph. Add an alias for the policy check; do not name it literally `audit` if doing so
would shadow the real `cargo-audit` binary (confirmed: it does, silently, under a deprecation
warning rather than a hard error) — name it for what it checks (e.g. `dep-policy`) and leave
`cargo audit` as the unaliased binary, so the KAT's two named commands stay two distinguishable
things. Do not copy comments pointing to template-local NOTES.
Acceptance/KAT: `cargo deny check` and `cargo audit` pass under the same policy.
Do not suppress advisories, widen trusted sources, or update dependencies silently to get green.
Necessary dependency remediation is a discovery requiring a separately scoped decision (this is
what B.3a was for; by the time B.3 runs, the graph is already clean).

### B.4 — Root toolchain and formatting
Replace rho's floating-stable file with one root rust-toolchain.toml pinning the verified stable
release. Declare rustfmt, clippy, and llvm-tools-preview; do not trigger unapproved installations.
Add rustfmt.toml with max_width = 100 and the `format` / `fmt-check` aliases.
Acceptance/KAT: root and member-directory toolchain discovery resolve the intended release absent
external overrides; no member toolchain files remain; `cargo fmt-check` and VERIFY pass.
Repeat relevant baselines if the selected release differs from the one initially measured.
Broad formatting changes are not authorized by this config fileset; surface and shard them separately.

### B.5 — Workspace lints
Define workspace Rust lints (unsafe_code forbid; missing_docs warn) and Clippy groups (all deny;
pedantic warn, explicit group priorities). Replace existing member duplicates with inheritance and
opt every member in. Add the workspace-wide strict `lint` alias.
Acceptance/KAT: all members inherit the policy; `cargo lint` passes with warnings denied; VERIFY passes.
If baseline findings require source/test remediation, stop for measured remediation shards before
adopting the gate. No broad lint suppressions, hidden behavioral changes, or catch-all cleanup.

### B.6 — Coverage
Measure workspace coverage under the pinned toolchain and record command, denominator, and result.
The ratified starting floor is 80%; any recalibration needs an explicit evidence-backed decision,
not automatic lowering to the observed result. Add the resulting `coverage` alias.
Acceptance/KAT: `cargo coverage` enforces the recorded floor and passes; VERIFY passes.
Do not execute ignored oracle tests or exclude arbitrary source files to raise the score.
Not running ignored tests does not automatically remove their otherwise-uncovered code from the
coverage denominator. Record the actual semantics.

### B.7 — Contributor guide and template handoff
Create docs/development.md from the now-verified commands and policies; make AGENTS point to it.
Explain setup, exact-toolchain upgrades, check commands, KAT-first coverage, and opt-in oracles.
Record the rDLP-side C-TemplateSeed handoff in NOTES: findings consumed, adapted conventions, and
concrete rust-template targets. Template execution remains separate.
Acceptance/KAT: documented commands agree with aliases and run successfully; relative links resolve;
the native assertions above and all adopted gates pass; planning-register anneal is clean.
Do not sweep unrelated mathematical prose, fix earlier implementation loose ends, or claim the
template repository has changed.

## Cross-session contracts

### C-Scope — prose-enforced; frozen by approved shard
Defined-in: SURVEY findings and this approval. Consumed-by: B.1–B.7.
Preserve C-Layout, mathematical behavior, oracle isolation, and topic-native durable prose.
Only enumerated configuration/manifest/contributor-guide work is authorized. Existing REFACTOR
residue remains CONSOLIDATE's responsibility. rust-template changes require their own execution plan.

### C-Baseline — prose/native-check-enforced; to be frozen at B.1
Defined-in: B.1. Consumed-by: B.2–B.7.
Verified workspace identity, exact-toolchain candidate, existing green regression commands, and
prospective-gate findings with missing data distinguished from success. This is readiness evidence,
not permission to treat failed new gates as adopted.

### C-License — native-check-enforced; to be frozen at B.2
Defined-in: B.2. Consumed-by: B.3, B.7.
All workspace members resolve GPL-3.0-or-later and the root license text is present.

### C-Audit — native-check-enforced; to be frozen at B.3
Defined-in: B.3. Consumed-by: subsequent gates and B.7.
The validated dependency policy and audit alias agree; exceptions require explicit justification.

### C-Toolchain — native-check/prose-enforced; to be frozen at B.4
Defined-in: B.4. Consumed-by: B.5–B.7.
One exact stable root toolchain, required development components, no member overrides, and
100-column rustfmt policy. This selects development tools; it does not establish an untested MSRV.

### C-Lints — compiler/native-check-enforced; to be frozen at B.5
Defined-in: B.5. Consumed-by: B.6–B.7.
Every member inherits the stated workspace policy and the strict workspace lint command is green.

### C-Coverage — native-check/prose-enforced; to be frozen at B.6
Defined-in: B.6. Consumed-by: B.7 and future test work.
Mathematical KATs remain primary; the measured, approved line floor is secondary. Command,
denominator, oracle exclusions, and threshold rationale are explicit.

### C-TemplateSeed — prose-enforced; rDLP side to be frozen at B.7
Defined-in: B.7. Consumed-by: a separate rust-template plan.
Seed root-toolchain policy, appropriately sized rolling-context conventions, opt-in oracle guidance,
and docs-reference discipline into the template. Do not prematurely freeze CONSOLIDATE's C-DocsLayer.
Expected template targets: rust-toolchain.toml, AGENTS.md, docs/development.md; confirm in that repo.
The full bidirectional contract remains pending until the template changes are verified.

## Progress ledger

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| B.1 | Record workspace quality-gate readiness | done | c1b0832 | C-Baseline (prose/native-check) |
| B.2 | Declare inherited workspace licensing | done | 0bfe20c | C-License (native-check) |
| B.3a | Remediate advisory-blocking transitive pins | pending | — | — |
| B.3 | Add workspace dependency audit policy | pending | — | — |
| B.4 | Normalize workspace formatting and toolchain scope | pending | — | — |
| B.5 | Enforce inherited workspace lint policy | pending | — | — |
| B.6 | Add the measured workspace coverage gate | pending | — | — |
| B.7 | ◆ Publish the contributor guide and template-seeding handoff | pending | — | — |

## Action-frame digest

### Preflight — 2026-09-18
Discovery/flex: `/plan-run` loop start found the working tree dirty — `.cargo/config.toml` and
`deny.toml` present as untracked files, matching B.3's expected fileset but never committed (`git
log --all` shows zero history for either path). `Cargo.lock` still resolves `crossbeam-epoch`
0.9.18, confirming B.3a's bump has not run — these are orphaned scratch artifacts from the earlier
B.3 dispatch that stopped mid-session at the RUSTSEC discovery recorded below, not a post-B.3a
completion left uncommitted.
Affected: none — no frozen contract is touched, and B.3's own ledger row is still `pending`, so
nothing here is a completed-but-unrecorded session.
Deferred: yes — a human must choose whether to `git clean` the stale pre-B.3a draft so B.3a
dispatches into a clean tree (B.3 gets freshly (re)dispatched afterward against the bumped lockfile,
per B.3's own "by the time B.3 runs, the graph is already clean" provision), or to inspect the draft
first for salvageable content. Precondition 5 (clean tree at loop start) fails until resolved; the
driver halted rather than choosing, per its "do not stomp uncommitted work" mandate.
Texture: this is the first `/plan-run` invocation to reach the loop since B.3's discovery-adjudication
reshard landed B.3a ahead of B.3 — the aborted dispatch's leftovers were never cleaned up before that
commit, so cold-resume here needed one extra disambiguating check (`git log --all` on the two paths)
before it was safe to conclude "orphaned scratch," not "an uncommitted finished session."

### B.3 — 2026-09-18
Discovery/flex: cargo-deny 0.20.2 and cargo-audit 0.22.2 got installed this session (approved
tool-installation exception); a live `cargo deny check`/`cargo audit` run against the actual graph
found RUSTSEC-2026-0204 (`crossbeam-epoch` 0.9.18, reached via `rayon` 1.12.0, a direct non-dev
dependency of `rho`) as a hard-failing vulnerability, plus RUSTSEC-2026-0190 (`anyhow` 1.0.102,
"unsound", reachable only via a wasm-target proc-macro path) as a cargo-audit-only warning that
cargo-deny's advisories check doesn't surface under schema defaults (confirmed: separate, both-
current advisory-db clones — a design difference, not staleness). B.3's own Acceptance/KAT could not
be reached without a dependency bump its own prose puts out of scope, so the implementer correctly
stopped and reported rather than resolving it. Forked `@plan-juncture` (discovery-adjudication);
verdict `additive-reshard` — inserted B.3a (Cargo.lock-only bump of both crates) ahead of B.3, user-
authorized (run has no `may-reshard`, so the fork's recommendation was surfaced, not auto-applied).
Also endorsed the implementer's `dep-policy` alias naming (naming it literally `audit` would silently
shadow the real `cargo-audit` binary — a genuine rust-template defect, carried to B.7's handoff).
Affected: C-Audit (not invalidated — its "agree; exceptions require explicit justification" wording
anticipates a live advisory; freezing is deferred until B.3a lands and B.3 re-runs clean). C-Scope,
C-Baseline, C-License all confirmed intact.
Deferred: yes — B.4 should treat `rho/Cargo.lock`'s independent copy of the same vulnerable pins as
load-bearing residue removal (a stale, vulnerable second lockfile), not merely cosmetic toolchain-
scope cleanup, per the B.1 baseline's original deferral. B.7 should carry the audit-alias-shadowing
hazard forward as a concrete rust-template target in the C-TemplateSeed handoff.
Texture: the juncture fork independently corrected the implementer's framing of the crossbeam-epoch
fix as possibly upstream-blocked ("or waiting on rayon-core") — it verified from the lockfile that
`rayon-core`'s requirement is an ordinary caret range, making the fix a pure lockfile bump available
today, not a wait on an upstream release. Worth remembering: a subagent's own diagnosis of "blocked
on upstream" is worth an independent recheck before it's treated as a hard constraint.

### B.1 — 2026-09-17
Discovery/flex: root toolchain candidate is rustc/cargo 1.95.0; `cargo fmt --all -- --check` fails
workspace-wide (no rustfmt.toml, never formatted, ~1000+ hunks across ~170 files); `cargo clippy
--workspace --all-targets` fails outright — the three manifests already opted into `[lints.clippy]
all = "deny"` (`shor`, `shared/padic`, `shared/gf2m`) hit hard compile errors under the installed
clippy 0.1.95, independent of any workspace-wide lint-adoption decision; `cargo-deny`, `cargo-audit`,
`cargo-llvm-cov`/`llvm-tools-preview`, `cargo-tarpaulin` are all absent (not checked, not failed); no
LICENSE file anywhere, every package resolves `license: null`; `rho/rust-toolchain.toml` (floating
`stable`) and a second independent `rho/Cargo.lock` are pre-workspace residue, currently non-
conflicting.
Affected: C-Baseline (produced, frozen as above); none of this contradicts C-Scope (the only other
frozen contract) — all findings are within the PLAN's own pre-anticipated "internal-continue" class
(Discoveries & risks: "coverage is unmeasured... cargo-deny and cargo-llvm-cov were not found...
recheck at B.1").
Deferred: yes — B.5's session detail already anticipates exactly the already-broken-clippy-policy
case ("If baseline findings require source/test remediation, stop for measured remediation shards
before adopting the gate"); the next juncture/executor touching B.5 should re-examine whether fixing
`shor`/`shared/padic`/`shared/gf2m`'s existing lint violations is a prerequisite shard before B.5
can adopt a workspace-wide policy. B.4 should re-examine the `rho/rust-toolchain.toml` +
`rho/Cargo.lock` residue when normalizing toolchain scope.
Texture: format/lint readiness is much rougher than "probe and adopt" — B.4's formatting row is
realistically a one-time whole-workspace reformat, not incremental; B.5 inherits a lint policy that
is already red on 3/9 crates before any workspace-wide widening is even considered.

## Discoveries & risks

- Confirmed: prior REFACTOR rows are done; durable discoveries are already in ROADMAP/NOTES.
- Confirmed survey evidence: coverage is unmeasured; lint inheritance currently covers only three
  members. cargo-deny and cargo-llvm-cov were not found in the surveyed environment; recheck at B.1.
- Internal-continue: command/schema adaptation and documented missing optional tools at baseline.
- Additive-reshard: measured formatting/lint/dependency remediation exceeds the config filesets.
  With halt-at-boundaries, surface the proposed shards; do not insert and execute them silently.
- Destructive-HALT: a proposed change invalidates frozen layout, behavior, oracle, or licensing
  intent and cannot be reconciled. Do not weaken gates to close a row.
- Actioned additive-reshard (B.3, 2026-09-18): a live `cargo deny check`/`cargo audit` run found a
  real, current advisory (RUSTSEC-2026-0204, `crossbeam-epoch` via `rayon`, a direct non-dev
  dependency of `rho`) blocking B.3's own stated Acceptance/KAT, plus a lower-severity, unreachable-
  in-this-workspace advisory (RUSTSEC-2026-0190, `anyhow`). Neither is a frozen-contract conflict —
  C-Audit's "agree; exceptions require explicit justification" wording anticipates exactly this case
  — so `@plan-juncture` returned `additive-reshard`, not `destructive-HALT`. User authorized inserting
  B.3a (Cargo.lock-only bump of both `crossbeam-epoch` and `anyhow`) ahead of B.3; see B.3a's session
  detail. The `rho/Cargo.lock` residue (deferred to B.4 at the B.1 baseline) now also carries the
  same vulnerable pins — its removal at B.4 is load-bearing (a stale, vulnerable second lockfile),
  not merely cosmetic scope cleanup.
- Free finding for C-TemplateSeed (B.3, 2026-09-18): naming a `.cargo/config.toml` alias literally
  `audit` silently shadows the real `cargo-audit` binary (a cargo deprecation warning today, a hard
  error in some future cargo release) — this is a genuine rust-template defect, already inside B.7's
  existing handoff scope. B.3 named its alias `dep-policy` instead; B.7 should carry the shadowing
  hazard itself forward as a concrete rust-template target, not just the alias name.
- Runtime/plugin loading, live attribution, accounting, and output persistence remain unverified.
  Deterministic tests recorded as passing in the config plan are not live integration evidence.
- There is no estimate that this entire graph fits within $10. Unfinished rows remain pending.

## Notes for executors

Five-input sizing: moderate ambient complexity; inheritance supplies the broad-row irreducible
floor; gate design-error cost and correctness-criticality favor separate adoption rows; strong
behavioral KATs but unknown lint/coverage readiness require baseline discovery first.
Outlier check: no broad source/test row, no repeated per-file cleanup, no catch-all final fixer.

REGISTER: state the property/reason/invariant, never the plan coordinate, in source/tests/human docs.
B.1 seeds this rule in AGENTS; sharding does not write AGENTS in architect mode.
Anneal only ALIGN-owned durable files, excluding rolling context and AGENTS. Denylist:
`\bB\.[1-7][a-z]?\b`, `sub-track`, `frozen contract`, `PLAN\.md`, `ROADMAP\.md`, `action-frame`,
`◆`, `\bCat [ABC]\b`, `plan-(run|shard|juncture)`. Do not sweep pre-existing unrelated prose.
(Widened to `[a-z]?` at the B.3 additive-reshard to also catch the inserted B.3a token.)

Pilot authorization: user approved real rDLP execution instead of disposable copies and $10 TOTAL
for the pilot, not per session. Start accounting at the fresh pilot-preflight session; include all
enrolled parents, workers, fixes, committers, reviews, and restarts. This sharding session is excluded.
There is no separate template budget. No pushes or automatic reversions.

Before live execution, @build runs the config plan's deterministic gate and checks runtime/plugin
loading. Establish a workable external spend-cap arrangement; the plugin cannot enforce the cap.
If only best-effort monitoring is available, disclose that and obtain approval before live calls.
Missing or zero-priced model accounting cannot certify remaining budget.

Use fresh observe-only startup first. First /plan-run dispatch authorization is B.1 ONLY: run its
normal verification, session commit, and ledger commit, then stop before B.2. This explicit pilot
stop is not a CLI flag and does not falsely mark B.1 as a sub-track-final row.
Review observe evidence; restart a fresh advisory-mode session before authorizing subsequent rows.
Carry cumulative spend across the restart. Do not restart or silently enroll existing sessions.

Suggested invocation from @plan-admin: `/plan-run docs/PLAN.md halt-at-boundaries`, accompanied by
the B.1-only observe-run instruction above. Resume remaining rows only after the pilot checkpoint.
Retain existing commit, regression, discovery, and scope gates. Budget exhaustion never authorizes
a false green commit; preserve/report incomplete work and request further authorization.

Record session IDs, modes, accounting comparisons, notice overhead, output preservation, failures,
and observed strategy changes in the opencode-config pilot report, not shipped rDLP documentation.
Task-boundary reports need no threshold crossing. Worker notices do; if none fires naturally,
mark live notice delivery unverified rather than inflate workload or change thresholds covertly.
