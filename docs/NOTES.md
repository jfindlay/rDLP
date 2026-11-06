# rDLP — Notes

Rolling-context durable framings, decisions, and mental models. Distinct from `docs/PLAN.md`
(current sub-track) and `docs/ROADMAP.md` (project-lifetime view). Captured framings that prove
durable graduate to the ROADMAP Discoveries log at a sub-track boundary.

---

## Transfer-style attacks shard into transfer / structure / solve (2026-06-15, E.H shard)

The small-characteristic index-calculus attack on binary-curve ECDLP decomposes into three
sub-tracks at clean contract seams: **E.H** (GHS/Weil descent) *transfers* an ECDLP on `E/GF(2^m)`
to a DLP on `Jac(C)/GF(2^l)`; **E.J** builds the Semaev summation polynomials; **E.K** runs the
index-calculus *solve*.

## D. REFACTOR freeze record (2026-06-21, D.6 ◆)

Captured at the D.6 ◆ freeze. These are the durable structural facts that bind CONSOLIDATE and
EXTEND.

### C-Layout freeze (compiler-enforced, FROZEN at D.6 ◆)

The final crate layout after D.1's wrapper collapse:

- **Removed:** `rho/src/field/mod.rs`, `rho/src/field/monty.rs` (dead code — never declared in
  `mod.rs`), `rho/src/field/naive.rs` (dead code — never declared in `mod.rs`),
  `rho/src/util/batch_inv.rs` (re-export shim + duplicated `batch_invert` tests),
  `rho/src/util/mp.rs` (empty stub), `rho/src/util/mod.rs` (collapsed away).
- **`rho` now imports directly:** `shared_field::{Fp, FpNaive4, FpMonty4}` and
  `shared_bigint::batch_invert` at call-sites (no intermediate re-export layer).
- **Preserved:** `F: Fp<4>` trait bounds throughout `rho` (the trait is genuinely generic; only the
  local aliases were contrivance). No crate split (open-Q 1 ratified: Track-E attack modules stay
  co-located with the Pollard-rho baseline). No other `shared/*` dedup (F-D5-03).
- **Proof:** `cargo check --workspace` and `cargo test --workspace` green at D.5 terminus
  (`3ef9aca`). The compiler-enforced contract is the proof.
- **Binds:** CONSOLIDATE (docs reference the final layout) and EXTEND (new code sits in it).

### `shor/Cargo.toml` dev-dep addition (structural necessity of D.1)

Before the collapse, `shor/tests/ecdlp_kat.rs` imported `rho::field::FpMonty` (the re-export
wrapper). Once `rho::field` was removed, the compiler forced a direct `shared-field` dev-dep in
`shor/Cargo.toml`. This is a mechanical consequence of the collapse, not a REFACTOR-scope
violation. Recorded here because it is a structural fact about the workspace's dependency graph.

### C-Coherence code-half freeze (FROZEN at D.6 ◆, with ~1% residue for CONSOLIDATE)

The code-depth de-provenancing catalog (F-D9-01…06, ~524 tokens) is consumed:

- **Scrubbed across D.2–D.5** (`89f1e22`, `84b993c`, `f2302a4`, `3ef9aca`): `Phase N` planning
  labels → optimization-layer names; `Track D` → topic-native language; sub-track IDs + contract-
  label tokens (`C-AnomalousLift`, `C-IndexCalc`, `C-Pairing`) → topic-native references or removal.
- **Identifier renames (D.3):** `sub_track_close_curve_axioms_intact` → `binary_curve_axioms_intact`
  (in `rho/tests/binary_curve_kat.rs`); `c_ek_relation_round_trip` → `relation_round_trip` (in
  `rho/tests/index_calculus_kat.rs`, contract-label prefix removal).
- **Three load-bearing carve-outs held:** algorithmic `Phase 1/2/3` step-indices in
  `rho/src/ecdlp/walk.rs` and `rho/src/index_calculus/solve.rs` preserved (not planning residue);
  quantum-phase-estimation `Phase N/M` tokens in `shor/tests/factor_kat.rs` preserved (math term).
- **Interface-contract vocabulary preserved:** `C-StateVec`, `C-PointAdd`, `C-ModExp`, `C-QFT`,
  `C-Factor`, `C-OrderFind`, `C-Sparse`, `C-Relation`, `C-Ideal`, `C-Res`, `C-NF`, `C-AlgSqrt`,
  `C-Padic`, `C-Hensel`, `C-PadicLog`, `C-HyperCurve`, `C-MovBridge` — these are arc-1 design
  vocabulary naming living interfaces, not planning-frame residue.
- **~1% residue for CONSOLIDATE:** five tokens escaped the sweep — `rho/src/curve/test_curves.rs:36`
  (`Pohlig–Hellman (E.A.2)`), `:139` (`Pohlig–Hellman (E.A)`), `rho/benches/attacks.rs:20,26`
  (`E.W table` / `E.W.2 chapter`), `rho/tests/hyperelliptic_kat.rs:6` (`## E.I.2`). All
  behaviourally inert; recommended for CONSOLIDATE (which will be in these files' neighbourhood).
- **KAT suite identically green throughout** (inertness invariant held across all four sessions).

---

## Transfer-style attacks shard into transfer / structure / solve (2026-06-15, E.H shard)

The load-bearing observation: E.H's correctness signal is the **homomorphism + logarithm-preservation
relation** (`D_{P+Q} = D_P + D_Q`, and `D_h = k·D_g` for known `k` via the frozen Cantor
`scalar_mul`) — the transfer is verified by *relationship-preservation*, with the solve delegated
downstream. This is the same "build the structure, delegate the attack" idiom seen across the
project:

- **E.I** built the Jacobian group law (Cantor on Mumford pairs) and left the DLP solver to its
  consumer.
- **MOV / SSA** (E.C / E.E) transferred a prime-field ECDLP to a different DLP and verified the
  transfer by the preserved log, calling a solver (NFS-DL / p-adic log) rather than re-deriving one.

Consequence for sharding: a transfer attack freezes the *relocated problem* (here C-GHSDescent)
without ever solving it; its sub-track-close KAT is relationship-preservation, not an end-to-end
break. Bundling the solve into the transfer sub-track couples two independently-delicate designs —
the coupling the inflection-juncture discipline warns against. This framing should inform E.J/E.K
sharding and the E.W cross-attack synthesis.

---

## SURVEY durable framings (2026-06-21, A.7 ◆)

Captured at the A.7 ◆ findings-ledger freeze. These framings proved durable across all six audit
sessions and are load-bearing for the five downstream campaigns.

### Scope-routing invariant

The SURVEY split execution by artifact register, and the routing is binding for all five campaigns:

> **Inline `//!`/`///`/`//` doc-comments + identifier/file/dir/bench/test renames → REFACTOR
> (compiler blast radius); human prose (PEDAGOGY/MATHEMATICS/README) → CONSOLIDATE; config/
> manifest/CI → ALIGN; CADO sidecar build → ORACLE; new algorithm/chapter → EXTEND.**

A finding routed to the wrong campaign is the primary cross-campaign defocus mode. The boundary
is sharp: inline doc-comments are code-adjacent (REFACTOR); human prose docs are CONSOLIDATE's.
The `PLAN.md E.D.2` reference in `shared/padic/tests/hensel_kat.rs` is a test file doc-comment
(code-adjacent) → REFACTOR, not CONSOLIDATE.

### The findings-vs-fixes line (the central SURVEY defocus risk)

A session that edits a `.rs`, manifest, config, or doc *content* (vs the SURVEY findings file)
has broken scope. **SURVEY records what should change; the campaigns change it.** Any code/config/
refactor edit is a HALT-and-surface defocus signal — the audit-substrate analogue of "Z.1/T.Z
must not write code." This line held across all six audit sessions.

### C-Testing-Philosophy doctrine (frozen at A.3, ratified at A.7)

The resolved doctrine for "what is the right testing target for toy-scale compute-heavy KAT-driven
pedagogical code":

- **Primary target:** mathematical-behavior KAT coverage. Field axioms, group laws, and
  algorithm-level correctness are the primary signals. Line coverage is secondary — a floor, not
  a ceiling.
- **Coverage gate:** 80% line coverage as a dead-code floor (not a correctness ceiling). Gate value
  is a recommendation; ALIGN calibrates against a measured baseline. The `#[ignore]`-gated oracle
  tests are excluded from the gate.
- **Two-tier pattern:** inline `#[test]` for substrate verification; external `tests/*_kat.rs` for
  algorithm-level KATs. New tests (CONSOLIDATE, F-EXTEND) follow this split.
- **Test-length norm:** one file per algorithm family, however long. Splitting is not warranted by
  length alone.
- **Oracle-gate norm:** `#[ignore]` for all external-tool and compute-heavy tests. The arc-1
  dev-oracle policy; new tests honor it.

The key insight: the entire test suite is KAT-based (no proptest/quickcheck); the compute-heavy
paths are already exercised at toy scale; the two-tier structure (inline substrate + external KAT)
is healthy and intentional. A line-coverage gate that misses a mathematical property is a weaker
correctness signal than a KAT-completeness check.

### C-Coherence classification defaults (frozen at A.7)

The artifact-stands-on-its-own-terms principle, applied to the SURVEY findings:

- **Code-half:** doc-comment-heavy, identifier-light. Exactly one planning-frame token in a code
  identifier (`sub_track_close_curve_axioms_intact`). All other tokens are in `//!`/`///`/`//`
  comments. Identifier renames are cheap (compiler-checked); doc-comment scrubbing is the bulk of
  REFACTOR's de-provenancing work.
- **Prose-half:** ~185 tokens across 6 files. ~110 are "Track X" tokens (grouping-coincides-with-
  topic — the five tracks map onto five real mathematical families). ~75 are pure residue (Phase N,
  session IDs, ◆ marks, fine-grained sub-track IDs).
- **Borderline default (open-Q 4):** for coarse groupings where track ≈ topic (S.A/S.B/S.C Shor
  sub-tracks; G.B–G.W GNFS pipeline stages), the adjudicator recommendation is
  **preserve-under-topic-label** — the coarse groupings map onto real mathematical families; only
  the fine-grained S.X.Y / G.X.Y session IDs are pure residue and dissolve.
- **The key distinction:** "Track ρ" → "Pollard rho" is a label change (grouping survives);
  "Phase 5" → "distinguished-point search" is a re-anchor (planning label replaced by topic name);
  "E.K.3 ◆" → removed (pure residue with no topic-native replacement).
## D. REFACTOR freeze record (2026-06-21, D.6 ◆)

Captured at the D.6 ◆ C-Layout + C-Coherence code-half freeze. These are the durable structural
facts that bind CONSOLIDATE and EXTEND.

### C-Layout freeze (compiler-enforced, FROZEN at D.6 ◆)

Final crate layout after the D.1 wrapper collapse:

- **Removed:** `rho/src/field/mod.rs`, `rho/src/field/monty.rs` (dead code — never declared in
  `mod.rs`), `rho/src/field/naive.rs` (dead code — same), `rho/src/util/batch_inv.rs` (re-export
  shim + duplicated `batch_invert` tests), `rho/src/util/mp.rs` (empty stub),
  `rho/src/util/mod.rs`.
- **`rho` now imports directly:** `shared_field::{Fp, FpNaive4, FpMonty4}` and
  `shared_bigint::batch_invert` at call-sites. The `F: Fp<4>` trait bounds are preserved (the trait
  is genuinely generic; only the local aliases and duplicated tests were contrivance).
- **No crate split** (open-Q 1 ratified): the 8 Track-E attack modules stay co-located with the
  Pollard-rho baseline in `rho`. The E.W cross-attack bench measures attacks against the rho baseline
  and they share `rho::curve` types directly — "attacks live with the baseline they're measured
  against" is load-bearing pedagogy.
- **No other `shared/*` dedup** (F-D5-03): the six `shared/*` crates have distinct, non-overlapping
  concerns.
- **Proof:** `cargo check --workspace` and `cargo test --workspace` green at D.5 terminus (`3ef9aca`).

**Structural necessity of the D.1 collapse:** `shor/tests/ecdlp_kat.rs` imported `rho::field::FpMonty`
(the re-export wrapper). Once `rho::field` was removed, the compiler forced a direct `shared-field`
dev-dep in `shor/Cargo.toml`. This is a mechanical consequence of the collapse, not a REFACTOR-scope
violation.

### C-Coherence code-half freeze (≈99% completeness, ~1% residue for CONSOLIDATE)

All ~524 cataloged F-D9-01…06 code-depth tokens re-anchored across D.2–D.5 (`89f1e22`, `84b993c`,
`f2302a4`, `3ef9aca`). The KAT suite stayed identically green throughout (inertness invariant held).

**Three load-bearing carve-outs held (verified):**
- Algorithmic `Phase 1/2/3` step-indices in `rho/src/ecdlp/walk.rs` (batched-inversion loop steps)
  and `Phase 1/2` in `rho/src/index_calculus/solve.rs` (Gaussian-elim steps) — these are numbered
  steps of one algorithm, not planning-frame session labels. Preserved.
- Quantum-phase-estimation `Phase N/M` tokens in `shor/tests/factor_kat.rs` — a mathematical term
  ("Phase 64/256 = 1/4"), not planning residue. Preserved.

**Codebase interface-contract vocabulary preserved:** `C-StateVec`, `C-PointAdd`, `C-ModExp`,
`C-QFT`, `C-Factor`, `C-OrderFind`, `C-Sparse`, `C-Relation`, `C-Ideal`, `C-Res`, `C-NF`,
`C-AlgSqrt`, `C-Padic`, `C-Hensel`, `C-PadicLog`, `C-HyperCurve`, `C-MovBridge` — these are arc-1
design vocabulary naming living interfaces, not planning-frame residue. Scrubbing them would be the
defocus failure mode.

**~1% residue for CONSOLIDATE (not a blocker):** Five pure-residue sub-track-ID tokens escaped the
D.2–D.5 sweep, all behaviourally inert doc-comment text in lower-density peripheral files:
- `rho/src/curve/test_curves.rs:36` — `Pohlig–Hellman (E.A.2)`
- `rho/src/curve/test_curves.rs:139` — `Pohlig–Hellman (E.A)`
- `rho/benches/attacks.rs:20,26` — `E.W table` / `E.W.2 chapter` (the line-1 `E.W.1` was scrubbed
  at D.3; interior refs were missed)
- `rho/tests/hyperelliptic_kat.rs:6` — `## E.I.2` section header

Recommended for CONSOLIDATE (which owns the prose-half coherence sweep and will be in these files'
neighbourhood). The code-half freeze stands with this residue logged.

### Identifier renames in D. REFACTOR

Two identifier renames were made (both compiler- and test-verified):
- `sub_track_close_curve_axioms_intact` → `binary_curve_axioms_intact` in
  `rho/tests/binary_curve_kat.rs:582` (D.3 scope — sub-track-ID removal from test name).
- `c_ek_relation_round_trip` → `relation_round_trip` in `rho/tests/index_calculus_kat.rs` (D.3
  scope — contract-label prefix `c_ek_` removed).

---

## Workspace verification-tooling baseline (2026-09-17)

A measured, honest snapshot of this workspace's build/test health and its readiness for format,
lint, dependency-audit, and coverage tooling, taken against one specific toolchain so that later
tool-adoption decisions have real evidence instead of assumptions. "Not checked" (tool absent) is
recorded separately from "checked and failed" (tool ran, did not pass) throughout.

### Toolchain and workspace identity — checked and passed

The active default stable toolchain, installed via `rustup` (`rustc`/`cargo` 1.95.0, with matching
`rustfmt` 1.9.0-stable and `clippy` 0.1.95 shipped alongside it), is verified working end to end and
is the recommended exact release for the workspace's root toolchain pin. `cargo metadata --no-deps
--format-version 1` resolves exactly the nine intended workspace members (`gnfs`, `rho`, `shor`,
`shared-bigint`, `shared-field`, `shared-gf2m`, `shared-numfield`, `shared-numth`, `shared-padic`).
`cargo check --workspace` and `cargo test --workspace` are both green: every test that runs passes;
the only non-passing outcome anywhere in the suite is the pre-existing `#[ignore]` gate on tests
that need an external PARI or CADO-NFS oracle, which is deliberate and orthogonal to this baseline.

### Format readiness — checked and failed

`cargo fmt --all -- --check` fails. No `rustfmt.toml` exists anywhere in the tree, so this ran under
rustfmt's built-in defaults, and the codebase appears never to have been passed through rustfmt at
all: the diff touches essentially every source and test file in every one of the nine crates (on
the order of a thousand hunks across ~170 files), not a handful of stragglers. Adopting a formatting
gate here is a one-time whole-workspace reformat, not an incremental cleanup.

### Lint readiness — checked and failed for the existing policy; mixed at default level elsewhere

`cargo clippy --workspace --all-targets` fails outright — nonzero exit even without `-D warnings`.
Exactly three manifests opt into a lint policy today (`shor`, `shared/padic`, `shared/gf2m`, each
with `[lints.clippy] all = "deny"`), and all three now produce hard compile errors, not warnings,
under the installed clippy: `shor` fails on several `needless_range_loop` and
`manual_is_multiple_of` sites, `shared/gf2m` on `needless_range_loop`, and `shared/padic` on a
`redundant_guard`-shaped lint.
The policy already committed today does not build clean against this environment's clippy version —
that is independent of any separate decision to widen the policy workspace-wide. The six crates with
no opt-in policy (`gnfs`, `rho`, `shared-bigint`, `shared-field`, `shared-numfield`, `shared-numth`)
only ever *warn* at clippy's default level; two of them (`rho`, `shared-field`) are already fully
clean, the rest carry a handful to a few dozen warnings apiece, concentrated in the largest crate
(`gnfs`). None of this is a required gate today — it is prospective-adoption evidence only.

### Audit and coverage readiness — not checked (tool absent)

`cargo-deny`, `cargo-audit`, `cargo-llvm-cov`, and `cargo-tarpaulin` are all absent from this
environment — `cargo <name>` reports "no such command" for each — and the `llvm-tools-preview`
`rustup` component that `cargo-llvm-cov` additionally needs is not installed either. None of these
tools were run; their eventual pass/fail is genuinely unmeasured, not a recorded failure, and
installing any of them is a separate decision from this baseline.

### Licensing — observed fact, not this baseline's gate

No `LICENSE` file exists at the repository root or in any member directory, and no member manifest
declares a `license` or `license-file` key: `cargo metadata` reports a null license for every one of
the nine packages.

### Toolchain-scope residue found while surveying

`rho/rust-toolchain.toml` (pinning the floating `stable` channel rather than an exact release) and
`rho/Cargo.lock` (a second, independent lockfile alongside the workspace root's) both date to the
repository's very first commit, before the Cargo workspace existed, and were never removed once
`rho` joined the workspace. Neither currently conflicts with the root `Cargo.toml`/`Cargo.lock` —
every command above ran clean — but both read as pre-workspace leftovers rather than an intentional
per-member override, worth resolving whenever the workspace's toolchain scope is normalized.

---

## Workspace coverage measurement (2026-09-22)

`cargo-llvm-cov` 0.9.1 is now installed against the pinned 1.95.0 toolchain (its
`llvm-tools-preview` requirement was already satisfied — see the toolchain baseline above), closing
the "not checked" gap the 2026-09-17 baseline recorded for coverage tooling specifically. This is
the first actual run.

**Command:** `cargo llvm-cov --workspace --fail-under-lines 80`, now also reachable as the
`cargo coverage` alias in `.cargo/config.toml`. Both invocations were run and both exit 0.

**Result:** measured line coverage **92.78%** (18594 total lines, 1342–1343 missed across two
successive runs — a single-line rounding-level fluctuation between builds, not a regression).
Region coverage measured 94.20–94.21%, function coverage 94.55%. The ratified 80% floor
(C-Testing-Philosophy, frozen at A.7: line coverage is a dead-code floor, not a correctness
ceiling) holds with comfortable margin, so no recalibration decision was needed or made this
session — the floor is recorded as-is.

**Denominator semantics, stated precisely:** `cargo llvm-cov --workspace` runs the same
`cargo test --workspace` invocation the workspace already gates on, so the `#[ignore]`-gated
PARI/CADO-NFS oracle tests are excluded from the coverage *run* exactly as they already are from
the plain test gate — no `--include-ignored` flag was added and no oracle binary was invoked. That
is a fact about which tests execute, not about which source lines count: their otherwise-uncovered
code is not excluded from the coverage *denominator*, and is counted as uncovered wherever an
oracle-only test would have been the sole thing exercising a line. No source file, module, or crate
was excluded from the denominator to inflate the score; `--ignore-filename-regex` was not used
anywhere.

**Shape of the uncovered residue (for whoever next revisits this floor):** two files show a
genuine 0.00% line coverage, both CLI-entry-point binaries never invoked as processes by the test
suite — `rho/src/ecdlp/cli.rs` (84 lines) and `rho/src/factor/cli.rs` (22 lines), thin
argument-parsing wrappers over already-tested library entry points. The weakest non-zero library
modules by line coverage: `shared/gf2m/src/normal.rs` (55.36%), `rho/src/index_calculus/mod.rs`
(69.44%), `gnfs/src/filter/merge.rs` (69.91%), `gnfs/src/dl/descent/node.rs` (72.45%),
`rho/src/ssa/mod.rs` (79.73%), `gnfs/src/polyselect/mod.rs` (80.18%), `shor/src/shor/mod.rs`
(80.54%), `shor/src/gates/mod.rs` (80.63%). None of this shape was used to trim the run or widen
the gate — recorded here as the honest per-crate picture behind the aggregate number.

**Not done this session:** no new test was written to close any gap, no ignored oracle test was
executed, and the 80% floor was not lowered — the measured number cleared it, so none of the three
forbidden moves were needed.

---

## C-TemplateSeed handoff — rDLP side (2026-09-22)

Records the rDLP-side half of C-TemplateSeed: what this campaign's SURVEY findings and its own
execution produced that is worth seeding into the `rust-template` repository. Template execution is
a separate decision made in that repository, by whoever owns it, on its own schedule; **no file in
`rust-template` has been read or changed by this session or any session in this campaign.**

**Findings consumed (now complete at rDLP):** F-D1-01 through F-D1-09 — `deny.toml`, `rustfmt.toml`,
`rust-toolchain.toml`, root `LICENSE`, workspace `[lints]` tables (with per-member
`[lints] workspace = true` inheritance replacing three previously-duplicated per-crate blocks),
the `cargo coverage` gate, `docs/development.md`, `.cargo/config.toml` aliases, and this campaign's
`AGENTS.md` edits. All nine are now real, verified state in this repository (see the coverage and
toolchain baselines above, and `docs/development.md` for the consolidated command reference).

**Adapted conventions mapped to C-TemplateSeed's three clauses (still pending in the template repo,
not evaluated there):**

- *Root-toolchain policy* — the exact-release pin-with-components pattern itself
  (`rust-toolchain.toml` pinning a specific release, e.g. `1.95.0`, rather than a floating channel,
  with `rustfmt`/`clippy`/`llvm-tools-preview` declared as required components) is a direct template
  target on its own, independent of the three F-D1-1x findings below.
- *Appropriately sized rolling-context conventions* (F-D1-10, multisession docs layering) — this
  repository's three-tier split (`docs/PLAN.md` for the live sub-track, `docs/ROADMAP.md` for the
  project-lifetime view, `docs/NOTES.md` for rolling dated framings that either graduate to ROADMAP
  or stay local) is the concrete shape to offer the template, sized down as needed for a template's
  typically smaller scope.
- *Opt-in oracle guidance* (F-D1-11, dev-oracle policy) — the `#[ignore]`-by-default convention for
  tests needing an external tool, named so the gate is discoverable (`kat_pari_*`, `kat_c_cado_*`
  here), run only via `cargo test -- --ignored`, and never promoted onto the default test path.
- *Docs-reference discipline* (F-D1-12, docs-register / three-layer reference discipline) — the
  register split this repository enforces between code-adjacent doc-comments, durable human prose
  (`README.md`, a contributor guide), and rolling context notes, with each layer stating which of
  the others it may reference and which it may not (durable prose never points at rolling notes;
  code never points at either).

**Concrete defects/patterns discovered while executing this plan, not from the original SURVEY
digest, also worth carrying into the template:**

- **Alias-naming shadowing hazard:** naming a `cargo-deny`-check alias literally `audit` in
  `.cargo/config.toml` silently shadows the real `cargo-audit` binary — cargo resolves an alias
  before a same-named installed binary, so `cargo audit` would silently run the alias instead of
  `cargo-audit`, today only a deprecation warning but a hard error in a future cargo release. This
  repository's alias is named `dep-policy` specifically to avoid the collision; a template
  bootstrapping both `cargo-deny` and `cargo-audit` should avoid the `audit` alias name from the
  start rather than rename it later.
- **The "declared-but-unadopted gate" pattern:** used twice in this plan (`cargo fmt-check` and
  `cargo lint` in this repository — see `docs/development.md`) for a gate whose alias/infrastructure
  is correct and runnable but whose workspace-wide conformance requires a separately authorized
  remediation pass on a pre-existing, previously-unformatted or previously-unlinted codebase. Worth
  naming as a reusable, honest intermediate state for any template bootstrapping onto an existing
  codebase rather than a fresh one: wire the alias and document it as declared-but-not-yet-enforced,
  rather than either skipping it entirely or claiming false conformance.

**Verification run this session (all commands actually executed, not assumed):**
`cargo check --workspace` and `cargo test --workspace` green (0 failures across all nine crates'
lib/integration/doc tests); `cargo dep-policy` and `cargo audit` both pass cleanly; `cargo coverage`
passes (92.7–92.8% measured, consistent with the 2026-09-22 coverage baseline above, single-line
fluctuation between runs as previously noted); `cargo fmt-check` and `cargo lint` both still fail,
unchanged from their recorded baseline — neither was touched, and neither is claimed as a gate in
`docs/development.md` or `AGENTS.md`.

---

## Correctness-first reassessment (2026-11-06)

Source inspection and a separate verification pass changed the next-work recommendation from
CONSOLIDATE to a bounded correctness campaign, beginning with GNFS. The user approved that priority
and a design-gated PLAN. No tests, builds, benchmarks, or oracles were run during this assessment.
The working tree was clean at `3cd6fe5` before the approved rolling-document edits.

### Confirmed facts and mathematical implications

- **GNFS mapping mismatch.** `gnfs/src/sqrt/algebraic.rs:154-169` returns `abs(Norm(beta)) mod N`;
  `docs/MATHEMATICS.md:1740-1744` specifies `beta(m) mod N`. The norm is multiplicative but is not
  the evaluation homomorphism required here. **Inferred:** the implemented result does not generally
  give the rational/algebraic square congruence. A source-level counterexample already appears in
  `gnfs/tests/sqrt_algebraic_kat.rs:91-134`: f=x²−2, N=7, m=3, gamma=9, beta=3. The test expects
  Y=2 and checks Y²=Norm(gamma) mod 7; the required rational product is 9 mod 7=2, whereas Y²=4.
  The test checks the wrong identity for the pipeline. No execution was needed for this arithmetic.
- **GNFS integration gap.** Inspected algebraic/factor fixtures use b=0 relations; the ordinary
  congruence tests use degree-one fields or special values. The full sieve-to-factor oracle test
  in `gnfs/tests/factor_end_to_end_kat.rs:271-286` is `todo!`. Existing stage tests are useful but
  are not evidence that the complete pipeline composes.
- **Reconstruction needs design review.** `algebraic.rs:99-129,252-302` uses heuristic root signs,
  constant-coefficient normalization across primes, fixed prime count, and local modular checks;
  there is no final exact beta²=gamma check. `gnfs/src/linalg/qc.rs` computes norm characters
  rather than individual root characters. **Open:** the justified reconstruction method, supported
  domain, obstruction semantics, and non-monic/basis corrections. Do not equate the confirmed final
  map defect with a complete repair specification.
- **NFS-DL placeholder.** `gnfs/src/dl/descent/solve.rs:483-501` discards the smoothing exponent and
  returns `Ok(0)` when its context-free k=1 frontier is empty. The alternative `solve_dl_full`
  exists and must not be conflated with that incomplete wrapper. Its value-bearing tests use
  supplied virtual logs and smooth targets. Structural descent tests do use sieve-produced
  relations, but mock logs do not verify the resulting logarithm. **Open:** sign conventions,
  normalization, Schirokauer terms, and target-kind handling require mathematical adjudication.
- **Fallback-masked success.** `shor/src/ecdlp/mod.rs:393-408` and
  `rho/src/index_calculus/solve.rs:135-147` enumerate candidate logarithms after their primary
  paths fail. Final-answer tests do not require the primary path to succeed. These fallbacks are
  documented locally; their presence does not prove the primary paths never work. It does mean
  final-answer success alone cannot demonstrate the advertised mechanism.
- **Oracle placeholders.** Besides the GNFS full-pipeline `todo!`,
  `gnfs/tests/merge_kat.rs:347` is `unimplemented!` and `rho/tests/mov_kat.rs:263` is an unconditional
  assertion. Ignoring external tests is correct policy; placeholder bodies are not completed oracle
  tests. Installing oracle binaries alone will not make them meaningful.

### Decision and boundaries

GNFS first: mathematical contract review, verified square reconstruction/map, dependency integration,
then a genuine toy sieve-to-factor test. A suitable fixture and repair algorithm remain undecided.
Implementation boundaries are provisional until the design review and user approval. Exact
verification may reject nonsquare kernel candidates; finite character checks are not assumed to
prove that every candidate is a field square. Domain restriction requires explicit approval.

NFS-DL follows as a separately scoped candidate; Shor/index-calculus algorithm-path evidence remains
queued. External-oracle comparisons should exercise repaired boundaries. Broad consolidation follows
the evidence, while narrow claim corrections accompany each repair. No re-opening of layout work,
automatic extension, or hidden format/lint remediation is authorized.

Historical campaign closure and the 92.78% line-coverage baseline remain valid records of their own
scope. They do not certify algorithm fidelity. Preserve the KAT-first doctrine, but strengthen which
mathematical contracts those KATs actually test.

### Capture candidate — algorithm-specific evidence

In a pedagogical algorithm library, a correct final answer does not demonstrate the named algorithm
when a fallback can supply it independently. Test the mechanism's contribution: require primary-path
success and verify intermediate identities, while testing fallback behavior separately. High line
coverage can coexist with this gap because executing a mechanism is not evidence that it produced
the answer. Candidate recorded for review; no global guidance was changed.

---

## GNFS mathematical contract and canonical reconstruction design (2026-11-06, G.1)

### Decision, evidence, and authorization

The user selected **canonical direct p-adic reconstruction plus the published split-prime CRT
fallback**, not a new reconstruction technique and not permanent inert-prime-only coverage. The
library's purpose is to reimplement, assemble, and ultimately unify established algorithms. Missing
substrate is a reason to implement the substrate, not to silently substitute a different algorithm.

The user approved recording this design, the proposed existing-workspace `shared-padic` dependency,
typed API changes, affine-collector boundary, and bounded integration discovery for review. This
authorized documentation writes only. Implementation, checks, installs, and commits remain blocked
until the revised PLAN is reviewed. No source, test, build, benchmark, or external oracle was run.
PLAN, NOTES, and ROADMAP already had uncommitted changes at session start; those were preserved.

**Confirmed** below means source or reference inspection. **Derived** labels mathematical reasoning
from the stated preconditions, not runtime verification. Fixture feasibility remains unmeasured.

### Canonical sources and what is actually selected

- Emmanuel Thomé, *Square Root Algorithms for the Number Field Sieve* (2012), introduction and
  §§1–4: <https://members.loria.fr/EThome/files/nfs-sqrt.pdf>. The introduction supplies the
  non-monic generator change and derivative correction; §1 gives direct lifting and coefficient
  bounds; §2 distinguishes classical Couveignes; §4 gives split-prime lifting and CRT with sign
  recombination. The complete paper was read through a text rendering.
- CADO-NFS source at `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`, inspected read-only:
  <https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sqrt>.
  `sqrt.cpp` is the scripted primary path. `FindSuitableModP` searches below 10^6 for an
  irreducible reduction, excluding primes dividing N or the leading coefficient. It uses
  extension-field Tonelli–Shanks and inverse-square-root Newton lifting. Its target precision is
  estimated, not certified; final square congruence is checked before gcd extraction. Do not copy
  unchecked inversions or infer an exact global-square certificate from that congruence alone.
- CADO's `crtalgsqrt.cpp` implements the published split-prime family separately, using
  meet-in-the-middle sign recombination. The ordinary Python driver does not automatically invoke
  it. Build configuration is conditional: `EXCLUDE_FROM_ALL` in a git checkout without
  `files.dist`, not an unconditional exclusion from release builds. Its README has a stale filename.
- Historical msieve mirror at `c8727d91305bdbe0972d160ef0ce61dd02ce9193` (2011), **not evidence
  about the latest msieve release**:
  <https://github.com/radii/msieve/tree/c8727d91305bdbe0972d160ef0ce61dd02ce9193/gnfs/sqrt>.
  `sqrt_a.c` describes and implements direct reciprocal Newton lifting, monic transformation,
  derivative correction, and evaluation. `sqrt.c` checks the final modular congruence. This version
  also has heuristic handling when no inert prime is found; that is not the selected fallback.

The direct method is the primary implementation: it is simpler than ideal-based reconstruction and
avoids a multi-prime sign problem, but alone excludes fields without inert primes. Thomé's method
closes that gap, at the cost of split-prime search and sign recombination. Full finite-field
polynomial factorization is **not** needed for the primary: polynomial arithmetic, gcd/powmod,
irreducibility testing, and a square root in an irreducible residue field suffice.

Classical Couveignes uses odd degree and inert primes; a known norm selects consistent signs across
primes. That norm is not the final map to Z/NZ. Montgomery/Nguyen is another established family,
requiring substantially more ideal/lattice machinery. Neither is enrolled as another implementation
variant in this campaign. The current code's method must not retain the Couveignes name.

### Confirmed seams beyond the original assessment

- `polyselect/mod.rs:219-261` already constructs the correct monic polynomial, whose generator is
  theta = c alpha. `sqrt/algebraic.rs:177-190` nevertheless multiplies a-b theta, not a-b alpha.
  Non-monic inputs therefore have an element mismatch before the final norm-map error.
- `factor_end_to_end_kat.rs:119-143` uses f=x²-2, m=6, N=35, although f(m)=34. It bypasses
  `PolyPair::verify`. Its gcd result is not evidence for the required polynomial homomorphism.
- `linalg/qc.rs:76-85` uses characters of homogeneous norms, rather than characters at individual
  roots. Once the recorded norm exponents have even sums, these columns retain at most norm-sign
  information; they do not distinguish the relevant unit square classes.
- `filter/merge.rs:249-284` replaces three incident rows by their total XOR. A shared column then
  remains set, and dependencies can be destroyed. Three copies of row [1] have a two-dimensional
  kernel, whereas the resulting single row [1] has none. This is not correct column elimination.
- `filter/matrix.rs:101-124` uses provenance union. Current destructive merges keep provenance
  disjoint, but correct weight-three elimination creates overlapping provenance. Both columns and
  provenance must then use symmetric difference.
- `sieve/mod.rs:314-346` verifies norm reconstruction but not the full ideal-label/index contract.
  `FactorBase::new` records bad-prime flags that the GF(2) builder does not use. Projective factors
  without affine roots are already unrepresented and their relations rejected.
- No full orchestrator exists. Public stages can be composed in an integration test. The square-root
  public functions have only GNFS assembly/test callers in the inspected workspace; rho uses NFS-DL,
  not these functions. Changing `PolyPair::number_field` globally would affect NFS-DL and is
  avoided.

### Common mathematical contract — derived

Use a polynomial f of degree d >= 1, leading coefficient c != 0, g=x-m, and N>1 with f(m)=0 mod N.
Normalize the sign of f locally if convenient; do not change the meaning of the caller's generator.
Require gcd(c,N)=1 for the unscaled map. A proper gcd is an incidental factor, not square-pipeline
success; gcd(c,N)=N is unsupported by this map. Validate these conditions before using `verify`,
whose current implementation is not safe for arbitrary N=0 inputs.

Define F(t)=c^(d-1) f(t/c), theta=c alpha, and A=Q[t]/(F). Require F squarefree exactly.
When f is irreducible, A is the number field K. Otherwise A is a product of number fields, not a
field. Exact quotient arithmetic and the following proof still apply; never claim that a split
prime proves irreducibility over Q. The genuine GNFS fixture must independently certify
irreducibility. An irreducible reduction at a good prime is a sufficient certificate.

For the nonempty lifted original relation set S, let k=|S| and h=ceil(k/2). Define:

    P     = product_S (a-bm) in Z
    gamma = product_S (a-b theta/c) in A
    eta   = c^(2h-k) product_S (ca-b theta) = c^(2h) gamma in Z[theta]
    H     = F'(theta)^2 eta in Z[theta]

All representatives are reduced to degree < d. The rational side requires P>0 and an exact square
root x0. Reject zero relation norms/products. For reducible F, require every relation factor to be
a unit in A; a zero divisor is unsupported and must be detected before prime search. This avoids
searching forever for a prime where an identically zero component becomes a unit.

If beta²=gamma, then w=c^h beta satisfies w²=eta and is integral, componentwise. The standard
derivative correction gives z=F'(theta)w in Z[theta]. Thus the reconstructed integer polynomial must
satisfy z²=H exactly modulo F. Recover beta=z/(c^h F'(theta)) over Q, reduce it, and independently
verify beta²=gamma. F' is invertible because F is squarefree; verify its inverse identity rather
than trusting `NumberFieldElement::inv` on arbitrary nonunits.

**Why derivative correction is sufficient.** Let theta_j be the d distinct complex roots of F and
w_j the corresponding values of w. Interpolation gives:

    z(t) = sum_j w_j F(t)/(t-theta_j)
    z_i  = Tr_A/Q(w h_i(theta))
    h_i(t) = sum_{ell=i+1}^d F_ell t^(ell-1-i)

Each w h_i(theta) is integral, so its trace is an integer. This proves integrality in the power
basis without assuming that Z[theta] is the maximal order or computing an integral basis.

**Certified coefficient bound.** Compute entirely with exact integers:

    R   = 1 + max_{i<d} |F_i|
    E   = sum_{i<d} |eta_i| R^i
    U   = ceil(sqrt(E))
    B_i = d U sum_{ell=i+1}^d |F_ell| R^(ell-1-i)
    B   = max_i B_i

Cauchy's bound gives |theta_j| <= R; hence |w_j| <= U. The trace formula gives |z_i| <= B_i for
every possible square root. A centered reconstruction modulus M **strictly greater than 2B** is
sufficient. Do not use floating point or the existing exact-square-only `isqrt` as if it returned a
floor square root for nonsquares. The bound is deliberately conservative: simpler certification
costs extra precision compared with certified numerical embedding bounds.

**Map.** Write beta=sum b_i theta^i with normalized rational b_i. Check each denominator is a unit
mod N, then evaluate at theta=cm to obtain y. This is a homomorphism on the localized order because
F(cm)=c^(d-1)f(m)=0 mod N. Therefore:

    y² = P mod N = (x0 mod N)².

Check that congruence before either gcd. There is no requirement to choose a positive real
embedding: beta and -beta give the same unordered gcd pair. Totally imaginary fields are allowed.

Checking actual normalized denominators is less restrictive than requiring gcd(disc(F),N)=1 or
blindly inverting F'(cm). A proper denominator gcd may reveal an incidental factor; gcd=N is an
unsupported mapping, not a nonsquare certificate. This contract does not claim that every possible
modular image is recoverable when this coefficientwise map fails.

The inspected msieve implementation instead scales both congruence sides, avoiding derivative
inversion. For even k, it uses X=c^(k/2)F'(cm)x0 and Y=z(cm). CADO normalizes back to the unscaled
image. The selected API exposes the unscaled contract above, with checked denominators and exact
identities. These are normalization choices within established reconstruction, not new algorithms.

### Reconstruction algorithms — selected design

**Direct Newton (primary).** Find an odd prime p, coprime to cN, for which F mod p is irreducible
and H mod p is nonzero. Use deterministic prime validation and an exact polynomial irreducibility
test. Compute s²=H in F_(p^d), returning NotSquare on a proved nonresidue, and set u=s^(-1).
In (Z/p^e Z)[t]/F, lift the inverse square root by:

    u <- u (3-Hu²)/2

Each step doubles precision. The invariant is Hu²=1 at the current precision. At M=p^e>2B, center
z=Hu mod (F,M) and check z²=H exactly. Only a global sign ambiguity exists at an inert prime; either
branch yields the same exact-root test. Failure at the proved bound certifies NotSquare. Failure to
find an inert prime within budget does not prove that none exists: dispatch the CRT fallback.

**Thomé CRT lifting (fallback).** Choose a bounded number ell>=1 of odd completely split primes,
coprime to cN, with distinct roots of F and nonzero H at every root. A nonresidue at any such root
is a valid NotSquare certificate. Use the existing simple-root Hensel lift to lift roots r_ij of F
to Q_i=p_i^e and then square roots s_ij of the fixed integer H(r_ij) modulo Q_i. Check both local
identities. Select e with M=product_i Q_i>2B. Taking ell=1 is an allowed published parameter choice,
not a special polynomial domain and not the primary algorithm.

For each local root, form its Lagrange polynomial L_ij modulo Q_i and the CRT idempotent weight
W_i=(M/Q_i) inverse(M/Q_i mod Q_i). The unknown global root has coefficients represented by:

    z = sum_{i,j} epsilon_ij W_i s_ij L_ij mod M, epsilon_ij in {-1,+1}.

Resolve signs using the published knapsack/meet-in-the-middle structure, with exact integer modular
interval tests. Fix one sign to remove the irrelevant global +/- symmetry. Split the remaining
signed contributions to one coefficient into two lists; match sums whose centered residue lies
within that coefficient's bound. Handle wraparound intervals exactly. For every match, reconstruct
all coefficients, check their bounds, and verify z²=H exactly. Do not accept a close leading
coefficient as proof. A wide interval may generate many matches; budget both lists and matches.
Exhausting all matching branches without a root certifies NotSquare; budget exhaustion does not.

This is Thomé §4, with exact verification and conservative bounds. It does not require floating
point, LLL, a new root-finding idea, or a full ideal factorization. Larger ell can reduce per-prime
precision but increases sign-search dimension. No production performance equivalence is claimed.

**Shared limits.** Bound prime search, coefficient/precision sizes, local-root search, and sign
recombination before expensive work. Fixed machine-prime range exhaustion is a resource outcome.
Use BigInt/BigUint for field cardinalities p^d and lifting moduli; never overflow a u64 exponent or
modular addition. Both methods return through the same exact verifier and modular mapper.

### Dependency constraints and supported collection boundary

Rational prime parity plus rational sign, on independently validated relations, gives the integer
square candidate P. Good affine ideal parity is also necessary for a field square. Even parity of
all principal-ideal valuations, when actually known, still need not give an element square: units,
real signatures, and nonprincipal ideal square roots can obstruct it. Sampled characters add
necessary conditions but do not certify completeness or choose reconstruction signs.

Use characters indexed by (q,s), with q odd, q not dividing c*disc(F), f(s)=0 mod q, and character
Legendre(a-bs,q). Select pairs for which all retained original relations evaluate nonzero; otherwise
skip that pair. Zero is not a valid character bit. A completely split q is sufficient but not
necessary. Compute merged-row characters by XOR over original provenance, not by norm characters.

Keep Relation's current full homogeneous-norm factor records, preserving NFS-DL consumers. In the
checked GNFS matrix path, classify p dividing c*disc(F) conservatively as exceptional and omit their
algebraic parity constraints rather than misrepresenting them as individual ideal valuations.
Retain zero columns if that preserves layout. Validate labels at good primes, including
a=br mod p. Keep rational parity/sign constraints for all rational factors.

This weakens the candidate filter, not successful-square verification. No complete exceptional or
projective ideal factorization is claimed. The collector still accepts only primitive affine pairs
whose two recorded norms factor completely over its represented bases; projective-only factors
remain outside collection. Do not silently change that boundary to obtain a fixture. A complete
projective/bad-ideal collector would need its own canonical design and fileset.

For non-monic inputs, Norm(gamma)=product F_hom(a,b)/c^k. Consequently, never add a check that the
raw homogeneous-norm product must be an integer square for every genuine field square. The exact
gamma test, not that shortcut, is authoritative.

Correct weight-three elimination chooses a pivot r0 and replaces r0,r1,r2 with r0 XOR r1 and
r0 XOR r2. The eliminated column vanishes from both rows, and their lift spans exactly the old
even-incidence combinations. XOR provenance permits overlap and cancellation. Row pruning is an
intentional loss of candidates; do not confuse it with the kernel-preservation obligation of merge.

Validate matrix shape, sorted unique indices, provenance bounds, per-row lifted parity, and kernel
membership before product construction. Deduplicate identical relation pairs in composition, without
adding a new special-q/lattice algorithm campaign. Try bounded deterministic combinations of
computed kernel vectors, not just basis vectors: two nonsquare candidates can combine to a square.

### Failure/API contract

Use typed, distinguishable outcomes for invalid polynomial/relations/matrix/dependency, rational
nonsquare, algebraic NotSquare, unsupported zero divisor or modular map, resource exhaustion,
verified root, valid-but-trivial congruence, and verified nontrivial factor. No unverified root can
be constructed as success. Invalid input/invariant failure is not silently treated as another
ordinary nonsquare candidate. Candidate exhaustion reports whether resource-limited attempts remain.

`rational_sqrt` and `algebraic_sqrt` become fallible; the algebraic result exposes an owned reduced
rational polynomial and its checked modular image, not a field element borrowing a temporary field.
`factor` returns a report/result distinguishing unsuccessful candidate classes instead of bare
Option. `factor_from_congruence` validates N and the congruence before gcd extraction. This is an
intentional API correction; do not preserve an unsafe panic wrapper to imitate compatibility.
Freeze concrete Rust names/signatures in the context unit before downstream implementation.

### Independently checkable fixtures — derived, not executed

- **Wrong map:** f=x²-2, N=7, m=3, gamma=9. beta=3 gives y²=2 mod 7; the old norm return
  y=2 gives y²=4. This is an algebra/map test, not a valid collected primitive relation.
- **Nonconstant root and assembly:** f=x²-2, N=527, m=23, relation (3,-2). gamma=3+2 alpha
  =(1+alpha)²; P=49, x0=7, y=24; gcd(7-24,527)=17 and gcd(7+24,527)=31. The negative b
  makes this an auxiliary seam fixture, not the positive-b sieve demonstration.
- **Non-monic, odd count:** f=2x²-1, N=7, m=2, relation (3,-4). theta=2 alpha, F=t²-2,
  gamma=3+2 theta, beta=1+theta. eta=12+8 theta, z=8+4 theta, y=5, y²=P=4 mod 7.
  P=11 is not an integer square: this tests algebraic reconstruction/map, not full assembly.
- **Root outside the order:** F=t²-8, eta=2, w=theta/2. Then z=F'w=8 and H=64. A search
  restricted to w in Z[theta] would falsely reject this square. Algebra-only fixture.
- **Local false positive:** F=t²-3, gamma=2+theta is a totally positive norm-one unit but not
  a square in Q(sqrt(3)). At q=23, roots 7 and 16 give values 9 and 18, both squares. Exact
  bounded reconstruction must reject it; at q=11 there is already a local nonresidue obstruction.
- **Sign ambiguity:** F=t²-2, gamma=9, p=7. In beta coordinates the local roots include
  +/-3 and +/-theta, but only +/-3 are global roots. Squaring locally cannot select signs.
- **No real embedding:** F=t²+1, beta=2+theta, gamma=3+4 theta. Exercise without Newton
  real-root finding or a real-sign convention.
- **No inert prime:** F=t⁴-10t²+1 is the minimal polynomial of sqrt(2)+sqrt(3), with Galois
  group V4. It has no inert primes. beta=1+theta+theta² squares to
  2 theta+13 theta²+2 theta³. Require CRT dispatch and exact reconstruction; do not replace this
  with an artificially exhausted primary search in a field that actually has inert primes.
- **Filter:** three identical one-column rows have three nonzero pair dependencies. The two
  pivot-relative rows must preserve them under provenance lifting, including overlap cancellation.

Also exercise precision just below/above the proved bound, exhausted budgets, invalid f(m), repeated
polynomial factors, zero character values, noninvertible mapping denominators, malformed provenance,
and trivial gcds. Root assertions allow global sign; factor assertions require a proper divisor
and the intermediate identities. Include method agreement on a field where both methods apply.

### Unselected material and remaining risks

The initially proposed single-split-prime enumeration is a specialization of Thomé §4. It is not
novel research and does not displace the primary algorithm. The trace bound is preserved above as
supporting derivation without a novelty claim. No genuinely new technique was established or
enrolled. Montgomery/Nguyen and other root variants remain optional future canonical coverage, not
prerequisites to completing this repair.

The substrate survey also found a possible bad-prime Dedekind-criterion defect: an integer
polynomial gcd is used where a finite-field gcd appears required. It was not independently
adjudicated here. The selected path does not rely on that criterion; record and scope its repair
separately rather than fixing all NFS-DL/ideal arithmetic while adding the needed finite-field
tools.

Runtime, fixture yield, useful dependency count, and performance of the certified bound remain
unverified. These are implementation/discovery risks, not excuses to relax the proof obligations.
