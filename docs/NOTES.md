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
