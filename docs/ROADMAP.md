# rDLP — Roadmap (Arc 2: Survey, Consolidation & Extension)

The durable, project-lifetime view of the work. Updated only at sub-track boundaries.
For the current sub-track's actionable detail, see `docs/PLAN.md`.

For the planning philosophy this document embodies, see
`~/.config/opencode/multisession/multi-session-planning.md` (the two-frame model).

**Arc 1 (the original implementation) is complete** — five tracks (ρ, G, D, E, S), six shared
crates, ~83 commit-shaped sessions, closed at `201df39` ("rDLP complete"). The arc-1 roadmap that
governed that work is preserved in git history at the commit prior to this document's introduction.
This is the **arc-2 roadmap**: a meta-campaign over the now-complete codebase — survey, revision,
consolidation, refactoring, and addendum/appendix work — plus a gated algorithmic-extension slot.

---

## Design intent

**Take the complete arc-1 rDLP — a self-consistent, pedagogically clear Rust reference library for
discrete-logarithm and integer-factorisation algorithms — and bring it to a finished, disciplined,
internally-coherent standard: aligned with its project seed, audited across docs/layout/testing,
equipped with a live correctness oracle, with its mathematical exposition reading as a continuous
textbook rather than a companion reference, and with its algorithmic spectrum confirmed complete (or
deliberately extended) across the integer / finite-field / elliptic-curve × classical / quantum
DLP matrix.**

**The artifact must stand on its own terms, not as the realization of an executed plan.** A reader
who never saw the roadmap should encounter a library organized by *topic* — the mathematics of
discrete logarithms — not by *track* and *phase*, the order in which the work was built. The git
history is where the construction narrative belongs; the code, docs, and names should cohere on the
subject's own principles. Planning-frame vocabulary (track/phase/sub-track ids, `◆` boundary marks,
letter-id labels) that leaked into the finished artifact is provenance residue to be re-anchored on
the mathematics — **except where a grouping genuinely coincides with a mathematical topic, in which
case the grouping survives and only the planning label changes** (e.g. "Track ρ" → "Pollard rho, the
birthday-bound baseline").

The arc-2 discipline is **survey-first**: the project's actual current shape is *discovered by audit*
before any consolidation or refactoring commits to scope. The anti-defocus anchor is this: arc 2 does
not add algorithmic mass for its own sake — every campaign serves either *coherence* (the library
reads and composes as one finished artifact) or *evidenced completeness* (a gap the survey actually
found, not a gap imagined at planning time). "While I'm here" rewrites and speculative extensions are
explicit risks, named and decided at sub-track boundaries, never absorbed mid-session.

Re-read this at every sub-track boundary. The question it answers: *is this work making the library
more finished and more coherent, or just bigger?*

---

## Sub-track DAG

Arc 2 is a survey substrate feeding five campaign slots. The substrate's only deliverable is a
findings ledger; every campaign is **gated behind survey findings** — the roadmap commits to the
*slot*, not to the specific changes, which are determined by what SURVEY discovers and fixed at
`/plan-shard` time per campaign.

```
                          ┌─────────────────────────────────────────┐
                          │  A. SURVEY  (substrate, findings-only)   │
                          │  audits all dimensions → Discoveries log │
                          └────────────────────┬────────────────────┘
                                               │ (every campaign depends on SURVEY)
        ┌──────────────┬─────────────────┬─────┴───────────┬──────────────────┐
        ▼              ▼                 ▼                 ▼                  ▼
  ┌───────────┐  ┌───────────┐    ┌───────────┐    ┌────────────┐    ┌──────────────┐
  │ B. ALIGN  │  │ C. ORACLE │    │ D.REFACTOR│───▶│E.CONSOLIDATE│   │  F. EXTEND   │
  │ template  │  │ CADO-NFS  │    │ crate     │    │ docs-layer  │   │ spectrum     │
  │ ↔ rDLP   │  │ sidecar   │    │ layout +  │    │ + math      │   │ additions    │
  │formalities│  │ build     │    │ dedup     │    │ continuity  │   │ (content TBD) │
  └───────────┘  └───────────┘    └─────┬─────┘    └──────┬──────┘   └──────┬───────┘
                                        │  (layout settles │ before docs    │
                                        │   reference it)   │ reference code │
                                        └──────────────────┘                │
                                                          (REFACTOR before ──┘
                                                           EXTEND: new code
                                                           sits in final layout)
```

**Edges (dependencies):**

| Sub-track | Depends on | Why the edge |
|-----------|-----------|--------------|
| A. SURVEY | — (substrate) | The root. Nothing commits scope before audit. |
| B. ALIGN | SURVEY | Backport direction is a survey finding (what's missing where). |
| C. ORACLE | SURVEY | Sidecar design depends on the survey's CADO-usage + testing-balance findings. |
| D. REFACTOR | SURVEY | Layout changes are gated on the layout/dedup audit. |
| E. CONSOLIDATE | SURVEY, **REFACTOR** | Docs reference code structure; let layout settle before rewriting the doc-layer + math continuity to point at it. |
| F. EXTEND | SURVEY, **REFACTOR** | New algorithmic code must sit in the final crate layout; its *content* is gated on the spectrum-gap audit + 4-reference distillation. |

**Per-sub-track charge and category mix:**

- **A. SURVEY** *(substrate — the irreducible foundational sub-track)*. Charge: audit every
  dimension and emit a findings/recommendations ledger; **no design, no code, no refactor**.
  Dimensions: (1) rust-template alignment, both directions; (2) docs-layer discipline (inline /
  human-code-ref / agent-docs reference directions + allowances); (3) project layout (crate peering,
  organisation, duplicative code); (4) testing balance, coverage, unit/integration mix, over-long
  tests — **and the meta-question of whether full-coverage is even the right target for
  toy-scale-but-compute-heavy pedagogical code**; (5) CADO-NFS sidecar usage and build needs;
  (6) math-exposition continuity (textbook-vs-reference, dir/file/notation/voice/audience
  consistency, code↔math cross-reference health); (7) algorithmic-spectrum completeness across the
  integer/field/EC × classical/quantum matrix; (8) distillation of the 4 quantum-DLP references into
  candidate inclusions; (9) **provenance leakage & topic-vs-track organization** — catalog every
  planning-frame token (track/phase/sub-track ids, `◆` marks, letter-id labels) at every depth (prose
  docs, code identifiers, module/file/dir names, benchmark labels, test names) and classify each as
  pure residue (re-anchor on topic) or grouping-that-coincides-with-a-topic (keep grouping, change
  label). Each dimension's findings land as Discoveries-log entries that shape the campaigns.
  Category: **substrate (audit)** — high-judgment, wide design surface, front-loaded.

- **B. ALIGN**. Charge: reconcile rDLP with its `~/Source/rust-template` seed in **both
  directions** — backport template formalities rDLP lacks (per survey: `deny.toml`, `rustfmt.toml`,
  `rust-toolchain.toml`, `LICENSE`, the `[lints]` table, a coverage gate, `development.md`), and
  forward rDLP-originated discipline worth seeding back into the template (multisession docs
  layering, the dev-oracle policy, the docs-register contracts). Category: **integrative / substrate
  (formalities)**. Largely independent of the other campaigns.

- **C. ORACLE**. Charge: specify and build the **dynamic CADO-NFS sidecar** — automated build for
  regression and comparison testing, available both automatically (CI/test harness) and on-demand
  (user-invoked). Honors the arc-1 dev-oracle policy (CADO as opt-in validation sidecar, never part
  of how rDLP computes). Category: **substrate (test infrastructure) + integrative**. Largely
  independent.

- **D. REFACTOR**. Charge: restructure crate peering/organisation and eliminate duplicative code,
  **strictly gated on the SURVEY layout audit** — no speculative restructuring. **Also owns the
  code-depth half of the de-provenancing (D9):** renaming planning-frame code identifiers, module/
  file/dir names, benchmark labels, and test names onto topic terms, where SURVEY classified them as
  residue (a rename is a layout change with compiler blast radius, so it lives here). Category:
  **substrate (layout)**. Must precede CONSOLIDATE and EXTEND so docs and new code target the final
  shape.

- **E. CONSOLIDATE**. Charge: the exposition campaign — (a) enforce the three-docs-layer discipline
  (inline / human-code-ref / agent-docs) with correct reference directions, (b) rewrite the
  mathematical exposition so it **reads as a continuous textbook with a suggested path**, internally
  consistent in dir/file names, notation, audience, voice, thoroughness — readable independently yet
  cross-referring to code where useful, and (c) **own the prose-depth half of the de-provenancing
  (D9):** re-anchor planning-frame vocabulary in all human-facing prose (PEDAGOGY/MATHEMATICS/README)
  onto topic terms, so the docs read as a treatment of the mathematics, not a record of an executed
  plan. Category: **integrative**. Depends on REFACTOR (docs point at final layout *and* at the
  already-renamed identifiers).

- **F. EXTEND** *(content gated behind survey findings)*. Charge: add the algorithmic components the
  SURVEY's spectrum-gap audit + 4-reference distillation actually justify — across integer / finite-
  field / elliptic-curve, classical / quantum. **The roadmap commits to this slot, not to its
  algorithms.** Possible content surfaced by intent (to be confirmed/refuted by survey): additional
  DLP subject areas; newer quantum resource-estimate methods; reference-derived techniques. Category:
  **algorithm + optimization** (decomposes per arc-1's substrate/algorithm/opt pattern once content
  is known). Depends on REFACTOR.

---

## Cross-track contracts (coarse)

Sketched at roadmap grain; frozen precisely at `/plan-shard` time per campaign and at inflection
junctures.

| Contract | Flavour | Spans | Sketch |
|----------|---------|-------|--------|
| **C-Findings** | prose | SURVEY → all campaigns | The survey findings ledger is the authoritative scope source for every campaign. A campaign that exceeds its findings is defocus; one that ignores a finding is rigidity. Frozen at SURVEY close. |
| **C-Testing-Philosophy** | prose | SURVEY → ALIGN, CONSOLIDATE, F | The resolved answer to "what is the right testing target for toy-scale compute-heavy KAT-driven pedagogical code" (coverage target, unit/integration balance, test-length norms). Frozen at SURVEY close; ALIGN's coverage-gate decision and any new tests honor it. |
| **C-Layout** | compiler | REFACTOR → CONSOLIDATE, EXTEND | The final crate-peering/module structure. Once REFACTOR freezes it, docs reference it and new code sits in it. Compiler-enforced (workspace members, module paths). |
| **C-DocsLayer** | prose | CONSOLIDATE → all future docs | The three-layer reference-direction discipline (inline ↔ human-code-ref ↔ agent-docs: who may reference whom). Frozen at CONSOLIDATE; governs all subsequent doc edits including the math textbook. |
| **C-MathSpine** | prose | CONSOLIDATE → EXTEND | The textbook's suggested-path spine, notation, voice, audience. Any EXTEND chapter must slot into this spine, not append incoherently. Frozen at CONSOLIDATE. |
| **C-Coherence** | prose | SURVEY → REFACTOR, CONSOLIDATE, EXTEND | The artifact-stands-on-its-own-terms principle: the finished library is organized and named by mathematical topic, not by planning track/phase; provenance residue is re-anchored, groupings that coincide with a topic are kept under topic labels. Frozen at SURVEY (the residue catalog + the keep/rename classification); enforced by REFACTOR (identifiers) and CONSOLIDATE (prose), and binding on all EXTEND naming. |
| **C-Oracle** | test | ORACLE → regression suite | The CADO-NFS sidecar build + invocation interface (automatic + on-demand). Downstream regression/comparison tests consume it. Frozen at ORACLE close. |
| **C-TemplateSeed** | prose | ALIGN ↔ rust-template | The agreed bidirectional formalities set: what flows template→rDLP and what flows rDLP→template. Recorded so future template-derived projects inherit the reconciled discipline. |

---

## Sequencing & scope

**Proposed order:** `A. SURVEY` → then `B. ALIGN` ∥ `C. ORACLE` ∥ `D. REFACTOR` (independent of each
other) → `E. CONSOLIDATE` (after REFACTOR) → `F. EXTEND` (after REFACTOR; content gated on SURVEY).

Reasoning:
- **SURVEY first, always.** It is the substrate; nothing else has its scope until it runs. This is the
  whole point of the survey-first discipline.
- **ALIGN / ORACLE / REFACTOR are mutually independent** after SURVEY and may run in any order or be
  interleaved by fatigue/appetite. REFACTOR should land before CONSOLIDATE and EXTEND.
- **CONSOLIDATE after REFACTOR** so the docs (and especially the math textbook's code cross-references)
  target the final layout rather than chasing a moving structure.
- **EXTEND last and conditional.** If SURVEY finds the spectrum complete, EXTEND may reduce to a
  recommendations note (a deliberate non-extension) rather than new tracks. That is a *successful*
  outcome, not a gap.

**Cost-of-wrong texture (informs commit size + `juncture-tier` at shard time):**

| Sub-track | Correctness-criticality | Design-error cost | Leaning |
|-----------|------------------------|-------------------|---------|
| A. SURVEY | low (findings, no code) — but **high judgment**; a missed finding mis-scopes everything downstream | high (mis-survey propagates) | **opus juncture** at close; the findings ledger is the highest-leverage artifact in arc 2 |
| B. ALIGN | medium (formalities, CI gates) | low (reversible config) | sonnet juncture; strong template reference reduces design risk |
| C. ORACLE | medium (test infra correctness) | medium (sidecar build is fiddly) | sonnet juncture; opus only if build automation proves intricate |
| D. REFACTOR | **high** (C-Layout is compiler-enforced; rework propagates) | **high** (wrong layout invalidates CONSOLIDATE + EXTEND) | **opus juncture**; this is the arc's rigidity/blast-radius node |
| E. CONSOLIDATE | medium (prose; the inner test loop doesn't guard prose) | medium (incoherent textbook is costly to re-thread) | opus juncture at the math-spine freeze (C-MathSpine) |
| F. EXTEND | high (new algorithmic correctness — KAT-gated) | high if content is non-obvious | **opus** for any substrate session; per arc-1's Category-A discipline |

**Commit-size leaning:** arc-1's codebase is a *garden* (coherent structure, strong KAT discipline,
trustworthy tests) — lever 5 (inner-loop bandwidth) is strong, which licenses smaller, denser commits
and *opting juncture-tier down to sonnet* where correctness-criticality is also low (ALIGN, ORACLE).
The two opus-leaning nodes (REFACTOR, F-EXTEND-substrate) are where blast radius and abductive design
weight concentrate — hold the default there.

---

## Discoveries & open questions

Seeded at construction; the reconcile step and `/plan-shard` boundary reconciliations append here as
arc 2 runs.

### Construction-time seeds (from intent + substrate read)

- **D1 — Template-formalities gap (confirmed, partial).** rDLP workspace root lacks `deny.toml`,
  `rustfmt.toml`, `rust-toolchain.toml`, `LICENSE`, the `[lints]` table, a coverage gate, and a
  `development.md` — all present in `rust-template`. SURVEY must confirm scope (per-crate vs
  workspace-level lints) and decide each item's flow direction. *Feeds ALIGN.*
- **D2 — Testing-philosophy meta-question (open).** Is full/100% line coverage the right target for a
  toy-scale-but-compute-heavy, KAT-driven pedagogical math library? Or does KAT + property-test
  coverage of *mathematical* behavior matter more than line coverage of compute kernels? Resolution
  becomes C-Testing-Philosophy. *Feeds ALIGN, CONSOLIDATE, F.* (This is a named design tension, not a
  rhetorical aside.)
- **D3 — Docs-layer reference directions (open).** The three layers (inline / human-code-ref
  `PEDAGOGY.md` / agent-docs) have prescribed reference directions and allowances; SURVEY audits
  whether the current ~257KB MATHEMATICS.md + per-crate PEDAGOGY.md + README set honors them.
  *Feeds CONSOLIDATE → C-DocsLayer.*
- **D4 — Math-exposition continuity (open).** MATHEMATICS.md is a 12-chapter survey assembled
  chapter-by-chapter at track boundaries (T.0/T.G/T.D/T.E/T.S/T.Z). Does it read as a *continuous
  textbook with a suggested path*, or as a stitched companion reference? Audit notation/voice/
  audience/dir-file-name consistency and code↔math cross-reference health. *Feeds CONSOLIDATE →
  C-MathSpine.*
- **D5 — Crate layout / dedup (open).** Nine workspace members; `rho` hosts both Track ρ and the
  large Track E (10+ attack modules). Is `rho` overloaded? Should Track E peer out? Is there
  duplicative code across `shared/*`? *Feeds REFACTOR → C-Layout.*
- **D6 — CADO-NFS sidecar dynamic build (open, design-heavy).** Specify automated + on-demand build
  for regression/comparison testing, honoring the opt-in-oracle policy. Open: build trigger model
  (CI vs lazy), version pinning, what comparisons it gates. *Feeds ORACLE → C-Oracle.*
- **D7 — Algorithmic-spectrum completeness (open, abductive).** The intent's discriminant: full
  spectrum from integer DLP → finite-field DLP → EC-over-finite-field DLP, classical and quantum,
  constrained to mathematically-necessary algorithms (excluding hardware/distributed specifics). Are
  there omitted DLP subject areas? SURVEY maps the current matrix against this charge. *Feeds F.*
- **D8 — Four quantum-DLP references, undistilled (open).** To read in SURVEY:
  (1) Google QuantumAI cryptocurrency whitepaper
  `https://quantumai.google/static/site-assets/downloads/cryptocurrency-whitepaper.pdf`;
  (2) `arxiv.org/abs/2603.28627`; (3) `arxiv.org/abs/2606.02235`;
  (4) `github.com/ecdsafail/ecdsafail-challenge`. Distill what belongs in a pedagogical/reference DLP
  library across its domains and strategize inclusion. *Feeds F (content), possibly CONSOLIDATE
  (new exposition).* **Note:** arxiv IDs `2603.*`/`2606.*` are dated beyond construction time
  (2026-06) — SURVEY must verify the identifiers resolve; if not, surface for correction.
- **D9 — Provenance leakage / topic-vs-track organization (open, coherence).** The arc-1 codebase and
  docs bear the fingerprints of their own construction: "Track ρ/G/D/E/S", "Phase α…ζ", sub-track
  letter-ids (E.W, G.A.3), `◆` boundary marks, and planning vocabulary appear in prose
  (PEDAGOGY/MATHEMATICS/README), and likely in code identifiers, module/file/dir names, benchmark
  labels, and test names. The artifact should read as a DLP library organized by *topic*, not as the
  realization of an executed plan — the git log carries the construction narrative; the code coheres
  on its own terms. SURVEY catalogs every occurrence at every depth and classifies each as **pure
  residue** (re-anchor on the mathematics) or **grouping-coincides-with-topic** (keep the grouping,
  change only the planning label — e.g. the five "tracks" map onto five real mathematical families, so
  the grouping survives as "Pollard rho / GNFS / NFS-DL / algebraic-ECDLP attacks / Shor"). Split
  execution: *prose fixes → CONSOLIDATE; identifier/file/dir/bench/test renames → REFACTOR* (a rename
  carries compiler blast radius). *Feeds REFACTOR + CONSOLIDATE → C-Coherence; binding on EXTEND
  naming.*

### Open questions for the human — RESOLVED (ratified 2026-06-21 at D. REFACTOR shard)

All four open-Qs surfaced at SURVEY close were ratified by the human at the D. REFACTOR shard,
accepting the A.7 adjudicator recommendations:

- **Q1 — `rho`-hosts-Track-E overload (D5): RATIFIED keep-with-baseline (no crate split).** The 8
  Track-E attack modules stay co-located with the Pollard-rho baseline. Rationale: the E.W
  cross-attack bench measures the attacks against the rho baseline; shared `rho::curve`/`rho::field`
  types; "attacks live with the baseline they're measured against" pedagogy is load-bearing. Binds
  C-Layout (frozen at REFACTOR D.6).
- **Q2 — Coverage gate value (D2): RATIFIED 80% floor, ALIGN re-tunes to measured baseline.** The
  doctrine shape (math-behavior KAT primary; line coverage a dead-code floor) is frozen
  (C-Testing-Philosophy); ALIGN calibrates the number against a real `cargo llvm-cov` baseline.
- **Q3 — F-EXTEND scope ceiling (D7): RATIFIED defer to arc-3.** The one genuine gap (hyperelliptic
  Jacobian DLP via Gaudry index calculus) is real but self-contained; arc-2's spine is
  consolidation/alignment. F-EXTEND in arc-2 reduces to a recommendations note (a deliberate
  non-extension); the gap becomes an arc-3 roadmap candidate.
- **Q4 — C-Coherence borderline default (D9): RATIFIED preserve-under-topic-label.** Coarse
  groupings that map onto real mathematical families (Track-X labels; coarse G.*/S.* pipeline/Shor
  stages) keep the grouping under a topic-native label; only fine-grained S.X.Y / G.X.Y session IDs
  are pure residue and dissolve. Binds REFACTOR (code-half, D.4/D.5) and CONSOLIDATE (prose-half).

### D. REFACTOR discoveries (folded up at D.6 ◆, 2026-06-21)

- **D10 — `rho/src/field/` had 3 files, not 1 (D.1 discovery).** `rho/src/field/mod.rs` was the
  25-line re-export wrapper; `monty.rs` and `naive.rs` were dead code never declared in `mod.rs`.
  All three were removed with the wrapper collapse. No layout-premise change — the extras are
  mechanical consequences of the collapse.
- **D11 — `shor/Cargo.toml` required a `shared-field` dev-dep (D.1 collapse-forced).** Before the
  collapse, `shor/tests/ecdlp_kat.rs` imported `rho::field::FpMonty` (the re-export wrapper). Once
  `rho::field` was removed, the compiler forced a direct `shared-field` dev-dep in `shor/Cargo.toml`.
  This is a structural necessity of the collapse, not a REFACTOR-scope violation. Recorded as a
  mechanical consequence of C-Layout.
- **D12 — C-Layout is now FROZEN (compiler-enforced, D.6 ◆).** Final crate layout: `rho/src/field/`
  (`mod.rs`, `monty.rs`, `naive.rs`) and `rho/src/util/{batch_inv,mp,mod}.rs` removed; `rho` imports
  `shared_field::{Fp, FpNaive4, FpMonty4}` and `shared_bigint::batch_invert` directly; `F: Fp<4>`
  trait bounds preserved; no crate split (open-Q 1 ratified); no other `shared/*` dedup (F-D5-03).
  Verified: `cargo check --workspace` and `cargo test --workspace` green at D.5 terminus (`3ef9aca`).
  Binds CONSOLIDATE (docs reference final layout) and EXTEND (new code sits in it).
- **D13 — C-Coherence code-half is FROZEN (D.6 ◆, with ~1% residue for CONSOLIDATE).** All ~524
  cataloged F-D9-01…06 code-depth tokens re-anchored across D.2–D.5 (`89f1e22`, `84b993c`,
  `f2302a4`, `3ef9aca`). Completeness ≈ 99%: five residue tokens escaped the sweep —
  `rho/src/curve/test_curves.rs:36` (E.A.2), `:139` (E.A), `rho/benches/attacks.rs:20,26`
  (E.W table / E.W.2 chapter), `rho/tests/hyperelliptic_kat.rs:6` (## E.I.2). All behaviourally
  inert; recommended for CONSOLIDATE (which will be in these files' neighbourhood). The three
  load-bearing carve-outs held: algorithmic `Phase 1/2/3` step-indices in `rho/src/ecdlp/walk.rs`
  and `rho/src/index_calculus/solve.rs` preserved; quantum-phase-estimation `Phase N/M` tokens in
  `shor/tests/factor_kat.rs` preserved. KAT suite identically green throughout (inertness invariant).
- **D14 — `c_ek_relation_round_trip` renamed to `relation_round_trip` (D.3 scope).** The single
  contract-label prefix (`c_ek_`) was removed from the test identifier in `rho/tests/
  index_calculus_kat.rs`, within D.3's sub-track-ID + contract-label re-anchoring scope. Compiler-
  and test-verified (the renamed test still passes). This is the only identifier rename in REFACTOR
  (the `sub_track_close_curve_axioms_intact` → `binary_curve_axioms_intact` rename was also in D.3).

---

### B. ALIGN discoveries (folded up at B.7 ◆, 2026-09-22)

- **D15 — Seven planned rows closed plus one inserted row; two additive-reshards, all
  human-authorized (B.7 ◆).** The rDLP portion ran readiness baseline (B.1), licensing (B.2), an
  inserted advisory-remediation row (B.3a), dependency-audit policy (B.3), root-toolchain/formatting
  (B.4), workspace lints (B.5), measured coverage (B.6), and the contributor-guide/template handoff
  (B.7). B.3a and the B.4/B.5 additive-reshards were each surfaced by `@plan-juncture`
  discovery-adjudication and explicitly authorized by the human, never applied silently.
- **D16 — Five contracts frozen in FULL: C-Scope, C-Baseline, C-License, C-Audit, C-Coverage
  (B.7 ◆).** The adopted, genuinely-green gates are `cargo check --workspace`,
  `cargo test --workspace`, `cargo dep-policy` (= `cargo deny check`), `cargo audit`, and
  `cargo coverage` (measured 92.78% line coverage against the ratified 80% floor).
- **D17 — Two contracts frozen SCOPED, not FULL: C-Toolchain and C-Lints (B.7 ◆).** Both declare
  and verify their policy while leaving workspace-wide conformance unadopted. C-Toolchain: one exact
  root toolchain (1.95.0) with required components (`rustfmt`, `clippy`, `llvm-tools-preview`), no
  member-level overrides, and a declared 100-column rustfmt policy are frozen — `cargo fmt-check` is
  red across 169 files workspace-wide, so that gate is not adopted. C-Lints: a single-source
  workspace lint policy (`unsafe_code` forbidden, `missing_docs` warn, `clippy::all` denied,
  `clippy::pedantic` warned with explicit group priorities) with universal member inheritance is
  frozen — `cargo lint` is red (234 measured `clippy::all`-group hits across 25 distinct lint kinds,
  plus a larger pre-existing `pedantic`/`missing_docs` backlog), so that gate is not adopted either.
  Both scoped freezes execute this campaign's own standing rule — recorded below in "ALIGN boundary
  decisions" and in the PLAN's Discoveries log — that measured remediation exceeding a config-only
  row's fileset is sharded separately, not hidden inside the row.
- **D18 — C-TemplateSeed frozen on the rDLP side only (B.7 ◆).** The root-toolchain
  pin-with-components pattern, a cargo-alias-naming-shadowing hazard (naming a `cargo-deny`-check
  alias literally `audit` would silently shadow the real `cargo-audit` binary), and the
  declared-but-unadopted-gate pattern are recorded in `docs/NOTES.md` as concrete rust-template
  targets. No `rust-template` file was read or changed this sub-track; the bidirectional contract
  stays OPEN.
- **D19 — Two sized but unscheduled follow-on items, neither owned by any plan yet.** A
  whole-workspace reformat (needs a `#[rustfmt::skip]` policy decision for numeric constant tables,
  plus resync of the exact source-line citations in `shared/numth/docs/PEDAGOGY.md` and this
  document — check whether those citations are already stale from REFACTOR before assuming the
  reformat is the only desync source); and a lint-remediation pass (the `clippy::all` floor plus the
  pedantic/`missing_docs` mass, including a suppression-policy decision on `many_single_char_names`
  against this project's register discipline for a mathematics library). Both are design-bearing,
  not mechanical. Whichever plan picks these up owns `fmt-check`/`lint` gate adoption respectively.
- **D20 — A durable-to-rolling reference newly introduced this sub-track, flagged for CONSOLIDATE
  (C-DocsLayer not yet frozen).** `.cargo/config.toml`'s comments and `docs/development.md` both
  point at `docs/NOTES.md` for measurement provenance — in tension with the three-layer reference
  discipline ALIGN itself just recorded as a rust-template seed target (durable prose should not
  point at rolling notes). This may or may not be a defect — pointing at where measurement
  provenance lives could be preferable to duplicating it — but CONSOLIDATE should adjudicate it as
  newly-authored this week, not misread it as pre-existing residue from before REFACTOR.

### ALIGN boundary decisions (2026-09-17)

- The rDLP portion is seven sessions: readiness baseline, licensing, dependency audit, root
  toolchain/formatting, workspace lints, measured coverage, and contributor-guide/template handoff.
- One exact stable root toolchain replaces the member-local floating-stable file. Select and verify
  the release during baseline work; forward-seed the root-toolchain policy into rust-template.
- Workspace-wide gate adoption may expose source remediation. Measure first and shard actual
  remediation separately; do not hide it inside a config-only row.
- B.7 closes only the rDLP portion. Bidirectional ALIGN remains open until template changes land.

## Status

All not-started at construction. `/plan-shard` marks in-progress; the reconcile step marks done.

| Sub-track | Status |
|-----------|--------|
| A. SURVEY | **done** (A.7 ◆ closed 2026-06-21; findings ledger frozen in docs/SURVEY.md; contracts frozen in docs/PLAN.md) |
| B. ALIGN | **rDLP portion done** (B.7 ◆ closed 2026-09-22; bidirectional ALIGN remains open — rust-template forward-seeds need their own execution plan) |
| C. ORACLE | not-started |
| D. REFACTOR | **done** (D.6 ◆ closed 2026-06-21; C-Layout + C-Coherence code-half frozen; discoveries folded up) |
| E. CONSOLIDATE | not-started |
| F. EXTEND | not-started (content gated on SURVEY; scope-ceiling open-Q for human) |

---

## Updates to this document

This document is rewritten **only at sub-track boundaries**. Day-to-day session work is captured in
`docs/PLAN.md`. Discoveries that affect this document are queued in the Discoveries log above and
integrated at the next sub-track boundary. A discovery severe enough to require immediate roadmap
revision triggers an inflection-point Opus juncture, not an ad-hoc edit by a Sonnet session.

- **2026-06-17 — Arc-2 construction.** This roadmap replaces the arc-1 roadmap (preserved in git
  history). Constructed via `/roadmap-construct` from the four-cluster intent seed. Six sub-tracks:
  one survey substrate + five gated campaign slots. Survey-first discipline; extension gated behind
  findings.
- **2026-06-17 — Coherence audit added (D9 / C-Coherence).** Folded in a ninth SURVEY dimension —
  provenance leakage & topic-vs-track organization — so the finished artifact stands on the
  mathematics' own terms rather than as the realization of an executed plan. Split execution across
  REFACTOR (identifiers/files/benches) and CONSOLIDATE (prose); no new sub-track. Design intent
  amended to state the standalone-artifact principle.
- **2026-06-20 — A. SURVEY sharded (sub-track-boundary reconciliation).** First arc-2 sub-track
  opened. `/plan-shard` from `@architect` decomposed the survey substrate into 7 findings-only
  sessions (6 dimension audits + 1 integrative ◆ ledger-freeze), the nine audit dimensions D1–D9
  grouped by shared substrate (D3/D4/D9-prose → the doc corpus; D5/D9-code → the crate tree).
  `juncture-tier: opus`, opt-down unavailable (lever-5 weak: no test inner-loop catches a wrong
  finding; the ◆ findings-ledger review is the only correctness gate). No prior arc-2 sub-track to
  mark done; no discoveries folded up (arc-1's completion already recorded at construction). The
  stale arc-1 *terminus* PLAN.md (Z.1+T.Z, done at `201df39`) is superseded by the SURVEY PLAN.
- **2026-06-21 — A. SURVEY closed at A.7 ◆ (findings ledger frozen).** All nine audit dimensions
  confirmed audited across A.1–A.6. Key findings folded up from `docs/SURVEY.md` and `docs/PLAN.md`:
  - **D1 (ALIGN scope):** 8 absent root formalities items confirmed (F-D1-01…09); 3 rDLP→template
    forward-seeds recorded (F-D1-10…12). All 12 findings feed ALIGN.
  - **D2 (C-Testing-Philosophy frozen):** math-behavior KAT coverage is the primary correctness
    target; 80% line gate as a dead-code floor (gate value open-Q for human); two-tier test pattern
    (inline substrate + external KAT); test-length norm (one file per algorithm family, however
    long); oracle-gate norm (`#[ignore]` for external-tool and compute-heavy tests).
  - **D5 (REFACTOR scope):** `rho` hosts 8 Track-E attack modules alongside 5 Track-ρ modules —
    layout overload, open-Q for human (peer out vs keep-with-baseline; adjudicator recommendation:
    keep-with-baseline). `rho/src/field/` and `rho/src/util/` are thin re-export wrappers with
    duplicated `batch_invert` tests — dedup candidate. No other `shared/*` dedup found.
  - **D6 (ORACLE scope):** sidecar trigger = lazy/on-demand + opt-in CI flag; pin = `git-2.0.1`
    (`88c4751`) on GitHub mirror; 4 tolerance-bounded comparisons (existing stubs); policy fidelity
    confirmed (validation-only, never on green path).
  - **D3/D4 (CONSOLIDATE scope):** 3 docs-layer defects found (bare Rust paths in MATHEMATICS.md;
    stale forward-ref in shor/PEDAGOGY.md; missing maths-first header in numfield/PEDAGOGY.md).
    Math-spine: chapters 8+9 heading inconsistency (`#` vs `##`); ToC complete; voice/notation
    consistent. One citation-style inconsistency.
  - **D7 (F-EXTEND scope):** exactly one genuine spectrum gap — **hyperelliptic Jacobian DLP
    (Gaudry index calculus)**, the deferred downstream solve of the GHS transfer (Ch. 10.3).
    Substrate exists (`rho/src/hyperelliptic/`: Mumford divisors, Cantor group law); missing layer
    is the index calculus on the Jacobian + KATs + MATHEMATICS.md treatment. All other absent
    algorithms (trial division, QS, BSGS, classical index calculus for 𝔽ₚ*) are deliberate
    non-extensions. EXTEND scope-ceiling is an open-Q for human (adjudicator recommendation:
    defer to arc-3).
  - **D8 (CONSOLIDATE citation updates):** all 4 ROADMAP-listed references resolve (risk flag not
    triggered). None adds a new mathematical algorithm. Three 2026 papers (Google QuantumAI
    whitepaper, Cain et al., Schrottenloher) feed CONSOLIDATE citation updates in §11.4.4 and
    §11.5.1; they do not scope F-EXTEND algorithmic work.
  - **D9 (C-Coherence frozen):** code-half: doc-comment-heavy, identifier-light (exactly 1
    identifier hit: `sub_track_close_curve_axioms_intact`). Prose-half: ~185 tokens; ~110
    grouping-coincides-with-topic (Track X labels); ~75 pure residue (Phase N, session IDs, ◆
    marks). S.A/S.B/S.C borderline — adjudicator recommendation: preserve-under-topic-label.
  - **Four human open-Qs surfaced:** (1) `rho`-overload split (adjudicator: keep-with-baseline);
    (2) coverage gate value (adjudicator: 80% floor, ALIGN re-tunes to measured baseline);
    (3) F-EXTEND scope ceiling (adjudicator: defer to arc-3); (4) C-Coherence borderline default
    (adjudicator: preserve-under-topic-label). These ride to human review; they do not block the
    freeze.
  - **Scope-routing invariant (load-bearing):** inline doc-comments + identifier/file/dir/bench/
    test renames → REFACTOR; human prose → CONSOLIDATE; config/manifest/CI → ALIGN; CADO sidecar
    → ORACLE; new algorithm/chapter → EXTEND. A finding routed to the wrong campaign is the
    primary cross-campaign defocus mode.
  A. SURVEY status: **done**. B–F campaigns: not-started (gated on SURVEY findings now frozen).
- **2026-06-21 — D. REFACTOR sharded (sub-track-boundary reconciliation).** Second arc-2 sub-track
  opened; first *campaign* (SURVEY was the substrate). `/plan-shard` from `@architect` decomposed
  REFACTOR into 6 sessions: D.1 (Opus) collapses the contrived `rho` re-export wrappers + freezes
  C-Layout (compiler-enforced); D.2–D.5 (Sonnet) execute the code-depth de-provenancing
  (F-D9-01…06, ~524 doc-comment tokens) by crate weight (rho 55% / gnfs 21% / shor+shared 24%),
  rho split at the Phase-N-vs-sub-track-ID token-class seam; D.6 ◆ (Opus @architect) freezes
  C-Layout + C-Coherence code-half and folds up. `juncture-tier: opus` **held** despite the lever-5
  opt-down being available (strong test loop + low criticality, the first arc-2 sub-track where it
  is) — honoring REFACTOR's pre-committed blast-radius-node status; the D.1 wrapper-collapse proved
  a genuine structural change. All four SURVEY open-Qs ratified (above). Design finding folded up:
  the `rho` wrappers shadowed an existing upstream alias (`shared_field::FpMonty4`) and preserved
  no parameterization the generic `Fp<L>` trait does not already give — the collapse removes
  contrivance, not capability. D. REFACTOR status: **in progress** (at shard time; closed at D.6 ◆, see below).
- **2026-06-21 — D. REFACTOR closed at D.6 ◆ (C-Layout + C-Coherence code-half frozen).** All six
  sessions complete. C-Layout FROZEN (compiler-enforced): `rho/src/field/` and `rho/src/util/
  {batch_inv,mp,mod}.rs` removed; `rho` imports `shared_field::{Fp, FpNaive4, FpMonty4}` and
  `shared_bigint::batch_invert` directly; `F: Fp<4>` bounds preserved; no crate split; no other
  `shared/*` dedup. C-Coherence code-half FROZEN: ~524 cataloged F-D9-01…06 tokens re-anchored
  (≈99% completeness; ~5 residue tokens for CONSOLIDATE). KAT suite identically green throughout.
  Durable discoveries D10–D14 folded into Discoveries log above. D. REFACTOR status: **done**.
