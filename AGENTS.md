# rDLP — agent guide

Orientation for anyone (human or agent) working in this workspace: the commands that are actually
verified to work, how this project treats tests that need an external oracle, and the vocabulary
discipline its durable text follows.

## Workspace

A Cargo workspace of nine members: `gnfs`, `rho`, `shor`, and six `shared/*` substrate crates
(`field`, `bigint`, `numth`, `numfield`, `padic`, `gf2m`). See `README.md` for the algorithmic map
and `docs/PEDAGOGY.md` for the section-by-section code tour.

## Commands that work today

See [`docs/development.md`](docs/development.md) for the full contributor guide: setup, toolchain
upgrades, and a precise enforced-vs-declared breakdown of every check command. Summary:

- Enforced (must stay green), each on its own cadence — `docs/development.md`'s "Verification
  cadence" is the authority: `cargo check --workspace` and `cargo test --workspace` once on the final
  inputs of each completed code change; `cargo coverage` (80% line-coverage floor, currently
  measuring around 92–93%) at the end of each planned group of related changes and at final delivery,
  earlier after substantial test removal; `cargo dep-policy` (dependency license/advisory/source
  policy) and `cargo audit` (RustSec advisory cross-check) when dependency inputs change, on the
  first active run of each day, and before a release. A result is reusable only if it was run on the
  same relevant inputs (the tree's content, not `HEAD` alone), is backed by a readable log, and, for
  the security checks, is fresh as that document defines. Report a check not yet due as deferred
  with its next checkpoint, never as passed.
- `cargo dep-policy` and `cargo audit` also depend on the date: a new advisory can turn them red
  against an unchanged tree. Treat that as a dependency-graph finding to triage, following
  `docs/development.md`'s "Dependency maintenance" section, not as a regression from your change.
- Declared and runnable, but not yet required to pass workspace-wide: `cargo fmt-check` and
  `cargo lint`. Both currently fail against the full codebase (a never-formatted tree and a large
  pre-existing lint backlog, respectively). Do not assume either is green, and do not "fix"
  formatting or lints as a drive-by side effect of an unrelated change — that rewrite belongs to its
  own dedicated change.
- `cargo build --release --workspace`, `cargo bench --bench <name>`, and the `rho-factor`/`rho-dlog`
  CLIs all work — see `README.md`'s build/bench/CLI sections for the exact invocations.

## Oracle policy

Tests that need an external number-theory tool — PARI/GP for algebraic cross-checks, CADO-NFS for
full sieving cross-checks — are `#[ignore]`-gated and excluded from `cargo test --workspace` by
design; their names say so (`kat_pari_*`, `kat_c_cado_*`, and similar). Run them manually with
`cargo test -- --ignored` once the corresponding oracle binary is on `PATH`. Never promote an oracle
test onto the default test path, and never substitute a benchmark run for a test run.

## Solver outcomes

Every solver success is exactly verified; failure classes stay distinct (budget exhaustion is never
a negative or a success); a fallback must be another published method through the same verifier,
never exhaustive search over the answer space. Full rules: `docs/development.md`, "Solver outcomes
and fallbacks".

## Register discipline

Source comments, identifiers, and this project's durable text (this file, `README.md`) carry only
mathematical and topic-native vocabulary — algorithm names, group/field/curve terms — never a
planning-session label, ticket ID, or phase number. `docs/NOTES.md` is where rolling, dated,
point-in-time framings and baselines belong instead; it is not itself a durable contract and code
never refers to it.
