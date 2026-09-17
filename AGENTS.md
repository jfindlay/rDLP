# rDLP — agent guide

Orientation for anyone (human or agent) working in this workspace: the commands that are actually
verified to work, how this project treats tests that need an external oracle, and the vocabulary
discipline its durable text follows.

## Workspace

A Cargo workspace of nine members: `gnfs`, `rho`, `shor`, and six `shared/*` substrate crates
(`field`, `bigint`, `numth`, `numfield`, `padic`, `gf2m`). See `README.md` for the algorithmic map
and `docs/PEDAGOGY.md` for the section-by-section code tour.

## Commands that work today

- `cargo check --workspace` and `cargo test --workspace` are the only enforced gates; keep both
  green before finishing any change.
- `cargo build --release --workspace`, `cargo bench --bench <name>`, and the `rho-factor`/`rho-dlog`
  CLIs all work — see `README.md`'s build/bench/CLI sections for the exact invocations.
- `cargo fmt --all -- --check` and `cargo clippy --workspace --all-targets` both run (rustfmt and
  clippy ship with the active stable toolchain), but neither is a required gate yet, and neither is
  currently clean. Do not assume either is green, and do not "fix" formatting or lints as a
  drive-by side effect of an unrelated change — that rewrite belongs to its own dedicated change.
- `cargo deny check`, `cargo audit`, and `cargo llvm-cov` are not installed in this environment.
  Treat any claim about what they would report as unverified until the tool is actually present.

## Oracle policy

Tests that need an external number-theory tool — PARI/GP for algebraic cross-checks, CADO-NFS for
full sieving cross-checks — are `#[ignore]`-gated and excluded from `cargo test --workspace` by
design; their names say so (`kat_pari_*`, `kat_c_cado_*`, and similar). Run them manually with
`cargo test -- --ignored` once the corresponding oracle binary is on `PATH`. Never promote an oracle
test onto the default test path, and never substitute a benchmark run for a test run.

## Register discipline

Source comments, identifiers, and this project's durable text (this file, `README.md`) carry only
mathematical and topic-native vocabulary — algorithm names, group/field/curve terms — never a
planning-session label, ticket ID, or phase number. `docs/NOTES.md` is where rolling, dated,
point-in-time framings and baselines belong instead; it is not itself a durable contract and code
never refers to it.
