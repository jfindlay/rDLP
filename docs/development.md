# Contributing to rDLP

A practical reference for building, checking, and testing this workspace. See
[`../README.md`](../README.md) for the algorithmic map, [`PEDAGOGY.md`](PEDAGOGY.md) for the
section-by-section code tour, and [`MATHEMATICS.md`](MATHEMATICS.md) for the maths-first textbook.

## Setup

The root [`../rust-toolchain.toml`](../rust-toolchain.toml) pins the workspace to an exact stable
release (currently 1.95.0) with the `rustfmt`, `clippy`, and `llvm-tools-preview` components.
`rustup` resolves this file automatically on any `cargo`/`rustc` invocation inside the workspace, so
no manual `rustup toolchain install` step is needed as long as that release is already present on
the machine; if it is not, the first `cargo` invocation in the workspace prompts `rustup` to fetch
it.

`rustfmt.toml` at the root declares `max_width = 100` explicitly, matching the workspace's general
line-length convention. `.cargo/config.toml` defines the command aliases used throughout this guide.

## Upgrading the pinned toolchain

Bumping the toolchain is a deliberate, verified step, not a one-line edit:

1. Edit the `channel` value in `rust-toolchain.toml` to the new exact release.
2. Re-run every command in the "Check commands" section below against the new toolchain and record
   the results in `docs/NOTES.md`, the same way the current baseline was measured. A different
   release can shift `clippy`/`rustfmt` diagnostics or dependency-audit results, so the existing
   readiness baselines only apply to the toolchain they were measured against — repeat the
   measurement rather than assuming the old numbers still hold.
3. Only update the components list if the new release changes a component name (for example, a
   future `rustup` may drop the `llvm-tools-preview` alias in favor of `llvm-tools`; the pin already
   uses the legacy alias name because that is what the currently pinned release accepts).

## Dependency maintenance

`cargo dep-policy` and `cargo audit` differ in kind from the other enforced gates: they depend on
the date as well as the tree. A new RustSec advisory can turn either red against an unchanged
`Cargo.lock`, so a failure there is a finding about the dependency graph, not necessarily a
regression from the change in hand. Triage it separately from the change being made.

No scheduler or CI runs these checks. Run `cargo audit` and `cargo dep-policy` periodically and
before cutting a release.

When an advisory fires, work through these in order:

1. **In-range fix.** `cargo update --dry-run`, then `cargo update -p <crate>` (or `--precise`).
   A requirement such as `"0.5"` means `>=0.5.0, <0.6.0`, so this only reaches versions the
   manifest already allows. It edits `Cargo.lock` only.
2. **Breaking upgrade.** If the patched release is outside the allowed range (for `0.x` crates
   the minor version is the breaking boundary), raise the requirement in every manifest that
   declares it, and in every crate that pins the same dependency transitively (dev-dependencies
   count, since `cargo audit` scans the whole lockfile). `cargo info <crate>` shows the current
   and latest version. `cargo update --breaking` exists but is unstable (nightly only), so edit
   the manifests directly. Do this on its own branch, then repeat every command in "Check
   commands" below.
3. **No fix available.** Only then add an entry to `deny.toml`'s `ignore` list, with the reason
   and the condition that would remove it. `deny.toml` carries no ignores today.

Dependencies that share a type across crate boundaries must move together. `crypto-bigint`'s
`rand_core` feature fixes which `rand_core` major version the workspace's `rand`/`rand_chacha`
must match, and `k256`/`p256` must match the `crypto-bigint` major version.

## Check commands

Run from the workspace root.

**Enforced — must stay green:**

| Command | What it checks |
|---|---|
| `cargo check --workspace` | Compiles every member without producing artifacts |
| `cargo test --workspace` | Runs every non-`#[ignore]`-gated test across all nine crates |
| `cargo dep-policy` (= `cargo deny check`) | Dependency licenses, advisories, and source policy |
| `cargo audit` | Cross-checks `Cargo.lock` against the RustSec advisory database |
| `cargo coverage` (= `cargo llvm-cov --workspace --fail-under-lines 80`) | Line-coverage floor |

**Declared and runnable, but not yet required to pass workspace-wide:**

- `cargo fmt-check` (= `cargo fmt --all -- --check`) reports diffs across most of the codebase —
  the workspace has never been passed through rustfmt as a whole. Do not treat this as a green
  gate, and do not fix formatting as a side effect of an unrelated change; a whole-workspace
  reformat is its own dedicated piece of work.
- `cargo lint` (= `cargo clippy --workspace --all-targets -- -D warnings`) fails outright: several
  hundred `clippy::all`-group hits plus a larger `pedantic`/doc-hygiene backlog. `cargo check`,
  `cargo test`, and plain `cargo clippy` (without `-D warnings`) stay usable in the meantime because
  `missing_docs` is only a `warn` at the workspace level; `cargo lint` promotes it (and the rest of
  `clippy::all`) to a hard error on top of that everyday loop. Do not attempt to close this backlog
  as a side effect of an unrelated change; it is a dedicated remediation effort.

`cargo format` (= `cargo fmt --all`) applies rustfmt in place and is unaffected by the
declared-but-unadopted distinction above — it is just not run automatically by any gate.

Also available and useful but not gates in themselves: `cargo build --release --workspace`,
`cargo bench --bench <name>`, and the `rho-factor`/`rho-dlog` CLIs — see `../README.md`'s build,
bench, and CLI sections for exact invocations.

## KAT-first coverage

The test suite is built around known-answer tests (KATs): field axioms, group laws, and
algorithm-level correctness checks against hand-verified or independently-computable results. This
is the primary correctness signal for the whole workspace — there is no proptest/quickcheck layer.
Two tiers exist by convention: inline `#[test]` functions for substrate-level verification
(arithmetic, group laws), and external `tests/*_kat.rs` files for algorithm-level KATs (one file per
algorithm family).

The 80% line-coverage floor enforced by `cargo coverage` is a secondary, complementary signal: it
catches dead or unexercised code paths. It is not a correctness ceiling, and clearing it (the
workspace currently measures around 92–93%) does not by itself establish that an algorithm is
correct — only its KATs do that. A change that raises coverage without adding or strengthening a
KAT has not, by itself, improved confidence in correctness.

## Opt-in oracles

Some tests need an external number-theory tool that is not part of this workspace's own dependency
graph: PARI/GP for algebraic cross-checks, CADO-NFS for full sieving cross-checks. These are
`#[ignore]`-gated by design, and their names say so (`kat_pari_*`, `kat_c_cado_*`, and similar
prefixes). They are excluded from `cargo test --workspace` on purpose.

To run them:

```sh
cargo test -- --ignored
```

once the corresponding oracle binary is on `PATH`. Never promote an oracle test onto the default
test path, and never substitute a benchmark run for a test run — a benchmark measures wall-clock
performance, not correctness.
