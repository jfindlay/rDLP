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
| 8 | Project the generic rho-dlog curve onto its largest prime-order subgroup | B | T1 | C-DlogOutcome (prior arc) | `rho/src/curve/{generic,mod}.rs`, `rho/src/ecdlp/cli.rs`, `rho/tests/ecdlp_kat.rs`, `README.md`, `docs/PEDAGOGY.md` §6 |
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

### 8 — Project the generic rho-dlog curve onto its largest prime-order subgroup

User-approved revision (2026-11-10): retain the curve and accept its small subgroup; do not search
for a replacement curve or remove `generic`. Deliverable per C-GenericSubgroup: with the original
P₀ = (3, 821_487_384_573_098_969), store G′ = c·P₀ =
(3_569_091_914_278_015_591, 7_073_522_384_094_063_206), c = 8_218_965_179_160, and prime order
n′ = 561_103 in `generic_curve()` (`generic::N` names the order). The CLI uses the stored order;
its separate `GENERIC_N` and "Computed by Sage" provenance go. No solver algorithm changes.

Rustdoc records the original point and multiplier, with PARI/GP discovery provenance: #E =
9_223_372_037_844_426_960 and ord P₀ = 4_611_686_018_922_213_480. These full orders and the claim
that n′ is the largest prime factor are externally computed, not proved by the subgroup KATs.
The in-tree proof establishes G′'s exact prime order and its relation to P₀ independently of those
full-order claims. Keep the original reference multiples and arithmetic-consistency tests on P₀
built explicitly, not on the changed `generator()`.

KATs: G′ on-curve, G′ ≠ ∞, n′·G′ = ∞, `is_prime(n′)`, c·P₀ = G′; Brent's input check accepts
(G′, n′), and `solve_brent` recovers a fixed nonzero k from (G′, k·G′), with a recorded seed.
Exercise the existing rejection paths with P₀: as base point it gives `WrongGeneratorOrder`;
as target with base G′ it gives `TargetOutsideSubgroup`. Do not silently project caller-supplied
CLI points. If c·P₀ ≠ G′ or n′·G′ ≠ ∞ in-tree: destructive-HALT; no constant adjustment or search.

README, PEDAGOGY §6 and `curve/mod.rs` describe a 63-bit field with a 20-bit subgroup of order
561103 (approximately 2¹⁹), suitable for runnable demonstrations, not a large-group performance
baseline. Explain that supplied base/target points must belong to that subgroup and provide a
working CLI example. Keep the warning about default runs that can exhaust their budget
(secp-toy Brent: per-attempt cap 2²⁹ − 1 vs ≈ 2.7·10⁹ expected steps, arithmetic only).

### 9 — Retarget the ECDLP benches to groups that finish

Deliverable per C-BenchFixture: every bench in `rho/benches/ecdlp.rs` asserts its solve succeeds
once before timing (as `bench_ssa` does) on a group sized so one iteration finishes in seconds in
release; reuse an existing prime-order fixture if one fits, else add one (prime order proved
in-tree; generator and order found by a bounded, recorded search — an external tool may find,
never certify). Budgets per C-DpBudget. BENCHMARKS: replace the ecdlp-bench setup description;
mark old numbers as not re-measured (re-measurement belongs to the tooling scope).

Choose fixtures independently of row 8's CLI repair: `generic` may be reused where it fits, but
its small group does not by itself establish representative DP/parallel performance and does not
support GLV. State the actual field/subgroup sizes and relevant tuning parameters; do not infer
algorithmic speedups from unlike groups or smoke-test completion. Replace BENCHMARKS's `generic`
"unknown†" row and footnote with the corrected subgroup description and provenance.

KATs: the fixture's prime-order test; `cargo bench -p rho --bench ecdlp -- --test` completes
(smoke, recorded; not a test substitute).

### 10 ◆ — Correct the follow-up exposition

Deliverable: (a) MOV R = G: `DegeneratePairing`'s rustdoc, the `vertical_at` comment
(`miller.rs:235–249`) and MATHEMATICS §10.5 / PEDAGOGY §10 say the evaluation point lies on the
Miller divisor's support (the value is undefined; the tangent-line numerator vanishes), reported
as `DegeneratePairing`. (b) Sweep superseded claims across all durable surfaces these solvers
touch, including unchanged passages: the SSA halving and "`LiftFailed` means bad input"; GHS
"valid only on the fixture"; MOV `debug_assert`-only domain; ECM "B1 ≳ 177 wraps"; `generic`
composite or wrong order (9_223_372_037_218_517_353 is neither #E nor a multiple of ord P₀),
including BENCHMARKS's `generic` row/footnote; rho factoring "probable prime in [ψ₁₂, 2^80) is
`Unsupported`"; `pohlig.rs`'s ECM mention. Keep field size distinct from subgroup size and external
full-order provenance distinct from in-tree subgroup proof. (c) Anneal the arc's durable files
against the denylist.

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

Compiler + test-enforced. **Frozen at 1 by the inflection juncture (2026-11-10).** Consumed by 2,
10. Row 2 adds no variant, renames nothing, and changes no payload; it only narrows the docs.

```rust
// rho/src/ssa/mod.rs
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SsaError {
    InvalidInput(SsaInputDefect),
    NotAnomalous,                    // proven negative: exact count on a valid curve, #E ≠ p
    Unsupported(&'static str),       // unchanged: p, a or b above one limb, or p < 5
    Inconsistent { k: u64 },         // the formal-log quotient mod p fails k·G == Q
    InvariantFailure(&'static str),  // a breach with no candidate k
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SsaInputDefect { FieldNotPrime, SingularCurve, GeneratorAtInfinity,
                          GeneratorNotOnCurve, TargetNotOnCurve, OrderNotFieldPrime }

// rho/src/ssa/lift.rs (lift helpers get their own error type, not SsaError)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LiftError { PointAtInfinity, Hensel(HenselError), Padic(ZpError) }
pub fn lift_point<F: Fp<4>>(..) -> Result<(Zp, Zp), LiftError>;   // row 2 may change the args
pub fn check_lift_on_curve(..) -> Result<bool, ZpError>;
```

**Order of checks in `ssa_solve`** (the first failure is returned): (1) `p`, `a`, `b` fit one limb
and p ≥ 5, else `Unsupported`; (2) p prime (`shared_numth::is_prime`, deterministic on u64), else
`FieldNotPrime`; (3) 4a³ + 27b² ≢ 0 mod p, else `SingularCurve`; (4) G ≠ ∞; (5) G on the curve;
(6) Q on the curve; (7) Q = ∞ → `Ok(0)`; (8) exact count, `NotAnomalous`; (9) n = p, else
`OrderNotFieldPrime`. Then the reduction runs. Checks 1–3 are one shared curve check, so
`verify_anomalous` also returns `InvalidInput(FieldNotPrime | SingularCurve)` and never `Ok(_)`
for a composite modulus or a singular cubic. Why this order: the count is a proven negative only on
a nonsingular curve over a prime field (a split-node cubic counts exactly p). Q = ∞ comes before the
count and n is checked after it because the existing `kat_ssa_solve_infinity_is_zero` (non-anomalous
curve, n = 12 → `Ok(0)`) and `kat_ssa_solve_not_anomalous` (n = 12 → `NotAnomalous`) must stay
unchanged. Before the count, n = #E is not an error. With this order every existing solve KAT and
the `bench_ssa` `Ok(3)` precheck keep passing unedited.

**Dropped from the draft:** `WrongGeneratorOrder` and `TargetOutsideSubgroup`. Once checks 1–9
pass, E(F_p) has prime order p, so every finite on-curve point generates it: n·G = ∞, n·Q = ∞ and
Q ∈ ⟨G⟩ all follow, and the variants would be dead. Row 1 records that implication in rustdoc and
does not run those scalar multiplications.

**Not reachable after validation, so removed from `SsaError`:** `LiftFailed` and `Padic`, along with
both `From<HenselError>` and `From<ZpError>` for `SsaError`. After checks 1–9: G ≠ ∞ and Q ≠ ∞;
y₀² ≡ x₀³ + ax₀ + b; y₀ ≠ 0 because an odd prime order has no 2-torsion; p ≥ 5; k ≥ 1; and one prime
throughout. So `hensel_lift`, `Zp::new`/`add`, and `padic_log`'s `Zp` arm cannot fail. Any such
error inside `ssa_solve` is mapped explicitly with `map_err` to `InvariantFailure`, never with `?`
through a `From`. The direct `lift_point` KATs keep their coverage through `LiftError`:
`kat_lift_2torsion_errors_non_simple_root` matches `LiftError::Hensel(NonSimpleRoot{..})`.
`PointAtInfinity` replaces the made-up `NotARoot` that was built with `r0 = 0` for the identity.

**Invariant class.** `InvariantFailure(msg)`, with one distinct `&'static str` per breach site that
names the violated property, covers these sites:
- in `elliptic_log_proj`: v_p(t) < 1, k_t < 1, and `ConvergenceViolation`;
- a zero residue at working precision, which means the valuation is at least the precision. This
  is treated as a breach and is never read as valuation 0 (`padic_valuation_bigint` returns 0 for
  zero);
- in `fp_divide_logs`: a g-log that is not a unit after the common p^v is divided out, or a zero
  g-log;
- any lift, `Zp` or `padic_log` error on validated input.

`Inconsistent { k }` stays the class for a candidate that fails verification. The field is renamed
from `k_raw` to `k`. In row 1 it is the raw quotient after the halved candidate also failed; from
row 2 on it is the only candidate.

Interim honesty: until row 2 makes the lift satisfy valuation 1, the rustdoc of both invariant
variants keeps today's caveat: on a degenerate lift of the stored canonical (a = 0) curve they can
occur on valid input. Row 2 removes the caveat.

**Seams for row 1's KAT:** a unit test calls `elliptic_log_proj` with a point outside the kernel
of reduction (X and Y both units at precision 8) and gets `InvariantFailure`. A unit test calls
`fp_divide_logs` with v(ψ_G) > v(ψ_Q) and also gets `InvariantFailure`. `Inconsistent` keeps the
`resolve_candidate_*` seam in row 1. Row 2, which deletes `resolve_candidate`, keeps a private
verification seam that is unit-tested for `Inconsistent { k }`.

**Second-lift failure (row 2): invariant class, with no exhaustion variant.** After `p·G̃ ∈ E₂` on
shift s₁, row 2 retries with one s₂ ≢ s₁ (mod p), s₂ ≢ 0. If p·G̃ is still in E₂, the result is
`InvariantFailure("p·G̃ in E₂ on two shifted lifts")`. The reasoning has three steps:
1. Inferred from Serre–Tate theory, recalled rather than re-read; row 2 cites a source or labels
   it. p·P̃ ∈ E₂ for some P ≠ O, which by cyclicity means for all P, exactly when the reduction
   sequence of Ẽ(Z/p²) splits. That happens exactly when Ẽ mod p² is the canonical lift. Only one
   isomorphism class mod p² is canonical.
2. Confirmed algebra. Short-Weierstrass lifts (a + pα, b + pβ) are isomorphic over Z/p² exactly
   when they differ by a multiple of (4a, 6b). The shift direction (1, 0) is parallel to that only
   if b ≡ 0. An anomalous curve has b ≢ 0, because b ≡ 0 puts (0, 0) in E[2] and makes #E even. So
   distinct s mod p give distinct classes, and at most one s mod p lands in E₂.
3. Corroborated by the stored measurements: every a = 0 lift (all b̃) is in E₂, which is one class;
   a = 7, 14, 21 give valuation 1.

Test obligation for row 2: on each SSA fixture, a KAT sweeps every s mod p and asserts that at most
one class puts p·G̃ in E₂. That turns this classification into a test-enforced property. If the
sweep finds two classes, that is the existing destructive-HALT discovery ("shifted lifts still land
in E₂").

**Display / source obligations.**
- Display texts:
  - `InvalidInput(d)` → "invalid SSA input: {d}". `SsaInputDefect` implements `Display`, one phrase
    per variant: field modulus not prime, curve singular (4a³ + 27b² ≡ 0 mod p), G at infinity, G
    not on curve, Q not on curve, n ≠ p (the order of every anomalous curve).
  - `NotAnomalous` and `Unsupported` keep their current text.
  - `Inconsistent { k }` → the quotient {k} mod p fails k·G = Q, an invariant failure. The text must
    not mention a halved candidate, so it is true in both rows.
  - `InvariantFailure(why)` → "invariant failure: {why}".
- `impl std::error::Error for SsaError` has `source() == None` for every variant, since no wrapped
  error survives.
- `LiftError` implements `Display`, plus `Error` with `source` returning the wrapped
  `HenselError`/`ZpError`.
- Tests match variants, not Display text.

### C-SsaLift

Test + prose-enforced. To be frozen at 2. Consumed by 10. Each attempt computes one lifted curve
(a ↦ a + p·s) and uses it for both point lifting and multiplication by p. If p·G̃ ∈ E₂, there is at
most one retry, with s₂ ≢ s₁, and a second E₂ result is `InvariantFailure` (per C-SsaOutcome). k is
the raw quotient on a lift with p·G̃ ∈ E₁ \ E₂. No halved candidate exists.

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

Test-enforced. Revised with user approval on 2026-11-10; to be frozen after implementation at 8.
Consumed by 9, 10. The stored `generic` generator G′ = 8_218_965_179_160·P₀ has exact prime order
n′ = 561_103, proved in-tree by on-curve/nonidentity checks, annihilation and primality; its
relation to P₀ is also checked in-tree. The CLI uses the stored order and never silently projects
supplied points. Original-point arithmetic KATs remain on P₀. The fixture is a runnable small
subgroup over a 63-bit field, not a promised large-group or representative parallel benchmark.

### C-BenchFixture

Test + prose-enforced. To be frozen at 9. Consumed by 10. Every ECDLP bench solves its instance
before timing; any fixture it adds has in-tree prime-order proof and recorded provenance.

## Progress ledger

Bindings: PLAN = `docs/PLAN.md`; VERIFY_TYPES = `cargo check --workspace`; VERIFY_TEST =
`cargo test --workspace`; boundary check = `cargo coverage`; freshness checks = `cargo
dep-policy`, `cargo audit` (separate commands, never combined). Cadence per "Verify gate".

Gate evidence (one line per check: check, outcome passed/reused/deferred, command and working
directory, exit status, time, log reference, input identity, environment; deferred → next
checkpoint): 

- `cargo dep-policy` — passed — `/home/justin/Source/rDLP` — exit 0 — 2026-11-10T08:33:39-05:00 —
  log `/tmp/opencode/rdlp/dep-policy.log` — inputs: base `a7f72ba`, clean tree — advisory db
  `~/.cargo/advisory-dbs/advisory-db-3157b0e258782691`, refreshed this run by `cargo deny`.
- `cargo audit` — passed — `/home/justin/Source/rDLP` — exit 0 — 2026-11-10T08:33:39-05:00 — log
  `/tmp/opencode/rdlp/audit.log` — inputs: base `a7f72ba`, clean tree — advisory db
  `~/.cargo/advisory-db` @ `8a4faeaaa8535304bc640bc34c8764d960b15850`, fetched this run (1439
  advisories).

- Row 1 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T08:45:14-05:00 — log `/tmp/opencode/rdlp/ssa-typed-check.log` — input identity
  `ef6f15f8…7169` = sha256 of (`git diff HEAD -- . ':!docs/PLAN.md'` + sha256sum of untracked
  files) over base `a7f72ba`, recomputed by the driver and equal — rustc/cargo 1.95.0, dev profile,
  default features.
- Row 1 `cargo test --workspace` — reused (worker record, validated, same identity) — exit 0 —
  2026-11-10T08:48:53-05:00 — log `/tmp/opencode/rdlp/ssa-typed-test.log` (no failures) — same
  environment.
- Row 2 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T09:03:05-05:00 — log `/tmp/opencode/rdlp/check_final.log` — input identity
  `bfb0eab9…5159` (same recipe, base `63381d0`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features.
- Row 2 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T09:07:30-05:00 — log
  `/tmp/opencode/rdlp/test_final.log` (86 `test result: ok`, no failures) — same environment.
  Smoke (not a gate): `cargo bench -p rho --bench attacks -- --test ssa` exit 0, log
  `/tmp/opencode/rdlp/bench_smoke.log`.
- Row 2 `cargo coverage` — deferred → ◆ boundary (row 10). Judged not "substantial test removal":
  3 unit tests deleted, 13 added.
- Row 3 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T09:23:05-05:00 — log `/tmp/opencode/rdlp/ghs-check-final3.log` — input
  identity `28d7ac9d…6d39` (same recipe, base `3c1c9dd`; recomputed by the driver and equal) —
  rustc/cargo 1.95.0, dev profile, default features.
- Row 3 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T09:27:38-05:00 — log
  `/tmp/opencode/rdlp/ghs-test-final3.log` (86 `test result: ok`, no failures) — same environment.
- Row 3 `cargo coverage` — deferred → ◆ boundary (row 10); no tests deleted (6 added).
- Row 4 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T09:51:35-05:00 — log `/tmp/opencode/rdlp/mov-check-workspace.log` — input
  identity `446638c0…4677` (same recipe, base `b0065a4`; recomputed by the driver and equal) —
  rustc/cargo 1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 4 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T09:55:40-05:00 — log
  `/tmp/opencode/rdlp/mov-test-workspace.log` (86 `test result: ok`, no failures) — same
  environment. Smoke (not a gate): `cargo test -p rho --release --test mov_kat
  mov_reduce_rejects_field_above_one_limb` exit 0, 1 passed (worker-reported, no log).
- Row 4 `cargo coverage` — deferred → ◆ boundary (row 10); no tests deleted (2 added).
- Row 5 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T10:00:00-05:00 — log `/tmp/opencode/rdlp/ecm-check.log` — input identity
  `5543ad0a…6559` (same recipe, base `2748ef5`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 5 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T10:03:58-05:00 — log
  `/tmp/opencode/rdlp/ecm-test.log` (86 `test result: ok`, no failures) — same environment.
- Row 5 `cargo coverage` — deferred → ◆ boundary (row 10); no tests deleted (4 added).
- Row 6 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T10:07:27-05:00 — log `/tmp/opencode/rdlp/inc-check.log` — input identity
  `21d74488…875d` (same recipe, base `347dae2`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 6 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T10:11:24-05:00 — log
  `/tmp/opencode/rdlp/inc-test.log` (86 `test result: ok`, no failures) — same environment.
- Row 6 `cargo coverage` — deferred → ◆ boundary (row 10); 1 test deleted, 17 added.
- Row 7 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T10:14:51-05:00 — log `/tmp/opencode/rdlp/rhob-check.log` — input identity
  `b9df6552…7de1` (same recipe, base `afd07af`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 7 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T10:19:14-05:00 — log
  `/tmp/opencode/rdlp/rhob-test.log` (86 `test result: ok`, no failures) — same environment.
- Row 7 `cargo coverage` — deferred → ◆ boundary (row 10). Judged not "substantial test removal":
  1 test and 1 branch + 1 constant deleted, 2 tests added, 1 extended.
- Row 8 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T12:31:53-05:00 — log `/tmp/opencode/rdlp/s8-check.log` — input identity
  `5aa7e47a…af36` (same recipe, base `3db95b1`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 8 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T12:35:10-05:00 — log
  `/tmp/opencode/rdlp/s8-test.log` (86 `test result: ok`, no failures) — same environment. Manual
  check (not a gate): release `rho-dlog --curve generic --q 1000932015531192845,9048581680910854061`
  → k = 424242 verified (worker-reported, no log).
- Row 8 `cargo coverage` — deferred → ◆ boundary (row 10); no tests deleted (3 added, helper
  signature changed).
- Row 9 `cargo check --workspace` — reused (worker record, validated) — `/home/justin/Source/rDLP` —
  exit 0 — 2026-11-10T13:25:54-05:00 — log `/tmp/opencode/rdlp/check-final.log` — input identity
  `aee7cdff…a09a` (same recipe, base `3270821`; recomputed by the driver and equal) — rustc/cargo
  1.95.0, dev profile, default features. No manifest or lockfile changed.
- Row 9 `cargo test --workspace` — reused (same identity) — exit 0 — 2026-11-10T13:29:43-05:00 — log
  `/tmp/opencode/rdlp/test-final.log` (86 `test result: ok`, no failures) — same environment. Smoke
  (not a gate, not a test substitute): `cargo bench -p rho --bench ecdlp -- --test` exit 0, release
  bench profile, 2026-11-10T13:30:58-05:00, log `/tmp/opencode/rdlp/smoke-final.log` (8 `Success`
  lines).
- Row 9 `cargo coverage` — deferred → ◆ boundary (row 10); one test added, no test deleted (one bench
  function, `bench_glv_vs_batch`, removed — benches are outside the coverage measure).
- Row 1 `cargo coverage` — deferred → ◆ boundary (row 10) or earlier after substantial test removal
  (row 2 or 7).

Pending checks: `cargo coverage` — due at the ◆ boundary (row 10) and final delivery, earlier after
substantial test removal (rows 2, 7); `cargo dep-policy`/`cargo audit` — due when dependency inputs
change and on first active run of each new day.

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Type SSA input and invariant failures | done | 829caaa | C-SsaOutcome (frozen by juncture; per-unit gate passed, coverage deferred to ◆) |
| 2 | Lift SSA curves non-canonically and drop the halved candidate | done | 0aa2b01 | C-SsaLift (test-enforced sweep; allowed extra `formal_log.rs`; coverage deferred to ◆) |
| 3 | Identify the GHS conorm subfield with `restrict` | done | ae1f1cc | C-GhsSubfieldIdentification (coverage deferred to ◆; `m > 63` has no variant, unreachable by type) |
| 4 | Guard MOV's single-limb field domain in release builds | done | 32452cd | C-MovDomain (coverage deferred to ◆) |
| 5 | Apply ECM stage-1 prime powers one ladder call at a time | done | a4a040f | C-EcmStage1 (coverage deferred to ◆) |
| 6 | Report non-dividing factor candidates as `Inconsistent` | done | 7ce7340 | C-FactorInconsistent (allowed extra `rho/src/factor/cli.rs`; coverage deferred to ◆) |
| 7 | Certify rho-factoring primes below 2^80 with the 13-base bound | done | 53b072e | C-RhoPrimeBound (allowed extras: `factor/{mod,cli}.rs`, both PEDAGOGY files; coverage deferred to ◆) |
| 8 | Project the generic rho-dlog curve onto its largest prime-order subgroup | done | e1764e2 | C-GenericSubgroup (user-approved revision; c·P₀ = G′ and 561103·G′ = ∞ held in-tree, no HALT; coverage deferred to ◆) |
| 9 | Retarget the ECDLP benches to groups that finish | done | f656f03 | C-BenchFixture (new 25-bit fixture `prime_order_25`; GLV bench dropped; coverage deferred to ◆) |
| 10 | Correct the follow-up exposition ◆ | pending | | |
| 11 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

### 1 — 2026-11-10
Discovery/flex: Inflection juncture self-cleared (`design-confident`) but narrowed the row's draft
enum: `WrongGeneratorOrder`/`TargetOutsideSubgroup` dropped (dead once checks pass), `FieldNotPrime`/
`SingularCurve` added, `LiftFailed`/`Padic` removed from `SsaError` (new `LiftError` for
`lift_point`), `k_raw` renamed `k`; second-lift double-E₂ failure is `InvariantFailure`, no
exhaustion variant.
Affected: C-SsaOutcome (frozen at 1), C-SsaLift (retry text updated by the juncture).
Deferred: yes — row 2 must (a) cite or label the Serre–Tate "only the canonical lift lands in E₂"
claim (recalled, not re-read) and (b) add a KAT sweeping every s mod p per fixture asserting at
most one lands in E₂; two classes would be the destructive-HALT discovery. Row 1 also made a zero
q-log in `fp_divide_logs` an `InvariantFailure` (not `Inconsistent{k:0}`) under the zero-residue
rule; row 2 should keep that when it deletes the halving.
Texture (row 1): Rustdoc on both invariant variants carries a caveat that they can occur on valid input
via the a = 0 lift; row 2 must remove it. Committer bodies came back unwrapped twice (baseline and
row 1); future dispatches should say to pass a hard-wrapped body via `-F`.

### 2 — 2026-11-10
Discovery/flex: Premise confirmed, no contract flexed: shifts s = 1..6 on `anomalous_toy` all give
valuation 1 and exactly one shift class (the stored a = 0 lift) is degenerate per fixture; a second
fixture (y² = x³ + x + 5 over F₁₁, G = (0, 4)) found by a recorded search agrees, and an
exploratory scan of 481 anomalous curves (5 ≤ p ≤ 59) also found one degenerate class each.
Working precision dropped from 8 to `WORKING_PRECISION` = 4 (k = 3 sufficed on the whole scan).
Affected: C-SsaLift (frozen at 2); `reduce::kernel_filtration_level` added as a public diagnostic.
Deferred: yes — the ◆ stale-claim sweep should check MATHEMATICS §10.2 / PEDAGOGY §11 and any
rustdoc for surviving "halved candidate", "p·G̃ = O in E(Z/p²)", "precision 8" or "LiftFailed means
bad input" wording. Smart 1999 was read: it states the failure on the canonical lift and the remedy
(another lift) but does not derive the mechanism; the docs label the Serre–Tate "only the canonical
lift lands in E₂" explanation as unverified, backed by the sweep. `lift.rs` comments still carry
"Principle-4" wording (plan vocabulary, left out of scope by the worker) — anneal candidate.

### 3 — 2026-11-10
Discovery/flex: The independent hand derivation (confirmed with `gp`) of the two transfer constants
agrees with `restrict` and disagrees with the old raw-bits output: the old constants were wrong
(raw reading swaps β and β² on GF(4) ⊂ GF(2⁶); D_P = [X+1, 2], D_{G+P} = [X+1, 3]). The transfer
itself is not suspect, so the planned destructive-HALT condition (disagreement with both) did not
fire. Homomorphism/log-preservation KATs pass under both readings because the fixture curve has F₂
coefficients and is conjugation-symmetric; a new KAT with a = β discriminates (fails on raw bits).
Affected: C-GhsSubfieldIdentification (frozen at 3). `GhsError` gained `SubfieldPolyMismatch` (and
`SubfieldDivisibility`); a pub `check_subfield_identification` was added in `ghs/mod.rs`.
Deferred: yes — the ◆ sweep should check that no durable text still calls the transfer "valid only
on the fixture" (worker fixed two PEDAGOGY lines outside §12) and that GHS statements elsewhere
(README, MATHEMATICS summaries) match. No `m > 63` variant exists: `Uint<1>` caps m at 63 by type,
replaced by a const assertion — the contract's "typed error for m > 63" is vacuous, not unmet.
Texture: one over-width line added to the PEDAGOGY test table (pre-existing over-width row);
ledger it for the anneal sweep.

### 4 — 2026-11-10
Discovery/flex: The inferred truncation at `tate.rs:147` was a live wrong-pairing path: with p ≈ 2^36,
ℓ = 3, k = 2 (y² = x³+1, fixture built with `gp`) the truncated exponent put `reduced_tate` outside
μ_ℓ. It is now exact (`BigUint`, panics only above 256 bits, which also removes the p^k overflow
for k ≥ 3). No existing fixture is affected (`pairing_toy` p = 47). The truncation is reachable
only through `reduced_tate` directly: `mov_reduce` end-to-end on a p > 2^32 field never finishes
because the out-of-scope k = 2 log engine's `find_irreducible_degree2` loops over all r in 0..p.
Affected: C-MovDomain (frozen at 4); no downstream contract.
Deferred: yes — the ◆ exposition row must not claim `mov_reduce` handles k = 2 fields near 2^64:
it accepts them but the log engine's irreducible search is O(p) (additive-reshard class, out of
scope). Also the new exponent test lives in `mov.rs` unit tests (`FpExt::coeffs` private).
Texture: one over-width row added to the PEDAGOGY outcome table (matches neighbours); the worker's
15-minute end-to-end attempt was removed; no stray processes remained.

### 5 — 2026-11-10
Discovery/flex: Premise confirmed: Z ≡ 0 mod a factor persists through later ladder calls; stage 2
does not depend on stage-1 point semantics (a proper factor returns earlier; Z ≡ 0 mod N is
`Inconclusive` before stage 2). The old product first exceeds 2^256 at B1 = 179 (not "≳ 177": the
latter was 256·ln 2); docs now say 179. Per-prime gcd chosen (a factor is caught before a second
prime of N also reaches the identity). Fixture n = (2^61−1)·1009, σ = 15, B1 = 250 (recorded search,
< 1 s, first hit).
Affected: C-EcmStage1 (frozen at 5); `stage1_scalar` removed; no downstream contract change.
Deferred: yes — the ◆ sweep should look for surviving "B1 ≳ 177 wraps" text outside the three edited
files: `docs/ROADMAP.md` ~:186 still describes the old scalar (ROADMAP is outside the durable-file
anneal exclusions but is a rolling doc — note, do not edit unless the row-10 scope covers it) and
`pohlig.rs`'s ECM mention. Row 6 edits the same `ecm.rs` (adds `EcmError::Inconsistent`).
Texture: worker used `python3` for a read-only calc and to splice one test into `ecm.rs`
(tool-discipline slip, no effect on correctness); the ladder fn's own doc comment was left alone —
re-read it at ◆.

### 6 — 2026-11-10
Discovery/flex: No contract flex. Both solvers now route candidates through a private
`judge_candidate` + search loop (`rho.rs`: `search_constants`; `ecm.rs`: `search_curves`) and
`verified_divisor` is gone from both. Neither error enum has `Display`; the property is in variant
doc comments and the CLI message. A compiler-forced `Inconsistent` arm was added to
`rho/src/factor/cli.rs`.
Affected: C-FactorInconsistent (frozen at 6); no downstream contract.
Deferred: yes — two human docs omit the new variants and are outside row 6's files:
`shared/numth/docs/PEDAGOGY.md` (~:317) and `docs/PEDAGOGY.md` (~:767) list the factor/ECM outcomes;
row 10's sweep must add `Inconsistent`. Row 7 edits `rho.rs` and `factor_kat.rs` again.
Texture: ECM and rho unit tests live in-module (new `mod tests` in `rho.rs`).

### 7 — 2026-11-10
Discovery/flex: The recalled constants were read and confirmed: Sorenson–Webster (arXiv:1509.00864,
Math. Comp. 2017) Theorem 1.1, ψ₁₂ = 318665857834031151167461, ψ₁₃ = 3317044064679887385961981;
2^80 < ψ₁₃; `is_prime` uses 12 bases below ψ₁₂ and the first 13 primes on [ψ₁₂, ψ₁₃). The
additive-reshard fallback for an unreadable paper did not fire. `MILLER_RABIN_PROOF_BOUND` already
existed; `PRIME_PROOF_BOUND` and the `Unsupported`-prime branch are gone; compile-time
`MAX_SUPPORTED ≤ MILLER_RABIN_PROOF_BOUND` assertion added.
Affected: C-RhoPrimeBound (frozen at 7); no downstream contract.
Deferred: yes — `docs/ROADMAP.md` ~:206–209 still says rho certifies only below ψ₁₂ and that the
constants are "recalled, not re-read" (rolling doc, outside any row's files; ◆ sweep decides whether
to touch it or leave it for the roadmap owner). Worker fixed `docs/PEDAGOGY.md` :777 and numth
PEDAGOGY :465 outside the listed files — row 10's sweep should re-read them.
Texture: none.

### 8 — 2026-11-10 (historical HALT; resolved by the approval below)
Discovery/flex: The row-8 worker (read-only, `gp`) found the plan's n′ = 77_507_328_043_853_087 is
prime but is NOT a factor of #E for `generic` (y² = x³ − 3x + 1 over p = 2⁶³ − 25, P₀ = (3,
821_487_384_573_098_969)); n′·(119·P₀) ≠ ∞. The juncture fork verified it independently: #E =
9_223_372_037_844_426_960 = 2⁴·3²·5·1493·1697·9011·561103; ord P₀ = 2³·3²·5·1493·1697·9011·561103;
stored `GENERIC_N` = 7·17·n′ is neither #E nor a multiple of ord P₀ (the "Computed by Sage" note is
false; it merely lies inside the Hasse interval). Today `--curve generic` reports `OrderNotPrime`
(`outcome.rs:149–154` tests primality before n·P = ∞), hiding the wrong order. No `Ok` is wrong; no
test, bench or solver run depends on `generic`; docs citing the composite order: README :231,
PEDAGOGY :877, BENCHMARKS :73–80, ROADMAP :194 (rolling), `curve/mod.rs:17`. The largest prime
subgroup is ℓ = 561_103 (≈ 2¹⁹), cofactor c = 8_218_965_179_160, G′ = c·P₀ =
(3569091914278015591, 7073522384094063206), ℓ·G′ = ∞ (gp).
Affected: C-GenericSubgroup (not yet frozen; text as written is unmeetable); C-BenchFixture and
row 9 unchanged. Row 8 made no edits; tree clean.
Verdict: juncture `additive-reshard` — redefine row 8 as "Project the generic rho-dlog curve onto
its largest prime-order subgroup": `generic_curve()` stores G′ above and n = 561_103 (`generic::N`);
`cli.rs` drops `GENERIC_N`/Sage note; rustdoc records P₀, c, #E, ord P₀ ("found with PARI/GP;
certified in-tree only by the KATs"); README and PEDAGOGY §6 describe a 63-bit field with a
prime-order subgroup of order 561103 and keep the secp-toy Brent exhaustion statement. KATs: G′ ≠ ∞,
561103·G′ = ∞, `is_prime(561103)`, c·P₀ = G′, Brent's input check accepts (G′, 561103), and
`solve_brent` recovers a fixed k from (G′, k·G′). Existing P₀ reference multiples/linearity tests
stay on P₀ built explicitly (not `generator()`). Files: the original row-8 list plus
`rho/src/curve/mod.rs:17` doc line. If in-tree c·P₀ ≠ G′ or 561103·G′ ≠ ∞: destructive-HALT (curve
arithmetic suspect; no search, no constant adjustment). C-GenericSubgroup new text: "stored
generator G′ = 8_218_965_179_160·P₀ has proven prime order n′ = 561_103 (all in-tree); the CLI uses
the stored order". Row 9 entry unchanged (it may reuse the 2¹⁹ subgroup, `tiny_a`, `tiny_b`; its
BENCHMARKS rewrite replaces the `generic` row "unknown†" and footnote). Row 10 (b) becomes
"`generic` composite or wrong order — 9_223_372_037_218_517_353 is neither #E nor a multiple of
ord P₀, including the BENCHMARKS `generic` row/footnote".
Not applied: the run was invoked without `may-reshard`, and the verdict changes what the user chose
(the plan's ≈2⁵⁶ group). User options: (a) accept the 19-bit `generic` (the spec above — resume with
`/plan-run docs/PLAN.md may-reshard`, or edit PLAN's row 8/contract/row 10 yourself and resume
without the flag); (b) replace `generic` with a curve over 2⁶³ − 25, a = −3 and a new b whose order
has a ~2⁵⁰–2⁵⁶ prime factor, found by a bounded recorded b-search and proved in-tree (follow-up
item; regenerates the reference KATs); (c) drop `generic`.
Deferred: yes — see the verdict spec; ◆ sweep must also catch the `generic` composite/wrong-order
claim in the five places listed.
Texture: the worker correctly refused to adjust constants; the fork corrected the worker's claim
that today's CLI fails with `WrongGeneratorOrder` (it is `OrderNotPrime`).

### 8 — 2026-11-10 (user-approved revision; implementation pending)

Decision: accept the 561103-order subgroup, retaining the curve and the original-point reference
tests. The row-8 deliverable, C-GenericSubgroup, row-9 benchmark caveats and row-10 stale-claim
sweep above now supersede the unapplied verdict/options in the historical halt entry. Only
`curve/mod.rs` is added to row 8's expected files; no solver, dependency or algorithm scope grows.
Confirmed by source inspection: the CLI holds the incorrect separate order constant; existing
`generic` KATs exercise arithmetic, and the ECDLP benches use `secp_k1_toy`, not `generic`.
Judgment: a runnable small subgroup serves this repair better than a replacement-curve search;
it sacrifices large-group experiments but preserves independent arithmetic references. A future
large-group fixture needs a separate requirement, not an assumed obligation of this repair.
Safeguards: prove the subgroup in-tree, preserve rejection of supplied points outside it, and
separate field size (63 bits), subgroup size (20 bits, approximately 2¹⁹), and benchmark suitability.
The external full orders are discovery provenance, not certified by the subgroup tests.
Resume from row 8 with the revised PLAN after review; `may-reshard` is not needed for this approved
change and was not granted. No implementation, tests or commits were performed for this revision.

### 8 — 2026-11-10 (implemented)
Discovery/flex: The destructive-HALT condition did not fire: c·P₀ = G′ and 561103·G′ = ∞ both hold
in-tree. `generic::N` is the order constant; the worker also exposed public constants for the field
prime, G′, P₀ and c in `generic.rs`. Brent recovers k = 424242 (seed 7) from (G′, k·G′); the CLI
example `--curve generic --q 1000932015531192845,9048581680910854061` was run in release and verified.
P₀ as base → `WrongGeneratorOrder`, P₀ as target → `TargetOutsideSubgroup` (tested); no projection of
supplied points.
Affected: C-GenericSubgroup (frozen at 8); no downstream contract flex.
Deferred: yes — `docs/BENCHMARKS.md` :78 and `docs/ROADMAP.md` :195 still cite the old composite order
9 223 372 037 218 517 353 and "group order is not prime". Row 9 rewrites the BENCHMARKS `generic` row
and footnote (subgroup ≈ 2¹⁹ over a 63-bit field; keep it distinct from field size and from the
external full-order provenance); row 10's sweep should re-read BENCHMARKS and PEDAGOGY :877 (§6 was
edited by row 8). ROADMAP is a rolling doc outside the durable-file list: note, do not edit unless the
row-10 scope covers it.
Texture: the `check_scalar_mul` test helper now takes the base point as a parameter.

### 9 — 2026-11-10
Discovery/flex: `generic` (≈ 2¹⁹ subgroup) was not reused; a new fixture `prime_order_25`
(y² = x³ − 3x + 72 over p = 16_777_259, whole-group prime order n = 16_777_729, G = (1, 6_299_712);
PARI/GP `ellcard` search a = −3, b = 1..200, first prime #E ≠ p at b = 72) serves every bench, with
the prime-order and Hasse-bound test in `ecdlp_kat.rs`. `bench_glv_vs_batch` was dropped: GLV exists
only on `secp_k1_toy` (≈ 10⁹ steps per solve), and the new fixture has no endomorphism. Tuning: θ = 4,
step budget 2²², Brent 20 attempts, batch 16, fixed seed, target `12_345_678_901 mod n`.
Affected: C-BenchFixture (frozen at 9); C-DpBudget and C-GenericSubgroup unchanged.
Deferred: yes — the worker found a SOLVER-BEHAVIOUR anomaly outside any row's scope (no solver edited):
`solve_dp_negmap` needs 3·10⁴–7·10⁴ steps on the 25-bit group vs ≈ 2·10³–10⁴ for `solve_dp`, did not
finish in 60–120 s on a 31-bit group (p = 2_147_483_659, n = 2_147_490_413) for θ ∈ {1,2,3,4,6,8}
while Brent/`solve_dp`/`solve_dp_batch` finished in seconds, and 1 of 6 seeds on the 25-bit group timed
out at 30 s (step counts cluster at multiples of the walker-restart threshold: walkers look trapped in
cycles). `solve_dp_batch` at 31 bits also took 1.1M steps at θ = 4 and 6.3M at θ = 6 vs ≈ 6·10⁴
expected. No contract is contradicted (budgets still bound every call; outcomes stay typed), so this
is an additive-reshard-class finding to SURFACE to the user at the end, never absorbed: negation-map
cycle handling is a new-algorithm concern outside this arc's scope. The ◆ fork and sweep should make
sure no durable text claims a negation-map speedup or says the DP solvers perform at the expected
cost; BENCHMARKS :502 already records the observation. BENCHMARKS :78 `generic` row and footnote were
rewritten (561_103 subgroup, in-tree proof, PARI/GP orders as uncertified provenance); `docs/ROADMAP.md`
:195 still cites the old composite order (rolling doc, note only).
Texture: the worker ran a 15-minute full bench on a 40-bit group while sizing before moving smaller
(bounded by timeout, no stray process); sizing examples were deleted; over-width lines added are
BENCHMARKS table rows only.

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
  decision). **Fired at row 8 (2026-11-10):** n′·G′ ≠ ∞ and `GENERIC_N` is not even a multiple of
  ord P₀. **Resolved by user approval 2026-11-10:** row 8 now uses the 561103-order subgroup;
  implementation and in-tree proof remain pending. If the specified c·P₀ or n′·G′ checks fail:
  destructive-HALT, not another constant adjustment. See the row-8 approval digest.
- **Inferred, enrolled 5:** Z ≡ 0 modulo a factor persists through later ladder calls, so the
  final gcd still detects it; if stage 2 depends on stage 1's point semantics, internal-continue.
- **Surfaced at row 9 (additive-reshard class, not absorbed):** `solve_dp_negmap` and `solve_dp_batch`
  take far more steps than expected on 25–31-bit groups, and `solve_dp_negmap` did not finish in
  60–120 s at 31 bits (see the row-9 digest). Needs its own investigation outside this arc.
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
