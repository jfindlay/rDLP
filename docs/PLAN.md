# rDLP — Shor ECDLP and elliptic index-calculus fidelity

Execution graph for the ROADMAP scope "Shor and elliptic index-calculus fidelity". The previous
(completed) Dedekind plan is preserved in git at `517f9a3:docs/PLAN.md`.

## Purpose (design intent)

Make the two fallback-masked solvers succeed by their advertised mechanism or fail with a typed,
honest outcome. Today `shor::ecdlp::solve_ecdlp` and `rho::index_calculus::index_calculus_dlp` both
end in an exhaustive search over the answer space and report the result as ordinary success, so
every end-to-end test passes whether or not the primary mechanism works. Static reading suggests
it does not: the Shor extraction applies the subgroup relation a′+b′k≡0 where measured outcomes
satisfy the dual relation d≡k·c, without rescaling the M=16 register to r=13; the index-calculus
kernel is a right kernel over factor-base columns whose result is discarded.

Follow the repository's GNFS failure contract: an alternative is admissible only if it is itself a
published method through the same exact verifier; resource exhaustion, mechanism negatives,
invariant failures and verified results are distinct; no unverified answer is success. Here that
means: remove both exhaustive fallbacks; Shor retries with fresh measurements (Shor 1997 §6), and
index calculus collects more relations (standard index-calculus practice) within a bound.

Scope boundary: the two solvers, their tests, benches and the exposition of these two mechanisms.
Shor factoring's register width, `order_from_phase` minimality and tuned seeds are ROADMAP
findings, not this arc. No new algorithms, no crate edges, no extension-field index calculus.

## Verify gate

- VERIFY_TYPES: `cargo check --workspace`
- VERIFY_TEST: `cargo test --workspace`
- Also mandatory: `cargo dep-policy`, `cargo audit`, `cargo coverage` (80% line floor).
- Not gates: `cargo fmt-check`, `cargo lint` (unadopted; do not fix drive-by).
- Cost-gated (not oracles): `cargo test -p shor --test ecdlp_kat -- --ignored` runs the all-k
  sweeps; run it at rows 2 and 4. Circuit runs cost ~10 s each in debug.

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 | Recover the ECDLP exponent from finite-precision Fourier outcomes | B | T1 | (existing circuit) | `shor/src/ecdlp/mod.rs`, `shor/tests/ecdlp_kat.rs` |
| 2 | Report Shor ECDLP measurement evidence without exhaustive fallback | B | T1 | C-ShorExtraction | `shor/src/ecdlp/mod.rs`, `shor/tests/ecdlp_kat.rs` |
| 3 | Recover the index-calculus log from a relation combination | B | T1 | (existing relations) | `rho/src/index_calculus/{linalg,solve,collect,mod}.rs`, `rho/tests/index_calculus_kat.rs`, `rho/benches/attacks.rs` |
| 4 ◆ | Correct the Shor ECDLP and index-calculus exposition | I | T1 (◆ fork T0) | C-ShorExtraction, C-ShorSolve, C-IndexCalcRecovery | `docs/MATHEMATICS.md`, `docs/PEDAGOGY.md`, `shor/docs/PEDAGOGY.md`, `docs/BENCHMARKS.md` |
| 5 | Anneal cleanup (skippable) | I | — | none | durable files touched by rows 1–4 |

Rows 1–2 (shor) and 3 (rho) are independent; run in row order. No `@architect` inflection rows:
the interfaces are frozen below and the only judgment-heavy derivation (row 1) is checked by an
exact, independent distribution KAT.

## Session detail

### 1 — Recover the ECDLP exponent from finite-precision Fourier outcomes

Deliverable: replace `extract_k_from_measurement` with the finite-precision dual-lattice rule and
add the exact outcome distribution, per C-ShorExtraction. Derivation to implement and state in
rustdoc: after H on both registers the state is Σ_{a,b<M} |a,b⟩|(a+bk mod r)·G⟩; outcome (c,d) has
probability P(c,d) = Σ_s |M⁻² Σ_{a+bk≡s} ω^{∓(ac+bd)}|², ω=e^{2πi/M}. At c/M=j/r, d/M=l/r the phase
is constructive iff l≡kj (mod r), independent of the transform's sign. So j=round(r·c/M) mod r,
l=round(r·d/M) mod r, and k=l·j⁻¹ mod r when j≢0. Exactly one candidate per measurement.

`solve_ecdlp` keeps its fallback in this row (row 2 removes it) but uses the new extraction.
Translate the module's leaked labels (`C-PointAdd`, `C-ECDLPSolve`) into plain prose.

KATs:
- **Independent distribution check:** `exponent_outcome_distribution` (from circuit amplitudes,
  marginalised over the point register) equals the classical Fourier-sum formula above within
  1e-9, sums to 1, for k ∈ {3, 7} on the default path and every k ∈ 2..13 under
  `#[ignore = "cost: ..."]`. A mismatch that vanishes under bit reversal of c or d is a register
  convention defect: decide against `qft` module KATs, fix the convention, record in the digest.
- `primary_success_probability(k)` > 0 for every tested k; record the values in the digest.
  Compute the old rule's success probability once and record it (do not test it).
- Synthetic extraction: exact peak pairs recover k; j≡0 → `None`; inputs ≥ M rejected.

Subtleties: never widen the candidate set (neighbouring roundings) to raise success probability.
At r=13, nine candidates would be almost exhaustive search. If the probabilities look too low for
a sane retry budget, raising `EXP_BITS` is internal-continue only if `qubit_budget_within_ceiling`
holds; record the cost change. Cache the distribution per k (e.g. `OnceLock`) in tests.

### 2 — Report Shor ECDLP measurement evidence without exhaustive fallback

Deliverable: `solve_ecdlp`/`solve_ecdlp_with_retries` return `Result<EcdlpSolution, EcdlpError>`
per C-ShorSolve. Delete the brute-force loop and its comment ("the measurement biases the search").
Each attempt uses a fresh seed; every accepted k is verified (k·G == Q) and carries the measured pair.
Correct the module and function rustdoc. Replace the vacuous end-to-end tests: each now asserts
`Measurement` evidence and `extract_k_from_measurement(a′, b′) == Some(k)`. Adapt the
`rho_crosscheck_*` tests to the new API. They stay ignored for cost: state the cost reason, and do
not call them oracle tests.

KATs: k ∈ {3, 5, 7, 9, 12} at fixed seeds → `Measurement` evidence consistent with the extraction;
Q=∞ and Q=G → `Trivial`; a budget and seed known to fail → `AttemptsExhausted { attempts }`, the
exact count; G=∞ → `InvalidGenerator`; all-k sweep (ignored, cost) succeeds within the default budget.

Subtleties: choose the default budget from row 1's recorded probabilities. Seeds may need
re-selection: re-select by a bounded recorded search, never by injecting outcomes.

### 3 — Recover the index-calculus log from a relation combination

Deliverable: per C-IndexCalcRecovery. `solve_ek_linalg` builds the transposed matrix (rows =
factor-base points, columns = relations) and returns the full left-kernel basis from the existing
`block_wiedemann_fl` engine: combinations u with Σᵢ uᵢeᵢ ≡ 0. From relations aᵢ+bᵢx ≡ Σⱼ eᵢⱼlⱼ
(mod ℓ), any u with uᵀb ≢ 0 gives x ≡ −uᵀa/uᵀb.

`index_calculus_dlp` collects relations, finds such a u, and checks the following:
- **Consistency:** every basis vector with uᵀb ≢ 0 yields the same x.
- **Verification:** x·G_ℓ == Q_ℓ.

It returns the solution with the u used. If no basis vector has uᵀb ≢ 0, it collects further
relations deterministically up to the full (a,b) enumeration bound. Then:
- exhaustion of that bound → `Undetermined`;
- disagreement or a failed verification → `Inconsistent` (an invariant failure, never retried).

Delete Phase 2 enumeration, the local `gauss_elim_f_ell` and its helpers, the unused
`RecoveryFailed`. Adapt `rho/benches/attacks.rs` to the new type. Correct the rustdoc ("kernel +
relation provenance", "consistency check") and translate `Principle-4` labels, including the
`principle4_annotation` test name, into the stated property: over E(F_p) the mechanism is
demonstrated, not the asymptotic advantage.

KATs:
- Left-kernel vectors annihilate the exponent matrix (replaces `kernel_correctness`).
- `adapter_fidelity` holds for the transposed layout.
- For every k ∈ 0..60, x = k mod 5 with uᵀb ≢ 0 and the identity Σuᵢ(aᵢ+bᵢx) ≡ 0 checked
  explicitly. This is cheap: keep it on the default path.
- A hand-built relation set with uᵀb ≡ 0 for all u → `Undetermined`.
- A corrupted relation → `Inconsistent`.
- Delete `index_calculus_dlp_smoke` (vacuous).

Subtleties: logs reduce mod ℓ because ℓ ∥ n=60, so the ℓ-Sylow subgroup is cyclic and every point's
ℓ-component lies in ⟨G_ℓ⟩. State this as the prose invariant. The six-file footprint is one
mechanism (kernel orientation, recovery formula, outcome type are inseparable): the irreducible
floor, not a split candidate. If `block_wiedemann_fl` does not return a full basis, the
`Undetermined` claim is not certified: compute the basis completely or rename the outcome to a
resource outcome, and record in the digest.

### 4 ◆ — Correct the Shor ECDLP and index-calculus exposition

Deliverable: correct the human docs against the frozen contracts.
- MATHEMATICS's Shor ECDLP section (~3954–4015): finite-precision derivation, the dual relation,
  all-qubit measurement, the stale `shor/src/ecdlp.rs` path, no general continued-fraction claim.
- MATHEMATICS's index-calculus section (~3209–3248): left-kernel recovery, the ℓ-Sylow invariant.
- PEDAGOGY's index-calculus tour (~1088–1167, ~1235–1237; `decompose`'s actual return type).
- shor PEDAGOGY (~304–308, ~463).
- BENCHMARKS (~473–514, ~699–733): drop "measurement biases" and success claims that the
  fallback produced.

KAT substitute (the row is prose): every quantitative claim published — success probabilities,
outcome masses, x-recovery coverage — must cite the test that computes it; a claim no test computes
is removed or labelled as an estimate. Flagged as the row's contract.

Integrative duty: sweep superseded claims (measurement-biased search, kernel provenance, Phase 2,
"KAT-verified … Complete", exhaustive-is-honest) across all durable surfaces these two solvers
touch, including unchanged passages. Anneal the arc's durable files against the denylist.

### 5 — Anneal cleanup (skippable)

`/plan-anneal` over the arc's durable files: translate leaked plan coordinates and rewrap ledgered
line-length overruns. The driver may mark it `skipped` if the ◆ anneal ran clean and no overruns
were ledgered.

## Cross-session contracts

### C-ShorExtraction

Compiler + test-enforced. Defined in 1. Consumed by 2, 4.

```rust
/// Candidate k from one measured exponent pair (c, d) ∈ [0, 2^EXP_BITS)²:
/// j = round(r·c/M) mod r, l = round(r·d/M) mod r, k = l·j⁻¹ mod r; `None` iff j ≡ 0.
pub fn extract_k_from_measurement(c: u64, d: u64, r: u64) -> Option<u64>;
/// Exact P(c, d) from circuit amplitudes, marginalised over the point register;
/// index `c + (d << EXP_BITS)`, length 2^(2·EXP_BITS).
pub fn exponent_outcome_distribution(g: Point, q: Point) -> Vec<f64>;
/// Σ P(c, d) over outcomes whose extraction returns log_g(q).
pub fn primary_success_probability(g: Point, q: Point) -> f64;
```

Prose invariant: one candidate per measurement; no widening of the candidate set.

### C-ShorSolve

Compiler + test-enforced. Defined in 2. Consumed by 4.

```rust
pub fn solve_ecdlp(g: Point, q: Point, seed: u64) -> Result<EcdlpSolution, EcdlpError>;
pub fn solve_ecdlp_with_retries(g: Point, q: Point, seed: u64, max_attempts: usize)
    -> Result<EcdlpSolution, EcdlpError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EcdlpSolution { pub k: u64, pub evidence: EcdlpEvidence }
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EcdlpEvidence {
    Trivial,                                                   // Q = ∞ or Q = G
    Measurement { seed: u64, attempt: usize, c: u64, d: u64 }, // extraction(c, d) == Some(k)
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EcdlpError {
    InvalidGenerator,                     // G = ∞
    AttemptsExhausted { attempts: usize } // resource outcome, not a mathematical negative
}
```

Prose invariant: every `Ok` has passed k·G == Q; no code path enumerates candidate exponents.

### C-IndexCalcRecovery

Compiler + test + prose-enforced. Defined in 3. Consumed by 4.

```rust
/// Left-kernel basis: each u ∈ F_ℓ^{#relations} satisfies Σᵢ uᵢ·eᵢ ≡ 0.
pub fn solve_ek_linalg(relations: &[Relation], strategy: &IndexCalcStrategy)
    -> Result<Vec<Vec<FpNaive>>, IndexCalcError>;
pub fn index_calculus_dlp(g, q, strategy: &IndexCalcStrategy)
    -> Result<IndexCalcSolution, IndexCalcError>;

pub struct IndexCalcSolution {
    pub x: u64,                      // log_G(Q) mod ℓ, verified x·G_ℓ == Q_ℓ
    pub combination: Vec<FpNaive>,   // the u used; uᵀb ≢ 0
    pub relations: Vec<Relation>,    // the system u combines
}
// IndexCalcError gains Undetermined { relations: usize } (bound exhausted, no u with uᵀb ≢ 0)
// and Inconsistent (combinations disagree or verification fails); RecoveryFailed is removed.
```

Prose invariants: relations hold mod ℓ because ℓ ∥ n; x comes only from a relation combination;
`Inconsistent` is an invariant failure, never retried or masked.

## Progress ledger

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Recover the ECDLP exponent from finite-precision Fourier outcomes | pending | | |
| 2 | Report Shor ECDLP measurement evidence without exhaustive fallback | pending | | |
| 3 | Recover the index-calculus log from a relation combination | pending | | |
| 4 | Correct the Shor ECDLP and index-calculus exposition ◆ | pending | | |
| 5 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

*(none yet)*

## Discoveries & risks

- **Inferred (sharding), enrolled 1:** wrong extraction orientation and no r/M rescaling. If the
  exact distribution shows the old rule was in fact correct under this circuit's conventions:
  internal-continue (the distribution is the contract, not the inference); record it.
- **Inferred, enrolled 3:** the kernel is a right kernel over factor-base columns; Phase 1 neither
  detects inconsistent rows nor verifies x.
- **Confirmed, enrolled 2–3:** every end-to-end test of both solvers would pass on the fallback
  alone; `index_calculus_dlp_smoke` asserts only `is_ok()`.
- **Internal-continue:** `EXP_BITS` raise within the qubit ceiling; seed re-selection; relation
  collection beyond fb_size+1 within the (a,b) bound; full-basis computation if the engine's is
  partial.
- **Additive-reshard (surface, do not absorb):** any edit outside `shor/src/ecdlp`, `rho/src/
  index_calculus`, their tests, `rho/benches/attacks.rs` and the row-4 docs; Shor factoring
  fidelity; deduplicating `iqft_on_qubits`; a new crate edge.
- **Destructive-HALT:** the circuit distribution disagrees with the Fourier-sum formula in a way
  no register convention explains; or the primary success probability is 0 for some k under every
  admissible `EXP_BITS` (the mechanism cannot work at this scale, a scope decision for the user).
- **Preflight:** PLAN and ROADMAP will be dirty after sharding; commit them before `/plan-run`.

## Notes for executors

**Sizing.**
- (1) Ambient complexity is moderate: two small modules in different crates.
- (2) The floor is moderate. Row 3's mechanism is atomic; row 1's derivation is one object.
- (3) Cost of a design error is moderate: public API breaks, but no consumers exist outside the
  tests and one bench.
- (4) Correctness-criticality is high: the arc exists because evidence lied. This drove the
  1|2 split at the extraction contract.
- (5) Inner-loop bandwidth is high for rho and for Shor's distribution check, but low for Shor
  end-to-end runs (~10 s each). Keep default-path circuit runs few.

No row uses F/C: the rows are small and tightly steered.

**Routing.** T1 implements every row; T0 runs only the ◆ boundary fork (intent check, claim
adjudication).

**Outcome contract.** `docs/development.md`, "Solver outcomes and fallbacks", is the authority for
rows 1–3; the Purpose section applies it to these two solvers.

**Register.** State the property or invariant, never the plan coordinate. Mathematical exposition
goes in the sections named in row 4; no new chapters.

**Anneal denylist** (grep enrolled durable files at ◆, excluding PLAN/ROADMAP/NOTES/SURVEY/AGENTS):

    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat [ABCI]\b
    \bC-[A-Z][A-Za-z]+\b
    \b[A-H]\.[A-Z0-9]+\b|[Pp]rinciple-?[0-9]+
    \bTrack [A-Z]\b|\bPhase [0-9]+\b
    (PLAN|ROADMAP|NOTES|SURVEY)(\.md)?\b|plan-(run|shard|juncture)

KAT, algorithmic step numbers, group notation (S₃, S4) and quantum phase terminology are domain
terms: classify collisions, do not rename.

**First run:** after committing the planning docs, from `@plan-admin`:

    /plan-run docs/PLAN.md halt-at-boundaries
