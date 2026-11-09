# rDLP — Solver-outcome conformance

Execution graph for the ROADMAP scope "Solver-outcome audit". The previous (completed) Shor ECDLP
and index-calculus plan is preserved in git at `01b4df9:docs/PLAN.md`.

## Purpose (design intent)

Bring every remaining solver under `docs/development.md`'s "Solver outcomes and fallbacks": success
returned only after an exact check independent of the method; invalid input, unsupported domain,
resource exhaustion, invariant failure and mathematical negatives kept distinct; no unbounded
search; provenance whenever more than one path can succeed. A static audit at sharding found
nearly every rho solver returning an unchecked `Option<u64>` whose `None` merges several failure
classes, distinguished-point walks with no budget, and one probable wrong-answer path (plain
distinguished-point collisions keyed by x alone). The repaired Shor ECDLP
(`shor::ecdlp::{EcdlpSolution, EcdlpError}`) is the reference shape.

Scope boundary: the rho ECDLP family (Brent, distinguished-point variants, binary-field, Koblitz),
Pohlig–Hellman, MOV, SSA, the GHS transfer's input guard, rho factoring, ECM, Shor factoring, their
tests, benches, CLIs and exposition. Out of scope: NFS-DL (including MOV's k=2 log engine), Shor
factoring's register width (documented, not changed), new algorithms, new crate edges.

## Verify gate

- VERIFY_TYPES: `cargo check --workspace`
- VERIFY_TEST: `cargo test --workspace`
- Also mandatory: `cargo dep-policy`, `cargo audit`, `cargo coverage` (80% line floor).
- Not gates: `cargo fmt-check`, `cargo lint` (unadopted; do not fix drive-by).
- Cost-gated (not oracles): `cargo test -p shor --test ecdlp_kat -- --ignored` at rows 1–2 (the
  rho cross-checks call the changed API); `cargo test -p shor --test factor_kat -- --ignored` at
  row 11.

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 @architect | Define a verified rho discrete-log outcome and adopt it in Brent rho | A | T1 (juncture T0) | (existing walks) | `rho/src/ecdlp/{outcome.rs (new),mod.rs,cli.rs,pohlig.rs}`, `rho/tests/{ecdlp_kat,ssa_kat}.rs`, `rho/benches/ecdlp.rs`, `shor/tests/ecdlp_kat.rs` |
| 2 | Bound the distinguished-point rho solvers and verify their results | B | T1 | C-DlogOutcome | `rho/src/ecdlp/{mod.rs,coordinator.rs,dp.rs,cli.rs}`, `rho/tests/ecdlp_kat.rs`, `rho/benches/ecdlp.rs`, `shor/tests/ecdlp_kat.rs` |
| 3 | Recover the log from negated distinguished-point collisions | B | T1 | C-DlogOutcome, C-DpBudget | `rho/src/ecdlp/{coordinator.rs,dp.rs,mod.rs}`, `rho/tests/ecdlp_kat.rs` |
| 4 | Verify Pohlig–Hellman results and name the small-subgroup method | B | T1 | C-DlogOutcome | `rho/src/ecdlp/pohlig.rs`, `rho/src/curve/test_curves.rs`, `rho/tests/{ecdlp_kat,index_calculus_kat}.rs`, `rho/benches/attacks.rs`, `docs/development.md` |
| 5 | Verify binary-field rho and Koblitz results | B | T1 | C-DlogOutcome | `rho/src/binary_ecdlp/{mod.rs,koblitz.rs}`, `rho/tests/binary_ecdlp_kat.rs` |
| 6 | Verify MOV results and type its failure classes | B | T1 | C-OutcomeShape | `rho/src/pairing/mov.rs`, `rho/tests/mov_kat.rs`, `rho/benches/attacks.rs` |
| 7 | Separate SSA failure classes and handle zero logs and large primes | B | T1 | C-OutcomeShape | `rho/src/ssa/{mod.rs,reduce.rs}`, `rho/tests/ssa_kat.rs`, `rho/benches/attacks.rs` |
| 8 | Reject off-curve and non-subfield points in the GHS transfer | B | T1 | C-OutcomeShape | `rho/src/ghs/{transfer.rs,mod.rs,reduce.rs}`, `rho/tests/ghs_kat.rs`, `rho/benches/attacks.rs` |
| 9 | Type rho factoring outcomes and bound each attempt | B | T1 | C-OutcomeShape | `rho/src/factor/{rho.rs,mod.rs,cli.rs}`, `rho/tests/factor_kat.rs`, `rho/benches/factor.rs` |
| 10 | Type ECM factoring outcomes | B | T1 | C-OutcomeShape | `shared/numth/src/{ecm.rs,lib.rs}` |
| 11 | Report Shor factoring provenance and return the minimal order | B | T1 | C-OutcomeShape | `shor/src/shor/mod.rs`, `shor/src/lib.rs`, `shor/tests/factor_kat.rs` |
| 12 ◆ | Correct the solver-outcome exposition | I | T1 (◆ fork T0) | all | `docs/MATHEMATICS.md`, `docs/PEDAGOGY.md`, `docs/BENCHMARKS.md`, `README.md`, `shor/docs/PEDAGOGY.md`, `docs/development.md` |
| 13 | Anneal cleanup (skippable) | I | — | none | durable files touched by rows 1–12 |

Rows 1→2→3 are sequential (shared walk and coordinator code); 4 and 5 consume row 1's frozen
contract; 6–11 are mutually independent and follow row order. Row 1 is the inflection point: the
juncture fork checks the drafted C-DlogOutcome against the walk, coordinator, GLV and Koblitz code
before row 1 implements it.

## Session detail

### 1 @architect — Define a verified rho discrete-log outcome and adopt it in Brent rho

Deliverable: `rho/src/ecdlp/outcome.rs` with C-DlogOutcome (re-exported from `rho::ecdlp`);
`solve_brent` returns `Result<DlogSolution<RhoEvidence>, DlogError>`. Validate input first
(`check_prime_order_input`: G ≠ ∞, G and Q on the curve, n prime via `shared-numth`, n·G = ∞,
n·Q = ∞), handle Q = ∞ and Q = G as `Trivial`, and verify every collision-derived k with
`verify_log` (k·G == Q) before returning. Degenerate collisions consume retries; spending all
retries or hitting the per-attempt step cap is `BudgetExhausted` (fix the rustdoc that calls the
cap "all retries were degenerate"). A non-degenerate collision that fails verification is
`Inconsistent`, never retried. The `rho-dlog` CLI drops its own re-check and reports each error
class distinctly.

The footprint is one atomic unit (irreducible floor): the signature change forces every
`solve_brent` caller in the same commit — `pohlig.rs` (interim adaptation only; row 4 rewrites
it), `ssa_kat.rs`, `shor/tests/ecdlp_kat.rs`, the bench and the CLI.

KATs:
- Fixed-seed fixtures return `Collision` evidence whose identity
  a1·G + b1·Q = ±(a2·G + b2·Q) the test rechecks from the evidence alone.
- Q = ∞ → `Trivial` with k = 0; Q = G → `Trivial` with k = 1.
- Each `InputDefect` variant is reachable (off-curve point, composite n, wrong n, G = ∞).
- `max_retries = 0` → `BudgetExhausted { attempts: 0, .. }`.
- `verify_log` rejects k + 1.

Subtleties: `check_prime_order_input` costs two scalar multiplications plus a primality test per
call; fine at these sizes. n·Q = ∞ does not prove Q ∈ ⟨G⟩ when n² divides #E; the residual case
surfaces as `Inconsistent`, and the rustdoc says so.

### 2 — Bound the distinguished-point rho solvers and verify their results

Deliverable per C-DpBudget: `solve_dp`, `solve_dp_negmap`, `solve_dp_batch` and `solve_dp_glv`
return the C-DlogOutcome types and take an explicit total-step budget. Walkers stop and the call
returns `BudgetExhausted` when it is spent, so no input can hang the solver. The coordinator
returns collision evidence, not a bare k; the solver verifies. `solve_dp_glv` returns
`Unsupported` off the curve its β/λ constants belong to (check the curve parameters, not a name).
Walker panic or channel loss is an invariant failure reported as such, not `None`. Delete the
"should not happen" `None` paths. The single-walker negmap limitation becomes a budget outcome.

KATs: the existing fixed-seed success fixtures assert evidence identities; a tiny budget →
`BudgetExhausted` for each variant; GLV on a non-GLV curve → `Unsupported`; replace the
wall-clock test `glv_fewer_steps_than_negmap` with a step-count comparison if the evidence
carries step counts (otherwise leave it and record).

Subtleties: pick default-path budgets with headroom so fixed-seed tests cannot become flaky;
record the chosen values in the digest.

### 3 — Recover the log from negated distinguished-point collisions

Deliverable per C-DpSignedCollision. `DpRecord` carries enough to tell P from −P (a y-parity bit)
and the full x-coordinate (not only the low word). A collision at equal x with opposite sign
yields k from a1 + b1·k ≡ −(a2 + b2·k) and evidence `negated: true`; unequal full x is not a
collision.

KATs:
- Synthetic coordinator test: records for a known P and −P with hand-chosen (a, b) recover the
  right k (deterministic; the primary evidence).
- A bounded, recorded seed search on a test curve for a plain `solve_dp` run that hits a negated
  collision; if found, before the fix it yields `Inconsistent` (record that) and after the fix
  `Ok` with `negated: true`.
- Two points sharing their low x word but differing in full x do not collide.

Subtleties: the bug is inferred from reading the code. If the plain walk provably never meets −P,
keep the synthetic test, record the reason, continue.

### 4 — Verify Pohlig–Hellman results and name the small-subgroup method

Deliverable per C-PohligEvidence. First read the primary source (Pohlig–Hellman 1978; the
Silver–Pohlig–Hellman table formulation) and record whether the order-p sub-step there is a
table/enumeration over the order-p subgroup. If yes, name the sub-method after the source, keep
the `p ≤ 64` size dispatch as the method's declared sub-step, and add one clarifying sentence to
`docs/development.md` rule 5: a size-selected sub-step that is itself the published method's step
is not a fallback. Return `DlogSolution<PohligEvidence>` with the per-digit method recorded; verify
k·G == Q after CRT (`Inconsistent` otherwise); check n·G = ∞ (`InvalidInput`); `n = 0` no longer
panics. Add a composite-order test curve with a prime factor > 64 so the rho branch runs
end to end. Adapt `attacks.rs` and the oracle use in `index_calculus_kat.rs`.

KATs: every k on the new fixture (or a sampled set if costly) with evidence showing both
sub-methods; the identity residue ≡ k (mod p^e) for each component; existing `composite_toy`
cases; wrong n → `InvalidInput`.

Subtleties: if the literature does not support the claim, HALT (see Discoveries). The fixture's
order must be factorable by the existing `factor_order` (u64 trial division).

### 5 — Verify binary-field rho and Koblitz results

Deliverable: `binary_ecdlp::{solve_brent, solve}` and `koblitz::solve_koblitz` return
`Result<DlogSolution<RhoEvidence>, DlogError>`, with binary-curve input checks mirroring row 1 and
k·G == Q verified. Koblitz evidence carries the τ-power (`lambda_power`) and sign of the class
equivalence. Remove or use the unused `_t` parameter. Budget and degenerate retries as in row 1.

KATs: existing binary and Koblitz fixtures assert rechecked evidence identities; the τ eigenvalue
from `compute_tau_eigenvalues` satisfies τ(P) = λ·P on the fixture; a corrupted eigenvalue gives
`Inconsistent`, not a wrong `Ok`.

### 6 — Verify MOV results and type its failure classes

Deliverable per C-MovOutcome: `mov_reduce` returns `Result<MovSolution, MovError>`; verify
k·G == Q in E[ℓ] (independent of the pairing) before returning; the degenerate pairing becomes
`DegeneratePairing` instead of a fabricated `DescentFailed { stuck_prime: 0 }`; check ℓ prime,
ℓ | p^k − 1 and G, Q ∈ E[ℓ] as `InvalidInput`; carry the pairing values as evidence. The k=2 log
engine stays as it is (NFS-DL scope); its errors pass through as `MovError::Dl`.

KATs: existing k ∈ {1, 2} fixtures assert evidence (g_pair^k == h_pair) and verification;
degenerate R → `DegeneratePairing`; a point outside E[ℓ] → `InvalidInput`.

### 7 — Separate SSA failure classes and handle zero logs and large primes

Deliverable per C-SsaOutcome: Q = ∞ returns `Ok(0)`; a verification failure after the curve was
proven anomalous is `Inconsistent`, not `NotAnomalous`; `verify_anomalous` counts exactly over
its supported range (u128 arithmetic) and returns `Unsupported` beyond it rather than a false
negative. State the mathematical source of the k_raw versus k_raw/2 correction in the rustdoc; if
both candidates can succeed, record which one did (rule 6).

KATs: Q = ∞ → `Ok(0)`; a large-p anomalous curve no longer reports `NotAnomalous` (or reports
`Unsupported`, per the chosen range); the existing exact-k and non-anomalous cases.

Subtleties: if the factor-2 correction turns out to mask a lifting defect rather than follow from
the formal-log derivation, surface it (additive), do not absorb it.

### 8 — Reject off-curve and non-subfield points in the GHS transfer

Deliverable per C-GhsTransferGuard: `transfer_point` returns `Result<_, GhsError>`, checks the
input is on E and that the conorm coordinates lie in GF(2^l) before reinterpreting them; remove
`GhsError::PointAtInfinity` if it remains unreachable. Adapt `ghs_descend`, `ghs_kat.rs` and the
bench.

KATs: off-curve input → `NotOnCurve`; a non-Galois-invariant conorm → `NotInSubfield`; existing
homomorphism and log-preservation KATs unchanged.

### 9 — Type rho factoring outcomes and bound each attempt

Deliverable per C-RhoFactorOutcome: `factor` returns `Result<Factorization, FactorError>`, checks
1 < d < n and n % d == 0 independently, carries `Even` or `Rho { c }` evidence, reports n ≤ 1 as
`InvalidInput`, a prime n as `Prime` only where the primality test is proven on the domain,
n ≥ 2^80 as `Unsupported` (the `mulmod` domain), and spent constants as `BudgetExhausted`.
`factor(2)` no longer returns 2. The primitives (`floyd`, `brent`, `brent_batched`) take an
explicit step bound and keep a documented single-failure `Option` ("this attempt failed"). The
`rho-factor` CLI distinguishes prime from exhaustion.

KATs: existing composites with evidence; n ∈ {0, 1} → `InvalidInput`; n ∈ {2, 3, a large prime} →
`Prime`; n ≥ 2^80 → `Unsupported`; `max_c = 0` → `BudgetExhausted`.

### 10 — Type ECM factoring outcomes

Deliverable per C-EcmOutcome: `ecm_factor` returns `Result<Uint<4>, EcmError>` with
`InvalidInput`, `Prime` (only where proven) and `CurvesExhausted { curves }` kept distinct, and
checks the divisor independently; n = 2 is no longer a trivial-divisor success. Update the
`shared/numth` re-export. Inline tests live in `ecm.rs`.

KATs: a semiprime factors with a checked divisor; n = 2 → `Prime`; `max_curves = 0` →
`CurvesExhausted { curves: 0 }`.

### 11 — Report Shor factoring provenance and return the minimal order

Deliverable per C-ShorFactorOutcome: `factor` returns `Result<ShorFactorization, ShorFactorError>`
with `Even`, `PrimePower`, `CommonFactor` (if the driver takes that path) or
`OrderFinding { a, seed, s, r }` evidence, p·q == n checked; prime n → `NotComposite` before any
circuit runs; inputs beyond the simulator ceiling → `Unsupported`; spent budget →
`AttemptsExhausted`. `order_from_phase` returns the multiplicative order: reduce a convergent
denominator d with a^d ≡ 1 by dividing out prime factors while the identity holds. Fix its rustdoc
(`≤ N` bound; it returned multiples). Correct the stale seed comments (s = 2389 cannot occur in a
6-bit register) and make `find_order_returns_valid_order` non-vacuous. The register width stays
n_bits(N); the rustdoc states the textbook 2^t ≥ N² bound and that sufficiency here is empirical.

KATs: each fixture asserts `OrderFinding` evidence and rechecks a^r ≡ 1 and the gcd step;
`cf_order_from_phase_n15_a7_s32` now expects 4; prime n → `NotComposite` without a circuit run;
even n → `Even`.

Subtleties: circuit runs are slow; keep default-path additions few.

### 12 ◆ — Correct the solver-outcome exposition

Deliverable: correct human docs against the frozen contracts: MATHEMATICS's rho, Pohlig–Hellman,
MOV, SSA and Shor-factoring passages; PEDAGOGY's rho/Pohlig/attack tours (including the
nonexistent `rho/tests/pohlig_kat.rs` reference near line 1361); BENCHMARKS success claims;
README CLI descriptions; shor PEDAGOGY's factoring chapter; check that `docs/development.md`'s
rule-5 clarification from row 4 reads as a standalone rule.

KAT substitute (the row is prose): every published behaviour claim (failure classes, budgets,
verification) cites the test that exercises it; a claim no test computes is removed or labelled.
Flagged as the row's contract.

Integrative duty: sweep superseded claims ("should not happen", "None if all retries were
degenerate", unverified "`Some(k)` with k·G = Q", "sufficient … verified for all four targets")
across all durable surfaces these solvers touch, including unchanged passages. Anneal the arc's
durable files against the denylist.

### 13 — Anneal cleanup (skippable)

`/plan-anneal` over the arc's durable files: translate leaked plan coordinates and rewrap ledgered
line-length overruns. The driver may mark it `skipped` if the ◆ anneal ran clean and no overruns
were ledgered.

## Cross-session contracts

### C-OutcomeShape

Prose + test-enforced. Defined here. Consumed by 1–11. Every repaired solver follows
`docs/development.md` rules 1–6 in the shape of `shor::ecdlp::{EcdlpSolution, EcdlpError}`: a
solution struct (answer plus evidence wherever more than one path can succeed) returned only after
an exact check independent of the method; a typed error whose variants separate invalid input,
unsupported domain, resource exhaustion, invariant failure (`Inconsistent`, never retried) and
proven negatives. A documented single-failure `Option` remains acceptable for primitives. No
crate edges are added to share types: rho, shor and shared-numth each own theirs.

### C-DlogOutcome

Compiler + test-enforced. *To be frozen at 1 (inflection juncture checks this draft against the
code).* Consumed by 2, 3, 4, 5, 12.

```rust
// rho/src/ecdlp/outcome.rs, re-exported from rho::ecdlp
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DlogSolution<E> { pub k: u64, pub evidence: E }  // k·G == Q checked by the solver

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DlogError {
    InvalidInput(InputDefect),
    Unsupported(&'static str),                         // outside the implementation's domain
    BudgetExhausted { attempts: usize, steps: u64 },   // resource outcome; proves nothing
    Inconsistent { k: u64 },                           // non-degenerate collision, k·G ≠ Q
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputDefect {
    GeneratorAtInfinity, PointNotOnCurve, OrderNotPrime, WrongGeneratorOrder,
    TargetOutsideSubgroup,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RhoEvidence {
    Trivial,  // Q = ∞ (k = 0) or Q = G (k = 1)
    // a1·G + b1·Q = (−1)^negated · λ^lambda_power · (a2·G + b2·Q); λ is the endomorphism
    // eigenvalue (GLV, Koblitz τ), lambda_power = 0 for plain walks.
    Collision { attempt: usize, a1: u64, b1: u64, a2: u64, b2: u64,
                negated: bool, lambda_power: u32 },
}

pub fn verify_log<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>, q: &AffinePoint<F>, k: u64) -> bool;
pub fn check_prime_order_input<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>,
    q: &AffinePoint<F>, n: u64) -> Result<(), DlogError>;
pub fn solve_brent<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>, q: &AffinePoint<F>, n: u64,
    seed: u64, max_retries: usize) -> Result<DlogSolution<RhoEvidence>, DlogError>;
```

Prose invariants: every `Ok` passed `verify_log`; `Inconsistent` is never retried or masked;
`lambda_power` and `negated` are over-specified for rows 2, 3 and 5.

### C-DpBudget

Compiler + test-enforced. *To be frozen at 2.* Consumed by 3, 12. The four distinguished-point
solvers return `Result<DlogSolution<RhoEvidence>, DlogError>` and take an explicit total-step
budget parameter; no input can make them run unboundedly; `solve_dp_glv` returns `Unsupported`
off its curve.

### C-DpSignedCollision

Test + prose-enforced. *To be frozen at 3.* Consumed by 12. A distinguished-point collision
requires equal full x; opposite y-parity yields the negated relation and `negated: true`.

### C-PohligEvidence

Compiler + test-enforced. *To be frozen at 4.* Consumed by 12.

```rust
pub fn solve_ecdlp_composite<F: Fp<4>>(curve: &Curve, g: &AffinePoint<F>, q: &AffinePoint<F>,
    n: u64) -> Result<DlogSolution<PohligEvidence>, DlogError>;
pub struct PohligEvidence { pub components: Vec<PrimePowerComponent> }
pub struct PrimePowerComponent { pub p: u64, pub e: u32, pub residue: u64,
                                 pub digits: Vec<DigitMethod> }
pub enum DigitMethod { SubgroupTable, Rho(RhoEvidence) }  // name per the primary source
```

Prose invariant: the small-subgroup step is the published method's own size-selected sub-step,
documented with its source; `docs/development.md` rule 5 says so in one sentence.

### C-MovOutcome

Compiler + test-enforced. *To be frozen at 6.* Consumed by 12. `mov_reduce -> Result<MovSolution,
MovError>`; `MovSolution { k_mod_ell, ell, g_pair, h_pair }` verified k·G == Q;
`MovError { InvalidInput, DegeneratePairing, Dl(gnfs::dl::SolveDlError), Inconsistent }`.

### C-SsaOutcome

Compiler + test-enforced. *To be frozen at 7.* Consumed by 12. `ssa_solve` keeps
`Result<u64, SsaError>`; `SsaError` gains `Inconsistent` and `Unsupported`; `NotAnomalous` only on
an exact count.

### C-GhsTransferGuard

Compiler + test-enforced. *To be frozen at 8.* Consumed by 12.
`transfer_point -> Result<_, GhsError>` with `NotOnCurve` and `NotInSubfield`.

### C-RhoFactorOutcome

Compiler + test-enforced. *To be frozen at 9.* Consumed by 12.

```rust
pub fn factor(n: u128, max_c: u64, batch_size: usize) -> Result<Factorization, FactorError>;
pub struct Factorization { pub d: u128, pub evidence: FactorEvidence } // 1 < d < n, n % d == 0
pub enum FactorEvidence { Even, Rho { c: u64 } }
pub enum FactorError { InvalidInput, Prime, Unsupported, BudgetExhausted { constants: u64 } }
```

### C-EcmOutcome

Compiler + test-enforced. *To be frozen at 10.* Consumed by 12.
`ecm_factor -> Result<Uint<4>, EcmError>` with `InvalidInput`, `Prime`, `CurvesExhausted { curves }`.

### C-ShorFactorOutcome

Compiler + test-enforced. *To be frozen at 11.* Consumed by 12.

```rust
pub fn factor(n: u64, seed: u64) -> Result<ShorFactorization, ShorFactorError>;
pub struct ShorFactorization { pub p: u64, pub q: u64, pub evidence: ShorEvidence } // p·q == n
pub enum ShorEvidence { Even, PrimePower { prime: u64, exponent: u32 }, CommonFactor { a: u64 },
                        OrderFinding { a: u64, seed: u64, s: u64, r: u64 } }
pub enum ShorFactorError { NotComposite, Unsupported { qubits: usize },
                           AttemptsExhausted { attempts: usize } }
```

Prose invariant: `order_from_phase` returns the multiplicative order, not a multiple.

## Progress ledger

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Define a verified rho discrete-log outcome and adopt it in Brent rho | pending | | |
| 2 | Bound the distinguished-point rho solvers and verify their results | pending | | |
| 3 | Recover the log from negated distinguished-point collisions | pending | | |
| 4 | Verify Pohlig–Hellman results and name the small-subgroup method | pending | | |
| 5 | Verify binary-field rho and Koblitz results | pending | | |
| 6 | Verify MOV results and type its failure classes | pending | | |
| 7 | Separate SSA failure classes and handle zero logs and large primes | pending | | |
| 8 | Reject off-curve and non-subfield points in the GHS transfer | pending | | |
| 9 | Type rho factoring outcomes and bound each attempt | pending | | |
| 10 | Type ECM factoring outcomes | pending | | |
| 11 | Report Shor factoring provenance and return the minimal order | pending | | |
| 12 | Correct the solver-outcome exposition ◆ | pending | | |
| 13 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

*(none yet)*

## Discoveries & risks

- **Inferred (sharding), enrolled 3:** plain `solve_dp` keys distinguished points by the low x
  word only, so a P/−P collision may return a wrong k. If the walk provably never meets −P:
  internal-continue (keep the synthetic KAT; record why).
- **Inferred, enrolled 4:** the Pohlig–Hellman/Silver order-p step is a table over the subgroup.
  If the primary source does not support this: **destructive-HALT** — the user's choice to keep
  the enumeration rests on it; the alternatives (rho on tiny groups, baby-step giant-step) are a
  scope decision.
- **Inferred, enrolled 5:** a wrong τ eigenvalue would give a silently wrong Koblitz k. If
  verification fails on a fixture: internal-continue fix within `binary_ecdlp`; record.
- **Inferred, enrolled 7:** the SSA k_raw/2 correction's mathematical source is unstated. If it
  masks a lifting defect: additive-reshard (surface, do not absorb).
- **Confirmed:** Pohlig–Hellman's end-to-end tests never reach rho (all primes ≤ 5); Shor
  factoring's `order_from_phase` returns multiples (its own test expects 8 for order 4); no rho
  solver verifies k·G == Q internally except SSA.
- **Internal-continue:** call-site-only adaptation in a file outside a row's expected files when
  the compiler forces it; budget defaults and seed re-selection by a bounded recorded search;
  replacing the wall-clock GLV test with a step comparison.
- **Additive-reshard (surface, do not absorb):** MOV's k=2 log engine and NFS-DL's `solve_dl`;
  Shor factoring register width; deduplicating `iqft_on_qubits`; principle-N label survivors
  outside a row's files; any new crate edge.
- **Destructive-HALT:** verification shows a solver returning wrong answers on existing fixtures
  in a way not fixable inside its module; or a contract needs a crate edge to hold.
- **Preflight:** PLAN and ROADMAP will be dirty after sharding; commit them before `/plan-run`.

## Notes for executors

**Sizing.**
- (1) Ambient complexity is moderate-high: `rho/src/ecdlp/mod.rs` is 1432 lines with threaded
  walkers; the other modules are small.
- (2) The floor is set by row 1's signature change (all callers at once) and row 2's shared
  coordinator; neither is a split candidate.
- (3) Cost of a design error is high for row 1 (consumed by rows 2–5), which put an inflection
  juncture on it; moderate elsewhere.
- (4) Correctness-criticality is high: the arc exists because results go unchecked. This drove
  the 2|3 split.
- (5) Inner-loop bandwidth is high for rho and numth (fast tests), low for Shor factoring
  (circuit runs).

Outlier check: rows 6, 7, 8 and 10 are below the usual band; they stay separate because they are
independent modules (no catch-all row). Row 1's eight files are one atomic signature change, mostly
call-site edits. No row uses F/C: the rows are small and tightly steered.

**Routing.** T1 implements every row; T0 runs the row-1 inflection juncture and the ◆ boundary
fork.

**Outcome contract.** `docs/development.md`, "Solver outcomes and fallbacks", is the authority;
`shor/src/ecdlp/mod.rs` is the worked reference shape.

**Register.** State the property or invariant, never the plan coordinate. Exposition goes in the
sections row 12 names; no new chapters.

**Anneal denylist** (grep enrolled durable files at ◆, excluding PLAN/ROADMAP/NOTES/SURVEY/AGENTS):

    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat [ABCI]\b
    \bC-[A-Z][A-Za-z]+\b
    \b[A-H]\.[A-Z0-9]+\b|[Pp]rinciple-?[0-9]+
    \bTrack [A-Z]\b|\bPhase [0-9]+\b
    (PLAN|ROADMAP|NOTES|SURVEY)(\.md)?\b|plan-(run|shard|juncture)

KAT, algorithmic step numbers, group notation (S₃, S4), bibliography keys ([S09]) and quantum
phase terminology are domain terms: classify collisions, do not rename.

**First run:** after committing the planning docs, from `@plan-admin`:

    /plan-run docs/PLAN.md halt-at-boundaries
