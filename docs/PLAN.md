# rDLP — Dedekind index criterion repair

Execution graph for the ROADMAP scope "Dedekind index criterion". The previous (completed) GNFS
square-root plan is preserved in git at `dd3aca4:docs/PLAN.md`.

## Purpose (design intent)

Make the Dedekind machinery in `shared/numfield` mathematically correct and honest about its scope.
`dedekind_criterion` returns wrong answers in both directions (x²−12 at 2 says p ∤ index; x²−3 at 3
says p | index) because it computes F from the reduced polynomial instead of the original integer f,
takes the radical with an integer gcd, and mishandles inseparable reductions (f' ≡ 0 mod p). The
public factorisation result returns `(p, α)` for inert primes, which is the unit ideal, not (p).

Scope boundary (DLP topicality): which degree-one primes (p, α − r) lie above p, with
multiplicity, and whether that bookkeeping is complete at p, is core GNFS/NFS-DL algebra and is in
scope. The full O_K factorisation at index-dividing primes (Round 2, Montes, distinct-degree
factorisation) is a neighbouring topic: mention it, do not implement or expound it. Do not wire
new consumers (the GNFS factor base, NFS-DL); that adoption belongs to the NFS-DL design scope.
Correct directly affected claims together with the code.

## Verify gate

- VERIFY_TYPES: `cargo check --workspace`
- VERIFY_TEST: `cargo test --workspace`
- Also mandatory: `cargo dep-policy`, `cargo audit`, `cargo coverage` (80% line floor).
- Not gates: `cargo fmt-check`, `cargo lint` (unadopted; do not fix drive-by).
- Opt-in oracle: `cargo test -p shared-numfield --test dedekind_pari_oracle_kat -- --ignored`
  (row 4; requires `gp` on PATH; must fail, not skip, when gp is missing).

## Session list

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 | Square-free decomposition over prime fields | A | Sonnet | C-FpPoly (existing) | `shared/numfield/src/fp_poly.rs`, `shared/numfield/tests/fp_poly_kat.rs`, `shared/numfield/docs/PEDAGOGY.md` |
| 2 | Canonical Dedekind index criterion | B | Sonnet | C-SquarefreeDecomposition | `shared/numfield/src/dedekind.rs`, `shared/numfield/tests/dedekind_kat.rs`, `shared/numfield/src/lib.rs`, `shared/numfield/docs/PEDAGOGY.md`, `docs/MATHEMATICS.md` |
| 3 | Honest degree-one prime report above p | B | Sonnet | C-SquarefreeDecomposition, C-IndexTest | `shared/numfield/src/dedekind.rs`, `shared/numfield/tests/dedekind_kat.rs`, `shared/numfield/src/lib.rs`, `shared/numfield/docs/PEDAGOGY.md` |
| 4 ◆ | Cross-check the index criterion and prime report against PARI/GP | I | Sonnet (◆ fork: Opus) | C-IndexTest, C-DedekindReport | `shared/numfield/tests/dedekind_pari_oracle_kat.rs` (new), `docs/development.md`, `shared/numfield/docs/PEDAGOGY.md` |
| 5 | Anneal cleanup (skippable) | I | — | none | durable files touched by rows 1–4 |

No `@architect` inflection rows: interfaces are frozen concretely below and no production code
consumes the changed API (only `is_bad_prime`/`discriminant`, which stay unchanged).

## Session detail

### 1 — Square-free decomposition over prime fields

Deliverable: `FpPoly::squarefree_decomposition` and `FpPoly::radical` per C-SquarefreeDecomposition,
including the characteristic-p case (f' = 0 ⇒ f = g(x^p) = g(x)^p, since Frobenius fixes F_p, so
the p-th root takes the coefficients at indices divisible by p). Short PEDAGOGY §12 subsection.

KATs (independently derived by hand): x² mod 2 → [(x, 2)]; x³ + 1 mod 3 → [(x+1, 3)];
x²(x+1)³(x+2)⁶ mod 3 → [(x,2),(x+1,3),(x+2,6)]; x(x+1)⁵ mod 5; a squarefree polynomial →
[(monic f, 1)]; nonzero constant → []; zero → `ZeroPolynomial`; (x+1)² over p = 2⁶¹−1 (no p-th
root path). Property check: lc · ∏ sᵢ^eᵢ = f, each sᵢ squarefree, pairwise gcd 1.

Subtleties: multiplicities can be multiples of p combined with non-multiples (the x+2 factor above);
a Yun-only implementation fails these. Deferrals: none.

### 2 — Canonical Dedekind index criterion

Deliverable: public `dedekind_index_test` per C-IndexTest, implementing the canonical criterion
(Cohen, *A Course in Computational Algebraic Number Theory*, §6.1, Dedekind's criterion): reduce
f̄ = ∏ tᵢ^eᵢ; ḡ = radical(f̄), h̄ = f̄/ḡ; take any monic integer lifts g, h; F = (f − g·h)/p
computed from the **original integer f**; T = gcd(F̄, ḡ, h̄) monic; p | [O_K : ℤ[α]] iff deg T > 0.
`dedekind_factor_extended` routes through it. Delete the old helpers it replaces
(`poly_exact_div_mod`, integer-gcd use). Correct rustdoc, numfield PEDAGOGY §9 (steps, the ℤ[√2]
example's reasoning, and the Round 2 passage's false "can be excluded from the factor base" claim;
keep Round 2 at mention level) and MATHEMATICS ≈387–392 (the code's "bad prime" means p | disc(f);
the prose means p | index; make the two terms distinct and consistent).

KATs: execute the two static traces first and record what the old code returned. Then: x²−12 at 2
→ divides, T = x; x²−3 at 3 → does not; x²−2 at 2 → does not; x²+3 at 2 → divides, T = x+1;
x³−10 at 3 → divides, T = x+2 (10 ≡ 1 mod 9); x³−2 at 3 and at 2 → does not (ℤ[∛2] is maximal);
a good prime → T = 1. Lift independence: shifting lift coefficients by p leaves T unchanged.
`NotMonic` and `ConstantPolynomial` errors.

Subtleties: confirm Cohen's statement that deg T = log_p [O′ : ℤ[α]] for the one-step enlargement
O′ = ℤ[α] + (U(α)/p)ℤ[α]. If it can't be confirmed from the reference, drop
`enlargement_exponent` and record the change in the digest (additive contract narrowing, internal).
Leave `is_bad_prime` and `discriminant` semantics unchanged (GNFS consumes them). Deferrals: the
report type (row 3).

### 3 — Honest degree-one prime report above p

Deliverable: replace `dedekind_factor`, `dedekind_factor_extended` and `DedekindResult` with a
single `dedekind_factor` returning `DedekindReport` per C-DedekindReport. No sentinel: an inert
prime has no degree-one primes and `nonlinear_degree = deg f`. Roots come from `FpPoly::roots` on
the radical under caller limits, so exhaustion is `RootSearchExhausted`, never a claim;
multiplicities come from the square-free decomposition. Correct numfield PEDAGOGY §7 ("what the
two-element form can and cannot represent"), §8 (remove the sentinel section), the §10 Dedekind
consumer passage (the factor base does its own root search and uses only `is_bad_prime`; state
that, without plan labels), the §11 KAT table, and the `lib.rs` item list.

KATs: x²−2 at 3 → no degree-one primes, nonlinear_degree 2 (regression for the unit-ideal
sentinel); at 7 → r ∈ {3, 4}, multiplicity 1; at 2 → r = 0, multiplicity 2, index test negative;
x³−x−1 at 5 → r = 2, nonlinear_degree 2; x³−10 at 3 → r = 1, multiplicity 3, index test positive;
truncated budget → `RootSearchExhausted`; invariant Σ multiplicity + nonlinear_degree = deg f for
every fixture.

Subtleties: when the index test is positive, the report claims only "prime ideals of ℤ[α] of norm
p" and multiplicities in f̄; never O_K ramification. Deferrals: factor-base adoption (NFS-DL
design), non-monic support.

### 4 ◆ — Cross-check the index criterion and prime report against PARI/GP

Deliverable: `#[ignore]`-gated `kat_pari_dedekind_*` tests. PARI computes v_p(index) =
(valuation(poldisc(f), p) − valuation(nfdisc(f), p))/2. Assert `divides_index ⇔ v_p > 0` and (if
kept) `enlargement_exponent ≤ v_p`. Where p ∤ index, `idealprimedec(nfinit(f), p)` must agree with
the report: degree-one primes matched by `idealval(nf, x − r, 𝔭) > 0` with e = multiplicity, and
Σ e·f over f > 1 = `nonlinear_degree`. Fixtures are chosen for coverage, not a fixed count:
inseparable reductions, divisible and non-divisible cases, ramified-but-maximal cases, p ∈ {2, 3,
5, 7}, degrees 2–4, including the generated GNFS cubic fixture's polynomial where monic. Check
`polisirreducible` in gp and fail loudly otherwise. A missing gp, a gp error or a mismatch fails
the test. Add the targeted invocation to `docs/development.md`; update the §11 KAT table.

Integrative duty: sweep superseded claims (the sentinel, "T = 1 because PID", the factor-base
consumer and exclusion claims, the old criterion steps) across numfield rustdoc, numfield PEDAGOGY,
MATHEMATICS and gnfs PEDAGOGY's bad-prime passages. Correct contradictions in the enrolled files
and record others in the digest. Anneal the enrolled durable files against the denylist.

Subtleties: duplicate the small gp-invocation helper locally (cross-crate test sharing would need a
new crate); record it as a minimalism candidate rather than refactoring.

### 5 — Anneal cleanup (skippable)

`/plan-anneal` over the arc's durable files: translate leaked plan coordinates and rewrap ledgered
line-length overruns. The driver may mark it `skipped` if the ◆ anneal ran clean and no overruns
were ledgered.

## Cross-session contracts

### C-FpPoly — existing checked prime-field polynomials (input)

Compiler-enforced. Defined in `shared/numfield/src/fp_poly.rs` (prior arc). Consumed by 1–3.
Additive extension only; no existing signature changes.

### C-SquarefreeDecomposition

Compiler + test-enforced. Defined in 1. Consumed by 2, 3.

```rust
impl FpPoly {
    /// f = lc(f) · ∏ sᵢ^eᵢ with each sᵢ monic, positive-degree, squarefree, pairwise coprime;
    /// returned in strictly ascending eᵢ. A nonzero constant returns `[]`.
    pub fn squarefree_decomposition(&self) -> Result<Vec<(FpPoly, usize)>, FpPolyError>;
    /// Monic product of the distinct irreducible factors; `one` for a nonzero constant.
    pub fn radical(&self) -> Result<FpPoly, FpPolyError>;
}
// Errors: FpPolyError::ZeroPolynomial for the zero polynomial.
```

### C-IndexTest

Compiler + test-enforced. Defined in 2. Consumed by 3, 4.

```rust
pub fn dedekind_index_test(f: &IntPoly, p: PrimeModulus) -> Result<IndexTest, DedekindError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexTest { /* private: prime, obstruction */ }
impl IndexTest {
    pub fn prime(&self) -> PrimeModulus;
    /// T = gcd(F̄, ḡ, h̄), monic.
    pub fn obstruction(&self) -> &FpPoly;
    /// p | [O_K : ℤ[α]] iff deg T > 0.
    pub fn divides_index(&self) -> bool;
    /// deg T = log_p [O′ : ℤ[α]]; may be dropped in 2 if the reference does not confirm it.
    pub fn enlargement_exponent(&self) -> usize;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DedekindError {
    NotMonic,
    ConstantPolynomial,
    Field(FpPolyError),
    /// Used from 3 on: resource exhaustion, not a mathematical statement.
    RootSearchExhausted { distinct_root_count: usize, evaluations: u64 },
}
```

Prose invariant: f is monic (irreducibility is the caller's, as for `NumberField`); F is formed
from the original integer f; the result is independent of the lifts chosen.

### C-DedekindReport

Compiler + test + prose-enforced. Defined in 3. Consumed by 4.

```rust
pub fn dedekind_factor<'a>(field: &'a NumberField, p: PrimeModulus, limits: RootSearchLimits)
    -> Result<DedekindReport<'a>, DedekindError>;

pub struct DedekindReport<'a> {
    pub degree_one: Vec<DegreeOnePrime<'a>>, // ascending r
    pub nonlinear_degree: usize,
    pub is_bad_prime: bool,                  // p | disc(f), unchanged meaning
    pub index: IndexTest,
}
pub struct DegreeOnePrime<'a> { pub ideal: Ideal<'a>, pub multiplicity: usize }
```

Prose invariants: Σ multiplicity + nonlinear_degree = deg f. If `!index.divides_index()`, each
listed ideal is a prime of O_K with residue degree 1 and ramification index = multiplicity, and
nonlinear_degree = Σ eᵢfᵢ over the primes with fᵢ > 1. If it divides the index, each listed ideal is
a prime of ℤ[α] of norm p and nothing is claimed about pO_K. No sentinel ideals.

## Progress ledger

Bindings: PLAN = `docs/PLAN.md`; VERIFY_TYPES = `cargo check --workspace`; VERIFY_TEST = `cargo test
--workspace`; mandatory extras per row: `cargo dep-policy`, `cargo audit`, `cargo coverage`. Run
config: no flags (self-continue at junctures; no `may-reshard`, no `halt-at-boundaries`).

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Square-free decomposition over prime fields | done | adf138e | C-SquarefreeDecomposition |
| 2 | Canonical Dedekind index criterion | done | 0073ab8 | C-IndexTest (without `enlargement_exponent`) |
| 3 | Honest degree-one prime report above p | done | 17b1604 | C-DedekindReport (report types carry no derives) |
| 4 | Cross-check against PARI/GP ◆ | done | a2e0d39 | none (test-only); boundary check pending |
| 5 | Anneal cleanup (skippable) | pending | | |

## Action-frame digest

### 2 — 2026-11-08
Discovery/flex: `enlargement_exponent` dropped from C-IndexTest: Cohen §6.1 Theorem 6.1.4 could not
be read to confirm deg T = log_p [O′ : ℤ[α]] (ring-ness of O′), as the plan's drop clause allowed.
Affected: C-IndexTest (additive narrowing; internal). Row 4's `enlargement_exponent ≤ v_p` assertion
is "if kept" and therefore void.
Deferred: yes — if a later scope wants it back, confirm against the printed theorem; the module index
p^{deg T} was derived independently by the implementer but not the ring property.
Texture: old-code traces reproduced: x²−12@2 wrongly "not divisible", x²−3@3 wrongly "divisible".
`DedekindError` did not previously exist (created, not reconciled); `RootSearchExhausted` is declared
but unconstructed until row 3. `dedekind_factor_extended` now panics for non-prime p (documented).
Lift-independence is tested as a unit test in `dedekind.rs` (private lift entry point), not in
`dedekind_kat.rs`. Committer wrote the commit body unwrapped (off the 72-col convention).

### 3 — 2026-11-08
Discovery/flex: `DedekindReport`/`DegreeOnePrime` carry no `Debug`/`Clone`/`PartialEq`/`Eq`, because
`Ideal` derives nothing and `ideal.rs` is outside the row's expected files; tests use `.err()` and
`matches!`. Reducible monic fixtures are used for mixed-multiplicity cases (`NumberField::new` does
not check irreducibility).
Affected: C-DedekindReport (derive set only; compiler-visible signature unchanged).
Deferred: yes — row 4 must not assume `enlargement_exponent` (dropped in 2; the C-IndexTest text in
this PLAN still lists it) and cannot `assert_eq!` on reports; boundary adjudicator should also check
that durable prose no longer lists `enlargement_exponent`.
Texture: no callers of the removed functions outside `shared/numfield` (confirmed by the
implementer; `gnfs` factor base does its own root search and uses only `is_bad_prime`). Truncated
budget on x²−2 at 7: budgets 0–4 exhaust, 5 succeeds; a rootless reduction needs no budget.

### 4 — 2026-11-08
Discovery/flex: none to the contracts: PARI/GP agreed with the criterion and the report on all 21
monic fixtures × p ∈ {2, 3, 5, 7} (84 cases; 5 `#[ignore]` tests, ~4–8 s with gp 2.17.3). Includes
the generated GNFS cubic x³+15x²+29x+8. The sweep found the superseded claims fixed in numfield and
MATHEMATICS, but two surviving conflations of "bad prime" with "p divides the index" outside the
enrolled files: `gnfs/docs/PEDAGOGY.md:881-883` ("At bad primes, Dedekind's theorem does not apply
directly") and `gnfs/src/sieve/factor_base.rs:31-32` (same claim in a module comment). Bad prime can
be maximal (x²−2 at 2).
Affected: none (contracts); durable prose outside the row's expected files.
Deferred: yes — the two gnfs hits are outside the enrolled scope (edit outside `shared/numfield` is
"additive-reshard: surface, do not absorb"); the boundary adjudicator rules on them. Also
minimalism candidate: the gp-invocation helper is duplicated in the new test file and in several
gnfs/rho oracle tests. Pre-existing plan labels in numfield PEDAGOGY §3/§4/§10/§11 (C-NF, G.A.1a,
the §10 contract table, D.A, "sub-track") are outside the edited passages: project-wide anneal scope.
Texture: committer wrote one commit body unwrapped in row 2 (0073ab8); not amended.

## Discoveries & risks

- **Confirmed (ROADMAP), enrolled 2:** wrong-result criterion; static traces x²−12@2 and x²−3@3 not
  yet executed. If a trace does not reproduce as described: internal-continue (the canonical
  criterion is the contract, not the trace).
- **Inferred during sharding, enrolled 3:** the inert sentinel (p, α) is the unit ideal
  (p ∤ f(0) and f(0) ∈ (α)); the rustdoc and PEDAGOGY §8 claim it equals (p).
- **Confirmed, enrolled 3:** PEDAGOGY §10 claims the GNFS factor base consumes `dedekind_factor`;
  it does its own root search and uses only `is_bad_prime`.
- **Inferred, not enrolled:** for non-monic f, `discriminant` returns lc(f)·disc(f) (Res(f, f')
  without division by lc), so `is_bad_prime` does flag p | c. The `factor_base.rs` comment that
  p | c "can be invisible" describes a different flag. Do not change `is_bad_prime`; leave it to the
  project-wide claim/anneal scope.
- **Destructive-HALT:** PARI disagrees with the canonical criterion or report on any fixture after
  the gp derivation is re-checked (contract-level contradiction).
- **Additive-reshard (surface, do not absorb):** any edit outside `shared/numfield` or the listed
  docs; any consumer wiring; non-monic support; Round 2/distinct-degree factorisation (out of scope
  by the purpose boundary).
- **Preflight:** PLAN and ROADMAP will be dirty after sharding; commit them before `/plan-run`.

## Notes for executors

**Sizing.** (1) Ambient complexity is moderate: one crate, a well-tested `FpPoly` substrate.
(2) The floor is low: each row is one mathematical object. (3) Cost of a design error is low to
moderate, since the public API break has no production consumers. (4) Correctness-criticality is
high, which pushes toward the 1|2 split. (5) Inner-loop bandwidth is high: fast exact KATs. No row
uses F/C: the rows are small and tightly steered.

**Routing.** Sonnet implements every row; Opus runs only the ◆ boundary fork (intent check and
claim adjudication).

**Register.** State the property or invariant, never the plan coordinate. Mathematical exposition
belongs in the PEDAGOGY/MATHEMATICS sections named per row; no new chapters.

**Anneal denylist** (grep enrolled durable files at ◆, excluding PLAN/ROADMAP/NOTES/SURVEY):

    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat [ABCI]\b
    \bC-[A-Z][A-Za-z]+\b
    \bG\.[A-Z0-9]|\bD\.A\b|principle-[0-9]+
    (PLAN|ROADMAP|NOTES|SURVEY)(\.md)?\b|plan-(run|shard|juncture)

KAT, algorithmic step numbers and group notation (S4) are domain terms: classify the collision, do
not rename. Pre-existing §10 plan labels outside the edited passage are project-wide anneal scope.

**First run:** after committing the planning docs, from `@plan-admin`:

    /plan-run docs/PLAN.md halt-at-boundaries
