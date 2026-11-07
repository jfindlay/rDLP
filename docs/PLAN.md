# GNFS correctness — canonical reconstruction and verified composition

**G.1 design accepted; sharding approved on 2026-11-06.** The user accepted the mathematical design,
then approved fourteen implementation rows, their contract boundaries and tier routing, deterministic
u64 prime validation, and budgeted root enumeration. The existing `gnfs -> shared-padic` path edge,
typed public API correction, affine-collection boundary, and integration-discovery limits are part
of that acceptance. Execution is a separate `/plan-run` action; this session starts no implementation
workers, runtime checks, installs, or commits. All implementation rows remain pending.

The prior completed ALIGN plan remains at `3cd6fe5:docs/PLAN.md`; its adopted gates remain binding.
ROADMAP retains unrelated deferred work. The dated G.1 NOTES entry records its earlier design-only
authorization historically; this acceptance supersedes that restriction, not its technical content.

## Purpose (design intent)

Demonstrate factorization through generated sieve relations, filtering, computed dependencies, exact
square reconstruction, and checked gcd assembly. Neither a manufactured intermediate answer nor an
incidental gcd counts as evidence for the advertised pipeline.

Implement the established direct Newton method first, with the published Thomé CRT method covering
fields without inert primes. Build reusable finite-field arithmetic rather than choosing the main
algorithm to fit currently missing substrate. This is not a new-algorithm campaign. Do not copy
heuristic precision estimates or unchecked inversions from production reference implementations.

Sources, exact derivations, independently checkable fixtures, and source-level findings are in
NOTES under **GNFS mathematical contract and canonical reconstruction design (2026-11-06, G.1)**.
That section is required handoff context until its contracts land in tests and focused human docs.

## Verify gate

Bindings are confirmed from `AGENTS.md`, `.cargo/config.toml`, and `docs/development.md`, not from
fresh execution. Before the first source change, re-establish every adopted baseline command.
Every delivery row must finish with all five gates green:

- `VERIFY_TYPES`: `cargo check --workspace`
- `VERIFY_TEST`: `cargo test --workspace`
- Additional mandatory gates: `cargo dep-policy`, `cargo audit`, `cargo coverage` (80% line floor).

`cargo fmt-check` and `cargo lint` remain declared but unadopted workspace gates. No drive-by
backlog cleanup, coverage exclusions, lowered thresholds, weakened mathematical assertions, or
benchmark substitution. Correct wrong expected values with independent mathematical justification.
External PARI/CADO checks remain opt-in and ignored; they validate results, never compute the
implementation's answer. Existing oracle placeholders are not evidence. Check tool availability
without authorizing an automatic installation; stop for a missing tool or failing adopted baseline.

## Session list

All workers are T1/Sonnet `@build`. `@architect` marks a T0/Opus inflection-design fork **before**
implementation, not an Opus implementation worker. `◆` marks delivery-group boundaries, not proof
that the entire campaign is repaired. Rows 1–7 deliver reconstruction; rows 8–14 deliver checked
composition. Each row freezes the corresponding contract listed below. File braces denote exactly
the named paths, not permission to edit an entire directory. New paths are identified in the detail.

| # | Session | Cat | Tier | Consumes | Expected files |
|---|---------|-----|------|----------|----------------|
| 1 | Add checked finite-field polynomial arithmetic — @architect | A | T1 + T0 design | C-Design | `shared/numfield/src/{fp_poly,lib}.rs`, `shared/numfield/tests/fp_poly_kat.rs`, `shared/numfield/docs/PEDAGOGY.md` |
| 2 | Add square roots in irreducible residue fields | B | T1 | C-Design, C-FpPoly | `shared/numfield/src/{residue_field,lib}.rs`, `shared/numfield/tests/residue_field_kat.rs`, `shared/numfield/docs/PEDAGOGY.md` |
| 3 | Define checked square-root contexts and outcomes — @architect | A | T1 + T0 design | C-Design | `gnfs/src/sqrt/{context,mod}.rs`, `gnfs/tests/sqrt_context_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 4 | Construct normalized targets and certified bounds | B | T1 | C-Design, C-SqrtContext | `gnfs/src/sqrt/context/normalize.rs`, `gnfs/src/sqrt/{context,mod}.rs`, `gnfs/tests/sqrt_normalize_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 5 | Reconstruct algebraic roots by p-adic Newton lifting | B | T1 | C-Design, C-FpPoly, C-ResidueField, C-SqrtContext, C-NormalizedTarget | `gnfs/src/sqrt/{padic,mod}.rs`, `gnfs/tests/sqrt_padic_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 6 | Construct lifted split-prime root contributions — @architect | A | T1 + T0 design | C-Design, C-FpPoly, C-ResidueField, C-SqrtContext, C-NormalizedTarget | `gnfs/src/sqrt/{crt,mod}.rs`, `gnfs/tests/sqrt_crt_kat.rs`, `gnfs/Cargo.toml`, `Cargo.lock`, `gnfs/docs/PEDAGOGY.md` |
| 7 | Recombine CRT signs and dispatch verified roots ◆ | I | T1; T0 boundary | C-Design, C-SqrtContext, C-NormalizedTarget, C-DirectRoot, C-SplitContributions | `gnfs/src/sqrt/{crt,mod,context}.rs`, `gnfs/tests/sqrt_crt_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 8 | Validate relations and classify algebraic parity — @architect | A | T1 + T0 design | C-Design, C-SqrtContext | `gnfs/src/sieve/{checked,mod,factor_base}.rs`, `gnfs/src/filter/mod.rs`, `gnfs/tests/checked_relations_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 9 | Preserve kernels through elimination and XOR lifts | B | T1 | C-Design, C-CheckedRelations | `gnfs/src/filter/{matrix,merge}.rs`, `gnfs/tests/merge_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 10 | Filter candidates with individual root characters | B | T1 | C-Design, C-FpPoly, C-CheckedRelations, C-XorLift | `gnfs/src/linalg/{qc,mod}.rs`, `gnfs/tests/linalg_substrate_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 11 | Validate dependencies against original relations — @architect | A | T1 + T0 design | C-Design, C-SqrtContext, C-NormalizedTarget, C-CheckedRelations, C-XorLift, C-RootCharacters | `gnfs/src/sqrt/{dependency,mod}.rs`, `gnfs/tests/sqrt_dependency_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 12 | Enumerate bounded combinations of kernel vectors | B | T1 | C-Design, C-Dependency | `gnfs/src/linalg/{kernel,mod}.rs`, `gnfs/tests/kernel_combinations_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 13 | Migrate factor assembly to verified root reports | I | T1 | C-Design, C-SqrtContext, C-NormalizedTarget, C-VerifiedRoot, C-CheckedRelations, C-Dependency, C-KernelSearch | `gnfs/src/sqrt/{algebraic,rational,assembly,mod}.rs`, `gnfs/src/lib.rs`, `gnfs/tests/{sqrt_algebraic_kat,sqrt_rational_kat,factor_end_to_end_kat}.rs`, `gnfs/docs/PEDAGOGY.md`, `docs/{MATHEMATICS,BENCHMARKS}.md` |
| 14 | Demonstrate factorization from generated relations ◆ | I | T1; T0 boundary | C-Design, C-VerifiedRoot, C-CheckedRelations, C-XorLift, C-RootCharacters, C-Dependency, C-KernelSearch, C-FactorReport | `gnfs/tests/factor_pipeline_kat.rs`, `README.md`, `gnfs/docs/PEDAGOGY.md`, `docs/MATHEMATICS.md` |

## Session detail

Each row is one green, commit-shaped unit, including substantive KATs and directly affected
explanatory text. No F/C pairs are selected. Module names may be refined at their interface juncture
with an explicit fileset update; an executor must not silently widen its fileset. Documentation
edits stay within the row's topic, not entire chapters. The public legacy path is not declared
repaired before rows 13–14 and the final independent review.

### 1 — Finite-field polynomials (G.2a)

New: `shared/numfield/src/fp_poly.rs`, `shared/numfield/tests/fp_poly_kat.rs`.
Deliver: polynomial arithmetic, division/gcd, powmod, and a canonical irreducibility test over a
validated prime modulus. No full factorization engine. Accept independent reducible/irreducible
fixtures, repeated factors, exact degree handling, and overflow-safe arithmetic.
Use a deterministic u64 primality algorithm with a documented valid witness/range argument, not
probable primality or a caller-trusted assertion. Cover composite/pseudoprime inputs. No new
numfield dependency is approved. Freeze canonical coefficient representation, checked modulus and
error types, arithmetic signatures, and resource-accounted exhaustive polynomial-root enumeration
for the later split-prime and character searches. Machine-prime range exhaustion is a resource
outcome; arbitrary-precision field sizes and lifting moduli are not restricted to u64.

### 2 — Residue-field square roots (G.2b)

New: `shared/numfield/src/residue_field.rs`, `shared/numfield/tests/residue_field_kat.rs`.
Deliver: checked residue-field context, inverses, quadratic-residue test, and Tonelli–Shanks with
bounded nonresidue search. Cardinalities/exponents use arbitrary precision. Accept degree >1 roots,
nonresidues, zero handling, invalid moduli, and explicit resource outcomes.

### 3 — Checked context and result types (G.2c, interface half)

New: `gnfs/src/sqrt/context.rs`, `gnfs/tests/sqrt_context_kat.rs`.
Deliver: checked polynomial/map context and owned root-witness/result types with non-forgeable success
construction, failure classes, and resource options. Validate N>1, g=x-m, f(m)=0 mod N, nonzero
leading coefficient coprime to N, and squarefree F. Define the input boundary for rejecting nonunit
relation factors before reconstruction. Accept invalid roots/moduli, repeated factors, and
incidental leading-coefficient gcd provenance. Do not call panic-prone legacy validation on invalid
N, alter `PolyPair::number_field`, or introduce placeholder reconstruction bodies.
Checkpoint: freeze concrete Rust signatures, witness ownership, resource accounting, and the planned
normalization/verifier interface before rows 4–7. Final checked dependency signatures freeze at 11.

### 4 — Normalization, bounds, and exact verification (G.2c, arithmetic half)

New: `gnfs/src/sqrt/context/normalize.rs` (child module of the context module, per the row-3
non-forgeability design; fileset refined at the row-3 inflection, 2026-11-07),
`gnfs/tests/sqrt_normalize_kat.rs`. `gnfs/src/sqrt/context.rs` is edited only to declare the child
module and drop the success seal's dead-code allowance.
Deliver: gamma/eta/H construction, certified integer coefficient bounds, exact z and beta identity
checks, verified derivative inverse, and normalized-denominator modular mapping. Reject zero/unit
violations before prime search. Accept the NOTES non-monic odd-cardinality and outside-order
fixtures, zero divisors, the wrong-norm-map counterexample, and noninvertible denominators with typed
outcomes. Use exact floor/ceiling integer square roots for bounds, not an exact-square-only helper.
This unit makes the verifier usable before either reconstruction algorithm exists.

### 5 — Canonical direct p-adic reconstruction (G.2d)

New: `gnfs/src/sqrt/padic.rs`, `gnfs/tests/sqrt_padic_kat.rs`.
Deliver: bounded inert-prime selection, residue-field inverse-square-root seed, Newton doubling in
the quotient ring, centered reconstruction, exact verification. Accept nonconstant roots, even
degree, no real embedding, bound-boundary tests, local nonsquare and global nonsquare certificates,
and prime/precision resource limits. No real-sign heuristic or fixed CRT-prime count.
Quotient-ring arithmetic modulo prime powers is new and remains private to this implementation
unless the interface checkpoint justifies sharing it. Existing scalar `Zp` is not this substrate.

### 6 — Split-prime local lifting (G.2e, local contribution half)

New: `gnfs/src/sqrt/crt.rs`, `gnfs/tests/sqrt_crt_kat.rs`.
Deliver: bounded split-prime/root selection, existing scalar `shared_padic::hensel_lift` reuse,
prime-power interpolation, CRT weights, and signed local contribution records without selecting
global signs. Manifest/lockfile changes are only for the existing `gnfs -> shared-padic` path edge.
Accept independently checked lifted F roots, local square identities, Lagrange basis identities,
CRT idempotents, and local nonresidue/resource outcomes. Moduli/interpolation use arbitrary
precision, not the old u64-only helpers. Freeze the contribution representation and bound-aware
matching interface at the inflection checkpoint, with list/match accounting for row 7.

### 7 — Exact sign recombination and method dispatch (G.2e, global half) ◆

Fileset refined at the row-6 inflection (2026-11-07): `gnfs/src/sqrt/context.rs` is added, only to
drop the stale dead-code allowance on the nonsquare-certificate constructor and to correct the CRT
limits doc comment ("each split prime", not "each prime power"). No other widening.

Deliver: exact meet-in-the-middle sign matching, complete verification, and bounded
primary-to-fallback dispatch, consuming signed local contributions rather than local candidate roots
misrepresented as global witnesses.
Support ell=1 and a small multi-prime case; never use arbitrary cross-prime sign normalization.
Accept the V4 quartic in NOTES, local false positives, cross-method agreement, and list/match budget
exhaustion. Do not silently return NotSquare after incomplete sign search. No new external package.
Handle modular interval wraparound exactly. Check all coefficient bounds and z²=H after each match.
This boundary reviews reconstruction fidelity and anneals its durable text; it does not certify the
still-unmigrated public factor path. The V4 fixture must exercise a genuinely no-inert-prime field.

### 8 — Validated relations and honest algebraic parity (G.3a)

Fileset refined at the row-8 inflection (2026-11-07): one new file, `gnfs/src/sieve/checked.rs`, holds
the checked relation types so the legacy collector in `sieve/mod.rs` stays untouched. No manifest
or `lib.rs` change. Doc claims about `Relation` that are false (sign through a real embedding;
characters guaranteeing a square) are corrected within this row's files and topic.

New: `gnfs/tests/checked_relations_kat.rs`.
Deliver: checked indices/exponents, exact norm/sign reconstruction, good affine root labels, and
explicit exceptional-column classification in the GNFS builder. Keep the full Relation record and
existing NFS-DL behavior. Do not impose raw homogeneous-norm-square tests on non-monic candidates.
Accept ordinary and exceptional prime cases with independent arithmetic; document the inherited
projective-collection limitation. Expansion to complete bad/projective ideals requires review.
Freeze a separate checked relation context and matrix layout; do not add fields to public Relation
or SparseMatrix just to carry validation state. Validation and column classification stay together
because the contract is that the new builder accepts only the validated record interpretation.

### 9 — Correct elimination and overlapping provenance (G.3c)

Deliver: pivot-relative weight-three elimination and XOR provenance. Accept exhaustive small-matrix
kernel/lift comparisons, the three-identical-rows counterexample, overlapping provenance
cancellation, and column-weight consistency. Do not mistake arbitrary row pruning for lossless
elimination. Consume the matrix layout frozen at 8; test disjoint and overlapping provenance.

### 10 — Individual quadratic characters (G.3b)

Deliver: (q,s) characters of a-bs, bounded selection of valid pairs, explicit zero handling, and
correct XOR through provenance. Character count counts pairs, not rational primes. Norm characters
are removed from this path. Accept the norm-one nonsquare unit in NOTES and merged-row consistency.
Characters may be populated after filtering only if the complete provenance identity is preserved.
Select valid q not dividing c*disc(F), with f(s)=0 mod q and nonzero evaluations on all retained
original relations. Keep the legacy `select_qc_primes` helper compiling while old algebraic callers
still need it; it is not used as the new character contract. No permanent unsafe wrapper is added.

### 11 — Checked original-row dependency context (G.3d, validation half)

New: `gnfs/src/sqrt/dependency.rs`, `gnfs/tests/sqrt_dependency_kat.rs`.
Deliver: shared validated original-row context, checked matrix shape/indices/provenance, per-row
lifted parity, kernel membership, relation-pair uniqueness, and exact P/gamma products. Duplicate
relations must be reconciled before matrix construction, never by silently renumbering provenance.
Accept malformed indices/provenance, overlap cancellation, invalid kernel membership, exact product
fixtures, and non-monic products without the false raw homogeneous-norm-square precondition.
Freeze constructor signatures, matrix-layout binding, owned validated candidate output, and the
public factor-attempt/report interface that row 13 will implement. Old public APIs remain unchanged
in this additive row. Fake exponent matrices belong only in isolated unit tests.

### 12 — Bounded computed-kernel combinations (G.3d, enumeration half)

Fileset refined at the row-11 inflection (2026-11-07): `gnfs/src/linalg/mod.rs` is added, only to
re-export the new public kernel-search items (the `kernel` module is private). Row 12 must emit a
`Vec<KernelVector>` plus a `KernelSearchStatus` carrying a `Complete` variant and `is_complete()`;
the dependency validation layer never trusts these and never calls `KernelVector::verify` or
`expand_provenance` (the latter panics on out-of-range rows).

New: `gnfs/tests/kernel_combinations_kat.rs`.
Deliver: deterministic budgeted enumeration of distinct nonzero XOR combinations of computed kernel
vectors, with explicit complete-versus-resource-limited exhaustion. Freeze ordering, deduplication,
work/output limits, and resumability semantics if any; do not implement an unbounded 2^r scan behind
an output-only cap. KATs cover combinations beyond basis vectors, dependent/duplicate basis inputs,
zero cancellation, deterministic order, and budget limits. Synthetic vectors are isolated tests of
enumeration, never manufactured evidence for the checked factor path.

### 13 — Checked dependency-to-factor API migration (G.3d, assembly half)

Deliver: fallible rational/algebraic roots, verified reconstruction/map, checked modular congruences,
typed factor-attempt reports, and retry over row 12's bounded candidates. Replace the old norm map,
real-sign heuristic, and fixed-prime routine; migrate every caller and re-export in one atomic commit.
Accept malformed-input errors, nonsquare-to-success retry, successful combinations beyond basis
vectors, trivial gcds, separate incidental-factor provenance, and the nonconstant 527 seam fixture.
Fix the invalid f(m) fixture using independently justified input, not weakened assertions. Ignored
oracle callers must compile and remain ignored. Do not retain a panic or unverified-root wrapper.
The fileset is intentionally larger than the default: source signatures and all callers must move
together to keep a green tree. Focused claims corrections belong here, including removal of the
misleading Couveignes name for the replaced routine. Unrelated benchmark-tooling defects are not
enrolled. Implementation reasoning remains coupled, so this is not an F/C compliance-tail split.

### 14 — Generated sieve-to-factor demonstration (G.4) ◆

New: `gnfs/tests/factor_pipeline_kat.rs`; keep test helpers local unless an explicit fileset change
is approved. Update only pipeline-evidence claims in the named documentation files.
Compose actual `line_sieve` -> validated/deduplicated relations -> matrix -> singleton/pruning/merge
-> root characters -> existing kernel solver -> checked lifted candidates -> square roots -> gcd.
Test-only composition is sufficient; do not add a new public orchestration framework by default.

Predeclared seeds:

1. N=45113, m=31, f=x³+15x²+29x+8, from the published worked example in
   <https://www.cs.umd.edu/~gasarch/TOPICS/factoring/NFSmadeeasy.pdf>. Adapt sign conventions, never
   import its relation list or dependency. Check irreducibility and f(m)=N independently.
2. N=10403, degree 3, using the existing `select_base_m` generator. Require valid non-monic output;
   an unexpected domain result is recorded, not silently replaced with a convenient polynomial.

For each seed, try these profiles, alternating seeds at each profile:

| Profile | B_rat | B_alg | A | b_max |
|---------|-------|-------|---|-------|
| 1 | 30 | 90 | 400 | 41 |
| 2 | 60 | 180 | 600 | 82 |
| 3 | 120 | 360 | 800 | 128 |

Use a in [-A,A], b in [1,b_max], eight valid root-character pairs, and deterministic solver seeds.
Before running discovery, pin all remaining options in its harness: primary prime search, fallback
prime count/search, precision limits, and recombination limits. These must instantiate, not weaken,
the reviewed contract. Changing those options after observing results counts as a new configuration.

Limits: at most six configurations, 10^6 lattice points total, 4096 distinct nonzero computed-kernel
combinations per configuration, and 30 minutes aggregate discovery execution (excluding compilation
and adopted workspace checks). Stop at the first qualifying success or the first exhausted limit.
Record attempted/unattempted configurations, yields, validation failures, ranks/dependencies,
reconstruction outcomes, and factor provenance. No unrecorded parameter/fixture search.

Closure requires one deterministic runtime-generated success with irreducible degree >1 and nonzero
b relations materially participating. Assert relation validity, lifted parity, exact square
identities, correct modular images, congruence, and a proper factor. Pinned parameters are allowed;
injected relations, handwritten kernels, brute-force factor fallbacks, and oracle-supplied answers
are not. Core KATs separately cover non-monic and no-inert-prime cases; a single pipeline fixture
does not establish all-input success. Budget exhaustion is an obstruction report, not closure.

Discovery and its qualifying fixture remain one unit: do not reset budgets between sessions, gate
retries, or resumed workers. Record the winning configuration as a deterministic ordinary KAT; its
subsequent verification runs are not further parameter search. At closure, anneal the enrolled
durable surface and state only the supported scope established by the tests.

### Final independent correctness and claim review (G.5; boundary, not a delivery row)

Owner: an independent review context at the row-14 T0 boundary, after success and all adopted gates.
Do not dispatch this as a no-code commit row or fabricate a KAT to satisfy a row-shaped review.
The driver supplies the implemented contracts, canonical source pointers, KAT evidence, discovery
record, and claims diff to the independent reviewer; implementation self-report is not sufficient.
Review canonical-method fidelity, proof preconditions, fixture independence, and honest tested
scope.
Then review external-oracle readiness and whether to shard NFS-DL next. Do not carry forward claims
that finite QC guarantees a field square, that global sign needs a real root, or that the former
fixed-prime routine differed from production only by a scale knob.
Record the review verdict and unresolved findings in the action-frame digest. A green row-14
implementation alone is not campaign acceptance; unresolved correctness findings block closure.

## Cross-session contracts

Mathematical semantics in C-Design are accepted now. Concrete interfaces below are **to be frozen
at their producer row**, with inflection-design decisions written here before that row executes.
Each freeze records actual signatures/representations and named KATs for cold downstream consumers.
Producer dispatches receive both the contracts they consume and their own contract to establish.

### C-Design — accepted mathematical and scope obligations

**Enforcement:** prose now; compiler/test enforcement accumulates in producer rows.
**Defined-in:** accepted G.1 design and NOTES derivation, 2026-11-06. **Consumed-by:** 1–14 and reviews.

- **Normalization:** c=lc(f), theta=c alpha, F(t)=c^(d-1)f(t/c), h=ceil(|S|/2).
  Reconstruct z²=H where eta=c^(2h-|S|) product(ca-b theta) and H=F'(theta)² eta.
  Recover beta=z/(c^h F'(theta)) over Q and verify beta²=product(a-b alpha) exactly.
- **Map:** evaluate normalized rational beta coefficients at theta=cm, checking denominators
  modulo N. The returned y squares to P=product(a-bm) mod N, not to Norm(gamma). The rational
  root comes from an exact integer square P>0. Check x²=y² mod N before gcd extraction.
- **Domain:** N>1, g=x-m, f(m)=0 mod N, nonzero leading coefficient coprime to N, squarefree F,
  and nonzero/unit relation factors in Q[t]/F. Non-monic, even-degree, totally imaginary, and
  no-inert-prime fields are not excluded. Reducible squarefree quotients are not called fields;
  a genuine GNFS demonstration must certify irreducibility. Unsupported maps remain explicit.
- **Reconstruction:** primary = finite-extension Tonelli–Shanks plus inverse Newton lifting;
  fallback = Thomé split-prime lifting/CRT with exact meet-in-the-middle sign recombination.
  Both use the NOTES coefficient bound, modulus M>2B, and the same exact verifier/map.
- **Failures:** distinguish invalid input/invariants, rational nonsquare, algebraic nonsquare,
  unsupported map/domain, resource exhaustion, trivial gcd, and verified factor. Only a local
  obstruction or complete branch exhaustion beyond a proved bound certifies NotSquare. Search
  exhaustion never proves absence of inert/split primes. Incidental factors have separate
  provenance. Invalid invariants must not be silently retried as ordinary nonsquares.
- **Candidates:** rational sign/prime parity and good ideal/root-character constraints are filters,
  not proof that every kernel vector is an element square. Try bounded combinations of computed
  kernel vectors. Recompute lifted parity from validated original rows.
- **Collection:** preserve Relation's norm records and NFS-DL consumers. The GNFS checked builder
  omits unproved algebraic parity constraints at primes dividing c*disc(F). The existing primitive
  affine collector remains; complete projective/exceptional-ideal collection is not claimed.
- **Provenance:** every merged row is the XOR of its original rows. Both column sets and provenance
  use symmetric difference. Correct elimination preserves the restricted kernel; intentional row
  pruning is a separately reported loss of candidates.
- **Compatibility:** no global reinterpretation of `PolyPair::number_field`, no unsafe compatibility
  wrappers returning unverified roots. Freeze owned witness types, error/result shapes, resource
  options, and checked dependency-context signatures before consuming them downstream.
- **Scope:** no broad formatting/lint cleanup, unrelated NFS-DL fixes, new algorithm families,
  external runtime/oracle dependency, or unapproved dependency expansion. Durable text uses topic
  names, not these session labels.

### C-FpPoly — checked prime-polynomial arithmetic

**Enforcement:** compiler + test. **Defined-in:** 1. **Consumed-by:** 2, 5, 6, 10.
**Status:** **frozen at row 1 (commit a5e168b, 2026-11-06)**, designed at the row-1 inflection and
implemented with no signature or semantic deviation (21 named KATs plus two extras pass). Record any
later deviation from these signatures here; do not adopt it silently. Fileset unchanged: no `Cargo.toml` edit. Existing `num-bigint` (including
`BigUint`), `num-integer`, and `num-traits` suffice. Module `shared_numfield::fp_poly`; lib.rs adds
`pub mod fp_poly;` and re-exports every public item below from the crate root. No new dependency,
no full factorization engine, no Berlekamp, Cantor–Zassenhaus or randomized splitting, no panics.

Signatures (bodies elided; every public item gets rustdoc with `# Errors` where fallible):

    pub fn is_prime_u64(n: u64) -> bool;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct PrimeModulus { /* private: p: u64, proved prime at construction */ }
    impl PrimeModulus {
        pub fn new(p: u64) -> Result<Self, FpPolyError>;         // NotPrime
        pub fn get(self) -> u64;
        pub fn to_biguint(self) -> BigUint;
        pub fn reduce(self, a: u64) -> u64;
        pub fn reduce_i64(self, a: i64) -> u64;                  // least nonnegative residue
        pub fn reduce_bigint(self, a: &BigInt) -> u64;           // least nonnegative residue
        pub fn add(self, a: u64, b: u64) -> u64;
        pub fn sub(self, a: u64, b: u64) -> u64;
        pub fn neg(self, a: u64) -> u64;
        pub fn mul(self, a: u64, b: u64) -> u64;
        pub fn pow(self, a: u64, e: u64) -> u64;                 // pow(_, 0) = 1
        pub fn inv(self, a: u64) -> Result<u64, FpPolyError>;    // DivisionByZero if a ≡ 0
    }

    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub struct FpPoly { /* private: modulus: PrimeModulus, coeffs: Vec<u64> */ }
    impl FpPoly {
        pub fn zero(m: PrimeModulus) -> Self;
        pub fn one(m: PrimeModulus) -> Self;
        pub fn constant(m: PrimeModulus, c: u64) -> Self;
        pub fn monomial(m: PrimeModulus, c: u64, degree: usize) -> Self;   // c·x^degree
        pub fn from_coeffs(m: PrimeModulus, coeffs: &[u64]) -> Self;
        pub fn from_int_poly(f: &IntPoly, m: PrimeModulus) -> Self;
        pub fn to_int_poly(&self) -> IntPoly;                    // representatives in [0, p)
        pub fn modulus(&self) -> PrimeModulus;
        pub fn coeffs(&self) -> &[u64];
        pub fn degree(&self) -> Option<usize>;
        pub fn leading_coeff(&self) -> Option<u64>;
        pub fn is_zero(&self) -> bool;
        pub fn is_monic(&self) -> bool;
        pub fn eval(&self, x: u64) -> u64;
        pub fn derivative(&self) -> Self;
        pub fn neg(&self) -> Self;
        pub fn scale(&self, c: u64) -> Self;
        pub fn add(&self, rhs: &Self) -> Result<Self, FpPolyError>;
        pub fn sub(&self, rhs: &Self) -> Result<Self, FpPolyError>;
        pub fn mul(&self, rhs: &Self) -> Result<Self, FpPolyError>;
        pub fn div_rem(&self, divisor: &Self) -> Result<(Self, Self), FpPolyError>;
        pub fn rem(&self, divisor: &Self) -> Result<Self, FpPolyError>;
        pub fn to_monic(&self) -> Result<Self, FpPolyError>;     // ZeroPolynomial
        pub fn gcd(&self, rhs: &Self) -> Result<Self, FpPolyError>;
        pub fn extended_gcd(&self, rhs: &Self) -> Result<ExtendedGcd, FpPolyError>;
        pub fn mul_mod(&self, rhs: &Self, modulus: &Self) -> Result<Self, FpPolyError>;
        pub fn pow_mod(&self, exp: &BigUint, modulus: &Self) -> Result<Self, FpPolyError>;
        pub fn is_squarefree(&self) -> bool;
        pub fn is_irreducible(&self) -> bool;
        pub fn distinct_root_count(&self) -> Result<usize, FpPolyError>; // ZeroPolynomial
        pub fn roots(&self, limits: RootSearchLimits) -> Result<RootEnumeration, FpPolyError>;
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct ExtendedGcd { pub gcd: FpPoly, pub s: FpPoly, pub t: FpPoly } // s·a + t·b = gcd

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct RootSearchLimits { pub max_evaluations: u64 }

    #[derive(Debug, Clone, PartialEq, Eq)]
    #[must_use]
    pub enum RootEnumeration {
        Complete { roots: Vec<u64>, evaluations: u64 },
        Truncated { found: Vec<u64>, distinct_root_count: usize, evaluations: u64 },
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum FpPolyError {
        NotPrime { candidate: u64 },
        ModulusMismatch { left: u64, right: u64 },
        DivisionByZero,
        ZeroPolynomial,
    }
    // Display and std::error::Error are implemented by hand: thiserror is not a numfield dependency.

Semantics consumers may rely on:

- **Canonical representation.** Coefficients are least-significant first, each in [0, p), trailing
  zeros trimmed, zero polynomial = empty slice. Fields are private, so every value is canonical.
  Equality and hashing include the modulus. Every constructor reduces its input, using least
  nonnegative residues for negative integers. Reducing a non-monic integer polynomial may lower
  its degree. Consumers that require preservation check `degree()`: F is monic, and f keeps
  degree d at q not dividing c.
- **Overflow safety.** Scalar sums and products widen to u128 and reduce immediately. Polynomial
  products reduce after every multiply-add and never accumulate unreduced u128 sums. No wrapping
  or unchecked arithmetic is allowed. Results are exact for every prime p < 2^64, including
  2^64−59. `derivative` multiplies by i mod p, so derivative of x^p is zero.
- **Checked operations.** A binary operation on different moduli returns `ModulusMismatch`. A
  zero divisor or modulus polynomial returns `DivisionByZero`. No input panics. `div_rem` satisfies
  a = q·b + r with r = 0 or deg r < deg b. `gcd` is monic, except gcd(0,0) = 0. `extended_gcd`
  returns the same canonical gcd with s·a + t·b = gcd exactly. `mul_mod` and `pow_mod` reduce
  their results modulo `modulus`; a nonzero constant modulus gives zero.
  `pow_mod(_, 0, m)` = 1 rem m. Arbitrary-precision exponents admit (p^d−1)/2 and p^d without u64
  overflow; inputs need not already be reduced.
- **Deterministic u64 primality.** Return false for n < 2. Trial-divide by the first twelve primes
  2..37; n equal to one of them is prime. Otherwise run a strong-probable-prime test to all twelve
  bases, using u128 multiplication. This is a proof, not probable primality, for all u64. The least
  strong pseudoprime to all twelve bases is psi_12 = 318665857834031151167461 > 2^64 (Sorenson and
  Webster, "Strong pseudoprimes to twelve prime bases", Math. Comp. 86 (2017) 985–1003). Fewer
  prime bases are insufficient: psi_9 = psi_10 = psi_11 = 3825123056546413051 < 2^64 (Jiang and
  Deng, Math. Comp. 83 (2014)). Both bounds match the table in `shared/numth/src/prime.rs`.
  Calling that crate would be a new numfield dependency, so this small duplicate is intentional; a
  later consolidation may remove it. A BigInt candidate above u64::MAX fails `u64::try_from`.
  Consumers report that machine-prime range exhaustion as a resource outcome, not invalid input.
- **Squarefreeness and irreducibility.** `is_squarefree`: zero is false, a nonzero constant is
  true, otherwise gcd(f, f') = 1. A positive-degree p-th power has f' = 0 and is not squarefree.
  `is_irreducible` uses Rabin's exact test. Degree < 1 is false, because zero and units are not
  irreducible; degree 1 is true. For g = to_monic(f) of degree n, g is irreducible iff
  x^(p^n) ≡ x mod g and gcd(x^(p^(n/r)) − x, g) = 1 for each prime r dividing n. Compute it by n
  successive p-th powers modulo g. The first condition rejects repeated factors because
  x^(p^n) − x is squarefree. A root test never substitutes for this test.
- **Exact root count.** `distinct_root_count` = deg gcd(f, (x^p mod f) − x). It uses
  `pow_mod`, which costs O(log p) products of degree < deg f. That cost is polynomial and not
  budgeted. Complete splitting into distinct linear factors holds iff `is_squarefree()` and
  `distinct_root_count() == degree()`. Consumers decide splitting this way, never from a list.
- **Resource-accounted exhaustive enumeration.** `roots` returns `ZeroPolynomial` for zero, because
  every element would be a root. It computes the exact count first. A count of 0 returns
  `Complete { roots: [], evaluations: 0 }`. Otherwise it evaluates a = 0, 1, 2, ..., in ascending
  order, one evaluation per candidate. It stops at the count-th root with `Complete`, giving
  distinct ascending roots, or at `max_evaluations` with `Truncated`. In `Truncated`, `found` lists
  exactly the roots in [0, evaluations), count is exact, and `found.len() < count`. It never
  exceeds the budget or uses randomness, and it offers no resume API; retry with a larger budget.
  The budget makes a large machine prime a resource outcome rather than an O(p) hang.
- **Error classes for consumers.** `NotPrime` is invalid input. `ModulusMismatch` is a caller
  invariant violation. `DivisionByZero` and `ZeroPolynomial` are invalid arguments. None is a
  mathematical certificate. Truncated enumeration is resource exhaustion, never proof of absence.

KATs, named in `shared/numfield/tests/fp_poly_kat.rs`, with independent comments justifying each
fixture: `kat_prime_modulus_accepts_small_and_boundary_primes` (2, 3, 5, 7, 2^31−1, 2^61−1,
2^64−59); `kat_prime_modulus_rejects_units_and_composites` (0, 1, 4, 9, 561, 2^64−1);
`kat_prime_modulus_rejects_strong_pseudoprimes` (2047, 1373653, 25326001, 3215031751,
2152302898747, 3474749660383, 341550071728321, 3825123056546413051; each test supplies a
nontrivial divisor and checks it by product, independently of Miller–Rabin);
`kat_is_prime_u64_matches_trial_division_below_100000`;
`kat_scalar_arithmetic_near_u64_max` (p=2^64−59: (p−1)+(p−1), (p−1)², inverse, Fermat, i64::MIN
and negative BigInt reduction); `kat_canonical_representation` (reduction, trimming, negative
coefficients, degree drop such as 2x²+x+1 mod 2, modulus-sensitive equality);
`kat_modulus_mismatch_is_error`; `kat_div_rem_identity_and_degree`; `kat_gcd_and_extended_gcd`
(gcd((x−1)(x−2),(x−1)(x−3)) = x−1 mod 7, coprime gives 1, gcd(0,0)=0, and Bézout);
`kat_pow_mod_matches_repeated_multiplication_and_frobenius` (x^(p^d) ≡ x mod irreducible f; a
BigUint exponent above 2^64 agrees with its reduction mod p^d−1; exponent 0; zero modulus);
`kat_derivative_characteristic_p`; `kat_irreducible_fixtures` (x²+1/F_3, x²−2/F_3 and F_5,
x³−x−1/F_2 and F_3, x⁴+x+1/F_2, x⁶+x+1/F_2, x²+1/F_(2^61−1) with p ≡ 3 mod 4, linear);
`kat_reducible_fixtures_without_roots` ((x²+1)(x²+x+2) = x⁴+x³+x+2 over F_3; t⁴−10t²+1 ≡
(t²+2)(t²+3) over F_5; (x³+x+1)(x³+x²+1) over F_2); `kat_repeated_factor_fixtures`
((x²+x+1)² over F_2; t⁴−10t²+1 ≡ (t²+1)² over F_3; x³ over F_3);
`kat_reducible_with_roots` (x²+1/F_5, x³−x−1/F_5, x²+1/F_(2^64−59) with p ≡ 1 mod 4);
`kat_degree_and_zero_handling`; `kat_distinct_root_count` (x⁵−x/F_5 → 5, (x−1)²(x−2)/F_5 → 2,
t⁴+1/F_5 → 0, x²+1/F_(2^64−59) → 2); `kat_root_enumeration_complete` (x²−2/F_7 gives [3,4]
with 5 evaluations; x⁵−x/F_5 gives [0..4]; (x−1)²(x−2)/F_5 gives [1,2] with 3 evaluations;
t⁴+1/F_5 and constants give [] with 0 evaluations; exhaustive cross-check for small p);
`kat_root_enumeration_truncated` (x²−2/F_7 with budgets 4 and 0; x²+1/F_(2^64−59) with budget
1000 is truncated, has count 2 and 1000 evaluations, and every found value is a root).

### C-ResidueField — irreducible residue fields

**Enforcement:** compiler + test. **Defined-in:** 2. **Consumed-by:** 5, 6.
**Status:** **frozen at row 2 (commit 237d353, 2026-11-06)**; signatures transcribed from the
worker's return, module `shared_numfield::residue_field`, crate-root re-exports, no new dependency.
Checked field context, inversion, quadratic-residue test and Tonelli–Shanks. Cardinality/exponents
are arbitrary precision; nonresidue search is bounded. 24 KATs (zero, degree >1 roots, nonsquares,
invalid field inputs, resource exhaustion, v2(q-1)=62 case, exhaustive small-field cross-checks).

    pub enum ResidueFieldError { Poly(FpPolyError), ZeroModulus, ConstantModulus,
        ReducibleModulus /* incl. repeated factors */, DivisionByZero, InvariantViolation }
    impl From<FpPolyError> for ResidueFieldError;   // hand-written Display/Error
    pub struct NonresidueSearchLimits { pub max_candidates: u64 }
    #[must_use] pub enum NonresidueSearch { Found { element: FpPoly, candidates_tried: u64 },
        Absent /* characteristic 2 certificate */, BudgetExhausted { candidates_tried: u64 } }
    #[must_use] pub enum SqrtOutcome { Root(FpPoly) /* verified by squaring; other root is -r */,
        NotSquare /* proved by Euler criterion */, ResourceExhausted { candidates_tried: u64 } }
    pub struct ResidueField { /* private: monic modulus, characteristic, degree, order = p^d */ }
    impl ResidueField {
        pub fn new(f: &FpPoly) -> Result<Self, ResidueFieldError>;
        pub fn from_int_poly(f: &IntPoly, p: u64) -> Result<Self, ResidueFieldError>;
        pub fn characteristic(&self) -> PrimeModulus;  pub fn degree(&self) -> usize;
        pub fn modulus_poly(&self) -> &FpPoly;  pub fn order(&self) -> &BigUint;   // q = p^d
        pub fn zero(&self) -> FpPoly;  pub fn one(&self) -> FpPoly;
        pub fn element(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError>;   // a rem f
        pub fn add|sub|mul(&self, a: &FpPoly, b: &FpPoly) -> Result<FpPoly, ResidueFieldError>;
        pub fn neg(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError>;
        pub fn pow(&self, a: &FpPoly, exp: &BigUint) -> Result<FpPoly, ResidueFieldError>;
        pub fn inv(&self, a: &FpPoly) -> Result<FpPoly, ResidueFieldError>;
        pub fn is_square(&self, a: &FpPoly) -> Result<bool, ResidueFieldError>;   // zero is a square
        pub fn find_nonresidue(&self, limits: NonresidueSearchLimits)
            -> Result<NonresidueSearch, ResidueFieldError>;
        pub fn sqrt(&self, a: &FpPoly, limits: NonresidueSearchLimits)
            -> Result<SqrtOutcome, ResidueFieldError>;
    }

Consumer notes: elements are owned `FpPoly` values, reduced mod f on every operation; `from_int_poly`
reduces mod p so field degree drops if p divides the leading coefficient; characteristic 2 is
handled (every element a square, `Absent`) so row 5 must itself restrict to odd p and reject zero H
before `sqrt`; `sqrt` rebuilds the nonresidue per call (no cache); the root's sign is deterministic
but unspecified; for q ≡ 3 mod 4 no nonresidue search runs, so budget 0 still yields a root.
`NonresidueSearch::BudgetExhausted` / `SqrtOutcome::ResourceExhausted` are resource outcomes only,
never evidence about the input.

### C-SqrtContext — owned context, outcomes, and limits

**Enforcement:** compiler + test + prose. **Defined-in:** 3. **Consumed-by:** 4–8, 11, 13.
**Status:** **frozen at row 3 (commit e82b2eb, 2026-11-07)**, designed at the row-3 inflection.
Implemented field-for-field; 15 KATs (12 named plus 3 extras), two `compile_fail` doctests and unit
tests pass. **Implementation deviations (benign):** (a) `NonSquareCertificate::new` carries a second
`#[allow(dead_code)]` with a why-comment (no non-test caller until reconstruction lands; the first
caller can drop it if its fileset includes `context.rs`, otherwise it may stay); (b) limit structs
derive `Debug, Clone, Copy, PartialEq, Eq` but not `Hash` because `shared_numfield::RootSearchLimits`
lacks `Hash`; (c) the pair (0,0) always yields `ZeroRationalFactor` (rational check runs first), so
only degree 1 triggers `ZeroAlgebraicFactor`. The planned row 4–7 surface below is binding intent
for those rows; each freezes its own part. Final checked dependency signatures freeze at 11.

**Fileset.** Row 3 keeps its fileset exactly. `gnfs/src/sqrt/mod.rs` adds `pub mod context;`,
re-exports every public item below at `gnfs::sqrt::*`, and adds one module-list bullet. Legacy
claims are not rewritten (row 13). `gnfs/src/lib.rs` is untouched, so no crate-root re-export
exists until row 13. No `Cargo.toml` edit is needed. gnfs already depends on `shared-numfield`,
`shared-bigint`, `num-bigint`, `num-integer`, `num-rational`, and `num-traits`. No `thiserror`, so
`Display` and `Error` are hand-written. **Row 4 refinement, which the driver must apply to the
session list:** the verifier becomes the *child* module `sqrt::context::normalize` in
`gnfs/src/sqrt/context/normalize.rs`, not `gnfs/src/sqrt/normalize.rs`. Row 4 edits
`gnfs/src/sqrt/context.rs` only to declare `pub mod normalize;` and drop the seal's `dead_code`
allowance. Its other files are unchanged. Rows 5–7 are unchanged.

**Why the nesting (non-forgeability).** Rust private items are visible only in their module and its
descendants. `VerifiedRoot` has private fields and a private `seal`. The compiler therefore excludes
every module outside the `sqrt::context` subtree from minting success. That includes
`sqrt::{padic,crt,assembly,algebraic}` and all external crates and tests. Within the subtree, the
review rule is that only `NormalizedTarget::verify_candidate` calls `seal`. That one call site is
greppable. A sibling `normalize.rs` would need `pub(super)` visibility, which `padic`/`crt` could
also forge. Rejected alternative: `#[path = "normalize.rs"]` would keep the planned path, but it is
unidiomatic. Cost of the chosen design: the trusted kernel (context and normalize) shares private
access, and the module tree is one level deeper.

Signatures (row 3 implements everything in this block; bodies elided; rustdoc with `# Errors`):

    // ── classification vocabulary (C-Design failure classes) ──
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum FailureClass { InvalidInput, InvariantViolation, RationalNonSquare,
        AlgebraicNonSquare, Unsupported, IncidentalFactor, ResourceExhausted, TrivialCongruence }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum InvalidInput {
        ModulusNotGreaterThanOne { n: BigInt },
        ConstantPolynomial,                                  // f zero or degree 0
        DegreeFieldMismatch { field_degree: usize, poly_degree: usize },   // PolyPair only
        RationalPolyNotLinear { expected: IntPoly, actual: IntPoly },      // g != x - m
        RootMismatch { f_of_m_mod_n: BigInt },               // least nonnegative, nonzero
        RepeatedFactor,                                      // disc(F) = 0
        EmptyRelationSet,
        ZeroRationalFactor { index: usize },                 // a - b m = 0
        ZeroAlgebraicFactor { index: usize },                // a - b theta/c = 0 in A
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum UnsupportedDomain {
        LeadingCoefficientNotUnit,                           // gcd(c, N) = N
        ZeroDivisorFactor { index: usize },                  // nonzero nonunit relation factor
        MapDenominatorNotUnit,                               // produced by row 4 only
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum IncidentalSource { LeadingCoefficient, MapDenominator }   // rows 11/13 may add
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct IncidentalFactor { /* private: divisor, cofactor: BigInt, source */ }
    impl IncidentalFactor {
        pub fn new(n: &BigInt, divisor: BigInt, source: IncidentalSource) -> Option<Self>;
        pub fn divisor(&self) -> &BigInt;  pub fn cofactor(&self) -> &BigInt;
        pub fn source(&self) -> IncidentalSource;
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum InputRejection { Invalid(InvalidInput), Unsupported(UnsupportedDomain),
        Incidental(IncidentalFactor) }
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct InvariantViolation { pub what: &'static str }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum SqrtError { Invalid(InvalidInput), Invariant(InvariantViolation) }
    // `class(&self) -> FailureClass` on InputRejection, SqrtError; Display + Error on all errors.

    // ── checked polynomial/map context ──
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct SqrtContext { /* private: n, m, f, c, monic F, degree, theta_image, disc(F) */ }
    impl SqrtContext {
        pub fn new(f: &IntPoly, m: &BigInt, n: &BigInt) -> Result<Self, InputRejection>;
        pub fn from_poly_pair(pair: &PolyPair) -> Result<Self, InputRejection>;
        pub fn modulus(&self) -> &BigInt;          pub fn root(&self) -> &BigInt;      // m, exact
        pub fn algebraic_poly(&self) -> &IntPoly;  pub fn leading_coeff(&self) -> &BigInt;
        pub fn monic_poly(&self) -> &IntPoly;      pub fn degree(&self) -> usize;
        pub fn theta_image(&self) -> &BigInt;      // c·m mod N in [0, N)
        pub fn discriminant(&self) -> &BigInt;     // disc(F), nonzero
    }

    // ── relation-factor input boundary ──
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct RelationFactors { /* private: context: SqrtContext, pairs: Vec<(BigInt, BigInt)> */ }
    impl RelationFactors {
        pub fn new(ctx: &SqrtContext, pairs: &[(BigInt, BigInt)]) -> Result<Self, InputRejection>;
        pub fn context(&self) -> &SqrtContext;  pub fn pairs(&self) -> &[(BigInt, BigInt)];
        pub fn relation_count(&self) -> usize;  // k >= 1
    }

    // ── resource options (all Copy, plain pub fields, deliberately no Default) ──
    pub struct PrimeSearchLimits { pub first_candidate: u64, pub max_primes_examined: u64 }
    pub struct PrecisionLimits { pub max_modulus_bits: u64 }
    pub struct DirectLimits { pub prime_search: PrimeSearchLimits,
        pub nonresidue_search: NonresidueSearchLimits }
    pub struct CrtLimits { pub prime_search: PrimeSearchLimits, pub max_primes: u32,
        pub root_search: RootSearchLimits, pub nonresidue_search: NonresidueSearchLimits }
    pub struct RecombinationLimits { pub max_list_entries: u64, pub max_matches: u64 }
    pub struct SqrtLimits { pub precision: PrecisionLimits, pub direct: DirectLimits,
        pub crt: CrtLimits, pub recombination: RecombinationLimits }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum ResourceKind { InertPrimeSearch, SplitPrimeSearch, MachinePrimeRange, RootSearch,
        NonresidueSearch, Precision, RecombinationList, RecombinationMatches }
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ResourceExhaustion { pub kind: ResourceKind, pub consumed: u64, pub limit: u64 }

    // ── owned outcomes and witnesses ──
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum ReconstructionMethod { DirectNewton { prime: u64, exponent: u32 },
        SplitPrimeCrt { primes: Vec<u64>, exponent: u32 } }  // row-6 inflection may refine fields
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum NonSquareEvidence {
        InertLocalNonresidue { prime: u64 },                 // row 5
        DirectBoundExhausted { prime: u64, exponent: u32 },  // row 5: ±z at M > 2B rejected
        SplitLocalNonresidue { prime: u64, root: u64 },      // row 6
        SignSearchExhausted { primes: Vec<u64>, exponent: u32 },   // row 7: complete matching
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct NonSquareCertificate { /* private: NonSquareEvidence */ }
    impl NonSquareCertificate {
        pub(super) fn new(evidence: NonSquareEvidence) -> Self;   // sqrt-module-internal
        pub fn evidence(&self) -> &NonSquareEvidence;
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct VerifiedRoot { /* private: context, relation_count, rational_product P,
        integral_root z: IntPoly, beta: RatPoly, image y: BigInt, method */ }
    impl VerifiedRoot {
        fn seal(parts: RootParts) -> Self;   // private; RootParts is a private plain struct
        pub fn context(&self) -> &SqrtContext;     pub fn relation_count(&self) -> usize;
        pub fn rational_product(&self) -> &BigInt; // P = product(a - b m) over Z, exact m
        pub fn integral_root(&self) -> &IntPoly;   // z, deg < d, z² ≡ H mod F
        pub fn beta(&self) -> &RatPoly;            // theta-basis, deg < d, beta² = gamma
        pub fn modular_image(&self) -> &BigInt;    // y = beta(cm) mod N in [0, N)
        pub fn method(&self) -> &ReconstructionMethod;
    }
    #[must_use]
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum AlgebraicSqrtOutcome {
        Verified(VerifiedRoot), NotSquare(NonSquareCertificate), Unsupported(UnsupportedDomain),
        Incidental(IncidentalFactor), ResourceExhausted(Vec<ResourceExhaustion>),  // non-empty
    }
    // AlgebraicSqrtOutcome::class(&self) -> Option<FailureClass> (None for Verified).

Semantics consumers may rely on:

- **Validation order.** Checks stop at the first failure, in this order. (1) N>1. This check comes
  first and never reaches legacy `PolyPair::verify`, which panics for N=0 through `%`. (2) deg f ≥ 1.
  (3) For `from_poly_pair` only: the degree field, then g = x − m. (4) f(m) mod N = 0, using
  `mod_floor` after (1). (5) F squarefree, exactly when `shared_numfield::discriminant(F) != 0`; F
  is monic, so this is exact. (6) g0 = gcd(|c|, N): 1 → accept; 1<g0<N → `Incidental` with
  `LeadingCoefficient`; g0 = N → `Unsupported(LeadingCoefficientNotUnit)`. All invalid-input checks
  precede domain checks, so a malformed input is never reported as an incidental factor.
  `from_poly_pair` reads the pair's public fields and never calls `verify`, `new`,
  `number_field`, or `monic_f`. The context computes F = c^(d−1) f(t/c) itself.
- **Context meaning.** f is stored as given; no sign normalization is applied, so c may be negative,
  while F is monic regardless. m is exact, not reduced mod N, because P needs the integer m.
  Reducible squarefree F is accepted, and docs call it the quotient algebra A, not a field. The
  context makes no irreducibility claim. `PolyPair::number_field` is untouched.
- **Relation boundary.** Every factor is validated before any prime search. An empty list →
  `EmptyRelationSet`. For each index in order: a − b·m = 0 → `ZeroRationalFactor`. If
  f_hom(a,b) = Σ f_i a^i b^(d−i) = 0, then ca − bθ vanishes at a root of F, so it is not a unit in
  A. That gives `ZeroAlgebraicFactor` when the factor is the zero element: (a,b)=(0,0), or d=1.
  Otherwise it is `ZeroDivisorFactor`. The pairs form a multiset; uniqueness belongs to row 11.
  Negative and zero b are allowed. `RelationFactors` owns a clone of its context, so the binding
  cannot be mismatched.
- **Ownership.** No public type has a lifetime parameter, and none borrows `NumberField`. Witnesses,
  targets, and outcomes own `IntPoly`/`RatPoly`/`BigInt` data plus a cloned `SqrtContext`.
  `VerifiedRoot` may be cloned, which is not forgery. It has no `Default`, no public constructor,
  and no deserialization. `NonSquareCertificate` is minted only inside `gnfs::sqrt`. Its guarantee
  is crate-external, not compiler-enforced against sibling modules; rows 5–7 mint only the evidence
  their contracts name. `IncidentalFactor::new` is public because it is self-certifying: it returns
  `None` unless 1 < divisor < N and divisor | N.
- **Error vs outcome.** `Err(SqrtError)` is invalid input or an internal invariant violation, which
  is never retried as a nonsquare. `Ok(outcome)` carries mathematical/resource results. Only
  `NotSquare` certifies that no square root exists. `ResourceExhausted` and every truncated search
  are never evidence. `ResourceExhaustion` lists are chronological, so the row-7 dispatcher reports
  primary and fallback exhaustion together.
- **Resource accounting.** Prime search examines primes p ≥ max(`first_candidate`, 3) in ascending
  order, using `is_prime_u64`, deterministically. Every prime examined counts, including primes
  skipped because p | cN, F is not irreducible or not completely split mod p, or H ≡ 0. The u64 end
  → `MachinePrimeRange`. `max_primes` is ℓ's upper bound (≥ 1; ℓ = 1 allowed). Row 6 selects ℓ.
  Precision: row 4 returns `ResourceExhausted{Precision, consumed: bits(2B+1), limit}` before any
  prime search when bits(2B+1) > `max_modulus_bits`. Rows 5/6 lift to the least e with p^e > 2B.
  Newton uses the ceil-halving precision chain, so the lifting modulus has at most
  bits(2B+1)+bits(p) bits. Root and nonresidue budgets pass through unchanged to C-FpPoly and
  C-ResidueField. Row 14 pins one `SqrtLimits` value before discovery.
- **Limits routing.** Row 4 receives `PrecisionLimits`, row 5 `DirectLimits`, row 6 `CrtLimits`,
  row 7 `RecombinationLimits`, and the row-7 dispatcher the whole `SqrtLimits`. All are passed by
  value.

**Planned verifier interface (row 4 implements and freezes in C-NormalizedTarget), in
`gnfs::sqrt::context::normalize`, re-exported at `gnfs::sqrt::*`:**

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct NormalizedTarget { /* private: factors (owning context), k, h, P, gamma: RatPoly,
        eta: IntPoly, H: IntPoly, F': IntPoly, verified inverse of F': RatPoly, B_i, B */ }
    #[must_use] pub enum TargetOutcome { Ready(NormalizedTarget), ResourceExhausted(ResourceExhaustion) }
    #[must_use] pub enum CandidateVerdict { Verified(VerifiedRoot), Rejected(CandidateRejection),
        Incidental(IncidentalFactor), Unsupported(UnsupportedDomain) }
    pub enum CandidateRejection { DegreeTooHigh, CoefficientOutOfBound { index: usize },
        NotSquareRootOfTarget }   // a candidate rejection is never itself a certificate
    impl NormalizedTarget {
        pub fn new(factors: RelationFactors, limits: PrecisionLimits)
            -> Result<TargetOutcome, SqrtError>;
        pub fn verify_candidate(&self, z: &IntPoly, method: ReconstructionMethod)
            -> Result<CandidateVerdict, SqrtError>;   // the only caller of VerifiedRoot::seal
        // accessors: context, relation_count, half_count, rational_product, gamma, eta,
        // target (H), derivative, derivative_inverse, coefficient_bounds, bound, minimum_modulus
    }

`verify_candidate` checks in this order. (1) deg z < d. (2) |z_i| ≤ B_i. (3) z² ≡ H mod F over Z.
These three are rejections. (4) β = z·(c^h F')^(−1), using the verified inverse; β² = γ over Q
exactly. (5) Every reduced denominator of β is a unit mod N. A proper gcd gives `Incidental` with
`MapDenominator`; gcd = N gives `Unsupported(MapDenominatorNotUnit)`. Both of those mean a root
exists but the map fails. (6) y = β(cm) mod N. (7) y² ≡ P mod N. Failure at (4) or (7) after (3)
passed is `Err(Invariant)`, never a rejection. Rows 5–7 propose centered z, call `verify_candidate`,
and wrap `Verified` into `AlgebraicSqrtOutcome::Verified`. Row 5: at an inert prime, a rejected
centered root at M > 2B yields `DirectBoundExhausted`. Planned, and frozen at the rows themselves:
`sqrt::padic::direct_root(&NormalizedTarget, DirectLimits)` and the row-7
`sqrt::reconstruct_root(&NormalizedTarget, SqrtLimits)`. Both return
`Result<AlgebraicSqrtOutcome, SqrtError>`. Row 6's contribution/outcome types are designed at its
inflection. The row-13 rational root should take `root.rational_product()` directly, so x and y
are bound to the same P.

KATs in `gnfs/tests/sqrt_context_kat.rs`, each with an independent arithmetic comment:
`kat_context_accepts_monic_and_non_monic` (x²−2, m=23, N=527; 2x²−1, m=2, N=7 → F=t²−2,
θ↦4, disc 8); `kat_context_rejects_modulus_le_one` (N = 1, 0, −5, no panic);
`kat_context_rejects_constant_and_zero_polynomial`; `kat_context_rejects_root_mismatch`
(x²−2, m=6, N=35, f(m)=34); `kat_from_poly_pair_rejects_degree_field_and_rational_poly`
(literal PolyPair, including N=0, no panic); `kat_context_rejects_repeated_factor` ((x−3)², m=3, N=7;
non-monic (2x−1)², m=4, N=7); `kat_context_incidental_leading_coefficient_factor` (proper gcd(c,N),
divisor·cofactor = N, source LeadingCoefficient); `kat_context_invalid_precedes_incidental` (proper
gcd(c,N) but f(m) ≢ 0 → Invalid); `kat_context_leading_coefficient_not_unit_unsupported`
(5x²+x, m=5, N=5: f(5)=130, F=t²+t, gcd(5,5)=N);
`kat_context_reducible_squarefree_accepted` (x²−1, m=4, N=15);
`kat_relation_factors_boundary` (empty; a−bm=0; (1,1) on x²−1 → ZeroDivisorFactor; degree-1
zero element, for example f=x−5, m=12, N=7, pair (5,1) → ZeroAlgebraicFactor; first failing index
reported); `kat_failure_classes_distinct` (class() mapping). Add a `compile_fail` rustdoc example
on `VerifiedRoot` that shows external construction fails, and a `context.rs` unit test exercising
`seal` and the accessors. Row 3 adds no reconstruction body, so no `VerifiedRoot` is produced
outside that unit test until row 4.

### C-NormalizedTarget — exact target, bound, verifier, and map

**Enforcement:** compiler + test + prose. **Defined-in:** 4. **Consumed-by:** 5–7, 11, 13.
**Status:** **frozen at row 4 (commit e636477, 2026-11-07)**; module
`gnfs::sqrt::context::normalize` (child of `context`; re-exported at `gnfs::sqrt::*`). 12 KATs plus
unit tests (the unreachable invariant paths are exercised by unit tests that corrupt private state).
Exactly one non-test `VerifiedRoot::seal` call site, in `NormalizedTarget::verify_candidate`.
Integral eta/H, rational gamma, inverse of F', coefficient bounds, and the only successful
root-verification/map path. Public surface as implemented (owned values; no lifetimes):

    pub struct NormalizedTarget { /* private; owns RelationFactors */ }
    #[must_use] pub enum TargetOutcome { Ready(NormalizedTarget), ResourceExhausted(ResourceExhaustion) }
    #[must_use] pub enum CandidateVerdict { Verified(VerifiedRoot), Rejected(CandidateRejection),
        Incidental(IncidentalFactor), Unsupported(UnsupportedDomain) }
    pub enum CandidateRejection { DegreeTooHigh, CoefficientOutOfBound { index: usize },
        NotSquareRootOfTarget }                                  // Debug, Clone, Copy, PartialEq, Eq
    impl NormalizedTarget {
        pub fn new(factors: RelationFactors, limits: PrecisionLimits)
            -> Result<TargetOutcome, SqrtError>;
        pub fn verify_candidate(&self, z: &IntPoly, method: ReconstructionMethod)
            -> Result<CandidateVerdict, SqrtError>;
        pub fn factors(&self) -> &RelationFactors;   // additive over the plan; rows 11/13 need the pairs
        pub fn context(&self) -> &SqrtContext;  pub fn relation_count(&self) -> usize;
        pub fn half_count(&self) -> usize;  pub fn rational_product(&self) -> &BigInt;  // P
        pub fn gamma(&self) -> &RatPoly;  pub fn eta(&self) -> &IntPoly;
        pub fn target(&self) -> &IntPoly;  /* H */  pub fn derivative(&self) -> &IntPoly;
        pub fn derivative_inverse(&self) -> &RatPoly;
        pub fn coefficient_bounds(&self) -> &[BigInt];  pub fn bound(&self) -> &BigInt;  // B
        pub fn minimum_modulus(&self) -> &BigInt;       // 2B+1
    }

Consumer notes: `new` returns `ResourceExhausted{Precision,..}` before computing the inverse when
bits(2B+1) exceeds `max_modulus_bits`; `Rejected` is never a nonsquare certificate (rows 5/7 decide
bound-exhaustion themselves); `z` of degree >= d is rejected, not reduced (callers reduce mod F
first); both ±z verify; z²=H and bound failures are `Rejected`, but beta²=gamma and y²≡P failures
after z²=H passes are `Err(Invariant)`; a proper gcd on a non-unit map denominator gives
`Incidental(MapDenominator)` (preferred over an earlier unsupported denominator), gcd = N gives
`Unsupported(MapDenominatorNotUnit)`; P>0 / exact-square is NOT checked here (row 13). Bound
formulas, with reduced eta and monic F of degree d, use exact integers:

    R = 1 + max_{i<d}|F_i|; E = sum_{i<d}|eta_i| R^i; U = ceil(sqrt(E))
    B_i = d U sum_{ell=i+1}^d |F_ell| R^(ell-1-i); B = max_i B_i

M must be strictly greater than 2B. Check z²=H modulo F over Z and beta²=gamma over Q, then normalized
denominator invertibility and evaluation at cm. KATs include odd-cardinality non-monic input,
outside-order roots, zero divisors, the norm-map counterexample, and unsupported denominator maps.

### C-DirectRoot — inverse Newton reconstruction

**Enforcement:** test + prose. **Defined-in:** 5. **Consumed-by:** 7.
**Status:** **frozen at row 5 (commit 7898119, 2026-11-07)**; module `gnfs::sqrt::padic`, one
public item re-exported at `gnfs::sqrt::*` (17 integration KATs plus 8 unit tests; the per-precision
Newton invariant H·u² ≡ 1 is a unit test because the quotient ring is private):

    pub fn direct_root(target: &NormalizedTarget, limits: DirectLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;

Outcomes: `Verified` with method `DirectNewton{prime, exponent}`; `NotSquare(InertLocalNonresidue
{prime})`; `NotSquare(DirectBoundExhausted{prime, exponent})` (centered reduced z rejected at
M > 2B); `Incidental`/`Unsupported` (a root exists, the map to Z/N failed — the dispatcher must
return these, not retry); `ResourceExhausted` with a single-entry list (`InertPrimeSearch`,
`MachinePrimeRange`, `NonresidueSearch`, `Precision`). Precision limits are enforced earlier, at
`NormalizedTarget::new`, before any prime search; `direct_root` takes no precision argument. The
prime search skips even p, p | cN, F reducible mod p and H ≡ 0 mod p; every prime examined counts
against the budget. Lift exponent e = least with p^e >= 2B+1; the exponent chain has ceil(e/2)
steps. `DegreeTooHigh` from the verifier is reported as an invariant error. Row 7 consumption:
build the target with `limits.precision`, call `direct_root(&target, limits.direct)`, fall back to
the split-prime method only on `ResourceExhausted` of `InertPrimeSearch`/`MachinePrimeRange`; never
fall back on `NotSquare`, `Incidental` or `Unsupported`. Prime-search exhaustion does not certify
anything. Original statement:
bounded selection of odd p coprime to cN with irreducible F mod p and nonzero H;
u <- u(3-Hu²)/2 doubles precision with Hu²=1. Center z=Hu at M>2B and use C-NormalizedTarget.
Local obstruction or failed full-bound reconstruction certifies nonsquare; prime-search exhaustion
does not. KATs cover nonconstant/even-degree/totally imaginary roots and precision boundaries.

### C-SplitContributions — lifted roots and signed CRT contributions

**Enforcement:** compiler + test. **Defined-in:** 6. **Consumed-by:** 7.
**Status:** **frozen at row 6 (commit 5da4883, 2026-11-07)**, designed at the row-6 inflection and
implemented with **no deviation** from the signatures and record layout below (15 integration KATs
plus 7 unit tests; fixtures independently re-derived by the worker with a separate script).
Implementation notes: `root_search` is a per-prime budget, not a total; `max_primes` is an upper
bound (if the budget ends after at least one usable prime, those primes are used and no exhaustion
is reported; `max_primes == 0` gives `{SplitPrimeSearch, 0, 0}`); `p | cN` skips are redundant for
d >= 2 (zero H at root 0 already skips them) but kept explicitly; `Precision` from
`split_contributions` is reachable only on u32 exponent overflow, so user-visible precision
exhaustion comes from `NormalizedTarget::new`; the stale `dead_code` allowance on
`NonSquareCertificate::new` remains. Ground truth for row 7's KATs (verified by the worker): of the
2^n sign vectors exactly the true root and its negation assemble to ±z (527 fixture ℓ=1 and ℓ=2;
V4 quartic `F = t⁴−10t²+1`, m=4, N=97, relation pair (1,−1) twice, `z = −4−20θ+20θ²+4θ³`,
β = 1+θ, at 23⁴ and 23²·47²). A false-positive fixture for `SignSearchExhausted`: `t²−3`, pair
(2,−1), `first_candidate` 14 gives `Ready` at prime 23 with M = 529 and all 4 sign vectors are
`Rejected` by `verify_candidate`. Original statement: distinct completely split odd primes coprime to cN,
nonzero H at each root, lifted F roots and local square roots, Lagrange polynomials and CRT
idempotents, and owned signed contributions. Select Q_i=p_i^e with M=product Q_i>2B. Check both
local identities. A local nonresidue proves nonsquare.

**Fileset (confirmed, no widening).** New `gnfs/src/sqrt/crt.rs` and `gnfs/tests/sqrt_crt_kat.rs`.
`gnfs/src/sqrt/mod.rs` adds `pub mod crt;`, re-exports every public item below at `gnfs::sqrt::*`,
and adds one module-list bullet. `gnfs/Cargo.toml` adds exactly
`shared-padic = { path = "../shared/padic" }`, and `Cargo.lock` gains only that edge. shared-padic
depends on bigint and numfield, not gnfs, so no cycle is formed. Its `thiserror` 1 is already
locked. `gnfs/docs/PEDAGOGY.md` gains one split-prime section. `context.rs`, `padic.rs`,
`normalize.rs`, and `algebraic.rs` are **not** touched. `crt.rs` therefore keeps its own private
copies of the prime cursor, which reports `SplitPrimeSearch`, plus the lifting exponent and
centering. That ~50-line duplication of `padic.rs` privates is a later consolidation candidate,
not this row's work. The stale `#[allow(dead_code)]` on `NonSquareCertificate::new` remains. Row 6
is another caller, through `pub(super)`, which is visible to `sqrt::crt`. Never call the legacy
u64 `garner_crt` or other `algebraic.rs` helpers.
Row 7's listed fileset (`crt.rs`, `mod.rs`, `sqrt_crt_kat.rs`, PEDAGOGY) is consistent with this
design and needs no new file: sign matching and the split-prime entry point live in `crt.rs`, the
method-neutral dispatcher in `mod.rs`, and the KATs extend `sqrt_crt_kat.rs`.

Signatures (row 6; bodies elided; rustdoc with `# Errors`; no lifetimes; all owned):

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum ContributionSign { Plus, Minus }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct LocalPrime { /* private: prime: u64, modulus Q = p^e: BigInt, weight W: BigInt */ }
    impl LocalPrime {
        pub fn prime(&self) -> u64;
        pub fn modulus(&self) -> &BigInt;     // Q_i = p_i^e
        pub fn weight(&self) -> &BigInt;      // W_i in [0, M): ≡ 1 mod Q_i, ≡ 0 mod Q_k (k ≠ i)
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct Contribution { /* private */ }
    impl Contribution {
        pub fn prime_index(&self) -> usize;   // index into ContributionSet::locals()
        pub fn prime(&self) -> u64;
        pub fn root(&self) -> u64;            // r in [0, p), F(r) ≡ 0 mod p
        pub fn lifted_root(&self) -> &BigInt; // R in [0, Q), R ≡ r mod p, F(R) ≡ 0 mod Q
        pub fn local_root(&self) -> &BigInt;  // s in [0, Q), s² ≡ H(R) mod Q (sign unspecified)
        pub fn lagrange(&self) -> &IntPoly;   // L mod Q, coeffs in [0, Q), degree exactly d − 1
        pub fn coefficients(&self) -> &[BigInt];  // C = W·s·L mod M; exactly d entries in [0, M)
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct ContributionSet { /* private: target: NormalizedTarget (owned clone), locals,
        exponent: u32, modulus M: BigInt, contributions */ }
    impl ContributionSet {
        pub fn target(&self) -> &NormalizedTarget;     // binding cannot be mismatched
        pub fn locals(&self) -> &[LocalPrime];         // ascending primes, ell >= 1
        pub fn exponent(&self) -> u32;                 // common e
        pub fn modulus(&self) -> &BigInt;              // M = product Q_i >= 2B + 1, odd
        pub fn contributions(&self) -> &[Contribution];   // ell·d entries, canonical order
        pub fn method(&self) -> ReconstructionMethod;  // SplitPrimeCrt { primes, exponent }
        pub fn assemble(&self, signs: &[ContributionSign]) -> Result<IntPoly, SqrtError>;
    }

    #[must_use]
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum SplitOutcome {
        Ready(ContributionSet),
        NotSquare(NonSquareCertificate),              // SplitLocalNonresidue only
        ResourceExhausted(ResourceExhaustion),        // single entry, see below
    }

    pub fn split_contributions(target: &NormalizedTarget, limits: CrtLimits)
        -> Result<SplitOutcome, SqrtError>;

**Frozen C-SqrtContext shapes are unchanged.** `ReconstructionMethod::SplitPrimeCrt { primes:
Vec<u64>, exponent: u32 }` is adequate because the exponent is common to all primes. The open
refinement is therefore declined, and no frozen enum, limit, or evidence shape changes.

Semantics row 7 and the KATs may rely on:

- **Selection.** Primes are examined in ascending order from `max(first_candidate, 3)`. Each prime
  counts against `prime_search.max_primes_examined`. Prime p is usable iff it is odd, p ∤ c·N,
  F mod p is squarefree (equivalently p ∤ disc F), `distinct_root_count() == d`, the complete root
  list fits `root_search`, and H(r) ≢ 0 mod p at every root. Otherwise p is skipped. Per-prime
  order: split test; enumerate roots (`Truncated` returns `ResourceExhausted{RootSearch,
  consumed: evaluations, limit: max_evaluations}`); zero-test H at all roots, skipping p if any is
  zero; then take residue square roots for ascending roots. A `SqrtOutcome::NotSquare` at root r
  returns `NotSquare(SplitLocalNonresidue { prime, root: r })` immediately. This certificate is
  valid even for reducible A, because z ∈ Z[t] and H(r) ≡ z(r)² mod p.
  `SqrtOutcome::ResourceExhausted` returns `{NonresidueSearch, consumed: candidates_tried, limit}`.
  Residue square roots use `ResidueField` with the linear modulus t over F_p (C-ResidueField
  reuse) and constant elements.
- **ell.** `max_primes` is an upper bound. Search stops once `max_primes` usable primes are found.
  If the prime budget or machine range ends after at least one usable prime, proceed with those
  found primes, so ell < max_primes, and record no exhaustion. With zero usable primes, return
  `ResourceExhausted{SplitPrimeSearch | MachinePrimeRange, consumed: primes examined, limit:
  max_primes_examined}`. `max_primes == 0` returns `{SplitPrimeSearch, consumed: 0, limit: 0}`
  before any search. *Tradeoff:* when fewer primes exist in the budget, the whole remaining prime
  budget is spent, and ell becomes budget-dependent. In exchange, a resource failure is reported
  only when nothing usable was found. Sign dimension grows with ell·d, so pinned limits should
  prefer `max_primes = 1`. ell = 1 is a valid published parameter.
- **Precision.** e is the least value with (product p_i)^e ≥ 2B+1 (`target.minimum_modulus()`).
  If e does not fit `u32`, return `{Precision, consumed: minimum.bits(), limit: u32::MAX}`, as row 5
  does. The target already enforced `max_modulus_bits`, so M has at most
  bits(2B+1)+bits(product p_i) bits.
- **Lifting (arbitrary precision; existing scalar substrate).** Use
  `shared_padic::hensel_lift(f: &IntPoly, r0: &BigInt, p: &BigInt, k: u32) -> Result<BigInt,
  HenselError>`, which returns a value in [0, p^k). Signatures match: `IntPoly` is the same
  `shared_numfield` type, and e is `u32`. Call it with `(monic F, r, p, e)` for R. Compute
  v = H(R) mod Q with `IntPoly::eval_mod`. Then call it with `(x² − v, s0, p, e)` for s, where s0 is
  the residue root. Preconditions are established beforehand: r is simple because F is squarefree
  mod p, and 2·s0 ≢ 0 because p is odd and v ≢ 0. Never pass k = 0, which returns an unlifted value.
  Every `HenselError` is `Err(SqrtError::Invariant)`. `hensel_lift` does not verify its result, so
  row 6 itself checks F(R) ≡ 0, R ≡ r, s² ≡ H(R) mod Q, and s ≡ s0 mod p. A failure is
  `Err(Invariant)`.
- **Interpolation and CRT (in `crt.rs`, BigInt throughout).** L_ij = product_{k≠j}(t − R_ik) ·
  (product_{k≠j}(R_ij − R_ik))^(−1) mod Q_i. The denominators are units because the roots are
  distinct mod p. M_i = M/Q_i and W_i = M_i·(M_i^(−1) mod Q_i) mod M. Use checked inverses: either
  `num_integer::Integer::extended_gcd` or `shared_padic::Zp::inv`, and verify the product is 1.
  ell = 1 gives W = 1. C_ij = W_i·s_ij·L_ij mod M, padded to exactly d coefficients. Row 6 checks
  W_i ≡ 1 mod Q_i at runtime. The other identities are KAT obligations.
- **Canonical order and global sign.** Contributions are ordered by ascending prime, then ascending
  root r. Index 0 is the first prime's least root. This order is part of the contract because row 7
  fixes index 0 to `Plus` as its only global-sign normalization.
- **assemble.** `signs.len()` must equal the contribution count; otherwise the result is
  `Err(Invariant)`. It returns the centered residue of sum ε_k C_k mod M, coefficientwise, as an
  `IntPoly` of degree < d. That value can be passed directly to `verify_candidate`. It selects no
  signs. Why it is sound: for any square root z, z(R_ij) ≡ ±s_ij mod Q_i, because H(R_ij) is a unit
  and Z/p^e has only ±s as unit square roots. Also F ≡ product(t − R_ij) mod Q_i. Hence z equals
  the assembly of some sign vector, and centered at M > 2B it equals z exactly.
- **No success or sign selection here.** Row 6 never calls `verify_candidate` and never returns
  `Verified`. Local squares at every root are not evidence of a global square. That is the local
  false-positive case.

**Planned row-7 consumption interface (binding intent; frozen in C-VerifiedRoot at row 7):**

    // crt.rs
    pub fn recombine_signs(set: &ContributionSet, limits: RecombinationLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;
    pub fn split_prime_root(target: &NormalizedTarget, crt: CrtLimits,
        recombination: RecombinationLimits) -> Result<AlgebraicSqrtOutcome, SqrtError>;
    // mod.rs (method-neutral dispatcher, planned name from C-SqrtContext)
    pub fn reconstruct_root(target: &NormalizedTarget, limits: SqrtLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;

Row 7 owns all of the following. Bound-aware matching uses coefficient κ = argmin_i B_i from
`set.target().coefficient_bounds()`, with ties going to the lowest i, and b = B_κ. Let
n = contributions count and c_k = C_k[κ]. Fix ε_0 = Plus. The left half is indices
1..=⌊(n−1)/2⌋. The right half is the remaining ⌈(n−1)/2⌉ indices. The left list holds
(c_0 + signed left sums) mod M, enumerated by ascending sign mask. The right list holds signed
right sums mod M, sorted by (residue, mask). For each left x, match y in the exact residue
interval [lo, hi], where lo = (−x−b) mod M and hi = (−x+b) mod M. If lo ≤ hi, it is one range.
Otherwise, scan [lo, M−1] and then [0, hi]. Because 2b+1 ≤ M, the interval never double-counts.
Matches are examined in that order. A match assembles all coefficients through `assemble` and calls
`verify_candidate`. `Verified`, `Incidental`, and `Unsupported` return immediately; a root
exists. `Rejected(DegreeTooHigh)` is `Err(Invariant)`. Other rejections continue.

Accounting:
- Before allocating, if 2^⌈(n−1)/2⌉ > `max_list_entries`, return
  `ResourceExhausted[{RecombinationList, consumed: required (saturating u64), limit}]`.
- Before examining match `max_matches + 1`, return
  `ResourceExhausted[{RecombinationMatches, consumed: max_matches, limit}]`.
- Zero total matches is complete, even if `max_matches = 0`. Complete exhaustion returns
  `NotSquare(SignSearchExhausted { primes, exponent })`. This is sound because one of ±z has
  ε_0 = Plus and |z_κ| ≤ B_κ, so the true vector is always a match.

Dispatcher mapping:
- `direct_root` falls back only on `ResourceExhausted` of `InertPrimeSearch` or
  `MachinePrimeRange`.
- `SplitOutcome::NotSquare` becomes `AlgebraicSqrtOutcome::NotSquare`.
- `SplitOutcome::ResourceExhausted(e)` and recombination exhaustion append after the primary entry,
  in chronological order.
- `Ready` goes to `recombine_signs`.

KATs (row 6, `gnfs/tests/sqrt_crt_kat.rs`, each with independent arithmetic comments; fixtures
derived here, worker re-verifies):
- `kat_split_selection_skips_unusable_primes`: x²−2, m=23, N=527=17·31, relation (3,−2), and
  H=24+16θ. With `first_candidate` 0, the primes examined are 3, 5, 7. Primes 3 and 5 are inert, and
  the first usable prime is 7, with roots 3 and 4, H(3) ≡ 2, and H(4) ≡ 4. With `max_primes` 2,
  the primes are [7, 23]. 11, 13, and 19 are inert, 17 divides N, and 8 primes are examined. At 23,
  the roots are 5 and 18, with H ≡ 12 and 13, both squares.
- `kat_lifted_roots_satisfy_polynomial_mod_prime_power`
- `kat_local_roots_satisfy_target_identity`
- `kat_lagrange_basis_identities`: L_ij(R_ik) ≡ δ_jk, sum_j L_ij ≡ 1 mod Q_i, and deg L = d−1.
- `kat_crt_idempotents`: W_i ≡ δ_ik mod Q_k, ΣW ≡ 1, and W² ≡ W mod M. Under ell = 1, W = 1.
- `kat_every_sign_vector_is_a_local_square`: for small n, every assembly satisfies z² ≡ H
  mod (F, M). Check with `IntPoly::mul` and `pseudo_div_rem`, not `crt.rs` code.
- `kat_known_root_sign_vector_assembles`: derive ε_ij independently from z(R_ij) ≡ ±s_ij, then
  require `assemble(ε) = ±z`. Use the 527 fixture with z = 4+2θ under ell = 1 and ell = 2, and use
  a V4 target, F = t⁴−10t²+1. The V4 case may use m = 4 and N = 97 with duplicated pairs such as
  (1,−1) twice, so β = 1+θ.
- `kat_local_nonresidue_certifies_not_square`: x²−3, m=4, N=13, and relation (2,−1). Then
  γ = 2+θ and H = 24+12θ. Prime 3 is not squarefree, while 5 and 7 are inert. At 11, the roots are
  5 and 6, and H(5) ≡ 7 is a nonresidue, giving `SplitLocalNonresidue { prime: 11, root: 5 }`.
- `kat_local_false_positive_is_not_a_certificate`: use the same target with `first_candidate`
  14. 17 and 19 are inert. At 23, the roots are 7 and 16, with H ≡ 16 and 9, both squares, so the
  result is `Ready`. Row 7 must later reject it with `SignSearchExhausted`.
- `kat_split_resource_outcomes`: cover `SplitPrimeSearch` with budget 1, `max_primes` 0,
  `RootSearch` with max_evaluations 0 at a split prime, and `NonresidueSearch` with budget 0 at a
  p ≡ 1 mod 4 split prime. Also cover a partial ell when the budget ends after one usable prime.
- `kat_degree_one_every_prime_splits`: L = 1 and n = ell.
- `kat_assemble_rejects_wrong_sign_count`

The V4 field genuinely has no inert prime, so it requires the fallback. The NOTES V4 β = 1+θ+θ²
is not obviously a product of relation factors (a − bθ). Row 7's V4 KAT should therefore use a
relation-pair square in the V4 field. It may check the NOTES β identity separately as an exact
algebra comment.

**Flagged recommendations (driver surfaces; none blocks row 6):**
- (a) Add `gnfs/src/sqrt/context.rs` to row 7's fileset only to delete the stale `dead_code`
  allowance and its comment on `NonSquareCertificate::new`. While there, correct the `CrtLimits`
  doc "modulo each prime power", which should say "modulo each split prime". This is a 3-line
  anneal. Otherwise, the item remains a deferred note.
- (b) Later consolidation of the duplicated prime cursor and helpers in `padic.rs` and `crt.rs`.

### C-VerifiedRoot — bounded canonical dispatch and exact sign matching

**Enforcement:** compiler + test + prose. **Defined-in:** 7. **Consumed-by:** 13, 14.
**Status:** **frozen at row 7 (commit 01e1907, 2026-11-07)**; 30 integration KATs in
`sqrt_crt_kat.rs` (15 from row 6, 15 new) plus 15 unit tests in `crt.rs`. Public surface, all
re-exported at `gnfs::sqrt::*`, owned values:

    // crt.rs
    pub fn recombine_signs(set: &ContributionSet, limits: RecombinationLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;
    pub fn split_prime_root(target: &NormalizedTarget, crt: CrtLimits,
        recombination: RecombinationLimits) -> Result<AlgebraicSqrtOutcome, SqrtError>;
    // mod.rs
    pub fn reconstruct_root(target: &NormalizedTarget, limits: SqrtLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;
    pub fn reconstruct_root_from_factors(factors: RelationFactors, limits: SqrtLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;   // addition over the plan: builds the target
                                                      // under limits.precision, then dispatches

Semantics: matching is on kappa = argmin B_i (lowest index on ties); the only sign normalisation is
eps_0 = Plus; the left half is indices `1..=floor((n-1)/2)`, the left list is indexed by mask and the
right list sorted by (residue, mask); a wrapped modular interval is scanned as the two exact ranges
`[lo, M-1]` and `[0, hi]` by binary search; list size `2^ceil((n-1)/2)` is checked before allocation
(saturating, `try_reserve`) and matches are counted before examining match `max_matches+1`; every
match is assembled and passed through `NormalizedTarget::verify_candidate`; `Incidental` and
`Unsupported` return from the first match that yields them; `Rejected(DegreeTooHigh)` is an
invariant error; `NotSquare(SignSearchExhausted{primes, exponent})` only after a complete search
with every match rejected (a zero-match complete search counts, even with `max_matches = 0`);
list/match exhaustion is `ResourceExhausted`. Dispatcher: `direct_root` first; the split-prime
method runs ONLY if direct returns `ResourceExhausted` consisting solely of
`InertPrimeSearch`/`MachinePrimeRange`; `NotSquare`, `Incidental`, `Unsupported` and
`NonresidueSearch`/`Precision` exhaustion pass through unchanged; if both are exhausted the lists
are concatenated primary-first; `limits.precision` is checked against bits(2B+1) before any search.
Consumer notes: cost is `2^ceil((l*d-1)/2)` list entries so pinned limits should keep
`max_primes = 1`; the order of matches depends on the stored sign of each local root, so KATs pin
match counts, not order. Row 7 also dropped the stale dead-code allowance on
`NonSquareCertificate::new` and fixed the `CrtLimits` doc in `context.rs`. Original statement:
primary-to-fallback dispatch through the common exact verifier. Match
sum epsilon_ij W_i s_ij L_ij modulo M using exact modular intervals and meet-in-the-middle lists.
Fix only one global sign; no independent cross-prime sign normalization. Bound lists and matches.
Absent a local obstruction, only complete matching-branch exhaustion at the certified bound proves
nonsquare. KATs include ell=1, multiple primes, V4 quartic, local false positives, and agreement of
the two methods.

### C-CheckedRelations — validated records and matrix interpretation

**Enforcement:** compiler + test + prose. **Defined-in:** 8. **Consumed-by:** 9–11, 13, 14.
**Status:** **frozen at row 8 (commit eba0a56, 2026-11-07)**, designed at the row-8 inflection. The
fileset refinement below was applied to the session list. 15 KATs (9 named plus 6 extras) pass; the
diff to `sieve/mod.rs`, `factor_base.rs` and `filter/mod.rs` is additive or doc-comment only
(`Relation`, `SparseMatrix`, `build_matrix`, `line_sieve` and NFS-DL untouched; the legacy matrix
path is pinned by a KAT). **Implementation deviations (benign):** `MatrixLayout` carries one extra
private cached field `num_cols`; `classify_prime` returns `Exceptional` for p < 2 (never panics);
`ExponentSide` gains `Display` and `CheckedRelationsError` gains `Error::source`; check ordering
left open by the design is whole-list passes in variant order for the factor base and per-list
(rational then algebraic) for exponents; `build_checked_matrix` reads rows through a private
infallible helper. Consumer notes: `SparseMatrix` has no `PartialEq` (compare `rows` and
`col_weights` field by field); `excess()` undercounts by `omitted_columns().len()`; row 10 should
call `classify_prime(ctx, q)` for the q ∤ c·disc(F) rule and additionally needs q odd; `(a,b)` and
`(−a,−b)` are both accepted and `b <= 0` is accepted, so the (3,−2) seam fixture is usable; kernel
vectors of the checked matrix are candidates only (fixture {(2,1),(5,2)} on f = 3x²+x−1 is a kernel
vector with Norm(γ)=117 not a rational square, and exact verification returns `NotSquare`); the
doc text in `filter/matrix.rs` saying `num_cols = FactorBase::matrix_width()` holds only for the
legacy builder and row 9 should correct it. Original statement: separate checked relation context, valid
index/exponent/sign/norm/affine-label identities, matrix column layout, and conservative
exceptional-prime classification. Preserve public Relation/SparseMatrix shapes and NFS-DL consumers.
Retain rational parity/sign; omit unproved algebraic constraints at p dividing c*disc(F).

**Fileset refinement, which the driver must apply to the session list:** one new file,
`gnfs/src/sieve/checked.rs`. Row 8 becomes `gnfs/src/sieve/{checked,mod,factor_base}.rs`,
`gnfs/src/filter/mod.rs`, `gnfs/tests/checked_relations_kat.rs`, `gnfs/docs/PEDAGOGY.md`. No other
widening and no manifest edit. gnfs already depends on `shared-numfield` (`is_prime_u64`,
`IntPoly::eval_mod`). *Why a new file:* the legacy collector, `Relation`, and its tests stay
byte-for-byte in `sieve/mod.rs`, so NFS-DL risk is confined to imports. The checked code (~400
lines with unit tests) would roughly double that module. *Cost:* one more path, plus an intra-crate
`sieve -> sqrt` module edge, which is acceptable in Rust because `sqrt -> sieve` already exists.
- `sieve/checked.rs`: validation and prime classification, kept together per the row intent.
- `sieve/mod.rs`: adds `pub mod checked;`, re-exports every public item below at `gnfs::sieve::*`,
  and adds a module-list bullet. It also corrects two `Relation` doc sentences. Algebraic sign is
  **not** resolved through a real embedding, and character columns are necessary filters that do
  **not** guarantee that a square root exists. No code change touches `Relation`, `ExponentVector`,
  `RelationError`, `Relation::new`, `Relation::verify`, or the `*_row_gf2` methods.
- `factor_base.rs`: **documentation only**. Make the same QC-claim correction on
  `obstruction_count`. Restate the bad-prime section: `is_bad_prime` is informational and not consulted by the checked
  path, which classifies by exact divisibility of c·disc(F). Also state the inherited projective
  limitation: for p | c and p | b, no affine root r satisfies a ≡ b·r, so `Relation::new` returns
  `None` and such relations are never collected.
- `filter/mod.rs`: `ColumnClass`, `MatrixLayout`, `ColumnCountOverflow`, and `build_checked_matrix`
  sit beside the unchanged legacy `build_matrix`. Add a doc note: the legacy builder treats every
  factor-base ideal column as a parity constraint. `gnfs/src/lib.rs` is untouched, so callers use
  module paths `gnfs::sieve::*` and `gnfs::filter::*`, following the row-3 precedent.

**Preservation (compiler-checkable).** No field, variant, or signature changes on `Relation`,
`ExponentVector`, `FactorBase`, `AlgebraicPrime`, `MatrixRow`, `SparseMatrix`, `build_matrix`, or
`line_sieve`. NFS-DL (`gnfs::dl`) consumes only `Relation`/`FactorBase`, and rho's index calculus
has its own `Relation` type. Both compile unchanged. The checked builder never reads or mutates
`fb.obstruction_count` or `AlgebraicPrime::is_bad_prime`; it snapshots the factor base it needs.

Signatures (bodies elided; rustdoc with `# Errors`; owned values; no lifetimes; hand-written
`Display` + `Error` on every error type; no panics on any input):

    // ── gnfs/src/sieve/checked.rs ──
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum PrimeClass { Good, Exceptional }          // Exceptional iff p | c·disc(F)
    pub fn classify_prime(ctx: &SqrtContext, p: u64) -> PrimeClass;   // exact BigInt divisibility

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum ExponentSide { Rational, Algebraic }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum FactorBaseRejection {                     // `position` = slice position
        RationalNotStrictlyAscending { position: usize },
        RationalNotPrime { position: usize, value: u64 },
        IdealsNotStrictlyAscending { position: usize },        // by (p, r)
        IdealIndexMismatch { position: usize, index: usize },  // index != position
        IdealNotPrime { position: usize, p: u64 },
        IdealRootOutOfRange { position: usize, p: u64, r: u64 },
        IdealNotRoot { position: usize, p: u64, r: u64 },      // f(r) ≢ 0 mod p, f from ctx
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum RelationRejection {
        NotCoprime { gcd: BigInt },
        MalformedExponents { side: ExponentSide, position: usize }, // index not strictly
                                                                    // ascending, or exponent 0
        IndexOutOfRange { side: ExponentSide, index: usize, size: usize },
        ZeroRationalNorm,                              // a − b·m = 0
        SignMismatch { expected: bool, stored: bool },
        RationalNormMismatch { norm: BigInt },         // Π p^e ≠ |a − b·m|
        ZeroAlgebraicNorm,                             // f_hom(a, b) = 0
        AffineLabelMismatch { position: usize, p: u64, r: u64 },  // a ≢ b·r mod p
        AlgebraicNormMismatch { norm: BigInt },        // Π p^e ≠ |f_hom(a, b)|
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum CheckedRelationsError {
        Relation { index: usize, rejection: RelationRejection },
        DuplicatePair { first: usize, second: usize },
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct RelationChecker { /* private: context: SqrtContext, rational_primes: Vec<u64>,
        ideals: Vec<(u64, u64)>, classes: Vec<PrimeClass> */ }
    impl RelationChecker {
        pub fn new(ctx: &SqrtContext, fb: &FactorBase) -> Result<Self, FactorBaseRejection>;
        pub fn context(&self) -> &SqrtContext;
        pub fn rational_primes(&self) -> &[u64];        // column order
        pub fn ideals(&self) -> &[(u64, u64)];          // (p, r), column order
        pub fn ideal_class(&self, index: usize) -> Option<PrimeClass>;
        pub fn check(&self, relation: &Relation) -> Result<CheckedRelation, RelationRejection>;
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct CheckedRelation { /* private; minted only by RelationChecker::check */ }
    impl CheckedRelation {
        pub fn a(&self) -> &BigInt;  pub fn b(&self) -> &BigInt;    // gcd(a, b) = 1; b any sign
        pub fn rational_norm(&self) -> &BigInt;       // a − b·m, signed, nonzero
        pub fn algebraic_norm(&self) -> &BigInt;      // f_hom(a,b) = Σ f_i a^i b^(d−i), signed, ≠ 0
        pub fn rational_sign(&self) -> bool;          // rational_norm < 0
        pub fn rational_exponents(&self) -> &[(usize, u32)];   // strictly ascending, e ≥ 1
        pub fn algebraic_exponents(&self) -> &[(usize, u32)];  // strictly ascending, e ≥ 1
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct CheckedRelations { /* private: checker: RelationChecker, relations */ }
    impl CheckedRelations {
        pub fn new(checker: &RelationChecker, relations: &[Relation])
            -> Result<Self, CheckedRelationsError>;   // clones checker; may be empty
        pub fn checker(&self) -> &RelationChecker;
        pub fn context(&self) -> &SqrtContext;
        pub fn relations(&self) -> &[CheckedRelation]; // position i == input i == provenance i
        pub fn len(&self) -> usize;  pub fn is_empty(&self) -> bool;
    }

    // ── gnfs/src/filter/mod.rs ──
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum ColumnClass {
        RationalPrime { prime: u64 },                  // [0, R)
        GoodIdeal { prime: u64, root: u64 },           // [R, R+A), constrained
        ExceptionalIdeal { prime: u64, root: u64 },    // [R, R+A), omitted: never set
        RationalSign,                                  // R + A
        RootCharacter { slot: usize },                 // R + A + 1 + slot; zero at build
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ColumnCountOverflow;                    // R + A + 1 + K overflowed usize

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct MatrixLayout { /* private: checker: RelationChecker (owned clone), K */ }
    impl MatrixLayout {
        pub fn checker(&self) -> &RelationChecker;     // the binding: compare by ==
        pub fn rational_size(&self) -> usize;          // R
        pub fn algebraic_size(&self) -> usize;         // A
        pub fn character_count(&self) -> usize;        // K
        pub fn sign_column(&self) -> usize;            // R + A == obstruction_col_start
        pub fn obstruction_col_start(&self) -> usize;  // R + A
        pub fn obstruction_count(&self) -> usize;      // 1 + K
        pub fn num_cols(&self) -> usize;               // R + A + 1 + K
        pub fn character_column(&self, slot: usize) -> Option<usize>;   // None if slot ≥ K
        pub fn class(&self, column: usize) -> Option<ColumnClass>;     // None if ≥ num_cols
        pub fn omitted_columns(&self) -> Vec<usize>;   // ascending ExceptionalIdeal columns
        pub fn row_columns(&self, relation: &CheckedRelation) -> Option<Vec<usize>>;
        pub fn shape_matches(&self, matrix: &SparseMatrix) -> bool;
    }
    pub fn build_checked_matrix(relations: &CheckedRelations, character_count: usize)
        -> Result<(SparseMatrix, MatrixLayout), ColumnCountOverflow>;

Semantics consumers may rely on:

- **Context binding.** `RelationChecker` owns a `SqrtContext` clone. Through it, the checker
  supplies N, the exact m, f, c, and nonzero disc(F). `CheckedRelations` and `MatrixLayout` own checker
  clones, so a matrix, layout, relation set, and context are bound by equality, never by convention.
  For every nonempty subset of `relations()`, `RelationFactors::new(relations.context(), pairs)`
  succeeds. Each pair has a − b·m ≠ 0 and f_hom ≠ 0, and the latter is exactly the unit criterion.
- **Factor-base validation (`RelationChecker::new`).** The checks follow the variant order in the
  enum, and the first failure is reported. Rational primes must be strictly ascending and prime
  under `is_prime_u64`. Ideals must be strictly ascending by (p, r), with `index == position`, p
  prime, r < p, and f(r) ≡ 0 mod p. Completeness is not required, because exact norm reconstruction
  catches any missing ideal. `is_bad_prime` and `obstruction_count` are ignored.
- **Classification.** `classify_prime` is `Exceptional` iff p | c or p | disc(F), using the context
  values and exact BigInt remainders. It is not derived from `is_bad_prime`. For d = 2,
  disc(F) = disc(f), so p | c can be invisible to a discriminant-only flag. Every ideal above an
  exceptional prime is exceptional, even when its root is simple. This whole-prime rule is
  conservative and easy to prove; its cost is that some genuine constraints are dropped. For
  example, (23, 3) on x³−x−1 is a simple root but is still omitted. 2 is `Good` when 2 ∤ c·disc(F).
  Root characters additionally require odd q; that is row 10's rule. Row 10 should reuse
  `classify_prime` for "q ∤ c·disc(F)" rather than re-deriving it.
- **Why good columns are honest (derived).** For p ∤ c·disc(F), p ∤ [O : Z[θ]], so Kummer–Dedekind
  gives a degree-1 prime 𝔭 = (p, θ − c·r) for each simple root r of f mod p. Since c is a p-unit,
  v_𝔭(a − bα) = v_𝔭(ca − bθ). For coprime (a, b), p ∤ b whenever p | f_hom, because otherwise
  f_hom ≡ c·a^d ≢ 0. Then r ≡ a/b is the unique root, so v_𝔭(a − bα) = v_p(f_hom(a, b)), and no
  other ideal above p divides the element. At exceptional primes, the recorded exponent is only a
  norm-valuation attribution, not an ideal valuation, so its parity is not a valid square
  constraint. In particular, for non-monic f, Norm(γ) = Π f_hom / c^k, so p | c parity is not
  necessary.
- **Relation check order (`check`).** The first failure is reported in this order:
  1. gcd(|a|, |b|) = 1. Any sign of b is accepted.
  2. Both exponent lists are well formed, rational first: strictly ascending indices with e ≥ 1,
     then index < size.
  3. Rational side: a − b·m ≠ 0; the stored sign matches; then exact reconstruction.
  4. Algebraic side: f_hom ≠ 0; every entry, good or exceptional, satisfies the affine label
     a ≡ b·r mod p for its ideal (p, r); then exact reconstruction.

  Reconstruction divides |norm| by each p exactly e times and requires remainder 1. It never forms
  p^e from an untrusted u32, so no input causes unbounded work. Success means each listed exponent
  is the exact p-adic valuation. Rational primes are distinct; for labels, a ≡ b·r with distinct r
  at one p would force p | b and p | a. Thus success also gives one algebraic entry per prime.
  `CheckedRelations::new` checks relations in index order. After relation i passes, an exact
  equality of (a, b) with an earlier pair returns `DuplicatePair { first, second: i }`. Duplicates
  are rejected, never silently removed or renumbered; composition deduplicates beforehand.
  (a, b) and (−a, −b) are distinct records: they have opposite rational signs and the product
  −(a − bα)². They are not rejected, and collectors with b ≥ 1 never produce both.
- **Why b ≤ 0 is accepted.** It keeps the NOTES seam fixture (3, −2) on x²−2, m = 23, N = 527,
  available through the checked path for rows 11 and 13. In that fixture, a − bm = 49 and f_hom = 1.
  Its canonical-b form (−3, 2) changes γ by −1, so it is not equivalent. *Cost:* duplicate detection
  is exact-pair, not ideal-level.
- **Column layout (frozen; rows 9–11 consume it cold).** Let R = `rational_primes().len()`, A =
  `ideals().len()`, and K = `character_count`. Columns are [0, R) rational primes, [R, R+A)
  factor-base ideals in index order, R+A the rational sign, and R+A+1+slot root characters. This is
  exactly the legacy `build_matrix` ordering, with `obstruction_col_start = R+A` and
  `obstruction_count = 1+K`. Singleton removal, merging, and the existing character convention
  (sign first, then characters) therefore read it unchanged. Exceptional ideal columns are
  **retained as permanently zero columns**, as NOTES permits, rather than compacted. Indices thus
  coincide with factor-base indices, and a zero column stays zero under any XOR. *Cost:*
  `SparseMatrix::excess()` undercounts by `omitted_columns().len()`. That is conservative: pruning
  keeps more rows.
- **Rows (`row_columns`, `build_checked_matrix`).** Row i is built from
  `relations.relations()[i]`, with provenance `[i]`. Its columns, ascending, are:
  - the rational indices with odd e;
  - R + index for every **good** ideal with odd e;
  - the sign column iff a − b·m < 0.

  Exceptional entries and character slots contribute nothing. There is no algebraic-sign column and
  no raw homogeneous-norm parity or square test. `col_weights` is exact, and an empty row is kept,
  because it is itself a dependency. `num_cols = R+A+1+K`. `row_columns` returns `None` only for a
  relation whose indices exceed this layout's sizes, meaning a foreign checker. `shape_matches`
  compares `num_cols`, `obstruction_col_start`, `obstruction_count`, and `col_weights.len()` only;
  rows and provenance are row 11's checks. The builder uses checked addition for R+A+1+K.
- **Downstream use.**
  - Row 9 relies on the block boundaries and the always-zero omitted columns.
  - Row 10 selects K pairs, possibly from `relations` before building. It sets
    `character_column(slot)` either at build time, through XOR over complete provenance, or after
    filtering through the same XOR.
  - Row 11 requires `layout.checker() == relations.checker()` and `layout.shape_matches(&matrix)`.
    It recomputes each original row as `row_columns(...)` plus its own character bits, and takes P
    and γ pairs from `relations()`. Pair uniqueness is already an invariant of `CheckedRelations`.
- **Inherited limitation (documented, not repaired).** The primitive affine collector stays as is.
  Relations with a projective factor (p | c, p | b, p | f_hom) are rejected by `Relation::new` and,
  if hand-built, by `AffineLabelMismatch` or `AlgebraicNormMismatch`. No complete
  projective/exceptional-ideal factorization is claimed. Omitting exceptional parity weakens the
  candidate filter but does not affect exact square verification. Expanding the collector requires
  review.

KATs in `gnfs/tests/checked_relations_kat.rs`, each with an independent arithmetic comment and
fixtures re-derived by the worker:
- **Fixture T (non-monic).** f = 3x²+x−1, m = 3, N = 29, since f(3) = 29. Here c = 3,
  F = t²+t−3, and disc(F) = 13. `FactorBase::new(&f, 7, 17)` gives rational primes [2, 3, 5, 7]
  and ideals (3,1)#0, (13,2)#1, (17,4)#2, (17,7)#3.
  - Mod 3, f ≡ x−1. Mod 13, f has the double root 2. Mod 17, the roots are 4 and 7.
  - 2, 5, 7, and 11 have no roots.
- **Fixture C (monic cubic).** x³−x−1, m = 2, N = 5, and disc = −23. Mod 23,
  x³−x−1 ≡ (x−3)(x−10)².

Named KATs:
- `kat_classify_primes_by_leading_coefficient_and_discriminant`: under T, 3 (p | c, p ∤ 13) and 13
  are `Exceptional`; 2, 5, and 17 are `Good`. Under C, both (23,3) and (23,10) are exceptional, and
  5 is good.
- `kat_checker_rejects_malformed_factor_base`: one case per variant.
- `kat_good_and_exceptional_relations_accepted`: T pairs with hand norms:

  | Pair | a − 3b | f_hom | Labels |
  |---|---|---|---|
  | (4,1) | 1 | 51 = 3·17 | (3,1), (17,4) |
  | (7,1) | 4 | 153 = 3²·17 | (3,1), (17,7) |
  | (2,1) | −1 | 13 | (13,2) |
  | (−3,2) | −9 | 17 | (17,7) |
  | (−9,2) | −15 | 221 | (13,2), (17,4) |
  | (5,2) | −1 | 81 = 3⁴ | (3,1) |

- `kat_checked_rows_omit_exceptional_columns`: under T with K = 2, `num_cols` = 11, the sign column
  is 8, and the character columns are 9 and 10. Rows are (4,1) [6], (7,1) [7], (2,1) [8],
  (−3,2) [7,8], (−9,2) [1,2,6,8], and (5,2) [8]. Columns 4 and 5 have weight 0. Legacy
  `build_matrix` differs exactly at those columns: (4,1) [4,6], (2,1) [5,8], and (−9,2)
  [1,2,5,6,8]. (4,1) has odd v₃(f_hom) at 3 | c and no bit. {(2,1),(5,2)} XORs to zero only in the
  checked matrix, yet Norm(γ) = (13/3)(81/3) = 117 is not a rational square. So omission admits a
  candidate that exact verification must reject.
- `kat_relation_rejections_are_typed`: one case per `RelationRejection` variant:
  - `NotCoprime`: (4,2).
  - `ZeroRationalNorm`: (3,1) under T.
  - `SignMismatch`: a flipped sign.
  - `MalformedExponents`: an unsorted index, a duplicate index, or exponent 0.
  - `IndexOutOfRange`: an index past the end of either list.
  - `RationalNormMismatch`: 2¹ claimed for a norm of 4.
  - `AlgebraicNormMismatch`.
  - `AffineLabelMismatch`: (17,7) claimed on (4,1).
  - `ZeroAlgebraicNorm`: x²−1, m = 4, N = 15, pair (1,1).
  - Projective `AffineLabelMismatch`: (1,3) under T has a − 3b = −8 and f_hom = −3 with 3 | b,
    and `Relation::new` returns `None` for it.
- `kat_duplicate_pairs_rejected`: [(4,1), (7,1), (4,1)] gives `DuplicatePair{0,2}`, while (a,b)
  and (−a,−b) are both accepted.
- `kat_negative_b_seam_relation_accepted`: x²−2, m = 23, N = 527, FB(7,7), and pair (3,−2). The
  rational entry is {7: 2}, there is no algebraic entry, and the row is empty.
- `kat_line_sieve_relations_validate`: C, FB(30,30), `LineSieveConfig::new(10,3)`. Every collected
  relation checks, and `CheckedRelations::new` succeeds. Checked rows equal legacy rows (with
  `fb.obstruction_count` = 1, K = 0) except at the two (23,·) columns. This is a positive control,
  not parameter search.
- `kat_layout_accessors_and_binding`: covers `class` for every column kind,
  `character_column(K) == None`, and `shape_matches` true and false. It also checks that
  `row_columns` returns `None` for a relation from a larger foreign checker, and that
  `RelationFactors::new` succeeds on subsets.

Existing `merge_kat`, `line_sieve_kat`, `linalg_substrate_kat`, `sqrt_algebraic_kat`, and `dl_*`
KATs stay unedited and green; that is the preservation evidence.

### C-XorLift — kernel-preserving elimination

**Enforcement:** test + prose. **Defined-in:** 9. **Consumed-by:** 10, 11, 14.
**Status:** **frozen at row 9 (commit 6f304e9, 2026-11-07)**; 12 new KATs in
`gnfs/tests/merge_kat.rs` (all 4×4 and 5×3 matrices exhaustively against brute-force Gray-code
kernels, sampled 8×7 and 10×8, three-identical-rows, overlapping provenance, pruning loss
2^19 of 2^20, checked-matrix output with an empty row and omitted columns). Public API after the
change (module `gnfs::filter::merge` and `gnfs::filter::matrix`; the new items are NOT re-exported
from `filter/mod.rs` because that file was outside the row's fileset — use the module paths):

    // matrix.rs
    MatrixRow::xor_merge        // now symmetric difference of cols AND provenance (was union)
    SparseMatrix::lift_dependency(&self, selected: &[usize]) -> Option<Vec<usize>>
        // sorted XOR of selected rows' provenance; repeats cancel; None if an index is out of range
    // merge.rs
    pub fn merge_columns(SparseMatrix) -> SparseMatrix;          // signature unchanged
    pub fn prune_cliques(SparseMatrix) -> SparseMatrix;          // unchanged; delegates below
    pub fn prune_cliques_reported(SparseMatrix) -> (SparseMatrix, PruneReport);
    pub struct PruneReport { pub removed: Vec<MatrixRow>, pub excess_before: isize,
        pub excess_after: isize }  // + is_lossless()
    pub fn eliminate_column(&mut SparseMatrix, col: usize) -> Result<(), EliminationRefusal>;
    pub enum EliminationRefusal { OutOfRange, Obstruction, UnsupportedWeight(u32),
        InconsistentWeights }    // Display + Error; the matrix is unchanged on refusal

Rules: weight two → the two rows become r0⊕r1; weight three → the pivot is the lightest incident
row (lowest index on ties) and the replacements are p⊕a then p⊕b with a, b ascending; replacements
are appended; `col_weights` equals the true weights after every step; obstruction columns are never
eliminated; omitted zero columns have weight 0 and are never touched; `merge_columns` runs a
weight-2 pass then a weight-3 pass in ascending column order and is not a fixpoint. Each reduced
row equals the XOR of the original rows its provenance names, in every column including character
columns, so row 10 may fill characters after filtering by XOR over complete provenance and row 11
lifts with `lift_dependency`. Lifting is a bijection onto the kernel of the starting matrix
restricted to surviving rows; pruning and singleton removal only drop rows. `linalg/kernel.rs`
already XORs provenance. Anything that unions provenance would now be wrong. Original statement:
for incident rows r0,r1,r2, replace with r0 XOR r1 and r0 XOR r2; apply symmetric
difference to columns and provenance. Exact lifted dependencies are preserved in the restricted
kernel. Pruning is separately reported candidate loss. KATs exhaust small kernels and overlapping
provenance, including the three-identical-rows counterexample and column-weight consistency.

### C-RootCharacters — individual affine-root filters

**Enforcement:** compiler + test + prose. **Defined-in:** 10. **Consumed-by:** 11, 14.
**Status:** **frozen at row 10 (commit 8100dd3, 2026-11-07)**; items are `gnfs::linalg::*` (module
`qc` stays private); 9 new KATs plus one legacy pin in `gnfs/tests/linalg_substrate_kat.rs`; no
Cargo change; the legacy norm-character helpers `populate_qc_columns`, `select_qc_primes` and
`DEFAULT_NUM_QC` are unchanged in behaviour (docs relabelled) and old algebraic callers still use
them until the migration row. Public API after the change:

    pub struct RootCharacter;  // (q, s): prime(), root(), evaluate(&BigInt, &BigInt) -> CharacterValue
    pub enum CharacterValue { Square /* bit 0 */, NonSquare /* bit 1 */, Zero /* no bit */ }
        // bit() -> Option<bool>; a-b·s reduced mod q, Legendre by Euler's criterion
    pub struct RootCharacters;  // owns a RelationChecker clone; minted only via new
        // new(&RelationChecker, &[(q, s)]) -> Result<_, RootCharacterError>
        //   rejection order: Even, NotPrime, Exceptional (classify_prime), RootOutOfRange, NotRoot,
        //   Duplicate{first}; each as RootCharacterError::Pair{position, prime, root, reason}
        // checker(), characters(), len() /* pairs, not primes */, is_empty()
        // relation_bits(&CheckedRelations, index) -> Result<Vec<bool>, RootCharacterError>
        // bit_table(&CheckedRelations) -> Result<Vec<Vec<bool>>, _>   // [relation][slot]
        // fill_columns(&mut SparseMatrix, &MatrixLayout, &CheckedRelations) -> Result<(), _>
        //   sets each row's character columns by XOR over its COMPLETE original provenance; zero
        //   values and malformed input refused before anything changes; refuses a second fill
        //   (CharacterColumnNotZero)
        // verify_columns(&SparseMatrix, &MatrixLayout, &CheckedRelations) -> Result<(), _>
        //   CharacterMismatch{row, slot} on a provenance-XOR difference
    pub fn select_root_characters(&CheckedRelations, count: usize, RootCharacterLimits)
        -> Result<RootCharacterSelection, RootCharacterError>;
    pub struct RootCharacterLimits { pub prime_search: PrimeSearchLimits, pub root_search: RootSearchLimits }
        // Default: start 3, 1000 primes, 100 000 evaluations per prime
    pub struct RootCharacterSelection;  // characters(), report(), into_parts()
    pub struct SelectionReport;  // rejected_pairs (vanish on a retained relation), exceptional_primes,
                                 // truncated_primes (a truncated prime still offers the roots found)
    pub enum RootCharacterError { Pair{..}, Exhausted{ resource: CharacterResource::{PrimeBudget{limit}
        | PrimeRange}, requested, found, report: Box<SelectionReport> }, RelationsCheckerMismatch,
        LayoutCheckerMismatch, LayoutCharacterCount{layout, characters}, MatrixShapeMismatch,
        RelationOutOfRange{index, len}, ZeroValue{slot, prime, root, relation},
        MalformedProvenance{row}, ProvenanceOutOfRange{row, index, relations},
        CharacterColumnNotZero{row, column}, CharacterMismatch{row, slot}, Arithmetic(FpPolyError) }
        // all with Display + Error

Consumer notes: the elimination pivot is order-dependent, so a matrix filled before merging and one
filled after can differ in structure — each satisfies the provenance-XOR identity, so do not assume
the two orders give identical matrices; zero rejection depends on the relation set (a pair that
vanishes on one retained relation is rejected for the whole set); row 11's "own character bits" call
is `relation_bits`. Fixture evidence (independent script): x²−3, m=7, N=23, u=2+α has bits 00 at the
q=23 roots 7 and 16 (so those characters cannot see it is a nonsquare) and bit 1 at pairs (11,6),
(13,4), (13,9); the square 13+4α has bit 0 everywhere; exact reconstruction certifies u a nonsquare
and 13+4α a root. Original statement: (q,s) selection, Legendre(a-bs,q), pair-count semantics, zero
rejection, and XOR over original provenance. q is odd and coprime to c*disc(F); f(s)=0 mod q. KATs distinguish the
norm-one nonsquare unit from a certified root and check merged-row consistency. Characters are
necessary filters, not square certificates or sign-selection rules.

### C-Dependency — checked original-row lifting and products

**Enforcement:** compiler + test + prose. **Defined-in:** 11. **Consumed-by:** 12–14.
**Status:** **frozen at row 11 (commit e70fa8f, 2026-11-07)**, designed at the row-11 inflection
(`design-confident`). Implemented field-for-field in `gnfs::sqrt::dependency` (five public items
re-exported at `gnfs::sqrt::*`); 11 contract KATs plus one extra
(`kat_context_accepts_characters_filled_before_or_after_filtering`) and 9 unit tests pass; every
contract hand value re-derived independently and correct. **Implementation deviations (benign):**
(a) `DependencyContext` derives `Debug, Clone` only, because `SparseMatrix` has no `PartialEq`;
(b) a `None` from `layout.row_columns` for a bound relation, which the contract calls impossible,
maps to `LayoutRelationsMismatch` instead of panicking; (c) every unreachable failure inside
`validate` is `DependencyError::Invariant`, reached in unit tests by corrupting private state.
Consumer notes: `DependencyContext::new` always takes a `RootCharacters` (use
`RootCharacters::new(&checker, &[])` for K = 0); `CheckedDependency::target(limits)` is independent of
`validate`, so a precision exhaustion appears only there. The factor-attempt/report signatures below
remain frozen-but-unimplemented until row 13. Original statement: matrix/context binding, sorted unique indices, provenance bounds, per-row lifted parity
and kernel membership, exact original-row products, and final checked factor-attempt/report
signatures. KATs reject malformed contexts and verify overlap cancellation and non-monic products.
Deduplicate pairs before building a matrix; never alter provenance indexing after validation.

**Fileset.** Row 11's fileset is unchanged: `gnfs/src/sqrt/{dependency,mod}.rs`,
`gnfs/tests/sqrt_dependency_kat.rs`, `gnfs/docs/PEDAGOGY.md`. No `context.rs`, `lib.rs`, `filter/*`,
`linalg/*` or manifest edit. `sqrt/mod.rs` adds `pub mod dependency;`, re-exports the five public
items below at `gnfs::sqrt::*`, and adds one module-list bullet. Its legacy Background text is not
rewritten here; row 13 does that. Row 11 does **not** implement the factor-attempt/report types.
Their signatures are frozen below and implemented by row 13, so no producer-less public API lands
early. Intra-crate edges `sqrt -> {sieve, filter, linalg}` already exist through the legacy modules.

Signatures (module `gnfs::sqrt::dependency`; owned values, no lifetimes; hand-written `Display` +
`Error` on every error type; no panic on any input, including a `KernelVector` struct literal):

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum DependencyContextError {            // first failure, checked in this order
        LayoutRelationsMismatch,                 // layout.checker() != relations.checker()
        CharactersRelationsMismatch,             // characters.checker() != relations.checker()
        CharacterCountMismatch { layout: usize, characters: usize },
        MatrixShapeMismatch,                     // !layout.shape_matches(&matrix)
        Characters(RootCharacterError),          // bit_table: a character vanishes (ZeroValue)
        // then per row, ascending row index:
        MalformedColumns { row: usize },         // cols not strictly ascending
        ColumnOutOfRange { row: usize, column: usize },        // column >= num_cols
        EmptyProvenance { row: usize },
        MalformedProvenance { row: usize },      // not strictly ascending
        ProvenanceOutOfRange { row: usize, index: usize, relations: usize },
        LiftedParityMismatch { row: usize, column: usize },    // smallest differing column
        // after all rows:
        ColumnWeightMismatch { column: usize },  // col_weights[c] != true weight
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct DependencyContext { /* private: relations: CheckedRelations, layout: MatrixLayout,
        characters: RootCharacters, matrix: SparseMatrix,
        original_columns: Vec<Vec<usize>>   // per original relation, incl. character columns */ }
    impl DependencyContext {
        pub fn new(relations: CheckedRelations, matrix: SparseMatrix, layout: MatrixLayout,
            characters: RootCharacters) -> Result<Self, DependencyContextError>;
        pub fn context(&self) -> &SqrtContext;     pub fn relations(&self) -> &CheckedRelations;
        pub fn layout(&self) -> &MatrixLayout;     pub fn characters(&self) -> &RootCharacters;
        pub fn matrix(&self) -> &SparseMatrix;     // the filtered matrix candidates index
        pub fn original_columns(&self, relation: usize) -> Option<&[usize]>;
        pub fn validate(&self, selection: &[usize]) -> Result<CheckedDependency, DependencyError>;
        pub fn validate_kernel_vector(&self, vector: &KernelVector)
            -> Result<CheckedDependency, DependencyError>;   // = validate(&vector.row_indices)
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum DependencyRejection {               // first failure, checked in this order
        EmptySelection,
        MalformedSelection { position: usize },  // selection[position] <= selection[position-1]
        RowOutOfRange { position: usize, row: usize, rows: usize },
        NotInKernel { column: usize },           // smallest column of odd selected weight
        EmptyLift,                               // provenance cancels to the empty set
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum DependencyError { Rejected(DependencyRejection), Invariant(InvariantViolation) }
    // DependencyError::class(&self) -> FailureClass: InvalidInput | InvariantViolation.

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct CheckedDependency { /* private; minted only by DependencyContext::validate:
        selection: Vec<usize>, relations: Vec<usize>, factors: RelationFactors,
        rational_product: BigInt, gamma: RatPoly */ }
    impl CheckedDependency {
        pub fn context(&self) -> &SqrtContext;
        pub fn selection(&self) -> &[usize];       // filtered-matrix rows, as validated
        pub fn relations(&self) -> &[usize];       // lifted original indices, ascending, non-empty
        pub fn factors(&self) -> &RelationFactors; // pairs in `relations()` order; owns context
        pub fn relation_count(&self) -> usize;     // k = relations().len() >= 1
        pub fn rational_product(&self) -> &BigInt; // P = prod(a - b m) over Z, exact m
        pub fn gamma(&self) -> &RatPoly;           // prod(a - b theta/c) mod F, theta basis
        pub fn target(&self, limits: PrecisionLimits) -> Result<TargetOutcome, SqrtError>;
            // = NormalizedTarget::new(self.factors().clone(), limits)
    }

Semantics consumers may rely on:

- **Binding.** The context owns the relation set, layout, characters, and filtered matrix. They
  belong together by checker equality and the character count. The matrix belongs to them by the
  per-row lifted-parity identity. `MatrixLayout` is minted only by `build_checked_matrix`.
  Relation sets that share a checker can produce layouts that compare equal. The identity check,
  not the layout, ties the matrix to *these* relations. A mismatch that slipped through could only
  produce a candidate that exact verification rejects, never a false root.
- **Per-row lifted parity (construction).** `original_columns(i)` is
  `layout.row_columns(relations()[i])`, which always succeeds for bound relations, together with
  `character_column(slot)` for every set bit of `characters.relation_bits(relations, i)`, sorted.
  Each matrix row's `cols` must equal the XOR of `original_columns` over its provenance. This
  subsumes `verify_columns`. With K > 0 it rejects unfilled character columns. It also forces
  omitted exceptional columns to stay zero. Filling characters before or after filtering both
  pass. Empty provenance is rejected because correct XOR elimination and row drops keep
  provenance vectors linearly independent, so a genuine row never has empty provenance
  (inferred). A duplicated genuine row passes construction; it can only produce `EmptyLift`.
- **Candidates are untrusted.** `validate` re-derives everything from the selection and never
  calls `KernelVector::verify`/`expand_provenance` (the latter panics out of range). Kernel
  membership is checked in the owned filtered matrix. The lift uses symmetric difference
  (`SparseMatrix::lift_dependency` or equivalent), so overlapping provenance cancels and nothing
  unions. Lifted parity is then recomputed from `original_columns`, not from the matrix. A
  failure after membership passed is `Invariant`, never a rejection, because construction makes it
  unreachable. So is a `RelationFactors::new` failure: every nonempty subset of `CheckedRelations`
  is a valid factor set.
- **Derived consequences of a validated dependency (inferred, not re-proved as preconditions).**
  Every rational factor-base column and the sign column have even lifted weight, and rational norms
  factor completely over the validated base. So P > 0 is an exact integer square. Good-ideal and
  character parities are necessary filters only. gamma need not be a square: `CheckedDependency`
  is a candidate, not a certificate.
- **Products.** P = Π(a − b·m) over Z with the exact m. gamma is computed directly as
  Π(c·a − b·θ) mod F over Z, scaled by 1/c^k, reduced, in the θ basis. That is the same
  representation as `NormalizedTarget::gamma`, so a KAT can check equality across the two
  implementations. No raw homogeneous-norm product is formed or tested. For non-monic f, Π f_hom
  = c^k·Norm(gamma) need not be an integer square even when gamma is.
- **Uniqueness.** Pair uniqueness is the `CheckedRelations::new` invariant (`DuplicatePair`).
  Composition removes duplicates from the `Relation` list *before* `CheckedRelations::new`, keeping
  first occurrences in order. Provenance then indexes the deduplicated list. No API accepts an
  index map or renumbers provenance after validation. (a, b) and (−a, −b) stay distinct records.
- **Fake matrices.** Fully synthetic `SparseMatrix` values appear only in `#[cfg(test)]` unit tests
  in `dependency.rs`, which also cover the unreachable invariant paths by corrupting private state.
  Integration KATs start from genuine `build_checked_matrix` output. Rejection KATs apply explicit,
  named corruptions to that output; a corrupted matrix is never success evidence.

**Final checked factor-attempt/report interface (frozen here; row 13 implements it in
`sqrt/{rational,algebraic,assembly}.rs`, replacing the legacy signatures of the same names, and adds
crate-root re-exports).** It consumes only C-Dependency, C-VerifiedRoot, and two row-12 items:
`KernelVector` and `KernelSearchStatus`, defined below under C-KernelSearch.

    // rational.rs
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct RationalRoot { /* private: product P > 0, root x0 > 0 with x0² = P exactly,
        image x = x0 mod N in [0, N) */ }    // product(), root(), modular_image()
    #[must_use] #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum RationalSqrtOutcome { Root(RationalRoot), NonSquare { product: BigInt } } // P <= 0 too
    pub fn rational_sqrt(dependency: &CheckedDependency) -> RationalSqrtOutcome;
    // algebraic.rs
    pub fn algebraic_sqrt(dependency: &CheckedDependency, limits: SqrtLimits)
        -> Result<AlgebraicSqrtOutcome, SqrtError>;
        // exactly reconstruct_root_from_factors(dependency.factors().clone(), limits)
    // assembly.rs
    #[must_use] #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum CongruenceOutcome { Proper { divisor: BigInt, cofactor: BigInt }, Trivial }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum CongruenceError { ModulusNotGreaterThanOne { n: BigInt }, NotCongruent }
    pub fn factor_from_congruence(x: &BigInt, y: &BigInt, n: &BigInt)
        -> Result<CongruenceOutcome, CongruenceError>;   // N > 1, x² ≡ y² (mod N), then gcds
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct VerifiedFactor { /* private, minted only in assembly.rs: divisor, cofactor,
        rational: RationalRoot, algebraic: VerifiedRoot, relations: Vec<usize> */ }
        // divisor(), cofactor(), rational_root(), algebraic_root(), relations()
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum AttemptOutcome {
        Factor(VerifiedFactor), TrivialCongruence { x: BigInt, y: BigInt },
        RationalNonSquare { product: BigInt }, AlgebraicNonSquare(NonSquareCertificate),
        Unsupported(UnsupportedDomain), Incidental(IncidentalFactor),
        ResourceExhausted(Vec<ResourceExhaustion>),
    }   // class(&self) -> Option<FailureClass>, None for Factor
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct AttemptRecord { pub candidate: usize, pub relations: Vec<usize>,
        pub outcome: AttemptOutcome }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum FactorError {
        InvalidCandidate { candidate: usize, error: DependencyError },
        Sqrt(SqrtError), Invariant(InvariantViolation),
    }
    pub fn attempt_factor(dependency: &CheckedDependency, limits: SqrtLimits)
        -> Result<AttemptOutcome, FactorError>;
    #[must_use] #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum FactorOutcome { Factor(VerifiedFactor), IncidentalOnly(IncidentalFactor), Exhausted }
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct FactorReport { /* private: outcome, attempts: Vec<AttemptRecord>,
        search: KernelSearchStatus */ }
        // outcome(), attempts(), search(), resource_limited()
    pub fn factor(context: &DependencyContext, candidates: &[KernelVector],
        search: KernelSearchStatus, limits: SqrtLimits) -> Result<FactorReport, FactorError>;

Report semantics row 13 must implement:

- `attempt_factor` runs in this order. (1) `rational_sqrt`; `NonSquare` gives `RationalNonSquare`,
  which is typed even though it is unreachable for a validated dependency. (2) `algebraic_sqrt`;
  every non-`Verified` outcome passes through unchanged into the matching variant. (3) It requires
  `root.rational_product() == rational.product()` and `root.context() == dependency.context()`;
  a mismatch is `Invariant`. This binds x and y to the same P. (4) `factor_from_congruence(x, y,
  N)`. `NotCongruent` or a modulus error is `Invariant`, because the context and both verified
  roots make them impossible. `Trivial` gives `TrivialCongruence`, and `Proper` mints
  `VerifiedFactor`.
- **Retry-loop invariant.** `AttemptOutcome::AlgebraicNonSquare` is produced only by passing through
  the certificate returned by `algebraic_sqrt`, hence only by `reconstruct_root_from_factors` →
  `reconstruct_root`. No other code path builds or relabels a nonsquare. Invalid candidates,
  invariant failures, rejections, and resource exhaustion are never recorded as nonsquare.
- `factor` validates and attempts candidates in slice order. A `DependencyError` aborts with
  `InvalidCandidate{candidate: position}`. Computed-kernel candidates are kernel vectors by
  construction, so such a failure is an upstream defect and is never retried. The loop stops at
  the first `Factor`. `Incidental` is recorded and the search continues; `IncidentalOnly` is
  returned only when no attempt factored, so pipeline evidence never rests on an incidental gcd.
  `resource_limited()` is true iff the search status is not complete or any attempt is
  `ResourceExhausted`. `Exhausted` with `resource_limited() == false` is a complete negative
  result for the supplied candidates only, not for the kernel.
- No panic, no unverified-root wrapper, and no legacy-shaped `Option<BigInt>` `factor` retained.

**Settled at this inflection.**

1. **`NonSquareCertificate::new` stays `pub(super)`; it is not tightened.** In the `sqrt` module
   tree, Rust visibility can name only ancestors. Restricting minting to `padic`/`crt` would mean
   moving both under a new parent module, which changes frozen module paths, or adding a
   source-grepping test; both cost more than they protect. Forgery can only produce a false
   negative (a skipped candidate), never a false root. Success stays compiler-sealed. Review
   invariant, greppable at the row-14 boundary: non-test `NonSquareCertificate::new(` call sites
   exist only in `padic.rs` and `crt.rs`. `dependency.rs`, `rational.rs`, `algebraic.rs`, and
   `assembly.rs` never call it. *Cost:* it remains crate-sibling-forgeable; the final review should
   confirm that this is acceptable.
2. **Retry-loop invariant:** as above. NotSquare comes only from `reconstruct_root` /
   `reconstruct_root_from_factors`, through `algebraic_sqrt`.
3. **Row-12 candidates:** row 12 emits `Vec<KernelVector>`, each with sorted unique filtered-matrix
   row indices, plus `KernelSearchStatus`. Validation trusts neither: it rechecks ordering, range,
   and kernel membership in the context's own matrix, lifts by XOR, and recomputes parity from the
   original rows.
4. **Fileset:** row 11 is unchanged (above). Flag for the driver: row 12's `kernel` module is
   private (`mod kernel;` in `linalg/mod.rs`), so its new public items need a re-export line. Row 12
   therefore needs `gnfs/src/linalg/mod.rs` added to its fileset as an explicit refinement. Row 13's
   fileset suffices, since it only reads `gnfs::linalg::{KernelVector, KernelSearchStatus}`.

KATs in `gnfs/tests/sqrt_dependency_kat.rs`, each with an independent arithmetic comment. The
worker re-derives every hand value with an independent script. Fixture T is f = 3x²+x−1, m = 3,
N = 29, with relations in the order (4,1), (7,1), (2,1), (−3,2), (−9,2), (5,2) as rows 0–5, and
F = t²+t−3, so θ² = 3 − θ. The seam fixture is x²−2, m = 23, N = 527, with (3,−2).

- `kat_context_binds_relations_layout_characters`: T with K = 0 and K = 2, where the two pairs are
  valid characters derived by the worker and filled. Covers accessors, and checks
  `original_columns` against `row_columns` plus the character bits.
- `kat_context_rejects_binding_mismatches`: a layout and characters from another checker, a K
  mismatch, and a truncated `col_weights`.
- `kat_context_rejects_malformed_rows`: named corruptions of a genuine matrix. These are unsorted or
  duplicate cols, a column ≥ num_cols, empty provenance, unsorted or duplicate provenance,
  provenance ≥ len, and a wrong column weight. Each yields its variant and the first row.
- `kat_context_rejects_lifted_parity_mismatch`: a flipped bit in one row with the weight adjusted,
  unfilled character columns at K > 0, and a character vanishing on a relation (`Characters`).
- `kat_validate_rejects_malformed_selection`: empty, unsorted, duplicate, and out-of-range
  selections, including a `KernelVector { row_indices: vec![5, 2] }` literal with no panic. Also a
  duplicated genuine row, where selecting both copies gives `EmptyLift`.
- `kat_validate_rejects_non_kernel_selection`: T rows {0, 1} give `NotInKernel` at column 6.
- `kat_overlap_cancellation_lifts_symmetric_difference`: a pinned genuine checked matrix and one
  `eliminate_column` weight-three elimination. A kernel selection containing both replacement
  rows p⊕a and p⊕b lifts to a set without the pivot. `relations()` and P equal the hand values
  for the symmetric difference, not the union.
- `kat_exact_products_seam_monic`: selection [0] gives P = 49, gamma = 3 + 2θ, and
  `target(..).gamma()` equal to `gamma()`.
- `kat_exact_products_non_monic_without_norm_square`: T at K = 0, with two cases.
  - Selection {2, 5}: P = 1 and gamma = (96 − 29θ)/9 = 32/3 − (29/9)θ. Π f_hom = 13·81 = 1053 is
    not a square, and Norm(gamma) = 117. Validation succeeds.
  - Selection {1, 2, 3} (k = 3, h = 2): P = 4·(−9)·(−1) = 36 and gamma = (−993 − 62θ)/27. Here
    Norm(−993 − 62θ) = Π Norm(3a − bθ) = Π 3·f_hom = 3³·33813 = 912951, and Π f_hom = 33813 =
    153·17·13 is not a square.
  - Both cases check `target(..).gamma()` equality. Neither asserts squareness.
- `kat_duplicates_reconciled_before_matrix`: [(4,1), (7,1), (4,1)] gives `DuplicatePair{0,2}`.
  After deduplication with first occurrences kept, the context builds and provenance indexes the
  deduplicated list. (a,b) and (−a,−b) are both accepted.
- `kat_validate_computed_kernel_vectors`: the pinned line-sieve matrix (x³−x−1, m = 2, N = 5,
  FB(30,30), `LineSieveConfig::new(10,3)`) after `remove_singletons` and `merge_columns`, run
  through `block_lanczos` with a fixed seed. Every returned vector that passes
  `KernelVector::verify` also validates, and its lifted parity holds. If any such vector fails, that
  is a discovery to surface, not a reason to weaken the KAT.

### C-KernelSearch — deterministic bounded kernel combinations

**Enforcement:** compiler + test. **Defined-in:** 12. **Consumed-by:** 13, 14.
To be frozen at 12: enumeration ordering and work/output bounds, distinct nonzero XOR candidates,
and exhaustion accounting. Candidates originate in the computed kernel; validation remains at the
checked dependency boundary. KATs include non-basis combinations and degenerate/duplicate inputs.
**Output shape fixed by the row-11 inflection.** Row 13's frozen `factor` consumes these, so row 12
must provide them. Candidates are emitted as `Vec<KernelVector>`, each with sorted unique
filtered-matrix row indices, in the frozen deterministic order. Row 12 also provides
`gnfs::linalg::KernelSearchStatus` (`Debug, Clone, PartialEq, Eq`) with a `Complete` variant, at
least one resource-limited variant whose fields are frozen at 12, and `is_complete(&self) -> bool`.
`kernel` is a private module, so the re-export needs `gnfs/src/linalg/mod.rs`, an explicit fileset
refinement the driver applies. Row 12 does not validate candidates against relations; C-Dependency
does.

### C-FactorReport — verified dependency-to-factor boundary

**Enforcement:** compiler + test + prose. **Defined-in:** 13. **Consumed-by:** 14.
To be frozen at 13: fallible rational_sqrt/algebraic_sqrt, checked factor_from_congruence, and factor
report with deterministic candidate retry. Require P>0 and exact integer square root, verified beta
and modular image, x²=y² mod N, and proper divisor checks. Incidental gcds have distinct provenance.
Candidate exhaustion reports unresolved resource-limited attempts. KATs include the nonconstant
527 fixture, nonsquare-to-success combinations, invalid inputs, trivial gcd, and API migration.
**Signatures frozen at the row-11 inflection.** See C-Dependency, "Final checked
factor-attempt/report interface". Row 13 implements them verbatim, records any deviation here, and
freezes this contract with its KAT evidence.

### C-PipelineEvidence — generated sieve-to-factor evidence

**Enforcement:** test + prose. **Defined-in:** 14. **Consumed-by:** final independent review.
To be frozen at 14: one qualifying deterministic ordinary KAT satisfying every discovery/closure
requirement in row 14, plus the attempted/unattempted configuration record. Actual generated rows,
computed kernels, irreducibility, nonzero b participation, exact square identities, modular images,
and proper factor provenance are mandatory. No manufactured relations or incidental-gcd shortcut.
The frozen evidence is a tested scope, not an all-input or production-performance claim.

## Progress ledger

Preflight bindings: VERIFY_TYPES=`cargo check --workspace`; VERIFY_TEST=`cargo test --workspace`;
additional mandatory gates=`cargo dep-policy`, `cargo audit`, `cargo coverage`. Run mode: default
self-continue (chosen by the user at run start over `halt-at-boundaries`). Baselines re-established
green at 0e974ab on 2026-11-06 (check, test, dep-policy, audit, coverage 92.76% lines). C-Design is
accepted input, not an uncompleted producer row. Concrete interfaces await their listed freezes.
Final G.5 independent review: **pending**. Row 1 allowed extras: none; unrelated dirty
`.cargo/config.toml` and `Cargo.toml` (build-tuning edits not made by any worker) were left
unstaged and uncommitted.

| # | Session | Status | Commit | Froze |
|---|---------|--------|--------|-------|
| 1 | Checked finite-field polynomials | done | a5e168b | C-FpPoly |
| 2 | Residue-field square roots | done | 237d353 | C-ResidueField |
| 3 | Checked square-root contexts | done | e82b2eb | C-SqrtContext (dependency signatures freeze at 11) |
| 4 | Normalized targets and bounds | done | e636477 | C-NormalizedTarget |
| 5 | Direct p-adic reconstruction | done | 7898119 | C-DirectRoot |
| 6 | Split-prime local contributions | done | 5da4883 | C-SplitContributions; manifest edge gnfs -> shared-padic |
| 7 | CRT recombination and dispatch ◆ | done | 01e1907 | C-VerifiedRoot (boundary: still-on-intent; anneal clean) |
| 8 | Checked relations and algebraic parity | done | eba0a56 | C-CheckedRelations |
| 9 | Elimination and XOR provenance | done | 6f304e9 | C-XorLift |
| 10 | Individual root characters | done | 8100dd3 | C-RootCharacters |
| 11 | Checked dependencies | done | e70fa8f | C-Dependency (context/validation types; factor-attempt/report signatures stay frozen for 13) |
| 12 | Bounded kernel combinations | pending | — | — |
| 13 | Verified factor API migration | pending | — | — |
| 14 | Generated sieve-to-factor evidence ◆ | pending | — | — |

## Action-frame digest

### 1 — 2026-11-06
Discovery/flex: Inflection fork returned `design-confident` (no halt); implemented with no deviation
from the designed C-FpPoly signatures. Primality is deliberately duplicated from `shared-numth`
because numfield may take no new dependency (candidate for later consolidation).
Affected: C-FpPoly (now frozen).
Deferred: yes — rows 5/6 need their own big-integer prime-power quotient arithmetic (`PrimeModulus`
is u64-only by design); consumers must decide splitting from `is_squarefree()` plus exact
`distinct_root_count() == degree()`, never from a possibly `Truncated` root list; non-monic integer
polynomials can lose degree on reduction, so check `degree()` after `from_int_poly`.
Texture: `.cargo/config.toml` and `Cargo.toml` became dirty mid-run with build-tuning edits (job
count, test threads, dev line-tables-only debuginfo) not made by any worker; committer ran in
exact-path mode and excluded them. The gates ran partly under that config.

### 2 — 2026-11-06
Discovery/flex: Interface was left to the worker (no inflection); frozen as C-ResidueField from its
return. Nothing was missing from fp_poly. `sqrt` works in characteristic 2 and treats zero as a
square with root zero, so row 5's odd-p and nonzero-H conditions are the consumer's to enforce.
Affected: C-ResidueField (now frozen); consumers 5, 6.
Deferred: yes — row 5 must filter to odd p and reject H ≡ 0 before calling `sqrt`; `sqrt` re-searches
a nonresidue per call (cost, not correctness).
Texture: the human edited `.cargo/config.toml` again (RUST_TEST_THREADS 4→2) during the run;
committer first refused as scope drift in exact mode and succeeded after the mode was restated.
`.cargo/config.toml` remains dirty and unstaged.

### 3 — 2026-11-07
Discovery/flex: Inflection fork returned `design-confident`. It refined row 4's fileset so the
verifier is the child module `gnfs/src/sqrt/context/normalize.rs` (plus a one-line edit of
`context.rs`): success (`VerifiedRoot`) is non-forgeable only to modules inside the `sqrt::context`
subtree, enforced by the compiler and pinned by two `compile_fail` doctests. Session list and row-4
detail were updated by the driver. `NonSquareCertificate` is forgeable from sibling modules inside
the crate (weaker than success); only the crate boundary protects it.
Affected: C-SqrtContext (frozen), C-NormalizedTarget (row 4 fileset refined).
Deferred: yes — row 4 must declare `pub mod normalize;` in `context.rs`, drop the seal's
`dead_code` allowance, and be the only `seal` caller (greppable invariant: only
`NormalizedTarget::verify_candidate`); row 5's fileset excludes `context.rs`, so the second
`dead_code` allowance on `NonSquareCertificate::new` may remain; the juncture re-examines the
sibling-forgeable certificate at row 7/11 boundaries; row 11 freezes final dependency signatures.
Texture: the human committed their `.cargo/config.toml` edit (79b9531) after the second surfacing; no
dirty files remain besides PLAN's own edits.

### 4 — 2026-11-07
Discovery/flex: Clean dispatch under the refined fileset; planned verifier entry points proved
adequate. Two small additions/decisions: `NormalizedTarget::factors()` accessor (rows 11/13 need the
pairs) and a denominator tie-break (a proper factor beats an earlier unsupported denominator).
`Rejected` carries no nonsquare evidence, so rows 5 and 7 own the decision that a failed centered
root at M > 2B means bound/sign-search exhaustion. `z` must be pre-reduced mod F.
Affected: C-NormalizedTarget (frozen); consumers 5–7, 11, 13.
Deferred: yes — row 5/7 must reduce candidates mod F and only certify NotSquare from local
obstruction or full-bound failure; row 13 must check P>0 and the exact integer square itself.

### 5 — 2026-11-07
Discovery/flex: Clean dispatch; `direct_root` kept the planned signature. The precision limit is
enforced at target construction rather than in `direct_root` (the planned `DirectLimits` has no
precision field). The `NonSquareCertificate::new` dead-code allowance in `context.rs` is now stale
(row 5 is its first caller) but `context.rs` was outside row 5's fileset, so it remains.
Affected: C-DirectRoot (frozen); consumer row 7.
Deferred: yes — remove the stale `#[allow(dead_code)]` on `NonSquareCertificate::new` in a row whose
fileset includes `context.rs` (none currently; flag for the row-7 boundary anneal or an explicit
fileset update); the V4 quartic x⁴−10x²+1 has no inert prime and needs the row-6/7 fallback.
Texture: the per-precision Newton-invariant check is a unit test, not an integration KAT, because the
quotient ring is private by design.

### 6 — 2026-11-07
Discovery/flex: Inflection fork returned `design-confident`; row 6 builds only the local pieces
(`split_contributions`, `ContributionSet`, `assemble`) and never chooses signs; row 7 owns
`recombine_signs`, `split_prime_root` and the `reconstruct_root` dispatcher. Frozen shapes
(`CrtLimits`, `ReconstructionMethod::SplitPrimeCrt{primes, exponent}`, `ResourceKind`) were left
unchanged. Implemented with no deviation; the only manifest change is the approved path edge
`gnfs -> shared-padic` (Cargo.lock +1 line).
Affected: C-SplitContributions (frozen); consumer row 7 (C-VerifiedRoot).
Deferred: yes — (a) row 7's fileset needs `gnfs/src/sqrt/context.rs` added (~3 lines: drop the
stale `#[allow(dead_code)]` on `NonSquareCertificate::new`, and correct the `CrtLimits` doc from
"modulo each prime power" to "each split prime"); the driver applies this explicit fileset update
before dispatching row 7. (b) `crt.rs` and `padic.rs` privately duplicate the prime-cursor and
lifting-exponent/centering helpers (consolidation candidate, not enrolled). (c) Pinned limits at
row 14 should use `max_primes = 1` because the sign search grows with ℓ·d. (d) Row 7's V4 test must
use a relation-pair square in the same field (e.g. (1,−1) twice, m=4, N=97); NOTES' β=1+θ+θ² is
not obviously buildable from relation pairs. (e) Row 7 matching is on the coefficient with the
smallest bound, fixes the first contribution's sign, splits wrapped intervals exactly, checks list
size before allocation and counts matches before examining each; zero matches after complete
examination certifies NotSquare.
Texture: the committer wrote its own commit body and it reads "contributions are kept unsigned
here", which loosely contradicts "signed contributions" (meaning: the records carry a ± slot applied
at assembly; no global sign is chosen). Left unamended: unpushed, wording is imprecise not false.

### 7 — 2026-11-07
Discovery/flex: Implemented per the planned interface with one additive public function,
`reconstruct_root_from_factors` (the plan's dispatcher takes a prebuilt target, but the dispatch
brief builds it under the precision limits), and the module header of `crt.rs` was rewritten because
it formerly said the module never chooses signs. The explicit fileset refinement (adding
`context.rs`) was used exactly as scoped: stale dead-code allowance dropped, `CrtLimits` doc fixed.
The worker found NOTES' false-positive fixture (`t²−3`, first 14) has zero matches rather than
rejected ones, so it used a second fixture `x²−2`, relation (−5,0) (1 match at ℓ=1, 2 at ℓ=2, all
rejected) for the all-matches-rejected case; a mutation check (dropping half the wrapped interval)
failed four KATs. Both methods agree up to sign (cross-method KAT).
Affected: C-VerifiedRoot (frozen); consumers 13, 14.
Deferred: yes — rows 13/14 must pin `max_primes = 1` (list size 2^ceil((ℓd−1)/2)); the V4 quartic
fixture runs through the dispatcher fallback; certificates other than success remain
sibling-forgeable inside the crate (the boundary review should confirm that is acceptable);
`crt.rs`/`padic.rs` still duplicate small private helpers (consolidation candidate, unenrolled).

### 7 ◆ boundary — 2026-11-07
Discovery/flex: Anneal gate (driver, mechanical): the plan's denylist over `gnfs/src/sqrt`,
`shared/numfield/{src,tests}` and `gnfs/tests/sqrt_*_kat.rs` gave zero hits; the two PEDAGOGY files
have three hits, all predating the run (git blame, June 2026: gnfs PEDAGOGY "of a frozen contract",
numfield PEDAGOGY "frozen in G.A ... sub-tracks" and "sub-track"), left as unrelated historical text.
Boundary fork (run mode default self-continue; the plan had recommended `halt-at-boundaries`)
returned `still-on-intent`: both methods follow the canonical design and share one bound/verifier;
NotSquare is issued only at four places (inert local nonresidue, split local nonresidue,
`DirectBoundExhausted`, `SignSearchExhausted` after complete search); budgets give
`ResourceExhausted`; `VerifiedRoot` has one production `seal` call site (normalize.rs); the V4
quartic exhausts direct's 50-prime budget and is solved only by the fallback.
Affected: C-VerifiedRoot, C-SqrtContext.
Deferred: yes — (1) rows 13/14: keep `max_primes = 1`; the precision limit bounds bits(2B+1), not the
working modulus, which can be ~ℓ·bits(p) larger (padic up to +bits(p)), so pin precision with
headroom. (2) Row 13 must correct legacy text still saying "Y = |Norm(β)| mod N" and "Couveignes'
CRT algorithm" (`gnfs/src/sqrt/mod.rs` background; legacy PEDAGOGY around lines 2905, 3044–3480,
4447–4624) and state which entry points use verified reconstruction. (3) Row 11 inflection decides
whether to tighten `NonSquareCertificate::new` (pub(super)); until then the invariant is that
production callers stay in `padic.rs` and `crt.rs`, and row 13's retry loop takes `NotSquare` only
from `reconstruct_root`/`reconstruct_root_from_factors`. A forged certificate can only cause a false
negative, never a false root. (4) Tested-scope gaps for row 14 / final review to close or name: no
non-monic (c ≠ 1) fixture through split primes and sign recombination; cross-method agreement only
at d ≤ 2; no reducible squarefree F through sign recombination; V4 and (−5,0) match counts are
observed regression values cross-checked only against brute force built on the code's own
`assemble`. (5) NOTES' fixture `t²−3`, first 14 has zero matches, so the all-rejected case uses
`x²−2`, (−5,0). (6) `crt.rs`/`padic.rs` duplicate `PrimeCursor`, `centered` and the exponent search
(consolidation candidate, unenrolled). (7) One unnecessary-text gap: nothing yet states that
`factor`/`algebraic_sqrt` do not use the new code (row 13 owns it).
Texture: the plan had recommended `halt-at-boundaries` for the first run; the user chose default
self-continue, so this first ◆ passed without a human review stop.

### 8 — 2026-11-07
Discovery/flex: Inflection fork returned `design-confident` and refined row 8's fileset with one new
file, `gnfs/src/sieve/checked.rs` (session list and row-8 detail updated by the driver), so the
legacy collector stays untouched. Exceptional prime ideals stay in the layout as permanently-zero
columns (column indices coincide with factor-base indices); the whole prime is exceptional when
p | c·disc(F), even at a simple root, which drops some genuine constraints (a conservative
tradeoff). The worker re-verified every hand-derived fixture number with an independent script and
found one more foreign-checker fixture ((1,4) on the 3x²+x−1 fixture). Projective factors
(p | c and p | b) are rejected, which is the documented inherited collection limitation.
Affected: C-CheckedRelations (frozen); consumers 9–11, 13, 14.
Deferred: yes — row 9 must correct the `filter/matrix.rs` doc claim `num_cols =
FactorBase::matrix_width()` for the checked matrix (`layout.num_cols()`) and note the `excess()`
undercount; row 10 reuses `classify_prime` and `layout.character_column(slot)`; rows 11/13 may
reuse the {(2,1),(5,2)} fixture as one where exact verification must reject a kernel candidate;
the sibling-forgeable `NonSquareCertificate` question moves to the row-11 inflection as already
recorded.
Texture: the committer was given a verbatim title and body this time and kept both.

### 9 — 2026-11-07
Discovery/flex: Clean dispatch; the discovered bug was in the legacy merge itself — weight-three
elimination replaced three rows by their total XOR (three identical rows lost two dependencies) and
provenance merged by union, so overlapping provenance did not cancel. Both are fixed; the existing
merge tests' assertions needed no change (disjoint-provenance results are unchanged); one scratch
doc comment narrating the wrong rule was replaced. The worker added API beyond the minimum
(`lift_dependency`, `PruneReport`, `eliminate_column`) as hooks for rows 10, 11 and 14, and fixed
stale `SparseMatrix` field docs in the linear-algebra chapter (a small spill within the same
PEDAGOGY file). New items are reachable only through `gnfs::filter::merge::*` and
`gnfs::filter::matrix::*` since `filter/mod.rs` was out of fileset.
Affected: C-XorLift (frozen); consumers 10, 11, 14.
Deferred: yes — row 10 fills character columns by XOR over complete provenance after filtering and
reuses `classify_prime`; row 11 lifts dependencies with `lift_dependency` and must not union
provenance anywhere; row 13/14 should not depend on `merge_columns` being a fixpoint; any later
row wanting crate-level re-exports of the new merge items needs `filter/mod.rs` in its fileset.

### 10 — 2026-11-07
Discovery/flex: Clean dispatch. Additions beyond the plan: `relation_bits`/`bit_table`, and a
truncated root enumeration still offers the roots found while reporting the prime (never proof of
absence). Legacy norm-character helpers remain callable (docs relabelled); the new filter is not
wired into the legacy end-to-end path, which is expected until the migration row. Verified the
norm-one nonsquare unit fixture independently: it passes the q=23 characters yet is not a square, and
exact reconstruction rejects it.
Affected: C-RootCharacters (frozen); consumers 11, 14.
Deferred: yes — row 11 must use `relation_bits`/`verify_columns` for provenance identities and must
treat characters as necessary filters only; row 13 migrates callers off the legacy
`populate_qc_columns`/`select_qc_primes` norm characters; row 14 must pin `RootCharacterLimits`
(the `Default` is not a pinned configuration) and use eight pairs as planned; fill order
(before vs after merging) gives structurally different but identity-satisfying matrices.

### 11 — 2026-11-07
Discovery/flex: The row-11 inflection fork returned `design-confident` (committed as e32eca3 before
dispatch, so self-cleared with no halt); the worker implemented it with no semantic deviation. The
pinned line-sieve KAT produced 14 relations, 7 rows after singleton removal and merging, and 9
`block_lanczos` vectors (seed 42); all 9 pass `KernelVector::verify`, validate, and have even lifted
parity, so no discovery against rows 12–14. The overlap-cancellation fixture confirms symmetric
difference (relations {0,1,3}, P = 225) versus a union ({0,1,2,3}, P = −1125).
Affected: C-Dependency (frozen); consumers 12–14.
Deferred: yes — (1) row 12 must add `gnfs/src/linalg/mod.rs` to its fileset (already in the
session list) and emit sorted-unique row indices; validation never trusts them. (2) Row 13 implements
the frozen factor-attempt/report interface and must take NotSquare only from
`reconstruct_root`/`reconstruct_root_from_factors`; `NonSquareCertificate::new` stays `pub(super)`
(settled at the inflection; the final review confirms that is acceptable). (3) Row 14 builds its
`DependencyContext` with `RootCharacters::new(&checker, &[])` only if K = 0, but the plan pins eight
pairs. (4) `DependencyContext` has no `PartialEq`, so row 13/14 tests compare outputs, not contexts.
Texture: `@build` is not dispatchable as a subagent in this environment, so `@general` (Sonnet) ran the
implementation, as the plan allows. Shell permission was rejected intermittently, including once for
`@committer`, which committed on a plain retry (not a scope or secret refusal).

## Discoveries & risks

- **Confirmed; internal-continue:** existing finite-field and prime-power quotient arithmetic is
  insufficient. Rows 1–7 explicitly supply it. No numfield dependency is approved; machine-prime
  validation and bounded root enumeration are the reviewed scope, not a hidden factorization engine.
- **Confirmed; internal-continue:** scalar Hensel lifting exists; the approved `gnfs -> shared-padic`
  edge introduces no dependency cycle. Existing number-field inversions are not trusted blindly on
  reducible quotients; exact identity checks remain mandatory.
- **Confirmed; internal-continue:** public sqrt callers are local GNFS APIs/tests. The old QC prime
  helper is still consumed by the old algebraic implementation until migration. Keep that path
  compiling at intermediate commits without presenting it as corrected. No GNFS benchmark target
  exists in the surveyed manifest; unrelated benchmark documentation repairs are deferred.
- **Unverified; destructive-HALT if limits exhaust:** certified-bound cost, useful dependency yield,
  and generated-fixture feasibility. The six-configuration/point/time budgets are not resettable.
  Exhaustion produces an obstruction record, never a weakened fixture or a completion claim.
- **Inferred from D21–D23; destructive-HALT if encountered:** a mathematical contradiction or
  capability reduction invalidating a downstream contract. Coverage and final-factor-only tests
  cannot adjudicate canonical-method fidelity; independent identities and review are required.
- **Additive-reshard:** a necessary independently-green helper or caller beyond the expected files,
  while preserving approved semantics. Surface to the juncture; do not quietly absorb it. A new
  dependency/install, expanded collector, new algorithm family, or changed domain requires human
  approval, not merely a fileset update.
- **Confirmed; execution preflight blocker:** PLAN, NOTES, and ROADMAP were dirty before sharding.
  This session changes only PLAN and ROADMAP and preserves NOTES. Reconcile and commit the reviewed
  documentation explicitly before launching `/plan-run`; its PLAN-only initialization exception
  does not authorize sweeping up unrelated dirty files. No commit permission is granted here.
- **Deferred, not enrolled:** NFS-DL/descent correctness; Shor/index-calculus fallback masking;
  complete projective/exceptional-ideal collection; possible Dedekind-criterion gcd defect;
  Montgomery/Nguyen variants; broad CONSOLIDATE, format/lint adoption, template handoff, and EXTEND.

## Notes for executors

### Sizing and routing

Five tuning inputs: (1) multi-crate ambient coupling pushes smaller; (2) irreducible mathematical
algorithms and atomic public API migration set a larger floor; (3) interface design errors have high
downstream cost, favouring early freezes; (4) mathematical correctness is critical, favouring
independent KAT-sized checkpoints; (5) existing substrate tests permit small increments, but weak
prior pipeline assurance prevents treating green coverage as proof. Baseline/gate costs are not yet
remeasured. The ~150–400 LOC / 2–4 file band is a default, not a cap.

Splits preserve sharp green boundaries: context versus normalized verifier (3/4), local CRT
contributions versus sign recombination (6/7), and dependency validation, pure candidate enumeration,
and public migration (11/12/13). Elimination precedes character tests so they consume correct XOR
provenance. The larger API migration is indivisible without broken callers; generated discovery is
kept with its fixture and shared budget. More green checkpoints cost more repeated gate/commit
overhead but reduce the rollback radius of a mathematical error.

No row uses F/C: after these splits, no material separable compliance tail is established. Do not
cold-transfer unresolved mathematical reasoning merely because a worker runs long. Sonnet implements;
Opus adjudicates interfaces at 1, 3, 6, 8, 11 and delivery boundaries at 7 and 14. The final boundary
also owns G.5 independent correctness review; no automatic NFS-DL enrollment follows it.

### Scope, source context, and run handoff

Required source context is the dated G.1 derivation in NOTES, including canonical references,
proofs, and independently checkable fixtures; dispatch the relevant excerpt with each row. Freeze
actual Rust signatures in the producer contract, not just a prose promise. Do not widen `Relation`,
reinterpret `PolyPair::number_field`, lower a gate, or revise an unsuccessful fixture silently.
Do not install missing tools/packages. Stop on a failing adopted baseline or an unreconciled
contract/scope change. Integration-budget exhaustion always requires human review.

After reviewing and committing the planning documents, switch to `@plan-admin` and invoke:

    /plan-run docs/PLAN.md halt-at-boundaries

This forces review at interface inflections and both delivery boundaries for the first corrective
run. This sharding session does not invoke it. Later execution authorizes only the chain's own
commits under its usual preconditions, never a push or an unrelated dirty-file commit.

### Register and anneal gate

State the property, reason, or invariant in durable text, never its plan coordinate. The root
`AGENTS.md` already provides this rule; it is unchanged. At each ◆, grep the enrolled durable files
(source, tests, and affected human docs), excluding PLAN/ROADMAP/NOTES and git metadata, for:

    \bG\.[1-5]([a-z])?(\.[0-9]+)?\b
    \bS[0-9]+\b
    sub-track|frozen contract|action-frame|◆|\bCat [ABCI]\b
    (PLAN|ROADMAP|NOTES)\.md|plan-(run|shard|juncture)

Check identifiers, prose, and filenames; translate actual planning leaks into standalone topic
language. KAT, algorithmic step numbers, and mathematical group notation such as S4 are legitimate
domain terms: classify these collisions rather than mechanically renaming them. The integrative
rows must leave the enrolled surface clean of actual plan coordinates, without broad cleanup of
unrelated historical text. Record evidence and any scope obstruction in the digest.
