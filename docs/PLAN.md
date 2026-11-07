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
| 7 | Recombine CRT signs and dispatch verified roots ◆ | I | T1; T0 boundary | C-Design, C-SqrtContext, C-NormalizedTarget, C-DirectRoot, C-SplitContributions | `gnfs/src/sqrt/{crt,mod}.rs`, `gnfs/tests/sqrt_crt_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 8 | Validate relations and classify algebraic parity — @architect | A | T1 + T0 design | C-Design, C-SqrtContext | `gnfs/src/sieve/{mod,factor_base}.rs`, `gnfs/src/filter/mod.rs`, `gnfs/tests/checked_relations_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 9 | Preserve kernels through elimination and XOR lifts | B | T1 | C-Design, C-CheckedRelations | `gnfs/src/filter/{matrix,merge}.rs`, `gnfs/tests/merge_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 10 | Filter candidates with individual root characters | B | T1 | C-Design, C-FpPoly, C-CheckedRelations, C-XorLift | `gnfs/src/linalg/{qc,mod}.rs`, `gnfs/tests/linalg_substrate_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 11 | Validate dependencies against original relations — @architect | A | T1 + T0 design | C-Design, C-SqrtContext, C-NormalizedTarget, C-CheckedRelations, C-XorLift, C-RootCharacters | `gnfs/src/sqrt/{dependency,mod}.rs`, `gnfs/tests/sqrt_dependency_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
| 12 | Enumerate bounded combinations of kernel vectors | B | T1 | C-Design, C-Dependency | `gnfs/src/linalg/kernel.rs`, `gnfs/tests/kernel_combinations_kat.rs`, `gnfs/docs/PEDAGOGY.md` |
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
To be frozen at 4: integral eta/H, rational gamma, inverse of F', coefficient bounds, and the only
successful root-verification/map path. With reduced eta and monic F of degree d, use exact integers:

    R = 1 + max_{i<d}|F_i|; E = sum_{i<d}|eta_i| R^i; U = ceil(sqrt(E))
    B_i = d U sum_{ell=i+1}^d |F_ell| R^(ell-1-i); B = max_i B_i

M must be strictly greater than 2B. Check z²=H modulo F over Z and beta²=gamma over Q, then normalized
denominator invertibility and evaluation at cm. KATs include odd-cardinality non-monic input,
outside-order roots, zero divisors, the norm-map counterexample, and unsupported denominator maps.

### C-DirectRoot — inverse Newton reconstruction

**Enforcement:** test + prose. **Defined-in:** 5. **Consumed-by:** 7.
To be frozen at 5: bounded selection of odd p coprime to cN with irreducible F mod p and nonzero H;
u <- u(3-Hu²)/2 doubles precision with Hu²=1. Center z=Hu at M>2B and use C-NormalizedTarget.
Local obstruction or failed full-bound reconstruction certifies nonsquare; prime-search exhaustion
does not. KATs cover nonconstant/even-degree/totally imaginary roots and precision boundaries.

### C-SplitContributions — lifted roots and signed CRT contributions

**Enforcement:** compiler + test. **Defined-in:** 6. **Consumed-by:** 7.
To be frozen at 6: distinct completely split odd primes coprime to cN, nonzero H at each root,
lifted F roots and local square roots, Lagrange polynomials and CRT idempotents, and owned signed
contributions. Select Q_i=p_i^e with M=product Q_i>2B. Check both local identities; a local
nonresidue proves nonsquare. KATs check interpolation, lifts, CRT weights, and resource limits.

### C-VerifiedRoot — bounded canonical dispatch and exact sign matching

**Enforcement:** compiler + test + prose. **Defined-in:** 7. **Consumed-by:** 13, 14.
To be frozen at 7: primary-to-fallback dispatch through the common exact verifier. Match
sum epsilon_ij W_i s_ij L_ij modulo M using exact modular intervals and meet-in-the-middle lists.
Fix only one global sign; no independent cross-prime sign normalization. Bound lists and matches.
Absent a local obstruction, only complete matching-branch exhaustion at the certified bound proves
nonsquare. KATs include ell=1, multiple primes, V4 quartic, local false positives, and agreement of
the two methods.

### C-CheckedRelations — validated records and matrix interpretation

**Enforcement:** compiler + test + prose. **Defined-in:** 8. **Consumed-by:** 9–11, 13, 14.
To be frozen at 8: separate checked relation context, valid index/exponent/sign/norm/affine-label
identities, matrix column layout, and conservative exceptional-prime classification. Preserve
public Relation/SparseMatrix shapes and NFS-DL consumers. Retain rational parity/sign; omit unproved
algebraic constraints at p dividing c*disc(F). KATs independently cover good and exceptional primes.

### C-XorLift — kernel-preserving elimination

**Enforcement:** test + prose. **Defined-in:** 9. **Consumed-by:** 10, 11, 14.
To be frozen at 9: for incident rows r0,r1,r2, replace with r0 XOR r1 and r0 XOR r2; apply symmetric
difference to columns and provenance. Exact lifted dependencies are preserved in the restricted
kernel. Pruning is separately reported candidate loss. KATs exhaust small kernels and overlapping
provenance, including the three-identical-rows counterexample and column-weight consistency.

### C-RootCharacters — individual affine-root filters

**Enforcement:** compiler + test + prose. **Defined-in:** 10. **Consumed-by:** 11, 14.
To be frozen at 10: (q,s) selection, Legendre(a-bs,q), pair-count semantics, zero rejection, and XOR
over original provenance. q is odd and coprime to c*disc(F); f(s)=0 mod q. KATs distinguish the
norm-one nonsquare unit from a certified root and check merged-row consistency. Characters are
necessary filters, not square certificates or sign-selection rules.

### C-Dependency — checked original-row lifting and products

**Enforcement:** compiler + test. **Defined-in:** 11. **Consumed-by:** 12–14.
To be frozen at 11: matrix/context binding, sorted unique indices, provenance bounds, per-row lifted
parity and kernel membership, exact original-row products, and final checked factor-attempt/report
signatures. KATs reject malformed contexts and verify overlap cancellation and non-monic products.
Deduplicate pairs before building a matrix; never alter provenance indexing after validation.

### C-KernelSearch — deterministic bounded kernel combinations

**Enforcement:** compiler + test. **Defined-in:** 12. **Consumed-by:** 13, 14.
To be frozen at 12: enumeration ordering and work/output bounds, distinct nonzero XOR candidates,
and exhaustion accounting. Candidates originate in the computed kernel; validation remains at the
checked dependency boundary. KATs include non-basis combinations and degenerate/duplicate inputs.

### C-FactorReport — verified dependency-to-factor boundary

**Enforcement:** compiler + test + prose. **Defined-in:** 13. **Consumed-by:** 14.
To be frozen at 13: fallible rational_sqrt/algebraic_sqrt, checked factor_from_congruence, and factor
report with deterministic candidate retry. Require P>0 and exact integer square root, verified beta
and modular image, x²=y² mod N, and proper divisor checks. Incidental gcds have distinct provenance.
Candidate exhaustion reports unresolved resource-limited attempts. KATs include the nonconstant
527 fixture, nonsquare-to-success combinations, invalid inputs, trivial gcd, and API migration.

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
| 4 | Normalized targets and bounds | pending | — | — |
| 5 | Direct p-adic reconstruction | pending | — | — |
| 6 | Split-prime local contributions | pending | — | — |
| 7 | CRT recombination and dispatch ◆ | pending | — | — |
| 8 | Checked relations and algebraic parity | pending | — | — |
| 9 | Elimination and XOR provenance | pending | — | — |
| 10 | Individual root characters | pending | — | — |
| 11 | Checked dependencies | pending | — | — |
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
