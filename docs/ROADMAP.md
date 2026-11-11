# rDLP — outstanding work

Reviewed 2026-11-08 against current source, recorded verification and mathematical references.
Reconciled 2026-11-11 at the solver-outcome follow-up boundary.
This is the single project-level register of known outstanding work. `PLAN.md` is reserved for the
currently authorised execution graph; it is not a second project backlog. Completed work and
superseded decisions belong in git history, not an accumulating status log here.

## Purpose and evidence

Make rDLP a mathematically correct, independently checked, pedagogically coherent reference library
for established factoring and discrete-logarithm algorithms. Subject coverage, algorithm fidelity,
correct final answers, and production performance are different claims; establish each separately.
Neither a green test suite nor a completed implementation graph proves all four.

**Scope boundary.** rDLP surveys DLP strategies, algorithms and phenomena; DLP topicality sets its
scope, not exhaustiveness in mathematics, code or exposition. Mathematics always offers a nearby
reason to include more, so draw the boundary between the umbra and the penumbra of the topic:
neighbouring results may be mentioned and linked, but not expounded or implemented.

**Current baseline, confirmed by source inspection:** checked GNFS reconstruction, modular mapping,
dependency validation, XOR provenance and elimination are implemented. An ordinary KAT factors
45113 from runtime-generated relations. Real opt-in PARI/GP comparisons cover reconstruction and
assembly; the completion ledger records them and the adopted workspace gates passing. This review
ran no tests, builds, benchmarks or external mathematical oracles. It found no demonstrated
wrong-result defect remaining in that checked factoring path; this is not an exhaustive proof.

Evidence labels below distinguish **confirmed** source/document facts from **inferred** mathematical
consequences or correctness concerns. Static counterexamples are identified as such, not represented
as executed regressions. **Deferred** describes scheduling, not confidence. Speculation does not
become an implementation requirement without adjudication.

Preserve these constraints while planning the work:

- Mathematical-behaviour tests are primary; the 80% line-coverage floor is complementary. Property
  tests exist in `shared/field/tests/sqrt_legendre_kat.rs`; KAT-first does not mean KAT-only.
- Keep `cargo check --workspace`, `cargo test --workspace`, `cargo dep-policy`, `cargo audit` and
  `cargo coverage` green. Re-establish their baseline before implementation. Dependency advisories
  can change independently of the source. Formatting and lint adoption remain separate work below.
- External tools validate rDLP's answers, never supply its implementation's answers. Tests requiring
  them stay opt-in. Do not treat every ignored test as an implemented oracle.
- Follow `docs/development.md`'s "Solver outcomes and fallbacks" in new and repaired code.
- Keep the existing crate layout; attacks can share their baseline's crate and arithmetic types.
  New algorithms, domain reductions, dependency edges and broad refactors require explicit review.
- Correct directly affected mathematical claims with the code. Artifact type must not force a
  correctness repair and its explanation into unrelated campaigns.
- Apply the four-surface minimalism principle below to new work as well as cleanup; do not create
  fresh explanation or machinery debt while retiring the existing backlog.

## Work register and sequencing

These are planning scopes, not an executable session graph. No new source work is authorised by
their presence here. Shard a selected scope into `PLAN.md` only after its contracts are settled.

| Scope | Current state | Dependency / priority |
|---|---|---|
| NFS-DL composition | Placeholder and unresolved maths | Design before implementation |
| Solver-outcome audit | Done 2026-11-10; residual findings below | NFS-DL and MOV k=2 via its design |
| DP walk step counts | Negation-map/batch far above expectation | Diagnose before any speedup claim |
| Oracles / CADO sidecar | Early GNFS slice selected | Design separately; parallel with repairs |
| GNFS assurance | Coverage limits, not proven defects | After urgent repairs |
| Mathematical exposition | Non-NFS-DL corrections executing in PLAN | In progress; NFS-DL claims via its design |
| Project-wide anneal | Legacy surfaces need a full pass | Audit broadly; repair in bounded units |
| Format / lint / measurements | Gates unadopted; timings stale | Dedicated changes |
| Four-surface minimalism | Explanation and machinery audit | Preserve evidence; simplify by topic |
| Template conventions | Bidirectional review | Parallel; use verified exemplars |
| Document retirement | SURVEY transfers pending | Other deletions need approval |
| Algorithm spectrum | Completeness unsupported | Reassess; extensions deferred |

The main dependency is correctness contracts → implementation and genuine composition evidence →
independent comparison and honest exposition. Oracle infrastructure can be prepared in parallel,
but cannot freeze meaningful comparisons against an undefined algorithm contract. Documentation
cleanup and template review need not wait for every algorithm repair. The user selected an early
CADO GNFS oracle slice alongside correctness repairs, not after the entire roadmap. Its design and
sharding belong in a separate session; NFS-DL comparisons wait for its repaired contract.
Extensions are not a way to avoid the outstanding correctness work.

Establish a whole-project anneal/minimalism baseline without delaying correctness repairs or the
early CADO slice. Execute the resulting cleanup in separately reviewed, coherent units, not a
repository-wide rewrite hidden in an algorithm change. Where explanation cleanup precedes structural
refactoring, do it first; refactor dependent test scaffolding after interfaces settle. Correct false
claims with their repairs. Use targeted annealing at delivery boundaries to prevent recurrence.

## Correctness

### NFS-DL contract and composition

**Confirmed:** `gnfs/src/dl/descent/solve.rs::solve_dl` discards the smoothing exponent and returns
`Ok(0)` when its k=1 frontier is empty. That is an explicit placeholder success. An empty frontier
means g^e*h=1, not that log_g(h)=0. The context-bearing `solve_dl_full` is a different API and must
not be conflated with this wrapper.

Establish one mathematical contract across collection, relation equations, linear algebra, virtual
logs and descent before applying local fixes. Inspect together:

- `dl/descent/recurse.rs::build_children`, `solve.rs::assemble_log`, and
  `dl/linalg/blockvec_fl.rs::build_fl_matrix`. **Confirmed:** assembly sums children positively,
  the matrix gives both sides positive exponents, and child construction skips all algebraic ideals
  above q by prime alone, without preserving the target exponent. **Inferred:** signs, target
  coefficient inversion and ideal identity do not consistently implement the stated equations.
- Rational frontier targets versus algebraic special-q rewriting, including the target root r,
  repeated factors, relation signs and the Schirokauer contribution.
- `dl/schirokauer.rs`, MATHEMATICS's Schirokauer exponent and claimed properties, and normalised
  virtual logs versus log_g of rational primes. These remain **unadjudicated**;
  do not present the present construction as independently established.
- Generator/subgroup, field encoding, normalisation and failure contracts. In
  `dl/ext/descent.rs::solve_dl_ext`, checking g^ell=1 alone admits g=1; it does not establish
  generation of a nontrivial prime-order subgroup.

**Confirmed scope mismatch:** the k=2 path used by MOV builds an exhaustive toy log table and
searches powers of g. It is not demonstrated NFS-DL composition. Value-bearing k=1 tests use
supplied virtual logs and smooth targets; structural descent tests with mock logs do not validate
the answer.

Closure requires value-checked composition from generated relations through computed logs and a
nontrivial descent, including backing out a nonzero smoothing exponent. Verify the final group
identity and intermediate relation equations independently. Remove placeholder success and correct
the stale k>1-unsupported and full-k=2-pipeline claims in rustdoc, MATHEMATICS and PEDAGOGY. Any
intentional capability restriction must be explicit and approved, not disguised as successful NFS.

### Shor factoring register width

**Confirmed:** `shor/src/shor/mod.rs::exp_len_for` sizes the exponent register at n_bits(N) qubits,
not the textbook ≈2·n_bits(N) that guarantees continued-fraction recovery (2^t ≥ N²); sufficiency
for the four fixture targets is empirical. Widening it pushes N=91 past the simulator's qubit
ceiling, so it is a target-selection decision, not a local fix.

### Solver-outcome follow-ups

`docs/development.md`'s "Solver outcomes and fallbacks" is the standard. Every solver in the
original audit (rho ECDLP in all four forms, Pohlig–Hellman, MOV, SSA, rho factoring, ECM,
Shor factoring, and the GHS transfer's input guard) now returns a verified, typed outcome; the
source and tests are the record. Shor ECDLP and elliptic index calculus conformed already. GHS
is a transfer, not a solver; the hyperelliptic and Semaev modules are substrates with no solver
path. Pohlig–Hellman's small-prime enumeration is the method's own size-selected sub-step
(Pohlig and Hellman 1978, §IV), not a fallback; its large-prime step is rho, not the paper's
time–memory step. A bare `Option` is not a violation when it documents a single failure class.

**Outside this scope:** NFS-DL's `solve_dl` returns `Ok(0)` for an empty frontier, and MOV's
k=2 discrete log runs through its exhaustive toy table; both belong to the NFS-DL contract design.

All follow-ups decided on 2026-11-09 landed on 2026-11-10; source, tests and git history are the
record. Residual findings, each outside what that work was allowed to change. Record new findings
here with an evidence label, a destination and a closure condition; delete an item when its fix
lands.

- **MOV k=2 engine has no budget (confirmed).** `mov_reduce` accepts single-limb k=2 fields, but
  the toy log engine's `find_irreducible_degree2` scans every r in 0..p, so a field near 2^64
  hangs instead of returning a typed outcome. Destination: the NFS-DL/MOV k=2 design. Closure: a
  budgeted or replaced k=2 engine with a typed exhaustion outcome.
- **Miller vertical line at x_Q = x_T (inferred, untested).** `vertical_at` returns 1 there, so for
  R ∈ ⟨G⟩ other than G nothing guarantees that a degenerate evaluation is detected. Every `Ok` is
  verified, so exposure is `Inconsistent` or `DegeneratePairing`, not a wrong logarithm. Closure:
  a derivation or a test over R ∈ ⟨G⟩.
- **Default CLI curve (confirmed, documented).** The default `secp-toy` Brent run cannot finish
  within its per-attempt cap (2²⁹ − 1 against ≈ 2.7·10⁹ expected steps; arithmetic, not run).
  The docs say so; decide whether the default should change.

### Distinguished-point walk step counts

**Confirmed observation, cause unknown:** on the 25-bit bench fixture `solve_dp_negmap` needs
3·10⁴–7·10⁴ steps against ≈ 2·10³–10⁴ for `solve_dp`, and 1 of 6 seeds timed out at 30 s; on a
31-bit group (p = 2_147_483_659, n = 2_147_490_413) it did not finish in 60–120 s for θ ∈ {1, 2,
3, 4, 6, 8}, and `solve_dp_batch` took 1.1·10⁶–6.3·10⁶ steps against ≈ 6·10⁴ expected. Step
counts cluster at multiples of the walker-restart threshold (**inferred:** walkers trapped in
fruitless cycles). Budgets bound every call and outcomes stay typed; no wrong answer is
implicated. The ECDLP benches' fixed seed masks the stall. Closure: a diagnosis, then a repair by
a published fruitless-cycle technique (read the source; recalled candidates: Bos–Kleinjung–Lenstra
2010, Bernstein–Lange–Schwabe 2011), before any negation-map or batch speedup is claimed or
measured.

## Oracles and assurance

### Implement comparisons, not just tool installation

**Confirmed:** meaningful `kat_pari_*` tests exist in `gnfs/tests/sqrt_pari_oracle_kat.rs`,
`factor_end_to_end_kat.rs` and `shared/numfield/tests/dedekind_pari_oracle_kat.rs`. Keep them opt-in and independently computed. The remaining inventory is:

- CADO placeholders in `gnfs/tests/{line_sieve,merge,lanczos,factor_end_to_end}_kat.rs`: a count
  compared with itself, `unimplemented!`, a hand-built matrix with no CADO invocation, and `todo!`.
- NFS-DL PARI placeholders in `gnfs/tests/dl_{individual_log,end_to_end,relation,ext}_kat.rs`:
  panic bodies or vacuous assertions, not algebraic comparisons.
- MOV, hyperelliptic and GHS placeholders in `rho/tests/{mov,hyperelliptic,ghs}_kat.rs`.
- Hard-coded, self-contained checks labelled as PARI tests in `shared/padic/tests/log_kat.rs` and
  `rho/tests/{ssa,semaev}_kat.rs`; distinguish oracle-derived known answers from live tool calls.
  Review their manual GP snippets rather than preserving unvalidated instructions.
- Ignored in-workspace cross-checks and genuinely expensive circuit tests in `shor/tests`:
  classify them by actual dependency/cost, not by the presence of `#[ignore]`.

Implement each meaningful comparison or remove its placeholder and retain the intended work here.
Do not replace an unimplemented check with an unconditional pass. Until resolved, descriptions must
say unimplemented, not merely tool unavailable. Fix the overstated oracle coverage in BENCHMARKS,
PEDAGOGY, AGENTS and development guidance. Document targeted invocations of implemented tests;
`cargo test -- --ignored` also encounters placeholders and expensive non-oracle tests.
An explicitly requested oracle run must fail on a missing tool, tool error or mismatch, not silently
pass or skip the comparison.

**Selected near-term delivery: one complete CADO GNFS comparison.** Design the independent evidence
and interchange boundary, then deliver a reproducible build, pinned inputs, artifact parsing and a
meaningful relation/dependency comparison together. An installation-only milestone does not close
the oracle gap. Use the existing checked GNFS boundary; do not freeze questionable NFS-DL behaviour
as expected output. Broader comparisons can follow their corresponding correctness work.

Early delivery wins on interoperability feedback and realistic-input evidence, but adds build/parser
complexity and competes for effort with direct repairs. Dedekind and NFS-DL repairs do not depend on
CADO: independent regressions and PARI remain appropriate for their algebraic contracts. Neither
installation nor implementation is authorised by this scheduling decision.

Relation-count tolerances and filtered matrix dimensions are heuristic diagnostics, not correctness
checks. A matching final factor is self-certifying and does not prove either pipeline's mechanism.
Prefer independently validated relations and checked dependencies, with cross-ingestion of a defined
relation/matrix format where useful. Define matched polynomial, sieve domain and large-prime policy
before comparing counts. Select comparisons by their properties, not a fixed quota; calibrate
statistical diagnostics separately from correctness assertions.

Then choose a reproducible upstream revision, build/cache layout, invocation and failure contract.
Verify the current upstream location and build requirements before choosing an immutable revision.
Lazy/on-demand building is a candidate, not a frozen design. No CI integration exists to attach an
opt-in flag to; choose that separately if CI is introduced.
The sidecar must remain optional and outside rDLP's computation and default test path.

### GNFS assurance extensions and local hardening

These are **confirmed limits of existing evidence**, not demonstrated failures of verified roots:

- Generated pipeline evidence is one monic cubic fixture succeeding at the first candidate. Retry,
  kernel combinations and CRT fallback have separate KATs, not generated-pipeline demonstrations.
- Broaden non-monic testing beyond degree two and single-prime recombination, and check agreement
  across methods in more domains, including reducible squarefree quotients where supported. Select
  independent expectations rather than treating regression-pinned counts as derivations.
- The PARI 45113 fixture uses the recorded dependency rather than regenerating the sieve. A genuine
  class-group-obstruction `NonSquareCertificate` fixture is still missing; PARI's Q(√−5) ideal
  facts alone do not exercise that crate outcome. Check incidental outcomes on their own contracts,
  not by requiring a mapped root when the map is intentionally unavailable.
- `NonSquareCertificate::new` is visible to sibling modules. The source-site test documents that it
  misses aliases and macros. Consider stronger construction confinement if justified; this is
  hardening against false negatives, not a discovered route to forging a verified root.
- Precision limits bound the certified minimum modulus, not every working allocation. Retain the
  documented lifting headroom and bounded sign-search accounting in new fixtures and measurements.

Keep fixture discovery bounded and recorded; never inject relations or an answer to make a generated
test succeed. Preserve the distinction between exact per-instance validation and all-input coverage.

## Mathematical exposition and documentation

### Correct false or unsupported claims

**Confirmed:** the following need correction, independently of CADO availability:

- MATHEMATICS's GNFS complexity derivation changes a leading constant by absorbing it into o(1),
  then substitutes another smoothness estimate. Re-derive the norm-size, Dickman and L-notation
  calculation and its optimisation together against a canonical reference; do not patch one number.
- README's GNFS introduction still says linear algebra produces a congruence of squares. It produces
  candidates; exact algebraic reconstruction remains necessary.
- `gnfs/src/sqrt/assembly.rs::AttemptOutcome::TrivialCongruence` says both gcds are trivial when
  x≡±y. For N=35 and x=y=5, gcd(x+y,N)=5. The detailed explanation correctly states the deliberate
  policy of not reporting this incidental divisor; align the short rustdoc with it.
- Residual non-monic prose calls b^d*f(a/b) the norm without the leading-coefficient factor.
  MATHEMATICS also uses `[O_K : Z[alpha]]` where alpha need not be integral. Use the normalised
  generator/order or a precise localisation; carry F_hom=c*Norm(a−b*alpha) consistently.
- NFS-DL equations, capability descriptions and full-implementation verdicts must follow the
  correctness adjudication above. Remove blanket claims that every described algorithm is already
  implemented and verified. A completed-history pass verdict is not current evidence.
- Four `gnfs` passages conflate p | disc(f) ("bad prime") with p | [O_K : ℤ[α]] or misstate
  consumers: `gnfs/docs/PEDAGOGY.md`'s bad-prime and `PolyPair::number_field()` passages, and
  `gnfs/src/sieve/factor_base.rs`'s module comment and residue-degree note. Some bad primes are
  maximal (x²−2 at 2); the factor base does its own root search, and `number_field()` serves only
  NFS-DL relation code. The same comment's "p | c can be invisible" describes a different flag:
  for non-monic f, `discriminant` returns lc(f)·disc(f), so `is_bad_prime` does flag p | c.
- `docs/development.md` incorrectly says there is no proptest layer. Describe the actual KAT-first
  mix and mathematical assertions, not stale counts or unconditional correctness guarantees.
- `shor/docs/PEDAGOGY.md` gives about 40 tests for `pointadd_kat.rs` (45 exist); several other
  orientation-doc counts have drifted. Prefer removing drifting counts to updating them.
- `docs/PEDAGOGY.md`'s GLV and negation-map passages (roughly the rho-optimisation tour) state √3
  and √6 walk-length reductions, some without the word "expected". **Confirmed:** no test asserts
  them; `glv_fewer_steps_than_negmap` checks only a 2× no-regression bound, and the Koblitz √(2m)
  saving is likewise unmeasured (the toy group has order 4, so orbits have at most two points).
  Label them "expected, not measured", or measure them on a fixture where the orbit is large.

Sweep related unchanged passages as well as edited lines. A clean planning-token grep does not
establish that the mathematics or capability claims are correct.

### Quantum-resource references

Rework MATHEMATICS's quantum resource discussion from primary sources with explicit evidence scope.
**Confirmed from the referenced abstracts:** logical versus physical qubits, per-point-addition
versus full-algorithm gate counts, and architecture/runtime assumptions must not be mixed.

- Babbush et al., *Securing Elliptic Curve Cryptocurrencies against Quantum Vulnerabilities*,
  <https://arxiv.org/abs/2603.28846v2>: distinguish full-Shor estimates from point-addition counts.
  State the withheld-circuit/zero-knowledge-proof evidence boundary and use the corrected version.
  The version of the separate Google-hosted PDF was not established.
- Cain et al., <https://arxiv.org/abs/2603.28627v1>: distinguish hardware-specific neutral-atom
  physical-qubit estimates from logical circuit resources; do not generalise the abstract's
  cryptographic-scale claim into a particular RSA resource claim without reading its derivation.
- Schrottenloher, <https://arxiv.org/abs/2606.02235v1>: disclosed point-addition circuits compare
  against Babbush et al.; explain that relationship rather than presenting independent estimates of
  identical scope. Distinguish secp256k1 results from generic prime-field variants.
- Re-check the Proos–Zalka baseline and Roetteler et al. comparison before publishing any
  attribution; the unsupported approximately 768-qubit claim is gone from MATHEMATICS and PEDAGOGY
  but SURVEY still carries it. Correct the Shor journal citation to SIAM J. Comput.

Put validated exposition and bibliography entries in MATHEMATICS. Competition harnesses and
engineering demonstrations do not automatically imply a new mathematical algorithm or implementation
requirement. Read the relevant full-text passages before publishing quantitative comparisons.

## Project-wide anneal

Annealing was introduced after this project began. A changed-lines-only check therefore cannot
establish that older surfaces meet the same standard. Audit application code, test code, project
configuration/infrastructure and all documentation, including untouched files and examples.

Check topic-native names, interfaces, references, obsolete scaffolding and construction residue.
Check claim accuracy and consistency between implementation, tests and exposition. The artifact should
stand independently of its construction history. Record actionable findings with a destination and
closure condition, not another permanent inventory. The known items below seed the audit; they are
not its scope ceiling.

Keep ownership clear: mathematical defects go to the correctness scopes, unnecessary mechanisms and
repeated explanation to minimalism, format/lint adoption to tooling, and obsolete-file removal to
document retirement. Do not duplicate those tasks here. Separate behaviour-preserving cleanup from
semantic repairs. Broad inspection catches legacy residue, but bounded edits limit review risk and
merge conflicts. Close the pass when findings are resolved or explicitly dispositioned, references
and claims agree, and affected checks pass; a clean keyword scan alone is not acceptance.

### Known navigation, layering and register findings

Keep valid code citations: human documentation is allowed to reference Rust paths and identifiers.
Retain useful citations and correct actual stale paths:
README and PEDAGOGY still name the removed `rho/src/field` and `rho/src/util` wrappers; a curve
rustdoc link uses `crate::field::Fp`; Shor citations name files now represented by `mod.rs` modules.
README also misplaces the root toolchain and the Pohlig–Hellman module.

Resolve these remaining documentation scopes without assuming a wholesale textbook rewrite:

- Replace to-be-written/appended references in GNFS and Shor PEDAGOGY. Add useful maths-first
  entry pointers in the numfield and GNFS guides.
- Normalise chapter-heading depth in MATHEMATICS, including nested headings. Choose an unambiguous
  section-reference convention: chapter-local GNFS §5 currently collides with whole-book chapter 5.
  Update all consumers together. Give the new checked-square-root sections stable navigable anchors
  without casually renumbering widely cited GNFS section ranges. A short reading path is optional.
- Track/phase labels are gone from README, PEDAGOGY and BENCHMARKS (MATHEMATICS's residual matches
  are bibliography keys). Replace the remainder in source/test comments, e.g. principle-N labels in
  gnfs linalg/sieve, `rho/src/semaev` and `rho/tests/ghs_kat.rs`. Preserve real technique
  groupings, algorithmic step numbers and quantum phase terminology. Replace undefined principle-N
  and contract labels with their actual invariant; do not erase useful interface explanations.
  **Confirmed by grep:** about 300 survivors remain in roughly 70 files outside the solver files
  already cleaned: `gnfs/` (including `gnfs/docs/PEDAGOGY.md`, the largest), `shared/padic`,
  `shared/gf2m`, `shared/numth` outside its guide, `rho/src/{semaev,hyperelliptic,binary_curve,
  pairing,ghs}` and `rho/src/ssa/lift.rs`, `shor/src/{sparse,qft,statevec,arith}`,
  `shor/tests/qft_kat.rs`, plus "C1"/ROADMAP mentions in `shared/numfield/docs/PEDAGOGY.md` and
  `gnfs/tests/{schirokauer_kat,special_q_kat}.rs`. The sweep regex must also cover the variants the
  original missed: `\bC-[A-Z][A-Za-z0-9]+\b`, `\bCat-?[ABCI]\b`, `[Pp]rinciple[- ]?[0-9]+`,
  `\bPhase [α-ω]`, `[α-ω]\.[0-9]+` and `\b(the|per the) plan\b`.
- One-time repository housekeeping, to delete once done: the first solver-repair commit carries a
  malformed `Co-authored-by` trailer (`@build <build@opencode.local>` instead of the identity format
  in `AGENTS-ROSTER.md`) and one later ledger commit has an unwrapped body, both fixable with
  `@git-editor`; a stale `git stash` entry ("WIP on master: d50ab68", a discarded experiment already
  superseded by committed work) can be dropped.
- Left for human triage by the last anneal pass (2026-11-10): "track"/family-letter vocabulary
  ("cross-track", "all five tracks" against six family letters, "α-Substrate" chapter titles,
  `docs/PEDAGOGY.md`'s "α"/"NF" labels, MATHEMATICS's "cross-track bridge"); "substrate
  over-specifies discipline" in two PEDAGOGY files; BENCHMARKS's "Tests added" headings;
  MATHEMATICS's Markup note on "the earlier undecided … choice"; README's "Two complementary
  artifacts" listing three and its stale `DynResidue` layout line; "frozen contract", "option A/B"
  and "Principle 1/3/4" in `rho/src/pairing/{ecext,mod}.rs`. Decide whether the track vocabulary
  is a reader-facing taxonomy (keep, define once) or construction residue (translate).
- Remove bare ROADMAP references as well as PLAN.md/NOTES.md/SURVEY.md references from durable
  prose and inline comments. The old token sweep missed extensionless names. Replace historical
  pass-verdict sections with honest scope statements where they contain useful substance.
- Move any needed measurement method/provenance into development documentation and remove its and
  `.cargo/config.toml`'s dependence on rolling notes. Prefer symbol links over brittle source-line
  numbers; remove drifting test counts and coverage percentages from orientation documents.
- Replace the textbook's construction-contract preamble with concise audience, notation, citation
  and depth conventions. History and completed chapter-pairing corrections need no replacement.

## Tooling and measurements

### Adopt formatting and lint gates deliberately

`cargo fmt-check` and `cargo lint` are declared but not required to pass workspace-wide. Measure
their current failures before sharding; the old backlog counts predate substantial source changes.
Choose a policy for numeric-table formatting and for justified mathematical identifiers before
remediation. Preserve topic-native notation rather than mechanically renaming mathematics to fit a
simplistic exemplar. Explicit exceptions need reasons; they are not permission to suppress defects.

Make formatting and lint adoption dedicated changes. Refresh affected links and documentation, then
adopt each gate only when green. Do not lower coverage or exclude files to compensate.

Re-measure published benchmarks whose Montgomery/RNG substrate or reconstruction method changed.
Record toolchain, hardware, input, command and method. Keep compared measurements in one controlled
run. BENCHMARKS must distinguish measurements from estimates and unsupported scale claims.
Do not invent new GNFS benchmark targets merely because old prose named nonexistent ones.

## Relentless minimalism across four surfaces

> Minimise the independent concepts, mechanisms and explanations that must be kept consistent,
> not merely the number of lines.

Give each explanation or policy one authoritative home. Repetition earns its place only through a
distinct reader need or an independent verification role. Audit all four surfaces, not documentation
alone:

| Surface | Look for |
|---|---|
| Application code | Unnecessary abstractions, wrappers, alternate paths and duplicate mechanisms |
| Test code | Repeated scaffolding or assertions adding no distinct evidence |
| Config/infrastructure | Repeated policy, redundant orchestration and unnecessary machinery |
| Documentation | Competing accounts, repeated explanations and construction narrative |

Explanation ownership is role-specific:

- Public rustdoc states interface contracts; inline comments explain necessary local rationale,
  hazards and non-obvious consequences, rather than narrating obvious code.
- Mathematical and pedagogical exposition owns derivations, conceptual motivation and worked
  examples. Preserve teaching value rather than replacing it with terse implementation summaries.
- Project docs own navigation, usage and operation; link to deeper accounts instead of copying them.
- Future commit messages explain why the change was needed, not a second copy of its implementation
  or design document. This work does not authorise rewriting existing git history.

Minimalism reduces maintenance and drift, but indiscriminate deletion or sharing can damage local
comprehension and test independence. Preserve independent expected-value derivations and oracles,
necessary runtime invariant checks, readable mathematical notation and useful reader context. Do not
replace independent verification with the implementation's helper merely to remove duplication, or
introduce a generic framework just to merge similar-looking code. No line-count reduction quota.

For each proposed simplification, identify what distinct value the removed material supplied, where
that value remains, and how correctness and explanation are preserved. Audit globally, then shard by
coherent topic and consumer boundary. Keep one live finding here rather than separate app/test/doc
copies of the same task. The candidates below are an initial set, not an exhaustive inventory.

### Initial maintenance and simplification candidates

- Review unused `proptest` dev-dependencies in gnfs, rho, shor and gf2m; keep the dependency where
  it is actually used in shared-field. Adding arbitrary property tests is not a dependency-cleanup
  requirement.
- Remove or justify the empty public `shared/bigint/src/mp.rs` placeholder; update rho's package
  description to reflect its actual algorithm scope.
- Consider consolidating the small prime-cursor, centering and exponent helpers in GNFS square-root
  code and repeated test/GP scaffolding. Weigh readability and independent checks against shared
  bugs and extra coupling. The deliberate primality-code duplicate avoids a crate edge; it is not
  automatically a deduplication defect.
- The `gp(script) -> Vec<String>` PARI/GP invocation helper (spawn `gp`, fail on a missing binary,
  stderr or a non-zero exit, return output lines) is copied in three oracle test files:
  `gnfs/tests/factor_end_to_end_kat.rs`, `gnfs/tests/sqrt_pari_oracle_kat.rs` and
  `shared/numfield/tests/dedekind_pari_oracle_kat.rs`. Sharing it across crates needs a small
  test-support crate (or a dev-dependency edge), so first decide whether that coupling is worth
  removing about thirty duplicated lines. Keep each file's own GP scripts and fixtures independent
  either way; dedupe only the process-spawning and error-handling wrapper.
- Reassess `shared/numth/src/smooth.rs::factor_base_up_to`'s large-bound sieve TODO against current
  primality implementation and measured need. Remove future-session prose regardless; optimisation
  is optional and its stated complexity rationale must be checked first.
- `gnfs/src/sieve/factor_base.rs` finds degree-one primes above p by its own root search, while
  numfield's `dedekind_factor` was designed for that role and has no production consumer. Adoption
  is blocked because `NumberField` is monic-only and GNFS polynomials need not be; whether the
  factor base or NFS-DL adopts the corrected report is for the NFS-DL contract design to decide.
- `iqft_on_qubits` is duplicated privately in `shor/src/ecdlp/mod.rs` and `shor/src/shor/mod.rs`
  beside the public `qft` module.
- Make independent-fixture provenance reproducible where tests cite unavailable Python/numpy/sympy
  scripts: retain the derivation or generator, or narrow the claim to what a reader can check.

## rust-template: conventions and exemplars in both directions

Treat `rust-template` as a deliberately simple baseline, not an exhaustive standard or a ceiling on
rDLP's design. Local project needs can expose missing reusable conventions. Conversely, complexity
and successful execution alone do not make an rDLP pattern exemplary.

Review both repositories and produce two kinds of output:

1. **Lightweight reusable conventions** suitable for template guidance: independently justified
   expectations; coverage versus correctness; meaningful opt-in external comparisons; checked input
   boundaries; semantic outcomes distinct from errors and exhausted budgets; precise statements of
   what private constructors and witness types actually guarantee. Include scale-appropriate
   separation of project backlog, current plan and transient notes, plus documentation registers
   and reference directions. Simple projects need not acquire every planning file. Human docs may
   cite code, but durable prose and inline docs must not depend on rolling files.
2. **Richer exemplars retained in rDLP**, with stable symbol references from the template if useful:
   explicit search limits and truncation reports; `SqrtContext`, `VerifiedRoot` and
   `NormalizedTarget::verify_candidate`; bounded fixture discovery; independent PARI comparisons.
   Explain when the design is warranted rather than transplanting a mathematical framework into a
   small application template.

For each candidate, verify that implementation, assertions and prose agree. Do not promote
placeholder tests, fallback-masked algorithms, source-scanning workarounds, duplicated scaffolding
or unresolved mathematical code as general conventions. A general rule, a linked exemplar, and an
intentionally local exception are different outcomes of this review.

**Confirmed concrete template issue:** `/home/justin/Source/rust-template/.cargo/config.toml` still
defines `audit = "deny check"`, which shadows cargo-audit when installed. Carry the alias-separation
correction to a separately approved template change. Its chosen 100% coverage floor need not become
rDLP's 80% floor, or vice versa. Lint/config differences need reasons, not automatic uniformity.
No template edits are authorised by this roadmap reconciliation.

This approach keeps the template teachable but puts more responsibility on example selection and
reference maintenance than copying every convention into one starter project.

## Retire SURVEY; review other rolling documents

The user selected eventual deletion of `docs/SURVEY.md`, not preservation as a historical findings
ledger. Its completed findings, stale inventories and superseded decisions need no migration.
Git history preserves them. The active scopes above replace it as the source of outstanding work.
SURVEY also quotes labels and comments since renamed in README, PEDAGOGY, BENCHMARKS and the
index-calculus module; those references now dangle, one more reason to transfer and delete it.

Before deletion, check and transfer only useful material that is not already present:

| Material | Destination |
|---|---|
| Still-live tasks and unresolved scope choices | This ROADMAP only |
| Current execution contracts and acceptance detail | A newly selected PLAN, only while needed |
| Test/oracle invocation and evidence conventions | `docs/development.md`; concise rules in AGENTS |
| Mathematical derivations, limitations and references | MATHEMATICS and the relevant PEDAGOGY |
| API invariants and implementation rationale | Rustdoc/comments and executable tests |
| General conventions or complex exemplars | Template review above, with separate approval |

The testing guide already carries the coverage floor and two-tier test pattern. Preserve the useful
test-length principle without old counts: split by coherent subject, not an arbitrary line limit.
Add mechanism-specific evidence and oracle-honesty guidance where absent. Keep real algorithmic
phase terminology while removing planning residue.

**Proposal requiring separate approval:** retire the completed GNFS PLAN (superseded in `PLAN.md`;
recover it with `git show dd3aca4:docs/PLAN.md`) after checking that its still-useful contracts
and derivations are represented in code, tests and focused human docs. Prune
or delete NOTES after the same transfer rather than maintaining a museum of old decisions. In
particular, the superseded norm-character argument must not remain authoritative: the corrected
account recognises retained norm-sign information and the limits of a sampled character set. Use the
current mathematical account, not the old contract body. Prefer transferring useful content over
maintaining obsolete contract bodies. This does not reopen the repaired GNFS implementation.

Check incoming links and bare document/contract references when removing files. Do not create a new
archive or duplicate backlog. These deletions and transfers require an implementation/docs session;
this reconciliation changes ROADMAP only.

## Algorithm-spectrum reassessment and deferred capabilities

**Confirmed assurance gap:** mathematical topic coverage does not establish algorithm fidelity.
Map theory, implementation, primary-path tests and independent oracles separately before claiming
there is only one remaining gap.
In particular, the inspected Shor code implements factoring and ECDLP, not the finite-field
multiplicative-group DLP pipeline. A theoretical generalisation is not its implementation.
Classify such omissions as scope decisions or extension candidates; do not silently implement them.

The hyperelliptic Jacobian DLP remains a **deferred candidate**, with its algorithmic regime still
to settle. Fixed-genus Gaudry cost O(q^(2−2/g)) is not subexponential in log(q). GHS over a small
field with growing genus calls for a different regime, such as Enge–Gaudry/ADH. The existing
hyperelliptic substrate is characteristic two and the GHS transfer fixture is genus one.
Higher-genus composition and reuse of mod-ell linear algebra need evidence and dependency review.
Recheck the stated GHS genus formula before using it to select a downstream solver; it was not
adjudicated in this review.
Settle the target regime and domain before selecting Gaudry as the implementation.

Other omitted families, including interval DLP, function-field/small-characteristic methods and
additional factoring methods, need classification against the intended scope, not automatic
enrolment. Claims such as GNFS strictly dominating QS at every size are not valid justifications for
omission. Circuit-resource papers and competition harnesses do not by themselves authorise new
algorithm families.

Complete projective/exceptional-ideal collection is also deferred: the current GNFS collector has
an explicit affine boundary and the checked matrix leaves exceptional algebraic columns
unconstrained. Exact root verification protects accepted results; collection completeness is not
claimed. Additional square-root families such as Montgomery/Nguyen are optional coverage, not
prerequisites to the correctness repairs above.

The Dedekind obstruction degree as log_p of the one-step enlargement index (Cohen §6.1) is not
exposed: the ring property was not confirmed against the printed theorem. Restore it only after
that confirmation.

## Maintenance rule

At a work boundary, update this register from actual source and evidence, not a session's title or
the completion of a checklist. Remove completed items once their useful content has landed. Do not
append campaign history, reopen accepted layout work by default, or convert an unverified hypothesis
into a frozen requirement. Keep future execution detail in the current PLAN and long-arc scope here.
