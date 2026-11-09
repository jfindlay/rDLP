# Shor's Algorithm: A Code-Tour Chapter

> **Code-first treatment.** This chapter is the code-tour sibling of `docs/MATHEMATICS.md` ch. 11.
> It assumes the reader knows the mathematics — period-finding, the QFT phase-estimation argument,
> the order-finding → factoring reduction, and the two-register hidden-subgroup construction for
> ECDLP — and focuses on *how the code realises* those ideas. For the mathematical foundations, see
> `docs/MATHEMATICS.md` ch. 11. For performance facts and qubit-budget tables, see
> `docs/BENCHMARKS.md` (the Shor simulator, factoring, and ECDLP sections).

---

## §1. At a Glance: Three Parts

The `shor` crate has three parts, each building on the last. The table below names the module
surface, the interfaces each part exposes, the canonical toy fixture, and the BENCHMARKS section.

| Part | Modules | Interfaces | Toy fixture | BENCHMARKS |
|------|---------|------------|-------------|------------|
| **Simulator** — state-vector circuits | `statevec`, `gates`, `sparse`, `measure`, `qft` | dense register, gate set, sparse register, Born-rule measurement, QFT | Bell state, GHZ-3, QFT on \|000⟩ | `docs/BENCHMARKS.md`, Shor simulator section |
| **Factoring** | `arith`, `shor` | modular exponentiation, order-finding, `factor` | 15 → {3,5}, 21 → {3,7}, 35 → {5,7}, 91 → {7,13} | `docs/BENCHMARKS.md`, Shor factoring section |
| **ECDLP** | `curve`, `ecc`, `ecdlp` | controlled point addition, `solve_ecdlp` | toy curve *y*² = *x*³ + 3 mod 7, *r* = 13; solves *k*·*G* = *Q* for *k* ∈ {3,5,7,9,12} | `docs/BENCHMARKS.md`, Shor ECDLP section |

**One gate set serves every circuit.** Every circuit in the factoring and ECDLP parts is assembled
entirely from the simulator's gate set. No gate is added to the simulator for the algorithms built
on it. This is the load-bearing substrate-reuse property: factoring and ECDLP are pure consumers of
the dense register.

---

## §2. The State-Vector Simulator

### §2.1 What it realises

The simulator is a classical state-vector quantum-circuit simulator: the substrate on which Shor's
algorithm runs. It is not a quantum computer; it is a classical program that faithfully simulates
the evolution of a quantum register by maintaining the full `2^n`-amplitude array. The simulator
is correct at any qubit count; only memory makes large *n* unreachable (see §6 on the resource
ceiling).

### §2.2 Module surface

**`statevec` — the dense register.**
`StateVec` holds `2^n` complex amplitudes as a `Vec<Complex<f64>>`. The basis-indexing convention
is **little-endian**: qubit 0 is the least-significant bit of the basis index. For an *n*-qubit
register, basis state |*q*_{*n*−1} … *q*₁ *q*₀⟩ maps to index
*i* = *q*₀ + 2·*q*₁ + 4·*q*₂ + … + 2^{*n*−1}·*q*_{*n*−1}.
This convention is fixed and consumed by every downstream module; a silent flip is a
wrong-answer bug. Constructors: `StateVec::zero(n)` (the |0…0⟩ state), `StateVec::basis(n, k)`
(the basis state |*k*⟩), `StateVec::from_amplitudes(n, amps)` (arbitrary normalised input). The
normalization invariant Σ|*aᵢ*|² = 1 is enforced at construction and preserved by all gates.

**`gates` — the universal gate set.**
Every gate is applied **in-place** by iterating the amplitude pairs it couples and updating them
with the gate's 2×2 unitary matrix. This is the standard O(2^*n*) state-vector method — no matrix
materialisation. The gate set:

- *Single-qubit:* X, Y, Z, H, S, T, phase(θ), arbitrary unitary U(*a*,*b*,*c*,*d*)
- *Two-qubit:* CNOT, controlled-phase(θ), SWAP
- *Multi-qubit:* Toffoli (CCX), multi-controlled-X (generalised Toffoli), multi-controlled unitary

The multi-controlled surface is more than the simulator's own KATs exercise (they use only single-
and two-qubit gates) but it is carried for the factoring part (modular exponentiation) and the
ECDLP part (point-addition circuit): the substrate is specified once, for all consumers.

**`sparse` — the sparse register.**
`SparseStateVec` stores only nonzero amplitudes in a `HashMap<usize, Complex<f64>>`. Gate
application iterates only the nonzero entries, applying the same unitary transformations as the
dense path. The sparse path is a no-regression mirror: a gate sequence on the sparse register
yields the same amplitudes as on the dense register (the sparse-dense agreement KAT). Sparsity
helps only while the state is sparse — a Hadamard on every qubit makes the state fully dense and
the sparse path matches the dense cost. This is a scale annotation: sparsity is state-dependent,
not a universal speedup. Dense↔sparse conversion: `SparseStateVec::from_dense` and
`SparseStateVec::to_dense`.

**`measure` — Born-rule measurement.**
`measure_all` samples a basis state *k* with probability |*aₖ*|² (the Born rule) and collapses
the register to |*k*⟩. `measure_qubit` measures a single qubit and renormalises the remaining
register conditioned on the outcome. The sampler uses `rand_chacha::ChaCha8Rng` seeded with a
`u64` seed, making KATs reproducible: the same seed + same state always produces the same
measurement outcome. `sample_counts` provides a multi-shot statistical helper.

**`qft` — the Quantum Fourier Transform.**
`qft` applies the standard Hadamard + controlled-phase ladder (O(*n*²) gates), built entirely from
the simulator's gate set. The **bit-reversal convention** is load-bearing: the standard QFT
circuit (Nielsen & Chuang) is designed for big-endian qubit ordering. In the little-endian
convention, the H + controlled-phase ladder maps |*j*⟩ → QFT|*j*_reversed⟩ instead of QFT|*j*⟩
— a silent wrong-answer bug for the period extraction in factoring and ECDLP. The fix is a single
**input bit-reversal step** (swap qubit *i* with qubit *n*−1−*i*) before the ladder; after this, the
ladder correctly computes QFT|*j*⟩ in natural little-endian order with no output bit-reversal
needed. `iqft` applies the inverse QFT (inverse controlled-phase + Hadamard ladder, reversed,
then output bit-reversal). The identity KAT `iqft(qft(sv)) == sv` is verified in `qft_kat.rs`.

### §2.3 Toy KAT

The canonical simulator KATs are in `shor/tests/statevec_kat.rs` and `shor/tests/qft_kat.rs`.

**Bell state.** Starting from |00⟩, apply H to qubit 0 then CNOT(0→1):

```rust
let mut sv = StateVec::zero(2);
gates::h(&mut sv, 0);
gates::cnot(&mut sv, 0, 1);
// Result: amplitudes[0] = 1/√2, amplitudes[3] = 1/√2, others = 0.
// (|00⟩ + |11⟩)/√2 in little-endian: index 0 = |00⟩, index 3 = |11⟩)
```

**QFT on |000⟩.** The QFT of the all-zero state is the uniform superposition:

```rust
let mut sv = StateVec::zero(3);
qft(&mut sv);
// All 8 amplitudes = 1/√8 (published: QFT|0⟩ = (1/√N)Σ|k⟩)
```

**QFT on |*j*⟩.** For *n* = 3, *j* = 1, the QFT amplitude at output index *k* is
(1/√8)·e^{2πi·*k*/8}, verified against Python-computed reference values.

### §2.4 Cross-reference

For the mathematical content — the QFT as a Fourier transform over ℤ/2^*n*ℤ, the phase-estimation
argument, and the Born-rule derivation — see `docs/MATHEMATICS.md` ch. 11. For the qubit-scaling
wall and the dense-vs-sparse performance comparison, see `docs/BENCHMARKS.md` (Shor simulator
section).

---

## §3. Shor's Factoring Algorithm

### §3.1 What it realises

This part runs the complete Shor-factoring algorithm end-to-end on the simulator: the reversible
modular-arithmetic quantum circuit builders, the order-finding circuit orchestration, the
classical continued-fraction period extraction, and the `factor(N)` driver. The algorithm factors
15 → {3,5}, 21 → {3,7}, 35 → {5,7}, and 91 → {7,13} — the canonical toy targets — using 8, 10, 12,
and 14 qubits respectively.

### §3.2 Module surface

**`arith` — reversible modular arithmetic.**
The `arith` module builds the controlled modular-exponentiation circuit
|*x*⟩|*y*⟩ → |*x*⟩|*y*·*aˣ* mod *N*⟩ from the simulator's gate set. The circuit hierarchy:

1. `controlled_add_mod(sv, ctrl, c, N, work_qubits)` — |*x*⟩ → |*x* + *c* mod *N*⟩ when ctrl = |1⟩.
2. `controlled_mult_mod(sv, ctrl, c, N, work_qubits)` — |*x*⟩ → |*c*·*x* mod *N*⟩ when ctrl = |1⟩.
3. `controlled_mod_exp(sv, a, layout)` — |*x*⟩|*y*⟩ → |*x*⟩|*y*·*aˣ* mod *N*⟩.

**Implementation approach: permutation synthesis.** Each modular-arithmetic operation is a
classical permutation on basis states. The permutation is synthesised as a product of
transpositions (selection-sort order), each implemented via a Gray-code path of single-bit-
difference transpositions, each implemented as a multi-controlled-X gate. This approach is:
ancilla-free (no scratch register beyond the work register), reversible (every permutation is
its own inverse composed with itself), and correct on superpositions (the gate sequence is
unitary). It is a demonstration-fidelity circuit, not gate-count-optimised (see §6.2).

**Register layout.** Qubits [0, *t*) are the exponent register (*t* bits);
qubits [*t*, *t*+*n*) are the work register (*n* bits, where *N* < 2^*n*). The work register
starts at |1⟩ for the mod-exp circuit (so *y*·*aˣ* mod *N* with *y* = 1 gives *aˣ* mod *N*).
`ModExpLayout::standard(N, exp_len)` constructs the standard layout. Total qubits: *t* + *n*.

| *N* | *n*_bits(*N*) | exp_len (*t*) | work (*n*) | Total qubits |
|-----|--------------|--------------|-----------|--------------|
| 15  | 4            | 4            | 4         | 8            |
| 21  | 5            | 5            | 5         | 10           |
| 35  | 6            | 6            | 6         | 12           |
| 91  | 7            | 7            | 7         | 14           |

Here $n_{\mathrm{bits}}(N) = \lceil \log_2 N \rceil$, the exponent register has
$t = n_{\mathrm{bits}}(N)$ qubits, and the work register has $n_{\mathrm{bits}}(N + 1)$ qubits
(`exp_len_for`,
`ModExpLayout::standard`; checked by `qubit_budgets_all_targets`). The textbook analysis takes
$2^t \geq N^2$, twice as many qubits, to bound the failure probability of continued-fraction
recovery; $t = \lceil \log_2 N \rceil$ is half that width, and its sufficiency here is empirical,
not proved. The only evidence is the fixtures: 15, 21, 35 and 91 each recover the order and factor
at a recorded seed (`factor_15_seed0`, `factor_21_seed0`, `factor_35_seed7`, `factor_91_seed1`).

**`shor` — order-finding and factoring driver** (`shor/src/shor/mod.rs`).
The `shor` module orchestrates the complete Shor-factoring algorithm:

1. **`run_order_finding_circuit(a, layout, seed)`** — allocates the register, puts the exponent
   register in uniform superposition (H on every exponent qubit), applies `controlled_mod_exp`,
   applies iQFT to the exponent register via `iqft_on_qubits`, measures with `measure_all_seeded`,
   and returns the measured phase numerator *s* (the lower `exp_len` bits of the measurement).

2. **`order_from_phase(s, t, a, N)`** — expands *s*/2^*t* as a continued fraction, takes the
   first convergent denominator *d* ≤ *N* with *a*^*d* ≡ 1 mod *N*, and reduces *d* to the
   multiplicative order *r*. A convergent denominator is a *multiple* of *r* (it can be a proper
   multiple, because the measured phase need not be in lowest terms), so the reduction divides out
   each prime *f* | *d* for as long as *a*^(*d*/*f*) ≡ 1 still holds. The result is the least
   *r* > 0 with *a*^*r* ≡ 1 mod *N*, never a multiple of it (`cf_order_from_phase_n15_a7_s32`:
   *s* = 32, *t* = 8, *a* = 7, *N* = 15 has first denominator 8 and order 4;
   `cf_order_from_phase_is_always_the_least_order`: every phase and every coprime base at the four
   register widths the targets use). It returns `None` when no convergent works, including
   *s* = 0 (`cf_order_from_phase_zero_phase`). This is the classical post-processing step that
   recovers the order *r* from the measured phase. See `docs/MATHEMATICS.md` ch. 11 for the
   continued-fraction approximation argument.

3. **`find_order(a, N, seed)`** — calls `run_order_finding_circuit` then `order_from_phase`; it is
   `None` for a measurement that carries no usable phase
   (`find_order_none_on_uninformative_measurement`) and otherwise the least order
   (`find_order_returns_valid_order`).

4. **`factor(N, seed)`** — the end-to-end driver, returning
   `Result<ShorFactorization, ShorFactorError>`. It classifies *N* classically before any circuit
   runs, in this order:

   | Check | Result | Exercised by |
   |-------|--------|--------------|
   | *N* ≤ 2 | `NotComposite` (0 and 1 have no factorisation; 2 is prime) | `factor_prime_is_not_composite` |
   | *N* even | `Ok`, *p* = 2, evidence `Even` | `factor_even_short_circuit` |
   | *N* = *p*^*e*, *e* ≥ 2 | `Ok`, *p* = the prime, evidence `PrimePower { prime, exponent }` | `factor_prime_power_short_circuit` |
   | circuit needs more than `MAX_QUBITS` = 25 qubits (*N* ≥ 4097) | `Unsupported { qubits }`, with `qubits` = *n*_bits(*N*) + *n*_bits(*N* + 1) | `factor_beyond_the_simulator_is_unsupported` |
   | *N* prime (trial division) | `NotComposite` | `factor_prime_is_not_composite` |

   Because the qubit ceiling is checked before the primality test, a large odd prime is
   `Unsupported`, not `NotComposite`: an odd *N* beyond the ceiling that is not a prime power is
   not classified further (`factor_beyond_the_simulator_is_unsupported`, which includes 1000003).
   Only then does the quantum order-finding loop run: for each of the first `MAX_BASES` = 20 bases
   *a* coprime to *N* (in increasing order), up to `SEEDS_PER_BASE` = 8 circuit runs at seed
   `seed + 1000·(base index) + (run index)`; each measurement *s* gives *r* by `order_from_phase`,
   and if *r* is even and *a*^(*r*/2) ≢ −1 mod *N*, then gcd(*a*^(*r*/2) ± 1, *N*) gives a
   factor. The first one found is returned with evidence `OrderFinding { a, seed, s, r }`, which
   a test rechecks from its fields: *r* is the least order of *a*, is even, and *p* is one of the
   two gcds (`assert_factor`, run by `factor_15_seed0`, `factor_21_seed0`, `factor_35_seed7`,
   `factor_91_seed1`). The returned *p*·*q* = *N* is checked by division. The bases are drawn
   coprime to *N*, so no base shares a factor with *N* and there is no common-factor outcome:
   taking `gcd(a, N)` as a shortcut would turn the base scan into trial division, which is not
   the algorithm. The loop is bounded: spending the budget is
   `AttemptsExhausted { attempts }`, with `attempts` the circuit runs actually made (fewer than
   160 for small *N*, which has fewer than 20 coprime bases), a resource outcome that says
   nothing about *N* (the unit test `spent_budget_is_attempts_exhausted_with_the_runs_made` in
   `shor/src/shor/mod.rs` exhausts a budget of one run, succeeds at a budget of two, and makes
   no run for a budget of zero). The same seed gives the same result, evidence included
   (`factor_reproducible`).

### §3.3 Toy KAT

The canonical factoring KATs are in `shor/tests/modexp_kat.rs` and `shor/tests/factor_kat.rs`.

**Modular exponentiation.** On basis state |*x*⟩|1⟩ with *a* = 2, *N* = 15, *x* = 3:

```rust
let layout = ModExpLayout::standard(15, 4); // 8 qubits total
let mut sv = make_modexp_state(&layout, 3, 1); // |exp=3⟩|work=1⟩
controlled_mod_exp(&mut sv, 2, &layout);
// work register reads 2^3 mod 15 = 8. Published: ord₂(15) = 4.
```

**End-to-end factoring.** The four canonical targets, each at a recorded seed:

```rust
factor(15, 0)  // → Ok(p: 3, q: 5,  evidence: OrderFinding { a: 2, seed: 0, s: 4,  r: 4 })
factor(21, 0)  // → Ok(p: 3, q: 7,  evidence: OrderFinding { a: 2, seed: 0, s: 27, r: 6 })
factor(35, 7)  // → Ok(p: 5, q: 7,  evidence: OrderFinding { a: 2, seed: 7, s: 59, r: 12 })
factor(91, 1)  // → Ok(p: 7, q: 13, evidence: OrderFinding { a: 2, seed: 1, s: 53, r: 12 })
```

These use 8, 10, 12 and 14 qubits. In each case the evidence shows that the first circuit run for
base 2 recovered the order (`factor_15_seed0`, `factor_21_seed0`, `factor_35_seed7`,
`factor_91_seed1`, which recheck the evidence from its fields). The seeds are ones that succeed
on the first run — a fixture choice, not a success rate: at seed 8, *N* = 35 does not recover an
order on its first run, and the ignored `find_seed_for_91` is the seed finder for 91. *N* = 91 is
the largest fixture; the simulator ceiling is `MAX_QUBITS` = 25, which admits *N* ≤ 4096, so 91 is
not at the ceiling.

### §3.4 Cross-reference

For the mathematical content — the order-finding → factoring reduction, the continued-fraction
approximation theorem, and the probability analysis — see `docs/MATHEMATICS.md` ch. 11. For the
qubit-budget table and the resource-scale annotation, see `docs/BENCHMARKS.md` (Shor factoring
section).

---

## §4. Shor's ECDLP Algorithm

### §4.1 What it realises

This part runs Shor's ECDLP algorithm end-to-end on the state-vector simulator: the toy elliptic
curve, the reversible controlled point-addition circuit built by permutation synthesis, the
two-register period-finding circuit, the classical dual-lattice discrete-log extraction, and the
`solve_ecdlp` driver. The algorithm solves the ECDLP *Q* = *k*·*G* on a 4-bit toy curve
(*r* = 13) using 17 qubits. Every success is verified by *k*·*G* = *Q* and carries the measured
pair that produced it; there is no classical search over *k*.

### §4.2 Module surface

**`curve` — the toy elliptic curve.**
The `curve` module defines the short-Weierstrass curve *y*² = *x*³ + 3 mod 7 with plain `u64`
arithmetic. Parameters: field prime *p* = 7, coefficients *a* = 0, *b* = 3, generator *G* = (1, 2),
prime group order *r* = 13. The 13 group elements (including ∞) are:

```text
0·G = ∞,      1·G = (1,2),  2·G = (6,3),  3·G = (2,2),
4·G = (4,5),  5·G = (3,3),  6·G = (5,3),  7·G = (5,4),
8·G = (3,4),  9·G = (4,2),  10·G = (2,5), 11·G = (6,4),
12·G = (1,5), 13·G = ∞
```

The classical group law (`add`, `double`, `negate`, `scalar_mul`) is the reference for the
circuit's KATs. The identity ∞ is encoded as the coordinate pair (0, 0) — a reserved sentinel
that does not appear on the curve (x = 0 gives y² = 3 mod 7, which has no solution).

`PointAddLayout` documents the qubit layout consumed by the point-addition circuit:
- qubit 0: control qubit
- qubits [1, 4): x-coordinate (3 bits, ⌈log₂ 7⌉ = 3)
- qubits [4, 7): y-coordinate (3 bits)
- qubits [7, 10): λ scratch register (3 bits, allocated but unused at runtime — see below)

Total: 10 qubits for the point-addition sub-circuit.

**`ecc` — the reversible controlled point-addition circuit.**
`controlled_point_add(sv, cg, layout)` applies |*P*⟩ → |*P* + *cG*⟩ when the control qubit is
|1⟩, for a classically-fixed point *cG*.

**Construction choice: permutation synthesis.** Since *cG* is a classical constant (precomputed),
this circuit performs *constant-point addition* — adding a fixed classical point to a quantum
register holding the running point *P*. The map (*x*_P, *y*_P) → (*x*_{P+cG}, *y*_{P+cG}) is
computed classically for all 13 group elements (including ∞) inside `build_point_add_permutation`,
then synthesised as a product of transpositions on the combined (*x*, *y*) register using the same
`apply_controlled_permutation` primitive from `arith`. The λ scratch register is allocated in the
layout for documentation and for the two-register circuit's layout, but is **not consumed at
runtime** — the permutation synthesis is ancilla-free. This mirrors the construction used for
modular arithmetic in the factoring part.

All exceptional cases (P = ∞, P = *cG*, P = −*cG*) are handled automatically by the permutation
synthesis — the classical group law is computed for every input, so no special-case branching is
needed in the circuit.

`controlled_point_add_inv(sv, cg, layout)` applies the inverse (adds −*cG*), used for
reversibility verification.

**`ecdlp` — the two-register ECDLP circuit and driver.**
The `ecdlp` module implements Shor's ECDLP algorithm as a two-register hidden-subgroup problem.
The register layout (17 qubits total):

```text
qubits [0, 4)   — a-register (4 bits, t = ⌈log₂ 13⌉ = 4)
qubits [4, 8)   — b-register (4 bits)
qubits [8, 11)  — x-coordinate of running point (3 bits)
qubits [11, 14) — y-coordinate of running point (3 bits)
qubits [14, 17) — λ scratch register (3 bits, allocated but unused)
Total: 17 qubits
```

**`run_period_finding_circuit(G, Q, seed)`** — the two-register circuit:

1. Allocate 17 qubits; initialise point register to |∞⟩ = |(0,0)⟩.
2. Apply H to every qubit of the a-register and b-register (uniform superposition).
3. For each bit *j* of the a-register: apply `controlled_point_add` with *cG* = 2^*j*·*G*,
   controlled on a-register qubit *j*.
4. For each bit *j* of the b-register: apply `controlled_point_add` with *cG* = 2^*j*·*Q*,
   controlled on b-register qubit *j*.
5. Apply iQFT to the a-register, then iQFT to the b-register.
6. Measure all qubits with `measure_all_seeded`.
7. Extract the pair (*c*, *d*) from the exponent register bits.

Because the circuit measures every qubit, the exponent pair (*c*, *d*) is read from the exponent
bits of the one measured basis state; the point and scratch bits are discarded.

**`extract_k_from_measurement(c, d, r)`** — the dual-lattice extraction. With *M* = 2^`EXP_BITS`
= 16 (not a multiple of *r* = 13), measured pairs concentrate near (*j*/*r*, *l*/*r*) with
*l* ≡ *k*·*j* (mod *r*), whichever sign convention the transform uses. The rule rounds each
measured value to the nearest multiple of 1/*r* (ties up): *j* = round(*r*·*c*/*M*) mod *r*,
*l* = round(*r*·*d*/*M*) mod *r*, and *k* = *l*·*j*⁻¹ mod *r*. It returns exactly one candidate
per measurement, or `None` iff *j* ≡ 0 (an uninformative measurement), and it panics when *c* or
*d* is ≥ *M*, which is invalid input rather than an uninformative measurement. The candidate set
is never widened to neighbouring roundings, and no continued-fraction step is involved. See
`docs/MATHEMATICS.md` §11.4 for the derivation.

**`exponent_outcome_distribution(G, Q)` and `primary_success_probability(G, Q)`** — the exact
distribution *P*(*c*, *d*) of the exponent pair, computed from the circuit's amplitudes, and its
total over the outcomes whose extraction returns the true logarithm. They are used by the tests
as the exact reference; `primary_success_probability` scans 0..*r* classically to find the
logarithm, so it is an oracle for tests, not part of the solver.

**`solve_ecdlp(G, Q, seed)` and `solve_ecdlp_with_retries(G, Q, seed, max_attempts)`** — the
end-to-end driver. `solve_ecdlp` uses `DEFAULT_MAX_ATTEMPTS` (8). The driver:

1. returns `EcdlpError::InvalidGenerator` for *G* = ∞;
2. handles the trivial cases (*Q* = ∞ → *k* = 0; *Q* = *G* → *k* = 1) with
   `EcdlpEvidence::Trivial`;
3. for each attempt *n*, runs `run_period_finding_circuit` at seed `seed.wrapping_add(n)` and
   takes the single candidate from `extract_k_from_measurement`;
4. accepts the candidate only if *k*·*G* = *Q*, returning `EcdlpSolution { k, evidence }` with
   `EcdlpEvidence::Measurement { seed, attempt, c, d }`, so a test can recheck the measured pair
   with the extraction rule;
5. otherwise measures again; after `max_attempts` it returns
   `EcdlpError::AttemptsExhausted { attempts }`.

There is no classical fallback: no exponent is enumerated. Repeating the same randomized
measurement with fresh randomness is the algorithm itself, and exhausting the budget is a
resource outcome that says nothing about whether a logarithm exists. A *Q* outside the group
generated by *G* fails verification on every attempt and also exhausts the budget.

### §4.3 Toy KAT

The canonical ECDLP KATs are in `shor/tests/pointadd_kat.rs` and `shor/tests/ecdlp_kat.rs`.

**Classical group law.** The full scalar-multiple table is verified against hand-computed values:

```rust
assert_eq!(scalar_mul(4, G), Point::Affine { x: 4, y: 5 }); // 4·G = (4,5)
assert_eq!(scalar_mul(13, G), Point::Infinity);               // order check: 13·G = ∞
```

**Point-addition circuit.** For each group element *P* and each *cG*, the circuit gives *P* + *cG*:

```rust
let layout = PointAddLayout::standard(); // 10 qubits
let mut sv = make_point_state(1, G, &layout); // |ctrl=1⟩|P=G⟩
controlled_point_add(&mut sv, G, &layout);
// Point register reads 2·G = (6,3). Exceptional cases (P=∞, P=−cG) also verified.
```

**Extraction and distribution.** `ecdlp_kat.rs` checks, without running the circuit, that the
nearest-integer representatives of exact lattice points recover *k* for every *k* ∈ 0..13 and
every *j* ∈ 1..13 (`extract_exact_peak_pairs_recover_k`), that only *c* = 0 yields `None` when
*M* = 16 and *r* = 13 (`extract_j_zero_returns_none`), and that out-of-range *c* or *d* panics
(`extract_rejects_c_out_of_range`, `extract_rejects_d_out_of_range`). It also checks that the
circuit's exact distribution equals the classical Fourier-sum formula to within 1e-9 and sums
to 1, for *k* = 3 and *k* = 7 on the default path (`distribution_matches_fourier_formula_k3`,
`distribution_matches_fourier_formula_k7`) and for every *k* ∈ 2..13 in an `#[ignore]`-gated test
(`distribution_matches_fourier_formula_all_k`). The one-candidate rule has positive success
probability for those *k* (`primary_success_probability_positive_k3_k7`, and
`primary_success_probability_positive_all_k` ignored for cost).

**End-to-end ECDLP.** The five canonical instances, each at seed 0:

```rust
solve_ecdlp(G, scalar_mul(3, G), 0)   // Ok(EcdlpSolution { k: 3, evidence: Measurement { .. } })
solve_ecdlp(G, scalar_mul(5, G), 0)   // k = 5
solve_ecdlp(G, scalar_mul(7, G), 0)   // k = 7
solve_ecdlp(G, scalar_mul(9, G), 0)   // k = 9
solve_ecdlp(G, scalar_mul(12, G), 0)  // k = 12
```

(`ecdlp_k3_seed0`, `ecdlp_k5_seed0`, `ecdlp_k7_seed0`, `ecdlp_k9_seed0`, `ecdlp_k12_seed0`.)
Each test asserts the reported *k*, that it succeeds on its first attempt with the recorded pair
(*c*, *d*) = (10, 13), (12, 13), (9, 12), (8, 13), (4, 12) for *k* = 3, 5, 7, 9, 12, and that the
`Measurement` evidence is itself consistent: the pair is in range and
`extract_k_from_measurement(c, d, r) == Some(k)`, so the test rechecks the mechanism and not
only the self-certifying *k*·*G* = *Q*. Seed 0 was kept for these fixtures because it extracts
correctly on the first attempt; that is a fixture choice, not a success rate. The other solver
outcomes have their own tests: `ecdlp_trivial_q_infinity` and `ecdlp_trivial_q_equals_g`
(`Trivial` evidence), `ecdlp_invalid_generator`, `ecdlp_budget_exhausted_on_wrong_candidate` (a
measurement whose single candidate is wrong, with a budget of one attempt, is
`AttemptsExhausted`), `ecdlp_retry_uses_fresh_seed` (a failed attempt is followed by an attempt
at the next seed), and `ecdlp_zero_budget_is_exhausted`.

**Success probability (estimate).** `primary_success_probability` computes the chance that one
measurement yields the true *k*. No test asserts its value. A single run of the ignored test
printed values around 0.64–0.66 for *k* = 2..12, so the probability that all eight default
attempts fail is roughly 0.36⁸ ≈ 3·10⁻⁴. Both numbers are estimates from that printout, not
tested facts.

### §4.4 Cross-reference

For the mathematical content — the two-register hidden-subgroup construction, the dual relation
and the finite-precision extraction rule, and the reduction from ECDLP to period-finding — see
`docs/MATHEMATICS.md` ch. 11. For the qubit-budget table and the resource-scale annotation, see
`docs/BENCHMARKS.md` (Shor ECDLP section).

---

## §5. How the Parts Reuse One Another

The crate has a clean substrate-reuse story. The diagram below shows the dependency chain:

```text
dense register + gate set (`statevec`, `gates`)
  ↓ gate set consumed by
sparse register, QFT, measurement (`sparse`, `qft`, `measure`)
  ↓ QFT + measure consumed by
modular exponentiation (`arith`)  ←── permutation primitives re-consumed by ──→  point addition (`ecc`)
  ↓ consumed by                                                                    ↓ consumed by
order-finding + factoring (`shor`)                                          ECDLP circuit + driver (`ecdlp`)
```

**The dense register is the foundation.** The dense register and gate set are the only substrate.
Every circuit in the factoring and ECDLP parts is assembled from `gates::h`,
`gates::controlled_phase`, `gates::swap`, `gates::multi_controlled_x`, and
`gates::multi_controlled_unitary`. No gate is added for the algorithms built on top.

**The QFT and measurement are reused by factoring and ECDLP.** The `iqft_on_qubits` helper (a
partial iQFT on a sub-register) appears in both `shor::shor` (for the exponent register in
factoring) and `shor::ecdlp` (for the a- and b-registers in ECDLP). The `measure_all_seeded`
function is the measurement primitive for both. The QFT bit-reversal convention is consumed
silently by both — a wrong convention would produce wrong-answer bugs in both.

**The `arith` primitives are re-consumed by point addition.** The
`apply_controlled_permutation` function from `arith` (the permutation-synthesis primitive) is
imported directly by `ecc::controlled_point_add`. Carrying the full permutation-synthesis
machinery for modular arithmetic pays off in the ECDLP part: the point-addition circuit reuses
the same primitive without any new gate or new abstraction.

**Point addition feeds the ECDLP circuit.** The `controlled_point_add` function is the only
circuit primitive consumed by `ecdlp::run_period_finding_circuit`. The ECDLP circuit is a loop
over the a- and b-register bits, each iteration calling `controlled_point_add` with a
classically-precomputed power of *G* or *Q*.

**The gate set is never extended for the algorithms.** This is the observable consequence of the
chain above: the `shor` crate's `gates` module is complete before the factoring and ECDLP
circuits are built, and they are pure consumers of it.

---

## §6. What the Implementation Claims

### §6.1 Algorithmic content

**Shor's factoring algorithm** is implemented head-on: the full quantum order-finding circuit
(exponent-register superposition → controlled-mod-exp → iQFT → measure), the classical
continued-fraction period extraction, and the even-order factor extraction via
gcd(*a*^(*r*/2) ± 1, *N*). The four canonical targets (15, 21, 35, 91) are all factored
correctly at recorded seeds, and every returned factorisation is checked by division. No step is
elided or approximated, with one stated difference from the textbook: the exponent register has
half the textbook width (§3.2), so the continued-fraction recovery is not covered by the textbook
failure bound, and its sufficiency is empirical.

**Shor's ECDLP algorithm** is implemented head-on: the full two-register period-finding circuit
(*a*·*G* + *b*·*Q* in superposition → iQFT on both registers → measure every qubit), the classical
dual-lattice extraction *k* = *l*·*j*⁻¹ (mod *r*) with *j* and *l* the nearest multiples of 1/*r*
to the measured *c*/*M* and *d*/*M*, and the verification *k*·*G* = *Q*, repeated with a fresh seed
for up to a fixed number of attempts. The five canonical instances (*k* = 3, 5, 7, 9, 12) are
solved at seed 0 (§4.3). The extraction rule is analysed only at *M* = 16, *r* = 13; nothing here
claims it recovers *k* for general *r*.

### §6.2 No engineering optimisation

**The point-addition circuit uses the straightforward permutation-synthesis construction, not a
qubit-optimised one.** The group law is computed classically inside `build_point_add_permutation`
for all 13 group elements; the resulting permutation is synthesised as a product of transpositions
via Gray-code paths of multi-controlled-X gates. This is a demonstration-fidelity circuit: it is
correct, reversible, and ancilla-free, but it is not gate-count-optimised. A production
implementation (e.g., Roetteler et al. 2017) would use an explicit reversible affine formula with
a dedicated λ register and Brent–Kung adder trees. The λ register is allocated in `PointAddLayout`
for documentation and future-compatibility, but is not consumed at runtime by the permutation
synthesis.

Similarly, the modular-exponentiation circuit uses direct permutation synthesis (not the
Vedral–Barenco–Ekert adder-based construction). The gate count is O(*N*·*t*) transpositions per
mod-exp stage, not the O(*n*²·*t*) of an adder-based circuit — but the permutation synthesis is
also not the asymptotically optimal construction. The choice is deliberate: the code exhibits the
algorithm's logic, not its engineering optimisation.

### §6.3 The ~25-qubit resource ceiling

**The ~25-qubit ceiling is a resource-scale wall, not a mathematical one.** The dense register
holds 2^*n* complex amplitudes. At *n* = 25 that is 2^25 ≈ 33 M entries (≈ 512 MiB of f64
pairs). This is the practical ceiling on a laptop. The factoring driver enforces it as the constant
`MAX_QUBITS` = 25 (*n*_bits(*N*) + *n*_bits(*N* + 1) qubits, which admits *N* ≤ 4096) and reports a
larger circuit as `Unsupported` before running anything
(`factor_beyond_the_simulator_is_unsupported`).

The ceiling is annotated consistently across all three parts:

- **Simulator** (`docs/BENCHMARKS.md`, Shor simulator section, `### Science↔engineering note`):
  the simulator demonstrates Shor's mathematics correctly at toy scale; the QFT produces the
  correct Fourier amplitudes and measurement samples from the correct Born-rule distribution. The
  ceiling is purely engineering.

- **Factoring** (`docs/BENCHMARKS.md`, Shor factoring section, `### Science↔engineering note`):
  the order-finding circuit demonstrates Shor's factoring mathematics correctly. This project's
  layout uses 2·*n*_bits(*N*) qubits (the qubit-budget table above), which is about 4100 qubits
  for a 2048-bit *N* — arithmetic on the layout, not a published estimate — requiring a
  2^4100-entry array, far beyond any classical computer.

- **ECDLP** (`docs/BENCHMARKS.md`, Shor ECDLP section, `### Science↔engineering note`): the
  two-register circuit demonstrates Shor's ECDLP mathematics correctly. This project's layout uses
  2·*t* + 3·⌈log₂ *p*⌉ qubits, which is 2·256 + 3·256 = 1280 for a 256-bit curve — again arithmetic
  on the layout, not a published estimate — requiring a 2^1280-entry array.

The simulator demonstrates Shor's mathematics correctly at toy scale; it does not claim quantum
speedup, which requires real quantum hardware out of scope by construction.

---

## §7. KAT Summary

The crate has six KAT files in `shor/tests/`. The table below summarises the test corpus.

| File | Tests | What is verified |
|------|-------|-----------------|
| `statevec_kat.rs` | 32 | Dense register: unitarity (X, Y, Z, H, S, T, CNOT, controlled-phase, SWAP, Toffoli, multi-controlled-X); Bell state (|00⟩ + |11⟩)/√2; GHZ states (3, 4, 5 qubits); normalization invariant; gate identities (HH=I, XX=I, S²=Z, T²=S) |
| `qft_kat.rs` | 36 | QFT on |0…0⟩ (*n* = 1,2,3,4): uniform superposition; QFT on basis states (*n* = 3,4): published Fourier amplitudes; QFT∘iQFT = identity; measurement distribution (Born-rule frequencies, seeded sampler); sparse-dense agreement (all gate types + round-trip conversion + the dense-state scale annotation) |
| `modexp_kat.rs` | 41 | Permutation correctness (controlled-add-mod, controlled-mult-mod, full sweeps); modular exponentiation correctness (2^3 mod 15 = 8, 7^2 mod 15 = 4, etc.); reversibility (forward + inverse = identity); ancilla-clean (full state restored); control-off no-op |
| `factor_kat.rs` | ~30 | Order KATs (ord₂(15)=4, ord₇(15)=4, ord₂(21)=6, ord₂(35)=12, ord₂(91)=12); continued-fraction KATs (known phase → known order, the least order for every phase and coprime base, deterministic classical); end-to-end factoring (15→3×5, 21→3×7, 35→5×7, 91→7×13) with the `OrderFinding` evidence rechecked; the `Even`, `PrimePower`, `NotComposite` and `Unsupported` outcomes; qubit-budget verification; `#[ignore]`-gated seed-finder for N=91 |
| `pointadd_kat.rs` | ~40 | Classical group-law correctness (all 13 scalar multiples, negation, commutativity, associativity, order check 13·G=∞); permutation correctness (full sweeps over all *P* for *cG* = G, 2G, 6G; exceptional cases ∞+G, G+(−G)); reversibility (full sweeps); scratch-clean (λ register = 0 after circuit); control-off no-op |
| `ecdlp_kat.rs` | 34 (7 `#[ignore]`) | Order/group-law KATs; dual-lattice extraction KATs (exact peak pairs → *k* for every *k* ∈ 0..13 and *j* ∈ 1..13, `None` iff *j* ≡ 0, out-of-range inputs rejected); exact outcome distribution equal to the Fourier-sum formula (*k* = 3, 7 on the default path); positive one-candidate success probability; end-to-end ECDLP (*k* = 3,5,7,9,12 at seed 0 with `Measurement` evidence; trivial cases; invalid generator; exhausted and zero budgets; retry at a fresh seed); qubit budget (17 qubits, within ceiling). `#[ignore]`-gated: `distribution_matches_fourier_formula_all_k`, `primary_success_probability_positive_all_k`, `ecdlp_all_k_values_within_default_budget`, `ecdlp_reproducible` (cost: each runs full circuit simulations), and three `rho_crosscheck_*` tests |

**The `#[ignore]`-gated tests.** Two categories:

1. *Slow full-sweep tests* (`distribution_matches_fourier_formula_all_k`,
   `primary_success_probability_positive_all_k`, `ecdlp_all_k_values_within_default_budget`,
   `ecdlp_reproducible`): each runs one or more full circuit simulations per *k* (roughly 10 s
   per simulation in debug mode, per the ignore reasons in the test file). The targeted fixture
   tests (*k* = 3, 5, 7, 9, 12; distribution for *k* = 3, 7) cover the same ground in the normal
   test suite.

2. *`rho` cross-check tests* (`rho_crosscheck_k3`, `rho_crosscheck_k7`,
   `rho_crosscheck_agrees_with_shor_k5`): compare the Shor solver against
   `rho::ecdlp::solve_brent` on the equivalent curve. These need no external tool: `rho` is a
   workspace dev-dependency and the comparison is in-workspace. They are kept off the default
   path for cost — two are cheap on their own, and the third runs a full Shor solve — and are
   not part of the green path.

---

## §8. Cross-References and Further Reading

**The maths-first sibling.** `docs/MATHEMATICS.md` ch. 11 is the mathematical counterpart to
this code-tour. It covers the QFT as a Fourier transform over ℤ/2^*n*ℤ, the phase-estimation
argument, the order-finding → factoring reduction (the period-finding payoff proof), the
two-register hidden-subgroup construction for ECDLP with its dual relation and finite-precision
extraction rule, and the PQC context (Shor's algorithm as the quantum threat to RSA and ECC).

**BENCHMARKS.** The performance facts cited in this document — qubit-budget tables, wall-clock
timings, the dense-vs-sparse comparison, and the `### Science↔engineering note` annotations —
are in the Shor sections of `docs/BENCHMARKS.md`. This code-tour cites those sections; it does
not extend them.

**Primary literature.**

- Shor, P.W. (1994). "Algorithms for quantum computation: discrete logarithms and factoring."
  *Proceedings of the 35th Annual Symposium on Foundations of Computer Science (FOCS 1994).*
  The original paper presenting both the factoring and discrete-logarithm algorithms.

- Proos, J., Zalka, C. (2003). "Shor's discrete logarithm quantum algorithm for elliptic curves."
  *Quantum Information and Computation* 3(4), 317–344.
  The two-register ECDLP construction that the `ecdlp` module implements.

- Nielsen, M.A., Chuang, I.L. (2000). *Quantum Computation and Quantum Information.*
  Cambridge University Press. §5.3 (order-finding and factoring), §5.4 (discrete logarithms).
  The standard textbook reference for the QFT, phase estimation, and Shor's algorithm.

- Roetteler, M., Naehrig, M., Svore, K.M., Lauter, K. (2017). "Quantum resource estimates for
  computing elliptic curve discrete logarithms." *ASIACRYPT 2017.*
  The production-scale ECDLP circuit design (the qubit-optimised construction that the `ecdlp`
  module deliberately does not implement, see §6.2).
