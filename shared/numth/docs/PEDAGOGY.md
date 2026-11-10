# The α-Substrate: A Code-Tour Chapter

> **Maths-first treatment.** For the mathematical foundations — Miller–Rabin primality (Fermat
> witnesses, strong pseudoprimes, the 2-Sylow structure argument), smooth-number theory
> (Canfield–Erdős–Pomerance density theorem), ECM (Lenstra's elliptic-curve method, the
> group-order smoothness argument, Montgomery-form curves), and Tonelli–Shanks (square roots mod
> $p$) — see `docs/MATHEMATICS.md §The α-Substrate: Primality, Smoothness, and ECM`. This chapter
> is the code-tour sibling: it assumes the reader knows the mathematics and focuses on the
> implementation in Rust.

This chapter explains the three shared crates `shared-field`, `shared-bigint`, and
`shared-numth`, which together form the arithmetic and number-theoretic substrate that the
algorithm crates build on. Like the number-field chapter (`shared/numfield/docs/PEDAGOGY.md`), it is
organised by the implementation, not by mathematical abstraction: each section states the
mathematics, then shows how it is realised in code. A reader who has read the `rho` PEDAGOGY chapter
can read this one without consulting the source.

The chapter covers prime-field arithmetic and its quadratic-residue machinery (the Legendre symbol
and the Tonelli–Shanks square root); Miller–Rabin primality testing; B-smoothness detection and the
`SmoothWitness` structure; Lenstra's Elliptic Curve Method (ECM) for factoring; and Montgomery's
batched field inversion.

These crates carry thorough code-level docstrings, but those do not give an integrative account.
The substrate's mathematical content — ECM in particular — is substantial enough to deserve one;
this chapter supplies it.

---

## 1. Introduction and Motivation

### What the α-substrate provides

The discrete-logarithm and factoring algorithms in this project all reduce, at the bottom, to three
kinds of computation:

1. **Field arithmetic** over GF(p) — addition, multiplication, inversion, powering, and the
   quadratic-residue operations (Legendre symbol, square root). This is `shared-field`.
2. **Number-theoretic primitives** — deciding primality, detecting smoothness, and factoring small
   cofactors. This is `shared-numth`.
3. **Cross-cutting arithmetic helpers** that do not belong to either — chiefly Montgomery's batched
   inversion, which amortises the cost of many field inversions. This is `shared-bigint`.

Each crate is consumed by multiple solvers. `shared-numth`'s smoothness predicate was designed
with three consumers in mind: NFS sieving, NFS discrete-log relation collection, and ECDLP index
calculus. ECM is a standalone factoring primitive; no solver in the workspace calls it
(Pohlig–Hellman's `factor_order` uses trial division instead). The `Fp` trait's square root is used
by pairing-friendly field constructions and anywhere a curve point must be recovered from its
x-coordinate.

### The crate names

In `Cargo.toml` the crates are `shared-field`, `shared-bigint`, and `shared-numth`; in Rust `use`
paths the hyphens become underscores (`shared_field::`, `shared_bigint::`, `shared_numth::`). Each
module is also `pub mod`, so both `shared_numth::is_prime` and `shared_numth::prime::is_prime` are
valid.

### Coefficient and limb conventions

Integer arithmetic is built on `crypto_bigint::Uint<L>`, a fixed-width unsigned integer of `L`
64-bit limbs. The substrate is generic over `L` where it can be (the `Fp<L>` trait), and fixed at
`Uint<4>` (256 bits) where the imminent consumers do not yet need more — see the C1 resolution in
the number-field chapter, which confirmed 256 bits suffices for toy-scale GNFS norms. The use of
`crypto-bigint` for the underlying big-integer arithmetic is a correctness-oracle dependency in the
same spirit as `num-bigint` in the number-field crate: the wide-integer arithmetic is delegated to a
well-tested library, while every algorithm above it is first-party.

---

## 2. Prime-Field Arithmetic: the `Fp<L>` Trait

### The trait and its two implementations

A prime field GF(p) is the set {0, 1, …, p−1} with arithmetic mod p, where p is prime so that every
non-zero element has a multiplicative inverse. The `Fp<L>` trait abstracts this, generic over the
limb count `L`, and is implemented twice: `FpNaive<L>` (schoolbook reduction after every operation)
and `FpMonty<L>` (Montgomery form, via `crypto-bigint`'s residue type). Both pass the same
known-answer tests; the `rho` chapter explains the Montgomery-form speedup in detail.

The trait's required surface is the standard field vocabulary: `zero`, `one`, `from_u64`,
`from_uint`, `to_uint`, `add`, `sub`, `neg`, `mul`, `square`, `pow`, and `inv`. Inversion is by
Fermat's little theorem (`a^(p−2) mod p`), which is why the trait assumes a prime modulus.

**Code reference:** `shared/field/src/lib.rs:66` (the `Fp<L>` trait); `FpNaive4` and `FpMonty4`
type aliases at lines 40 and 43.

### The Legendre symbol via Euler's criterion

For an odd prime p and an integer a not divisible by p, the **Legendre symbol** (a/p) is +1 if a is
a quadratic residue mod p (a square), −1 if it is a non-residue, and 0 if p | a. **Euler's
criterion** computes it directly:

> (a/p) ≡ a^((p−1)/2) (mod p)

The right side is +1 for a residue, p−1 (≡ −1) for a non-residue, and 0 when p | a. The
implementation evaluates this power and maps the three outcomes to the `i8` values {1, −1, 0}:

```rust
fn legendre(&self, p: &Uint<L>) -> i8;  // a^((p−1)/2) mod p → {0, 1, −1}
```

This is a **default method** on the `Fp` trait — written once in terms of `pow`, inherited by both
`FpNaive` and `FpMonty` with no per-implementation code. The cost is one modular exponentiation.

**Code reference:** `shared/field/src/lib.rs:140`.

### The Tonelli–Shanks square root

The companion operation is the modular square root: given a quadratic residue a, find r with
r² ≡ a (mod p). The implementation provides this as `sqrt`, returning `Option<Self>` — `Some(r)` for
a residue, `None` for a non-residue (the Legendre symbol decides which).

There are two paths, chosen by p mod 4:

- **p ≡ 3 (mod 4): the shortcut.** Here r = a^((p+1)/4) is a square root, because
  r² = a^((p+1)/2) = a · a^((p−1)/2) = a · (a/p) = a for a residue. One exponentiation, no search.
- **p ≡ 1 (mod 4): the full Tonelli–Shanks loop.** Write p − 1 = Q · 2^S with Q odd. The algorithm
  maintains a tuple (M, c, t, r) and repeatedly squares t to find the multiplicative order of the
  "error" term, correcting r by a precomputed power of a quadratic non-residue z at each step until
  t ≡ 1. Finding z is by trial: test n = 2, 3, 5, … until (n/p) = −1. For every prime this project
  uses, z is found within about ten trials.

```rust
fn sqrt(&self, p: &Uint<L>) -> Option<Self>;
// p ≡ 3 (mod 4): a^((p+1)/4)
// p ≡ 1 (mod 4): factor p−1 = Q·2^S, find QNR z, iterate Tonelli–Shanks
```

Like `legendre`, `sqrt` is a default method assuming an odd prime modulus; behaviour is unspecified
for composite p. Neither checks that p is prime: Tonelli–Shanks is correct only for an odd prime
modulus, so a caller that needs that guarantee certifies it first with `shared_numth::is_prime`.

**Code reference:** `shared/field/src/lib.rs:166`.

**A subtlety worth stating.** The trial search for a quadratic non-residue z is the one place the
square root is not a closed-form computation. It terminates quickly because non-residues have
density 1/2 — the probability that the first ten candidates {2, 3, 5, 7, 11, …} are *all* residues
is about 2⁻¹⁰. There is no known fast deterministic way to produce a non-residue (it is related to
open problems about the least quadratic non-residue under GRH), so the trial loop is the standard
practical choice. This is a small instance of a recurring pattern: an algorithm that is "almost"
closed-form but contains a search whose termination rests on a density argument.

---

## 3. Miller–Rabin Primality Testing

### The strong pseudoprime test

To decide whether an integer n is prime, the substrate uses the **Miller–Rabin** strong
pseudoprime test. Write n − 1 = 2^s · d with d odd. For a witness base a, compute x = a^d mod n,
then square it up to s − 1 times. If the sequence a^d, a^(2d), …, a^(2^(s−1)·d) ever hits n − 1
(or starts at 1), the witness a is *satisfied* and offers no evidence of compositeness; if it never
does, a is a *witness to compositeness* and n is definitely composite.

The mathematical content is that for a prime n, every a ∈ (1, n−1) is satisfied; for a composite n,
at least 3/4 of the bases are witnesses to compositeness. So a composite that passes k random bases
slips through with probability at most 4⁻ᵏ.

```rust
pub fn miller_rabin(n: &Uint<4>, witnesses: &[u64]) -> bool;
// n−1 = 2^s · d; for each a: x = a^d mod n, square up to s−1 times seeking n−1
```

**Code reference:** `shared/numth/src/prime.rs`, `miller_rabin`. Modular exponentiation routes
through `FpNaive::<4>::pow`.

### From probabilistic to deterministic

Miller–Rabin is probabilistic in general, but for bounded n it can be made *deterministic* by using
a fixed witness set known to have no exceptions below a threshold. The substrate uses the published
tables (Pomerance–Selfridge–Wagstaff, with refinements by Jaeschke and Sorenson–Webster, 2015). The
twelve prime bases {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} are a deterministic test for all
n < ψ₁₂ = 318,665,857,834,031,151,167,461 (about 2^78.1); adding the base 41, the thirteen bases are
a deterministic test for all n < ψ₁₃ = 3,317,044,064,679,887,385,961,981 (about 2^81.5), exported as
`MILLER_RABIN_PROOF_BOUND`. The constants are Theorem 1.1 of Sorenson and Webster, "Strong
pseudoprimes to twelve prime bases", *Math. Comp.* 86 (2017), arXiv:1509.00864, which was read for
them; the tests verify the boundary, not the exhaustive search behind it: ψ₁₂ =
399,165,290,221 · 798,330,580,441 is composite, the twelve bases accept it, and base 41 rejects it
(`psi_12_is_the_product_it_claims`, `twelve_bases_accept_psi_12_but_base_41_rejects_it`,
`psi_12_rejected`).

The public entry point `is_prime` selects the minimal sufficient witness set for the magnitude of n,
using the thirteen-base set from ψ₁₂ up. For n at or above ψ₁₃ (which still fits in `Uint<4>`) the
same thirteen bases give only a *probable-prime* verdict with no proven error bound. A "composite"
verdict is always a proof; a "prime" verdict is a proof only below ψ₁₃. Callers that need a proof
therefore relate their domain to that bound. `ecm_factor` compares a prime verdict against it and
reports `Unsupported` at or above it (§5). Rho factoring's supported range ends at 2^80, below ψ₁₃,
so every prime verdict there is a proof, and a compile-time assertion keeps the range from growing
past the bound.

```rust
pub fn is_prime(n: &Uint<4>) -> bool;
// selects the minimal deterministic witness set for n's magnitude
pub const MILLER_RABIN_PROOF_BOUND: u128; // ψ₁₃: below this, both verdicts are proofs
```

**Code reference:** `shared/numth/src/prime.rs`, `is_prime`; the witness table
`DETERMINISTIC_WITNESSES`.

The KAT coverage cross-checks `is_prime` against brute-force trial division for all n up to 1000
(no false negatives), recognises the first 100 primes and several Mersenne primes, accepts the prime
2^79 − 67 between the two bounds (`prime_above_twelve_base_bound_accepted`), and rejects the
classic base-2 pseudoprimes (e.g. 2047 = 23·89, which passes base 2 but fails base 23 — a test that
directly exercises the witness-set logic).

---

## 4. B-Smoothness and the `SmoothWitness`

### What smoothness is, and why it is central

An integer is **B-smooth** if all of its prime factors are ≤ B. Smoothness is the engine of every
index-calculus and sieve algorithm in this project: the whole strategy of NFS, NFS-DL, and ECDLP
index calculus is to collect integers (or norms, or curve points) that happen to be smooth over a
fixed factor base, because smooth objects factor completely over a small fixed set of primes and so
become rows in a linear system.

The substrate detects smoothness by trial division over a sorted prime factor base and records the
result in a `SmoothWitness`:

```rust
pub struct SmoothWitness {
    pub factors: Vec<(u64, u32)>,  // sorted (prime, exponent) pairs
    pub cofactor: Uint<4>,         // unfactored remainder; == 1 iff fully B-smooth
}
```

The `factors` list is the factorisation over the base; the `cofactor` is whatever is left after
dividing all base primes out. The integer is fully B-smooth exactly when `cofactor == 1`; a
`cofactor > 1` records a *partial* factorisation, where the remainder is a single large prime (or
product) outside the base — the data the large-prime variations need.

**Code reference:** `shared/numth/src/smooth.rs:47` (the struct); `is_smooth` (line 56, tests
`cofactor == 1`), `product` (line 63, reconstructs `cofactor · ∏ pᵉ`), and `verify` (line 77, the
round-trip check used in tests).

### The witness as a shared interface

The `SmoothWitness` shape is deliberately general: it was designed to serve three consumers —
integer smoothness with a prime factor base (NFS sieving, NFS discrete-log relation collection)
and, structurally similar though semantically different, smoothness of curve points via Semaev
polynomials (ECDLP index calculus). All three inspect `cofactor` to decide full versus partial
smoothness; the `factors` list becomes the exponent vector for the linear-algebra stage. The
witness was shaped for the index-calculus consumer even though that consumer is built much later
than the substrate.

### Trial division and the factor base

```rust
pub fn trial_smooth(n: &Uint<4>, factor_base: &[u64]) -> SmoothWitness;
pub fn factor_base_up_to(bound: u64) -> Vec<u64>;  // {p prime : p ≤ B}
```

`trial_smooth` divides n by each base prime as many times as it goes, accumulating exponents;
`factor_base_up_to` builds the base by testing each candidate with `is_prime`. The smoothness KAT
cross-checks `trial_smooth` against a brute-force reference for all n in 1..1000, and exercises the
partial-smoothness path (e.g. 77 = 7·11 over the base {2,3,5,7} leaves cofactor 11).

**Code reference:** `shared/numth/src/smooth.rs:97` (`trial_smooth`), `:145` (`factor_base_up_to`).

**Scale limit (under-exposed at toy scale).** `factor_base_up_to` builds the base by calling
`is_prime` on every candidate, which is O(B·√B) in the worst case. For the factor-base bounds used
at toy scale (B ≤ ~10⁴ in tests) this is invisible. At NFS scale, where B reaches the millions, a
sieve of Eratosthenes (O(B log log B)) is the correct algorithm, and the trial-division approach
would dominate setup time. The code carries this as an explicit TODO at the function and module
level. The *mathematics* of the factor base is unchanged by the choice of construction algorithm;
only the *performance* shifts, and only at a scale this project does not run at, so the gap is
invisible in the tests. The honest annotation is: correct as written, not production-scale, and
the gap is engineering, not mathematics.

---

## 5. Lenstra's Elliptic Curve Method (ECM)

ECM is the most mathematically substantial component of the α-substrate, and the reason this chapter
exists. It is a sub-exponential factoring method whose running time depends on the size of the
*smallest* prime factor rather than the size of N — which makes it the right tool for splitting off
small-to-medium prime factors. No other solver in the workspace calls it; it is exercised by its own
tests.

### The idea: Pollard's p−1, lifted onto a curve

Pollard's p−1 method factors N by exploiting that if p | N and p − 1 is B-smooth, then for a base a,
a^(k) ≡ 1 (mod p) once k = lcm(1, …, B) is a multiple of the order of a mod p. Then
gcd(a^k − 1, N) reveals p. The weakness is that it only works when p − 1 happens to be smooth — a
property of the single fixed group (ℤ/pℤ)*.

ECM (Lenstra 1987) removes that dependence on luck by replacing the fixed group (ℤ/pℤ)* with the
group of points on a *random* elliptic curve E over ℤ/Nℤ. The order |E(F_p)| varies with the curve
choice and ranges over an interval around p (Hasse's bound). If the chosen curve has B-smooth order
mod some prime p | N, then multiplying a starting point P by k = ∏ qᵉ (prime powers ≤ B) drives
k·P to the identity mod p while it stays non-trivial mod N. The group law's modular inversion then
fails with a denominator sharing a factor with N, and gcd(denominator, N) = p falls out.

The key advantage: if one curve's order is not smooth, you simply try another curve. ECM converts a
question of luck (is p − 1 smooth?) into a question you can re-roll (is *some* curve's order smooth?).

### Montgomery-form curves and the (X:Z) ladder

ECM uses **Montgomery-form curves** B·v² = u³ + A·u² + u, with the **Suyama parameterization** to
generate a valid (curve, starting point) pair from a single integer σ. Montgomery form admits
projective (X:Z) coordinates that omit the Y coordinate entirely: scalar multiplication needs only
the x-coordinate, computed by the **Montgomery ladder** using *differential addition* (given P, Q,
and P − Q, compute P + Q) and doubling. The ladder is uniform and efficient, and the Z = 0 point
represents the identity.

The Suyama parameterization computes, from σ, the curve constant A24 = (A+2)/4 and a starting point,
returning early with a factor if any of the constructions hits a non-invertible denominator mod N
(itself a lucky factor find).

**Code references** (the public surface, then the internals a reader will want):

| Item | Location | Role |
|------|----------|------|
| `ecm_factor(n, b1, b2, max_curves) -> Result<Uint<4>, EcmError>` | `ecm.rs` | top-level: classify `n`, then try up to `max_curves` Suyama curves (σ from 6), return a checked divisor or an `EcmError` (below) |
| `ecm_one_curve(n, sigma, b1, b2) -> EcmResult` | `ecm.rs` | one curve; `EcmResult::Factor(p)` or `Inconclusive` |
| `suyama_param(sigma, n)` | `ecm.rs` | Suyama parameterization → A24 and starting point |
| `montgomery_ladder(p, k, a24, n)` | `ecm.rs` | k·P via the (X:Z) ladder |
| `diff_add` / `point_double` | `ecm.rs` | Montgomery differential addition and doubling |
| `run_stage1` / `run_stage2` | `ecm.rs` | the two stages (below) |

**What `ecm_factor` returns.** `EcmError` is `InvalidInput`, `Prime`, `Unsupported`,
`CurvesExhausted { curves }`, or `Inconsistent { d }`. `ecm_factor` classifies `n` first, in this
order: `n ≤ 1` is `InvalidInput` (`zero_and_one_are_invalid_input`); `n = 2` is `Prime`
(`two_is_prime`); an even `n > 2` yields the divisor 2 with no curve run
(`even_n_yields_divisor_two`); an odd `n` is tested for primality, and a prime is `Prime` below
`MILLER_RABIN_PROOF_BOUND` (`odd_primes_are_prime`) but `Unsupported` at or above it, where the
verdict is only "probably prime" (`2^89 − 1` in `prime_above_proof_bound_is_unsupported`). The
composite ψ₁₂ is never certified prime (`composite_psi_12_is_never_prime`). Only an odd composite
reaches the curves. A divisor is returned only after `1 < d < n` and `n mod d = 0` are checked in
`Uint` arithmetic, independently of the gcd that produced it
(`judge_candidate_classifies_trivial_and_proper`). A curve that yields no candidate, or only `1`
or `n` (a gcd that is not proper), is a failed curve that consumes one of `max_curves`
(`search_trivial_candidates_consume_curves`). Any other candidate that does not divide `n` is
`Inconsistent { d }`: a gcd with `n` always divides `n`, so this is a breach of the divisor
invariant, reported at once and never retried (`judge_candidate_flags_non_dividing_candidates`,
`search_surfaces_inconsistent_candidate_without_retrying`,
`search_inconsistent_after_failed_curves`); the tests drive it through a stubbed curve, since no
input reaches it. `CurvesExhausted { curves }` counts the curves that ran, is a resource outcome
and proves nothing about `n` (`zero_max_curves_is_exhausted_at_zero`,
`exhaustion_counts_curves_run`). `ecm_factor` does not test the divisor it finds for primality.
All of these are unit tests in `shared/numth/src/ecm.rs`.

### Stage 1 and Stage 2

**Stage 1** computes Q = k·P with k = ∏ { q^⌊log_q(B1)⌋ : q prime ≤ B1 } — every prime power up to
B1. After each ladder call, gcd(Z-coordinate, N) is tested; a non-trivial gcd is the factor. Stage 1
succeeds when |E(F_p)| is fully B1-smooth.

*No product scalar.* The product k = lcm(1, …, B1) ≈ e^B1 exceeds 2^256 once B1 reaches 179
(256·ln 2 ≈ 177; the product first passes 256 bits when the prime 179 is multiplied in), so it
cannot live in a fixed-width integer at larger B1. `run_stage1` therefore never forms it: it applies
each prime power q = p^⌊log_p B1⌋ to the running point with its own `montgomery_ladder` call, one
prime at a time, and (q₁·q₂·…)·P is the same point however the factors are grouped. The test
`stage1_per_prime_power_equals_single_product_ladder` checks that for B1 up to 178 (every fifth
bound, plus 177 and 178) the result equals the single-scalar ladder (compared projectively,
X₁Z₂ = X₂Z₁), and
`single_product_scalar_overflows_at_179` pins where the single product stops fitting. The test
`stage1_prime_power_above_overflow_point_finds_factor` (N = (2^61 − 1)·1009, σ = 15, B1 = 250)
finds 1009 with the per-prime ladder while the product truncated mod 2^256 does not.

Testing gcd(Z, N) after every prime rather than once at the end is safe because once Z ≡ 0 (mod r)
for a prime r | N it stays 0: doubling and differential addition of points with Z ≡ 0 return
Z ≡ 0. It is also better, since the gcd is a proper divisor only while Z vanishes mod some primes of
N and not others; testing per prime reports the factor before every prime of N has reached the
identity.

**Stage 2** is the standard refinement for the common case where |E(F_p)| is B1-smooth *except for
one* prime factor q in (B1, B2]. Rather than raising B1 (which re-does all the work), stage 2 walks
the multiples q·Q for primes q ∈ (B1, B2] using differential addition with a fixed step (a
baby-step/giant-step structure with step size 2), accumulates the product of their Z-coordinates,
and takes a *single* gcd with N at the end. This catches the one-large-prime case cheaply. The KAT
`stage2_finds_factor_stage1_misses` (N = 1009·3541, B1 = 20, B2 = 200) asserts only that stage 1
plus stage 2 factors N; its own comment notes that stage 1 alone may or may not succeed on the
curves tried, so the test does not show that stage 1 alone fails.

**Code reference:** `shared/numth/src/ecm.rs` — `ecm_factor`, `EcmError`, `EcmResult`.

**Scale disconnect (under-exposed).** The ECM module carries three explicit TODOs:
the Brent–Suyama parameterization (an improved σ formula that raises the smoothness-hit rate), an
FFT-based stage-2 continuation, and stage-2 step-size optimisation. All three are *scale-only*
optimisations: their mathematical content is real, but their payoff only appears when ECM is run
at the scale where it factors large cofactors at speed. At the toy scale of the KATs (semiprimes up
to a few hundred thousand, B1 in the tens-to-hundreds) the simple Suyama parameterization and
linear stage-2 walk are entirely sufficient, and the optimised versions would demonstrate no
measurable difference. They are therefore documented and deferred, with the disconnect noted: the
omission is engineering-at-scale, not missing mathematics. ECM's *algorithmic* content — the
curve-group idea, the two-stage structure, the Montgomery ladder — is complete.

---

## 6. Batched Field Inversion

### Why batch inversions

Modular inversion is the most expensive field operation — typically tens of multiplications' worth.
Many algorithms (the affine elliptic-curve walks in the `rho` crate, future sieve and
linear-algebra steps) need *many* independent inversions at once. **Montgomery's batched inversion
trick** computes n inverses with a single inversion and 3(n−1) multiplications, trading n expensive
inversions for one inversion plus cheap multiplications.

The method is a forward/backward pass over prefix products:

1. **Forward:** build prefix products p[i] = x[0]·x[1]·…·x[i].
2. **Invert once:** compute (p[n−1])⁻¹, the only inversion.
3. **Backward:** recover each x[i]⁻¹ = (running inverse) · p[i−1], updating the running inverse by
   multiplying in x[i] as you walk back down.

```rust
pub fn batch_invert<const L: usize, F: Fp<L>>(xs: &mut [F], p: &Uint<L>);
// 1 inversion + 3(n−1) multiplications; inverts xs in place
```

It is generic over any `Fp<L>` implementation, lives in `shared-bigint` because it belongs to
neither the field-arithmetic nor the number-theory crate cleanly, and is the only public item that
crate currently re-exports. The empty-slice case is a no-op; a zero element has no inverse and is a
precondition violation. The `rho` chapter discusses how this same trick makes affine-coordinate
walks competitive with Jacobian ones.

**Code reference:** `shared/bigint/src/batch_inv.rs:25`.

### The `mp` stub

`shared/bigint/src/mp.rs` is an intentional placeholder — module docstring only, no code — reserved
for multi-precision helpers beyond `crypto-bigint` if gaps appear in that API.
It is noted here so a reader is not surprised by an empty module; it carries no mathematics yet.

**Code reference:** `shared/bigint/src/mp.rs:1`.

---

## 7. KAT Summary

The known-answer and property tests below pin the mathematical content of the α-substrate. (Field
QR tests are in `shared/field/tests/sqrt_legendre_kat.rs`; the `numth` tests are inline
`#[cfg(test)]` modules in each source file.)

| Test | Module | Mathematical fact verified |
|------|--------|---------------------------|
| `legendre_p5` … `legendre_p17` | field | (a/p) correct for all residues/non-residues mod small primes, both p ≡ 1 and p ≡ 3 mod 4 |
| `legendre_p1048517` | field | 2 is a non-residue when p ≡ 5 (mod 8) |
| `prop_legendre_multiplicative` | field | (ab/p) = (a/p)(b/p) exhaustively mod 17 |
| `sqrt_residues_p7` | field | Tonelli–Shanks p ≡ 3 (mod 4) shortcut path round-trips |
| `sqrt_residues_p1009` | field | exhaustive: every QR mod 1009 (full loop, p ≡ 1 mod 4) round-trips |
| `proptest_sqrt_roundtrip` | field | sqrt(x²)² ≡ x² for random x |
| `proptest_sqrt_none_for_qnr` | field | sqrt returns `None` exactly when (a/p) = −1 |
| `naive_monty_agree_p17` | field | `FpNaive` and `FpMonty` agree on legendre/sqrt for all a mod 17 |
| `first_100_primes`, `no_false_negatives_up_to_1000` | numth/prime | `is_prime` matches brute-force trial division below 1000 |
| `miller_rabin_base2_accepts_2047` / `…_base23_rejects_2047` | numth/prime | 2047 = 23·89 is a base-2 pseudoprime caught by base 23 (witness-set logic) |
| `mersenne_primes_recognised` | numth/prime | Mersenne primes recognised |
| `smooth_60_over_235` | numth/smooth | 60 = 2²·3·5 fully smooth over {2,3,5}; exponents correct |
| `partial_77_over_2357` | numth/smooth | 77 = 7·11 leaves cofactor 11 (partial smoothness) |
| `matches_brute_force_up_to_1000` | numth/smooth | `trial_smooth` matches a brute-force reference for all n < 1000 |
| `round_trip_verify`, `product_with_cofactor` | numth/smooth | `SmoothWitness::product` reconstructs the original integer |
| `ecm_factors_143` … `ecm_factors_near_2_32` | numth/ecm | ECM splits semiprimes from 11·13 up to 211·223 |
| `stage2_finds_factor_stage1_misses` | numth/ecm | N = 1009·3541, B1 = 20, B2 = 200 — stage 1 plus stage 2 factors N (does not show stage 1 alone fails) |
| `stage2_extends_reach` | numth/ecm | N = 101·103, B1 = 10, B2 = 200 — a tight B1 with stage 2 still factors N |
| `zero_and_one_are_invalid_input`, `two_is_prime`, `odd_primes_are_prime` | numth/ecm | `InvalidInput`, and `Prime` for 2 and for odd primes below the proof bound |
| `even_n_yields_divisor_two` | numth/ecm | an even n > 2 gives the divisor 2 with no curve |
| `prime_above_proof_bound_is_unsupported`, `composite_psi_12_is_never_prime` | numth/ecm | 2^89 − 1 is `Unsupported`; the composite ψ₁₂ is never certified prime |
| `zero_max_curves_is_exhausted_at_zero`, `exhaustion_counts_curves_run`, `judge_candidate_classifies_trivial_and_proper`, `search_trivial_candidates_consume_curves` | numth/ecm | `CurvesExhausted` counts the curves run; a divisor is returned only after the `1 < d < n`, `n mod d = 0` check; the trivial candidates `1` and `n` consume a curve |
| `judge_candidate_flags_non_dividing_candidates`, `search_surfaces_inconsistent_candidate_without_retrying`, `search_inconsistent_after_failed_curves` | numth/ecm | a candidate other than `1` or `n` that does not divide `n` is `Inconsistent { d }`, ends the search at once and is never retried (stubbed curves; no input reaches it) |
| `psi_12_is_the_product_it_claims`, `twelve_bases_accept_psi_12_but_base_41_rejects_it`, `psi_12_rejected` | numth/prime | ψ₁₂ = 399165290221·798330580441 is accepted by the twelve bases, rejected by base 41 and by `is_prime` |
| `prime_above_twelve_base_bound_accepted` | numth/prime | the prime 2^79 − 67, between ψ₁₂ and ψ₁₃, is accepted |
| `batch_inv_matches_individual` | bigint | batched inverse of a 10-element batch matches individual `inv` calls |
| `batch_inv_empty` | bigint | empty batch is a no-op, not a panic |

---

## 8. Downstream Consumption

| Component | Consumed by |
|-----------|-------------|
| `Fp<L>` arithmetic | everything with field arithmetic — `rho`, number fields (which build on `Fp`), the NFS discrete-log solvers |
| `Fp::legendre`, `Fp::sqrt` | pairing-friendly prime construction (needs p ≡ 1 mod 4 square roots), curve point recovery |
| `is_prime`, `miller_rabin` | `factor_base_up_to`, Pohlig–Hellman group-order analysis, polynomial-selection screening |
| `trial_smooth`, `SmoothWitness` | NFS sieving, NFS discrete-log relation collection, ECDLP index calculus |
| `ecm_factor` | no caller in the workspace other than its own tests (Pohlig–Hellman's `factor_order` uses trial division) |
| `MILLER_RABIN_PROOF_BOUND` | the bound `ecm_factor` compares a prime verdict against; rho factoring asserts at compile time that its supported limit 2^80 does not exceed it |
| `batch_invert` | `rho` affine walks; future sieve/linear-algebra inner loops |

The smoothness witness is the substrate's most load-bearing interface across solvers; the
number-field chapter (§10) tracks its evolution. Its current `Uint<4>` / `(u64, u32)` shape was
confirmed sufficient for toy-scale norms in the number-field arithmetic.

---

## Further Reading

1. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective.* 2nd ed.
   Springer. The standard computational reference: Miller–Rabin and deterministic witness sets
   (§3.4–3.5), Tonelli–Shanks (§2.3), smoothness and the factor base (§3.2), and ECM (§7.4).

2. **Lenstra, H. W. Jr. (1987).** "Factoring integers with elliptic curves." *Annals of
   Mathematics*, 126(3), 649–673. The original ECM paper.

3. **Montgomery, P. L. (1987).** "Speeding the Pollard and elliptic curve methods of factorization."
   *Mathematics of Computation*, 48(177), 243–264. Montgomery-form curves, the (X:Z) ladder and
   differential addition, and the batched-inversion trick.

4. **Suyama, H. (1985).** Informal note on the parameterization used to generate ECM curves with
   guaranteed group-order divisibility; described in Crandall–Pomerance §7.4.

5. **Pomerance, C., Selfridge, J. L., and Wagstaff, S. S. Jr. (1980).** "The pseudoprimes to
   25·10⁹." *Mathematics of Computation*, 35(151), 1003–1026. The witness-set tables behind
   deterministic Miller–Rabin, with later refinements by Jaeschke (1993) and Sorenson–Webster
   (2017, "Strong pseudoprimes to twelve prime bases", *Mathematics of Computation*, 86; its
   Theorem 1.1 gives ψ₁₂ and ψ₁₃).

6. **Bach, E., and Shallit, J. (1996).** *Algorithmic Number Theory, Vol. 1.* MIT Press. The
   quadratic-residue machinery (Legendre/Jacobi symbols, Euler's criterion) and the least-non-residue
   questions behind the Tonelli–Shanks trial search.
