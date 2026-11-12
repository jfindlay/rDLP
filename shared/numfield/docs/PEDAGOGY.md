# Number-Field Substrate: A Code-Tour Chapter

This chapter explains the `shared/numfield` crate — the algebraic substrate that the General Number
Field Sieve (GNFS) and related algorithms sit on. It is organised by the implementation, not by
mathematical abstraction: each section introduces the mathematics, then shows how it is realised in
code. A reader who has read the `rho` PEDAGOGY chapter can read this one without consulting the
source.

The chapter covers polynomial arithmetic, number fields and their elements, norm and trace, the
resultant and subresultant GCD, ideal representation, and Dedekind factorisation — including bad
primes and the Dedekind criterion.

---

## 1. Introduction and Motivation

### Why number fields arise in GNFS

The General Number Field Sieve factors a large integer N by exploiting a homomorphism between two
rings: the rational integers ℤ and the ring of integers of a number field K = ℚ(α). The key idea
is to find many pairs (a, b) ∈ ℤ² such that both a + bm (the rational norm, for a suitable integer
m) and the algebraic value, Norm_{K/ℚ}(a + bα) for monic f and c times that norm for leading
coefficient c (see below), are B-smooth — divisible only by primes up to a bound B. When enough
such pairs are collected, linear algebra over 𝔽₂ finds candidate dependencies; a verified square
root of one gives a congruence of squares mod N, whose gcd with N is a non-trivial factor unless
the congruence is trivial.

The algebraic side of the sieve requires computing Norm_{K/ℚ}(a + bα) for each sieve pair (a, b).
For monic f this norm is an integer. For a leading coefficient c ≠ 1 the integer is instead
c · Norm_{K/ℚ}(a + bα) = (−b)^d f(−a/b), because the norm itself can have a denominator. Its prime
factorisation determines which algebraic prime ideals divide the ideal (a + bα) in ℤ[α]. The factor base on the algebraic side consists of prime ideals above
rational primes p ≤ B — specifically, the prime ideals (p, α − r) where r is a root of f mod p.

This is why the `shared/numfield` crate exists: it provides the arithmetic of K = ℚ(α), the norm
map, and the ideal-theoretic machinery that the sieve consumes.

### What ℤ[α] gives us that ℤ does not

The rational integers ℤ have unique factorisation. The ring ℤ[α] — polynomials in α with integer
coefficients, where α satisfies a monic irreducible f ∈ ℤ[x] — is in general not a unique
factorisation domain. However, the ring of integers ℤ_K of K = ℚ(α) is a Dedekind domain: every
non-zero ideal factors uniquely into prime ideals. The GNFS exploits this ideal-theoretic unique
factorisation on the algebraic side.

For the toy-scale GNFS implemented in this project, ℤ[α] is used directly (not the full ring of
integers ℤ_K, which may be larger). The Dedekind criterion (§9) detects when ℤ[α] ≠ ℤ_K at a
given prime p; for the polynomials used in practice, this distinction rarely matters at toy scale.

### Coefficient types: a correctness-oracle dependency

All polynomial and element arithmetic in this crate uses `num_bigint::BigInt` for integer
coefficients and `num_rational::BigRational` for rational coefficients. This is a deliberate
design choice, analogous to the role of CADO-NFS as a correctness oracle for the sieve: the
underlying integer and rational arithmetic is delegated to a well-tested third-party library, while
the number-field abstraction and all algorithms above it are first-party. The alternative —
implementing arbitrary-precision arithmetic from scratch — would be an engineering optimisation
orthogonal to the mathematical content of the crate.

---

## 2. Polynomial Arithmetic: `IntPoly` and `RatPoly`

### Data representation

Both polynomial types store coefficients in a `Vec` with the **least-significant coefficient
first**: `coeffs[i]` is the coefficient of xⁱ. The zero polynomial has an empty `coeffs` vector.
The **trailing-zero invariant** is maintained throughout: `from_coeffs` trims any trailing zeros
before storing, and all arithmetic operations produce trimmed results. This means the degree of a
non-zero polynomial is always `coeffs.len() - 1`, and the zero polynomial has `degree() == None`.

```rust
// f = x² − 2: coefficients [−2, 0, 1] (constant term first)
let f = IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)]);
assert_eq!(f.degree(), Some(2));
assert_eq!(f.leading_coeff(), Some(&bi(1)));
```

The trailing-zero invariant is not merely cosmetic: it ensures that degree comparisons and leading
coefficient lookups are O(1), and that equality is coefficient-wise without a normalisation step.

### `IntPoly` — polynomials over ℤ

`IntPoly` provides the arithmetic of ℤ[x]: addition, subtraction, negation, multiplication, and
scalar multiplication. All operations return new polynomials; the inputs are not mutated. Evaluation
uses Horner's method for efficiency.

The most important non-trivial operation is **pseudo-division** (`pseudo_div_rem`). Over ℤ[x],
ordinary long division may require dividing by the leading coefficient of the divisor, which is not
always an integer. Pseudo-division avoids this by pre-multiplying the dividend by a power of the
leading coefficient of the divisor:

> Given f and g with leading coefficient lc(g), pseudo-division computes
> lc(g)^e · f = q · g + r, where e = deg(f) − deg(g) + 1.

The exponent e is chosen so that every division step in the long-division loop is exact over ℤ.
This is the building block for the subresultant GCD algorithm (§6).

```rust
// f = 2x² + 3x + 1, g = 2x + 1.
// lc(g) = 2, e = 2, multiplier = 4.
// 4f = 8x² + 12x + 4 = (4x + 4)(2x + 1) + 0.
let (q, r) = f.pseudo_div_rem(&g);
// q = 4x + 4, r = 0
```

`IntPoly` also provides `to_rat_poly`, which embeds ℤ[x] into ℚ[x] by converting each `BigInt`
coefficient to a `BigRational`.

### `RatPoly` — polynomials over ℚ

`RatPoly` provides the arithmetic of ℚ[x], with the same storage convention. The key operation
that `IntPoly` lacks is **exact polynomial long division** (`div_rem`): over ℚ, the leading
coefficient of the divisor is always a unit, so division is always possible without pre-scaling.

```rust
// (x² − 1) / (x − 1) = x + 1, remainder 0.
let (q, r) = dividend.div_rem(&divisor);
```

The `rem` method is a convenience wrapper that discards the quotient.

### Why both types are needed

`IntPoly` is used for:
- The defining polynomial f ∈ ℤ[x] (stored in `NumberField`).
- The resultant computation (§6), which stays in ℤ throughout via the Bareiss algorithm.
- The discriminant and Dedekind criterion (§8, §9), which work mod p over ℤ.

`RatPoly` is used for:
- Element arithmetic in K = ℚ(α): elements have rational coefficients in general.
- Reduction mod f (the `rem` operation), which requires exact division.
- The extended Euclidean algorithm for inversion (§3).

The embedding `IntPoly::to_rat_poly` bridges the two: the defining polynomial f lives in ℤ[x] but
must be treated as an element of ℚ[x] when reducing element products mod f.

---

## 3. Number Fields and Elements

### The number field K = ℚ(α)

A number field K = ℚ(α) is determined by a monic irreducible polynomial f ∈ ℤ[x] of degree d ≥ 1.
The field is the quotient ring ℚ[x]/(f(x)): elements of K are equivalence classes of polynomials
in ℚ[x] under the relation g ≡ h iff f | (g − h). The primitive element α is the image of x in
this quotient; it satisfies f(α) = 0 by construction.

`NumberField` wraps the defining polynomial:

```rust
pub struct NumberField {
    pub f: IntPoly,  // monic, irreducible, degree ≥ 1
}
```

The constructor panics if f is not monic or has degree < 1. It does **not** verify irreducibility —
this is the caller's responsibility, enforced by KAT coverage rather than a runtime check (which
would require a polynomial factorisation algorithm).

The extension degree [K : ℚ] = deg(f) is the dimension of K as a ℚ-vector space. The standard
basis is {1, α, α², …, α^{d−1}}.

### Elements of K

`NumberFieldElement<'a>` represents an element of K as a polynomial in α of degree strictly less
than d, with rational coefficients:

```rust
pub struct NumberFieldElement<'a> {
    pub field: &'a NumberField,
    pub poly: RatPoly,  // degree < field.degree()
}
```

The lifetime `'a` ties the element to its ambient field. The degree bound `poly.degree() < d` is
an invariant maintained by all arithmetic operations.

Constructing elements:

```rust
let k = NumberField::new(IntPoly::from_coeffs(vec![bi(-2), bi(0), bi(1)])); // x² − 2
let alpha = k.alpha();          // the element α (polynomial x mod f)
let one = k.from_int(bi(1));    // the element 1 (constant polynomial)
let beta = one.add(&alpha);     // 1 + α
```

Addition and subtraction cannot increase the degree beyond d − 1 (since both operands have degree
< d), so they require no reduction. Negation is coefficient-wise.

### The Mul canonicalisation contract

Multiplication is the critical operation. The product of two polynomials of degree < d has degree
up to 2d − 2, which exceeds the degree bound. `mul` therefore **eagerly reduces mod f** after
polynomial multiplication:

```rust
pub fn mul(&self, rhs: &Self) -> Self {
    let product = self.poly.mul(&rhs.poly);
    Self::reduced(self.field, product)  // reduces product mod f via RatPoly::rem
}
```

The `reduced` helper calls `poly.rem(&f_rat)`, which performs exact polynomial long division over
ℚ and returns the remainder. This remainder has degree < d by the definition of polynomial
division.

The consequence of eager reduction is that **equality is coefficient-wise**: two elements are equal
if and only if their `poly` fields are equal. No re-reduction is needed before comparison. This is
the canonical-form invariant of `Mul`: it is what makes `PartialEq` for `NumberFieldElement`
correct and cheap.

```rust
let alpha_sq = alpha.square();
// α² ≡ 2 mod (x² − 2), so alpha_sq.poly == [2/1] (the constant 2)
assert_eq!(alpha_sq, k.from_int(bi(2)));
```

`square` and `pow` (square-and-multiply) both delegate to `mul` and therefore inherit the
canonicalisation contract.

### Inversion via extended Euclidean in ℚ[x]

Since f is irreducible over ℚ, for any non-zero β ∈ K the polynomials β.poly and f are coprime in
ℚ[x]. The extended Euclidean algorithm finds s, t ∈ ℚ[x] such that:

> s · β.poly + t · f = 1

The Bézout coefficient s is then the inverse of β.poly modulo f. The implementation uses the
standard iterative extended GCD over ℚ[x], normalising the GCD to be monic:

```rust
pub fn inv(&self) -> Self {
    let f_rat = self.field.f.to_rat_poly();
    let (_, s, _) = extended_gcd_rat(&self.poly, &f_rat);
    Self::reduced(self.field, s)
}
```

The final `reduced` call is defensive: s should already have degree < d (since it is a Bézout
coefficient for polynomials of degree < d and d respectively), but the reduction ensures the
invariant is maintained regardless.

---

## 4. Norm via Resultant

### The algebraic definition

The field norm N_{K/ℚ}(β) of an element β ∈ K is the product of all conjugates of β:

> N_{K/ℚ}(β) = ∏ᵢ β(αᵢ)

where α₁, …, αd are the d roots of f (in some algebraic closure). For β represented as a
polynomial g(α) ∈ ℚ[x], this product equals the resultant of f and g with respect to x:

> N_{K/ℚ}(g(α)) = Res_x(f, g)

This identity holds because the resultant of f and g is, by definition, the product ∏ᵢ g(αᵢ)
times the leading coefficient of f raised to deg(g) — and since f is monic, the leading coefficient
factor is 1.

### The Sylvester matrix

The resultant Res(f, g) is the determinant of the **Sylvester matrix**, a (d + e) × (d + e) matrix
where d = deg(f) and e = deg(g). The matrix is built from two blocks of shifted coefficient rows:
the first e rows are shifts of f (by x⁰, x¹, …, x^{e−1}), and the next d rows are shifts of g
(by x⁰, x¹, …, x^{d−1}). Coefficients appear in descending order of degree within each row.

For the norm computation, deg(g) < deg(f) = d, so the Sylvester matrix is (d + e) × (d + e) with
e < d. The implementation in `element.rs` builds this matrix with `BigRational` entries and
computes its determinant via Gaussian elimination with partial pivoting:

```rust
fn sylvester_resultant(f: &IntPoly, g: &RatPoly) -> BigRational {
    let d = f.degree().unwrap();
    let e = g.degree().unwrap_or(0);
    let n = d + e;
    // Build the (d+e)×(d+e) Sylvester matrix...
    det_rational(&mut mat, n)
}
```

Since the arithmetic is over ℚ (exact), Gaussian elimination is numerically exact. The result is
always a rational number, but for algebraic integers (elements of ℤ[α]) it is always an integer.

### The sign convention

A common source of confusion is the sign of the norm for linear elements. For monic f of degree d
and the element β = α − c (for c ∈ ℤ), the norm is:

> N_{K/ℚ}(α − c) = Res_x(f, x − c) = (−1)^d · f(c)

The factor (−1)^d arises from the Sylvester matrix construction: the resultant of f (degree d) and
g = x − c (degree 1) equals (−1)^{d·1} times the product of (αᵢ − c) over all roots αᵢ of f,
which is (−1)^d · f(c) for monic f.

The sign is easy to get wrong: a first draft of the norm KAT stated Norm(α − 1) = −1 for
f = x³ − x − 1, but the correct value is:

> N(α − 1) = (−1)³ · f(1) = −1 · (1 − 1 − 1) = −1 · (−1) = 1

The implementation is mathematically correct; the KAT spec was corrected to match. The formula
Norm(α − c) = (−1)^d · f(c) is the reliable reference.

### The numerical-embedding alternative

An alternative approach to computing the norm is via numerical embeddings: embed K into ℂ^d by
sending α to each of the d roots of f (computed numerically), evaluate β at each embedding, and
multiply the results. This approach is simpler to implement but introduces floating-point error,
which must be controlled by working with sufficient precision and rounding to the nearest integer.
For the GNFS use case, where norms can be large (up to ~2^{150} for toy-scale parameters), the
rounding step is non-trivial. The resultant approach is exact and is therefore the implementation
chosen here.

---

## 5. Trace

The field trace Tr_{K/ℚ}(β) is the sum of all conjugates of β:

> Tr_{K/ℚ}(β) = ∑ᵢ β(αᵢ)

For β = g(α) with g ∈ ℚ[x], the trace equals the trace of the matrix g(C), where C is the
**companion matrix** of f. The companion matrix is the d × d matrix whose characteristic polynomial
is f; its eigenvalues are exactly the roots α₁, …, αd of f.

The companion matrix of a monic polynomial f = xd + c_{d−1} x^{d−1} + … + c₀ is:

```
C = [ 0   0  …  0  −c₀  ]
    [ 1   0  …  0  −c₁  ]
    [ 0   1  …  0  −c₂  ]
    [ …   …  …  …   …   ]
    [ 0   0  …  1  −c_{d−1} ]
```

(subdiagonal of 1s, last column of negated lower coefficients). The trace of g(C) is computed by
evaluating the polynomial g at the matrix C using Horner's method, then summing the diagonal
entries.

For the standard basis element α (g = x), the trace is the sum of the roots of f, which equals
−c_{d−1} by Vieta's formulas. For the constant element 1 (g = 1), the trace is d (the degree of
the extension).

The trace is used less frequently than the norm in the GNFS context, but it appears in the
computation of the discriminant (via the trace form) and in lattice-based algorithms.

---

## 6. Resultant and Subresultant GCD

The `resultant.rs` module provides a standalone, general-purpose resultant and GCD over ℤ[x],
separate from the norm-specific Sylvester computation in `element.rs`. This module is consumed by
the discriminant computation (§8) and by the Dedekind criterion (§9).

### Resultant via Sylvester matrix + Bareiss determinant

The public `resultant(f, g)` function builds the Sylvester matrix of f and g as a flat
row-major `Vec<BigInt>` and computes its determinant using the **Bareiss fraction-free algorithm**.

Bareiss elimination avoids rational arithmetic entirely. The key identity is:

> M[i][j] ← (M[i][j] · pivot − M[i][col] · M[col][j]) / prev_pivot

where `prev_pivot` is the pivot from the previous column. This division is always exact over ℤ
(this is the content of the Bareiss theorem), so the entire computation stays in ℤ without
introducing fractions. The result is the determinant of the Sylvester matrix, which equals
Res(f, g) ∈ ℤ.

The Bareiss algorithm is preferred over Gaussian elimination with rational arithmetic for the
integer resultant because it avoids the coefficient explosion that would occur if intermediate
results were kept as fractions.

### Subresultant GCD via pseudo-remainder sequence

The `subresultant_gcd(f, g)` function computes the primitive GCD of f and g over ℤ[x] — a
polynomial proportional to gcd(f, g) over ℚ[x], with positive leading coefficient.

The naive approach — Euclidean algorithm using `pseudo_div_rem` at each step — suffers from
**coefficient explosion**: the coefficients of intermediate remainders grow exponentially in the
degree. The subresultant PRS (pseudo-remainder sequence) avoids this by reducing each remainder to
its primitive part before the next step:

```rust
loop {
    let (_, r) = a.pseudo_div_rem(&b);
    if r.degree().is_none() {
        return primitive_part(&b);  // b is the GCD
    }
    a = b;
    b = primitive_part(&r);  // remove content before next step
}
```

The `primitive_part` function divides all coefficients by their GCD (the content of the
polynomial) and normalises the leading coefficient to be positive. This keeps coefficients small
throughout the sequence. The result is the primitive part of the last non-zero remainder, which is
the primitive GCD of f and g.

The `pseudo_div_rem` building block (in `poly.rs`) is the engine of this algorithm. It computes
lc(g)^e · f = q · g + r over ℤ, where the pre-multiplication by lc(g)^e ensures all divisions in
the long-division loop are exact.

---

## 7. Ideal Representation

### Two-element primary form

An ideal I in ℤ[α] is represented in **two-element primary form**: I = (p, α − r), the ideal
generated by the rational prime p and the element α − r for some r ∈ ℤ. This representation is
standard in the NFS literature and is sufficient for the prime ideals that appear in the algebraic
factor base.

```rust
pub struct Ideal<'a> {
    pub field: &'a NumberField,
    pub p: BigInt,   // rational prime (positive)
    pub r: BigInt,   // integer r; second generator is α − r
}
```

The invariant `p > 0` is enforced by the constructor. The ideal (p, α − r) is a prime ideal of
ℤ[α] above the rational prime p when r is a root of f mod p — that is, when f(r) ≡ 0 (mod p).

### Ideal norm

The norm of the ideal I = (p, α − r) is defined as the index [ℤ[α] : I] — the number of cosets
of I in ℤ[α] as an abelian group. For a prime ideal above p with **residue degree 1** (meaning
the quotient ℤ[α]/I is isomorphic to 𝔽_p), the norm is simply p:

```rust
pub fn norm(&self) -> BigInt {
    self.p.clone()
}
```

The residue degree 1 convention is the standard NFS assumption: the factor base consists only of
prime ideals (p, α − r) where r is a root of f mod p, which have residue degree 1. Prime ideals
above p arising from irreducible factors of f mod p of degree > 1 have higher residue degree and
are not included in the factor base at toy scale.

### Ideal multiplication via CRT

The product of two prime ideals I = (p₁, α − r₁) and J = (p₂, α − r₂) above distinct primes
p₁ ≠ p₂ is represented as (p₁p₂, α − r) where r is determined by the **Chinese Remainder
Theorem**:

> r ≡ r₁ (mod p₁),   r ≡ r₂ (mod p₂)

The CRT solution is:

> r = r₁ + p₁ · ((r₂ − r₁) · p₁⁻¹ mod p₂)

reduced into [0, p₁p₂). This is implemented in `Ideal::mul`:

```rust
pub fn mul(&self, rhs: &Self) -> Self {
    let p = &self.p * &rhs.p;
    let diff = &rhs.r - &self.r;
    let p1_inv_mod_p2 = mod_inverse(&self.p, &rhs.p).expect("p₁ and p₂ must be coprime");
    let t = mod_reduce(&(diff * &p1_inv_mod_p2), &rhs.p);
    let r = mod_reduce(&(&self.r + &self.p * t), &p);
    Self { field: self.field, p, r }
}
```

The norm-multiplicativity identity N(IJ) = N(I) · N(J) follows immediately: N(IJ) = p₁p₂ =
N(I) · N(J).

### What the two-element form can and cannot represent

The two-element form (p, α − r) is the ideal of ℤ[α] generated by p and α − r. When r is a root of
f mod p, it is a prime ideal of ℤ[α] of norm p: ℤ[α]/(p, α − r) ≅ 𝔽_p. It is not the prime of
residue degree 1 for an arbitrary pair (p, r): the form carries no check that f(r) ≡ 0 (mod p), and
when f has no root mod p, (p, α − r) is the unit ideal for every r (it contains f(r), because
α − r divides f(α) − f(r) = −f(r), and p ∤ f(r), so it contains gcd(p, f(r)) = 1). The form
cannot represent:
- Prime ideals of residue degree > 1 (arising from irreducible factors of f mod p of degree > 1).
  Their second generator is a lift of that factor evaluated at α, not α − r for an integer r.
- Ideals generated by something other than a rational integer and one linear element, in general.
- Products of ideals above the same prime p (the CRT formula requires p₁ ≠ p₂).

It also says nothing about the ring of integers ℤ_K: a prime of ℤ[α] of norm p is a prime of ℤ_K
only when p does not divide the index [ℤ_K : ℤ[α]] (§9).

For a complete ideal-theoretic treatment, one would use the **Hermite Normal Form (HNF) basis**
representation, which can represent any ideal as a ℤ-module. This is noted as future work; the
two-element form is sufficient for degree-one primes.

---

## 8. Dedekind Factorisation

### Dedekind's theorem

Given a monic polynomial f ∈ ℤ[x] and a rational prime p that does not divide the index
[ℤ_K : ℤ[α]] (in particular any p ∤ disc(f)), **Dedekind's theorem** states that the prime ideals
of ℤ_K above p correspond bijectively to the irreducible factors of f mod p. Specifically:

- If f ≡ ∏ᵢ gᵢ(x)^{eᵢ} (mod p) with gᵢ irreducible over 𝔽_p, then the prime ideals above p are
  𝔭ᵢ = (p, gᵢ(α)).
- For a linear factor gᵢ = x − r, the prime ideal is 𝔭ᵢ = (p, α − r) in two-element form.
- The ramification index eᵢ and residue degree fᵢ = deg(gᵢ) satisfy ∑ eᵢ fᵢ = d = deg(f).

If p ∤ disc(f), f is squarefree mod p (all eᵢ = 1), so the factorisation is unramified. A bad
prime may still ramify, with some eᵢ > 1, and still not divide the index (§9).

### The degree-one report

`dedekind_factor(field, p, limits)` reports the degree-one part of this factorisation. It reduces f
mod p, takes the square-free decomposition f̄ = ∏ sₑ^e, and enumerates the roots of the radical
∏ sₑ with `FpPoly::roots` under the caller's `RootSearchLimits`. Each root r lies in exactly one
sₑ, and that e is its multiplicity as a root of f̄. The result is a `DedekindReport`:

```rust
pub struct DedekindReport<'a> {
    pub degree_one: Vec<DegreeOnePrime<'a>>, // (p, α − r) with multiplicity, ascending r
    pub nonlinear_degree: usize,             // deg f − Σ multiplicity
    pub is_bad_prime: bool,                  // p | disc(f)
    pub index: IndexTest,                    // does p divide [ℤ_K : ℤ[α]]?
}
```

The invariant is Σ multiplicity + nonlinear_degree = deg f. Irreducible factors of degree > 1 have
no two-element form (§7), so they appear only through `nonlinear_degree`: when p ∤ index it equals
Σ eᵢfᵢ over the primes with fᵢ > 1.

How to read `degree_one` depends on the index test.

- If p does not divide the index, each listed ideal is a prime of ℤ_K of residue degree 1, and its
  multiplicity is its ramification index.
- If p divides the index, each listed ideal is only a prime of ℤ[α] of norm p, and its
  multiplicity is the multiplicity of r in f̄, not a ramification index. Nothing is claimed about
  the factorisation of (p) in ℤ_K.

### Inert primes and exhausted searches

If f has no roots mod p — all irreducible factors of f mod p have degree > 1 — the report has an
empty `degree_one` and `nonlinear_degree = deg f`. No ideal stands in for it: a placeholder such as
(p, α − 0) would be the unit ideal whenever f has no root mod p (§7), which is not a prime at all,
and could not be told apart from a real degree-one prime above p without re-evaluating f.

The root count is computed exactly before any scan (`deg gcd(f̄, x^p − x)`), so "no roots" is a
proved statement and costs no budget. When roots exist but the evaluation budget runs out before
all are found, `dedekind_factor` returns `DedekindError::RootSearchExhausted`, carrying the exact
distinct-root count and the evaluations spent, rather than a report whose multiplicities and
`nonlinear_degree` would be wrong.

---

## 9. Bad Primes and the Dedekind Criterion

### Two notions: bad primes and primes dividing the index

Two distinct conditions on a prime p are easy to conflate:

- p is **bad** for f if p | disc(f). This is what `is_bad_prime` tests for monic f; for a non-monic
  f the `discriminant` function returns lc(f)·disc(f), so `is_bad_prime` also flags a prime that
  divides only the leading coefficient (this is how the GNFS factor base uses it). The index
  statements below concern monic f. The discriminant of a monic polynomial f of degree d is:

  > disc(f) = (−1)^{d(d−1)/2} · Res(f, f')

  where f' is the formal derivative. The `discriminant` function computes this via the standalone
  `resultant` from `resultant.rs`:

  ```rust
  pub fn discriminant(f: &IntPoly) -> BigInt {
      let f_prime = formal_derivative(f);
      let res = resultant(f, &f_prime);
      let sign_exp = d * (d - 1) / 2;
      if sign_exp % 2 == 0 { res } else { -res }
  }
  ```

- p **divides the index** if p | [ℤ_K : ℤ[α]], that is, if ℤ[α] is strictly smaller than the full
  ring of integers ℤ_K at p. This is what the Dedekind criterion decides.

The two are linked by disc(f) = [ℤ_K : ℤ[α]]² · disc(K): every prime dividing the index divides
disc(f), so is bad. The converse fails. A bad prime is a necessary but not sufficient condition for
dividing the index.

Dedekind's theorem (§8) reads the factorisation of (p) off f mod p exactly when p does not divide
the index. For a prime dividing the index, f mod p must be non-squarefree (if it were squarefree,
h̄ = 1 and T = 1 in the criterion below), though the converse fails, and the primes of ℤ[α] above p
need not be the primes of ℤ_K.

### The Dedekind criterion

The **Dedekind criterion** (also called the index criterion; Cohen, *A Course in Computational
Algebraic Number Theory*, §6.1.2) decides, for any prime p, whether p divides the index. The
algorithm, as implemented by `dedekind_index_test`, is:

1. Factor f̄ = f mod p as ∏ tᵢ^eᵢ over 𝔽_p, with tᵢ distinct monic irreducibles.
2. Let ḡ = ∏ tᵢ (the radical of f̄) and h̄ = f̄ / ḡ.
3. Lift ḡ and h̄ to monic integer polynomials g and h. Any lifts will do.
4. Form F = (f − g·h) / p, computed from the **original integer polynomial f**. The numerator is
   divisible by p coefficient-wise because g·h ≡ f (mod p) by construction.
5. Compute T = gcd(F̄, ḡ, h̄) in 𝔽_p[x], normalised monic.
6. p divides [ℤ_K : ℤ[α]] if and only if deg T > 0.

Two points are easy to get wrong.

- **F comes from f, not from f mod p.** F̄ depends on f modulo p², but reducing f first keeps only
  f mod p, which does not determine it: x² − 2 and x² − 12 both reduce to x² mod 2 yet have
  different obstructions, as the example below shows.
- **The lifts do not matter.** Replacing g by g + p·u changes F by −u·h, and replacing h by
  h + p·v changes F by −g·v; after reduction these are multiples of h̄ and ḡ respectively, so the
  gcd with ḡ and h̄ is unaffected. The implementation lifts with representatives in [0, p) and a
  unit test shifts the lifts by multiples of p to confirm T is unchanged.

The factorisation in step 1 is only needed up to the radical, so the code uses the square-free
decomposition of f̄ (`FpPoly::radical`) rather than a full irreducible factorisation. This
handles multiplicities that are multiples of p, where the characteristic-0 shortcut f̄ / gcd(f̄, f̄')
fails: a factor tᵉ with p | e has zero derivative, so it is divided out entirely and is missed.

### The ℤ[√2] example

The polynomial f = x² − 2 has discriminant 8, so p = 2 is a bad prime. The Dedekind criterion
gives T = 1 at p = 2:

- f̄ = x² (since −2 ≡ 0 mod 2), so ḡ = x and h̄ = x.
- Lift g = h = x, so g·h = x².
- F = (f − g·h)/2 = (−2)/2 = −1, and F̄ = 1.
- T = gcd(1, x, x) = 1.

So p = 2 does not divide [ℤ_K : ℤ[√2]]. This is correct: ℤ[√2] is already the full ring of
integers of ℚ(√2), whose discriminant is 8 = disc(f), so the index is 1. The prime 2 is totally
ramified (2 = (√2)² up to a unit) yet ℤ[√2] is maximal at 2. The prime 2 is bad (it divides the
discriminant) but does not divide the index.

Contrast f = x² − 12 at p = 2, where α = 2√3. Again f̄ = x², ḡ = h̄ = x, g = h = x, but now
F = (f − x²)/2 = −6, so F̄ = 0 and T = gcd(0, x, x) = x has degree 1. Here 2 divides the index:
ℤ[2√3] has index 2 in ℤ[√3]. A computation that reduced f mod 2 to x² before forming F would
get F = 0 for *both* polynomials and could not tell them apart.

### The `IndexTest` type

`dedekind_index_test(f, p)` takes a `PrimeModulus` and returns an `IndexTest`, which records the
prime and the obstruction T. `IndexTest::divides_index` is true iff deg T > 0. Errors are
`DedekindError::NotMonic` and `DedekindError::ConstantPolynomial` for invalid input, and
`DedekindError::Field` for a failure inside the prime-field arithmetic.

`dedekind_factor` (§8) stores this test as the `index` of its report, next to `is_bad_prime`.
`index.divides_index()` implies `is_bad_prime`, never the other way round.

### Cross-checking against PARI/GP

The opt-in oracle suite `dedekind_pari_oracle_kat.rs` (run with `cargo test -p shared-numfield
--test dedekind_pari_oracle_kat -- --ignored`; it fails if `gp` is missing) compares the criterion
and the report with PARI/GP on irreducible fixtures of degrees 2–4 at p ∈ {2, 3, 5, 7}. PARI
computes v_p([ℤ_K : ℤ[α]]) = (v_p(disc f) − v_p(disc K)) / 2 from `poldisc` and `nfdisc`, and
`divides_index` must hold exactly when it is positive. Where p does not divide the index,
`idealprimedec` gives the true factorisation of (p) in ℤ_K: each degree-one prime is matched to a
report entry by the residue r with `idealval(nf, x − r, 𝔭) > 0`, and its ramification index must be
the entry's multiplicity, while Σ eᵢfᵢ over the primes with fᵢ > 1 must be `nonlinear_degree`.
Whether or not p divides the index, the listed roots and multiplicities must be those of `factormod`
applied to f, since that is all the report claims at an index-dividing prime.

### Beyond the criterion: Round 2

When `index.divides_index()` is true, ℤ[α] is strictly smaller than ℤ_K at p, the listed ideals are
primes of ℤ[α] rather than of ℤ_K, and the correct factorisation of (p) in ℤ_K requires
enlarging the order, which is the job of the **Round 2 algorithm** (also called Zassenhaus's
algorithm, built on HNF computation over ℤ-modules). It is mentioned here for context and is not
implemented in this crate: the reason is engineering scale, not mathematical omission.

Primes dividing the index cannot simply be excluded from a factor base: the prime ideals of ℤ_K
above p still occur in the factorisations of norms, and the ideals of ℤ[α] above p do not
describe them. Such primes must be handled by Round 2 or reported as unresolved.

---

## 10. Downstream Use of the Interfaces

### Which components consume each interface

The interfaces of this crate are consumed by several other components:

| Interface | Consumed by |
|-----------|-------------|
| `NumberField` / `NumberFieldElement` arithmetic | ideals, resultants, Dedekind factorisation, polynomial selection, sieving, NFS-DL, the p-adic substrate |
| prime ideals `(p, α − r)` | Dedekind factorisation, sieving, NFS-DL |
| resultants | the norm, polynomial selection, NFS-DL |
| Dedekind factorisation | no solver yet |

**Polynomial selection** uses `IntPoly` and the resultant to compute discriminants and
polynomial scores. The defining polynomial f is selected to have small coefficients and a large
number of roots mod small primes — both properties that depend on the resultant and discriminant
machinery.

**Sieving** uses this crate less than one might expect. For each sieve pair (a, b) it evaluates
the homogeneous form F(a, b) = b^d·f(a/b) = c·Norm(a − bα) directly from the coefficients of the
original f (`gnfs/src/sieve/norms.rs`), not `NumberFieldElement::norm()`, and factors that value
over the algebraic factor base. The factor base does not call `dedekind_factor`: it finds the
roots of f mod p by its own search and takes only the `is_bad_prime` test from this crate. The
degree-one report of §8 is therefore not an input to the factor base today.

**NFS-DL** builds the field elements a − bα with `NumberField` (`gnfs/src/dl/relation.rs`), which
sieving does not, and otherwise shares the factor-base and ideal machinery with sieving; the main
difference is in the linear algebra step.

**Pairings** do **not** use `NumberField` or `NumberFieldElement` for extension-field
arithmetic.  Pairing target fields are `F_{p^k}` — finite fields of **characteristic p**, not
char-0 `ℚ[x]/(f)` number fields in the sense of this crate.  The degree-2 extension `F_{47^2}`
used by the toy Weil/Tate fixture (and the degree-12 / degree-6 extensions of BN/BLS12 curves at
crypto scale) are `F_{p^k}`, not `ℚ`-extensions.  The char-0 `NumberFieldElement` arithmetic is
**not** directly reusable: the coefficient field differs (`F_p` vs `ℚ`).  The real char-p
substrate is `rho/src/pairing/fpext.rs` (`FpExt<F>`, a polynomial quotient `F_p[u]/(m(u))`
built over the `shared::field` `Fp<4>` prime field).  (The p-adic substrate is separate from
pairings.)

### `Uint<4>` sufficiency

Is the `shared::numth` crate's `Uint<4>` (256-bit integers) sufficient for GNFS algebraic values at
toy scale? It is:

For a toy 60–80 bit semiprime N with factor base bound B ≤ 10⁴:
- Polynomial degree d = 3–5, coefficients C ~ N^{1/d} ~ 2^{20–27}.
- Sieve bound M ~ 10⁴ = 2^{13}.
- Norm bound: |Norm| ≤ (d+1)^{(d+1)/2} · (M · C)^d ≈ 2^{120–150}.

**256 bits suffices.** The `shared::numth` crate remains at `Uint<4>` with no widening needed.
For larger N (e.g., 128-bit semiprimes) the bound would exceed 256 bits and `Uint<4>` would need
to be widened — the `Fp` trait is designed to support this widening without breaking downstream
consumers.

---

## 11. KAT Summary

The following table lists the key known-answer tests across the crate, with the mathematical fact
each one verifies. All tests are in `shared/numfield/tests/`.

| Test | File | Mathematical fact verified |
|------|------|---------------------------|
| `kat1_product_1_plus_alpha_times_1_minus_alpha` | `numfield_kat.rs` | (1+α)(1−α) = −1 in ℚ(√2); Mul canonicalisation and reduction mod f = x²−2 |
| `kat2_norm_1_plus_alpha_in_sqrt2` | `numfield_kat.rs` | N_{ℚ(√2)/ℚ}(1+√2) = −1; Sylvester resultant for norm |
| `kat3a_norm_alpha_cubic` | `numfield_kat.rs` | N(α) = 1 for f = x³−x−1; Norm(α) = (−1)^d · f(0) = (−1)³·(−1) = 1 |
| `kat3b_norm_alpha_minus_1_cubic` | `numfield_kat.rs` | N(α−1) = 1 for f = x³−x−1; sign convention Norm(α−c) = (−1)^d · f(c) |
| `alpha_squared_equals_2_in_sqrt2` | `numfield_kat.rs` | α² ≡ 2 mod (x²−2); eager reduction in `square` |
| `inv_round_trip_in_sqrt2` | `numfield_kat.rs` | β · β⁻¹ = 1; extended Euclidean inversion in ℚ[x] |
| `trace_of_alpha_is_zero_in_sqrt2` | `numfield_kat.rs` | Tr(√2) = 0; companion matrix trace computation |
| `trace_of_one_is_degree_in_sqrt2` | `numfield_kat.rs` | Tr(1) = [K:ℚ] = 2; trace of the identity element |
| `cubic_alpha_cubed_equals_alpha_plus_1` | `numfield_kat.rs` | α³ = α+1 in ℚ(α) for f = x³−x−1; `pow` via square-and-multiply |
| `norm_kat_prime_ideal` | `ideal_kat.rs` | N((5, α−2)) = 5; ideal norm = p for residue degree 1 |
| `norm_multiplicativity_kat` | `ideal_kat.rs` | N(IJ) = N(I)·N(J) = 35 for I=(5,α−2), J=(7,α−3) |
| `crt_consistency_kat` | `ideal_kat.rs` | CRT: r ≡ 2 (mod 5), r ≡ 3 (mod 7) gives r = 17 in [0,35) |
| `kat_resultant_1_shared_root` | `resultant_kat.rs` | Res(x²−1, x−1) = 0; shared root implies zero resultant |
| `kat_resultant_2_coprime_quadratics` | `resultant_kat.rs` | Res(x²−2, x²−3) = 1; coprime polynomials have non-zero resultant |
| `kat_resultant_3_norm_consistency` | `resultant_kat.rs` | Res(x²−2, x+1) = −1; cross-check with the field norm |
| `kat_gcd_1_shared_factor` | `resultant_kat.rs` | gcd(x²−1, x−1) ∝ x−1; subresultant PRS finds shared factor |
| `kat_gcd_2_coprime_polynomials` | `resultant_kat.rs` | gcd(x²−2, x²−3) is constant; coprime polynomials have trivial GCD |
| `kat1_inert_prime_has_no_degree_one_primes` | `dedekind_kat.rs` | x²−2 at 3: no roots, so no degree-one primes and nonlinear degree 2; no ideal such as (3, α), which is the unit ideal, is listed |
| `kat2_split_prime` | `dedekind_kat.rs` | x²−2 at 7: roots 3 and 4, each of multiplicity 1, ascending |
| `kat3_cubic_partial_split` | `dedekind_kat.rs` | x³−x−1 at 5: only root 2 (the other factor is an irreducible quadratic), nonlinear degree 2 |
| `kat4_norm_product` | `dedekind_kat.rs` | ∏ N(𝔭)^e = p^d = 49 for totally split p=7 in a degree-2 field |
| `kat5_discriminant_quadratic` | `dedekind_kat.rs` | disc(x²−2) = 8; discriminant via resultant |
| `kat6_discriminant_cubic` | `dedekind_kat.rs` | disc(x³−x−1) = −23; standard reference value |
| `kat7_is_bad_prime` | `dedekind_kat.rs` | 2 | disc(x²−2) = 8 (bad); 3 ∤ 8 (good) |
| `kat8_ramified_prime_with_maximal_order` | `dedekind_kat.rs` | x²−2 at 2: root 0 of multiplicity 2, bad, index test negative: ℤ[√2] is maximal so 2 is totally ramified |
| `kat9_degree_one_report_when_the_index_is_divisible` | `dedekind_kat.rs` | x³−10 at 3: root 1 of multiplicity 3 and index test positive; the multiplicity is that of f̄, not a ramification index |
| `kat_divisible_index_with_a_double_root` | `dedekind_kat.rs` | x²−12 at 2: root 0 of multiplicity 2 with positive index test |
| `kat_bad_prime_with_a_maximal_order` | `dedekind_kat.rs` | x²−3 at 3: root 0 of multiplicity 2, bad prime, index test negative |
| `kat_totally_ramified_cyclotomic_prime` | `dedekind_kat.rs` | Φ₅ at 5: root 1 of multiplicity 4, nonlinear degree 0, ℤ[ζ₅] maximal |
| `kat_irreducible_cubic_reduction_is_all_nonlinear` | `dedekind_kat.rs` | x³−x−1 at 2: x³+x+1 has no roots, so nonlinear degree 3 |
| `kat_roots_ascend_across_multiplicity_classes` | `dedekind_kat.rs` | roots from different squarefree parts are merged in ascending order; irreducible quadratic leaves nonlinear degree 2 |
| `kat_multiplicity_equal_to_the_characteristic` | `dedekind_kat.rs` | (x−1)³(x+1) mod 3: multiplicity 3 = p comes from the p-th-power branch of the decomposition |
| `kat_truncated_budget_is_an_error_not_a_claim` | `dedekind_kat.rs` | x²−2 at 7: budgets 0–4 give `RootSearchExhausted` with the exact root count 2 and the evaluations spent; 5 succeeds |
| `kat_rootless_reduction_needs_no_budget` | `dedekind_kat.rs` | the root count is exact before scanning, so an inert prime is reported under a zero budget |
| `kat_multiplicities_plus_nonlinear_degree_equal_the_degree` | `dedekind_kat.rs` | Σ multiplicity + nonlinear degree = deg f over thirteen (f, p) fixtures; roots ascend and are roots mod p |
| `report_index_matches_the_standalone_index_test` | `dedekind_kat.rs` | the report's `index` equals `dedekind_index_test` on the same inputs |
| `index_divisible_x2_minus_12_at_2` | `dedekind_kat.rs` | F = −6 from the integer f gives T = x: 2 divides the index of ℤ[2√3] |
| `index_not_divisible_x2_minus_3_at_3` | `dedekind_kat.rs` | F = −1 gives T = 1: ℤ[√3] is 3-maximal; reducing f mod 3 first would wrongly give T = x |
| `index_not_divisible_x2_minus_2_at_2` | `dedekind_kat.rs` | T = 1: ℤ[√2] is 2-maximal though 2 divides disc |
| `index_divisible_x2_plus_3_at_2` | `dedekind_kat.rs` | T = x+1: 2 divides the index of ℤ[√−3] |
| `index_divisible_x3_minus_10_at_3` | `dedekind_kat.rs` | T = x+2 (10 ≡ 1 mod 9): multiplicity 3 = p combined with a non-trivial obstruction |
| `index_not_divisible_x3_minus_2_at_3`, `…_at_2` | `dedekind_kat.rs` | ℤ[∛2] is maximal: T = 1 at 3 and at 2 |
| `index_not_divisible_cyclotomic_5_at_5` | `dedekind_kat.rs` | totally ramified f̄ = (x−1)⁴ yet T = 1: ℤ[ζ₅] is maximal |
| `index_not_divisible_at_good_primes` | `dedekind_kat.rs` | p ∤ disc(f) gives T = 1 |
| `index_test_rejects_invalid_polynomials` | `dedekind_kat.rs` | non-monic input is `NotMonic`; zero or constant input is `ConstantPolynomial` |
| `kat_pari_dedekind_fixtures_are_irreducible_and_monic` | `dedekind_pari_oracle_kat.rs` | (ignored; needs `gp`) every oracle fixture is monic and `polisirreducible` in PARI, so the criterion's precondition holds |
| `kat_pari_dedekind_index_criterion_matches_the_index_valuation` | `dedekind_pari_oracle_kat.rs` | (ignored; needs `gp`) `divides_index` ⇔ v_p([ℤ_K : ℤ[α]]) > 0, with the valuation from `poldisc`/`nfdisc`; `is_bad_prime` ⇔ p ∣ disc(f); both index-dividing and bad-but-maximal primes occur |
| `kat_pari_dedekind_report_lists_the_roots_of_the_reduction` | `dedekind_pari_oracle_kat.rs` | (ignored; needs `gp`) the report's roots and multiplicities are the linear factors of `factormod(f, p)`, and Σ multiplicity + `nonlinear_degree` = deg f, at every tested p |
| `kat_pari_dedekind_report_matches_idealprimedec_when_p_does_not_divide_the_index` | `dedekind_pari_oracle_kat.rs` | (ignored; needs `gp`) degree-one primes (matched by `idealval(nf, x − r, 𝔭) > 0`) have e = multiplicity, and Σ e·f over f > 1 = `nonlinear_degree`, including ramified-but-maximal primes |
| `kat_pari_dedekind_report_claims_only_the_reduction_when_p_divides_the_index` | `dedekind_pari_oracle_kat.rs` | (ignored; needs `gp`) at primes dividing the index the report asserts only f̄'s linear factors and the degree sum, including inseparable f̄ and a repeated irreducible quadratic factor |
| `kat_prime_modulus_accepts_small_and_boundary_primes` | `fp_poly_kat.rs` | 2, 3, 5, 7, 2³¹−1, 2⁶¹−1 and 2⁶⁴−59 are accepted as moduli |
| `kat_prime_modulus_rejects_units_and_composites` | `fp_poly_kat.rs` | 0, 1, 4, 9, 561 and 2⁶⁴−1 are rejected with `NotPrime` |
| `kat_prime_modulus_rejects_strong_pseudoprimes` | `fp_poly_kat.rs` | each listed strong pseudoprime equals an explicit product of factors, so it is composite independently of Miller–Rabin |
| `kat_is_prime_u64_matches_trial_division_below_100000` | `fp_poly_kat.rs` | `is_prime_u64` agrees with a sieve of Eratosthenes on every n < 100000 |
| `kat_scalar_arithmetic_near_u64_max` | `fp_poly_kat.rs` | p = 2⁶⁴−59: (p−1)+(p−1), (p−1)², inverses, Fermat, `i64::MIN` and negative `BigInt` reduction do not overflow |
| `kat_canonical_representation` | `fp_poly_kat.rs` | reduction, trimming, least nonnegative residues, degree drop of 2x²+x+1 mod 2, modulus-sensitive equality |
| `kat_modulus_mismatch_is_error` | `fp_poly_kat.rs` | every binary operation on different moduli returns `ModulusMismatch` |
| `kat_div_rem_identity_and_degree` | `fp_poly_kat.rs` | a = q·b + r with r = 0 or deg r < deg b; division by zero is an error |
| `kat_gcd_and_extended_gcd` | `fp_poly_kat.rs` | gcd((x−1)(x−2),(x−1)(x−3)) = x−1 mod 7; gcd(0,0) = 0; s·a + t·b = gcd exactly |
| `kat_pow_mod_matches_repeated_multiplication_and_frobenius` | `fp_poly_kat.rs` | x^(p^d) ≡ x mod an irreducible f; a > 2⁶⁴ exponent reduces mod p^d−1 |
| `kat_derivative_characteristic_p` | `fp_poly_kat.rs` | d/dx x^p = 0; p·x^(p−1) vanishes mod p |
| `kat_irreducible_fixtures` | `fp_poly_kat.rs` | Rabin's test accepts independently irreducible quadratics, cubics, GF(16) and GF(64) moduli, and linear polynomials |
| `kat_reducible_fixtures_without_roots` | `fp_poly_kat.rs` | products of two irreducibles (no roots) are squarefree but not irreducible |
| `kat_repeated_factor_fixtures` | `fp_poly_kat.rs` | squares, cubes and p-th powers are not squarefree and not irreducible |
| `kat_reducible_with_roots` | `fp_poly_kat.rs` | x²+1 mod 5, x³−x−1 mod 5 and x²+1 mod 2⁶⁴−59 are reducible with the stated root counts |
| `kat_degree_and_zero_handling` | `fp_poly_kat.rs` | zero and unit polynomials are not irreducible; zero has no root count |
| `kat_distinct_root_count` | `fp_poly_kat.rs` | deg gcd(f, x^p−x) counts distinct roots, including repeated and root-free cases |
| `kat_root_enumeration_complete` | `fp_poly_kat.rs` | ascending complete enumeration with exact evaluation counts; agrees with brute force over small fields |
| `kat_root_enumeration_truncated` | `fp_poly_kat.rs` | a spent budget returns the exact root count and exactly the roots found so far |
| `kat_valid_fields_report_characteristic_degree_and_order` | `residue_field_kat.rs` | x²+1 over F₃, x³−2 over F₇, a linear modulus and x⁴+x+1 over F₂ are fields of order 9, 343, p and 16 |
| `kat_non_monic_modulus_is_scaled_to_monic` | `residue_field_kat.rs` | 2x²+2 and x²+1 over F₃ define the same field |
| `kat_cardinality_exceeds_machine_words` | `residue_field_kat.rs` | p² for p = 2⁶⁴−59 and p³ for p = 2⁶¹−1 are exact `BigUint` orders beyond `u64` |
| `kat_invalid_moduli_are_rejected_with_the_right_error` | `residue_field_kat.rs` | zero, constant, rooted-reducible, rootless-reducible and repeated-factor moduli each get their own error |
| `kat_integer_polynomial_entry_checks_the_prime_and_the_reduction` | `residue_field_kat.rs` | non-prime characteristics are rejected; a degree drop mod p gives a lower-degree field; a reduction to zero is rejected |
| `kat_error_display_and_source` | `residue_field_kat.rs` | every error variant has a message; only the wrapped prime-field error has a source |
| `kat_elements_over_another_characteristic_are_rejected` | `residue_field_kat.rs` | every operation rejects an element over a different prime |
| `kat_arithmetic_reduces_modulo_the_modulus` | `residue_field_kat.rs` | in F₃[x]/(x²+1): x² = −1, (1+x)² = 2x, unreduced spellings reduce, u⁸ = 1 |
| `kat_inverse_by_multiplication` | `residue_field_kat.rs` | (1+x)⁻¹ = 4+x in F₅[x]/(x²−2) by hand and by product; large-prime inverses |
| `kat_zero_has_no_inverse_but_is_a_square_with_root_zero` | `residue_field_kat.rs` | zero (and the modulus as an unreduced zero) is uninvertible, a square, and has root zero without any search |
| `kat_exhaustive_even_degree_fields` | `residue_field_kat.rs` | `is_square`, `sqrt` and `inv` agree with enumeration of all of F₉, F₂₅, F₄₉, F₁₂₁, F₁₆₉, F₂₈₉ |
| `kat_exhaustive_odd_degree_and_prime_fields` | `residue_field_kat.rs` | the same cross-check for F₂₇, F₃₄₃, F₁₇, F₉₇, F₂₅₇ (long Tonelli–Shanks loops, 2-adic valuation up to 8) |
| `kat_exhaustive_characteristic_two` | `residue_field_kat.rs` | every element of GF(2), GF(8), GF(16) is a square with a verified root |
| `kat_roots_of_degree_two_squares_are_plus_or_minus_the_witness` | `residue_field_kat.rs` | the root of w² is ±w, including over F_(p²) with p = 2⁶⁴−59 |
| `kat_deep_two_adic_valuation_over_a_machine_prime` | `residue_field_kat.rs` | F_(p²), p = 2⁶¹−1: v₂(q−1) = 62; squareness agrees with the norm criterion and roots square back |
| `kat_odd_degree_over_a_machine_prime` | `residue_field_kat.rs` | F_(p³), p = 2⁶¹−1: roots of squares; scalar nonresidues of F_p stay nonsquares |
| `kat_large_prime_field_roots_and_nonsquares` | `residue_field_kat.rs` | p = 2⁶⁴−59: Euler-criterion agreement, 2 is a nonsquare, −1 has a root |
| `kat_unreduced_inputs_are_accepted` | `residue_field_kat.rs` | an unreduced element is treated as its class |
| `kat_characteristic_two_has_unique_roots_and_no_nonresidues` | `residue_field_kat.rs` | GF(16): the root of a is a⁸, there is no nonresidue, and a zero budget still succeeds |
| `kat_search_order_and_budget_boundary` | `residue_field_kat.rs` | odd-degree search order: constants, first nonresidue 3 in F₁₇ and F₇[x]/(x³−2); budgets below it are exhausted |
| `kat_even_degree_search_skips_constants` | `residue_field_kat.rs` | even-degree search order: monic linear polynomials; F_(p²) finds a nonsquare quickly despite p−1 square constants |
| `kat_budget_exhaustion_is_distinct_from_nonsquare` | `residue_field_kat.rs` | in F₂₅ a zero budget exhausts exactly the 9 squares that need a nonresidue; nonsquares and zero are unaffected |
| `kat_too_small_budget_reports_candidates_spent` | `residue_field_kat.rs` | in F₁₇ budgets 0, 1, 2 report that many candidates tried and budget 3 succeeds |
| `kat_three_mod_four_needs_no_nonresidue` | `residue_field_kat.rs` | q ≡ 3 mod 4: every root is found with a zero budget |

---

## 12. Polynomials over a Prime Field

`fp_poly` is the characteristic-p counterpart of `IntPoly`. It supplies the arithmetic that residue
fields F_p[x]/(f) and root-based prime searches are built on. Its design choices are about making
wrong answers unrepresentable rather than about speed.

### Canonical values and checked operations

`PrimeModulus` can only be built from a prime: `PrimeModulus::new` runs `is_prime_u64` and returns
`NotPrime` otherwise, so downstream code never re-validates. `FpPoly` stores `(modulus, coeffs)` with
private fields. Coefficients are least-significant first, each in `[0, p)`, trailing zeros trimmed,
and the zero polynomial is the empty slice. Every constructor reduces its input (negative integers
take their least nonnegative residue), so structural equality is mathematical equality, and
equality includes the modulus. Reducing a non-monic integer polynomial can lower its degree
(2x²+x+1 is x+1 over F₂), so a caller that needs a degree preserved checks `degree()`.

Operations return `Result`: mixed moduli give `ModulusMismatch`, a zero divisor or modulus gives
`DivisionByZero`, and nothing panics. These are invalid-input signals, never mathematical
certificates.

### Overflow safety

The moduli range over every prime below 2⁶⁴, so a product of two residues needs 128 bits. Scalar
sums and products widen to `u128` and reduce at once. Polynomial products reduce after every
multiply-add instead of accumulating unreduced sums. The tests exercise p = 2⁶⁴−59 (the largest
prime below 2⁶⁴), where (p−1)·(p−1) and (p−1)+(p−1) both overflow `u64`.

### Deterministic primality

`is_prime_u64` is a proof, not probable primality. After trial division by 2, 3, …, 37 it applies
the strong-probable-prime test to the first twelve primes as bases. The smallest composite that
passes all twelve bases is ψ₁₂ = 318665857834031151167461 (Sorenson and Webster, 2017), which exceeds
2⁶⁴, so no `u64` composite passes. Eleven bases are not enough: 3825123056546413051 < 2⁶⁴ is a
strong pseudoprime to the first nine, ten and eleven prime bases (Jiang and Deng, 2014). The
pseudoprime tests pin that boundary from the composite side by factoring each listed number
explicitly. `numth::prime` carries the same table for `Uint<4>`; calling it from here would add a
dependency to this crate, so the small `u64` version is duplicated deliberately.

### Irreducibility and root counting

Rabin's test decides irreducibility exactly. For monic g of degree n over F_p, g is irreducible iff
x^(pⁿ) ≡ x (mod g) and gcd(x^(p^(n/r)) − x, g) = 1 for every prime r dividing n. The powers are n
successive p-th powers modulo g. The first condition rejects repeated factors because x^(pⁿ) − x is
squarefree. A root search never replaces this test: a reducible quartic like (x²+1)(x²+x+2) over F₃
has no roots at all.

`distinct_root_count` is deg gcd(f, x^p − x), since x^p − x is the product of (x − a) over F_p. It
costs O(log p) polynomial products regardless of p, so it is exact and unbudgeted. A polynomial
splits completely into distinct linear factors iff it is squarefree and its count equals its degree.
Consumers decide splitting that way and never from an enumerated root list.

### Budgeted enumeration

`roots` evaluates candidates 0, 1, 2, … in order, one evaluation each, and stops at the count-th
root or at the budget. The outcome is a `RootEnumeration`: `Complete` carries all roots and
`Truncated` carries exactly the roots below the evaluations performed, together with the exact
count. Exhausting the budget is a resource outcome, never evidence that roots are absent, and there
is no resume: retry with a larger budget. For primes near 2⁶⁴ the budget turns an O(p) scan into a
reported resource limit instead of a hang.

### Square-free decomposition

`squarefree_decomposition` writes f = lc(f) · ∏ sᵢ^eᵢ with each sᵢ monic, squarefree, of positive
degree and pairwise coprime, listed with strictly ascending eᵢ; sᵢ is the product of every
irreducible factor of f whose multiplicity is exactly eᵢ. `radical` is the product of the sᵢ.

Over ℚ, Yun's algorithm alone suffices: c = gcd(f, f′) carries every factor to one power less, and
repeatedly splitting w = f/c against c peels off the sᵢ in order. Over F_p it can miss factors,
because the derivative of a factor of multiplicity divisible by p vanishes: if q^p divides f then
q^p stays entirely inside c, with no trace in w. Yun's loop therefore returns only the parts whose
multiplicity is prime to p, and what remains of c is a p-th power. Since Frobenius fixes F_p,
h(x)^p = h(x^p) for every h over F_p, so that remainder is c = h(x^p) and h is read off by taking
the coefficients of c at indices divisible by p. Decomposing h recursively and multiplying its
multiplicities by p finishes the job. The same step handles f′ = 0 outright (then c = f).

The mixed case is the one a Yun-only implementation gets wrong: x²(x+1)³(x+2)⁶ over F₃ yields x at
multiplicity 2 from the main loop, and the remainder ((x+1)(x+2)²)³ yields x+1 at 3 and x+2 at 6
from the recursion. Multiplicities from the main loop are never divisible by p, and those from the
p-th-root levels have strictly increasing p-adic valuation, so no two levels produce the same
multiplicity; entries are still merged by multiplicity defensively.

---

## 13. Residue Fields and Square Roots

`residue_field` builds the field F_q = F_p[x]/(f), q = pᵈ, on top of `fp_poly` and takes square
roots in it. The checked constructor is the whole validity story: `ResidueField::new` rejects the
zero polynomial, a nonzero constant (the zero ring), and any `f` for which `is_irreducible` is false
(a root test would miss rootless reducibles such as a product of two quadratics). It scales `f` to
monic, so a non-monic modulus and its monic multiple are the same field. After that, no operation
can silently work in a ring with zero divisors.

### Elements and exponents

An element is a plain `FpPoly`. Every operation first reduces its operands modulo `f`, so an
unreduced polynomial is just another spelling of its class, and a polynomial over a different prime
is a `ModulusMismatch`. Nothing borrows the field, which lets a caller keep elements after the field
is gone. The cardinality q and the exponents derived from it — (q−1)/2, the odd part t of q−1,
q/2 — are `BigUint`, because already q = p² exceeds a machine word for p near 2⁶⁴. `inv` runs the
extended Euclidean algorithm against `f` and returns `DivisionByZero` for zero.

### Squareness and Tonelli–Shanks

In odd characteristic the Euler criterion is exact: a nonzero a is a square iff a^((q−1)/2) = 1, and
otherwise the power is −1. `is_square` counts zero as a square. Equivalently a is a square in F_q iff
its norm to F_p is a square in F_p, which the tests use as an independent oracle.

Write q−1 = 2ˢ·t with t odd. Tonelli–Shanks starts from r = a^((t+1)/2) and tt = aᵗ, and repeatedly
multiplies r by a power of c = zᵗ, for a nonresidue z, until tt = 1. Each pass lowers m (initially
s) strictly, so there are at most s passes: large 2-adic valuations of q−1 are slow but terminate.
Two consequences:

- A nonresidue is needed only if aᵗ ≠ 1. For q ≡ 3 mod 4 (s = 1) it never is, so those roots need
  no search at all.
- In characteristic 2 the Frobenius map is a bijection, so every element is a square with the
  unique root a^(q/2).

Every root is squared and compared with the input before it is returned, so an unverified root is
unrepresentable; a failure there is reported as an invariant violation rather than a wrong root.

### The bounded nonresidue search

`find_nonresidue` tests candidates in a fixed order under `NonresidueSearchLimits`, with no
randomness. The order is chosen so that each list provably contains a nonsquare and finds one fast:

- For odd d the candidates are the constants 1, 2, …: c^((q−1)/2) = Legendre(c)^(1+p+…+pᵈ⁻¹) and the
  exponent is odd, so a nonresidue of F_p stays a nonresidue.
- For even d every constant is a square, so a plain digit-by-digit enumeration would burn the first
  p−1 candidates on constants. The candidates are instead the monic polynomials of degree 1 … d−1.
  Any nonsquare a = c·m with monic m forces m to be a nonsquare, so the list is complete.

The three outcomes are separate types of fact. `Found` carries a verified nonsquare;
`Absent` is a certificate (characteristic 2); `BudgetExhausted` is a resource outcome that says
nothing about the field. `sqrt` reports the same distinction as `NotSquare` (a proof, by the Euler
criterion), `Root` (verified by squaring) and `ResourceExhausted` (the search ran dry, not a claim
about the input). The tests cross-check every element of fields up to 289 elements against
brute-force squaring, and exercise q−1 with 2-adic valuation 62 using the norm criterion.

---

## Further Reading

1. **Cohen, H. (1993).** *A Course in Computational Algebraic Number Theory.* Springer GTM 138.
   The standard reference for algorithmic number theory: polynomial arithmetic, resultants,
   discriminants, Dedekind's theorem, and the Round 2 algorithm (§6.1).

2. **Crandall, R., and Pomerance, C. (2005).** *Prime Numbers: A Computational Perspective.*
   2nd ed. Springer. Chapter 6 covers the Number Field Sieve, including the algebraic side of the
   sieve and the role of norms and ideals.

3. **Lenstra, A. K., and Lenstra, H. W. Jr. (eds.) (1993).** *The Development of the Number Field
   Sieve.* Springer LNM 1554. The original papers on GNFS, including Buhler, Lenstra, and Pomerance
   on the algorithm and Lenstra on the algebraic side.

4. **Pohst, M., and Zassenhaus, H. (1989).** *Algorithmic Algebraic Number Theory.* Cambridge
   University Press. The Round 2 / HNF algorithm for computing the maximal order is described in
   detail here.

5. **Neukirch, J. (1999).** *Algebraic Number Theory.* Springer. The mathematical foundation:
   Dedekind domains, ideal factorisation, discriminants, and ramification theory.
