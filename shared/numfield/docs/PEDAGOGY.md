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
m) and Norm_{K/ℚ}(a + bα) (the algebraic norm) are B-smooth — divisible only by primes up to a
bound B. When enough such pairs are collected, linear algebra over 𝔽₂ finds candidate dependencies;
a verified square root of one gives a congruence of squares mod N, whose gcd with N is a
non-trivial factor unless the congruence is trivial.

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
the Mul canonicalisation contract (C-NF): it is what makes `PartialEq` for `NumberFieldElement`
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

This sign tripped up the G.A.1a KAT specification: the initial spec stated Norm(α − 1) = −1 for
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

The two-element form (p, α − r) represents exactly the prime ideals of residue degree 1 above p.
It cannot represent:
- Prime ideals of residue degree > 1 (arising from irreducible factors of f mod p of degree > 1).
- Non-principal ideals in general (though all prime ideals of residue degree 1 are principal in
  ℤ[α] when ℤ[α] is a PID).
- Products of ideals above the same prime p (the CRT formula requires p₁ ≠ p₂).

For a complete ideal-theoretic treatment, one would use the **Hermite Normal Form (HNF) basis**
representation, which can represent any ideal as a ℤ-module. This is noted as future work; the
two-element form is sufficient for the NFS factor base at toy scale.

---

## 8. Dedekind Factorisation

### Dedekind's theorem

Given a monic polynomial f ∈ ℤ[x] and a rational prime p with p ∤ disc(f), **Dedekind's theorem**
states that the prime ideals of ℤ[α] above p correspond bijectively to the irreducible factors of
f mod p. Specifically:

- If f ≡ ∏ᵢ gᵢ(x)^{eᵢ} (mod p) with gᵢ irreducible over 𝔽_p, then the prime ideals above p are
  𝔭ᵢ = (p, gᵢ(α)).
- For a linear factor gᵢ = x − r, the prime ideal is 𝔭ᵢ = (p, α − r) in two-element form.
- The ramification index eᵢ and residue degree fᵢ = deg(gᵢ) satisfy ∑ eᵢ fᵢ = d = deg(f).

The condition p ∤ disc(f) ensures that f is squarefree mod p (all eᵢ = 1), so the factorisation
is unramified.

### The linear-factor case

The implementation in `dedekind.rs` handles the **linear-factor case**: it finds all roots r ∈
{0, …, p−1} with f(r) ≡ 0 (mod p) by brute-force evaluation, and returns the corresponding prime
ideals (p, α − r):

```rust
pub fn dedekind_factor<'a>(field: &'a NumberField, p: &BigInt) -> Vec<Ideal<'a>> {
    let f = &field.f;
    let mut roots = Vec::new();
    for r in 0..p_to_usize(p) {
        let r_big = BigInt::from(r);
        if mod_reduce(&f.eval(&r_big), p).is_zero() {
            roots.push(r_big);
        }
    }
    // ...
}
```

This is correct and sufficient for the NFS factor base, which consists only of prime ideals of
residue degree 1. Irreducible factors of f mod p of degree > 1 contribute prime ideals of higher
residue degree, which are not included in the factor base at toy scale.

### The inert-prime sentinel convention

If f has no roots mod p — that is, all irreducible factors of f mod p have degree > 1 — then (p)
is **inert** in ℤ[α] (or factors only into higher-degree prime ideals). In this case,
`dedekind_factor` returns a single **sentinel ideal** (p, α − 0):

```rust
if roots.is_empty() {
    vec![Ideal::new(field, p.clone(), BigInt::zero())]
} else {
    roots.into_iter().map(|r| Ideal::new(field, p.clone(), r)).collect()
}
```

The sentinel convention avoids an ambiguous empty `Vec` return. Callers that need to distinguish
inert primes from split primes should check whether the returned ideal's r satisfies f(r) ≡ 0
(mod p). The downstream implication for G.C factor-base construction is that inert primes are
excluded from the algebraic factor base: the sentinel signals "no linear-factor ideals above p."

---

## 9. Bad Primes and the Dedekind Criterion

### What makes a prime "bad"

A prime p is **bad** for f (equivalently, **ramified** in the discriminant) if p | disc(f). The
discriminant of a monic polynomial f of degree d is:

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

For bad primes, Dedekind's theorem does not apply directly: f may be non-squarefree mod p (some
eᵢ > 1), and the ring ℤ[α] may not be the full ring of integers ℤ_K at p.

### The Dedekind criterion

The **Dedekind criterion** (also called the index criterion) tests whether p divides the index
[ℤ_K : ℤ[α]] — that is, whether ℤ[α] is strictly smaller than the full ring of integers at p.
The algorithm is:

1. Compute g = squarefree part of f mod p: g = f / gcd(f, f') mod p.
2. Compute h = f / g mod p.
3. Compute t = (g·h − f) / p mod p. (The numerator is divisible by p coefficient-wise because
   g·h ≡ f mod p by construction.)
4. Compute T = gcd(g, gcd(h, t)) mod p.
5. If T = 1 in 𝔽_p[x] (constant polynomial with value 1 mod p), then p ∤ [ℤ_K : ℤ[α]] and
   Dedekind's theorem applies normally. If T ≠ 1, then p | [ℤ_K : ℤ[α]] and ℤ[α] is not the
   maximal order at p.

The GCDs in steps 4 and 5 are computed via `subresultant_gcd` from `resultant.rs`, followed by
reduction mod p.

### The ℤ[√2] example

The polynomial f = x² − 2 has discriminant 8, so p = 2 is a bad prime. However, the Dedekind
criterion gives T = 1 at p = 2:

- f mod 2 = x² (since −2 ≡ 0 mod 2).
- f' mod 2 = 0 (since 2x ≡ 0 mod 2).
- gcd(x², 0) = x², so g = x²/x² = 1, h = x²/1 = x².
- g·h − f_mod = x² − x² = 0, so t = 0.
- T = gcd(1, gcd(x², 0)) = gcd(1, x²) = 1.

The criterion returns false (T = 1), meaning p = 2 does not divide [ℤ_K : ℤ[√2]]. This is
correct: ℤ[√2] is already the full ring of integers of ℚ(√2) — it is a PID, and the prime 2 is
totally ramified (2 = −(√2)²) but ℤ[√2] is already maximal. The prime 2 is bad (it divides the
discriminant) but the index is not divisible by 2.

This distinction — **bad prime** (p | disc) vs. **index divisible** (p | [ℤ_K : ℤ[α]]) — is
pedagogically important. A bad prime is a necessary but not sufficient condition for index
divisibility.

### The `DedekindResult` struct

The extended function `dedekind_factor_extended` returns a `DedekindResult` that carries both the
prime ideals and the bad-prime metadata:

```rust
pub struct DedekindResult<'a> {
    pub ideals: Vec<Ideal<'a>>,
    pub is_bad_prime: bool,
    pub index_divisible: bool,
}
```

The `index_divisible` flag is always false for good primes (p ∤ disc implies p ∤ [ℤ_K : ℤ[α]]
automatically). For bad primes, it reflects the Dedekind criterion result.

### The Round 2 / HNF algorithm

When `index_divisible` is true — that is, when ℤ[α] is strictly smaller than ℤ_K at p — the
correct factorisation of (p) requires the **Round 2 algorithm** (also called the HNF-basis
algorithm or Zassenhaus's algorithm). This algorithm computes a ℤ-basis for the full ring of
integers ℤ_K at p, expressed as an HNF (Hermite Normal Form) matrix, and then factors (p) using
this basis.

This algorithm is described but not implemented in this crate. The reason is engineering scale, not
mathematical omission: the Round 2 algorithm requires HNF computation over ℤ-modules, which is a
substantial implementation effort orthogonal to the core GNFS substrate. For the toy-scale
polynomials used in this project, index-divisible bad primes are rare and can be excluded from the
factor base without loss of correctness at demonstration fidelity.

---

## 10. Cross-Track Implications

### How the contracts are consumed downstream

The contracts frozen in G.A are consumed by multiple downstream sub-tracks:

| Contract | Frozen at | Consumed by |
|----------|-----------|-------------|
| C-NF | G.A.1a | G.A.1b, G.A.2, G.A.3, G.B, G.C, D.A, E.D |
| C-Ideal | G.A.1b | G.A.3, G.C, D.A |
| C-Res | G.A.2 | G.A.1a (norm), G.B, D.A |
| C-Dedekind | G.A.3 | G.C (factor-base construction) |

**G.B (polynomial selection)** uses `IntPoly` and the resultant to compute discriminants and
polynomial scores. The defining polynomial f is selected to have small coefficients and a large
number of roots mod small primes — both properties that depend on the resultant and discriminant
machinery.

**G.C (sieving)** is the primary consumer of the algebraic substrate. For each sieve pair (a, b),
it evaluates `NumberFieldElement::norm()` to obtain the algebraic norm, then factors it over the
algebraic factor base using the prime ideals from `dedekind_factor`. The inert-prime sentinel
convention directly affects which primes are included in the factor base: primes for which
`dedekind_factor` returns only the sentinel are excluded.

**D.A (NFS-DL)** uses the same algebraic substrate as G.C, adapted for the discrete logarithm
setting. The norm and ideal machinery are identical; the difference is in the linear algebra step.

**E.B (pairings)** does **not** use `NumberField` or `NumberFieldElement` for extension-field
arithmetic.  Pairing target fields are `F_{p^k}` — finite fields of **characteristic p**, not
char-0 `ℚ[x]/(f)` number fields in the sense of this crate.  The degree-2 extension `F_{47^2}`
used by the toy Weil/Tate fixture (and the degree-12 / degree-6 extensions of BN/BLS12 curves at
crypto scale) are `F_{p^k}`, not `ℚ`-extensions.  The char-0 `NumberFieldElement` arithmetic is
**not** directly reusable: the coefficient field differs (`F_p` vs `ℚ`).  The real char-p
substrate is `rho/src/pairing/fpext.rs` (`FpExt<F>`, a polynomial quotient `F_p[u]/(m(u))`
built over the `shared::field` `Fp<4>` prime field).  (Note: E.D is p-adic arithmetic, not
pairings.)

### The C1 resolution: `Uint<4>` sufficiency

A key question deferred from Phase α was whether the `shared::numth` crate's `Uint<4>` (256-bit
integers) is sufficient for GNFS norms. The G.A.1a session resolved this:

For a toy 60–80 bit semiprime N with factor base bound B ≤ 10⁴:
- Polynomial degree d = 3–5, coefficients C ~ N^{1/d} ~ 2^{20–27}.
- Sieve bound M ~ 10⁴ = 2^{13}.
- Norm bound: |Norm| ≤ (d+1)^{(d+1)/2} · (M · C)^d ≈ 2^{120–150}.

**256 bits suffices.** The `shared::numth` crate remains at `Uint<4>` with no widening needed.
If a future sub-track targets larger N (e.g., 128-bit semiprimes), the norm bound would exceed 256
bits and `Uint<4>` would need to be widened — but the `Fp` trait was designed at α.2 to support
this widening without breaking downstream consumers.

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
| `kat_resultant_3_norm_consistency` | `resultant_kat.rs` | Res(x²−2, x+1) = −1; cross-check with G.A.1a norm |
| `kat_gcd_1_shared_factor` | `resultant_kat.rs` | gcd(x²−1, x−1) ∝ x−1; subresultant PRS finds shared factor |
| `kat_gcd_2_coprime_polynomials` | `resultant_kat.rs` | gcd(x²−2, x²−3) is constant; coprime polynomials have trivial GCD |
| `kat1_inert_prime` | `dedekind_kat.rs` | p=3 is inert in ℤ[√2]; sentinel ideal (3, α−0) returned |
| `kat2_split_prime` | `dedekind_kat.rs` | p=7 splits in ℤ[√2]; ideals (7,α−3) and (7,α−4) |
| `kat3_cubic_partial_split` | `dedekind_kat.rs` | p=5 partially splits in ℤ[α] for f=x³−x−1; linear factor at r=2 |
| `kat4_norm_product` | `dedekind_kat.rs` | ∏ N(𝔭ᵢ) = p^d = 49 for totally split p=7 in degree-2 field |
| `kat5_discriminant_quadratic` | `dedekind_kat.rs` | disc(x²−2) = 8; discriminant via resultant |
| `kat6_discriminant_cubic` | `dedekind_kat.rs` | disc(x³−x−1) = −23; standard reference value |
| `kat7_is_bad_prime` | `dedekind_kat.rs` | 2 | disc(x²−2) = 8 (bad); 3 ∤ 8 (good) |
| `kat8_bad_prime_extended` | `dedekind_kat.rs` | p=2 is bad for x²−2 but index_divisible=false; ℤ[√2] is maximal order |
| `kat9_good_prime_extended` | `dedekind_kat.rs` | p=7 is good for x²−2; extended interface agrees with basic interface |
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
