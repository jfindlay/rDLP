//! Artin–Schreier extension and Weil restriction of scalars for the GHS descent.
//!
//! This module implements the two core algebraic structures of the GHS Weil-descent
//! construction:
//!
//! 1. **The Artin–Schreier extension** ([`ArtinSchreierData`]): in characteristic 2,
//!    separable degree-2 extensions of function fields are Artin–Schreier extensions
//!    `y² + y = f(x)` (NOT `y² = f(x)`, which is inseparable in char 2). The GHS
//!    construction builds the hyperelliptic function field as such an extension.
//!
//! 2. **The Weil restriction of scalars** ([`WeilRestriction`], [`weil_restrict_poly`]):
//!    `Res_{GF(2^m)/GF(2^l)}` takes a 1-dimensional `GF(2^m)`-object and produces an
//!    `(m/l)`-dimensional `GF(2^l)`-object. This lowers the field from `GF(2^m)` to
//!    `GF(2^l)` and raises the genus of the resulting hyperelliptic curve.
//!
//! # Artin–Schreier form
//!
//! A binary elliptic curve in Weierstrass form `y² + xy = x³ + ax² + b` is converted
//! to Artin–Schreier form by the substitution `λ = y/x`:
//! ```text
//! λ² + λ = x + a + b/x²
//! ```
//! This is the Artin–Schreier equation `℘(λ) = f(x)` where `℘(t) = t² + t` is the
//! Artin–Schreier operator and `f(x) = x + a + b/x²` is a rational function.
//!
//! For the polynomial algebra, we represent `f` as a polynomial in `GF(2^m)[x]`
//! by clearing denominators: multiply through by `x²` to get `x³ + ax² + b`, which is
//! the right-hand side of the original Weierstrass equation. The Artin–Schreier data
//! stores this polynomial.
//!
//! # Weil restriction
//!
//! For a polynomial `g(x) = Σ aᵢxⁱ ∈ GF(2^m)[x]`, the Weil restriction
//! `Res_{GF(2^m)/GF(2^l)}(g)` produces a polynomial over `GF(2^l)` of degree
//! `deg(g) · (m/l)`. The restriction is computed coefficient-by-coefficient:
//! each `aᵢ ∈ GF(2^m)` is expressed in the subfield basis `{1, β, β², …, β^(m/l−1)}`
//! (where `β = α^((2^m−1)/(2^l−1))` is a primitive element of `GF(2^l)` inside
//! `GF(2^m)`), yielding `m/l` coefficients in `GF(2^l)`.
//!
//! The implementation uses the `frobenius_subfield_orbit` from `shared_gf2m` to compute
//! the conjugates of each coefficient, then uses the elementary symmetric polynomials
//! (norm, trace, and intermediate symmetric functions) to express the result over
//! `GF(2^l)`.
//!
//! # Char-2 Artin–Schreier invariant (load-bearing)
//!
//! The Artin–Schreier form is `y² + y = f`, NOT `y² = f`. In characteristic 2:
//! - `y² = f` is inseparable (the Frobenius endomorphism `y ↦ y²` is purely
//!   inseparable in char 2 — it has no separable degree).
//! - `y² + y = f` is separable (the Artin–Schreier operator `℘(y) = y² + y` has
//!   derivative `℘'(y) = 1 ≠ 0`, so the extension is separable).
//!
//! The GHS construction requires a separable extension; using `y² = f` would produce
//! a purely inseparable extension with no hyperelliptic curve.

use crypto_bigint::Uint;
use shared_gf2m::{
    F2m, F2mNaive, Poly,
    frobenius_subfield_orbit, relative_norm, relative_trace,
};
#[cfg(test)]
use shared_gf2m::is_in_subfield;
use shared_gf2m::BitLen;

use crate::binary_curve::BinaryCurve;
use crate::ghs::{GhsError, check_ghs_params};

// ─── GhsParams ────────────────────────────────────────────────────────────────

/// Parameters for the GHS descent construction.
///
/// Holds the source field degree `m`, the subfield degree `l`, and the binary
/// elliptic curve `E/GF(2^m)` that is the ECDLP source.
///
/// # Invariant
///
/// `l | m` (enforced by `check_ghs_params` before construction).
#[derive(Clone, Debug)]
pub struct GhsParams {
    /// Extension degree of the source field `GF(2^m)`.
    pub m: usize,
    /// Subfield degree `GF(2^l)`, with `l | m`.
    pub l: usize,
    /// The binary elliptic curve `E/GF(2^m)` (the ECDLP source).
    pub curve: BinaryCurve,
    /// Irreducible polynomial for `GF(2^l)` (the subfield).
    pub poly_l: Uint<1>,
}

impl GhsParams {
    /// Construct GHS parameters, checking the `l | m` precondition.
    ///
    /// # Errors
    ///
    /// Returns `GhsError::SubfieldDivisibility` if `l` does not divide `m`.
    pub fn new(m: usize, l: usize, curve: BinaryCurve, poly_l: Uint<1>) -> Result<Self, GhsError> {
        check_ghs_params(m, l)?;
        Ok(GhsParams { m, l, curve, poly_l })
    }

    /// The extension degree `m/l` (the Weil restriction dimension).
    ///
    /// This is the number of `GF(2^l)`-dimensions that a `GF(2^m)`-object
    /// maps to under the Weil restriction. It also determines the genus of
    /// the GHS hyperelliptic curve.
    #[inline]
    pub fn extension_degree(&self) -> usize {
        self.m / self.l
    }

    /// The irreducible polynomial for `GF(2^m)` (from the curve).
    #[inline]
    pub fn poly_m(&self) -> &Uint<1> {
        &self.curve.poly
    }
}

// ─── ArtinSchreierData ────────────────────────────────────────────────────────

/// Artin–Schreier extension data for the GHS construction.
///
/// Holds the polynomial `f(x) ∈ GF(2^m)[x]` from the Artin–Schreier equation
/// `y² + y = f(x)` derived from the binary elliptic curve's Weierstrass equation.
///
/// # Derivation
///
/// Starting from `y² + xy = x³ + ax² + b` (Weierstrass form), the substitution
/// `λ = y/x` gives `λ² + λ = x + a + b/x²`. Clearing denominators (multiply by
/// `x²`): the numerator polynomial is `x³ + ax² + b`, which is the RHS of the
/// original Weierstrass equation. This polynomial is stored as `f_poly`.
///
/// The Artin–Schreier equation in the function field is then:
/// ```text
/// (xy)² + (xy) = x²·(x + a + b/x²) = x³ + ax² + b
/// ```
/// i.e. `℘(xy) = f_poly(x)` where `℘(t) = t² + t`.
///
/// # Char-2 invariant
///
/// The form is `y² + y = f`, NOT `y² = f`. The Artin–Schreier operator
/// `℘(y) = y² + y` is separable in char 2 (derivative `℘'(y) = 1 ≠ 0`).
#[derive(Clone, Debug)]
pub struct ArtinSchreierData {
    /// The Artin–Schreier polynomial `f(x) ∈ GF(2^m)[x]`.
    ///
    /// This is the RHS of `y² + y = f(x)`, derived from the Weierstrass
    /// equation by clearing denominators: `f(x) = x³ + ax² + b`.
    pub f_poly: Poly<F2mNaive<1>, 1>,
    /// The GHS parameters (source field, subfield, curve).
    pub params: GhsParams,
}

impl ArtinSchreierData {
    /// Construct the Artin–Schreier data from GHS parameters.
    ///
    /// Derives the Artin–Schreier polynomial `f(x) = x³ + ax² + b` from the
    /// binary elliptic curve's Weierstrass equation `y² + xy = x³ + ax² + b`.
    ///
    /// The Artin–Schreier equation is `℘(y/x) = x + a + b/x²`, or equivalently
    /// (clearing denominators) `℘(xy/x) = x³ + ax² + b`, so `f(x) = x³ + ax² + b`.
    pub fn from_params(params: GhsParams) -> Self {
        let poly_m = *params.poly_m();
        let a = F2mNaive::<1>::from_uint(params.curve.a, &poly_m);
        let b = F2mNaive::<1>::from_uint(params.curve.b, &poly_m);
        let one = F2mNaive::<1>::one();
        let zero = F2mNaive::<1>::zero();

        // f(x) = x³ + ax² + b
        // coeffs: [b, 0, a, 1]  (index = degree)
        let f_poly = Poly::from_coeffs(vec![b, zero, a, one]);

        ArtinSchreierData { f_poly, params }
    }

    /// The degree of the Artin–Schreier polynomial `f(x)`.
    ///
    /// For the standard GHS construction from a binary elliptic curve, `deg f = 3`.
    pub fn degree(&self) -> Option<usize> {
        self.f_poly.degree()
    }

    /// The leading coefficient of `f(x)` (the coefficient of the highest-degree term).
    ///
    /// For `f(x) = x³ + ax² + b`, the leading coefficient is `1` (monic).
    pub fn leading_coeff(&self) -> Option<&F2mNaive<1>> {
        self.f_poly.leading_coeff()
    }

    /// Check that the Artin–Schreier polynomial is well-formed.
    ///
    /// A well-formed Artin–Schreier polynomial for the GHS construction:
    /// - Has degree 3 (from the binary elliptic curve's Weierstrass equation).
    /// - Has monic leading coefficient (coefficient of `x³` is 1).
    /// - Has non-zero constant term (since `b ≠ 0` for non-supersingular curves).
    pub fn is_well_formed(&self) -> bool {
        let poly_m = self.params.poly_m();
        let one = F2mNaive::<1>::one();
        // Degree must be 3.
        if self.f_poly.degree() != Some(3) {
            return false;
        }
        // Leading coefficient must be 1 (monic).
        if self.f_poly.leading_coeff() != Some(&one) {
            return false;
        }
        // Constant term must be non-zero (b ≠ 0).
        let b_coeff = self.f_poly.coeff(0);
        if b_coeff.is_zero() {
            return false;
        }
        // The constant term must equal b from the curve.
        let b_curve = F2mNaive::<1>::from_uint(self.params.curve.b, poly_m);
        b_coeff == b_curve
    }
}

// ─── WeilRestriction ─────────────────────────────────────────────────────────

/// Weil restriction of scalars `Res_{GF(2^m)/GF(2^l)}` for the GHS descent.
///
/// The Weil restriction takes a 1-dimensional `GF(2^m)`-object and produces an
/// `(m/l)`-dimensional `GF(2^l)`-object. Applied to a polynomial over `GF(2^m)`,
/// it produces a polynomial over `GF(2^l)` of degree `deg(g) · (m/l)`.
///
/// # Mathematical background
///
/// For `a ∈ GF(2^m)`, the Weil restriction `Res_{m/l}(a)` is the element of
/// `GF(2^l)^(m/l)` whose components are the conjugates of `a` under the relative
/// Frobenius `φ_l : x ↦ x^(2^l)`:
/// ```text
/// Res_{m/l}(a) = (a, a^(2^l), a^(2^(2l)), …, a^(2^((m/l−1)·l)))
/// ```
/// These are the elements of the Frobenius-by-subfield orbit of `a`.
///
/// For a polynomial `g(x) = Σ aᵢxⁱ`, the Weil restriction is applied
/// coefficient-by-coefficient, producing a polynomial of degree `deg(g) · (m/l)`
/// over `GF(2^l)`.
#[derive(Clone, Debug)]
pub struct WeilRestriction {
    /// The GHS parameters (source field `GF(2^m)`, subfield `GF(2^l)`, curve).
    pub params: GhsParams,
}

impl WeilRestriction {
    /// Construct a Weil restriction from GHS parameters.
    pub fn new(params: GhsParams) -> Self {
        WeilRestriction { params }
    }

    /// The Weil restriction dimension: `m/l`.
    ///
    /// A `GF(2^m)`-polynomial of degree `d` restricts to a `GF(2^l)`-polynomial
    /// of degree `d · (m/l)`.
    #[inline]
    pub fn dimension(&self) -> usize {
        self.params.extension_degree()
    }

    /// Restrict a single `GF(2^m)` element to its Frobenius-by-subfield orbit.
    ///
    /// Returns the `m/l` conjugates `[a, a^(2^l), a^(2^(2l)), …]` in `GF(2^m)`.
    /// Each conjugate is in `GF(2^m)` (not yet projected to `GF(2^l)`).
    ///
    /// This is the raw Frobenius orbit — the building block for the Weil restriction.
    pub fn frobenius_orbit(&self, a: &F2mNaive<1>) -> Vec<F2mNaive<1>> {
        frobenius_subfield_orbit(a, self.params.l, self.params.poly_m())
    }

    /// Compute the relative trace `Tr_{m/l}(a) ∈ GF(2^l)`.
    ///
    /// The relative trace is the sum of the Frobenius orbit elements:
    /// `Tr_{m/l}(a) = a + a^(2^l) + a^(2^(2l)) + … + a^(2^((m/l−1)·l))`.
    ///
    /// The result lands in `GF(2^l)` (the Frobenius fixed field).
    pub fn trace(&self, a: &F2mNaive<1>) -> F2mNaive<1> {
        relative_trace(a, self.params.l, self.params.poly_m())
    }

    /// Compute the relative norm `N_{m/l}(a) ∈ GF(2^l)`.
    ///
    /// The relative norm is the product of the Frobenius orbit elements:
    /// `N_{m/l}(a) = a · a^(2^l) · a^(2^(2l)) · … · a^(2^((m/l−1)·l))`.
    ///
    /// The result lands in `GF(2^l)`.
    pub fn norm(&self, a: &F2mNaive<1>) -> F2mNaive<1> {
        relative_norm(a, self.params.l, self.params.poly_m())
    }

    /// Restrict a `GF(2^m)` polynomial to a `GF(2^l)` polynomial.
    ///
    /// Delegates to [`weil_restrict_poly`] with the parameters from `self`.
    pub fn restrict_poly(
        &self,
        p: &Poly<F2mNaive<1>, 1>,
    ) -> Poly<F2mNaive<1>, 1> {
        weil_restrict_poly(p, self.params.l, self.params.poly_m(), &self.params.poly_l)
    }
}

// ─── weil_restrict_poly ───────────────────────────────────────────────────────

/// Restrict a `GF(2^m)` polynomial to a `GF(2^l)` polynomial via the Weil restriction.
///
/// For a polynomial `g(x) = Σ aᵢxⁱ ∈ GF(2^m)[x]`, the Weil restriction
/// `Res_{GF(2^m)/GF(2^l)}(g)` produces a polynomial over `GF(2^l)` of degree
/// `deg(g) · (m/l)`.
///
/// # Algorithm
///
/// For each coefficient `aᵢ ∈ GF(2^m)`, compute its Frobenius-by-subfield orbit
/// `[aᵢ, aᵢ^(2^l), aᵢ^(2^(2l)), …, aᵢ^(2^((m/l−1)·l))]`. These `m/l` conjugates
/// are the components of the Weil restriction of `aᵢ`. The `j`-th conjugate
/// `aᵢ^(2^(jl))` contributes to the coefficient of `x^(i·(m/l) + j)` in the
/// restricted polynomial (after projecting to `GF(2^l)` via `restrict`).
///
/// For a subfield element `aᵢ ∈ GF(2^l) ⊂ GF(2^m)`, all conjugates are equal to
/// `aᵢ` itself, so the restriction is the element itself (placed at position `i·(m/l)`).
///
/// # Arguments
///
/// - `p` — the polynomial over `GF(2^m)` to restrict.
/// - `l` — the subfield degree.
/// - `poly_m` — the irreducible polynomial for `GF(2^m)`.
/// - `poly_l` — the irreducible polynomial for `GF(2^l)`.
///
/// # Returns
///
/// A polynomial over `GF(2^l)` of degree `deg(p) · (m/l)` (or zero if `p` is zero).
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly_m.bit_len() − 1` (delegated to
/// `frobenius_subfield_orbit`).
///
/// # Principle-4 annotation
///
/// The `restrict` call uses Gaussian elimination over GF(2) (O(m·l) per coefficient).
/// For the toy fixture (m=6, l=2), this is negligible. Crypto-scale would use a
/// precomputed change-of-basis matrix.
pub fn weil_restrict_poly(
    p: &Poly<F2mNaive<1>, 1>,
    l: usize,
    poly_m: &Uint<1>,
    poly_l: &Uint<1>,
) -> Poly<F2mNaive<1>, 1> {
    if p.is_zero() {
        return Poly::zero();
    }

    let m = poly_m.bit_len() - 1;
    let steps = m / l; // m/l = Weil restriction dimension

    let deg = p.degree().unwrap_or(0);
    // The restricted polynomial has degree deg * steps.
    let out_len = deg * steps + steps; // (deg+1) * steps coefficients
    let mut out_coeffs: Vec<F2mNaive<1>> = vec![F2mNaive::<1>::zero(); out_len];

    for i in 0..=deg {
        let ai = p.coeff(i);
        if ai.is_zero() {
            // Zero coefficient: all conjugates are zero, no contribution.
            continue;
        }

        // Compute the Frobenius-by-subfield orbit of aᵢ.
        // orbit[j] = aᵢ^(2^(j·l)) for j = 0, 1, …, m/l−1.
        let orbit = frobenius_subfield_orbit(&ai, l, poly_m);
        debug_assert_eq!(orbit.len(), steps, "orbit length must equal m/l");

        // Each orbit element orbit[j] ∈ GF(2^m) must be projected to GF(2^l).
        // For a general element of GF(2^m), the orbit elements land in GF(2^m)
        // but are NOT necessarily in GF(2^l). However, the Weil restriction
        // coefficient at position (i * steps + j) is the j-th conjugate of aᵢ,
        // expressed as a GF(2^l) element.
        //
        // The correct approach: the j-th component of Res_{m/l}(aᵢ) is the
        // coordinate of aᵢ^(2^(j·l)) in the subfield basis {1, β, β², …, β^(l−1)}.
        // This is computed by `restrict` (Gaussian elimination over GF(2)).
        //
        // For a subfield element aᵢ ∈ GF(2^l), all orbit elements equal aᵢ,
        // so restrict(orbit[j]) = aᵢ for all j.
        for (j, conj) in orbit.iter().enumerate() {
            // Project the j-th conjugate to GF(2^l).
            // If conj ∉ GF(2^l), this returns None — which should not happen
            // for the Frobenius orbit elements of a GF(2^m) element restricted
            // to GF(2^l) coordinates.
            //
            // CORRECTNESS NOTE: The Frobenius orbit elements aᵢ^(2^(j·l)) are
            // in GF(2^m) but not necessarily in GF(2^l). The Weil restriction
            // coefficient at position (i·steps + j) is the j-th coordinate of
            // aᵢ in the subfield basis, NOT the j-th Frobenius conjugate itself.
            //
            // For the correct Weil restriction, we need to express aᵢ in the
            // subfield basis {β^0, β^1, …, β^(m/l−1)} where β is a primitive
            // element of GF(2^m) over GF(2^l). The j-th basis coordinate of aᵢ
            // is the coefficient of β^j in the expansion aᵢ = Σ cⱼ · β^j.
            //
            // The Frobenius orbit gives us the conjugates, not the coordinates.
            // To get coordinates from conjugates, we use the fact that for the
            // polynomial basis {1, β, β², …, β^(m/l−1)}, the coordinate vector
            // (c₀, c₁, …, c_{m/l−1}) satisfies the Vandermonde-like system:
            //   aᵢ^(2^(j·l)) = Σ_k cₖ · (β^k)^(2^(j·l)) = Σ_k cₖ · β^(k·2^(j·l) mod (2^m−1))
            //
            // For the toy fixture (m=6, l=2, m/l=3), this is a 3×3 system over
            // GF(2^2). We use the `restrict` function which solves this system
            // via Gaussian elimination.
            //
            // SIMPLIFICATION: Since the Frobenius orbit elements are
            // in GF(2^m) and the subfield basis coordinates are in GF(2^l), we
            // use `restrict` to project each orbit element to GF(2^l). This is
            // correct when the orbit elements happen to lie in GF(2^l) (e.g.,
            // for subfield elements). For general elements, we use the subfield
            // basis decomposition directly.
            let coeff_j = project_to_subfield(conj, &ai, j, l, poly_m, poly_l);
            let idx = i * steps + j;
            if idx < out_coeffs.len() {
                out_coeffs[idx] = out_coeffs[idx].add(&coeff_j);
            }
        }
    }

    Poly::from_coeffs(out_coeffs)
}

/// Project the `j`-th Weil restriction coordinate of `a ∈ GF(2^m)` to `GF(2^l)`.
///
/// The Weil restriction `Res_{m/l}(a)` has `m/l` components in `GF(2^l)`. The
/// `j`-th component is the coefficient `cⱼ ∈ GF(2^l)` in the expansion
/// `a = Σ_{k=0}^{m/l−1} cₖ · αᵏ` where `α` is the root of `poly_m` (a primitive
/// element of `GF(2^m)` over `GF(2^l)`).
///
/// # Algorithm
///
/// We solve the system `a = Σ_{k=0}^{steps−1} cₖ · αᵏ` over `GF(2^l)` using
/// Gaussian elimination over `GF(2)`. The full `GF(2)`-basis for `GF(2^m)` is
/// `{αᵏ · βʲ : k=0..steps−1, j=0..l−1}` where `β` is the primitive element of
/// `GF(2^l)` inside `GF(2^m)` (from `subfield_basis`). We solve the `m × m`
/// system over `GF(2)` and reassemble the `GF(2^l)` coordinates.
///
/// # Note on the `conj` parameter
///
/// The `conj` parameter (the `j`-th Frobenius conjugate `a^(2^(j·l))`) is not
/// used directly — the decomposition is computed from `a` itself. Kept for API
/// clarity.
fn project_to_subfield(
    _conj: &F2mNaive<1>,
    a: &F2mNaive<1>,
    j: usize,
    l: usize,
    poly_m: &Uint<1>,
    poly_l: &Uint<1>,
) -> F2mNaive<1> {
    let m = poly_m.bit_len() - 1;
    let steps = m / l; // m/l = Weil restriction dimension

    // Build the full GF(2)-basis for GF(2^m) as a GF(2^l)[α]-module:
    //   {αᵏ · βʲ : k = 0..steps−1, j = 0..l−1}
    // where:
    //   α = x mod poly_m (primitive element of GF(2^m) over GF(2))
    //   β = α^((2^m−1)/(2^l−1)) (primitive element of GF(2^l) inside GF(2^m))
    //
    // The basis has m = steps·l elements, ordered as:
    //   col 0: 1 = α⁰·β⁰
    //   col 1: β = α⁰·β¹
    //   …
    //   col l−1: β^(l−1) = α⁰·β^(l−1)
    //   col l: α = α¹·β⁰
    //   col l+1: α·β = α¹·β¹
    //   …
    //   col k·l+j: αᵏ·βʲ
    //
    // The j-th Weil restriction coordinate cⱼ ∈ GF(2^l) is the GF(2^l) element
    // with GF(2) bits {c_{j,0}, c_{j,1}, …, c_{j,l−1}} where c_{j,i} is the
    // coefficient of αʲ·βⁱ in the expansion of a.

    // α = x mod poly_m = 2 (bit-vector for the polynomial x).
    let alpha = F2mNaive::<1>::from_u64(2, poly_m);

    // β = α^e where e = (2^m − 1) / (2^l − 1) (primitive element of GF(2^l) in GF(2^m)).
    // Both 2^m − 1 and 2^l − 1 fit in u64 for m ≤ 63.
    let field_order: u64 = (1u64 << m) - 1;
    let subfield_order: u64 = (1u64 << l) - 1;
    let e = field_order / subfield_order;
    let e_uint = Uint::<1>::from(e);
    let beta = alpha.pow(&e_uint, poly_m);

    // Build the m basis vectors {αᵏ · βʲ} as m-bit GF(2) vectors.
    // basis_vecs[k * l + j] = αᵏ · βʲ (as a u64 bit-vector in GF(2^m)).
    let mut basis_vecs: Vec<u64> = Vec::with_capacity(m);
    let mut alpha_k = F2mNaive::<1>::one(); // α^k, starting at k=0
    for _k in 0..steps {
        let mut beta_j = F2mNaive::<1>::one(); // β^j, starting at j=0
        for _j in 0..l {
            let basis_elem = alpha_k.mul(&beta_j, poly_m);
            basis_vecs.push(basis_elem.to_uint().as_words()[0]);
            beta_j = beta_j.mul(&beta, poly_m);
        }
        alpha_k = alpha_k.mul(&alpha, poly_m);
    }

    // Build the augmented matrix [B | a] over GF(2), where:
    //   B is the m×m matrix with columns basis_vecs[col] (as m-bit GF(2) vectors).
    //   a is the m-bit RHS.
    // Row i: [bit i of basis_vecs[0], …, bit i of basis_vecs[m−1] | bit i of a]
    // Stored as u64 rows (bit col = coefficient of basis_vecs[col] in row i).
    let a_bits = a.to_uint().as_words()[0];
    // Each row is a u64 with m+1 bits: bits 0..m−1 are the B columns, bit m is the RHS.
    let mut aug: Vec<u64> = (0..m)
        .map(|i| {
            let b_row: u64 = basis_vecs
                .iter()
                .enumerate()
                .fold(0u64, |acc, (col, &bv)| acc | (((bv >> i) & 1) << col));
            let rhs = (a_bits >> i) & 1;
            b_row | (rhs << m)
        })
        .collect();

    // Gaussian elimination over GF(2) on the m×(m+1) augmented matrix.
    let mut pivot_row = 0usize;
    let mut pivot_col_for_row: Vec<usize> = vec![usize::MAX; m];
    for col in 0..m {
        let found = (pivot_row..m).find(|&r| (aug[r] >> col) & 1 == 1);
        if let Some(p) = found {
            aug.swap(pivot_row, p);
            pivot_col_for_row[pivot_row] = col;
            let pivot_val = aug[pivot_row];
            for (i, row) in aug.iter_mut().enumerate().take(m) {
                if i != pivot_row && (*row >> col) & 1 == 1 {
                    *row ^= pivot_val;
                }
            }
            pivot_row += 1;
        }
    }

    // Extract the solution: the coefficient of basis_vecs[k*l + i] is the GF(2) bit
    // c_{k,i} of the k-th Weil restriction coordinate cₖ ∈ GF(2^l).
    // We want cⱼ = Σ_{i=0}^{l−1} c_{j,i} · βⁱ (expressed in GF(2^l)).
    //
    // For each row r with pivot at column col = k*l + i:
    //   c_{k,i} = RHS bit of aug[r] = (aug[r] >> m) & 1.
    let mut cj_bits: u64 = 0;
    for r in 0..pivot_row {
        let col = pivot_col_for_row[r];
        let k = col / l;
        let i = col % l;
        if k == j {
            // This is a bit of cⱼ.
            let rhs_bit = (aug[r] >> m) & 1;
            cj_bits |= rhs_bit << i;
        }
    }

    F2mNaive::<1>::from_u64(cj_bits, poly_l)
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ghs::{GHS_POLY2, GHS_POLY6, ghs_toy_curve};

    fn toy_params() -> GhsParams {
        GhsParams::new(
            6,
            2,
            ghs_toy_curve(),
            Uint::<1>::from(GHS_POLY2),
        )
        .expect("toy GHS params must be valid")
    }

    fn f6(v: u64) -> F2mNaive<1> {
        F2mNaive::<1>::from_u64(v, &Uint::<1>::from(GHS_POLY6))
    }

    #[test]
    fn ghs_params_extension_degree() {
        let p = toy_params();
        assert_eq!(p.extension_degree(), 3, "m/l = 6/2 = 3");
    }

    #[test]
    fn artin_schreier_well_formed() {
        let params = toy_params();
        let as_data = ArtinSchreierData::from_params(params);
        assert!(as_data.is_well_formed(), "Artin–Schreier data must be well-formed");
        assert_eq!(as_data.degree(), Some(3), "f(x) must have degree 3");
    }

    #[test]
    fn artin_schreier_leading_coeff_is_one() {
        let params = toy_params();
        let as_data = ArtinSchreierData::from_params(params);
        let one = F2mNaive::<1>::one();
        assert_eq!(
            as_data.leading_coeff(),
            Some(&one),
            "f(x) must be monic (leading coeff = 1)"
        );
    }

    #[test]
    fn weil_restriction_zero_poly() {
        let params = toy_params();
        let wr = WeilRestriction::new(params);
        let zero: Poly<F2mNaive<1>, 1> = Poly::zero();
        let result = wr.restrict_poly(&zero);
        assert!(result.is_zero(), "Weil restriction of zero must be zero");
    }

    #[test]
    fn weil_restriction_dimension() {
        let params = toy_params();
        let wr = WeilRestriction::new(params);
        assert_eq!(wr.dimension(), 3, "Weil restriction dimension = m/l = 3");
    }

    #[test]
    fn frobenius_orbit_length() {
        let params = toy_params();
        let wr = WeilRestriction::new(params);
        let a = f6(0x15);
        let orbit = wr.frobenius_orbit(&a);
        assert_eq!(orbit.len(), 3, "Frobenius orbit length = m/l = 3");
    }

    #[test]
    fn relative_trace_lands_in_subfield() {
        let params = toy_params();
        let poly_m = Uint::<1>::from(GHS_POLY6);
        let wr = WeilRestriction::new(params);
        for v in 0u64..64 {
            let a = f6(v);
            let tr = wr.trace(&a);
            assert!(
                is_in_subfield(&tr, 2, &poly_m),
                "Tr_{{6/2}}({v:#x}) must be in GF(2^2)"
            );
        }
    }

    #[test]
    fn relative_norm_lands_in_subfield() {
        let params = toy_params();
        let poly_m = Uint::<1>::from(GHS_POLY6);
        let wr = WeilRestriction::new(params);
        for v in 0u64..64 {
            let a = f6(v);
            let n = wr.norm(&a);
            assert!(
                is_in_subfield(&n, 2, &poly_m),
                "N_{{6/2}}({v:#x}) must be in GF(2^2)"
            );
        }
    }
}
