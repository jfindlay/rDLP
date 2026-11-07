//! Subfield substrate GF(2^l) ⊂ GF(2^m): embedding, relative trace/norm, Frobenius-by-subfield.
//!
//! This module implements the relative-field substrate the GHS/Weil descent stands on.
//! The subfield GF(2^l) is the Frobenius fixed field: `{a ∈ GF(2^m) : a^(2^l) = a}`.
//! All operations require `l | m`; this is asserted loudly.
//!
//! # Design: per-call poly threading
//!
//! All operations that depend on the field identity take `poly: &Uint<L>` as an explicit
//! parameter, mirroring the C-F2m idiom.  The field degree `m = poly.bit_len() − 1` is
//! recovered at runtime.
//!
//! # Relative trace and norm
//!
//! The relative trace `Tr_{m/l}` and norm `N_{m/l}` iterate the **l-th Frobenius power**
//! (steps of `l`), NOT the absolute Frobenius (steps of 1).  The absolute `trace` on the
//! `F2m` trait lands in GF(2); the relative trace lands in GF(2^l).  These are distinct maps.
//!
//! # GHS fixture
//!
//! The canonical GHS fixture is `m = 6, l = 2, m/l = 3` (odd).  Odd `m/l` keeps the
//! descent in the imaginary/ramified hyperelliptic model that the frozen C-HyperCurve handles.
//! Even `m/l` risks the real/split model.  The KATs assert non-degeneracy at this fixture.
//!
//! # Subfield basis
//!
//! The polynomial basis for GF(2^l) inside GF(2^m) is `{1, β, β², …, β^(l−1)}` where
//! `β = α^((2^m − 1)/(2^l − 1))` is a primitive element of GF(2^l) expressed in GF(2^m)
//! (here `α` is the root of `poly`, i.e. the primitive element of GF(2^m)).

use crate::BitLen;
use crypto_bigint::Uint;

use crate::naive::F2mNaive;
use crate::F2m;

// ── Precondition helper ───────────────────────────────────────────────────────

/// Assert that `l` divides `m`, panicking with a clear message if not.
///
/// GF(2^l) ⊆ GF(2^m) if and only if `l | m`.  All subfield operations require this.
#[inline]
fn assert_l_divides_m(l: usize, m: usize) {
    assert!(
        l > 0 && m % l == 0,
        "subfield: l must divide m (l = {l}, m = {m}, m % l = {}); \
         GF(2^l) ⊆ GF(2^m) iff l | m",
        m % l
    );
}

// ── Subfield membership ───────────────────────────────────────────────────────

/// Test whether `a ∈ GF(2^l)` inside GF(2^m).
///
/// GF(2^l) is the Frobenius fixed field: `a ∈ GF(2^l)` iff `a^(2^l) = a`.
/// This is checked by applying `frobenius` `l` times and comparing.
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly.bit_len() − 1`.
#[must_use]
pub fn is_in_subfield<const L: usize>(a: &F2mNaive<L>, l: usize, poly: &Uint<L>) -> bool
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    assert_l_divides_m(l, m);
    // Apply Frobenius l times: a → a^(2^l).
    let mut cur = a.clone();
    for _ in 0..l {
        cur = cur.frobenius(poly);
    }
    cur == *a
}

// ── Subfield embedding and restriction ───────────────────────────────────────

/// Embed a GF(2^l) element into GF(2^m).
///
/// Maps `a ∈ GF(2^l)` (expressed in the polynomial basis of GF(2^l) with root `α_l` of
/// `poly_l`) to its image in GF(2^m).  The embedding sends the primitive element `α_l` of
/// GF(2^l) to `β = α_m^((2^m−1)/(2^l−1))`, the primitive element of GF(2^l) inside GF(2^m).
///
/// Concretely: `embed(Σ aᵢ · α_l^i) = Σ aᵢ · β^i` where `β^i` are the subfield basis
/// vectors computed by `subfield_basis`.
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly_m.bit_len() − 1`.
/// Panics if `m > 63` (see `subfield_basis`).
#[must_use]
pub fn embed<const L: usize>(
    a: &F2mNaive<L>,
    poly_l: &Uint<L>,
    poly_m: &Uint<L>,
) -> F2mNaive<L>
where
    F2mNaive<L>: F2m<L>,
{
    let l = poly_l.bit_len() - 1; // degree of GF(2^l)
    let basis = subfield_basis(l, poly_m); // {1, β, β², …, β^(l−1)} in GF(2^m)
    let a_bits = a.to_uint();
    // Compute Σ aᵢ · basis[i] (XOR in char 2).
    let mut result = F2mNaive::<L>::zero();
    for (i, b) in basis.iter().enumerate().take(l) {
        if a_bits.bit_at(i) {
            result = result.add(b);
        }
    }
    result
}

/// Restrict a GF(2^m) element to GF(2^l), returning `Some` if it lies in the subfield.
///
/// Returns `Some(c)` where `c` is the GF(2^l) representation of `a` (i.e. the coefficients
/// in the polynomial basis `{1, β, β², …, β^(l−1)}` of GF(2^l) inside GF(2^m)), or `None`
/// if `a ∉ GF(2^l)`.
///
/// The returned element is the GF(2^l) bit-vector: bit `i` is the coefficient of `α_l^i`
/// where `α_l` is the primitive element of GF(2^l).  This is the inverse of `embed`.
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly_m.bit_len() − 1`.
/// Panics if `m > 63` (see `subfield_basis`).
#[must_use]
pub fn restrict<const L: usize>(
    a: &F2mNaive<L>,
    l: usize,
    poly_m: &Uint<L>,
) -> Option<F2mNaive<L>>
where
    F2mNaive<L>: F2m<L>,
{
    if !is_in_subfield(a, l, poly_m) {
        return None;
    }
    // Recover the GF(2^l) coordinates by solving a over the subfield basis.
    // The basis is {b₀=1, b₁=β, b₂=β², …, b_{l-1}=β^(l-1)} in GF(2^m).
    // We need to find bits c₀, c₁, …, c_{l-1} ∈ GF(2) such that
    //   a = Σ cᵢ · bᵢ   (in GF(2^m)).
    // This is a linear system over GF(2) with m equations and l unknowns.
    // Since the basis is linearly independent (it spans GF(2^l)), the system has a unique
    // solution.  We use Gaussian elimination on the m×l augmented matrix.
    //
    // We work with the low 64 bits of each basis vector (valid for m ≤ 63).
    let m = poly_m.bit_len() - 1;
    let basis = subfield_basis(l, poly_m);
    let a_bits = a.to_uint().as_words()[0];
    let basis_bits: Vec<u64> = basis.iter().map(|b| b.to_uint().as_words()[0]).collect();

    // Build the augmented matrix [A | b] over GF(2), where:
    //   A[j][i] = bit j of basis[i]   (coefficient of cᵢ in equation j)
    //   b[j]    = bit j of a           (RHS of equation j)
    // Row j: bits 0..l-1 are the A columns, bit l is the RHS.
    // We have m rows (one per bit position of GF(2^m) elements).
    let mut aug: Vec<u64> = (0..m)
        .map(|j| {
            let row: u64 = basis_bits
                .iter()
                .enumerate()
                .fold(0u64, |acc, (i, &b)| acc | (((b >> j) & 1) << i));
            let rhs = (a_bits >> j) & 1;
            row | (rhs << l)
        })
        .collect();

    // Gaussian elimination on the m×(l+1) augmented matrix.
    // We pivot on the l unknown columns (0..l).
    let mut pivot_row = 0usize;
    let mut pivot_cols = vec![usize::MAX; l];
    for col in 0..l {
        let found = aug[pivot_row..].iter().position(|r| (r >> col) & 1 == 1);
        if let Some(p) = found {
            aug.swap(pivot_row, pivot_row + p);
            pivot_cols[pivot_row] = col;
            let pivot_val = aug[pivot_row];
            for (i, row) in aug.iter_mut().enumerate().take(m) {
                if i != pivot_row && (*row >> col) & 1 == 1 {
                    *row ^= pivot_val;
                }
            }
            pivot_row += 1;
        }
    }

    // Extract solution: bit pivot_cols[r] of the result = RHS bit of aug[r].
    let mut result: u64 = 0;
    for r in 0..pivot_row {
        let col = pivot_cols[r];
        let rhs = (aug[r] >> l) & 1;
        result |= rhs << col;
    }

    Some(F2mNaive::<L>::from_u64(result, poly_m))
}

// ── Frobenius-by-subfield orbit ───────────────────────────────────────────────

/// Compute the Frobenius-by-subfield orbit of `a`: `[a, a^(2^l), a^(2^(2l)), …]`.
///
/// The orbit has length `m/l`.  Each step applies the l-th Frobenius power (i.e. applies
/// `frobenius` `l` times from the previous entry).  The orbit is the set of conjugates of
/// `a` under the relative Frobenius `φ_l : x ↦ x^(2^l)`.
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly.bit_len() − 1`.
#[must_use]
pub fn frobenius_subfield_orbit<const L: usize>(
    a: &F2mNaive<L>,
    l: usize,
    poly: &Uint<L>,
) -> Vec<F2mNaive<L>>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    assert_l_divides_m(l, m);
    let steps = m / l; // orbit length = m/l
    let mut orbit = Vec::with_capacity(steps);
    let mut cur = a.clone();
    for _ in 0..steps {
        orbit.push(cur.clone());
        // Advance by applying Frobenius l times: cur → cur^(2^l).
        for _ in 0..l {
            cur = cur.frobenius(poly);
        }
    }
    orbit
}

// ── Relative trace ────────────────────────────────────────────────────────────

/// Compute the relative trace `Tr_{m/l}(a) = Σ_{i=0}^{m/l−1} a^(2^(il))`.
///
/// Sums `m/l` terms, each obtained by applying the l-th Frobenius power from the previous.
/// The result lands in GF(2^l) (the Frobenius fixed field).
///
/// **Char-2 note:** addition is XOR.  The relative trace is NOT the absolute trace
/// (which sums `m` terms stepping by 1 and lands in GF(2)).
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly.bit_len() − 1`.
#[must_use]
pub fn relative_trace<const L: usize>(
    a: &F2mNaive<L>,
    l: usize,
    poly: &Uint<L>,
) -> F2mNaive<L>
where
    F2mNaive<L>: F2m<L>,
{
    let orbit = frobenius_subfield_orbit(a, l, poly);
    // Sum all orbit elements (XOR in char 2).
    let mut acc = F2mNaive::<L>::zero();
    for term in &orbit {
        acc = acc.add(term);
    }
    acc
}

// ── Relative norm ─────────────────────────────────────────────────────────────

/// Compute the relative norm `N_{m/l}(a) = Π_{i=0}^{m/l−1} a^(2^(il))`.
///
/// Products `m/l` terms, each obtained by applying the l-th Frobenius power from the
/// previous.  The result lands in GF(2^l).  The norm is multiplicative:
/// `N_{m/l}(a·b) = N_{m/l}(a)·N_{m/l}(b)`.
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly.bit_len() − 1`.
#[must_use]
pub fn relative_norm<const L: usize>(
    a: &F2mNaive<L>,
    l: usize,
    poly: &Uint<L>,
) -> F2mNaive<L>
where
    F2mNaive<L>: F2m<L>,
{
    let orbit = frobenius_subfield_orbit(a, l, poly);
    // Product of all orbit elements.
    let mut acc = F2mNaive::<L>::one();
    for term in &orbit {
        acc = acc.mul(term, poly);
    }
    acc
}

// ── Subfield basis ────────────────────────────────────────────────────────────

/// Compute the polynomial basis `{1, β, β², …, β^(l−1)}` for GF(2^l) inside GF(2^m).
///
/// `β = α^((2^m − 1)/(2^l − 1))` is a primitive element of GF(2^l) expressed in GF(2^m),
/// where `α` is the root of `poly_m` (the primitive element of GF(2^m)).  The basis has
/// `l` elements.
///
/// This basis is over-specified for the GHS/Weil descent and carried for index calculus over
/// GF(2^l) (see `rho::index_calculus`).
///
/// # Panics
///
/// Panics if `l` does not divide `m = poly_m.bit_len() − 1`.
/// Panics if `m > 63` (the exponent `(2^m − 1)/(2^l − 1)` would overflow `u64`).
#[must_use]
pub fn subfield_basis<const L: usize>(
    l: usize,
    poly_m: &Uint<L>,
) -> Vec<F2mNaive<L>>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly_m.bit_len() - 1;
    assert_l_divides_m(l, m);
    assert!(
        m <= 63,
        "subfield_basis: m = {m} > 63; exponent (2^m − 1)/(2^l − 1) overflows u64"
    );

    // α is the root of poly_m, represented as the polynomial x (bit-vector 0b10 = 2).
    let alpha = F2mNaive::<L>::from_u64(2, poly_m);

    // β = α^e where e = (2^m − 1) / (2^l − 1).
    // Both 2^m − 1 and 2^l − 1 fit in u64 for m ≤ 63.
    let field_order: u64 = (1u64 << m) - 1; // 2^m − 1
    let subfield_order: u64 = (1u64 << l) - 1; // 2^l − 1
    let e = field_order / subfield_order; // exact since l | m
    let e_uint = Uint::<L>::from(e);
    let beta = alpha.pow(&e_uint, poly_m);

    // Build {1, β, β², …, β^(l−1)}.
    let mut basis = Vec::with_capacity(l);
    let mut cur = F2mNaive::<L>::one();
    for _ in 0..l {
        basis.push(cur.clone());
        cur = cur.mul(&beta, poly_m);
    }
    basis
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // GF(2^6) with irreducible x⁶ + x + 1 = 0x43.
    fn poly6() -> Uint<1> {
        Uint::<1>::from(0x43u64)
    }

    // GF(2^2) with irreducible x² + x + 1 = 0x7.
    fn poly2() -> Uint<1> {
        Uint::<1>::from(0x7u64)
    }

    type F = F2mNaive<1>;

    fn f6(v: u64) -> F {
        F::from_u64(v, &poly6())
    }

    #[test]
    fn is_in_subfield_zero_and_one() {
        let p6 = poly6();
        // 0 and 1 are in every subfield.
        assert!(is_in_subfield(&F::zero(), 2, &p6));
        assert!(is_in_subfield(&F::one(), 2, &p6));
        assert!(is_in_subfield(&F::zero(), 3, &p6));
        assert!(is_in_subfield(&F::one(), 3, &p6));
    }

    #[test]
    fn gf2_subfield_has_4_elements() {
        // GF(2^2) inside GF(2^6): exactly 4 elements satisfy a^4 = a.
        let p6 = poly6();
        let count = (0u64..64).filter(|&v| is_in_subfield(&f6(v), 2, &p6)).count();
        assert_eq!(count, 4, "GF(2^2) should have exactly 4 elements inside GF(2^6)");
    }

    #[test]
    fn gf3_subfield_has_8_elements() {
        // GF(2^3) inside GF(2^6): exactly 8 elements satisfy a^8 = a.
        let p6 = poly6();
        let count = (0u64..64).filter(|&v| is_in_subfield(&f6(v), 3, &p6)).count();
        assert_eq!(count, 8, "GF(2^3) should have exactly 8 elements inside GF(2^6)");
    }

    #[test]
    fn relative_trace_lands_in_subfield_l2() {
        // Tr_{6/2}(a) ∈ GF(2^2) for all a ∈ GF(2^6).
        let p6 = poly6();
        for v in 0u64..64 {
            let a = f6(v);
            let tr = relative_trace(&a, 2, &p6);
            assert!(
                is_in_subfield(&tr, 2, &p6),
                "Tr_{{6/2}}({v:#x}) = {tr:?} not in GF(2^2)"
            );
        }
    }

    #[test]
    fn relative_trace_lands_in_subfield_l3() {
        // Tr_{6/3}(a) ∈ GF(2^3) for all a ∈ GF(2^6).
        let p6 = poly6();
        for v in 0u64..64 {
            let a = f6(v);
            let tr = relative_trace(&a, 3, &p6);
            assert!(
                is_in_subfield(&tr, 3, &p6),
                "Tr_{{6/3}}({v:#x}) = {tr:?} not in GF(2^3)"
            );
        }
    }

    #[test]
    fn relative_norm_lands_in_subfield_l2() {
        // N_{6/2}(a) ∈ GF(2^2) for all a ∈ GF(2^6).
        let p6 = poly6();
        for v in 0u64..64 {
            let a = f6(v);
            let n = relative_norm(&a, 2, &p6);
            assert!(
                is_in_subfield(&n, 2, &p6),
                "N_{{6/2}}({v:#x}) = {n:?} not in GF(2^2)"
            );
        }
    }

    #[test]
    fn frobenius_orbit_length() {
        let p6 = poly6();
        let a = f6(0x15);
        let orbit = frobenius_subfield_orbit(&a, 2, &p6);
        assert_eq!(orbit.len(), 3, "orbit length should be m/l = 6/2 = 3");
    }

    #[test]
    fn embed_restrict_round_trip() {
        // For elements in GF(2^2), restrict(embed(c)) = Some(c).
        let p6 = poly6();
        let p2 = poly2();
        for v in 0u64..4 {
            let c = F::from_u64(v, &p2);
            let embedded = embed(&c, &p2, &p6);
            let restricted = restrict(&embedded, 2, &p6);
            assert!(restricted.is_some(), "restrict(embed({v:#x})) should be Some");
            assert_eq!(
                restricted.unwrap().to_uint(),
                c.to_uint(),
                "embed∘restrict round-trip failed for {v:#x}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "l must divide m")]
    fn l_not_dividing_m_panics() {
        let p6 = poly6();
        let a = f6(1);
        // l=4 does not divide m=6 → should panic.
        let _ = is_in_subfield(&a, 4, &p6);
    }
}
