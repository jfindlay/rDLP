//! Polynomial-basis ↔ normal-basis change-of-basis isomorphism for GF(2^m).
//!
//! A **normal basis** for GF(2^m) over GF(2) is the set
//! `{β, β², β⁴, …, β^(2^(m−1))}` where β is a *normal element* — an element
//! whose Frobenius orbit spans the whole field.  Every element `a ∈ GF(2^m)`
//! can be written uniquely as
//!
//! ```text
//! a = a₀·β + a₁·β² + a₂·β⁴ + … + a_{m−1}·β^(2^(m−1))
//! ```
//!
//! where each `aᵢ ∈ GF(2)`.  The normal-basis representation stores the
//! coefficient vector `(a₀, a₁, …, a_{m−1})` as a bit-vector.
//!
//! # Algorithms
//!
//! - **`frobenius_orbit`**: compute the m basis vectors `β^(2^i)` in polynomial
//!   basis by repeated squaring.
//! - **`is_normal_element`**: check that the m×m matrix whose rows are the
//!   Frobenius orbit has rank m over GF(2) (Gaussian elimination).
//! - **`find_normal_element`**: search GF(2^m) for the smallest normal element.
//! - **`normal_to_poly`**: given normal-basis coordinates, XOR the corresponding
//!   basis vectors to recover the polynomial-basis representation.
//! - **`poly_to_normal`**: solve the linear system M·x = a over GF(2) (where M
//!   has rows β^(2^i)) via Gaussian elimination to find the normal-basis
//!   coordinates.
//!
//! # Storage contract
//!
//! All `Uint<L>` values are polynomial coefficient bit-vectors (bit `i` = coeff
//! of `x^i`).  Only XOR, shift, and bit operations are meaningful.

use crate::BitLen;
use crypto_bigint::Uint;

use crate::naive::F2mNaive;
use crate::F2m;

// ── Frobenius orbit ───────────────────────────────────────────────────────────

/// Compute the Frobenius orbit of `beta` in GF(2^m): `[β, β², β⁴, …, β^(2^(m−1))]`.
///
/// Each entry is in polynomial-basis representation.  The orbit has exactly `m`
/// elements when `beta` is a normal element.
///
/// `poly` is the irreducible polynomial defining GF(2^m); `m = poly.bit_len() − 1`.
#[must_use]
pub fn frobenius_orbit<const L: usize>(beta: &Uint<L>, poly: &Uint<L>) -> Vec<Uint<L>>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1; // degree of the field
    let mut orbit = Vec::with_capacity(m);
    let mut cur = F2mNaive::<L>::from_uint(*beta, poly);
    for _ in 0..m {
        orbit.push(cur.to_uint());
        cur = cur.square(poly);
    }
    orbit
}

// ── Normal-element check ──────────────────────────────────────────────────────

/// Return `true` if `beta` is a normal element of GF(2^m) over GF(2).
///
/// `beta` is normal iff its Frobenius orbit `{β, β², β⁴, …, β^(2^(m−1))}`
/// consists of exactly `m` linearly independent vectors over GF(2).
///
/// Equivalently, the m×m matrix M whose rows are the orbit vectors (as
/// polynomial-basis bit-vectors, taking only the low `m` bits) has rank `m`
/// over GF(2).
///
/// `poly` is the irreducible polynomial; `m = poly.bit_len() − 1`.
#[must_use]
pub fn is_normal_element<const L: usize>(beta: &Uint<L>, poly: &Uint<L>) -> bool
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    let orbit = frobenius_orbit(beta, poly);
    if orbit.len() != m {
        return false;
    }
    gf2_matrix_rank_m(&orbit, m) == m
}

/// Compute the rank of an m×m GF(2) matrix given as a slice of `m` row vectors.
///
/// Each row is a `Uint<L>` whose low `m` bits are the row entries.
/// Uses Gaussian elimination over GF(2) (XOR-based).
fn gf2_matrix_rank_m<const L: usize>(rows: &[Uint<L>], m: usize) -> usize {
    // Work with u64 for small m (m ≤ 64 covers all toy fields).
    // Extract the low 64 bits of each row.
    let mut mat: Vec<u64> = rows
        .iter()
        .map(|r| r.as_words()[0]) // low 64-bit word
        .collect();

    let mut rank = 0usize;
    for col in (0..m).rev() {
        // Find a pivot row at or below `rank` with bit `col` set.
        let pivot = mat[rank..].iter().position(|r| (r >> col) & 1 == 1);
        if let Some(p) = pivot {
            mat.swap(rank, rank + p);
            let pivot_row = mat[rank];
            // Eliminate this column from all other rows.
            for (i, row) in mat.iter_mut().enumerate() {
                if i != rank && (*row >> col) & 1 == 1 {
                    *row ^= pivot_row;
                }
            }
            rank += 1;
        }
    }
    rank
}

// ── Find a normal element ─────────────────────────────────────────────────────

/// Find the smallest non-zero normal element of GF(2^m) over GF(2).
///
/// Searches elements `1, 2, 3, …` in order and returns the first one whose
/// Frobenius orbit spans the field.
///
/// # Panics
///
/// Panics if no normal element is found (which cannot happen for a valid
/// irreducible polynomial, since normal bases always exist over finite fields).
#[must_use]
pub fn find_normal_element<const L: usize>(poly: &Uint<L>) -> Uint<L>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    let field_size: u64 = 1u64 << m; // 2^m
    for v in 1u64..field_size {
        let beta = Uint::<L>::from(v);
        if is_normal_element(&beta, poly) {
            return beta;
        }
    }
    panic!("find_normal_element: no normal element found — irreducible polynomial may be invalid");
}

// ── Basis conversion ──────────────────────────────────────────────────────────

/// Convert a normal-basis representation to polynomial basis.
///
/// Given normal-basis coordinates `a_normal` (bit `i` = coefficient `aᵢ` of
/// `β^(2^i)`) and the normal element `beta` (in polynomial basis), compute
///
/// ```text
/// a_poly = Σ aᵢ · β^(2^i)   (sum in GF(2^m), i.e. XOR)
/// ```
///
/// The result is the polynomial-basis representation of the same field element.
///
/// `poly` is the irreducible polynomial; `beta` must be a normal element.
#[must_use]
pub fn normal_to_poly<const L: usize>(
    a_normal: &Uint<L>,
    beta: &Uint<L>,
    poly: &Uint<L>,
) -> Uint<L>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    let orbit = frobenius_orbit(beta, poly);
    debug_assert_eq!(orbit.len(), m, "normal_to_poly: beta is not a normal element");

    // XOR the orbit vectors corresponding to set bits in a_normal.
    let mut result = Uint::<L>::ZERO;
    for (i, &basis_vec) in orbit.iter().enumerate().take(m) {
        if a_normal.bit_at(i) {
            result ^= basis_vec;
        }
    }
    result
}

/// Convert a polynomial-basis representation to normal basis.
///
/// Given `a_poly` (polynomial-basis representation) and the normal element
/// `beta`, find the normal-basis coordinates `a_normal` such that
///
/// ```text
/// a_poly = Σ aᵢ · β^(2^i)
/// ```
///
/// This is done by solving the linear system over GF(2).  The system has one
/// equation per bit position `j` (0 ≤ j < m):
///
/// ```text
/// Σᵢ xᵢ · orbit[i][j] = a_poly[j]
/// ```
///
/// where `orbit[i][j]` is bit `j` of `β^(2^i)` in polynomial basis.
///
/// The augmented matrix has `m` rows (one per bit position `j`) and `m+1`
/// columns: column `i` is the coefficient of `xᵢ` in equation `j`, and the
/// last column is the RHS `a_poly[j]`.
///
/// Uses Gaussian elimination over GF(2).
///
/// # Panics
///
/// Panics if `beta` is not a normal element (the system has no unique solution).
#[must_use]
pub fn poly_to_normal<const L: usize>(
    a_poly: &Uint<L>,
    beta: &Uint<L>,
    poly: &Uint<L>,
) -> Uint<L>
where
    F2mNaive<L>: F2m<L>,
{
    let m = poly.bit_len() - 1;
    let orbit = frobenius_orbit(beta, poly);
    debug_assert_eq!(orbit.len(), m, "poly_to_normal: beta is not a normal element");

    // Build the augmented matrix [A | b] over GF(2), where:
    //   A[j][i] = bit j of orbit[i]   (coefficient of xᵢ in equation j)
    //   b[j]    = bit j of a_poly      (RHS of equation j)
    //
    // We work with u64 for small m (m ≤ 64).
    // Each augmented row is a u64: bits [m-1..0] are the A columns (xᵢ coefficients),
    // bit m is the RHS.
    //
    // Row j: A[j][i] = (orbit[i].as_words()[0] >> j) & 1 for i in 0..m
    //        b[j]    = (a_poly.as_words()[0] >> j) & 1
    let a_bits = a_poly.as_words()[0];
    let orbit_bits: Vec<u64> = orbit.iter().map(|r| r.as_words()[0]).collect();

    let mut aug: Vec<u64> = (0..m)
        .map(|j| {
            // Build row j: bit i of the row = A[j][i] = bit j of orbit[i].
            let row: u64 = orbit_bits
                .iter()
                .enumerate()
                .take(m)
                .fold(0u64, |acc, (i, &ob)| acc | (((ob >> j) & 1) << i));
            // Append RHS as bit m.
            let rhs = (a_bits >> j) & 1;
            row | (rhs << m)
        })
        .collect();

    // Gaussian elimination over GF(2) on the augmented matrix.
    // Pivot on columns 0..m (the xᵢ columns).
    let mut pivot_row = 0usize;
    let mut pivot_cols = vec![usize::MAX; m]; // pivot_cols[r] = column of pivot in row r

    for col in 0..m {
        // Find a row at or below pivot_row with bit `col` set.
        let found = aug[pivot_row..].iter().position(|r| (r >> col) & 1 == 1);
        if let Some(p) = found {
            aug.swap(pivot_row, pivot_row + p);
            pivot_cols[pivot_row] = col;
            let pivot_val = aug[pivot_row];
            // Eliminate this column from all other rows.
            for (i, row) in aug.iter_mut().enumerate().take(m) {
                if i != pivot_row && (*row >> col) & 1 == 1 {
                    *row ^= pivot_val;
                }
            }
            pivot_row += 1;
        }
    }

    debug_assert_eq!(
        pivot_row, m,
        "poly_to_normal: beta is not a normal element — matrix is singular"
    );

    // Back-substitute: each row r has a unique pivot column pivot_cols[r].
    // The solution bit for that column is the RHS bit (bit m of aug[r]).
    let mut result: u64 = 0;
    for r in 0..m {
        let col = pivot_cols[r];
        let rhs = (aug[r] >> m) & 1;
        result |= rhs << col;
    }

    Uint::<L>::from(result)
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64) // x⁴ + x + 1
    }

    fn poly8() -> Uint<1> {
        Uint::<1>::from(0x11bu64) // x⁸ + x⁴ + x³ + x + 1 (AES)
    }

    #[test]
    fn find_normal_gf4() {
        let p = poly4();
        let beta = find_normal_element(&p);
        assert!(
            is_normal_element(&beta, &p),
            "find_normal_element returned non-normal element 0x{:x}",
            beta
        );
    }

    #[test]
    fn find_normal_gf8() {
        let p = poly8();
        let beta = find_normal_element(&p);
        assert!(
            is_normal_element(&beta, &p),
            "find_normal_element returned non-normal element 0x{:x}",
            beta
        );
    }

    #[test]
    fn round_trip_poly_normal_poly_gf4() {
        let p = poly4();
        let beta = find_normal_element(&p);
        for v in 0u64..16 {
            let a_poly = Uint::<1>::from(v);
            let a_norm = poly_to_normal(&a_poly, &beta, &p);
            let a_back = normal_to_poly(&a_norm, &beta, &p);
            assert_eq!(
                a_poly, a_back,
                "round-trip failed for v={v:#x}: poly→normal→poly"
            );
        }
    }

    #[test]
    fn round_trip_normal_poly_normal_gf4() {
        let p = poly4();
        let beta = find_normal_element(&p);
        for v in 0u64..16 {
            let a_norm = Uint::<1>::from(v);
            let a_poly = normal_to_poly(&a_norm, &beta, &p);
            let a_back = poly_to_normal(&a_poly, &beta, &p);
            assert_eq!(
                a_norm, a_back,
                "round-trip failed for v={v:#x}: normal→poly→normal"
            );
        }
    }
}
