//! Pohlig–Hellman ECDLP: order factorization, prime-subgroup projection, and composite-order
//! ECDLP reduction.
//!
//! This module provides the full Pohlig–Hellman reduction:
//!
//! - [`factor_order`] — prime-power decomposition of a u64 group order.
//! - [`project_to_subgroup`] — map `(G, Q)` to the order-`p^e` subgroup via
//!   `[n / p^e]`-scalar-multiplication.
//! - [`solve_ecdlp_composite`] — composite-order ECDLP via prime-power lift + CRT combine.
//!
//! # Algorithm
//!
//! [`solve_ecdlp_composite`] decomposes the DLP `Q = k·G` in a group of composite order
//! `n = ∏ pᵢ^{eᵢ}` into independent DLPs in each prime-power subgroup, solves each by
//! handing off to the frozen rho solvers (which require a **prime** order), and reassembles
//! the answer by the Chinese Remainder Theorem.
//!
//! The prime-power lift (for `eᵢ > 1`) recovers `k mod pᵢ^{eᵢ}` digit-by-digit in base `pᵢ`:
//! at each digit step, the current remainder is projected to the order-`pᵢ` sub-subgroup and
//! solved by a single rho call with prime order `pᵢ`.
//!
//! # Scope
//!
//! Toy-order-scoped (u64 group orders). `factor_order` uses trial division via
//! `shared_numth::trial_smooth` over a `√n` factor base, with an `is_prime`
//! short-circuit on the cofactor. For orders beyond u64 (not needed here),
//! `shared_numth::ecm_factor` is available as a fallback — noted as a
//! principle-4 annotation, not wired.

use crypto_bigint::Uint;

use shared_numth::{factor_base_up_to, is_prime, trial_smooth};

use shared_field::Fp;

use crate::curve::{AffinePoint, Curve, JacobianPoint};
use crate::ecdlp::solve_brent;

// ── factor_order ──────────────────────────────────────────────────────────────

/// Compute the prime-power decomposition of a u64 group order.
///
/// Returns a sorted `Vec<(prime, exponent)>` such that `∏ pᵢ^{eᵢ} = n`.
/// Uses trial division over a `√n` factor base (`factor_base_up_to(isqrt(n))`)
/// with an `is_prime` short-circuit on the remaining cofactor.
///
/// Toy-order-scoped (u64). For orders that require ECM, `shared_numth::ecm_factor`
/// is available as a fallback (principle-4 annotation — not wired here).
///
/// # Panics
///
/// Panics if `n == 0`.
pub fn factor_order(n: u64) -> Vec<(u64, u32)> {
    assert!(n != 0, "factor_order: n must be non-zero");

    if n == 1 {
        return vec![];
    }

    // Build a factor base up to √n for trial division.
    let sqrt_n = isqrt(n);
    let base = factor_base_up_to(sqrt_n);

    let n_uint = Uint::<4>::from(n);
    let witness = trial_smooth(&n_uint, &base);

    let mut factors = witness.factors.clone();

    // The cofactor is whatever trial_smooth could not factor.
    // If cofactor > 1, it must be prime (since the base covers all primes ≤ √n,
    // and any composite cofactor would have a factor ≤ √cofactor ≤ √n).
    let cofactor_uint = witness.cofactor;
    if cofactor_uint > Uint::<4>::ONE {
        // Sanity: the cofactor should be prime.
        debug_assert!(
            is_prime(&cofactor_uint),
            "factor_order: cofactor is composite — base too small?"
        );
        let cofactor_u64 = cofactor_uint.as_words()[0];
        factors.push((cofactor_u64, 1));
        factors.sort_unstable_by_key(|&(p, _)| p);
    }

    factors
}

// ── isqrt ─────────────────────────────────────────────────────────────────────

/// Integer square root: largest `k` such that `k² ≤ n`.
fn isqrt(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut x = (n as f64).sqrt() as u64;
    // Correct for floating-point rounding.
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

// ── project_to_subgroup ───────────────────────────────────────────────────────

/// Project `(G, Q)` to the order-`p^e` subgroup of a curve of order `n`.
///
/// Given a curve of order `n` with `p^e | n`, maps both the base point and the
/// target to the unique subgroup of order `p^e` by multiplying by the cofactor
/// `n / p^e`. Returns `(G', Q')` where `G' = [n/p^e]·G` and `Q' = [n/p^e]·Q`.
///
/// `G'` generates the order-`p^e` subgroup (provided `G` has full order `n`).
/// The DLP `Q' = k·G'` in the subgroup gives `k mod p^e`.
///
/// # Preconditions
///
/// - `p^e` must divide `n` exactly.
/// - `G` must have full order `n` (not a subgroup generator).
pub fn project_to_subgroup<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    p_power: u64,
) -> (AffinePoint<F>, AffinePoint<F>) {
    debug_assert!(n % p_power == 0, "project_to_subgroup: p^e does not divide n");
    let cofactor = n / p_power;
    let cofactor_uint = Uint::<4>::from(cofactor);
    let g_sub = curve.scalar_mul(g, &cofactor_uint);
    let q_sub = curve.scalar_mul(q, &cofactor_uint);
    (g_sub, q_sub)
}

// ── solve_ecdlp_composite ─────────────────────────────────────────────────────

/// Solve `Q = k·G` on a curve whose generator has composite order `n`.
///
/// Decomposes the DLP into independent DLPs in each prime-power subgroup via
/// Pohlig–Hellman, solves each with the frozen rho solver (`solve_brent`), and
/// reassembles the answer by the Chinese Remainder Theorem.
///
/// # Algorithm
///
/// For each prime power `pᵢ^{eᵢ}` in the factorization of `n`:
/// 1. Project `(G, Q)` to the order-`pᵢ^{eᵢ}` subgroup.
/// 2. Recover `k mod pᵢ^{eᵢ}` digit-by-digit in base `pᵢ` (the prime-power lift).
///    Each digit is solved by a rho call with prime order `pᵢ` — never `pᵢ^{eᵢ}`,
///    because `inv_mod_prime` inside the rho solver uses Fermat's little theorem and
///    is silently wrong on composite moduli.
/// 3. CRT-combine the per-prime-power residues into `k mod n`.
///
/// # Returns
///
/// `Some(k)` with `k·G = Q`, or `None` if any per-subgroup rho call fails.
///
/// # Preconditions
///
/// - `G` must have full order `n` (not a subgroup generator).
/// - `n` must equal the product of the prime powers returned by `factor_order(n)`.
pub fn solve_ecdlp_composite<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
) -> Option<u64> {
    // Handle the trivial case: Q = ∞ means k = 0 (or n, but 0 is canonical).
    if q.is_infinity() {
        return Some(0);
    }

    let factors = factor_order(n);

    // Collect (x_i mod p_i^e_i, p_i^e_i) for each prime power.
    let mut residues: Vec<(u64, u64)> = Vec::with_capacity(factors.len());

    for &(p, e) in &factors {
        let p_power = p.pow(e); // pᵢ^{eᵢ}
        let x_pe = solve_prime_power(curve, g, q, n, p, e, p_power)?;
        residues.push((x_pe, p_power));
    }

    // CRT-combine all residues into k mod n.
    Some(crt_combine(&residues, n))
}

/// Solve `d · G = Q` by brute force for small prime order `p` (p ≤ 64).
///
/// Iterates `d = 1, 2, …, p-1`, computing `d·G` by repeated addition.
/// Returns `Some(d)` when `d·G = Q`, or `None` if no solution exists in `[0, p)`.
///
/// Used instead of `solve_brent` when `p` is small enough that the rho walk
/// table (20 entries) would exceed the group size, causing excessive degeneration.
fn solve_small_dlog<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    p: u64,
) -> Option<u64> {
    // d=0: G·0 = ∞; caller already handles the ∞ case before calling here.
    let mut acc = g.clone(); // 1·G
    for d in 1..p {
        if &acc == q {
            return Some(d);
        }
        // acc ← acc + G
        let acc_jac = JacobianPoint::from_affine(&acc, &curve.p);
        acc = curve.add_mixed(&acc_jac, g).to_affine(&curve.p);
    }
    // d = p: p·G = ∞ (group order), which is the identity — not a valid solution
    // for a non-identity Q (caller checks Q ≠ ∞ before calling).
    None
}

/// Recover `k mod p^e` by the digit-by-digit Pohlig–Hellman prime-power lift.
///
/// Projects `(G, Q)` to the order-`p^e` subgroup, then iteratively recovers each
/// base-`p` digit of `k mod p^e` by projecting the current remainder to the
/// order-`p` sub-subgroup and solving a prime-order DLP with `solve_brent`.
///
/// Every rho call receives the prime order `p`, never `p^e` — the prime-only
/// precondition of `inv_mod_prime` (Fermat) is honoured throughout.
fn solve_prime_power<F: Fp<4>>(
    curve: &Curve,
    g: &AffinePoint<F>,
    q: &AffinePoint<F>,
    n: u64,
    p: u64,
    e: u32,
    p_power: u64, // p^e
) -> Option<u64> {
    // Project (G, Q) to the order-p^e subgroup.
    let (g_pe, q_pe) = project_to_subgroup(curve, g, q, n, p_power);

    // γ = p^(e-1)·G_pe has order p — the generator of the order-p sub-subgroup.
    let p_e_minus_1 = p.pow(e - 1); // p^(e-1)
    let gamma = curve.scalar_mul(&g_pe, &Uint::<4>::from(p_e_minus_1));

    // Accumulate x mod p^e digit by digit.
    let mut x_pe: u64 = 0;
    // Q_curr tracks the "remaining" target after subtracting known digits.
    let mut q_curr = q_pe.clone();

    for j in 0..e {
        // Project Q_curr to the order-p sub-subgroup: rhs = p^(e-1-j) · Q_curr.
        // For j = e-1, the exponent is 0 and rhs = Q_curr (already order-p).
        let exp = p.pow(e - 1 - j); // p^(e-1-j)
        let rhs = curve.scalar_mul(&q_curr, &Uint::<4>::from(exp));

        // Solve d_j · γ = rhs with prime order p.
        // γ has order p (prime), so solve_brent is valid here.
        // For small p (≤ 64), brute-force is faster and avoids rho walk degeneration
        // on tiny groups where the walk table (20 entries) exceeds the group size.
        let d_j = if rhs.is_infinity() {
            // rhs = ∞ means d_j = 0.
            0u64
        } else if p <= 64 {
            solve_small_dlog(curve, &gamma, &rhs, p)?
        } else {
            // Interim adaptation: this function still reports a single failure class.
            solve_brent(curve, &gamma, &rhs, p, 0, 50).ok().map(|s| s.k)?
        };

        // Accumulate: x_pe += d_j · p^j.
        x_pe += d_j * p.pow(j);

        // Update Q_curr: subtract d_j · p^j · G_pe.
        // Q_curr ← Q_curr − d_j · p^j · G_pe
        if d_j != 0 {
            let scalar = d_j * p.pow(j); // d_j · p^j, fits in u64 since < p^e ≤ n
            let correction = curve.scalar_mul(&g_pe, &Uint::<4>::from(scalar));
            let neg_correction = curve.negate(&correction);
            let q_curr_jac = JacobianPoint::from_affine(&q_curr, &curve.p);
            q_curr = curve.add_mixed(&q_curr_jac, &neg_correction).to_affine(&curve.p);
        }
        // If d_j == 0, Q_curr is unchanged (subtracting 0·G_pe = ∞ is a no-op).
    }

    Some(x_pe)
}

/// Combine residues `{(xᵢ, mᵢ)}` into `x mod M` by the Chinese Remainder Theorem.
///
/// `mᵢ` are pairwise coprime (they are distinct prime powers), and `M = ∏ mᵢ`.
/// Uses the standard constructive CRT formula:
/// `x = Σ xᵢ · Mᵢ · (Mᵢ⁻¹ mod mᵢ) mod M`, where `Mᵢ = M / mᵢ`.
fn crt_combine(residues: &[(u64, u64)], modulus: u64) -> u64 {
    let mut x: u64 = 0;
    for &(xi, mi) in residues {
        let big_mi = modulus / mi; // Mᵢ = M / mᵢ
        // Mᵢ⁻¹ mod mᵢ: since gcd(Mᵢ, mᵢ) = 1, the inverse exists.
        let mi_inv = mod_inv(big_mi % mi, mi);
        // Contribution: xᵢ · Mᵢ · (Mᵢ⁻¹ mod mᵢ) mod M.
        // Use u128 arithmetic to avoid overflow (all values ≤ n ≤ u64::MAX).
        let term = ((xi as u128 * (big_mi as u128 % modulus as u128) % modulus as u128)
            * mi_inv as u128)
            % modulus as u128;
        x = ((x as u128 + term) % modulus as u128) as u64;
    }
    x
}

/// Modular inverse of `a` mod `m` via the extended Euclidean algorithm.
///
/// Requires `gcd(a, m) = 1`. Returns `r` such that `a · r ≡ 1 (mod m)`.
///
/// # Panics
///
/// Panics if `a == 0` or `m == 0`.
fn mod_inv(a: u64, m: u64) -> u64 {
    assert!(a != 0 && m != 0, "mod_inv: zero input");
    // Extended Euclidean algorithm over i128 to handle signed intermediates.
    let mut old_r = a as i128;
    let mut r = m as i128;
    let mut old_s: i128 = 1;
    let mut s: i128 = 0;

    while r != 0 {
        let q = old_r / r;
        let tmp_r = r;
        r = old_r - q * r;
        old_r = tmp_r;
        let tmp_s = s;
        s = old_s - q * s;
        old_s = tmp_s;
    }

    // old_r is gcd; old_s is the Bézout coefficient for a.
    debug_assert_eq!(old_r, 1, "mod_inv: gcd({a}, {m}) ≠ 1");
    ((old_s % m as i128 + m as i128) % m as i128) as u64
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_field::FpMonty4 as FpMonty;

    use crate::curve::test_curves::{composite_toy, COMPOSITE_TOY_FACTORS, COMPOSITE_TOY_N};

    // ── factor_order KATs ─────────────────────────────────────────────────────

    /// factor_order reproduces the recorded factorization of the composite fixture.
    #[test]
    fn factor_order_composite_toy() {
        let factors = factor_order(COMPOSITE_TOY_N);
        assert_eq!(
            factors, COMPOSITE_TOY_FACTORS,
            "factor_order({COMPOSITE_TOY_N}) = {factors:?}, expected {COMPOSITE_TOY_FACTORS:?}"
        );
    }

    /// The product of p^e for all factors equals COMPOSITE_TOY_N.
    #[test]
    fn factor_order_product_equals_n() {
        let factors = factor_order(COMPOSITE_TOY_N);
        let product: u64 = factors.iter().map(|&(p, e)| p.pow(e)).product();
        assert_eq!(
            product, COMPOSITE_TOY_N,
            "∏ pᵢ^eᵢ = {product}, expected {COMPOSITE_TOY_N}"
        );
    }

    /// factor_order handles small known values correctly.
    #[test]
    fn factor_order_small_cases() {
        assert_eq!(factor_order(1), vec![]);
        assert_eq!(factor_order(2), vec![(2, 1)]);
        assert_eq!(factor_order(4), vec![(2, 2)]);
        assert_eq!(factor_order(12), vec![(2, 2), (3, 1)]);
        assert_eq!(factor_order(60), vec![(2, 2), (3, 1), (5, 1)]);
        assert_eq!(factor_order(360), vec![(2, 3), (3, 2), (5, 1)]);
    }

    /// factor_order handles prime inputs correctly.
    #[test]
    fn factor_order_prime() {
        assert_eq!(factor_order(7), vec![(7, 1)]);
        assert_eq!(factor_order(97), vec![(97, 1)]);
    }

    // ── project_to_subgroup KATs ──────────────────────────────────────────────

    /// Projecting to the order-4 subgroup (p=2, e=2) gives a point of order 4.
    #[test]
    fn project_order_4_subgroup() {
        let curve = composite_toy();
        let g: AffinePoint<FpMonty> = curve.generator();
        let n = COMPOSITE_TOY_N;
        // p^e = 4 (prime 2, exponent 2)
        let (g_sub, _) = project_to_subgroup(&curve, &g, &g, n, 4);
        // g_sub should have order 4: 4*g_sub = infinity
        let four_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(4u64));
        assert!(four_g_sub.is_infinity(), "project order-4: 4*G' should be ∞");
        // And 2*g_sub != infinity (order is exactly 4, not 2)
        let two_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(2u64));
        assert!(!two_g_sub.is_infinity(), "project order-4: 2*G' should not be ∞");
    }

    /// Projecting to the order-3 subgroup gives a point of order 3.
    #[test]
    fn project_order_3_subgroup() {
        let curve = composite_toy();
        let g: AffinePoint<FpMonty> = curve.generator();
        let n = COMPOSITE_TOY_N;
        // p^e = 3
        let (g_sub, _) = project_to_subgroup(&curve, &g, &g, n, 3);
        let three_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(3u64));
        assert!(three_g_sub.is_infinity(), "project order-3: 3*G' should be ∞");
        let one_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(1u64));
        assert!(!one_g_sub.is_infinity(), "project order-3: G' should not be ∞");
    }

    /// Projecting to the order-5 subgroup gives a point of order 5.
    #[test]
    fn project_order_5_subgroup() {
        let curve = composite_toy();
        let g: AffinePoint<FpMonty> = curve.generator();
        let n = COMPOSITE_TOY_N;
        // p^e = 5
        let (g_sub, _) = project_to_subgroup(&curve, &g, &g, n, 5);
        let five_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(5u64));
        assert!(five_g_sub.is_infinity(), "project order-5: 5*G' should be ∞");
        let one_g_sub = curve.scalar_mul(&g_sub, &Uint::<4>::from(1u64));
        assert!(!one_g_sub.is_infinity(), "project order-5: G' should not be ∞");
    }

    // ── isqrt unit tests ──────────────────────────────────────────────────────

    /// isqrt returns the correct integer square root for small values.
    #[test]
    fn isqrt_small() {
        assert_eq!(isqrt(0), 0);
        assert_eq!(isqrt(1), 1);
        assert_eq!(isqrt(3), 1);
        assert_eq!(isqrt(4), 2);
        assert_eq!(isqrt(8), 2);
        assert_eq!(isqrt(9), 3);
        assert_eq!(isqrt(60), 7);
        assert_eq!(isqrt(100), 10);
    }
}
