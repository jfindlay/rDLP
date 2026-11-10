//! Tate pairing and reduced-Tate pairing with final exponentiation.
//!
//! The **Tate pairing** is a bilinear, non-degenerate map
//!
//! ```text
//! t_ℓ : E[ℓ] × E(F_{p^k}) / ℓ·E(F_{p^k}) → F_{p^k}* / (F_{p^k}*)^ℓ
//! ```
//!
//! computed as a single Miller evaluation:
//!
//! ```text
//! t_ℓ(P, Q) = f_{ℓ,P}(Q)
//! ```
//!
//! The **reduced Tate pairing** (also called the Tate–Lichtenbaum pairing) maps
//! into the unique subgroup `μ_ℓ ⊂ F_{p^k}*` of ℓ-th roots of unity via the
//! **final exponentiation**:
//!
//! ```text
//! τ_ℓ(P, Q) = f_{ℓ,P}(Q)^{(p^k − 1)/ℓ}
//! ```
//!
//! The exponent `(p^k − 1)/ℓ` is the **final exponentiation exponent**.  It
//! divides evenly because `ℓ | p^k − 1` (the embedding-degree condition).
//!
//! # Tate vs Weil
//!
//! | Property | Weil | Tate (reduced) |
//! |----------|------|----------------|
//! | Miller calls | **two** (`f_{ℓ,P}(Q)` and `f_{ℓ,Q}(P)`) | **one** (`f_{ℓ,P}(Q)`) |
//! | Final exponentiation | none (lands in `μ_ℓ` directly via ratio) | `^{(p^k−1)/ℓ}` (projects into `μ_ℓ`) |
//! | Sign factor | `(−1)^ℓ` | none |
//!
//! The Tate pairing is computationally cheaper (one Miller call vs two) and is
//! the standard choice in pairing-based cryptography.  The Weil pairing is
//! mathematically cleaner (no final exponentiation) but requires two Miller
//! evaluations.
//!
//! # Final exponentiation
//!
//! For the toy fixture (`p = 47`, `k = 2`, `ℓ = 3`):
//! - `p^k − 1 = 47² − 1 = 2208 = 3 · 736`
//! - Final exponentiation exponent: `(p^k − 1)/ℓ = 736`
//! - `f_{3,P}(Q)^{736} ∈ μ_3 ⊂ F_{47^2}*`
//!
//! The divisibility `ℓ | p^k − 1` is the embedding-degree condition and holds
//! by construction of the fixture.

use crypto_bigint::Uint;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use shared_field::Fp;

use crate::curve::Curve;
use crate::pairing::ecext::PairingPoint;
use crate::pairing::fpext::{FpExt, IrreducibleModulus};
use crate::pairing::miller::miller_loop;

// ── Public entry points ───────────────────────────────────────────────────────

/// Compute the raw Tate pairing `t_ℓ(P, Q) = f_{ℓ,P}(Q)`.
///
/// Returns the Miller function `f_{ℓ,P}` evaluated at `Q`.  This is the raw
/// (unreduced) Tate pairing value — an element of `F_{p^k}*` that represents a
/// coset in `F_{p^k}* / (F_{p^k}*)^ℓ`.  To obtain a well-defined element of
/// `μ_ℓ`, apply the final exponentiation via [`reduced_tate`].
///
/// # Parameters
///
/// - `curve` — the short-Weierstrass curve `y² = x³ + ax + b`.
/// - `modulus` — the irreducible polynomial defining `F_{p^k}`.
/// - `p_point` — the first pairing argument `P ∈ E[ℓ]`.
/// - `q_point` — the second pairing argument `Q ∈ E(F_{p^k})`, linearly
///   independent of `P` for a non-degenerate result.
/// - `ell` — the torsion order; must be ≥ 2.
///
/// # Panics
///
/// Panics if `ell < 2`.
pub fn tate_pairing<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    p_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    ell: u64,
) -> FpExt<F> {
    assert!(ell >= 2, "tate_pairing: ell must be >= 2");
    // One Miller call — the key distinction from the Weil pairing (which uses two).
    miller_loop(curve, modulus, p_point, q_point, ell)
}

/// Compute the reduced Tate pairing `τ_ℓ(P, Q) = f_{ℓ,P}(Q)^{(p^k − 1)/ℓ}`.
///
/// Applies the **final exponentiation** `^{(p^k − 1)/ℓ}` to the raw Tate
/// pairing value, projecting it into the unique subgroup `μ_ℓ ⊂ F_{p^k}*` of
/// ℓ-th roots of unity.  The result satisfies `result^ℓ = 1`.
///
/// # Final exponentiation
///
/// The exponent `(p^k − 1)/ℓ` must divide evenly — this is guaranteed by the
/// embedding-degree condition `ℓ | p^k − 1`.  For the toy fixture:
/// `(47² − 1)/3 = 2208/3 = 736`.
///
/// # Parameters
///
/// - `curve` — the short-Weierstrass curve `y² = x³ + ax + b`.
/// - `modulus` — the irreducible polynomial defining `F_{p^k}`.
/// - `p_point` — the first pairing argument `P ∈ E[ℓ]`.
/// - `q_point` — the second pairing argument `Q ∈ E(F_{p^k})`, linearly
///   independent of `P` for a non-degenerate result.
/// - `ell` — the torsion order; must be ≥ 2.
///
/// # Panics
///
/// Panics if `ell < 2`, if `ell` does not divide `p^k − 1` (the
/// embedding-degree condition is violated), or if `(p^k − 1)/ℓ` does not fit 256 bits.
pub fn reduced_tate<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    p_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    ell: u64,
) -> FpExt<F> {
    assert!(ell >= 2, "reduced_tate: ell must be >= 2");

    let p = &curve.p;
    let k = modulus.degree();

    // f = f_{ℓ,P}(Q) — one Miller call.
    let f = tate_pairing(curve, modulus, p_point, q_point, ell);

    // p^k − 1 exactly: for the fixture p = 47, k = 2 → 2208. Even for a one-limb p the exponent
    // (p^k − 1)/ℓ can exceed one limb (and p^k can exceed 128 bits), so no fixed-width integer
    // is used here.
    let p_big = p
        .as_words()
        .iter()
        .rev()
        .fold(BigUint::zero(), |acc, &w| (acc << 64u32) + BigUint::from(w));
    let k = u32::try_from(k).expect("extension degree fits u32");
    let pk_minus_1 = p_big.pow(k) - BigUint::one();
    let ell_big = BigUint::from(ell);

    // Precondition: ℓ | p^k − 1 (the embedding-degree condition).
    assert!(
        (&pk_minus_1 % &ell_big).is_zero(),
        "reduced_tate: ell = {ell} must divide p^k − 1 = {pk_minus_1} \
         (embedding-degree condition violated)"
    );

    let exp_big = pk_minus_1 / ell_big;

    // Convert the exponent to Uint<4> for FpExt::pow, refusing one that does not fit.
    let digits = exp_big.to_u64_digits();
    assert!(digits.len() <= 4, "reduced_tate: final exponent exceeds 256 bits");
    let mut words = [0u64; 4];
    words[..digits.len()].copy_from_slice(&digits);
    let exp = Uint::<4>::from_words(words);

    f.pow(&exp, modulus, p)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use shared_field::FpNaive;

    use crate::pairing::test_curves::pairing_toy;

    /// `tate_pairing(P, Q, ℓ)` returns a non-trivial element of `F_{p^k}*`.
    ///
    /// Verifies the raw Miller value is non-trivial before final exponentiation.
    #[test]
    fn tate_non_trivial() {
        let (curve, modulus, ell, p_point, q_point) = pairing_toy();
        let p = curve.p;
        let f = tate_pairing::<FpNaive<4>>(&curve, &modulus, &p_point, &q_point, ell);
        assert!(
            !f.is_one(&p),
            "tate_pairing(P, Q, ℓ) should be non-trivial (≠ 1)"
        );
    }

    /// `reduced_tate(Q, P, ℓ)^ℓ = 1` — the result lands in `μ_ℓ`.
    ///
    /// Uses Q as the first (Miller-base) argument and P as the second.  See the
    /// `test_curves` module docstring for why this argument order is required.
    #[test]
    fn reduced_tate_lands_in_mu_ell() {
        let (curve, modulus, ell, p_point, q_point) = pairing_toy();
        let p = curve.p;
        // Q as first argument, P as second — the non-degenerate direction.
        let tau = reduced_tate::<FpNaive<4>>(&curve, &modulus, &q_point, &p_point, ell);
        let tau_ell = tau.pow_u64(ell, &modulus, &p);
        assert!(
            tau_ell.is_one(&p),
            "reduced_tate(Q, P, ℓ)^ℓ should be 1 (result is in μ_ℓ)"
        );
    }

    /// `reduced_tate(Q, P, ℓ) ≠ 1` — the reduced Tate pairing is non-degenerate.
    ///
    /// Uses Q as the first (Miller-base) argument and P as the second.
    /// `P ∈ E(F_p)[ℓ]` is NOT in `ℓ·E(F_{p^k})`, so the Tate pairing
    /// `t_ℓ(Q, P)` is non-degenerate.  The reverse `t_ℓ(P, Q) = 1` because
    /// `Q ∈ ℓ·E(F_{p^k})` for this fixture.
    #[test]
    fn reduced_tate_non_degenerate() {
        let (curve, modulus, ell, p_point, q_point) = pairing_toy();
        let p = curve.p;
        // Q as first argument, P as second — the non-degenerate direction.
        let tau = reduced_tate::<FpNaive<4>>(&curve, &modulus, &q_point, &p_point, ell);
        assert!(
            !tau.is_one(&p),
            "reduced_tate(Q, P, ℓ) should be non-trivial (≠ 1) for independent Q, P"
        );
    }
}
