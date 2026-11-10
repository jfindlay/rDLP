//! MOV bridge and MOV/Frey–Rück reduction.
//!
//! # MOV bridge: `F_{p^k}` element → base-p `BigInt` encoding
//!
//! This module is the rho-side half of the MOV bridge. It converts a pairing output
//! [`FpExt<F>`] into the base-p [`BigInt`] encoding that [`gnfs::dl::solve_dl`] consumes
//! at k=2.
//!
//! The encoding convention belongs to gnfs. The rho side produces only the coefficient
//! `Vec<BigInt>` and delegates the encoding to the gnfs-side helper
//! [`gnfs::dl::ext::target::fpext_coeffs_to_dl_target`]. The rho side does NOT re-derive
//! the base-p encoding itself.
//!
//! # MOV/Frey–Rück reduction: `mov_reduce`
//!
//! [`mov_reduce`] is the MOV/Frey–Rück reduction entry. Given an ECDLP instance
//! `(curve, G, Q, R, ell)` — where `Q = k·G` in the order-ℓ subgroup and `R` is a
//! μ_ℓ-generator with `e(G, R) ≠ 1` — it transports the ECDLP into a discrete log in
//! `F_{p^k}*` via the pairing's bilinearity, encodes both pairing outputs through the
//! MOV bridge, and calls `gnfs::dl::solve_dl` to recover `k mod ℓ`.
//!
//! The mathematical shape: `e(Q, R) = e(k·G, R) = e(G, R)^k`. Writing `g := e(G, R)` and
//! `h := e(Q, R)`, both in `μ_ℓ ⊂ F_{p^k}*`, the ECDLP scalar `k` is exactly `log_g(h)`.
//!
//! # Modulus-consistency guard
//!
//! The gnfs-side helper asserts that the pairing's irreducible modulus matches what
//! `find_irreducible_degree2(p)` returns. This guard is mandatory: without it, a prime p
//! where `find_irreducible_degree2` picks a different irreducible than the pairing's modulus
//! would compute a DL in a different F_{p²} and return a wrong discrete log with no error.
//!
//! # Scale limits
//!
//! The `Uint<4> → BigInt` conversion assumes p < 2^64 (single limb), the domain of the
//! finite-field log engine behind it. [`mov_reduce`] rejects any larger `p` with
//! [`MovError::Unsupported`] before doing arithmetic on it, so no conversion on its path
//! truncates silently. Supporting larger fields would need the full limb vector.
//!
//! The recovered log is `k mod ℓ` (the subgroup log). A full composite-order lift via
//! Pohlig–Hellman is possible but not wired here.

use crypto_bigint::Uint;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::ToPrimitive;
use shared_field::Fp;
use shared_numth::is_prime;

use crate::curve::Curve;
use crate::pairing::ecext::PairingPoint;
use crate::pairing::fpext::{FpExt, IrreducibleModulus};
use crate::pairing::tate::reduced_tate;

/// Convert a `FpExt<F>` pairing output to the base-p `BigInt` encoding `solve_dl` consumes.
///
/// Extracts the coefficient vector from `elem` via `to_uint_vec`, converts each `Uint<4>` to
/// `BigInt` (toy-only: single-limb, p < 2^64), then delegates the base-p encoding and
/// modulus-consistency guard to the gnfs-side helper
/// [`gnfs::dl::ext::target::fpext_coeffs_to_dl_target`].
///
/// The single-limb check on the coefficients is a `debug_assert`, so it is absent from release
/// builds. [`mov_reduce`] does not rely on it: it rejects `p ≥ 2^64` with [`MovError::Unsupported`]
/// first, which keeps every reduced coefficient below `2^64`. A direct caller passing coefficients
/// with a non-zero upper limb gets a silently truncated encoding in a release build.
///
/// # Arguments
///
/// - `elem`: The pairing output element in F_{p^k}.
/// - `modulus`: The irreducible modulus the pairing used (carried by the pairing fixture).
/// - `p`: The base-field prime as `u64`.
///
/// # Panics
///
/// Panics if the pairing's modulus does not match what `find_irreducible_degree2(p)` returns
/// (forwarded from the gnfs-side modulus-consistency guard).
pub fn fpext_to_bigint<F: Fp<4>>(
    elem: &FpExt<F>,
    modulus: &IrreducibleModulus<F>,
    p: u64,
) -> BigInt {
    let p_big = BigInt::from(p);

    // Step 1: extract coefficient vector as Vec<Uint<4>>.
    let uint_vec: Vec<Uint<4>> = elem.to_uint_vec();

    // Step 2: convert each Uint<4> to BigInt.
    // Toy-scale limit: assumes p < 2^64 (single limb). A crypto-scale p needs the full
    // limb vector.
    let coeffs: Vec<BigInt> = uint_vec
        .iter()
        .map(|u| {
            let words = u.as_words();
            debug_assert!(
                words[1] == 0 && words[2] == 0 && words[3] == 0,
                "fpext_to_bigint: Uint<4> has non-zero upper limbs — p >= 2^64 is not supported \
                 (single-limb toy-scale limit)"
            );
            BigInt::from(words[0])
        })
        .collect();

    // Step 3: convert the pairing's IrreducibleModulus<F> coefficients to Vec<BigInt>.
    // The modulus is passed to the gnfs-side helper for the modulus-consistency guard.
    let modulus_bigint: Vec<BigInt> = modulus
        .coeffs
        .iter()
        .map(|c| {
            let uint_val = c.to_uint();
            let words = uint_val.as_words();
            debug_assert!(
                words[1] == 0 && words[2] == 0 && words[3] == 0,
                "fpext_to_bigint: modulus coefficient Uint<4> has non-zero upper limbs"
            );
            BigInt::from(words[0])
        })
        .collect();

    // Step 4: delegate to the gnfs-side helper (encoding convention + modulus-consistency guard).
    gnfs::dl::ext::target::fpext_coeffs_to_dl_target(coeffs, &p_big, &modulus_bigint)
}

// ─── MOV/Frey–Rück reduction ──────────────────────────────────────────────────

/// A verified MOV/Frey–Rück solution: `k·G == Q` in `E[ℓ]`, with the pairing values behind it.
///
/// `#[non_exhaustive]` keeps the fields readable everywhere while only this crate can construct a
/// solution, so an unverified answer cannot be built by a caller.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct MovSolution<F: Clone + PartialEq + Eq> {
    /// The logarithm modulo ℓ, in `[0, ℓ)`, satisfying `k_mod_ell·G == Q`.
    pub k_mod_ell: u64,
    /// The prime subgroup order ℓ.
    pub ell: u64,
    /// `g_pair = e(G, R) ∈ μ_ℓ`, the base of the finite-field log.
    pub g_pair: FpExt<F>,
    /// `h_pair = e(Q, R) ∈ μ_ℓ`, the target of the finite-field log: `g_pair^k = h_pair`.
    pub h_pair: FpExt<F>,
}

/// A precondition of the MOV instance that does not hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MovInputDefect {
    /// ℓ is not prime.
    EllNotPrime,
    /// ℓ does not divide `p^k − 1` for the extension degree `k` of the modulus.
    EllDoesNotDivideFieldOrder,
    /// The generator `G` is the point at infinity.
    GeneratorAtInfinity,
    /// The generator `G` does not satisfy the curve equation.
    GeneratorNotOnCurve,
    /// `ℓ·G` is not the point at infinity, so `G ∉ E[ℓ]`.
    GeneratorNotInTorsion,
    /// The target `Q` does not satisfy the curve equation.
    TargetNotOnCurve,
    /// `ℓ·Q` is not the point at infinity, so `Q ∉ E[ℓ]`.
    TargetNotInTorsion,
    /// The pairing point `R` does not satisfy the curve equation.
    PairingPointNotOnCurve,
    /// `ℓ·R` is not the point at infinity, so `R ∉ E[ℓ]`.
    PairingPointNotInTorsion,
}

/// Why [`mov_reduce`] did not return a verified logarithm.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MovError {
    /// The instance does not satisfy the reduction's preconditions.
    InvalidInput(MovInputDefect),
    /// A pairing value is unusable: `e(G, R) = 1`, or a Miller evaluation returned 0. A 0 arises
    /// when the evaluation point lies on the support of the Miller function's divisor, where the
    /// value is undefined: for `R = G` the first tangent-line numerator vanishes at the
    /// evaluation point (`mov_reduce_degenerate_r_is_degenerate_pairing`). `R` is then not a
    /// μ_ℓ-generator for `G` and the finite-field log carries no information about `k`. A property
    /// of the supplied `R`, not a solver failure.
    DegeneratePairing,
    /// The instance lies outside the reduction's supported domain: the field characteristic `p`
    /// has a non-zero word above the lowest (`p ≥ 2^64`), which the single-limb bridge to the
    /// finite-field log engine cannot represent. Reported before any arithmetic on `p`; says
    /// nothing about whether the instance is valid.
    Unsupported(&'static str),
    /// The finite-field log engine failed; its error is passed through unchanged.
    Dl(gnfs::dl::SolveDlError),
    /// The log engine returned a value that fails `k·G == Q` in `E[ℓ]`.
    ///
    /// An invariant failure (the pairing values and the points disagree): surfaced, never
    /// retried.
    Inconsistent {
        /// The candidate, reduced modulo ℓ, that failed verification.
        k_mod_ell: u64,
    },
}

impl From<gnfs::dl::SolveDlError> for MovError {
    fn from(e: gnfs::dl::SolveDlError) -> Self {
        MovError::Dl(e)
    }
}

/// `base^exp mod m` for `m ≥ 1`, without overflow (products of residues below 2^64 fit `u128`).
fn pow_mod(base: u64, exp: usize, m: u64) -> u64 {
    let m = u128::from(m);
    let mut acc: u128 = 1 % m;
    for _ in 0..exp {
        acc = acc * (u128::from(base) % m) % m;
    }
    u64::try_from(acc).expect("residue below a u64 modulus fits u64")
}

/// Validate the MOV instance: ℓ prime, `ℓ | p^k − 1`, and `G`, `Q`, `R` points of `E[ℓ]` (on the
/// curve with `ℓ·P = ∞`).
///
/// Checks run in the order of [`MovInputDefect`]'s variants; the first failure is returned. `G`
/// may not be ∞; `Q` and `R` may be (`Q = ∞` is `k = 0`; `R = ∞` is a degenerate pairing).
fn check_mov_input<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    g_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    r_point: &PairingPoint<F>,
    ell: u64,
) -> Result<(), MovError> {
    let defect = |d: MovInputDefect| Err(MovError::InvalidInput(d));
    let p = &curve.p;
    let k = modulus.degree();

    if !is_prime(&Uint::<4>::from(ell)) {
        return defect(MovInputDefect::EllNotPrime);
    }
    // ℓ | p^k − 1  ⇔  p^k ≡ 1 (mod ℓ). Reading `p` as its lowest word is exact: `mov_reduce`
    // rejects p ≥ 2^64 before calling this.
    debug_assert!(p.as_words()[1..].iter().all(|&w| w == 0), "check_mov_input: p >= 2^64");
    if pow_mod(p.as_words()[0] % ell, k, ell) != 1 {
        return defect(MovInputDefect::EllDoesNotDivideFieldOrder);
    }
    if g_point.is_infinity() {
        return defect(MovInputDefect::GeneratorAtInfinity);
    }

    let a = FpExt::from_base(F::from_uint(curve.a, p), k, p);
    let b = FpExt::from_base(F::from_uint(curve.b, p), k, p);
    let points = [
        (g_point, MovInputDefect::GeneratorNotOnCurve, MovInputDefect::GeneratorNotInTorsion),
        (q_point, MovInputDefect::TargetNotOnCurve, MovInputDefect::TargetNotInTorsion),
        (
            r_point,
            MovInputDefect::PairingPointNotOnCurve,
            MovInputDefect::PairingPointNotInTorsion,
        ),
    ];
    for (point, off_curve, off_torsion) in points {
        if !point.is_on_curve(&a, &b, modulus, p) {
            return defect(off_curve);
        }
        if !point.scalar_mul(ell, &a, modulus, p).is_infinity() {
            return defect(off_torsion);
        }
    }
    Ok(())
}

/// Reduce an ECDLP to a finite-field DLP via the MOV/Frey–Rück reduction.
///
/// Given the ECDLP `Q = k·G` in the order-ℓ subgroup of `E(F_p)`, and a μ_ℓ-generator
/// `R ∈ E(F_{p^k})[ℓ]` with `e(G, R) ≠ 1`, this function:
///
/// 1. Validates the instance (see *Errors*): ℓ prime, `ℓ | p^k − 1`, and `G`, `Q`, `R` points of
///    `E[ℓ]`. The reduced Tate pairing and its final exponentiation are only meaningful then.
/// 2. Computes `g_pair = reduced_tate(R, G, ell)` — the pairing `e(G, R) ∈ μ_ℓ` — and rejects
///    `g_pair ∈ {0, 1}` (R is not a μ_ℓ-generator with respect to G, so the DL is undefined).
/// 3. Computes `h_pair = reduced_tate(R, Q, ell)` — the pairing `e(Q, R) ∈ μ_ℓ` — rejecting 0.
/// 4. Encodes both pairing outputs as base-p `BigInt` via [`fpext_to_bigint`] (the MOV bridge).
/// 5. Calls `gnfs::dl::solve_dl(&g_big, &h_big, &p_big, 2, &ell_big)` to recover `k mod ℓ`.
/// 6. Reduces the result modulo ℓ and checks `k·G == Q` directly in `E[ℓ]`, independently of the
///    pairing and of the log engine, before returning it with the pairing values.
///
/// # Argument order for non-degeneracy
///
/// The reduced Tate pairing is non-degenerate when the **second** argument is not in
/// `ℓ·E(F_{p^k})`. For the `pairing_toy` fixture, `P ∈ E(F_p)[ℓ]` is not in `ℓ·E(F_{p^k})`,
/// so `reduced_tate(R, G, ell)` (R first, G second) is non-degenerate. The caller must supply
/// R such that this holds.
///
/// # Recovered log
///
/// The returned value is `k mod ℓ` — the discrete log in the order-ℓ subgroup, so the check
/// `k·G == Q` is exact in `E[ℓ]` with `k` reduced modulo ℓ. A full composite-order lift via
/// Pohlig–Hellman is possible but not wired here.
///
/// # Arguments
///
/// - `curve` — the short-Weierstrass curve `y² = x³ + ax + b`.
/// - `modulus` — the irreducible polynomial defining `F_{p^k}` (must match what
///   `gnfs::dl::ext::find_irreducible_degree2(p)` returns — enforced by the modulus-consistency
///   guard in [`fpext_to_bigint`]).
/// - `g_point` — the ECDLP generator `G ∈ E[ℓ]` (typically a base-field point).
/// - `q_point` — the ECDLP target `Q = k·G ∈ E[ℓ]`.
/// - `r_point` — a μ_ℓ-generator `R ∈ E(F_{p^k})[ℓ] \ E(F_p)` with `e(G, R) ≠ 1`.
/// - `ell` — the subgroup order ℓ (a prime dividing `p^k − 1`).
///
/// # Errors
///
/// - [`MovError::Unsupported`] if `p ≥ 2^64` (checked before anything else).
/// - [`MovError::InvalidInput`] with the first [`MovInputDefect`] found: ℓ not prime;
///   `ℓ ∤ p^k − 1`; `G = ∞`; or `G`, `Q`, `R` off the curve or outside `E[ℓ]`.
/// - [`MovError::DegeneratePairing`] if `R = ∞`, `e(G, R) = 1`, or a pairing evaluation is 0.
/// - [`MovError::Dl`] with the error of `solve_dl` unchanged (e.g. `InitSmoothingFailed`,
///   `DescentFailed`, or `Unsupported` for an extension degree above 2).
/// - [`MovError::Inconsistent`] if the log engine's value fails `k·G == Q`.
///
/// # Panics
///
/// Panics if the pairing's modulus does not match `find_irreducible_degree2(p)` (forwarded
/// from the modulus-consistency guard in [`fpext_to_bigint`]).
pub fn mov_reduce<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    g_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    r_point: &PairingPoint<F>,
    ell: u64,
) -> Result<MovSolution<F>, MovError> {
    // The domain guard comes first: every later single-limb conversion of `p` (input checks, the
    // final exponent, the bridge to the log engine) relies on it, in release builds as well.
    let p_uint = curve.p;
    let p_words = p_uint.as_words();
    if p_words[1..].iter().any(|&w| w != 0) {
        return Err(MovError::Unsupported("field characteristic p >= 2^64 (single-limb only)"));
    }
    let p_u64: u64 = p_words[0];

    check_mov_input(curve, modulus, g_point, q_point, r_point, ell)?;

    let k = modulus.degree();

    // e(R, ∞) = 1 and the Miller loop is undefined at a point at infinity, so handle R = ∞ and
    // Q = ∞ without evaluating the pairing.
    if r_point.is_infinity() {
        return Err(MovError::DegeneratePairing);
    }

    // g_pair = e(G, R) = reduced_tate(R, G, ell).
    // Argument order: R first (extension-field point), G second (base-field point).
    // This is the non-degenerate direction for the toy fixture (see `test_curves` module).
    //
    // A usable pairing value is a unit of F_{p^k}. The value 1 is the trivial pairing; the value
    // 0 is what the Miller evaluation returns when the evaluation point lies on the support of
    // the Miller function's divisor (R = G), where the pairing is undefined. Neither carries a
    // log.
    let g_pair: FpExt<F> = reduced_tate(curve, modulus, r_point, g_point, ell);
    if g_pair.is_one(&p_uint) || g_pair.is_zero(&p_uint) {
        return Err(MovError::DegeneratePairing);
    }

    // h_pair = e(Q, R) = reduced_tate(R, Q, ell); the pairing at Q = ∞ is the identity.
    let h_pair: FpExt<F> = if q_point.is_infinity() {
        FpExt::one(k, &p_uint)
    } else {
        reduced_tate(curve, modulus, r_point, q_point, ell)
    };
    if h_pair.is_zero(&p_uint) {
        return Err(MovError::DegeneratePairing);
    }

    // Encode both pairing outputs as base-p BigInt via the MOV bridge. The modulus-consistency
    // guard in fpext_to_bigint ensures the pairing's modulus matches what solve_dl uses
    // internally — the silent-wrong-field defense.
    let g_big: BigInt = fpext_to_bigint(&g_pair, modulus, p_u64);
    let h_big: BigInt = fpext_to_bigint(&h_pair, modulus, p_u64);

    let p_big = BigInt::from(p_u64);
    let ell_big = BigInt::from(ell);

    // Recover log_{g_pair}(h_pair) with the extension-field solver (degree 2 at toy scale).
    let k_big = gnfs::dl::solve_dl(&g_big, &h_big, &p_big, 2, &ell_big)?;

    // The subgroup log is defined modulo ℓ; reduce whatever representative the engine returned.
    let k_mod_ell: u64 = k_big
        .mod_floor(&ell_big)
        .to_u64()
        .expect("a residue modulo a u64 ℓ fits u64");

    if !log_matches(curve, modulus, g_point, q_point, k_mod_ell) {
        return Err(MovError::Inconsistent { k_mod_ell });
    }

    Ok(MovSolution { k_mod_ell, ell, g_pair, h_pair })
}

/// Check `k·G == Q` directly on the curve, independently of the pairing and the log engine.
fn log_matches<F: Fp<4>>(
    curve: &Curve,
    modulus: &IrreducibleModulus<F>,
    g_point: &PairingPoint<F>,
    q_point: &PairingPoint<F>,
    k: u64,
) -> bool {
    let p = &curve.p;
    let a = FpExt::from_base(F::from_uint(curve.a, p), modulus.degree(), p);
    g_point.scalar_mul(k, &a, modulus, p) == *q_point
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use num_bigint::BigInt;
    use shared_field::FpNaive;

    use crate::pairing::fpext::IrreducibleModulus;
    use crate::pairing::test_curves::{PAIRING_TOY_P, pairing_toy};

    // ── Fixture helpers ───────────────────────────────────────────────────────

    fn p47() -> Uint<4> {
        Uint::<4>::from(PAIRING_TOY_P)
    }

    fn modulus_k2() -> IrreducibleModulus<FpNaive<4>> {
        let p = p47();
        IrreducibleModulus::new(
            vec![
                FpNaive::<4>::from_u64(1, &p),
                FpNaive::<4>::from_u64(0, &p),
                FpNaive::<4>::from_u64(1, &p),
            ],
            &p,
        )
    }

    fn fp2(c0: u64, c1: u64) -> FpExt<FpNaive<4>> {
        let p = p47();
        FpExt { coeffs: vec![FpNaive::<4>::from_u64(c0, &p), FpNaive::<4>::from_u64(c1, &p)] }
    }

    // ── KAT 1: Round-trip ─────────────────────────────────────────────────────

    /// KAT: a known μ_ℓ element (23, 6) at p=47 → bridge → base-p BigInt = 23 + 6*47 = 305.
    ///
    /// This matches the dl_ext_kat.rs encoding: ζ = (23, 6) → BigInt(305).
    #[test]
    fn kat_round_trip_zeta() {
        let elem = fp2(23, 6);
        let modulus = modulus_k2();
        let result = fpext_to_bigint(&elem, &modulus, PAIRING_TOY_P);
        assert_eq!(
            result,
            BigInt::from(23 + 6 * 47i64),
            "round-trip: (23, 6) at p=47 should encode to 23 + 6*47 = 305"
        );
    }

    /// KAT: ζ² = (23, 41) at p=47 → bridge → base-p BigInt = 23 + 41*47 = 1950.
    ///
    /// Matches dl_ext_kat.rs: zeta_sq() = BigInt(1950).
    #[test]
    fn kat_round_trip_zeta_sq() {
        let elem = fp2(23, 41);
        let modulus = modulus_k2();
        let result = fpext_to_bigint(&elem, &modulus, PAIRING_TOY_P);
        assert_eq!(
            result,
            BigInt::from(23 + 41 * 47i64),
            "round-trip: (23, 41) at p=47 should encode to 23 + 41*47 = 1950"
        );
    }

    /// KAT: identity element (1, 0) at p=47 → bridge → BigInt(1).
    #[test]
    fn kat_round_trip_identity() {
        let elem = fp2(1, 0);
        let modulus = modulus_k2();
        let result = fpext_to_bigint(&elem, &modulus, PAIRING_TOY_P);
        assert_eq!(result, BigInt::from(1i64), "identity (1, 0) should encode to 1");
    }

    // ── KAT 2: Modulus-consistency guard fires ────────────────────────────────

    /// KAT: a wrong modulus causes a panic (modulus-consistency guard).
    ///
    /// The correct modulus for p=47 is u²+1 = [1, 0, 1]. We supply [1, 1, 1] (u²+u+1),
    /// which differs from what find_irreducible_degree2(47) returns, so the guard must panic.
    #[test]
    fn kat_modulus_guard_fires_on_wrong_modulus() {
        let p = p47();
        // Wrong modulus: u²+u+1 = [1, 1, 1] (not what find_irreducible_degree2(47) returns).
        let wrong_modulus = IrreducibleModulus::new(
            vec![
                FpNaive::<4>::from_u64(1, &p),
                FpNaive::<4>::from_u64(1, &p),
                FpNaive::<4>::from_u64(1, &p),
            ],
            &p,
        );
        let elem = fp2(23, 6);
        let result = std::panic::catch_unwind(|| {
            fpext_to_bigint(&elem, &wrong_modulus, PAIRING_TOY_P)
        });
        assert!(
            result.is_err(),
            "fpext_to_bigint with wrong modulus should panic (modulus-consistency guard)"
        );
    }

    // ── KAT 3: Uint<4> → BigInt step ─────────────────────────────────────────

    /// KAT: Uint<4> with value 23 converts to BigInt::from(23).
    ///
    /// Verifies the single-limb extraction step used in fpext_to_bigint.
    #[test]
    fn kat_uint4_to_bigint_single_limb() {
        let u = Uint::<4>::from(23u64);
        let words = u.as_words();
        // Upper limbs must be zero for toy-scale p < 2^64.
        assert_eq!(words[1], 0, "upper limb 1 should be zero");
        assert_eq!(words[2], 0, "upper limb 2 should be zero");
        assert_eq!(words[3], 0, "upper limb 3 should be zero");
        let b = BigInt::from(words[0]);
        assert_eq!(b, BigInt::from(23i64), "Uint<4>(23) should convert to BigInt(23)");
    }

    /// KAT: Uint<4> with value 0 converts to BigInt::from(0).
    #[test]
    fn kat_uint4_to_bigint_zero() {
        let u = Uint::<4>::from(0u64);
        let words = u.as_words();
        let b = BigInt::from(words[0]);
        assert_eq!(b, BigInt::from(0i64), "Uint<4>(0) should convert to BigInt(0)");
    }

    /// KAT: Uint<4> with value 46 (= p-1 for p=47) converts to BigInt::from(46).
    #[test]
    fn kat_uint4_to_bigint_p_minus_1() {
        let u = Uint::<4>::from(46u64);
        let words = u.as_words();
        let b = BigInt::from(words[0]);
        assert_eq!(b, BigInt::from(46i64), "Uint<4>(46) should convert to BigInt(46)");
    }

    // ── Verification and input checks ─────────────────────────────────────────

    /// `log_matches` accepts exactly the true residue of `k` modulo ℓ, so a wrong log-engine
    /// value would be reported as inconsistent.
    #[test]
    fn log_matches_accepts_only_the_true_log() {
        let (curve, modulus, _ell, g, _r) = pairing_toy();
        let p = curve.p;
        let a = FpExt::from_base(FpNaive::<4>::from_u64(1, &p), 2, &p);
        let q = g.scalar_mul(2, &a, &modulus, &p);
        assert!(log_matches(&curve, &modulus, &g, &q, 2));
        assert!(!log_matches(&curve, &modulus, &g, &q, 0));
        assert!(!log_matches(&curve, &modulus, &g, &q, 1));
    }

    /// `pow_mod` agrees with direct evaluation, including the degenerate modulus 1.
    #[test]
    fn pow_mod_matches_direct_power() {
        assert_eq!(pow_mod(47, 2, 3), 2209 % 3);
        assert_eq!(pow_mod(47, 2, 5), 2209 % 5);
        assert_eq!(pow_mod(7, 0, 5), 1);
        assert_eq!(pow_mod(7, 3, 1), 0);
        // Residues near 2^64 must not overflow.
        assert_eq!(pow_mod(u64::MAX - 1, 2, u64::MAX), 1);
    }

    // ── Final exponent above one limb ─────────────────────────────────────────

    /// `p = 68719476767`, `p ≡ 11 (mod 12)` and just above `2^36`, so `(p² − 1)/3 > 2^64`: the
    /// final-exponentiation exponent does not fit one limb.
    const LARGE_P: u64 = 68_719_476_767;

    /// The curve `y² = x³ + 1` over `F_LARGE_P` (supersingular, `#E = p + 1`, divisible by 3) with
    /// `F_{p²} = F_p[u]/(u² + 1)`, `G = (0, 1) ∈ E(F_p)[3]`, and
    /// `R = (x₀ω, 2ω + 1) ∈ E(F_{p²})[3] \ E(F_p)`.
    ///
    /// Here `ω` is a primitive cube root of unity and `x₀³ = −4`, so `R` has `y² = −3 = x³ + 1`
    /// and `x(R)` is a root of the 3-division polynomial `x(x³ + 4)`. Coordinates computed offline
    /// and checked with `3·R = ∞`.
    fn large_exponent_fixture() -> (
        Curve,
        IrreducibleModulus<FpNaive<4>>,
        PairingPoint<FpNaive<4>>,
        PairingPoint<FpNaive<4>>,
    ) {
        let p = Uint::<4>::from(LARGE_P);
        let c = |v: u64| FpNaive::<4>::from_u64(v, &p);
        let curve = Curve {
            p,
            a: Uint::<4>::from(0u64),
            b: Uint::<4>::from(1u64),
            n: Uint::<4>::from(3u64),
            gx: Uint::<4>::from(0u64),
            gy: Uint::<4>::from(1u64),
        };
        let modulus = IrreducibleModulus::new(vec![c(1), c(0), c(1)], &p);
        let g = PairingPoint::new(
            FpExt::from_base(c(0), 2, &p),
            FpExt::from_base(c(1), 2, &p),
        );
        let r = PairingPoint::new(
            FpExt { coeffs: vec![c(55_629_843_007), c(11_542_864_339)] },
            FpExt { coeffs: vec![c(0), c(32_645_870_166)] },
        );
        (curve, modulus, g, r)
    }

    /// With a final exponent above `2^64` the reduced Tate pairing lands in `μ_3` and is bilinear.
    ///
    /// A final exponent cut down to its low word lands outside `μ_3` almost surely, so
    /// `e³ = 1` and `e(R, 2G) = e(R, G)²` fail unless the whole exponent is applied.
    #[test]
    fn reduced_tate_with_exponent_above_one_limb_lands_in_mu_ell() {
        let (curve, modulus, g, r) = large_exponent_fixture();
        let p = curve.p;
        let ell = 3;
        assert!(
            (u128::from(LARGE_P) * u128::from(LARGE_P) - 1) / u128::from(ell)
                > u128::from(u64::MAX),
            "fixture assumption: the final exponent exceeds 2^64"
        );
        let a = FpExt::from_base(FpNaive::<4>::from_u64(0, &p), 2, &p);
        let b = FpExt::from_base(FpNaive::<4>::from_u64(1, &p), 2, &p);
        assert!(g.is_on_curve(&a, &b, &modulus, &p) && r.is_on_curve(&a, &b, &modulus, &p));
        assert!(g.scalar_mul(ell, &a, &modulus, &p).is_infinity());
        assert!(r.scalar_mul(ell, &a, &modulus, &p).is_infinity());

        let e_g = reduced_tate::<FpNaive<4>>(&curve, &modulus, &r, &g, ell);
        assert!(!e_g.is_one(&p) && !e_g.is_zero(&p), "e(G, R) is a nontrivial value");
        assert!(e_g.pow_u64(ell, &modulus, &p).is_one(&p), "e(G, R)^3 = 1");

        let g2 = g.scalar_mul(2, &a, &modulus, &p);
        let e_g2 = reduced_tate::<FpNaive<4>>(&curve, &modulus, &r, &g2, ell);
        assert_eq!(e_g2, e_g.square(&modulus, &p), "e(2G, R) = e(G, R)^2");
    }

    // ── KAT 4: pairing_toy fixture integration ────────────────────────────────

    /// KAT: the bridge is consistent with the pairing_toy fixture modulus.
    ///
    /// Uses the actual pairing_toy() modulus (not a hand-constructed one) to verify
    /// the bridge works end-to-end with the fixture.
    #[test]
    fn kat_pairing_toy_modulus_consistent() {
        let (_curve, modulus, _ell, _p_point, _q_point) = pairing_toy();
        let elem = fp2(23, 6);
        let result = fpext_to_bigint(&elem, &modulus, PAIRING_TOY_P);
        // ζ = (23, 6) → 23 + 6*47 = 305.
        assert_eq!(
            result,
            BigInt::from(305i64),
            "pairing_toy modulus: (23, 6) should encode to 305"
        );
    }
}
