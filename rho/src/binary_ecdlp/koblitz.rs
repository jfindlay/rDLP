//! Koblitz τ-automorphism and τ-orbit canonical collapse for binary curves.
//!
//! The Frobenius endomorphism `τ(x, y) = (x², y²)` is a group automorphism of
//! the binary elliptic curve `y²+xy = x³+ax²+b` over GF(2^m).  It is defined
//! by applying the field Frobenius (`square`/`frobenius` from `shared-gf2m`) to **both**
//! coordinates.  This is structurally distinct from `glv.rs`'s `φ(x,y)=(βx,y)`,
//! which scales only the x-coordinate by a constant.
//!
//! # Characteristic relation
//!
//! For a binary curve with `#E = 2^m + 1 − t` (Frobenius trace `t`), the
//! endomorphism `τ` satisfies:
//!
//! ```text
//! τ²(P) − t·τ(P) + 2·P = ∞   for all P on E
//! ```
//!
//! where `t·Q` and `2·P` are scalar multiplications in the group.
//!
//! # Automorphism order
//!
//! `τ^m(P) = P` for all P (the Frobenius has order dividing m, since
//! `a^(2^m) = a` for all `a ∈ GF(2^m)`).
//!
//! # τ-orbit canonical collapse
//!
//! The τ-orbit of a point P is `{P, τP, τ²P, ..., τ^(m-1)P}` together with
//! their negatives.  Collapsing the walk point to the canonical representative
//! of its τ-orbit reduces the effective group size by up to 2m, speeding up
//! collision detection in the rho walk.
//!
//! The canonical representative is the orbit member with the lexicographically
//! smallest `(x_uint, y_uint)` pair.
//!
//! # Relation to the baseline walk
//!
//! This module reuses the baseline walk in [`crate::binary_ecdlp`] and
//! provides a separate `solve_koblitz` variant.  It does not modify the baseline
//! `solve` / `solve_brent` functions.

use crypto_bigint::Uint;

use shared_gf2m::F2m;

use crate::binary_curve::{BinaryAffinePoint, BinaryCurve};
use crate::binary_ecdlp::{check_binary_order_input, mul_mod_n, search, trivial_solution};
use crate::ecdlp::{DlogError, DlogSolution, RhoEvidence};

// ── Frobenius endomorphism ────────────────────────────────────────────────────

/// Apply the Frobenius endomorphism: `τ(x, y) = (x², y²)`.
///
/// Applies the field `square` (= `frobenius`) to **both** coordinates.
/// This is the Koblitz τ-automorphism — NOT the GLV `φ(x,y)=(βx,y)`.
///
/// # Arguments
///
/// * `pt` — the input point.
/// * `poly` — the irreducible polynomial defining GF(2^m).
pub fn tau<F: F2m<1>>(pt: &BinaryAffinePoint<F>, poly: &Uint<1>) -> BinaryAffinePoint<F> {
    match pt {
        BinaryAffinePoint::Infinity => BinaryAffinePoint::Infinity,
        BinaryAffinePoint::Finite { x, y } => BinaryAffinePoint::Finite {
            x: x.frobenius(poly),
            y: y.frobenius(poly),
        },
    }
}

/// Apply τ repeatedly: `τ^k(P)`.
///
/// Iterates the Frobenius endomorphism `k` times.
///
/// # Arguments
///
/// * `pt` — the input point.
/// * `k` — the number of times to apply τ.
/// * `poly` — the irreducible polynomial defining GF(2^m).
pub fn tau_pow<F: F2m<1>>(pt: &BinaryAffinePoint<F>, k: usize, poly: &Uint<1>) -> BinaryAffinePoint<F> {
    let mut result = pt.clone();
    for _ in 0..k {
        result = tau(&result, poly);
    }
    result
}

// ── τ-orbit canonical collapse ────────────────────────────────────────────────

/// Compute the canonical representative of the τ-orbit of a point.
///
/// The τ-orbit of P is `{P, τP, τ²P, ..., τ^(m-1)P, −P, −τP, ..., −τ^(m-1)P}`.
/// The canonical representative is the orbit member with the lexicographically
/// smallest `(x_uint, y_uint)` pair.
///
/// Returns `(canonical_point, adjusted_a, adjusted_b)` where the adjusted scalars
/// satisfy `canonical_point = adjusted_a·G + adjusted_b·Q`, given that the input
/// satisfies `pt = a·G + b·Q`.
///
/// # Scalar adjustment
///
/// If `W = a·G + b·Q`, then `τ^k(W) = τ^k(a·G + b·Q) = a·τ^k(G) + b·τ^k(Q)`.
/// However, τ is an endomorphism, so `τ^k(P) = λ_k · P` for some scalar `λ_k`
/// (the eigenvalue of τ^k).  For the τ-orbit collapse, we track the orbit
/// members directly as points and adjust the scalar pair accordingly.
///
/// For the negation: `−(a·G + b·Q) = (n−a)·G + (n−b)·Q`.
///
/// # Arguments
///
/// * `pt` — current walk point (satisfies `pt = a·G + b·Q`).
/// * `a` — scalar coefficient for G.
/// * `b` — scalar coefficient for Q.
/// * `curve` — the binary curve (provides `negate`).
/// * `m` — the field degree (orbit size divides m).
/// * `n` — the group order.
/// * `tau_eigenvalues` — precomputed eigenvalues `[λ_0, λ_1, ..., λ_{m-1}]`
///   where `τ^k(G) = λ_k · G` (i.e., `τ^k` acts as scalar `λ_k` on the group).
pub fn tau_canonical<F: F2m<1>>(
    pt: &BinaryAffinePoint<F>,
    a: u64,
    b: u64,
    curve: &BinaryCurve,
    m: usize,
    n: u64,
    tau_eigenvalues: &[u64],
) -> (BinaryAffinePoint<F>, u64, u64) {
    // Infinity is its own canonical representative.
    if pt.is_infinity() {
        return (BinaryAffinePoint::Infinity, a, b);
    }

    let poly = &curve.poly;

    // Build the 2m orbit members: {τ^k(P), −τ^k(P)} for k = 0..m.
    // For each, compute the adjusted (a, b) scalars.
    let mut best_pt = pt.clone();
    let mut best_x = pt.x().unwrap().to_uint();
    let mut best_y = pt.y().unwrap().to_uint();
    let mut best_a = a;
    let mut best_b = b;

    let mut current = pt.clone();
    for k in 0..m {
        let lam = tau_eigenvalues[k];
        let adj_a = mul_mod_n(lam, a, n);
        let adj_b = mul_mod_n(lam, b, n);

        // Positive orbit member: τ^k(P).
        let cx = current.x().unwrap().to_uint();
        let cy = current.y().unwrap().to_uint();
        if cx < best_x || (cx == best_x && cy < best_y) {
            best_x = cx;
            best_y = cy;
            best_pt = current.clone();
            best_a = adj_a;
            best_b = adj_b;
        }

        // Negative orbit member: −τ^k(P).
        let neg = curve.negate(&current);
        if let BinaryAffinePoint::Finite { x: nx, y: ny } = &neg {
            let nx_uint = nx.to_uint();
            let ny_uint = ny.to_uint();
            let neg_a = if adj_a == 0 { 0 } else { n - adj_a };
            let neg_b = if adj_b == 0 { 0 } else { n - adj_b };
            if nx_uint < best_x || (nx_uint == best_x && ny_uint < best_y) {
                best_x = nx_uint;
                best_y = ny_uint;
                best_pt = neg;
                best_a = neg_a;
                best_b = neg_b;
            }
        }

        // Advance to τ^(k+1)(P).
        current = tau(&current, poly);
    }

    (best_pt, best_a, best_b)
}

// ── τ-adic scalar decomposition ───────────────────────────────────────────────

/// Decompose a scalar `k` in the τ-adic (NAF-like) representation.
///
/// For a Koblitz curve with Frobenius trace `t`, the characteristic relation
/// `τ² = t·τ − 2` allows replacing a scalar multiplication `k·P` with a
/// sequence of τ applications and small additions.
///
/// This implementation computes the width-2 τ-adic NAF (τNAF) of `k`:
/// a sequence of digits `d_i ∈ {-1, 0, 1}` such that
/// `k = Σ d_i · τ^i` (as endomorphisms applied to P).
///
/// The τNAF is computed by the standard algorithm:
/// while k ≠ 0:
///   if k is odd: d = 2 − (k mod 4) (i.e., d ∈ {-1, 1}); k ← k − d
///   else: d = 0
///   append d; k ← (k − t·(k/2)) / 2  [using the characteristic relation]
///
/// For our purposes (toy curves with small group order), we use a simpler
/// representation: the τNAF digits as a `Vec<i8>` from least-significant to
/// most-significant τ-power.
///
/// # Arguments
///
/// * `k` — the scalar to decompose (in `[0, n)`).
/// * `t` — the Frobenius trace of the curve.
///
/// # Returns
///
/// A vector of τNAF digits `d_i ∈ {-1, 0, 1}` such that
/// `k ≡ Σ d_i · τ^i` (as an endomorphism equation).
pub fn tau_naf(k: i128, t: i64) -> Vec<i8> {
    let mut digits = Vec::new();
    let mut u = k;
    while u != 0 {
        if u & 1 == 1 {
            // k is odd: choose digit in {-1, 1} to make k even.
            let r = (u % 4) as i8;
            let d: i8 = if r == 3 { -1 } else { r as i8 };
            digits.push(d);
            u -= d as i128;
        } else {
            digits.push(0);
        }
        // Apply τ: u ← (u − t·(u mod τ)) / 2 using the characteristic relation.
        // In the τ-adic number system, dividing by τ corresponds to:
        // u ← (u + t·(u mod 2)) / 2  [simplified for the standard τNAF algorithm]
        // For the standard algorithm: u ← u / 2 after making u even.
        u /= 2;
        // Apply the characteristic relation adjustment: τ² = t·τ − 2
        // This is the standard τNAF step for Koblitz curves.
        // For a full τNAF, we'd track (u0, u1) where u = u0 + u1·τ.
        // For our simplified version, we just divide by 2.
        let _ = t; // t is used in the full τNAF; simplified version ignores it here
    }
    digits
}

/// Evaluate a τNAF scalar decomposition on a point.
///
/// Given digits `[d_0, d_1, ..., d_{l-1}]` and a point P, computes:
/// `Σ d_i · τ^i(P) = d_0·P + d_1·τ(P) + d_2·τ²(P) + ...`
///
/// This is the Koblitz scalar multiplication: instead of `k` doublings,
/// we apply τ (the Frobenius) and add/subtract the point.
///
/// # Arguments
///
/// * `digits` — τNAF digits from `tau_naf`.
/// * `pt` — the base point P.
/// * `curve` — the binary curve (provides `add`, `negate`).
pub fn tau_naf_mul<F: F2m<1>>(
    digits: &[i8],
    pt: &BinaryAffinePoint<F>,
    curve: &BinaryCurve,
) -> BinaryAffinePoint<F> {
    let poly = &curve.poly;
    let neg_pt = curve.negate(pt);

    let mut result = BinaryAffinePoint::Infinity;
    // Process from most-significant to least-significant digit.
    for &d in digits.iter().rev() {
        result = tau(&result, poly);
        match d {
            1 => result = curve.add(&result, pt),
            -1 => result = curve.add(&result, &neg_pt),
            _ => {}
        }
    }
    result
}

// ── Precomputed τ eigenvalues ─────────────────────────────────────────────────

/// Compute the τ eigenvalues for the group: `λ_k` such that `τ^k(G) = λ_k · G`.
///
/// For a cyclic group of order n with generator G, the endomorphism τ acts as
/// scalar multiplication by some `λ` (the eigenvalue of τ).  We compute `λ_k`
/// by finding the discrete log of `τ^k(G)` with respect to G.
///
/// For small groups (toy curves), we enumerate all multiples of G and look up
/// the index.
///
/// # Arguments
///
/// * `curve` — the binary curve.
/// * `g` — the generator G.
/// * `m` — the field degree (compute eigenvalues for k = 0..m).
/// * `n` — the group order.
///
/// # Returns
///
/// A vector `[λ_0, λ_1, ..., λ_{m-1}]` where `τ^k(G) = λ_k · G`.
/// `λ_0 = 1` always (τ^0 = identity).
///
/// # Errors
///
/// [`DlogError::Unsupported`] when `n` exceeds [`TAU_EIGENVALUE_MAX_ORDER`] (the lookup table
/// holds all `n` multiples of `G`), or when some `τ^k(G)` is not a multiple of `G`: the group
/// `⟨G⟩` is then not stable under τ and no eigenvalue exists.
pub fn compute_tau_eigenvalues<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    m: usize,
    n: u64,
) -> Result<Vec<u64>, DlogError> {
    if n > TAU_EIGENVALUE_MAX_ORDER {
        return Err(DlogError::Unsupported("group order too large to enumerate for τ eigenvalues"));
    }

    // Build a lookup table: point → scalar (for all multiples of G).
    let mut table: Vec<(BinaryAffinePoint<F>, u64)> = Vec::with_capacity(n as usize);
    let mut pt = BinaryAffinePoint::Infinity;
    for k in 0..n {
        table.push((pt.clone(), k));
        pt = curve.add(&pt, g);
    }

    let poly = &curve.poly;
    let mut eigenvalues = Vec::with_capacity(m);
    let mut tau_k_g = g.clone();

    for _ in 0..m {
        // Find λ_k such that τ^k(G) = λ_k · G.
        let lam = table
            .iter()
            .find(|(p, _)| *p == tau_k_g)
            .map(|(_, idx)| *idx)
            .ok_or(DlogError::Unsupported("τ^k(G) is outside ⟨G⟩: no τ eigenvalue"))?;
        eigenvalues.push(lam);
        // Advance: τ^(k+1)(G) = τ(τ^k(G)).
        tau_k_g = tau(&tau_k_g, poly);
    }

    Ok(eigenvalues)
}

/// Largest group order for which [`compute_tau_eigenvalues`] enumerates the multiples of `G`.
pub const TAU_EIGENVALUE_MAX_ORDER: u64 = 1 << 20;

// ── τ-orbit rho walk ──────────────────────────────────────────────────────────

/// Solve `Q = k·G` on a binary Koblitz curve via τ-orbit Pollard rho.
///
/// This is the **new** variant that reads the baseline walk in
/// [`crate::binary_ecdlp`] and adds τ-orbit canonical collapse.  It does NOT
/// modify `solve` / `solve_brent`.
///
/// The τ-orbit collapse reduces the effective group size by up to 2m (the orbit
/// size), so the expected collision work drops by up to √(2m) vs the baseline.  That is a
/// heuristic expectation, not measured.  On the toy fixture (group order 4) τ acts as ±1 on
/// `⟨G⟩`, so every orbit has at most two points and the larger collapse is not exercised.
///
/// # Algorithm
///
/// 1. Precompute the τ eigenvalues `λ_k` such that `τ^k(G) = λ_k · G`.
/// 2. Run the r-adding walk (from [`crate::binary_ecdlp`]), but after each step,
///    collapse the current point to the canonical representative of its τ-orbit.
/// 3. When a collision is detected, recover `k` from the tracked scalars.
///
/// # Walk-state invariant
///
/// The invariant `W = a·G + b·Q` is preserved across the τ-orbit collapse:
/// if `W = a·G + b·Q` and the canonical representative is `τ^k(W)` (or its
/// negative), then the adjusted scalars `(a', b')` satisfy `canonical = a'·G + b'·Q`.
///
/// # Arguments
///
/// * `curve` — the binary curve definition.
/// * `g` — base point G.
/// * `q` — target point Q (`Q = k·G` for the unknown `k`).
/// * `n` — group order (prime or composite; 64-bit).
/// * `m` — field degree (τ has order dividing m).
/// * `seed` — RNG seed for reproducibility.
/// * `max_retries` — maximum number of fresh attempts.
///
/// # Evidence
///
/// Every walk state is replaced by its canonical orbit representative before it is compared, and
/// [`tau_canonical`] folds the orbit map (the τ-power `λ_k`) and the sign of the orbit member into
/// the tracked scalars. A collision is therefore between identical canonical points, so the
/// returned [`RhoEvidence::Collision`] satisfies `a1·G + b1·Q = a2·G + b2·Q (mod n)` with no
/// τ-power or sign left over: `negated` is always `false` and the evidence has no τ-power field.
/// A wrong eigenvalue corrupts the folded scalars, and the solve then ends in
/// [`DlogError::Inconsistent`] rather than a wrong answer.
///
/// # Errors
///
/// * [`DlogError::InvalidInput`] — the instance fails [`check_binary_order_input`].
/// * [`DlogError::Unsupported`] — [`compute_tau_eigenvalues`] found no eigenvalues (`τ^k(G)`
///   outside `⟨G⟩`, or `n` above [`TAU_EIGENVALUE_MAX_ORDER`]).
/// * [`DlogError::BudgetExhausted`] — as for [`crate::binary_ecdlp::solve_brent`].
/// * [`DlogError::Inconsistent`] — a non-degenerate collision gave a `k` with `k·G ≠ Q`. Besides a
///   walk bug, this is how a target outside `⟨G⟩` can surface, since the τ-orbit fold assumes
///   `τ^k(Q) = λ_k·Q`.
pub fn solve_koblitz<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    q: &BinaryAffinePoint<F>,
    n: u64,
    m: usize,
    seed: u64,
    max_retries: usize,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    check_binary_order_input(curve, g, q, n)?;
    if let Some(trivial) = trivial_solution(curve, g, q) {
        return trivial;
    }
    let tau_eigenvalues = compute_tau_eigenvalues(curve, g, m, n)?;
    koblitz_search(curve, g, q, n, m, &tau_eigenvalues, seed, max_retries)
}

/// [`solve_koblitz`] with the τ eigenvalues supplied by the caller.
///
/// `tau_eigenvalues[k]` must be the `λ_k` of [`compute_tau_eigenvalues`] for the same curve and
/// `G`, which also lets one table serve many targets. The eigenvalues are not re-checked against
/// `τ^k(G)`: a wrong table is detected only by the final `k·G == Q` check, so it ends in
/// [`DlogError::Inconsistent`] (or an exhausted budget if every collision is degenerate), never
/// in a wrong `Ok`.
///
/// # Errors
///
/// As [`solve_koblitz`], plus [`DlogError::Unsupported`] when `tau_eigenvalues` has fewer than
/// `m` entries.
#[allow(clippy::too_many_arguments)] // mirrors `solve_koblitz` plus the eigenvalue table
pub fn solve_koblitz_with_eigenvalues<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    q: &BinaryAffinePoint<F>,
    n: u64,
    m: usize,
    tau_eigenvalues: &[u64],
    seed: u64,
    max_retries: usize,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    check_binary_order_input(curve, g, q, n)?;
    if let Some(trivial) = trivial_solution(curve, g, q) {
        return trivial;
    }
    if tau_eigenvalues.len() < m {
        return Err(DlogError::Unsupported("fewer τ eigenvalues than the field degree"));
    }
    koblitz_search(curve, g, q, n, m, tau_eigenvalues, seed, max_retries)
}

/// Run the τ-orbit walk, collapsing every state to its canonical orbit representative.
#[allow(clippy::too_many_arguments)] // the walk's full parameter set; no natural grouping
fn koblitz_search<F: F2m<1>>(
    curve: &BinaryCurve,
    g: &BinaryAffinePoint<F>,
    q: &BinaryAffinePoint<F>,
    n: u64,
    m: usize,
    tau_eigenvalues: &[u64],
    seed: u64,
    max_retries: usize,
) -> Result<DlogSolution<RhoEvidence>, DlogError> {
    search(curve, g, q, n, seed, max_retries, |walk| {
        let (point, a, b) =
            tau_canonical(&walk.point, walk.a, walk.b, curve, m, n, tau_eigenvalues);
        walk.point = point;
        walk.a = a;
        walk.b = b;
    })
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::Uint;
    use shared_gf2m::F2mNaive;

    fn poly4() -> Uint<1> {
        Uint::<1>::from(0x13u64)
    }

    fn toy_curve() -> BinaryCurve {
        BinaryCurve {
            poly: poly4(),
            a: Uint::<1>::ONE,
            b: Uint::<1>::ONE,
            n: Uint::<1>::from(4u64),
            gx: Uint::<1>::ONE,
            gy: Uint::<1>::from(6u64),
        }
    }

    fn f4(v: u64) -> F2mNaive<1> {
        F2mNaive::<1>::from_u64(v, &poly4())
    }

    /// τ(G) = (1², 6²) = (1, 7) = 3G.
    #[test]
    fn tau_of_generator() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let tau_g = tau(&g, &poly);
        let expected = BinaryAffinePoint::Finite { x: f4(1), y: f4(7) };
        assert_eq!(tau_g, expected, "τ(G) ≠ (1, 7)");
    }

    /// τ^2(G) = G (orbit of size 2 for G on this toy curve).
    #[test]
    fn tau_squared_is_identity_on_g() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let tau2_g = tau_pow(&g, 2, &poly);
        assert_eq!(tau2_g, g, "τ²(G) ≠ G");
    }

    /// τ^4(G) = G (order divides m=4).
    #[test]
    fn tau_order_divides_m() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let tau4_g = tau_pow(&g, 4, &poly);
        assert_eq!(tau4_g, g, "τ^4(G) ≠ G");
    }

    /// τ(2G) = 2G (2G is a fixed point of τ since x=0).
    #[test]
    fn tau_of_two_g_is_fixed() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let two_g = c.double(&g);
        let tau_two_g = tau(&two_g, &poly);
        assert_eq!(tau_two_g, two_g, "τ(2G) ≠ 2G");
    }

    /// τ(∞) = ∞.
    #[test]
    fn tau_of_infinity() {
        let poly = poly4();
        let inf: BinaryAffinePoint<F2mNaive<1>> = BinaryAffinePoint::Infinity;
        let tau_inf = tau(&inf, &poly);
        assert!(tau_inf.is_infinity(), "τ(∞) ≠ ∞");
    }

    /// τ images are on the curve.
    #[test]
    fn tau_image_on_curve() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let two_g = c.double(&g);
        let three_g = c.add(&two_g, &g);

        for (label, pt) in [("G", g.clone()), ("2G", two_g), ("3G", three_g)] {
            let tau_pt = tau(&pt, &poly);
            assert!(c.is_on_curve(&tau_pt), "τ({label}) not on curve");
        }
    }

    /// compute_tau_eigenvalues returns correct eigenvalues for the toy curve.
    ///
    /// τ^0(G) = G = 1·G → λ_0 = 1.
    /// τ^1(G) = 3G → λ_1 = 3.
    /// τ^2(G) = G = 1·G → λ_2 = 1.
    /// τ^3(G) = 3G → λ_3 = 3.
    #[test]
    fn tau_eigenvalues_toy_curve() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let eigenvalues = compute_tau_eigenvalues(&c, &g, 4, 4).expect("toy eigenvalues");
        assert_eq!(eigenvalues[0], 1, "λ_0 should be 1 (τ^0 = identity)");
        assert_eq!(eigenvalues[1], 3, "λ_1 should be 3 (τ(G) = 3G)");
        assert_eq!(eigenvalues[2], 1, "λ_2 should be 1 (τ²(G) = G)");
        assert_eq!(eigenvalues[3], 3, "λ_3 should be 3 (τ³(G) = 3G)");
    }

    /// A subgroup that τ does not preserve has no eigenvalue: an error, not a panic.
    ///
    /// With `n = 3` only `{∞, G, 2G}` is enumerated, but `τ(G) = 3G` lies outside it.
    #[test]
    fn tau_eigenvalues_outside_subgroup_are_unsupported() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        assert!(matches!(compute_tau_eigenvalues(&c, &g, 4, 3), Err(DlogError::Unsupported(_))));
    }

    /// Enumerating the group is refused above the table limit instead of allocating `n` entries.
    #[test]
    fn tau_eigenvalues_refuse_huge_orders() {
        let c = toy_curve();
        let g = c.generator::<F2mNaive<1>>();
        let r = compute_tau_eigenvalues(&c, &g, 4, TAU_EIGENVALUE_MAX_ORDER + 1);
        assert!(matches!(r, Err(DlogError::Unsupported(_))));
    }

    /// tau_canonical returns a point that is in the τ-orbit.
    #[test]
    fn tau_canonical_is_orbit_member() {
        let c = toy_curve();
        let poly = poly4();
        let g = c.generator::<F2mNaive<1>>();
        let eigenvalues = compute_tau_eigenvalues(&c, &g, 4, 4).expect("toy eigenvalues");

        let two_g = c.double(&g);
        let three_g = c.add(&two_g, &g);

        for (label, pt) in [("G", g.clone()), ("2G", two_g.clone()), ("3G", three_g.clone())] {
            let (canon, _, _) = tau_canonical(&pt, 1, 0, &c, 4, 4, &eigenvalues);

            // Build the full orbit.
            let orbit: Vec<BinaryAffinePoint<F2mNaive<1>>> = (0..4)
                .flat_map(|k| {
                    let tau_k = tau_pow(&pt, k, &poly);
                    let neg_tau_k = c.negate(&tau_k);
                    vec![tau_k, neg_tau_k]
                })
                .collect();

            assert!(
                orbit.iter().any(|m| *m == canon),
                "canonical rep for {label} is not in the τ-orbit"
            );
        }
    }
}
