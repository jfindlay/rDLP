//! Point decomposition via the Semaev polynomial for the index-calculus ECDLP solver.
//!
//! This module implements the Semaev-based point-decomposition step: given a point `Q`,
//! find a sum of `m` factor-base points equal to `Q` by finding roots of the Semaev
//! polynomial `S_{m+1}` specialised at `Q`'s x-coordinate.
//!
//! # Algorithm (m = 2, the ratified arity)
//!
//! 1. Compute `S_3(x_0, x_1, x_2)` — the 3-variable Semaev polynomial.
//! 2. Specialise at `Q`'s x-coordinate: `S_3(Q.x, x_1, x_2)` via `partial_eval`.
//! 3. For each factor-base point `P_i` (x-coordinate `x_i`), substitute `x_1 = x_i`
//!    to get a 1-variable polynomial in `x_2`.
//! 4. For each factor-base point `P_j`, evaluate the 1-variable polynomial at `x_j`.
//!    If it evaluates to 0, then `S_3(Q.x, P_i.x, P_j.x) = 0`, meaning there exist
//!    y-values such that `Q' + P_i + P_j = ∞` for some point `Q'` with x-coordinate
//!    `Q.x` (either `Q` or `-Q`).
//! 5. Verify via the group law: check whether `P_i + P_j = Q` (or `P_i + P_j = -Q`).
//!    Accept only the case `P_i + P_j = Q`.
//! 6. Return the first valid decomposition, or `None` if none exists.
//!
//! # Semaev convention
//!
//! `S_3(x_0, x_1, x_2) = 0` iff there exist y_i such that `P_0 + P_1 + P_2 = ∞`.
//! This is an existential condition on y-coordinates — the polynomial depends only on
//! x-coordinates. A vanishing triple `(Q.x, P_i.x, P_j.x)` means either `Q + P_i + P_j = ∞`
//! (i.e., `P_i + P_j = -Q`) or `-Q + P_i + P_j = ∞` (i.e., `P_i + P_j = Q`). The group-law
//! check disambiguates.
//!
//! # Green path vs msolve sidecar
//!
//! The green path uses native enumeration over the factor base (O(|FB|²) for m = 2).
//! Native enumeration does not scale past toy `m`.
//! An optional msolve cross-check is `#[ignore]`-gated (not implemented here — the
//! green path is the deliverable).
//!
//! # Scale limits
//!
//! Native enumeration is O(|FB|^m) — acceptable for toy `m = 2` and `|FB| = 6`.
//! Crypto-scale decomposition requires Gröbner-basis solvers (msolve, F4/F5). Over `E(F_p)` the
//! mechanism is demonstrated, not the asymptotic advantage.

use shared_field::{Fp, FpNaive4 as FpNaive};

use crate::curve::{AffinePoint, JacobianPoint};
use crate::index_calculus::strategy::{FbPoint, IndexCalcStrategy};
use crate::semaev::semaev_poly;

/// Attempt to decompose `q` as a sum of `strategy.m` factor-base points.
///
/// Returns `Some(fb_points)` where `fb_points` is a list of `m` factor-base points
/// (with repetition allowed) whose sum equals `q` via the frozen group law. Returns
/// `None` if no such decomposition exists over the current factor base.
///
/// Green path: native enumeration via Semaev polynomial root-finding (no live oracle).
/// The msolve cross-check is `#[ignore]`-gated.
///
/// # Panics (debug)
///
/// Panics if `q` is the point at infinity (decomposition is undefined for ∞).
pub fn decompose(q: AffinePoint<FpNaive>, strategy: &IndexCalcStrategy) -> Option<Vec<FbPoint>> {
    debug_assert!(!q.is_infinity(), "decompose: q must be a finite point");

    // Only m = 2 is ratified; the algorithm is specialised for this arity.
    // SCALE: generalising to m > 2 requires an m-nested loop or a Gröbner solver.
    assert_eq!(strategy.m, 2, "decompose: only m = 2 is implemented (native enumeration does not scale)");

    let curve = &strategy.curve;
    let p_uint = &curve.p;
    let p_u64 = p_uint.as_words()[0];
    let a_u64 = curve.a.as_words()[0];
    let b_u64 = curve.b.as_words()[0];

    // Extract Q's x-coordinate as u64 (toy-scale: p = 47 fits in u64).
    let q_x_u64 = match &q {
        AffinePoint::Finite { x, .. } => x.to_uint().as_words()[0],
        AffinePoint::Infinity => unreachable!("guarded by debug_assert above"),
    };

    // Step 1: Compute S_3(x_0, x_1, x_2) — the 3-variable Semaev polynomial.
    // S_3 has num_vars = 3; variables are indexed 0, 1, 2.
    let s3 = semaev_poly(3, a_u64, b_u64, p_u64)
        .expect("semaev_poly(3) should not fail for valid curve parameters");

    // Step 2: Specialise at Q's x-coordinate: fix variable 0 = Q.x.
    // partial_eval takes &[Option<u64>]: Some(val) fixes the variable, None leaves it free.
    // After fixing x_0 = Q.x, the result is a 2-variable polynomial in (x_1, x_2).
    let s3_at_q = s3
        .partial_eval(&[Some(q_x_u64), None, None])
        .expect("partial_eval: arity 3, fixing variable 0");
    // s3_at_q has num_vars = 2 (the two free variables, renumbered 0 and 1).

    // Steps 3–5: Enumerate factor-base pairs (P_i, P_j).
    // For each P_i, substitute x_0 = P_i.x into s3_at_q → 1-variable polynomial in x_1.
    // For each P_j, evaluate at x_1 = P_j.x → scalar. If 0, check group law.
    for fb_i in &strategy.factor_base {
        let x_i_u64 = fb_point_x_u64(fb_i);

        // Step 3: Substitute x_0 = P_i.x into the 2-variable polynomial.
        // s3_at_q has variables (x_1, x_2) renumbered as (0, 1).
        // Fix variable 0 = P_i.x → 1-variable polynomial in x_2 (renumbered to 0).
        let s3_at_q_pi = s3_at_q
            .partial_eval(&[Some(x_i_u64), None])
            .expect("partial_eval: arity 2, fixing variable 0");
        // s3_at_q_pi has num_vars = 1.

        for fb_j in &strategy.factor_base {
            let x_j_u64 = fb_point_x_u64(fb_j);

            // Step 4: Evaluate the 1-variable polynomial at P_j.x.
            let val = s3_at_q_pi
                .eval(&[x_j_u64])
                .expect("eval: arity 1");

            if val != 0 {
                // S_3(Q.x, P_i.x, P_j.x) ≠ 0 — not a Semaev triple.
                continue;
            }

            // S_3(Q.x, P_i.x, P_j.x) = 0: there exist y-values such that some point
            // with x = Q.x, P_i, P_j sum to ∞. Check whether P_i + P_j = Q via the
            // group law (the Semaev condition is existential over y-values of Q).
            //
            // Step 5: Verify P_i + P_j = Q using the frozen group law.
            let sum = add_affine(curve, &fb_i.point, &fb_j.point, p_uint);
            if sum == q {
                return Some(vec![fb_i.clone(), fb_j.clone()]);
            }
            // If P_i + P_j = -Q (the other y-branch), this is not a valid decomposition
            // for Q — skip. The Semaev polynomial is symmetric in all three arguments,
            // so the (Q, P_i, P_j) triple may correspond to -Q + P_i + P_j = ∞.
        }
    }

    None
}

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Extract the x-coordinate of a factor-base point as `u64` (toy-scale).
fn fb_point_x_u64(fb: &FbPoint) -> u64 {
    match &fb.point {
        AffinePoint::Finite { x, .. } => x.to_uint().as_words()[0],
        AffinePoint::Infinity => panic!("factor-base point is the point at infinity"),
    }
}

/// Add two affine points using the frozen `Curve` group law.
///
/// Converts to Jacobian, adds, converts back to affine.
fn add_affine(
    curve: &crate::curve::Curve,
    p1: &AffinePoint<FpNaive>,
    p2: &AffinePoint<FpNaive>,
    p: &crypto_bigint::Uint<4>,
) -> AffinePoint<FpNaive> {
    let j1 = JacobianPoint::from_affine(p1, p);
    let j2 = JacobianPoint::from_affine(p2, p);
    curve.add_jacobian(&j1, &j2).to_affine(p)
}

// ─── unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index_calculus::strategy::IndexCalcStrategy;

    /// Smoke test: `decompose` does not panic on a factor-base point.
    ///
    /// A factor-base point P_0 + P_1 = some point Q; decompose(Q) should return
    /// Some([P_0, P_1]) or Some([P_1, P_0]) (or another valid pair).
    #[test]
    fn decompose_smoke_no_panic() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let curve = &strategy.curve;
        let p = &curve.p;

        // Compute P_0 + P_1 directly.
        let p0 = &strategy.factor_base[0].point;
        let p1 = &strategy.factor_base[1].point;
        let j0 = JacobianPoint::from_affine(p0, p);
        let j1 = JacobianPoint::from_affine(p1, p);
        let q = curve.add_jacobian(&j0, &j1).to_affine(p);

        if q.is_infinity() {
            // P_0 + P_1 = ∞ — skip (decomposition undefined for ∞).
            return;
        }

        // decompose should not panic.
        let _result = decompose(q, &strategy);
    }

    /// Decomposition correctness: if `decompose(Q)` returns `Some(pts)`, then
    /// `Σ pts = Q` via the frozen group law.
    #[test]
    fn decompose_result_sums_to_q() {
        let strategy = IndexCalcStrategy::toy().expect("toy strategy should build");
        let curve = &strategy.curve;
        let p = &curve.p;

        // Try all pairs of factor-base points to find a decomposable Q.
        for fb_i in &strategy.factor_base {
            for fb_j in &strategy.factor_base {
                let ji = JacobianPoint::from_affine(&fb_i.point, p);
                let jj = JacobianPoint::from_affine(&fb_j.point, p);
                let q = curve.add_jacobian(&ji, &jj).to_affine(p);

                if q.is_infinity() {
                    continue;
                }

                if let Some(decomp) = decompose(q.clone(), &strategy) {
                    // Verify the decomposition sums to Q.
                    assert_eq!(decomp.len(), 2, "m = 2: decomposition should have 2 points");
                    let jd0 = JacobianPoint::from_affine(&decomp[0].point, p);
                    let jd1 = JacobianPoint::from_affine(&decomp[1].point, p);
                    let sum = curve.add_jacobian(&jd0, &jd1).to_affine(p);
                    assert_eq!(
                        sum, q,
                        "decompose returned a decomposition that does not sum to Q"
                    );
                    return; // Found and verified one decomposition — test passes.
                }
            }
        }
        // If no decomposable Q was found among all pairs, the test is vacuously passing
        // (the factor base is too small or the curve has no decomposable points).
        // This should not happen for the toy fixture.
        panic!("no decomposable Q found among all factor-base pairs — unexpected for toy fixture");
    }
}
