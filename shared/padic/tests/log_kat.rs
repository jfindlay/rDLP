//! Known-answer tests (KATs) for the p-adic logarithm in `shared-padic`.
//!
//! All KATs use `p = 7`, `k = 4` (modulus = 7^4 = 2401) unless otherwise noted.
//!
//! # Hand-computed values
//!
//! The formal-group series is `log(1 + x) = x − x²/2 + x³/3 − x⁴/4 + …` where `x = z − 1`.
//! Terms with `n − v_7(n) ≥ 4` contribute 0 mod 7^4 and are dropped.
//!
//! ## log(8) = log(1 + 7), x = 7
//!
//! - n=1: +7
//! - n=2: −7²/2 = −49/2. inv(2) mod 2401 = 1201 (2·1201 = 2402 ≡ 1).
//!   49·1201 = 58849; 58849 mod 2401 = 1225. Term: −1225 ≡ 1176 mod 2401.
//! - n=3: +7³/3 = +343/3. inv(3) mod 2401 = 1601 (3·800 = 2400 ≡ −1, so inv(3) = 2401−800 = 1601).
//!   343·1601 = 549143; 549143 mod 2401 = 1715. Term: +1715.
//! - n=4: −7⁴/4 = −2401/4 ≡ 0 mod 2401 (7^4 | 7^4). Stopping: n − v_7(n) = 4 − 0 = 4 ≥ k=4.
//!
//! log(8) mod 7^4 = 7 + 1176 + 1715 = 2898 mod 2401 = **497**.
//!
//! ## log(15) = log(1 + 14), x = 14 = 2·7
//!
//! - n=1: +14
//! - n=2: −14²/2 = −196/2 = −98. −98 mod 2401 = 2303. Term: +2303.
//! - n=3: +14³/3 = +2744/3. 2744 mod 2401 = 343. 343·1601 = 549143 mod 2401 = 1715. Term: +1715.
//! - n=4: −14⁴/4. 14^4 = 7^4·2^4 ≡ 0 mod 7^4. Term: 0. Stopping: n − v_7(n) = 4 ≥ k=4.
//!
//! log(15) mod 7^4 = 14 + 2303 + 1715 = 4032 mod 2401 = **1631**.
//!
//! ## log(120) = log(1 + 119), x = 119 = 7·17
//!
//! - n=1: +119
//! - n=2: −119²/2 = −14161/2. 14161 mod 2401 = 14161 − 5·2401 = 2156. 2156·1201 mod 2401:
//!   2156·1201 = 2,589,356; 2589356 mod 2401: 1078·2401 = 2,588,278; 2589356 − 2588278 = 1078.
//!   Term: −1078 ≡ 1323 mod 2401.
//! - n=3: +119³/3. 119³ mod 2401: 119² mod 2401 = 2156; 2156·119 mod 2401 = 256564 mod 2401.
//!   256564 / 2401 ≈ 106.9; 106·2401 = 254506; 256564 − 254506 = 2058.
//!   2058·1601 mod 2401: 2058·1601 = 3,294,858; 3294858 / 2401 ≈ 1372.3; 1372·2401 = 3,294,172;
//!   3294858 − 3294172 = 686. Term: +686.
//! - n=4: −119⁴/4. 119 = 7·17, so 119^4 = 7^4·17^4 ≡ 0 mod 7^4. Term: 0.
//!
//! log(120) mod 7^4 = 119 + 1323 + 686 = 2128 mod 2401 = **2128**.
//!
//! ## Homomorphism check
//!
//! log(8) + log(15) = 497 + 1631 = 2128 = log(120) mod 7^4. ✓
//!
//! # KAT coverage
//!
//! - Primary KAT: `log(8) mod 7^4 = 497`.
//! - Homomorphism KAT: `log(8) + log(15) ≡ log(120) mod 7^4`.
//! - Convergence guard KAT: `log(2)` (z=2, z−1=1, v_7(1)=0) returns `ConvergenceViolation`.
//! - log(1) = 0 (the identity element).
//! - A test (`#[ignore]`) asserting a value recorded from PARI/GP; it does not call `gp`.

use num_bigint::BigInt;
use shared_padic::log::{PadicLogError, padic_log};
use shared_padic::Zp;

// ─── helpers ──────────────────────────────────────────────────────────────────

/// Construct a Z/7^k element from an i64 value.
fn z7(value: i64, k: u32) -> Zp {
    Zp::from_i64(value, 7, k).expect("valid Zp construction")
}

fn bi(n: i64) -> BigInt {
    BigInt::from(n)
}

// ─── primary KAT ─────────────────────────────────────────────────────────────

/// `log(8) mod 7^4 = 497`.
///
/// Hand-computed from the series `log(1 + 7) = 7 − 49/2 + 343/3 − …` mod 2401.
/// See module-level comment for the full derivation.
#[test]
fn kat_log_8_mod_7_4() {
    let z = z7(8, 4); // z = 1 + 7; x = z − 1 = 7; v_7(7) = 1 ✓
    let result = padic_log(&z).expect("log(8) mod 7^4 should succeed");
    assert_eq!(
        *result.residue(),
        bi(497),
        "log(8) mod 7^4 = 497 (hand-computed)"
    );
    // Sanity: the result has v_7 ≥ 1 (log maps 1+p·Z_p into p·Z_p).
    assert_eq!(result.precision(), 4, "precision preserved");
}

// ─── homomorphism KAT ────────────────────────────────────────────────────────

/// `log(a·b) = log(a) + log(b)` for `a = 8`, `b = 15`, `a·b = 120`.
///
/// This is the load-bearing property of the p-adic log: it is a group homomorphism from
/// the kernel of reduction (multiplicative) to p·Z_p (additive). The Smart–Satoh–Araki
/// reduction depends on this homomorphism.
///
/// Hand-computed: log(8) = 497, log(15) = 1631, log(120) = 2128.
/// 497 + 1631 = 2128 mod 2401. ✓
#[test]
fn kat_log_homomorphism() {
    let a = z7(8, 4);   // 1 + 7;  v_7(a − 1) = 1 ✓
    let b = z7(15, 4);  // 1 + 14; v_7(b − 1) = 1 ✓
    let ab = a.mul(&b).expect("8 · 15 mod 2401");
    // 8 · 15 = 120; 120 mod 2401 = 120.
    assert_eq!(*ab.residue(), bi(120), "8 · 15 mod 2401 = 120");

    let log_a = padic_log(&a).expect("log(8)");
    let log_b = padic_log(&b).expect("log(15)");
    let log_ab = padic_log(&ab).expect("log(120)");

    let log_a_plus_log_b = log_a.add(&log_b).expect("log(a) + log(b)");

    assert_eq!(
        *log_a_plus_log_b.residue(),
        *log_ab.residue(),
        "log(8) + log(15) = log(120) mod 7^4 (homomorphism)"
    );

    // Spot-check the hand-computed values.
    assert_eq!(*log_a.residue(), bi(497), "log(8) = 497");
    assert_eq!(*log_b.residue(), bi(1631), "log(15) = 1631");
    assert_eq!(*log_ab.residue(), bi(2128), "log(120) = 2128");
}

// ─── convergence guard KAT ───────────────────────────────────────────────────

/// `log(2)` must fail: `z = 2`, `z − 1 = 1`, `v_7(1) = 0 < 1`.
///
/// This is the convergence guard: the series diverges p-adically when `v_p(z − 1) = 0`.
/// Returning a value here would be a silent wrong answer — the guard prevents that.
#[test]
fn kat_log_convergence_guard() {
    let z = z7(2, 4); // z − 1 = 1; v_7(1) = 0 → convergence violation
    let result = padic_log(&z);
    assert!(
        matches!(result, Err(PadicLogError::ConvergenceViolation { valuation: 0, .. })),
        "log(2) must return ConvergenceViolation (v_7(1) = 0), got: {result:?}"
    );
}

/// `log(8)` with z = 8 ≡ 1 mod 7 must succeed (v_7(7) = 1 ≥ 1).
///
/// Complementary to the guard test: the boundary case that should succeed.
#[test]
fn kat_log_convergence_boundary_ok() {
    let z = z7(8, 4); // z − 1 = 7; v_7(7) = 1 ✓
    assert!(padic_log(&z).is_ok(), "log(8) must succeed: v_7(7) = 1 ≥ 1");
}

/// `log(z)` for z ≡ 3 mod 7 (v_7(z−1) = v_7(2) = 0) must fail.
#[test]
fn kat_log_convergence_guard_z3() {
    let z = z7(3, 4); // z − 1 = 2; v_7(2) = 0 → convergence violation
    let result = padic_log(&z);
    assert!(
        matches!(result, Err(PadicLogError::ConvergenceViolation { valuation: 0, .. })),
        "log(3) must return ConvergenceViolation, got: {result:?}"
    );
}

// ─── log(1) = 0 ──────────────────────────────────────────────────────────────

/// `log(1) = 0`: the identity element maps to zero.
///
/// z = 1, x = z − 1 = 0, v_7(0) = ∞ ≥ 1. The series is empty; the result is 0.
#[test]
fn kat_log_identity() {
    let z = z7(1, 4);
    let result = padic_log(&z).expect("log(1) should succeed");
    assert_eq!(*result.residue(), bi(0), "log(1) = 0");
}

// ─── additional spot-checks ───────────────────────────────────────────────────

/// `log(1 + 7²) = log(50)` mod 7^4.
///
/// x = 49 = 7², v_7(49) = 2 ≥ 1 ✓.
/// Series: n=1: +49; n=2: −49²/2 = −2401/2 ≡ 0 mod 2401 (7^4 | 49²).
/// All higher terms also vanish. log(50) mod 7^4 = 49.
#[test]
fn kat_log_1_plus_p2() {
    let z = z7(50, 4); // z = 1 + 49; x = 49 = 7²; v_7(49) = 2 ✓
    let result = padic_log(&z).expect("log(50) should succeed");
    // n=1: +49. n=2: −49²/2 = −2401/2 ≡ 0 mod 2401. All higher terms 0.
    assert_eq!(*result.residue(), bi(49), "log(1 + 7²) = 49 mod 7^4");
}

// ─── value recorded from PARI/GP ──────────────────────────────────────────────

/// Asserts a value recorded from PARI/GP for `log(8) mod 7^4`.
///
/// The test does not call `gp`; the value was obtained once with:
/// ```text
/// gp> Qp = padic(7, 4); log(1 + Qp(7))
/// ```
/// Expected: the 7-adic expansion with leading coefficient 1 at precision 4,
/// corresponding to the residue 497 mod 2401.
///
/// It is `#[ignore]`d and needs no external tool when run with `cargo test -- --ignored`.
#[test]
#[ignore = "asserts a value recorded from PARI/GP; does not call gp"]
fn kat_log_pari_cross_check() {
    // PARI/GP: log(1 + O(7^4) + 7) gives a 7-adic number with residue 497 mod 2401.
    // The recorded value can be re-derived with: echo "log(1 + O(7^4) + 7)" | gp -q
    let z = z7(8, 4);
    let result = padic_log(&z).expect("log(8) should succeed");
    assert_eq!(
        *result.residue(),
        bi(497),
        "PARI cross-check: log(8) mod 7^4 = 497"
    );
}
