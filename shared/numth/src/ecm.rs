//! Lenstra's Elliptic Curve Method (ECM) for integer factorization.
//!
//! ECM (Lenstra 1987) factors a composite ``N`` by choosing a random elliptic
//! curve ``E`` over ``Z/NZ`` and computing a ``B1``-smooth scalar multiple of a
//! starting point ``P``.  The key insight: if ``p | N`` and the group order
//! ``|E(F_p)|`` is ``B1``-smooth, then computing ``k * P`` (for ``k`` the
//! product of prime powers ≤ ``B1``) drives the point to the identity mod ``p``
//! while it remains non-identity mod ``N``.  A failed modular inversion in the
//! group law then reveals ``gcd(denominator, N) = p``.
//!
//! # Curve form
//!
//! ECM uses **Montgomery-form curves** (Suyama parameterization):
//! ``B * v² = u³ + A * u² + u`` over ``Z/NZ``.  Montgomery projective
//! coordinates ``(X:Z)`` omit the ``Y`` coordinate; only the ``x``-coordinate
//! of scalar multiples is needed.  The Montgomery ladder computes ``k * P``
//! using only differential addition and doubling, which are efficient and
//! numerically stable.
//!
//! # Note on curve form choice
//!
//! ECM uses Montgomery-form curves (Suyama parameterization), which are
//! categorically different from the short Weierstrass curves in ``rho::curve``.
//! Reusing ``rho``'s curve machinery would add complexity without benefit;
//! the Montgomery group law is implemented directly here.
//!
//! # Stage 1 and Stage 2
//!
//! **Stage 1**: compute ``k * P`` where ``k = ∏ { q^⌊log_q(B1)⌋ : q prime, q ≤ B1 }``, by
//! applying each prime power ``q`` to the running point with its own ladder call (the product
//! ``k ≈ e^B1`` is never formed).  A factor is found if ``gcd(Z-coordinate, N)`` is a proper
//! divisor after any of those calls.
//!
//! **Stage 2**: given the stage-1 result ``Q = k * P``, check primes ``q`` in
//! ``(B1, B2]``.  If ``q * Q`` is the identity mod ``p`` for some ``p | N``,
//! then ``gcd(∏ x(j*Q) - x(r*Q), N) > 1``.  Uses a baby-step giant-step
//! approach with step size 2.
//!
//! # Deferred
//!
//! - Brent–Suyama parameterization (improved sigma formula) — TODO.
//! - FFT continuation for stage 2 — TODO.
//! - Stage 2 step-size optimization — TODO.

use crypto_bigint::{NonZero, Odd, Uint};

use crate::prime::{is_prime, MILLER_RABIN_PROOF_BOUND};
use crate::smooth::factor_base_up_to;

// ── Modular arithmetic helpers ────────────────────────────────────────────────
//
// ECM operates mod N (a composite), so we cannot use the Fp<L> trait (which
// assumes the modulus is prime and uses Fermat's little theorem for inversion).
// All arithmetic is implemented directly on Uint<4>.

/// Compute ``(a * b) mod n`` for ``a, b < n``.
///
/// Widens to ``Uint<8>`` for the full product, then reduces via ``rem``.
/// This mirrors the approach in ``FpNaive::mul``.
#[inline]
fn mul_mod(a: &Uint<4>, b: &Uint<4>, n: &Uint<4>) -> Uint<4> {
    let (lo, hi) = a.widening_mul(b);
    // Assemble the 512-bit product: value = hi * 2^256 + lo.
    let wide = Uint::<8>::from((lo, hi));
    // Embed n into Uint<8> as a low-half value for rem.
    let n_wide = Uint::<8>::from((*n, Uint::<4>::ZERO));
    let nz = NonZero::new(n_wide).expect("modulus must be non-zero");
    let rem: Uint<8> = wide.rem(&nz);
    // rem < n < 2^256, so the high 4 limbs are zero; extract the low 4.
    let words = rem.as_words();
    let mut lo_words = [0u64; 4];
    lo_words.copy_from_slice(&words[..4]);
    Uint::<4>::from_words(lo_words)
}

/// Compute ``(a + b) mod n`` for ``a, b < n``.
#[inline]
fn add_mod(a: &Uint<4>, b: &Uint<4>, n: &Uint<4>) -> Uint<4> {
    a.add_mod(b, &NonZero::new(*n).expect("modulus must be non-zero"))
}

/// Compute ``(a - b) mod n`` for ``a, b < n``.
#[inline]
fn sub_mod(a: &Uint<4>, b: &Uint<4>, n: &Uint<4>) -> Uint<4> {
    a.sub_mod(b, &NonZero::new(*n).expect("modulus must be non-zero"))
}

/// Binary GCD (Stein's algorithm) for ``Uint<4>``.
///
/// Returns ``gcd(a, b)``.  Handles the case where either input is zero.
fn gcd(mut a: Uint<4>, mut b: Uint<4>) -> Uint<4> {
    if a == Uint::<4>::ZERO {
        return b;
    }
    if b == Uint::<4>::ZERO {
        return a;
    }

    // Factor out common powers of 2.
    let a_tz = a.trailing_zeros();
    let b_tz = b.trailing_zeros();
    let common_twos = a_tz.min(b_tz);

    a >>= a_tz;
    b >>= b_tz;

    loop {
        // Both a and b are now odd.
        // Ensure a <= b.
        if a > b {
            let tmp = a;
            a = b;
            b = tmp;
        }
        // b = b - a (b - a >= 0 since a <= b).
        b = b.wrapping_sub(&a);
        if b == Uint::<4>::ZERO {
            break;
        }
        // Remove trailing zeros from b.
        b >>= b.trailing_zeros();
    }

    // Restore common factor of 2.
    a <<= common_twos;
    a
}

/// Attempt to compute the modular inverse of ``a`` mod ``n`` (``n`` must be odd).
///
/// Returns ``Ok(a^{-1} mod n)`` if ``gcd(a, n) == 1``, or ``Err(g)`` where
/// ``g = gcd(a, n)`` if ``g > 1``.  The ``Err`` case signals a factor of ``n``
/// has been found.
///
/// Uses ``crypto_bigint``'s ``inv_odd_mod``, which implements the constant-time
/// binary extended GCD.  When the inverse does not exist, we fall back to a
/// binary GCD to recover the actual factor.
fn try_inv(a: &Uint<4>, n: &Uint<4>) -> Result<Uint<4>, Uint<4>> {
    if *a == Uint::<4>::ZERO {
        // gcd(0, n) = n — degenerate, return n as the "factor" (caller discards).
        return Err(*n);
    }
    let odd_n = Odd::new(*n).expect("ECM modulus must be odd");
    if let Some(inv) = a.invert_odd_mod(&odd_n).into_option() {
        Ok(inv)
    } else {
        // No inverse: gcd(a, n) > 1.  Compute the actual gcd.
        let g = gcd(*a, *n);
        Err(g)
    }
}

// ── Montgomery projective point ───────────────────────────────────────────────

/// A point on a Montgomery curve in projective ``(X:Z)`` coordinates.
///
/// The affine ``x``-coordinate is ``X * Z^{-1} mod N``.  The ``Y`` coordinate
/// is not tracked; ECM only needs ``x``-coordinates for the scalar multiple.
/// The point at infinity is represented as ``Z = 0``.
#[derive(Clone, Debug)]
struct MontPoint {
    x: Uint<4>,
    z: Uint<4>,
}

impl MontPoint {
    /// The point at infinity: ``(1:0)``.
    fn infinity() -> Self {
        MontPoint { x: Uint::<4>::ONE, z: Uint::<4>::ZERO }
    }

    /// Return ``true`` if this point is the identity (``Z == 0``).
    fn is_infinity(&self) -> bool {
        self.z == Uint::<4>::ZERO
    }
}

// ── Montgomery curve differential addition and doubling ──────────────────────
//
// Reference: Bernstein–Birkner–Joye–Lange–Peters, "Twisted Edwards Curves",
// and the standard Montgomery ladder literature.
//
// The curve parameter is stored as ``A24 = (A + 2) / 4 mod N``.
// This avoids repeated division by 4 in the inner loop.

/// Differential addition: given ``P``, ``Q``, and ``P - Q``, compute ``P + Q``.
///
/// All points are in projective ``(X:Z)`` coordinates.  The formula requires
/// ``P - Q`` (the difference) as an additional input, which is always available
/// in the Montgomery ladder.
///
/// Formula (from Montgomery 1987, Algorithm 3):
/// ```text
/// U = (Xp - Zp)(Xq + Zq) + (Xp + Zp)(Xq - Zq)
/// V = (Xp - Zp)(Xq + Zq) - (Xp + Zp)(Xq - Zq)
/// X_{P+Q} = Z_{P-Q} * U^2
/// Z_{P+Q} = X_{P-Q} * V^2
/// ```
fn diff_add(p: &MontPoint, q: &MontPoint, diff: &MontPoint, n: &Uint<4>) -> MontPoint {
    let u = mul_mod(
        &add_mod(&p.x, &p.z, n),
        &sub_mod(&q.x, &q.z, n),
        n,
    );
    let v = mul_mod(
        &sub_mod(&p.x, &p.z, n),
        &add_mod(&q.x, &q.z, n),
        n,
    );
    let add_uv = add_mod(&u, &v, n);
    let sub_uv = sub_mod(&u, &v, n);
    let x = mul_mod(&diff.z, &mul_mod(&add_uv, &add_uv, n), n);
    let z = mul_mod(&diff.x, &mul_mod(&sub_uv, &sub_uv, n), n);
    MontPoint { x, z }
}

/// Point doubling: given ``P`` and the curve parameter ``A24 = (A+2)/4``, compute ``2P``.
///
/// Formula (from Montgomery 1987):
/// ```text
/// T1 = (X + Z)^2
/// T2 = (X - Z)^2
/// X_{2P} = T1 * T2
/// T3 = T1 - T2
/// Z_{2P} = T3 * (T2 + A24 * T3)
/// ```
fn point_double(p: &MontPoint, a24: &Uint<4>, n: &Uint<4>) -> MontPoint {
    let xpz = add_mod(&p.x, &p.z, n);
    let xmz = sub_mod(&p.x, &p.z, n);
    let t1 = mul_mod(&xpz, &xpz, n);
    let t2 = mul_mod(&xmz, &xmz, n);
    let x_new = mul_mod(&t1, &t2, n);
    let t3 = sub_mod(&t1, &t2, n);
    let z_new = mul_mod(&t3, &add_mod(&t2, &mul_mod(a24, &t3, n), n), n);
    MontPoint { x: x_new, z: z_new }
}

// ── Montgomery ladder ─────────────────────────────────────────────────────────

/// Compute ``k * P`` on a Montgomery curve using the Montgomery ladder.
///
/// Processes bits of ``k`` from MSB to LSB.  Maintains the invariant that
/// ``R1 - R0 = P`` throughout.  Returns the resulting point ``k * P``.
///
/// The ladder never requires a modular inversion; factors are detected by the
/// caller inspecting the ``Z`` coordinate of the result.
fn montgomery_ladder(p: &MontPoint, k: &Uint<4>, a24: &Uint<4>, n: &Uint<4>) -> MontPoint {
    if *k == Uint::<4>::ZERO {
        return MontPoint::infinity();
    }

    // Find the most significant bit of k.
    let bit_len = k.bits();
    if bit_len == 0 {
        return MontPoint::infinity();
    }

    // Initialise: R0 = P, R1 = 2P.  Invariant: R1 - R0 = P.
    let mut r0 = p.clone();
    let mut r1 = point_double(p, a24, n);

    // Process bits from MSB-1 down to 0 (the MSB itself is always 1, so we
    // start the ladder with R0=P, R1=2P and process the remaining bits).
    for i in (0..bit_len - 1).rev() {
        if bool::from(k.bit(i)) {
            // Bit is 1: (R0, R1) ← (R0 + R1, 2*R1).
            r0 = diff_add(&r0, &r1, p, n);
            r1 = point_double(&r1, a24, n);
        } else {
            // Bit is 0: (R0, R1) ← (2*R0, R0 + R1).
            r1 = diff_add(&r0, &r1, p, n);
            r0 = point_double(&r0, a24, n);
        }
    }

    r0
}

// ── Suyama parameterization ───────────────────────────────────────────────────

/// Result of the Suyama parameterization: a Montgomery curve and starting point.
enum SuyamaResult {
    /// Curve and starting point generated successfully.
    Ok { a24: Uint<4>, p: MontPoint },
    /// A non-trivial factor of N was found during parameterization.
    Factor(Uint<4>),
}

/// Generate a Montgomery curve and starting point from the Suyama parameter ``sigma``.
///
/// The Suyama parameterization (Suyama 1985) generates a curve and a starting
/// point with a known rational point, which ensures the starting point is on
/// the curve.  The formulas are:
///
/// ```text
/// u = sigma^2 - 5  (mod N)
/// v = 4 * sigma    (mod N)
/// A24 = (v - u)^3 * (3*u + v) / (4 * u^3 * v) - 2   (mod N)
/// ```
///
/// where ``A24 = (A + 2) / 4``.  The starting point is ``(u^3 : v^3)``.
///
/// If the denominator ``4 * u^3 * v`` shares a factor with ``N``, that factor
/// is returned directly.
fn suyama_param(sigma: u64, n: &Uint<4>) -> SuyamaResult {
    let sigma_u = Uint::<4>::from(sigma);
    let sigma_mod = {
        let nz = NonZero::new(*n).expect("n must be non-zero");
        sigma_u.rem(&nz)
    };

    // u = sigma^2 - 5 (mod N).
    let sigma_sq = mul_mod(&sigma_mod, &sigma_mod, n);
    let five = Uint::<4>::from(5u64);
    // sigma^2 - 5: use sub_mod, but sigma^2 may be < 5, so use wrapping sub with mod.
    let u = if sigma_sq >= five {
        sigma_sq.wrapping_sub(&five)
    } else {
        // sigma^2 < 5: result is sigma^2 - 5 + N.
        add_mod(&sigma_sq, &n.wrapping_sub(&five), n)
    };

    // v = 4 * sigma (mod N).
    let four = Uint::<4>::from(4u64);
    let v = mul_mod(&four, &sigma_mod, n);

    // u^3 mod N.
    let u2 = mul_mod(&u, &u, n);
    let u3 = mul_mod(&u2, &u, n);

    // v^3 mod N.
    let v2 = mul_mod(&v, &v, n);
    let v3 = mul_mod(&v2, &v, n);

    // Denominator: 4 * u^3 * v.
    let denom = mul_mod(&four, &mul_mod(&u3, &v, n), n);

    // Attempt to invert the denominator.  If gcd(denom, N) > 1, we found a factor.
    let denom_inv = match try_inv(&denom, n) {
        Ok(inv) => inv,
        Err(g) => {
            if g > Uint::<4>::ONE && g < *n {
                return SuyamaResult::Factor(g);
            }
            // g == 1 (shouldn't happen since try_inv returned Err) or g == N (degenerate).
            // Fall through: treat as inconclusive by returning a dummy curve.
            // This branch is unreachable in practice.
            return SuyamaResult::Factor(*n);
        }
    };

    // Numerator: (v - u)^3 * (3*u + v).
    let vmu = sub_mod(&v, &u, n);
    let vmu2 = mul_mod(&vmu, &vmu, n);
    let vmu3 = mul_mod(&vmu2, &vmu, n);
    let three = Uint::<4>::from(3u64);
    let three_u_plus_v = add_mod(&mul_mod(&three, &u, n), &v, n);
    let numer = mul_mod(&vmu3, &three_u_plus_v, n);

    // A24 = numer * denom_inv - 2 (mod N).
    let a24_plus_2 = mul_mod(&numer, &denom_inv, n);
    let two = Uint::<4>::from(2u64);
    let a24 = sub_mod(&a24_plus_2, &two, n);

    // Starting point: (u^3 : v^3).
    let p = MontPoint { x: u3, z: v3 };

    SuyamaResult::Ok { a24, p }
}

// ── Stage 1 ───────────────────────────────────────────────────────────────────

/// Return the largest power of the prime ``p`` that is ≤ ``b1``: ``p^⌊log_p(B1)⌋``.
///
/// Requires ``2 ≤ p ≤ b1``.  The power always fits in a ``u64`` because it is ≤ ``b1``.
fn largest_prime_power_le(p: u64, b1: u64) -> u64 {
    let mut q = p;
    // `q <= b1 / p` is `q * p <= b1` without the overflow.
    while q <= b1 / p {
        q *= p;
    }
    q
}

/// Run stage 1 of ECM on a single curve.
///
/// Applies ``q = p^⌊log_p(B1)⌋`` to the point, one prime ``p ≤ B1`` at a time, each with its own
/// ladder call.  The composite scalar ``∏ q = lcm(1, …, B1)`` grows like ``e^B1`` and never fits
/// a fixed-width integer for large ``B1``; applying the prime powers separately yields the same
/// multiple ``(∏ q) * P`` at every ``B1`` without forming that product.  If ``|E(F_r)|`` is
/// ``B1``-power-smooth for a prime ``r | N``, the point becomes the identity mod ``r``.
///
/// After each ladder call ``gcd(Z, N)`` is tested.  Once ``Z ≡ 0 (mod r)`` the point stays at
/// the identity mod ``r`` (every doubling and differential addition of a point with ``Z ≡ 0``
/// has ``Z ≡ 0``), so the test could be deferred to the end; testing per prime instead returns
/// the factor as soon as it appears, before a second prime of ``N`` also reaches the identity
/// and makes the gcd with ``N`` improper.
///
/// Returns ``Ok(point)`` if stage 1 completed without finding a proper factor, or
/// ``Err(factor)`` for a non-trivial factor of ``N``.  If ``Z ≡ 0 (mod N)`` the identity is
/// returned immediately: every later multiple is the identity as well.
fn run_stage1(
    p: &MontPoint,
    a24: &Uint<4>,
    b1: u64,
    n: &Uint<4>,
) -> Result<MontPoint, Uint<4>> {
    let mut q = p.clone();
    for prime in factor_base_up_to(b1) {
        let power = Uint::<4>::from(largest_prime_power_le(prime, b1));
        q = montgomery_ladder(&q, &power, a24, n);

        let g = gcd(q.z, *n);
        if g == *n {
            return Ok(q);
        }
        if g > Uint::<4>::ONE {
            return Err(g);
        }
    }

    Ok(q)
}

// ── Stage 2 ───────────────────────────────────────────────────────────────────

/// Run stage 2 of ECM (standard continuation).
///
/// Given the stage-1 result ``Q = k * P``, checks primes ``q`` in ``(B1, B2]``
/// by accumulating the product of Z-coordinates of ``j * Q`` for each prime
/// ``j`` in ``(B1, B2]``, then taking a GCD with ``N``.
///
/// If ``Z_{jQ} ≡ 0 (mod p)`` for some prime ``p | N``, then ``j * Q`` is the
/// identity in ``E(F_p)``, meaning ``|E(F_p)| | j``.  Since ``j`` is prime and
/// ``j > B1``, this is the stage-2 detection event.
///
/// Walks odd multiples of ``Q`` using repeated differential addition with step
/// ``2Q``:
/// - ``j = 1``: ``Q``
/// - ``j = 3``: ``3Q = diff_add(Q, 2Q, -Q)``  (``-Q`` has same x as ``Q``)
/// - ``j = 5``: ``5Q = diff_add(3Q, 2Q, Q)``
/// - ``j = 2k+1``: ``diff_add((2k-1)Q, 2Q, (2k-3)Q)``
///
/// ``stage2_primes`` is a pre-sorted list of primes in ``(B1, B2]``, computed
/// once by the caller to avoid redundant primality testing across curves.
///
/// Returns ``Ok(())`` if no factor found, ``Err(factor)`` if a factor is found.
fn run_stage2(
    q: &MontPoint,
    a24: &Uint<4>,
    stage2_primes: &[u64],
    n: &Uint<4>,
) -> Result<(), Uint<4>> {
    if stage2_primes.is_empty() {
        return Ok(());
    }

    let b2 = *stage2_primes.last().unwrap();

    // Precompute 2*Q for stepping.
    let two_q = point_double(q, a24, n);

    // Walk odd multiples of Q: j = 1, 3, 5, ...
    // Track cur = j*Q and prev_cur = (j-2)*Q for differential addition.
    //
    // For j=3: diff_add(Q, 2Q, -Q).  Since -Q has the same x-coordinate as Q
    // in Montgomery x-only coordinates, we use Q as the difference argument.
    // This is correct: diff_add(P, Q, P-Q) with P=Q, Q_arg=2Q, P-Q_arg=-Q≡Q.
    let mut prev_cur = q.clone(); // (j-2)*Q; for j=3 this is -Q ≡ Q (same x)
    let mut cur = diff_add(q, &two_q, q, n); // j=3: 3Q
    let mut j: u64 = 3;

    // Accumulate product of Z-coordinates for prime j values in (B1, B2].
    let mut acc = Uint::<4>::ONE;

    // Index into the sorted stage2_primes list for efficient lookup.
    let mut prime_idx = 0usize;

    loop {
        if j > b2 || prime_idx >= stage2_primes.len() {
            break;
        }

        // Advance prime_idx to the current j (primes are sorted).
        while prime_idx < stage2_primes.len() && stage2_primes[prime_idx] < j {
            prime_idx += 1;
        }

        // Accumulate Z_{jQ} if j is a stage-2 prime.
        if prime_idx < stage2_primes.len() && stage2_primes[prime_idx] == j {
            acc = mul_mod(&acc, &cur.z, n);
            prime_idx += 1;
        }

        // Step to j+2: next = diff_add(cur, 2Q, prev_cur).
        let next = diff_add(&cur, &two_q, &prev_cur, n);
        prev_cur = cur;
        cur = next;
        j += 2;
    }

    // Check gcd(acc, N).
    let g = gcd(acc, *n);
    if g > Uint::<4>::ONE && g < *n {
        return Err(g);
    }

    Ok(())
}

// ── Public interface ──────────────────────────────────────────────────────────

/// Result of an ECM run: either a factor was found or the run was inconclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EcmResult {
    /// A proper factor of N was found.
    Factor(Uint<4>),
    /// ECM completed without finding a factor (try another curve).
    Inconclusive,
}

/// Run ECM on ``n`` with a single curve parameterized by ``sigma``.
///
/// Stage 1 uses bound ``b1``; stage 2 uses bound ``b2`` (set ``b2 = b1`` to
/// skip stage 2).  Returns ``EcmResult::Factor(p)`` if a proper factor ``p``
/// is found, or ``EcmResult::Inconclusive`` if not.
///
/// For repeated calls with the same bounds (e.g., from ``ecm_factor``), prefer
/// the internal ``ecm_one_curve_with_primes`` variant to avoid recomputing the
/// stage-2 prime list on every call.
///
/// :param n: The integer to factor.  Must be odd and composite.
/// :param sigma: Suyama curve parameter.  Different values give different curves.
/// :param b1: Stage-1 smoothness bound.
/// :param b2: Stage-2 smoothness bound.  Set equal to ``b1`` to skip stage 2.
/// :returns: ``EcmResult::Factor(p)`` or ``EcmResult::Inconclusive``.
pub fn ecm_one_curve(n: &Uint<4>, sigma: u64, b1: u64, b2: u64) -> EcmResult {
    // Precompute stage-2 primes for this call.
    let stage2_primes: Vec<u64> = if b2 > b1 {
        factor_base_up_to(b2).into_iter().filter(|&p| p > b1).collect()
    } else {
        Vec::new()
    };
    ecm_one_curve_inner(n, sigma, b1, &stage2_primes)
}

/// Internal: run ECM on ``n`` with a single curve, using a precomputed stage-2 prime list.
///
/// ``stage2_primes`` must be a sorted list of primes in ``(B1, B2]``.  This
/// avoids recomputing the prime list on every curve when called from
/// ``ecm_factor``.
fn ecm_one_curve_inner(
    n: &Uint<4>,
    sigma: u64,
    b1: u64,
    stage2_primes: &[u64],
) -> EcmResult {
    // Trivial cases.
    if *n <= Uint::<4>::ONE {
        return EcmResult::Inconclusive;
    }

    // Generate the curve and starting point via Suyama parameterization.
    let (a24, p) = match suyama_param(sigma, n) {
        SuyamaResult::Ok { a24, p } => (a24, p),
        SuyamaResult::Factor(g) => {
            if g > Uint::<4>::ONE && g < *n {
                return EcmResult::Factor(g);
            }
            return EcmResult::Inconclusive;
        }
    };

    // Stage 1.
    let q = match run_stage1(&p, &a24, b1, n) {
        Ok(q) => q,
        Err(g) => {
            if g > Uint::<4>::ONE && g < *n {
                return EcmResult::Factor(g);
            }
            return EcmResult::Inconclusive;
        }
    };

    // If stage 1 produced the identity, the curve is degenerate for this N.
    if q.is_infinity() {
        return EcmResult::Inconclusive;
    }

    // Stage 2.
    if !stage2_primes.is_empty() {
        match run_stage2(&q, &a24, stage2_primes, n) {
            Ok(()) => {}
            Err(g) => {
                if g > Uint::<4>::ONE && g < *n {
                    return EcmResult::Factor(g);
                }
            }
        }
    }

    EcmResult::Inconclusive
}

/// Why [`ecm_factor`] did not return a divisor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EcmError {
    /// ``n ≤ 1`` has no non-trivial factorization to ask for.
    InvalidInput,
    /// ``n`` is prime, certified by a primality test that is proven on ``n``: ``n = 2``, or odd
    /// and below [`MILLER_RABIN_PROOF_BOUND`].
    Prime,
    /// ``n`` is odd, at least [`MILLER_RABIN_PROOF_BOUND`], and "probably prime", which cannot
    /// be certified.
    Unsupported,
    /// Every curve tried found no divisor or only a trivial one (``1`` or ``n``).  This is a
    /// resource outcome and proves nothing about ``n``.
    CurvesExhausted {
        /// Number of curves that were run and failed.
        curves: u32,
    },
    /// A curve produced a candidate ``d`` that is neither ``1`` nor ``n`` yet does not divide
    /// ``n``.  A gcd with ``n`` always divides ``n``, so this is a breach of the divisor
    /// invariant in the arithmetic, not a failed curve; it is reported at once and never retried.
    Inconsistent {
        /// The candidate that failed ``n mod d == 0``.
        d: Uint<4>,
    },
}

/// Judge one curve's result: a proper divisor of ``n``, a failed curve, or an invariant breach.
///
/// ``Inconclusive`` and the trivial candidates ``1`` and ``n`` are failed curves (``Ok(None)``).
/// Any other candidate is checked as ``1 < d < n`` and ``n mod d == 0``, independently of how it
/// was found; one that fails is [`EcmError::Inconsistent`], never a failed curve.
fn judge_candidate(n: &Uint<4>, result: EcmResult) -> Result<Option<Uint<4>>, EcmError> {
    let EcmResult::Factor(d) = result else {
        return Ok(None);
    };
    if d == Uint::<4>::ONE || d == *n {
        return Ok(None);
    }
    let divides = NonZero::new(d)
        .into_option()
        .is_some_and(|nonzero| d < *n && n.rem(&nonzero) == Uint::<4>::ZERO);
    if divides {
        Ok(Some(d))
    } else {
        Err(EcmError::Inconsistent { d })
    }
}

/// Try curves ``sigma = 6, 7, …`` (at most ``max_curves``) through ``run_curve`` and return the
/// first verified divisor of ``n``.
///
/// A curve that finds nothing, or only a trivial candidate, counts as a failed curve; an
/// inconsistent candidate ends the search with [`EcmError::Inconsistent`].  ``n`` must be odd
/// and composite.
fn search_curves<C>(n: &Uint<4>, max_curves: u32, mut run_curve: C) -> Result<Uint<4>, EcmError>
where
    C: FnMut(u64) -> EcmResult,
{
    let mut curves = 0u32;
    for i in 0..max_curves {
        let sigma = 6u64 + u64::from(i);
        if let Some(d) = judge_candidate(n, run_curve(sigma))? {
            return Ok(d);
        }
        curves += 1;
    }
    Err(EcmError::CurvesExhausted { curves })
}

/// Attempt to factor ``n`` by running ECM with multiple curves.
///
/// Classifies ``n`` first, in this order: ``n ≤ 1`` is invalid; ``2`` is prime; an even
/// ``n > 2`` yields the divisor ``2`` with no curve run; an odd ``n`` is tested for primality,
/// where a "prime" verdict is [`EcmError::Prime`] below [`MILLER_RABIN_PROOF_BOUND`] (the
/// test is proven there) and [`EcmError::Unsupported`] at or above it.  An odd composite is then
/// attacked with up to ``max_curves`` curves, ``sigma`` starting from 6 (``sigma`` 0 through 5
/// are degenerate in the Suyama parameterization).  A divisor is returned only after it passes
/// ``1 < d < n`` and ``n mod d == 0`` in ``Uint`` arithmetic.  A curve that yields no candidate,
/// or only ``1`` or ``n``, counts as a failed curve; any other candidate that fails the check
/// is [`EcmError::Inconsistent`] and ends the search.
///
/// The stage-2 prime list is precomputed once and reused across all curves.
///
/// :param n: The integer to factor.
/// :param b1: Stage-1 smoothness bound.
/// :param b2: Stage-2 smoothness bound.  Set equal to ``b1`` to skip stage 2.
/// :param max_curves: Maximum number of curves to try.
/// :returns: A proper divisor of ``n``.
///
/// # Errors
///
/// * [`EcmError::InvalidInput`] — ``n ≤ 1``.
/// * [`EcmError::Prime`] — ``n`` is prime and the primality verdict is proven.
/// * [`EcmError::Unsupported`] — ``n`` is odd, at least [`MILLER_RABIN_PROOF_BOUND`], and
///   "probably prime".
/// * [`EcmError::CurvesExhausted`] — every curve failed; ``curves`` counts the curves that ran
///   (``0`` when ``max_curves = 0``).  This proves nothing about ``n``.
/// * [`EcmError::Inconsistent`] — a curve produced a candidate other than ``1`` or ``n`` that
///   does not divide ``n``; reported at once, without trying further curves.
pub fn ecm_factor(n: &Uint<4>, b1: u64, b2: u64, max_curves: u32) -> Result<Uint<4>, EcmError> {
    let two = Uint::<4>::from(2u64);
    if *n <= Uint::<4>::ONE {
        return Err(EcmError::InvalidInput);
    }
    if *n == two {
        return Err(EcmError::Prime);
    }
    // An even n > 2 is divisible by 2 (bit 0 is clear), which is proper because n > 2.
    if !bool::from(n.bit(0)) {
        return Ok(two);
    }
    if is_prime(n) {
        // "Composite" verdicts are always proofs; "prime" verdicts only below the bound.
        let proven = *n < Uint::<4>::from(MILLER_RABIN_PROOF_BOUND);
        return Err(if proven { EcmError::Prime } else { EcmError::Unsupported });
    }

    // Precompute stage-2 primes once for all curves.
    let stage2_primes: Vec<u64> = if b2 > b1 {
        factor_base_up_to(b2).into_iter().filter(|&p| p > b1).collect()
    } else {
        Vec::new()
    };

    search_curves(n, max_curves, |sigma| ecm_one_curve_inner(n, sigma, b1, &stage2_primes))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn u(v: u64) -> Uint<4> {
        Uint::<4>::from(v)
    }

    /// Run ``ecm_factor`` on ``n`` and require a proper prime divisor of it.
    fn prime_divisor_of(n: &Uint<4>, b1: u64, b2: u64, max_curves: u32) -> Uint<4> {
        let d = ecm_factor(n, b1, b2, max_curves)
            .unwrap_or_else(|e| panic!("ecm_factor failed on {n:?}: {e:?}"));
        assert!(d > u(1) && d < *n, "factor must be proper");
        assert_eq!(n.rem(&NonZero::new(d).unwrap()), u(0), "factor must divide N");
        assert!(is_prime(&d), "factor must be prime");
        d
    }

    // ── Unit tests: modular arithmetic helpers ────────────────────────────────

    /// Verify ``mul_mod`` for small values.
    #[test]
    fn mul_mod_small() {
        let n = u(13);
        // 5 * 6 = 30 ≡ 4 (mod 13)
        assert_eq!(mul_mod(&u(5), &u(6), &n), u(4));
        // 12 * 12 = 144 ≡ 1 (mod 13)
        assert_eq!(mul_mod(&u(12), &u(12), &n), u(1));
    }

    /// Verify ``gcd`` against known values.
    #[test]
    fn gcd_known_values() {
        assert_eq!(gcd(u(0), u(7)), u(7));
        assert_eq!(gcd(u(7), u(0)), u(7));
        assert_eq!(gcd(u(12), u(8)), u(4));
        assert_eq!(gcd(u(35), u(14)), u(7));
        assert_eq!(gcd(u(17), u(13)), u(1));
        // gcd(143, 11) = 11
        assert_eq!(gcd(u(143), u(11)), u(11));
    }

    /// Verify ``try_inv`` returns the correct inverse when gcd = 1.
    #[test]
    fn try_inv_prime_modulus() {
        let n = u(13);
        // 5^{-1} mod 13 = 8 (since 5*8 = 40 ≡ 1 mod 13)
        let inv = try_inv(&u(5), &n).expect("inverse should exist");
        assert_eq!(mul_mod(&u(5), &inv, &n), u(1));
    }

    /// Verify ``try_inv`` returns a factor when gcd > 1.
    #[test]
    fn try_inv_finds_factor() {
        // N = 11 * 13 = 143.  gcd(11, 143) = 11.
        let n = u(143);
        let result = try_inv(&u(11), &n);
        assert!(result.is_err(), "should find a factor");
        let g = result.unwrap_err();
        assert!(g > u(1) && g < n, "factor should be proper");
        assert_eq!(u(143).rem(&NonZero::new(g).unwrap()), u(0), "factor should divide N");
    }

    // ── KAT: stage 1 alone ────────────────────────────────────────────────────

    /// ECM factors N = 143 = 11 * 13 with B1 = 11.
    ///
    /// With B1 = 11, the stage-1 scalar includes 11 as a prime power.
    /// The group order |E(F_11)| for a suitable curve is 11-smooth, so
    /// stage 1 should find a factor.
    #[test]
    fn ecm_factors_143() {
        prime_divisor_of(&u(143), 11, 11, 20);
    }

    /// ECM factors N = 28891 = 167 * 173 with B1 = 200.
    #[test]
    fn ecm_factors_28891() {
        prime_divisor_of(&u(28891), 200, 200, 50); // 167 * 173
    }

    /// ECM factors a semiprime with two 3-digit prime factors.
    ///
    /// N = 101 * 103 = 10403.  Both factors are prime and near 100.
    /// With B1 = 100 and B2 = 200, ECM should find a factor.
    #[test]
    fn ecm_factors_medium_semiprime() {
        // 101 and 103 are both prime.
        prime_divisor_of(&u(101u64 * 103u64), 100, 200, 100);
    }

    /// ECM factors a semiprime with two 3-digit prime factors.
    ///
    /// N = 211 * 223 = 47053.  Both factors are prime.
    /// With B1 = 200 and B2 = 400, ECM should find a factor.
    #[test]
    fn ecm_factors_near_2_32() {
        // 211 and 223 are both prime.
        prime_divisor_of(&u(211u64 * 223u64), 200, 400, 100);
    }

    // ── KAT: stage 1 + stage 2 ───────────────────────────────────────────────

    /// Stage 1 alone fails but stage 1 + stage 2 succeeds.
    ///
    /// We engineer an instance where |E(F_p)| has a prime factor just above B1
    /// but ≤ B2.
    ///
    /// Construction: N = 1009 * 3541 = 3,572,869.  With B1 = 20 and B2 = 200,
    /// stage 2 catches curves whose group order has a prime factor in (20, 200].
    /// Since the group order of E(F_1009) is roughly 1009 ± O(sqrt(1009)) ≈ 1009,
    /// many curves will have a prime factor in (20, 200] in their group order.
    ///
    /// The test verifies that stage 1+2 succeeds (the key property).  Stage 1
    /// alone may or may not succeed depending on which curves are tried.
    #[test]
    fn stage2_finds_factor_stage1_misses() {
        // N = 1009 * 3541 = 3,572,869.
        let n = u(1009u64 * 3541u64);
        // B1 = 20: stage 1 only catches curves where |E(F_p)| is 20-smooth.
        // B2 = 200: stage 2 catches primes up to 200.  Stage 1+2 must succeed.
        prime_divisor_of(&n, 20, 200, 50);
    }

    /// Verify that stage 2 with a tight B1 and wider B2 finds factors.
    ///
    /// N = 101 * 103 = 10403.  With B1 = 10 (very tight), stage 1 alone rarely
    /// succeeds.  With B2 = 200, stage 2 should find a factor.
    #[test]
    fn stage2_extends_reach() {
        // 101 and 103 are both prime.
        let n = u(101u64 * 103u64);
        // B1 = 10: very tight — stage 1 only catches 10-smooth group orders.
        // B2 = 200: stage 2 catches primes up to 200.
        prime_divisor_of(&n, 10, 200, 50);
    }

    // ── Stage 1: one ladder call per prime power ──────────────────────────────

    /// Model of a single-product stage-1 scalar accumulated in a fixed 256-bit integer.
    ///
    /// Returns ``Some(k)`` with ``k = ∏ q^⌊log_q(B1)⌋`` only while the exact product fits in
    /// ``Uint<4>``, and ``None`` once it would not.  `wrapping` selects the overflowing variant
    /// (the product reduced mod 2^256), returned as ``Some`` unconditionally.
    fn single_product_scalar(b1: u64, wrapping: bool) -> Option<Uint<4>> {
        let mut k = Uint::<4>::ONE;
        for p in factor_base_up_to(b1) {
            let q = Uint::<4>::from(largest_prime_power_le(p, b1));
            let (lo, hi) = k.widening_mul(&q);
            if hi != Uint::<4>::ZERO && !wrapping {
                return None;
            }
            k = lo;
        }
        Some(k)
    }

    /// ``true`` when ``a`` and ``b`` are the same projective point mod ``n``: ``X1 Z2 = X2 Z1``.
    fn projectively_equal(a: &MontPoint, b: &MontPoint, n: &Uint<4>) -> bool {
        mul_mod(&a.x, &b.z, n) == mul_mod(&b.x, &a.z, n)
    }

    /// ``largest_prime_power_le`` returns ``p^⌊log_p(B1)⌋`` on known values and at the extremes.
    #[test]
    fn largest_prime_power_known_values() {
        assert_eq!(largest_prime_power_le(2, 100), 64);
        assert_eq!(largest_prime_power_le(3, 100), 81);
        assert_eq!(largest_prime_power_le(5, 100), 25);
        assert_eq!(largest_prime_power_le(7, 100), 49);
        assert_eq!(largest_prime_power_le(97, 100), 97);
        // The power equals b1 exactly, and the product would overflow a naive `q * p` check.
        assert_eq!(largest_prime_power_le(2, 64), 64);
        assert_eq!(largest_prime_power_le(3, u64::MAX), 12_157_665_459_056_928_801);
    }

    /// For ``B1`` up to 178 (every fifth bound, plus 177 and 178) the per-prime-power point
    /// equals the single-product ladder point (as a projective point), because the product
    /// still fits in 256 bits there.
    ///
    /// ``n`` is a product of two 61- and 89-bit Mersenne primes, so no factor can appear at
    /// these tiny bounds and ``run_stage1`` always returns its point.
    #[test]
    fn stage1_per_prime_power_equals_single_product_ladder() {
        let n = ((Uint::<4>::ONE << 61u32).wrapping_sub(&Uint::<4>::ONE))
            .widening_mul(&(Uint::<4>::ONE << 89u32).wrapping_sub(&Uint::<4>::ONE))
            .0;
        for sigma in [6u64, 7, 123] {
            let (a24, p) = match suyama_param(sigma, &n) {
                SuyamaResult::Ok { a24, p } => (a24, p),
                SuyamaResult::Factor(_) => panic!("unexpected factor in curve construction"),
            };
            for b1 in (2..=178u64).step_by(5).chain([177, 178]) {
                let k = single_product_scalar(b1, false)
                    .unwrap_or_else(|| panic!("product must fit in 256 bits at B1 = {b1}"));
                let reference = montgomery_ladder(&p, &k, &a24, &n);
                let q = run_stage1(&p, &a24, b1, &n)
                    .unwrap_or_else(|g| panic!("no factor expected at B1 = {b1}, got {g:?}"));
                assert!(
                    projectively_equal(&q, &reference, &n),
                    "sigma = {sigma}, B1 = {b1}: per-prime point differs from single product"
                );
            }
        }
    }

    /// The exact product of prime powers stops fitting in 256 bits at ``B1 = 179``, where the
    /// prime 179 first multiplies it past ``2^256``.
    #[test]
    fn single_product_scalar_overflows_at_179() {
        assert!(single_product_scalar(178, false).is_some());
        assert!(single_product_scalar(179, false).is_none());
    }

    /// A factor that only appears once prime powers above the single-product overflow point are
    /// applied.
    ///
    /// ``n = (2^61 − 1) · 1009`` with Suyama ``sigma = 15``.  The search that found it ran
    /// primes ``r ∈ (1000, 3000]`` in increasing order, ``n = (2^61 − 1) · r``,
    /// ``sigma ∈ [6, 300)``, and ``B1 ∈ {200, 250, 300, 400}``, keeping the first triple where
    /// stage 1 (i) finds nothing at ``B1 = 178`` (the last bound whose scalar fits in 256 bits),
    /// (ii) finds exactly ``r`` at ``B1``, and (iii) finds nothing when the scalar is the product
    /// reduced mod 2^256.  The first hit is ``r = 1009``, ``sigma = 15``, ``B1 = 250``.
    #[test]
    fn stage1_prime_power_above_overflow_point_finds_factor() {
        let n = ((Uint::<4>::ONE << 61u32).wrapping_sub(&Uint::<4>::ONE))
            .widening_mul(&u(1009))
            .0;
        let (sigma, b1) = (15u64, 250u64);
        let (a24, p) = match suyama_param(sigma, &n) {
            SuyamaResult::Ok { a24, p } => (a24, p),
            SuyamaResult::Factor(_) => panic!("unexpected factor in curve construction"),
        };

        // Within the exact single-product range stage 1 does not reach the identity mod 1009.
        assert!(run_stage1(&p, &a24, 178, &n).is_ok());

        // The truncated product loses the factor...
        let wrapped = single_product_scalar(b1, true).unwrap();
        let old = montgomery_ladder(&p, &wrapped, &a24, &n);
        let g_old = gcd(old.z, n);
        assert!(g_old == u(1) || g_old == n, "truncated scalar must not find a factor");

        // ...and the per-prime-power ladder finds 1009, through the public entry points too.
        assert_eq!(run_stage1(&p, &a24, b1, &n).err(), Some(u(1009)));
        assert_eq!(ecm_one_curve(&n, sigma, b1, b1), EcmResult::Factor(u(1009)));
    }

    // ── Property tests ────────────────────────────────────────────────────────

    /// For a range of small semiprimes, ecm_factor always finds a proper prime factor.
    #[test]
    fn ecm_factors_small_semiprimes() {
        // Small semiprimes: p * q where p, q are small primes.
        let semiprimes: &[(u64, u64)] = &[
            (11, 13),
            (17, 19),
            (23, 29),
            (31, 37),
            (41, 43),
            (53, 59),
            (97, 101),
            (127, 131),
            (251, 257),
        ];
        for &(p, q) in semiprimes {
            let d = prime_divisor_of(&u(p * q), 300, 300, 50);
            assert!(d == u(p) || d == u(q), "divisor of {p} * {q} must be a prime factor");
        }
    }

    // ── Classification of inputs ──────────────────────────────────────────────

    /// ``n ≤ 1`` has no factorization to ask for.
    #[test]
    fn zero_and_one_are_invalid_input() {
        for v in [0u64, 1] {
            assert_eq!(ecm_factor(&u(v), 100, 100, 10), Err(EcmError::InvalidInput));
        }
    }

    /// ``n = 2`` is prime; it is not a trivial-divisor success.
    #[test]
    fn two_is_prime() {
        assert_eq!(ecm_factor(&u(2), 100, 100, 10), Err(EcmError::Prime));
    }

    /// An even ``n > 2`` yields the divisor 2 even when no curve may run.
    #[test]
    fn even_n_yields_divisor_two() {
        for v in [4u64, 10, 1 << 40] {
            assert_eq!(ecm_factor(&u(v), 100, 100, 0), Ok(u(2)));
        }
    }

    /// Odd primes below the proven bound are certified ``Prime`` and no curve is spent.
    #[test]
    fn odd_primes_are_prime() {
        for v in [3u64, 5, 1_000_003, 2_305_843_009_213_693_951] {
            assert_eq!(ecm_factor(&u(v), 100, 100, 10), Err(EcmError::Prime), "n = {v}");
        }
    }

    /// 2^89 − 1 is a Mersenne prime above the proven bound: it tests prime but that verdict
    /// is unproven, so the input is unsupported rather than ``Prime``.
    #[test]
    fn prime_above_proof_bound_is_unsupported() {
        let n = (Uint::<4>::ONE << 89u32).wrapping_sub(&Uint::<4>::ONE);
        assert!(n >= Uint::<4>::from(MILLER_RABIN_PROOF_BOUND));
        assert_eq!(ecm_factor(&n, 100, 100, 10), Err(EcmError::Unsupported));
    }

    /// ψ₁₂ is composite, so ``ecm_factor`` must never certify it ``Prime``; any divisor it
    /// does return is checked.
    #[test]
    fn composite_psi_12_is_never_prime() {
        let n = Uint::<4>::from(399_165_290_221u128 * 798_330_580_441u128);
        match ecm_factor(&n, 100, 100, 3) {
            Ok(d) => assert_eq!(n.rem(&NonZero::new(d).unwrap()), u(0)),
            Err(e) => assert_eq!(e, EcmError::CurvesExhausted { curves: 3 }),
        }
    }

    /// With no curves allowed, an odd composite is exhausted at zero curves.
    #[test]
    fn zero_max_curves_is_exhausted_at_zero() {
        assert_eq!(
            ecm_factor(&u(143), 11, 11, 0),
            Err(EcmError::CurvesExhausted { curves: 0 })
        );
    }

    /// Exhaustion counts the curves that ran, and is distinct from every other outcome.
    #[test]
    fn exhaustion_counts_curves_run() {
        // Two 20-bit primes with B1 = B2 = 2: no curve can reach either group order.
        let n = u(1_000_003u64 * 1_000_033u64);
        assert_eq!(ecm_factor(&n, 2, 2, 3), Err(EcmError::CurvesExhausted { curves: 3 }));
    }

    /// No candidate and the trivial candidates ``1`` and ``n`` are failed curves; proper
    /// divisors pass.
    #[test]
    fn judge_candidate_classifies_trivial_and_proper() {
        let n = u(143);
        assert_eq!(judge_candidate(&n, EcmResult::Inconclusive), Ok(None));
        assert_eq!(judge_candidate(&n, EcmResult::Factor(u(1))), Ok(None));
        assert_eq!(judge_candidate(&n, EcmResult::Factor(n)), Ok(None));
        assert_eq!(judge_candidate(&n, EcmResult::Factor(u(11))), Ok(Some(u(11))));
        assert_eq!(judge_candidate(&n, EcmResult::Factor(u(13))), Ok(Some(u(13))));
    }

    /// Any other candidate is an invariant breach carrying the candidate, including ``0`` and
    /// values above ``n``.
    #[test]
    fn judge_candidate_flags_non_dividing_candidates() {
        let n = u(143);
        for v in [0u64, 2, 7, 12, 142, 144, u64::MAX] {
            assert_eq!(
                judge_candidate(&n, EcmResult::Factor(u(v))),
                Err(EcmError::Inconsistent { d: u(v) }),
                "d = {v}"
            );
        }
    }

    /// A non-dividing candidate ends the curve search at once, without running more curves.
    #[test]
    fn search_surfaces_inconsistent_candidate_without_retrying() {
        let mut runs = 0u32;
        let result = search_curves(&u(143), 10, |_| {
            runs += 1;
            EcmResult::Factor(u(7))
        });
        assert_eq!(result, Err(EcmError::Inconsistent { d: u(7) }));
        assert_eq!(runs, 1);
    }

    /// An inconsistent candidate after failed curves still ends the search; sigma starts at 6.
    #[test]
    fn search_inconsistent_after_failed_curves() {
        let mut sigmas = Vec::new();
        let result = search_curves(&u(143), 10, |sigma| {
            sigmas.push(sigma);
            if sigma == 8 { EcmResult::Factor(u(5)) } else { EcmResult::Inconclusive }
        });
        assert_eq!(result, Err(EcmError::Inconsistent { d: u(5) }));
        assert_eq!(sigmas, vec![6, 7, 8]);
    }

    /// A candidate equal to ``n`` (or ``1``) consumes a curve like ``Inconclusive`` does and does
    /// not stop the search.
    #[test]
    fn search_trivial_candidates_consume_curves() {
        let n = u(143);
        assert_eq!(
            search_curves(&n, 4, |_| EcmResult::Factor(n)),
            Err(EcmError::CurvesExhausted { curves: 4 })
        );
        assert_eq!(
            search_curves(&n, 4, |s| if s % 2 == 0 {
                EcmResult::Factor(u(1))
            } else {
                EcmResult::Factor(n)
            }),
            Err(EcmError::CurvesExhausted { curves: 4 })
        );
        assert_eq!(
            search_curves(&n, 4, |_| EcmResult::Inconclusive),
            Err(EcmError::CurvesExhausted { curves: 4 })
        );
    }

    /// A proper divisor from a curve is returned after the failed curves before it.
    #[test]
    fn search_returns_verified_divisor() {
        let n = u(143);
        let result = search_curves(&n, 10, |s| {
            EcmResult::Factor(if s == 7 { u(13) } else { n })
        });
        assert_eq!(result, Ok(u(13)));
    }
}
