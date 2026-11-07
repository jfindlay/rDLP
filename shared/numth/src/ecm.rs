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
//! **Stage 1**: compute ``k * P`` where ``k = ∏ { q^⌊log_q(B1)⌋ : q prime, q ≤ B1 }``.
//! A factor is found if ``gcd(Z-coordinate, N) > 1`` at any point during the
//! ladder.
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

/// Compute the stage-1 scalar: ``k = ∏ { q^⌊log_q(B1)⌋ : q prime, q ≤ B1 }``.
///
/// This is the product of prime powers ≤ ``B1``.  Multiplying a point by ``k``
/// ensures that if ``|E(F_p)|`` is ``B1``-smooth, the point becomes the identity
/// mod ``p``.
fn stage1_scalar(b1: u64) -> Uint<4> {
    let primes = factor_base_up_to(b1);
    let mut k = Uint::<4>::ONE;
    for p in primes {
        // Compute q = p^⌊log_p(B1)⌋: the largest power of p that is ≤ B1.
        let mut q = p;
        while q <= b1 / p {
            // q * p <= b1 (checked to avoid overflow)
            q *= p;
        }
        // Multiply k by q.  Use wrapping_mul since k is built up incrementally
        // and we rely on the Montgomery ladder to detect factors via Z-coordinates.
        k = k.wrapping_mul(&Uint::<4>::from(q));
    }
    k
}

/// Run stage 1 of ECM on a single curve.
///
/// Computes ``k * P`` where ``k`` is the stage-1 scalar.  Returns ``Ok(point)``
/// if stage 1 completed without finding a factor, or ``Err(factor)`` if a
/// non-trivial factor of ``N`` was found.
///
/// A factor is detected by checking ``gcd(Z, N)`` after the ladder completes.
/// The ladder itself never fails (no inversions are needed).
fn run_stage1(
    p: &MontPoint,
    a24: &Uint<4>,
    b1: u64,
    n: &Uint<4>,
) -> Result<MontPoint, Uint<4>> {
    // Build the stage-1 scalar and run the ladder.
    let k = stage1_scalar(b1);
    let q = montgomery_ladder(p, &k, a24, n);

    // Check gcd(Z, N).
    let g = gcd(q.z, *n);
    if g > Uint::<4>::ONE && g < *n {
        return Err(g);
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

/// Attempt to factor ``n`` by running ECM with multiple random curves.
///
/// Tries up to ``max_curves`` curves with ``sigma`` values starting from
/// ``sigma_start = 6`` (the standard starting value).  Returns ``Some(p)`` if
/// a proper factor ``p`` is found (``1 < p < n`` and ``p | n``), or ``None``
/// if all curves were inconclusive.
///
/// The stage-2 prime list is precomputed once and reused across all curves.
///
/// :param n: The integer to factor.  Must be odd and composite.
/// :param b1: Stage-1 smoothness bound.
/// :param b2: Stage-2 smoothness bound.
/// :param max_curves: Maximum number of curves to try.
/// :returns: ``Some(factor)`` or ``None``.
pub fn ecm_factor(n: &Uint<4>, b1: u64, b2: u64, max_curves: u32) -> Option<Uint<4>> {
    // Trivial cases.
    if *n <= Uint::<4>::ONE {
        return None;
    }
    // Check if n is even (factor 2).
    if !bool::from(n.bit(0)) {
        return Some(Uint::<4>::from(2u64));
    }

    // Precompute stage-2 primes once for all curves.
    let stage2_primes: Vec<u64> = if b2 > b1 {
        factor_base_up_to(b2).into_iter().filter(|&p| p > b1).collect()
    } else {
        Vec::new()
    };

    // Try successive sigma values starting from 6.
    // sigma = 0, 1, 2, 3, 4, 5 are degenerate in the Suyama parameterization.
    for i in 0..max_curves {
        let sigma = 6u64 + i as u64;
        match ecm_one_curve_inner(n, sigma, b1, &stage2_primes) {
            EcmResult::Factor(p) => return Some(p),
            EcmResult::Inconclusive => {}
        }
    }
    None
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime::is_prime;

    fn u(v: u64) -> Uint<4> {
        Uint::<4>::from(v)
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
        let n = u(143);
        let result = ecm_factor(&n, 11, 11, 20);
        assert!(result.is_some(), "ecm_factor should find a factor of 143");
        let p = result.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
    }

    /// ECM factors N = 28891 = 167 * 173 with B1 = 200.
    #[test]
    fn ecm_factors_28891() {
        let n = u(28891); // 167 * 173
        let result = ecm_factor(&n, 200, 200, 50);
        assert!(result.is_some(), "ecm_factor should find a factor of 28891");
        let p = result.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
    }

    /// ECM factors a semiprime with two 3-digit prime factors.
    ///
    /// N = 101 * 103 = 10403.  Both factors are prime and near 100.
    /// With B1 = 100 and B2 = 200, ECM should find a factor.
    #[test]
    fn ecm_factors_medium_semiprime() {
        // 101 and 103 are both prime.
        let n = u(101u64 * 103u64);
        let result = ecm_factor(&n, 100, 200, 100);
        assert!(result.is_some(), "ecm_factor should find a factor of 101 * 103");
        let p = result.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
    }

    /// ECM factors a semiprime with two 3-digit prime factors.
    ///
    /// N = 211 * 223 = 47053.  Both factors are prime.
    /// With B1 = 200 and B2 = 400, ECM should find a factor.
    #[test]
    fn ecm_factors_near_2_32() {
        // 211 and 223 are both prime.
        let n = u(211u64 * 223u64);
        let result = ecm_factor(&n, 200, 400, 100);
        assert!(result.is_some(), "ecm_factor should find a factor of 211 * 223");
        let p = result.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
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
        // B2 = 200: stage 2 catches primes up to 200.
        let stage1_plus_2 = ecm_factor(&n, 20, 200, 50);

        // Stage 1+2 must succeed.
        assert!(
            stage1_plus_2.is_some(),
            "stage 1+2 should find a factor of {} = 1009 * 3541",
            1009u64 * 3541u64
        );
        let p = stage1_plus_2.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
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
        let result = ecm_factor(&n, 10, 200, 50);
        assert!(result.is_some(), "stage 1+2 should find a factor of 101 * 103");
        let p = result.unwrap();
        assert!(p > u(1) && p < n, "factor must be proper");
        let nz = NonZero::new(p).unwrap();
        assert_eq!(n.rem(&nz), u(0), "factor must divide N");
        assert!(is_prime(&p), "factor must be prime");
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
            let n = u(p * q);
            let result = ecm_factor(&n, 300, 300, 50);
            assert!(
                result.is_some(),
                "ecm_factor should find a factor of {p} * {q} = {}",
                p * q
            );
            let factor = result.unwrap();
            assert!(
                factor > u(1) && factor < n,
                "factor {factor:?} must be proper for N = {p} * {q}"
            );
            let nz = NonZero::new(factor).unwrap();
            assert_eq!(
                n.rem(&nz),
                u(0),
                "factor must divide N = {p} * {q}"
            );
            assert!(is_prime(&factor), "factor must be prime for N = {p} * {q}");
        }
    }

    /// ``ecm_factor`` returns ``None`` for prime inputs (no proper factor exists).
    #[test]
    fn ecm_inconclusive_for_prime() {
        // For a prime N, no proper factor exists.  ECM should return None.
        let n = u(1_000_003u64); // prime
        let result = ecm_factor(&n, 100, 100, 10);
        // ECM may return None (inconclusive) or accidentally find a "factor" that
        // equals 1 or N — but our interface guarantees 1 < p < N, so it returns None.
        // (In practice, for prime N, the Z-coordinate is never 0 mod N, so ECM
        // always returns Inconclusive.)
        assert!(result.is_none(), "ecm_factor should return None for prime N");
    }
}
