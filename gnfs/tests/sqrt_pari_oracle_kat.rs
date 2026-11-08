//! PARI/GP cross-checks of the algebraic facts behind the exact square root.
//!
//! Every test here is `#[ignore]`d: it needs `gp` (GP/PARI) on `PATH`, which the default test run
//! does not require. Run them by hand with
//!
//! ```text
//! cargo test -p gnfs --test sqrt_pari_oracle_kat -- --ignored
//! ```
//!
//! A test fails (it does not skip) when `gp` is missing.
//!
//! **Independence.** Each test computes one side from the crate's own public API (a verified
//! root, a nonsquare certificate, a modular image, a discriminant) and the other side from `gp -q`
//! run on the fixture inputs alone: the polynomial `f`, the relation pairs `(a, b)`, `m` and `N`.
//! PARI re-derives `γ = Π(a − b·α)`, `P = Π(a − b·m)`, the discriminant and the irreducibility of
//! `f` for itself. The crate never reads a PARI answer to produce its own, and PARI is never told
//! a result it is asked to "confirm", with one deliberate exception: the verified `β` is handed to
//! PARI as an input to be *checked* (`β² = γ` in `ℚ[y]/(f)`, and `β(m) mod N` recomputed), and
//! PARI's two roots of `t² − γ` are compared with it.
//!
//! **Fixtures** (all reused from the crate's own known-answer tests):
//!
//! - **S** (seam): `x² − 2`, `m = 23`, `N = 527`, relation (3, −2): `β = ±(1 + α)`, `y = ±24`,
//!   `Norm(β) = −1`, `24² = 576 ≡ 49 = P`.
//! - **R135**: the same field and modulus with relations (−3, 2), (−2, 1), (2, −1):
//!   `γ = 34 + 24α = (4 + 3α)²`.
//! - **NM** (non-monic): `2x² − 2x − 5`, `m = 5`, `N = 35`, relations (9, −11), (−1, −1):
//!   `β = ±(4 + θ/2)`.
//! - **C2** (non-monic): `2x² − 1`, `m = 2`, `N = 7`, relation (3, −4): `β = ±(1 + θ)`.
//! - **Cubic**: `x³ − 2`, `m = 7`, `N = 341`, relations (3, 1), (−1, −1).
//! - **V4**: `x⁴ − 10x² + 1`, `m = 4`, `N = 97`, relation (1, −1) twice: `β = ±(1 + α)`; the field
//!   `ℚ(√2, √3)` has no inert prime, so no reduction of `f` is irreducible.
//! - **45113 cubic**: `x³ + 15x² + 29x + 8`, `m = 31`, `N = 45113`, and the 22 relations of the
//!   dependency that the generated sieve-to-factor run selects.
//!
//! The nonsquare fixtures are the certified ones: `1 − α` and `3 − α` in `ℚ(√2)` at `m = 4`,
//! `N = 7`; `2 + α`, `10 + 7α` and `−5` in `ℚ(√2)` at `m = 23`; `2 + α` in `ℚ(√3)` (a norm-one
//! unit that is a square modulo many primes); and `2 − α` over the non-monic `3x² + x + 2`
//! (the field `ℚ(√−23)`).

use std::io::Write;
use std::process::{Command, Stdio};

use gnfs::sqrt::{
    AlgebraicSqrtOutcome, CrtLimits, DirectLimits, PrecisionLimits, PrimeSearchLimits,
    RecombinationLimits, RelationFactors, SqrtContext, SqrtLimits, VerifiedRoot,
    reconstruct_root_from_factors,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use shared_numfield::{
    FpPoly, IntPoly, NonresidueSearchLimits, PrimeModulus, RatPoly, RootSearchLimits, discriminant,
    is_prime_u64, resultant,
};

const GP_MISSING: &str = "PARI/GP (`gp`) must be on PATH to run the oracle tests";

// ─── fixtures ────────────────────────────────────────────────────────────────

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

/// A polynomial, a root `m` of it modulo `N`, and relation pairs, exactly as the crate takes them.
struct Fixture {
    name: &'static str,
    f: Vec<i64>,
    m: i64,
    n: i64,
    pairs: Vec<(i64, i64)>,
}

fn fixture(name: &'static str, f: &[i64], m: i64, n: i64, pairs: &[(i64, i64)]) -> Fixture {
    Fixture { name, f: f.to_vec(), m, n, pairs: pairs.to_vec() }
}

impl Fixture {
    fn context(&self) -> SqrtContext {
        SqrtContext::new(&poly(&self.f), &big(self.m), &big(self.n)).expect("valid context")
    }

    fn factors(&self) -> RelationFactors {
        let pairs: Vec<(BigInt, BigInt)> =
            self.pairs.iter().map(|&(a, b)| (big(a), big(b))).collect();
        RelationFactors::new(&self.context(), &pairs).expect("valid relation factors")
    }

    /// GP definitions of the fixture, derived by PARI from the inputs alone: the field
    /// polynomial `ff` in `y`, the monic `Ft = c^(d−1) ff(y/c)` that PARI's number field code
    /// needs, `γ = gg` and the rational product `pp`.
    fn gp_preamble(&self) -> String {
        let coefficients: Vec<String> = self.f.iter().map(i64::to_string).collect();
        let pairs: Vec<String> = self.pairs.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
        format!(
            "ff = Polrev([{}], 'y); mm = {}; nn = {}; rel = [{}];\n\
             cc = pollead(ff); dd = poldegree(ff); Ft = cc^(dd-1)*subst(ff, y, y/cc);\n\
             gg = prod(i=1, #rel, rel[i][1] - rel[i][2]*y);\n\
             pp = prod(i=1, #rel, rel[i][1] - rel[i][2]*mm);\n",
            coefficients.join(","),
            self.m,
            self.n,
            pairs.join(","),
        )
    }
}

fn fixture_s() -> Fixture {
    fixture("S", &[-2, 0, 1], 23, 527, &[(3, -2)])
}

fn fixture_r135() -> Fixture {
    fixture("R135", &[-2, 0, 1], 23, 527, &[(-3, 2), (-2, 1), (2, -1)])
}

fn fixture_nm() -> Fixture {
    fixture("NM", &[-5, -2, 2], 5, 35, &[(9, -11), (-1, -1)])
}

fn fixture_c2() -> Fixture {
    fixture("C2", &[-1, 0, 2], 2, 7, &[(3, -4)])
}

fn fixture_cubic() -> Fixture {
    fixture("x^3-2", &[-2, 0, 0, 1], 7, 341, &[(3, 1), (-1, -1)])
}

fn fixture_v4() -> Fixture {
    fixture("V4", &[1, 0, -10, 0, 1], 4, 97, &[(1, -1), (1, -1)])
}

/// The 22 relations of the dependency selected by the generated sieve-to-factor run on this
/// polynomial; each is a coprime pair with `b ≥ 1`.
fn fixture_45113() -> Fixture {
    fixture(
        "45113",
        &[8, 29, 15, 1],
        31,
        45113,
        &[
            (-61, 1),
            (-32, 1),
            (-3, 1),
            (2, 1),
            (-33, 2),
            (-25, 2),
            (-1, 2),
            (5, 2),
            (-17, 3),
            (-2, 3),
            (-19, 4),
            (-14, 5),
            (-15, 7),
            (7, 9),
            (-33, 17),
            (-5, 17),
            (1, 19),
            (16, 25),
            (-375, 29),
            (-1, 33),
            (-78, 37),
            (-5, 41),
        ],
    )
}

/// Fixtures with a verified algebraic root.
fn root_fixtures() -> Vec<Fixture> {
    vec![
        fixture_s(),
        fixture_r135(),
        fixture_nm(),
        fixture_c2(),
        fixture_cubic(),
        fixture_v4(),
        fixture_45113(),
    ]
}

/// Fixtures whose product is certified not to be a square in the field.
fn nonsquare_fixtures() -> Vec<Fixture> {
    vec![
        fixture("1-a", &[-2, 0, 1], 4, 7, &[(1, 1)]),
        fixture("3-a", &[-2, 0, 1], 4, 7, &[(3, 1)]),
        fixture("2+a", &[-2, 0, 1], 23, 527, &[(2, -1)]),
        fixture("10+7a", &[-2, 0, 1], 23, 527, &[(-3, 2), (-2, 1)]),
        fixture("-5", &[-2, 0, 1], 23, 527, &[(-5, 0)]),
        fixture("sqrt3 2+a", &[-3, 0, 1], 4, 13, &[(2, -1)]),
        fixture("sqrt-23 2-a", &[2, 1, 3], 7, 13, &[(2, 1)]),
    ]
}

// ─── limits ──────────────────────────────────────────────────────────────────

/// Direct reconstruction from the first prime, then split primes one at a time; the precision cap
/// is far above what any fixture needs.
fn ample() -> SqrtLimits {
    SqrtLimits {
        precision: PrecisionLimits { max_modulus_bits: 65_536 },
        direct: DirectLimits {
            prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 100 },
            nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
        },
        crt: CrtLimits {
            prime_search: PrimeSearchLimits { first_candidate: 0, max_primes_examined: 200 },
            max_primes: 1,
            root_search: RootSearchLimits { max_evaluations: 100_000 },
            nonresidue_search: NonresidueSearchLimits { max_candidates: 1000 },
        },
        recombination: RecombinationLimits { max_list_entries: 1 << 20, max_matches: 1 << 20 },
    }
}

fn verified_root(fixture: &Fixture) -> VerifiedRoot {
    match reconstruct_root_from_factors(fixture.factors(), ample()).expect("no invariant error") {
        AlgebraicSqrtOutcome::Verified(root) => root,
        other => panic!("{}: expected a verified root, got {other:?}", fixture.name),
    }
}

// ─── PARI plumbing ───────────────────────────────────────────────────────────

/// Run `gp -q` on `script` (fed on standard input) and return its output lines. Any text on
/// standard error (a GP error or warning) or a failing exit status is a test failure.
fn gp(script: &str) -> Vec<String> {
    let mut child = Command::new("gp")
        .arg("-q")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("{GP_MISSING}: {e}"));
    child
        .stdin
        .take()
        .expect("piped stdin")
        .write_all(script.as_bytes())
        .expect("write the script to gp");
    let output = child.wait_with_output().expect("wait for gp");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success() && stderr.trim().is_empty(),
        "gp failed ({}):\n{stderr}\nscript:\n{script}",
        output.status
    );
    String::from_utf8(output.stdout)
        .expect("gp output is UTF-8")
        .lines()
        .map(|l| l.trim().to_string())
        .collect()
}

fn gp_rational(q: &BigRational) -> String {
    format!("({}/{})", q.numer(), q.denom())
}

/// `Polrev([...], 'y)` for coefficients lowest degree first.
fn gp_poly(coeffs: &[BigRational]) -> String {
    let items: Vec<String> = coeffs.iter().map(gp_rational).collect();
    format!("Polrev([{}], 'y)", items.join(","))
}

/// A root `β` found in the `θ = c·α` basis, re-expressed in the `α` basis: the coefficient of
/// `θ^i` times `c^i`.
fn in_alpha_basis(beta: &RatPoly, c: &BigInt) -> Vec<BigRational> {
    let mut power = BigInt::one();
    let mut coeffs = Vec::with_capacity(beta.coeffs.len());
    for coeff in &beta.coeffs {
        coeffs.push(coeff * BigRational::from_integer(power.clone()));
        power *= c;
    }
    coeffs
}

fn parse_int(line: &str) -> BigInt {
    line.parse().unwrap_or_else(|e| panic!("`{line}` is not an integer: {e}"))
}

/// A prime `p` not dividing the leading coefficient at which `f` keeps its degree and is
/// irreducible: a certificate that `f` is irreducible over `Q`.
fn irreducibility_certificate(f: &IntPoly) -> Option<u64> {
    let degree = f.degree()?;
    let lead = f.leading_coeff()?;
    (3u64..1000).filter(|&p| is_prime_u64(p)).find(|&p| {
        let modulus = PrimeModulus::new(p).expect("odd prime modulus");
        let reduced = FpPoly::from_int_poly(f, modulus);
        !(lead % BigInt::from(p)).is_zero()
            && reduced.degree() == Some(degree)
            && reduced.is_irreducible()
    })
}

// ─── polynomials ─────────────────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_confirms_irreducibility_and_discriminant_of_the_fixture_polynomials() {
    let fixtures = [fixture_nm(), fixture_v4(), fixture_cubic(), fixture_45113(), fixture_c2()];
    for fixture in &fixtures {
        let context = fixture.context();
        let lines = gp(&format!(
            "{}print(polisirreducible(ff));\n\
             print(poldisc(Ft));\n\
             print(poldisc(Ft) == cc^((dd-1)*(dd-2))*poldisc(ff));\n",
            fixture.gp_preamble()
        ));
        assert_eq!(lines.len(), 3, "{}: {lines:?}", fixture.name);
        // PARI proves f irreducible over Q, and F(t) = c^(d−1) f(t/c) has the crate's
        // discriminant; disc(F) and disc(f) differ by the expected power of the leading
        // coefficient.
        assert_eq!(lines[0], "1", "{}: f is irreducible", fixture.name);
        assert_eq!(parse_int(&lines[1]), *context.discriminant(), "{}", fixture.name);
        assert_eq!(parse_int(&lines[1]), discriminant(context.monic_poly()), "{}", fixture.name);
        assert_eq!(lines[2], "1", "{}: disc(F) = c^((d−1)(d−2))·disc(f)", fixture.name);

        // The crate certifies irreducibility only through an inert prime. The V4 field has none
        // below 1000 (its Galois group has exponent two), so PARI is the only certificate there.
        let certificate = irreducibility_certificate(context.algebraic_poly());
        if fixture.name == "V4" {
            assert_eq!(certificate, None, "V4 is irreducible but reducible modulo every prime");
        } else {
            assert!(certificate.is_some(), "{}: the crate certifies f irreducible", fixture.name);
        }
    }
}

// ─── roots ───────────────────────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_confirms_verified_roots_square_to_their_targets() {
    for fixture in &root_fixtures() {
        let root = verified_root(fixture);
        let beta = in_alpha_basis(root.beta(), fixture.context().leading_coeff());
        let lines = gp(&format!(
            "{}bb = {};\n\
             print(polisirreducible(ff));\n\
             print(Mod(bb, ff)^2 == Mod(gg, ff));\n\
             nf = nfinit(Ft); gt = subst(gg, y, y/cc); bt = subst(bb, y, y/cc);\n\
             print(nfeltmul(nf, nfalgtobasis(nf, bt), nfalgtobasis(nf, bt)) == nfalgtobasis(nf, gt));\n\
             rts = nfroots(nf, x^2 - Mod(gt, Ft));\n\
             print(#rts);\n\
             print(#select(r -> nfalgtobasis(nf, r) == nfalgtobasis(nf, bt), rts));\n",
            fixture.gp_preamble(),
            gp_poly(&beta),
        ));
        // f irreducible; β² = γ in Q[y]/(f); the same through nfeltmul in the maximal order's
        // number field; and β is one of the two roots of t² − γ that PARI finds for itself.
        assert_eq!(lines, ["1", "1", "1", "2", "1"], "{}", fixture.name);
    }
}

// ─── non-squares ─────────────────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_confirms_certified_nonsquares_have_no_root_in_the_field() {
    for fixture in &nonsquare_fixtures() {
        let outcome = reconstruct_root_from_factors(fixture.factors(), ample())
            .expect("no invariant error");
        let AlgebraicSqrtOutcome::NotSquare(certificate) = outcome else {
            panic!("{}: expected a nonsquare certificate, got {outcome:?}", fixture.name);
        };
        let lines = gp(&format!(
            "{}print(polisirreducible(ff));\n\
             nf = nfinit(Ft); gt = subst(gg, y, y/cc);\n\
             print(#nfroots(nf, x^2 - Mod(gt, Ft)));\n",
            fixture.gp_preamble()
        ));
        assert_eq!(
            lines,
            ["1", "0"],
            "{}: PARI finds a root of t² − γ although the crate certified {certificate:?}",
            fixture.name
        );
    }
}

// ─── the map ─────────────────────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_recomputes_the_modular_image_and_the_congruence_from_each_root() {
    for fixture in &root_fixtures() {
        let root = verified_root(fixture);
        let n = big(fixture.n);
        let beta = in_alpha_basis(root.beta(), fixture.context().leading_coeff());
        let lines = gp(&format!(
            "{}bb = {};\n\
             yv = Mod(1, nn)*subst(bb, y, mm);\n\
             print(yv^2 == Mod(pp, nn));\n\
             print(lift(yv));\n\
             print(pp);\n",
            fixture.gp_preamble(),
            gp_poly(&beta),
        ));
        assert_eq!(lines.len(), 3, "{}: {lines:?}", fixture.name);
        // PARI evaluates β at α = m modulo N (every denominator is a unit) and finds y² ≡ P.
        assert_eq!(lines[0], "1", "{}: y² ≡ P (mod N)", fixture.name);
        // PARI's y and P are the crate's modular image and rational product.
        assert_eq!(&parse_int(&lines[1]), root.modular_image(), "{}", fixture.name);
        assert_eq!(&parse_int(&lines[2]), root.rational_product(), "{}", fixture.name);
        let y = root.modular_image();
        assert_eq!((y * y - root.rational_product()) % &n, BigInt::zero(), "{}", fixture.name);
    }
}

// ─── norm identity and the seam ──────────────────────────────────────────────

/// `F(a, b) = b^d f(a/b) = Σ f_i a^i b^(d−i)`.
fn homogeneous_value(f: &[i64], a: i64, b: i64) -> BigInt {
    let d = f.len() - 1;
    f.iter()
        .enumerate()
        .map(|(i, &c)| big(c) * big(a).pow(i as u32) * big(b).pow((d - i) as u32))
        .sum()
}

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_confirms_the_norm_identity_and_that_the_norm_is_not_the_image() {
    // F(a, b) = c·Norm(a − b·α), including the non-monic fixtures.
    for fixture in [fixture_c2(), fixture_nm(), fixture_s()] {
        let lines = gp(&format!(
            "{}for(i = 1, #rel, print(cc*norm(Mod(rel[i][1] - rel[i][2]*y, ff))));\n",
            fixture.gp_preamble()
        ));
        assert_eq!(lines.len(), fixture.pairs.len(), "{}", fixture.name);
        for (line, &(a, b)) in lines.iter().zip(&fixture.pairs) {
            assert_eq!(
                parse_int(line),
                homogeneous_value(&fixture.f, a, b),
                "{}: c·Norm({a} − {b}α)",
                fixture.name
            );
        }
    }

    // The seam: β = ±(1 + α) has Norm(β) = −1, the image is y = ±24, and 24² = 576 ≡ 49 = P
    // (mod 527); the norm would give 1² = 1 ≠ 49.
    let fixture = fixture_s();
    let root = verified_root(&fixture);
    let beta = in_alpha_basis(root.beta(), fixture.context().leading_coeff());
    assert!(beta.iter().all(|q| q.is_integer()), "the seam root is an algebraic integer");
    let beta_int = IntPoly::from_coeffs(beta.iter().map(|q| q.numer().clone()).collect());
    let crate_norm = resultant(&poly(&fixture.f), &beta_int);
    let lines = gp(&format!(
        "{}bb = {};\n\
         print(norm(Mod(bb, ff)));\n\
         yv = Mod(1, nn)*subst(bb, y, mm);\n\
         print(lift(yv));\n\
         print(yv^2 == Mod(pp, nn));\n\
         print(yv^2 == Mod(norm(Mod(bb, ff)), nn));\n",
        fixture.gp_preamble(),
        gp_poly(&beta),
    ));
    assert_eq!(lines.len(), 4, "{lines:?}");
    assert_eq!(parse_int(&lines[0]), big(-1), "Norm(β) = −1 in PARI");
    assert_eq!(parse_int(&lines[0]), crate_norm, "the crate's resultant agrees");
    assert_eq!(&parse_int(&lines[1]), root.modular_image());
    assert!([big(24), big(503)].contains(root.modular_image()));
    assert_eq!(lines[2], "1", "y² ≡ P (mod N)");
    assert_eq!(lines[3], "0", "y² is not congruent to Norm(β)");
}

// ─── the obstruction group ───────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_confirms_the_obstruction_group_of_two_quadratic_fields() {
    // ℚ(√3): trivial class group, unit rank 1, torsion {±1}. The group V of elements with all
    // valuations even, modulo squares, has dimension r₁ + r₂ + rk₂ Cl = 2 + 0 + 0.
    // 2 + √3 is the fundamental unit, a nonsquare of norm 1. That is the certified nonsquare
    // `2 + α` of the `x² − 3` fixture.
    let nonsquare_unit = nonsquare_fixtures().into_iter().nth(5).expect("six fixtures");
    assert_eq!(nonsquare_unit.name, "sqrt3 2+a");
    let outcome = reconstruct_root_from_factors(nonsquare_unit.factors(), ample())
        .expect("no invariant error");
    assert!(matches!(outcome, AlgebraicSqrtOutcome::NotSquare(_)), "{outcome:?}");
    let lines = gp(&format!(
        "{}bnf = bnfinit(ff, 1); nf = bnf.nf;\n\
         print(bnf.sign[1], \" \", bnf.sign[2]);\n\
         print(#bnf.fu);\n\
         print(bnf.tu[1]);\n\
         print(#select(c -> c % 2 == 0, bnf.cyc));\n\
         u = Mod(gg, ff);\n\
         print(norm(u));\n\
         print(#nfroots(nf, x^2 - u));\n\
         print(bnfisunit(bnf, u)[1] % 2);\n",
        nonsquare_unit.gp_preamble()
    ));
    assert_eq!(lines, ["2 0", "1", "2", "0", "1", "0", "1"]);
    assert_eq!(obstruction_dimension(&lines), 2, "dim V = r₁ + r₂ + rk₂ Cl for ℚ(√3)");

    // ℚ(√−5): class group Z/2, unit rank 0, torsion {±1}, so dim V = 0 + 1 + 1 = 2. The prime
    // 2 ramifies, (2) = J² with J = (2, 1 + √−5), and J is not principal.
    let lines = gp("ff = y^2 + 5;\n\
         bnf = bnfinit(ff, 1); nf = bnf.nf;\n\
         print(bnf.sign[1], \" \", bnf.sign[2]);\n\
         print(#bnf.fu);\n\
         print(bnf.tu[1]);\n\
         print(#select(c -> c % 2 == 0, bnf.cyc));\n\
         pr = idealprimedec(nf, 2);\n\
         print(#pr, \" \", pr[1].e, \" \", pr[1].f);\n\
         J = idealhnf(nf, 2, 1 + y);\n\
         print(J == idealhnf(nf, pr[1]));\n\
         print(idealpow(nf, J, 2) == idealhnf(nf, 2));\n\
         print(bnfisprincipal(bnf, J, 0) != [0]~);\n");
    assert_eq!(lines, ["0 1", "0", "2", "1", "1 2 1", "1", "1", "1"]);
    assert_eq!(obstruction_dimension(&lines), 2, "dim V = r₁ + r₂ + rk₂ Cl for ℚ(√−5)");
}

/// `dim V = r₁ + r₂ + rk₂ Cl` from PARI's printed signature (line 0), unit rank (line 1) and
/// 2-rank of the class group (line 3), after checking that the unit rank is `r₁ + r₂ − 1`.
fn obstruction_dimension(lines: &[String]) -> i64 {
    let signature: Vec<i64> =
        lines[0].split_whitespace().map(|t| t.parse().expect("signature")).collect();
    let (unit_rank, class_two_rank): (i64, i64) =
        (lines[1].parse().expect("unit rank"), lines[3].parse().expect("2-rank"));
    let places = signature[0] + signature[1];
    assert_eq!(unit_rank, places - 1, "unit rank is r₁ + r₂ − 1");
    places + class_two_rank
}
