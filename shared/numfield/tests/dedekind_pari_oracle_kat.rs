//! PARI/GP cross-checks of the Dedekind index criterion and the degree-one prime report.
//!
//! Every test here is `#[ignore]`d: it needs `gp` (GP/PARI) on `PATH`, which the default test run
//! does not require. Run them by hand with
//!
//! ```text
//! cargo test -p shared-numfield --test dedekind_pari_oracle_kat -- --ignored
//! ```
//!
//! A test fails (it does not skip) when `gp` is missing, prints to standard error, or disagrees.
//!
//! **Independence.** The crate side is `dedekind_index_test` and `dedekind_factor`; the PARI side
//! is `gp -q` run on the polynomial and the prime alone.
//!
//! - The index valuation is `v_p([O_K : Z[α]]) = (v_p(poldisc f) − v_p(nfdisc f)) / 2`, because
//!   `disc f = [O_K : Z[α]]² · disc K`. The criterion must report `divides_index` exactly when
//!   `v_p > 0`.
//! - Where `p` does not divide the index, `idealprimedec(nfinit(f), p)` is the decomposition of
//!   `(p)` in `O_K`. Each degree-one prime is matched to a report entry by the residue `r` for
//!   which `idealval(nf, x − r, 𝔭) > 0`, and its ramification index must equal the entry's
//!   multiplicity. The remaining primes (residue degree above one) must account for
//!   `nonlinear_degree` as `Σ e·f`.
//! - Whether or not `p` divides the index, the report claims only facts about `f mod p`: the
//!   listed `r` are exactly the roots of `f̄` with their multiplicities (PARI's `factormod`), and
//!   the multiplicities plus `nonlinear_degree` sum to `deg f`. Nothing about `pO_K` is asserted
//!   there, since the report claims nothing.
//!
//! **Fixtures** are chosen for coverage and each is checked irreducible by `polisirreducible`;
//! every one is tried at `p ∈ {2, 3, 5, 7}`. They span degrees 2–4 and include:
//!
//! - `p` dividing the index, with `f̄` inseparable (`x² − 5`, `x² − 8`, `x⁴ − 5`), with a repeated
//!   and a simple root (`x³ + x² − 2x + 8` at 2, where `2` splits completely in `O_K`), with a
//!   repeated irreducible quadratic factor (`x⁴ + 2x³ + 3x² + 2x + 5` at 2), and at an odd prime
//!   (`x² − 18`, `x³ − 10` at 3);
//! - `p` bad but not dividing the index, so ramified yet maximal (`x² − 2`, `x² − 3`, `x² + 1`,
//!   `x³ − 2`, `x³ − 3x + 1`, `x⁴ + 1`, the fifth cyclotomic quartic at 5), including
//!   `x² − 18` at 2, where `f̄ = x²` is inseparable but `Z[α]` is maximal at 2, and a ramified
//!   prime of residue degree 2 (`x⁴ + 2x³ + 3x² + 2x + 3` at 2);
//! - `p ∤ disc f`: inert, split and partially split primes, among them the monic cubic
//!   `x³ + 15x² + 29x + 8` used by the factoring pipeline tests.

use std::io::Write;
use std::process::{Command, Stdio};

use num_bigint::BigInt;
use shared_numfield::{
    dedekind_factor, dedekind_index_test, is_bad_prime, IntPoly, NumberField, PrimeModulus,
    RootSearchLimits,
};

const GP_MISSING: &str = "PARI/GP (`gp`) must be on PATH to run the oracle tests";

/// A budget far above any root search over the small primes used here.
const AMPLE: RootSearchLimits = RootSearchLimits { max_evaluations: 1_000 };

const PRIMES: [u64; 4] = [2, 3, 5, 7];

// ─── fixtures ────────────────────────────────────────────────────────────────

/// A monic polynomial with its coefficients lowest degree first.
struct Fixture {
    name: &'static str,
    coeffs: &'static [i64],
}

const FIXTURES: &[Fixture] = &[
    Fixture { name: "x^2-2", coeffs: &[-2, 0, 1] },
    Fixture { name: "x^2-3", coeffs: &[-3, 0, 1] },
    Fixture { name: "x^2-5", coeffs: &[-5, 0, 1] },
    Fixture { name: "x^2-17", coeffs: &[-17, 0, 1] },
    Fixture { name: "x^2+1", coeffs: &[1, 0, 1] },
    Fixture { name: "x^2+3", coeffs: &[3, 0, 1] },
    Fixture { name: "x^2-8", coeffs: &[-8, 0, 1] },
    Fixture { name: "x^2-18", coeffs: &[-18, 0, 1] },
    Fixture { name: "x^3-2", coeffs: &[-2, 0, 0, 1] },
    Fixture { name: "x^3-10", coeffs: &[-10, 0, 0, 1] },
    Fixture { name: "x^3-x-1", coeffs: &[-1, -1, 0, 1] },
    Fixture { name: "x^3+15x^2+29x+8", coeffs: &[8, 29, 15, 1] },
    Fixture { name: "x^3+x^2-2x+8", coeffs: &[8, -2, 1, 1] },
    Fixture { name: "x^3-3x+1", coeffs: &[1, -3, 0, 1] },
    Fixture { name: "x^4+1", coeffs: &[1, 0, 0, 0, 1] },
    Fixture { name: "x^4-10x^2+1", coeffs: &[1, 0, -10, 0, 1] },
    Fixture { name: "x^4-2", coeffs: &[-2, 0, 0, 0, 1] },
    Fixture { name: "x^4-5", coeffs: &[-5, 0, 0, 0, 1] },
    Fixture { name: "x^4+x^3+x^2+x+1", coeffs: &[1, 1, 1, 1, 1] },
    Fixture { name: "x^4+2x^3+3x^2+2x+3", coeffs: &[3, 2, 3, 2, 1] },
    Fixture { name: "x^4+2x^3+3x^2+2x+5", coeffs: &[5, 2, 3, 2, 1] },
];

impl Fixture {
    fn poly(&self) -> IntPoly {
        IntPoly::from_coeffs(self.coeffs.iter().map(|&c| BigInt::from(c)).collect())
    }

    fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// The GP definition `f = …` of the fixture polynomial in the variable `x`.
    fn gp_definition(&self) -> String {
        let coefficients: Vec<String> = self.coeffs.iter().map(i64::to_string).collect();
        format!("f = Polrev([{}], 'x);", coefficients.join(","))
    }
}

fn prime(p: u64) -> PrimeModulus {
    PrimeModulus::new(p).expect("small prime")
}

/// Every fixture at every tested prime.
fn cases() -> Vec<(&'static Fixture, u64)> {
    FIXTURES.iter().flat_map(|fx| PRIMES.iter().map(move |&p| (fx, p))).collect()
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

/// What PARI says about `f` at `p`.
struct Oracle {
    /// `v_p([O_K : Z[α]])`.
    index_valuation: i64,
    /// `v_p(disc f)`.
    discriminant_valuation: i64,
    /// `(r, e)` for each degree-one prime of `O_K` above `p`, from `idealprimedec`.
    ideal_degree_one: Vec<(u64, usize)>,
    /// `Σ e·f` over the primes of `O_K` above `p` with residue degree above one.
    ideal_nonlinear_weight: usize,
    /// `(r, e)` for each linear factor `(x − r)^e` of `f mod p`, from `factormod`.
    reduction_roots: Vec<(u64, usize)>,
}

/// Ask PARI about the fixture at `p`, failing unless PARI proves the polynomial irreducible.
///
/// The script is one statement per line: the number field is only built once irreducibility
/// is confirmed, so a reducible fixture fails on its own assertion rather than inside `nfinit`.
fn oracle(fx: &Fixture, p: u64) -> Oracle {
    let script = format!(
        "{def} p = {p};\n\
         print(polisirreducible(f));\n\
         if(!polisirreducible(f), quit);\n\
         nf = nfinit(f);\n\
         print((valuation(poldisc(f),p) - valuation(nfdisc(f),p))/2);\n\
         print(valuation(poldisc(f),p));\n\
         pd = idealprimedec(nf,p);\n\
         for(i=1, #pd, if(pd[i].f==1, my(rs=select(r->idealval(nf,x-r,pd[i])>0, [0..p-1])); \
         print(\"D \", #rs, \" \", if(#rs==1, rs[1], -1), \" \", pd[i].e), \
         print(\"N \", pd[i].e*pd[i].f)));\n\
         fm = factormod(f,p);\n\
         for(i=1, matsize(fm)[1], my(g=lift(fm[i,1])); if(poldegree(g)==1, \
         print(\"R \", (-polcoeff(g,0))%p, \" \", fm[i,2])));\n",
        def = fx.gp_definition(),
    );
    let lines = gp(&script);
    assert_eq!(
        lines.first().map(String::as_str),
        Some("1"),
        "{} must be irreducible over Q: {lines:?}",
        fx.name
    );
    let mut oracle = Oracle {
        index_valuation: parse(&lines[1], fx, p),
        discriminant_valuation: parse(&lines[2], fx, p),
        ideal_degree_one: Vec::new(),
        ideal_nonlinear_weight: 0,
        reduction_roots: Vec::new(),
    };
    for line in &lines[3..] {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.as_slice() {
            ["D", count, r, e] => {
                // A degree-one prime contains exactly one linear form x − r with r in 0..p.
                assert_eq!(*count, "1", "{} at {p}: degree-one prime without one r", fx.name);
                oracle.ideal_degree_one.push((parse(r, fx, p), parse(e, fx, p)));
            }
            ["N", weight] => oracle.ideal_nonlinear_weight += parse::<usize>(weight, fx, p),
            ["R", r, e] => oracle.reduction_roots.push((parse(r, fx, p), parse(e, fx, p))),
            _ => panic!("{} at {p}: unrecognised gp line `{line}`", fx.name),
        }
    }
    oracle.ideal_degree_one.sort_unstable();
    oracle.reduction_roots.sort_unstable();
    oracle
}

fn parse<T: std::str::FromStr>(text: &str, fx: &Fixture, p: u64) -> T {
    text.parse().unwrap_or_else(|_| panic!("{} at {p}: `{text}` is not a number", fx.name))
}

// ─── crate side ──────────────────────────────────────────────────────────────

/// The report's degree-one entries as `(r, multiplicity)`, ascending in `r`, and its
/// `nonlinear_degree`. Also checks that every listed ideal lies over `p`.
fn report_summary(fx: &Fixture, p: u64) -> (Vec<(u64, usize)>, usize) {
    let field = NumberField::new(fx.poly());
    let report = dedekind_factor(&field, prime(p), AMPLE).expect("root search within budget");
    let entries = report
        .degree_one
        .iter()
        .map(|d| {
            assert_eq!(d.ideal.p, BigInt::from(p), "{} at {p}: wrong prime", fx.name);
            let r = u64::try_from(&d.ideal.r).expect("r is a residue");
            (r, d.multiplicity)
        })
        .collect();
    (entries, report.nonlinear_degree)
}

// ─── tests ───────────────────────────────────────────────────────────────────

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_dedekind_fixtures_are_irreducible_and_monic() {
    for fx in FIXTURES {
        assert_eq!(fx.coeffs.last(), Some(&1), "{} is monic", fx.name);
        let lines = gp(&format!("{} print(polisirreducible(f));\n", fx.gp_definition()));
        assert_eq!(lines, ["1"], "{} must be irreducible over Q", fx.name);
    }
}

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_dedekind_index_criterion_matches_the_index_valuation() {
    let (mut dividing, mut bad_but_maximal) = (0, 0);
    for (fx, p) in cases() {
        let pari = oracle(fx, p);
        let test = dedekind_index_test(&fx.poly(), prime(p)).expect("monic, nonconstant");
        assert_eq!(
            test.divides_index(),
            pari.index_valuation > 0,
            "{} at {p}: v_p(index) = {} but obstruction is {:?}",
            fx.name,
            pari.index_valuation,
            test.obstruction()
        );
        // Every prime dividing the index divides disc f; the converse fails.
        assert_eq!(is_bad_prime(&fx.poly(), &BigInt::from(p)), pari.discriminant_valuation > 0);
        assert!(pari.index_valuation == 0 || pari.discriminant_valuation > 0);
        if test.divides_index() {
            dividing += 1;
        } else if pari.discriminant_valuation > 0 {
            bad_but_maximal += 1;
        }
    }
    // Guard the fixture set: both interesting classes must stay represented.
    assert!(dividing >= 6, "only {dividing} cases with p dividing the index");
    assert!(bad_but_maximal >= 10, "only {bad_but_maximal} bad primes that are maximal");
}

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_dedekind_report_lists_the_roots_of_the_reduction() {
    for (fx, p) in cases() {
        let pari = oracle(fx, p);
        let (entries, nonlinear) = report_summary(fx, p);
        assert_eq!(entries, pari.reduction_roots, "{} at {p}: roots of f mod p", fx.name);
        let linear: usize = entries.iter().map(|&(_, e)| e).sum();
        assert_eq!(linear + nonlinear, fx.degree(), "{} at {p}: degrees add up", fx.name);
    }
}

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_dedekind_report_matches_idealprimedec_when_p_does_not_divide_the_index() {
    let mut ramified_degree_one = 0;
    let mut nonlinear_primes = 0;
    for (fx, p) in cases() {
        let pari = oracle(fx, p);
        if pari.index_valuation > 0 {
            continue;
        }
        let (entries, nonlinear) = report_summary(fx, p);
        assert_eq!(entries, pari.ideal_degree_one, "{} at {p}: degree-one primes", fx.name);
        assert_eq!(nonlinear, pari.ideal_nonlinear_weight, "{} at {p}: Σ e·f, f > 1", fx.name);
        ramified_degree_one += entries.iter().filter(|&&(_, e)| e > 1).count();
        nonlinear_primes += usize::from(nonlinear > 0);
    }
    // Guard the fixture set: ramified degree-one primes and nonlinear contributions both occur.
    assert!(ramified_degree_one >= 10, "only {ramified_degree_one} ramified degree-one primes");
    assert!(nonlinear_primes >= 20, "only {nonlinear_primes} cases with nonlinear primes");
}

#[test]
#[ignore = "requires PARI/GP (`gp`) on PATH; run manually with --ignored"]
fn kat_pari_dedekind_report_claims_only_the_reduction_when_p_divides_the_index() {
    let mut dividing = 0;
    for (fx, p) in cases() {
        let pari = oracle(fx, p);
        if pari.index_valuation == 0 {
            continue;
        }
        dividing += 1;
        let field = NumberField::new(fx.poly());
        let report = dedekind_factor(&field, prime(p), AMPLE).expect("root search within budget");
        assert!(report.index.divides_index(), "{} at {p}", fx.name);
        let linear: usize = report.degree_one.iter().map(|d| d.multiplicity).sum();
        assert_eq!(linear + report.nonlinear_degree, fx.degree(), "{} at {p}", fx.name);
        // The listed ideals are (p, α − r) for roots r of f̄; PARI's factorisation of f mod p
        // contains the linear factor x − r with that exponent.
        let listed: Vec<(u64, usize)> = report
            .degree_one
            .iter()
            .map(|d| (u64::try_from(&d.ideal.r).expect("r is a residue"), d.multiplicity))
            .collect();
        assert_eq!(listed, pari.reduction_roots, "{} at {p}", fx.name);
    }
    assert!(dividing >= 6, "only {dividing} cases with p dividing the index");
}
