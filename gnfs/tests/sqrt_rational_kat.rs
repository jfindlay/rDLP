//! Known-answer tests for the rational square root of a checked dependency.
//!
//! The rational root of a dependency is the exact integer `x₀` with `x₀² = P`, where
//! `P = ∏(a − b·m)` over the lifted relations, reduced modulo `N`. `P > 0` and `P` a perfect
//! square are checked by `rational_sqrt` itself.
//!
//! **Fixtures.** Every relation set is built through the checked pipeline (factor base, checked
//! relations, `build_checked_matrix`, `DependencyContext`), and each dependency is a validated
//! kernel selection, so the sign and rational-prime columns already force `P` to be a positive
//! square. A non-square `P` is therefore not reachable through a validated dependency; the
//! non-square, zero, and negative branches are exercised on the exact-root routine itself by unit
//! tests beside `rational_sqrt`.
//!
//! - **T** (non-monic): `f = 3x² + x − 1`, `m = 3`, `N = 29`, factor base `FB(7, 17)`; relations
//!   (4,1), (7,1), (2,1), (−3,2), (−9,2), (5,2) as rows 0–5. Rational norms `a − 3b`: 1, 4, −1,
//!   −9, −15, −1. Its kernel is spanned by rows {1,2,3} and {2,5} (the third nonzero vector is
//!   {1,3,5}).
//! - **S** (seam): `f = x² − 2`, `m = 23`, `N = 527 = 17·31` (`f(23) = 527`), relation (3,−2):
//!   `a − 23b = 3 + 46 = 49`.
//! - **R** (eight relations on `x² − 2`, `m = 23`, `N = 527`, `FB(13, 13)`): (10,−1), (−3,2),
//!   (−2,−1), (−2,1), (4,3), (2,−1), (10,1), (1,0) as rows 0–7. Rational norms `a − 23b`: 33,
//!   −49, 21, −25, −65, 25, −13, 1.
//!
//! Hand values were recomputed with an independent exact-arithmetic script (Python `fractions`,
//! no Rust), and the arithmetic is written beside each use.

use gnfs::filter::build_checked_matrix;
use gnfs::linalg::{
    KernelSearchLimits, KernelSearchStatus, KernelVector, RootCharacters, search_kernel_combinations,
};
use gnfs::sieve::{
    CheckedRelations, ExponentVector, FactorBase, Relation, RelationChecker,
};
use gnfs::sqrt::{CheckedDependency, DependencyContext, RationalSqrtOutcome, SqrtContext};
use gnfs::{RationalRoot, rational_sqrt};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use shared_numfield::IntPoly;

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

/// Divide every power of `p` out of `n` and return the exponent.
fn strip(n: &mut BigInt, p: u64) -> u32 {
    let p = BigInt::from(p);
    let mut e = 0;
    while !n.is_zero() && (&*n % &p).is_zero() {
        *n /= &p;
        e += 1;
    }
    e
}

/// The relation `(a, b)` with exponent vectors read off its two norms by trial division over the
/// factor base; the checker re-verifies every norm, so a wrong factorization cannot pass.
fn relation_of(fb: &FactorBase, f: &IntPoly, m: i64, a: i64, b: i64) -> Relation {
    let (a_big, b_big) = (big(a), big(b));
    let mut rational_norm = &a_big - &b_big * big(m);
    let rational_sign = rational_norm.is_negative();
    rational_norm = rational_norm.abs();
    let rational: Vec<(usize, u32)> = fb
        .rational_primes
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| {
            let e = strip(&mut rational_norm, p);
            (e > 0).then_some((i, e))
        })
        .collect();
    assert_eq!(rational_norm, big(1), "the rational norm of ({a},{b}) is smooth over the base");

    let d = f.coeffs.len() - 1;
    let mut algebraic_norm = BigInt::zero();
    for (i, c) in f.coeffs.iter().enumerate() {
        algebraic_norm += c * a_big.pow(i as u32) * b_big.pow((d - i) as u32);
    }
    algebraic_norm = algebraic_norm.abs();
    let algebraic: Vec<(usize, u32)> = fb
        .algebraic_ideals
        .iter()
        .filter_map(|ideal| {
            let label = (&a_big - &b_big * big(ideal.r as i64)) % BigInt::from(ideal.p);
            if !label.is_zero() {
                return None;
            }
            let e = strip(&mut algebraic_norm, ideal.p);
            (e > 0).then_some((ideal.index, e))
        })
        .collect();
    assert_eq!(algebraic_norm, big(1), "the algebraic norm of ({a},{b}) is smooth over the base");

    Relation {
        a: a_big,
        b: b_big,
        rational_exponents: ExponentVector { entries: rational },
        algebraic_exponents: ExponentVector { entries: algebraic },
        rational_sign,
    }
}

/// A dependency context over genuine checked relations, with no root-character columns.
fn dependency_context(
    f: &[i64],
    m: i64,
    n: i64,
    bounds: (u64, u64),
    pairs: &[(i64, i64)],
) -> DependencyContext {
    let f = poly(f);
    let ctx = SqrtContext::new(&f, &big(m), &big(n)).expect("valid context");
    let fb = FactorBase::new(&f, bounds.0, bounds.1);
    let checker = RelationChecker::new(&ctx, &fb).expect("valid factor base");
    let relations: Vec<Relation> =
        pairs.iter().map(|&(a, b)| relation_of(&fb, &f, m, a, b)).collect();
    let checked = CheckedRelations::new(&checker, &relations).expect("valid relations");
    let (matrix, layout) = build_checked_matrix(&checked, 0).expect("column count fits");
    let characters = RootCharacters::new(checked.checker(), &[]).expect("no characters");
    DependencyContext::new(checked, matrix, layout, characters).expect("genuine parts")
}

fn context_t() -> DependencyContext {
    dependency_context(
        &[-1, 1, 3],
        3,
        29,
        (7, 17),
        &[(4, 1), (7, 1), (2, 1), (-3, 2), (-9, 2), (5, 2)],
    )
}

fn context_s() -> DependencyContext {
    dependency_context(&[-2, 0, 1], 23, 527, (13, 13), &[(3, -2)])
}

fn context_r() -> DependencyContext {
    dependency_context(
        &[-2, 0, 1],
        23,
        527,
        (13, 13),
        &[(10, -1), (-3, 2), (-2, -1), (-2, 1), (4, 3), (2, -1), (10, 1), (1, 0)],
    )
}

fn validated(context: &DependencyContext, rows: &[usize]) -> CheckedDependency {
    context.validate(rows).unwrap_or_else(|e| panic!("{rows:?} is not a dependency: {e}"))
}

fn root_of(dependency: &CheckedDependency) -> RationalRoot {
    match rational_sqrt(dependency) {
        RationalSqrtOutcome::Root(root) => root,
        other => panic!("expected a rational root, got {other:?}"),
    }
}

#[test]
fn kat_rational_root_is_the_exact_square_root_of_the_product() {
    // T rows {2,5}: P = (2 − 3)(5 − 6) = (−1)(−1) = 1, so x₀ = 1.
    let t = context_t();
    let root = root_of(&validated(&t, &[2, 5]));
    assert_eq!(root.product(), &big(1));
    assert_eq!(root.root(), &big(1));
    assert_eq!(root.modular_image(), &big(1));

    // T rows {1,2,3}: P = (7−3)(2−3)(−3−6) = 4·(−1)·(−9) = 36, so x₀ = 6 and 6 < 29.
    let root = root_of(&validated(&t, &[1, 2, 3]));
    assert_eq!((root.product(), root.root(), root.modular_image()), (&big(36), &big(6), &big(6)));

    // T rows {1,3,5}, the third nonzero kernel vector: P = 4·(−9)·(−1) = 36 as well.
    let root = root_of(&validated(&t, &[1, 3, 5]));
    assert_eq!((root.product(), root.root()), (&big(36), &big(6)));

    // S: P = 3 − (−2)·23 = 49 = 7², so x₀ = 7.
    let root = root_of(&validated(&context_s(), &[0]));
    assert_eq!((root.product(), root.root(), root.modular_image()), (&big(49), &big(7), &big(7)));

    // R rows {1,3,5}: P = (−49)(−25)(25) = 30625 = 175², and 175 < 527.
    let r = context_r();
    let root = root_of(&validated(&r, &[1, 3, 5]));
    assert_eq!(root.product(), &big(30_625));
    assert_eq!(root.root(), &big(175));
    assert_eq!(root.modular_image(), &big(175));

    // R row {7} is (1,0): P = 1, a unit product.
    let root = root_of(&validated(&r, &[7]));
    assert_eq!((root.product(), root.root()), (&big(1), &big(1)));
}

#[test]
fn kat_rational_root_binds_to_the_dependency_product() {
    // The root is computed from the product the dependency carries, so the same P is available
    // to the algebraic side: x₀² = P exactly, x₀ > 0, and the image is x₀ mod N in [0, N).
    let r = context_r();
    for rows in [&[1usize, 3][..], &[5], &[1, 3, 5], &[7]] {
        let dependency = validated(&r, rows);
        let root = root_of(&dependency);
        assert_eq!(root.product(), dependency.rational_product());
        assert_eq!(&(root.root() * root.root()), root.product());
        assert!(root.root().is_positive());
        let n = dependency.context().modulus();
        assert_eq!(root.modular_image(), &(root.root() % n));
        assert!(!root.modular_image().is_negative() && root.modular_image() < n);
    }
}

#[test]
fn kat_validated_dependencies_never_give_a_rational_nonsquare() {
    // Positive control for the module note: over every nonzero kernel combination of T (rank 2)
    // and R (rank 3), the rational root exists. T's kernel is spanned by {1,2,3} and {2,5}; R's
    // by {1,3}, {5} and {7}.
    let limits = KernelSearchLimits { max_xor_operations: 1 << 20, max_candidates: 1 << 20 };
    for (context, basis, expected) in [
        (context_t(), vec![vec![1, 2, 3], vec![2, 5]], 3usize),
        (context_r(), vec![vec![1, 3], vec![5], vec![7]], 7),
    ] {
        let vectors: Vec<KernelVector> =
            basis.into_iter().map(|row_indices| KernelVector { row_indices }).collect();
        for vector in &vectors {
            assert!(vector.verify(context.matrix()), "{vector:?} is a kernel vector");
        }
        let search = search_kernel_combinations(&vectors, limits);
        assert_eq!(search.status, KernelSearchStatus::Complete);
        assert_eq!(search.candidates.len(), expected);
        for candidate in &search.candidates {
            let dependency = context.validate_kernel_vector(candidate).unwrap();
            let root = root_of(&dependency);
            assert_eq!(&(root.root() * root.root()), dependency.rational_product());
        }
    }
}
