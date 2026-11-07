//! Known-answer tests for the bounded enumeration of XOR combinations of kernel vectors.
//!
//! A kernel vector is a set of filtered-matrix row indices, and the sum of two kernel vectors is
//! their symmetric difference. The enumeration reads only the row indices it is given. Hand-built
//! vectors here isolate the enumeration logic; they say nothing about any matrix. Only the last
//! test uses vectors computed from a real matrix, and only to confirm that each combination is
//! again a kernel vector of that matrix.
//!
//! **Fixture A** has three independent vectors `v1 = {0,1}`, `v2 = {1,2}`, `v3 = {3}`. Their
//! largest indices 1, 2, 3 are distinct, so reducing them costs no XOR. The seven nonzero
//! combinations, by subset size and then lexicographic order, are (recomputed with an independent
//! Python script, `symmetric_difference` over `set`s)
//!
//! | subset | combination | XORs |
//! |--------|-------------|------|
//! | v1 | {0,1} | 0 |
//! | v2 | {1,2} | 0 |
//! | v3 | {3} | 0 |
//! | v1⊕v2 | {0,1}△{1,2} = {0,2} | 1 |
//! | v1⊕v3 | {0,1,3} | 1 |
//! | v2⊕v3 | {1,2,3} | 1 |
//! | v1⊕v2⊕v3 | {0,2}△{3} = {0,2,3} | 2 |
//!
//! so a complete enumeration costs 0+0+0+1+1+1+2 = 5 XOR operations.

use std::collections::HashSet;

use gnfs::filter::{build_checked_matrix, merge_columns, remove_singletons};
use gnfs::linalg::{
    block_lanczos, search_kernel_combinations, KernelSearch, KernelSearchLimits,
    KernelSearchPhase, KernelSearchResource, KernelSearchStatus, KernelVector, MatrixOperator,
};
use gnfs::polyselect::PolyPair;
use gnfs::sieve::{
    line_sieve, CheckedRelations, FactorBase, LineSieveConfig, RelationChecker,
};
use gnfs::sqrt::SqrtContext;
use num_bigint::BigInt;
use shared_numfield::IntPoly;

fn kv(rows: &[usize]) -> KernelVector {
    KernelVector { row_indices: rows.to_vec() }
}

fn kvs(vectors: &[&[usize]]) -> Vec<KernelVector> {
    vectors.iter().map(|rows| kv(rows)).collect()
}

fn rows(search: &KernelSearch) -> Vec<Vec<usize>> {
    search.candidates.iter().map(|c| c.row_indices.clone()).collect()
}

fn limits(max_xor_operations: u64, max_candidates: u64) -> KernelSearchLimits {
    KernelSearchLimits { max_xor_operations, max_candidates }
}

const AMPLE: KernelSearchLimits = KernelSearchLimits { max_xor_operations: 1 << 40, max_candidates: 1 << 40 };

fn fixture_a() -> Vec<KernelVector> {
    kvs(&[&[0, 1], &[1, 2], &[3]])
}

/// Fixture A's combinations in the frozen order.
fn fixture_a_sequence() -> Vec<Vec<usize>> {
    vec![
        vec![0, 1],
        vec![1, 2],
        vec![3],
        vec![0, 2],
        vec![0, 1, 3],
        vec![1, 2, 3],
        vec![0, 2, 3],
    ]
}

fn limited(
    resource: KernelSearchResource,
    phase: KernelSearchPhase,
    xor_operations: u64,
    candidates: u64,
) -> KernelSearchStatus {
    KernelSearchStatus::ResourceLimited { resource, phase, xor_operations, candidates }
}

#[test]
fn kat_combinations_include_sums_beyond_basis_vectors() {
    let search = search_kernel_combinations(&fixture_a(), AMPLE);
    assert_eq!(search.status, KernelSearchStatus::Complete);
    assert_eq!(search.rank, 3);
    assert_eq!(rows(&search), fixture_a_sequence());
    // v1 ⊕ v2 = {0,1} △ {1,2} = {0,2} is its own candidate, absent from the inputs.
    assert!(rows(&search).contains(&vec![0, 2]));
    assert!(!fixture_a().contains(&kv(&[0, 2])));
    // 2^3 − 1 distinct nonzero combinations.
    assert_eq!(search.candidates.len(), 7);
}

#[test]
fn kat_dependent_and_duplicate_inputs_reduce_to_independent_set() {
    // Inputs in order: a = {0,1}; a again; b = {1,2}; c = {0,2} = a ⊕ b; the zero vector; {4,4},
    // which cancels to zero; d = {5}. Only a, b, d are independent, so the rank is 3 and a complete
    // search has 2^3 − 1 = 7 candidates, with a ⊕ b = {0,2} emitted once (it equals c).
    //
    // Reduction XORs (echelon keyed by the largest index): the second a meets pivot a (top 1):
    // one XOR, to the empty set. c = {0,2} meets pivot b (top 2): {0,2}△{1,2} = {0,1}, then pivot a
    // (top 1): {0,1}△{0,1} = {}: two XORs. Reduction costs 3; enumeration costs 5 as in Fixture A.
    let inputs = kvs(&[&[0, 1], &[0, 1], &[1, 2], &[0, 2], &[], &[4, 4], &[5]]);
    let search = search_kernel_combinations(&inputs, AMPLE);
    assert_eq!(search.status, KernelSearchStatus::Complete);
    assert_eq!(search.rank, 3);
    assert_eq!(
        rows(&search),
        vec![
            vec![0, 1],
            vec![1, 2],
            vec![5],
            vec![0, 2],
            vec![0, 1, 5],
            vec![1, 2, 5],
            vec![0, 2, 5],
        ]
    );
    // The work budget accounts for the reduction: 3 + 5 = 8 is exactly enough, 7 is not.
    assert!(search_kernel_combinations(&inputs, limits(8, 100)).status.is_complete());
    assert_eq!(
        search_kernel_combinations(&inputs, limits(7, 100)).status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Enumeration, 7, 6)
    );
}

#[test]
fn kat_large_dependent_input_set_has_rank_based_count() {
    // e_i = {i} for i in 0..10 are independent. The extra vectors are sums of those:
    // {0,1} = e0 ⊕ e1, {2,3,4} = e2 ⊕ e3 ⊕ e4, and {0,1,2,3,4,5,6,7,8,9} = Σ e_i. The rank is 10,
    // so a complete search has 2^10 − 1 = 1023 candidates, all distinct, sorted, and nonzero.
    let mut inputs: Vec<KernelVector> = (0..10).map(|i| kv(&[i])).collect();
    inputs.push(kv(&[0, 1]));
    inputs.push(kv(&[2, 3, 4]));
    inputs.push(kv(&(0..10).collect::<Vec<_>>()));
    let search = search_kernel_combinations(&inputs, AMPLE);
    assert_eq!(search.status, KernelSearchStatus::Complete);
    assert_eq!(search.rank, 10);
    assert_eq!(search.candidates.len(), 1023);
    let distinct: HashSet<&Vec<usize>> = search.candidates.iter().map(|c| &c.row_indices).collect();
    assert_eq!(distinct.len(), 1023, "no combination is emitted twice");
    for candidate in &search.candidates {
        assert!(!candidate.row_indices.is_empty());
        assert!(candidate.row_indices.windows(2).all(|w| w[0] < w[1]), "{candidate:?}");
    }
}

#[test]
fn kat_zero_cancellation_never_emitted() {
    // v ⊕ v is empty: a duplicated input contributes one candidate, not the empty sum.
    let doubled = search_kernel_combinations(&kvs(&[&[0, 1], &[0, 1]]), AMPLE);
    assert_eq!(rows(&doubled), vec![vec![0, 1]]);
    assert_eq!(doubled.rank, 1);
    assert!(doubled.status.is_complete());

    // a ⊕ b ⊕ c = {0,1} △ {1,2} △ {0,2} = {} (each index appears twice), so c is dependent and the
    // empty triple sum is never a candidate. The rank is 2: 3 candidates, not 7.
    let cycle = search_kernel_combinations(&kvs(&[&[0, 1], &[1, 2], &[0, 2]]), AMPLE);
    assert_eq!(rows(&cycle), vec![vec![0, 1], vec![1, 2], vec![0, 2]]);
    assert_eq!(cycle.rank, 2);
    assert!(cycle.status.is_complete());

    // No vectors, and vectors that are zero or cancel to zero: nothing to emit, and complete,
    // because 2^0 − 1 = 0.
    for inputs in [kvs(&[]), kvs(&[&[]]), kvs(&[&[7, 7], &[]])] {
        let search = search_kernel_combinations(&inputs, limits(0, 0));
        assert!(search.candidates.is_empty());
        assert_eq!(search.rank, 0);
        assert!(search.status.is_complete());
    }
}

#[test]
fn kat_malformed_input_rows_do_not_panic() {
    // Unsorted and repeated indices are read as GF(2) sets: an index listed an even number of times
    // cancels. {3,1,2,1} → {2,3}; {9,9,9} → {9}; {2,2,5} → {5}.
    let inputs = kvs(&[&[3, 1, 2, 1], &[9, 9, 9], &[2, 2, 5]]);
    let search = search_kernel_combinations(&inputs, AMPLE);
    assert_eq!(
        rows(&search),
        vec![vec![2, 3], vec![9], vec![5], vec![2, 3, 9], vec![2, 3, 5], vec![5, 9], vec![2, 3, 5, 9]]
    );
    assert!(search.status.is_complete());

    // Out-of-range values are labels, carried through: this layer has no matrix to range-check.
    let huge = search_kernel_combinations(&kvs(&[&[usize::MAX, 0], &[usize::MAX]]), AMPLE);
    assert_eq!(rows(&huge), vec![vec![0, usize::MAX], vec![usize::MAX], vec![0]]);
    assert!(huge.status.is_complete());
}

#[test]
fn kat_order_is_deterministic_and_pinned() {
    let first = search_kernel_combinations(&fixture_a(), AMPLE);
    let second = search_kernel_combinations(&fixture_a(), AMPLE);
    assert_eq!(first, second);
    assert_eq!(rows(&first), fixture_a_sequence());

    // Four independent singletons: subsets by size, then lexicographic over input position.
    // Size 1: 0 1 2 3; size 2: 01 02 03 12 13 23; size 3: 012 013 023 123; size 4: 0123.
    let singles: Vec<KernelVector> = (0..4).map(|i| kv(&[i])).collect();
    let search = search_kernel_combinations(&singles, AMPLE);
    let expected: Vec<Vec<usize>> = [
        &[0][..],
        &[1],
        &[2],
        &[3],
        &[0, 1],
        &[0, 2],
        &[0, 3],
        &[1, 2],
        &[1, 3],
        &[2, 3],
        &[0, 1, 2],
        &[0, 1, 3],
        &[0, 2, 3],
        &[1, 2, 3],
        &[0, 1, 2, 3],
    ]
    .iter()
    .map(|s| s.to_vec())
    .collect();
    assert_eq!(rows(&search), expected);

    // The order follows input order of the independent vectors, not their index values.
    let swapped = search_kernel_combinations(&kvs(&[&[3], &[1, 2], &[0, 1]]), AMPLE);
    assert_eq!(&rows(&swapped)[..3], &[vec![3], vec![1, 2], vec![0, 1]]);
    assert_eq!(rows(&swapped)[3], vec![1, 2, 3]);
}

#[test]
fn kat_xor_budget_limits_enumeration() {
    let full = rows(&search_kernel_combinations(&fixture_a(), AMPLE));

    // Exactly sufficient: 5 XORs are needed and 5 are allowed.
    let exact = search_kernel_combinations(&fixture_a(), limits(5, 100));
    assert_eq!(exact.status, KernelSearchStatus::Complete);
    assert_eq!(rows(&exact), full);

    // One short: the three pairs cost 3, then the triple's first XOR is the 4th, its second would
    // be the 5th. Six candidates are out, four XORs done.
    let short = search_kernel_combinations(&fixture_a(), limits(4, 100));
    assert_eq!(
        short.status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Enumeration, 4, 6)
    );
    assert!(!short.status.is_complete());
    // A smaller budget's output is a prefix of a larger budget's: retrying reproduces it.
    assert_eq!(rows(&short), full[..6].to_vec());

    // Zero budget: single vectors need no XOR, the first pair does.
    let zero = search_kernel_combinations(&fixture_a(), limits(0, 100));
    assert_eq!(
        zero.status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Enumeration, 0, 3)
    );
    assert_eq!(rows(&zero), full[..3].to_vec());
}

#[test]
fn kat_xor_budget_can_end_during_reduction() {
    // Inputs a, a, b, c = a ⊕ b, d = {5}. The second a costs 1 XOR; c costs 2 (see the dependent
    // input test). With 2 XORs allowed, c's second XOR is refused: the reduction is cut short with
    // a and b found, nothing emitted, and the status says the reduction was running.
    let inputs = kvs(&[&[0, 1], &[0, 1], &[1, 2], &[0, 2], &[5]]);
    let search = search_kernel_combinations(&inputs, limits(2, 100));
    assert_eq!(
        search.status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Reduction, 2, 0)
    );
    assert!(search.candidates.is_empty());
    assert_eq!(search.rank, 2);
}

#[test]
fn kat_candidate_budget_limits_output() {
    let full = rows(&search_kernel_combinations(&fixture_a(), AMPLE));

    // Exactly the 7 combinations: complete, because nothing remains to refuse.
    let exact = search_kernel_combinations(&fixture_a(), limits(1000, 7));
    assert_eq!(exact.status, KernelSearchStatus::Complete);
    assert_eq!(rows(&exact), full);

    // 6 allowed: the seventh is refused. Work done: the 3 pair XORs.
    let short = search_kernel_combinations(&fixture_a(), limits(1000, 6));
    assert_eq!(
        short.status,
        limited(KernelSearchResource::Candidates, KernelSearchPhase::Enumeration, 3, 6)
    );
    assert_eq!(rows(&short), full[..6].to_vec());

    // Zero output allowed with a nonempty search space: limited, not complete.
    let none = search_kernel_combinations(&fixture_a(), limits(1000, 0));
    assert_eq!(
        none.status,
        limited(KernelSearchResource::Candidates, KernelSearchPhase::Enumeration, 0, 0)
    );
    assert!(none.candidates.is_empty());

    // When both budgets are spent at one step, the output budget is the one reported: after 3
    // singles with 0 XORs allowed and 3 candidates allowed.
    let both = search_kernel_combinations(&fixture_a(), limits(0, 3));
    assert_eq!(
        both.status,
        limited(KernelSearchResource::Candidates, KernelSearchPhase::Enumeration, 0, 3)
    );
}

#[test]
fn kat_large_input_set_with_small_budget_does_not_scan_the_subset_space() {
    // 64 independent singletons have 2^64 − 1 combinations; the call must return at once with a
    // limited status. All 64 singles come first, then pairs in lexicographic order.
    let singles: Vec<KernelVector> = (0..64).map(|i| kv(&[i])).collect();

    // C(64,2) = 2016 pairs. With 100 candidates allowed: 64 singles + 36 pairs = 100, and each pair
    // costs one XOR, so 36 XORs were done.
    let capped = search_kernel_combinations(&singles, limits(u64::MAX, 100));
    assert_eq!(
        capped.status,
        limited(KernelSearchResource::Candidates, KernelSearchPhase::Enumeration, 36, 100)
    );
    assert_eq!(capped.rank, 64);
    assert_eq!(rows(&capped)[63], vec![63]);
    assert_eq!(rows(&capped)[64], vec![0, 1]);
    // Pairs start {0,1},{0,2},…,{0,63},{1,2},…: the 36th pair is {0,36}, at position 64 + 35 = 99.
    assert_eq!(rows(&capped)[99], vec![0, 36]);

    // The same space with a tiny XOR budget: 64 singles, 10 pairs, then the 11th pair is refused.
    let starved = search_kernel_combinations(&singles, limits(10, u64::MAX));
    assert_eq!(
        starved.status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Enumeration, 10, 74)
    );

    // More vectors than any integer width could index as a subset mask still returns promptly:
    // 300 singles, then 50 pairs (one XOR each), then the 51st pair is refused.
    let many: Vec<KernelVector> = (0..300).map(|i| kv(&[i])).collect();
    let wide = search_kernel_combinations(&many, limits(50, 400));
    assert_eq!(wide.rank, 300);
    assert_eq!(
        wide.status,
        limited(KernelSearchResource::XorOperations, KernelSearchPhase::Enumeration, 50, 350)
    );
}

#[test]
fn kat_status_is_complete_only_for_complete() {
    assert!(KernelSearchStatus::Complete.is_complete());
    for resource in [KernelSearchResource::XorOperations, KernelSearchResource::Candidates] {
        for phase in [KernelSearchPhase::Reduction, KernelSearchPhase::Enumeration] {
            assert!(!limited(resource, phase, 0, 0).is_complete());
        }
    }
}

// ─── computed kernel vectors ─────────────────────────────────────────────────

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn poly(coeffs: &[i64]) -> IntPoly {
    IntPoly::from_coeffs(coeffs.iter().map(|&c| big(c)).collect())
}

#[test]
fn kat_combinations_of_computed_vectors_are_kernel_vectors() {
    // The pinned line-sieve matrix: x³ − x − 1, m = 2, N = 5 (f(2) = 5), FB(30, 30),
    // LineSieveConfig::new(10, 3), after singleton removal and column merging. Lanczos with a fixed
    // seed supplies the vectors; the vectors that verify are the input, and every combination of
    // them must again be a left-nullspace vector, because the kernel is closed under XOR.
    let f = poly(&[-1, -1, 0, 1]);
    let pair = PolyPair::new(f.clone(), poly(&[-2, 1]), big(2), big(5));
    pair.verify().expect("f(2) = 5 = N");
    let fb = FactorBase::new(&f, 30, 30);
    let collected = line_sieve(&pair, &fb, &LineSieveConfig::new(10, 3));
    let ctx = SqrtContext::from_poly_pair(&pair).unwrap();
    let checker = RelationChecker::new(&ctx, &fb).unwrap();
    let relations = CheckedRelations::new(&checker, &collected).unwrap();
    let (matrix, _layout) = build_checked_matrix(&relations, 0).unwrap();
    let matrix = merge_columns(remove_singletons(matrix));

    let operator = MatrixOperator::new(&matrix);
    let computed: Vec<KernelVector> =
        block_lanczos(&operator, 42).into_iter().filter(|v| v.verify(&matrix)).collect();
    assert!(!computed.is_empty(), "positive control: the pinned matrix has a nontrivial kernel");

    let search = search_kernel_combinations(&computed, AMPLE);
    assert!(search.status.is_complete());
    // Non-vacuous: at least two independent vectors, so some candidate is a genuine sum.
    assert!(search.rank >= 2, "rank {}", search.rank);
    assert_eq!(search.candidates.len(), (1usize << search.rank) - 1);
    let distinct: HashSet<&Vec<usize>> = search.candidates.iter().map(|c| &c.row_indices).collect();
    assert_eq!(distinct.len(), search.candidates.len());
    for candidate in &search.candidates {
        assert!(!candidate.is_empty());
        assert!(candidate.verify(&matrix), "{:?} is not in the kernel", candidate.row_indices);
    }
    // Determinism on computed input.
    assert_eq!(search, search_kernel_combinations(&computed, AMPLE));
    // Every independent computed vector is emitted as itself, before any sum.
    for (position, vector) in search.candidates.iter().take(search.rank).enumerate() {
        assert!(computed.contains(vector), "candidate {position} is a computed vector");
    }
}
