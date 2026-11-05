//! M2: every multiplier agrees with Schoolbook (RUST_PLAN.md §7 rule 7), on
//! random inputs and on the edge cases that broke the 2021 C++.

use aks2::poly::{CyclicPoly, Kronecker, KroneckerBytes, PolyMul, Ring, Schoolbook};
use proptest::prelude::*;
use rug::Integer;
use rug::integer::Order;

/// n ≥ 2 of up to `max_limbs` 64-bit limbs, biased towards awkward shapes.
fn modulus(max_limbs: usize) -> impl Strategy<Value = Integer> {
    prop_oneof![
        (2u64..1000).prop_map(Integer::from),
        any::<u64>().prop_map(|x| Integer::from(x.max(2))),
        proptest::collection::vec(any::<u64>(), 1..=max_limbs)
            .prop_map(|limbs| Integer::from_digits(&limbs, Order::Lsf).max(Integer::from(2))),
        // 2^k and 2^k ± 1: slot widths land exactly on byte boundaries.
        (1u32..=64 * max_limbs as u32, -1i32..=1)
            .prop_map(|(k, d)| { (Integer::from(1) << k) + d })
            .prop_filter("n >= 2", |n| *n >= 2),
    ]
}

/// One polynomial of the given ring: random, or one of the edge-case shapes.
fn poly(ring: Ring) -> impl Strategy<Value = (Ring, CyclicPoly)> {
    let r = ring.r();
    let n = ring.n().clone();
    let limbs = n.significant_bits().div_ceil(64) as usize + 1;
    let random = proptest::collection::vec(proptest::collection::vec(any::<u64>(), limbs), r)
        .prop_map(|cs| {
            cs.into_iter()
                .map(|l| Integer::from_digits(&l, Order::Lsf))
                .collect::<Vec<_>>()
        });
    let nm1 = Integer::from(&n - 1u32);
    let shape = prop_oneof![
        Just(vec![Integer::new()]),                        // zero
        any::<u64>().prop_map(|c| vec![Integer::from(c)]), // constant
        Just(vec![Integer::from(1); r]),                   // all ones
        Just(vec![nm1.clone(); r]),                        // all n - 1: worst-case slots
        (0..r).prop_map(move |k| {
            // sparse: n - 1 at x^k
            let mut v = vec![Integer::new(); k + 1];
            v[k] = Integer::from(&n - 1u32);
            v
        }),
    ];
    prop_oneof![2 => random, 1 => shape].prop_map(move |cs| {
        let p = ring.from_coeffs(cs);
        (ring.clone(), p)
    })
}

fn ring_and_two(
    max_limbs: usize,
    max_r: usize,
) -> impl Strategy<Value = (Ring, CyclicPoly, CyclicPoly)> {
    (modulus(max_limbs), 1..=max_r)
        .prop_flat_map(|(n, r)| {
            let ring = Ring::new(n, r);
            (poly(ring.clone()), poly(ring))
        })
        .prop_map(|((ring, a), (_, b))| (ring, a, b))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn kronecker_mul_matches_schoolbook((ring, a, b) in ring_and_two(8, 64)) {
        prop_assert_eq!(Kronecker.mul(&ring, &a, &b), Schoolbook.mul(&ring, &a, &b));
    }

    #[test]
    fn kronecker_square_matches_schoolbook((ring, a, _b) in ring_and_two(8, 64)) {
        prop_assert_eq!(Kronecker.square(&ring, &a), Schoolbook.mul(&ring, &a, &a));
    }

    #[test]
    fn kronecker_bytes_matches_schoolbook((ring, a, b) in ring_and_two(8, 64)) {
        prop_assert_eq!(KroneckerBytes.mul(&ring, &a, &b), Schoolbook.mul(&ring, &a, &b));
        prop_assert_eq!(KroneckerBytes.square(&ring, &a), Schoolbook.mul(&ring, &a, &a));
    }

    #[test]
    fn square_into_reuses_any_buffer((ring, a, b) in ring_and_two(8, 64)) {
        // The output buffer's old contents (here: b) must not leak into a².
        let mut out = b;
        Kronecker.square_into(&ring, &a, &mut out);
        prop_assert_eq!(out, Schoolbook.mul(&ring, &a, &a));
    }

    #[test]
    fn pow_x_minus_a_matches_pow((ring, _a, _b) in ring_and_two(4, 32), a in any::<u64>(), e in 0u64..5000) {
        let base = ring.from_coeffs([-Integer::from(a), Integer::from(1)]);
        let e = Integer::from(e);
        prop_assert_eq!(ring.pow_x_minus_a(&Kronecker, a, &e), ring.pow(&Schoolbook, &base, &e));
    }

    #[test]
    fn ring_laws((ring, a, b) in ring_and_two(4, 32)) {
        let m = Kronecker;
        prop_assert_eq!(m.mul(&ring, &a, &b), m.mul(&ring, &b, &a));
        prop_assert_eq!(m.mul(&ring, &a, &ring.one()), a.clone());
        // (a + b)(a - b) = a² - b²
        let lhs = m.mul(&ring, &ring.add(&a, &b), &ring.sub(&a, &b));
        let rhs = ring.sub(&m.square(&ring, &a), &m.square(&ring, &b));
        prop_assert_eq!(lhs, rhs);
    }

    #[test]
    fn pow_matches_repeated_multiplication((ring, a, _b) in ring_and_two(2, 16), e in 0u32..40) {
        let mut want = ring.one();
        for _ in 0..e {
            want = Schoolbook.mul(&ring, &want, &a);
        }
        prop_assert_eq!(ring.pow(&Kronecker, &a, &Integer::from(e)), want);
    }
}

/// The plan's full range (n up to 512 bits, r in the thousands) is too slow for
/// Schoolbook in debug builds; CI runs this in release mode.
#[test]
#[ignore = "release-mode only: cargo test --release -- --include-ignored"]
fn kronecker_matches_schoolbook_large_r() {
    let mut runner = proptest::test_runner::TestRunner::new(ProptestConfig::with_cases(16));
    runner
        .run(&ring_and_two(8, 4_096), |(ring, a, b)| {
            prop_assert_eq!(Kronecker.mul(&ring, &a, &b), Schoolbook.mul(&ring, &a, &b));
            prop_assert_eq!(Kronecker.square(&ring, &a), Schoolbook.mul(&ring, &a, &a));
            prop_assert_eq!(
                KroneckerBytes.square(&ring, &a),
                Schoolbook.mul(&ring, &a, &a)
            );
            Ok(())
        })
        .unwrap();
}
