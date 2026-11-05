//! Kronecker (M5), the M2 byte-slot Kronecker and Schoolbook by r, for n of the sizes in RUST_PLAN.md §2.
//!
//!     cargo bench --bench poly
//!
//! Kronecker squaring should grow almost linearly in r (GMP's multiply is
//! close to O(m log m)), not by 4× per doubling as the 2021 code did.

use std::hint::black_box;

use aks2::poly::{Kronecker, KroneckerBytes, PolyMul, Ring, Schoolbook};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use rug::Integer;
use rug::rand::RandState;

fn random_poly(ring: &Ring, rng: &mut RandState<'_>) -> aks2::poly::CyclicPoly {
    let n = ring.n().clone();
    ring.from_coeffs((0..ring.r()).map(|_| Integer::from(n.random_below_ref(rng))))
}

fn squaring(c: &mut Criterion) {
    let mut rng = RandState::new();
    for n in [
        "1000000007",
        "1003026954441971",
        "689960931088884849033689023336009222695077",
    ] {
        let n: Integer = n.parse().unwrap();
        let mut group = c.benchmark_group(format!("square/{}-digit n", n.to_string().len()));
        for r in [256usize, 512, 1024, 2048, 4096, 8192, 19051] {
            let ring = Ring::new(n.clone(), r);
            let a = random_poly(&ring, &mut rng);
            group.bench_with_input(BenchmarkId::new("kronecker", r), &a, |bch, a| {
                bch.iter(|| Kronecker.square(&ring, black_box(a)))
            });
            group.bench_with_input(BenchmarkId::new("kronecker-bytes", r), &a, |bch, a| {
                bch.iter(|| KroneckerBytes.square(&ring, black_box(a)))
            });
            if r <= 1024 {
                group.bench_with_input(BenchmarkId::new("schoolbook", r), &a, |bch, a| {
                    bch.iter(|| Schoolbook.mul(&ring, black_box(a), black_box(a)))
                });
            }
        }
        group.finish();
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(20);
    targets = squaring
}
criterion_main!(benches);
