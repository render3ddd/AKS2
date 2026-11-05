//! Schoolbook multiplication: O(r²) coefficient products. Slow but obviously
//! correct, so it is the oracle every other back-end is tested against.

use rug::Integer;

use super::{CyclicPoly, PolyMul, Ring};

/// The classical O(r²) algorithm, accumulating straight into the r cyclic slots.
#[derive(Debug, Clone, Copy, Default)]
pub struct Schoolbook;

impl PolyMul for Schoolbook {
    fn mul(&self, ring: &Ring, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly {
        ring.check(a);
        ring.check(b);
        let r = ring.r();
        let mut acc = vec![Integer::new(); r];
        for (i, x) in a.coeffs.iter().enumerate().filter(|(_, x)| !x.is_zero()) {
            for (j, y) in b.coeffs.iter().enumerate() {
                let k = if i + j >= r { i + j - r } else { i + j };
                acc[k] += x * y; // mpz_addmul, no temporary
            }
        }
        ring.fold_reduce(acc)
    }

    fn name(&self) -> &'static str {
        "schoolbook"
    }
}
