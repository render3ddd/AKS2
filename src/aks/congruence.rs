//! The congruence tests: (x − a)ⁿ ≡ x^(n mod r) − a in ℤₙ[x]/(xʳ − 1).
//!
//! The same form as the corrected C++ (`CongruenceZnx.h`), so the failing *a*
//! reported for a composite is directly comparable.

use rayon::prelude::*;
use rug::Integer;

use crate::poly::{CyclicPoly, PolyMul, Ring};

/// Precomputed pieces shared by every *a* for one n and r.
#[derive(Debug, Clone)]
pub struct Congruence {
    ring: Ring,
    /// x^(n mod r): xⁿ reduced mod xʳ − 1 needs no power at all.
    x_to_n: CyclicPoly,
}

impl Congruence {
    /// Sets up the ring ℤₙ[x]/(xʳ − 1).
    pub fn new(n: &Integer, r: u64) -> Self {
        let ring = Ring::new(n.clone(), usize::try_from(r).expect("r fits usize"));
        let x_to_n = ring.x_pow(n);
        Congruence { ring, x_to_n }
    }

    /// The ring the tests run in.
    pub fn ring(&self) -> &Ring {
        &self.ring
    }

    /// Whether (x − a)ⁿ ≡ x^(n mod r) − a, using multiplier `m` for the squarings.
    pub fn holds<M: PolyMul + ?Sized>(&self, m: &M, a: u64) -> bool {
        let lhs = self.ring.pow_x_minus_a(m, a, self.ring.n());
        let rhs = self
            .ring
            .sub(&self.x_to_n, &self.ring.monomial(0, Integer::from(a)));
        lhs == rhs
    }

    /// The smallest a in 1..=`a_max` for which the congruence fails, if any.
    pub fn first_failure<M: PolyMul + ?Sized>(&self, m: &M, a_max: u64) -> Option<u64> {
        (1..=a_max).find(|&a| !self.holds(m, a))
    }

    /// [`Congruence::first_failure`] across the current rayon thread pool.
    ///
    /// `find_first` returns the *smallest* failing a, whatever order the
    /// threads finish in, and stops work on larger a once one fails. Each
    /// test owns its polynomials, so nothing is shared mutably: the data race
    /// in the 2021 `CongruenceZnxHPC.h` (REFACTOR.md §1.6) cannot happen.
    pub fn first_failure_par<M: PolyMul + ?Sized>(&self, m: &M, a_max: u64) -> Option<u64> {
        (1..=a_max)
            .into_par_iter()
            .find_first(|&a| !self.holds(m, a))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{Kronecker, Schoolbook};

    #[test]
    fn primes_pass_and_composites_fail() {
        let prime = Integer::from(10007);
        let c = Congruence::new(&prime, 173);
        assert_eq!(c.first_failure(&Kronecker, 20), None);

        // 91 = 7 * 13 with r = 5 (deliberately not the AKS r): some a fails.
        let comp = Integer::from(91);
        let c = Congruence::new(&comp, 5);
        let k = c.first_failure(&Kronecker, 90);
        assert!(k.is_some());
        assert_eq!(c.first_failure(&Schoolbook, 90), k);
    }

    #[test]
    fn parallel_reports_the_smallest_failing_a() {
        // Small deliberately-wrong r values make many a fail; the parallel
        // search must still report the same, smallest one as the serial one.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap();
        for (n, r) in [(91u64, 5u64), (561, 7), (1_022_117, 11), (10007 * 10009, 3)] {
            let c = Congruence::new(&Integer::from(n), r);
            let a_max = 200.min(n - 1);
            let serial = c.first_failure(&Kronecker, a_max);
            for _ in 0..5 {
                assert_eq!(
                    pool.install(|| c.first_failure_par(&Kronecker, a_max)),
                    serial,
                    "n={n} r={r}"
                );
            }
        }
    }
}
