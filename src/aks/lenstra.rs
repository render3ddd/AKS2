//! The whole test: Lenstra's (2002) version of AKS, as in the corrected C++.
//!
//! 1. If n is a perfect power, COMPOSITE.
//! 2. Find the smallest r with ord_r(n) > ⌊log₂ n⌋², returning COMPOSITE on
//!    any 1 < gcd(r′, n) < n on the way, and PRIME if no such r is below n.
//! 3. For 1 ≤ a ≤ bound, if (x − a)ⁿ ≢ x^(n mod r) − a, COMPOSITE.
//! 4. PRIME (certified only when the bound is the theorem's).

use std::fmt;
use std::time::{Duration, Instant};

use rug::Integer;

use super::congruence::Congruence;
use super::params::{ARange, Params, RSearch, a_bound, find_r};
use crate::num::perfect_power;
use crate::poly::{Kronecker, PolyMul, Schoolbook};

/// Which [`PolyMul`] back-end runs the squarings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Multiplier {
    /// Kronecker substitution (binary segmentation): the default.
    #[default]
    Kronecker,
    /// O(r²) schoolbook multiplication: for testing.
    Schoolbook,
}

impl Multiplier {
    /// Every back-end, for tests and benchmarks.
    pub const ALL: [Multiplier; 2] = [Multiplier::Kronecker, Multiplier::Schoolbook];

    /// The back-end itself.
    pub fn get(self) -> &'static dyn PolyMul {
        match self {
            Multiplier::Kronecker => &Kronecker,
            Multiplier::Schoolbook => &Schoolbook,
        }
    }

    /// Short name, as used on the command line.
    pub fn name(self) -> &'static str {
        self.get().name()
    }
}

impl std::str::FromStr for Multiplier {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        Multiplier::ALL
            .into_iter()
            .find(|m| m.name() == s)
            .ok_or_else(|| format!("unknown multiplier {s:?} (expected kronecker or schoolbook)"))
    }
}

/// How to run the test.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Options {
    /// Range of *a* for the congruence tests.
    pub range: ARange,
    /// Polynomial multiplier.
    pub multiplier: Multiplier,
    /// Worker threads for the congruence tests: `None` uses rayon's global
    /// pool (one thread per core, or `RAYON_NUM_THREADS`); `Some(1)` runs
    /// serially on the calling thread.
    pub threads: Option<usize>,
}

/// Why the test stopped where it did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Reason {
    /// n = 2 or 3.
    Small,
    /// n = base^exp with exp ≥ 2.
    PerfectPower {
        /// The smallest base.
        base: Integer,
        /// The largest exponent.
        exp: u32,
    },
    /// 1 < gcd(r′, n) < n for some r′ ≤ r.
    Divisor(u64),
    /// No r below n with ord_r(n) > ⌊log₂ n⌋², and n has no factor below n.
    NoRBelowN,
    /// All congruences 1 ≤ a ≤ bound hold.
    CongruencesHold,
    /// The congruence for this *a* fails.
    CongruenceFails(u64),
}

/// Everything a run found out, for printing, logging and tests.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trace {
    /// The input.
    pub n: Integer,
    /// The verdict.
    pub prime: bool,
    /// Why.
    pub reason: Reason,
    /// The options used.
    pub options: Options,
    /// r, ord_r(n), φ(r), λ(r), when the run got that far.
    pub params: Option<Params>,
    /// The upper bound on *a*, when the run got that far.
    pub a_max: Option<u64>,
    /// Threads the congruence tests ran on (1 if the run stopped before them).
    pub threads: usize,
    /// Wall-clock time for the whole test.
    pub elapsed: Duration,
}

impl Trace {
    /// Whether a PRIME verdict is backed by the theorem. `false` for composites.
    pub fn certified(&self) -> bool {
        self.prime && (self.reason != Reason::CongruencesHold || self.options.range.is_certified())
    }
}

/// n < 2: neither prime nor composite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NotApplicable(pub Integer);

impl fmt::Display for NotApplicable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} is neither prime nor composite", self.0)
    }
}

impl std::error::Error for NotApplicable {}

/// Runs the AKS test on n ≥ 2.
pub fn aks(n: &Integer, options: Options) -> Result<Trace, NotApplicable> {
    if *n < 2 {
        return Err(NotApplicable(n.clone()));
    }
    let start = Instant::now();
    let done = |prime, reason, params, a_max, threads| Trace {
        n: n.clone(),
        prime,
        reason,
        options,
        params,
        a_max,
        threads,
        elapsed: start.elapsed(),
    };

    if *n <= 3 {
        return Ok(done(true, Reason::Small, None, None, 1));
    }
    if let Some((base, exp)) = perfect_power(n) {
        return Ok(done(
            false,
            Reason::PerfectPower { base, exp },
            None,
            None,
            1,
        ));
    }
    let params = match find_r(n) {
        RSearch::Divisor(d) => return Ok(done(false, Reason::Divisor(d), None, None, 1)),
        RSearch::NoRBelowN => return Ok(done(true, Reason::NoRBelowN, None, None, 1)),
        RSearch::Found(p) => p,
    };
    let a_max = a_bound(options.range, n, &params);
    let congruence = Congruence::new(n, params.r);
    let m = options.multiplier.get();
    let (failure, threads) = match options.threads {
        Some(1) => (congruence.first_failure(m, a_max), 1),
        Some(t) => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .expect("building the thread pool");
            (pool.install(|| congruence.first_failure_par(m, a_max)), t)
        }
        None => (
            congruence.first_failure_par(m, a_max),
            rayon::current_num_threads(),
        ),
    };
    Ok(match failure {
        None => done(
            true,
            Reason::CongruencesHold,
            Some(params),
            Some(a_max),
            threads,
        ),
        Some(a) => done(
            false,
            Reason::CongruenceFails(a),
            Some(params),
            Some(a_max),
            threads,
        ),
    })
}

/// Convenience: the verdict alone, with default options. Panics for n < 2.
pub fn is_prime(n: &Integer) -> bool {
    aks(n, Options::default()).expect("n >= 2").prime
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(n: u64) -> Trace {
        aks(&Integer::from(n), Options::default()).unwrap()
    }

    #[test]
    fn reasons() {
        assert_eq!(run(2).reason, Reason::Small);
        assert_eq!(
            run(64).reason,
            Reason::PerfectPower {
                base: Integer::from(2),
                exp: 6
            }
        );
        assert_eq!(run(91).reason, Reason::Divisor(7));
        assert_eq!(run(5).reason, Reason::NoRBelowN);
        let t = run(10007);
        assert_eq!(t.reason, Reason::CongruencesHold);
        assert_eq!(t.params.unwrap().r, 173);
        assert_eq!(t.a_max, Some(174));
        assert!(t.prime && t.certified());
        assert!(aks(&Integer::from(1), Options::default()).is_err());
    }

    #[test]
    fn thread_counts_agree() {
        // A semiprime whose factors exceed r, so a congruence must fail.
        let n = Integer::from(1009u64 * 1013);
        let runs: Vec<Trace> = [Some(1), Some(2), Some(5), None]
            .into_iter()
            .map(|threads| {
                aks(
                    &n,
                    Options {
                        threads,
                        ..Options::default()
                    },
                )
                .unwrap()
            })
            .collect();
        assert!(matches!(runs[0].reason, Reason::CongruenceFails(_)));
        assert!(runs.iter().all(|t| t.reason == runs[0].reason));
        assert_eq!(runs[2].threads, 5);
    }

    #[test]
    fn uncertified_modes() {
        let o = Options {
            range: ARange::LambdaMinusOne,
            ..Options::default()
        };
        let t = aks(&Integer::from(10007), o).unwrap();
        assert!(t.prime && !t.certified());
        assert_eq!(t.a_max, Some(171));
    }

    #[test]
    fn multipliers_parse() {
        assert_eq!("schoolbook".parse(), Ok(Multiplier::Schoolbook));
        assert!("fft".parse::<Multiplier>().is_err());
    }
}
