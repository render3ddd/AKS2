//! Choosing r and the range of *a* (REFACTOR.md §1.1, §1.2, §2.4).

use std::fmt;

use rug::Integer;

use crate::num::{gcd, log2_f64, log2_floor, mult_order, totients};

/// How far the congruence tests (x + a)ⁿ ≡ xⁿ + a run, as 1 ≤ a ≤ bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ARange {
    /// ⌊√φ(r) · log₂ n⌋, capped at n − 1: the bound the Lenstra/AKS proof needs.
    #[default]
    Theorem,
    /// φ(r) − 1: the dissertation's Algorithm 2.4.2. **Not certified.**
    PhiMinusOne,
    /// λ(r) − 1: Egan (2005), dissertation §5.2.4. Can fall below the
    /// theorem's bound, so a PRIME verdict is **not certified**.
    LambdaMinusOne,
}

impl ARange {
    /// Whether a PRIME verdict in this mode is backed by the theorem.
    pub const fn is_certified(self) -> bool {
        matches!(self, ARange::Theorem)
    }
}

impl fmt::Display for ARange {
    /// The same wording as the C++ reference's `ARangeName`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ARange::Theorem => "theorem: floor(sqrt(phi(r)) log2 n)",
            ARange::PhiMinusOne => "phi(r) - 1 [not certified]",
            ARange::LambdaMinusOne => "lambda(r) - 1 [not certified]",
        })
    }
}

/// The parameters AKS needs once a suitable r is found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Params {
    /// The smallest r with gcd(r, n) = 1 and ord_r(n) > ⌊log₂ n⌋².
    pub r: u64,
    /// ord_r(n), the multiplicative order of n modulo r.
    pub order: u64,
    /// φ(r).
    pub phi: u64,
    /// λ(r).
    pub lambda: u64,
}

/// The outcome of the search for r (steps 2–4 of the algorithm).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RSearch {
    /// A suitable r exists below n, and gcd(a, n) = 1 for all a ≤ r.
    Found(Params),
    /// Some r′ ≤ r shares the factor `d` with n, 1 < d < n: n is COMPOSITE.
    Divisor(u64),
    /// No suitable r below n, and gcd(a, n) = 1 for every a < n: n is PRIME
    /// (the algorithm's "if n ≤ r, output PRIME").
    NoRBelowN,
}

/// Searches r = 2, 3, … for the smallest r with ord_r(n) > ⌊log₂ n⌋², checking
/// gcd(r, n) on the way, exactly as the corrected C++ `FindR` loop does.
///
/// Orders are only computed for r > ⌊log₂ n⌋² + 1 with λ(r) > ⌊log₂ n⌋²,
/// since ord_r(n) ≤ λ(r) < r. That skips almost all the multiplications the C++
/// does, without changing the answer.
///
/// # Panics
/// If n < 2, or if r would exceed 2³² (n would need about 2¹⁶ⁿ digits).
pub fn find_r(n: &Integer) -> RSearch {
    assert!(*n >= 2, "find_r needs n >= 2");
    let log2n = u64::from(log2_floor(n));
    let k = log2n * log2n; // need ord_r(n) > k

    let mut r: u64 = 2;
    while *n > r {
        let r32 = u32::try_from(r).expect("r exceeds 2^32");
        let n_mod_r = u64::from(n.mod_u(r32));
        let g = gcd(r, n_mod_r);
        if g != 1 {
            return RSearch::Divisor(g);
        }
        if r > k + 1 {
            let (phi, lambda) = totients(r);
            if lambda > k {
                let order = mult_order(n_mod_r, r).expect("gcd(n, r) = 1");
                if order > k {
                    return RSearch::Found(Params {
                        r,
                        order,
                        phi,
                        lambda,
                    });
                }
            }
        }
        r += 1;
    }
    RSearch::NoRBelowN
}

/// The upper bound on *a* for the chosen mode.
///
/// The theorem bound is ⌊√φ · log₂ n⌋ in `f64`, as in the C++ reference, capped
/// at n − 1 because *a* must stay below n.
pub fn a_bound(mode: ARange, n: &Integer, params: &Params) -> u64 {
    match mode {
        ARange::Theorem => {
            let bound = ((params.phi as f64).sqrt() * log2_f64(n)).floor() as u64;
            if *n <= bound {
                n.to_u64().expect("n <= bound fits u64") - 1
            } else {
                bound
            }
        }
        ARange::PhiMinusOne => params.phi - 1,
        ARange::LambdaMinusOne => params.lambda - 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn found(n: &str) -> Params {
        match find_r(&n.parse().unwrap()) {
            RSearch::Found(p) => p,
            other => panic!("{n}: {other:?}"),
        }
    }

    #[test]
    fn rust_plan_values() {
        assert_eq!(found("10007").r, 173);
        assert_eq!(found("23456789").r, 599);
        assert_eq!(found("1000000007").r, 857);
        assert_eq!(found("689960931088884849033689023336009222695077").r, 19051);
    }

    #[test]
    fn small_cases() {
        assert_eq!(find_r(&Integer::from(2)), RSearch::NoRBelowN);
        assert_eq!(find_r(&Integer::from(5)), RSearch::NoRBelowN);
        assert_eq!(find_r(&Integer::from(91)), RSearch::Divisor(7));
        assert_eq!(find_r(&Integer::from(4)), RSearch::Divisor(2));
    }

    #[test]
    fn found_params_satisfy_the_definition() {
        for n in [31u64, 10007, 23_456_789, 1_000_000_007, 9_999_954_997] {
            let p = found(&n.to_string());
            let k = u64::from(log2_floor(&Integer::from(n))).pow(2);
            assert!(p.order > k);
            assert_eq!(gcd(p.r, n % p.r), 1);
            // r is the smallest: every smaller coprime r' has ord <= k.
            for r in 2..p.r {
                assert_eq!(gcd(r, n % r), 1);
                assert!(mult_order(n % r, r).unwrap() <= k, "n = {n}, r' = {r}");
            }
        }
    }

    #[test]
    fn a_bound_modes() {
        let n = Integer::from(10007);
        let p = found("10007");
        assert_eq!(a_bound(ARange::Theorem, &n, &p), 174);
        assert_eq!(a_bound(ARange::PhiMinusOne, &n, &p), 171);
        assert_eq!(a_bound(ARange::LambdaMinusOne, &n, &p), 171);
        assert!(ARange::Theorem.is_certified() && !ARange::LambdaMinusOne.is_certified());
    }
}
