//! Word-sized and big-integer number theory used to choose the AKS parameters.
//!
//! Everything about `r` lives in `u64`: in AKS, r is at most about
//! (log₂ n)⁵, so it fits in 32 bits for any n anyone will test. Only `n` itself
//! is a big [`Integer`].

use rug::Integer;
use rug::ops::Pow;

/// ⌊log₂ n⌋ for n ≥ 1, from the bit length (never from a float logarithm;
/// the 2021 code used the natural log here, see REFACTOR.md §1.1).
///
/// # Panics
/// If `n < 1`.
pub fn log2_floor(n: &Integer) -> u32 {
    assert!(*n >= 1, "log2_floor needs n >= 1");
    n.significant_bits() - 1
}

/// log₂ n as an `f64`, accurate to about one ulp for any size of n.
///
/// Used only for the theorem's bound on *a*, which the C++ reference also
/// computes in `double`.
///
/// # Panics
/// If `n < 1`.
pub fn log2_f64(n: &Integer) -> f64 {
    assert!(*n >= 1, "log2_f64 needs n >= 1");
    // n = m · 2^e with 0.5 ≤ m < 1, so log₂ n = e + log₂ m.
    let (m, e) = n.to_f64_exp();
    f64::from(e) + m.log2()
}

/// Greatest common divisor.
pub const fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Least common multiple.
///
/// # Panics
/// On overflow (debug builds), which cannot happen for λ(r) since λ(r) < r.
pub const fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 {
        0
    } else {
        a / gcd(a, b) * b
    }
}

/// a · b mod m without overflow.
pub const fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// aᵉ mod m by square-and-multiply. Returns 0 when m = 1.
pub const fn pow_mod(mut a: u64, mut e: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut acc = 1;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod(acc, a, m);
        }
        a = mul_mod(a, a, m);
        e >>= 1;
    }
    acc
}

/// Prime factorisation by trial division, as `(prime, exponent)` pairs in
/// increasing order. `factor_u64(1)` is empty. Fast for the r values AKS
/// uses (r < 2³²); correct, if slow, for any u64.
///
/// # Panics
/// If `n == 0`.
pub fn factor_u64(mut n: u64) -> Vec<(u64, u32)> {
    assert!(n > 0, "factor_u64(0) is undefined");
    let mut out = Vec::new();
    let mut take = |n: &mut u64, p: u64| {
        let mut e = 0;
        while *n % p == 0 {
            *n /= p;
            e += 1;
        }
        if e > 0 {
            out.push((p, e));
        }
    };
    take(&mut n, 2);
    take(&mut n, 3);
    // Remaining primes are 6k ± 1.
    let mut p = 5;
    while p <= n / p {
        take(&mut n, p);
        take(&mut n, p + 2);
        p += 6;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

/// Euler's totient φ(r) and Carmichael's function λ(r), from one factorisation.
///
/// λ(r) is the exponent of the unit group (ℤ/rℤ)ˣ: the lcm of λ(pᵏ) over the
/// prime powers of r, where λ(pᵏ) = φ(pᵏ) except λ(2ᵏ) = 2ᵏ⁻² for k ≥ 3.
/// Both are 1 for r = 1.
///
/// # Panics
/// If `r == 0`.
pub fn totients(r: u64) -> (u64, u64) {
    let mut phi = 1;
    let mut lambda = 1;
    for (p, k) in factor_u64(r) {
        let pk1 = p.pow(k - 1); // p^(k-1)
        let phi_pk = (p - 1) * pk1;
        let lambda_pk = if p == 2 && k >= 3 { phi_pk / 2 } else { phi_pk };
        phi *= phi_pk;
        lambda = lcm(lambda, lambda_pk);
    }
    (phi, lambda)
}

/// Euler's totient φ(r). See [`totients`].
pub fn euler_phi(r: u64) -> u64 {
    totients(r).0
}

/// Carmichael's function λ(r). See [`totients`].
pub fn carmichael_lambda(r: u64) -> u64 {
    totients(r).1
}

/// The multiplicative order of a modulo r: the least k ≥ 1 with aᵏ ≡ 1 (mod r).
///
/// Returns `None` if gcd(a, r) ≠ 1 (no such k), and `Some(1)` for r = 1.
/// Works down from λ(r), which every order divides, so it costs a handful of
/// modular exponentiations rather than up to r multiplications.
pub fn mult_order(a: u64, r: u64) -> Option<u64> {
    if r == 1 {
        return Some(1);
    }
    let a = a % r;
    if gcd(a, r) != 1 {
        return None;
    }
    let mut ord = carmichael_lambda(r);
    for (p, _) in factor_u64(ord) {
        while ord % p == 0 && pow_mod(a, ord / p, r) == 1 {
            ord /= p;
        }
    }
    Some(ord)
}

/// Whether ord_r(a) > k, for gcd(a, r) = 1. Cheap when λ(r) ≤ k, since every
/// order divides λ(r).
pub fn order_exceeds(a: u64, r: u64, k: u64) -> bool {
    carmichael_lambda(r) > k && mult_order(a, r).is_some_and(|o| o > k)
}

/// If n = aᵇ with b ≥ 2, returns `(a, b)` with the **largest** such b (so the
/// smallest base). Exact for every size of n, unlike the 2021 `PerfectPower.h`
/// (REFACTOR.md §1.3). Returns `None` for non-powers and for n < 4
/// (0 and 1 are trivial powers, but AKS only asks about n > 1).
pub fn perfect_power(n: &Integer) -> Option<(Integer, u32)> {
    if *n < 4 || !n.is_perfect_power() {
        return None;
    }
    // n ≥ 4, so ⌊log₂ n⌋ ≥ 2 bounds the exponent; try the largest first.
    (2..=log2_floor(n)).rev().find_map(|b| {
        let (root, rem) = n.clone().root_rem(Integer::new(), b);
        (rem == 0).then(|| {
            debug_assert_eq!(root.clone().pow(b), *n);
            (root, b)
        })
    })
}

/// Whether n = aᵇ for some integers a, b ≥ 2.
pub fn is_perfect_power(n: &Integer) -> bool {
    *n >= 4 && n.is_perfect_power()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int(n: u64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn log2_floor_is_bit_length_minus_one() {
        assert_eq!(log2_floor(&int(1)), 0);
        assert_eq!(log2_floor(&int(2)), 1);
        assert_eq!(log2_floor(&int(3)), 1);
        assert_eq!(log2_floor(&int(1023)), 9);
        assert_eq!(log2_floor(&int(1024)), 10);
        let big = Integer::from(1) << 1000;
        assert_eq!(log2_floor(&big), 1000);
        assert_eq!(log2_floor(&(big - 1u32)), 999);
    }

    #[test]
    fn log2_f64_matches_std_and_scales() {
        for n in [1u64, 2, 3, 10, 10007, 23_456_789, u64::MAX] {
            let want = (n as f64).log2();
            assert!((log2_f64(&int(n)) - want).abs() < 1e-12, "n = {n}");
        }
        let big = Integer::from(1) << 5000;
        assert_eq!(log2_f64(&big), 5000.0);
    }

    #[test]
    fn factorisation_round_trips() {
        for n in 1..5000u64 {
            let f = factor_u64(n);
            assert_eq!(f.iter().map(|&(p, k)| p.pow(k)).product::<u64>(), n);
            assert!(f.windows(2).all(|w| w[0].0 < w[1].0));
        }
        assert_eq!(factor_u64(4_294_967_291), vec![(4_294_967_291, 1)]); // largest 32-bit prime
        assert_eq!(factor_u64(1 << 40), vec![(2, 40)]);
    }

    #[test]
    fn totients_brute_force() {
        for r in 1..=1000u64 {
            let units: Vec<u64> = (1..=r).filter(|&a| gcd(a, r) == 1).collect();
            let phi = units.len() as u64;
            // λ(r): least k with a^k ≡ 1 for every unit (r = 1: the empty group, 1).
            let lambda = (1..=r)
                .find(|&k| units.iter().all(|&a| pow_mod(a, k, r) == 1 % r))
                .unwrap();
            assert_eq!(totients(r), (phi, lambda), "r = {r}");
        }
    }

    #[test]
    fn mult_order_brute_force() {
        for r in 1..=300u64 {
            for a in 0..r {
                let want =
                    (gcd(a, r) == 1).then(|| (1..=r).find(|&k| pow_mod(a, k, r) == 1 % r).unwrap());
                assert_eq!(mult_order(a, r), want, "a = {a}, r = {r}");
            }
        }
    }

    #[test]
    fn order_exceeds_agrees_with_mult_order() {
        for r in 2..=300u64 {
            for a in (1..r).filter(|&a| gcd(a, r) == 1) {
                let o = mult_order(a, r).unwrap();
                for k in [0, o - 1, o, o + 1, r] {
                    assert_eq!(order_exceeds(a, r, k), o > k, "a={a} r={r} k={k}");
                }
            }
        }
    }

    #[test]
    fn perfect_powers() {
        assert_eq!(perfect_power(&int(64)), Some((int(2), 6)));
        assert_eq!(perfect_power(&int(36)), Some((int(6), 2)));
        assert_eq!(perfect_power(&int(4)), Some((int(2), 2)));
        for n in [0u64, 1, 2, 3, 5, 6, 12, 10007] {
            assert_eq!(perfect_power(&int(n)), None, "n = {n}");
        }
        // Bases beyond 2^64, which the 2021 PerfectPower.h missed.
        let base: Integer = (Integer::from(1) << 70) + 25u32;
        let n = base.clone().pow(3);
        assert_eq!(perfect_power(&n), Some((base, 3)));
        assert!(!is_perfect_power(&(n + 1u32)));
        // Exhaustive against brute force for small n.
        for n in 2..20_000u64 {
            let brute = (2..=14u32).any(|b| {
                let r = (n as f64).powf(1.0 / f64::from(b)).round() as u64;
                (r.saturating_sub(1)..=r + 1).any(|x| x >= 2 && x.checked_pow(b) == Some(n))
            });
            assert_eq!(is_perfect_power(&int(n)), brute, "n = {n}");
            assert_eq!(perfect_power(&int(n)).is_some(), brute, "n = {n}");
        }
    }
}
