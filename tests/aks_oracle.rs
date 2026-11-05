//! M3: verdicts against GMP's probable-prime test, which is deterministic for
//! n < 2⁶⁴ (BPSW has no known counterexample, and none exist below 2⁶⁴), and
//! composite families chosen to reach the congruence step.

use aks2::aks::{Multiplier, Options, Reason, aks, is_prime};
use rug::Integer;
use rug::integer::IsPrime;
use rug::ops::Pow;

fn oracle(n: &Integer) -> bool {
    n.is_probably_prime(30) != IsPrime::No
}

fn sweep(range: std::ops::RangeInclusive<u64>) {
    for n in range {
        let n = Integer::from(n);
        assert_eq!(is_prime(&n), oracle(&n), "n = {n}");
    }
}

#[test]
fn sweep_to_3000() {
    sweep(2..=3000);
}

/// CI release job: every n up to 2·10⁴.
#[test]
#[ignore = "release-mode only: cargo test --release -- --include-ignored"]
fn sweep_to_20000() {
    sweep(3001..=20_000);
}

/// The plan's 10⁵ sweep: ~9,600 primes at a fraction of a second each. Opt-in.
#[test]
#[ignore = "slow: cargo test --release -- --ignored slow_"]
fn slow_sweep_to_100000() {
    sweep(20_001..=100_000);
}

fn primes_from(start: u64, count: usize) -> Vec<Integer> {
    let mut p = Integer::from(start);
    (0..count)
        .map(|_| {
            p.next_prime_mut();
            p.clone()
        })
        .collect()
}

/// Products of two primes larger than r: no gcd catches them, so only the
/// congruences can, and every multiplier must find the same smallest failing a.
#[test]
fn semiprimes_fail_a_congruence() {
    let ps = primes_from(1000, 4);
    for (i, p) in ps.iter().enumerate() {
        for q in &ps[i + 1..] {
            let n = Integer::from(p * q);
            let runs: Vec<_> = Multiplier::ALL
                .iter()
                .map(|&multiplier| {
                    aks(
                        &n,
                        Options {
                            threads: Some(1),
                            multiplier,
                            ..Options::default()
                        },
                    )
                    .unwrap()
                })
                .collect();
            let t = &runs[0];
            assert!(!t.prime, "n = {n}");
            assert!(
                matches!(t.reason, Reason::CongruenceFails(_)),
                "n = {n}: {:?}",
                t.reason
            );
            assert!(
                t.params.unwrap().r < 1000,
                "n = {n}: r too big for this test"
            );
            assert!(
                runs.iter().all(|u| u.reason == t.reason),
                "n = {n}: multipliers disagree"
            );
        }
    }
}

/// Carmichael numbers fool Fermat's test; AKS must not be fooled.
#[test]
fn carmichael_numbers() {
    let small = [
        561u64, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 62745, 63973, 75361, 101101,
    ];
    for n in small {
        assert!(!is_prime(&Integer::from(n)), "n = {n}");
    }
    // (6k+1)(12k+1)(18k+1) with all three prime is Carmichael (Chernick).
    // k = 370: factors 2221, 4441, 6661 all exceed r, so only a congruence can catch it.
    let k = Integer::from(370u32);
    let f = |m: u32| Integer::from(&k * m) + 1u32;
    let (a, b, c) = (f(6), f(12), f(18));
    assert!(
        [&a, &b, &c].iter().all(|x| oracle(x)),
        "k = 370 must give three primes"
    );
    let n = a * b * c;
    let t = aks(&n, Options::default()).unwrap();
    assert!(!t.prime, "n = {n}: {:?}", t.reason);
    assert!(
        matches!(t.reason, Reason::CongruenceFails(_)),
        "n = {n}: {:?}",
        t.reason
    );
}

/// Prime powers, including bases beyond 2⁶⁴ that the 2021 PerfectPower.h missed.
#[test]
fn prime_powers() {
    for (p, e) in [(2u64, 10u32), (3, 7), (10007, 2), (65_537, 3)] {
        let n = Integer::from(p).pow(e);
        let t = aks(&n, Options::default()).unwrap();
        assert!(matches!(t.reason, Reason::PerfectPower { .. }), "{p}^{e}");
    }
    let p = primes_from(u64::MAX, 1).remove(0);
    for e in [2u32, 3, 5] {
        let n = p.clone().pow(e);
        let t = aks(&n, Options::default()).unwrap();
        assert_eq!(
            t.reason,
            Reason::PerfectPower {
                base: p.clone(),
                exp: e
            }
        );
    }
}
