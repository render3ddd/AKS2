//! M3: verdicts and parameters against the C++ v1.1 reference data.

mod common;

use aks2::aks::{Options, Reason, aks};
use rug::Integer;

fn check(row: &common::CppRow) {
    let t = aks(&Integer::from(row.n), Options::default()).unwrap();
    let ctx = format!("n = {}: Rust {t:?}\nC++ {row:?}", row.n);
    assert_eq!(t.prime, row.prime, "{ctx}");
    let (reason, witness) = match t.reason {
        Reason::Small => ("small", None),
        Reason::PerfectPower { .. } => ("perfect_power", None),
        Reason::Divisor(d) => ("divisor", Some(d)),
        Reason::NoRBelowN => ("n_le_r", None),
        Reason::CongruencesHold => ("congruence", None),
        Reason::CongruenceFails(a) => ("congruence", Some(a)),
    };
    assert_eq!(reason, row.reason, "{ctx}");
    assert_eq!(witness, row.witness, "{ctx}");
    assert_eq!(t.params.map(|p| p.r), row.r, "{ctx}");
    assert_eq!(t.params.map(|p| p.phi), row.phi, "{ctx}");
    assert_eq!(t.params.map(|p| p.lambda), row.lambda, "{ctx}");
    assert_eq!(t.a_max, row.a_max, "{ctx}");
    assert!(t.certified() || !t.prime);
}

#[test]
fn matches_cpp_up_to_1000() {
    common::cpp_reference()
        .iter()
        .filter(|r| r.n <= 1000)
        .for_each(check);
}

/// Every n ≤ 10⁴ and the Appendix B primes up to 7 digits. Run by CI's release job.
#[test]
#[ignore = "release-mode only: cargo test --release -- --include-ignored"]
fn matches_cpp_up_to_7_digits() {
    common::cpp_reference()
        .iter()
        .filter(|r| r.n > 1000 && r.n < 10_000_000)
        .for_each(check);
}

/// The Appendix B primes with 8–10 digits: minutes in serial. Opt-in.
#[test]
#[ignore = "slow: cargo test --release -- --ignored slow_"]
fn slow_matches_cpp_appendix_b_8_to_10_digits() {
    common::cpp_reference()
        .iter()
        .filter(|r| r.n >= 10_000_000)
        .for_each(check);
}
