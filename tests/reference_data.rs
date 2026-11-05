//! Sanity checks on the C++ reference data itself, so later milestones can trust it.

mod common;

fn is_prime_trial(n: u64) -> bool {
    n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0)
}

#[test]
fn covers_every_n_up_to_ten_thousand() {
    let rows = common::cpp_reference();
    let small: Vec<u64> = rows.iter().map(|r| r.n).filter(|&n| n <= 10_000).collect();
    assert_eq!(small, (2..=10_000).collect::<Vec<_>>());
}

#[test]
fn verdicts_match_trial_division() {
    for row in common::cpp_reference() {
        assert_eq!(row.prime, is_prime_trial(row.n), "n = {}", row.n);
    }
}

#[test]
fn rows_are_internally_consistent() {
    for row in common::cpp_reference() {
        match row.reason.as_str() {
            "small" => assert!(row.n == 2 || row.n == 3),
            "perfect_power" | "divisor" => assert!(!row.prime, "n = {}", row.n),
            "n_le_r" => assert!(row.prime && row.r.is_none(), "n = {}", row.n),
            "congruence" => {
                let (r, phi, lambda, a) = (
                    row.r.unwrap(),
                    row.phi.unwrap(),
                    row.lambda.unwrap(),
                    row.a_max.unwrap(),
                );
                assert!(r < row.n && phi < r && phi % lambda == 0, "n = {}", row.n);
                assert!(a < row.n, "a_max must be capped below n (n = {})", row.n);
            }
            other => panic!("unknown reason {other:?} for n = {}", row.n),
        }
        if row.reason == "divisor" {
            let d = row.witness.unwrap();
            assert!(d > 1 && d < row.n && row.n % d == 0, "n = {}", row.n);
        }
    }
}

#[test]
fn spot_values_from_rust_plan() {
    let rows = common::cpp_reference();
    let get = |n| rows.iter().find(|r| r.n == n).unwrap();
    assert_eq!(get(23_456_789).r, Some(599));
    assert_eq!(get(1_000_000_007).r, Some(857));
}
