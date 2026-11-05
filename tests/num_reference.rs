//! M1: the number theory layer against independent and C++ reference data.

mod common;

use aks2::aks::{ARange, RSearch, a_bound, find_r};
use aks2::num::{is_perfect_power, totients};
use rug::Integer;

#[test]
fn totients_match_sympy() {
    let text = include_str!("data/totients.csv");
    let mut count = 0;
    for line in text.lines().skip(1) {
        let v: Vec<u64> = line.split(',').map(|f| f.parse().unwrap()).collect();
        assert_eq!(totients(v[0]), (v[1], v[2]), "r = {}", v[0]);
        count += 1;
    }
    assert_eq!(count, 20_000);
}

/// Every row of the C++ v1.1 export: perfect-power detection, the r search
/// (divisor found, r found, or none below n), φ(r), λ(r) and the theorem bound.
#[test]
fn parameters_match_cpp_reference() {
    for row in common::cpp_reference() {
        let n = Integer::from(row.n);
        let ctx = format!("n = {}, C++ row {row:?}", row.n);
        match row.reason.as_str() {
            "small" => continue,
            "perfect_power" => {
                assert!(is_perfect_power(&n), "{ctx}");
                continue;
            }
            _ => assert!(!is_perfect_power(&n), "{ctx}"),
        }
        match (find_r(&n), row.reason.as_str()) {
            (RSearch::Divisor(d), "divisor") => assert_eq!(Some(d), row.witness, "{ctx}"),
            (RSearch::NoRBelowN, "n_le_r") => {}
            (RSearch::Found(p), "congruence") => {
                assert_eq!(Some(p.r), row.r, "{ctx}");
                assert_eq!(Some(p.phi), row.phi, "{ctx}");
                assert_eq!(Some(p.lambda), row.lambda, "{ctx}");
                assert_eq!(Some(a_bound(ARange::Theorem, &n, &p)), row.a_max, "{ctx}");
            }
            (got, _) => panic!("{ctx}: Rust gave {got:?}"),
        }
    }
}
