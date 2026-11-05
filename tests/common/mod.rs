//! Shared test helpers.

/// One row of `tests/data/cpp_reference.csv`, exported from the corrected C++
/// (render3d/AKS tag `v1.1`) by `tools/cpp_reference.py`.
#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)] // not every test binary reads every field
pub struct CppRow {
    pub n: u64,
    pub prime: bool,
    pub reason: String,
    pub witness: Option<u64>,
    pub r: Option<u64>,
    pub phi: Option<u64>,
    pub lambda: Option<u64>,
    pub a_max: Option<u64>,
}

fn opt(s: &str) -> Option<u64> {
    (!s.is_empty()).then(|| s.parse().expect("integer field"))
}

/// Reads the whole reference file (it is ~10,000 rows, so this is cheap).
#[allow(dead_code)]
pub fn cpp_reference() -> Vec<CppRow> {
    let text = include_str!("../data/cpp_reference.csv");
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("n,prime,reason,witness,r,phi,lambda,a_max"),
        "unexpected header"
    );
    lines
        .map(|line| {
            let f: Vec<&str> = line.split(',').collect();
            assert_eq!(f.len(), 8, "bad row: {line}");
            CppRow {
                n: f[0].parse().expect("n"),
                prime: f[1] == "1",
                reason: f[2].to_owned(),
                witness: opt(f[3]),
                r: opt(f[4]),
                phi: opt(f[5]),
                lambda: opt(f[6]),
                a_max: opt(f[7]),
            }
        })
        .collect()
}
