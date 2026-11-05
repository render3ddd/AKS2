# aks2

A Rust rewrite of my [C++ AKS project](https://github.com/render3d/AKS): the Lenstra (2002) variant of the Agrawal–Kayal–Saxena deterministic primality test, originally written for my 2021 University of Bath MSc dissertation, *Implementing alternative polynomial multiplication algorithms and the AKS primality test in parallel*.

The corrected C++ (tag [`v1.1`](https://github.com/render3d/AKS/releases/tag/v1.1)) is the numerical and performance reference for this port. [`RUST_PLAN.md`](RUST_PLAN.md) sets out the design and milestones.

## Status

| Milestone | Contents | State |
|---|---|---|
| M0 | Crate setup, licence, CI, C++ reference data | done |
| M1 | Number theory: log₂, perfect powers, multiplicative order, φ, λ, choice of r | done |
| M2 | Polynomial ring Zₙ[x]/(xʳ − 1), schoolbook and Kronecker (binary segmentation) multiplication | done |
| M3 | Serial Lenstra AKS and the `aks` command | done |
| M4 | Parallel congruence tests (rayon) | done |
| M5 | Performance | — |

## Building

Needs a Rust toolchain (1.85 or newer) and a C compiler: the [`rug`](https://crates.io/crates/rug) crate builds GMP from source on the first build, which takes a few minutes. On macOS the Xcode Command Line Tools are enough.

```bash
cargo build --release
cargo test
```

## Usage

```bash
cargo run --release -- 23456789               # r, parameters, verdict, time
echo 1000000007 | cargo run --release -- -q   # n from stdin; verdict only
cargo run --release -- --lambda 10007         # dissertation's λ(r) − 1 bound (not certified)
cargo run --release -- --mul schoolbook --log logs/run.csv 10007
cargo run --release -- -j 4 1000000007         # 4 threads (default: one per core)
```

The exit status is 0 for prime and 1 otherwise, as with the C++ programs.

## Reference data

`tests/data/cpp_reference.csv` holds the corrected C++'s verdict and parameters (r, φ(r), λ(r), and the theorem's bound on *a*) for every n ≤ 10,000 and for the dissertation's Appendix B primes up to 10 digits. To regenerate it from a built C++ checkout:

```bash
python3 tools/cpp_reference.py --bin ../AKS/devVR/LenstraZnx.out --out tests/data/cpp_reference.csv
tools/cpp_baselines.sh ../AKS     # C++ timings for RUST_PLAN.md §2
```

## Citing

If you use this code, please cite it: see [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button for it).

## Licence

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option. Both require the copyright notice, and so the author's credit, to be kept in copies. The [`rug`](https://crates.io/crates/rug) dependency and GMP are LGPL-3.0-or-later; anyone distributing binaries must also meet the LGPL's terms for those libraries.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in this work, as defined in the Apache-2.0 licence, shall be dual licensed as above, without any additional terms or conditions.
