# M3: serial Lenstra AKS

`aks2::aks::aks(n, Options)` runs the whole test and returns a `Trace`: the verdict, the reason it stopped, r, ord_r(n), φ(r), λ(r), the bound on *a*, and the time taken. The `aks` binary wraps it with the corrected C++'s command line and CSV log format.

Steps, as in the corrected C++:

1. Perfect power?
2. Search for r, with a gcd check at each step.
3. Congruences (x − a)ⁿ ≡ x^(n mod r) − a for 1 ≤ a ≤ bound.

There is no separate "Z version" and "Zn[x] version" any more. Both are the same algorithm with a different `PolyMul`: `--mul kronecker` (binary segmentation, the default) or `--mul schoolbook`.

The one structural speed-up over the C++: (x − a)ⁿ is computed by square-and-multiply. Each multiply step is by x − a, so it is done as a rotation minus a scalar multiple, in O(r). Only the squarings use the multiplier.

## Correctness

| Check | Scope | Where |
|---|---|---|
| Verdict, reason, witness, r, φ, λ, a-bound = C++ v1.1 | all 10,028 reference rows (n ≤ 10⁴ and the 30 Appendix B primes ≤ 10 digits) | `tests/aks_reference.rs` |
| Verdict = GMP `is_probably_prime(30)` | every n ≤ 20,000 (CI); every n ≤ 10⁵ (opt-in `slow_`) | `tests/aks_oracle.rs` |
| Composites that only a congruence can catch | semiprimes p·q with p, q > r; Chernick Carmichael number 2221·4441·6661 | `tests/aks_oracle.rs` |
| Prime powers | including (2⁶⁴ + 13)ᵉ, beyond the 2021 `PerfectPower.h` | `tests/aks_oracle.rs` |

Debug builds run n ≤ 1000 against the reference and n ≤ 3000 against the oracle. The CI release job adds n ≤ 10⁷ against the reference and n ≤ 20,000 against the oracle. The 8–10-digit Appendix B primes (90 s serially) and the 10⁵ sweep are opt-in: `cargo test --release -- --ignored slow_`.

## Speed (serial, Linux x86-64 container)

| n | r | Rust, 1 thread | C++ Zn[x] (NTL), 2 threads | C++ Z, 2 threads |
|---|---:|---:|---:|---:|
| 23456789 | 599 | 2.1 s | 1.4 s | 2.0 s |
| 1000000007 | 857 | 6.5 s | 3.5 s | 5.9 s |

On one thread, Rust nearly matches the C++ on two.
