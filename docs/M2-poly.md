# M2: the cyclic polynomial ring and its multipliers

`src/poly/` implements ℤₙ[x]/(xʳ − 1) (`Ring`, `CyclicPoly`) and two multipliers behind the `PolyMul` trait:

- **`Schoolbook`**: O(r²) cyclic convolution with GMP multiply-accumulate. The oracle for every test.
- **`Kronecker`**: binary segmentation, the dissertation's Z version. Unsigned byte-aligned slots with 2ᵇ > min(U, V)·(n − 1)², one GMP multiply (or square), then unpack, fold and reduce mod n. Zero operands are handled up front.

Folding happens before reduction ("lazy reduction"), so each product does r reductions mod n rather than 2r − 1.

## Tests

`tests/poly_props.rs` (proptest) checks Kronecker multiply and square against Schoolbook, the ring laws, and `pow` against repeated multiplication. Inputs:

- n up to 512 bits, including 2ᵏ and 2ᵏ ± 1;
- random polynomials, plus the edge cases that broke the 2021 C++: zero, constants, all-ones, all-(n − 1) and sparse.

The debug suite uses r ≤ 64. A release-only test (`--include-ignored`) goes to r ≤ 4096.

## Benchmark

`cargo bench --bench poly`. Linux x86-64 container, one thread, median time for one squaring:

| r | 10-digit n (Kronecker) | 10-digit n (Schoolbook) | 16-digit n (Kronecker) | 42-digit n (Kronecker) |
|---:|---:|---:|---:|---:|
| 256 | 0.08 ms | 0.77 ms | 0.11 ms | 0.30 ms |
| 512 | 0.17 ms | 3.2 ms | 0.25 ms | 0.68 ms |
| 1024 | 0.47 ms | 14 ms | 0.54 ms | 1.4 ms |
| 2048 | 0.88 ms | — | 1.1 ms | 3.0 ms |
| 4096 | 1.6 ms | — | 2.2 ms | 6.3 ms |
| 8192 | 3.6 ms | — | 5.0 ms | 14 ms |
| 19051 | 8.8 ms | — | 11 ms | 35 ms |

Kronecker roughly doubles per doubling of r; Schoolbook quadruples. The 2021 code's 4× per doubling was a symptom of its O(D²·b) unpacking loop, not of the method.

A rough estimate for 1000000007 (r = 857, a ≤ 874, about 30 squarings per power): 874 × 30 × 0.4 ms ≈ 10 s serial, with multiplications by x − a extra. The corrected C++ Z version takes 5.9 s on 2 threads. Two things should close the gap: multiplying by x − a in O(r) as a shift and subtract (M3), and parallelism (M4). Coefficient reduction and per-coefficient allocation are the next targets (M5).
