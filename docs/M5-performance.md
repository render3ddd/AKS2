# M5: performance, first round

## Where the time went

Callgrind on `aks -j 1 1000003` with the M2 multiplier (r = 379):

| Cost | Share |
|---|---:|
| `mpz_import` / `mpz_export` (byte-slot packing and unpacking) | 34% |
| `malloc` / `free` (a fresh `Integer` per coefficient per product) | 27% |
| Reductions mod n (`tdiv_r`, `divrem_1`) | ≈ 10% |
| The GMP squaring itself | 9% |

The arithmetic AKS is actually about was under a tenth of the run.

## Changes

1. **Bit-packed limb slots** (`Kronecker`). Slots are exactly ⌈log₂(min(U,V)·(n−1)²)⌉ bits, with no rounding up to bytes. Packing and unpacking shift 64-bit limbs, and GMP imports and exports native limbs with a plain copy.
2. **Fold on limbs, then reduce once.** Slot i and slot i + r are added as limb arrays before any big-integer work. For n < 2⁶⁴ the reduction is word arithmetic (u128 Horner). Otherwise it reuses the output coefficient's allocation.
3. **No allocation in the power loop.** `PolyMul::square_into` writes into an existing polynomial, and `pow_x_minus_a` ping-pongs between two buffers. Multiplication by x − a is done in place, with a word-arithmetic path for n < 2⁶⁴.

After the change, the same profile executes 4× fewer instructions, and the GMP squaring is the largest single cost (28%). That is where a Kronecker implementation should spend its time.

The M2 version is kept as `KroneckerBytes` (`--mul kronecker-bytes`) for comparison and as a second test oracle.

## Results (Linux x86-64 container)

Whole test, best of three:

| n | r | M3/M4, 1 thread | **M5, 1 thread** | **M5, 2 threads** | C++ v1.1 Zn[x] (NTL), 2 threads | C++ v1.1 Z, 2 threads |
|---|---:|---:|---:|---:|---:|---:|
| 23456789 | 599 | 2.12 s | **0.56 s** | **0.30 s** | 1.37 s | 2.0 s |
| 1000000007 | 857 | 6.39 s | **2.23 s** | **1.17 s** | 3.50 s | 5.9 s |
| 1003026954441971 | 2441 | — | — | **38 s** | ≈ 100 s | — |

Test suite: the 10⁵ sweep went from 866 s to 179 s, and the 8–10-digit Appendix B primes run in 14 s.

One squaring (`cargo bench --bench poly`):

| r | 10-digit: M5 | 10-digit: M2 bytes | 16-digit: M5 | 16-digit: M2 bytes | 42-digit: M5 | 42-digit: M2 bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 30 µs | 70 µs | 48 µs | 108 µs | 171 µs | 278 µs |
| 1024 | 172 µs | 347 µs | 273 µs | 495 µs | 810 µs | 1.27 ms |
| 4096 | 0.91 ms | 1.48 ms | 1.26 ms | 2.06 ms | 4.05 ms | 5.43 ms |
| 19051 | 5.3 ms | 7.9 ms | 8.3 ms | 11.3 ms | 23.5 ms | 32.3 ms |

Squaring alone is 1.4–2.3× faster. The whole test gains more (3–4×) because the power loop no longer allocates and the x − a steps use word arithmetic.

## Against the plan's goals

- *"Kronecker at least matching the corrected C++ Z version, and within about 1.5× of NTL."* Exceeded: 2.4–6.6× faster than the C++ Z version and 2.6–4.5× faster than the C++ NTL build on the same machine and thread count.
- *Karatsuba and NTT multipliers.* Not needed for speed on the CPU: Kronecker hands one large multiplication to GMP, which already switches between Toom–Cook and FFT by size. An NTT is still the natural GPU kernel (M7).

## Further ideas (not done)

- Montgomery or Barrett reduction for n ≥ 2⁶⁴ (today: GMP `tdiv_r` per coefficient).
- Packing straight from the previous product's limbs, skipping the `Integer` coefficients in between.
- A word-sized coefficient representation (`Vec<u64>`) for n < 2⁶⁴, which would remove the remaining `Integer` overhead entirely.
- Reusing the GMP product's buffer between squarings (needs `unsafe` or a `rug` API for limb access; the crate is `#![forbid(unsafe_code)]`).
