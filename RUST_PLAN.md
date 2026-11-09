# AKS v2: Plan for re-implementing the C++ AKS project in Rust

This plan covers re-implementing [render3d/AKS](https://github.com/render3d/AKS) (C++/NTL/GMP) in this repository. The aim is to meet the same objectives as the 2021 MSc dissertation, *Implementing alternative polynomial multiplication algorithms and the AKS primality test in parallel*, while correcting the problems recorded in the C++ repo's `REFACTOR.md`.

**Status (November 2026).** The C++ correction phase is complete: every issue in `REFACTOR.md` is fixed on the C++ repo's `main` (`fef9c34`, tagged **`v1.1`**; the 2021 submission is tagged `v1.0`). That corrected C++ is now the **numerical and performance reference** for this port: same parameters (r, φ(r), λ(r), range of *a*), same verdicts, and timings to beat (§6, M0 and M5). **M0 is done** (§6).

---

## 1. Objectives carried over from the dissertation

| Dissertation objective (§4.1–4.2) | Rust equivalent | Done when |
|---|---|---|
| 1. Implement the Lenstra (2002) variant of AKS in the **Zn[x]** form, using fast polynomial multiplication | `aks::lenstra` with a library multiplier (schoolbook → Karatsuba → NTT) | Agrees with a probable-prime oracle for all n ≤ 10⁶ and for every prime in Appendix B |
| 2. Implement the **Z** form (Crandall & Papadopoulos), with polynomial multiplication by binary segmentation | Kronecker-substitution multiplier on GMP big integers via `rug` | Same correctness bar; within ~3× of the Zn[x] multiplier (the corrected C++ achieves about 1–3× of NTL) |
| 3. Profile binary segmentation against the library methods | `criterion` benchmarks by degree, r, and digits of n | Log–log plots of time against log n, as the supervisor suggested (p. 45) |
| 4. Run congruence tests in parallel on multi-core CPUs (the NTL thread-pool experiment) | `rayon` over the range of *a*, with early termination | Near-linear scaling on the M3 Pro's performance cores, and the smallest failing *a* reported deterministically |
| 5. Evaluate the Carmichael λ(r) versus Euler φ(r) bound on *a* | Selectable `ARange` mode; the theorem bound is the default | λ/φ experiments re-run with a **valid** r and clearly labelled "not certified" |
| 6. Explore GPGPU (originally CUDA) | Stretch milestone M7: GPU NTT for the polynomial squaring | A feasibility report and prototype; not required for v1 |
| 7. Assess practicality against other tests (Miller–Rabin, ECPP) | Comparison against `rug`'s `is_probably_prime` (GMP BPSW + Miller–Rabin) in benchmarks | A table in the README |

Out of scope for v1: AKS variants beyond Lenstra (2002), MIRACL/FLINT parity, and primality certificates.

---

## 2. Starting point

**This repository (after M0)**

- `Cargo.toml`: package **`aks2`**, edition **2024**, `rust-version = "1.85"` (the edition-2024 minimum, and `rug` 1.30's own MSRV); one dependency, `rug = "1.30"` with only its `integer` feature (no MPFR/MPC build). `num-bigint` is removed.
- A library (`src/lib.rs`) plus the `aks` binary (`src/bin/aks.rs`), both placeholders until M1/M3.
- Licence **MIT OR Apache-2.0**, with a `CITATION.cff` (software + the dissertation as preferred citation). This replaces CC BY 4.0, which is not designed for software.
- GitHub Actions (`.github/workflows/ci.yml`): fmt; clippy `-D warnings` and tests on `ubuntu-latest` and `macos-latest`; a release-mode run that includes `#[ignore]` tests except those named `slow_*`; and an MSRV `cargo check` on 1.85.
- `tests/data/cpp_reference.csv` (C++ v1.1 output for 10,028 inputs) and the scripts in `tools/` that produce it and the timings below.

**The C++ reference (AKS repo, `main` at `fef9c34`)**

- Both programs share one command line: `LenstraZnx.out | LenstraZ.out [--theorem | --phi | --lambda] [n]`; *n* from standard input if omitted; exit status 0 for prime, 1 otherwise. Each run appends r, φ(r), λ(r), the range of *a* and the time to a CSV in `logs/`.
- `make test`: 2,322 binary segmentation checks against NTL, plus an AKS sweep of n = 2…1000 against NTL's `ProbPrime`.
- Measured baselines (`tools/cpp_baselines.sh`, best of three; the C++ uses every hardware thread):

  | n | r | Zn[x] (NTL), M3 Pro, 11 threads | Z (binary segmentation), M3 Pro, 11 threads | Zn[x], Linux container, 2 threads | Z, Linux container, 2 threads |
  |---|---|---:|---:|---:|---:|
  | 23456789 | 599 | **0.23 s** | **1.03 s** | ≈ 1.4–1.9 s | ≈ 2.0–2.5 s |
  | 1000000007 | 857 | **0.53 s** | **3.11 s** | ≈ 3.5–4.8 s | ≈ 5.9–7 s |
  | 1003026954441971 (16 digits) | 2441 | not run (`LONG=1`) | — | ≈ 100 s | — |

  The M3 Pro has 11 hardware threads: 5 performance cores and 6 efficiency cores. The C++ cannot be limited to fewer threads, so like-for-like comparisons with Rust use `aks -j 11` (§6, M4).

- Profile of the Z version: big-integer squaring ≈ 25% of run time, general multiplies ≈ 10%, and **reducing the unpacked coefficients mod n ≈ 29%**, the largest single cost.

---

## 3. Key technical decisions

### 3.1 Big integers: use `rug` and drop `num-bigint`

`rug` wraps GMP, the same backend NTL used in the C++ project, so performance comparisons stay fair. The Z version's whole premise is "one fast big-integer multiply", so GMP's multiplication is essential. `rug` also provides what the C++ code had to hand-roll, or got wrong:

| Need | `rug::Integer` | Replaces (C++) |
|---|---|---|
| Perfect-power test | `is_perfect_power()` (exact) | `PerfectPower.h`, which misses bases ≥ 2⁶³ and can hang |
| ⌊log₂ n⌋ | `significant_bits() - 1` | `log(n)`, which was the natural log |
| Packing/unpacking for Kronecker substitution | `from_digits::<u64>()`, `to_digits::<u64>()`, `keep_bits()` | `evaluate()` plus the `%`/`/` loop, which was O(D²·b) |
| gcd, modular powers, roots | `gcd`, `pow_mod`, `root_rem` | assorted |
| Test oracle | `is_probably_prime(reps)` | — |

`gmp-mpfr-sys` builds GMP from source on first compile, so no Homebrew is needed (expect a few minutes for the first build). On macOS it can link a system GMP with the `use-system-libs` feature instead.

`num-bigint` was removed in M0. If a pure-Rust or WASM build ever becomes a goal, it could return behind a `pure-rust` cargo feature and a `BigInt` trait.

### 3.2 Polynomial arithmetic: write our own cyclic ring

Rust has no mature NTL equivalent, but AKS needs only one ring: **Z_n[x]/(xʳ − 1)**. That ring is small enough to own outright.

- `struct CyclicPoly { coeffs: Vec<Integer> /* len r, each in [0, n) */ }`. Reduction mod xʳ − 1 is just folding coefficient i + r onto i, which is O(r) and needs no polynomial division.
- A `PolyMul` trait with interchangeable back-ends, chosen per run:
  - `Schoolbook`: the reference implementation, used as the oracle in property tests.
  - `Karatsuba`: for mid-size r.
  - `Kronecker`: **the Z version, done properly**. Pack the coefficients into fixed-width bit slots of b ≥ ⌈log₂(min(U,V)·(n−1)²)⌉ bits using `from_digits`, do one GMP multiply, then unpack with `to_digits`. Coefficients are non-negative mod n, so no sign bit is needed. The C++ fix keeps signed recovery only for `ZZX`. This is proven in the corrected C++ (`ZZpXmultiply`, byte-aligned slots): within about 1–3× of NTL, against 119–3,197× for the 2021 code.
  - `Ntt` (M5): multi-prime number-theoretic transform over 64-bit primes with CRT. This is the Zn[x] "library-class" path, analogous to NTL's FFT multiplication. It is also the natural kernel for a GPU (M7).
- Dedicated `square()` for every back-end, since powermod is mostly squarings (§5.2.2 of the dissertation). In the C++ this made squares 1.4–1.65× faster but the whole Z program only about 8% faster, because coefficient reduction dominates (§2).
- **Coefficient reduction is the next bottleneck.** Unpacking yields 2·r − 1 coefficients of up to b bits, each reduced mod n with a full division. Plan for: lazy reduction (fold mod xʳ − 1 on unreduced coefficients first, then reduce r values instead of 2r − 1), and Barrett or Montgomery reduction with a precomputed reciprocal of n (`rug` has no built-in Barrett, so implement it on `Integer` or on limbs).

### 3.3 Parallelism: `rayon`

```rust
// smallest failing a, found in parallel, deterministic result
let witness = (1..=a_max).into_par_iter().find_first(|&a| !congruence_holds(n, r, a));
```

`find_first` cancels outstanding work once a lower index has failed. That fixes both the data race in `CongruenceZnxHPC.h` and the lack of early termination. Each worker owns its buffers, so nothing is shared mutably.

### 3.4 Toolchain on the M3 Pro (macOS 27 Golden Gate)

- Stable Rust via `rustup`, target `aarch64-apple-darwin`.
- Xcode Command Line Tools are needed for the C compiler that builds GMP.
- `cargo fmt`, `cargo clippy -D warnings`, `cargo test`, `cargo bench`.
- Optional: `cargo-flamegraph` or Instruments for profiling.

---

## 4. Proposed crate layout

```
aks2/
├── Cargo.toml              # [lib] + [[bin]] aks; deps: rug, rayon, clap, csv/serde; dev: proptest, criterion
├── src/
│   ├── lib.rs
│   ├── num.rs              # log2_floor, is_perfect_power, mult_order(n mod r, r) on u64, factor_u64,
│   │                       # euler_phi, carmichael_lambda (from one factorisation)
│   ├── poly/
│   │   ├── mod.rs          # CyclicPoly over Z_n[x]/(x^r - 1), add/sub/fold/eq
│   │   ├── schoolbook.rs
│   │   ├── karatsuba.rs
│   │   ├── kronecker.rs    # the "Z version" multiplier
│   │   └── ntt.rs          # M5
│   ├── aks/
│   │   ├── params.rs       # find_r (true ord_r(n) > floor(log2 n)^2), a_bound(ARange)
│   │   ├── congruence.rs   # (x - a)^n == x^(n mod r) - a  in Z_n[x]/(x^r - 1); serial + rayon
│   │   └── lenstra.rs      # orchestration -> Verdict { Prime, Composite(reason) } + Trace
│   └── bin/aks.rs          # CLI
├── tests/                  # integration: sweeps, golden vectors, Appendix B primes
├── benches/                # criterion: mul/sqr by degree, powmod, full test by digits
└── docs/                   # results, plots, comparison with the C++ build
```

CLI sketch, a superset of the corrected C++ programs' command line:

```
aks [--theorem | --phi | --lambda] [n]                      # same flags as the C++; n from stdin if omitted
aks 23456789 --mul kronecker|schoolbook|karatsuba|ntt \
             --threads 8 --log logs/run.csv --trace
```

- Exit status 0 for prime, 1 otherwise, as in the C++ (so the C++ `test/aksSweep.sh` style of check works unchanged).
- The CSV output keeps the C++ column order (`Int, Cores, Prime (T/F), Time (ms), Comments`), with r, φ(r), λ(r) and the range mode in the comments, so the 2021 logs, the corrected C++ and new runs can be plotted together.

---

## 5. C++ → Rust mapping

| C++ (AKS repo) | Rust | Notes |
|---|---|---|
| `devVR/LenstraZnx.cpp`, `LenstraZ.cpp` | `aks::lenstra` + `bin/aks.rs` | One code path; the variant is chosen with `--mul`; same `--theorem/--phi/--lambda` flags and exit status |
| `PerfectPower.h` | `num::is_perfect_power` → `rug` | Exact |
| `Euler.h`, `Carmichael.h`, `Totients.h` | `num::{euler_phi, carmichael_lambda}` | One trial-division factorisation up to √r, as `Totients.h::PhiAndLambda` does |
| r search in `Lenstra*.cpp`, `FindR.h` | `aks::params::find_r` | log₂, a true order test (`OrderExceeds`), u64 arithmetic, and the n ≤ r ⇒ PRIME step |
| `ARange.h` | `aks::params::ARange`, `a_bound` | Same three modes and default |
| `CongruenceZnx.h`, `CongruenceZnxHPC.h`, `Cores.h` | `aks::congruence` | Right-hand side is x^(n mod r) − a with no powermod; rayon `find_first` replaces the atomic minimum and interleaving |
| `CongruenceZ.h`, `biSegMultiplyZZpX.h` | `poly::kronecker` + `CyclicPoly::pow` | Packing as in `ZZpXmultiply`; `ZZpXsqr` → `square()`; `FoldCyclic` → `CyclicPoly::fold` |
| `biSegMultiplyZZX.h` (signed, over Z) | `poly::kronecker::mul_signed` (test-only) | Keeps the corrected sign-bit algorithm as a documented reference |
| `test/unitTest.cpp`, `biSegCorrectnessTest.cpp` | `tests/`, `proptest` | Every multiplier is checked against `Schoolbook`; port the edge cases from `biSegCorrectnessTest.cpp` |
| `test/aksSweep.sh`, `test/probPrime.cpp`, `Makefile` | `tests/` (sweep vs `is_probably_prime`), `cargo test` | — |
| `test/perflog*.h`, `chronoTest.cpp` | `benches/` (criterion) | Replaces hand-timed loops |
| `CUDA/`, `.cu` | M7 (GPU NTT) | See §6 |
| `logs/` CSVs | `docs/results/` | Import the 2021 data for comparison |

---

## 6. Milestones

Each milestone ends with green `cargo test` / `clippy` and a short `docs/` note.

### M0 — Housekeeping and reference data (done)

- [x] `.gitignore` and this plan committed; `num-bigint` removed; `rug` bumped to 1.30 (integer feature only); crate renamed `aks2`, edition 2024.
- [x] GitHub Actions: `macos-latest` (arm64) and `ubuntu-latest`, fmt, clippy, tests, MSRV.
- [x] Licence: MIT OR Apache-2.0, plus `CITATION.cff` (§9).
- [x] **Reference data from the corrected C++** (tag `v1.1`): `tools/cpp_reference.py` runs `LenstraZnx.out --theorem` for every n in 2…10⁴ and the 30 Appendix B primes up to 10 digits, and parses its output into `tests/data/cpp_reference.csv` (columns n, prime, reason, witness, r, phi, lambda, a_max). 10,028 rows, 1,258 primes; every verdict agrees with trial division (`tests/reference_data.rs`). Parsing standard output rather than the CSV log keeps the witness (divisor or failing *a*) and the exact bound printed. Of the rows, 8,646 end at a gcd divisor, 124 are perfect powers, 7 end at "n ≤ r", and 1,249 run the congruence tests.
- [x] **Re-measure the C++ baselines** (§2) on the M3 Pro with `tools/cpp_baselines.sh ../AKS`: 0.23 s and 0.53 s for NTL, 1.03 s and 3.11 s for the Z version, all on 11 threads.

Note for M1: the C++ computes the theorem bound in `double`, as `floor(sqrt(φ) · ln n / ln 2)`. That is exact unless the product falls within rounding error of an integer. The Rust version uses the same `f64` formula, with log₂ n taken from the top 64 bits of n so it works beyond `f64` range, and the reference data checks for parity. If a mismatch ever turns up, it will be investigated, not hidden; a higher-precision bound is an option for M5.

### M1 — Number theory layer (done)

- [x] `src/num.rs`: `log2_floor`, `log2_f64`, `gcd`/`lcm`/`mul_mod`/`pow_mod`, `factor_u64`, `totients` (φ and λ from one factorisation), `mult_order` (works down from λ(r)), `order_exceeds`, `perfect_power` (exact, with the witness) and `is_perfect_power`.
- [x] `src/aks/params.rs`: `ARange` (with the C++ wording and `is_certified`), `Params { r, order, phi, lambda }`, `find_r` → `RSearch::{Found, Divisor, NoRBelowN}`, and `a_bound`.
- [x] `find_r` only computes an order when r > ⌊log₂ n⌋² + 1 and λ(r) > ⌊log₂ n⌋², since ord_r(n) ≤ λ(r) < r. Same answer as the C++, far fewer multiplications: the 42-digit prime's r = 19051 takes milliseconds.
- [x] Tests:
  - Unit tests against brute force: φ and λ for r ≤ 1000, orders for r ≤ 300, perfect powers for n < 20,000, plus bases beyond 2⁶⁴.
  - φ and λ for r ≤ 20,000 against SymPy (`tests/data/totients.csv`, from `tools/totients.py`). The plan said 10⁵; 20,000 already covers every r used for n up to 42 digits and keeps the file to 300 KB.
  - `find_r`: 10007 → 173, 23456789 → 599, 1000000007 → 857, and the 42-digit prime 689960931088884849033689023336009222695077 → 19051.
  - Every one of the 10,028 rows of `tests/data/cpp_reference.csv` agrees: perfect powers, divisor found, "n ≤ r", and for the 1,249 congruence rows r, φ(r), λ(r) and the theorem bound (so the `f64` bound matches the C++ bit for bit on these inputs).

### M2 — Cyclic polynomial ring and multipliers (done)

- [x] `Ring` / `CyclicPoly` (`src/poly/mod.rs`), with construction that folds and reduces, add/sub, lazy fold-then-reduce, and left-to-right `pow`. The `PolyMul` trait has a default `square`.
- [x] `Schoolbook` and `Kronecker` (mul plus a dedicated square; unsigned byte-aligned slots, 2ᵇ > min(U,V)·(n−1)²).
- [x] Property tests (`tests/poly_props.rs`): `Kronecker == Schoolbook` for n up to 512 bits, including 2ᵏ ± 1, with the edge cases (zero, constants, all-ones, all-(n−1), sparse), plus the ring laws and `pow`. The debug suite uses r ≤ 64. A release-mode `#[ignore]` test goes to r ≤ 4096: Schoolbook at r = 20,000 would take minutes per case, so the planned r up to 20,000 is covered by the AKS runs in M3 instead.
- [x] Benchmark (`benches/poly.rs`, results in `docs/M2-poly.md`): Kronecker squaring about doubles per doubling of r (8.8 ms at r = 19,051 for a 10-digit n), against 4× for Schoolbook.

### M3 — Serial Lenstra AKS, both variants (done)

- [x] `ARange::{Theorem, PhiMinusOne, LambdaMinusOne}`, defaulting to `Theorem` = ⌊√φ(r)·log₂ n⌋, capped below n (from M1).
- [x] `aks()` returns a `Trace`: the verdict, `Reason` (small / perfect power with witness / divisor / no r below n / congruences hold / failing *a*), `Params`, a_max, the elapsed time, and `certified()`. n < 2 gives `NotApplicable`. The "Z" and "Zn[x]" versions are now the same driver with a different `Multiplier`.
- [x] (x − a)ⁿ uses an O(r) multiply by x − a for every multiply step; only the squarings go through `PolyMul`.
- [x] The `aks` binary: the C++ flags, n from stdin, the same exit codes, `--mul`, `--quiet`, and `--log` in the C++ CSV format. Uncertified PRIME verdicts say "(not certified)".
- [x] Tests (`tests/aks_reference.rs`, `tests/aks_oracle.rs`; see `docs/M3-aks.md`):
  - All 10,028 C++ reference rows agree on verdict, reason, witness and parameters, including every Appendix B prime up to 10 digits.
  - Sweep against `is_probably_prime(30)`: n ≤ 3000 in debug, 20,000 in the release CI job, and 10⁵ opt-in. The plan's 10⁶ in CI would take hours even with M4, so it is not run.
  - Composite suites: semiprimes above r, a Chernick Carmichael number above r (both must fail a congruence, with the same *a* for every multiplier), and prime powers with bases above 2⁶⁴.
  - Long-running tests follow one convention: `#[ignore]` means the release CI job runs it, and a `slow_` prefix means opt-in only.
- Serial speed: 23456789 takes 2.1 s and 1000000007 takes 6.5 s on one thread, against 1.4 s and 3.5 s for the C++ NTL build on two threads.

### M4 — Parallel congruence testing (done)

- [x] `rayon` `find_first` over 1 ≤ a ≤ bound (`Congruence::first_failure_par`), with `-j/--threads` (`Options::threads`); the default is rayon's global pool. The CSV "Cores" column records the threads used.
- [x] Determinism test: on small r where many *a* fail, the parallel search reports the same smallest failing *a* as the serial one, run after run. The multiplier comparison test pins `threads: Some(1)`.
- [x] Scaling: `tools/scaling.sh` (`docs/M4-parallel.md`).
  - On the M3 Pro, 1000000007 scales at 92–93% efficiency up to 4 threads, all on performance cores. Efficiency falls to 72% at 8 threads and 59% at 11, once the efficiency cores join in: 1.11 s on 1 thread, 0.17 s on 11 (6.4×).
  - With all 11 threads, Rust is **3.1× faster than the C++ NTL build** (0.17 s against 0.53 s) and **18× faster than the C++ Z version** (3.11 s). It already matches NTL's 11-thread time on 2 threads (0.60 s).
  - In the 2-vCPU container it gives 1.7× on 2 threads.
- The opt-in 10⁵ sweep (`slow_sweep_to_100000`) passes: every n ≤ 10⁵ agrees with GMP, in 14 minutes on 2 threads.

### M5 — Performance work (first round done; open-ended)

Measured, not guessed: callgrind showed byte-level import/export and per-coefficient allocation at over 60% of the run, and the GMP squaring at 9% (`docs/M5-performance.md`).

- [x] **Coefficient handling first.** Bit-packed limb slots; folding on limbs, then one reduction per coefficient (word arithmetic for n < 2⁶⁴); `square_into` and an in-place multiply by x − a, so the power loop allocates nothing per coefficient. The M2 version stays as `KroneckerBytes`.
- [x] Specialised squaring (GMP's square path), unchanged from M2.
- [x] **Goal met and exceeded:** with 2 threads, 0.30 s for 23456789, 1.17 s for 1000000007 and 38 s for the 16-digit prime. That is 2.6–4.5× faster than the corrected C++ NTL build and 2.4–6.6× faster than the C++ Z version on the same machine.
- [ ] `Karatsuba` / `Ntt` back-ends: not needed for CPU speed, because GMP already chooses Toom/FFT inside the one big multiply. NTT remains the plan for the GPU (M7).
- [ ] Further ideas, in `docs/M5-performance.md`: Montgomery/Barrett for n ≥ 2⁶⁴, a `Vec<u64>` coefficient representation for small n, and packing straight from limbs.

### M6 — Experiments and write-up (2–3 days)

- Re-run the dissertation's comparisons (Figures 5.3–5.6, A.4–A.12) with **valid parameters**: time against log n on log–log axes, Kronecker against NTT, and λ against φ.
- Publish the tables and plots in `docs/`, with the 2021 numbers alongside for contrast.

### M7 — Stretch: GPU exploration

The 2021 plan of one GPU thread per value of *a*, each doing big-integer polynomial arithmetic, was blocked by the lack of device-side bignum libraries (Ch. 3.4). A better route is to keep the bignums on the CPU and offload the **NTT squaring** to the GPU. NTT over 32/64-bit primes needs only word arithmetic.

- Apple Silicon: `wgpu` compute shaders (Metal backend), or `metal-rs`.
- NVIDIA on Linux: `cudarc`.
- Deliverable: a prototype NTT kernel and a measured break-even degree. If it is a win, a GPU `PolyMul` back-end.

---

## 7. Correctness rules carried over from `REFACTOR.md`

1. ⌊log₂ n⌋ comes from bit length, never from a float logarithm.
2. r is the smallest r with gcd(r, n) = 1 and ord_r(n) > ⌊log₂ n⌋², computed exactly in u64. If some 1 < gcd(a, n) < n for a ≤ r, the answer is COMPOSITE. If n ≤ r, the answer is PRIME.
3. The default range of *a* is the theorem's. φ and λ modes print a "not certified" warning.
4. Perfect powers are detected exactly.
5. The multiplier's size bound uses magnitudes, and slot widths leave room for a sign bit wherever coefficients can be negative.
6. No shared mutable state in parallel code; use `find_first` / `find_map_first`.
7. Every multiplier is property-tested against `Schoolbook`. Every AKS verdict in the test suite is checked against an oracle.

---

## 8. Risks

| Risk | Mitigation |
|---|---|
| `rug`/GMP build friction on a new macOS | Default static build via `gmp-mpfr-sys` (needs only the Xcode CLT); fall back to `use-system-libs` with Homebrew `gmp` |
| Kronecker memory use for large r (the product is an integer of 2·r·b bits; ≈ 1.4 MB at r ≈ 19,051 for the 42-digit prime) | Fine on CPU; reuse buffers |
| Large inputs are slow even when correct (a 42-digit prime needs r ≈ 19,051 and about 19,000 values of *a*; the corrected C++ takes about 100 s for a 16-digit prime) | Long-running tests are opt-in; benchmarks are sampled |
| Rust and C++ disagree on a parameter or verdict | The C++ reference data (M0) pins expected values; investigate any mismatch in both before "fixing" either |
| Copyleft: `rug` is LGPL-3.0+ | Fine to depend on; to distribute binaries, either allow relinking (dynamic GMP) or ship the source, as the LGPL requires |
| Scope creep into GPU work | M7 is explicitly a stretch goal and only starts after M6 |

---

## 9. Open decisions

1. **The C++ repo's licence** is still CC BY 4.0. Switching it to MIT OR Apache-2.0 to match would be consistent; it is the author's call.

Resolved:

- The C++ fixes landed first (AKS `main`, `fef9c34`, tag `v1.1`), so the C++ is the numerical reference for this port.
- Licence: **MIT OR Apache-2.0** (the Rust ecosystem's default; permissive, and both require the copyright notice to be kept, which gives the attribution the author wants), plus `CITATION.cff` for academic citation. Compatible with depending on LGPL `rug`.
- `num-bigint`: **removed**; no milestone needs it.
- Edition **2024**, `rust-version = "1.85"`.
- Crate name **`aks2`** (crate names conventionally avoid dashes, which become underscores in code).

---

## 10. First steps on the M3 Pro

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh   # once
xcode-select --install                                           # once, if not already done
rustup update stable                                             # edition 2024 needs Rust 1.85+
cd AKS2
cargo build            # first build compiles GMP via gmp-mpfr-sys (a few minutes)
cargo test

# reference build of the corrected C++ (for M0's baselines; reference data is already committed)
cd ../AKS && git switch --detach v1.1 && make && make test
cd ../AKS2 && tools/cpp_baselines.sh ../AKS
```
