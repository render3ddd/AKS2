# M4: parallel congruence tests

The congruences for different *a* are independent, so `Congruence::first_failure_par` hands 1 ≤ a ≤ bound to rayon's `find_first`. Two properties matter for AKS:

- **The smallest failing *a* is reported**, whatever order the threads finish in, so a composite's witness is the same as in a serial run. This is tested repeatedly on deliberately small r, where many values of *a* fail.
- **Early termination.** Once some *a* fails, work on larger *a* stops. Each test owns its polynomials, so nothing is shared mutably, unlike the 2021 `CongruenceZnxHPC.h` (REFACTOR.md §1.6).

`aks -j N` (or `Options::threads`) selects the thread count:

- the default is one thread per core, or `RAYON_NUM_THREADS`;
- `-j 1` runs serially on the calling thread;
- the CSV log records the count in the C++ "Cores" column.

## Scaling

`tools/scaling.sh [n] [max-threads]`, best of three runs.

### M3 Pro (5 performance + 6 efficiency cores, 11 threads)

n = 1000000007, measured after M5, so with the tuned multiplier:

| threads | time | speed-up | efficiency |
|---:|---:|---:|---:|
| 1 | 1.108 s | 1.00× | 100% |
| 2 | 0.595 s | 1.86× | 93% |
| 4 | 0.300 s | 3.69× | 92% |
| 8 | 0.192 s | 5.77× | 72% |
| 11 | 0.172 s | 6.44× | 59% |

The curve shows the core types: near-linear while every thread has a performance core (up to 5), then flatter as efficiency cores are added. Each efficiency core adds between a quarter and two-fifths of a performance core's throughput here.

Against the corrected C++ on the same Mac (`tools/cpp_baselines.sh`, which always uses all 11 threads):

| n = 1000000007 | time | Rust is |
|---|---:|---:|
| C++ Zn[x] (NTL), 11 threads | 0.525 s | 3.1× faster |
| C++ Z (binary segmentation), 11 threads | 3.109 s | 18× faster |
| Rust, 11 threads | 0.172 s | — |
| Rust, 2 threads | 0.595 s | about level with NTL on 11 |

### Linux x86-64 container with 2 vCPUs

Measured before M5:

| n | 1 thread | 2 threads | speed-up |
|---|---:|---:|---:|
| 23456789 | 2.12 s | 1.26 s | 1.68× |
| 1000000007 | 6.39 s | 3.71 s | 1.72× |

The two vCPUs are probably hyperthreads of one core, which would explain the efficiency below 2×.

With 2 threads the Rust build matches the corrected C++ NTL build on the same machine (1.4 s and 3.5 s) and beats the C++ Z version (2.0 s and 5.9 s). The M5 multiplier then makes it roughly 3× faster again.
