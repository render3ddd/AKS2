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

Linux x86-64 container with 2 vCPUs:

| n | 1 thread | 2 threads | speed-up |
|---|---:|---:|---:|
| 23456789 | 2.12 s | 1.26 s | 1.68× |
| 1000000007 | 6.39 s | 3.71 s | 1.72× |

The two vCPUs are probably hyperthreads of one core, which would explain the efficiency below 2×. On the M3 Pro, `tools/scaling.sh 1000000007` should scale almost linearly up to the performance-core count. After that the efficiency cores join and the curve flattens. **M3 Pro figures to be added by the author.**

With 2 threads the Rust build matches the corrected C++ NTL build on the same machine (1.4 s and 3.5 s) and beats the C++ Z version (2.0 s and 5.9 s). The M5 multiplier then makes it roughly 3× faster again.
