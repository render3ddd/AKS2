#!/usr/bin/env python3
"""Export reference data from the corrected C++ AKS (render3d/AKS, tag v1.1).

Runs `LenstraZnx.out --theorem n` for every n in 2..=LIMIT plus the Appendix B
primes of the dissertation up to 10 digits, parses what the program prints, and
writes one CSV row per n. The Rust test suite (milestones M1 and M3) checks its
parameters and verdicts against this file.

Usage (from anywhere; needs the C++ repo built with `make`):

    python3 tools/cpp_reference.py --bin ../AKS/devVR/LenstraZnx.out \
        --out tests/data/cpp_reference.csv [--limit 10000] [--jobs 8]

Columns
    n        the input
    prime    1 if the C++ said prime (exit status 0), else 0
    reason   small          n = 2 or 3, answered directly
             perfect_power  n = a^b, b > 1
             divisor        1 < gcd(r', n) < n for some r' searched before r
             n_le_r         no r < n with ord_r(n) > floor(log2 n)^2: PRIME
             congruence     all congruences checked (prime) or one failed
    witness  the divisor (reason=divisor) or the failing a (congruence, composite)
    r, phi, lambda, a_max
             the AKS parameters: r, phi(r), lambda(r) and the theorem's upper
             bound on a, floor(sqrt(phi(r)) * log2 n) capped at n - 1.
             Empty when the run ended before r was found.
"""

import argparse
import concurrent.futures as cf
import csv
import re
import subprocess
import sys

APPENDIX_B = [
    347, 11701, 23456789, 1000000007, 10133, 32561, 160583, 268069, 2969023,
    7741009, 12434839, 86457079, 100982933, 879673117, 2492813497, 88903,
    99397, 107339, 777619, 1000099, 4740623, 26933611, 10012333, 100041703,
    103736293, 1001772091, 2958347047, 104729, 200000033, 9999954997,
]

RE_R = re.compile(r"^r = (\d+)$", re.M)
RE_PL = re.compile(r"^Phi\((\d+)\) = (\d+), Lambda\(\d+\) = (\d+)$", re.M)
RE_A = re.compile(r"^Testing 1 <= a <= (\d+) \(theorem", re.M)
RE_DIV = re.compile(r"^(\d+) is a divisor\.$", re.M)
RE_FAIL = re.compile(r"^The a which fails is (\d+)$", re.M)


def run(binary, n):
    p = subprocess.run([binary, "--theorem", str(n)], capture_output=True,
                       text=True, check=False)
    out = p.stdout
    prime = {0: 1, 1: 0}.get(p.returncode)
    if prime is None:
        raise RuntimeError(f"n={n}: exit status {p.returncode}\n{out}{p.stderr}")
    row = {"n": n, "prime": prime, "reason": "", "witness": "",
           "r": "", "phi": "", "lambda": "", "a_max": ""}
    m_r, m_pl, m_a = RE_R.search(out), RE_PL.search(out), RE_A.search(out)
    if m_r:
        row["r"] = m_r.group(1)
    if m_pl:
        row["phi"], row["lambda"] = m_pl.group(2), m_pl.group(3)
    if m_a:
        row["a_max"] = m_a.group(1)

    if n in (2, 3):
        row["reason"] = "small"
    elif "is a perfect power" in out:
        row["reason"] = "perfect_power"
    elif m := RE_DIV.search(out):
        row["reason"], row["witness"] = "divisor", m.group(1)
    elif m_r:
        row["reason"] = "congruence"
        if m := RE_FAIL.search(out):
            row["witness"] = m.group(1)
    else:
        row["reason"] = "n_le_r"

    # Sanity: every congruence row carries the full parameter set.
    if row["reason"] == "congruence" and not (row["phi"] and row["a_max"]):
        raise RuntimeError(f"n={n}: could not parse parameters\n{out}")
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bin", required=True, help="path to LenstraZnx.out")
    ap.add_argument("--out", required=True, help="CSV file to write")
    ap.add_argument("--limit", type=int, default=10_000)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    ns = sorted(set(range(2, args.limit + 1)) | set(APPENDIX_B))
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        rows = list(ex.map(lambda n: run(args.bin, n), ns))

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    primes = sum(r["prime"] for r in rows)
    print(f"wrote {len(rows)} rows ({primes} prime) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
