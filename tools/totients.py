#!/usr/bin/env python3
"""Write tests/data/totients.csv: r, phi(r), lambda(r) for 1 <= r <= LIMIT, from SymPy.

    python3 tools/totients.py [--limit 20000] > tests/data/totients.csv

An independent check on aks2::num::totients. 20,000 covers every r that AKS
picks for n up to about 42 digits (r = 19051 for the dissertation's largest prime).
"""
import argparse
try:  # SymPy >= 1.13
    from sympy.functions.combinatorial.numbers import reduced_totient, totient
except ImportError:
    from sympy.ntheory import reduced_totient, totient

ap = argparse.ArgumentParser()
ap.add_argument("--limit", type=int, default=20_000)
args = ap.parse_args()
print("r,phi,lambda")
for r in range(1, args.limit + 1):
    print(f"{r},{totient(r)},{reduced_totient(r)}")
