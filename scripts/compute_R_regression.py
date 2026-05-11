#!/usr/bin/env python3
"""
compute_R_regression.py

For each value of k specified, compute

    R(k) = pi(x; 18k) / [ x / (log x)^2 ]

empirically up to x = xmax, to be compared with the Hardy-Littlewood
prediction

    R_HL(k) = 2 * C_2 * prod_{p odd, p|18k} (p-1)/(p-2)

The ratio R / R_HL must tend to 1 as x -> infinity. The regression
form

    R(k) = C + a/k + O(1/k^2)

applies to the subset of k for which 18k has the same odd-prime
radical (e.g., k = 5, 10, 15, 20 share rad = 30, so the HL factor is
constant); for mixed families R_HL(k) varies with the prime factors
of k.

The output CSV contains one row per (k, x_checkpoint) combination,
with the empirical count, R_emp, R_HL, and their ratio.

Usage:
    python3 compute_R_regression.py --xmax 50_000_000 \\
        --k 5 7 10 11 12 13 15 17 20

Dependencies: numpy, scipy.
"""

from __future__ import annotations
import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad


TWIN_PRIME_CONSTANT = 0.6601618158468695739278121100145557784326233602847334133194484233354056423


def sieve_of_eratosthenes(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return sieve


def odd_prime_factors(n: int) -> list[int]:
    """Returns the distinct odd prime factors of n."""
    factors = []
    # First strip all factors of 2 (we want only odd primes)
    while n % 2 == 0:
        n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        factors.append(n)
    return factors


def hl_constant_for_gap(g: int) -> float:
    odd_factors = odd_prime_factors(g)
    product = 1.0
    for p in odd_factors:
        product *= (p - 1) / (p - 2)
    return 2 * TWIN_PRIME_CONSTANT * product


def li2(x: float, xmin: float = 2.0) -> float:
    """Secondary logarithmic integral li_2(x) = int_{xmin}^x dt/(log t)^2."""
    val, _ = quad(lambda t: 1.0 / (math.log(t) ** 2), xmin, x)
    return val


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xmax",
        type=lambda s: int(s.replace("_", "")),
        default=50_000_000,
        help="Maximum threshold X (default: 5*10^7)",
    )
    parser.add_argument(
        "--k",
        type=int,
        nargs="+",
        default=[5, 7, 10, 11, 12, 13, 15, 17, 20],
        help="Values of k to analyze",
    )
    parser.add_argument(
        "--out",
        default="../data/R_regression.csv",
        help="Output CSV file",
    )
    args = parser.parse_args()

    max_g = 18 * max(args.k)
    sieve_max = args.xmax + max_g
    print(f"[1/3] Sieving primes up to {sieve_max:,}...")
    t0 = time.time()
    sieve = sieve_of_eratosthenes(sieve_max)
    print(f"      done in {time.time() - t0:.1f}s, "
          f"memory: {sieve.nbytes / 1024**2:.1f} MB")

    # Logarithmic checkpoints
    checkpoints = []
    x = 1_000_000
    while x <= args.xmax:
        checkpoints.append(x)
        x *= 2
    if checkpoints[-1] != args.xmax:
        checkpoints.append(args.xmax)

    rows = []
    print(f"\n[2/3] Computing R(k) for k = {args.k}\n")
    print(f"{'k':>3s} {'g':>4s} {'X':>12s} {'count':>10s} "
          f"{'R_emp':>8s} {'R_HL':>8s} {'ratio':>8s}")
    print("-" * 66)

    for k in args.k:
        g = 18 * k
        H = hl_constant_for_gap(g)
        # Identify indices of primes p such that (p, p+g) are both prime
        p_indices = np.arange(5, args.xmax + 1)
        is_pair = sieve[p_indices] & sieve[p_indices + g]
        pair_p = p_indices[is_pair]

        for X in checkpoints:
            count = int(np.sum(pair_p <= X))
            # Empirical R: R_emp = count * (log X)^2 / X
            # (raw asymptotic form, does not use li_2)
            R_emp = count * (math.log(X) ** 2) / X
            # Empirical R via li_2 (more accurate)
            R_emp_li2 = count / li2(X, xmin=5.0)
            R_HL = H
            ratio = R_emp_li2 / R_HL
            row = {
                "k": k,
                "g": g,
                "X": X,
                "count": count,
                "R_emp_naive": R_emp,
                "R_emp_li2": R_emp_li2,
                "R_HL": R_HL,
                "ratio_emp_to_HL": ratio,
                "log_X": math.log(X),
            }
            rows.append(row)
            print(f"{k:>3d} {g:>4d} {X:>12,} {count:>10,} "
                  f"{R_emp_li2:>8.4f} {R_HL:>8.4f} {ratio:>8.4f}")
        print()

    print(f"[3/3] Saving CSV...")
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"      saved to {out_path}")

    # Summary at x_max: R_emp_li2 vs R_HL for each k
    print(f"\nSummary at X = {args.xmax:,}:")
    print(f"{'k':>3s} {'g':>4s} {'R_emp':>8s} {'R_HL':>8s} "
          f"{'ratio':>8s} {'dev_pct':>9s}")
    print("-" * 50)
    for k in args.k:
        last = [r for r in rows if r["k"] == k and r["X"] == args.xmax][0]
        dev = (last["ratio_emp_to_HL"] - 1) * 100
        print(f"{k:>3d} {last['g']:>4d} "
              f"{last['R_emp_li2']:>8.4f} {last['R_HL']:>8.4f} "
              f"{last['ratio_emp_to_HL']:>8.4f} {dev:>+8.2f}%")


if __name__ == "__main__":
    main()
