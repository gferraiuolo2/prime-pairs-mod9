#!/usr/bin/env python3
"""
compute_twins_dr9.py

Counts twin prime pairs (p, p+2) with p <= X, classified by digital
root dr(p) = p mod 9. The only admissible classes (for p > 3) are
a in {2, 5, 8}.

Output: a CSV file with cumulative counts at increasing thresholds,
suitable for plotting the convergence toward 1/3.

Usage:
    python3 compute_twins_dr9.py --xmax 10_000_000
    python3 compute_twins_dr9.py --xmax 100_000_000   # slower

Dependencies: numpy (for the sieve).
"""

from __future__ import annotations
import argparse
import csv
import time
from pathlib import Path

import numpy as np


def sieve_of_eratosthenes(n: int) -> np.ndarray:
    """Returns a boolean array where sieve[i] = True iff i is prime."""
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return sieve


def count_twins_by_dr(xmax: int, checkpoints: list[int]) -> list[dict]:
    """
    For each checkpoint X <= xmax, count the twin prime pairs (p, p+2)
    with p <= X, split by dr(p) mod 9.

    Only p > 3 is considered (the pair (3,5) is excluded as a unique
    modular exception, lying outside the standard unit classes).
    """
    print(f"[1/3] Sieving primes up to {xmax + 2:,}...")
    t0 = time.time()
    sieve = sieve_of_eratosthenes(xmax + 2)
    print(f"      sieve done in {time.time() - t0:.1f}s")

    print(f"[2/3] Identifying twin pairs and classifying by dr...")
    t0 = time.time()

    # Indices of primes p such that (p, p+2) are both prime and p <= xmax.
    # We exclude p=3 (the pair (3,5)) for modular uniformity.
    p_indices = np.arange(5, xmax + 1)  # candidates p >= 5
    is_twin = sieve[p_indices] & sieve[p_indices + 2]
    twin_p = p_indices[is_twin]

    print(f"      found {len(twin_p):,} twin pairs in [5, {xmax:,}]")
    print(f"      classification done in {time.time() - t0:.1f}s")

    print(f"[3/3] Computing counts at checkpoints...")
    dr = twin_p % 9  # dr(p) for primes not divisible by 3 is p mod 9 in {1,2,4,5,7,8}
    # For twins (p, p+2) with p > 3, p mod 3 = 2, hence p mod 9 in {2,5,8}.

    rows: list[dict] = []
    for X in checkpoints:
        mask = twin_p <= X
        sub_dr = dr[mask]
        total = len(sub_dr)
        if total == 0:
            continue
        c2 = int(np.sum(sub_dr == 2))
        c5 = int(np.sum(sub_dr == 5))
        c8 = int(np.sum(sub_dr == 8))
        # Sanity check: no pairs in inadmissible classes
        c_other = total - (c2 + c5 + c8)
        rows.append(
            {
                "X": X,
                "total": total,
                "count_dr2": c2,
                "count_dr5": c5,
                "count_dr8": c8,
                "frac_dr2": c2 / total,
                "frac_dr5": c5 / total,
                "frac_dr8": c8 / total,
                "delta_dr2_pp": (c2 / total - 1 / 3) * 100,
                "delta_dr5_pp": (c5 / total - 1 / 3) * 100,
                "delta_dr8_pp": (c8 / total - 1 / 3) * 100,
                "non_admissible": c_other,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xmax",
        type=lambda s: int(s.replace("_", "")),
        default=10_000_000,
        help="Maximum threshold X (default: 10^7)",
    )
    parser.add_argument(
        "--out",
        default="../data/twins_dr9.csv",
        help="Output CSV file",
    )
    args = parser.parse_args()

    # Logarithmically spaced checkpoints for the convergence plot
    checkpoints = []
    x = 10_000
    while x <= args.xmax:
        checkpoints.append(x)
        x *= 2
    if checkpoints[-1] != args.xmax:
        checkpoints.append(args.xmax)

    rows = count_twins_by_dr(args.xmax, checkpoints)

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nFinal result at X = {args.xmax:,}:")
    last = rows[-1]
    print(f"  total twin pairs: {last['total']:,}")
    print(f"  dr=2: {last['count_dr2']:,} ({last['frac_dr2']*100:.4f}%, "
          f"Delta = {last['delta_dr2_pp']:+.4f} pp)")
    print(f"  dr=5: {last['count_dr5']:,} ({last['frac_dr5']*100:.4f}%, "
          f"Delta = {last['delta_dr5_pp']:+.4f} pp)")
    print(f"  dr=8: {last['count_dr8']:,} ({last['frac_dr8']*100:.4f}%, "
          f"Delta = {last['delta_dr8_pp']:+.4f} pp)")
    print(f"  inadmissible (sanity check, must be 0): {last['non_admissible']}")
    print(f"\nCSV saved to: {out_path}")


if __name__ == "__main__":
    main()
