#!/usr/bin/env python3
"""
compute_permutation_families.py

Counts prime pairs (p, p+g) with g = 18k for user-specified k values,
classified by the digital root dr(p) = p mod 9. Since g is divisible
by 9, we have dr(p+g) = dr(p), and the 6 admissible unit classes are
{1, 2, 4, 5, 7, 8}.

For each k it returns:
  - empirical count per DR class
  - empirical fraction per DR class
  - Hardy-Littlewood prediction (singular series) per DR class
  - relative deviation (empirical - HL) / HL

Usage:
    python3 compute_permutation_families.py --xmax 10_000_000 \
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


# Twin-prime constant C_2 = prod_{p>2} (1 - 1/(p-1)^2)
# High-precision numerical value (Wrench, 1961; OEIS A005597)
TWIN_PRIME_CONSTANT = 0.6601618158468695739278121100145557784326233602847334133194484233354056423


def sieve_of_eratosthenes(n: int) -> np.ndarray:
    """In-memory segmented sieve of Eratosthenes. For n ~ 10^9 needs ~125 MB."""
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
    # n is now odd. Trial division by odd numbers.
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
    """
    Hardy-Littlewood constant for the pair (p, p+g):
        H(g) = 2 * C_2 * prod_{p odd, p|g} (p-1)/(p-2)
    such that pi(x; g) ~ H(g) * x / (log x)^2.
    """
    odd_factors = odd_prime_factors(g)
    product = 1.0
    for p in odd_factors:
        product *= (p - 1) / (p - 2)
    return 2 * TWIN_PRIME_CONSTANT * product


def hl_prediction_per_dr(g: int, x_max: int, x_min: int = 2) -> float:
    """
    Asymptotic Hardy-Littlewood prediction restricted to one arithmetic
    progression modulo 9, for a pair (p, p+g) with g divisible by 9.

    For g divisible by 9 the admissible DR classes are the 6 units of
    Z/9Z = {1,2,4,5,7,8}, and each receives 1/6 of the total density
    (by symmetry of the local factor at p=3).

    Uses the secondary logarithmic integral li_2(x) = int_2^x dt/(log t)^2
    for a more accurate estimate than x/(log x)^2.
    """
    H = hl_constant_for_gap(g)
    # Approximation: integral of 1/(log t)^2 from x_min to x_max
    # via adaptive quadrature
    from scipy.integrate import quad

    integral, _ = quad(lambda t: 1.0 / (math.log(t) ** 2), x_min, x_max)
    return H * integral / 6.0


def count_pairs_by_dr(
    xmax: int, k: int, sieve: np.ndarray
) -> dict:
    """
    For a given k, count the pairs (p, p+18k) with p in [5, xmax-18k]
    and p, p+18k both prime, classified by dr(p) mod 9.
    """
    g = 18 * k
    if xmax + g >= len(sieve):
        raise ValueError(f"Sieve too small: need sieve up to {xmax + g}")

    p_indices = np.arange(5, xmax + 1)
    is_pair = sieve[p_indices] & sieve[p_indices + g]
    pair_p = p_indices[is_pair]
    dr = pair_p % 9

    total = len(pair_p)
    counts = {a: int(np.sum(dr == a)) for a in [1, 2, 4, 5, 7, 8]}
    fractions = {a: counts[a] / total if total > 0 else 0.0 for a in counts}

    # Sanity check: no pairs in inadmissible classes
    non_admissible = total - sum(counts.values())
    assert non_admissible == 0, f"k={k}: found {non_admissible} inadmissible pairs"

    return {
        "k": k,
        "g": g,
        "total": total,
        "counts": counts,
        "fractions": fractions,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xmax",
        type=lambda s: int(s.replace("_", "")),
        default=10_000_000,
        help="Maximum threshold X for p (default: 10^7)",
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
        default="../data/permutation_families.csv",
        help="Output CSV file",
    )
    args = parser.parse_args()

    max_g = 18 * max(args.k)
    sieve_max = args.xmax + max_g
    print(f"[1/3] Sieving primes up to {sieve_max:,}...")
    t0 = time.time()
    sieve = sieve_of_eratosthenes(sieve_max)
    print(f"      done in {time.time() - t0:.1f}s, "
          f"sieve memory: {sieve.nbytes / 1024**2:.1f} MB")

    results = []
    print(f"\n[2/3] Counting pairs for k = {args.k}")
    print(f"      Twin constant C_2 = {TWIN_PRIME_CONSTANT:.10f}\n")

    for k in args.k:
        t0 = time.time()
        r = count_pairs_by_dr(args.xmax, k, sieve)
        H = hl_constant_for_gap(r["g"])
        # HL prediction per DR class (1/6 of the total)
        from scipy.integrate import quad
        integral, _ = quad(lambda t: 1.0 / (math.log(t) ** 2), 5, args.xmax)
        pred_per_class = H * integral / 6.0

        print(f"k={k:2d} (g={r['g']:3d}): total {r['total']:>7,} pairs, "
              f"H(g)={H:.4f}, predicted/class={pred_per_class:.0f}, "
              f"computed in {time.time() - t0:.1f}s")
        for a in [1, 2, 4, 5, 7, 8]:
            c = r["counts"][a]
            f = r["fractions"][a] * 100
            dev_vs_uniform = (r["fractions"][a] - 1/6) * 100
            dev_vs_hl = (c - pred_per_class) / pred_per_class * 100
            row = {
                "k": k,
                "g": r["g"],
                "dr": a,
                "count": c,
                "fraction_pct": f,
                "delta_uniform_pp": dev_vs_uniform,
                "hl_predicted": pred_per_class,
                "deviation_vs_hl_pct": dev_vs_hl,
                "H_g": H,
                "xmax": args.xmax,
            }
            results.append(row)
            print(f"     dr={a}: {c:>6,} ({f:6.3f}%, "
                  f"Delta_unif={dev_vs_uniform:+.3f} pp, "
                  f"dev/HL={dev_vs_hl:+6.2f}%)")
        print()

    print(f"[3/3] Saving CSV...")
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(f"      saved to {out_path}")


if __name__ == "__main__":
    main()
