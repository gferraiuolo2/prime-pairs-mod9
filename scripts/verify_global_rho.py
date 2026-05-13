#!/usr/bin/env python3
"""
verify_global_rho.py

Computes the aggregate ratio
    rho_glob(X) := pi_2(X) / [ H(2) * li_2(X) ],
with
    H(2)   = 2 * C_2   (Hardy-Littlewood constant for gap g = 2),
    li_2(X) = int_{x0}^{X} dt / (log t)^2,
for X = 10^n, n = 7..18.

Empirical pi_2(X) values are taken from OEIS A007508 (Number of twin
prime pairs below 10^n), tabulated in Oliveira e Silva,
"Tables of values of pi(x) and of pi2(x)",
https://sweet.ua.pt/tos/primes.html

The script reproduces Table "tab:rho_global" of the paper.

Usage:
    python3 verify_global_rho.py [--x0 2.0] [--out global_rho.csv]

Dependencies: scipy.
"""

from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

from scipy.integrate import quad


# Twin-prime constant C_2 (Wrench 1961; OEIS A005597), 25 digits.
C2 = 0.6601618158468695739278121
H2 = 2 * C2   # ~ 1.3203236316937391478556242


# OEIS A007508: pi_2(10^n) INCLUDING the pair (3,5).
# n -> pi_2(10^n)
PI2_OEIS: dict[int, int] = {
    7:                 58_980,
    8:                440_312,
    9:              3_424_506,
   10:             27_412_679,
   11:            224_376_048,
   12:          1_870_585_220,
   13:         15_834_664_872,
   14:        135_780_321_665,
   15:      1_177_209_242_304,
   16:     10_304_195_697_298,
   17:     90_948_839_353_159,
   18:    808_675_888_577_436,
}


def li2(x: float, x0: float = 2.0) -> float:
    """Secondary logarithmic integral li_2(x) = int_{x0}^x dt/(log t)^2."""
    val, _ = quad(
        lambda t: 1.0 / (math.log(t) ** 2),
        x0, x,
        limit=200, epsabs=1e-3, epsrel=1e-14,
    )
    return val


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--x0", type=float, default=2.0,
        help="Lower bound for li_2 (default: 2.0)",
    )
    parser.add_argument(
        "--out", default="../data/global_rho.csv",
        help="Output CSV file (default: ../data/global_rho.csv, "
             "matching the layout of the other compute_* scripts)",
    )
    args = parser.parse_args()

    print(f"H(2) = 2 * C_2 = {H2:.16f}")
    print(f"li_2 lower bound x0 = {args.x0}\n")

    print(f"{'n':>3s} {'X = 10^n':>22s} {'pi_2(X) (OEIS)':>22s} "
          f"{'H(2)*li_2(X)':>22s} {'rho':>10s} {'(rho-1) ppm':>14s}")
    print("-" * 100)

    rows = []
    for n in sorted(PI2_OEIS):
        X = 10 ** n
        pi2 = PI2_OEIS[n]
        L = li2(X, x0=args.x0)
        pred = H2 * L
        rho = pi2 / pred
        ppm = (rho - 1) * 1e6
        print(f"{n:>3d} {X:>22,d} {pi2:>22,d} {pred:>22,.1f} "
              f"{rho:>10.7f} {ppm:>+13.1f}")
        rows.append({
            "n": n,
            "X": X,
            "pi2_OEIS": pi2,
            "li2": L,
            "H2_times_li2": pred,
            "rho_glob": rho,
            "rho_minus_1_ppm": ppm,
        })

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nCSV saved to: {out_path}")


if __name__ == "__main__":
    main()
