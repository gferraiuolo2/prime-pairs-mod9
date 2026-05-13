#!/usr/bin/env python3
"""
analysis.py

Two statistical analyses on the empirical data, cited in §3.1 and §3.2
of the paper:

  (A) Power-law fit on the aggregate ratio
        rho_glob(X) = pi_2(X) / [2*C_2 * li_2(X)]
      for twin primes, using OEIS A007508 values tabulated in
      ../data/global_rho.csv (produced by verify_global_rho.py).

      Three models are fitted on n = 8..18 (the asymptotic regime):
        (A1) log|rho-1|             = a + b log X
        (A2) log|rho-1|             = a + b log X + c log(log X)
        (A3) log(|rho-1| / log X)   = a + b log X

  (B) Chi-square goodness-of-fit on the DR-mod-9 partition of twin
      primes at every logarithmic checkpoint in ../data/twins_dr9.csv
      (produced by compute_twins_dr9.py).

Usage:
    python3 analysis.py

The script prints results to stdout; no files are written. The numbers
reported here are quoted verbatim in §3.1 and §3.2 of paper_v6_en.tex.

Dependencies: numpy, scipy.
"""

from __future__ import annotations
import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"


def main() -> None:
    # ------------------------------------------------------------------
    # (A) Power-law fit on aggregate ratio
    # ------------------------------------------------------------------
    print("=" * 72)
    print("(A) Slope fit:  log|rho_glob - 1|  vs  log X")
    print("=" * 72)

    rows = []
    with open(DATA / "global_rho.csv") as f:
        for r in csv.DictReader(f):
            n = int(r["n"])
            X = float(r["X"])
            rho = float(r["rho_glob"])
            rows.append((n, X, rho, abs(rho - 1)))

    print(f"\n{'n':>3s} {'X':>22s} {'|rho-1|':>14s}")
    for n, X, _rho, e in rows:
        print(f"{n:>3d} {X:>22,.0f} {e:>14.3e}")

    data_asymp = [(X, e) for (n, X, _r, e) in rows if n >= 8 and e > 0]
    logX = np.array([math.log(X) for X, _ in data_asymp])
    logE = np.array([math.log(e) for _, e in data_asymp])
    loglogX = np.log(logX)

    print(f"\nFit on n = 8..18 ({len(data_asymp)} points)")
    print("-" * 72)

    s1, i1, r1, p1, se1 = stats.linregress(logX, logE)
    print(f"\n(A1) log|rho-1| = a + b log X")
    print(f"     slope b = {s1:+.4f}  (SE {se1:.4f}, 95% CI "
          f"[{s1 - 1.96 * se1:+.4f}, {s1 + 1.96 * se1:+.4f}])")
    print(f"     intercept a = {i1:+.3f}")
    print(f"     R^2 = {r1 ** 2:.5f}")
    print(f"     distance from theoretical -0.5: "
          f"{(s1 - (-0.5)) / se1:+.2f} sigma")

    M = np.column_stack([np.ones_like(logX), logX, loglogX])
    coef, *_ = np.linalg.lstsq(M, logE, rcond=None)
    a2, b2, c2 = coef
    pred2 = M @ coef
    ss_res = float(np.sum((logE - pred2) ** 2))
    ss_tot = float(np.sum((logE - logE.mean()) ** 2))
    r2_2 = 1 - ss_res / ss_tot
    dof = len(logE) - 3
    mse = ss_res / dof
    cov = mse * np.linalg.inv(M.T @ M)
    se_b2 = math.sqrt(cov[1, 1])
    se_c2 = math.sqrt(cov[2, 2])
    print(f"\n(A2) log|rho-1| = a + b log X + c log(log X)")
    print(f"     slope b = {b2:+.4f}  (SE {se_b2:.4f})")
    print(f"     coef c  = {c2:+.4f}  (SE {se_c2:.4f})")
    print(f"     R^2 = {r2_2:.5f}")
    print(f"     (theoretical X^(-1/2) log X => b = -0.5, c = +1)")

    logE_minus = logE - loglogX
    s3, i3, r3, p3, se3 = stats.linregress(logX, logE_minus)
    print(f"\n(A3) log(|rho-1| / log X) = a + b log X")
    print(f"     slope b = {s3:+.4f}  (SE {se3:.4f}, 95% CI "
          f"[{s3 - 1.96 * se3:+.4f}, {s3 + 1.96 * se3:+.4f}])")
    print(f"     intercept a = {i3:+.3f}")
    print(f"     R^2 = {r3 ** 2:.5f}")
    print(f"     distance from theoretical -0.5: "
          f"{(s3 - (-0.5)) / se3:+.2f} sigma")

    # ------------------------------------------------------------------
    # (B) Chi-square goodness-of-fit on DR classes
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("(B) Chi-square test: DR mod 9 equidistribution for twin primes")
    print("=" * 72)

    twin_rows = []
    with open(DATA / "twins_dr9.csv") as f:
        for r in csv.DictReader(f):
            twin_rows.append({
                "X":     int(r["X"]),
                "total": int(r["total"]),
                "c2":    int(r["count_dr2"]),
                "c5":    int(r["count_dr5"]),
                "c8":    int(r["count_dr8"]),
            })

    print(f"\n{'X':>14s} {'N':>10s} {'c=2':>10s} {'c=5':>10s} {'c=8':>10s} "
          f"{'chi2':>8s} {'p-value':>10s}")
    print("-" * 76)
    for r in twin_rows:
        obs = np.array([r["c2"], r["c5"], r["c8"]])
        N = r["total"]
        exp = np.array([N / 3.0] * 3)
        chi2 = float(np.sum((obs - exp) ** 2 / exp))
        pval = stats.chi2.sf(chi2, df=2)
        print(f"{r['X']:>14,d} {N:>10,d} {obs[0]:>10,d} {obs[1]:>10,d} "
              f"{obs[2]:>10,d} {chi2:>8.3f} {pval:>10.4f}")

    print("\n" + "-" * 76)
    print("Interpretation at the three main thresholds (10^7, 10^8, 10^9)")
    print("-" * 76)
    crit_05 = stats.chi2.ppf(0.95, df=2)
    for r in twin_rows:
        if r["X"] in {10_000_000, 100_000_000, 1_000_000_000}:
            obs = np.array([r["c2"], r["c5"], r["c8"]])
            N = r["total"]
            exp = np.array([N / 3.0] * 3)
            chi2 = float(np.sum((obs - exp) ** 2 / exp))
            pval = stats.chi2.sf(chi2, df=2)
            verdict = "PASS" if chi2 < crit_05 else "REJECT (0.05)"
            print(f"\nX = 10^{int(math.log10(r['X']))}: "
                  f"chi^2 = {chi2:.3f}, df=2, p = {pval:.4f}, "
                  f"critical(0.05) = {crit_05:.2f}  =>  {verdict}")


if __name__ == "__main__":
    main()
