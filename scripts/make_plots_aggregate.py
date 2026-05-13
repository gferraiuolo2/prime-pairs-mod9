#!/usr/bin/env python3
"""
make_plots_aggregate.py

Renders Figure 2 of the paper: log-log decay of |rho_glob(X) - 1|
for twin primes, using OEIS A007508 data tabulated in
../data/global_rho.csv (produced by verify_global_rho.py).

Output: ../figures/decadimento_aggregato.{pdf,png}

The figure shows the empirical decay over 11 orders of magnitude
(X = 10^7 .. 10^18), with a dashed guideline proportional to
X^{-1/2} * log X normalised at X = 10^10, as cited in §3.2 of the
paper. The slope quantification (-0.45 +- 0.03) is reported by
analysis.py.

Usage:
    python3 make_plots_aggregate.py

Dependencies: numpy, matplotlib.
"""

from __future__ import annotations
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rcParams.update({
    "font.family":  "serif",
    "font.serif":   ["Times New Roman", "DejaVu Serif"],
    "font.size":    10,
    "axes.labelsize":  11,
    "axes.titlesize":  11,
    "legend.fontsize": 9,
    "figure.figsize":  (6, 4),
    "lines.linewidth": 1.4,
    "axes.grid":     True,
    "grid.alpha":    0.3,
    "grid.linestyle": "--",
})


def main() -> None:
    base = Path(__file__).resolve().parent.parent
    data_dir = base / "data"
    out_dir  = base / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(data_dir / "global_rho.csv") as f:
        for r in csv.DictReader(f):
            r_X   = float(r["X"])
            r_rho = float(r["rho_glob"])
            rows.append((r_X, abs(r_rho - 1.0)))

    Xs   = [r[0] for r in rows]
    errs = [r[1] for r in rows]

    fig, ax = plt.subplots()
    ax.loglog(Xs, errs, "o-", color="C0", markersize=6,
              label=r"$|\rho_{\mathrm{glob}}(X) - 1|$")

    # Guideline X^{-1/2} log X, normalised at the second point (X=10^10)
    X_guide = np.logspace(8, 18.2, 80)
    X_ref, err_ref = Xs[1], errs[1]
    guide = err_ref * (X_ref / X_guide) ** 0.5 * (np.log(X_guide) / math.log(X_ref))
    ax.loglog(X_guide, guide, "--", color="grey", linewidth=1.1,
              label=r"$\propto X^{-1/2}\log X$  (normalised at $X=10^{10}$)")

    ax.set_xlabel(r"Threshold $X$")
    ax.set_ylabel(r"$|\rho_{\mathrm{glob}}(X) - 1|$")
    ax.set_title(r"Aggregate decay of HL deviation for twin primes")
    ax.set_xlim(5e7, 5e18)
    ax.set_ylim(1e-9, 1e-3)
    ax.legend(loc="upper right", framealpha=0.95)
    ax.text(0.02, 0.02,
            r"$\pi_2(X)$ from OEIS A007508 (Oliveira e Silva)",
            transform=ax.transAxes, fontsize=8, color="dimgrey",
            verticalalignment="bottom")

    fig.tight_layout()
    fig.savefig(out_dir / "decadimento_aggregato.pdf")
    fig.savefig(out_dir / "decadimento_aggregato.png", dpi=150)
    print(f"  wrote decadimento_aggregato.{{pdf,png}} in {out_dir}")


if __name__ == "__main__":
    main()
