#!/usr/bin/env python3
"""
make_plots.py

Generates the paper figures from the CSV files in ../data/.
Output: ../figures/{twins_convergence,dr_equipartition,HL_confirmation}.{pdf,png}

Usage:
    python3 make_plots.py
"""

from __future__ import annotations
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl


# Uniform publication style
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "figure.figsize": (6, 4),
    "lines.linewidth": 1.4,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def load_csv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (ValueError, TypeError):
                pass
    return rows


def plot_twins_convergence(data_dir: Path, out_dir: Path) -> None:
    rows = load_csv(data_dir / "twins_dr9.csv")
    Xs = [r["X"] for r in rows]

    fig, ax = plt.subplots()
    for cls, color, marker in [(2, "C0", "o"), (5, "C1", "s"), (8, "C2", "^")]:
        deltas = [r[f"delta_dr{cls}_pp"] for r in rows]
        ax.plot(Xs, deltas, color=color, marker=marker,
                markersize=4, label=f"$a={cls}$")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"Threshold $X$")
    ax.set_ylabel(r"$\pi_2(X;9,a)/\pi_2(X) - 1/3$  [percentage points]")
    ax.set_title("Convergence of DR fractions for twin primes")
    ax.legend(title="class", loc="best", framealpha=0.9)
    ax.set_ylim(-0.4, 0.4)

    fig.tight_layout()
    # Keep the same filename used in the LaTeX source
    fig.savefig(out_dir / "convergenza_gemelli.pdf")
    fig.savefig(out_dir / "convergenza_gemelli.png", dpi=150)
    plt.close(fig)
    print(f"  wrote convergenza_gemelli.{{pdf,png}}")


def plot_equipartition(data_dir: Path, out_dir: Path) -> None:
    rows = load_csv(data_dir / "permutation_families.csv")
    ks = sorted({int(r["k"]) for r in rows})
    drs = [1, 2, 4, 5, 7, 8]

    fig, ax = plt.subplots(figsize=(7, 4))
    width = 0.13
    x_pos = list(range(len(ks)))

    for i, dr in enumerate(drs):
        deltas = []
        for k in ks:
            r = next(r for r in rows if int(r["k"]) == k and int(r["dr"]) == dr)
            deltas.append(r["delta_uniform_pp"])
        offsets = [x + (i - 2.5) * width for x in x_pos]
        ax.bar(offsets, deltas, width=width, label=f"$a={dr}$")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"$k={k}$" for k in ks])
    ax.set_ylabel(r"DR fraction $-1/6$  [pp]")
    ax.set_title("Equidistribution of DR classes for pairs $(p,p+18k)$")
    ax.legend(title="class", loc="best", ncol=3, framealpha=0.9)
    ax.set_ylim(-0.25, 0.25)

    fig.tight_layout()
    fig.savefig(out_dir / "equipartizione_dr.pdf")
    fig.savefig(out_dir / "equipartizione_dr.png", dpi=150)
    plt.close(fig)
    print(f"  wrote equipartizione_dr.{{pdf,png}}")


def plot_HL_confirmation(data_dir: Path, out_dir: Path) -> None:
    rows = load_csv(data_dir / "R_regression.csv")
    ks = sorted({int(r["k"]) for r in rows})
    Xmax = max(r["X"] for r in rows)

    fig, ax = plt.subplots()
    ratios_at_xmax = []
    for k in ks:
        r = next(r for r in rows if int(r["k"]) == k and r["X"] == Xmax)
        ratios_at_xmax.append((k, r["ratio_emp_to_HL"]))

    ks_arr = [t[0] for t in ratios_at_xmax]
    ratios = [t[1] for t in ratios_at_xmax]
    ax.plot(ks_arr, ratios, "o-", color="C0", markersize=6)
    ax.axhline(1, color="black", linewidth=0.8, linestyle="--",
               label="exact HL")
    ax.axhspan(0.99, 1.01, color="C0", alpha=0.15,
               label=r"$\pm 1\%$")

    ax.set_xlabel("$k$")
    ax.set_ylabel(r"$R_{\mathrm{emp}}(k) / R_{\mathrm{HL}}(k)$")
    ax.set_title(f"Hardy--Littlewood confirmation at $X = {int(Xmax):,}$")
    ax.legend(loc="best", framealpha=0.9)
    ax.set_ylim(0.985, 1.015)

    fig.tight_layout()
    fig.savefig(out_dir / "conferma_HL.pdf")
    fig.savefig(out_dir / "conferma_HL.png", dpi=150)
    plt.close(fig)
    print(f"  wrote conferma_HL.{{pdf,png}}")


def main() -> None:
    base = Path(__file__).parent.parent
    data_dir = base / "data"
    out_dir = base / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Generating figures:")
    plot_twins_convergence(data_dir, out_dir)
    plot_equipartition(data_dir, out_dir)
    plot_HL_confirmation(data_dir, out_dir)
    print(f"\nFigures saved in {out_dir}")


if __name__ == "__main__":
    main()
