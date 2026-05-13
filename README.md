# Paper: Modulo 9 classification of prime pairs $(p, p+18k)$

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20117270.svg)](https://doi.org/10.5281/zenodo.20117270)
[![License: MIT](https://img.shields.io/badge/Code-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20%26%20Paper-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**GitHub**: <https://github.com/gferraiuolo2/prime-pairs-mod9>
**Zenodo (Concept DOI)**: <https://doi.org/10.5281/zenodo.20117270>
**Contact**: g.ferraiuolo2@gmail.com

> **Note on DOI.** Always cite the *Concept DOI* `10.5281/zenodo.20117270`
> above, which resolves to the latest version. Each release produces a
> distinct *Version DOI* (e.g. `…20117271` for v1.2.0, etc.); these
> should be used only when a specific snapshot must be pinned.

---

LaTeX source, Python computation scripts and datasets accompanying the
paper "Modulo 9 classification of prime pairs $(p,p+g)$: an empirical
verification of the Hardy–Littlewood conjecture for $g=18k$"
by Giovanni Ferraiuolo.

## What's new in v1.2.1 (May 2026)

- **Twin primes extended to $X = 10^9$** for the DR-mod-9 partition:
  three admissible classes within $\pm 0.025$ pp of $1/3$ (one-$\sigma$
  width), chi-square test $\chi^2 = 0.82$ (df $=2$, $p = 0.67$).
- **Aggregate Hardy–Littlewood benchmark up to $X = 10^{18}$** via
  Oliveira e Silva / OEIS A007508: $|\rho_{\rm glob} - 1| \le
  3\cdot 10^{-8}$ at $X = 10^{18}$. Log–log fit of the deviation gives a
  slope of $-0.45 \pm 0.03$ on $\log(|\rho-1|/\log X)$ vs $\log X$,
  $\approx 2\sigma$ from the theoretical $-1/2$, in line with the
  function-field analogues of Sawin–Shusterman and Keating–Rudnick.
- **Memory bug fix** in `compute_twins_dr9.py`: the previous version
  allocated `np.arange(5, xmax+1)` (7.5 GiB at $X=10^9$). The patched
  version uses `np.flatnonzero(sieve_mask)`, peak memory $\sim 2$ GB.
- **New scripts**: `verify_global_rho.py` (aggregate OEIS benchmark)
  and `analysis.py` (slope fit + chi-square goodness-of-fit).
- **New figure**: `decadimento_aggregato.pdf` (log–log decay over 11
  orders of magnitude).

See `CHANGELOG.md` for the full version history.

## Repository layout

```
prime-pairs-mod9/
├── paper.tex           # LaTeX source (compile with pdflatex)
├── paper.pdf           # Compiled PDF (output)
├── README.md                 # This file
├── CHANGELOG.md              # Version history
├── LICENSE                   # MIT (code)
├── LICENSE-data              # CC-BY-4.0 (data, figures, paper)
├── .zenodo.json              # Zenodo metadata
├── scripts/
│   ├── compute_twins_dr9.py             # Twin primes by DR mod 9 (patched)
│   ├── compute_permutation_families.py  # Pairs (p, p+18k) by DR
│   ├── compute_R_regression.py          # R_emp vs R_HL
│   ├── verify_global_rho.py             # NEW: OEIS aggregate benchmark
│   ├── analysis.py                      # NEW: slope fit + chi-square
│   ├── make_plots.py                    # Figures 1, 3, 4
│   └── make_plots_aggregate.py          # NEW: figure 2 (decay)
├── data/
│   ├── twins_dr9.csv             # 18 checkpoints, X from 1e4 to 1e9
│   ├── permutation_families.csv  # 9 values of k, DR breakdown
│   ├── R_regression.csv          # ratio R_emp / R_HL
│   └── global_rho.csv            # NEW: aggregate rho_glob from OEIS
└── figures/
    ├── convergenza_gemelli.pdf      # Twin DR convergence
    ├── decadimento_aggregato.pdf    # NEW: aggregate decay (log-log)
    ├── equipartizione_dr.pdf        # DR equidistribution for g=18k
    └── conferma_HL.pdf              # HL ratio rho(k) at X=5e7
```

## Dependencies

- Python ≥ 3.10
- numpy ≥ 1.20
- scipy ≥ 1.7
- matplotlib ≥ 3.5
- pdflatex (texlive-latex-base + texlive-fonts-recommended)

Minimal installation on Ubuntu/Debian:
```bash
pip install numpy scipy matplotlib
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

## Reproduction pipeline

The five computation scripts are independent. Run them in the
following order to reproduce all the data and figures in the paper.

### 1. Twin primes classified by DR mod 9

```bash
cd scripts
python3 compute_twins_dr9.py --xmax 1_000_000_000
```

- Runtime: ~30 s on a single CPU core
- Memory: ~2 GB peak (sieve + temporary bool mask)
- Output: `data/twins_dr9.csv` with 18 logarithmic checkpoints
- Populates Table 1 and Figure 1 of the paper

For a quick check (~12 s, ~100 MB) restrict to $X = 10^8$:
```bash
python3 compute_twins_dr9.py --xmax 100_000_000
```

### 2. Aggregate Hardy–Littlewood benchmark up to $X = 10^{18}$

```bash
python3 verify_global_rho.py --out ../data/global_rho.csv
```

- Runtime: <1 s (no sieve; uses OEIS A007508 + scipy.integrate)
- Output: `data/global_rho.csv` with $\rho_{\rm glob}(10^n)$ for $n=7,\ldots,18$
- Populates Table 2 and Figure 2 of the paper

### 3. Gap family $g = 18k$, DR equidistribution

```bash
python3 compute_permutation_families.py --xmax 10_000_000 \
    --k 5 7 10 11 12 13 15 17 20
```

- Runtime: ~5 s
- Memory: ~20 MB
- Output: `data/permutation_families.csv`
- Populates Table 3 and Figure 3 of the paper

### 4. Quantitative HL confirmation for $g = 18k$

```bash
python3 compute_R_regression.py --xmax 50_000_000 \
    --k 5 7 10 11 12 13 15 17 20
```

- Runtime: ~30 s
- Memory: ~60 MB
- Output: `data/R_regression.csv`
- Populates Table 4 and Figure 4 of the paper

### 5. Statistical analysis (slope fit + chi-square)

```bash
python3 analysis.py
```

- Runtime: <1 s
- Reads `data/global_rho.csv` and `data/twins_dr9.csv`
- Prints:
  - log–log slope of $|\rho_{\rm glob}-1|$ vs $X$ on $n=8,\ldots,18$
  - chi-square test of equidistribution at every checkpoint
- Numbers reported in §3.1–§3.2 of the paper

### 6. Figure generation

```bash
python3 make_plots.py            # convergenza, equipartizione, conferma_HL
python3 make_plots_aggregate.py  # decadimento_aggregato (NEW)
```

- Generates `.pdf` and `.png` versions in `figures/`
- Uses Times New Roman serif (falls back to DejaVu Serif if missing)

### 7. LaTeX compilation

```bash
cd ..
pdflatex paper.tex
pdflatex paper.tex   # second pass to resolve cross-references
```

## Numerical verification of the HL constant

To check that `H(g)` is computed correctly (an earlier draft had a
bug in the factorization: factors of 2 were not stripped before trial
division, which led to spurious "anomalies"):

```bash
python3 -c "
from compute_R_regression import odd_prime_factors, hl_constant_for_gap
for k in [5,7,10,11,12,13,15,17,20]:
    g = 18*k
    print(f'k={k:2d} g={g:3d}  odd primes={odd_prime_factors(g)}  H={hl_constant_for_gap(g):.4f}')
"
```

Expected output:
```
k= 5 g= 90  odd primes=[3, 5]   H=3.5209
k= 7 g=126  odd primes=[3, 7]   H=3.1688
k=10 g=180  odd primes=[3, 5]   H=3.5209
k=11 g=198  odd primes=[3, 11]  H=2.9341
k=12 g=216  odd primes=[3]      H=2.6406
k=13 g=234  odd primes=[3, 13]  H=2.8807
k=15 g=270  odd primes=[3, 5]   H=3.5209
k=17 g=306  odd primes=[3, 17]  H=2.8167
k=20 g=360  odd primes=[3, 5]   H=3.5209
```

## Cross-check against OEIS A007508

`compute_twins_dr9.py` excludes the pair $(3,5)$ for modular
uniformity, so its `total` column equals OEIS A007508 minus 1 for
$X = 10^n$, $n \ge 2$. Verified through $n = 9$:

| n | A007508 (incl. (3,5)) | script output | Δ |
|---:|---:|---:|---:|
| 7 | 58,980 | 58,979 | 1 ✓ |
| 8 | 440,312 | 440,311 | 1 ✓ |
| 9 | 3,424,506 | 3,424,505 | 1 ✓ |

## Methodological notes

1. **Sanity check**: all scripts internally verify that no pair falls
   in inadmissible DR classes (`assert` in
   `compute_permutation_families.py`).
2. **Replacing $x/(\log x)^2$ with $\operatorname{li}_2(x)$**: the
   regression and aggregate scripts use the secondary logarithmic
   integral to reduce the truncation error of the raw asymptotic.
   The difference is significant (~5% at $X = 10^7$).
3. **$\operatorname{li}_2$ lower bound**: `verify_global_rho.py` uses
   $x_0 = 2$ (configurable via `--x0`); `compute_R_regression.py` uses
   $x_0 = 5$. At the scales tabulated in the paper the choice is
   immaterial ($< 4$ for $X \ge 10^8$); the v1.2.1 paper makes the
   choice explicit in §3.2.
4. **Constant $C_2$**: used to 25 decimal places (Wrench, 1961;
   OEIS A005597) to avoid round-off in $H(g) = 2 C_2 \prod (p-1)/(p-2)$.
5. **Scaling limits**: the scripts are designed for $X \le 10^9$
   (bool in-memory sieve, ~1 GB at $10^9$). For larger $X$ a
   disk-segmented sieve is needed (e.g. `primesieve`); for $g = 2$,
   $X$ up to $10^{18}$ is accessible *aggregatedly* via OEIS A007508
   (no per-class breakdown is publicly tabulated).
6. **Statistical interpretation**: `analysis.py` reports both a
   linear regression on $\log X$ and a chi-square goodness-of-fit at
   each checkpoint. At $X = 10^9$ the chi-square is 0.82 (df $= 2$,
   $p = 0.67$), well above any conventional significance threshold.

## Citation

If you use this material, please cite the Zenodo Concept DOI:

> Ferraiuolo, G. (2026). *Modulo 9 classification of prime pairs
> $(p, p+g)$: an empirical verification of the Hardy–Littlewood
> conjecture for $g = 18k$.* Zenodo.
> <https://doi.org/10.5281/zenodo.20117270>

BibTeX entry:

```bibtex
@software{ferraiuolo_2026_primes_mod9,
  author       = {Ferraiuolo, Giovanni},
  title        = {Modulo 9 classification of prime pairs (p, p+g):
                  an empirical verification of the Hardy--Littlewood
                  conjecture for g = 18k},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20117270},
  url          = {https://doi.org/10.5281/zenodo.20117270}
}
```

## Licensing

This repository uses a **dual-license** scheme that is standard for
research compendia:

| Content | License | File |
|---------|---------|------|
| Python source code (`scripts/*.py`) | MIT License | `LICENSE` |
| Datasets (`data/*.csv`), figures (`figures/*`), paper (`paper.{tex,pdf}`) | Creative Commons Attribution 4.0 International (CC-BY-4.0) | `LICENSE-data` |

In short:

- You may use, modify, and redistribute the **code** under MIT terms
  (attribution appreciated, not legally required for derivative code).
- You may share and adapt the **data, figures, and paper text** under
  CC-BY-4.0, provided you give appropriate credit (citation, link to
  license, indication of changes).

## Acknowledgments

The author thanks Will Sawin (Columbia University) for valuable
comments on MathOverflow ([Q511187](https://mathoverflow.net/q/511187)),
in particular for the pointer to Sawin–Shusterman's *Annals* 2022
paper as the function-field analogue, and for suggesting the
extension of the empirical benchmark to the twin-prime case at
$X = 10^{18}$.

This work made use of AI assistance (Anthropic Claude) for code
review, LaTeX formatting, and bibliographic search; all mathematical
content, experimental design, and intellectual decisions are those of
the author.

## Contact

Giovanni Ferraiuolo — g.ferraiuolo2@gmail.com
