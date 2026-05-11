# Paper: Modulo 9 classification of prime pairs $(p, p+18k)$

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20117270.svg)](https://doi.org/10.5281/zenodo.20117270)
[![License: MIT](https://img.shields.io/badge/Code-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20%26%20Paper-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**GitHub**: <https://github.com/gferraiuolo2/prime-pairs-mod9>  
**Zenodo (DOI)**: <https://doi.org/10.5281/zenodo.20117271>  
**Contact**: g.ferraiuolo2@gmail.com

---

LaTeX source, Python computation scripts and datasets accompanying the
paper "Modulo 9 classification of prime pairs $(p,p+g)$: an empirical
verification of the Hardy–Littlewood conjecture for $g=18k$"
by Giovanni Ferraiuolo.

## Repository layout

```
paper_v6/
├── paper_v6_en.tex       # LaTeX source (compile with pdflatex)
├── paper_v6_en.pdf       # Compiled PDF (output)
├── README.md             # This file
├── scripts/
│   ├── compute_twins_dr9.py             # Counts twin primes by DR mod 9
│   ├── compute_permutation_families.py  # Counts (p, p+18k) pairs by DR
│   ├── compute_R_regression.py          # Compares R_emp vs R_HL
│   └── make_plots.py                    # Generates the PDF/PNG figures
├── data/
│   ├── twins_dr9.csv             # Output of compute_twins_dr9.py
│   ├── permutation_families.csv  # Output of compute_permutation_families.py
│   └── R_regression.csv          # Output of compute_R_regression.py
└── figures/
    ├── convergenza_gemelli.pdf
    ├── equipartizione_dr.pdf
    └── conferma_HL.pdf
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

The three computation scripts are independent. Run them in the
following order to reproduce all the data and figures in the paper.

### 1. Twin primes classified by DR mod 9

```bash
cd scripts
python3 compute_twins_dr9.py --xmax 100_000_000
```

- Runtime: ~12 s on a single CPU core
- Memory: ~100 MB (bool sieve over 10⁸)
- Output: `data/twins_dr9.csv` with counts at logarithmic checkpoints
- Populates Table 1 and Figure 1 of the paper

For X = 10⁹ (a stronger test, if you have ≥1 GB RAM):
```bash
python3 compute_twins_dr9.py --xmax 1_000_000_000
```

### 2. Gap family g = 18k, DR equidistribution

```bash
python3 compute_permutation_families.py --xmax 10_000_000 \
    --k 5 7 10 11 12 13 15 17 20
```

- Runtime: ~5 s
- Memory: ~20 MB
- Output: `data/permutation_families.csv`
- Populates Tables 2 and 3 and Figure 2 of the paper

To extend to k = 23, 25, 30: add values to the `--k` flag. The script
raises an exception if the sieve is insufficient for the largest k
(it requires a sieve up to `xmax + 18*max(k)`).

### 3. Quantitative HL confirmation

```bash
python3 compute_R_regression.py --xmax 50_000_000 \
    --k 5 7 10 11 12 13 15 17 20
```

- Runtime: ~30 s
- Memory: ~60 MB
- Output: `data/R_regression.csv`
- Populates Table 4 and Figure 3 of the paper

For X = 10⁹, use `--xmax 1_000_000_000` (requires ~1.2 GB RAM).

### 4. Figure generation

```bash
python3 make_plots.py
```

- Generates all files in `figures/` (PDF + high-resolution PNG)
- Uses Times New Roman serif (falls back to DejaVu Serif if missing)

### 5. LaTeX compilation

```bash
cd ..
pdflatex paper_v6_en.tex
pdflatex paper_v6_en.tex   # second pass to resolve cross-references
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

The values of H are independent of k within the same sub-family that
shares the same odd-prime set (e.g. {3,5} for k=5,10,15,20), as
predicted by the formula H(g) = 2*C2 * prod_{p | g, p>2} (p-1)/(p-2).

## Methodological notes

1. **Sanity check**: all scripts internally verify that no pair falls
   in inadmissible DR classes (assertion in
   `compute_permutation_families.py`).
2. **Replacing x/(log x)² with li₂(x)**: the regression script uses
   the secondary logarithmic integral to reduce the truncation error
   of the raw asymptotic. The difference is significant (~5% at
   x = 10⁷, ~3% at x = 10⁹).
3. **Constant C₂**: used to 70 decimal places (Wrench, 1961;
   OEIS A005597) to avoid round-off errors.
4. **Scaling limits**: the scripts are designed for X ≤ 10⁹ (bool
   in-memory sieve). For larger X, a disk-segmented sieve is needed
   (e.g. primesieve).

## Citation

If you use this material, please cite the Zenodo archive:

> Ferraiuolo, G. (2026). *Modulo 9 classification of prime pairs
> (p, p+g): an empirical verification of the Hardy–Littlewood
> conjecture for g = 18k.* Zenodo.
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
| Datasets (`data/*.csv`), figures (`figures/*`), paper (`paper_v6_en.{tex,pdf}`) | Creative Commons Attribution 4.0 International (CC-BY-4.0) | `LICENSE-data` |

In short:

- You may use, modify, and redistribute the **code** under MIT terms
  (attribution appreciated, not legally required for derivative code).
- You may share and adapt the **data, figures, and paper text** under
  CC-BY-4.0, provided you give appropriate credit (citation, link to
  license, indication of changes).

Suggested citation for any reuse:

> Ferraiuolo, G. (2026). *Modulo 9 classification of prime pairs
> (p, p+g): an empirical verification of the Hardy–Littlewood
> conjecture for g = 18k.* Zenodo.
> <https://doi.org/10.5281/zenodo.20117271>

## Contact

Giovanni Ferraiuolo — g.ferraiuolo2@gmail.com
