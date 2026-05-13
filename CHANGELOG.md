# Changelog

All notable changes to this project are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [v1.2.1] — 2026-05-12

### Added
- **§3.2 of the paper**: new subsection *Aggregate verification
  against published tables*, comparing $\pi_2(X)$ from OEIS A007508
  (Oliveira e Silva) against $H(2)\,\operatorname{li}_2(X)$ for
  $X = 10^n$, $n = 7,\ldots,18$. Table `tab:rho_global` and figure
  `fig:rho_global`.
- **Figure 2** (`decadimento_aggregato.pdf`): log–log decay of
  $|\rho_{\rm glob}(X)-1|$ over 11 orders of magnitude, with
  guideline $\propto X^{-1/2}\log X$.
- **`verify_global_rho.py`**: standalone script that computes the
  aggregate ratio from OEIS A007508 (no sieve needed, <1 s runtime).
- **`analysis.py`**: standalone script with the two statistical
  analyses now cited in the paper:
  - least-squares slope of $\log(|\rho-1|/\log X)$ vs $\log X$,
    yielding $-0.45 \pm 0.03$ ($R^2 = 0.97$);
  - chi-square goodness-of-fit on the three DR classes for twin
    primes at every checkpoint (max $\chi^2 = 3.5$, min $p = 0.17$).
- **`make_plots_aggregate.py`**: renders the new Figure 2.
- **Sawin–Shusterman 2022** (Annals 196, 457–506) explicitly cited as
  the rigorous function-field analogue, following Sawin's
  recommendation on MathOverflow Q511187.
- **Oliveira e Silva** and **OEIS A007508** added as bibliographic
  sources for the aggregate $\pi_2(X)$ data.

### Changed
- **Twin-primes DR-mod-9 computation extended to $X = 10^9$** (Table 1
  now has three columns: $10^7, 10^8, 10^9$). All three admissible
  classes lie within $\pm 0.025$ pp of $1/3$ at the new scale.
- **Figure 1** (`convergenza_gemelli.pdf`) regenerated with the
  extended dataset: now 18 logarithmic checkpoints from $10^4$ to
  $10^9$ (was 11 up to $10^8$).
- **§3.1** rewritten to document the chi-square goodness-of-fit
  test ($\chi^2 = 0.82$, df $=2$, $p = 0.67$ at $X = 10^9$).
- **§6.3** updated to cross-reference §3.2 for the quantitative slope
  estimate rather than repeating the fit details.
- **Abstract** and **Conclusion** rephrased to reflect the new
  empirical ranges and the $X^{-1/2}\log X$ scaling argument.
- **Zenodo DOI** in the manuscript and README aligned to the
  Concept DOI `10.5281/zenodo.20117270` (was incorrectly pointing to
  the v1.2.0 Version DOI `…20117271`).

### Fixed
- **Memory bug in `compute_twins_dr9.py`**: the previous
  implementation materialised `np.arange(5, xmax + 1)` as an int64
  array, requiring 7.5 GiB of RAM at $X = 10^9$ (would crash on any
  machine with <8 GB free). Replaced with
  `np.flatnonzero(sieve[5:xmax+1] & sieve[7:xmax+3]) + 5`, which has
  identical output and peak memory ~2 GB.
- **Bibliographic typo in `Nicely11`**: removed the spurious
  reference to OEIS A001097 (which is the list of twin primes, not
  the counting function). The counting function is OEIS A007508,
  now cited separately.

### Verified
- Cross-check of the twin-prime totals against OEIS A007508 at
  $X = 10^7, 10^8, 10^9$: agreement modulo the exclusion of the
  pair $(3,5)$ (script output equals OEIS minus 1).

## [v1.2.0] — 2026-04-25

### Added
- §6.3 *On the magnitude of the empirical error* introduces the
  function-field context: Keating–Rudnick variance and
  Sawin–Shusterman power-saving twin primes over $\mathbb{F}_q[T]$.
- Cramér heuristic and Maier's *Primes in short intervals* cited as
  classical references on the error term.
- Acknowledgments section including disclosure of AI assistance and
  thanks to Will Sawin (MathOverflow Q511187).
- Empirical benchmark for twin primes at $X = 10^{18}$ documented in
  §6.3, with explicit values $\pi_2(10^{18}) = 808{,}675{,}888{,}577{,}436$
  and prediction $\approx 808{,}675{,}901{,}493{,}606$.

### Changed
- Bibliography expanded to include Maier (1985), Keating–Rudnick
  (IMRN 2014), Sawin–Shusterman (Annals 2022), and the MathOverflow
  question.

## [v1.1] — 2026-03-14

### Added
- Family $g = 18k$ extended to $k \in \{5,7,10,11,12,13,15,17,20\}$
  (was a smaller subset).
- §5 *Quantitative confirmation of Hardy–Littlewood* with
  Table 4 / Figure 3.

### Fixed
- `odd_prime_factors`: strip factors of 2 *before* trial division by
  odd primes. The previous version produced incorrect $H(g)$ values
  for some $g$, which had led to spurious "anomalies" in earlier
  drafts.

## [v1.0] — 2026-02-01

Initial release accompanying the first preprint.

[v1.2.1]: https://github.com/gferraiuolo2/prime-pairs-mod9/releases/tag/v1.2.1
[v1.2.0]: https://github.com/gferraiuolo2/prime-pairs-mod9/releases/tag/v1.2.0
[v1.1]:   https://github.com/gferraiuolo2/prime-pairs-mod9/releases/tag/v1.1
[v1.0]:   https://github.com/gferraiuolo2/prime-pairs-mod9/releases/tag/v1.0
