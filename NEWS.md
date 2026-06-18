# OHG 1.0.0

First stable release. The public API and the result schema are now stable.

## Bug fixes

* `NLES` is now reported under the default `max_cutoff_frac = 0.25`. The leading-edge
  magnitude null was being poisoned by permutation draws whose hits all fell below the
  cutoff cap `L`: those draws have no leading edge, but they were folded into the null
  as `NA`, which silently gated `NLES` to `NA` for every pathway — even with a clean,
  well-spread `weight` — and emitted a warning that wrongly blamed the weight vector.
  The null now conditions on a non-empty leading edge, the same event the observed
  score conditions on, so `NLES` is reported whenever the null is genuinely
  informative. The stability gate still returns `NA` only when there are truly too few
  usable permutations or the null spread is near zero.
* The caller's global random-number state is now restored on exit even when
  `seed = NULL` (a `seed = NULL` run could previously leave `.Random.seed` altered).

## Breaking changes

* Removed three result columns that were exact duplicates of others:
  `leading_edge_size` (identical to `cutoff_rank`), `n_leading_edge` (identical to
  `overlap`), and `NES_OHG` (identical to `NLES`). Use the canonical columns.

## Documentation

* `ohg_enrichment()` now documents every column of the returned tibble.

# OHG 0.1.1

Documentation and packaging only — no change to the analysis or results. Adds a
Zenodo DOI badge and `CITATION.cff`, clarifies the adaptive permutation budget,
sizes the comparison figures, and refines the OHG-vs-GSEA framing.

# OHG 0.1.0

First release.

Ordered / minimum-hypergeometric (mHG) pathway enrichment for ranked gene lists:
scan every cutoff, keep the most enriched, calibrate it against permutations, and
return an honest p-value, the exact leading-edge genes, and a separate effect size.

## Features

* `ohg_enrichment()` — one ranked gene list against a gene-set collection; one row
  per pathway, most enriched first.
* **Permutation-calibrated significance**, not the raw mHG statistic, so it is
  honest and BH/FDR-correctable.
* **Adaptive Besag–Clifford resampling** (default `method = "adaptive"`) removes the
  `1/(n_perm + 1)` p-value floor by spending extra permutations only where the
  significance decision is in doubt; `n_perm_used` / `resolution_limited` report what
  each pathway received. A fixed-`n_perm` `method = "permutation"` is also available.
* **Three separate axes, never pre-mixed:** significance (`p_value`, `p_adjust`),
  localization (`cutoff_rank`, `leading_edge_fraction`, `hits`), and magnitude
  (`NLES` / `NLES_signed`, a robust median/MAD Normalized Leading-Edge Score, gated
  to `NA` when the null is too thin to size reliably).
* **XL-mHG `L`/`X` restriction** (`max_cutoff_frac`, `min_hits`) keeps the optimal
  cutoff near the top, so a leading edge stays a short, interpretable gene list; the
  same restriction is applied to the null, so calibration stays valid.
* **Direction-aware:** `"up"`, `"down"`, or `"both"` (corrected together), with
  `collapse_both` for one row per pathway. Direction comes from the run and the sign
  of `rank_stat`, never from the effect weight.
* **Flexible gene-set inputs:** a named list, a `.gmt` file (MSigDB / Enrichr /
  g:Profiler), or a `GSEABase::GeneSetCollection`.
* **Helpers:** `read_gmt()`, `ohg_shrink_lfc()` (adaptive shrinkage via `ashr`),
  `ohg_winsorize()` (no-SE tail capping), and `plot_ohg_leading_edge()`.
* **Parallel and reproducible:** pass `n_cores > 1` (with `future` / `furrr`); OHG
  owns and restores the backend itself, and results are identical across core counts
  for a fixed `seed`.

## Recommended inputs

* Rank by signed significance, `rank_stat = sign(log2FC) * -log10(p)` — ordering by
  evidence, with magnitude left to `NLES`.
* `weight = abs(log2FC)`, cleaned first with `ohg_shrink_lfc()` (with a standard
  error) or `ohg_winsorize()`.

## Performance

A continuous metric has no ties, so the mHG kernel takes an O(*m*) fast path. A full
run on 19,484 genes against 3,759 GO:BP gene sets completes in about 1.5 minutes on 8
cores, with an adaptive permutation budget — a baseline of `n_perm = 2000`, escalating
to roughly 150,000 for near-significant pathways so their p-values resolve precisely.
