# OHG

**Ordered Hypergeometric (minimum-hypergeometric, mHG) enrichment** for one
ranked gene list — e.g., genes ordered by a differential expression analysis
(by signed significance, fold change, or a test statistic) — against a collection
of gene sets (pathways).

OHG learns the optimal cutoff from the ranking itself and reports three
**orthogonal, never-pre-mixed** axes:

| Axis | Question | Columns |
|---|---|---|
| **Significance** | Is the top-of-list overlap more than chance? | `p_value`, `p_adjust` |
| **Localization** | Where is the leading edge, and how much of the set is in it? | `cutoff_rank`, `leading_edge_fraction` |
| **Magnitude** | How hard did the leading-edge genes move? | `NLES`, `NLES_signed` |

The test stays **count-based on purpose**: magnitude never enters the
hypergeometric p-value — that is what `NLES` (and, in a meta-tool, fGSEA) is for.
There is **no composite/blended ranking column**; combining axes is a one-line
downstream choice.

## Installation

```r
# install.packages("remotes")
remotes::install_github("MEladawi/OHG")    # or install_local() on a clone
```

Imports: `stats`, `dplyr`, `purrr`, `tibble`, `stringr`, `methods`.
Optional (Suggests): `furrr`/`future` (parallel), `ggplot2` (plot), `GSEABase`
(`GeneSetCollection` input), `ashr` (LFC shrinkage).

## Quick start

```r
library(OHG)

# A ranked list (most important first) and a few gene sets (pathways).
ranked <- paste0("g", 1:500)
sets <- list(
  PATHWAY_A = ranked[c(1:18, 40, 75)],   # strong top enrichment
  PATHWAY_B = sample(ranked, 30),
  PATHWAY_C = sample(ranked, 22)
)

res <- ohg_enrichment(ranked, sets, n_perm = 2000L, seed = 1)
res[, c("pathway", "set_size", "cutoff_rank", "leading_edge_fraction",
        "p_value", "p_adjust")]
```

Gene sets may also be supplied as a **`.gmt` file** (MSigDB / Enrichr /
g:Profiler) or a `GSEABase::GeneSetCollection` — all three funnel through
`coerce_gene_sets()`:

```r
res <- ohg_enrichment(ranked, "path/to/sets.gmt", n_perm = 2000L, seed = 1)
```

## `rank_stat` vs `weight` — two inputs, two jobs

These are **different** and easy to conflate:

- **`rank_stat`** decides *where genes sit*. It is **ordering-only** — never
  multiplied into the statistic (contrast GSEA). Its **sign** separates the
  up/down tails and sets the default `direction`. A finite signed statistic
  (`t`, Wald `z`, moderated-`t`, or signed significance `sign(LFC)·-log10(p)`)
  is ideal because it stays finite where a p-value underflows.
- **`weight`** measures *effect size* for `NLES`, on a **non-negative magnitude**
  (`abs(weight)` is used; the sign is ignored). Any finite choice works:
  `|LFC|`, `|LFC|·-log10(p)`, or `-log10(p)`.

> **Build any `-log10(p)` term in log space** from the finite log-p your pipeline
> already carries — never as `-log10(0) = Inf`. OHG hard-errors on non-finite
> `weight`.

```r
lfc   <- rnorm(500)
neglp <- -log10(runif(500))                 # built from finite values
res <- ohg_enrichment(
  ranked,
  sets,
  rank_stat = sort(sign(lfc) * neglp, decreasing = TRUE),  # ordering + sign
  weight    = abs(lfc) * neglp,                            # magnitude (|π-value|)
  seed = 1
)
```

## Directionality

The test scans only the **top** of the list it is given, so `direction` points it
at the right tail:

- `"up"` — list as-is (`NLES_signed = NLES`).
- `"down"` — runs on the reversed list (`dir_sign = -1`).
- `"both"` — both runs; each `(pathway, direction)` is a separate hypothesis and
  **all rows enter one pooled `p_adjust`**. `NLES_signed > 0` means enriched among
  up-regulated genes, `< 0` among down-regulated.

When `direction` is not supplied it is **inferred**: a signed `rank_stat`
(crosses zero) ⇒ `"both"`; a non-negative or absent `rank_stat` ⇒ `"up"`.
`collapse_both = TRUE` keeps the more significant direction per pathway with a
×2 Bonferroni penalty.

## Scoring vs interpreting a term

For a per-term **score** to rank significant pathways:

```r
library(dplyr)
res |>
  filter(p_adjust < 0.05) |>     # gate on significance
  arrange(desc(NLES_signed))     # rank on magnitude
```

`NLES_signed` is built to be comparable across terms and datasets, so it is
sufficient as the single ranking score. `leading_edge_fraction = overlap /
set_size` is **not** a score — it is an **interpretation** (coverage) column.
The two measure loosely-correlated axes:

- a **low-fraction / high-`NLES`** term is a strong *driver sub-module* of a big
  pathway (a few genes, large effect);
- a **high-fraction / low-`NLES`** term is coherent whole-pathway involvement at
  modest magnitude.

`leading_edge_fraction` is what stops "pathway X is enriched" being over-read as
whole-pathway when it is 4 of 80 genes.

## Ties

When `rank_stat` has ties, cutting *inside* a tie block would make the statistic
depend on arbitrary within-block order. OHG restricts cutoffs to **tie-block
boundaries** (a block is included whole or not at all), so the result is
invariant to within-block permutation. **Tip:** rank by a finer statistic to
shrink large tie blocks.

## Effect-size stability (`NLES = NA`)

`NLES` is gated to `NA` (with a specific warning) — but the pathway is **never
dropped** and its `p_value`/`p_adjust` stay valid — when the permutation null is
too thin to size the effect: `mad(E_b) ≈ 0`, fewer than `min_nles_support`
distinct null magnitudes, or `B < min_perm_nles`. This gate is **independent of
`min_set_size`**: a small significant set is reported as "enriched, but too small
to size the effect reliably".

## Optional LFC shrinkage

If your DE table carries a standard error, shrink before ranking — agnostic to DE
tool:

```r
shrunk <- ohg_shrink_lfc(lfc, se)     # needs 'ashr' (Suggests)
```

Shrinkage is a function of `(LFC, SE)`, not `(LFC, p)`; OHG never
reverse-engineers an SE from a p-value. If only LFC + p are available, winsorize
`|LFC|` (cap near the 99th percentile) rather than faking shrinkage.

## Parallelism (HPC)

OHG parallelizes the permutation nulls over **distinct set sizes** via `furrr`
when `n_cores > 1`. It **never calls `plan()` itself** — set the backend yourself:

```r
library(future)
plan(multicore)     # Unix/HPC: forked workers share the weight via copy-on-write
res <- ohg_enrichment(ranked, sets, n_perm = 5000L, seed = 1, n_cores = 8L)
```

`plan(multicore)` is recommended for Unix/HPC batch jobs (avoid in RStudio;
falls back to `multisession` on Windows). Results are reproducible across core
counts for a fixed `seed` (per-size RNG streams keyed to `m`). In a batch script,
read the worker count from the scheduler:

```r
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
```

## Output schema

One row per `(pathway[, direction])`, sorted by `p_value`: `pathway`,
`direction` (when `direction != "up"`), `set_size`, `cutoff_rank`,
`leading_edge_size`, `overlap`, `leading_edge_fraction`, `neg_log10_mHG`,
`mHG_stat`, `p_value`, `p_adjust`, `p_adjust_method`, `E_obs`, `NLES`,
`NLES_signed`, `n_leading_edge`, `hits`. (`NES_OHG` is a deprecated alias for
`NLES`.)

## License

MIT © 2026 Mahmoud Eladawi.
