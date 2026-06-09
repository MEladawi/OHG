# OHG — Ordered Hypergeometric enrichment

Pathway enrichment for a **ranked gene list** — for example, all genes ordered from
most up-regulated to most down-regulated by a differential expression analysis — that
doesn't make you guess where to draw the "significant genes" cutoff.

## Background & motivation

Classic over-representation analysis (ORA) — tools like **Enrichr**, DAVID, or g:Profiler —
starts by splitting your genes into a "significant" set and everything else, using a
threshold: usually an FDR or p-value cut (e.g. FDR < 0.05), often combined with a
fold-change filter. It then asks whether each pathway is over-represented in that
significant set. Two things make this lossy:

- **The threshold is arbitrary.** FDR < 0.05 vs < 0.10, with or without a fold-change
  filter, give different gene sets and different enrichment calls.
- **The ranking is thrown away.** Inside the significant set, genes are an unordered bag —
  the test is *pure overlap* (a hypergeometric / Fisher count), so a gene that barely cleared
  the cut counts exactly the same as your strongest hit, and everything below the cut counts
  for nothing.

**OHG uses the whole ranking and removes the threshold.** Instead of one fixed cut, it walks
down your ranked list and tries *every* cutoff, scoring each with a hypergeometric test
("how surprising is this much overlap in the top-k?"), and keeps the single most surprising
one. Because it works from positions rather than an in-or-out label, a gene's rank matters:
the genes above the chosen cutoff are the **leading edge** — the part of the pathway actually
driving the signal.

Because OHG peeked at many cutoffs and kept the best, that raw score is over-optimistic. So
it **calibrates**: it reruns the same procedure on thousands of random rankings and reports
where your real score falls among them. The result is an honest, comparable p-value — the
`minimum-hypergeometric (mHG)` test.

## What you get back: three separate answers

OHG reports three things that are easy to confuse, so it keeps them in **separate columns
and never blends them into one score**:

| Axis | The question it answers | Columns |
|---|---|---|
| **Significance** | Is the top-of-list overlap more than chance? | `p_value`, `p_adjust` |
| **Localization** | Where is the cutoff, and how much of the pathway sits above it? | `cutoff_rank`, `leading_edge_fraction` |
| **Magnitude** | How strongly did the leading-edge genes actually move? | `NLES`, `NLES_signed` |

These measure different things. A pathway can be highly significant on only a few of its
genes, or modestly significant across most of them — keeping the axes apart lets you tell
those apart. Combining them into a ranking is a one-liner *you* control (see
[Ranking your results](#ranking-your-results)); the package never does it for you.

## Installation

```r
# install.packages("remotes")
remotes::install_github("MEladawi/OHG")    # or install_local() on a clone
```

Required: `stats`, `dplyr`, `purrr`, `tibble`, `methods`.
Optional (only for extra features): `furrr`/`future` (parallel runs), `ggplot2` (plotting),
`GSEABase` (`GeneSetCollection` input), `ashr` (fold-change shrinkage).

## Quick start

```r
library(OHG)

# `ranked` is your gene list, most important first.
# `sets` maps each pathway name to its member genes.
ranked <- paste0("g", 1:500)
sets <- list(
  PATHWAY_A = ranked[c(1:18, 40, 75)],   # clustered near the top -> enriched
  PATHWAY_B = sample(ranked, 30),        # scattered -> not enriched
  PATHWAY_C = sample(ranked, 22)
)

res <- ohg_enrichment(ranked, sets, n_perm = 2000L, seed = 1)

# The headline columns:
res[, c("pathway", "set_size", "cutoff_rank", "leading_edge_fraction",
        "p_value", "p_adjust")]
```

Each row is one pathway, sorted by `p_value` (most enriched first). `set_size` is how many
of the pathway's genes are in your list, `cutoff_rank` is where OHG drew the line, and
`leading_edge_fraction` is the share of the pathway above it.

Your pathways can also come from a **`.gmt` file** (the standard MSigDB / Enrichr /
g:Profiler format) or a `GSEABase::GeneSetCollection` — just pass it instead of the list:

```r
res <- ohg_enrichment(ranked, "path/to/sets.gmt", n_perm = 2000L, seed = 1)
```

## Two optional inputs, two different jobs: `rank_stat` and `weight`

By default OHG just uses the order you handed it. Two optional arguments let you be more
precise, and they do **different** jobs — this is the one thing worth getting right:

- **`rank_stat` — *where each gene sits*.** The numeric values you ranked by (e.g. a signed
  test statistic). OHG uses it only for ordering and to handle ties; it is **never** plugged
  into the test arithmetic. Its **sign** is what separates up- from down-regulated genes.
- **`weight` — *how big each gene's effect is*.** Used only to compute the magnitude score
  (`NLES`). OHG uses its size (`abs(weight)`), never its sign.

A handy way to remember it: **`rank_stat` decides the order; `weight` measures the push.**

```r
lfc   <- rnorm(500)               # log fold-changes
neglp <- -log10(runif(500))       # significance, built from finite p-values

res <- ohg_enrichment(
  ranked, sets,
  rank_stat = sort(sign(lfc) * neglp, decreasing = TRUE),  # order + direction
  weight    = abs(lfc) * neglp,                            # effect magnitude
  seed = 1
)
```

> **Tip:** if a weight involves `-log10(p)`, build it from the small log-p your DE tool
> already provides — never as `-log10(0) = Inf`. OHG stops with an error on non-finite
> weights rather than silently mangling them.

## Up, down, or both

OHG only ever scans the **top** of the list it's given, so `direction` just points it at the
end you care about:

- `"up"` — your list as-is (enrichment among the most up-regulated genes).
- `"down"` — the reversed list (enrichment among the most down-regulated genes).
- `"both"` — runs each end separately and reports both, correcting them together. Here
  `NLES_signed > 0` means the pathway is enriched among up-regulated genes and `< 0` among
  down-regulated.

If you don't set `direction`, OHG infers it: a `rank_stat` that crosses zero (has both signs)
defaults to `"both"`; one that's all non-negative (or absent) defaults to `"up"`. Set
`collapse_both = TRUE` to keep only the stronger direction per pathway (with a ×2 penalty for
having looked at both).

## Ranking your results

To turn the table into a ranked shortlist, first keep the significant pathways, then sort by
magnitude:

```r
library(dplyr)
res |>
  filter(p_adjust < 0.05) |>   # 1. keep what's statistically real
  arrange(desc(NLES_signed))   # 2. rank those by effect strength
```

`NLES_signed` is designed to be comparable across pathways and datasets, so it's enough as a
single ranking score. Keep `leading_edge_fraction` around to **interpret**, not to rank — it
tells you *coverage*:

- **few genes, large effect** (low fraction, high `NLES`): a strong sub-module of a big pathway;
- **most genes, modest effect** (high fraction, low `NLES`): broad, gentle involvement.

It's what stops "pathway X is enriched" from being over-read as the whole pathway when it's
really 4 of 80 genes.

## Notes for specific situations

**Ties in `rank_stat`.** If many genes share the same value, the order within that block is
meaningless, so OHG only ever cuts *between* tied blocks — never inside one. Results are
therefore unaffected by how ties happen to be ordered. If you have large tied blocks, ranking
by a finer statistic will sharpen the result.

**When `NLES` comes back `NA`.** The magnitude score needs enough spread in the random
background to be trustworthy. For very small or degenerate cases OHG sets `NLES = NA` and warns
— but it **never drops the pathway**: its `p_value`/`p_adjust` are still valid. In short:
"enriched, but too small to size the effect reliably."

**Fold-change shrinkage (optional).** If your DE table has a standard error, you can shrink
noisy fold-changes before ranking:

```r
shrunk <- ohg_shrink_lfc(lfc, se)   # needs the 'ashr' package
```

This needs `(LFC, SE)`, not p-values — OHG never tries to reverse-engineer a standard error
from a p-value. If you only have LFC and p, cap extreme `|LFC|` values instead.

**Big jobs (HPC).** OHG can run the permutation step in parallel over distinct pathway sizes.
It never picks a backend for you — you set one with `future::plan()`:

```r
library(future)
plan(multicore)   # Unix/HPC; forked workers share memory efficiently
res <- ohg_enrichment(ranked, sets, n_perm = 5000L, seed = 1, n_cores = 8L)
```

Results are reproducible across any number of cores for a fixed `seed`. On a cluster, read the
core count from the scheduler: `as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))`.

## Full output columns

One row per pathway (or per pathway × direction), sorted by `p_value`:

`pathway`, `direction` (only when not `"up"`), `set_size`, `cutoff_rank`,
`leading_edge_size`, `overlap`, `leading_edge_fraction`, `neg_log10_mHG`, `mHG_stat`,
`p_value`, `p_adjust`, `p_adjust_method`, `E_obs`, `NLES`, `NLES_signed`, `n_leading_edge`,
`hits`. (`NES_OHG` is a deprecated alias for `NLES`.)

## License

MIT © 2026 Mahmoud Eladawi.
