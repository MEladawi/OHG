# OHG — Ordered Hypergeometric enrichment

Pathway enrichment for a **ranked gene list** — for example, all genes ordered from
most up-regulated to most down-regulated by a differential expression analysis — that
doesn't make you guess where to draw the "significant genes" cutoff.

## Background & motivation

Classic over-representation analysis (ORA) — tools like **Enrichr**, DAVID, or g:Profiler —
starts by splitting your genes into a "significant" set and everything else, using a
threshold: usually an FDR or p-value cut (e.g. FDR < 0.05), often combined with a
fold-change filter. It then asks whether each pathway is over-represented in that set. Two
things make this lossy:

- **The threshold is arbitrary.** FDR < 0.05 vs < 0.10, with or without a fold-change
  filter, give different gene sets and different enrichment calls.
- **The ranking is thrown away.** Inside the significant set, genes are an unordered bag —
  the test is *pure overlap*, so a gene that barely cleared the cut counts exactly the same as
  your strongest hit, and everything below the cut counts for nothing.

**OHG keeps the ranking and removes the threshold.** It walks down your ranked list, tries
*every* cutoff, scores each with a hypergeometric test, and keeps the most surprising one. The
genes above that cutoff are the **leading edge** — the part of the pathway driving the signal.
Since it tried many cutoffs, OHG then **calibrates** the score against thousands of random
rankings to give an honest, comparable p-value (the *minimum-hypergeometric*, or mHG, test).

## Installation

```r
# install.packages("remotes")
remotes::install_github("MEladawi/OHG")
```

## A quick example

**The scenario:** you ran a differential expression analysis, so you have a results table
like this:

| gene   | log2FC | pvalue |
|--------|-------:|-------:|
| GENE1  |   4.2  | 1e-8   |
| GENE2  |  −3.8  | 3e-7   |
| GENE3  |   3.1  | 5e-6   |
| …      |   …    | …      |

**Step 1 — rank the genes**, most important first. A standard metric is *signed significance*,
which combines the p-value's strength with the fold-change's direction:

```r
library(OHG)

# Recommended: log2FC * -log10(p) -- the sign and size of the fold change set the
# direction and spread, -log10(p) weights it by evidence. Alternatives: signed
# significance sign(log2FC) * -log10(p), log2FC alone, or a moderated-t / Wald z.
# (In real use, clean the fold change first -- see "Cleaning the fold change" below.)
de$rank_stat <- de$log2FC * -log10(de$pvalue)
de <- de[order(de$rank_stat, decreasing = TRUE), ]
```

**Step 2 — run OHG** on that list, your pathways, and a per-gene effect size. A pathway is
just a named set of genes (or a `.gmt` file):

```r
pathways <- list(
  CELL_CYCLE = c("GENE1", "GENE3", "GENE12", "GENE41"),
  APOPTOSIS  = c("GENE2", "GENE9",  "GENE57")
)

res <- ohg_enrichment(
  ranked    = de$gene,
  gene_sets = pathways,
  rank_stat = de$rank_stat,     # log2FC * -log10(p) recommended (alt: log2FC alone, moderated-t, Wald z)
  weight    = abs(de$log2FC),   # effect size -> NLES; |LFC| recommended (alt: abs(log2FC)*-log10(p), -log10(p))
  direction = "up",             # test the up-regulated end; omit to test both ends
  seed      = 1
)
```

**What you get back** — one row per pathway, most enriched first:

| pathway    | set_size | overlap | leading_edge_fraction | p_value | p_adjust | NLES |
|------------|---------:|--------:|----------------------:|--------:|---------:|-----:|
| CELL_CYCLE |       45 |      18 |                  0.40 |  0.0005 |   0.0015 |  2.8 |
| APOPTOSIS  |       30 |       4 |                  0.13 |   0.21  |    0.31  |  0.4 |

That's the whole idea: **a ranked gene list in, ranked pathways out.** `CELL_CYCLE` is enriched
near the top of your list (small `p_adjust`) *and* its genes moved strongly (high `NLES`);
`APOPTOSIS` is neither. The sections below explain every column and input.

## Reading the output

OHG answers three different questions and keeps them in **separate columns** (it never blends
them into one score):

| Axis | The question | Key columns |
|---|---|---|
| **Significance** — is it real? | Is the top-of-list overlap more than chance? | `p_value`, `p_adjust` |
| **Localization** — where is it? | Where is the cutoff, and how much of the pathway is above it? | `cutoff_rank`, `leading_edge_fraction` |
| **Magnitude** — how strong? | How hard did the leading-edge genes move? | `NLES`, `NLES_signed` |

A pathway can be significant on just a few of its genes, or modest across most of them —
keeping the axes apart lets you tell those apart. (Multiple-testing correction is
Benjamini–Hochberg/FDR by default.)

## The two inputs that shape the result: `rank_stat` and `weight`

The quick example passed two optional inputs. They use **different** parts of your DE table and
do **different** jobs — confusing them is the one easy mistake:

- **`rank_stat` decides the order.** The values you ranked by. OHG uses it for ties and, via
  its **sign**, to tell up- from down-regulated genes — it is never plugged into the test
  arithmetic. The recommended choice is `log2FC * -log10(p)`: the fold change sets the direction
  and spread, `-log10(p)` weights it by evidence, and the product is continuous (few ties).
  Alternatives: signed significance `sign(log2FC) * -log10(p)`, `log2FC` alone, or a moderated-t /
  Wald z statistic. Omit it and OHG simply trusts the order you gave.
- **`weight` measures the push.** Each gene's effect size, used only for the magnitude score
  (`NLES`); OHG uses its size, never its sign. The recommended choice is `abs(log2FC)` on a
  **cleaned** fold change — shrunken with `ohg_shrink_lfc()` if you have a standard error, else
  winsorized with `ohg_winsorize()` — so noisy low-count genes don't dominate. Alternatives:
  `abs(log2FC) * -log10(p)` (the "π-value") or `-log10(p)`. Omit it and `NLES` is `NA` (the overlap
  and p-value columns still come back).

Mnemonic: **`rank_stat` decides where a gene sits; `weight` measures how hard it moved.**

Keep the two complementary: significance already lives in the `p_value` column, so a pure
effect size like `abs(log2FC)` keeps the magnitude axis separate. Folding significance into the
weight (`abs(log2FC) * -log10(p)`) is allowed but overlaps `p_value`; for noisy
large-fold-change genes, clean the `|LFC|` first — `ohg_shrink_lfc()` (below) using the standard
error, or `ohg_winsorize()` when you only have LFC and p.

> When a weight uses `-log10(p)`, build it from the p-values your DE tool reports — never from a
> p that has underflowed to 0 (`-log10(0) = Inf`). OHG errors on non-finite weights.

Pathways can also be a standard **`.gmt` file** (MSigDB / Enrichr / g:Profiler) or a
`GSEABase::GeneSetCollection`. Nothing else changes — only the `gene_sets` argument:

```r
res <- ohg_enrichment(
  ranked    = de$gene,
  gene_sets = "path/to/hallmark.gmt",   # <- the only difference vs. the quick example
  rank_stat = de$rank_stat,
  weight    = abs(de$log2FC),
  direction = "up",
  seed      = 1
)
```

## Up, down, or both

OHG scans only the **top** of the list it's given, so `direction` points it at the end you
care about:

- `"up"` — your list as-is; enrichment among the **most up-regulated** genes.
- `"down"` — the reversed list; enrichment among the **most down-regulated** genes.
- `"both"` — both ends, reported separately and corrected together (`NLES_signed > 0` = up,
  `< 0` = down).

If you don't set it, OHG infers it from `rank_stat`: a metric with both signs (like signed
significance) defaults to `"both"`; an all-positive or absent one defaults to `"up"`.

## Ranking and interpreting your results

To get a shortlist, keep the significant pathways, then sort by magnitude:

```r
library(dplyr)
res |>
  filter(p_adjust < 0.05) |>   # 1. keep what's significant
  arrange(desc(NLES_signed))   # 2. rank those by effect strength
```

Use `leading_edge_fraction` to **interpret**, not to rank — it tells you *coverage*:

- **few genes, large effect** (low fraction, high `NLES`) — a strong driver sub-module of a big pathway;
- **most genes, modest effect** (high fraction, low `NLES`) — broad, gentle involvement.

It's what stops "pathway X is enriched" from being over-read as the whole pathway when the
signal is really 4 of its 80 genes.

## Notes for specific situations

**Ties in `rank_stat`.** Genes with the *exact same* `rank_stat` value form a "tied block", and
their order inside it is arbitrary — just however your sort happened to break the tie. OHG never
puts the leading-edge cutoff *inside* a tied block; it takes the whole block or none of it. So if
those tied genes were shuffled, you'd get the identical result — the cutoff, p-value, and
leading edge don't move.

The trade-off is resolution. A metric with large tied blocks — e.g. ranking by `log2FC` when
many genes round to the same value, or by p-values that are all exactly 1 — gives OHG only a few
places to draw the cutoff, so the leading edge is coarse. A finer, more continuous metric (the
recommended `log2FC * -log10(p)`, or a moderated-t statistic) breaks those ties and lets OHG place
the cutoff more precisely.

**When `NLES` is `NA`.** The magnitude score needs enough spread in the random background. For
small or degenerate cases OHG sets `NLES = NA` and warns, but **never drops the pathway** — its
`p_value`/`p_adjust` stay valid ("enriched, but too small to size the effect reliably").

**Cleaning the fold change (recommended).** Both `log2FC * -log10(p)` and `abs(log2FC)` inherit the
low-count problem: a low-expression gene can show a wild, unreliable fold change that pulls it high
in the ranking *and* hands it a large weight. Clean the LFC **once**, upstream, then derive both
inputs from the same cleaned value so the two axes stay consistent:

```r
# Prefer shrinkage when you have a standard error; otherwise winsorize the tails.
clean_lfc <- if ("se" %in% names(de)) {
  ohg_shrink_lfc(de$log2FC, de$se)   # adaptive shrinkage (needs the 'ashr' package)
} else {
  ohg_winsorize(de$log2FC)           # cap the tails, no SE needed (default p = 0.99)
}

res <- ohg_enrichment(
  de$gene, pathways,
  rank_stat = clean_lfc * -log10(de$pvalue),  # signed ordering metric
  weight    = abs(clean_lfc)                  # magnitude for NLES
)
```

Use one or the other, not both. The OHG core never transforms your weight — it only *warns* about
pathological ones — so the cleaning stays visible in your own script.

**Big jobs (HPC).** Pass `n_cores > 1` and OHG runs the permutation work in parallel, setting up
(and restoring) the `future` backend itself — you only need the optional `future` and `furrr`
packages installed:

```r
res <- ohg_enrichment(de$gene, pathways, n_perm = 5000L, seed = 1, n_cores = 8L)
```

Results are reproducible across any number of cores for a fixed `seed`.

## Reference

### Inputs

| Input | Shape | Notes |
|---|---|---|
| `ranked_genes` | character vector, most important first | Must be unique (duplicates dropped with a warning). This list **is** the background — no separate background argument. |
| `gene_sets` | named list, `.gmt` path, or `GeneSetCollection` | Each pathway is intersected with `ranked_genes`, so `set_size` counts only genes present in your list. |
| `rank_stat` | numeric, same length & order as `ranked_genes` | Must be sorted non-increasing (OHG checks). Optional. |
| `weight` | numeric, same length & order as `ranked_genes`, finite | Per-gene effect magnitude for `NLES`; `abs()` is used. Optional. |

### Arguments of `ohg_enrichment()`

| Argument | Default | What it does |
|---|---|---|
| `ranked_genes`, `gene_sets` | required | Your ranked gene list and the pathways to test. |
| `rank_stat` | `NULL` | Values you ordered by; sets ties and the default direction. `NULL` trusts the given order. |
| `weight` | `NULL` | Effect magnitude feeding `NLES`. `NULL` ⇒ `NLES` is `NA` (overlap results still reported). |
| `direction` | `NULL` | `"up"`, `"down"`, or `"both"`; `NULL` infers it from the `rank_stat` sign. |
| `p_adjust_method` | `"BH"` | Correction; any of `stats::p.adjust.methods` (`"BH"`, `"BY"`, `"bonferroni"`, `"none"`, …). |
| `n_perm` | `2000` | Calibration permutations. Smallest possible p-value is `1/(n_perm + 1)`; raise for very significant hits. |
| `min_set_size` | `3` | Drop pathways with fewer than this many genes in your list. Inclusion only — never gates `NLES`. |
| `min_perm_nles` | `1000` | Minimum permutations before `NLES` is reported (else `NA`). |
| `min_nles_support` | `10` | Minimum distinct background magnitudes for `NLES` (else `NA`). |
| `robust_nles` | `TRUE` | `NLES` uses median/MAD (`TRUE`) or mean/SD (`FALSE`). |
| `collapse_both` | `FALSE` | With `"both"`, keep only the stronger direction per pathway (×2 penalty). |
| `method` | `"permutation"` | Calibration method; only `"permutation"` is implemented. |
| `seed` | `NULL` | Integer for reproducibility; the caller's RNG state is restored afterwards. |
| `n_cores` | `1` | Parallel workers (needs `furrr`/`future` + a `future::plan()`). Identical results for a fixed `seed`. |

Helpers: `read_gmt()`, `ohg_shrink_lfc()`, `ohg_winsorize()`, `plot_ohg_leading_edge()`.
See `?ohg_enrichment`.

### Output columns

One row per pathway (or pathway × direction), sorted by `p_value`:

| Column | Meaning |
|---|---|
| `pathway` | gene-set name |
| `direction` | `"up"`/`"down"` (only when not `"up"` alone) |
| `set_size` | pathway genes present in your list |
| `cutoff_rank` / `leading_edge_size` | rank of the chosen cutoff |
| `overlap` / `n_leading_edge` | pathway genes above the cutoff |
| `leading_edge_fraction` | `overlap / set_size` (coverage) |
| `neg_log10_mHG` / `mHG_stat` | the raw mHG statistic (not a p-value) |
| `p_value` / `p_adjust` | calibrated significance, raw / corrected |
| `p_adjust_method` | the correction used |
| `E_obs` | observed leading-edge magnitude (median `abs(weight)`) |
| `NLES` / `NLES_signed` | magnitude score, unsigned / direction-signed |
| `hits` | leading-edge gene names |

(`NES_OHG` is a deprecated alias for `NLES`.)

## License

MIT © 2026 Mahmoud Eladawi.
