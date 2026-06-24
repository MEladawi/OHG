toy <- function() {
  set.seed(11)
  ranked <- paste0("g", 1:300)
  list(
    ranked = ranked,
    sets = list(
      HIT = ranked[c(1:12, 50, 80)], # strong top enrichment
      MEH = sample(ranked, 15),
      MEH2 = sample(ranked, 20)
    )
  )
}

test_that("default p_adjust is BH and recorded per row", {
  d <- toy()
  res <- ohg_enrichment_quiet(d$ranked, d$sets, n_perm = 500L, seed = 1)
  expect_true(all(res$p_adjust_method == "BH"))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_equal(res$p_adjust, stats::p.adjust(res$p_value, "BH"))
  expect_false(is.unsorted(res$p_value)) # default sort by p_value
})

test_that("changing method changes p_adjust; pooled once", {
  d <- toy()
  res_bonf <- ohg_enrichment_quiet(d$ranked, d$sets,
    p_adjust_method = "bonferroni",
    n_perm = 500L, seed = 1
  )
  expect_equal(res_bonf$p_adjust, stats::p.adjust(res_bonf$p_value, "bonferroni"))
  expect_true(all(res_bonf$p_adjust_method == "bonferroni"))
})

test_that("schema columns present; overlap>=1; no duplicate columns", {
  d <- toy()
  res <- ohg_enrichment_quiet(d$ranked, d$sets,
    weight = rev(seq_along(d$ranked)),
    n_perm = 1200L, seed = 1
  )
  needed <- c(
    "pathway", "set_size", "cutoff_rank", "overlap",
    "leading_edge_fraction", "neg_log10_mHG", "mHG_stat", "p_value",
    "p_adjust", "p_adjust_method", "E_obs", "NLES", "NLES_signed",
    "hits"
  )
  expect_true(all(needed %in% names(res)))
  expect_true(all(res$overlap >= 1L))
  expect_equal(res$leading_edge_fraction, res$overlap / res$set_size)
  # the removed duplicate columns must be gone
  expect_false(any(c("leading_edge_size", "n_leading_edge", "NES_OHG") %in% names(res)))
})

test_that("BH family counts tested pathways with no eligible leading edge", {
  set.seed(7)
  N <- 1500L
  ranked <- paste0("g", seq_len(N))
  # Small random sets with a tight cap: most have no hit in the top L = 5% and so
  # are tested (p = 1) but never emitted. Plant one genuine top hit.
  sets <- stats::setNames(
    lapply(seq_len(150), function(i) sample(ranked, 12)), paste0("S", seq_len(150))
  )
  sets[["TOP"]] <- ranked[1:12]

  res <- ohg_enrichment_quiet(
    ranked, sets,
    method = "permutation", n_perm = 800L, max_cutoff_frac = 0.05, seed = 1
  )
  pruned <- prune_gene_sets(coerce_gene_sets(sets), ranked, 3L, N)
  n_family <- length(pruned$kept)

  # More hypotheses tested than rows emitted (overlap == 0 sets are dropped).
  expect_gt(n_family, nrow(res))
  # Adjustment must use the full tested family, not the emitted subset.
  expect_equal(res$p_adjust, stats::p.adjust(res$p_value, "BH", n = n_family))
  expect_false(isTRUE(all.equal(
    res$p_adjust, stats::p.adjust(res$p_value, "BH")
  )))
})

test_that("adding tested-but-null pathways makes survivors more conservative", {
  N <- 1500L
  ranked <- paste0("g", seq_len(N))
  base_sets <- list(TOP = ranked[1:12], TOP2 = ranked[c(1:6, 40:45)])
  set.seed(3)
  # Filler drawn from the BOTTOM of the list: with L = 5% none have a top-L hit,
  # so they are tested (p = 1) and dropped, changing only the family size.
  filler <- stats::setNames(
    lapply(seq_len(200), function(i) sample(ranked[200:N], 10)), paste0("F", seq_len(200))
  )
  small <- ohg_enrichment_quiet(
    ranked, base_sets,
    method = "permutation", n_perm = 800L, max_cutoff_frac = 0.05, seed = 1
  )
  big <- ohg_enrichment_quiet(
    ranked, c(base_sets, filler),
    method = "permutation", n_perm = 800L, max_cutoff_frac = 0.05, seed = 1
  )
  # TOP's raw p-value is unchanged (null depends only on N, m, boundaries); only
  # the family size grew, so its adjusted p-value must rise.
  expect_equal(
    small$p_value[small$pathway == "TOP"], big$p_value[big$pathway == "TOP"]
  )
  expect_gt(
    big$p_adjust[big$pathway == "TOP"], small$p_adjust[small$pathway == "TOP"]
  )
})

test_that("collapse_both uses the pathway count as the BH family size", {
  set.seed(9)
  N <- 1000L
  ranked <- paste0("g", seq_len(N))
  rs <- c(seq(3, 0.01, length.out = N / 2), seq(-0.01, -3, length.out = N / 2))
  sets <- stats::setNames(
    lapply(seq_len(60), function(i) sample(ranked, 12)), paste0("S", seq_len(60))
  )
  sets[["UP"]] <- ranked[1:12]
  sets[["DOWN"]] <- ranked[(N - 11):N]

  res <- ohg_enrichment_quiet(
    ranked, sets,
    rank_stat = rs, direction = "both", collapse_both = TRUE,
    method = "permutation", n_perm = 600L, max_cutoff_frac = 0.1, seed = 1
  )
  pruned <- prune_gene_sets(coerce_gene_sets(sets), ranked, 3L, N)
  n_family <- length(pruned$kept) # collapse: one reported hypothesis per pathway

  expect_lte(nrow(res), n_family)
  expect_equal(res$p_adjust, stats::p.adjust(res$p_value, "BH", n = n_family))
})
