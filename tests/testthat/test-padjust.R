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
