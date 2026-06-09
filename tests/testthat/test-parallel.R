test_that("identical result for n_cores = 1 vs > 1 with the same seed", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  ranked <- paste0("g", 1:250)
  sets <- list(
    A = ranked[1:12],
    B = ranked[c(2:9, 100, 150)],
    C = ranked[c(1:6, 200:215)]
  )

  r1 <- ohg_enrichment(ranked, sets, n_perm = 600L, seed = 99, n_cores = 1L)
  r2 <- ohg_enrichment(ranked, sets, n_perm = 600L, seed = 99, n_cores = 2L)

  expect_equal(r1$pathway, r2$pathway)
  expect_equal(r1$p_value, r2$p_value)
  expect_equal(r1$neg_log10_mHG, r2$neg_log10_mHG)
})

test_that("build_nulls is reproducible and keyed by m, not order", {
  b1 <- build_nulls(c(5L, 8L), N = 100L, B = 200L, boundaries = seq_len(100L), seed = 7)
  b2 <- build_nulls(c(8L, 5L), N = 100L, B = 200L, boundaries = seq_len(100L), seed = 7)
  expect_equal(b1[["5"]]$log_stat_b, b2[["5"]]$log_stat_b)
  expect_equal(b1[["8"]]$log_stat_b, b2[["8"]]$log_stat_b)
})
