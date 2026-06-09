test_that("golden example (§2.4): N=10, T={g1,g2,g5}", {
  ranked <- paste0("g", 1:10)
  res <- ohg_statistic(ranked, T_eff = c("g1", "g2", "g5"), N = 10L)

  expect_equal(exp(res$log_stat), 1 / 15)
  expect_equal(res$log_stat, log(1 / 15))
  expect_identical(res$cutoff, 2L)
  expect_identical(res$overlap, 2L)
  expect_identical(ranked[res$le_idx], c("g1", "g2"))
})

test_that("tail values match the worked table", {
  tails <- stats::phyper(c(0, 1, 2), m = 3, n = 7, k = c(1, 2, 5), lower.tail = FALSE)
  expect_equal(tails, c(3 / 10, 1 / 15, 1 / 12))
})

test_that("no hits => log_stat 0 (log 1), overlap 0, cutoff NA", {
  res <- ohg_statistic(paste0("g", 1:5), T_eff = "zzz", N = 5L)
  expect_identical(res$log_stat, 0)
  expect_identical(res$overlap, 0L)
  expect_true(is.na(res$cutoff))
})

test_that("passing precomputed boundaries reduces to the distinct-rank path", {
  ranked <- paste0("g", 1:10)
  res <- ohg_statistic(ranked, c("g1", "g2", "g5"), 10L, boundaries = 1:10)
  expect_identical(res$cutoff, 2L)
  expect_identical(res$overlap, 2L)
})
