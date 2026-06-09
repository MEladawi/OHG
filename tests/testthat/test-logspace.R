test_that("log path stays finite where exp(log_stat) underflows to 0", {
  N <- 2000L
  ranked <- paste0("g", seq_len(N))
  T_eff <- ranked[1:50] # perfect top-50 enrichment: an extreme signal
  res <- ohg_statistic(ranked, T_eff, N)

  expect_true(is.finite(res$log_stat))
  expect_lt(res$log_stat, -100) # far beyond the double-precision exp range
  expect_equal(exp(res$log_stat), 0) # underflows, as expected

  nl10 <- -res$log_stat / log(10)
  expect_true(is.finite(nl10) && nl10 > 40)
})
