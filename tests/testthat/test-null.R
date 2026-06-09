test_that("p_emp in [0,1] and every permutation has overlap >= 1", {
  set.seed(42)
  N <- 200L
  m <- 10L
  null <- ohg_permutation_null(N = N, m = m, B = 500L, boundaries = seq_len(N))
  expect_length(null$log_stat_b, 500L)
  expect_true(all(null$overlap_b >= 1L)) # R_b subset U => never zero overlap

  log_stat_obs <- min(stats::phyper(
    0:(m - 1), m, N - m, 1:m,
    lower.tail = FALSE, log.p = TRUE
  ))
  p_emp <- (1 + sum(null$log_stat_b <= log_stat_obs)) / (1 + length(null$log_stat_b))
  expect_gte(p_emp, 0)
  expect_lte(p_emp, 1)
})

test_that("CALIBRATION: random ranking => p_emp ~ Uniform (the correction works)", {
  set.seed(7)
  N <- 300L
  m <- 15L
  B <- 999L
  null <- ohg_permutation_null(N = N, m = m, B = B, boundaries = seq_len(N))

  # Under H0 a random placement of m hits is exchangeable with the null draws,
  # so the empirical p-values must be ~ Uniform(0, 1).
  reps <- 400L
  p_vals <- vapply(seq_len(reps), function(i) {
    pos <- sort(sample.int(N, m))
    log_pv <- stats::phyper(seq_len(m) - 1L, m, N - m, pos,
      lower.tail = FALSE, log.p = TRUE
    )
    log_obs <- min(log_pv)
    (1 + sum(null$log_stat_b <= log_obs)) / (1 + B)
  }, numeric(1))

  ks <- suppressWarnings(stats::ks.test(p_vals, "punif"))
  expect_gt(ks$p.value, 0.01) # do NOT reject uniformity at alpha = 0.01
  expect_lt(abs(mean(p_vals) - 0.5), 0.06)
})
