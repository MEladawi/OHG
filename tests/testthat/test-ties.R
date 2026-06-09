test_that("tie_boundaries gives every position for distinct ranks", {
  expect_identical(tie_boundaries(c(5, 4, 3, 2, 1)), 1:5)
  expect_identical(tie_boundaries(NULL, n = 4L), 1:4)
})

test_that("tie_boundaries collapses tie blocks to their last position", {
  # ranks: 5,5 | 4 | 3,3,3 | 1  -> boundaries at 2, 3, 6, 7
  expect_identical(tie_boundaries(c(5, 5, 4, 3, 3, 3, 1)), c(2L, 3L, 6L, 7L))
})

test_that("within-tie-block shuffles leave stat/cutoff/overlap unchanged", {
  set.seed(1)
  base_rank <- c(5, 5, 5, 4, 4, 3, 2, 2, 1, 1)
  ranked <- paste0("g", 1:10)
  T_eff <- c("g2", "g7", "g9")
  b0 <- tie_boundaries(base_rank)
  r0 <- ohg_statistic(ranked, T_eff, 10L, boundaries = b0)

  perm <- ranked
  for (block in list(1:3, 4:5, 6, 7:8, 9:10)) {
    perm[block] <- sample(ranked[block])
  }
  r1 <- ohg_statistic(perm, T_eff, 10L, boundaries = tie_boundaries(base_rank))

  expect_identical(r1$log_stat, r0$log_stat)
  expect_identical(r1$cutoff, r0$cutoff)
  expect_identical(r1$overlap, r0$overlap)
})

test_that("two hits in one tie block => one boundary, q jumps by 2, evals < m", {
  base_rank <- c(5, 5, 4, 3, 2, 1)
  ranked <- paste0("g", 1:6)
  T_eff <- c("g1", "g2", "g5") # g1, g2 share the first tie block
  b <- tie_boundaries(base_rank)
  is_hit <- ranked %in% T_eff
  q_all <- cumsum(is_hit)[b]
  keep <- q_all >= 1L & c(TRUE, diff(q_all) > 0L)

  expect_lt(sum(keep), length(T_eff)) # strictly fewer evaluations than m = 3
  expect_identical(q_all[keep][1L], 2L) # first kept boundary already carries 2 hits
})
