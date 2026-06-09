make_null_le <- function(N, m, B, seed = 1) {
  set.seed(seed)
  ohg_permutation_null(N, m, B, boundaries = seq_len(N))$le_idx_b
}

test_that("NLES uses abs(weight) and one code path for any magnitude", {
  N <- 100L
  m <- 8L
  B <- 1500L
  le_idx_b <- make_null_le(N, m, B)
  w <- runif(N, 0.1, 5)
  le_obs <- 1:6

  res_pos <- compute_effect(le_obs, w, le_idx_b, 1, TRUE, 1000L, 10L)
  res_neg <- compute_effect(le_obs, -w, le_idx_b, 1, TRUE, 1000L, 10L) # sign ignored

  expect_equal(res_pos$E_obs, stats::median(w[le_obs]))
  expect_equal(res_pos$NLES, res_neg$NLES)
  expect_true(is.finite(res_pos$NLES))
})

test_that("NLES_signed = dir_sign * NLES", {
  N <- 100L
  m <- 8L
  le_idx_b <- make_null_le(N, m, 1500L)
  w <- runif(N, 0.1, 5)
  up <- compute_effect(1:6, w, le_idx_b, 1, TRUE, 1000L, 10L)
  dn <- compute_effect(1:6, w, le_idx_b, -1, TRUE, 1000L, 10L)
  expect_equal(up$NLES_signed, up$NLES)
  expect_equal(dn$NLES_signed, -dn$NLES)
})

test_that("NLES = NA when weight missing", {
  res <- compute_effect(1:6, NULL, list(), 1, TRUE, 1000L, 10L)
  expect_true(is.na(res$NLES))
  expect_true(is.na(res$E_obs))
})

test_that("NLES = NA + warning when null spread ~ 0 (all-equal weight)", {
  N <- 100L
  m <- 8L
  le_idx_b <- make_null_le(N, m, 1500L)
  w <- rep(2, N) # constant => mad(E_b) = 0
  expect_warning(
    res <- compute_effect(1:6, w, le_idx_b, 1, TRUE, 1000L, 10L),
    "near-zero spread"
  )
  expect_true(is.na(res$NLES))
})

test_that("NLES = NA when B < min_perm_nles", {
  N <- 100L
  m <- 8L
  le_small <- make_null_le(N, m, 200L) # B too small
  w <- runif(N, 0.1, 5)
  expect_warning(r1 <- compute_effect(1:6, w, le_small, 1, TRUE, 1000L, 10L))
  expect_true(is.na(r1$NLES))
})

test_that("mean/SD path available behind robust = FALSE", {
  N <- 100L
  m <- 8L
  le_idx_b <- make_null_le(N, m, 1500L)
  w <- runif(N, 0.1, 5)
  le_obs <- 1:6
  res <- compute_effect(le_obs, w, le_idx_b, 1, FALSE, 1000L, 10L)
  expect_equal(res$E_obs, mean(w[le_obs]))
})
