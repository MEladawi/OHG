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

test_that("split null-summary + per-pathway helpers reproduce compute_effect", {
  # The hot path hoists the null-side summary (computed once per set size) out of
  # the per-pathway loop. The split must be numerically identical to the
  # single-shot compute_effect for every gating outcome.
  N <- 100L
  m <- 8L
  B <- 1500L
  le_idx_b <- make_null_le(N, m, B)
  le_obs <- 1:6

  # (a) ordinary case: finite NLES
  w <- runif(N, 0.1, 5)
  ce <- compute_effect(le_obs, w, le_idx_b, -1, TRUE, 1000L, 10L)
  summ <- .nles_null_summary(le_idx_b, w, TRUE, 1000L, 10L)
  sp <- .nles_from_summary(le_obs, w, -1, summ)
  expect_equal(sp$E_obs, ce$E_obs)
  expect_equal(sp$NLES, ce$NLES)
  expect_equal(sp$NLES_signed, ce$NLES_signed)

  # (b) gated case: degenerate spread => NA NLES but real E_obs, same as before
  wc <- rep(2, N)
  ce2 <- suppressWarnings(compute_effect(le_obs, wc, le_idx_b, 1, TRUE, 1000L, 10L))
  summ2 <- suppressWarnings(.nles_null_summary(le_idx_b, wc, TRUE, 1000L, 10L))
  sp2 <- .nles_from_summary(le_obs, wc, 1, summ2)
  expect_equal(sp2$E_obs, ce2$E_obs)
  expect_true(is.na(sp2$NLES))

  # (c) the null summary warns once (not per pathway) when gated
  expect_warning(.nles_null_summary(le_idx_b, wc, TRUE, 1000L, 10L), "near-zero spread")
})

test_that(".eb_from_leidx is the single E_b definition and abs's signed weights", {
  set.seed(7)
  null <- ohg_permutation_null(N = 120L, m = 8L, B = 200L)
  w_signed <- rnorm(120)
  # E_b computed via the helper ...
  eb_helper <- .eb_from_leidx(null$le_idx_b, w_signed, robust = TRUE)
  # ... must equal the explicit abs+median formula it replaces.
  center <- stats::median
  eb_manual <- vapply(null$le_idx_b, function(idx) center(abs(w_signed)[idx]), numeric(1))
  expect_identical(eb_helper, eb_manual)
})
