# XL-mHG L (largest cutoff) / X (minimum prefix hits) restriction.

test_that("L makes the statistic ignore enrichment outside the top L", {
  ranked <- paste0("g", 1:100)
  bottom <- paste0("g", 51:60) # all hits in the bottom half

  full <- ohg_statistic(ranked, bottom, N = 100L) # L = N: finds the deep bunch
  restr <- ohg_statistic(ranked, bottom, N = 100L, L = 25L) # top 25 only

  expect_lt(full$log_stat, 0) # enriched somewhere deep
  expect_equal(restr$log_stat, 0) # nothing in the top 25 -> not enriched
  expect_equal(restr$overlap, 0L)
  expect_length(restr$le_idx, 0L)
})

test_that("cutoff never exceeds L for a top set", {
  ranked <- paste0("g", 1:100)
  top <- paste0("g", c(1, 2, 3, 30, 31, 32))

  restr <- ohg_statistic(ranked, top, N = 100L, L = 10L)
  expect_lte(restr$cutoff, 10L)
  expect_gte(restr$overlap, 1L)
})

test_that("X requires at least min_hits before a cutoff is eligible", {
  ranked <- paste0("g", 1:100)
  lone <- paste0("g", c(1, 80, 81, 82, 83)) # one lone top gene, rest deep

  x1 <- ohg_statistic(ranked, lone, N = 100L, X = 1L)
  x3 <- ohg_statistic(ranked, lone, N = 100L, X = 3L)

  expect_equal(x1$cutoff, 1L) # a single top gene can drive cutoff 1
  expect_gte(x3$overlap, 3L) # X = 3 forces at least 3 hits in the prefix
  expect_gte(x3$cutoff, 81L) # ... so the cutoff sits at the 3rd hit or deeper
})

test_that("ohg_enrichment caps cutoff_rank at L and records L_used/X_used", {
  set.seed(1)
  ranked <- paste0("g", 1:400)
  sets <- list(TOP = ranked[1:15], SPREAD = ranked[round(seq(1, 400, length.out = 25))])
  res <- ohg_enrichment_quiet(ranked, sets,
    max_cutoff_frac = 0.25, min_hits = 2L, n_perm = 500L, seed = 1
  )
  expect_true(all(res$cutoff_rank <= ceiling(0.25 * 400)))
  expect_equal(attr(res, "L_used"), ceiling(0.25 * 400))
  expect_equal(attr(res, "X_used"), 2L)
})

test_that("a bottom-only set is dropped under L but kept when unrestricted (up)", {
  set.seed(2)
  ranked <- paste0("g", 1:400)
  sets <- list(TOP = ranked[1:15], BOTTOM = ranked[360:385])

  restr <- ohg_enrichment_quiet(ranked, sets,
    direction = "up", max_cutoff_frac = 0.25, n_perm = 500L, seed = 1
  )
  full <- ohg_enrichment_quiet(ranked, sets,
    direction = "up", max_cutoff_frac = 1, n_perm = 500L, seed = 1
  )
  expect_false("BOTTOM" %in% restr$pathway) # no hits in the top 25% -> dropped
  expect_true("BOTTOM" %in% full$pathway) # unrestricted finds the deep bunch
  expect_true("TOP" %in% restr$pathway) # the concentrated set survives
})

test_that("CALIBRATION under L/X: observed and null share L => valid p (centred)", {
  # The §A.4 trap: if the null does not use the same L/X as the observed path, the
  # p-value is biased. Matched L/X must give a valid (level-controlling) p-value.
  set.seed(7)
  N <- 300L
  m <- 15L
  B <- 1500L
  L <- 75L
  null <- ohg_permutation_null(N, m, B, boundaries = seq_len(N), L = L, X = 1L)
  reps <- 500L
  p_vals <- vapply(seq_len(reps), function(i) {
    is_hit <- logical(N)
    is_hit[sample.int(N, m)] <- TRUE
    obs <- .mhg_core(is_hit, m, N, seq_len(N), L = L, X = 1L)$log_stat
    (1 + sum(null$log_stat_b <= obs)) / (1 + B)
  }, numeric(1))
  # mHG is discrete and L adds a small conservative atom, so test validity, not
  # exact uniformity: not anti-conservative, and roughly centred.
  expect_lte(mean(p_vals <= 0.05), 0.08)
  expect_gt(mean(p_vals), 0.45)
  expect_lt(mean(p_vals), 0.60)
})

test_that("max_cutoff_frac = 1 reproduces the unrestricted scan", {
  ranked <- paste0("g", 1:100)
  deep <- paste0("g", 40:55)
  a <- ohg_statistic(ranked, deep, N = 100L) # default L = N
  b <- ohg_statistic(ranked, deep, N = 100L, L = 100L) # explicit L = N
  expect_equal(a$log_stat, b$log_stat)
  expect_equal(a$cutoff, b$cutoff)
})
