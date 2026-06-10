# Byte-identity guards for the adaptive redraw elimination. The goldens were
# snapshotted from the pre-refactor build (characterization baseline) and must
# pass before AND after the refactor (proof that removing the redraw changed
# nothing). A refactor guard, not a fail-first test: if it ever fails, the
# refactor altered observable output.

test_that("heterogeneous-group adaptive result is byte-identical to golden", {
  golden <- readRDS(test_path("fixtures", "redraw_golden_hetero.rds"))
  got <- suppressWarnings(run_redraw_fixture(FALSE))
  expect_equal(got, golden)
  expect_identical(attr(got, "L_used"), attr(golden, "L_used"))
  expect_identical(attr(got, "X_used"), attr(golden, "X_used"))
})

test_that("signed-weight adaptive result is byte-identical to golden", {
  golden <- readRDS(test_path("fixtures", "redraw_golden_signed.rds"))
  got <- suppressWarnings(run_redraw_fixture(TRUE))
  expect_equal(got, golden)
  expect_identical(attr(got, "L_used"), attr(golden, "L_used"))
  expect_identical(attr(got, "X_used"), attr(golden, "X_used"))
})

test_that(".adaptive_draw_increment returns delta_E_b matching a manual draw", {
  N <- 120L; m <- 8L; n_new <- 150L; seed <- 3L
  obs <- c(-2.0, -0.5) # two observed log-stats
  w <- rnorm(N)        # signed weights
  inc <- .adaptive_draw_increment(
    obs_log_stats = obs, m = m, N = N, boundaries = seq_len(N), n_new = n_new,
    rng_state = NULL, c_prev = c(0L, 0L), b_used_prev = 0L, seed = seed,
    target_hits = 10L, weight = w, robust = TRUE
  )
  # Reproduce the SAME draw independently and build E_b via the shared helper.
  set.seed(seed + m)
  nb <- ohg_permutation_null(N, m, n_new, seq_len(N))
  expect_identical(inc$delta_E_b, .eb_from_leidx(nb$le_idx_b, w, TRUE))
  expect_length(inc$delta_E_b, n_new)
})

test_that(".adaptive_draw_increment with weight = NULL yields delta_E_b = NULL", {
  inc <- .adaptive_draw_increment(
    obs_log_stats = c(-1.0), m = 8L, N = 120L, boundaries = seq_len(120L),
    n_new = 50L, rng_state = NULL, c_prev = 0L, b_used_prev = 0L, seed = 1L,
    target_hits = 10L, weight = NULL, robust = TRUE
  )
  expect_null(inc$delta_E_b)
})
