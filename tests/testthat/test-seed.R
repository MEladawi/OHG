test_that("large integer seeds do not overflow set.seed", {
  ranked <- paste0("g", 1:200)
  sets <- list(A = ranked[1:12], B = ranked[5:25])
  big <- .Machine$integer.max - 3L # seed + m would overflow integer for m >= 4

  expect_no_error(
    ohg_enrichment(
      ranked, sets,
      n_perm = 200L, method = "permutation", seed = big
    )
  )
  # suppressWarnings: a tiny null can hit the permutation cap and emit the
  # legitimate "conservative lower-bound" warning; this test only asserts the
  # large seed does not abort the run.
  expect_no_error(
    suppressWarnings(
      ohg_enrichment(
        ranked, sets,
        n_perm = 200L, method = "adaptive", n_perm_max = 2000L, seed = big
      )
    )
  )
})

test_that("a keyed seed near the integer ceiling is still reproducible", {
  ranked <- paste0("g", 1:200)
  sets <- list(A = ranked[1:12], B = ranked[5:25])
  big <- .Machine$integer.max - 3L

  r1 <- ohg_enrichment(ranked, sets, n_perm = 300L, method = "permutation", seed = big)
  r2 <- ohg_enrichment(ranked, sets, n_perm = 300L, method = "permutation", seed = big)
  expect_equal(r1$p_value, r2$p_value)
})
