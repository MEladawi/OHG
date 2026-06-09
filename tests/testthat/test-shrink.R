test_that("ohg_shrink_lfc shrinks from (LFC, SE) only", {
  skip_if_not_installed("ashr")
  set.seed(1)
  lfc <- c(3, 0.1, -2, 0.05, 5)
  se <- c(0.2, 2.0, 0.3, 3.0, 0.1)
  out <- ohg_shrink_lfc(lfc, se)
  expect_length(out, length(lfc))
  expect_true(all(is.finite(out)))
  # The high-SE estimate shrinks proportionally more than the low-SE estimate.
  expect_lt(abs(out[2]) / abs(lfc[2]), abs(out[1]) / abs(lfc[1]))
})

test_that("clean error when ashr is absent", {
  if (requireNamespace("ashr", quietly = TRUE)) {
    skip("ashr installed")
  }
  expect_error(ohg_shrink_lfc(c(1, 2), c(1, 1)), "ashr")
})

test_that("length-mismatch errors", {
  expect_error(ohg_shrink_lfc(c(1, 2, 3), c(1, 1)), "same length")
})
