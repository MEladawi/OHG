test_that("ohg_winsorize caps both tails symmetrically at quantile(abs(x), p)", {
  set.seed(1)
  x <- c(rnorm(100), 50, -60) # two wild tails
  out <- ohg_winsorize(x, p = 0.95)
  cap <- attr(out, "cap")

  expect_equal(cap, unname(stats::quantile(abs(x), 0.95)))
  expect_length(out, length(x))
  expect_lte(max(abs(out)), cap + 1e-9) # nothing exceeds the cap
  expect_equal(out[101], cap) # +50 -> +cap
  expect_equal(out[102], -cap) # -60 -> -cap (sign kept, not abs)
  expect_equal(attr(out, "n_capped"), sum(abs(x) > cap))
})

test_that("ohg_winsorize keeps sign and leaves the body unchanged at p = 1", {
  x <- c(-3, -0.1, 0.2, 2)
  out <- ohg_winsorize(x, p = 1) # cap = max(abs(x)) => nothing clipped
  expect_equal(as.numeric(out), x)
  expect_equal(attr(out, "n_capped"), 0L)
})

test_that("ohg_winsorize tames Inf and validates its arguments", {
  out <- ohg_winsorize(c(1, 2, Inf, -3), p = 0.75)
  expect_true(all(is.finite(out))) # Inf clipped to the cap
  expect_error(ohg_winsorize(c(1, 2), p = 0), "in \\(0, 1\\]")
  expect_error(ohg_winsorize(c(1, 2), p = 1.5), "in \\(0, 1\\]")
  expect_error(ohg_winsorize("a"), "numeric")
})
