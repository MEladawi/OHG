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
