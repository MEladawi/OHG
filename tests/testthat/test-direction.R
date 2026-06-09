test_that("infer_direction: signed => both, non-negative/NULL => up", {
  expect_identical(infer_direction(c(3, 1, -1, -4)), "both")
  expect_identical(infer_direction(c(5, 4, 3, 1)), "up")
  expect_identical(infer_direction(NULL), "up")
})

test_that("explicit direction always overrides inference", {
  expect_identical(infer_direction(c(3, -1), supplied = "up"), "up")
  expect_identical(infer_direction(c(5, 4), supplied = "both"), "both")
})

test_that("down == up on the reversed list; both pools one adjustment", {
  set.seed(3)
  ranked <- paste0("g", 1:200)
  rs <- seq(100, -99) # signed, strictly decreasing
  sets <- list(TOP = ranked[1:12], BOTTOM = ranked[189:200])

  dn <- ohg_enrichment(ranked, sets,
    rank_stat = rs, direction = "down",
    n_perm = 500L, seed = 1
  )
  # The down run is the up run on reversed inputs: BOTTOM is strongly enriched down.
  expect_lt(dn$p_value[dn$pathway == "BOTTOM"][1], 0.05)

  both <- ohg_enrichment(ranked, sets,
    rank_stat = rs, direction = "both",
    n_perm = 500L, seed = 1
  )
  expect_true(all(c("up", "down") %in% both$direction))
  expect_equal(both$p_adjust, stats::p.adjust(both$p_value, "BH"))
})

test_that("inferred direction: signed => direction column, non-negative => up only", {
  set.seed(5)
  ranked <- paste0("g", 1:120)
  res_both <- ohg_enrichment(ranked, list(S = ranked[1:10]),
    rank_stat = seq(60, -59), n_perm = 300L, seed = 1
  )
  expect_true("direction" %in% names(res_both))

  res_up <- ohg_enrichment(ranked, list(S = ranked[1:10]),
    rank_stat = seq(120, 1), n_perm = 300L, seed = 1
  )
  expect_false("direction" %in% names(res_up))
})
