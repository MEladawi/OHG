test_that("infer_direction: signed => both, non-negative/NULL => up", {
  expect_identical(infer_direction(c(3, 1, -1, -4)), "both")
  expect_identical(infer_direction(c(5, 4, 3, 1)), "up")
  expect_identical(infer_direction(NULL), "up")
})

test_that("infer_direction: all-negative => down (signal sits at the bottom)", {
  # A signed statistic that never crosses zero (every gene down-regulated) puts
  # its strongest genes at the bottom of the descending sort, so the meaningful
  # tail is the bottom -- defaulting to "up" would test the weakest genes.
  expect_identical(infer_direction(c(-1, -2, -4)), "down")
  expect_identical(infer_direction(c(0, -1, -3)), "down") # zeros don't count as positive
  expect_identical(infer_direction(c(0, 0, 0)), "up") # no negatives => still "up"
})

test_that("inferred all-negative direction finds bottom enrichment by default", {
  ranked <- paste0("g", 1:200)
  rs <- seq(-1, -200) # strictly decreasing, all negative => default "down"
  sets <- list(TOP = ranked[1:12], BOTTOM = ranked[189:200])
  res <- ohg_enrichment_quiet(ranked, sets, rank_stat = rs, n_perm = 500L, seed = 1)
  # default inference must pick the down tail, where BOTTOM is strongly enriched
  expect_false("up" %in% res$direction)
  expect_lt(res$p_value[res$pathway == "BOTTOM"][1], 0.05)
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

  dn <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "down",
    n_perm = 500L, seed = 1
  )
  # The down run is the up run on reversed inputs: BOTTOM is strongly enriched down.
  expect_lt(dn$p_value[dn$pathway == "BOTTOM"][1], 0.05)

  both <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "both",
    n_perm = 500L, seed = 1
  )
  expect_true(all(c("up", "down") %in% both$direction))
  # One pooled BH adjustment over BOTH tails, with the family counting every
  # tested hypothesis (each kept pathway x each direction). TOP enriches up but
  # not down, and BOTTOM the reverse, so the TOP-down and BOTTOM-up tails are
  # tested null results (p = 1) dropped from the output -- they still belong in
  # the family, hence n = 2 * kept, not the emitted-row count.
  kept <- length(prune_gene_sets(coerce_gene_sets(sets), ranked, 3L, length(ranked))$kept)
  expect_equal(both$p_adjust, stats::p.adjust(both$p_value, "BH", n = 2L * kept))
})

test_that("both-direction rows match separate up/down runs (shared-null safety)", {
  ranked <- paste0("g", 1:200)
  rs <- seq(100, -99) # distinct => up and down share boundaries
  sets <- list(TOP = ranked[1:12], BOTTOM = ranked[189:200])

  both <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "both", weight = abs(rs), n_perm = 1200L, seed = 1
  )
  up <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "up", weight = abs(rs), n_perm = 1200L, seed = 1
  )
  dn <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "down", weight = abs(rs), n_perm = 1200L, seed = 1
  )

  bu <- both[both$direction == "up", ]
  bd <- both[both$direction == "down", ]
  expect_equal(bu$p_value[order(bu$pathway)], up$p_value[order(up$pathway)])
  expect_equal(bd$p_value[order(bd$pathway)], dn$p_value[order(dn$pathway)])
  # effect axis too: null reuse must not corrupt direction-specific E_obs / NLES
  expect_equal(bd$E_obs[order(bd$pathway)], dn$E_obs[order(dn$pathway)])
  expect_equal(bd$NLES_signed[order(bd$pathway)], dn$NLES_signed[order(dn$pathway)])
})

test_that("collapse_both keeps the more significant direction with a x2 penalty", {
  set.seed(8)
  ranked <- paste0("g", 1:200)
  rs <- seq(100, -99)
  sets <- list(TOP = ranked[1:12], BOTTOM = ranked[189:200])

  both <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "both", n_perm = 500L, seed = 1
  )
  coll <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, direction = "both", collapse_both = TRUE,
    n_perm = 500L, seed = 1
  )

  # one row per pathway after collapsing
  expect_equal(nrow(coll), length(unique(coll$pathway)))
  expect_setequal(coll$pathway, c("TOP", "BOTTOM"))

  # collapsed p_value is the per-pathway minimum, x2-penalized and capped at 1
  best <- dplyr::summarise(both, p = min(p_value), .by = pathway)
  expected <- pmin(1, best$p[match(coll$pathway, best$pathway)] * 2)
  expect_equal(coll$p_value, expected)
})

test_that("inferred direction: signed => direction column, non-negative => up only", {
  set.seed(5)
  ranked <- paste0("g", 1:120)
  res_both <- ohg_enrichment_quiet(ranked, list(S = ranked[1:10]),
    rank_stat = seq(60, -59), n_perm = 300L, seed = 1
  )
  expect_true("direction" %in% names(res_both))

  res_up <- ohg_enrichment_quiet(ranked, list(S = ranked[1:10]),
    rank_stat = seq(120, 1), n_perm = 300L, seed = 1
  )
  expect_false("direction" %in% names(res_up))
})
