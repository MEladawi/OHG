test_that("de-dup keeps best rank and warns", {
  expect_warning(
    v <- validate_inputs(c("g1", "g2", "g1"), rank_stat = c(3, 2, 1), weight = NULL),
    "duplicate"
  )
  expect_identical(v$ranked_genes, c("g1", "g2"))
  expect_identical(v$rank_stat, c(3, 2))
})

test_that("non-increasing rank_stat is asserted", {
  expect_error(
    validate_inputs(c("g1", "g2"), rank_stat = c(1, 5), weight = NULL),
    "non-increasing|sorted"
  )
})

test_that("weight guard: hard error on non-finite, soft warn on extreme range", {
  expect_error(validate_inputs(c("a", "b"), NULL, weight = c(1, Inf)), "finite")
  expect_error(validate_inputs(c("a", "b"), NULL, weight = c(1, NA)), "finite")
  expect_warning(
    validate_inputs(paste0("g", 1:4), NULL, weight = c(1e9, 1, 1, 1)),
    "dynamic range|extreme"
  )
})

test_that("p_adjust_method validated against stats::p.adjust.methods", {
  expect_error(
    validate_inputs(c("a", "b"), NULL, NULL, p_adjust_method = "nope"),
    "p_adjust_method"
  )
})

test_that("set_size m = |T intersect U|; min_set_size gates inclusion not NLES", {
  v <- validate_inputs(paste0("g", 1:10), NULL, NULL)
  sets <- list(BIG = paste0("g", 1:5), TINY = c("g1", "zzz"))
  pruned <- prune_gene_sets(sets, universe = v$ranked_genes, min_set_size = 3L, N = 10L)
  expect_named(pruned$kept, "BIG")
  expect_identical(pruned$kept$BIG$m, 5L)
  expect_true("TINY" %in% pruned$dropped$pathway) # m=1 < 3 dropped from inclusion
})
