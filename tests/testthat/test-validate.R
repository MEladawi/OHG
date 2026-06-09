test_that("de-dup keeps best rank and warns", {
  expect_warning(
    v <- validate_inputs(c("g1", "g2", "g1"), rank_stat = c(3, 2, 1), weight = NULL),
    "duplicate"
  )
  expect_identical(v$ranked_genes, c("g1", "g2"))
  expect_identical(v$rank_stat, c(3, 2))
})

test_that("unsorted rank_stat is sorted descending; genes and weight follow", {
  v <- validate_inputs(
    c("g1", "g2", "g3"),
    rank_stat = c(1, 5, 3), weight = c(10, 20, 30)
  )
  expect_identical(v$ranked_genes, c("g2", "g3", "g1")) # by rank_stat 5, 3, 1
  expect_identical(v$rank_stat, c(5, 3, 1))
  expect_identical(v$weight, c(20, 30, 10)) # weight tracks the gene order
})

test_that("tie fraction is reported when rank_stat has ties", {
  expect_message(
    validate_inputs(paste0("g", 1:4), rank_stat = c(5, 5, 3, 1), weight = NULL),
    "ties"
  )
  expect_silent(
    validate_inputs(paste0("g", 1:4), rank_stat = c(5, 4, 3, 1), weight = NULL)
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
