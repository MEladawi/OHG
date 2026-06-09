test_that("read_gmt parses name, description, and genes; drops blanks", {
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c(
    "SET_A\tdescription A\tg1\tg2\tg3",
    "SET_B\thttp://example.org\tg2\tg4"
  ), gmt)

  res <- read_gmt(gmt)

  expect_type(res, "list")
  expect_named(res, c("SET_A", "SET_B"))
  expect_identical(res$SET_A, c("g1", "g2", "g3"))
  expect_identical(res$SET_B, c("g2", "g4"))
  expect_identical(attr(res, "descriptions")[["SET_A"]], "description A")
})

test_that("read_gmt tolerates trailing empty fields and blank lines", {
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c("SET_A\tdesc\tg1\tg2\t\t", "", "SET_B\tdesc\tg3"), gmt)
  res <- read_gmt(gmt)
  expect_identical(res$SET_A, c("g1", "g2"))
  expect_named(res, c("SET_A", "SET_B"))
})

test_that("coerce_gene_sets unifies list, .gmt path, and GeneSetCollection", {
  lst <- list(SET_A = c("g1", "g2", "g3"), SET_B = c("g2", "g4"))
  expect_identical(coerce_gene_sets(lst), lapply(lst, unique))

  gmt <- tempfile(fileext = ".gmt")
  writeLines(c("SET_A\tdesc\tg1\tg2\tg3", "SET_B\tdesc\tg2\tg4"), gmt)
  from_path <- coerce_gene_sets(gmt)
  expect_identical(unname(from_path[["SET_A"]]), c("g1", "g2", "g3"))

  skip_if_not_installed("GSEABase")
  gs_a <- GSEABase::GeneSet(c("g1", "g2", "g3"), setName = "SET_A")
  gs_b <- GSEABase::GeneSet(c("g2", "g4"), setName = "SET_B")
  gsc <- GSEABase::GeneSetCollection(list(gs_a, gs_b))
  from_gsc <- coerce_gene_sets(gsc)
  expect_named(from_gsc, c("SET_A", "SET_B"))
  expect_identical(sort(from_gsc[["SET_A"]]), sort(c("g1", "g2", "g3")))
})

test_that("coerce_gene_sets errors on unnamed list and bad type", {
  expect_error(coerce_gene_sets(list(c("g1"), c("g2"))), "named")
  expect_error(coerce_gene_sets(42), "named list|GeneSetCollection|.gmt")
})
