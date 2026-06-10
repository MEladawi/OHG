# Byte-identity guard for the no-ties mHG kernel fast path. The golden was
# snapshotted from the pre-refactor build. The O(m) fast path (tie-free) and the
# O(N) tie-aware fallback must both reproduce it exactly -- every log_stat, overlap,
# and leading-edge index vector. If this fails, the refactor changed observable
# kernel output.
test_that("ohg_permutation_null is byte-identical to the pre-refactor kernel golden", {
  golden <- readRDS(test_path("fixtures", "kernel_golden.rds"))
  got <- kernel_fixture()
  expect_identical(got$tie_free, golden$tie_free) # the fast path
  expect_identical(got$tied, golden$tied)         # the tie-aware fallback
})

test_that(".mhg_core_pos matches .mhg_core on random tie-free inputs (incl. X, edges)", {
  set.seed(11L)
  N <- 400L
  for (rep in seq_len(40L)) {
    m <- sample.int(150L, 1L)
    L <- sample(c(40L, 100L, N), 1L)
    X <- sample(1:4, 1L)
    pos <- sort(sample.int(N, m))
    is_hit <- logical(N)
    is_hit[pos] <- TRUE
    a <- .mhg_core(is_hit, m, N, seq_len(N), L = L, X = X)
    b <- .mhg_core_pos(pos, m, N, L = L, X = X)
    expect_identical(b, a)
  }
  # explicit edge: no hits in the top L -> log_stat 0, empty leading edge
  pos <- c(300L, 350L, 399L)
  is_hit <- logical(N); is_hit[pos] <- TRUE
  expect_identical(
    .mhg_core_pos(pos, 3L, N, L = 50L, X = 1L),
    .mhg_core(is_hit, 3L, N, seq_len(N), L = 50L, X = 1L)
  )
  # explicit edge: empty hit set
  expect_identical(
    .mhg_core_pos(integer(0), 0L, N, L = 100L, X = 1L),
    .mhg_core(logical(N), 0L, N, seq_len(N), L = 100L, X = 1L)
  )
})
