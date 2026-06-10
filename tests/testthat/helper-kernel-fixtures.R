# Deterministic kernel fixtures for the no-ties fast-path byte-identity golden.
# Exercises ohg_permutation_null on BOTH the tie-free case (boundaries == seq_len(N),
# the O(m) fast path) and a tie-aware case (the O(N) fallback), across several set
# sizes and X values plus edge cases (a set with no hits in the top L). The golden
# is snapshotted on the pre-refactor build; the fast path must reproduce it byte for
# byte. Standalone calls use the default RNG kind; set.seed before each draw set.
kernel_fixture <- function() {
  N <- 500L
  Lq <- as.integer(ceiling(0.25 * N)) # 125

  # --- tie-free (fast path): boundaries == seq_len(N) ---
  tf_cfg <- list(
    c(m = 10L, X = 1L), c(m = 47L, X = 1L), c(m = 47L, X = 3L),
    c(m = 120L, X = 1L), c(m = 3L, X = 1L), c(m = 250L, X = 5L)
  )
  tie_free <- lapply(tf_cfg, function(p) {
    set.seed(100L + p[["m"]] + p[["X"]])
    ohg_permutation_null(N, p[["m"]], 200L, seq_len(N), L = Lq, X = p[["X"]])
  })
  names(tie_free) <- vapply(tf_cfg, function(p) paste(p[["m"]], p[["X"]], sep = "_"), "")

  # tiny-L edge: many draws place all hits below L -> log_stat 0, le_idx empty
  set.seed(999L)
  tie_free[["tinyL"]] <- ohg_permutation_null(N, 8L, 200L, seq_len(N), L = 20L, X = 1L)

  # --- tie-aware (fallback): boundaries from a heavily-tied rank vector ---
  set.seed(7L)
  rs <- sort(sample.int(60L, N, replace = TRUE), decreasing = TRUE) # many ties
  tb <- tie_boundaries(rs, n = N)
  tied <- lapply(c(10L, 47L, 120L), function(m) {
    set.seed(200L + m)
    ohg_permutation_null(N, m, 200L, tb, L = Lq, X = 1L)
  })
  names(tied) <- as.character(c(10L, 47L, 120L))

  list(tie_free = tie_free, tied = tied, boundaries_tied = tb)
}
