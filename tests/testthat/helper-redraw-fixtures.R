# Deterministic fixtures for the adaptive redraw-elimination byte-identity
# goldens. Shared by the golden-generation step (run on the pre-refactor build)
# and by test-redraw-byte-identical.R (run on every build).

# Build the inputs. `signed_weight = TRUE` injects negative weights to guard the
# abs contract in .eb_from_leidx(); FALSE uses magnitudes.
redraw_inputs <- function(signed_weight = FALSE) {
  set.seed(20260610L)
  N <- 240L
  genes <- paste0("g", seq_len(N))
  rank_stat <- rnorm(N) # continuous -> few ties
  weight <- if (signed_weight) rnorm(N) else abs(rnorm(N))
  top <- genes[order(rank_stat, decreasing = TRUE)] # best-ranked first
  sets <- list(
    HITS_TOP    = top[1:10],            # size 10, strongly enriched up -> caps
    RANDOM_10   = sample(genes, 10),    # size 10, ~random -> resolves early
    HITS_BOTTOM = top[(N - 9):N],       # size 10, enriched only in 'down'
    MIXED_15    = c(top[1:3], sample(genes, 12)) # size 15, moderate
  )
  list(genes = genes, sets = sets, rank_stat = rank_stat, weight = weight)
}

# Run the adaptive engine on a fixture. Small n_perm_max forces a cap cheaply;
# both goldens MUST use identical run parameters on the current and refactored
# builds so the diff varies only the implementation.
run_redraw_fixture <- function(signed_weight = FALSE) {
  fx <- redraw_inputs(signed_weight)
  ohg_enrichment(
    fx$genes, fx$sets,
    rank_stat = fx$rank_stat, weight = fx$weight,
    direction = "both", method = "adaptive",
    n_perm = 50L, min_perm_nles = 50L, min_nles_support = 10L,
    n_perm_max = 200L, target_hits = 10L,
    seed = 99L, n_cores = 1L
  )
}
