#' Permutation null for the mHG statistic
#'
#' Draws `B` random position sets `R_b` of size `m` uniformly from the universe
#' `seq_len(N)` (a random placement of `m` hits among `N` positions) and evaluates
#' the identical tie-aware mHG statistic on each. Because every draw is a size-`m`
#' subset of the universe, each permutation has exactly `m` hits and a leading edge
#' of size at least 1 — zero-overlap permutations cannot occur. Position draws and
#' `log_stat_b` depend only on `(N, m, boundaries)` and so are portable across
#' datasets; the per-permutation leading-edge indices are returned for the
#' call-specific effect-size null.
#'
#' @param N Universe size.
#' @param m Number of hits (the set size).
#' @param B Number of permutations.
#' @param boundaries Integer tie-block boundaries; use `seq_len(N)` for distinct ranks.
#'
#' @return A list with `log_stat_b` (length-`B` numeric), `overlap_b` (length-`B`
#'   integer), and `le_idx_b` (length-`B` list of leading-edge position vectors).
#'
#' @examples
#' ohg_permutation_null(N = 100L, m = 8L, B = 200L)
#'
#' @export
ohg_permutation_null <- function(N, m, B, boundaries = seq_len(N)) {
  stat_one <- function(pos) {
    is_hit <- logical(N)
    is_hit[pos] <- TRUE
    .mhg_core(is_hit, m, N, boundaries) # shared kernel (see statistic.R)
  }

  draws <- lapply(seq_len(B), function(b) stat_one(sort(sample.int(N, m))))
  list(
    log_stat_b = vapply(draws, `[[`, numeric(1), "log_stat"),
    overlap_b = vapply(draws, function(d) as.integer(d$overlap), integer(1)),
    le_idx_b = lapply(draws, `[[`, "le_idx")
  )
}

#' Build one permutation null per distinct set size, reproducibly
#'
#' Constructs `ohg_permutation_null()` once per distinct `m` (descending, for load
#' balance). When `n_cores > 1` and both `furrr` and `future` are installed, the
#' builds run in parallel; each `m` seeds an independent stream keyed to `m` (not
#' execution order), so results are identical to the sequential path for the same
#' `seed`.
#'
#' @param ms Integer set sizes (duplicates collapsed internally).
#' @param N Universe size.
#' @param B Permutations per `m`.
#' @param boundaries Integer tie-block boundaries.
#' @param n_cores Worker count; `> 1` parallelizes when furrr/future are available.
#' @param seed Integer base seed or `NULL`.
#'
#' @return A named list of nulls, keyed by `as.character(m)`.
#'
#' @keywords internal
#' @export
build_nulls <- function(ms, N, B, boundaries, n_cores = 1L, seed = NULL) {
  ms <- sort(unique(ms), decreasing = TRUE)
  one <- function(mm) {
    if (!is.null(seed)) set.seed(seed + mm) # stream keyed to m, not order
    ohg_permutation_null(N, mm, B, boundaries)
  }
  use_par <- n_cores > 1L &&
    requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)
  res <- if (use_par) {
    furrr::future_map(ms, one, .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(ms, one)
  }
  stats::setNames(res, as.character(ms))
}
