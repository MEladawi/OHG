#' Permutation null for the mHG statistic
#'
#' Draws `B` random position sets `R_b` of size `m` uniformly from the universe
#' `seq_len(N)` (a random placement of `m` hits among `N` positions) and evaluates
#' the identical tie-aware mHG statistic on each. Because every draw is a size-`m`
#' subset of the universe, each permutation has exactly `m` hits and a leading edge
#' of size at least 1 -- zero-overlap permutations cannot occur. Position draws and
#' `log_stat_b` depend only on `(N, m, boundaries)` and so are portable across
#' datasets; the per-permutation leading-edge indices are returned for the
#' call-specific effect-size null.
#'
#' @param N Universe size.
#' @param m Number of hits (the set size).
#' @param B Number of permutations.
#' @param boundaries Integer tie-block boundaries; use `seq_len(N)` for distinct ranks.
#' @param L,X XL-mHG restriction passed to the shared kernel (largest cutoff `L`,
#'   minimum prefix hits `X`); must match the observed path for a valid p-value.
#'
#' @return A list with `log_stat_b` (length-`B` numeric), `overlap_b` (length-`B`
#'   integer), and `le_idx_b` (length-`B` list of leading-edge position vectors).
#'
#' @examples
#' ohg_permutation_null(N = 100L, m = 8L, B = 200L)
#'
#' @export
ohg_permutation_null <- function(N, m, B, boundaries = seq_len(N), L = N, X = 1L) {
  stat_one <- function(pos) {
    is_hit <- logical(N)
    is_hit[pos] <- TRUE
    .mhg_core(is_hit, m, N, boundaries, L = L, X = X) # shared kernel (statistic.R)
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
#' @param L,X XL-mHG restriction forwarded to [ohg_permutation_null()].
#'
#' @return A named list of nulls, keyed by `as.character(m)`.
#'
#' @keywords internal
#' @export
build_nulls <- function(ms, N, B, boundaries, n_cores = 1L, seed = NULL,
                        L = N, X = 1L) {
  ms <- sort(unique(ms), decreasing = TRUE)
  one <- function(mm) {
    if (!is.null(seed)) set.seed(seed + mm) # stream keyed to m, not order
    ohg_permutation_null(N, mm, B, boundaries, L = L, X = X)
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

# Besag-Clifford (1991) sequential Monte-Carlo p-value for ONE pathway, read off a
# stream of null log-statistics. `obs` is the observed log-statistic; an exceedance
# is a draw at least as enriched (`log_stat_b <= obs`). If the h-th exceedance is
# reached at draw L, p_hat = h / L (resolved); otherwise the stream is capped and
# p_hat = (c + 1) / (b_max + 1), a conservative bound flagged resolution_limited.
# Do NOT replace this with (1 + c)/(1 + B_used) under adaptive stopping -- that
# estimator is only valid for a fixed, pre-chosen B. Internal kernel helper.
.bc_finalize <- function(log_stat_b, obs, h, b_max) {
  cs <- cumsum(log_stat_b <= obs)
  c_final <- cs[length(cs)]
  if (c_final >= h) {
    list(
      p_hat = h / which(cs >= h)[1L],
      n_exceed = as.integer(c_final), resolution_limited = FALSE
    )
  } else {
    list(
      p_hat = (c_final + 1) / (b_max + 1),
      n_exceed = as.integer(c_final), resolution_limited = TRUE
    )
  }
}

# Adaptive sequential null for ALL pathways of one (set size, orientation). Draws a
# baseline `n_perm0`, then expands the SAME resumable RNG stream in geometric
# batches (doubling) until the strongest candidate accumulates `target_hits`
# exceedances or the stream reaches `n_perm_max`; every pathway reads its own
# (c, L) off the shared stream. Reusing one stream makes the result independent of
# how the draws are batched. Returns per-pathway p_hat/n_exceed/resolution_limited,
# the shared `n_perm_used`, and the accumulated `le_idx_b` buffer for NLES.
# Internal kernel helper (cf. .mhg_core); seeding/scoping is owned by the caller.
.adaptive_null <- function(obs_log_stats, m, N, boundaries,
                           n_perm0, n_perm_max, target_hits, seed, L = N, X = 1L) {
  if (!is.null(seed)) set.seed(seed + m) # resumable stream keyed to m, not order
  log_stat_b <- numeric(0)
  le_idx_b <- list()
  b_used <- 0L
  repeat {
    batch <- if (b_used == 0L) n_perm0 else b_used # geometric doubling
    batch <- min(batch, n_perm_max - b_used)
    draw <- ohg_permutation_null(N, m, batch, boundaries, L = L, X = X)
    log_stat_b <- c(log_stat_b, draw$log_stat_b)
    le_idx_b <- c(le_idx_b, draw$le_idx_b)
    b_used <- length(log_stat_b)
    c_i <- vapply(obs_log_stats, function(o) sum(log_stat_b <= o), integer(1))
    if (b_used >= n_perm_max || all(c_i >= target_hits)) break
  }
  fin <- lapply(obs_log_stats, function(o) {
    .bc_finalize(log_stat_b, o, target_hits, n_perm_max)
  })
  list(
    p_hat = vapply(fin, `[[`, numeric(1), "p_hat"),
    n_exceed = vapply(fin, `[[`, integer(1), "n_exceed"),
    resolution_limited = vapply(fin, `[[`, logical(1), "resolution_limited"),
    n_perm_used = b_used,
    le_idx_b = le_idx_b
  )
}
