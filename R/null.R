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
  # Tie-free (continuous metric): the sampled hit positions ARE the eligible
  # cutoffs, so the O(m) kernel takes them directly -- never building the length-N
  # is_hit. Decided once here, not per draw. Tied metrics fall back to .mhg_core.
  # Both paths consume the identical sample.int() stream, so output is byte-for-byte
  # unchanged and reproducible across cores for a fixed seed.
  if (identical(boundaries, seq_len(N))) {
    draws <- lapply(seq_len(B), function(b) {
      .mhg_core_pos(sort(sample.int(N, m)), m, N, L = L, X = X)
    })
  } else {
    stat_one <- function(pos) {
      is_hit <- logical(N)
      is_hit[pos] <- TRUE
      .mhg_core(is_hit, m, N, boundaries, L = L, X = X) # shared kernel (statistic.R)
    }
    draws <- lapply(seq_len(B), function(b) stat_one(sort(sample.int(N, m))))
  }
  list(
    log_stat_b = vapply(draws, `[[`, numeric(1), "log_stat"),
    overlap_b = vapply(draws, function(d) as.integer(d$overlap), integer(1)),
    le_idx_b = lapply(draws, `[[`, "le_idx")
  )
}

# Per-size RNG seed keyed to m, overflow-safe. `seed + m` on two integers can
# exceed .Machine$integer.max and silently overflow to NA, which then makes
# set.seed(NA) abort ("supplied seed is not a valid integer") for large seeds.
# Add in double precision and wrap into the valid integer range first; the wrap
# is deterministic, so the per-m streams stay reproducible. Internal helper.
.keyed_seed <- function(seed, m) {
  as.integer((as.double(seed) + as.double(m)) %% .Machine$integer.max)
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
    if (!is.null(seed)) set.seed(.keyed_seed(seed, mm)) # stream keyed to m, not order
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

# Draw the NEXT `n_new` permutations of the size-`m` null, continuing one resumable
# RNG stream (seeded `seed + m` on the first call, then resumed from `rng_state`) so
# the rounds never redraw work already done -- the cumulative stream is identical to
# a single `set.seed(seed + m); draw(b_used)`, which is what the finalize pass uses
# to recover leading edges. Returns each pathway's exceedance increment `delta_c`,
# the global draw index at which it first reaches `target_hits` (`L_hit`, given the
# carried `c_prev`/`b_used_prev`, NA if not yet), and the advanced `rng_state`. All
# payloads are tiny integer vectors, so a worker ships them back cheaply. Internal.
.adaptive_draw_increment <- function(obs_log_stats, m, N, boundaries, n_new,
                                     rng_state, c_prev, b_used_prev, seed,
                                     target_hits, L = N, X = 1L,
                                     weight = NULL, robust = TRUE) {
  if (is.null(rng_state)) {
    set.seed(.keyed_seed(seed, m)) # first draw: open the stream keyed to m
  } else {
    assign(".Random.seed", rng_state, envir = .GlobalEnv) # resume where we stopped
  }
  nb <- ohg_permutation_null(N, m, n_new, boundaries, L = L, X = X)
  lsb <- nb$log_stat_b
  new_state <- get(".Random.seed", envir = .GlobalEnv)
  delta_c <- integer(length(obs_log_stats))
  L_hit <- rep(NA_integer_, length(obs_log_stats))
  for (i in seq_along(obs_log_stats)) {
    hits <- lsb <= obs_log_stats[i]
    delta_c[i] <- sum(hits)
    w <- which(c_prev[i] + cumsum(hits) >= target_hits)[1L]
    if (!is.na(w)) L_hit[i] <- as.integer(b_used_prev + w)
  }
  # Group-level NLES null increment: center(abs(weight)[le_idx]) per fresh draw,
  # via the single shared E_b definition. NULL weight -> no NLES, skip the work.
  delta_E_b <- if (is.null(weight)) {
    NULL
  } else {
    .eb_from_leidx(nb$le_idx_b, weight, robust)
  }
  list(
    delta_c = delta_c, L_hit = L_hit, rng_state = new_state,
    delta_E_b = delta_E_b
  )
}
