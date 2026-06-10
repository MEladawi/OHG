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
    set.seed(seed + m) # first draw: open the stream keyed to m
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

# One-sided Clopper-Pearson (Beta) LOWER bound on a pathway's exceedance rate given
# `c` exceedances in `n` draws -- the smallest true p consistent with the draws at
# the simultaneous confidence carried in `gamma`. This is the "best case" a pathway
# could still resolve to: if even this optimistic p is non-significant after BH, no
# further draws can rescue it (DECIDED-ACCEPT). `c == 0` -> 0 (no exceedances yet
# could still be highly significant; never auto-accepted, ref amendment section 7).
# `gamma` is the per-bound error budget, split simultaneously across all m pathways
# and the lower side: `gamma = gate_conf / (2 * m)`. Internal helper.
.p_lower_cp <- function(c, n, gamma) {
  ifelse(c == 0L, 0, stats::qbeta(gamma, c, pmax(n - c + 1L, 1L)))
}

# Each pathway's contribution to the global BH multiset used by the accept gate,
# plus the mask of pathways still ramping. `resolved` -> h / L (locked once the
# target_hits-th exceedance is reached); `accepted`/`capped` -> their fixed reported
# estimate; a still-`ramping` pathway contributes its Clopper-Pearson lower bound
# (`.p_lower_cp`) -- its best achievable p. Feeding the lower bound to the BH step-up
# makes "stop ramping" the rigorous DECIDED-ACCEPT test: BH-adjusting this multiset
# and accepting where padj(pLo) > gate_alpha is exactly `pLo > tau_hi` (the highest
# threshold the future could produce), so acceptance is safe. Internal helper.
.group_p_opt <- function(g, target_hits, b_max, gamma) {
  st <- g$status
  p <- numeric(g$n)
  res <- st == "resolved"
  acc <- st == "accepted"
  cap <- st == "capped"
  rmp <- st == "ramping"
  p[res] <- target_hits / g$L_hit[res]
  p[acc] <- (g$c[acc] + 1) / (g$b_used + 1)
  p[cap] <- (g$c[cap] + 1) / (b_max + 1)
  p[rmp] <- .p_lower_cp(g$c[rmp], g$b_used, gamma)
  list(p_opt = p, ramping = rmp)
}

# The multiplicity gate. Given the global provisional p multiset (best case for the
# ramping pathways) and which entries are ramping, BH-adjust and stop ramping any
# pathway whose best-case adjusted p already exceeds `alpha`. Because every ramping
# pathway is scored at its floor here, each one's adjusted p is the smallest it can
# ever attain (own p minimal, competition maximal), so gating on `padj > alpha`
# never drops a pathway that could become significant. Internal helper.
.gate_decisions <- function(p_opt, ramping, alpha) {
  padj <- stats::p.adjust(p_opt, method = "BH")
  ramping & (padj > alpha)
}
