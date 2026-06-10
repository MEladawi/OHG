#' Minimum-hypergeometric (mHG) statistic in log space
#'
#' Computes the most-enriched prefix of `ranked_genes` against `T_eff`, working
#' entirely in natural-log space. The minimum upper-tail hypergeometric
#' probability is attained immediately after a hit, so evaluation is restricted
#' to hit-carrying tie-block boundaries (see [tie_boundaries()]); with distinct
#' ranks this is one evaluation per hit. Raw probabilities are never formed, so
#' the path is underflow-safe.
#'
#' @param ranked_genes Character vector, most important first (length `N`).
#' @param T_eff Character vector of pathway genes, already intersected with the
#'   universe and de-duplicated.
#' @param N Universe size (equal to `length(ranked_genes)`).
#' @param boundaries Optional integer tie-block boundaries from [tie_boundaries()].
#'   When `NULL`, distinct ranks are assumed (every position is a boundary).
#' @param L Largest cutoff considered (XL-mHG `L`): boundaries deeper than `L` are
#'   not eligible for the optimum, so `cutoff` never exceeds `L`. Default `N` (no
#'   restriction). A set with no hits at or above rank `L` returns `log_stat = 0`.
#' @param X Minimum prefix hits before a cutoff is eligible (XL-mHG `X`): only
#'   boundaries with cumulative hits `>= X` compete. Default `1` (off).
#'
#' @return A list with components `log_stat` (natural log of the mHG statistic),
#'   `cutoff` (optimal prefix size, or `NA`), `overlap` (hits inside the leading
#'   edge), and `le_idx` (integer positions of the leading-edge hits).
#'
#' @examples
#' ohg_statistic(paste0("g", 1:10), c("g1", "g2", "g5"), N = 10L)
#'
#' @export
ohg_statistic <- function(ranked_genes, T_eff, N, boundaries = NULL, L = N, X = 1L) {
  m <- length(T_eff)
  is_hit <- ranked_genes %in% T_eff
  if (is.null(boundaries)) boundaries <- seq_len(N) # distinct ranks
  .mhg_core(is_hit, m, N, boundaries, L = L, X = X)
}

# Internal mHG kernel shared by ohg_statistic() and ohg_permutation_null().
# Given a length-N hit indicator and the candidate cutoff `boundaries`, evaluates
# the upper-tail hypergeometric only at boundaries that add a new hit (constant q
# with larger k can only worsen the tail; statistic plan section 2.3) and returns
# the most-enriched prefix in natural-log space. The single source of truth for
# the statistic; both the observed and permutation paths call it. The XL-mHG
# restriction (k <= L, q >= X) is applied HERE so observed and null pass through
# the identical code path -- the calibration depends on that (amendment A.4).
.mhg_core <- function(is_hit, m, N, boundaries, L = N, X = 1L) {
  pos <- which(is_hit)
  if (length(pos) == 0L) {
    return(list(log_stat = 0, cutoff = NA_integer_, overlap = 0L, le_idx = integer(0)))
  }
  q_all <- cumsum(is_hit)[boundaries]
  # Eligible: adds a new hit, has at least X hits, and sits no deeper than L.
  keep <- q_all >= X & c(TRUE, diff(q_all) > 0L) & boundaries <= L
  if (!any(keep)) {
    # Nothing enriched within the top L (or never X hits) -> not a top hit.
    return(list(log_stat = 0, cutoff = NA_integer_, overlap = 0L, le_idx = integer(0)))
  }
  k <- boundaries[keep]
  q <- q_all[keep]

  log_pv <- stats::phyper(
    q = q - 1L, m = m, n = N - m, k = k,
    lower.tail = FALSE, log.p = TRUE
  )
  log_stat <- min(log_pv)
  j <- max(which(log_pv == log_stat))
  cutoff <- k[j]

  list(
    log_stat = log_stat,
    cutoff = as.integer(cutoff),
    overlap = as.integer(q[j]),
    le_idx = pos[pos <= cutoff]
  )
}

# O(m) tie-free twin of .mhg_core. When every rank is its own tie-block
# (boundaries == seq_len(N)), the eligible cutoffs are exactly the hit positions in
# the top L, so the length-N keep/q_all/is_hit machinery collapses to the sorted hit
# positions `pos`. Byte-identical to .mhg_core on tie-free input; the permutation
# null passes its already-sorted sample here directly, never materializing is_hit.
# `pos` must be the ascending hit positions; `m == length(pos)`. Internal kernel.
.mhg_core_pos <- function(pos, m, N, L = N, X = 1L) {
  none <- list(
    log_stat = 0, cutoff = NA_integer_, overlap = 0L, le_idx = integer(0)
  )
  if (length(pos) == 0L) {
    return(none)
  }
  jL <- sum(pos <= L) # hits within the top L (pos is sorted ascending)
  if (jL < X) {
    return(none) # fewer than X hits in the top L -> not a top hit
  }
  q <- X:jL # candidate overlaps; the cutoffs are the hit positions at those ranks
  k <- pos[q]
  log_pv <- stats::phyper(
    q = q - 1L, m = m, n = N - m, k = k,
    lower.tail = FALSE, log.p = TRUE
  )
  log_stat <- min(log_pv)
  overlap <- q[max(which(log_pv == log_stat))]
  list(
    log_stat = log_stat,
    cutoff = as.integer(pos[overlap]),
    overlap = as.integer(overlap),
    le_idx = pos[seq_len(overlap)]
  )
}
