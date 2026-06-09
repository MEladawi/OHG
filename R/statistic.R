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
#'
#' @return A list with components `log_stat` (natural log of the mHG statistic),
#'   `cutoff` (optimal prefix size, or `NA`), `overlap` (hits inside the leading
#'   edge), and `le_idx` (integer positions of the leading-edge hits).
#'
#' @examples
#' ohg_statistic(paste0("g", 1:10), c("g1", "g2", "g5"), N = 10L)
#'
#' @export
ohg_statistic <- function(ranked_genes, T_eff, N, boundaries = NULL) {
  m <- length(T_eff)
  is_hit <- ranked_genes %in% T_eff
  if (is.null(boundaries)) boundaries <- seq_len(N) # distinct ranks
  .mhg_core(is_hit, m, N, boundaries)
}

# Internal mHG kernel shared by ohg_statistic() and ohg_permutation_null().
# Given a length-N hit indicator and the candidate cutoff `boundaries`, evaluates
# the upper-tail hypergeometric only at boundaries that add a new hit (constant q
# with larger k can only worsen the tail; statistic plan section 2.3) and returns
# the most-enriched prefix in natural-log space. The single source of truth for
# the statistic; both the observed and permutation paths call it.
.mhg_core <- function(is_hit, m, N, boundaries) {
  pos <- which(is_hit)
  if (length(pos) == 0L) {
    return(list(log_stat = 0, cutoff = NA_integer_, overlap = 0L, le_idx = integer(0)))
  }
  q_all <- cumsum(is_hit)[boundaries]
  keep <- q_all >= 1L & c(TRUE, diff(q_all) > 0L)
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
