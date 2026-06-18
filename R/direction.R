#' Tie-block boundary positions
#'
#' Given the ordering statistic, returns the last position of each maximal run of
#' equal values. Cutoffs in the mHG scan are restricted to these boundaries so a
#' tie block is included whole or not at all (within a tie block, gene order is
#' arbitrary). With distinct ranks this reduces to `seq_len(n)`.
#'
#' @param rank_stat Numeric ordering statistic (already sorted non-increasing),
#'   or `NULL` for a fully-resolved order.
#' @param n Length of the ranked list; used when `rank_stat` is `NULL`.
#'
#' @return Integer vector of boundary positions, strictly increasing, ending at `n`.
#'
#' @examples
#' tie_boundaries(c(5, 5, 4, 3, 3, 3, 1))
#'
#' @export
tie_boundaries <- function(rank_stat, n = length(rank_stat)) {
  if (is.null(rank_stat)) {
    # length(NULL) == 0, so the default n would silently yield integer(0). Require
    # an explicit n: the caller alone knows the list length when no rank_stat is
    # given to measure it from.
    if (missing(n)) {
      stop("`n` is required when `rank_stat` is NULL.", call. = FALSE)
    }
    return(seq_len(n))
  }
  n <- length(rank_stat)
  if (n == 0L) {
    return(integer(0))
  }
  changed <- c(rank_stat[-n] != rank_stat[-1L], TRUE)
  which(changed)
}

#' Infer the default test direction from `rank_stat`
#'
#' A signed `rank_stat` (one that crosses zero) defaults to `"both"` -- both ends
#' are meaningful. A non-negative or absent `rank_stat` defaults to `"up"` -- only
#' the top is meaningful. A user-supplied `direction` always wins.
#'
#' @param rank_stat Numeric ordering statistic, or `NULL`.
#' @param supplied A user-supplied direction (`"up"`/`"down"`/`"both"`) or `NULL`.
#'
#' @return One of `"up"`, `"down"`, `"both"`.
#'
#' @examples
#' infer_direction(c(3, 1, -1, -4))
#' infer_direction(c(5, 4, 3, 1))
#'
#' @export
infer_direction <- function(rank_stat, supplied = NULL) {
  if (!is.null(supplied)) {
    return(match.arg(supplied, c("up", "down", "both")))
  }
  if (is.null(rank_stat)) {
    return("up")
  }
  if (any(rank_stat < 0) && any(rank_stat > 0)) "both" else "up"
}
