#' Normalized leading-edge score (NLES) and observed magnitude
#'
#' Standardizes the observed leading-edge magnitude against the permutation null
#' of leading-edge magnitudes for the same set size. The magnitude axis is always
#' `abs(weight)`; the sign of `weight` is never read (direction comes from
#' `dir_sign`). The column -- not the pathway -- is gated to `NA` when the null
#' spread is degenerate or the null is too thin to standardize reliably.
#'
#' @param le_obs Integer positions of the observed leading-edge hits.
#' @param weight Numeric per-gene magnitude aligned to the ranked list, or `NULL`.
#' @param le_idx_b List of per-permutation leading-edge position vectors.
#' @param dir_sign `+1` (up run) or `-1` (down run).
#' @param robust If `TRUE`, use median/MAD; if `FALSE`, mean/SD.
#' @param min_perm_nles Minimum number of permutations required for a reported NLES.
#' @param min_nles_support Minimum number of distinct null magnitudes required.
#'
#' @return A list with `E_obs`, `NLES`, and `NLES_signed`.
#'
#' @examples
#' null <- ohg_permutation_null(N = 100L, m = 8L, B = 1500L)
#' compute_effect(1:6, runif(100), null$le_idx_b, 1, TRUE, 1000L, 10L)
#'
#' @export
compute_effect <- function(le_obs, weight, le_idx_b, dir_sign, robust,
                           min_perm_nles, min_nles_support) {
  if (is.null(weight)) {
    return(list(E_obs = NA_real_, NLES = NA_real_, NLES_signed = NA_real_))
  }

  w <- abs(weight)
  center <- if (robust) stats::median else mean
  scale_fn <- if (robust) stats::mad else stats::sd

  E_obs <- center(w[le_obs])
  E_b <- vapply(le_idx_b, function(idx) center(w[idx]), numeric(1))

  B <- length(E_b)
  spread <- scale_fn(E_b)
  eps <- .Machine$double.eps^0.5
  n_distinct <- length(unique(E_b))

  degenerate_spread <- !is.finite(spread) || spread < eps
  if (B < min_perm_nles || degenerate_spread || n_distinct < min_nles_support) {
    msg <- if (degenerate_spread) {
      paste0(
        "NLES skipped: the permutation null has near-zero spread (mad(E_b) ~= 0). ",
        "This usually means the `weight` vector has a high density of identical or ",
        "exact-zero values, so random leading edges yield near-constant median ",
        "magnitudes."
      )
    } else if (B < min_perm_nles) {
      sprintf(
        paste0(
          "NLES skipped: only %d permutation(s) available (< min_perm_nles = %d), ",
          "too few to size the effect reliably. Increase `n_perm`."
        ),
        B, min_perm_nles
      )
    } else {
      sprintf(
        paste0(
          "NLES skipped: the permutation null has only %d distinct leading-edge ",
          "magnitude(s) (< min_nles_support = %d), too thin to standardize against."
        ),
        n_distinct, min_nles_support
      )
    }
    warning(msg, call. = FALSE)
    return(list(E_obs = E_obs, NLES = NA_real_, NLES_signed = NA_real_))
  }

  nles <- (E_obs - center(E_b)) / spread
  list(E_obs = E_obs, NLES = nles, NLES_signed = dir_sign * nles)
}
