#' Optional symmetric winsorizing for a signed effect (LFC) vector
#'
#' Caps the tails of a SIGNED vector so a few wild values cannot dominate -- a
#' fallback for cleaning a log-fold-change when no standard error is available for
#' [ohg_shrink_lfc()]. Works in signed space: any value whose magnitude exceeds the
#' `p` quantile of `abs(x)` is clipped to `+cap` or `-cap`
#' (`cap = quantile(abs(x), p)`), both tails symmetrically. It does **not** take
#' `abs()` -- the sign is kept so the result can drive
#' `rank_stat = clean_lfc * -log10(p)`; take `abs()` yourself for the weight. The
#' OHG core never winsorizes; this is a visible, opt-in step in your own script.
#'
#' @param x Numeric signed vector (e.g. log-fold-changes).
#' @param p Upper quantile of `abs(x)` to cap at, a single number in `(0, 1]`. This
#'   is a dataset decision, not a universal constant: raise it when the tail is
#'   genuinely large-effect rather than low-count noise. Default `0.99`.
#'
#' @return `x` with both tails capped, carrying attributes `cap` (the cap value)
#'   and `n_capped` (how many values were clipped) so the transform is inspectable.
#'
#' @examples
#' ohg_winsorize(c(0.1, -0.3, 8, -9, 0.2), p = 0.8)
#'
#' @export
ohg_winsorize <- function(x, p = 0.99) {
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p > 1) {
    stop("`p` must be a single number in (0, 1].", call. = FALSE)
  }
  ax <- abs(x)
  cap <- stats::quantile(ax[is.finite(ax)], probs = p, names = FALSE)
  hi <- !is.na(x) & x > cap
  lo <- !is.na(x) & x < -cap
  out <- x
  out[hi] <- cap
  out[lo] <- -cap
  attr(out, "cap") <- cap
  attr(out, "n_capped") <- sum(hi) + sum(lo)
  out
}
