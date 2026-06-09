#' Optional, tool-agnostic LFC shrinkage
#'
#' Shrinks log-fold-changes from `(LFC, SE)` under a unimodal prior via `ashr`.
#' Shrinkage is a function of `(LFC, SE)`, not `(LFC, p-value)`, so it cannot be
#' recovered from a table carrying only LFC and p. The OHG core never shrinks and
#' never reverse-engineers an SE from a p-value; this is a user-facing convenience
#' only, agnostic to the DE tool that produced the estimates.
#'
#' @param lfc Numeric log-fold-changes.
#' @param se Numeric standard errors, the same length as `lfc`.
#' @param method Currently only `"ashr"`.
#'
#' @return Numeric vector of posterior-mean shrunken LFCs.
#'
#' @examples
#' \dontrun{
#' ohg_shrink_lfc(c(3, 0.1, -2), c(0.2, 2.0, 0.3))
#' }
#'
#' @export
ohg_shrink_lfc <- function(lfc, se, method = "ashr") {
  method <- match.arg(method, "ashr")
  if (length(lfc) != length(se)) {
    stop("`lfc` and `se` must be the same length.", call. = FALSE)
  }
  if (any(!is.finite(lfc)) || any(!is.finite(se))) {
    stop("`lfc` and `se` must be finite (no NA/NaN/Inf).", call. = FALSE)
  }
  if (any(se < 0)) {
    stop("`se` must be non-negative (it is a standard error).", call. = FALSE)
  }
  if (!requireNamespace("ashr", quietly = TRUE)) {
    stop(
      "`ohg_shrink_lfc()` needs the 'ashr' package (in Suggests). ",
      "Install it, or pass an already-shrunken LFC.",
      call. = FALSE
    )
  }
  fit <- ashr::ash(betahat = lfc, sebetahat = se, mixcompdist = "normal")
  ashr::get_pm(fit)
}
