#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods is
#' @importFrom stats phyper p.adjust p.adjust.methods median mad
## usethis namespace: end
NULL

# Bare aes() column names in plot_ohg_leading_edge() (avoids an rlang import).
utils::globalVariables(c("rank", "cum_hits"))
