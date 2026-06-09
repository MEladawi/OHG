#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods is
#' @importFrom stats phyper p.adjust p.adjust.methods median mad
## usethis namespace: end
NULL

# Bare column names used in ggplot2::aes() and dplyr verbs (avoids an rlang import).
utils::globalVariables(c("rank", "cum_hits", "p_value", "pathway"))
