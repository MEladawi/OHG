#' Coerce a gene-set collection to a named list of character vectors
#'
#' Normalizes the three accepted `gene_sets` inputs at the package boundary: a
#' named `list` of character vectors, a path to a `.gmt` file, or a
#' `GSEABase::GeneSetCollection`. Each set is de-duplicated. The `.gmt` path is
#' parsed by [read_gmt()]; the `GeneSetCollection` branch needs the suggested
#' GSEABase package.
#'
#' @param gene_sets A named list of character vectors, a single `.gmt` file path,
#'   or a `GSEABase::GeneSetCollection`.
#'
#' @return A named list of unique character vectors.
#'
#' @examples
#' coerce_gene_sets(list(SET_A = c("g1", "g2"), SET_B = c("g2", "g3")))
#'
#' @export
coerce_gene_sets <- function(gene_sets) {
  if (is.character(gene_sets) && length(gene_sets) == 1L && file.exists(gene_sets)) {
    return(coerce_gene_sets(read_gmt(gene_sets)))
  }
  if (methods::is(gene_sets, "GeneSetCollection")) {
    if (!requireNamespace("GSEABase", quietly = TRUE)) {
      stop(
        "A `GeneSetCollection` was supplied but the 'GSEABase' package is not ",
        "installed. Install GSEABase or pass a named list / .gmt path instead.",
        call. = FALSE
      )
    }
    ids <- GSEABase::geneIds(gene_sets)
    return(lapply(ids, function(g) unique(as.character(g))))
  }
  if (is.list(gene_sets)) {
    nm <- names(gene_sets)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("`gene_sets` list must be fully named (one name per gene set).", call. = FALSE)
    }
    return(lapply(gene_sets, function(g) unique(as.character(g))))
  }
  stop(
    "`gene_sets` must be a named list, a `.gmt` file path, or a ",
    "`GSEABase::GeneSetCollection`.",
    call. = FALSE
  )
}
