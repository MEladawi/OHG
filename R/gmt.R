#' Read a GMT gene-set file
#'
#' Parses an MSigDB / Enrichr / g:Profiler `.gmt` file in which each line is
#' `set_name <tab> description <tab> gene1 <tab> gene2 ...`. The description
#' field is read and dropped from the values, but retained on the
#' `"descriptions"` attribute of the returned list. Base R only.
#'
#' @param path Path to a `.gmt` file.
#'
#' @return A named list of unique character vectors (genes per set), with a named
#'   character vector of descriptions on `attr(x, "descriptions")`.
#'
#' @examples
#' gmt <- tempfile(fileext = ".gmt")
#' writeLines("SET_A\tdescription\tg1\tg2\tg3", gmt)
#' read_gmt(gmt)
#'
#' @export
read_gmt <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("`path` must be a single existing file path to a .gmt file.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0L) {
    stop("GMT file is empty: ", path, call. = FALSE)
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  set_names <- vapply(fields, function(f) trimws(f[[1L]]), character(1))
  descriptions <- vapply(
    fields,
    function(f) if (length(f) >= 2L) trimws(f[[2L]]) else "",
    character(1)
  )
  genes <- lapply(fields, function(f) {
    g <- if (length(f) >= 3L) f[-(1:2)] else character(0)
    g <- trimws(g)
    unique(g[nzchar(g)])
  })

  names(genes) <- set_names
  names(descriptions) <- set_names
  attr(genes, "descriptions") <- descriptions
  genes
}
