#' Validate and normalize OHG inputs
#'
#' Enforces the input contract: a non-empty unique character ranking, an aligned
#' finite non-increasing `rank_stat` (if present), and a `weight` that is finite
#' (hard error) but otherwise the user's responsibility (soft warnings on
#' pathological-but-finite values). De-duplicates `ranked_genes`, keeping the best
#' (earliest) rank.
#'
#' @param ranked_genes Character vector, most important first.
#' @param rank_stat Numeric ordering statistic aligned to `ranked_genes`, or `NULL`.
#' @param weight Numeric magnitude aligned to `ranked_genes`, or `NULL`.
#' @param p_adjust_method One of `stats::p.adjust.methods`.
#'
#' @return A list with the de-duplicated `ranked_genes`, aligned `rank_stat` and
#'   `weight`, `N`, `tie_fraction`, and the validated `p_adjust_method`.
#'
#' @examples
#' validate_inputs(c("g1", "g2", "g3"), rank_stat = c(3, 2, 1), weight = NULL)
#'
#' @export
validate_inputs <- function(ranked_genes, rank_stat, weight, p_adjust_method = "BH") {
  if (!is.character(ranked_genes) || length(ranked_genes) == 0L) {
    stop("`ranked_genes` must be a non-empty character vector.", call. = FALSE)
  }
  if (!p_adjust_method %in% stats::p.adjust.methods) {
    stop(
      "`p_adjust_method` must be one of: ",
      paste(stats::p.adjust.methods, collapse = ", "), ".",
      call. = FALSE
    )
  }

  dup <- duplicated(ranked_genes)
  if (any(dup)) {
    warning(
      sum(dup), " duplicate gene(s) dropped, keeping best (earliest) rank.",
      call. = FALSE
    )
    keep <- !dup
    ranked_genes <- ranked_genes[keep]
    if (!is.null(rank_stat)) rank_stat <- rank_stat[keep]
    if (!is.null(weight)) weight <- weight[keep]
  }
  N <- length(ranked_genes)

  tie_fraction <- 0
  if (!is.null(rank_stat)) {
    if (!is.numeric(rank_stat) || length(rank_stat) != N || any(!is.finite(rank_stat))) {
      stop("`rank_stat` must be finite numeric aligned to `ranked_genes`.", call. = FALSE)
    }
    if (is.unsorted(rev(rank_stat))) {
      stop(
        "`rank_stat` must be sorted non-increasing (most important first).",
        call. = FALSE
      )
    }
    tie_fraction <- mean(duplicated(rank_stat))
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight) || length(weight) != N) {
      stop("`weight` must be numeric aligned to `ranked_genes`.", call. = FALSE)
    }
    if (any(!is.finite(weight))) {
      stop(
        "`weight` has non-finite values (NA/NaN/Inf). Build any -log10(p) term ",
        "in log space from the finite log-p, never as -log10(0).",
        call. = FALSE
      )
    }
    aw <- abs(weight)
    pos <- aw[aw > 0]
    med <- if (length(pos)) stats::median(pos) else 0
    if (is.finite(med) && med > 0 && max(aw) / med > 1e6) {
      warning(
        "`weight` has extreme dynamic range (max/median > 1e6); consider ",
        "winsorizing. Proceeding with values as given.",
        call. = FALSE
      )
    }
    if (mean(aw == 0) > 0.5) {
      warning("`weight` is >50% exact zeros; NLES may be NA.", call. = FALSE)
    }
    if (length(unique(weight)) == 1L) {
      warning("`weight` is all-equal; NLES will be NA.", call. = FALSE)
    }
  }

  list(
    ranked_genes = ranked_genes, rank_stat = rank_stat, weight = weight,
    N = N, tie_fraction = tie_fraction, p_adjust_method = p_adjust_method
  )
}

#' Intersect, de-duplicate, and size-filter gene sets against the universe
#'
#' Computes `T_eff = unique(T) intersect U` and `m = |T_eff|` for each set, then
#' drops sets with `m < min_set_size` or `m > N`. `min_set_size` controls
#' significance inclusion only; it never gates the NLES column.
#'
#' @param gene_sets Named list of character vectors (post-coercion).
#' @param universe The ranked gene vector.
#' @param min_set_size Minimum `m` for significance inclusion.
#' @param N Universe size.
#'
#' @return A list with `kept` (named list of `list(genes, m)`) and `dropped`
#'   (a tibble of `pathway` and `reason`).
#'
#' @keywords internal
#' @export
prune_gene_sets <- function(gene_sets, universe, min_set_size, N) {
  kept <- list()
  dropped <- list()
  for (nm in names(gene_sets)) {
    t_eff <- intersect(unique(gene_sets[[nm]]), universe)
    m <- length(t_eff)
    if (m < min_set_size) {
      dropped[[length(dropped) + 1L]] <- tibble::tibble(
        pathway = nm,
        reason = sprintf("m=%d < min_set_size=%d", m, min_set_size)
      )
    } else if (m > N) {
      dropped[[length(dropped) + 1L]] <- tibble::tibble(pathway = nm, reason = "m > N")
    } else {
      kept[[nm]] <- list(genes = t_eff, m = m)
    }
  }
  list(
    kept = kept,
    dropped = if (length(dropped)) {
      dplyr::bind_rows(dropped)
    } else {
      tibble::tibble(pathway = character(0), reason = character(0))
    }
  )
}
