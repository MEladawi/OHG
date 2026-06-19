#' Validate and normalize OHG inputs
#'
#' Enforces the input contract: a non-empty character ranking, an aligned finite
#' `rank_stat` (if present), and a `weight` that is finite (hard error) but
#' otherwise the user's responsibility (soft warnings on pathological-but-finite
#' values). When `rank_stat` is given, **OHG sorts the genes by it (descending)**
#' so callers need not pre-sort; ties keep their given order. De-duplicates
#' `ranked_genes`, keeping the best-ranked copy.
#'
#' @param ranked_genes Character vector. When `rank_stat` is supplied the order is
#'   irrelevant (OHG sorts); with `rank_stat = NULL` it is taken as the ranking.
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

  n0 <- length(ranked_genes)
  if (!is.null(rank_stat) &&
    (!is.numeric(rank_stat) || length(rank_stat) != n0 || any(!is.finite(rank_stat)))) {
    stop("`rank_stat` must be finite numeric aligned to `ranked_genes`.", call. = FALSE)
  }
  if (!is.null(weight) && (!is.numeric(weight) || length(weight) != n0)) {
    stop("`weight` must be numeric aligned to `ranked_genes`.", call. = FALSE)
  }

  # OHG orders the list itself -- callers pass genes and rank_stat in any order.
  # Sort by rank_stat (descending; stable, so tie blocks keep their given order).
  if (!is.null(rank_stat)) {
    ord <- order(rank_stat, decreasing = TRUE)
    ranked_genes <- ranked_genes[ord]
    rank_stat <- rank_stat[ord]
    if (!is.null(weight)) weight <- weight[ord]
  }

  # De-duplicate, keeping the first (now best-ranked) occurrence of each gene.
  dup <- duplicated(ranked_genes)
  if (any(dup)) {
    warning(
      sum(dup), " duplicate gene(s) dropped, keeping the best-ranked copy.",
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
    tie_fraction <- mean(duplicated(rank_stat))
    if (tie_fraction > 0) {
      message(sprintf(
        paste0(
          "rank_stat has ties (%.1f%% of genes share a value); using tie-aware ",
          "evaluation. Rank by a finer statistic to shrink tie blocks."
        ),
        100 * tie_fraction
      ))
    }
  }

  if (!is.null(weight)) {
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
#' drops sets with `m < min_set_size` (`m` cannot exceed `N`: `T_eff` is an
#' intersection with the universe). `min_set_size` controls significance inclusion
#' only; it never gates the NLES column.
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
  # m = |T_eff| is bounded above by N because T_eff is an intersection with the
  # universe, so only the lower size filter can ever fire.
  for (nm in names(gene_sets)) {
    t_eff <- intersect(unique(gene_sets[[nm]]), universe)
    m <- length(t_eff)
    if (m < min_set_size) {
      dropped[[length(dropped) + 1L]] <- tibble::tibble(
        pathway = nm,
        reason = sprintf("m=%d < min_set_size=%d", m, min_set_size)
      )
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
