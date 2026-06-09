#' Plot the running leading-edge curve for one gene set
#'
#' Shows the cumulative hit count of one gene set along the ranked list and marks
#' the optimal mHG cutoff. Requires the suggested `ggplot2` package.
#'
#' @param ranked_genes Character vector. When `rank_stat` is supplied the order is
#'   irrelevant (the plot sorts by it); with `rank_stat = NULL` it is the ranking.
#' @param gene_set Character vector of pathway genes.
#' @param rank_stat Optional numeric ordering statistic for a tie-aware cutoff.
#' @param max_cutoff_frac,min_hits XL-mHG restriction for the marked cutoff, matching
#'   [ohg_enrichment()] (`L = ceil(max_cutoff_frac * N)`, `X = min_hits`). Default
#'   `max_cutoff_frac = 1` (unrestricted); set `0.25` to mirror the default analysis.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' genes <- paste0("g", 1:50)
#' lfc <- rnorm(50)
#' p <- runif(50)
#' # same recipe as ohg_enrichment(): clean the LFC, rank by clean_lfc * -log10(p)
#' plot_ohg_leading_edge(
#'   genes, genes[c(1, 2, 3, 20)],
#'   rank_stat = ohg_winsorize(lfc) * -log10(p)
#' )
#' }
#'
#' @export
plot_ohg_leading_edge <- function(ranked_genes, gene_set, rank_stat = NULL,
                                  max_cutoff_frac = 1, min_hits = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "`plot_ohg_leading_edge()` needs the 'ggplot2' package (in Suggests).",
      call. = FALSE
    )
  }
  N <- length(ranked_genes)
  if (!is.null(rank_stat) && length(rank_stat) != N) {
    stop("`rank_stat` must be the same length as `ranked_genes`.", call. = FALSE)
  }
  # Order by rank_stat (descending) so the caller need not pre-sort, like ohg_enrichment().
  if (!is.null(rank_stat)) {
    ord <- order(rank_stat, decreasing = TRUE)
    ranked_genes <- ranked_genes[ord]
    rank_stat <- rank_stat[ord]
  }
  L <- as.integer(ceiling(max_cutoff_frac * N))
  t_eff <- intersect(unique(gene_set), ranked_genes)
  boundaries <- tie_boundaries(rank_stat, n = N)
  stat <- ohg_statistic(
    ranked_genes, t_eff, N,
    boundaries = boundaries, L = L, X = as.integer(min_hits)
  )

  df <- tibble::tibble(
    rank = seq_len(N),
    cum_hits = cumsum(ranked_genes %in% t_eff)
  )
  ggplot2::ggplot(df, ggplot2::aes(x = rank, y = cum_hits)) +
    ggplot2::geom_step() +
    ggplot2::geom_vline(
      xintercept = stat$cutoff, linetype = "dashed", colour = "firebrick"
    ) +
    ggplot2::labs(
      x = "Rank position", y = "Cumulative hits",
      title = "OHG leading edge",
      subtitle = sprintf("cutoff = %s, overlap = %d", stat$cutoff, stat$overlap)
    ) +
    ggplot2::theme_minimal()
}
