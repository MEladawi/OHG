#' Plot the running leading-edge curve for one gene set
#'
#' Shows the cumulative hit count of one gene set along the ranked list and marks
#' the optimal mHG cutoff. Requires the suggested `ggplot2` package.
#'
#' @param ranked_genes Character vector, most important first.
#' @param gene_set Character vector of pathway genes.
#' @param rank_stat Optional numeric ordering statistic for a tie-aware cutoff.
#' @param ... Reserved for future use.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_ohg_leading_edge(paste0("g", 1:50), paste0("g", c(1, 2, 3, 20)))
#' }
#'
#' @export
plot_ohg_leading_edge <- function(ranked_genes, gene_set, rank_stat = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "`plot_ohg_leading_edge()` needs the 'ggplot2' package (in Suggests).",
      call. = FALSE
    )
  }
  N <- length(ranked_genes)
  t_eff <- intersect(unique(gene_set), ranked_genes)
  boundaries <- tie_boundaries(rank_stat, n = N)
  stat <- ohg_statistic(ranked_genes, t_eff, N, boundaries = boundaries)

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
