#' Ordered hypergeometric (mHG) enrichment
#'
#' Runs the ordered / minimum-hypergeometric enrichment test for one ranked gene
#' list against a collection of gene sets. The reported `p_value` is a
#' permutation-calibrated p-value (significance); localization is reported by
#' `cutoff_rank` and `leading_edge_fraction`; magnitude by `NLES`. These three
#' axes are kept separate and never pre-mixed — the package emits no composite
#' ranking column.
#'
#' @param ranked_genes Unique character vector, most important first. The universe
#'   is the ranked list; there is no separate universe-size argument.
#' @param gene_sets A named list of character vectors, a `.gmt` file path, or a
#'   `GSEABase::GeneSetCollection`.
#' @param rank_stat Numeric ordering statistic (ordering only — never multiplied
#'   into the test); `NULL` assumes a fully-resolved order.
#' @param weight Numeric effect magnitude (`abs()` is used; sign ignored); `NULL`
#'   yields `NLES = NA` with overlap-based outputs still reported.
#' @param direction `NULL` (inferred from the `rank_stat` sign), or one of
#'   `"up"`, `"down"`, `"both"`.
#' @param p_adjust_method Any `stats::p.adjust.methods`; default `"BH"` (FDR).
#' @param n_perm Permutations per distinct set size.
#' @param min_set_size Significance-inclusion floor; does not gate `NLES`.
#' @param min_perm_nles,min_nles_support `NLES` stability gates.
#' @param robust_nles Median/MAD (default) versus mean/SD.
#' @param collapse_both If `TRUE` with `direction = "both"`, keep the more
#'   significant direction per pathway with a x2 Bonferroni penalty before adjustment.
#' @param method `"permutation"` (the only implemented method).
#' @param seed Integer RNG seed or `NULL`.
#' @param n_cores Worker count; `> 1` parallelizes the per-size nulls via furrr.
#'
#' @return A tibble, one row per `(pathway[, direction])`, sorted by `p_value`.
#'
#' @examples
#' ranked <- paste0("g", 1:200)
#' sets <- list(HIT = ranked[1:12], MEH = sample(ranked, 15))
#' ohg_enrichment(ranked, sets, n_perm = 500L, seed = 1)
#'
#' @export
ohg_enrichment <- function(ranked_genes, gene_sets, rank_stat = NULL, weight = NULL,
                           direction = NULL, p_adjust_method = "BH", n_perm = 2000L,
                           min_set_size = 3L, min_perm_nles = 1000L,
                           min_nles_support = 10L, robust_nles = TRUE,
                           collapse_both = FALSE, method = "permutation",
                           seed = NULL, n_cores = 1L) {
  method <- match.arg(method, c("permutation", "exact"))
  if (method == "exact") {
    stop("method = 'exact' is not implemented yet.", call. = FALSE)
  }

  v <- validate_inputs(ranked_genes, rank_stat, weight, p_adjust_method)
  sets <- coerce_gene_sets(gene_sets)
  dir <- infer_direction(v$rank_stat, supplied = direction)
  dirs <- if (dir == "both") c("up", "down") else dir

  pruned <- prune_gene_sets(sets, v$ranked_genes, min_set_size, v$N)
  if (length(pruned$kept) == 0L) {
    stop("No gene set passed min_set_size = ", min_set_size, ".", call. = FALSE)
  }

  if (!is.null(seed)) {
    # Use L'Ecuyer-CMRG for reproducible per-size streams, but leave the caller's
    # global RNG kind and seed exactly as we found them once we return.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    on.exit(
      {
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
  }

  rows <- list()
  for (run_dir in dirs) {
    dir_sign <- if (run_dir == "up") 1 else -1
    if (run_dir == "up") {
      rg <- v$ranked_genes
      rs <- v$rank_stat
      w <- v$weight
    } else {
      rg <- rev(v$ranked_genes)
      rs <- if (is.null(v$rank_stat)) NULL else rev(v$rank_stat)
      w <- if (is.null(v$weight)) NULL else rev(v$weight)
    }
    boundaries <- tie_boundaries(rs, n = v$N)

    ms <- vapply(pruned$kept, `[[`, integer(1), "m")
    nulls <- build_nulls(ms, v$N, n_perm, boundaries, n_cores = n_cores, seed = seed)

    rows <- c(rows, purrr::imap(pruned$kept, function(set, nm) {
      stat <- ohg_statistic(rg, set$genes, v$N, boundaries = boundaries)
      if (stat$overlap == 0L) {
        return(NULL)
      }
      nb <- nulls[[as.character(set$m)]]
      p_emp <- (1 + sum(nb$log_stat_b <= stat$log_stat)) / (1 + n_perm)
      eff <- compute_effect(
        stat$le_idx, w, nb$le_idx_b, dir_sign, robust_nles,
        min_perm_nles, min_nles_support
      )
      tibble::tibble(
        pathway = nm, direction = run_dir, set_size = set$m,
        cutoff_rank = stat$cutoff, leading_edge_size = stat$cutoff,
        overlap = stat$overlap, leading_edge_fraction = stat$overlap / set$m,
        neg_log10_mHG = -stat$log_stat / log(10), mHG_stat = exp(stat$log_stat),
        p_value = p_emp, E_obs = eff$E_obs, NLES = eff$NLES,
        NLES_signed = eff$NLES_signed, n_leading_edge = stat$overlap,
        hits = paste(rg[stat$le_idx], collapse = " ")
      )
    }))
  }

  out <- purrr::list_rbind(rows)
  if (nrow(out) == 0L) {
    stop("No pathway had overlap >= 1.", call. = FALSE)
  }

  if (isTRUE(collapse_both) && dir == "both") {
    out <- out |>
      dplyr::mutate(p_value = pmin(1, p_value * 2)) |>
      dplyr::slice_min(p_value, n = 1L, by = pathway, with_ties = FALSE)
  }

  out$p_adjust <- stats::p.adjust(out$p_value, method = p_adjust_method)
  out$p_adjust_method <- p_adjust_method
  out$NES_OHG <- out$NLES

  col_order <- c(
    "pathway", "direction", "set_size", "cutoff_rank", "leading_edge_size",
    "overlap", "leading_edge_fraction", "neg_log10_mHG", "mHG_stat",
    "p_value", "p_adjust", "p_adjust_method", "E_obs", "NLES",
    "NLES_signed", "NES_OHG", "n_leading_edge", "hits"
  )
  if (dir == "up") {
    out$direction <- NULL
    col_order <- setdiff(col_order, "direction")
  }
  out <- out[, col_order]
  out[order(out$p_value), ]
}
