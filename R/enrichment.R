#' Ordered hypergeometric (mHG) enrichment
#'
#' Runs the ordered / minimum-hypergeometric enrichment test for one ranked gene
#' list against a collection of gene sets. The reported `p_value` is a
#' permutation-calibrated p-value (significance); localization is reported by
#' `cutoff_rank` and `leading_edge_fraction`; magnitude by `NLES`. These three
#' axes are kept separate and never pre-mixed -- the package emits no composite
#' ranking column.
#'
#' @param ranked_genes Character vector of unique gene ids; the universe is this
#'   list (no separate universe-size argument). When `rank_stat` is supplied OHG
#'   sorts the genes by it, so the order you pass does not matter; with
#'   `rank_stat = NULL` the given order is taken as the ranking (best first).
#' @param gene_sets A named list of character vectors, a `.gmt` file path, or a
#'   `GSEABase::GeneSetCollection`.
#' @param rank_stat Numeric ordering statistic (ordering only -- never multiplied
#'   into the test); `NULL` assumes a fully-resolved order. A signed
#'   `log2FC * -log10(p)` is recommended (continuous, so few ties).
#' @param weight Numeric effect magnitude (`abs()` is used; sign ignored); `NULL`
#'   yields `NLES = NA` with overlap-based outputs still reported. `abs(log2FC)` is
#'   recommended, ideally on a shrunken fold change (see [ohg_shrink_lfc()]).
#' @param direction `NULL` (inferred from the `rank_stat` sign), or one of
#'   `"up"`, `"down"`, `"both"`.
#' @param p_adjust_method Any `stats::p.adjust.methods`; default `"BH"` (FDR).
#' @param n_perm Baseline permutations per distinct set size (`B0`). Also the
#'   fixed permutation count when `method = "permutation"`.
#' @param target_hits Exceedances (`h`) a pathway must accumulate before its
#'   sequential p-value is resolved as `h / L`. Only used when
#'   `method = "adaptive"`.
#' @param n_perm_max Permutation cap (`B_max`). `NULL` (default) uses
#'   `max(1e5, ceil(n_tests / alpha))` so a rank-1 pathway can still clear `alpha`
#'   after correction. Only used when `method = "adaptive"`.
#' @param max_cutoff_frac Largest cutoff the scan may reach, as a fraction of the
#'   list: `L = ceil(max_cutoff_frac * N)`. Keeps the optimal cutoff (and the
#'   leading edge) near the top instead of placing it deep in the list on a
#'   diffuse whole-list shift. Default `0.25` (top quarter); lower it (e.g. `0.10`)
#'   to demand sharper signals, set `1` for the unrestricted scan. The same `L` is
#'   applied to the permutation null, so the p-value stays calibrated. The resolved
#'   `L` is attached to the result as `attr(x, "L_used")`.
#' @param min_hits Minimum leading-edge genes before a cutoff is eligible (XL-mHG
#'   `X`). Default `1` (off); `3`-`5` ignores "enrichment" resting on one or two top
#'   genes. Attached as `attr(x, "X_used")`.
#' @param alpha Significance level used only to set the default `n_perm_max` and
#'   to decide whether to warn about cap-limited pathways.
#' @param min_set_size Significance-inclusion floor; does not gate `NLES`.
#' @param min_perm_nles,min_nles_support `NLES` stability gates.
#' @param robust_nles Median/MAD (default) versus mean/SD.
#' @param collapse_both If `TRUE` with `direction = "both"`, keep the more
#'   significant direction per pathway with a x2 Bonferroni penalty before adjustment.
#' @param method One of `"adaptive"` (default; Besag-Clifford sequential
#'   Monte-Carlo, removing the `1/(n_perm + 1)` floor by spending extra
#'   permutations only on near-floor pathways), `"permutation"` (fixed-`n_perm`
#'   empirical p-value), or `"exact"` (the mHG dynamic program, not yet
#'   implemented).
#' @param seed Integer RNG seed or `NULL`.
#' @param n_cores Worker count. With `> 1` (and `future`/`furrr` installed) the
#'   per-size work -- null draws and pathway scoring -- runs in parallel. The
#'   package sets up the `future` plan itself (multicore where supported, else
#'   multisession) and restores the caller's plan on exit; just pass `n_cores`.
#'
#' @return A tibble, one row per `(pathway[, direction])`, sorted by `p_value`.
#'   Resampling diagnostics `n_perm_used` (draws spent on this size), `n_exceed`
#'   (exceedance count), and `resolution_limited` (`TRUE` when the cap was hit
#'   before resolving, so `p_value` is a conservative bound) are included.
#'
#' @examples
#' # Minimal: a fully-ordered list, no rank_stat/weight.
#' ranked <- paste0("g", 1:200)
#' sets <- list(HIT = ranked[1:12], MEH = sample(ranked, 15))
#' ohg_enrichment(ranked, sets, n_perm = 500L, seed = 1)
#'
#' # Recommended recipe from a DE table: clean the log-fold-change once, then
#' # rank by clean_lfc * -log10(p) and weight by abs(clean_lfc). OHG sorts by
#' # rank_stat for you -- pass the genes in any order.
#' set.seed(1)
#' lfc <- rnorm(200, sd = 2) # toy log-fold-changes
#' p <- runif(200) # toy p-values
#' clean_lfc <- ohg_winsorize(lfc) # no SE here -> winsorize; else ohg_shrink_lfc()
#' ohg_enrichment(
#'   ranked, sets,
#'   rank_stat = clean_lfc * -log10(p),
#'   weight = abs(clean_lfc),
#'   n_perm = 500L, seed = 1
#' )
#'
#' @export
ohg_enrichment <- function(ranked_genes, gene_sets, rank_stat = NULL, weight = NULL,
                           direction = NULL, p_adjust_method = "BH", n_perm = 2000L,
                           target_hits = 10L, n_perm_max = NULL,
                           max_cutoff_frac = 0.25, min_hits = 1L,
                           alpha = 0.05, min_set_size = 3L, min_perm_nles = 1000L,
                           min_nles_support = 10L, robust_nles = TRUE,
                           collapse_both = FALSE, method = "adaptive",
                           seed = NULL, n_cores = 1L) {
  method <- match.arg(method, c("adaptive", "permutation", "exact"))
  if (method == "exact") {
    stop("method = 'exact' is not implemented yet.", call. = FALSE)
  }
  adaptive <- method == "adaptive"

  v <- validate_inputs(ranked_genes, rank_stat, weight, p_adjust_method)
  sets <- coerce_gene_sets(gene_sets)
  dir <- infer_direction(v$rank_stat, supplied = direction)
  dirs <- if (dir == "both") c("up", "down") else dir

  pruned <- prune_gene_sets(sets, v$ranked_genes, min_set_size, v$N)
  if (length(pruned$kept) == 0L) {
    stop("No gene set passed min_set_size = ", min_set_size, ".", call. = FALSE)
  }

  # XL-mHG restriction: the optimal cutoff may not sit deeper than L = the top
  # `max_cutoff_frac` of the list (so a "leading edge" stays near the top, not a
  # diffuse whole-list shift), and needs at least X = min_hits genes. Applied
  # identically to the observed statistic and the permutation null (statistic.R).
  if (length(max_cutoff_frac) != 1L || !is.finite(max_cutoff_frac) ||
    max_cutoff_frac <= 0 || max_cutoff_frac > 1) {
    stop("`max_cutoff_frac` must be a single number in (0, 1].", call. = FALSE)
  }
  if (length(min_hits) != 1L || !is.finite(min_hits) || min_hits < 1) {
    stop("`min_hits` must be a single integer >= 1.", call. = FALSE)
  }
  L <- as.integer(ceiling(max_cutoff_frac * v$N))
  X <- as.integer(min_hits)

  # Adaptive (Besag-Clifford) resampling removes the fixed-B floor 1/(n_perm+1).
  # Baseline draws B0 = max(n_perm, min_perm_nles) feed both the p-value and NLES;
  # the cap B_max defaults so a rank-1 pathway can still clear alpha after BH
  # (floor 1/(B_max+1) <= alpha / n_tests). n_perm0 == n_perm when adaptive = FALSE.
  n_tests <- length(pruned$kept) * length(dirs)
  n_perm0 <- if (adaptive) max(n_perm, min_perm_nles) else n_perm
  b_max <- if (is.null(n_perm_max)) {
    max(100000L, as.integer(ceiling(n_tests / alpha)))
  } else {
    as.integer(n_perm_max)
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

  # Parallel backend: the package owns plan setup so callers only pass n_cores.
  # Prefer multicore (fork) where supported -- it shares memory copy-on-write, so
  # the ranked-list/gene-set globals are never serialized to workers and forked
  # workers inherit the loaded namespace; otherwise fall back to multisession
  # (portable). The caller's own plan is restored on exit.
  use_par <- n_cores > 1L &&
    requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)
  if (n_cores > 1L && !use_par) {
    warning(
      "n_cores > 1 needs the 'future' and 'furrr' packages; running sequentially.",
      call. = FALSE
    )
  }
  if (use_par) {
    backend <- if (future::supportsMulticore()) {
      future::multicore
    } else {
      future::multisession
    }
    oplan <- future::plan(backend, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
  }

  # Per-direction context: the (possibly reversed) ranked list, weights, tie-block
  # boundaries and direction sign. Small relative to the permutation nulls.
  dir_ctx <- lapply(dirs, function(run_dir) {
    if (run_dir == "up") {
      rg <- v$ranked_genes
      rs <- v$rank_stat
      w <- v$weight
    } else {
      rg <- rev(v$ranked_genes)
      rs <- if (is.null(v$rank_stat)) NULL else rev(v$rank_stat)
      w <- if (is.null(v$weight)) NULL else rev(v$weight)
    }
    list(
      label = run_dir, dir_sign = if (run_dir == "up") 1 else -1,
      rg = rg, w = w, boundaries = tie_boundaries(rs, n = v$N)
    )
  })

  # One parallel unit per distinct set size: a worker draws that size's null and
  # scores every pathway (and direction) of that size, so the heavy work (null
  # draws + NLES standardization) runs on the worker and only small result rows
  # return. RNG is keyed to m (set.seed(seed + m)), so results are identical
  # sequentially or in parallel and across worker counts. With adaptive = TRUE the
  # sequential null is decided per (size, orientation) -- re-seeded per orientation
  # so a `both` run's draws match the separate up/down runs exactly. With
  # adaptive = FALSE the fixed-B null is reused across tails when boundaries match.
  by_size <- split(pruned$kept, vapply(pruned$kept, `[[`, integer(1), "m"))
  sizes <- as.integer(names(by_size))

  score_size <- function(sets_m, m) {
    fixed_cache <- list() # fixed-B path: null draws keyed by boundaries
    rows_m <- list()
    for (ctx in dir_ctx) {
      stats <- lapply(sets_m, function(set) {
        ohg_statistic(ctx$rg, set$genes, v$N, boundaries = ctx$boundaries, L = L, X = X)
      })
      # unname so imap() below indexes p_info by position, not by pathway name
      keep <- unname(which(vapply(stats, function(s) s$overlap >= 1L, logical(1))))
      if (length(keep) == 0L) next
      obs <- vapply(stats[keep], `[[`, numeric(1), "log_stat")

      if (adaptive) {
        an <- .adaptive_null(
          obs, m, v$N, ctx$boundaries, n_perm0, b_max, target_hits, seed,
          L = L, X = X
        )
        le_idx_b <- an$le_idx_b
        p_info <- list(
          p = an$p_hat, c = an$n_exceed,
          used = an$n_perm_used, rl = an$resolution_limited
        )
      } else {
        nb <- NULL
        for (nc in fixed_cache) {
          if (identical(nc$boundaries, ctx$boundaries)) {
            nb <- nc$nb
            break
          }
        }
        if (is.null(nb)) {
          if (!is.null(seed)) set.seed(seed + m) # stream keyed to m, not order
          nb <- ohg_permutation_null(v$N, m, n_perm, ctx$boundaries, L = L, X = X)
          fixed_cache <- c(
            fixed_cache, list(list(boundaries = ctx$boundaries, nb = nb))
          )
        }
        le_idx_b <- nb$le_idx_b
        c_exc <- vapply(obs, function(o) sum(nb$log_stat_b <= o), integer(1))
        p_info <- list(
          p = (1 + c_exc) / (1 + n_perm), c = c_exc,
          used = n_perm, rl = rep(FALSE, length(keep))
        )
      }

      # NLES summary uses this direction's weights -> direction-specific; hoisted
      # once per (size, direction) from the (possibly expanded) null buffer.
      summ <- if (is.null(v$weight)) {
        NULL
      } else {
        .nles_null_summary(
          le_idx_b, ctx$w, robust_nles, min_perm_nles, min_nles_support
        )
      }

      rows_m <- c(rows_m, purrr::imap(keep, function(j, i) {
        st <- stats[[j]]
        set <- sets_m[[j]]
        eff <- .nles_from_summary(st$le_idx, ctx$w, ctx$dir_sign, summ)
        tibble::tibble(
          pathway = names(sets_m)[j], direction = ctx$label, set_size = set$m,
          cutoff_rank = st$cutoff, leading_edge_size = st$cutoff,
          overlap = st$overlap, leading_edge_fraction = st$overlap / set$m,
          neg_log10_mHG = -st$log_stat / log(10), mHG_stat = exp(st$log_stat),
          p_value = p_info$p[i], n_perm_used = p_info$used,
          n_exceed = p_info$c[i], resolution_limited = p_info$rl[i],
          E_obs = eff$E_obs, NLES = eff$NLES, NLES_signed = eff$NLES_signed,
          n_leading_edge = st$overlap, hits = paste(ctx$rg[st$le_idx], collapse = " ")
        )
      }))
    }
    rows_m
  }

  results <- if (use_par) {
    furrr::future_map2(
      by_size, sizes, score_size,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    purrr::map2(by_size, sizes, score_size)
  }

  out <- purrr::list_rbind(purrr::list_flatten(results))
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

  # Honesty flag: a pathway capped at b_max without resolving reports a conservative
  # lower-bound p (it can only be called LESS significant than truth). Warn once if
  # any such pathway is at or near significance -- there the imprecision matters.
  near_sig <- out$resolution_limited & out$p_adjust <= alpha * 1.5
  if (any(near_sig, na.rm = TRUE)) {
    warning(sprintf(
      paste0(
        "%d pathway(s) hit the permutation cap (n_perm_max = %d) before resolving; ",
        "their p-values are conservative lower-bounds. Raise n_perm_max (toward ",
        "ceil(target_hits * n_tests / alpha) = %d) if you need them resolved."
      ),
      sum(near_sig), b_max, as.integer(ceiling(target_hits * n_tests / alpha))
    ), call. = FALSE)
  }

  col_order <- c(
    "pathway", "direction", "set_size", "cutoff_rank", "leading_edge_size",
    "overlap", "leading_edge_fraction", "neg_log10_mHG", "mHG_stat",
    "p_value", "p_adjust", "p_adjust_method", "n_perm_used", "n_exceed",
    "resolution_limited", "E_obs", "NLES",
    "NLES_signed", "NES_OHG", "n_leading_edge", "hits"
  )
  if (dir == "up") {
    out$direction <- NULL
    col_order <- setdiff(col_order, "direction")
  }
  out <- out[order(out$p_value), col_order]
  # Resolved cutoff restriction, recorded as metadata (cutoff_rank never exceeds L).
  attr(out, "L_used") <- L
  attr(out, "X_used") <- X
  out
}
