#' Normalized leading-edge score (NLES) and observed magnitude
#'
#' Standardizes the observed leading-edge magnitude against the permutation null
#' of leading-edge magnitudes for the same set size. The magnitude axis is always
#' `abs(weight)`; the sign of `weight` is never read (direction comes from
#' `dir_sign`). The null is conditioned on a non-empty leading edge: permutation
#' draws whose hits all fall below the cutoff cap have no leading-edge magnitude
#' and are excluded, matching the event the observed statistic conditions on. The
#' column -- not the pathway -- is gated to `NA` when the null spread is degenerate
#' or there are too few non-empty draws to standardize reliably.
#'
#' @param le_obs Integer positions of the observed leading-edge hits.
#' @param weight Numeric per-gene magnitude aligned to the ranked list, or `NULL`.
#' @param le_idx_b List of per-permutation leading-edge position vectors.
#' @param dir_sign `+1` (up run) or `-1` (down run).
#' @param robust If `TRUE`, use median/MAD; if `FALSE`, mean/SD.
#' @param min_perm_nles Minimum number of permutations required for a reported NLES.
#' @param min_nles_support Minimum number of distinct null magnitudes required.
#'
#' @return A list with `E_obs`, `NLES`, and `NLES_signed`.
#'
#' @examples
#' null <- ohg_permutation_null(N = 100L, m = 8L, B = 1500L)
#' compute_effect(1:6, runif(100), null$le_idx_b, 1, TRUE, 1000L, 10L)
#'
#' @export
compute_effect <- function(le_obs, weight, le_idx_b, dir_sign, robust,
                           min_perm_nles, min_nles_support) {
  if (is.null(weight)) {
    return(list(E_obs = NA_real_, NLES = NA_real_, NLES_signed = NA_real_))
  }
  summary <- .nles_null_summary(
    le_idx_b, weight, robust, min_perm_nles, min_nles_support
  )
  .nles_from_summary(le_obs, weight, dir_sign, summary)
}

# Single definition of E_b -- the per-permutation null leading-edge
# magnitudes -- from leading-edge indices. `weight` (= ctx$w) is SIGNED, so
# this helper owns the abs; the magnitude axis is always abs(weight). Both the
# fixed-B path (.nles_null_summary) and the adaptive path
# (.adaptive_draw_increment) call this, so their E_b match by construction.
# Internal kernel helper (cf. .mhg_core).
.eb_from_leidx <- function(le_idx_b, weight, robust) {
  w <- abs(weight)
  center <- if (robust) stats::median else mean
  # A null draw whose hits all fall below the cutoff cap L has no top-L leading
  # edge (le_idx = integer(0)) and therefore no leading-edge magnitude. Drop it
  # so the null is the distribution of leading-edge magnitude GIVEN a non-empty
  # leading edge -- the same event the observed statistic conditions on. Coercing
  # an empty draw to center(numeric(0)) = NA instead would poison E_b with NAs,
  # collapsing the spread to NA and silently gating NLES off for every pathway.
  le_idx_b <- le_idx_b[lengths(le_idx_b) > 0L]
  vapply(le_idx_b, function(idx) center(w[idx]), numeric(1))
}

# Summarize a precomputed null E_b vector for ONE set size+direction: spread,
# distinct-count, the stability gates, and the single warning. Depends only
# on (E_b, robust, gates) -- never on the observed pathway -- so every pathway
# of a given size+direction shares it. Splitting this out lets the adaptive
# engine feed an E_b accumulated across rounds (no redraw) while the fixed-B
# path feeds an E_b built fresh from le_idx_b. E_b holds one entry per null draw
# WITH a non-empty leading edge (empty draws are dropped upstream in
# .eb_from_leidx), so B = length(E_b) counts the informative draws and the B-gate
# fires precisely when too few of them survive. It stays correct under any later
# n_perm/min_perm_nles default change. Internal kernel helper (cf. .mhg_core).
.nles_summary_from_Eb <- function(E_b, robust, min_perm_nles, min_nles_support) {
  center <- if (robust) stats::median else mean
  scale_fn <- if (robust) stats::mad else stats::sd

  B <- length(E_b)
  spread <- scale_fn(E_b)
  eps <- .Machine$double.eps^0.5
  n_distinct <- length(unique(E_b))

  degenerate_spread <- !is.finite(spread) || spread < eps
  gated <- B < min_perm_nles || degenerate_spread ||
    n_distinct < min_nles_support
  if (gated) {
    # Check the count gate before the spread gate: when too few non-empty draws
    # survive, B (and an all-empty E_b) drives spread to NA, which would otherwise
    # mis-attribute the cause to the weight vector. Report the true cause first.
    msg <- if (B < min_perm_nles) {
      sprintf(
        paste0(
          "NLES skipped: only %d permutation(s) with a non-empty leading edge ",
          "(< min_perm_nles = %d), too few to size the effect reliably. Increase ",
          "`n_perm`."
        ),
        B, min_perm_nles
      )
    } else if (degenerate_spread) {
      paste0(
        "NLES skipped: the permutation null has near-zero spread (mad(E_b) ~= 0). ",
        "This usually means the `weight` vector has a high density of identical or ",
        "exact-zero values, so random leading edges yield near-constant median ",
        "magnitudes."
      )
    } else {
      sprintf(
        paste0(
          "NLES skipped: the permutation null has only %d distinct leading-edge ",
          "magnitude(s) (< min_nles_support = %d), too thin to standardize against."
        ),
        n_distinct, min_nles_support
      )
    }
    warning(msg, call. = FALSE)
  }

  list(center = center, center_b = center(E_b), spread = spread, gated = gated)
}

# Summarize the permutation null of leading-edge magnitudes for ONE set size
# from a fresh le_idx_b (fixed-B path / compute_effect). Builds E_b via the
# shared .eb_from_leidx() helper, then delegates to .nles_summary_from_Eb().
# The hot path (ohg_enrichment, method = "permutation") computes this once per
# size and reuses it across all pathways of that size. Internal kernel helper
# (cf. .mhg_core).
.nles_null_summary <- function(le_idx_b, weight, robust,
                               min_perm_nles, min_nles_support) {
  E_b <- .eb_from_leidx(le_idx_b, weight, robust)
  .nles_summary_from_Eb(E_b, robust, min_perm_nles, min_nles_support)
}

# Standardize ONE observed leading edge against a precomputed null summary. Cheap
# per-pathway step: just the observed center and a divide. Mirrors the tail of the
# old compute_effect exactly. Internal kernel helper (cf. .mhg_core).
.nles_from_summary <- function(le_obs, weight, dir_sign, summary) {
  if (is.null(weight)) {
    return(list(E_obs = NA_real_, NLES = NA_real_, NLES_signed = NA_real_))
  }
  E_obs <- summary$center(abs(weight)[le_obs])
  if (summary$gated) {
    return(list(E_obs = E_obs, NLES = NA_real_, NLES_signed = NA_real_))
  }
  nles <- (E_obs - summary$center_b) / summary$spread
  list(E_obs = E_obs, NLES = nles, NLES_signed = dir_sign * nles)
}
