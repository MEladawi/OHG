# Adaptive (Besag-Clifford 1991 sequential Monte-Carlo) p-values.

test_that("Besag-Clifford estimator: h/L when resolved, (c+1)/(B_max+1) when capped", {
  # exceedance = log_stat_b <= obs (smaller log = more enriched); obs = 0 here.
  stream <- c(1, -1, 1, -1, -1, 1, 1, 1) # <= 0 at draws 2, 4, 5 -> cumsum hits 2 at draw 4
  r <- .bc_finalize(stream, obs = 0, h = 2L, b_max = 100L)
  expect_false(r$resolution_limited)
  expect_equal(r$p_hat, 2 / 4) # 2nd exceedance lands at draw 4
  expect_equal(r$n_exceed, 3L)

  # only 3 exceedances but h = 5 -> capped, conservative (c + 1)/(B_max + 1)
  r2 <- .bc_finalize(stream, obs = 0, h = 5L, b_max = 100L)
  expect_true(r2$resolution_limited)
  expect_equal(r2$n_exceed, 3L)
  expect_equal(r2$p_hat, (3 + 1) / (100 + 1))
})

test_that("adaptive removes the 1/(B0+1) floor for a strong pathway", {
  ranked <- paste0("g", 1:400)
  fillers <- setNames(
    lapply(seq_len(50), function(i) ranked[sample.int(400, 15)]),
    paste0("F", seq_len(50))
  )
  sets <- c(list(STRONG = ranked[1:15]), fillers) # STRONG = exact top block

  fixed <- ohg_enrichment(ranked, sets, n_perm = 2000L, method = "permutation", seed = 1)
  adap <- suppressWarnings(ohg_enrichment(ranked, sets,
    n_perm = 2000L, method = "adaptive",
    target_hits = 10L, n_perm_max = 20000L, seed = 1
  ))
  pf <- fixed$p_value[fixed$pathway == "STRONG"]
  pa <- adap$p_value[adap$pathway == "STRONG"]
  expect_lt(pa, pf) # adaptive resolves below the fixed-B value
  expect_lt(pa, 1 / (2000 + 1)) # ... and below the fixed-B floor itself
})

test_that("adaptive spends draws only where needed (per m-group)", {
  ranked <- paste0("g", 1:400)
  # different sizes => different m-groups: STRONG (m=15, top block) expands;
  # NULLBIG (m=40, spread across the whole list => not enriched) does not.
  sets <- list(
    STRONG = ranked[1:15],
    NULLBIG = ranked[round(seq(2, 399, length.out = 40))]
  )
  res <- suppressWarnings(ohg_enrichment(ranked, sets,
    n_perm = 1000L, method = "adaptive",
    target_hits = 10L, n_perm_max = 20000L, seed = 1
  ))
  expect_equal(res$n_perm_used[res$pathway == "NULLBIG"], 1000L) # finalized at B0
  expect_gt(res$n_perm_used[res$pathway == "STRONG"], 1000L) # expanded
})

test_that("same-m pathways share one draw sequence (one n_perm_used)", {
  ranked <- paste0("g", 1:400)
  sets <- list(STRONG = ranked[1:15], MILD = ranked[c(1:5, 200:209)]) # both m = 15
  res <- suppressWarnings(ohg_enrichment(ranked, sets,
    n_perm = 1000L, method = "adaptive",
    target_hits = 10L, n_perm_max = 50000L, seed = 1
  ))
  expect_equal(
    res$n_perm_used[res$pathway == "STRONG"],
    res$n_perm_used[res$pathway == "MILD"]
  )
})

test_that("capped pathway is flagged resolution_limited and warns near significance", {
  ranked <- paste0("g", 1:400)
  fillers <- setNames(
    lapply(seq_len(30), function(i) ranked[sample.int(400, 15)]),
    paste0("F", seq_len(30))
  )
  sets <- c(list(STRONG = ranked[1:15]), fillers)
  expect_warning(
    res <- ohg_enrichment(ranked, sets,
      n_perm = 500L, method = "adaptive", target_hits = 10L,
      n_perm_max = 2000L, seed = 1
    ),
    "permutation cap"
  )
  expect_true(res$resolution_limited[res$pathway == "STRONG"])
})

test_that("CALIBRATION: adaptive p_hat ~ uniform under a random ranking", {
  set.seed(11)
  ranked <- paste0("g", 1:200)
  rs <- sort(rnorm(200), decreasing = TRUE) # random sets vs any ranking => H0
  sets <- setNames(
    lapply(seq_len(80), function(i) ranked[sample.int(200, 15)]),
    paste0("S", seq_len(80))
  )
  res <- ohg_enrichment(ranked, sets,
    rank_stat = rs, method = "adaptive", n_perm = 500L,
    target_hits = 10L, n_perm_max = 3000L, seed = 1
  )
  # The correct H0 invariant for a sequential Monte-Carlo p-value is VALIDITY
  # (level control), not exact uniformity: the mHG statistic is discrete on a
  # small list and Besag-Clifford is mildly conservative, so the right check is
  # that the p-values are not anti-conservative and are roughly centred.
  expect_lte(mean(res$p_value <= 0.05), 0.10) # ~alpha (with Monte-Carlo slack)
  expect_lte(mean(res$p_value <= 0.10), 0.17)
  expect_lte(mean(res$p_value <= 0.20), 0.28)
  expect_gt(mean(res$p_value), 0.45)
  expect_lt(mean(res$p_value), 0.62)
})

test_that("adaptive is reproducible across worker counts", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  ranked <- paste0("g", 1:300)
  sets <- list(A = ranked[1:15], B = ranked[c(1:8, 200:206)], C = ranked[50:69])
  r1 <- suppressWarnings(ohg_enrichment(ranked, sets,
    method = "adaptive", n_perm = 800L,
    n_perm_max = 8000L, seed = 3, n_cores = 1L
  ))
  r2 <- suppressWarnings(ohg_enrichment(ranked, sets,
    method = "adaptive", n_perm = 800L,
    n_perm_max = 8000L, seed = 3, n_cores = 2L
  ))
  expect_equal(r1$p_value[order(r1$pathway)], r2$p_value[order(r2$pathway)])
  expect_equal(r1$n_perm_used[order(r1$pathway)], r2$n_perm_used[order(r2$pathway)])
})

test_that("method 'permutation' reproduces the fixed-B p-value exactly", {
  ranked <- paste0("g", 1:250)
  sets <- list(A = ranked[1:12], B = ranked[c(2:9, 100, 150)])
  off <- ohg_enrichment(ranked, sets, method = "permutation", n_perm = 600L, seed = 99)
  expect_true(all(off$n_perm_used == 600L))
  expect_true(all(!off$resolution_limited))
  expect_equal(off$p_value, (1 + off$n_exceed) / (1 + 600L))
})
