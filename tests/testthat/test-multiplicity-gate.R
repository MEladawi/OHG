# Multiplicity-aware ramp gate. The adaptive sampler must stop spending
# permutations on pathways that cannot clear the BH threshold even in their best
# achievable case, while never gating one that could still be significant. The
# headline "spends ~10x fewer permutations" is a large-M property verified on real
# data; here we pin the gate's decision logic deterministically.

test_that("gate decision: hopeless ramping pathways are gated, a strong one is not", {
  # 20 ramping pathways whose best-case p is 0.5 cannot clear BH; one resolved
  # pathway at a tiny p is not ramping and is never gated.
  p_opt <- c(strong = 1e-8, rep(0.5, 20))
  ramping <- c(FALSE, rep(TRUE, 20))
  gated <- .gate_decisions(p_opt, ramping, alpha = 0.05)
  expect_false(gated[1])
  expect_true(all(gated[-1]))
})

test_that("gate decision: a ramping pathway that could still clear alpha is kept", {
  # Few tests => BH lenient; a best-case p of 0.01 across 5 tests can still reach
  # padj < 0.05, so the gate must keep it alive.
  p_opt <- c(1e-8, 0.01, 0.9, 0.9, 0.9)
  ramping <- c(FALSE, TRUE, TRUE, TRUE, TRUE)
  gated <- .gate_decisions(p_opt, ramping, alpha = 0.05)
  expect_false(gated[2]) # padj(0.01) = 0.01 * 5 / 2 = 0.025 <= 0.05 -> keep
  expect_true(all(gated[3:5]))
})

test_that("gate decision: never gates a finalized (non-ramping) pathway", {
  expect_false(any(.gate_decisions(c(0.9, 0.9), c(FALSE, FALSE), alpha = 0.05)))
})

test_that("Clopper-Pearson lower bound: 0 at c=0, rises with c, below the estimate", {
  gamma <- 0.05 / (2 * 200) # gate_conf / (2m)
  # c = 0 cannot be ruled out as significant -> exactly 0 (never auto-accepted).
  expect_equal(.p_lower_cp(0L, 200L, gamma), 0)
  # c > 0: a positive lower bound, below the point estimate (c + 1)/(n + 1).
  lb <- .p_lower_cp(8L, 200L, gamma)
  expect_gt(lb, 0)
  expect_lt(lb, (8 + 1) / (200 + 1))
  # tightens (rises) as more draws accumulate at the same exceedance rate.
  expect_gt(.p_lower_cp(80L, 2000L, gamma), .p_lower_cp(8L, 200L, gamma))
})

test_that("the gate accepts a decided ramping pathway, not an undecided one", {
  b_max <- 10000L
  th <- 10L
  gamma <- 0.05 / (2 * 203) # gate_conf / (2m), m = total rows below
  # One resolved strong, one undecided (c = 0, lower bound 0), one decided ramping
  # (c = 8 of 200 -> lower bound well above the threshold set by a bulk of large-p
  # finalized rows).
  g <- list(
    n = 3L, status = c("resolved", "ramping", "ramping"),
    c = c(60L, 0L, 8L), b_used = 200L, L_hit = c(30L, NA_integer_, NA_integer_)
  )
  bulk <- list(
    n = 200L, status = rep("accepted", 200L), c = rep(100L, 200L),
    b_used = 200L, L_hit = rep(NA_integer_, 200L)
  )
  a <- .group_p_opt(g, th, b_max, gamma)
  b <- .group_p_opt(bulk, th, b_max, gamma)
  accept <- .gate_decisions(c(a$p_opt, b$p_opt), c(a$ramping, b$ramping), 0.05)
  expect_false(accept[1]) # resolved -> not ramping
  expect_false(accept[2]) # c = 0 -> lower bound 0 -> could be significant -> keep
  expect_true(accept[3]) # c = 8 decided non-significant among the bulk -> accepted
})

test_that("the gate never touches a strongly significant pathway", {
  ranked <- paste0("g", 1:600)
  fillers <- stats::setNames(
    lapply(seq_len(40), function(i) ranked[sample.int(600, 12L)]),
    paste0("F", seq_len(40))
  )
  sets <- c(list(STRONG = ranked[1:14]), fillers)
  # alpha = 1 disables the gate (padj can never exceed 1); 0.05 enables it.
  gate <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 1000L, target_hits = 10L,
    n_perm_max = 20000L, alpha = 0.05, seed = 7
  )
  nogate <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 1000L, target_hits = 10L,
    n_perm_max = 20000L, alpha = 1, seed = 7
  )
  g <- gate[gate$pathway == "STRONG", ]
  n <- nogate[nogate$pathway == "STRONG", ]
  expect_equal(g$p_value, n$p_value)
  expect_equal(g$resolution_limited, n$resolution_limited)
  expect_equal(g$n_exceed, n$n_exceed)
})

test_that("gate preserves the significant set and never spends more permutations", {
  ranked <- paste0("g", 1:600)
  fillers <- stats::setNames(
    lapply(seq_len(40), function(i) ranked[sample.int(600, 12L + (i %% 15L))]),
    paste0("F", seq_len(40))
  )
  sets <- c(list(STRONG = ranked[1:14]), fillers)
  gate <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 1000L, target_hits = 10L,
    n_perm_max = 20000L, alpha = 0.05, seed = 9
  )
  nogate <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 1000L, target_hits = 10L,
    n_perm_max = 20000L, alpha = 1, seed = 9
  )
  expect_equal(
    sort(gate$pathway[gate$p_adjust <= 0.05]),
    sort(nogate$pathway[nogate$p_adjust <= 0.05])
  )
  dg <- sum(unique(gate[, c("set_size", "n_perm_used")])$n_perm_used)
  dn <- sum(unique(nogate[, c("set_size", "n_perm_used")])$n_perm_used)
  expect_lte(dg, dn)
})

test_that("gated run stays reproducible across worker counts", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  ranked <- paste0("g", 1:800)
  set.seed(22)
  sets <- c(
    list(A = ranked[1:15], B = ranked[c(1:8, 400:406)]),
    stats::setNames(
      lapply(1:20, function(k) ranked[sample.int(800, 20L + k)]),
      paste0("N", 1:20)
    )
  )
  r1 <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 800L, n_perm_max = 12000L, seed = 3, n_cores = 1L
  )
  r2 <- ohg_enrichment_quiet(ranked, sets,
    method = "adaptive", n_perm = 800L, n_perm_max = 12000L, seed = 3, n_cores = 2L
  )
  expect_equal(r1$p_value[order(r1$pathway)], r2$p_value[order(r2$pathway)])
  expect_equal(r1$n_perm_used[order(r1$pathway)], r2$n_perm_used[order(r2$pathway)])
})

test_that("CALIBRATION holds with the gate active (level control under H0)", {
  set.seed(404)
  ranked <- paste0("g", 1:200)
  rs <- sort(rnorm(200), decreasing = TRUE)
  sets <- stats::setNames(
    lapply(seq_len(120), function(i) ranked[sample.int(200, 15)]),
    paste0("S", seq_len(120))
  )
  res <- ohg_enrichment_quiet(ranked, sets,
    rank_stat = rs, method = "adaptive", n_perm = 500L,
    target_hits = 10L, n_perm_max = 3000L, alpha = 0.05, seed = 1
  )
  # Sequential-MC validity: not anti-conservative, roughly centred (the gate only
  # coarsens hopeless rows, so it cannot inflate the false-positive rate).
  expect_lte(mean(res$p_value <= 0.05), 0.10)
  expect_lte(mean(res$p_value <= 0.10), 0.18)
  expect_gt(mean(res$p_value), 0.45)
})
