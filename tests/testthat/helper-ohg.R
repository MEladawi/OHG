# Several tests use strong synthetic gene sets that cap the adaptive sampler at
# n_perm_max and so emit the honest "permutation cap" warning. Those tests assert
# direction / BH / parallel / schema relationships unrelated to that warning, so
# they call through this quiet wrapper to keep the suite output clean. Tests that
# are specifically about the warning use expect_warning() on ohg_enrichment().
ohg_enrichment_quiet <- function(...) suppressWarnings(ohg_enrichment(...))
