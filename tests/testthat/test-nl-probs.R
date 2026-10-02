# One individual's nested logit probabilities: nl_individual_probs() works on
# raw arrays in per-thread buffers, and must reproduce bit for bit the
# Armadillo expressions it replaced, kept as a reference in
# src/kernel_test_exports.cpp.

nl_probs_fields <- c("P_i", "P_j_given_k", "P_k", "log_I_k", "log_P_i",
                     "log_P_outside")

expect_nl_probs_identical <- function(V, nest0, lambda, ioo, what) {
  r <- choicer:::test_nl_individual_probs(V, nest0, lambda, ioo)
  for (f in nl_probs_fields) {
    expect_true(identical(r[[f]], r[[paste0(f, "_ref")]], num.eq = FALSE),
                label = sprintf("%s (outside option %s): %s", what, ioo, f))
  }
}

test_that("nl_individual_probs() reproduces its Armadillo reference", {
  set.seed(20261002)
  for (rep in seq_len(60)) {
    n_nests <- sample(c(1:6, 17L, 25L), 1L)
    m <- sample(c(1:8, 15:18, 40L), 1L)
    # Nests drawn from a subset, so some nests are empty and some singletons
    used <- sort(sample.int(n_nests, sample.int(n_nests, 1L)))
    nest0 <- as.integer(used[sample.int(length(used), m, replace = TRUE)] - 1L)
    lambda <- runif(n_nests, 0.05, 1)
    lambda[sample.int(n_nests, max(1L, n_nests %/% 3L))] <- 1
    V <- rnorm(m, sd = sample(c(0.5, 5, 50), 1L))
    for (ioo in c(FALSE, TRUE)) {
      expect_nl_probs_identical(V, nest0, lambda, ioo, sprintf("draw %d", rep))
    }
  }
})

test_that("nl_individual_probs() matches its reference, every nest occupied", {
  # Many alternatives per nest and moderate utilities: the regime in which the
  # nest sum's order (two accumulators, the odd tail, then the outside term)
  # changes the last bit, so a reordering would show here.
  set.seed(20261003)
  for (rep in seq_len(60)) {
    n_nests <- sample(c(2:9, 16:18), 1L)
    m <- n_nests * sample(2:4, 1L)
    nest0 <- sample(c(seq_len(n_nests) - 1L,
                      sample.int(n_nests, m - n_nests, replace = TRUE) - 1L))
    lambda <- runif(n_nests, 0.3, 1)
    V <- rnorm(m, sd = 0.5)
    for (ioo in c(FALSE, TRUE)) {
      expect_nl_probs_identical(V, nest0, lambda, ioo,
                                sprintf("dense draw %d", rep))
    }
  }
})

test_that("nl_individual_probs() matches its reference on edge inputs", {
  edge <- list(
    list("one alternative", 1.3, 0L, 0.7),
    list("lambda near 0", c(0.4, -0.2, 1.1), c(0L, 0L, 1L), c(1e-3, 1)),
    list("a -Inf utility", c(-Inf, 0.3, 0.1), c(0L, 0L, 1L), c(0.6, 0.9)),
    list("a whole nest at -Inf", c(-Inf, -Inf, 0.2), c(0L, 0L, 1L),
         c(0.6, 0.9)),
    list("every utility -Inf", c(-Inf, -Inf), c(0L, 1L), c(0.5, 1)),
    list("a +Inf utility", c(Inf, 0.3, 0.1), c(0L, 0L, 1L), c(0.6, 0.9)),
    list("a NaN utility", c(NaN, 0.3, 0.1), c(0L, 1L, 1L), c(0.6, 0.9)),
    list("large utilities", c(800, 790, -700, 760), c(0L, 1L, 1L, 2L),
         c(0.3, 0.8, 1)),
    list("empty nests", c(0.1, 0.2), c(3L, 3L), c(0.5, 0.6, 0.7, 0.8)),
    list("signed zeros", c(0, -0, 0), c(0L, 0L, 1L), c(1, 1)),
    # Three singleton nests whose sums change in the last bit when the nest
    # sum uses one accumulator, adds the outside term first, or puts the odd
    # term in the second accumulator (found on arm64)
    list("three singletons, a", c(1.01, -0.07, -1.14), 0:2, c(1, 1, 1)),
    list("three singletons, b", c(0.35, 1.17, -0.48), 0:2, c(1, 1, 1)),
    list("three singletons, c", c(0.19, -0.03, 0.47), 0:2, c(1, 1, 1))
  )
  for (e in edge) {
    for (ioo in c(FALSE, TRUE)) {
      expect_nl_probs_identical(e[[2]], e[[3]], e[[4]], ioo, e[[1]])
    }
  }
})

test_that("nl_blp_contraction() checks nest codes and lambda before its loop", {
  # It takes lambda directly and never runs nl_parse_theta(), so these inputs
  # used to fail inside the parallel loop (a terminate under OpenMP).
  d <- create_nl_inputs()
  J <- length(d$nest_idx)
  n_nests <- max(d$nest_idx)
  args <- list(delta = rep(0, J), target_shares = rep(1 / J, J), X = d$X,
               beta = c(0.1, -0.1), lambda = c(1, 0.7, 0.8)[seq_len(n_nests)],
               alt_idx = d$alt_idx, nest_idx = d$nest_idx, M = d$M,
               weights = d$weights)
  bad_code <- args
  bad_code$nest_idx[2L] <- 0L
  expect_error(do.call(nl_blp_contraction, bad_code),
               "nest_idx must use 1-based nest indices (found 0).",
               fixed = TRUE)
  # Also for an alternative that no choice situation offers, as in the other
  # nested logit kernels (it used to be ignored)
  unoffered <- args
  unoffered$nest_idx <- c(args$nest_idx, 0L)
  expect_error(do.call(nl_blp_contraction, unoffered),
               "nest_idx must use 1-based nest indices (found 0).",
               fixed = TRUE)
  for (n_lambda in c(n_nests - 1L, n_nests + 1L)) {
    wrong <- args
    wrong$lambda <- rep(0.8, n_lambda)
    expect_error(do.call(nl_blp_contraction, wrong),
                 sprintf(paste("lambda must have one entry per nest of",
                               "nest_idx (%d); it has %d."),
                         n_nests, n_lambda), fixed = TRUE)
  }
})
