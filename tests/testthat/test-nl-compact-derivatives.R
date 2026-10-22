# The nested logit BHHH kernel on designs with many alternatives and sparse
# choice sets. A situation's score is nonzero only on the coefficients, the
# dissimilarities of the nests it offers and the constants of the
# alternatives it offers, so the kernel forms w s s' in that compact block
# and adds it into the n x n result. These designs (30 alternatives in eight
# nests of three and six singletons, choice sets of 3-8, dissimilarities on
# both sides of one) exercise the slot maps: nests offered in part or not at
# all, singleton nests (their dissimilarity fixed at one), many small nests,
# one large nest, a situation listing an alternative twice, the outside
# option, no constants.

# N situations, each offering `m_range` distinct alternatives drawn from
# 1..J, plus one repeated row when `dup`; weights vary across situations.
ncd_fixture <- function(name, seed, N = 60L, J = 30L, m_range = c(3L, 8L),
                        K = 2L, nest_idx = as.integer(c(rep(1:8, each = 3L), 9:14)),
                        include_outside_option = FALSE, use_asc = TRUE,
                        dup = FALSE, shuffle = FALSE) {
  set.seed(seed)
  alts <- lapply(seq_len(N), function(t) {
    m <- sample(m_range[1L]:m_range[2L], 1L)
    a <- sort(sample.int(J, min(m, J)))
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    if (shuffle) a <- a[sample.int(length(a))]
    a
  })
  M <- lengths(alts)
  n_rows <- sum(M)
  choice_idx <- vapply(M, function(m) {
    if (include_outside_option) sample(0:m, 1L) else sample.int(m, 1L)
  }, integer(1L))
  K_l <- sum(tabulate(nest_idx) > 1L)
  n_delta <- if (!use_asc) 0L else if (include_outside_option) J else J - 1L
  list(name = name, X = matrix(stats::rnorm(n_rows * K), n_rows, K),
       alt_idx = as.integer(unlist(alts)), choice_idx = as.integer(choice_idx),
       nest_idx = as.integer(nest_idx), M = as.integer(M),
       weights = stats::runif(N, 0.5, 2),
       theta = c(stats::rnorm(K, sd = 0.4), stats::runif(K_l, 0.55, 1.2),
                 stats::rnorm(n_delta, sd = 0.4)),
       use_asc = use_asc, include_outside_option = include_outside_option,
       J = J, K = K, K_l = K_l, N = N)
}

ncd_cells <- function() {
  list(
    ncd_fixture("base", 3001),
    ncd_fixture("outside option", 3002, include_outside_option = TRUE),
    ncd_fixture("repeated alternative", 3003, dup = TRUE),
    ncd_fixture("repeated alternative, shuffled rows", 3004, dup = TRUE,
                shuffle = TRUE),
    ncd_fixture("no constants", 3005, use_asc = FALSE),
    ncd_fixture("one coefficient", 3006, K = 1L),
    ncd_fixture("many small nests", 3007, nest_idx = rep(1:15, each = 2L)),
    ncd_fixture("one large nest and singletons", 3008,
                nest_idx = c(rep(1L, 20L), 2:11)),
    ncd_fixture("wide sets", 3009, J = 60L, m_range = c(10L, 40L), N = 30L,
                nest_idx = rep(1:20, each = 3L))
  )
}

ncd_call <- function(kernel, f, theta = f$theta, weights = f$weights, ...) {
  ns <- asNamespace("choicer")
  a <- list(theta = theta, X = f$X, alt_idx = f$alt_idx,
            choice_idx = f$choice_idx, nest_idx = f$nest_idx, M = f$M,
            weights = weights, use_asc = f$use_asc,
            include_outside_option = f$include_outside_option)
  switch(kernel,
    gradient = do.call(ns$nl_loglik_gradient_parallel, c(a, list(...))),
    bhhh = do.call(ns$nl_bhhh_parallel, c(a, list(...))),
    scores = do.call(ns$nl_scores_parallel, c(a[-7L], list(...))))
}

# The dense BHHH update, situation by situation: B + w_i s_i s_i'.
ncd_bhhh_ref <- function(S, w) {
  B <- matrix(0, ncol(S), ncol(S))
  for (i in seq_len(nrow(S))) B <- B + w[i] * outer(S[i, ], S[i, ])
  B
}

# Equal non-finite patterns (missing values, infinities and their signs), and
# finite entries within `tol`; off CRAN also NA and NaN apart.
ncd_same_pattern <- function(a, b, what, tol = 1e-12) {
  a <- unclass(a)
  expect_identical(is.na(a), is.na(b), label = paste(what, "NA or NaN"))
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    expect_identical(is.nan(a), is.nan(b), label = paste(what, "NaN"))
  }
  expect_identical(is.infinite(a), is.infinite(b), label = paste(what, "Inf"))
  expect_identical(a[is.infinite(a)], b[is.infinite(b)],
                   label = paste(what, "infinities and their signs"))
  fin <- is.finite(b)
  if (any(fin)) {
    expect_equal(a[fin], b[fin], tolerance = tol, label = paste(what, "finite entries"))
  }
}

test_that("the compact NL BHHH is the weighted cross-product of the scores", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in ncd_cells()) {
    for (nt in 1:2) {
      set_omp_threads(nt)
      S <- ncd_call("scores", f)
      B <- ncd_call("bhhh", f)
      expect_identical(B, t(B), label = paste(f$name, "symmetric"))
      expect_equal(unclass(B), crossprod(sqrt(f$weights) * S), tolerance = 1e-12,
                   label = sprintf("[%s, %d threads]", f$name, nt))
    }
  }
})

test_that("at one thread the compact NL BHHH is the dense update bit for bit", {
  # Each entry is the product of the two scores, times the weight, added in
  # situation order, rounded at each step as R rounds the reference.
  skip_on_cran()  # a compiler that fuses the dense update's terms differs
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (f in ncd_cells()) {
    S <- ncd_call("scores", f)
    for (w in list(f$weights, f$weights^2)) {
      expect_identical(unclass(ncd_call("bhhh", f, weights = w)),
                       ncd_bhhh_ref(S, w), label = f$name)
    }
  }
})

test_that("relabeling alternatives and nests permutes the NL BHHH", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f <- ncd_fixture("relabel", 3011, dup = TRUE)
  set.seed(3111)
  # alternatives: code c becomes perm[c] (the reference stays first), its
  # nest and constant with it
  perm <- c(1L, 1L + sample.int(f$J - 1L))
  g <- f
  g$alt_idx <- perm[f$alt_idx]
  g$nest_idx[perm] <- f$nest_idx
  g$theta[f$K + f$K_l + perm[-1L] - 1L] <- f$theta[f$K + f$K_l + seq_len(f$J - 1L)]
  idx <- c(seq_len(f$K + f$K_l), f$K + f$K_l + perm[-1L] - 1L)
  expect_identical(ncd_call("bhhh", g)[idx, idx], ncd_call("bhhh", f))
  # nests: nest code k becomes sigma[k] among the non-singleton nests, so the
  # dissimilarities' order in theta changes with it (and the probabilities'
  # sum over the nests their order, hence the tolerance)
  sigma <- c(sample.int(8L), 9:14)
  h <- f
  h$nest_idx <- sigma[f$nest_idx]
  h$theta[f$K + sigma[1:8]] <- f$theta[f$K + 1:8]
  idx <- c(seq_len(f$K), f$K + sigma[1:8], f$K + f$K_l + seq_len(f$J - 1L))
  expect_equal(unclass(ncd_call("bhhh", h))[idx, idx], unclass(ncd_call("bhhh", f)),
               tolerance = 1e-12)
})

test_that("non-finite scores and weights spread over the NL BHHH as the dense update did", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f0 <- ncd_fixture("non-finite", 3014)
  row_off <- c(0L, cumsum(f0$M))
  i <- 3L
  r_other <- row_off[i] + setdiff(seq_len(f0$M[i]), f0$choice_idx[i])[1L]
  cells <- list(
    # utilities that overflow: NaN probabilities and scores
    overflow = function(f) { f$X[r_other, 1L] <- 1e308; f$theta[1L] <- 2; f },
    inf_weight = function(f) { f$weights[i] <- Inf; f },
    nan_weight = function(f) { f$weights[i] <- NaN; f },
    na_weight = function(f) { f$weights[i] <- NA_real_; f })
  for (nm in names(cells)) {
    f <- cells[[nm]](f0)
    S <- ncd_call("scores", f)
    B <- ncd_call("bhhh", f)
    expect_true(any(!is.finite(B)), label = paste(nm, "has non-finite entries"))
    ncd_same_pattern(B, ncd_bhhh_ref(S, f$weights), paste(nm, "BHHH"))
  }
})

test_that("a shared NL BHHH result equals the per-thread accumulators' sum", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in ncd_cells()[c(1L, 2L, 3L, 5L, 7L, 9L)]) {
    set_omp_threads(2L)
    shared <- ncd_call("bhhh", f, acc_bytes = 0)
    expect_identical(shared, t(shared), label = paste(f$name, "symmetric"))
    expect_equal(unclass(shared), unclass(ncd_call("bhhh", f)), tolerance = 1e-12,
                 label = paste(f$name, "shared vs per-thread"))
    # one thread adds into the result itself, whatever the budget
    set_omp_threads(1L)
    expect_identical(ncd_call("bhhh", f, acc_bytes = 0), ncd_call("bhhh", f),
                     label = paste(f$name, "one thread"))
  }
  set_omp_threads(2L)
  f0 <- ncd_fixture("non-finite, shared", 3015)
  row_off <- c(0L, cumsum(f0$M))
  r_other <- row_off[5L] + setdiff(seq_len(f0$M[5L]), f0$choice_idx[5L])[1L]
  f <- f0
  f$X[r_other, 1L] <- 1e308
  f$theta[1L] <- 2
  cells <- list(overflow = f)
  for (w_bad in c(Inf, NA_real_)) {
    f <- f0
    f$weights[5L] <- w_bad
    cells[[paste("weight", w_bad)]] <- f
  }
  for (nm in names(cells)) {
    ncd_same_pattern(ncd_call("bhhh", cells[[nm]], acc_bytes = 0),
                     unclass(ncd_call("bhhh", cells[[nm]])), paste(nm, "shared"))
  }
})
