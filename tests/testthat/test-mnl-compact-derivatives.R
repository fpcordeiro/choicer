# The multinomial logit Hessian, BHHH and cluster-meat kernels on designs
# with many alternatives and sparse choice sets. A situation's score and
# Hessian are nonzero only on the coefficients and the free constants of the
# alternatives it offers, so the kernels form them in that compact block and
# add it into the n x n result (the meat a cluster's summed score, over the
# parameters its situations touch). These designs (30 alternatives, choice
# sets of 3-8) make the blocks small next to n and exercise the slot maps:
# the normalized reference, the outside option, a situation listing an
# alternative twice (adjacent or not), the reference never offered,
# situations offering only the reference, choice sets of one alternative, no
# constants.

# N situations, each offering `m_range` distinct alternatives drawn from 1..J
# (from 2..J when `absent_ref`), plus one repeated row when `dup`; weights
# vary across situations.
mcd_fixture <- function(name, seed, N = 60L, J = 30L, m_range = c(3L, 8L),
                        K = 2L, include_outside_option = FALSE, use_asc = TRUE,
                        dup = FALSE, absent_ref = FALSE, ref_only = FALSE,
                        shuffle = FALSE, single = FALSE) {
  set.seed(seed)
  pool <- if (absent_ref) 2:J else seq_len(J)
  alts <- lapply(seq_len(N), function(t) {
    # `ref_only`: about a third of the situations offer only the reference
    if (ref_only && stats::runif(1L) < 1 / 3) return(1L)
    # `single`: about a third offer one alternative
    if (single && stats::runif(1L) < 1 / 3) return(pool[sample.int(length(pool), 1L)])
    m <- sample(m_range[1L]:m_range[2L], 1L)
    a <- sort(pool[sample.int(length(pool), min(m, length(pool)))])
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    # `shuffle`: rows in no particular order (repeats not adjacent)
    if (shuffle) a <- a[sample.int(length(a))]
    a
  })
  M <- lengths(alts)
  n_rows <- sum(M)
  choice_idx <- vapply(M, function(m) {
    if (include_outside_option) sample(0:m, 1L) else sample.int(m, 1L)
  }, integer(1L))
  n_delta <- if (!use_asc) 0L else if (include_outside_option) J else J - 1L
  list(name = name, X = matrix(stats::rnorm(n_rows * K), n_rows, K),
       alt_idx = as.integer(unlist(alts)), choice_idx = as.integer(choice_idx),
       M = as.integer(M), weights = stats::runif(N, 0.5, 2),
       theta = c(stats::rnorm(K, sd = 0.4), stats::rnorm(n_delta, sd = 0.4)),
       use_asc = use_asc, include_outside_option = include_outside_option,
       J = J, K = K, N = N)
}

mcd_cells <- function() {
  list(
    mcd_fixture("base", 2001),
    mcd_fixture("outside option", 2002, include_outside_option = TRUE),
    mcd_fixture("repeated alternative", 2003, dup = TRUE),
    mcd_fixture("repeated alternative, shuffled rows", 2004, dup = TRUE,
                shuffle = TRUE),
    mcd_fixture("reference never offered", 2005, absent_ref = TRUE),
    mcd_fixture("situations with the reference only", 2006, ref_only = TRUE),
    mcd_fixture("single-alternative sets", 2007, single = TRUE),
    mcd_fixture("no constants", 2008, use_asc = FALSE),
    mcd_fixture("one coefficient", 2009, K = 1L),
    mcd_fixture("wide sets", 2010, J = 60L, m_range = c(10L, 40L), N = 30L)
  )
}

mcd_call <- function(kernel, f, theta = f$theta, weights = f$weights, ...) {
  ns <- asNamespace("choicer")
  a <- list(theta = theta, X = f$X, alt_idx = f$alt_idx,
            choice_idx = f$choice_idx, M = f$M, weights = weights,
            use_asc = f$use_asc, include_outside_option = f$include_outside_option)
  switch(kernel,
    gradient = do.call(ns$mnl_loglik_gradient_parallel, c(a, list(...))),
    hessian = do.call(ns$mnl_loglik_hessian_parallel, c(a, list(...))),
    bhhh = do.call(ns$mnl_bhhh_parallel, c(a, list(...))),
    scores = do.call(ns$mnl_scores_parallel, c(a[-6L], list(...))),
    meat = do.call(ns$mnl_cluster_meat_parallel, c(a, list(...))))
}

# The dense BHHH update, situation by situation: B + w_i s_i s_i'.
mcd_bhhh_ref <- function(S, w) {
  B <- matrix(0, ncol(S), ncol(S))
  for (i in seq_len(nrow(S))) B <- B + w[i] * outer(S[i, ], S[i, ])
  B
}

# Equal non-finite patterns (missing values, infinities and their signs), and
# finite entries within `tol`. Off CRAN also NA and NaN apart: R's NA keeps
# its payload through the kernels' arithmetic and exp() here, as it does
# through R's, which not every platform's libm guarantees.
mcd_same_pattern <- function(a, b, what, tol = 1e-12) {
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

# The free constant of each row (its position in theta), or NA.
mcd_free <- function(f) {
  if (!f$use_asc) return(rep(NA_integer_, length(f$alt_idx)))
  j <- if (f$include_outside_option) f$alt_idx else f$alt_idx - 1L
  ifelse(j >= 1L, f$K + j, NA_integer_)
}

# The Hessian of the negated log-likelihood as the dense kernel formed it.
# Situation by situation, with p the probabilities, mu_b = sum_a p_a x_a and
# mu_d the probabilities summed by constant (zero off the situation's
# constants): BB = sum_a p_a x_a x_a' - mu_b mu_b', BD = sum_a p_a x_a e_a' -
# mu_b mu_d' and DD = diag(mu_d) - mu_d mu_d', over every constant, times
# w_i. A non-finite mean therefore fills its parameter's row and column, and
# a non-finite weight every entry (w * 0).
mcd_dense_hessian <- function(f, weights = f$weights) {
  n <- length(f$theta); K <- f$K; J_d <- n - K
  b <- seq_len(K); d <- K + seq_len(J_d)
  fr <- mcd_free(f)
  row_off <- c(0L, cumsum(f$M))
  H <- matrix(0, n, n)
  for (i in seq_len(f$N)) {
    r <- (row_off[i] + 1L):row_off[i + 1L]
    Xi <- f$X[r, , drop = FALSE]
    v <- drop(Xi %*% f$theta[b]) + ifelse(is.na(fr[r]), 0, f$theta[fr[r]])
    if (f$include_outside_option) v <- c(0, v)
    e <- exp(v - if (all(is.na(v))) -Inf else max(v, na.rm = TRUE))
    p <- e / sum(e)
    if (f$include_outside_option) p <- p[-1L]
    BB <- crossprod(Xi * p, Xi)
    BB[lower.tri(BB)] <- t(BB)[lower.tri(BB)]
    mu_b <- colSums(Xi * p)
    H[b, b] <- H[b, b] + (BB - outer(mu_b, mu_b)) * weights[i]
    if (J_d > 0L) {
      BD <- matrix(0, K, J_d)
      mu_d <- numeric(J_d)
      for (a in which(!is.na(fr[r]))) {
        j <- fr[r[a]] - K
        BD[, j] <- BD[, j] + p[a] * Xi[a, ]
        mu_d[j] <- mu_d[j] + p[a]
      }
      BD <- BD - outer(mu_b, mu_d)
      H[b, d] <- H[b, d] + BD * weights[i]
      H[d, b] <- H[d, b] + t(BD) * weights[i]
      H[d, d] <- H[d, d] + (diag(mu_d, J_d) - outer(mu_d, mu_d)) * weights[i]
    }
  }
  H
}

mcd_scale_dev <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))

test_that("the compact Hessian matches numDeriv of the gradient on sparse designs", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  for (f in mcd_cells()) {
    num <- numDeriv::jacobian(function(th) mcd_call("gradient", f, theta = th)$gradient,
                              f$theta)
    expect_lt(mcd_scale_dev(unclass(mcd_call("hessian", f)), num), 1e-6,
              label = f$name)
  }
})

test_that("the compact Hessian is the dense one", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in mcd_cells()) {
    for (nt in 1:2) {
      set_omp_threads(nt)
      expect_equal(unclass(mcd_call("hessian", f)), mcd_dense_hessian(f),
                   tolerance = 1e-12, label = sprintf("[%s, %d threads]", f$name, nt))
    }
  }
})

test_that("the Hessian is exactly symmetric, linear in the weights and permuted with the labels", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (f in mcd_cells()) {
    H <- mcd_call("hessian", f)
    expect_identical(H, t(H), label = f$name)
    expect_identical(mcd_call("hessian", f, weights = 2 * f$weights), 2 * H,
                     label = f$name)
  }
  # Relabeling the non-reference alternatives permutes the constants' block
  f <- mcd_fixture("relabel", 2011, dup = TRUE)
  set.seed(2111)
  perm <- c(1L, 1L + sample.int(f$J - 1L))  # code c becomes perm[c]
  g <- f
  g$alt_idx <- perm[f$alt_idx]
  g$theta[f$K + perm[-1L] - 1L] <- f$theta[f$K + seq_len(f$J - 1L)]
  idx <- c(seq_len(f$K), f$K + perm[-1L] - 1L)
  expect_identical(mcd_call("hessian", g)[idx, idx], mcd_call("hessian", f))
})

test_that("non-finite means and weights spread over the Hessian as the dense products did", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f0 <- mcd_fixture("non-finite", 2014)
  row_off <- c(0L, cumsum(f0$M))
  i <- 3L
  r_other <- row_off[i] + setdiff(seq_len(f0$M[i]), f0$choice_idx[i])[1L]
  cells <- list(
    # utilities that overflow: every probability NaN, so every mean
    overflow = function(f) { f$X[r_other, 1L] <- 1e308; f$theta[1L] <- 2; f },
    # a missing covariate: missing probabilities
    na_x = function(f) { f$X[r_other, 1L] <- NA_real_; f },
    # p x x' past the largest double: an infinite entry, finite means
    huge_x = function(f) { f$X[r_other, 1L] <- 1e200; f$theta[1L] <- 1e-300; f },
    inf_weight = function(f) { f$weights[i] <- Inf; f },
    nan_weight = function(f) { f$weights[i] <- NaN; f },
    na_weight = function(f) { f$weights[i] <- NA_real_; f })
  for (nm in names(cells)) {
    f <- cells[[nm]](f0)
    H <- mcd_call("hessian", f)
    expect_true(any(!is.finite(H)), label = paste(nm, "has non-finite entries"))
    mcd_same_pattern(H, mcd_dense_hessian(f), paste(nm, "Hessian"))
  }
})

test_that("the compact BHHH is the weighted cross-product of the scores", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in mcd_cells()) {
    for (nt in 1:2) {
      set_omp_threads(nt)
      S <- mcd_call("scores", f)
      B <- mcd_call("bhhh", f)
      expect_identical(B, t(B), label = paste(f$name, "symmetric"))
      expect_equal(unclass(B), crossprod(sqrt(f$weights) * S), tolerance = 1e-12,
                   label = sprintf("[%s, %d threads]", f$name, nt))
    }
  }
})

test_that("at one thread the compact BHHH is the dense update bit for bit", {
  # Each entry is the product of the two scores, times the weight, added in
  # situation order, rounded at each step as R rounds the reference.
  skip_on_cran()  # a compiler that fuses the dense update's terms differs
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (f in mcd_cells()) {
    S <- mcd_call("scores", f)
    for (w in list(f$weights, f$weights^2)) {
      expect_identical(unclass(mcd_call("bhhh", f, weights = w)),
                       mcd_bhhh_ref(S, w), label = f$name)
    }
  }
})

test_that("relabeling the alternatives permutes the BHHH's constants' block", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f <- mcd_fixture("relabel", 2011, dup = TRUE)
  set.seed(2111)
  perm <- c(1L, 1L + sample.int(f$J - 1L))  # code c becomes perm[c]
  g <- f
  g$alt_idx <- perm[f$alt_idx]
  # rows within a situation stay in their order; the constants follow
  g$theta[f$K + perm[-1L] - 1L] <- f$theta[f$K + seq_len(f$J - 1L)]
  idx <- c(seq_len(f$K), f$K + perm[-1L] - 1L)
  expect_identical(mcd_call("bhhh", g)[idx, idx], mcd_call("bhhh", f))
})

test_that("non-finite scores and weights spread over the BHHH as the dense update did", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f0 <- mcd_fixture("non-finite", 2014)
  row_off <- c(0L, cumsum(f0$M))
  i <- 3L
  r_ch <- row_off[i] + f0$choice_idx[i]
  r_other <- row_off[i] + setdiff(seq_len(f0$M[i]), f0$choice_idx[i])[1L]
  cells <- list(
    # utilities that overflow: every probability, and the score, NaN
    overflow = function(f) { f$X[r_other, 1L] <- 1e308; f$theta[1L] <- 2; f },
    # a missing covariate: missing probabilities
    na_x = function(f) { f$X[r_other, 1L] <- NA_real_; f },
    # finite utilities, an infinite score: covariates of -1e308 (chosen) and
    # 1e308 times a coefficient of 1e-306
    infinite_score = function(f) {
      f$X[r_ch, 1L] <- -1e308; f$X[r_other, 1L] <- 1e308; f$theta[1L] <- 1e-306; f
    },
    inf_weight = function(f) { f$weights[i] <- Inf; f },
    nan_weight = function(f) { f$weights[i] <- NaN; f },
    na_weight = function(f) { f$weights[i] <- NA_real_; f })
  for (nm in names(cells)) {
    f <- cells[[nm]](f0)
    S <- mcd_call("scores", f)
    if (nm == "infinite_score") {
      expect_true(is.infinite(S[i, 1L]) && all(is.finite(S[i, -1L])))
    }
    B <- mcd_call("bhhh", f)
    expect_true(any(!is.finite(B)), label = paste(nm, "has non-finite entries"))
    mcd_same_pattern(B, mcd_bhhh_ref(S, f$weights), paste(nm, "BHHH"))
  }
})

# --- Cluster meat ---------------------------------------------------------

# The cluster meat kernel with labels coded as the R route codes them, and the
# R route's meat, crossprod(rowsum(w * S, cluster)).
mcd_meat <- function(f, cl, ...) {
  lab <- as.character(cl)
  mcd_call("meat", f, cluster = match(lab, unique(lab)), ...)
}
mcd_meat_ref <- function(f, cl) {
  crossprod(rowsum(f$weights * mcd_call("scores", f), group = as.character(cl)))
}

test_that("the cluster meat equals crossprod(rowsum(w * S, cluster))", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in mcd_cells()) {
    set.seed(2100)
    labs <- list(singletons = seq_len(f$N), one = rep(1L, f$N),
                 coarse = sample(1:5, f$N, replace = TRUE),
                 few = rep(c("b", "a"), length.out = f$N),
                 shuffled = sample(seq_len(f$N)) * 10.5)
    for (nm in names(labs)) {
      ref <- mcd_meat_ref(f, labs[[nm]])
      for (nt in 1:2) {
        set_omp_threads(nt)
        M <- mcd_meat(f, labs[[nm]])
        expect_identical(M, t(M), label = paste(f$name, nm, "symmetric"))
        expect_equal(unclass(M), ref, tolerance = 1e-12,
                     label = sprintf("[%s, %s, %d threads]", f$name, nm, nt))
      }
    }
  }
})

test_that("unused cluster codes and their numbering do not change the meat", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f <- mcd_fixture("codes", 2012)
  set.seed(2101)
  g <- sample(1:7, f$N, replace = TRUE)
  M1 <- mcd_call("meat", f, cluster = g)
  # the same grouping with unused codes between: the same sums in the same order
  expect_identical(mcd_call("meat", f, cluster = c(1L, 4L, 9L, 12L, 20L, 21L, 40L)[g]), M1)
  # codes in the reverse order: the clusters' updates in another order
  expect_equal(mcd_call("meat", f, cluster = 8L - g), M1, tolerance = 1e-13)
})

test_that("the cluster meat checks its codes after the other inputs", {
  f <- mcd_fixture("codes", 2013)
  expect_error(mcd_call("meat", f, cluster = rep(1L, f$N - 1L)),
               "cluster length \\([0-9]+\\) does not match the number of choice situations")
  bad <- rep(1L, f$N); bad[3L] <- 0L
  expect_error(mcd_call("meat", f, cluster = bad),
               "cluster codes must be positive integers \\(situation 3 has 0\\)")
  bad[3L] <- NA_integer_
  expect_error(mcd_call("meat", f, cluster = bad), "situation 3 has NA")
  g <- f; g$choice_idx[2L] <- 99L
  expect_error(mcd_call("meat", g, cluster = bad), "Invalid chosen alternative index")
})

test_that("non-finite scores, weights and sums spread over the meat as crossprod() did", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f0 <- mcd_fixture("non-finite", 2014)
  row_off <- c(0L, cumsum(f0$M))
  other <- function(i) row_off[i] + setdiff(seq_len(f0$M[i]), f0$choice_idx[i])[1L]
  chosen <- function(i) row_off[i] + f0$choice_idx[i]
  cl <- rep(1:6, length.out = f0$N)  # situations 3 and 9 in cluster 3
  cells <- list(
    # utilities that overflow: a NaN score
    overflow = function(f) { f$X[other(3L), 1L] <- 1e308; f$theta[1L] <- 2; f },
    # finite utilities, an infinite score
    infinite_score = function(f) {
      f$X[chosen(3L), 1L] <- -1e308; f$X[other(3L), 1L] <- 1e308
      f$theta[1L] <- 1e-306; f
    },
    # two finite weighted scores near -1e308 in one cluster: their sum overflows
    overflowing_sum = function(f) {
      for (i in c(3L, 9L)) { f$X[chosen(i), 1L] <- 0; f$X[other(i), 1L] <- 1e308 }
      f$theta[1L] <- 1e-306; f$weights[c(3L, 9L)] <- 1; f
    },
    inf_weight = function(f) { f$weights[3L] <- Inf; f },
    na_weight = function(f) { f$weights[3L] <- NA_real_; f })
  for (nm in names(cells)) {
    f <- cells[[nm]](f0)
    if (nm == "overflowing_sum") {
      S <- mcd_call("scores", f)
      expect_true(all(is.finite(S)) && !is.finite(sum(S[c(3L, 9L), 1L])))
    }
    M <- mcd_meat(f, cl)
    expect_true(any(!is.finite(M)), label = paste(nm, "has non-finite entries"))
    mcd_same_pattern(M, mcd_meat_ref(f, cl), paste(nm, "meat"))
  }
})

test_that("shared results equal the per-thread accumulators' sums", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in mcd_cells()[c(1L, 2L, 3L, 6L, 8L, 10L)]) {
    codes <- (seq_len(f$N) - 1L) %/% 3L + 1L
    for (k in c("hessian", "bhhh", "meat")) {
      what <- paste(f$name, k)
      kcall <- function(...) {
        if (k == "meat") mcd_call(k, f, cluster = codes, ...) else mcd_call(k, f, ...)
      }
      set_omp_threads(2L)
      shared <- kcall(acc_bytes = 0)
      expect_identical(shared, t(shared), label = paste(what, "symmetric"))
      expect_equal(unclass(shared), unclass(kcall()), tolerance = 1e-12,
                   label = paste(what, "shared vs per-thread"))
      # one thread adds into the result itself, whatever the budget
      set_omp_threads(1L)
      expect_identical(kcall(acc_bytes = 0), kcall(), label = paste(what, "one thread"))
    }
  }
  # Non-finite entries spread over a shared result as over per-thread ones
  set_omp_threads(2L)
  f0 <- mcd_fixture("non-finite, shared", 2015)
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
  codes <- rep(1:5, length.out = f0$N)
  for (nm in names(cells)) {
    for (k in c("hessian", "bhhh", "meat")) {
      extra <- if (k == "meat") list(cluster = codes) else list()
      mcd_same_pattern(do.call(mcd_call, c(list(k, cells[[nm]], acc_bytes = 0), extra)),
                       unclass(do.call(mcd_call, c(list(k, cells[[nm]]), extra))),
                       paste(nm, k, "shared"))
    }
  }
})
