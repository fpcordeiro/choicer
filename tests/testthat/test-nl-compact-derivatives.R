# The nested logit Hessian, BHHH and cluster-meat kernels on designs with
# many alternatives and sparse choice sets. A situation's score and Hessian
# are nonzero only on the coefficients, the dissimilarities of the nests it
# offers and the constants of the alternatives it offers, so the kernels form
# them in that compact block and add it into the n x n result (the meat a
# cluster's summed score, over the parameters its situations touch). These
# designs (30 alternatives in eight nests of three and six singletons, choice
# sets of 3-8, dissimilarities on both sides of one) exercise the slot maps:
# nests offered in part or not at all, singleton nests (their dissimilarity
# fixed at one), many small nests, one large nest, a situation listing an
# alternative twice, the outside option, no constants, and twenty
# coefficients (Armadillo keeps vectors of more than 16 on the heap).

# N situations, each offering `m_range` distinct alternatives drawn from
# 1..J, plus one repeated row when `dup`; weights vary across situations.
ncd_fixture <- function(name, seed, N = 60L, J = 30L, m_range = c(3L, 8L),
                        K = 2L, nest_idx = as.integer(c(rep(1:8, each = 3L), 9:14)),
                        include_outside_option = FALSE, use_asc = TRUE,
                        dup = FALSE, shuffle = FALSE, runs = FALSE) {
  set.seed(seed)
  alts <- lapply(seq_len(N), function(t) {
    m <- sample(m_range[1L]:m_range[2L], 1L)
    a <- sort(sample.int(J, min(m, J)))
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    if (shuffle) a <- a[sample.int(length(a))]
    a
  })
  # `runs`: most situations repeat the previous one's rows (the kernels keep
  # a slot map while the codes repeat)
  if (runs) {
    for (t in seq_along(alts)[-1L]) {
      if (stats::runif(1L) < 0.7) alts[[t]] <- alts[[t - 1L]]
    }
  }
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
                nest_idx = rep(1:20, each = 3L)),
    ncd_fixture("twenty coefficients", 3010, K = 20L),
    # Products X_i' v in the shapes Armadillo dispatches differently (more
    # than 16 covariates and choices, a 4 x 4 X_i), and many nests of which a
    # situation offers a few (the loops over the nests it offers)
    ncd_fixture("wide covariates and sets", 3021, J = 40L, m_range = c(10L, 30L),
                K = 17L, N = 30L, nest_idx = c(rep(1:10, each = 3L), 11:20)),
    ncd_fixture("four alternatives, four covariates", 3022, J = 4L,
                m_range = c(4L, 4L), K = 4L, nest_idx = c(1L, 1L, 2L, 2L)),
    ncd_fixture("many nests, few offered", 3023, J = 120L, m_range = c(2L, 6L),
                nest_idx = c(rep(1:50, each = 2L), 51:70), shuffle = TRUE),
    ncd_fixture("runs of equal sets, unsorted, repeated", 3024, runs = TRUE,
                dup = TRUE, shuffle = TRUE),
    ncd_fixture("runs of equal sets, outside option", 3025, runs = TRUE,
                include_outside_option = TRUE)
  )
}

# The free constant of each row (its position in theta), NA where none.
ncd_free <- function(f) {
  if (!f$use_asc) return(rep(NA_integer_, length(f$alt_idx)))
  j <- if (f$include_outside_option) f$alt_idx else f$alt_idx - 1L
  ifelse(j >= 1L, f$K + f$K_l + j, NA_integer_)
}

# Each nest's dissimilarity at theta (1 for a singleton).
ncd_lambda <- function(f, theta = f$theta) {
  n_k <- tabulate(f$nest_idx)
  lam <- rep(1, length(n_k))
  lam[n_k > 1L] <- theta[f$K + seq_len(f$K_l)]
  lam
}

# NL probabilities of one situation's rows (nests: their nests) and of the
# outside option, in R.
ncd_probs <- function(v, nests, lambda, outside) {
  ks <- unique(nests)
  iv <- vapply(ks, function(k) {
    u <- v[nests == k] / lambda[k]
    max(u) + log(sum(exp(u - max(u))))
  }, numeric(1L))
  a <- lambda[ks] * iv
  top <- max(c(if (outside) 0, a))
  den <- (if (outside) exp(-top) else 0) + sum(exp(a - top))
  p_k <- exp(a - top) / den
  pos <- match(nests, ks)
  list(p = exp(v / lambda[nests] - iv[pos]) * p_k[pos],
       p0 = if (outside) exp(-top) / den else 0)
}

ncd_call <- function(kernel, f, theta = f$theta, weights = f$weights, ...) {
  ns <- asNamespace("choicer")
  a <- list(theta = theta, X = f$X, alt_idx = f$alt_idx,
            choice_idx = f$choice_idx, nest_idx = f$nest_idx, M = f$M,
            weights = weights, use_asc = f$use_asc,
            include_outside_option = f$include_outside_option)
  switch(kernel,
    gradient = do.call(ns$nl_loglik_gradient_parallel, c(a, list(...))),
    hessian = do.call(ns$nl_loglik_hessian_parallel, c(a, list(...))),
    bhhh = do.call(ns$nl_bhhh_parallel, c(a, list(...))),
    scores = do.call(ns$nl_scores_parallel, c(a[-7L], list(...))),
    meat = do.call(ns$nl_cluster_meat_parallel, c(a, list(...))))
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

# The parameters situation i touches (its compact block): the coefficients,
# the dissimilarities of the non-singleton nests it offers, the constants of
# its alternatives.
ncd_block <- function(f, i) {
  row_off <- c(0L, cumsum(f$M))
  alts <- f$alt_idx[(row_off[i] + 1L):row_off[i + 1L]]
  multi <- which(tabulate(f$nest_idx) > 1L)
  lam <- f$K + match(intersect(multi, f$nest_idx[alts]), multi)
  j <- if (!f$use_asc) integer(0) else if (f$include_outside_option) alts else alts[alts > 1L] - 1L
  sort(unique(c(seq_len(f$K), lam, f$K + f$K_l + j)))
}

# f without situation i
ncd_drop <- function(f, i) {
  row_off <- c(0L, cumsum(f$M))
  rows <- (row_off[i] + 1L):row_off[i + 1L]
  f$X <- f$X[-rows, , drop = FALSE]
  f$alt_idx <- f$alt_idx[-rows]
  f$choice_idx <- f$choice_idx[-i]
  f$M <- f$M[-i]
  f$weights <- f$weights[-i]
  f$N <- f$N - 1L
  f
}

# f with situation i only
ncd_keep <- function(f, i) {
  row_off <- c(0L, cumsum(f$M))
  rows <- (row_off[i] + 1L):row_off[i + 1L]
  f$X <- f$X[rows, , drop = FALSE]
  f$alt_idx <- f$alt_idx[rows]
  f$choice_idx <- f$choice_idx[i]
  f$M <- f$M[i]
  f$weights <- f$weights[i]
  f$N <- 1L
  f
}

# The budgets (acc_bytes) of two threads sharing the result: 0 keeps the
# coefficients' rows private, the bytes of the threads' copies of the
# coefficients' and dissimilarities' rows keep those (their packed triangles
# do not fit; without constants the two coincide).
ncd_budgets <- function(f) {
  n <- length(f$theta)
  s <- f$K + f$K_l
  c(coefficients = 0, dissimilarities = 2 * 8 * (n * s - s * (s - 1) / 2))
}

ncd_scale_dev <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))

test_that("the compact NL Hessian matches numDeriv of the gradient", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  for (f in ncd_cells()) {
    num <- numDeriv::jacobian(function(th) ncd_call("gradient", f, theta = th)$gradient,
                              f$theta)
    expect_lt(ncd_scale_dev(unclass(ncd_call("hessian", f)), num), 1e-6, label = f$name)
  }
})

test_that("the NL Hessian is exactly symmetric, linear in the weights and permuted with the labels", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (f in ncd_cells()) {
    H <- ncd_call("hessian", f)
    expect_identical(H, t(H), label = f$name)
    expect_identical(ncd_call("hessian", f, weights = 2 * f$weights), 2 * H,
                     label = f$name)
  }
  f <- ncd_fixture("relabel", 3011, dup = TRUE)
  set.seed(3111)
  perm <- c(1L, 1L + sample.int(f$J - 1L))
  g <- f
  g$alt_idx <- perm[f$alt_idx]
  g$nest_idx[perm] <- f$nest_idx
  g$theta[f$K + f$K_l + perm[-1L] - 1L] <- f$theta[f$K + f$K_l + seq_len(f$J - 1L)]
  idx <- c(seq_len(f$K + f$K_l), f$K + f$K_l + perm[-1L] - 1L)
  expect_identical(ncd_call("hessian", g)[idx, idx], ncd_call("hessian", f))
  sigma <- c(sample.int(8L), 9:14)
  h <- f
  h$nest_idx <- sigma[f$nest_idx]
  h$theta[f$K + sigma[1:8]] <- f$theta[f$K + 1:8]
  idx <- c(seq_len(f$K), f$K + sigma[1:8], f$K + f$K_l + seq_len(f$J - 1L))
  expect_equal(unclass(ncd_call("hessian", h))[idx, idx],
               unclass(ncd_call("hessian", f)), tolerance = 1e-12)
})

test_that("an NL situation whose chosen probability is not finite is left out", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f <- ncd_fixture("skip", 3016)
  row_off <- c(0L, cumsum(f$M))
  # the chosen alternative's utility at -Inf, so its log-probability
  i <- 3L
  r_chosen <- row_off[i] + f$choice_idx[i]
  g <- f
  g$X[r_chosen, 1L] <- -1e308
  g$theta[1L] <- 2
  d <- ncd_drop(f, i)
  d$theta[1L] <- 2
  H <- unclass(ncd_call("hessian", g))
  expect_true(all(is.finite(H)))
  expect_equal(H, unclass(ncd_call("hessian", d)), tolerance = 1e-14)
})

test_that("NaN probabilities in an NL situation zero the entries their terms reach", {
  # An unchosen alternative's utility at +Inf, or a missing covariate, leaves
  # the chosen probability finite (the alternative's nest drops out) and the
  # nest's other probabilities not numbers. The situation's terms that use
  # them are not finite, and every entry they reach is set to 0 at the end,
  # as the dense code set it; its other terms are added as usual.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (seed in c(3016, 3018)) {
    f <- ncd_fixture("NaN probabilities", seed)
    f$theta[1L] <- 2
    i <- 3L
    row_off <- c(0L, cumsum(f$M))
    r_other <- row_off[i] + setdiff(seq_len(f$M[i]), f$choice_idx[i])[1L]
    rest <- unclass(ncd_call("hessian", ncd_drop(f, i)))
    blk <- ncd_block(f, i)
    for (v in c(1e308, NA_real_, NaN)) {
      what <- paste(seed, "covariate", v)
      g <- f
      g$X[r_other, 1L] <- v
      H <- unclass(ncd_call("hessian", g))
      own <- unclass(ncd_call("hessian", ncd_keep(g, i)))
      expect_true(all(is.finite(H)), label = what)
      zero <- H == 0 & rest != 0
      expect_true(any(zero), label = paste(what, "zeroes entries"))
      expect_true(all(zero[-blk, ] == FALSE) && all(zero[, -blk] == FALSE),
                  label = paste(what, "only in its block"))
      expect_true(all(own[zero] == 0), label = paste(what, "where its terms are not finite"))
      expect_equal(H[!zero], (rest + own)[!zero], tolerance = 1e-12,
                   label = paste(what, "elsewhere"))
    }
  }
})

test_that("a non-finite NL weight zeroes exactly the entries its situation touches", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f <- ncd_fixture("weights", 3017)
  i <- 4L
  blk <- ncd_block(f, i)
  ref <- unclass(ncd_call("hessian", ncd_drop(f, i)))
  inside <- matrix(FALSE, length(f$theta), length(f$theta))
  inside[blk, blk] <- TRUE
  for (w_bad in c(Inf, NaN, NA_real_)) {
    g <- f
    g$weights[i] <- w_bad
    H <- unclass(ncd_call("hessian", g))
    expect_true(all(H[inside] == 0), label = paste("weight", w_bad, "block"))
    expect_equal(H[!inside], ref[!inside], tolerance = 1e-14,
                 label = paste("weight", w_bad, "elsewhere"))
  }
})

test_that("the NL Hessian's sums overflow only where the triangles' sums do", {
  # One situation, one nest of two alternatives (lambda = 2) and a covariate
  # near sqrt(DBL_MAX): a term of about 0.56 DBL_MAX, the entry about 0.09
  # DBL_MAX. Each triangle's running sum stays finite, and so does the sum
  # of their halves; the sum of the two terms would not.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  ns <- asNamespace("choicer")
  hess <- function(w) {
    ns$nl_loglik_hessian_parallel(c(1e-200, 0.1, 2, 0.05),
                                  cbind(c(sqrt(1.7e308), 0), c(0.3, -0.2)),
                                  1:2, 1L, c(1L, 1L), 2L, w)
  }
  H <- hess(1.5)
  expect_true(all(is.finite(H)) && H[1L, 1L] > 1.5e307)
  expect_identical(H, 2 * hess(0.75))
})

test_that("a shared NL Hessian equals the per-thread accumulators' sum", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in ncd_cells()[c(1L, 2L, 3L, 5L, 7L, 9L, 10L)]) {
    for (b in names(ncd_budgets(f))) {
      what <- paste(f$name, b, "private")
      set_omp_threads(2L)
      shared <- ncd_call("hessian", f, acc_bytes = ncd_budgets(f)[[b]])
      expect_identical(shared, t(shared), label = paste(what, "symmetric"))
      expect_equal(unclass(shared), unclass(ncd_call("hessian", f)), tolerance = 1e-12,
                   label = paste(what, "shared vs per-thread"))
    }
    set_omp_threads(1L)
    expect_identical(ncd_call("hessian", f, acc_bytes = 0), ncd_call("hessian", f),
                     label = paste(f$name, "one thread"))
  }
})

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

test_that("the NL gradient, scores, predictions and shares agree with R on every design", {
  # The kernels form each situation's probabilities, nest sums and products
  # X_i' v in per-thread arrays, and loop over the nests the situation
  # offers: the gradient against numDeriv of the objective (the
  # dissimilarities' entries included) and the weighted sum of the scores,
  # its OPG diagonal against the scores', the probabilities and utilities
  # against R's, the shares against the weighted sums of the probabilities.
  skip_if_not_installed("numDeriv")
  on.exit(set_omp_threads(2L), add = TRUE)
  ns <- asNamespace("choicer")
  for (f in ncd_cells()) {
    fr <- ncd_free(f)
    lam <- ncd_lambda(f)
    v_ref <- drop(f$X %*% f$theta[seq_len(f$K)]) + ifelse(is.na(fr), 0, f$theta[fr])
    row_off <- c(0L, cumsum(f$M))
    p_ref <- numeric(length(v_ref))
    shares_ref <- numeric(f$J + f$include_outside_option)
    for (i in seq_len(f$N)) {
      r <- (row_off[i] + 1L):row_off[i + 1L]
      p <- ncd_probs(v_ref[r], f$nest_idx[f$alt_idx[r]], lam, f$include_outside_option)
      p_ref[r] <- p$p
      if (f$include_outside_option) {
        shares_ref[1L] <- shares_ref[1L] + f$weights[i] * p$p0
      }
      a <- f$alt_idx[r] + f$include_outside_option
      for (k in seq_along(r)) shares_ref[a[k]] <- shares_ref[a[k]] + f$weights[i] * p$p[k]
    }
    shares_ref <- shares_ref / sum(f$weights)
    for (nt in 1:2) {
      set_omp_threads(nt)
      what <- sprintf("[%s, %d threads]", f$name, nt)
      g <- ncd_call("gradient", f, opg_diag = TRUE)
      S <- ncd_call("scores", f)
      expect_equal(drop(g$gradient), -colSums(f$weights * S), tolerance = 1e-12,
                   label = what)
      expect_equal(drop(g$opg_diag), colSums(f$weights * S^2), tolerance = 1e-12,
                   label = what)
      pr <- ns$nl_predict(f$theta, f$X, f$alt_idx, f$M, f$nest_idx, f$use_asc,
                          f$include_outside_option)
      expect_equal(drop(pr$utility), v_ref, tolerance = 1e-14,
                   label = paste(what, "utilities"))
      expect_equal(drop(pr$choice_prob), p_ref, tolerance = 1e-12,
                   label = paste(what, "probabilities"))
      # one share per code up to the largest offered (the outside option first)
      sh <- drop(ns$nl_predict_shares(f$theta, f$X, f$alt_idx, f$M, f$weights,
                                      f$nest_idx, f$use_asc, f$include_outside_option))
      expect_equal(sh, shares_ref[seq_along(sh)], tolerance = 1e-12,
                   label = paste(what, "shares"))
      expect_true(all(shares_ref[-seq_along(sh)] == 0), label = paste(what, "unoffered"))
    }
    set_omp_threads(1L)
    num <- numDeriv::grad(function(th) ncd_call("gradient", f, theta = th)$objective,
                          f$theta)
    expect_lt(ncd_scale_dev(drop(ncd_call("gradient", f)$gradient), num), 1e-6,
              label = f$name)
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

test_that("the NL score matrix refuses to exceed max_bytes", {
  f <- ncd_fixture("guard", 3018, N = 5L)
  expect_error(ncd_call("scores", f, max_bytes = 100),
               "The score matrix would take 0.00 GiB (5 choice situations x",
               fixed = TRUE)
  expect_equal(dim(ncd_call("scores", f)), c(5L, length(f$theta)))
})

# --- Cluster meat ---------------------------------------------------------

# The cluster meat kernel with labels coded as the R route codes them, and the
# R route's meat, crossprod(rowsum(w * S, cluster)).
ncd_meat <- function(f, cl, ...) {
  lab <- as.character(cl)
  ncd_call("meat", f, cluster = match(lab, unique(lab)), ...)
}
ncd_meat_ref <- function(f, cl) {
  crossprod(rowsum(f$weights * ncd_call("scores", f), group = as.character(cl)))
}

test_that("the NL cluster meat equals crossprod(rowsum(w * S, cluster))", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in ncd_cells()) {
    set.seed(3100)
    labs <- list(singletons = seq_len(f$N), one = rep(1L, f$N),
                 coarse = sample(1:5, f$N, replace = TRUE),
                 few = rep(c("b", "a"), length.out = f$N))
    for (nm in names(labs)) {
      ref <- ncd_meat_ref(f, labs[[nm]])
      for (nt in 1:2) {
        set_omp_threads(nt)
        M <- ncd_meat(f, labs[[nm]])
        expect_identical(M, t(M), label = paste(f$name, nm, "symmetric"))
        expect_equal(unclass(M), ref, tolerance = 1e-12,
                     label = sprintf("[%s, %s, %d threads]", f$name, nm, nt))
      }
    }
  }
})

test_that("the NL cluster meat checks its codes after the other inputs and before its memory", {
  f <- ncd_fixture("codes", 3013)
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  expect_error(ncd_call("meat", f, cluster = rep(1L, f$N - 1L)),
               "cluster length \\([0-9]+\\) does not match the number of choice situations")
  bad <- rep(1L, f$N); bad[3L] <- 0L
  expect_error(ncd_call("meat", f, cluster = bad),
               "cluster codes must be positive integers \\(situation 3 has 0\\)")
  bad[3L] <- NA_integer_
  expect_error(ncd_call("meat", f, cluster = bad), "situation 3 has NA")
  g <- f; g$choice_idx[2L] <- 99L
  expect_error(ncd_call("meat", g, cluster = bad), "Invalid chosen alternative index")
  expect_error(ncd_call("meat", f, cluster = rep(1L, f$N)), "needs about")
})

test_that("non-finite scores, weights and sums spread over the NL meat as crossprod() did", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  f0 <- ncd_fixture("non-finite", 3014)
  row_off <- c(0L, cumsum(f0$M))
  other <- function(i) row_off[i] + setdiff(seq_len(f0$M[i]), f0$choice_idx[i])[1L]
  chosen <- function(i) row_off[i] + f0$choice_idx[i]
  cl <- rep(1:6, length.out = f0$N)  # situations 3 and 9 in cluster 3
  cells <- list(
    # utilities that overflow: NaN scores
    overflow = function(f) { f$X[other(3L), 1L] <- 1e308; f$theta[1L] <- 2; f },
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
      S <- ncd_call("scores", f)
      expect_true(all(is.finite(S)) && !is.finite(sum(S[c(3L, 9L), 1L])))
    }
    M <- ncd_meat(f, cl)
    expect_true(any(!is.finite(M)), label = paste(nm, "has non-finite entries"))
    ncd_same_pattern(M, ncd_meat_ref(f, cl), paste(nm, "meat"))
  }
})

test_that("shared NL BHHH and meat results equal the per-thread accumulators' sums", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(set_omp_threads(2L), add = TRUE)
  for (f in ncd_cells()[c(1L, 2L, 3L, 5L, 7L, 9L, 10L)]) {
    codes <- (seq_len(f$N) - 1L) %/% 3L + 1L
    for (k in c("bhhh", "meat")) {
      what <- paste(f$name, k)
      kcall <- function(...) {
        if (k == "meat") ncd_call(k, f, cluster = codes, ...) else ncd_call(k, f, ...)
      }
      set_omp_threads(2L)
      for (b in names(ncd_budgets(f))) {
        shared <- kcall(acc_bytes = ncd_budgets(f)[[b]])
        expect_identical(shared, t(shared), label = paste(what, b, "symmetric"))
        expect_equal(unclass(shared), unclass(kcall()), tolerance = 1e-12,
                     label = paste(what, b, "shared vs per-thread"))
      }
      # one thread adds into the result itself, whatever the budget
      set_omp_threads(1L)
      expect_identical(kcall(acc_bytes = 0), kcall(), label = paste(what, "one thread"))
    }
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
  codes <- rep(1:5, length.out = f0$N)
  for (nm in names(cells)) {
    for (k in c("bhhh", "meat")) {
      extra <- if (k == "meat") list(cluster = codes) else list()
      ncd_same_pattern(do.call(ncd_call, c(list(k, cells[[nm]], acc_bytes = 0), extra)),
                       unclass(do.call(ncd_call, c(list(k, cells[[nm]]), extra))),
                       paste(nm, k, "shared"))
    }
  }
})

test_that("a kept NL slot map follows which nests' inclusive values are finite", {
  # Three situations with one set of rows, nests (1, 1, 2, 2): the second
  # sends nest 2's utilities to -Inf. Its log I_2 = -Inf leaves nest 2's
  # dissimilarity out of its block (a zero score there) while its other
  # scores turn NaN, so the BHHH and the meat take the dense update, whose
  # NaN reach every entry but those the excluded zero meets with a finite
  # score. The map built for the first situation, where nest 2 is finite,
  # must not serve the second, nor the second's the third.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  set.seed(3201)
  X <- matrix(stats::rnorm(12L * 2L), 12L, 2L)
  X[7:8, 1L] <- -1e308  # theta[1] = 2: V = -Inf in nest 2 of situation 2
  f <- list(name = "finiteness", X = X, alt_idx = rep(1:4, 3L),
            choice_idx = c(1L, 2L, 1L), nest_idx = c(1L, 1L, 2L, 2L),
            M = rep(4L, 3L), weights = c(0.7, 1.3, 1.1),
            theta = c(2, -0.4, 0.6, 0.8, 0.3, -0.2, 0.5),
            use_asc = TRUE, include_outside_option = FALSE, J = 4L, K = 2L,
            K_l = 2L, N = 3L)
  S <- ncd_call("scores", f)
  expect_identical(S[2L, 4L], 0)  # nest 2's dissimilarity: out of the block
  expect_true(any(is.nan(S[2L, ])) && all(is.finite(S[-2L, ])))
  B <- ncd_call("bhhh", f)
  ncd_same_pattern(B, ncd_bhhh_ref(S, f$weights), "BHHH")
  expect_true(all(is.finite(B[4L, c(3L, 5L)])), label = "excluded score meets finite ones")
  cl <- c(1L, 1L, 2L)
  ncd_same_pattern(ncd_meat(f, cl), ncd_meat_ref(f, cl), "meat")
})
