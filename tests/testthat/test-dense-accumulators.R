# The elasticity and diversion-ratio kernels of every model add each choice
# situation's terms, by alternative, into a J x J sum (DenseAcc): the result
# itself with one thread, the threads' copies while they fit in acc_bytes,
# and above that the shared result, added into atomically (with row and
# column 0 kept per thread when they are the outside option's). They check
# the count of alternatives and their memory before allocating, and the
# methods label the matrix in place.

# N situations over J alternatives: the first offers every one, the others
# m_range of them, every fifth listing one twice (with `unoffered`, the last
# alternative is offered nowhere); K covariates, nests of three.
dac_fixture <- function(seed, ioo = FALSE, N = 40L, J = 9L, m_range = c(2L, 6L),
                        K = 2L, unoffered = FALSE) {
  set.seed(seed)
  pool <- if (unoffered) seq_len(J - 1L) else seq_len(J)
  alts <- lapply(seq_len(N), function(i) {
    a <- sort(sample(pool, sample(m_range[1L]:m_range[2L], 1L)))
    if (i %% 5L == 0L) a <- sort(c(a, a[1L]))
    a
  })
  alts[[1L]] <- pool
  M <- lengths(alts)
  rows <- sum(M)
  nest_idx <- as.integer((seq_len(J) - 1L) %/% 3L + 1L)
  K_l <- sum(tabulate(nest_idx) > 1L)
  n_delta <- if (ioo) J else J - 1L
  list(X = matrix(stats::rnorm(rows * K), rows, K),
       W = matrix(stats::rnorm(rows * 2L), rows, 2L),
       alt_idx = as.integer(unlist(alts)), M = as.integer(M),
       choice_idx = vapply(M, function(m) sample.int(m, 1L), integer(1)),
       weights = stats::runif(N, 0.5, 2), nest_idx = nest_idx, ioo = ioo,
       beta = stats::rnorm(K, sd = 0.5), lambda = stats::runif(K_l, 0.5, 0.9),
       delta = stats::rnorm(n_delta, sd = 0.4), mxl = stats::rnorm(K + 4L, sd = 0.3),
       eta = get_halton_normals(5L, N, 2L), J = J, K = K)
}

# The six kernels on fixture f; ... reaches the kernel (acc_bytes).
dac_kernels <- c("mnl_elas", "mnl_dr", "nl_elas", "nl_dr", "mxl_elas", "mxl_dr",
                 "mxl_elas_rc", "mxl_dr_rc")

dac_call <- function(kernel, f, ...) {
  ns <- asNamespace("choicer")
  th_mnl <- c(f$beta, f$delta)
  th_nl <- c(f$beta, f$lambda, f$delta)
  # beta (K), mu (2), the Cholesky diagonal (2), delta
  th_mxl <- c(f$mxl, f$delta)
  rc <- endsWith(kernel, "_rc")
  switch(sub("_rc$", "", kernel),
    mnl_elas = ns$mnl_elasticities_parallel(th_mnl, f$X, f$alt_idx, f$choice_idx,
                                            f$M, f$weights, 2L, TRUE, f$ioo, ...),
    mnl_dr = ns$mnl_diversion_ratios_parallel(th_mnl, f$X, f$alt_idx, f$M,
                                              f$weights, TRUE, f$ioo, ...),
    nl_elas = ns$nl_elasticities_parallel(th_nl, f$X, f$alt_idx, f$choice_idx,
                                          f$nest_idx, f$M, f$weights, 2L, TRUE,
                                          f$ioo, ...),
    nl_dr = ns$nl_diversion_ratios_parallel(th_nl, f$X, f$alt_idx, f$nest_idx,
                                            f$M, f$weights, TRUE, f$ioo, ...),
    mxl_elas = ns$mxl_elasticities_parallel(th_mxl, f$X, f$W, f$alt_idx,
                                            f$choice_idx, f$M, f$weights, f$eta,
                                            c(0L, 0L), if (rc) 1L else 2L, rc,
                                            FALSE, TRUE, TRUE, f$ioo, ...),
    mxl_dr = ns$mxl_diversion_ratios_parallel(th_mxl, f$X, f$W, f$alt_idx, f$M,
                                              f$weights, f$eta, c(0L, 0L),
                                              if (rc) 1L else 2L, rc, FALSE, TRUE,
                                              TRUE, f$ioo, ...))
}

test_that("at one thread the J x J sums do not depend on acc_bytes", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (ioo in c(FALSE, TRUE)) {
    f <- dac_fixture(5101 + ioo, ioo = ioo)
    for (k in dac_kernels) {
      expect_identical(dac_call(k, f, acc_bytes = 0), dac_call(k, f),
                       label = paste(k, "outside option:", ioo))
    }
  }
})

test_that("shared, per-thread and one-thread sums agree", {
  # Two threads add into their copies (the default budget) or into the
  # shared result (acc_bytes = 0), with row and column 0 per thread when they
  # are the outside option's; the sums differ from one thread's in the order
  # of their additions only.
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_omp_threads(2L), add = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    f <- dac_fixture(5111 + ioo, ioo = ioo, N = 120L)
    for (k in dac_kernels) {
      what <- paste(k, "outside option:", ioo)
      set_omp_threads(1L)
      one <- dac_call(k, f)
      set_omp_threads(2L)
      expect_equal(dac_call(k, f), one, tolerance = 1e-12,
                   label = paste(what, "copies"))
      expect_equal(dac_call(k, f, acc_bytes = 0), one, tolerance = 1e-12,
                   label = paste(what, "shared"))
    }
  }
})

test_that("one thread's J x J sums are +0 where a fused term rounds to -0", {
  # Alternative 2's probability is subnormal, and its weighted
  # cross-elasticity's exact value lies below the smallest subnormal: a fused
  # multiply-add from +0 rounds it to -0. The threads' sums were always added
  # into a zero matrix, which made it +0; with one thread too.
  on.exit(set_omp_threads(2L), add = TRUE)
  ns <- asNamespace("choicer")
  for (nt in 1:2) {
    set_omp_threads(nt)
    x <- ns$mnl_elasticities_parallel(c(1, -746), matrix(c(0, 1), 2L), c(1L, 2L),
                                      1L, 2L, 0.25, 1L, TRUE, FALSE)
    expect_identical(1 / x[1L, 2L], Inf, label = sprintf("%d thread(s)", nt))
  }
})

test_that("MNL diversion ratios match their formula, with zero diagonal and unoffered columns", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    f <- dac_fixture(5121 + ioo, ioo = ioo, unoffered = TRUE)
    n <- f$J + ioo
    num <- matrix(0, n, n)
    den <- numeric(n)
    v <- drop(f$X %*% f$beta) +
      if (ioo) f$delta[f$alt_idx] else c(0, f$delta)[f$alt_idx]
    row_off <- c(0L, cumsum(f$M))
    for (i in seq_along(f$M)) {
      r <- (row_off[i] + 1L):row_off[i + 1L]
      u <- if (ioo) c(0, v[r]) else v[r]
      p <- exp(u - max(u)) / sum(exp(u - max(u)))
      g <- if (ioo) c(1L, f$alt_idx[r] + 1L) else f$alt_idx[r]
      for (j in seq_along(g)) {
        den[g[j]] <- den[g[j]] + f$weights[i] * p[j] * (1 - p[j])
        for (k in seq_along(g)[-j]) {
          num[g[k], g[j]] <- num[g[k], g[j]] + f$weights[i] * p[j] * p[k]
        }
      }
    }
    dr <- sweep(num, 2L, den, "/")
    dr[, den <= 1e-15] <- 0
    diag(dr) <- 0
    for (nt in 1:2) {
      set_omp_threads(nt)
      for (acc in c(2^31, 0)) {
        got <- dac_call("mnl_dr", f, acc_bytes = acc)
        expect_equal(unclass(got), dr, tolerance = 1e-12,
                     label = sprintf("outside option %s, %d threads, acc %g", ioo, nt, acc))
        expect_true(all(diag(got) == 0) && all(got[, n] == 0) &&
                      !any(1 / got[, n] < 0), label = "zero diagonal and column, +0")
      }
    }
  }
})

test_that("the J x J kernels check the count of alternatives, then their memory", {
  ns <- asNamespace("choicer")
  # alternative code 2^31 - 1 with an outside option: J + 1 does not fit an int
  X <- matrix(c(0.5, -0.2, 0.1), 3L, 1L)
  alt <- c(1L, .Machine$integer.max, 2L)
  for (call in list(
    function() ns$mnl_elasticities_parallel(0.3, X, alt, 1L, 3L, 1, 1L, FALSE, TRUE),
    function() ns$mnl_diversion_ratios_parallel(0.3, X, alt, 3L, 1, FALSE, TRUE))) {
    expect_error(call(), "more alternatives than an int can count")
  }

  f <- dac_fixture(5131, ioo = TRUE)
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  for (k in dac_kernels) {
    expect_error(dac_call(k, f),
                 "matrix needs about [0-9.]+ kB even with one thread .* more than the 0 kB of options\\(choicer.max_memory\\)",
                 label = k)
  }
  # the input checks come first
  expect_error(ns$mnl_elasticities_parallel(c(f$beta, f$delta), f$X, f$alt_idx,
                                            f$choice_idx, f$M, f$weights, 7L,
                                            TRUE, TRUE),
               "elast_var_idx is out of bounds")
})

test_that("the J x J kernels name the threads that fit", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_omp_threads(2L), add = TRUE)
  f <- dac_fixture(5141, J = 60L, N = 80L)
  for (k in c("mnl_elas", "nl_dr", "mxl_elas")) {
    # one thread's estimate, then a limit between it and two threads' (the
    # copies of the J x J sums do not fit twice)
    set_omp_threads(1L)
    old <- options(choicer.max_memory = 1)
    msg <- tryCatch(dac_call(k, f), error = conditionMessage)
    options(old)
    m <- regmatches(msg, regexec("needs about ([0-9.]+) (kB|MB)", msg))[[1]]
    one <- as.numeric(m[2L]) * c(kB = 1e3, MB = 1e6)[[m[3L]]]
    set_omp_threads(2L)
    old <- options(choicer.max_memory = one * 1.02)
    expect_error(dac_call(k, f), "with 2 threads .* Use set_num_threads\\(1\\)",
                 label = k)
    options(old)
  }
})

test_that("elasticities() and diversion_ratios() label their matrices with the alternatives", {
  set.seed(5151)
  N <- 60L; J <- 4L
  dt <- data.table(id = rep(seq_len(N), each = J),
                   alt = factor(rep(c("bus", "car", "rail", "air"), N),
                                levels = c("bus", "car", "rail", "air")),
                   x1 = rnorm(N * J), x2 = rnorm(N * J))
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt[, nest := ifelse(alt %in% c("bus", "rail"), "ground", "other")]
  fits <- suppressMessages(list(
    mnl = run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2")),
    nl = run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")))
  for (m in names(fits)) {
    labs <- as.character(fits[[m]]$alt_mapping[[2]])
    for (res in list(elasticities(fits[[m]], "x1"), diversion_ratios(fits[[m]]))) {
      expect_identical(dimnames(res), list(labs, labs), label = m)
    }
  }
})
