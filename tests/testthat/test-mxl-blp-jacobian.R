# The share Jacobian on the estimation draws (mxl_blp_jacobian_parallel()),
# for run_mxlogit(concentrate_out_asc = "blp"), against the brute-force
# reference of helper-blp-reference.R; and spd_solve(), its solves.

# The panel kernels' cells (helper-mxl-panel.R), with constants
blpj_cells <- function() {
  lapply(mxlp_cell_configs(), function(cfg) {
    cfg$use_asc <- TRUE
    do.call(mxlp_fixture, c(list(name = mxlp_cell_name(cfg)), cfg))
  })
}

test_that("the share Jacobian is the derivative of the predicted counts", {
  skip_if_not_installed("numDeriv")
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in c(blpj_cells(), list(mxlp_cross_section(blpj_cells()[[2L]])))) {
    n <- length(fx$theta)
    K_c <- n - length(blpr_free(fx))
    fr <- blpr_free(fx)
    r <- blpr_call(fx, variance = TRUE)
    # the objective and the gradient are the gradient kernel's, bit for bit
    g <- choicer:::mxl_loglik_gradient_parallel(
      fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
      fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
      Ti = fx[["Ti"]])
    expect_identical(r$objective, g$objective, label = fx$name)
    expect_identical(r$gradient, g$gradient, label = fx$name)
    # the counts, every inside alternative's
    ref <- blpr_counts(fx$theta, fx)
    expect_equal(as.numeric(r$counts), ref, tolerance = 1e-13, label = fx$name)
    # the Jacobian of the free constants' counts: C_c, C_d; a zero
    # continuous block
    jn <- numDeriv::jacobian(function(th) blpr_counts(th, fx)[fr], fx$theta)
    d_idx <- K_c + seq_along(fr)
    expect_equal(r$jacobian[d_idx, ], jn, tolerance = 1e-7, label = fx$name)
    expect_equal(r$jacobian[seq_len(K_c), d_idx], t(jn[, seq_len(K_c)]),
                 tolerance = 1e-7, label = fx$name)
    expect_true(all(r$jacobian[seq_len(K_c), seq_len(K_c)] == 0))
    expect_identical(r$jacobian, t(r$jacobian))
    # K, the units' across-draw covariance of their predicted counts
    expect_equal(r$K, blpr_K(fx$theta, fx), tolerance = 1e-12, label = fx$name)
    # without the variance mode: the same matrices to rounding (the
    # probabilities read from the score's residuals), no K
    r0 <- blpr_call(fx)
    expect_null(r0$K)
    expect_equal(r0$jacobian, r$jacobian, tolerance = 1e-14, label = fx$name)
    expect_equal(as.numeric(r0$counts), ref, tolerance = 1e-13, label = fx$name)
    expect_identical(r0$objective, r$objective)
  }
})

test_that("the share Jacobian on sparse designs", {
  # helper-mxl-compact.R's designs: 30 alternatives, situations listing an
  # alternative twice or only the reference, the reference never offered,
  # trailing constants that no situation offers, one decision maker
  skip_if_not_installed("numDeriv")
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in Filter(function(fx) fx$use_asc, cdx_cells())) {
    fr <- blpr_free(fx)
    K_c <- length(fx$theta) - length(fr)
    r <- blpr_call(fx, variance = TRUE)
    r0 <- blpr_call(fx)
    g <- choicer:::mxl_loglik_gradient_parallel(
      fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
      fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
      Ti = fx[["Ti"]])
    expect_identical(r0$objective, g$objective, label = fx$name)
    expect_identical(r0$gradient, g$gradient, label = fx$name)
    ref <- blpr_counts(fx$theta, fx)
    expect_length(r0$counts, fx$J)
    expect_equal(as.numeric(r0$counts), ref, tolerance = 1e-13, label = fx$name)
    expect_equal(as.numeric(r$counts), ref, tolerance = 1e-13, label = fx$name)
    expect_equal(r$K, blpr_K(fx$theta, fx), tolerance = 1e-12, label = fx$name)
    # the Jacobian against numDeriv of the kernel's own counts (checked
    # against the reference above)
    jn <- numDeriv::jacobian(function(th) as.numeric(blpr_call(fx, theta = th)$counts)[fr],
                             fx$theta)
    d_idx <- K_c + seq_along(fr)
    expect_equal(r0$jacobian[d_idx, ], jn, tolerance = 1e-7, label = fx$name)
    expect_equal(r$jacobian, r0$jacobian, tolerance = 1e-14, label = fx$name)
  }
})

test_that("split draw batches, the shared accumulator and draws formed on the fly", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- blpj_cells()[[2L]]   # correlated, outside option, person weights
  base <- blpr_call(fx, variance = TRUE)
  # a unit's draws in batches of 5 forms its probabilities again: the same
  # matrices to rounding, and the gradient kernel's objective and gradient
  # at the same batches, bit for bit
  b <- blpr_call(fx, draw_batch = 5L, variance = TRUE)
  g <- choicer:::mxl_loglik_gradient_parallel(
    fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
    fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, TRUE, Ti = fx$Ti,
    draw_batch = 5L)
  expect_identical(b$objective, g$objective)
  expect_identical(b$gradient, g$gradient)
  expect_equal(b$jacobian, base$jacobian, tolerance = 1e-14)
  expect_equal(b$K, base$K, tolerance = 1e-14)
  expect_equal(blpr_call(fx, draw_batch = 5L)$jacobian, base$jacobian,
               tolerance = 1e-14)
  # two threads, each with its triangles (the default) or adding into the
  # shared result: the same results to rounding, the counts and the
  # objective's sum included
  set_omp_threads(2L)
  for (ab in list(NULL, 0)) {
    two <- if (is.null(ab)) {
      blpr_call(fx, variance = TRUE)
    } else {
      blpr_call(fx, acc_bytes = ab, variance = TRUE)
    }
    lab <- if (is.null(ab)) "triangles" else "shared"
    expect_equal(two$jacobian, base$jacobian, tolerance = 1e-13, label = lab)
    expect_equal(two$K, base$K, tolerance = 1e-13, label = lab)
    expect_equal(two$counts, base$counts, tolerance = 1e-14, label = lab)
    expect_equal(two$objective, base$objective, tolerance = 1e-14, label = lab)
    expect_equal(two$gradient, base$gradient, tolerance = 1e-13, label = lab)
  }
  set_omp_threads(1L)
  # generate mode with identity scrambling forms get_halton_normals()'
  # draws, the store-mode cube's
  gen <- blpr_call(fx, eta = array(0, c(fx$K_w, 0L, 0L)), gen_seed = 0L,
                   gen_scramble = 0L, gen_S = fx$S, variance = TRUE)
  expect_equal(gen$jacobian, base$jacobian, tolerance = 1e-13)
  expect_equal(gen$K, base$K, tolerance = 1e-13)
})

test_that("a non-finite weight makes the matrices and the counts NaN", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- blpj_cells()[[1L]]
  w <- fx$weights
  w[seq_len(fx$Ti[1L])] <- NaN   # the first decision maker's weight
  r <- blpr_call(fx, weights = w, variance = TRUE)
  expect_true(all(is.nan(r$jacobian)))
  expect_true(all(is.nan(r$K)))
  expect_true(all(is.nan(r$counts)))
  expect_error(choicer:::mxl_blp_jacobian_parallel(
    fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
    fx$rc_dist, fx$rc_correlation, fx$rc_mean, FALSE, FALSE, Ti = fx$Ti),
    "needs the constants")
})

test_that("the share Jacobian stops before allocating more than choicer.max_memory", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  # 400 alternatives and 200 draws: the variance mode's K accumulator and
  # each unit's draws of its counts are a large part of the estimate
  fx <- mxlp_fixture("memory", 3104, J = 400L, S = 200L, use_asc = TRUE)
  one_thread <- function(variance) {
    msg <- tryCatch(blpr_call(fx, variance = variance), error = conditionMessage)
    m <- regmatches(msg, regexec(
      "share Jacobian needs about ([0-9.]+) (kB|MB|GB) even with one thread", msg))[[1L]]
    as.numeric(m[2L]) * c(kB = 1e3, MB = 1e6, GB = 1e9)[[m[3L]]]
  }
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  expect_error(blpr_call(fx),
               "The share Jacobian needs about [0-9.]+ [kMG]B even with one thread .* more than the 0 kB of options\\(choicer.max_memory\\)")
  expect_gt(one_thread(TRUE), one_thread(FALSE))
  # Inf lifts the check, and a large limit does not bind
  options(choicer.max_memory = NULL)
  ref <- blpr_call(fx, variance = TRUE)
  for (lim in list(Inf, 1e15)) {
    options(choicer.max_memory = lim)
    expect_identical(blpr_call(fx, variance = TRUE), ref, label = format(lim))
  }
})

test_that("spd_solve() solves a positive definite system and reports why not", {
  set.seed(5)
  A <- crossprod(matrix(rnorm(400), 20)) + diag(20)
  A[1, ] <- A[1, ] * 1e6
  A[, 1] <- A[, 1] * 1e6
  B <- matrix(rnorm(60), 20)
  D <- choicer:::.jacobi_scale(A)
  s <- choicer:::spd_solve(A, B, D)
  expect_identical(s$info, 0L)
  expect_equal(s$x, solve(A, B), tolerance = 1e-12)
  # only the upper triangle is read
  A2 <- A
  A2[lower.tri(A2)] <- NaN
  expect_identical(choicer:::spd_solve(A2, B, D)$x, s$x)
  # not positive definite: the leading minor's order, zeros
  N <- A
  N[3, 3] <- -1
  n <- choicer:::spd_solve(N, B, rep(1, 20))
  expect_identical(n$info, 3L)
  expect_true(all(n$x == 0))
  # positive definite but singular to working precision: n + 1
  Sg <- diag(c(1, 1e-17, rep(1, 18)))
  expect_identical(choicer:::spd_solve(Sg, B, rep(1, 20))$info, 21L)
  # a value that is not finite: -1
  expect_identical(choicer:::spd_solve(A, replace(B, 2L, Inf), D)$info, -1L)
  # no right-hand side, no system
  expect_identical(dim(choicer:::spd_solve(A, B[, 0, drop = FALSE], D)$x), c(20L, 0L))
  expect_error(choicer:::spd_solve(A, B[-1L, ], D), "must be a double matrix with 20 rows")
  # the trailing block after `offset` rows and columns, read in place
  big <- matrix(NaN, 23, 23)
  big[4:23, 4:23] <- A
  expect_identical(choicer:::spd_solve(big, B, D, 3L), s)
  expect_error(choicer:::spd_solve(big, B, D, 24L), "between 0 and 23")
  expect_error(choicer:::spd_solve(big, B, D, 2L), "with 21 rows")
})
