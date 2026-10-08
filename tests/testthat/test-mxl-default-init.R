# Tests for the default theta_init in run_mxlogit().
#
# When theta_init = NULL the start is in the columns' units
# (.mxl_default_start()):
#   - zeros on the beta, mu and ASC blocks and the Cholesky off-diagonal;
#   - on the Cholesky diagonal, log L_pp = -log(h_p) for a normal random
#     coefficient, h_p the typical step of its column (design_column_step():
#     the distance between a two-valued column's values, the mean distance
#     from a value filling more than half of the column, else the sd), and
#     log(0.5) for a log-normal one (the sd of the logarithm of its random
#     part);
#   - moved into lower/upper where it falls outside them.
# A start, not an estimator: the first two tests guard that the default
# reaches the same maximum as an explicit zero start.

test_that("default theta_init reaches the same MLE as zeros (uncorrelated)", {
  skip_on_cran()
  skip_on_ci()
  # outside_option=FALSE so J=3 inside alts -> J-1 = 2 free ASCs.
  sim <- simulate_mxl_data(N = 400L, J = 3L, seed = 31L, outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)

  common <- list(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    S = 30L,
    rc_mean = FALSE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    scale_vars = "none",
    control = list(xtol_rel = 1e-12, maxeval = 5000L)
  )

  # n_params = K_x + 0 + K_w + (J - 1) = 2 + 2 + 2 = 6
  n_params <- 6L

  fit_default <- suppressMessages(do.call(run_mxlogit, common))
  fit_zeros   <- suppressMessages(do.call(
    run_mxlogit,
    c(common, list(theta_init = rep(0, n_params)))
  ))

  # Both should converge so the comparison is meaningful.
  expect_true(fit_default$convergence > 0)
  expect_true(fit_zeros$convergence > 0)

  expect_equal(coef(fit_default), coef(fit_zeros), tolerance = 1e-5)
  expect_equal(
    as.numeric(logLik(fit_default)),
    as.numeric(logLik(fit_zeros)),
    tolerance = 1e-8
  )
})

test_that("default theta_init reaches the same MLE as zeros (correlated, rc_mean)", {
  skip_on_cran()
  skip_on_ci()
  # Cover the rc_correlation=TRUE branch where the Cholesky diagonal positions
  # within param_map$sigma are cumsum(seq_len(K_w)) rather than every entry,
  # and rc_mean=TRUE adds a mu block (which the default leaves at zero).
  sim <- simulate_mxl_data(N = 400L, J = 3L, seed = 17L, outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)

  common <- list(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    S = 30L,
    rc_mean = TRUE,
    rc_correlation = TRUE,
    use_asc = TRUE,
    scale_vars = "none",
    control = list(xtol_rel = 1e-12, maxeval = 5000L)
  )

  # n_params = K_x + K_w + K_w(K_w+1)/2 + (J - 1) = 2 + 2 + 3 + 2 = 9
  n_params <- 9L

  fit_default <- suppressMessages(do.call(run_mxlogit, common))
  fit_zeros   <- suppressMessages(do.call(
    run_mxlogit,
    c(common, list(theta_init = rep(0, n_params)))
  ))

  expect_true(fit_default$convergence > 0)
  expect_true(fit_zeros$convergence > 0)

  # MLE invariance: log-likelihood at the two optima must agree to the
  # convergence tolerance. This is the canonical regression guard for the
  # default-init change in the correlated+rc_mean case — coefficients along
  # flat directions of the Cholesky surface (off-diagonals and diagonals near
  # zero RC variance) can drift between starts while the likelihood is
  # identical, so logLik is the load-bearing invariant here.
  expect_equal(
    as.numeric(logLik(fit_default)),
    as.numeric(logLik(fit_zeros)),
    tolerance = 1e-8
  )
})

test_that("design_column_step() is a column's typical step", {
  set.seed(1)
  n <- 1000L
  dense <- cbind(rnorm(n), 1e3 * runif(n), 1e-6 * rexp(n), 1e8 + rnorm(n))
  expect_equal(choicer:::design_column_step(dense), unname(apply(dense, 2, stats::sd)),
               tolerance = 1e-13)
  # Two values: the distance between them, whatever their shares (a dummy's
  # step is 1, balanced or not). A value filling more than half of a column
  # with more values: the mean distance from it over the other rows.
  dummy <- as.numeric(seq_len(n) %% 40L == 0L)
  balanced <- rep(c(0, 1), length.out = n)
  sparse <- ifelse(seq_len(n) %% 5L == 0L, 3 + runif(n), 0)
  W <- cbind(dummy, balanced, 1 - dummy, 7 * dummy - 2, sparse)
  expect_equal(choicer:::design_column_step(W),
               c(1, 1, 1, 7, mean(sparse[sparse != 0])), tolerance = 1e-13)
  # Three balanced levels: the standard deviation
  three <- rep(c(0, 1, 2), length.out = n)
  expect_equal(choicer:::design_column_step(cbind(three)), stats::sd(three),
               tolerance = 1e-13)
  # No variation relative to the level, a step that is not finite, no rows: 0
  flat <- rep(3, n)
  nudged <- flat
  nudged[7] <- 3 * (1 + 4 * .Machine$double.eps)
  inf <- dense[, 1]
  inf[5] <- Inf
  na <- dense[, 1]
  na[9] <- NA
  inf_gap <- c(rep(0, 800), Inf, seq_len(n - 801))
  expect_identical(choicer:::design_column_step(cbind(flat, nudged, inf, na,
                                                      1e8 + 1e-6 * dense[, 1],
                                                      inf_gap)),
                   rep(0, 6))
  expect_identical(choicer:::design_column_step(dense[0, , drop = FALSE]), c(0, 0, 0, 0))
  expect_identical(choicer:::design_column_step(dense[1, , drop = FALSE]), c(0, 0, 0, 0))
  expect_identical(choicer:::design_column_step(dense[, 0, drop = FALSE]), numeric())
  expect_error(choicer:::design_column_step(matrix(1:4, 2)), "double matrix")
  expect_error(choicer:::design_column_step(c(1, 2)), "double matrix")
})

test_that("the default start is in the columns' units", {
  sim <- simulate_mxl_data(N = 60L, J = 3L, seed = 5L, outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  dt[, w1 := w1 * 1000]  # a column in other units
  dt[, w3 := as.numeric(x1 == max(x1)), by = id]  # a dummy, one per situation
  d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                        c("w1", "w2", "w3"))
  sdW <- unname(apply(d$W, 2, stats::sd))
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- theta_init
    stop("captured")
  })
  start <- function(rnd = c("w1", "w2"), ...) {
    expect_error(run_mxlogit(dt, "id", "alt", "choice", c("x1", "x2"), rnd,
                             S = 5L, scale_vars = "none", ...),
                 "captured")
    seen
  }
  # x1, x2, L_11, L_22, ASC_2, ASC_3
  expect_equal(start(), c(0, 0, -log(sdW[1:2]), 0, 0), tolerance = 1e-14)
  # A dummy's step is 1: L_33 = 1, where 1 / sd(w3) would be above 2
  expect_gt(1 / sdW[3], 2)
  expect_equal(start(c("w1", "w2", "w3")), c(0, 0, -log(sdW[1:2]), 0, 0, 0),
               tolerance = 1e-14)
  # x1, x2, Mu_w1, Mu_w2, L_11, L_21, L_22, ASC_2, ASC_3
  expect_equal(start(rc_correlation = TRUE, rc_mean = TRUE),
               c(0, 0, 0, 0, -log(sdW[1]), 0, -log(sdW[2]), 0, 0),
               tolerance = 1e-14)
  # A log-normal coefficient: L_pp = 0.5, the sd of the logarithm of its
  # random part (alone, and in a correlated factor)
  expect_equal(start(rc_dist = c(1L, 0L), rc_mean = TRUE),
               c(0, 0, 0, 0, log(0.5), -log(sdW[2]), 0, 0), tolerance = 1e-14)
  expect_equal(start(rc_dist = c(1L, 0L), rc_correlation = TRUE),
               c(0, 0, log(0.5), 0, -log(sdW[2]), 0, 0), tolerance = 1e-14)
  # Moved into the bounds, with a message; a supplied start is used as given
  expect_lt(-log(sdW[1]), -1)
  expect_gt(-log(sdW[2]), -1)
  expect_message(start(lower = c(L_11 = -1, x1 = -5), upper = c(L_22 = -1)),
                 "default start of L_11, L_22 lies outside")
  expect_identical(seen, c(0, 0, -1, -1, 0, 0))
  # A parameter pinned by lower == upper starts there without a word
  expect_message(start(lower = c(x1 = 0.3, L_11 = -1), upper = c(x1 = 0.3)),
                 "default start of L_11 lies outside")
  expect_identical(seen[1:3], c(0.3, 0, -1))
  expect_identical(start(theta_init = rep(0.1, 6), lower = c(L_11 = 1)),
                   rep(0.1, 6))
  # Bounds that hold no start
  for (b in list(list(lower = c(L_11 = 1), upper = c(L_11 = 0)),
                 list(lower = c(x1 = Inf)), list(upper = c(ASC_2 = -Inf)))) {
    expect_error(do.call(run_mxlogit, c(list(dt, "id", "alt", "choice",
                                             c("x1", "x2"), c("w1", "w2"),
                                             S = 5L), b)),
                 "Invalid bounds for (L_11|x1|ASC_2)")
  }
  # The start moves with the column's units
  pm <- list(beta = 1:2, sigma = 3:4, asc = 5:6)
  d2 <- d
  d2$W <- d$W[, 1:2]
  d2$W[, 1] <- d2$W[, 1] * 2^-7
  expect_equal(
    choicer:::.mxl_default_start(6L, pm, d2$W, c(0L, 0L), FALSE)[3:4],
    c(-log(sdW[1]) + 7 * log(2), -log(sdW[2])), tolerance = 1e-14)
  # An integer W (a hand-built input_data) is read as double
  Wi <- matrix(c(1L, 5L, 2L, 8L, 3L, 3L, 1L, 4L), 4)
  expect_identical(choicer:::.mxl_default_start(6L, pm, Wi, c(0L, 0L), FALSE),
                   choicer:::.mxl_default_start(6L, pm, Wi + 0, c(0L, 0L), FALSE))
  # No variation, or a value that is not finite: L_pp = 1; tiny units clamp
  Wz <- cbind(rep(4, 6), c(1, 2, Inf, 4, 5, 6))
  expect_identical(choicer:::.mxl_default_start(6L, pm, Wz, c(0L, 0L), FALSE),
                   rep(0, 6))
  Wt <- cbind(1e-152 * (1:6), 1:6)  # a unit below 2^-500, squares still normal
  expect_equal(choicer:::.mxl_default_start(6L, pm, Wt, c(0L, 0L), FALSE)[3],
               500 * log(2), tolerance = 1e-15)
})

test_that("from the default start a fit reaches the same maximum in any units", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_mxl_data(N = 1000L, J = 4L, beta = c(0.8, -0.5),
                           Sigma = diag(c(1, 0.6)), outside_option = FALSE,
                           seed = 1L)
  dt <- data.table::as.data.table(sim$data)
  fit <- function(d) suppressMessages(run_mxlogit(
    d, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"), S = 20L,
    draws = "generate", seed = 3L))
  a <- fit(dt)
  # The random coefficients' columns in other units: from L_pp = 0.5 in these
  # units (the 0.2.x start) the optimizer failed about 1,900
  # log-likelihood units short of the maximum.
  b <- fit(data.table::copy(dt)[, `:=`(w1 = w1 * 1e3, w2 = w2 * 1e-3)])
  expect_true(b$convergence %in% 1:4)
  expect_equal(b$loglik, a$loglik, tolerance = 1e-8)
  # The same model: row p of the Cholesky factor divided by the column's
  # factor (its log diagonal shifted by -log of it), the rest unchanged.
  shift <- c(L_11 = log(1e3), L_22 = log(1e-3))
  th_b <- coef(b)
  th_b[names(shift)] <- th_b[names(shift)] + shift
  expect_lt(max(abs(th_b - coef(a)) / a$se), 1e-3)
})

test_that("a fit whose default start lies outside the bounds starts at the bound", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_mxl_data(N = 400L, J = 3L, beta = c(0.8, -0.5),
                           Sigma = diag(c(1, 0.6)), outside_option = FALSE,
                           seed = 9L)
  dt <- data.table::as.data.table(sim$data)
  dt[, w1 := w1 * 1000]
  # The default start of L_11 is -log(sd(w1)), about -6.4, below the bound;
  # the maximum, near log(1 / 1000), lies below it too, so the bound binds.
  msgs <- capture_messages(fit <- run_mxlogit(
    dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"), S = 20L,
    draws = "generate", seed = 3L, lower = c(L_11 = -5)))
  expect_true(any(grepl("default start of L_11 lies outside", msgs)))
  expect_true(fit$convergence %in% 1:4)
  expect_equal(coef(fit)[["L_11"]], -5, tolerance = 1e-8)
})
