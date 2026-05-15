# Tests for the default theta_init in run_mxlogit().
#
# When theta_init = NULL the cold start is:
#   - zeros on the beta, mu, and ASC blocks
#   - log(0.5) on the Cholesky diagonal positions in param_map$sigma
# Starting from log(0.5) on the diagonal (RC variance 0.25) is a better-
# conditioned cold start than log(1) = 0 (unit RC variance), but it MUST
# reach the same MLE as the explicit zero start. This test is the regression
# guard that the new default is just a starting-point change, not a different
# optimum.

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
