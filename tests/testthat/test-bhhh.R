# Tests for BHHH (OPG) information matrix for Mixed Logit

test_that("mxl_bhhh_parallel returns a symmetric matrix with correct dims", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)  # beta + L (diag) + ASCs (J-1)

  set.seed(123)
  theta <- runif(n_params, -0.3, 0.3)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("mxl_bhhh_parallel is positive semi-definite", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)

  set.seed(7)
  theta <- runif(n_params, -0.3, 0.3)
  eta <- get_halton_normals(S = 40, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )

  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-8))
})

test_that("run_mxlogit(se_method = 'bhhh') returns a valid choicer_mxl fit", {
  skip_on_cran()

  dt <- create_small_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 40L, se_method = "bhhh"
  )

  expect_s3_class(fit, "choicer_mxl")
  expect_identical(fit$se_method, "bhhh")
  expect_true(!is.null(fit$vcov))
  expect_true(!is.null(fit$se))
  expect_true(all(is.finite(fit$se)))
  expect_false(any(is.na(fit$se)))
})

test_that("BHHH and Hessian SEs have the same order of magnitude at convergence", {
  skip_on_cran()

  # Use a well-posed DGP with real signal so the MLE is identified.
  # The generic `create_small_mxl_data` helper picks choices uniformly at
  # random, which leaves the likelihood nearly flat and yields pathological
  # SEs for both estimators; that is not informative for comparing BHHH to
  # the analytical Hessian.
  set.seed(101)
  N <- 200; J <- 3
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), N)
  )
  dt[, x1 := rnorm(.N)]
  dt[, w1 := rnorm(.N)]
  U <- 1.0 * dt$x1 + 0.5 * dt$w1 + rnorm(N * J)
  dt[, U := U]
  dt[, choice := as.integer(U == max(U)), by = id]
  dt[, U := NULL]

  fit_hess <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    S = 50L
  )

  # Re-fit starting from converged point and use BHHH SEs.
  fit_bhhh <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    S = 50L, se_method = "bhhh",
    theta_init = unname(coef(fit_hess))
  )

  hess_se <- fit_hess$se
  bhhh_se <- fit_bhhh$se

  expect_true(all(is.finite(hess_se)))
  expect_true(all(is.finite(bhhh_se)))
  expect_true(all(hess_se > 0))
  expect_true(all(bhhh_se > 0))

  # Same order of magnitude. The information matrix equality holds only
  # asymptotically at the true MLE; with simulation noise the tolerance
  # must be generous.
  expect_true(all(abs(log(bhhh_se / hess_se)) < 1.0))
})

test_that("summary prints BHHH label without error", {
  skip_on_cran()

  dt <- create_small_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 40L, se_method = "bhhh"
  )

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("BHHH", out)))
})
