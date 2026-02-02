# Tests for Hessian computations

# --- MNL Hessian tests ---

test_that("mnl_loglik_hessian_parallel matches numerical Hessian", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  set.seed(42)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Analytical Hessian
  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  # Numerical Hessian
  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

test_that("mnl_loglik_hessian_parallel is symmetric", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  H <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  expect_equal(H, t(H), tolerance = 1e-12)
})

test_that("mnl_loglik_numeric_hessian matches analytical", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.2, 0.2)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  H_num <- mnl_loglik_numeric_hessian(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

test_that("mnl_loglik_hessian at zero parameters", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

# --- MXL Hessian tests ---

test_that("mxl_hessian_parallel matches numerical Hessian", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 111)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.2, 0.2),
    log(runif(K_w, 0.4, 0.7)),
    runif(J - 1, -0.1, 0.1)
  )

  H_anal <- mxl_hessian_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rep(0L, K_w),
      rc_correlation = FALSE, rc_mean = FALSE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

test_that("mxl_hessian_parallel is symmetric", {
  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 20

  eta_draws <- get_halton_normals(S, N, K_w)
  theta <- c(0.1, log(0.5), log(0.4), 0.1, -0.1)

  H <- mxl_hessian_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("mxl_loglik_numeric_hessian returns correct dimensions", {
  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 20

  eta_draws <- get_halton_normals(S, N, K_w)
  n_params <- K_x + K_w + J - 1
  theta <- rep(0, n_params)

  H <- mxl_loglik_numeric_hessian(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_true(all(is.finite(H)))
})

# --- NL Hessian tests ---

test_that("nl_loglik_numeric_hessian returns correct dimensions", {
  inputs <- create_nl_inputs()

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  n_params <- K_x + K_l + J - 1
  theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))

  H <- nl_loglik_numeric_hessian(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_true(all(is.finite(H)))
})

test_that("nl_loglik_numeric_hessian matches numDeriv", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 222)

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  set.seed(42)
  theta <- c(runif(K_x, -0.2, 0.2), runif(K_l, 0.3, 0.7), runif(J - 1, -0.1, 0.1))

  H_pkg <- nl_loglik_numeric_hessian(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  obj_fn <- function(th) {
    nl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_pkg, H_num, tolerance = TOL_HESS)
})

# Note: create_nl_inputs is defined in setup.R
