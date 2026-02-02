# Tests for Mixed Logit likelihood and gradient functions

test_that("mxl_loglik_gradient_parallel returns correct structure", {
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

  # theta: beta (K_x) + L_params (K_w) + delta (J-1)
  theta <- rep(0, K_x + K_w + J - 1)

  result <- mxl_loglik_gradient_parallel(
    theta = theta,
    X = inputs$X,
    W = inputs$W,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("objective", "gradient"))
  expect_length(result$gradient, length(theta))
  expect_true(is.finite(result$objective))
  expect_true(all(is.finite(result$gradient)))
})

test_that("mxl_loglik_gradient has accurate gradient (uncorrelated)", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 123)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 30

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.3, 0.3),       # beta
    log(runif(K_w, 0.3, 0.8)),   # L_params (log scale)
    runif(J - 1, -0.2, 0.2)      # delta
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

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("mxl_loglik_gradient has accurate gradient (correlated)", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 456)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = TRUE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  L_size <- K_w * (K_w + 1) / 2
  S <- 30

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.3, 0.3),      # beta
    runif(L_size, -0.3, 0.3),   # L_params (includes off-diagonal)
    runif(J - 1, -0.2, 0.2)     # delta
  )

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rep(0L, K_w),
      rc_correlation = TRUE, rc_mean = FALSE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = TRUE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("mxl_loglik_gradient handles log-normal distribution", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 789)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 30

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.3, 0.3),      # beta
    log(runif(K_w, 0.3, 0.8)),  # L_params
    runif(J - 1, -0.2, 0.2)     # delta
  )

  # rc_dist = 1 for log-normal
  rc_dist <- rep(1L, K_w)

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rc_dist,
      rc_correlation = FALSE, rc_mean = FALSE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rc_dist,
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("mxl_loglik_gradient handles rc_mean = TRUE", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 321)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 30

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  # theta: beta (K_x) + mu (K_w) + L_params (K_w) + delta (J-1)
  theta <- c(
    runif(K_x, -0.3, 0.3),      # beta
    runif(K_w, -0.3, 0.3),      # mu (random coefficient means)
    log(runif(K_w, 0.3, 0.8)),  # L_params
    runif(J - 1, -0.2, 0.2)     # delta
  )

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rep(0L, K_w),
      rc_correlation = FALSE, rc_mean = TRUE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = TRUE,
    use_asc = TRUE, include_outside_option = FALSE
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("mxl_loglik_gradient is deterministic", {
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
  theta <- c(0.5, log(0.5), log(0.3), 0.1, -0.1)

  results <- replicate(5, {
    mxl_loglik_gradient_parallel(
      theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rep(0L, K_w),
      rc_correlation = FALSE, rc_mean = FALSE,
      use_asc = TRUE, include_outside_option = FALSE
    )
  }, simplify = FALSE)

  # Allow tiny FP differences with parallel reduction
  objectives <- sapply(results, `[[`, "objective")
  expect_equal(max(objectives) - min(objectives), 0, tolerance = 1e-12)
})

test_that("mxl_loglik_gradient returns finite values with large parameters", {
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

  # Large parameter values
  theta <- c(5, log(2), log(2), 3, -3)

  result <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  expect_true(is.finite(result$objective))
  expect_true(all(is.finite(result$gradient)))
})
