# Tests for MNL likelihood and gradient functions

test_that("mnl_loglik_gradient_parallel returns correct structure", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  result <- mnl_loglik_gradient_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("objective", "gradient"))
  expect_length(result$gradient, length(theta))
  expect_true(is.finite(result$objective))
  expect_true(all(is.finite(result$gradient)))
})

test_that("mnl_loglik_gradient returns known value for equal probabilities", {
  # With zero parameters and zero covariates, all alternatives have equal probability
  dt <- create_equal_prob_data(N = 100, J = 3)
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", "x1")

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)  # All zeros

  result <- mnl_loglik_gradient_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  # With equal probabilities P = 1/J, negative log-likelihood = -N * log(1/J) = N * log(J)
  expected_nll <- 100 * log(3)
  expect_equal(result$objective, expected_nll, tolerance = TOL_LOGLIK)
})

test_that("mnl_loglik_gradient has accurate gradient (numDeriv comparison)", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Test at random parameter values
  set.seed(42)
  theta <- runif(K + J - 1, -0.5, 0.5)

  # Objective function wrapper
  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  # Analytical gradient
  analytic_grad <- mnl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )$gradient

  # Numerical gradient
  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(drop(analytic_grad), numeric_grad, tolerance = TOL_GRAD)
})

test_that("mnl_loglik_gradient is accurate at zero parameters", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  analytic_grad <- mnl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(drop(analytic_grad), numeric_grad, tolerance = TOL_GRAD)
})

test_that("mnl_loglik_gradient works without ASCs", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  theta <- rep(0, K)  # Only beta, no ASCs

  result <- mnl_loglik_gradient_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = FALSE,
    include_outside_option = FALSE
  )

  expect_length(result$gradient, K)
  expect_true(is.finite(result$objective))
})

test_that("mnl_loglik_gradient handles observation weights", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Equal weights
  result_equal <- mnl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  # Double the weights
  double_weights <- inputs$weights * 2
  result_double <- mnl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, double_weights, TRUE, FALSE
  )

  # Objective should double, gradient should double

  expect_equal(result_double$objective, result_equal$objective * 2, tolerance = TOL_LOGLIK)
  expect_equal(result_double$gradient, result_equal$gradient * 2, tolerance = TOL_LOGLIK)
})

test_that("mnl_loglik_gradient is deterministic (OpenMP stability)", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  # Run multiple times
  results <- replicate(5, {
    mnl_loglik_gradient_parallel(
      theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )
  }, simplify = FALSE)

  # All objectives should be very close (tiny FP differences possible with parallel reduction)
  objectives <- sapply(results, `[[`, "objective")
  expect_equal(max(objectives) - min(objectives), 0, tolerance = 1e-12)

  # All gradients should be very close
  for (i in 2:length(results)) {
    expect_equal(results[[1]]$gradient, results[[i]]$gradient, tolerance = 1e-12)
  }
})

test_that("mnl_loglik_gradient handles large parameter values gracefully", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Large positive parameters
  theta_large <- rep(10, K + J - 1)

  result <- mnl_loglik_gradient_parallel(
    theta_large, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  # Should not produce NaN or Inf
  expect_true(is.finite(result$objective))
  expect_true(all(is.finite(result$gradient)))
})

test_that("mnl_loglik_gradient handles single individual", {
  dt <- data.table(
    id = rep(1L, 3),
    alt = 1:3,
    choice = c(1L, 0L, 0L),
    x1 = c(0.5, -0.3, 0.1),
    x2 = c(0.2, 0.8, -0.5)
  )
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  result <- mnl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  expect_true(is.finite(result$objective))
  expect_length(result$gradient, length(theta))
})
