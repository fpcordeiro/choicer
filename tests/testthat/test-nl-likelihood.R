# Tests for Nested Logit likelihood and gradient functions
# Note: create_nl_inputs is defined in setup.R

test_that("nl_loglik_gradient_parallel returns correct structure", {
  inputs <- create_nl_inputs()

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)  # 2 non-singleton nests

  # theta: beta (K_x) + lambda (K_l) + delta (J-1)
  theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))

  result <- nl_loglik_gradient_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    nest_idx = inputs$nest_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = inputs$include_outside_option
  )

  expect_type(result, "list")
  expect_named(result, c("objective", "gradient"))
  expect_length(result$gradient, length(theta))
  expect_true(is.finite(result$objective))
  expect_true(all(is.finite(result$gradient)))
})

test_that("nl_loglik_gradient has accurate gradient", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 456)

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.3, 0.3),       # beta
    runif(K_l, 0.3, 0.8),        # lambda (must be positive)
    runif(J - 1, -0.2, 0.2)      # delta
  )

  obj_fn <- function(th) {
    nl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )$objective
  }

  analytic_grad <- nl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("nl_loglik_gradient is accurate at zero beta", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 789)

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  # Zero betas, non-zero lambdas and deltas
  theta <- c(rep(0, K_x), 0.6, 0.4, runif(J - 1, -0.1, 0.1))

  obj_fn <- function(th) {
    nl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )$objective
  }

  analytic_grad <- nl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("nl_loglik_gradient works without ASCs", {
  inputs <- create_nl_inputs()

  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  # theta: only beta and lambda (no delta)
  theta <- c(rep(0, K_x), rep(0.5, K_l))

  result <- nl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = FALSE, include_outside_option = inputs$include_outside_option
  )

  expect_length(result$gradient, K_x + K_l)
  expect_true(is.finite(result$objective))
})

test_that("nl_loglik_gradient is deterministic", {
  inputs <- create_nl_inputs()

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  theta <- c(0.3, -0.2, 0.6, 0.4, rep(0.1, J - 1))

  results <- replicate(5, {
    nl_loglik_gradient_parallel(
      theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )
  }, simplify = FALSE)

  # Allow tiny FP differences with parallel reduction
  objectives <- sapply(results, `[[`, "objective")
  expect_equal(max(objectives) - min(objectives), 0, tolerance = 1e-12)
})

test_that("nl_loglik_gradient handles lambda near boundary", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 321)

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  # Lambda near zero (but positive) and near one
  theta <- c(rep(0.1, K_x), 0.05, 0.95, rep(0, J - 1))

  obj_fn <- function(th) {
    nl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )$objective
  }

  analytic_grad <- nl_loglik_gradient_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )$gradient

  numeric_grad <- numDeriv::grad(obj_fn, theta, method = "Richardson")

  expect_equal(analytic_grad, numeric_grad, tolerance = TOL_GRAD)
})

test_that("nl_loglik_gradient returns finite with various lambda values", {
  inputs <- create_nl_inputs()

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  # Test various lambda combinations
  lambda_vals <- list(
    c(0.1, 0.1),
    c(0.5, 0.5),
    c(0.9, 0.9),
    c(0.2, 0.8)
  )

  for (lv in lambda_vals) {
    theta <- c(rep(0, K_x), lv, rep(0, J - 1))

    result <- nl_loglik_gradient_parallel(
      theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )

    expect_true(is.finite(result$objective),
                label = paste("lambda =", paste(lv, collapse = ",")))
    expect_true(all(is.finite(result$gradient)),
                label = paste("lambda =", paste(lv, collapse = ",")))
  }
})
