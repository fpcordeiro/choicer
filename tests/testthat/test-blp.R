# Tests for BLP contraction mapping functions

test_that("blp_contraction returns correct length", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- rep(1/J, J)  # Equal shares
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  expect_length(result, J)
})

test_that("blp_contraction converges to target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  K <- ncol(inputs$X)

  # Use observed market shares as target
  target_shares <- inputs$alt_mapping$MKT_SHARE

  beta <- c(0.3, -0.2)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  # Verify shares match target
  # Construct theta: beta + delta (excluding first delta for identification)
  theta <- c(beta, delta_result[2:J])

  pred_shares <- mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  )
  # Convert from matrix to vector if needed
  pred_shares <- as.vector(pred_shares)

  expect_equal(pred_shares, target_shares, tolerance = 1e-6)
})

test_that("blp_contraction handles equal target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- rep(1/J, J)
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  # Should converge without error
  expect_true(all(is.finite(delta_result)))
})

test_that("blp_contraction returns finite values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- c(0.5, 0.3, 0.2)
  beta <- c(1.0, -0.5)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  expect_true(all(is.finite(delta_result)))
})

# MXL BLP contraction tests - skip if function not available
test_that("mxl_blp_contraction returns correct length", {
  skip_if_not(exists("mxl_blp_contraction"),
              "mxl_blp_contraction not exported")

  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  target_shares <- rep(1/J, J)
  beta <- 0.3
  L_params <- c(log(0.5), log(0.4))
  delta_init <- rep(0, J)

  result <- mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    W = inputs$W,
    beta = beta,
    mu = rep(0, K_w),
    L_params = L_params,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_length(result, J)
})

test_that("mxl_blp_contraction returns finite values", {
  skip_if_not(exists("mxl_blp_contraction"),
              "mxl_blp_contraction not exported")

  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  J <- nrow(inputs$alt_mapping)
  K_w <- ncol(inputs$W)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  target_shares <- c(0.4, 0.35, 0.25)
  beta <- 0.5
  L_params <- c(log(0.5), log(0.4))
  delta_init <- rep(0, J)

  result <- mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    W = inputs$W,
    beta = beta,
    mu = rep(0, K_w),
    L_params = L_params,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_true(all(is.finite(result)))
})

test_that("blp_contraction respects convergence tolerance", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- inputs$alt_mapping$MKT_SHARE
  beta <- c(0.3, -0.2)
  delta_init <- rep(0, J)

  # Tighter tolerance
  delta_tight <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-12
  )

  # Looser tolerance
  delta_loose <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-6
  )

  # Both should be close but tight should be more accurate
  theta_tight <- c(beta, delta_tight[2:J])
  theta_loose <- c(beta, delta_loose[2:J])

  shares_tight <- as.vector(mnl_predict_shares(
    theta_tight, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))
  shares_loose <- as.vector(mnl_predict_shares(
    theta_loose, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))

  error_tight <- max(abs(shares_tight - target_shares))
  error_loose <- max(abs(shares_loose - target_shares))

  # Tight tolerance should give smaller error
  expect_true(error_tight <= error_loose)
})

test_that("blp_contraction handles skewed target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)

  # Very skewed shares
  target_shares <- c(0.9, 0.08, 0.02)
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_true(all(is.finite(delta_result)))

  # Verify shares are recovered
  theta <- c(beta, delta_result[2:J])
  pred_shares <- as.vector(mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))

  expect_equal(pred_shares, target_shares, tolerance = 1e-5)
})
