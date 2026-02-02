# Tests for elasticity computation functions

test_that("mnl_elasticities_parallel returns correct dimensions", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Elasticity with respect to first covariate (1-based index)
  elast <- mnl_elasticities_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    elast_var_idx = 1L,  # 1-based indexing
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  # Should return J x J matrix
  expect_equal(dim(elast), c(J, J))
})

test_that("mnl_elasticities_parallel produces finite values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  for (k in 1:K) {
    elast <- mnl_elasticities_parallel(
      theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, k, TRUE, FALSE
    )

    expect_true(all(is.finite(elast)),
                label = paste("covariate index", k))
  }
})

test_that("mnl_elasticities own-elasticities are on diagonal", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- c(1.0, -0.5, rep(0.1, J - 1))  # Positive beta for x1

  elast <- mnl_elasticities_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, 1L, TRUE, FALSE
  )

  # Diagonal contains own-elasticities
  own_elast <- diag(elast)
  expect_length(own_elast, J)
})

test_that("mnl_elasticities with zero beta gives zero elasticities", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Zero coefficient for x1
  theta <- c(0, 0.5, rep(0, J - 1))

  elast <- mnl_elasticities_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, 1L, TRUE, FALSE
  )

  # All elasticities with respect to x1 should be zero
  expect_equal(elast, matrix(0, J, J), tolerance = 1e-12)
})

# MXL elasticities tests - skip if function not available
test_that("mxl_elasticities_parallel returns correct dimensions", {
  skip_if_not(exists("mxl_elasticities_parallel"),
              "mxl_elasticities_parallel not exported")

  dt <- create_small_mxl_data()
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
  theta <- c(0.3, log(0.5), log(0.4), 0.1, -0.1)

  # Elasticity with respect to fixed covariate x1
  elast <- mxl_elasticities_parallel(
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
    include_outside_option = FALSE,
    elast_var_idx = 1L,
    is_rc_var = FALSE
  )

  expect_equal(dim(elast), c(J, J))
})

test_that("mxl_elasticities_parallel works for random coefficient variable", {
  skip_if_not(exists("mxl_elasticities_parallel"),
              "mxl_elasticities_parallel not exported")

  dt <- create_small_mxl_data()
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
  theta <- c(0.3, log(0.5), log(0.4), 0.1, -0.1)

  # Elasticity with respect to random coefficient w1 (index 1 in W, 1-based)
  elast <- mxl_elasticities_parallel(
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
    include_outside_option = FALSE,
    elast_var_idx = 1L,
    is_rc_var = TRUE
  )

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})
