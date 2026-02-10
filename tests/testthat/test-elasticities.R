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
    is_random_coef = FALSE
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
    is_random_coef = TRUE
  )

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})


# =============================================================================
# MNL Diversion Ratio Tests
# =============================================================================

test_that("mnl_diversion_ratios_parallel returns correct dimensions", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_equal(dim(dr), c(J, J))
})

test_that("mnl_diversion_ratios_parallel produces finite non-negative values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_true(all(is.finite(dr)))
  expect_true(all(dr >= 0))
})

test_that("mnl_diversion_ratios_parallel has zero diagonal", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_equal(diag(dr), rep(0, J))
})

test_that("mnl_diversion_ratios column sums equal 1 (no outside option, homogeneous choice sets)", {
  # With homogeneous choice sets and no outside option,
  # column sums of DR matrix should be 1 (all diverted demand stays inside)
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  col_sums <- colSums(dr)
  expect_equal(col_sums, rep(1, J), tolerance = 1e-10)
})

test_that("mnl_diversion_ratios IIA property: DR proportional to market share", {
  # Under MNL with IIA and identical choice sets, DR(j->k) = s_k / (1 - s_j)
  # This holds exactly when all individuals face the same covariates,
  # so we use zero X (only ASCs drive shares).
  set.seed(77)
  N <- 50
  J <- 3
  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1 = 0  # constant covariate -> ASCs fully determine shares
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", "x1")

  K <- ncol(inputs$X)
  J_out <- nrow(inputs$alt_mapping)
  # theta: beta for x1, then J-1 ASCs with different values
  theta <- c(0.5, 0.3, -0.2)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  shares <- as.vector(mnl_predict_shares(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  ))

  # For each column j, DR(j->k) = s_k / (1 - s_j)
  for (j in seq_len(J_out)) {
    expected_dr <- shares / (1 - shares[j])
    expected_dr[j] <- 0
    expect_equal(dr[, j], expected_dr, tolerance = 1e-10)
  }
})

test_that("mnl_diversion_ratios with outside option has correct dimensions and column sums", {
  set.seed(42)
  N <- 30
  J <- 4

  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J-1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  inputs <- prepare_mnl_data(
    dt, "id", "j", "choice", c("x1", "x2"),
    outside_opt_label = 0L,
    include_outside_option = TRUE
  )

  K <- ncol(inputs$X)
  # alt_mapping includes the outside option row; ASCs are for inside alts only
  J_asc <- nrow(inputs$alt_mapping) - 1
  theta <- runif(K + J_asc, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = TRUE
  )

  # Should be (J_asc + 1) x (J_asc + 1) matrix (inside alts + outside option)
  J_total <- J_asc + 1
  expect_equal(dim(dr), c(J_total, J_total))

  # Column sums should be 1
  col_sums <- colSums(dr)
  expect_equal(col_sums, rep(1, J_total), tolerance = 1e-10)

  # Diagonal should be 0
  expect_equal(diag(dr), rep(0, J_total))
})
