# Tests for prediction functions

test_that("mnl_predict returns correct structure", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
  )

  expect_type(preds, "list")
  expect_true("choice_prob" %in% names(preds))
  expect_true("utility" %in% names(preds))
  expect_length(preds$choice_prob, nrow(inputs$X))
  expect_length(preds$utility, nrow(inputs$X))
})

test_that("mnl_predict probabilities sum to 1 within individuals", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
  )

  # Add probabilities back to data
  dt$prob <- preds$choice_prob

  # Sum probabilities within each individual
  prob_sums <- dt[, sum(prob), by = id]$V1

  expect_equal(prob_sums, rep(1, length(prob_sums)), tolerance = TOL_PROB)
})

test_that("mnl_predict probabilities are between 0 and 1", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Test with various parameter values
  theta_sets <- list(
    rep(0, K + J - 1),
    runif(K + J - 1, -1, 1),
    c(rep(5, K), rep(-5, J - 1))
  )

  for (theta in theta_sets) {
    preds <- mnl_predict(
      theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
    )

    expect_true(all(preds$choice_prob >= 0))
    expect_true(all(preds$choice_prob <= 1))
  }
})

test_that("mnl_predict with zero parameters gives equal probabilities", {
  dt <- create_equal_prob_data(N = 50, J = 4)
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", "x1")

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
  )

  # All probabilities should be 1/J = 0.25
  unique_probs <- unique(as.vector(round(preds$choice_prob, 8)))
  expect_equal(unique_probs, 1/J, tolerance = TOL_PROB)
})

test_that("mnl_predict_shares sum to 1", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  shares <- as.vector(mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M, inputs$weights, TRUE, FALSE
  ))

  expect_length(shares, J)
  expect_equal(sum(shares), 1, tolerance = TOL_PROB)
})

test_that("mnl_predict_shares are non-negative", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -2, 2)

  shares <- as.vector(mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M, inputs$weights, TRUE, FALSE
  ))

  expect_true(all(shares >= 0))
})

test_that("mnl_predict_shares with equal weights match average probabilities", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Get shares
  shares <- as.vector(mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M, inputs$weights, TRUE, FALSE
  ))

  # Calculate manually from individual probabilities
  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
  )

  dt$prob <- preds$choice_prob
  manual_shares <- dt[, mean(prob), by = alt][order(alt)]$V1

  expect_equal(shares, manual_shares, tolerance = TOL_PROB)
})

test_that("mnl_predict handles large parameter values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Large positive parameters
  theta <- rep(10, K + J - 1)

  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, TRUE, FALSE
  )

  # Should still produce valid probabilities
  expect_true(all(is.finite(preds$choice_prob)))
  expect_true(all(preds$choice_prob >= 0))
  expect_true(all(preds$choice_prob <= 1))

  dt$prob <- preds$choice_prob
  prob_sums <- dt[, sum(prob), by = id]$V1
  expect_equal(prob_sums, rep(1, length(prob_sums)), tolerance = TOL_PROB)
})

test_that("mnl_predict without ASCs has correct parameter count", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  theta <- runif(K, -0.3, 0.3)  # Only beta, no ASCs

  preds <- mnl_predict(
    theta, inputs$X, inputs$alt_idx, inputs$M, FALSE, FALSE
  )

  expect_length(preds$choice_prob, nrow(inputs$X))
  expect_true(all(is.finite(preds$choice_prob)))
})
