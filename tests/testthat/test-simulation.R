# Tests for R/simulation.R DGP functions.

test_that("simulate_mnl_data is reproducible across calls at fixed seed", {
  a <- simulate_mnl_data(seed = 123)
  b <- simulate_mnl_data(seed = 123)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "mnl")
})

test_that("simulate_mxl_data is reproducible across calls at fixed seed", {
  a <- simulate_mxl_data(seed = 123)
  b <- simulate_mxl_data(seed = 123)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "mxl")
})

test_that("simulate_nl_data is reproducible across calls at fixed seed", {
  a <- simulate_nl_data(seed = 123)
  b <- simulate_nl_data(seed = 123)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "nl")
})

test_that("simulate_mnp_data is reproducible across calls at fixed seed", {
  a <- simulate_mnp_data(N = 200, seed = 123)
  b <- simulate_mnp_data(N = 200, seed = 123)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "mnp")
})

test_that("simulate_mnp_data produces balanced sets and identified truth", {
  Sigma <- matrix(c(4, 1, 1, 6), 2, 2)   # sigma_11 = 4 -> non-unit DGP scale
  beta <- c(0.8, -0.6)
  delta <- c(0.5, -0.5)
  sim <- simulate_mnp_data(N = 100, J = 3, beta = beta, delta = delta,
                           Sigma = Sigma, seed = 7L)

  counts <- sim$data[, .N, by = id]
  expect_true(all(counts$N == 3L))
  expect_equal(sim$data[, sum(choice), by = id]$V1, rep(1L, 100))

  # true_params are reported on the identified (sigma_11 = 1) scale
  expect_equal(sim$true_params$Sigma[1, 1], 1)
  expect_equal(sim$true_params$beta, beta / 2)
  expect_equal(sim$true_params$delta, delta / 2)
  expect_equal(sim$true_params$Sigma, Sigma / 4)

  expect_error(simulate_mnp_data(N = 10, J = 4, Sigma = diag(2)), "Sigma")
  expect_error(simulate_mnp_data(N = 10, J = 3, delta = 1:3), "delta")
})

test_that("simulate_mxl_data supports price_cols and rc_dist extensions", {
  sim <- simulate_mxl_data(
    N = 200, J = 4, seed = 1L,
    beta = c(0.5, -0.4),
    mu = c(0.1, 0.2),
    Sigma = matrix(c(0.8, 0.1, 0.1, 0.6), nrow = 2),
    rc_dist = c(1L, 0L),
    price_cols = "w1"
  )
  # Every non-zero w1 value should be negative (price column).
  w1_nonzero <- sim$data[sim$data$w1 != 0, ]$w1
  expect_true(all(w1_nonzero < 0))
  expect_equal(sim$true_params$rc_dist, c(1L, 0L))
  expect_equal(sim$true_params$mu, c(0.1, 0.2))
  expect_true(sim$true_params$rc_correlation)
})

test_that("simulate_mnl_data honors outside_option = FALSE and vary_choice_set = FALSE", {
  sim <- simulate_mnl_data(
    N = 50, J = 4, seed = 7L,
    outside_option = FALSE,
    vary_choice_set = FALSE
  )
  counts <- sim$data[, .N, by = id]
  expect_true(all(counts$N == 4L))
  expect_true(!any(sim$data$alt == 0L))
})
