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

test_that("simulate_mnl_data matches legacy helpers.R bit-for-bit at seed = 123", {
  # Digest snapshot captured from the pre-refactor `inst/simulations/helpers.R`
  # (sha256 of `data.table` with columns id, alt, x1, x2, choice). Hash is
  # platform-sensitive (data.table attribute order, integer representation),
  # so this regression check runs locally only.
  skip_on_cran()
  skip_if_not_installed("digest")
  sim <- simulate_mnl_data(seed = 123)
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "e15122bc6da53d1c53dcce5c4a93f77aab8ab3cace14414f60598e3eca893024"
  )
})

test_that("simulate_mxl_data matches legacy helpers.R bit-for-bit at seed = 123", {
  skip_on_cran()
  skip_if_not_installed("digest")
  sim <- simulate_mxl_data(seed = 123)
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "b797e867fa78eac3acad445f963f91e47d05349e32523c88523a3d2c087f1a52"
  )
})

test_that("simulate_nl_data matches legacy helpers.R bit-for-bit at seed = 123", {
  skip_on_cran()
  skip_if_not_installed("digest")
  sim <- simulate_nl_data(seed = 123)
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "3a64b533d3c514590e26fe5d9d7e6a3e17ca7a5a503ad361348146954a4a7e6a"
  )
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
