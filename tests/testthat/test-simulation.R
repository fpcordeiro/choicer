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
  # (sha256 of `data.table` with columns id, alt, x1, x2, choice). Any change
  # to the default DGP output will invalidate this snapshot.
  skip_if_not_installed("digest")
  sim <- simulate_mnl_data(seed = 123)
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "09e81a884419d936055720dda5133fe240f8b359b517de61492fade96a315183"
  )
})

test_that("simulate_mxl_data matches legacy helpers.R bit-for-bit at seed = 123", {
  skip_if_not_installed("digest")
  sim <- simulate_mxl_data(seed = 123)
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "dc793012b3a80ec7ee943d3ca6a260f5efe0745e9249699eeaba8103dfbed91d"
  )
})

test_that("simulate_nl_data matches legacy helpers.R bit-for-bit at seed = 123", {
  skip_if_not_installed("digest")
  sim <- simulate_nl_data(seed = 123)
  # NOTE: updated after L4 fix (nest column is now retained in the returned
  # data). If the DGP changes again, recompute with:
  #   digest::digest(simulate_nl_data(seed = 123)$data, algo = "sha256")
  expect_equal(
    digest::digest(sim$data, algo = "sha256"),
    "68771f00440cc776d32215892f657780729a4d82f92c76726d8eabe787bb4f33"
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
