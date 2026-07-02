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

# --- HB DGPs (simulate_hmnl_data / simulate_hmnp_data) ---

test_that("simulate_hmnl_data is reproducible and carries the HB truth fields", {
  a <- simulate_hmnl_data(N = 30, T = 3, J = 4, seed = 123)
  b <- simulate_hmnl_data(N = 30, T = 3, J = 4, seed = 123)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "hmnl")

  tp <- a$true_params
  expect_named(tp, c("beta", "W", "theta", "sigma_d", "delta", "xi", "Z",
                     "rc_dist"))
  expect_length(tp$delta, 4L)
  expect_length(tp$xi, 4L)
  expect_equal(dim(tp$Z), c(4L, 2L))          # intercept + z1
  # the realized delta obeys the mean function delta = Z theta + xi
  expect_equal(tp$delta, drop(tp$Z %*% tp$theta) + tp$xi)
  expect_equal(tp$rc_dist, c(0L, 0L))

  # row counts: outside option adds exactly one alt = 0 row per task
  n_tasks <- 30L * 3L
  expect_equal(nrow(a$data), n_tasks * 5L)
  expect_equal(a$data[alt == 0L, .N], n_tasks)
  # z covariates constant within alternative, zero on the outside rows
  expect_equal(a$data[, uniqueN(z1), by = alt][["V1"]], rep(1L, 5L))
  expect_true(all(a$data[alt == 0L, z1 == 0]))
  # one (argmax) choice per task
  expect_equal(a$data[, sum(choice), by = task][["V1"]], rep(1L, n_tasks))
})

test_that("simulate_hmnl_data honors include_outside, vary_choice_set, rc_dist", {
  sim <- simulate_hmnl_data(N = 20, T = 2, J = 5, seed = 7,
                            rc_dist = c(1L, 0L),
                            include_outside = FALSE, vary_choice_set = TRUE)
  expect_true(!any(sim$data$alt == 0L))
  counts <- sim$data[, .N, by = task]
  expect_equal(nrow(counts), 40L)
  expect_true(all(counts$N >= 2L & counts$N <= 5L))
  expect_equal(sim$data[, sum(choice), by = task][["V1"]], rep(1L, 40L))
  # log-normal flag stored on the chain (log) scale convention
  expect_equal(sim$true_params$rc_dist, c(1L, 0L))
  expect_equal(sim$true_params$beta, c(0.8, -0.6))

  # input validation
  expect_error(simulate_hmnl_data(N = 5, T = 1, J = 3, W = diag(3)),
               "`W` must be a 2 x 2 matrix")
  expect_error(
    simulate_hmnl_data(N = 5, T = 1, J = 3, theta = c(0.5, -0.4),
                       Z = matrix(0, 3, 2)),
    "`Z` must be a 3 x 1 matrix"
  )
})

test_that("simulate_hmnp_data reports truth on the identified scale", {
  sigma <- 2
  a <- simulate_hmnp_data(N = 25, T = 2, J = 3, seed = 11, sigma = sigma)
  b <- simulate_hmnp_data(N = 25, T = 2, J = 3, seed = 11, sigma = sigma)
  expect_equal(a$data, b$data)
  expect_equal(a$true_params, b$true_params)
  expect_s3_class(a, "choicer_sim")
  expect_equal(a$model, "hmnp")

  tp <- a$true_params
  expect_named(tp, c("beta", "W", "theta", "sigma_d", "delta", "xi", "Z"))
  expect_equal(tp$beta, c(0.8, -0.6) / sigma)
  expect_equal(tp$W, diag(0.5, 2) / sigma^2)
  expect_equal(tp$theta, c(0.5, -0.4) / sigma)
  expect_equal(tp$sigma_d, 0.5 / sigma)
  # the identified-scale delta still obeys the mean function
  expect_equal(tp$delta, drop(tp$Z %*% tp$theta) + tp$xi)
  expect_equal(a$settings$sigma, sigma)

  # rows: 50 tasks x (3 inside + outside)
  expect_equal(nrow(a$data), 50L * 4L)
  expect_equal(a$data[alt == 0L, .N], 50L)
  expect_equal(a$data[, sum(choice), by = task][["V1"]], rep(1L, 50L))

  expect_error(simulate_hmnp_data(N = 5, T = 1, J = 3, sigma = 0),
               "`sigma` must be a positive number")
})
