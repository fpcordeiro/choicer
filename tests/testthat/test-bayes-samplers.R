# Tests for the from-scratch RNG core (rng.h) and the MCMC sampling
# primitives (bayes_samplers.h), reached via the internal *_cpp wrappers.

test_that("rnorm_cpp draws standard normals", {
  x <- choicer:::rnorm_cpp(2e5, 7)
  expect_equal(mean(x), 0, tolerance = 0.02)
  expect_equal(sd(x), 1, tolerance = 0.02)
  expect_equal(mean(x^3), 0, tolerance = 0.05)   # skewness
  expect_equal(mean(x^4), 3, tolerance = 0.15)   # kurtosis
})

test_that("rnorm_cpp is reproducible by seed", {
  expect_identical(choicer:::rnorm_cpp(100, 42), choicer:::rnorm_cpp(100, 42))
  expect_false(identical(choicer:::rnorm_cpp(100, 42),
                         choicer:::rnorm_cpp(100, 43)))
})

test_that("rgamma_cpp matches gamma moments incl. non-integer shape", {
  for (a in c(0.5, 1.7, 2.5)) {
    g <- as.numeric(choicer:::rgamma_cpp(2e5, a, 11))
    expect_true(all(g > 0))
    expect_equal(mean(g), a, tolerance = 0.03)
    expect_equal(var(g), a, tolerance = 0.08)
  }
})

test_that("rtruncnorm_cpp: one-sided lower truncation matches theory", {
  x <- choicer:::rtruncnorm_cpp(2e5, 0, 1, 1, Inf, 3)
  expect_true(all(x > 1))
  expect_equal(mean(x), dnorm(1) / pnorm(1, lower.tail = FALSE),
               tolerance = 0.01)
})

test_that("rtruncnorm_cpp: far tail (a = 8) is finite and exact", {
  # The naive qnorm(runif(pnorm(a), 1)) approach returns Inf here; this
  # exercises the Robert (1995) exponential-rejection branch.
  x <- choicer:::rtruncnorm_cpp(2e5, 0, 1, 8, Inf, 3)
  expect_true(all(is.finite(x)))
  expect_true(all(x > 8))
  expect_equal(mean(x), dnorm(8) / pnorm(8, lower.tail = FALSE),
               tolerance = 0.001)
})

test_that("rtruncnorm_cpp: two-sided and mirrored intervals", {
  x <- choicer:::rtruncnorm_cpp(2e5, 0, 1, -1, 0.5, 3)
  expect_true(all(x > -1 & x < 0.5))
  m <- (dnorm(-1) - dnorm(0.5)) / (pnorm(0.5) - pnorm(-1))
  expect_equal(mean(x), m, tolerance = 0.01)

  # Interval entirely in the lower tail (mirror branch)
  y <- choicer:::rtruncnorm_cpp(2e5, 0, 1, -Inf, -8, 3)
  expect_true(all(is.finite(y)))
  expect_true(all(y < -8))
  expect_equal(mean(y), -dnorm(8) / pnorm(8, lower.tail = FALSE),
               tolerance = 0.001)
})

test_that("rtruncnorm_cpp respects location and scale", {
  x <- choicer:::rtruncnorm_cpp(2e5, 2, 3, 2, Inf, 5)
  expect_true(all(x > 2))
  expect_equal(mean(x), 2 + 3 * dnorm(0) / 0.5, tolerance = 0.02)
})

test_that("rmvnorm_cpp recovers mean and covariance", {
  Sg <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
  mu <- c(1, -1)
  x <- choicer:::rmvnorm_cpp(5e4, mu, Sg, 5)
  expect_equal(colMeans(x), mu, tolerance = 0.03)
  expect_equal(cov(x), Sg, tolerance = 0.05)
})

test_that("rwishart_cpp has mean df * S and draws symmetric PD matrices", {
  Sg <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
  df <- 7
  draws <- lapply(1:1000, function(i) choicer:::rwishart_cpp(df, Sg, i))
  avg <- Reduce(`+`, draws) / length(draws)
  expect_equal(avg, df * Sg, tolerance = 0.1)
  for (W in draws[1:20]) {
    expect_equal(W, t(W))
    expect_true(all(eigen(W, only.values = TRUE)$values > 0))
  }
})

test_that("riwishart_cpp has mean V / (df - p - 1) and draws PD matrices", {
  Sg <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
  df <- 8   # p = 2 -> prior mean Sg / 5
  draws <- lapply(1:2000, function(i) choicer:::riwishart_cpp(df, Sg, i))
  avg <- Reduce(`+`, draws) / length(draws)
  expect_equal(avg, Sg / (df - 2 - 1), tolerance = 0.1)
  expect_true(all(eigen(draws[[1]], only.values = TRUE)$values > 0))
})

test_that("rwishart_cpp rejects df < p", {
  expect_error(choicer:::rwishart_cpp(1, diag(2), 1), "degrees of freedom")
})
