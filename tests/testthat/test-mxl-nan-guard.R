# Tests for the NaN-safe likelihood guard in mxl_loglik_gradient_parallel().
#
# When the C++ likelihood evaluation produces a non-finite objective, the
# function returns objective = 1e10 and gradient = zeros (a sentinel that lets
# nloptr/optim continue searching without poisoning the line search). When the
# objective is finite but some gradient entries are non-finite, those entries
# are zeroed.
#
# We exercise the guard by evaluating at a deliberately pathological theta
# (rep(1e6, n_params)) where the random-coefficient utilities explode and the
# log-sum-exp returns Inf / NaN in the un-guarded path.

test_that("mxl_loglik_gradient_parallel returns sentinel at pathological theta", {
  sim <- simulate_mxl_data(N = 200L, J = 4L, seed = 11L)
  dt <- data.table::as.data.table(sim$data)

  inputs <- prepare_mxl_data(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    rc_correlation = FALSE
  )

  K_w <- ncol(inputs$W)
  K_x <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  rc_mean <- FALSE
  rc_correlation <- inputs$rc_correlation
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- J - 1
  n_params <- K_x + mu_size + L_size + n_asc

  eta_draws <- get_halton_normals(S = 30L, N = inputs$N, K_w = K_w)

  bad_theta <- rep(1e6, n_params)
  result <- mxl_loglik_gradient_parallel(
    theta = bad_theta,
    X = inputs$X,
    W = inputs$W,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    use_asc = TRUE,
    include_outside_option = inputs$include_outside_option
  )

  # Sentinel objective: must be finite (NOT NaN/Inf) and exactly 1e10.
  expect_true(is.finite(result$objective))
  expect_equal(result$objective, 1e10)
  expect_true(result$overflow)

  # Gradient must be entirely finite.
  expect_true(all(is.finite(result$gradient)))
  expect_length(result$gradient, n_params)
})

test_that("the overflow sentinel stays above every objective seen", {
  # Finite objectives are unbounded (probabilities are handled in log space),
  # so after a start above 1e10 the kernel's fixed sentinel must be lifted, or
  # the line search would take an overflowing trial point for an improvement.
  vals <- c(2e10, 1e10, 5, 1e10, 3e10, 1e10)
  i <- 0L
  f <- choicer:::.lift_sentinel(function(theta) {
    i <<- i + 1L
    list(objective = vals[i], gradient = 0, overflow = i %% 2L == 0L)
  })
  got <- vapply(seq_along(vals), function(k) f(0)$objective, numeric(1))
  expect_equal(got, c(2e10, 2e11, 5, 2e11, 3e10, 3e11))
  # From a start below 1e10 the sentinel is unchanged.
  i <- 0L
  vals <- c(100, 1e10)
  g <- choicer:::.lift_sentinel(function(theta) {
    i <<- i + 1L
    list(objective = vals[i], gradient = 0, overflow = i == 2L)
  })
  expect_equal(c(g(0)$objective, g(0)$objective), c(100, 1e10))
})

test_that("a finite objective equal to the sentinel is preserved", {
  kernel <- function(beta) {
    mxl_loglik_gradient_parallel(
      theta = c(beta, 0), X = matrix(c(-1, 0), 2L), W = matrix(0, 2L),
      alt_idx = 1:2, choice_idx = 1L, M = 2L, weights = 1,
      eta_draws = array(0, c(1L, 1L, 1L)), rc_dist = 0L,
      rc_correlation = FALSE, rc_mean = FALSE, use_asc = FALSE,
      include_outside_option = FALSE)
  }
  f <- choicer:::.lift_sentinel(kernel)
  expect_equal(f(2e10)$objective, 2e10)
  raw <- kernel(1e10)
  expect_identical(raw$objective, 1e10)
  expect_false(raw$overflow)
  res <- f(1e10)
  expect_identical(res$objective, raw$objective)
  expect_identical(res$gradient, raw$gradient)
  expect_equal(drop(res$gradient), c(1, 0))
  expect_named(res, c("objective", "gradient"))
  # A real overflow is still penalized above the largest valid objective.
  expect_true(kernel(Inf)$overflow)
  bad <- f(Inf)
  expect_equal(bad$objective, 2e11)
  expect_true(all(bad$gradient == 0))
})
