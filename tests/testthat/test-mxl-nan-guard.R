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

  # Gradient must be entirely finite.
  expect_true(all(is.finite(result$gradient)))
  expect_length(result$gradient, n_params)
})
