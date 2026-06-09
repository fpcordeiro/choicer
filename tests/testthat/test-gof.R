# Tests for gof(): McFadden pseudo R-squared and hit rate

# =============================================================================
# Equal-shares null on balanced, unweighted MNL
# =============================================================================

test_that("gof equal-shares null is exact on balanced unweighted MNL", {
  dt <- create_small_mnl_data()  # N = 20, J = 3
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )

  g <- gof(fit)

  expect_s3_class(g, "choicer_gof")
  expect_identical(g$null, "equal_shares")
  expect_equal(g$loglik_null, -20 * log(3), tolerance = TOL_LOGLIK)
  expect_equal(g$loglik, fit$loglik)
  expect_equal(g$mcfadden_r2, 1 - g$loglik / g$loglik_null,
               tolerance = TOL_LOGLIK)
  expect_equal(g$mcfadden_r2_adj,
               1 - (g$loglik - g$n_params) / g$loglik_null,
               tolerance = TOL_LOGLIK)
  # The equal-shares null (beta = 0, ASC = 0) is nested in the model
  expect_gte(g$mcfadden_r2, 0)
  expect_lt(g$mcfadden_r2, 1)
  expect_lt(g$mcfadden_r2_adj, g$mcfadden_r2)
  expect_identical(g$nobs, fit$nobs)
  expect_identical(g$n_params, fit$n_params)
})

# =============================================================================
# Hit rate
# =============================================================================

test_that("gof hit rate equals an explicit per-individual loop", {
  sim <- simulate_mnl_data(
    N = 300, J = 3, beta = c(0.8, -0.6), seed = 99,
    outside_option = FALSE, vary_choice_set = FALSE
  )
  fit <- run_mnlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )

  g <- gof(fit)

  d <- fit$data
  p <- predict(fit, type = "probabilities")$choice_prob
  ends <- cumsum(d$M)
  starts <- ends - d$M + 1L
  hits <- numeric(length(d$M))
  for (i in seq_along(d$M)) {
    pb <- p[starts[i]:ends[i]]
    hits[i] <- as.numeric(which.max(pb) == d$choice_idx[i])
  }
  expect_equal(g$hit_rate, sum(d$weights * hits) / sum(d$weights),
               tolerance = TOL_LOGLIK)
})

test_that("gof hit rate is high on a strongly separated DGP", {
  sim <- simulate_mnl_data(
    N = 500, J = 3, beta = c(6, -6), seed = 7,
    outside_option = FALSE, vary_choice_set = FALSE
  )
  fit <- run_mnlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )

  g <- gof(fit)

  expect_gt(g$hit_rate, 0.6)
  expect_lte(g$hit_rate, 1)
  expect_gt(g$mcfadden_r2, 0.3)
})

# =============================================================================
# Weighted fit
# =============================================================================

test_that("gof equal-shares null matches the manual weighted formula", {
  dt <- create_small_mnl_data()  # ids 1..20 in ascending order, J = 3
  set.seed(1)
  w <- runif(20, 0.5, 2)
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), weights = w
  )

  g <- gof(fit)

  expect_equal(g$loglik_null, -sum(fit$data$weights * log(fit$data$M)),
               tolerance = TOL_LOGLIK)
  expect_equal(g$loglik_null, -sum(w * log(3)), tolerance = TOL_LOGLIK)
  expect_true(is.finite(g$hit_rate))
})

# =============================================================================
# Market-shares null
# =============================================================================

test_that("gof market-shares null matches the manual formula (balanced)", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )

  g <- gof(fit, null = "market_shares")

  am <- fit$alt_mapping
  keep <- am$N_CHOICES > 0
  expect_identical(g$null, "market_shares")
  expect_equal(g$loglik_null,
               sum(am$N_CHOICES[keep] * log(am$MKT_SHARE[keep])),
               tolerance = TOL_LOGLIK)
  # The market-shares null fits at least as well as equal shares
  g_eq <- gof(fit, null = "equal_shares")
  expect_gte(g$loglik_null, g_eq$loglik_null - TOL_LOGLIK)
})

test_that("gof market-shares null errors on varying choice sets", {
  sim <- simulate_mnl_data(
    N = 200, J = 4, seed = 5,
    outside_option = FALSE, vary_choice_set = TRUE
  )
  fit <- run_mnlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )

  expect_error(gof(fit, null = "market_shares"), "balanced")
  expect_error(gof(fit, null = "market_shares"), "ASC-only")

  # equal_shares remains exact for the unbalanced design
  g <- gof(fit)
  expect_equal(g$loglik_null, -sum(log(fit$data$M)), tolerance = TOL_LOGLIK)
})

test_that("gof market-shares null errors with non-uniform weights", {
  dt <- create_small_mnl_data()
  set.seed(2)
  w <- runif(20, 0.5, 2)
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), weights = w
  )

  expect_error(gof(fit, null = "market_shares"), "uniform weights")
})

# =============================================================================
# keep_data = FALSE
# =============================================================================

test_that("gof without stored data returns NA fields with a message", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), keep_data = FALSE
  )

  expect_message(g <- gof(fit), "keep_data = TRUE")

  expect_s3_class(g, "choicer_gof")
  expect_equal(g$loglik, fit$loglik)
  expect_true(is.na(g$loglik_null))
  expect_true(is.na(g$mcfadden_r2))
  expect_true(is.na(g$mcfadden_r2_adj))
  expect_true(is.na(g$hit_rate))

  # print must not error on NA fields
  expect_output(print(g), "Goodness of fit")
  expect_output(print(g), "keep_data")

  # footer helper no-ops gracefully on NA fields
  expect_silent(print_gof_lines(g))
})

# =============================================================================
# include_outside_option
# =============================================================================

test_that("gof with outside option uses M + 1 and handles outside wins", {
  set.seed(11)
  N <- 80
  J <- 4  # alternatives 0..3; j = 0 is the outside option
  dt <- data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J - 1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  # Make the outside option dominant (~70% of choices) so its predicted
  # probability is the largest for most individuals.
  dt[, choice := {
    pick <- sample(c(rep(0L, 7), 1L, 2L, 3L), 1)
    as.integer(j == pick)
  }, by = id]

  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "j", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    outside_opt_label = 0L, include_outside_option = TRUE
  )

  g <- gof(fit)

  # M counts inside alternatives only (3); the outside good adds + 1
  expect_equal(unique(fit$data$M), 3L)
  expect_equal(g$loglik_null, -N * log(4), tolerance = TOL_LOGLIK)

  # Manual hit-rate loop with outside probability = 1 - sum(p_inside)
  d <- fit$data
  p <- predict(fit, type = "probabilities")$choice_prob
  ends <- cumsum(d$M)
  starts <- ends - d$M + 1L
  preds_i <- integer(length(d$M))
  hits <- numeric(length(d$M))
  for (i in seq_along(d$M)) {
    pb <- p[starts[i]:ends[i]]
    po <- 1 - sum(pb)
    preds_i[i] <- if (po > max(pb)) 0L else which.max(pb)
    hits[i] <- as.numeric(preds_i[i] == d$choice_idx[i])
  }
  expect_equal(g$hit_rate, sum(d$weights * hits) / sum(d$weights),
               tolerance = TOL_LOGLIK)

  # The outside-wins branch is actually exercised
  expect_gt(mean(preds_i == 0L), 0.5)
  expect_gt(g$hit_rate, 0.5)

  # market-shares null includes the outside-option row of alt_mapping
  g_ms <- gof(fit, null = "market_shares")
  am <- fit$alt_mapping
  keep <- am$N_CHOICES > 0
  expect_equal(g_ms$loglik_null,
               sum(am$N_CHOICES[keep] * log(am$MKT_SHARE[keep])),
               tolerance = TOL_LOGLIK)
})

# =============================================================================
# MXL and NL smoke tests
# =============================================================================

test_that("gof returns finite numbers for a fitted MXL", {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 50L)
  )

  g <- gof(fit)

  expect_s3_class(g, "choicer_gof")
  expect_true(is.finite(g$loglik_null))
  expect_equal(g$loglik_null, -30 * log(3), tolerance = TOL_LOGLIK)
  expect_true(is.finite(g$mcfadden_r2))
  expect_true(is.finite(g$mcfadden_r2_adj))
  expect_true(g$hit_rate >= 0 && g$hit_rate <= 1)
})

test_that("gof returns finite numbers for a fitted NL", {
  dt <- create_small_nl_data()
  fit <- run_nestlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), nest_col = "nest",
    control = list(maxeval = 50L)
  )

  g <- gof(fit)

  expect_s3_class(g, "choicer_gof")
  expect_true(is.finite(g$loglik_null))
  expect_equal(g$loglik_null, -30 * log(6), tolerance = TOL_LOGLIK)
  expect_true(is.finite(g$mcfadden_r2))
  expect_true(g$hit_rate >= 0 && g$hit_rate <= 1)
})

# =============================================================================
# print methods
# =============================================================================

test_that("print.choicer_gof and print_gof_lines print fitted measures", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2")
  )
  g <- gof(fit)

  expect_output(print(g), "Goodness of fit")
  expect_output(print(g), "McFadden R2")
  expect_output(print(g), "equal shares")
  expect_output(print(g), "Hit rate")

  expect_output(print_gof_lines(g), "McFadden R2: \\d\\.\\d{3} \\(adj: ")
  expect_output(print_gof_lines(g), "\\| Hit rate: \\d\\.\\d{3}")

  capture.output(res <- withVisible(print(g)))
  expect_false(res$visible)
})
