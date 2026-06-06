# Tests for run_mnlogit(scale_vars = ...): pre-estimation column scaling of the
# fixed-coefficient design matrix X with a delta-method back-transform, so
# reported quantities (coefficients, vcov, se, loglik) are in natural units.
#
# The fit object also exposes two new fields:
#   - scale_vars: character (e.g., "none", "sd", "mad", "iqr")
#   - sX:         named numeric vector, length K_x (column scales for X)
# Stored data$X is in natural (original) units regardless of the scaling path.

# =============================================================================
# Test 1: Invariance across sd / mad / iqr — natural-scale coefficients, SEs,
# vcov, and log-likelihood agree with the unscaled fit.
# =============================================================================

fit_mnl_pair <- function(dt, scale_choice) {
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)
  common <- list(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    use_asc = TRUE,
    control = ctrl
  )
  list(
    none   = suppressMessages(do.call(run_mnlogit, c(common, list(scale_vars = "none")))),
    scaled = suppressMessages(do.call(run_mnlogit, c(common, list(scale_vars = scale_choice))))
  )
}

for (scale_choice in c("sd", "mad", "iqr")) {
  test_that(paste0("scale_vars is invariant: scale_vars='", scale_choice, "'"), {
    skip_on_cran()
    skip_on_ci()
    sim <- simulate_mnl_data(N = 1500L, J = 4L, seed = 2026L,
                             outside_option = FALSE)
    dt <- data.table::as.data.table(sim$data)

    fits <- fit_mnl_pair(dt, scale_choice)
    f1 <- fits$none; f2 <- fits$scaled

    expect_true(f1$convergence > 0)
    expect_true(f2$convergence > 0)

    expect_lt(max(abs(coef(f1) - coef(f2))), 1e-5)
    expect_false(anyNA(f1$se))
    expect_false(anyNA(f2$se))
    expect_lt(max(abs(f1$se - f2$se), na.rm = TRUE), 1e-5)
    expect_lt(abs(f1$loglik - f2$loglik), 1e-8)
    expect_lt(max(abs(f1$vcov - f2$vcov), na.rm = TRUE), 1e-5)
  })
}

# =============================================================================
# Test 2: Fields — scale_vars, sX, and natural-scale storage of data$X.
# =============================================================================

test_that("scale_vars='sd' records sX and stores natural-scale X", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mnl_data(N = 800L, J = 3L, seed = 17L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  # Rescale x1 so its SD is far from 1 (the situation scaling addresses).
  dt[, x1 := x1 * 10]
  expected_sd_x1 <- stats::sd(dt$x1)
  expected_sd_x2 <- stats::sd(dt$x2)

  fit <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), use_asc = TRUE, scale_vars = "sd",
    control = list(maxeval = 2000L)
  ))

  expect_equal(fit$scale_vars, "sd")
  expect_named(fit$sX, c("x1", "x2"))
  expect_equal(unname(fit$sX["x1"]), expected_sd_x1, tolerance = 1e-10)
  expect_equal(unname(fit$sX["x2"]), expected_sd_x2, tolerance = 1e-10)

  # Stored design matrix is in natural (user-input) units.
  expect_equal(stats::sd(fit$data$X[, "x1"]), expected_sd_x1, tolerance = 1e-10)
  expect_equal(stats::sd(fit$data$X[, "x2"]), expected_sd_x2, tolerance = 1e-10)
})

# =============================================================================
# Test 3: Post-estimation — predict(), summary(), elasticities() agree between
# scaled and unscaled fits.
# =============================================================================

test_that("scale_vars='sd' post-estimation matches unscaled fit", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mnl_data(N = 1000L, J = 4L, seed = 31L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)

  f_none <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), use_asc = TRUE, scale_vars = "none",
    control = ctrl
  ))
  f_sd <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), use_asc = TRUE, scale_vars = "sd",
    control = ctrl
  ))

  p_none <- predict(f_none, type = "shares")
  p_sd   <- predict(f_sd, type = "shares")
  expect_equal(p_none, p_sd, tolerance = 1e-5)

  # summary() must run end-to-end on a scaled fit.
  expect_silent(s <- summary(f_sd))
  expect_s3_class(s, "summary.choicer_mnl")

  e_none <- elasticities(f_none, elast_var = "x1")
  e_sd   <- elasticities(f_sd, elast_var = "x1")
  expect_equal(e_none, e_sd, tolerance = 1e-5)
})

# =============================================================================
# Test 4: Error path — a (near-)constant covariate triggers the scale error.
# =============================================================================

test_that("scale_vars='sd' errors on a near-constant covariate", {
  skip_on_cran()
  sim <- simulate_mnl_data(N = 400L, J = 3L, seed = 5L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  dt[, x2 := 1]  # constant column -> SD == 0

  expect_error(
    suppressMessages(run_mnlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), use_asc = TRUE, scale_vars = "sd",
      control = list(maxeval = 100L)
    )),
    "scale <"
  )
})

# =============================================================================
# Test 5: Lazy-recompute vcov parity — the eagerly-stored (scaled-then-back-
# transformed) vcov matches a fresh lazy recompute from the stored natural-scale
# data via compute_hessian()/invert_hessian().
# =============================================================================

test_that("scale_vars='sd' eager vcov matches lazy recompute from stored data", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mnl_data(N = 1500L, J = 4L, seed = 2026L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)

  fit_sd <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), use_asc = TRUE, scale_vars = "sd",
    keep_data = TRUE,
    control = list(xtol_rel = 1e-12, maxeval = 5000L)
  ))

  recomputed <- choicer:::invert_hessian(
    choicer:::compute_hessian(fit_sd)
  )$vcov

  expect_equal(unname(fit_sd$vcov), unname(recomputed), tolerance = 1e-6)
})

# =============================================================================
# Test 6: Advanced pathway — pre-built input_data + scale_vars is invariant
# (coef / se / loglik) against the unscaled fit on the same prepared data.
# =============================================================================

test_that("scale_vars is invariant via the advanced input_data pathway", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mnl_data(N = 1500L, J = 4L, seed = 2026L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)

  prep <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  f_none <- suppressMessages(run_mnlogit(
    input_data = prep, use_asc = TRUE, scale_vars = "none", control = ctrl
  ))
  f_sd <- suppressMessages(run_mnlogit(
    input_data = prep, use_asc = TRUE, scale_vars = "sd", control = ctrl
  ))

  expect_true(f_none$convergence > 0)
  expect_true(f_sd$convergence > 0)

  expect_lt(max(abs(coef(f_none) - coef(f_sd))), 1e-5)
  expect_false(anyNA(f_none$se))
  expect_false(anyNA(f_sd$se))
  expect_lt(max(abs(f_none$se - f_sd$se)), 1e-5)
  expect_lt(abs(f_none$loglik - f_sd$loglik), 1e-8)
})
