# Tests for run_mxlogit(scale_vars = "sd"): pre-estimation column scaling of
# X and W design matrices with delta-method back-transform so reported
# quantities (coefficients, vcov, se, sigma) are in the user's natural units.
#
# The fit object also exposes three new fields:
#   - scale_vars: character (e.g., "none" or "sd")
#   - sX:         named numeric vector, length K_x (column scales for X)
#   - sW:         named numeric vector, length K_w (column scales for W;
#                  log-normal columns get sW = 1 by design)
# Stored data$X and data$W are in natural (original) units regardless of
# the scaling pathway taken during optimization.

# =============================================================================
# Test 1: Invariance — across all (rc_correlation, rc_mean) combinations.
# Each variant exercises a different Jacobian-block path:
#   rc_mean = TRUE/FALSE         -> mu block present/absent
#   rc_correlation = TRUE/FALSE  -> off-diagonal L params present/absent
#                                   (the diagonal-only branch is the FALSE arm)
# =============================================================================

fit_pair <- function(dt, rc_correlation, rc_mean, scale_choice) {
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)
  common <- list(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    S = 80L,
    rc_mean = rc_mean,
    rc_correlation = rc_correlation,
    use_asc = TRUE,
    control = ctrl
  )
  list(
    none   = suppressMessages(do.call(run_mxlogit, c(common, list(scale_vars = "none")))),
    scaled = suppressMessages(do.call(run_mxlogit, c(common, list(scale_vars = scale_choice))))
  )
}

# The Jacobian back-transform is identical across sd / mad / iqr — only the
# per-column denominator differs. Looping over all three exercises the scale
# selection switch and confirms each robust SD-equivalent scale yields the same
# natural-scale coefficients, vcov, and log-likelihood as the unscaled fit.
for (scale_choice in c("sd", "mad", "iqr")) {
  for (rc_correlation in c(TRUE, FALSE)) {
    for (rc_mean in c(TRUE, FALSE)) {
      label <- sprintf(
        "scale_vars='%s', rc_correlation=%s, rc_mean=%s",
        scale_choice, rc_correlation, rc_mean
      )
      test_that(paste("scale_vars is invariant:", label), {
        skip_on_cran()
        skip_on_ci()
        sim <- simulate_mxl_data(N = 1200L, J = 4L, seed = 2026)
        dt <- data.table::as.data.table(sim$data)

        fits <- fit_pair(
          dt,
          rc_correlation = rc_correlation,
          rc_mean = rc_mean,
          scale_choice = scale_choice
        )
        f1 <- fits$none; f2 <- fits$scaled

        # Both fits must converge for the comparison to be meaningful
        expect_true(f1$convergence > 0)
        expect_true(f2$convergence > 0)

        # Coefficients in natural scale agree across pathways
        expect_lt(max(abs(coef(f1) - coef(f2))), 1e-5)

        # Standard errors back-transformed via delta method agree
        expect_lt(max(abs(f1$se - f2$se), na.rm = TRUE), 1e-5)

        # Log-likelihood is scale-invariant
        expect_lt(abs(f1$loglik - f2$loglik), 1e-8)

        # Full vcov matrix agrees
        expect_lt(max(abs(f1$vcov - f2$vcov), na.rm = TRUE), 1e-5)
      })
    }
  }
}

# =============================================================================
# Test 2: Log-normal pass-through — sW = 1 for log-normal RC columns
# =============================================================================

test_that("scale_vars='sd' passes log-normal RC columns through unchanged", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mxl_data(N = 800L, J = 3L, seed = 7L)
  dt <- data.table::as.data.table(sim$data)
  expected_sd_w2 <- stats::sd(dt$w2)
  expected_sd_x1 <- stats::sd(dt$x1)
  expected_sd_x2 <- stats::sd(dt$x2)

  # rc_dist = c(1L, 0L) -> w1 is log-normal, w2 is normal.
  # We expect a one-line message announcing log-normal pass-through.
  # capture_messages() returns the messages; the fit assignment escapes the
  # block via lexical scoping in the test_that frame.
  fit <- NULL
  msgs <- capture_messages({
    fit <- run_mxlogit(
      data = dt,
      id_col = "id",
      alt_col = "alt",
      choice_col = "choice",
      covariate_cols = c("x1", "x2"),
      random_var_cols = c("w1", "w2"),
      S = 50L,
      rc_dist = c(1L, 0L),
      rc_mean = TRUE,
      rc_correlation = FALSE,
      use_asc = TRUE,
      scale_vars = "sd",
      control = list(maxeval = 100L)
    )
  })
  expect_true(any(grepl("log-normal", msgs)))

  expect_equal(fit$scale_vars, "sd")

  # sW: log-normal column (w1) gets pass-through scale = 1; normal column
  # (w2) gets its sample SD.
  expect_named(fit$sW, c("w1", "w2"))
  expect_equal(unname(fit$sW["w1"]), 1)
  expect_equal(unname(fit$sW["w2"]), expected_sd_w2, tolerance = 1e-10)

  # sX: both fixed-coef columns get their sample SDs.
  expect_named(fit$sX, c("x1", "x2"))
  expect_equal(unname(fit$sX["x1"]), expected_sd_x1, tolerance = 1e-10)
  expect_equal(unname(fit$sX["x2"]), expected_sd_x2, tolerance = 1e-10)
})

# =============================================================================
# Test 3: Natural-scale storage — fit$data holds original-units matrices
# =============================================================================

test_that("scale_vars='sd' stores natural-scale X/W and post-estimation methods work", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mxl_data(N = 600L, J = 3L, seed = 31L)
  dt <- data.table::as.data.table(sim$data)

  # Deliberately rescale x1 so its SD is far from 1; this is the situation
  # scale_vars='sd' is designed to handle.
  dt[, x1 := x1 * 10]
  expected_sd_x1 <- stats::sd(dt$x1)
  expected_sd_w1 <- stats::sd(dt$w1)

  fit <- suppressMessages(run_mxlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    S = 50L,
    rc_mean = TRUE,
    rc_correlation = TRUE,
    use_asc = TRUE,
    scale_vars = "sd",
    control = list(maxeval = 500L)
  ))

  # Stored design matrix is in natural (user-input) units: compare to the
  # input SD directly, not to fit$sX which by construction equals it under
  # this code path (the assertion against the input SD is the real invariant).
  expect_equal(stats::sd(fit$data$X[, "x1"]), expected_sd_x1, tolerance = 1e-10)
  expect_equal(stats::sd(fit$data$W[, "w1"]), expected_sd_w1, tolerance = 1e-10)

  # predict() consumes fit$data, so it must work end-to-end and return
  # finite probabilities in [0, 1].
  preds <- predict(fit, type = "probabilities")
  expect_true(all(is.finite(preds$choice_prob)))
  expect_true(all(preds$choice_prob >= 0 & preds$choice_prob <= 1))

  # elasticities() likewise must consume natural-scale data correctly.
  elast <- elasticities(fit, elast_var = "x1", is_random_coef = FALSE)
  expect_true(is.matrix(elast))
  expect_true(all(is.finite(elast)))
})

# =============================================================================
# Test 4: Natural-scale storage holds even for log-normal RC columns —
# combines the pass-through carve-out with the natural-units storage invariant.
# Without an explicit test, a regression could substitute the scaled W back in
# for the stored matrix and the existing tests would not catch it for the
# log-normal column (since sW[k] == 1 there masks the difference).
# =============================================================================

test_that("scale_vars='sd' stores log-normal W column in natural units", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mxl_data(N = 500L, J = 3L, seed = 19L)
  dt <- data.table::as.data.table(sim$data)
  # Force w1 onto a non-unit scale so a regression that stored the *scaled*
  # matrix would produce sd(stored$W[,"w1"]) far from sd(dt$w1).
  dt[, w2 := w2 * 5]  # w2 is the normal RC column (rc_dist[2] = 0)
  expected_sd_w1 <- stats::sd(dt$w1)
  expected_sd_w2 <- stats::sd(dt$w2)

  fit <- suppressMessages(run_mxlogit(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = c("w1", "w2"),
    S = 50L,
    rc_dist = c(1L, 0L),  # w1 log-normal, w2 normal
    rc_mean = TRUE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    scale_vars = "sd",
    control = list(maxeval = 300L)
  ))

  # w1 is log-normal so sW["w1"] == 1 by carve-out; storage must still be
  # in natural units (sd matches dt$w1 directly, not 1).
  expect_equal(unname(fit$sW[["w1"]]), 1)
  expect_equal(stats::sd(fit$data$W[, "w1"]), expected_sd_w1, tolerance = 1e-10)

  # w2 is normal, sW["w2"] != 1; storage still in natural units.
  expect_gt(unname(fit$sW[["w2"]]), 1)
  expect_equal(stats::sd(fit$data$W[, "w2"]), expected_sd_w2, tolerance = 1e-10)
})

# =============================================================================
# Test 5: theta_init is interpreted in natural-scale units regardless of the
# scaling mode. Warm-starting from a converged fit's coefficients should reach
# the same optimum quickly under both scale_vars settings.
# =============================================================================

test_that("theta_init is natural-scale under scale_vars='sd' (warm-start works)", {
  skip_on_cran()
  skip_on_ci()
  sim <- simulate_mxl_data(N = 900L, J = 4L, seed = 41L)
  dt <- data.table::as.data.table(sim$data)
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)

  # Warm-up fit produces natural-scale coefficients.
  warmup <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 60L, rc_mean = TRUE, rc_correlation = TRUE, use_asc = TRUE,
    scale_vars = "none", control = ctrl
  ))
  theta_natural <- unname(coef(warmup))

  # Re-fit with the natural-scale theta_init under both modes. If the forward
  # transform inside run_mxlogit is correct, both restarts should converge to
  # the same point as the warm-up.
  fit_none <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 60L, rc_mean = TRUE, rc_correlation = TRUE, use_asc = TRUE,
    scale_vars = "none", theta_init = theta_natural, control = ctrl
  ))
  fit_sd <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 60L, rc_mean = TRUE, rc_correlation = TRUE, use_asc = TRUE,
    scale_vars = "sd",   theta_init = theta_natural, control = ctrl
  ))

  expect_lt(max(abs(coef(fit_sd) - coef(warmup))), 1e-5)
  expect_lt(max(abs(coef(fit_sd) - coef(fit_none))), 1e-5)
  expect_lt(abs(fit_sd$loglik - warmup$loglik), 1e-8)
})
