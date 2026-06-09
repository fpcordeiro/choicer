# Tests for wtp(): willingness-to-pay ratios with delta-method standard errors

# =============================================================================
# Fixtures (memoized so each fit runs at most once per test file)
# =============================================================================

.wtp_fits <- new.env(parent = emptyenv())

# MNL on simulated data with known betas; x2 plays the role of price.
get_fit_wtp_mnl <- function() {
  if (is.null(.wtp_fits$mnl)) {
    sim <- simulate_mnl_data(
      N = 2000, J = 4, beta = c(0.8, -0.6), seed = 42,
      outside_option = FALSE, vary_choice_set = FALSE
    )
    .wtp_fits$mnl <- run_mnlogit(
      data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2")
    )
  }
  .wtp_fits$mnl
}

# MXL with fixed price (x1 in X) and one normal random coefficient (w1),
# rc_mean = TRUE so the mu block exists.
get_fit_wtp_mxl_normal <- function() {
  if (is.null(.wtp_fits$mxl_normal)) {
    dt <- create_identified_mxl_data(seed = 101, N = 200)
    .wtp_fits$mxl_normal <- run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = "w1",
      rc_mean = TRUE, S = 50L
    )
  }
  .wtp_fits$mxl_normal
}

# Same data, rc_mean = FALSE: no mu block.
get_fit_wtp_mxl_nomean <- function() {
  if (is.null(.wtp_fits$mxl_nomean)) {
    dt <- create_identified_mxl_data(seed = 101, N = 200)
    .wtp_fits$mxl_nomean <- run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = "w1",
      rc_mean = FALSE, S = 50L
    )
  }
  .wtp_fits$mxl_nomean
}

# MXL with a log-normal random coefficient and a fixed price variable (x2).
get_fit_wtp_mxl_lognormal <- function() {
  if (is.null(.wtp_fits$mxl_lognormal)) {
    sim <- simulate_mxl_data(
      N = 400, J = 3, beta = c(0.5, -1.0), mu = 0.3,
      Sigma = matrix(0.25, 1, 1), rc_dist = 1L, seed = 7,
      outside_option = FALSE, vary_choice_set = FALSE
    )
    .wtp_fits$mxl_lognormal <- run_mxlogit(
      data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = "w1",
      rc_dist = 1L, rc_mean = TRUE, S = 50L
    )
  }
  .wtp_fits$mxl_lognormal
}

get_fit_wtp_nl <- function() {
  if (is.null(.wtp_fits$nl)) {
    dt <- create_small_nl_data()
    .wtp_fits$nl <- run_nestlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), nest_col = "nest"
    )
  }
  .wtp_fits$nl
}

# =============================================================================
# MNL: exact ratio identity, MC recovery, structure
# =============================================================================

test_that("wtp.choicer_fit returns -coef[k]/coef[p] exactly (MNL)", {
  fit <- get_fit_wtp_mnl()
  cf <- coef(fit)

  w <- wtp(fit, price_var = "x2")

  expect_s3_class(w, "choicer_wtp")
  expect_s3_class(w, "data.frame")
  expect_equal(
    colnames(w),
    c("Estimate", "Std_Error", "z_value", "CI_lower", "CI_upper")
  )
  expect_true("x1" %in% rownames(w))
  expect_equal(w["x1", "Estimate"], unname(-cf["x1"] / cf["x2"]),
               tolerance = 1e-12)
  expect_identical(attr(w, "price_var"), "x2")
  expect_identical(attr(w, "level"), 0.95)
})

test_that("wtp recovers the true ratio on a known DGP (MNL)", {
  fit <- get_fit_wtp_mnl()
  w <- wtp(fit, price_var = "x2")

  # True WTP = -beta_x1 / beta_x2 = -0.8 / (-0.6) = 4/3
  true_wtp <- -0.8 / (-0.6)
  expect_lt(abs(w["x1", "Estimate"] - true_wtp), 0.2)
  expect_true(is.finite(w["x1", "Std_Error"]))
  # CI should cover the estimate
  expect_lt(w["x1", "CI_lower"], w["x1", "Estimate"])
  expect_gt(w["x1", "CI_upper"], w["x1", "Estimate"])
})

test_that("wtp of an ASC equals -asc/beta_price", {
  fit <- get_fit_wtp_mnl()
  cf <- coef(fit)

  w <- wtp(fit, price_var = "x2", attr_vars = c("x1", "ASC_2"))

  expect_equal(w["ASC_2", "Estimate"], unname(-cf["ASC_2"] / cf["x2"]),
               tolerance = 1e-12)
})

# =============================================================================
# Delta-method SE vs numDeriv
# =============================================================================

test_that("wtp delta-method SE matches numDeriv jacobian (MNL)", {
  skip_if_not_installed("numDeriv")
  fit <- get_fit_wtp_mnl()
  th <- coef(fit)
  V <- vcov(fit)
  skip_if(is.null(V), "vcov unavailable")

  k <- which(names(th) == "x1")
  p <- which(names(th) == "x2")
  g <- function(t) -t[k] / t[p]
  Jg <- numDeriv::jacobian(g, th)
  se_num <- sqrt(as.numeric(Jg %*% V %*% t(Jg)))

  w <- wtp(fit, price_var = "x2")
  expect_equal(w["x1", "Std_Error"], se_num, tolerance = TOL_GRAD)
})

test_that("wtp log-normal median WTP and SE match numDeriv (MXL)", {
  skip_if_not_installed("numDeriv")
  fit <- get_fit_wtp_mxl_lognormal()
  th <- coef(fit)
  V <- vcov(fit)
  skip_if(is.null(V), "vcov unavailable")

  mu_idx <- fit$param_map$mu[1]
  p_idx <- which(names(th) == "x2")

  w <- wtp(fit, price_var = "x2")
  lbl <- "w1"
  expect_true(lbl %in% rownames(w))
  expect_true(lbl %in% attr(w, "median_rows"))

  # Median WTP point estimate for the shifted log-normal
  # beta = exp(mu) + exp(L eta): -(exp(mu) + 1)/beta_p
  expect_equal(w[lbl, "Estimate"],
               unname(-(exp(th[mu_idx]) + 1) / th[p_idx]),
               tolerance = 1e-12)

  # Delta-method SE vs full numerical jacobian
  g <- function(t) -(exp(t[mu_idx]) + 1) / t[p_idx]
  Jg <- numDeriv::jacobian(g, th)
  se_num <- sqrt(as.numeric(Jg %*% V %*% t(Jg)))
  expect_equal(w[lbl, "Std_Error"], se_num, tolerance = TOL_GRAD)
})

test_that("wtp log-normal RC with rc_mean = FALSE has median WTP -1/beta_p", {
  skip_if_not_installed("numDeriv")
  sim <- simulate_mxl_data(
    N = 400, J = 3, beta = c(0.5, -1.0), mu = 0,
    Sigma = matrix(0.25, 1, 1), rc_dist = 1L, seed = 11,
    outside_option = FALSE, vary_choice_set = FALSE
  )
  fit <- run_mxlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = "w1",
    rc_dist = 1L, rc_mean = FALSE, S = 50L
  )
  th <- coef(fit)
  p_idx <- which(names(th) == "x2")

  w <- wtp(fit, price_var = "x2")
  expect_true("w1" %in% rownames(w))
  expect_true("w1" %in% attr(w, "median_rows"))
  # beta = exp(L eta) has median 1: median WTP = -1/beta_p
  expect_equal(w["w1", "Estimate"], unname(-1 / th[p_idx]),
               tolerance = 1e-12)

  V <- vcov(fit)
  if (!is.null(V) && is.finite(V[p_idx, p_idx])) {
    g <- function(t) -1 / t[p_idx]
    Jg <- numDeriv::jacobian(g, th)
    se_num <- sqrt(as.numeric(Jg %*% V %*% t(Jg)))
    expect_equal(w["w1", "Std_Error"], se_num, tolerance = TOL_GRAD)
  }
})

# =============================================================================
# MXL: mu-block handling and random-price guard
# =============================================================================

test_that("wtp includes normal-RC mean row -mu/beta_p by default (MXL)", {
  fit <- get_fit_wtp_mxl_normal()
  th <- coef(fit)

  w <- wtp(fit, price_var = "x1")

  expect_true("Mu_w1" %in% rownames(w))
  expect_equal(
    w["Mu_w1", "Estimate"],
    unname(-th[fit$param_map$mu[1]] / th[fit$param_map$beta[1]]),
    tolerance = 1e-12
  )
  # Normal RC rows are mean WTP, not median
  expect_false("Mu_w1" %in% attr(w, "median_rows"))
})

test_that("wtp errors on a random price coefficient (MXL)", {
  fit <- get_fit_wtp_mxl_normal()

  expect_error(wtp(fit, price_var = "w1"), "[Rr]andom price")
  expect_error(wtp(fit, price_var = "w1"), "WTP space")
})

test_that("wtp excludes normal RCs when rc_mean = FALSE (MXL)", {
  fit <- get_fit_wtp_mxl_nomean()

  w <- wtp(fit, price_var = "x1")

  # x1 is the only fixed coefficient and w1 is a mean-zero normal RC,
  # so no default attribute rows remain
  expect_equal(nrow(w), 0L)
  expect_false(any(grepl("w1", rownames(w))))
  expect_match(attr(w, "rc_note"), "rc_mean = FALSE")

  # ASCs can still be requested explicitly
  cf <- coef(fit)
  w2 <- wtp(fit, price_var = "x1", attr_vars = "ASC_2")
  expect_equal(w2["ASC_2", "Estimate"], unname(-cf["ASC_2"] / cf["x1"]),
               tolerance = 1e-12)

  # Asking for the random coefficient by name is an informative error
  expect_error(wtp(fit, price_var = "x1", attr_vars = "w1"), "rc_mean")
})

# =============================================================================
# NL: dispatches through the choicer_fit method
# =============================================================================

test_that("wtp works for nested logit via the choicer_fit method", {
  fit <- get_fit_wtp_nl()
  cf <- coef(fit)

  w <- wtp(fit, price_var = "x1")

  expect_s3_class(w, "choicer_wtp")
  expect_true("x2" %in% rownames(w))
  expect_equal(w["x2", "Estimate"], unname(-cf["x2"] / cf["x1"]),
               tolerance = 1e-12)
  # Lambda (nest) parameters are not WTP attributes
  expect_false(any(grepl("^Lambda", rownames(w))))
})

# =============================================================================
# Confidence level and input validation
# =============================================================================

test_that("wtp level = 0.90 gives a narrower CI than 0.95", {
  fit <- get_fit_wtp_mnl()

  w95 <- wtp(fit, price_var = "x2", level = 0.95)
  w90 <- wtp(fit, price_var = "x2", level = 0.90)

  width95 <- w95["x1", "CI_upper"] - w95["x1", "CI_lower"]
  width90 <- w90["x1", "CI_upper"] - w90["x1", "CI_lower"]
  expect_lt(width90, width95)
  expect_identical(attr(w90, "level"), 0.90)
})

test_that("wtp validates price_var, attr_vars, and level", {
  fit <- get_fit_wtp_mnl()

  expect_error(wtp(fit, price_var = "nonexistent"), "not found")
  expect_error(wtp(fit, price_var = "ASC_2"), "not found")
  expect_error(wtp(fit, price_var = "x2", attr_vars = "nonexistent"),
               "not found")
  expect_error(wtp(fit, price_var = "x2", attr_vars = "x2"),
               "must not include the price variable")
  expect_error(wtp(fit, price_var = "x2", level = 1.5), "level")
  expect_error(wtp(fit, price_var = "x2", level = 0), "level")
})

# =============================================================================
# print method
# =============================================================================

test_that("print.choicer_wtp prints header with WTP and price_var", {
  fit <- get_fit_wtp_mnl()
  w <- wtp(fit, price_var = "x2")

  expect_output(print(w), "WTP")
  expect_output(print(w), "x2")
  expect_output(print(w), "Estimate")

  # print returns its argument invisibly
  capture.output(res <- withVisible(print(w)))
  expect_false(res$visible)
  expect_identical(res$value, w)
})

test_that("print.choicer_wtp flags median rows for log-normal RCs", {
  fit <- get_fit_wtp_mxl_lognormal()
  w <- wtp(fit, price_var = "x2")

  expect_output(print(w), "median WTP")
  expect_output(print(w), "w1")
})
