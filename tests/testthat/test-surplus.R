# Tests for logsum() and consumer_surplus()
# - MNL logsum vs manual log-sum-exp (with and without outside option)
# - MNL mean-CS delta-method SE vs numDeriv
# - NL with all lambda = 1 collapses to the plain MNL logsum
# - mxl_logsum kernel vs brute-force R simulation (normal and shifted
#   log-normal random coefficients)
# - Jensen sanity: simulated E[logsum] >= logsum of draw-averaged utilities
# - consumer_surplus identities, newdata round trip, direction, validation
# - print.choicer_cs smoke test

# =============================================================================
# Helpers
# =============================================================================

# Reference blockwise log-sum-exp in plain R
manual_logsum <- function(V, M, ioo = FALSE) {
  ends <- cumsum(M)
  starts <- ends - M + 1L
  vapply(seq_along(M), function(i) {
    v <- V[starts[i]:ends[i]]
    if (ioo) v <- c(v, 0)
    log(sum(exp(v)))
  }, numeric(1))
}

# MNL dataset with an explicit outside option (alt = 0, zero covariates);
# mirrors the fixture used in test-predict-newdata.R.
create_mnl_outside_data <- function(seed = 7, N = 40, J_inside = 3) {
  set.seed(seed)
  dt <- data.table(
    id = rep(1:N, each = J_inside + 1),
    alt = rep(0:J_inside, N),
    x1 = rnorm(N * (J_inside + 1)),
    x2 = runif(N * (J_inside + 1), -1, 1)
  )
  dt[alt == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := {
    pick <- sample.int(J_inside + 1, 1) - 1L
    as.integer(alt == pick)
  }, by = id]
  dt[]
}

# Memoized fixtures so each fit runs at most once per test file
.surplus_fits <- new.env(parent = emptyenv())

get_fit_surplus_mnl <- function() {
  if (is.null(.surplus_fits$mnl)) {
    dt <- create_small_mnl_data()
    .surplus_fits$mnl <- run_mnlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2")
    )
  }
  .surplus_fits$mnl
}

get_fit_surplus_mnl_outside <- function() {
  if (is.null(.surplus_fits$mnl_outside)) {
    .surplus_fits$mnl_outside <- run_mnlogit(
      data = create_mnl_outside_data(), id_col = "id", alt_col = "alt",
      choice_col = "choice", covariate_cols = c("x1", "x2"),
      outside_opt_label = 0L, include_outside_option = TRUE
    )
  }
  .surplus_fits$mnl_outside
}

# MNL on simulated data with known betas (x2 plays the role of price, with a
# negative coefficient so -alpha > 0); used for the direction tests.
get_fit_surplus_mnl_sim <- function() {
  if (is.null(.surplus_fits$mnl_sim)) {
    sim <- simulate_mnl_data(
      N = 500, J = 3, beta = c(0.8, -0.6), seed = 42,
      outside_option = FALSE, vary_choice_set = FALSE
    )
    .surplus_fits$mnl_sim <- list(
      fit = run_mnlogit(
        data = sim$data, id_col = "id", alt_col = "alt",
        choice_col = "choice", covariate_cols = c("x1", "x2")
      ),
      dt = sim$data
    )
  }
  .surplus_fits$mnl_sim
}

get_fit_surplus_mxl <- function() {
  if (is.null(.surplus_fits$mxl)) {
    dt <- create_small_mxl_data()
    .surplus_fits$mxl <- run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = c("w1", "w2"),
      S = 10L, control = list(maxeval = 50L)
    )
  }
  .surplus_fits$mxl
}

get_fit_surplus_nl <- function() {
  if (is.null(.surplus_fits$nl)) {
    dt <- create_small_nl_data()
    .surplus_fits$nl <- run_nestlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), nest_col = "nest"
    )
  }
  .surplus_fits$nl
}

# =============================================================================
# 1. MNL logsum vs manual R log-sum-exp
# =============================================================================

test_that("MNL logsum matches the manual log(sum(exp(V))) per block", {
  fit <- get_fit_surplus_mnl()
  V <- predict(fit, type = "probabilities")$utility

  expect_equal(logsum(fit), manual_logsum(V, fit$data$M, ioo = FALSE),
               tolerance = TOL_LOGLIK)
})

test_that("MNL logsum with outside option includes the exp(0) = 1 term", {
  fit <- get_fit_surplus_mnl_outside()
  V <- predict(fit, type = "probabilities")$utility
  M <- fit$data$M

  ls <- logsum(fit)
  expect_equal(ls, manual_logsum(V, M, ioo = TRUE), tolerance = TOL_LOGLIK)
  # The +1 term matters: dropping it gives a strictly smaller logsum
  expect_true(all(ls > manual_logsum(V, M, ioo = FALSE)))
})

test_that("logsum errors without stored data and mirrors predict's message", {
  fit <- get_fit_surplus_mnl()
  fit_slim <- fit
  fit_slim$data <- NULL
  expect_error(logsum(fit_slim), "keep_data")
})

# =============================================================================
# 2. MNL mean-CS delta-method SE vs numDeriv
# =============================================================================

# Weighted mean CS as a pure-R function of theta (recomputes utilities,
# logsums, and CS from the stored design matrix).
cs_mean_fn <- function(theta, fit) {
  d <- fit$data
  pm <- fit$param_map
  ioo <- fit$include_outside_option
  V <- as.numeric(d$X %*% theta[pm$beta])
  if (fit$use_asc && !is.null(pm$asc)) {
    asc_full <- if (ioo) theta[pm$asc] else c(0, theta[pm$asc])
    V <- V + asc_full[d$alt_idx]
  }
  ls <- manual_logsum(V, d$M, ioo = ioo)
  price_idx <- match("x2", names(fit$coefficients))
  cs <- ls / (-theta[price_idx])
  sum(d$weights * cs) / sum(d$weights)
}

test_that("MNL mean-CS delta SE matches numDeriv (fit with ASCs)", {
  skip_if_not_installed("numDeriv")
  fit <- get_fit_surplus_mnl()
  expect_true(!is.null(fit$param_map$asc))  # ASC gradient slots exercised

  out <- consumer_surplus(fit, price_var = "x2")

  g <- numDeriv::grad(cs_mean_fn, coef(fit), fit = fit)
  se_num <- sqrt(as.numeric(t(g) %*% vcov(fit) %*% g))

  expect_true(is.finite(out$se_mean_cs))
  expect_equal(out$se_mean_cs, se_num, tolerance = TOL_GRAD)
})

test_that("MNL mean-CS delta SE matches numDeriv with an outside option", {
  skip_if_not_installed("numDeriv")
  fit <- get_fit_surplus_mnl_outside()

  out <- consumer_surplus(fit, price_var = "x2")

  g <- numDeriv::grad(cs_mean_fn, coef(fit), fit = fit)
  se_num <- sqrt(as.numeric(t(g) %*% vcov(fit) %*% g))

  expect_true(is.finite(out$se_mean_cs))
  expect_equal(out$se_mean_cs, se_num, tolerance = TOL_GRAD)
})

# =============================================================================
# 3. NL with all lambda = 1 equals the plain MNL logsum
# =============================================================================

test_that("NL logsum with all lambda = 1 collapses to log(sum(exp(V)))", {
  fit <- get_fit_surplus_nl()
  fit1 <- fit
  fit1$coefficients[fit1$param_map$lambda] <- 1

  V <- predict(fit1, type = "probabilities")$utility  # V does not depend on lambda
  expect_equal(logsum(fit1), manual_logsum(V, fit1$data$M, ioo = FALSE),
               tolerance = TOL_LOGLIK)
})

test_that("NL logsum matches a direct R implementation of the nested formula", {
  fit <- get_fit_surplus_nl()
  d <- fit$data
  V <- predict(fit, type = "probabilities")$utility

  n_nests <- max(d$nest_idx)
  lambda_full <- rep(1, n_nests)
  lambda_full[which(tabulate(d$nest_idx, n_nests) > 1)] <-
    coef(fit)[fit$param_map$lambda]
  nests <- d$nest_idx[d$alt_idx]

  ends <- cumsum(d$M)
  starts <- ends - d$M + 1L
  manual <- vapply(seq_along(d$M), function(i) {
    rows <- starts[i]:ends[i]
    terms <- vapply(unique(nests[rows]), function(b) {
      lambda_full[b] * log(sum(exp(V[rows][nests[rows] == b] / lambda_full[b])))
    }, numeric(1))
    log(sum(exp(terms)))
  }, numeric(1))

  expect_equal(logsum(fit), manual, tolerance = TOL_LOGLIK)
})

# =============================================================================
# 4. mxl_logsum kernel vs brute-force R simulation
# =============================================================================

test_that("mxl_logsum matches brute force: normal RCs, diagonal L, ASCs", {
  set.seed(2024)
  N <- 5; J <- 3; S <- 4; K_w <- 2
  M <- rep(J, N)
  alt_idx <- rep(1:J, N)
  X <- matrix(rnorm(N * J), ncol = 1)
  W <- matrix(rnorm(N * J * K_w), ncol = K_w)
  # theta = [beta(1), L_params(2, log-diagonal), asc(J - 1 = 2)]
  theta <- c(0.5, log(0.7), log(0.4), 0.3, -0.2)
  eta <- get_halton_normals(S, N, K_w)

  ls <- as.numeric(mxl_logsum(
    theta, X, W, alt_idx, M, eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  ))

  # Brute force: rc_correlation = FALSE => L = diag(exp(L_params))
  L <- diag(exp(theta[2:3]))
  delta <- c(0, theta[4:5])
  manual <- vapply(seq_len(N), function(i) {
    rows <- ((i - 1) * J + 1):(i * J)
    base <- as.numeric(X[rows, , drop = FALSE] %*% theta[1]) + delta[alt_idx[rows]]
    mean(vapply(seq_len(S), function(s) {
      v <- base + as.numeric(W[rows, ] %*% (L %*% eta[, s, i]))
      log(sum(exp(v)))
    }, numeric(1)))
  }, numeric(1))

  expect_equal(ls, manual, tolerance = 1e-8)
})

test_that("mxl_logsum matches brute force: shifted log-normal RC + outside", {
  set.seed(2025)
  N <- 5; J <- 3; S <- 4; K_w <- 1
  M <- rep(J, N)
  alt_idx <- rep(1:J, N)
  X <- matrix(rnorm(N * J), ncol = 1)
  W <- matrix(runif(N * J, 0, 1), ncol = 1)
  # theta = [beta(1), mu(1), L_params(1)] (rc_mean = TRUE, use_asc = FALSE)
  theta <- c(-0.4, 0.2, log(0.5))
  eta <- get_halton_normals(S, N, K_w)

  ls <- as.numeric(mxl_logsum(
    theta, X, W, alt_idx, M, eta, rc_dist = 1L,
    rc_correlation = FALSE, rc_mean = TRUE,
    use_asc = FALSE, include_outside_option = TRUE
  ))

  # Shifted log-normal: beta_k = exp(mu) + exp(L * eta)
  L <- exp(theta[3])
  mu_final <- exp(theta[2])
  manual <- vapply(seq_len(N), function(i) {
    rows <- ((i - 1) * J + 1):(i * J)
    base <- as.numeric(X[rows, , drop = FALSE] %*% theta[1]) +
      W[rows, 1] * mu_final
    mean(vapply(seq_len(S), function(s) {
      gamma_s <- exp(L * eta[1, s, i])
      v <- c(0, base + W[rows, 1] * gamma_s)  # outside option's exp(0) slot
      log(sum(exp(v)))
    }, numeric(1)))
  }, numeric(1))

  expect_equal(ls, manual, tolerance = 1e-8)
})

# =============================================================================
# 5. Jensen sanity: E[logsum] >= logsum of draw-averaged utilities
# =============================================================================

test_that("MXL logsum dominates the logsum of mxl_predict's averaged utility", {
  fit <- get_fit_surplus_mxl()

  ls_sim <- logsum(fit)
  V_avg <- predict(fit, type = "probabilities")$utility
  ls_naive <- manual_logsum(V_avg, fit$data$M, ioo = FALSE)

  # log-sum-exp is convex, so E_s[lse(V_s)] >= lse(E_s[V_s]) (Jensen);
  # strict for a nondegenerate Sigma.
  expect_true(all(ls_sim >= ls_naive - 1e-10))
  expect_gt(min(ls_sim - ls_naive), 0)
})

# =============================================================================
# 6. consumer_surplus: identities, round trip, direction, validation
# =============================================================================

test_that("cs equals logsum / (-alpha) for all model classes", {
  fit_mnl <- get_fit_surplus_mnl()
  cs_mnl <- consumer_surplus(fit_mnl, price_var = "x2")
  expect_s3_class(cs_mnl, "choicer_cs")
  expect_equal(cs_mnl$cs, logsum(fit_mnl) / (-coef(fit_mnl)[["x2"]]),
               tolerance = 1e-12)
  expect_equal(cs_mnl$mean_cs, mean(cs_mnl$cs), tolerance = 1e-12)
  expect_identical(cs_mnl$n, length(fit_mnl$data$M))

  fit_mxl <- get_fit_surplus_mxl()
  cs_mxl <- consumer_surplus(fit_mxl, price_var = "x1")
  expect_equal(cs_mxl$cs, logsum(fit_mxl) / (-coef(fit_mxl)[["x1"]]),
               tolerance = 1e-12)
  expect_identical(cs_mxl$se_mean_cs, NA_real_)
  expect_identical(cs_mxl$ci, c(NA_real_, NA_real_))

  fit_nl <- get_fit_surplus_nl()
  cs_nl <- consumer_surplus(fit_nl, price_var = "x2")
  expect_equal(cs_nl$cs, logsum(fit_nl) / (-coef(fit_nl)[["x2"]]),
               tolerance = 1e-12)
  expect_identical(cs_nl$se_mean_cs, NA_real_)
  expect_identical(cs_nl$ci, c(NA_real_, NA_real_))
})

test_that("newdata = original data reproduces the stored-data CS", {
  m <- get_fit_surplus_mnl_sim()

  cs0 <- consumer_surplus(m$fit, price_var = "x2")
  cs1 <- consumer_surplus(m$fit, price_var = "x2", newdata = m$dt)

  expect_equal(cs1$cs, cs0$cs, tolerance = 1e-12)
  expect_equal(cs1$mean_cs, cs0$mean_cs, tolerance = 1e-12)
  expect_equal(cs1$se_mean_cs, cs0$se_mean_cs, tolerance = 1e-12)
})

test_that("improving an attribute with a positive coefficient raises mean CS", {
  m <- get_fit_surplus_mnl_sim()
  expect_gt(coef(m$fit)[["x1"]], 0)  # DGP: beta_x1 = 0.8
  expect_lt(coef(m$fit)[["x2"]], 0)  # DGP: beta_x2 = -0.6 (price)

  cs0 <- consumer_surplus(m$fit, price_var = "x2")
  dt_cf <- copy(m$dt)[alt == 2, x1 := x1 + 1]
  cs1 <- consumer_surplus(m$fit, price_var = "x2", newdata = dt_cf)
  expect_gt(cs1$mean_cs, cs0$mean_cs)

  # ... and a price increase lowers it
  dt_price <- copy(m$dt)[alt == 2, x2 := x2 + 1]
  cs2 <- consumer_surplus(m$fit, price_var = "x2", newdata = dt_price)
  expect_lt(cs2$mean_cs, cs0$mean_cs)
})

test_that("price_var validation: unknown variable and random MXL price error", {
  fit_mnl <- get_fit_surplus_mnl()
  expect_error(consumer_surplus(fit_mnl, price_var = "nope"),
               "not found among fixed-coefficient")
  expect_error(consumer_surplus(fit_mnl, price_var = 1L),
               "single variable name")

  fit_mxl <- get_fit_surplus_mxl()
  expect_error(consumer_surplus(fit_mxl, price_var = "w1"),
               "Random price coefficients are not supported")

  expect_error(consumer_surplus(fit_mnl, price_var = "x2", level = 1.2),
               "'level' must be a single number")
})

# =============================================================================
# 7. print.choicer_cs smoke test
# =============================================================================

test_that("print.choicer_cs prints the summary and returns invisibly", {
  fit <- get_fit_surplus_mnl()
  cs <- consumer_surplus(fit, price_var = "x2")

  expect_output(print(cs), "Consumer surplus, price variable: 'x2'")
  expect_output(print(cs), "Mean CS:")
  expect_output(print(cs), "SE \\(delta method\\):")
  expect_output(print(cs), "95% CI:")

  cs_na <- consumer_surplus(get_fit_surplus_nl(), price_var = "x2")
  expect_output(print(cs_na), "SE: NA")

  expect_invisible(print(cs))
})
