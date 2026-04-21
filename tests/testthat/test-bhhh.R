# Tests for BHHH (OPG) information matrix for Mixed Logit

test_that("mxl_bhhh_parallel returns a symmetric matrix with correct dims", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)  # beta + L (diag) + ASCs (J-1)

  set.seed(123)
  theta <- runif(n_params, -0.3, 0.3)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("mxl_bhhh_parallel is positive semi-definite", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)

  set.seed(7)
  theta <- runif(n_params, -0.3, 0.3)
  eta <- get_halton_normals(S = 40, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )

  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-8))
})

test_that("mxl_bhhh_parallel returns PSD under rc_mean = TRUE", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + K_w + (J - 1)  # beta + mu + L (diag) + ASCs

  set.seed(11)
  theta <- runif(n_params, -0.2, 0.2)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = TRUE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_equal(H, t(H), tolerance = 1e-10)
  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-8))
})

test_that("mxl_bhhh_parallel returns PSD under rc_correlation = TRUE", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"),
                        rc_correlation = TRUE)

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  L_size <- K_w * (K_w + 1) / 2
  n_params <- K_x + L_size + (J - 1)  # beta + L (lower-tri) + ASCs

  set.seed(13)
  theta <- runif(n_params, -0.2, 0.2)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = TRUE, rc_mean = FALSE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_equal(H, t(H), tolerance = 1e-10)
  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-8))
})

test_that("mxl_bhhh_parallel returns PSD with include_outside_option = TRUE", {
  # Create a small dataset where some ids have no inside choice (outside option).
  set.seed(17)
  N <- 40; J <- 3
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), N),
    x1 = rnorm(N * J),
    w1 = rnorm(N * J)
  )
  dt[, choice := 0L]
  # ~60% of ids choose an inside alternative, the rest take the outside option
  inside_ids <- sample(seq_len(N), size = 25)
  dt[id %in% inside_ids, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  d <- prepare_mxl_data(
    dt, "id", "alt", "choice",
    covariate_cols = "x1",
    random_var_cols = "w1",
    include_outside_option = TRUE
  )

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J_inside <- nrow(d$alt_mapping) - 1  # outside option occupies one row
  n_params <- K_x + K_w + J_inside      # ASCs for all inside alts (outside = 0)

  set.seed(19)
  theta <- runif(n_params, -0.2, 0.2)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights,
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = TRUE
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_equal(H, t(H), tolerance = 1e-10)
  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-8))
})

test_that("mxl_bhhh_parallel scales linearly with uniform weights", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = c("w1", "w2"))

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)

  set.seed(23)
  theta <- runif(n_params, -0.2, 0.2)
  eta <- get_halton_normals(S = 30, N = d$N, K_w = K_w)

  H_w1 <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = rep(1, d$N),
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )
  H_w3 <- mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = rep(3, d$N),
    eta_draws = eta, rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE
  )

  # Weight enters exactly once per individual: w=3 scales the matrix by 3, not 9.
  expect_equal(H_w3, 3 * H_w1, tolerance = 1e-10)
})

test_that("mxl_bhhh_parallel errors on mismatched weights length", {
  dt <- create_small_mxl_data()
  d <- prepare_mxl_data(dt, "id", "alt", "choice",
                        covariate_cols = "x1",
                        random_var_cols = "w1")

  K_x <- ncol(d$X)
  K_w <- ncol(d$W)
  J <- nrow(d$alt_mapping)
  n_params <- K_x + K_w + (J - 1)

  theta <- rep(0, n_params)
  eta <- get_halton_normals(S = 10, N = d$N, K_w = K_w)

  expect_error(
    mxl_bhhh_parallel(
      theta = theta, X = d$X, W = d$W,
      alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      M = d$M, weights = rep(1, d$N + 1),  # wrong length
      eta_draws = eta, rc_dist = rep(0L, K_w),
      rc_correlation = FALSE, rc_mean = FALSE
    ),
    regexp = "weights length"
  )
})

test_that("run_mxlogit(se_method = 'bhhh') returns a valid choicer_mxl fit", {
  skip_on_cran()

  dt <- create_small_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 40L, se_method = "bhhh"
  )

  expect_s3_class(fit, "choicer_mxl")
  expect_identical(fit$se_method, "bhhh")
  expect_true(!is.null(fit$vcov))
  expect_true(!is.null(fit$se))
  expect_true(all(is.finite(fit$se)))
  expect_false(any(is.na(fit$se)))
})

test_that("ensure_vcov lazily recomputes BHHH for a bhhh-method fit", {
  skip_on_cran()

  dt <- create_identified_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    S = 40L, se_method = "bhhh"
  )

  original_vcov <- fit$vcov
  original_se <- fit$se

  # Clear vcov to force lazy recomputation via ensure_vcov()
  fit$vcov <- NULL
  fit$se <- NULL

  fit_restored <- choicer:::ensure_vcov(fit)

  expect_true(!is.null(fit_restored$vcov))
  expect_true(!is.null(fit_restored$se))
  # Must match the BHHH path, not the Hessian path: compare to the original.
  expect_equal(fit_restored$vcov, original_vcov, tolerance = 1e-10)
  expect_equal(fit_restored$se, original_se, tolerance = 1e-10)
})

test_that("BHHH and Hessian SEs have the same order of magnitude at convergence", {
  skip_on_cran()

  # Well-posed DGP: the generic `create_small_mxl_data` picks choices uniformly,
  # leaving the likelihood near-flat and producing pathological SEs for both
  # estimators. See `create_identified_mxl_data()` in setup.R.
  dt <- create_identified_mxl_data()

  fit_hess <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    S = 50L
  )

  # Re-fit starting from converged point and use BHHH SEs.
  fit_bhhh <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    S = 50L, se_method = "bhhh",
    theta_init = unname(coef(fit_hess))
  )

  hess_se <- fit_hess$se
  bhhh_se <- fit_bhhh$se

  expect_true(all(is.finite(hess_se)))
  expect_true(all(is.finite(bhhh_se)))
  expect_true(all(hess_se > 0))
  expect_true(all(bhhh_se > 0))

  # Same order of magnitude. The information matrix equality holds only
  # asymptotically at the true MLE; with simulation noise the tolerance
  # must be generous.
  expect_true(all(abs(log(bhhh_se / hess_se)) < 1.0))
})

test_that("summary prints BHHH label without error", {
  skip_on_cran()

  dt <- create_small_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 40L, se_method = "bhhh"
  )

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("BHHH", out)))
})
