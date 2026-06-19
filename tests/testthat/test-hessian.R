# Tests for Hessian computations

# --- MNL Hessian tests ---

test_that("mnl_loglik_hessian_parallel matches numerical Hessian", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  set.seed(42)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Analytical Hessian
  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  # Numerical Hessian
  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

test_that("mnl_loglik_hessian_parallel is symmetric", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  H <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  expect_equal(H, t(H), tolerance = 1e-12)
})

test_that("mnl_loglik_hessian at zero parameters", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- rep(0, K + J - 1)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, FALSE
  )

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})


# --- MNL Hessian regression tests for block-decomposition branches ---

# S3: use_asc = FALSE — only the BB block is populated; J_asc == 0,
# so the delta/BD/DD sub-matrices are entirely absent.
test_that("mnl Hessian block-decomp S3: use_asc=FALSE matches numerical oracle", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  # use_asc = FALSE: theta is just the K beta parameters (no ASC block)
  set.seed(42)
  theta <- runif(K, -0.3, 0.3)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, FALSE, FALSE
  )

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, FALSE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_true(all(is.finite(H_anal)),
              label = "Hessian (use_asc=FALSE) is finite")
  expect_equal(H_anal, t(H_anal), tolerance = 1e-10,
               label = "Hessian (use_asc=FALSE) is symmetric")
  expect_equal(H_anal, H_num, tolerance = TOL_HESS,
               label = "Hessian (use_asc=FALSE) matches numerical oracle")
})

test_that("mnl Hessian block-decomp S3: use_asc=FALSE yields finite non-negative SEs", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"), use_asc = FALSE)
  ses <- fit$se
  expect_true(all(is.finite(ses)),
              label = "SEs are finite when use_asc=FALSE")
  expect_true(all(ses >= 0),
              label = "SEs are non-negative when use_asc=FALSE")
})

# S4: include_outside_option = TRUE — the outside option contributes a zero
# utility row; the ASC indexing shifts so every inside alternative has a free
# delta (J_asc == J, not J-1), and choice_idx == 0 signals outside-option
# choosers.  This exercises the `delta_pos = alt_idx0_i[a]` branch in the
# block-decomposition Hessian (vs. `a_id - 1` in the baseline).
test_that("mnl Hessian block-decomp S4: include_outside_option=TRUE matches numerical oracle", {
  skip_if_not_installed("numDeriv")

  set.seed(77)
  N <- 30; J <- 3
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    x2  = runif(N * J, -1, 1)
  )
  # Mark ~3/4 of individuals as inside-option choosers; rest are outside-option choosers
  dt[, choice := 0L]
  set.seed(77)
  chosen_ids <- sample(1:N, ceiling(3 * N / 4))
  for (id_i in chosen_ids) {
    rows_i <- which(dt$id == id_i)
    dt[rows_i[sample(J, 1)], choice := 1L]
  }

  inputs <- prepare_mnl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    include_outside_option = TRUE
  )

  K <- ncol(inputs$X)
  J_inside <- nrow(inputs$alt_mapping) - 1L  # exclude outside-option row
  # With include_outside_option=TRUE all J inside alts have free ASCs
  n_params <- K + J_inside
  set.seed(42)
  theta <- runif(n_params, -0.3, 0.3)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, TRUE, TRUE
  )

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, TRUE, TRUE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_true(all(is.finite(H_anal)),
              label = "Hessian (include_outside_option=TRUE) is finite")
  expect_equal(H_anal, t(H_anal), tolerance = 1e-10,
               label = "Hessian (include_outside_option=TRUE) is symmetric")
  expect_equal(H_anal, H_num, tolerance = TOL_HESS,
               label = "Hessian (include_outside_option=TRUE) matches numerical oracle")
})

test_that("mnl Hessian block-decomp S4: include_outside_option=TRUE yields finite non-negative SEs", {
  set.seed(77)
  N <- 30; J <- 3
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    x2  = runif(N * J, -1, 1)
  )
  dt[, choice := 0L]
  set.seed(77)
  chosen_ids <- sample(1:N, ceiling(3 * N / 4))
  for (id_i in chosen_ids) {
    rows_i <- which(dt$id == id_i)
    dt[rows_i[sample(J, 1)], choice := 1L]
  }
  fit <- run_mnlogit(
    dt, "id", "alt", "choice", c("x1", "x2"),
    include_outside_option = TRUE
  )
  ses <- fit$se
  expect_true(all(is.finite(ses)),
              label = "SEs are finite when include_outside_option=TRUE")
  expect_true(all(ses >= 0),
              label = "SEs are non-negative when include_outside_option=TRUE")
})

# S6: non-uniform per-individual weights (w_i != 1) — exercises the weight
# accumulation path in the block-decomposition Hessian.  The analytical
# Hessian must still match the jacobian-of-gradient numerical oracle, because
# the numerical oracle also uses the same weighted objective.
test_that("mnl Hessian block-decomp S6: non-uniform weights matches numerical oracle", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  set.seed(42)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Replace unit weights with random positive weights
  N <- inputs$N
  set.seed(7)
  w_rand <- runif(N, 0.5, 2.5)

  H_anal <- mnl_loglik_hessian_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, w_rand, TRUE, FALSE
  )

  obj_fn <- function(th) {
    mnl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, w_rand, TRUE, FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_true(all(is.finite(H_anal)),
              label = "Hessian (non-uniform weights) is finite")
  expect_equal(H_anal, t(H_anal), tolerance = 1e-10,
               label = "Hessian (non-uniform weights) is symmetric")
  expect_equal(H_anal, H_num, tolerance = TOL_HESS,
               label = "Hessian (non-uniform weights) matches numerical oracle")
})

test_that("mnl Hessian block-decomp S6: non-uniform weights yields finite non-negative SEs", {
  set.seed(42)
  N <- 20; J <- 3
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    x2  = runif(N * J, -1, 1)
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
  set.seed(7)
  w_rand <- runif(N, 0.5, 2.5)
  fit <- run_mnlogit(
    dt, "id", "alt", "choice", c("x1", "x2"),
    weights = w_rand
  )
  ses <- fit$se
  expect_true(all(is.finite(ses)),
              label = "SEs are finite with non-uniform weights")
  expect_true(all(ses >= 0),
              label = "SEs are non-negative with non-uniform weights")
})

# --- MXL Hessian tests ---

test_that("mxl_hessian_parallel matches numerical Hessian", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 111)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.2, 0.2),
    log(runif(K_w, 0.4, 0.7)),
    runif(J - 1, -0.1, 0.1)
  )

  H_anal <- mxl_hessian_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rep(0L, K_w),
      rc_correlation = FALSE, rc_mean = FALSE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

test_that("mxl_hessian_parallel is symmetric", {
  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 20

  eta_draws <- get_halton_normals(S, N, K_w)
  theta <- c(0.1, log(0.5), log(0.4), 0.1, -0.1)

  H <- mxl_hessian_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("mxl_hessian_parallel matches numerical Hessian with mixed rc_dist", {
  skip_if_not_installed("numDeriv")

  dt <- create_small_mxl_data(seed = 333)
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  J <- nrow(inputs$alt_mapping)
  S <- 30

  eta_draws <- get_halton_normals(S, N, K_w)
  rc_dist <- c(0L, 1L)  # mix of normal + log-normal

  set.seed(42)
  theta <- c(
    runif(K_x, -0.2, 0.2),
    runif(K_w, -0.1, 0.1),           # mu (rc_mean = TRUE)
    log(runif(K_w, 0.4, 0.7)),       # L diagonal on log scale
    runif(J - 1, -0.1, 0.1)
  )

  H_anal <- mxl_hessian_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_draws,
    rc_dist = rc_dist,
    rc_correlation = FALSE, rc_mean = TRUE,
    use_asc = TRUE, include_outside_option = FALSE
  )

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      th, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, eta_draws,
      rc_dist = rc_dist,
      rc_correlation = FALSE, rc_mean = TRUE,
      use_asc = TRUE, include_outside_option = FALSE
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_anal, H_num, tolerance = TOL_HESS)
})

# --- NL Hessian tests ---

test_that("nl_loglik_numeric_hessian returns correct dimensions", {
  inputs <- create_nl_inputs()

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  n_params <- K_x + K_l + J - 1
  theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))

  H <- nl_loglik_numeric_hessian(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(dim(H), c(n_params, n_params))
  expect_true(all(is.finite(H)))
})

test_that("nl_loglik_numeric_hessian matches numDeriv", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 222)

  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  set.seed(42)
  theta <- c(runif(K_x, -0.2, 0.2), runif(K_l, 0.3, 0.7), runif(J - 1, -0.1, 0.1))

  H_pkg <- nl_loglik_numeric_hessian(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  obj_fn <- function(th) {
    nl_loglik_gradient_parallel(
      th, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )$objective
  }

  H_num <- numDeriv::hessian(obj_fn, theta, method = "Richardson")

  expect_equal(H_pkg, H_num, tolerance = TOL_HESS)
})

# Note: create_nl_inputs is defined in setup.R
