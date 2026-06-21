# Tests for nl_loglik_hessian_parallel (analytical Hessian of negated NL log-likelihood)
# Builder unit tests — driven by spec.md, not test-spec.md.

# ---------------------------------------------------------------------------
# Helper: build a small NL dataset with 4 alternatives in 2 nests
# (delegates to the shared create_nl_inputs from setup.R)
# ---------------------------------------------------------------------------

# ============================================================
# Block 1: structure checks
# ============================================================

test_that("nl_loglik_hessian_parallel returns a finite, square, symmetric matrix", {
  inputs <- create_nl_inputs()

  J    <- nrow(inputs$alt_mapping)
  K_x  <- ncol(inputs$X)
  K_l  <- sum(table(inputs$nest_idx) > 1)   # non-singleton nests
  P    <- K_x + K_l + (J - 1)               # total params (with ASCs)

  theta <- c(rep(0.1, K_x), rep(0.6, K_l), rep(0, J - 1))

  H <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE,
    include_outside_option = inputs$include_outside_option
  )

  expect_true(is.matrix(H))
  expect_equal(nrow(H), P)
  expect_equal(ncol(H), P)
  expect_true(all(is.finite(H)))
  # Symmetry (enforced by the function itself)
  expect_equal(H, t(H), tolerance = 1e-12)
})

test_that("nl_loglik_hessian_parallel without ASCs has dimensions K_x + K_l", {
  inputs <- create_nl_inputs()

  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  theta <- c(rep(0.1, K_x), rep(0.6, K_l))

  H <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = FALSE,
    include_outside_option = inputs$include_outside_option
  )

  expect_equal(nrow(H), K_x + K_l)
  expect_equal(ncol(H), K_x + K_l)
  expect_true(all(is.finite(H)))
})

# ============================================================
# Block 2: numerical accuracy — match numeric oracle
# ============================================================

test_that("nl_loglik_hessian_parallel matches numeric oracle at random theta", {
  skip_if_not_installed("numDeriv")

  inputs <- create_nl_inputs(seed = 456)

  J   <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  set.seed(42)
  theta <- c(
    runif(K_x, -0.3, 0.3),
    runif(K_l, 0.3, 0.8),
    runif(J - 1, -0.2, 0.2)
  )

  H_an <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  H_num <- nl_loglik_numeric_hessian(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(H_an, H_num, tolerance = TOL_HESS)
})

test_that("nl_loglik_hessian_parallel matches numeric oracle at zero beta", {
  inputs <- create_nl_inputs(seed = 789)

  J   <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  theta <- c(rep(0, K_x), 0.6, 0.4, rep(0, J - 1))

  H_an <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  H_num <- nl_loglik_numeric_hessian(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(H_an, H_num, tolerance = TOL_HESS)
})

test_that("nl_loglik_hessian_parallel matches numeric oracle without ASCs", {
  inputs <- create_nl_inputs(seed = 321)

  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  set.seed(7)
  theta <- c(runif(K_x, -0.4, 0.4), runif(K_l, 0.2, 0.9))

  H_an <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = FALSE, include_outside_option = inputs$include_outside_option
  )

  H_num <- nl_loglik_numeric_hessian(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = FALSE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(H_an, H_num, tolerance = TOL_HESS)
})

test_that("nl_loglik_hessian_parallel is accurate near lambda = 1 (MNL limit)", {
  inputs <- create_nl_inputs(seed = 123)

  J   <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  # Lambda very close to 1 -> approaches MNL
  theta <- c(0.3, -0.2, 0.99, 0.99, rep(0, J - 1))

  H_an <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  H_num <- nl_loglik_numeric_hessian(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(H_an, H_num, tolerance = TOL_HESS)
})

test_that("nl_loglik_hessian_parallel is accurate at small lambda (0.2, 0.3)", {
  inputs <- create_nl_inputs(seed = 654)

  J   <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)

  theta <- c(0.1, -0.1, 0.2, 0.3, rep(0, J - 1))

  H_an <- nl_loglik_hessian_parallel(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  H_num <- nl_loglik_numeric_hessian(
    theta,
    inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$nest_idx, inputs$M, inputs$weights,
    use_asc = TRUE, include_outside_option = inputs$include_outside_option
  )

  expect_equal(H_an, H_num, tolerance = TOL_HESS)
})

# ============================================================
# Block 3: positive semi-definiteness at the MLE
# ============================================================

test_that("nl_loglik_hessian_parallel is positive definite at MLE", {
  skip_on_cran()  # optimizer-dependent

  set.seed(2024)
  N <- 400; J <- 4
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    x2  = rnorm(N * J),
    nest = ifelse(rep(1:J, N) <= 2, "A", "B")
  )
  dt[, util := 0.5 * x1 - 0.3 * x2 + c(0, 0.3, -0.2, -0.4)[alt], by = seq_len(nrow(dt))]
  dt[, choice := {
    r  <- runif(1)
    pr <- exp(util) / sum(exp(util))
    cs <- cumsum(pr)
    as.integer(alt == alt[which(cs >= r)[1]])
  }, by = id]

  fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                       se_method = "hessian")

  d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
  theta <- coef(fit)
  H <- nl_loglik_hessian_parallel(
    theta,
    d$X, d$alt_idx, d$choice_idx, d$nest_idx, d$M, d$weights
  )

  eigvals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigvals > -1e-6),
    label = paste("Min eigenvalue:", min(eigvals)))
})

# ============================================================
# Block 4: determinism
# ============================================================

test_that("nl_loglik_hessian_parallel is deterministic across repeated calls", {
  inputs <- create_nl_inputs()

  J   <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_l <- sum(table(inputs$nest_idx) > 1)
  theta <- c(0.3, -0.2, 0.6, 0.4, rep(0.1, J - 1))

  results <- replicate(5, {
    nl_loglik_hessian_parallel(
      theta,
      inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$nest_idx, inputs$M, inputs$weights,
      use_asc = TRUE, include_outside_option = inputs$include_outside_option
    )
  }, simplify = FALSE)

  for (i in 2:5) {
    expect_equal(results[[i]], results[[1]], tolerance = 1e-12)
  }
})

# ============================================================
# Block 5: se_method dispatch in run_nestlogit
# ============================================================

test_that("run_nestlogit stores se_method field in fitted object", {
  skip_on_cran()

  set.seed(77)
  N <- 150; J <- 4
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    nest = ifelse(rep(1:J, N) <= 2, "A", "B")
  )
  dt[, choice := {
    pr <- runif(J); pr <- pr / sum(pr)
    as.integer(alt == sample(1:J, 1, prob = pr))
  }, by = id]

  fit_h <- run_nestlogit(dt, "id", "alt", "choice", "x1", "nest",
                         se_method = "hessian")
  fit_n <- run_nestlogit(dt, "id", "alt", "choice", "x1", "nest",
                         se_method = "numeric")

  expect_equal(fit_h$se_method, "hessian")
  expect_equal(fit_n$se_method, "numeric")
})

test_that("run_nestlogit default se_method is 'hessian'", {
  skip_on_cran()

  set.seed(88)
  N <- 100; J <- 4
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    nest = ifelse(rep(1:J, N) <= 2, "A", "B")
  )
  dt[, choice := {
    pr <- runif(J); pr <- pr / sum(pr)
    as.integer(alt == sample(1:J, 1, prob = pr))
  }, by = id]

  fit <- run_nestlogit(dt, "id", "alt", "choice", "x1", "nest")
  expect_equal(fit$se_method, "hessian")
})

test_that("run_nestlogit rejects invalid se_method", {
  set.seed(99)
  N <- 100; J <- 4
  dt <- data.table::data.table(
    id  = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1  = rnorm(N * J),
    nest = ifelse(rep(1:J, N) <= 2, "A", "B")
  )
  dt[, choice := {
    pr <- runif(J); pr <- pr / sum(pr)
    as.integer(alt == sample(1:J, 1, prob = pr))
  }, by = id]

  expect_error(
    run_nestlogit(dt, "id", "alt", "choice", "x1", "nest",
                  se_method = "not-a-method"),
    "'arg' should be one of"
  )
})
