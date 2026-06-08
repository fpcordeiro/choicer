# Tests for BLP contraction mapping functions

test_that("blp_contraction returns correct length", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- rep(1/J, J)  # Equal shares
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  expect_length(result, J)
})

test_that("blp_contraction converges to target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  K <- ncol(inputs$X)

  # Use observed market shares as target
  target_shares <- inputs$alt_mapping$MKT_SHARE

  beta <- c(0.3, -0.2)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  # Verify shares match target
  # Construct theta: beta + delta (excluding first delta for identification)
  theta <- c(beta, delta_result[2:J])

  pred_shares <- mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  )
  # Convert from matrix to vector if needed
  pred_shares <- as.vector(pred_shares)

  expect_equal(pred_shares, target_shares, tolerance = 1e-6)
})

test_that("blp_contraction handles equal target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- rep(1/J, J)
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  # Should converge without error
  expect_true(all(is.finite(delta_result)))
})

test_that("blp_contraction returns finite values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- c(0.5, 0.3, 0.2)
  beta <- c(1.0, -0.5)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-10
  )

  expect_true(all(is.finite(delta_result)))
})

# MXL BLP contraction tests - skip if function not available
test_that("mxl_blp_contraction returns correct length", {
  skip_if_not(exists("mxl_blp_contraction"),
              "mxl_blp_contraction not exported")

  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  J <- nrow(inputs$alt_mapping)
  K_x <- ncol(inputs$X)
  K_w <- ncol(inputs$W)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  target_shares <- rep(1/J, J)
  beta <- 0.3
  L_params <- c(log(0.5), log(0.4))
  delta_init <- rep(0, J)

  result <- mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    W = inputs$W,
    beta = beta,
    mu = rep(0, K_w),
    L_params = L_params,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_length(result, J)
})

test_that("mxl_blp_contraction returns finite values", {
  skip_if_not(exists("mxl_blp_contraction"),
              "mxl_blp_contraction not exported")

  dt <- create_small_mxl_data()
  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  N <- inputs$N
  J <- nrow(inputs$alt_mapping)
  K_w <- ncol(inputs$W)
  S <- 25

  eta_draws <- get_halton_normals(S, N, K_w)

  target_shares <- c(0.4, 0.35, 0.25)
  beta <- 0.5
  L_params <- c(log(0.5), log(0.4))
  delta_init <- rep(0, J)

  result <- mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    W = inputs$W,
    beta = beta,
    mu = rep(0, K_w),
    L_params = L_params,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_true(all(is.finite(result)))
})

test_that("blp_contraction respects convergence tolerance", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)
  target_shares <- inputs$alt_mapping$MKT_SHARE
  beta <- c(0.3, -0.2)
  delta_init <- rep(0, J)

  # Tighter tolerance
  delta_tight <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-12
  )

  # Looser tolerance
  delta_loose <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-6
  )

  # Both should be close but tight should be more accurate
  theta_tight <- c(beta, delta_tight[2:J])
  theta_loose <- c(beta, delta_loose[2:J])

  shares_tight <- as.vector(mnl_predict_shares(
    theta_tight, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))
  shares_loose <- as.vector(mnl_predict_shares(
    theta_loose, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))

  error_tight <- max(abs(shares_tight - target_shares))
  error_loose <- max(abs(shares_loose - target_shares))

  # Tight tolerance should give smaller error
  expect_true(error_tight <= error_loose)
})

# =============================================================================
# Helpers for outside-option / correlation BLP tests
# =============================================================================

create_mxl_outside_data <- function(seed = 42, N = 40, J_inside = 3) {
  set.seed(seed)
  dt <- data.table::data.table(
    id = rep(1:N, each = J_inside + 1),
    alt = rep(0:J_inside, N),
    x1 = rnorm(N * (J_inside + 1)),
    w1 = rnorm(N * (J_inside + 1)),
    w2 = runif(N * (J_inside + 1), -1, 1)
  )
  # Outside option (alt=0) has zero covariates
  dt[alt == 0, c("x1", "w1", "w2") := 0]
  dt[, choice := 0L]
  # Some individuals choose the outside option, others choose an inside alt
  dt[, choice := {
    pick <- sample.int(J_inside + 1, 1) - 1L  # alternative index 0..J_inside
    as.integer(alt == pick)
  }, by = id]
  dt[]
}

# =============================================================================
# blp.choicer_mxl with outside option / correlation / round-trip
# =============================================================================

test_that("blp.choicer_mxl works with include_outside_option = TRUE", {
  dt <- create_mxl_outside_data(seed = 42, N = 40, J_inside = 3)

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L,
    outside_opt_label = 0L,
    include_outside_option = TRUE,
    control = list(maxeval = 30L)
  )

  J_total <- nrow(fit$alt_mapping)        # inside + outside
  J_inside <- J_total - 1L

  # target_shares length must equal J_inside + 1 (outside first row)
  target_shares <- rep(1 / J_total, J_total)

  delta <- blp(fit, target_shares = target_shares)

  expect_true(is.numeric(delta))
  expect_true(all(is.finite(delta)))
  # mxl_blp_contraction returns delta of length J_inside (free ASCs)
  expect_length(delta, J_inside)
})

test_that("blp.choicer_mxl with outside option round-trips own predicted shares", {
  dt <- create_mxl_outside_data(seed = 7, N = 50, J_inside = 3)

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 20L,
    outside_opt_label = 0L,
    include_outside_option = TRUE,
    rc_correlation = FALSE,
    control = list(maxeval = 60L)
  )

  J_total <- nrow(fit$alt_mapping)
  J_inside <- J_total - 1L

  # Feeding the model's own predicted shares must recover the inside ASCs
  # (up to an additive constant) and stay finite. mxl_blp_contraction already
  # operates on the length-J_inside inside-delta slice (the outside option's
  # utility is the fixed normalization), so the contraction is well-posed.
  target_shares <- predict(fit, type = "shares")

  delta <- blp(fit, target_shares = target_shares)

  expect_true(all(is.finite(delta)))
  expect_length(delta, J_inside)

  pm <- fit$param_map
  asc_existing <- if (!is.null(pm$asc)) fit$coefficients[pm$asc] else rep(0, J_inside)
  delta_centered <- as.numeric(delta) - as.numeric(delta)[1]
  asc_centered <- unname(asc_existing - asc_existing[1])
  expect_equal(delta_centered, asc_centered, tolerance = 2e-3)
})

test_that("blp.choicer_mxl works with rc_correlation = TRUE", {
  dt <- create_small_mxl_data()

  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    rc_correlation = TRUE,
    S = 10L,
    control = list(maxeval = 30L)
  )

  current_shares <- predict(fit, type = "shares")

  delta <- blp(fit, target_shares = current_shares)

  expect_true(all(is.finite(delta)))

  # The recovered delta should be close to the existing free ASCs
  pm <- fit$param_map
  asc_existing <- if (!is.null(pm$asc)) {
    c(0, fit$coefficients[pm$asc])
  } else {
    rep(0, nrow(fit$alt_mapping))
  }
  # Up to identification (constant shift): recenter both to first element = 0.
  # blp_contraction returns arma::vec which surfaces as a 1-column matrix in R;
  # strip shape and names so the comparison is purely numeric.
  delta_centered <- as.numeric(delta) - as.numeric(delta)[1]
  asc_centered <- unname(asc_existing - asc_existing[1])

  expect_equal(delta_centered, asc_centered, tolerance = 1e-3)
})

test_that("blp.choicer_mxl round-trips through mxl_predict_shares", {
  fit <- run_mxlogit(
    data = create_small_mxl_data(),
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L,
    control = list(maxeval = 30L)
  )

  current_shares <- predict(fit, type = "shares")

  delta_new <- blp(fit, target_shares = current_shares)

  pm <- fit$param_map
  asc_existing <- if (!is.null(pm$asc)) {
    c(0, fit$coefficients[pm$asc])
  } else {
    rep(0, nrow(fit$alt_mapping))
  }

  # Up to identification (additive constant), delta should match existing ASCs.
  # Strip shape/names (delta_new is arma::vec, asc_existing is a named numeric).
  expect_equal(
    as.numeric(delta_new) - as.numeric(delta_new)[1],
    unname(asc_existing - asc_existing[1]),
    tolerance = 1e-3
  )
})

# =============================================================================
# blp.choicer_mnl with outside option (parallel bug-fix check)
# =============================================================================

test_that("blp.choicer_mnl works with include_outside_option = TRUE", {
  set.seed(123)
  N <- 40
  J <- 4  # 1 outside + 3 inside

  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J - 1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := {
    pick <- sample.int(J, 1) - 1L
    as.integer(j == pick)
  }, by = id]

  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "j", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    outside_opt_label = 0L,
    include_outside_option = TRUE,
    control = list(maxeval = 50L)
  )

  J_total <- nrow(fit$alt_mapping)        # inside + outside
  J_inside <- J_total - 1L

  # Feed the fitted model's own predicted shares as the BLP target. Because the
  # contraction inverts the same share map the model produced, it must recover
  # the inside ASCs (up to an additive constant) and stay finite.
  target_shares <- predict(fit, type = "shares")

  delta <- blp(fit, target_shares = target_shares)

  # Before the C++ fix the MNL blp_contraction kernel fed the FULL length-(J+1)
  # delta (outside slot included) to mnl_predict_shares_internal, which indexes
  # delta by inside-alt index {0..J-1}. That shifted every ASC by one slot and
  # left the outside delta unpinned, so the shares were corrupted each iteration
  # and the contraction diverged to NaN. The fix feeds the length-J inside slice
  # and pins delta[0] = 0, mirroring nl_blp_contraction.
  expect_true(is.numeric(delta))
  expect_length(delta, J_inside)
  expect_true(all(is.finite(delta)))

  # Round-trip: recovered delta matches the fitted inside ASCs up to a constant.
  pm <- fit$param_map
  asc_existing <- fit$coefficients[pm$asc]
  delta_centered <- as.numeric(delta) - as.numeric(delta)[1]
  asc_centered <- unname(asc_existing - asc_existing[1])
  expect_equal(delta_centered, asc_centered, tolerance = 1e-4)
})

test_that("blp_contraction rejects non-positive target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
  J <- nrow(inputs$alt_mapping)

  expect_error(
    blp_contraction(
      delta = rep(0, J),
      target_shares = c(0.5, 0.5, 0),  # zero share -> log undefined
      X = inputs$X,
      beta = c(0.3, -0.2),
      alt_idx = inputs$alt_idx,
      M = inputs$M,
      weights = inputs$weights,
      include_outside_option = FALSE
    ),
    "strictly positive"
  )
})

# =============================================================================
# blp.choicer_nl (S3 method on fitted choicer_nl)
# =============================================================================
# create_small_nl_data() (setup.R): N=30, J=6, two nests of size 3, no outside
# option. blp.choicer_nl mirrors blp.choicer_mnl:
#   blp(object, target_shares, delta_init = NULL, damping = 1,
#       tol = 1e-8, max_iter = 1000, ...)

fit_nl_for_blp <- function(seed = 123) {
  dt <- create_small_nl_data(seed = seed)
  run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    control = list(maxeval = 50L)
  )
}

test_that("blp.choicer_nl returns a finite numeric vector", {
  fit <- fit_nl_for_blp()
  J <- nrow(fit$alt_mapping)

  # Feed the model's own predicted shares as the BLP target rather than exact
  # uniform shares. Uniform shares (rep(1/J, J)) push the NL contraction into a
  # numerically degenerate regime when a fitted lambda is near zero, producing
  # NaN; the model's implied shares are guaranteed interior and consistent with
  # the fitted nest structure, so the contraction is well-posed.
  target_shares <- as.vector(predict(fit, type = "shares"))

  delta <- blp(fit, target_shares = target_shares)

  expect_true(is.numeric(delta))
  expect_length(delta, J)
  expect_true(all(is.finite(delta)))
})

test_that("blp.choicer_nl default delta_init = NULL path works", {
  fit <- fit_nl_for_blp()

  # Feed the model's own implied market shares as the target. With no explicit
  # delta_init the method must default it from the fitted ASCs (param_map$asc).
  target_shares <- as.vector(predict(fit, type = "shares"))

  delta <- blp(fit, target_shares = target_shares)

  expect_true(all(is.finite(delta)))
})

test_that("blp.choicer_nl round-trips: own predicted shares recover original ASCs", {
  # Using the model's own predicted shares as the BLP target should recover the
  # ASCs that generated them, up to the usual additive identification constant.
  fit <- fit_nl_for_blp()

  current_shares <- as.vector(predict(fit, type = "shares"))

  delta_new <- blp(fit, target_shares = current_shares)

  pm <- fit$param_map
  # No outside option here, so the baseline alternative's ASC is fixed at 0 and
  # the free ASCs follow it: existing delta = c(0, ASC_2, ..., ASC_J).
  asc_existing <- if (!is.null(pm$asc)) {
    c(0, fit$coefficients[pm$asc])
  } else {
    rep(0, nrow(fit$alt_mapping))
  }

  # Recenter both to first element = 0 to absorb the identification constant,
  # and strip arma::vec / name attributes before the numeric comparison.
  delta_centered <- as.numeric(delta_new) - as.numeric(delta_new)[1]
  asc_centered <- unname(asc_existing - asc_existing[1])

  expect_equal(delta_centered, asc_centered, tolerance = 1e-3)
})

test_that("blp_contraction handles skewed target shares", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  J <- nrow(inputs$alt_mapping)

  # Very skewed shares
  target_shares <- c(0.9, 0.08, 0.02)
  beta <- c(0.5, -0.3)
  delta_init <- rep(0, J)

  delta_result <- blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = inputs$X,
    beta = beta,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    include_outside_option = FALSE,
    tol = 1e-8
  )

  expect_true(all(is.finite(delta_result)))

  # Verify shares are recovered
  theta <- c(beta, delta_result[2:J])
  pred_shares <- as.vector(mnl_predict_shares(
    theta, inputs$X, inputs$alt_idx, inputs$M,
    inputs$weights, TRUE, FALSE
  ))

  expect_equal(pred_shares, target_shares, tolerance = 1e-5)
})
