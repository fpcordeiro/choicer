# Tests for elasticity computation functions

test_that("mnl_elasticities_parallel returns correct dimensions", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  # Elasticity with respect to first covariate (1-based index)
  elast <- mnl_elasticities_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    elast_var_idx = 1L,  # 1-based indexing
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  # Should return J x J matrix
  expect_equal(dim(elast), c(J, J))
})

test_that("mnl_elasticities_parallel produces finite values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  for (k in 1:K) {
    elast <- mnl_elasticities_parallel(
      theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
      inputs$M, inputs$weights, k, TRUE, FALSE
    )

    expect_true(all(is.finite(elast)),
                label = paste("covariate index", k))
  }
})

test_that("mnl_elasticities own-elasticities are on diagonal", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- c(1.0, -0.5, rep(0.1, J - 1))  # Positive beta for x1

  elast <- mnl_elasticities_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, 1L, TRUE, FALSE
  )

  # Diagonal contains own-elasticities
  own_elast <- diag(elast)
  expect_length(own_elast, J)
})

test_that("mnl_elasticities with zero beta gives zero elasticities", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)

  # Zero coefficient for x1
  theta <- c(0, 0.5, rep(0, J - 1))

  elast <- mnl_elasticities_parallel(
    theta, inputs$X, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, 1L, TRUE, FALSE
  )

  # All elasticities with respect to x1 should be zero
  expect_equal(elast, matrix(0, J, J), tolerance = 1e-12)
})

# MXL elasticities tests - skip if function not available
test_that("mxl_elasticities_parallel returns correct dimensions", {
  skip_if_not(exists("mxl_elasticities_parallel"),
              "mxl_elasticities_parallel not exported")

  dt <- create_small_mxl_data()
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
  theta <- c(0.3, log(0.5), log(0.4), 0.1, -0.1)

  # Elasticity with respect to fixed covariate x1
  elast <- mxl_elasticities_parallel(
    theta = theta,
    X = inputs$X,
    W = inputs$W,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    use_asc = TRUE,
    include_outside_option = FALSE,
    elast_var_idx = 1L,
    is_random_coef = FALSE
  )

  expect_equal(dim(elast), c(J, J))
})

test_that("mxl_elasticities_parallel works for random coefficient variable", {
  skip_if_not(exists("mxl_elasticities_parallel"),
              "mxl_elasticities_parallel not exported")

  dt <- create_small_mxl_data()
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
  theta <- c(0.3, log(0.5), log(0.4), 0.1, -0.1)

  # Elasticity with respect to random coefficient w1 (index 1 in W, 1-based)
  elast <- mxl_elasticities_parallel(
    theta = theta,
    X = inputs$X,
    W = inputs$W,
    alt_idx = inputs$alt_idx,
    choice_idx = inputs$choice_idx,
    M = inputs$M,
    weights = inputs$weights,
    eta_draws = eta_draws,
    rc_dist = rep(0L, K_w),
    rc_correlation = FALSE,
    rc_mean = FALSE,
    use_asc = TRUE,
    include_outside_option = FALSE,
    elast_var_idx = 1L,
    is_random_coef = TRUE
  )

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})


# =============================================================================
# MNL Diversion Ratio Tests
# =============================================================================

test_that("mnl_diversion_ratios_parallel returns correct dimensions", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_equal(dim(dr), c(J, J))
})

test_that("mnl_diversion_ratios_parallel produces finite non-negative values", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.5, 0.5)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_true(all(is.finite(dr)))
  expect_true(all(dr >= 0))
})

test_that("mnl_diversion_ratios_parallel has zero diagonal", {
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  expect_equal(diag(dr), rep(0, J))
})

test_that("mnl_diversion_ratios column sums equal 1 (no outside option, homogeneous choice sets)", {
  # With homogeneous choice sets and no outside option,
  # column sums of DR matrix should be 1 (all diverted demand stays inside)
  dt <- create_small_mnl_data()
  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  K <- ncol(inputs$X)
  J <- nrow(inputs$alt_mapping)
  theta <- runif(K + J - 1, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  col_sums <- colSums(dr)
  expect_equal(col_sums, rep(1, J), tolerance = 1e-10)
})

test_that("mnl_diversion_ratios IIA property: DR proportional to market share", {
  # Under MNL with IIA and identical choice sets, DR(j->k) = s_k / (1 - s_j)
  # This holds exactly when all individuals face the same covariates,
  # so we use zero X (only ASCs drive shares).
  set.seed(77)
  N <- 50
  J <- 3
  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1 = 0  # constant covariate -> ASCs fully determine shares
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  inputs <- prepare_mnl_data(dt, "id", "alt", "choice", "x1")

  K <- ncol(inputs$X)
  J_out <- nrow(inputs$alt_mapping)
  # theta: beta for x1, then J-1 ASCs with different values
  theta <- c(0.5, 0.3, -0.2)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  )

  shares <- as.vector(mnl_predict_shares(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = FALSE
  ))

  # For each column j, DR(j->k) = s_k / (1 - s_j)
  for (j in seq_len(J_out)) {
    expected_dr <- shares / (1 - shares[j])
    expected_dr[j] <- 0
    expect_equal(dr[, j], expected_dr, tolerance = 1e-10)
  }
})

# =============================================================================
# NL elasticities (S3 method on fitted choicer_nl)
# =============================================================================
# create_small_nl_data() (setup.R): N=30, J=6, two nests of size 3
#   nest 1 = alternatives {1, 2, 3}, nest 2 = alternatives {4, 5, 6}.

fit_nl_for_elast <- function(seed = 123) {
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

test_that("elasticities.choicer_nl returns labeled J x J matrix of finite values", {
  fit <- fit_nl_for_elast()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = "x1")

  expect_true(is.matrix(elast))
  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
  expect_true(!is.null(rownames(elast)))
  expect_true(!is.null(colnames(elast)))
  expect_equal(rownames(elast), colnames(elast))
})

test_that("elasticities.choicer_nl accepts a 1-based integer index into X", {
  fit <- fit_nl_for_elast()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = 1L)

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})

test_that("elasticities.choicer_nl own-elasticity sign opposes a positive beta", {
  fit <- fit_nl_for_elast()
  pm <- fit$param_map

  # Force a positive coefficient on x1 (index 1 of X) so the own-elasticity
  # sign is determined. choicer's elasticity convention (identical to the
  # shipped elasticities.choicer_mnl) reports the own-elasticity DIAGONAL as
  # NEGATIVE for a positive beta: it measures the proportional change in an
  # alternative's own choice probability, which for a desirable good is dampened
  # relative to the raw utility gain, yielding a negative diagonal entry.
  fit$coefficients[pm$beta[1]] <- 1.0

  elast <- elasticities(fit, elast_var = "x1")
  own_elast <- diag(elast)

  expect_length(own_elast, nrow(elast))
  expect_true(all(own_elast < 0))
})

test_that("elasticities.choicer_nl breaks IIA: within-nest vs cross-nest cross-elasticities differ", {
  # This is the defining NL property. For a price change in alternative j,
  # the proportional substitution toward same-nest alternatives differs from
  # substitution toward other-nest alternatives whenever lambda != 1.
  # Under MNL/IIA every off-diagonal entry in a column would be identical.
  fit <- fit_nl_for_elast()
  pm <- fit$param_map

  # Pin parameters to a clearly non-MNL configuration: a sizeable beta and
  # nest dissimilarity parameters well away from 1.
  fit$coefficients[pm$beta[1]] <- 1.0
  fit$coefficients[pm$lambda] <- 0.4

  elast <- elasticities(fit, elast_var = "x1")

  # Column 1 corresponds to a change in alternative 1, which lives in nest 1
  # ({1,2,3}). Same-nest off-diagonal rows: {2, 3}. Cross-nest rows: {4,5,6}.
  col <- elast[, 1]
  within_nest <- col[c(2, 3)]
  cross_nest <- col[c(4, 5, 6)]

  # The off-diagonal entries are NOT all equal (IIA would force equality).
  off_diag <- col[-1]
  expect_false(isTRUE(all.equal(
    rep(off_diag[1], length(off_diag)), off_diag,
    tolerance = 1e-6
  )))

  # Concretely, within-nest substitution differs from cross-nest substitution.
  expect_false(isTRUE(all.equal(
    mean(within_nest), mean(cross_nest),
    tolerance = 1e-6
  )))
})

test_that("elasticities.choicer_nl with near-zero beta gives near-zero elasticities", {
  fit <- fit_nl_for_elast()
  pm <- fit$param_map

  # Zero out the x1 coefficient: an attribute with no utility weight cannot
  # move any choice probability, so all x1-elasticities should be ~0.
  fit$coefficients[pm$beta[1]] <- 0

  elast <- elasticities(fit, elast_var = "x1")

  expect_equal(elast, matrix(0, nrow(elast), ncol(elast)),
               tolerance = 1e-8, ignore_attr = TRUE)
})

# =============================================================================
# NL diversion ratios (S3 method on fitted choicer_nl)
# =============================================================================
# diversion_ratios.choicer_nl mirrors the MNL signature: NO variable argument.

test_that("diversion_ratios.choicer_nl returns labeled J x J non-negative matrix", {
  fit <- fit_nl_for_elast()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_true(is.matrix(dr))
  expect_equal(dim(dr), c(J, J))
  expect_true(all(is.finite(dr)))
  expect_true(all(dr >= 0))
  expect_true(!is.null(rownames(dr)))
  expect_true(!is.null(colnames(dr)))
})

test_that("diversion_ratios.choicer_nl has zero diagonal", {
  fit <- fit_nl_for_elast()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_equal(unname(diag(dr)), rep(0, J))
})

test_that("diversion_ratios.choicer_nl column sums equal 1 (no outside option)", {
  # With homogeneous choice sets and no outside option, all demand diverted
  # away from j is recaptured by the remaining inside alternatives, so each
  # column of the diversion matrix sums to 1.
  fit <- fit_nl_for_elast()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  col_sums <- unname(colSums(dr))
  expect_equal(col_sums, rep(1, J), tolerance = 1e-8)
})

test_that("diversion_ratios.choicer_nl column sums equal 1 (with outside option)", {
  # create_nl_inputs() builds 6 alternatives in 3 nests where j = 0 is a
  # singleton with zero covariates; treat it as the outside option here.
  # We rebuild the dataset (rather than reuse the input object) so we can pass
  # outside_opt_label / include_outside_option through prepare_nl_data().
  set.seed(123)
  N <- 30
  J <- 6
  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J - 1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
  dt[, nest := data.table::fifelse(j == 0, "A",
                data.table::fifelse(j <= 2, "B", "C"))]

  fit <- run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "j",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    outside_opt_label = 0L,
    include_outside_option = TRUE,
    control = list(maxeval = 50L)
  )

  dr <- diversion_ratios(fit)

  J_total <- nrow(fit$alt_mapping)  # inside alts + outside option
  expect_equal(dim(dr), c(J_total, J_total))
  expect_equal(unname(diag(dr)), rep(0, J_total))

  col_sums <- unname(colSums(dr))
  expect_equal(col_sums, rep(1, J_total), tolerance = 1e-8)
})

# =============================================================================
# NL -> MNL equivalence sanity check
# =============================================================================

test_that("NL with all lambda = 1 reproduces MNL elasticities and diversion ratios", {
  # When every nest dissimilarity parameter lambda_g = 1, the nested logit
  # collapses exactly to the multinomial logit. An NL fit on the same data,
  # with its lambdas pinned to 1 and betas/ASCs aligned to the MNL fit, should
  # therefore yield the same elasticity and diversion matrices as MNL.
  #
  # create_small_nl_data() has 2 multi-member nests (so K_l = 2 lambdas). It is
  # not all-singleton, but lambda = 1 is the analytic MNL limit, so this is the
  # cleaner equivalence path than building a singleton structure.
  dt <- create_small_nl_data(seed = 123)

  fit_mnl <- run_mnlogit(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 50L)
  )

  fit_nl <- run_nestlogit(
    data = dt,
    id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), nest_col = "nest",
    control = list(maxeval = 50L)
  )

  # Align the NL coefficients to the MNL fit and set every lambda to 1.
  pm_nl <- fit_nl$param_map
  mnl_coef <- coef(fit_mnl)
  fit_nl$coefficients[pm_nl$beta] <- mnl_coef[c("x1", "x2")]
  fit_nl$coefficients[pm_nl$lambda] <- 1
  if (!is.null(pm_nl$asc)) {
    asc_names <- names(coef(fit_mnl))[fit_mnl$param_map$asc]
    fit_nl$coefficients[pm_nl$asc] <- mnl_coef[asc_names]
  }

  strip <- function(m) {
    out <- unclass(m)
    dimnames(out) <- NULL
    out
  }

  elast_mnl <- elasticities(fit_mnl, elast_var = "x1")
  elast_nl  <- elasticities(fit_nl, elast_var = "x1")
  expect_equal(strip(elast_nl), strip(elast_mnl), tolerance = 1e-5)

  dr_mnl <- diversion_ratios(fit_mnl)
  dr_nl  <- diversion_ratios(fit_nl)
  expect_equal(strip(dr_nl), strip(dr_mnl), tolerance = 1e-5)
})

test_that("mnl_diversion_ratios with outside option has correct dimensions and column sums", {
  set.seed(42)
  N <- 30
  J <- 4

  dt <- data.table::data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J-1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  inputs <- prepare_mnl_data(
    dt, "id", "j", "choice", c("x1", "x2"),
    outside_opt_label = 0L,
    include_outside_option = TRUE
  )

  K <- ncol(inputs$X)
  # alt_mapping includes the outside option row; ASCs are for inside alts only
  J_asc <- nrow(inputs$alt_mapping) - 1
  theta <- runif(K + J_asc, -0.3, 0.3)

  dr <- mnl_diversion_ratios_parallel(
    theta = theta,
    X = inputs$X,
    alt_idx = inputs$alt_idx,
    M = inputs$M,
    weights = inputs$weights,
    use_asc = TRUE,
    include_outside_option = TRUE
  )

  # Should be (J_asc + 1) x (J_asc + 1) matrix (inside alts + outside option)
  J_total <- J_asc + 1
  expect_equal(dim(dr), c(J_total, J_total))

  # Column sums should be 1
  col_sums <- colSums(dr)
  expect_equal(col_sums, rep(1, J_total), tolerance = 1e-10)

  # Diagonal should be 0
  expect_equal(diag(dr), rep(0, J_total))
})
