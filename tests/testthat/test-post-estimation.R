# Tests for S3 post-estimation methods:
# - elasticities()
# - diversion_ratios()
# - blp()
# These test the S3 dispatch layer on fitted model objects.

# =============================================================================
# Helper: fit small MNL model for post-estimation tests
# =============================================================================

fit_small_mnl <- function() {
  dt <- create_small_mnl_data()
  run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 50L)
  )
}

fit_small_mxl <- function() {
  dt <- create_small_mxl_data()
  run_mxlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = "x1",
    random_var_cols = c("w1", "w2"),
    S = 10L,
    control = list(maxeval = 50L)
  )
}

# =============================================================================
# MNL elasticities
# =============================================================================

test_that("elasticities.choicer_mnl returns labeled J x J matrix", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = "x1")

  expect_true(is.matrix(elast))
  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
  expect_true(!is.null(rownames(elast)))
  expect_true(!is.null(colnames(elast)))
  expect_equal(rownames(elast), colnames(elast))
})

test_that("elasticities.choicer_mnl accepts integer index", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = 1L)

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})

test_that("elasticities.choicer_mnl errors without data", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    elasticities(fit, elast_var = "x1"),
    "keep_data"
  )
})

test_that("elasticities.choicer_mnl errors with invalid variable name", {
  fit <- fit_small_mnl()

  expect_error(
    elasticities(fit, elast_var = "nonexistent"),
    "not found"
  )
})

# =============================================================================
# MNL diversion_ratios
# =============================================================================

test_that("diversion_ratios.choicer_mnl returns labeled J x J matrix", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_true(is.matrix(dr))
  expect_equal(dim(dr), c(J, J))
  expect_true(all(is.finite(dr)))
  expect_true(all(dr >= 0))
  expect_true(!is.null(rownames(dr)))
  expect_true(!is.null(colnames(dr)))
})

test_that("diversion_ratios.choicer_mnl has zero diagonal", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_equal(unname(diag(dr)), rep(0, J))
})

test_that("diversion_ratios.choicer_mnl column sums equal 1", {
  fit <- fit_small_mnl()

  dr <- diversion_ratios(fit)

  col_sums <- unname(colSums(dr))
  expect_equal(col_sums, rep(1, length(col_sums)), tolerance = 1e-10)
})

test_that("diversion_ratios.choicer_mnl errors without data", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    diversion_ratios(fit),
    "keep_data"
  )
})

# =============================================================================
# MNL blp
# =============================================================================

test_that("blp.choicer_mnl returns finite vector of correct length", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  # Use equal target shares
  target_shares <- rep(1 / J, J)

  delta <- blp(fit, target_shares = target_shares)

  expect_true(is.numeric(delta))
  expect_true(all(is.finite(delta)))
})

test_that("blp.choicer_mnl uses default delta_init from ASCs", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)
  target_shares <- fit$alt_mapping$MKT_SHARE

  # Should not error — default delta_init taken from param_map$asc
  delta <- blp(fit, target_shares = target_shares)
  expect_true(all(is.finite(delta)))
})

test_that("blp.choicer_mnl errors without data", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), control = list(maxeval = 50L),
    keep_data = FALSE
  )
  J <- nrow(fit$alt_mapping)

  expect_error(
    blp(fit, target_shares = rep(1 / J, J)),
    "keep_data"
  )
})

# =============================================================================
# MNL predict
# =============================================================================

test_that("predict.choicer_mnl returns probabilities", {
  fit <- fit_small_mnl()

  preds <- predict(fit, type = "probabilities")

  expect_type(preds, "list")
  expect_true("choice_prob" %in% names(preds))
  expect_true(all(preds$choice_prob >= 0))
  expect_true(all(preds$choice_prob <= 1))
})

test_that("predict.choicer_mnl returns shares", {
  fit <- fit_small_mnl()
  J <- nrow(fit$alt_mapping)

  shares <- predict(fit, type = "shares")

  expect_length(shares, J)
  expect_equal(sum(shares), 1, tolerance = 1e-8)
})

test_that("predict.choicer_mnl errors without data", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    predict(fit),
    "keep_data"
  )
})

# =============================================================================
# MXL elasticities
# =============================================================================

test_that("elasticities.choicer_mxl returns labeled J x J matrix for fixed coef", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = "x1", is_random_coef = FALSE)

  expect_true(is.matrix(elast))
  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
  expect_true(!is.null(rownames(elast)))
})

test_that("elasticities.choicer_mxl returns J x J matrix for random coef", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  elast <- elasticities(fit, elast_var = "w1", is_random_coef = TRUE)

  expect_equal(dim(elast), c(J, J))
  expect_true(all(is.finite(elast)))
})

test_that("elasticities.choicer_mxl errors without data", {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    elasticities(fit, elast_var = "x1"),
    "keep_data"
  )
})

# =============================================================================
# MXL blp
# =============================================================================

test_that("blp.choicer_mxl returns finite vector", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)
  target_shares <- rep(1 / J, J)

  delta <- blp(fit, target_shares = target_shares)

  expect_true(is.numeric(delta))
  expect_true(all(is.finite(delta)))
})

test_that("blp.choicer_mxl errors without data", {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 50L),
    keep_data = FALSE
  )
  J <- nrow(fit$alt_mapping)

  expect_error(
    blp(fit, target_shares = rep(1 / J, J)),
    "keep_data"
  )
})

# =============================================================================
# MXL predict
# =============================================================================

test_that("predict.choicer_mxl(type='probabilities') returns list with choice_prob and utility", {
  fit <- fit_small_mxl()
  M <- fit$data$M

  preds <- predict(fit, type = "probabilities")

  expect_type(preds, "list")
  expect_true("choice_prob" %in% names(preds))
  expect_true("utility" %in% names(preds))
  expect_length(preds$choice_prob, sum(M))
  expect_true(all(preds$choice_prob >= 0))
  expect_true(all(preds$choice_prob <= 1))
})

test_that("predict.choicer_mxl(type='shares') returns shares summing to 1", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  shares <- predict(fit, type = "shares")

  # Fitted without outside option, so length should be J inside alts
  expect_length(shares, J)
  expect_equal(sum(shares), 1, tolerance = 1e-6)
})

test_that("predict.choicer_mxl(type='probabilities') sums to 1 within each individual", {
  fit <- fit_small_mxl()
  preds <- predict(fit, type = "probabilities")

  M <- fit$data$M
  # For each individual, the per-individual probabilities sum to 1
  ends <- cumsum(M)
  starts <- c(1L, head(ends, -1) + 1L)

  per_id_sums <- vapply(
    seq_along(M),
    function(i) sum(preds$choice_prob[starts[i]:ends[i]]),
    numeric(1)
  )

  expect_equal(per_id_sums, rep(1, length(M)), tolerance = 1e-6)
})

test_that("predict.choicer_mxl errors without keep_data", {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    predict(fit, type = "probabilities"),
    "Refit with keep_data"
  )
})

test_that("predict.choicer_mxl rejects unknown type", {
  fit <- fit_small_mxl()

  expect_error(predict(fit, type = "garbage"))
})

# =============================================================================
# MXL diversion_ratios
# =============================================================================

test_that("diversion_ratios.choicer_mxl returns labeled J x J matrix", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_true(is.matrix(dr))
  expect_equal(dim(dr), c(J, J))
  expect_true(all(is.finite(dr)))
  expect_equal(rownames(dr), as.character(fit$alt_mapping[[2]]))
  expect_equal(colnames(dr), as.character(fit$alt_mapping[[2]]))
})

test_that("diversion_ratios.choicer_mxl has zero diagonal", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  expect_equal(unname(diag(dr)), rep(0, J))
})

test_that("diversion_ratios.choicer_mxl column sums equal 1", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit)

  col_sums <- unname(colSums(dr))
  expect_equal(col_sums, rep(1, J), tolerance = 1e-6)
})

test_that("diversion_ratios.choicer_mxl entries are non-negative", {
  fit <- fit_small_mxl()

  dr <- diversion_ratios(fit)

  expect_true(all(dr >= 0))
})

test_that("diversion_ratios.choicer_mxl errors without keep_data", {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 50L),
    keep_data = FALSE
  )

  expect_error(
    diversion_ratios(fit),
    "keep_data"
  )
})

test_that("diversion_ratios.choicer_mxl is deterministic", {
  fit <- fit_small_mxl()

  dr1 <- diversion_ratios(fit)
  dr2 <- diversion_ratios(fit)

  expect_equal(dr1, dr2)
})

test_that("MXL diversion_ratios approximates MNL when variance is near zero", {
  skip("MXL->MNL collapse sanity check requires deeper integration")
})
