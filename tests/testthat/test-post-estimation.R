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

# Fit a small NL model (N=30, J=6, 2 nests of size 3) via the convenience
# pathway. create_small_nl_data() returns a raw data.table with columns
# id, alt, nest, x1, x2, choice.
fit_small_nl <- function(seed = 123, keep_data = TRUE) {
  dt <- create_small_nl_data(seed = seed)
  run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    control = list(maxeval = 50L),
    keep_data = keep_data
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

  dr <- diversion_ratios(fit, wrt_var = "x1")

  expect_true(is.matrix(dr))
  expect_equal(dim(dr), c(J, J))
  expect_true(all(is.finite(dr)))
  expect_equal(rownames(dr), as.character(fit$alt_mapping[[2]]))
  expect_equal(colnames(dr), as.character(fit$alt_mapping[[2]]))
})

test_that("diversion_ratios.choicer_mxl has zero diagonal", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit, wrt_var = "x1")

  expect_equal(unname(diag(dr)), rep(0, J))
})

test_that("diversion_ratios.choicer_mxl column sums equal 1", {
  fit <- fit_small_mxl()
  J <- nrow(fit$alt_mapping)

  dr <- diversion_ratios(fit, wrt_var = "x1")

  col_sums <- unname(colSums(dr))
  expect_equal(col_sums, rep(1, J), tolerance = 1e-6)
})

test_that("diversion_ratios.choicer_mxl entries are non-negative (fixed-coef variable)", {
  # For a fixed-coef variable, beta cancels in the ratio, so DR >= 0 always
  fit <- fit_small_mxl()

  dr <- diversion_ratios(fit, wrt_var = "x1")

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
    diversion_ratios(fit, wrt_var = "x1"),
    "keep_data"
  )
})

test_that("diversion_ratios.choicer_mxl is deterministic", {
  fit <- fit_small_mxl()

  dr1 <- diversion_ratios(fit, wrt_var = "x1")
  dr2 <- diversion_ratios(fit, wrt_var = "x1")

  expect_equal(dr1, dr2)
})

test_that("diversion_ratios.choicer_mxl depends on which random-coef variable is perturbed", {
  # With random coefficients beta_{ik}^s varies across draws and does not
  # cancel in the ratio. The diversion matrix should differ between w1 and w2.
  fit <- fit_small_mxl()

  dr_w1 <- diversion_ratios(fit, wrt_var = "w1", is_random_coef = TRUE)
  dr_w2 <- diversion_ratios(fit, wrt_var = "w2", is_random_coef = TRUE)

  # Both still satisfy the column-sum identity
  expect_equal(unname(colSums(dr_w1)), rep(1, nrow(dr_w1)), tolerance = 1e-6)
  expect_equal(unname(colSums(dr_w2)), rep(1, nrow(dr_w2)), tolerance = 1e-6)

  # But the matrices themselves should differ (not equal up to numerical noise)
  expect_false(isTRUE(all.equal(dr_w1, dr_w2, tolerance = 1e-3)))
})

test_that("MXL diversion_ratios collapses to MNL when sigma is near zero", {
  # When Sigma -> 0, beta_{ik}^s = mu_k* + gamma_{ik}^{s*} -> mu_k* (a constant
  # across draws and individuals), so beta cancels in the ratio and the MXL
  # diversion matrix should match the MNL diversion matrix computed on the
  # same dataset with the same effective coefficients.
  dt <- create_small_mxl_data()

  fit_mnl <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "w1"),
    control = list(maxeval = 50L)
  )

  fit_mxl <- run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1",
    rc_correlation = FALSE, rc_mean = TRUE,
    S = 20L, control = list(maxeval = 20L)
  )

  # Force MXL parameters to match MNL with near-zero variance.
  # MXL layout (rc_mean=TRUE, rc_correlation=FALSE, 1 random coef):
  #   theta = [beta_x1, mu_w1, L_11_param, ASC_2, ASC_3]
  # MNL layout: [x1, w1, ASC_2, ASC_3].
  pm <- fit_mxl$param_map
  mnl_coef <- coef(fit_mnl)
  fit_mxl$coefficients[pm$beta] <- mnl_coef["x1"]
  fit_mxl$coefficients[pm$mu] <- mnl_coef["w1"]
  fit_mxl$coefficients[pm$sigma] <- log(1e-8)  # L_diag = 1e-8 -> Sigma ~= 1e-16
  fit_mxl$coefficients[pm$asc] <- mnl_coef[c("ASC_2", "ASC_3")]

  dr_mnl <- diversion_ratios(fit_mnl)
  dr_mxl_x1 <- diversion_ratios(fit_mxl, wrt_var = "x1")
  dr_mxl_w1 <- diversion_ratios(fit_mxl, wrt_var = "w1", is_random_coef = TRUE)

  # Strip dimnames for the numeric comparison (both have the same labels but
  # all.equal is sensitive to attributes)
  strip <- function(m) {
    out <- unclass(m)
    dimnames(out) <- NULL
    out
  }

  # Fixed-coef wrt_var: beta cancels exactly, so should match to machine eps
  expect_equal(strip(dr_mxl_x1), strip(dr_mnl), tolerance = 1e-6)
  # Random-coef wrt_var with Sigma ~ 0: beta_{ik}^s ~ mu_w1, also cancels
  expect_equal(strip(dr_mxl_w1), strip(dr_mnl), tolerance = 1e-6)
})

# =============================================================================
# NL S3 dispatch
# =============================================================================
# The four post-estimation generics must dispatch on class choicer_nl and
# return objects shaped like their MNL counterparts. Detailed numeric checks
# live in test-predictions.R, test-elasticities.R, and test-blp.R; here we
# confirm dispatch, the keep_data = FALSE guard, and bad-variable handling.

test_that("post-estimation generics dispatch on a fitted choicer_nl", {
  fit <- fit_small_nl()
  expect_s3_class(fit, "choicer_nl")

  J <- nrow(fit$alt_mapping)

  preds <- predict(fit, type = "probabilities")
  expect_type(preds, "list")
  expect_true("choice_prob" %in% names(preds))

  elast <- elasticities(fit, elast_var = "x1")
  expect_true(is.matrix(elast))
  expect_equal(dim(elast), c(J, J))

  dr <- diversion_ratios(fit)
  expect_true(is.matrix(dr))
  expect_equal(dim(dr), c(J, J))

  delta <- blp(fit, target_shares = rep(1 / J, J))
  expect_true(is.numeric(delta))
})

# --- keep_data = FALSE guards ------------------------------------------------

test_that("predict.choicer_nl errors without data", {
  fit <- fit_small_nl(keep_data = FALSE)

  expect_error(
    predict(fit, type = "probabilities"),
    "keep_data"
  )
})

test_that("elasticities.choicer_nl errors without data", {
  fit <- fit_small_nl(keep_data = FALSE)

  expect_error(
    elasticities(fit, elast_var = "x1"),
    "keep_data"
  )
})

test_that("diversion_ratios.choicer_nl errors without data", {
  fit <- fit_small_nl(keep_data = FALSE)

  expect_error(
    diversion_ratios(fit),
    "keep_data"
  )
})

test_that("blp.choicer_nl errors without data", {
  fit <- fit_small_nl(keep_data = FALSE)
  J <- nrow(fit$alt_mapping)

  expect_error(
    blp(fit, target_shares = rep(1 / J, J)),
    "keep_data"
  )
})

# --- bad-variable handling ---------------------------------------------------

test_that("elasticities.choicer_nl errors with invalid variable name", {
  fit <- fit_small_nl()

  expect_error(
    elasticities(fit, elast_var = "nonexistent"),
    "not found"
  )
})

test_that("diversion_ratios.choicer_nl takes no variable argument (MNL-style)", {
  # diversion_ratios.choicer_nl mirrors the MNL signature: diversion_ratios(fit).
  # It must not require a wrt_var/elast_var argument.
  fit <- fit_small_nl()

  dr <- diversion_ratios(fit)
  expect_true(is.matrix(dr))
})
