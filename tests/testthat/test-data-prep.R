# Tests for data preparation functions:
# - prepare_mnl_data()
# - prepare_mxl_data()
# - check_collinearity()
# - remove_nullspace_cols()

# --- prepare_mnl_data tests ---

test_that("prepare_mnl_data validates required columns", {
  dt <- data.table(
    id = rep(1:3, each = 2),
    alt = rep(1:2, 3),
    choice = rep(c(1L, 0L), 3),
    x1 = rnorm(6)
  )

  # Should work with correct columns
  expect_no_error(
    prepare_mnl_data(dt, "id", "alt", "choice", "x1")
  )

  # Should error with missing column
  expect_error(
    prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2")),
    "Missing columns"
  )

  # Should error with wrong id column
  expect_error(
    prepare_mnl_data(dt, "wrong_id", "alt", "choice", "x1"),
    "Missing columns"
  )
})

test_that("prepare_mnl_data handles missing values with warning", {
  dt <- create_small_mnl_data()
  dt[1, x1 := NA]

  expect_warning(
    prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2")),
    "choice situations containing missing values"
  )
})

test_that("prepare_mnl_data errors when all data has NA", {
  dt <- data.table(
    id = rep(1:2, each = 2),
    alt = rep(1:2, 2),
    choice = rep(c(1L, 0L), 2),
    x1 = NA_real_
  )

  expect_error(
    prepare_mnl_data(dt, "id", "alt", "choice", "x1"),
    "All choice situations removed"
  )
})

test_that("prepare_mnl_data validates choice column format", {
  # Multiple choices per individual
  dt <- data.table(
    id = rep(1:3, each = 2),
    alt = rep(1:2, 3),
    choice = c(1L, 1L, 1L, 0L, 0L, 1L),  # id 1 has two choices
    x1 = rnorm(6)
  )

  expect_error(
    prepare_mnl_data(dt, "id", "alt", "choice", "x1", include_outside_option = FALSE),
    "exactly one chosen"
  )

  # No choice for an individual
  dt2 <- data.table(
    id = rep(1:3, each = 2),
    alt = rep(1:2, 3),
    choice = c(0L, 0L, 1L, 0L, 0L, 1L),  # id 1 has no choice
    x1 = rnorm(6)
  )

  expect_error(
    prepare_mnl_data(dt2, "id", "alt", "choice", "x1", include_outside_option = FALSE),
    "exactly one chosen"
  )
})

test_that("prepare_mnl_data validates numeric covariates", {
  dt <- data.table(
    id = rep(1:3, each = 2),
    alt = rep(1:2, 3),
    choice = rep(c(1L, 0L), 3),
    x1 = c("a", "b", "c", "d", "e", "f")  # Character covariate
  )

  expect_error(
    prepare_mnl_data(dt, "id", "alt", "choice", "x1"),
    "must be numeric"
  )
})

test_that("prepare_mnl_data returns correct output structure", {
  dt <- create_small_mnl_data()
  result <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  # Check required outputs exist
  expect_true("X" %in% names(result))
  expect_true("alt_idx" %in% names(result))
  expect_true("choice_idx" %in% names(result))
  expect_true("M" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("alt_mapping" %in% names(result))
  expect_true("N" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$X), nrow(dt))
  expect_equal(ncol(result$X), 2)  # x1 and x2
  expect_equal(length(result$M), result$N)
  expect_equal(length(result$choice_idx), result$N)
})

test_that("prepare_mnl_data handles outside option correctly", {
  dt <- data.table(
    id = rep(1:5, each = 3),
    alt = rep(c(0L, 1L, 2L), 5),  # 0 is outside option
    choice = rep(c(0L, 1L, 0L), 5),
    x1 = rnorm(15)
  )

  result <- prepare_mnl_data(
    dt, "id", "alt", "choice", "x1",
    outside_opt_label = 0L,
    include_outside_option = TRUE
  )

  # Outside option should be first in alt_mapping
  expect_equal(result$alt_mapping$alt[1], 0L)
  expect_true(result$include_outside_option)
})

# --- prepare_mxl_data tests ---

test_that("prepare_mxl_data validates required columns", {
  dt <- create_small_mxl_data()

  # Should work with correct columns
  expect_no_error(
    prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
  )

  # Should error with missing random var column
  expect_error(
    prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w3")),
    "Missing columns"
  )
})

test_that("prepare_mxl_data returns correct output structure", {
  dt <- create_small_mxl_data()
  result <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  # Check required outputs
  expect_true("X" %in% names(result))
  expect_true("W" %in% names(result))
  expect_true("rc_correlation" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$X), nrow(dt))
  expect_equal(ncol(result$X), 1)  # x1 only
  expect_equal(ncol(result$W), 2)  # w1 and w2
  expect_false(result$rc_correlation)
})

test_that("prepare_mxl_data handles rc_correlation flag", {
  dt <- create_small_mxl_data()

  result_uncorr <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  result_corr <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = TRUE
  )

  expect_false(result_uncorr$rc_correlation)
  expect_true(result_corr$rc_correlation)
})

# --- check_collinearity tests ---

test_that("check_collinearity detects perfectly collinear columns", {
  # Create matrix with collinear columns: col3 = 2*col1
  X <- matrix(c(
    1, 2, 2,
    2, 3, 4,
    3, 4, 6,
    4, 5, 8
  ), nrow = 4, byrow = TRUE)
  colnames(X) <- c("a", "b", "c")

  result <- check_collinearity(X)

  # Should drop one of the collinear columns
  expect_equal(ncol(result$mat), 2)
  expect_true(length(result$dropped) == 1)
  expect_true("c" %in% result$dropped || "a" %in% result$dropped)
})

test_that("check_collinearity keeps independent columns", {
  X <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow = 4, byrow = TRUE)
  colnames(X) <- c("a", "b", "c")

  result <- check_collinearity(X)

  expect_equal(ncol(result$mat), 3)
  expect_equal(length(result$dropped), 0)
})

test_that("check_collinearity handles single column", {
  X <- matrix(1:4, ncol = 1)
  colnames(X) <- "a"

  result <- check_collinearity(X)

  expect_equal(ncol(result$mat), 1)
  expect_equal(length(result$dropped), 0)
})
