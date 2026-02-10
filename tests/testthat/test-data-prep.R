# Tests for data preparation functions:
# - prepare_mnl_data()
# - prepare_mxl_data()
# - prepare_nl_data()
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

  # Warns about removed rows before erroring
  expect_error(
    suppressWarnings(prepare_mnl_data(dt, "id", "alt", "choice", "x1")),
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

# --- S3 class checks for prepare_*_data() ---

test_that("prepare_mnl_data returns choicer_data_mnl class", {
  dt <- create_small_mnl_data()
  result <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  expect_s3_class(result, "choicer_data_mnl")
  expect_true(is.list(result))
  expect_true(!is.null(result$data_spec))
  expect_equal(result$data_spec$id_col, "id")
  expect_equal(result$data_spec$alt_col, "alt")
})

test_that("prepare_mxl_data returns choicer_data_mxl class", {
  dt <- create_small_mxl_data()
  result <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )

  expect_s3_class(result, "choicer_data_mxl")
  expect_true(is.list(result))
  expect_true(!is.null(result$data_spec))
  expect_equal(result$data_spec$random_var_cols, c("w1", "w2"))
})

test_that("prepare_nl_data returns choicer_data_nl class", {
  dt <- create_small_nl_data()
  result <- prepare_nl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    nest_col = "nest"
  )

  expect_s3_class(result, "choicer_data_nl")
  expect_true(is.list(result))
  expect_true(!is.null(result$nest_idx))
  expect_true(!is.null(result$data_spec))
  expect_equal(result$data_spec$nest_col, "nest")
})

# --- prepare_nl_data tests ---

test_that("prepare_nl_data returns correct structure", {
  dt <- create_small_nl_data()
  result <- prepare_nl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    nest_col = "nest"
  )

  expect_true("X" %in% names(result))
  expect_true("alt_idx" %in% names(result))
  expect_true("choice_idx" %in% names(result))
  expect_true("M" %in% names(result))
  expect_true("nest_idx" %in% names(result))
  expect_true("alt_mapping" %in% names(result))

  # nest_idx has length J (number of alternatives)
  J <- nrow(result$alt_mapping)
  expect_length(result$nest_idx, J)

  # 2 nests in test data
  expect_equal(length(unique(result$nest_idx)), 2)
})

test_that("prepare_nl_data validates missing nest_col", {
  dt <- create_small_nl_data()

  expect_error(
    prepare_nl_data(
      dt, "id", "alt", "choice", c("x1", "x2"),
      nest_col = "nonexistent"
    ),
    "Missing column"
  )
})

test_that("prepare_nl_data validates alternatives in multiple nests", {
  dt <- create_small_nl_data()
  # Create conflicting nest assignments: alt 1 in both nests
  dt[alt == 1 & id == 1, nest := 2L]

  expect_error(
    prepare_nl_data(
      dt, "id", "alt", "choice", c("x1", "x2"),
      nest_col = "nest"
    ),
    "multiple nests"
  )
})

test_that("prepare_nl_data validates at least 2 nests", {
  dt <- create_small_nl_data()
  dt[, nest := 1L]  # All in same nest

  expect_error(
    prepare_nl_data(
      dt, "id", "alt", "choice", c("x1", "x2"),
      nest_col = "nest"
    ),
    "At least 2 nests"
  )
})

test_that("prepare_nl_data validates no NA nest assignments", {
  dt <- create_small_nl_data()
  dt[alt == 1, nest := NA_integer_]

  expect_error(
    prepare_nl_data(
      dt, "id", "alt", "choice", c("x1", "x2"),
      nest_col = "nest"
    ),
    "Missing nest assignments"
  )
})

test_that("prepare_nl_data output is compatible with run_nestlogit()", {
  dt <- create_small_nl_data()
  nl_data <- prepare_nl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    nest_col = "nest"
  )

  # Should have all fields needed by run_nestlogit(input_data = ...)
  expect_true(!is.null(nl_data$X))
  expect_true(!is.null(nl_data$alt_idx))
  expect_true(!is.null(nl_data$choice_idx))
  expect_true(!is.null(nl_data$nest_idx))
  expect_true(!is.null(nl_data$M))
  expect_true(!is.null(nl_data$weights))
  expect_true(!is.null(nl_data$alt_mapping))
  expect_true(!is.null(nl_data$include_outside_option))
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
