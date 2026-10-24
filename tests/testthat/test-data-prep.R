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

# --- the preps copy only the columns they use ---

prep_all <- function(d) {
  list(mnl = prepare_mnl_data(d, "id", "alt", "choice", c("x1", "x2")),
       mxl = prepare_mxl_data(d, "id", "alt", "choice", "x1", "x2"),
       nl = prepare_nl_data(d, "id", "alt", "choice", c("x1", "x2"), "nest"))
}

test_that("preparing data leaves the caller's data unchanged", {
  df <- as.data.frame(create_small_nl_data())
  set.seed(1)
  df <- df[sample(nrow(df)), ]         # shuffled, so the preps reorder rows
  dt <- data.table::as.data.table(df)
  data.table::setkeyv(dt, "x1")        # a key the preps do not sort by
  data.table::setindexv(dt, "alt")
  df0 <- data.table::copy(df)
  dt0 <- data.table::copy(dt)
  prep_all(df)
  prep_all(dt)
  expect_identical(df, df0)
  expect_identical(dt, dt0)            # values, row order, key and index
})

test_that("columns a model does not use do not affect its preparation", {
  df <- as.data.frame(create_small_nl_data())
  extra <- df
  extra$note <- c(NA, "a")                      # NAs in an unused column
  extra$lst <- as.list(seq_len(nrow(extra)))    # a list column
  extra$x_unused <- NA_real_
  expect_silent(p_extra <- prep_all(extra))
  expect_identical(p_extra, prep_all(df))
  # Inputs outside the direct route keep the full-copy route.
  expect_identical(prep_all(as.list(df)), prep_all(df))
  with_matrix <- df
  with_matrix$m <- matrix(1, nrow(df), 2)
  expect_identical(prep_all(with_matrix), prep_all(df))
})

test_that("a missing value in a used column drops its whole choice situation", {
  df <- as.data.frame(create_small_nl_data())
  df$x2[df$id == 3][2] <- NA
  df$x1[df$id == 7][1] <- NA
  expect_warning(p <- prepare_mxl_data(df, "id", "alt", "choice", "x1", "x2"),
                 "Removed 2 choice situations containing missing values")
  expect_false(any(p$situation_ids %in% c(3, 7)))
  expect_equal(p$N, 28)
  # x2 is not used here, so situation 3 stays.
  expect_warning(q <- prepare_mnl_data(df, "id", "alt", "choice", "x1"),
                 "Removed 1 choice situations containing missing values")
  expect_false(7 %in% q$situation_ids)
  expect_true(3 %in% q$situation_ids)
})

test_that("an infinite covariate value drops its choice situation, naming the column", {
  df <- as.data.frame(create_small_nl_data())
  df$person <- (df$id - 1L) %/% 3L + 1L
  df$x2[df$id == 3][2] <- Inf
  df$x1[df$id == 7][1] <- -Inf
  kept <- df[!df$id %in% c(3, 7), ]
  msg <- function(n, cols) {
    paste0("Removed ", n, " choice situations containing non-finite ",
           "covariate values (in ", cols, ").")
  }
  # Each preparation drops them, naming the covariates that hold the values,
  # and returns what the data without those situations give
  expect_warning(p <- prepare_mnl_data(df, "id", "alt", "choice", c("x1", "x2")),
                 msg(2, "x1, x2"), fixed = TRUE)
  expect_identical(p, prepare_mnl_data(kept, "id", "alt", "choice", c("x1", "x2")))
  expect_warning(p <- prepare_nl_data(df, "id", "alt", "choice", c("x1", "x2"),
                                      "nest"),
                 msg(2, "x1, x2"), fixed = TRUE)
  expect_identical(p, prepare_nl_data(kept, "id", "alt", "choice", c("x1", "x2"),
                                      "nest"))
  # Fixed and random coefficients, in the cross-section and in a panel
  expect_warning(p <- prepare_mxl_data(df, "id", "alt", "choice", "x1", "x2"),
                 msg(2, "x1, x2"), fixed = TRUE)
  expect_identical(p, prepare_mxl_data(kept, "id", "alt", "choice", "x1", "x2"))
  expect_warning(p <- prepare_mxl_data(df, "id", "alt", "choice", "x1", "x2",
                                       person_col = "person"),
                 msg(2, "x1, x2"), fixed = TRUE)
  expect_identical(p, prepare_mxl_data(kept, "id", "alt", "choice", "x1", "x2",
                                       person_col = "person"))
  # Only the covariates the model uses count
  expect_warning(p <- prepare_mnl_data(df, "id", "alt", "choice", "x1"),
                 msg(1, "x1"), fixed = TRUE)
  expect_true(3 %in% p$situation_ids)
  # A NaN is a missing value, dropped as before
  df$x1[df$id == 9][1] <- NaN
  expect_warning(
    expect_warning(prepare_mnl_data(df, "id", "alt", "choice", "x1"),
                   "Removed 1 choice situations containing missing values",
                   fixed = TRUE),
    msg(1, "x1"), fixed = TRUE)
})

test_that("an infinite value in a row the preparation removes does not count", {
  df <- as.data.frame(create_small_nl_data())
  # The outside option's own rows are removed before the covariates are read
  df$x1[df$alt == 6][1:2] <- Inf
  expect_no_warning(p <- prepare_mnl_data(
    df, "id", "alt", "choice", "x1", include_outside_option = TRUE,
    outside_opt_label = 6L))
  expect_equal(p$N, 30)
})

test_that("a preparation that would drop every choice situation stops", {
  dt <- create_small_mnl_data()
  dt[alt == 1, x1 := Inf]
  expect_error(
    suppressWarnings(prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))),
    "All choice situations removed due to non-finite covariate values.",
    fixed = TRUE)
  expect_error(
    suppressWarnings(prepare_mxl_data(dt, "id", "alt", "choice", "x2", "x1")),
    "All choice situations removed due to non-finite covariate values.",
    fixed = TRUE)
  # A base label whose alternative only a dropped situation offered
  d <- data.frame(id = rep(1:4, each = 2L), alt = c(1, 2, 1, 2, 1, 2, 1, 3),
                  choice = rep(c(1L, 0L), 4L), x1 = c(1:7, Inf))
  expect_error(
    suppressWarnings(prepare_mnl_data(d, "id", "alt", "choice", "x1",
                                      outside_opt_label = 3)),
    paste("that remain after removing the 1 choice situation with non-finite",
          "covariate values"), fixed = TRUE)
})

test_that(".rows_not_finite() flags the rows in use and names their columns", {
  x <- data.frame(a = c(1, Inf, 3), b = c(-Inf, 2, 3), c = c(1, 2, 3))
  # b's -Inf is in a row no longer in use
  expect_identical(.rows_not_finite(x, 1:3, 2:3),
                   list(rows = c(TRUE, FALSE), cols = "a"))
  expect_identical(.rows_not_finite(x, 2:3, 2:3),
                   list(rows = c(FALSE, FALSE), cols = character(0)))
  expect_identical(.rows_not_finite(x, 1:3, 1:3),
                   list(rows = c(TRUE, TRUE, FALSE), cols = c("a", "b")))
  # A column whose other values are finite is not scanned for its NAs
  x$d <- c(1, NA, 3)
  expect_identical(.rows_not_finite(x, 3:4, 1:3),
                   list(rows = c(FALSE, FALSE, FALSE), cols = character(0)))
})

test_that("situations dropped for missing values, then infinite ones, keep the rows aligned", {
  df <- as.data.frame(create_small_nl_data())
  df$x1[df$id == 4][2] <- NA
  df$x2[df$id == 11][1] <- Inf
  df$x1[df$id == 20][3] <- -Inf
  kept <- df[!df$id %in% c(4, 11, 20), ]
  expect_identical(
    suppressWarnings(prepare_mnl_data(df, "id", "alt", "choice", c("x1", "x2"))),
    prepare_mnl_data(kept, "id", "alt", "choice", c("x1", "x2")))
  expect_identical(
    suppressWarnings(prepare_mxl_data(df, "id", "alt", "choice", "x1", "x2")),
    prepare_mxl_data(kept, "id", "alt", "choice", "x1", "x2"))
  dt <- create_small_mnl_data()
  dt[id == 2 & alt == 3, x1 := NA]
  dt[id == 6 & alt == 1, x2 := Inf]   # the base alternative's row
  expect_identical(
    suppressWarnings(prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))),
    prepare_mnp_data(dt[!id %in% c(2, 6)], "id", "alt", "choice", c("x1", "x2")))
})

test_that("a covariate that is both fixed and random is named once", {
  df <- as.data.frame(create_small_nl_data())
  df$x1[df$id == 3][1] <- Inf
  expect_warning(
    p <- prepare_mxl_data(df, "id", "alt", "choice", c("x1", "x2"), "x1"),
    "Removed 1 choice situations containing non-finite covariate values (in x1).",
    fixed = TRUE)
  expect_false(3 %in% p$situation_ids)
})

test_that("the warning names the first ten covariates", {
  expect_warning(
    .warn_not_finite(4L, paste0("x", 1:12)),
    paste("Removed 4 choice situations containing non-finite covariate values",
          "(in x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 and 2 more)."),
    fixed = TRUE)
})

test_that("a panel's situation ids are checked before infinite values drop them", {
  # Ids that restart within decision makers: dropping id 1 would also drop
  # the second decision maker's first situation
  p <- data.frame(person = rep(1:2, each = 4L), id = c(1, 1, 2, 2, 1, 1, 3, 3),
                  alt = rep(1:2, 4L), choice = rep(c(1L, 0L), 4L),
                  x1 = c(Inf, 2:8), w1 = c(0.3, -1, 0.2, 1, 2, -0.5, 0.1, 0.7))
  expect_error(prepare_mxl_data(p, "id", "alt", "choice", "x1", "w1",
                                person_col = "person"),
               "Each 'id' must belong to exactly one 'person'", fixed = TRUE)
})

test_that("the base label's error counts the situations dropped for each reason", {
  # Alternative 3 is offered only in situations 4 and 5, both dropped
  d <- data.frame(id = rep(1:5, each = 2L),
                  alt = c(1, 2, 1, 2, 1, 2, 1, 3, 1, 3),
                  choice = rep(c(1L, 0L), 5L), x1 = as.numeric(1:10),
                  w1 = as.numeric(1:10))
  prep <- function(d) {
    suppressWarnings(prepare_mxl_data(d, "id", "alt", "choice", "x1", "w1",
                                      outside_opt_label = 3))
  }
  removed <- function(why) {
    paste("that remain after removing the 2 choice situations with", why)
  }
  both <- d
  both$x1[8] <- Inf
  both$w1[10] <- NA
  expect_error(prep(both),
               removed("missing values or non-finite covariate values"),
               fixed = TRUE)
  inf <- d
  inf$x1[c(8, 10)] <- Inf
  expect_error(prep(inf), removed("non-finite covariate values"), fixed = TRUE)
  na <- d
  na$w1[c(8, 10)] <- NA
  expect_error(prep(na), removed("missing values:"), fixed = TRUE)
})

test_that("design matrices are double, with the covariate names as column names", {
  df <- as.data.frame(create_small_nl_data())
  df$k1 <- rep_len(0:3, nrow(df))      # integer covariates
  df$k2 <- rep_len(c(1L, 0L, 0L), nrow(df))
  num <- function(cols) {
    matrix(as.numeric(unlist(df[cols])), ncol = length(cols),
           dimnames = list(NULL, cols))
  }
  p <- prepare_mxl_data(df, "id", "alt", "choice", c("k1", "k2"), c("k1", "x1"))
  expect_identical(p$X, num(c("k1", "k2")))
  expect_identical(p$W, num(c("k1", "x1")))
  expect_identical(prepare_mnl_data(df, "id", "alt", "choice", c("x1", "k2"))$X,
                   num(c("x1", "k2")))
  expect_identical(prepare_nl_data(df, "id", "alt", "choice", "k1", "nest")$X,
                   num("k1"))
  # A named covariate vector leaves no names on the column names, and a
  # repeated covariate gives a repeated column, as as.matrix() did.
  expect_identical(prepare_mnl_data(df, "id", "alt", "choice", c(a = "x1", b = "k2"))$X,
                   num(c("x1", "k2")))
  expect_identical(prepare_mxl_data(df, "id", "alt", "choice", "x1", c(w = "k1"))$W,
                   num("k1"))
  expect_identical(.gather_matrix(df, c("x1", "x1"), seq_len(nrow(df))),
                   num(c("x1", "x1")))
  cc <- character(0)
  expect_identical(.gather_matrix(df, cc, seq_len(nrow(df))),
                   as.matrix(data.table::as.data.table(df)[, ..cc]))
})

test_that("prep_gather_design gathers rows exactly and checks its inputs", {
  cols <- list(c(1.5, NA, 3.25, -1), c(4L, NA, 6L, 7L))
  expect_identical(prep_gather_design(cols, c(4L, 1L, 2L)),
                   matrix(c(-1, 1.5, NA, 7, 4, NA), 3, 2))
  expect_identical(prep_gather_design(cols, integer(0)), matrix(numeric(0), 0, 2))
  expect_error(prep_gather_design(cols, c(1L, 5L)), "out of range")
  expect_error(prep_gather_design(cols, c(0L, 1L)), "out of range")
  expect_error(prep_gather_design(cols, c(NA_integer_, 1L)), "out of range")
  expect_error(prep_gather_design(list(letters[1:4]), 1L), "neither integer nor double")
})

test_that("gathered design rows follow the prepared order through filters and sorts", {
  set.seed(3)
  df <- as.data.frame(create_small_nl_data())
  df$k <- rep_len(0:4, nrow(df))                   # an integer covariate
  df$person <- (df$id - 1L) %/% 3L + 1L            # three situations per decision maker
  outside <- df[!duplicated(df$id), ]
  outside$alt <- 0L
  outside$choice <- 0L
  full <- rbind(df, outside)                       # physical outside-option rows
  full$x2[full$id == 5 & full$alt == 2] <- NA      # drops situation 5
  full <- full[sample(nrow(full)), ]               # shuffled
  inside <- full[full$alt != 0 & full$id != 5, ]
  expected <- function(cols, ...) {
    m <- as.matrix(inside[order(...), cols, drop = FALSE])
    storage.mode(m) <- "double"
    dimnames(m) <- list(NULL, cols)
    m
  }
  expect_warning(
    p <- prepare_mnl_data(full, "id", "alt", "choice", c("x1", "k", "x2"),
                          outside_opt_label = 0L, include_outside_option = TRUE),
    "Removed 1 choice situations")
  expect_identical(p$X, expected(c("x1", "k", "x2"), inside$id, inside$alt))
  expect_warning(
    q <- prepare_mxl_data(full[full$alt != 0, ], "id", "alt", "choice",
                          c("x1", "k"), "x2", person_col = "person"),
    "Removed 1 choice situations")
  expect_identical(q$X, expected(c("x1", "k"), inside$person, inside$id, inside$alt))
  expect_identical(q$W, expected("x2", inside$person, inside$id, inside$alt))
})

test_that("situation-level columns are read correctly for ids named V1 or N", {
  df <- as.data.frame(create_small_nl_data())      # ids 1..30, six rows each
  df$person <- (df$id - 1L) %/% 3L + 1L            # three situations per person
  df$wt <- df$id / 10
  df$cl <- df$id %% 4L
  df$pw <- df$person / 10                          # constant within person
  df$pcl <- df$person %% 4L
  for (nm in c("V1", "N")) {
    d <- df
    names(d)[names(d) == "id"] <- nm
    q <- prepare_mnl_data(d, nm, "alt", "choice", "x1", weights_col = "wt",
                          cluster_col = "cl")
    expect_identical(q$weights, (1:30) / 10)
    expect_identical(q$cluster, (1:30) %% 4L)
    expect_identical(q$M, rep(6L, 30))
    p <- prepare_mxl_data(d, nm, "alt", "choice", "x1", "x2", weights_col = "pw",
                          cluster_col = "pcl", person_col = "person")
    expect_identical(p$weights, rep((1:10) / 10, each = 3))
    expect_identical(p$cluster, rep((1:10) %% 4L, each = 3))
    expect_identical(p$Ti, rep(3L, 10))
    expect_identical(p$M, rep(6L, 30))
  }
})

test_that("predict(newdata) builds the same choice sets for ids named V1 or N", {
  df <- as.data.frame(create_small_nl_data())      # ids 1..30, six rows each
  df$person <- (df$id - 1L) %/% 3L + 1L
  for (nm in c("V1", "N")) {
    d <- df
    names(d)[names(d) == "id"] <- nm
    fits <- suppressMessages(list(
      mnl = run_mnlogit(d, nm, "alt", "choice", "x1"),
      mxl = run_mxlogit(d, nm, "alt", "choice", "x1", "x2", S = 20L,
                        person_col = "person")
    ))
    for (f in fits) {
      expect_identical(prepare_newdata(f, d)$M, rep(6L, 30))
      expect_equal(predict(f, newdata = d), predict(f), tolerance = 1e-12)
    }
  }
})

test_that("covariates may carry the names of the working columns", {
  # alt_int and idx_in_group are columns of the preps' working table; the
  # covariates are read from the data, so they no longer collide.
  df <- as.data.frame(create_small_nl_data())      # sorted by id, alt; no NA
  df$alt_int <- 2 * df$x1
  df$idx_in_group <- df$x2 + 1
  cols <- c("alt_int", "idx_in_group")
  ref <- as.matrix(df[cols])
  expect_identical(prepare_mnl_data(df, "id", "alt", "choice", cols)$X, ref)
  expect_identical(prepare_nl_data(df, "id", "alt", "choice", cols, "nest")$X,
                   ref)
  expect_identical(prepare_mxl_data(df, "id", "alt", "choice", "x1", cols)$W,
                   ref)
})

test_that("a repeated covariate stops MNL and NL preparation, as before", {
  df <- as.data.frame(create_small_nl_data())
  expect_error(prepare_mnl_data(df, "id", "alt", "choice", c("x1", "x1")),
               "names a column more than once: x1")
  expect_error(prepare_nl_data(df, "id", "alt", "choice", c("x1", "x1"), "nest"),
               "names a column more than once: x1")
})

test_that("a covariate that is also an index column stops MNL, MXL and NL preparation", {
  # It left their working table before they read it again: they stopped
  # with an unrelated error, or counted another object named like the choice
  # column. They now say why.
  dt <- create_small_nl_data()
  dt[, `:=`(w = 1 + id %% 3, cl = id %% 4, pid = (id + 1) %/% 2)]
  msg <- function(col) {
    paste0("A covariate cannot also be the id, alternative, choice, weight, ",
           "cluster or decision-maker column: '", col, "'.")
  }
  for (col in c("id", "alt", "choice", "w", "cl")) {
    expect_error(prepare_mnl_data(dt, "id", "alt", "choice", c("x1", col),
                                  weights_col = "w", cluster_col = "cl"),
                 msg(col), fixed = TRUE)
    expect_error(prepare_nl_data(dt, "id", "alt", "choice", c(col, "x1"),
                                 "nest", weights_col = "w", cluster_col = "cl"),
                 msg(col), fixed = TRUE)
  }
  for (col in c("id", "alt", "choice", "pid")) {
    expect_error(prepare_mxl_data(dt, "id", "alt", "choice", c("x1", col),
                                  "x2", person_col = "pid"),
                 msg(col), fixed = TRUE)
    expect_error(prepare_mxl_data(dt, "id", "alt", "choice", "x1", col,
                                  person_col = "pid"),
                 msg(col), fixed = TRUE)
  }
  # A weight or cluster column the call does not use is an ordinary covariate
  expect_identical(colnames(prepare_mnl_data(dt, "id", "alt", "choice",
                                             c("x1", "w", "cl"))$X),
                   c("x1", "w", "cl"))
})

test_that("a character or factor choice column is reported as not 0/1", {
  # "0"/"1" pass a %in% c(0, 1) check by coercion, but the choices cannot be
  # counted; the preparations and the WESML helpers now say so instead of
  # failing in sum().
  dt <- create_small_mnl_data()
  dt[, nest := ifelse(alt == 1L, "a", "b")]
  for (as_type in list(as.character, factor)) {
    d <- copy(dt)[, choice := as_type(choice)]
    msg <- "`choice` must contain only 0 and 1."
    expect_error(prepare_mnl_data(d, "id", "alt", "choice", "x1"), msg,
                 fixed = TRUE)
    expect_error(prepare_mxl_data(d, "id", "alt", "choice", "x1", "x2"), msg,
                 fixed = TRUE)
    expect_error(prepare_nl_data(d, "id", "alt", "choice", "x1", "nest"), msg,
                 fixed = TRUE)
    expect_error(prepare_mnp_data(d, "id", "alt", "choice", "x1"), msg,
                 fixed = TRUE)
    expect_error(prepare_hmnl_data(d, "id", "alt", "choice", "x1"), msg,
                 fixed = TRUE)
    expect_error(prepare_hmnp_data(d, "id", "alt", "choice", "x1"), msg,
                 fixed = TRUE)
    expect_error(wesml_weights(d, "id", "alt", "choice",
                               Q = c(`1` = 0.3, `2` = 0.3, `3` = 0.4)),
                 msg, fixed = TRUE)
    expect_error(sample_by_choice(d, "id", "alt", "choice", n_per_alt = 2L),
                 msg, fixed = TRUE)
  }
  # Logical choices are counted as before
  d <- copy(dt)[, choice := as.logical(choice)]
  expect_identical(prepare_mnl_data(d, "id", "alt", "choice", "x1")$choice_idx,
                   prepare_mnl_data(dt, "id", "alt", "choice", "x1")$choice_idx)
})

test_that("attributes of the choice column do not reach alt_mapping", {
  # A pdata.frame's or a labelled data set's columns carry attributes. sum()
  # drops them from the choice counts, where data.table's GForce sum would
  # keep them, so the counts in alt_mapping are summed by sum().
  dt <- create_small_mnl_data()
  dt[, nest := ifelse(alt == 1L, "a", "b")]
  d <- copy(dt)[, choice := structure(choice, note = "chosen")]
  preps <- list(
    function(x) prepare_mnl_data(x, "id", "alt", "choice", "x1"),
    function(x) prepare_mxl_data(x, "id", "alt", "choice", "x1", "x2"),
    function(x) prepare_nl_data(x, "id", "alt", "choice", "x1", "nest"),
    function(x) prepare_mnp_data(x, "id", "alt", "choice", "x1"),
    function(x) prepare_hmnl_data(x, "id", "alt", "choice", "x1")
  )
  for (prep in preps) {
    expect_identical(prep(d), prep(dt))
  }
})


test_that("prep_gather_design differences rows against base rows in one pass", {
  cols <- list(c(1.5, NA, 3.25, -1), c(4L, NA, 6L, 7L))
  r <- c(4L, 1L, 3L)
  b <- c(1L, 1L, 2L)
  # As R's arithmetic forms the differences, NA included, in double.
  expect_identical(prep_gather_design(cols, r, base = b),
                   cbind(cols[[1]][r] - cols[[1]][b], cols[[2]][r] - cols[[2]][b]))
  expect_identical(prep_gather_design(cols, integer(0), base = integer(0)),
                   matrix(numeric(0), 0, 2))
  expect_error(prep_gather_design(cols, 1:2, base = 1L), "one index per element")
  expect_error(prep_gather_design(cols, 1:2, base = c(1L, 5L)), "out of range")
  expect_error(prep_gather_design(cols, 1:2, base = c(1L, NA)), "out of range")
  expect_error(prep_gather_design(cols, 1:2, base = c(1, 2)), "NULL or an integer vector")
  # Past the OpenMP threshold.
  set.seed(1)
  n <- 200001L
  big <- list(rnorm(n), sample.int(1000L, n, TRUE))
  r <- sample.int(n)
  b <- sample.int(n)
  expect_identical(prep_gather_design(big, r, base = b),
                   cbind(big[[1]][r] - big[[1]][b], big[[2]][r] - big[[2]][b]))
})

test_that("column arguments that are not names are reported as missing columns", {
  # The reserved-name checks run before the data are read; an argument that
  # is NA or a number is left to the check for missing columns.
  dt <- create_small_mnl_data()
  dt[, nest := ifelse(alt == 1L, "a", "b")]
  missing_na <- "Missing columns: NA"
  expect_error(prepare_mnl_data(dt, "id", "alt", "choice", c("x1", NA)),
               missing_na, fixed = TRUE)
  expect_error(prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                                weights_col = NA_character_),
               missing_na, fixed = TRUE)
  expect_error(prepare_mnl_data(dt, 1, 2, 5, 3), "Missing columns: 1, 2, 5, 3",
               fixed = TRUE)
  expect_error(prepare_mnl_data(dt, "id", c("alt", "zz"), "choice", "x1"),
               "Missing columns: zz", fixed = TRUE)
  expect_error(prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2",
                                person_col = NA_character_),
               missing_na, fixed = TRUE)
  expect_error(prepare_nl_data(dt, "id", "alt", "choice", c("x1", NA), "nest"),
               missing_na, fixed = TRUE)
  expect_error(prepare_mnp_data(dt, "id", "alt", "choice", c("x1", NA)),
               missing_na, fixed = TRUE)
  expect_error(prepare_hmnl_data(dt, "id", "alt", "choice", c("x1", NA)),
               missing_na, fixed = TRUE)
  expect_error(prepare_hmnp_data(dt, "id", "alt", "choice", "x1",
                                 person_col = NA_character_),
               missing_na, fixed = TRUE)
})

test_that("an outside_opt_label that codes no alternative, or leaves some uncoded, is an error", {
  # Without the outside option the label names the base, the alternative
  # coded first. A label absent from the data took code 1 with no rows (a
  # phantom base): fits with constants stopped in the kernels with an error
  # about the constants' block, nested logit fits without them with one
  # about the nests, and multinomial and mixed logit fits without them
  # carried an alternative that is never offered.
  dt <- create_small_nl_data()  # alternatives 1..6 in nests {1, 2, 3}, {4, 5, 6}
  msg <- paste0("`outside_opt_label` (99) names none of the alternatives in ",
                "'alt': 1, 2, 3, 4, 5, 6. Without the outside option it names ",
                "the base alternative, coded first: name one of these, or ",
                "leave `outside_opt_label` NULL to make the first in sort order ",
                "the base.")
  expect_error(prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                                outside_opt_label = 99L), msg, fixed = TRUE)
  expect_error(prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                               "nest", outside_opt_label = 99L),
               msg, fixed = TRUE)
  expect_error(prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2",
                                outside_opt_label = 99L), msg, fixed = TRUE)
  for (ua in c(TRUE, FALSE)) {
    expect_error(run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                             use_asc = ua, outside_opt_label = 99L),
                 msg, fixed = TRUE)
    expect_error(run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                               "nest", use_asc = ua, outside_opt_label = 99L),
                 msg, fixed = TRUE)
    expect_error(run_mxlogit(dt, "id", "alt", "choice", "x1", "x2", S = 5L,
                             use_asc = ua, outside_opt_label = 99L),
                 msg, fixed = TRUE)
  }

  # Labels that coded the alternatives still do, whatever their class: the
  # check reads the codes, not the labels as text
  base <- prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                           outside_opt_label = 2L)
  expect_identical(base$alt_mapping$alt[1], 2L)
  for (lab in list(2, "2")) {
    p <- prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                          outside_opt_label = lab)
    expect_identical(p$alt_idx, base$alt_idx)
    expect_equal(p$alt_mapping, base$alt_mapping)
  }
  # Round codes that as.character() prints in scientific notation: an integer
  # label with a double column, and (as factor() matches numbers as numbers
  # from R 4.5.0) a double label with an integer column
  big <- data.table::copy(dt)[, alt := alt * 1e5]
  expect_identical(prepare_mnl_data(big, "id", "alt", "choice", "x1",
                                    outside_opt_label = 200000L)$alt_idx,
                   base$alt_idx)
  if (getRversion() >= "4.5.0") {
    big_int <- data.table::copy(dt)[, alt := alt * 100000L]
    expect_identical(prepare_mnl_data(big_int, "id", "alt", "choice", "x1",
                                      outside_opt_label = 2e5)$alt_idx,
                     base$alt_idx)
  }

  # A label of another class than the column can leave alternatives
  # without a code (the preparations returned them, and fits stopped)
  dates <- data.table::copy(dt)[, alt := as.Date("2020-01-01") + alt]
  expect_error(prepare_mnl_data(dates, "id", "alt", "choice", "x1",
                                outside_opt_label = "2020-01-03"),
               paste0("`outside_opt_label` (\"2020-01-03\", of class ",
                      "character) cannot be combined with the alternatives in ",
                      "'alt' (of class Date) without leaving some of them ",
                      "uncoded. Give the label in the class of 'alt'."),
               fixed = TRUE)
  expect_identical(prepare_mnl_data(dates, "id", "alt", "choice", "x1",
                                    outside_opt_label = as.Date("2020-01-03")
                                    )$alt_idx, base$alt_idx)
  letters_dt <- data.table::copy(dt)[, alt := letters[alt]]
  expect_error(prepare_mnl_data(letters_dt, "id", "alt", "choice", "x1",
                                outside_opt_label = factor("b")),
               "(\"b\", of class factor) cannot be combined with the",
               fixed = TRUE)
  if (requireNamespace("bit64", quietly = TRUE)) {
    i64 <- data.table::copy(dt)[, alt := bit64::as.integer64(alt)]
    expect_error(prepare_mnl_data(i64, "id", "alt", "choice", "x1",
                                  outside_opt_label = 2L),
                 "(of class integer64) without leaving some", fixed = TRUE)
  }

  # A factor's labels, a long list shortened, and situations removed for
  # missing values
  dt_f <- data.table::copy(dt)[, alt := factor(letters[alt])]
  expect_error(prepare_mnl_data(dt_f, "id", "alt", "choice", "x1",
                                outside_opt_label = "z"),
               paste0("(\"z\") names none of the alternatives in 'alt': ",
                      "\"a\", \"b\", \"c\", \"d\", \"e\", \"f\"."),
               fixed = TRUE)
  expect_identical(as.character(prepare_mnl_data(
    dt_f, "id", "alt", "choice", "x1", outside_opt_label = "c"
  )$alt_mapping$alt[1]), "c")
  dt12 <- data.table::data.table(id = rep(1:4, each = 12L),
                                 alt = rep(1:12, 4L), x1 = rnorm(48L))
  dt12[, choice := as.integer(alt == (id %% 12L) + 1L)]
  expect_error(prepare_mnl_data(dt12, "id", "alt", "choice", "x1",
                                outside_opt_label = 0L),
               "'alt': 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 and 2 more.", fixed = TRUE)
  extra <- data.table::data.table(id = 1:2, alt = 7L, x1 = NA_real_,
                                  x2 = 0, nest = 2L, choice = 0L)
  with_7 <- data.table::rbindlist(list(dt, extra), use.names = TRUE)
  expect_error(suppressWarnings(prepare_mnl_data(
    with_7, "id", "alt", "choice", "x1", outside_opt_label = 7L)),
    paste0("(7) names none of the alternatives in 'alt' that remain after ",
           "removing the 2 choice situations with missing values: 1, 2, 3, 4, ",
           "5, 6."), fixed = TRUE)

  # An NA label names no alternative either (factor() dropped it from the
  # levels, so the first alternative silently became the base); among
  # others, the absent ones are reported
  expect_error(prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                                outside_opt_label = NA),
               "(NA) names none of the alternatives in 'alt'", fixed = TRUE)
  expect_error(prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                                outside_opt_label = c(NA, 3L, 99L)),
               "(NA, 99) names none of the alternatives in 'alt'", fixed = TRUE)

  # With the outside option the label only names the outside option's row of
  # alt_mapping, which the data need not contain
  p <- prepare_mnl_data(dt, "id", "alt", "choice", "x1",
                        include_outside_option = TRUE, outside_opt_label = 99L)
  expect_identical(p$alt_mapping$alt[1], 99L)
})
