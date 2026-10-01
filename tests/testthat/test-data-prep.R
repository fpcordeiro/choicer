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
  expect_identical(.gather_matrix(df, c("x1", "x1"), seq_len(nrow(df)), "X"),
                   num(c("x1", "x1")))
  cc <- character(0)
  expect_identical(.gather_matrix(df, cc, seq_len(nrow(df)), "X"),
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

test_that("design matrices past 2^32 - 1 values are refused before they are built", {
  expect_error(.check_design_size(2^31, c("a", "b"), "The design matrix X"),
               "X would have 2,147,483,648 rows and 2 columns.*more than 2\\^32 - 1")
  expect_silent(.check_design_size(2^31 - 1, c("a", "b"), "The design matrix X"))
  # Wired into .gather_matrix() ahead of any allocation (seq_len() is compact).
  df <- as.data.frame(create_small_nl_data())
  expect_error(.gather_matrix(df, c("x1", "x2", "id"), seq_len(.Machine$integer.max), "X"),
               "more than 2\\^32 - 1")
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
