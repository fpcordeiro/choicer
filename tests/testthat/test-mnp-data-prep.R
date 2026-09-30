# Tests for prepare_mnp_data(): differencing, ASC block, base alternative,
# choice coding, and the balanced-choice-set requirement.

test_that("prepare_mnp_data differences covariates against the base", {
  dt <- data.table(
    id = c(1, 1, 1, 2, 2, 2),
    alt = c("a", "b", "c", "a", "b", "c"),
    x1 = c(1, 2, 4, 0, 1, -1),
    choice = c(0L, 1L, 0L, 1L, 0L, 0L)
  )
  d <- prepare_mnp_data(dt, "id", "alt", "choice", "x1")

  expect_s3_class(d, "choicer_data_mnp")
  expect_equal(d$J, 3L)
  expect_equal(d$p, 2L)
  expect_equal(d$N, 2L)
  expect_equal(d$base_alt, "a")          # first alternative in sort order

  # Rows are (b - a), (c - a) per id
  expect_equal(unname(d$X[, "x1"]), c(1, 3, 1, -1))

  # ASC block: p x p identity tiled N times
  expect_equal(unname(d$X[, "ASC_b"]), c(1, 0, 1, 0))
  expect_equal(unname(d$X[, "ASC_c"]), c(0, 1, 0, 1))

  # y coding: id 1 chose "b" -> 1; id 2 chose base "a" -> 0
  expect_equal(d$y, c(1L, 0L))

  expect_equal(d$param_map$beta, 1L)
  expect_equal(d$param_map$asc, c(2L, 3L))
})

test_that("prepare_mnp_data honours base_alt", {
  dt <- data.table(
    id = c(1, 1, 1, 2, 2, 2),
    alt = c("a", "b", "c", "a", "b", "c"),
    x1 = c(1, 2, 4, 0, 1, -1),
    choice = c(0L, 1L, 0L, 1L, 0L, 0L)
  )
  d <- prepare_mnp_data(dt, "id", "alt", "choice", "x1", base_alt = "c")

  expect_equal(d$base_alt, "c")
  # Levels are (c, a, b): rows are (a - c), (b - c) per id
  expect_equal(unname(d$X[, "x1"]), c(-3, -2, 1, 2))
  expect_equal(colnames(d$X), c("x1", "ASC_a", "ASC_b"))
  # id 1 chose "b" (3rd level) -> y = 2; id 2 chose "a" (2nd level) -> y = 1
  expect_equal(d$y, c(2L, 1L))

  expect_error(
    prepare_mnp_data(dt, "id", "alt", "choice", "x1", base_alt = "zzz"),
    "not one of the alternatives"
  )
})

test_that("prepare_mnp_data without ASCs has covariate columns only", {
  dt <- create_small_mnl_data()
  d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"),
                        use_asc = FALSE)
  expect_equal(colnames(d$X), c("x1", "x2"))
  expect_null(d$param_map$asc)
})

test_that("prepare_mnp_data errors on unbalanced choice sets", {
  dt <- create_small_mnl_data()
  drop_row <- which(dt$choice == 0L)[1]    # keep every id's chosen row intact
  expect_error(
    prepare_mnp_data(dt[-drop_row], "id", "alt", "choice", c("x1", "x2")),
    "balanced choice sets"
  )

  # Duplicated (id, alt) pair with the same row count also errors
  dt2 <- copy(dt)
  dt2[2, alt := 1L]
  expect_error(
    prepare_mnp_data(dt2, "id", "alt", "choice", c("x1", "x2")),
    "balanced choice sets"
  )
})

test_that("alternative-invariant covariates are dropped against the ASCs", {
  dt <- create_small_mnl_data()
  dt[, z := rnorm(1), by = id]   # constant within each choice situation
  expect_message(
    d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "z")),
    "collinearity"
  )
  expect_false("z" %in% colnames(d$X))
  expect_equal(d$dropped_cols, "z")
  expect_equal(d$param_map$beta, 1L)
})

test_that("prepare_mnp_data enforces one choice per situation", {
  dt <- create_small_mnl_data()
  dt[1:3, choice := 0L]
  expect_error(
    prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2")),
    "exactly one chosen"
  )
})

# --- the prep copies only the columns it uses ---

test_that("prepare_mnp_data leaves the caller's data unchanged", {
  dt <- create_small_mnl_data()
  set.seed(1)
  dt <- dt[sample(.N)]                 # shuffled, so the prep reorders rows
  data.table::setkeyv(dt, "x1")        # a key the prep does not sort by
  data.table::setindexv(dt, "alt")
  df <- as.data.frame(dt)
  dt0 <- data.table::copy(dt)
  df0 <- data.table::copy(df)
  prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  prepare_mnp_data(df, "id", "alt", "choice", c("x1", "x2"))
  expect_identical(dt, dt0)            # values, row order, key and index
  expect_identical(df, df0)
})

test_that("columns prepare_mnp_data does not use do not affect it", {
  df <- as.data.frame(create_small_mnl_data())
  extra <- df
  extra$note <- c(NA, "a")                     # NAs in an unused column
  extra$lst <- as.list(seq_len(nrow(extra)))   # a list column
  extra$x_unused <- NA_real_
  expect_silent(p <- prepare_mnp_data(extra, "id", "alt", "choice", c("x1", "x2")))
  expect_identical(p, prepare_mnp_data(df, "id", "alt", "choice", c("x1", "x2")))
})

test_that("a missing value in a used column drops its choice situation", {
  dt <- create_small_mnl_data()
  dt[, x3 := 1]
  dt[id == 3 & alt == 2, x2 := NA]
  dt[id == 7 & alt == 1, x3 := NA]     # x3 is not used
  expect_warning(d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2")),
                 "Removed 1 choice situations containing missing values")
  expect_equal(d$N, 19L)
  expect_identical(d, prepare_mnp_data(dt[id != 3], "id", "alt", "choice",
                                       c("x1", "x2")))
})
