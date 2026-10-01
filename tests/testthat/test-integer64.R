# integer64 (bit64) inputs in every design builder
#
# integer64 keeps 64-bit integers in the storage of a double vector, and
# is.numeric() accepts it, so code that reads the storage directly
# (as.matrix(), the C++ kernels) used to see the raw bit patterns: 1 as
# 4.9e-324, 2^40 as 5.4e-312, small negative values as NaN; and base R's
# is.na() misreads them unless bit64 is loaded. Every builder now reads
# integer64 covariates and weights as their values, as as.double() converts
# them, so an integer64 input must give the same prepared object, fit and
# prediction as the same values stored as double.

# Long-format choice data. `big` holds integer values that an integer64
# column typically carries (negative ones, whose bit patterns are NaNs, and
# ones beyond 2^31 - 1); `cnt` small counts; `z` an alternative-level
# covariate; `w` a weight constant within each choice situation. With
# `as_int64`, those four columns are integer64, otherwise double.
int64_choice_data <- function(as_int64, N = 40L, J = 3L, seed = 11L) {
  set.seed(seed)
  dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  dt[, x1 := rnorm(.N)]
  dt[, big := sample(c(-7, -1, 0, 3, 2^31, 2^40 + 1), .N, replace = TRUE)]
  dt[, cnt := as.numeric(sample(0:9, .N, replace = TRUE))]
  dt[, z := 3 * alt - 1]
  dt[, w := rep(as.numeric(sample(1:4, N, replace = TRUE)), each = J)]
  dt[, nest := ifelse(alt == 1L, "a", "b")]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1L))), by = id]
  i64_cols <- c("big", "cnt", "z", "w")
  if (as_int64) {
    dt[, (i64_cols) := lapply(.SD, bit64::as.integer64), .SDcols = i64_cols]
  }
  dt[]
}

# An integer64 matrix with the (integer) values and shape of `m`
int64_matrix <- function(m) {
  out <- bit64::as.integer64(as.vector(m))
  dim(out) <- dim(m)
  dimnames(out) <- dimnames(m)
  out
}

# The messages of the warnings `expr` gives, which are muffled
warning_messages <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  msgs
}

# The warning for integer64 values from 2^53 up, and the error without bit64
rounded_msg <- function(what) {
  paste0(what, " has integer64 values of magnitude 2^53 or more; they ",
         "were rounded to the nearest double.")
}
no_bit64_msg <- function(what) {
  paste0(what, " is of class integer64; install the bit64 package to read ",
         "it.")
}

test_that("integer64 covariates enter the MNL, NL and MXL designs as numbers", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)

  mnl <- prepare_mnl_data(d64, "id", "alt", "choice", c("x1", "big", "cnt"))
  expect_identical(
    mnl, prepare_mnl_data(dbl, "id", "alt", "choice", c("x1", "big", "cnt"))
  )
  expect_setequal(mnl$X[, "big"], c(-7, -1, 0, 3, 2^31, 2^40 + 1))

  expect_identical(
    prepare_nl_data(d64, "id", "alt", "choice", c("x1", "big"),
                    nest_col = "nest"),
    prepare_nl_data(dbl, "id", "alt", "choice", c("x1", "big"),
                    nest_col = "nest")
  )

  # integer64 in the fixed (X) and in the random-coefficient (W) design
  mxl <- prepare_mxl_data(d64, "id", "alt", "choice", c("x1", "big"),
                          random_var_cols = "cnt")
  expect_identical(
    mxl, prepare_mxl_data(dbl, "id", "alt", "choice", c("x1", "big"),
                          random_var_cols = "cnt")
  )
  expect_true(is.double(mxl$W))
  expect_identical(
    prepare_mxl_data(d64, "id", "alt", "choice", "x1", random_var_cols = "big"),
    prepare_mxl_data(dbl, "id", "alt", "choice", "x1", random_var_cols = "big")
  )

  # the caller's data is left as it was
  expect_s3_class(d64$big, "integer64")
})

test_that("a missing integer64 value drops its choice situation", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)
  set(d64, 5L, "big", bit64::NA_integer64_)
  set(dbl, 5L, "big", NA_real_)

  expect_warning(
    mnl <- prepare_mnl_data(d64, "id", "alt", "choice", c("x1", "big")),
    "Removed 1 choice situations containing missing values"
  )
  expect_warning(
    ref <- prepare_mnl_data(dbl, "id", "alt", "choice", c("x1", "big")),
    "Removed 1 choice situations containing missing values"
  )
  expect_identical(mnl, ref)
  expect_identical(mnl$N, 39L)
})

test_that("integer64 covariates enter the MNP differenced design", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)

  mnp <- prepare_mnp_data(d64, "id", "alt", "choice", c("x1", "big"))
  expect_identical(mnp,
                   prepare_mnp_data(dbl, "id", "alt", "choice", c("x1", "big")))
  expect_true(all(is.finite(mnp$X)))
})

test_that("integer64 structural, alternative-level and control-function columns enter the HB designs", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)

  hmnl <- prepare_hmnl_data(d64, "id", "alt", "choice", c("x1", "big"),
                            alt_covariate_cols = "z", cf_residual_col = "cnt")
  expect_identical(
    hmnl, prepare_hmnl_data(dbl, "id", "alt", "choice", c("x1", "big"),
                            alt_covariate_cols = "z", cf_residual_col = "cnt")
  )
  expect_identical(unname(hmnl$Z[, "z"]), c(2, 5, 8))

  expect_identical(
    prepare_hmnp_data(d64, "id", "alt", "choice", c("x1", "big"),
                      alt_covariate_cols = "z", cf_residual_col = "cnt"),
    prepare_hmnp_data(dbl, "id", "alt", "choice", c("x1", "big"),
                      alt_covariate_cols = "z", cf_residual_col = "cnt")
  )
})

test_that("integer64 weights are read as numbers, as a column or a vector", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)

  mnl <- prepare_mnl_data(d64, "id", "alt", "choice", "x1", weights_col = "w")
  expect_identical(
    mnl, prepare_mnl_data(dbl, "id", "alt", "choice", "x1", weights_col = "w")
  )
  expect_identical(mnl$weights, dbl$w[dbl$alt == 1L])
  expect_identical(
    prepare_mxl_data(d64, "id", "alt", "choice", "x1", random_var_cols = "cnt",
                     weights_col = "w"),
    prepare_mxl_data(dbl, "id", "alt", "choice", "x1", random_var_cols = "cnt",
                     weights_col = "w")
  )

  # positional weights, one per choice situation in ascending-id order
  w <- dbl$w[dbl$alt == 1L]
  expect_identical(
    prepare_mnl_data(dbl, "id", "alt", "choice", "x1",
                     weights = bit64::as.integer64(w)),
    prepare_mnl_data(dbl, "id", "alt", "choice", "x1", weights = w)
  )
  expect_identical(
    prepare_mxl_data(dbl, "id", "alt", "choice", "x1", random_var_cols = "cnt",
                     weights = bit64::as.integer64(w)),
    prepare_mxl_data(dbl, "id", "alt", "choice", "x1", random_var_cols = "cnt",
                     weights = w)
  )
})

test_that("fits and counterfactuals on integer64 data match those on doubles", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)
  fit_mnl <- function(data) {
    run_mnlogit(data, "id", "alt", "choice", c("x1", "cnt"),
                control = list(maxeval = 50L))
  }
  fit <- fit_mnl(dbl)
  fit64 <- fit_mnl(d64)
  # identical kernel inputs; the fits themselves agree only to the round-off
  # of the kernels' multi-threaded reductions, which vary from run to run
  expect_identical(fit64$data, fit$data)
  expect_equal(coef(fit64), coef(fit), tolerance = 1e-10)
  expect_equal(vcov(fit64), vcov(fit), tolerance = 1e-10)

  expect_identical(predict(fit, newdata = d64), predict(fit, newdata = dbl))
  w <- dbl$w[dbl$alt == 1L]
  expect_equal(  # shares are summed across threads, like the fits above
    predict(fit, type = "shares", newdata = d64,
            weights = bit64::as.integer64(w)),
    predict(fit, type = "shares", newdata = dbl, weights = w),
    tolerance = 1e-12
  )
  expect_identical(logsum(fit, newdata = d64), logsum(fit, newdata = dbl))
  expect_identical(consumer_surplus(fit, "cnt", newdata = d64),
                   consumer_surplus(fit, "cnt", newdata = dbl))

  # MXL: newdata fills both X and W
  mxl <- run_mxlogit(dbl, "id", "alt", "choice", "x1",
                     random_var_cols = "cnt", S = 10L,
                     control = list(maxeval = 30L))
  nd <- prepare_newdata(mxl, d64)
  expect_identical(nd, prepare_newdata(mxl, dbl))
  expect_true(is.double(nd$W))
})

test_that("the list form of newdata takes integer64 X and W matrices", {
  skip_if_not_installed("bit64")
  dbl <- int64_choice_data(FALSE)
  fit <- run_mnlogit(dbl, "id", "alt", "choice", c("x1", "cnt"),
                     control = list(maxeval = 50L))
  d <- fit$data
  X <- round(10 * d$X)
  expect_identical(
    validate_newdata_list(fit, list(X = int64_matrix(X), alt_idx = d$alt_idx,
                                    M = d$M)),
    validate_newdata_list(fit, list(X = X, alt_idx = d$alt_idx, M = d$M))
  )
  expect_identical(
    predict(fit, newdata = list(X = int64_matrix(X), alt_idx = d$alt_idx,
                                M = d$M)),
    predict(fit, newdata = list(X = X, alt_idx = d$alt_idx, M = d$M))
  )

  mxl <- run_mxlogit(dbl, "id", "alt", "choice", "x1",
                     random_var_cols = "cnt", S = 10L,
                     control = list(maxeval = 30L))
  dm <- mxl$data
  X <- round(10 * dm$X)
  expect_identical(
    validate_newdata_list(mxl, list(X = int64_matrix(X), W = int64_matrix(dm$W),
                                    alt_idx = dm$alt_idx, M = dm$M)),
    validate_newdata_list(mxl, list(X = X, W = dm$W, alt_idx = dm$alt_idx,
                                    M = dm$M))
  )
})

test_that("hierarchical Bayes prediction converts integer64 newdata", {
  skip_if_not_installed("bit64")
  sim <- simulate_hmnl_data(N = 30, T = 2, J = 4, seed = 42)
  dbl <- sim$data
  dbl[, x1 := round(3 * x1)]
  dbl[, z1 := round(4 * z1)]
  d <- prepare_hmnl_data(dbl, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, mcmc = list(R = 60, burn = 20, seed = 1))
  )

  # a new alternative (posterior-predictive delta) with its own z1
  new_rows <- dbl[alt == 1L][, `:=`(alt = 9L, z1 = 7)]
  nd <- rbind(dbl, new_rows)
  nd64 <- copy(nd)[, c("x1", "z1") := lapply(.SD, bit64::as.integer64),
                   .SDcols = c("x1", "z1")]
  rd <- .hb_resolve_newdata(fit, nd64)
  expect_identical(rd, .hb_resolve_newdata(fit, nd))
  expect_identical(unname(rd$z_new[, "z1"]), 7)

  local_mocked_bindings(.bit64_available = function() FALSE)
  expect_error(.hb_resolve_newdata(fit, nd64), no_bit64_msg("Column 'x1'"),
               fixed = TRUE)
})

test_that("integer64 values from 2^53 up are rounded, with a warning naming the column", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  # 2^54 + 1: the doubles around it are 4 apart, so the nearest is 2^54
  set(d64, 1L, "big", bit64::as.integer64("18014398509481985"))
  expect_warning(
    mnl <- prepare_mnl_data(d64, "id", "alt", "choice", c("x1", "big")),
    "Column 'big' has integer64 values of magnitude 2^53 or more",
    fixed = TRUE
  )
  expect_identical(unname(mnl$X[1L, "big"]), 2^54)

  # the probit differences the rounded values; id 1's base row is row 1
  expect_warning(
    mnp <- prepare_mnp_data(d64, "id", "alt", "choice", c("x1", "big")),
    "Column 'big' has integer64 values of magnitude 2^53 or more",
    fixed = TRUE
  )
  big <- int64_choice_data(FALSE)$big
  expect_identical(unname(mnp$X[1:2, "big"]), big[2:3] - 2^54)

  # values below 2^53 never warn, wherever they sit in the column
  expect_no_warning(
    prepare_mnl_data(int64_choice_data(TRUE), "id", "alt", "choice",
                     c("x1", "big"))
  )
})

test_that("prep_gather_design() reads integer64 storage as bit64's as.double()", {
  skip_if_not_installed("bit64")
  x <- bit64::as.integer64(c("-9007199254740993", "-7", "-1", "0", "1",
                             "2147483648", "1099511627777",
                             "9007199254740993", "9223372036854775807", NA))
  ref <- suppressWarnings(as.double(x))
  rows <- c(10L, 3L, 1L, 8L, 9L, 2L, 4L, 5L, 6L, 7L)
  g <- prep_gather_design(list(x), rows)
  expect_identical(g[, 1], ref[rows])
  expect_identical(attr(g, "int64_big"), TRUE)

  # differenced (MNP) and alongside other columns: only integer64 columns
  # with values from 2^53 up are flagged
  small <- x[2:7]
  d <- prep_gather_design(list(1:6, small, x[1:6]), 1:6, base = rep(4L, 6))
  expect_identical(d[, 2], as.double(small) - as.double(small[4]))
  expect_identical(d[, 3], ref[1:6] - ref[4])
  expect_identical(attr(d, "int64_big"), c(FALSE, FALSE, TRUE))
  expect_null(attributes(prep_gather_design(list(small), 1:6))[["int64_big"]])

  # the multi-threaded branch (more than 100000 rows) agrees, NA included
  n <- 200001L
  set.seed(5)
  v <- bit64::as.integer64(sample(-1e6:1e6, n, replace = TRUE))
  v[c(7L, n)] <- NA
  idx <- sample.int(n)
  expect_identical(prep_gather_design(list(v), idx)[, 1],
                   as.double(v)[idx])
  expect_identical(
    prep_gather_design(list(v), idx, base = rev(idx))[, 1],
    as.double(v)[idx] - as.double(v)[rev(idx)]
  )

  # .gather_matrix() turns the flag into one warning naming the column
  df <- data.frame(a = 1:10)
  df$big <- x
  expect_warning(
    m <- .gather_matrix(df, c("a", "big", "big"), rows, "X"),
    "Column 'big' has integer64 values of magnitude 2^53 or more",
    fixed = TRUE
  )
  expect_null(attr(m, "int64_big"))
  expect_identical(dimnames(m), list(NULL, c("a", "big", "big")))
})

test_that("each object converted warns once about values from 2^53 up", {
  skip_if_not_installed("bit64")
  big <- bit64::as.integer64("18014398509481985")  # 2^54 + 1; nearest 2^54
  d64 <- int64_choice_data(TRUE)
  # set(), not d64[..., := big], where `big` would be the column of that name
  set(d64, 1L, "big", big)
  set(d64, which(d64$alt == 1L), "z", big)  # constant within alternative 1
  set(d64, which(d64$id == 1L), "w", big)   # one weight per choice situation

  # a column in both X and W, in both X and Z, or twice in Z: one warning
  expect_identical(
    warning_messages(prepare_mxl_data(d64, "id", "alt", "choice",
                                      c("x1", "big"), random_var_cols = "big")),
    rounded_msg("Column 'big'")
  )
  expect_identical(
    warning_messages(prepare_hmnl_data(d64, "id", "alt", "choice",
                                       c("x1", "z"), alt_covariate_cols = "z")),
    rounded_msg("Column 'z'")
  )
  expect_identical(
    warning_messages(hb <- suppressMessages(
      prepare_hmnl_data(d64, "id", "alt", "choice", "x1",
                        alt_covariate_cols = c("z", "z"))
    )),
    rounded_msg("Column 'z'")
  )
  expect_identical(unname(hb$Z[1L, "z"]), 2^54)

  # the paths converted in R: weights, newdata and its list form
  expect_identical(
    warning_messages(mnl <- prepare_mnl_data(d64, "id", "alt", "choice", "x1",
                                             weights_col = "w")),
    rounded_msg("Column 'w'")
  )
  expect_identical(mnl$weights[1L], 2^54)
  expect_identical(
    warning_messages(prepare_mnl_data(
      int64_choice_data(FALSE), "id", "alt", "choice", "x1",
      weights = c(big, bit64::as.integer64(rep(1, 39)))
    )),
    rounded_msg("`weights`")
  )
  expect_identical(warning_messages(.validate_pred_weights(big, 1L)),
                   rounded_msg("'weights'"))

  fit <- run_mnlogit(int64_choice_data(FALSE), "id", "alt", "choice",
                     c("x1", "cnt"), control = list(maxeval = 50L))
  nd <- int64_choice_data(TRUE)
  set(nd, 1L, "cnt", big)
  expect_identical(warning_messages(prepare_newdata(fit, nd)),
                   rounded_msg("Column 'cnt'"))
  X <- bit64::as.integer64(as.vector(round(10 * fit$data$X)))
  X[1L] <- big
  dim(X) <- dim(fit$data$X)
  expect_identical(
    warning_messages(validate_newdata_list(
      fit, list(X = X, alt_idx = fit$data$alt_idx, M = fit$data$M)
    )),
    rounded_msg("newdata$X")
  )
})

test_that("bit64 is loaded before the missing-value scan reads the columns", {
  skip_if_not_installed("bit64")
  loaded <- FALSE
  calls <- 0L
  loaded_at_scan <- logical(0)
  rows_with_na <- .rows_with_na
  local_mocked_bindings(
    .bit64_available = function() {
      loaded <<- TRUE
      calls <<- calls + 1L
      TRUE
    },
    .rows_with_na = function(x, cols = seq_along(x)) {
      loaded_at_scan <<- c(loaded_at_scan, loaded)
      rows_with_na(x, cols)
    }
  )
  d64 <- int64_choice_data(TRUE)
  preps <- list(
    function() prepare_mnl_data(d64, "id", "alt", "choice", c("x1", "big")),
    function() prepare_mxl_data(d64, "id", "alt", "choice", "x1",
                                random_var_cols = "big"),
    function() prepare_mnp_data(d64, "id", "alt", "choice", c("x1", "big")),
    function() prepare_hmnl_data(d64, "id", "alt", "choice", c("x1", "big")),
    # only an index column is integer64; the scan reads it too
    function() {
      dbl <- int64_choice_data(FALSE)
      dbl[, id64 := bit64::as.integer64(id)]
      prepare_mnl_data(dbl, "id64", "alt", "choice", "x1")
    }
  )
  for (prep in preps) {
    loaded <- FALSE
    prep()
  }
  expect_identical(loaded_at_scan, rep(TRUE, length(preps)))

  # and bit64 is left alone when nothing is integer64
  calls <- 0L
  prepare_mnl_data(int64_choice_data(FALSE), "id", "alt", "choice", "x1")
  prepare_hmnl_data(int64_choice_data(FALSE), "id", "alt", "choice", "x1")
  expect_identical(calls, 0L)
})

test_that("without bit64, every path that meets an integer64 column says so", {
  skip_if_not_installed("bit64")
  d64 <- int64_choice_data(TRUE)
  dbl <- int64_choice_data(FALSE)
  d_nest <- copy(dbl)[, nest := bit64::as.integer64(ifelse(alt == 1L, 1, 2))]
  fit <- run_mnlogit(dbl, "id", "alt", "choice", c("x1", "cnt"),
                     control = list(maxeval = 50L))
  X64 <- int64_matrix(round(10 * fit$data$X))
  local_mocked_bindings(.bit64_available = function() FALSE)

  for (prep in list(prepare_mnl_data, prepare_mnp_data, prepare_hmnl_data)) {
    expect_error(prep(d64, "id", "alt", "choice", c("x1", "big")),
                 no_bit64_msg("Column 'big'"), fixed = TRUE)
  }
  expect_error(
    prepare_mxl_data(d64, "id", "alt", "choice", "x1", random_var_cols = "big"),
    no_bit64_msg("Column 'big'"), fixed = TRUE
  )
  expect_error(
    prepare_nl_data(d_nest, "id", "alt", "choice", "x1", nest_col = "nest"),
    no_bit64_msg("Column 'nest'"), fixed = TRUE
  )
  expect_error(
    prepare_mnl_data(dbl, "id", "alt", "choice", "x1",
                     weights = bit64::as.integer64(rep(1, 40))),
    no_bit64_msg("`weights`"), fixed = TRUE
  )
  expect_error(prepare_newdata(fit, d64), no_bit64_msg("Column 'cnt'"),
               fixed = TRUE)
  expect_error(
    validate_newdata_list(fit, list(X = X64, alt_idx = fit$data$alt_idx,
                                    M = fit$data$M)),
    no_bit64_msg("newdata$X"), fixed = TRUE
  )
})
