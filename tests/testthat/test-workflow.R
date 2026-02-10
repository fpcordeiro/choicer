# Tests for dual convenience/advanced workflow in run_*logit() functions
# Verifies both pathways produce valid choicer_fit objects

# =============================================================================
# MNL Workflow Tests
# =============================================================================

test_that("run_mnlogit convenience workflow produces valid object", {
  dt <- create_small_mnl_data()

  fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 10L)
  )

  expect_s3_class(fit, "choicer_mnl")
  expect_s3_class(fit, "choicer_fit")
  expect_true(is.finite(fit$loglik))
  expect_true(length(fit$coefficients) > 0)
  expect_true(!is.null(fit$vcov))
  expect_true(!is.null(fit$data))
  expect_true(!is.null(fit$data_spec))
})

test_that("run_mnlogit advanced workflow produces valid object", {
  dt <- create_small_mnl_data()
  input_data <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  fit <- run_mnlogit(
    input_data = input_data,
    control = list(maxeval = 10L)
  )

  expect_s3_class(fit, "choicer_mnl")
  expect_true(is.finite(fit$loglik))
  expect_true(length(fit$coefficients) > 0)
})

test_that("run_mnlogit errors when both data and input_data provided", {
  dt <- create_small_mnl_data()
  input_data <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))

  expect_error(
    run_mnlogit(
      data = dt,
      id_col = "id",
      alt_col = "alt",
      choice_col = "choice",
      covariate_cols = c("x1", "x2"),
      input_data = input_data
    ),
    "not both"
  )
})

test_that("run_mnlogit errors when neither data nor input_data provided", {
  expect_error(
    run_mnlogit(),
    "Supply either"
  )
})

test_that("run_mnlogit convenience workflow requires column names", {
  dt <- create_small_mnl_data()

  expect_error(
    run_mnlogit(data = dt, id_col = "id"),
    "requires"
  )
})

test_that("run_mnlogit keep_data = FALSE omits data", {
  dt <- create_small_mnl_data()

  fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 10L),
    keep_data = FALSE
  )

  expect_null(fit$data)
})

# =============================================================================
# MXL Workflow Tests
# =============================================================================

test_that("run_mxlogit convenience workflow produces valid object", {
  dt <- create_small_mxl_data()

  # suppressWarnings: low maxeval can produce NaN in Hessian inversion
  fit <- suppressWarnings(run_mxlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = "x1",
    random_var_cols = c("w1", "w2"),
    S = 10L,
    control = list(maxeval = 10L)
  ))

  expect_s3_class(fit, "choicer_mxl")
  expect_s3_class(fit, "choicer_fit")
  expect_true(is.finite(fit$loglik))
  expect_true(length(fit$coefficients) > 0)
  expect_true(!is.null(fit$draws_info))
  expect_true(!is.null(fit$data))
  expect_true(!is.null(fit$data_spec))
  expect_equal(fit$draws_info$S, 10L)
})

test_that("run_mxlogit advanced workflow produces valid object", {
  dt <- create_small_mxl_data()
  input_data <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2"),
    rc_correlation = FALSE
  )
  eta_draws <- get_halton_normals(10L, input_data$N, ncol(input_data$W))

  fit <- suppressWarnings(run_mxlogit(
    input_data = input_data,
    eta_draws = eta_draws,
    control = list(maxeval = 10L)
  ))

  expect_s3_class(fit, "choicer_mxl")
  expect_true(is.finite(fit$loglik))
})

test_that("run_mxlogit errors when both data and input_data provided", {
  dt <- create_small_mxl_data()
  input_data <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2")
  )

  expect_error(
    run_mxlogit(
      data = dt,
      id_col = "id",
      alt_col = "alt",
      choice_col = "choice",
      covariate_cols = "x1",
      random_var_cols = c("w1", "w2"),
      input_data = input_data,
      eta_draws = array(0, c(2, 10, 30))
    ),
    "not both"
  )
})

test_that("run_mxlogit advanced workflow requires eta_draws", {
  dt <- create_small_mxl_data()
  input_data <- prepare_mxl_data(
    dt, "id", "alt", "choice", "x1", c("w1", "w2")
  )

  expect_error(
    run_mxlogit(input_data = input_data),
    "eta_draws"
  )
})

# =============================================================================
# NL Workflow Tests
# =============================================================================

test_that("run_nestlogit convenience workflow produces valid object", {
  dt <- create_small_nl_data()

  fit <- suppressWarnings(run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    control = list(maxeval = 10L)
  ))

  expect_s3_class(fit, "choicer_nl")
  expect_s3_class(fit, "choicer_fit")
  expect_true(is.finite(fit$loglik))
  expect_true(length(fit$coefficients) > 0)
  expect_true(!is.null(fit$lambda))
  expect_true(!is.null(fit$data))
  expect_true(!is.null(fit$data_spec))
  expect_equal(fit$data_spec$nest_col, "nest")
})

test_that("run_nestlogit advanced workflow produces valid object", {
  dt <- create_small_nl_data()
  nl_data <- prepare_nl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    nest_col = "nest"
  )

  fit <- suppressWarnings(run_nestlogit(
    input_data = nl_data,
    control = list(maxeval = 10L)
  ))

  expect_s3_class(fit, "choicer_nl")
  expect_true(is.finite(fit$loglik))
})

test_that("run_nestlogit errors when both data and input_data provided", {
  dt <- create_small_nl_data()
  nl_data <- prepare_nl_data(
    dt, "id", "alt", "choice", c("x1", "x2"),
    nest_col = "nest"
  )

  expect_error(
    run_nestlogit(
      data = dt,
      id_col = "id",
      alt_col = "alt",
      choice_col = "choice",
      covariate_cols = c("x1", "x2"),
      nest_col = "nest",
      input_data = nl_data
    ),
    "not both"
  )
})

test_that("run_nestlogit convenience workflow requires nest_col", {
  dt <- create_small_nl_data()

  expect_error(
    run_nestlogit(
      data = dt,
      id_col = "id",
      alt_col = "alt",
      choice_col = "choice",
      covariate_cols = c("x1", "x2")
    ),
    "nest_col"
  )
})

# =============================================================================
# S3 method tests on fitted objects
# =============================================================================

test_that("coef, logLik, nobs, vcov work on MNL fit", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 10L)
  )

  expect_true(length(coef(fit)) > 0)
  expect_s3_class(logLik(fit), "logLik")
  expect_true(nobs(fit) > 0)
  expect_true(is.matrix(vcov(fit)))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("summary and print work on MNL fit", {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 10L)
  )

  s <- summary(fit)
  expect_s3_class(s, "summary.choicer_mnl")
  expect_true("coefficients" %in% names(s))

  # print should not error
  expect_output(print(fit))
  expect_output(print(s))
})

test_that("summary and print work on MXL fit", {
  dt <- create_small_mxl_data()
  fit <- suppressWarnings(run_mxlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = "x1",
    random_var_cols = c("w1", "w2"),
    S = 10L,
    control = list(maxeval = 10L)
  ))

  s <- suppressWarnings(summary(fit))
  expect_s3_class(s, "summary.choicer_mxl")
  expect_output(print(s))
})

test_that("summary and print work on NL fit", {
  dt <- create_small_nl_data()
  fit <- suppressWarnings(run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    control = list(maxeval = 10L)
  ))

  s <- summary(fit)
  expect_s3_class(s, "summary.choicer_nl")
  expect_output(print(s))
})
