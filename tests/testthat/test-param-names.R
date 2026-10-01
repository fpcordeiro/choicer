# Parameter names must be unique: coefficients, vcov(), summary(), wtp()
# and run_mxlogit()'s named bounds look parameters up by name, so a
# covariate may not take a name the model generates (ASC_<label>,
# Mu_<variable>, L_<i><j>, Lambda_<k>) or a label summary() prints in its
# place (Sigma_<i><j>, exp(Mu_<variable>)).

make_clash_data <- function(seed = 21) {
  set.seed(seed)
  dt <- data.table(id = rep(1:40, each = 3), alt = rep(c("a", "b", "c"), 40))
  dt[, `:=`(x1 = rnorm(.N), w = rnorm(.N), w2 = rnorm(.N),
            pw = runif(.N, 0.5, 2), ASC_b = rnorm(.N), Mu_w = rnorm(.N),
            Mu_pw = rnorm(.N), `exp(Mu_pw)` = rnorm(.N), L_11 = rnorm(.N),
            L_21 = rnorm(.N), Sigma_11 = rnorm(.N), Sigma_21 = rnorm(.N),
            Lambda_1 = rnorm(.N))]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt[, nest := ifelse(alt == "a", "A", "B")]   # one lambda: Lambda_1
  dt[]
}

test_that("run_mnlogit() rejects a covariate named like an ASC", {
  dt <- make_clash_data()
  expect_error(
    run_mnlogit(dt, "id", "alt", "choice", c("x1", "ASC_b")),
    "Parameter names must be unique; repeated: 'ASC_b'. Covariates may not",
    fixed = TRUE
  )
  # Without ASCs the name is free
  fit <- suppressMessages(run_mnlogit(dt, "id", "alt", "choice",
                                      c("x1", "ASC_b"), use_asc = FALSE))
  expect_identical(names(coef(fit)), c("x1", "ASC_b"))

  # The advanced workflow is checked too
  d <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "w"))
  colnames(d$X) <- c("x", "x")
  expect_error(run_mnlogit(input_data = d), "repeated: 'x'.", fixed = TRUE)
})

test_that("run_nestlogit() rejects covariates named like a lambda or an ASC", {
  dt <- make_clash_data()
  nl <- function(...) {
    suppressMessages(run_nestlogit(dt, "id", "alt", "choice",
                                   nest_col = "nest", ...))
  }
  expect_error(nl(covariate_cols = c("x1", "Lambda_1")),
               "repeated: 'Lambda_1'.", fixed = TRUE)
  expect_error(nl(covariate_cols = c("ASC_b", "x1")),
               "repeated: 'ASC_b'.", fixed = TRUE)
  # The check runs before the optimizer, on supplied names too
  expect_error(nl(covariate_cols = "x1",
                  param_names = c("x1", "x1", "ASC_b", "ASC_c")),
               "repeated: 'x1'.", fixed = TRUE)
  # Unique supplied names free the generated ones
  fit <- nl(covariate_cols = c("x1", "Lambda_1"),
            param_names = c("x1", "Lambda_1", "lambda_B", "ASC_b", "ASC_c"))
  expect_identical(names(coef(fit)),
                   c("x1", "Lambda_1", "lambda_B", "ASC_b", "ASC_c"))
})

test_that("run_mxlogit() rejects covariates named like its parameters", {
  dt <- make_clash_data()
  mxl <- function(cov, rc = "w", ...) {
    suppressMessages(run_mxlogit(dt, "id", "alt", "choice", cov,
                                 random_var_cols = rc, S = 10L, ...))
  }
  expect_error(mxl(c("x1", "ASC_b")), "repeated: 'ASC_b'.", fixed = TRUE)
  expect_error(mxl(c("x1", "Mu_w"), rc_mean = TRUE), "repeated: 'Mu_w'.",
               fixed = TRUE)
  expect_error(mxl(c("x1", "L_11")), "repeated: 'L_11'.", fixed = TRUE)
  # Sigma_11 is no parameter name, but summary() prints L_11 as Sigma_11
  expect_error(mxl(c("x1", "Sigma_11")), "repeated: 'Sigma_11'.",
               fixed = TRUE)
  expect_error(mxl(c("L_11", "Sigma_11", "ASC_b")),
               "repeated: 'L_11', 'Sigma_11', 'ASC_b'.", fixed = TRUE)
  # Correlated coefficients add off-diagonal names
  expect_error(mxl(c("x1", "L_21"), rc = c("w", "w2"), rc_correlation = TRUE),
               "repeated: 'L_21'.", fixed = TRUE)
  expect_error(mxl(c("x1", "Sigma_21"), rc = c("w", "w2"),
                   rc_correlation = TRUE),
               "repeated: 'Sigma_21'.", fixed = TRUE)
  # A log-normal mean is Mu_pw in the coefficients, exp(Mu_pw) in summary()
  expect_error(mxl(c("x1", "Mu_pw"), rc = "pw", rc_dist = 1L, rc_mean = TRUE),
               "repeated: 'Mu_pw'.", fixed = TRUE)
  expect_error(mxl(c("x1", "exp(Mu_pw)"), rc = "pw", rc_dist = 1L,
                   rc_mean = TRUE),
               "repeated: 'exp(Mu_pw)'.", fixed = TRUE)

  # Names the configuration does not generate are free: Mu_w needs
  # rc_mean = TRUE, L_21 needs correlated coefficients
  fit <- mxl(c("x1", "Mu_w"))
  expect_identical(names(coef(fit)), c("x1", "Mu_w", "L_11", "ASC_b", "ASC_c"))
  expect_identical(rownames(summary(fit, gof = FALSE)$coefficients),
                   c("x1", "Mu_w", "Sigma_11", "ASC_b", "ASC_c"))
  fit <- mxl(c("x1", "L_21"), rc = c("w", "w2"))
  expect_identical(names(coef(fit)),
                   c("x1", "L_21", "L_11", "L_22", "ASC_b", "ASC_c"))
})

test_that("covariance labels follow the Cholesky parameter order", {
  expect_identical(.mxl_cov_names("L", 3L, TRUE),
                   c("L_11", "L_21", "L_22", "L_31", "L_32", "L_33"))
  expect_identical(.mxl_cov_names("Sigma", 3L, FALSE),
                   c("Sigma_11", "Sigma_22", "Sigma_33"))
  expect_identical(.mxl_cov_names("L", 0L, TRUE), character(0))
})
