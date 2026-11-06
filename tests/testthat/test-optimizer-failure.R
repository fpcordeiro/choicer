# An optimizer's failure code (nloptr's NLOPT_FAILURE, -1, and its kin): the
# code alone cannot tell the start values from the maximum, so the wrappers
# certify the point (a Newton phase of no steps where no polish was asked for)
# or warn (class choicer_optimizer_failed). The optimizer is mocked where the
# path matters, so that these do not depend on where L-BFGS's line search
# gives up; the one real failure is a first step one gradient long on badly
# scaled columns.

# The real nloptr adapter (read once, before any mock), relabelled: code -1 at
# whatever point it stops
real_nloptr <- choicer:::run_nloptr
fail_after_run <- function(theta_init, eval_f, lower, upper, control) {
  r <- real_nloptr(theta_init, eval_f, lower, upper, control)
  r$convergence <- -1L
  r$message <- "NLOPT_FAILURE: Generic failure code."
  r
}

# Code -1 at the start values, as the line search leaves them
fail_at_start <- function(theta_init, eval_f, lower, upper, control) {
  list(par = theta_init, value = eval_f(theta_init)$objective, convergence = -1L,
       message = "NLOPT_FAILURE: Generic failure code.", iterations = 12L, raw = list())
}

# Long-format choices, J = 4 and a column of large values (a first trial step
# one gradient long lands far past the maximum)
ofail_data <- function(N = 400L, scale = 1e4, seed = 11L) {
  set.seed(seed)
  J <- 4L
  x <- rnorm(N * J) * scale
  z <- rnorm(N * J)
  v <- 0.5 * x / scale - 0.5 * z + rep(c(0, 0.3, -0.2, 0.1), N)
  u <- v - log(-log(runif(N * J)))
  id <- rep(seq_len(N), each = J)
  data.table::data.table(id = id, alt = rep(seq_len(J), N), x = x, z = z,
                         nest = rep(c(1L, 1L, 2L, 2L), N),
                         choice = as.integer(u == ave(u, id, FUN = max)))
}

test_that(".optimizer_failed() reads each optimizer's failure codes", {
  failed <- choicer:::.optimizer_failed
  for (code in c(-1L, -2L, -3L, -4L, -5L)) {
    expect_true(failed(list(convergence = code), "nloptr"))
    expect_true(failed(list(convergence = code), NULL))
  }
  # successes, and the caps
  for (code in 1:6) expect_false(failed(list(convergence = code), "nloptr"))
  for (code in c(10L, 51L, 52L)) expect_true(failed(list(convergence = code), "optim"))
  for (code in c(0L, 1L)) expect_false(failed(list(convergence = code), "optim"))
  # a custom optimizer's codes follow no convention
  expect_false(failed(list(convergence = -1L), function(...) NULL))
  expect_false(failed(list(convergence = NA_integer_), "nloptr"))
})

test_that("run_nestlogit() warns where L-BFGS fails at the start, and the coordinates fix it", {
  skip_on_cran()  # the optimizer's real path
  d <- ofail_data()
  # the fit's own default, "none": a first step one gradient long
  expect_warning(
    fit <- suppressMessages(run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest")),
    class = "choicer_optimizer_failed")
  expect_identical(fit$convergence, -1L)
  expect_match(fit$message, "NLOPT_FAILURE", fixed = TRUE)
  start <- c(0, 0, 0.5, 0.5, 0, 0, 0)
  expect_equal(unname(coef(fit)), start, tolerance = 1e-12)
  # the check's record: no steps, nothing certified
  pol <- fit$optimizer$polish
  expect_identical(pol$iterations, 0L)
  expect_false(pol$converged)
  expect_identical(pol$optimizer$convergence, -1L)
  expect_false(choicer:::.is_converged(fit))
  # scale_vars = "bhhh" and a polish both reach the same maximum, without a warning
  ok <- suppressMessages(run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest",
                                       scale_vars = "bhhh", polish = "newton"))
  expect_true(ok$optimizer$polish$converged)
  expect_gt(as.numeric(logLik(ok)), as.numeric(logLik(fit)) + 1)
  pol2 <- suppressMessages(run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest",
                                         polish = "newton"))
  expect_equal(as.numeric(logLik(pol2)), as.numeric(logLik(ok)), tolerance = 1e-8)
  expect_true(pol2$optimizer$polish$converged)
  expect_identical(pol2$optimizer$polish$optimizer$convergence, -1L)
})

test_that("the warning says what happened at the start values", {
  d <- ofail_data(N = 200L, scale = 1)
  local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
  for (run in list(
    function(...) run_mnlogit(d, "id", "alt", "choice", c("x", "z"), ...),
    function(...) run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest", ...))) {
    w <- NULL
    fit <- withCallingHandlers(
      suppressMessages(run()),
      choicer_optimizer_failed = function(c) {
        w <<- conditionMessage(c)
        invokeRestart("muffleWarning")
      })
    expect_match(w, "NLOPT_FAILURE: Generic failure code.", fixed = TRUE)
    expect_match(w, "returned the start values", fixed = TRUE)
    expect_match(w, "polish = \"newton\"", fixed = TRUE)
    expect_identical(fit$convergence, -1L)
    # the variance is the one at the point returned
    expect_false(is.null(fit$vcov))
    expect_false(choicer:::.is_converged(fit))
  }
})

test_that("a failure code at the maximum is certified, not warned of", {
  d <- ofail_data(N = 200L, scale = 1)
  plain <- suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "z")))
  local_mocked_bindings(run_nloptr = fail_after_run, .package = "choicer")
  expect_message(
    fit <- expect_no_warning(run_mnlogit(d, "id", "alt", "choice", c("x", "z"))),
    "a Newton phase of no steps", fixed = TRUE)
  # the estimates and the variance are the optimizer's point's (two threads:
  # two runs of one fit agree to rounding only)
  expect_equal(coef(fit), coef(plain), tolerance = 1e-8)
  expect_equal(fit$vcov, plain$vcov, tolerance = 1e-8)
  expect_identical(fit$convergence, 1L)
  expect_true(choicer:::.is_converged(fit))
  expect_match(fit$message, "^Newton polish converged after 0 iterations")
  expect_match(fit$message, "NLOPT_FAILURE", fixed = TRUE)
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$optimizer$polish$iterations, 0L)
  expect_identical(fit$optimizer$polish$optimizer$convergence, -1L)

  nl <- suppressMessages(expect_no_warning(
    run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest",
                  scale_vars = "bhhh", control = list(xtol_rel = 1e-14))))
  expect_identical(nl$convergence, 1L)
  expect_identical(nl$optimizer$polish$iterations, 0L)
})

test_that("a polish that does not certify the point keeps the code and warns", {
  d <- ofail_data(N = 200L, scale = 1)
  local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
  # certified from the start values: the code is the polish's
  fit <- expect_no_warning(suppressMessages(
    run_mnlogit(d, "id", "alt", "choice", c("x", "z"), polish = "newton")))
  expect_identical(fit$convergence, 1L)
  expect_gt(fit$optimizer$polish$iterations, 0L)
  # one step is not enough: moved, so not "the start values"
  expect_warning(
    short <- suppressMessages(
      run_mnlogit(d, "id", "alt", "choice", c("x", "z"), polish = "newton",
                  polish_control = list(max_iter = 1L))),
    "not certified as a maximum", class = "choicer_optimizer_failed")
  expect_identical(short$convergence, -1L)
  expect_identical(short$optimizer$polish$iterations, 1L)
})

test_that("a fit without a failure code takes no check", {
  d <- ofail_data(N = 200L, scale = 1)
  fit <- expect_no_warning(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "z"))))
  expect_null(fit$optimizer$polish)
  expect_true(fit$convergence %in% 1:4)
})

test_that("run_mxlogit() checks the optimizer's stop, the profile likelihood and \"blp\" too", {
  d <- ofail_data(N = 150L, scale = 1)
  args <- list(data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
               covariate_cols = "z", random_var_cols = "x", S = 20L, seed = 3L)
  fit_with <- function(...) suppressMessages(do.call("run_mxlogit", c(args, list(...))))

  ok <- fit_with(polish = "newton")
  expect_true(ok$optimizer$polish$converged)
  th <- unname(coef(ok))

  # at the start values, far from the maximum: warned
  local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
  expect_warning(far <- fit_with(), "returned the start values",
                 class = "choicer_optimizer_failed")
  expect_identical(far$convergence, -1L)
  # at the maximum: certified, and the estimates the start values there
  at_max <- expect_no_warning(fit_with(theta_init = th))
  expect_identical(at_max$convergence, 1L)
  expect_identical(at_max$optimizer$polish$iterations, 0L)
  expect_equal(unname(coef(at_max)), th, tolerance = 1e-12)

  # the profile likelihood's joint point
  prof <- expect_no_warning(fit_with(theta_init = th,
                                     concentrate_out_asc = "profile_likelihood"))
  expect_identical(prof$convergence, 1L)

  # "blp" takes no polish; its outer decrement judges the point
  local_mocked_bindings(run_nloptr = fail_after_run, .package = "choicer")
  blp_ok <- expect_no_warning(fit_with(concentrate_out_asc = "blp"))
  expect_null(blp_ok$optimizer$polish)
  expect_lt(blp_ok$optimizer$concentration$decrement, 1e-8)
  expect_identical(blp_ok$convergence, 1L)
  expect_match(blp_ok$message, "NLOPT_FAILURE", fixed = TRUE)

  local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
  expect_warning(blp_far <- fit_with(concentrate_out_asc = "blp"), "returned the start values",
                 class = "choicer_optimizer_failed")
  expect_identical(blp_far$convergence, -1L)
})

test_that("optim's abnormal line search is checked as well", {
  d <- ofail_data(N = 200L, scale = 1)
  fail_optim <- function(theta_init, eval_f, lower, upper, control) {
    list(par = theta_init, value = eval_f(theta_init)$objective, convergence = 52L,
         message = "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH", iterations = 3L, raw = list())
  }
  local_mocked_bindings(run_optim = fail_optim, .package = "choicer")
  expect_warning(
    fit <- suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "z"), optimizer = "optim")),
    "ABNORMAL_TERMINATION_IN_LNSRCH", class = "choicer_optimizer_failed")
  expect_identical(fit$convergence, 52L)
  # a custom optimizer's code is its own
  custom <- function(theta_init, eval_f, lower, upper) {
    list(par = theta_init, value = eval_f(theta_init)$objective, convergence = -1L)
  }
  expect_no_warning(suppressMessages(
    run_mnlogit(d, "id", "alt", "choice", c("x", "z"), optimizer = custom)))
})
