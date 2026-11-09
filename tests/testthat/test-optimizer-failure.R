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

test_that("run_nestlogit() warns where L-BFGS fails at the start, and its default avoids it", {
  skip_on_cran()  # the optimizer's real path
  d <- ofail_data()
  # "none", the parameters' own units: a first step one gradient long
  expect_warning(
    fit <- suppressMessages(run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest",
                                          scale_vars = "none")),
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
  # the default, "bhhh", and a polish from "none" reach the same maximum,
  # without a warning
  expect_identical(eval(formals(run_nestlogit)$scale_vars)[1L], "bhhh")
  ok <- suppressMessages(expect_no_warning(
    run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest", polish = "newton")))
  expect_identical(ok$scale_vars, "bhhh")
  expect_true(ok$optimizer$polish$converged)
  expect_gt(as.numeric(logLik(ok)), as.numeric(logLik(fit)) + 1)
  pol2 <- suppressMessages(run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest",
                                         scale_vars = "none", polish = "newton"))
  expect_equal(as.numeric(logLik(pol2)), as.numeric(logLik(ok)), tolerance = 1e-8)
  expect_true(pol2$optimizer$polish$converged)
  expect_identical(pol2$optimizer$polish$optimizer$convergence, -1L)
})

test_that("the NL's scales do not freeze a dissimilarity started near its bound", {
  skip_on_cran()  # the optimizer's real path
  set.seed(11)
  N <- 800L
  J <- 5L
  id <- rep(seq_len(N), each = J)
  x1 <- rnorm(N * J)
  x2 <- 100 * rnorm(N * J)
  u <- 0.8 * x1 - 0.01 * x2 + rep(c(0, 0.3, -0.2, 0.5, 0.1), N) - log(-log(runif(N * J)))
  d <- data.table::data.table(id = id, alt = rep(seq_len(J), N), x1 = x1, x2 = x2,
                              nest = rep(c(1L, 1L, 2L, 2L, 3L), N),
                              choice = as.integer(u == ave(u, id, FUN = max)))
  nl <- function(...) suppressMessages(
    run_nestlogit(d, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest", ...))
  best <- nl(polish = "newton")
  expect_true(best$optimizer$polish$converged)
  m <- coef(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x1", "x2"))))
  lam <- best$param_map$lambda
  start <- numeric(length(coef(best)))
  start[best$param_map$beta] <- m[c("x1", "x2")]
  start[best$param_map$asc] <- m[grep("^ASC", names(m))]
  for (low in c(1e-2, 1e-3, 1e-4)) {
    start[lam] <- low
    fit <- expect_no_warning(nl(theta_init = start))
    expect_equal(fit$loglik, best$loglik, tolerance = 1e-6)
  }
  # below that no coordinates save the fit, and it says so before the optimizer
  start[lam] <- 1e-6
  expect_warning(nl(theta_init = start), "A nest parameter starts below 1e-4", fixed = TRUE)
  # a start well inside is read as it is
  expect_equal(nl(theta_init = unname(coef(best)))$loglik, best$loglik, tolerance = 1e-6)
})

test_that("the NL's BHHH pass leaves the scales at 1 where the start overflows", {
  d <- ofail_data(N = 200L, scale = 1)
  real <- choicer:::nl_loglik_gradient_parallel
  local_mocked_bindings(nl_loglik_gradient_parallel = function(...) {
    r <- real(...)
    if (isTRUE(list(...)$opg_diag)) r$overflow <- TRUE
    r
  }, .package = "choicer")
  expect_message(
    fit <- run_nestlogit(d, "id", "alt", "choice", c("x", "z"), nest_col = "nest"),
    "overflows at the start values", fixed = TRUE)
  expect_true(all(fit$param_scale == 1))
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
  # (a success code, whatever the real optimizer's last step)
  local_mocked_bindings(run_nloptr = function(theta_init, eval_f, lower, upper, control) {
    r <- real_nloptr(theta_init, eval_f, lower, upper, control)
    r$convergence <- 3L
    r
  }, .package = "choicer")
  fit <- expect_no_warning(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "z"))))
  expect_null(fit$optimizer$polish)
  expect_identical(fit$convergence, 3L)
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

test_that("a \"blp\" point is not certified where the outer curvature is not positive definite", {
  d <- ofail_data(N = 150L, scale = 1)
  # at this start the curvature along the share constraint is indefinite, where
  # the decrement's formula can read negative
  local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
  expect_warning(
    fit <- suppressMessages(run_mxlogit(
      d, "id", "alt", "choice", "z", "x", S = 20L, seed = 3L, concentrate_out_asc = "blp",
      theta_init = c(-1.879, 0.367, 0, 0, 0))),
    class = "choicer_optimizer_failed")
  expect_identical(fit$convergence, -1L)
  expect_identical(fit$optimizer$concentration$decrement, NA_real_)
})

# Stationary points of "blp" that are minima of the log-likelihood along the
# share constraint, with the sample's own shares. "normal": one normal random
# coefficient, its standard deviation near zero, where the BHHH matrix of
# se_method = "bhhh" is positive definite too (its decrement, about 1e-27,
# certified the point). "lognormal": one log-normal coefficient, where the
# curvature along the constraint that the variance takes is positive definite
# under the Hessian too (+8.1, against the objective's -5.3): it omits the
# constants' score times the curvature of delta(theta_c), and that score is not
# zero where the shares match. Returns the fit from a start (`fit_at`), the
# point (`start`), the negated log-likelihood's second derivative along the
# constraint there (`curv`, by differences of the outer gradient) and the
# variance's curvature along the constraint (`It`, under the Hessian).
blp_stationary_min <- function(design = c("normal", "lognormal")) {
  design <- match.arg(design)
  lognormal <- design == "lognormal"
  set.seed(if (lognormal) 20L else 1L)
  N <- if (lognormal) 40L else 300L
  J <- 3L
  w <- rnorm(N * J)
  b <- if (lognormal) exp(rnorm(N, sd = 1)) else rnorm(N, sd = 2)
  u <- rep(b, each = J) * w + (if (lognormal) rep(c(0, 1, -1), N) else 0) -
    log(-log(runif(N * J)))
  id <- rep(seq_len(N), each = J)
  d <- data.table::data.table(id = id, alt = rep(seq_len(J), N), w = w,
                              choice = as.integer(u == ave(u, id, FUN = max)))
  fit_at <- function(start = NULL, se = "hessian", control = list()) {
    suppressMessages(choicer::run_mxlogit(
      d, "id", "alt", "choice", NULL, "w", S = 20L, draws = "generate", seed = 17L,
      rc_dist = as.integer(lognormal), concentrate_out_asc = "blp", theta_init = start,
      se_method = se, scale_vars = "none", control = control))
  }
  f <- suppressWarnings(fit_at(control = list(maxeval = 1L)))
  gp <- choicer:::.mxl_gen_params(f$draws_info, N = N)
  inner <- choicer:::.mxl_blp_setup(
    f$data, isTRUE(f$include_outside_option), gp$eta_draws,
    list(seed = gp$gen_seed, scramble = gp$gen_scramble, S = gp$gen_S),
    f$rc_dist, f$rc_correlation, f$rc_mean, f$param_map, f$n_params,
    f$share_targets, choicer:::.concentrate_control(list(tol_shares = 1e-13), "blp"))
  solve_at <- function(x) inner$solve(x, unname(coef(f)[f$param_map$asc]), 4000L)
  x <- stats::uniroot(function(x) solve_at(x)$gradient[1],
                      if (lognormal) c(1.5, 1.75) else c(-6, -4), tol = 1e-12)$root
  curv <- (solve_at(x + 1e-4)$gradient[1] - solve_at(x - 1e-4)$gradient[1]) / 2e-4
  f$coefficients[] <- c(x, solve_at(x)$delta)
  list(fit_at = fit_at, start = unname(coef(f)), curv = curv,
       It = choicer:::.blp_fit_vcov(f, "hessian")$It[1, 1])
}

test_that("a \"blp\" stationary point that is not a maximum is not certified under any se_method", {
  for (design in c("normal", "lognormal")) {
    pt <- blp_stationary_min(design)
    expect_lt(pt$curv, 0)
    # (the variance's curvature along the constraint reads positive at the
    # log-normal's minimum)
    if (design == "lognormal") expect_gt(pt$It, 0)
    local_mocked_bindings(run_nloptr = fail_at_start, .package = "choicer")
    for (se in c("hessian", "bhhh")) {
      expect_warning(fit <- pt$fit_at(pt$start, se), class = "choicer_optimizer_failed")
      expect_identical(fit$convergence, -1L)
      expect_identical(fit$optimizer$concentration$decrement, NA_real_)
    }
  }
})

test_that("a \"blp\" maximum is certified by the objective's own curvature under any se_method", {
  d <- ofail_data(N = 150L, scale = 1)
  fit_with <- function(...) {
    suppressMessages(run_mxlogit(d, "id", "alt", "choice", "z", "x", S = 20L, seed = 3L,
                                 concentrate_out_asc = "blp", keep_data = TRUE, ...))
  }
  local_mocked_bindings(run_nloptr = fail_after_run, .package = "choicer")
  hess <- expect_no_warning(fit_with())
  bhhh <- expect_no_warning(fit_with(se_method = "bhhh"))
  expect_identical(hess$convergence, 1L)
  expect_identical(bhhh$convergence, 1L)
  expect_lt(hess$optimizer$concentration$decrement, 1e-8)
  # the same check whatever the standard errors
  expect_equal(bhhh$optimizer$concentration$decrement,
               hess$optimizer$concentration$decrement, tolerance = 1e-6)
  # and where the variance fails
  local_mocked_bindings(.blp_fit_vcov = function(...) stop("no variance here"),
                        .package = "choicer")
  expect_warning(novc <- fit_with(), "no variance here", class = "choicer_vcov_failed")
  expect_identical(novc$convergence, 1L)
  expect_equal(novc$optimizer$concentration$decrement,
               hess$optimizer$concentration$decrement, tolerance = 1e-6)
})

test_that(".blp_outer_hessian() differences the outer gradient", {
  # gradient A theta + theta^2 / 2 (and a constant's entry, dropped): Hessian
  # A + diag(theta), which central differences of a quadratic reach exactly
  A <- matrix(c(4, 1, 0, 1, 3, -1, 0, -1, 2), 3L)
  th0 <- c(0.3, -1, 2)
  calls <- 0L
  solve <- function(theta_c, delta0, max_iter, best) {
    calls <<- calls + 1L
    stopifnot(identical(best$theta_c, th0), identical(delta0, 0.5))
    list(converged = TRUE, gradient = c(A %*% theta_c + theta_c^2 / 2, 99))
  }
  oh <- choicer:::.blp_outer_hessian(solve, th0, list(delta = 0.5, extra = NULL),
                                     rep(1e-3, 3L), 10L)
  expect_null(oh$why)
  expect_equal(oh$hessian, A + diag(th0), tolerance = 1e-10)
  expect_identical(calls, 6L)
  # a solve that fails leaves no Hessian, and says why
  fails <- function(...) list(converged = FALSE, message = "the share inversion stopped")
  oh <- choicer:::.blp_outer_hessian(fails, th0, list(delta = 0.5), rep(1e-3, 3L), 10L)
  expect_null(oh$hessian)
  expect_identical(oh$why, "the share inversion stopped")
})

test_that("the start values are named only where the estimates are them", {
  warned <- function(par, start, scale = NULL) {
    w <- tryCatch(
      choicer:::.warn_optimizer_failed(list(par = par, convergence = -1L, message = "M"),
                                       "nloptr", start, NULL, "bhhh", FALSE, scale = scale),
      warning = function(c) conditionMessage(c))
    grepl("returned the start values", w, fixed = TRUE)
  }
  expect_true(warned(c(0, 0.5), c(0, 0.5)))
  # the coordinates' round trip, within rounding of the parameter's scale
  expect_true(warned(c(1e-25, 0.5 * (1 + 1e-15)), c(0, 0.5), scale = c(1e-12, 1)))
  # a small move on a parameter of small scale is a move
  expect_false(warned(c(7e-13, 0.5), c(0, 0.5), scale = c(1e-12, 1)))
  # without coordinates, any move is one
  expect_false(warned(c(1e-20, 0.5), c(0, 0.5)))
  # a point that is not finite is not the start
  expect_false(warned(c(NaN, NaN), c(0, 0.5)))
})

test_that("the warning reads without an optimizer's message", {
  w <- tryCatch(
    choicer:::.warn_optimizer_failed(list(par = 1, convergence = 10L, message = ""),
                                     "optim", 0, NULL, "bhhh", FALSE),
    warning = function(c) conditionMessage(c))
  expect_match(w, "a failure (code 10)", fixed = TRUE)
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
