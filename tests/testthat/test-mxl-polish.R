# run_mxlogit(polish = ): the Newton phase on the mixed logit. A converged
# polish certifies the maximum the optimizer approached (a decrement below
# tol_decrement) and reports the optimizer's success code, keeping the
# optimizer's own status in fit$optimizer$polish; polish = "none" (the
# default) changes nothing.

# A J = 30 cross-section; simulate_mxl_data()'s own outside rows are dropped
# and the preparation adds the outside option. Every alternative is chosen,
# none in every situation that offers it.
mxp_data <- function(N = 600L, seed = 11L) {
  sim <- simulate_mxl_data(N = N, J = 30L, beta = c(0.8, -0.6),
                           delta = seq(-1, 1, length.out = 30L),
                           Sigma = diag(c(1, 0.5)), outside_option = TRUE,
                           seed = seed)
  data.table::as.data.table(sim$data)[alt != 0L]
}

mxp_fit <- function(dt, ...) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 20L, include_outside_option = TRUE, outside_opt_label = 0L, ...))
}

# A fit's distance from a reference: the largest difference in a
# coefficient, over its absolute value or its standard error, whichever is
# larger
mxp_dist <- function(fit, ref) {
  max(abs(coef(fit) - coef(ref)) / pmax(abs(coef(ref)), ref$se))
}

MXP_TIGHT <- list(xtol_rel = 1e-14, maxeval = 5000L)
MXP_LOOSE <- list(maxeval = 6L)
# Enough steps for a polish expected to converge to do so on any platform
MXP_STEPS <- list(max_iter = 100L)
MXP_CERT <- list(max_iter = 100L, tol_decrement = 1e-16)

# The fit without its call and run time, for bitwise comparisons
mxp_strip <- function(fit) {
  fit$call <- NULL
  fit$optimizer$elapsed_time <- NULL
  fit
}

test_that("polish = \"none\" is the fit without a polish, bit for bit", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  f0 <- mxp_fit(dt, control = MXP_LOOSE)
  f1 <- mxp_fit(dt, control = MXP_LOOSE, polish = "none",
                polish_control = list(max_iter = 3L))
  expect_true(identical(mxp_strip(f1), mxp_strip(f0), num.eq = FALSE,
                        single.NA = FALSE))
  expect_named(f0$optimizer, c("name", "control", "elapsed_time", "iterations"))
})

test_that("a loose fit, polished, reaches the tightly converged maximum", {
  skip_on_cran()
  dt <- mxp_data()
  ref <- mxp_fit(dt, control = MXP_TIGHT)
  # the reference's own certificate: max_iter = 0 does not move it (a second
  # fit at two threads differs from the first in the last bits)
  cert <- mxp_fit(dt, control = MXP_TIGHT, polish = "newton",
                  polish_control = list(max_iter = 0L))
  expect_true(cert$optimizer$polish$converged)
  expect_lt(cert$optimizer$polish$decrement, 1e-10)
  expect_equal(coef(cert), coef(ref), tolerance = 1e-10)
  expect_identical(cert$optimizer$polish$iterations, 0L)

  for (pm in c("newton", "bhhh")) {
    fit <- mxp_fit(dt, control = MXP_LOOSE, polish = pm,
                   polish_control = list(max_iter = 100L,
                                         tol_decrement = 1e-16))
    p <- fit$optimizer$polish
    expect_true(p$converged, label = pm)
    expect_identical(p$method, pm)
    expect_gt(p$iterations, 0L)
    expect_lt(p$decrement, 1e-16)
    # the optimizer's own status: nloptr stopped at maxeval
    expect_identical(p$optimizer$convergence, 5L)
    expect_identical(fit$convergence, 1L)
    expect_match(fit$message, "polish converged after")
    expect_match(fit$message, "NLOPT_MAXEVAL_REACHED", fixed = TRUE)
    expect_identical(p$message, paste("the", if (pm == "newton") "Newton" else
                                        "BHHH", "decrement is below tol_decrement"))
    expect_gt(p$improvement, 0)
    expect_lt(mxp_dist(fit, ref), 1e-6)
    expect_equal(fit$loglik, ref$loglik, tolerance = 1e-12)
    expect_equal(fit$se, ref$se, tolerance = 1e-5)
  }
})

test_that("the record, print() and summary() report the polish", {
  dt <- mxp_data()
  fit <- mxp_fit(dt, control = MXP_LOOSE, polish = "newton",
                 polish_control = MXP_STEPS)
  p <- fit$optimizer$polish
  expect_named(p, c("method", "converged", "iterations", "damped",
                    "decrement", "lambda", "message", "evaluations", "active",
                    "held", "flat", "improvement", "error", "elapsed_time",
                    "optimizer"))
  expect_true(p$converged)
  expect_null(p$error)
  expect_length(p$active, 0L)
  expect_length(p$held, 0L)
  expect_length(p$flat, 0L)
  expect_named(p$evaluations, c("objective", "curvature"))
  expect_named(p$optimizer, c("convergence", "message", "iterations",
                              "elapsed_time"))
  expect_gte(p$elapsed_time, 0)
  # the fit's run time is the optimizer's and the polish's
  expect_equal(fit$optimizer$elapsed_time,
               p$optimizer$elapsed_time + p$elapsed_time)
  expect_output(print(fit), "Newton polish converged")
  expect_output(print(summary(fit)), "Newton polish converged")
})

test_that("a polish that does not converge keeps the optimizer's code", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  f0 <- mxp_fit(dt, control = MXP_LOOSE)
  suppressMessages(expect_message(
    f1 <- run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      S = 20L, include_outside_option = TRUE, outside_opt_label = 0L,
      control = MXP_LOOSE, polish = "newton",
      polish_control = list(max_iter = 0L)),
    paste0("^Newton polish: not converged after 0 iterations \\(Newton ",
           "decrement [^)]+\\): the iteration limit \\(max_iter\\) is reached, ")))
  p <- f1$optimizer$polish
  expect_false(p$converged)
  expect_identical(f1$convergence, f0$convergence)
  expect_match(f1$message,
               "Newton polish not converged after 0 iterations: the iteration limit")
  # max_iter = 0 only evaluates the decrement: the fit is the optimizer's
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$loglik, f0$loglik)
  expect_identical(f1$vcov, f0$vcov)
  expect_gt(p$decrement, 1e-8)
})

test_that("a binding bound stays active, and the polish certifies the rest", {
  skip_on_cran()
  dt <- mxp_data()
  ref <- mxp_fit(dt, upper = c(x1 = 0.5), control = MXP_TIGHT)
  fit <- mxp_fit(dt, upper = c(x1 = 0.5), control = MXP_LOOSE,
                 polish = "newton", polish_control = MXP_CERT)
  p <- fit$optimizer$polish
  expect_true(p$converged)
  expect_identical(p$active, "x1")
  expect_identical(coef(fit)[["x1"]], 0.5)
  expect_lt(mxp_dist(fit, ref), 1e-6)
})

test_that("constants without a finite estimate are held at the optimizer's", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  # alternative 30 is never chosen
  dt[alt == 30L, choice := 0L]
  expect_warning(f0 <- mxp_fit(dt, control = list(maxeval = 15L)),
                 class = "choicer_unbounded_constants")
  for (pm in c("newton", "bhhh")) {
    expect_warning(fit <- mxp_fit(dt, control = list(maxeval = 15L),
                                  polish = pm, polish_control = MXP_STEPS),
                   class = "choicer_unbounded_constants")
    p <- fit$optimizer$polish
    expect_true(p$converged, label = pm)
    expect_identical(p$held, "ASC_30")
    expect_identical(coef(fit)[["ASC_30"]], coef(f0)[["ASC_30"]])
  }
})

test_that("a variance that vanishes neither stalls Newton nor runs BHHH away", {
  # In mode_choice the travel-time coefficient's variance goes to zero: its
  # log-Cholesky diagonal heads toward -Inf, where BHHH's curvature falls
  # like exp(4 L_11) and Newton's like exp(2 L_11). BHHH's step in L_11 runs
  # past ten maximum steps: L_11 is held, so that it neither drifts to where
  # its information underflows (every standard error would then be NA) nor,
  # shortened, shortens the other parameters' steps (they would barely move)
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  mc <- data.table::as.data.table(mode_choice)
  fit_mc <- function(...) {
    suppressMessages(run_mxlogit(
      data = mc, id_col = "id", alt_col = "mode", choice_col = "choice",
      covariate_cols = c("vcost", "wait"), random_var_cols = "travel",
      S = 50L, control = list(maxeval = 10L), ...))
  }
  fn <- fit_mc(polish = "newton", polish_control = MXP_STEPS)
  expect_true(fn$optimizer$polish$converged)
  expect_lt(coef(fn)[["L_11"]], -10)
  fb <- fit_mc(polish = "bhhh", polish_control = MXP_STEPS)
  pb <- fb$optimizer$polish
  expect_false(pb$converged)
  expect_identical(pb$flat, "L_11")
  expect_match(pb$message, "a curvature all but vanished", fixed = TRUE)
  expect_gt(coef(fb)[["L_11"]], -30)
  # the other parameters reach the maximum given L_11, within 1e-2 of
  # Newton's log-likelihood
  expect_gt(pb$improvement, 0)
  expect_lt(fn$loglik - fb$loglik, 1e-2)
  expect_false(anyNA(fb$se))
})

test_that("\"bhhh\" with fewer likelihood units than parameters stops", {
  sim <- simulate_mxl_data(N = 3L, T = 5L, J = 3L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 5L)
  dt <- data.table::as.data.table(sim$data)
  # 3 decision makers, 4 parameters without constants
  expect_error(suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 10L, use_asc = FALSE, person_col = "pid", polish = "bhhh")),
    "the BHHH matrix of 3 decision makers has rank at most 3")
})

test_that("the BHHH rank check counts only the parameters the polish moves", {
  sim <- simulate_mxl_data(N = 3L, T = 5L, J = 3L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 5L)
  dt <- data.table::as.data.table(sim$data)
  # x1 pinned by lower == upper: 3 parameters move, one per decision maker
  fit <- suppressWarnings(suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 10L, use_asc = FALSE, person_col = "pid", polish = "bhhh",
    lower = c(x1 = 0.8), upper = c(x1 = 0.8))))
  expect_identical(fit$optimizer$polish$method, "bhhh")
  expect_null(fit$optimizer$polish$error)
  expect_identical(coef(fit)[["x1"]], 0.8)
})

# Half of g' C^-1 g at the fit's estimates, over the mean weight, with g the
# gradient and C the Hessian (or BHHH) computed afresh from the fit's stored
# data and draws
mxp_decrement <- function(fit, curvature) {
  d <- fit$data
  gp <- .mxl_gen_params(fit$draws_info, N = length(.unit_first(d)))
  g <- mxl_loglik_gradient_parallel(
    theta = fit$coefficients, X = d$X, W = d$W, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
    rc_dist = fit$rc_dist, rc_correlation = fit$rc_correlation,
    rc_mean = fit$rc_mean, eta_draws = gp$eta_draws, use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option,
    gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
    Ti = d$Ti)$gradient
  fit$se_method <- if (curvature == "bhhh") "bhhh" else "hessian"
  0.5 * sum(g * solve(compute_hessian(fit), g)) / mean(d$weights)
}

test_that("the polish's curvature takes the fit's draws, weights and panel", {
  # With max_iter = 0 the polish only evaluates the decrement at the
  # optimizer's solution. A curvature on other draws (another mode or
  # scramble), on squared weights or per choice situation would give
  # another decrement there.
  sim <- simulate_mxl_data(N = 150L, T = 4L, J = 5L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 7L)
  dt <- data.table::as.data.table(sim$data)
  set.seed(8)
  wt <- stats::runif(150L, 0.5, 2)
  dt[, w := wt[pid]]
  fit_p <- function(...) {
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      S = 20L, person_col = "pid", weights_col = "w",
      se_method = "sandwich", ...))
  }
  modes <- list(list(draws = "generate", seed = 3L),
                list(draws = "generate", seed = 3L, scramble = "none"),
                list(draws = "store"))
  for (dm in modes) {
    for (pm in c("newton", "bhhh")) {
      lab <- paste(c(dm$draws, dm$scramble, pm), collapse = " ")
      fit <- do.call(fit_p, c(dm, list(control = list(maxeval = 10L),
                                       polish = pm,
                                       polish_control = list(max_iter = 0L))))
      p <- fit$optimizer$polish
      expect_identical(p$lambda, 0, label = lab)
      expect_gt(p$decrement, 1e-6, label = lab)
      expect_equal(p$decrement, mxp_decrement(fit, pm), tolerance = 1e-8,
                   label = lab)
    }
  }
  # and without constants
  for (pm in c("newton", "bhhh")) {
    fit <- fit_p(draws = "generate", seed = 3L, use_asc = FALSE,
                 control = list(maxeval = 10L), polish = pm,
                 polish_control = list(max_iter = 0L))
    p <- fit$optimizer$polish
    expect_null(p$error, label = pm)
    expect_identical(p$lambda, 0, label = pm)
    expect_equal(p$decrement, mxp_decrement(fit, pm), tolerance = 1e-8,
                 label = pm)
  }

  skip_on_cran()
  ref <- fit_p(draws = "generate", seed = 3L, control = MXP_TIGHT)
  fit <- fit_p(draws = "generate", seed = 3L, control = MXP_LOOSE,
               polish = "newton", polish_control = MXP_CERT)
  expect_true(fit$optimizer$polish$converged)
  expect_lt(mxp_dist(fit, ref), 1e-6)
})

test_that("the polish works in natural units, whatever the optimizer's coordinates", {
  skip_on_cran()
  dt <- mxp_data()
  ref <- mxp_fit(dt, control = MXP_TIGHT)
  for (sv in c("none", "sd")) {
    fit <- mxp_fit(dt, control = MXP_LOOSE, scale_vars = sv,
                   polish = "newton", polish_control = MXP_CERT)
    expect_true(fit$optimizer$polish$converged, label = sv)
    expect_lt(mxp_dist(fit, ref), 1e-6)
  }
})

test_that("keep_data = FALSE and the advanced pathway polish the same fit", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  fit <- mxp_fit(dt, control = MXP_LOOSE, polish = "newton")
  fk <- mxp_fit(dt, control = MXP_LOOSE, polish = "newton", keep_data = FALSE)
  expect_null(fk$data)
  expect_identical(coef(fk), coef(fit))
  expect_identical(fk$vcov, fit$vcov)
  inp <- prepare_mxl_data(dt, id_col = "id", alt_col = "alt",
                          choice_col = "choice", covariate_cols = c("x1", "x2"),
                          random_var_cols = c("w1", "w2"),
                          include_outside_option = TRUE, outside_opt_label = 0L)
  fa <- suppressMessages(run_mxlogit(
    input_data = inp, eta_draws = get_halton_normals(20L, length(inp$M), 2L),
    control = MXP_LOOSE, polish = "newton"))
  expect_identical(coef(fa), coef(fit))
  expect_identical(fa$optimizer$polish$decrement,
                   fit$optimizer$polish$decrement)
  expect_identical(fa$vcov, fit$vcov)
})

test_that("the polish gets the fit's coordinates, bounds and settings", {
  dt <- mxp_data()
  polish_fn <- run_polish
  seen <- NULL
  local_mocked_bindings(run_polish = function(theta, eval_f, curv_f, ...) {
    seen <<- list(...)
    polish_fn(theta, eval_f, curv_f, ...)
  })
  fit <- mxp_fit(dt, control = MXP_LOOSE, upper = c(x1 = 2), polish = "bhhh",
                 polish_control = list(max_iter = 3L))
  n <- length(coef(fit))
  expect_identical(unname(seen$scale), unname(fit$param_scale))
  expect_false(all(seen$scale == 1))
  expect_identical(unname(seen$lower), rep(-Inf, n))
  expect_identical(unname(seen$upper), c(2, rep(Inf, n - 1L)))
  expect_identical(seen$control, .polish_control(list(max_iter = 3L)))
  expect_identical(seen$method, "bhhh")
  expect_identical(seen$fixed, integer())
})

test_that("correlated, shifted and log-normal random coefficients are polished", {
  skip_on_cran()
  sim <- simulate_mxl_data(N = 600L, J = 10L, beta = c(0.8, -0.6),
                           mu = c(0.5, -0.5),
                           Sigma = matrix(c(1, 0.3, 0.3, 0.5), 2L),
                           rc_dist = c(0L, 1L), rc_correlation = TRUE,
                           outside_option = TRUE, seed = 13L)
  dt <- data.table::as.data.table(sim$data)[alt != 0L]
  fit_l <- function(...) {
    mxp_fit(dt, rc_dist = c(0L, 1L), rc_correlation = TRUE, rc_mean = TRUE,
            ...)
  }
  ref <- fit_l(control = MXP_TIGHT)
  fit <- fit_l(control = MXP_LOOSE, polish = "newton",
               polish_control = MXP_CERT)
  expect_true(fit$optimizer$polish$converged)
  expect_lt(mxp_dist(fit, ref), 1e-6)
})

test_that("optim and a custom optimizer hand their solution to the polish", {
  skip_on_cran()
  dt <- mxp_data()
  ref <- mxp_fit(dt, control = MXP_TIGHT)
  fit <- mxp_fit(dt, optimizer = "optim", control = list(maxit = 3L),
                 polish = "newton", polish_control = MXP_CERT)
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$convergence, 0L)
  expect_identical(fit$optimizer$polish$optimizer$convergence, 1L)
  expect_lt(mxp_dist(fit, ref), 1e-6)

  # a custom optimizer that returns its start, without a message: the
  # polish starts there, and the fit's message is the polish's alone
  hold <- function(theta_init, eval_f, ...) {
    list(par = theta_init, value = eval_f(theta_init)$objective,
         convergence = 9L)
  }
  fit <- mxp_fit(dt, optimizer = hold, polish = "newton",
                 polish_control = MXP_CERT)
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$convergence, 0L)
  expect_identical(fit$optimizer$polish$optimizer$convergence, 9L)
  expect_identical(fit$optimizer$polish$optimizer$message, "")
  expect_match(fit$message,
               "^Newton polish converged after [0-9]+ iterations \\(Newton decrement [^)]+\\)$")
  expect_lt(mxp_dist(fit, ref), 1e-6)
})

test_that("an optimizer's message of any length leads into the polish's", {
  dt <- mxp_data(N = 200L)
  msg_opt <- function(msg) {
    function(theta_init, eval_f, ...) {
      list(par = theta_init, value = eval_f(theta_init)$objective,
           convergence = 0L, message = msg)
    }
  }
  fit <- suppressWarnings(mxp_fit(dt, optimizer = msg_opt(character(0)),
                                  polish = "newton",
                                  polish_control = list(max_iter = 0L)))
  expect_match(fit$message, "^Newton polish not converged after 0 iterations")
  fit <- suppressWarnings(mxp_fit(dt, optimizer = msg_opt(c("one", "two")),
                                  polish = "newton",
                                  polish_control = list(max_iter = 0L)))
  expect_match(fit$message, "^one; two; Newton polish not converged")
})

# A weighted panel with clusters of decision makers, for the variance after a
# polish: unit weights would hide B(w) standing in for the robust meat B(w^2)
mxp_panel <- function() {
  sim <- simulate_mxl_data(N = 120L, T = 3L, J = 5L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 21L)
  dt <- data.table::as.data.table(sim$data)
  set.seed(22)
  wt <- stats::runif(120L, 0.5, 2)
  dt[, w := wt[pid]]
  dt[, cl := pid %% 15L]
  dt[]
}
# Non-uniform weights make hessian/bhhh fits warn that they are not WESML
# corrections; that advice is beside the point here, and only it is muffled
mxp_panel_fit <- function(dt, se, ...) {
  withCallingHandlers(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      S = 20L, person_col = "pid", weights_col = "w", se_method = se,
      cluster_col = if (se == "cluster") "cl", ...)),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "Non-uniform weights detected")) {
        invokeRestart("muffleWarning")
      }
    })
}

test_that("the polish's decrement rule does not depend on the weights' scale", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_panel()
  dt16 <- data.table::copy(dt)[, w := 16 * w]
  hold <- function(theta_init, eval_f, ...) {
    list(par = theta_init, value = eval_f(theta_init)$objective,
         convergence = 0L)
  }
  # the polish from the same loose point at both scales, in the columns'
  # coordinates, which do not move with the weights
  start <- unname(coef(mxp_panel_fit(dt, "hessian", control = list(maxeval = 25L))))
  run <- function(d, ctl) {
    mxp_panel_fit(d, "hessian", optimizer = hold, theta_init = start,
                  scale_vars = "sd", polish = "newton", polish_control = ctl)
  }
  # the decrement over the mean weight, bit for bit at a power of two
  nu <- run(dt, list(max_iter = 0L))$optimizer$polish$decrement
  expect_true(is.finite(nu) && nu > 0)
  expect_identical(run(dt16, list(max_iter = 0L))$optimizer$polish$decrement, nu)
  # a tolerance twice that certifies at once at both scales, where one in
  # raw log-likelihood units would not at 16 times the weights
  ctl <- list(tol_decrement = 2 * nu)
  a <- run(dt, ctl)$optimizer$polish
  b <- run(dt16, ctl)$optimizer$polish
  expect_true(a$converged)
  expect_identical(a$iterations, 0L)
  expect_identical(b[c("converged", "iterations", "decrement")],
                   a[c("converged", "iterations", "decrement")])
})

test_that("the polish keeps its curvature only for a variance that needs it", {
  dt <- mxp_panel()
  polish_fn <- run_polish
  kept <- NULL
  local_mocked_bindings(run_polish = function(...) {
    kept <<- list(...)$keep_curvature
    polish_fn(...)
  })
  for (pm in c("newton", "bhhh")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      kept <- NULL
      mxp_panel_fit(dt, se, control = list(maxeval = 10L), polish = pm,
                    polish_control = list(max_iter = 1L))
      want <- if (pm == "newton") se != "bhhh" else se == "bhhh"
      expect_identical(kept, want, label = paste(pm, se))
    }
  }
})

test_that("the variance reuses the polish's last curvature where se_method needs it", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_panel()
  # the real kernels, counted
  real_h <- mxl_hessian_parallel
  real_b <- mxl_bhhh_parallel
  n_h <- 0L
  n_b <- 0L
  local_mocked_bindings(
    mxl_hessian_parallel = function(...) {
      n_h <<- n_h + 1L
      real_h(...)
    },
    mxl_bhhh_parallel = function(...) {
      n_b <<- n_b + 1L
      real_b(...)
    })
  for (pm in c("newton", "bhhh")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      lab <- paste(pm, se)
      n_h <- 0L
      n_b <- 0L
      fit <- mxp_panel_fit(dt, se, control = list(maxeval = 10L), polish = pm,
                           polish_control = list(max_iter = 3L))
      k <- fit$optimizer$polish$evaluations[["curvature"]]
      expect_gt(k, 1L, label = lab)
      # Newton's Hessian is the bread of every variance but "bhhh"'s; BHHH's
      # matrix (weights unsquared) is the "bhhh" variance's, never the
      # sandwich's meat (weights squared, computed); the cluster meat has its
      # own kernel
      reuse <- if (pm == "newton") se != "bhhh" else se == "bhhh"
      want_h <- (pm == "newton") * k + (se != "bhhh" && !reuse)
      want_b <- (pm == "bhhh") * k + (se == "sandwich") +
        (se == "bhhh" && !reuse)
      expect_identical(c(n_h, n_b), as.integer(c(want_h, want_b)), label = lab)
      type <- switch(se, sandwich = "robust", se)
      expect_identical(suppressMessages(vcov(fit, type = type)), fit$vcov,
                       label = lab)
      expect_false(any(vapply(fit$optimizer$polish, is.matrix, TRUE)),
                   label = lab)
    }
  }
})

test_that("a polish that takes no step leaves the fit's variance as it was", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_panel()
  for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
    f0 <- mxp_panel_fit(dt, se, control = list(maxeval = 10L))
    for (pm in c("newton", "bhhh")) {
      lab <- paste(pm, se)
      f1 <- mxp_panel_fit(dt, se, control = list(maxeval = 10L), polish = pm,
                          polish_control = list(max_iter = 0L))
      expect_identical(coef(f1), coef(f0), label = lab)
      expect_identical(f1$vcov, f0$vcov, label = lab)
      expect_identical(f1$se, f0$se, label = lab)
    }
  }
})

test_that("a polish whose curvature fails leaves the variance its own", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_panel()
  real_h <- mxl_hessian_parallel
  n_h <- 0L
  # the Hessian fails at the polish's second evaluation: the polish ends at
  # the point it reached, whose curvature it never had (the failing call
  # cleared the last one; that run_polish() keeps no curvature from another
  # point is test-polish.R's)
  local_mocked_bindings(mxl_hessian_parallel = function(...) {
    n_h <<- n_h + 1L
    if (n_h == 2L) stop("cannot allocate the Hessian")
    real_h(...)
  })
  expect_warning(fit <- mxp_panel_fit(dt, "hessian",
                                      control = list(maxeval = 10L),
                                      polish = "newton"),
                 class = "choicer_polish_failed")
  expect_identical(fit$optimizer$polish$iterations, 1L)
  expect_identical(n_h, 3L)
  expect_identical(suppressMessages(vcov(fit, type = "hessian")), fit$vcov)
})

test_that("a polish that hands back no curvature leaves the variance its own", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_panel()
  real_h <- mxl_hessian_parallel
  real_polish <- run_polish
  n_h <- 0L
  # as when the polish's final guard returns the start
  local_mocked_bindings(
    mxl_hessian_parallel = function(...) {
      n_h <<- n_h + 1L
      real_h(...)
    },
    run_polish = function(...) {
      r <- real_polish(...)
      r$curvature_at_par <- NULL
      r
    })
  fit <- mxp_panel_fit(dt, "hessian", control = list(maxeval = 10L),
                       polish = "newton", polish_control = list(max_iter = 2L))
  k <- fit$optimizer$polish$evaluations[["curvature"]]
  expect_identical(n_h, k + 1L)
  expect_identical(suppressMessages(vcov(fit, type = "hessian")), fit$vcov)
})

test_that("a reused curvature whose inversion fails takes the fallback", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  f0 <- mxp_fit(dt, control = MXP_LOOSE, polish = "newton")
  real_h <- mxl_hessian_parallel
  n_h <- 0L
  local_mocked_bindings(
    equilibrated_inverse = function(A, D) {
      stop("cannot allocate a working copy of 2.3 Gb")
    },
    mxl_hessian_parallel = function(...) {
      n_h <<- n_h + 1L
      real_h(...)
    })
  expect_warning(f1 <- mxp_fit(dt, control = MXP_LOOSE, polish = "newton"),
                 "cannot allocate a working copy", class = "choicer_vcov_failed")
  # the variance evaluated no Hessian of its own: it inverted the polish's
  expect_identical(n_h, f1$optimizer$polish$evaluations[["curvature"]])
  expect_identical(coef(f1), coef(f0))
  expect_null(f1$vcov)
  expect_match(f1$variance_error, "cannot allocate a working copy")
})

test_that("an error in the curvature keeps the optimizer's estimates", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  # the variance takes the BHHH kernel, so only the polish fails
  f0 <- mxp_fit(dt, control = MXP_LOOSE, se_method = "bhhh")
  local_mocked_bindings(mxl_hessian_parallel = function(...) {
    stop("cannot allocate the Hessian")
  })
  expect_warning(f1 <- mxp_fit(dt, control = MXP_LOOSE, se_method = "bhhh",
                               polish = "newton"),
                 "Newton polish stopped with an error",
                 class = "choicer_polish_failed")
  p <- f1$optimizer$polish
  expect_identical(p$error, "cannot allocate the Hessian")
  expect_false(p$converged)
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$convergence, f0$convergence)
  expect_identical(f1$vcov, f0$vcov)
})

test_that("an error after accepted steps keeps the last accepted point", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- mxp_data()
  f0 <- mxp_fit(dt, control = MXP_LOOSE, se_method = "bhhh")
  # one step, then the iteration limit: the point the failing polish reaches
  f1 <- mxp_fit(dt, control = MXP_LOOSE, se_method = "bhhh", polish = "newton",
                polish_control = list(max_iter = 1L))
  hessian <- mxl_hessian_parallel
  calls <- 0L
  local_mocked_bindings(mxl_hessian_parallel = function(...) {
    calls <<- calls + 1L
    if (calls > 1L) stop("cannot allocate the Hessian")
    hessian(...)
  })
  expect_warning(f2 <- mxp_fit(dt, control = MXP_LOOSE, se_method = "bhhh",
                               polish = "newton"),
                 class = "choicer_polish_failed")
  p <- f2$optimizer$polish
  expect_identical(calls, 2L)
  expect_identical(p$iterations, 1L)
  expect_identical(p$error, "cannot allocate the Hessian")
  expect_false(p$converged)
  expect_identical(coef(f2), coef(f1))
  expect_identical(f2$loglik, f1$loglik)
  expect_gt(f2$loglik, f0$loglik)
  expect_identical(f2$vcov, f1$vcov)
})

test_that("a start where the objective is not finite is returned as it is", {
  dt <- mxp_data()
  # a custom optimizer that returns a NaN coefficient
  nan_opt <- function(theta_init, eval_f, ...) {
    theta_init[1L] <- NaN
    list(par = theta_init, value = 1, convergence = 9L)
  }
  fit <- suppressWarnings(mxp_fit(dt, optimizer = nan_opt, polish = "newton"))
  p <- fit$optimizer$polish
  expect_false(p$converged)
  expect_identical(p$message, "the objective is not finite at the start")
  expect_identical(p$iterations, 0L)
  # the optimizer's solution and value are kept
  expect_true(is.nan(coef(fit)[[1L]]))
  expect_identical(fit$loglik, -1)
  expect_identical(fit$convergence, 9L)
})

test_that("constants that diverge together leave the polish unconverged", {
  # the outside option never chosen: every inside constant rises together,
  # and none is held
  sim <- simulate_mxl_data(N = 600L, J = 30L, beta = c(0.8, -0.6),
                           delta = seq(-1, 1, length.out = 30L),
                           Sigma = diag(c(1, 0.5)), outside_option = TRUE,
                           seed = 11L)
  dt <- data.table::as.data.table(sim$data)
  dt <- dt[id %in% dt[alt != 0L & choice == 1L, id] & alt != 0L]
  fit_o <- function(...) {
    expect_warning(fit <- mxp_fit(dt, control = list(maxeval = 200L),
                                  polish = "newton", ...),
                   class = "choicer_unbounded_constants")
    fit
  }
  # a decrement below tol_decrement at the optimizer's solution
  fit <- fit_o(polish_control = list(max_iter = 0L, tol_decrement = 1e3))
  p <- fit$optimizer$polish
  expect_false(p$converged)
  expect_length(p$held, 0L)
  expect_match(p$message, paste0(
    "the Newton decrement is below tol_decrement, but the constants diverge ",
    "together (the outside option is never chosen), so there is no maximum ",
    "to certify"), fixed = TRUE)
  expect_identical(fit$convergence, p$optimizer$convergence)
  expect_false(fit_o(polish_control = MXP_STEPS)$optimizer$polish$converged)
  # an upper bound on the constant of an alternative that is chosen stops
  # their rise: the constrained maximum is certified
  fit <- fit_o(upper = c(ASC_15 = 2), polish_control = MXP_STEPS)
  p <- fit$optimizer$polish
  expect_true(p$converged)
  expect_identical(p$active, "ASC_15")
  expect_identical(p$message, "the Newton decrement is below tol_decrement")
  # so does holding the constant of an alternative chosen wherever offered
  # (a decrement below tol_decrement at the optimizer's solution)
  dt_held <- dt[!(alt == 30L & !id %in% dt[alt == 30L & choice == 1L, id])]
  expect_warning(fit <- mxp_fit(dt_held, control = list(maxeval = 6L),
                                polish = "newton",
                                polish_control = list(max_iter = 0L,
                                                      tol_decrement = 1e3)),
                 class = "choicer_unbounded_constants")
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$optimizer$polish$held, "ASC_30")
  # a bound on the constant of an alternative never chosen does not: the
  # others rise past it at no cost
  dt <- dt[!id %in% dt[alt == 30L & choice == 1L, id]]
  fit <- fit_o(lower = c(ASC_30 = -5), upper = c(ASC_30 = -5),
               polish_control = list(max_iter = 0L, tol_decrement = 1e3))
  expect_false(fit$optimizer$polish$converged)
  expect_match(fit$optimizer$polish$message,
               "(the outside option is never chosen)", fixed = TRUE)

  # the base alternative never chosen, without an outside option
  sim <- simulate_mxl_data(N = 600L, J = 6L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 12L)
  dt <- data.table::as.data.table(sim$data)
  dt <- dt[!id %in% dt[alt == 1L & choice == 1L, id]]
  expect_warning(fit <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 20L, control = list(maxeval = 100L), polish = "newton",
    polish_control = list(max_iter = 0L, tol_decrement = 1e3))),
    class = "choicer_unbounded_constants")
  expect_false(fit$optimizer$polish$converged)
  expect_match(fit$optimizer$polish$message,
               "(the base alternative is never chosen)", fixed = TRUE)

  # the base alternative chosen wherever it is offered: the other constants
  # fall together, and no bound is taken to stop them
  sim <- simulate_mxl_data(N = 900L, J = 6L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE,
                           seed = 13L)
  dt <- data.table::as.data.table(sim$data)
  dt <- dt[!id %in% dt[alt == 1L, id] | id %in% dt[alt == 1L & choice == 1L, id]]
  expect_warning(fit <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 20L, control = list(maxeval = 6L), upper = c(ASC_2 = 5),
    polish = "newton",
    polish_control = list(max_iter = 0L, tol_decrement = 1e3))),
    class = "choicer_unbounded_constants")
  expect_false(fit$optimizer$polish$converged)
  expect_match(fit$optimizer$polish$message,
               "(the base alternative is chosen wherever it is offered)",
               fixed = TRUE)
})

test_that("polish and polish_control are checked before the data", {
  expect_error(run_mxlogit(polish = "lbfgs"), "should be one of")
  expect_error(run_mxlogit(polish_control = list(maxit = 2)),
               "Unknown `polish_control` setting")
  expect_error(run_mxlogit(polish_control = list(max_iter = -1)),
               "whole number")
})
