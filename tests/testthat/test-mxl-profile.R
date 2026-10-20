# run_mxlogit(concentrate_out_asc = "profile_likelihood"): the constants
# concentrated out of the likelihood. The profile log-likelihood
# l*(theta_c) = max_delta l(theta_c, delta) has the joint maximizer, so the
# estimator does not change: the outer optimizer works on the continuous
# parameters theta_c alone, each evaluation solving the constants' first-order
# condition by Newton steps (.mxl_profile_solver()), its gradient by the
# envelope theorem; the variance is the full Hessian's at
# (theta_c_hat, delta_hat).

# A J = 12 cross-section with the outside option (simulate_mxl_data()'s own
# outside rows dropped; the preparation adds the outside option), or, with
# outside = FALSE, without it (the first alternative the base). Every
# alternative is chosen, none in every situation that offers it.
prf_data <- function(N = 400L, J = 12L, seed = 31L, T = 1L, outside = TRUE) {
  sim <- simulate_mxl_data(N = N, J = J, T = T, beta = c(0.8, -0.6),
                           delta = seq(-1, 1, length.out = J),
                           Sigma = diag(c(1, 0.5)), outside_option = outside,
                           seed = seed)
  dt <- data.table::as.data.table(sim$data)
  if (outside) dt <- dt[alt != 0L]
  dt[]
}

prf_fit <- function(dt, ..., outside = TRUE) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    S = 20L, include_outside_option = outside,
    outside_opt_label = if (outside) 0L, ...))
}

# A panel of 300 decision makers with four situations each, whose random
# coefficients' variances are interior (on the cross-section the second one's
# maximum can sit at zero, its log-Cholesky diagonal running off), for the
# comparisons of standard errors and BHHH's certificate
prf_panel <- function() prf_data(N = 300L, T = 4L, seed = 35L)
prf_panel_fit <- function(dt, ...) prf_fit(dt, person_col = "pid", ...)

PRF <- list(concentrate_out_asc = "profile_likelihood")
# A maximum both routes certify (a Newton polish to a tiny decrement)
PRF_CERT <- list(polish = "newton",
                 polish_control = list(max_iter = 100L, tol_decrement = 1e-14))
PRF_TIGHT <- list(xtol_rel = 1e-12, maxeval = 5000L)

# A fit's distance from a reference: the largest difference in a
# coefficient, over its absolute value or its standard error, whichever is
# larger
prf_dist <- function(fit, ref) {
  max(abs(coef(fit) - coef(ref)) / pmax(abs(coef(ref)), ref$se))
}

# The fit's kernels' gradient at its estimates, from its stored data and draws
prf_gradient <- function(fit) {
  d <- fit$data
  gp <- .mxl_gen_params(fit$draws_info, N = length(.unit_first(d)))
  mxl_loglik_gradient_parallel(
    unname(coef(fit)), d$X, d$W, d$alt_idx, d$choice_idx, d$M, d$weights,
    gp$eta_draws, fit$rc_dist, fit$rc_correlation, fit$rc_mean, TRUE,
    fit$include_outside_option, gp$gen_seed, gp$gen_scramble, gp$gen_S, d$Ti)
}

test_that("the profile fit is the joint maximum, constants and all", {
  dt <- prf_data()
  ref <- do.call(prf_fit, c(list(dt, control = PRF_TIGHT), PRF_CERT))
  expect_true(ref$optimizer$polish$converged)
  fit <- do.call(prf_fit, c(list(dt), PRF, PRF_CERT))
  expect_true(fit$optimizer$polish$converged)
  expect_lt(prf_dist(fit, ref), 1e-6)
  expect_lt(abs(fit$loglik / ref$loglik - 1), 1e-9)
  # without the polish, the profile's own solution
  f0 <- do.call(prf_fit, c(list(dt), PRF))
  expect_lt(prf_dist(f0, ref), 1e-4)
  expect_lt(abs(f0$loglik - ref$loglik), 1e-6)
  # the same structure as a joint fit, with the concentration's record
  expect_identical(names(coef(f0)), names(coef(ref)))
  expect_identical(f0$param_map, ref$param_map)
  rec <- f0$optimizer$concentration
  expect_named(rec, c("method", "evaluations", "solves", "failed",
                      "inner_iterations", "passes", "max_residual", "residual",
                      "decrement", "elapsed_time"))
  expect_identical(rec$method, "profile_likelihood")
  expect_identical(rec$failed, 0L)
  expect_gt(rec$solves, 1L)
  expect_gte(rec$passes, rec$solves)
  expect_named(f0$optimizer, c("name", "control", "elapsed_time", "iterations",
                               "concentration"))
})

test_that("the profile fit is the joint maximum across designs and optimizers", {
  skip_on_cran()
  base <- prf_data(outside = FALSE)
  panel <- prf_data(N = 150L, T = 3L, seed = 32L)
  set.seed(33)
  wt <- stats::runif(150L, 0.5, 2)
  panel[, w := wt[pid]]
  cases <- list(
    base = list(dt = base, outside = FALSE),
    panel_generate_weighted = list(dt = panel, person_col = "pid",
                                   draws = "generate", seed = 5L,
                                   weights_col = "w", se_method = "sandwich"),
    none = list(dt = prf_data(), scale_vars = "none"),
    sd = list(dt = prf_data(), scale_vars = "sd"),
    optim = list(dt = prf_data(), optimizer = "optim"),
    bound = list(dt = prf_data(), upper = c(x1 = 0.5)),
    correlated_means = list(dt = prf_data(), rc_correlation = TRUE, rc_mean = TRUE)
  )
  for (nm in names(cases)) {
    a <- cases[[nm]]
    dt <- a$dt
    a$dt <- NULL
    ctl <- if (identical(a$optimizer, "optim")) list(maxit = 2000L) else PRF_TIGHT
    ref <- do.call(prf_fit, c(list(dt, control = ctl), a, PRF_CERT))
    fit <- do.call(prf_fit, c(list(dt), a, PRF, PRF_CERT))
    expect_true(fit$optimizer$polish$converged, label = nm)
    expect_lt(prf_dist(fit, ref), 1e-6, label = nm)
    expect_lt(abs(fit$loglik / ref$loglik - 1), 1e-9, label = nm)
    if (nm == "bound") expect_identical(unname(coef(fit)["x1"]), 0.5)
  }
  # a custom optimizer gets the profile in the optimizer's coordinates
  dt <- prf_data()
  my_opt <- function(theta_init, eval_f, lower, upper, control) {
    o <- stats::optim(theta_init, function(x) eval_f(x)$objective,
                      function(x) eval_f(x)$gradient, method = "BFGS",
                      control = list(maxit = 500L, reltol = 1e-14))
    list(par = o$par, value = o$value, convergence = o$convergence)
  }
  ref <- do.call(prf_fit, c(list(dt, control = PRF_TIGHT), PRF_CERT))
  fit <- do.call(prf_fit, c(list(dt, optimizer = my_opt), PRF, PRF_CERT))
  expect_lt(prf_dist(fit, ref), 1e-6)
  # the advanced pathway, with the same draws, is the same fit (at one
  # thread, where the kernels' sums are reproducible)
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  inp <- prepare_mxl_data(dt, id_col = "id", alt_col = "alt",
                          choice_col = "choice", covariate_cols = c("x1", "x2"),
                          random_var_cols = c("w1", "w2"),
                          include_outside_option = TRUE, outside_opt_label = 0L)
  conv <- do.call(prf_fit, c(list(dt), PRF))
  adv <- suppressMessages(run_mxlogit(
    input_data = inp, eta_draws = get_halton_normals(20L, length(inp$M), 2L),
    concentrate_out_asc = "profile_likelihood"))
  expect_identical(coef(adv), coef(conv))
  expect_identical(adv$vcov, conv$vcov)
})

test_that("the profile's log-likelihood and score are the kernels' at its estimates", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data()
  fit <- do.call(prf_fit, c(list(dt), PRF))
  g <- prf_gradient(fit)
  expect_identical(fit$loglik, -g$objective)
  # the constants' first-order condition holds: each score is small next to
  # its chosen count (the residual rule) or the decrement is (the Newton rule)
  rec <- fit$optimizer$concentration
  expect_true(rec$residual < 1e-10 || isTRUE(rec$decrement < 1e-12))
  expect_lt(max(abs(g$gradient[fit$param_map$asc])), 1e-4)
})

test_that("the outer gradient is the profile's, by the envelope theorem", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data()
  inp <- prepare_mxl_data(dt, id_col = "id", alt_col = "alt",
                          choice_col = "choice", covariate_cols = c("x1", "x2"),
                          random_var_cols = c("w1", "w2"),
                          include_outside_option = TRUE, outside_opt_label = 0L)
  eta <- get_halton_normals(20L, length(inp$M), 2L)
  kernel <- function(th) {
    mxl_delta_hessian_parallel(th, inp$X, inp$W, inp$alt_idx, inp$choice_idx,
                               inp$M, inp$weights, eta, c(0L, 0L), FALSE,
                               FALSE, TRUE, TRUE)
  }
  J <- nrow(inp$alt_mapping) - 1L
  c_idx <- 1:4
  d_idx <- 4L + seq_len(J)
  cc <- .chosen_counts(inp, J)
  solve <- .mxl_profile_solver(kernel, c_idx, d_idx, cc$obs, cc$ref, 1,
                               .concentrate_control(NULL))
  th_c <- c(0.5, -0.4, log(0.8), log(0.6))
  s <- solve(th_c, log(cc$obs / cc$ref), 100L)
  expect_true(s$converged)
  f <- function(tc) solve(tc, s$delta, 25L)$value
  expect_equal(s$gradient[c_idx], numDeriv::grad(f, th_c), tolerance = 1e-6)
})

test_that("the profile's variance is the full Hessian's at its estimates", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_panel()
  dt[, cl := pid %% 20L]
  for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
    fit <- do.call(prf_panel_fit, c(list(dt, se_method = se,
                                         cluster_col = if (se == "cluster") "cl"),
                                    PRF))
    type <- switch(se, sandwich = "robust", se)
    expect_identical(suppressMessages(vcov(fit, type = type)), fit$vcov,
                     label = se)
  }
  skip_on_cran()
  # and the joint fit's at the same maximum, not the continuous block's
  # alone: the Schur complement's inverse is the full inverse's block
  ref <- do.call(prf_panel_fit, c(list(dt, control = PRF_TIGHT), PRF_CERT))
  fit <- do.call(prf_panel_fit, c(list(dt), PRF, PRF_CERT))
  expect_true(all(coef(fit)[c("L_11", "L_22")] > -3))
  expect_equal(fit$se, ref$se, tolerance = 1e-4)
  # treating the constants as known would understate the standard errors:
  # the Schur complement's inverse dominates the continuous block's inverse
  c_idx <- setdiff(seq_along(coef(fit)), fit$param_map$asc)
  H <- compute_hessian(fit)
  se_cc <- sqrt(diag(solve(H[c_idx, c_idx])))
  expect_true(all(fit$se[c_idx] >= se_cc * (1 - 1e-10)))
})

test_that("a polish after the profile certifies it and reuses its curvature", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_panel()
  real_h <- mxl_hessian_parallel
  n_h <- 0L
  local_mocked_bindings(mxl_hessian_parallel = function(...) {
    n_h <<- n_h + 1L
    real_h(...)
  })
  fit <- do.call(prf_panel_fit, c(list(dt, polish = "newton"), PRF))
  p <- fit$optimizer$polish
  expect_true(p$converged)
  expect_identical(n_h, p$evaluations[["curvature"]])
  expect_identical(fit$convergence, 1L)
  fb <- do.call(prf_panel_fit, c(list(dt, polish = "bhhh",
                                      polish_control = list(max_iter = 100L)),
                                 PRF))
  expect_true(fb$optimizer$polish$converged)
  expect_identical(fb$optimizer$concentration$method, "profile_likelihood")
})

# Data whose constants lack a finite estimate (as test-unbounded-constants.R
# builds them): N situations over J alternatives, the chosen one cycling
# through `chosen`
prf_ucon_data <- function(N = 40L, J = 4L, chosen = seq_len(J), seed = 1L) {
  set.seed(seed)
  pick <- chosen[(seq_len(N) - 1L) %% length(chosen) + 1L]
  dt <- data.table::data.table(id = rep(seq_len(N), each = J),
                               alt = rep(seq_len(J), N))
  dt[, choice := as.integer(alt == pick[id])]
  dt[, `:=`(x1 = stats::rnorm(.N), w1 = stats::rnorm(.N))]
  dt[]
}
prf_ucon_fit <- function(dt, ...) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1", S = 10L,
    concentrate_out_asc = "profile_likelihood", ...))
}

test_that("constants without a finite estimate stop the profile before the draws", {
  n_draws <- 0L
  n_delta <- 0L
  real_d <- mxl_delta_hessian_parallel
  local_mocked_bindings(
    get_halton_normals = function(...) {
      n_draws <<- n_draws + 1L
      .halton_cube(...)
    },
    mxl_delta_hessian_parallel = function(...) {
      n_delta <<- n_delta + 1L
      real_d(...)
    })
  never <- prf_ucon_data(chosen = 1:3)
  always <- prf_ucon_data()
  always[id <= 3L, choice := 0L]
  always <- rbind(always, data.table::data.table(
    id = 1:3, alt = 5L, choice = 1L, x1 = c(0.1, -0.2, 0.3), w1 = c(0.4, 0.5, -0.6)))
  data.table::setorderv(always, c("id", "alt"))
  dup <- rbind(always, data.table::data.table(id = 1L, alt = 5L, choice = 0L,
                                             x1 = 0.7, w1 = -0.1))
  data.table::setorderv(dup, c("id", "alt"))
  out_never <- prf_ucon_data()
  out_always <- data.table::copy(out_never)[, choice := 0L]
  base_never <- prf_ucon_data(chosen = 2:4)
  base_always <- prf_ucon_data(chosen = 2:4)
  base_always <- base_always[!(alt == 1L & id > 5L)]
  base_always[id <= 5L, choice := as.integer(alt == 1L)]
  cases <- list(
    never = list(dt = never, msg = "Never chosen in the data"),
    always = list(dt = always, msg = "Chosen in every choice situation that offers them"),
    outside_never = list(dt = out_never, msg = "The outside option is never chosen",
                         include_outside_option = TRUE),
    outside_always = list(dt = out_always, msg = "The outside option is chosen in every",
                          include_outside_option = TRUE),
    base_never = list(dt = base_never, msg = "The base alternative 1"),
    base_always = list(dt = base_always, msg = "is chosen in every choice situation that offers it"),
    # listed twice where it is chosen: alt_mapping counts its rows, 4 of them
    # for 3 choices, which the warning's counts miss
    always_listed_twice = list(dt = dup, msg = "Chosen in every choice situation that offers them"))
  for (nm in names(cases)) {
    a <- cases[[nm]]
    err <- expect_error(do.call(prf_ucon_fit, c(list(a$dt), a[setdiff(names(a), c("dt", "msg"))])),
                        class = "choicer_unbounded_constants", label = nm)
    expect_match(conditionMessage(err), "needs every alternative-specific constant to have a finite",
                 fixed = TRUE, label = nm)
    expect_match(conditionMessage(err), a$msg, fixed = TRUE, label = nm)
    expect_s3_class(err, "error")
  }
  expect_identical(c(n_draws, n_delta), c(0L, 0L))
  # the warning stands without the profile
  expect_warning(suppressMessages(run_mxlogit(
    data = never, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1", S = 10L)),
    class = "choicer_unbounded_constants")
})

test_that("the profile's other inputs are checked before the draws", {
  dt <- prf_data(N = 100L, J = 6L)
  n_draws <- 0L
  local_mocked_bindings(get_halton_normals = function(...) {
    n_draws <<- n_draws + 1L
    .halton_cube(...)
  })
  expect_error(do.call(prf_fit, c(list(dt, use_asc = FALSE), PRF)),
               "so it needs them: use_asc = TRUE", fixed = TRUE)
  expect_error(do.call(prf_fit, c(list(dt, lower = c(ASC_2 = -5)), PRF)),
               "without bounds, but `lower`/`upper` bound ASC_2", fixed = TRUE)
  expect_error(do.call(prf_fit, c(list(dt, lower = c(ASC_3 = 0), upper = c(ASC_3 = 0)),
                                  PRF)),
               "bound ASC_3", fixed = TRUE)
  expect_error(do.call(prf_fit, c(list(dt, concentrate_control = list(max_dense = 3L)),
                                  PRF)),
               "up to concentrate_control$max_dense = 3 constants; this model has 6",
               fixed = TRUE)
  # theta_init, under any coordinates
  for (sv in c("bhhh", "none")) {
    expect_error(do.call(prf_fit, c(list(dt, theta_init = c(0.5, -0.5),
                                          scale_vars = sv), PRF)),
                 "`theta_init` must have one value per parameter (10); got 2.",
                 fixed = TRUE)
    expect_error(do.call(prf_fit, c(list(dt, theta_init = c(NA, rep(0, 9)),
                                          scale_vars = sv), PRF)),
                 "`theta_init` must be finite.", fixed = TRUE)
  }
  expect_identical(n_draws, 0L)
  # the reserved method and the settings, before the data are read
  expect_error(run_mxlogit(data = "not data", concentrate_out_asc = "blp"),
               "\"blp\" is reserved", fixed = TRUE)
  expect_error(run_mxlogit(data = "not data",
                           concentrate_control = list(max_iter = 0L)),
               "`concentrate_control$max_iter` must be a whole number >= 1.",
               fixed = TRUE)
  expect_error(run_mxlogit(data = "not data",
                           concentrate_control = list(tol_residual = -1)),
               "`concentrate_control$tol_residual` must be a positive number.",
               fixed = TRUE)
  expect_error(run_mxlogit(data = "not data", concentrate_control = list(foo = 1)),
               "Unknown `concentrate_control` setting(s): foo", fixed = TRUE)
  expect_error(run_mxlogit(data = "not data", concentrate_control = 5),
               "`concentrate_control` must be a list.", fixed = TRUE)
  expect_identical(.concentrate_control(list(max_iter = 7, tol_decrement = NULL)),
                   list(max_iter = 7L, tol_residual = 1e-10, tol_decrement = 1e-12,
                        max_dense = 10000L))
})

test_that("a warm start comes from the best converged point, never a trial", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data(N = 200L, J = 8L)
  real_solver <- .mxl_profile_solver
  solves <- list()
  local_mocked_bindings(.mxl_profile_solver = function(...) {
    s <- real_solver(...)
    function(theta_c, delta0, max_iter) {
      r <- s(theta_c, delta0, max_iter)
      solves[[length(solves) + 1L]] <<- list(theta_c = theta_c, delta0 = delta0,
                                             delta = r$delta, converged = r$converged,
                                             max_iter = max_iter)
      r
    }
  })
  # x0 (the start: the cache), a far trial, a point a hair uphill of x0 (so
  # that x0 stays the best: a cache by all.equal() would take it for x0), x0
  # again
  path <- NULL
  opt_path <- function(theta_init, eval_f, ...) {
    x0 <- theta_init
    e0 <- eval_f(x0)
    far <- x0 + c(3, -3, 2, 2)
    near <- x0 + 1e-10 * e0$gradient / sqrt(sum(e0$gradient^2))
    f <- c(e0$objective,
           vapply(list(far, near, x0), function(x) eval_f(x)$objective, 0))
    path <<- f
    list(par = x0, value = f[4L], convergence = 0L)
  }
  # the far trial's solve may fail: the optimizer then sees an overflow
  fit <- suppressWarnings(do.call(prf_fit, c(list(dt, optimizer = opt_path), PRF)),
                          classes = "choicer_profile_unconverged")
  # the start's solve, then the far trial and the near point; x0 twice from
  # the cache
  expect_length(solves, 3L)
  # the first solve, from the cold start, may take four times max_iter
  expect_equal(vapply(solves, `[[`, 0, "max_iter"), c(100, 25, 25))
  d0 <- solves[[1L]]$delta
  expect_identical(solves[[2L]]$delta0, d0)
  expect_identical(solves[[3L]]$delta0, d0)
  expect_identical(path[4L], path[1L])
  expect_gt(path[3L], path[1L])
  expect_identical(fit$optimizer$concentration$evaluations, 4L)
  # the constants returned are the start's, bit for bit
  expect_identical(unname(coef(fit)[fit$param_map$asc]), d0)
})

test_that("an inner solve that fails is an overflow to the optimizer", {
  dt <- prf_data(N = 200L, J = 8L)
  real_solver <- .mxl_profile_solver
  n <- 0L
  local_mocked_bindings(.mxl_profile_solver = function(...) {
    s <- real_solver(...)
    function(theta_c, delta0, max_iter) {
      n <<- n + 1L
      r <- s(theta_c, delta0, max_iter)
      if (n == 2L) {
        r$converged <- FALSE
        r$error <- "the kernel ran out of memory"
      }
      r
    }
  })
  seen <- NULL
  opt_two <- function(theta_init, eval_f, ...) {
    f1 <- eval_f(theta_init)$objective
    f2 <- eval_f(theta_init + 0.1)$objective
    seen <<- c(f1, f2)
    list(par = theta_init, value = f1, convergence = 0L)
  }
  # after the best point (the start): the optimizer may have stopped on it
  w <- expect_warning(fit <- do.call(prf_fit, c(list(dt, optimizer = opt_two), PRF)),
                      "1 of the inner solves after the optimizer's best point did not converge",
                      class = "choicer_profile_unconverged")
  expect_match(conditionMessage(w), "the last error: the kernel ran out of memory",
               fixed = TRUE)
  expect_gte(seen[2L], 1e10)
  expect_lt(seen[1L], 1e10)
  expect_identical(fit$optimizer$concentration$failed, 1L)
  expect_match(fit$message, "profile likelihood: 1 of the inner solves")
  # before a better point: the line search's ordinary backtracking, counted
  # without a warning
  n <- 0L
  f0 <- NULL
  opt_back <- function(theta_init, eval_f, ...) {
    e <- eval_f(theta_init)
    f0 <<- e$objective
    eval_f(theta_init + 0.1)
    x <- theta_init - 0.01 * e$gradient / sqrt(sum(e$gradient^2))
    list(par = x, value = eval_f(x)$objective, convergence = 0L)
  }
  expect_no_warning(fit <- do.call(prf_fit, c(list(dt, optimizer = opt_back), PRF)))
  expect_identical(fit$optimizer$concentration$failed, 1L)
  expect_lt(-fit$loglik, f0)
  expect_false(grepl("profile likelihood", fit$message, fixed = TRUE))
  # a first solve that fails stops the fit
  local_mocked_bindings(.mxl_profile_solver = function(...) {
    s <- real_solver(...)
    function(theta_c, delta0, max_iter) {
      r <- s(theta_c, delta0, max_iter)
      r$converged <- FALSE
      r
    }
  })
  expect_error(do.call(prf_fit, c(list(dt), PRF)),
               "the constants could not be solved at the start values",
               class = "choicer_profile_unconverged")
  # an error there (a kernel's memory check, say) is reported as itself
  local_mocked_bindings(mxl_delta_hessian_parallel = function(...) {
    stop("The Hessian of the constants needs about 9 GB")
  })
  err <- expect_error(do.call(prf_fit, c(list(dt), PRF)),
                      "stopped with an error: The Hessian of the constants needs about 9 GB",
                      fixed = TRUE)
  expect_false(inherits(err, "choicer_profile_unconverged"))
})

test_that("the constants are solved at the point the optimizer returns", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data(N = 200L, J = 8L)
  # a custom optimizer that evaluates its start (the cache) and returns a
  # point downhill that it never evaluated: the constants are solved there
  opt_away <- function(theta_init, eval_f, ...) {
    e <- eval_f(theta_init)
    list(par = theta_init - 0.05 * e$gradient / sqrt(sum(e$gradient^2)),
         value = e$objective, convergence = 0L)
  }
  fit <- do.call(prf_fit, c(list(dt, optimizer = opt_away), PRF))
  rec <- fit$optimizer$concentration
  expect_identical(rec$evaluations, 1L)
  expect_identical(rec$solves, 2L)  # at the start and at the point returned
  g <- prf_gradient(fit)
  expect_identical(fit$loglik, -g$objective)
  expect_true(rec$residual < 1e-10 || isTRUE(rec$decrement < 1e-12))
  expect_lt(max(abs(g$gradient[fit$param_map$asc])), 1e-4)
  # when that solve fails, the best point solved (here the start) is returned
  real_solver <- .mxl_profile_solver
  n <- 0L
  local_mocked_bindings(.mxl_profile_solver = function(...) {
    s <- real_solver(...)
    function(theta_c, delta0, max_iter) {
      n <<- n + 1L
      r <- s(theta_c, delta0, max_iter)
      if (n == 1L) {
        s0 <<- r
        t0 <<- theta_c
      }
      if (n == 2L) r$converged <- FALSE
      r
    }
  })
  s0 <- NULL
  t0 <- NULL
  expect_warning(f2 <- do.call(prf_fit, c(list(dt, optimizer = opt_away), PRF)),
                 "could not be solved at the optimizer's solution, so the best point",
                 class = "choicer_profile_unconverged")
  expect_identical(n, 2L)
  expect_identical(unname(coef(f2)[f2$param_map$asc]), s0$delta)
  expect_identical(unname(coef(f2)[-f2$param_map$asc]), unname(t0))
  expect_identical(f2$loglik, -s0$value)
  expect_identical(f2$loglik, -prf_gradient(f2)$objective)
  expect_identical(f2$optimizer$concentration$failed, 0L)
  expect_match(f2$message, "profile likelihood: the constants could not be solved")
})

test_that("the constants' counts are weighted, and the residual rule scale-free", {
  # .chosen_counts(): the summed weights of the situations choosing each
  # free constant's alternative, and the reference's
  inp <- list(M = c(3L, 2L, 3L, 2L), choice_idx = c(2L, 0L, 3L, 1L),
              alt_idx = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 3L, 1L, 3L),
              weights = c(0.5, 2, 1.5, 4), include_outside_option = TRUE)
  cc <- .chosen_counts(inp, 3L)
  expect_identical(cc$obs, c(4, 0.5, 1.5))
  expect_identical(cc$ref, 2)
  inp$include_outside_option <- FALSE
  inp$choice_idx <- c(2L, 1L, 3L, 1L)
  cc <- .chosen_counts(inp, 2L)
  expect_identical(cc$obs, c(2.5, 1.5))
  expect_identical(cc$ref, 4)
  # the solve stops at the start when every score, and the reference's (minus
  # their sum), is below tol_residual relative to its count plus the mean
  # weight, whatever the weights' scale, and goes on otherwise; a quadratic
  # objective in delta (H = diag(h)) whose curvature the stub reports as
  # h / k (k = 1: exact)
  quad <- function(delta_star, h, k = 1, nonfinite = FALSE) {
    function(th) {
      dl <- th[3:4]
      g <- h * (dl - delta_star)
      list(objective = 0.5 * sum(h * (dl - delta_star)^2), gradient = c(0, 0, g),
           overflow = FALSE, nonfinite = nonfinite, hessian = diag(h / k))
    }
  }
  ctl <- .concentrate_control(list(tol_decrement = 1e-300))
  for (scale in c(1, 1e-6, 1e6)) {
    obs <- c(10, 30) * scale
    h <- c(5, 20) * scale
    s <- .mxl_profile_solver(quad(c(1, 2), h), 1:2, 3:4, obs, 1e3 * scale,
                             scale, ctl)
    # scores 0.5e-10 (obs + w_bar), just inside the rule
    eps <- 0.5e-10 * (obs + scale) / h
    r <- s(c(0, 0), c(1, 2) + eps, 25L)
    expect_true(r$converged, label = scale)
    expect_identical(r$iterations, 0L, label = scale)
    expect_identical(r$message, "every constant's score is below tol_residual")
    expect_identical(r$decrement, NA_real_)
    # twice the rule's threshold: a Newton step, then the rule, which leaves
    # no decrement at the point it stops
    r <- s(c(0, 0), c(1, 2) + 4 * eps, 25L)
    expect_identical(r$iterations, 1L, label = scale)
    expect_identical(r$message, "every constant's score is below tol_residual")
    expect_identical(r$decrement, NA_real_)
    expect_lt(r$residual, 1e-10)
  }
  eps <- 0.5e-10 * (c(10, 30) + 1) / c(5, 20)
  # the same scores against a reference chosen once: its score, minus their
  # sum, is outside the rule
  s <- .mxl_profile_solver(quad(c(1, 2), c(5, 20)), 1:2, 3:4, c(10, 30), 1, 1, ctl)
  expect_gt(s(c(0, 0), c(1, 2) + eps, 25L)$iterations, 0L)
  # a constant whose alternative was chosen at a weight of 0.1, the mean
  # weight 1: inside the rule, which adds the mean weight to the count
  eps <- 0.5e-10 * (c(0.1, 30) + 1) / c(5, 20)
  s <- .mxl_profile_solver(quad(c(1, 2), c(5, 20)), 1:2, 3:4, c(0.1, 30), 1e3, 1, ctl)
  expect_identical(s(c(0, 0), c(1, 2) + eps, 25L)$iterations, 0L)
  # steps nine times too long take four trials each: the budget of three
  # passes a step, not one, lets 20 steps converge, and 10 steps' budget
  # (31 passes) runs out
  s <- .mxl_profile_solver(quad(c(1, 2), c(5, 20), k = 9), 1:2, 3:4, c(10, 30),
                           1e3, 1, ctl)
  r <- s(c(0, 0), c(3, -1), 20L)
  expect_true(r$converged)
  expect_gt(r$passes, 21L)
  r <- s(c(0, 0), c(3, -1), 10L)
  expect_false(r$converged)
  expect_identical(r$message, "the solve used its 31 passes")
  expect_null(r$error)
  # the decrement rule over the mean weight: steps nine times too long shrink
  # the decrement 64-fold per step, so a rule in raw units would stop the
  # run whose weights are scaled by 2^-30 (objective, gradient, curvature,
  # counts and mean weight alike) about five steps early
  ctl_d <- .concentrate_control(list(tol_residual = 1e-300))
  run_at <- function(c) {
    s <- .mxl_profile_solver(quad(c(1, 2), c(5, 20) * c, k = 9), 1:2, 3:4,
                             c(10, 30) * c, 1e3 * c, c, ctl_d)
    s(c(0, 0), c(3, -1), 50L)
  }
  r1 <- run_at(1)
  r2 <- run_at(2^-30)
  expect_true(r1$converged)
  expect_identical(r1$message, "the Newton decrement is below tol_decrement")
  expect_identical(r2[c("iterations", "passes", "delta", "message")],
                   r1[c("iterations", "passes", "delta", "message")])
  # a point whose gradient has entries that are not finite fails at once
  s <- .mxl_profile_solver(quad(c(1, 2), c(5, 20), nonfinite = TRUE), 1:2, 3:4,
                           c(10, 30), 1e3, 1, ctl)
  r <- s(c(0, 0), c(3, -1), 20L)
  expect_false(r$converged)
  expect_identical(r$passes, 1L)
})

test_that("the inner solve's stopping rules do not depend on the weights' scale", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data(N = 200L, J = 8L)
  inp <- prepare_mxl_data(dt, id_col = "id", alt_col = "alt",
                          choice_col = "choice", covariate_cols = c("x1", "x2"),
                          random_var_cols = c("w1", "w2"),
                          include_outside_option = TRUE, outside_opt_label = 0L)
  eta <- get_halton_normals(20L, length(inp$M), 2L)
  set.seed(36)
  w <- stats::runif(length(inp$M), 0.5, 2)
  J <- nrow(inp$alt_mapping) - 1L
  # weights scaled by 2^-40 scale the objective, gradient and block exactly
  # (2^-20 would not do: the solve's decrements happen to skip the window a
  # rule in raw units would open)
  solve_at <- function(wt) {
    d <- inp
    d$weights <- wt
    kernel <- function(th) {
      mxl_delta_hessian_parallel(th, d$X, d$W, d$alt_idx, d$choice_idx, d$M,
                                 d$weights, eta, c(0L, 0L), FALSE, FALSE, TRUE,
                                 TRUE)
    }
    cc <- .chosen_counts(d, J)
    s <- .mxl_profile_solver(kernel, 1:4, 4L + seq_len(J), cc$obs, cc$ref,
                             mean(wt), .concentrate_control(NULL))
    s(c(0.5, -0.4, log(0.8), log(0.6)), log(cc$obs / cc$ref), 100L)
  }
  r1 <- solve_at(w)
  r2 <- solve_at(w * 2^-40)
  expect_true(r1$converged)
  expect_identical(r2$delta, r1$delta)
  expect_identical(r2[c("iterations", "passes", "residual", "decrement",
                        "message")],
                   r1[c("iterations", "passes", "residual", "decrement",
                        "message")])
  expect_identical(r2$value, r1$value * 2^-40)
})

test_that("concentrate_out_asc = \"none\" is the fit without it, bit for bit", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- prf_data(N = 200L, J = 8L)
  strip <- function(fit) {
    fit$call <- NULL
    fit$optimizer$elapsed_time <- NULL
    fit
  }
  f0 <- prf_fit(dt, control = list(maxeval = 20L))
  f1 <- prf_fit(dt, control = list(maxeval = 20L), concentrate_out_asc = "none",
                concentrate_control = list(max_iter = 3L))
  expect_true(identical(strip(f1), strip(f0), num.eq = FALSE, single.NA = FALSE))
  expect_named(f0$optimizer, c("name", "control", "elapsed_time", "iterations"))
})

test_that("print() and summary() show a profile fit", {
  dt <- prf_data(N = 200L, J = 8L)
  fit <- do.call(prf_fit, c(list(dt), PRF))
  expect_output(print(fit), "Convergence")
  expect_output(print(summary(fit)), "ASC_")
})
