# run_mxlogit(use_asc = FALSE): the mixed logit without alternative-specific
# constants. Its parameters are the fixed coefficients, the random
# coefficients' means (rc_mean) and their Cholesky factor, and it has no
# `asc` block, as in run_mnlogit() and run_nestlogit(). It is the
# use_asc = TRUE model with every constant held at 0, which these tests use
# as the reference: at the same point the two give the same likelihood,
# information, predictions and post-estimation results.

# Data from simulate_mxl_data() without constants in the DGP: x1, x2 fixed,
# w1, w2 random. `outside_option = TRUE` adds an alt = 0 row of zero
# covariates per situation.
uasc_sim <- function(seed, outside_option = FALSE, T = 1L, N = 200L) {
  simulate_mxl_data(N = N, J = 3L, T = T, beta = c(0.8, -0.6),
                    Sigma = diag(c(1, 0.5)), delta = rep(0, 3L),
                    outside_option = outside_option, seed = seed)
}

uasc_fit <- function(data, ..., use_asc = FALSE, quiet = TRUE) {
  fit <- function() {
    run_mxlogit(
      data = data, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      S = 20L, use_asc = use_asc, ...
    )
  }
  res <- if (quiet) suppressMessages(fit()) else fit()
  # the variance computed at fit time (one that fails only warns)
  expect_null(res$variance_error)
  res
}

# A custom optimizer that returns its start: the fit is held there.
uasc_hold <- function(theta_init, eval_f, ...) {
  list(par = theta_init, value = eval_f(theta_init)$objective,
       convergence = 0L)
}

# The use_asc = TRUE model held at a no-constant fit's estimates, every
# constant at 0, on the same data and draws (`...` repeats the fit's
# arguments).
uasc_held <- function(fit0, data, ...) {
  n_asc <- nrow(fit0$alt_mapping) - 1L
  uasc_fit(data, ..., use_asc = TRUE, scale_vars = "none",
           theta_init = c(unname(coef(fit0)), rep(0, n_asc)),
           optimizer = uasc_hold)
}

test_that("run_mxlogit(use_asc = FALSE) fits a model without constants", {
  # Releases since 0.2.0 stopped here ("Theta vector too long: 6 parameters
  # given but the model expects 4": J - 1 constants were counted whatever
  # use_asc said); 0.1.0 returned them, held at 0, with every SE NA
  dt <- uasc_sim(21L)$data
  fits <- lapply(c(bhhh = "bhhh", none = "none", sd = "sd"),
                 function(sv) uasc_fit(dt, scale_vars = sv))
  for (fit in fits) {
    expect_named(coef(fit), c("x1", "x2", "L_11", "L_22"))
    expect_equal(fit$n_params, 4)
    expect_null(fit$param_map$asc)
    expect_false(fit$use_asc)
    expect_equal(dim(fit$vcov), c(4L, 4L))
    expect_named(fit$se, names(coef(fit)))
    expect_true(all(is.finite(fit$se)))
    expect_false(any(grepl("ASC", rownames(summary(fit)$coefficients))))
  }
  # one model: the three coordinates reach the same maximum
  expect_equal(fits$none$loglik, fits$bhhh$loglik, tolerance = 1e-8)
  expect_equal(fits$sd$loglik, fits$bhhh$loglik, tolerance = 1e-8)
  expect_equal(coef(fits$none), coef(fits$bhhh), tolerance = 1e-4)

  # With the outside option, alt_mapping's first row is the outside good:
  # still no constant
  sim <- uasc_sim(22L, outside_option = TRUE)
  fit <- uasc_fit(sim$data[alt != 0L], include_outside_option = TRUE,
                  outside_opt_label = 0L)
  expect_equal(nrow(fit$alt_mapping), 4L)
  expect_named(coef(fit), c("x1", "x2", "L_11", "L_22"))
  expect_null(fit$param_map$asc)
  # ... the same likelihood as the outside good as a base row of zero
  # covariates, whose utility is 0 either way
  base_row <- uasc_fit(sim$data, outside_opt_label = 0L)
  expect_equal(fit$loglik, base_row$loglik, tolerance = 1e-8)

  # Without constants no parameter is named ASC_<label>, so a covariate may be
  # (with them, the names clash)
  dt2 <- data.table::copy(dt)
  data.table::setnames(dt2, "x2", "ASC_2")
  args <- list(data = dt2, id_col = "id", alt_col = "alt",
               choice_col = "choice", covariate_cols = c("x1", "ASC_2"),
               random_var_cols = c("w1", "w2"), S = 20L)
  fit <- suppressMessages(do.call(run_mxlogit, c(args, use_asc = FALSE)))
  expect_named(coef(fit), c("x1", "ASC_2", "L_11", "L_22"))
  expect_equal(unname(coef(fit)), unname(coef(fits$bhhh)))
  expect_error(do.call(run_mxlogit, c(args, use_asc = TRUE)),
               "Parameter names must be unique; repeated: 'ASC_2'")
})

test_that("start values and bounds have one entry per parameter", {
  dt <- uasc_sim(23L, N = 60L)$data
  W <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                        c("w1", "w2"))$W
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- list(theta_init = theta_init, lower = lower, upper = upper)
    stop("captured")
  })
  start <- function(...) {
    expect_error(uasc_fit(dt, scale_vars = "none", ...), "captured")
    seen
  }
  # x1, x2, L_11, L_22: the default start in the columns' units
  s <- start()
  expect_equal(s$theta_init, c(0, 0, -log(choicer:::design_column_step(W))),
               tolerance = 1e-14)
  expect_identical(s$lower, rep(-Inf, 4L))
  expect_identical(s$upper, rep(Inf, 4L))
  # the use_asc = TRUE start without its constants
  expect_error(uasc_fit(dt, scale_vars = "none", use_asc = TRUE), "captured")
  expect_identical(s$theta_init, seen$theta_init[1:4])
  # named bounds address the same layout, and move the default start
  # (L_22 starts at -log(sd(w2)), about 0.55 for a uniform w2 on (-1, 1))
  suppressMessages(expect_message(
    s <- start(lower = c(L_22 = 1), upper = c(x2 = 1), quiet = FALSE),
    "default start of L_22 lies outside"
  ))
  expect_identical(s$lower, c(-Inf, -Inf, -Inf, 1))
  expect_identical(s$upper, c(Inf, 1, Inf, Inf))
  expect_identical(s$theta_init[4], 1)
  # a supplied start is used as given
  expect_identical(start(theta_init = c(0.1, -0.1, 0, 0))$theta_init,
                   c(0.1, -0.1, 0, 0))
  # the constants' names and slots do not exist
  expect_error(uasc_fit(dt, lower = c(ASC_2 = 0)),
               "Unknown parameter name\\(s\\) in `lower`: ASC_2")
  expect_error(uasc_fit(dt, upper = rep(1, 6L)),
               "must have length n_params \\(4\\); got 6")
  expect_error(uasc_fit(dt, theta_init = rep(0, 6L)),
               "one value per parameter \\(4\\); got 6")
})

test_that("a fit without constants is the model with its constants at 0", {
  sim <- uasc_sim(24L, outside_option = TRUE)
  inside <- sim$data[alt != 0L]
  designs <- list(
    cross = list(data = uasc_sim(25L)$data),
    outside = list(data = inside, include_outside_option = TRUE,
                   outside_opt_label = 0L),
    base_row = list(data = sim$data, outside_opt_label = 0L),
    panel = list(data = uasc_sim(26L, T = 3L, N = 80L)$data,
                 person_col = "pid"),
    generate = list(data = inside, include_outside_option = TRUE,
                    outside_opt_label = 0L, draws = "generate", seed = 7L)
  )
  for (nm in names(designs)) {
    a <- designs[[nm]]
    fit0 <- do.call(uasc_fit, a)
    fitA <- do.call(uasc_held, c(list(fit0), a))
    idx <- seq_len(fit0$n_params)
    expect_true(all(startsWith(names(coef(fitA))[-idx], "ASC_")), label = nm)
    expect_equal(fit0$loglik, fitA$loglik, tolerance = 1e-12, label = nm)
    # Each unit's score and the Hessian are the full model's without the
    # constants' entries, so every variance is built from the full model's
    # blocks (compared before inversion: a small sample can leave a variance
    # on its zero plateau, where the information is singular)
    expect_equal(choicer:::compute_scores(fit0),
                 choicer:::compute_scores(fitA)[, idx, drop = FALSE],
                 tolerance = 1e-12, ignore_attr = TRUE, label = nm)
    expect_equal(choicer:::.compute_bread(fit0),
                 choicer:::.compute_bread(fitA)[idx, idx],
                 tolerance = 1e-12, ignore_attr = TRUE, label = nm)
    # Prediction and post-estimation
    expect_equal(predict(fit0), predict(fitA), tolerance = 1e-12, label = nm)
    nd <- data.table::copy(a$data)
    nd[, x1 := x1 + 0.5]
    cf0 <- suppressMessages(predict(fit0, newdata = nd, type = "shares"))
    expect_equal(cf0, suppressMessages(predict(fitA, newdata = nd,
                                               type = "shares")),
                 tolerance = 1e-12, label = nm)
    expect_equal(elasticities(fit0, "x1"), elasticities(fitA, "x1"),
                 tolerance = 1e-12, label = nm)
    expect_equal(elasticities(fit0, "w1", is_random_coef = TRUE),
                 elasticities(fitA, "w1", is_random_coef = TRUE),
                 tolerance = 1e-12, label = nm)
    expect_equal(diversion_ratios(fit0, "w2", is_random_coef = TRUE),
                 diversion_ratios(fitA, "w2", is_random_coef = TRUE),
                 tolerance = 1e-12, label = nm)
    expect_equal(logsum(fit0), logsum(fitA), tolerance = 1e-12, label = nm)
    expect_equal(conditional_tastes(fit0)[c("mean", "sd")],
                 conditional_tastes(fitA)[c("mean", "sd")],
                 tolerance = 1e-12, label = nm)
    # The constants that would give a counterfactual's shares, from zeros
    expect_equal(blp(fit0, target_shares = cf0),
                 blp(fitA, target_shares = cf0), tolerance = 1e-10,
                 label = nm)
  }
})

test_that("post-estimation runs on a fit without constants", {
  sim <- uasc_sim(27L, outside_option = TRUE)
  fit <- uasc_fit(sim$data[alt != 0L], include_outside_option = TRUE,
                  outside_opt_label = 0L, rc_mean = TRUE)
  p <- fit$n_params
  expect_named(coef(fit), c("x1", "x2", "Mu_w1", "Mu_w2", "L_11", "L_22"))
  cl <- stats::setNames(rep(seq_len(40L), length.out = fit$nobs),
                        fit$data$situation_ids)
  for (v in list(vcov(fit, type = "bhhh"), vcov(fit, type = "robust"),
                 vcov(fit, type = "cluster", cluster = cl))) {
    expect_equal(dim(v), c(p, p))
    expect_true(all(is.finite(diag(v)) & diag(v) > 0))
  }
  expect_length(suppressWarnings(wesml_vcov(fit, type = "se")), p)
  w <- wtp(fit, price_var = "x2")
  expect_identical(rownames(w), c("x1", "Mu_w1", "Mu_w2"))
  expect_error(wtp(fit, price_var = "x2", attr_vars = "ASC_1"),
               "not found")
  expect_true(is.finite(consumer_surplus(fit, price_var = "x2")$mean_cs))
  expect_true(is.finite(gof(fit)$loglik))
  rec <- recovery_table(fit, sim)
  expect_setequal(unique(rec$group), c("beta", "mu", "sigma"))
})
