# blp() on fits: the checks it makes in R before the solver, the start, the
# method and its settings, and the result's row names and "convergence"
# attribute (R/blp.R: .blp_common()).

bm_mnl_fit <- function(outside = FALSE, seed = 7L) {
  set.seed(seed)
  N <- 160L
  labs <- c("bus", "car", "train", "walk")
  dt <- data.table::data.table(id = rep(seq_len(N), each = 4L), alt = rep(labs, N),
                               x1 = stats::rnorm(4L * N), x2 = stats::runif(4L * N))
  # with an outside option, a situation with no inside choice chose it
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N + outside, 1L)), by = id]
  if (outside) {
    suppressWarnings(run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                                 include_outside_option = TRUE))
  } else {
    run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
  }
}

# The reproducer of the review of PR #56 as a fit of each model: three choice
# situations (each repeated) in which alternatives 2, 3 and 4 each dominate
# one, with beta = 1, the nested logit's lambdas 1 and the mixed logit's
# random design zero, so that each is the same multinomial logit.
bm_repro_x <- c(0, 3, -3, -3, 0, -3, 3, -3, 0, -3, -3, 3)
bm_repro_fit <- function(model) {
  R <- 10L
  set.seed(5L)
  dt <- data.table::data.table(id = rep(seq_len(3L * R), each = 4L),
                               alt = rep(1:4, 3L * R), x1 = rep(bm_repro_x, R),
                               nest = rep(c(1L, 1L, 2L, 2L), 3L * R),
                               w1 = stats::rnorm(12L * R))
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  fit <- suppressWarnings(switch(model,
    mnl = run_mnlogit(dt, "id", "alt", "choice", "x1"),
    nl = run_nestlogit(dt, "id", "alt", "choice", "x1", "nest",
                       control = list(maxeval = 20L)),
    mxl = run_mxlogit(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
                      covariate_cols = "x1", random_var_cols = "w1", S = 5L,
                      control = list(maxeval = 20L))))
  pm <- fit$param_map
  fit$coefficients[pm$beta] <- 1
  if (model == "nl") fit$coefficients[pm$lambda] <- 1
  if (model == "mxl") fit$data$W[] <- 0
  fit
}

test_that("every method solves the review's reproducer through every model's blp()", {
  # The diagonal step alone diverged here (test-blp-engine.R)
  ds <- c(0, 1, 1, 1)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, ds[-1L]), matrix(bm_repro_x),
                                               rep(1:4, 3L), rep(4L, 3L), rep(1, 3L)))
  for (model in c("mnl", "nl", "mxl")) {
    fit <- bm_repro_fit(model)
    # the fit is the reproducer: its solution meets the targets at once
    at <- blp(fit, target_shares = S, delta_init = ds)
    expect_identical(attr(at, "convergence")$evaluations, 1L, label = model)
    for (method in c("simple", "original", "anderson", "squarem")) {
      lab <- sprintf("%s, %s", model, method)
      d <- blp(fit, target_shares = S, delta_init = rep(0, 4L), method = method)
      cv <- attr(d, "convergence")
      expect_true(cv$converged, label = lab)
      expect_identical(cv$method, method, label = lab)
      expect_lte(cv$evaluations, 10L, label = lab)
      expect_equal(as.numeric(d), ds, tolerance = 1e-8, label = lab)
    }
  }
})

test_that("every method recovers known constants on each model", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dn <- create_small_nl_data()
  nl <- run_nestlogit(dn, "id", "alt", "choice", c("x1", "x2"), "nest",
                      control = list(maxeval = 50L))
  nl$coefficients[nl$param_map$lambda] <- c(0.7, 0.85)
  mxl <- suppressWarnings(run_mxlogit(
    data = create_small_mxl_data(), id_col = "id", alt_col = "alt",
    choice_col = "choice", covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    S = 10L, control = list(maxeval = 30L)))
  fits <- list(mnl = bm_mnl_fit(), mnl_outside = bm_mnl_fit(outside = TRUE), nl = nl,
               mxl = mxl)
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    asc <- fit$param_map$asc
    outside <- isTRUE(fit$include_outside_option)
    set.seed(match(nm, names(fits)))
    fit$coefficients[asc] <- stats::rnorm(length(asc))
    ds <- if (outside) fit$coefficients[asc] else c(0, fit$coefficients[asc])
    S <- as.numeric(predict(fit, type = "shares"))
    for (method in c("simple", "original", "anderson", "squarem")) {
      lab <- sprintf("%s, %s", nm, method)
      d <- blp(fit, target_shares = S, delta_init = rep(0, length(ds)), tol = 1e-10,
               method = method)
      expect_true(attr(d, "convergence")$converged, label = lab)
      expect_identical(attr(d, "convergence")$method, method, label = lab)
      expect_equal(as.numeric(d), as.numeric(ds), tolerance = 1e-8, label = lab)
    }
    # a method's settings, through control
    d <- blp(fit, target_shares = S, delta_init = rep(0, length(ds)), tol = 1e-10,
             method = "anderson", control = list(memory = 3L, max_condition = 1e6))
    expect_equal(as.numeric(d), as.numeric(ds), tolerance = 1e-8, label = nm)
  }
})

test_that("method = \"original\" is Berry, Levinsohn and Pakes's contraction", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # (identical results at one thread)
  for (outside in c(FALSE, TRUE)) {
    fit <- bm_mnl_fit(outside)
    # (far from the fitted shares, about 1/4 each)
    S <- if (outside) c(0.3, 0.25, 0.2, 0.15, 0.1) else c(0.4, 0.3, 0.2, 0.1)
    d <- blp(fit, target_shares = S, method = "original", tol = 1e-8)
    pm <- fit$param_map
    start <- if (outside) fit$coefficients[pm$asc] else c(0, fit$coefficients[pm$asc])
    # the engine's plain map: delta + log S - log s(delta), normalized (with
    # blp()'s pass)
    r <- choicer:::mnl_blp_solve(as.numeric(start), S, fit$data$X,
                                 fit$coefficients[pm$beta], fit$data$alt_idx, fit$data$M,
                                 fit$data$weights, outside,
                                 choicer:::.blp_control(map = "plain", patience = 0L,
                                                        pass = "exp"))
    expect_identical(as.numeric(d), as.numeric(r$delta))
    expect_identical(attr(d, "convergence")$evaluations, r$evaluations)
    # one step from the start, by hand: two evaluations return the better of
    # the start and the point one step on
    resid <- function(dl) {
      s <- as.numeric(choicer:::mnl_predict_shares(
        c(fit$coefficients[pm$beta], if (outside) dl else dl[-1L]), fit$data$X,
        fit$data$alt_idx, fit$data$M, fit$data$weights, TRUE, outside))
      r <- log(S) - log(s)
      list(step = if (outside) r[-1L] else r - r[1L],
           res = max(abs(if (outside) r else r - r[1L])))
    }
    start <- as.numeric(start)
    on1 <- start + resid(start)$step
    best <- if (resid(on1)$res < resid(start)$res) on1 else start
    d2 <- suppressWarnings(blp(fit, target_shares = S, method = "original", max_iter = 2))
    expect_equal(as.numeric(d2), best, tolerance = 1e-12)
  }
  # a nested logit's damping multiplies Berry's step too
  dt <- create_small_nl_data()
  nl <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                      control = list(maxeval = 50L))
  nl$coefficients[nl$param_map$lambda] <- c(0.7, 0.85)
  S <- rev(as.numeric(predict(nl, type = "shares")))
  o1 <- blp(nl, target_shares = S, method = "original")
  o5 <- blp(nl, target_shares = S, method = "original", damping = 0.5)
  expect_equal(as.numeric(o5), as.numeric(o1), tolerance = 1e-6)
  expect_gt(attr(o5, "convergence")$evaluations, attr(o1, "convergence")$evaluations)
  expect_identical(attr(o5, "convergence")$method, "original")
})

test_that("blp() labels its constants and records its convergence", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # (identical results at one thread)
  for (outside in c(FALSE, TRUE)) {
    fit <- bm_mnl_fit(outside)
    am <- fit$alt_mapping
    inside <- am$alt_int > 0L
    labs <- as.character(am[[2L]][inside])[order(am$alt_int[inside])]
    S <- as.numeric(predict(fit, type = "shares"))
    d <- blp(fit, target_shares = rev(S))
    expect_true(is.matrix(d))
    expect_identical(dim(d), c(length(labs), 1L))
    expect_identical(rownames(d), labs)
    cv <- attr(d, "convergence")
    expect_identical(names(cv), c("converged", "status", "iterations", "evaluations",
                                  "residual", "tol", "method"))
    expect_true(cv$converged)
    expect_identical(cv$status, "converged")
    expect_identical(cv$method, "anderson")
    expect_identical(cv$tol, 1e-12)
    expect_gte(cv$evaluations, 1L)
    expect_lte(cv$residual, 1e-12)
    # the same constants as the export's
    pm <- fit$param_map
    start <- if (outside) fit$coefficients[pm$asc] else c(0, fit$coefficients[pm$asc])
    ex <- blp_contraction(as.numeric(start), rev(S), fit$data$X,
                          fit$coefficients[pm$beta], fit$data$alt_idx, fit$data$M,
                          fit$data$weights, outside)
    expect_identical(as.numeric(d), as.numeric(ex))
  }
})

test_that("blp() that stops early warns with the record it attaches", {
  fit <- bm_mnl_fit()
  S <- as.numeric(predict(fit, type = "shares"))
  w <- expect_warning(d <- blp(fit, target_shares = rev(S), tol = 1e-14, max_iter = 2),
                      class = "choicer_blp_unconverged")
  cv <- attr(d, "convergence")
  expect_identical(w$convergence, cv)
  expect_false(cv$converged)
  expect_identical(cv$status, "max_iter")
  expect_identical(cv$evaluations, 2L)
  # max_iter = 0: the fitted constants come back, unevaluated
  w0 <- expect_warning(d0 <- blp(fit, target_shares = rev(S), max_iter = 0),
                       class = "choicer_blp_unconverged")
  expect_identical(as.numeric(d0), c(0, as.numeric(fit$coefficients[fit$param_map$asc])))
  expect_identical(attr(d0, "convergence")$evaluations, 0L)
  expect_true(is.na(attr(d0, "convergence")$residual))
})

test_that("blp() checks its targets before the solver", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # (identical results at one thread)
  fit <- bm_mnl_fit()
  S <- as.numeric(predict(fit, type = "shares"))
  expect_error(blp(fit, target_shares = S[-1]), "has 3 entries; it must have 4")
  expect_error(blp(fit, target_shares = replace(S, 2L, NA)),
               "must be finite and positive \\(entry 2")
  expect_error(blp(fit, target_shares = replace(S, 3L, 0)),
               "must be finite and positive \\(entry 3")
  expect_error(blp(fit, target_shares = S * 1.01), "must sum to 1; they sum to")
  expect_error(blp(fit, target_shares = cbind(S, S)), "one-column matrix")
  expect_error(blp(fit, target_shares = as.character(S)), "must be a numeric vector")
  expect_error(blp(fit, target_shares = as.list(S)), "must be a numeric vector")
  expect_silent(d <- blp(fit, target_shares = matrix(S, ncol = 1L)))
  expect_identical(d, blp(fit, target_shares = S))
  fit_oo <- bm_mnl_fit(outside = TRUE)
  So <- as.numeric(predict(fit_oo, type = "shares"))
  expect_error(blp(fit_oo, target_shares = So[-1]),
               "the outside option's share, then one per inside alternative")
  expect_error(blp(fit_oo, target_shares = So * 0.9),
               "must sum to 1, the outside option's share included")
  # the engine's own checks still apply: here, a target the situations
  # offering the alternative cannot reach
  set.seed(3L)
  dt <- data.table::data.table(id = rep(1:60, each = 3L), alt = rep(1:3, 60L),
                               x1 = stats::rnorm(180L))
  dt <- dt[!(alt == 3L & id > 6L)]  # alternative 3 in the first 6 situations only
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  fit3 <- suppressWarnings(run_mnlogit(dt, "id", "alt", "choice", "x1"))
  expect_error(blp(fit3, target_shares = c(0.4, 0.3, 0.3)),
               "not below the weight of the choice situations that offer it")
})

test_that("blp() checks its start, method and settings", {
  fit <- bm_mnl_fit()
  S <- as.numeric(predict(fit, type = "shares"))
  expect_error(blp(fit, target_shares = S, delta_init = rep(0, 3L)),
               "must be a numeric vector of 4 constants")
  expect_error(blp(fit, target_shares = S, delta_init = c(0, NA, 0, 0)),
               "`delta_init` must be finite")
  d <- blp(fit, target_shares = rev(S), delta_init = rep(0, 4L))
  expect_true(attr(d, "convergence")$converged)
  expect_error(blp(fit, target_shares = S, method = "newton"),
               "must be \"anderson\" (the default: Anderson acceleration", fixed = TRUE)
  expect_error(blp(fit, target_shares = S, method = "simple", control = list(m = 5)),
               "`control` takes no settings with method \"simple\" \\(it has `m`\\)")
  expect_error(blp(fit, target_shares = S, control = list(m = 5)),
               "`control` takes `memory` and `max_condition` with method \"anderson\"")
  expect_error(blp(fit, target_shares = S, method = "original", control = list(m = 5)),
               "`control` takes no settings with method \"original\"")
  expect_error(blp(fit, target_shares = S, method = "anderson", control = list(m = 5)),
               "`control` takes `memory` and `max_condition` with method \"anderson\"")
  expect_error(blp(fit, target_shares = S, method = "anderson",
                   control = list(memory = 0)), "`control\\$memory` must be a whole")
  expect_error(blp(fit, target_shares = S, method = "squarem",
                   control = list(scheme = 4)), "`control\\$scheme` must be 1, 2 or 3")
  expect_error(blp(fit, target_shares = S, control = "x"), "must be a list")
  expect_error(blp(fit, target_shares = S, tol = -1), "`tol` must be a single")
  expect_error(blp(fit, target_shares = S, max_iter = 2.5), "`max_iter` must be a whole")
})

test_that("blp() on a nested logit labels, damps and checks damping", {
  dt <- create_small_nl_data()
  fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                       control = list(maxeval = 50L))
  S <- as.numeric(predict(fit, type = "shares"))
  # (the step alone, whose halved steps cost evaluations)
  d1 <- blp(fit, target_shares = rev(S), method = "simple")
  d2 <- blp(fit, target_shares = rev(S), damping = 0.5, method = "simple")
  expect_identical(rownames(d1), as.character(1:6))
  expect_equal(as.numeric(d2), as.numeric(d1), tolerance = 1e-6)
  expect_gt(attr(d2, "convergence")$evaluations, attr(d1, "convergence")$evaluations)
  expect_error(blp(fit, target_shares = S, damping = 0), "`damping` must be a single")
})

test_that("blp() and the exports evaluate the shares from kept exponentials", {
  # the pass and sums blp() and the exports hand the models' solvers, and
  # whether their evaluations read kept exponentials (all but the mixed
  # logit's export, which keeps nothing); the internal callers' defaults stay
  # the usual pass and plain sums
  expect_identical(choicer:::.blp_control()[c("pass", "compensate")],
                   list(pass = "log", compensate = "off"))
  seen <- character()
  read <- logical()
  spy <- function(real) {
    force(real)
    function(...) {
      ctl <- Filter(function(a) is.list(a) && !is.null(a$pass), list(...))[[1L]]
      seen <<- c(seen, paste(ctl$pass, ctl$compensate))
      r <- real(...)
      read <<- c(read, isTRUE(r$exp_kept) && r$exp_evals > 0)
      r
    }
  }
  X <- matrix(bm_repro_x)
  alt <- rep(1:4, 3L)
  M <- rep(4L, 3L)
  w <- rep(1, 3L)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, 1, 1, 1), X, alt, M, w))
  eta <- array(stats::rnorm(15L), c(1L, 5L, 3L))
  with_mocked_bindings({
    blp(bm_repro_fit("mnl"), target_shares = S)
    blp_contraction(rep(0, 4L), S, X, 1, alt, M, w)
  }, mnl_blp_solve = spy(choicer:::mnl_blp_solve))
  with_mocked_bindings({
    blp(bm_repro_fit("nl"), target_shares = S)
    nl_blp_contraction(rep(0, 4L), S, X, 1, c(1, 1), alt, c(1L, 1L, 2L, 2L), M, w)
  }, nl_blp_solve = spy(choicer:::nl_blp_solve))
  with_mocked_bindings({
    blp(bm_repro_fit("mxl"), target_shares = S)
    mxl_blp_contraction(rep(0, 4L), S, X, matrix(0, 12L, 1L), 1, 0, log(0.5), alt, M,
                        w, eta, 0L, FALSE, TRUE, FALSE)
  }, mxl_blp_solve = spy(choicer:::mxl_blp_solve))
  expect_identical(seen, rep("exp auto", 6L))
  expect_identical(read, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("blp() on a mixed logit: a free start, and what its budget keeps", {
  # (bit for bit at one thread: above it, sums over situations combine in a
  # varying order)
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dt <- create_small_mxl_data()
  fit <- suppressWarnings(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 10L,
    control = list(maxeval = 30L)))
  S <- as.numeric(predict(fit, type = "shares"))
  d0 <- blp(fit, target_shares = rev(S), keep_draws_bytes = 0)
  # the draws alone (8 K_w S N bytes) change nothing
  dG <- blp(fit, target_shares = rev(S),
            keep_draws_bytes = 8 * 2 * 10 * length(fit$data$M))
  expect_identical(dG, d0)
  # every situation's exponentials too: the shares agree to rounding, and
  # the constants to about tol
  dInf <- blp(fit, target_shares = rev(S), keep_draws_bytes = Inf)
  expect_equal(as.numeric(dInf), as.numeric(d0), tolerance = 1e-10)
  expect_true(attr(dInf, "convergence")$converged)
  expect_identical(rownames(d0), as.character(1:3))
  expect_true(attr(d0, "convergence")$converged)
  # without an outside option the start may omit the first constant's zero
  dfree <- blp(fit, target_shares = rev(S), delta_init = c(0.2, -0.1))
  dfull <- blp(fit, target_shares = rev(S), delta_init = c(0, 0.2, -0.1))
  expect_identical(dfree, dfull)
  expect_error(blp(fit, target_shares = S, delta_init = 0),
               "of 3 constants, one per inside alternative \\(or 2")
})
