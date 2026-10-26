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
    expect_identical(cv$method, "simple")
    expect_identical(cv$tol, 1e-8)
    expect_gte(cv$evaluations, 1L)
    expect_lte(cv$residual, 1e-8)
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
  expect_error(blp(fit, target_shares = S, method = "anderson"), "must be \"simple\"")
  expect_error(blp(fit, target_shares = S, control = list(m = 5)),
               "`control` takes no settings with method \"simple\" \\(it has `m`\\)")
  expect_error(blp(fit, target_shares = S, control = "x"), "must be a list")
  expect_error(blp(fit, target_shares = S, tol = -1), "`tol` must be a single")
  expect_error(blp(fit, target_shares = S, max_iter = 2.5), "`max_iter` must be a whole")
})

test_that("blp() on a nested logit labels, damps and checks damping", {
  dt <- create_small_nl_data()
  fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                       control = list(maxeval = 50L))
  S <- as.numeric(predict(fit, type = "shares"))
  d1 <- blp(fit, target_shares = rev(S))
  d2 <- blp(fit, target_shares = rev(S), damping = 0.5)
  expect_identical(rownames(d1), as.character(1:6))
  expect_equal(as.numeric(d2), as.numeric(d1), tolerance = 1e-6)
  expect_gt(attr(d2, "convergence")$evaluations, attr(d1, "convergence")$evaluations)
  expect_error(blp(fit, target_shares = S, damping = 0), "`damping` must be a single")
})

test_that("blp() on a mixed logit: a free start, and kept draws change nothing", {
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
  dInf <- blp(fit, target_shares = rev(S), keep_draws_bytes = Inf)
  expect_identical(d0, dInf)
  expect_identical(rownames(d0), as.character(1:3))
  expect_true(attr(d0, "convergence")$converged)
  # without an outside option the start may omit the first constant's zero
  dfree <- blp(fit, target_shares = rev(S), delta_init = c(0.2, -0.1))
  dfull <- blp(fit, target_shares = rev(S), delta_init = c(0, 0.2, -0.1))
  expect_identical(dfree, dfull)
  expect_error(blp(fit, target_shares = S, delta_init = 0),
               "of 3 constants, one per inside alternative \\(or 2")
})
