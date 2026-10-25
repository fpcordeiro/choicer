# run_mnlogit(polish = ) and run_nestlogit(polish = ): the Newton phase on the
# multinomial and nested logits. A converged polish certifies the maximum the
# optimizer approached (a decrement below tol_decrement) and reports the
# optimizer's success code, keeping the optimizer's own status in
# fit$optimizer$polish; polish = "none" (the default) changes nothing.

# N situations over J alternatives in sets of lo to hi, x1 and x2, a weight
# and a cluster per situation; choices from the model at beta = (0.8, -0.5)
# and constants drawn once, an MNL, or an NL over `nests` with dissimilarity
# `lambda` in every nest of two or more; with `ioo`, the outside option can
# be chosen; with `alone`, each situation offers exactly one of those
# alternatives.
lp_data <- function(seed, N = 400L, J = 8L, lo = 3L, hi = 6L, ioo = FALSE,
                    nests = NULL, lambda = 0.6, alone = NULL) {
  set.seed(seed)
  M <- sample(lo:hi, N, replace = TRUE)
  pool <- setdiff(seq_len(J), alone)
  alt <- unlist(lapply(M, function(m) {
    if (length(alone)) sort(c(sample(pool, m - 1L), sample(alone, 1L)))
    else sort(sample.int(J, m))
  }))
  X <- matrix(stats::rnorm(sum(M) * 2L), ncol = 2L)
  delta <- c(if (!ioo) 0, stats::rnorm(J - (!ioo), sd = 0.5))
  start <- c(0L, cumsum(M))
  pr <- if (is.null(nests)) {
    V <- drop(X %*% c(0.8, -0.5)) + delta[alt]
    unlist(lapply(seq_len(N), function(i) {
      v <- c(if (ioo) 0, V[(start[i] + 1L):start[i + 1L]])
      p <- exp(v - max(v)); p <- p / sum(p)
      if (ioo) p[-1L] else p
    }))
  } else {
    K_l <- sum(tabulate(nests) > 1L)
    theta <- c(0.8, -0.5, rep_len(lambda, K_l), if (ioo) delta else delta[-1L])
    nl_predict(theta, X, as.integer(alt), as.integer(M), as.integer(nests),
               TRUE, ioo)$choice_prob
  }
  choice <- integer(length(alt))
  for (i in seq_len(N)) {
    p <- pr[(start[i] + 1L):start[i + 1L]]
    p <- c(if (ioo) max(0, 1 - sum(p)), p)
    k <- sample.int(length(p), 1L, prob = p) - ioo
    if (k >= 1L) choice[start[i] + k] <- 1L
  }
  dt <- data.table::data.table(id = rep(seq_len(N), M), alt = alt, choice = choice,
                               x1 = X[, 1L], x2 = X[, 2L])
  dt[, w := stats::runif(N, 0.5, 2)[id]]
  dt[, clus := (id - 1L) %/% 20L]
  if (!is.null(nests)) dt[, nest := nests[alt]]
  dt[]
}

LP_NESTS <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L)

# A fit, its messages muffled. Non-uniform weights make hessian/bhhh fits
# warn that they are not WESML corrections; that advice is beside the point
# here, and only it is muffled.
lp_run <- function(model, dt, ...) {
  if (model == "mnl") {
    run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"), ...)
  } else {
    run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest", ...)
  }
}
lp_fit <- function(model, dt, ...) {
  withCallingHandlers(
    suppressMessages(lp_run(model, dt, ...)),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "Non-uniform weights detected")) {
        invokeRestart("muffleWarning")
      }
    })
}

lp_dat <- function(model, seed, ...) {
  lp_data(seed, nests = if (model == "nl") LP_NESTS, ...)
}

# A fit's distance from a reference: the largest difference in a
# coefficient, over its absolute value or its standard error, whichever is
# larger
lp_dist <- function(fit, ref) {
  max(abs(coef(fit) - coef(ref)) / pmax(abs(coef(ref)), ref$se))
}

LP_TIGHT <- list(xtol_rel = 1e-14, maxeval = 10000L)
LP_LOOSE <- list(maxeval = 6L)
LP_STEPS <- list(max_iter = 100L)
LP_CERT <- list(max_iter = 100L, tol_decrement = 1e-16)

# The fit without its call and run time, for bitwise comparisons
lp_strip <- function(fit) {
  fit$call <- NULL
  fit$optimizer$elapsed_time <- NULL
  fit
}

# Half of g' C^-1 g at the fit's estimates, over the mean weight, with g the
# gradient and C the Hessian (or BHHH) from the fit's stored data
lp_decrement <- function(model, fit, curvature) {
  d <- fit$data
  g <- if (model == "mnl") {
    mnl_loglik_gradient_parallel(fit$coefficients, d$X, d$alt_idx, d$choice_idx,
                                 d$M, d$weights, fit$use_asc,
                                 fit$include_outside_option)$gradient
  } else {
    nl_loglik_gradient_parallel(fit$coefficients, d$X, d$alt_idx, d$choice_idx,
                                d$nest_idx, d$M, d$weights, fit$use_asc,
                                fit$include_outside_option)$gradient
  }
  fit$se_method <- if (curvature == "bhhh") "bhhh" else "hessian"
  0.5 * sum(g * solve(compute_hessian(fit), g)) / mean(d$weights)
}

test_that("polish = \"none\" ignores polish_control and records no polish", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 11L)
    f0 <- lp_fit(model, dt, control = LP_LOOSE)
    f1 <- lp_fit(model, dt, control = LP_LOOSE, polish = "none",
                 polish_control = list(max_iter = 3L))
    expect_true(identical(lp_strip(f1), lp_strip(f0), num.eq = FALSE,
                          single.NA = FALSE), label = model)
    expect_named(f0$optimizer, c("name", "control", "elapsed_time", "iterations"))
  }
})

test_that("a loose fit, polished, reaches the tightly converged maximum", {
  skip_on_cran()
  cases <- list(
    mnl = list("mnl", lp_dat("mnl", 21L)),
    mnl_outside = list("mnl", lp_data(22L, ioo = TRUE), include_outside_option = TRUE),
    mnl_weights = list("mnl", lp_dat("mnl", 23L), weights_col = "w"),
    nl = list("nl", lp_dat("nl", 24L)),
    nl_outside = list("nl", lp_data(25L, ioo = TRUE, nests = LP_NESTS),
                      include_outside_option = TRUE),
    nl_far = list("nl", lp_data(26L, N = 800L, nests = LP_NESTS, lambda = 0.3)),
    nl_weights = list("nl", lp_dat("nl", 27L), weights_col = "w"))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    model <- cs[[1L]]; dt <- cs[[2L]]; extra <- cs[-(1:2)]
    ref <- do.call(lp_fit, c(list(model, dt, control = LP_TIGHT), extra))
    for (pm in c("newton", "bhhh")) {
      fit <- do.call(lp_fit, c(list(model, dt, control = LP_LOOSE, polish = pm,
                                    polish_control = LP_CERT), extra))
      lab <- paste(nm, pm)
      p <- fit$optimizer$polish
      expect_true(p$converged, label = lab)
      expect_lt(lp_dist(fit, ref), 1e-6)
      expect_equal(fit$loglik, ref$loglik, tolerance = 1e-10, label = lab)
      expect_identical(fit$convergence, 1L, label = lab)
    }
  }
})

test_that("the record, print() and summary() report the polish", {
  # a curvature slowed down, so that the polish's run time shows in the sum
  slow <- function(f) {
    force(f)
    function(...) {
      Sys.sleep(0.02)
      f(...)
    }
  }
  local_mocked_bindings(
    mnl_loglik_hessian_parallel = slow(mnl_loglik_hessian_parallel),
    nl_loglik_hessian_parallel = slow(nl_loglik_hessian_parallel))
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 31L)
    fit <- lp_fit(model, dt, control = LP_LOOSE, polish = "newton",
                  polish_control = LP_STEPS)
    p <- fit$optimizer$polish
    expect_named(fit$optimizer, c("name", "control", "elapsed_time", "iterations",
                                  "polish"))
    expect_named(p, c("method", "converged", "iterations", "damped", "decrement",
                      "lambda", "message", "evaluations", "active", "held", "flat",
                      "improvement", "error", "elapsed_time", "optimizer"))
    expect_named(p$optimizer, c("convergence", "message", "iterations",
                                "elapsed_time"))
    expect_true(p$converged, label = model)
    expect_null(p$error)
    expect_identical(p$method, "newton")
    expect_named(p$evaluations, c("objective", "curvature"))
    expect_match(fit$message, "^Newton polish converged after [0-9]+ iterations? ")
    expect_match(fit$message, "the optimizer: ", fixed = TRUE)
    # the fit's run time is the optimizer's and the polish's
    expect_gte(p$elapsed_time, 0.02)
    expect_equal(fit$optimizer$elapsed_time,
                 p$optimizer$elapsed_time + p$elapsed_time)
    expect_output(print(fit), "Newton polish converged")
    expect_output(print(summary(fit)), "Newton polish converged")
    # one line on the console
    msgs <- character()
    withCallingHandlers(
      lp_run(model, dt, control = LP_LOOSE, polish = "bhhh"),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    expect_identical(sum(grepl("^BHHH polish: (not )?converged after", msgs)), 1L,
                     label = model)
  }
})

test_that("a polish that does not converge keeps the optimizer's code", {
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 41L)
    fit <- lp_fit(model, dt, control = list(maxeval = 2L), polish = "newton",
                  polish_control = list(max_iter = 1L, tol_decrement = 1e-30))
    p <- fit$optimizer$polish
    expect_false(p$converged, label = model)
    expect_identical(fit$convergence, p$optimizer$convergence)
    expect_match(fit$message, "Newton polish not converged after 1 iteration: ",
                 fixed = TRUE)
  }
})

test_that("the polish's curvature takes the fit's data and weights", {
  # With max_iter = 0 the polish only evaluates the decrement at the
  # optimizer's solution: a curvature on other data, or on squared weights,
  # would give another decrement there. (The stops are close enough to the
  # maximum for a positive definite Hessian, and far enough for a decrement
  # that is not negligible.)
  for (model in c("mnl", "nl")) for (oo in c(FALSE, TRUE)) {
    dt <- lp_data(51L + oo, ioo = oo, nests = if (model == "nl") LP_NESTS)
    for (pm in c("newton", "bhhh")) {
      lab <- paste(model, "oo =", oo, pm)
      fit <- lp_fit(model, dt, include_outside_option = oo, weights_col = "w",
                    control = list(maxeval = if (model == "mnl") 3L else 10L),
                    polish = pm, polish_control = list(max_iter = 0L))
      p <- fit$optimizer$polish
      expect_identical(p$iterations, 0L, label = lab)
      expect_identical(p$lambda, 0, label = lab)
      expect_gt(p$decrement, 1e-6)
      expect_equal(p$decrement, lp_decrement(model, fit, pm), tolerance = 1e-8,
                   label = lab)
    }
  }
  # and without constants
  for (model in c("mnl", "nl")) for (pm in c("newton", "bhhh")) {
    dt <- lp_dat(model, 53L)
    fit <- lp_fit(model, dt, use_asc = FALSE, control = list(maxeval = 5L),
                  polish = pm, polish_control = list(max_iter = 0L))
    expect_equal(fit$optimizer$polish$decrement, lp_decrement(model, fit, pm),
                 tolerance = 1e-8, label = paste(model, pm))
  }
})

test_that("the polish's decrement rule does not depend on the weights' scale", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 61L)
    dt16 <- data.table::copy(dt)[, w := 16 * w]
    # the polish from the same loose point at both scales, in the
    # parameters' own coordinates, which do not move with the weights
    start <- unname(coef(lp_fit(model, dt, weights_col = "w", scale_vars = "none",
                                control = LP_LOOSE)))
    hold <- function(theta_init, eval_f, ...) {
      list(par = start, value = eval_f(start)$objective, convergence = 0L)
    }
    run <- function(d, ctl) {
      lp_fit(model, d, weights_col = "w", scale_vars = "none", optimizer = hold,
             polish = "newton", polish_control = ctl)$optimizer$polish
    }
    # the decrement over the mean weight, bit for bit at a power of two
    nu <- run(dt, list(max_iter = 0L))$decrement
    expect_true(is.finite(nu) && nu > 0, label = model)
    expect_identical(run(dt16, list(max_iter = 0L))$decrement, nu, label = model)
    # a tolerance twice that certifies at once at both scales, where one in
    # raw log-likelihood units would not at 16 times the weights
    ctl <- list(tol_decrement = 2 * nu)
    a <- run(dt, ctl)
    b <- run(dt16, ctl)
    expect_true(a$converged, label = model)
    expect_identical(a$iterations, 0L)
    expect_identical(b[c("converged", "iterations", "decrement")],
                     a[c("converged", "iterations", "decrement")], label = model)
  }
})

test_that("the polish works in natural units, whatever the optimizer's coordinates", {
  skip_on_cran()
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 71L)
    fits <- lapply(c("none", "sd", "bhhh"), function(sv) {
      lp_fit(model, dt, scale_vars = sv, control = LP_LOOSE, polish = "newton",
             polish_control = LP_CERT)
    })
    for (f in fits) expect_true(f$optimizer$polish$converged, label = model)
    expect_lt(lp_dist(fits[[1L]], fits[[3L]]), 1e-7)
    expect_lt(lp_dist(fits[[2L]], fits[[3L]]), 1e-7)
  }
})

test_that("under scale_vars = \"none\" the steps' scales come from a BHHH pass", {
  skip_on_cran()
  # Columns in units 10^6 apart: in the parameters' own units (scales of 1)
  # the steps that the coefficient of the small-unit column needs, hundreds
  # of units, run past the stall rule's 100 and are held; at the BHHH
  # scales the polish certifies the maximum
  ones <- FALSE
  real <- .polish_scales
  n <- 0L
  local_mocked_bindings(.polish_scales = function(model_f, theta) {
    n <<- n + 1L
    if (ones) rep(1, length(theta)) else real(model_f, theta)
  })
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 81L, N = 1500L)
    dt[, `:=`(x1 = x1 * 1e3, x2 = x2 * 1e-3)]
    ref <- lp_fit(model, dt, scale_vars = "bhhh", control = LP_TIGHT,
                  polish = "newton", polish_control = LP_CERT)
    ones <- FALSE
    fit <- lp_fit(model, dt, scale_vars = "none", control = list(maxeval = 20L),
                  polish = "newton", polish_control = list(tol_decrement = 1e-16))
    expect_true(fit$optimizer$polish$converged, label = model)
    expect_lt(lp_dist(fit, ref), 1e-6)
    ones <- TRUE
    fit <- lp_fit(model, dt, scale_vars = "none", control = list(maxeval = 20L),
                  polish = "newton")
    p <- fit$optimizer$polish
    expect_false(p$converged, label = model)
    expect_identical(p$flat, "x2", label = model)
  }
  # the pass runs only where the fit's coordinates have no scales
  ones <- FALSE
  dt <- lp_dat("mnl", 82L)
  n <- 0L
  lp_fit("mnl", dt, scale_vars = "bhhh", polish = "newton")
  lp_fit("mnl", dt, scale_vars = "sd", polish = "newton")
  lp_fit("mnl", dt, scale_vars = "none")
  expect_identical(n, 0L)
  lp_fit("mnl", dt, scale_vars = "none", polish = "newton")
  expect_identical(n, 1L)
})

test_that("constants without a finite estimate are held at the optimizer's", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 91L)
    dt[alt == 8L, choice := 0L]
    dt <- dt[id %in% dt[choice == 1L, id]]
    fit_h <- function(...) {
      expect_warning(f <- lp_fit(model, dt, control = LP_LOOSE, ...),
                     class = "choicer_unbounded_constants")
      f
    }
    f0 <- fit_h()
    fit <- fit_h(polish = "newton", polish_control = LP_STEPS)
    p <- fit$optimizer$polish
    expect_identical(p$held, "ASC_8", label = model)
    expect_true(p$converged, label = model)
    # the optimizer's value, the others certified given it
    expect_identical(coef(fit)[["ASC_8"]], coef(f0)[["ASC_8"]], label = model)
    expect_gt(fit$loglik, f0$loglik)
  }
})

test_that("a never-chosen constant is held only while its nest's dissimilarity is at most 1", {
  # mode_choice with a bus never chosen in the public nest, whose
  # dissimilarity is about 1.95 at the maximum
  data(mode_choice, package = "choicer", envir = environment())
  mc <- data.table::as.data.table(mode_choice)
  b2 <- mc[mode == "bus"]
  b2[, `:=`(mode = "bus2", choice = 0L, gcost = gcost + 5)]
  mc2 <- rbind(mc, b2)[order(id, mode)]
  mc2[, nest := ifelse(mode == "car", "private", "public")]
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fit_mc <- function(...) {
    suppressMessages(suppressWarnings(run_nestlogit(
      mc2, "id", "mode", "choice", c("gcost", "wait"), nest_col = "nest", ...),
      classes = "choicer_unbounded_constants"))
  }
  # from a stop with lambda above 1: the constant is moved, not held (here
  # it has no finite estimate either, and falls until its derivatives no
  # longer move the decrement)
  f0 <- fit_mc(control = list(maxeval = 500L))
  expect_gt(coef(f0)[["Lambda_1"]], 1)
  fit <- fit_mc(control = list(maxeval = 500L), polish = "newton",
                polish_control = list(max_iter = 50L))
  p <- fit$optimizer$polish
  expect_length(p$held, 0L)
  expect_true(p$converged)
  expect_lt(coef(fit)[["ASC_bus2"]], coef(f0)[["ASC_bus2"]])
  # from a stop with lambda at most 1, held; the polish takes lambda above
  # 1, where the held constant can have a finite estimate: not certified
  f0 <- fit_mc(control = list(maxeval = 3L))
  expect_lte(coef(f0)[["Lambda_1"]], 1)
  fit <- fit_mc(control = list(maxeval = 3L), polish = "newton",
                polish_control = list(max_iter = 50L))
  p <- fit$optimizer$polish
  expect_identical(p$held, "ASC_bus2")
  expect_identical(coef(fit)[["ASC_bus2"]], coef(f0)[["ASC_bus2"]])
  expect_gt(coef(fit)[["Lambda_1"]], 1)
  expect_false(p$converged)
  expect_match(p$message, paste(
    "the Newton decrement is below tol_decrement, but the dissimilarity of",
    "the nest of the never-chosen alternative whose constant (ASC_bus2) was",
    "held ended above 1, where such a constant can have a finite estimate"),
    fixed = TRUE)
  # with that bus the base alternative, the reference never chosen: the
  # other constants move against it, and the same point is certified
  fit_b <- fit_mc(control = list(maxeval = 500L), polish = "newton",
                  polish_control = list(max_iter = 50L), outside_opt_label = "bus2")
  expect_true(fit_b$optimizer$polish$converged)
  expect_equal(fit_b$loglik, fit_mc(control = list(maxeval = 500L), polish = "newton",
                                    polish_control = list(max_iter = 50L))$loglik,
               tolerance = 1e-8)
})

test_that("a moved constant whose nest ends at or below 1 is not certified", {
  # alternative 8, in the third nest, never chosen; the optimizer stops with
  # that nest's dissimilarity at 1.5, so the polish moves the constant, and
  # ends with the dissimilarity near the data's 0.6
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- lp_dat("nl", 91L)
  dt[alt == 8L, choice := 0L]
  dt <- dt[id %in% dt[choice == 1L, id]]
  fit_n <- function(...) {
    suppressWarnings(lp_fit("nl", dt, ...), classes = "choicer_unbounded_constants")
  }
  f0 <- fit_n(control = LP_LOOSE)
  th <- unname(coef(f0))
  th[f0$param_map$lambda[3L]] <- 1.5
  stop_at <- function(theta_init, eval_f, lower, upper, control) {
    list(par = th, value = eval_f(th)$objective, convergence = 0L)
  }
  fit <- fit_n(optimizer = stop_at, polish = "newton", polish_control = LP_STEPS)
  p <- fit$optimizer$polish
  expect_length(p$held, 0L)
  expect_lte(coef(fit)[["Lambda_3"]], 1)
  expect_false(p$converged)
  expect_match(p$message, paste(
    "the dissimilarity of the nest of the never-chosen alternative whose",
    "constant (ASC_8) was moved ended at or below 1, where such a constant has",
    "no finite estimate"), fixed = TRUE)
  # a never-chosen alternative alone in its nest, and one chosen wherever
  # offered, are held whatever the dissimilarities
  nests <- c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L)
  dt <- lp_data(92L, nests = nests)
  dt[alt == 6L, choice := 0L]
  dt <- dt[id %in% dt[choice == 1L, id]]
  dt <- dt[!(id %in% dt[alt == 4L, id] & !id %in% dt[alt == 4L & choice == 1L, id])]
  f0 <- fit_n(control = LP_LOOSE)
  th <- unname(coef(f0))
  th[f0$param_map$lambda] <- 1.5
  fit <- fit_n(optimizer = stop_at, polish = "newton", polish_control = list(max_iter = 0L))
  expect_setequal(fit$optimizer$polish$held, c("ASC_4", "ASC_6"))
})

test_that("constants that diverge together leave the polish unconverged", {
  # the outside option never chosen: every inside constant rises together
  dt <- lp_data(101L, ioo = TRUE)
  dt <- dt[id %in% dt[choice == 1L, id]]
  fit_o <- function(d, ...) {
    expect_warning(f <- lp_fit("mnl", d, include_outside_option = TRUE,
                               control = list(maxeval = 200L), polish = "newton",
                               ...),
                   class = "choicer_unbounded_constants")
    f
  }
  fit <- fit_o(dt, polish_control = list(max_iter = 0L, tol_decrement = 1e3))
  p <- fit$optimizer$polish
  expect_false(p$converged)
  expect_length(p$held, 0L)
  expect_match(p$message, paste0(
    "the Newton decrement is below tol_decrement, but the constants diverge ",
    "together (the outside option is never chosen), so there is no maximum ",
    "to certify"), fixed = TRUE)
  # holding the constant of an alternative chosen wherever offered, which
  # the situations that choose it offer with the others, stops their rise:
  # certified
  dt_held <- dt[!(id %in% dt[alt == 8L, id] & !id %in% dt[alt == 8L & choice == 1L, id])]
  fit <- fit_o(dt_held, polish_control = list(max_iter = 0L, tol_decrement = 1e3))
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$optimizer$polish$held, "ASC_8")
  # held where it is offered alone, it stops nothing
  alone <- data.table::data.table(id = max(dt$id) + seq_len(40L), alt = 8L,
                                  choice = 1L, x1 = stats::rnorm(40L),
                                  x2 = stats::rnorm(40L), w = 1, clus = 0L)
  dt_alone <- rbind(dt[alt != 8L], alone)
  dt_alone <- dt_alone[id %in% dt_alone[choice == 1L, id]]
  fit <- fit_o(dt_alone, polish_control = list(max_iter = 0L, tol_decrement = 1e3))
  expect_identical(fit$optimizer$polish$held, "ASC_8")
  expect_false(fit$optimizer$polish$converged)
  expect_match(fit$optimizer$polish$message, "the constants diverge together", fixed = TRUE)
})

test_that("\"bhhh\" with fewer choice situations than parameters stops", {
  # J alternatives offered in each of N situations, situation i choosing
  # alternative chosen[i]
  lp_small <- function(J, chosen) {
    N <- length(chosen)
    set.seed(J + N)
    dt <- data.table::data.table(id = rep(seq_len(N), each = J),
                                 alt = rep(seq_len(J), N),
                                 x1 = stats::rnorm(N * J), x2 = stats::rnorm(N * J))
    dt[, choice := as.integer(alt == chosen[id])]
    dt[, nest := c(1L, 1L, 2L, 2L, 3L, 3L)[alt]]
    dt[]
  }
  # 4 situations, 2 coefficients and 3 constants, every one with an estimate
  dt <- lp_small(4L, 1:4)
  for (model in c("mnl", "nl")) {
    expect_error(lp_fit(model, dt, polish = "bhhh"), paste0(
      "polish = \"bhhh\" needs at least as many likelihood units as the ",
      "parameters it moves: the BHHH matrix of 4 choice situations has rank at ",
      "most 4, below the ", if (model == "mnl") 5 else 7, " parameters"),
      fixed = TRUE)
  }
  # held constants are not moved, so they do not count: 5 situations, 2
  # coefficients and 5 constants, two of them never chosen
  dt <- lp_small(6L, c(1:4, 1L))
  expect_warning(f <- lp_fit("mnl", dt, polish = "bhhh",
                             polish_control = list(max_iter = 0L)),
                 class = "choicer_unbounded_constants")
  expect_identical(f$optimizer$polish$held, c("ASC_5", "ASC_6"))
  expect_null(f$optimizer$polish$error)
  # in the nested logit a never-chosen constant in a nest of two or more
  # counts as moving, since the polish moves it where its nest's
  # dissimilarity exceeds 1: 9 situations, 2 coefficients, 3 dissimilarities
  # and 5 constants
  dt <- lp_small(6L, c(1:5, 1:3, 5L))
  expect_error(suppressWarnings(lp_fit("nl", dt, polish = "bhhh"),
                                classes = "choicer_unbounded_constants"),
               "the BHHH matrix of 9 choice situations has rank at most 9, below the 10 parameters",
               fixed = TRUE)
})

test_that("keep_data = FALSE and the advanced pathway polish the same fit", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- lp_dat("nl", 121L)
  f1 <- lp_fit("nl", dt, control = LP_LOOSE, polish = "newton")
  f2 <- lp_fit("nl", dt, control = LP_LOOSE, polish = "newton", keep_data = FALSE)
  inp <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
  f3 <- suppressMessages(run_nestlogit(input_data = inp, control = LP_LOOSE,
                                       polish = "newton"))
  expect_identical(coef(f2), coef(f1))
  expect_identical(coef(f3), coef(f1))
  expect_identical(f2$vcov, f1$vcov)
  expect_null(f2$data)
})

test_that("optim and a custom optimizer hand their solution to the polish", {
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 131L)
    f_o <- lp_fit(model, dt, optimizer = "optim", control = list(maxit = 3L),
                  polish = "newton", polish_control = LP_STEPS)
    expect_true(f_o$optimizer$polish$converged, label = model)
    expect_identical(f_o$convergence, 0L)
    custom <- function(theta_init, eval_f, lower, upper, control) {
      list(par = theta_init, value = eval_f(theta_init)$objective, convergence = 7L)
    }
    f_c <- lp_fit(model, dt, optimizer = custom, polish = "newton",
                  polish_control = LP_STEPS)
    expect_true(f_c$optimizer$polish$converged, label = model)
    expect_identical(f_c$convergence, 0L)
    expect_identical(f_c$optimizer$polish$optimizer$convergence, 7L)
  }
  # an optimizer that ignores the dissimilarities' bound (optim's BFGS drops
  # it): the polish starts from the solution projected onto the bounds, and
  # its steps stay there
  dt <- lp_dat("nl", 132L)
  ref <- lp_fit("nl", dt, control = LP_TIGHT)
  # (the steps' scales, under scale_vars = "none", from the projected point)
  seen <- NULL
  real_ps <- .polish_scales
  local_mocked_bindings(.polish_scales = function(model_f, theta) {
    seen <<- theta
    real_ps(model_f, theta)
  })
  k <- ref$param_map$lambda[2L]
  for (bad in c(-0.5, 0)) {
    th <- unname(coef(ref))
    th[k] <- bad
    outside <- function(theta_init, eval_f, lower, upper, control) {
      list(par = th, value = NA_real_, convergence = 0L)
    }
    seen <- NULL
    f <- lp_fit("nl", dt, optimizer = outside, polish = "newton")
    expect_identical(seen[k], 1e-16)
    expect_true(all(coef(f)[f$param_map$lambda] >= 1e-16))
    expect_null(f$optimizer$polish$error)
  }
})

test_that("a BHHH pass for the steps' scales that fails leaves them at 1", {
  polish_fn <- run_polish
  seen <- NULL
  called <- FALSE
  local_mocked_bindings(
    .polish_scales = function(...) stop("the pass failed"),
    run_polish = function(theta, eval_f, curv_f, ...) {
      called <<- TRUE
      seen <<- list(...)$scale
      polish_fn(theta, eval_f, curv_f, ...)
    })
  dt <- lp_dat("nl", 133L)
  fit <- lp_fit("nl", dt, control = LP_LOOSE, polish = "newton",
                polish_control = LP_STEPS)
  p <- fit$optimizer$polish
  # (no scales: run_polish() takes 1 for each)
  expect_true(called)
  expect_null(seen)
  expect_null(p$error)
  expect_true(p$converged)
})

test_that("the variance reuses the polish's last curvature where se_method needs it", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  n_h <- n_b <- n_num <- 0L
  count <- function(f, counter) {
    force(f)
    function(...) {
      assign(counter, get(counter, envir = cnt) + 1L, envir = cnt)
      f(...)
    }
  }
  cnt <- environment()
  local_mocked_bindings(
    mnl_loglik_hessian_parallel = count(mnl_loglik_hessian_parallel, "n_h"),
    mnl_bhhh_parallel = count(mnl_bhhh_parallel, "n_b"),
    nl_loglik_hessian_parallel = count(nl_loglik_hessian_parallel, "n_h"),
    nl_bhhh_parallel = count(nl_bhhh_parallel, "n_b"),
    nl_loglik_numeric_hessian = count(nl_loglik_numeric_hessian, "n_num"))
  # the Hessian and BHHH calls beyond the polish's: the variance's own, none
  # where it takes the polish's last curvature (the robust meat squares the
  # weights, so "sandwich" evaluates the BHHH kernel; "cluster" has its own)
  cases <- list(
    list("newton", "hessian", h = 0L, b = 0L),
    list("newton", "sandwich", h = 0L, b = 1L),
    list("newton", "cluster", h = 0L, b = 0L),
    list("newton", "bhhh", h = 0L, b = 1L),
    list("bhhh", "bhhh", h = 0L, b = 0L),
    list("bhhh", "hessian", h = 1L, b = 0L),
    list("bhhh", "sandwich", h = 1L, b = 1L))
  for (model in c("mnl", "nl")) {
    dt <- lp_dat(model, 141L)
    for (cs in cases) {
      pm <- cs[[1L]]; se <- cs[[2L]]
      lab <- paste(model, pm, "/", se)
      n_h <- 0L; n_b <- 0L
      fit <- lp_fit(model, dt, control = LP_LOOSE, polish = pm, se_method = se,
                    cluster_col = if (se == "cluster") "clus",
                    polish_control = LP_STEPS)
      p <- fit$optimizer$polish
      expect_true(p$converged, label = lab)
      cc <- p$evaluations[["curvature"]]
      expect_identical(n_h, cs$h + if (pm == "newton") cc else 0L, label = lab)
      expect_identical(n_b, cs$b + if (pm == "bhhh") cc else 0L, label = lab)
      # the variance is the one computed afresh at the estimates
      v <- suppressMessages(vcov(fit, type = switch(se, sandwich = "robust", se)))
      expect_identical(fit$vcov, v, label = lab)
    }
  }
  # "numeric" takes its own numeric Hessian, once, never the polish's
  # analytic one
  n_h <- 0L; n_num <- 0L
  fit <- lp_fit("nl", dt, control = LP_LOOSE, polish = "newton", se_method = "numeric",
                polish_control = LP_STEPS)
  expect_identical(n_h, fit$optimizer$polish$evaluations[["curvature"]])
  expect_identical(n_num, 1L)
  expect_identical(fit$vcov, suppressMessages(vcov(fit, type = "hessian")))
})

test_that("an error in the curvature keeps the optimizer's estimates", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- lp_dat("mnl", 151L)
  ref <- lp_fit("mnl", dt, control = LP_LOOSE)
  # the polish's first curvature fails, the variance's does not
  n <- 0L
  real <- mnl_loglik_hessian_parallel
  local_mocked_bindings(mnl_loglik_hessian_parallel = function(...) {
    n <<- n + 1L
    if (n == 1L) stop("the curvature failed")
    real(...)
  })
  expect_warning(fit <- lp_fit("mnl", dt, control = LP_LOOSE, polish = "newton"),
                 class = "choicer_polish_failed")
  p <- fit$optimizer$polish
  expect_identical(p$error, "the curvature failed")
  expect_false(p$converged)
  expect_identical(coef(fit), coef(ref))
  # the variance is the fit's own, computed afresh
  expect_identical(fit$vcov, ref$vcov)
})

test_that("a trial point whose objective is not the likelihood fails", {
  # the kernel's overflow flag set at every point but the optimizer's
  # solution: no step is accepted
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- lp_dat("mnl", 161L)
  ref <- lp_fit("mnl", dt, control = LP_LOOSE)
  real <- mnl_loglik_gradient_parallel
  local_mocked_bindings(mnl_loglik_gradient_parallel = function(theta, ...) {
    r <- real(theta, ...)
    if (!isTRUE(all.equal(unname(theta), unname(coef(ref)), tolerance = 0))) {
      r$overflow <- TRUE
    }
    r
  })
  fit <- lp_fit("mnl", dt, control = LP_LOOSE, polish = "newton")
  p <- fit$optimizer$polish
  expect_false(p$converged)
  expect_identical(p$iterations, 0L)
  expect_identical(coef(fit), coef(ref))
})

test_that("a nest that no situation offers twice leaves its parameter unidentified", {
  # alternatives 7 and 8, the third nest of two, are never offered together:
  # the likelihood does not depend on that nest's dissimilarity. At the
  # default start of 0.5 its terms vanish exactly (1/lambda is exact there),
  # so the polish leaves it where the optimizer did, and BHHH reports it
  # without information. (The Hessian's diagonal entry, zero in exact
  # arithmetic, carries a rounding residue, so the Newton polish may certify
  # the rest.)
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  nests <- c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L)
  dt <- lp_data(171L, N = 300L, nests = nests, alone = 7:8)
  f0 <- lp_fit("nl", dt, control = LP_LOOSE)
  for (pm in c("newton", "bhhh")) {
    fit <- lp_fit("nl", dt, control = LP_LOOSE, polish = pm,
                  polish_control = LP_STEPS)
    p <- fit$optimizer$polish
    expect_identical(coef(fit)[["Lambda_3"]], coef(f0)[["Lambda_3"]], label = pm)
    expect_null(p$error)
  }
  expect_false(p$converged)
  expect_identical(p$flat, "Lambda_3")
  expect_match(p$message, "Lambda_3 has no usable curvature at this point (no information)",
               fixed = TRUE)
  # from 0.6 its gradient and curvature are rounding noise: BHHH still
  # leaves it out of its steps and does not certify
  start <- c(0, 0, 0.5, 0.5, 0.6, rep(0, 7))
  fit <- lp_fit("nl", dt, theta_init = start, control = LP_LOOSE, polish = "bhhh",
                polish_control = LP_STEPS)
  p <- fit$optimizer$polish
  expect_null(p$error)
  expect_false(p$converged)
  expect_true("Lambda_3" %in% p$flat)
})

test_that("the anchors against a common rise are formed before the optimizer", {
  local_mocked_bindings(.polish_anchors = function(...) stop("the anchors failed"))
  ran <- FALSE
  watch <- function(theta_init, eval_f, lower, upper, control) {
    ran <<- TRUE
    list(par = theta_init, value = eval_f(theta_init)$objective, convergence = 0L)
  }
  for (model in c("mnl", "nl")) {
    expect_error(lp_fit(model, lp_dat(model, 181L), optimizer = watch, polish = "newton"),
                 "the anchors failed", fixed = TRUE)
    expect_false(ran, label = model)
  }
})

test_that("polish and polish_control are checked before the data", {
  expect_error(run_mnlogit(polish = "lbfgs"), "should be one of")
  expect_error(run_nestlogit(polish = "lbfgs"), "should be one of")
  expect_error(run_mnlogit(polish_control = list(steps = 3L)),
               "Unknown `polish_control` setting(s): steps", fixed = TRUE)
  expect_error(run_nestlogit(polish_control = list(max_iter = -1L)),
               "`polish_control$max_iter` must be a whole number >= 0.", fixed = TRUE)
})
