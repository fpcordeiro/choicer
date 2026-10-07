# scale_vars as the optimizer's coordinates: theta = m * theta_t + c. The
# kernels evaluate the model at theta; the optimizer works on theta_t, with
# the gradient times m. "bhhh" takes m from the BHHH diagonal at the start
# values, in powers of two.

test_that("column scales give the map main's fits built", {
  pm <- list(beta = 1:2, mu = 3:4, sigma = 5:7, asc = 8:9)
  sX <- c(3, 0.25); sW <- c(7, 0.5)
  m <- choicer:::.coordinate_map("sd", pm, 9L, sX = sX, sW = sW,
                                 rc_correlation = TRUE)
  # Cholesky rows (row-major lower: L11, L21, L22): row i divided by sW[i],
  # the log diagonal shifted by -log(sW[i])
  expect_identical(m$scale, c(1 / 3, 4, 1 / 7, 2, 1, 2, 1, 1, 1))
  expect_identical(m$shift, c(0, 0, 0, 0, -log(7), 0, -log(0.5), 0, 0))
  pu <- list(beta = 1:2, sigma = 3:4, asc = 5L)
  mu <- choicer:::.coordinate_map("mad", pu, 5L, sX = sX, sW = sW)
  expect_identical(mu$scale, c(1 / 3, 4, 1, 1, 1))
  expect_identical(mu$shift, c(0, 0, -log(7), -log(0.5), 0))
  # K_w = 3: L11, L21, L22, L31, L32, L33
  p3 <- list(beta = 1L, sigma = 2:7)
  m3 <- choicer:::.coordinate_map("iqr", p3, 7L, sX = 3, sW = c(7, 0.5, 4),
                                  rc_correlation = TRUE)
  expect_identical(m3$scale, c(1 / 3, 1, 2, 1, 0.25, 0.25, 1))
  expect_identical(m3$shift, c(0, -log(7), 0, -log(0.5), 0, 0, -log(4)))
  # MNL and NL: beta only.
  pn <- list(beta = 1:2, lambda = 3L, asc = 4:5)
  mn <- choicer:::.coordinate_map("iqr", pn, 5L, sX = sX)
  expect_identical(mn$scale, c(1 / 3, 4, 1, 1, 1))
  expect_identical(mn$shift, rep(0, 5))
  none <- choicer:::.coordinate_map("none", pm, 9L)
  expect_identical(none, list(scale = rep(1, 9), shift = rep(0, 9)))
})

test_that("BHHH scales are floored, clamped powers of two", {
  B <- c(4, 1e-6, 2.5e7, 0, NaN, Inf, 1e-20, 16)
  m <- choicer:::.bhhh_scales(B)
  expect_identical(log2(m), round(log2(m)))
  med <- stats::median(B[is.finite(B) & B > 0])
  expect_identical(m[c(1L, 8L)], c(0.5, 0.25))
  # zero, NaN, Inf and an entry below 1e-12 times the median take the median
  expect_identical(m[4:7], rep(2^round(log2(1 / sqrt(med))), 4))
  expect_identical(m[2:3], c(2^10, 2^-12))
  # within the floor, but beyond the clamp
  expect_identical(choicer:::.bhhh_scales(c(1e-300, 1e-295, 1e-290)),
                   rep(2^400, 3))
  expect_identical(choicer:::.bhhh_scales(c(1e300, 1e295, 1e290)),
                   rep(2^-400, 3))
  expect_message(none <- choicer:::.bhhh_scales(c(NA, 0, -1)),
                 "no finite, positive entry")
  expect_identical(none, rep(1, 3))
})

test_that("the map round-trips, exactly for powers of two", {
  set.seed(1)
  th <- c(stats::rnorm(6) * 10^(-3:2), -0, 0)
  mb <- list(scale = 2^c(-30, -3, 0, 4, 9, 40, 2, -2), shift = rep(0, 8))
  expect_identical(choicer:::.from_coordinates(choicer:::.to_coordinates(th, mb), mb), th)
  # -0 stays -0 without a shift
  expect_identical(1 / choicer:::.from_coordinates(-0, list(scale = 2, shift = 0)), -Inf)
  ms <- list(scale = c(1 / 3, 1 / 7, 1, 1, 0.1, 1, 1, 1),
             shift = c(0, 0, -log(7), 0, 0, 0, 0, 0))
  expect_equal(choicer:::.from_coordinates(choicer:::.to_coordinates(th, ms), ms),
               th, tolerance = 1e-15)
  # infinite bounds stay infinite
  expect_identical(choicer:::.to_coordinates(c(-Inf, Inf), list(scale = c(4, 0.5),
                                                                 shift = c(1, 0))),
                   c(-Inf, Inf))
})

test_that("the coordinate objective's gradient is the chain rule's", {
  skip_if_not_installed("numDeriv")
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  fx <- mxlp_fixture("coord", 41L, rc_dist = c(0L, 1L), rc_mean = TRUE,
                     rc_correlation = TRUE, include_outside_option = TRUE,
                     weight_type = "person")
  f <- function(theta) mxlp_call("gradient", fx, theta = theta)
  B <- diag(mxlp_call("bhhh", fx))
  map <- list(scale = choicer:::.bhhh_scales(B), shift = rep(0, length(B)))
  ft <- choicer:::.coordinate_objective(f, map)
  th_t <- choicer:::.to_coordinates(fx$theta, map)
  r <- ft(th_t)
  expect_identical(r$objective, f(fx$theta)$objective)
  num <- numDeriv::grad(function(x) ft(x)$objective, th_t)
  expect_equal(as.numeric(r$gradient), num, tolerance = 1e-6)
})

test_that("the coordinate objective is the scaled design's objective", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  for (cfg in list(list(rc_dist = c(0L, 0L), corr = TRUE, mean = FALSE, panel = TRUE),
                   list(rc_dist = c(0L, 1L), corr = FALSE, mean = TRUE, panel = FALSE))) {
    fx <- mxlp_fixture("eqv", 45L, rc_dist = cfg$rc_dist,
                       rc_correlation = cfg$corr, rc_mean = cfg$mean,
                       weight_type = "person")
    if (!cfg$panel) fx <- mxlp_cross_section(fx, 46L)
    fx$X <- fx$X %*% diag(c(100, 0.025))
    fx$W <- fx$W %*% diag(c(25, 1))
    pm <- list(beta = 1:2)
    pos <- 2L
    if (cfg$mean) { pm$mu <- pos + 1:2; pos <- pos + 2L }
    L_size <- if (cfg$corr) 3L else 2L
    pm$sigma <- pos + seq_len(L_size)
    n <- length(fx$theta)
    for (method in c("sd", "mad", "iqr")) {
      sX <- choicer:::.column_scales(fx$X, method)
      sW <- choicer:::.column_scales(fx$W, method)
      sW[cfg$rc_dist == 1L] <- 1
      map <- choicer:::.coordinate_map(method, pm, n, sX = sX, sW = sW,
                                       rc_correlation = cfg$corr)
      th_t <- choicer:::.to_coordinates(fx$theta, map)
      for (generate in c(FALSE, TRUE)) {
        coord <- choicer:::.coordinate_objective(
          function(theta) mxlp_call("gradient", fx, theta = theta,
                                    generate = generate), map)(th_t)
        fs <- fx
        fs$X <- sweep(fx$X, 2, sX, "/")
        fs$W <- sweep(fx$W, 2, sW, "/")
        scaled <- mxlp_call("gradient", fs, theta = th_t, generate = generate)
        # Two roundings of one function: far below a wrong map's 1e-2..1
        what <- sprintf("[%s, %s, %s]", method,
                        if (cfg$panel) "panel" else "cross-section",
                        if (generate) "generate" else "store")
        mxlp_expect_close(coord$objective, scaled$objective, 1e-11,
                          paste("objective", what))
        mxlp_expect_close(coord$gradient, scaled$gradient, 1e-11,
                          paste("gradient", what))
      }
    }
  }
})

# Simulated mixed logit data with badly scaled columns and decision-maker
# weights.
sc_data <- function(N, T, seed, ioo = FALSE) {
  sim <- simulate_mxl_data(N = N, J = 4L, T = T, beta = c(0.8, -0.5),
                           Sigma = diag(c(1, 0.6)), outside_option = ioo,
                           vary_choice_set = TRUE, seed = seed)
  dt <- data.table::as.data.table(sim$data)
  if (ioo) dt <- dt[alt != 0]
  dt[, `:=`(x1 = x1 * 10, x2 = x2 / 40, w1 = w1 * 5, w2 = abs(w2) + 0.1)]
  unit <- if (T > 1L) dt$pid else dt$id
  u <- match(unit, unique(unit))
  dt[, w := 0.5 + (u %% 5L) / 4]
  dt
}

sc_fit <- function(dt, panel = FALSE, ...) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    person_col = if (panel) "pid", S = 20L, seed = 3L,
    control = list(maxeval = 500L, xtol_rel = 1e-10), ...))
}

test_that("\"bhhh\" fits reach the unscaled optimum", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  cases <- list(
    list(panel = FALSE, ioo = FALSE, args = list(rc_correlation = TRUE)),
    list(panel = TRUE, ioo = TRUE,
         args = list(draws = "generate", rc_dist = c(0L, 1L), rc_mean = TRUE,
                     include_outside_option = TRUE, weights_col = "w",
                     se_method = "sandwich")))
  first <- NULL
  for (cs in cases) {
    dt <- sc_data(if (cs$panel) 80L else 200L, if (cs$panel) 3L else 1L,
                  51L + cs$panel, cs$ioo)
    none <- do.call(sc_fit, c(list(dt, cs$panel, scale_vars = "none"), cs$args))
    bh <- do.call(sc_fit, c(list(dt, cs$panel, scale_vars = "bhhh"), cs$args))
    # The unscaled optimizer polishing from the "bhhh" estimates: the optimum
    # both reach (unscaled from its cold start may stop short of it).
    ref <- do.call(sc_fit, c(list(dt, cs$panel, scale_vars = "none",
                                  theta_init = coef(bh)), cs$args))
    what <- if (cs$panel) "panel" else "cross-section"
    expect_gte(bh$loglik, none$loglik - 1e-6, label = what)
    expect_equal(bh$loglik, ref$loglik, tolerance = 1e-8, label = what)
    expect_true(all(is.finite(bh$se)))
    expect_lt(max(abs(coef(bh) - coef(ref)) / bh$se), 1e-3)
    expect_identical(log2(bh$param_scale), round(log2(bh$param_scale)))
    expect_identical(unname(bh$param_shift), rep(0, length(coef(bh))))
    expect_identical(unname(bh$sX), c(1, 1))
    expect_identical(names(bh$sW), c("w1", "w2"))
    type <- if (identical(bh$se_method, "sandwich")) "robust" else bh$se_method
    expect_identical(suppressMessages(vcov(bh, type = type)), bh$vcov)
    if (is.null(first)) first <- list(dt = dt, none = none, bh = bh, ref = ref)
  }
  # A warm start (the unscaled estimates) reaches the same optimum in fewer
  # iterations than the cold start.
  warm <- sc_fit(first$dt, rc_correlation = TRUE, scale_vars = "bhhh",
                 theta_init = coef(first$none))
  expect_lt(warm$optimizer$iterations, first$bh$optimizer$iterations)
  expect_equal(warm$loglik, first$ref$loglik, tolerance = 1e-8)
})

test_that("\"bhhh\" leaves the parameters unscaled where the likelihood overflows", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sc_data(60L, 1L, 61L)
  th <- c(4e306, 0, log(0.5), log(0.5), 0, 0, 0)  # x1's coefficient overflows
  msgs <- capture_messages(
    fit <- run_mxlogit(dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"),
                       S = 10L, scale_vars = "bhhh", theta_init = th,
                       control = list(maxeval = 20L)))
  expect_true(any(grepl("the likelihood overflows at the start values", msgs)))
  expect_identical(unname(fit$param_scale), rep(1, 7))
})

test_that("\"none\" hands the optimizer the natural problem", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sc_data(60L, 1L, 71L)
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- list(theta_init = theta_init, lower = lower, upper = upper,
                  r = eval_f(theta_init))
    stop("captured")
  })
  expect_error(run_mxlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                           c("w1", "w2"), S = 10L, scale_vars = "none"),
               "captured")
  expect_identical(seen$theta_init, c(0, 0, log(0.5), log(0.5), 0, 0, 0))
  expect_identical(seen$lower, rep(-Inf, 7))
  expect_identical(seen$upper, rep(Inf, 7))
  d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"))
  k <- mxl_loglik_gradient_parallel(
    theta = seen$theta_init, X = d$X, W = d$W, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
    eta_draws = get_halton_normals(10L, d$N, 2L), rc_dist = c(0L, 0L),
    rc_correlation = FALSE, rc_mean = FALSE)
  expect_identical(seen$r$objective, k$objective)
  expect_identical(seen$r$gradient, k$gradient)
  expect_null(seen$r$overflow)  # .lift_sentinel() strips it
})

test_that("\"bhhh\" hands the optimizer mapped start values and bounds", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sc_data(60L, 1L, 71L)
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- list(theta_init = theta_init, lower = lower, upper = upper,
                  r = eval_f(theta_init))
    stop("captured")
  })
  # No scale_vars: the default is "bhhh"
  expect_error(run_mxlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                           c("w1", "w2"), S = 10L,
                           lower = c(x2 = -30, L_11 = -3),
                           upper = c(x1 = 0.5)), "captured")
  th0 <- c(0, 0, log(0.5), log(0.5), 0, 0, 0)
  d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"))
  k <- mxl_loglik_gradient_parallel(
    theta = th0, X = d$X, W = d$W, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
    eta_draws = get_halton_normals(10L, d$N, 2L), rc_dist = c(0L, 0L),
    rc_correlation = FALSE, rc_mean = FALSE, opg_diag = TRUE)
  m <- choicer:::.bhhh_scales(k$opg_diag)
  expect_identical(log2(m), round(log2(m)))
  expect_identical(seen$theta_init, th0 / m)
  expect_identical(unname(seen$lower), c(-Inf, -30, -3, -Inf, -Inf, -Inf, -Inf) / m)
  expect_identical(unname(seen$upper), c(0.5, Inf, Inf, Inf, Inf, Inf, Inf) / m)
  expect_identical(seen$r$objective, k$objective)
  expect_identical(seen$r$gradient, k$gradient * m)
  expect_null(seen$r$overflow)
})

test_that("\"sd\" hands the optimizer the column map's start values and bounds", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sc_data(60L, 1L, 71L)
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- list(theta_init = theta_init, lower = lower,
                  r = eval_f(theta_init))
    stop("captured")
  })
  expect_error(run_mxlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                           c("w1", "w2"), S = 10L, scale_vars = "sd",
                           lower = c(L_11 = -3)), "captured")
  d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"))
  map <- choicer:::.coordinate_map(
    "sd", list(beta = 1:2, sigma = 3:4, asc = 5:7), 7L,
    sX = choicer:::.column_scales(d$X, "sd"),
    sW = choicer:::.column_scales(d$W, "sd"))
  th0 <- c(0, 0, log(0.5), log(0.5), 0, 0, 0)
  expect_identical(seen$theta_init, (th0 - map$shift) / map$scale)
  # The log diagonal's shift moves the bound: L_11 >= -3 natural
  expect_identical(unname(seen$lower), c(-Inf, -Inf, -3 - map$shift[3],
                                         rep(-Inf, 4)))
  k <- mxl_loglik_gradient_parallel(
    theta = choicer:::.from_coordinates(seen$theta_init, map), X = d$X,
    W = d$W, alt_idx = d$alt_idx, choice_idx = d$choice_idx, M = d$M,
    weights = d$weights, eta_draws = get_halton_normals(10L, d$N, 2L),
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE)
  expect_identical(seen$r$objective, k$objective)
  expect_identical(seen$r$gradient, k$gradient * map$scale)
})

test_that("a \"bhhh\" fit keeps its natural bounds", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sc_data(200L, 1L, 51L)
  free <- sc_fit(dt, scale_vars = "bhhh")
  ub <- coef(free)[["x2"]] - 1  # binds: one unit below the free estimate
  th <- coef(free)
  th[["x2"]] <- ub - 1           # a feasible start
  bnd <- sc_fit(dt, scale_vars = "bhhh", upper = c(x2 = ub),
                theta_init = unname(th))
  expect_lte(coef(bnd)[["x2"]], ub)
  expect_equal(coef(bnd)[["x2"]], ub, tolerance = 1e-8)
  expect_lt(bnd$loglik, free$loglik)
})

test_that("multinomial scaled fits reach the unscaled optimum", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 81, outside_option = TRUE)
  dm <- data.table::as.data.table(sim$data)
  dm[, x1 := x1 * 100]
  fit <- function(sv) suppressMessages(run_mnlogit(
    dm, "id", "alt", "choice", c("x1", "x2"), outside_opt_label = 0L,
    include_outside_option = TRUE, scale_vars = sv,
    control = list(xtol_rel = 1e-10)))
  none <- fit("none")
  fits <- list(bhhh = fit("bhhh"), sd = fit("sd"))
  for (sv in names(fits)) {
    f <- fits[[sv]]
    expect_equal(f$loglik, none$loglik, tolerance = 1e-10, label = sv)
    expect_lt(max(abs(coef(f) - coef(none)) / none$se), 1e-4, label = sv)
    expect_identical(suppressMessages(vcov(f, type = "hessian")), f$vcov,
                     label = sv)
  }
  bh <- fits$bhhh
  expect_identical(log2(bh$param_scale), round(log2(bh$param_scale)))
  expect_identical(unname(bh$sX), c(1, 1))
  fs <- fits$sd
  expect_identical(unname(fs$param_scale), c(1 / unname(fs$sX), rep(1, 4)))
  expect_identical(unname(fs$param_shift), rep(0, 6))
})

test_that("a scaled fit's kernels evaluate the natural design", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  ns <- asNamespace("choicer")
  kernels <- c("mxl_loglik_gradient_parallel", "mxl_hessian_parallel",
               "mxl_bhhh_parallel", "mxl_cluster_meat_parallel")
  seen <- list()
  mocks <- lapply(stats::setNames(nm = kernels), function(name) {
    real <- get(name, envir = ns)
    function(..., X, W) {
      seen[[length(seen) + 1L]] <<- list(kernel = name, X = X, W = W)
      real(..., X = X, W = W)
    }
  })
  local_mocked_bindings(!!!mocks)
  dt <- sc_data(120L, 1L, 91L)
  dt[, grp := (match(id, unique(id)) - 1L) %% 7L]
  d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"))
  # The variance's kernels for each se_method; the gradient runs throughout.
  variance <- list(sandwich = c("mxl_hessian_parallel", "mxl_bhhh_parallel"),
                   cluster = c("mxl_hessian_parallel", "mxl_cluster_meat_parallel"))
  for (sv in c("sd", "bhhh")) {
    for (se in names(variance)) {
      what <- paste0("[", sv, ", ", se, "]")
      seen <- list()
      fit <- sc_fit(dt, scale_vars = sv, se_method = se,
                    cluster_col = if (se == "cluster") "grp")
      ran <- vapply(seen, `[[`, "", "kernel")
      expect_true(all(variance[[se]] %in% ran),
                  label = paste("the variance's kernels ran", what))
      expect_gt(sum(ran == "mxl_loglik_gradient_parallel"), 1L,
                label = paste("gradient calls", what))
      natural <- vapply(seen, function(s) {
        identical(s$X, d$X) && identical(s$W, d$W)
      }, TRUE)
      expect_true(all(natural), label = paste("natural design", what))
    }
  }
})

# Simulated nested logit data: two nests of 2 and 3 alternatives and a
# singleton nest holding j = 0 (an inside alternative here, the constants'
# reference), X in large units.
nl_sc_data <- function(N, seed) {
  sim <- simulate_nl_data(N = N, seed = seed)
  dt <- data.table::as.data.table(sim$data)
  dt[, X := X * 100]
  dt
}

test_that("nested logit scale_vars fits reach the unscaled optimum", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- nl_sc_data(600L, 21L)
  fit <- function(sv, ...) suppressMessages(run_nestlogit(
    dt, "id", "j", "choice", c("X", "W"), "nest", scale_vars = sv,
    control = list(xtol_rel = 1e-10, maxeval = 1000L), ...))
  bh <- fit("bhhh")
  sd <- fit("sd")
  # The unscaled optimizer polishing from the "bhhh" estimates: the optimum
  # both reach (unscaled from its cold start may stop short of it).
  ref <- fit("none", theta_init = coef(bh))
  for (f in list(bh, sd)) {
    expect_equal(f$loglik, ref$loglik, tolerance = 1e-10)
    expect_lt(max(abs(coef(f) - coef(ref)) / ref$se), 1e-3)
    expect_identical(suppressMessages(vcov(f, type = "hessian")), f$vcov)
  }
  expect_identical(log2(bh$param_scale), round(log2(bh$param_scale)))
  expect_identical(unname(bh$sX), c(1, 1))
  # Column scales reach the coefficients only.
  pm <- sd$param_map
  expect_identical(unname(sd$param_scale[pm$beta]), 1 / unname(sd$sX))
  expect_identical(unname(sd$param_scale[c(pm$lambda, pm$asc)]),
                   rep(1, length(c(pm$lambda, pm$asc))))
  expect_identical(unname(sd$param_shift), rep(0, length(coef(sd))))
  expect_identical(ref$scale_vars, "none")
  expect_identical(unname(ref$param_scale), rep(1, length(coef(ref))))
})

test_that("nested logit \"bhhh\" maps the start values and the nest parameters' bound", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- nl_sc_data(200L, 22L)
  th0 <- c(0.005, -0.3, 0.7, 0.4, 0, 0.1, -0.1, 0.2, 0.1)
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- list(theta_init = theta_init, lower = lower,
                  r = eval_f(theta_init))
    stop("captured")
  })
  expect_error(run_nestlogit(dt, "id", "j", "choice", c("X", "W"), "nest",
                             theta_init = th0, scale_vars = "bhhh"),
               "captured")
  d <- prepare_nl_data(dt, "id", "j", "choice", c("X", "W"), "nest")
  g <- nl_loglik_gradient_parallel(th0, d$X, d$alt_idx, d$choice_idx,
                                   d$nest_idx, d$M, d$weights, TRUE, FALSE,
                                   opg_diag = TRUE)
  m <- choicer:::.bhhh_scales(g$opg_diag)
  expect_true(all(m[3:4] != 1))  # the bound's mapping is visible
  expect_identical(seen$theta_init, th0 / m)
  expect_identical(seen$lower, c(-Inf, -Inf, 1e-16 / m[3:4], rep(-Inf, 5)))
  expect_identical(seen$r$objective, g$objective)
  expect_identical(seen$r$gradient, g$gradient * m)
})

test_that("a scaled fit's start values must have one value per parameter", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  # The default ("bhhh") checks them before its gradient pass.
  dm <- sc_data(60L, 1L, 71L)
  fit_x <- function(th) run_mxlogit(dm, "id", "alt", "choice", c("x1", "x2"),
                                    c("w1", "w2"), S = 10L, theta_init = th)
  expect_error(fit_x(rep(0, 8L)), "one value per parameter (7); got 8",
               fixed = TRUE)
  expect_error(fit_x(c(0, NA, 0, 0, 0, 0, 0)), "`theta_init` must be finite",
               fixed = TRUE)
})

test_that("a scaled nested logit's start values must have one value per parameter", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- nl_sc_data(200L, 22L)
  expect_error(suppressMessages(run_nestlogit(
    dt, "id", "j", "choice", c("X", "W"), "nest", scale_vars = "sd",
    theta_init = rep(0.1, 8L))), "one value per parameter (9); got 8",
    fixed = TRUE)
})

test_that("\"bhhh\" is the mixed and multinomial logits' default, \"none\" the nested logit's", {
  default <- function(f) eval(formals(f)$scale_vars)[1L]
  expect_identical(default(run_mxlogit), "bhhh")
  expect_identical(default(run_mnlogit), "bhhh")
  expect_identical(default(run_nestlogit), "none")
})
