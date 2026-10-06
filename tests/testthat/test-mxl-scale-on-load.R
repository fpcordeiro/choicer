# run_mxlogit(scale_vars = ) gives the estimation kernels the column scales,
# and the kernels divide each decision maker's rows of X and W by them as they
# load them, instead of the fit holding a sweep()ed copy of the design. Both
# are one IEEE division per value, so every kernel result, and so every fit,
# must be bit-identical to the sweep() path. Compared at one thread: with
# more, the threads' partial sums combine in varying order.

sol_bitwise <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE)

# Column scales as a fit passes them: log-normal W columns get 1.
sol_scales <- function(fx, kind) {
  col_sd <- function(M) vapply(seq_len(ncol(M)), function(k) stats::sd(M[, k]), 0)
  sX <- switch(kind,
               odd = c(1 / 3, 7321.5, 0.0029)[seq_len(ncol(fx$X))],
               sd = col_sd(fx$X))
  sW <- switch(kind,
               odd = c(0.0029, 3, 1 / 7)[seq_len(ncol(fx$W))],
               sd = col_sd(fx$W))
  sW[fx$rc_dist == 1L] <- 1
  list(sX = sX, sW = sW)
}

# The fixture with X and W divided by sweep(), as run_mxlogit() formed them
# before the kernels took the scales.
sol_swept <- function(fx, sc) {
  fx$X <- sweep(fx$X, 2, sc$sX, "/")
  fx$W <- sweep(fx$W, 2, sc$sW, "/")
  fx
}

test_that("the kernels a fit calls scale X and W on load exactly as sweep() did", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  for (fx in mxlp_build_cells()) {
    for (panel in c(TRUE, FALSE)) {
      f <- if (panel) fx else mxlp_cross_section(fx)
      for (kind in c("odd", "sd")) {
        sc <- sol_scales(f, kind)
        sw <- sol_swept(f, sc)
        for (gen in c(FALSE, TRUE)) {
          for (k in mxlp_weighted) {
            what <- sprintf("%s [%s, %s, %s scales, %s]", k, f$name,
                            if (panel) "panel" else "cross-section", kind,
                            if (gen) "generate" else "store")
            expect_true(sol_bitwise(
              mxlp_call(k, f, generate = gen, sX = sc$sX, sW = sc$sW),
              mxlp_call(k, sw, generate = gen)), label = what)
          }
        }
      }
    }
  }
})

test_that("scaling on load holds across draw batches and for NA in a log-normal column", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  cells <- mxlp_build_cells()
  # Cell 5 (normal + log-normal, alternative-level W) and cell 7 (two normal
  # columns and a log-normal one, row-aligned W).
  for (fx in cells[c(5L, 7L)]) {
    sc <- sol_scales(fx, "odd")
    sw <- sol_swept(fx, sc)
    for (k in mxlp_weighted) {
      expect_true(sol_bitwise(
        mxlp_call(k, fx, draw_batch = 5L, sX = sc$sX, sW = sc$sW),
        mxlp_call(k, sw, draw_batch = 5L)),
        label = sprintf("%s in batches of 5 draws [%s]", k, fx$name))
    }
    # NA in a log-normal column (scale 1): the results match sweep()'s. They
    # cannot show that the kernels divide that column by 1 rather than copy
    # it, since every result passes through arithmetic, which quiets R's
    # signaling-NaN NA either way; the code divides, as sweep() did.
    ln <- which(fx$rc_dist == 1L)[1L]
    fx$W[2L, ln] <- NA
    sw <- sol_swept(fx, sc)
    for (k in mxlp_weighted) {
      expect_true(sol_bitwise(mxlp_call(k, fx, sX = sc$sX, sW = sc$sW),
                              mxlp_call(k, sw)),
                  label = sprintf("%s with NA in a log-normal column [%s]", k,
                                  fx$name))
    }
  }
})

test_that("the scales of X and of W each apply to their own matrix", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  cells <- mxlp_build_cells()
  # Cell 3 (alternative-level W) and cell 7 (row-aligned W).
  for (fx in cells[c(3L, 7L)]) {
    sc <- sol_scales(fx, "odd")
    only_x <- fx
    only_x$X <- sweep(fx$X, 2, sc$sX, "/")
    only_w <- fx
    only_w$W <- sweep(fx$W, 2, sc$sW, "/")
    for (k in mxlp_weighted) {
      expect_true(sol_bitwise(mxlp_call(k, fx, sX = sc$sX), mxlp_call(k, only_x)),
                  label = sprintf("%s with sX alone [%s]", k, fx$name))
      expect_true(sol_bitwise(mxlp_call(k, fx, sW = sc$sW), mxlp_call(k, only_w)),
                  label = sprintf("%s with sW alone [%s]", k, fx$name))
    }
  }
})

test_that("the kernels check the column scales after their other inputs", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  fx <- mxlp_build_cells()[[1L]]  # K_x = 2, K_w = 2
  for (k in mxlp_weighted) {
    expect_error(mxlp_call(k, fx, sX = 1),
                 "sX must hold one scale per column of X (2); got 1.", fixed = TRUE)
    expect_error(mxlp_call(k, fx, sW = c(1, 2, 3)),
                 "sW must hold one scale per column of W (2); got 3.", fixed = TRUE)
  }
  g <- function(...) mxlp_call("gradient", fx, ...)
  expect_error(g(sX = numeric(0)),
               "sX must hold one scale per column of X (2); got 0.", fixed = TRUE)
  expect_error(g(sX = c(1, NA)),
               "sX must be finite and positive; sX[2] is NA.", fixed = TRUE)
  expect_error(g(sX = c(NaN, 1)),
               "sX must be finite and positive; sX[1] is NaN.", fixed = TRUE)
  expect_error(g(sX = c(Inf, 1)),
               "sX must be finite and positive; sX[1] is Inf.", fixed = TRUE)
  expect_error(g(sW = c(0, 1)),
               "sW must be finite and positive; sW[1] is 0.", fixed = TRUE)
  expect_error(g(sW = c(1, -2.5)),
               "sW must be finite and positive; sW[2] is -2.5.", fixed = TRUE)
  expect_error(g(sX = c("a", "b")), "Not compatible with requested type")
  # Integer scales convert exactly; NULL (the default) reads the matrices as
  # they are.
  expect_true(sol_bitwise(g(sX = c(2L, 4L)), g(sX = c(2, 4))))
  expect_true(sol_bitwise(
    mxl_loglik_gradient_parallel(
      fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights,
      fx$eta, fx$rc_dist, fx$rc_correlation, fx$rc_mean, fx$use_asc,
      fx$include_outside_option, Ti = fx$Ti, sX = NULL, sW = NULL),
    g()))
  # An earlier check keeps its precedence over a bad scale.
  expect_error(g(eta = fx$eta[1L, , , drop = FALSE], sX = 1),
               "eta_draws 1st dimension (1) does not match K_w (2)", fixed = TRUE)
})

# --- run_mxlogit() end to end ------------------------------------------------

sol_kernels <- c("mxl_loglik_gradient_parallel", "mxl_hessian_parallel",
                 "mxl_bhhh_parallel", "mxl_cluster_meat_parallel")

# Evaluate `expr` with the four kernels a fit calls replaced by recorders:
# each records the design and scales it received and calls route(f), f being
# the real kernel, taken before the bindings are replaced. Returns the value
# and the record.
sol_with_kernels <- function(expr, route) {
  ns <- asNamespace("choicer")
  seen <- list()
  mocks <- lapply(stats::setNames(nm = sol_kernels), function(name) {
    real <- get(name, envir = ns)
    f <- route(real)
    function(..., X, W, sX = NULL, sW = NULL) {
      seen[[length(seen) + 1L]] <<- list(kernel = name, X = X, W = W,
                                          sX = sX, sW = sW)
      f(..., X = X, W = W, sX = sX, sW = sW)
    }
  })
  local_mocked_bindings(!!!mocks)
  value <- expr
  list(value = value, seen = seen)
}

sol_new <- function(f) f

sol_fit <- function(args, route) {
  r <- sol_with_kernels(suppressMessages(do.call("run_mxlogit", args)), route)
  r$value$optimizer$elapsed_time <- NULL
  r
}

# Badly scaled columns, a positive log-normal regressor (w2), and weights and
# cluster labels constant within decision maker.
sol_data <- function(N, T, ioo, seed) {
  sim <- simulate_mxl_data(N = N, J = 4L, T = T, beta = c(0.8, -0.5),
                           Sigma = diag(c(1, 0.6)), outside_option = ioo,
                           vary_choice_set = TRUE, seed = seed)
  dt <- data.table::as.data.table(sim$data)
  # The simulated outside option is an explicit alt = 0 row; the implicit
  # outside option needs those rows dropped.
  if (ioo) dt <- dt[alt != 0]
  dt[, `:=`(x1 = x1 * 100, x2 = x2 / 40, w1 = w1 * 25, w2 = abs(w2) + 0.1)]
  unit <- if (T > 1L) dt$pid else dt$id
  u <- match(unit, unique(unit))
  dt[, `:=`(grp = (u - 1L) %% 7L, w = 0.5 + (u %% 5L) / 4)]
  dt
}

sol_args <- function(dt, panel, ...) {
  c(list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
         person_col = if (panel) "pid", S = 20L, seed = 3L,
         control = list(maxeval = 300L)),
    list(...))
}

# The scaling is the optimizer's coordinates: every kernel call a scaled fit
# makes, the objective's and the variance's, gets the natural design and no
# column scales.
test_that("a scaled fit's kernels evaluate the natural design", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  cases <- list(
    list(panel = FALSE, ioo = FALSE, sv = "sd", draws = "store",
         rc_dist = c(0L, 0L), rc_mean = FALSE, corr = TRUE, weights = FALSE,
         se = c("hessian", "bhhh", "sandwich", "cluster")),
    list(panel = FALSE, ioo = TRUE, sv = "mad", draws = "generate",
         rc_dist = c(0L, 1L), rc_mean = TRUE, corr = FALSE, weights = FALSE,
         se = "hessian"),
    list(panel = TRUE, ioo = TRUE, sv = "iqr", draws = "store",
         rc_dist = c(1L, 1L), rc_mean = TRUE, corr = FALSE, weights = FALSE,
         se = c("hessian", "cluster")),
    list(panel = TRUE, ioo = FALSE, sv = "bhhh", draws = "generate",
         rc_dist = c(0L, 0L), rc_mean = FALSE, corr = TRUE, weights = TRUE,
         se = "sandwich"))
  eager <- list(hessian = "mxl_hessian_parallel", bhhh = "mxl_bhhh_parallel",
                sandwich = c("mxl_hessian_parallel", "mxl_bhhh_parallel"),
                cluster = c("mxl_hessian_parallel", "mxl_cluster_meat_parallel"))
  for (cs in cases) {
    dt <- sol_data(if (cs$panel) 60L else 180L, if (cs$panel) 3L else 1L,
                   cs$ioo, 300L + cs$panel + 2L * cs$ioo)
    for (se in cs$se) {
      what <- sprintf("[%s, %s, %s, rc_dist = %s, se_method = %s]",
                      if (cs$panel) "panel" else "cross-section", cs$sv,
                      cs$draws, paste(cs$rc_dist, collapse = ""), se)
      args <- sol_args(dt, cs$panel, rc_dist = cs$rc_dist,
                       rc_mean = cs$rc_mean, rc_correlation = cs$corr,
                       include_outside_option = cs$ioo, draws = cs$draws,
                       scale_vars = cs$sv, se_method = se,
                       cluster_col = if (se == "cluster") "grp",
                       weights_col = if (cs$weights) "w")
      new <- sol_fit(args, sol_new)
      kernels <- vapply(new$seen, `[[`, "", "kernel")
      expect_true(all(eager[[se]] %in% kernels),
                  label = paste("the variance kernels ran", what))
      expect_gt(sum(kernels == "mxl_loglik_gradient_parallel"), 1L)
      fit <- new$value
      # The optimizer's coordinates are scaled, the data are not: every call
      # gets the natural design and no column scales.
      got <- vapply(new$seen, function(s) {
        identical(s$X, fit$data$X) && identical(s$W, fit$data$W) &&
          is.null(s$sX) && is.null(s$sW)
      }, TRUE)
      expect_true(all(got),
                  label = paste("every call got the natural design and no scales", what))
    }
  }
})

test_that("unscaled fits and post-hoc variances pass no scales", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  dt <- sol_data(180L, 1L, FALSE, 310L)
  args <- sol_args(dt, FALSE, rc_correlation = TRUE)
  no_scales <- function(seen) {
    all(vapply(seen, function(s) is.null(s$sX) && is.null(s$sW), TRUE))
  }
  none <- sol_fit(c(args, list(scale_vars = "none")), sol_new)
  expect_true(no_scales(none$seen))
  scaled <- sol_fit(c(args, list(scale_vars = "sd")), sol_new)
  expect_true(no_scales(scaled$seen))
  # Post-hoc variances work with natural coefficients and data.
  fit <- scaled$value
  cl <- dt[!duplicated(id), stats::setNames(grp, id)]
  post <- sol_with_kernels(suppressMessages(list(
    hessian = vcov(fit, type = "hessian"),
    bhhh = vcov(fit, type = "bhhh"),
    robust = vcov(fit, type = "robust"),
    cluster = vcov(fit, type = "cluster", cluster = cl))), sol_new)
  expect_gte(length(post$seen), 4L)
  expect_true(no_scales(post$seen))
})
