# run_mxlogit(concentrate_out_asc = "blp"): the constants solve the share
# equations on the estimation draws (the sample's own shares, supplied
# targets, or a choice-based sample's Q) and the other parameters maximize
# the likelihood given them. The kernels and the variance's algebra are
# test-mxl-blp-jacobian.R's and test-mxl-blp-vcov.R's; here, the estimator.

# A design with an outside option (or not), a cross-section or a panel, one
# or two markets (each decision maker in one): utility
#   V = x1 + (-1 + 0.8 e_n) w1 + delta_j + Gumbel,  delta = seq(-1, 1),
# x1 ~ N(0, 1), w1 ~ N(0, 1.5^2), e_n ~ N(0, 1) per decision maker, choice
# sets at least half of the market's alternatives (two at least).
blpc_data <- function(seed, N = 400L, J = 4L, T = 1L, outside = TRUE,
                      markets = 1L, weights = FALSE) {
  set.seed(seed)
  n_sit <- N * T
  pid <- rep(seq_len(N), each = T)
  mk_alt <- rep(seq_len(markets), each = J %/% markets)
  mk_sit <- rep(sample.int(markets, N, replace = TRUE), each = T)
  alts <- lapply(seq_len(n_sit), function(t) {
    pool <- which(mk_alt == mk_sit[t])
    k <- max(2L, ceiling(length(pool) / 2)):length(pool)
    k <- k[sample.int(length(k), 1L)]
    sort(pool[sample.int(length(pool), k)])
  })
  dt <- data.table::data.table(id = rep(seq_len(n_sit), lengths(alts)),
                               alt = unlist(alts))
  dt[, pid := pid[id]]
  dt[, x1 := stats::rnorm(.N)]
  dt[, w1 := stats::rnorm(.N, sd = 1.5)]
  e <- stats::rnorm(N)
  delta <- seq(-1, 1, length.out = J)
  dt[, p := exp(x1 + (-1 + 0.8 * e[pid]) * w1 + delta[alt])]
  dt[, p := p / (sum(p) + outside), by = id]
  r <- stats::runif(n_sit)
  dt[, choice := {
    cp <- cumsum(p)
    as.integer(r[id[1L]] <= cp & r[id[1L]] > cp - p)
  }, by = id]
  if (weights) dt[, wt := stats::runif(1L, 0.5, 2), by = pid]
  dt[, p := NULL]
  dt[]
}

# The fit, without the advice on weights with se_method "hessian" or
# "bhhh" (these tests' weights are frequencies); other warnings pass
blpc_fit <- function(dt, ..., outside = TRUE, panel = FALSE, S = 20L,
                     covariate_cols = c("x1", "w1"), random_var_cols = "w1") {
  withCallingHandlers(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = covariate_cols, random_var_cols = random_var_cols, S = S,
      include_outside_option = outside, outside_opt_label = if (outside) 0L,
      person_col = if (panel) "pid", ...)),
    warning = function(w) {
      if (grepl("Non-uniform weights detected", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
}

# The fit's predicted counts on its estimation draws, every inside
# alternative's
blpc_counts <- function(fit) {
  d <- fit$data
  gp <- choicer:::.mxl_gen_params(fit$draws_info,
                                  N = length(choicer:::.unit_first(d)))
  as.numeric(choicer:::mxl_blp_jacobian_parallel(
    coef(fit), d$X, d$W, d$alt_idx, d$choice_idx, d$M, d$weights, gp$eta_draws,
    fit$rc_dist, fit$rc_correlation, fit$rc_mean, TRUE,
    isTRUE(fit$include_outside_option), gp$gen_seed, gp$gen_scramble,
    gp$gen_S, Ti = d$Ti)$counts)
}

# The reference's fixture (helper-blp-reference.R) for a fit: its data and
# estimation draws at its estimate, and each likelihood unit's cluster
blpc_fx <- function(fit) {
  d <- fit$data
  first_u <- choicer:::.unit_first(d)
  gp <- choicer:::.mxl_gen_params(fit$draws_info, N = length(first_u))
  outside <- isTRUE(fit$include_outside_option)
  list(theta = unname(coef(fit)), X = d$X, W = d$W, alt_idx = d$alt_idx,
       choice_idx = d$choice_idx, M = d$M, weights = d$weights, Ti = d$Ti,
       eta = gp$eta_draws, rc_dist = fit$rc_dist,
       rc_correlation = fit$rc_correlation, rc_mean = fit$rc_mean, use_asc = TRUE,
       include_outside_option = outside, K_w = ncol(d$W),
       J = length(fit$param_map$asc) + !outside,
       cluster = if (!is.null(d$cluster)) d$cluster[first_u])
}

test_that("the ML fit's own shares as targets give the ML estimates", {
  skip_on_cran()  # two fits per design
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  cfgs <- list(list(T = 1L, outside = TRUE), list(T = 3L, outside = TRUE),
               list(T = 1L, outside = FALSE), list(T = 3L, outside = TRUE, weights = TRUE))
  for (cfg in cfgs) {
    w <- isTRUE(cfg$weights)
    dt <- blpc_data(601, N = if (cfg$T > 1L) 150L else 450L, T = cfg$T,
                    outside = cfg$outside, weights = w)
    args <- list(outside = cfg$outside, panel = cfg$T > 1L,
                 weights_col = if (w) "wt")
    ml <- do.call(blpc_fit, c(list(dt, polish = "newton",
                                   polish_control = list(tol_decrement = 1e-14)), args))
    # the ML fit's shares on its estimation draws (one market: within it)
    s <- blpc_counts(ml) / sum(ml$data$weights)
    fb <- do.call(blpc_fit, c(list(dt, concentrate_out_asc = "blp", target_shares = s),
                              args))
    lab <- paste0("T = ", cfg$T, if (!cfg$outside) ", no outside option",
                  if (w) ", weights")
    expect_lt(max(abs(coef(fb) - coef(ml)) / ml$se), 1e-4, label = lab)
    expect_equal(fb$loglik, ml$loglik, tolerance = 1e-9, label = lab)
    expect_identical(fb$share_targets$source, "supplied", label = lab)
    # the fit matches its targets
    expect_equal(blpc_counts(fb) / sum(fb$data$weights), s, tolerance = 1e-10,
                 label = lab)
  }
})

test_that("correlated, located and log-normal coefficients: ML's shares give ML", {
  # (with rc_correlation the parameter map's positions are doubles)
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(619, N = 400L)
  dt[, nw1 := -w1]
  cfgs <- list(
    correlated = list(rc_correlation = TRUE),
    `log-normal, with its location` = list(covariate_cols = "x1",
                                           random_var_cols = "nw1",
                                           rc_mean = TRUE, rc_dist = 1L))
  for (nm in names(cfgs)) {
    ml <- do.call(blpc_fit, c(list(dt, polish = "newton",
                                   polish_control = list(tol_decrement = 1e-14)),
                              cfgs[[nm]]))
    s <- blpc_counts(ml) / sum(ml$data$weights)
    fb <- do.call(blpc_fit, c(list(dt, concentrate_out_asc = "blp", target_shares = s),
                              cfgs[[nm]]))
    expect_lt(max(abs(coef(fb) - coef(ml)) / ml$se), 1e-4, label = nm)
    expect_equal(fb$loglik, ml$loglik, tolerance = 1e-9, label = nm)
  }
})

test_that("the sample's own shares are matched, and the fit is certified", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(602, N = 400L)
  fb <- blpc_fit(dt, concentrate_out_asc = "blp")
  d <- fb$data
  first <- c(0, cumsum(d$M))[seq_along(d$M)]
  ch <- d$choice_idx > 0L
  obs <- tabulate(d$alt_idx[first[ch] + d$choice_idx[ch]], 4L)
  expect_equal(blpc_counts(fb), obs, tolerance = 1e-10)
  expect_identical(fb$share_targets$source, "sample")
  rec <- fb$optimizer$concentration
  expect_identical(rec$method, "blp")
  expect_lt(rec$max_residual, 1e-11)
  expect_lt(rec$decrement, 1e-10)
  expect_identical(names(rec$lambda), names(coef(fb))[fb$param_map$asc])
  expect_identical(rec$solves, rec$jacobians)
  # predict() and blp() on the fit return the targets and the constants
  # (a cross-section, store-mode draws)
  S <- fb$share_targets$engine
  expect_equal(as.numeric(predict(fb, type = "shares")), S, tolerance = 1e-10)
  expect_equal(as.numeric(blp(fb, target_shares = S)),
               unname(coef(fb)[fb$param_map$asc]), tolerance = 1e-8)
  expect_output(print(summary(fb)), "share-constrained")
  # the log-likelihood is the likelihood's at the estimate (not the
  # optimizer's corrected objective), the gradient kernel's at one thread
  gp <- choicer:::.mxl_gen_params(fb$draws_info, N = length(d$M))
  ll <- choicer:::mxl_loglik_gradient_parallel(
    coef(fb), d$X, d$W, d$alt_idx, d$choice_idx, d$M, d$weights, gp$eta_draws,
    fb$rc_dist, fb$rc_correlation, fb$rc_mean, TRUE, TRUE)$objective
  expect_identical(fb$loglik, -ll)
  # without the data kept: the same fit, variance and targets
  fk <- blpc_fit(dt, concentrate_out_asc = "blp", keep_data = FALSE)
  expect_null(fk$data)
  expect_identical(coef(fk), coef(fb))
  expect_identical(fk$vcov, fb$vcov)
  expect_identical(fk$share_targets, fb$share_targets)
  expect_output(print(summary(fk)), "share-constrained")
})

test_that("without heterogeneity the estimates are the logit's", {
  # W = 0: the share equations of the sample's shares are the logit's
  # first-order conditions for the constants, and L moves nothing
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(603, N = 500L)
  d <- suppressMessages(prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "w1"),
                                         "w1", include_outside_option = TRUE,
                                         outside_opt_label = 0L))
  d$W[] <- 0
  fb <- suppressMessages(suppressWarnings(run_mxlogit(
    input_data = d, eta_draws = get_halton_normals(20L, length(d$M), 1L),
    concentrate_out_asc = "blp", scale_vars = "none", theta_init = rep(0, 7L),
    control = list(xtol_rel = 1e-12, ftol_rel = 1e-14, maxeval = 2000L))))
  mn <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "w1"), include_outside_option = TRUE,
    outside_opt_label = 0L))
  keep <- c(fb$param_map$beta, fb$param_map$asc)
  expect_equal(unname(coef(fb)[keep]), unname(coef(mn)), tolerance = 1e-5)
})

test_that("the implicit gradient and the corrected objective", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (T in c(1L, 3L)) {
    dt <- blpc_data(604, N = if (T > 1L) 100L else 300L, T = T)
    # targets the micro data disagree with, so that the constants' likelihood
    # score, and the plain objective's first-order error, are not small
    fb <- blpc_fit(dt, concentrate_out_asc = "blp", panel = T > 1L,
                   target_shares = c(0.10, 0.14, 0.20, 0.26))
    d <- fb$data
    U <- length(choicer:::.unit_first(d))
    eta <- get_halton_normals(20L, U, 1L)
    gen <- list(seed = -1L, scramble = 1L, S = 0L)
    setup <- function(tol) {
      choicer:::.mxl_blp_setup(d, TRUE, eta, gen, fb$rc_dist, fb$rc_correlation,
                               fb$rc_mean, fb$param_map, length(coef(fb)),
                               fb$share_targets,
                               choicer:::.concentrate_control(list(tol_shares = tol),
                                                              "blp"))
    }
    exact <- setup(1e-14)
    c_idx <- seq_len(length(coef(fb)) - length(fb$param_map$asc))
    th <- coef(fb)[c_idx] + c(0.05, -0.04, 0.03)   # away from the optimum
    d0 <- unname(coef(fb)[fb$param_map$asc])
    s <- exact$solve(th, d0, 1000L)
    expect_true(s$converged)
    num <- numDeriv::grad(function(x) exact$solve(x, d0, 1000L)$value, th)
    expect_equal(s$gradient[c_idx], num, tolerance = 1e-6, label = paste("T =", T))
    # a loose solve: the corrected objective's error is second order, the
    # plain objective's first
    loose <- setup(1e-3)$solve(th, d0 - 0.2, 1000L)
    plain <- exact$jac_f(c(th, loose$delta))$objective
    expect_gt(abs(plain - s$value), 1e-6)
    expect_lt(abs(loose$value - s$value), 1e-2 * abs(plain - s$value),
              label = paste("T =", T))
  }
})

test_that("supplied targets: recovery and the constants' precision", {
  skip_on_cran()
  set_omp_threads(2L)
  # the DGP's population shares within each market (a large simulation)
  pop <- blpc_data(605, N = 100000L, T = 1L)
  set.seed(1)
  pop[, e := stats::rnorm(1L), by = id]
  pop[, p := exp(x1 + (-1 + 0.8 * e) * w1 + seq(-1, 1, length.out = 4L)[alt])]
  pop[, p := p / (sum(p) + 1), by = id]
  s_pop <- pop[, .(s = sum(p) / data.table::uniqueN(pop$id)), by = alt][order(alt)]$s
  dt <- blpc_data(606, N = 2000L)
  truth <- c(1, -1, log(0.8), seq(-1, 1, length.out = 4L))
  fs <- blpc_fit(dt, concentrate_out_asc = "blp")
  fp <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = s_pop)
  for (f in list(fs, fp)) {
    expect_true(all(abs(coef(f) - truth) / f$se < 4))
  }
  # exact targets carry information the sample's own shares do not
  asc <- fp$param_map$asc
  expect_true(all(fp$se[asc] < fs$se[asc]))
  # a target sample of n choices adds its covariance
  V_s <- (diag(s_pop) - s_pop %o% s_pop) / 5000
  fv <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = s_pop,
                 target_vcov = V_s)
  expect_equal(coef(fv), coef(fp), tolerance = 1e-8)   # (two threads)
  expect_true(all(fv$se[asc] > fp$se[asc]))
})

test_that("the variance by every route, at one thread", {
  skip_on_cran()  # eight fits
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(607, N = 120L, T = 3L, markets = 2L, J = 4L, weights = TRUE)
  # clusters of decision makers, within markets (alternatives 1-2 and 3-4)
  dt[, cl := (pid %% 7L) * 10L + (alt[1L] > 2L), by = id]
  s_m <- c(0.3, 0.4, 0.35, 0.3)   # within each market's two alternatives
  for (tg in list(NULL, s_m)) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      fit <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE,
                      weights_col = "wt", cluster_col = "cl", se_method = se,
                      target_shares = tg)
      lab <- paste(se, if (is.null(tg)) "sample" else "supplied")
      expect_false(is.null(fit$vcov), label = lab)
      expect_identical(unname(vcov(fit, type = choicer:::.vcov_type(se))),
                       unname(fit$vcov), label = lab)
      lazy <- fit
      lazy["vcov"] <- list(NULL)
      expect_identical(unname(vcov(lazy)), unname(fit$vcov), label = lab)
    }
  }
})

test_that("the fit's variance is the stacked equations' (a numerical reference)", {
  # blpr_vcov_bf(): G^{-1} Omega G^{-T} with G by numDeriv, from the fit's
  # data, draws, targets and clusters, for each kind of target
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  check <- function(fit, type, tg = NULL, Vs = NULL, W_free = NULL, lab) {
    fx <- blpc_fx(fit)
    ref <- blpr_vcov_bf(fx, type, tg, cluster = fx$cluster, Vs = Vs, W_free = W_free)
    expect_equal(unname(fit$vcov), ref, tolerance = 1e-6, label = lab)
  }
  # the sample's shares: a cross-section with weights and clusters
  dt <- blpc_data(616, N = 90L, weights = TRUE)
  dt[, cl := id %% 9L]
  for (se in c("hessian", "sandwich", "cluster")) {
    fit <- blpc_fit(dt, concentrate_out_asc = "blp", weights_col = "wt",
                    cluster_col = "cl", se_method = se, S = 10L)
    check(fit, se, lab = paste("sample,", se))
  }
  # supplied targets with their covariance, without an outside option (the
  # base's row and column dropped), a panel with weights (the sandwich
  # squares them)
  dp <- blpc_data(617, N = 40L, T = 3L, outside = FALSE, weights = TRUE)
  s <- c(0.15, 0.2, 0.3, 0.35)
  V_s <- (diag(s) - s %o% s) / 2000
  for (se in c("hessian", "sandwich")) {
    fit <- blpc_fit(dp, concentrate_out_asc = "blp", outside = FALSE, panel = TRUE,
                    weights_col = "wt", target_shares = s, target_vcov = V_s,
                    se_method = se, S = 10L)
    W <- sum(fit$data$weights)
    check(fit, se, list(s = s[-1L], market_of_alt = rep(1L, 4L)), V_s[-1L, -1L],
          rep(W, 3L), lab = paste("supplied, no outside option,", se))
  }
  # the covariance given to four digits: singular by construction, its
  # smallest eigenvalue then -8e-5 of its largest, but not on the block the
  # variance takes (the base's row and column dropped)
  fr <- blpc_fit(dp, concentrate_out_asc = "blp", outside = FALSE, panel = TRUE,
                 weights_col = "wt", target_shares = s, target_vcov = signif(V_s, 4L),
                 se_method = "sandwich", S = 10L)
  expect_equal(fr$vcov, fit$vcov, tolerance = 1e-3)
  # two markets' targets, weights, clusters within markets
  d2 <- blpc_data(618, N = 60L, T = 2L, markets = 2L, weights = TRUE)
  d2[, cl := (pid %% 5L) * 10L + (alt[1L] > 2L), by = id]
  s_m <- c(0.3, 0.4, 0.35, 0.3)
  for (se in c("hessian", "cluster")) {
    fit <- blpc_fit(d2, concentrate_out_asc = "blp", panel = TRUE, weights_col = "wt",
                    cluster_col = "cl", se_method = se, target_shares = s_m, S = 10L)
    check(fit, se, list(s = s_m, market_of_alt = fit$share_targets$market),
          lab = paste("two markets,", se))
  }
})

test_that("a T = 1 panel is the cross-section, and the pathways agree", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(608, N = 300L)
  f1 <- blpc_fit(dt, concentrate_out_asc = "blp")
  fp <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE)
  expect_equal(coef(fp), coef(f1), tolerance = 1e-10)
  expect_equal(fp$vcov, f1$vcov, tolerance = 1e-10)
  d <- suppressMessages(prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "w1"),
                                         "w1", include_outside_option = TRUE,
                                         outside_opt_label = 0L))
  fa <- suppressMessages(run_mxlogit(
    input_data = d, eta_draws = get_halton_normals(20L, length(d$M), 1L),
    concentrate_out_asc = "blp"))
  expect_identical(coef(fa), coef(f1))
  expect_identical(fa$vcov, f1$vcov)
})

test_that("the checks before the draws", {
  dt <- blpc_data(609, N = 200L)
  n_draws <- 0L
  local_mocked_bindings(get_halton_normals = function(...) {
    n_draws <<- n_draws + 1L
    .halton_cube(...)
  })
  blp_err <- function(..., regexp) {
    expect_error(blpc_fit(dt, concentrate_out_asc = "blp", ...), regexp, fixed = TRUE)
  }
  blp_err(use_asc = FALSE, regexp = "so it needs them: use_asc = TRUE")
  blp_err(lower = c(rep(-Inf, 3L), -5, rep(-Inf, 3L)),
          regexp = "without bounds, but `lower`/`upper` bound")
  blp_err(polish = "newton", regexp = "polish = \"newton\" would move them")
  blp_err(concentrate_control = list(max_dense = 2L),
          regexp = "up to concentrate_control$max_dense = 2 constants; this model has 4")
  blp_err(target_shares = c(0.2, 0.2, 0.2),
          regexp = "`target_shares` must have one share per inside alternative (4)")
  blp_err(target_shares = c(a = 0.2, b = 0.2, c = 0.2, d = 0.2),
          regexp = "`target_shares` names must be the inside alternatives' labels")
  blp_err(target_shares = c(0.2, 0.2, -0.1, 0.2),
          regexp = "`target_shares` must be finite and positive (entry 3 is -0.1)")
  blp_err(target_shares = c(`4` = 0.2, `3` = -0.1, `2` = 0.2, `1` = 0.2),
          regexp = "`target_shares` must be finite and positive (entry 2, \"3\", is -0.1)")
  blp_err(target_shares = c(`0` = 0.2, `1` = 0.2, `2` = 0.2, `3` = 0.2, `4` = 0.2),
          regexp = paste0("`target_shares` names must be the inside alternatives' ",
                          "labels, each once (4; the outside option's share follows"))
  blp_err(target_shares = c(0.3, 0.3, 0.3, 0.3),
          regexp = paste0("target_shares, as targets, fail the share inversion's ",
                          "checks. The target shares leave no share to the outside option"))
  # above the weight of the situations that offer it, by label
  blp_err(target_shares = c(0.9, 0.02, 0.02, 0.02),
          regexp = "Alternative 1's target share (0.9) is not below")
  blp_err(target_shares = c(0.9, 0.02, 0.02, 0.02),
          regexp = "(In its numbering, 1 is \"1\".)")
  blp_err(target_vcov = diag(4),
          regexp = "`target_vcov` is the sampling covariance of supplied `target_shares`")
  blp_err(target_shares = rep(0.2, 4L), target_vcov = diag(3),
          regexp = "`target_vcov` must be a finite symmetric 4 x 4 matrix")
  not_cov <- paste0("`target_vcov` must be a covariance, positive semi-definite to ",
                    "within rounding (checked on the inside alternatives")
  blp_err(target_shares = rep(0.2, 4L), target_vcov = diag(c(1, 1, 1, -1)),
          regexp = not_cov)
  v_cor <- diag(4)
  v_cor[1L, 2L] <- v_cor[2L, 1L] <- 2   # a correlation of 2
  blp_err(target_shares = rep(0.2, 4L), target_vcov = v_cor, regexp = not_cov)
  v_zero <- diag(c(0, 1, 1, 1))
  v_zero[1L, 2L] <- v_zero[2L, 1L] <- 0.1   # covariance without variance
  blp_err(target_shares = rep(0.2, 4L), target_vcov = v_zero, regexp = not_cov)
  blp_err(target_shares = rep(0.2, 4L),
          target_vcov = matrix(0, 4L, 4L, dimnames = list(letters[1:4], letters[1:4])),
          regexp = "`target_vcov`'s row and column names, when it has them, must both be")
  blp_err(target_shares = c(0.3, 0.2, 0.2, 0.2, 0.2),
          regexp = "`target_shares` with the outside option's share first must sum to 1")
  # the sample's shares: an alternative never chosen
  dz <- data.table::copy(dt)[alt != 4L | choice == 0L]
  expect_error(blpc_fit(dz, concentrate_out_asc = "blp"),
               class = "choicer_unbounded_constants")
  # ... the outside option never chosen in one of two markets (alternatives 3
  # and 4), though chosen in the other
  d3 <- blpc_data(615, N = 200L, markets = 2L)
  d3[, choice := if (alt[1L] > 2L && sum(choice) == 0L) as.integer(alt == max(alt)) else choice,
     by = id]
  e3 <- expect_error(blpc_fit(d3, concentrate_out_asc = "blp"),
                     class = "choicer_unbounded_constants")
  expect_match(conditionMessage(e3), "no one chooses it in the market of 3, 4", fixed = TRUE)
  # ... without an outside option, an alternative chosen only where offered
  # alone (once, and in two situations that list it twice)
  d4 <- blpc_data(620, N = 200L, outside = FALSE)
  d4[, choice := if (any(alt == 4L & choice == 1L)) as.integer(alt == min(alt)) else choice,
     by = id]
  n4 <- max(d4$id)
  d4 <- rbind(d4, data.table::data.table(id = n4 + c(1L, 2L, 2L, 3L, 3L), alt = 4L,
                                         pid = n4 + c(1L, 2L, 2L, 3L, 3L),
                                         x1 = c(0.1, 0.2, 0.3, -0.1, 0.4),
                                         w1 = c(0.2, -0.3, 0.1, 0.5, 0.2),
                                         choice = c(1L, 1L, 0L, 1L, 0L)))
  e4 <- expect_error(blpc_fit(d4, concentrate_out_asc = "blp", outside = FALSE),
                     class = "choicer_unbounded_constants")
  expect_match(conditionMessage(e4), "4 is chosen only where offered alone", fixed = TRUE)
  # two markets without an outside option
  d2 <- blpc_data(610, N = 200L, outside = FALSE, markets = 2L)
  expect_error(blpc_fit(d2, concentrate_out_asc = "blp", outside = FALSE),
               "split them into 2 groups that never meet", fixed = TRUE)
  expect_identical(n_draws, 0L)
  # the start's solve, after the draws, with too few evaluations of the shares
  e <- expect_error(blpc_fit(dt, concentrate_out_asc = "blp",
                             concentrate_control = list(max_iter = 1L)),
                    class = "choicer_blp_unconverged")
  expect_match(conditionMessage(e),
               "at the start values, the share inversion stopped (max_iter)", fixed = TRUE)
  expect_identical(n_draws, 1L)
})

test_that("solves failing after the best point warn, and the best point solved stands", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(621, N = 200L)
  # every solve after the third (the start's and the optimizer's first two)
  # fails
  real_setup <- .mxl_blp_setup
  n_solves <- 0L
  local_mocked_bindings(.mxl_blp_setup = function(...) {
    s <- real_setup(...)
    solve <- s$solve
    s$solve <- function(theta_c, delta0, max_iter, best = NULL) {
      n_solves <<- n_solves + 1L
      if (n_solves <= 3L) return(solve(theta_c, delta0, max_iter, best))
      list(converged = FALSE, value = NA_real_, delta = delta0, gradient = NULL,
           iterations = 0L, passes = 0L, residual = NA_real_,
           decrement = NA_real_, message = NULL, error = "a mocked failure")
    }
    s
  })
  # (the solves it saw as overflows can also leave a failure code, its own warning)
  w <- suppressWarnings(
    expect_warning(fit <- blpc_fit(dt, concentrate_out_asc = "blp"),
                   class = "choicer_blp_unconverged"),
    classes = "choicer_optimizer_failed")
  expect_match(conditionMessage(w), paste0(
    "of the inner solves after the optimizer's best point did not converge (it ",
    "saw them as overflows, so its success code can mislead; the last error: a ",
    "mocked failure)"), fixed = TRUE)
  expect_gt(fit$optimizer$concentration$failed, 0L)
  # the constants solve the share equations at the returned point
  expect_equal(blpc_counts(fit), fit$share_targets$engine[-1L] * sum(fit$data$weights),
               tolerance = 1e-10)
  expect_false(is.null(fit$vcov))
})

test_that("the constants are solved where the optimizer stops, or the best point stands", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(622, N = 200L)
  # a custom optimizer that evaluates its start (the cache) and returns a
  # point downhill that it never evaluated: the constants are solved there
  opt_away <- function(theta_init, eval_f, ...) {
    e <- eval_f(theta_init)
    list(par = theta_init - 0.05 * e$gradient / sqrt(sum(e$gradient^2)),
         value = e$objective, convergence = 0L)
  }
  ll <- function(fit) {
    d <- fit$data
    gp <- choicer:::.mxl_gen_params(fit$draws_info, N = length(d$M))
    -choicer:::mxl_loglik_gradient_parallel(
      coef(fit), d$X, d$W, d$alt_idx, d$choice_idx, d$M, d$weights, gp$eta_draws,
      fit$rc_dist, fit$rc_correlation, fit$rc_mean, TRUE, TRUE)$objective
  }
  fit <- blpc_fit(dt, concentrate_out_asc = "blp", optimizer = opt_away)
  rec <- fit$optimizer$concentration
  expect_identical(rec$evaluations, 1L)
  expect_identical(rec$solves, 2L)   # at the start and at the point returned
  expect_identical(fit$loglik, ll(fit))
  expect_equal(blpc_counts(fit), fit$share_targets$engine[-1L] * sum(fit$data$weights),
               tolerance = 1e-10)
  # when that solve fails, the best point solved (here the start) is returned
  real_setup <- .mxl_blp_setup
  n <- 0L
  s0 <- NULL
  local_mocked_bindings(.mxl_blp_setup = function(...) {
    s <- real_setup(...)
    solve <- s$solve
    s$solve <- function(theta_c, delta0, max_iter, best = NULL) {
      n <<- n + 1L
      r <- solve(theta_c, delta0, max_iter, best)
      if (n == 1L) s0 <<- c(r, list(theta_c = theta_c))
      if (n == 2L) {
        r$converged <- FALSE
        r$error <- "the kernel ran out of memory"
      }
      r
    }
    s
  })
  w <- expect_warning(f2 <- blpc_fit(dt, concentrate_out_asc = "blp", optimizer = opt_away),
                      class = "choicer_blp_unconverged")
  expect_match(conditionMessage(w), paste0(
    "the share equations could not be solved at the optimizer's solution (the ",
    "kernel ran out of memory), so the best point solved is returned"), fixed = TRUE)
  expect_identical(n, 2L)
  expect_identical(unname(coef(f2)[f2$param_map$asc]), s0$delta)
  expect_identical(unname(coef(f2)[-f2$param_map$asc]), unname(s0$theta_c))
  expect_identical(f2$loglik, -s0$extra$objective)
  expect_identical(f2$loglik, ll(f2))
  expect_identical(f2$optimizer$concentration$failed, 0L)
  expect_match(f2$message, "share-matching constants: the share equations could not be solved",
               fixed = TRUE)
})

test_that("market-level targets: names, clusters and decision makers within markets", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(611, N = 300L, markets = 2L)
  s_m <- c(`1` = 0.3, `2` = 0.4, `3` = 0.35, `4` = 0.3)
  # their covariance, from 500 and 800 choices in the two markets (named
  # by label, as diag() and %o% name it)
  V_m <- (diag(s_m) - s_m %o% s_m) / c(500, 500, 800, 800)
  V_m[1:2, 3:4] <- V_m[3:4, 1:2] <- 0
  V_m <- (V_m + t(V_m)) / 2
  f1 <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = s_m,
                 target_vcov = V_m)
  # target_vcov in target_shares' order, named or not, or by its names
  f2 <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = rev(s_m),
                 target_vcov = unname(V_m[4:1, 4:1]))
  f3 <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = rev(s_m),
                 target_vcov = V_m)
  expect_identical(coef(f2), coef(f1))
  expect_identical(f2$vcov, f1$vcov)
  expect_identical(f3$vcov, f1$vcov)
  # (which the order matters to)
  f4 <- blpc_fit(dt, concentrate_out_asc = "blp", target_shares = s_m,
                 target_vcov = unname(V_m[4:1, 4:1]))
  expect_gt(max(abs(f4$vcov - f1$vcov)), 1e-3 * max(abs(f1$vcov)))
  # the fit matches each market's targets within it
  mk <- f1$share_targets$market
  W_m <- f1$share_targets$market_weight[mk]
  expect_equal(blpc_counts(f1) / W_m, unname(s_m), tolerance = 1e-10)
  # the outside option's share first is for one market only
  expect_error(blpc_fit(dt, concentrate_out_asc = "blp",
                        target_shares = c(0.4, 0.15, 0.15, 0.15, 0.15)),
               "taken for one market only; these data have 2", fixed = TRUE)
  # clusters across markets, at the fit and after it
  dt[, cl := id %% 5L]
  expect_error(blpc_fit(dt, concentrate_out_asc = "blp", target_shares = s_m,
                        cluster_col = "cl", se_method = "cluster"),
               "needs clusters within markets", fixed = TRUE)
  cl <- dt[, .(cl = cl[1L]), by = id]
  expect_error(vcov(f1, type = "cluster", cluster = stats::setNames(cl$cl, cl$id)),
               "need clusters within markets", fixed = TRUE)
  # a decision maker in two markets
  dp <- blpc_data(612, N = 100L, T = 2L, markets = 2L)
  dp[id %% 2L == 0L, pid := pid + 1L]   # moves situations between persons
  dp <- dp[order(pid, id)]
  expect_error(blpc_fit(dp, concentrate_out_asc = "blp", target_shares = s_m,
                        panel = TRUE),
               "each decision maker's choice situations in one market", fixed = TRUE)
})

test_that("a choice-based sample's population shares are the default targets", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(613, N = 600L)
  Q <- c(`0` = 0.3, `1` = 0.1, `2` = 0.15, `3` = 0.2, `4` = 0.25)
  dt_w <- wesml_weights(dt, "id", "alt", "choice", Q = Q, attach = TRUE,
                        outside_opt_label = 0L, include_outside_option = TRUE)
  expect_message(
    fit <- withCallingHandlers(run_mxlogit(
      data = dt_w, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "w1"), random_var_cols = "w1", S = 20L,
      include_outside_option = TRUE, outside_opt_label = 0L,
      se_method = "sandwich", concentrate_out_asc = "blp"),
      message = function(m) {
        if (!grepl("population shares (its Q)", conditionMessage(m), fixed = TRUE)) {
          invokeRestart("muffleMessage")
        }
      }),
    "population shares (its Q), taken as exact", fixed = TRUE)
  expect_identical(fit$share_targets$source, "wesml")
  expect_equal(blpc_counts(fit) / sum(fit$data$weights), unname(Q[-1L]),
               tolerance = 1e-10)
  expect_false(is.null(fit$vcov))
  # Q failing the share inversion's checks: the message on Q comes first, and
  # the error names Q (the check made to fail: the weights make the weighted
  # sample's shares Q, which then pass)
  real_check <- blp_check_targets
  local_mocked_bindings(blp_check_targets = function(..., target_shares = NULL) {
    if (is.null(target_shares)) return(real_check(...))
    stop("Alternative 2's target share (0.15) is not below the weight of the ",
         "choice situations that offer it (0.1): no constants attain it.")
  })
  seen <- character()
  expect_error(withCallingHandlers(
    run_mxlogit(
      data = dt_w, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "w1"), random_var_cols = "w1", S = 20L,
      include_outside_option = TRUE, outside_opt_label = 0L,
      se_method = "sandwich", concentrate_out_asc = "blp"),
    message = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }),
    paste0("the choice-based sample's population shares Q, as targets, fail the ",
           "share inversion's checks. Alternative 2's target share (0.15) is not ",
           "below the weight of the choice situations that offer it (0.1): no ",
           "constants attain it. (In its numbering, 2 is \"2\".)"), fixed = TRUE)
  expect_true(any(grepl("population shares (its Q), taken as exact", seen, fixed = TRUE)))
  # target_vcov describes supplied targets only
  expect_error(blpc_fit(dt, concentrate_out_asc = "blp",
                        target_vcov = diag(4)),
               "sampling covariance of supplied `target_shares`", fixed = TRUE)
})

test_that("draws formed on the fly, a start of the user's and scaled coordinates", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- blpc_data(614, N = 120L, T = 3L)
  # generate mode: the share inversion forms each decision maker's Halton
  # block as the kernels do
  fg <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE, draws = "generate",
                 seed = 9L)
  expect_identical(fg$share_targets$source, "sample")
  d <- fg$data
  first <- c(0, cumsum(d$M))[seq_along(d$M)]
  ch <- d$choice_idx > 0L
  expect_equal(blpc_counts(fg), tabulate(d$alt_idx[first[ch] + d$choice_idx[ch]], 4L),
               tolerance = 1e-10)
  expect_identical(unname(vcov(fg, type = "hessian")), unname(fg$vcov))
  # a supplied start (its constants start the first solve) and column
  # scales reach the same estimates
  f0 <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE)
  f1 <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE,
                 theta_init = unname(coef(f0)) + 0.05)
  f2 <- blpc_fit(dt, concentrate_out_asc = "blp", panel = TRUE, scale_vars = "sd")
  expect_lt(max(abs(coef(f1) - coef(f0)) / f0$se), 1e-3)
  expect_lt(max(abs(coef(f2) - coef(f0)) / f0$se), 1e-3)
})
