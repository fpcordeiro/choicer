# _validation/mxl_mc_helpers.R
#
# Helpers for the Mixed Logit Monte Carlo validation suite.
# Consumed by _validation/mxl_monte_carlo.R.
#
# Scope:
#   - scenario_spec(): scenario DSL returning sim_fun / fit_fun / meta
#   - fit_with_both_ses(): run_mxlogit() + BHHH SE augmentation
#   - restricted_fit(): single-parameter equality restriction via nloptr bounds
#   - hessian_agreement(): analytical vs numerical Hessian at theta_hat
#   - natural_scale_recovery(): delta-method-transformed recovery table
#   - plotting helpers (base R graphics only; no ggplot)
#   - build_report(): composes REPORT.md

suppressPackageStartupMessages({
  library(choicer)
  library(data.table)
})

# ---- Scenario DSL ----------------------------------------------------------

# Every scenario is a list:
#   id       character scalar identifier (e.g. "A", "B", "D")
#   purpose  one-line description
#   sim_fun  function(seed, arg = ...) -> choicer_sim
#   fit_fun  function(sim, ...) -> choicer_mxl fit augmented with bhhh_se, etc.
#   meta     list with N, S, K_w, rc_dist, rc_correlation, rc_mean, R, etc.
#   N_grid   optional vector of N values (Scenario A only)
#   S_grid   optional vector of S values (Scenario D only)

# Halton-normal draw count convention used across the study.
.halton_S <- function(N) {
  S <- 50 * ceiling((1.5 * sqrt(N)) / 50)
  max(S, 50)
}

# Sigma fixtures from the plan's "DGP fixtures" section.
.sigma_uncorrelated  <- diag(c(1.0, 1.5))
.sigma_correlated_2  <- matrix(c(1.0, 0.5, 0.5, 1.5), nrow = 2)
.sigma_correlated_3  <- matrix(
  c(1.0, 0.4, 0.2,
    0.4, 1.5, 0.3,
    0.2, 0.3, 0.8),
  nrow = 3, byrow = TRUE
)
.sigma_weak          <- diag(c(0.1, 1.5))

# Factory. `quick` halves the default R per scenario and reduces Scenario A's
# N-grid to {1k, 5k}.
scenario_spec <- function(id, quick = FALSE) {
  beta_true <- c(0.8, -0.6)
  J <- 8L

  mk_mxl_fit_fun <- function(K_w, rc_dist, rc_correlation, rc_mean,
                             random_var_cols, S_override = NULL) {
    force(K_w); force(rc_dist); force(rc_correlation); force(rc_mean)
    force(random_var_cols); force(S_override)
    function(sim, S = NULL) {
      if (is.null(S)) {
        S <- if (!is.null(S_override)) S_override else .halton_S(sim$settings$N)
      }
      # With outside_option = TRUE in the DGP, alt = 0 is present in the
      # data. We pass outside_opt_label = 0L, include_outside_option = FALSE
      # so prepare_mxl_data() re-levels alt_int with the outside first; the
      # C++ fixes alt_int = 1 (alt = 0) ASC to zero. The remaining J ASCs
      # match truth$delta (length J) exactly: no drop-first-ASC fixup is
      # triggered in recovery_table(). This matches the convention used by
      # inst/simulations/mxl_simulation.R and lines up with the plan's
      # "all J inside ASCs are free, truth has J ASCs" spirit -- the
      # outside's delta is identified only up to the zero reference, so
      # truth$delta has J entries and the fit has J free ASCs.
      input <- prepare_mxl_data(
        data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
        covariate_cols = c("x1", "x2"),
        random_var_cols = random_var_cols,
        outside_opt_label = 0L,
        include_outside_option = FALSE,
        rc_correlation = rc_correlation
      )
      eta_draws <- get_halton_normals(S, input$N, K_w)
      fit <- run_mxlogit(
        input_data = input, eta_draws = eta_draws,
        rc_dist = rc_dist, rc_correlation = rc_correlation,
        rc_mean = rc_mean, use_asc = TRUE,
        se_method = "hessian",
        control = list(print_level = 0L, maxeval = 400L)
      )
      # Augment with BHHH standard errors at theta_hat.
      fit <- .augment_bhhh(fit, input, eta_draws)
      # Stash hyperparameters for downstream helpers.
      fit$._scenario <- list(
        id = id, input = input, eta_draws = eta_draws,
        K_w = K_w, rc_dist = rc_dist, rc_correlation = rc_correlation,
        rc_mean = rc_mean, random_var_cols = random_var_cols, S = S
      )
      fit
    }
  }

  mk_sim_fun <- function(N, Sigma, mu = NULL, rc_dist) {
    force(N); force(Sigma); force(mu); force(rc_dist)
    function(seed) {
      K_w <- ncol(Sigma)
      simulate_mxl_data(
        N = N, J = J,
        beta = beta_true,
        Sigma = Sigma,
        mu = mu, rc_dist = rc_dist,
        seed = seed,
        outside_option = TRUE,
        vary_choice_set = TRUE
      )
    }
  }

  switch(id,
    A = {
      N_grid <- if (quick) c(1000L, 5000L) else c(1000L, 2500L, 5000L, 10000L)
      R_map  <- if (quick) setNames(rep(50L, length(N_grid)), N_grid)
                else setNames(c(1000L, 1000L, 1000L, 2000L), c(1000, 2500, 5000, 10000))
      list(
        id = "A", purpose = "Consistency across N",
        N_grid = N_grid,
        R_map  = R_map,
        sim_fun_for_N = function(N) mk_sim_fun(N, .sigma_uncorrelated,
                                              rc_dist = c(0L, 0L)),
        fit_fun = mk_mxl_fit_fun(
          K_w = 2L, rc_dist = c(0L, 0L),
          rc_correlation = FALSE, rc_mean = FALSE,
          random_var_cols = c("w1", "w2")
        ),
        meta = list(K_w = 2L, rc_dist = c(0L, 0L), rc_correlation = FALSE,
                    rc_mean = FALSE, purpose = "Consistency across N",
                    Sigma = .sigma_uncorrelated, beta = beta_true,
                    lr_restriction = list(idx = 2L, value = beta_true[2]))
        )
    },
    B = {
      N <- 5000L
      R <- if (quick) 50L else 1000L
      list(
        id = "B", purpose = "Correlated Sigma, K_w = 2",
        N = N, R = R,
        sim_fun = mk_sim_fun(N, .sigma_correlated_2, rc_dist = c(0L, 0L)),
        fit_fun = mk_mxl_fit_fun(
          K_w = 2L, rc_dist = c(0L, 0L),
          rc_correlation = TRUE, rc_mean = FALSE,
          random_var_cols = c("w1", "w2")
        ),
        meta = list(K_w = 2L, rc_dist = c(0L, 0L), rc_correlation = TRUE,
                    rc_mean = FALSE, Sigma = .sigma_correlated_2,
                    beta = beta_true)
      )
    },
    B2 = {
      N <- 5000L
      R <- if (quick) 50L else 1000L
      list(
        id = "B2", purpose = "Correlated Sigma, K_w = 3",
        N = N, R = R,
        sim_fun = mk_sim_fun(N, .sigma_correlated_3, rc_dist = rep(0L, 3)),
        fit_fun = mk_mxl_fit_fun(
          K_w = 3L, rc_dist = rep(0L, 3),
          rc_correlation = TRUE, rc_mean = FALSE,
          random_var_cols = c("w1", "w2", "w3")
        ),
        meta = list(K_w = 3L, rc_dist = rep(0L, 3), rc_correlation = TRUE,
                    rc_mean = FALSE, Sigma = .sigma_correlated_3,
                    beta = beta_true)
      )
    },
    C = {
      N <- 5000L
      R <- if (quick) 50L else 1000L
      mu_true <- c(-0.2, 0.4)
      list(
        id = "C", purpose = "Log-normal RC with mu recovery",
        N = N, R = R,
        sim_fun = mk_sim_fun(N, .sigma_uncorrelated,
                            mu = mu_true, rc_dist = c(1L, 0L)),
        fit_fun = mk_mxl_fit_fun(
          K_w = 2L, rc_dist = c(1L, 0L),
          rc_correlation = FALSE, rc_mean = TRUE,
          random_var_cols = c("w1", "w2")
        ),
        meta = list(K_w = 2L, rc_dist = c(1L, 0L), rc_correlation = FALSE,
                    rc_mean = TRUE, Sigma = .sigma_uncorrelated,
                    mu = mu_true, beta = beta_true)
      )
    },
    D = {
      N <- 5000L
      R <- if (quick) 50L else 300L
      S_grid <- if (quick) c(25L, 100L) else c(25L, 50L, 100L, 250L, 500L)
      list(
        id = "D", purpose = "Simulation bias O(1 / S)",
        N = N, R = R, S_grid = S_grid,
        sim_fun = mk_sim_fun(N, .sigma_uncorrelated, rc_dist = c(0L, 0L)),
        fit_fun_for_S = function(S_override) mk_mxl_fit_fun(
          K_w = 2L, rc_dist = c(0L, 0L),
          rc_correlation = FALSE, rc_mean = FALSE,
          random_var_cols = c("w1", "w2"),
          S_override = S_override
        ),
        meta = list(K_w = 2L, rc_dist = c(0L, 0L), rc_correlation = FALSE,
                    rc_mean = FALSE, Sigma = .sigma_uncorrelated,
                    beta = beta_true)
      )
    },
    F = {
      N <- 5000L
      R <- if (quick) 50L else 1000L
      list(
        id = "F", purpose = "Weak identification (small sigma)",
        N = N, R = R,
        sim_fun = mk_sim_fun(N, .sigma_weak, rc_dist = c(0L, 0L)),
        fit_fun = mk_mxl_fit_fun(
          K_w = 2L, rc_dist = c(0L, 0L),
          rc_correlation = FALSE, rc_mean = FALSE,
          random_var_cols = c("w1", "w2")
        ),
        meta = list(K_w = 2L, rc_dist = c(0L, 0L), rc_correlation = FALSE,
                    rc_mean = FALSE, Sigma = .sigma_weak,
                    beta = beta_true)
      )
    },
    stop("Unknown scenario id: '", id, "'")
  )
}

# ---- BHHH augmentation -----------------------------------------------------

.augment_bhhh <- function(fit, input, eta_draws) {
  bhhh <- tryCatch(
    mxl_bhhh_parallel(
      theta = fit$coefficients,
      X = input$X, W = input$W,
      alt_idx = input$alt_idx, choice_idx = input$choice_idx,
      M = input$M, weights = input$weights,
      eta_draws = eta_draws,
      rc_dist = fit$rc_dist,
      rc_correlation = fit$rc_correlation,
      rc_mean = fit$rc_mean,
      use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option
    ),
    error = function(e) NULL
  )
  if (is.null(bhhh)) {
    fit$bhhh_se <- rep(NA_real_, length(fit$coefficients))
    fit$bhhh_vcov <- NULL
    return(fit)
  }
  vv <- tryCatch(solve(bhhh), error = function(e) NULL)
  if (is.null(vv)) {
    fit$bhhh_se <- rep(NA_real_, length(fit$coefficients))
    fit$bhhh_vcov <- NULL
    return(fit)
  }
  se_bhhh <- sqrt(pmax(diag(vv), 0))
  names(se_bhhh) <- names(fit$coefficients)
  fit$bhhh_vcov <- vv
  fit$bhhh_se   <- se_bhhh
  fit
}

# ---- Restricted fit (LR test, Scenario A) ----------------------------------

# Impose a single equality restriction beta[idx] = value via nloptr bounds
# (lower = upper = value). Returns fitted value + restricted loglik.
restricted_fit <- function(sim, scenario, restriction) {
  K_w <- scenario$meta$K_w
  rc_dist <- scenario$meta$rc_dist
  rc_correlation <- scenario$meta$rc_correlation
  rc_mean <- scenario$meta$rc_mean

  input <- prepare_mxl_data(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    random_var_cols = paste0("w", seq_len(K_w)),
    outside_opt_label = 0L,
    include_outside_option = FALSE,
    rc_correlation = rc_correlation
  )
  eta_draws <- get_halton_normals(.halton_S(sim$settings$N), input$N, K_w)

  J <- nrow(input$alt_mapping)
  K_x <- ncol(input$X)
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- J - 1
  n_params <- K_x + mu_size + L_size + n_asc

  theta_init <- rep(0, n_params)
  theta_init[restriction$idx] <- restriction$value
  lower <- rep(-Inf, n_params); upper <- rep(Inf, n_params)
  lower[restriction$idx] <- restriction$value
  upper[restriction$idx] <- restriction$value

  eval_f <- function(theta) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input$X, W = input$W,
      alt_idx = input$alt_idx, choice_idx = input$choice_idx,
      M = input$M, weights = input$weights,
      rc_dist = rc_dist, rc_correlation = rc_correlation,
      rc_mean = rc_mean, eta_draws = eta_draws,
      use_asc = TRUE, include_outside_option = input$include_outside_option
    )
  }

  opts <- list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-8,
               maxeval = 400L, print_level = 0L)
  raw <- tryCatch(
    nloptr::nloptr(x0 = theta_init, eval_f = eval_f,
                   lb = lower, ub = upper, opts = opts),
    error = function(e) NULL
  )
  if (is.null(raw)) return(list(loglik_restricted = NA_real_, convergence = NA))
  list(loglik_restricted = -raw$objective, convergence = raw$status)
}

# ---- Hessian analytical-vs-numerical agreement -----------------------------

hessian_agreement <- function(fit, sim = NULL) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    return(list(max_rel_err = NA_real_, note = "numDeriv not installed"))
  }
  s <- fit$._scenario
  if (is.null(s)) {
    return(list(max_rel_err = NA_real_, note = "fit missing _scenario"))
  }
  input <- s$input; eta_draws <- s$eta_draws

  # Scalar negative-loglik for numDeriv. The parallel entry point returns a
  # list with $objective and $gradient so we wrap with a $objective accessor.
  eval_scalar <- function(theta) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input$X, W = input$W,
      alt_idx = input$alt_idx, choice_idx = input$choice_idx,
      M = input$M, weights = input$weights,
      rc_dist = fit$rc_dist,
      rc_correlation = fit$rc_correlation,
      rc_mean = fit$rc_mean,
      eta_draws = eta_draws,
      use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option
    )$objective
  }

  H_num <- tryCatch(
    numDeriv::hessian(eval_scalar, fit$coefficients, method = "Richardson"),
    error = function(e) NULL
  )
  H_ana <- tryCatch(
    mxl_hessian_parallel(
      theta = fit$coefficients,
      X = input$X, W = input$W,
      alt_idx = input$alt_idx, choice_idx = input$choice_idx,
      M = input$M, weights = input$weights,
      eta_draws = eta_draws,
      rc_dist = fit$rc_dist, rc_correlation = fit$rc_correlation,
      rc_mean = fit$rc_mean, use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option
    ),
    error = function(e) NULL
  )
  if (is.null(H_num) || is.null(H_ana)) {
    return(list(max_rel_err = NA_real_, note = "Hessian evaluation failed"))
  }
  denom <- max(abs(H_num))
  if (!is.finite(denom) || denom <= 0) {
    return(list(max_rel_err = NA_real_, note = "H_num ~= 0"))
  }
  err <- max(abs(H_ana - H_num)) / denom
  list(max_rel_err = err, note = NA_character_)
}

# ---- Natural-scale recovery (delta method on Cholesky / log-normal) --------

# Internal delta-method transformer. Duplicates the core of
# R/methods.R:apply_mxl_delta_method() so the validation suite does not depend
# on an unexported internal. Keeps Step 6 of the plan's orchestration (factor
# a .mxl_delta_sigma helper) deferred: the existing summary method is left
# untouched.
.mxl_delta_sigma <- function(theta, vcov_mat, param_map, K_w, rc_correlation,
                             rc_dist, rc_mean) {
  est <- theta
  se_out <- sqrt(pmax(diag(vcov_mat), 0))

  if (rc_mean && !is.null(param_map$mu)) {
    idx_mu <- param_map$mu
    for (k in seq_len(K_w)) {
      if (rc_dist[k] == 1L) {
        ci <- idx_mu[k]
        mu_hat <- est[ci]
        est[ci] <- exp(mu_hat)
        se_out[ci] <- exp(mu_hat) * se_out[ci]
      }
    }
  }
  if (!is.null(param_map$sigma)) {
    idx_s <- param_map$sigma
    L_params_hat <- theta[idx_s]
    Sigma_hat <- build_var_mat(L_params_hat, K_w, rc_correlation)
    if (rc_correlation) {
      est[idx_s] <- Sigma_hat[lower.tri(Sigma_hat, diag = TRUE)]
    } else {
      est[idx_s] <- diag(Sigma_hat)
    }
    J_mat <- jacobian_vech_Sigma(L_params_hat, K_w, rc_correlation)
    V_sub <- vcov_mat[idx_s, idx_s, drop = FALSE]
    V_sigma <- J_mat %*% V_sub %*% t(J_mat)
    se_out[idx_s] <- sqrt(pmax(diag(V_sigma), 0))
  }
  list(estimate = est, se = se_out)
}

# Build a recovery_table-shaped data.table on the natural scale. Raw
# parameters become exp(mu) and vech(Sigma) (or diag(Sigma) when uncorrelated).
natural_scale_recovery <- function(fit, sim, level = 0.95) {
  pm <- fit$param_map
  vv <- fit$vcov
  if (is.null(vv)) return(NULL)

  transformed <- .mxl_delta_sigma(
    theta = fit$coefficients,
    vcov_mat = vv,
    param_map = pm,
    K_w = fit$._scenario$K_w,
    rc_correlation = fit$rc_correlation,
    rc_dist = fit$rc_dist,
    rc_mean = fit$rc_mean
  )

  est_nat <- transformed$estimate
  se_nat  <- transformed$se

  truth <- sim$true_params
  z_crit <- stats::qnorm(1 - (1 - level) / 2)

  # Natural-scale truth block:
  #   beta: as is
  #   mu  : exp(mu) elementwise where rc_dist == 1
  #   sigma: vech(Sigma) if correlated, else diag(Sigma)
  #   asc : as is
  rows <- list()

  idx_beta <- pm$beta
  rows[[length(rows) + 1]] <- data.table(
    parameter = names(fit$coefficients)[idx_beta], group = "beta",
    true = as.numeric(truth$beta),
    estimate = as.numeric(est_nat[idx_beta]),
    se = as.numeric(se_nat[idx_beta])
  )

  if (!is.null(pm$mu) && !is.null(truth$mu)) {
    idx_mu <- pm$mu
    mu_true_nat <- truth$mu
    for (k in seq_along(idx_mu)) {
      if (fit$rc_dist[k] == 1L) mu_true_nat[k] <- exp(truth$mu[k])
    }
    rows[[length(rows) + 1]] <- data.table(
      parameter = names(fit$coefficients)[idx_mu], group = "mu",
      true = as.numeric(mu_true_nat),
      estimate = as.numeric(est_nat[idx_mu]),
      se = as.numeric(se_nat[idx_mu])
    )
  }

  if (!is.null(pm$sigma)) {
    idx_s <- pm$sigma
    Sigma_true <- truth$Sigma
    if (fit$rc_correlation) {
      sigma_true_nat <- Sigma_true[lower.tri(Sigma_true, diag = TRUE)]
    } else {
      sigma_true_nat <- diag(Sigma_true)
    }
    rows[[length(rows) + 1]] <- data.table(
      parameter = names(fit$coefficients)[idx_s], group = "sigma",
      true = as.numeric(sigma_true_nat),
      estimate = as.numeric(est_nat[idx_s]),
      se = as.numeric(se_nat[idx_s])
    )
  }

  if (!is.null(pm$asc)) {
    idx_a <- pm$asc
    asc_true <- truth$delta[seq_along(idx_a)]  # OO case: length matches
    rows[[length(rows) + 1]] <- data.table(
      parameter = names(fit$coefficients)[idx_a], group = "asc",
      true = as.numeric(asc_true),
      estimate = as.numeric(est_nat[idx_a]),
      se = as.numeric(se_nat[idx_a])
    )
  }

  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, `:=`(
    bias         = estimate - true,
    rel_bias_pct = 100 * (estimate - true) / true,
    z_vs_true    = (estimate - true) / se,
    lower_ci     = estimate - z_crit * se,
    upper_ci     = estimate + z_crit * se,
    covers       = (true >= estimate - z_crit * se) &
                   (true <= estimate + z_crit * se)
  )]
  data.table::setcolorder(out, c("parameter", "group", "true", "estimate",
                                 "se", "bias", "rel_bias_pct", "z_vs_true",
                                 "lower_ci", "upper_ci", "covers"))
  out
}

# ---- Augmented per-rep record ---------------------------------------------

# Glue: combine raw-scale recovery rows with natural-scale recovery rows and
# attach BHHH SEs on the raw scale. Returns a long data.table tagged with
# scale in $scale.
augment_replication <- function(fit, sim) {
  raw <- recovery_table(fit, sim)
  raw[, scale := "raw"]
  if (!is.null(fit$bhhh_se)) {
    pm <- fit$param_map
    raw[, se_bhhh := NA_real_]
    # For each row, pull the bhhh se by parameter name
    se_by_name <- fit$bhhh_se
    raw[, se_bhhh := se_by_name[parameter]]
  }
  nat <- tryCatch(natural_scale_recovery(fit, sim),
                  error = function(e) NULL)
  if (!is.null(nat)) {
    nat[, scale := "natural"]
    nat[, se_bhhh := NA_real_]  # BHHH transform not tracked on natural scale
    out <- rbindlist(list(raw, nat), use.names = TRUE, fill = TRUE)
  } else {
    out <- raw
  }
  out
}

# ---- Plotting helpers (base R graphics) ------------------------------------

.png_device <- function(path, width = 7, height = 6, res = 150) {
  grDevices::png(path, width = width, height = height, units = "in", res = res)
}

plot_qq <- function(reps, scenario_id, outfile) {
  # Faceted QQ of z_vs_true vs N(0, 1), one panel per parameter.
  params <- unique(reps$parameter)
  n_params <- length(params)
  ncol <- min(4L, max(1L, ceiling(sqrt(n_params))))
  nrow <- ceiling(n_params / ncol)

  .png_device(outfile, width = 2 * ncol + 1.5, height = 2 * nrow + 1)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
  on.exit(par(op), add = TRUE, after = FALSE)
  for (p in params) {
    z <- reps[parameter == p, z_vs_true]
    z <- z[is.finite(z)]
    if (!length(z)) { plot.new(); title(p); next }
    stats::qqnorm(z, main = p, cex.main = 0.85, pch = 16, cex = 0.4)
    stats::qqline(z, col = "red", lwd = 1)
  }
  invisible(NULL)
}

plot_consistency <- function(consistency_tbl, outfile) {
  # Expects columns: N, parameter, group, abs_bias, rmse.
  .png_device(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE, after = FALSE)

  params <- unique(consistency_tbl$parameter)
  cols <- grDevices::rainbow(length(params))

  # log|bias| vs log N with fitted slope for each parameter.
  x_range <- log(range(consistency_tbl$N))
  y_range <- log(range(consistency_tbl$abs_bias[consistency_tbl$abs_bias > 0], na.rm = TRUE))
  plot(NA, xlim = x_range, ylim = y_range,
       xlab = "log N", ylab = "log |bias|",
       main = "Consistency: |bias| vs N")
  for (i in seq_along(params)) {
    d <- consistency_tbl[parameter == params[i] & abs_bias > 0]
    points(log(d$N), log(d$abs_bias), col = cols[i], pch = 16)
    if (nrow(d) >= 2) {
      lines(log(d$N), log(d$abs_bias), col = cols[i])
    }
  }
  abline(lm(y ~ x, data = data.frame(
    x = log(consistency_tbl$N),
    y = log(pmax(consistency_tbl$abs_bias, 1e-12))
  )), col = "grey40", lty = 2)

  y_range2 <- log(range(consistency_tbl$rmse[consistency_tbl$rmse > 0], na.rm = TRUE))
  plot(NA, xlim = x_range, ylim = y_range2,
       xlab = "log N", ylab = "log RMSE", main = "Consistency: RMSE vs N")
  for (i in seq_along(params)) {
    d <- consistency_tbl[parameter == params[i] & rmse > 0]
    points(log(d$N), log(d$rmse), col = cols[i], pch = 16)
    if (nrow(d) >= 2) lines(log(d$N), log(d$rmse), col = cols[i])
  }
  abline(lm(y ~ x, data = data.frame(
    x = log(consistency_tbl$N),
    y = log(pmax(consistency_tbl$rmse, 1e-12))
  )), col = "grey40", lty = 2)
  invisible(NULL)
}

plot_coverage <- function(coverage_tbl, outfile) {
  # Expects columns: N, parameter, cov95, cov95_lower, cov95_upper.
  .png_device(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE, after = FALSE)

  params <- unique(coverage_tbl$parameter)
  cols <- grDevices::rainbow(length(params))
  x_range <- range(coverage_tbl$N)
  plot(NA, xlim = x_range, ylim = c(0.8, 1),
       xlab = "N", ylab = "empirical coverage", log = "x",
       main = "95 pct Wald coverage with Wilson band")
  abline(h = 0.95, col = "grey40", lty = 2)
  for (i in seq_along(params)) {
    d <- coverage_tbl[parameter == params[i]]
    points(d$N, d$cov95, col = cols[i], pch = 16)
    segments(d$N, d$cov95_lower, d$N, d$cov95_upper, col = cols[i])
  }
  legend("bottomright", legend = params, col = cols, pch = 16,
         cex = 0.7, bty = "n")
  invisible(NULL)
}

plot_sim_bias <- function(simbias_tbl, outfile) {
  # Expects columns: S, parameter, abs_bias.
  .png_device(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE, after = FALSE)

  params <- unique(simbias_tbl$parameter)
  cols <- grDevices::rainbow(length(params))
  x <- 1 / simbias_tbl$S
  plot(NA, xlim = range(x), ylim = range(simbias_tbl$abs_bias, na.rm = TRUE),
       xlab = "1 / S", ylab = "|bias|",
       main = "Simulation bias vs 1 / S")
  for (i in seq_along(params)) {
    d <- simbias_tbl[parameter == params[i]]
    points(1 / d$S, d$abs_bias, col = cols[i], pch = 16)
    if (nrow(d) >= 2) {
      fit <- lm(abs_bias ~ I(1 / S), data = d)
      abline(fit, col = cols[i], lty = 2)
    }
  }
  legend("topleft", legend = params, col = cols, pch = 16,
         cex = 0.7, bty = "n")
  invisible(NULL)
}

plot_se_ratios <- function(se_tbl, outfile) {
  # Expects: scenario, parameter, hess_ratio, bhhh_ratio.
  .png_device(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mar = c(5, 4, 2, 1))
  on.exit(par(op), add = TRUE, after = FALSE)

  # Stack by scenario x parameter for a grouped bar-ish scatter.
  se_tbl <- copy(se_tbl)
  se_tbl[, label := paste0(scenario, ":", parameter)]
  n <- nrow(se_tbl)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, 1.5),
       xlab = "", ylab = "SE / sd_emp",
       main = "Information-matrix equality check",
       xaxt = "n")
  abline(h = c(0.9, 1.0, 1.1), col = c("grey50", "black", "grey50"),
         lty = c(2, 1, 2))
  points(seq_len(n) - 0.1, se_tbl$hess_ratio, col = "blue", pch = 16)
  points(seq_len(n) + 0.1, se_tbl$bhhh_ratio, col = "red", pch = 17)
  axis(1, at = seq_len(n), labels = se_tbl$label, las = 2, cex.axis = 0.6)
  legend("topright", legend = c("Hessian", "BHHH"),
         col = c("blue", "red"), pch = c(16, 17), bty = "n", cex = 0.8)
  invisible(NULL)
}

plot_lr_cdf <- function(lr_values, outfile, df = 1L) {
  lr_values <- lr_values[is.finite(lr_values) & lr_values >= 0]
  if (!length(lr_values)) return(invisible(NULL))
  .png_device(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE, after = FALSE)

  x_grid <- seq(0, max(quantile(lr_values, 0.99), qchisq(0.999, df)),
                length.out = 500)
  plot(ecdf(lr_values), main = sprintf("LR empirical CDF vs chi2_%d", df),
       xlab = "LR", xlim = range(x_grid), cex = 0.3,
       pch = 16, col = "black")
  lines(x_grid, pchisq(x_grid, df = df), col = "red", lwd = 2)
  legend("bottomright", legend = c("empirical", sprintf("chi2_%d", df)),
         col = c("black", "red"), lty = 1, bty = "n", cex = 0.8)
  invisible(NULL)
}

# ---- Report builder --------------------------------------------------------

build_report <- function(scenario_results, outfile) {
  # scenario_results: named list keyed by scenario id. Each element contains
  #   $asymptotics_raw (choicer_mc_asymptotics)
  #   $asymptotics_natural (choicer_mc_asymptotics or NULL)
  #   $meta (from scenario_spec)
  #   $extras  list of narrative bits (lr_ks_p, lr_size, hessian_err, ...)
  lines <- character(0)
  push <- function(...) lines <<- c(lines, paste0(...))

  push("# MXL Monte Carlo Validation Report")
  push("")
  push("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  push("")
  if (isTRUE(getOption("choicer.quick_mode", FALSE))) {
    push("**QUICK mode**: reduced R per scenario and Scenario A grid trimmed. ",
         "This is a pre-flight sanity run, not the publication-quality ",
         "certification pass. See `_validation/README.md` for the full-run ",
         "recipe.")
    push("")
  }

  push("## Pass / fail matrix")
  push("")
  header <- paste0("| scenario | parameter | group | ",
                   "pass_bias | pass_se_ratio | pass_cov95 | pass_skew | pass_kurt |")
  sep <- paste0("|----------|-----------|-------|",
                "-----------|----------------|------------|-----------|-----------|")
  push(header); push(sep)

  for (nm in names(scenario_results)) {
    sr <- scenario_results[[nm]]
    a <- sr$asymptotics_raw
    if (is.null(a) || !nrow(a)) next
    for (i in seq_len(nrow(a))) {
      push(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                   nm, a$parameter[i], a$group[i],
                   a$pass_bias[i], a$pass_se_ratio[i], a$pass_cov95[i],
                   a$pass_skew[i], a$pass_kurt[i]))
    }
  }
  push("")

  push("## Per-scenario narratives")
  push("")
  for (nm in names(scenario_results)) {
    sr <- scenario_results[[nm]]
    push("### Scenario ", nm, ": ", sr$meta$purpose %||% "")
    push("")

    a <- sr$asymptotics_raw
    if (!is.null(a) && nrow(a)) {
      pass_cols <- c("pass_bias", "pass_se_ratio", "pass_cov95",
                     "pass_skew", "pass_kurt")
      all_pass <- rep(TRUE, nrow(a))
      for (col in pass_cols) all_pass <- all_pass & as.logical(a[[col]])
      n_all <- sum(all_pass, na.rm = TRUE)
      push(sprintf("- raw scale: %d / %d parameters pass all 5 flags",
                   n_all, nrow(a)))
      bias_range <- range(abs(a$bias_mc_se), na.rm = TRUE)
      se_range <- range(a$se_ratio, na.rm = TRUE)
      push(sprintf("- |bias / MC-SE| range: [%.2f, %.2f]",
                   bias_range[1], bias_range[2]))
      push(sprintf("- se_ratio (Hessian) range: [%.3f, %.3f]",
                   se_range[1], se_range[2]))
      push(sprintf("- cov95 range: [%.3f, %.3f]",
                   min(a$cov95, na.rm = TRUE),
                   max(a$cov95, na.rm = TRUE)))
    }

    ab <- sr$asymptotics_raw_bhhh
    if (!is.null(ab) && nrow(ab)) {
      se_range_b <- range(ab$se_ratio, na.rm = TRUE)
      cov95_range_b <- range(ab$cov95, na.rm = TRUE)
      push(sprintf("- se_ratio (BHHH) range: [%.3f, %.3f]",
                   se_range_b[1], se_range_b[2]))
      push(sprintf("- cov95 (BHHH) range: [%.3f, %.3f]",
                   cov95_range_b[1], cov95_range_b[2]))
    }

    nat <- sr$asymptotics_natural
    if (!is.null(nat) && nrow(nat)) {
      pass_cols <- c("pass_bias", "pass_se_ratio", "pass_cov95",
                     "pass_skew", "pass_kurt")
      all_pass <- rep(TRUE, nrow(nat))
      for (col in pass_cols) all_pass <- all_pass & as.logical(nat[[col]])
      n_all <- sum(all_pass, na.rm = TRUE)
      push(sprintf("- natural scale: %d / %d parameters pass all 5 flags",
                   n_all, nrow(nat)))
    }

    ex <- sr$extras
    if (!is.null(ex)) {
      if (!is.null(ex$hessian_err)) {
        push(sprintf("- analytical-vs-numerical Hessian max rel err: %.2e",
                     ex$hessian_err))
      }
      if (!is.null(ex$lr_ks_p)) {
        push(sprintf("- LR KS p-value vs chi2_1: %.4f", ex$lr_ks_p))
      }
      if (!is.null(ex$lr_size)) {
        push(sprintf("- LR nominal-5pct size: %.3f (Wilson band %.3f -- %.3f)",
                     ex$lr_size, ex$lr_size_lower, ex$lr_size_upper))
      }
      if (!is.null(ex$conv_rate)) {
        push(sprintf("- convergence rate: %.3f", ex$conv_rate))
      }
      if (!is.null(ex$simbias_slopes) && nrow(ex$simbias_slopes)) {
        push("- claim 5: OLS slope of |bias| on 1/S (positive slope + t > 0 ",
             "=> bias decays toward zero as S grows):")
        sb <- ex$simbias_slopes
        for (i in seq_len(nrow(sb))) {
          push(sprintf("  - %s: slope=%.4f (se=%.4f, t=%.2f, p=%.3f)",
                       sb$parameter[i],
                       sb$slope[i], sb$slope_se[i],
                       sb$slope_t[i], sb$slope_p[i]))
        }
      }
    }
    push("")
  }

  writeLines(lines, outfile)
  invisible(outfile)
}
