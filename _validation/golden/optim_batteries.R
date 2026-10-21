# ============================================================================
# 6D follow-up 2 batteries (run_optim() evaluates the model once per point),
# sourced by refactor_identity_check.R after start_batteries.R. Main's optim
# adapter runs the kernel for fn and again for gr; the branch keeps the last
# point's evaluation (.one_point_cache()). Every entry must be identical
# (class I): optim sees the same numbers.
#   * optim_ro : run_optimizer("optim") (and .optimize_in_coordinates() with a
#                "bhhh" map) on logged MNL / MXL / NL kernel objectives: BFGS,
#                L-BFGS-B (free and binding bounds; also forced through
#                control$method), CG, Nelder-Mead, parscale/fnscale, the MXL
#                overflow sentinel on a path that steps into the overflow,
#                generated draws on a panel. Recorded: run_optimizer()'s
#                result, optim's counts included, and the points the kernel
#                evaluated with consecutive repeats dropped (main's fn/gr
#                pairs collapse to the branch's single evaluation).
#   * optim_fit: run_mnlogit / run_mxlogit / run_nestlogit with optimizer =
#                "optim" over scale_vars, se_method, weights, clusters, the
#                outside option, bounds, panels, generated draws, log-normal
#                coefficients, the sentinel's path and control variants:
#                the whole fit (timings stripped) with its conditions.
# Kernel calls per entry differ by design: they are printed and saved next to
# this script (optim_calls<ref>_<lib>.rds), not compared.
# All calls run at one thread (the harness sets it).
# ============================================================================

OC_CALLS <- list()
oc_counter <- new.env()
oc_counter$n <- 0L

oc_collapse <- function(pts) {
  if (length(pts) < 2L) return(pts)
  repeated <- vapply(seq_along(pts)[-1L], function(i) {
    identical(pts[[i]], pts[[i - 1L]], num.eq = FALSE, single.NA = FALSE)
  }, TRUE)
  pts[c(TRUE, !repeated)]
}

# An objective that logs each point it evaluates.
oc_logged <- function(f) {
  env <- new.env()
  env$pts <- list()
  env$f <- function(theta) {
    env$pts[[length(env$pts) + 1L]] <- theta
    f(theta)
  }
  env
}

# run_optimizer("optim") (or, with a map, .optimize_in_coordinates()) on a
# logged objective.
oc_run <- function(key, f, theta0, lower = NULL, upper = NULL,
                   control = list(), map = NULL) {
  obj <- oc_logged(f)
  val <- p2_val(if (is.null(map)) {
    choicer:::run_optimizer("optim", theta0, obj$f, lower = lower,
                            upper = upper, control = control)
  } else {
    choicer:::.optimize_in_coordinates(map, "optim", theta0, obj$f,
                                       lower = lower, upper = upper,
                                       control = control)
  })
  OC_CALLS[[key]] <<- length(obj$pts)
  list(result = val, points = p2_leaf(oc_collapse(obj$pts)))
}

# A whole fit with its conditions (messages without the run time); kernel
# calls counted through the traced gradient kernels.
oc_fit_record <- function(key, expr) {
  n0 <- oc_counter$n
  warns <- character(); msgs <- character(); err <- NA_character_
  fit <- tryCatch(
    withCallingHandlers(expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      }),
    error = function(e) {
      err <<- sprintf("%s | %s", class(e)[1], conditionMessage(e))
      NULL
    })
  OC_CALLS[[key]] <<- oc_counter$n - n0
  list(fit = if (is.null(fit)) NULL else p2_leaf(p2_strip(fit)),
       warnings = warns,
       messages = msgs[!startsWith(msgs, "Optimization run time")],
       error = err)
}

oc_mxl_lognormal_data <- function() {
  sim <- simulate_mxl_data(N = 150L, J = 3L, seed = 62L, outside_option = FALSE,
                           rc_dist = c(0L, 1L), mu = c(0, -1),
                           Sigma = diag(c(1, 0.3)))
  dt <- data.table::as.data.table(sim$data)
  dt[, w2 := abs(w2) + 0.1]
  dt
}

# ---------------------------------------------------------------------------
# optim_ro
# ---------------------------------------------------------------------------

run_optim_ro <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  # MNL: BFGS by default (no bounds)
  for (ioo in c(TRUE, FALSE)) {
    sim <- simulate_mnl_data(N = 300, J = 4, seed = 11000 + ioo,
                             outside_option = ioo)
    d <- if (ioo) {
      prepare_mnl_data(sim$data, "id", "alt", "choice", c("x1", "x2"),
                       outside_opt_label = 0L, include_outside_option = TRUE)
    } else {
      prepare_mnl_data(sim$data, "id", "alt", "choice", c("x1", "x2"))
    }
    f <- function(theta, ...) choicer:::mnl_loglik_gradient_parallel(
      theta, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, TRUE,
      d$include_outside_option, ...)
    n <- 2L + if (ioo) 4L else 3L
    th0 <- rep(0, n)
    ctls <- list(
      bfgs = list(),
      bfgs_tight = list(reltol = 1e-12, maxit = 500L),
      lbfgsb = list(method = "L-BFGS-B"),
      cg = list(method = "CG", maxit = 300L),
      nm = list(method = "Nelder-Mead", maxit = 400L),
      scaled = list(parscale = seq(0.5, 2, length.out = n), fnscale = 3))
    for (nm in names(ctls)) {
      key <- sprintf("ro|mnl|oo=%d|%s", ioo, nm)
      out[[key]] <- oc_run(key, f, th0, control = ctls[[nm]])
    }
    # The "bhhh" map through .optimize_in_coordinates()
    map <- list(scale = choicer:::.bhhh_scales(f(th0, opg_diag = TRUE)$opg_diag),
                shift = rep(0, n))
    key <- sprintf("ro|mnl|oo=%d|bhhh_map", ioo)
    out[[key]] <- oc_run(key, f, th0, map = map)
  }

  # MXL: L-BFGS-B, as run_mxlogit() passes bounds; the sentinel wraps the
  # kernel as in run_mxlogit().
  sim <- simulate_mxl_data(N = 150L, J = 3L, T = 2L, seed = 11010L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  for (panel in c(FALSE, TRUE)) {
    d <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                          c("w1", "w2"), person_col = if (panel) "pid")
    U <- if (panel) length(unique(dt$pid)) else d$N
    for (draws in c("store", "generate")) {
      eta <- if (draws == "store") get_halton_normals(10L, U, 2L) else
        array(0, c(0L, 0L, 0L))
      for (corr in c(FALSE, TRUE)) {
        f <- choicer:::.lift_sentinel(function(theta, ...) {
          a <- list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
                    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
                    rc_dist = c(0L, 0L), rc_correlation = corr,
                    rc_mean = FALSE, eta_draws = eta, use_asc = TRUE,
                    include_outside_option = FALSE,
                    gen_seed = if (draws == "generate") 5L else -1L,
                    gen_scramble = 1L,
                    gen_S = if (draws == "generate") 10L else 0L, ...)
          if (!is.null(d$Ti)) a$Ti <- d$Ti
          do.call(choicer:::mxl_loglik_gradient_parallel, a)
        })
        n <- 2L + (if (corr) 3L else 2L) + 2L
        th0 <- c(0, 0, rep(log(0.5), if (corr) 3L else 2L), 0, 0)
        if (corr) th0[4L] <- 0
        base <- sprintf("ro|mxl|%s|%s|corr=%d", if (panel) "panel" else "xsec",
                        draws, corr)
        free <- list(lower = rep(-Inf, n), upper = rep(Inf, n))
        bind <- list(lower = rep(-Inf, n), upper = c(0.5, rep(Inf, n - 1L)))
        out[[paste0(base, "|free")]] <- do.call(oc_run, c(
          list(paste0(base, "|free"), f, th0), free))
        out[[paste0(base, "|bind")]] <- do.call(oc_run, c(
          list(paste0(base, "|bind"), f, th0), bind))
        if (!panel && draws == "store") {
          out[[paste0(base, "|bfgs")]] <- oc_run(paste0(base, "|bfgs"), f, th0)
          out[[paste0(base, "|cg")]] <- oc_run(
            paste0(base, "|cg"), f, th0, control = list(method = "CG",
                                                        maxit = 200L))
          g0 <- f(th0)  # a map from the start's gradient
          map <- list(scale = 2^round(log2(1 / pmax(abs(g0$gradient), 1e-3))),
                      shift = rep(0, n))
          out[[paste0(base, "|map")]] <- do.call(oc_run, c(
            list(paste0(base, "|map"), f, th0, map = map), free))
        }
      }
    }
  }
  # The sentinel on a path that steps into the overflow (a log-normal
  # coefficient's mean and log-sd started high): x1, x2, Mu_w1, Mu_w2, L_11,
  # L_22, ASC_2, ASC_3.
  dl <- oc_mxl_lognormal_data()
  d <- prepare_mxl_data(dl, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"))
  eta <- get_halton_normals(10L, d$N, 2L)
  ov <- new.env(); ov$n <- 0L
  f <- choicer:::.lift_sentinel(function(theta) {
    r <- choicer:::mxl_loglik_gradient_parallel(
      theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
      choice_idx = d$choice_idx, M = d$M, weights = d$weights,
      rc_dist = c(0L, 1L), rc_correlation = FALSE, rc_mean = TRUE,
      eta_draws = eta, use_asc = TRUE, include_outside_option = FALSE)
    if (isTRUE(r$overflow)) ov$n <- ov$n + 1L
    r
  })
  for (th in list(c(0, 0, 0, 6, 0, 3, 0, 0), c(0, 0, 0, 5, 0, 4.5, 0, 0),
                  c(0, 0, 0, 3, 0, 5.5, 0, 0))) {
    key <- sprintf("ro|mxl|sentinel|mu=%g|l=%g", th[4], th[6])
    n0 <- ov$n
    out[[key]] <- oc_run(key, f, th, lower = rep(-Inf, 8), upper = rep(Inf, 8))
    out[[key]]$overflow_seen <- ov$n > n0
    cat(sprintf("  %s: %d overflowing evaluations\n", key, ov$n - n0))
  }

  # NL: L-BFGS-B with the nest parameters' lower bound.
  simn <- simulate_nl_data(N = 300, seed = 11020)
  dn <- prepare_nl_data(simn$data, "id", "j", "choice", c("X", "W"), "nest")
  f <- function(theta, ...) choicer:::nl_loglik_gradient_parallel(
    theta, dn$X, dn$alt_idx, dn$choice_idx, dn$nest_idx, dn$M, dn$weights,
    TRUE, dn$include_outside_option, ...)
  # X, W, two nest parameters (nests of two or more), five constants
  n <- 2L + 2L + (length(unique(dn$alt_idx)) - 1L)
  th0 <- c(0, 0, rep(0.5, 2L), rep(0, n - 4L))
  lb <- c(-Inf, -Inf, rep(1e-16, 2L), rep(-Inf, n - 4L))
  out[["ro|nl|lbfgsb"]] <- tryCatch(
    oc_run("ro|nl|lbfgsb", f, th0, lower = lb),
    error = function(e) list(error = conditionMessage(e)))
  out[["ro|nl|lbfgsb_tight"]] <- tryCatch(
    oc_run("ro|nl|lbfgsb_tight", f, th0, lower = lb,
           control = list(factr = 1e3, pgtol = 1e-10)),
    error = function(e) list(error = conditionMessage(e)))
  out[["ro|nl|meta"]] <- list(n = n)
  out
}

# ---------------------------------------------------------------------------
# optim_fit
# ---------------------------------------------------------------------------

run_optim_fit <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  kernels <- c("mnl_loglik_gradient_parallel", "mxl_loglik_gradient_parallel",
               "nl_loglik_gradient_parallel")
  for (k in kernels) {
    suppressMessages(trace(k, tracer = quote(assign("n", oc_counter$n + 1L, envir = oc_counter)),
                           where = ns, print = FALSE))
  }
  on.exit(for (k in kernels) suppressMessages(untrace(k, where = ns)), add = TRUE)

  # MNL (BFGS)
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 11100, outside_option = TRUE)
  dm <- data.table::as.data.table(sim$data)
  dm[, x1 := x1 * 100]
  dm[, w := 0.5 + (id %% 4L) / 2]
  dm[, clus := (id - 1L) %/% 10L]
  mnl <- function(...) list(dm, "id", "alt", "choice", c("x1", "x2"),
                            outside_opt_label = 0L,
                            include_outside_option = TRUE,
                            optimizer = "optim", ...)
  for (sv in c("bhhh", "none", "sd", "mad", "iqr")) {
    key <- sprintf("fit|mnl|%s|hessian", sv)
    out[[key]] <- oc_fit_record(key, do.call("run_mnlogit", mnl(scale_vars = sv)))
  }
  for (se in c("bhhh", "sandwich", "cluster")) {
    key <- sprintf("fit|mnl|bhhh|%s|w", se)
    out[[key]] <- oc_fit_record(key, do.call("run_mnlogit", mnl(
      se_method = se, weights_col = "w",
      cluster_col = if (se == "cluster") "clus")))
  }
  ctls <- list(tight = list(reltol = 1e-12, maxit = 500L),
               lbfgsb = list(method = "L-BFGS-B"),
               cg = list(method = "CG", maxit = 300L))
  for (nm in names(ctls)) {
    key <- sprintf("fit|mnl|none|%s", nm)
    out[[key]] <- oc_fit_record(key, do.call("run_mnlogit", mnl(
      scale_vars = "none", control = ctls[[nm]])))
  }
  sim0 <- simulate_mnl_data(N = 300, J = 4, seed = 11101, outside_option = FALSE)
  out[["fit|mnl|noasc"]] <- oc_fit_record("fit|mnl|noasc", do.call(
    "run_mnlogit", list(sim0$data, "id", "alt", "choice", c("x1", "x2"),
                        use_asc = FALSE, optimizer = "optim")))

  # MXL (L-BFGS-B)
  dx <- p6b_mxl_data(200L, 1L, FALSE, 11110L)
  mxl <- function(...) list(data = dx, id_col = "id", alt_col = "alt",
                            choice_col = "choice",
                            covariate_cols = c("x1", "x2"),
                            random_var_cols = c("w1", "w2"), S = 10L,
                            optimizer = "optim", ...)
  for (sv in c("bhhh", "none", "sd", "mad")) {
    key <- sprintf("fit|mxl|xsec|%s|hessian", sv)
    out[[key]] <- oc_fit_record(key, do.call("run_mxlogit", mxl(scale_vars = sv)))
  }
  out[["fit|mxl|xsec|bhhh|corr|sandwich|w"]] <- oc_fit_record(
    "fit|mxl|xsec|bhhh|corr|sandwich|w", do.call("run_mxlogit", mxl(
      rc_correlation = TRUE, se_method = "sandwich", weights_col = "w")))
  out[["fit|mxl|xsec|bhhh|cluster"]] <- oc_fit_record(
    "fit|mxl|xsec|bhhh|cluster", do.call("run_mxlogit", mxl(
      se_method = "cluster", cluster_col = "cl")))
  out[["fit|mxl|xsec|bhhh|upper"]] <- oc_fit_record(
    "fit|mxl|xsec|bhhh|upper", do.call("run_mxlogit", mxl(
      upper = c(x1 = 0.005), lower = c(L_11 = -6))))
  for (nm in c("factr", "bfgs")) {
    ctl <- if (nm == "factr") list(factr = 1e5, maxit = 300L) else
      list(method = "BFGS")
    key <- sprintf("fit|mxl|xsec|none|%s", nm)
    out[[key]] <- oc_fit_record(key, do.call("run_mxlogit", mxl(
      scale_vars = "none", control = ctl)))
  }
  dp <- p6b_mxl_data(80L, 3L, TRUE, 11120L)
  for (sv in c("bhhh", "none")) {
    key <- sprintf("fit|mxl|panel_gen|%s", sv)
    out[[key]] <- oc_fit_record(key, do.call("run_mxlogit", list(
      data = dp, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      person_col = "pid", rc_dist = c(0L, 1L), rc_mean = TRUE, S = 10L,
      draws = "generate", seed = 3L, include_outside_option = TRUE,
      scale_vars = sv, optimizer = "optim")))
  }
  # The sentinel's path (ro|mxl|sentinel|mu=6|l=3 as a fit)
  dl <- oc_mxl_lognormal_data()
  for (sv in c("none", "bhhh")) {
    key <- sprintf("fit|mxl|sentinel|%s", sv)
    out[[key]] <- oc_fit_record(key, do.call("run_mxlogit", list(
      dl, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"),
      rc_dist = c(0L, 1L), rc_mean = TRUE, S = 10L,
      theta_init = c(0, 0, 0, 6, 0, 3, 0, 0), scale_vars = sv,
      optimizer = "optim")))
  }

  # NL (L-BFGS-B: the nest parameters' lower bound)
  simn <- simulate_nl_data(N = 400, seed = 11130)
  dn <- data.table::as.data.table(simn$data)
  dn[, clus := (id - 1L) %/% 10L]
  nl <- function(...) list(dn, "id", "j", "choice", c("X", "W"), "nest",
                           optimizer = "optim", ...)
  for (sv in c("none", "bhhh", "sd")) {
    key <- sprintf("fit|nl|%s|hessian", sv)
    out[[key]] <- oc_fit_record(key, do.call("run_nestlogit", nl(scale_vars = sv)))
  }
  for (se in c("numeric", "bhhh", "sandwich", "cluster")) {
    key <- sprintf("fit|nl|none|%s", se)
    out[[key]] <- oc_fit_record(key, do.call("run_nestlogit", nl(
      se_method = se, cluster_col = if (se == "cluster") "clus")))
  }
  out[["fit|nl|none|tight"]] <- oc_fit_record("fit|nl|none|tight", do.call(
    "run_nestlogit", nl(control = list(factr = 1e3, pgtol = 1e-10))))
  out
}

oc_save_calls <- function() {
  path <- file.path(.out_dir, sprintf("optim_calls%s_%s_%s.rds", ref_name,
                                         basename(.lib), mode))
  saveRDS(OC_CALLS, path)
  cat(sprintf("  optim kernel calls: %d over %d entries (saved %s)\n",
              sum(unlist(OC_CALLS)), length(OC_CALLS), path))
}
