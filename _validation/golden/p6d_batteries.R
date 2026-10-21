# ============================================================================
# Phase 6D batteries (scale audit: scale_vars as a change of the optimizer's
# coordinates, the BHHH-diagonal option), sourced by refactor_identity_check.R
# after p6b_batteries.R. Version-aware: main computes each reference the way
# main can, a 6D build the new way; p6d/compare.R holds each commit's
# expected class per entry (main has no tolerance group for MNL/NL).
#   * p6d_diag : the gradient kernels' OPG diagonal (formal `opg_diag`)
#                against diag() of the BHHH kernel with the same weights
#                (main), and the objective/gradient a flag-on call returns
#                (main: the flag-off call). MXL store/gen0/gen1 x panel /
#                cross-section x outside option x normal/mixed x
#                rc_mean/rc_corr, forced draw batches, two non-finite cells;
#                MNL/NL x outside option x use_asc, a wide design, repeated
#                alternatives; finite non-uniform weights.
#   * p6d_coord: the optimizer's objective and gradient at a fixed
#                theta-tilde for sd/mad/iqr: main's scaled evaluation (MXL
#                kernels with sX/sW, MNL on sweep()ed X) against a 6D build's
#                natural kernel at theta = m * theta-tilde + c with the
#                gradient times m (.coordinate_objective()); and the map
#                (m, c): main's run_mxlogit() arithmetic against
#                .coordinate_map().
#   * p6d_vcov : small fits per model x se_method x scale: coefficients,
#                fit$vcov, se, vcov(fit, type = ), the lazy ensure_vcov()
#                route (vcov dropped), wesml_vcov() on weighted fits, their
#                identity flags, param_scale / param_shift; badly scaled
#                designs whose natural Hessian is singular to solve().
#   * p6d_new  : "bhhh" fits (MXL, MNL, NL) and NL scale_vars fits: errors on
#                main (match.arg, unused argument), values on a 6D build.
# All calls run at one thread (the harness sets it).
# ============================================================================

p6d_has_formal <- function(fn, arg) arg %in% names(formals(fn))
p6d_has_coord <- function() {
  exists(".coordinate_objective", envir = asNamespace("choicer"), inherits = FALSE)
}
p6d_num <- function(x) as.numeric(x)

# ---------------------------------------------------------------------------
# p6d_diag
# ---------------------------------------------------------------------------

# One gradient call with the OPG diagonal: the 6D kernel's own output, or on
# main the flag-off gradient plus diag() of the BHHH kernel with the same
# weights (B = sum_u w_u s_u s_u').
p6d_diag_call <- function(grad_fn, bhhh_fn, args) {
  if (p6d_has_formal(grad_fn, "opg_diag")) {
    g <- do.call(grad_fn, c(args, list(opg_diag = TRUE)))
    dg <- p6d_num(g$opg_diag)
    g$opg_diag <- NULL
  } else {
    g <- do.call(grad_fn, args)
    dg <- p6d_num(diag(do.call(bhhh_fn, args)))
  }
  list(diag = dg, grad = g)
}

p6d_mxl_args <- function(theta, d, W, md, rc_dist, rc_corr, rc_mean, ioo,
                         db = NULL) {
  a <- list(theta = theta, X = d$X, W = W, alt_idx = d$alt_idx,
            choice_idx = d$choice_idx, M = d$M, weights = d$weights,
            eta_draws = md$eta, rc_dist = rc_dist, rc_correlation = rc_corr,
            rc_mean = rc_mean, use_asc = TRUE, include_outside_option = ioo,
            gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS)
  if (!is.null(d$Ti)) a$Ti <- d$Ti
  if (!is.null(db)) a$draw_batch <- as.integer(db)
  a
}

run_p6d_diag <- function() {
  out <- list()
  J <- 8L; K_x <- 2L; K_w <- 2L; S <- 6L
  dists <- list(normal = c(0L, 0L), mixed = c(0L, 1L))
  rcs <- list(c(mean = TRUE, corr = TRUE), c(mean = FALSE, corr = FALSE),
              c(mean = TRUE, corr = FALSE))
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 10L else 24L
      d <- p4_design(9900 + 2 * panel + ioo, U, J, 3L, 6L, K_x, ioo, panel)
      set.seed(9910 + 2 * panel + ioo)
      W <- matrix(rnorm(d$total * K_w), d$total, K_w)
      modes <- p4_modes(S, U, K_w)
      for (mode_name in names(modes)) {
        for (dist_name in names(dists)) {
          for (rc in rcs) {
            theta <- p4_theta(9920 + rc[["mean"]] + 2 * rc[["corr"]], K_x, K_w,
                              rc[["corr"]], rc[["mean"]],
                              n_delta_free(J, TRUE, ioo))
            key <- sprintf("diag|mxl|%s|oo=%d|%s|%s|mean=%d|corr=%d",
                           if (panel) "panel" else "xsec", ioo, mode_name,
                           dist_name, rc[["mean"]], rc[["corr"]])
            out[[key]] <- p6d_diag_call(
              mxl_loglik_gradient_parallel, mxl_bhhh_parallel,
              p6d_mxl_args(theta, d, W, modes[[mode_name]], dists[[dist_name]],
                           rc[["corr"]], rc[["mean"]], ioo))
          }
        }
      }
      # Forced draw batches of 2 draws.
      theta <- p4_theta(9930, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, ioo))
      for (mode_name in names(modes)) {
        key <- sprintf("diag|mxl|%s|oo=%d|%s|batches",
                       if (panel) "panel" else "xsec", ioo, mode_name)
        out[[key]] <- p6d_diag_call(
          mxl_loglik_gradient_parallel, mxl_bhhh_parallel,
          p6d_mxl_args(theta, d, W, modes[[mode_name]], c(0L, 1L), TRUE, TRUE,
                       ioo, db = 2L))
      }
    }
  }
  # Non-finite cells: a huge covariate (overflowing utilities) and a huge
  # log-normal mean. Recorded with their conditions, never aborting the run.
  d <- p4_design(9940, 24L, J, 3L, 6L, K_x, FALSE, FALSE)
  set.seed(9941)
  W <- matrix(rnorm(d$total * K_w), d$total, K_w)
  md <- p4_modes(S, 24L, K_w)$store
  theta <- p4_theta(9942, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, FALSE))
  d_big <- d
  d_big$X[5L, 1L] <- 1e308
  theta_big <- theta
  theta_big[1L] <- 4  # 4e308 overflows to Inf: the kernel's sentinel
  out[["diag|mxl|edge|bigx"]] <- p2_val(p6d_diag_call(
    mxl_loglik_gradient_parallel, mxl_bhhh_parallel,
    p6d_mxl_args(theta_big, d_big, W, md, c(0L, 0L), TRUE, TRUE, FALSE)))
  theta_ln <- theta
  theta_ln[K_x + 2L] <- 800  # exp(800) overflows
  out[["diag|mxl|edge|bigmu"]] <- p2_val(p6d_diag_call(
    mxl_loglik_gradient_parallel, mxl_bhhh_parallel,
    p6d_mxl_args(theta_ln, d, W, md, c(0L, 1L), TRUE, TRUE, FALSE)))

  # MNL and NL.
  mnl_bhhh <- choicer:::mnl_bhhh_parallel
  nl_bhhh <- choicer:::nl_bhhh_parallel
  mnl_designs <- list(small = list(N = 40, J = 5, K = 3),
                      wide = list(N = 60, J = 30, K = 3))
  for (dn in names(mnl_designs)) {
    g <- mnl_designs[[dn]]
    for (ioo in c(FALSE, TRUE)) {
      dat <- p2_dat(9950 + 10 * (dn == "wide") + ioo, g$N, g$J, g$K, ioo)
      for (use_asc in c(TRUE, FALSE)) {
        set.seed(9960 + ioo + 2 * use_asc)
        theta <- rnorm(g$K + n_delta_free(g$J, use_asc, ioo), sd = 0.4)
        a <- list(theta = theta, X = dat$X, alt_idx = dat$alt_idx,
                  choice_idx = dat$choice_idx, M = dat$M,
                  weights = dat$weights, use_asc = use_asc,
                  include_outside_option = ioo)
        out[[sprintf("diag|mnl|%s|asc=%d|oo=%d", dn, use_asc, ioo)]] <-
          p6d_diag_call(mnl_loglik_gradient_parallel, mnl_bhhh, a)
      }
    }
  }
  # Repeated alternatives within a situation (one row duplicated per set).
  for (ioo in c(FALSE, TRUE)) {
    d <- p4_design(9970 + ioo, 30L, 8L, 3L, 6L, 2L, ioo, FALSE, dup = TRUE)
    set.seed(9972 + ioo)
    theta <- rnorm(2L + n_delta_free(8L, TRUE, ioo), sd = 0.4)
    a <- list(theta = theta, X = d$X, alt_idx = d$alt_idx,
              choice_idx = d$choice_idx, M = d$M, weights = d$weights,
              use_asc = TRUE, include_outside_option = ioo)
    out[[sprintf("diag|mnl|dup|oo=%d", ioo)]] <-
      p2_val(p6d_diag_call(mnl_loglik_gradient_parallel, mnl_bhhh, a))
    nest <- c(1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L)
    th_nl <- p2_nl_theta(2L, nest, 8L, TRUE, ioo, seed = 9974 + ioo)
    a_nl <- list(theta = th_nl, X = d$X, alt_idx = d$alt_idx,
                 choice_idx = d$choice_idx, nest_idx = nest, M = d$M,
                 weights = d$weights, use_asc = TRUE,
                 include_outside_option = ioo)
    out[[sprintf("diag|nl|dup|oo=%d", ioo)]] <-
      p2_val(p6d_diag_call(nl_loglik_gradient_parallel, nl_bhhh, a_nl))
  }
  nl_designs <- list(
    small = list(N = 40, J = 5, K = 3, nest = c(1L, 1L, 2L, 2L, 3L)),
    nests20 = list(N = 60, J = 40, K = 3,
                   nest = as.integer(c(rep(1:15, each = 2), 16L, 16L, 16L,
                                       17L, 17L, 18L, 18L, 19L, 19L, 20L))))
  for (dn in names(nl_designs)) {
    g <- nl_designs[[dn]]
    stopifnot(length(g$nest) == g$J)
    for (ioo in c(FALSE, TRUE)) {
      dat <- p2_dat(9980 + 10 * (dn == "nests20") + ioo, g$N, g$J, g$K, ioo)
      for (use_asc in c(TRUE, FALSE)) {
        theta <- p2_nl_theta(g$K, g$nest, g$J, use_asc, ioo,
                             seed = 9990 + ioo + 2 * use_asc)
        a <- list(theta = theta, X = dat$X, alt_idx = dat$alt_idx,
                  choice_idx = dat$choice_idx, nest_idx = g$nest, M = dat$M,
                  weights = dat$weights, use_asc = use_asc,
                  include_outside_option = ioo)
        out[[sprintf("diag|nl|%s|asc=%d|oo=%d", dn, use_asc, ioo)]] <-
          p6d_diag_call(nl_loglik_gradient_parallel, nl_bhhh, a)
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# p6d_coord
# ---------------------------------------------------------------------------

# Main's map, as run_mxlogit() / run_mnlogit() build bt_mult and bt_shift.
p6d_main_map <- function(param_map, n_params, sX, sW, rc_corr, K_w) {
  m <- rep(1, n_params); sh <- rep(0, n_params)
  if (length(sX)) m[param_map$beta] <- 1 / sX
  if (length(param_map$mu)) m[param_map$mu] <- 1 / sW
  if (K_w > 0L) {
    if (rc_corr) {
      idx <- 1L
      for (i in seq_len(K_w)) {
        for (j in seq_len(i)) {
          pos <- param_map$sigma[idx]
          if (i == j) sh[pos] <- -log(sW[i]) else m[pos] <- 1 / sW[i]
          idx <- idx + 1L
        }
      }
    } else {
      for (i in seq_len(K_w)) sh[param_map$sigma[i]] <- -log(sW[i])
    }
  }
  list(scale = m, shift = sh)
}

p6d_map <- function(method, param_map, n_params, sX, sW, rc_corr, K_w) {
  if (p6d_has_coord()) {
    choicer:::.coordinate_map(method, param_map, n_params, sX = sX, sW = sW,
                              rc_correlation = rc_corr)
  } else {
    p6d_main_map(param_map, n_params, sX, sW, rc_corr, K_w)
  }
}

p6d_mxl_param_map <- function(K_x, K_w, rc_corr, rc_mean, n_delta) {
  L_size <- if (rc_corr) (K_w * (K_w + 1L)) %/% 2L else K_w
  pm <- list(beta = seq_len(K_x))
  pos <- K_x
  if (rc_mean) { pm$mu <- pos + seq_len(K_w); pos <- pos + K_w }
  pm$sigma <- pos + seq_len(L_size); pos <- pos + L_size
  pm$asc <- pos + seq_len(n_delta)
  pm
}

run_p6d_coord <- function() {
  out <- list()
  J <- 8L; K_x <- 2L; K_w <- 2L; S <- 6L
  dists <- list(normal = c(0L, 0L), mixed = c(0L, 1L))
  rcs <- list(c(mean = TRUE, corr = TRUE), c(mean = FALSE, corr = FALSE),
              c(mean = TRUE, corr = FALSE))
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 10L else 24L
      d <- p4_design(10000 + 2 * panel + ioo, U, J, 3L, 6L, K_x, ioo, panel)
      d$X <- d$X %*% diag(c(100, 0.025))  # units the scales must undo
      colnames(d$X) <- c("x1", "x2")
      set.seed(10010 + 2 * panel + ioo)
      W <- matrix(rnorm(d$total * K_w), d$total, K_w) %*% diag(c(25, 0.3))
      W[, 2L] <- abs(W[, 2L]) + 0.1
      colnames(W) <- c("w1", "w2")
      modes <- p4_modes(S, U, K_w)
      n_delta <- n_delta_free(J, TRUE, ioo)
      for (mode_name in names(modes)) {
        md <- modes[[mode_name]]
        for (dist_name in names(dists)) {
          rc_dist <- dists[[dist_name]]
          for (rc in rcs) {
            rc_corr <- rc[["corr"]]; rc_mean <- rc[["mean"]]
            pm <- p6d_mxl_param_map(K_x, K_w, rc_corr, rc_mean, n_delta)
            n_params <- max(unlist(pm))
            th_t <- p4_theta(10020 + rc_mean + 2 * rc_corr, K_x, K_w, rc_corr,
                             rc_mean, n_delta)
            for (method in c("sd", "mad", "iqr")) {
              sX <- choicer:::.column_scales(d$X, method)
              sW <- choicer:::.column_scales(W, method)
              sW[rc_dist == 1L] <- 1
              map <- p6d_map(method, pm, n_params, sX, sW, rc_corr, K_w)
              a <- p6d_mxl_args(NULL, d, W, md, rc_dist, rc_corr, rc_mean, ioo)
              a$theta <- NULL; a$X <- NULL; a$W <- NULL
              natural <- function(th) {
                do.call(mxl_loglik_gradient_parallel,
                        c(list(theta = th, X = d$X, W = W), a))
              }
              val <- if (p6d_has_coord()) {
                choicer:::.coordinate_objective(natural, map)(th_t)
              } else {
                do.call(mxl_loglik_gradient_parallel,
                        c(list(theta = th_t, X = d$X, W = W), a,
                          list(sX = sX, sW = sW)))
              }
              key <- sprintf("coord|mxl|%s|oo=%d|%s|%s|mean=%d|corr=%d|%s",
                             if (panel) "panel" else "xsec", ioo, mode_name,
                             dist_name, rc_mean, rc_corr, method)
              out[[key]] <- list(value = val, scale = map$scale,
                                 shift = map$shift)
            }
          }
        }
      }
    }
  }
  # MNL: main's run_mnlogit() divides X by sweep(); ASCs unscaled.
  for (ioo in c(FALSE, TRUE)) {
    dat <- p2_dat(10040 + ioo, 60, 6, 3, ioo)
    dat$X <- dat$X %*% diag(c(100, 0.025, 7))
    for (use_asc in c(TRUE, FALSE)) {
      set.seed(10050 + ioo + 2 * use_asc)
      n_delta <- n_delta_free(6, use_asc, ioo)
      th_t <- rnorm(3 + n_delta, sd = 0.4)
      pm <- list(beta = 1:3)
      if (n_delta) pm$asc <- 3L + seq_len(n_delta)
      for (method in c("sd", "mad", "iqr")) {
        sX <- choicer:::.column_scales(dat$X, method)
        map <- p6d_map(method, pm, 3L + n_delta, sX, NULL, FALSE, 0L)
        natural <- function(th) {
          mnl_loglik_gradient_parallel(th, dat$X, dat$alt_idx, dat$choice_idx,
                                       dat$M, dat$weights, use_asc, ioo)
        }
        val <- if (p6d_has_coord()) {
          choicer:::.coordinate_objective(natural, map)(th_t)
        } else {
          mnl_loglik_gradient_parallel(th_t, sweep(dat$X, 2, sX, "/"),
                                       dat$alt_idx, dat$choice_idx, dat$M,
                                       dat$weights, use_asc, ioo)
        }
        out[[sprintf("coord|mnl|asc=%d|oo=%d|%s", use_asc, ioo, method)]] <-
          list(value = val, scale = map$scale, shift = map$shift)
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# p6d_vcov
# ---------------------------------------------------------------------------

p6d_se_type <- function(se) {
  switch(se, sandwich = "robust", numeric = "hessian", se)
}

# Fields main's fit objects lack; p2_strip() drops them from whole-fit leaves.
p6d_new_fields <- function(fit) {
  f <- c("param_scale", "param_shift")
  if (inherits(fit, "choicer_nl")) f <- c(f, "scale_vars", "sX")
  f
}

# A fit's variance record: the fit's own vcov/se, the post-hoc route of its
# se_method, the lazy route (vcov dropped, then vcov()), wesml_vcov() for a
# weighted fit, whether they are identical to the fit's, and the new fields.
p6d_vcov_record <- function(fit, wesml = FALSE) {
  rec <- list(coef = fit$coefficients, loglik = fit$loglik,
              convergence = fit$convergence,
              iterations = fit$optimizer$iterations,
              vcov = fit$vcov, se = fit$se,
              names = setdiff(names(fit), p6d_new_fields(fit)),
              new = list(param_scale = fit$param_scale,
                         param_shift = fit$param_shift,
                         scale_vars = fit$scale_vars, sX = fit$sX))
  if (!is.null(fit[["data"]])) {
    ph <- p2_val(vcov(fit, type = p6d_se_type(fit$se_method)))
    lazy_fit <- fit
    lazy_fit$vcov <- NULL
    lazy_fit$se <- NULL
    lz <- p2_val(vcov(lazy_fit))
    rec$posthoc <- ph
    rec$lazy <- lz
    rec$id_posthoc <- identical(unname(ph$value), unname(fit$vcov))
    rec$id_lazy <- identical(unname(lz$value), unname(fit$vcov))
    if (wesml) rec$wesml <- p2_val(wesml_vcov(fit))
  }
  rec
}

# A fit with its conditions (messages without the run time); the record is
# taken from the fit as returned, before any stripping.
p6d_fit_record <- function(expr, wesml = FALSE) {
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
  msgs <- msgs[!startsWith(msgs, "Optimization run time")]
  list(record = if (!is.null(fit)) p2_rt(p6d_vcov_record(fit, wesml)),
       warnings = warns, messages = msgs, error = err)
}

run_p6d_vcov <- function() {
  out <- list()
  # MXL cross-section, stored draws, weights, clusters.
  dt <- p6b_mxl_data(300L, 1L, FALSE, 10100L)
  for (sv in c("none", "sd")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      args <- list(data = dt, id_col = "id", alt_col = "alt",
                   choice_col = "choice", covariate_cols = c("x1", "x2"),
                   random_var_cols = c("w1", "w2"), rc_correlation = TRUE,
                   S = 20L, scale_vars = sv, se_method = se,
                   weights_col = "w")
      if (se == "cluster") args$cluster_col <- "cl"
      out[[sprintf("vc|mxl|xsec_store|%s|%s", sv, se)]] <-
        p6d_fit_record(do.call("run_mxlogit", args), wesml = TRUE)
    }
  }
  # MXL panel, generated draws, outside option, a log-normal coefficient.
  dtp <- p6b_mxl_data(120L, 3L, TRUE, 10110L)
  for (sv in c("none", "mad")) {
    for (se in c("hessian", "sandwich")) {
      args <- list(data = dtp, id_col = "id", alt_col = "alt",
                   choice_col = "choice", covariate_cols = c("x1", "x2"),
                   random_var_cols = c("w1", "w2"), person_col = "pid",
                   rc_dist = c(0L, 1L), rc_mean = TRUE, S = 20L,
                   draws = "generate", seed = 3L,
                   include_outside_option = TRUE, scale_vars = sv,
                   se_method = se)
      out[[sprintf("vc|mxl|panel_gen|%s|%s", sv, se)]] <-
        p6d_fit_record(do.call("run_mxlogit", args))
    }
  }
  # Advanced workflow with draws that are not Halton (post hoc regenerates
  # Halton draws; the fit's own vcov uses the caller's).
  d_adv <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                            c("w1", "w2"))
  set.seed(10120)
  eta <- array(rnorm(2L * 20L * d_adv$N), dim = c(2L, 20L, d_adv$N))
  for (sv in c("none", "sd")) {
    out[[sprintf("vc|mxl|adv_eta|%s", sv)]] <- p6d_fit_record(
      run_mxlogit(input_data = d_adv, eta_draws = eta, scale_vars = sv))
  }
  # keep_data = FALSE keeps a NULL data slot.
  for (sv in c("none", "sd")) {
    out[[sprintf("vc|mxl|nodata|%s", sv)]] <- p6d_fit_record(do.call(
      "run_mxlogit",
      list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
           covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
           S = 20L, scale_vars = sv, keep_data = FALSE)))
  }
  # MNL.
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 10130, outside_option = FALSE)
  dm <- data.table::as.data.table(sim$data)
  dm[, x1 := x1 * 100]
  dm[, clus := (id - 1L) %/% 10L]
  set.seed(10131)
  dm[, w := runif(.N, 0.5, 2)[match(id, unique(id))]]
  for (sv in c("none", "sd")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      args <- list(dm, "id", "alt", "choice", c("x1", "x2"), se_method = se,
                   scale_vars = sv, weights_col = "w",
                   cluster_col = if (se == "cluster") "clus")
      out[[sprintf("vc|mnl|%s|%s", sv, se)]] <-
        p6d_fit_record(do.call("run_mnlogit", args), wesml = TRUE)
    }
  }
  # NL (no scale_vars on main).
  simn <- simulate_nl_data(N = 400, seed = 10140)
  dn <- data.table::as.data.table(simn$data)
  dn[, clus := (id - 1L) %/% 10L]
  for (se in c("hessian", "numeric", "bhhh", "sandwich", "cluster")) {
    out[[sprintf("vc|nl|none|%s", se)]] <- p6d_fit_record(do.call(
      "run_nestlogit",
      list(dn, "id", "j", "choice", c("X", "W"), "nest", se_method = se,
           cluster_col = if (se == "cluster") "clus")))
  }
  # Badly scaled designs: the natural Hessian is singular to solve().
  simb <- simulate_mxl_data(N = 300, J = 4L, beta = c(0.8, -0.5),
                            Sigma = diag(c(1, 0.6)), seed = 10150)
  db <- data.table::as.data.table(simb$data)
  db[, x1 := x1 * 1e7]
  db[, x2 := x2 * 1e-7]
  for (sv in c("none", "sd")) {
    for (se in c("hessian", "sandwich")) {
      out[[sprintf("vc|mxl|bad|%s|%s", sv, se)]] <- p6d_fit_record(do.call(
        "run_mxlogit",
        list(data = db, id_col = "id", alt_col = "alt", choice_col = "choice",
             covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
             S = 20L, scale_vars = sv, se_method = se,
             control = list(maxeval = 300L))))
    }
  }
  simbm <- simulate_mnl_data(N = 300, J = 4, seed = 10160, outside_option = FALSE)
  dbm <- data.table::as.data.table(simbm$data)
  dbm[, x1 := x1 * 1e7]
  dbm[, x2 := x2 * 1e-7]
  for (sv in c("none", "sd")) {
    out[[sprintf("vc|mnl|bad|%s|hessian", sv)]] <- p6d_fit_record(
      do.call("run_mnlogit", list(dbm, "id", "alt", "choice", c("x1", "x2"),
                                  scale_vars = sv)))
  }
  out
}

# ---------------------------------------------------------------------------
# p6d_new
# ---------------------------------------------------------------------------

run_p6d_new <- function() {
  out <- list()
  dt <- p6b_mxl_data(300L, 1L, FALSE, 10200L)
  for (se in c("hessian", "sandwich")) {
    out[[sprintf("new|mxl|xsec_store|bhhh|%s", se)]] <- p6d_fit_record(do.call(
      "run_mxlogit",
      list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
           covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
           rc_correlation = TRUE, S = 20L, scale_vars = "bhhh",
           se_method = se, weights_col = "w")))
  }
  dtp <- p6b_mxl_data(120L, 3L, TRUE, 10210L)
  out[["new|mxl|panel_gen|bhhh|hessian"]] <- p6d_fit_record(do.call(
    "run_mxlogit",
    list(data = dtp, id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
         person_col = "pid", rc_dist = c(0L, 1L), rc_mean = TRUE, S = 20L,
         draws = "generate", seed = 3L, include_outside_option = TRUE,
         scale_vars = "bhhh")))
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 10220, outside_option = TRUE)
  dm <- data.table::as.data.table(sim$data)
  dm[, x1 := x1 * 100]
  for (se in c("hessian", "bhhh")) {
    out[[sprintf("new|mnl|bhhh|%s", se)]] <- p6d_fit_record(do.call(
      "run_mnlogit",
      list(dm, "id", "alt", "choice", c("x1", "x2"), se_method = se,
           scale_vars = "bhhh", outside_opt_label = 0L,
           include_outside_option = TRUE)))
  }
  simn <- simulate_nl_data(N = 400, seed = 10230)
  dn <- data.table::as.data.table(simn$data)
  dn[, X := X * 100]
  for (sv in c("none", "sd", "mad", "iqr", "bhhh")) {
    out[[sprintf("new|nl|%s|hessian", sv)]] <- p6d_fit_record(do.call(
      "run_nestlogit",
      list(dn, "id", "j", "choice", c("X", "W"), "nest", scale_vars = sv)))
  }
  out
}
