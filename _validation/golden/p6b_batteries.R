# ============================================================================
# Phase 6B batteries (scale audit: scale_vars applied on load), sourced by
# refactor_identity_check.R after p4_batteries.R. Every key exists on main and
# on the Phase 6B builds; every entry is compared with identical().
#   * scale_kern: the four kernels a fit calls (gradient, Hessian, BHHH,
#                 cluster meat) at a scaled-space theta. On builds whose
#                 kernels take column scales (formal `sX`) they get the
#                 natural X and W plus sX/sW; on main, X and W divided by
#                 sweep(), as main's run_mxlogit() formed them. Grid:
#                 store/gen0/gen1, cross-section and panel, +/- outside
#                 option, normal / mixed / all-log-normal (log-normal columns
#                 scaled by 1, as the fit passes them), row-aligned and
#                 alternative-level W, rc_mean x rc_corr, data-derived and
#                 awkward scales; extras: forced draw batches, powers of two,
#                 extreme scales, K_x = 0, NA in a log-normal column.
#   * scale_fit : run_mxlogit() fits with scale_vars = sd/mad/iqr on badly
#                 scaled data: the stripped fit (coefficients, vcov, se,
#                 loglik, optimizer iterations, sX, sW, sigma, stored data),
#                 messages (without the run time) and warnings; every
#                 se_method, panels, store and generate, the outside option,
#                 weights, the advanced workflow (also with an
#                 alternative-level W), post-hoc vcov(type =) of one fit.
#   * scale_cols: .column_scales() on assorted matrices (value, names,
#                 errors): apply() on main, a column loop on C2.
#   * scale_mnl : run_mnlogit() fits with scale_vars (C2's MNL path).
# ============================================================================

p6b_has_scales <- function() {
  "sX" %in% names(formals(choicer:::mxl_loglik_gradient_parallel))
}

# One kernel call at a scaled-space theta: the natural X/W and the column
# scales on builds that take them, X/W divided by sweep() otherwise.
p6b_call <- function(fn, X, W, sX, sW, args) {
  if (p6b_has_scales()) {
    do.call(fn, c(list(X = X, W = W), args, list(sX = sX, sW = sW)))
  } else {
    do.call(fn, c(list(X = sweep(X, 2, sX, "/"), W = sweep(W, 2, sW, "/")),
                  args))
  }
}

# The four fit-time kernels on one cell; `db` forces draw batches. With
# `wrap`, each call's value, warnings, messages and error are recorded.
p6b_kernels <- function(theta, d, X, W, sc, md, rc_dist, rc_corr, rc_mean, ioo,
                        cl, db = NULL, wrap = FALSE) {
  a <- list(theta = theta, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
            M = d$M, weights = d$weights, eta_draws = md$eta,
            rc_dist = rc_dist, rc_correlation = rc_corr, rc_mean = rc_mean,
            use_asc = TRUE, include_outside_option = ioo,
            gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS)
  if (!is.null(d$Ti)) a$Ti <- d$Ti
  if (!is.null(db)) a$draw_batch <- as.integer(db)
  run <- function(fn, args) {
    if (wrap) p2_val(p6b_call(fn, X, W, sc$sX, sc$sW, args))
    else p6b_call(fn, X, W, sc$sX, sc$sW, args)
  }
  list(
    grad = run(mxl_loglik_gradient_parallel, a),
    hess = run(mxl_hessian_parallel, a),
    bhhh = run(mxl_bhhh_parallel, a),
    meat = run(choicer:::mxl_cluster_meat_parallel, c(a, list(cluster = cl)))
  )
}

p6b_col_sd <- function(M) {
  vapply(seq_len(ncol(M)), function(k) stats::sd(M[, k]), numeric(1))
}

# Column scales as the fit passes them: log-normal W columns get 1.
p6b_scales <- function(kind, X, W, rc_dist) {
  pick <- function(v, n) v[seq_len(n)]
  sX <- switch(kind,
               sd = p6b_col_sd(X),
               odd = pick(c(1 / 3, 7321.5, 0.0029), ncol(X)),
               pow2 = pick(c(0.5, 4, 2), ncol(X)),
               extreme = pick(c(1e-7, 1e7, 3e-5), ncol(X)))
  sW <- switch(kind,
               sd = p6b_col_sd(W),
               odd = pick(c(0.0029, 3, 1 / 7), ncol(W)),
               pow2 = pick(c(2, 0.25, 8), ncol(W)),
               extreme = pick(c(1e-3, 1e3, 0.37), ncol(W)))
  sW[rc_dist == 1L] <- 1
  list(sX = sX, sW = sW)
}

run_scale_kern <- function() {
  out <- list()
  J <- 8L; K_x <- 2L; K_w <- 2L; S <- 6L
  dists <- list(normal = c(0L, 0L), mixed = c(0L, 1L), lognormal = c(1L, 1L))
  rcs <- list(c(mean = TRUE, corr = TRUE), c(mean = FALSE, corr = FALSE),
              c(mean = TRUE, corr = FALSE))
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 10L else 24L
      d <- p4_design(9600 + 2 * panel + ioo, U, J, 3L, 6L, K_x, ioo, panel)
      set.seed(9610 + 2 * panel + ioo)
      W_row <- matrix(rnorm(d$total * K_w), d$total, K_w) %*% diag(c(5, 0.2))
      W_alt <- matrix(rnorm(J * K_w), J, K_w)
      cl <- as.integer((seq_len(U) - 1L) %% 4L + 1L)
      modes <- p4_modes(S, U, K_w)
      for (mode_name in names(modes)) {
        for (dist_name in names(dists)) {
          for (wl in c("row", "alt")) {
            W <- if (wl == "row") W_row else W_alt
            for (rc in rcs) {
              theta <- p4_theta(9620 + rc[["mean"]] + 2 * rc[["corr"]], K_x, K_w,
                                rc[["corr"]], rc[["mean"]],
                                n_delta_free(J, TRUE, ioo))
              for (sk in c("sd", "odd")) {
                sc <- p6b_scales(sk, d$X, W, dists[[dist_name]])
                key <- sprintf("kern|%s|oo=%d|%s|%s|W=%s|mean=%d|corr=%d|s=%s",
                               if (panel) "panel" else "xsec", ioo, mode_name,
                               dist_name, wl, rc[["mean"]], rc[["corr"]], sk)
                out[[key]] <- p6b_kernels(theta, d, d$X, W, sc, modes[[mode_name]],
                                          dists[[dist_name]], rc[["corr"]],
                                          rc[["mean"]], ioo, cl)
              }
            }
          }
        }
      }
    }
  }

  # Extras on a cross-section and a panel without the outside option.
  for (panel in c(FALSE, TRUE)) {
    U <- if (panel) 10L else 24L
    d <- p4_design(9640 + panel, U, J, 3L, 6L, K_x, FALSE, panel)
    set.seed(9650 + panel)
    W <- matrix(rnorm(d$total * K_w), d$total, K_w)
    cl <- as.integer((seq_len(U) - 1L) %% 3L + 1L)
    modes <- p4_modes(S, U, K_w)
    theta <- p4_theta(9660, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, FALSE))
    tag <- if (panel) "panel" else "xsec"
    for (mode_name in names(modes)) {
      md <- modes[[mode_name]]
      for (sk in c("pow2", "extreme")) {
        sc <- p6b_scales(sk, d$X, W, c(0L, 1L))
        out[[sprintf("extra|%s|%s|s=%s", tag, mode_name, sk)]] <-
          p6b_kernels(theta, d, d$X, W, sc, md, c(0L, 1L), TRUE, TRUE, FALSE, cl)
      }
      # Forced draw batches (2 draws per batch).
      sc <- p6b_scales("odd", d$X, W, c(0L, 1L))
      out[[sprintf("extra|%s|%s|batches", tag, mode_name)]] <-
        p6b_kernels(theta, d, d$X, W, sc, md, c(0L, 1L), TRUE, TRUE, FALSE, cl,
                    db = 2L)
      # NA in a log-normal W column: divided by its scale of 1, as sweep() did.
      W_na <- W
      W_na[3L, 2L] <- NA
      out[[sprintf("extra|%s|%s|na_lognormal", tag, mode_name)]] <-
        p6b_kernels(theta, d, d$X, W_na, sc, md, c(0L, 1L), TRUE, TRUE, FALSE, cl)
    }
    # K_x = 0: no fixed coefficients, an empty sX. The kernels refuse K_x = 0
    # (parse_mxl_theta); the error must come before any scale check.
    d0 <- p4_design(9670 + panel, U, J, 3L, 6L, 0L, FALSE, panel)
    set.seed(9680 + panel)
    W0 <- matrix(rnorm(d0$total * K_w), d0$total, K_w)
    theta0 <- p4_theta(9690, 0L, K_w, TRUE, FALSE, n_delta_free(J, TRUE, FALSE))
    modes0 <- p4_modes(S, U, K_w)
    for (mode_name in names(modes0)) {
      sc <- p6b_scales("odd", d0$X, W0, c(0L, 0L))
      out[[sprintf("extra|%s|%s|kx0", tag, mode_name)]] <-
        p6b_kernels(theta0, d0, d0$X, W0, sc, modes0[[mode_name]], c(0L, 0L),
                    TRUE, FALSE, FALSE, cl, wrap = TRUE)
    }
  }
  out
}

# A fit's recorded outcome: the stripped value, warnings and messages without
# the optimizer's run time (it varies run to run). Fits are called by name, so
# that the stored call holds the symbol, not the function (whose body differs
# between builds).
p6b_fit_val <- function(expr) {
  v <- p2_val(expr)
  v$messages <- v$messages[!startsWith(v$messages, "Optimization run time")]
  v
}

# Simulated MXL data with badly scaled columns, weights and cluster labels
# constant within decision maker (or situation).
p6b_mxl_data <- function(N, T, ioo, seed) {
  sim <- simulate_mxl_data(N = N, J = 5L, T = T, beta = c(0.8, -0.5),
                           Sigma = diag(c(1.0, 0.6)),
                           delta = seq(-0.4, 0.4, length.out = 5L),
                           outside_option = ioo, vary_choice_set = TRUE,
                           seed = seed)
  dt <- data.table::as.data.table(sim$data)
  # The simulated outside option is an explicit alt = 0 row; the kernels'
  # implicit outside option needs those rows dropped.
  if (ioo) dt <- dt[alt != 0]
  dt[, x1 := x1 * 100]
  dt[, x2 := x2 / 40]
  dt[, w1 := w1 * 25]
  dt[, w2 := abs(w2) + 0.1]  # positive: a reasonable log-normal regressor
  unit <- if (T > 1L) dt$pid else dt$id
  u <- unique(unit)
  set.seed(seed + 1L)
  dt[, w := runif(length(u), 0.5, 2)[match(unit, u)]]
  dt[, cl := sample(1:9, length(u), replace = TRUE)[match(unit, u)]]
  dt
}

run_scale_fit <- function() {
  out <- list()
  cfgs <- list(
    list(name = "xsec_store", N = 400L, T = 1L, ioo = FALSE, draws = "store",
         rc_dist = c(0L, 0L), corr = TRUE, mean = FALSE, w = FALSE,
         scales = c("sd", "mad", "iqr"),
         se = c("hessian", "bhhh", "sandwich", "cluster")),
    list(name = "xsec_gen_mixed", N = 400L, T = 1L, ioo = FALSE,
         draws = "generate", rc_dist = c(0L, 1L), corr = FALSE, mean = TRUE,
         w = TRUE, scales = c("sd", "mad", "iqr"), se = "hessian"),
    list(name = "panel_store", N = 120L, T = 3L, ioo = FALSE, draws = "store",
         rc_dist = c(0L, 0L), corr = TRUE, mean = TRUE, w = TRUE,
         scales = c("sd", "mad"), se = c("hessian", "sandwich", "cluster")),
    list(name = "panel_gen_lognormal", N = 120L, T = 3L, ioo = TRUE,
         draws = "generate", rc_dist = c(1L, 1L), corr = FALSE, mean = TRUE,
         w = FALSE, scales = "sd", se = c("hessian", "bhhh")),
    list(name = "xsec_oo_store", N = 400L, T = 1L, ioo = TRUE, draws = "store",
         rc_dist = c(0L, 0L), corr = TRUE, mean = TRUE, w = FALSE,
         scales = "iqr", se = c("bhhh", "cluster")))
  for (cfg in cfgs) {
    dt <- p6b_mxl_data(cfg$N, cfg$T, cfg$ioo, 9700L + cfg$N + cfg$T + cfg$ioo)
    for (sv in cfg$scales) {
      for (se in cfg$se) {
        args <- list(data = dt, id_col = "id", alt_col = "alt",
                     choice_col = "choice", covariate_cols = c("x1", "x2"),
                     random_var_cols = c("w1", "w2"),
                     person_col = if (cfg$T > 1L) "pid",
                     rc_dist = cfg$rc_dist, rc_correlation = cfg$corr,
                     rc_mean = cfg$mean, S = 30L, draws = cfg$draws, seed = 7L,
                     include_outside_option = cfg$ioo, scale_vars = sv,
                     se_method = se, control = list(maxeval = 300L))
        if (cfg$w) args$weights_col <- "w"
        if (se == "cluster") args$cluster_col <- "cl"
        key <- sprintf("fit|%s|%s|%s", cfg$name, sv, se)
        res <- list(fit = p6b_fit_val(do.call("run_mxlogit", args)))
        if (cfg$name == "xsec_store" && sv == "sd" && se == "hessian" &&
            !is.null(res$fit$value)) {
          f <- suppressMessages(do.call("run_mxlogit", args))
          res$posthoc <- list(
            robust = p2_val(vcov(f, type = "robust")),
            bhhh = p2_val(vcov(f, type = "bhhh")),
            cluster = p2_val(vcov(f, type = "cluster",
                                  cluster = dt[!duplicated(id), setNames(cl, id)])))
        }
        out[[key]] <- res
      }
    }
  }

  # Advanced workflow: prepared input and a store-mode cube; also with an
  # alternative-level W (one row per global alternative).
  dt <- p6b_mxl_data(400L, 1L, FALSE, 9790L)
  d_adv <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                            c("w1", "w2"), rc_correlation = TRUE)
  eta <- get_halton_normals(30L, d_adv$N, 2L)
  for (sv in c("sd", "mad")) {
    out[[sprintf("fit|adv|%s", sv)]] <- list(fit = p6b_fit_val(
      run_mxlogit(input_data = d_adv, eta_draws = eta, scale_vars = sv,
                  rc_mean = TRUE, control = list(maxeval = 300L))))
  }
  d_alt <- d_adv
  set.seed(9791)
  d_alt$W <- matrix(rnorm(nrow(d_adv$alt_mapping) * 2L) * c(9, 0.3),
                    nrow(d_adv$alt_mapping), 2L,
                    dimnames = list(NULL, c("w1", "w2")))
  out[["fit|adv_altW|sd"]] <- list(fit = p6b_fit_val(
    run_mxlogit(input_data = d_alt, eta_draws = eta, scale_vars = "sd",
                se_method = "bhhh", control = list(maxeval = 300L))))
  out
}

run_scale_cols <- function() {
  out <- list()
  set.seed(9800)
  base <- matrix(rnorm(60 * 3) * rep(c(1, 50, 0.01), each = 60), 60, 3)
  nm3 <- c("a", "b", "c")
  mats <- list(
    plain = base,
    colnames = `colnames<-`(base, nm3),
    rownames = `rownames<-`(base, paste0("r", 1:60)),
    both = `dimnames<-`(base, list(paste0("r", 1:60), nm3)),
    one_row = `colnames<-`(base[1, , drop = FALSE], nm3),
    one_row_both = `dimnames<-`(base[1, , drop = FALSE], list("r1", nm3)),
    zero_row = `colnames<-`(base[0, , drop = FALSE], nm3),
    one_col = base[, 2, drop = FALSE],
    integer = matrix(sample.int(9L, 120L, replace = TRUE), 60, 2,
                     dimnames = list(NULL, c("i1", "i2"))),
    na = local({ m <- base; m[5, 2] <- NA; m }),
    inf = local({ m <- base; m[7, 3] <- Inf; m }),
    constant = cbind(base[, 1], 4),
    dummy = cbind(d = rbinom(60, 1, 0.2), x = base[, 1]))
  for (nm in names(mats)) {
    for (m in c("sd", "mad", "iqr")) {
      out[[sprintf("cols|%s|%s", nm, m)]] <-
        p2_val(choicer:::.column_scales(mats[[nm]], m))
    }
  }
  out
}

run_scale_mnl <- function() {
  out <- list()
  for (ioo in c(TRUE, FALSE)) {
    sim <- simulate_mnl_data(N = 300, J = 4, seed = 9810 + ioo,
                             outside_option = ioo)
    d_use <- data.table::as.data.table(sim$data)
    d_use[, x1 := x1 * 100]
    d_use[, clus := (id - 1L) %/% 10L]
    for (sv in c("sd", "mad", "iqr")) {
      for (se in c("hessian", "sandwich", "cluster")) {
        args <- list(d_use, "id", "alt", "choice", c("x1", "x2"),
                     se_method = se, scale_vars = sv,
                     cluster_col = if (se == "cluster") "clus")
        if (ioo) args <- c(args, list(outside_opt_label = 0L,
                                      include_outside_option = TRUE))
        out[[sprintf("mnl|oo=%d|%s|%s", ioo, sv, se)]] <-
          list(fit = p6b_fit_val(do.call("run_mnlogit", args)))
      }
    }
  }
  out
}
