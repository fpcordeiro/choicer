# ============================================================================
# Phase 4 batteries (scale audit: compact Hessian / BHHH, the cluster meat
# kernel, the MXL score-based variance routes), sourced by
# refactor_identity_check.R after p3_batteries.R. Every key exists on main
# and on the Phase 4 builds: where a build lacks the meat kernel, the meat is
# computed as main's R route computes it (crossprod over the score matrix).
#   * mxl_asc : the five estimation kernels on many-alternative designs
#               (J = 30, choice sets of 5-10), the shapes the compact blocks
#               are about: ± outside option, cross-section and panel, store
#               and generate mode, normal and log-normal, forced draw
#               batches, S = 1, use_asc = FALSE, duplicate alternatives, the
#               reference alternative absent from a unit, extra trailing
#               deltas, a single unit.
#   * mxl_meat: robust meat (BHHH with w^2; main's vcov path:
#               crossprod(w * S)) and cluster meats (main:
#               crossprod(rowsum(w * S, cl))) for several label sets.
#   * mxl_vcov: small run_mxlogit() fits (coef identical; vcov / se under
#               tolerance), eager se_method and post-hoc vcov(type = ).
#   * mxl_edge_asc: non-finite cases on many-alternative data: overflowed
#               units (BHHH / score NaN rows), zero-weight -Inf draws (NaN
#               Hessian with a finite likelihood), a huge covariate (Inf
#               inside the block). Non-finite patterns must match.
# Slots hess / bhhh / meat / bread / vcov / se are the tolerance group with
# --tol-slots=hess,bhhh,meat,bread,vcov,se; grad / scores / tastes must then be
# identical.
# ============================================================================

p4_has_meat <- function() {
  exists("mxl_cluster_meat_parallel", envir = asNamespace("choicer"),
         inherits = FALSE)
}

# A many-alternative design. N situations (or sum(Ti) with U units), J
# alternatives, sets of m_lo..m_hi distinct alternatives drawn without the
# reference when `absent_ref` (alternative 1 then never appears), with one
# duplicated row per situation when `dup`, codes drawn from 1..J_codes.
p4_design <- function(seed, U, J, m_lo, m_hi, K_x, ioo, panel,
                      dup = FALSE, absent_ref = FALSE, J_codes = J) {
  set.seed(seed)
  Ti <- if (panel) sample(c(1L, 2L, 3L, 5L), U, replace = TRUE) else rep(1L, U)
  N <- sum(Ti)
  pool <- if (absent_ref) 2:J_codes else seq_len(J_codes)
  alts <- lapply(seq_len(N), function(t) {
    m <- sample(m_lo:m_hi, 1L)
    a <- sort(pool[sample.int(length(pool), min(m, length(pool)))])
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    a
  })
  M <- lengths(alts)
  alt_idx <- as.integer(unlist(alts))
  total <- sum(M)
  X <- matrix(rnorm(total * K_x), total, K_x)
  lo <- if (ioo) 0L else 1L
  choice_idx <- vapply(M, function(m) sample(lo:m, 1L), integer(1))
  w_u <- runif(U, 0.5, 2)
  list(X = X, alt_idx = alt_idx, choice_idx = as.integer(choice_idx),
       M = as.integer(M), weights = rep(w_u, Ti),
       Ti = if (panel) as.integer(Ti) else NULL, U = U, N = N,
       total = total, J = J, K = K_x)
}

p4_theta <- function(seed, K_x, K_w, rc_corr, rc_mean, n_delta) {
  set.seed(seed)
  L_size <- if (rc_corr) (K_w * (K_w + 1L)) %/% 2L else K_w
  c(rnorm(K_x, sd = 0.4), if (rc_mean) rnorm(K_w, sd = 0.3),
    rnorm(L_size, sd = 0.3) - 0.5, rnorm(n_delta, sd = 0.4))
}

p4_modes <- function(S, n_units, K_w) {
  list(store = list(eta = get_halton_normals(S, n_units, K_w), gs = -1L,
                    sc = 1L, gS = 0L),
       gen0 = list(eta = array(0, dim = c(K_w, 0L, 0L)), gs = 13L, sc = 0L,
                   gS = S),
       gen1 = list(eta = array(0, dim = c(K_w, 0L, 0L)), gs = 13L, sc = 1L,
                   gS = S))
}

# The five estimation kernels on one cell. `db` forces draw batches.
p4_kernels <- function(theta, d, W, md, rc_dist, rc_corr, rc_mean, use_asc,
                       ioo, db = NULL) {
  extra <- list(gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS)
  if (!is.null(d$Ti)) extra$Ti <- d$Ti
  if (!is.null(db)) extra$draw_batch <- as.integer(db)
  wt <- list(theta, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, md$eta,
             rc_dist, rc_corr, rc_mean, use_asc, ioo)
  nw <- wt[-7L]
  list(
    grad = do.call(mxl_loglik_gradient_parallel, c(wt, extra)),
    hess = do.call(mxl_hessian_parallel, c(wt, extra)),
    bhhh = do.call(mxl_bhhh_parallel, c(wt, extra)),
    scores = do.call(mxl_scores_parallel, c(nw, extra)),
    tastes = do.call(mxl_conditional_tastes_parallel, c(nw, extra))
  )
}

run_mxl_asc <- function() {
  out <- list()
  J <- 30L; K_x <- 2L; K_w <- 2L; S <- 12L
  dists <- list(normal = c(0L, 0L), mixed = c(0L, 1L))
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 14L else 40L
      d <- p4_design(9100 + 2 * panel + ioo, U, J, 5L, 10L, K_x, ioo, panel)
      set.seed(9200 + 2 * panel + ioo)
      W <- matrix(rnorm(d$total * K_w), d$total, K_w)
      modes <- p4_modes(S, U, K_w)
      for (mode_name in names(modes)) {
        for (dist_name in names(dists)) {
          for (rc_mean in c(TRUE, FALSE)) {
            for (rc_corr in c(TRUE, FALSE)) {
              theta <- p4_theta(9300 + rc_mean + 2 * rc_corr, K_x, K_w, rc_corr,
                                rc_mean, n_delta_free(J, TRUE, ioo))
              key <- sprintf("asc|%s|oo=%d|%s|%s|mean=%d|corr=%d",
                             if (panel) "panel" else "xsec", ioo, mode_name,
                             dist_name, rc_mean, rc_corr)
              out[[key]] <- p4_kernels(theta, d, W, modes[[mode_name]],
                                       dists[[dist_name]], rc_corr, rc_mean,
                                       TRUE, ioo)
            }
          }
        }
      }
      # Forced draw batches (pass 2 forms one situation's rows at a time).
      theta <- p4_theta(9301, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, ioo))
      for (db in c(1L, 5L)) {
        key <- sprintf("asc|%s|oo=%d|store|batch=%d",
                       if (panel) "panel" else "xsec", ioo, db)
        out[[key]] <- p4_kernels(theta, d, W, modes$store, dists$mixed, TRUE,
                                 TRUE, TRUE, ioo, db = db)
      }
      # use_asc = FALSE: the compact block is the whole parameter vector.
      theta <- p4_theta(9302, K_x, K_w, TRUE, TRUE, 0L)
      key <- sprintf("asc|%s|oo=%d|store|noasc", if (panel) "panel" else "xsec",
                     ioo)
      out[[key]] <- p4_kernels(theta, d, W, modes$store, dists$mixed, TRUE,
                               TRUE, FALSE, ioo)
    }
  }
  # Shapes: S = 1, duplicate alternatives, the reference absent, extra
  # trailing deltas (codes below J), a single unit, units < threads.
  shapes <- list(
    s1 = list(S = 1L), dup = list(dup = TRUE), absent = list(absent_ref = TRUE),
    trailing = list(J_codes = 24L), one_unit = list(U = 1L),
    two_units = list(U = 2L))
  for (sh in names(shapes)) {
    for (panel in c(FALSE, TRUE)) {
      for (ioo in c(FALSE, TRUE)) {
        cfg <- shapes[[sh]]
        U <- cfg$U %||% (if (panel) 10L else 25L)
        S_c <- cfg$S %||% S
        d <- p4_design(9400 + 2 * panel + ioo, U, J, 5L, 10L, K_x, ioo, panel,
                       dup = isTRUE(cfg$dup), absent_ref = isTRUE(cfg$absent_ref),
                       J_codes = cfg$J_codes %||% J)
        set.seed(9500 + 2 * panel + ioo)
        W <- matrix(rnorm(d$total * K_w), d$total, K_w)
        theta <- p4_theta(9303, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, ioo))
        md <- p4_modes(S_c, U, K_w)
        for (mode_name in c("store", "gen1")) {
          key <- sprintf("shape|%s|%s|oo=%d|%s", sh,
                         if (panel) "panel" else "xsec", ioo, mode_name)
          out[[key]] <- p4_kernels(theta, d, W, md[[mode_name]], c(0L, 1L),
                                   TRUE, TRUE, TRUE, ioo)
        }
      }
    }
  }
  out
}

# Cluster meat: the kernel where it exists, else main's R route.
p4_meat_cluster <- function(theta, d, W, md, rc_dist, ioo, cl) {
  codes <- match(as.character(cl), unique(as.character(cl)))
  extra <- list(gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS)
  if (!is.null(d$Ti)) extra$Ti <- d$Ti
  if (p4_has_meat()) {
    return(do.call(choicer:::mxl_cluster_meat_parallel, c(
      list(theta, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, codes,
           md$eta, rc_dist, TRUE, TRUE, TRUE, ioo), extra)))
  }
  S <- do.call(mxl_scores_parallel, c(
    list(theta, d$X, W, d$alt_idx, d$choice_idx, d$M, md$eta, rc_dist, TRUE,
         TRUE, TRUE, ioo), extra))
  u <- if (is.null(d$Ti)) seq_len(d$N) else cumsum(c(1L, d$Ti[-length(d$Ti)]))
  crossprod(rowsum(d$weights[u] * S, group = as.character(cl)))
}

run_mxl_meat <- function() {
  out <- list()
  J <- 30L; K_x <- 2L; K_w <- 2L; S <- 12L
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 16L else 48L
      d <- p4_design(9600 + 2 * panel + ioo, U, J, 5L, 10L, K_x, ioo, panel)
      set.seed(9700 + 2 * panel + ioo)
      W <- matrix(rnorm(d$total * K_w), d$total, K_w)
      theta <- p4_theta(9304, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, ioo))
      md_all <- p4_modes(S, U, K_w)
      set.seed(9800 + 2 * panel + ioo)
      labels <- list(singletons = seq_len(U), one = rep(1L, U),
                     coarse = sample(1:5, U, replace = TRUE),
                     few = rep(c("b", "a"), length.out = U),
                     shuffled = sample(seq_len(U)) * 10.5)
      for (mode_name in c("store", "gen1")) {
        md <- md_all[[mode_name]]
        key <- sprintf("meat|%s|oo=%d|%s", if (panel) "panel" else "xsec", ioo,
                       mode_name)
        extra <- list(gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS)
        if (!is.null(d$Ti)) extra$Ti <- d$Ti
        res <- list()
        # Robust meat: main's vcov path forms crossprod(w * S); the Phase 4
        # route calls the BHHH kernel with w^2.
        res$meat <- list(robust = if (p4_has_meat()) {
          do.call(mxl_bhhh_parallel, c(list(theta, d$X, W, d$alt_idx,
            d$choice_idx, d$M, d$weights^2, md$eta, c(0L, 1L), TRUE, TRUE,
            TRUE, ioo), extra))
        } else {
          S_u <- do.call(mxl_scores_parallel, c(list(theta, d$X, W, d$alt_idx,
            d$choice_idx, d$M, md$eta, c(0L, 1L), TRUE, TRUE, TRUE, ioo),
            extra))
          u <- if (is.null(d$Ti)) seq_len(d$N) else
            cumsum(c(1L, d$Ti[-length(d$Ti)]))
          crossprod(d$weights[u] * S_u)
        })
        for (lb in names(labels)) {
          res$meat[[lb]] <- p4_meat_cluster(theta, d, W, md, c(0L, 1L), ioo,
                                            labels[[lb]])
        }
        out[[key]] <- res
      }
    }
  }
  out
}

# Small fits: the variance routes of run_mxlogit() and vcov(type = ).
run_mxl_vcov <- function() {
  out <- list()
  # Well-conditioned fits (bread condition numbers about 20-30 on main), so
  # that vcov is held to the same 1e-10 as its bread and meat.
  for (cfg in list(
    list(name = "xsec", N = 1500L, T = 1L, ioo = FALSE, draws = "store",
         corr = TRUE),
    list(name = "xsec_gen", N = 1500L, T = 1L, ioo = FALSE, draws = "generate",
         corr = FALSE),
    list(name = "panel", N = 150L, T = 3L, ioo = FALSE, draws = "store",
         corr = TRUE),
    list(name = "panel_oo_gen", N = 300L, T = 3L, ioo = TRUE, draws = "generate",
         corr = TRUE))) {
    sim <- simulate_mxl_data(N = cfg$N, J = 8L,
                             T = cfg$T, beta = c(0.8, -0.5),
                             Sigma = if (cfg$T > 1L) diag(c(1.0, 0.8)) else
                               diag(c(1.5, 1.0)),
                             delta = seq(-0.5, 0.5, length.out = 8L),
                             outside_option = cfg$ioo, vary_choice_set = TRUE,
                             seed = 9900L)
    dt <- data.table::as.data.table(sim$data)
    # The simulated outside option is an explicit alt = 0 row; the kernels'
    # implicit outside option needs those rows dropped.
    if (cfg$ioo) dt <- dt[alt != 0]
    person <- if (cfg$T > 1L) "pid" else NULL
    set.seed(9901)
    ids <- unique(dt$id)
    grp <- if (cfg$T > 1L) {
      p <- unique(dt[, c("id", "pid"), with = FALSE])
      setNames(sample(1:12, length(unique(p$pid)), replace = TRUE)[
        match(p$pid, unique(p$pid))], p$id)
    } else setNames(sample(1:12, length(ids), replace = TRUE), ids)
    dt[, cl := grp[as.character(id)]]
    w_id <- if (cfg$T > 1L) {
      pw <- runif(length(unique(dt$pid)), 0.5, 2)
      pw[match(dt$pid, unique(dt$pid))]
    } else runif(length(ids), 0.5, 2)[match(dt$id, ids)]
    dt[, w := w_id]
    common <- list(data = dt, id_col = "id", alt_col = "alt",
                   choice_col = "choice", covariate_cols = c("x1", "x2"),
                   random_var_cols = c("w1", "w2"), person_col = person,
                   rc_correlation = cfg$corr,
                   weights_col = "w", S = 50L, draws = cfg$draws, seed = 7L,
                   include_outside_option = cfg$ioo)
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      args <- common
      args$se_method <- se
      if (se == "cluster") args$cluster_col <- "cl"
      fit <- suppressMessages(suppressWarnings(do.call(run_mxlogit, args)))
      key <- sprintf("vcov|%s|%s", cfg$name, se)
      res <- list(coef = unname(fit$coefficients), vcov = unname(fit$vcov),
                  se = unname(fit$se), loglik = fit$loglik)
      if (se == "hessian") {
        cl <- dt[!duplicated(id), setNames(cl, id)]
        res$posthoc <- list(
          hessian = unname(vcov(fit, type = "hessian")),
          bhhh = unname(vcov(fit, type = "bhhh")),
          robust = unname(vcov(fit, type = "robust")),
          cluster = unname(vcov(fit, type = "cluster", cluster = cl)))
        res$bread <- unname(choicer:::.compute_bread(fit))
      }
      out[[key]] <- res
    }
  }
  out
}

# Non-finite cases on many-alternative data (store mode, one unit edited).
run_mxl_edge_asc <- function() {
  out <- list()
  J <- 30L; K_x <- 2L; K_w <- 2L; S <- 12L
  for (panel in c(FALSE, TRUE)) {
    for (ioo in c(FALSE, TRUE)) {
      U <- if (panel) 12L else 30L
      d0 <- p4_design(9950 + 2 * panel + ioo, U, J, 5L, 10L, K_x, ioo, panel)
      set.seed(9960 + 2 * panel + ioo)
      W0 <- matrix(rnorm(d0$total * K_w), d0$total, K_w)
      theta <- p4_theta(9305, K_x, K_w, TRUE, TRUE, n_delta_free(J, TRUE, ioo))
      row0 <- c(0L, cumsum(d0$M))
      first <- if (is.null(d0$Ti)) seq_len(d0$N) else
        cumsum(c(1L, d0$Ti[-length(d0$Ti)]))
      u <- 3L
      t <- first[u]
      if (d0$choice_idx[t] == 0L) d0$choice_idx[t] <- 1L
      chosen_row <- row0[t] + d0$choice_idx[t]
      md <- p4_modes(S, U, K_w)$store
      # (1) the chosen fixed utility overflows: the unit's likelihood is 0.
      d <- d0; th <- theta; th[1L] <- 2
      d$X[chosen_row, 1L] <- -1e308
      out[[sprintf("edge|overflow|%s|oo=%d", if (panel) "panel" else "xsec",
                   ioo)]] <- p4_kernels(th, d, W0, md, c(0L, 0L), TRUE, TRUE,
                                        TRUE, ioo)
      # (2) zero-weight -Inf draws: draws 1-5 of unit u give its choice
      # probability 0, the others do not (NaN Hessian, finite likelihood).
      d <- d0; W <- W0
      t_end <- t + (if (is.null(d$Ti)) 1L else d$Ti[u])
      rows_u <- (row0[t] + 1L):row0[t_end]
      W[rows_u, 1L] <- 0
      W[chosen_row, 1L] <- 10
      md2 <- md
      md2$eta[1L, 1:5, u] <- -1e308
      out[[sprintf("edge|zero_weight_inf|%s|oo=%d",
                   if (panel) "panel" else "xsec", ioo)]] <-
        p4_kernels(theta, d, W, md2, c(0L, 0L), TRUE, TRUE, TRUE, ioo)
      # (3) a huge covariate with a tiny coefficient: finite utility, Inf
      # second moments inside the unit's block.
      d <- d0; th <- theta; th[2L] <- 1e-200
      d$X[chosen_row, 2L] <- 1e200
      out[[sprintf("edge|huge_x|%s|oo=%d", if (panel) "panel" else "xsec",
                   ioo)]] <- p4_kernels(th, d, W0, md, c(0L, 0L), TRUE, TRUE,
                                        TRUE, ioo)
      # (4) a non-finite weight (direct calls only; R rejects it).
      d <- d0
      d$weights[t] <- NaN
      if (!is.null(d$Ti)) d$weights[first[u] - 1L + seq_len(d$Ti[u])] <- NaN
      out[[sprintf("edge|nan_weight|%s|oo=%d", if (panel) "panel" else "xsec",
                   ioo)]] <- p4_kernels(theta, d, W0, md, c(0L, 0L), TRUE,
                                        TRUE, TRUE, ioo)[c("hess", "bhhh")]
    }
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
