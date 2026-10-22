# ============================================================================
# Phase 8A batteries (scale audit, the MNL/NL variance kernels: compact
# Hessian and BHHH blocks, the cluster-meat kernels, the variance routes),
# sourced by refactor_identity_check.R after p5b_batteries.R. Every key
# exists on main and on the Phase 8 builds: where a build lacks a cluster-meat
# kernel, the meat is computed as main's R route computes it,
# crossprod(rowsum(w * S, cluster)).
#   * p8_mnl : mnl_loglik_hessian_parallel(), mnl_bhhh_parallel() (weights w
#              and w^2), mnl_scores_parallel(), the gradient with its OPG
#              diagonal, and cluster meats for three label sets, on
#              many-alternative designs (J 6-200, sets of 1-40 straddling 16,
#              K = 1 and K > 1), +- outside option, +- constants, a repeated
#              alternative, the reference never offered.
#   * p8_nl  : the NL kernels likewise, the analytic Hessian included, on
#              nest structures with singletons, many nests and nests a
#              situation only partly offers, dissimilarities below and above 1.
#   * p8_edge: non-finite cells, one kind of non-finite value per call: the
#              chosen utility at -Inf, an unchosen one at +Inf, a huge
#              covariate times a tiny coefficient, an NA / NaN / Inf / zero /
#              negative weight, NA in X, and for NL a NaN dissimilarity.
#   * p8_fit : run_mnlogit() and run_nestlogit() fits on many-alternative data,
#              recorded field by field (estimates, log-likelihood, the fit's
#              variance, every post-hoc vcov(type = ), the lazy route,
#              wesml_vcov()), with keep_data = FALSE twins.
# Compare with the tooling's p8/compare.R (per-entry expected classes).
# ============================================================================

p8_ns <- asNamespace("choicer")
p8_fn <- function(name) get0(name, envir = p8_ns, inherits = FALSE)

p8_theta <- function(seed, K, n_delta, mid = numeric()) {
  set.seed(seed)
  c(rnorm(K, sd = 0.4), mid, rnorm(n_delta, sd = 0.4))
}

p8_nl_theta <- function(seed, K, nest_idx, n_delta) {
  set.seed(seed)
  K_l <- sum(tabulate(nest_idx) > 1L)
  c(rnorm(K, sd = 0.4), runif(K_l, 0.5, 1.5), rnorm(n_delta, sd = 0.4))
}

# The weight-free U x n score matrix of either model.
p8_scores <- function(model, theta, d, asc, ioo, nest_idx = NULL) {
  if (model == "mnl") {
    choicer:::mnl_scores_parallel(theta, d$X, d$alt_idx, d$choice_idx, d$M,
                                  asc, ioo)
  } else {
    choicer:::nl_scores_parallel(theta, d$X, d$alt_idx, d$choice_idx,
                                 nest_idx, d$M, asc, ioo)
  }
}

# Cluster meat sum_g s_g s_g', s_g the cluster's weighted score sum: the
# kernel where the build has one, main's R route otherwise.
p8_meat_cluster <- function(model, theta, d, asc, ioo, cl, nest_idx = NULL) {
  kern <- p8_fn(paste0(model, "_cluster_meat_parallel"))
  if (!is.null(kern)) {
    codes <- match(as.character(cl), unique(as.character(cl)))
    args <- list(theta = theta, X = d$X, alt_idx = d$alt_idx,
                 choice_idx = d$choice_idx)
    if (model == "nl") args$nest_idx <- nest_idx
    args <- c(args, list(M = d$M, weights = d$weights, cluster = codes,
                         use_asc = asc, include_outside_option = ioo))
    return(do.call(kern, args))
  }
  S <- p8_scores(model, theta, d, asc, ioo, nest_idx)
  crossprod(rowsum(d$weights * S, group = as.character(cl)))
}

p8_labels <- function(seed, N) {
  set.seed(seed)
  list(singletons = seq_len(N), coarse = sample(1:5, N, replace = TRUE),
       few = rep(c("b", "a"), length.out = N))
}

p8_mnl_cells <- function(theta, d, asc, ioo, labels) {
  args <- list(theta, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, asc, ioo)
  args_w2 <- args
  args_w2[[6L]] <- d$weights^2
  res <- list(
    grad = p2_val(do.call(choicer:::mnl_loglik_gradient_parallel,
                          c(args, list(opg_diag = TRUE)))),
    hess = p2_val(do.call(choicer:::mnl_loglik_hessian_parallel, args)),
    bhhh = p2_val(do.call(choicer:::mnl_bhhh_parallel, args)),
    bhhh_w2 = p2_val(do.call(choicer:::mnl_bhhh_parallel, args_w2)),
    scores = p2_val(p8_scores("mnl", theta, d, asc, ioo)))
  for (lb in names(labels)) {
    res[[paste0("meat_", lb)]] <-
      p2_val(p8_meat_cluster("mnl", theta, d, asc, ioo, labels[[lb]]))
  }
  res
}

p8_nl_cells <- function(theta, d, nest_idx, asc, ioo, labels) {
  args <- list(theta, d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights,
               asc, ioo)
  args_w2 <- args
  args_w2[[7L]] <- d$weights^2
  res <- list(
    grad = p2_val(do.call(choicer:::nl_loglik_gradient_parallel,
                          c(args, list(opg_diag = TRUE)))),
    hess = p2_val(do.call(choicer:::nl_loglik_hessian_parallel, args)),
    bhhh = p2_val(do.call(choicer:::nl_bhhh_parallel, args)),
    bhhh_w2 = p2_val(do.call(choicer:::nl_bhhh_parallel, args_w2)),
    scores = p2_val(p8_scores("nl", theta, d, asc, ioo, nest_idx)))
  for (lb in names(labels)) {
    res[[paste0("meat_", lb)]] <- p2_val(
      p8_meat_cluster("nl", theta, d, asc, ioo, labels[[lb]], nest_idx))
  }
  res
}

# Designs: J alternatives, U situations with sets of lo..hi alternatives
# (lo < hi: p4_design() draws the set size with sample(lo:hi, 1)).
P8_MNL_DESIGNS <- list(
  sparse = list(U = 60L, J = 30L, lo = 3L, hi = 8L, K = 2L, dup = TRUE),
  wide = list(U = 80L, J = 200L, lo = 10L, hi = 40L, K = 3L, dup = TRUE),
  k1 = list(U = 60L, J = 30L, lo = 3L, hi = 8L, K = 1L, dup = FALSE),
  single = list(U = 60L, J = 12L, lo = 1L, hi = 3L, K = 2L, dup = FALSE),
  nearfull = list(U = 40L, J = 6L, lo = 5L, hi = 6L, K = 3L, dup = FALSE))

run_p8_mnl <- function() {
  out <- list()
  i <- 0L
  for (nm in names(P8_MNL_DESIGNS)) {
    g <- P8_MNL_DESIGNS[[nm]]
    for (ioo in c(FALSE, TRUE)) {
      for (asc in c(TRUE, FALSE)) {
        for (dup in if (g$dup) c(FALSE, TRUE) else FALSE) {
          i <- i + 1L
          d <- p4_design(12000L + i, g$U, g$J, g$lo, g$hi, g$K, ioo, FALSE,
                         dup = dup)
          theta <- p8_theta(12100L + i, g$K, n_delta_free(g$J, asc, ioo))
          key <- sprintf("mnl|%s|oo=%d|asc=%d|dup=%d", nm, ioo, asc, dup)
          out[[key]] <- p8_mnl_cells(theta, d, asc, ioo,
                                     p8_labels(12200L + i, d$N))
        }
      }
    }
  }
  # The reference alternative (code 1) never offered, with constants.
  for (ioo in c(FALSE, TRUE)) {
    d <- p4_design(12300L + ioo, 50L, 20L, 3L, 7L, 2L, ioo, FALSE,
                   absent_ref = TRUE)
    theta <- p8_theta(12310L + ioo, 2L, n_delta_free(20L, TRUE, ioo))
    out[[sprintf("mnl|absent_ref|oo=%d", ioo)]] <-
      p8_mnl_cells(theta, d, TRUE, ioo, p8_labels(12320L + ioo, d$N))
  }
  out
}

# Nest structures over J alternatives (codes 1..J): `nests` holds each
# alternative's nest.
P8_NL_DESIGNS <- list(
  mix = list(U = 60L, J = 30L, lo = 3L, hi = 8L, K = 2L,
             nests = c(rep(1:8, each = 3L), 9:14), dup = TRUE),
  many = list(U = 80L, J = 60L, lo = 10L, hi = 30L, K = 3L,
              nests = c(rep(1:20, each = 2L), rep(21:25, each = 4L)),
              dup = FALSE),
  partly = list(U = 60L, J = 40L, lo = 3L, hi = 6L, K = 2L,
                nests = rep(1:10, each = 4L), dup = FALSE),
  k1 = list(U = 60L, J = 30L, lo = 3L, hi = 8L, K = 1L,
            nests = c(rep(1:8, each = 3L), 9:14), dup = FALSE),
  single = list(U = 60L, J = 12L, lo = 1L, hi = 3L, K = 2L,
                nests = c(rep(1:3, each = 3L), 4:6), dup = FALSE))

run_p8_nl <- function() {
  out <- list()
  i <- 0L
  for (nm in names(P8_NL_DESIGNS)) {
    g <- P8_NL_DESIGNS[[nm]]
    nest_idx <- as.integer(g$nests)
    for (ioo in c(FALSE, TRUE)) {
      for (asc in c(TRUE, FALSE)) {
        for (dup in if (g$dup) c(FALSE, TRUE) else FALSE) {
          i <- i + 1L
          d <- p4_design(13000L + i, g$U, g$J, g$lo, g$hi, g$K, ioo, FALSE,
                         dup = dup)
          theta <- p8_nl_theta(13100L + i, g$K, nest_idx,
                               n_delta_free(g$J, asc, ioo))
          key <- sprintf("nl|%s|oo=%d|asc=%d|dup=%d", nm, ioo, asc, dup)
          out[[key]] <- p8_nl_cells(theta, d, nest_idx, asc, ioo,
                                    p8_labels(13200L + i, d$N))
        }
      }
    }
  }
  out
}

# Non-finite cells on one sparse design per model and outside-option
# setting; situation 3's rows are edited (its choice moved off the outside
# option first).
run_p8_edge <- function() {
  out <- list()
  K <- 2L; J <- 30L
  nest_idx <- as.integer(c(rep(1:8, each = 3L), 9:14))
  for (model in c("mnl", "nl")) {
    for (ioo in c(FALSE, TRUE)) {
      seed <- 14000L + 10L * (model == "nl") + ioo
      d0 <- p4_design(seed, 40L, J, 3L, 8L, K, ioo, FALSE)
      th0 <- if (model == "mnl") {
        p8_theta(seed + 100L, K, n_delta_free(J, TRUE, ioo))
      } else {
        p8_nl_theta(seed + 100L, K, nest_idx, n_delta_free(J, TRUE, ioo))
      }
      t <- 3L
      if (d0$choice_idx[t] == 0L) d0$choice_idx[t] <- 1L
      row0 <- c(0L, cumsum(d0$M))
      chosen_row <- row0[t] + d0$choice_idx[t]
      other_row <- row0[t] + (if (d0$choice_idx[t] == 1L) 2L else 1L)
      cells <- list(
        chosen_minus_inf = function(d, th) {
          d$X[chosen_row, 1L] <- -1e308; th[1L] <- 2; list(d, th)
        },
        unchosen_plus_inf = function(d, th) {
          d$X[other_row, 1L] <- 1e308; th[1L] <- 2; list(d, th)
        },
        huge_x = function(d, th) {
          d$X[chosen_row, 2L] <- 1e200; th[2L] <- 1e-200; list(d, th)
        },
        na_weight = function(d, th) { d$weights[t] <- NA_real_; list(d, th) },
        nan_weight = function(d, th) { d$weights[t] <- NaN; list(d, th) },
        inf_weight = function(d, th) { d$weights[t] <- Inf; list(d, th) },
        zero_weight = function(d, th) { d$weights[t] <- 0; list(d, th) },
        neg_weight = function(d, th) { d$weights[t] <- -1; list(d, th) },
        na_x = function(d, th) { d$X[other_row, 1L] <- NA_real_; list(d, th) })
      if (model == "nl") {
        cells$nan_lambda <- function(d, th) { th[K + 1L] <- NaN; list(d, th) }
      }
      labels <- p8_labels(seed + 200L, d0$N)["coarse"]
      for (nm in names(cells)) {
        e <- cells[[nm]](d0, th0)
        key <- sprintf("edge|%s|%s|oo=%d", model, nm, ioo)
        out[[key]] <- if (model == "mnl") {
          p8_mnl_cells(e[[2L]], e[[1L]], TRUE, ioo, labels)
        } else {
          p8_nl_cells(e[[2L]], e[[1L]], nest_idx, TRUE, ioo, labels)
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# p8_fit: fits on many-alternative data, recorded field by field
# ---------------------------------------------------------------------------

# Long data with choices drawn from the model at a known parameter: N
# situations, sets of lo..hi of J alternatives (codes 1..J), x1 and x2, a
# per-situation weight and cluster. With `nests`, an NL model (dissimilarity
# 0.6 in multi-alternative nests), else an MNL; with `ioo` the outside
# option (utility 0, its own nest) can be chosen and has no rows.
p8_frame <- function(seed, N, J, lo, hi, ioo, nests = NULL) {
  set.seed(seed)
  M <- sample(lo:hi, N, replace = TRUE)
  alt <- unlist(lapply(M, function(m) sort(sample.int(J, m))))
  id <- rep(seq_len(N), M)
  X <- matrix(rnorm(sum(M) * 2L), ncol = 2L)
  delta <- c(if (!ioo) 0, rnorm(J - (!ioo), sd = 0.5))
  beta <- c(0.8, -0.5)
  if (is.null(nests)) {
    V <- drop(X %*% beta) + delta[alt]
    choice <- integer(length(alt))
    start <- c(0L, cumsum(M))
    for (i in seq_len(N)) {
      v <- c(if (ioo) 0, V[(start[i] + 1L):start[i + 1L]])
      p <- exp(v - max(v)); p <- p / sum(p)
      k <- sample.int(length(p), 1L, prob = p) - ioo
      if (k >= 1L) choice[start[i] + k] <- 1L
    }
  } else {
    K_l <- sum(tabulate(nests) > 1L)
    theta <- c(beta, rep(0.6, K_l), if (ioo) delta else delta[-1L])
    pr <- choicer:::nl_predict(theta, X, as.integer(alt), as.integer(M),
                               as.integer(nests), TRUE, ioo)$choice_prob
    choice <- integer(length(alt))
    start <- c(0L, cumsum(M))
    for (i in seq_len(N)) {
      p <- pr[(start[i] + 1L):start[i + 1L]]
      p <- c(if (ioo) max(0, 1 - sum(p)), p)
      k <- sample.int(length(p), 1L, prob = p) - ioo
      if (k >= 1L) choice[start[i] + k] <- 1L
    }
  }
  w_i <- runif(N, 0.5, 2)
  dt <- data.table::data.table(id = id, alt = alt, choice = choice,
                               x1 = X[, 1L], x2 = X[, 2L], w = w_i[id],
                               clus = (id - 1L) %/% 15L)
  if (!is.null(nests)) dt[, nest := nests[alt]]
  dt[]
}

# A fit's variance record: estimates and the fit's variance, every post-hoc
# type (the cluster type with the fit's own labels, or the frame's cluster
# per situation), the lazy route (vcov and se slots set to NULL), whether the
# post-hoc route of the fit's se_method and the lazy route equal the fit's
# variance, and wesml_vcov() for a weighted fit.
p8_vcov_record <- function(fit, cl = NULL, wesml = FALSE) {
  rec <- list(coef = fit$coefficients, loglik = fit$loglik,
              convergence = fit$convergence,
              iterations = fit$optimizer$iterations, vcov = fit$vcov,
              se = fit$se, names = names(fit),
              variance_error = fit$variance_error)
  if (!is.null(fit[["data"]])) {
    for (ty in c("hessian", "bhhh", "robust")) {
      rec[[paste0("ph_", ty)]] <- p2_val(vcov(fit, type = ty))
    }
    rec$ph_cluster <- if (identical(fit$se_method, "cluster")) {
      p2_val(vcov(fit, type = "cluster"))
    } else if (!is.null(cl)) {
      p2_val(vcov(fit, type = "cluster", cluster = cl))
    }
    lazy_fit <- fit
    lazy_fit["vcov"] <- list(NULL)
    lazy_fit["se"] <- list(NULL)
    rec$lazy <- p2_val(vcov(lazy_fit))
    own <- rec[[paste0("ph_", p6d_se_type(fit$se_method))]]$value
    rec$id_posthoc <- identical(unname(own), unname(fit$vcov))
    rec$id_lazy <- identical(unname(rec$lazy$value), unname(fit$vcov))
    if (wesml) rec$wesml <- p2_val(wesml_vcov(fit))
  }
  rec
}

p8_fit_record <- function(expr, cl = NULL, wesml = FALSE) {
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
  list(record = if (!is.null(fit)) p2_rt(p8_vcov_record(fit, cl, wesml)),
       warnings = warns, messages = msgs, error = err)
}

# Per-situation cluster labels named by situation id.
p8_situation_labels <- function(dt) {
  first <- !duplicated(dt$id)
  stats::setNames(dt$clus[first], as.character(dt$id[first]))
}

run_p8_fit <- function() {
  out <- list()
  for (ioo in c(FALSE, TRUE)) {
    dm <- p8_frame(15000L + ioo, 600L, 30L, 4L, 8L, ioo)
    cl <- p8_situation_labels(dm)
    base <- list(dm, "id", "alt", "choice", c("x1", "x2"), weights_col = "w")
    if (ioo) base$include_outside_option <- TRUE
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      args <- c(base, list(se_method = se,
                           cluster_col = if (se == "cluster") "clus"))
      out[[sprintf("fit|mnl|oo=%d|%s", ioo, se)]] <- p8_fit_record(
        do.call("run_mnlogit", args), cl = cl, wesml = se == "sandwich")
    }
    for (se in c("hessian", "sandwich", "cluster")) {
      args <- c(base, list(se_method = se, keep_data = FALSE,
                           cluster_col = if (se == "cluster") "clus"))
      out[[sprintf("fit|mnl|oo=%d|%s|nodata", ioo, se)]] <-
        p8_fit_record(do.call("run_mnlogit", args))
    }
    out[[sprintf("fit|mnl|oo=%d|hessian|noasc", ioo)]] <- p8_fit_record(
      do.call("run_mnlogit", c(base, list(use_asc = FALSE))), cl = cl)
  }
  nests <- as.integer(c(rep(1:8, each = 3L), 9:14))
  for (ioo in c(FALSE, TRUE)) {
    # Seeds whose fits converge at interior dissimilarities (with the
    # outside option, many seeds stop at nloptr's generic failure or a
    # dissimilarity at its lower bound: weakly identified nests).
    dn <- p8_frame(if (ioo) 15106L else 15100L, 800L, 30L, 4L, 8L, ioo,
                   nests = nests)
    cl <- p8_situation_labels(dn)
    base <- list(dn, "id", "alt", "choice", c("x1", "x2"), "nest",
                 weights_col = "w")
    if (ioo) base$include_outside_option <- TRUE
    ses <- if (ioo) c("hessian", "sandwich") else
      c("hessian", "numeric", "bhhh", "sandwich", "cluster")
    for (se in ses) {
      args <- c(base, list(se_method = se,
                           cluster_col = if (se == "cluster") "clus"))
      out[[sprintf("fit|nl|oo=%d|%s", ioo, se)]] <- p8_fit_record(
        do.call("run_nestlogit", args), cl = cl, wesml = se == "sandwich")
    }
    if (!ioo) {
      for (se in c("hessian", "cluster")) {
        args <- c(base, list(se_method = se, keep_data = FALSE,
                             cluster_col = if (se == "cluster") "clus"))
        out[[sprintf("fit|nl|oo=%d|%s|nodata", ioo, se)]] <-
          p8_fit_record(do.call("run_nestlogit", args))
      }
    }
  }
  out
}
