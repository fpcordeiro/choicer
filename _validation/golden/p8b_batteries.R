# ============================================================================
# Phase 8B batteries (scale audit: the multinomial and nested logit kernels'
# working arrays, their nest loops and slot maps, every model's elasticity
# and diversion-ratio accumulators, the cluster labels' coding), sourced by
# refactor_identity_check.R after p8_batteries.R (whose helpers they use, with
# p2's and p3's). Every call exists with the same R signature on main and on
# the Phase 8B builds.
#   * p8b_mnl : every MNL kernel (gradient with and without its OPG diagonal,
#               Hessian, BHHH with w and w^2, scores, cluster meat, predict,
#               shares, BLP, elasticities of the first and last covariate,
#               diversion) on designs whose temporaries reach the heap (K = 17,
#               sets of 10-30), with sets of exactly 4 and K = 4 (gemv's tiny
#               square path), K = 1 (the operand swap), runs of identical
#               choice sets, unsorted and repeated alternatives; +- outside
#               option, +- constants.
#   * p8b_nl  : the NL kernels likewise, the analytic Hessian and BLP with
#               damping included, on 60 nests of 2 (situations touching a few),
#               odd and even nest counts, singletons, K = 17 and K = 4 designs.
#   * p8b_edge: a zero weight on a situation whose X' d overflows (BLAS and
#               tiny-square paths), infinite utilities inside runs of equal
#               choice sets, NA covariates, and dissimilarities near 1e-16,
#               +Inf and NaN, through every kernel; nl_blp_contraction() with
#               lambda <= 0 (no check there).
#   * p8b_mxl : the MXL elasticity and diversion kernels on wide sets with
#               repeated alternatives, +- outside option, store and generate
#               draws.
#   * p8b_cl  : vcov(type = "cluster") of MNL, NL and MXL (cross-section and
#               panel) fits with cluster labels of every type (integer,
#               logical, character, factor with reordered, unused and NA
#               levels, a factor with duplicated levels, doubles that print
#               alike, dates), and fit-time cluster_col of several types.
# Compare dumps bitwise (the tooling's p8b/compare_bits.R): identical() with
# its defaults equates +0 with -0 and ignores NaN payloads.
# ============================================================================

# N situations over J alternatives (codes 1..J), sets of lo..hi, K
# covariates. `runs`: consecutive situations share their set, in runs of 1-5;
# `unsorted`: each set's rows in random order; `dup`: one alternative listed
# twice. The first situation offers every alternative.
p8b_design <- function(seed, N, J, lo, hi, K, ioo, runs = FALSE,
                       unsorted = FALSE, dup = FALSE, sd_x = 1) {
  set.seed(seed)
  alts <- vector("list", N)
  i <- 1L
  while (i <= N) {
    m <- sample(lo:hi, 1L)
    a <- sort(sample.int(J, min(m, J)))
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    if (unsorted) a <- a[sample.int(length(a))]
    r <- if (runs) sample(1:5, 1L) else 1L
    for (k in seq_len(min(r, N - i + 1L))) {
      alts[[i]] <- a
      i <- i + 1L
    }
  }
  alts[[1L]] <- if (unsorted) sample.int(J) else seq_len(J)
  M <- lengths(alts)
  total <- sum(M)
  lo_c <- if (ioo) 0L else 1L
  list(X = matrix(rnorm(total * K, sd = sd_x), total, K),
       alt_idx = as.integer(unlist(alts)), M = as.integer(M),
       choice_idx = vapply(M, function(m) sample(lo_c:m, 1L), integer(1)),
       weights = runif(N, 0.5, 2), N = N, J = J, K = K)
}

# Every MNL kernel's output on one design.
p8b_mnl_cells <- function(theta, d, asc, ioo, labels, K) {
  args <- list(theta, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, asc, ioo)
  args_w2 <- args
  args_w2[[6L]] <- d$weights^2
  beta <- theta[seq_len(K)]
  res <- list(
    grad = p2_val(do.call(choicer:::mnl_loglik_gradient_parallel, args)),
    grad_opg = p2_val(do.call(choicer:::mnl_loglik_gradient_parallel,
                              c(args, list(opg_diag = TRUE)))),
    hess = p2_val(do.call(choicer:::mnl_loglik_hessian_parallel, args)),
    bhhh = p2_val(do.call(choicer:::mnl_bhhh_parallel, args)),
    bhhh_w2 = p2_val(do.call(choicer:::mnl_bhhh_parallel, args_w2)),
    scores = p2_val(p8_scores("mnl", theta, d, asc, ioo)),
    pred = p2_val(choicer:::mnl_predict(theta, d$X, d$alt_idx, d$M, asc, ioo)),
    shares = p2_val(choicer:::mnl_predict_shares(theta, d$X, d$alt_idx, d$M,
                                                 d$weights, asc, ioo)),
    elas1 = p2_val(choicer:::mnl_elasticities_parallel(
      theta, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, 1L, asc, ioo)),
    elasK = p2_val(choicer:::mnl_elasticities_parallel(
      theta, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, K, asc, ioo)),
    dr = p2_val(choicer:::mnl_diversion_ratios_parallel(
      theta, d$X, d$alt_idx, d$M, d$weights, asc, ioo)))
  if (asc) {
    J <- d$J
    target <- rep(1 / (J + ioo), J + ioo)
    # (uniform targets some designs cannot attain: the engine's check stops
    # them, as the leaf expects from its inputs)
    res$blp <- golden_blp("blp_contraction", list(
      rep(0, J), target, d$X, beta, d$alt_idx, d$M, d$weights, ioo,
      max_iter = 25L))
  }
  for (lb in names(labels)) {
    res[[paste0("meat_", lb)]] <-
      p2_val(p8_meat_cluster("mnl", theta, d, asc, ioo, labels[[lb]]))
  }
  res
}

# Every NL kernel's output on one design (lam: the full dissimilarity vector
# for the BLP contraction).
p8b_nl_cells <- function(theta, d, nest_idx, asc, ioo, labels, K) {
  args <- list(theta, d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights,
               asc, ioo)
  args_w2 <- args
  args_w2[[7L]] <- d$weights^2
  K_l <- sum(tabulate(nest_idx) > 1L)
  lam <- rep(1, max(nest_idx))
  lam[tabulate(nest_idx, max(nest_idx)) > 1L] <- theta[K + seq_len(K_l)]
  res <- list(
    grad = p2_val(do.call(choicer:::nl_loglik_gradient_parallel, args)),
    grad_opg = p2_val(do.call(choicer:::nl_loglik_gradient_parallel,
                              c(args, list(opg_diag = TRUE)))),
    hess = p2_val(do.call(choicer:::nl_loglik_hessian_parallel, args)),
    bhhh = p2_val(do.call(choicer:::nl_bhhh_parallel, args)),
    bhhh_w2 = p2_val(do.call(choicer:::nl_bhhh_parallel, args_w2)),
    scores = p2_val(p8_scores("nl", theta, d, asc, ioo, nest_idx)),
    pred = p2_val(choicer:::nl_predict(theta, d$X, d$alt_idx, d$M, nest_idx,
                                       asc, ioo)),
    shares = p2_val(choicer:::nl_predict_shares(theta, d$X, d$alt_idx, d$M,
                                                d$weights, nest_idx, asc, ioo)),
    elas1 = p2_val(choicer:::nl_elasticities_parallel(
      theta, d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights, 1L, asc,
      ioo)),
    elasK = p2_val(choicer:::nl_elasticities_parallel(
      theta, d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights, K, asc,
      ioo)),
    dr = p2_val(choicer:::nl_diversion_ratios_parallel(
      theta, d$X, d$alt_idx, nest_idx, d$M, d$weights, asc, ioo)))
  if (asc) {
    J <- d$J
    target <- rep(1 / (J + ioo), J + ioo)
    beta <- theta[seq_len(K)]
    # a dissimilarity that is NaN or infinite (the lambda cells), or a
    # uniform target a design cannot attain, stops in the engine's checks, as
    # the leaf expects from its inputs
    res$blp <- golden_blp("nl_blp_contraction", list(
      rep(0, J), target, d$X, beta, lam, d$alt_idx, nest_idx, d$M, d$weights,
      ioo, max_iter = 25L))
    res$blp_damped <- golden_blp("nl_blp_contraction", list(
      rep(0, J), target, d$X, beta, lam, d$alt_idx, nest_idx, d$M, d$weights,
      ioo, damping = 0.5, max_iter = 25L))
  }
  for (lb in names(labels)) {
    res[[paste0("meat_", lb)]] <- p2_val(
      p8_meat_cluster("nl", theta, d, asc, ioo, labels[[lb]], nest_idx))
  }
  res
}

P8B_MNL_DESIGNS <- list(
  k17 = list(N = 60L, J = 40L, lo = 10L, hi = 30L, K = 17L, runs = TRUE,
             unsorted = FALSE, dup = FALSE),
  k4sq = list(N = 80L, J = 4L, lo = 4L, hi = 4L, K = 4L, runs = TRUE,
              unsorted = FALSE, dup = FALSE),
  k4sq_unsorted = list(N = 80L, J = 4L, lo = 4L, hi = 4L, K = 4L, runs = TRUE,
                       unsorted = TRUE, dup = FALSE),
  k1 = list(N = 60L, J = 20L, lo = 3L, hi = 20L, K = 1L, runs = FALSE,
            unsorted = TRUE, dup = FALSE),
  k3dup = list(N = 60L, J = 25L, lo = 2L, hi = 8L, K = 3L, runs = TRUE,
               unsorted = TRUE, dup = TRUE))

run_p8b_mnl <- function() {
  out <- list()
  i <- 0L
  for (nm in names(P8B_MNL_DESIGNS)) {
    g <- P8B_MNL_DESIGNS[[nm]]
    for (ioo in c(FALSE, TRUE)) {
      for (asc in c(TRUE, FALSE)) {
        i <- i + 1L
        d <- p8b_design(16000L + i, g$N, g$J, g$lo, g$hi, g$K, ioo,
                        runs = g$runs, unsorted = g$unsorted, dup = g$dup,
                        sd_x = 0.5)
        theta <- p8_theta(16100L + i, g$K, n_delta_free(g$J, asc, ioo))
        key <- sprintf("mnl|%s|oo=%d|asc=%d", nm, ioo, asc)
        out[[key]] <- p8b_mnl_cells(theta, d, asc, ioo,
                                    p8_labels(16200L + i, d$N)["coarse"], g$K)
      }
    }
  }
  out
}

# Nest structures over J alternatives (codes 1..J): `nests` holds each
# alternative's nest.
P8B_NL_DESIGNS <- list(
  n60 = list(N = 80L, J = 120L, lo = 3L, hi = 8L, K = 3L, runs = TRUE,
             unsorted = TRUE, nests = rep(1:60, each = 2L)),
  odd7 = list(N = 60L, J = 21L, lo = 3L, hi = 8L, K = 2L, runs = FALSE,
              unsorted = FALSE, nests = rep(1:7, each = 3L)),
  even8 = list(N = 60L, J = 24L, lo = 3L, hi = 8L, K = 2L, runs = TRUE,
               unsorted = FALSE, nests = rep(1:8, each = 3L)),
  singles = list(N = 60L, J = 20L, lo = 2L, hi = 7L, K = 2L, runs = TRUE,
                 unsorted = TRUE, nests = c(rep(1:6, each = 3L), 7L, 8L)),
  k17 = list(N = 50L, J = 40L, lo = 10L, hi = 30L, K = 17L, runs = TRUE,
             unsorted = FALSE, nests = rep(1:5, each = 8L)),
  k4sq = list(N = 80L, J = 4L, lo = 4L, hi = 4L, K = 4L, runs = TRUE,
              unsorted = TRUE, nests = c(1L, 1L, 2L, 2L)))

run_p8b_nl <- function() {
  out <- list()
  i <- 0L
  for (nm in names(P8B_NL_DESIGNS)) {
    g <- P8B_NL_DESIGNS[[nm]]
    nest_idx <- as.integer(g$nests)
    for (ioo in c(FALSE, TRUE)) {
      for (asc in c(TRUE, FALSE)) {
        i <- i + 1L
        d <- p8b_design(17000L + i, g$N, g$J, g$lo, g$hi, g$K, ioo,
                        runs = g$runs, unsorted = g$unsorted, sd_x = 0.5)
        theta <- p8_nl_theta(17100L + i, g$K, nest_idx,
                             n_delta_free(g$J, asc, ioo))
        key <- sprintf("nl|%s|oo=%d|asc=%d", nm, ioo, asc)
        out[[key]] <- p8b_nl_cells(theta, d, nest_idx, asc, ioo,
                                   p8_labels(17200L + i, d$N)["coarse"], g$K)
      }
    }
  }
  out
}

# Non-finite and boundary cells.
run_p8b_edge <- function() {
  out <- list()
  # A zero weight on a situation whose X' d overflows, on BLAS's path (sets of
  # 6, K = 3: alpha = 0 returns before the product) and on Armadillo's tiny
  # square path (sets of 4, K = 4: 0 * Inf); the column's coefficient keeps
  # the utilities finite.
  for (shape in c("blas", "tinysq")) {
    K <- if (shape == "blas") 3L else 4L
    J <- if (shape == "blas") 8L else 4L
    m <- if (shape == "blas") 6L else 4L
    for (model in c("mnl", "nl")) {
      for (ioo in c(FALSE, TRUE)) {
        seed <- 18000L + 10L * (shape == "tinysq") + 2L * (model == "nl") + ioo
        d <- p8b_design(seed, 30L, J, m, m, K, ioo)
        t <- 4L
        rows <- c(0L, cumsum(d$M))[t] + seq_len(d$M[t])
        d$X[rows, K] <- 1e308
        d$weights[t] <- 0
        if (model == "mnl") {
          theta <- p8_theta(seed + 100L, K, n_delta_free(J, TRUE, ioo))
          theta[K] <- 1e-300
          out[[sprintf("w0overflow|%s|mnl|oo=%d", shape, ioo)]] <-
            p8b_mnl_cells(theta, d, TRUE, ioo, list(), K)
        } else {
          nest_idx <- as.integer(rep(seq_len(J / 2L), each = 2L))
          theta <- p8_nl_theta(seed + 100L, K, nest_idx, n_delta_free(J, TRUE, ioo))
          theta[K] <- 1e-300
          out[[sprintf("w0overflow|%s|nl|oo=%d", shape, ioo)]] <-
            p8b_nl_cells(theta, d, nest_idx, TRUE, ioo, list(), K)
        }
      }
    }
  }
  # Infinite utilities inside a run of equal choice sets (not its last
  # situation), and an NA covariate, through every kernel.
  nest_idx <- as.integer(rep(1:6, each = 2L))
  for (model in c("mnl", "nl")) {
    for (ioo in c(FALSE, TRUE)) {
      seed <- 18100L + 2L * (model == "nl") + ioo
      d0 <- p8b_design(seed, 40L, 12L, 4L, 4L, 2L, ioo)
      # situations 2-6 share situation 2's set
      for (t in 3:6) {
        r0 <- c(0L, cumsum(d0$M))
        d0$alt_idx[r0[t] + seq_len(d0$M[t])] <- d0$alt_idx[r0[2L] + seq_len(d0$M[2L])]
      }
      r0 <- c(0L, cumsum(d0$M))
      th0 <- if (model == "mnl") {
        p8_theta(seed + 100L, 2L, n_delta_free(12L, TRUE, ioo))
      } else {
        p8_nl_theta(seed + 100L, 2L, nest_idx, n_delta_free(12L, TRUE, ioo))
      }
      cells <- list(
        plus_inf = function(d, th) {
          d$X[r0[4L] + 2L, 1L] <- 1e308; th[1L] <- 3; list(d, th)
        },
        minus_inf = function(d, th) {
          d$X[r0[4L] + 1L, 1L] <- -1e308; th[1L] <- 3; list(d, th)
        },
        na_x = function(d, th) { d$X[r0[4L] + 2L, 2L] <- NA_real_; list(d, th) })
      for (nm in names(cells)) {
        e <- cells[[nm]](d0, th0)
        key <- sprintf("run|%s|%s|oo=%d", model, nm, ioo)
        out[[key]] <- if (model == "mnl") {
          p8b_mnl_cells(e[[2L]], e[[1L]], TRUE, ioo, list(), 2L)
        } else {
          p8b_nl_cells(e[[2L]], e[[1L]], nest_idx, TRUE, ioo, list(), 2L)
        }
      }
    }
  }
  # Dissimilarities at the edges: near 1e-16 (the fits' lower bound), +Inf and
  # NaN, on 60 nests of 2.
  g <- P8B_NL_DESIGNS$n60
  nest_idx <- as.integer(g$nests)
  for (ioo in c(FALSE, TRUE)) {
    d <- p8b_design(18200L + ioo, g$N, g$J, g$lo, g$hi, g$K, ioo, runs = TRUE,
                    sd_x = 0.5)
    th0 <- p8_nl_theta(18210L + ioo, g$K, nest_idx, n_delta_free(g$J, TRUE, ioo))
    lam_vals <- list(tiny = 1e-16, inf = Inf, nan = NaN)
    for (nm in names(lam_vals)) {
      th <- th0
      th[g$K + c(1L, 37L)] <- lam_vals[[nm]]
      out[[sprintf("lambda|%s|oo=%d", nm, ioo)]] <-
        p8b_nl_cells(th, d, nest_idx, TRUE, ioo, list(), g$K)
    }
  }
  # nl_blp_contraction() takes lambda as given: zero and negative entries.
  d <- p8b_design(18300L, 40L, 12L, 3L, 6L, 2L, FALSE)
  nest_idx <- as.integer(rep(1:6, each = 2L))
  lam_vals <- list(zero = 0, negative = -0.5)
  for (nm in names(lam_vals)) {
    lam <- rep(0.7, 6L)
    lam[c(2L, 5L)] <- lam_vals[[nm]]
    out[[sprintf("blp_lambda|%s", nm)]] <- golden_blp("nl_blp_contraction", list(
      rep(0, 12L), rep(1 / 12, 12L), d$X, c(0.3, -0.2), lam, d$alt_idx,
      nest_idx, d$M, d$weights, FALSE, max_iter = 10L), ref = FALSE)
  }
  out
}

# The MXL elasticity and diversion kernels on wide sets (straddling 16, with
# repeated alternatives), +- outside option, store and generate draws.
run_p8b_mxl <- function() {
  out <- list()
  i <- 0L
  for (ioo in c(FALSE, TRUE)) {
    for (gen in c(FALSE, TRUE)) {
      for (dup in c(FALSE, TRUE)) {
        i <- i + 1L
        set.seed(19000L + i)
        M <- sample(10:30, 40L, replace = TRUE)
        f <- p3_case(19100L + i, M, J = 40L, K_x = 3L, K_w = 2L, S = 7L,
                     dup = dup, ioo = ioo, rc_corr = TRUE, rc_mean = TRUE,
                     rc_dist = c(0L, 1L), gen = gen)
        f$var_r <- 2L
        key <- sprintf("mxl|oo=%d|gen=%d|dup=%d", ioo, gen, dup)
        out[[key]] <- lapply(stats::setNames(nm = c("elas_fixed", "elas_random",
                                                    "dr_fixed", "dr_random")),
                             function(k) p2_val(p3_call(k, f)))
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# p8b_cl: cluster labels of every type through vcov(type = "cluster")
# ---------------------------------------------------------------------------

# The cluster labels of a frame's situations (ids 1..N), one per situation,
# named by id: integer clusters of 15 situations, in other types.
p8b_label_sets <- function(n_sit, seed) {
  set.seed(seed)
  g <- (seq_len(n_sit) - 1L) %/% 15L
  G <- max(g) + 1L
  sets <- list(
    integer = g,
    logical = g %% 2L == 0L,
    character = sprintf("c%03d", g),
    factor = factor(sprintf("c%03d", g)),
    factor_rev = factor(sprintf("c%03d", g), levels = rev(sprintf("c%03d", 0:(G - 1L)))),
    factor_unused = factor(sprintf("c%03d", g), levels = c("zz", sprintf("c%03d", 0:(G - 1L)), "yy")),
    factor_na = addNA(factor(sprintf("c%03d", g))),
    factor_dup = structure(as.integer(g %% 4L) + 1L, levels = c("a", "a", "b", "c"),
                           class = "factor"),
    double = c(0.1 + 0.2, 0.3, 1.5, 2.25)[g %% 4L + 1L],
    date = as.Date("2020-01-01") + g)
  lapply(sets, function(x) stats::setNames(x, as.character(seq_len(n_sit))))
}

p8b_cl_record <- function(fit, sets) {
  lapply(sets, function(cl) p2_val(vcov(fit, type = "cluster", cluster = cl)))
}

run_p8b_cl <- function() {
  out <- list()
  sets <- NULL
  for (model in c("mnl", "nl")) {
    nests <- if (model == "nl") as.integer(c(rep(1:8, each = 3L), 9:14))
    dt <- p8_frame(if (model == "nl") 15100L else 15000L, 300L, 30L, 4L, 8L,
                   FALSE, nests = nests)
    sets <- p8b_label_sets(300L, 19500L)
    fit <- suppressWarnings(if (model == "nl") {
      run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
    } else {
      run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
    })
    out[[paste0("posthoc|", model)]] <- p8b_cl_record(fit, sets)
    # Fit-time cluster_col of several types
    for (ty in c("integer", "character", "factor", "double")) {
      dtc <- data.table::copy(dt)
      v <- sets[[ty]][dtc$id]
      dtc[, clus := v]
      out[[sprintf("fittime|%s|%s", model, ty)]] <- p8_fit_record(
        if (model == "nl") {
          run_nestlogit(dtc, "id", "alt", "choice", c("x1", "x2"), "nest",
                        se_method = "cluster", cluster_col = "clus")
        } else {
          run_mnlogit(dtc, "id", "alt", "choice", c("x1", "x2"),
                      se_method = "cluster", cluster_col = "clus")
        })
    }
  }
  # Mixed logit: a cross-section, and a panel whose clusters nest persons.
  dt <- p8_frame(15200L, 240L, 12L, 3L, 6L, FALSE)
  dt[, person := (id - 1L) %/% 3L + 1L]
  for (panel in c(FALSE, TRUE)) {
    fit <- suppressWarnings(run_mxlogit(
      dt, "id", "alt", "choice", "x1", "x2", S = 20L, draws = "generate",
      seed = 7L, person_col = if (panel) "person"))
    s <- p8b_label_sets(240L, 19600L)
    if (panel) {
      # constant within each person (persons are 3 consecutive situations)
      s <- lapply(s, function(x) {
        y <- x[rep(seq(1L, 240L, by = 3L), each = 3L)]
        stats::setNames(y, names(x))
      })
    }
    out[[sprintf("posthoc|mxl|panel=%d", panel)]] <- p8b_cl_record(fit, s)
  }
  out
}
