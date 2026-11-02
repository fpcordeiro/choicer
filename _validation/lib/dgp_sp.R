# ============================================================================
# D1 sp_panel: a stated-preference conjoint panel (Revelt & Train 1998), and
# D1-L (variant "labeled"), its labeled-products design for the hierarchical
# Bayes samplers (F9). Definitions only, except the registration of DGP "D1"
# (dgp_register()).
#
# Design (both variants). N persons x T tasks (T = 8). Every attribute takes
# discrete levels, balanced within each person and position: the T tasks of
# one person at one position carry a random permutation of the attribute's
# level set repeated to length T (each of 4 levels twice in 8 tasks, each of
# 2 levels four times), drawn independently per attribute, position and
# person. Every level appears equally often across a person's tasks and
# positions, and the attributes are independent of one another.
#
# D1 (variant "default"). Each task offers 3 unlabeled alternatives (codes
# 1, 2, 3) and an opt-out, the outside option (no rows, utility 0;
# include_outside_option = TRUE). Attributes and levels:
#   price    {1, 2, 3, 4} price units     fixed coefficient b_p
#   time     {0.5, 1, 1.5, 2} hours       enters as neg_time = -time
#   quality  {1, 2, 3, 4}, a rating
#   irr      {0, 1}, irrelevant           fixed coefficient 0 (F4's LR test)
# Utility of alternative j in task t for person i:
#   delta + b_p price + 0 irr + b_ti neg_time + b_qi quality + e_ijt,
#   b_ti = exp(mu_t) + exp(g_1i),  b_qi = mu_q + g_2i,  g_i = L eta_i,
# choicer's shifted log-normal and normal coefficients (rc_dist = c(1, 0),
# rc_mean = TRUE) with a correlated L (rc_correlation = TRUE), and one eta_i
# per person for all their tasks (the panel likelihood's model). The three
# position constants are equal in truth (delta: inside against the
# opt-out); choicer estimates three (ASC_1, ASC_2, ASC_3).
#
# theta0 (D1). Literals chosen 2026-11-02 (SP_TRUE): b_p = -1 per price
# unit; mu_q = 1.5, so that the mean WTP for one quality level is
# mu_q / |b_p| = 1.5 price units; sd(g_2) = 1, so sigma / |mu| = 0.67 for
# quality; exp(mu_t) = 0.5 and sd(g_1) = 0.7, so the time coefficient has
# mean 0.5 + exp(0.7^2 / 2) = 1.778 per hour, SD
# sqrt((exp(0.49) - 1) exp(0.49)) = 1.016 (SD / mean 0.57) and median 1.5
# (the median VTTS: 1.5 price units per hour). The underlying normals
# (g_1, g_2) correlate at rho = 0.4; the coefficients then correlate at
# rho sd(g_1) / sqrt(exp(sd(g_1)^2) - 1) = 0.35: people who mind time more
# also value quality more. L = chol(Sigma): L_11 = log 0.7, L_21 = 0.4,
# L_22 = log(sqrt(1 - 0.16)). delta was calibrated by sp_calibrate() the same
# day: the expected opt-out share 0.15 at N = 50,000 persons, seed 1.
#
# D1-L (variant "labeled"). J = 6 labeled products P1..P6, all offered in
# every task, with an implicit outside option (no rows). Each product has a
# fixed characteristic z_j (SPL_Z); the tasks vary price {1, 2, 3, 4} and
# feature {0, 1} (e.g. on display) by product. The hierarchical model of
# run_hmnlogit() / run_hmnprobit():
#   U_ijt = delta_j + price_ijt b_pi + feature_ijt b_fi + e_ijt,
#   (b_pi, b_fi) ~ N(b, W),  delta_j = theta_0 + theta_z z_j + xi_j,
#   xi_j ~ N(0, sigma_d^2),
# e Gumbel (shock = "logit", HMNL's model) or N(0, 1) on every product and
# the opt-out (shock = "probit", HMNP's model at sigma = 1, its identified
# scale). theta0 literals (2026-11-02, SPL_TRUE): b = (-0.8, 0.6), W with
# SDs 0.4 and 0.5 and correlation -0.3, theta_z = 0.6, sigma_d^2 = 0.25,
# xi = 0.5 times a fixed standard normal draw (SPL_XI_STD, rnorm(6) at seed
# 20261102), theta_0 calibrated by sp_calibrate("labeled") to the expected
# opt-out share 0.15 (N = 50,000, seed 1). J = 6 = 3P with one z (P = 2), the
# least that avoids run_hmnlogit()'s few-alternatives warning. Data in
# prepare_hmnl_data()'s layout: id_col "id", alt_col "alt", choice_col
# "choice", covariate_cols c("price", "feature"), person_col "person",
# alt_covariate_cols "z". The truth is named as choicer's draws: b by
# covariate, W[<row>,<col>] (vech, row by row), (Intercept) and z for theta,
# sigma_d2, and delta_<label>; meta$hb holds them as a list, with xi.
#
# Randomness, from `seed` in this order: the attributes (one T x J*N block
# of uniforms per attribute, in SP_LEVELS' order), the tastes' normals eta
# (2 x N), then the choices (one uniform per task by inverse CDF; probit:
# one normal per row, then one per task for the opt-out). Choice
# probabilities come from lib/simulate.R's engines: the logit given each
# person's tastes, mxl_cond_probs().
# ============================================================================

SP_VARIANTS <- c("default", "labeled")
SP_T <- 8L
SP_ALTS <- list(default = 1:3, labeled = sprintf("P%d", 1:6))
SP_X <- c("price", "irr")
SP_W <- c("neg_time", "quality")
SPL_X <- c("price", "feature")
SP_LEVELS <- list(default = list(price = c(1, 2, 3, 4), time = c(0.5, 1, 1.5, 2),
                                 quality = c(1, 2, 3, 4), irr = c(0, 1)),
                  labeled = list(price = c(1, 2, 3, 4), feature = c(0, 1)))

SP_TRUE <- list(
  beta = c(price = -1, irr = 0),
  mu = c(neg_time = log(0.5), quality = 1.5),
  s = c(neg_time = 0.7, quality = 1),          # sd of g = L eta
  rho = 0.4,                                   # corr(g_1, g_2)
  delta = 2.4059591418534101)                  # each position's constant (sp_calibrate())

SPL_Z <- c(P1 = -1, P2 = -0.6, P3 = -0.2, P4 = 0.2, P5 = 0.6, P6 = 1)
SPL_XI_STD <- c(P1 = -0.91804940550063097, P2 = 0.73987621094601408, P3 = 0.87236004768931952,
                P4 = -0.034328405463817449, P5 = 1.273575005071049, P6 = -1.0947943654440757)
SPL_TRUE <- list(
  b = c(price = -0.8, feature = 0.6),
  W = matrix(c(0.16, -0.06, -0.06, 0.25), 2L, dimnames = list(SPL_X, SPL_X)),
  theta = c("(Intercept)" = 1.2599756978744163, z = 0.6),   # theta_0: sp_calibrate("labeled")
  sigma_d2 = 0.25)

# The engines' specification (lib/simulate.R) of a variant: D1's MXL; D1-L's
# conditional logit given the tastes, written as an MXL without fixed
# coefficients (normal tastes with means b, L = chol(W), a constant per
# product).
sp_spec <- function(variant = "default") {
  variant <- match.arg(variant, SP_VARIANTS)
  if (variant == "default") {
    return(list(beta = SP_X, alts = SP_ALTS$default, outside = TRUE, random = SP_W,
                rc_dist = c(1L, 0L), rc_mean = TRUE, correlated = TRUE))
  }
  list(beta = character(0), alts = SP_ALTS$labeled, outside = TRUE, random = SPL_X,
       rc_dist = c(0L, 0L), rc_mean = TRUE, correlated = TRUE)
}

# L from the tastes' SDs and the correlation of their normals.
.sp_L <- function(s, rho) matrix(c(s[1], rho * s[2], 0, s[2] * sqrt(1 - rho^2)), 2L)

# A variant's truth: `theta`, the engines' (and for D1 choicer's) parameter
# vector; `truth`, the named truth families score against (D1: theta;
# D1-L: the hierarchical parameters named as choicer's draws); D1-L's `hb`
# list (b, W, theta, sigma_d2, xi, delta, Z). `override` replaces values:
# D1 a named numeric in choicer's names; D1-L a list with any of b, W,
# theta, sigma_d2, xi (xi defaults to sqrt(sigma_d2) * SPL_XI_STD).
sp_truth <- function(variant = "default", override = NULL) {
  variant <- match.arg(variant, SP_VARIANTS)
  spec <- sp_spec(variant)
  if (variant == "default") {
    p <- SP_TRUE
    th <- pack_theta("mxl", spec, list(
      beta = p$beta, mu = p$mu, L = .sp_L(unname(p$s[SP_W]), p$rho),
      asc = stats::setNames(rep(p$delta, 3L), SP_ALTS$default)))
    if (!is.null(override)) {
      bad <- setdiff(names(override), names(th))
      if (is.null(names(override)) || length(bad)) {
        stop("D1: `theta` must be named by choicer's parameters (",
             paste(names(th), collapse = ", "), ")", call. = FALSE)
      }
      th[names(override)] <- override
    }
    return(list(theta = th, truth = th))
  }
  if (!is.null(override) &&
      (!is.list(override) || is.null(names(override)) || !all(nzchar(names(override))) ||
       length(setdiff(names(override), c("b", "W", "theta", "sigma_d2", "xi"))))) {
    stop("D1-L: `theta` must be a named list with any of b, W, theta, sigma_d2, xi",
         call. = FALSE)
  }
  hb <- utils::modifyList(SPL_TRUE, as.list(override)[intersect(names(override), names(SPL_TRUE))])
  hb$b <- stats::setNames(as.numeric(hb$b), SPL_X)
  hb$W <- matrix(as.numeric(hb$W), 2L, dimnames = list(SPL_X, SPL_X))
  hb$theta <- stats::setNames(as.numeric(hb$theta), c("(Intercept)", "z"))
  hb$xi <- stats::setNames(as.numeric(override$xi %||% (sqrt(hb$sigma_d2) * SPL_XI_STD)),
                           SP_ALTS$labeled)
  hb$Z <- cbind("(Intercept)" = 1, z = SPL_Z)
  hb$delta <- stats::setNames(drop(hb$Z %*% hb$theta) + hb$xi, SP_ALTS$labeled)
  if (!isSymmetric(hb$W) || any(eigen(hb$W, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("D1-L: `W` must be symmetric positive definite", call. = FALSE)
  }
  th <- pack_theta("mxl", spec, list(mu = hb$b, L = t(chol(hb$W)), asc = hb$delta))
  i <- rep(1:2, 1:2)
  j <- sequence(1:2)
  truth <- c(hb$b, stats::setNames(hb$W[cbind(i, j)], sprintf("W[%s,%s]", SPL_X[i], SPL_X[j])),
             hb$theta, sigma_d2 = hb$sigma_d2,
             stats::setNames(hb$delta, paste0("delta_", SP_ALTS$labeled)))
  list(theta = th, truth = truth, hb = hb)
}

# --- Design ----------------------------------------------------------------------

# T x n levels: each column a random permutation of `levels` repeated to
# length T (one uniform per cell, columns in order); T must be a multiple of
# the number of levels, so that every column is balanced.
.sp_balanced <- function(levels, T, n) {
  if (T %% length(levels) != 0L) {
    stop("D1: T = ", T, " is not a multiple of ", length(levels),
         " levels, so the design cannot be balanced", call. = FALSE)
  }
  o <- order(rep(seq_len(n), each = T), stats::runif(T * n))
  matrix(rep_len(levels, T)[(o - 1L) %% T + 1L], T, n)
}

# One attribute for N persons x T tasks x J positions, in row order (person,
# task, position); columns of the balanced block person-major.
.sp_attribute <- function(levels, T, J, N) {
  as.vector(aperm(array(.sp_balanced(levels, T, J * N), c(T, J, N)), c(2L, 1L, 3L)))
}

# The long design without choices: person, id (task, unique over persons),
# alt, the attributes (and neg_time for D1, z for D1-L), in choicer's row
# order.
sp_design <- function(N, T, variant = "default") {
  alts <- SP_ALTS[[variant]]
  J <- length(alts)
  dt <- data.table::data.table(person = rep(seq_len(N), each = T * J),
                               id = rep(seq_len(N * T), each = J), alt = rep(alts, N * T))
  lev <- SP_LEVELS[[variant]]
  for (a in names(lev)) data.table::set(dt, j = a, value = .sp_attribute(lev[[a]], T, J, N))
  if (variant == "default") {
    dt[, neg_time := -time]
  } else {
    dt[, z := unname(SPL_Z[alt])]
  }
  dt
}

# Probit choices: the maximum of V + N(0, 1) over the task's rows and the
# opt-out (0 + N(0, 1)); one normal per row, then one per task.
.sp_probit_choice <- function(d, par, tastes) {
  N <- .n_sit(d$sit)
  u <- .utility(d, par) + .random_part(d$W, t(tastes), .unit_of_row(d)) +
    stats::rnorm(length(d$sit))
  u0 <- stats::rnorm(N)
  best <- u == .group_max(u, d$sit, N)[d$sit] & u > u0[d$sit]
  chosen <- integer(N)
  chosen[d$sit[best]] <- which(best)
  list(chosen = chosen, choice = as.integer(best))
}

# --- Simulation ------------------------------------------------------------------------

# A dataset: list(data, truth, meta). `data` (person, id, alt, choice, the
# attributes) in choicer's row order; `truth` sp_truth()'s; `meta`: variant,
# shock, spec, N, T, the engines' theta, D1-L's hb list, eta (2 x N), the
# persons' realized coefficients (`tastes`, N x 2; D1: b_t and b_q, D1-L:
# price and feature), the expected opt-out share given the tastes (logit),
# and with keep_probs the probabilities the choices were drawn from.
sp_simulate <- function(N, seed, variant = "default", T = SP_T, theta = NULL,
                        shock = c("logit", "probit"), keep_probs = FALSE) {
  variant <- match.arg(variant, SP_VARIANTS)
  shock <- match.arg(shock)
  if (variant == "default" && shock != "logit") {
    stop("D1: shock = \"probit\" is D1-L's", call. = FALSE)
  }
  N <- as.integer(N)
  T <- as.integer(T)
  tr <- sp_truth(variant, theta)
  spec <- sp_spec(variant)
  par <- unpack_theta("mxl", spec, tr$theta)
  pr <- NULL
  vs_with_seed(seed, {
    data <- sp_design(N, T, variant)
    eta <- matrix(stats::rnorm(2L * N), 2L)
    tastes <- t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean))
    colnames(tastes) <- spec$random
    d <- design(data, spec, person = "person")
    ch <- if (shock == "logit") {
      pr <- mxl_cond_probs(d, par, tastes)
      draw_choice(pr$p, pr$p0, d$sit)
    } else {
      .sp_probit_choice(d, par, tastes)
    }
  })
  data[, choice := ch$choice]
  data.table::setcolorder(data, c("person", "id", "alt", "choice"))
  list(data = data, truth = tr$truth,
       meta = list(variant = variant, shock = shock, spec = spec, N = N, T = T,
                   theta = tr$theta, hb = tr$hb, eta = eta, tastes = tastes,
                   outside_expected = if (!is.null(pr)) mean(pr$p0) else NA_real_,
                   probs = if (keep_probs) pr))
}

# --- Moments and calibration -----------------------------------------------------------

# D1's coefficient distribution at theta (analytic): mean, SD and SD / |mean|
# of b_t = exp(mu_t) + exp(g_1) and b_q = mu_q + g_2, their correlation, the
# mean WTP for quality (mu_q / |b_p|) and the mean and median VTTS (price
# units per hour).
sp_coef_moments <- function(theta = sp_truth("default")$theta) {
  par <- unpack_theta("mxl", sp_spec("default"), theta)
  Sig <- par$L %*% t(par$L)
  bp <- abs(par$beta[["price"]])
  m_t <- exp(par$mu[["neg_time"]]) + exp(Sig[1, 1] / 2)
  sd_t <- sqrt((exp(Sig[1, 1]) - 1) * exp(Sig[1, 1]))
  m_q <- par$mu[["quality"]]
  sd_q <- sqrt(Sig[2, 2])
  c(mean_b_time = m_t, sd_b_time = sd_t, cv_time = sd_t / m_t, mean_b_quality = m_q,
    sd_b_quality = sd_q, cv_quality = sd_q / abs(m_q),
    corr_b = Sig[2, 1] / (sqrt(Sig[2, 2]) * sqrt(exp(Sig[1, 1]) - 1)),
    wtp_quality = m_q / bp, vtts_mean = m_t / bp,
    vtts_median = (exp(par$mu[["neg_time"]]) + 1) / bp)
}

# Calibration moments of a dataset. Both variants: the realized and the
# expected opt-out share (targeted 10-20%). D1 from the persons' realized
# coefficients: mean, SD and SD / |mean| of b_t and b_q (targeted 0.5-1),
# their correlation, the mean WTP for quality (targeted about 1.5) and the
# mean VTTS; the positions' shares (equal in truth). D1-L: the products'
# shares, the tastes' means, SDs and correlation.
sp_moments <- function(sim) {
  m <- sim$meta
  dt <- sim$data
  n_sit <- m$N * m$T
  ch <- dt[choice == 1L]
  out <- c(outside_share = 1 - nrow(ch) / n_sit, outside_expected = m$outside_expected)
  tb <- m$tastes
  if (m$variant == "default") {
    bp <- abs(m$theta[["price"]])
    bt <- tb[, "neg_time"]
    bq <- tb[, "quality"]
    pos <- tabulate(ch$alt, 3L) / n_sit
    return(c(out, mean_b_time = mean(bt), sd_b_time = stats::sd(bt),
             cv_time = stats::sd(bt) / mean(bt), mean_b_quality = mean(bq),
             sd_b_quality = stats::sd(bq), cv_quality = stats::sd(bq) / abs(mean(bq)),
             corr_b = stats::cor(bt, bq), wtp_quality = mean(bq) / bp, vtts_mean = mean(bt) / bp,
             share_pos1 = pos[1], share_pos2 = pos[2], share_pos3 = pos[3]))
  }
  sh <- table(factor(ch$alt, levels = SP_ALTS$labeled)) / n_sit
  c(out, stats::setNames(as.numeric(sh), paste0("share_", names(sh))),
    mean_b_price = mean(tb[, "price"]), sd_b_price = stats::sd(tb[, "price"]),
    mean_b_feature = mean(tb[, "feature"]), sd_b_feature = stats::sd(tb[, "feature"]),
    corr_b = stats::cor(tb[, "price"], tb[, "feature"]))
}

SP_TARGETS <- data.frame(
  moment = c("outside_share", "cv_time", "cv_quality", "wtp_quality"),
  lo = c(0.10, 0.5, 0.5, 1.4),
  hi = c(0.20, 1.0, 1.0, 1.6),
  stringsAsFactors = FALSE)

# Calibrates the constants' level to an expected opt-out share `target`: D1's
# delta (every position's constant), D1-L's theta_0 (every product's delta),
# by the contraction level += logit(1 - target) - logit(1 - s0), s0 the mean
# opt-out probability over N persons' tasks given their tastes (one seed:
# common draws across iterations). Not run when the file is sourced; its
# result is frozen in SP_TRUE$delta and SPL_TRUE$theta.
sp_calibrate <- function(variant = "default", N = 50000L, seed = 1L, target = 0.15, T = SP_T,
                         tol = 1e-12, max_iter = 100L) {
  variant <- match.arg(variant, SP_VARIANTS)
  spec <- sp_spec(variant)
  th <- sp_truth(variant)$theta
  asc <- theta_blocks("mxl", spec)$asc
  vs_with_seed(seed, {
    data <- sp_design(as.integer(N), as.integer(T), variant)
    eta <- matrix(stats::rnorm(2L * N), 2L)
  })
  d <- design(data, spec, person = "person")
  par <- unpack_theta("mxl", spec, th)
  tastes <- t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean))
  shift <- 0
  for (it in seq_len(max_iter)) {
    par$delta[] <- unpack_theta("mxl", spec, replace(th, asc, th[asc] + shift))$delta
    s0 <- mean(mxl_cond_probs(d, par, tastes)$p0)
    step <- stats::qlogis(1 - target) - stats::qlogis(1 - s0)
    shift <- shift + step
    if (abs(step) < tol) break
  }
  level <- if (variant == "default") SP_TRUE$delta else SPL_TRUE$theta[["(Intercept)"]]
  list(level = level + shift, outside = s0, iterations = it)
}

# --- Registration ------------------------------------------------------------------------

dgp_register(list(
  id = "D1",
  name = paste("stated-preference conjoint panel: N persons x 8 tasks, 3 unlabeled",
               "alternatives + opt-out (D1-L: 6 labeled products)"),
  theta0 = list(default = sp_truth("default")$truth, labeled = sp_truth("labeled")$truth),
  simulate = function(N, seed, variant = "default", ...) sp_simulate(N, seed, variant, ...),
  moments = sp_moments, targets = SP_TARGETS, variants = SP_VARIANTS))
