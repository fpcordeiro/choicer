# ============================================================================
# D2 mode_rp: revealed-preference intercity mode choice, a population built
# from the shipped `mode_choice` (Greene & Hensher 1997; 210 travellers, all
# of air, bus, car and train offered). That sample is choice-based, so the
# population resamples traveller type i with probability
#   pi_i = Q_{y_i} / n_{y_i}            (its WESML weight Q/H over the sample size),
# Q the population mode shares of Greene's treatment (air .14, train .13, bus
# .09, car .64) and n_j the sample's choosers of j: the population's joint
# covariate distribution. Utility, air the reference (the smallest label):
#   V_ij = b_w wait_ij + b_t travel_ij + b_c vcost_ij + b_ia income_i 1{j = air} + delta_j.
# Variants (theta0 frozen below):
#   "default"  MNL; theta0 = the WESML MNL fit of mode_choice at Q. With a
#              full set of constants the WESML first-order conditions make the
#              population's expected shares at theta0 equal Q exactly.
#   "nl"       NL (RU2), nests ground = {bus, car, train}, fly = {air} (a
#              singleton, lambda 1); theta0 = the WESML NL fit.
#   "mxl"      MXL, a normal coefficient on travel (Mu_travel, L_11 = log sd).
#              The free WESML MXL fit puts sigma at 0 (sigma / SE(sigma) =
#              0.05), so sigma = |Mu_travel| (coefficient of variation 1): the
#              WESML MXL fit with sigma held at |Mu_travel| (a fixed point),
#              then the constants re-solved so that the population's expected
#              shares equal Q (the MXL's first-order conditions do not match
#              shares). theta_star (D2_THETA_STAR) is the MNL's pseudo-true
#              value under it (F11's truth), exact.
# The MXL's probabilities integrate the travel coefficient by Gauss-Hermite
# with D2_GH_NODES nodes: sigma * travel spans about 15 utils per standard
# deviation, so the logit is steep in the draw and the rule needs about 2000
# nodes for 1e-12 (1000 give 2e-9).
# Registers DGP "D2" at source time; definitions otherwise. d2_calibrate()
# reproduces the frozen values (it fits choicer; sourcing never runs it).
# ============================================================================

D2_VARIANTS <- c("default", "nl", "mxl")
D2_MODEL <- c(default = "mnl", nl = "nl", mxl = "mxl")
D2_ALTS <- c("air", "bus", "car", "train")                      # choicer's codes 1..4
D2_Q <- c(air = 0.14, bus = 0.09, car = 0.64, train = 0.13)    # Greene's population shares
D2_X <- c("wait", "travel", "vcost", "income_air")
D2_NEST <- c(air = "fly", bus = "ground", car = "ground", train = "ground")
D2_GH_NODES <- 2000L                   # Gauss-Hermite nodes of the travel coefficient's integral

# The suite's specification of a variant (lib/simulate.R).
d2_spec <- function(variant = "default") {
  variant <- match.arg(variant, D2_VARIANTS)
  switch(variant,
         default = list(beta = D2_X, alts = D2_ALTS),
         nl = list(beta = D2_X, alts = D2_ALTS, nests = unname(D2_NEST[D2_ALTS])),
         mxl = list(beta = setdiff(D2_X, "travel"), alts = D2_ALTS, random = "travel",
                    rc_dist = 0L, rc_mean = TRUE, correlated = FALSE))
}

# theta0 per variant (choicer's names and order). Calibrated by
# d2_calibrate() on 2026-11-02 with choicer 0.2.1.9000 at one thread: WESML
# fits of mode_choice at D2_Q, polished by Newton (decrements below 1e-13);
# "mxl" by generate-mode draws (S = 5000, seed 1), sigma held at |Mu_travel|
# (secant fixed point, 9 fits), then the constants matched to Q.
D2_THETA0 <- list(
  default = c(wait = -0.13231547689285264, travel = -0.0031859149548581498,
              vcost = -0.011245984858241102, income_air = 0.00094640773966060091,
              ASC_bus = -2.2302828402765487, ASC_car = -5.578990684385257,
              ASC_train = -1.981287797982572),
  nl = c(wait = -0.086759510523331085, travel = -0.0029301066098455049,
         vcost = -0.0096616325615906402, income_air = 0.0028399069475784234,
         Lambda_1 = 0.54501549923023007, ASC_bus = -0.29872311865278933,
         ASC_car = -2.7225232345814097, ASC_train = -0.11500412268264219),
  mxl = c(wait = -0.17142211847528663, vcost = -0.032668227218943094,
          income_air = 0.030013278155533726, Mu_travel = -0.011097269455489046,
          L_11 = -4.5010561917560485, ASC_bus = 4.0815189109041725,
          ASC_car = -0.91502768607287199, ASC_train = 4.6552996700932505))

# The MNL's pseudo-true value under the "mxl" population (F11's truth):
# argmax_theta sum_i pi_i sum_j P_true(j | i) log P_MNL(j | i; theta), exact
# up to the fit's convergence (d2_theta_star_fit(), same date and build).
D2_THETA_STAR <- c(wait = -0.09400971537680089, travel = -0.006341521292787246,
                   vcost = -0.019261853546803032, income_air = 0.0043542798622291179,
                   ASC_bus = 0.25813988130947924, ASC_car = -2.0247088214200413,
                   ASC_train = 0.63443814324373859)

.d2_cache <- new.env(parent = emptyenv())

# --- The population ----------------------------------------------------------

# The 210 traveller types: `long` (4 rows per type in choicer's order: type,
# alt, wait, travel, vcost, income, income_air, nest), `y` (each type's mode
# in mode_choice), `n_y` (the sample's choosers per mode), `pi` (population
# weights, summing to 1) and `n`. Cached.
d2_types <- function() {
  if (!is.null(.d2_cache$types)) return(.d2_cache$types)
  long <- d2_mode_data()
  long[, type := match(id, unique(id))]
  n <- max(long$type)
  if (nrow(long) != 4L * n || !identical(long$alt, rep(D2_ALTS, n)) ||
      !all(tabulate(long$type[long$choice == 1L], n) == 1L)) {
    stop("mode_choice is not 4 modes per traveller with one choice each")
  }
  y <- long$alt[long$choice == 1L]
  n_y <- table(factor(y, levels = D2_ALTS))
  pi <- unname(D2_Q[y] / as.numeric(n_y[y]))
  long <- long[, list(type, alt, wait, travel, vcost, income, income_air, nest)]
  .d2_cache$types <- list(long = long, y = y, n_y = n_y, pi = pi, n = n)
}

# Gauss-Hermite rule for E f(z), z ~ N(0, 1) (Golub-Welsch): nodes z,
# weights w. Cached by n.
d2_gh <- function(n = D2_GH_NODES) {
  key <- paste0("gh_", n)
  if (!is.null(.d2_cache[[key]])) return(.d2_cache[[key]])
  jac <- matrix(0, n, n)
  jac[cbind(1:(n - 1L), 2:n)] <- jac[cbind(2:n, 1:(n - 1L))] <- sqrt(1:(n - 1L) / 2)
  eg <- eigen(jac, symmetric = TRUE)
  .d2_cache[[key]] <- list(z = sqrt(2) * eg$values, w = eg$vectors[1L, ]^2)
}

# Each type's choice probabilities at theta (a types x 4 matrix, columns
# D2_ALTS): MNL and NL by the engines, MXL integrated over the travel
# coefficient by Gauss-Hermite.
d2_type_probs <- function(variant = "default", theta = D2_THETA0[[variant]],
                          nodes = D2_GH_NODES) {
  ty <- d2_types()
  model <- D2_MODEL[[variant]]
  spec <- d2_spec(variant)
  par <- unpack_theta(model, spec, theta)
  if (model != "mxl") {
    d <- design(ty$long, spec, id = "type")
    p <- if (model == "mnl") mnl_probs(d, par)$p else nl_probs(d, par)$p
    return(matrix(p, ty$n, 4L, byrow = TRUE, dimnames = list(NULL, D2_ALTS)))
  }
  # The types replicated once per node (node-major situations), cached by n.
  gh <- d2_gh(nodes)
  key <- paste0("gh_design_", nodes)
  if (is.null(.d2_cache[[key]])) {
    rep_long <- ty$long[rep(seq_len(nrow(ty$long)), nodes)]
    rep_long[, rep_id := rep(seq_len(ty$n * nodes), each = 4L)]
    .d2_cache[[key]] <- design(rep_long, spec, id = "rep_id")
  }
  b <- rc_coefs(par$mu, par$L, matrix(gh$z, 1L), par$rc_dist, par$rc_mean)
  p <- mxl_cond_probs(.d2_cache[[key]], par, matrix(rep(b, each = ty$n), ncol = 1L))$p
  p <- matrix(p, 4L * ty$n, nodes) %*% gh$w
  matrix(p, ty$n, 4L, byrow = TRUE, dimnames = list(NULL, D2_ALTS))
}

# The population's expected shares at theta: sum_i pi_i P_i. Cached for the
# frozen theta0.
d2_pop_shares <- function(variant = "default", theta = NULL, nodes = D2_GH_NODES) {
  if (is.null(theta)) {
    key <- paste0("shares_", variant)
    if (is.null(.d2_cache[[key]])) {
      .d2_cache[[key]] <- drop(d2_types()$pi %*% d2_type_probs(variant))
    }
    return(.d2_cache[[key]])
  }
  drop(d2_types()$pi %*% d2_type_probs(variant, theta, nodes))
}

# --- Simulation ----------------------------------------------------------------

# N situations of the variant: types drawn from pi, then (MXL) each
# traveller's travel coefficient, then the choices by inverse CDF, all from
# `seed` in that order. data: id, alt, choice, wait, travel, vcost, income,
# income_air, nest, type (choicer's row order, no outside rows). meta: the
# variant, model and spec, N, Q, the population's expected shares at theta0
# (`pop_shares`), the probabilities the choices were drawn from (`p`, per
# row; MXL: given the tastes), `type`, MXL `tastes` (N x 1) and theta_star.
d2_simulate <- function(N, seed, variant = "default") {
  variant <- match.arg(variant, D2_VARIANTS)
  ty <- d2_types()
  model <- D2_MODEL[[variant]]
  spec <- d2_spec(variant)
  theta <- D2_THETA0[[variant]]
  par <- unpack_theta(model, spec, theta)
  N <- as.integer(N)
  tastes <- NULL
  vs_with_seed(seed, {
    type <- sample.int(ty$n, N, replace = TRUE, prob = ty$pi)
    rows <- (rep(type, each = 4L) - 1L) * 4L + rep(1:4, N)  # outside [ ]: `type` is a column
    data <- ty$long[rows]
    data[, id := rep(seq_len(N), each = 4L)]
    d <- design(data, spec)
    p <- switch(model, mnl = mnl_probs(d, par), nl = nl_probs(d, par), mxl = {
      tastes <- draw_tastes(N, par$mu, par$L, par$rc_dist, par$rc_mean)
      colnames(tastes) <- spec$random
      mxl_cond_probs(d, par, tastes)
    })
    data[, choice := draw_choice(p$p, NULL, d$sit)$choice]
  })
  data.table::setcolorder(data, c("id", "alt", "choice", "wait", "travel", "vcost", "income",
                                  "income_air", "nest", "type"))
  list(data = data, truth = theta,
       meta = list(variant = variant, model = model, spec = spec, N = N, Q = D2_Q,
                   pop_shares = d2_pop_shares(variant), p = p$p, type = type, tastes = tastes,
                   theta_star = if (variant == "mxl") D2_THETA_STAR))
}

# Calibration moments: the simulated mode shares; the population's exact
# expected shares at theta0 and their largest gap from Q (`pop_gap`: 0 up
# to the WESML fit's convergence for "default"); the mean VTTS
# (b_travel / b_vcost, per hour; MXL: over the travellers' tastes); the air
# share by income tercile.
d2_moments <- function(sim) {
  dt <- sim$data
  ch <- dt[choice == 1L]
  sh <- table(factor(ch$alt, levels = D2_ALTS)) / nrow(ch)
  th <- sim$truth
  b_t <- if (sim$meta$model == "mxl") sim$meta$tastes[, "travel"] else th[["travel"]]
  inc <- dt$income[dt$alt == "air"]
  terc <- findInterval(inc, stats::quantile(inc, c(1, 2) / 3), left.open = TRUE) + 1L
  air <- dt$choice[dt$alt == "air"]
  gap <- max(abs(sim$meta$pop_shares - D2_Q[D2_ALTS]))
  out <- c(stats::setNames(as.numeric(sh), paste0("share_", D2_ALTS)),
           stats::setNames(sim$meta$pop_shares, paste0("pop_share_", D2_ALTS)),
           stats::setNames(gap, if (sim$meta$variant == "default") "pop_gap"
                                else paste0("pop_gap_", sim$meta$variant)),
           vtts_per_hour = 60 * mean(b_t) / th[["vcost"]],
           stats::setNames(tapply(air, factor(terc, 1:3), mean), paste0("air_share_inc_t", 1:3)))
  out
}

D2_TARGETS <- data.frame(
  moment = c(paste0("share_", D2_ALTS), "pop_gap"),
  lo = c(D2_Q[D2_ALTS] - 0.01, 0), hi = c(D2_Q[D2_ALTS] + 0.01, 1e-8),
  stringsAsFactors = FALSE)

# --- Calibration (choicer fits; never run at source time) ------------------------

# The shipped mode_choice in the DGP's columns (id, alt, choice, wait,
# travel, vcost, income, income_air, nest), choicer's row order; the
# references (lib/references.R) fit it too.
d2_mode_data <- function() {
  e <- new.env()
  utils::data("mode_choice", package = "choicer", envir = e)
  mc <- data.table::as.data.table(e$mode_choice)
  dt <- data.table::data.table(id = mc$id, alt = as.character(mc$mode), choice = mc$choice,
                               wait = mc$wait, travel = mc$travel, vcost = mc$vcost,
                               income = mc$income)
  dt[, income_air := income * (alt == "air")]
  dt[, nest := unname(D2_NEST[alt])]
  data.table::setorder(dt, id, alt)
  dt
}

# mode_choice with WESML weights at D2_Q attached (`.wesml_weight`).
d2_wesml_data <- function() {
  choicer::wesml_weights(d2_mode_data(), "id", "alt", "choice", Q = D2_Q, attach = TRUE)
}

# Weighted fits warn that weights may be sampling weights; calibration knows.
.d2_quiet <- function(expr) {
  withCallingHandlers(suppressMessages(expr), warning = function(w) {
    if (startsWith(conditionMessage(w), "Non-uniform weights detected")) {
      invokeRestart("muffleWarning")
    }
  })
}

# A weighted MNL, or NL with the ground nest (D2's specification, polished
# by Newton; the NL in the "bhhh" coordinates), on the 210 types with one
# situation per (type i, mode j), choice j, weight W[i, j] (a types x 4
# matrix; cells of weight 0 left out): the log-likelihood of a population
# (W = pi_i P(j | i); its Hessian and BHHH routes are then the population's
# H^-1 and B^-1 per traveller) or of a sample (W = its (type, mode) counts:
# exactly the log-likelihood of its travellers). `...` goes to the fit.
d2_pseudo_fit <- function(W, model = c("mnl", "nl"), ...) {
  model <- match.arg(model)
  ty <- d2_types()
  pseudo <- ty$long[rep(seq_len(nrow(ty$long)), 4L)]
  k <- rep(1:4, each = nrow(ty$long))                      # the pseudo-situation's choice
  pseudo[, `:=`(id = (k - 1L) * ty$n + type, choice = as.integer(alt == D2_ALTS[k]))]
  pseudo[, w := W[cbind(type, k)]]
  pseudo <- pseudo[w > 0]
  data.table::setorder(pseudo, id, alt)
  .d2_quiet(if (model == "mnl") {
    choicer::run_mnlogit(pseudo, "id", "alt", "choice", D2_X, weights_col = "w",
                         polish = "newton", ...)
  } else {
    choicer::run_nestlogit(pseudo, "id", "alt", "choice", D2_X, nest_col = "nest",
                           weights_col = "w", scale_vars = "bhhh", polish = "newton", ...)
  })
}

# The (type, mode) distribution of a Bernoulli choice-based sample from the
# variant's population, design shares H against the population's shares Q:
# W[i, j] = pi_i P_ij f_j / sum, f_j = H_j / Q_j (the sum is 1).
d2_cb_weights <- function(variant = "default", H = rep(0.25, 4L)) {
  P <- d2_type_probs(variant)
  W <- d2_types()$pi * P
  sweep(W, 2L, H / colSums(W), `*`)
}

# The MNL's pseudo-true value under the "mxl" truth, exactly: the
# population is a mixture over the 210 types, so theta_star maximizes
# sum_i pi_i sum_j P_true(j | i) log P_MNL(j | i; theta), P_true by
# Gauss-Hermite: d2_pseudo_fit() at W = pi_i P_true(j | i). Returns the fit.
d2_theta_star_fit <- function(theta_mxl = D2_THETA0$mxl, nodes = D2_GH_NODES) {
  d2_pseudo_fit(d2_types()$pi * d2_type_probs("mxl", theta_mxl, nodes),
                control = list(xtol_rel = 1e-14))
}

# The "mxl" theta with its constants re-solved so that the population's
# expected shares equal Q: delta_j += log Q_j - log s_j (air's held at 0),
# until the log-share gaps are below `tol`.
d2_match_shares <- function(theta, nodes = D2_GH_NODES, tol = 1e-13, max_iter = 1000L) {
  asc <- paste0("ASC_", D2_ALTS[-1L])
  for (it in seq_len(max_iter)) {
    r <- log(D2_Q[D2_ALTS]) - log(d2_pop_shares("mxl", theta, nodes))
    if (max(abs(r - r[1L])) < tol) return(theta)
    theta[asc] <- theta[asc] + (r[-1L] - r[1L])
  }
  stop("d2_match_shares(): no convergence in ", max_iter, " iterations")
}

# Reproduces D2_THETA0 and D2_THETA_STAR: the WESML fits of mode_choice at
# D2_Q (MNL, NL; MXL with S draws, generate mode, seed 1), sigma's
# determination in the free MXL fit, and, when sigma / SE(sigma) < 2, the
# WESML MXL fit with sigma held at cv |Mu_travel| (the fixed point of sigma
# -> cv |Mu_travel(sigma)|, by secant steps), its constants then matched to
# Q (d2_match_shares()), and theta_star under it. Returns the values, the
# fixed point's path and the fits. About 3 minutes.
d2_calibrate <- function(S = 5000L, cv = 1, nodes = D2_GH_NODES) {
  mc <- d2_wesml_data()
  fit <- function(f, ...) .d2_quiet(f(mc, "id", "alt", "choice", weights_col = ".wesml_weight",
                                      se_method = "sandwich", polish = "newton", ...))
  mnl <- fit(choicer::run_mnlogit, covariate_cols = D2_X)
  nl <- fit(choicer::run_nestlogit, covariate_cols = D2_X, nest_col = "nest")
  sp <- d2_spec("mxl")
  mxl_fit <- function(...) {
    fit(choicer::run_mxlogit, covariate_cols = sp$beta, random_var_cols = sp$random,
        rc_dist = sp$rc_dist, rc_mean = TRUE, S = S, draws = "generate", seed = 1L, ...)
  }
  mxl_free <- mxl_fit()
  sigma_z <- 1 / mxl_free$se[["L_11"]]          # sigma / SE(sigma) by the delta method
  mxl <- mxl_free
  path <- NULL
  if (!is.finite(sigma_z) || sigma_z < 2) {
    held <- function(sig) {
      mxl <<- mxl_fit(theta_init = replace(coef(mxl), "L_11", log(sig)),
                      lower = c(L_11 = log(sig)), upper = c(L_11 = log(sig)))
      cv * abs(coef(mxl)[["Mu_travel"]]) - sig
    }
    s0 <- cv * abs(coef(mxl_free)[["Mu_travel"]])
    g0 <- held(s0)
    s1 <- s0 + g0
    for (it in 1:30) {
      g1 <- held(s1)
      path <- rbind(path, c(sigma = s1, gap = g1))
      if (abs(g1) < 1e-8 * s1) break
      s2 <- s1 - g1 * (s1 - s0) / (g1 - g0)
      s0 <- s1
      g0 <- g1
      s1 <- s2
    }
  }
  theta_mxl <- d2_match_shares(coef(mxl), nodes)
  star <- d2_theta_star_fit(theta_mxl, nodes)
  list(theta0 = list(default = coef(mnl), nl = coef(nl), mxl = theta_mxl),
       theta_star = coef(star), sigma_z = sigma_z, sigma_path = path,
       fits = list(mnl = mnl, nl = nl, mxl_free = mxl_free, mxl = mxl, star = star))
}

dgp_register(list(
  id = "D2", name = "mode_rp: intercity mode choice (mode_choice population, Greene's Q)",
  theta0 = D2_THETA0, simulate = d2_simulate, moments = d2_moments, targets = D2_TARGETS,
  variants = D2_VARIANTS))
