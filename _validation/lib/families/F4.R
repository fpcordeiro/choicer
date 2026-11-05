# ============================================================================
# F4: mixed logit on a stated-preference panel (plan §5; claims C1, C2, C4)
# on D1 (lib/dgp_sp.R): N persons x 8 tasks, 3 unlabeled alternatives and an
# opt-out, a fixed price, an irrelevant attribute (true coefficient 0), and
# correlated random coefficients, shifted log-normal on -time and normal on
# quality (rc_mean = TRUE, rc_correlation = TRUE). Cells: persons in {200,
# 800, 3.2k}, the panel at S = 500 generate-mode draws (seeded from the
# rep's seed), the composite at S_xsec = 800 at 200 and 800 persons only
# (S1, 2026-11-02: sqrt(U)/S = 0.1 at 800 persons' 6,400 tasks); "200.q"
# (200 persons, S = 200 for both) for the quick tier.
# Methods (polish = "newton"; the main fit from choicer's default start and
# from theta0, the better certified one kept, fit_two_starts(), its record
# in diag$start; multi-start on reps 10, 20, ...; diag$plateau flags an L
# diagonal driven to -Inf):
#   panel  run_mxlogit(person_col = "person"): the panel likelihood (Revelt &
#          Train 1998), one draw block per person; routes hessian, bhhh,
#          robust; lr: the LR test of irr = 0, the fit without irr on the
#          same draws, started at the estimate without irr (df = 1).
#   xsec   the composite likelihood: the same model without person_col, one
#          draw block per task, cluster_col = "person"; routes hessian and
#          cluster; joint Wald by route (hessian, cluster).
# Each task's marginal is the mixed logit, so both are consistent for
# theta0; the composite is inefficient and its Hessian SEs are too small,
# its cluster SEs valid (Varin, Reid & Firth 2011).
# Derived, both: the mean WTP for quality, Mu_quality / -price (wtp()'s
# delta SE, truth wtp() at theta0), and Sigma_11, Sigma_21, Sigma_22, the
# covariance of g = L eta (sigma_vech(), delta method); their SEs come from
# the Hessian variance in both methods (metrics read a derived SE as the
# hessian route).
# Gates: panel at 3.2k (sqrt(U)/S = 0.11), bias, Hessian SE ratio and tails
# of every parameter, joint Wald size on the common parameters (beta, mu,
# L), LR size, Sigma's tails; xsec at 800, bias, cluster SE ratio and tails
# (gates$cell named by method). Findings: xsec Hessian SE ratio < 1; the
# efficiency loss ESE_xsec / ESE_panel at 200 and 800 persons, the panel at
# S = 500 against the composite at S = 800; rates over N.
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

f4_simulate <- function(cell, seed) {
  get_dgp("D1")$simulate(N = cell$N, seed = seed, variant = "default", T = cell$T)
}

# One D1 fit with S draws seeded from the rep's seed; `x` the fixed
# covariates.
f4_run <- function(sim, seed, S, x = SP_X, ...) {
  choicer::run_mxlogit(sim$data, "id", "alt", "choice", x, SP_W, rc_dist = c(1L, 0L),
                       rc_mean = TRUE, rc_correlation = TRUE, include_outside_option = TRUE,
                       S = S, draws = "generate", seed = seed + 1L, polish = "newton", ...)
}

# The same, certified (require_converged()).
f4_fit <- function(sim, seed, S, what, ...) require_converged(f4_run(sim, seed, S, ...), what)

# The mean WTP for quality and Sigma's vech (sigma_vech()), SEs from
# fit$vcov (the caller passes the composite's fit with its Hessian variance
# there: metrics read a derived SE as the hessian route).
f4_estimands <- function(fit, sim) {
  w <- choicer::wtp(fit, price_var = "price", attr_vars = "quality")
  w0 <- choicer::wtp(at_truth(fit, sim$truth), price_var = "price", attr_vars = "quality")
  sg <- names(coef(fit))[fit$param_map$sigma]
  s_hat <- sigma_vech(coef(fit)[sg], 2L, TRUE, if (!is.null(fit$vcov)) fit$vcov[sg, sg])
  s_hat$truth <- sigma_vech(sim$truth[sg], 2L, TRUE)$est
  rbind(data.frame(name = "wtp_quality", est = w$Estimate, se = w$Std_Error, truth = w0$Estimate,
                   stringsAsFactors = FALSE), s_hat)
}

# The common parameters (everything but the constants), by name.
.f4_common <- function(fit) names(coef(fit))[-as.integer(fit$param_map$asc)]

# A fit's diagnostics: fit_two_starts()'s record, the multi-start on reps
# 10, 20, ..., and the plateau: min_L_diag, the smaller of choicer's
# log-diagonal parameters L_11 and L_22 (truth -0.36 and -0.09), and plateau
# = min_L_diag < -5. L_11 to -Inf drives the time coefficient's SD sd(g_1)
# to 0; L_22 to -Inf makes Sigma singular, |corr(g_1, g_2)| -> 1 with
# sd(g_2) = |L_21|. The score vanishes along such a ridge, so the polish can
# certify a point there that is not the maximum.
.f4_diag <- function(fit, refit, cell, seed) {
  diag <- fit_diag(fit)
  diag$start <- attr(fit, "start")
  diag$min_L_diag <- min(coef(fit)[c("L_11", "L_22")])
  diag$plateau <- diag$min_L_diag < -5
  if (vs_seed_rep("F4", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = refit, seed = seed)
  }
  diag
}

f4_panel <- function(sim, cell, seed) {
  what <- "the panel fit"
  fit <- fit_two_starts(function(theta_init) {
    f4_run(sim, seed, cell$S, person_col = "person", theta_init = theta_init)
  }, sim$truth, what)
  refit <- function(theta_init) {
    f4_fit(sim, seed, cell$S, what, person_col = "person", theta_init = theta_init)
  }
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  cf <- coef(fit)
  fit_r <- f4_fit(sim, seed, cell$S, "the LR test's restricted fit", x = setdiff(SP_X, "irr"),
                  person_col = "person", theta_init = unname(cf[names(cf) != "irr"]))
  diag <- .f4_diag(fit, refit, cell, seed)
  diag$lr_restricted <- unlist(fit_diag(fit_r))
  list(est = est_table(fit, sim$truth, V = V),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(cf[.f4_common(fit)], sim$truth, V$hessian),
       lr = lr_test(fit$loglik, fit_r$loglik, 1),
       derived = f4_estimands(fit, sim),
       diag = diag)
}

f4_xsec <- function(sim, cell, seed) {
  what <- "the composite fit"
  fit <- fit_two_starts(function(theta_init) {
    f4_run(sim, seed, cell$S_xsec, cluster_col = "person", theta_init = theta_init)
  }, sim$truth, what)
  refit <- function(theta_init) {
    f4_fit(sim, seed, cell$S_xsec, what, cluster_col = "person", theta_init = theta_init)
  }
  V <- fit_vcovs(fit, c("hessian", "cluster"))
  common <- coef(fit)[.f4_common(fit)]
  fit_h <- fit
  fit_h$vcov <- V$hessian
  list(est = est_table(fit, sim$truth, V = V),
       vcov_c = vcov_common(fit, route = "cluster", V = V$cluster),
       joint_wald = list(hessian = joint_wald(common, sim$truth, V$hessian),
                         cluster = joint_wald(common, sim$truth, V$cluster)),
       derived = f4_estimands(fit_h, sim),
       diag = .f4_diag(fit, refit, cell, seed))
}

family <- list(
  id = "F4",
  claims = c("C1", "C2", "C4"),
  dgp = "D1",
  cells = data.frame(cell = c("200", "800", "3.2k", "200.q"), N = c(200L, 800L, 3200L, 200L),
                     T = 8L, S = c(500L, 500L, 500L, 200L),
                     S_xsec = c(800L, 800L, NA_integer_, 200L), stringsAsFactors = FALSE),
  methods = list(panel = f4_panel, xsec = f4_xsec),
  methods_by_cell = list("200" = c("panel", "xsec"), "800" = c("panel", "xsec"),
                         "3.2k" = "panel", "200.q" = c("panel", "xsec")),
  R = list(quick = c("200.q" = 100), mc = c("200" = 1000, "800" = 1000, "3.2k" = 1000)),
  simulate = f4_simulate,
  estimands = f4_estimands,
  gates = list(cell = list(panel = "3.2k", xsec = "800"),
               metrics = list(panel = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                        "tail_upper:hessian", "joint_wald_size", "lr_size",
                                        paste0(c("tail_lower:hessian@", "tail_upper:hessian@"),
                                               rep(c("Sigma_11", "Sigma_21", "Sigma_22"),
                                                   each = 2L))),
                              xsec = c("bias_std", "se_ratio:cluster", "tail_lower:cluster",
                                       "tail_upper:cluster"))),
  # The quick tier at 200.q (S1; 200 persons, S = 200 for both methods), R = 300 on the family's
  # seeds, 2026-11-05 (MCSE 0.058 for bias, 0.04 for an SE ratio; no failed rep; the composite's
  # default start lost in 239 of 300 fits, and 20 of its kept fits ended on a plateau); reps 1-100,
  # Q3's own, pass every row. Kept: bias of beta and mu in both methods, and of the constants in the
  # panel (|value| at most 0.099); the SE ratio of price and mu, Hessian in the panel (0.993, 1.03,
  # 0.959), cluster in the composite (0.963, 1.03, 0.966); the panel's WTP SE ratio (1.01). Dropped,
  # all at this rung (C5): irr's SE ratio, 0.906 +- 0.037 (panel), 0.925 +- 0.038 (cluster), which
  # alone fail about 5% of R = 100 runs; the constants' SE ratio (0.93); every tail (price 2.0 /
  # 4.0%, irr 3.3 / 4.7%, Mu_neg_time 5.7 / 0%: one-sided errors O(N^-1/2)); the L block (L_22's bias
  # -0.327 +- 0.058, L_21's SE ratio 0.805, L_11's lower tail 9.3%) and Sigma (SE ratio 0.814, tails
  # to 10.3%); the joint Wald size, 22.0% panel, 40.3% composite; the LR size, 8.0% +- 1.6 (mean LR /
  # q 1.22). Q2 checks the L block at 3.2k persons.
  # False failures at R = 100, 1,000 bootstrap samples of the 300: 2.2%, recentered 0.7%.
  gates_quick = list(cell = "200.q",
                     metrics = list(panel = c("bias_std@beta", "bias_std@mu", "bias_std@asc",
                                              "se_ratio:hessian@price", "se_ratio:hessian@mu",
                                              "se_ratio:hessian@wtp_quality"),
                                    xsec = c("bias_std@beta", "bias_std@mu",
                                             "se_ratio:cluster@price", "se_ratio:cluster@mu")))
)
