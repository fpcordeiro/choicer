# ============================================================================
# F3: mixed logit, cross-section (plan §5; claims C1, C4, C5) on D4-40:
# school choice among each student's 15 nearest of one market's 40 schools
# (40 constants) and an outside option, a log-normal taste on -distance
# beside a fixed distance coefficient and a mean-zero normal taste on
# quality (rc_mean = FALSE, uncorrelated; lib/dgp_geo.R's header). S1
# (2026-11-02): 40 schools, not D4-small's 160, so that the gates test C1
# and not the incidental constants (C6): with 160, the distance
# coefficient's bias/SD was -0.163 J/sqrt(N) (-0.74 at 1k, -0.41 at 4k),
# vanishing when the constants were known. Here J/sqrt(N) = 1.26, 0.63, 0.32.
# Cells N in {1k, 4k, 16k} students, generate-mode draws seeded from the
# rep's seed, S = 500, 500, 520 (sqrt(N)/S <= 0.25), and "1k.q" (N = 1000,
# S = 200) for the quick tier.
# Every fit runs from choicer's default start and from theta0 and keeps the
# certified fit with the higher log-likelihood (fit_two_starts(), lib/fit.R;
# record in diag$start, with the smaller log-diagonal L_kk and the plateau
# flag L_kk < -5: from the default start a few fits at 1k ended with
# sigma_d -> 0, certified by the polish). Methods:
#   ml   run_mxlogit(polish = "newton"), every cell; routes hessian, bhhh,
#        robust; multi-start on reps 10, 20, ...;
#   blp  concentrate_out_asc = "blp" with the sample's shares (polish
#        "none", no bounds), at 16k only (R = 500, the same reps' data);
#        route hessian. Its theta0 start is theta0 whole: the constants'
#        part only starts the first share inversion, whose solution
#        delta(theta_c) does not depend on it.
# Derived, both methods: Sigma_11 and Sigma_22, the variances of L eta
# (sigma_vech(), delta method from the sigma block's covariance), and
# mwtp_q_ses = b_qs / (1 - b_d), the median over students of the WTP for
# quality x SES in km: the distance coefficient is b_d - exp(g), whose median
# is b_d - 1, so the WTP's median is at the median taste (delta SE from the
# fit's Hessian covariance; wtp() refuses a random price, and here the
# random coefficient is on the numeraire). Never-chosen alternatives are
# dropped in every rep (drop_unbounded_alts(); "blp" requires it; counts in
# diag$drop).
# Gates at 16k, ml: bias, Hessian SE ratio and tails for beta and L; joint
# Wald size on beta and L; Sigma's tails. Findings: rates over N, "blp"
# against ML (ESE ratio), the median WTP, the global-max rate.
# ============================================================================

f3_simulate <- function(cell, seed) {
  geo_drop(get_dgp("D4-40")$simulate(N = cell$N, seed = seed, variant = "default"))
}

# One fit, as it comes (fit_two_starts() certifies it); f3_fit() certifies
# a single fit (the multi-start's refits).
f3_run <- function(sim, cell, seed, ...) {
  choicer::run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W, rc_dist = c(1L, 0L),
                       rc_mean = FALSE, include_outside_option = TRUE, S = cell$S,
                       draws = "generate", seed = seed + 1L, ...)
}

f3_fit <- function(sim, cell, seed, ...) {
  require_converged(f3_run(sim, cell, seed, ...), "the MXL fit")
}

f3_estimands <- function(fit, sim) {
  cf <- coef(fit)
  V <- fit$vcov
  sg <- names(cf)[fit$param_map$sigma]
  K_w <- fit$draws_info$K_w
  s_hat <- sigma_vech(cf[sg], K_w, fit$rc_correlation, if (!is.null(V)) V[sg, sg])
  s_hat$truth <- sigma_vech(sim$truth[sg], K_w, fit$rc_correlation)$est
  mwtp <- function(b) b[["q_ses"]] / (1 - b[["dist"]])
  g <- c(q_ses = 1 / (1 - cf[["dist"]]), dist = cf[["q_ses"]] / (1 - cf[["dist"]])^2)
  se <- if (is.null(V)) NA_real_ else sqrt(drop(g %*% V[names(g), names(g)] %*% g))
  rbind(s_hat, data.frame(name = "mwtp_q_ses", est = mwtp(cf), se = se,
                          truth = mwtp(sim$truth), stringsAsFactors = FALSE))
}

# The common parameters (everything but the constants), by name.
.f3_common <- function(fit) names(coef(fit))[-as.integer(fit$param_map$asc)]

# fit_two_starts()'s record, the smaller log-diagonal parameter and the
# plateau flag, the counts dropped.
.f3_diag <- function(fit, sim) {
  diag <- fit_diag(fit)
  diag$start <- attr(fit, "start")
  diag$min_L_diag <- min(coef(fit)[names(coef(fit))[fit$param_map$sigma]])
  diag$plateau <- diag$min_L_diag < -5
  diag$drop <- sim$meta$drop
  diag
}

f3_ml <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f3_run(sim, cell, seed, polish = "newton", theta_init = theta_init)
  }, sim$truth, "the ML fit")
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  diag <- .f3_diag(fit, sim)
  if (vs_seed_rep("F3", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = function(theta_init) {
      f3_fit(sim, cell, seed, polish = "newton", theta_init = theta_init)
    }, seed = seed)
  }
  list(est = est_table(fit, sim$truth, V = V, gated = c("beta", "sigma")),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[.f3_common(fit)], sim$truth, V$hessian),
       derived = f3_estimands(fit, sim),
       diag = diag)
}

f3_blp <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f3_run(sim, cell, seed, concentrate_out_asc = "blp", theta_init = theta_init)
  }, sim$truth, "the blp fit")
  V <- fit_vcovs(fit, "hessian")
  conc <- fit$optimizer$concentration
  diag <- .f3_diag(fit, sim)
  diag$concentration <- unlist(conc[c("evaluations", "solves", "failed", "max_residual")])
  list(est = est_table(fit, sim$truth, V = V, gated = FALSE),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[.f3_common(fit)], sim$truth, V$hessian),
       derived = f3_estimands(fit, sim),
       diag = diag)
}

family <- list(
  id = "F3",
  claims = c("C1", "C4", "C5"),
  dgp = "D4-40",
  cells = data.frame(cell = c("1k", "4k", "16k", "1k.q"), N = c(1000L, 4000L, 16000L, 1000L),
                     S = c(500L, 500L, 520L, 200L), stringsAsFactors = FALSE),
  methods = list(ml = f3_ml, blp = f3_blp),
  methods_by_cell = list("1k" = "ml", "4k" = "ml", "16k" = c("ml", "blp"), "1k.q" = "ml"),
  R = list(quick = c("1k.q" = 100), mc = c("1k" = 1000, "4k" = 1000, "16k" = 1000)),
  R_by_method = list("16k" = c(blp = 500)),
  simulate = f3_simulate,
  estimands = f3_estimands,
  gates = list(cell = "16k",
               metrics = list(ml = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                     "tail_upper:hessian", "joint_wald_size",
                                     paste0(c("tail_lower:hessian@", "tail_upper:hessian@"),
                                            rep(c("Sigma_11", "Sigma_22"), each = 2L)))))
)
