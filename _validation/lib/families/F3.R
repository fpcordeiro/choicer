# ============================================================================
# F3: mixed logit, cross-section (plan §5; claims C1, C4, C5, C6) on
# D4-small: school choice among each student's 15 nearest schools of 4
# markets x 40 (160 constants) and an outside option, a log-normal taste on
# -distance beside a fixed distance coefficient and a mean-zero normal taste
# on quality (rc_mean = FALSE, uncorrelated; lib/dgp_geo.R's header).
# Cells N in {1k, 4k, 16k} students at S = 500 generate-mode draws (seeded
# from the rep's seed), and "1k.q" (N = 1000, S = 200) for the quick tier.
# Methods:
#   ml   run_mxlogit(polish = "newton"), every cell; routes hessian, bhhh,
#        robust; multi-start on reps 10, 20, ...;
#   blp  concentrate_out_asc = "blp" with the sample's shares (polish
#        "none", no bounds), at 16k only (R = 500, the same reps' data);
#        route hessian.
# Derived, both methods: Sigma_11 and Sigma_22, the variances of L eta
# (Sigma = L L', delta method from the sigma block's covariance), and
# mwtp_q_ses = b_qs / (1 - b_d), the median over students of the WTP for
# quality x SES in km: the distance coefficient is b_d - exp(g), whose median
# is b_d - 1, so the WTP's median is at the median taste (delta SE from the
# fit's Hessian covariance; wtp() refuses a random price, and here the
# random coefficient is on the numeraire). Never-chosen alternatives are
# dropped in every rep (drop_unbounded_alts(); "blp" requires it; counts in
# diag$drop). Every fit, multi-start refits included, must converge
# (require_converged(), lib/fit.R: ml's polish certifies it, blp's optimizer
# code).
# Gates at 16k, ml: bias, Hessian SE ratio and tails for beta and L; joint
# Wald size on beta and L; Sigma's tails. Findings: rates over N, "blp"
# against ML (ESE ratio), the median WTP, the global-max rate.
# ============================================================================

f3_simulate <- function(cell, seed) {
  geo_drop(get_dgp("D4")$simulate(N = cell$N, seed = seed, variant = "default"))
}

f3_fit <- function(sim, cell, seed, ...) {
  require_converged(choicer::run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W,
                                         rc_dist = c(1L, 0L), rc_mean = FALSE,
                                         include_outside_option = TRUE, S = cell$S,
                                         draws = "generate", seed = seed + 1L, ...),
                    "the MXL fit")
}

# vech of Sigma = L L' (row by row: Sigma_11, Sigma_21, Sigma_22, ...; the
# diagonal only for an uncorrelated L) from choicer's sigma block, with the
# delta-method SEs from V, the block's covariance (NA without one).
f3_sigma <- function(sigma, K_w, correlated, V = NULL) {
  L <- L_from_theta(unname(sigma), K_w, correlated)
  pij <- if (correlated) cbind(rep(seq_len(K_w), seq_len(K_w)), sequence(seq_len(K_w)))
         else cbind(seq_len(K_w), seq_len(K_w))
  Sig <- L %*% t(L)
  Jac <- matrix(0, nrow(pij), nrow(pij))
  for (p in seq_len(nrow(pij))) {
    dL <- matrix(0, K_w, K_w)
    dL[pij[p, 1L], pij[p, 2L]] <- if (pij[p, 1L] == pij[p, 2L]) L[pij[p, 1L], pij[p, 1L]] else 1
    Jac[, p] <- (dL %*% t(L) + L %*% t(dL))[pij]
  }
  se <- if (is.null(V)) rep(NA_real_, nrow(pij))
        else sqrt(pmax(diag(Jac %*% V %*% t(Jac)), 0))
  data.frame(name = sprintf("Sigma_%d%d", pij[, 1L], pij[, 2L]), est = Sig[pij], se = se,
             stringsAsFactors = FALSE)
}

f3_estimands <- function(fit, sim) {
  cf <- coef(fit)
  V <- fit$vcov
  sg <- names(cf)[fit$param_map$sigma]
  K_w <- fit$draws_info$K_w
  s_hat <- f3_sigma(cf[sg], K_w, fit$rc_correlation, if (!is.null(V)) V[sg, sg])
  s_hat$truth <- f3_sigma(sim$truth[sg], K_w, fit$rc_correlation)$est
  mwtp <- function(b) b[["q_ses"]] / (1 - b[["dist"]])
  g <- c(q_ses = 1 / (1 - cf[["dist"]]), dist = cf[["q_ses"]] / (1 - cf[["dist"]])^2)
  se <- if (is.null(V)) NA_real_ else sqrt(drop(g %*% V[names(g), names(g)] %*% g))
  rbind(s_hat, data.frame(name = "mwtp_q_ses", est = mwtp(cf), se = se,
                          truth = mwtp(sim$truth), stringsAsFactors = FALSE))
}

# The common parameters (everything but the constants), by name.
.f3_common <- function(fit) names(coef(fit))[-as.integer(fit$param_map$asc)]

f3_ml <- function(sim, cell, seed) {
  refit <- function(theta_init) f3_fit(sim, cell, seed, polish = "newton", theta_init = theta_init)
  fit <- refit(NULL)
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  diag <- c(fit_diag(fit), list(drop = sim$meta$drop))
  if (vs_seed_rep("F3", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = refit, seed = seed)
  }
  list(est = est_table(fit, sim$truth, V = V, gated = c("beta", "sigma")),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[.f3_common(fit)], sim$truth, V$hessian),
       derived = f3_estimands(fit, sim),
       diag = diag)
}

f3_blp <- function(sim, cell, seed) {
  fit <- f3_fit(sim, cell, seed, concentrate_out_asc = "blp")
  V <- fit_vcovs(fit, "hessian")
  conc <- fit$optimizer$concentration
  list(est = est_table(fit, sim$truth, V = V, gated = FALSE),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[.f3_common(fit)], sim$truth, V$hessian),
       derived = f3_estimands(fit, sim),
       diag = c(fit_diag(fit), list(
         drop = sim$meta$drop,
         concentration = unlist(conc[c("evaluations", "solves", "failed", "max_residual")]))))
}

family <- list(
  id = "F3",
  claims = c("C1", "C4", "C5", "C6"),
  dgp = "D4",
  cells = data.frame(cell = c("1k", "4k", "16k", "1k.q"), N = c(1000L, 4000L, 16000L, 1000L),
                     S = c(500L, 500L, 500L, 200L), stringsAsFactors = FALSE),
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
