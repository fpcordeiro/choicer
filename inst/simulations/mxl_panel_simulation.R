# Panel Mixed Logit - Parameter and Taste Recovery Simulation
#
# Estimand: the population distribution of tastes f(gamma | theta) when each
# decision maker keeps one taste vector across repeated choices (Revelt and
# Train 1998), plus what a decision maker's own choices reveal about theirs.
#   1. A panel fit (person_col) recovers the fixed coefficients, the ASCs and
#      the Cholesky factor of the random-coefficient covariance.
#   2. Conditional tastes E[gamma_n | y_n] (Revelt and Train 2000) track the
#      realized tastes of the DGP, shrunk toward the population mean.
#   3. The same data fitted as a cross-section with decision-maker-clustered
#      standard errors: what the panel likelihood buys in precision.
#
# Run from package root: Rscript inst/simulations/mxl_panel_simulation.R

library(choicer)

# 1) DGP ======================================================================
# N decision makers x T choice situations each. The utility of inside
# alternative j in situation t of decision maker n is
#   U_ntj = delta_j + x_ntj'beta + w_ntj'gamma_n + eps_ntj,
# with fixed beta = (0.8, -0.6) on (x1, x2), ASCs delta = (0.5, -0.5, ...),
# zero-mean random coefficients gamma_n ~ N(0, Sigma) on (w1, w2) with
# Sigma = [1, 0.5; 0.5, 1.5], and iid EV1 shocks. Every choice set also holds
# a physical outside row (alt = 0, covariates 0, utility eps_nt0).
#
# Beyond beta and delta, common to everyone, the DGP holds only gamma_n fixed
# within a person: one draw per decision maker, kept in all T of their
# situations. The choice set (2 to J inside alternatives), the attributes
# and the shocks are redrawn in every situation, so a person's choices are
# dependent only through their tastes - the dependence the panel likelihood
# models and a cross-section ignores.
sim <- simulate_mxl_data(N = 2000, T = 10, J = 5, seed = 123)
print(sim)

# 2) Panel fit ================================================================
# person_col = "pid": decision maker n contributes one likelihood term that
# integrates the product of their T logit probabilities over a single taste
# draw,
#   L_n(theta) = int prod_t P_nt(j_nt | gamma) f(gamma | theta) d gamma,
# so the within-person correlation of choices identifies Sigma most sharply
# (the marginals identify it too, but only through the shape of the mixture).
#
# Draws: S Halton draws per decision maker, generated on the fly from
# `seed`. MSL matches ML asymptotically only if S grows faster than sqrt(N):
# the simulation bias is O(1/S), the sampling error O(1/sqrt(N)). The bias
# comes from taking the log of a simulated average, a gap that widens with
# the spread of tastes, so it typically pulls the variance estimates toward
# zero. A product of T probabilities is sharply peaked in gamma and needs
# far more draws than a single choice probability: S = 1000 >> sqrt(N) = 45.
#
# rc_mean = FALSE (the default) imposes the DGP's zero taste means. Alt 0 is
# the first alternative, so ASC_j is utility relative to the outside good.
S <- 1000L
fit <- run_mxlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  random_var_cols        = c("w1", "w2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  rc_correlation         = TRUE,
  person_col             = "pid",
  S                      = S,
  draws                  = "generate",
  seed                   = 2026L
)

cat("\n")
summary(fit)

# 3) Parameter recovery =======================================================
# The sigma block compares the Cholesky parameters the optimizer works in
# (log diagonal, raw off-diagonal): truth (0, 0.5, log(sqrt(1.25))). The
# summary above reports Sigma itself, by the delta method.
cat("\n--- Parameter Recovery (panel fit) ---\n")
rt_panel <- recovery_table(fit, sim$true_params)
print(rt_panel)

# 4) Conditional tastes =======================================================
# E[gamma_n | y_n, x_n; theta_hat] (printed as E[beta | observed choices]):
# the mean taste among decision makers who face n's attributes and make n's
# T choices. It is a posterior mean, not n's own gamma_n: it tracks the
# realized tastes but shrinks toward the population mean, so its SD across
# decision makers falls short of the SD of tastes (population SDs 1 and
# sqrt(1.5) = 1.22). The shortfall in variance is the taste variance that T
# choices leave unresolved (`residual`). Because E[gamma | y] is the best
# predictor of gamma given y, Cov(E[gamma | y], gamma) = Var(E[gamma | y]),
# so at the true theta
#   cor(E[gamma | y], gamma) = sd(E[gamma | y]) / sd(gamma):
# how well the conditional means track tastes and how much they shrink are
# one number. The last two columns below should agree up to sampling,
# simulation and estimation error in theta_hat.
ct <- conditional_tastes(fit)
cat("\n")
print(ct)

gamma_true <- sim$true_params$gamma_i[, colnames(ct$mean), drop = FALSE]
taste_fit <- data.frame(
  sd_cond  = apply(ct$mean, 1, stats::sd),
  sd_gamma = apply(gamma_true, 1, stats::sd),
  cor      = vapply(rownames(ct$mean), function(k) {
    stats::cor(ct$mean[k, ], gamma_true[k, ])
  }, numeric(1))
)
taste_fit$sd_ratio <- taste_fit$sd_cond / taste_fit$sd_gamma
taste_fit <- taste_fit[c("sd_cond", "sd_gamma", "sd_ratio", "cor")]
cat("\n--- Conditional means vs realized tastes (", ncol(ct$mean),
    " decision makers) ---\n", sep = "")
print(round(taste_fit, 3))

# 5) Contrast: the same data as a cross-section ===============================
# Without person_col every choice situation gets its own taste draw. On
# panel data that is a composite (pseudo-) likelihood: each situation's
# marginal choice probability is still a correct mixed logit probability, so
# it targets the same theta, but it discards the within-person correlation
# of choices and identifies Sigma only through the shape of the marginal
# probabilities. Its standard errors must be clustered by decision maker
# (cluster_col = "pid" selects se_method = "cluster"). Its log-likelihood
# sums marginal log-probabilities over situations rather than joint ones
# over decision makers, so it - and its AIC/BIC - is not comparable with the
# panel fit's. The same S = 1000 draws per choice situation also clear
# sqrt(N * T) = 141, the rule for this fit's N * T likelihood units.
fit_cs <- run_mxlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  random_var_cols        = c("w1", "w2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  rc_correlation         = TRUE,
  cluster_col            = "pid",
  S                      = S,
  draws                  = "generate",
  seed                   = 2026L
)
cat("Cross-sectional fit convergence:", fit_cs$convergence,
    "(", fit_cs$message, ")\n")

cat("\n--- Parameter Recovery (cross-sectional fit, clustered SEs) ---\n")
rt_cs <- recovery_table(fit_cs, sim$true_params)
print(rt_cs)

# SE ratio, cross-sectional (clustered) / panel (Hessian). Its square,
# n_factor, is roughly how many times as many decision makers the
# cross-sectional fit would need to match the panel's precision. The
# Cholesky block is where the panel likelihood pays off. The mean-utility
# parameters (beta, ASCs) are pinned down by marginal choice probabilities,
# which both fits use, so they gain little.
se_ratio <- data.frame(
  group    = rt_panel$group,
  se_panel = rt_panel$se,
  se_cross = rt_cs$se[match(rt_panel$parameter, rt_cs$parameter)],
  row.names = rt_panel$parameter
)
se_ratio$ratio <- se_ratio$se_cross / se_ratio$se_panel
se_ratio$n_factor <- se_ratio$ratio^2
num <- c("se_panel", "se_cross", "ratio", "n_factor")
se_ratio[num] <- lapply(se_ratio[num], round, 3)
cat("\n--- SE ratio: cross-sectional / panel ---\n")
print(se_ratio)
