# Bayesian Multinomial Probit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mnp_simulation.R

library(choicer)

# DGP =========================================================================
# Latent utility differences against the base alternative (alt 1) with
# correlated errors. The default Sigma has sigma_11 = 1, so the DGP scale is
# already the identified scale that run_mnprobit() reports.
sim <- simulate_mnp_data(N = 5000, J = 3, seed = 123)
print(sim)

# Estimation ==================================================================
# Gibbs sampler with data augmentation (Albert-Chib 1993; McCulloch-Rossi
# 1994). The chain runs on the non-identified parameterization; reported
# quantities are normalized per draw by sigma_11. The master RNG seed is
# drawn from R's RNG, so set.seed() makes the whole run reproducible - and
# draws are bitwise invariant to the number of OpenMP threads.
set.seed(42)
fit <- run_mnprobit(
  data           = sim$data,
  id_col         = "id",
  alt_col        = "alt",
  choice_col     = "choice",
  covariate_cols = c("x1", "x2"),
  use_asc        = TRUE,
  mcmc           = list(R = 20000, burn = 5000, thin = 5, trace = 5000)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
# Posterior means / SDs of the identified draws against the identified truth.
# The Sigma_11 row is the identification normalization and is exact by
# construction (se = 0); lower_ci / upper_ci are normal-approximation
# credible intervals. Expect the betas to recover tightly. The ASC and Sigma
# posteriors are the weakly identified part of the MNP (the data carry only
# ordering information about the latent scale), so at this N they can sit
# 1-2 posterior SDs from truth for a single realization; they concentrate on
# the truth as N grows. Single-realization `covers` is one Bernoulli draw -
# coverage statements need a Monte Carlo over replications.
cat("\n--- Parameter Recovery ---\n")
print(recovery_table(fit, sim))

# Chain Diagnostics ===========================================================
# Effective sample size from the kept (thinned) draws, using the standard
# initial-positive-sequence truncation of the autocorrelation sum:
# ESS = n / (1 + 2 * sum(rho_k)). Gibbs draws are autocorrelated, so ESS,
# not the number of kept draws, measures the information in the chain;
# covariance elements (Sigma) typically mix slower than the coefficients.
ess <- function(x) {
  n <- length(x)
  rho <- as.numeric(stats::acf(x, lag.max = min(200L, n - 1L), plot = FALSE)$acf)[-1]
  cut <- which(rho < 0)[1]
  if (!is.na(cut)) rho <- rho[seq_len(cut - 1L)]
  n / (1 + 2 * sum(rho))
}

draws <- cbind(fit$draws$beta, fit$draws$sigma[, -1L, drop = FALSE])  # Sigma_11 == 1
diagnostics <- data.frame(
  parameter = colnames(draws),
  acf_lag1  = round(apply(draws, 2, function(x) stats::acf(x, lag.max = 1, plot = FALSE)$acf[2]), 3),
  ess       = round(apply(draws, 2, ess)),
  row.names = NULL
)
cat("\n--- Chain diagnostics (kept draws:", fit$mcmc$R_keep, ") ---\n")
print(diagnostics)

# Identified covariance ========================================================
# Posterior mean of Sigma / sigma_11 against the truth. Covariance parameters
# are the weakly identified part of the MNP: expect visibly wider posteriors
# than for beta, tightening only slowly with N.
cat("\n--- Posterior mean Sigma (identified) vs truth ---\n")
cat("Posterior mean:\n")
print(round(fit$sigma, 3))
cat("Truth:\n")
print(sim$true_params$Sigma)
