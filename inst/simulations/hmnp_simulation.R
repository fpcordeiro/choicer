# Hierarchical Bayesian Multinomial Probit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/hmnp_simulation.R

library(choicer)

# DGP =========================================================================
# Same hierarchical structure as the HMNL demo but with iid N(0, sigma^2)
# utility shocks (including on the implicit outside option). The probit
# likelihood identifies parameters only up to the common scale sigma, so
# simulate_hmnp_data() reports truth on the identified scale (beta/sigma,
# W/sigma^2, theta/sigma, sigma_d/sigma, delta/sigma) - with sigma = 1 the
# scales coincide.
sim <- simulate_hmnp_data(
  N = 1000, T = 12, J = 6,
  beta = c(0.8, -0.6), theta = c(0.3, -0.5), sigma_d = 0.5, sigma = 1,
  seed = 123
)
print(sim)

# Estimation ==================================================================
# Fully conjugate Albert-Chib Gibbs in un-differenced utility space: TN
# sweep over inside rows plus the stochastic outside latent, conjugate
# beta_i / delta_j / hierarchy draws, and a free (non-identified) sigma^2
# chain whose every kept draw is normalized away (PX-DA). No Metropolis
# steps, so no acceptance rates to tune.
set.seed(42)
fit <- run_hmnprobit(
  data               = sim$data,
  id_col             = "task",
  alt_col            = "alt",
  choice_col         = "choice",
  covariate_cols     = c("x1", "x2"),
  person_col         = "pid",
  alt_covariate_cols = "z1",
  chains             = 2,
  mcmc               = list(R = 8000, burn = 2000, thin = 2, trace = 2000)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
cat("\n--- Parameter Recovery (identified scale) ---\n")
print(recovery_table(fit, sim))

# Chain Diagnostics ===========================================================
ess <- function(x) {
  n <- length(x)
  rho <- as.numeric(stats::acf(x, lag.max = min(200L, n - 1L),
                               plot = FALSE)$acf)[-1]
  cut <- which(rho < 0)[1]
  if (!is.na(cut)) rho <- rho[seq_len(cut - 1L)]
  n / (1 + 2 * sum(rho))
}
draws <- cbind(fit$draws$b, fit$draws$theta,
               sigma_d2 = fit$draws$sigma_d2)
cat("\n--- Effective sample sizes (kept draws:", nrow(fit$draws$b), ") ---\n")
print(round(vapply(as.data.frame(draws), ess, numeric(1))))
cat("\n--- split-R-hat (2 chains, identified draws) ---\n")
print(round(fit$rhat, 3))

# The raw sigma^2 trace: non-identified, so it wanders (parameter
# expansion); only ratios like b/sigma are stationary.
cat("\nraw sigma^2 trace quantiles:",
    round(stats::quantile(fit$draws$sigma2, c(0.05, 0.5, 0.95)), 3), "\n")

# Scale relation to the HMNL ==================================================
# EV1 shocks have standard deviation pi/sqrt(6) ~ 1.28; iid probit shocks
# have (identified) standard deviation 1. Fitting BOTH models to their own
# DGPs with the same structural coefficients, the HMNL posterior b is
# therefore ~ pi/sqrt(6) times the HMNP posterior b. (Both models carry a
# stochastic outside option, which is what makes the comparison
# meaningful.)
cat("\n--- HMNP b scaled by pi/sqrt(6) (approximates an HMNL b on EV1 data) ---\n")
print(round(rbind(hmnp_b = coef(fit),
                  scaled = coef(fit) * pi / sqrt(6)), 3))

# Counterfactual via the quadrature engine ====================================
# Choice probabilities use the exact 1-D Gauss-Hermite representation of
# the iid-probit integral (no simulation noise in the probability
# evaluation itself).
cf <- data.table::copy(sim$data)
cf[cf$alt == 1L, "x1"] <- cf[cf$alt == 1L, ][["x1"]] + 0.25
set.seed(7)
base_share <- predict(fit, n_draws = 100)
set.seed(7)
cf_share <- predict(fit, newdata = cf, n_draws = 100)
cat("\n--- Subsidy counterfactual: alternative 1 share ---\n")
cat("baseline:", round(base_share[base_share$alternative == "1", ][["share"]], 4),
    "-> subsidized:", round(cf_share[cf_share$alternative == "1", ][["share"]], 4), "\n")

# Posterior-predictive check ==================================================
set.seed(7)
cat("\n--- Posterior-predictive share check ---\n")
print(ppc_shares(fit, n_draws = 100))
