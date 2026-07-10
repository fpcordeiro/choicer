# Hierarchical Bayesian Multinomial Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/hmnl_simulation.R

library(choicer)

# DGP =========================================================================
# Panel logit with respondent-level random coefficients beta_i ~ N(b, W) and
# a BLP-style alternative-level random effect delta_j = z_j'theta + xi_j,
# xi_j ~ N(0, sigma_d^2), against an implicit outside option (systematic
# utility 0 plus its own EV1 shock). The outside good anchors the level of
# delta, so delta_j is mean utility relative to the outside option.
sim <- simulate_hmnl_data(
  N = 1000, T = 12, J = 6,
  beta = c(0.8, -0.6), theta = c(0.3, -0.5), sigma_d = 0.5,
  seed = 123
)
print(sim)

# Estimation ==================================================================
# Adaptive RW-Metropolis-within-Gibbs: parallel per-respondent beta_i
# updates (proposal precision H_i + W^-1 at the pooled MLE), a strictly
# serial delta_j sweep (the delta conditionals are coupled through the
# softmax denominators), and conjugate hierarchy draws. set.seed() makes
# the whole run reproducible given a fixed number of OpenMP threads (across
# thread counts, draws agree only up to floating-point reduction-order
# round-off). Two chains feed the split-R-hat table.
set.seed(42)
fit <- run_hmnlogit(
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
# Posterior means / SDs against the DGP truth: the population means b, the
# random-coefficient variances diag(W), the delta mean function theta, the
# alternative-effect scale sigma_d, and the realized delta_j ladder. The
# delta LEVEL is identified by the outside-option share, so it is the
# weakly-identified direction when the outside share is small; the delta
# SHAPE (cross-alternative contrasts) recovers tightly.
cat("\n--- Parameter Recovery ---\n")
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
cat("\n--- split-R-hat (2 chains) ---\n")
print(round(fit$rhat, 3))
cat("\nMean acceptance: beta", round(fit$accept$mean_beta, 3),
    "| delta", round(fit$accept$mean_delta, 3), "\n")

# Cross-check vs the frequentist mixed logit =================================
# The same DGP through run_mxlogit (random coefficients on x1/x2, ASCs
# absorbing delta): the MSLE mu should agree with the posterior mean b to
# within sampling noise. The frequentist model has no partial pooling on
# the ASCs, so the comparison is on the structural coefficients.
cat("\n--- Frequentist cross-check (run_mxlogit mu vs posterior b) ---\n")
# run_mxlogit requires at least one fixed covariate; a pure-noise column
# (true coefficient 0) keeps x1/x2 random and exactly comparable.
mxl_data <- data.table::copy(sim$data)
mxl_data$noise <- stats::rnorm(nrow(mxl_data))
mxl <- run_mxlogit(
  data = mxl_data, id_col = "task", alt_col = "alt", choice_col = "choice",
  covariate_cols = "noise", random_var_cols = c("x1", "x2"), rc_mean = TRUE,
  use_asc = TRUE, include_outside_option = TRUE, S = 200
)
comp <- rbind(hmnl_b = coef(fit), mxl_mu = coef(mxl)[c("Mu_x1", "Mu_x2")])
print(round(comp, 3))

# Policy counterfactual (compensating variation) =============================
# A subsidy that improves x1 by +0.25 for alternative 1 only. CV is the
# posterior of the logsum difference scaled by the marginal utility of
# income (here the negative x2 coefficient), summed over choice situations.
cf <- data.table::copy(sim$data)
cf[cf$alt == 1L, "x1"] <- cf[cf$alt == 1L, ][["x1"]] + 0.25
set.seed(7)
base_share <- predict(fit, n_draws = 200)
set.seed(7)
cf_share <- predict(fit, newdata = cf, n_draws = 200)
cat("\n--- Subsidy counterfactual: alternative 1 share ---\n")
cat("baseline:", round(base_share[base_share$alternative == "1", ][["share"]], 4),
    "-> subsidized:", round(cf_share[cf_share$alternative == "1", ][["share"]], 4), "\n")
set.seed(7)
cs <- consumer_surplus(fit, price_var = "x2", newdata = cf, n_draws = 200)
cat("Compensating variation (posterior quantiles, utility-of-income units):\n")
print(round(attr(cs, "cv"), 2))

# Entry counterfactual (posterior-predictive delta) ==========================
# A new alternative never seen in estimation, with alternative-level
# covariate z1 = 0.4: its delta comes from the posterior predictive
# N(z_new'theta, sigma_d^2), so its share carries the full entry
# uncertainty (a wider interval than the incumbents').
entrant <- sim$data[sim$data$alt == 1L, ]
entrant$alt <- 99L
entrant$z1 <- 0.4
entrant$choice <- 0L
entry_data <- rbind(sim$data, entrant)
set.seed(7)
p_entry <- predict(fit, newdata = entry_data, n_draws = 200)
cat("\n--- Entry counterfactual: shares with a new alternative (99) ---\n")
print(p_entry)

# Posterior-predictive check ==================================================
set.seed(7)
cat("\n--- Posterior-predictive share check ---\n")
print(ppc_shares(fit, n_draws = 200))
