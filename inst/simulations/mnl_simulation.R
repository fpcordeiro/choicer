# Multinomial Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mnl_simulation.R

library(choicer)

# DGP =========================================================================
sim <- simulate_mnl_data(N = 50000, J = 500, seed = 123)
print(sim)

# Estimation ==================================================================
fit <- run_mnlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  use_asc                = TRUE,
  control                = list(print_level = 1L)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
cat("\n--- Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))

# Post-estimation =============================================================
# x2 plays the role of price (true beta_x2 = -0.6 < 0, so -alpha > 0).

# Goodness of fit: McFadden R2 (vs. the equal-shares null, exact for the
# varying choice sets used by this DGP) and the in-sample hit rate. The same
# figures appear in the summary() footer above.
cat("\n--- Goodness of fit ---\n")
print(gof(fit))

# Predicted market shares
cat("\n--- Predicted shares ---\n")
print(predict(fit, type = "shares"))

# Own- and cross-elasticities of choice probabilities w.r.t. x2
cat("\n--- Elasticities (x2) ---\n")
print(elasticities(fit, elast_var = "x2"))

# Diversion ratios: share of demand lost by j that flows to k
cat("\n--- Diversion ratios ---\n")
print(diversion_ratios(fit))

# Willingness to pay for x1 (and the ASCs) in units of x2, with
# delta-method standard errors. True WTP_x1 = -0.8 / (-0.6) = 1.33.
cat("\n--- Willingness to pay ---\n")
print(wtp(fit, price_var = "x2"))

# Expected consumer surplus (logsum / -alpha). Levels depend on the ASC
# normalization; differences across scenarios are the meaningful quantity.
cs0 <- consumer_surplus(fit, price_var = "x2")
cat("\n--- Consumer surplus (baseline) ---\n")
print(cs0)

# Counterfactual: raise x2 ("price") by 0.5 for alternative 2 only =============
dt_cf <- data.table::copy(sim$data)[alt == 2, x2 := x2 + 0.5]

cat("\n--- Counterfactual shares (x2 + 0.5 for alt 2) ---\n")
shares_cf  <- predict(fit, type = "shares", newdata = dt_cf)
shares_cmp <- rbind(baseline       = drop(predict(fit, type = "shares")),
                    counterfactual = drop(shares_cf))
colnames(shares_cmp) <- as.character(fit$alt_mapping[[2]])
print(shares_cmp)

# Welfare change: negative, since the price of one alternative rose
cs1 <- consumer_surplus(fit, price_var = "x2", newdata = dt_cf)
cat("\nDelta mean consumer surplus:", round(cs1$mean_cs - cs0$mean_cs, 4),
    "(negative: the price increase lowers expected surplus)\n")
