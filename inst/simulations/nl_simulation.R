# Nested Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/nl_simulation.R

library(choicer)

# DGP =========================================================================
sim <- simulate_nl_data(N = 5e4, seed = 123)
print(sim)

# `sim$data` already contains a `nest` column (integer; 0L = outside);
# run_nestlogit() accepts integer nest labels directly.

# Estimation ==================================================================
fit <- run_nestlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "j",
  choice_col             = "choice",
  covariate_cols         = c("X", "W"),
  nest_col               = "nest",
  use_asc                = TRUE,
  include_outside_option = TRUE,
  outside_opt_label      = 0L,
  control                = list(print_level = 1L, check_derivatives = TRUE)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
cat("\n--- Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))

# Post-estimation =============================================================
# W plays the role of price (true beta_W = -0.8 < 0, so -alpha > 0).

# Goodness of fit. This DGP is balanced (every alternative in every choice
# set), so the closed-form market-shares (ASC-only) null is also valid.
cat("\n--- Goodness of fit ---\n")
print(gof(fit))
print(gof(fit, null = "market_shares"))

# Predicted market shares (includes the outside option implicitly:
# inside shares sum to 1 - outside share)
cat("\n--- Predicted shares ---\n")
print(predict(fit, type = "shares"))

# Elasticities and diversion ratios under the nested substitution pattern:
# within-nest diversion exceeds cross-nest diversion (no IIA across nests)
cat("\n--- Elasticities (W) ---\n")
print(elasticities(fit, elast_var = "W"))
cat("\n--- Diversion ratios ---\n")
print(diversion_ratios(fit))

# Willingness to pay for X in units of W, with delta-method SEs.
# True WTP_X = -1.5 / (-0.8) = 1.875. Lambdas are not WTP attributes.
cat("\n--- Willingness to pay ---\n")
print(wtp(fit, price_var = "W"))

# Expected consumer surplus from the nested logsum (point estimate for NL)
cs0 <- consumer_surplus(fit, price_var = "W")
cat("\n--- Consumer surplus (baseline) ---\n")
print(cs0)

# Counterfactual: raise W ("price") by 0.5 for alternative 2 only =============
dt_cf <- data.table::copy(sim$data)[j == 2, W := W + 0.5]

cat("\n--- Counterfactual shares (W + 0.5 for alt 2) ---\n")
shares_cf  <- predict(fit, type = "shares", newdata = dt_cf)
shares_cmp <- rbind(baseline       = drop(predict(fit, type = "shares")),
                    counterfactual = drop(shares_cf))
colnames(shares_cmp) <- as.character(fit$alt_mapping[[2]])
print(shares_cmp)

# Welfare change: negative, since the price of one alternative rose
cs1 <- consumer_surplus(fit, price_var = "W", newdata = dt_cf)
cat("\nDelta mean consumer surplus:", round(cs1$mean_cs - cs0$mean_cs, 4),
    "(negative: the price increase lowers expected surplus)\n")
