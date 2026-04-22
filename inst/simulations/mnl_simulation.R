# Multinomial Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mnl_simulation.R

library(choicer)

# DGP =========================================================================
sim <- simulate_mnl_data(N = 50000, J = 5, seed = 123)
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
