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
