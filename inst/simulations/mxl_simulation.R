# Mixed Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mxl_simulation.R

library(choicer)

# DGP =========================================================================
sim <- simulate_mxl_data(N = 5e3, J = 8, seed = 123)
print(sim)

# Data preparation ============================================================
mxl_inputs <- prepare_mxl_data(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  random_var_cols        = c("w1", "w2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  rc_correlation         = TRUE
)

# Halton draws: S > sqrt(N) for MSL consistency
S <- 50 * ceiling((1.5 * sqrt(sim$settings$N)) / 50)
eta_draws <- get_halton_normals(S, mxl_inputs$N, sim$settings$K_w)
cat("Simulation draws: S =", S, "\n\n")

# Estimation ==================================================================
fit <- run_mxlogit(
  input_data = mxl_inputs,
  eta_draws  = eta_draws,
  use_asc    = TRUE,
  control    = list(print_level = 1L, check_derivatives = TRUE)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
cat("\n--- Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))
