# Multinomial Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mnl_simulation.R

source("inst/simulations/helpers.R")

# DGP =========================================================================
dgp <- simulate_mnl_data(N = 50000, J = 5, seed = 123)
cat("MNL simulation: N =", dgp$settings$N, ", J =", dgp$settings$J, "\n\n")

# Estimation ==================================================================
fit <- run_mnlogit(
  data           = dgp$data,
  id_col         = "id",
  alt_col        = "alt",
  choice_col     = "choice",
  covariate_cols = c("x1", "x2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  use_asc = TRUE,
  control = list(print_level = 1L)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
theta <- coef(fit)
K_x   <- dgp$settings$K_x
beta_est  <- theta[1:K_x]
delta_est <- theta[(K_x + 1):length(theta)]

cat("\n--- Beta Recovery ---\n")
compare_params(dgp$true_params$beta, beta_est, c("x1", "x2"))

cat("\n--- Delta (ASC) Recovery ---\n")
compare_params(
  dgp$true_params$delta, delta_est,
  paste0("alt_", seq_along(dgp$true_params$delta))
)
