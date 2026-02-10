# Mixed Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/mxl_simulation.R

source("inst/simulations/helpers.R")

# DGP =========================================================================
dgp <- simulate_mxl_data(N = 5e+4, J = 8, seed = 123)
cat("MXL simulation: N =", dgp$settings$N, ", J =", dgp$settings$J, "\n\n")

# Data preparation ============================================================
mxl_inputs <- prepare_mxl_data(
  data             = dgp$data,
  id_col           = "id",
  alt_col          = "alt",
  choice_col       = "choice",
  covariate_cols   = c("x1", "x2"),
  random_var_cols  = c("w1", "w2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  rc_correlation = TRUE
)

# Halton draws: S > sqrt(N) for MSL consistency
S <- 50 * ceiling((1.5 * sqrt(dgp$settings$N)) / 50)
eta_draws <- get_halton_normals(S, mxl_inputs$N, dgp$settings$K_w)
cat("Simulation draws: S =", S, "\n\n")

# Estimation ==================================================================
fit <- run_mxlogit(
  input_data = mxl_inputs,
  eta_draws  = eta_draws,
  use_asc    = TRUE,
  control = list(print_level = 1L)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
theta  <- coef(fit)
K_x    <- dgp$settings$K_x
K_w    <- dgp$settings$K_w
L_size <- K_w * (K_w + 1) / 2

beta_est    <- theta[1:K_x]
L_params_est <- theta[(K_x + 1):(K_x + L_size)]
delta_est   <- theta[(K_x + L_size + 1):length(theta)]
Sigma_est   <- build_var_mat(L_params_est, K_w, TRUE)

cat("\n--- Beta Recovery ---\n")
compare_params(dgp$true_params$beta, beta_est, c("x1", "x2"))

cat("\n--- Sigma Recovery ---\n")
cat("True Sigma:\n")
print(dgp$true_params$Sigma)
cat("Estimated Sigma:\n")
print(round(Sigma_est, 4))

cat("\n--- Delta (ASC) Recovery ---\n")
compare_params(
  dgp$true_params$delta, delta_est,
  paste0("alt_", seq_along(dgp$true_params$delta))
)
