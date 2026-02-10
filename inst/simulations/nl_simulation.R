# Nested Logit - Parameter Recovery Simulation
# Run from package root: Rscript inst/simulations/nl_simulation.R

source("inst/simulations/helpers.R")

# DGP =========================================================================
dgp <- simulate_nl_data(N = 5e+4, seed = 123)
cat("NL simulation: N =", dgp$settings$N,
    ", J_inside =", dgp$settings$J_inside, "\n\n")

# Data preparation ============================================================
# Add nest column to data: outside option (j=0) -> "outside", inside alts -> "nest_g"
dt <- dgp$data
dt[, nest := "outside"]
for (g in seq_along(dgp$nest_structure)) {
  dt[j %in% dgp$nest_structure[[g]], nest := paste0("nest_", g)]
}

# Estimation via convenience workflow ==========================================
fit <- run_nestlogit(
  data           = dt,
  id_col         = "id",
  alt_col        = "j",
  choice_col     = "choice",
  covariate_cols = c("X", "W"),
  nest_col       = "nest",
  use_asc        = TRUE,
  include_outside_option = TRUE,
  outside_opt_label = 0L,
  control = list(print_level = 1L)
)

cat("\n")
summary(fit)

# Parameter Recovery ==========================================================
theta <- coef(fit)
K_x   <- 2L
K_l   <- length(dgp$true_params$lambdas)

beta_est   <- theta[1:K_x]
lambda_est <- theta[(K_x + 1):(K_x + K_l)]
delta_est  <- theta[(K_x + K_l + 1):length(theta)]

cat("\n--- Beta Recovery ---\n")
compare_params(dgp$true_params$beta, beta_est, c("X", "W"))

cat("\n--- Lambda Recovery ---\n")
compare_params(
  dgp$true_params$lambdas, lambda_est,
  paste0("lambda_", seq_along(dgp$true_params$lambdas))
)

cat("\n--- Delta (ASC) Recovery ---\n")
compare_params(
  dgp$true_params$delta, delta_est,
  paste0("alt_", seq_along(dgp$true_params$delta))
)
