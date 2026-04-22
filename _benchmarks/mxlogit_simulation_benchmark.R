# MXL cross-package benchmark
# Compares choicer against logitr (and optionally mixl/gmnl/apollo) on a common
# simulated DGP that matches the `_benchmarks/` preset:
#   - log-normal on w1 (NEGATIVE price), normal on w2
#   - correlated random coefficients via `Sigma_true`
#   - mean shifters `mu_true`
# Run from package root: Rscript _benchmarks/mxlogit_simulation_benchmark.R

rm(list = ls(all.names = TRUE))
gc()

devtools::load_all()
source("_benchmarks/bench_helpers.R")

# Simulated DGP ----------------------------------------------------------------
sim <- simulate_mxl_data(
  N          = 1e4,
  J          = 6,
  beta       = c(0.8, -0.6),
  mu         = c(0.3, 0.7),
  Sigma      = matrix(c(1.5, 0.5,
                        0.5, 1.0), nrow = 2),
  rc_dist    = c(1L, 0L),       # log-normal on w1, normal on w2
  price_cols = "w1",            # w1 drawn negative
  seed       = 123
)
print(sim)

# Run benchmarks ---------------------------------------------------------------
# logitr is the only reliably installed MXL reference; extend as the environment
# gains `mixl`, `gmnl`, `apollo`.
res <- benchmark_fit(
  sim,
  packages = c("choicer", "logitr")
)

print(res)

# Parameter recovery (choicer) -------------------------------------------------
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
S <- 50 * ceiling((1.5 * sqrt(sim$settings$N)) / 50)
eta_draws <- get_halton_normals(S, mxl_inputs$N, sim$settings$K_w)

fit <- run_mxlogit(
  input_data = mxl_inputs,
  eta_draws  = eta_draws,
  rc_dist    = sim$true_params$rc_dist,
  rc_mean    = TRUE,
  use_asc    = TRUE,
  control    = list(print_level = 0L)
)
cat("\n--- choicer Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))
