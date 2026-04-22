# MNL cross-package benchmark
# Compares choicer against mlogit and logitr on a common simulated DGP.
# Run from package root: Rscript _benchmarks/mnlogit_simulation_benchmark.R

rm(list = ls(all.names = TRUE))
gc()

devtools::load_all()
source("_benchmarks/bench_helpers.R")

# Simulated DGP ----------------------------------------------------------------
sim <- simulate_mnl_data(
  N    = 1e4,
  J    = 50,
  beta = c(0.8, -0.6),
  seed = 123
)
print(sim)

# Run benchmarks ---------------------------------------------------------------
res <- benchmark_fit(
  sim,
  packages = c("choicer", "mlogit", "logitr")
)

print(res)

# Parameter recovery (choicer rows only) ---------------------------------------
fit <- run_mnlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  use_asc                = TRUE,
  control                = list(print_level = 0L)
)
cat("\n--- choicer Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))
