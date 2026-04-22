# NL cross-package benchmark
# Compares choicer against mlogit on a common simulated DGP.
# Run from package root: Rscript _benchmarks/nestlogit_benchmark.R

rm(list = ls(all.names = TRUE))
gc()

devtools::load_all()
source("_benchmarks/bench_helpers.R")

# Simulated DGP ----------------------------------------------------------------
sim <- simulate_nl_data(N = 5e4, seed = 123)
print(sim)

# Run benchmarks ---------------------------------------------------------------
res <- benchmark_fit(
  sim,
  packages = c("choicer", "mlogit")
)

print(res)

# Parameter recovery (choicer) -------------------------------------------------
dt <- copy(sim$data)
dt[, nest := "outside"]
for (g in seq_along(sim$settings$nest_structure)) {
  dt[j %in% sim$settings$nest_structure[[g]], nest := paste0("nest_", g)]
}

fit <- run_nestlogit(
  data                   = dt,
  id_col                 = "id",
  alt_col                = "j",
  choice_col             = "choice",
  covariate_cols         = c("X", "W"),
  nest_col               = "nest",
  use_asc                = TRUE,
  include_outside_option = TRUE,
  outside_opt_label      = 0L,
  control                = list(print_level = 0L)
)

cat("\n--- choicer Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))
