#!/usr/bin/env Rscript
# Analysis: fit log-log slopes (empirical scaling exponents) per profiled
# block over the N sweep (fixed J) and the J sweep (fixed N), then
# extrapolate each block's per-iteration cost to the target scale
# (N=30000, J=200) and to the full feasibility budget (R=10000, burn=5000,
# 2 SEQUENTIAL chains -- see R/hmnlogit_utils.R / hmnprobit_utils.R's
# lapply-based chain loop).
#
# NOT executed as part of Phase 1 (the full grid has not been run yet). Kept
# here, ready to run once the grid in configs/proposed_grid.csv is approved
# and 20_harness.R has populated results/hb_profiling_results.csv.
#
# Usage:
#   Rscript _benchmarks/hb_profiling/30_analyze.R \
#     --results _benchmarks/hb_profiling/results/hb_profiling_results.csv

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}
this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file) == 1) dirname(this_file) else "."
results_path <- get_arg("--results",
  file.path(script_dir, "results", "hb_profiling_results.csv"))

stopifnot(file.exists(results_path))
res <- fread(results_path)
res <- res[status == "OK"]

block_cols <- grep("^prof_t_", names(res), value = TRUE)
message("Profiled block columns found: ", paste(block_cols, collapse = ", "))

# --- Log-log slope helper: fit log(block_time / R) ~ log(x) -----------------
loglog_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  if (sum(keep) < 2) return(c(slope = NA_real_, intercept = NA_real_, n = sum(keep)))
  fit <- lm(log(y[keep]) ~ log(x[keep]))
  c(slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]),
    n = sum(keep))
}

fit_sweep <- function(dt, sweep_var, id_prefix, model_name) {
  sub <- dt[model == model_name & grepl(paste0("^", id_prefix), config_id)]
  if (nrow(sub) < 2) return(NULL)
  out <- rbindlist(lapply(block_cols, function(bc) {
    per_iter <- sub[[bc]] / sub$R
    s <- loglog_slope(sub[[sweep_var]], per_iter)
    data.table(model = model_name, sweep = sweep_var, block = bc,
              slope = s["slope"], intercept = s["intercept"], n = s["n"])
  }))
  out
}

slopes <- rbindlist(c(
  lapply(c("hmnl", "hmnp"), function(m) fit_sweep(res, "N", "nsweep_", m)),
  lapply(c("hmnl", "hmnp"), function(m) fit_sweep(res, "J", "jsweep_", m))
), fill = TRUE)
slopes <- slopes[!sapply(slopes$slope, is.null)]

cat("\n=== Empirical log-log scaling exponents (per-iteration block cost) ===\n")
print(slopes)

# --- Extrapolate to the target scale + full feasibility budget --------------
# Target: N = 30000, J = 200, R = 10000, burn = 5000, 2 SEQUENTIAL chains.
target_N <- 30000; target_J <- 200
target_R <- 10000; target_chains <- 2

anchor <- res[grepl("^target_anchor_", config_id)]
extrapolate_from_anchor <- function(model_name) {
  a <- anchor[model == model_name]
  if (nrow(a) == 0) {
    message("No target_anchor_", model_name, " row in results yet; skipping.")
    return(NULL)
  }
  per_iter_cols <- block_cols
  per_iter <- sapply(per_iter_cols, function(bc) a[[bc]][1] / a$R[1])
  total_per_iter <- sum(per_iter, na.rm = TRUE)
  full_run_sec <- total_per_iter * target_R * target_chains
  data.table(model = model_name,
            anchor_N = a$N[1], anchor_J = a$J[1], anchor_R = a$R[1],
            sec_per_iter_at_anchor = total_per_iter,
            extrapolated_full_budget_sec = full_run_sec,
            extrapolated_full_budget_hours = full_run_sec / 3600)
}
extrap <- rbindlist(lapply(c("hmnl", "hmnp"), extrapolate_from_anchor), fill = TRUE)

cat("\n=== Extrapolation to N=30000, J=200, R=10000, burn=5000, 2 sequential chains ===\n")
print(extrap)

cat("\nNOTE: this extrapolation uses the target-scale anchor config directly ",
   "(N=30000, J=200 already run at small R) when available; the N/J sweep ",
   "slopes above are provided to sanity-check that anchor extrapolation ",
   "against the independently-fit power law.\n", sep = "")
