#!/usr/bin/env Rscript
# Analysis of the full HB profiling grid (Phase 2): per-block log-log
# empirical scaling exponents over the N sweep and J sweep (including the
# HMNP fused-pass and HMNL beta-phase thread-summed sub-splits), a
# threads={1,5} Amdahl decomposition, the keep_beta_i=2 memory-probe vs the
# analytic K*N*R_keep*8B formula, predicted-vs-actual wall time per config,
# and extrapolation to the target scale + full feasibility budget
# (N=30000, J=200, K=10, P=10, R=10000, burn=5000, 2 SEQUENTIAL chains,
# 5 threads).
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
pred_path <- get_arg("--predictions",
  file.path(script_dir, "configs", "predicted_grid_times.csv"))

stopifnot(file.exists(results_path))
res <- fread(results_path)
res_ok <- res[status == "OK"]

wall_block_cols <- grep("^prof_t_", names(res_ok), value = TRUE)
cpu_split_cols <- grep("^prof_(fused|beta_(propconstruct|loglik|acceptadapt))_", names(res_ok), value = TRUE)
message("Wall-clock block columns: ", paste(wall_block_cols, collapse = ", "))
message("CPU-sum sub-split columns: ", paste(cpu_split_cols, collapse = ", "))

# =============================================================================
# 1. Log-log empirical scaling exponents (N sweep, J sweep)
# =============================================================================
loglog_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  if (sum(keep) < 2) return(c(slope = NA_real_, intercept = NA_real_, n = sum(keep)))
  fit <- lm(log(y[keep]) ~ log(x[keep]))
  c(slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]), n = sum(keep))
}

fit_sweep <- function(dt, sweep_var, id_prefix, model_name, cols) {
  sub <- dt[model == model_name & grepl(paste0("^", id_prefix), config_id)]
  if (nrow(sub) < 2) return(NULL)
  rbindlist(lapply(cols, function(bc) {
    if (!bc %in% names(sub)) return(NULL)
    per_iter <- sub[[bc]] / sub$R
    s <- loglog_slope(sub[[sweep_var]], per_iter)
    data.table(model = model_name, sweep = sweep_var, block = bc,
              slope = s["slope"], intercept = s["intercept"], n_points = s["n"])
  }), fill = TRUE)
}

all_cols <- c(wall_block_cols, cpu_split_cols)
slopes <- rbindlist(c(
  lapply(c("hmnl", "hmnp"), function(m) fit_sweep(res_ok, "N", "nsweep_", m, all_cols)),
  lapply(c("hmnl", "hmnp"), function(m) fit_sweep(res_ok, "J", "jsweep_", m, all_cols))
), fill = TRUE)
slopes <- slopes[!is.na(slope)]
slopes[, slope := round(slope, 3)]

cat("\n=== 1. Empirical log-log scaling exponents (per-iteration block cost) ===\n")
cat("(slope ~ 1 => linear in the sweep variable; ~ 2 => quadratic; etc.)\n")
print(slopes[order(model, sweep, block)], nrows = 100)

# =============================================================================
# 2. Threads {1, 5} Amdahl decomposition (N=4000, J=100)
# =============================================================================
cat("\n=== 2. Threads {1, 5} Amdahl decomposition (N=4000, J=100) ===\n")
amdahl <- function(model_name) {
  t1 <- res_ok[config_id == paste0("threads_", model_name, "_t1_N4000_J100")]
  t5 <- res_ok[config_id == paste0("threads_", model_name, "_t5_N4000_J100")]
  if (nrow(t1) == 0 || nrow(t5) == 0) return(NULL)

  if (model_name == "hmnl") {
    serial_col <- "prof_t_delta_serial"    # never parallelized
    parallel_cols <- c("prof_t_cache_rebuild", "prof_t_beta_mh")
  } else {
    serial_col <- "prof_t_hierarchy"       # master-only; recording negligible
    parallel_cols <- c("prof_t_fused_pass", "prof_t_delta_parallel", "prof_t_rss")
  }
  serial_t1 <- t1[[serial_col]] / t1$R
  serial_t5 <- t5[[serial_col]] / t5$R
  parallel_t1 <- sum(sapply(parallel_cols, function(c) t1[[c]] / t1$R))
  parallel_t5 <- sum(sapply(parallel_cols, function(c) t5[[c]] / t5$R))
  total_t1 <- t1$prof_t_total / t1$R
  total_t5 <- t5$prof_t_total / t5$R

  speedup_parallel_blocks <- parallel_t1 / parallel_t5
  speedup_total <- total_t1 / total_t5
  serial_share_t5 <- serial_t5 / total_t5
  # Amdahl ceiling at threads -> infinity, holding the measured serial
  # fraction (at 5 threads) fixed: max possible speedup vs 1 thread is
  # 1 / serial_fraction_at_1_thread.
  serial_fraction_1thread <- serial_t1 / total_t1
  amdahl_ceiling <- 1 / serial_fraction_1thread

  data.table(model = model_name,
            serial_block = serial_col,
            sec_per_iter_1thread = total_t1, sec_per_iter_5thread = total_t5,
            measured_speedup_total = speedup_total,
            measured_speedup_parallel_blocks_only = speedup_parallel_blocks,
            serial_fraction_at_1thread = serial_fraction_1thread,
            serial_share_of_iter_at_5threads = serial_share_t5,
            amdahl_ceiling_speedup_at_infinite_threads = amdahl_ceiling)
}
amdahl_dt <- rbindlist(lapply(c("hmnl", "hmnp"), amdahl), fill = TRUE)
print(amdahl_dt)

# =============================================================================
# 3. keep_beta_i = "draws" memory probe vs analytic K*N*R_keep*8B
# =============================================================================
cat("\n=== 3. Memory probe (keep_beta_i=draws, N=8000, J=50, R=400) vs analytic cube ===\n")
mem_probe <- function(model_name) {
  probe <- res_ok[config_id == paste0("mem_probe_", model_name)]
  # Same N, J, K, P but keep_beta_i = "means" and different R -- isolates the
  # base (N,J)-dependent footprint (X, Z, per-row caches) from the R_keep x
  # K x N draw cube added by keep_beta_i = "draws".
  base_cfg <- res_ok[config_id == paste0("nsweep_", model_name, "_N8000_J50")]
  if (nrow(probe) == 0 || nrow(base_cfg) == 0) return(NULL)

  K <- probe$K[1]; N <- probe$N[1]
  R_keep <- floor((probe$R[1] - probe$burn[1]) / 1) + 1L   # thin = 1
  # (R - burn + thin - 1) %/% thin, thin = 1:
  R_keep <- (probe$R[1] - probe$burn[1] + 1L - 1L) %/% 1L
  analytic_cube_bytes <- 8 * as.numeric(K) * N * R_keep

  data.table(model = model_name,
            probe_peak_rss_bytes = probe$peak_rss_bytes[1],
            base_peak_rss_bytes = base_cfg$peak_rss_bytes[1],
            measured_delta_bytes = probe$peak_rss_bytes[1] - base_cfg$peak_rss_bytes[1],
            analytic_cube_bytes = analytic_cube_bytes,
            R_keep = R_keep,
            ratio_measured_to_analytic = (probe$peak_rss_bytes[1] - base_cfg$peak_rss_bytes[1]) / analytic_cube_bytes)
}
mem_dt <- rbindlist(lapply(c("hmnl", "hmnp"), mem_probe), fill = TRUE)
print(mem_dt)

# =============================================================================
# 4. Predicted-vs-actual wall time per config (flag > 3x deviation)
# =============================================================================
cat("\n=== 4. Predicted-vs-actual wall time (flag ratio > 3x either bound) ===\n")
if (file.exists(pred_path)) {
  pred <- fread(pred_path)
  cmp <- merge(res_ok[model %in% c("hmnl", "hmnp"),
                      .(config_id, model, N, J, R, prof_t_total, wall_outer_sec)],
              pred[, .(config_id, pred_total_low_sec, pred_total_high_sec)],
              by = "config_id", all.x = TRUE)
  cmp[, ratio_vs_high := prof_t_total / pred_total_high_sec]
  cmp[, ratio_vs_low := prof_t_total / pred_total_low_sec]
  cmp[, flag_gt_3x := (ratio_vs_high > 3) | (ratio_vs_low > 3) |
                      (ratio_vs_high < 1/3 & ratio_vs_low < 1/3)]
  print(cmp[order(-ratio_vs_high)])
  n_flagged <- sum(cmp$flag_gt_3x, na.rm = TRUE)
  cat("\nConfigs flagged (>3x or <1/3x prediction):", n_flagged, "\n")
} else {
  message("No prediction file at ", pred_path, "; skipping comparison.")
  cmp <- NULL
}

# =============================================================================
# 5. Extrapolation to target scale (measured anchor) + full feasibility budget
# =============================================================================
cat("\n=== 5. Target-scale anchor (MEASURED) + full-budget extrapolation ===\n")
target_R <- 10000; target_chains <- 2

anchor <- res_ok[grepl("^target_anchor_", config_id)]
extrapolate_from_anchor <- function(model_name) {
  a <- anchor[model == model_name]
  if (nrow(a) == 0) return(NULL)
  cols <- if (model_name == "hmnl") {
    c("prof_t_cache_rebuild", "prof_t_beta_mh", "prof_t_delta_serial",
      "prof_t_hierarchy", "prof_t_recording")
  } else {
    c("prof_t_startup_gram", "prof_t_fused_pass", "prof_t_delta_parallel",
      "prof_t_rss", "prof_t_hierarchy", "prof_t_recording")
  }
  per_iter <- sapply(cols, function(bc) a[[bc]][1] / a$R[1])
  # startup_gram is one-time (not per iteration); exclude from the per-iter
  # sum, add back once at the end.
  onetime <- if ("prof_t_startup_gram" %in% cols) per_iter["prof_t_startup_gram"] * a$R[1] else 0
  per_iter_iter_only <- per_iter[setdiff(names(per_iter), "prof_t_startup_gram")]
  total_per_iter <- sum(per_iter_iter_only, na.rm = TRUE)
  per_chain_sec <- total_per_iter * target_R + onetime
  full_budget_sec <- per_chain_sec * target_chains
  data.table(model = model_name,
            anchor_N = a$N[1], anchor_J = a$J[1], anchor_R = a$R[1],
            anchor_peak_rss_bytes = a$peak_rss_bytes[1],
            sec_per_iter_at_anchor = total_per_iter,
            per_chain_sec = per_chain_sec, per_chain_min = per_chain_sec / 60,
            full_budget_sec = full_budget_sec, full_budget_min = full_budget_sec / 60)
}
extrap <- rbindlist(lapply(c("hmnl", "hmnp"), extrapolate_from_anchor), fill = TRUE)
print(extrap)

cat("\nNOTE: extrapolation is now MEASURED (target_anchor_hmnl/hmnp actually ran\n",
   "at N=30000, J=200) times target_R/anchor_R -- linear-in-R extrapolation\n",
   "of an already-measured per-iteration cost at the true target N and J, not\n",
   "a cross-scale power-law guess. This is the load-bearing number for the\n",
   "feasibility verdict.\n", sep = "")

# --- Save all analysis tables for the report -------------------------------
out_dir <- file.path(script_dir, "results")
fwrite(slopes, file.path(out_dir, "analysis_scaling_slopes.csv"))
fwrite(amdahl_dt, file.path(out_dir, "analysis_amdahl.csv"))
fwrite(mem_dt, file.path(out_dir, "analysis_memory_probe.csv"))
if (!is.null(cmp)) fwrite(cmp, file.path(out_dir, "analysis_predicted_vs_actual.csv"))
fwrite(extrap, file.path(out_dir, "analysis_target_extrapolation.csv"))
message("\nWrote analysis_*.csv to ", out_dir)
