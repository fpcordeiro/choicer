#!/usr/bin/env Rscript
# Predicts per-config wall time for configs/proposed_grid*.csv by extrapolating
# from the Phase-1 smoke measurements (N=500, J=20, K=10, P=10, threads=5,
# R=100) in results/hb_profiling_results.csv, using block-level scaling laws
# read directly off the kernel code (src/hmnlogit.cpp, src/hmnprobit.cpp):
#
#   cache_rebuild      (HMNL) : O(N*J) parallel   -> scales with N*J
#   beta_mh            (HMNL) : per-respondent cost is O(J) [task loglik] +
#                                O(K^3) [Cholesky, K fixed at 10] -- AMBIGUOUS
#                                which term dominates as J grows; reported as
#                                a [low, high] = [scales-with-N, scales-with-N*J]
#                                range. This ambiguity is exactly why a real
#                                J-sweep (not just this one smoke point) is
#                                needed -- flagged in harness-summary.md.
#   delta_serial       (HMNL) : O(N*J), single thread            -> N*J
#   hierarchy       (HMNL/HMNP): O(N*K^2) dominates O(J*P^2) at target scale
#                                (N*K^2 = 3e6 vs J*P^2 = 2e4 at N=30000,J=200)
#                                -> scales with N
#   recording       (HMNL/HMNP): O(J) per recorded iter, tiny magnitude -> N/A
#   fused_pass         (HMNP) : split via its OWN thread-summed sub-blocks
#                                (this is exactly why that instrumentation was
#                                added): tn_sweep O(N*J), beta_draw O(N*K^2)
#                                [K fixed], mu_x_refresh O(N*J); each scaled
#                                independently, then divided by thread count.
#   delta_parallel     (HMNP) : O(N*J) work spread over J-parallel loop -> N*J
#   rss                (HMNP) : O(N*J) parallel over respondents       -> N*J
#   startup_gram       (HMNP) : ONE-TIME (not per iteration); O(N*J), parallel
#
# All "parallel, work scales with X" blocks are assumed thread-count-invariant
# in TOTAL work; wall time is assumed to scale as (work) / threads when
# threads changes between the smoke baseline (5) and a grid config.
#
# This is a PRE-GRID analytical estimate from a SINGLE (N, J) data point --
# it cannot distinguish an N-only law from an N*J law with certainty (see the
# beta_mh range above). The actual N-sweep / J-sweep grid run (pending
# approval) will fit the true exponent empirically via 30_analyze.R.

suppressPackageStartupMessages(library(data.table))

this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file) == 1) dirname(this_file) else "."
root <- normalizePath(file.path(script_dir, ".."))

smoke <- fread(file.path(root, "results", "hb_profiling_results.csv"))
smoke <- smoke[status == "OK"]
N0 <- 500; J0 <- 20; threads0 <- 5

hmnl0 <- smoke[model == "hmnl"][1]
hmnp0 <- smoke[model == "hmnp"][1]

predict_hmnl <- function(N, J, R, threads) {
  rNJ <- (N * J) / (N0 * J0)
  rN  <- N / N0
  thr_ratio <- threads0 / threads   # wall time ~ work / threads

  cache_rebuild <- hmnl0$prof_t_cache_rebuild / hmnl0$prof_R * rNJ * thr_ratio
  beta_low  <- hmnl0$prof_t_beta_mh / hmnl0$prof_R * rN  * thr_ratio
  beta_high <- hmnl0$prof_t_beta_mh / hmnl0$prof_R * rNJ * thr_ratio
  delta_serial <- hmnl0$prof_t_delta_serial / hmnl0$prof_R * rNJ  # single thread, no thr_ratio
  hierarchy <- hmnl0$prof_t_hierarchy / hmnl0$prof_R * rN          # master-only, no thr_ratio
  recording <- hmnl0$prof_t_recording / hmnl0$prof_R * rN

  per_iter_low  <- cache_rebuild + beta_low  + delta_serial + hierarchy + recording
  per_iter_high <- cache_rebuild + beta_high + delta_serial + hierarchy + recording
  list(per_iter_low = per_iter_low, per_iter_high = per_iter_high,
      total_low = per_iter_low * R, total_high = per_iter_high * R)
}

predict_hmnp <- function(N, J, R, threads) {
  rNJ <- (N * J) / (N0 * J0)
  rN  <- N / N0
  thr_ratio <- threads0 / threads

  tn_per_iter <- hmnp0$prof_fused_tn_sweep_cpu_sum / hmnp0$prof_R
  bd_per_iter <- hmnp0$prof_fused_beta_draw_cpu_sum / hmnp0$prof_R
  mx_per_iter <- hmnp0$prof_fused_mu_x_refresh_cpu_sum / hmnp0$prof_R
  fused <- (tn_per_iter * rNJ + bd_per_iter * rN + mx_per_iter * rNJ) / threads

  delta_parallel <- hmnp0$prof_t_delta_parallel / hmnp0$prof_R * rNJ * thr_ratio
  rss <- hmnp0$prof_t_rss / hmnp0$prof_R * rNJ * thr_ratio
  hierarchy <- hmnp0$prof_t_hierarchy / hmnp0$prof_R * rN
  recording <- hmnp0$prof_t_recording / hmnp0$prof_R * rN
  startup_once <- hmnp0$prof_t_startup_gram * rNJ * thr_ratio

  per_iter <- fused + delta_parallel + rss + hierarchy + recording
  list(per_iter_low = per_iter, per_iter_high = per_iter,
      total_low = per_iter * R + startup_once,
      total_high = per_iter * R + startup_once)
}

grid <- fread(file.path(script_dir, "proposed_grid.csv"))
grid[, `:=`(pred_total_low_sec = NA_real_, pred_total_high_sec = NA_real_)]
for (i in seq_len(nrow(grid))) {
  cfg <- grid[i]
  pred <- if (cfg$model == "hmnl") {
    predict_hmnl(cfg$N, cfg$J, cfg$R, cfg$threads)
  } else {
    predict_hmnp(cfg$N, cfg$J, cfg$R, cfg$threads)
  }
  grid[i, pred_total_low_sec := pred$total_low]
  grid[i, pred_total_high_sec := pred$total_high]
}
grid[, pred_total_low_min := pred_total_low_sec / 60]
grid[, pred_total_high_min := pred_total_high_sec / 60]
grid[, flag_over_10min := pred_total_high_min > 10]

out_path <- file.path(script_dir, "predicted_grid_times.csv")
fwrite(grid, out_path)
cat("Wrote predictions to", out_path, "\n\n")
print(grid[, .(config_id, model, N, J, R, threads,
              pred_total_low_min = round(pred_total_low_min, 3),
              pred_total_high_min = round(pred_total_high_min, 3),
              flag_over_10min)])

# --- Full feasibility-budget extrapolation (R=10000, burn=5000, 2 SEQUENTIAL
# chains -- run_h{mnl,mnp}ogit's chains are NOT parallel, see
# R/hmnlogit_utils.R / hmnprobit_utils.R's lapply-based extra_chains loop) --
target <- predict_hmnl(30000, 200, 10000, 5)
cat("\nHMNL full budget (N=30000, J=200, R=10000, 2 sequential chains):\n")
cat(sprintf("  per-chain: [%.1f, %.1f] min | both chains: [%.1f, %.1f] min\n",
           target$total_low / 60, target$total_high / 60,
           2 * target$total_low / 60, 2 * target$total_high / 60))

target_p <- predict_hmnp(30000, 200, 10000, 5)
cat("\nHMNP full budget (N=30000, J=200, R=10000, 2 sequential chains):\n")
cat(sprintf("  per-chain: %.1f min | both chains: %.1f min\n",
           target_p$total_low / 60, 2 * target_p$total_low / 60))
