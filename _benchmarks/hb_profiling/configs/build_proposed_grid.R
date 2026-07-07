#!/usr/bin/env Rscript
# Builds the PROPOSED profiling grid (see request.md / harness-summary.md).
# NOT executed automatically -- writes configs/proposed_grid.csv for the user
# to review/approve before any full-grid run (Phase 1 hard scope limit: only
# the smoke_grid.csv configs may actually be executed in this phase).

suppressPackageStartupMessages(library(data.table))

mk <- function(config_id, model, N, J, K = 10L, P = 10L, R, burn, thin = 1L,
              threads = 5L, keep_beta_i = "means", seed = 1L) {
  data.table(config_id = config_id, model = model, N = N, J = J, K = K, P = P,
            R = R, burn = burn, thin = thin, threads = threads,
            keep_beta_i = keep_beta_i, seed = seed)
}

rows <- list()

# --- N sweep at J=50, K=10, P=10, R=200, burn=100, threads=5, both models ---
for (N in c(1000L, 2000L, 4000L, 8000L)) {
  for (model in c("hmnl", "hmnp")) {
    rows[[length(rows) + 1]] <- mk(
      sprintf("nsweep_%s_N%d_J50", model, N), model, N = N, J = 50L,
      R = 200L, burn = 100L
    )
  }
}

# --- J sweep at N=4000, K=10, P=10, R=200, burn=100, threads=5, both models -
for (J in c(25L, 50L, 100L, 200L)) {
  for (model in c("hmnl", "hmnp")) {
    rows[[length(rows) + 1]] <- mk(
      sprintf("jsweep_%s_N4000_J%d", model, J), model, N = 4000L, J = J,
      R = 200L, burn = 100L
    )
  }
}

# --- Thread points at N=4000, J=100: threads in {1, 5} ----------------------
for (th in c(1L, 5L)) {
  for (model in c("hmnl", "hmnp")) {
    rows[[length(rows) + 1]] <- mk(
      sprintf("threads_%s_t%d_N4000_J100", model, th), model, N = 4000L,
      J = 100L, R = 200L, burn = 100L, threads = th
    )
  }
}

# --- Target-scale anchor: N=30000, J=200, R=30, burn=10, both models --------
for (model in c("hmnl", "hmnp")) {
  rows[[length(rows) + 1]] <- mk(
    sprintf("target_anchor_%s", model), model, N = 30000L, J = 200L,
    R = 30L, burn = 10L
  )
}

# --- Memory probe: keep_beta_i = draws at N=8000, J=50, R=400, burn=200 -----
for (model in c("hmnl", "hmnp")) {
  rows[[length(rows) + 1]] <- mk(
    sprintf("mem_probe_%s", model), model, N = 8000L, J = 50L, R = 400L,
    burn = 200L, keep_beta_i = "draws"
  )
}

grid <- rbindlist(rows)
out_path <- file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "proposed_grid.csv")
if (length(out_path) == 0 || is.na(out_path) || out_path == "") {
  out_path <- "proposed_grid.csv"
}
fwrite(grid, out_path)
message("Wrote ", nrow(grid), " HMNL/HMNP config(s) to ", out_path)

# --- run_mnprobit contrast: separate grid (wall-time only, not profiled, not
# routed through 10_run_config.R's HMNL/HMNP path) --------------------------
mnp_grid <- data.table(
  config_id = sprintf("mnprobit_contrast_J%d", c(4L, 8L, 12L)),
  model = "mnprobit", N = 2000L, J = c(4L, 8L, 12L), K = 10L, P = NA_integer_,
  R = 200L, burn = 40L, thin = 1L, threads = 5L, keep_beta_i = "none",
  seed = 1L
)
mnp_out <- file.path(dirname(out_path), "proposed_grid_mnprobit.csv")
fwrite(mnp_grid, mnp_out)
message("Wrote ", nrow(mnp_grid), " run_mnprobit contrast config(s) to ", mnp_out)
