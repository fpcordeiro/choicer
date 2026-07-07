#!/usr/bin/env Rscript
# Config-driven harness: for each row of a config grid, launch
# `_benchmarks/hb_profiling/10_run_config.R` as a FRESH subprocess wrapped in
# `/usr/bin/time -l` (macOS: reports "maximum resident set size" in BYTES),
# parse the worker's JSON result and the peak-RSS line, and append one row to
# the results CSV.
#
# Usage:
#   Rscript _benchmarks/hb_profiling/20_harness.R \
#     --grid _benchmarks/hb_profiling/configs/smoke_grid.csv \
#     --lib-dir /path/to/templib \
#     --results _benchmarks/hb_profiling/results/hb_profiling_results.csv
#
# The grid CSV has one row per config with columns matching 10_run_config.R's
# key=value args: config_id,model,N,J,K,P,R,burn,thin,threads,keep_beta_i,seed
#
# Fixed seeds throughout (the `seed` column); nothing here draws from R's
# global RNG, so re-running the same grid is deterministic (the kernels'
# own draws are additionally deterministic per hb_internal.h's per-
# (iteration, unit) RNG-stream design).

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}

this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file) == 1) dirname(this_file) else "."

grid_path    <- get_arg("--grid", file.path(script_dir, "configs", "smoke_grid.csv"))
lib_dir      <- get_arg("--lib-dir",
  "/Users/fernando/.claude/plugins/data/statsclaw-statsclaw/workspace/choicer/tmp/hb_prof_build/templib")
results_path <- get_arg("--results",
  file.path(script_dir, "results", "hb_profiling_results.csv"))
worker_path  <- file.path(script_dir, "10_run_config.R")
tmp_dir      <- file.path(script_dir, "results", "tmp")

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

grid <- fread(grid_path)
message("Loaded ", nrow(grid), " config(s) from ", grid_path)

parse_max_rss <- function(time_stderr_lines) {
  # macOS `/usr/bin/time -l` line: "  1234567  maximum resident set size"
  ln <- grep("maximum resident set size", time_stderr_lines, value = TRUE)
  if (length(ln) == 0) return(NA_real_)
  as.numeric(trimws(sub("maximum resident set size.*$", "", ln[1])))
}

run_one <- function(cfg) {
  config_id <- as.character(cfg$config_id)
  out_json <- file.path(tmp_dir, paste0(config_id, ".json"))
  err_log  <- file.path(tmp_dir, paste0(config_id, ".time.err"))
  if (file.exists(out_json)) file.remove(out_json)

  worker_args <- c(
    sprintf("model=%s", cfg$model),
    sprintf("N=%d", cfg$N), sprintf("J=%d", cfg$J),
    sprintf("K=%d", cfg$K), sprintf("P=%d", cfg$P),
    sprintf("R=%d", cfg$R), sprintf("burn=%d", cfg$burn),
    sprintf("thin=%d", cfg$thin %||% 1L),
    sprintf("threads=%d", cfg$threads),
    sprintf("keep_beta_i=%s", cfg$keep_beta_i %||% "means"),
    sprintf("seed=%d", cfg$seed %||% 1L),
    sprintf("lib_dir=%s", lib_dir),
    sprintf("config_id=%s", config_id),
    sprintf("out=%s", out_json)
  )
  cmd <- sprintf(
    "/usr/bin/time -l Rscript %s %s > %s 2>&1",
    shQuote(worker_path), paste(shQuote(worker_args), collapse = " "),
    shQuote(err_log)
  )
  message("Running: ", config_id, " (", cfg$model, ", N=", cfg$N,
          ", J=", cfg$J, ", R=", cfg$R, ")")
  t0 <- Sys.time()
  status <- system(cmd)
  t1 <- Sys.time()
  wall_outer <- as.numeric(difftime(t1, t0, units = "secs"))

  log_lines <- readLines(err_log, warn = FALSE)
  peak_rss <- parse_max_rss(log_lines)

  if (status != 0 || !file.exists(out_json)) {
    warning("Config ", config_id, " FAILED (status ", status, "); see ",
            err_log)
    return(data.table(
      config_id = config_id, model = cfg$model, N = cfg$N, J = cfg$J,
      K = cfg$K, P = cfg$P, R = cfg$R, burn = cfg$burn, threads = cfg$threads,
      keep_beta_i = cfg$keep_beta_i %||% "means", status = "FAILED",
      wall_outer_sec = wall_outer, peak_rss_bytes = peak_rss
    ))
  }

  res <- jsonlite::fromJSON(out_json)
  prof <- res$profile
  row <- data.table(
    config_id = config_id, model = cfg$model, N = cfg$N, J = cfg$J,
    K = cfg$K, P = cfg$P, R = cfg$R, burn = cfg$burn, threads = cfg$threads,
    keep_beta_i = cfg$keep_beta_i %||% "means", status = "OK",
    t_setup_sec = res$t_setup_sec,
    wall_outer_sec = wall_outer,
    peak_rss_bytes = peak_rss
  )
  if (!is.null(prof)) {
    for (nm in names(prof)) {
      if (is.numeric(prof[[nm]]) && length(prof[[nm]]) == 1) {
        row[[paste0("prof_", nm)]] <- prof[[nm]]
      }
    }
    if (!is.null(prof$t_total) && cfg$R > 0) {
      row[["sec_per_iter"]] <- prof$t_total / cfg$R
    }
  }
  row
}

all_rows <- rbindlist(lapply(seq_len(nrow(grid)), function(i) run_one(grid[i])),
                      fill = TRUE)

if (file.exists(results_path)) {
  prev <- fread(results_path)
  all_rows <- rbindlist(list(prev, all_rows), fill = TRUE)
}
fwrite(all_rows, results_path)
message("Wrote ", nrow(all_rows), " total row(s) to ", results_path)
print(all_rows)
