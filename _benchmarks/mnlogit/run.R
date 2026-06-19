#!/usr/bin/env Rscript

root <- local({
  path <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) break
    parent <- dirname(path)
    if (identical(parent, path)) stop("Run this script from inside the choicer repository.")
    path <- parent
  }
  path
})
setwd(root)

suppressPackageStartupMessages({
  library(data.table)
})

source("_benchmarks/benchmark_run.R")
source("_benchmarks/benchmark_metadata.R")
source("_benchmarks/benchmark_metrics.R")
source("_benchmarks/benchmark_plotting.R")
source("_benchmarks/benchmark_tables.R")
source("_benchmarks/mnlogit/config.R")
source("_benchmarks/mnlogit/dgp.R")
source("_benchmarks/mnlogit/dispatchers.R")

config <- mnlogit_benchmark_config()
config <- bench_parse_cli(defaults = config)
config <- mnlogit_normalize_config(config)

output_dir <- bench_output_dir(root, "mnlogit", config$tag)
message("Writing benchmark outputs to: ", output_dir)

bench_load_choicer(root, clean_dll = config$clean_dll)
if (exists("set_num_threads", mode = "function")) {
  set_num_threads(config$n_threads)
}

specs <- bench_spec_grid(config)
raw_rows <- vector("list", nrow(specs) * config$n_runs * length(config$packages))
coef_rows <- vector("list", length(raw_rows))
idx <- 1L

for (i in seq_len(nrow(specs))) {
  spec <- specs[i]
  for (run in seq_len(config$n_runs)) {
    seed <- config$seed + i * 10000L + run
    sim <- simulate_mnlogit_benchmark_data(
      N = spec$N,
      J = spec$J,
      beta = config$beta,
      seed = seed
    )
    for (pkg in config$packages) {
      message(sprintf(
        "[mnlogit] spec=%s run=%d package=%s",
        spec$spec_id, run, pkg
      ))
      res <- mnlogit_run_package(pkg, sim$data, sim$x_cols, spec, run, config)
      raw_rows[[idx]] <- res$raw
      coef_rows[[idx]] <- res$coef
      idx <- idx + 1L
    }
  }
}

raw_results <- data.table::rbindlist(raw_rows, use.names = TRUE, fill = TRUE)
coef_results <- data.table::rbindlist(coef_rows, use.names = TRUE, fill = TRUE)
coef_results <- benchmark_compute_coef_checks(coef_results, raw_results, config)
summary_results <- benchmark_summarise_results(raw_results, coef_results, config)

raw_path <- file.path(output_dir, "raw_results.csv")
coef_path <- file.path(output_dir, "coef_results.csv")
summary_path <- file.path(output_dir, "summary_results.csv")
table_path <- file.path(output_dir, "summary_table.md")
session_path <- file.path(output_dir, "session_info.txt")

data.table::fwrite(raw_results, raw_path)
data.table::fwrite(coef_results, coef_path)
data.table::fwrite(summary_results, summary_path)
table_path <- benchmark_write_summary_table(summary_results, output_dir)
plot_paths <- benchmark_write_runtime_plot(summary_results, output_dir, metric = config$plot_metric)
bench_write_session_info(session_path)

metadata <- benchmark_capture_metadata(
  root = root,
  model = "mnlogit",
  config = config,
  packages = unique(c(config$packages, vapply(config$packages, mnlogit_dependency, character(1L)))),
  output_dir = output_dir,
  output_files = c(raw_path, coef_path, summary_path, table_path, session_path, plot_paths)
)
benchmark_write_metadata(metadata, output_dir)

message("Benchmark complete.")
message("Summary table: ", table_path)
if (length(plot_paths)) message("Runtime plot: ", plot_paths[[1L]])
