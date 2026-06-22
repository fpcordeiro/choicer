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
source("_benchmarks/mnlogit/process.R")

config <- mnlogit_benchmark_config()
config <- bench_parse_cli(defaults = config)
config <- mnlogit_normalize_config(config)

output_dir <- bench_output_dir(root, "mnlogit", config$tag)
message("Writing benchmark outputs to: ", output_dir)
mnlogit_prepare_attempt_dirs(output_dir)

bench_load_choicer(root, clean_dll = config$clean_dll)
choicer_thread_state <- bench_set_choicer_threads(config$n_threads)
message(
  "choicer OpenMP threads: requested=", config$n_threads,
  ", reported_max=", choicer_thread_state$omp_get_max_threads
)

specs <- bench_spec_grid(config)
attempt_status <- mnlogit_attempt_grid(specs, config)
status_path <- file.path(output_dir, "run_status.csv")
mnlogit_write_run_status(attempt_status, status_path)

message("Benchmark attempts: ", nrow(attempt_status))
if (is.finite(config$max_run_sec)) {
  message("Per-attempt child-process timeout: ", config$max_run_sec, " seconds")
}

for (i in seq_len(nrow(attempt_status))) {
  attempt <- attempt_status[i]
  message(sprintf(
    "[mnlogit] attempt=%d/%d spec=%s run=%d package=%s",
    i, nrow(attempt_status), attempt$spec_id, attempt$run, attempt$package
  ))
  mnlogit_run_child_attempt(
    attempt = attempt,
    root = root,
    output_dir = output_dir,
    config = config,
    status = attempt_status,
    status_path = status_path
  )
  mnlogit_reconcile_partial_outputs(output_dir, status_path, config)
}

final_outputs <- mnlogit_reconcile_partial_outputs(output_dir, status_path, config)
if (is.null(final_outputs)) {
  stop("No benchmark results were produced.")
}
raw_results <- final_outputs$raw_results
coef_results <- final_outputs$coef_results
summary_results <- final_outputs$summary_results
raw_path <- final_outputs$raw_path
coef_path <- final_outputs$coef_path
summary_path <- final_outputs$summary_path
table_path <- file.path(output_dir, "summary_table.md")
session_path <- file.path(output_dir, "session_info.txt")

table_path <- benchmark_write_summary_table(summary_results, output_dir)
plot_paths <- benchmark_write_runtime_plot(summary_results, output_dir, metric = config$plot_metric)
bench_write_session_info(session_path)

metadata <- benchmark_capture_metadata(
  root = root,
  model = "mnlogit",
  config = config,
  packages = unique(c(config$packages, vapply(config$packages, mnlogit_dependency, character(1L)))),
  choicer_thread_state = choicer_thread_state,
  output_dir = output_dir,
  output_files = c(
    raw_path, coef_path, summary_path, status_path, table_path, session_path,
    plot_paths, file.path(output_dir, "partials"), file.path(output_dir, "logs")
  )
)
benchmark_write_metadata(metadata, output_dir)

message("Benchmark complete.")
message("Summary table: ", table_path)
if (length(plot_paths)) message("Runtime plot: ", plot_paths[[1L]])
