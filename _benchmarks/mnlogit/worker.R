#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
request_arg <- grep("^--request=", args, value = TRUE)
if (length(request_arg) != 1L) {
  stop("Worker requires exactly one --request=<path> argument.")
}
request_path <- sub("^--request=", "", request_arg[[1L]])
request <- readRDS(request_path)

root <- request$root
setwd(root)

suppressPackageStartupMessages({
  library(data.table)
})

source("_benchmarks/benchmark_run.R")
source("_benchmarks/mnlogit/config.R")
source("_benchmarks/mnlogit/dgp.R")
source("_benchmarks/mnlogit/dispatchers.R")

attempt <- request$attempt
config <- request$config

if (identical(attempt$package, "choicer")) {
  bench_load_choicer(root, clean_dll = FALSE)
  invisible(bench_set_choicer_threads(config$n_threads))
}

spec <- data.table::data.table(
  sweep = attempt$sweep,
  dimension = attempt$dimension,
  dimension_value = attempt$dimension_value,
  spec_id = attempt$spec_id,
  N = as.integer(attempt$N),
  J = as.integer(attempt$J)
)

sim <- simulate_mnlogit_benchmark_data(
  N = spec$N,
  J = spec$J,
  beta = config$beta,
  seed = as.integer(attempt$seed)
)

res <- mnlogit_run_package(
  package = attempt$package,
  dt = sim$data,
  x_cols = sim$x_cols,
  spec = spec,
  run = as.integer(attempt$run),
  config = config
)

bench_atomic_fwrite(res$raw, request$paths$raw)
if (nrow(res$coef)) {
  bench_atomic_fwrite(res$coef, request$paths$coef)
}
