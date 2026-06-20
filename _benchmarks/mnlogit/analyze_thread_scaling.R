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

defaults <- list(
  input_dir = file.path(root, "_benchmarks", "mnlogit", "output"),
  output_dir = NULL,
  pattern = "choicer_threads_",
  package = "choicer",
  metric = "mean_fit_time_sec",
  metrics = c("mean_fit_time_sec", "mean_total_time_sec", "mean_vcov_time_sec"),
  require_verified_threads = TRUE,
  within_best = 1.10,
  tag = "thread_scaling_analysis"
)

config <- bench_parse_cli(defaults = defaults)
config$input_dir <- normalizePath(config$input_dir[[1L]], winslash = "/", mustWork = TRUE)
config$pattern <- config$pattern[[1L]]
config$package <- config$package[[1L]]
config$metric <- config$metric[[1L]]
config$metrics <- unique(c(config$metrics, config$metric))
config$require_verified_threads <- isTRUE(config$require_verified_threads[[1L]])
config$within_best <- as.numeric(config$within_best[[1L]])
config$tag <- config$tag[[1L]]

if (is.null(config$output_dir)) {
  output_dir <- bench_output_dir(root, "mnlogit", config$tag)
} else {
  output_dir <- normalizePath(config$output_dir[[1L]], winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
}

message("Writing thread-scaling analysis to: ", output_dir)

parse_thread_from_text <- function(x, pattern = "choicer_threads_") {
  if (is.null(x) || !length(x) || is.na(x)) return(NA_integer_)
  rx <- paste0(".*", gsub("([\\W])", "\\\\\\1", pattern), "([0-9]+).*")
  if (!grepl(rx, x)) return(NA_integer_)
  as.integer(sub(rx, "\\1", x))
}

read_json_metadata <- function(path) {
  if (!file.exists(path) || !requireNamespace("jsonlite", quietly = TRUE)) {
    return(list())
  }
  tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) list()
  )
}

read_run_metadata <- function(run_dir) {
  meta <- read_json_metadata(file.path(run_dir, "metadata.json"))
  run_tag <- meta$run_tag %||% basename(run_dir)
  n_threads <- suppressWarnings(as.integer(meta$config$n_threads %||% NA_integer_))
  if (is.na(n_threads)) {
    n_threads <- parse_thread_from_text(run_tag, config$pattern)
  }
  if (is.na(n_threads)) {
    n_threads <- parse_thread_from_text(basename(run_dir), config$pattern)
  }
  reported_threads <- suppressWarnings(as.integer(
    meta$choicer_threads$omp_get_max_threads %||% NA_integer_
  ))
  openmp_enabled <- meta$choicer_threads$openmp_enabled %||% NA
  thread_state_verified <- isTRUE(openmp_enabled) &&
    !is.na(reported_threads) &&
    !is.na(n_threads) &&
    identical(as.integer(reported_threads), as.integer(n_threads))

  data.table(
    source_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    run_folder = basename(run_dir),
    run_tag = as.character(run_tag %||% NA_character_),
    n_threads = n_threads,
    reported_threads = reported_threads,
    openmp_enabled = openmp_enabled,
    thread_state_verified = thread_state_verified,
    timestamp = as.character(meta$timestamp %||% NA_character_),
    command = as.character(meta$command %||% NA_character_)
  )
}

run_dirs <- list.dirs(config$input_dir, full.names = TRUE, recursive = FALSE)
metadata <- rbindlist(lapply(run_dirs, read_run_metadata), use.names = TRUE, fill = TRUE)
metadata <- metadata[
  grepl(config$pattern, run_tag, fixed = TRUE) |
    grepl(config$pattern, run_folder, fixed = TRUE)
]
metadata <- metadata[!is.na(n_threads)]

if (!nrow(metadata)) {
  stop("No thread benchmark folders matched pattern '", config$pattern, "' in ", config$input_dir, ".")
}
if (isTRUE(config$require_verified_threads)) {
  verified_idx <- metadata$thread_state_verified %in% TRUE
  unverified <- metadata[!verified_idx]
  metadata <- metadata[verified_idx]
  if (nrow(unverified)) {
    message(
      "Ignoring ", nrow(unverified),
      " unverified thread run folder(s). Re-run benchmarks with the fixed run.R ",
      "or pass --require-verified-threads=false to include them."
    )
  }
  if (!nrow(metadata)) {
    stop(
      "No matched run folders have verified OpenMP thread metadata. ",
      "Re-run the thread benchmark with the fixed _benchmarks/mnlogit/run.R."
    )
  }
}

read_raw_results <- function(meta_row) {
  path <- file.path(meta_row$source_dir, "raw_results.csv")
  if (!file.exists(path)) return(data.table())
  dt <- data.table::fread(path)
  dt[, `:=`(
    source_dir = meta_row$source_dir,
    run_folder = meta_row$run_folder,
    run_tag = meta_row$run_tag,
    n_threads = as.integer(meta_row$n_threads)
  )]
  if (!"vcov_time_sec" %in% names(dt)) dt[, vcov_time_sec := NA_real_]
  if (!"vcov_timing_status" %in% names(dt)) dt[, vcov_timing_status := NA_character_]
  if (!"requested_compute_vcov" %in% names(dt)) dt[, requested_compute_vcov := NA]
  if (!"compute_vcov" %in% names(dt)) dt[, compute_vcov := NA]
  dt[]
}

read_summary_results <- function(meta_row) {
  path <- file.path(meta_row$source_dir, "summary_results.csv")
  if (!file.exists(path)) return(data.table())
  dt <- data.table::fread(path)
  dt[, `:=`(
    source_dir = meta_row$source_dir,
    run_folder = meta_row$run_folder,
    run_tag = meta_row$run_tag,
    n_threads = as.integer(meta_row$n_threads)
  )]
  if (!"mean_vcov_time_sec" %in% names(dt)) dt[, mean_vcov_time_sec := NA_real_]
  if (!"sd_vcov_time_sec" %in% names(dt)) dt[, sd_vcov_time_sec := NA_real_]
  if (!"median_vcov_time_sec" %in% names(dt)) dt[, median_vcov_time_sec := NA_real_]
  if (!"vcov_timing_status" %in% names(dt)) dt[, vcov_timing_status := NA_character_]
  dt[]
}

raw_results <- rbindlist(
  lapply(seq_len(nrow(metadata)), function(i) read_raw_results(metadata[i])),
  use.names = TRUE,
  fill = TRUE
)

summary_inputs <- rbindlist(
  lapply(seq_len(nrow(metadata)), function(i) read_summary_results(metadata[i])),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(raw_results)) {
  raw_results <- raw_results[package == config$package]
  raw_results[, ok := status == "ok"]
  thread_summary <- raw_results[, .(
    n_attempts = .N,
    n_ok = sum(ok, na.rm = TRUE),
    success_rate = mean(ok, na.rm = TRUE),
    convergence_rate = mean(converged %in% TRUE, na.rm = TRUE),
    mean_setup_time_sec = mean(setup_time_sec, na.rm = TRUE),
    sd_setup_time_sec = stats::sd(setup_time_sec, na.rm = TRUE),
    median_setup_time_sec = stats::median(setup_time_sec, na.rm = TRUE),
    mean_fit_time_sec = mean(fit_time_sec, na.rm = TRUE),
    sd_fit_time_sec = stats::sd(fit_time_sec, na.rm = TRUE),
    median_fit_time_sec = stats::median(fit_time_sec, na.rm = TRUE),
    mean_vcov_time_sec = mean(vcov_time_sec, na.rm = TRUE),
    sd_vcov_time_sec = stats::sd(vcov_time_sec, na.rm = TRUE),
    median_vcov_time_sec = stats::median(vcov_time_sec, na.rm = TRUE),
    mean_total_time_sec = mean(total_time_sec, na.rm = TRUE),
    sd_total_time_sec = stats::sd(total_time_sec, na.rm = TRUE),
    median_total_time_sec = stats::median(total_time_sec, na.rm = TRUE),
    vcov_timing_status = paste(sort(unique(stats::na.omit(vcov_timing_status))), collapse = ",")
  ), by = .(n_threads, sweep, dimension, dimension_value, spec_id, N, J, package)]
} else if (nrow(summary_inputs)) {
  thread_summary <- summary_inputs[package == config$package]
  if (!"n_attempts" %in% names(thread_summary)) {
    thread_summary[, n_attempts := n_runs %||% n_ok]
  }
} else {
  stop("No raw_results.csv or summary_results.csv files were found in matched run folders.")
}

numeric_cols <- names(thread_summary)[vapply(thread_summary, is.numeric, logical(1L))]
for (col in numeric_cols) {
  data.table::set(thread_summary, which(is.nan(thread_summary[[col]])), col, NA_real_)
}

thread_summary[, spec_label := sprintf("%s: N=%s, J=%s", sweep, N, J)]
setorder(thread_summary, sweep, dimension_value, n_threads)

available_metrics <- config$metrics[
  config$metrics %in% names(thread_summary) &
    vapply(config$metrics, function(m) any(is.finite(thread_summary[[m]])), logical(1L))
]
if (!config$metric %in% available_metrics) {
  stop("Primary metric '", config$metric, "' is not available in the collected results.")
}

analyze_metric <- function(dt, metric, within_best) {
  metric_name <- metric
  work <- copy(dt[is.finite(get(metric)) & n_ok > 0])
  if (!nrow(work)) return(list(by_thread = data.table(), by_spec = data.table()))

  spec_cols <- c("sweep", "dimension", "dimension_value", "spec_id", "N", "J", "package")
  min_thread <- min(work$n_threads, na.rm = TRUE)
  baseline <- work[n_threads == min_thread, c(spec_cols, metric), with = FALSE]
  setnames(baseline, metric, "baseline_time_sec")

  best <- work[, {
    idx <- which.min(get(metric))
    .(
      best_thread = n_threads[idx],
      best_time_sec = get(metric)[idx]
    )
  }, by = spec_cols]

  out <- merge(work, best, by = spec_cols, all.x = TRUE, sort = FALSE)
  out <- merge(out, baseline, by = spec_cols, all.x = TRUE, sort = FALSE)
  out[, `:=`(
    metric = metric,
    time_sec = get(metric),
    rel_to_best = get(metric) / best_time_sec,
    speedup_vs_min_thread = baseline_time_sec / get(metric),
    efficiency_vs_min_thread = (baseline_time_sec / get(metric)) / (n_threads / min_thread)
  )]

  by_spec <- out[rel_to_best <= within_best, {
    idx <- which.min(n_threads)
    .(
      metric = metric_name,
      recommended_threads = n_threads[idx],
      recommended_time_sec = time_sec[idx],
      best_thread = best_thread[idx],
      best_time_sec = best_time_sec[idx],
      rel_to_best = rel_to_best[idx],
      speedup_vs_min_thread = speedup_vs_min_thread[idx]
    )
  }, by = spec_cols]

  by_thread <- out[, .(
    metric = metric_name,
    n_specs = .N,
    n_best = sum(n_threads == best_thread, na.rm = TRUE),
    geo_mean_rel_to_best = exp(mean(log(rel_to_best), na.rm = TRUE)),
    median_rel_to_best = stats::median(rel_to_best, na.rm = TRUE),
    p90_rel_to_best = as.numeric(stats::quantile(rel_to_best, 0.90, na.rm = TRUE, names = FALSE)),
    geo_mean_speedup_vs_min_thread = exp(mean(log(speedup_vs_min_thread), na.rm = TRUE)),
    median_speedup_vs_min_thread = stats::median(speedup_vs_min_thread, na.rm = TRUE),
    mean_efficiency_vs_min_thread = mean(efficiency_vs_min_thread, na.rm = TRUE),
    min_success_rate = min(success_rate, na.rm = TRUE)
  ), by = .(n_threads)]
  setorder(by_thread, n_threads)

  list(detail = out, by_thread = by_thread, by_spec = by_spec)
}

metric_results <- lapply(available_metrics, function(m) {
  res <- analyze_metric(thread_summary, m, config$within_best)
  res$detail[, metric := m]
  res$by_thread[, metric := m]
  res$by_spec[, metric := m]
  res
})
names(metric_results) <- available_metrics

thread_detail <- rbindlist(lapply(metric_results, `[[`, "detail"), use.names = TRUE, fill = TRUE)
thread_overall <- rbindlist(lapply(metric_results, `[[`, "by_thread"), use.names = TRUE, fill = TRUE)
spec_recommendations <- rbindlist(lapply(metric_results, `[[`, "by_spec"), use.names = TRUE, fill = TRUE)

choose_overall <- function(overall_dt, metric, within_best) {
  metric_value <- metric
  dt <- copy(overall_dt[overall_dt$metric == metric_value & is.finite(geo_mean_rel_to_best)])
  if (!nrow(dt)) return(data.table())
  fastest <- dt[which.min(geo_mean_rel_to_best)]
  eligible <- dt[geo_mean_rel_to_best <= within_best]
  recommended <- if (nrow(eligible)) eligible[which.min(n_threads)] else fastest
  recommended[, `:=`(
    fastest_threads = fastest$n_threads,
    fastest_geo_mean_rel_to_best = fastest$geo_mean_rel_to_best,
    within_best = within_best
  )]
  recommended[]
}

overall_recommendation <- rbindlist(
  lapply(available_metrics, function(m) choose_overall(thread_overall, m, config$within_best)),
  use.names = TRUE,
  fill = TRUE
)

sweep_recommendation <- thread_detail[, .(
  geo_mean_rel_to_best = exp(mean(log(rel_to_best), na.rm = TRUE)),
  median_rel_to_best = stats::median(rel_to_best, na.rm = TRUE),
  geo_mean_speedup_vs_min_thread = exp(mean(log(speedup_vs_min_thread), na.rm = TRUE)),
  min_success_rate = min(success_rate, na.rm = TRUE)
), by = .(metric, sweep, n_threads)]
sweep_recommendation <- sweep_recommendation[, {
  fastest <- .SD[which.min(geo_mean_rel_to_best)]
  eligible <- .SD[geo_mean_rel_to_best <= config$within_best]
  rec <- if (nrow(eligible)) eligible[which.min(n_threads)] else fastest
  rec[, `:=`(
    fastest_threads = fastest$n_threads,
    fastest_geo_mean_rel_to_best = fastest$geo_mean_rel_to_best
  )]
  rec
}, by = .(metric, sweep)]
setorder(sweep_recommendation, metric, sweep)

write_plot <- function(detail, metric, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(character())
  metric_value <- metric
  dt <- copy(detail[detail$metric == metric_value & is.finite(time_sec)])
  if (!nrow(dt)) return(character())
  dt[, n_threads_factor := factor(n_threads, levels = sort(unique(n_threads)))]

  p_time <- ggplot2::ggplot(
    dt,
    ggplot2::aes(x = n_threads, y = time_sec, color = spec_label, group = spec_id)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ sweep, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(dt$n_threads))) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Threads",
      y = paste0(metric, " (seconds, log10 scale)"),
      color = "Spec"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())

  p_speedup <- ggplot2::ggplot(
    dt,
    ggplot2::aes(x = n_threads, y = speedup_vs_min_thread, color = spec_label, group = spec_id)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ sweep, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(dt$n_threads))) +
    ggplot2::labs(
      x = "Threads",
      y = "Speedup vs lowest tested thread count",
      color = "Spec"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())

  stub <- gsub("[^A-Za-z0-9_.-]", "_", metric)
  time_png <- file.path(output_dir, paste0("thread_scaling_", stub, ".png"))
  time_pdf <- file.path(output_dir, paste0("thread_scaling_", stub, ".pdf"))
  speed_png <- file.path(output_dir, paste0("thread_speedup_", stub, ".png"))
  speed_pdf <- file.path(output_dir, paste0("thread_speedup_", stub, ".pdf"))
  ggplot2::ggsave(time_png, p_time, width = 9, height = 6, dpi = 160)
  ggplot2::ggsave(time_pdf, p_time, width = 9, height = 6)
  ggplot2::ggsave(speed_png, p_speedup, width = 9, height = 6, dpi = 160)
  ggplot2::ggsave(speed_pdf, p_speedup, width = 9, height = 6)
  c(time_png, time_pdf, speed_png, speed_pdf)
}

plot_paths <- unlist(lapply(available_metrics, function(m) write_plot(thread_detail, m, output_dir)))

metadata_path <- file.path(output_dir, "thread_scaling_run_folders.csv")
summary_path <- file.path(output_dir, "thread_scaling_summary.csv")
detail_path <- file.path(output_dir, "thread_scaling_detail.csv")
overall_path <- file.path(output_dir, "thread_scaling_overall.csv")
spec_path <- file.path(output_dir, "thread_scaling_spec_recommendations.csv")
verdict_path <- file.path(output_dir, "thread_scaling_verdict.md")

fwrite(metadata, metadata_path)
fwrite(thread_summary, summary_path)
fwrite(thread_detail, detail_path)
fwrite(thread_overall, overall_path)
fwrite(spec_recommendations, spec_path)

primary_rec <- overall_recommendation[metric == config$metric]
primary_sweeps <- sweep_recommendation[metric == config$metric]

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_threads <- function(x) {
  sprintf("%d %s", x, ifelse(x == 1L, "thread", "threads"))
}

verdict_lines <- c(
  "# Choicer Thread-Scaling Verdict",
  "",
  paste("Input pattern:", config$pattern),
  paste("Output directory:", output_dir),
  paste("Primary metric:", config$metric),
  paste("Practical threshold:", sprintf("smallest thread count within %.1f%% of the per-spec best geometric mean", (config$within_best - 1) * 100)),
  "",
  "## Verdict",
  ""
)

if (nrow(primary_rec)) {
  verdict_lines <- c(
    verdict_lines,
    sprintf(
      "Use **%s** as the practical default for this benchmark on this machine.",
      fmt_threads(primary_rec$n_threads)
    ),
    "",
    sprintf(
      "The fastest geometric-mean result used **%s**. The recommended setting has geometric mean relative runtime %.3fx the spec-wise best and geometric mean speedup %.3fx versus the lowest tested thread count.",
      fmt_threads(primary_rec$fastest_threads),
      primary_rec$geo_mean_rel_to_best,
      primary_rec$geo_mean_speedup_vs_min_thread
    ),
    ""
  )
}

if (nrow(primary_sweeps)) {
  verdict_lines <- c(verdict_lines, "## By Sweep", "")
  for (i in seq_len(nrow(primary_sweeps))) {
    row <- primary_sweeps[i]
    verdict_lines <- c(
      verdict_lines,
      sprintf(
        "- `%s` sweep: recommend **%s**; fastest tested was **%s**; geometric mean relative runtime %.3fx.",
        row$sweep,
        fmt_threads(row$n_threads),
        fmt_threads(row$fastest_threads),
        row$geo_mean_rel_to_best
      )
    )
  }
  verdict_lines <- c(verdict_lines, "")
}

verdict_lines <- c(
  verdict_lines,
  "## Overall Metrics",
  "",
  "| metric | recommended_threads | fastest_threads | geo_mean_rel_to_best | median_rel_to_best | p90_rel_to_best | geo_mean_speedup_vs_min_thread |",
  "|:--|--:|--:|--:|--:|--:|--:|"
)
for (i in seq_len(nrow(overall_recommendation))) {
  row <- overall_recommendation[i]
  verdict_lines <- c(
    verdict_lines,
    sprintf(
      "| %s | %d | %d | %s | %s | %s | %s |",
      row$metric,
      row$n_threads,
      row$fastest_threads,
      fmt_num(row$geo_mean_rel_to_best),
      fmt_num(row$median_rel_to_best),
      fmt_num(row$p90_rel_to_best),
      fmt_num(row$geo_mean_speedup_vs_min_thread)
    )
  )
}

verdict_lines <- c(
  verdict_lines,
  "",
  "## Notes",
  "",
  "- The recommendation is intentionally the smallest thread count close to the fastest result, not necessarily the absolute fastest count.",
  "- `mean_fit_time_sec` is usually the right primary metric for algorithmic scaling.",
  "- `mean_total_time_sec` includes setup overhead and is useful when assessing the whole per-package benchmark call.",
  "- `mean_vcov_time_sec` is reported only when the underlying benchmark run recorded a separately timed covariance step.",
  "",
  "## Files",
  "",
  paste("- Run folders:", metadata_path),
  paste("- Combined summary:", summary_path),
  paste("- Normalized detail:", detail_path),
  paste("- Overall ranking:", overall_path),
  paste("- Spec recommendations:", spec_path)
)
if (length(plot_paths)) {
  verdict_lines <- c(verdict_lines, paste("- Plot:", plot_paths))
}

writeLines(verdict_lines, verdict_path)

message("Analysis complete.")
message("Verdict: ", verdict_path)
if (nrow(primary_rec)) {
  message("Recommended threads for ", config$metric, ": ", primary_rec$n_threads)
}
