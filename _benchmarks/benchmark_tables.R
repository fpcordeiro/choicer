benchmark_write_summary_table <- function(summary_results, output_dir) {
  table_dt <- data.table::copy(summary_results)
  keep <- c(
    "sweep", "N", "J", "package", "n_ok", "n_runs",
    "mean_fit_time_sec", "sd_fit_time_sec", "median_fit_time_sec",
    "convergence_rate", "estimates_within_tolerance",
    "loglik_within_tolerance", "max_coef_abs_diff", "max_loglik_rel_diff"
  )
  table_dt <- table_dt[, keep, with = FALSE]
  table_dt[, runs := paste0(n_ok, "/", n_runs)]
  table_dt[, c("n_ok", "n_runs") := NULL]

  round_cols <- c(
    "mean_fit_time_sec", "sd_fit_time_sec", "median_fit_time_sec",
    "convergence_rate", "max_coef_abs_diff", "max_loglik_rel_diff"
  )
  for (col in intersect(round_cols, names(table_dt))) {
    table_dt[, (col) := round(get(col), 6)]
  }
  data.table::setcolorder(table_dt, c("sweep", "N", "J", "package", "runs"))

  path <- file.path(output_dir, "summary_table.md")
  if (requireNamespace("knitr", quietly = TRUE)) {
    txt <- knitr::kable(table_dt, format = "markdown")
    writeLines(txt, path)
  } else {
    con <- file(path, open = "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(table_dt, con, sep = " | ", row.names = FALSE, quote = FALSE)
  }
  path
}
