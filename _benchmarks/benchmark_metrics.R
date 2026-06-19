benchmark_compute_coef_checks <- function(coef_results, raw_results, config) {
  if (!nrow(coef_results)) return(coef_results)

  key_cols <- c("spec_id", "run", "parameter")
  ref <- coef_results[
    package == "choicer" & status == "ok",
    c(key_cols, "estimate"),
    with = FALSE
  ]
  data.table::setnames(ref, "estimate", "reference_estimate")

  out <- merge(coef_results, ref, by = key_cols, all.x = TRUE, sort = FALSE)
  out[, abs_diff := abs(estimate - reference_estimate)]
  out[, rel_diff := abs_diff / pmax(abs(reference_estimate), .Machine$double.eps)]
  out[, within_tolerance := data.table::fifelse(
    is.na(abs_diff),
    NA,
    abs_diff <= config$coef_abs_tol | rel_diff <= config$coef_rel_tol
  )]

  raw_ref <- raw_results[
    package == "choicer" & status == "ok",
    .(spec_id, run, reference_loglik = loglik)
  ]
  out <- merge(out, raw_ref, by = c("spec_id", "run"), all.x = TRUE, sort = FALSE)
  out[, loglik_abs_diff := abs(loglik - reference_loglik)]
  out[, loglik_rel_diff := loglik_abs_diff / pmax(abs(reference_loglik), .Machine$double.eps)]
  out[, loglik_within_tolerance := data.table::fifelse(
    is.na(loglik_rel_diff),
    NA,
    loglik_rel_diff <= config$loglik_rel_tol
  )]
  out[]
}

benchmark_summarise_results <- function(raw_results, coef_checks, config) {
  raw <- data.table::copy(raw_results)
  raw[, ok := status == "ok"]
  summary <- raw[, .(
    n_runs = .N,
    n_ok = sum(ok, na.rm = TRUE),
    n_converged = sum(converged %in% TRUE, na.rm = TRUE),
    success_rate = mean(ok, na.rm = TRUE),
    convergence_rate = mean(converged %in% TRUE, na.rm = TRUE),
    mean_fit_time_sec = mean(fit_time_sec, na.rm = TRUE),
    sd_fit_time_sec = stats::sd(fit_time_sec, na.rm = TRUE),
    median_fit_time_sec = stats::median(fit_time_sec, na.rm = TRUE),
    mean_total_time_sec = mean(total_time_sec, na.rm = TRUE),
    sd_total_time_sec = stats::sd(total_time_sec, na.rm = TRUE),
    median_total_time_sec = stats::median(total_time_sec, na.rm = TRUE),
    mean_loglik = mean(loglik, na.rm = TRUE)
  ), by = .(sweep, dimension, dimension_value, spec_id, N, J, package)]

  if (nrow(coef_checks)) {
    coef_summary <- coef_checks[status == "ok", .(
      estimates_within_tolerance = all(within_tolerance %in% TRUE),
      max_coef_abs_diff = max(abs_diff, na.rm = TRUE),
      max_coef_rel_diff = max(rel_diff, na.rm = TRUE),
      loglik_within_tolerance = all(loglik_within_tolerance %in% TRUE),
      max_loglik_abs_diff = max(loglik_abs_diff, na.rm = TRUE),
      max_loglik_rel_diff = max(loglik_rel_diff, na.rm = TRUE)
    ), by = .(sweep, dimension, dimension_value, spec_id, N, J, package)]
    summary <- merge(
      summary, coef_summary,
      by = c("sweep", "dimension", "dimension_value", "spec_id", "N", "J", "package"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    summary[, `:=`(
      estimates_within_tolerance = NA,
      max_coef_abs_diff = NA_real_,
      max_coef_rel_diff = NA_real_,
      loglik_within_tolerance = NA,
      max_loglik_abs_diff = NA_real_,
      max_loglik_rel_diff = NA_real_
    )]
  }

  numeric_cols <- names(summary)[vapply(summary, is.numeric, logical(1L))]
  for (col in numeric_cols) {
    data.table::set(summary, which(is.nan(summary[[col]])), col, NA_real_)
  }
  data.table::setorder(summary, sweep, dimension_value, package)
  summary[]
}
