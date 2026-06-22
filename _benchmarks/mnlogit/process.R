mnlogit_safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

mnlogit_output_path <- function(output_dir, path) {
  if (grepl("^([A-Za-z]:)?[\\/]", path)) return(path)
  file.path(output_dir, path)
}

mnlogit_attempt_grid <- function(specs, config) {
  n_total <- nrow(specs) * config$n_runs * length(config$packages)
  rows <- vector("list", n_total)
  idx <- 1L
  timeout_sec <- if (is.finite(config$max_run_sec)) config$max_run_sec else NA_real_

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i]
    for (run in seq_len(config$n_runs)) {
      seed <- config$seed + i * 10000L + run
      for (pkg in config$packages) {
        attempt_id <- mnlogit_safe_filename(sprintf(
          "%04d_%s_run%03d_%s",
          idx, spec$spec_id, run, pkg
        ))
        dep <- mnlogit_dependency(pkg)
        rows[[idx]] <- data.table::data.table(
          benchmark = "mnlogit",
          attempt_index = idx,
          attempt_id = attempt_id,
          sweep = spec$sweep,
          dimension = spec$dimension,
          dimension_value = spec$dimension_value,
          spec_id = spec$spec_id,
          N = as.integer(spec$N),
          J = as.integer(spec$J),
          run = as.integer(run),
          seed = as.integer(seed),
          package = pkg,
          package_dependency = dep,
          package_version = bench_package_version(dep),
          status = "pending",
          started_at = NA_character_,
          finished_at = NA_character_,
          elapsed_wall_sec = NA_real_,
          exit_status = NA_integer_,
          timeout_sec = timeout_sec,
          raw_path = file.path("partials", paste0(attempt_id, "_raw.csv")),
          coef_path = file.path("partials", paste0(attempt_id, "_coef.csv")),
          stdout_path = file.path("logs", paste0(attempt_id, ".out")),
          stderr_path = file.path("logs", paste0(attempt_id, ".err")),
          error = NA_character_
        )
        idx <- idx + 1L
      }
    }
  }

  data.table::rbindlist(rows, use.names = TRUE)
}

mnlogit_prepare_attempt_dirs <- function(output_dir) {
  dirs <- file.path(output_dir, c("partials", "logs", "requests"))
  for (dir in dirs) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

mnlogit_write_run_status <- function(status, status_path) {
  bench_atomic_fwrite(status, status_path)
}

mnlogit_set_attempt_status <- function(status, status_path, attempt_id, values) {
  row_idx <- match(attempt_id, status$attempt_id)
  if (is.na(row_idx)) stop("Unknown benchmark attempt: ", attempt_id)
  for (nm in names(values)) {
    data.table::set(status, i = row_idx, j = nm, value = values[[nm]])
  }
  mnlogit_write_run_status(status, status_path)
  invisible(status)
}

mnlogit_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
}

mnlogit_log_tail <- function(path, n = 40L) {
  if (!file.exists(path)) return(NA_character_)
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (!length(txt)) return(NA_character_)
  paste(utils::tail(txt, n), collapse = "\n")
}

mnlogit_clean_process_warnings <- function(warnings) {
  if (!length(warnings)) return(character())
  vapply(warnings, function(w) {
    if (grepl("timed out after", w, ignore.case = TRUE)) {
      return(sub(".*timed out after ([^[:space:]]+).*", "Worker process timed out after \\1.", w))
    }
    w
  }, character(1L), USE.NAMES = FALSE)
}

mnlogit_process_error_text <- function(stderr_path, stdout_path, process_warnings = character()) {
  pieces <- c(
    mnlogit_clean_process_warnings(process_warnings),
    mnlogit_log_tail(stderr_path),
    mnlogit_log_tail(stdout_path, n = 20L)
  )
  pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
  if (!length(pieces)) return(NA_character_)
  paste(pieces, collapse = "\n")
}

mnlogit_write_failure_partial <- function(attempt, output_dir, status, elapsed_wall_sec,
                                          exit_status, error, process_warnings = character()) {
  timeout_sec <- if (is.finite(attempt$timeout_sec)) attempt$timeout_sec else NA_real_
  raw <- mnlogit_raw_row(
    package = attempt$package,
    spec = attempt,
    run = attempt$run,
    status = status,
    error = error,
    warnings = if (length(process_warnings)) paste(process_warnings, collapse = " | ") else NA_character_,
    timeout_sec = timeout_sec,
    process_elapsed_sec = elapsed_wall_sec,
    worker_exit_status = exit_status
  )
  bench_atomic_fwrite(raw, mnlogit_output_path(output_dir, attempt$raw_path))
  invisible(raw)
}

mnlogit_update_partial_process_metadata <- function(raw_path, timeout_sec,
                                                    elapsed_wall_sec, exit_status) {
  raw <- data.table::fread(raw_path)
  if (!"timeout_sec" %in% names(raw)) raw[, timeout_sec := NA_real_]
  if (!"process_elapsed_sec" %in% names(raw)) raw[, process_elapsed_sec := NA_real_]
  if (!"worker_exit_status" %in% names(raw)) raw[, worker_exit_status := NA_integer_]
  timeout_value <- timeout_sec
  elapsed_value <- elapsed_wall_sec
  exit_value <- as.integer(exit_status)
  raw[, `:=`(
    timeout_sec = timeout_value,
    process_elapsed_sec = elapsed_value,
    worker_exit_status = exit_value
  )]
  bench_atomic_fwrite(raw, raw_path)
  raw[]
}

mnlogit_reconcile_partial_outputs <- function(output_dir, status_path, config) {
  if (!file.exists(status_path)) return(NULL)
  status <- data.table::fread(status_path)

  read_outputs <- function(paths) {
    paths <- paths[!is.na(paths) & nzchar(paths)]
    paths <- vapply(paths, mnlogit_output_path, character(1L), output_dir = output_dir)
    paths <- paths[file.exists(paths)]
    if (!length(paths)) return(data.table::data.table())
    data.table::rbindlist(lapply(paths, data.table::fread), use.names = TRUE, fill = TRUE)
  }

  raw_results <- read_outputs(status$raw_path)
  if (!nrow(raw_results)) return(NULL)

  coef_results <- read_outputs(status$coef_path)
  coef_results <- benchmark_compute_coef_checks(coef_results, raw_results, config)
  summary_results <- benchmark_summarise_results(raw_results, coef_results, config)

  raw_path <- file.path(output_dir, "raw_results.csv")
  coef_path <- file.path(output_dir, "coef_results.csv")
  summary_path <- file.path(output_dir, "summary_results.csv")

  bench_atomic_fwrite(raw_results, raw_path)
  bench_atomic_fwrite(coef_results, coef_path)
  bench_atomic_fwrite(summary_results, summary_path)

  list(
    raw_results = raw_results,
    coef_results = coef_results,
    summary_results = summary_results,
    raw_path = raw_path,
    coef_path = coef_path,
    summary_path = summary_path
  )
}

mnlogit_run_child_attempt <- function(attempt, root, output_dir, config,
                                      status, status_path) {
  raw_path <- mnlogit_output_path(output_dir, attempt$raw_path)
  coef_path <- mnlogit_output_path(output_dir, attempt$coef_path)
  stdout_path <- mnlogit_output_path(output_dir, attempt$stdout_path)
  stderr_path <- mnlogit_output_path(output_dir, attempt$stderr_path)
  request_path <- file.path(output_dir, "requests", paste0(attempt$attempt_id, ".rds"))
  timeout_arg <- if (is.finite(config$max_run_sec)) config$max_run_sec else 0
  started_at <- mnlogit_timestamp()

  mnlogit_set_attempt_status(status, status_path, attempt$attempt_id, list(
    status = "running",
    started_at = started_at,
    finished_at = NA_character_,
    elapsed_wall_sec = NA_real_,
    exit_status = NA_integer_,
    error = NA_character_
  ))

  request <- list(
    root = root,
    attempt = as.list(attempt),
    config = config,
    paths = list(raw = raw_path, coef = coef_path)
  )
  bench_atomic_save_rds(request, request_path)
  on.exit(if (file.exists(request_path)) unlink(request_path), add = TRUE)

  process_warnings <- character()
  process_error <- NA_character_
  elapsed_start <- proc.time()[["elapsed"]]
  worker_path <- file.path(root, "_benchmarks", "mnlogit", "worker.R")
  rscript <- file.path(R.home("bin"), "Rscript")
  exit_status <- tryCatch(
    withCallingHandlers(
      system2(
        rscript,
        c("--vanilla", worker_path, paste0("--request=", request_path)),
        stdout = stdout_path,
        stderr = stderr_path,
        timeout = timeout_arg
      ),
      warning = function(w) {
        process_warnings <<- c(process_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      process_error <<- conditionMessage(e)
      1L
    }
  )
  elapsed_wall_sec <- as.numeric(proc.time()[["elapsed"]] - elapsed_start)
  exit_status <- as.integer(exit_status[[1L]] %||% 0L)
  process_warnings <- mnlogit_clean_process_warnings(process_warnings)
  timed_out <- is.finite(config$max_run_sec) && (
    identical(exit_status, 124L) ||
      any(grepl("timed out|timeout", process_warnings, ignore.case = TRUE)) ||
      (!file.exists(raw_path) && elapsed_wall_sec >= config$max_run_sec)
  )

  timeout_sec <- if (is.finite(config$max_run_sec)) config$max_run_sec else NA_real_
  if (identical(exit_status, 0L) && file.exists(raw_path)) {
    raw <- mnlogit_update_partial_process_metadata(
      raw_path = raw_path,
      timeout_sec = timeout_sec,
      elapsed_wall_sec = elapsed_wall_sec,
      exit_status = exit_status
    )
    attempt_status <- unique(raw$status)
    attempt_status <- if (length(attempt_status) == 1L) attempt_status[[1L]] else "ok"
    attempt_error <- raw$error[!is.na(raw$error) & nzchar(raw$error)]
    attempt_error <- if (length(attempt_error)) paste(unique(attempt_error), collapse = " | ") else NA_character_
  } else {
    attempt_status <- if (timed_out) "timeout" else "error"
    attempt_error <- if (!is.na(process_error)) {
      process_error
    } else if (timed_out) {
      sprintf("Worker process exceeded max_run_sec = %s.", config$max_run_sec)
    } else if (!file.exists(raw_path)) {
      "Worker process did not produce a raw partial result."
    } else {
      sprintf("Worker process exited with status %s.", exit_status)
    }
    log_error <- if (timed_out) {
      NA_character_
    } else {
      mnlogit_process_error_text(stderr_path, stdout_path, process_warnings)
    }
    if (!is.na(log_error)) attempt_error <- paste(c(attempt_error, log_error), collapse = "\n")
    mnlogit_write_failure_partial(
      attempt = attempt,
      output_dir = output_dir,
      status = attempt_status,
      elapsed_wall_sec = elapsed_wall_sec,
      exit_status = exit_status,
      error = attempt_error,
      process_warnings = process_warnings
    )
  }

  mnlogit_set_attempt_status(status, status_path, attempt$attempt_id, list(
    status = attempt_status,
    finished_at = mnlogit_timestamp(),
    elapsed_wall_sec = elapsed_wall_sec,
    exit_status = exit_status,
    error = attempt_error
  ))

  invisible(attempt_status)
}
