benchmark_git_value <- function(root, args) {
  out <- tryCatch(
    system2("git", c("-C", root, args), stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  if (!length(out)) "" else paste(out, collapse = "\n")
}

benchmark_capture_metadata <- function(root, model, config, packages,
                                       output_dir, output_files = character()) {
  pkg_versions <- stats::setNames(
    vapply(packages, bench_package_version, character(1L)),
    packages
  )
  missing_packages <- names(pkg_versions)[is.na(pkg_versions)]
  ext <- extSoftVersion()
  ext_value <- function(name) {
    if (name %in% names(ext)) unname(ext[[name]]) else NA_character_
  }

  list(
    benchmark = model,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    run_tag = config$tag %||% NA_character_,
    command = paste(commandArgs(FALSE), collapse = " "),
    root = root,
    output_dir = output_dir,
    output_files = as.list(output_files),
    git = list(
      sha = benchmark_git_value(root, c("rev-parse", "HEAD")),
      branch = benchmark_git_value(root, c("rev-parse", "--abbrev-ref", "HEAD")),
      status_short = benchmark_git_value(root, c("status", "--short"))
    ),
    r = list(
      version = R.version.string,
      platform = R.version$platform,
      blas = ext_value("BLAS"),
      lapack = ext_value("LAPACK")
    ),
    config = config,
    packages = as.list(pkg_versions),
    missing_packages = as.list(missing_packages)
  )
}

benchmark_write_metadata <- function(metadata, output_dir) {
  json_path <- file.path(output_dir, "metadata.json")
  txt_path <- file.path(output_dir, "metadata.txt")

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(metadata, json_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    dput(metadata, file = json_path)
  }

  con <- file(txt_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    paste("Benchmark:", metadata$benchmark),
    paste("Timestamp:", metadata$timestamp),
    paste("Tag:", metadata$run_tag),
    paste("Repository:", metadata$root),
    paste("Output:", metadata$output_dir),
    paste("Git SHA:", metadata$git$sha),
    paste("Git branch:", metadata$git$branch),
    "",
    "Command:",
    metadata$command,
    "",
    "Package versions:"
  ), con)
  for (nm in names(metadata$packages)) {
    writeLines(sprintf("  %s: %s", nm, metadata$packages[[nm]] %||% NA_character_), con)
  }
  if (length(metadata$missing_packages)) {
    writeLines(c("", "Missing packages:", paste(" ", unlist(metadata$missing_packages))), con)
  }
  if (nzchar(metadata$git$status_short %||% "")) {
    writeLines(c("", "Git status:", metadata$git$status_short), con)
  }

  invisible(c(json_path, txt_path))
}
