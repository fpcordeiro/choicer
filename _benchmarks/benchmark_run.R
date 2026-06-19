# Shared benchmark execution helpers. These files are sourced by scripts under
# _benchmarks/ and are intentionally not part of the installed package.

`%||%` <- function(a, b) if (is.null(a)) b else a

bench_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc)) {
      first <- readLines(desc, n = 1L, warn = FALSE)
      if (length(first) && grepl("^Package:\\s+choicer\\s*$", first)) return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Could not locate the choicer repository root from ", start, ".")
}

bench_source <- function(path, root = bench_repo_root()) {
  source(file.path(root, path), chdir = TRUE)
}

bench_parse_cli <- function(args = commandArgs(trailingOnly = TRUE),
                            defaults = list()) {
  out <- defaults
  if (!length(args)) return(out)

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("CLI arguments must use --name=value form; got: ", arg)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    value <- if (length(pieces) == 1L) "TRUE" else paste(pieces[-1L], collapse = "=")

    old <- out[[key]]
    if (is.null(old)) {
      out[[key]] <- value
    } else if (is.logical(old)) {
      out[[key]] <- tolower(value) %in% c("true", "t", "1", "yes", "y")
    } else if (is.integer(old)) {
      out[[key]] <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    } else if (is.numeric(old)) {
      out[[key]] <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
    } else if (is.character(old)) {
      out[[key]] <- strsplit(value, ",", fixed = TRUE)[[1L]]
    } else {
      out[[key]] <- value
    }
  }
  out
}

bench_timestamp <- function(time = Sys.time()) {
  format(time, "%Y%m%d_%H%M%S")
}

bench_output_dir <- function(root, model, tag = NULL, timestamp = bench_timestamp()) {
  suffix <- if (!is.null(tag) && nzchar(tag)) paste0("_", gsub("[^A-Za-z0-9_.-]", "-", tag)) else ""
  out <- file.path(root, "_benchmarks", model, "output", paste0(timestamp, suffix))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

bench_load_choicer <- function(root = bench_repo_root(), clean_dll = TRUE) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  options(pkg.build_extra_flags = FALSE)

  if (requireNamespace("pkgload", quietly = TRUE)) {
    if (isTRUE(clean_dll) && requireNamespace("pkgbuild", quietly = TRUE)) {
      try(pkgbuild::clean_dll(root), silent = TRUE)
    }
    pkgload::load_all(root, quiet = TRUE, export_all = FALSE)
    return(invisible(TRUE))
  }

  if (requireNamespace("devtools", quietly = TRUE)) {
    if (isTRUE(clean_dll) && requireNamespace("pkgbuild", quietly = TRUE)) {
      try(pkgbuild::clean_dll(root), silent = TRUE)
    }
    devtools::load_all(root, quiet = TRUE, export_all = FALSE)
    return(invisible(TRUE))
  }

  suppressPackageStartupMessages(library(choicer))
  invisible(TRUE)
}

bench_spec_grid <- function(config) {
  grid_mode <- config$grid_mode %||% "two_sweeps"
  if (!grid_mode %in% c("two_sweeps", "full_grid")) {
    stop("config$grid_mode must be 'two_sweeps' or 'full_grid'.")
  }

  if (grid_mode == "two_sweeps") {
    n_sweep <- data.table::data.table(
      sweep = "N",
      dimension = "N",
      dimension_value = as.integer(config$n_vec),
      N = as.integer(config$n_vec),
      J = as.integer(config$fixed_j)
    )
    j_sweep <- data.table::data.table(
      sweep = "J",
      dimension = "J",
      dimension_value = as.integer(config$j_vec),
      N = as.integer(config$fixed_n),
      J = as.integer(config$j_vec)
    )
    grid <- data.table::rbindlist(list(n_sweep, j_sweep), use.names = TRUE)
  } else {
    grid <- data.table::CJ(
      N = as.integer(config$n_vec),
      J = as.integer(config$j_vec),
      sorted = FALSE
    )
    grid[, `:=`(
      sweep = "full_grid",
      dimension = "N_J",
      dimension_value = seq_len(.N)
    )]
    data.table::setcolorder(grid, c("sweep", "dimension", "dimension_value", "N", "J"))
  }

  grid[, spec_id := sprintf("%s_N%s_J%s", sweep, N, J)]
  grid[]
}

bench_time <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  list(value = value, elapsed = as.numeric(elapsed))
}

bench_quiet <- function(expr) {
  value <- NULL
  messages <- character()
  warnings <- character()
  output <- utils::capture.output({
    value <- withCallingHandlers(
      force(expr),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })
  list(value = value, output = output, messages = messages, warnings = warnings)
}

bench_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

bench_write_session_info <- function(path) {
  utils::capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}
