# Shared helpers for the population-scale benchmarks (development only; not part
# of the package). Sourced by the scripts in this directory, which run one
# configuration per R process:
#   R_LIBS=<library holding one choicer build> Rscript _benchmarks/scale/<script>.R ...
#
# - load_choicer() loads choicer from the first R_LIBS entry and stops if it came
#   from elsewhere (a site library's release would otherwise shadow the build).
# - gen_regime() builds synthetic claims, census or 4-alternative data directly in
#   the prepared (stacked) layout with gen_data.cpp, deterministic per unit.
# - The mc_*() and rss_*() helpers read the macOS malloc-counting shim
#   (memcount.c, preloaded by the Rmc front end that make_rmc.sh writes); without
#   the shim they return NA.

.script_file <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) != 1L) stop("Run the benchmark scripts with Rscript or R --file=.")
  normalizePath(f, mustWork = TRUE)
}
# The scripts' directory: from the running script's --file=, unless the caller
# set SCALE_DIR before sourcing this file (as `Rscript -e` code must)
if (!exists("SCALE_DIR", inherits = FALSE)) SCALE_DIR <- dirname(.script_file())
`%||%` <- function(a, b) if (is.null(a)) b else a

# Command-line options key=value or --key=value; a bare --flag is "1".
scale_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- sub("^--", "", args)
  kv <- strsplit(args, "=", fixed = TRUE)
  setNames(lapply(kv, function(x) if (length(x) > 1L) paste(x[-1L], collapse = "=") else "1"),
           vapply(kv, `[[`, "", 1L))
}

# The output root: output/ next to the scripts (ignored by git), unless given.
scale_output <- function(...) {
  d <- file.path(Sys.getenv("CHOICER_SCALE_OUT", file.path(SCALE_DIR, "output")), ...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

load_choicer <- function() {
  lib <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1]]
  lib <- lib[nzchar(lib)]
  if (!length(lib)) stop("Set R_LIBS to the library that holds the choicer build to measure.")
  lib <- normalizePath(lib[1L], mustWork = TRUE)
  suppressMessages(library(choicer, lib.loc = lib))
  p <- getNamespaceInfo("choicer", "path")
  if (normalizePath(dirname(p)) != lib) stop("choicer loaded from ", p, ", not from ", lib)
  invisible(lib)
}

# Physical cores this process may use, and how they were counted. On Linux
# parallel::detectCores() counts logical CPUs whatever its `logical` argument,
# so count the distinct (core, socket) pairs among the CPUs of the process's
# affinity mask, from lscpu or else from sysfs; elsewhere detectCores(logical =
# FALSE). The count must lie between 1 and the logical CPUs.
physical_cores <- function() {
  logical <- parallel::detectCores()
  ok <- function(n) length(n) == 1L && is.finite(n) && n >= 1 && (is.na(logical) || n <= logical)
  if (file.exists("/proc/self/status")) {
    allowed <- NULL
    a <- grep("^Cpus_allowed_list:", readLines("/proc/self/status"), value = TRUE)
    if (length(a)) {
      ranges <- strsplit(trimws(sub("^Cpus_allowed_list:", "", a)), ",", fixed = TRUE)[[1]]
      allowed <- unlist(lapply(ranges, function(r) {
        b <- as.integer(strsplit(r, "-", fixed = TRUE)[[1]])
        seq(b[1L], b[length(b)])
      }))
    }
    cpus <- tryCatch(suppressWarnings(system("lscpu -p=CPU,CORE,SOCKET", intern = TRUE,
                                             ignore.stderr = TRUE)), error = function(e) character())
    cpus <- cpus[nzchar(cpus) & !startsWith(cpus, "#")]
    if (length(cpus)) {
      f <- tryCatch(utils::read.csv(text = cpus, header = FALSE, colClasses = "character",
                                    na.strings = c("", "-"), fill = TRUE),
                    error = function(e) NULL)
      if (!is.null(f) && ncol(f) >= 2L) {
        if (ncol(f) < 3L) f[[3L]] <- NA_character_
        f[[3L]][is.na(f[[3L]])] <- "0"  # no socket field: one socket
        if (!is.null(allowed)) f <- f[as.integer(f[[1L]]) %in% allowed, , drop = FALSE]
        if (nrow(f) && !anyNA(f[[2L]])) {
          n <- nrow(unique(f[, 2:3]))
          if (ok(n)) return(list(n = n, method = "lscpu (core, socket) pairs"))
        }
      }
    }
    # sysfs: one sibling set (the logical CPUs of a core) per physical core;
    # core_id alone restarts within each cluster on some ARM machines
    dirs <- Sys.glob("/sys/devices/system/cpu/cpu[0-9]*/topology")
    if (length(dirs)) {
      id <- as.integer(sub(".*/cpu([0-9]+)/topology$", "\\1", dirs))
      if (!is.null(allowed)) dirs <- dirs[id %in% allowed]
      rd <- function(d, f) {
        x <- tryCatch(readLines(file.path(d, f), warn = FALSE), error = function(e) character())
        if (length(x)) x[1L] else NA_character_
      }
      sib <- vapply(dirs, function(d) {
        x <- rd(d, "core_cpus_list")
        if (is.na(x)) rd(d, "thread_siblings_list") else x
      }, "")
      if (length(sib) && !anyNA(sib)) {
        n <- length(unique(sib))
        if (ok(n)) return(list(n = n, method = "sysfs (distinct sibling sets)"))
      }
    }
  }
  n <- parallel::detectCores(logical = FALSE)
  linux <- file.exists("/proc/self/status")  # where logical = FALSE is ignored
  if (ok(n)) {
    return(list(n = n, method = if (linux) "detectCores() (logical CPUs; no topology found)" else
                                 "detectCores(logical = FALSE)"))
  }
  list(n = if (ok(logical)) logical else 1L, method = "detectCores() (logical CPUs)")
}

loadavg <- function() {
  if (file.exists("/proc/loadavg")) return(paste(readLines("/proc/loadavg", 1L), collapse = " "))
  paste(tryCatch(system("sysctl -n vm.loadavg", intern = TRUE), error = function(e) NA), collapse = " ")
}

# The malloc shim (macOS): heap blocks, bytes moved and peak live heap, plus a
# sampled peak RSS.
shim_loaded <- function() nzchar(Sys.getenv("DYLD_INSERT_LIBRARIES")) && nzchar(Sys.getenv("MC_SHIM"))
if (shim_loaded()) dyn.load(Sys.getenv("MC_SHIM"))
mc_reset <- function() if (shim_loaded()) invisible(.C("mc_reset"))
mc_disable <- function() if (shim_loaded()) invisible(.C("mc_disable"))
mc_get <- function() {
  nm <- c("n_alloc", "bytes_alloc", "live", "peak_live", "n_free", "bytes_freed")
  if (!shim_loaded()) return(setNames(rep(NA_real_, 6), nm))
  setNames(.C("mc_get", out = double(6))$out, nm)
}
rss_start <- function(ms = 2) if (shim_loaded()) invisible(.C("mc_rss_start", as.double(ms)))
rss_stop <- function() {
  nm <- c("peak_rss", "peak_fp", "rss", "fp")
  if (!shim_loaded()) return(setNames(rep(NA_real_, 4), nm))
  setNames(.C("mc_rss_stop", out = double(4))$out, nm)
}
rss_now <- function() {
  if (!shim_loaded()) return(c(rss = NA_real_, fp = NA_real_))
  setNames(.C("mc_rss_now", out = double(2))$out, c("rss", "fp"))
}

gen_env <- new.env()
load_gen <- function() {
  if (!exists("gen_mxl_data", envir = gen_env)) {
    Rcpp::sourceCpp(file.path(SCALE_DIR, "gen_data.cpp"), env = gen_env,
                    cacheDir = scale_output("cpp-cache"))
  }
  get("gen_mxl_data", envir = gen_env)
}

K_X <- 10L; K_W <- 3L
true_par <- function(K_x = K_X, K_w = K_W) {
  beta <- 0.5 * seq(-1, 1, length.out = K_x)
  mu <- c(-1, 0.5, 1)[seq_len(K_w)]
  L <- diag(c(0.8, 0.6, 0.5)[seq_len(K_w)], K_w)
  if (K_w >= 2) L[2, 1] <- 0.2
  if (K_w >= 3) { L[3, 1] <- -0.1; L[3, 2] <- 0.15 }
  list(beta = beta, mu = mu, L = L)
}

# Regimes:
#   claims    patients with a Pareto(1.1) number of visits (up to 500), each
#             choosing among 20-200 of 500 hospitals in their region;
#   census    students with mostly one choice (up to three), among 10-100 of
#             3,000 schools in their region;
#   census20k the census design with 20,000 schools (about 17,000 parameters);
#   small     a cross-section of 4 alternatives, one choice per unit.
gen_regime <- function(regime, U, seed = 1, fat = FALSE, threads = 11L,
                       K_x = K_X, K_w = K_W) {
  f <- load_gen()
  tp <- true_par(K_x, K_w)
  cfg <- switch(regime,
    claims = list(J_total = 500L, n_regions = 50L, m_min = 20L, m_max = 200L,
                  keep_p = 0.9, t_dist = "pareto", alpha = 1.1, t_min = 2L,
                  T_max = 500L, p_single = 0),
    census = list(J_total = 3000L, n_regions = 150L, m_min = 10L, m_max = 100L,
                  keep_p = 1.0, t_dist = "short", alpha = 1, t_min = 1L,
                  T_max = 3L, p_single = 0.85),
    census20k = list(J_total = 20000L, n_regions = 1000L, m_min = 10L, m_max = 100L,
                     keep_p = 1.0, t_dist = "short", alpha = 1, t_min = 1L,
                     T_max = 3L, p_single = 0.85),
    small = list(J_total = 4L, n_regions = 1L, m_min = 4L, m_max = 4L,
                 keep_p = 1.0, t_dist = "short", alpha = 1, t_min = 1L,
                 T_max = 1L, p_single = 1.0),
    stop("unknown regime ", regime))
  d <- f(U, K_x, K_w, cfg$J_total, cfg$n_regions, cfg$m_min, cfg$m_max,
         cfg$keep_p, cfg$t_dist, cfg$alpha, cfg$t_min, cfg$T_max, cfg$p_single,
         tp$beta, tp$mu, tp$L, 0.5, if (fat) 1 else 0, if (fat) 3 else 1,
         seed, threads)
  d$weights <- rep(1, length(d$M))
  d$regime <- regime
  d$U <- U
  # True theta in choicer's order: beta, mu, L (row-major lower triangle, log
  # diagonal), delta_2..J
  Lp <- unlist(lapply(seq_len(K_w), function(i) sapply(seq_len(i), function(j)
    if (i == j) log(tp$L[i, i]) else tp$L[i, j])))
  d$theta_true <- c(tp$beta, tp$mu, Lp, d$delta[-1])
  d$K_x <- K_x; d$K_w <- K_w
  d
}

# The 0.2.x cold start (zeros, log 0.5 on the Cholesky diagonal), used to time
# the kernels at a fixed point.
theta_start <- function(d) {
  th <- rep(0, length(d$theta_true))
  th[d$K_x + d$K_w + cumsum(seq_len(d$K_w))] <- log(0.5)
  th
}

# Rows of each likelihood unit (decision maker; a situation without Ti).
unit_rows <- function(d) {
  if (is.null(d$Ti)) return(as.numeric(d$M))
  diff(c(0, cumsum(as.numeric(d$M))[cumsum(as.numeric(d$Ti))]))
}

describe <- function(d) {
  T_q <- stats::quantile(d$Ti, c(0.5, 0.9, 0.99, 1))
  cat(sprintf("[%s] U=%d N=%.0f rows=%.3g J=%d n_params=%d | T: median %g p90 %g p99 %g max %g | M: mean %.1f max %d | max unit rows %.0f | X+W %.2f GB\n",
              d$regime, d$U, d$N, d$rows, d$J, length(d$theta_true), T_q[1], T_q[2], T_q[3], T_q[4],
              mean(d$M), max(d$M), max(unit_rows(d)), 8 * d$rows * (d$K_x + d$K_w) / 1e9))
}
