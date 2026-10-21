#!/usr/bin/env Rscript
# A sweep of the mixed logit kernels' draw-batch budget on this machine.
#
# Each estimation kernel forms a decision maker's rows-by-draws matrices a batch
# of draws at a time, so that a thread's working memory stays within a fixed
# budget, MXL_BATCH_BYTES (src/mxlogit.cpp; 4 MiB): a unit of R rows is split
# into batches when n_mats * 8 * R * S bytes exceed it (n_mats = 2 for the
# gradient and BHHH kernels, 1 for the Hessian and the conditional tastes). A
# batch that stays in a core's cache is fastest; 4 MiB lies in the range that
# ran fastest (1-16 MiB) on the 11-core laptop where it was tuned. Where many
# cores share a cache, another budget may be faster. The budget is fixed when
# the package is compiled, so this script builds one library per budget from a
# source tree (-DCHOICER_MXL_BATCH_BYTES=<bytes>), times the same kernels on the
# same synthetic data with each, and reports times relative to 4 MiB.
#
#   Rscript _benchmarks/scale/budget_sweep.R --src=<choicer source tree> \
#     [--budgets=262144,1048576,2097152,4194304,8388608,16777216,67108864] \
#     [--threads=<half>,<all physical cores>] \
#     [--cases=heavy,claims,claims_s500,claims_s1000,census_s1000,census,small] \
#     [--kernels=grad,tastes,bhhh] [--rounds=3] [--reps=3] [--out=<dir>] [--quick] [--no-build]
#
# --quick runs every case at a tenth of its size or less, with two budgets
# unless --budgets is given: a check of the setup, not a measurement.
# --no-build reuses the libraries already in --out. A run resumes, with
# --no-build, in an --out that holds results of the same configuration (more
# --rounds add rounds); it refuses another configuration, and refuses to
# rebuild into an --out that holds results. The package's tests assume the
# default budget (some check its batch thresholds and the memory check's
# messages): do not run them against these libraries.
#
# Libraries go to <out>/lib-<budget>, each built from a fresh copy of the tree
# with the machine's own Makevars (tools::makevars_user(), copied in) plus the
# define, and checked for its namespace path, OpenMP and its budget (the
# per-thread memory the kernels' memory check reports includes it). One R
# process per (round, case, budget), the budgets in an order rotated by round
# and case and reversed in alternate rounds, with the BLAS single-threaded (the
# kernels call it inside their parallel regions). A budget that batches every
# unit of a case exactly as another does is run once for that case and
# reported with that budget's times. Results: <out>/machine.txt,
# <out>/raw.csv (every timed call), <out>/summary.csv and a table on the
# console. Decision rule (README.md): change the budget only if one beats
# 4 MiB by more than 5% at the largest thread count in every round on each
# long panel (claims_s500, claims_s1000, census_s1000) where it batches at
# least a fifth of the rows differently, and loses no more than 2% wherever
# else it batches a twentieth of the rows or more differently.
#
# On a multi-socket server consider `numactl --interleave=all Rscript ...`: the
# synthetic data are written by one thread, so they sit on one node otherwise.
# A ~/.Renviron that sets R_LIBS overrides the R_LIBS this script passes to its
# processes; the library checks would then fail, by design.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
opt <- scale_args()
REF <- 4194304

CASES <- list(
  heavy        = list(regime = "claims", U = 1500L,   S = 1000L, quick_U = 300L,   quick_S = 300L),
  claims       = list(regime = "claims", U = 15000L,  S = 100L,  quick_U = 1500L,  quick_S = 100L),
  claims_s500  = list(regime = "claims", U = 15000L,  S = 500L,  quick_U = 1500L,  quick_S = 500L),
  claims_s1000 = list(regime = "claims", U = 15000L,  S = 1000L, quick_U = 1500L,  quick_S = 1000L),
  census_s1000 = list(regime = "census", U = 200000L, S = 1000L, quick_U = 20000L, quick_S = 1000L),
  census       = list(regime = "census", U = 200000L, S = 100L,  quick_U = 20000L, quick_S = 100L),
  small        = list(regime = "small",  U = 250000L, S = 100L,  quick_U = 25000L, quick_S = 100L))
# The long panels whose results decide. heavy (1,500 patients) is reported but
# not decisive: at many threads its few largest patients set the wall time.
LONG <- c("claims_s500", "claims_s1000", "census_s1000")
N_MATS <- c(grad = 2L, bhhh = 2L, tastes = 1L, hess = 1L)

case_spec <- function(name, quick) {
  cs <- CASES[[name]]
  if (is.null(cs)) stop("unknown case ", name, "; cases: ", paste(names(CASES), collapse = ", "))
  if (quick) { cs$U <- cs$quick_U; cs$S <- cs$quick_S }
  cs
}

# The batches of each unit of R rows (mxl_batch_count() in src/mxlogit.cpp).
batch_count <- function(R, S, n_mats, budget) {
  b_max <- budget / (8 * n_mats * pmax(R, 1))
  ifelse(b_max >= S, 1L, as.integer(ceiling(S / pmax(1, floor(b_max)))))
}

# ---------------------------------------------------------------------------
# Worker: one (round, case, budget), in its own process.
if (!is.null(opt$worker)) {
  load_choicer()
  ns <- asNamespace("choicer")
  quick <- opt$quick == "1"
  cs <- case_spec(opt$case, quick)
  budget <- as.numeric(opt$budget)
  threads <- as.integer(strsplit(opt$threads, ",")[[1]])
  kernels <- strsplit(opt$kernels, ",")[[1]]
  reps <- as.integer(opt$reps)
  d <- gen_regime(cs$regime, cs$U, threads = max(threads))
  args <- list(theta = theta_start(d), X = d$X, W = d$W, alt_idx = d$alt_idx,
               choice_idx = d$choice_idx, M = d$M, eta_draws = array(0, dim = c(d$K_w, 0L, 0L)),
               rc_dist = rep(0L, d$K_w), rc_correlation = TRUE, rc_mean = TRUE, use_asc = TRUE,
               include_outside_option = FALSE, gen_seed = 17L, gen_scramble = 1L, gen_S = cs$S,
               Ti = d$Ti)
  wargs <- c(args, list(weights = d$weights))
  call_kernel <- function(k) switch(k,
    grad = do.call(ns$mxl_loglik_gradient_parallel, wargs),
    bhhh = do.call(ns$mxl_bhhh_parallel, wargs),
    hess = do.call(ns$mxl_hessian_parallel, wargs),
    tastes = do.call(ns$mxl_conditional_tastes_parallel, args),
    stop("unknown kernel ", k))
  check_value <- function(x) {
    if (is.list(x) && !is.null(x$objective)) return(x$objective)
    sum(vapply(if (is.list(x)) x else list(x), function(e) sum(as.numeric(e)), 0))
  }
  # Resident memory (Linux): the high-water mark is reset before each kernel,
  # and the kernel's peak is reported as its increase over the resident size
  # right after the reset (heap kept from earlier calls is resident already)
  vm_status <- function(key) {
    if (!file.exists("/proc/self/status")) return(NA_real_)
    v <- grep(paste0("^", key, ":"), readLines("/proc/self/status"), value = TRUE)
    if (!length(v)) NA_real_ else as.numeric(gsub("[^0-9]", "", v)) * 1024
  }
  reset_hwm <- function() isTRUE(tryCatch({ writeLines("5", "/proc/self/clear_refs"); TRUE },
                                          error = function(e) FALSE, warning = function(w) FALSE))
  R_u <- unit_rows(d)
  rows <- list()
  for (k in kernels) for (nt in threads) {
    set_num_threads(nt)
    ti <- thread_info()
    invisible(gc())
    reset_ok <- reset_hwm()
    rss0 <- vm_status("VmRSS")
    reset_ok <- reset_ok && isTRUE(vm_status("VmHWM") <= rss0 + 4 * 2^20)
    cv <- check_value(call_kernel(k))
    tt <- vapply(seq_len(reps), function(r) system.time(call_kernel(k), gcFirst = FALSE)[["elapsed"]], 0)
    vm_add <- vm_status("VmHWM") - rss0
    bc <- batch_count(R_u, cs$S, N_MATS[[k]], budget)
    rows[[length(rows) + 1L]] <- data.frame(
      round = as.integer(opt$round), case = opt$case, regime = cs$regime, U = cs$U, S = cs$S,
      quick = quick, kernel = k, threads = nt, team = ti$omp_get_num_threads,
      budget = budget, rep = seq_len(reps), sec = tt, check = cv, split_units = mean(bc > 1L),
      split_rows = sum(R_u[bc > 1L]) / sum(R_u), max_batches = max(bc),
      vm_add = if (reset_ok) vm_add else NA_real_, load = loadavg(), stringsAsFactors = FALSE)
    cat(sprintf("  round %s %-12s %-6s %3d threads, budget %9.0f: median %.4f s (units split %.3f)\n",
                opt$round, opt$case, k, nt, budget, stats::median(tt), mean(bc > 1L)))
  }
  saveRDS(do.call(rbind, rows), opt$out)
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# Driver
quick <- !is.null(opt$quick)
no_build <- !is.null(opt[["no-build"]])
src <- opt$src
if (is.null(src) && !no_build) stop("--src=<choicer source tree> is required, or --no-build with --out")
if (!is.null(src)) src <- normalizePath(src, mustWork = TRUE)
budget_txt <- strsplit(opt$budgets %||% if (quick) "1048576,4194304" else
  "262144,1048576,2097152,4194304,8388608,16777216,67108864", ",")[[1]]
budgets <- suppressWarnings(as.numeric(budget_txt))
if (anyNA(budgets)) stop("budgets are numbers of bytes; could not read: ",
                         paste(budget_txt[is.na(budgets)], collapse = ", "))
budgets <- sort(unique(c(budgets, REF)))
if (any(budgets < 1024 | budgets > 2^34)) stop("budgets are in bytes, from 1024 to 2^34 (16 GiB)")
pc <- physical_cores()
cores <- pc$n
threads <- opt$threads %||% paste(unique(c(max(1L, cores %/% 2L), cores)), collapse = ",")
top <- max(as.integer(strsplit(threads, ",")[[1]]))
cases <- strsplit(opt$cases %||% paste(names(CASES), collapse = ","), ",")[[1]]
for (name in cases) case_spec(name, quick)
kernels <- opt$kernels %||% "grad,tastes,bhhh"
kern_v <- strsplit(kernels, ",")[[1]]
if (!all(kern_v %in% names(N_MATS))) stop("kernels are among ", paste(names(N_MATS), collapse = ", "))
rounds <- as.integer(opt$rounds %||% if (quick) "1" else "3")
reps <- as.integer(opt$reps %||% if (quick) "1" else "3")
out <- opt$out %||% scale_output(format(Sys.time(), "budget-%Y%m%d-%H%M%S"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out)
rscript <- file.path(R.home("bin"), "Rscript")
self <- .script_file()

# The configuration a resumed run must share. A run resumes only with
# --no-build, on the libraries already built: rebuilding into an --out that
# holds results could average two builds (another source tree, compiler or
# Makevars) without a trace.
stamp <- file.path(out, "config.rds")
old <- if (file.exists(stamp)) readRDS(stamp) else NULL
has_results <- length(list.files(out, pattern = "^res-r.*\\.rds$")) > 0L
if (has_results && !no_build) {
  stop(out, " holds results: resume them with --no-build, or use another --out")
}
if (no_build && is.null(old)) stop("--no-build needs an --out from an earlier run (", stamp, " is missing)")
config <- list(budgets = budgets, threads = threads, cases = cases, kernels = kernels, reps = reps,
               quick = quick, sizes = lapply(setNames(cases, cases), function(n) case_spec(n, quick)[c("U", "S")]))
if (!is.null(old) && !identical(old, config) && has_results) {
  stop(out, " holds results of another configuration (", stamp, "); use another --out")
}
saveRDS(config, stamp)

# Each process: the library first in R_LIBS, the caller's libraries after it
# (for choicer's dependencies); the BLAS single-threaded; no ~/.Rprofile
user_libs <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1]]
lib_of <- function(b) file.path(out, sprintf("lib-%.0f", b))
env_for <- function(b) {
  c(paste0("R_LIBS=", shQuote(paste(c(lib_of(b), user_libs[nzchar(user_libs)]),
                                    collapse = .Platform$path.sep))),
    "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1", "FLEXIBLAS_NUM_THREADS=1",
    "BLIS_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1",
    paste0("CHOICER_SCALE_OUT=", shQuote(Sys.getenv("CHOICER_SCALE_OUT", file.path(SCALE_DIR, "output")))))
}

# The machine and the build's provenance
sh <- function(cmd) tryCatch(suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE)),
                             error = function(e) character())
omp <- grep("^(OMP_|GOMP_|KMP_)", names(Sys.getenv()), value = TRUE)
info <- c(paste("date:", format(Sys.time())), paste("R:", R.version.string),
          paste("BLAS:", extSoftVersion()[["BLAS"]]), paste("LAPACK:", La_library()),
          paste("physical cores available:", cores, "(", pc$method, ") | logical CPUs:",
                parallel::detectCores(), "| nproc:", paste(sh("nproc"), collapse = " ")),
          paste("cgroup CPU quota (v2 cpu.max; v1 cfs quota and period):",
                paste(c(sh("cat /sys/fs/cgroup/cpu.max"), sh("cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us"),
                        sh("cat /sys/fs/cgroup/cpu/cpu.cfs_period_us")), collapse = " ")),
          paste("load:", loadavg()), if (length(omp)) paste0(omp, "=", Sys.getenv(omp)) else "no OMP_/GOMP_/KMP_ variables set",
          paste("numa_balancing:", paste(sh("cat /proc/sys/kernel/numa_balancing"), collapse = " ")),
          paste("budgets:", paste(budgets, collapse = ", ")), paste("threads:", threads),
          paste("cases:", paste(cases, collapse = ", ")), paste("kernels:", kernels),
          paste("rounds:", rounds, "reps:", reps, if (quick) "(quick)" else ""),
          paste("source:", if (no_build) "the libraries already in --out (see the first run's lines above)" else src),
          if (!no_build) c(
            paste("source revision:", paste(sh(paste("git -C", shQuote(src), "describe --always --dirty")), collapse = " ")),
            paste("choicer version:", read.dcf(file.path(src, "DESCRIPTION"), fields = "Version")[1, 1])),
          paste("CXX17:", paste(sh(paste(shQuote(file.path(R.home("bin"), "R")), "CMD config CXX17")), collapse = " ")),
          paste("CXX17FLAGS:", paste(sh(paste(shQuote(file.path(R.home("bin"), "R")), "CMD config CXX17FLAGS")), collapse = " ")),
          paste("user Makevars:", paste(tools::makevars_user(), collapse = " ")),
          "--- lscpu / sysctl", sh("lscpu"),
          sh(paste("sysctl machdep.cpu.brand_string hw.ncpu hw.l1dcachesize hw.l2cachesize",
                   "hw.l3cachesize hw.perflevel0.l2cachesize hw.perflevel1.l2cachesize")),
          "--- numactl -H", sh("numactl -H"))
mt <- file.path(out, "machine.txt")
if (no_build && file.exists(mt)) {
  cat(c(paste("--- resumed with --no-build:", format(Sys.time())), info), file = mt, sep = "\n", append = TRUE)
} else {
  writeLines(info, mt)
}
cat(paste(info[1:16], collapse = "\n"), "\n")

# One library per budget, each from a fresh copy of the source tree (make keeps
# object files whose sources did not change, whatever the flags). The user's
# Makevars (OpenMP flags, compilers) is copied into each build's Makevars.
build_lib <- function(b) {
  s <- file.path(out, sprintf("src-%.0f", b))
  unlink(s, recursive = TRUE); dir.create(s)
  for (x in c("DESCRIPTION", "NAMESPACE", "R", "src", "inst", "data", "man")) {
    if (file.exists(file.path(src, x))) file.copy(file.path(src, x), s, recursive = TRUE)
  }
  unlink(list.files(file.path(s, "src"), pattern = "\\.(o|so|dll)$", full.names = TRUE))
  mk <- file.path(out, sprintf("Makevars-%.0f", b))
  user <- tools::makevars_user()
  writeLines(c(if (length(user)) readLines(user[1L]),
               sprintf("CPPFLAGS += -DCHOICER_MXL_BATCH_BYTES=%.0f", b)), mk)
  lib <- lib_of(b)
  unlink(lib, recursive = TRUE); dir.create(lib)
  log <- file.path(out, sprintf("install-%.0f.log", b))
  st <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-docs", "--no-multiarch",
                "--preclean", paste0("--library=", shQuote(lib)), shQuote(s)), stdout = log, stderr = log,
                env = paste0("R_MAKEVARS_USER=", shQuote(mk)))
  if (st != 0L) stop("building the ", b, "-byte library failed; see ", log)
  unlink(s, recursive = TRUE)
}
# Namespace path, OpenMP, and the per-thread bytes that the memory check reports
# for a one-situation design (they include the budget).
probe_lib <- function(b) {
  code <- paste(
    "lib <- strsplit(Sys.getenv('R_LIBS'), .Platform$path.sep, fixed = TRUE)[[1]][1]",
    "suppressMessages(library(choicer, lib.loc = lib))",
    "stopifnot(normalizePath(dirname(getNamespaceInfo('choicer', 'path'))) == normalizePath(lib))",
    "op <- options(choicer.max_memory = 1)",
    paste0("m <- tryCatch(choicer:::mxl_loglik_gradient_parallel(rep(0, 5), matrix(c(1, 0, 0, 1), 4), ",
           "matrix(1, 4, 1), 1:4, 1L, 4L, 1, array(0, c(1, 0, 0)), 0L, TRUE, FALSE, TRUE, FALSE, ",
           "gen_seed = 1L, gen_S = 10L), error = conditionMessage)"),
    # No backslashes: R's front end passes -e code through echo, which may
    # interpret them (dash; macOS /bin/sh)
    "need <- sub(' per thread', '', regmatches(m, regexpr('[0-9.]+ [kMG]B per thread', m)))",
    "writeLines(paste('PROBE', thread_info()$openmp_enabled, if (length(need)) need else 'none'))", sep = "; ")
  r <- system2(rscript, c("--no-init-file", "-e", shQuote(code)), stdout = TRUE, stderr = TRUE,
               env = env_for(b))
  p <- grep("^PROBE ", r, value = TRUE)
  f <- if (length(p)) strsplit(trimws(p[length(p)]), " ")[[1]] else character()
  unit <- c(kB = 1e3, MB = 1e6, GB = 1e9)[if (length(f) == 4L) f[4] else NA_character_]
  # Half the last digit the message prints: 1 kB below 1 MB, 0.1 MB below 1 GB,
  # then 0.01 GB
  half <- c(kB = 500, MB = 5e4, GB = 5e6)[if (length(f) == 4L) f[4] else NA_character_]
  if (is.na(unit)) stop("probing the ", b, "-byte library failed:\n", paste(r, collapse = "\n"))
  list(openmp = identical(f[2], "TRUE"), per_thread = as.numeric(f[3]) * unname(unit),
       half = unname(half), raw = paste(r, collapse = " "))
}
if (!no_build) {
  if (!file.exists(file.path(src, "src", "mxlogit.cpp")) ||
      !any(grepl("CHOICER_MXL_BATCH_BYTES", readLines(file.path(src, "src", "mxlogit.cpp")))))
    stop(src, " is not a choicer source tree whose budget can be set (CHOICER_MXL_BATCH_BYTES)")
  for (b in budgets) { cat(sprintf("building the %.0f-byte library\n", b)); build_lib(b) }
  cc <- grep("mxlogit\\.cpp", readLines(file.path(out, sprintf("install-%.0f.log", REF))), value = TRUE)
  writeLines(c("mxlogit.cpp compile line (the 4 MiB build):", cc), file.path(out, "compile.txt"))
}
probes <- lapply(budgets, probe_lib)
need <- vapply(probes, `[[`, 0, "per_thread")
cat(sprintf("library %10.0f B: OpenMP %s, per-thread need %.0f B\n", budgets,
            vapply(probes, `[[`, NA, "openmp"), need), sep = "")
if (!all(vapply(probes, `[[`, NA, "openmp"))) {
  stop("a library was built without OpenMP; see its install log:\n",
       paste(vapply(probes, `[[`, "", "raw"), collapse = "\n"))
}
# Each library's need less its budget must equal the 4 MiB library's to within
# the two messages' rounding.
half <- vapply(probes, `[[`, 0, "half")
ref <- which(budgets == REF)
if (anyNA(need) || anyNA(half) ||
    any(abs((need - budgets) - (need[ref] - REF)) > half + half[ref])) {
  stop("the libraries' per-thread needs do not follow their budgets:\n",
       paste(vapply(probes, `[[`, "", "raw"), collapse = "\n"))
}

# The cases' rows per decision maker, the budgets whose batching differs, and
# for each budget the share of rows it batches differently from 4 MiB
units_of <- function(name) {
  cs <- case_spec(name, quick)
  f <- file.path(out, sprintf("units-%s-%d.rds", name, cs$U))
  if (!file.exists(f)) {
    code <- sprintf("SCALE_DIR <- %s; source(file.path(SCALE_DIR, 'scale_common.R')); saveRDS(unit_rows(gen_regime('%s', %dL)), %s)",
                    deparse(SCALE_DIR), cs$regime, cs$U, deparse(f))
    r <- system2(rscript, c("--no-echo", "--no-restore", "--no-init-file", "-e", shQuote(code)),
                 stdout = TRUE, stderr = TRUE, env = env_for(REF))
    if (!file.exists(f)) stop("generating the ", name, " case failed:\n", paste(r, collapse = "\n"))
  }
  readRDS(f)
}
plan <- list(); rep_of <- list(); diff_rows <- list()
for (name in cases) {
  cs <- case_spec(name, quick)
  R_u <- units_of(name)
  parts <- lapply(budgets, function(b) lapply(unique(N_MATS[kern_v]), function(n)
    batch_count(R_u, cs$S, n, b)))
  first <- vapply(seq_along(budgets), function(i)
    Position(function(p) identical(p, parts[[i]]), parts), 0L)
  ref_i <- which(budgets == REF)
  first[first == first[ref_i]] <- ref_i  # 4 MiB represents the budgets that batch like it
  run_b <- budgets[first == seq_along(budgets)]
  plan[[name]] <- run_b
  rep_of[[name]] <- setNames(budgets[first], budgets)
  diff_rows[[name]] <- do.call(rbind, lapply(kern_v, function(k) {
    ref_bc <- batch_count(R_u, cs$S, N_MATS[[k]], REF)
    data.frame(case = name, kernel = k, budget = budgets,
               diff_rows = vapply(budgets, function(b)
                 sum(R_u[batch_count(R_u, cs$S, N_MATS[[k]], b) != ref_bc]) / sum(R_u), 0))
  }))
  cat(sprintf(paste0("case %-12s %s U = %d, S = %d, max unit rows %.0f (%.1f%% of the rows; one thread's ",
                     "share at %d threads is %.1f%%): units split at 4 MiB %.4f; budgets run: %s\n"),
              name, cs$regime, cs$U, cs$S, max(R_u), 100 * max(R_u) / sum(R_u), top, 100 / top,
              mean(batch_count(R_u, cs$S, 2L, REF) > 1L), paste(run_b, collapse = ", ")))
}
diff_rows <- do.call(rbind, diff_rows)

# Rounds: every case, its budgets rotated by round and case, reversed in
# alternate rounds, so that no budget keeps its place in the order
log_dir <- file.path(out, "logs"); dir.create(log_dir, showWarnings = FALSE)
for (r in seq_len(rounds)) {
  for (ci in seq_along(cases)) {
    name <- cases[ci]
    b_r <- plan[[name]]
    b_r <- b_r[(seq_along(b_r) + r + ci - 2L) %% length(b_r) + 1L]
    if (r %% 2L == 0L) b_r <- rev(b_r)
    for (b in b_r) {
      f <- file.path(out, sprintf("res-r%d-%s-%.0f.rds", r, name, b))
      if (file.exists(f)) next
      st <- system2(rscript, c("--no-echo", "--no-restore", "--no-init-file", shQuote(self), "--worker",
                               paste0("--case=", name), paste0("--budget=", b), paste0("--round=", r),
                               paste0("--threads=", threads), paste0("--kernels=", kernels),
                               paste0("--reps=", reps), paste0("--quick=", if (quick) "1" else "0"),
                               paste0("--out=", shQuote(f))),
                    stdout = "", stderr = file.path(log_dir, sprintf("r%d-%s-%.0f.err", r, name, b)),
                    env = env_for(b))
      if (st != 0L) warning("round ", r, " ", name, " budget ", b, " failed; see logs/")
    }
  }
}

# Summary: per round, each budget's median over its reps against 4 MiB's in the
# same round; across rounds, the median of those ratios and the worst one. A
# budget not run for a case takes the times of the budget that batches it the
# same way (same_as).
files <- list.files(out, pattern = "^res-r.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no results in ", out)
raw <- do.call(rbind, lapply(files, readRDS))
expected <- do.call(rbind, lapply(cases, function(name) expand.grid(
  round = seq_len(rounds), case = name, budget = plan[[name]], stringsAsFactors = FALSE)))
got <- unique(raw[c("round", "case", "budget")])
missing <- expected[!do.call(paste, expected) %in% do.call(paste, got), ]
if (nrow(missing)) {
  warning("results missing for ", nrow(missing), " (round, case, budget): ",
          paste(apply(missing, 1L, paste, collapse = "/"), collapse = ", "))
}
raw$same_as <- NA_real_
fill <- do.call(rbind, lapply(cases, function(name) {
  m <- rep_of[[name]]
  do.call(rbind, lapply(names(m)[m != as.numeric(names(m))], function(b) {
    x <- raw[raw$case == name & raw$budget == m[[b]], ]
    if (!nrow(x)) return(NULL)
    x$same_as <- x$budget; x$budget <- as.numeric(b)
    x
  }))
}))
utils::write.csv(raw, file.path(out, "raw.csv"), row.names = FALSE)
all_rows <- rbind(raw, fill)
med <- stats::aggregate(sec ~ round + case + kernel + threads + budget, all_rows, stats::median)
ref_t <- med[med$budget == REF, c("round", "case", "kernel", "threads", "sec")]
names(ref_t)[5] <- "sec_ref"
med <- merge(med, ref_t)
med$ratio <- med$sec / med$sec_ref
chk <- stats::aggregate(check ~ case + kernel + budget, raw, `[`, 1L)
chk_ref <- chk[chk$budget == REF, c("case", "kernel", "check")]
names(chk_ref)[3] <- "check_ref"
chk <- merge(chk, chk_ref)
chk$rel_diff <- abs(chk$check - chk$check_ref) / pmax(abs(chk$check_ref), 1)
split <- stats::aggregate(cbind(split_units, split_rows, vm_add) ~ case + kernel + threads + budget,
                          all_rows, stats::median, na.action = na.pass)
summ <- stats::aggregate(cbind(sec, ratio) ~ case + kernel + threads + budget, med, stats::median)
worst <- stats::aggregate(ratio ~ case + kernel + threads + budget, med, max)
names(worst)[5] <- "ratio_worst"
sa <- unique(all_rows[c("case", "budget", "same_as")])
summ <- merge(merge(merge(merge(summ, worst), split), diff_rows, all.x = TRUE),
              chk[, c("case", "kernel", "budget", "rel_diff")], all.x = TRUE)
summ <- merge(summ, sa, all.x = TRUE)
summ <- summ[order(summ$case, summ$kernel, summ$threads, summ$budget), ]
utils::write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)
op <- options(width = 220)
print(format(transform(summ, budget_MiB = budget / 2^20, sec = round(sec, 4), ratio = round(ratio, 3),
                       ratio_worst = round(ratio_worst, 3), split_units = round(split_units, 4),
                       diff_rows = round(diff_rows, 4), vm_add_GB = round(vm_add / 1e9, 2),
                       same_as_MiB = same_as / 2^20, rel_diff = signif(rel_diff, 2))[
  , c("case", "kernel", "threads", "budget_MiB", "same_as_MiB", "sec", "ratio", "ratio_worst",
      "split_units", "diff_rows", "vm_add_GB", "rel_diff")]), row.names = FALSE)
options(op)
if (any(chk$rel_diff > 1e-10, na.rm = TRUE)) {
  warning("results differ across budgets by more than 1e-10 relative; see summary.csv (rel_diff)")
}

# The decision rule, at the largest thread count: a budget beats 4 MiB if, on
# every long panel and kernel where it batches at least a fifth of the rows
# differently, it is more than 5% faster in every round (and there is at least
# one such), and wherever else it batches at least a twentieth of the rows
# differently it is no more than 2% slower (median over rounds). Where it
# batches fewer rows differently, any effect is within the timings' noise, and
# where it batches every row as 4 MiB does its ratio is 1 by construction;
# both are left out. Every result must be present, every check value agree
# (its representative's, for a budget run once with another), and every run
# have the thread team asked.
m_top <- med[med$threads == top & med$budget != REF, ]
m_top <- merge(m_top, diff_rows, all.x = TRUE)
s_top <- summ[summ$threads == top & summ$budget != REF, ]
cand <- unique(m_top$budget)
complete <- nrow(missing) == 0L
teams_ok <- all(raw$team == raw$threads)
if (!teams_ok) warning("some runs got fewer threads than asked (OMP_THREAD_LIMIT, OMP_DYNAMIC?); see raw.csv")
long_run <- any(cases %in% LONG)
wins <- vapply(cand, function(b) {
  dec <- m_top[m_top$budget == b & m_top$case %in% LONG & m_top$diff_rows >= 0.2, ]
  rest <- s_top[s_top$budget == b & s_top$diff_rows >= 0.05 & !(s_top$case %in% LONG & s_top$diff_rows >= 0.2), ]
  # check values of the budget's own runs, or of the budget it borrows times from
  agree <- all(vapply(cases, function(name) {
    rb <- rep_of[[name]][[as.character(b)]]
    all(chk$rel_diff[chk$case == name & chk$budget == rb] <= 1e-10, na.rm = TRUE)
  }, NA))
  agree && nrow(dec) > 0 && all(dec$ratio < 0.95) && all(rest$ratio <= 1.02)
}, NA)
verdict <- if (!complete) "undecided (results missing)" else
  if (!teams_ok) "undecided (thread teams smaller than asked)" else
  if (!long_run) "undecided (no long panel among the cases)" else
  if (any(wins)) paste(cand[wins], collapse = ", ") else "none (keep 4 MiB)"
cat(sprintf(paste0("\nAt %d threads, budgets more than 5%% faster than 4 MiB in every round on the long panels ",
                   "where they batch a fifth of the rows differently, and at most 2%% slower elsewhere: %s\n"),
            top, verdict))
cat("Results in", out, "\n")
