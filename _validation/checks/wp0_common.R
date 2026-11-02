# ============================================================================
# Self-check of the suite's foundation (WP0a): the build guard, arguments,
# seeds, rep files, the output guard, metadata, threads, memory, the family
# loader and its contract, run_unit() on the T0 fixture, the fit helpers and
# at_truth()'s identity for MNL, NL and generate-mode MXL fits. Prints PASS,
# FAIL or SKIP per check and exits 1 on any FAIL. Temporary files go to
# tempdir(); the output-guard check creates and removes one directory under
# _validation/output/ (ignored by git).
#
# Usage (under 2 minutes; the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 Rscript _validation/checks/wp0_common.R
# ============================================================================

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                               value = TRUE)[1]))))
lib_dir <- file.path(root, "lib")
source(file.path(lib_dir, "families.R"), local = (boot <- new.env()))
L <- boot$vs_source_lib(lib_dir)
attach(L, name = "vs_lib", warn.conflicts = FALSE)
build <- vs_load_choicer()
vs_threads(1L)
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- file.path(tempdir(), "wp0_check")
dir.create(tmp, showWarnings = FALSE)

results <- logical()
notes <- character()
# Details printed under the check's verdict (one line per element).
note <- function(...) notes <<- c(notes, paste0("      ", strsplit(paste0(...), "\n")[[1]]))
# A check passes when `expr` is TRUE, is skipped when it is NA, and fails
# otherwise or on an error (whose message is printed).
check <- function(name, expr) {
  notes <<- character()
  res <- tryCatch(expr, error = function(e) structure(FALSE, msg = conditionMessage(e)))
  verdict <- if (identical(res, NA)) "SKIP" else if (isTRUE(res)) "PASS" else "FAIL"
  cat(sprintf("%s  %s\n", verdict, name))
  if (!is.null(attr(res, "msg"))) note("error: ", attr(res, "msg"))
  if (length(notes)) cat(notes, sep = "\n")
  results[name] <<- verdict != "FAIL"
}
# The message of the error `expr` raises ("" when it raises none).
err_msg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
# Runs an R script in a subprocess under `env`; returns its output and status.
run_r <- function(script_args, env) {
  out <- suppressWarnings(system2(rscript, script_args, stdout = TRUE, stderr = TRUE, env = env))
  list(out = out, status = attr(out, "status") %||% 0L)
}

cat(sprintf("choicer %s from %s\n\n", build$version, build$path))

# --- Build guard ------------------------------------------------------------

check("build guard: the dev build loads from R_LIBS[1]", {
  note("version ", build$version, ", lib ", build$lib)
  package_version(build$version) >= "0.2.1.9000" &&
    startsWith(getNamespaceInfo("choicer", "path"), build$lib)
})

guard_script <- file.path(tmp, "guard.R")
writeLines(c(sprintf("source(%s, local = (b <- new.env()))",
                     deparse(file.path(lib_dir, "families.R"))),
             sprintf("L <- b$vs_source_lib(%s)", deparse(lib_dir)),
             "L$vs_load_choicer()", "cat('LOADED\\n')"), guard_script)

check("build guard: a site library first (choicer 0.2.1) is refused", {
  site <- NULL
  for (lp in setdiff(.libPaths(), build$lib)) {
    d <- file.path(lp, "choicer", "DESCRIPTION")
    if (file.exists(d) && package_version(read.dcf(d, "Version")[1]) < "0.2.1.9000") {
      site <- lp
      break
    }
  }
  if (is.null(site)) {
    note("no library on .libPaths() holds an older choicer")
    NA
  } else {
    r <- run_r(shQuote(guard_script), c(paste0("R_LIBS=", site), "OMP_NUM_THREADS=1"))
    msg <- grep("older than", r$out, value = TRUE)
    note("R_LIBS=", site, ": exit ", r$status)
    note(if (length(msg)) msg[1] else paste(r$out, collapse = " | "))
    r$status != 0 && length(msg) > 0 && !any(r$out == "LOADED")
  }
})

check("build guard: R_LIBS unset is refused", {
  r <- run_r(shQuote(guard_script), c("R_LIBS=", "OMP_NUM_THREADS=1"))
  msg <- grep("R_LIBS is not set", r$out, value = TRUE)
  note("exit ", r$status, ": ", if (length(msg)) msg[1] else paste(r$out, collapse = " | "))
  r$status != 0 && length(msg) > 0
})

check("vs_require: present packages pass, a missing one is named", {
  vs_require(c("survival", "mlogit", "bayesm", "data.table", "Rcpp", "numDeriv"))
  m <- err_msg(vs_require(c("survival", "noSuchPkgWp0")))
  note(m)
  grepl("noSuchPkgWp0", m) && !grepl("survival", m)
})

# --- Arguments --------------------------------------------------------------

check("vs_args: values take their defaults' types; unknown keys list the known", {
  d <- list(N = 2e5, seed = 1L, run = NULL, verbose = FALSE, only = "all")
  a <- vs_args(c("--N=500", "--seed=7", "--run=x1", "--verbose", "--only=Q0,Q1.3"), d)
  m1 <- err_msg(vs_args("--bogus=1", d))
  m2 <- err_msg(vs_args("--seed=7.5", d))
  m3 <- err_msg(vs_args("positional", d))
  note(m1)
  identical(a, list(N = 500, seed = 7L, run = "x1", verbose = TRUE, only = "Q0,Q1.3")) &&
    grepl("--N, --seed, --run, --verbose, --only", m1, fixed = TRUE) && nzchar(m2) && nzchar(m3) &&
    identical(vs_args(character(), d), d)
})

# --- Seeds ------------------------------------------------------------------

check("seeds: FNV-1a test vectors and a known answer", {
  s <- vs_seed("T0", "200", 1:3)
  note("vs_seed(\"T0\", \"200\", 1:3) = ", paste(s, collapse = ", "))
  L$.vs_fnv1a32("a") == 0xe40c292c && L$.vs_fnv1a32("foobar") == 0xbf9cf968 &&
    L$.vs_fnv1a32("") == 0x811c9dc5 && identical(s, c(1660404319L, 1660404320L, 1660404321L))
})

check("seeds: distinct, whole and in [1, 2^31 - 2] over 10^6 reps of one cell", {
  s <- vs_seed("F3", "16k", seq_len(1e6))
  !anyDuplicated(s) && all(s == floor(s)) && min(s) >= 1 && max(s) <= 2^31 - 2 &&
    !anyNA(as.integer(s))
})

check("seeds: collisions across 20 cells x 10^4 reps", {
  cells <- expand.grid(family = paste0("F", 1:4), cell = c("500", "2k", "8k", "16k", "3.2k"),
                       stringsAsFactors = FALSE)
  s <- unlist(Map(function(f, c) vs_seed(f, c, seq_len(1e4)), cells$family, cells$cell))
  n_coll <- length(s) - length(unique(s))
  note(length(s), " seeds, ", n_coll, " collisions")
  n_coll == 0
})

check("seeds: vs_seed_rep() inverts vs_seed()", {
  identical(vs_seed_rep("F1", "8k", vs_seed("F1", "8k", 1:1000)), 1:1000)
})

# --- Rep files --------------------------------------------------------------

rec1 <- new_rep_record("T0", "200", "ml", 1, seed = 1)
out1 <- file.path(tmp, "reps")
unlink(out1, recursive = TRUE)

check("rep files: atomic write, rep_exists(), read back, no temporary left", {
  p <- rep_path(out1, "T0", "200", "ml", 1)
  before <- rep_exists(p)
  write_rep(rec1, p)
  left <- list.files(dirname(p), all.files = TRUE, no.. = TRUE)
  note(sub(paste0(out1, "/"), "", p, fixed = TRUE), "; directory holds: ",
       paste(left, collapse = ", "))
  !before && rep_exists(p) && identical(readRDS(p), rec1) && identical(left, "rep_00001.rds") &&
    basename(p) == "rep_00001.rds"
})

check("rep files: rep_path() refuses reps above 99999 and unsafe names", {
  nzchar(err_msg(rep_path(out1, "T0", "200", "ml", 100000))) &&
    nzchar(err_msg(rep_path(out1, "T0", "a/b", "ml", 1))) &&
    nzchar(err_msg(rep_path(out1, "T0", "200", "ml", 0)))
})

check("rep files: a corrupt file is reported, not fatal", {
  writeBin(as.raw(1:20), rep_path(out1, "T0", "200", "ml", 2))
  cls <- NULL
  recs <- withCallingHandlers(read_reps(out1, "T0"), warning = function(w) {
    cls <<- class(w)[1]
    invokeRestart("muffleWarning")
  })
  idx <- rep_index(out1)
  note(paste(capture.output(print(idx[c("family", "cell", "method", "rep", "ok")],
                                  row.names = FALSE)), collapse = "\n"))
  length(recs) == 1L && identical(cls, "vs_unreadable_rep") &&
    basename(attr(recs, "unreadable")) == "rep_00002.rds" &&
    identical(idx$ok, c(TRUE, NA)) && identical(idx$rep, 1:2)
})

# --- Output guard -----------------------------------------------------------

check("output guard: refuses a directory git does not ignore (_validation/lib)", {
  m <- err_msg(vs_out_dir(file.path(root, "lib"), create = FALSE))
  note(substr(m, 1, 200))
  grepl("^Refusing", m)
})

check("output guard: accepts _validation/output/<run> and creates it", {
  had_output <- dir.exists(file.path(root, "output"))
  run <- paste0(".wp0_check_", Sys.getpid())
  d <- vs_out_dir(run = run)
  ok <- dir.exists(d) && d == file.path(normalizePath(root), "output", run)
  unlink(d, recursive = TRUE)
  o <- file.path(root, "output")
  if (!had_output && !length(list.files(o, all.files = TRUE, no.. = TRUE))) {
    unlink(o, recursive = TRUE)
  }
  note(d, " (removed)")
  ok
})

check("output guard: accepts a directory outside any package tree", {
  d <- vs_out_dir(file.path(tmp, "outside", "run1"))
  dir.exists(d)
})

# --- Metadata, threads, memory ----------------------------------------------

check("metadata: provenance fields and the meta file", {
  m <- vs_metadata(ref = "wp0")
  f <- write_metadata(tmp, ref = "wp0")
  note("cpu: ", m$cpu, "; git ", substr(m$git_sha, 1, 12), " dirty ", m$git_dirty, "; ",
       basename(f))
  note("BLAS: ", m$blas, "; LAPACK: ", m$lapack)
  note("git_dirty_paths (", length(m$git_dirty_paths), "): ",
       paste(m$git_dirty_paths, collapse = ", "))
  identical(m$choicer_version, build$version) && grepl("^[0-9a-f]{40}$", m$git_sha) &&
    !is.na(m$cpu) && inherits(m$session, "sessionInfo") && file.exists(f) &&
    identical(readRDS(f)$ref, "wp0") && grepl("^meta_.+_[0-9]+\\.rds$", basename(f)) &&
    is.character(m$git_dirty_paths) && length(m$git_dirty_paths) <= 20L &&
    identical(m$git_dirty, length(m$git_dirty_paths) > 0L)
})

check("threads: vs_threads() sets choicer's threads and returns the previous count", {
  prev <- vs_threads(1L)
  is.numeric(prev) && choicer::thread_info()$omp_get_max_threads == 1L
})

check("memory: heap peak and RSS in MB", {
  m0 <- vs_mem_reset()
  x <- numeric(5e6)
  m <- vs_mem_peak(m0)
  rm(x)
  note(sprintf("heap %.1f MB, rss %.1f MB (%s)", m$heap_mb, m$rss_mb, m$rss_kind))
  m$heap_mb >= 38 && is.finite(m$rss_mb) && !is.na(m$rss_kind)
})

# --- Families ---------------------------------------------------------------

fam <- load_family("T0", L)

check("families: load_family(\"T0\") validates; fixtures listed only on request", {
  identical(fam$id, "T0") && !"T0" %in% list_families(lib_dir) &&
    "T0" %in% list_families(lib_dir, fixtures = TRUE) &&
    identical(cell_methods(fam, "800"), c("ml", "fail")) &&
    identical(family_cell(fam, "800")$N, 800L)
})

check("families: a broken copy of T0.R (R$quick names an unknown cell) is rejected", {
  bdir <- file.path(tmp, "families")
  dir.create(bdir, showWarnings = FALSE)
  src <- readLines(file.path(lib_dir, "families", "T0.R"))
  hit <- grep('quick = c("200" = 4', src, fixed = TRUE)
  src[hit] <- sub('quick = c("200" = 4, "800" = 4)', 'quick = c("200" = 4, "900" = 4)', src[hit],
                  fixed = TRUE)
  writeLines(src, file.path(bdir, "T0.R"))
  m <- err_msg(load_family("T0", L, dir = bdir))
  note(m)
  length(hit) == 1L && grepl("`R$quick`", m, fixed = TRUE) && grepl("900", m)
})

check("families: validate_family() names the offending field", {
  breaks <- list(
    "`cells$cell`"           = function(f) { f$cells$cell[2] <- "200"; f },
    "`methods$ml`"           = function(f) { f$methods$ml <- function(sim, cell) NULL; f },
    "`methods_by_cell$200`"  = function(f) { f$methods_by_cell <- list("200" = "nope"); f },
    "`gates$cell`"           = function(f) { f$gates$cell <- "16k"; f },
    "`gates`"                = function(f) { f$gates <- NULL; f },
    "`simulate`"             = function(f) { f$simulate <- function(N, seed) NULL; f },
    "`id`"                   = function(f) { f$id <- "T1"; f })
  ok <- vapply(names(breaks), function(field) {
    m <- err_msg(validate_family(breaks[[field]](fam), id = "T0"))
    hit <- grepl(field, m, fixed = TRUE)
    if (!hit) note(field, " -> ", m)
    hit
  }, logical(1))
  note(sum(ok), " of ", length(ok), " breakages named their field")
  all(ok)
})

check("families: gates name several cells and per-method metrics; bad shapes name the field", {
  multi <- fam
  multi$gates <- list(cell = c("200", "800"), metrics = c("bias_std@beta", "se_ratio:hessian"))
  breaks <- list(
    "`gates$metrics`"    = function(f) { f$gates$metrics <- list(nope = "bias_std"); f },
    "`gates$metrics$ml`" = function(f) { f$gates$metrics <- list(ml = c("bias_std", "")); f },
    "`gates$cell`"       = function(f) { f$gates$cell <- character(); f },
    "`gates$metrics` "   = function(f) { f$gates$metrics <- c("bias_std", NA); f })
  ok <- vapply(names(breaks), function(field) {
    m <- err_msg(validate_family(breaks[[field]](fam), id = "T0"))
    hit <- grepl(trimws(field), m, fixed = TRUE)
    if (!hit) note(field, " -> ", m)
    hit
  }, logical(1))
  note("T0 gates: cell ", fam$gates$cell, "; metrics$ml = ",
       paste(fam$gates$metrics$ml, collapse = ", "))
  note(sum(ok), " of ", length(ok), " bad shapes named their field")
  identical(validate_family(multi, id = "T0")$gates$cell, c("200", "800")) &&
    identical(names(fam$gates$metrics), "ml") && all(ok)
})

check("families: R_by_method (optional) is checked against the cells, methods and R$mc", {
  rbm <- function(x) { f <- fam; f$R_by_method <- x; f }
  breaks <- list(
    "`R_by_method`"      = list("900" = c(ml = 2)),
    "`R_by_method$800` " = list("800" = c(nope = 2)),
    "`R_by_method$800`"  = list("800" = c(fail = 5)),
    " `R_by_method$800`" = list("800" = c(fail = 1.5)))
  ok <- vapply(names(breaks), function(field) {
    m <- err_msg(validate_family(rbm(breaks[[field]]), id = "T0"))
    hit <- grepl(trimws(field), m, fixed = TRUE)
    note(trimws(field), " -> ", m)
    hit
  }, logical(1))
  good <- validate_family(rbm(list("800" = c(fail = 2))), id = "T0")
  identical(good$R_by_method, list("800" = c(fail = 2))) && all(ok)
})

check("families: gate cells per method and gates_quick are checked (methods, cells, runs)", {
  # fail runs only at 200 here, so fail cannot be gated at 800.
  f0 <- fam
  f0$methods_by_cell <- list("800" = "ml")
  with_gates <- function(cell, quick = NULL) {
    f <- f0
    f$gates$cell <- cell
    f$gates_quick <- quick
    f
  }
  breaks <- list(
    "`gates$cell`"       = with_gates(list(nope = "800")),
    "`gates$cell$fail`"  = with_gates(list(ml = "800", fail = "800")),
    "`gates$cell$ml`"    = with_gates(list(ml = c("800", "999"))),
    "`gates_quick$cell`" = with_gates("800", list(cell = "999", metrics = "bias_std")),
    "`gates_quick`"      = with_gates("800", list(cell = "800")))
  ok <- vapply(names(breaks), function(field) {
    m <- err_msg(validate_family(breaks[[field]], id = "T0"))
    note(field, " -> ", m)
    grepl(field, m, fixed = TRUE)
  }, logical(1))
  good <- validate_family(with_gates(list(ml = "800", fail = "200"),
                                     list(cell = list(ml = "200"), metrics = "bias_std")), id = "T0")
  is.list(good$gates$cell) && identical(good$gates_quick$cell, list(ml = "200")) && all(ok)
})

check("families: findings is function(records); gates = NULL only with findings", {
  no_gates <- fam
  no_gates["gates"] <- list(NULL)
  m1 <- err_msg(validate_family(no_gates, id = "T0"))
  bad_fn <- fam
  bad_fn$findings <- function(recs) list()
  m2 <- err_msg(validate_family(bad_fn, id = "T0"))
  fo <- no_gates
  fo$findings <- function(records) list(n = data.frame(reps = length(records)))
  m3 <- err_msg(validate_family(fo, id = "T0"))
  note("gates NULL, no findings -> ", m1)
  note("findings(recs) -> ", m2)
  note("gates NULL with findings -> ", if (nzchar(m3)) m3 else "accepted")
  grepl("`gates`", m1, fixed = TRUE) && grepl("findings-only", m1) &&
    grepl("`findings`", m2, fixed = TRUE) && !nzchar(m3)
})

check("families: pooled_groups is NULL or a function of one argument", {
  pg <- function(x) { f <- fam; f$pooled_groups <- x; f }
  m1 <- err_msg(validate_family(pg(function(records, extra) NULL), id = "T0"))
  m2 <- err_msg(validate_family(pg("terciles"), id = "T0"))
  m3 <- err_msg(validate_family(pg(function(records) c(ASC_2 = "low")), id = "T0"))
  note("two arguments -> ", m1)
  note("a string -> ", m2)
  note("function(records) -> ", if (nzchar(m3)) m3 else "accepted")
  grepl("`pooled_groups`", m1, fixed = TRUE) && grepl("`pooled_groups`", m2, fixed = TRUE) &&
    !nzchar(m3)
})

check("DGP registry: an unregistered DGP fails with the registered list", {
  m <- err_msg(get_dgp("NOPE"))
  note(m)
  grepl("^no DGP 'NOPE' registered; registered: ", m)
})

# --- run_unit() on T0 ---------------------------------------------------------

recs <- suppressMessages(lapply(c("ml", "fail"), function(m) {
  lapply(1:4, function(r) run_unit(fam, "200", m, r, ref = "wp0"))
}))
names(recs) <- c("ml", "fail")
# A record is valid when new_rep_record() accepts its fields as they stand.
valid <- function(rec) {
  isTRUE(tryCatch({ do.call(new_rep_record, rec); TRUE }, error = function(e) FALSE))
}

check("run_unit: T0 ml, reps 1-4: valid records, same seeds as vs_seed()", {
  r <- recs$ml
  note("ok: ", paste(sapply(r, `[[`, "ok"), collapse = " "), "; warnings: ",
       paste(sapply(r, function(x) nrow(x$warnings)), collapse = " "))
  all(sapply(r, valid)) && all(sapply(r, `[[`, "ok")) &&
    identical(sapply(r, `[[`, "seed"), vs_seed("T0", "200", 1:4)) &&
    all(sapply(r, function(x) nrow(x$est) == 4L && nrow(x$warnings) == 0L)) &&
    all(sapply(r, function(x) identical(x$meta$ref, "wp0")))
})

check("run_unit: T0 fail, reps 1-4: errors on even reps kept, warnings on odd reps captured", {
  r <- recs$fail
  for (x in r) {
    note("rep ", x$rep, ": ok ", x$ok, "; error: ", x$error %||% "-", "; warnings: ",
         if (nrow(x$warnings)) paste(x$warnings$class, x$warnings$stage, sep = "@") else "-")
  }
  ok <- sapply(r, `[[`, "ok")
  all(sapply(r, valid)) && identical(ok, c(TRUE, FALSE, TRUE, FALSE)) &&
    all(grepl("deliberate failure on even rep", sapply(r[!ok], `[[`, "error"))) &&
    all(sapply(r[ok], function(x) identical(x$warnings$class, "vs_fixture_warning") &&
                 identical(x$warnings$stage, "method"))) &&
    all(sapply(r[!ok], function(x) is.null(x$est) && is.finite(x$diag$unit_sec)))
})

check("run_unit: a cached sim gives the same estimates (methods share the rep's data)", {
  seed <- vs_seed("T0", "200", 1)
  sim <- fam$simulate(family_cell(fam, "200"), seed)
  r <- suppressMessages(run_unit(fam, "200", "ml", 1, sim = sim))
  identical(r$est, recs$ml[[1]]$est) && identical(r$derived, recs$ml[[1]]$derived) &&
    is.na(r$diag$sim_sec)
})

check("run_unit never throws: simulate failure, unknown cell or method, bad result", {
  f2 <- fam
  f2$simulate <- function(cell, seed) { warning("sim warns"); stop("boom") }
  f3 <- fam
  f3$methods$ml <- function(sim, cell, seed) list(est = NULL, bogus = 1)
  f4 <- fam
  f4$methods$ml <- function(sim, cell, seed) list(est = data.frame(name = "a"))
  r <- list(run_unit(f2, "200", "ml", 1), run_unit(fam, "999", "ml", 1),
            run_unit(fam, "200", "nope", 1), suppressMessages(run_unit(f3, "200", "ml", 1)),
            suppressMessages(run_unit(f4, "200", "ml", 1)), run_unit(NULL, "200", "ml", 1))
  for (x in r) note("ok ", x$ok, ": ", x$error)
  all(!sapply(r, `[[`, "ok")) && grepl("^simulate: boom", r[[1]]$error) &&
    identical(r[[1]]$warnings$stage, "simulate") && grepl("bogus", r[[4]]$error) &&
    grepl("lacks column", r[[5]]$error) && all(sapply(r[1:5], valid))
})

check("run_unit: messages muffled, timing lines dropped, the rest counted in diag$messages", {
  ml0 <- fam$methods$ml
  f5 <- fam
  f5$methods$ml <- function(sim, cell, seed) {
    message("Optimization run time 0h:0m:1s")
    message("start projected into the bounds")
    message("  start projected into the bounds  ")
    message("another note")
    ml0(sim, cell, seed)
  }
  printed <- capture.output(err <- capture.output(r <- run_unit(f5, "200", "ml", 1),
                                                  type = "message"))
  note("diag$messages: ", paste(names(r$diag$messages), r$diag$messages, sep = " = ",
                                collapse = "; "), "; printed: ", length(c(printed, err)), " lines")
  identical(r$diag$messages, c("start projected into the bounds" = 2L, "another note" = 1L)) &&
    !length(printed) && !length(err) && r$ok && identical(recs$ml[[1]]$diag$messages, integer(0))
})

check("run_unit records round-trip through rep files and rep_index()", {
  out2 <- file.path(tmp, "reps2")
  unlink(out2, recursive = TRUE)
  for (x in c(recs$ml, recs$fail)) write_rep(x, rep_path(out2, x$family, x$cell, x$method, x$rep))
  idx <- rep_index(out2)
  back <- read_reps(out2, "T0", method = "fail")
  n_ok <- tapply(idx$ok, idx$method, sum)
  note(nrow(idx), " files; ok by method: ", paste(names(n_ok), n_ok, collapse = ", "))
  nrow(idx) == 8L && sum(idx$ok) == 6L && length(back) == 4L &&
    identical(back[[2]]$error, recs$fail[[2]]$error)
})

# --- Fit helpers --------------------------------------------------------------

sim_mnl <- simulate_mnl_data(N = 400, J = 3, seed = 11, outside_option = FALSE)
fit_mnl <- suppressMessages(run_mnlogit(sim_mnl$data, "id", "alt", "choice", c("x1", "x2")))
sim_nl <- simulate_nl_data(N = 600, seed = 12)
fit_nl <- suppressMessages(run_nestlogit(sim_nl$data, "id", "j", "choice", c("X", "W"),
                                         nest_col = "nest", outside_opt_label = 0L,
                                         include_outside_option = TRUE))
sim_mxl <- simulate_mxl_data(N = 800, seed = 13, outside_option = FALSE)
fit_mxl <- suppressMessages(run_mxlogit(sim_mxl$data, "id", "alt", "choice", c("x1", "x2"),
                                        c("w1", "w2"), S = 50L, draws = "generate", seed = 7L))

check("at_truth: predict() at coef(fit), S_mult = 1, is identical (MNL, NL, MXL generate)", {
  same <- function(f) {
    identical(predict(at_truth(f, coef(f), S_mult = 1)), predict(f)) &&
      identical(predict(at_truth(f, coef(f), S_mult = 1), type = "shares"),
                predict(f, type = "shares"))
  }
  fit_store <- suppressMessages(run_mxlogit(sim_mxl$data, "id", "alt", "choice", c("x1", "x2"),
                                            c("w1", "w2"), S = 50L, rc_correlation = TRUE))
  r <- c(mnl = same(fit_mnl), nl = same(fit_nl), mxl = same(fit_mxl), mxl_store = same(fit_store),
         sigma = identical(at_truth(fit_store, coef(fit_store))$sigma, fit_store$sigma),
         lambda = identical(at_truth(fit_nl, coef(fit_nl))$lambda, fit_nl$lambda))
  note(paste(names(r), r, collapse = ", "))
  all(r)
})

check("at_truth: S_mult = 4 moves the MXL prediction slightly; another theta moves it", {
  p1 <- predict(fit_mxl)$choice_prob
  p4 <- predict(at_truth(fit_mxl, coef(fit_mxl), S_mult = 4))$choice_prob
  th <- coef(fit_mxl)
  th[["x1"]] <- th[["x1"]] + 0.5
  pt <- predict(at_truth(fit_mxl, th, S_mult = 1))$choice_prob
  d4 <- max(abs(p4 - p1))
  note(sprintf("max |P(S = 200) - P(S = 50)| = %.2e; S in draws_info: %d; theta moved: %.2e", d4,
               at_truth(fit_mxl, coef(fit_mxl))$draws_info$S, max(abs(pt - p1))))
  d4 > 0 && d4 < 0.05 && max(abs(pt - p1)) > 1e-3
})

check("at_truth: wtp() at theta0 is -theta_a / theta_p with no variance computed", {
  th <- c(x1 = 0.8, x2 = -0.6, ASC_2 = 0, ASC_3 = 0)
  a <- at_truth(fit_mnl, th)
  w <- wtp(a, price_var = "x2", attr_vars = "x1")
  nl <- at_truth(fit_nl, replace(coef(fit_nl), "Lambda_1", 0.5))
  note("wtp = ", w$Estimate, ", se = ", w$Std_Error, "; NL lambda field ",
       paste(nl$lambda, collapse = ", "))
  isTRUE(all.equal(w$Estimate, 0.8 / 0.6)) && is.na(w$Std_Error) && is.null(a$vcov) &&
    is.na(a$loglik) && nl$lambda[["Lambda_1"]] == 0.5 &&
    nzchar(err_msg(at_truth(fit_mnl, th[1:3]))) && nzchar(err_msg(at_truth(fit_mnl, c(th, z = 1))))
})

check("est_table: routes, truth by name, gating, failed route as NA plus a warning", {
  tr <- c(x1 = 0.8, x2 = -0.6, ASC_2 = 0.5, ASC_3 = -0.5)
  cls <- character()
  tab <- withCallingHandlers(est_table(fit_mnl, tr, routes = c("hessian", "robust", "cluster"),
                                       gated = c("beta", "ASC_2")),
                             warning = function(w) {
                               cls <<- c(cls, class(w)[1])
                               invokeRestart("muffleWarning")
                             })
  note(paste(capture.output(print(tab, row.names = FALSE, digits = 4)), collapse = "\n"))
  m1 <- err_msg(est_table(fit_mnl, tr[-1]))
  m2 <- err_msg(est_table(fit_mnl, c(tr, extra = 1)))
  m3 <- err_msg(est_table(fit_mnl, tr, routes = "sandwich"))
  note(m1, " | ", m2)
  identical(names(tab), names(VS_EST_COLS)) &&
    isTRUE(all.equal(tab$se_hessian, unname(fit_mnl$se))) && all(is.na(tab$se_bhhh)) &&
    all(is.finite(tab$se_robust)) && all(is.na(tab$se_cluster)) &&
    identical(cls, "vs_route_failed") && identical(tab$gated, c(TRUE, TRUE, TRUE, FALSE)) &&
    grepl("missing: x1", m1) && grepl("extra: extra", m2) && grepl("robust", m3)
})

check("vcov_common, joint_wald, lr_test", {
  vc <- vcov_common(fit_nl)
  V <- vcov(fit_nl)
  b <- c("X", "W", "Lambda_1", "Lambda_2")
  jw <- joint_wald(c(a = 1, b = 2), c(a = 0, b = 0, c = 9),
                   matrix(c(1, 0, 0, 4), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  lr <- lr_test(-100, -101.5, 1)
  identical(dim(vc$common), c(4L, 4L)) && identical(rownames(vc$common), b) &&
    isTRUE(all.equal(vc$common, V[b, b])) && identical(names(vc$asc_var), paste0("ASC_", 1:5)) &&
    isTRUE(all.equal(jw$W, 2)) && jw$df == 2L && isTRUE(all.equal(lr$LR, 3)) &&
    isTRUE(all.equal(lr$p, pchisq(3, 1, lower.tail = FALSE)))
})

check("multistart (NL) and fit_diag", {
  refit <- function(theta_init) {
    suppressMessages(run_nestlogit(sim_nl$data, "id", "j", "choice", c("X", "W"),
                                   nest_col = "nest", outside_opt_label = 0L,
                                   include_outside_option = TRUE, theta_init = theta_init))
  }
  ms <- multistart(fit_nl, refit, n = 2, seed = 5)
  d <- fit_diag(fit_nl)
  note(sprintf("dlogL %.2e, max |dtheta|/SE %.2e, failed %d; code %d, converged %s, %d iterations",
               ms$dlogL, ms$max_dtheta_se, ms$n_failed, d$code, d$converged, d$iterations))
  ms$n_failed == 0L && abs(ms$dlogL) < 1e-4 && ms$max_dtheta_se < 0.05 && isTRUE(d$converged) &&
    is.na(d$decrement)
})

fake_dgp <- c(
  "dgp_register(list(",
  "  id = 'Z9', name = 'fake logit shares', theta0 = c(a = 0.5),",
  "  simulate = function(N, seed, variant = 'default', ...) {",
  "    vs_set_seed(seed)",
  "    list(data = data.frame(y = stats::rbinom(N, 1, plogis(0.5))), truth = c(a = 0.5),",
  "         meta = list())",
  "  },",
  "  moments = function(sim) c(share_y = mean(sim$data$y), n = nrow(sim$data), untargeted = 1),",
  "  targets = data.frame(moment = c('share_y', 'n'), lo = c(0.55, 1), hi = c(0.70, 10)),",
  "  variants = 'default'))")

check("DGP registry: dgp_register() validates and names the field", {
  before <- list_dgps()                 # the lib's own DGPs, registered at source time
  ok_entry <- list(id = "Z8", name = "x", theta0 = c(a = 1),
                   simulate = function(N, seed, variant = "default") NULL,
                   moments = function(sim) NULL,
                   targets = data.frame(moment = "m", lo = 0, hi = 1), variants = "default")
  breaks <- list(
    "`variants`" = function(e) { e$variants <- "alt"; e },
    "`targets`"  = function(e) { e$targets$lo <- 2; e },
    "`theta0`"   = function(e) { e$theta0 <- c(1, 2); e },
    "`simulate`" = function(e) { e$simulate <- function(N) NULL; e },
    "missing field(s) moments" = function(e) { e$moments <- NULL; e })
  ok <- vapply(names(breaks), function(field) {
    m <- err_msg(dgp_register(breaks[[field]](ok_entry)))
    hit <- grepl(field, m, fixed = TRUE)
    if (!hit) note(field, " -> ", m)
    hit
  }, logical(1))
  note(sum(ok), " of ", length(ok), " breakages named their field")
  all(ok) && identical(list_dgps(), before)       # the failed registrations left it as it was
})

check("calibrate.R: a registered DGP prints its moments against the targets", {
  cdir <- file.path(tmp, "lib_copy")
  dir.create(cdir, showWarnings = FALSE)
  file.copy(list.files(lib_dir, pattern = "\\.R$", full.names = TRUE), cdir, overwrite = TRUE)
  writeLines(fake_dgp, file.path(cdir, "dgp_zfake.R"))
  r <- run_r(c(shQuote(file.path(cdir, "calibrate.R")), "--dgp=Z9", "--N=4000", "--seed=3"),
             c(paste0("R_LIBS=", build$lib), "OMP_NUM_THREADS=1"))
  note(paste(r$out, collapse = "\n"))
  r$status == 0L && any(grepl("^ +share_y +0\\.6[0-9]* +0\\.55 +0\\.70? +TRUE", r$out)) &&
    any(grepl("untargeted.*NA", r$out)) && any(grepl("1 of 2 targeted moments within", r$out))
})

check("calibrate.R: an unregistered DGP fails cleanly", {
  r <- run_r(c(shQuote(file.path(lib_dir, "calibrate.R")), "--dgp=NOPE"),
             c(paste0("R_LIBS=", build$lib), "OMP_NUM_THREADS=1"))
  note("exit ", r$status, ": ", paste(r$out, collapse = " | "))
  r$status == 1L && any(grepl("^calibrate: no DGP 'NOPE' registered; registered: ", r$out))
})

# --- The record's structure, for the contract ---------------------------------

cat("\nstr() of the T0 record (cell 200, method fail, rep 1):\n")
str(recs$fail[[1]], vec.len = 2, give.attr = FALSE)

unlink(tmp, recursive = TRUE)
n_fail <- sum(!results)
cat(sprintf("\n%d checks: %d passed or skipped, %d failed\n", length(results), sum(results),
            n_fail))
quit(save = "no", status = if (n_fail) 1L else 0L)
