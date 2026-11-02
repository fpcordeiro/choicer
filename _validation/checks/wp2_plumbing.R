# ============================================================================
# Self-check of the runners and the report (WP2; plan §8 "Plumbing") on the
# T0 fixture: mc.R at R = 2 per cell renders through report.R (and so do a
# run whose reps all failed, SBC ranks, a block of many constants and gated
# checks without a value); a partial run reports its absent units, their
# gate entries as na and its missing reps; run_rep()'s records equal
# run_unit()'s with its own simulation (warnings, messages and failures of
# the simulation included); two shards write exactly the units of one run;
# a shard killed mid-run resumes without duplicates and leaves its earlier
# records as they were, and a file cut short is rerun; failures and
# warnings are counted; --retry-failed reruns only the failed units;
# R_by_method runs a method at fewer paired reps (a copy of the tree whose
# T0 sets it); quick.R runs Q3's T0 mini Monte Carlo on two workers and
# report.R --tier=quick gates its rep files as quick.R did; run_shards.sh
# runs shards, stops them with their R processes, and starts none without
# --family or while a shard of the run is live; seed ranges are checked;
# argument errors stop. Prints PASS or FAIL per
# check and exits 1 on any FAIL. Outputs go to --out (default a directory
# under tempdir()), removed at the end unless --keep.
#
# Usage (under 5 minutes; the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 Rscript _validation/checks/wp2_plumbing.R
#       [--out=DIR] [--keep]
# ============================================================================

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                               value = TRUE)[1]))))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
sys.source(file.path(root, "mc.R"), envir = globalenv(), keep.source = FALSE)  # mc_seed_overlaps()
args <- vs_args(defaults = list(out = NULL, keep = FALSE))
vs_threads(1L)
t_start <- proc.time()[["elapsed"]]
tmp <- args$out %||% file.path(tempdir(), "wp2_plumbing")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
tmp <- normalizePath(tmp)
rscript <- file.path(R.home("bin"), "Rscript")
env <- c(paste0("R_LIBS=", b$build$lib), "OMP_NUM_THREADS=1", "OMP_THREAD_LIMIT=1")
fam <- load_family("T0", b$L)

results <- logical()
notes <- character()
note <- function(...) notes <<- c(notes, paste0("      ", strsplit(paste0(...), "\n")[[1]]))
check <- function(name, expr) {
  notes <<- character()
  res <- tryCatch(expr, error = function(e) structure(FALSE, msg = conditionMessage(e)))
  verdict <- if (isTRUE(res)) "PASS" else "FAIL"
  cat(sprintf("%s  %s\n", verdict, name))
  if (!is.null(attr(res, "msg"))) note("error: ", attr(res, "msg"))
  if (length(notes)) cat(notes, sep = "\n")
  results[name] <<- isTRUE(res)
}

# Runs a script of _validation/ (or any file) with the build's environment.
run <- function(script, ...) {
  out <- suppressWarnings(system2(rscript, c(shQuote(script), ...), stdout = TRUE,
                                  stderr = TRUE, env = env))
  list(out = out, status = attr(out, "status") %||% 0L)
}
mc <- function(...) run(file.path(root, "mc.R"), "--family=T0", "--fixtures",
                        paste0("--out=", shQuote(tmp)), ...)
report <- function(run_tag, ...) run(file.path(root, "report.R"),
                                     paste0("--run=", shQuote(file.path(tmp, run_tag))), ...)

# The rep files of a run, relative to it, and their md5 sums.
files_of <- function(run_tag) {
  d <- file.path(tmp, run_tag)
  rel <- list.files(d, pattern = "^rep_[0-9]{5}\\.rds$", recursive = TRUE)
  stats::setNames(unname(tools::md5sum(file.path(d, rel))), rel)
}
# The deterministic part of a record (no times, memory or meta).
core <- function(r) {
  c(r[c("family", "cell", "method", "rep", "seed", "ok", "error", "warnings", "est", "vcov_c",
        "joint_wald", "lr", "derived")],
    r$diag[c("code", "converged", "iterations", "decrement", "messages")])
}
same_runs <- function(a, b, rel = intersect(names(files_of(a)), names(files_of(b)))) {
  all(vapply(rel, function(f) {
    identical(core(readRDS(file.path(tmp, a, f))), core(readRDS(file.path(tmp, b, f))))
  }, logical(1)))
}
pid_alive <- function(pid) isTRUE(suppressWarnings(tools::pskill(pid, 0L)))
# Polls `cond` every 0.05 s until it holds or `timeout` seconds pass.
wait_for <- function(cond, timeout = 60) {
  t0 <- proc.time()[["elapsed"]]
  while (!cond()) {
    if (proc.time()[["elapsed"]] - t0 > timeout) return(FALSE)
    Sys.sleep(0.05)
  }
  TRUE
}

cat(sprintf("choicer %s from %s; metrics from %s\n\n", b$build$version, b$build$path,
            file.path(root, "lib", "metrics.R")))

# --- (1) A run at R = 2 per cell renders ------------------------------------------

check("(1) mc.R at R = 2 per cell writes 8 rep files and report.R renders them", {
  r <- mc("--reps=1:2", "--run=full")
  note(paste(utils::tail(r$out, 1), collapse = "\n"))
  f <- files_of("full")
  rp <- report("full")
  note(paste(rp$out, collapse = "\n"))
  md <- readLines(file.path(tmp, "full", "REPORT.md"))
  r$status == 0L && length(f) == 8L && rp$status == 0L &&
    all(c("## T0", "### T0 / cell 800 / method ml") %in% md) &&
    file.exists(file.path(tmp, "full", "tables", "T0_metrics.csv")) &&
    file.exists(file.path(tmp, "full", "plots", "T0_qq_800_ml.png")) &&
    file.exists(file.path(tmp, "full", "plots", "T0_ladder_ml_N.png"))
})

check("(1) report.R renders a run whose reps all failed", {
  r <- mc("--methods=fail", "--reps=1:4", "--shard=1/2", "--run=allfail")
  rp <- report("allfail")
  idx <- rep_index(file.path(tmp, "allfail"))
  md <- readLines(file.path(tmp, "allfail", "REPORT.md"))
  note(sprintf("%d records, %d ok; report exit %d", nrow(idx), sum(idx$ok), rp$status))
  r$status == 0L && nrow(idx) == 4L && !any(idx$ok) && rp$status == 0L &&
    any(grepl("failure_rate", md, fixed = TRUE))
})

check("(1) report.R draws SBC ECDF plots and lists the metrics' sbc_ecdf rows (diag$sbc)", {
  # Records of a made-up method carrying uniform ranks in diag$sbc.
  vs_with_seed(1, for (r in 1:50) {
    rec <- new_rep_record("T0", "800", "sbc", r, vs_seed("T0", "800", r),
                          diag = list(sbc = list(L = 99L, ranks = c(beta = sample(0:99, 1),
                                                                    loglik = sample(0:99, 1)))))
    write_rep(rec, rep_path(file.path(tmp, "sbc"), "T0", "800", "sbc", r))
  })
  rp <- report("sbc")
  md <- readLines(file.path(tmp, "sbc", "REPORT.md"))
  m <- utils::read.csv(file.path(tmp, "sbc", "tables", "T0_metrics.csv"))
  sb <- m[m$metric == "sbc_ecdf", c("estimand", "R", "value", "p", "verdict")]
  rows <- grep("^\\| 800 \\| sbc \\| (beta|loglik) \\| 50 \\|", md, value = TRUE)
  note(paste(utils::capture.output(print(sb, row.names = FALSE)), collapse = "\n"))
  note(paste(rows, collapse = "\n"))
  rp$status == 0L && file.exists(file.path(tmp, "sbc", "plots", "T0_sbc_800_sbc.png")) &&
    identical(sb$estimand, c("beta", "loglik")) && all(sb$R == 50L) && all(sb$p > 0.001) &&
    any(md == "### SBC ECDF tests") && length(rows) == 2L &&
    !file.exists(file.path(tmp, "sbc", "tables", "T0_sbc.csv"))
})

check("(1) run_rep()'s records equal run_unit()'s with its own simulation", {
  # Every unit of the run, against run_unit() simulating itself.
  ok1 <- all(vapply(names(files_of("full")), function(f) {
    rec <- readRDS(file.path(tmp, "full", f))
    identical(core(rec), core(run_unit(fam, rec$cell, rec$method, rec$rep)))
  }, logical(1)))
  # A simulation that warns and messages, and one that fails.
  noisy <- fam
  noisy$simulate <- function(cell, seed) {
    warning(warningCondition("sim warning", class = "vs_test_warning"))
    message("sim note")
    fam$simulate(cell, seed)
  }
  broken <- fam
  broken$simulate <- function(cell, seed) stop("no data today")
  cmp <- function(f) {
    recs <- run_rep(f, "200", 1L, c("ml", "fail"))
    all(vapply(names(recs), function(m) identical(core(recs[[m]]), core(run_unit(f, "200", m, 1L))),
               logical(1)))
  }
  rn <- run_rep(noisy, "200", 1L, "fail")$fail
  rb <- run_rep(broken, "200", 1L, "ml")$ml
  note("noisy rep: warnings ", paste(rn$warnings$stage, rn$warnings$class, collapse = ", "),
       "; messages ", paste(names(rn$diag$messages), rn$diag$messages, collapse = ", "))
  note("broken rep: ok ", rb$ok, ", error '", rb$error, "'")
  ok1 && cmp(noisy) && cmp(broken) && !rb$ok && rb$error == "simulate: no data today" &&
    identical(rn$warnings$stage, c("simulate", "method"))
})

check("(1) report.R summarizes a block of many constants; every row stays in the CSV", {
  # 40 reps of a made-up method with 2 slopes and 15 constants.
  nm <- c("b1", "b2", sprintf("ASC_%d", 2:16))
  vs_with_seed(2, for (r in 1:40) {
    est <- data.frame(name = nm, block = c("beta", "beta", rep("asc", 15)),
                      est = stats::rnorm(17, 0, 0.1), truth = 0, se_hessian = 0.1,
                      se_bhhh = NA_real_, se_robust = NA_real_, se_cluster = NA_real_,
                      gated = TRUE, stringsAsFactors = FALSE)
    write_rep(new_rep_record("T0", "200", "many", r, vs_seed("T0", "200", r), est = est),
              rep_path(file.path(tmp, "many"), "T0", "200", "many", r))
  })
  rp <- report("many")
  md <- readLines(file.path(tmp, "many", "REPORT.md"))
  tab <- utils::read.csv(file.path(tmp, "many", "tables", "T0_metrics.csv"))
  blk <- grep("^Block asc, 15 estimands", md, value = TRUE)
  pooled <- grep("^\\| pooled_", md, value = TRUE)
  note(blk, "\n", paste(pooled, collapse = "\n"))
  rp$status == 0L && length(blk) == 1L && length(pooled) == 3L &&
    any(md == "Its worst rows:") && !any(grepl("^\\| ASC_[0-9]+ \\| hessian", md)) &&
    any(grepl("^\\| b1 \\| hessian", md)) &&
    setequal(unique(tab$estimand[tab$block %in% "asc" & tab$metric == "bias_std"]), nm[-(1:2)])
})

# Every gated row of 800/ml (T0's gate cell and gated method) is na, at
# least T0's four gate entries, and REPORT.md and stdout say so in bold.
na_gated_lines <- function(rp, run_tag) {
  m <- utils::read.csv(file.path(tmp, run_tag, "tables", "T0_metrics.csv"))
  g <- m[m$gated, ]
  k <- nrow(g)
  md <- readLines(file.path(tmp, run_tag, "REPORT.md"))
  hits <- c(grep(sprintf("^\\*\\*%d gated check\\(s\\) without a value", k), md, value = TRUE),
            grep(sprintf("^\\*\\*Gated rows without a value \\(na\\): T0 \\(%d\\)", k), md,
                 value = TRUE))
  note(paste(c(grep("T0:", rp$out, value = TRUE), hits,
               paste("gated:", paste(g$metric, g$route, collapse = ", "))), collapse = "\n"))
  rp$status == 0L && k >= 4L && all(g$verdict == "na") && all(g$cell == 800L & g$method == "ml") &&
    all(c("bias_std", "se_ratio", "tail_lower", "tail_upper") %in% g$metric) &&
    length(hits) == 2L &&
    any(startsWith(rp$out, sprintf("  T0: %d gated rows, 0 pass, 0 fail, %d na", k, k)))
}

check("(1) gated checks without a value (no row: all reps failed) are reported as na", {
  for (r in 1:2) {
    write_rep(new_rep_record("T0", "800", "ml", r, vs_seed("T0", "800", r), ok = FALSE,
                             error = "made-up failure"),
              rep_path(file.path(tmp, "nagate"), "T0", "800", "ml", r))
  }
  rp <- report("nagate")
  na_gated_lines(rp, "nagate")
})

check("(1, 2) a partial run: absent units in bold, their gate entries na, missing reps counted", {
  # Cell 200, method ml, reps 1-2 only: T0's gate cell 800 has no rep file.
  r <- mc("--cells=200", "--methods=ml", "--reps=1:2", "--run=partial")
  rp <- report("partial")
  md <- readLines(file.path(tmp, "partial", "REPORT.md"))
  u <- utils::read.csv(file.path(tmp, "partial", "tables", "T0_units.csv"))
  s <- utils::read.csv(file.path(tmp, "partial", "tables", "summary.csv"))
  miss <- grep("^\\*\\*Reps missing against the design", md, value = TRUE)
  na <- na_gated_lines(rp, "partial")
  note(miss)
  note(paste(utils::capture.output(print(u[c("cell", "method", "status", "files", "expected",
                                             "missing")], row.names = FALSE)), collapse = "\n"))
  # 4 units x R$mc = 4 expected reps; 2 present, the 2 at the end of 200/ml's range missing too.
  r$status == 0L && rp$status == 0L && na &&
    sum(u$status == "**absent**") == 3L && s$missing == 14L && s$absent == 3L &&
    identical(miss, "**Reps missing against the design (an incomplete run): T0 (14).**") &&
    any(grepl("| **14** |", md, fixed = TRUE)) &&
    any(md == "**Absent from the run: no rep file.**")
})

check("(3) resume reruns a rep file that does not read (cut short)", {
  r1 <- mc("--cells=800", "--reps=1:1", "--run=cut")
  f <- file.path(tmp, "cut", "T0", "800", "ml", "rep_00001.rds")
  bytes <- readBin(f, "raw", file.size(f))
  writeBin(bytes[seq_len(length(bytes) %/% 2)], f)
  cut <- is.null(.vs_read_rep(f))
  r2 <- mc("--cells=800", "--reps=1:1", "--run=cut")
  note(paste(grep("^T0 |^done", r2$out, value = TRUE), collapse = "\n"))
  r1$status == 0L && cut && r2$status == 0L &&
    any(r2$out == "T0 800 rep 1: unreadable rep file(s) of ml rerun") &&
    identical(core(readRDS(f)), core(readRDS(file.path(tmp, "full", "T0", "800", "ml",
                                                       "rep_00001.rds"))))
})

# --- (2) Shards -------------------------------------------------------------------

check("(2) shards 1/2 and 2/2 write exactly the units of the unsharded run", {
  r1 <- mc("--reps=1:2", "--run=sharded", "--shard=1/2")
  r2 <- mc("--reps=1:2", "--run=sharded", "--shard=2/2")
  a <- files_of("full")
  s <- files_of("sharded")
  seeds <- function(t, rel) vapply(rel, function(f) readRDS(file.path(tmp, t, f))$seed, 1)
  note("shard 1/2: ", paste(grep("^T0 ", r1$out, value = TRUE), collapse = "; "))
  note("shard 2/2: ", paste(grep("^T0 ", r2$out, value = TRUE), collapse = "; "))
  r1$status == 0L && r2$status == 0L && setequal(names(a), names(s)) &&
    identical(seeds("full", names(a)), seeds("sharded", names(a))) && same_runs("full", "sharded")
})

# --- (3) Killed and resumed ----------------------------------------------------------

check("(3) a shard killed mid-run resumes without duplicates; earlier records unchanged", {
  cmd <- c(env, shQuote(rscript), shQuote(file.path(root, "mc.R")),
           "--family=T0", "--fixtures", "--cells=800", "--reps=1:200", "--shard=1/2",
           "--run=killed", paste0("--out=", shQuote(tmp)))
  pid <- as.integer(system(paste(paste(cmd, collapse = " "), ">", shQuote(file.path(tmp, "killed.log")),
                                 "2>&1 & echo $!"), intern = TRUE))
  appeared <- wait_for(function() length(files_of("killed")) > 0L)
  tools::pskill(pid)
  gone <- wait_for(function() !pid_alive(pid), 30)
  before <- files_of("killed")
  tmpfiles <- list.files(file.path(tmp, "killed"), pattern = "\\.tmp$", recursive = TRUE,
                         all.files = TRUE)
  res <- mc("--cells=800", "--reps=1:200", "--shard=1/2", "--run=killed")
  after <- files_of("killed")
  idx <- rep_index(file.path(tmp, "killed"))
  expected <- as.vector(outer(c("fail", "ml"), seq(2, 200, 2), function(m, r) {
    sprintf("T0/800/%s/rep_%05d.rds", m, r)
  }))
  # The records written before the kill against the same units computed here.
  earlier <- vapply(names(before), function(f) {
    rec <- readRDS(file.path(tmp, "killed", f))
    identical(core(rec), core(run_rep(fam, rec$cell, rec$rep, rec$method)[[1]]))
  }, logical(1))
  done <- grep("^done:", res$out, value = TRUE)
  note(sprintf("pid %d killed after %d rep file(s) (stopped: %s; %d temporary file(s) left)",
               pid, length(before), gone, length(tmpfiles)))
  note("resume: ", done)
  appeared && gone && length(before) < length(expected) && res$status == 0L &&
    setequal(names(after), expected) && length(after) == length(expected) &&
    identical(after[names(before)], before) && all(earlier) && !anyNA(idx$ok) &&
    grepl(sprintf("%d skipped", length(before)), done)
})

# --- (4) Failures and warnings counted ----------------------------------------------

check("(4) the fail method's failures enter failure_rate and its warnings are counted by class", {
  tab <- utils::read.csv(file.path(tmp, "full", "tables", "T0_metrics.csv"))
  fr <- tab[tab$metric == "failure_rate", c("cell", "method", "value", "n_failed")]
  w <- utils::read.csv(file.path(tmp, "full", "tables", "T0_warnings.csv"))
  fl <- utils::read.csv(file.path(tmp, "full", "tables", "T0_failures.csv"))
  note(paste(utils::capture.output(print(fr, row.names = FALSE)), collapse = "\n"))
  note(paste(utils::capture.output(print(w[c("cell", "method", "stage", "class", "count")],
                                         row.names = FALSE)), collapse = "\n"))
  note(paste(utils::capture.output(print(fl, row.names = FALSE)), collapse = "\n"))
  f <- fr[fr$method == "fail", ]
  m <- fr[fr$method == "ml", ]
  nrow(f) == 2L && all(f$value == 0.5) && all(f$n_failed == 1L) && all(m$value == 0) &&
    nrow(w) == 2L && all(w$class == "vs_fixture_warning") && all(w$count == 1L) &&
    all(w$method == "fail") && nrow(fl) == 2L && all(fl$count == 1L) &&
    all(fl$error == "T0 fixture: deliberate failure on even rep 2")
})

# --- (5) Retrying failed units -------------------------------------------------------

check("(5) --retry-failed reruns only the failed units", {
  before <- files_of("full")
  r <- mc("--reps=1:2", "--run=full", "--retry-failed")
  after <- files_of("full")
  lines <- grep("^T0 ", r$out, value = TRUE)
  failed <- c("T0/200/fail/rep_00002.rds", "T0/800/fail/rep_00002.rds")
  changed <- names(after)[after != before[names(after)]]
  note(paste(lines, collapse = "\n"))
  note("rewritten: ", paste(changed, collapse = ", "))
  r$status == 0L && setequal(names(after), names(before)) && setequal(changed, failed) &&
    identical(sub(" [(][0-9.]+ s[)]$", "", lines),
              c("T0 200 rep 2: fail FAIL", "T0 800 rep 2: fail FAIL"))
})

# --- (5') Methods at different R (R_by_method) ---------------------------------------

check("(5') R_by_method runs a method at fewer reps on the same paired reps; --reps overrides", {
  # A copy of the scripts and the library whose T0 runs fail at 2 reps in cell 800.
  tree <- file.path(tmp, "tree", "_validation")
  dir.create(tree, recursive = TRUE)
  file.copy(file.path(root, c("mc.R", "report.R", "lib")), tree, recursive = TRUE)
  cat('\nfamily$R_by_method <- list("800" = c(fail = 2))\n',
      file = file.path(tree, "lib", "families", "T0.R"), append = TRUE)
  m1 <- run(file.path(tree, "mc.R"), "--family=T0", "--fixtures", "--run=rbm",
            paste0("--out=", shQuote(tmp)))
  m2 <- run(file.path(tree, "mc.R"), "--family=T0", "--fixtures", "--run=rbm3", "--reps=1:3",
            paste0("--out=", shQuote(tmp)))
  rp <- run(file.path(tree, "report.R"), paste0("--run=", shQuote(file.path(tmp, "rbm"))))
  u <- utils::read.csv(file.path(tmp, "rbm", "tables", "T0_units.csv"))
  f <- names(files_of("rbm"))
  n <- table(sub("/rep_.*$", "", f))
  note(paste(names(n), n, sep = ": ", collapse = "; "))
  note(paste(grep("^T0 800", m1$out, value = TRUE), collapse = "\n"))
  m1$status == 0L && m2$status == 0L && rp$status == 0L &&
    identical(as.integer(n[c("T0/200/fail", "T0/200/ml", "T0/800/fail", "T0/800/ml")]),
              c(4L, 4L, 2L, 4L)) &&
    u$expected[u$cell == 800 & u$method == "fail"] == 2L &&
    all(u$expected[!(u$cell == 800 & u$method == "fail")] == 4L) &&
    length(files_of("rbm3")) == 12L &&
    same_runs("full", "rbm", rel = names(files_of("full")))
})

# --- (6) quick.R's Q3 ---------------------------------------------------------------

check("(6) quick.R --only=Q3:T0 --fixtures --workers=2 runs T0's mini Monte Carlo", {
  q <- run(file.path(root, "quick.R"), "--only=Q3:T0", "--fixtures", "--workers=2",
           paste0("--out=", shQuote(file.path(tmp, "quick"))))
  note(paste(q$out, collapse = "\n"))
  csv <- list.files(file.path(tmp, "quick"), pattern = "_Q3_T0\\.csv$", full.names = TRUE)
  g <- if (length(csv)) utils::read.csv(csv[1]) else NULL
  g <- g[g$gated, ]
  reps_dir <- list.files(file.path(tmp, "quick"), pattern = "_reps$", full.names = TRUE)[1]
  reps <- rep_index(reps_dir)
  # T0's gate entries: 4 metrics x 2 slopes x 2 cells (route_failure_rate rows besides).
  est_rows <- g$metric %in% c("bias_std", "se_ratio", "tail_lower", "tail_upper")
  any(grepl("^Q3:T0 gate table", q$out)) && !is.null(g) && sum(est_rows) == 16L &&
    nrow(reps) == 16L && q$status == as.integer(any(g$verdict %in% c("fail", "na")))
})

check("(6) report.R --tier=quick on Q3's rep files gates as quick.R did (cells, tolerances)", {
  reps_dir <- list.files(file.path(tmp, "quick"), pattern = "_reps$", full.names = TRUE)[1]
  q <- utils::read.csv(list.files(file.path(tmp, "quick"), pattern = "_Q3_T0\\.csv$",
                                  full.names = TRUE)[1])
  rp <- run(file.path(root, "report.R"), paste0("--run=", shQuote(reps_dir)), "--tier=quick")
  m <- utils::read.csv(file.path(reps_dir, "tables", "T0_metrics.csv"))
  s <- utils::read.csv(file.path(reps_dir, "tables", "summary.csv"))
  key <- function(d) paste(d$cell, d$method, d$metric, d$estimand, d$route, d$verdict)
  note(grep("T0:", rp$out, value = TRUE))
  note("gated cells: report ", paste(sort(unique(m$cell[m$gated])), collapse = ", "),
       "; quick.R ", paste(sort(unique(q$cell[q$gated])), collapse = ", "))
  rp$status == 0L && setequal(key(m[m$gated, ]), key(q[q$gated, ])) &&
    setequal(m$cell[m$gated], c(200L, 800L)) && s$missing == 0L && s$absent == 0L &&
    any(readLines(file.path(reps_dir, "REPORT.md")) == "### Units (tier quick)")
})

# --- (7) run_shards.sh ----------------------------------------------------------------

shards <- function(...) {
  out <- suppressWarnings(system2("sh", c(shQuote(file.path(root, "run_shards.sh")), "--family=T0",
                                          "--fixtures", paste0("--out=", shQuote(tmp)), ...),
                                  stdout = TRUE, stderr = TRUE, env = env))
  list(out = out, status = attr(out, "status") %||% 0L)
}
shard_pids <- function(run_tag) {
  p <- list.files(file.path(tmp, run_tag, "logs"), pattern = "\\.pid$", full.names = TRUE)
  stats::setNames(as.integer(vapply(p, readLines, "")), basename(p))
}

check("(7) run_shards.sh runs two shards that write the units of the unsharded run", {
  s <- shards("--reps=1:2", "--shards=2", "--run=viash")
  pids <- shard_pids("viash")
  done <- wait_for(function() !any(vapply(pids, pid_alive, logical(1))), 120)
  logs <- list.files(file.path(tmp, "viash", "logs"), pattern = "\\.log$", full.names = TRUE)
  tails <- vapply(logs, function(l) utils::tail(readLines(l), 1), "")
  note(paste(grep("^(shard|run directory)", s$out, value = TRUE), collapse = "\n"))
  note(paste(tails, collapse = "\n"))
  s$status == 0L && length(pids) == 2L && done && all(grepl("exit status 0", tails)) &&
    setequal(names(files_of("viash")), names(files_of("full"))) && same_runs("full", "viash")
})

check("(7) stop_shards.sh stops the shards and their R processes", {
  s <- shards("--cells=800", "--reps=1:400", "--shards=2", "--run=stopped")
  pids <- shard_pids("stopped")
  started <- wait_for(function() length(files_of("stopped")) > 0L)
  kids <- as.integer(unlist(lapply(pids, function(p) {
    suppressWarnings(system2("pgrep", c("-P", p), stdout = TRUE))
  })))
  st <- suppressWarnings(system2("sh", shQuote(file.path(tmp, "stopped", "logs", "stop_shards.sh")),
                                 stdout = TRUE, stderr = TRUE))
  gone <- wait_for(function() !any(vapply(c(pids, kids), pid_alive, logical(1))), 30)
  note("wrappers ", paste(pids, collapse = ", "), "; R processes ", paste(kids, collapse = ", "))
  note(paste(st, collapse = "\n"))
  s$status == 0L && started && length(kids) == 2L && length(st) == 2L && gone &&
    length(files_of("stopped")) < 800L
})

check("(7) run_shards.sh starts no shard without --family or while a shard of the run is live", {
  # No --family: refused before the warm pass.
  nf <- suppressWarnings(system2("sh", c(shQuote(file.path(root, "run_shards.sh")), "--fixtures",
                                         "--shards=2", "--run=nofam",
                                         paste0("--out=", shQuote(tmp))),
                                 stdout = TRUE, stderr = TRUE, env = env))
  # Shard 2's pid file names a live process (this one): shard 1 must not start either.
  dir.create(file.path(tmp, "live", "logs"), recursive = TRUE)
  writeLines(as.character(Sys.getpid()), file.path(tmp, "live", "logs", "shard_2_of_2.pid"))
  lv <- shards("--reps=1:2", "--shards=2", "--run=live")
  note(paste(grep("run_shards.sh:", c(nf, lv$out), value = TRUE), collapse = "\n"))
  !identical(attr(nf, "status") %||% 0L, 0L) && any(grepl("--family=F is required", nf)) &&
    !dir.exists(file.path(tmp, "nofam")) && lv$status != 0L &&
    any(grepl("running already .* no shard started", lv$out)) &&
    !file.exists(file.path(tmp, "live", "logs", "shard_1_of_2.pid")) && !length(files_of("live"))
})

# --- (8) Seeds and arguments ----------------------------------------------------------

check("(8) --check-seeds: T0 has no overlap; overlapping ranges are found exactly", {
  r <- mc("--check-seeds")
  M <- 2^31 - 2
  D <- (vs_seed("T0", "800", 1) - vs_seed("T0", "200", 1)) %% M
  fake <- function(n200) {
    f <- fam
    f$R$mc <- c("200" = n200, "800" = 1)
    list(T0 = f)
  }
  near <- mc_seed_overlaps(fake(D + 1))$overlaps
  apart <- mc_seed_overlaps(fake(D))$overlaps
  note(paste(r$out, collapse = "\n"))
  note(sprintf("cell 800's first seed lies %.0f seeds after cell 200's", D))
  r$status == 0L && any(r$out == "no overlapping seed ranges") && NROW(near) == 1L &&
    is.null(apart)
})

check("(8) argument errors stop before any unit runs", {
  a <- run(file.path(root, "mc.R"), "--family=T0", paste0("--out=", shQuote(tmp)), "--run=x")
  s <- run(file.path(root, "mc.R"), "--family=T0", "--fixtures", "--shard=1/2",
           paste0("--out=", shQuote(tmp)))
  u <- mc("--cells=999", "--run=x")
  q <- run(file.path(root, "quick.R"), "--only=Q3:T0")
  msg <- function(r) grep("Error", r$out, value = TRUE)[1]
  for (r in list(a, s, u, q)) note(msg(r))
  a$status != 0L && grepl("plumbing fixture", msg(a)) && s$status != 0L &&
    grepl("--shard needs --run", msg(s), fixed = TRUE) && u$status != 0L &&
    grepl("unknown cell", msg(u)) && q$status != 0L && grepl("add --fixtures", msg(q)) &&
    !dir.exists(file.path(tmp, "x", "T0"))
})

cat("\nREPORT.md of the R = 2 run (excerpt):\n")
md <- readLines(file.path(tmp, "full", "REPORT.md"))
cat(md[seq_len(min(length(md), which(md == "### T0 / cell 200 / method fail")[1] - 1L))], sep = "\n")

if (!args$keep) unlink(tmp, recursive = TRUE)
n_fail <- sum(!results)
cat(sprintf("\n%d checks: %d passed, %d failed (%.0f s)\n", length(results), sum(results), n_fail,
            proc.time()[["elapsed"]] - t_start))
quit(save = "no", status = if (n_fail) 1L else 0L)
