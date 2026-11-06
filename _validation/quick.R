# ============================================================================
# The quick tier (plan §5; 15 minutes at 8 workers): a registry of tests
# Q0 ... Q4, each a function(ctx, entries) returning check rows
# data.frame(test, check, value, threshold, pass), pass TRUE or FALSE. A row
# without a threshold and with pass NA is information, printed and kept apart; a
# row with a threshold and pass NA (a comparison that met a NaN) fails.
#   Q0    every DGP against choicer at theta0 (lib/quick_checks.R)
#   Q1.1  MNL on mode_choice against survival::clogit and mlogit (ref_q11())
#   Q1.2  NL (GROUND nest) against mlogit (ref_q12())
#   Q1.3  MXL, S = 4000, against exact Gauss-Hermite ML (ref_q13(), about 110 s)
#   Q1.4  MNP against bayesm::rmnpGibbs (ref_q14(), about 50 s)
#   Q2    one large-N fit per family F1-F4, F6, on workers (lib/quick_checks.R)
#   Q3    a mini Monte Carlo of the bottom rungs of F1-F4 and F8 on workers
#         through run_rep() (the long tier's paired reps and seeds), then
#         metrics_cell(tol_set = "quick") and gate_quick() (lib/metrics.R), gated
#         by the family's gates_quick (the checks theory supports at that rung,
#         confirmed at R = 300), else at the entry's cells with the family's
#         gate metrics; every entry also bounds its failed reps, and reports how
#         often the default start failed or lost (an information row)
#   Q4    the mini-SBC of MNP, HMNL, HMNP on workers: q4_check() (lib/dgp_bayes.R),
#         through gate_quick() too, at F9's alpha; its failed reps bounded alike
# A reference that is missing or fails is a failed check, printed at once and
# in the summary; a test without a check row fails too.
#
# Usage (the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/quick.R
#       [--only=Q0,Q1.3,Q2:F1,Q3:F2] [--workers=8] [--out=DIR] [--ref=NAME]
#       [--reps=N] [--perturb=x] [--fixtures]
#   --only      the tests to run (default all); Q2:<entry> and Q3:<entry> run
#               one of their entries
#   --reps      every cell of Q3 and Q4 at N reps instead of its R$quick (the
#               checks and the budget's timing at R = 2; the gates hold at 100)
#   --perturb   adds x to theta0 on choicer's side of Q0: a failure injection
#               for checks/wp9_quick.R
#   --fixtures  adds Q3's T0 plumbing entry (never run by default)
# Writes into --out (default _validation/output/quick/): quick_<stamp>.csv
# (every check row, with the test's seconds and `run`, what the run was:
# q_run_kind(), "tier" only for the whole tier unshortened), quick_<stamp>_info.csv
# (the information rows), quick_<stamp>_Q3_<entry>.csv and quick_<stamp>_Q4.csv
# (the metrics), and the reps of Q3 and Q4 under quick_<stamp>_reps/ with
# q3_plan.rds, the units they ran, Q4's included. Render them with this script's
# verdicts (both sides call gate_quick()) by
#   Rscript _validation/report.R --run=<out>/quick_<stamp>_reps --tier=quick
# Prints a summary and the gate tables; exits 1 if any check fails.
# ============================================================================

# --- Q3's entries and Q4's ----------------------------------------------------------

# An entry: a family and, optionally, `cells`, `methods` and `R` (default: the
# cells of the family's R$quick with their reps, every method of the cell).
Q3_ENTRIES <- list(
  F1 = list(family = "F1"),
  F2 = list(family = "F2"),
  F3 = list(family = "F3"),
  F4 = list(family = "F4"),
  F8 = list(family = "F8", methods = "mnl_unweighted")   # plan: F8-MNL, under quotas
)
Q3_FIXTURES <- list(T0 = list(family = "T0"))
Q4_ENTRY <- list(F9 = list(family = "F9"))

q3_entries <- function(fixtures) c(Q3_ENTRIES, if (fixtures) Q3_FIXTURES)

# --- On a worker --------------------------------------------------------------------

# One rep of one (family, cell), every method of the entry (families loaded
# once per worker). Errors come back as a string.
q_rep_task <- function(task) {
  tryCatch({
    fam <- get0(task$family, envir = .q_fams, inherits = FALSE)
    if (is.null(fam)) {
      fam <- load_family(task$family, .q_L)
      assign(task$family, fam, envir = .q_fams)
    }
    run_rep(fam, task$cell, task$rep, task$methods, ref = task$ref)
  }, error = function(e) conditionMessage(e))
}

# One unit of Q2; an error is its failed row.
q2_task <- function(k) {
  tryCatch(q2_unit(k, .q_L), error = function(e) {
    q_row("Q2", sprintf("%s: %s", k, conditionMessage(e)), pass = FALSE)
  })
}

# Runs on each worker: the master's R_LIBS, the build guard and the library
# (vs_boot()), the same choicer as the master's, one thread.
q_worker_init <- function(root, rlibs, path) {
  Sys.setenv(R_LIBS = rlibs, OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  .libPaths(c(strsplit(rlibs, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
  sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
  b <- boot$vs_boot(root)
  if (getNamespaceInfo("choicer", "path") != path) {
    stop("a worker loaded choicer from ", getNamespaceInfo("choicer", "path"), ", not ", path)
  }
  assign(".q_L", b$L, envir = globalenv())
  assign(".q_fams", new.env(), envir = globalenv())
  vs_threads(1L)
  invisible(NULL)
}

# --- The mini Monte Carlos ----------------------------------------------------------

# Runs the reps of `ents` on the workers (largest N first, one per dispatch:
# parLapplyLB's default sends two chunks per worker), writes their rep files
# and adds the units to q3_plan.rds (report.R's name for it; Q4's units are there
# too). Returns list(plan = entry -> (fam, cells), recs, rows): a family that
# cannot be loaded and a worker's error are failed rows.
q_mc_run <- function(ctx, test, ents) {
  rows <- list()
  plan <- list()
  tasks <- list()
  units <- list()
  for (k in names(ents)) {
    e <- ents[[k]]
    id <- paste0(test, ":", k)
    fam <- tryCatch(load_family(e$family, ctx$L), error = function(err) err)
    if (inherits(fam, "error")) {
      rows[[k]] <- q_row(id, paste("not run:", conditionMessage(fam)), pass = FALSE)
      next
    }
    rq <- fam$R$quick
    cells <- e$cells %||% names(rq)[rq > 0]
    if (!length(cells)) {
      rows[[k]] <- q_row(id, "not run: the family has no quick cells (R$quick)", pass = FALSE)
      next
    }
    for (cl in cells) {
      ms <- intersect(cell_methods(fam, cl), e$methods %||% names(fam$methods))
      R <- if (ctx$reps > 0L) ctx$reps else e$R %||% rq[[cl]]
      if (length(ms)) {
        units[[length(units) + 1L]] <- data.frame(family = fam$id, cell = cl, method = ms,
                                                  R = as.integer(R), stringsAsFactors = FALSE)
      }
      for (r in seq_len(R)) {
        tasks[[length(tasks) + 1L]] <- list(entry = k, family = fam$id, cell = cl, rep = r,
                                            methods = ms, ref = ctx$ref,
                                            N = fam$cells$N[match(cl, fam$cells$cell)])
      }
    }
    plan[[k]] <- list(fam = fam, cells = cells)
  }
  if (!length(tasks)) return(list(plan = plan, recs = list(), rows = rows))
  tasks <- tasks[order(-vapply(tasks, `[[`, 1, "N"))]
  cat(sprintf("%s: %d reps of %s on %d workers\n", test, length(tasks),
              paste(names(plan), collapse = ", "), min(ctx$workers, length(tasks))))
  flush(stdout())
  res <- parallel::parLapplyLB(ctx$cluster(length(tasks)), tasks, q_rep_task, chunk.size = 1L)
  dir.create(ctx$reps_dir, showWarnings = FALSE, recursive = TRUE)
  write_metadata(ctx$reps_dir, ctx$ref)
  plan_file <- file.path(ctx$reps_dir, "q3_plan.rds")
  saveRDS(rbind(if (file.exists(plan_file)) readRDS(plan_file), do.call(rbind, units)), plan_file)
  recs <- list()
  for (i in seq_along(tasks)) {
    if (is.character(res[[i]])) {
      rows[[length(rows) + 1L]] <- q_row(paste0(test, ":", tasks[[i]]$entry),
                                         sprintf("worker error at %s/%s rep %d: %s",
                                                 tasks[[i]]$family, tasks[[i]]$cell,
                                                 tasks[[i]]$rep, res[[i]]), pass = FALSE)
      next
    }
    for (rec in res[[i]]) {
      write_rep(rec, rep_path(ctx$reps_dir, rec$family, rec$cell, rec$method, rec$rep))
      recs[[length(recs) + 1L]] <- rec
    }
  }
  list(plan = plan, recs = recs, rows = rows)
}

# One row per gated failure_rate row of gate_quick()'s table `tab`, with its
# verdict: a (cell, method)'s share of failed reps (an error, or a fit no start
# certified) against the quick tier's tolerance on failures, 2%, whatever the
# surviving reps say: metrics read those alone; the message of the first failure
# when there is one. At R = 300 no entry failed a rep. T0's `fail` method fails
# half on purpose, so the fixture has no such row and gets no check here.
q_failed_rows <- function(id, recs, tab) {
  f <- tab[tab$metric == "failure_rate" & tab$gated, , drop = FALSE]
  key <- vapply(recs, function(r) paste0(r$cell, "/", r$method), "")
  do.call(rbind, lapply(seq_len(nrow(f)), function(i) {
    rr <- recs[key == paste0(f$cell[i], "/", f$method[i])]
    bad <- !vapply(rr, function(r) isTRUE(r$ok), logical(1))
    q_row(id, sprintf("%s/%s failed reps: %d of %d%s", f$cell[i], f$method[i], sum(bad), length(rr),
                      if (any(bad)) paste0(" (", rr[bad][[1L]]$error, ")") else ""),
          f$value[i], f$tol[i], f$verdict[i] == "pass")
  }))
}

# One information row per (cell, method) of records with fit_two_starts()'s
# record (diag$start): how many reps' default start failed or lost. theta0's fit
# takes over silently, so a dead default start shows nowhere else. Value: their
# share. NULL when no record has one.
q_start_rows <- function(id, recs) {
  recs <- Filter(function(r) !is.null(r$diag$start), recs)
  key <- vapply(recs, function(r) paste0(r$cell, "/", r$method), "")
  do.call(rbind, lapply(unique(key), function(s) {
    st <- lapply(recs[key == s], function(r) r$diag$start)
    failed <- vapply(st, function(x) !is.null(x$default_error), logical(1))
    lost <- vapply(st, function(x) identical(x$winner, "theta0"), logical(1)) & !failed
    q_row(id, sprintf("%s two starts: of %d reps the default start failed in %d, lost in %d", s,
                      length(st), sum(failed), sum(lost)), mean(failed | lost))
  }))
}

# The gated rows of one entry as check rows; "na" (gated, not computable) fails.
q_gate_rows <- function(id, g) {
  q_row(id, paste0(g$cell, "/", g$method, " ", g$metric,
                   ifelse(is.na(g$route), "", paste0(":", g$route)),
                   ifelse(is.na(g$estimand), "", paste0(" ", g$estimand))),
        g$value, g$tol, unname(c(pass = TRUE, fail = FALSE, na = FALSE)[g$verdict]))
}

q3_run <- function(ctx, entries) {
  ents <- q3_entries(ctx$fixtures)
  if (!identical(entries, "*")) ents <- ents[entries]
  m <- q_mc_run(ctx, "Q3", ents)
  rows <- m$rows

  # Each entry gated over the entry (one family): by the family's gates_quick
  # where it has one, else at the entry's cells with the family's gate metrics.
  for (k in names(m$plan)) {
    p <- m$plan[[k]]
    id <- paste0("Q3:", k)
    gates <- p$fam$gates_quick %||%
      (if (!is.null(p$fam$gates)) list(cell = p$cells, metrics = p$fam$gates$metrics))
    if (is.null(gates)) {
      rows[[length(rows) + 1L]] <- q_row(id, "no gates: a findings-only family", pass = FALSE)
      next
    }
    recs <- Filter(function(r) r$family == p$fam$id && r$cell %in% p$cells, m$recs)
    key <- vapply(recs, function(r) paste(r$cell, r$method, sep = "\r"), "")
    tab <- do.call(rbind, lapply(unique(key), function(s) {
      metrics_cell(recs[key == s], gates = gates, tol_set = "quick")
    }))
    if (is.null(tab)) {
      rows[[length(rows) + 1L]] <- q_row(id, "no records came back", pass = FALSE)
      next
    }
    tab <- gate_quick(tab, p$fam)
    utils::write.csv(tab, file.path(ctx$out, sprintf("quick_%s_Q3_%s.csv", ctx$stamp, k)),
                     row.names = FALSE)
    g <- tab[tab$gated & tab$metric != "failure_rate", , drop = FALSE]
    gc <- gates$cell
    cat(sprintf("\nQ3:%s gate table (%scells %s; tolerances quick):\n", k,
                if (is.null(p$fam$gates_quick)) "" else "gates_quick, ",
                if (is.list(gc)) paste(names(gc), vapply(gc, paste, "", collapse = ","),
                                       sep = ": ", collapse = "; ")
                else paste(gc, collapse = ", ")))
    rows[[length(rows) + 1L]] <- q_failed_rows(id, recs, tab)
    rows[[length(rows) + 1L]] <- q_start_rows(id, recs)
    if (!nrow(g)) {
      cat("  no gated rows\n")
      rows[[length(rows) + 1L]] <- q_row(id, "no gated rows: the entry checks nothing",
                                         pass = FALSE)
      next
    }
    pm <- ifelse(is.na(g$value), "NA", sprintf("%.4g +- %.2g", g$value, g$mcse))
    print(data.frame(cell = g$cell, method = g$method, metric = g$metric,
                     estimand = ifelse(is.na(g$estimand), "", g$estimand),
                     route = ifelse(is.na(g$route), "", g$route), value = pm, target = g$target,
                     tol = g$tol, p_holm = signif(g$p_holm, 2), verdict = g$verdict),
          row.names = FALSE)
    rows[[length(rows) + 1L]] <- q_gate_rows(id, g)
  }
  do.call(rbind, unname(rows))
}

q4_run <- function(ctx, entries) {
  m <- q_mc_run(ctx, "Q4", Q4_ENTRY)
  p <- m$plan$F9
  if (is.null(p) || !length(m$recs)) return(do.call(rbind, unname(m$rows)))
  chk <- q4_check(p$fam, m$recs)
  utils::write.csv(chk$tab, file.path(ctx$out, sprintf("quick_%s_Q4.csv", ctx$stamp)),
                   row.names = FALSE)
  r <- chk$rows[order(chk$rows$value), , drop = FALSE]
  cat(sprintf("\nQ4: %d quantities; smallest Holm p (threshold %g):\n", nrow(r), r$threshold[1]))
  print(utils::head(r[c("check", "value", "pass")], 5L), row.names = FALSE)
  rbind(do.call(rbind, unname(m$rows)), q_failed_rows("Q4", m$recs, chk$tab), chk$rows)
}

# --- Q2 -----------------------------------------------------------------------------

# The units selected, the costliest first, on the workers. F6's port is
# compiled here first, so that workers only load it.
q2_run <- function(ctx, entries) {
  keys <- if (identical(entries, "*")) names(Q2_UNITS) else entries
  if (any(vapply(Q2_UNITS[keys], function(u) u$family == "F6", logical(1)))) geo_gen_compile()
  keys <- keys[order(-vapply(Q2_UNITS[keys], `[[`, 1, "cost"))]
  cat(sprintf("Q2: %s on %d workers\n", paste(keys, collapse = ", "),
              min(ctx$workers, length(keys))))
  flush(stdout())
  do.call(rbind, parallel::parLapplyLB(ctx$cluster(length(keys)), as.list(keys), q2_task,
                                       chunk.size = 1L))
}

# --- The registry ------------------------------------------------------------------

QUICK_TESTS <- list(
  Q0     = function(ctx, entries) q0_checks(ctx$perturb),
  `Q1.1` = function(ctx, entries) ref_q11(),
  `Q1.2` = function(ctx, entries) ref_q12(),
  `Q1.3` = function(ctx, entries) ref_q13(),
  `Q1.4` = function(ctx, entries) ref_q14(),
  Q2     = q2_run,
  Q3     = q3_run,
  Q4     = q4_run
)

# --only as a list named by test: "*" for every entry, else entry names
# (Q3:F2,Q3:F3), checked against Q2's and Q3's entries; other tests have none.
q_select <- function(only, fixtures) {
  if (is.null(only)) {
    return(stats::setNames(rep(list("*"), length(QUICK_TESTS)), names(QUICK_TESTS)))
  }
  tok <- trimws(strsplit(only, ",", fixed = TRUE)[[1]])
  tok <- tok[nzchar(tok)]
  id <- sub(":.*$", "", tok)
  ent <- ifelse(grepl(":", tok, fixed = TRUE), sub("^[^:]*:", "", tok), "*")
  bad <- setdiff(id, names(QUICK_TESTS))
  if (length(bad)) {
    stop("unknown test(s) ", paste(bad, collapse = ", "), "; tests: ",
         paste(names(QUICK_TESTS), collapse = ", "), call. = FALSE)
  }
  sel <- list()
  for (i in seq_along(id)) {
    sel[[id[i]]] <- if ("*" %in% c(sel[[id[i]]], ent[i])) "*" else c(sel[[id[i]]], ent[i])
  }
  bad <- setdiff(sel$Q3, c("*", names(q3_entries(fixtures))))
  if (length(bad)) {
    stop("Q3 has no entry ", paste(bad, collapse = ", "), "; entries: ",
         paste(names(q3_entries(fixtures)), collapse = ", "),
         if (any(bad %in% names(Q3_FIXTURES))) " (a fixture: add --fixtures)", call. = FALSE)
  }
  bad <- setdiff(sel$Q2, c("*", names(Q2_UNITS)))
  if (length(bad)) {
    stop("Q2 has no entry ", paste(bad, collapse = ", "), "; entries: ",
         paste(names(Q2_UNITS), collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(names(sel)[!vapply(sel, identical, logical(1), "*")], c("Q2", "Q3"))
  if (length(bad)) {
    stop("only Q2 and Q3 have entries; ", paste(bad, collapse = ", "), " has none", call. = FALSE)
  }
  sel[intersect(names(QUICK_TESTS), names(sel))]
}

# --- Main -------------------------------------------------------------------------

t_start <- proc.time()[["elapsed"]]
options(width = 200L)
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE)[1])))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
args <- vs_args(defaults = list(only = NULL, workers = 8L, out = NULL, ref = NULL,
                                reps = 0L, perturb = 0, fixtures = FALSE))
if (args$workers < 1L) stop("--workers must be at least 1", call. = FALSE)
if (args$reps < 0L) stop("--reps must be a positive number of reps", call. = FALSE)
sel <- q_select(args$only, args$fixtures)
run_kind <- q_run_kind(args$only, args$reps, args$perturb, args$fixtures)
vs_threads(1L)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out <- if (is.null(args$out)) vs_out_dir(run = "quick") else vs_out_dir(out = args$out)
workers <- new.env()
ctx <- list(
  L = b$L, out = out, stamp = stamp, ref = args$ref, fixtures = args$fixtures,
  workers = args$workers, reps = args$reps, perturb = args$perturb,
  reps_dir = file.path(out, sprintf("quick_%s_reps", stamp)),
  # The workers: as many as asked and the test has tasks; a later test with
  # more tasks gets a larger cluster.
  cluster = function(n_tasks) {
    n <- min(args$workers, n_tasks)
    if (!is.null(workers$cl) && length(workers$cl) < n) {
      parallel::stopCluster(workers$cl)
      workers$cl <- NULL
    }
    if (is.null(workers$cl)) {
      Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
      workers$cl <- parallel::makePSOCKcluster(n)
      parallel::clusterCall(workers$cl, q_worker_init, root, Sys.getenv("R_LIBS"),
                            getNamespaceInfo("choicer", "path"))
    }
    workers$cl
  })
cat(sprintf("quick.R: choicer %s from %s; tests %s\n", b$build$version, b$build$path,
            paste(names(sel), collapse = ", ")))

rows <- list()
info <- list()
for (id in names(sel)) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(QUICK_TESTS[[id]](ctx, sel[[id]]), error = function(e) {
    q_row(id, paste("error:", conditionMessage(e)), pass = FALSE)
  })
  sec <- proc.time()[["elapsed"]] - t0
  if (!is.null(r)) r$pass[is.na(r$pass) & !is.na(r$threshold)] <- FALSE   # met a NaN
  if (is.null(r) || all(is.na(r$pass))) {
    r <- rbind(r, q_row(id, "the test returned no check row", pass = FALSE))
  }
  r$seconds <- sec
  info[[id]] <- r[is.na(r$pass), , drop = FALSE]
  rows[[id]] <- r[!is.na(r$pass), , drop = FALSE]
  bad <- r[r$pass %in% FALSE, , drop = FALSE]
  cat(sprintf("%s: %d checks, %d failed, %d information rows (%.1f s)\n", id, nrow(rows[[id]]),
              nrow(bad), nrow(info[[id]]), sec))
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  FAILED %s: %s = %s (threshold %s)\n", bad$test[i], bad$check[i],
                format(bad$value[i], digits = 4), format(bad$threshold[i], digits = 4)))
  }
  flush(stdout())
}
if (!is.null(workers$cl)) parallel::stopCluster(workers$cl)
res <- do.call(rbind, unname(rows))
inf <- do.call(rbind, unname(info))
res$run <- run_kind
csv <- file.path(out, sprintf("quick_%s.csv", stamp))
utils::write.csv(res, csv, row.names = FALSE)
if (nrow(inf)) {
  inf$run <- run_kind
  utils::write.csv(inf, file.path(out, sprintf("quick_%s_info.csv", stamp)), row.names = FALSE)
  cat("\nInformation rows:\n")
  print(inf[c("test", "check", "value")], row.names = FALSE, right = FALSE)
}

# Summary: per test, the checks passed and failed; then the failures.
summ <- do.call(rbind, lapply(split(res, factor(res$test, unique(res$test))), function(d) {
  data.frame(test = d$test[1], checks = nrow(d), pass = sum(d$pass %in% TRUE),
             fail = sum(d$pass %in% FALSE), seconds = round(d$seconds[1], 1))
}))
cat("\nSummary:\n")
print(summ, row.names = FALSE)
bad <- res[res$pass %in% FALSE, , drop = FALSE]
if (nrow(bad)) {
  cat("\nFailed checks:\n")
  print(bad[c("test", "check", "value", "threshold")], row.names = FALSE, right = FALSE)
}
if (dir.exists(ctx$reps_dir)) {
  cat(sprintf("\nThe mini Monte Carlos' report: Rscript %s --run=%s --tier=quick\n",
              file.path(root, "report.R"), ctx$reps_dir))
}
cat(sprintf("\n%s\ntotal %.1f s; %s; run: %s\n", csv, proc.time()[["elapsed"]] - t_start,
            if (nrow(bad)) sprintf("%d check(s) FAILED", nrow(bad)) else "no check failed",
            if (run_kind == "tier") "the whole quick tier"
            else paste0(run_kind, " (not the quick tier)")))
quit(save = "no", status = if (nrow(bad)) 1L else 0L)
