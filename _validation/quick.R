# ============================================================================
# The quick tier (plan §5; 15 minutes at 8 workers): a registry of tests
# Q0 ... Q4, each a function(ctx, entries) returning check rows
# data.frame(test, check, value, threshold, pass); pass NA marks a row that
# checked nothing (a test not implemented yet, a family not written yet).
# WP9 fills the tests in. Q3's mechanism is generic: a mini Monte Carlo of
# its entries (family, cells, R) on PSOCK workers through run_rep() (the
# long tier's paired reps and seeds), then metrics_cell(...,
# tol_set = "quick") and gate_family(), gated by the family's gates_quick
# where it has one, else at the entry's cells with the family's gate metrics.
#
# Usage (the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/quick.R
#       [--only=Q0,Q1.3,Q3:F2] [--workers=8] [--out=DIR] [--ref=NAME] [--fixtures]
#   --only      the tests to run (default all); Q3:<entry> runs one Q3 entry
#   --fixtures  adds Q3's T0 plumbing entry (never run by default)
# Writes into --out (default _validation/output/quick/): quick_<stamp>.csv
# (every check row, with the test's seconds), quick_<stamp>_Q3_<entry>.csv
# (the entry's metrics), and Q3's rep files under quick_<stamp>_reps/ with
# q3_plan.rds, the units Q3 ran. Render them, gated as here, with
#   Rscript _validation/report.R --run=<out>/quick_<stamp>_reps --tier=quick
# Prints a summary and the Q3 gate tables; exits 1 if any check fails.
# ============================================================================

# --- Q3: the mini Monte Carlo -------------------------------------------------------

# Q3's entries (plan §5): a family and, optionally, `cells`, `methods` and
# `R` (default: the cells of the family's R$quick with their reps, every
# method of the cell). They are gated by the family's gates_quick, else at
# the entry's cells with the family's gate metrics, with the quick
# tolerances.
Q3_ENTRIES <- list(
  F1 = list(family = "F1"),
  F2 = list(family = "F2"),
  F3 = list(family = "F3"),
  F4 = list(family = "F4"),
  F8 = list(family = "F8")   # WP9: methods = the MNL methods (plan: F8-MNL)
)
Q3_FIXTURES <- list(T0 = list(family = "T0"))

# A worker's task: one rep of one (family, cell), every method of the
# entry (families loaded once per worker). Errors come back as a string.
q3_task <- function(task) {
  tryCatch({
    fam <- get0(task$family, envir = .q_fams, inherits = FALSE)
    if (is.null(fam)) {
      fam <- load_family(task$family, .q_L)
      assign(task$family, fam, envir = .q_fams)
    }
    run_rep(fam, task$cell, task$rep, task$methods, ref = task$ref)
  }, error = function(e) conditionMessage(e))
}

q3_entries <- function(fixtures) c(Q3_ENTRIES, if (fixtures) Q3_FIXTURES)

q3_run <- function(ctx, entries) {
  ents <- q3_entries(ctx$fixtures)
  if (!identical(entries, "*")) ents <- ents[entries]
  rows <- list()
  plan <- list()
  tasks <- list()
  units <- list()   # q3_plan.rds: the units run, for report.R --tier=quick
  for (k in names(ents)) {
    e <- ents[[k]]
    fam <- tryCatch(load_family(e$family, ctx$L), error = function(err) err)
    if (inherits(fam, "error")) {
      rows[[k]] <- q_row(paste0("Q3:", k), paste("not run:", conditionMessage(fam)))
      next
    }
    rq <- fam$R$quick
    cells <- e$cells %||% names(rq)[rq > 0]
    if (!length(cells)) {
      rows[[k]] <- q_row(paste0("Q3:", k), "not run: the family has no quick cells (R$quick)")
      next
    }
    for (cl in cells) {
      ms <- intersect(cell_methods(fam, cl), e$methods %||% names(fam$methods))
      R <- e$R %||% rq[[cl]]
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
  if (!length(tasks)) return(do.call(rbind, rows))

  # The reps, largest N first and one per dispatch (parLapplyLB's default
  # sends two chunks per worker), then their rep files.
  tasks <- tasks[order(-vapply(tasks, `[[`, 1, "N"))]
  cat(sprintf("Q3: %d reps of %s on %d workers\n", length(tasks),
              paste(names(plan), collapse = ", "), min(ctx$workers, length(tasks))))
  res <- parallel::parLapplyLB(ctx$cluster(length(tasks)), tasks, q3_task, chunk.size = 1L)
  dir.create(ctx$reps_dir, showWarnings = FALSE, recursive = TRUE)
  write_metadata(ctx$reps_dir, ctx$ref)
  saveRDS(do.call(rbind, units), file.path(ctx$reps_dir, "q3_plan.rds"))
  recs <- list()
  for (i in seq_along(tasks)) {
    if (is.character(res[[i]])) {
      rows[[length(rows) + 1L]] <- q_row(paste0("Q3:", tasks[[i]]$entry),
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

  # Each entry gated over the entry (one family): by the family's
  # gates_quick where it has one, else at the entry's cells with the
  # family's gate metrics; a findings-only family (gates NULL) is not gated.
  key <- vapply(recs, function(r) paste(r$family, r$cell, r$method, sep = "\r"), "")
  for (k in names(plan)) {
    p <- plan[[k]]
    gates <- p$fam$gates_quick %||%
      (if (!is.null(p$fam$gates)) list(cell = p$cells, metrics = p$fam$gates$metrics))
    if (is.null(gates)) {
      rows[[length(rows) + 1L]] <- q_row(paste0("Q3:", k), "no gates: a findings-only family")
      next
    }
    sel <- unique(key[vapply(recs, function(r) r$family == p$fam$id && r$cell %in% p$cells,
                             logical(1))])
    tab <- do.call(rbind, lapply(sel, function(s) {
      metrics_cell(recs[key == s], gates = gates, tol_set = "quick")
    }))
    if (is.null(tab)) next
    tab <- gate_family(tab)
    utils::write.csv(tab, file.path(ctx$out, sprintf("quick_%s_Q3_%s.csv", ctx$stamp, k)),
                     row.names = FALSE)
    g <- tab[tab$gated, , drop = FALSE]
    gc <- gates$cell
    cat(sprintf("\nQ3:%s gate table (%scells %s; tolerances quick):\n", k,
                if (is.null(p$fam$gates_quick)) "" else "gates_quick, ",
                if (is.list(gc)) paste(names(gc), vapply(gc, paste, "", collapse = ","),
                                       sep = ": ", collapse = "; ")
                else paste(gc, collapse = ", ")))
    if (!nrow(g)) {
      cat("  no gated rows\n")
      next
    }
    pm <- ifelse(is.na(g$value), "NA", sprintf("%.4g +- %.2g", g$value, g$mcse))
    print(data.frame(cell = g$cell, method = g$method, metric = g$metric,
                     estimand = ifelse(is.na(g$estimand), "", g$estimand),
                     route = ifelse(is.na(g$route), "", g$route), value = pm, target = g$target,
                     tol = g$tol, p_holm = signif(g$p_holm, 2), verdict = g$verdict),
          row.names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      test = paste0("Q3:", k),
      check = paste0(g$cell, "/", g$method, " ", g$metric, ifelse(is.na(g$route), "",
                                                                   paste0(":", g$route)),
                     ifelse(is.na(g$estimand), "", paste0(" ", g$estimand))),
      value = g$value, threshold = g$tol,
      pass = c(pass = TRUE, fail = FALSE, na = FALSE)[g$verdict], stringsAsFactors = FALSE)
  }
  do.call(rbind, unname(rows))
}

# --- The registry ------------------------------------------------------------------

# A check row; pass NA checks nothing.
q_row <- function(test, check, value = NA_real_, threshold = NA_real_, pass = NA) {
  data.frame(test = test, check = check, value = value, threshold = threshold, pass = pass,
             stringsAsFactors = FALSE)
}

# A test not implemented yet (WP9): one row, pass NA.
q_stub <- function(test, what) {
  force(test)
  force(what)
  function(ctx, entries) q_row(test, paste("not implemented (WP9):", what))
}

QUICK_TESTS <- list(
  Q0     = q_stub("Q0", "the suite's probabilities at theta0 against predict() (MNL, NL, MXL)"),
  `Q1.1` = q_stub("Q1.1", "MNL on mode_choice against survival::clogit and mlogit"),
  `Q1.2` = q_stub("Q1.2", "NL (GROUND nest) against mlogit, un.nest.el = TRUE"),
  `Q1.3` = q_stub("Q1.3", "MXL, S = 4000, against Gauss-Hermite ML"),
  `Q1.4` = q_stub("Q1.4", "MNP against bayesm::rmnpGibbs"),
  Q2     = q_stub("Q2", "one large-N fit per family F1-F4, F6"),
  Q3     = q3_run,
  Q4     = q_stub("Q4", "mini-SBC of MNP, HMNL, HMNP")
)

# --only as a list named by test: "*" for every entry, else entry names
# (Q3:F2,Q3:F3), checked against Q3's entries.
q_select <- function(only, fixtures) {
  if (is.null(only)) return(stats::setNames(rep(list("*"), length(QUICK_TESTS)), names(QUICK_TESTS)))
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
  sel[intersect(names(QUICK_TESTS), names(sel))]
}

# --- Workers ----------------------------------------------------------------------

# Runs on each worker: the master's R_LIBS, the build guard and the library
# (vs_boot()), the same choicer as the master's, one thread. Returns the
# worker's pid.
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
  Sys.getpid()
}

# --- Main -------------------------------------------------------------------------

t_start <- proc.time()[["elapsed"]]
options(width = 200L)
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE)[1])))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
args <- vs_args(defaults = list(only = NULL, workers = 8L, out = NULL, ref = NULL,
                                fixtures = FALSE))
if (args$workers < 1L) stop("--workers must be at least 1", call. = FALSE)
sel <- q_select(args$only, args$fixtures)
vs_threads(1L)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out <- if (is.null(args$out)) vs_out_dir(run = "quick") else vs_out_dir(out = args$out)
workers <- new.env()
ctx <- list(
  root = root, L = b$L, out = out, stamp = stamp, ref = args$ref, fixtures = args$fixtures,
  workers = args$workers, reps_dir = file.path(out, sprintf("quick_%s_reps", stamp)),
  cluster = function(n_tasks) {
    if (is.null(workers$cl)) {
      n <- min(args$workers, n_tasks)
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
for (id in names(sel)) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(QUICK_TESTS[[id]](ctx, sel[[id]]), error = function(e) {
    q_row(id, paste("error:", conditionMessage(e)), pass = FALSE)
  })
  if (!is.null(r) && nrow(r)) {
    r$seconds <- proc.time()[["elapsed"]] - t0
    rows[[id]] <- r
  }
}
if (!is.null(workers$cl)) parallel::stopCluster(workers$cl)
res <- do.call(rbind, unname(rows))
csv <- file.path(out, sprintf("quick_%s.csv", stamp))
utils::write.csv(res, csv, row.names = FALSE)

# Summary: per test, the checks passed, failed and not run; then the failures.
summ <- do.call(rbind, lapply(split(res, factor(res$test, unique(res$test))), function(d) {
  data.frame(test = d$test[1], checks = nrow(d), pass = sum(d$pass %in% TRUE),
             fail = sum(d$pass %in% FALSE), not_run = sum(is.na(d$pass)),
             seconds = round(d$seconds[1], 1))
}))
cat("\nSummary:\n")
print(summ, row.names = FALSE)
bad <- res[res$pass %in% FALSE, , drop = FALSE]
if (nrow(bad)) {
  cat("\nFailed checks:\n")
  print(bad[c("test", "check", "value", "threshold")], row.names = FALSE)
}
if (dir.exists(ctx$reps_dir)) {
  cat(sprintf("\nQ3's report: Rscript %s --run=%s --tier=quick\n", file.path(root, "report.R"),
              ctx$reps_dir))
}
cat(sprintf("\n%s\ntotal %.1f s; %s\n", csv, proc.time()[["elapsed"]] - t_start,
            if (nrow(bad)) sprintf("%d check(s) FAILED", nrow(bad)) else "no check failed"))
quit(save = "no", status = if (nrow(bad)) 1L else 0L)
