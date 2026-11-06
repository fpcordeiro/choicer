# ============================================================================
# Self-check of WP9: quick.R on small settings, through its command line.
#   entries  Q3's entries and gates_quick, Q2's units, against the families: the
#            bottom rungs, R = 100, gate entries that parse, no "blp" in Q2; Q2's
#            SE tolerance; the run marker; the row builders of Q1.3 and Q1.4.
#   q0       --only=Q0 passes with no information row; --perturb fails it (exit
#            1) on every comparison, its names rows standing.
#   q1       --only=Q1.1 and Q1.2 pass.
#   q2       --only=Q2:F1,Q2:F6.mnl runs those units alone; unknown selections
#            stop before anything runs.
#   q3, q4   --only=Q3:F1,Q3:F3,Q3:F8,Q4 --reps=3: the entries' rep files (F8:
#            mnl_unweighted at quota alone), gates_quick's table, every row a
#            verdict, the two-start information row, every quantity of Q4's three
#            cells; q4_check()'s floor on a quantity's reps.
#   report   report.R --tier=quick on those rep files gates as quick.R did.
#   verdicts report.R --tier=quick against quick.R's q3_run() and q4_run() on the
#            same fixture records, no fits: F9 with 10 of 100 failed MNP records,
#            F9 with a Holm p between 0.001 and 0.01, F1 with 5 of 100 failed, and
#            controls: every gated verdict, by row, on both sides. A findings-only
#            family gets no gated row.
#   inject   failures, each from a tree copy, each exit 1 with the row printed: a
#            comparison that met a NaN; a NaN BHHH SE of a constant in Q2; a gated
#            row "na"; a missing reference; a test without a check row; a family
#            that does not load; failed reps beyond the 2% tolerance.
# No check row of any run may have pass NA, and no row with a threshold may be
# filed as information (no_na()).
#
# Usage (about a minute; one OpenMP thread, the dev build first on R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp9_quick.R [--only=entries,q0,...] [--out=DIR]
#       [--workers=2]
# Exits with status 1 on any FAIL.
# ============================================================================

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                               value = TRUE)[1]))))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
args <- vs_args(defaults = list(only = "entries,q0,q1,q2,q3,q4,report,verdicts,inject", out = NULL,
                                workers = 2L))
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
out <- if (is.null(args$out)) tempfile("wp9_") else args$out
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out)
options(width = 160L)
t_start <- proc.time()[["elapsed"]]
cat("choicer", b$build$version, "from", b$build$path, "\n")

rscript <- file.path(R.home("bin"), "Rscript")
env <- c(paste0("R_LIBS=", Sys.getenv("R_LIBS")), "OMP_NUM_THREADS=1", "OMP_THREAD_LIMIT=1")
results <- list()
notes <- character()
note <- function(...) notes <<- c(notes, paste0("      ", strsplit(paste0(...), "\n")[[1]]))
check <- function(name, expr) {
  notes <<- character()
  res <- tryCatch(expr, error = function(e) structure(FALSE, msg = conditionMessage(e)))
  cat(sprintf("%s  %s\n", if (isTRUE(res)) "PASS" else "FAIL", name))
  if (!is.null(attr(res, "msg"))) note("error: ", attr(res, "msg"))
  if (length(notes)) cat(notes, sep = "\n")
  results[name] <<- isTRUE(res)
}

# quick.R in tree `at` (the suite itself by default), each run with its own
# output directory: its exit status, printed lines, seconds, check rows and
# information rows.
runs <- new.env()
quick <- function(tag, ..., at = root) {
  d <- file.path(out, tag)
  t0 <- proc.time()[["elapsed"]]
  o <- suppressWarnings(system2(rscript, c(shQuote(file.path(at, "quick.R")),
                                           paste0("--out=", shQuote(d)), ...),
                                stdout = TRUE, stderr = TRUE, env = env))
  csv <- list.files(d, pattern = "^quick_[0-9-]+\\.csv$", full.names = TRUE)
  inf <- list.files(d, pattern = "^quick_[0-9-]+_info\\.csv$", full.names = TRUE)
  r <- list(out = o, status = attr(o, "status") %||% 0L, sec = proc.time()[["elapsed"]] - t0,
            dir = d, rows = if (length(csv)) utils::read.csv(csv[1], stringsAsFactors = FALSE),
            info = if (length(inf)) utils::read.csv(inf[1], stringsAsFactors = FALSE))
  note(sprintf("%s: exit %d, %.1f s, %d rows, %d information rows; run: %s", tag, r$status, r$sec,
               NROW(r$rows), NROW(r$info), paste(unique(r$rows$run), collapse = " | ")))
  runs[[tag]] <- r
  r
}

# Every row of a run has its verdict: the check rows TRUE or FALSE, and what was
# filed as information has neither a verdict nor a threshold.
no_na <- function(r) {
  !is.null(r$rows) && nrow(r$rows) > 0L && !anyNA(r$rows$pass) && !anyNA(r$rows$run) &&
    all(is.na(r$info$pass)) && all(is.na(r$info$threshold))
}
printed <- function(r, s) any(grepl(s, r$out, fixed = TRUE))

# A copy of the suite's scripts and library, patched by appending to lib files.
tree <- function(tag, patches) {
  d <- file.path(out, paste0("tree_", tag), "_validation")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(root, c("quick.R", "report.R", "lib")), d, recursive = TRUE)
  for (f in names(patches)) {
    cat("\n", patches[[f]], "\n", file = file.path(d, "lib", f), append = TRUE, sep = "")
  }
  d
}

# Definitions of quick.R, parsed from the script (its top level runs a tier).
quick_defs <- function(keep) {
  defs <- new.env(parent = globalenv())
  for (e in parse(file.path(root, "quick.R"), keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% keep) {
      eval(e, defs)
    }
  }
  defs
}

# --- entries -------------------------------------------------------------------------

if ("entries" %in% only) {
  defs <- quick_defs(c("Q3_ENTRIES", "Q3_FIXTURES", "Q4_ENTRY"))
  check("Q3's entries are the plan's bottom rungs at R = 100, with parsing gates_quick", {
    want <- list(F1 = "500", F2 = "1k", F3 = "1k.q", F4 = "200.q", F8 = "quota")
    ok <- setequal(names(defs$Q3_ENTRIES), names(want)) &&
      identical(defs$Q3_ENTRIES$F8$methods, "mnl_unweighted")
    for (k in names(want)) {
      fam <- load_family(k, b$L)
      ms <- intersect(cell_methods(fam, want[[k]]),
                      defs$Q3_ENTRIES[[k]]$methods %||% names(fam$methods))
      rq <- fam$R$quick
      gq <- fam$gates_quick
      for (m in ms) .parse_gates(gq, m)
      gated <- if (is.list(gq$cell)) unique(unlist(gq$cell)) else gq$cell
      note(sprintf("%s: cells %s, R %s, methods %s, gates_quick at %s: %s", k,
                   paste(names(rq), collapse = ","), paste(rq, collapse = ","),
                   paste(ms, collapse = ","), paste(gated, collapse = ","),
                   paste(unlist(gq$metrics), collapse = " ")))
      ok <- ok && identical(names(rq), want[[k]]) && all(rq == 100) && !is.null(gq) &&
        identical(gated, want[[k]])
    }
    ok && "bias_std@asc" %in% load_family("F8", b$L)$gates_quick$metrics$mnl_unweighted
  })
  check("Q2's units: families and cells exist, S where the family draws, no \"blp\" fit", {
    ok <- all(c("F1", "F2", "F3", "F4", "F6") %in% vapply(Q2_UNITS, `[[`, "", "family"))
    for (k in names(Q2_UNITS)) {
      u <- Q2_UNITS[[k]]
      fam <- load_family(u$family, b$L)
      cell <- family_cell(fam, u$base)
      ok <- ok && !is.null(u$cost) &&
        !any(grepl("blp", deparse(body(u$fit)), fixed = TRUE)) && "blp" != u$base
      note(sprintf("%s: %s at %s, N = %d%s", k, u$family, u$base, u$N %||% cell$N,
                   if (!is.null(cell$S) && !is.na(cell$S)) sprintf(", S = %d", cell$S) else ""))
    }
    ok
  })
  check("Q2's SE tolerance: 10%, for a constant chosen n times at least 2 / sqrt(n)", {
    fit <- list(coefficients = c(b = 1, ASC_x = 0, ASC_y = 0, ASC_z = 0),
                param_map = list(beta = 1L, asc = 2:4))
    data <- data.frame(alt = c(rep("x", 25L), rep("y", 900L), "z"), choice = c(rep(1L, 925L), 0L))
    tol <- .q2_se_tol(fit, data)
    note(paste(names(tol), signif(tol, 3), collapse = "; "))
    identical(names(tol), names(fit$coefficients)) && all(tol[c("b", "ASC_y")] == 0.10) &&
      tol[["ASC_x"]] == 0.4 && is.infinite(tol[["ASC_z"]])
  })
  check("the run marker: \"tier\" only for the whole tier unshortened", {
    k <- list(q_run_kind(), q_run_kind(only = "Q0"), q_run_kind(reps = 2L),
              q_run_kind(perturb = 1e-6), q_run_kind(fixtures = TRUE))
    note(paste(unlist(k), collapse = " | "))
    identical(k[[1]], "tier") && !any(unlist(k[-1]) == "tier")
  })
  check("the comparisons of Q1.3 and Q1.4 go through .ref_row(), so that a NaN fails", {
    files <- file.path(root, "lib", c("references_quad.R", "dgp_bayes.R"))
    src <- unlist(lapply(files, readLines))
    note(sprintf(".ref_row calls: %d", sum(grepl(".ref_row(", src, fixed = TRUE))))
    !any(grepl("data.frame(test", src, fixed = TRUE)) &&
      sum(grepl(".ref_row(", src, fixed = TRUE)) >= 5L &&
      .ref_row("x", "a NaN", NaN, 1)$pass %in% FALSE
  })
}

# --- q0 ------------------------------------------------------------------------------

if ("q0" %in% only) {
  check("Q0 passes: every comparison below its tolerance, no information row", {
    r <- quick("q0", "--only=Q0")
    v <- r$rows[grepl("max \\||, rel$", r$rows$check), ]
    note(sprintf("%d rows, largest: %s", nrow(r$rows),
                 paste(sprintf("%.1e", tapply(v$value, v$threshold, max)), collapse = ", ")))
    r$status == 0L && no_na(r) && all(r$rows$pass) && nrow(v) >= 21L &&
      all(v$threshold == 1e-12) && all(v$value < v$threshold) && is.null(r$info) &&
      any(grepl("^D2 mxl: the DGP's own tastes", r$rows$check)) &&
      all(r$rows$run == "only=Q0") && printed(r, "(not the quick tier)")
  })
  check("Q0 with --perturb=1e-6 (theta on choicer's side) fails every comparison, exit 1", {
    r <- quick("q0_perturbed", "--only=Q0", "--perturb=1e-6")
    v <- r$rows[grepl("max \\||, rel$", r$rows$check), ]
    n <- r$rows[grepl("names\\(coef", r$rows$check), ]
    note(sprintf("comparisons failed: %d of %d; smallest value %.1e", sum(!v$pass), nrow(v),
                 min(v$value)))
    r$status == 1L && no_na(r) && !any(v$pass) && all(n$pass) && printed(r, "FAILED Q0") &&
      printed(r, "check(s) FAILED") && all(grepl("perturb=", r$rows$run, fixed = TRUE))
  })
}

# --- q1 ------------------------------------------------------------------------------

if ("q1" %in% only) {
  check("Q1.1 and Q1.2 pass", {
    r <- quick("q1", "--only=Q1.1,Q1.2")
    note(paste(sprintf("%s %.1e < %.0e", r$rows$test, r$rows$value, r$rows$threshold),
               collapse = "; "))
    r$status == 0L && no_na(r) && all(r$rows$pass) &&
      setequal(r$rows$test, c("Q1.1", "Q1.2")) && all(r$rows$value < r$rows$threshold)
  })
}

# --- q2 ------------------------------------------------------------------------------

if ("q2" %in% only) {
  check("--only=Q2:F1,Q2:F6.mnl runs those units alone and passes", {
    r <- quick("q2", "--only=Q2:F1,Q2:F6.mnl", paste0("--workers=", args$workers))
    keys <- unique(sub(" .*$", "", r$rows$check))
    for (i in seq_len(nrow(r$rows))) note(sprintf("%s: %.3g", r$rows$check[i], r$rows$value[i]))
    r$status == 0L && no_na(r) && all(r$rows$pass) && setequal(keys, c("F1", "F6.mnl")) &&
      nrow(r$rows) == 8L && NROW(r$info) == 2L
  })
  check("unknown selections stop before anything runs, naming the entries", {
    cases <- list(c("--only=Q2:F9", "Q2 has no entry F9"),
                  c("--only=Q3:F6", "Q3 has no entry F6"),
                  c("--only=Q9", "unknown test"),
                  c("--only=Q1.1:x", "only Q2 and Q3 have entries"),
                  c("--only=Q3:T0", "add --fixtures"))
    ok <- TRUE
    for (k in seq_along(cases)) {
      r <- quick(paste0("sel", k), cases[[k]][1])
      ok <- ok && r$status != 0L && printed(r, cases[[k]][2]) && is.null(r$rows)
    }
    ok
  })
}

# --- q3, q4, report ------------------------------------------------------------------

if (any(c("q3", "q4", "report") %in% only)) {
  w <- paste0("--workers=", args$workers)
  check("--only=Q3:F1,Q3:F3,Q3:F8,Q4 --reps=3: the entries' rep files; every row a verdict", {
    r <- quick("mc", "--only=Q3:F1,Q3:F3,Q3:F8,Q4", "--reps=3", w)
    reps <- rep_index(list.files(r$dir, pattern = "_reps$", full.names = TRUE)[1])
    u <- unique(paste(reps$family, reps$cell, reps$method))
    note(paste(u, collapse = "; "))
    note(paste(names(table(r$rows$test)), table(r$rows$test), collapse = "; "))
    fr <- r$rows[grepl("failed reps", r$rows$check, fixed = TRUE), ]
    st <- r$info[grepl("two starts", r$info$check, fixed = TRUE), ]
    note(paste(st$check, collapse = "; "))
    r$status == as.integer(any(!r$rows$pass)) && no_na(r) && nrow(reps) == 3L * 6L &&
      setequal(u, c("F1 500 ml", "F3 1k.q ml", "F8 quota mnl_unweighted", "F9 q_sbc_mnp mnp",
                    "F9 q_sbc_hmnl hmnl", "F9 q_sbc_hmnp hmnp")) &&
      printed(r, "Q3:F1 gate table (gates_quick") && printed(r, "Q3:F8 gate table (gates_quick") &&
      all(r$rows$pass) && all(c("Q3:F1", "Q3:F3", "Q3:F8", "Q4") %in% r$rows$test) &&
      all(fr$threshold == tolerances("quick")[["failure"]]) &&
      any(grepl("failed reps: 0 of 3", fr$check, fixed = TRUE)) &&
      sum(r$rows$test == "Q4" & r$rows$threshold == 0.001) > 10L &&
      any(grepl("q_sbc_mnp .* \\(R = 3 of 3\\)", r$rows$check)) &&
      nrow(st) == 1L && all(st$test == "Q3:F3") &&
      all(r$rows$run == "only=Q3:F1,Q3:F3,Q3:F8,Q4 reps=3")
  })
  check("q4_check(): a quantity whose ranks cover under 95% of the reps fails", {
    r <- runs$mc
    recs <- read_reps(list.files(r$dir, pattern = "_reps$", full.names = TRUE)[1], "F9")
    fam <- load_family("F9", b$L)
    i <- which(vapply(recs, `[[`, "", "cell") == "q_sbc_mnp")[1L]
    q <- names(recs[[i]]$diag$sbc$ranks)[1L]
    recs[[i]]$diag$sbc$ranks <- recs[[i]]$diag$sbc$ranks[-1L]
    rows <- q4_check(fam, recs)$rows
    hit <- rows[startsWith(rows$check, paste0("q_sbc_mnp ", q, " (R = ")), ]
    note(sprintf("%s: Holm p %.3f, pass %s", hit$check, hit$value, hit$pass))
    nrow(hit) == 1L && !hit$pass && hit$value > hit$threshold &&
      grepl("(R = 2 of 3)", hit$check, fixed = TRUE) && sum(!rows$pass) == 1L
  })
  if ("report" %in% only) check("report.R --tier=quick gates the rep files as quick.R did", {
    r <- runs$mc
    rd <- list.files(r$dir, pattern = "_reps$", full.names = TRUE)[1]
    o <- suppressWarnings(system2(rscript, c(shQuote(file.path(root, "report.R")),
                                             paste0("--run=", shQuote(rd)), "--tier=quick"),
                                  stdout = TRUE, stderr = TRUE, env = env))
    note(sprintf("report exit %s", attr(o, "status") %||% 0L))
    same <- vapply(c("F1", "F3", "F8"), function(k) {
      q <- utils::read.csv(list.files(r$dir, pattern = sprintf("_Q3_%s\\.csv$", k),
                                      full.names = TRUE)[1])
      p <- utils::read.csv(file.path(rd, "tables", sprintf("%s_metrics.csv", k)))
      key <- function(t) {
        sort(paste(t$cell, t$method, t$metric, t$estimand, t$route, t$verdict)[t$gated])
      }
      identical(key(q), key(p))
    }, TRUE)
    note(paste(names(same), same, collapse = "; "))
    is.null(attr(o, "status")) && all(same)
  })
}

# --- verdicts ------------------------------------------------------------------------

# The report against the runner on fixture records, no fits. The runner is quick.R's own
# q3_run() or q4_run(), with q_mc_run() replaced by one that hands them the records; the
# report is report.R --tier=quick, in a subprocess, on those records' rep files. Each gated
# row of its metrics CSV is labelled as the runner labels its check row, and `agree()` asks
# for the same rows with the same verdicts: every row, not the counts.
if ("verdicts" %in% only) {
  verdicts <- function(tag, recs) {
    d <- file.path(out, paste0("verdicts_", tag))
    dir.create(d, showWarnings = FALSE)
    for (r in recs) write_rep(r, rep_path(d, r$family, r$cell, r$method, r$rep))
    fam <- load_family(recs[[1]]$family, b$L)
    cells <- vapply(recs, `[[`, "", "cell")
    plan <- unique(data.frame(family = fam$id, cell = cells,
                              method = vapply(recs, `[[`, "", "method"), R = 100L,
                              stringsAsFactors = FALSE))
    saveRDS(plan, file.path(d, "q3_plan.rds"))
    defs <- quick_defs(c("Q3_ENTRIES", "Q3_FIXTURES", "Q4_ENTRY", "q3_entries", "q3_run", "q4_run",
                         "q_failed_rows", "q_start_rows", "q_gate_rows"))
    defs$q_mc_run <- function(ctx, test, ents) {
      list(plan = stats::setNames(list(list(fam = fam, cells = unique(cells))), fam$id),
           recs = recs, rows = list())
    }
    ctx <- list(fixtures = FALSE, out = out, stamp = tag)
    utils::capture.output(rows <- if (fam$id == "F9") defs$q4_run(ctx, "*") else
      defs$q3_run(ctx, "*"))
    o <- suppressWarnings(system2(rscript, c(shQuote(file.path(root, "report.R")),
                                             paste0("--run=", shQuote(d)), "--tier=quick"),
                                  stdout = TRUE, stderr = TRUE, env = env))
    tab <- utils::read.csv(file.path(d, "tables", paste0(fam$id, "_metrics.csv")),
                           stringsAsFactors = FALSE)
    g <- tab[tab$gated, , drop = FALSE]
    lab <- defs$q_gate_rows("report", g)$check
    f <- g$metric == "failure_rate"
    n <- g$R[f][match(paste(g$cell, g$method), paste(g$cell, g$method)[f])]
    s <- g$metric == "sbc_ecdf"
    lab[s] <- sprintf("%s %s (R = %d of %d)", g$cell[s], g$estimand[s], g$R[s], n[s])
    lab[f] <- sprintf("%s/%s failed reps: %d of %d", g$cell[f], g$method[f], g$n_failed[f], g$R[f])
    run <- rows[!is.na(rows$pass), , drop = FALSE]
    list(runner = stats::setNames(run$pass, sub("( failed reps: [0-9]+ of [0-9]+) \\(.*\\)$", "\\1",
                                                run$check)),
         report = stats::setNames(g$verdict == "pass", lab), rows = run, g = g,
         out = o, status = attr(o, "status") %||% 0L)
  }
  # The same rows with the same verdicts, `n` of them, failing exactly `fails`.
  agree <- function(v, n, fails = character()) {
    note(sprintf("runner: %d rows, %d fail; report: %d rows, %d fail; report exit %d",
                 length(v$runner), sum(!v$runner), length(v$report), sum(!v$report), v$status))
    for (k in names(v$runner)[!v$runner]) note("runner fails: ", k)
    for (k in names(v$report)[!v$report]) note("report fails: ", k)
    v$status == 0L && length(v$runner) == n && !anyDuplicated(names(v$runner)) &&
      setequal(names(v$runner), names(v$report)) &&
      identical(v$runner[names(v$report)], v$report) && setequal(names(v$runner)[!v$runner], fails)
  }

  # F9's records: `ranks(cell)` gives the 100 ranks of the quantity q, NA a failed rep. The
  # ranks are in 0..99; 0:99 is the even spread.
  f9 <- load_family("F9", b$L)
  f9_cells <- names(f9$R$quick)[f9$R$quick > 0]
  f9_recs <- function(ranks) {
    unlist(lapply(f9_cells, function(cl) {
      rk <- ranks(cl)
      lapply(1:100, function(r) {
        if (is.na(rk[r])) {
          new_rep_record("F9", cl, cell_methods(f9, cl), r, seed = r, ok = FALSE,
                         error = "injected")
        } else {
          new_rep_record("F9", cl, cell_methods(f9, cl), r, seed = r,
                         diag = list(sbc = list(L = 99L, ranks = c(q = rk[r]))))
        }
      })
    }), recursive = FALSE)
  }
  mnp <- "q_sbc_mnp"
  check("verdicts, F9 with 10 of 100 MNP records failed: the runner and the report fail MNP", {
    v <- verdicts("f9_failed", f9_recs(function(cl) {
      if (cl == mnp) c(rep(NA, 10L), floor(0:89 * 100 / 90)) else 0:99
    }))
    cov <- v$rows[grepl("(R = 90 of 100)", v$rows$check, fixed = TRUE), ]
    note(sprintf("%s: Holm p %.3f, above %g", cov$check, cov$value, cov$threshold))
    agree(v, 6L, c("q_sbc_mnp q (R = 90 of 100)", "q_sbc_mnp/mnp failed reps: 10 of 100")) &&
      cov$value > cov$threshold && printed(v, "F9: 6 gated rows, 4 pass, 2 fail, 0 na")
  })
  check("verdicts, F9 with complete ranks and a Holm p in (0.001, 0.01): both pass at Q4's alpha", {
    skew <- floor(100 * ((1:100 - 0.5) / 100)^1.6)
    v <- verdicts("f9_alpha", f9_recs(function(cl) if (cl == mnp) skew else 0:99))
    p <- v$g$p_holm[v$g$metric == "sbc_ecdf" & v$g$cell == mnp]
    note(sprintf("the skewed cell's Holm p: %.5f", p))
    agree(v, 6L) && p > 0.001 && p < 0.01 && printed(v, "F9: 6 gated rows, 6 pass, 0 fail, 0 na")
  })

  # F1's records at cell 500: its coefficients with exact SEs, each rep's z and LR drawn on a
  # frozen seed, otherwise unremarkable. The first `failed` reps fail; the others are the same
  # draws whichever `failed` is.
  th <- D2_THETA0$default
  f1_recs <- function(failed = 0L) {
    beta <- !startsWith(names(th), "ASC_")
    vs_with_seed(20261106L, lapply(1:100, function(r) {
      e <- th + 0.1 * stats::rnorm(length(th))
      lr <- stats::rchisq(1L, 1)
      if (r <= failed) {
        return(new_rep_record("F1", "500", "ml", r, seed = r, ok = FALSE, error = "injected"))
      }
      est <- data.frame(name = names(th), block = ifelse(beta, "beta", "asc"), est = unname(e),
                        truth = unname(th), se_hessian = 0.1, se_bhhh = NA_real_,
                        se_robust = NA_real_, se_cluster = NA_real_, gated = TRUE,
                        stringsAsFactors = FALSE)
      new_rep_record("F1", "500", "ml", r, seed = r, est = est,
                     joint_wald = joint_wald(e[beta], th[beta], diag(0.01, sum(beta))),
                     lr = lr_test(lr / 2, 0, 1))
    }))
  }
  check("verdicts, F1 with 5 of 100 failed: both fail the unit, and agree on the other 31 rows", {
    v <- verdicts("f1_failed", f1_recs(5L))
    agree(v, 32L, "500/ml failed reps: 5 of 100") &&
      printed(v, "F1: 32 gated rows, 31 pass, 1 fail")
  })
  check("verdicts, a control where nothing is wrong: every verdict passes on both sides", {
    f <- verdicts("f1_clean", f1_recs())
    q <- verdicts("f9_clean", f9_recs(function(cl) 0:99))
    agree(f, 32L) && agree(q, 6L)
  })
  check("verdicts: a findings-only family gets no gated row in the quick tier", {
    nm <- c("price/ref", "price/S100/g1", "price/S200/g1")
    v <- lapply(1:4, function(d) {
      new_rep_record("F5", "F4", "msl", d, seed = d, ok = d > 1L, error = if (d == 1L) "injected",
                     est = if (d > 1L) data.frame(
                       name = nm, block = "beta", est = c(0, 0.2, 0.1) * d, truth = 0,
                       se_hessian = NA_real_, se_bhhh = NA_real_, se_robust = NA_real_,
                       se_cluster = NA_real_, gated = FALSE, stringsAsFactors = FALSE),
                     diag = if (d > 1L) list(se_ref = c(price = 1), S_ref = 400L, units = 800L,
                                             src = c(S_src = 100L)))
    })
    d <- file.path(out, "verdicts_f5")
    dir.create(d, showWarnings = FALSE)
    for (r in v) write_rep(r, rep_path(d, r$family, r$cell, r$method, r$rep))
    o <- suppressWarnings(system2(rscript, c(shQuote(file.path(root, "report.R")),
                                             paste0("--run=", shQuote(d)), "--tier=quick"),
                                  stdout = TRUE, stderr = TRUE, env = env))
    tab <- utils::read.csv(file.path(d, "tables", "F5_metrics.csv"), stringsAsFactors = FALSE)
    fr <- tab[tab$metric == "failure_rate", ]
    note(sprintf("%d rows, %d gated; failure_rate %.2f, verdict %s", nrow(tab), sum(tab$gated),
                 fr$value, fr$verdict))
    is.null(attr(o, "status")) && nrow(tab) > 0L && !any(tab$gated) && nrow(fr) == 1L &&
      fr$value == 0.25 && fr$verdict == "finding"
  })
}

# --- inject --------------------------------------------------------------------------

if ("inject" %in% only) {
  check("a comparison that met a NaN fails; it is never filed as information", {
    d <- tree("nan", list(references.R = paste(
      'ref_q11 <- function() {',
      '  rbind(q_row("Q1.1", "a comparison that met a NaN", NaN, 1e-4, pass = NaN < 1e-4),',
      '        q_row("Q1.1", "an information row", 1))',
      '}',
      'ref_q12 <- function() .ref_row("Q1.2", "a NaN through the constructor", NaN, 1e-4)',
      sep = "\n")))
    r <- quick("nan", "--only=Q1.1,Q1.2", at = d)
    note(paste(r$rows$check, r$rows$pass, collapse = "; "))
    note(paste(r$info$check, collapse = "; "))
    r$status == 1L && no_na(r) && nrow(r$rows) == 2L && !any(r$rows$pass) &&
      identical(r$info$check, "an information row") && printed(r, "FAILED Q1.1") &&
      printed(r, "FAILED Q1.2")
  })
  check("a NaN BHHH SE on the rarest constant fails Q2's check", {
    d <- tree("nan_se", list(fit.R = paste(
      'fit_vcovs <- local({',
      '  f <- fit_vcovs',
      '  function(fit, routes = "hessian", cluster = NULL) {',
      '    V <- f(fit, routes, cluster)',
      '    a <- as.integer(fit$param_map$asc)',
      '    i <- a[which.max(diag(V$bhhh)[a])]',
      '    V$bhhh[i, i] <- NaN',
      '    V',
      '  }',
      '})', sep = "\n")))
    r <- quick("nan_se", "--only=Q2:F6.mnl", "--workers=1", at = d)
    k <- r$rows[grepl("SE_bhhh", r$rows$check, fixed = TRUE), ]
    rb <- r$rows[grepl("SE_robust", r$rows$check, fixed = TRUE), ]
    note(paste(k$check, k$value, k$pass))
    r$status == 1L && no_na(r) && nrow(k) == 1L && !k$pass && is.na(k$value) && all(rb$pass) &&
      printed(r, "FAILED Q2")
  })
  check("a gated row \"na\" fails the entry: exit 1, the row printed", {
    d <- tree("na", list(`families/T0.R` = paste0(
      'family$gates_quick <- list(cell = "200", metrics = list(ml = "lr_size"))')))
    r <- quick("na", "--only=Q3:T0", "--fixtures", "--workers=1", at = d)
    k <- r$rows[grepl("lr_size", r$rows$check, fixed = TRUE), ]
    note(paste(k$check, k$value, k$pass))
    r$status == 1L && no_na(r) && nrow(k) >= 1L && !any(k$pass) && printed(r, "FAILED Q3:T0")
  })
  check("a missing reference is a failed row, printed; the other tests still run", {
    d <- tree("ref", list(
      references.R = 'ref_q11 <- function() { vs_require("no.such.package"); NULL }'))
    r <- quick("ref", "--only=Q1.1,Q1.2", at = d)
    k <- r$rows[r$rows$test == "Q1.1", ]
    note(k$check)
    r$status == 1L && no_na(r) && nrow(k) == 1L && !k$pass &&
      grepl("Missing package(s): no.such.package", k$check, fixed = TRUE) &&
      printed(r, "FAILED Q1.1") && all(r$rows$pass[r$rows$test == "Q1.2"])
  })
  check("a test that returns no check row, or information alone, fails", {
    d <- tree("empty", list(references.R = paste(
      'ref_q11 <- function() NULL',
      'ref_q12 <- function() q_row("Q1.2", "information alone", 1)', sep = "\n")))
    r <- quick("empty", "--only=Q1.1,Q1.2", at = d)
    note(paste(r$rows$check, collapse = "; "))
    r$status == 1L && no_na(r) && sum(grepl("no check row", r$rows$check, fixed = TRUE)) == 2L
  })
  check("a family that does not load fails its entry, the others still run", {
    d <- tree("nofam", list())
    unlink(file.path(d, "lib", "families", "F8.R"))
    r <- quick("nofam", "--only=Q3:F8,Q3:F1", "--reps=2", "--workers=1", at = d)
    k <- r$rows[r$rows$test == "Q3:F8", ]
    note(k$check)
    r$status == 1L && no_na(r) && nrow(k) == 1L && !k$pass &&
      grepl("not run: no family file", k$check, fixed = TRUE) && any(r$rows$test == "Q3:F1")
  })
  check("failed reps beyond the bound fail the entry, whatever the survivors say", {
    d <- tree("failed", list(`families/F1.R` = paste(
      'family$methods$ml <- local({',
      '  ml <- family$methods$ml',
      '  function(sim, cell, seed) {',
      '    if (vs_seed_rep("F1", cell$cell, seed) %% 2 == 0) stop("injected failure")',
      '    ml(sim, cell, seed)',
      '  }',
      '})', sep = "\n")))
    r <- quick("failed", "--only=Q3:F1", "--reps=4", "--workers=1", at = d)
    k <- r$rows[grepl("failed reps", r$rows$check, fixed = TRUE), ]
    note(k$check)
    r$status == 1L && no_na(r) && nrow(k) == 1L && !k$pass && k$value == 0.5 &&
      grepl("injected failure", k$check, fixed = TRUE)
  })
}

# --- Summary -------------------------------------------------------------------------

v <- unlist(results)
cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(v), sum(!v),
            proc.time()[["elapsed"]] - t_start))
quit(save = "no", status = if (any(!v)) 1L else 0L)
