# ============================================================================
# The long tier's runner. Runs the units (family, cell, method, rep) of one
# family, or of several, into rep files, one per unit:
#   <out>/<run>/<family>/<cell>/<method>/rep_NNNNN.rds
# Each rep is simulated once and its data go to every method of the rep
# (paired comparisons; run_rep() in lib/fit.R). A unit whose rep file exists
# and reads is skipped (resume; an unreadable file is rerun), so a run is
# resumed or extended by running it again with the same --run;
# --retry-failed also reruns the units whose record failed. Shard k of K
# takes the reps with rep %% K == k - 1. Prints one
# line per rep. When sourced (the plumbing check does) it only defines
# functions.
#
# Usage (the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/mc.R --family=F3
#       [--cells=a,b] [--methods=x,y] [--reps=a:b] [--shard=k/K] [--threads=T]
#       [--run=TAG] [--out=DIR] [--retry-failed] [--ref=NAME] [--fixtures]
#   --family=all        every family (T* plumbing fixtures only with --fixtures)
#   --reps=a:b          reps a..b of every selected cell and method (default
#                       1:R_by_method[[cell]][method] where the family gives
#                       it, else 1:R$mc[cell])
#   --run=TAG           the run directory <out>/<TAG> (default a time stamp;
#                       required with --shard); --out defaults to
#                       _validation/output/, through the output guard
#   --warm              source the library, run its compile hooks (lib
#                       functions named *_compile, e.g. geo_gen_compile()),
#                       check the selection, print RUN_DIR=<dir>, and exit:
#                       run_shards.sh runs it once before starting shards
#   --check-seeds       list the pairs of cells (of --family, else of every
#                       family) whose seed ranges overlap; exit 1 if any
# ============================================================================

# --- Selecting units ----------------------------------------------------------

.mc_split <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",", fixed = TRUE)[[1]])

# --reps=a:b (or a single rep) as an integer vector within 1..99999.
.mc_reps <- function(x) {
  if (is.null(x)) return(NULL)
  ab <- suppressWarnings(as.integer(strsplit(x, ":", fixed = TRUE)[[1]]))
  if (length(ab) == 1L) ab <- c(ab, ab)
  if (length(ab) != 2L || anyNA(ab) || ab[1] < 1L || ab[2] < ab[1] || ab[2] > 99999L) {
    stop("--reps=", x, " must be a:b with 1 <= a <= b <= 99999", call. = FALSE)
  }
  seq.int(ab[1], ab[2])
}

# --shard=k/K as c(k, K).
.mc_shard <- function(x) {
  if (is.null(x)) return(c(1L, 1L))
  kK <- suppressWarnings(as.integer(strsplit(x, "/", fixed = TRUE)[[1]]))
  if (length(kK) != 2L || anyNA(kK) || kK[1] < 1L || kK[1] > kK[2]) {
    stop("--shard=", x, " must be k/K with 1 <= k <= K", call. = FALSE)
  }
  kK
}

# The families named by --family (a list named by id): "all" is every
# family file, the T* fixtures only with fixtures = TRUE.
mc_families <- function(L, family, fixtures = FALSE) {
  ids <- .mc_split(family)
  if (identical(ids, "all")) ids <- list_families(fixtures = fixtures)
  if (!length(ids)) stop("no families selected", call. = FALSE)
  fix <- ids[startsWith(ids, "T")]
  if (length(fix) && !fixtures) {
    stop("family ", paste(fix, collapse = ", "), " is a plumbing fixture: add --fixtures",
         call. = FALSE)
  }
  stats::setNames(lapply(ids, load_family, lib_env = L), ids)
}

# The default reps of a unit (mc tier): 1..R_by_method[[cell]][method] where
# the family gives it (methods of one cell run at different R on the same
# paired reps), else 1..R$<tier>[cell]; none when the family names neither.
.mc_default_reps <- function(fam, cell, method, tier = "mc") {
  n <- if (tier == "mc") fam$R_by_method[[cell]][method]
  if (is.null(n) || is.na(n)) n <- fam$R[[tier]][cell]
  if (is.null(n) || is.na(n)) integer(0) else seq_len(n)
}

# The selected units as a data.frame(family, cell, rep, method), in the order
# they run: family, cell, rep, then the cell's methods. --cells and
# --methods must each match something.
mc_units <- function(fams, cells = NULL, methods = NULL, reps = NULL, shard = c(1L, 1L)) {
  all_cells <- unlist(lapply(fams, function(f) f$cells$cell))
  all_methods <- unlist(lapply(fams, function(f) names(f$methods)))
  if (length(bad <- setdiff(cells, all_cells))) {
    stop("unknown cell(s) ", paste(bad, collapse = ", "), "; cells: ",
         paste(unique(all_cells), collapse = ", "), call. = FALSE)
  }
  if (length(bad <- setdiff(methods, all_methods))) {
    stop("unknown method(s) ", paste(bad, collapse = ", "), "; methods: ",
         paste(unique(all_methods), collapse = ", "), call. = FALSE)
  }
  out <- list()
  for (fam in fams) {
    for (cl in fam$cells$cell) {
      if (!is.null(cells) && !cl %in% cells) next
      ms <- cell_methods(fam, cl)
      if (!is.null(methods)) ms <- intersect(ms, methods)
      u <- do.call(rbind, lapply(ms, function(m) {
        rr <- reps %||% .mc_default_reps(fam, cl, m)
        rr <- rr[rr %% shard[2] == shard[1] - 1L]
        if (length(rr)) data.frame(rep = rr, method = m, k = match(m, ms))
      }))
      if (is.null(u)) next
      u <- u[order(u$rep, u$k), , drop = FALSE]
      out[[length(out) + 1L]] <- data.frame(family = fam$id, cell = cl, rep = u$rep,
                                            method = u$method, stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(family = character(), cell = character(), rep = integer(),
                      method = character(), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# --- Seed ranges ---------------------------------------------------------------

# The pairs of (family, cell) whose seed ranges overlap. Seeds are
# 1 + (h + rep) mod (2^31 - 2), so a cell's reps a..b take the cyclic
# interval of seeds starting at vs_seed(family, cell, a), of length b - a + 1;
# two intervals [o_i, o_i + n_i) and [o_j, o_j + n_j) mod M overlap when
# (o_j - o_i) mod M < n_i or (o_i - o_j) mod M < n_j.
mc_seed_overlaps <- function(fams, reps = NULL) {
  M <- 2^31 - 2
  rng <- do.call(rbind, lapply(fams, function(f) do.call(rbind, lapply(f$cells$cell, function(cl) {
    n <- if (is.null(reps)) f$R$mc[cl] else length(reps)
    if (is.null(n) || is.na(n) || n < 1) return(NULL)
    a <- if (is.null(reps)) 1 else min(reps)
    data.frame(family = f$id, cell = cl, first = a, last = a + n - 1,
               o = vs_seed(f$id, cl, a) - 1, n = as.numeric(n), stringsAsFactors = FALSE)
  }))))
  hits <- list()
  for (i in seq_len(max(0L, NROW(rng) - 1L))) {
    for (j in (i + 1L):nrow(rng)) {
      d <- (rng$o[j] - rng$o[i]) %% M
      if (d < rng$n[i] || (M - d) %% M < rng$n[j]) {
        hits[[length(hits) + 1L]] <- data.frame(a = i, b = j)
      }
    }
  }
  list(ranges = rng, overlaps = if (length(hits)) do.call(rbind, hits) else NULL)
}

# --- Main -------------------------------------------------------------------------

mc_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  t_start <- proc.time()[["elapsed"]]
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  root <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
  b <- boot$vs_boot(root)
  args <- vs_args(argv, defaults = list(
    family = NULL, cells = NULL, methods = NULL, reps = NULL, shard = NULL, threads = 1L,
    run = NULL, out = NULL, `retry-failed` = FALSE, ref = NULL, warm = FALSE,
    fixtures = FALSE, `check-seeds` = FALSE))
  reps <- .mc_reps(args$reps)
  shard <- .mc_shard(args$shard)

  if (args$`check-seeds`) {
    fams <- mc_families(b$L, args$family %||% "all", args$fixtures)
    s <- mc_seed_overlaps(fams, reps)
    rng <- s$ranges
    lab <- function(i) sprintf("%s/%s (reps %d-%d, seeds %d+)", rng$family[i], rng$cell[i],
                               rng$first[i], rng$last[i], as.integer(rng$o[i] + 1))
    cat(sprintf("seed ranges: %d families, %d cells\n", length(fams), NROW(rng)))
    if (is.null(s$overlaps)) {
      cat("no overlapping seed ranges\n")
      return(0L)
    }
    for (k in seq_len(nrow(s$overlaps))) {
      cat("overlap:", lab(s$overlaps$a[k]), "and", lab(s$overlaps$b[k]), "\n")
    }
    return(1L)
  }

  # The selection is checked before anything is compiled or written.
  if (is.null(args$family)) stop("--family is required (an id, a list a,b, or all)", call. = FALSE)
  if (shard[2] > 1L && is.null(args$run)) {
    stop("--shard needs --run, so that the shards write into one run", call. = FALSE)
  }
  fams <- mc_families(b$L, args$family, args$fixtures)
  units <- mc_units(fams, .mc_split(args$cells), .mc_split(args$methods), reps, shard)
  run <- args$run %||% format(Sys.time(), "%Y%m%d-%H%M%S")
  run_dir <- if (is.null(args$out)) vs_out_dir(run = run)
             else vs_out_dir(out = file.path(args$out, run))

  if (args$warm) {
    hooks <- Filter(function(f) is.function(get(f, envir = b$L)), ls(b$L, pattern = "_compile$"))
    for (h in hooks) {
      cat("warm: ", h, "()\n", sep = "")
      get(h, envir = b$L)()
    }
    if (!length(hooks)) cat("warm: the library defines no compile hook\n")
    if (!is.null(units)) {
      cat(sprintf("warm: %d units selected (%s)\n", nrow(units), paste(names(fams), collapse = ", ")))
    }
    cat("RUN_DIR=", run_dir, "\n", sep = "")
    return(0L)
  }

  vs_threads(args$threads)
  write_metadata(run_dir, args$ref)
  cat(sprintf("mc.R: choicer %s from %s; %d thread(s); shard %d/%d\n", b$build$version,
              b$build$path, choicer::thread_info()$omp_get_max_threads, shard[1], shard[2]))
  cat(sprintf("run %s: %d units (%s)\n", run_dir, nrow(units), paste(names(fams), collapse = ", ")))
  flush(stdout())

  # A unit runs when its rep file is missing or unreadable (a file cut short
  # by a kill), and with --retry-failed when its record failed.
  n_run <- n_failed <- n_skipped <- n_unreadable <- 0L
  groups <- split(seq_len(nrow(units)), paste(units$family, units$cell, units$rep, sep = "\r"))
  groups <- groups[order(vapply(groups, `[`, 1L, 1L))]
  for (g in groups) {
    fam <- fams[[units$family[g[1]]]]
    cl <- units$cell[g[1]]
    r <- units$rep[g[1]]
    paths <- vapply(units$method[g], function(m) rep_path(run_dir, fam$id, cl, m, r), "")
    state <- vapply(paths, function(p) {
      if (!file.exists(p)) return("missing")
      rec <- .vs_read_rep(p)
      if (is.null(rec)) "unreadable" else if (rec$ok) "ok" else "failed"
    }, "")
    todo <- state %in% c("missing", "unreadable") | (args$`retry-failed` & state == "failed")
    if (any(state == "unreadable")) {
      cat(sprintf("%s %s rep %d: unreadable rep file(s) of %s rerun\n", fam$id, cl, r,
                  paste(units$method[g][state == "unreadable"], collapse = ", ")))
      n_unreadable <- n_unreadable + sum(state == "unreadable")
    }
    n_skipped <- n_skipped + sum(!todo)
    if (!any(todo)) next
    t0 <- proc.time()[["elapsed"]]
    recs <- run_rep(fam, cl, r, units$method[g][todo], ref = args$ref)
    for (k in seq_along(recs)) write_rep(recs[[k]], paths[todo][k])
    ok <- vapply(recs, function(x) x$ok, logical(1))
    n_run <- n_run + length(recs)
    n_failed <- n_failed + sum(!ok)
    cat(sprintf("%s %s rep %d: %s (%.1f s)\n", fam$id, cl, r,
                paste(names(recs), ifelse(ok, "ok", "FAIL"), sep = " ", collapse = ", "),
                proc.time()[["elapsed"]] - t0))
    flush(stdout())
    rm(recs)
  }
  cat(sprintf(paste("done: %d units run (%d failed, %d unreadable files rerun),",
                    "%d skipped (rep files present), %.1f s\n"),
              n_run, n_failed, n_unreadable, n_skipped, proc.time()[["elapsed"]] - t_start))
  0L
}

if (sys.nframe() == 0L) quit(save = "no", status = mc_main())
