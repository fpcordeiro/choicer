# ============================================================================
# The validation suite's shared infrastructure: the build guard, command-line
# arguments, the output guard, seeds, rep files, run metadata, threads and
# memory. Definitions only: sourcing this file has no side effects. Scripts
# source it through vs_source_lib() (lib/families.R), which also sets
# `.vs_root`, the `_validation/` directory, in the library's environment.
# ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Build guard ------------------------------------------------------------

# Loads choicer from the first entry of R_LIBS and checks that the namespace
# came from there and is a development build: the site library's 0.2.1
# lacks person_col, polish and concentrate_out_asc, and shadows a dev build
# unless the scratch library comes first. Returns the path and version.
vs_load_choicer <- function(min_version = "0.2.1.9000") {
  how <- paste0("Run as R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript <script>, where <lib> holds ",
                "a development build of choicer (>= ", min_version, ").")
  lib <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1]]
  lib <- lib[nzchar(lib)]
  if (!length(lib)) stop("R_LIBS is not set, so the build under test is not explicit. ", how,
                         call. = FALSE)
  if (!dir.exists(lib[1])) stop("R_LIBS[1] = ", lib[1], " does not exist. ", how, call. = FALSE)
  lib <- normalizePath(lib[1])
  if (!file.exists(file.path(lib, "choicer", "DESCRIPTION"))) {
    stop("R_LIBS[1] = ", lib, " holds no choicer. ", how, call. = FALSE)
  }
  suppressMessages(library(choicer, lib.loc = lib))
  path <- getNamespaceInfo("choicer", "path")
  if (normalizePath(dirname(path)) != lib) {
    stop("choicer was loaded from ", path, ", not from R_LIBS[1] = ", lib,
         " (another copy was loaded first). ", how, call. = FALSE)
  }
  version <- getNamespaceVersion("choicer")[["version"]]
  if (package_version(version) < min_version) {
    stop("choicer ", version, " at ", path, " is older than ", min_version,
         ": it lacks person_col, polish and concentrate_out_asc (the site library's ",
         "release?). ", how, call. = FALSE)
  }
  invisible(list(path = path, version = version, lib = lib))
}

# Stops naming the reference packages that are not installed: the quick tier
# needs survival, mlogit and bayesm, the mc tier data.table and Rcpp.
vs_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing package(s): ", paste(missing, collapse = ", "),
         ". Install them into a library on R_LIBS or .libPaths().", call. = FALSE)
  }
  invisible(TRUE)
}

# --- Arguments --------------------------------------------------------------

# Parses --key=value and bare --flag arguments against `defaults`, a named
# list: unknown keys are an error listing the known ones. A value takes the
# type of its default (logical, integer, double; character otherwise or for a
# NULL default); a bare flag is TRUE and needs a logical or NULL default.
vs_args <- function(args = commandArgs(trailingOnly = TRUE), defaults = list()) {
  known <- names(defaults)
  if (length(defaults) && (is.null(known) || any(!nzchar(known)))) {
    stop("vs_args(): every default needs a name")
  }
  out <- defaults
  for (a in args) {
    if (!startsWith(a, "--") || nchar(a) < 3L) {
      stop("Unexpected argument '", a, "': use --key=value or --flag", call. = FALSE)
    }
    kv <- substring(a, 3L)
    eq <- regexpr("=", kv, fixed = TRUE)
    key <- if (eq > 0L) substr(kv, 1L, eq - 1L) else kv
    if (!key %in% known) {
      stop("Unknown argument --", key, "; known: ",
           if (length(known)) paste0("--", known, collapse = ", ") else "none", call. = FALSE)
    }
    out[key] <- list(if (eq > 0L) .vs_arg_value(substring(kv, eq + 1L), defaults[[key]], key)
                     else .vs_arg_flag(defaults[[key]], key))
  }
  out
}

.vs_arg_value <- function(value, default, key) {
  if (is.null(default) || is.character(default) || length(default) != 1L) return(value)
  v <- if (is.logical(default)) {
    switch(tolower(value), "true" = , "yes" = , "1" = TRUE, "false" = , "no" = , "0" = FALSE, NA)
  } else if (is.integer(default)) {
    if (grepl("^-?[0-9]+$", value)) suppressWarnings(as.integer(value)) else NA_integer_
  } else if (is.double(default)) {
    suppressWarnings(as.numeric(value))
  } else {
    value
  }
  if (is.na(v)) {
    stop("--", key, "=", value, " is not a valid ", typeof(default), call. = FALSE)
  }
  v
}

.vs_arg_flag <- function(default, key) {
  if (is.null(default) || is.logical(default)) return(TRUE)
  stop("--", key, " needs a value (--", key, "=...)", call. = FALSE)
}

# --- Output guard -----------------------------------------------------------

# The `_validation/` directory, set by vs_source_lib().
vs_root <- function() {
  root <- get0(".vs_root", envir = environment(vs_root))
  if (is.null(root)) stop("`.vs_root` is unset: source the library with vs_source_lib()")
  root
}

# An absolute path whether or not it exists yet (symbolic links resolved in
# the part that exists).
.vs_abs_path <- function(p) {
  p <- path.expand(p)
  if (!grepl("^/", p) && !grepl("^[A-Za-z]:", p)) p <- file.path(getwd(), p)
  rest <- character()
  while (!dir.exists(p)) {
    rest <- c(basename(p), rest)
    p <- dirname(p)
  }
  do.call(file.path, as.list(c(normalizePath(p), rest)))
}

# The run's output directory, created: `out` itself, or by default
# `_validation/output/<run>`. Inside a package tree (an ancestor holds a
# DESCRIPTION) only a directory git ignores is accepted, checked on a probe
# path below it with `git check-ignore -q` (exit 0 ignored, 1 not ignored,
# 128 not a repository): the suite never writes into tracked files.
vs_out_dir <- function(out = NULL, run = NULL, create = TRUE) {
  if (is.null(out)) {
    if (is.null(run) || !grepl("^[A-Za-z0-9._-]+$", run)) {
      stop("vs_out_dir(): name the run (letters, digits, '.', '_', '-') or pass `out`")
    }
    out <- file.path(vs_root(), "output", run)
  }
  out <- .vs_abs_path(out)
  tree <- out
  repeat {
    if (file.exists(file.path(tree, "DESCRIPTION"))) break
    if (dirname(tree) == tree) { tree <- NULL; break }
    tree <- dirname(tree)
  }
  if (!is.null(tree)) {
    anc <- out
    while (!dir.exists(anc)) anc <- dirname(anc)
    probe <- file.path(out, ".vs_probe")
    st <- suppressWarnings(tryCatch(
      system2("git", c("-C", shQuote(anc), "check-ignore", "-q", "--", shQuote(probe)),
              stdout = FALSE, stderr = FALSE),
      error = function(e) -1L))
    if (!identical(as.integer(st), 0L)) {
      stop("Refusing to write into ", out, ": it lies inside the package tree ", tree,
           " and git does not ignore it (git check-ignore exit ", st, "). ",
           "Use the default _validation/output/<run> or a directory outside the tree.",
           call. = FALSE)
    }
  }
  if (create) dir.create(out, showWarnings = FALSE, recursive = TRUE)
  out
}

# --- Seeds ------------------------------------------------------------------

# A replication's seed, shared by every method of the rep (paired
# comparisons: they fit the same simulated data):
#   seed = 1 + (h + rep) mod (2^31 - 2),   h = fmix32(fnv1a32("<family>|<cell>")),
# FNV-1a over the UTF-8 bytes, then MurmurHash3's 32-bit finalizer, in exact
# double arithmetic (no packages, the same on every platform). Seeds are
# integers in [1, 2^31 - 2]; reps 1..2^31 - 2 of one (family, cell) get
# distinct seeds. Vectorized over `rep`.
vs_seed <- function(family, cell, rep) {
  if (!is.character(family) || length(family) != 1L || !is.character(cell) ||
      length(cell) != 1L) {
    stop("vs_seed(): `family` and `cell` must be single strings")
  }
  if (!is.numeric(rep) || anyNA(rep) || any(rep < 1) || any(rep != floor(rep))) {
    stop("vs_seed(): `rep` must be positive whole numbers")
  }
  h <- .vs_fmix32(.vs_fnv1a32(paste0(family, "|", cell)))
  as.integer(1 + (h + rep) %% (2^31 - 2))
}

# The rep a seed belongs to: the inverse of vs_seed() for one (family, cell).
vs_seed_rep <- function(family, cell, seed) {
  as.integer((as.numeric(seed) - vs_seed(family, cell, 1)) %% (2^31 - 2) + 1)
}

# Unsigned 32-bit arithmetic on doubles: values below 2^32, every product
# split into 16-bit halves so that it stays below 2^53.
.vs_u32_xor <- function(a, b) {
  bitwXor(a %/% 65536, b %/% 65536) * 65536 + bitwXor(a %% 65536, b %% 65536)
}
.vs_u32_mul <- function(a, b) {
  al <- a %% 65536; ah <- a %/% 65536; bl <- b %% 65536; bh <- b %/% 65536
  (al * bl + ((ah * bl + al * bh) %% 65536) * 65536) %% 2^32
}
.vs_fnv1a32 <- function(s) {
  h <- 2166136261
  for (b in as.integer(charToRaw(enc2utf8(s)))) h <- .vs_u32_mul(.vs_u32_xor(h, b), 16777619)
  h
}
.vs_fmix32 <- function(h) {
  h <- .vs_u32_xor(h, h %/% 65536)
  h <- .vs_u32_mul(h, 2246822507)            # 0x85ebca6b
  h <- .vs_u32_xor(h, h %/% 8192)
  h <- .vs_u32_mul(h, 3266489909)            # 0xc2b2ae35
  .vs_u32_xor(h, h %/% 65536)
}

# Seeds R's generator with the default kinds named, so that a user's
# RNGkind() cannot change a replication.
vs_set_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
}

# Evaluates `expr` under `seed`, then restores the caller's generator state.
vs_with_seed <- function(seed, expr) {
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old, envir = globalenv())
  })
  vs_set_seed(seed)
  expr
}

# --- Rep files --------------------------------------------------------------

.vs_check_part <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || !grepl("^[A-Za-z0-9._+-]+$", x)) {
    stop(what, " must be a single name of letters, digits, '.', '_', '+' or '-'")
  }
  x
}

# <out>/<family>/<cell>/<method>/rep_NNNNN.rds
rep_path <- function(out, family, cell, method, rep) {
  if (!is.numeric(rep) || length(rep) != 1L || is.na(rep) || rep < 1 || rep > 99999 ||
      rep != floor(rep)) {
    stop("rep_path(): `rep` must be a whole number in 1..99999")
  }
  file.path(out, .vs_check_part(family, "family"), .vs_check_part(cell, "cell"),
            .vs_check_part(method, "method"), sprintf("rep_%05d.rds", as.integer(rep)))
}

# Saves `obj` to `path` atomically: a temporary file in the same directory,
# then a rename, so that a reader or a resumed run never sees a partial file.
.vs_save_atomic <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  ok <- tryCatch({ saveRDS(obj, tmp); file.rename(tmp, path) },
                 error = function(e) { unlink(tmp); stop(e) })
  if (!ok) {
    unlink(tmp)
    stop("could not rename ", tmp, " to ", path)
  }
  invisible(path)
}

write_rep <- function(record, path) .vs_save_atomic(record, path)

rep_exists <- function(path) file.exists(path)

# The rep files below `out`, parsed from their paths (nothing is read).
.vs_rep_files <- function(out) {
  rel <- list.files(out, pattern = "^rep_[0-9]{5}\\.rds$", recursive = TRUE)
  parts <- strsplit(rel, "/", fixed = TRUE)
  keep <- lengths(parts) == 4L
  parts <- parts[keep]
  data.frame(
    family = vapply(parts, `[`, "", 1L),
    cell   = vapply(parts, `[`, "", 2L),
    method = vapply(parts, `[`, "", 3L),
    rep    = as.integer(substr(vapply(parts, `[`, "", 4L), 5L, 9L)),
    path   = file.path(out, rel[keep]),
    stringsAsFactors = FALSE
  )
}

.vs_read_rep <- function(path) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.list(r) && is.logical(r$ok) && length(r$ok) == 1L) r else NULL
}

# The records of one family (optionally one cell and one method), ordered by
# cell, method and rep. Unreadable files (corrupt, partial, not a record) are
# skipped with one warning of class "vs_unreadable_rep" and listed in the
# result's "unreadable" attribute.
read_reps <- function(out, family, cell = NULL, method = NULL) {
  idx <- .vs_rep_files(out)
  sel <- idx$family == family & (is.null(cell) | idx$cell %in% cell) &
    (is.null(method) | idx$method %in% method)
  idx <- idx[sel, , drop = FALSE]
  idx <- idx[order(idx$cell, idx$method, idx$rep, method = "radix"), , drop = FALSE]
  recs <- lapply(idx$path, .vs_read_rep)
  bad <- vapply(recs, is.null, logical(1))
  if (any(bad)) {
    warning(warningCondition(
      paste0("read_reps(): ", sum(bad), " unreadable rep file(s): ",
             paste(idx$path[bad], collapse = ", ")), class = "vs_unreadable_rep"))
  }
  recs <- recs[!bad]
  attr(recs, "unreadable") <- idx$path[bad]
  recs
}

# One row per rep file below `out`: family, cell, method, rep, ok (NA for an
# unreadable file; not read when read = FALSE) and path.
rep_index <- function(out, read = TRUE) {
  idx <- .vs_rep_files(out)
  idx <- idx[order(idx$family, idx$cell, idx$method, idx$rep, method = "radix"), , drop = FALSE]
  idx$ok <- if (read) {
    vapply(idx$path, function(p) { r <- .vs_read_rep(p); if (is.null(r)) NA else r$ok },
           logical(1), USE.NAMES = FALSE)
  } else {
    rep(NA, nrow(idx))
  }
  rownames(idx) <- NULL
  idx[c("family", "cell", "method", "rep", "ok", "path")]
}

# --- Metadata ---------------------------------------------------------------

# A command's first output line, trimmed (NA on failure or no output); with
# all = TRUE every line (character(0) for none, NA on failure).
.vs_cmd <- function(cmd, args, all = FALSE) {
  out <- suppressWarnings(tryCatch(system2(cmd, args, stdout = TRUE, stderr = FALSE),
                                   error = function(e) NULL))
  if (is.null(out) || !is.null(attr(out, "status"))) return(NA_character_)
  if (all) return(out)
  if (!length(out)) NA_character_ else trimws(out[1])
}

# The paths that differ from HEAD: tracked files anywhere in the tree, and
# untracked files under `root` (`_validation/`) that git does not ignore,
# relative to the repository's root (porcelain status). NA when git fails.
.vs_git_dirty <- function(root) {
  st <- function(...) .vs_cmd("git", c("-C", shQuote(root), "status", "--porcelain", ...),
                              all = TRUE)
  tracked <- st("--untracked-files=no")
  here <- st("--untracked-files=all", "--", ".")
  if (anyNA(tracked) || anyNA(here)) return(NA_character_)
  unique(substring(c(tracked, here), 4L))
}

.vs_cpu <- function() {
  if (Sys.info()[["sysname"]] == "Darwin") {
    return(.vs_cmd("sysctl", c("-n", "machdep.cpu.brand_string")))
  }
  if (file.exists("/proc/cpuinfo")) {
    m <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
    if (length(m)) return(trimws(sub("^[^:]*:", "", m[1])))
  }
  NA_character_
}

# The run's provenance: choicer's path and version, the --ref, the session,
# BLAS/LAPACK, choicer's threads, the CPU, host, time, and the git SHA of the
# tree holding `_validation/`, dirty when a tracked file differs from it or
# an untracked, not ignored file lies under `_validation/` (new lib or family
# files), with the first 20 such paths.
vs_metadata <- function(ref = NULL) {
  root <- get0(".vs_root", envir = environment(vs_metadata))
  loaded <- isNamespaceLoaded("choicer")
  sha <- if (is.null(root)) NA_character_
         else .vs_cmd("git", c("-C", shQuote(root), "rev-parse", "HEAD"))
  paths <- if (is.na(sha)) NA_character_ else .vs_git_dirty(root)
  dirty <- if (anyNA(paths)) NA else length(paths) > 0L
  list(
    choicer_path    = if (loaded) getNamespaceInfo("choicer", "path") else NA_character_,
    choicer_version = if (loaded) getNamespaceVersion("choicer")[["version"]] else NA_character_,
    ref             = ref,
    session         = utils::sessionInfo(),
    blas            = extSoftVersion()[["BLAS"]],
    lapack          = La_library(),
    lapack_version  = La_version(),
    threads         = if (loaded) choicer::thread_info() else NULL,
    cpu             = .vs_cpu(),
    host            = Sys.info()[["nodename"]],
    pid             = Sys.getpid(),
    time            = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_sha         = sha,
    git_dirty       = dirty,
    git_dirty_paths = utils::head(paths, 20L)
  )
}

# Writes vs_metadata() to <out>/meta_<host>_<pid>.rds (one file per process,
# so shards never collide) and returns its path.
write_metadata <- function(out, ref = NULL) {
  m <- vs_metadata(ref)
  host <- gsub("[^A-Za-z0-9._-]", "_", m$host)
  .vs_save_atomic(m, file.path(out, sprintf("meta_%s_%d.rds", host, m$pid)))
}

# --- Threads ----------------------------------------------------------------

# Sets choicer's OpenMP threads; returns the previous count, invisibly. A
# build without OpenMP runs at one thread whatever is asked.
vs_threads <- function(n) {
  info <- choicer::thread_info()
  if (!isTRUE(info$openmp_enabled)) return(invisible(1L))
  prev <- info$omp_get_max_threads
  choicer::set_num_threads(as.integer(n))
  invisible(prev)
}

# --- Memory -----------------------------------------------------------------

# Starts a peak-memory window: resets R's "max used" (gc(reset = TRUE)) and,
# on Linux, the process's peak RSS (writing 5 to /proc/self/clear_refs;
# failure ignored). Returns whether the RSS peak was reset.
vs_mem_reset <- function() {
  gc(reset = TRUE)
  ok <- file.exists("/proc/self/clear_refs") &&
    isTRUE(tryCatch({ cat("5", file = "/proc/self/clear_refs"); TRUE },
                    error = function(e) FALSE, warning = function(w) FALSE))
  invisible(list(rss_reset = ok))
}

# The window's peak memory in MB: `heap`, R's heap peak (gc()'s "max used";
# blind to C++ allocations such as Armadillo's), and `rss`, the process's
# resident set: on Linux its peak (VmHWM: kind "peak" since the reset, or
# "lifetime_peak" where it failed), on macOS the current RSS from ps (kind
# "current": the OS keeps no resettable peak).
vs_mem_peak <- function(reset = NULL) {
  g <- gc()
  heap <- sum(g[, which(colnames(g) == "max used") + 1L])
  rss <- NA_real_
  kind <- NA_character_
  if (file.exists("/proc/self/status")) {
    st <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
    if (length(st)) {
      rss <- as.numeric(gsub("[^0-9]", "", st[1])) / 1024
      kind <- if (isTRUE(reset$rss_reset)) "peak" else "lifetime_peak"
    }
  } else {
    kb <- suppressWarnings(as.numeric(.vs_cmd("ps", c("-o", "rss=", "-p", Sys.getpid()))))
    if (!is.na(kb)) {
      rss <- kb / 1024
      kind <- "current"
    }
  }
  list(heap_mb = heap, rss_mb = rss, rss_kind = kind)
}
