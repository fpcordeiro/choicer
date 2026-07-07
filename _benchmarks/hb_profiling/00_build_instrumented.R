#!/usr/bin/env Rscript
# Build the HB-profiling instrumented copy of choicer.
#
# What this does (and does NOT do):
#   1. rsync's the current repo source tree into a scratch build directory
#      (never inside the repo -- default is the workspace tmp/ area used by
#      the profiling run; override with --build-dir).
#   2. Applies _benchmarks/hb_profiling/patches/hb_profile.patch to the COPY
#      ONLY. The patch adds `#ifdef HB_PROFILE` timer blocks to
#      src/hb_internal.h, src/hmnlogit.cpp, src/hmnprobit.cpp, and adds
#      `PKG_CPPFLAGS = -DHB_PROFILE` to src/Makevars(.win) in the copy.
#   3. R CMD INSTALLs the patched copy into a throwaway library
#      (--lib-dir, default <build-dir>/templib).
#
# The package source under version control is NEVER modified: every edit
# lives in the patch file, applied only to an rsync'd copy outside the repo.
# Verify with `git status --porcelain` in the repo after running this.
#
# Usage:
#   Rscript _benchmarks/hb_profiling/00_build_instrumented.R \
#     [--build-dir DIR] [--lib-dir DIR] [--repo-root DIR]
#
# The harness scripts (10_run_config.R) load the instrumented build via
# `library(choicer, lib.loc = <lib-dir>)` in a fresh subprocess, never in a
# session that also touches the ordinarily-installed choicer.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}

this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
default_repo_root <- if (length(this_file) == 1) {
  normalizePath(file.path(dirname(this_file), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
repo_root <- normalizePath(get_arg("--repo-root", default_repo_root),
                           mustWork = TRUE)

build_dir <- get_arg(
  "--build-dir",
  "/Users/fernando/.claude/plugins/data/statsclaw-statsclaw/workspace/choicer/tmp/hb_prof_build"
)
lib_dir <- get_arg("--lib-dir", file.path(build_dir, "templib"))
src_copy <- file.path(build_dir, "choicer_src")
patch_file <- file.path(repo_root, "_benchmarks", "hb_profiling", "patches",
                        "hb_profile.patch")

stopifnot(file.exists(patch_file))
if (startsWith(normalizePath(build_dir, mustWork = FALSE),
               normalizePath(repo_root, mustWork = TRUE))) {
  stop("build-dir must NOT be inside the repo tree: ", build_dir)
}

message("repo_root  = ", repo_root)
message("build_dir  = ", build_dir)
message("lib_dir    = ", lib_dir)

dir.create(src_copy, recursive = TRUE, showWarnings = FALSE)
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. rsync the repo source into the scratch copy (excluding VCS/build
# artifacts and the non-CRAN dev trees) ---------------------------------------
rsync_excludes <- c(
  ".git", "*.o", "*.so", ".Rproj.user", "_benchmarks", "_validation",
  "data_room", ".claude", "docs"
)
excl_args <- paste(sprintf("--exclude='%s'", rsync_excludes), collapse = " ")
rsync_cmd <- sprintf("rsync -a %s %s/ %s/", excl_args,
                     shQuote(repo_root, type = "sh"),
                     shQuote(src_copy, type = "sh"))
message("Syncing source: ", rsync_cmd)
status <- system(rsync_cmd)
if (status != 0) stop("rsync failed with status ", status)

# --- 2. Apply the instrumentation patch to the COPY only ---------------------
# Idempotency: skip if hmnlogit.cpp already has the HB_PROFILE marker (e.g.
# a repeat run without a fresh rsync target).
already_patched <- any(grepl(
  "HB_PROFILE",
  readLines(file.path(src_copy, "src", "hmnlogit.cpp"), warn = FALSE)
))
if (already_patched) {
  message("Copy already carries the HB_PROFILE patch; skipping `patch`.")
} else {
  patch_cmd <- sprintf("cd %s && patch -p0 < %s",
                       shQuote(src_copy, type = "sh"),
                       shQuote(patch_file, type = "sh"))
  message("Applying patch: ", patch_cmd)
  status <- system(patch_cmd)
  if (status != 0) stop("patch failed with status ", status)
}

# --- 3. R CMD INSTALL the instrumented copy into the throwaway library ------
install_cmd <- sprintf(
  "R CMD INSTALL -l %s --no-docs %s",
  shQuote(lib_dir, type = "sh"),
  shQuote(file.path(src_copy), type = "sh")
)
message("Installing: ", install_cmd)
status <- system(install_cmd)
if (status != 0) stop("R CMD INSTALL failed with status ", status)

message("Done. Load with: library(choicer, lib.loc = ", shQuote(lib_dir), ")")
