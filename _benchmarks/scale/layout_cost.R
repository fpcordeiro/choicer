#!/usr/bin/env Rscript
# The per-call cost of the mixed logit kernels' layout against a likelihood
# evaluation, for a server run at 10^9 rows.
#
# Every estimation kernel call builds the layout of the stacked design
# (MxlLayout, src/mxlogit.cpp): it fills N + 1 eight-byte row offsets and makes
# one pass over the rows and three serial passes over the situations. The
# serial passes grow with the number of situations, so their share of an
# evaluation is largest for many small choice sets, few draws and many threads.
# A layout cache would pay only if that share exceeds about 5% (README.md).
# This script measures it on a cross-section of 4-alternative choice sets with
# one fixed and one random coefficient (the layout does not depend on the
# number of covariates; X and W take 16 bytes a row, 16 GB at 10^9 rows): the
# layout as the gradient call that options(choicer.max_memory = 1) stops right
# after building it (choicer 0.3.0's memory check), against full gradient calls.
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/layout_cost.R [--units=25000000] \
#     [--threads=<all>] [--S=25,100] [--reps=5] [--panel] [--out=<file.rds>]
#
# --units is the number of choice situations (4 rows each; the default 2.5e7 is
# 10^8 rows; 2.5e8, 10^9 rows, peaks near 40 GB while the data are generated
# and holds about 27 GB after). --panel passes one decision maker per situation
# (Ti), as a panel fit of single choices would, which adds the units' offsets.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()
ns <- asNamespace("choicer")
if (!exists(".MXL_STORE_POINTS", envir = ns, inherits = FALSE)) {
  stop("this build has no memory check to stop at (choicer 0.3.0 or later is needed)")
}
opt <- scale_args()
U <- as.numeric(opt$units %||% "2.5e7")
threads <- as.integer(strsplit(opt$threads %||% as.character(physical_cores()$n), ",")[[1]])
S_v <- as.integer(strsplit(opt$S %||% "25,100", ",")[[1]])
reps <- as.integer(opt$reps %||% "5")
panel <- !is.null(opt$panel)
out <- opt$out %||% file.path(scale_output("layout"), sprintf("layout_U%.0f%s.rds", U, if (panel) "_panel" else ""))

cat(sprintf("== layout cost | %.0f situations (%.3g rows) | threads %s | S %s | load %s\n",
            U, 4 * U, paste(threads, collapse = ","), paste(S_v, collapse = ","), loadavg()))
t_gen <- system.time(d <- gen_regime("small", as.integer(U), K_x = 1L, K_w = 1L,
                                     threads = max(threads)))[["elapsed"]]
d$gamma <- NULL  # the realized tastes, not needed here
invisible(gc())
cat(sprintf("  generated %.0f situations, %.3g rows in %.0f s\n", length(d$M), sum(as.numeric(d$M)), t_gen))
args <- list(theta = theta_start(d), X = d$X, W = d$W, alt_idx = d$alt_idx,
             choice_idx = d$choice_idx, M = d$M, weights = d$weights,
             eta_draws = array(0, dim = c(1L, 0L, 0L)), rc_dist = 0L, rc_correlation = TRUE,
             rc_mean = TRUE, use_asc = TRUE, include_outside_option = FALSE, gen_seed = 17L,
             gen_scramble = 1L, Ti = if (panel) d$Ti else NULL)
layout_call <- function() {
  old <- options(choicer.max_memory = 1)
  on.exit(options(old))
  msg <- tryCatch({ do.call(ns$mxl_loglik_gradient_parallel, c(args, list(gen_S = 10L))); "no stop" },
                  error = conditionMessage)
  if (!grepl("needs about", msg, fixed = TRUE)) stop("the call did not stop at the memory check: ", msg)
}
timed <- function(f, n) vapply(seq_len(n), function(i) system.time(f(), gcFirst = FALSE)[["elapsed"]], 0)

rows <- list()
for (nt in threads) {
  set_num_threads(nt)
  invisible(gc())
  layout_call()
  t_lay <- timed(layout_call, reps)
  for (S in S_v) {
    grad <- function() do.call(ns$mxl_loglik_gradient_parallel, c(args, list(gen_S = S)))
    grad()
    t_grad <- timed(grad, max(1L, reps %/% 2L))
    rows[[length(rows) + 1L]] <- data.frame(
      threads = nt, S = S, situations = U, rows = 4 * U, layout_s = stats::median(t_lay),
      gradient_s = stats::median(t_grad), share = stats::median(t_lay) / stats::median(t_grad),
      offsets_gb = 8 * (U + 1) * (1 + panel) / 1e9, panel = panel, load = loadavg(),
      stringsAsFactors = FALSE)
    cat(sprintf("  %3d threads, S = %4d: layout %.4f s, gradient %.3f s, share %.2f%%\n",
                nt, S, stats::median(t_lay), stats::median(t_grad),
                100 * stats::median(t_lay) / stats::median(t_grad)))
  }
}
res <- do.call(rbind, rows)
attr(res, "info") <- list(lib = lib, t_gen = t_gen, R = R.version.string,
                          blas = extSoftVersion()[["BLAS"]], time = format(Sys.time()))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
saveRDS(res, out)
cat("saved", out, "\n")
