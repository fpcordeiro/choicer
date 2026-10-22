#!/usr/bin/env Rscript
# Multinomial and nested logit fits on synthetic population-scale data: one
# build (R_LIBS), one configuration, one process. For each se_method, one
# run_mnlogit() / run_nestlogit() call on prepared data (advanced pathway),
# timed in full and split into the optimizer's time (fit$optimizer$elapsed_time,
# the "bhhh" scale pass included) and the rest (the fit-time variance, mostly);
# then, on the first fit, each post-hoc vcov(fit, type = ) timed.
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/fit_mnl_nl.R tag=<label> \
#     model=<mnl|nl> regime=<claims|census|census20k|small> U=<units> threads=11 \
#     se=hessian,bhhh,sandwich,cluster[,numeric] posthoc=hessian,bhhh,robust,cluster \
#     [nests=50] [weights=unit|runif] [seed=1] [Kx=10] [out=<file.rds>]
#
# Data as bench_mnl_nl.R (gen_regime()'s design read as a cross-section; NL
# nests a partition of the alternatives), in long form with one cluster per
# decision maker. Alternatives never chosen, and the situations offering one
# chosen wherever it is offered, are not removed: the warning on constants
# without a finite estimate is recorded, not avoided (time is what is measured).
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()

opt <- scale_args()
get_opt <- function(k, default) opt[[k]] %||% default
tag <- get_opt("tag", basename(lib))
model <- get_opt("model", "mnl")
regime <- get_opt("regime", "census")
U <- as.integer(get_opt("U", "20000"))
nt <- as.integer(get_opt("threads", "11"))
ses <- strsplit(get_opt("se", "hessian,bhhh,sandwich,cluster"), ",")[[1]]
posthoc <- strsplit(get_opt("posthoc", "hessian,bhhh,robust,cluster"), ",")[[1]]
seed <- as.numeric(get_opt("seed", "1"))
Kx <- as.integer(get_opt("Kx", "10"))
n_nests <- as.integer(get_opt("nests", "50"))
wmode <- get_opt("weights", "unit")
out <- get_opt("out", file.path(scale_output("fits"),
                                sprintf("%s_%s_%s_U%d_t%d.rds", tag, model, regime, U, nt)))
stopifnot(model %in% c("mnl", "nl"), wmode %in% c("unit", "runif"))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
set_num_threads(nt)

g <- gen_regime(regime, U, seed = seed, K_x = Kx, K_w = 1L)
N <- length(g$M)
sit <- rep.int(seq_len(N), g$M)
chosen <- integer(length(sit))
first <- c(0L, cumsum(g$M))[seq_len(N)]
chosen[first + g$choice_idx] <- 1L
set.seed(seed + 1)
w_i <- if (wmode == "runif") runif(N, 0.5, 2) else rep(1, N)
dt <- data.table::data.table(id = sit, alt = g$alt_idx, choice = chosen,
                             w = w_i[sit], clus = rep.int(seq_along(g$Ti), g$Ti)[sit])
xn <- paste0("x", seq_len(Kx))
for (k in seq_len(Kx)) data.table::set(dt, j = xn[k], value = g$X[, k])
if (model == "nl") dt[, nest := (alt - 1L) %% n_nests + 1L]
rm(g); invisible(gc())
cat(sprintf("== %s %s fits | %s U=%d situations=%d rows=%.3g | threads %d | load %s\n",
            tag, model, regime, U, N, nrow(dt), nt, loadavg()))

t_prep <- system.time(inp <- if (model == "mnl") {
  prepare_mnl_data(dt, "id", "alt", "choice", xn, weights_col = "w", cluster_col = "clus")
} else {
  prepare_nl_data(dt, "id", "alt", "choice", xn, "nest", weights_col = "w",
                  cluster_col = "clus")
})[["elapsed"]]
rm(dt); invisible(gc())
cat(sprintf("  prepared in %.1fs: %d parameters\n", t_prep,
            ncol(inp$X) + (if (model == "nl") sum(tabulate(inp$nest_idx) > 1L) else 0L) +
              nrow(inp$alt_mapping) - 1L))

fit_one <- function(se) {
  warns <- character()
  t <- system.time(fit <- withCallingHandlers(
    if (model == "mnl") run_mnlogit(input_data = inp, se_method = se) else
      run_nestlogit(input_data = inp, se_method = se),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage")))[["elapsed"]]
  list(fit = fit, t = t, warns = warns)
}

rows <- list(); first_fit <- NULL
for (se in ses) {
  r <- fit_one(se)
  f <- r$fit
  t_opt <- f$optimizer$elapsed_time %||% NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "fit", what = se, sec = r$t, sec_opt = t_opt, sec_rest = r$t - t_opt,
    iterations = f$optimizer$iterations %||% NA, convergence = f$convergence,
    loglik = f$loglik, se_ok = !is.null(f$vcov), n_par = length(f$coefficients),
    warnings = length(r$warns), stringsAsFactors = FALSE)
  cat(sprintf("  fit %-9s %8.2fs (optimizer %.2fs, rest %.2fs) it %s conv %s loglik %.6f vcov %s\n",
              se, r$t, t_opt, r$t - t_opt, f$optimizer$iterations, f$convergence, f$loglik,
              if (is.null(f$vcov)) "NULL" else "ok"))
  if (is.null(first_fit)) first_fit <- f
  rm(f, r); invisible(gc())
}
for (ty in posthoc) {
  t <- system.time(v <- tryCatch(vcov(first_fit, type = ty), error = function(e) e))[["elapsed"]]
  err <- if (inherits(v, "error")) conditionMessage(v) else ""
  rows[[length(rows) + 1L]] <- data.frame(
    kind = "posthoc", what = ty, sec = t, sec_opt = NA_real_, sec_rest = NA_real_,
    iterations = NA, convergence = NA, loglik = NA_real_, se_ok = !nzchar(err),
    n_par = length(first_fit$coefficients), warnings = 0L, stringsAsFactors = FALSE)
  cat(sprintf("  vcov(type = %-10s %8.2fs%s\n", sprintf('"%s")', ty), t,
              if (nzchar(err)) paste(" ERROR:", substr(err, 1, 160)) else ""))
}
res <- do.call(rbind, rows)
attr(res, "info") <- list(tag = tag, lib = lib, model = model, regime = regime, U = U,
                          N = N, threads = nt, t_prep = t_prep, weights = wmode,
                          nests = if (model == "nl") n_nests, load_end = loadavg(),
                          finished = format(Sys.time()))
saveRDS(res, out)
cat(sprintf("  saved %s\n", out))
