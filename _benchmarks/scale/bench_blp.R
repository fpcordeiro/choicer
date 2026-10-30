#!/usr/bin/env Rscript
# The share inversion's methods on synthetic population-scale data: one build
# (R_LIBS), one configuration, one process (BLP engine plan, Phase B, R4).
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/bench_blp.R tag=<label> \
#     model=<mnl|nl|mxl> regime=<claims|census|census20k|small> U=<units> \
#     methods=simple,original,anderson,squarem threads=1,11 reps=3 \
#     [mode=time|mem] [tol=1e-8] [max_iter=1000] [S=20] [nests=50] [ioo=0] \
#     [shift=0.5] [start=fit|zero] [seed=1] [Kx=10] [pass=log|exp]
#     [compensate=off|on|auto] [out=<file.rds>]
#
# Data: gen_regime()'s design read as a cross-section of choice situations, as
# in bench_mnl_nl.R: X (Kx columns) and the generator's constants; the nested
# logit's nests as there (alternative j in nest (j - 1) %% nests + 1,
# dissimilarity 0.7 in every nest that holds two or more; at most J / 2 nests,
# so that some nest does); the mixed logit's
# three correlated normal random coefficients (gen_regime()'s W, mu and L),
# their draws formed on the fly (S per situation, seeded), kept across
# evaluations within blp()'s default budget of 2^32 bytes. The targets are
# counterfactual, the model's shares at delta* = delta + N(0, shift^2) (the
# generator's constants moved), and the inversion starts from the generator's
# constants (start=fit, as blp() starts from a fit's) or from zeros.
#
# Each method runs once untimed (warm), then reps times per thread count
# through the internal solvers (*_blp_solve() with .blp_control()): a row
# records the status, evaluations of the shares, iterations, Berry's
# contraction joining (race) and restarts, the residual, the largest error
# against delta* (each constant relative to the first's without an outside
# option), the time and the time per evaluation. mode=mem runs each method
# once under the malloc shim (heap blocks, bytes and peak live heap; see
# make_rmc.sh). Methods a build lacks (releases before Phase B have no
# "anderson" or "squarem") record an error instead of running another.
# pass=exp evaluates from the exp-space cache (Phase C; .blp_control(pass =));
# a build without that setting runs its own passes, and its rows say so
# (pass "log"); exp_kept says whether a model's pass kept the cache (NA: the
# model's pass has no exp-space form in this build). compensate= sets the
# passes' compensated share sums (Phase C, .blp_control(compensate =)) where
# the build has them.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()
ns <- asNamespace("choicer")

opt <- scale_args()
get_opt <- function(k, default) opt[[k]] %||% default
tag <- get_opt("tag", basename(lib))
model <- get_opt("model", "mnl")
regime <- get_opt("regime", "census")
U <- as.integer(get_opt("U", "20000"))
methods <- strsplit(get_opt("methods", "simple,original,anderson"), ",")[[1]]
threads <- as.integer(strsplit(get_opt("threads", "1,11"), ",")[[1]])
reps <- as.integer(get_opt("reps", "3"))
mode <- get_opt("mode", "time")
tol <- as.numeric(get_opt("tol", "1e-8"))
max_iter <- as.integer(get_opt("max_iter", "1000"))
S <- as.integer(get_opt("S", "20"))
n_nests <- as.integer(get_opt("nests", "50"))
ioo <- get_opt("ioo", "0") == "1"
shift <- as.numeric(get_opt("shift", "0.5"))
start_at <- get_opt("start", "fit")
seed <- as.numeric(get_opt("seed", "1"))
Kx <- as.integer(get_opt("Kx", "10"))
pass <- get_opt("pass", "log")
compensate <- get_opt("compensate", "off")
out <- get_opt("out", file.path(scale_output("blp"),
                                sprintf("%s_%s_%s_U%d_%s.rds", tag, model, regime, U, mode)))
stopifnot(model %in% c("mnl", "nl", "mxl"), mode %in% c("time", "mem"),
          start_at %in% c("fit", "zero"), pass %in% c("log", "exp"),
          compensate %in% c("off", "on", "auto"))
# (a build without a setting runs its own passes)
if (!"pass" %in% names(formals(ns$.blp_control))) pass <- "log"
if (!"compensate" %in% names(formals(ns$.blp_control))) compensate <- "off"

if (mode == "mem" && !shim_loaded()) stop("mode=mem needs the malloc shim (output/Rmc with MC_SHIM set)")
if (mode == "time" && shim_loaded()) stop("mode=time must run without the shim")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

cat(sprintf("== %s %s %s | %s U=%d Kx=%d%s%s | methods %s | threads %s | tol %g | pass %s | compensate %s | load %s\n",
            tag, model, mode, regime, U, Kx,
            switch(model, nl = sprintf(" nests=%d", n_nests), mxl = sprintf(" S=%d", S), ""),
            if (ioo) " outside option" else "", paste(methods, collapse = ","),
            paste(threads, collapse = ","), tol, pass, compensate, loadavg()))

K_w <- if (model == "mxl") 3L else 1L
t_gen <- system.time(g <- gen_regime(regime, U, seed = seed, K_x = Kx, K_w = K_w))[["elapsed"]]
d <- list(X = g$X, alt_idx = g$alt_idx, M = g$M, weights = rep(1, length(g$M)))
J <- g$J
tp <- true_par(Kx, K_w)
beta <- tp$beta
# the free constants of delta (J), as the kernels take them
free <- function(delta) if (ioo) delta else delta[-1L] - delta[1L]
delta0 <- if (ioo) g$delta else g$delta - g$delta[1L]
set.seed(seed + 2)
delta_star <- delta0 + stats::rnorm(J, 0, shift)
if (!ioo) delta_star <- delta_star - delta_star[1L]
start <- if (start_at == "fit") delta0 else rep(0, J)
nest_idx <- NULL; lambda_full <- NULL
if (model == "nl") {
  n_nests <- min(n_nests, max(1L, J %/% 2L))  # (some nest holds two or more)
  nest_idx <- as.integer((seq_len(J) - 1L) %% n_nests + 1L)
  big <- tabulate(nest_idx, n_nests) > 1L
  lambda_full <- ifelse(big, 0.7, 1)
}
if (model == "mxl") {
  Lp <- unlist(lapply(seq_len(K_w), function(i) sapply(seq_len(i), function(j)
    if (i == j) log(tp$L[i, i]) else tp$L[i, j])))
  empty <- array(0, c(K_w, 0L, 0L))
  gen <- list(gen_seed = as.integer(seed), gen_scramble = 1L, gen_S = S)
}
target <- as.numeric(switch(model,
  mnl = ns$mnl_predict_shares(c(beta, free(delta_star)), d$X, d$alt_idx, d$M, d$weights,
                              TRUE, ioo),
  nl = ns$nl_predict_shares(c(beta, lambda_full[big], free(delta_star)), d$X, d$alt_idx, d$M,
                            d$weights, nest_idx, TRUE, ioo),
  mxl = do.call(ns$mxl_predict_shares, c(list(
    c(beta, tp$mu, Lp, free(delta_star)), d$X, g$W, d$alt_idx, d$M, d$weights, empty,
    rep(0L, K_w), TRUE, TRUE, TRUE, ioo), gen))))
N <- length(d$M); rows <- length(d$alt_idx)
cat(sprintf("[%s] situations=%.0f rows=%.3g J=%d | M: mean %.1f max %d | generated in %.1fs\n",
            regime, N, rows, J, mean(d$M), max(d$M), t_gen))

has_method <- function(m) {
  m %in% c("simple", "original") ||
    (exists(".BLP_SETTINGS", envir = ns) && m %in% names(get(".BLP_SETTINGS", envir = ns)))
}
solve_blp <- function(m) {
  extra <- c(if (pass != "log") list(pass = pass),
             if (compensate != "off") list(compensate = compensate))
  ctl <- do.call(ns$.blp_control, c(list(tol, max_iter, method = m), extra))
  switch(model,
    mnl = ns$mnl_blp_solve(start, target, d$X, beta, d$alt_idx, d$M, d$weights, ioo, ctl),
    nl = ns$nl_blp_solve(start, target, d$X, beta, lambda_full, d$alt_idx, nest_idx, d$M,
                         d$weights, ioo, ctl),
    mxl = do.call(ns$mxl_blp_solve, c(list(
      start, target, d$X, g$W, beta, tp$mu, Lp, d$alt_idx, d$M, d$weights, empty,
      rep(0L, K_w), TRUE, TRUE, ioo, ctl), gen, list(cache_bytes = 2^32))))
}
delta_err <- function(delta) {
  delta <- as.numeric(delta)
  if (!ioo) delta <- delta - delta[1L]
  max(abs(delta - delta_star))
}

rows_out <- list()
add_row <- function(...) rows_out[[length(rows_out) + 1L]] <<- data.frame(
  tag = tag, model = model, regime = regime, U = U, N = N, rows = rows, J = J, S = S,
  ioo = ioo, tol = tol, max_iter = max_iter, start = start_at, pass = pass,
  compensate = compensate, ...,
  load = loadavg(),
  stringsAsFactors = FALSE)
record <- function(m, nt, rep, sec, r, a = NULL) {
  add_row(method = m, threads = nt, rep = rep, sec = sec,
          code = r$code, evals = r$evaluations, iterations = r$iterations,
          race = r$race %||% NA_integer_, restarts = r$restarts %||% NA_integer_,
          exp_kept = r$exp_kept %||% NA, exp_evals = r$exp_evals %||% NA_real_,
          exp_forms = r$exp_forms %||% NA_real_, exp_floors = r$exp_floors %||% NA_real_,
          exp_prefix = r$exp_prefix %||% NA_real_, exp_reads = r$exp_reads %||% NA_real_,
          passes = r$passes %||% NA_real_, reads = r$reads %||% NA_real_,
          residual = r$residual, err = delta_err(r$delta), sec_per_eval = sec / r$evaluations,
          n_alloc = if (is.null(a)) NA_real_ else a[["n_alloc"]],
          peak_live = if (is.null(a)) NA_real_ else a[["peak_live"]], error = "")
}
for (m in methods) {
  if (!has_method(m)) {
    add_row(method = m, threads = NA_integer_, rep = NA_integer_, sec = NA_real_, code = NA_integer_,
            evals = NA_integer_, iterations = NA_integer_, race = NA_integer_,
            restarts = NA_integer_, exp_kept = NA, exp_evals = NA_real_,
            exp_forms = NA_real_, exp_floors = NA_real_, exp_prefix = NA_real_,
            exp_reads = NA_real_, passes = NA_real_, reads = NA_real_,
            residual = NA_real_, err = NA_real_, sec_per_eval = NA_real_,
            n_alloc = NA_real_, peak_live = NA_real_, error = "method not in this build")
    cat(sprintf("  %-9s not in this build\n", m))
    next
  }
  for (nt in threads) {
    set_num_threads(nt)
    invisible(gc())
    if (mode == "mem") {
      mc_reset()
      t0 <- Sys.time()
      r <- suppressWarnings(solve_blp(m))
      sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      a <- mc_get(); mc_disable()
      record(m, nt, 0L, sec, r, a)
      cat(sprintf("  %-9s threads=%2d mem: %d evaluations %.2fs | %.0f blocks, peak live %.3f GB\n",
                  m, nt, r$evaluations, sec, a[["n_alloc"]], a[["peak_live"]] / 1e9))
      next
    }
    invisible(suppressWarnings(solve_blp(m)))  # warm
    tt <- numeric(0)
    for (rp in seq_len(reps)) {
      invisible(gc())
      t0 <- Sys.time()
      r <- suppressWarnings(solve_blp(m))
      tt[rp] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      record(m, nt, rp, tt[rp], r)
    }
    cat(sprintf("  %-9s threads=%2d code %d | %4d evaluations | time med %.3fs min %.3fs | %.4fs per evaluation | residual %.2e err %.2e\n",
                m, nt, r$code, r$evaluations, stats::median(tt), min(tt),
                min(tt) / r$evaluations, r$residual, delta_err(r$delta)))
  }
}
res <- do.call(rbind, rows_out)
saveRDS(res, out)
cat("saved", out, "\n")
