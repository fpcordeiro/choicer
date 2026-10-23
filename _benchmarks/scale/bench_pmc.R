#!/usr/bin/env Rscript
# Instructions retired, cycles and time per call of the multinomial, nested or
# mixed logit kernels at one thread, for one build (macOS: the process's own
# counters, pmc.c, no root). Instructions per situation (per decision maker in
# the mixed logit) tell more work from slower execution, where wall times on a
# loaded machine cannot; instructions per cycle tell the rest.
#
#   KMP_BLOCKTIME=0 R_LIBS=<lib> Rscript _benchmarks/scale/bench_pmc.R tag=<label> \
#     models=mnl,nl[,mxl] kernels=grad,bhhh[,grad_opg,hess,robust,predict,shares,elas,dr] \
#     [regime=small] [U=2500000] [Kx=10] [nests=2] [S=100] [draws=store|generate] \
#     [reps=5] [out=<file.rds>]
#
# KMP_BLOCKTIME=0 keeps libomp's idle workers from spinning into the count.
# Data as bench_mnl_nl.R (NL: alternative j in nest (j - 1) %% nests + 1,
# dissimilarity 0.7) and bench_kernels.R (MXL: K_w = 3 correlated normal
# coefficients with means, its 0.2.x start). The gradient, whose work no Phase 8
# change touched until 8B, is the control. Expect 1-2% between runs.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()
ns <- asNamespace("choicer")

so <- file.path(scale_output(), "pmc.so")
if (!file.exists(so)) {
  st <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "SHLIB", "-o", shQuote(so), shQuote(file.path(SCALE_DIR, "pmc.c"))),
                stdout = FALSE, stderr = FALSE)
  if (st != 0L || !file.exists(so)) stop("could not build pmc.so (macOS only)")
}
dyn.load(so)
pmc <- function() .Call("pmc_read")

opt <- scale_args()
get_opt <- function(k, default) opt[[k]] %||% default
tag <- get_opt("tag", basename(lib))
models <- strsplit(get_opt("models", "mnl,nl"), ",")[[1]]
kernels <- strsplit(get_opt("kernels", "grad,bhhh"), ",")[[1]]
regime <- get_opt("regime", "small")
U <- as.integer(get_opt("U", "2500000"))
Kx <- as.integer(get_opt("Kx", "10"))
n_nests <- as.integer(get_opt("nests", "2"))
S <- as.integer(get_opt("S", "100"))
draws <- get_opt("draws", "store")
reps <- as.integer(get_opt("reps", "5"))
out <- get_opt("out", file.path(scale_output("pmc"), sprintf("%s_%s_U%d_K%d.rds", tag, regime, U, Kx)))
stopifnot(draws %in% c("store", "generate"))
cat(sprintf("== %s | %s U=%d Kx=%d | models %s | kernels %s | load %s\n", tag, regime, U, Kx,
            paste(models, collapse = ","), paste(kernels, collapse = ","), loadavg()))

set_num_threads(1L)
rows <- list()
run_model <- function(model) {
  if (model == "mxl") {
    d <- gen_regime(regime, U, seed = 1, K_x = Kx, K_w = 3L)
    theta <- theta_start(d)
    n_units <- length(d$Ti)
    eta <- array(0, dim = c(d$K_w, 0L, 0L))
    if (draws == "store") eta <- get_halton_normals(S, n_units, d$K_w)
    a <- c(list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
                M = d$M, eta_draws = eta, rc_dist = rep(0L, d$K_w)),
           list(rc_correlation = TRUE, rc_mean = TRUE, use_asc = TRUE, include_outside_option = FALSE,
                gen_seed = if (draws == "store") -1L else 17L, gen_scramble = 1L,
                gen_S = if (draws == "store") 0L else S),
           list(Ti = d$Ti, weights = d$weights))
    calls <- list(grad = function() do.call(ns$mxl_loglik_gradient_parallel, a),
                  bhhh = function() do.call(ns$mxl_bhhh_parallel, a),
                  hess = function() do.call(ns$mxl_hessian_parallel, a))
    return(list(calls = calls, units = n_units, n_par = length(theta)))
  }
  g <- gen_regime(regime, U, seed = 1, K_x = Kx, K_w = 1L)
  d <- list(X = g$X, alt_idx = g$alt_idx, choice_idx = g$choice_idx, M = g$M,
            weights = rep(1, length(g$M)))
  J <- g$J
  beta <- 0.5 * seq(-1, 1, length.out = Kx)
  delta <- g$delta[-1]
  rm(g); invisible(gc())
  nest_idx <- NULL
  if (model == "nl") {
    nest_idx <- as.integer((seq_len(J) - 1L) %% n_nests + 1L)
    theta <- c(beta, rep(0.7, sum(tabulate(nest_idx) > 1L)), delta)
  } else {
    theta <- c(beta, delta)
  }
  w_args <- function(w) {
    a <- list(theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx)
    if (model == "nl") a$nest_idx <- nest_idx
    c(a, list(M = d$M, weights = w, use_asc = TRUE, include_outside_option = FALSE))
  }
  p_args <- function(w = NULL) {
    a <- list(theta = theta, X = d$X, alt_idx = d$alt_idx, M = d$M)
    if (model == "nl") a$nest_idx <- nest_idx
    if (!is.null(w)) a$weights <- w
    c(a, list(use_asc = TRUE, include_outside_option = FALSE))
  }
  kern <- function(name) get(paste0(model, "_", name), envir = ns)
  calls <- list(
    grad = function() do.call(kern("loglik_gradient_parallel"), w_args(d$weights)),
    grad_opg = function() do.call(kern("loglik_gradient_parallel"),
                                  c(w_args(d$weights), list(opg_diag = TRUE))),
    bhhh = function() do.call(kern("bhhh_parallel"), w_args(d$weights)),
    robust = function() do.call(kern("bhhh_parallel"), w_args(d$weights^2)),
    hess = function() do.call(kern("loglik_hessian_parallel"), w_args(d$weights)),
    predict = function() do.call(kern("predict"), p_args()),
    shares = function() do.call(kern("predict_shares"), p_args(d$weights)),
    elas = function() do.call(kern("elasticities_parallel"),
                              c(w_args(d$weights), list(elast_var_idx = 1L))),
    dr = function() do.call(kern("diversion_ratios_parallel"), p_args(d$weights)))
  list(calls = calls, units = length(d$M), n_par = length(theta))
}

for (model in models) {
  m <- run_model(model)
  for (k in kernels) {
    f <- m$calls[[k]]
    if (is.null(f)) next
    res <- f(); rm(res); invisible(gc()); Sys.sleep(0.3)  # warm
    for (r in seq_len(reps)) {
      invisible(gc())
      p0 <- pmc(); t0 <- proc.time()[["elapsed"]]
      res <- f()
      t1 <- proc.time()[["elapsed"]]; p1 <- pmc()
      chk <- if (is.list(res)) res[[1L]][[1L]] else sum(res)
      rm(res)
      rows[[length(rows) + 1L]] <- data.frame(
        tag = tag, model = model, kernel = k, rep = r, Kx = Kx, units = m$units,
        n_par = m$n_par, sec = t1 - t0, instr = p1[1] - p0[1], cycles = p1[2] - p0[2],
        check = chk)
    }
  }
}
res <- do.call(rbind, rows)
attr(res, "info") <- list(tag = tag, lib = lib, regime = regime, U = U, Kx = Kx, nests = n_nests,
                          S = S, draws = draws, load = loadavg(), finished = format(Sys.time()))
saveRDS(res, out)
agg <- aggregate(cbind(sec, instr, cycles) ~ tag + model + kernel + Kx + n_par + units,
                 data = res, FUN = stats::median)
agg$instr_unit <- agg$instr / agg$units
agg$cyc_unit <- agg$cycles / agg$units
agg$ipc <- agg$instr / agg$cycles
print(agg[, c("tag", "model", "kernel", "Kx", "n_par", "sec", "instr_unit", "cyc_unit", "ipc")],
      digits = 5, row.names = FALSE)
cat(sprintf("  saved %s\n", out))
