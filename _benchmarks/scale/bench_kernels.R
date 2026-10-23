#!/usr/bin/env Rscript
# One timing or memory run of the mixed logit's compiled kernels on synthetic
# population-scale data: one build (R_LIBS), one configuration, one process.
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/bench_kernels.R tag=<label> \
#     regime=<claims|census|census20k|small> U=<units> S=<draws> threads=1,11 reps=3 \
#     mode=<time|mem> draws=<generate|store> out=<file.rds> \
#     what=grad[,hess,bhhh,scores,tastes,predict,elas,elas_r,dr,dr_r,meat_cluster,meatr_robust,meatr_cluster,layout] \
#     [seed=1] [Kx=10] [Kw=3] [warm=1] [clusters=2000] [acc=<bytes>]
#
# mode=time runs without the shim: per (kernel, threads) a warm call (warm=1)
# and `reps` timed calls. mode=mem runs under the malloc shim (output/Rmc with
# MC_SHIM set; macOS): one call per (kernel, threads) with its heap blocks,
# bytes moved, peak live heap and sampled peak RSS.
#
# Draws: generate mode by default; draws=store builds the estimation cube
# (get_halton_normals(), one block per unit) only for estimation kernels, and
# store-mode prediction follows each build's own predict(): builds with
# store-point draws (0.3.0) form them on the fly, older ones take the
# per-situation cube. Kernels:
#   grad, hess, bhhh, scores, tastes, predict  the internal kernels behind the
#     log-likelihood and gradient, the Hessian, the BHHH matrix, the score
#     matrix, conditional_tastes() and predict();
#   elas, dr       the elasticities and diversion ratios of the first fixed
#     coefficient (elas_r, dr_r: of the first random one), with the prediction
#     draws; `acc` passes acc_bytes to builds that take it (0.3.0's J x J
#     accumulators; 0 forces the threads to share one result);
#   meat_cluster   the clustered meat kernel (0.3.0), `clusters` clusters of units;
#   meatr_robust, meatr_cluster  the meats as releases before 0.3.0 formed them:
#     the U x n score matrix, then crossprod() (and rowsum() by cluster);
#   layout         (0.3.0) the gradient call that options(choicer.max_memory = 1)
#     stops right after its layout is built: the per-call layout cost.
# Every kernel is evaluated at the 0.2.x cold start (theta_start()); the
# result's check value (the objective, or a sum) must agree across builds.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()
ns <- asNamespace("choicer")

opt <- scale_args()
get_opt <- function(k, default) opt[[k]] %||% default
tag <- get_opt("tag", basename(lib))
regime <- get_opt("regime", "claims")
U <- as.integer(get_opt("U", "2000"))
S <- as.integer(get_opt("S", "100"))
threads <- as.integer(strsplit(get_opt("threads", "11"), ",")[[1]])
reps <- as.integer(get_opt("reps", "3"))
mode <- get_opt("mode", "time")
what <- strsplit(get_opt("what", "grad"), ",")[[1]]
seed <- as.numeric(get_opt("seed", "1"))
Kx <- as.integer(get_opt("Kx", "10"))
Kw <- as.integer(get_opt("Kw", "3"))
warm <- get_opt("warm", "1") == "1"
G <- as.integer(get_opt("clusters", "2000"))
draws_mode <- get_opt("draws", "generate")
acc <- get_opt("acc", "")
out <- get_opt("out", file.path(scale_output("kernels"),
                                sprintf("%s_%s_U%d_S%d_%s.rds", tag, regime, U, S, mode)))
stopifnot(mode %in% c("time", "mem"), draws_mode %in% c("generate", "store"))
if (mode == "mem" && !shim_loaded()) stop("mode=mem needs the malloc shim (output/Rmc with MC_SHIM set)")
if (mode == "time" && shim_loaded()) stop("mode=time must run without the shim")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

store_points <- exists(".MXL_STORE_POINTS", envir = ns, inherits = FALSE)
cat(sprintf("== %s (%s) %s | %s U=%d S=%d Kx=%d Kw=%d draws=%s | threads %s | what %s | load %s\n",
            tag, if (store_points) "store points" else "cube", mode, regime, U, S, Kx, Kw, draws_mode,
            paste(threads, collapse = ","), paste(what, collapse = ","), loadavg()))

rss0 <- rss_now()
t_gen <- system.time(d <- gen_regime(regime, U, seed = seed, K_x = Kx, K_w = Kw))[["elapsed"]]
rss_data <- rss_now()
describe(d)
n_units <- length(d$Ti)
theta <- theta_start(d)
cl_u <- as.integer((seq_len(n_units) - 1L) %% G + 1L)

est_kernels <- c("grad", "hess", "bhhh", "scores", "tastes", "meat_cluster",
                 "meatr_robust", "meatr_cluster")
empty <- array(0, dim = c(d$K_w, 0L, 0L))
eta_est <- empty
t_draws <- NA_real_; draws_alloc <- NA_real_; draws_rss <- NA_real_
if (draws_mode == "store" && any(what %in% est_kernels)) {
  rss_b <- rss_now()[["rss"]]; mc_reset(); rss_start(5)
  t_draws <- system.time(eta_est <- get_halton_normals(S, n_units, d$K_w))[["elapsed"]]
  m_d <- rss_stop(); a_d <- mc_get(); mc_disable()
  draws_alloc <- a_d[["bytes_alloc"]]; draws_rss <- m_d[["peak_rss"]] - rss_b
  cat(sprintf("  estimation cube %d x %d x %d (%.2f GB) in %.1fs, %.2f GB allocated, peak RSS +%.2f GB\n",
              d$K_w, S, n_units, 8 * d$K_w * S * n_units / 1e9, t_draws, draws_alloc / 1e9, draws_rss / 1e9))
}
est_flags <- list(rc_correlation = TRUE, rc_mean = TRUE, use_asc = TRUE,
                  include_outside_option = FALSE,
                  gen_seed = if (draws_mode == "store") -1L else 17L, gen_scramble = 1L,
                  gen_S = if (draws_mode == "store") 0L else S)
pred <- list(eta_draws = empty, gen_seed = 17L, gen_scramble = 1L, gen_S = S)
t_pd <- NA_real_
if (draws_mode == "store" && any(c("predict", "elas", "elas_r", "dr", "dr_r") %in% what)) {
  if (store_points) {
    pred <- list(eta_draws = empty, gen_seed = 0L, gen_scramble = 2L, gen_S = S)
  } else {
    rss_b <- rss_now()[["rss"]]; mc_reset(); rss_start(5)
    t_pd <- system.time(eta_p <- get_halton_normals(S, length(d$M), d$K_w))[["elapsed"]]
    m_d <- rss_stop(); a_d <- mc_get(); mc_disable()
    pred <- list(eta_draws = eta_p, gen_seed = -1L, gen_scramble = 1L, gen_S = 0L)
    cat(sprintf("  per-situation cube %d x %d x %d (%.2f GB) in %.1fs, %.2f GB allocated, peak RSS +%.2f GB\n",
                d$K_w, S, length(d$M), 8 * d$K_w * S * length(d$M) / 1e9, t_pd,
                a_d[["bytes_alloc"]] / 1e9, (m_d[["peak_rss"]] - rss_b) / 1e9))
  }
}
est_args <- c(list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
                   choice_idx = d$choice_idx, M = d$M, eta_draws = eta_est,
                   rc_dist = rep(0L, d$K_w)), est_flags, list(Ti = d$Ti))
weighted <- c(est_args, list(weights = d$weights))
# The elasticity and diversion kernels (prediction draws), and acc_bytes for
# builds that take it
pred_kernel <- function(f, random, with_choice) {
  a <- list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx)
  if (with_choice) a$choice_idx <- d$choice_idx
  a <- c(a, list(M = d$M, weights = d$weights, eta_draws = pred$eta_draws,
                 rc_dist = rep(0L, d$K_w), elast_var_idx = 1L, is_random_coef = random),
         est_flags[c("rc_correlation", "rc_mean", "use_asc", "include_outside_option")],
         pred[c("gen_seed", "gen_scramble", "gen_S")])
  if (nzchar(acc) && "acc_bytes" %in% names(formals(f))) a$acc_bytes <- as.numeric(acc)
  do.call(f, a)
}
call_kernel <- function(k) {
  switch(k,
    grad = do.call(ns$mxl_loglik_gradient_parallel, weighted),
    hess = do.call(ns$mxl_hessian_parallel, weighted),
    bhhh = do.call(ns$mxl_bhhh_parallel, weighted),
    scores = do.call(ns$mxl_scores_parallel, est_args),
    tastes = do.call(ns$mxl_conditional_tastes_parallel, est_args),
    meat_cluster = do.call(ns$mxl_cluster_meat_parallel, c(weighted, list(cluster = cl_u))),
    meatr_robust = { Sm <- do.call(ns$mxl_scores_parallel, est_args); crossprod(Sm) },
    meatr_cluster = { Sm <- do.call(ns$mxl_scores_parallel, est_args); crossprod(rowsum(Sm, cl_u)) },
    predict = do.call(ns$mxl_predict,
                      c(list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx, M = d$M,
                             eta_draws = pred$eta_draws, rc_dist = rep(0L, d$K_w)),
                        est_flags[c("rc_correlation", "rc_mean", "use_asc", "include_outside_option")],
                        pred[c("gen_seed", "gen_scramble", "gen_S")])),
    elas = pred_kernel(ns$mxl_elasticities_parallel, FALSE, TRUE),
    elas_r = pred_kernel(ns$mxl_elasticities_parallel, TRUE, TRUE),
    dr = pred_kernel(ns$mxl_diversion_ratios_parallel, FALSE, FALSE),
    dr_r = pred_kernel(ns$mxl_diversion_ratios_parallel, TRUE, FALSE),
    layout = {
      if (!store_points) stop("layout needs the memory check of choicer 0.3.0")
      old <- options(choicer.max_memory = 1)
      on.exit(options(old), add = TRUE)
      msg <- tryCatch({ do.call(ns$mxl_loglik_gradient_parallel, weighted); "no stop" },
                      error = conditionMessage)
      if (!grepl("needs about", msg, fixed = TRUE)) stop("the layout probe did not stop at the memory check: ", msg)
      msg
    },
    stop("unknown kernel ", k))
}
check_value <- function(x) {
  if (is.character(x)) return(NA_real_)
  if (is.list(x) && !is.null(x$objective)) return(x$objective)
  if (is.numeric(x)) return(sum(x))
  sum(vapply(x, function(e) sum(as.numeric(e)), 0))
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
info <- list(tag = tag, lib = lib, store_points = store_points, mode = mode, regime = regime, acc = acc,
             U = U, N = d$N, rows = d$rows, J = d$J, n_par = length(theta), S = S, Kx = Kx,
             Kw = Kw, seed = seed, draws = draws_mode, t_gen = t_gen, t_draws = t_draws,
             draws_alloc = draws_alloc, draws_rss = draws_rss, t_pred_cube = t_pd,
             rss0 = rss0, rss_data = rss_data, max_unit_rows = max(unit_rows(d)),
             load_start = loadavg(), started = format(Sys.time()))

for (k in what) {
  for (nt in threads) {
    set_num_threads(nt)
    invisible(gc())
    if (mode == "mem") {
      mc_reset(); rss_start(2)
      t0 <- Sys.time()
      res <- tryCatch(call_kernel(k), error = function(e) e)
      t_meas <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      m <- rss_stop(); a <- mc_get(); mc_disable()
      err <- if (inherits(res, "error")) conditionMessage(res) else ""
      cv <- if (nzchar(err)) NA_real_ else check_value(res)
      rm(res); invisible(gc())
      add_row(what = k, threads = nt, rep = 0L, sec = t_meas, check = cv,
              n_alloc = a[["n_alloc"]], bytes_alloc = a[["bytes_alloc"]], peak_live = a[["peak_live"]],
              peak_rss = m[["peak_rss"]] - rss_data[["rss"]], error = err)
      cat(sprintf("  %-13s threads=%2d mem: %.3fs | %.0f blocks %.2f GB, peak live %.3f GB, peak RSS +%.2f GB%s\n",
                  k, nt, t_meas, a[["n_alloc"]], a[["bytes_alloc"]] / 1e9, a[["peak_live"]] / 1e9,
                  (m[["peak_rss"]] - rss_data[["rss"]]) / 1e9, if (nzchar(err)) paste(" ERROR:", err) else ""))
      next
    }
    cv <- NA_real_; err <- ""
    if (warm) {
      res <- tryCatch(call_kernel(k), error = function(e) e)
      if (inherits(res, "error")) err <- conditionMessage(res) else cv <- check_value(res)
      rm(res); invisible(gc())
    }
    tt <- numeric(0)
    if (!nzchar(err)) for (r in seq_len(reps)) {
      t0 <- Sys.time()
      res <- tryCatch(call_kernel(k), error = function(e) e)
      tt[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (inherits(res, "error")) { err <- conditionMessage(res); tt <- tt[-r]; break }
      if (is.na(cv)) cv <- check_value(res)
      rm(res)
      add_row(what = k, threads = nt, rep = r, sec = tt[r], check = cv, n_alloc = NA_real_,
              bytes_alloc = NA_real_, peak_live = NA_real_, peak_rss = NA_real_, error = "")
    }
    if (nzchar(err)) add_row(what = k, threads = nt, rep = NA_integer_, sec = NA_real_, check = NA_real_,
                             n_alloc = NA_real_, bytes_alloc = NA_real_, peak_live = NA_real_,
                             peak_rss = NA_real_, error = err)
    cat(sprintf("  %-13s threads=%2d time: med %.4fs [min %.4f, n %d] check %.10g%s\n", k, nt,
                if (length(tt)) stats::median(tt) else NA, if (length(tt)) min(tt) else NA,
                length(tt), cv, if (nzchar(err)) paste(" ERROR:", substr(err, 1, 200)) else ""))
  }
}

res <- if (length(rows)) do.call(rbind, rows) else data.frame()
info$load_end <- loadavg()
info$finished <- format(Sys.time())
attr(res, "info") <- info
saveRDS(res, out)
cat(sprintf("  saved %s | load %s -> %s\n", out, info$load_start, info$load_end))
