#!/usr/bin/env Rscript
# One timing or memory run of the multinomial or nested logit's compiled
# kernels on synthetic population-scale data: one build (R_LIBS), one
# configuration, one process. Output rows and options follow bench_kernels.R.
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/bench_mnl_nl.R tag=<label> \
#     model=<mnl|nl> regime=<claims|census|census20k|small> U=<units> \
#     threads=1,11 reps=3 mode=<time|mem> out=<file.rds> \
#     what=grad[,grad_opg,hess,bhhh,robust,scores,meat_cluster,meatr_robust,meatr_cluster,elas,dr,predict,shares,blp] \
#     [nests=50] [clusters=unit|<G>] [weights=unit|runif] [situations=<n>] [seed=1] [Kx=10] [warm=1] \
#     [ioo=0] [acc=<bytes>] [blp_iter=3]
#
# Data: gen_regime()'s design (the mixed logit's generator; its choices come
# from that model), read as a cross-section of choice situations: X (Kx
# columns) and a constant for every alternative but the first, which the
# design always offers some situation. NL adds `nests` nests, a partition of
# the alternatives (alternative j in nest (j - 1) %% nests + 1, so the
# alternatives of a region spread over many nests), dissimilarity 0.7 in every
# nest that holds two or more. The kernels are evaluated at the generator's
# coefficients and constants. `situations = n` keeps the first n situations
# (every parameter stays): for timing per situation the dense kernels of
# releases before 0.3.0, whose calls at 17,000 parameters would take days.
# Clusters: the decision maker of each situation (`unit`, the default), or
# `G` groups of decision makers.
#
# `ioo=1` adds the outside option (a constant for every alternative). `acc`
# passes acc_bytes to the elasticity and diversion kernels of builds that take
# it (0.3.0's J x J accumulators; 0 forces the threads to share one result).
#
# Kernels:
#   grad          the log-likelihood and gradient;
#   grad_opg      the same with the OPG diagonal (scale_vars = "bhhh"'s pass);
#   hess, bhhh    the Hessian and the BHHH matrix (weights w);
#   robust        the BHHH kernel with weights w^2 (the sandwich's meat);
#   scores        the N x n score matrix;
#   meat_cluster  the clustered meat kernel (0.3.0);
#   meatr_robust, meatr_cluster  the meats as releases before 0.3.0 formed them
#                 post hoc: the score matrix, then crossprod(w * S) (and
#                 rowsum() by cluster);
#   elas, dr      elasticities (of the first covariate) and diversion ratios;
#   predict, shares  probabilities and utilities, and predicted shares;
#   blp           blp_contraction() / nl_blp_contraction() at tol = 0, from
#                 zeros to the shares at theta: blp_iter evaluations of the
#                 shares (0.3.0's engine; before it, blp_iter iterations, so
#                 blp_iter + 1 passes).
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
threads <- as.integer(strsplit(get_opt("threads", "11"), ",")[[1]])
reps <- as.integer(get_opt("reps", "3"))
mode <- get_opt("mode", "time")
what <- strsplit(get_opt("what", "grad"), ",")[[1]]
seed <- as.numeric(get_opt("seed", "1"))
Kx <- as.integer(get_opt("Kx", "10"))
warm <- get_opt("warm", "1") == "1"
n_nests <- as.integer(get_opt("nests", "50"))
clusters <- get_opt("clusters", "unit")
wmode <- get_opt("weights", "unit")
n_sit <- as.numeric(get_opt("situations", "0"))
ioo <- get_opt("ioo", "0") == "1"
acc <- get_opt("acc", "")
blp_iter <- as.integer(get_opt("blp_iter", "3"))
out <- get_opt("out", file.path(scale_output("kernels"),
                                sprintf("%s_%s_%s_U%d_%s.rds", tag, model, regime, U, mode)))
stopifnot(model %in% c("mnl", "nl"), mode %in% c("time", "mem"),
          wmode %in% c("unit", "runif"))
if (mode == "mem" && !shim_loaded()) stop("mode=mem needs the malloc shim (output/Rmc with MC_SHIM set)")
if (mode == "time" && shim_loaded()) stop("mode=time must run without the shim")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

cat(sprintf("== %s %s %s | %s U=%d Kx=%d%s%s | weights %s clusters %s%s | threads %s | what %s | load %s\n",
            tag, model, mode, regime, U, Kx, if (model == "nl") sprintf(" nests=%d", n_nests) else "",
            if (ioo) " outside option" else "",
            wmode, clusters, if (n_sit > 0) sprintf(" | first %.0f situations", n_sit) else "",
            paste(threads, collapse = ","), paste(what, collapse = ","), loadavg()))

rss0 <- rss_now()
t_gen <- system.time(g <- gen_regime(regime, U, seed = seed, K_x = Kx, K_w = 1L))[["elapsed"]]
unit_of <- rep.int(seq_along(g$Ti), g$Ti)       # decision maker of each situation
d <- list(X = g$X, alt_idx = g$alt_idx, choice_idx = g$choice_idx, M = g$M)
set.seed(seed + 1)
d$weights <- if (wmode == "runif") runif(length(d$M), 0.5, 2) else rep(1, length(d$M))
if (n_sit > 0 && n_sit < length(d$M)) {
  keep_rows <- sum(as.numeric(d$M[seq_len(n_sit)]))
  d$X <- d$X[seq_len(keep_rows), , drop = FALSE]
  d$alt_idx <- d$alt_idx[seq_len(keep_rows)]
  d$choice_idx <- d$choice_idx[seq_len(n_sit)]
  d$M <- d$M[seq_len(n_sit)]
  d$weights <- d$weights[seq_len(n_sit)]
  unit_of <- unit_of[seq_len(n_sit)]
}
g$W <- NULL
J <- g$J
beta <- 0.5 * seq(-1, 1, length.out = Kx)
delta <- if (ioo) g$delta else g$delta[-1]
nest_idx <- NULL
if (model == "nl") {
  nest_idx <- as.integer((seq_len(J) - 1L) %% n_nests + 1L)
  K_l <- sum(tabulate(nest_idx) > 1L)
  theta <- c(beta, rep(0.7, K_l), delta)
  lambda_full <- ifelse(tabulate(nest_idx) > 1L, 0.7, 1)
} else {
  theta <- c(beta, delta)
}
G <- if (clusters == "unit") NA_integer_ else as.integer(clusters)
cl <- if (is.na(G)) unit_of else as.integer((unit_of - 1L) %% G + 1L)
n_cl <- length(unique(cl))
N <- length(d$M); rows <- length(d$alt_idx)
rss_data <- rss_now()
cat(sprintf("[%s] situations=%.0f rows=%.3g J=%d n_params=%d clusters=%d | M: mean %.1f max %d | X %.2f GB\n",
            regime, N, rows, J, length(theta), n_cl, mean(d$M), max(d$M), 8 * rows * Kx / 1e9))

kern <- function(name) get0(paste0(model, "_", name), envir = ns, inherits = FALSE)
w_args <- function(w) {
  a <- list(theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx)
  if (model == "nl") a$nest_idx <- nest_idx
  c(a, list(M = d$M, weights = w, use_asc = TRUE, include_outside_option = ioo))
}
nw_args <- function() {
  a <- w_args(NULL); a$weights <- NULL; a
}
p_args <- function(w = NULL) {
  a <- list(theta = theta, X = d$X, alt_idx = d$alt_idx, M = d$M)
  if (model == "nl") a$nest_idx <- nest_idx
  if (!is.null(w)) a$weights <- w
  c(a, list(use_asc = TRUE, include_outside_option = ioo))
}
# acc_bytes for the J x J kernels of builds that take it
acc_arg <- function(f) {
  if (nzchar(acc) && "acc_bytes" %in% names(formals(f))) list(acc_bytes = as.numeric(acc)) else list()
}
# The inversion's target: the shares at theta (attainable, unlike uniform
# shares on regional choice sets, which the engine's checks stop)
blp_target <- if ("blp" %in% what) as.numeric(do.call(kern("predict_shares"), p_args(d$weights)))
# (the step alone where the build has methods, so that builds before and
# after Phase B, whose default is Anderson acceleration, time the same pass)
blp_method <- if ("method" %in% names(formals(ns$blp_contraction))) {
  list(method = "simple")
} else {
  list()
}
call_blp <- function() {
  # (it stops at max_iter: the engine warns, builds before 0.3.0 printed)
  utils::capture.output(r <- suppressWarnings(if (model == "nl") {
    do.call(ns$nl_blp_contraction, c(list(rep(0, J), blp_target, d$X, beta, lambda_full,
                                          d$alt_idx, nest_idx, d$M, d$weights, ioo,
                                          tol = 0, max_iter = blp_iter), blp_method))
  } else {
    do.call(ns$blp_contraction, c(list(rep(0, J), blp_target, d$X, beta, d$alt_idx, d$M,
                                       d$weights, ioo, tol = 0, max_iter = blp_iter),
                                  blp_method))
  }))
  r
}
call_kernel <- function(k) {
  switch(k,
    grad = do.call(kern("loglik_gradient_parallel"), w_args(d$weights)),
    grad_opg = do.call(kern("loglik_gradient_parallel"), c(w_args(d$weights), list(opg_diag = TRUE))),
    hess = do.call(kern("loglik_hessian_parallel"), w_args(d$weights)),
    bhhh = do.call(kern("bhhh_parallel"), w_args(d$weights)),
    robust = do.call(kern("bhhh_parallel"), w_args(d$weights^2)),
    scores = do.call(kern("scores_parallel"), nw_args()),
    meat_cluster = {
      f <- kern("cluster_meat_parallel")
      if (is.null(f)) stop("no cluster-meat kernel in this build")
      do.call(f, c(w_args(d$weights), list(cluster = match(cl, unique(cl)))))
    },
    meatr_robust = { S <- do.call(kern("scores_parallel"), nw_args()); crossprod(d$weights * S) },
    meatr_cluster = {
      S <- do.call(kern("scores_parallel"), nw_args())
      crossprod(rowsum(d$weights * S, group = as.character(cl)))
    },
    elas = do.call(kern("elasticities_parallel"), c(w_args(d$weights), list(elast_var_idx = 1L),
                                                     acc_arg(kern("elasticities_parallel")))),
    dr = do.call(kern("diversion_ratios_parallel"), c(p_args(d$weights),
                                                     acc_arg(kern("diversion_ratios_parallel")))),
    blp = call_blp(),
    predict = do.call(kern("predict"), p_args()),
    shares = do.call(kern("predict_shares"), p_args(d$weights)),
    stop("unknown kernel ", k))
}
check_value <- function(x) {
  if (is.list(x) && !is.null(x$objective)) return(x$objective)
  if (is.numeric(x)) return(sum(x))
  sum(vapply(x, function(e) sum(as.numeric(e)), 0))
}

rows_out <- list()
add_row <- function(...) rows_out[[length(rows_out) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
info <- list(tag = tag, lib = lib, model = model, mode = mode, regime = regime, U = U, N = N,
             ioo = ioo, acc = acc, blp_iter = blp_iter,
             rows = rows, J = J, n_par = length(theta), nests = if (model == "nl") n_nests,
             Kx = Kx, seed = seed, weights = wmode, clusters = clusters, n_clusters = n_cl,
             situations = n_sit, t_gen = t_gen, rss0 = rss0, rss_data = rss_data,
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

res <- if (length(rows_out)) do.call(rbind, rows_out) else data.frame()
info$load_end <- loadavg()
info$finished <- format(Sys.time())
attr(res, "info") <- info
saveRDS(res, out)
cat(sprintf("  saved %s | load %s -> %s\n", out, info$load_start, info$load_end))
