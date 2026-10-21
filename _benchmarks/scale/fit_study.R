#!/usr/bin/env Rscript
# Fits of run_mxlogit() on synthetic population-scale data: the default fit,
# the Newton phase (polish) and the profile likelihood (concentrate_out_asc),
# with every compiled kernel's calls counted, so that fits can be compared by
# their work when wall times are unreliable (a loaded machine).
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/fit_study.R fit regime=census U=20000 cell=lbfgs \
#     [threads=11] [out=<dir>]
#   R_LIBS=<lib> Rscript _benchmarks/scale/fit_study.R costs regime=census U=20000 \
#     ref=<dir>/census_20000_reference.rds [threads=11] [rounds=3] [out=<file.rds>]
#   Rscript _benchmarks/scale/fit_study.R table dir=<dir> [costs=<dir>] [csv=<file.csv>]
#
# Data: gen_regime() made D7-clean (every constant has a finite estimate):
# alternatives never chosen are dropped, and so are the choice situations that
# offer an alternative chosen wherever it is offered, until neither remains;
# the most chosen alternative is the base; no outside option; store-mode draws,
# S = 50 (run_mxlogit()'s advanced pathway takes a cube); rc_mean = TRUE.
# Cells (fit):
#   reference         the default fit, then Newton to a decrement below 1e-10
#   lbfgs             the default fit
#   lbfgs_newton      + polish = "newton"
#   lbfgs_bhhh        + polish = "bhhh"
#   k<K>_<newton|bhhh>  L-BFGS stopped after K evaluations, then the polish
#   profile_lbfgs     concentrate_out_asc = "profile_likelihood"
#   profile_newton    the same + polish = "newton"
# Each cell saves <dir>/<regime>_<U>_<cell>.rds: times, the kernels' calls,
# the Newton and BHHH decrements at the endpoint and the gradient there (from
# the kernels, not counted), and the estimates and standard errors. costs
# times each kernel pass at the reference estimates (medians over rounds);
# table compares each cell with its design's reference and, with costs=,
# models each cell's time as its calls times those medians.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
cmd <- commandArgs(TRUE)[1]
opt <- scale_args(commandArgs(TRUE)[-1])

# Keep rows `keep_row` of a stacked design and rebuild M, choice_idx and Ti
# (the chosen row of every kept situation must be kept).
subset_rows <- function(d, keep_row) {
  sit <- rep(seq_along(d$M), d$M)
  off <- c(0, cumsum(as.numeric(d$M)))[seq_along(d$M)]
  chosen_row <- off + d$choice_idx
  M2 <- tabulate(sit[keep_row], length(d$M))
  keep_sit <- M2 > 0L
  stopifnot(all(keep_row[chosen_row[keep_sit]]))
  cs <- cumsum(keep_row)
  off_kept <- c(0, cumsum(as.numeric(M2)))[seq_along(M2)]
  choice2 <- as.integer(cs[chosen_row] - off_kept)
  Ti2 <- tabulate(rep(seq_along(d$Ti), d$Ti)[keep_sit], length(d$Ti))
  d$X <- d$X[keep_row, , drop = FALSE]
  d$W <- d$W[keep_row, , drop = FALSE]
  d$alt_idx <- d$alt_idx[keep_row]
  d$M <- as.integer(M2[keep_sit])
  d$choice_idx <- choice2[keep_sit]
  d$Ti <- as.integer(Ti2[Ti2 > 0L])
  d
}

clean_d7 <- function(d) {
  repeat {
    off <- c(0, cumsum(as.numeric(d$M)))[seq_along(d$M)]
    J <- max(d$alt_idx)
    n_ch <- tabulate(d$alt_idx[off + d$choice_idx], J)
    n_obs <- tabulate(d$alt_idx, J)
    never <- n_obs > 0 & n_ch == 0
    always <- n_obs > 0 & n_ch == n_obs
    if (!any(never) && !any(always)) return(d)
    sit <- rep(seq_along(d$M), d$M)
    bad_sit <- unique(sit[always[d$alt_idx]])
    cat(sprintf("  dropping %d never-chosen alternatives and %d situations offering an always-chosen one\n",
                sum(never), length(bad_sit)))
    d <- subset_rows(d, !never[d$alt_idx] & !(sit %in% bad_sit))
  }
}

study_input <- function(regime, U, threads, S = 50L, n_clusters = 50L) {
  d <- gen_regime(regime, U, threads = threads)
  describe(d)
  d <- clean_d7(d)
  off <- c(0, cumsum(as.numeric(d$M)))[seq_along(d$M)]
  present <- sort(unique(d$alt_idx))
  n_chosen <- tabulate(match(d$alt_idx[off + d$choice_idx], present), length(present))
  order_alt <- present[order(-n_chosen, present)]
  code <- match(d$alt_idx, order_alt)
  J <- length(order_alt)
  am <- data.table::data.table(alt_int = seq_len(J), alt = order_alt, N_OBS = tabulate(code, J),
                               N_CHOICES = tabulate(code[off + d$choice_idx], J))
  am[, `:=`(TAKE_RATE = N_CHOICES / N_OBS, MKT_SHARE = N_CHOICES / sum(N_CHOICES))]
  unit <- rep(seq_along(d$Ti), d$Ti)
  inp <- structure(list(
    X = d$X, W = d$W, alt_idx = code, choice_idx = d$choice_idx, M = d$M, N = length(d$M),
    weights = rep(1, length(d$M)), cluster = (unit - 1L) %% n_clusters + 1L,
    situation_ids = seq_along(d$M), Ti = d$Ti, person_ids = seq_along(d$Ti),
    include_outside_option = FALSE, rc_correlation = TRUE, alt_mapping = am,
    data_spec = list(id_col = "id", alt_col = "alt", choice_col = "choice")),
    class = c("choicer_data_mxl", "list"))
  stopifnot(is.null(choicer:::.unbounded_constants(am)))
  cat(sprintf("[%s %d] J %d; situations %d; units %d; parameters %d\n", regime, U, J,
              length(d$M), length(d$Ti), ncol(d$X) + ncol(d$W) * (ncol(d$W) + 3L) / 2L + J - 1L))
  list(inp = inp, eta = get_halton_normals(S, length(d$Ti), ncol(d$W)))
}

kernel_args <- function(x) {
  list(X = x$inp$X, W = x$inp$W, alt_idx = x$inp$alt_idx, choice_idx = x$inp$choice_idx,
       M = x$inp$M, weights = x$inp$weights, eta_draws = x$eta, rc_dist = rep(0L, ncol(x$inp$W)),
       rc_correlation = TRUE, rc_mean = TRUE, use_asc = TRUE, include_outside_option = FALSE,
       Ti = x$inp$Ti)
}

if (identical(cmd, "fit")) {
  lib <- load_choicer()
  ns <- asNamespace("choicer")
  regime <- opt$regime; U <- as.integer(opt$U); cell <- opt$cell
  threads <- as.integer(opt$threads %||% "11")
  out_dir <- opt$out %||% scale_output("fits")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(out_dir, sprintf("%s_%d_%s.rds", regime, U, cell))
  if (file.exists(out)) { cat("exists:", out, "\n"); quit(save = "no") }
  set_num_threads(threads)
  x <- study_input(regime, U, threads)
  args <- switch(sub("^k[0-9]+_", "k_", cell),
    reference = list(polish = "newton", polish_control = list(tol_decrement = 1e-10, max_iter = 100L)),
    lbfgs = list(),
    lbfgs_newton = list(polish = "newton"),
    lbfgs_bhhh = list(polish = "bhhh"),
    k_newton = , k_bhhh = list(control = list(maxeval = as.integer(sub("^k([0-9]+)_.*$", "\\1", cell))),
                               polish = sub("^k[0-9]+_", "", cell),
                               polish_control = list(max_iter = if (grepl("newton", cell)) 50L else 100L)),
    profile_lbfgs = list(concentrate_out_asc = "profile_likelihood"),
    profile_newton = list(concentrate_out_asc = "profile_likelihood", polish = "newton"),
    stop("unknown cell ", cell))

  # Count the kernels' calls during the fit; newton_step() is the Cholesky
  calls <- c(gradient = 0L, hessian = 0L, bhhh = 0L, meat = 0L, delta = 0L, cholesky = 0L)
  kern <- c(mxl_loglik_gradient_parallel = "gradient", mxl_hessian_parallel = "hessian",
            mxl_bhhh_parallel = "bhhh", mxl_cluster_meat_parallel = "meat",
            mxl_delta_hessian_parallel = "delta", newton_step = "cholesky")
  kern <- kern[names(kern) %in% ls(ns, all.names = TRUE)]
  for (k in names(kern)) {
    suppressMessages(trace(k, tracer = bquote(calls[[.(kern[[k]])]] <<- calls[[.(kern[[k]])]] + 1L),
                           where = ns, print = FALSE))
  }
  load0 <- loadavg()
  t <- system.time(f <- suppressMessages(do.call(run_mxlogit, c(list(
    input_data = x$inp, eta_draws = x$eta, rc_mean = TRUE), args))), gcFirst = FALSE)
  for (k in names(kern)) suppressMessages(untrace(k, where = ns))

  # Certificates at the endpoint, from the kernels (not counted)
  theta <- unname(coef(f))
  ka <- c(list(theta = theta), kernel_args(x))
  gr <- do.call(ns$mxl_loglik_gradient_parallel, ka)
  dec <- function(C) {
    r <- ns$newton_step(C, gr$gradient, ns$.jacobi_scale(C), seq_along(theta), 0)
    if (r$info != 0L) NA_real_ else -0.5 * sum(gr$gradient * r$step)
  }
  dec_newton <- dec(do.call(ns$mxl_hessian_parallel, ka)); invisible(gc(FALSE))
  dec_bhhh <- dec(do.call(ns$mxl_bhhh_parallel, ka)); invisible(gc(FALSE))
  p <- f$optimizer$polish
  res <- list(
    regime = regime, U = U, cell = cell, threads = threads, n_params = length(theta),
    J = nrow(x$inp$alt_mapping), situations = length(x$inp$M), units = length(x$inp$Ti),
    wall = t[["elapsed"]], fit_time = f$optimizer$elapsed_time,
    optimizer_time = p$optimizer$elapsed_time %||% f$optimizer$elapsed_time,
    polish_time = p$elapsed_time %||% 0, iterations = f$optimizer$iterations,
    convergence = f$convergence, message = f$message,
    polish = if (!is.null(p)) p[c("converged", "iterations", "damped", "decrement", "lambda",
                                  "message", "evaluations", "held", "flat", "improvement", "error")],
    concentration = f$optimizer$concentration, calls = calls, loglik = f$loglik, coef = coef(f),
    se = f$se, objective_at_end = gr$objective, gradient_at_end = gr$gradient,
    decrement_newton = dec_newton, decrement_bhhh = dec_bhhh, load_start = load0,
    load_end = loadavg(), blas = extSoftVersion()[["BLAS"]], lib = lib,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  saveRDS(res, out)
  cat(sprintf("[%s %d %s] wall %.1f s (polish %.1f) | LL %.6f | iterations %s | Newton decrement %.3g, BHHH %.3g | calls %s\n",
              regime, U, cell, res$wall, res$polish_time, res$loglik, res$iterations, dec_newton,
              dec_bhhh, paste(names(calls), calls, sep = "=", collapse = " ")))
} else if (identical(cmd, "costs")) {
  lib <- load_choicer()
  ns <- asNamespace("choicer")
  regime <- opt$regime; U <- as.integer(opt$U)
  threads <- as.integer(opt$threads %||% "11"); rounds <- as.integer(opt$rounds %||% "3")
  ref <- readRDS(opt$ref)
  out <- opt$out %||% file.path(dirname(opt$ref), sprintf("costs_%s_%d.rds", regime, U))
  set_num_threads(threads)
  x <- study_input(regime, U, threads)
  theta <- unname(ref$coef)
  ka <- c(list(theta = theta), kernel_args(x))
  K_c <- ncol(x$inp$X) + ncol(x$inp$W) * (ncol(x$inp$W) + 3L) / 2L
  d_idx <- (K_c + 1L):length(theta)
  call_k <- function(k) do.call(ns[[k]], ka)
  H <- call_k("mxl_hessian_parallel"); B <- call_k("mxl_bhhh_parallel")
  g <- call_k("mxl_loglik_gradient_parallel")$gradient
  Hd <- H[d_idx, d_idx]; gd <- g[d_idx]
  D <- ns$.jacobi_scale(H); Dd <- ns$.jacobi_scale(Hd)
  tasks <- list(
    gradient = function() call_k("mxl_loglik_gradient_parallel"),
    hessian = function() call_k("mxl_hessian_parallel"),
    bhhh = function() call_k("mxl_bhhh_parallel"),
    delta = function() call_k("mxl_delta_hessian_parallel"),
    cholesky_p = function() ns$newton_step(H, g, D, seq_along(theta), 0),
    cholesky_d = function() ns$newton_step(Hd, gd, Dd, seq_along(d_idx), 0),
    inverse_p = function() ns$invert_hessian(H),
    sandwich_p = function() ns$.sandwich_combine(H, B))
  times <- matrix(NA_real_, rounds, length(tasks), dimnames = list(NULL, names(tasks)))
  for (r in seq_len(rounds)) for (k in names(tasks)) {
    gc(FALSE)
    times[r, k] <- system.time(suppressMessages(tasks[[k]]()), gcFirst = FALSE)[["elapsed"]]
  }
  res <- list(regime = regime, U = U, threads = threads, n_params = length(theta),
              J_delta = length(d_idx), times = times, median = apply(times, 2, stats::median),
              load = loadavg(), blas = extSoftVersion()[["BLAS"]])
  saveRDS(res, out)
  print(round(res$median, 4))
} else if (identical(cmd, "table")) {
  files <- list.files(opt$dir, pattern = "\\.rds$", full.names = TRUE)
  cells <- Filter(function(x) !is.null(x$cell), lapply(files, readRDS))
  refs <- Filter(function(x) x$cell == "reference", cells)
  names(refs) <- vapply(refs, function(x) paste(x$regime, x$U), "")
  costs <- list()
  if (!is.null(opt$costs)) for (f in list.files(opt$costs, pattern = "^costs_.*\\.rds$", full.names = TRUE)) {
    x <- readRDS(f); costs[[paste(x$regime, x$U)]] <- x$median
  }
  num <- function(x) if (is.null(x)) NA_real_ else x
  rows <- lapply(cells, function(x) {
    key <- paste(x$regime, x$U); r <- refs[[key]]; m <- costs[[key]]
    n <- x$calls; p <- x$polish
    # Cholesky factorizations: of the constants' block inside the profile's
    # solves, of the information at p in the polish
    prof <- grepl("^profile", x$cell)
    n_chol_p <- if (prof) num(p$evaluations[["curvature"]]) %||% 0 else n[["cholesky"]]
    if (is.na(n_chol_p)) n_chol_p <- 0
    model <- if (is.null(m)) NA_real_ else sum(
      n[["gradient"]] * m[["gradient"]], n[["hessian"]] * m[["hessian"]], n[["bhhh"]] * m[["bhhh"]],
      n[["delta"]] * m[["delta"]], n_chol_p * m[["cholesky_p"]],
      (n[["cholesky"]] - n_chol_p) * m[["cholesky_d"]], m[["inverse_p"]], na.rm = TRUE)
    data.frame(
      design = key, cell = x$cell, n_params = x$n_params, wall = x$wall, model = model,
      iterations = x$iterations, certified = if (is.null(p)) NA else isTRUE(p$converged),
      ll_gap = if (is.null(r)) NA_real_ else r$loglik - x$loglik,
      dec_newton = x$decrement_newton, dec_bhhh = x$decrement_bhhh,
      dcoef_se = if (is.null(r)) NA_real_ else max(abs(x$coef - r$coef) / r$se, na.rm = TRUE),
      gradient = n[["gradient"]], hessian = n[["hessian"]], bhhh = n[["bhhh"]], delta = n[["delta"]],
      cholesky = n[["cholesky"]],
      inner = if (is.null(x$concentration)) NA_integer_ else x$concentration$inner_iterations,
      load = paste(x$load_start, "->", x$load_end), stringsAsFactors = FALSE)
  })
  t <- do.call(rbind, rows)
  ord <- c("reference", "lbfgs", "lbfgs_newton", "lbfgs_bhhh",
           paste0("k", rep(c(25, 50, 100, 200), each = 2), c("_newton", "_bhhh")),
           "profile_lbfgs", "profile_newton")
  t <- t[order(t$design, match(t$cell, ord)), ]
  op <- options(width = 250)
  print(format(t, digits = 3), row.names = FALSE)
  options(op)
  if (!is.null(opt$csv)) utils::write.csv(t, opt$csv, row.names = FALSE)
} else {
  stop("usage: fit_study.R fit|costs|table key=value ... (see the header)")
}
