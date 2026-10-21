# Phase 5B batteries.
# pdelta: mxl_delta_hessian_parallel(), the objective, gradient and the
#   constants' block of the information in one pass, on small designs: a
#   cross-section with the outside option, the base normalization, a panel,
#   generated draws, non-uniform weights, split draw batches, the shared
#   accumulator, an alternative-level W. Class D against main (a new kernel);
#   pinned for later phases.
# profile: run_mxlogit(concentrate_out_asc = "profile_likelihood") on the
#   same kind of data (outside option, base, panel with generated draws and
#   weights, scale_vars none/sd, optim, a bound on a continuous parameter, a
#   Newton polish after it, se_method bhhh/sandwich/cluster), plus "none".
#   Class D against main; on a build without the argument the records are
#   errors.

p5b_sim <- function(seed, N = 200L, J = 8L, T = 1L, outside = TRUE) {
  d <- data.table::as.data.table(simulate_mxl_data(
    N = N, J = J, T = T, beta = c(0.8, -0.6),
    delta = seq(-1, 1, length.out = J), Sigma = diag(c(1, 0.5)),
    outside_option = outside, seed = seed)$data)
  if (outside) d <- d[alt != 0L]
  d[]
}

p5b_prep <- function(dt, outside = TRUE, person_col = NULL) {
  prepare_mxl_data(dt, id_col = "id", alt_col = "alt", choice_col = "choice",
                   covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
                   include_outside_option = outside,
                   outside_opt_label = if (outside) 0L, person_col = person_col)
}

run_pdelta <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  kern <- get0("mxl_delta_hessian_parallel", envir = ns)
  call <- function(inp, eta, theta, ...) {
    p2_val(kern(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx, inp$M,
                inp$weights, eta, c(0L, 0L), FALSE, TRUE, TRUE,
                inp$include_outside_option, Ti = inp$Ti, ...))
  }
  if (is.null(kern)) {
    out[["pdelta|absent"]] <- p2_val(stop("no mxl_delta_hessian_parallel()"))
    return(out)
  }
  theta_of <- function(inp) {
    J <- nrow(inp$alt_mapping) - 1L
    set.seed(91)
    c(0.5, -0.4, 0.2, -0.1, log(0.9), log(0.6), stats::rnorm(J, sd = 0.5))
  }
  x <- p5b_prep(p5b_sim(81L))
  eta_x <- get_halton_normals(15L, length(x$M), 2L)
  th <- theta_of(x)
  out[["pdelta|cross"]] <- call(x, eta_x, th)
  out[["pdelta|cross|batches"]] <- call(x, eta_x, th, draw_batch = 4L)
  out[["pdelta|cross|shared"]] <- call(x, eta_x, th, acc_bytes = 0)
  out[["pdelta|cross|generate"]] <- call(x, array(0, dim = c(2L, 0L, 0L)), th,
                                         gen_seed = 3L, gen_scramble = 1L,
                                         gen_S = 15L)
  xw <- x
  set.seed(92)
  xw$weights <- stats::runif(length(x$M), 0.5, 2)
  out[["pdelta|cross|weights"]] <- call(xw, eta_x, th)
  xa <- x
  set.seed(93)
  xa$W <- matrix(stats::rnorm(2L * (nrow(x$alt_mapping) - 1L)), ncol = 2L)
  out[["pdelta|cross|alt_W"]] <- call(xa, eta_x, th)
  b <- p5b_prep(p5b_sim(82L, outside = FALSE), outside = FALSE)
  th_b <- theta_of(b)
  out[["pdelta|base"]] <- call(b, get_halton_normals(15L, length(b$M), 2L), th_b)
  pnl <- p5b_prep(p5b_sim(83L, N = 80L, T = 3L), person_col = "pid")
  out[["pdelta|panel"]] <- call(pnl, get_halton_normals(15L, length(pnl$Ti), 2L),
                                theta_of(pnl))
  out
}

run_profile <- function() {
  out <- list()
  fit <- function(dt, ...) {
    p5a_rec(run_mxlogit(data = dt, id_col = "id", alt_col = "alt",
                        choice_col = "choice", covariate_cols = c("x1", "x2"),
                        random_var_cols = c("w1", "w2"), S = 15L, ...))
  }
  dt_x <- p5b_sim(84L)
  dt_b <- p5b_sim(85L, outside = FALSE)
  dt_p <- p5b_sim(86L, N = 80L, T = 3L)
  set.seed(87)
  wt <- stats::runif(80L, 0.5, 2)
  dt_p[, w := wt[pid]]
  dt_x[, cl := id %% 10L]
  o <- list(include_outside_option = TRUE, outside_opt_label = 0L)
  P <- list(concentrate_out_asc = "profile_likelihood")
  cases <- list(
    cross = c(o, P),
    base = P,
    panel_gen_w = c(o, P, list(person_col = "pid", draws = "generate", seed = 5L,
                               weights_col = "w", se_method = "sandwich")),
    none = c(o, P, list(scale_vars = "none")),
    sd = c(o, P, list(scale_vars = "sd")),
    optim = c(o, P, list(optimizer = "optim")),
    bound = c(o, P, list(upper = c(x1 = 0.5))),
    newton = c(o, P, list(polish = "newton")),
    bhhh_se = c(o, P, list(se_method = "bhhh")),
    cluster = c(o, P, list(cluster_col = "cl")),
    off = c(o, list(concentrate_out_asc = "none"))
  )
  data_of <- list(base = dt_b, panel_gen_w = dt_p)
  for (nm in names(cases)) {
    dt <- data_of[[nm]] %||% dt_x
    out[[paste("profile", nm, sep = "|")]] <- do.call(fit, c(list(dt), cases[[nm]]))
  }
  out
}
