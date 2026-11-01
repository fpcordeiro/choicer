# BLP engine plan, Phase F: run_mxlogit(concentrate_out_asc = "blp").
# pf_units: the share inversion on a panel's estimation draws
#   (mxl_blp_solve(Ti =)): each situation reads its decision maker's draw
#   block, from a cube of one slice per decision maker or formed on the fly
#   (with and without the Gamma cache), with the usual pass and the
#   exp-space pass (its exponentials kept for every situation), with the
#   outside option and without; one situation per decision maker; and the
#   same draws repeated per situation without Ti, which main computes too.
#   Identity leaves, outside the solver tier (golden_blp() recomputes a
#   residual through the prediction kernels, which read one block per
#   situation): a change to the engine's iteration moves them, and Phases D
#   and E re-pin them. The entries with Ti are class D against main, whose
#   solver takes no Ti (its records are errors).

run_pf_units <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  bsolve <- get0("mxl_blp_solve", envir = ns)
  ctl <- function(pass) ns$.blp_control(1e-12, 1000L, method = "anderson",
                                        settings = ns$.BLP_SETTINGS$anderson,
                                        pass = pass, compensate = "off")
  rec <- function(inp, target, eta, pass, Ti, gen = c(-1L, 1L, 0L), cache = 0) {
    args <- list(rep(0, nrow(inp$alt_mapping) - inp$include_outside_option),
                 target, inp$X, inp$W, c(0.8, -0.6), c(0, 0), c(log(0.9), log(0.6)),
                 inp$alt_idx, inp$M, inp$weights, eta, c(0L, 0L), FALSE, FALSE,
                 inp$include_outside_option, ctl(pass), gen[1L], gen[2L], gen[3L],
                 cache)
    if (!is.null(Ti)) args$Ti <- Ti
    p2_val(do.call(bsolve, args))
  }
  # targets: each design's observed shares, outside first
  target_of <- function(inp) {
    J <- nrow(inp$alt_mapping) - inp$include_outside_option
    first <- c(0, cumsum(as.numeric(inp$M)))[seq_along(inp$M)]
    ch <- inp$choice_idx
    code <- integer(length(ch))
    code[ch > 0L] <- inp$alt_idx[first[ch > 0L] + ch[ch > 0L]]
    s <- tabulate(code[ch > 0L], J) / length(ch)
    if (inp$include_outside_option) c(1 - sum(s), s) else s
  }
  for (outside in c(TRUE, FALSE)) {
    pnl <- p5b_prep(p5b_sim(301L + outside, N = 60L, T = 3L, outside = outside),
                    outside = outside, person_col = "pid")
    tg <- target_of(pnl)
    U <- length(pnl$Ti)
    eta_u <- get_halton_normals(12L, U, 2L)
    tag <- if (outside) "pf_units|outside" else "pf_units|base"
    out[[paste(tag, "store", "log", sep = "|")]] <- rec(pnl, tg, eta_u, "log", pnl$Ti)
    out[[paste(tag, "store", "exp", sep = "|")]] <- rec(pnl, tg, eta_u, "exp", pnl$Ti,
                                                         cache = Inf)
    empty <- array(0, dim = c(2L, 0L, 0L))
    for (cache in c(0, Inf)) {
      out[[paste(tag, "generate", paste0("cache=", cache), sep = "|")]] <-
        rec(pnl, tg, empty, "log", pnl$Ti, gen = c(7L, 1L, 12L), cache = cache)
    }
    out[[paste(tag, "generate", "exp", sep = "|")]] <-
      rec(pnl, tg, empty, "exp", pnl$Ti, gen = c(7L, 1L, 12L), cache = Inf)
    # the same draws repeated per situation, without Ti
    out[[paste(tag, "repeated", sep = "|")]] <-
      rec(pnl, tg, eta_u[, , rep(seq_len(U), pnl$Ti), drop = FALSE], "log", NULL)
  }
  x <- p5b_prep(p5b_sim(303L, N = 80L))
  eta_x <- get_halton_normals(12L, length(x$M), 2L)
  out[["pf_units|one_per_unit"]] <- rec(x, target_of(x), eta_x, "log",
                                        rep(1L, length(x$M)))
  out[["pf_units|cross"]] <- rec(x, target_of(x), eta_x, "log", NULL)
  out
}

# pf_jac: mxl_blp_jacobian_parallel(), the objective, gradient, predicted
#   counts and share Jacobian (with variance = TRUE also K) on the estimation
#   draws, on small designs: a cross-section with the outside option, the
#   base normalization, a panel (also from pass 1's residuals), generated
#   draws, non-uniform weights, split draw batches, acc_bytes = 0 (one thread
#   here, so the same single result as without it), an alternative-level W, a
#   log-normal correlated design with rc_mean. Class D against main (a new
#   kernel); pinned for the later steps.
run_pf_jac <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  kern <- get0("mxl_blp_jacobian_parallel", envir = ns)
  if (is.null(kern)) {
    out[["pf_jac|absent"]] <- p2_val(stop("no mxl_blp_jacobian_parallel()"))
    return(out)
  }
  rec <- function(inp, eta, theta, rc_dist = c(0L, 0L), corr = FALSE, mean = FALSE,
                  ...) {
    p2_val(kern(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx, inp$M,
                inp$weights, eta, rc_dist, corr, mean, TRUE,
                inp$include_outside_option, Ti = inp$Ti, ...))
  }
  theta_of <- function(inp, n_head = 6L) {
    J <- nrow(inp$alt_mapping) - 1L
    set.seed(191)
    c(0.5, -0.4, 0.2, -0.1, log(0.9), log(0.6), stats::rnorm(J, sd = 0.5))[
      c(seq_len(n_head), 6L + seq_len(J))]
  }
  x <- p5b_prep(p5b_sim(181L))
  eta_x <- get_halton_normals(15L, length(x$M), 2L)
  th <- theta_of(x)
  # (theta: beta (2), L (2): the diagonal design has two continuous blocks)
  th <- th[-(3:4)]
  out[["pf_jac|cross"]] <- rec(x, eta_x, th)
  out[["pf_jac|cross|variance"]] <- rec(x, eta_x, th, variance = TRUE)
  out[["pf_jac|cross|batches"]] <- rec(x, eta_x, th, draw_batch = 4L, variance = TRUE)
  out[["pf_jac|cross|shared"]] <- rec(x, eta_x, th, acc_bytes = 0, variance = TRUE)
  out[["pf_jac|cross|generate"]] <- rec(x, array(0, dim = c(2L, 0L, 0L)), th,
                                         gen_seed = 3L, gen_scramble = 1L,
                                         gen_S = 15L, variance = TRUE)
  xw <- x
  set.seed(192)
  xw$weights <- stats::runif(length(x$M), 0.5, 2)
  out[["pf_jac|cross|weights"]] <- rec(xw, eta_x, th, variance = TRUE)
  xa <- x
  set.seed(193)
  xa$W <- matrix(stats::rnorm(2L * (nrow(x$alt_mapping) - 1L)), ncol = 2L)
  out[["pf_jac|cross|alt_W"]] <- rec(xa, eta_x, th, variance = TRUE)
  # log-normal second coefficient, correlated, with rc_mean
  th_l <- c(th[1:2], 0.1, -0.2, log(0.7), 0.2, log(0.5), th[-(1:4)])
  out[["pf_jac|cross|lognormal"]] <- rec(x, eta_x, th_l, rc_dist = c(0L, 1L),
                                          corr = TRUE, mean = TRUE, variance = TRUE)
  b <- p5b_prep(p5b_sim(182L, outside = FALSE), outside = FALSE)
  th_b <- theta_of(b)[-(3:4)]
  out[["pf_jac|base"]] <- rec(b, get_halton_normals(15L, length(b$M), 2L), th_b,
                               variance = TRUE)
  pnl <- p5b_prep(p5b_sim(183L, N = 80L, T = 3L), person_col = "pid")
  eta_p <- get_halton_normals(15L, length(pnl$Ti), 2L)
  out[["pf_jac|panel"]] <- rec(pnl, eta_p, theta_of(pnl)[-(3:4)], variance = TRUE)
  # the probabilities from pass 1's residuals, decision makers of several
  # situations each
  out[["pf_jac|panel|residuals"]] <- rec(pnl, eta_p, theta_of(pnl)[-(3:4)])
  out
}

# pf_meat: mxl_blp_meat_parallel(), the share-matching constants' meat (the
#   units' projected scores and share residuals; robust, by cluster, and for
#   supplied targets with the per-market sums, also by cluster), with D from
#   the share Jacobian at theta, on a cross-section with the outside option,
#   the base normalization, a panel, log-normal correlated coefficients with
#   rc_mean, an alternative-level W, in draw batches and on draws formed on
#   the fly. Class D against the steps before it (a new kernel); pinned for
#   the later steps.
run_pf_meat <- function() {
  out <- list()
  ns <- asNamespace("choicer")
  kern <- get0("mxl_blp_meat_parallel", envir = ns)
  if (is.null(kern)) {
    out[["pf_meat|absent"]] <- p2_val(stop("no mxl_blp_meat_parallel()"))
    return(out)
  }
  jac <- get("mxl_blp_jacobian_parallel", envir = ns)
  theta_of <- function(inp, n_head = 4L) {
    J <- nrow(inp$alt_mapping) - 1L
    set.seed(194)
    c(c(0.5, -0.4, log(0.9), log(0.6), 0.2, -0.1, log(0.7))[seq_len(n_head)],
      stats::rnorm(J, sd = 0.5))
  }
  # gen: the generator's arguments, for both kernels; spec: rc_dist, corr,
  # mean and the continuous parameters' count
  rec <- function(inp, eta, supplied = FALSE, weights = inp$weights, gen = list(),
                  spec = list(rc_dist = c(0L, 0L), corr = FALSE, mean = FALSE, K_c = 4L),
                  th = theta_of(inp, spec$K_c), ...) {
    p2_val({
      j <- do.call(jac, c(list(th, inp$X, inp$W, inp$alt_idx, inp$choice_idx, inp$M,
                               inp$weights, eta, spec$rc_dist, spec$corr, spec$mean,
                               TRUE, inp$include_outside_option), gen,
                          list(Ti = inp$Ti)))
      di <- spec$K_c + seq_len(length(th) - spec$K_c)
      D <- -solve(j$jacobian[di, di], j$jacobian[di, seq_len(spec$K_c)])
      do.call(kern, c(list(th, inp$X, inp$W, inp$alt_idx, inp$choice_idx, inp$M,
                           weights, eta, spec$rc_dist, D, supplied, spec$corr,
                           spec$mean, TRUE, inp$include_outside_option), gen,
                      list(Ti = inp$Ti, ...)))
    })
  }
  x <- p5b_prep(p5b_sim(184L))
  eta_x <- get_halton_normals(15L, length(x$M), 2L)
  N <- length(x$M)
  out[["pf_meat|cross|robust"]] <- rec(x, eta_x, weights = x$weights^2)
  out[["pf_meat|cross|cluster"]] <- rec(x, eta_x, cluster = rep(1:7, length.out = N))
  out[["pf_meat|cross|supplied"]] <- rec(x, eta_x, TRUE, alpha = x$weights,
                                         market = rep(1L, N))
  out[["pf_meat|cross|batches"]] <- rec(x, eta_x, weights = x$weights^2,
                                        draw_batch = 4L)
  # log-normal second coefficient, correlated, with rc_mean
  ln <- list(rc_dist = c(0L, 1L), corr = TRUE, mean = TRUE, K_c = 7L)
  out[["pf_meat|cross|lognormal"]] <- rec(x, eta_x, spec = ln, weights = x$weights^2)
  xa <- x
  set.seed(195)
  xa$W <- matrix(stats::rnorm(2L * (nrow(x$alt_mapping) - 1L)), ncol = 2L)
  out[["pf_meat|cross|alt_W"]] <- rec(xa, eta_x, weights = x$weights^2)
  b <- p5b_prep(p5b_sim(182L, outside = FALSE), outside = FALSE)
  out[["pf_meat|base"]] <- rec(b, get_halton_normals(15L, length(b$M), 2L),
                               weights = b$weights^2)
  pnl <- p5b_prep(p5b_sim(185L, N = 80L, T = 3L), person_col = "pid")
  eta_p <- get_halton_normals(15L, length(pnl$Ti), 2L)
  U <- length(pnl$Ti)
  w_u <- pnl$weights[c(0L, cumsum(pnl$Ti))[seq_len(U)] + 1L]
  cl_p <- rep(1:5, length.out = U)
  out[["pf_meat|panel|robust"]] <- rec(pnl, eta_p, weights = pnl$weights^2)
  out[["pf_meat|panel|cluster"]] <- rec(pnl, eta_p, cluster = cl_p)
  out[["pf_meat|panel|supplied"]] <- rec(pnl, eta_p, TRUE, alpha = w_u^2 * pnl$Ti,
                                         market = rep(1L, U))
  # by cluster with the per-market sums: kappa_g w_u
  kap <- tapply(w_u * pnl$Ti, cl_p, sum)[cl_p] * w_u
  out[["pf_meat|panel|supplied_cluster"]] <- rec(pnl, eta_p, TRUE, cluster = cl_p,
                                                 alpha = as.numeric(kap),
                                                 market = rep(1L, U))
  out[["pf_meat|panel|batches"]] <- rec(pnl, eta_p, TRUE, alpha = w_u^2 * pnl$Ti,
                                        market = rep(1L, U), draw_batch = 4L)
  out[["pf_meat|panel|generate"]] <- rec(pnl, array(0, dim = c(2L, 0L, 0L)),
                                         gen = list(gen_seed = 3L, gen_scramble = 1L,
                                                    gen_S = 15L))
  out
}

# pf_blpfit: run_mxlogit(concentrate_out_asc = "blp") on p5b_sim()'s data:
#   the sample's shares (outside option; the base normalization; a panel with
#   generated draws, weights and the sandwich), supplied targets within the
#   market (unnamed, named, with their covariance), se_method bhhh and
#   cluster. The inner settings are pinned here; the records move with the
#   share inversion's iteration, so a phase that changes it (BLP plan, Phases
#   D and E) re-pins them. Class D against the steps before it (on a build
#   without the method the records are errors).
run_pf_blpfit <- function() {
  out <- list()
  fit <- function(dt, ...) {
    p5a_rec(run_mxlogit(data = dt, id_col = "id", alt_col = "alt",
                        choice_col = "choice", covariate_cols = c("x1", "x2"),
                        random_var_cols = c("w1", "w2"), S = 15L,
                        concentrate_out_asc = "blp",
                        concentrate_control = list(max_iter = 1000L,
                                                   tol_shares = 1e-12), ...))
  }
  dt_x <- p5b_sim(186L)
  dt_b <- p5b_sim(187L, outside = FALSE)
  dt_p <- p5b_sim(188L, N = 80L, T = 3L)
  set.seed(189)
  wt <- stats::runif(80L, 0.5, 2)
  dt_p[, w := wt[pid]]
  dt_x[, cl := id %% 10L]
  # supplied targets: the sample's within-market shares pulled toward their
  # mean
  n_x <- data.table::uniqueN(dt_x$id)
  s <- tabulate(dt_x[choice == 1L]$alt, 8L) / n_x
  s_t <- 0.9 * s + 0.1 * mean(s)
  o <- list(include_outside_option = TRUE, outside_opt_label = 0L)
  cases <- list(
    cross = o,
    base = list(),
    panel_gen_w = c(o, list(person_col = "pid", draws = "generate", seed = 5L,
                            weights_col = "w", se_method = "sandwich")),
    supplied = c(o, list(target_shares = s_t)),
    supplied_named = c(o, list(target_shares = stats::setNames(rev(s_t),
                                                               as.character(8:1)))),
    supplied_vcov = c(o, list(target_shares = s_t,
                              target_vcov = (diag(s_t) - s_t %o% s_t) / 2000)),
    bhhh_se = c(o, list(se_method = "bhhh")),
    cluster = c(o, list(cluster_col = "cl"))
  )
  data_of <- list(base = dt_b, panel_gen_w = dt_p)
  for (nm in names(cases)) {
    dt <- data_of[[nm]] %||% dt_x
    out[[paste("pf_blpfit", nm, sep = "|")]] <- do.call(fit, c(list(dt), cases[[nm]]))
  }
  out
}
