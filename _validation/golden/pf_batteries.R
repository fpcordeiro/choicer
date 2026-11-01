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
