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
