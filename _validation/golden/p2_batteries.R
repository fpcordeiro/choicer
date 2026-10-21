# ============================================================================
# Phase 2 batteries (scale audit: index vectors read in place, the NL nest
# helper, HB/MNP integer widths, ARMA_64BIT_WORD), sourced by
# refactor_identity_check.R after the earlier batteries. Every call exists
# with the same R signature on main and on the Phase 2 builds.
#   * mnl_x  : MNL bhhh / scores / shares / BLP over asc x outside option,
#              small (J = 5) and wide (J = 30, sets straddling 16) designs.
#   * nl_x   : NL analytic Hessian / bhhh / scores and the rest of the NL
#              family, including a 20-nest design with sets straddling 16.
#   * idx_dbl: every MNL/NL kernel with double-typed index vectors, and
#              whether the result equals the integer-typed call.
#   * err_p2 : every MNL/NL validation path, per kernel: condition class,
#              message and whether a call is attached, plus warnings.
#   * hb     : hmnl_gibbs / hmnp_gibbs / mnp_gibbs and run_hmnlogit /
#              run_hmnprobit / run_mnprobit at 1 thread, fixed seeds.
#   * fits   : small run_mnlogit / run_nestlogit fits and their
#              post-estimation, timings removed.
# Inputs that terminate R on main (nl_blp_contraction with a bad nest code or
# lambda length) are not here: p2_blp_abort.R runs them in subprocesses.
# ============================================================================

p2_rt <- function(x) unserialize(serialize(x, NULL))

p2_leaf <- function(x) {
  x <- p2_rt(x)
  if (is.list(x)) attr(x, "golden_leaf") <- TRUE
  x
}

# Condition of an expression: "<no error>" or "<class> | <message> |
# call=<TRUE/FALSE>", followed by any warnings, as one character vector.
p2_cond <- function(expr) {
  warns <- character()
  out <- tryCatch(
    withCallingHandlers({
      force(expr)
      "<no error>"
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      sprintf("%s | %s | call=%s", class(e)[1], conditionMessage(e),
              !is.null(conditionCall(e)))
    })
  c(out, warns)
}

# Value of an expression plus its warnings and messages (or the error).
p2_val <- function(expr) {
  warns <- character(); msgs <- character()
  err <- NA_character_
  val <- tryCatch(
    withCallingHandlers(expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      }),
    error = function(e) {
      err <<- sprintf("%s | %s", class(e)[1], conditionMessage(e))
      NULL
    })
  list(value = if (is.null(val)) NULL else p2_leaf(p2_strip(val)),
       warnings = warns, messages = msgs, error = err)
}

# Drop wall-clock fields anywhere in a fitted object.
p2_strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("elapsed_time", "elapsed", "time_elapsed")
    # Phase 6D: fields main's fit objects lack (the p6d batteries record
    # them), so whole-fit leaves compare the fields both builds have.
    if (inherits(x, "choicer_fit")) {
      drop <- drop | names(x) %in% c("param_scale", "param_shift")
      if (inherits(x, "choicer_nl")) drop <- drop | names(x) %in% c("scale_vars", "sX")
    }
    if (any(drop)) {
      a <- attributes(x)
      x <- x[!drop]
      keep <- setdiff(names(a), c("names", "row.names"))
      for (nm in keep) attr(x, nm) <- a[[nm]]
    }
    for (i in seq_along(x)) {
      if (!is.null(x[[i]])) x[[i]] <- p2_strip(x[[i]])
    }
  }
  x
}

p2_dat <- function(seed, N, J, K, ioo) {
  make_choice_data(seed = seed, N = N, J = J, K = K,
                   include_outside_option = ioo)
}

# ---------------------------------------------------------------------------
# MNL extras
# ---------------------------------------------------------------------------

p2_mnl_outputs <- function(theta, dat, use_asc, ioo, with_blp) {
  res <- list()
  res$grad <- mnl_loglik_gradient_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights, use_asc, ioo)
  res$bhhh <- choicer:::mnl_bhhh_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights, use_asc, ioo)
  res$bhhh_w2 <- choicer:::mnl_bhhh_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights^2, use_asc, ioo)
  res$scores <- choicer:::mnl_scores_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, use_asc, ioo)
  res$pred <- mnl_predict(theta, dat$X, dat$alt_idx, dat$M, use_asc, ioo)
  res$shares <- mnl_predict_shares(
    theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
  res$hess <- mnl_loglik_hessian_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights, use_asc, ioo)
  res$elas2 <- mnl_elasticities_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
    elast_var_idx = 2L, use_asc, ioo)
  res$dr <- mnl_diversion_ratios_parallel(
    theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
  if (with_blp) {
    beta <- theta[seq_len(dat$K)]
    res$blp <- blp_contraction(rep(0, dat$J), res$shares, dat$X, beta,
                               dat$alt_idx, dat$M, dat$weights, ioo)
    res$blp_short <- suppressWarnings(utils::capture.output(
      v <- blp_contraction(rep(0.2, dat$J), res$shares, dat$X, beta,
                           dat$alt_idx, dat$M, dat$weights, ioo,
                           tol = 1e-14, max_iter = 3L)))
    res$blp_short_val <- v
  }
  res
}

run_mnl_x <- function() {
  out <- list()
  designs <- list(small = list(N = 40, J = 5, K = 3), wide = list(N = 60, J = 30, K = 3))
  for (dn in names(designs)) {
    g <- designs[[dn]]
    for (ioo in c(TRUE, FALSE)) {
      dat <- p2_dat(seed = 4100 + 10 * ioo + (dn == "wide"), N = g$N, J = g$J,
                    K = g$K, ioo = ioo)
      for (use_asc in c(TRUE, FALSE)) {
        for (rep in 1:2) {
          set.seed(4200 + 10 * rep + 2 * ioo + use_asc)
          theta <- rnorm(dat$K + n_delta_free(dat$J, use_asc, ioo), sd = 0.4)
          key <- sprintf("mnl_x|%s|asc=%d|oo=%d|rep=%d", dn, use_asc, ioo, rep)
          out[[key]] <- p2_mnl_outputs(theta, dat, use_asc, ioo,
                                       with_blp = use_asc && rep == 1)
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# NL extras
# ---------------------------------------------------------------------------

p2_nl_theta <- function(K, nest_idx, J, use_asc, ioo, seed) {
  n_nonsing <- sum(table(nest_idx) > 1)
  set.seed(seed)
  theta <- rnorm(K + n_nonsing + n_delta_free(J, use_asc, ioo), sd = 0.4)
  theta[K + seq_len(n_nonsing)] <- runif(n_nonsing, 0.55, 1.0)
  theta
}

p2_lambda_full <- function(theta, K, nest_idx) {
  n_nests <- max(nest_idx)
  counts <- tabulate(nest_idx, n_nests)
  lam <- rep(1, n_nests)
  ns <- which(counts > 1)
  lam[ns] <- theta[K + seq_along(ns)]
  lam
}

p2_nl_outputs <- function(theta, dat, nest_idx, use_asc, ioo, with_blp) {
  res <- list()
  res$grad <- nl_loglik_gradient_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M, dat$weights,
    use_asc, ioo)
  res$hess <- choicer:::nl_loglik_hessian_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M, dat$weights,
    use_asc, ioo)
  res$bhhh <- choicer:::nl_bhhh_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M, dat$weights,
    use_asc, ioo)
  res$scores <- choicer:::nl_scores_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M, use_asc, ioo)
  res$pred <- nl_predict(theta, dat$X, dat$alt_idx, dat$M, nest_idx, use_asc, ioo)
  res$shares <- nl_predict_shares(theta, dat$X, dat$alt_idx, dat$M, dat$weights,
                                  nest_idx, use_asc, ioo)
  res$elas2 <- nl_elasticities_parallel(
    theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M, dat$weights,
    elast_var_idx = 2L, use_asc, ioo)
  res$dr <- nl_diversion_ratios_parallel(
    theta, dat$X, dat$alt_idx, nest_idx, dat$M, dat$weights, use_asc, ioo)
  if (with_blp) {
    beta <- theta[seq_len(dat$K)]
    lam <- p2_lambda_full(theta, dat$K, nest_idx)
    res$blp <- nl_blp_contraction(rep(0, dat$J), res$shares, dat$X, beta, lam,
                                  dat$alt_idx, nest_idx, dat$M, dat$weights, ioo)
    res$blp_damped <- nl_blp_contraction(rep(0, dat$J), res$shares, dat$X, beta,
                                         lam, dat$alt_idx, nest_idx, dat$M,
                                         dat$weights, ioo, damping = 0.6)
  }
  res
}

run_nl_x <- function() {
  out <- list()
  designs <- list(
    small = list(N = 40, J = 5, K = 3, nest = c(1L, 1L, 2L, 2L, 3L)),
    # 20 nests over 40 alternatives: 15 pairs, 2 triples... and 3 singletons
    wide = list(N = 60, J = 40, K = 3,
                nest = as.integer(c(rep(1:15, each = 2), 16L, 16L, 16L,
                                    17L, 17L, 17L, 17L, 18L, 19L, 20L))))
  for (dn in names(designs)) {
    g <- designs[[dn]]
    stopifnot(length(g$nest) == g$J)
    for (ioo in c(TRUE, FALSE)) {
      dat <- p2_dat(seed = 5100 + 10 * ioo + (dn == "wide"), N = g$N, J = g$J,
                    K = g$K, ioo = ioo)
      for (use_asc in c(TRUE, FALSE)) {
        for (rep in 1:2) {
          theta <- p2_nl_theta(g$K, g$nest, g$J, use_asc, ioo,
                               seed = 5200 + 10 * rep + 2 * ioo + use_asc)
          key <- sprintf("nl_x|%s|asc=%d|oo=%d|rep=%d", dn, use_asc, ioo, rep)
          out[[key]] <- p2_nl_outputs(theta, dat, g$nest, use_asc, ioo,
                                      with_blp = use_asc && rep == 1)
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Double-typed index vectors
# ---------------------------------------------------------------------------

p2_as_dbl <- function(dat) {
  dat$alt_idx <- as.double(dat$alt_idx)
  dat$choice_idx <- as.double(dat$choice_idx)
  dat$M <- as.double(dat$M)
  dat
}

run_idx_dbl <- function() {
  out <- list()
  for (ioo in c(TRUE, FALSE)) {
    dat <- p2_dat(seed = 6100 + ioo, N = 50, J = 30, K = 3, ioo = ioo)
    for (use_asc in c(TRUE, FALSE)) {
      set.seed(6200 + 2 * ioo + use_asc)
      theta <- rnorm(dat$K + n_delta_free(dat$J, use_asc, ioo), sd = 0.4)
      a <- p2_mnl_outputs(theta, dat, use_asc, ioo, with_blp = use_asc)
      b <- p2_mnl_outputs(theta, p2_as_dbl(dat), use_asc, ioo, with_blp = use_asc)
      key <- sprintf("idx_dbl|mnl|asc=%d|oo=%d", use_asc, ioo)
      out[[key]] <- list(dbl = b, same_as_int = identical(a, b))
    }
    nest <- as.integer(c(rep(1:12, each = 2), 13L, 13L, 13L, 14L, 15L, 16L))
    dat <- p2_dat(seed = 6300 + ioo, N = 50, J = 30, K = 3, ioo = ioo)
    for (use_asc in c(TRUE, FALSE)) {
      theta <- p2_nl_theta(3, nest, 30, use_asc, ioo, seed = 6400 + 2 * ioo + use_asc)
      a <- p2_nl_outputs(theta, dat, nest, use_asc, ioo, with_blp = use_asc)
      b <- p2_nl_outputs(theta, p2_as_dbl(dat), as.double(nest), use_asc, ioo,
                         with_blp = use_asc)
      key <- sprintf("idx_dbl|nl|asc=%d|oo=%d", use_asc, ioo)
      out[[key]] <- list(dbl = b, same_as_int = identical(a, b))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Error and message battery
# ---------------------------------------------------------------------------

# Call one MNL kernel by name with the inputs in `f` (theta, X, alt_idx,
# choice_idx, M, weights, use_asc, ioo, plus extras for BLP/elasticities).
p2_mnl_call <- function(k, f) {
  switch(k,
    grad = mnl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                        f$M, f$weights, f$use_asc, f$ioo),
    bhhh = choicer:::mnl_bhhh_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                       f$M, f$weights, f$use_asc, f$ioo),
    scores = choicer:::mnl_scores_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                           f$M, f$use_asc, f$ioo),
    pred = mnl_predict(f$theta, f$X, f$alt_idx, f$M, f$use_asc, f$ioo),
    shares = mnl_predict_shares(f$theta, f$X, f$alt_idx, f$M, f$weights,
                                f$use_asc, f$ioo),
    hess = mnl_loglik_hessian_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                       f$M, f$weights, f$use_asc, f$ioo),
    elas = mnl_elasticities_parallel(f$theta, f$X, f$alt_idx, f$choice_idx, f$M,
                                     f$weights, f$var_idx %||% 1L, f$use_asc,
                                     f$ioo),
    dr = mnl_diversion_ratios_parallel(f$theta, f$X, f$alt_idx, f$M, f$weights,
                                       f$use_asc, f$ioo),
    blp = blp_contraction(f$delta %||% rep(0, f$J), f$target %||% rep(1 / (f$J + f$ioo), f$J + f$ioo),
                          f$X, f$theta[seq_len(ncol(f$X))], f$alt_idx, f$M,
                          f$weights, f$ioo, max_iter = 5L))
}

p2_nl_call <- function(k, f) {
  switch(k,
    grad = nl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                       f$nest_idx, f$M, f$weights, f$use_asc,
                                       f$ioo),
    hess = choicer:::nl_loglik_hessian_parallel(f$theta, f$X, f$alt_idx,
                                                f$choice_idx, f$nest_idx, f$M,
                                                f$weights, f$use_asc, f$ioo),
    num_hess = nl_loglik_numeric_hessian(f$theta, f$X, f$alt_idx, f$choice_idx,
                                         f$nest_idx, f$M, f$weights, f$use_asc,
                                         f$ioo),
    bhhh = choicer:::nl_bhhh_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                      f$nest_idx, f$M, f$weights, f$use_asc,
                                      f$ioo),
    scores = choicer:::nl_scores_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                          f$nest_idx, f$M, f$use_asc, f$ioo),
    pred = nl_predict(f$theta, f$X, f$alt_idx, f$M, f$nest_idx, f$use_asc, f$ioo),
    shares = nl_predict_shares(f$theta, f$X, f$alt_idx, f$M, f$weights,
                               f$nest_idx, f$use_asc, f$ioo),
    elas = nl_elasticities_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                    f$nest_idx, f$M, f$weights, f$var_idx %||% 1L,
                                    f$use_asc, f$ioo),
    dr = nl_diversion_ratios_parallel(f$theta, f$X, f$alt_idx, f$nest_idx, f$M,
                                      f$weights, f$use_asc, f$ioo),
    blp = nl_blp_contraction(f$delta %||% rep(0, f$J),
                             f$target %||% rep(1 / (f$J + f$ioo), f$J + f$ioo),
                             f$X, f$theta[seq_len(ncol(f$X))],
                             f$lambda %||% c(0.7, 0.8, 1), f$alt_idx, f$nest_idx,
                             f$M, f$weights, f$ioo, max_iter = 5L))
}

P2_MNL_K <- c("grad", "bhhh", "scores", "pred", "shares", "hess", "elas", "dr", "blp")
P2_NL_K <- c("grad", "hess", "num_hess", "bhhh", "scores", "pred", "shares",
             "elas", "dr", "blp")

p2_with <- function(f, ...) {
  a <- list(...)
  for (nm in names(a)) f[nm] <- list(a[[nm]])
  f
}

run_err_p2 <- function() {
  out <- list()
  J <- 4L; K <- 2L
  base <- function(ioo, use_asc = TRUE) {
    d <- p2_dat(seed = 7100 + ioo, N = 10, J = J, K = K, ioo = ioo)
    set.seed(7200 + ioo + 2 * use_asc)
    th <- rnorm(K + n_delta_free(J, use_asc, ioo), sd = 0.3)
    list(theta = th, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
         M = d$M, weights = d$weights, use_asc = use_asc, ioo = ioo, J = J)
  }
  f0 <- base(FALSE); fo <- base(TRUE); fn <- base(FALSE, FALSE)
  n <- sum(f0$M)
  alt_rep <- function(f, i, v) { a <- f$alt_idx; a[i] <- v; a }
  ch_rep <- function(f, i, v) { a <- f$choice_idx; a[i] <- v; a }
  empty <- function(f) p2_with(f, X = f$X[0, , drop = FALSE], alt_idx = integer(0),
                               choice_idx = integer(0), M = integer(0),
                               weights = numeric(0))
  cases <- list(
    ok = f0, ok_oo = fo, ok_noasc = fn,
    M_zero = p2_with(f0, M = replace(f0$M, 3L, 0L)),
    M_na = p2_with(f0, M = replace(f0$M, 3L, NA_integer_)),
    M_neg = p2_with(f0, M = replace(f0$M, 3L, -2L)),
    X_short = p2_with(f0, X = f0$X[-1L, , drop = FALSE]),
    M_long = p2_with(f0, M = c(f0$M, 2L)),
    alt_short = p2_with(f0, alt_idx = f0$alt_idx[-1L]),
    w_short = p2_with(f0, weights = f0$weights[-1L]),
    ch_short = p2_with(f0, choice_idx = f0$choice_idx[-1L]),
    alt_zero = p2_with(f0, alt_idx = alt_rep(f0, 5L, 0L)),
    alt_na = p2_with(f0, alt_idx = alt_rep(f0, 5L, NA_integer_)),
    alt_neg = p2_with(f0, alt_idx = alt_rep(f0, 5L, -1L)),
    alt_3e9 = p2_with(f0, alt_idx = replace(as.double(f0$alt_idx), 5L, 3e9)),
    alt_3e9_noasc = p2_with(fn, alt_idx = replace(as.double(fn$alt_idx), 5L, 3e9)),
    alt_big = p2_with(f0, alt_idx = alt_rep(f0, 5L, J + 1L)),
    alt_big_noasc = p2_with(fn, alt_idx = alt_rep(fn, 5L, J + 1L)),
    alt_frac = p2_with(f0, alt_idx = as.double(f0$alt_idx) + 0.25),
    alt_lgl = p2_with(f0, alt_idx = rep(TRUE, n)),
    alt_factor = p2_with(f0, alt_idx = factor(f0$alt_idx)),
    alt_null = p2_with(f0, alt_idx = NULL),
    alt_chr = p2_with(f0, alt_idx = as.character(f0$alt_idx)),
    ch_past = p2_with(f0, choice_idx = ch_rep(f0, 4L, f0$M[4L] + 1L)),
    ch_zero = p2_with(f0, choice_idx = ch_rep(f0, 4L, 0L)),
    ch_na = p2_with(f0, choice_idx = ch_rep(f0, 4L, NA_integer_)),
    ch_neg = p2_with(f0, choice_idx = ch_rep(f0, 4L, -1L)),
    ch_past_oo = p2_with(fo, choice_idx = ch_rep(fo, 4L, fo$M[4L] + 1L)),
    ch_na_oo = p2_with(fo, choice_idx = ch_rep(fo, 4L, NA_integer_)),
    ch_neg_oo = p2_with(fo, choice_idx = ch_rep(fo, 4L, -1L)),
    ch_null = p2_with(f0, choice_idx = NULL),
    theta_short = p2_with(f0, theta = f0$theta[1:2]),
    theta_long_noasc = p2_with(fn, theta = c(fn$theta, 0.1)),
    delta_short = p2_with(f0, theta = f0$theta[-length(f0$theta)]),
    empty_asc = empty(f0), empty_noasc = empty(fn),
    empty_oo = empty(fo),
    w_zero = p2_with(f0, weights = rep(0, length(f0$weights))),
    var_oob = p2_with(f0, var_idx = 3L),
    blp_target_len = p2_with(f0, target = rep(0.25, 3)),
    blp_target_zero = p2_with(f0, target = c(0.5, 0.5, 0, 0)))
  for (cn in names(cases)) {
    for (k in P2_MNL_K) {
      out[[sprintf("err|mnl|%s|%s", cn, k)]] <- p2_cond(p2_mnl_call(k, cases[[cn]]))
    }
  }
  # NL: J = 4 with nests {1, 2}, {3, 4} unless a case says otherwise.
  nest <- c(1L, 1L, 2L, 2L)
  nl_base <- function(ioo, use_asc = TRUE) {
    f <- base(ioo, use_asc)
    th <- f$theta
    f$theta <- c(th[seq_len(K)], 0.7, 0.8, th[-seq_len(K)])
    f$nest_idx <- nest
    f$lambda <- c(0.7, 0.8)
    f
  }
  g0 <- nl_base(FALSE); go <- nl_base(TRUE); gn <- nl_base(FALSE, FALSE)
  nl_cases <- list(
    ok = g0, ok_oo = go, ok_noasc = gn,
    M_zero = p2_with(g0, M = replace(g0$M, 3L, 0L)),
    M_na = p2_with(g0, M = replace(g0$M, 3L, NA_integer_)),
    X_short = p2_with(g0, X = g0$X[-1L, , drop = FALSE]),
    alt_short = p2_with(g0, alt_idx = g0$alt_idx[-1L]),
    w_short = p2_with(g0, weights = g0$weights[-1L]),
    ch_short = p2_with(g0, choice_idx = g0$choice_idx[-1L]),
    alt_zero = p2_with(g0, alt_idx = alt_rep(g0, 5L, 0L)),
    alt_na = p2_with(g0, alt_idx = alt_rep(g0, 5L, NA_integer_)),
    alt_neg = p2_with(g0, alt_idx = alt_rep(g0, 5L, -1L)),
    alt_3e9 = p2_with(g0, alt_idx = replace(as.double(g0$alt_idx), 5L, 3e9)),
    alt_big = p2_with(g0, alt_idx = alt_rep(g0, 5L, J + 1L)),
    alt_big_noasc = p2_with(gn, alt_idx = alt_rep(gn, 5L, J + 1L)),
    alt_frac = p2_with(g0, alt_idx = as.double(g0$alt_idx) + 0.25),
    alt_null = p2_with(g0, alt_idx = NULL),
    ch_past = p2_with(g0, choice_idx = ch_rep(g0, 4L, g0$M[4L] + 1L)),
    ch_zero = p2_with(g0, choice_idx = ch_rep(g0, 4L, 0L)),
    ch_na = p2_with(g0, choice_idx = ch_rep(g0, 4L, NA_integer_)),
    ch_neg = p2_with(g0, choice_idx = ch_rep(g0, 4L, -1L)),
    ch_na_oo = p2_with(go, choice_idx = ch_rep(go, 4L, NA_integer_)),
    ch_past_oo = p2_with(go, choice_idx = ch_rep(go, 4L, go$M[4L] + 1L)),
    nest_short = p2_with(g0, nest_idx = c(1L, 1L, 2L)),
    nest_empty = p2_with(g0, nest_idx = integer(0)),
    nest_zero = p2_with(g0, nest_idx = c(1L, 1L, 0L, 2L)),
    nest_na_last = p2_with(g0, nest_idx = c(1L, 1L, 2L, NA_integer_)),
    nest_neg_last = p2_with(g0, nest_idx = c(1L, 1L, 2L, -1L)),
    nest_dbl = p2_with(g0, nest_idx = as.double(nest)),
    nest_all_single = p2_with(g0, nest_idx = 1:4, theta = c(g0$theta[1:2], g0$theta[-(1:4)])),
    lambda_neg = p2_with(g0, theta = replace(g0$theta, K + 1L, -0.5)),
    theta_short = p2_with(g0, theta = g0$theta[1:3]),
    theta_long_noasc = p2_with(gn, theta = c(gn$theta, 0.1)),
    delta_short = p2_with(g0, theta = g0$theta[-length(g0$theta)]),
    empty_asc = empty(g0), empty_noasc = empty(gn),
    w_zero = p2_with(g0, weights = rep(0, length(g0$weights))),
    var_oob = p2_with(g0, var_idx = 3L),
    blp_target_len = p2_with(g0, target = rep(0.25, 3)),
    blp_target_zero = p2_with(g0, target = c(0.5, 0.5, 0, 0)))
  for (cn in names(nl_cases)) {
    for (k in P2_NL_K) {
      if (k == "blp" && cn %in% c("nest_zero", "nest_na_last", "nest_neg_last",
                                  "nest_short", "nest_empty", "nest_all_single")) {
        next # terminating or Armadillo-internal on main: see p2_blp_abort.R
      }
      out[[sprintf("err|nl|%s|%s", cn, k)]] <- p2_cond(p2_nl_call(k, nl_cases[[cn]]))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Hierarchical Bayes and MNP
# ---------------------------------------------------------------------------

run_hb <- function() {
  out <- list()
  mk <- function(sim_fun, prep_fun, person, seed) {
    sim <- sim_fun(N = 20, T = 2, J = 3, seed = seed, theta = c(0.5, -0.4),
                   sigma_d = 0.5)
    prep_fun(sim$data, "task", "alt", "choice", c("x1", "x2"),
             person_col = if (person) "pid" else NULL,
             alt_covariate_cols = "z1")
  }
  for (person in c(TRUE, FALSE)) {
    d <- mk(simulate_hmnl_data, prepare_hmnl_data, person, 42)
    for (kb in 0:2) {
      args <- list(
        X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
        include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
        rc_dist = as.integer(d$rc_dist), beta_pooled = rep(0, d$K_struct),
        delta_init = rep(0, d$J), theta_init = rep(0, d$P),
        b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
        nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
        theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
        sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
        R = 120, burn = 20, thin = 2, seed = 7, keep_beta_i = kb,
        s_init = 2.38 / sqrt(d$K_struct), accept_target = 0.234)
      out[[sprintf("hb|hmnl_gibbs|person=%d|kb=%d", person, kb)]] <-
        p2_leaf(do.call(choicer:::hmnl_gibbs, args))
    }
    d <- mk(simulate_hmnp_data, prepare_hmnp_data, person, 43)
    for (kb in 0:2) {
      args <- list(
        X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
        include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
        delta_init = rep(0, d$J), theta_init = rep(0, d$P),
        b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
        nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
        theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
        sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
        a0 = 3, s0 = 3, R = 120, burn = 20, thin = 2, seed = 7,
        keep_beta_i = kb)
      out[[sprintf("hb|hmnp_gibbs|person=%d|kb=%d", person, kb)]] <-
        p2_leaf(do.call(choicer:::hmnp_gibbs, args))
    }
  }
  # MNP Gibbs on a small balanced design
  set.seed(44)
  N <- 30; J <- 3
  dt <- data.table::data.table(id = rep(1:N, each = J), alt = rep(1:J, N),
                               x1 = rnorm(N * J), x2 = runif(N * J, -1, 1))
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
  dm <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  out[["hb|mnp_gibbs"]] <- p2_leaf(choicer:::mnp_gibbs(
    X = dm$X, y = dm$y, p = dm$p, beta_bar = rep(0, dm$K), A = 0.01 * diag(dm$K),
    nu = dm$p + 3, V = (dm$p + 3) * diag(dm$p), R = 150, burn = 30, thin = 2,
    seed = 123))
  # End-to-end wrappers (MNL warm start included), timings removed.
  sim <- simulate_hmnl_data(N = 20, T = 2, J = 3, seed = 45, theta = c(0.5, -0.4),
                            sigma_d = 0.5)
  for (kb in c("means", "draws")) {
    out[[sprintf("hb|run_hmnlogit|%s", kb)]] <- p2_val(suppressMessages(run_hmnlogit(
      sim$data, "task", "alt", "choice", c("x1", "x2"), person_col = "pid",
      alt_covariate_cols = "z1", mcmc = list(R = 200, burn = 50, thin = 2, seed = 9),
      chains = 2, keep_beta_i = kb)))
  }
  simp <- simulate_hmnp_data(N = 20, T = 2, J = 3, seed = 46, theta = c(0.5, -0.4),
                             sigma_d = 0.5)
  out[["hb|run_hmnprobit"]] <- p2_val(suppressMessages(run_hmnprobit(
    simp$data, "task", "alt", "choice", c("x1", "x2"), person_col = "pid",
    alt_covariate_cols = "z1", mcmc = list(R = 200, burn = 50, thin = 2, seed = 9))))
  out[["hb|run_mnprobit"]] <- p2_val(suppressMessages(run_mnprobit(
    dt, "id", "alt", "choice", c("x1", "x2"),
    mcmc = list(R = 300, burn = 50, thin = 2, seed = 9))))
  out
}

# ---------------------------------------------------------------------------
# Fits and their post-estimation
# ---------------------------------------------------------------------------

run_fits <- function() {
  out <- list()
  for (ioo in c(TRUE, FALSE)) {
    sim <- simulate_mnl_data(N = 300, J = 4, seed = 47 + ioo, outside_option = ioo)
    d_use <- sim$data
    d_use[, clus := (id - 1L) %/% 10L]
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      args <- list(d_use, "id", "alt", "choice", c("x1", "x2"), se_method = se,
                   cluster_col = if (se == "cluster") "clus" else NULL)
      if (ioo) args <- c(args, list(outside_opt_label = 0L, include_outside_option = TRUE))
      key <- sprintf("fits|mnl|oo=%d|se=%s", ioo, se)
      # By name: the stored call then holds the symbol, not the function
      # (whose body differs between builds; Phase 6D).
      fit <- suppressMessages(do.call("run_mnlogit", args))
      res <- list(fit = p2_leaf(p2_strip(fit)))
      if (se == "hessian") {
        res$pred_p <- p2_val(predict(fit, type = "probabilities"))
        res$pred_s <- p2_val(predict(fit, type = "shares"))
        res$pred_new <- p2_val(predict(fit, type = "probabilities", newdata = d_use))
        res$elas <- p2_val(elasticities(fit, "x1"))
        res$dr <- p2_val(diversion_ratios(fit))
        res$blp <- p2_val(blp(fit, target_shares = predict(fit, type = "shares")))
        res$logsum <- p2_val(logsum(fit))
        res$cs <- p2_val(consumer_surplus(fit, price_var = "x1"))
        res$wtp <- p2_val(wtp(fit, price_var = "x1", attr_vars = "x2"))
        res$gof <- p2_val(gof(fit))
        res$vcov_robust <- p2_val(vcov(fit, type = "robust"))
        res$vcov_bhhh <- p2_val(vcov(fit, type = "bhhh"))
      }
      out[[key]] <- res
    }
  }
  simn <- simulate_nl_data(N = 400, seed = 48)
  dn <- simn$data
  dn[, clus := (id - 1L) %/% 10L]
  for (se in c("hessian", "numeric", "bhhh", "sandwich", "cluster")) {
    key <- sprintf("fits|nl|se=%s", se)
    fit <- suppressMessages(run_nestlogit(
      dn, "id", "j", "choice", c("X", "W"), "nest", se_method = se,
      cluster_col = if (se == "cluster") "clus" else NULL))
    res <- list(fit = p2_leaf(p2_strip(fit)))
    if (se == "hessian") {
      res$pred_p <- p2_val(predict(fit, type = "probabilities"))
      res$pred_s <- p2_val(predict(fit, type = "shares"))
      res$elas <- p2_val(elasticities(fit, "X"))
      res$dr <- p2_val(diversion_ratios(fit))
      res$blp <- p2_val(blp(fit, target_shares = predict(fit, type = "shares")))
      res$logsum <- p2_val(logsum(fit))
      res$wtp <- p2_val(wtp(fit, price_var = "X", attr_vars = "W"))
      res$gof <- p2_val(gof(fit))
      res$vcov_robust <- p2_val(vcov(fit, type = "robust"))
      res$vcov_bhhh <- p2_val(vcov(fit, type = "bhhh"))
    }
    out[[key]] <- res
  }
  out
}
