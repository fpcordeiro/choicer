# The variance of run_mxlogit(concentrate_out_asc = "blp"): the meat kernel
# (mxl_blp_meat_parallel()) against sums formed from the units' scores and
# residuals, and .blp_vcov() against the brute force of
# helper-blp-reference.R (G by numDeriv of the stacked equations), the
# block formulas with solve(), and ML's variance in the logit limit.

# The panel kernels' cells with constants, without the ones where the
# constants absorb mu (an alternative-level W), whose curvature along the
# share constraint is singular
blpv_cells <- function() {
  cells <- Filter(function(cfg) !(cfg$rc_mean && cfg$W_layout == "alt"),
                  mxlp_cell_configs())
  # an alternative-level W without mu
  cells <- c(cells, list(list(seed = 111, rc_dist = c(0L, 1L), rc_correlation = TRUE,
                              rc_mean = FALSE, use_asc = TRUE,
                              include_outside_option = TRUE, W_layout = "alt",
                              weight_type = "person")))
  lapply(cells, function(cfg) {
    cfg$use_asc <- TRUE
    do.call(mxlp_fixture, c(list(name = mxlp_cell_name(cfg)), cfg))
  })
}

blpv_meat <- function(fx, D, supplied, weights = fx$weights, ...) {
  choicer:::mxl_blp_meat_parallel(
    fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, weights, fx$eta,
    fx$rc_dist, D, supplied, fx$rc_correlation, fx$rc_mean, TRUE,
    fx$include_outside_option, Ti = fx[["Ti"]], ...)
}

blpv_kernel <- function(fn, fx, weights = fx$weights, ...) {
  fn(fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, weights, fx$eta,
     fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
     Ti = fx[["Ti"]], ...)
}

# The pieces at fx$theta: index sets, the Jacobian (with K), D, the units'
# weights and their v_u = (m_u, c_u) from the scores kernel and the reference
blpv_parts <- function(fx, tg = NULL) {
  n <- length(fx$theta)
  J_d <- length(blpr_free(fx))
  n_c <- n - J_d
  ci <- seq_len(n_c)
  di <- n_c + seq_len(J_d)
  jac <- blpr_call(fx, variance = TRUE)
  D <- -solve(jac$jacobian[di, di], jac$jacobian[di, ci])
  sc <- choicer:::mxl_scores_parallel(
    fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$eta, fx$rc_dist,
    fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
    Ti = fx[["Ti"]])
  v <- cbind(sc[, ci, drop = FALSE] + sc[, di, drop = FALSE] %*% D,
             blpr_residuals(fx$theta, fx, tg))
  list(n = n, n_c = n_c, ci = ci, di = di, jac = jac, D = D,
       w = mxlp_unit_weights(fx), v = v)
}

test_that("the meat sums the units' projected scores and share residuals", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in c(blpv_cells(), list(mxlp_cross_section(blpv_cells()[[2L]])))) {
    p <- blpv_parts(fx)
    U <- length(p$w)
    # robust: the weights as given (here squared, as for the sandwich)
    m2 <- blpv_meat(fx, p$D, FALSE, weights = fx$weights^2)
    expect_equal(m2$meat, crossprod(p$w * p$v), tolerance = 1e-13, label = fx$name)
    expect_null(m2$A_c)
    # clustered: pairs of units
    cl <- rep(seq_len(ceiling(U / 2)), each = 2L)[seq_len(U)]
    expect_equal(blpv_meat(fx, p$D, FALSE, cluster = cl)$meat,
                 crossprod(rowsum(p$w * p$v, cl)), tolerance = 1e-13,
                 label = fx$name)
    # supplied targets (one market): the residuals' part on the units' own
    # alternatives, and the per-market sums; with the correction, the dense
    # residuals' meat
    J_d <- length(p$di)
    s <- seq_len(J_d) / (J_d * (J_d + 1))
    tg <- list(s = s, market_of_alt = rep(1L, fx$J))
    ps <- blpv_parts(fx, tg)
    T_u <- diff(mxlp_unit_offsets(fx[["Ti"]], length(fx$M)))
    ka <- choicer:::.blp_kappa(p$w, T_u, rep(1L, U), 1L, power = 2)
    ms <- blpv_meat(fx, p$D, TRUE, weights = fx$weights^2, alpha = ka$alpha,
                    market = rep(1L, U))
    vs <- ps$v
    vs[, p$di] <- vs[, p$di] + outer(T_u, s)   # the part on the units' own slots
    expect_equal(ms$meat, crossprod(p$w * vs), tolerance = 1e-13, label = fx$name)
    expect_equal(ms$A_c, matrix(colSums(ka$alpha * vs[, p$ci, drop = FALSE])),
                 tolerance = 1e-13, label = fx$name)
    expect_equal(ms$A_d, colSums(ka$alpha * vs[, p$di, drop = FALSE]),
                 tolerance = 1e-13, label = fx$name)
    O <- choicer:::.blp_meat_correct(
      ms, list(source = "supplied", s = s, market = rep(1L, J_d), b = ka$b), p$n_c)
    expect_equal(O, crossprod(p$w * ps$v), tolerance = 1e-13, label = fx$name)
  }
})

test_that("draw batches, the shared accumulator and draws formed on the fly", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- blpv_cells()[[2L]]   # correlated, outside option, person weights
  p <- blpv_parts(fx)
  U <- length(p$w)
  ka <- choicer:::.blp_kappa(p$w, fx$Ti, rep(1:2, length.out = U), 2L, power = 1)
  base <- blpv_meat(fx, p$D, TRUE, alpha = ka$alpha, market = rep(1:2, length.out = U))
  # units in batches of 5 draws form their probabilities again
  b <- blpv_meat(fx, p$D, TRUE, alpha = ka$alpha, market = rep(1:2, length.out = U),
                 draw_batch = 5L)
  expect_equal(b, base, tolerance = 1e-14)
  # two threads, the shared result
  set_omp_threads(2L)
  expect_equal(blpv_meat(fx, p$D, TRUE, alpha = ka$alpha,
                         market = rep(1:2, length.out = U), acc_bytes = 0),
               base, tolerance = 1e-13)
  cl <- rep(1:3, length.out = U)
  set_omp_threads(1L)
  mk <- rep(1:2, length.out = U)
  c1 <- blpv_meat(fx, p$D, TRUE, cluster = cl, alpha = ka$alpha, market = mk)
  set_omp_threads(2L)
  expect_equal(blpv_meat(fx, p$D, TRUE, cluster = cl, alpha = ka$alpha, market = mk),
               c1, tolerance = 1e-13)
  expect_equal(blpv_meat(fx, p$D, TRUE, cluster = cl, alpha = ka$alpha, market = mk,
                         acc_bytes = 0), c1, tolerance = 1e-13)
  set_omp_threads(1L)
  # a weight that is not finite spreads over the whole meat
  w <- fx$weights
  w[seq_len(fx$Ti[1L])] <- NaN
  expect_true(all(is.nan(blpv_meat(fx, p$D, FALSE, weights = w)$meat)))
  # identity-scrambled Halton draws formed on the fly are the cube's
  g <- choicer:::mxl_blp_meat_parallel(
    fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights,
    array(0, c(fx$K_w, 0L, 0L)), fx$rc_dist, p$D, TRUE, fx$rc_correlation,
    fx$rc_mean, TRUE, TRUE, gen_seed = 0L, gen_scramble = 0L, gen_S = fx$S,
    Ti = fx$Ti, alpha = ka$alpha, market = rep(1:2, length.out = U))
  expect_equal(g, base, tolerance = 1e-13)
})

test_that("the variance is the stacked equations' sandwich", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()  # numerical derivatives of the reference
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  cells <- blpv_cells()
  for (fx in list(cells[[1L]], cells[[2L]], cells[[3L]],
                  mxlp_cross_section(cells[[4L]]))) {
    p <- blpv_parts(fx)
    U <- length(p$w)
    I <- blpv_kernel(choicer:::mxl_hessian_parallel, fx)
    tg <- list(source = "sample")
    cl <- rep(seq_len(ceiling(U / 3)), each = 3L)[seq_len(U)]
    V_h <- suppressMessages(choicer:::.blp_vcov(I, p$jac, "hessian", tg, n_c = p$n_c))
    expect_equal(V_h$vcov, blpr_vcov_bf(fx, "hessian"), tolerance = 1e-7,
                 label = fx$name)
    # standard errors where the variance is positive (the fixtures' theta is
    # no estimate, so some need not be), NA where it is negative
    ok <- !is.na(V_h$se)
    expect_equal(V_h$se[ok]^2, diag(V_h$vcov)[ok])
    expect_true(all(diag(V_h$vcov)[!ok] < 0))
    V_r <- suppressMessages(choicer:::.blp_vcov(
      I, p$jac, "sandwich", tg, blpv_meat(fx, p$D, FALSE, weights = fx$weights^2),
      p$n_c))
    expect_equal(V_r$vcov, blpr_vcov_bf(fx, "sandwich"), tolerance = 1e-7,
                 label = fx$name)
    V_c <- suppressMessages(choicer:::.blp_vcov(
      I, p$jac, "cluster", tg, blpv_meat(fx, p$D, FALSE, cluster = cl), p$n_c))
    expect_equal(V_c$vcov, blpr_vcov_bf(fx, "cluster", cluster = cl),
                 tolerance = 1e-7, label = fx$name)
  }
})

test_that("supplied targets in two markets, with their sampling covariance", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (cfg in list(list(rc_correlation = TRUE, weight_type = "person"),
                   list(rc_dist = c(1L, 0L), rc_mean = TRUE, weight_type = "person"))) {
    fx <- do.call(blpr_two_markets, c(list(seed = 501), cfg))
    p <- blpv_parts(fx)
    s <- c(0.2, 0.3, 0.1, 0.25, 0.15, 0.35)   # within-market inside shares
    W_m <- as.numeric(tapply(fx$weights,
                             blpr_market_of_situations(fx, fx$market_of_alt), sum))
    W_free <- W_m[fx$market_of_alt]
    Vs <- matrix(0, 6, 6)
    for (m in 1:2) {
      j <- which(fx$market_of_alt == m)
      Vs[j, j] <- (diag(s[j]) - s[j] %o% s[j]) / 400
    }
    I <- blpv_kernel(choicer:::mxl_hessian_parallel, fx)
    tg <- function(ka) list(source = "supplied", s = s, market = fx$market_of_alt,
                            b = ka$b, vcov = Vs, W = W_free)
    k1 <- choicer:::.blp_kappa(p$w, fx$Ti, fx$unit_market, 2L, power = 1)
    k2 <- choicer:::.blp_kappa(p$w, fx$Ti, fx$unit_market, 2L, power = 2)
    cl <- rep(seq_len(ceiling(fx$U / 3)), each = 3L)[seq_len(fx$U)]
    cl <- match(cl + 100L * fx$unit_market, unique(cl + 100L * fx$unit_market))
    kc <- choicer:::.blp_kappa(p$w, fx$Ti, fx$unit_market, 2L, cluster = cl)
    meat <- function(w, ka, ...) {
      blpv_meat(fx, p$D, TRUE, weights = w, alpha = ka$alpha,
                market = fx$unit_market, ...)
    }
    tg_r <- list(s = s, market_of_alt = fx$market_of_alt)
    V_h <- suppressMessages(choicer:::.blp_vcov(I, p$jac, "hessian", tg(k1),
                                                meat(fx$weights, k1), p$n_c))
    expect_equal(V_h$vcov, blpr_vcov_bf(fx, "hessian", tg_r, Vs = Vs, W_free = W_free),
                 tolerance = 1e-7)
    V_r <- suppressMessages(choicer:::.blp_vcov(I, p$jac, "sandwich", tg(k2),
                                                meat(fx$weights^2, k2), p$n_c))
    expect_equal(V_r$vcov, blpr_vcov_bf(fx, "sandwich", tg_r, Vs = Vs, W_free = W_free),
                 tolerance = 1e-7)
    V_c <- suppressMessages(choicer:::.blp_vcov(I, p$jac, "cluster", tg(kc),
                                                meat(fx$weights, kc, cluster = cl),
                                                p$n_c))
    expect_equal(V_c$vcov, blpr_vcov_bf(fx, "cluster", tg_r, cluster = cl, Vs = Vs,
                                        W_free = W_free), tolerance = 1e-7)
    # the targets' covariance only adds variance
    V_0 <- suppressMessages(choicer:::.blp_vcov(
      I, p$jac, "hessian", replace(tg(k1), "vcov", list(NULL)),
      meat(fx$weights, k1), p$n_c))
    expect_true(all(diag(V_h$vcov) >= diag(V_0$vcov)))
  }
})

test_that("\"bhhh\" takes the BHHH matrix for the information", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- blpv_cells()[[2L]]
  p <- blpv_parts(fx)
  B <- blpv_kernel(choicer:::mxl_bhhh_parallel, fx)
  Cd <- p$jac$jacobian[p$di, p$di]
  Cc <- p$jac$jacobian[p$di, p$ci]
  It <- B[p$ci, p$ci] + B[p$ci, p$di] %*% p$D + t(p$D) %*% B[p$di, p$ci] +
    t(p$D) %*% B[p$di, p$di] %*% p$D
  V <- suppressMessages(choicer:::.blp_vcov(B, p$jac, "bhhh", list(source = "sample"),
                                            n_c = p$n_c))
  expect_equal(V$vcov, blpr_assemble(B, Cc, Cd, It, matrix(0, p$n_c, length(p$di)),
                                     Cd + p$jac$K), tolerance = 1e-12)
})

test_that("without heterogeneity the variance is the logit's", {
  # W = 0: every draw is the same, so the share equations of the sample's
  # shares are the constants' first-order conditions, K = 0, C = the
  # information's constants' rows, and each type gives ML's variance (L,
  # with no score or curvature, left out)
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  cfgs <- mxlp_cell_configs()[c(1L, 2L)]
  for (cfg in cfgs) {
    cfg$use_asc <- TRUE
    fx0 <- do.call(mxlp_fixture, c(list(name = mxlp_cell_name(cfg)), cfg))
    for (fx in list(fx0, mxlp_cross_section(fx0))) {
      fx$W[] <- 0
      p <- blpv_parts(fx)
      L_idx <- ncol(fx$X) +
        seq_len(if (fx$rc_correlation) fx$K_w * (fx$K_w + 1L) / 2L else fx$K_w)
      keep <- setdiff(seq_len(p$n), L_idx)
      sub <- function(x) x[keep, keep]
      expect_lt(max(abs(p$jac$K)), 1e-25)
      I <- blpv_kernel(choicer:::mxl_hessian_parallel, fx)
      expect_equal(p$jac$jacobian[p$di, ], I[p$di, ], tolerance = 1e-14)
      U <- length(p$w)
      cl <- rep(seq_len(ceiling(U / 2)), each = 2L)[seq_len(U)]
      jk <- list(jacobian = sub(p$jac$jacobian), K = p$jac$K)
      nk <- p$n_c - length(L_idx)
      tg <- list(source = "sample")
      Ii <- solve(sub(I))
      V_h <- suppressMessages(choicer:::.blp_vcov(sub(I), jk, "hessian", tg, n_c = nk))
      expect_equal(V_h$vcov, Ii, tolerance = 1e-13, label = fx$name)
      m2 <- blpv_meat(fx, p$D, FALSE, weights = fx$weights^2)
      V_r <- suppressMessages(choicer:::.blp_vcov(sub(I), jk, "sandwich", tg,
                                                  list(meat = sub(m2$meat)), nk))
      B2 <- sub(blpv_kernel(choicer:::mxl_bhhh_parallel, fx, weights = fx$weights^2))
      expect_equal(V_r$vcov, Ii %*% B2 %*% Ii, tolerance = 1e-13, label = fx$name)
      mc <- blpv_meat(fx, p$D, FALSE, cluster = cl)
      V_c <- suppressMessages(choicer:::.blp_vcov(sub(I), jk, "cluster", tg,
                                                  list(meat = sub(mc$meat)), nk))
      Bc <- sub(choicer:::mxl_cluster_meat_parallel(
        fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, cl,
        fx$eta, fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE,
        fx$include_outside_option, Ti = fx[["Ti"]]))
      expect_equal(V_c$vcov, Ii %*% Bc %*% Ii, tolerance = 1e-13, label = fx$name)
    }
  }
})

test_that("the variance reports what it cannot compute", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- blpv_cells()[[2L]]
  p <- blpv_parts(fx)
  I <- blpv_kernel(choicer:::mxl_hessian_parallel, fx)
  tg <- list(source = "sample")
  # a share Jacobian that is not positive definite
  bad <- p$jac
  bad$jacobian[p$di[1L], p$di[1L]] <- -1
  expect_message(r <- choicer:::.blp_vcov(I, bad, "hessian", tg, n_c = p$n_c),
                 "share Jacobian on the constants is not positive definite")
  expect_null(r$vcov)
  expect_true(all(is.na(r$se)))
  expect_length(r$se, p$n)
  # a value that is not finite
  bad <- p$jac
  bad$K[1L, 1L] <- NaN
  expect_message(r <- choicer:::.blp_vcov(I, bad, "hessian", tg, n_c = p$n_c),
                 "not finite")
  expect_null(r$vcov)
  # no curvature along the constraint
  I0 <- I * 0
  expect_message(r <- choicer:::.blp_vcov(I0, p$jac, "hessian", tg, n_c = p$n_c),
                 "singular")
  expect_null(r$vcov)
})

test_that("the meat kernel checks its inputs", {
  fx <- blpv_cells()[[2L]]
  p <- blpv_parts(fx)
  U <- length(p$w)
  expect_error(blpv_meat(fx, p$D[-1L, , drop = FALSE], FALSE),
               "D must be [0-9]+ x [0-9]+")
  expect_error(blpv_meat(fx, p$D, FALSE, alpha = rep(1, U)), "market is needed")
  expect_error(blpv_meat(fx, p$D, FALSE, alpha = rep(1, U - 1L), market = rep(1L, U - 1L)),
               "one entry per likelihood unit")
  expect_error(blpv_meat(fx, p$D, FALSE, alpha = rep(1, U), market = c(0L, rep(1L, U - 1L))),
               "market codes must be positive integers \\(unit 1 has 0\\)")
  expect_error(blpv_meat(fx, p$D, FALSE, cluster = rep(1L, U - 1L)),
               "cluster length")
  expect_error(blpv_meat(fx, p$D, FALSE, cluster = c(NA, rep(1L, U - 1L))),
               "unit 1 has NA")
  expect_error(choicer:::mxl_blp_meat_parallel(
    fx$theta[-(length(fx$theta))], fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M,
    fx$weights, fx$eta, fx$rc_dist, p$D, FALSE, fx$rc_correlation, fx$rc_mean,
    FALSE, TRUE, Ti = fx$Ti), "needs the constants")
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  expect_error(blpv_meat(fx, p$D, FALSE),
               "The share-matching meat needs about [0-9.]+ [kMG]B even with one thread")
})

test_that("the meat on sparse designs", {
  # helper-mxl-compact.R's designs: a situation listing an alternative twice,
  # units and situations without a free constant, the reference never
  # offered, trailing constants; D any J_d x K_c matrix (the meat is defined
  # for every D)
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in Filter(function(fx) fx$use_asc, cdx_cells())) {
    n <- length(fx$theta)
    J_d <- length(blpr_free(fx))
    n_c <- n - J_d
    set.seed(4)
    D <- matrix(stats::rnorm(J_d * n_c), J_d, n_c)
    sc <- choicer:::mxl_scores_parallel(
      fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$eta, fx$rc_dist,
      fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
      Ti = fx[["Ti"]])
    w <- mxlp_unit_weights(fx)
    for (supplied in c(FALSE, TRUE)) {
      tg <- if (supplied) list(s = rep(0.01, J_d), market_of_alt = rep(1L, fx$J))
      cu <- blpr_residuals(fx$theta, fx, tg)
      if (supplied) {
        T_u <- diff(mxlp_unit_offsets(fx[["Ti"]], length(fx$M)))
        cu <- cu + outer(T_u, rep(0.01, J_d))   # the part on the units' own slots
      }
      v <- cbind(sc[, seq_len(n_c), drop = FALSE] + sc[, -seq_len(n_c), drop = FALSE] %*% D,
                 cu)
      m <- blpv_meat(fx, D, supplied, weights = fx$weights^2)$meat
      expect_equal(m, crossprod(w * v), tolerance = 1e-13,
                   label = paste(fx$name, if (supplied) "supplied" else "sample"))
    }
  }
})

test_that("the model's Omega is exact under the simulated model", {
  # A decision maker's choices are distributed as its simulated likelihood
  # says, P~(y) = (1/S) sum_s prod_t P_ts(y_t). Summed over every choice
  # sequence y with those probabilities: E[c c'] = C_d,u + K_u (the model's
  # variance of the unit's chosen counts), E[s c'] = -C_u' (so that Omega_12
  # = [I D'] E[s c'] = 0), and E[s s'] = E[-d^2 log P~] (the information
  # equality), exactly.
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  fx <- mxlp_fixture("enumeration", 2207, J = 3L, S = 6L, rc_correlation = TRUE,
                     include_outside_option = TRUE, use_asc = TRUE)
  off <- mxlp_unit_offsets(fx$Ti, length(fx$M))
  row_off <- c(0L, cumsum(fx$M))
  n <- length(fx$theta)
  di <- n - 3L + seq_len(3L)
  units <- utils::head(which(fx$Ti <= 3L & fx$Ti >= 2L), 2L)
  expect_length(units, 2L)
  for (u in units) {
    sit <- (off[u] + 1L):off[u + 1L]
    rows <- unlist(lapply(sit, function(t) (row_off[t] + 1L):row_off[t + 1L]))
    fu <- fx
    fu$X <- fx$X[rows, , drop = FALSE]
    fu$W <- fx$W[rows, , drop = FALSE]
    fu$alt_idx <- fx$alt_idx[rows]
    fu$M <- fx$M[sit]
    fu$weights <- rep(1, length(sit))
    fu$Ti <- length(sit)
    fu$eta <- fx$eta[, , u, drop = FALSE]
    fu$choice_idx <- rep(0L, length(sit))
    jac <- blpr_call(fu, variance = TRUE)
    Q <- as.numeric(jac$counts)
    first <- c(0L, cumsum(fu$M))[seq_along(fu$M)]
    seqs <- as.matrix(expand.grid(lapply(fu$M, function(m) 0:m)))
    E_p <- 0
    E_cc <- E_sc <- E_ss <- E_H <- 0
    for (k in seq_len(nrow(seqs))) {
      fu$choice_idx <- as.integer(seqs[k, ])
      p <- exp(-blpv_kernel(choicer:::mxl_loglik_gradient_parallel, fu)$objective)
      ch <- fu$choice_idx > 0L
      cvec <- Q - tabulate(fu$alt_idx[first[ch] + fu$choice_idx[ch]], 3L)
      s <- as.numeric(choicer:::mxl_scores_parallel(
        fu$theta, fu$X, fu$W, fu$alt_idx, fu$choice_idx, fu$M, fu$eta, fu$rc_dist,
        fu$rc_correlation, fu$rc_mean, TRUE, TRUE, Ti = fu$Ti))
      E_p <- E_p + p
      E_cc <- E_cc + p * tcrossprod(cvec)
      E_sc <- E_sc + p * tcrossprod(s, cvec)
      E_ss <- E_ss + p * tcrossprod(s)
      E_H <- E_H + p * blpv_kernel(choicer:::mxl_hessian_parallel, fu)
    }
    expect_equal(E_p, 1, tolerance = 1e-12)
    expect_equal(E_cc, jac$jacobian[di, di] + jac$K, tolerance = 1e-10)
    expect_equal(E_sc, -t(jac$jacobian[di, ]), tolerance = 1e-10)
    expect_equal(E_ss, E_H, tolerance = 1e-10)
  }
})
