# The share inversion's engine (src/blp_engine.h, src/fixed_point.h): its
# passes against the prediction kernels and a finite-difference Jacobian, the
# two-level map and its safeguards (the stalls, and the race with Berry's
# contraction among them), components and the targets' checks,
# the stopping rule and the warning, through the internal solvers
# (*_blp_solve()) and the exported contractions.

# Choice sets over J alternatives in K regions (each situation takes m of its
# region's, plus one from another with probability `border`; the first
# situation of each region offers the whole region).
eng_design <- function(seed, N, J, m, K = 1L, border = 0, Kx = 2L, Kw = 0L) {
  set.seed(seed)
  reg <- ((seq_len(J) - 1L) %% K) + 1L
  sets <- lapply(seq_len(N), function(i) {
    k <- ((i - 1L) %% K) + 1L
    own <- which(reg == k)
    pick <- if (i <= K || length(own) <= m) own else sort(own[sample.int(length(own), m)])
    if (border > 0 && K > 1L && stats::runif(1) < border) {
      others <- which(reg != k)
      pick <- sort(unique(c(pick, others[sample.int(length(others), 1L)])))
    }
    pick
  })
  alt <- as.integer(unlist(sets))
  n <- length(alt)
  list(X = matrix(stats::rnorm(n * Kx), n, Kx), alt = alt, M = lengths(sets),
       w = stats::runif(N, 0.5, 2), J = J, reg = reg,
       W = if (Kw > 0L) matrix(stats::rnorm(n * Kw), n, Kw))
}

eng_free <- function(d, ioo) if (ioo) d else (d - d[1L])[-1L]
eng_beta <- c(0.8, -0.5)

mnl_s <- function(d, delta, ioo, w = d$w) {
  as.numeric(choicer:::mnl_predict_shares(c(eng_beta, eng_free(delta, ioo)), d$X,
                                          d$alt, d$M, w, TRUE, ioo))
}
# (theta holds the dissimilarities of the nests with more than one alternative)
nl_s <- function(d, delta, ioo, nest, lam, beta = eng_beta) {
  as.numeric(choicer:::nl_predict_shares(
    c(beta, lam[tabulate(nest) > 1L], eng_free(delta, ioo)), d$X, d$alt, d$M, d$w,
    nest, TRUE, ioo))
}
eng_ctl <- function(...) choicer:::.blp_control(...)
# The last pass an internal solver kept (control$keep_pass), as plain vectors
eng_pass <- function(r) lapply(r$pass, as.numeric)

# The largest absolute normalized log-share residual at delta: inside
# alternatives, relative to the first without an outside option (one
# component), and the outside share's with one.
eng_resid <- function(s, target, ioo) {
  r <- log(target) - log(s)
  max(abs(if (ioo) r else r - r[1L]))
}

test_that("a pass's shares are the prediction kernels' at the same delta", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- eng_design(1001L, 200L, 6L, 4L)
  nest <- c(1L, 1L, 2L, 2L, 3L, 3L)
  lam <- c(0.6, 0.8, 0.5)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(3L + ioo)
    delta <- stats::rnorm(6L)
    if (!ioo) delta <- delta - delta[1L]
    target <- rep(1 / (6 + ioo), 6 + ioo)
    ctl <- eng_ctl(max_iter = 1L, keep_pass = TRUE)
    r <- choicer:::mnl_blp_solve(delta, target, d$X, eng_beta, d$alt, d$M, d$w,
                                 ioo, ctl)
    s <- mnl_s(d, delta, ioo)
    expect_equal(eng_pass(r)$s, if (ioo) s[-1L] else s, tolerance = 1e-14)
    if (ioo) expect_equal(eng_pass(r)$s0, s[1L], tolerance = 1e-14)
    # the residual at the one evaluation, the outside share's included
    expect_equal(r$residual, eng_resid(s, target, ioo), tolerance = 1e-12)
    r <- choicer:::nl_blp_solve(delta, target, d$X, eng_beta, lam, d$alt, nest,
                                d$M, d$w, ioo, ctl)
    s <- nl_s(d, delta, ioo, nest, lam)
    expect_equal(eng_pass(r)$s, if (ioo) s[-1L] else s, tolerance = 1e-14)
  }
})

# The largest relative error of a against b, entry by entry (entries below
# 1e-8 of the largest measured against that floor: finite differences
# cannot resolve them)
max_rel <- function(a, b) {
  max(abs(a - b) / pmax(abs(b), 1e-8 * max(abs(b)), .Machine$double.xmin))
}

# The share Jacobian by central differences of the prediction kernel
fd_jacobian <- function(sfun, delta, J, ioo) {
  h <- 1e-6
  sapply(seq_len(J), function(k) {
    up <- delta; up[k] <- up[k] + h
    dn <- delta; dn[k] <- dn[k] - h
    (sfun(up) - sfun(dn)) / (2 * h)
  })
}

test_that("a pass's A_jj and co-choice with the base match the share Jacobian", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # repeated codes in some situations
  d <- eng_design(1002L, 150L, 5L, 3L)
  r0 <- c(0L, cumsum(d$M))
  for (t in 1:10) d$alt[r0[t] + 1L] <- d$alt[r0[t] + 2L]
  nest <- c(1L, 1L, 2L, 2L, 2L)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(11L + ioo)
    delta <- stats::rnorm(5L)
    if (!ioo) delta <- delta - delta[1L]
    target <- rep(1 / (5 + ioo), 5 + ioo)
    ctl <- eng_ctl(max_iter = 1L, keep_pass = TRUE)
    for (model in c("mnl", "nl_small", "nl_big")) {
      lam <- switch(model, mnl = NULL, nl_small = c(1e-3, 0.4), nl_big = c(1.4, 0.7))
      sfun <- function(dl) {
        s <- if (is.null(lam)) mnl_s(d, dl, ioo) else nl_s(d, dl, ioo, nest, lam)
        s
      }
      r <- if (is.null(lam)) {
        choicer:::mnl_blp_solve(delta, target, d$X, eng_beta, d$alt, d$M, d$w, ioo, ctl)
      } else {
        choicer:::nl_blp_solve(delta, target, d$X, eng_beta, lam, d$alt, nest, d$M,
                               d$w, ioo, ctl)
      }
      A <- fd_jacobian(sfun, delta, 5L, ioo)  # rows: outside first when ioo
      Ain <- if (ioo) A[-1L, , drop = FALSE] else A
      lab <- sprintf("%s, ioo %d", model, ioo)
      ps <- eng_pass(r)
      expect_lt(max_rel(ps$A, diag(Ain)), 1e-5, label = lab)
      # the co-choice with the base: minus the base's row of the Jacobian
      cb <- if (ioo) -A[1L, ] else -Ain[1L, ]
      idx <- if (ioo) seq_len(5L) else 2:5
      expect_lt(max_rel(ps$C[idx], cb[idx]), 1e-5, label = lab)
      if (identical(model, "nl_big")) {
        # a nest with lambda > 1: the pass also forms the diagonal and the
        # co-choice its probabilities give with lambda = 1, over each
        # situation's totals per code
        P <- as.numeric(choicer:::nl_predict(
          c(eng_beta, lam, eng_free(delta, ioo)), d$X, d$alt, d$M, nest, TRUE,
          ioo)$choice_prob)
        sit <- rep(seq_along(d$M), d$M)
        key <- paste(sit, d$alt)
        tot <- as.numeric(rowsum(P, key, reorder = FALSE))
        first <- !duplicated(key)
        sit_k <- sit[first]
        alt_k <- d$alt[first]
        Pb <- if (ioo) (1 - as.numeric(rowsum(P, sit)))[sit_k] else
          as.numeric(rowsum(P * (d$alt == 1L), sit))[sit_k]
        wk <- d$w[sit_k] / sum(d$w)
        expect_equal(ps$A1, as.numeric(rowsum(wk * tot * (1 - tot), alt_k)),
                     tolerance = 1e-13, label = lab)
        expect_equal(ps$C1, as.numeric(rowsum(wk * tot * Pb, alt_k)),
                     tolerance = 1e-13, label = lab)
      } else {
        expect_null(ps$A1)
      }
    }
  }
})

test_that("a mixed logit pass's A_jj and co-choice match the share Jacobian", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  set.seed(1017L)
  N <- 60L
  J <- 5L
  sets <- lapply(seq_len(N), function(i) sort(sample.int(J, sample(3:5, 1L))))
  # a code listed twice: the reference in some situations, others elsewhere
  for (i in 1:8) {
    st <- sets[[i]]
    sets[[i]] <- if (i <= 4L && 1L %in% st) sort(c(st, 1L)) else sort(c(st, st[length(st)]))
  }
  alt <- as.integer(unlist(sets))
  M <- lengths(sets)
  n <- length(alt)
  X <- matrix(stats::rnorm(2L * n), n, 2L)
  W <- matrix(stats::rnorm(2L * n), n, 2L)
  w <- stats::runif(N, 0.5, 2)
  eta <- array(stats::rnorm(2L * 7L * N), c(2L, 7L, N))
  mu <- c(0.3, -0.2)
  Lp <- c(log(0.7), 0.2, log(0.5))
  for (ioo in c(FALSE, TRUE)) {
    delta <- stats::rnorm(J)
    if (!ioo) delta <- delta - delta[1L]
    sh <- function(dl) {
      as.numeric(choicer:::mxl_predict_shares(c(eng_beta, mu, Lp, eng_free(dl, ioo)), X, W,
                                              alt, M, w, eta, c(0L, 0L), TRUE, TRUE, TRUE,
                                              ioo))
    }
    r <- choicer:::mxl_blp_solve(delta, rep(1 / (J + ioo), J + ioo), X, W, eng_beta, mu,
                                 Lp, alt, M, w, eta, c(0L, 0L), TRUE, TRUE, ioo,
                                 eng_ctl(max_iter = 1L, keep_pass = TRUE), -1L, 1L, 0L, 0)
    ps <- eng_pass(r)
    s0 <- sh(delta)
    expect_equal(ps$s, if (ioo) s0[-1L] else s0, tolerance = 1e-14)
    Jac <- fd_jacobian(sh, delta, J, ioo)
    Ain <- if (ioo) Jac[-1L, , drop = FALSE] else Jac
    lab <- sprintf("ioo %d", ioo)
    expect_lt(max_rel(ps$A, diag(Ain)), 1e-5, label = lab)
    cb <- if (ioo) -Jac[1L, ] else -Ain[1L, ]
    idx <- if (ioo) seq_len(J) else 2:J
    expect_lt(max_rel(ps$C[idx], cb[idx]), 1e-5, label = lab)
  }
})

test_that("a nested logit's diagonal keeps its digits at lambda = 1e-12", {
  # Nest-mates whose utilities differ by -38 to 38 times lambda: the textbook
  # 1 - P(j|g) of the dominant row loses every digit, and with it A_jj, whose
  # within-nest term is (1 - P(j|g)) / lambda. The reference: the same sums
  # in 90-digit decimal arithmetic.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  set.seed(5L)
  N <- 40L
  k <- seq(-38, 38, length.out = N)
  x <- numeric(3L * N)
  x[seq(2L, 3L * N, 3L)] <- k * 1e-12
  x[seq(3L, 3L * N, 3L)] <- stats::rnorm(N)
  w <- stats::runif(N, 0.5, 2)
  alt <- rep(1:3, N)
  M <- rep(3L, N)
  ref <- list(
    list(A = c(5.10359425975835705e+09, 5.10359425973838234e+09, 1.91240862205404621e-01),
         C = c(2.36856963774622162e+11, 5.10359425965274906e+09, 1.05607888115373735e-01),
         s = c(2.41960558034380546e-01, 1.66026894194309160e-01, 5.92012547771310293e-01)),
    list(A = c(3.73051456641028547e+09, 3.73051456638933754e+09, 1.95998161005084576e-01),
         C = c(5.38058247922503508e-02, 3.30698224526153814e-02, 9.79990805021454248e-02),
         s = c(1.56717273525368905e-01, 1.15033101468046184e-01, 4.56499250015452496e-01)))
  for (ioo in c(FALSE, TRUE)) {
    r <- choicer:::nl_blp_solve(c(0, 0, 0.3), rep(1 / (3 + ioo), 3 + ioo),
                                cbind(x, 0), c(1, 0), c(1e-12, 1), alt, c(1L, 1L, 2L), M,
                                w, ioo, eng_ctl(max_iter = 1L, keep_pass = TRUE))
    ps <- eng_pass(r)
    rf <- ref[[ioo + 1L]]
    lab <- sprintf("ioo %d", ioo)
    expect_lt(max_rel(ps$A, rf$A), 1e-12, label = lab)
    expect_lt(max_rel(ps$C, rf$C), 1e-12, label = lab)
    expect_lt(max_rel(ps$s, rf$s), 1e-12, label = lab)
  }
})

test_that("each model recovers a known delta from zeros, in both layouts", {
  # (the step alone; Anderson's and SQUAREM's legs on this design follow below)
  d <- eng_design(1003L, 400L, 8L, 5L, K = 2L, border = 0.1, Kw = 2L)
  nest <- rep(1:4, each = 2L)
  lam <- c(0.5, 0.7, 0.9, 0.6)
  eta <- get_halton_normals(20L, length(d$M), 2L)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(21L + ioo)
    ds <- stats::rnorm(8L)
    if (!ioo) ds <- ds - ds[1L]
    S <- mnl_s(d, ds, ioo)
    est <- blp_contraction(rep(0, 8L), S, d$X, eng_beta, d$alt, d$M, d$w, ioo,
                           tol = 1e-12, method = "simple")
    expect_equal(dim(est), c(8L, 1L))
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    if (!ioo) expect_identical(est[1L], 0)
    S <- nl_s(d, ds, ioo, nest, lam)
    est <- nl_blp_contraction(rep(0, 8L), S, d$X, eng_beta, lam, d$alt, nest, d$M,
                              d$w, ioo, tol = 1e-12, method = "simple")
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    for (rc in list(c(0L, 0L), c(0L, 1L))) {
      th <- c(eng_beta, 0.3, -0.2, log(0.7), 0.2, log(0.5))
      S <- as.numeric(choicer:::mxl_predict_shares(
        c(th, eng_free(ds, ioo)), d$X, d$W, d$alt, d$M, d$w, eta, rc, TRUE, TRUE,
        TRUE, ioo))
      for (start in if (ioo) list(rep(0, 8L)) else list(rep(0, 8L), rep(0, 7L))) {
        est <- mxl_blp_contraction(start, S, d$X, d$W, eng_beta, c(0.3, -0.2),
                                   c(log(0.7), 0.2, log(0.5)), d$alt, d$M, d$w, eta,
                                   rc, TRUE, TRUE, ioo, tol = 1e-12, method = "simple")
        expect_equal(dim(est), c(8L, 1L))
        expect_equal(as.numeric(est), ds, tolerance = 1e-9)
      }
    }
  }
})

test_that("draws formed on the fly give the cube's inversion, store and generate", {
  d <- eng_design(1004L, 120L, 4L, 4L, Kw = 2L)
  g <- list(gen_seed = 5L, gen_scramble = 1L, gen_S = 15L)
  empty <- array(0, c(2L, 0L, 0L))
  th <- c(eng_beta, log(0.8), log(0.6))
  ds <- c(0, 0.4, -0.3, 0.2)
  S <- as.numeric(do.call(choicer:::mxl_predict_shares, c(list(
    c(th, ds[-1L]), d$X, d$W, d$alt, d$M, d$w, empty, c(0L, 0L), FALSE, FALSE,
    TRUE, FALSE), g)))
  est <- do.call(mxl_blp_contraction, c(list(
    rep(0, 4L), S, d$X, d$W, eng_beta, c(0, 0), c(log(0.8), log(0.6)), d$alt, d$M,
    d$w, empty, c(0L, 0L), FALSE, FALSE, FALSE, tol = 1e-12), g))
  expect_equal(as.numeric(est), ds, tolerance = 1e-9)
  # an integer cube is read as doubles
  cube <- array(as.integer(round(get_halton_normals(10L, 120L, 2L))), c(2L, 10L, 120L))
  Sc <- as.numeric(choicer:::mxl_predict_shares(
    c(th, ds[-1L]), d$X, d$W, d$alt, d$M, d$w, cube * 1.0, c(0L, 0L), FALSE, FALSE,
    TRUE, FALSE))
  est <- mxl_blp_contraction(rep(0, 4L), Sc, d$X, d$W, eng_beta, c(0, 0),
                             c(log(0.8), log(0.6)), d$alt, d$M, d$w, cube,
                             c(0L, 0L), FALSE, FALSE, FALSE, tol = 1e-12)
  expect_equal(as.numeric(est), ds, tolerance = 1e-9)
})

test_that("the nested logit converges where Berry's contraction does not", {
  d <- eng_design(1005L, 300L, 9L, 9L)
  nest <- rep(1:3, each = 3L)
  for (lam in c(0.05, 0.3, 1, 1.4)) {
    for (ioo in c(FALSE, TRUE)) {
      set.seed(31L + ioo)
      ds <- stats::rnorm(9L)
      if (!ioo) ds <- ds - ds[1L]
      S <- nl_s(d, ds, ioo, nest, rep(lam, 3L))
      r <- choicer:::nl_blp_solve(rep(0, 9L), S, d$X, eng_beta, rep(lam, 3L), d$alt,
                                  nest, d$M, d$w, ioo, eng_ctl(tol = 1e-10))
      lab <- sprintf("lambda %g, ioo %d", lam, ioo)
      expect_identical(r$code, 0L, label = lab)
      expect_lte(r$residual, 1e-10)
      expect_equal(as.numeric(r$delta), ds, tolerance = 1e-7, label = lab)
      if (lam == 0.3 && !ioo) {
        # Berry's map overshoots: the plain step's Jacobian has eigenvalues
        # near 1 - 1/lambda < -1
        p <- choicer:::nl_blp_solve(rep(0, 9L), S, d$X, eng_beta, rep(lam, 3L),
                                    d$alt, nest, d$M, d$w, ioo,
                                    eng_ctl(tol = 1e-10, max_iter = 300L, map = "plain"))
        expect_false(identical(p$code, 0L))
      }
    }
  }
})

test_that("the common shift converges where the diagonal step alone diverges", {
  # The reproducer of the review of PR #56: three choice situations, in each
  # of which one of alternatives 2-4 dominates, and the reference, 1, is
  # rarely chosen. At the solution, a common shift of the free constants
  # moves every log-share residual against the base at the rate h = 1, while
  # the diagonal step takes it to move each at D = 0.355: alone, that step
  # multiplies the common shift's error by 1 - h / D = -1.82 and diverges,
  # where Berry's contraction converges.
  X <- matrix(c(0, 3, -3, -3, 0, -3, 3, -3, 0, -3, -3, 3), ncol = 1L)
  alt <- rep(1:4, 3L)
  M <- rep(4L, 3L)
  w <- rep(1, 3L)
  ds <- c(0, 1, 1, 1)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, ds[-1L]), X, alt, M, w))
  at <- eng_pass(choicer:::mnl_blp_solve(ds, S, X, 1, alt, M, w, FALSE,
                                         eng_ctl(tol = 0, max_iter = 1L, keep_pass = TRUE)))
  D <- at$A / at$s + at$C / at$s[1L]
  h <- at$C / at$s + sum(at$C[-1L]) / at$s[1L]
  expect_equal((1 - h / D)[-1L], rep(-1.8175, 3L), tolerance = 1e-4)
  # the three models' exports, from zeros
  expect_equal(as.numeric(blp_contraction(rep(0, 4L), S, X, 1, alt, M, w)), ds,
               tolerance = 1e-8)
  expect_equal(as.numeric(nl_blp_contraction(rep(0, 4L), S, X, 1, c(1, 1), alt,
                                             c(1L, 1L, 2L, 2L), M, w)), ds,
               tolerance = 1e-8)
  eta <- array(stats::rnorm(15L), c(1L, 5L, 3L))
  expect_equal(as.numeric(mxl_blp_contraction(rep(0, 4L), S, X, matrix(0, 12L, 1L), 1,
                                              0, log(0.5), alt, M, w, eta, 0L, FALSE,
                                              TRUE, FALSE)), ds, tolerance = 1e-8)
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE, eng_ctl())
  expect_identical(r$code, 0L)
  expect_lte(r$evaluations, 5L)
  # the diagonal step alone never improves on its start; a stall, its
  # residual having grown past twice the best, takes it back to the best
  # point and halves it, and then it contracts
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(map = "jacobi", patience = 0L))
  expect_identical(r$code, 1L)
  expect_identical(as.numeric(r$delta), rep(0, 4L))
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(map = "jacobi", patience = 30L))
  expect_identical(r$code, 0L)
  expect_gte(r$stalls, 1L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-8)
})

test_that("a stall shortens the steps without pulling a climb back", {
  # The reproducer of the re-review of PR #56. Far from the solution the
  # reference, 1, barely responds to the common shift, which sends the
  # constants of 2 and 3 to about -27, where their shares come only from the
  # two situations without the reference; the way back climbs for about 60
  # steps at a residual a little above the best so far. Sent back to that
  # best point at every stall with shorter steps, the iteration stayed there
  # for good (residual 0.0966 after 100,000 evaluations): from it, the
  # two-level step raises the residual however short. Going on from where it
  # is, at half the step, it converges; Berry's contraction takes 624.
  X <- matrix(c(2, -10, -4, -4, -8, 0, 3, 4, 0, 3, 3), ncol = 1L)
  alt <- c(1:4, 1:4, 2L, 2:3)
  M <- c(4L, 4L, 1L, 2L)
  w <- c(1.074, 0.940, 0.553, 0.889)
  ds <- c(0, 2, 0, 3)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, ds[-1L]), X, alt, M, w))
  resid <- function(d) {
    eng_resid(as.numeric(choicer:::mnl_predict_shares(c(1, d[-1L]), X, alt, M, w)), S,
              FALSE)
  }
  trap <- c(0, -26.63545, -28.54366, 2.83895)
  at <- eng_pass(choicer:::mnl_blp_solve(trap, S, X, 1, alt, M, w, FALSE,
                                         eng_ctl(tol = 0, max_iter = 1L, keep_pass = TRUE)))
  r <- log(S) - log(at$s)
  num <- (r - r[1L])[-1L]
  C <- at$C[-1L]
  s <- at$s[-1L]
  h <- C / s + sum(C) / at$s[1L]
  tau <- sum(s * num) / sum(s * h)
  step <- c(0, tau + (num - h * tau) / (at$A[-1L] / s + C / at$s[1L]))
  expect_gt(resid(trap + 1e-5 * step), resid(trap))
  expect_lt(resid(trap + 1e-3 * c(0, num)), resid(trap))  # Berry's step lowers it
  # the three models' exports, from zeros, with the step and its stalls (at
  # tol 1e-8 the constants are within about 1e-7)
  expect_equal(as.numeric(blp_contraction(rep(0, 4L), S, X, 1, alt, M, w, tol = 1e-8,
                                          method = "simple")), ds, tolerance = 1e-6)
  expect_equal(as.numeric(nl_blp_contraction(rep(0, 4L), S, X, 1, c(1, 1), alt,
                                             c(1L, 1L, 2L, 2L), M, w, tol = 1e-8,
                                             method = "simple")),
               ds, tolerance = 1e-6)
  eta <- array(stats::rnorm(20L), c(1L, 5L, 4L))
  expect_equal(as.numeric(mxl_blp_contraction(rep(0, 4L), S, X, matrix(0, 11L, 1L), 1,
                                              0, log(0.5), alt, M, w, eta, 0L, FALSE,
                                              TRUE, FALSE, tol = 1e-8, method = "simple")),
               ds, tolerance = 1e-6)
  # a stall on the climb (one here, 199 evaluations), and no need of Berry's
  # contraction
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE, eng_ctl())
  expect_identical(r$code, 0L)
  expect_gte(r$stalls, 1L)
  expect_identical(r$race, 0L)
  expect_lte(r$evaluations, 400L)
  expect_identical(choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                                           eng_ctl(patience = 0L, window = 0L))$code, 0L)
  expect_identical(choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                                           eng_ctl(method = "original"))$code, 0L)
})

test_that("Berry's contraction joins an iteration that stops halving its residual", {
  # Alternative 2 nearly owns the only situation that offers it (its
  # diagonal is 4.5e-5 at the solution), so the divided step amplifies its
  # small residuals to the 5-unit cap and wrecks the share of 4, its
  # neighbour there; the two oscillate, and the residual decays too slowly to
  # meet tol in 20,000 evaluations. Berry's contraction needs 6.
  X <- matrix(c(6, -5, 0, -3, 8, 5, 7, -6, -10), ncol = 1L)
  alt <- c(5L, 6L, 1L, 3L, 5L, 6L, 2L, 4L, 6L)
  M <- c(2L, 4L, 3L)
  w <- c(1.106, 0.545, 0.985)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, 0, 1, 3, 0, -1), X, alt, M, w))
  r <- choicer:::mnl_blp_solve(rep(0, 6L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(window = 0L, max_iter = 2000L))
  expect_identical(r$code, 1L)
  rb <- choicer:::mnl_blp_solve(rep(0, 6L), S, X, 1, alt, M, w, FALSE,
                                eng_ctl(method = "original"))
  expect_identical(rb$code, 0L)
  # with the race (the default), Berry's iteration from the start converges
  # first: the same point, after at most T + 2E evaluations
  r <- choicer:::mnl_blp_solve(rep(0, 6L), S, X, 1, alt, M, w, FALSE, eng_ctl())
  expect_identical(r$code, 0L)
  expect_true(r$berry)
  expect_gt(r$race, 0L)
  expect_lte(r$evaluations, r$race + 2L * rb$evaluations)
  expect_identical(as.numeric(r$delta), as.numeric(rb$delta))
  # a budget that ends where the race would begin, or just after it
  for (k in r$race + 0:2) {
    rk <- choicer:::mnl_blp_solve(rep(0, 6L), S, X, 1, alt, M, w, FALSE,
                                  eng_ctl(max_iter = k))
    expect_identical(rk$code, 1L)
    expect_identical(rk$evaluations, k)
    expect_identical(dim(rk$delta), c(6L, 1L))
    expect_true(all(is.finite(rk$delta)))
    expect_identical(rk$race, if (k > r$race) r$race else 0L)
  }
  # and the exports (at the solver's tol, with their pass)
  rx <- choicer:::mnl_blp_solve(rep(0, 6L), S, X, 1, alt, M, w, FALSE, eng_ctl(pass = "exp"))
  expect_identical(rx$race, r$race)
  expect_identical(blp_contraction(rep(0, 6L), S, X, 1, alt, M, w, tol = 1e-8,
                                   method = "simple"), rx$delta)
  eta <- array(stats::rnorm(15L), c(1L, 5L, 3L))
  d <- mxl_blp_contraction(rep(0, 6L), S, X, matrix(0, 9L, 1L), 1, 0, log(0.5), alt, M, w,
                           eta, 0L, FALSE, TRUE, FALSE, tol = 1e-8, method = "simple")
  expect_equal(as.numeric(d), as.numeric(r$delta), tolerance = 1e-8)
})

test_that("Berry's contraction goes on again when the step stops halving its residual", {
  # The reproducer of the third review of PR #56: a nested logit (lambda = 0.2
  # in the nest of 1 and 2) whose third situation offers only 2 and 4; their
  # shares elsewhere are small. From zeros the step carries the constants of 2
  # and 4 down together, to about -90, where that situation still gives them
  # their shares and the residual (7.8e-4) no longer responds to their common
  # level; on the way back to the solution, (-3, -5), it rises before it
  # falls. Berry's contraction, joining after evaluation 201, fell behind and
  # stopped, and the step stayed at that point for good (1,000,000
  # evaluations). The contraction now goes on whenever the step alone fails to
  # halve its best residual within 200 iterations, and converges as it would
  # alone (8,367 of its own evaluations; 19,720 in all, above T + 2E = 16,935
  # and within T + 4E = 33,669).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  X <- matrix(c(6, 0, 7, 1, -2, 5, 8, 0, -1, 10, 0, 3, 0, -2, -3), ncol = 1L)
  alt <- c(1:5, 1:5, 2L, 4L, 1:3)
  M <- c(5L, 5L, 2L, 3L)
  w <- c(1.035, 0.524, 2.255, 1.560)
  nest <- c(1L, 1L, 2L, 2L, 3L)
  lam <- c(0.2, 1, 1)
  sh <- function(d) {
    as.numeric(choicer:::nl_predict_shares(c(1, lam[1:2], d[-1L]), X, alt, M, w, nest,
                                           TRUE, FALSE))
  }
  S <- sh(c(0, -3, -5, -5, 1))
  resid <- function(d) eng_resid(sh(d), S, FALSE)
  trap <- c(0, -90.3615118, -4.9997695, -92.3622983, 1.1272149)
  up <- c(0, 1, 0, 1, 0)
  expect_lt(abs(resid(trap + 40 * up) - resid(trap)), 1e-12)
  expect_gt(resid(trap + 80 * up), resid(trap))
  r0 <- choicer:::nl_blp_solve(rep(0, 5L), S, X, 1, lam, alt, nest, M, w, FALSE,
                               eng_ctl(window = 0L, max_iter = 5000L))
  expect_identical(r0$code, 1L)
  expect_lt(max(r0$delta[c(2L, 4L)]), -50)
  rb <- choicer:::nl_blp_solve(rep(0, 5L), S, X, 1, lam, alt, nest, M, w, FALSE,
                               eng_ctl(method = "original", max_iter = 20000L))
  expect_identical(rb$code, 0L)
  r <- choicer:::nl_blp_solve(rep(0, 5L), S, X, 1, lam, alt, nest, M, w, FALSE,
                              eng_ctl(max_iter = 30000L))
  expect_identical(r$code, 0L)
  expect_true(r$berry)
  expect_gt(r$race, 0L)
  expect_gte(r$resumes, 1L)
  expect_lte(r$evaluations, r$race + 4L * rb$evaluations)
  expect_identical(r$delta, rb$delta)
  # and the export (at the solver's tol, with its pass)
  rx <- choicer:::nl_blp_solve(rep(0, 5L), S, X, 1, lam, alt, nest, M, w, FALSE,
                               eng_ctl(max_iter = 30000L, pass = "exp"))
  expect_identical(rx$code, 0L)
  expect_identical(nl_blp_contraction(rep(0, 5L), S, X, 1, lam, alt, nest, M, w,
                                      tol = 1e-8, max_iter = 30000L, method = "simple"),
                   rx$delta)
  # A logit from a sweep of small hard designs, whose first situation offers
  # only 2 and 3: the step's best residual stays at 0.0027 for over 46,000
  # evaluations while their constants come back down from about 2.6 too high
  # through larger residuals. The contraction, joining after evaluation 201,
  # halves its best residual over its first 100 evaluations (3.2 to 0.0042)
  # but not over its next 100 (0.0035, above the step's 0.0027), and stops at
  # evaluation 602. The step then has 200 iterations alone, after which the
  # contraction goes on, once: with a budget of 803 evaluations, not 802,
  # since it needs one left. Through the logit's and the mixed logit's exports
  # (a zero random design)
  X <- matrix(c(9, 2, 3, -8, -5, 0, 3, -9, -1, 4, 0, 7), ncol = 1L)
  alt <- c(2L, 3L, 1:5, 4L, 2L, 4L, 5L, 4L)
  M <- c(2L, 5L, 1L, 3L, 1L)
  w <- c(0.558, 1.128, 0.697, 0.536, 1.152)
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, -2, 0, 3, 2), X, alt, M, w))
  rb <- choicer:::mnl_blp_solve(rep(0, 5L), S, X, 1, alt, M, w, FALSE,
                                eng_ctl(method = "original", max_iter = 100000L))
  expect_identical(rb$code, 0L)
  r <- choicer:::mnl_blp_solve(rep(0, 5L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(max_iter = 200000L))
  expect_identical(r$code, 0L)
  expect_true(r$berry)
  expect_gte(r$resumes, 1L)
  expect_lte(r$evaluations, r$race + 4L * rb$evaluations)
  expect_identical(r$delta, rb$delta)
  for (k in c(802L, 803L)) {
    rk <- choicer:::mnl_blp_solve(rep(0, 5L), S, X, 1, alt, M, w, FALSE,
                                  eng_ctl(max_iter = k))
    expect_identical(c(rk$race, rk$resumes), c(r$race, if (k == 802L) 0L else 1L))
  }
  rx <- choicer:::mnl_blp_solve(rep(0, 5L), S, X, 1, alt, M, w, FALSE,
                                eng_ctl(max_iter = 200000L, pass = "exp"))
  expect_identical(rx$code, 0L)
  expect_identical(blp_contraction(rep(0, 5L), S, X, 1, alt, M, w, tol = 1e-8,
                                   max_iter = 200000L, method = "simple"), rx$delta)
  eta <- array(stats::rnorm(25L), c(1L, 5L, 5L))
  d <- mxl_blp_contraction(rep(0, 5L), S, X, matrix(0, 12L, 1L), 1, 0, log(0.5), alt, M,
                           w, eta, 0L, FALSE, TRUE, FALSE, tol = 1e-8, max_iter = 200000L,
                           method = "simple")
  expect_equal(as.numeric(d), as.numeric(r$delta), tolerance = 1e-8)
})

test_that("an alternative that complements the base has a diagonal of at least 1", {
  # Situations of two kinds: {1, 2, 3}, where 3 dominates, and {2, 3}, where 2
  # does; 1 and 2 share a nest with lambda = 2.5. Where 1 is offered, raising
  # delta_2 raises its nest's inclusive value and 1's share with it, so the
  # co-choice of 2 with the reference 1 is negative and the exact diagonal
  # nets out to about 0 (here -0.2): unbounded, its steps leave the iteration
  # at its start; floored at 1, it converges at least as fast as Berry's
  # contraction.
  N <- 200L
  sets <- lapply(seq_len(N), function(i) if (i %% 2L) 1:3 else 2:3)
  alt <- as.integer(unlist(sets))
  M <- lengths(sets)
  typ <- rep(seq_len(N) %% 2L, M)
  X <- cbind(ifelse(typ == 1L & alt == 3L, 5, ifelse(typ == 0L & alt == 2L, 2, 0)), 0)
  d <- list(X = X, alt = alt, M = M, w = rep(1, N))
  nest <- c(1L, 1L, 2L)
  lam <- c(2.5, 1)
  bt <- c(1, 0)
  ds <- c(0, 0.5, -0.2)
  S <- nl_s(d, ds, FALSE, nest, lam, bt)
  r <- choicer:::nl_blp_solve(rep(0, 3L), S, X, bt, lam, alt, nest, M, d$w, FALSE,
                              eng_ctl(tol = 1e-10, keep_pass = TRUE))
  expect_identical(r$code, 0L)
  expect_gt(r$fallbacks, 0L)
  expect_lt(eng_pass(r)$C[2L], 0)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-7)
  p <- choicer:::nl_blp_solve(rep(0, 3L), S, X, bt, lam, alt, nest, M, d$w, FALSE,
                              eng_ctl(tol = 1e-10, map = "plain"))
  expect_lte(r$evaluations, p$evaluations)
})

test_that("a nest with a large lambda keeps a logit's diagonal as a floor", {
  # With lambda = 3 in nests of ten and an outside option, the exact diagonal
  # is about a third of what a common shift of a nest's constants needs, so
  # its step would overshoot that shift threefold and diverge; Berry's
  # contraction stalls (its rate on the shifts within a nest is 1 - 1/lambda,
  # and slower with the outside option).
  set.seed(1016L)
  J <- 30L
  N <- 800L
  sets <- lapply(seq_len(N), function(i) sort(sample.int(J, sample(5:J, 1L))))
  alt <- as.integer(unlist(sets))
  d <- list(X = cbind(stats::rnorm(length(alt)), 0), alt = alt, M = lengths(sets),
            w = rep(1, N))
  nest <- ((seq_len(J) - 1L) %/% 10L) + 1L
  lam <- rep(3, 3L)
  bt <- c(1, 0)
  ds <- stats::rnorm(J)
  S <- nl_s(d, ds, TRUE, nest, lam, bt)
  r <- choicer:::nl_blp_solve(rep(0, J), S, d$X, bt, lam, alt, nest, d$M, d$w, TRUE,
                              eng_ctl(tol = 1e-10))
  expect_identical(r$code, 0L)
  expect_lt(r$evaluations, 200L)
  expect_gt(r$fallbacks, 0L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-7)
  p <- choicer:::nl_blp_solve(rep(0, J), S, d$X, bt, lam, alt, nest, d$M, d$w, TRUE,
                              eng_ctl(tol = 1e-10, map = "plain"))
  expect_false(identical(p$code, 0L))
})

test_that("the step cap rescues a far noisy start in a segmented market", {
  # (slow for every map: the regions' levels relative to one another meet
  # through few situations; about 3,650 evaluations here, of which Berry's
  # contraction, racing the step until it falls behind, takes about 600; its
  # own 5,500)
  d <- eng_design(1007L, 1200L, 24L, 4L, K = 4L, border = 0.02)
  set.seed(41L)
  ds <- stats::rnorm(24L)
  ds <- ds - ds[1L]
  S <- mnl_s(d, ds, FALSE)
  set.seed(42L)
  start <- ds + stats::rnorm(24L, 0, 3)
  r <- choicer:::mnl_blp_solve(start, S, d$X, eng_beta, d$alt, d$M, d$w, FALSE,
                               eng_ctl(tol = 1e-10, max_iter = 5000L))
  expect_identical(r$code, 0L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-6)
})

test_that("an outside share that underflows takes the step without its correction", {
  d <- eng_design(1008L, 200L, 5L, 5L)
  ds <- c(0.3, -0.2, 0.1, 0.4, -0.5)
  S <- mnl_s(d, ds, TRUE)
  # (Anderson acceleration takes that step unextrapolated until the share
  # returns: from this start it stayed there for 20,000 evaluations when it
  # extrapolated it)
  for (method in c("simple", "anderson", "squarem")) {
    r <- choicer:::mnl_blp_solve(rep(800, 5L), S, d$X, eng_beta, d$alt, d$M, d$w,
                                 TRUE, eng_ctl(tol = 1e-10, method = method))
    expect_identical(r$code, 0L, label = method)
    expect_gt(r$fallbacks, 0L)
    expect_equal(as.numeric(r$delta), ds, tolerance = 1e-8, label = method)
  }
})

test_that("components hold their references; their targets must fit their weight", {
  d <- eng_design(1009L, 300L, 6L, 3L, K = 2L, border = 0)
  ds <- c(0, 0.4, -0.2, 1.5, 0.3, -0.6)
  S <- mnl_s(d, ds, FALSE)
  start <- c(0, 0, 0, 2, 0, 0)  # component 2's reference, alternative 2, held at 0
  est <- blp_contraction(start, S, d$X, eng_beta, d$alt, d$M, d$w, tol = 1e-12)
  in1 <- d$reg == 1L
  expect_equal(as.numeric(est)[in1], ds[in1], tolerance = 1e-9)
  expect_identical(as.numeric(est)[2L], 0)
  expect_equal(as.numeric(est)[!in1] - est[2L], ds[!in1] - ds[2L], tolerance = 1e-9)
  bad <- S
  bad[in1] <- S[in1] * 1.1
  bad[!in1] <- S[!in1] * (1 - sum(bad[in1])) / sum(S[!in1])
  expect_error(blp_contraction(rep(0, 6L), bad, d$X, eng_beta, d$alt, d$M, d$w),
               "their choice situations weigh")
})

test_that("targets the data cannot attain stop before any pass", {
  d <- eng_design(1010L, 200L, 6L, 3L, K = 2L, border = 0.2)
  S <- mnl_s(d, rep(0, 6L), FALSE)
  expect_error(blp_contraction(rep(0, 6L), S * 1.01, d$X, eng_beta, d$alt, d$M, d$w),
               "target_shares must sum to 1; they sum to")
  expect_error(blp_contraction(rep(0, 6L), c(0.5, rep(0.1, 6L)), d$X, eng_beta,
                               d$alt, d$M, d$w, TRUE),
               "must sum to 1, the outside option's share included")
  # an alternative's target above the weight of the situations offering it
  offer <- tapply(d$w[rep(seq_along(d$M), d$M)], d$alt, sum) / sum(d$w)
  tg <- rep((1 - 1.1 * offer[6L]) / 5, 6L)
  tg[6L] <- 1.1 * offer[6L]
  expect_error(blp_contraction(rep(0, 6L), tg, d$X, eng_beta, d$alt, d$M, d$w),
               "is not below the weight of the choice situations that offer it")
  # below the weight of the situations offering it alone
  alt <- c(1L, 1L, 2L, 1L, 2L, 3L, 2L, 3L)
  M <- c(1L, 2L, 3L, 2L)
  X <- matrix(stats::rnorm(16L), 8L, 2L)
  expect_error(blp_contraction(rep(0, 3L), c(0.2, 0.4, 0.4), X, eng_beta, alt, M,
                               rep(1, 4L)),
               "is not above the weight of the choice situations that offer it alone")
  # an alternative offered only by a situation of zero weight
  expect_error(blp_contraction(rep(0, 3L), c(0.5, 0.3, 0.2), X, eng_beta, alt, M,
                               c(1, 1, 1, 0)),
               NA)
  expect_error(blp_contraction(rep(0, 4L), rep(0.25, 4L), rbind(X, 0), eng_beta,
                               c(alt, 4L), c(M[1:3], 3L), c(1, 1, 1, 0)),
               "Alternative 4 is offered in no choice situation of positive weight")
  # an alternative alone in its component (offered only by itself): its target
  # is its weight, and its constant is held
  alt2 <- c(1L, 2L, 1L, 2L, 3L)
  M2 <- c(2L, 2L, 1L)
  X2 <- matrix(stats::rnorm(10L), 5L, 2L)
  w2 <- c(1, 1, 2)
  s12 <- mnl_s(list(X = X2, alt = alt2, M = M2, w = w2), c(0, 0.3, 0), FALSE)
  est <- blp_contraction(c(0, 0, 5), s12, X2, eng_beta, alt2, M2, w2, tol = 1e-12)
  expect_equal(as.numeric(est), c(0, 0.3, 5), tolerance = 1e-9)
})

test_that("inputs that used to give NaN or drift are errors", {
  d <- eng_design(1011L, 100L, 4L, 4L)
  S <- rep(0.25, 4L)
  expect_error(blp_contraction(c(0, NA, 0, 0), S, d$X, eng_beta, d$alt, d$M, d$w),
               "delta must be finite")
  expect_error(blp_contraction(rep(0, 4L), c(0.25, NaN, 0.25, 0.25), d$X, eng_beta,
                               d$alt, d$M, d$w), "strictly positive", fixed = TRUE)
  expect_error(blp_contraction(rep(0, 4L), S, d$X, eng_beta, d$alt, d$M,
                               replace(d$w, 3L, -1)),
               "weights must be finite and non-negative")
  expect_error(blp_contraction(rep(0, 4L), S, d$X, eng_beta, d$alt, d$M, d$w,
                               tol = -1), "`tol` must be a single non-negative")
  expect_error(blp_contraction(rep(0, 4L), S, d$X, eng_beta, d$alt, d$M, d$w,
                               max_iter = 1.5), "`max_iter` must be a whole number")
  nest <- c(1L, 1L, 2L, 2L)
  for (lam in list(c(0.5, 0), c(-0.5, 1), c(Inf, 1), c(NaN, 1))) {
    expect_error(nl_blp_contraction(rep(0, 4L), S, d$X, eng_beta, lam, d$alt, nest,
                                    d$M, d$w), "lambda must be positive and finite")
  }
  for (damp in list(0, 1.5, -1, NA_real_, c(0.5, 0.5))) {
    expect_error(nl_blp_contraction(rep(0, 4L), S, d$X, eng_beta, c(0.5, 0.5),
                                    d$alt, nest, d$M, d$w, damping = damp),
                 "`damping` must be a single number in (0, 1]", fixed = TRUE)
  }
})

test_that("max_iter counts evaluations; an unconverged call warns with its record", {
  d <- eng_design(1012L, 300L, 6L, 6L)
  ds <- c(0, 0.5, -0.3, 0.2, 0.8, -0.4)
  S <- mnl_s(d, ds, FALSE)
  w0 <- expect_warning(
    est <- blp_contraction(c(1, 2, 3, 4, 5, 6), S, d$X, eng_beta, d$alt, d$M, d$w,
                           max_iter = 0L),
    class = "choicer_blp_unconverged")
  expect_identical(as.numeric(est), c(0, 1, 2, 3, 4, 5))  # the start, normalized
  cv <- w0$convergence
  expect_identical(cv$status, "max_iter")
  expect_identical(cv$evaluations, 0L)
  expect_true(is.na(cv$residual))
  w3 <- expect_warning(
    est3 <- blp_contraction(rep(0, 6L), S, d$X, eng_beta, d$alt, d$M, d$w,
                            tol = 1e-15, max_iter = 3L),
    class = "choicer_blp_unconverged")
  cv <- w3$convergence
  expect_identical(cv$evaluations, 3L)
  expect_false(cv$converged)
  # the residual reported is the returned point's (to rounding: a difference
  # of logs, so absolute)
  expect_lt(abs(cv$residual - eng_resid(mnl_s(d, as.numeric(est3), FALSE), S, FALSE)),
            1e-13)
  # converged: silent, within tol
  expect_silent(est <- blp_contraction(rep(0, 6L), S, d$X, eng_beta, d$alt, d$M,
                                       d$w, tol = 1e-10))
  expect_lte(eng_resid(mnl_s(d, as.numeric(est), FALSE), S, FALSE), 1e-10)
})

test_that("a start whose shares underflow stops with the warning", {
  d <- eng_design(1013L, 100L, 4L, 4L)
  S <- mnl_s(d, c(0, 0.2, 0.1, -0.3), FALSE)
  w <- expect_warning(
    est <- blp_contraction(c(0, 0, 0, -1000), S, d$X, eng_beta, d$alt, d$M, d$w),
    class = "choicer_blp_unconverged")
  expect_identical(w$convergence$status, "non_finite_start")
  expect_identical(as.numeric(est), c(0, 0, 0, -1000))
})

test_that("tiny targets and weights are met in log shares", {
  d <- eng_design(1014L, 300L, 6L, 6L)
  ds <- c(0, 0.5, -0.3, 0.2, 0.8, -24)
  S <- mnl_s(d, ds, FALSE)
  expect_lt(S[6L], 1e-9)
  est <- blp_contraction(rep(0, 6L), S, d$X, eng_beta, d$alt, d$M, d$w, tol = 1e-12)
  expect_equal(as.numeric(est), ds, tolerance = 1e-9)
  w0 <- replace(d$w, c(3L, 9L), 0)
  S0 <- mnl_s(d, ds, FALSE, w0)
  est <- blp_contraction(rep(0, 6L), S0, d$X, eng_beta, d$alt, d$M, w0, tol = 1e-12)
  expect_equal(as.numeric(est), ds, tolerance = 1e-9)
})

test_that("two threads agree with one to the tolerance", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled))
  on.exit(set_omp_threads(2L), add = TRUE)
  d <- eng_design(1015L, 3000L, 10L, 5L, K = 2L, border = 0.2)
  set.seed(51L)
  ds <- stats::rnorm(10L)
  ds <- ds - ds[1L]
  S <- mnl_s(d, ds, FALSE)
  for (method in c("simple", "anderson", "squarem")) {
    set_omp_threads(1L)
    e1 <- blp_contraction(rep(0, 10L), S, d$X, eng_beta, d$alt, d$M, d$w,
                          tol = 1e-12, method = method)
    set_omp_threads(2L)
    e2 <- blp_contraction(rep(0, 10L), S, d$X, eng_beta, d$alt, d$M, d$w,
                          tol = 1e-12, method = method)
    expect_equal(e2, e1, tolerance = 1e-10, label = method)
  }
})

test_that("Anderson acceleration recovers a known delta in each model and layout", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # (identical results at one thread)
  d <- eng_design(1003L, 400L, 8L, 5L, K = 2L, border = 0.1, Kw = 2L)
  nest <- rep(1:4, each = 2L)
  lam <- c(0.5, 0.7, 0.9, 0.6)
  eta <- get_halton_normals(20L, length(d$M), 2L)
  th <- c(eng_beta, 0.3, -0.2, log(0.7), 0.2, log(0.5))
  for (ioo in c(FALSE, TRUE)) {
    set.seed(21L + ioo)
    ds <- stats::rnorm(8L)
    if (!ioo) ds <- ds - ds[1L]
    S <- mnl_s(d, ds, ioo)
    est <- blp_contraction(rep(0, 8L), S, d$X, eng_beta, d$alt, d$M, d$w, ioo,
                           tol = 1e-12, method = "anderson")
    expect_equal(dim(est), c(8L, 1L))
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    if (!ioo) expect_identical(est[1L], 0)
    r <- choicer:::mnl_blp_solve(rep(0, 8L), S, d$X, eng_beta, d$alt, d$M, d$w, ioo,
                                 eng_ctl(tol = 1e-12, method = "anderson", pass = "exp"))
    expect_identical(r$code, 0L)
    expect_identical(r$race, 0L)
    expect_lte(r$residual, 1e-12)
    expect_identical(r$delta, est)
    S <- nl_s(d, ds, ioo, nest, lam)
    est <- nl_blp_contraction(rep(0, 8L), S, d$X, eng_beta, lam, d$alt, nest, d$M,
                              d$w, ioo, tol = 1e-12, method = "anderson")
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    S <- as.numeric(choicer:::mxl_predict_shares(
      c(th, eng_free(ds, ioo)), d$X, d$W, d$alt, d$M, d$w, eta, c(0L, 1L), TRUE, TRUE,
      TRUE, ioo))
    est <- mxl_blp_contraction(rep(0, 8L), S, d$X, d$W, eng_beta, c(0.3, -0.2),
                               c(log(0.7), 0.2, log(0.5)), d$alt, d$M, d$w, eta,
                               c(0L, 1L), TRUE, TRUE, ioo, tol = 1e-12,
                               method = "anderson")
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
  }
})

# The race reproducers above, and a sweep design, with Anderson acceleration.
# Its least squares amplify last-bit differences, so the paths (who wins the
# race, after how many evaluations) can differ across platforms: the first
# test holds anywhere, the second pins the development machine's paths and is
# skipped on CRAN and CI (as the package's other path snapshots are).
aa_cases <- function() {
  list(
    # the re-review's: the two-level step climbs back for about 60 steps
    climb = list(X = matrix(c(2, -10, -4, -4, -8, 0, 3, 4, 0, 3, 3), ncol = 1L),
                 alt = c(1:4, 1:4, 2L, 2:3), M = c(4L, 4L, 1L, 2L),
                 w = c(1.074, 0.940, 0.553, 0.889), ds = c(0, 2, 0, 3)),
    # alternative 2 nearly owns its situation, and the divided step oscillates
    oscillate = list(X = matrix(c(6, -5, 0, -3, 8, 5, 7, -6, -10), ncol = 1L),
                     alt = c(5L, 6L, 1L, 3L, 5L, 6L, 2L, 4L, 6L), M = c(2L, 4L, 3L),
                     w = c(1.106, 0.545, 0.985), ds = c(0, 0, 1, 3, 0, -1)),
    # a logit from a sweep of small hard designs (two situations)
    sweep = list(X = matrix(c(10, -3, 2, 2, 6, 7, -7, 10, 7, -9, -10), ncol = 1L),
                 alt = c(1:6, 1:3, 5:6), M = 6:5, w = c(0.521, 0.926),
                 ds = c(0, -2, 2, 1, 0, -1)))
}
aa_solve <- function(cs, ...) {
  S <- as.numeric(choicer:::mnl_predict_shares(c(1, cs$ds[-1L]), cs$X, cs$alt, cs$M,
                                               cs$w))
  choicer:::mnl_blp_solve(rep(0, length(cs$ds)), S, cs$X, 1, cs$alt, cs$M, cs$w, FALSE,
                          eng_ctl(...))
}
# The third review's nested logit (above): the two-level step carries the
# constants of 2 and 4 down to where the residual no longer responds to their
# level
aa_nl_trap <- function(...) {
  X <- matrix(c(6, 0, 7, 1, -2, 5, 8, 0, -1, 10, 0, 3, 0, -2, -3), ncol = 1L)
  alt <- c(1:5, 1:5, 2L, 4L, 1:3)
  M <- c(5L, 5L, 2L, 3L)
  w <- c(1.035, 0.524, 2.255, 1.560)
  nest <- c(1L, 1L, 2L, 2L, 3L)
  lam <- c(0.2, 1, 1)
  S <- as.numeric(choicer:::nl_predict_shares(c(1, lam[1:2], -3, -5, -5, 1), X, alt, M,
                                              w, nest, TRUE, FALSE))
  choicer:::nl_blp_solve(rep(0, 5L), S, X, 1, lam, alt, nest, M, w, FALSE, eng_ctl(...))
}

test_that("Anderson acceleration converges on the race reproducers", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (nm in names(aa_cases())) {
    r <- aa_solve(aa_cases()[[nm]], method = "anderson", max_iter = 20000L)
    # (the shares within tol; the oscillating design's constants are
    # determined only to about 0.02 by them, its conditioning being 2e6)
    expect_identical(r$code, 0L, label = nm)
    expect_lte(r$residual, 1e-8)
  }
  # (Berry's contraction converges in 8,367 evaluations, so the race
  # converges within T + 4E, T being the evaluations before it joins)
  expect_identical(aa_nl_trap(method = "anderson", max_iter = 40000L)$code, 0L)
})

test_that("Anderson acceleration's paths on the race reproducers", {
  skip_on_cran()  # (the development machine's paths: see aa_cases())
  skip_on_ci()
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  cs <- aa_cases()
  # Within T + 4E wherever the race runs (no trial point failed here)
  for (nm in names(cs)) {
    rb <- aa_solve(cs[[nm]], method = "original", max_iter = 20000L)
    r <- aa_solve(cs[[nm]], method = "anderson", max_iter = 20000L)
    if (r$race > 0L) expect_lte(r$evaluations, r$race + 4L * rb$evaluations)
  }
  rb <- aa_nl_trap(method = "original", max_iter = 20000L)
  r <- aa_nl_trap(method = "anderson", max_iter = 40000L)
  expect_gt(r$race, 0L)
  expect_lte(r$evaluations, r$race + 4L * rb$evaluations)
  # The climb: restarts that keep it, a few dozen evaluations (the two-level
  # step's stalls: about 200), and Phase R's restart from the best point
  # gets there too
  r <- aa_solve(cs$climb, method = "anderson")
  rs <- aa_solve(cs$climb)
  expect_gte(r$restarts, 1L)
  expect_identical(r$race, 0L)
  expect_lt(r$evaluations, rs$evaluations)
  expect_identical(aa_solve(cs$climb, method = "anderson", restart_best = TRUE)$code, 0L)
  # The oscillation: the two-level step needs Berry's contraction; Anderson's
  # iteration converges by itself, in fewer evaluations than the race;
  # restarting from its best point every time, it stays there without the
  # race
  rs <- aa_solve(cs$oscillate)
  r <- aa_solve(cs$oscillate, method = "anderson")
  expect_gt(rs$race, 0L)
  expect_identical(r$race, 0L)
  expect_lt(r$evaluations, rs$evaluations)
  rb <- aa_solve(cs$oscillate, method = "anderson", restart_best = TRUE, window = 0L)
  expect_identical(rb$code, 1L)
  expect_gt(rb$restarts, 100L)
  # The sweep design: extrapolating fails to make progress, and without the
  # race the iteration stays near a residual of 0.004; with it, Berry's
  # contraction joins after 201 evaluations and converges first, at its own
  # point, within T + 2E
  r0 <- aa_solve(cs$sweep, method = "anderson", window = 0L, max_iter = 5000L)
  expect_identical(r0$code, 1L)
  expect_gt(r0$residual, 1e-3)
  rb <- aa_solve(cs$sweep, method = "original", max_iter = 20000L)
  r <- aa_solve(cs$sweep, method = "anderson", max_iter = 20000L)
  expect_true(r$berry)
  expect_gt(r$race, 0L)
  expect_lte(r$evaluations, r$race + 2L * rb$evaluations)
  expect_identical(r$delta, rb$delta)
})

test_that("Anderson's least squares drop the differences they cannot use", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- eng_design(1016L, 600L, 12L, 4L, K = 3L, border = 0.05)
  set.seed(61L)
  ds <- stats::rnorm(12L)
  ds <- ds - ds[1L]
  S <- mnl_s(d, ds, FALSE)
  run <- function(method = "anderson", ...) {
    choicer:::mnl_blp_solve(rep(0, 12L), S, d$X, eng_beta, d$alt, d$M, d$w, FALSE,
                            eng_ctl(tol = 1e-12, max_iter = 5000L, method = method,
                                    settings = list(...)))
  }
  # A condition limit of 1 keeps only the newest difference (unless the
  # triangular factor's diagonal is constant): the iterates of memory 1,
  # which needs Berry's contraction here (three regions that share few
  # alternatives)
  r1 <- run(memory = 1L)
  rc <- run(max_condition = 1)
  expect_identical(r1$code, 0L)
  expect_identical(rc[c("delta", "evaluations", "restarts", "race")],
                   r1[c("delta", "evaluations", "restarts", "race")])
  # memory 10 (the default): a few dozen evaluations, where the two-level
  # step alone needs over a thousand
  r10 <- run()
  rs <- run("simple")
  expect_identical(r10$code, 0L)
  expect_identical(r10$race, 0L)
  expect_lt(10L * r10$evaluations, rs$evaluations)
  expect_lt(r10$evaluations, r1$evaluations)
  expect_equal(as.numeric(r10$delta), ds, tolerance = 1e-9)
  # Memory beyond the free constants (two alternatives: one) is cut to them
  ds2 <- c(0, 0.7)
  X2 <- matrix(c(0.5, -0.5, 1, 0.2), ncol = 1L)
  S2 <- as.numeric(choicer:::mnl_predict_shares(c(1, 0.7), X2, c(1:2, 1:2), c(2L, 2L),
                                                c(1, 1)))
  r2 <- choicer:::mnl_blp_solve(c(0, 0), S2, X2, 1, c(1:2, 1:2), c(2L, 2L), c(1, 1),
                                FALSE, eng_ctl(tol = 1e-12, method = "anderson"))
  expect_identical(r2$code, 0L)
  expect_equal(as.numeric(r2$delta), ds2, tolerance = 1e-10)
})

test_that("SQUAREM recovers a known delta in each model and layout", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # (identical results at one thread)
  d <- eng_design(1003L, 400L, 8L, 5L, K = 2L, border = 0.1, Kw = 2L)
  nest <- rep(1:4, each = 2L)
  lam <- c(0.5, 0.7, 0.9, 0.6)
  eta <- get_halton_normals(20L, length(d$M), 2L)
  th <- c(eng_beta, 0.3, -0.2, log(0.7), 0.2, log(0.5))
  for (ioo in c(FALSE, TRUE)) {
    set.seed(21L + ioo)
    ds <- stats::rnorm(8L)
    if (!ioo) ds <- ds - ds[1L]
    S <- mnl_s(d, ds, ioo)
    est <- blp_contraction(rep(0, 8L), S, d$X, eng_beta, d$alt, d$M, d$w, ioo,
                           tol = 1e-12, method = "squarem")
    expect_equal(dim(est), c(8L, 1L))
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    if (!ioo) expect_identical(est[1L], 0)
    r <- choicer:::mnl_blp_solve(rep(0, 8L), S, d$X, eng_beta, d$alt, d$M, d$w, ioo,
                                 eng_ctl(tol = 1e-12, method = "squarem", pass = "exp"))
    expect_identical(r$code, 0L)
    expect_lte(r$residual, 1e-12)
    expect_identical(r$delta, est)
    S <- nl_s(d, ds, ioo, nest, lam)
    est <- nl_blp_contraction(rep(0, 8L), S, d$X, eng_beta, lam, d$alt, nest, d$M,
                              d$w, ioo, tol = 1e-12, method = "squarem")
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    S <- as.numeric(choicer:::mxl_predict_shares(
      c(th, eng_free(ds, ioo)), d$X, d$W, d$alt, d$M, d$w, eta, c(0L, 1L), TRUE, TRUE,
      TRUE, ioo))
    est <- mxl_blp_contraction(rep(0, 8L), S, d$X, d$W, eng_beta, c(0.3, -0.2),
                               c(log(0.7), 0.2, log(0.5)), d$alt, d$M, d$w, eta,
                               c(0L, 1L), TRUE, TRUE, ioo, tol = 1e-12,
                               method = "squarem")
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
  }
})

test_that("SQUAREM's cycle: unit steps are the plain iteration, longer ones save evaluations", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- eng_design(1016L, 600L, 12L, 4L, K = 3L, border = 0.05)
  set.seed(61L)
  ds <- stats::rnorm(12L)
  ds <- ds - ds[1L]
  S <- mnl_s(d, ds, FALSE)
  run <- function(method, ...) {
    choicer:::mnl_blp_solve(rep(0, 12L), S, d$X, eng_beta, d$alt, d$M, d$w, FALSE,
                            eng_ctl(tol = 1e-10, max_iter = 5000L, method = method, ...))
  }
  # Bounds that stay at 1 (step_factor 1) make every extrapolation the double
  # step x2 = G(G(x0)): each turn is a damped base step, the plain iteration
  # without stalls, point for point
  r1 <- run("squarem", settings = list(step_factor = 1))
  rp <- run("simple", patience = 0L)
  expect_identical(r1[c("delta", "code", "evaluations", "iterations", "race")],
                   rp[c("delta", "code", "evaluations", "iterations", "race")])
  expect_identical(r1$restarts, 0L)
  # The default bounds grow while the extrapolations succeed: every scheme
  # converges, and in fewer evaluations than the plain iteration (698, 752
  # and 788 against 947 here; counts that another platform's rounding moves,
  # so not compared on CRAN)
  ev <- integer(3L)
  for (sc in 1:3) {
    r <- run("squarem", settings = list(scheme = sc))
    expect_identical(r$code, 0L, label = paste("scheme", sc))
    expect_equal(as.numeric(r$delta), ds, tolerance = 1e-8, label = paste("scheme", sc))
    ev[sc] <- r$evaluations
  }
  skip_on_cran()
  for (sc in 1:3) expect_lt(ev[sc], rp$evaluations, label = paste("scheme", sc))
})

test_that("the methods' settings are checked before any evaluation", {
  on.exit(set_omp_threads(2L), add = TRUE)
  d <- eng_design(1017L, 100L, 4L, 4L)
  S <- mnl_s(d, c(0, 0.2, 0.1, -0.3), FALSE)
  bc <- function(...) blp_contraction(rep(0, 4L), S, d$X, eng_beta, d$alt, d$M, d$w, ...)
  expect_error(bc(method = "newton"), "`method` must be \"anderson\" (the default",
               fixed = TRUE)
  expect_error(bc(method = NA_character_), "`method` must be")
  expect_error(bc(method = c("simple", "anderson")), "`method` must be")
  expect_error(bc(method = "simple", control = list(memory = 5L)),
               "`control` takes no settings with method \"simple\" (it has `memory`)",
               fixed = TRUE)
  # the default is Anderson acceleration, and takes its settings (at one
  # thread, where its path is the same from run to run)
  set_omp_threads(1L)
  expect_identical(bc(control = list(memory = 5L)),
                   bc(method = "anderson", control = list(memory = 5L)))
  expect_error(bc(method = "anderson", control = list(m = 5L)),
               "`control` takes `memory` and `max_condition` with method \"anderson\"; it has `m`",
               fixed = TRUE)
  expect_error(bc(method = "anderson", control = list(5L)), "must name each of its settings once")
  expect_error(bc(method = "anderson", control = list(memory = 5L, memory = 6L)),
               "must name each of its settings once")
  for (m in list(0L, 101L, 2.5, NA, "a", c(5L, 6L), NULL)) {
    expect_error(bc(method = "anderson", control = list(memory = m)),
                 "`control$memory` must be a whole number from 1 to 100", fixed = TRUE)
  }
  for (mc in list(0.5, NA_real_, "a", c(10, 100), NULL)) {
    expect_error(bc(method = "anderson", control = list(max_condition = mc)),
                 "`control$max_condition` must be a single number of at least 1",
                 fixed = TRUE)
  }
  expect_error(bc(method = "anderson", control = "memory"), "must be a list")
  expect_error(bc(method = "squarem", control = list(memory = 5L)),
               "`control` takes `scheme`, `step_min`, `step_max` and `step_factor` with method \"squarem\"; it has `memory`",
               fixed = TRUE)
  for (sc in list(0L, 4L, 2.5, NA, "3", NULL)) {
    expect_error(bc(method = "squarem", control = list(scheme = sc)),
                 "`control$scheme` must be 1, 2 or 3", fixed = TRUE)
  }
  for (b in list(list(step_min = 0), list(step_min = -1), list(step_min = 2),
                 list(step_max = Inf), list(step_max = NA_real_), list(step_min = "1"))) {
    expect_error(bc(method = "squarem", control = b),
                 "`control$step_min` and `control$step_max` must be finite numbers",
                 fixed = TRUE)
  }
  for (f in list(0.5, Inf, NA_real_, c(2, 4))) {
    expect_error(bc(method = "squarem", control = list(step_factor = f)),
                 "`control$step_factor` must be a finite number of at least 1",
                 fixed = TRUE)
  }
  expect_silent(bc(method = "squarem", control = list(scheme = 1L, step_max = 8)))
  ctl <- eng_ctl(method = "squarem", settings = list(step_factor = 2))
  expect_identical(ctl[c("iteration", "scheme", "step_min", "step_max", "step_factor")],
                   list(iteration = "squarem", scheme = 3L, step_min = 1, step_max = 1,
                        step_factor = 2))
  # the settings given, the others at their defaults
  expect_silent(e <- bc(method = "anderson", control = list(memory = 3L)))
  expect_silent(e <- bc(method = "anderson", control = list(max_condition = Inf)))
  ctl <- eng_ctl(method = "anderson", settings = list(memory = 3L))
  expect_identical(ctl[c("iteration", "memory", "max_condition", "patience", "window")],
                   list(iteration = "anderson", memory = 3L, max_condition = 1e8,
                        patience = 0L, window = 100L))
})

test_that("the exports keep their signatures and return n x 1 matrices", {
  # (method and control come last, after the arguments they always had)
  expect_identical(names(formals(blp_contraction)),
                   c("delta", "target_shares", "X", "beta", "alt_idx", "M", "weights",
                     "include_outside_option", "tol", "max_iter", "method", "control"))
  expect_identical(names(formals(nl_blp_contraction)),
                   c("delta", "target_shares", "X", "beta", "lambda", "alt_idx",
                     "nest_idx", "M", "weights", "include_outside_option", "damping",
                     "tol", "max_iter", "method", "control"))
  expect_identical(names(formals(mxl_blp_contraction)),
                   c("delta", "target_shares", "X", "W", "beta", "mu", "L_params",
                     "alt_idx", "M", "weights", "eta_draws", "rc_dist",
                     "rc_correlation", "rc_mean", "include_outside_option", "tol",
                     "max_iter", "gen_seed", "gen_scramble", "gen_S", "method",
                     "control"))
  for (f in list(blp_contraction, nl_blp_contraction, mxl_blp_contraction)) {
    expect_identical(formals(f)$method, "anderson")
    expect_identical(formals(f)$tol, 1e-12)
    expect_identical(formals(f)$control, quote(list()))
  }
  # blp()'s methods take the same defaults
  for (cls in c("choicer_mnl", "choicer_nl", "choicer_mxl")) {
    fm <- formals(utils::getS3method("blp", cls))
    expect_identical(fm$method, "anderson", label = cls)
    expect_identical(fm$tol, 1e-12, label = cls)
    expect_identical(fm$control, quote(list()), label = cls)
  }
  expect_identical(formals(blp_contraction)$max_iter, 1000L)
  expect_identical(formals(nl_blp_contraction)$damping, 1)
})
