# The share inversion's engine (src/blp_engine.h, src/fixed_point.h): its
# passes against the prediction kernels and a finite-difference Jacobian, the
# two-level map and its safeguards, components and the targets' checks,
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
                           tol = 1e-12)
    expect_equal(dim(est), c(8L, 1L))
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    if (!ioo) expect_identical(est[1L], 0)
    S <- nl_s(d, ds, ioo, nest, lam)
    est <- nl_blp_contraction(rep(0, 8L), S, d$X, eng_beta, lam, d$alt, nest, d$M,
                              d$w, ioo, tol = 1e-12)
    expect_equal(as.numeric(est), ds, tolerance = 1e-9)
    for (rc in list(c(0L, 0L), c(0L, 1L))) {
      th <- c(eng_beta, 0.3, -0.2, log(0.7), 0.2, log(0.5))
      S <- as.numeric(choicer:::mxl_predict_shares(
        c(th, eng_free(ds, ioo)), d$X, d$W, d$alt, d$M, d$w, eta, rc, TRUE, TRUE,
        TRUE, ioo))
      for (start in if (ioo) list(rep(0, 8L)) else list(rep(0, 8L), rep(0, 7L))) {
        est <- mxl_blp_contraction(start, S, d$X, d$W, eng_beta, c(0.3, -0.2),
                                   c(log(0.7), 0.2, log(0.5)), d$alt, d$M, d$w, eta,
                                   rc, TRUE, TRUE, ioo, tol = 1e-12)
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
  # the diagonal step alone never improves on its start; the restarts halve
  # it until it contracts
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(map = "jacobi", patience = 0L))
  expect_identical(r$code, 1L)
  expect_identical(as.numeric(r$delta), rep(0, 4L))
  r <- choicer:::mnl_blp_solve(rep(0, 4L), S, X, 1, alt, M, w, FALSE,
                               eng_ctl(map = "jacobi", patience = 30L))
  expect_identical(r$code, 0L)
  expect_gte(r$restarts, 1L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-8)
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
  # through few situations; about 3,050 evaluations here, Berry's 5,500)
  d <- eng_design(1007L, 1200L, 24L, 4L, K = 4L, border = 0.02)
  set.seed(41L)
  ds <- stats::rnorm(24L)
  ds <- ds - ds[1L]
  S <- mnl_s(d, ds, FALSE)
  set.seed(42L)
  start <- ds + stats::rnorm(24L, 0, 3)
  r <- choicer:::mnl_blp_solve(start, S, d$X, eng_beta, d$alt, d$M, d$w, FALSE,
                               eng_ctl(tol = 1e-10, max_iter = 4000L))
  expect_identical(r$code, 0L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-6)
})

test_that("an outside share that underflows takes the step without its correction", {
  d <- eng_design(1008L, 200L, 5L, 5L)
  ds <- c(0.3, -0.2, 0.1, 0.4, -0.5)
  S <- mnl_s(d, ds, TRUE)
  r <- choicer:::mnl_blp_solve(rep(800, 5L), S, d$X, eng_beta, d$alt, d$M, d$w,
                               TRUE, eng_ctl(tol = 1e-10))
  expect_identical(r$code, 0L)
  expect_gt(r$fallbacks, 0L)
  expect_equal(as.numeric(r$delta), ds, tolerance = 1e-8)
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
  set_omp_threads(1L)
  e1 <- blp_contraction(rep(0, 10L), S, d$X, eng_beta, d$alt, d$M, d$w, tol = 1e-12)
  set_omp_threads(2L)
  e2 <- blp_contraction(rep(0, 10L), S, d$X, eng_beta, d$alt, d$M, d$w, tol = 1e-12)
  expect_equal(e2, e1, tolerance = 1e-10)
})

test_that("the exports keep their signatures and return n x 1 matrices", {
  expect_identical(names(formals(blp_contraction)),
                   c("delta", "target_shares", "X", "beta", "alt_idx", "M", "weights",
                     "include_outside_option", "tol", "max_iter"))
  expect_identical(names(formals(nl_blp_contraction)),
                   c("delta", "target_shares", "X", "beta", "lambda", "alt_idx",
                     "nest_idx", "M", "weights", "include_outside_option", "damping",
                     "tol", "max_iter"))
  expect_identical(names(formals(mxl_blp_contraction)),
                   c("delta", "target_shares", "X", "W", "beta", "mu", "L_params",
                     "alt_idx", "M", "weights", "eta_draws", "rc_dist",
                     "rc_correlation", "rc_mean", "include_outside_option", "tol",
                     "max_iter", "gen_seed", "gen_scramble", "gen_S"))
  expect_identical(formals(blp_contraction)$max_iter, 1000L)
  expect_identical(formals(nl_blp_contraction)$damping, 1)
})
