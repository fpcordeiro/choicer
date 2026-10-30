# The share inversion's exp-space pass (src/blp_engine.h, BlpExpCache): an
# evaluation from exponentials kept at an origin against the prediction kernels
# and the usual (log-space) pass, the span and the floor that send an
# evaluation back to the usual pass, and the cache's budget.

xp_design <- function(seed, N, J, m, Kx = 2L) {
  set.seed(seed)
  sets <- lapply(seq_len(N), function(i) sort(sample.int(J, m)))
  alt <- as.integer(unlist(sets))
  n <- length(alt)
  list(X = matrix(stats::rnorm(n * Kx), n, Kx), alt = alt, M = lengths(sets),
       w = stats::runif(N, 0.5, 2), J = J)
}
xp_beta <- c(0.8, -0.5)
xp_norm <- function(d, ioo) if (ioo) d else d - d[1L]
xp_free <- function(d, ioo) if (ioo) d else (d - d[1L])[-1L]
xp_ctl <- function(...) choicer:::.blp_control(...)
xp_pass <- function(r) lapply(r$pass, as.numeric)

xp_mnl_s <- function(d, delta, ioo) {
  as.numeric(choicer:::mnl_predict_shares(c(xp_beta, xp_free(delta, ioo)), d$X,
                                          d$alt, d$M, d$w, TRUE, ioo))
}
xp_mnl <- function(d, start, S, ioo, ctl) {
  choicer:::mnl_blp_solve(start, S, d$X, xp_beta, d$alt, d$M, d$w, ioo, ctl)
}
# Targets every alternative can meet: the shares at some constants
xp_target <- function(d, ioo, seed) {
  set.seed(seed)
  xp_mnl_s(d, stats::rnorm(d$J), ioo)
}
# The floor F(d) of an evaluation at distance d from the cache's origin
xp_floor <- function(dist) 2^53 * .Machine$double.xmin * exp(2 * dist)

test_that("an evaluation from the cache gives the prediction kernel's shares", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2001L, 300L, 7L, 4L)
  # a code listed twice in some situations
  r0 <- c(0L, cumsum(d$M))
  for (t in 1:12) d$alt[r0[t] + 1L] <- d$alt[r0[t] + 2L]
  for (ioo in c(FALSE, TRUE)) {
    set.seed(5L + ioo)
    ds <- xp_norm(stats::rnorm(7L), ioo)
    d1 <- xp_norm(ds + stats::rnorm(7L, 0, 0.5), ioo)
    S <- xp_target(d, ioo, 7L)
    lab <- sprintf("ioo %d", ioo)
    rx <- xp_mnl(d, d1, S, ioo, xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp",
                                        origin = ds))
    rl <- xp_mnl(d, d1, S, ioo, xp_ctl(max_iter = 1L, keep_pass = TRUE))
    expect_true(rx$exp_kept, label = lab)
    expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(1, 1, 0),
                     label = lab)
    s1 <- xp_mnl_s(d, d1, ioo)
    px <- xp_pass(rx)
    expect_equal(px$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
    if (ioo) expect_equal(px$s0, s1[1L], tolerance = 1e-14, label = lab)
    # the diagonal and the co-choice as the usual pass forms them
    pl <- xp_pass(rl)
    expect_equal(px$A, pl$A, tolerance = 1e-13, label = lab)
    expect_equal(px$C, pl$C, tolerance = 1e-13, label = lab)
    expect_equal(rx$residual, rl$residual, tolerance = 1e-12, label = lab)
    # the usual pass reports no cache
    expect_null(rl$exp_evals)
  }
})

test_that("the exp-space pass solves as the usual pass does", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2002L, 600L, 9L, 5L)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(11L + ioo)
    ds <- xp_norm(stats::rnorm(9L, 0, 1.5), ioo)
    S <- xp_mnl_s(d, ds, ioo)
    for (method in c("simple", "original", "anderson", "squarem")) {
      lab <- sprintf("%s, ioo %d", method, ioo)
      rl <- xp_mnl(d, rep(0, 9L), S, ioo, xp_ctl(1e-12, 1000L, method = method))
      rx <- xp_mnl(d, rep(0, 9L), S, ioo, xp_ctl(1e-12, 1000L, method = method,
                                                   pass = "exp"))
      expect_identical(rx$code, 0L, label = lab)
      expect_identical(rx$code, rl$code, label = lab)
      expect_equal(as.numeric(rx$delta), ds, tolerance = 1e-10, label = lab)
      expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-11,
                   label = lab)
      # every evaluation from the cache, but the first, which forms it
      expect_identical(rx$exp_forms, 1, label = lab)
      expect_identical(rx$exp_evals, rx$evaluations - 1, label = lab)
    }
  }
  skip_if_not(isTRUE(thread_info()$openmp_enabled))
  set_omp_threads(2L)
  S <- xp_mnl_s(d, ds, TRUE)
  r1 <- xp_mnl(d, rep(0, 9L), S, TRUE, xp_ctl(1e-12, 1000L, method = "anderson",
                                               pass = "exp"))
  expect_equal(as.numeric(r1$delta), ds, tolerance = 1e-10)
})

test_that("an evaluation beyond the span is the usual pass, and moves the origin", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2003L, 200L, 5L, 3L)
  S <- xp_target(d, TRUE, 13L)
  ds <- c(0.3, -0.2, 0.5, 0.1, -0.4)
  for (far in c(49, 51)) {
    d1 <- replace(ds, 3L, ds[3L] + far)
    rx <- xp_mnl(d, d1, S, TRUE, xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp",
                                         origin = ds))
    rl <- xp_mnl(d, d1, S, TRUE, xp_ctl(max_iter = 1L, keep_pass = TRUE))
    if (far < 50) {
      expect_identical(c(rx$exp_evals, rx$exp_forms), c(1, 1))
      expect_equal(xp_pass(rx)$s, xp_pass(rl)$s, tolerance = 1e-14)
    } else {
      expect_identical(c(rx$exp_evals, rx$exp_forms), c(0, 2))
      expect_identical(xp_pass(rx), xp_pass(rl))  # the usual pass, bit for bit
    }
  }
})

test_that("an evaluation that gives a share below its floor is the usual pass", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2004L, 200L, 4L, 4L)
  S <- c(0.3, 0.3, 0.3, 0.1)
  # the origin's fourth share far below the others', but above F(50); 40 log
  # units lower it falls below F(40)
  ds <- c(0, 0.2, -0.1, -560)
  d1 <- replace(ds, 4L, -600)
  s0 <- xp_mnl_s(d, ds, FALSE)
  s1 <- xp_mnl_s(d, d1, FALSE)
  expect_gt(s0[4L], xp_floor(50))
  expect_lt(s1[4L], xp_floor(40))
  rx <- xp_mnl(d, d1, S, FALSE, xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp",
                                        origin = ds))
  rl <- xp_mnl(d, d1, S, FALSE, xp_ctl(max_iter = 1L, keep_pass = TRUE))
  expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(0, 2, 1))
  expect_identical(xp_pass(rx), xp_pass(rl))
  # A cache whose own pass gave a share below F(50) is not used: here the
  # fourth share underflows at the origin
  d2 <- replace(ds, 4L, -800)
  rx <- xp_mnl(d, replace(d2, 2L, 0.3), S, FALSE,
               xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = d2))
  expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(0, 2, 0))
  # From such a start the solver climbs out and then uses the cache
  start <- replace(rep(0, 4L), 4L, -700)
  S <- xp_mnl_s(d, c(0, 0.2, -0.1, -3), FALSE)
  rl <- xp_mnl(d, start, S, FALSE, xp_ctl(1e-12, 1000L, method = "anderson"))
  rx <- xp_mnl(d, start, S, FALSE, xp_ctl(1e-12, 1000L, method = "anderson",
                                           pass = "exp"))
  expect_identical(rx$code, rl$code)
  expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-10)
  expect_gt(rx$exp_forms, 1)
  expect_gt(rx$exp_evals, 0)
  # An evaluation accepted near its floor (a share near 1e-250 at distance
  # 45) is the usual pass's to rounding, share by share
  ds <- c(0, 0.2, -0.1, -530)
  d1 <- replace(ds, 4L, -575)
  expect_gt(xp_mnl_s(d, d1, FALSE)[4L], 10 * xp_floor(45))
  rx <- xp_mnl(d, d1, S, FALSE, xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp",
                                        origin = ds))
  rl <- xp_mnl(d, d1, S, FALSE, xp_ctl(max_iter = 1L, keep_pass = TRUE))
  expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(1, 1, 0))
  expect_lt(max(abs(log(xp_pass(rx)$s) - log(xp_pass(rl)$s))), 1e-13)
})

test_that("an outside share below its floor sends the evaluation to the usual pass", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2006L, 200L, 4L, 3L)
  S <- c(0.1, rep(0.225, 4L))
  # every inside constant 560 up: the outside share near e^-560, above F(50);
  # 40 more and it falls below F(40)
  ds <- rep(560, 4L)
  d1 <- rep(600, 4L)
  expect_gt(xp_mnl_s(d, ds, TRUE)[1L], xp_floor(50))
  expect_lt(xp_mnl_s(d, d1, TRUE)[1L], xp_floor(40))
  rx <- xp_mnl(d, d1, S, TRUE, xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp",
                                       origin = ds))
  rl <- xp_mnl(d, d1, S, TRUE, xp_ctl(max_iter = 1L, keep_pass = TRUE))
  expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(0, 2, 1))
  expect_identical(xp_pass(rx), xp_pass(rl))
  # an origin whose outside share underflows: the cache waits
  d2 <- rep(800, 4L)
  rx <- xp_mnl(d, d2 + c(0, 0.1, 0, 0), S, TRUE,
               xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = d2))
  expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(0, 2, 0))
})

test_that("the cache keeps each component's reference and outside share", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # two regions of three alternatives, no situation offering both
  set.seed(2007L)
  sets <- lapply(seq_len(400L), function(i) sort(3L * (i %% 2L) + sample.int(3L, 2L)))
  alt <- as.integer(unlist(sets))
  d <- list(X = matrix(stats::rnorm(2L * length(alt)), ncol = 2L), alt = alt,
            M = lengths(sets), w = stats::runif(400L, 0.5, 2), J = 6L)
  ds <- c(0, 0.4, -0.3, 0, 0.2, -0.5)  # (the references, 1 and 4, at 0)
  for (ioo in c(FALSE, TRUE)) {
    S <- xp_mnl_s(d, ds, ioo)
    for (method in c("simple", "anderson")) {
      lab <- sprintf("%s, ioo %d", method, ioo)
      rl <- xp_mnl(d, rep(0, 6L), S, ioo, xp_ctl(1e-12, 1000L, method = method))
      rx <- xp_mnl(d, rep(0, 6L), S, ioo, xp_ctl(1e-12, 1000L, method = method,
                                                  pass = "exp"))
      expect_identical(rx$code, 0L, label = lab)
      expect_equal(as.numeric(rx$delta), ds, tolerance = 1e-10, label = lab)
      expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-11,
                   label = lab)
      expect_gt(rx$exp_evals, 0)
    }
  }
})

test_that("a test origin must hold finite constants", {
  d <- xp_design(2008L, 50L, 4L, 3L)
  S <- xp_target(d, FALSE, 19L)
  for (o in list(c(0, 1), c(0, NA, 0, 0))) {
    expect_error(xp_mnl(d, rep(0, 4L), S, FALSE, xp_ctl(max_iter = 2L, pass = "exp",
                                                         origin = o)),
                 "origin must hold 4 finite constants")
  }
})

test_that("the cache is kept only within its budget", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2005L, 150L, 5L, 3L)
  n <- length(d$alt)
  N <- length(d$M)
  for (ioo in c(FALSE, TRUE)) {
    S <- xp_target(d, ioo, 17L)
    need <- 8 * (n + if (ioo) N else 0)
    rl <- xp_mnl(d, rep(0, 5L), S, ioo, xp_ctl(1e-10, 1000L))
    for (b in list(list(0, FALSE), list(NaN, FALSE), list(need - 1, FALSE),
                   list(need, TRUE), list(Inf, TRUE))) {
      lab <- sprintf("ioo %d, budget %s", ioo, format(b[[1]]))
      rx <- xp_mnl(d, rep(0, 5L), S, ioo, xp_ctl(1e-10, 1000L, pass = "exp",
                                                  cache_bytes = b[[1]]))
      expect_identical(rx$exp_kept, b[[2]], label = lab)
      if (!b[[2]]) {
        expect_identical(c(rx$exp_evals, rx$exp_forms), c(0, 0), label = lab)
        expect_identical(rx$delta, rl$delta, label = lab)
      }
    }
    # capped by half the memory limit, less what the call holds
    old <- options(choicer.max_memory = 1e4)
    on.exit(options(old), add = TRUE)
    rx <- xp_mnl(d, rep(0, 5L), S, ioo, xp_ctl(1e-10, 1000L, pass = "exp"))
    expect_false(rx$exp_kept)
    options(choicer.max_memory = Inf)
    rx <- xp_mnl(d, rep(0, 5L), S, ioo, xp_ctl(1e-10, 1000L, pass = "exp"))
    expect_true(rx$exp_kept)
    options(old)
    # one evaluation keeps nothing (no evaluation could use it)
    rx <- xp_mnl(d, rep(0, 5L), S, ioo, xp_ctl(1e-10, 1L, pass = "exp"))
    expect_false(rx$exp_kept)
  }
})

# -- The nested logit ---------------------------------------------------------

xp_nl_s <- function(d, delta, ioo, nest, lam) {
  as.numeric(choicer:::nl_predict_shares(
    c(xp_beta, lam[tabulate(nest) > 1L], xp_free(delta, ioo)), d$X, d$alt, d$M,
    d$w, nest, TRUE, ioo))
}
xp_nl <- function(d, start, S, ioo, nest, lam, ctl) {
  choicer:::nl_blp_solve(start, S, d$X, xp_beta, lam, d$alt, nest, d$M, d$w, ioo, ctl)
}

test_that("a nested logit's evaluation from the cache gives the prediction kernel's shares", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # three nests of three and a singleton; repeated codes in some situations
  d <- xp_design(2101L, 300L, 10L, 6L)
  nest <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L)
  r0 <- c(0L, cumsum(d$M))
  for (t in 1:12) d$alt[r0[t] + 1L] <- d$alt[r0[t] + 2L]
  for (lam_in in c(1e-3, 0.4, 1.4)) {
    lam <- c(lam_in, 0.7, lam_in, 1)
    for (ioo in c(FALSE, TRUE)) {
      set.seed(23L + ioo)
      ds <- xp_norm(stats::rnorm(10L), ioo)
      # well within the span, in units of lambda
      d1 <- xp_norm(ds + stats::rnorm(10L, 0, 0.1 * min(lam)), ioo)
      S <- xp_nl_s(d, xp_norm(stats::rnorm(10L), ioo), ioo, nest, lam)
      lab <- sprintf("lambda %g, ioo %d", lam_in, ioo)
      rx <- xp_nl(d, d1, S, ioo, nest, lam,
                  xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
      rl <- xp_nl(d, d1, S, ioo, nest, lam, xp_ctl(max_iter = 1L, keep_pass = TRUE))
      expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(1, 1, 0),
                       label = lab)
      s1 <- xp_nl_s(d, d1, ioo, nest, lam)
      px <- xp_pass(rx)
      pl <- xp_pass(rl)
      expect_equal(px$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
      if (ioo) expect_equal(px$s0, s1[1L], tolerance = 1e-14, label = lab)
      expect_equal(px$A, pl$A, tolerance = 1e-12, label = lab)
      expect_equal(px$C, pl$C, tolerance = 1e-12, label = lab)
      if (lam_in > 1) {
        expect_equal(px$A1, pl$A1, tolerance = 1e-12, label = lab)
        expect_equal(px$C1, pl$C1, tolerance = 1e-12, label = lab)
      }
    }
  }
})

test_that("a nested logit's span is in units of its dissimilarities", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2102L, 200L, 6L, 4L)
  nest <- c(1L, 1L, 1L, 2L, 2L, 2L)
  lam <- c(0.5, 1)
  S <- xp_nl_s(d, c(0, 0.2, -0.1, 0.3, 0, -0.2), FALSE, nest, lam)
  ds <- c(0, 0.2, -0.1, 0.3, 0, -0.2)
  for (move in c(24, 26)) {  # scaled: 48 and 52
    d1 <- replace(ds, 2L, ds[2L] + move)
    rx <- xp_nl(d, d1, S, FALSE, nest, lam,
                xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
    rl <- xp_nl(d, d1, S, FALSE, nest, lam, xp_ctl(max_iter = 1L, keep_pass = TRUE))
    if (move < 25) {
      expect_identical(c(rx$exp_evals, rx$exp_forms), c(1, 1))
      expect_equal(xp_pass(rx)$s, xp_pass(rl)$s, tolerance = 1e-14)
    } else {
      expect_identical(c(rx$exp_evals, rx$exp_forms), c(0, 2))
      expect_identical(xp_pass(rx), xp_pass(rl))
    }
  }
  # an alternative of the lambda = 1 nest moved by 40 is within the span (in
  # units of 0.5 it would not be)
  d1 <- replace(ds, 5L, ds[5L] + 40)
  rx <- xp_nl(d, d1, S, FALSE, nest, lam,
              xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
  rl <- xp_nl(d, d1, S, FALSE, nest, lam, xp_ctl(max_iter = 1L, keep_pass = TRUE))
  expect_identical(c(rx$exp_evals, rx$exp_forms), c(1, 1))
  expect_equal(xp_pass(rx)$s, xp_pass(rl)$s, tolerance = 1e-14)
})

test_that("rows negligible at the cache's origin are kept for when they dominate", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # each alternative far below its nest-mates at the origin (scaled distances
  # up to 40), then moved within the span until it dominates its nest
  d <- xp_design(2105L, 200L, 6L, 4L)
  nest <- c(1L, 1L, 1L, 2L, 2L, 2L)
  for (ioo in c(FALSE, TRUE)) {
    for (lam in c(0.3, 1)) {
      ds <- xp_norm(c(0, -8, 0.5, 0, -12 * lam, 0.3) , ioo)
      d1 <- xp_norm(replace(ds, c(2L, 5L), c(4, 4 * lam)), ioo)
      S <- xp_nl_s(d, ds, ioo, nest, rep(lam, 2L))
      lab <- sprintf("lambda %g, ioo %d", lam, ioo)
      rx <- xp_nl(d, d1, S, ioo, nest, rep(lam, 2L),
                  xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
      expect_identical(c(rx$exp_evals, rx$exp_forms), c(1, 1), label = lab)
      s1 <- xp_nl_s(d, d1, ioo, nest, rep(lam, 2L))
      expect_equal(xp_pass(rx)$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14,
                   label = lab)
      expect_lt(max(abs(log(xp_pass(rx)$s) - log(if (ioo) s1[-1L] else s1))), 1e-12,
                label = lab)
      # the same for the multinomial logit
      rx <- xp_mnl(d, d1, S[seq_along(S)], ioo,
                   xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
      s1 <- xp_mnl_s(d, d1, ioo)
      expect_lt(max(abs(log(xp_pass(rx)$s) - log(if (ioo) s1[-1L] else s1))), 1e-12,
                label = lab)
    }
  }
})

test_that("the nested logit's exp-space pass solves as its usual pass does", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2103L, 600L, 9L, 6L)
  nest <- rep(1:3, each = 3L)
  for (lam in c(0.2, 0.7, 1.5)) {
    for (ioo in c(FALSE, TRUE)) {
      set.seed(29L + ioo)
      ds <- xp_norm(stats::rnorm(9L), ioo)
      S <- xp_nl_s(d, ds, ioo, nest, rep(lam, 3L))
      for (method in c("simple", "anderson")) {
        lab <- sprintf("%s, lambda %g, ioo %d", method, lam, ioo)
        rl <- xp_nl(d, rep(0, 9L), S, ioo, nest, rep(lam, 3L),
                    xp_ctl(1e-12, 2000L, method = method))
        rx <- xp_nl(d, rep(0, 9L), S, ioo, nest, rep(lam, 3L),
                    xp_ctl(1e-12, 2000L, method = method, pass = "exp"))
        expect_identical(rx$code, 0L, label = lab)
        expect_equal(as.numeric(rx$delta), ds, tolerance = 1e-9, label = lab)
        expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-10,
                     label = lab)
        expect_identical(rx$exp_evals + rx$exp_forms, as.numeric(rx$evaluations),
                         label = lab)
        expect_gt(rx$exp_evals, 0)
      }
    }
  }
  skip_if_not(isTRUE(thread_info()$openmp_enabled))
  set_omp_threads(2L)
  S <- xp_nl_s(d, ds, TRUE, nest, rep(0.4, 3L))
  r2 <- xp_nl(d, rep(0, 9L), S, TRUE, nest, rep(0.4, 3L),
              xp_ctl(1e-12, 2000L, method = "anderson", pass = "exp"))
  expect_equal(as.numeric(r2$delta), ds, tolerance = 1e-9)
  expect_gt(r2$exp_evals, 0)
})

test_that("a nested logit with more nests than rows keeps its cache on the offered path", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # six nests of two, three rows per situation: the probabilities go over the
  # nests each situation offers (nl_individual_probs_offered())
  d <- xp_design(2104L, 300L, 12L, 3L)
  nest <- rep(1:6, each = 2L)
  lam <- c(0.3, 0.6, 1, 1.4, 0.8, 0.5)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(31L + ioo)
    ds <- xp_norm(stats::rnorm(12L), ioo)
    S <- xp_nl_s(d, xp_norm(stats::rnorm(12L), ioo), ioo, nest, lam)
    d1 <- xp_norm(ds + stats::rnorm(12L, 0, 0.05), ioo)
    lab <- sprintf("ioo %d", ioo)
    rx <- xp_nl(d, d1, S, ioo, nest, lam,
                xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
    expect_identical(c(rx$exp_evals, rx$exp_forms), c(1, 1), label = lab)
    s1 <- xp_nl_s(d, d1, ioo, nest, lam)
    expect_equal(xp_pass(rx)$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
    # beyond the span (30 utility units in a nest with lambda 0.3): the usual
    # pass, which forms the cache, bit for bit
    d2 <- replace(ds, 2L, ds[2L] + 30)
    rx <- xp_nl(d, d2, S, ioo, nest, lam,
                xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds))
    rl <- xp_nl(d, d2, S, ioo, nest, lam, xp_ctl(max_iter = 1L, keep_pass = TRUE))
    expect_identical(c(rx$exp_evals, rx$exp_forms), c(0, 2), label = lab)
    expect_identical(xp_pass(rx), xp_pass(rl), label = lab)
    # and the cache is kept only within its budget
    rx <- xp_nl(d, rep(0, 12L), S, ioo, nest, lam,
                xp_ctl(1e-10, 1000L, pass = "exp", cache_bytes = 8 * length(d$alt) - 1))
    expect_false(rx$exp_kept, label = lab)
    rl <- xp_nl(d, rep(0, 12L), S, ioo, nest, lam, xp_ctl(1e-10, 1000L))
    expect_identical(rx$delta, rl$delta, label = lab)
    rx <- xp_nl(d, rep(0, 12L), S, ioo, nest, lam,
                xp_ctl(1e-10, 1000L, pass = "exp", cache_bytes = 8 * length(d$alt)))
    expect_true(rx$exp_kept, label = lab)
  }
})

# -- Compensated share sums ---------------------------------------------------

test_that("compensated share sums keep the small terms a plain sum drops", {
  skip_on_cran()  # (two million rows)
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # one situation of weight 1, where alternative 1's probability is p (in
  # [2^-8, 2^-7)), then 2^20 of weight 2^-62 where it is 1/2: each small term,
  # 2^-63, is below half a unit of the running share, so a plain sum drops
  # every one; a plain total of the weights drops half of them (one of its two
  # accumulators starts at 1). The exact share is (p + 2^-43) / (1 + 2^-42),
  # p + 2^-43 and 1 + 2^-42 being exact.
  n_small <- 2^20
  N <- n_small + 1
  alt <- rep(1:2, N)
  M <- rep(2L, N)
  X <- matrix(0, 2 * N, 1L)
  X[2L, 1L] <- 5
  w <- c(1, rep(2^-62, n_small))
  run <- function(comp, X, alt, M, w) {
    choicer:::mnl_blp_solve(c(0, 0), c(0.5, 0.5), X, 1, alt, M, w, FALSE,
                            xp_ctl(max_iter = 1L, keep_pass = TRUE, compensate = comp))
  }
  # p as the pass forms it, from the first situation alone
  p <- xp_pass(run("off", X[1:2, , drop = FALSE], 1:2, 2L, 1))$s[1L]
  ref <- (p + 2^-43) / (1 + 2^-42)
  rk <- run("on", X, alt, M, w)
  rp <- run("off", X, alt, M, w)
  expect_true(rk$compensated)
  expect_null(rp$compensated)
  expect_lt(abs(xp_pass(rk)$s[1L] / ref - 1), 2^-52)
  expect_gt(abs(xp_pass(rp)$s[1L] / ref - 1), 1e-14)
  # "auto" compensates only above 2^22 situations (per thread for the passes)
  ra <- run("auto", X, alt, M, w)
  expect_false(ra$compensated)
  expect_identical(xp_pass(ra), xp_pass(rp))
})

test_that("\"auto\" compensates the sums above 2^22 situations per thread", {
  skip_on_cran()  # (eight million rows)
  on.exit(set_omp_threads(2L), add = TRUE)
  N <- 2^22 + 1
  alt <- rep(1:2, N)
  M <- rep(2L, N)
  X <- matrix(0, length(alt), 1L)
  w <- rep(1, N)
  run <- function(comp) {
    choicer:::mnl_blp_solve(c(0, 0), c(0.5, 0.5), X, 1, alt, M, w, FALSE,
                            xp_ctl(max_iter = 1L, compensate = comp))
  }
  set_omp_threads(1L)
  expect_true(run("auto")$compensated)
  ti <- thread_info()
  skip_if_not(isTRUE(ti$openmp_enabled) && isTRUE(ti$omp_get_thread_limit >= 2L))
  set_omp_threads(2L)
  expect_false(run("auto")$compensated)  # (2^21 + 1 situations per thread)
})

test_that("compensated sums bound each target by its group's weight", {
  # One situation of weight 1, then 2^16 of weight 2^-56, every one offering
  # alternatives 1 and 2: plain sums lose the small weights (1 + 2^-56 rounds
  # to 1) and compensated sums keep them (1 + 2^-40). The targets are
  # attained at the start; an offer weight summed plainly over the
  # compensated total (1 - 2^-40 of it) fell below alternative 1's (PR #58's
  # review). The offer weights now share the total's compensation and order,
  # so alternative 1's bound is its group's weight exactly, whatever the
  # weights' common scale.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  n <- 2^16 + 1
  alt <- rep(1:2, n)
  M <- rep(2L, n)
  X <- matrix(0, 2 * n, 1L)
  S <- c(1 - 5e-13, 5e-13)
  d0 <- c(0, log(S[2] / S[1]))
  ctl <- xp_ctl(1e-8, 1L, compensate = "on")
  for (s in c(1, 1.3)) {
    w <- s * c(1, rep(2^-56, n - 1))
    lab <- sprintf("scale %g", s)
    r <- choicer:::mnl_blp_solve(d0, S, X, 0, alt, M, w, FALSE, ctl)
    expect_identical(r$code, 0L, label = lab)
    r <- choicer:::nl_blp_solve(d0, S, X, 0, 1, alt, c(1L, 1L), M, w, FALSE, ctl)
    expect_identical(r$code, 0L, label = lab)
    r <- choicer:::mxl_blp_solve(d0, S, X, matrix(0, 2 * n, 1L), 0, 0, 0, alt, M, w,
                                 array(0, c(1L, 0L, 0L)), 0L, FALSE, FALSE, FALSE, ctl,
                                 1L, 1L, 2L, 0)
    expect_identical(r$code, 0L, label = lab)
    # with an outside option, its share and alternative 2's below 2^-40
    So <- c(3e-13, 1 - 6e-13, 3e-13)
    r <- choicer:::mnl_blp_solve(log(So[-1L] / So[1L]), So, X, 0, alt, M, w, TRUE, ctl)
    expect_identical(r$code, 0L, label = lab)
  }
  # and a second group, alternatives 3 and 4 in one situation of weight 1:
  # alternative 1's bound is its group's weight over the total
  alt2 <- c(alt, 3:4)
  w2 <- c(1, rep(2^-56, n - 1), 1)
  mA <- (1 + 2^-40) / (2 + 2^-40)
  S2 <- c(mA * S, (1 - mA) / 2, (1 - mA) / 2)
  r <- choicer:::mnl_blp_solve(c(d0, 0, 0), S2, rbind(X, 0, 0), 0, alt2, c(M, 2L), w2,
                               FALSE, ctl)
  expect_identical(r$code, 0L)
})

test_that("blp() and the exports accept targets near a bound above 2^22 situations", {
  skip_on_cran()  # (eight million rows)
  # PR #58's review: alternatives 1 and 2 offered in every one of 2^22 + 1
  # situations, so each offer bound is 1 exactly; the targets are attained
  # at the start. With weights 1.3 the plain offer weight over the
  # compensated total was 0.99999999992, below alternative 1's target.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  N <- 2^22 + 1
  S <- c(1 - 7e-11, 7e-11)
  X <- matrix(0, 2 * N, 1L)
  alt <- rep(1:2, N)
  M <- rep(2L, N)
  d0 <- c(0, log(S[2] / S[1]))
  for (wv in c(1.3, 1)) {
    w <- rep(wv, N)
    lab <- sprintf("weights %g", wv)
    expect_equal(as.numeric(blp_contraction(d0, S, X, 0, alt, M, w, tol = 1e-8,
                                            max_iter = 1L)), d0, label = lab)
    expect_equal(as.numeric(nl_blp_contraction(d0, S, X, 0, 1, alt, c(1L, 1L), M, w,
                                               tol = 1e-8, max_iter = 1L)), d0,
                 label = lab)
    expect_equal(as.numeric(mxl_blp_contraction(d0, S, X, matrix(0, 2 * N, 1L), 0, 0, 0,
                                                alt, M, w, array(0, c(1L, 0L, 0L)), 0L,
                                                FALSE, FALSE, tol = 1e-8, max_iter = 1L,
                                                gen_seed = 1L, gen_scramble = 1L,
                                                gen_S = 2L)), d0, label = lab)
  }
})

test_that("compensated sums keep the components' weights exact", {
  skip_on_cran()  # (two million rows)
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  # an outside option and two components: alternatives 1 and 2 in one
  # situation of weight 1 and 2^20 of weight 2^-62, alternatives 3 and 4 in 10
  # of weight 1. Without covariates every situation of a component has the
  # same probabilities, so the shares are each component's probabilities
  # times its weight, 1 + 2^-42 or 10, over 11 + 2^-42 (all exact), and the
  # outside option's share in a component (its weight less its inside
  # targets) needs that weight exactly.
  n_small <- 2^20
  alt <- c(rep(1:2, 1 + n_small), rep(3:4, 10))
  M <- rep(2L, 1 + n_small + 10)
  X <- matrix(0, length(alt), 1L)
  w <- c(1, rep(2^-62, n_small), rep(1, 10))
  ds <- c(0.5, -0.3, 0.2, 0.4)
  one <- function(d) {  # a component's probabilities, the outside option's first
    as.numeric(choicer:::mnl_predict_shares(c(0, d), matrix(0, 2, 1), 1:2, 2L, 1, TRUE, TRUE))
  }
  P1 <- one(ds[1:2])
  P2 <- one(ds[3:4])
  w1 <- 1 + 2^-42
  W <- 11 + 2^-42
  S <- c(P1[1L] * w1 + P2[1L] * 10, P1[2:3] * w1, P2[2:3] * 10) / W
  run <- function(comp) {
    choicer:::mnl_blp_solve(ds, S, X, 0, alt, M, w, TRUE,
                            xp_ctl(0, 1L, compensate = comp))
  }
  # (a component's outside target, its weight less its inside targets, loses
  # a few units to cancellation)
  expect_lt(run("on")$residual, 4e-15)
  expect_gt(run("off")$residual, 1e-13)
})

test_that("compensated sums give each model's prediction to rounding", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2201L, 300L, 8L, 5L)
  nest <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  lam <- c(0.5, 0.8, 1, 1.3)
  for (ioo in c(FALSE, TRUE)) {
    set.seed(41L + ioo)
    ds <- xp_norm(stats::rnorm(8L), ioo)
    d1 <- xp_norm(ds + stats::rnorm(8L, 0, 0.1), ioo)
    S <- xp_mnl_s(d, ds, ioo)
    S2 <- xp_nl_s(d, ds, ioo, nest, lam)
    s1 <- xp_mnl_s(d, d1, ioo)
    s2 <- xp_nl_s(d, d1, ioo, nest, lam)
    for (pass in c("log", "exp")) {
      lab <- sprintf("%s, ioo %d", pass, ioo)
      # (the exp pass from a cache formed at ds: one evaluation from it)
      ctl <- xp_ctl(max_iter = 1L, keep_pass = TRUE, compensate = "on", pass = pass,
                    origin = if (pass == "exp") ds)
      r <- xp_mnl(d, d1, S, ioo, ctl)
      expect_true(r$compensated, label = lab)
      if (pass == "exp") expect_identical(r$exp_evals, 1, label = lab)
      expect_equal(xp_pass(r)$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
      if (ioo) expect_equal(xp_pass(r)$s0, s1[1L], tolerance = 1e-14, label = lab)
      r <- xp_nl(d, d1, S2, ioo, nest, lam, ctl)
      expect_true(r$compensated, label = lab)
      if (pass == "exp") expect_identical(r$exp_evals, 1, label = lab)
      expect_equal(xp_pass(r)$s, if (ioo) s2[-1L] else s2, tolerance = 1e-14, label = lab)
      if (ioo) expect_equal(xp_pass(r)$s0, s2[1L], tolerance = 1e-14, label = lab)
    }
    # and they solve as the plain sums do, on either pass
    for (method in c("simple", "anderson")) {
      lab <- sprintf("%s, ioo %d", method, ioo)
      rl <- xp_mnl(d, rep(0, 8L), S, ioo, xp_ctl(1e-12, 1000L, method = method))
      rk <- xp_mnl(d, rep(0, 8L), S, ioo, xp_ctl(1e-12, 1000L, method = method,
                                                  compensate = "on", pass = "exp"))
      expect_identical(rk$code, 0L, label = lab)
      expect_gt(rk$exp_evals, 0)
      expect_equal(as.numeric(rk$delta), as.numeric(rl$delta), tolerance = 1e-11,
                   label = lab)
      rl <- xp_nl(d, rep(0, 8L), S2, ioo, nest, lam, xp_ctl(1e-12, 2000L, method = method))
      rk <- xp_nl(d, rep(0, 8L), S2, ioo, nest, lam,
                  xp_ctl(1e-12, 2000L, method = method, compensate = "on", pass = "exp"))
      expect_identical(rk$code, 0L, label = lab)
      expect_gt(rk$exp_evals, 0)
      expect_equal(as.numeric(rk$delta), as.numeric(rl$delta), tolerance = 1e-10,
                   label = lab)
    }
  }
  skip_if_not(isTRUE(thread_info()$openmp_enabled))
  set_omp_threads(2L)
  rk <- xp_mnl(d, rep(0, 8L), S, TRUE, xp_ctl(1e-12, 1000L, method = "anderson",
                                               compensate = "on", pass = "exp"))
  expect_equal(as.numeric(rk$delta), ds, tolerance = 1e-10)
})

test_that("the mixed logit's compensated sums keep its draws and its prediction", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_design(2202L, 200L, 6L, 4L)
  set.seed(43L)
  W <- matrix(stats::rnorm(2L * length(d$alt)), ncol = 2L)
  mu <- c(0.3, -0.2)
  Lp <- c(log(0.7), 0.2, log(0.5))
  S_ <- 15L
  gen <- list(gen_seed = 3L, gen_scramble = 1L, gen_S = S_)
  empty <- array(0, c(2L, 0L, 0L))
  for (ioo in c(FALSE, TRUE)) {
    ds <- xp_norm(c(0, 0.3, -0.2, 0.5, 0.1, -0.4), ioo)
    s1 <- as.numeric(do.call(choicer:::mxl_predict_shares, c(list(
      c(xp_beta, mu, Lp, xp_free(ds, ioo)), d$X, W, d$alt, d$M, d$w, empty,
      c(0L, 0L), TRUE, TRUE, TRUE, ioo), gen)))
    solve <- function(start, ctl, budget) {
      do.call(choicer:::mxl_blp_solve, c(list(
        start, s1, d$X, W, xp_beta, mu, Lp, d$alt, d$M, d$w, empty, c(0L, 0L), TRUE,
        TRUE, ioo, ctl), gen, list(cache_bytes = budget)))
    }
    lab <- sprintf("ioo %d", ioo)
    r <- solve(ds, xp_ctl(max_iter = 1L, keep_pass = TRUE, compensate = "on"), 0)
    expect_true(r$compensated, label = lab)
    expect_equal(xp_pass(r)$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
    if (ioo) expect_equal(xp_pass(r)$s0, s1[1L], tolerance = 1e-14, label = lab)
    rk <- solve(rep(0, 6L), xp_ctl(1e-12, 1000L, method = "anderson", compensate = "on"), Inf)
    rl <- solve(rep(0, 6L), xp_ctl(1e-12, 1000L, method = "anderson"), Inf)
    expect_identical(rk$code, 0L, label = lab)
    expect_equal(as.numeric(rk$delta), as.numeric(rl$delta), tolerance = 1e-11, label = lab)
    # the draws kept as before
    expect_true(rk$kept, label = lab)
    expect_identical(rk$reads, length(d$M) * (rk$passes - 1), label = lab)
  }
})

# -- The mixed logit's exp-space cache ------------------------------------------

xp_mxl_design <- function(seed, N, J, m) {
  d <- xp_design(seed, N, J, m)
  set.seed(seed + 1L)
  d$W <- matrix(stats::rnorm(2L * length(d$alt)), ncol = 2L)
  d
}
xp_mu <- c(0.3, -0.2)
xp_Lp <- c(log(0.7), 0.2, log(0.5))
# draws formed on the fly (seeded) or a cube in memory
xp_draws <- function(N, S = 12L) {
  list(generate = list(eta = array(0, c(2L, 0L, 0L)), gen_seed = 5L, gen_scramble = 1L,
                       gen_S = S),
       cube = list(eta = get_halton_normals(S, N, 2L), gen_seed = -1L, gen_scramble = 1L,
                   gen_S = 0L))
}
xp_mxl_s <- function(d, delta, ioo, g) {
  as.numeric(choicer:::mxl_predict_shares(
    c(xp_beta, xp_mu, xp_Lp, xp_free(delta, ioo)), d$X, d$W, d$alt, d$M, d$w, g$eta,
    c(0L, 0L), TRUE, TRUE, TRUE, ioo, g$gen_seed, g$gen_scramble, g$gen_S))
}
xp_mxl <- function(d, start, S, ioo, g, ctl, budget) {
  choicer:::mxl_blp_solve(start, S, d$X, d$W, xp_beta, xp_mu, xp_Lp, d$alt, d$M, d$w,
                          g$eta, c(0L, 0L), TRUE, TRUE, ioo, ctl, g$gen_seed,
                          g$gen_scramble, g$gen_S, budget)
}

test_that("a mixed logit's evaluation from the cache gives the prediction kernel's shares", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_mxl_design(2301L, 160L, 6L, 4L)
  N <- length(d$M)
  n <- length(d$alt)
  for (gn in c("generate", "cube")) {
    g <- xp_draws(N)[[gn]]
    S_ <- 12
    gamma_b <- if (gn == "generate") 8 * 2 * S_ * N else 0
    for (ioo in c(FALSE, TRUE)) {
      exp_b <- 8 * S_ * (n + ioo * N)
      set.seed(61L + ioo)
      ds <- xp_norm(stats::rnorm(6L), ioo)
      d1 <- xp_norm(ds + stats::rnorm(6L, 0, 0.3), ioo)
      S <- xp_mxl_s(d, xp_norm(stats::rnorm(6L), ioo), ioo, g)
      s1 <- xp_mxl_s(d, d1, ioo, g)
      rl <- xp_mxl(d, d1, S, ioo, g, xp_ctl(max_iter = 1L, keep_pass = TRUE), 0)
      for (bu in c("all", "half")) {
        budget <- if (bu == "all") Inf else gamma_b + exp_b / 2
        lab <- sprintf("%s, ioo %d, %s", gn, ioo, bu)
        rx <- xp_mxl(d, d1, S, ioo, g,
                     xp_ctl(max_iter = 1L, keep_pass = TRUE, pass = "exp", origin = ds),
                     budget)
        expect_true(rx$exp_kept, label = lab)
        expect_identical(c(rx$exp_evals, rx$exp_forms, rx$exp_floors), c(1, 1, 0),
                         label = lab)
        if (bu == "all") {
          expect_identical(rx$exp_prefix, as.numeric(N), label = lab)
        } else {
          expect_gt(rx$exp_prefix, 0)
          expect_lt(rx$exp_prefix, N)
        }
        expect_identical(rx$exp_reads, rx$exp_prefix, label = lab)
        px <- xp_pass(rx)
        expect_equal(px$s, if (ioo) s1[-1L] else s1, tolerance = 1e-14, label = lab)
        if (ioo) expect_equal(px$s0, s1[1L], tolerance = 1e-14, label = lab)
        pl <- xp_pass(rl)
        expect_equal(px$A, pl$A, tolerance = 1e-13, label = lab)
        expect_equal(px$C, pl$C, tolerance = 1e-13, label = lab)
      }
    }
  }
})

test_that("the mixed logit's exp-space pass solves as its usual pass does", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_mxl_design(2302L, 200L, 6L, 4L)
  N <- length(d$M)
  g <- xp_draws(N)$generate
  for (ioo in c(FALSE, TRUE)) {
    ds <- xp_norm(c(0, 0.6, -0.4, 0.9, 0.2, -0.7), ioo)
    S <- xp_mxl_s(d, ds, ioo, g)
    for (method in c("simple", "anderson")) {
      lab <- sprintf("%s, ioo %d", method, ioo)
      rl <- xp_mxl(d, rep(0, 6L), S, ioo, g, xp_ctl(1e-12, 1000L, method = method), Inf)
      rx <- xp_mxl(d, rep(0, 6L), S, ioo, g,
                   xp_ctl(1e-12, 1000L, method = method, pass = "exp"), Inf)
      expect_identical(rx$code, 0L, label = lab)
      expect_equal(as.numeric(rx$delta), ds, tolerance = 1e-10, label = lab)
      expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-11,
                   label = lab)
      expect_gt(rx$exp_evals, 0)
      # every situation of every pass after the first is served by a cache
      expect_true(rx$kept, label = lab)
      expect_identical(rx$reads + rx$exp_reads, N * (rx$passes - 1), label = lab)
      expect_identical(rx$exp_reads, rx$exp_prefix * (rx$exp_evals + rx$exp_floors),
                       label = lab)
    }
  }
  skip_if_not(isTRUE(thread_info()$openmp_enabled))
  set_omp_threads(2L)
  S <- xp_mxl_s(d, ds, TRUE, g)
  r2 <- xp_mxl(d, rep(0, 6L), S, TRUE, g, xp_ctl(1e-12, 1000L, method = "anderson",
                                                    pass = "exp"), Inf)
  expect_equal(as.numeric(r2$delta), ds, tolerance = 1e-10, label = "two threads")
  expect_gt(r2$exp_evals, 0)
  expect_identical(r2$reads + r2$exp_reads, N * (r2$passes - 1), label = "two threads")
})

test_that("a mixed logit's evaluation beyond the span, or below its floor, is the usual pass", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_mxl_design(2304L, 120L, 5L, 3L)
  N <- length(d$M)
  g <- xp_draws(N)$generate
  ds <- c(0, 0.3, -0.2, 0.5, 0.1)
  S <- xp_mxl_s(d, ds, FALSE, g)
  one <- function(d1, origin) {
    list(x = xp_mxl(d, d1, S, FALSE, g, xp_ctl(max_iter = 1L, keep_pass = TRUE,
                                                 pass = "exp", origin = origin), Inf),
         l = xp_mxl(d, d1, S, FALSE, g, xp_ctl(max_iter = 1L, keep_pass = TRUE), Inf))
  }
  # beyond the span: the cache formed at the origin, then the usual pass at d1
  r <- one(replace(ds, 3L, ds[3L] + 51), ds)
  expect_identical(c(r$x$exp_evals, r$x$exp_forms), c(0, 2))
  expect_identical(xp_pass(r$x), xp_pass(r$l))
  expect_identical(r$x$passes, 2)
  expect_identical(r$x$reads + r$x$exp_reads, N * (r$x$passes - 1))
  # below the floor at d1 (above it at the origin): the evaluation from the
  # cache, then the usual pass (the random coefficients' tails lift the share
  # some e^12 above the logit's, hence the deeper constants)
  d0 <- replace(ds, 5L, -575)
  d1 <- replace(ds, 5L, -624)
  expect_gt(xp_mxl_s(d, d0, FALSE, g)[5L], xp_floor(50))
  expect_lt(xp_mxl_s(d, d1, FALSE, g)[5L], xp_floor(49))
  r <- one(d1, d0)
  expect_identical(c(r$x$exp_evals, r$x$exp_forms, r$x$exp_floors), c(0, 2, 1))
  expect_identical(xp_pass(r$x), xp_pass(r$l))
  expect_identical(r$x$passes, 3)
  expect_identical(r$x$reads + r$x$exp_reads, N * (r$x$passes - 1))
})

test_that("the mixed logit's exp-space cache takes what Gamma leaves of the budget", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_mxl_design(2303L, 120L, 5L, 3L)
  N <- length(d$M)
  n <- length(d$alt)
  g <- xp_draws(N)$generate
  gamma_b <- 8 * 2 * 12 * N
  S <- xp_mxl_s(d, c(0, 0.3, -0.2, 0.5, 0.1), FALSE, g)
  rl <- xp_mxl(d, rep(0, 5L), S, FALSE, g, xp_ctl(1e-10, 1000L), Inf)
  for (b in list(list(0, FALSE, FALSE), list(NaN, FALSE, FALSE),
                 list(gamma_b, TRUE, FALSE), list(gamma_b + 8 * 12 * 2, TRUE, FALSE),
                 list(gamma_b + 8 * 12 * 3, TRUE, TRUE),
                 list(gamma_b + 8 * 12 * n, TRUE, TRUE), list(Inf, TRUE, TRUE))) {
    lab <- format(b[[1]])
    rx <- xp_mxl(d, rep(0, 5L), S, FALSE, g, xp_ctl(1e-10, 1000L, pass = "exp"), b[[1]])
    expect_identical(rx$kept, b[[2]], label = lab)
    expect_identical(rx$exp_kept, b[[3]], label = lab)
    if (!rx$exp_kept) {
      expect_identical(rx$exp_prefix, 0, label = lab)
      expect_identical(rx$exp_reads, 0, label = lab)
      expect_identical(rx$delta, rl$delta, label = lab)
    }
  }
  # the longest prefix the budget holds, in both layouts (three rows per
  # situation, twelve draws, and one outside option's value per draw)
  for (ioo in c(FALSE, TRUE)) {
    So <- xp_mxl_s(d, xp_norm(c(0, 0.3, -0.2, 0.5, 0.1), ioo), ioo, g)
    per <- 8 * 12 * (3 + ioo)
    for (b in list(c(per, 1), c(N * per - 1, N - 1), c(N * per, N))) {
      rx <- xp_mxl(d, rep(0, 5L), So, ioo, g, xp_ctl(1e-10, 1000L, pass = "exp"),
                   gamma_b + b[1L])
      expect_identical(rx$exp_prefix, b[2L], label = sprintf("ioo %d, %g", ioo, b[1L]))
    }
    # control's cache_bytes caps it too
    rx <- xp_mxl(d, rep(0, 5L), So, ioo, g,
                 xp_ctl(1e-10, 1000L, pass = "exp", cache_bytes = 0), Inf)
    expect_false(rx$exp_kept, label = sprintf("ioo %d", ioo))
  }
  # a budget the cap leaves no room under keeps Gamma (within the budget, as
  # before) and no exponentials
  old <- options(choicer.max_memory = 1e4)
  on.exit(options(old), add = TRUE)
  rx <- xp_mxl(d, rep(0, 5L), S, FALSE, g, xp_ctl(1e-10, 1000L, pass = "exp"), Inf)
  options(old)
  expect_true(rx$kept)
  expect_false(rx$exp_kept)
})

test_that("a cube's partial prefix solves with compensated sums", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- xp_mxl_design(2305L, 200L, 6L, 4L)
  N <- length(d$M)
  g <- xp_draws(N)$cube
  ds <- c(0, 0.6, -0.4, 0.9, 0.2, -0.7)
  S <- xp_mxl_s(d, ds, FALSE, g)
  rl <- xp_mxl(d, rep(0, 6L), S, FALSE, g, xp_ctl(1e-12, 1000L, method = "anderson"), Inf)
  rx <- xp_mxl(d, rep(0, 6L), S, FALSE, g,
               xp_ctl(1e-12, 1000L, method = "anderson", pass = "exp", compensate = "on"),
               8 * 12 * length(d$alt) / 2)
  expect_identical(rx$code, 0L)
  expect_true(rx$compensated)
  expect_false(rx$kept)  # (a cube is read in place)
  expect_gt(rx$exp_prefix, 0)
  expect_lt(rx$exp_prefix, N)
  expect_gt(rx$exp_evals, 0)
  expect_equal(as.numeric(rx$delta), as.numeric(rl$delta), tolerance = 1e-11)
  expect_identical(rx$exp_reads, rx$exp_prefix * (rx$exp_evals + rx$exp_floors))
})
