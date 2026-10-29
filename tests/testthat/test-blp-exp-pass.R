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
