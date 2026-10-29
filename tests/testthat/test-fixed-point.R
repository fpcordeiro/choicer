# The fixed-point solver (src/fixed_point.h) on a test map
# (test_fp_solve(), src/kernel_test_exports.cpp): the linear step
# f(x) = b - A x, residuals scripted for the first evaluations, evaluations
# that fail on a schedule, and a fallback step f(x) / 2; every evaluated point
# is recorded (trace_x, one column per evaluation).

fp_ctl <- function(...) choicer:::.blp_control(...)
fp_run <- function(A, b, x0, ctl, fail = integer(), script = numeric(),
                   fallback = FALSE) {
  choicer:::test_fp_solve(as.matrix(A), b, x0, ctl, as.integer(fail),
                          as.numeric(script), fallback)
}
fp_step <- function(A, b, x) b - drop(A %*% x)

test_that("Anderson acceleration solves a linear map in a few iterations", {
  A <- diag(c(1, 0.01, 0.001))
  xs <- c(1, 2, 3)
  b <- drop(A %*% xs)
  r <- fp_run(A, b, rep(0, 3L), fp_ctl(tol = 1e-10, max_iter = 50L, method = "anderson"))
  expect_identical(r$code, 0L)
  expect_lte(r$evaluations, 10L)
  expect_equal(r$x, xs, tolerance = 1e-8)
  # the plain iteration crawls along the slowest coordinate (rate 0.999)
  p <- fp_run(A, b, rep(0, 3L), fp_ctl(tol = 1e-10, max_iter = 500L, patience = 0L,
                                       window = 0L))
  expect_identical(p$code, 1L)
})

test_that("a base step whose trial point fails is halved, without a restart", {
  A <- diag(c(1, 0.5, 0.25))
  b <- drop(A %*% c(1, 2, 3))
  x0 <- rep(0, 3L)
  # the first step from the start (evaluation 2) fails: the next trial is half
  # of it, not the same point again
  r <- fp_run(A, b, x0, fp_ctl(tol = 1e-10, max_iter = 200L, method = "anderson"),
              fail = 2L)
  expect_identical(r$code, 0L)
  expect_identical(r$restarts, 0L)
  expect_identical(r$backtracks, 1L)
  expect_equal(r$trace_x[, 3L], x0 + 0.5 * fp_step(A, b, x0))
})

test_that("an extrapolated point that fails restarts from where the iteration is", {
  A <- diag(c(1, 0.5, 0.25))
  b <- drop(A %*% c(1, 2, 3))
  x0 <- rep(0, 3L)
  # evaluation 3 is the first extrapolated point (one difference kept); it
  # fails, and evaluation 4 is the step from evaluation 2's point
  r <- fp_run(A, b, x0, fp_ctl(tol = 1e-10, max_iter = 200L, method = "anderson"),
              fail = 3L)
  expect_identical(r$code, 0L)
  expect_identical(r$restarts, 1L)
  expect_identical(r$backtracks, 0L)
  x2 <- r$trace_x[, 2L]
  expect_equal(r$trace_x[, 4L], x2 + fp_step(A, b, x2))
  expect_equal(r$x, c(1, 2, 3), tolerance = 1e-8)
})

test_that("a stalled iteration restarts in place, or from its best point where it diverges", {
  A <- diag(c(1, 0.5, 0.25))
  b <- drop(A %*% c(1, 2, 3))
  x0 <- rep(0, 3L)
  ctl <- fp_ctl(tol = 1e-10, max_iter = 200L, method = "anderson")
  # residuals scripted: the best, 0.5, at evaluation 2, then five iterations
  # without one below 0.45; at the last the residual is 0.95 (at most twice
  # the best: in place) or 1.2 (past it: back to the best point)
  for (last in c(0.95, 1.2)) {
    r <- fp_run(A, b, x0, ctl, script = c(1, 0.5, 0.6, 0.7, 0.8, 0.9, last))
    expect_identical(r$code, 0L)
    expect_identical(r$restarts, 1L)
    from <- if (last > 1) r$trace_x[, 2L] else r$trace_x[, 7L]
    expect_equal(r$trace_x[, 8L], from + fp_step(A, b, from), label = format(last))
  }
  # Phase R's rule: from the best point whatever the residual
  rb <- fp_run(A, b, x0, fp_ctl(tol = 1e-10, max_iter = 200L, method = "anderson",
                                restart_best = TRUE),
               script = c(1, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95))
  expect_equal(rb$trace_x[, 8L], rb$trace_x[, 2L] + fp_step(A, b, rb$trace_x[, 2L]))
})

test_that("a restart whose step keeps failing ends the iteration, and the fallback takes over", {
  A <- diag(c(1, 0.5, 0.25))
  b <- drop(A %*% c(1, 2, 3))
  x0 <- rep(0, 3L)
  ctl <- fp_ctl(tol = 1e-10, max_iter = 2000L, method = "anderson")
  # the extrapolated point (3) and the restart's step and its ten halvings
  # (4 to 14) fail: the fallback starts from the start at evaluation 15 and
  # converges alone
  r <- fp_run(A, b, x0, ctl, fail = 3:14, fallback = TRUE)
  expect_identical(r$code, 0L)
  expect_identical(r$restarts, 1L)
  expect_identical(r$backtracks, 10L)
  expect_identical(r$race, 14L)
  expect_true(r$fallback)
  expect_equal(r$trace_x[, 15L], x0)
  expect_identical(r$trace_fb[15L], 1L)
  # without a fallback the iteration stops there
  r0 <- fp_run(A, b, x0, ctl, fail = 3:14)
  expect_identical(r0$code, 2L)
  expect_identical(r0$evaluations, 14L)
  # and the budget can end inside the halvings
  r6 <- fp_run(A, b, x0, fp_ctl(tol = 1e-10, max_iter = 6L, method = "anderson"),
               fail = 3:14, fallback = TRUE)
  expect_identical(r6$code, 1L)
  expect_identical(r6$evaluations, 6L)
  expect_identical(r6$race, 0L)
})

test_that("Anderson's least squares drop differences of steps that vanish", {
  # A constant step (A = 0): every difference of steps is zero, a zero
  # diagonal of the triangular factor, so none is kept and every step is the
  # plain one (the residual never falls, so the iteration restarts every five)
  r <- fp_run(matrix(0, 2L, 2L), c(1, -1), c(0, 0),
              fp_ctl(tol = 1e-12, max_iter = 31L, method = "anderson"))
  expect_identical(r$code, 1L)
  expect_identical(r$evaluations, 31L)
  expect_true(all(is.finite(r$trace_x)))
  expect_equal(r$trace_x[, 31L], c(30, -30))
  expect_identical(r$restarts, 5L)
  # A multiple of the identity: the first extrapolation is exact
  A <- diag(rep(0.3, 4L))
  r <- fp_run(A, drop(A %*% c(1, -1, 2, 0)), rep(0, 4L),
              fp_ctl(tol = 1e-12, max_iter = 100L, method = "anderson"))
  expect_identical(r$code, 0L)
  expect_identical(r$evaluations, 3L)
  expect_equal(r$x, c(1, -1, 2, 0), tolerance = 1e-12)
  # memory beyond the coordinates is cut to them
  r1 <- fp_run(diag(0.2, 1L), 0.2, 0, fp_ctl(tol = 1e-12, max_iter = 100L,
                                             method = "anderson",
                                             settings = list(memory = 50L)))
  expect_identical(r1$code, 0L)
  expect_equal(r1$x, 1, tolerance = 1e-10)
})

test_that("an infinite residual at the start counts no stagnation", {
  A <- diag(c(1, 0.5))
  b <- drop(A %*% c(1, 2))
  # residuals infinite for the first eight evaluations (as where an outside
  # share underflows at the start): no restart while there is no best point
  r <- fp_run(A, b, c(0, 0), fp_ctl(tol = 1e-10, max_iter = 200L, method = "anderson"),
              script = rep(Inf, 8L))
  expect_identical(r$code, 0L)
  expect_identical(r$restarts, 0L)
})
