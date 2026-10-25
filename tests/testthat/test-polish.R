# The Newton phase (the run_*logit() wrappers' polish): newton_step(), the
# damped Cholesky step on the free parameters, run_polish(), the projected,
# damped Newton iteration, and the wiring the wrappers share (.polish_held(),
# .polish_fit()), here on synthetic objectives with known optima.

# An objective in run_polish()'s form: f, its gradient, an overflow flag
pol_obj <- function(f, g, overflow = FALSE) {
  list(objective = f, gradient = g, overflow = overflow)
}

# The quadratic 1/2 x'Ax - b'x + c, minimized at solve(A, b)
pol_quad <- function(A, b, c = 0) {
  list(f = function(x) pol_obj(0.5 * sum(x * (A %*% x)) - sum(b * x) + c,
                               as.numeric(A %*% x - b)),
       C = function(x) A)
}

# Rosenbrock's function, minimized at (1, 1); its Hessian is indefinite where
# x2 > 3 x1^2 + 0.005
pol_rosen <- list(
  f = function(x) {
    pol_obj(100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2,
            c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
              200 * (x[2] - x[1]^2)))
  },
  C = function(x) {
    matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
             -400 * x[1], 200), 2)
  }
)

# A random symmetric positive definite matrix of order n
pol_spd <- function(n, seed) {
  set.seed(seed)
  X <- matrix(stats::rnorm(4 * n * n), 4 * n)
  crossprod(X) / n + diag(n)
}

test_that("newton_step() solves the scaled, damped system on the free set", {
  for (n in c(6L, 65L, 130L)) {
    H <- pol_spd(n, n)
    # spread the curvatures over orders of magnitude
    s <- 10^seq(-3, 3, length.out = n)
    H <- H * tcrossprod(s)
    set.seed(n + 1L)
    g <- stats::rnorm(n)
    D <- .jacobi_scale(H)
    r <- newton_step(H, g, D, seq_len(n), 0)
    expect_identical(r$info, 0L)
    expect_equal(r$step, -solve(H, g), tolerance = 1e-10)

    # a free subset: the others do not move
    fr <- sort(sample.int(n, n %/% 2L))
    r <- newton_step(H, g, D, fr, 0)
    expect_identical(r$info, 0L)
    expect_equal(r$step[fr], -solve(H[fr, fr], g[fr]), tolerance = 1e-10)
    expect_identical(r$step[-fr], rep(0, n - length(fr)))

    # lambda on the scaled block's diagonal: Marquardt's damping
    lam <- 0.37
    Ds <- diag(D[fr])
    ref <- Ds %*% solve(Ds %*% H[fr, fr] %*% Ds + lam * diag(length(fr)),
                        -D[fr] * g[fr])
    r <- newton_step(H, g, D, fr, lam)
    expect_equal(r$step[fr], as.numeric(ref), tolerance = 1e-10)
  }
})

test_that("newton_step() reports a block that is not positive definite", {
  H <- pol_spd(5L, 1L)
  g <- seq_len(5) / 5
  H[3, 3] <- -1
  r <- newton_step(H, g, .jacobi_scale(H), 1:5, 0)
  expect_gt(r$info, 0L)
  expect_identical(r$step, rep(0, 5))
  # enough damping makes it positive definite
  expect_identical(newton_step(H, g, .jacobi_scale(H), 1:5, 1e3)$info, 0L)
  # the indefinite entry outside the free set does not matter
  expect_identical(newton_step(H, g, .jacobi_scale(H), c(1L, 2L, 4L), 0)$info,
                   0L)
})

test_that("newton_step() flags values that are not finite, scaled or not", {
  H <- pol_spd(4L, 2L)
  g <- c(1, -1, 0.5, 2)
  D <- .jacobi_scale(H)
  Hn <- H
  Hn[2, 4] <- Hn[4, 2] <- NaN
  expect_identical(newton_step(Hn, g, D, 1:4, 0)$info, -1L)
  expect_identical(newton_step(H, replace(g, 3, Inf), D, 1:4, 0)$info, -1L)
  # outside the free set a non-finite value is not read
  expect_identical(newton_step(Hn, g, D, 1:3, 0)$info, 0L)
  # the scaling itself overflows: a huge off-diagonal entry times 2^1000
  Ho <- H
  Ho[1, 2] <- Ho[2, 1] <- 1e300
  expect_identical(newton_step(Ho, g, rep(2^500, 4), 1:4, 0)$info, -1L)
})

test_that("newton_step() treats a block singular to working precision as unusable", {
  # positive definite, its Cholesky factorization succeeds, but its
  # reciprocal condition number, about 2^-53, is below the machine epsilon:
  # the code is order + 1
  H <- matrix(c(1, 1 - 2^-52, 1 - 2^-52, 1), 2)
  expect_identical(newton_step(H, c(1, 1), c(1, 1), 1:2, 0)$info, 3L)
  # damping makes it usable
  expect_identical(newton_step(H, c(1, 1), c(1, 1), 1:2, 1e-3)$info, 0L)
})

test_that("newton_step() scales by the product of the scales, exactly", {
  # D = (2^500, 2^-500): the off-diagonal entry 2^600 times D_1 alone would
  # overflow; times D_1 D_2 = 1 it does not, so the indefinite block is
  # reported as such rather than as not finite
  H <- matrix(c(2^-1000, 2^600, 2^600, 2^1000), 2)
  D <- .jacobi_scale(H)
  expect_identical(D, c(2^500, 2^-500))
  expect_gt(newton_step(H, c(1, 1), D, 1:2, 0)$info, 0L)
})

test_that("newton_step() checks its arguments", {
  H <- pol_spd(3L, 3L)
  g <- c(1, 2, 3)
  D <- c(1, 1, 1)
  r <- newton_step(H, g, D, integer(), 0)
  expect_identical(r, list(step = c(0, 0, 0), info = 0L))
  expect_error(newton_step(H, g, D, c(2L, 1L), 0), "strictly increasing")
  expect_error(newton_step(H, g, D, c(1L, 1L), 0), "strictly increasing")
  expect_error(newton_step(H, g, D, c(0L, 1L), 0), "strictly increasing")
  expect_error(newton_step(H, g, D, c(1L, 4L), 0), "strictly increasing")
  expect_error(newton_step(H, g, D, c(1L, NA), 0), "strictly increasing")
  expect_error(newton_step(H, g, D, c(1, 2), 0), "integer vector")
  expect_error(newton_step(H, g, D, 1:3, -1), "nonnegative")
  expect_error(newton_step(H, g, D, 1:3, NaN), "nonnegative")
  expect_error(newton_step(H[, 1:2], g, D, 1:2, 0), "square")
  expect_error(newton_step(1:9, g, D, 1:3, 0), "double matrix")
  expect_error(newton_step(H, g[1:2], D, 1:3, 0), "length 3")
  expect_error(newton_step(H, g, D[1:2], 1:3, 0), "length 3")
})

test_that(".polish_control() merges settings with the defaults and checks them", {
  def <- .polish_control(NULL)
  expect_identical(def, list(max_iter = 20L, tol_decrement = 1e-8,
                             tol_step = 1e-12))
  expect_identical(.polish_control(list()), def)
  expect_identical(.polish_control(list(max_iter = 5)),
                   modifyList(def, list(max_iter = 5L)))
  expect_identical(.polish_control(list(tol_decrement = NULL)), def)
  expect_error(.polish_control(list(maxiter = 5)), "Unknown `polish_control`")
  expect_error(.polish_control(list(max_iter = 1, max_iter = 2)), "Duplicate")
  expect_error(.polish_control(list(5)), "named list")
  expect_error(.polish_control(c(max_iter = 5)), "must be a list")
  expect_error(.polish_control(list(max_iter = 1.5)), "whole number")
  expect_error(.polish_control(list(max_iter = -1)), "whole number")
  expect_error(.polish_control(list(max_iter = NA)), "whole number")
  expect_error(.polish_control(list(max_iter = TRUE)), "whole number")
  expect_error(.polish_control(list(tol_decrement = 0)), "positive number")
  expect_error(.polish_control(list(tol_step = c(1, 2))), "positive number")
  expect_error(.polish_control(list(tol_step = Inf)), "positive number")
})

test_that("a quadratic is solved in one step and certified", {
  q <- pol_quad(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  r <- run_polish(c(a = 5, b = -3), q$f, q$C)
  expect_true(r$converged)
  expect_identical(r$iterations, 1L)
  expect_identical(r$damped, 0L)
  expect_identical(r$lambda, 0)
  expect_lt(r$decrement, 1e-20)
  expect_named(r$par, c("a", "b"))
  expect_equal(unname(r$par), solve(matrix(c(4, 1, 1, 3), 2), c(1, 2)),
               tolerance = 1e-14)
  expect_identical(r$message, "the Newton decrement is below tol_decrement")
  expect_identical(r$evaluations, c(objective = 2L, curvature = 2L))
  expect_gt(r$improvement, 0)
  expect_null(r$error)
  expect_length(r$active, 0L)
  expect_length(r$held, 0L)
  expect_identical(r$method, "newton")
  # the curvature's name, in the result and the message
  r <- run_polish(c(a = 5, b = -3), q$f, q$C, method = "bhhh")
  expect_identical(r$method, "bhhh")
  expect_identical(r$message, "the BHHH decrement is below tol_decrement")
})

test_that("an indefinite start takes damped steps to the minimum", {
  # (0, 1): the Hessian's first diagonal entry is -398
  r <- run_polish(c(0, 1), pol_rosen$f, pol_rosen$C,
                  control = .polish_control(list(max_iter = 100L)))
  expect_true(r$converged)
  expect_gt(r$damped, 0L)
  expect_identical(r$lambda, 0)
  # a decrement below 1e-8 bounds the objective's gap, not the coordinates'
  expect_lt(r$decrement, 1e-8)
  expect_lt(r$value, 1e-8)
  expect_lt(max(abs(r$par - 1)), 1e-4)
  # without damping the indefinite start cannot be left
  r0 <- run_polish(c(0, 1), pol_rosen$f, pol_rosen$C,
                   control = .polish_control(list(max_iter = 0L)))
  expect_false(r0$converged)
  expect_gt(r0$lambda, 0)
  expect_identical(r0$par, c(0, 1))
  # the decrement reported is the undamped one, which an indefinite block
  # does not give: NA, not the damped step's
  expect_identical(r0$decrement, NA_real_)
})

test_that("bounds hold and an active bound is certified at its KKT point", {
  A <- matrix(c(4, 1, 1, 3), 2)
  q <- pol_quad(A, c(1, 2))
  # the minimum has a = 1/11; an upper bound below it binds
  r <- run_polish(c(a = 0, b = 0), q$f, q$C, upper = c(0.05, Inf))
  expect_true(r$converged)
  expect_identical(r$active, "a")
  expect_identical(unname(r$par[1]), 0.05)
  # b minimizes the quadratic with a held at the bound
  expect_equal(unname(r$par[2]), (2 - 0.05) / 3, tolerance = 1e-14)
  # the gradient of a points across the bound
  expect_lt(q$f(r$par)$gradient[1], 0)

  # strong coupling: the first coordinate reaches its bound during the polish
  Ac <- matrix(c(1, 0.95, 0.95, 1), 2)
  qc <- pol_quad(Ac, c(1, -2))
  x_free <- solve(Ac, c(1, -2))
  expect_gt(x_free[1], 2)
  r <- run_polish(c(0, 0), qc$f, qc$C, lower = c(-Inf, -Inf),
                  upper = c(2, Inf))
  expect_true(r$converged)
  expect_identical(r$par[1], 2)
  expect_equal(r$par[2], -2 - 0.95 * 2, tolerance = 1e-14)
  expect_identical(r$active, "1")
  expect_true(all(r$par <= c(2, Inf)))

  # a start outside the box is projected into it first
  r <- run_polish(c(a = 1, b = 0), q$f, q$C, upper = c(0.05, Inf))
  expect_identical(unname(r$par[1]), 0.05)
  expect_true(r$converged)
})

test_that("a parameter pulled across its bound does not slow the others", {
  # f = 1/2 |x - (100, 1)|^2 with x1 <= 0: x1 sits on its bound, its gradient
  # pointing across it; x2 takes its full Newton step at once
  f <- function(x) pol_obj(0.5 * sum((x - c(100, 1))^2), x - c(100, 1))
  r <- run_polish(c(0, 0), f, function(x) diag(2), upper = c(0, Inf))
  expect_true(r$converged)
  expect_identical(r$iterations, 1L)
  expect_identical(r$par, c(0, 1))
  expect_identical(r$active, "1")
})

test_that("an almost-active parameter is held on its bound (epsilon-active)", {
  # x1 starts 1e-6 above its lower bound with its gradient pointing across
  # it: the epsilon-active set moves it onto the bound in the first step,
  # where an exact active set would take the coupled step first
  A <- matrix(c(1, 0.95, 0.95, 1), 2)
  q <- pol_quad(A, c(-1, 1))
  r <- run_polish(c(1e-6, 0), q$f, q$C, lower = c(0, -Inf))
  expect_true(r$converged)
  expect_identical(r$iterations, 1L)
  expect_identical(r$evaluations, c(objective = 2L, curvature = 2L))
  expect_identical(r$par[1], 0)
})

test_that("a projected step that ascends at first is halved, not damped", {
  # Near a lower bound the projected path of the Newton step ascends at
  # alpha = 1 and descends for smaller alpha
  A <- matrix(c(1, -0.9, -0.9, 1), 2)
  x0 <- c(0.1, 0)
  b <- as.numeric(A %*% x0) - c(0.3, -0.2)
  q <- pol_quad(A, b)
  r <- run_polish(x0, q$f, q$C, lower = c(0, -Inf),
                  control = .polish_control(list(max_iter = 1L)))
  expect_identical(r$iterations, 1L)
  expect_identical(r$damped, 0L)
  expect_lt(r$value, q$f(x0)$objective)
})

test_that("a parameter whose step runs off is held and named", {
  # The curvature in x2 is a trillion times too small, as BHHH's is in the
  # log-Cholesky diagonal of a variance going to zero: its step would carry
  # it 1e11 away. It is held; x1 reaches its minimum, but the polish is not
  # certified
  f <- function(x) pol_obj(0.5 * sum((x - c(1, 0))^2), x - c(1, 0))
  r <- run_polish(c(a = 0, b = 0.5), f, function(x) diag(c(1, 1e-12)))
  expect_false(r$converged)
  expect_identical(r$flat, "b")
  expect_identical(unname(r$par), c(1, 0.5))
  expect_match(r$message, "b has no usable curvature")
  expect_match(r$message, paste("the Newton decrement over the other",
                               "parameters is below tol_decrement, but b"),
               fixed = TRUE)
  expect_match(r$message, "a curvature all but vanished", fixed = TRUE)
})

test_that("a trial point whose objective overflows is backtracked", {
  # f = (x - 1)^2 overflows beyond x = 3; a curvature ten times too small
  # makes the first Newton step 10 times too long
  calls <- numeric()
  f <- function(x) {
    calls <<- c(calls, x)
    if (x > 3) return(pol_obj(1e10, 0, overflow = TRUE))
    pol_obj((x - 1)^2, 2 * (x - 1))
  }
  r <- run_polish(-1, f, function(x) matrix(0.2),
                  control = .polish_control(list(max_iter = 200L)))
  expect_true(any(calls > 3))
  expect_lt(abs(r$par - 1), 1e-4)
  expect_lt(r$value, 1e-8)
  # no point it accepted overflowed
  expect_lte(r$par, 3)
})

test_that("a start whose objective overflows is not certified", {
  f <- function(x) pol_obj(1e10, rep(0, 2), overflow = TRUE)
  r <- run_polish(c(a = 1, b = 2), f, function(x) diag(2))
  expect_false(r$converged)
  expect_identical(r$par, c(a = 1, b = 2))
  expect_identical(r$value, NA_real_)
  expect_match(r$message, "objective is not finite at the start")
  expect_identical(r$evaluations, c(objective = 1L, curvature = 0L))
})

test_that("at the objective's noise floor the gradient's rule accepts a step", {
  # f = C + 1/2 x'Ax - b'x, plus 1e-4 at any point but the start: near the
  # optimum the Newton step's gain (1e-5) is below that noise, so Armijo
  # refuses it; the approximate Wolfe condition, which reads the exact
  # gradient, accepts it, and the next decrement certifies the point
  C0 <- 1e9
  A <- diag(c(2, 1))
  b <- c(1, 1)
  x_opt <- solve(A, b)
  x0 <- x_opt + c(sqrt(1e-5), 0)
  f <- function(x) {
    q <- 0.5 * sum(x * (A %*% x)) - sum(b * x) + C0
    pol_obj(q + if (identical(x, x0)) 0 else 1e-4, as.numeric(A %*% x - b))
  }
  r <- run_polish(x0, f, function(x) A)
  expect_true(r$converged)
  expect_equal(r$par, x_opt, tolerance = 1e-12)
  expect_identical(r$iterations, 1L)
})

test_that("the noise rule is not used where the predicted gain is large", {
  # h(x) = -0.2625 x^3 + 1.0375 x^2 - x rises from h(0) = 0 to h(2) = 0.05,
  # with h'(2) = 0. A curvature of 0.5 at 0 asks for a step to 2: uphill, but
  # by less than the allowance 1e-10 * 1e9 and with a flat gradient there.
  # Its predicted gain, 2, is far above the allowance, so only Armijo
  # applies, and the line search halves to x = 1, where h(1) = -0.225.
  C0 <- 1e9
  f <- function(x) {
    pol_obj(C0 - 0.2625 * x^3 + 1.0375 * x^2 - x,
            -0.7875 * x^2 + 2.075 * x - 1)
  }
  r <- run_polish(0, f, function(x) matrix(0.5),
                  control = .polish_control(list(max_iter = 1L)))
  expect_identical(r$iterations, 1L)
  expect_equal(r$par, 1, tolerance = 1e-14)
  expect_lt(r$value, C0)
})

test_that("a failed line search retries with a damped step", {
  # The objective overflows within 10 degrees of the Newton direction, so
  # the undamped step's line search fails; damping turns the step toward
  # the negative gradient, out of that cone
  A <- matrix(c(1, 0.9, 0.9, 1), 2)
  g0 <- c(-1, 0)
  dn <- -solve(A, g0)
  cone <- function(x) {
    d <- x
    nd <- sqrt(sum(d^2))
    nd > 0 && sum(d * dn) / (nd * sqrt(sum(dn^2))) > cos(10 * pi / 180)
  }
  f <- function(x) {
    if (cone(x)) return(pol_obj(1e10, c(0, 0), overflow = TRUE))
    pol_obj(0.5 * sum(x * (A %*% x)) + sum(g0 * x), as.numeric(A %*% x + g0))
  }
  r <- run_polish(c(0, 0), f, function(x) A,
                  control = .polish_control(list(max_iter = 1L)))
  expect_identical(r$iterations, 1L)
  expect_identical(r$damped, 1L)
  expect_lt(r$value, 0)
})

test_that("a line search that always fails, and damping that never helps, stop", {
  # every point but the start overflows
  x0 <- c(1, 1)
  f <- function(x) {
    if (!identical(x, x0)) return(pol_obj(1e10, c(0, 0), overflow = TRUE))
    pol_obj(1, c(1, 1))
  }
  r <- run_polish(x0, f, function(x) diag(2))
  expect_false(r$converged)
  expect_identical(r$message, "the line search failed")
  expect_identical(r$par, x0)
  # a curvature no damping up to 1e10 makes positive definite
  q <- pol_quad(diag(2), c(1, 1))
  r <- run_polish(c(0, 0), q$f, function(x) matrix(c(1, 1e12, 1e12, 1), 2))
  expect_false(r$converged)
  expect_identical(r$message, "no damping makes the curvature positive definite")
  expect_identical(r$decrement, NA_real_)
})

test_that("a gradient that is not finite stops the start and fails a trial", {
  f <- function(x) pol_obj(1, c(NaN, 0))
  r <- run_polish(c(0, 0), f, function(x) diag(2))
  expect_identical(r$message, "the gradient is not finite at the start")
  expect_false(r$converged)
  # a trial point with a gradient that is not finite is backtracked
  q <- pol_quad(diag(2), c(1, 1))
  seen <- list()
  f2 <- function(x) {
    seen[[length(seen) + 1L]] <<- x
    out <- q$f(x)
    if (length(seen) == 2L) out$gradient[1] <- Inf
    out
  }
  r <- run_polish(c(0, 0), f2, q$C)
  expect_true(r$converged)
  expect_equal(r$par, c(1, 1), tolerance = 1e-14)
  expect_identical(seen[[3L]], c(0.5, 0.5))
})

test_that("damping is warm-started from the previous step's", {
  # Rosenbrock from (0, 1): the first steps need damping; later ones, where
  # the diagonal is not positive, start from a hundredth of the last damping
  real <- newton_step
  lams <- numeric()
  local_mocked_bindings(newton_step = function(H, g, D, free_idx, lambda) {
    lams <<- c(lams, lambda)
    real(H, g, D, free_idx, lambda)
  })
  r <- run_polish(c(0, 1), pol_rosen$f, pol_rosen$C,
                  control = .polish_control(list(max_iter = 100L)))
  expect_true(r$converged)
  expect_gt(r$damped, 1L)
  # the first iteration (its diagonal not positive) starts at 1e-6 and grows
  # tenfold until a step is accepted, at 10; the next starts at 10 / 100
  expect_equal(lams[1:8], 10^(-6:1), tolerance = 1e-12)
  expect_equal(lams[9], 0.1, tolerance = 1e-12)
})

test_that("a last point worse than the start returns the start", {
  # f = C + 1/2 (x - 1)^2 with a curvature twice too large, plus 0.06 at
  # every evaluation after the start: each step gains about 1e-5, inside the
  # allowance 1e-10 * 1e9 = 0.1, and the gated rule accepts it; two steps
  # leave the point 0.12 above the start, beyond the allowance
  C0 <- 1e9
  k <- 0L
  f <- function(x) {
    k <<- k + 1L
    pol_obj(C0 + 0.5 * (x - 1)^2 + 0.06 * (k - 1L), x - 1)
  }
  r <- run_polish(1 - 0.005, f, function(x) matrix(2),
                  control = .polish_control(list(max_iter = 3L)))
  expect_false(r$converged)
  expect_identical(r$par, 1 - 0.005)
  expect_match(r$message, "worse than the start, which is returned")
})

test_that("max_iter = 0 certifies the start without moving it", {
  q <- pol_quad(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  x_opt <- solve(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  r <- run_polish(x_opt, q$f, q$C, control = .polish_control(list(max_iter = 0L)))
  expect_true(r$converged)
  expect_identical(r$par, x_opt)
  expect_identical(r$iterations, 0L)
  r <- run_polish(x_opt + 1, q$f, q$C,
                  control = .polish_control(list(max_iter = 0L)))
  expect_false(r$converged)
  expect_identical(r$par, x_opt + 1)
  expect_gt(r$decrement, 1)
  expect_match(r$message, "iteration limit")
})

test_that("a step below tol_step stops the polish", {
  q <- pol_quad(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  r <- run_polish(c(0, 0), q$f, q$C,
                  control = .polish_control(list(tol_step = 1e3)))
  expect_false(r$converged)
  expect_match(r$message, "step is below tol_step")
  expect_identical(r$par, c(0, 0))
})

test_that("a curvature that is not finite, or an error, ends the polish", {
  q <- pol_quad(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  r <- run_polish(c(0, 0), q$f, function(x) matrix(NaN, 2, 2))
  expect_false(r$converged)
  expect_match(r$message, "curvature is not finite")
  expect_identical(r$par, c(0, 0))

  # the error after the first step leaves the polish at the accepted point
  n <- 0L
  C_err <- function(x) {
    n <<- n + 1L
    if (n > 1L) stop("cannot allocate the Hessian")
    matrix(c(4, 1, 1, 3), 2)
  }
  r <- run_polish(c(0, 0), q$f, C_err)
  expect_false(r$converged)
  expect_identical(r$error, "cannot allocate the Hessian")
  expect_match(r$message, "^an error: cannot allocate the Hessian")
  expect_identical(r$iterations, 1L)
  expect_equal(r$par, solve(matrix(c(4, 1, 1, 3), 2), c(1, 2)),
               tolerance = 1e-14)
})

test_that("a held parameter's gradient does not widen the active margin", {
  # x1 is held far from its minimum (gradient -10); x2 starts 5e-4 above its
  # lower bound, its minimum at 2e-4, curvature 1.5: its margin comes from
  # its own projected gradient (4.5e-4), so it stays free and one Newton step
  # reaches 2e-4. A margin of 1e-3, set by the held parameter, would make it
  # epsilon-active, its scaled-gradient move overshooting to 5e-5
  w <- c(1, 1.5)
  f <- function(x) {
    pol_obj(0.5 * sum(w * (x - c(10, 2e-4))^2), w * (x - c(10, 2e-4)))
  }
  r <- run_polish(c(0, 5e-4), f, function(x) diag(w), lower = c(-Inf, 0),
                  fixed = 1L)
  expect_true(r$converged)
  expect_identical(r$iterations, 1L)
  expect_equal(r$par, c(0, 2e-4), tolerance = 1e-12)
})

test_that("held parameters stay put and the rest is certified given them", {
  A <- matrix(c(4, 1, 1, 3), 2)
  q <- pol_quad(A, c(1, 2))
  r <- run_polish(c(a = 0.5, b = 0), q$f, q$C, fixed = 1L)
  expect_true(r$converged)
  expect_identical(r$held, "a")
  expect_identical(unname(r$par[1]), 0.5)
  expect_equal(unname(r$par[2]), (2 - 0.5) / 3, tolerance = 1e-14)
})

test_that("a flat parameter is left alone and stops the certificate", {
  # x3 does not enter the objective: its column of the curvature and its
  # gradient are zero, which would leave the undamped system singular. It is
  # left out of the steps; the others reach their minimum, but the polish
  # does not certify a parameter that carries no information
  A <- matrix(c(4, 1, 0, 1, 3, 0, 0, 0, 0), 3)
  q <- pol_quad(A, c(1, 2, 0))
  r <- run_polish(c(a = 0, b = 0, c = 7), q$f, q$C)
  expect_false(r$converged)
  expect_identical(r$flat, "c")
  expect_match(r$message, "below tol_decrement, but c has no usable curvature")
  expect_match(r$message, "(no information)", fixed = TRUE)
  expect_identical(r$lambda, 0)
  expect_identical(unname(r$par[3]), 7)
  expect_equal(unname(r$par[1:2]), solve(A[1:2, 1:2], c(1, 2)),
               tolerance = 1e-14)
})

test_that("a step longer than 10 units of the coordinates is shortened", {
  # A curvature 12.5 times too small in x1 asks for a step of 25, 50 times
  # its scale 0.5 (between the maximum step and the 100 at which it would be
  # held): the first trial point moves x1 by 10 times its scale, the step's
  # direction kept
  trials <- list()
  f <- function(x) {
    trials[[length(trials) + 1L]] <<- x
    pol_obj(sum((x - c(1, 2))^2), 2 * (x - c(1, 2)))
  }
  C_bad <- function(x) diag(c(0.08, 2))
  r <- run_polish(c(0, 0), f, C_bad, scale = c(0.5, 1),
                  control = .polish_control(list(max_iter = 1L)))
  first <- trials[[2L]]
  expect_equal(abs(first[1]) / 0.5, 10, tolerance = 1e-14)
  # the direction is the Newton step's: x2's step is 2, x1's 25
  expect_equal(first[2] / first[1], 2 / 25, tolerance = 1e-10)
  expect_identical(r$iterations, 1L)
  expect_length(r$flat, 0L)
})

test_that("a step past ten maximum steps on a vanished curvature is held", {
  # x1's curvature, 1e6 times too small, asks for a step of 1e6, 2e6 times
  # its scale: held for the iteration, it does not shorten x2's step, which
  # reaches x2's minimum at once
  f <- function(x) pol_obj(sum((x - c(1, 2))^2), 2 * (x - c(1, 2)))
  r <- run_polish(c(a = 0, b = 0), f, function(x) diag(c(2e-6, 2)),
                  scale = c(0.5, 1))
  expect_false(r$converged)
  expect_identical(r$flat, "a")
  expect_identical(unname(r$par[1]), 0)
  expect_equal(unname(r$par[2]), 2, tolerance = 1e-14)
  expect_identical(r$iterations, 1L)
  expect_match(r$message, "a has no usable curvature")
})

test_that("a step that is not finite goes to damping, not round the stall test", {
  # A solver that returns NaN, undamped, for a parameter on a vanished
  # curvature: the stall test is NA there, and the step is damped instead
  f <- function(x) pol_obj(sum((x - c(1, 2))^2), 2 * (x - c(1, 2)))
  calls <- 0L
  local_mocked_bindings(newton_step = function(H, g, D, free_idx, lambda) {
    calls <<- calls + 1L
    if (calls > 50L) stop("the stall test looped")
    step <- numeric(length(g))
    step[free_idx] <- -g[free_idx] / diag(H)[free_idx]
    if (lambda == 0 && 1L %in% free_idx) step[1L] <- NaN
    list(step = step, info = 0L)
  })
  r <- run_polish(c(a = 0, b = 0), f, function(x) diag(c(2e-6, 2)),
                  scale = c(0.5, 1),
                  control = .polish_control(list(max_iter = 1L)))
  expect_null(r$error)
  expect_lt(calls, 50L)
  expect_gt(r$damped, 0L)
})

test_that("a vanished curvature with a short step is not held", {
  # x1's curvature is 1e6 times below x2's, but its Newton step, 5, is short
  # of ten maximum steps: it takes it, and the minimum is certified
  f <- function(x) {
    pol_obj(0.5e-6 * (x[1] - 5)^2 + 0.5 * (x[2] - 1)^2,
            c(1e-6 * (x[1] - 5), x[2] - 1))
  }
  r <- run_polish(c(a = 0, b = 0), f, function(x) diag(c(1e-6, 1)))
  expect_true(r$converged)
  expect_length(r$flat, 0L)
  expect_equal(unname(r$par), c(5, 1), tolerance = 1e-12)
})

test_that("a start far from the minimum is shortened, not held", {
  # Steps of a thousand times the scale on a sound curvature: each step is
  # shortened to the maximum step, its direction kept, until the minimum is
  # within reach
  A <- matrix(c(1, 0.5, 0.5, 1), 2)
  q <- pol_quad(A, c(1, 2))
  xs <- solve(A, c(1, 2))
  r <- run_polish(c(p = xs[1] + 1000, q = xs[2] - 700), q$f, q$C,
                  control = .polish_control(list(max_iter = 200L)))
  expect_true(r$converged)
  expect_length(r$flat, 0L)
  expect_gt(r$iterations, 50L)
  expect_equal(unname(r$par), xs, tolerance = 1e-12)
})

test_that("with every parameter flat, no decrement is claimed", {
  f <- function(x) pol_obj(1, c(0, 0))
  r <- run_polish(c(a = 1, b = 2), f, function(x) matrix(0, 2, 2))
  expect_false(r$converged)
  expect_identical(r$flat, c("a", "b"))
  expect_identical(r$message,
                   "a, b have no usable curvature at this point (no information)")
})

test_that("the last curvature comes back when it was evaluated at the point returned", {
  # A curvature that records where it was evaluated: the result's curvature
  # must be the one at `par`, for a caller (the fit's variance) that would
  # otherwise evaluate it there again
  tag <- function(curv) function(x) structure(curv(x), at = x)
  rp <- function(...) run_polish(..., keep_curvature = TRUE)
  here <- function(r, C, what) {
    expect_identical(r$curvature_at_par, C(r$par), label = what)
  }
  q <- pol_quad(matrix(c(4, 1, 1, 3), 2), c(1, 2))
  Cq <- tag(q$C)
  Cr <- tag(pol_rosen$C)

  r <- rp(c(a = 5, b = -3), q$f, Cq)
  expect_true(r$converged)
  here(r, Cq, "converged")
  # it is released with the polish unless the caller keeps it
  expect_null(run_polish(c(a = 5, b = -3), q$f, Cq)$curvature_at_par)
  r <- rp(c(0, 1), pol_rosen$f, Cr,
          control = .polish_control(list(max_iter = 2L)))
  expect_match(r$message, "iteration limit")
  expect_identical(r$iterations, 2L)
  here(r, Cr, "the iteration limit after steps")
  r <- rp(c(0, 1), pol_rosen$f, Cr,
          control = .polish_control(list(max_iter = 0L)))
  expect_identical(r$par, c(0, 1))
  here(r, Cr, "max_iter = 0")
  r <- rp(c(0, 0), q$f, Cq,
          control = .polish_control(list(tol_step = 1e3)))
  expect_match(r$message, "step is below tol_step")
  here(r, Cq, "a step below tol_step")
  # every point but the start overflows: the line search fails at the start
  x0 <- c(1, 1)
  f_over <- function(x) {
    if (!identical(x, x0)) return(pol_obj(1e10, c(0, 0), overflow = TRUE))
    pol_obj(1, c(1, 1))
  }
  C1 <- tag(function(x) diag(2))
  r <- rp(x0, f_over, C1)
  expect_identical(r$message, "the line search failed")
  here(r, C1, "a failed line search")
  # a flat parameter stops the certificate where its curvature was evaluated
  A3 <- matrix(c(4, 1, 0, 1, 3, 0, 0, 0, 0), 3)
  q3 <- pol_quad(A3, c(1, 2, 0))
  C3 <- tag(q3$C)
  r <- rp(c(a = 0, b = 0, c = 7), q3$f, C3)
  expect_identical(r$flat, "c")
  here(r, C3, "a flat parameter")
  # a curvature that is not finite is returned as evaluated
  Cnan <- tag(function(x) matrix(NaN, 2, 2))
  r <- rp(c(0, 0), q$f, Cnan)
  expect_match(r$message, "curvature is not finite")
  here(r, Cnan, "a curvature that is not finite")
  # an error in a trial evaluation ends at the last accepted point, whose
  # curvature was evaluated before the step was tried: with the curvature
  # twice too large each step goes halfway and is accepted, so the fourth
  # evaluation is the third step's trial
  C2A <- tag(function(x) 2 * matrix(c(4, 1, 1, 3), 2))
  k <- 0L
  f_err <- function(x) {
    k <<- k + 1L
    if (k == 4L) stop("cannot allocate the gradient")
    q$f(x)
  }
  r <- rp(c(0, 0), f_err, C2A)
  expect_identical(r$error, "cannot allocate the gradient")
  expect_identical(r$iterations, 2L)
  here(r, C2A, "an error in a trial evaluation")

  # No curvature at the point returned: the start returned by the final
  # guard (the curvature was evaluated after the steps it abandons) ...
  C0 <- 1e9
  n <- 0L
  f_up <- function(x) {
    n <<- n + 1L
    pol_obj(C0 + 0.5 * (x - 1)^2 + 0.06 * (n - 1L), x - 1)
  }
  C2 <- tag(function(x) matrix(2))
  r <- rp(1 - 0.005, f_up, C2,
          control = .polish_control(list(max_iter = 3L)))
  expect_match(r$message, "worse than the start, which is returned")
  expect_null(r$curvature_at_par)
  # ... an error in the curvature itself, at the start or after a step ...
  C_err <- function(at_call) {
    m <- 0L
    function(x) {
      m <<- m + 1L
      if (m == at_call) stop("cannot allocate the Hessian")
      Cq(x)
    }
  }
  r <- rp(c(0, 0), q$f, C_err(1L))
  expect_identical(r$error, "cannot allocate the Hessian")
  expect_identical(r$iterations, 0L)
  expect_null(r$curvature_at_par)
  r <- rp(c(0, 0), q$f, C_err(2L))
  expect_identical(r$iterations, 1L)
  expect_null(r$curvature_at_par)
  # ... and a start whose objective is not finite
  r <- rp(c(0, 0), function(x) pol_obj(Inf, c(0, 0)), Cq)
  expect_match(r$message, "objective is not finite at the start")
  expect_null(r$curvature_at_par)
})

test_that(".polish_held() holds the constants without an estimate and checks BHHH's rank", {
  # rows 2 to 4 of alt_mapping have the constants at positions 3 to 5: row 3
  # is never chosen, row 4 chosen wherever offered
  ub <- list(never = 3L, always = 4L)
  expect_identical(.polish_held("newton", ub, 3:5, 5L, 3L), c(4L, 5L))
  expect_identical(.polish_held("none", ub, 3:5, 5L, 3L), integer())
  expect_identical(.polish_held("newton", NULL, 3:5, 5L, 3L), integer())
  # the units are counted only for "bhhh"
  expect_identical(.polish_held("newton", ub, 3:5, 5L, stop("counted")), c(4L, 5L))
  # without bounds (NULL) every parameter not held moves
  expect_error(.polish_held("bhhh", NULL, integer(), 5L, 3L),
               "the BHHH matrix of 3 choice situations has rank at most 3, below the 5 parameters",
               fixed = TRUE)
  expect_identical(.polish_held("bhhh", NULL, integer(), 5L, 5L), integer())
  # held and pinned parameters do not move
  expect_identical(.polish_held("bhhh", ub, 3:5, 5L, 3L), c(4L, 5L))
  expect_identical(.polish_held("bhhh", NULL, integer(), 5L, 4L,
                                lower = c(0, rep(-Inf, 4)),
                                upper = c(0, rep(Inf, 4))), integer())
  expect_error(.polish_held("bhhh", NULL, integer(), 5L, 2L, units = "decision makers",
                            lower = c(0, rep(-Inf, 4))),
               "the BHHH matrix of 2 decision makers has rank at most 2, below the 5 parameters",
               fixed = TRUE)
})

test_that(".polish_fit() certifies a rise only where held or bounded constants anchor the rest", {
  # three constants, the quadratic's minimum at (1, 2, 3), for alternatives
  # 1 to 3; the outside option (row 1) never chosen. Situations: {1, 2, 3}
  # choosing 1, {1, 2, 3} choosing 2, {2, 3} choosing 3
  q <- pol_quad(diag(3), c(1, 2, 3))
  am <- data.frame(alt_int = 0:3, N_CHOICES = c(0, 1, 1, 1))
  ub <- list(never = integer(), always = integer(), ref = 1L, ref_never = TRUE,
             ref_always = FALSE, outside = TRUE)
  linked <- list(alt_idx = c(1L, 2L, 3L, 1L, 2L, 3L, 2L, 3L),
                 choice_idx = c(1L, 2L, 2L), M = c(3L, 3L, 2L))
  opt <- list(par = c(0, 0, 0), value = 0, convergence = 4L, message = "xtol",
              iterations = 7L)
  t_opt <- system.time(NULL)
  pf <- function(held, upper = NULL, unbounded = ub, situations = linked,
                 polish = "newton", ref_open = NULL) {
    suppressMessages(.polish_fit(
      polish, opt, t_opt, c("a", "b", "c"), q$f, q$C, lower = NULL,
      upper = upper, scale = rep(1, 3), control = .polish_control(NULL),
      held = held, keep_curvature = FALSE, unbounded = unbounded,
      asc_idx = 1:3, alt_mapping = am,
      anchors = .polish_anchors("newton", unbounded, held, upper, 1:3, am,
                                situations),
      weights = rep(2, 4), optimizer = NULL, ref_open = ref_open))
  }
  # a held constant whose alternative is chosen where the others are
  # offered: the rise stops, and the rest is certified given it
  r <- pf(1L)
  expect_true(r$record$converged)
  expect_identical(r$opt$convergence, 1L)
  expect_identical(r$opt$par, c(0, 2, 3))
  expect_identical(r$record$held, "a")
  expect_identical(r$record$elapsed_time, r$time[["elapsed"]])
  expect_identical(r$record$optimizer$elapsed_time, t_opt[["elapsed"]])
  # neither held nor bounded, without bounds or with infinite ones: no
  # maximum to certify
  for (up in list(NULL, rep(Inf, 3))) {
    r <- pf(integer(), up)
    expect_false(r$record$converged)
    expect_identical(r$opt$convergence, 4L)
    expect_match(r$record$message, paste(
      "the constants diverge together (the outside option is never chosen),",
      "so there is no maximum to certify"), fixed = TRUE)
  }
  # a constant bounded above stops it too
  expect_true(pf(integer(), c(Inf, 5, Inf))$record$converged)
  # the held alternative chosen only where it is offered alone stops nothing
  apart <- list(alt_idx = c(1L, 2L, 3L, 2L, 3L), choice_idx = c(1L, 1L, 2L),
                M = c(1L, 2L, 2L))
  r <- pf(1L, situations = apart)
  expect_false(r$record$converged)
  expect_match(r$record$message, "the constants diverge together", fixed = TRUE)
  # a caller's ref_open (the nested logit's never-chosen base in a nest whose
  # dissimilarity exceeds 1) sets the rule aside
  expect_true(pf(1L, situations = apart, ref_open = function(par) TRUE)$record$converged)
  # a base alternative chosen wherever offered: the constants fall together,
  # never certified
  ub_fall <- utils::modifyList(ub, list(ref_never = FALSE, ref_always = TRUE,
                                        outside = FALSE))
  r <- pf(1L, unbounded = ub_fall, ref_open = function(par) TRUE)
  expect_false(r$record$converged)
  expect_match(r$record$message, paste(
    "the constants diverge together (the base alternative is chosen",
    "wherever it is offered)"), fixed = TRUE)
  expect_error(pf(1L, polish = "none"))
})

test_that(".polish_anchors() anchors through the situations' choices", {
  am <- data.frame(alt_int = 0:4)
  ub <- list(never = integer(), always = integer(), ref = 1L, ref_never = TRUE,
             ref_always = FALSE, outside = TRUE)
  anc <- function(held, sit, upper = NULL, polish = "newton", unbounded = ub,
                  map = am) {
    .polish_anchors(polish, unbounded, held, upper, seq_len(nrow(map) - 1L), map, sit)
  }
  # {1, 2} choosing 1, {2, 3} choosing 2, {3, 4} choosing 3: anchored from 1
  # in three steps
  chain <- list(alt_idx = c(1L, 2L, 2L, 3L, 3L, 4L), choice_idx = c(1L, 1L, 1L),
                M = c(2L, 2L, 2L))
  expect_identical(anc(1L, chain), rep(TRUE, 4))
  # from 4, which no situation chooses, nothing else
  expect_identical(anc(4L, chain), c(FALSE, FALSE, FALSE, TRUE))
  # from 3: alternative 1 stays free
  expect_identical(anc(3L, chain), c(FALSE, FALSE, TRUE, TRUE))
  # a constant bounded above anchors as a held one does
  expect_identical(anc(integer(), chain, upper = c(2, Inf, Inf, Inf)), rep(TRUE, 4))
  # a situation that chooses the outside option (choice 0; never, where the
  # rule applies) is skipped
  out <- list(alt_idx = c(1L, 2L, 2L, 3L), choice_idx = c(1L, 0L), M = c(2L, 2L))
  expect_identical(anc(1L, out, map = data.frame(alt_int = 0:3)), c(TRUE, TRUE, FALSE))
  # needed only for a polish whose reference is never chosen
  expect_null(anc(1L, chain, polish = "none"))
  expect_null(anc(1L, chain, unbounded = NULL))
  expect_null(anc(1L, chain, unbounded = utils::modifyList(ub, list(ref_never = FALSE,
                                                                    ref_always = TRUE))))
})
