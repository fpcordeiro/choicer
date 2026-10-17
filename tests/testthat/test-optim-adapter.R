# optimizer = "optim": stats::optim() takes the objective and its gradient
# from separate functions, while one kernel pass computes both. The adapter
# keeps the last point's evaluation (.one_point_cache()), so the gradient at
# the point whose objective optim has just computed costs no kernel call. The
# path without the cache (.one_point_cache() mocked as the identity) is the
# reference: optim must see the same numbers, bit for bit on one thread.

# The points evaluated, consecutive repeats dropped.
oc_collapse <- function(pts) {
  if (length(pts) < 2L) return(pts)
  repeated <- vapply(seq_along(pts)[-1L], function(i) {
    identical(pts[[i]], pts[[i - 1L]], num.eq = FALSE, single.NA = FALSE)
  }, TRUE)
  pts[c(TRUE, !repeated)]
}

# Run `code` without the cache.
oc_uncached <- function(code) {
  with_mocked_bindings(code, .one_point_cache = function(eval_f) eval_f)
}

# Rosenbrock's function in four dimensions; `pts` logs the points evaluated.
oc_rosenbrock <- function() {
  env <- new.env()
  env$pts <- list()
  env$f <- function(theta) {
    env$pts[[length(env$pts) + 1L]] <- theta
    i <- seq_len(length(theta) - 1L)
    d <- theta[i + 1L] - theta[i]^2
    e <- 1 - theta[i]
    g <- numeric(length(theta))
    g[i] <- -400 * theta[i] * d - 2 * e
    g[i + 1L] <- g[i + 1L] + 200 * d
    list(objective = sum(100 * d^2 + e^2), gradient = g)
  }
  env
}

test_that("the cache answers its last point and evaluates any other", {
  pts <- list()
  f <- function(theta) {
    pts[[length(pts) + 1L]] <<- theta
    list(objective = sum(theta^2), gradient = 2 * theta)
  }
  g <- choicer:::.one_point_cache(f)
  a <- c(1, -2)
  b <- c(1, 2)
  ra <- g(a)
  expect_identical(ra, list(objective = 5, gradient = c(2, -4)))
  expect_identical(g(a), ra)
  expect_identical(pts, list(a))
  g(b)
  g(a)  # an earlier point is evaluated afresh
  expect_identical(pts, list(a, b, a))
  g(c(0, 1))
  g(c(-0, 1))  # bit for bit: -0 is not 0
  expect_length(pts, 5L)
  g(c(u = -0, v = 1))  # names included
  expect_length(pts, 6L)
  # An evaluation that fails leaves the last one in place.
  h <- choicer:::.one_point_cache(function(theta) {
    if (theta[1] > 5) stop("not here")
    f(theta)
  })
  rh <- h(a)
  expect_error(h(c(6, 0)), "not here")
  expect_error(h(c(6, 0)), "not here")
  expect_identical(h(a), rh)
  expect_length(pts, 7L)
})

test_that("optim's methods see the uncached numbers, no point evaluated twice in a row", {
  cases <- list(
    BFGS = list(),
    `L-BFGS-B` = list(lower = rep(-2, 4), upper = rep(0.9, 4)),  # binding
    `L-BFGS-B free` = list(lower = rep(-5, 4), upper = rep(5, 4)),
    CG = list(method = "CG"),
    `Nelder-Mead` = list(method = "Nelder-Mead"))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    control <- list(maxit = 500L)
    control$method <- cs$method
    run <- function() {
      obj <- oc_rosenbrock()
      res <- choicer:::run_optimizer(
        "optim", c(-1.2, 1, -1.2, 1), obj$f, lower = cs$lower,
        upper = cs$upper, control = control)
      list(res = res, pts = obj$pts)
    }
    cached <- run()
    plain <- oc_uncached(run())
    expect_identical(cached$res, plain$res, label = nm)
    # The kernel sees the uncached path's points, less its repeats
    expect_identical(cached$pts, oc_collapse(plain$pts), label = nm)
    n <- cached$res$raw$counts
    if (nm %in% c("BFGS", "L-BFGS-B", "L-BFGS-B free")) {
      # Each call of fn and gr was an evaluation. Every gradient is asked for
      # at the point just evaluated, so at most one evaluation per objective
      # (fewer where L-BFGS-B asks for a point twice in a row).
      expect_length(plain$pts, n[["function"]] + n[["gradient"]])
      expect_lte(length(cached$pts), n[["function"]])
    }
    if (startsWith(nm, "L-BFGS-B")) {
      expect_identical(n[["function"]], n[["gradient"]], label = nm)
    }
    if (nm == "CG") {
      # Gradients back at earlier points, evaluated afresh
      expect_gt(length(cached$pts), n[["function"]])
    }
  }
})

test_that("the sentinel and the coordinates keep their numbers through the cache", {
  # Overflow beyond theta_1 = 2, the minimum at 3: the line searches step
  # into the overflow. The objective starts at 9.0e10, so the sentinel's
  # value there, ten times the largest objective seen, depends on the path.
  model <- function() {
    env <- new.env()
    env$pts <- list()
    env$f <- function(theta) {
      env$pts[[length(env$pts) + 1L]] <- theta
      if (theta[1] > 2) {
        return(list(objective = 1e10, gradient = c(0, 0), overflow = TRUE))
      }
      list(objective = 1e7 * (1e3 * (theta[1] - 3)^2 + (theta[2] + 1)^2),
           gradient = 1e7 * c(2e3 * (theta[1] - 3), 2 * (theta[2] + 1)),
           overflow = FALSE)
    }
    env
  }
  map <- list(scale = c(2^-4, 2), shift = c(0, 0))
  run <- function(...) {
    obj <- model()
    sentinel <- choicer:::.lift_sentinel(obj$f)
    lifted <- numeric()
    eval_f <- function(theta) {
      r <- sentinel(theta)
      if (theta[1] > 2) lifted <<- c(lifted, r$objective)
      r
    }
    res <- choicer:::.optimize_in_coordinates(
      map = map, optimizer = "optim", theta_init = c(0, 0), eval_f = eval_f,
      ...)
    list(res = res, pts = obj$pts, lifted = lifted)
  }
  for (b in list(list(), list(lower = c(-10, -10), upper = c(10, 10)))) {
    what <- if (length(b)) "L-BFGS-B" else "BFGS"
    cached <- do.call(run, b)
    plain <- oc_uncached(do.call(run, b))
    # The overflow was reached and its objective lifted above the kernel's
    expect_gt(length(cached$lifted), 0L, label = what)
    expect_true(all(cached$lifted > 1e10), label = what)
    expect_identical(cached$res, plain$res, label = what)
    expect_identical(cached$pts, oc_collapse(plain$pts), label = what)
    n <- cached$res$raw$counts
    expect_length(plain$pts, n[["function"]] + n[["gradient"]])
    expect_lte(length(cached$pts), n[["function"]])
    expect_lte(cached$res$par[1], 2)
  }
})

# A fit with its gradient kernel's calls counted and optim's counts kept;
# without the cache, the points optim asks for are logged at the cache's
# place (in the optimizer's coordinates).
oc_fit <- function(fit, kernel, cache = TRUE) {
  ns <- asNamespace("choicer")
  real_kernel <- get(kernel, envir = ns)
  real_run_optim <- get("run_optim", envir = ns)
  calls <- 0L
  counts <- NULL
  pts <- list()
  mocks <- stats::setNames(list(function(...) {
    calls <<- calls + 1L
    real_kernel(...)
  }), kernel)
  mocks$run_optim <- function(...) {
    r <- real_run_optim(...)
    counts <<- r$raw$counts
    r
  }
  if (!cache) {
    mocks$.one_point_cache <- function(eval_f) {
      force(eval_f)  # run_optim() rebinds the name to the result
      function(theta) {
        pts[[length(pts) + 1L]] <<- theta
        eval_f(theta)
      }
    }
  }
  local_mocked_bindings(!!!mocks)
  value <- suppressMessages(fit())
  value$optimizer$elapsed_time <- NULL
  list(fit = value, calls = calls, counts = counts, pts = pts)
}

# The cached fit against the uncached one: the same fit and optim counts.
# Kernel calls beyond the B pass: one per call of fn and gr without the
# cache; with it, one per point asked for, less immediate repeats.
oc_expect_once <- function(fit, kernel, bpass, what) {
  cached <- oc_fit(fit, kernel)
  plain <- oc_fit(fit, kernel, cache = FALSE)
  expect_identical(cached$fit, plain$fit, label = what)
  expect_identical(cached$counts, plain$counts, label = what)
  n <- plain$counts
  expect_length(plain$pts, n[["function"]] + n[["gradient"]])
  expect_identical(plain$calls, bpass + length(plain$pts), label = what)
  expect_identical(cached$calls, bpass + length(oc_collapse(plain$pts)),
                   label = what)
  expect_lte(cached$calls, bpass + n[["function"]], label = what)
  invisible(cached)
}

test_that("multinomial logit fits by BFGS are the uncached path's", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 61, outside_option = TRUE)
  for (sv in c("bhhh", "none", "sd")) {
    fit <- function() run_mnlogit(
      sim$data, "id", "alt", "choice", c("x1", "x2"), outside_opt_label = 0L,
      include_outside_option = TRUE, optimizer = "optim", scale_vars = sv)
    r <- oc_expect_once(fit, "mnl_loglik_gradient_parallel",
                        bpass = as.integer(sv == "bhhh"), what = sv)
    expect_identical(r$fit$convergence, 0L, label = sv)
  }
})

test_that("mixed logit fits by L-BFGS-B are the uncached path's at one pass per point", {
  skip_on_cran()
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_mxl_data(N = 120L, J = 3L, T = 2L, seed = 62L,
                           outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  cases <- list(
    list(),  # "bhhh", the default
    list(scale_vars = "none", rc_correlation = TRUE),
    list(scale_vars = "sd", person_col = "pid", draws = "generate", seed = 4L,
         upper = c(x1 = 0.5)))  # binding
  for (cs in cases) {
    what <- if (is.null(cs$scale_vars)) "bhhh" else cs$scale_vars
    fit <- function() do.call("run_mxlogit", c(list(
      dt, "id", "alt", "choice", c("x1", "x2"), c("w1", "w2"), S = 10L,
      optimizer = "optim"), cs))
    r <- oc_expect_once(fit, "mxl_loglik_gradient_parallel",
                        bpass = as.integer(what == "bhhh"), what = what)
    n <- r$counts
    expect_identical(n[["function"]], n[["gradient"]], label = what)
    if (!is.null(cs$upper)) expect_lte(coef(r$fit)[["x1"]], 0.5)
  }
})

test_that("nested logit fits by L-BFGS-B are the uncached path's at one pass per point", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  sim <- simulate_nl_data(N = 300, seed = 63)
  for (sv in c("none", "bhhh")) {
    fit <- function() run_nestlogit(
      sim$data, "id", "j", "choice", c("X", "W"), "nest", optimizer = "optim",
      scale_vars = sv)
    r <- oc_expect_once(fit, "nl_loglik_gradient_parallel",
                        bpass = as.integer(sv == "bhhh"), what = sv)
    expect_identical(r$counts[["function"]], r$counts[["gradient"]],
                     label = sv)
  }
})
