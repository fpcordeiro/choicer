# A log-Cholesky diagonal gone to minus infinity: the polish's decrement cannot
# tell a false stop, where the draws' asymmetry makes a zero variance a local
# maximum below a higher one, from a variance whose maximum is at zero. The
# polish restarts once from a larger standard deviation, the higher likelihood
# stands, and a variance that ends at zero warns (choicer_variance_plateau).

PLATEAU <- "choicer_variance_plateau"

# Random coefficient on w, true SD `sd`; a fixed one on x
plat_data <- function(seed = 1L, N = 300L, J = 4L, sd = 1.5) {
  set.seed(seed)
  id <- rep(seq_len(N), each = J)
  w <- rnorm(N * J)
  w2 <- rnorm(N * J)
  x <- rnorm(N * J)
  b <- rep(rnorm(N, 1, sd), each = J)
  u <- b * w - 0.5 * x - 0.5 * w2 - log(-log(runif(N * J)))
  data.table(id = id, alt = rep(seq_len(J), N), w = w, w2 = w2, x = x,
             person = (id - 1L) %/% 3L,
             choice = as.integer(u == ave(u, id, FUN = max)))
}

# Logit data: every true variance is zero
plat_zero <- function() {
  set.seed(3)
  N <- 1500L
  J <- 4L
  id <- rep(seq_len(N), each = J)
  w <- rnorm(N * J)
  w2 <- rnorm(N * J)
  x <- rnorm(N * J)
  u <- w - 0.5 * w2 + 0.5 * x - log(-log(runif(N * J)))
  data.table(id = id, alt = rep(seq_len(J), N), w = w, w2 = w2, x = x,
             choice = as.integer(u == ave(u, id, FUN = max)))
}

# The MNL's estimates as an MXL start on w, with L at `L`
plat_start <- function(d, L = -16) {
  cf <- coef(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "w"))))
  c(cf[["x"]], cf[["w"]], L, unname(cf[grep("^ASC", names(cf))]))
}

plat_fit <- function(d, ...) {
  suppressMessages(choicer::run_mxlogit(d, "id", "alt", "choice", "x", "w", rc_mean = TRUE,
                                        S = 50L, ...))
}

test_that(".mxl_variance_plateau() reads the default start's level, not pinned parameters", {
  plat <- choicer:::.mxl_variance_plateau
  lev <- c(0, 2, log(0.5), 0)
  idx <- 2:3
  lo <- rep(-Inf, 4L)
  up <- rep(Inf, 4L)
  # below 1e-4 of the level (a column of step exp(-2); a log-normal's 0.5)
  expect_identical(plat(c(9, 2 + log(1e-5), log(0.5), 9), lev, idx, lo, up), 2L)
  expect_identical(plat(c(9, 2 + log(1e-3), log(0.5) + log(1e-5), 9), lev, idx, lo, up), 3L)
  expect_identical(plat(c(9, 2 + log(1e-3), log(0.5) + log(1e-3), 9), lev, idx, lo, up), integer())
  # a parameter pinned there, or bounded above below the line, is the user's
  expect_identical(plat(c(9, -30, -30, 9), lev, idx, replace(lo, 2L, -30), replace(up, 2L, -30)), 3L)
  expect_identical(plat(c(9, -30, -30, 9), lev, idx, lo, replace(up, 2L, -20)), 3L)
  # a point that is not finite is not on it
  expect_identical(plat(c(9, NaN, NaN, 9), lev, idx, lo, up), integer())
})

test_that("a false plateau: the restart reaches the higher maximum, and nothing warns", {
  skip_on_cran()  # the optimizer's real path
  d <- plat_data()
  expect_no_warning(
    fit <- plat_fit(d, theta_init = plat_start(d), polish = "newton"), class = PLATEAU)
  r <- fit$optimizer$polish$restart
  expect_identical(r$parameters, "L_11")
  expect_true(r$adopted)
  expect_gt(r$gain, 1)
  expect_true(r$converged)
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$convergence, 1L)
  expect_gt(coef(fit)[["L_11"]], -2)
  expect_match(fit$message, "a restart from a larger standard deviation of L_11 raised", fixed = TRUE)
  expect_true(all(is.finite(fit$se)))
  # the maximum the default start reaches, on the same draws
  plain <- plat_fit(d, polish = "newton")
  expect_equal(fit$loglik, plain$loglik, tolerance = 1e-6)
  expect_null(plain$optimizer$polish$restart)
  # the kept point's variance is its own: the Hessian route evaluated afresh agrees
  expect_equal(fit$vcov, vcov(fit, type = "hessian"), tolerance = 1e-8)
})

test_that("a zero variance: the restart comes back, the first point stands, and it warns", {
  skip_on_cran()
  d <- plat_zero()
  w <- expect_warning(
    fit <- suppressMessages(choicer::run_mxlogit(d, "id", "alt", "choice", "x", c("w", "w2"),
                                                 rc_mean = TRUE, S = 100L, polish = "newton")),
    class = PLATEAU)
  expect_match(conditionMessage(w), "numerically zero", fixed = TRUE)
  expect_match(conditionMessage(w), "found no higher likelihood", fixed = TRUE)
  r <- fit$optimizer$polish$restart
  expect_false(r$adopted)
  expect_lt(abs(r$gain), 1e-6)
  expect_true(fit$optimizer$polish$converged)
  expect_identical(fit$convergence, 1L)
  expect_match(fit$message, "found no higher log-likelihood", fixed = TRUE)
  # the variance is evaluated at the first point, where it stands
  expect_equal(fit$vcov, vcov(fit, type = "hessian"), tolerance = 1e-8)
})

test_that("a higher point that is not certified stands, with its verdict", {
  skip_on_cran()
  d <- plat_data()
  # from the plateau, no optimizer steps and no Newton steps certify it at
  # once; the restart's L-BFGS climbs, stopped short, and no step follows
  fit <- plat_fit(d, theta_init = plat_start(d), polish = "newton", control = list(maxeval = 8L),
                  polish_control = list(max_iter = 0L))
  r <- fit$optimizer$polish$restart
  expect_true(r$adopted)
  expect_false(r$converged)
  expect_false(fit$optimizer$polish$converged)
  expect_false(choicer:::.is_converged(fit))
})

test_that("the warning says what the restart can say, and only that", {
  say <- function(names, restart, can_restart = TRUE) {
    tryCatch(choicer:::.warn_variance_plateau(names, rep(-20, length(names)), restart, can_restart),
             choicer_variance_plateau = function(w) conditionMessage(w))
  }
  came_back <- list(parameters = "L_11", adopted = FALSE, gain = -1e-9)
  expect_match(say("L_11", came_back), "found no higher likelihood", fixed = TRUE)
  rose <- list(parameters = "L_11", adopted = TRUE, gain = 4.93)
  m <- say("L_11", rose)
  expect_match(m, "raised the likelihood by 4.93 and ended here again", fixed = TRUE)
  expect_no_match(m, "found no higher", fixed = TRUE)
  # a diagonal the restart did not raise, gone there during it
  m <- say("L_22", rose)
  expect_match(m, "went there during a restart from a larger standard deviation of L_11", fixed = TRUE)
  expect_no_match(m, "found no higher", fixed = TRUE)
  m <- say("L_11", list(parameters = "L_11", adopted = FALSE, gain = NA_real_, error = "boom"))
  expect_match(m, "stopped with an error (boom)", fixed = TRUE)
  expect_no_match(m, "found no higher", fixed = TRUE)
  expect_match(say("L_11", NULL), "polish = \"newton\"", fixed = TRUE)
  expect_no_match(say("L_11", NULL, can_restart = FALSE), "polish = \"newton\"", fixed = TRUE)
})

test_that("the restart raises a normal's level and a log-normal's to 1, in fresh coordinates", {
  d <- plat_data()
  m <- coef(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "w", "w2"))))
  start <- c(m[["x"]], m[["w"]], m[["w2"]], -16, -16, unname(m[grep("^ASC", names(m))]))
  real <- choicer:::.optimize_in_coordinates
  real_map <- choicer:::.coordinate_map
  starts <- list()
  maps <- list()
  opg <- list()
  local_mocked_bindings(
    .optimize_in_coordinates = function(map, optimizer, theta_init, ...) {
      starts[[length(starts) + 1L]] <<- theta_init
      maps[[length(maps) + 1L]] <<- map
      real(map, optimizer, theta_init, ...)
    },
    .coordinate_map = function(..., opg_diag = NULL) {
      opg[[length(opg) + 1L]] <<- opg_diag
      real_map(..., opg_diag = opg_diag)
    }, .package = "choicer")
  fit <- suppressWarnings(suppressMessages(choicer::run_mxlogit(
    d, "id", "alt", "choice", "x", c("w", "w2"), rc_dist = c(0L, 1L), rc_mean = TRUE, S = 50L,
    theta_init = start, polish = "newton", control = list(maxeval = 1L),
    polish_control = list(max_iter = 0L))),
    classes = "choicer_variance_plateau")
  expect_length(starts, 2L)
  L <- fit$param_map$sigma
  expect_equal(starts[[2L]][L[2L]], 0)
  expect_equal(starts[[2L]][L[1L]], -log(choicer:::design_column_step(fit$data$W)[1L]))
  # the coordinates read at the restart's start, not the first's
  expect_length(opg, 2L)
  expect_false(identical(opg[[1L]], opg[[2L]]))
  # the fit carries those of the point that stands
  r <- fit$optimizer$polish$restart
  expect_identical(unname(fit$param_scale), maps[[if (r$adopted) 2L else 1L]]$scale)
})

test_that("a first point that stands keeps its coordinates, and the restart can be turned off", {
  skip_on_cran()
  d <- plat_zero()
  real <- choicer:::.optimize_in_coordinates
  maps <- list()
  local_mocked_bindings(.optimize_in_coordinates = function(map, ...) {
    maps[[length(maps) + 1L]] <<- map
    real(map, ...)
  }, .package = "choicer")
  fit_z <- function(...) {
    suppressMessages(choicer::run_mxlogit(d, "id", "alt", "choice", "x", c("w", "w2"),
                                          rc_mean = TRUE, S = 100L, polish = "newton", ...))
  }
  fit <- suppressWarnings(fit_z(), classes = PLATEAU)
  r <- fit$optimizer$polish$restart
  expect_false(r$adopted)
  expect_length(maps, 2L)
  expect_identical(unname(fit$param_scale), maps[[1L]]$scale)
  maps <- list()
  expect_warning(off <- fit_z(polish_control = list(restart = FALSE)),
                 "fit with polish", class = PLATEAU)
  expect_length(maps, 1L)
  expect_null(off$optimizer$polish$restart)
  # (two threads: the same point to rounding, which moves a diagonal at -18 more)
  expect_equal(off$loglik, fit$loglik, tolerance = 1e-8)
  expect_error(fit_z(polish_control = list(restart = NA)),
               "`polish_control$restart` must be TRUE or FALSE", fixed = TRUE)
  # the setting is the mixed logit's
  expect_error(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "w"), polish = "newton",
                                            polish_control = list(restart = FALSE))),
               "Unknown `polish_control`", fixed = TRUE)
})

test_that("an adopted BHHH polish's curvature serves its variance", {
  skip_on_cran()
  d <- plat_data()
  fit <- suppressWarnings(
    plat_fit(d, theta_init = plat_start(d), polish = "bhhh", se_method = "bhhh"),
    classes = "choicer_variance_plateau")
  r <- fit$optimizer$polish$restart
  expect_true(r$adopted)
  expect_equal(fit$vcov, vcov(fit, type = "bhhh"), tolerance = 1e-8)
})

test_that("an error in the restart leaves the first point", {
  d <- plat_data()
  start <- plat_start(d)
  real <- choicer:::.optimize_in_coordinates
  calls <- 0L
  local_mocked_bindings(.optimize_in_coordinates = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("a mocked failure")
    real(...)
  }, .package = "choicer")
  w <- expect_warning(fit <- plat_fit(d, theta_init = start, polish = "newton"), class = PLATEAU)
  expect_identical(calls, 2L)
  r <- fit$optimizer$polish$restart
  expect_false(r$adopted)
  expect_identical(r$error, "a mocked failure")
  expect_lt(coef(fit)[["L_11"]], -9)
  expect_true(fit$optimizer$polish$converged)
  expect_match(conditionMessage(w), "stopped with an error (a mocked failure)", fixed = TRUE)
  expect_match(fit$message, "stopped with an error: a mocked failure", fixed = TRUE)
  expect_no_match(fit$message, "found no higher", fixed = TRUE)
})

test_that("without a polish, and under the profile likelihood, the fit warns without restarting", {
  d <- plat_data()
  start <- plat_start(d)
  # the optimizer at the plateau (one evaluation), no polish
  w <- expect_warning(fit <- plat_fit(d, theta_init = start, control = list(maxeval = 1L)),
                      class = PLATEAU)
  expect_match(conditionMessage(w), "polish = \"newton\"", fixed = TRUE)
  expect_null(fit$optimizer$polish)
  expect_equal(coef(fit)[["L_11"]], -16, tolerance = 1e-12)

  w <- expect_warning(
    prof <- plat_fit(d, theta_init = start, polish = "newton", control = list(maxeval = 1L),
                     concentrate_out_asc = "profile_likelihood"),
    class = PLATEAU)
  expect_null(prof$optimizer$polish$restart)
  expect_match(conditionMessage(w), "start elsewhere (theta_init)", fixed = TRUE)
})

test_that("a diagonal pinned low is the user's, and an ordinary fit is not touched", {
  d <- plat_data()
  start <- plat_start(d)
  pinned <- expect_no_warning(
    plat_fit(d, theta_init = start, polish = "newton", lower = c(L_11 = -16), upper = c(L_11 = -16)),
    class = PLATEAU)
  expect_null(pinned$optimizer$polish$restart)
  expect_identical(coef(pinned)[["L_11"]], -16)
})

test_that("a failure code's check of no steps does not restart", {
  d <- plat_data()
  start <- plat_start(d)
  local_mocked_bindings(run_nloptr = function(theta_init, eval_f, lower, upper, control) {
    list(par = theta_init, value = eval_f(theta_init)$objective, convergence = -1L,
         message = "NLOPT_FAILURE: Generic failure code.", iterations = 1L, raw = list())
  }, .package = "choicer")
  fit <- suppressWarnings(plat_fit(d, theta_init = start), classes = PLATEAU)
  expect_identical(fit$optimizer$polish$iterations, 0L)
  expect_null(fit$optimizer$polish$restart)
})

test_that("the Cholesky diagonal is read under correlated coefficients", {
  skip_on_cran()
  d <- plat_data()
  m <- coef(suppressMessages(run_mnlogit(d, "id", "alt", "choice", c("x", "w", "w2"))))
  start <- c(m[["x"]], m[["w"]], m[["w2"]], 0.2, 0, -16, unname(m[grep("^ASC", names(m))]))
  fit <- suppressWarnings(suppressMessages(choicer::run_mxlogit(
    d, "id", "alt", "choice", "x", c("w", "w2"), rc_mean = TRUE, rc_correlation = TRUE, S = 50L,
    theta_init = start, control = list(maxeval = 1L), polish = "newton",
    polish_control = list(max_iter = 0L))), classes = PLATEAU)
  r <- fit$optimizer$polish$restart
  expect_identical(r$parameters, "L_22")
})
