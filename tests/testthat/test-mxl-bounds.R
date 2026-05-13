# Tests for run_mxlogit(lower=, upper=): parameter-bound plumbing.
#
# The optimizer accepts bounds in three forms:
#   (1) NULL                                  -> unbounded (default).
#   (2) Unnamed length-n_params numeric       -> nloptr-native, used as-is.
#   (3) Named (partial) numeric               -> names must be a subset of
#                                                 param_names; missing entries
#                                                 default to +/- Inf.
# Bounds are interpreted in natural-scale units and forward-transformed
# internally when scale_vars != "none". These tests run with scale_vars="none"
# so we can check the bound is honored against the natural-units value the user
# passed in directly.

build_inputs <- function(seed = 2026) {
  # outside_option=FALSE keeps the parameter count tidy: J=3 inside alts
  # contribute exactly J-1 = 2 free ASCs (the first is normalized to 0).
  sim <- simulate_mxl_data(N = 400L, J = 3L, seed = seed, outside_option = FALSE)
  dt <- data.table::as.data.table(sim$data)
  list(
    dt = dt,
    common = list(
      data = dt,
      id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"),
      random_var_cols = c("w1", "w2"),
      S = 30L,
      rc_mean = FALSE,
      rc_correlation = FALSE,
      use_asc = TRUE,
      scale_vars = "none",
      control = list(maxeval = 200L)
    )
  )
}

# K_x = 2 (x1, x2); K_w = 2 (w1, w2); rc_mean=FALSE -> no mu block;
# rc_correlation=FALSE -> L block has K_w=2 diagonal entries (L_11, L_22);
# J=3 inside alts -> n_asc = J - 1 = 2 ASCs.
# n_params = 2 + 0 + 2 + 2 = 6.
# Index order: x1, x2, L_11, L_22, ASC_2, ASC_3.

test_that("full-length unnamed bounds clip the Cholesky diagonal", {
  ctx <- build_inputs()
  n_params <- 6L

  lo <- rep(-Inf, n_params); up <- rep(Inf, n_params)
  # Positions 3, 4 are the Cholesky diagonal entries L_11, L_22.
  lo[3:4] <- -0.3
  up[3:4] <-  0.3
  # Start inside the box; the default theta_init places log(0.5) ~= -0.693
  # on the diagonal which is outside the lower bound.
  init <- rep(0, n_params)

  fit <- suppressMessages(do.call(run_mxlogit, c(
    ctx$common,
    list(lower = lo, upper = up, theta_init = init)
  )))

  estimates <- coef(fit)
  diag_vals <- estimates[c("L_11", "L_22")]
  slack <- 1e-8
  expect_true(all(diag_vals >= -0.3 - slack))
  expect_true(all(diag_vals <=  0.3 + slack))
})

test_that("named partial bounds clip only the named parameter", {
  ctx <- build_inputs()
  n_params <- 6L

  fit <- suppressMessages(do.call(run_mxlogit, c(
    ctx$common,
    list(
      lower = c(L_11 = -0.3),
      upper = c(L_11 =  0.3),
      # Default theta_init has log(0.5) on the diagonal which would be outside
      # this tight a lower bound, so cold-start at zero (which is inside).
      theta_init = rep(0, n_params)
    )
  )))

  slack <- 1e-8
  l11 <- coef(fit)[["L_11"]]
  expect_true(l11 >= -0.3 - slack)
  expect_true(l11 <=  0.3 + slack)
})

test_that("unnamed bound vector of wrong length errors informatively", {
  ctx <- build_inputs()
  bad <- rep(0, 99L)  # n_params is 6, so 99 is clearly wrong.
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(ctx$common, list(lower = bad)))),
    "must have length"
  )
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(ctx$common, list(upper = bad)))),
    "must have length"
  )
})

test_that("named bound vector with unknown name errors and names the offender", {
  ctx <- build_inputs()
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(
      ctx$common,
      list(lower = c(bogus = -1))
    ))),
    "bogus"
  )
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(
      ctx$common,
      list(upper = c(bogus = 1))
    ))),
    "bogus"
  )
})

test_that("bound vector with NA/NaN errors", {
  ctx <- build_inputs()
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(
      ctx$common, list(lower = c(L_11 = NA_real_))
    ))),
    "NA/NaN"
  )
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(
      ctx$common, list(upper = c(L_11 = NaN))
    ))),
    "NA/NaN"
  )
})

test_that("bound vector with duplicate names errors", {
  ctx <- build_inputs()
  expect_error(
    suppressMessages(do.call(run_mxlogit, c(
      ctx$common, list(lower = c(L_11 = -1, L_11 = -2))
    ))),
    "Duplicate"
  )
})
