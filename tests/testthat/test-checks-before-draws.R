# Inputs that can fail are checked before the expensive work. run_mxlogit()
# checks a store-mode cube's size, then computes the column scales of
# scale_vars = "sd", "mad" and "iqr" (with their check) and checks a supplied
# start, all before it builds the cube; the three wrappers check that every
# parameter has a name before the optimizer runs (they used to fail when
# naming the estimates, after it).

cbd_data <- function(seed = 11L) {
  set.seed(seed)
  dt <- data.table::data.table(id = rep(1:120, each = 3L), alt = rep(1:3, 120L))
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  # A mostly-zero column: its median absolute deviation is 0
  dt[, w0 := ifelse(stats::runif(.N) < 0.8, 0, stats::rnorm(.N))]
  dt[, nest := c("A", "A", "B")[alt]]
  dt[, V := 0.8 * x1 - 0.5 * x2 + 0.3 * w0]
  dt[, choice := as.integer(alt == sample(alt, 1L, prob = exp(V) / sum(exp(V)))),
     by = id]
  dt[, V := NULL]
  dt[]
}

# Count the likelihood evaluations of the three models until the test ends
cbd_count_evals <- function(env = parent.frame()) {
  counter <- new.env()
  counter$n <- 0L
  count <- function(name) {
    real <- get(name, envir = asNamespace("choicer"))
    function(...) {
      counter$n <- counter$n + 1L
      real(...)
    }
  }
  local_mocked_bindings(
    mnl_loglik_gradient_parallel = count("mnl_loglik_gradient_parallel"),
    mxl_loglik_gradient_parallel = count("mxl_loglik_gradient_parallel"),
    nl_loglik_gradient_parallel = count("nl_loglik_gradient_parallel"),
    .env = env
  )
  counter
}

test_that("run_mxlogit() checks the column scales before it builds store-mode draws", {
  dt <- cbd_data()
  local_mocked_bindings(get_halton_normals = function(...) stop("draws built"))
  expect_error(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = c("x2", "w0"), S = 10L,
      draws = "store", scale_vars = "mad")),
    "scale_vars='mad': normal random-coefficient column\\(s\\) with scale"
  )
  # The same call with scales that pass reaches the draws
  expect_error(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = c("x2", "w0"), S = 10L,
      draws = "store", scale_vars = "sd")),
    "draws built"
  )
})

test_that("run_mxlogit() checks a supplied start before it builds store-mode draws", {
  dt <- cbd_data()
  local_mocked_bindings(get_halton_normals = function(...) stop("draws built"))
  for (sv in c("bhhh", "sd")) {
    expect_error(
      suppressMessages(run_mxlogit(
        data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
        covariate_cols = "x1", random_var_cols = "x2", S = 10L,
        draws = "store", scale_vars = sv, theta_init = c(0, 0))),
      "`theta_init` must have one value per parameter", info = sv
    )
    expect_error(
      suppressMessages(run_mxlogit(
        data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
        covariate_cols = "x1", random_var_cols = "x2", S = 10L,
        draws = "store", scale_vars = sv, theta_init = c(0, NA, 0, 0))),
      "`theta_init` must be finite", info = sv
    )
  }
  # Under "none" the kernels check the start's length, as before
  expect_error(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = "x2", S = 10L,
      draws = "store", scale_vars = "none", theta_init = c(0, 0))),
    "draws built"
  )
})

test_that("run_mxlogit() checks a store-mode cube's size before the column scales", {
  dt <- cbd_data()
  local_mocked_bindings(.halton_points_ok = function(S, N) FALSE,
                        .column_scales = function(...) stop("scales computed"))
  expect_error(
    suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = "x2", S = 10L,
      draws = "store", scale_vars = "sd")),
    "A store-mode draw cube needs S \\* N"
  )
})

test_that("every parameter needs a name before the optimizer runs", {
  dt <- cbd_data()
  counter <- cbd_count_evals()

  inp <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
  colnames(inp$X) <- NULL
  expect_error(suppressMessages(run_mnlogit(input_data = inp)),
               "`input_data\\$X` needs one column name per covariate \\(2 columns, 0 names\\)")
  expect_identical(counter$n, 0L)

  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  colnames(inp$X) <- NULL
  eta <- get_halton_normals(10L, length(inp$M), 1L)
  expect_error(suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta)),
               "`input_data\\$X` and `input_data\\$W` need one column name")
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("x2", "w0"))
  colnames(inp$W) <- NULL
  eta <- get_halton_normals(10L, length(inp$M), 2L)
  expect_error(suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta,
                                            rc_mean = TRUE)),
               "`input_data\\$X` and `input_data\\$W` need one column name")
  # One unnamed W column used to name its mean "Mu_"
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  colnames(inp$W) <- NULL
  eta <- get_halton_normals(10L, length(inp$M), 1L)
  expect_error(suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta,
                                            rc_mean = TRUE)),
               "`input_data\\$X` and `input_data\\$W` need one column name")
  expect_identical(counter$n, 0L)

  expect_error(
    suppressMessages(run_nestlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), nest_col = "nest",
      param_names = c("b1", "b2"))),
    "`param_names` has 2 names; the model has 5 parameters"
  )
  expect_identical(counter$n, 0L)
})

test_that("an empty block of parameters has no names", {
  # Every nest a single alternative: no dissimilarity parameter, so no
  # "Lambda_" name, and the kernel's own error comes first, as before
  dt <- cbd_data()
  dt[, nest1 := paste0("n", alt)]
  expect_error(
    suppressMessages(run_nestlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), nest_col = "nest1")),
    "No non-singleton nests found"
  )
})
