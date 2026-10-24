# Inputs that can fail are checked before the expensive work. run_mxlogit()
# checks a store-mode cube's size, then computes the column scales of
# scale_vars = "sd", "mad" and "iqr" (with their check) and checks a supplied
# start, all before it builds the cube; the three wrappers check that every
# parameter has a name before the optimizer runs (they used to fail when
# naming the estimates, after it). An eta_draws cube stored as integers
# reaches the mixed logit's kernels as doubles: they read the cube in place,
# and RcppArmadillo's conversion of an integer cube is a temporary that R may
# free while they read it; a cube of another type is an error. The three
# wrappers check the values of an input_data before any evaluation of the
# likelihood: that its weights are finite and strictly positive, with the
# preparations' messages (integer64 weights converted to double first), and
# that its X and W, and eta_draws, are finite (edited values used to leave a
# fit at its start values, reporting success).

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

test_that("an integer eta_draws cube reaches the kernels as doubles", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # so that the two fits agree bit for bit
  dt <- cbd_data()
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  eta <- get_halton_normals(20L, length(inp$M), 1L)
  eta_int <- array(as.integer(round(2 * eta)), dim(eta))
  eta_dbl <- array(as.double(eta_int), dim(eta))
  types <- character()
  real <- get("mxl_loglik_gradient_parallel", envir = asNamespace("choicer"))
  local_mocked_bindings(mxl_loglik_gradient_parallel = function(..., eta_draws) {
    types <<- c(types, typeof(eta_draws))
    real(..., eta_draws = eta_draws)
  })
  f_int <- suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta_int,
                                        se_method = "bhhh"))
  expect_true(length(types) > 0L)
  expect_true(all(types == "double"))
  f_dbl <- suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta_dbl,
                                        se_method = "bhhh"))
  expect_identical(coef(f_int), coef(f_dbl))
  expect_identical(f_int$loglik, f_dbl$loglik)
  expect_identical(f_int$vcov, f_dbl$vcov)
})

test_that("an eta_draws cube that is not numeric stops before the kernels", {
  dt <- cbd_data()
  counter <- cbd_count_evals()
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  eta <- get_halton_normals(10L, length(inp$M), 1L)
  eta_chr <- array(format(eta), dim(eta))
  expect_error(suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta_chr)),
               "`eta_draws` must be a numeric array")
  expect_identical(counter$n, 0L)
})

test_that("an integer64 eta_draws cube stops before the kernels", {
  skip_if_not_installed("bit64")
  dt <- cbd_data()
  counter <- cbd_count_evals()
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  eta <- get_halton_normals(10L, length(inp$M), 1L)
  eta_i64 <- bit64::as.integer64(round(eta))
  dim(eta_i64) <- dim(eta)
  expect_error(suppressMessages(run_mxlogit(input_data = inp, eta_draws = eta_i64)),
               "`eta_draws` must be a numeric array")
  expect_identical(counter$n, 0L)
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

test_that("input_data's weights must be finite and strictly positive", {
  dt <- cbd_data()
  counter <- cbd_count_evals()
  eta <- get_halton_normals(10L, 120L, 1L)
  preps <- list(
    mnl = prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2")),
    nl = prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest"),
    mxl = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  )
  fits <- list(
    mnl = function(d, ...) run_mnlogit(input_data = d, ...),
    nl = function(d, ...) run_nestlogit(input_data = d, ...),
    mxl = function(d, ...) run_mxlogit(input_data = d, eta_draws = eta, ...)
  )
  not_finite <- "Weights must be finite, but non-finite values (NA/NaN/Inf) were found."
  not_positive <- "Weights must be strictly positive, but values <= 0 were found."
  for (m in names(preps)) {
    for (bad in c(NA, NaN, Inf, -Inf, 0, -1)) {
      d <- preps[[m]]
      d$weights[3L] <- bad
      # Before the warning about non-uniform weights, which named another
      # cause
      expect_no_warning(expect_error(
        suppressMessages(fits[[m]](d)),
        if (is.finite(bad)) not_positive else not_finite,
        fixed = TRUE, info = paste(m, bad)))
    }
    # These used to fit at the start values, the robust variance's without a
    # warning
    d <- preps[[m]]
    d$weights[c(1L, 25L)] <- c(NA, NaN)
    expect_error(suppressMessages(fits[[m]](d, se_method = "sandwich")),
                 not_finite, fixed = TRUE, info = m)
  }
  d <- preps$mxl
  d$weights[1L] <- NA
  expect_error(
    suppressMessages(fits$mxl(d, concentrate_out_asc = "profile_likelihood")),
    not_finite, fixed = TRUE)
  # A hand-built list whose weights the wrappers reach by `$`'s partial
  # matching is checked as they read it
  d <- preps$mnl
  names(d)[names(d) == "weights"] <- "weights_raw"
  d$weights_raw[1L] <- NA
  expect_error(suppressMessages(fits$mnl(d)), not_finite, fixed = TRUE)
  # A panel: the first decision maker's weights
  dt[, person := (id - 1L) %/% 4L + 1L]
  pan <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2",
                          person_col = "person")
  pan$weights[1:4] <- NA
  expect_error(
    suppressMessages(run_mxlogit(input_data = pan,
                                 eta_draws = get_halton_normals(10L, 30L, 1L))),
    not_finite, fixed = TRUE)
  expect_identical(counter$n, 0L)
})

test_that("input_data's X and W, and eta_draws, must be finite", {
  dt <- cbd_data()
  counter <- cbd_count_evals()
  eta <- get_halton_normals(10L, 120L, 1L)
  mnl <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
  nl <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
  mxl <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  msg <- function(what) paste0("`", what, "` must contain only finite values.")
  for (bad in c(NA, NaN, Inf, -Inf)) {
    d <- mnl
    d$X[5L, 2L] <- bad
    expect_error(suppressMessages(run_mnlogit(input_data = d)),
                 msg("input_data$X"), fixed = TRUE, info = format(bad))
    d <- nl
    d$X[5L, 2L] <- bad
    expect_error(suppressMessages(run_nestlogit(input_data = d)),
                 msg("input_data$X"), fixed = TRUE, info = format(bad))
    d <- mxl
    d$X[5L, 1L] <- bad
    expect_error(suppressMessages(run_mxlogit(input_data = d, eta_draws = eta)),
                 msg("input_data$X"), fixed = TRUE, info = format(bad))
    d <- mxl
    d$W[5L, 1L] <- bad
    expect_error(suppressMessages(run_mxlogit(input_data = d, eta_draws = eta)),
                 msg("input_data$W"), fixed = TRUE, info = format(bad))
    e <- eta
    e[1L, 4L, 7L] <- bad
    expect_error(suppressMessages(run_mxlogit(input_data = mxl, eta_draws = e)),
                 msg("eta_draws"), fixed = TRUE, info = format(bad))
  }
  # A cube stored as integers is checked once converted
  e <- array(1L, dim(eta))
  e[2L] <- NA
  expect_error(suppressMessages(run_mxlogit(input_data = mxl, eta_draws = e)),
               msg("eta_draws"), fixed = TRUE)
  expect_identical(counter$n, 0L)
})

test_that(".all_finite() reads every numeric and logical type", {
  expect_true(.all_finite(NULL))
  # An empty x never reaches min() and max(), which warn on it
  expect_no_warning(expect_true(.all_finite(numeric(0))))
  expect_true(.all_finite(matrix(1:4, 2L)))
  expect_false(.all_finite(c(1L, NA)))
  expect_false(.all_finite(c(TRUE, NA)))
  for (bad in c(NA, NaN, Inf, -Inf)) {
    expect_false(.all_finite(c(bad, 1, 2)), info = format(bad))
    expect_false(.all_finite(c(1, 2, bad)), info = format(bad))
  }
  # A design of another type is left to the kernels' type errors
  d <- list(weights = 1, X = data.frame(x = NA_real_))
  expect_identical(.check_input_data(d), d)
})

test_that("integer64 weights in input_data are read as their values", {
  skip_if_not_installed("bit64")
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # so that the fits agree bit for bit
  dt <- cbd_data()
  eta <- get_halton_normals(10L, 120L, 1L)
  w <- rep(1:3, 40L)
  preps <- list(
    mnl = prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2")),
    nl = prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest"),
    mxl = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  )
  # The robust variance, which weighted fits call for, without its warning
  fits <- list(
    mnl = function(d) run_mnlogit(input_data = d, se_method = "sandwich"),
    nl = function(d) run_nestlogit(input_data = d, se_method = "sandwich"),
    mxl = function(d) run_mxlogit(input_data = d, eta_draws = eta,
                                  se_method = "sandwich")
  )
  for (m in names(preps)) {
    d <- preps[[m]]
    d$weights <- as.double(w)
    f_dbl <- suppressMessages(fits[[m]](d))
    d$weights <- bit64::as.integer64(w)
    f_i64 <- suppressMessages(fits[[m]](d))
    expect_identical(coef(f_i64), coef(f_dbl), info = m)
    expect_identical(f_i64$loglik, f_dbl$loglik, info = m)
    expect_identical(f_i64$vcov, f_dbl$vcov, info = m)
    expect_identical(f_i64$data$weights, as.double(w), info = m)
  }
  # A missing integer64 weight is missing, not its raw bits (those of -0)
  d <- preps$mnl
  d$weights <- bit64::as.integer64(c(NA, w[-1L]))
  expect_error(suppressMessages(run_mnlogit(input_data = d)),
               "Weights must be finite", fixed = TRUE)
  # Weights rounded to doubles warn, naming input_data's
  big <- bit64::as.integer64(c("9007199254740993", "1"))  # 2^53 + 1, 1
  expect_warning(.check_input_data(list(weights = big)),
                 "`input_data$weights` has integer64 values of magnitude 2^53",
                 fixed = TRUE)
  local_mocked_bindings(.bit64_available = function() FALSE)
  d$weights <- bit64::as.integer64(w)
  expect_error(suppressMessages(run_mnlogit(input_data = d)),
               "`input_data$weights` is of class integer64; install the bit64",
               fixed = TRUE)
})
