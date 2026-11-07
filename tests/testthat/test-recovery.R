# Tests for R/recovery.R (recovery_table, monte_carlo, summary.choicer_mc).

test_that("recovery_table.choicer_fit produces expected columns and row count for MNL", {
  sim <- simulate_mnl_data(N = 500, J = 3, seed = 42L)
  fit <- run_mnlogit(
    data                   = sim$data,
    id_col                 = "id",
    alt_col                = "alt",
    choice_col             = "choice",
    covariate_cols         = c("x1", "x2"),
    outside_opt_label      = 0L,
    include_outside_option = FALSE,
    use_asc                = TRUE,
    control                = list(print_level = 0L)
  )

  rt <- recovery_table(fit, sim)

  expect_s3_class(rt, "choicer_recovery")
  expect_setequal(
    names(rt),
    c("parameter", "group", "true", "estimate", "se", "bias",
      "rel_bias_pct", "z_vs_true", "lower_ci", "upper_ci", "covers")
  )
  # beta (2) + asc (3) = 5 rows (outside=0 present in data, treated as inside;
  # first inside alt's ASC is normalized to zero -> 3 free ASCs matches delta len 3).
  expect_equal(nrow(rt), 5L)
  expect_setequal(unique(rt$group), c("beta", "asc"))
  expect_true(all(!is.na(rt$estimate)))
  expect_true(all(rt$upper_ci >= rt$lower_ci))
})

test_that("recovery_table accepts either a choicer_sim or a bare true_params list", {
  sim <- simulate_mnl_data(N = 300, J = 3, seed = 11L)
  fit <- run_mnlogit(
    data                   = sim$data,
    id_col                 = "id",
    alt_col                = "alt",
    choice_col             = "choice",
    covariate_cols         = c("x1", "x2"),
    outside_opt_label      = 0L,
    include_outside_option = FALSE,
    use_asc                = TRUE,
    control                = list(print_level = 0L)
  )

  rt_sim  <- recovery_table(fit, sim)
  rt_list <- recovery_table(fit, sim$true_params)
  expect_equal(rt_sim, rt_list)
})

test_that("recovery_table handles normalized-ASC case (no outside in data/fit)", {
  # Reproducer for the length-mismatch crash: DGP omits the outside option and
  # the fit runs with include_outside_option = FALSE, so the estimator fixes
  # the first inside alt's ASC to zero. Truth has J deltas, fit has J-1 ASCs.
  sim <- simulate_mnl_data(N = 500, J = 4, seed = 42L,
                           outside_option = FALSE, vary_choice_set = FALSE)
  fit <- run_mnlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), use_asc = TRUE,
    control = list(print_level = 0L)
  )
  expect_silent(rt <- recovery_table(fit, sim))
  expect_s3_class(rt, "choicer_recovery")
  # beta (2) + asc (J - 1 = 3) = 5 rows; truth$delta[1] is dropped internally.
  expect_equal(nrow(rt), 5L)
  expect_equal(sum(rt$group == "asc"), 3L)
})

test_that("recovery_table works for MXL with correlated random coefficients", {
  skip_on_cran()
  sim <- simulate_mxl_data(
    N = 300, J = 3, seed = 7L,
    beta = c(0.5),
    Sigma = matrix(c(0.6, 0.2, 0.2, 0.4), nrow = 2),
    rc_correlation = TRUE,
    outside_option = FALSE, vary_choice_set = FALSE
  )
  input <- prepare_mxl_data(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    rc_correlation = TRUE
  )
  eta <- get_halton_normals(50L, input$N, sim$settings$K_w)
  fit <- run_mxlogit(
    input_data = input, eta_draws = eta,
    rc_mean = FALSE, use_asc = TRUE,
    control = list(print_level = 0L, maxeval = 200L)
  )
  rt <- recovery_table(fit, sim)
  expect_s3_class(rt, "choicer_recovery")
  # sigma block exercises the Cholesky (L_params) column-major packing.
  expect_true("sigma" %in% unique(rt$group))
  # Rows per block: beta (1) + sigma (3) + asc (J - 1 = 2) = 6.
  expect_equal(nrow(rt), 6L)
})

test_that("recovery_table works for NL (exercises lambda block)", {
  skip_on_cran()
  sim <- simulate_nl_data(N = 2000, seed = 123L)
  fit <- run_nestlogit(
    data = sim$data, id_col = "id", alt_col = "j", choice_col = "choice",
    covariate_cols = c("X", "W"), nest_col = "nest",
    use_asc = TRUE, include_outside_option = TRUE, outside_opt_label = 0L,
    control = list(print_level = 0L)
  )
  rt <- recovery_table(fit, sim)
  expect_s3_class(rt, "choicer_recovery")
  expect_true("lambda" %in% unique(rt$group))
  # 2 beta + 2 lambda + 5 asc (all inside alts) = 9.
  expect_equal(nrow(rt), 9L)
})

test_that("recovery_table works for MNP (posterior summaries + sigma block)", {
  sim <- simulate_mnp_data(N = 300, J = 3, seed = 42L)
  fit <- suppressMessages(
    run_mnprobit(
      data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), use_asc = TRUE,
      mcmc = list(R = 400, burn = 100, seed = 5)
    )
  )

  rt <- recovery_table(fit, sim)
  expect_s3_class(rt, "choicer_recovery")
  expect_setequal(
    names(rt),
    c("parameter", "group", "true", "estimate", "se", "bias",
      "rel_bias_pct", "z_vs_true", "lower_ci", "upper_ci", "covers")
  )
  # 2 beta + 2 asc + 3 sigma (lower triangle of the 2x2 identified Sigma) = 7.
  expect_equal(nrow(rt), 7L)
  expect_setequal(unique(rt$group), c("beta", "asc", "sigma"))

  # Sigma_11 row is the identification normalization: exact, zero posterior SD
  s11 <- rt[parameter == "Sigma_11"]
  expect_equal(s11$true, 1)
  expect_equal(s11$estimate, 1)
  expect_equal(s11$se, 0)
  expect_true(s11$covers)

  # Accepts a bare true_params list too
  expect_equal(rt, recovery_table(fit, sim$true_params))
})

test_that("monte_carlo MNL smoke test: coverage and bias in plausible ranges", {
  skip_on_cran()

  sim_fun <- function(seed) {
    simulate_mnl_data(N = 1000, J = 3, seed = seed)
  }
  fit_fun <- function(sim) {
    run_mnlogit(
      data                   = sim$data,
      id_col                 = "id",
      alt_col                = "alt",
      choice_col             = "choice",
      covariate_cols         = c("x1", "x2"),
      outside_opt_label      = 0L,
      include_outside_option = FALSE,
      use_asc                = TRUE,
      control                = list(print_level = 0L)
    )
  }

  mc <- monte_carlo(sim_fun, fit_fun, R = 20L, seed = 1L, progress = FALSE)
  expect_s3_class(mc, "choicer_mc")
  expect_equal(mc$meta$R, 20L)

  s <- summary(mc)
  expect_s3_class(s, "choicer_mc_summary")
  expect_gte(s$conv_rate, 0.9)

  tbl <- s$summary
  # One summary row per parameter; MNL with J=3 and ASCs has 2 betas + 3 ASCs.
  expect_gte(nrow(tbl), 5L)
  # Coverage should be within a generous band at R = 20.
  expect_true(all(tbl$coverage >= 0.6 & tbl$coverage <= 1.0))
  # Mean bias within 3 MC-standard-errors of zero.
  expect_true(all(abs(tbl$bias) <= 3 * tbl$sd_est / sqrt(tbl$R_success) + 1e-3))
})

# --- summary.choicer_mc(level =) ----------------------------------------------

test_that("summary.choicer_mc() computes coverage at `level`", {
  # Fake reps with known distances from the truth, in SEs: |z| = 0.1, 0.5, 1,
  # 1.5, 1.7, 1.9, 2.1, 2.3, 3 and an NA SE
  z <- c(0.1, -0.5, 1, -1.5, 1.7, -1.9, 2.1, -2.3, 3, 0)
  se <- rep(0.2, 10)
  se[10] <- NA_real_
  est <- 1 + z * 0.2
  reps <- data.table(
    rep_id = 1:10, seed = 1:10, parameter = "b", group = "beta", true = 1,
    estimate = est, se = se, bias = est - 1, rel_bias_pct = NA_real_,
    z_vs_true = z, lower_ci = NA_real_, upper_ci = NA_real_,
    covers = (1 >= est - qnorm(0.975) * se) & (1 <= est + qnorm(0.975) * se),
    loglik = -1, converged = TRUE, time_sec = 0, error = NA_character_)
  covers0 <- data.table::copy(reps$covers)
  mc <- structure(list(replications = reps,
                       meta = list(R = 10L, seed = 1L, seeds = 1:10, parallel = FALSE)),
                  class = "choicer_mc")
  cov_at <- function(level) summary(mc, level = level)$summary$coverage
  # the NA SE counts in neither the numerator nor the denominator
  expect_equal(cov_at(0.95), 6 / 9)      # |z| <= 1.96
  expect_equal(cov_at(0.50), 2 / 9)      # |z| <= 0.674
  expect_equal(cov_at(0.80), 3 / 9)      # |z| <= 1.282
  expect_equal(cov_at(0.99), 8 / 9)      # |z| <= 2.576
  # the default reads as the reps' own indicator, and recovery_table() agrees
  expect_equal(summary(mc)$summary$coverage, mean(reps$covers, na.rm = TRUE))
  expect_equal(recovery_table(mc, level = 0.5)$summary$coverage, 2 / 9)
  expect_identical(summary(mc, level = 0.5)$level, 0.5)
  # summary() leaves the object's reps as they were
  expect_identical(mc$replications$covers, covers0)
})

# --- monte_carlo() on the samplers' fits --------------------------------------

test_that("monte_carlo() runs on MNP and HB fits, which have no likelihood", {
  skip_on_cran()
  mnp <- suppressWarnings(monte_carlo(
    function(seed) simulate_mnp_data(N = 150, J = 3, seed = seed),
    function(sim) {
      run_mnprobit(sim$data, "id", "alt", "choice", c("x1", "x2"),
                   mcmc = list(R = 300, burn = 100))
    },
    R = 2L, seed = 5L, progress = FALSE))
  expect_s3_class(mnp, "choicer_mc")
  # a rep that fails is recorded with its error: none did
  expect_true(all(is.na(mnp$replications$error)))
  expect_true(all(is.na(mnp$replications$loglik)))
  expect_true(all(mnp$replications$converged))
  s <- summary(mnp)
  expect_equal(s$conv_rate, 1)
  expect_gt(nrow(s$summary), 0L)
  expect_true(all(s$summary$R_success == 2L))

  hb <- suppressWarnings(monte_carlo(
    function(seed) simulate_hmnl_data(N = 40, T = 2, J = 3, seed = seed),
    function(sim) {
      run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"), person_col = "pid",
                   alt_covariate_cols = "z1", mcmc = list(R = 200, burn = 50))
    },
    R = 2L, seed = 7L, progress = FALSE))
  expect_true(all(is.na(hb$replications$error)))
  expect_true(all(is.na(hb$replications$loglik)))
  # One row per alternative effect across the reps, whose truth differs by seed
  s_hb <- summary(hb)$summary
  delta <- s_hb[group == "delta"]
  expect_identical(anyDuplicated(delta$parameter), 0L)
  expect_true(all(delta$R_success == 2L))
  expect_true(all(is.finite(delta$sd_est)))
})
