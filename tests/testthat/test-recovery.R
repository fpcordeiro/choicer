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
