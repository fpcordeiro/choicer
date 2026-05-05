# Tests for mc_asymptotics() and print.choicer_mc_asymptotics().
#
# These tests build a choicer_mc fixture by hand (no real fitting): given
# known sampling distributions for each parameter across "replications",
# we assert that mc_asymptotics() recovers the correct bias, coverage,
# moments, SE ratio, and pass/fail flags. Testing the math, not the fit
# pipeline.

make_fake_choicer_mc <- function(parameter, group, true, estimates, ses,
                                 converged = NULL,
                                 level = 0.95) {
  R <- length(estimates)
  if (is.null(converged)) converged <- rep(TRUE, R)
  stopifnot(length(ses) == R, length(converged) == R)

  z_crit <- stats::qnorm(1 - (1 - level) / 2)
  lower_ci <- estimates - z_crit * ses
  upper_ci <- estimates + z_crit * ses
  covers <- (true >= lower_ci) & (true <= upper_ci)

  replications <- data.table::data.table(
    rep_id       = seq_len(R),
    seed         = seq_len(R),
    parameter    = parameter,
    group        = group,
    true         = true,
    estimate     = estimates,
    se           = ses,
    bias         = estimates - true,
    rel_bias_pct = 100 * (estimates - true) / true,
    z_vs_true    = (estimates - true) / ses,
    lower_ci     = lower_ci,
    upper_ci     = upper_ci,
    covers       = covers,
    loglik       = NA_real_,
    converged    = converged,
    time_sec     = NA_real_,
    error        = NA_character_
  )

  structure(
    list(
      replications = replications,
      meta = list(R = R, seed = 1L, seeds = seq_len(R),
                  parallel = FALSE, timestamp = Sys.time())
    ),
    class = "choicer_mc"
  )
}

# Build a fake MC for a single parameter with iid draws (est, se).
single_param_mc <- function(estimates, ses, true = 0, name = "beta_1",
                            group = "beta", converged = NULL) {
  R <- length(estimates)
  make_fake_choicer_mc(
    parameter = rep(name, R), group = rep(group, R),
    true = rep(true, R), estimates = estimates, ses = ses,
    converged = converged
  )
}

test_that("mc_asymptotics recovers bias, sd_emp, mean_se on canned normal draws", {
  set.seed(2026)
  R <- 2000
  true <- 1.0
  sd_true <- 0.1
  estimates <- stats::rnorm(R, mean = true, sd = sd_true)
  ses <- rep(sd_true, R)  # reported SE equals true SD

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  expect_s3_class(out, "choicer_mc_asymptotics")
  expect_equal(nrow(out), 1L)
  expect_equal(out$R_used, R)
  expect_equal(out$R_excluded, 0L)
  # bias within 3 MC-SE of zero
  expect_lt(abs(out$bias_mc_se), 3)
  # se_ratio ~ 1
  expect_true(abs(out$se_ratio - 1) < 0.1)
  # Passes all flags
  expect_true(out$pass_bias)
  expect_true(out$pass_se_ratio)
  expect_true(out$pass_cov95)
  expect_true(out$pass_skew)
  expect_true(out$pass_kurt)
})

test_that("mc_asymptotics detects large bias via pass_bias = FALSE", {
  set.seed(2027)
  R <- 500
  true <- 0
  # Shift mean well away from truth.
  estimates <- stats::rnorm(R, mean = 0.5, sd = 0.1)
  ses <- rep(0.1, R)

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  expect_false(out$pass_bias)
  expect_gt(abs(out$bias_mc_se), 3)
})

test_that("mc_asymptotics detects SE/SD mismatch via pass_se_ratio = FALSE", {
  set.seed(2028)
  R <- 500
  true <- 0
  sd_true <- 0.2
  estimates <- stats::rnorm(R, mean = true, sd = sd_true)
  # Report SEs that are half the true SD: se_ratio ~ 0.5.
  ses <- rep(sd_true / 2, R)

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  expect_false(out$pass_se_ratio)
  expect_lt(out$se_ratio, 0.9)
})

test_that("mc_asymptotics empirical coverage matches the reported-SE band", {
  set.seed(2029)
  R <- 2000
  true <- 0
  sd_true <- 0.15
  # Estimates drawn from N(true, sd_true^2); SEs match true SD exactly.
  estimates <- stats::rnorm(R, mean = true, sd = sd_true)
  ses <- rep(sd_true, R)

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  # Nominal 95 pct coverage should land near 0.95 with R=2000.
  expect_lt(abs(out$cov95 - 0.95), 0.03)
  expect_true(out$pass_cov95)
  # Wilson band contains the nominal level.
  expect_lte(out$cov95_lower, 0.95)
  expect_gte(out$cov95_upper, 0.95)

  # Lower nominal levels should give lower empirical coverage.
  expect_lt(out$cov90, out$cov95)
  expect_lt(out$cov95, out$cov99)
})

test_that("mc_asymptotics normality tests pass on normal z and flag skew/kurt", {
  set.seed(2030)
  R <- 1000
  true <- 0
  sd_true <- 1
  estimates <- stats::rnorm(R, mean = true, sd = sd_true)
  ses <- rep(sd_true, R)

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  expect_true(abs(out$mean_z) < 0.15)
  expect_true(abs(out$sd_z - 1) < 0.1)
  expect_true(abs(out$skew_z) < 0.3)
  expect_true(abs(out$kurt_excess_z) < 0.5)
  # p-values should generally be above 0.01 under the null.
  expect_gt(out$jb_p, 0.01)

  # Heavy-tailed input (t3): normality should fail, skew OK but kurtosis fails.
  set.seed(31)
  R2 <- 500
  heavy <- stats::rt(R2, df = 3) * sd_true + true
  ses2 <- rep(sd_true, R2)
  mc2 <- single_param_mc(heavy, ses2, true = true)
  out2 <- mc_asymptotics(mc2)
  # Excess kurtosis of t3 is undefined / very large in samples; pass_kurt
  # should fail.
  expect_false(out2$pass_kurt)
})

test_that("mc_asymptotics excludes non-converged reps and reports R_excluded", {
  set.seed(2032)
  R <- 100
  estimates <- stats::rnorm(R, mean = 0, sd = 0.1)
  ses <- rep(0.1, R)
  converged <- c(rep(TRUE, 80), rep(FALSE, 20))

  mc <- single_param_mc(estimates, ses, true = 0, converged = converged)
  out <- mc_asymptotics(mc)

  expect_equal(out$R_used, 80L)
  expect_equal(out$R_excluded, 20L)
})

test_that("mc_asymptotics handles multiple parameters independently", {
  set.seed(2033)
  R <- 300
  estimates <- c(
    stats::rnorm(R, mean = 1.0, sd = 0.1),
    stats::rnorm(R, mean = -0.5, sd = 0.15)
  )
  ses <- c(rep(0.1, R), rep(0.15, R))

  mc <- make_fake_choicer_mc(
    parameter = c(rep("beta_1", R), rep("beta_2", R)),
    group     = c(rep("beta", R), rep("beta", R)),
    true      = c(rep(1.0, R), rep(-0.5, R)),
    estimates = estimates,
    ses       = ses
  )
  # Replace rep_id so each parameter has aligned rep ids (both go 1..R).
  mc$replications[, rep_id := c(seq_len(R), seq_len(R))]
  out <- mc_asymptotics(mc)

  expect_equal(nrow(out), 2L)
  expect_setequal(out$parameter, c("beta_1", "beta_2"))
  # Both should pass bias and se_ratio flags.
  expect_true(all(out$pass_bias))
  expect_true(all(out$pass_se_ratio))
})

test_that("mc_asymptotics winsorized columns differ from raw under heavy tails", {
  set.seed(2034)
  R <- 500
  true <- 0
  # 95% clean + 5% outliers in one tail.
  inner <- stats::rnorm(floor(R * 0.95), mean = 0, sd = 0.1)
  outer <- stats::rnorm(R - length(inner), mean = 10, sd = 0.1)
  estimates <- c(inner, outer)
  ses <- rep(0.1, R)

  mc <- single_param_mc(estimates, ses, true = true)
  out <- mc_asymptotics(mc)

  # Winsorized bias should be smaller in magnitude.
  expect_lt(abs(out$bias_w), abs(out$bias))
  # Winsorized SD smaller too.
  expect_lt(out$sd_emp_w, out$sd_emp)
})

test_that("Wilson CI contains observed proportion and stays in [0, 1]", {
  # Hit the internal helper via :::; it is unexported.
  wilson <- choicer:::.wilson_ci
  ci <- wilson(0.5, 100, level = 0.95)
  expect_gte(ci$lower, 0); expect_lte(ci$upper, 1)
  expect_lt(ci$lower, 0.5); expect_gt(ci$upper, 0.5)
  # Boundary: p = 0.
  ci0 <- wilson(0, 100, level = 0.95)
  expect_equal(ci0$lower, 0)
  # n = 0: both NA.
  ci_empty <- wilson(0.3, 0, level = 0.95)
  expect_true(is.na(ci_empty$lower) && is.na(ci_empty$upper))
})

test_that("Jarque-Bera p-value is high for normal data, low for skewed data", {
  jb_p <- choicer:::.jarque_bera_p
  set.seed(2035)
  expect_gt(jb_p(stats::rnorm(1000)), 0.05)
  expect_lt(jb_p(stats::rexp(1000, rate = 1)), 0.01)
})

test_that("mc_asymptotics(se_col) selects alternative SE column", {
  # Fixture: est ~ N(true, sd_true). Two SE columns on the replications:
  # `se` is correctly calibrated (matches sd_true); `se_bhhh` is inflated
  # by 2x. We should see se_ratio ~ 1 with default se_col = "se" and
  # se_ratio ~ 2 with se_col = "se_bhhh".
  set.seed(2050)
  R <- 1000
  true <- 0.5
  sd_true <- 0.1
  estimates <- stats::rnorm(R, mean = true, sd = sd_true)
  ses <- rep(sd_true, R)
  mc <- single_param_mc(estimates, ses, true = true)
  mc$replications[, se_bhhh := 2 * ses]

  out_hess <- mc_asymptotics(mc, se_col = "se")
  out_bhhh <- mc_asymptotics(mc, se_col = "se_bhhh")

  expect_true(abs(out_hess$se_ratio - 1) < 0.1)
  expect_true(abs(out_bhhh$se_ratio - 2) < 0.1)

  # meta attribute records which column was used.
  expect_equal(attr(out_hess, "meta")$se_col, "se")
  expect_equal(attr(out_bhhh, "meta")$se_col, "se_bhhh")

  # cov95 with inflated SE should be > 0.95 (bands are too wide).
  expect_gt(out_bhhh$cov95, out_hess$cov95)

  # Error when se_col doesn't exist.
  expect_error(mc_asymptotics(mc, se_col = "does_not_exist"),
               "not found in mc\\$replications")
})

test_that("print.choicer_mc_asymptotics shows pass/fail matrix header", {
  set.seed(2036)
  R <- 300
  estimates <- stats::rnorm(R, mean = 0, sd = 0.1)
  ses <- rep(0.1, R)
  mc <- single_param_mc(estimates, ses, true = 0)
  out <- mc_asymptotics(mc)

  out_lines <- capture.output(print(out))
  expect_true(any(grepl("choicer_mc_asymptotics", out_lines)))
  expect_true(any(grepl("pass_bias", out_lines)))
  expect_true(any(grepl("pass_se_ratio", out_lines)))
  expect_true(any(grepl("pass_cov95", out_lines)))
  expect_true(any(grepl("all_pass", out_lines)))
})
