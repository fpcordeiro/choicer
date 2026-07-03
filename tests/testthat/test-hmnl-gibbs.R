# Tests for the HMNL Gibbs kernel (hmnl_gibbs) and the run_hmnlogit wrapper.
# CRAN-fast: small fixtures, short chains. Recovery accuracy lives in
# test-hmnl-recovery.R (skip_on_cran).

make_hmnl_fixture <- function(N = 20, T = 2, J = 3, seed = 42, ...) {
  sim <- simulate_hmnl_data(N = N, T = T, J = J, seed = seed,
                            theta = c(0.5, -0.4), sigma_d = 0.5, ...)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  list(sim = sim, d = d)
}

make_hmnl_args <- function(d, R = 200, burn = 50, thin = 1, seed = 7,
                           keep_beta_i = 1L) {
  list(
    X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
    include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
    rc_dist = as.integer(d$rc_dist), beta_pooled = rep(0, d$K_struct),
    delta_init = rep(0, d$J), theta_init = rep(0, d$P),
    b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
    nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
    theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
    sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
    R = R, burn = burn, thin = thin, seed = seed,
    keep_beta_i = keep_beta_i, s_init = 2.38 / sqrt(d$K_struct),
    accept_target = 0.234
  )
}

test_that("hmnl_gibbs returns correctly shaped, finite draws", {
  fx <- make_hmnl_fixture()
  d <- fx$d
  out <- do.call(hmnl_gibbs, make_hmnl_args(d, R = 200, burn = 50, thin = 3))

  expect_equal(out$R_keep, 50L)            # ceil(150 / 3)
  expect_equal(dim(out$bdraw), c(50L, d$K_struct))
  expect_equal(dim(out$wdraw), c(50L, d$K_struct * (d$K_struct + 1L) / 2L))
  expect_equal(dim(out$deltadraw), c(50L, d$J))
  expect_equal(dim(out$thetadraw), c(50L, d$P))
  expect_length(as.numeric(out$sigma_d2draw), 50L)
  expect_length(as.numeric(out$loglik_trace), 50L)

  expect_true(all(is.finite(out$bdraw)))
  expect_true(all(is.finite(out$wdraw)))
  expect_true(all(is.finite(out$deltadraw)))
  expect_true(all(is.finite(out$thetadraw)))
  expect_true(all(as.numeric(out$sigma_d2draw) > 0))
  expect_true(all(is.finite(out$loglik_trace)))

  # W draws: every kept diagonal is positive (vech row-major lower triangle;
  # note %/% binds tighter than *, hence the explicit parentheses)
  diag_idx <- vapply(seq_len(d$K_struct),
                     function(k) (k * (k + 1L)) %/% 2L, integer(1L))
  expect_true(all(out$wdraw[, diag_idx] > 0))

  # acceptance rates and adapted scales are valid
  expect_true(all(out$accept_rate_beta >= 0 & out$accept_rate_beta <= 1))
  expect_true(all(out$accept_rate_delta >= 0 & out$accept_rate_delta <= 1))
  expect_true(all(out$s_final >= 0.01 & out$s_final <= 10))
  expect_true(all(out$s_delta_final >= 0.01 & out$s_delta_final <= 10))

  # beta_i summaries (keep_beta_i = 1)
  expect_equal(dim(out$beta_i_mean), c(d$K_struct, d$N_persons))
  expect_equal(dim(out$beta_i_sd), c(d$K_struct, d$N_persons))
  expect_null(out$beta_i_draws)
  expect_true(all(is.finite(out$beta_i_mean)))

  # quality-ladder summaries always present
  expect_length(as.numeric(out$delta_mean), d$J)
  expect_length(as.numeric(out$xi_mean), d$J)
})

test_that("hmnl_gibbs draws are seed-reproducible and thread-invariant", {
  fx <- make_hmnl_fixture()
  args <- make_hmnl_args(fx$d)

  set_num_threads(1)
  out1 <- do.call(hmnl_gibbs, args)
  set_num_threads(4)
  out4 <- do.call(hmnl_gibbs, args)

  # Bitwise identical across thread counts: per-(iteration, unit) RNG
  # streams + fixed-order accumulation + the serial delta sweep.
  expect_identical(out1$bdraw, out4$bdraw)
  expect_identical(out1$wdraw, out4$wdraw)
  expect_identical(out1$deltadraw, out4$deltadraw)
  expect_identical(out1$thetadraw, out4$thetadraw)
  expect_identical(out1$sigma_d2draw, out4$sigma_d2draw)
  expect_identical(out1$loglik_trace, out4$loglik_trace)

  # Same seed twice -> identical; different seed -> different
  out_again <- do.call(hmnl_gibbs, args)
  expect_identical(out4$bdraw, out_again$bdraw)
  args2 <- args
  args2$seed <- 999
  expect_false(identical(do.call(hmnl_gibbs, args2)$bdraw, out1$bdraw))
})

test_that("hmnl_gibbs is thread-invariant in the J > N regime", {
  # 10 respondents, 25 alternatives: the master RNG tags (2N + J + 0..3)
  # sit beyond every delta tag only because the base is 2N + J, not 3N.
  sim <- simulate_hmnl_data(N = 10, T = 1, J = 25, seed = 5,
                            theta = c(0.2, -0.3), sigma_d = 0.4)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         alt_covariate_cols = "z1")
  expect_true(d$J > d$N_persons)
  args <- make_hmnl_args(d, R = 100, burn = 30)

  set_num_threads(1)
  out1 <- do.call(hmnl_gibbs, args)
  set_num_threads(4)
  out4 <- do.call(hmnl_gibbs, args)
  expect_identical(out1$bdraw, out4$bdraw)
  expect_identical(out1$deltadraw, out4$deltadraw)
  expect_identical(out1$sigma_d2draw, out4$sigma_d2draw)
})

test_that("hmnl_gibbs keep_beta_i modes are consistent", {
  fx <- make_hmnl_fixture()
  args0 <- make_hmnl_args(fx$d, keep_beta_i = 0L)
  args2 <- make_hmnl_args(fx$d, keep_beta_i = 2L)

  out0 <- do.call(hmnl_gibbs, args0)
  out2 <- do.call(hmnl_gibbs, args2)

  expect_null(out0$beta_i_mean)
  expect_null(out0$beta_i_draws)
  # identical population draws regardless of storage mode
  expect_identical(out0$bdraw, out2$bdraw)

  # cube means match the Welford means
  expect_equal(dim(out2$beta_i_draws),
               c(fx$d$K_struct, fx$d$N_persons, out2$R_keep))
  expect_equal(apply(out2$beta_i_draws, c(1, 2), mean),
               unname(as.matrix(out2$beta_i_mean)), tolerance = 1e-12)
})

test_that("hmnl_gibbs validates inputs", {
  fx <- make_hmnl_fixture()
  args <- make_hmnl_args(fx$d)

  a <- args; a$include_outside_option <- FALSE
  expect_error(do.call(hmnl_gibbs, a), "include_outside_option")
  a <- args; a$R <- 0L
  expect_error(do.call(hmnl_gibbs, a), "R must be")
  a <- args; a$burn <- args$R
  expect_error(do.call(hmnl_gibbs, a), "burn")
  a <- args; a$b_bar <- rep(0, 5)
  expect_error(do.call(hmnl_gibbs, a), "b_bar")
  a <- args; a$rc_dist <- c(0L, 2L)
  expect_error(do.call(hmnl_gibbs, a), "rc_dist")
  a <- args; a$delta_init <- rep(0, 2)
  expect_error(do.call(hmnl_gibbs, a), "delta_init")
  a <- args; a$seed <- -1
  expect_error(do.call(hmnl_gibbs, a), "seed")
  a <- args; a$sd_prior <- list(half_cauchy = TRUE, s_d = -1, c0 = 3, d0 = 3)
  expect_error(do.call(hmnl_gibbs, a), "sd_prior")
})

test_that("run_hmnlogit is reproducible end-to-end via set.seed", {
  fx <- make_hmnl_fixture()

  set.seed(123)
  fit1 <- suppressWarnings(run_hmnlogit(input_data = fx$d,
                                        mcmc = list(R = 150, burn = 50)))
  set.seed(123)
  fit2 <- suppressWarnings(run_hmnlogit(input_data = fx$d,
                                        mcmc = list(R = 150, burn = 50)))
  # The Gibbs kernel is bitwise reproducible (tested above), but the pooled
  # MLE initialization runs through the OpenMP frequentist likelihood, whose
  # reduction order is only machine-precision stable across calls — so the
  # end-to-end check is tolerance-based, not identical().
  expect_equal(fit1$draws$b, fit2$draws$b, tolerance = 1e-8)
  expect_equal(fit1$draws$delta, fit2$draws$delta, tolerance = 1e-8)

  set.seed(456)
  fit3 <- suppressWarnings(run_hmnlogit(input_data = fx$d,
                                        mcmc = list(R = 150, burn = 50)))
  expect_false(identical(fit1$draws$b, fit3$draws$b))
})

test_that("run_hmnlogit builds a complete classed fit", {
  fx <- make_hmnl_fixture()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 200, burn = 50),
                 chains = 2)
  )

  expect_s3_class(fit, c("choicer_hmnl", "choicer_hb"), exact = TRUE)
  expect_equal(fit$model, "hmnl")
  expect_named(fit$coefficients, colnames(fx$d$X))
  expect_equal(dim(fit$vcov), c(fx$d$K_struct, fx$d$K_struct))
  expect_equal(nobs(fit), fx$d$n_tasks)
  expect_equal(fit$n_persons, fx$d$N_persons)
  expect_equal(nrow(fit$delta), fx$d$J)
  expect_equal(nrow(fit$xi), fx$d$J)
  expect_true(all(is.finite(fit$draws$sigma_d2)) &&
                all(fit$draws$sigma_d2 > 0))
  expect_true(is.numeric(fit$rhat) && length(fit$rhat) >= fx$d$K_struct)

  # coef components
  expect_length(coef(fit), fx$d$K_struct)
  expect_length(coef(fit, "theta"), fx$d$P)
  expect_length(coef(fit, "delta"), fx$d$J)
  expect_length(coef(fit, "xi"), fx$d$J)

  # print/summary smoke
  expect_output(print(fit), "Hierarchical Bayesian Multinomial Logit")
  s <- summary(fit)
  expect_s3_class(s, "summary.choicer_hb")
  expect_output(print(s), "Quality ladder")
  expect_output(print(s), "split-R-hat")
})

test_that("run_hmnlogit input validation and guards", {
  fx <- make_hmnl_fixture()

  expect_error(run_hmnlogit(), "Supply either")
  expect_error(run_hmnlogit(data = fx$sim$data, input_data = fx$d),
               "not both")
  expect_error(run_hmnlogit(input_data = list(a = 1)), "choicer_data_hmnl")
  expect_error(run_hmnlogit(fx$sim$data, "task", "alt", "choice"),
               "Convenience workflow requires")

  # v1 requires the outside option (identification anchor)
  sim_no <- simulate_hmnl_data(N = 15, T = 2, J = 3, seed = 9,
                               include_outside = FALSE)
  d_no_out <- prepare_hmnl_data(sim_no$data, "task", "alt", "choice",
                                c("x1", "x2"), person_col = "pid",
                                include_outside_option = FALSE)
  expect_error(run_hmnlogit(input_data = d_no_out),
               "include_outside_option = TRUE")

  # prior dimension checks
  expect_error(
    run_hmnlogit(input_data = fx$d, prior = list(b_bar = rep(0, 5)),
                 mcmc = list(R = 100, burn = 20)),
    "b_bar"
  )
  expect_error(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 100, burn = 100)),
    "greater than"
  )
})
