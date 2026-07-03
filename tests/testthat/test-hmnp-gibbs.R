# Tests for the HMNP Gibbs kernel (hmnp_gibbs) and the run_hmnprobit wrapper.
# CRAN-fast: small fixtures, short chains. Recovery accuracy lives in
# test-hmnp-recovery.R (skip_on_cran).

make_hmnp_fixture <- function(N = 20, T = 2, J = 3, seed = 42, ...) {
  sim <- simulate_hmnp_data(N = N, T = T, J = J, seed = seed,
                            theta = c(0.5, -0.4), sigma_d = 0.5, ...)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  list(sim = sim, d = d)
}

make_hmnp_args <- function(d, R = 200, burn = 50, thin = 1, seed = 7,
                           keep_beta_i = 1L) {
  list(
    X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
    include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
    delta_init = rep(0, d$J), theta_init = rep(0, d$P),
    b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
    nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
    theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
    sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
    a0 = 3, s0 = 3,
    R = R, burn = burn, thin = thin, seed = seed, keep_beta_i = keep_beta_i
  )
}

test_that("hmnp_gibbs returns correctly shaped, finite draws", {
  fx <- make_hmnp_fixture()
  d <- fx$d
  out <- do.call(hmnp_gibbs, make_hmnp_args(d, R = 200, burn = 50, thin = 3))

  expect_equal(out$R_keep, 50L)            # ceil(150 / 3)
  expect_equal(dim(out$bdraw), c(50L, d$K_struct))
  expect_equal(dim(out$wdraw), c(50L, d$K_struct * (d$K_struct + 1L) / 2L))
  expect_equal(dim(out$deltadraw), c(50L, d$J))
  expect_equal(dim(out$thetadraw), c(50L, d$P))
  expect_length(as.numeric(out$sigma_d2draw), 50L)
  expect_length(as.numeric(out$sigma2draw), 50L)

  expect_true(all(is.finite(out$bdraw)))
  expect_true(all(is.finite(out$wdraw)))
  expect_true(all(is.finite(out$deltadraw)))
  expect_true(all(is.finite(out$thetadraw)))
  expect_true(all(as.numeric(out$sigma_d2draw) > 0))
  expect_true(all(as.numeric(out$sigma2draw) > 0))

  # W diagonals positive (vech row-major lower triangle; parentheses matter:
  # %/% binds tighter than *)
  diag_idx <- vapply(seq_len(d$K_struct),
                     function(k) (k * (k + 1L)) %/% 2L, integer(1L))
  expect_true(all(out$wdraw[, diag_idx] > 0))

  # beta_i summaries (identified scale)
  expect_equal(dim(out$beta_i_mean), c(d$K_struct, d$N_persons))
  expect_equal(dim(out$beta_i_sd), c(d$K_struct, d$N_persons))
  expect_null(out$beta_i_draws)
  expect_true(all(is.finite(out$beta_i_mean)))
  expect_length(as.numeric(out$delta_mean), d$J)
  expect_length(as.numeric(out$xi_mean), d$J)
})

test_that("hmnp_gibbs draws are seed-reproducible and thread-invariant", {
  fx <- make_hmnp_fixture()
  args <- make_hmnp_args(fx$d)

  set_num_threads(1)
  out1 <- do.call(hmnp_gibbs, args)
  set_num_threads(4)
  out4 <- do.call(hmnp_gibbs, args)

  # Bitwise identical across thread counts: per-(iteration, unit) RNG
  # streams, fixed-order per-j sufficient statistics, fixed-order RSS
  # reduction on master.
  expect_identical(out1$bdraw, out4$bdraw)
  expect_identical(out1$wdraw, out4$wdraw)
  expect_identical(out1$deltadraw, out4$deltadraw)
  expect_identical(out1$thetadraw, out4$thetadraw)
  expect_identical(out1$sigma_d2draw, out4$sigma_d2draw)
  expect_identical(out1$sigma2draw, out4$sigma2draw)

  # Same seed twice -> identical; different seed -> different
  out_again <- do.call(hmnp_gibbs, args)
  expect_identical(out4$bdraw, out_again$bdraw)
  args2 <- args
  args2$seed <- 999
  expect_false(identical(do.call(hmnp_gibbs, args2)$bdraw, out1$bdraw))
})

test_that("hmnp_gibbs is thread-invariant in the J > N regime", {
  sim <- simulate_hmnp_data(N = 10, T = 1, J = 25, seed = 5,
                            theta = c(0.2, -0.3), sigma_d = 0.4)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         alt_covariate_cols = "z1")
  expect_true(d$J > d$N_persons)
  args <- make_hmnp_args(d, R = 100, burn = 30)

  set_num_threads(1)
  out1 <- do.call(hmnp_gibbs, args)
  set_num_threads(4)
  out4 <- do.call(hmnp_gibbs, args)
  expect_identical(out1$bdraw, out4$bdraw)
  expect_identical(out1$deltadraw, out4$deltadraw)
  expect_identical(out1$sigma2draw, out4$sigma2draw)
})

test_that("hmnp_gibbs keep_beta_i modes are consistent", {
  fx <- make_hmnp_fixture()
  out0 <- do.call(hmnp_gibbs, make_hmnp_args(fx$d, keep_beta_i = 0L))
  out2 <- do.call(hmnp_gibbs, make_hmnp_args(fx$d, keep_beta_i = 2L))

  expect_null(out0$beta_i_mean)
  expect_null(out0$beta_i_draws)
  expect_identical(out0$bdraw, out2$bdraw)

  # the cube holds per-draw-normalized beta_i / sigma draws, so its means
  # match the Welford means on the identified scale
  expect_equal(dim(out2$beta_i_draws),
               c(fx$d$K_struct, fx$d$N_persons, out2$R_keep))
  expect_equal(apply(out2$beta_i_draws, c(1, 2), mean),
               unname(as.matrix(out2$beta_i_mean)), tolerance = 1e-12)
})

test_that("hmnp_gibbs validates inputs", {
  fx <- make_hmnp_fixture()
  args <- make_hmnp_args(fx$d)

  a <- args; a$include_outside_option <- FALSE
  expect_error(do.call(hmnp_gibbs, a), "include_outside_option")
  a <- args; a$R <- 0L
  expect_error(do.call(hmnp_gibbs, a), "R must be")
  a <- args; a$burn <- args$R
  expect_error(do.call(hmnp_gibbs, a), "burn")
  a <- args; a$b_bar <- rep(0, 5)
  expect_error(do.call(hmnp_gibbs, a), "b_bar")
  a <- args; a$delta_init <- rep(0, 2)
  expect_error(do.call(hmnp_gibbs, a), "delta_init")
  a <- args; a$a0 <- -1
  expect_error(do.call(hmnp_gibbs, a), "a0 and s0")
  a <- args; a$seed <- -1
  expect_error(do.call(hmnp_gibbs, a), "seed")
})

test_that("run_hmnprobit is exactly reproducible end-to-end via set.seed", {
  fx <- make_hmnp_fixture()

  # No optimizer init: the HMNP wrapper is deterministic given the kernel
  # seed, so end-to-end reproducibility is exact (contrast run_hmnlogit).
  set.seed(123)
  fit1 <- suppressWarnings(run_hmnprobit(input_data = fx$d,
                                         mcmc = list(R = 150, burn = 50)))
  set.seed(123)
  fit2 <- suppressWarnings(run_hmnprobit(input_data = fx$d,
                                         mcmc = list(R = 150, burn = 50)))
  expect_identical(fit1$draws$b, fit2$draws$b)
  expect_identical(fit1$draws$delta, fit2$draws$delta)

  set.seed(456)
  fit3 <- suppressWarnings(run_hmnprobit(input_data = fx$d,
                                         mcmc = list(R = 150, burn = 50)))
  expect_false(identical(fit1$draws$b, fit3$draws$b))
})

test_that("run_hmnprobit normalizes every identified draw by sigma", {
  fx <- make_hmnp_fixture()
  fit <- suppressWarnings(run_hmnprobit(input_data = fx$d,
                                        mcmc = list(R = 200, burn = 50)))

  s <- sqrt(fit$draws$sigma2)
  expect_equal(fit$draws$b, fit$draws$b_raw / s, ignore_attr = TRUE)
  expect_equal(fit$draws$delta, fit$draws$delta_raw / s,
               ignore_attr = TRUE)
  expect_equal(fit$draws$theta, fit$draws$theta_raw / s,
               ignore_attr = TRUE)
  expect_equal(fit$draws$w_vech, fit$draws$w_vech_raw / fit$draws$sigma2,
               ignore_attr = TRUE)
  expect_equal(fit$draws$sigma_d2, fit$draws$sigma_d2_raw / fit$draws$sigma2,
               ignore_attr = TRUE)
})

test_that("run_hmnprobit builds a complete classed fit", {
  fx <- make_hmnp_fixture()
  fit <- suppressWarnings(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = 200, burn = 50),
                  chains = 2)
  )

  expect_s3_class(fit, c("choicer_hmnp", "choicer_hb"), exact = TRUE)
  expect_equal(fit$model, "hmnp")
  expect_named(fit$coefficients, colnames(fx$d$X))
  expect_equal(nobs(fit), fx$d$n_tasks)
  expect_equal(nrow(fit$delta), fx$d$J)
  expect_true(all(fit$draws$sigma2 > 0))
  expect_true(is.numeric(fit$rhat))

  expect_length(coef(fit), fx$d$K_struct)
  expect_length(coef(fit, "delta"), fx$d$J)

  expect_output(print(fit), "Hierarchical Bayesian Multinomial Probit")
  s <- summary(fit)
  expect_s3_class(s, "summary.choicer_hb")
  expect_output(print(s), "Quality ladder")
  expect_output(print(s), "Raw shock variance")

  # v1 requires the outside option
  sim_no <- simulate_hmnp_data(N = 15, T = 2, J = 3, seed = 9,
                               include_outside = FALSE)
  d_no <- prepare_hmnp_data(sim_no$data, "task", "alt", "choice",
                            c("x1", "x2"), person_col = "pid",
                            include_outside_option = FALSE)
  expect_error(run_hmnprobit(input_data = d_no),
               "include_outside_option = TRUE")
})
