# Parameter-recovery tests for run_hmnlogit on simulated data. Not run on
# CRAN (multi-second chains); tolerances are loose by design -- these guard
# against systematic estimator bugs, not Monte Carlo noise.

test_that("run_hmnlogit recovers the DGP on a mid-sized panel", {
  skip_on_cran()

  sim <- simulate_hmnl_data(N = 200, T = 8, J = 5, beta = c(0.8, -0.6),
                            theta = c(0.5, -0.4), sigma_d = 0.5, seed = 42)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, mcmc = list(R = 2500, burn = 800, seed = 7))
  )
  tp <- sim$true_params

  # b within 3 posterior SDs (and sane absolute bias)
  b_hat <- coef(fit)
  b_sd <- fit$se
  expect_true(all(abs(b_hat - tp$beta) < pmax(3 * b_sd, 0.15)))

  # diag(W): posterior mean within 50% of truth or 3 posterior SDs
  # (parentheses matter: %/% binds tighter than *)
  diag_idx <- vapply(1:2, function(k) (k * (k + 1L)) %/% 2L, integer(1L))
  w_hat <- colMeans(fit$draws$w_vech)[diag_idx]
  w_sd <- apply(fit$draws$w_vech, 2, stats::sd)[diag_idx]
  w_true <- diag(tp$W)
  expect_true(all(abs(w_hat - w_true) < pmax(0.5 * w_true, 3 * w_sd)))

  # theta within 3 posterior SDs (J = 5 alternatives -> wide posteriors)
  th_hat <- colMeans(fit$draws$theta)
  th_sd <- apply(fit$draws$theta, 2, stats::sd)
  expect_true(all(abs(th_hat - tp$theta) < 3 * th_sd + 0.05))

  # sigma_d in the right order of magnitude
  sd_hat <- mean(sqrt(fit$draws$sigma_d2))
  expect_gt(sd_hat, 0.3 * tp$sigma_d)
  expect_lt(sd_hat, 3.0 * tp$sigma_d)

  # delta shape recovered (level is anchored by the outside share)
  expect_gt(cor(colMeans(fit$draws$delta), tp$delta), 0.9)

  # recovery_table runs and covers all five blocks
  rt <- recovery_table(fit, sim)
  expect_setequal(unique(rt$group),
                  c("beta", "w", "theta", "sigma_d", "delta"))
  expect_true(all(is.finite(rt$estimate)))
})

test_that("run_hmnlogit recovers log-normal coordinates on the chain scale", {
  skip_on_cran()

  sim <- simulate_hmnl_data(N = 250, T = 8, J = 4, beta = c(0.3, -0.6),
                            rc_dist = c(1L, 0L), theta = c(0.4, 0),
                            sigma_d = 0.4, seed = 11)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", rc_dist = c(1L, 0L))
  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, rc_dist = NULL,
                 mcmc = list(R = 2500, burn = 800, seed = 3))
  )

  # truth is on the chain (log) scale for the log-normal coordinate
  b_hat <- coef(fit)
  b_sd <- fit$se
  expect_true(all(abs(b_hat - sim$true_params$beta) < pmax(3 * b_sd, 0.2)))
})

test_that("degenerate DGP (no heterogeneity) concentrates at the truth", {
  skip_on_cran()

  # W ~ 0 and sigma_d = 0: plain MNL with delta = Z theta. Tight priors on
  # both variance blocks make the fit collapse to the pooled model, so b
  # must sit close to the DGP truth.
  sim <- simulate_hmnl_data(N = 400, T = 5, J = 4, beta = c(0.8, -0.6),
                            W = diag(1e-6, 2), theta = c(0.3, -0.2),
                            sigma_d = 0, seed = 21)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(run_hmnlogit(
    input_data = d,
    prior = list(nu = 500, V = 0.5 * diag(2),
                 sd_prior = list(half_cauchy = FALSE, c0 = 500, d0 = 0.5)),
    mcmc = list(R = 2000, burn = 600, seed = 13)
  ))

  expect_true(all(abs(coef(fit) - c(0.8, -0.6)) < 0.1))
  # posterior W is pinned near the tight prior mean (~0.001)
  diag_idx <- c(1L, 3L)
  expect_true(all(colMeans(fit$draws$w_vech)[diag_idx] < 0.01))
  # delta collapses to Z theta: xi posterior means ~ 0
  expect_true(all(abs(fit$xi$mean) < 0.1))
})

test_that("cross-sectional mode (T = 1, person_col = NULL) recovers b and theta", {
  skip_on_cran()

  sim <- simulate_hmnl_data(N = 1500, T = 1, J = 5, beta = c(0.8, -0.6),
                            theta = c(0.5, -0.4), sigma_d = 0.5, seed = 31)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         alt_covariate_cols = "z1")   # person_col = NULL
  expect_true(all(d$Ti == 1L))
  expect_equal(d$N_persons, 1500L)

  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, keep_beta_i = "none",
                 mcmc = list(R = 1500, burn = 500, seed = 17))
  )
  b_hat <- coef(fit)
  expect_true(all(abs(b_hat - sim$true_params$beta) < pmax(3 * fit$se, 0.2)))
  th_hat <- colMeans(fit$draws$theta)
  th_sd <- apply(fit$draws$theta, 2, stats::sd)
  expect_true(all(abs(th_hat - sim$true_params$theta) < 3 * th_sd + 0.05))
  expect_gt(cor(colMeans(fit$draws$delta), sim$true_params$delta), 0.85)
})
