# Parameter-recovery tests for run_hmnprobit on simulated data. Not run on
# CRAN. simulate_hmnp_data stores truth on the identified (/sigma) scale, so
# posterior summaries compare directly.

test_that("run_hmnprobit recovers the DGP on a mid-sized panel", {
  skip_on_cran()

  # J = 12 (not 8): the delta mean-function intercept theta_0 is the
  # weakly-identified block (regressing the J alternative effects on Z), and
  # at J = 8 its posterior mean sits right at the 3-SD tolerance for this
  # seed. The Gibbs chain is bitwise-reproducible across thread counts but
  # not across platform BLAS/LAPACK (the master-thread Wishart / inv_sympd
  # calls), so a borderline theta_0 that just passes on one platform's chain
  # realization can just fail on another's. J = 12 gives theta_0 enough
  # alternatives to concentrate near truth with margin on every platform.
  sim <- simulate_hmnp_data(N = 300, T = 5, J = 12, beta = c(0.8, -0.6),
                            theta = c(-0.5, -0.4), sigma_d = 0.5, sigma = 1,
                            seed = 42)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnprobit(input_data = d, mcmc = list(R = 2500, burn = 800, seed = 7))
  )
  tp <- sim$true_params

  b_hat <- coef(fit)
  expect_true(all(abs(b_hat - tp$beta) < pmax(3 * fit$se, 0.15)))

  w_hat <- diag(fit$W_mean)
  w_true <- diag(tp$W)
  diag_idx <- vapply(1:2, function(k) (k * (k + 1L)) %/% 2L, integer(1L))
  w_sd <- apply(fit$draws$w_vech, 2, stats::sd)[diag_idx]
  expect_true(all(abs(w_hat - w_true) < pmax(0.5 * w_true, 3 * w_sd)))

  th_hat <- colMeans(fit$draws$theta)
  th_sd <- apply(fit$draws$theta, 2, stats::sd)
  expect_true(all(abs(th_hat - tp$theta) < 3 * th_sd + 0.05))

  sd_hat <- mean(sqrt(fit$draws$sigma_d2))
  expect_gt(sd_hat, 0.3 * tp$sigma_d)
  expect_lt(sd_hat, 3.0 * tp$sigma_d)

  expect_gt(cor(colMeans(fit$draws$delta), tp$delta), 0.9)

  rt <- recovery_table(fit, sim)
  expect_setequal(unique(rt$group),
                  c("beta", "w", "theta", "sigma_d", "delta"))
  expect_true(all(is.finite(rt$estimate)))
})

test_that("truth is recovered on the identified scale when sigma != 1", {
  skip_on_cran()

  # DGP scale sigma = 2: the raw chains drift on a different scale, but the
  # per-draw-normalized posteriors must match the identified truth
  # (beta / sigma etc.) that simulate_hmnp_data reports.
  sim <- simulate_hmnp_data(N = 300, T = 5, J = 5, beta = c(1.2, -0.9),
                            theta = c(-0.6, 0.3), sigma_d = 0.6, sigma = 2,
                            seed = 11)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnprobit(input_data = d, mcmc = list(R = 2500, burn = 800, seed = 3))
  )
  tp <- sim$true_params
  expect_equal(unname(tp$beta), c(1.2, -0.9) / 2)   # identified-scale truth

  expect_true(all(abs(coef(fit) - tp$beta) < pmax(3 * fit$se, 0.15)))
  expect_gt(cor(colMeans(fit$draws$delta), tp$delta), 0.9)
})

test_that("J = 1 vs outside matches sqrt(2) x glm probit (degenerate DGP)", {
  skip_on_cran()

  # One inside alternative against the stochastic outside: a binary probit
  # whose composite error has variance 2 sigma^2, so the identified HMNP
  # coefficients equal sqrt(2) times the glm probit coefficients. Tight
  # priors pin W ~ 0 and sigma_d ~ 0 (matching the degenerate DGP); the
  # near-degenerate hierarchy mixes slowly (the documented funnel regime),
  # hence the long chain and the loose tolerance.
  sim <- simulate_hmnp_data(N = 3000, T = 1, J = 1, beta = c(0.9, -0.7),
                            theta = c(0.4), sigma_d = 0, sigma = 1,
                            W = diag(1e-8, 2), seed = 42)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"))
  fit <- suppressWarnings(run_hmnprobit(
    input_data = d, keep_beta_i = "none",
    prior = list(nu = 500, V = 0.5 * diag(2),
                 sd_prior = list(half_cauchy = FALSE, c0 = 500, d0 = 0.5)),
    mcmc = list(R = 20000, burn = 8000, seed = 7)
  ))
  glm_fit <- stats::glm(choice ~ x1 + x2,
                        family = stats::binomial(link = "probit"),
                        data = sim$data)

  expect_true(all(abs(coef(fit) -
                        sqrt(2) * coef(glm_fit)[c("x1", "x2")]) < 0.08))
  expect_lt(abs(coef(fit, "delta") - sqrt(2) * coef(glm_fit)[1]), 0.08)
})

test_that("cross-sectional mode (T = 1, person_col = NULL) recovers b", {
  skip_on_cran()

  sim <- simulate_hmnp_data(N = 1500, T = 1, J = 5, beta = c(0.8, -0.6),
                            theta = c(-0.4, -0.3), sigma_d = 0.5, seed = 31)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         alt_covariate_cols = "z1")   # person_col = NULL
  expect_true(all(d$Ti == 1L))

  fit <- suppressWarnings(
    run_hmnprobit(input_data = d, keep_beta_i = "none",
                  mcmc = list(R = 2000, burn = 600, seed = 17))
  )
  expect_true(all(abs(coef(fit) - sim$true_params$beta) <
                    pmax(3 * fit$se, 0.2)))
  expect_gt(cor(colMeans(fit$draws$delta), sim$true_params$delta), 0.85)
})
