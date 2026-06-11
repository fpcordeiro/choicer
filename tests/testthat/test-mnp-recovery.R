# Parameter recovery on simulated MNP data with correlated differenced
# errors. The correlated Sigma exercises the cross terms of the latent
# truncated-normal conditionals: a sign or indexing error in the
# conditional-mean formula would systematically distort both beta and Sigma.

test_that("MNP recovers identified parameters on simulated data", {
  skip_on_cran()

  sim <- create_mnp_sim_data(N = 800, J = 3,
                             beta = c(1.0, -0.5),
                             Sigma = matrix(c(1.0, 0.5, 0.5, 1.5), 2, 2),
                             seed = 99)
  fit <- suppressMessages(
    run_mnprobit(sim$data, "id", "alt", "choice", sim$covariate_cols,
                 use_asc = TRUE, mcmc = list(R = 4000, burn = 1000, seed = 7))
  )

  est <- coef(fit)

  # Covariate coefficients near the true identified values
  expect_lt(max(abs(est[sim$covariate_cols] - sim$beta_id)), 0.15)

  # ASCs near zero (true DGP has none)
  expect_lt(max(abs(est[c("ASC_2", "ASC_3")])), 0.25)

  # Identified Sigma: unit (1, 1) element by construction; covariance
  # structure near the truth
  expect_equal(unname(fit$sigma[1, 1]), 1)
  expect_lt(abs(fit$sigma[2, 1] - sim$Sigma_id[2, 1]), 0.2)
  expect_lt(abs(fit$sigma[2, 2] - sim$Sigma_id[2, 2]), 0.35)

  # Posterior uncertainty is sensible: positive SDs, finite draws
  expect_true(all(fit$se > 0))
  expect_true(all(is.finite(fit$draws$beta)))

  # Method smoke tests on a real fit
  s <- summary(fit)
  expect_equal(nrow(s$coefficients), 4L)
  expect_equal(nrow(s$sigma_table), 3L)
  expect_equal(nobs(fit), 800L)
  expect_equal(dim(vcov(fit)), c(4L, 4L))
})
