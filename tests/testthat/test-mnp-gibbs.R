# Tests for the mnp_gibbs() engine and the run_mnprobit() wrapper:
# output structure, reproducibility, thread invariance, validation, and the
# binary (J = 2) special case.

make_gibbs_args <- function(d, R = 200, burn = 50, thin = 1, seed = 123) {
  list(
    X = d$X, y = d$y, p = d$p,
    beta_bar = rep(0, d$K), A = 0.01 * diag(d$K),
    nu = d$p + 3, V = (d$p + 3) * diag(d$p),
    R = R, burn = burn, thin = thin, seed = seed
  )
}

vech_row_to_mat <- function(v, p) {
  S <- matrix(0, p, p)
  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      S[i, j] <- v[k]
      S[j, i] <- v[k]
      k <- k + 1L
    }
  }
  S
}

test_that("mnp_gibbs returns well-formed draws", {
  dt <- create_small_mnl_data()
  d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  out <- do.call(mnp_gibbs, make_gibbs_args(d))

  expect_equal(out$R_keep, 150L)
  expect_equal(dim(out$betadraw), c(150L, d$K))
  expect_equal(dim(out$sigmadraw), c(150L, d$p * (d$p + 1) / 2))
  expect_true(all(is.finite(out$betadraw)))
  expect_true(all(is.finite(out$sigmadraw)))
  expect_true(all(out$sigmadraw[, 1] > 0))

  # Every reconstructed Sigma draw is symmetric PD, and the storage order
  # round-trips through vech_row()
  for (r in c(1L, 50L, 150L)) {
    S <- vech_row_to_mat(out$sigmadraw[r, ], d$p)
    expect_true(all(eigen(S, only.values = TRUE)$values > 0))
    expect_equal(unname(choicer:::vech_row(S)), unname(out$sigmadraw[r, ]))
  }
})

test_that("mnp_gibbs honours thinning", {
  dt <- create_small_mnl_data()
  d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  out <- do.call(mnp_gibbs, make_gibbs_args(d, R = 200, burn = 50, thin = 3))
  expect_equal(out$R_keep, 50L)          # ceil(150 / 3)
  expect_equal(nrow(out$betadraw), 50L)
})

test_that("mnp_gibbs draws are seed-reproducible and thread-invariant", {
  dt <- create_small_mnl_data()
  d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  args <- make_gibbs_args(d)

  set_num_threads(1)
  out1 <- do.call(mnp_gibbs, args)
  set_num_threads(4)
  out4 <- do.call(mnp_gibbs, args)

  # Bitwise identical across thread counts: each (iteration, observation)
  # task has its own RNG stream, so scheduling cannot affect the draws
  expect_identical(out1$betadraw, out4$betadraw)
  expect_identical(out1$sigmadraw, out4$sigmadraw)

  # Same seed twice -> identical; different seed -> different
  out_again <- do.call(mnp_gibbs, args)
  expect_identical(out4$betadraw, out_again$betadraw)
  args2 <- args
  args2$seed <- 999
  expect_false(identical(do.call(mnp_gibbs, args2)$betadraw, out1$betadraw))
})

test_that("mnp_gibbs validates inputs", {
  dt <- create_small_mnl_data()
  d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
  args <- make_gibbs_args(d)

  bad <- args; bad$y[1] <- 99L
  expect_error(do.call(mnp_gibbs, bad), "outside")

  bad <- args; bad$nu <- 1
  expect_error(do.call(mnp_gibbs, bad), "nu")

  bad <- args; bad$burn <- 300L
  expect_error(do.call(mnp_gibbs, bad), "burn")

  bad <- args; bad$X <- args$X[-1, , drop = FALSE]
  expect_error(do.call(mnp_gibbs, bad), "rows")

  bad <- args; bad$beta_bar <- rep(0, d$K + 1)
  expect_error(do.call(mnp_gibbs, bad), "beta_bar")

  bad <- args; bad$seed <- -1
  expect_error(do.call(mnp_gibbs, bad), "seed")
})

test_that("run_mnprobit returns a well-formed choicer_mnp object", {
  dt <- create_small_mnl_data()
  fit <- suppressMessages(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 mcmc = list(R = 200, burn = 50, seed = 11))
  )

  expect_s3_class(fit, "choicer_mnp")
  expect_false(inherits(fit, "choicer_fit"))
  expect_named(fit$coefficients, c("x1", "x2", "ASC_2", "ASC_3"))
  expect_equal(fit$mcmc$R_keep, 150L)
  expect_equal(dim(fit$draws$beta), c(150L, 4L))

  # Identified normalization: sigma_11 draws are exactly 1 and the posterior
  # mean Sigma has unit (1, 1) element
  expect_true(all(abs(fit$draws$sigma[, "Sigma_11"] - 1) < 1e-12))
  expect_equal(unname(fit$sigma[1, 1]), 1)

  # Eager posterior vcov of the identified draws
  expect_equal(unname(vcov(fit)), unname(cov(fit$draws$beta)))
  expect_equal(nobs(fit), 20L)

  expect_output(print(fit), "Bayesian Multinomial Probit")
  s <- summary(fit)
  expect_s3_class(s, "summary.choicer_mnp")
  expect_named(s$coefficients, c("Mean", "SD", "CI_lo", "Median", "CI_hi"))
  expect_equal(nrow(s$coefficients), 4L)
  expect_equal(nrow(s$sigma_table), 3L)
  expect_output(print(s), "per-draw normalization by sigma_11")
})

test_that("run_mnprobit is reproducible via set.seed", {
  dt <- create_small_mnl_data()
  set.seed(31)
  f1 <- suppressMessages(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 mcmc = list(R = 100, burn = 10))
  )
  set.seed(31)
  f2 <- suppressMessages(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 mcmc = list(R = 100, burn = 10))
  )
  expect_identical(f1$draws$beta, f2$draws$beta)
  expect_identical(f1$mcmc$seed, f2$mcmc$seed)
})

test_that("run_mnprobit validates prior and mcmc settings", {
  dt <- create_small_mnl_data()
  expect_error(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 prior = list(nu = 1), mcmc = list(R = 100, burn = 10)),
    "nu"
  )
  expect_error(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 prior = list(A = diag(3)), mcmc = list(R = 100, burn = 10)),
    "prior\\$A"
  )
  expect_error(
    run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                 mcmc = list(R = 100, burn = 200)),
    "greater than"
  )
  expect_error(run_mnprobit(), "Supply either")
  expect_error(
    run_mnprobit(data = dt, input_data = list()),
    "not both"
  )
})

test_that("J = 2 reduces to binary probit", {
  set.seed(5)
  N <- 400
  dt <- data.table(id = rep(1:N, each = 2), alt = rep(1:2, N))
  dt[, x1 := rnorm(.N)]
  xd <- matrix(dt$x1, nrow = 2)            # 2 x N
  w <- (xd[2, ] - xd[1, ]) * 1.0 + rnorm(N)
  chosen <- ifelse(w > 0, 2L, 1L)
  dt[, choice := as.integer(alt == chosen[id])]

  fit <- suppressMessages(
    run_mnprobit(dt, "id", "alt", "choice", "x1", use_asc = FALSE,
                 mcmc = list(R = 2000, burn = 500, seed = 1))
  )

  # With p = 1, the identified Sigma draws are identically 1
  expect_true(all(abs(fit$draws$sigma - 1) < 1e-12))

  # Posterior mean under a diffuse prior is close to the probit MLE
  d <- prepare_mnp_data(dt, "id", "alt", "choice", "x1", use_asc = FALSE)
  g <- suppressWarnings(
    glm(d$y ~ 0 + d$X, family = binomial(link = "probit"))
  )
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 0.1)
})
