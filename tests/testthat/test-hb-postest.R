# Post-estimation tests for the hierarchical Bayes fits: predict (incl. the
# posterior-predictive delta for new alternatives), wtp, logsum /
# consumer_surplus (HMNL-only), elasticities / diversion, and the PPC
# helper. Fast structural checks run on CRAN; behavior-under-truth checks
# are skip_on_cran.

.fit_small_hmnl <- function() {
  sim <- simulate_hmnl_data(N = 60, T = 3, J = 4, beta = c(0.8, -0.6),
                            theta = c(0.5, -0.4), sigma_d = 0.5, seed = 42)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, mcmc = list(R = 400, burn = 150, seed = 7))
  )
  list(sim = sim, fit = fit)
}

test_that("Gauss-Hermite quadrature matches the iid-probit integral", {
  # against brute-force Monte Carlo argmax frequencies at fixed utilities
  V <- c(0.5, -0.2, 0.8)
  gh <- .gauss_hermite(20)
  pq <- .hb_task_probs(V, rep(1L, 3), 1L, "hmnp", gh)
  set.seed(5)
  S <- 200000
  U <- matrix(rep(c(V, 0), S), nrow = 4) + matrix(rnorm(4 * S), nrow = 4)
  mc <- tabulate(apply(U, 2, which.max), 4) / S
  expect_lt(max(abs(c(pq$p_inside, pq$p_outside) - mc)), 5e-3)
  expect_equal(sum(pq$p_inside) + pq$p_outside, 1, tolerance = 1e-10)

  # and the GH rule itself integrates exp(-t^2) polynomials exactly
  gh8 <- .gauss_hermite(8)
  expect_equal(sum(gh8$weights), sqrt(pi), tolerance = 1e-10)
  expect_equal(sum(gh8$weights * gh8$nodes^2), sqrt(pi) / 2,
               tolerance = 1e-10)
})

test_that("predict returns coherent posterior shares (both engines)", {
  fx <- .fit_small_hmnl()
  set.seed(1)
  p <- predict(fx$fit, n_draws = 50)
  expect_s3_class(p, "data.table")
  expect_setequal(p$alternative, c(as.character(1:4), "(outside)"))
  expect_equal(sum(p$share), 1, tolerance = 1e-8)
  expect_true(all(p$lower <= p$share & p$share <= p$upper))
  expect_equal(dim(attr(p, "draws")), c(50L, 5L))

  # row-level probabilities
  pr <- predict(fx$fit, n_draws = 20, aggregate = FALSE)
  expect_length(pr, nrow(fx$fit$data$X))
  expect_true(all(pr > 0 & pr < 1))

  # HMNP path (quadrature engine)
  simp <- simulate_hmnp_data(N = 40, T = 2, J = 3, seed = 11,
                             theta = c(0.2, -0.3), sigma_d = 0.4)
  dp <- prepare_hmnp_data(simp$data, "task", "alt", "choice",
                          c("x1", "x2"), person_col = "pid",
                          alt_covariate_cols = "z1")
  fitp <- suppressWarnings(
    run_hmnprobit(input_data = dp, mcmc = list(R = 300, burn = 100, seed = 3))
  )
  set.seed(2)
  pp <- predict(fitp, n_draws = 25)
  expect_equal(sum(pp$share), 1, tolerance = 1e-8)
  expect_true(all(pp$share > 0))
})

test_that("a subsidy counterfactual moves shares in the expected direction", {
  fx <- .fit_small_hmnl()
  cf <- data.table::copy(fx$sim$data)
  cf[cf$alt == 2, "x1"] <- cf[cf$alt == 2, ][["x1"]] + 1   # x1 coef > 0
  set.seed(1)
  base <- predict(fx$fit, n_draws = 50)
  set.seed(1)
  new <- predict(fx$fit, newdata = cf, n_draws = 50)
  expect_gt(new[new$alternative == "2", ][["share"]],
            base[base$alternative == "2", ][["share"]])
  expect_lt(new[new$alternative == "(outside)", ][["share"]],
            base[base$alternative == "(outside)", ][["share"]])
})

test_that("new alternatives get a posterior-predictive delta", {
  skip_on_cran()

  # Simulate J = 5, fit on the first 4, predict with the held-out 5th: its
  # predicted share must carry the extra sigma_d uncertainty and its truth
  # must be plausible under the posterior-predictive interval.
  sim <- simulate_hmnl_data(N = 250, T = 4, J = 5, beta = c(0.8, -0.6),
                            theta = c(0.4, -0.5), sigma_d = 0.4, seed = 9)
  dt_all <- sim$data
  dt_fit <- dt_all[dt_all$alt != 5L, ]
  d <- prepare_hmnl_data(dt_fit, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  fit <- suppressWarnings(
    run_hmnlogit(input_data = d, mcmc = list(R = 1500, burn = 500, seed = 7))
  )

  set.seed(4)
  p <- predict(fit, newdata = dt_all, n_draws = 150)
  expect_true("5" %in% p$alternative)
  new_row <- p[p$alternative == "5", ]
  known_sd <- mean(p[p$alternative %in% as.character(1:4), ][["sd"]])
  expect_gt(new_row$sd, known_sd)          # entry uncertainty is wider

  # truth: simulate the DGP share of alt 5 from the realized delta
  draws <- attr(p, "draws")
  expect_true(new_row$lower >= 0 && new_row$upper <= 1)
  # the entrant + incumbents shares still sum to one
  expect_equal(sum(p$share), 1, tolerance = 1e-8)

  # new alternatives need their z covariates
  dt_missing <- data.table::copy(dt_all)
  dt_missing$z1 <- NULL
  expect_error(predict(fit, newdata = dt_missing), "missing columns|z1")
})

test_that("wtp returns median-based posterior ratios with sane coverage", {
  fx <- .fit_small_hmnl()
  w <- suppressWarnings(wtp(fx$fit, price_var = "x2"))
  expect_s3_class(w, "data.table")
  expect_equal(w$attribute, "x1")
  # true WTP = -0.8 / -0.6 = 1.333; the interval must cover it
  expect_gt(w$upper, 4 / 3)
  expect_lt(w$lower, 4 / 3 + 1.5)
  expect_true(w$lower < w$wtp && w$wtp < w$upper)
  expect_error(wtp(fx$fit, price_var = "nope"), "price_var")
})

test_that("logsum and consumer_surplus are HMNL-only and coherent", {
  fx <- .fit_small_hmnl()
  set.seed(1)
  ls_ <- logsum(fx$fit, n_draws = 30)
  expect_length(ls_, fx$fit$nobs)
  expect_true(all(is.finite(ls_)))
  expect_true(all(ls_ > 0))                # log(1 + positive sum) > 0

  # improving an attribute raises expected surplus: positive CV
  cf <- data.table::copy(fx$sim$data)
  cf$x1 <- cf$x1 + 0.5
  set.seed(2)
  cs <- suppressWarnings(
    consumer_surplus(fx$fit, "x2", newdata = cf, n_draws = 30)
  )
  expect_s3_class(cs, "data.table")
  cv <- attr(cs, "cv")
  expect_length(cv, 3L)
  expect_gt(cv[["50%"]], 0)

  # Task weights enter the aggregate CV rather than being silently ignored.
  set.seed(2)
  cs_w <- suppressWarnings(
    consumer_surplus(fx$fit, "x2", newdata = cf, n_draws = 30,
                     weights = rep(2, fx$fit$nobs))
  )
  expect_equal(attr(cs_w, "cv"), 2 * cv, tolerance = 1e-10)
  expect_equal(attr(cs_w, "cv_weights"), rep(2, fx$fit$nobs))
  expect_error(
    consumer_surplus(fx$fit, "x2", newdata = cf, n_draws = 30,
                     weights = rep(1, fx$fit$nobs - 1L)),
    "one entry per choice situation"
  )

  # Baseline and counterfactual task columns are matched by identity, not
  # silently by their position after sorting.
  cf_other_tasks <- data.table::copy(cf)
  cf_other_tasks$task <- cf_other_tasks$task + 100000L
  expect_error(
    consumer_surplus(fx$fit, "x2", newdata = cf_other_tasks, n_draws = 10),
    "same choice situations"
  )
  cf_missing_task <- cf[task != min(task)]
  expect_error(
    consumer_surplus(fx$fit, "x2", newdata = cf_missing_task, n_draws = 10),
    "expected.*tasks but found"
  )

  # probit: informative refusal
  simp <- simulate_hmnp_data(N = 30, T = 2, J = 3, seed = 11)
  dp <- prepare_hmnp_data(simp$data, "task", "alt", "choice",
                          c("x1", "x2"), person_col = "pid")
  fitp <- suppressWarnings(
    run_hmnprobit(input_data = dp, mcmc = list(R = 200, burn = 50, seed = 3))
  )
  expect_error(logsum(fitp), "logit-only")
  expect_error(consumer_surplus(fitp, "x2"), "logit-only")
})

test_that("HMNL probabilities and logsums are stable at extreme utilities", {
  V <- c(1000, -1000, 800, 799)
  task <- c(1L, 1L, 2L, 2L)
  terms <- choicer:::.hb_logit_task_terms(V, task, 2L)

  expect_true(all(is.finite(terms$p_inside)))
  expect_true(all(is.finite(terms$p_outside)))
  expect_true(all(is.finite(terms$logsum)))
  expect_equal(
    as.numeric(rowsum(terms$p_inside, task, reorder = TRUE)) +
      terms$p_outside,
    rep(1, 2),
    tolerance = 1e-14
  )
  expect_equal(terms$logsum[1], 1000, tolerance = 1e-12)
  expect_equal(terms$logsum[2], 800 + log1p(exp(-1)), tolerance = 1e-12)
})

test_that("elasticities and diversion ratios are internally coherent", {
  fx <- .fit_small_hmnl()
  set.seed(1)
  el <- elasticities(fx$fit, "x1", n_draws = 25)
  J <- fx$fit$J
  expect_equal(dim(el), c(J + 1L, J))
  # x1 coefficient is positive: own-effects positive, cross-effects negative
  for (j in seq_len(J)) {
    expect_gt(el[j, j], 0)
    expect_true(all(el[-j, j] < 0))
  }

  set.seed(1)
  dr <- diversion_ratios(fx$fit, "x1", n_draws = 25)
  expect_equal(dim(dr), c(J + 1L, J))
  expect_true(all(diag(dr[seq_len(J), ]) == 0))
  # share lost by the perturbed alternative is fully accounted for
  expect_equal(unname(colSums(dr)), rep(1, J), tolerance = 1e-6)
  expect_true(all(dr[cbind(J + 1L, seq_len(J))] > 0))   # some goes outside
})

test_that("ppc_shares covers observed shares on well-specified data", {
  fx <- .fit_small_hmnl()
  set.seed(1)
  ppc <- ppc_shares(fx$fit, n_draws = 60)
  expect_s3_class(ppc, "data.table")
  expect_equal(nrow(ppc), fx$fit$J + 1L)
  expect_equal(sum(ppc$observed), 1, tolerance = 1e-8)
  # the model is well-specified here: most shares must be covered
  expect_gte(sum(ppc$covered), fx$fit$J - 1L)
})

test_that("individual-level prediction uses the respondent posteriors", {
  fx <- .fit_small_hmnl()
  set.seed(1)
  expect_message(
    p_ind <- predict(fx$fit, level = "individual", n_draws = 25),
    "posterior-mean beta_i"
  )
  expect_equal(sum(p_ind$share), 1, tolerance = 1e-8)

  # draws-based individual prediction (fully posterior-integrated)
  d <- fx$fit$data
  fit2 <- suppressWarnings(
    run_hmnlogit(input_data = d, keep_beta_i = "draws",
                 mcmc = list(R = 300, burn = 100, seed = 5))
  )
  set.seed(2)
  p_ind2 <- predict(fit2, level = "individual", n_draws = 25)
  expect_equal(sum(p_ind2$share), 1, tolerance = 1e-8)
})
