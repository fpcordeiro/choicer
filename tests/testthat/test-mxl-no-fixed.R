# A mixed logit whose every coefficient is random (no fixed covariates): the
# preparations' empty X, and the kernels against the same model with an extra
# fixed column held at 0, an oracle by which the likelihood, its derivatives
# and the post-estimation results must agree.

nofix_data <- function(N = 120L, J = 3L, seed = 2L) {
  set.seed(seed)
  id <- rep(seq_len(N), each = J)
  w1 <- rnorm(N * J)
  w2 <- rnorm(N * J)
  b <- rep(rnorm(N, 1, 0.7), each = J)
  u <- b * w1 - 0.5 * w2 - log(-log(runif(N * J)))
  data.table(id = id, alt = rep(seq_len(J), N), w1 = w1, w2 = w2, x = rnorm(N * J),
             person = (id - 1L) %/% 4L, nest = rep(c(1L, 1L, 2L), N),
             choice = as.integer(u == ave(u, id, FUN = max)))
}

test_that("the preparations return an empty X, row for row, without covariates", {
  d <- nofix_data()
  n <- nrow(d)
  for (cols in list(NULL, character(0))) {
    for (input in list(prepare_mxl_data(d, "id", "alt", "choice", cols, "w1"),
                       prepare_mnl_data(d, "id", "alt", "choice", cols),
                       prepare_nl_data(d, "id", "alt", "choice", cols, "nest"))) {
      expect_identical(dim(input$X), c(n, 0L))
      expect_type(input$X, "double")
      expect_length(input$alt_idx, n)
    }
  }
  # the multinomial probit's design holds the constants then
  mnp <- prepare_mnp_data(d, "id", "alt", "choice", NULL)
  expect_identical(ncol(mnp$X), 2L)
})

test_that("an all-random mixed logit equals the model with a fixed column held at 0", {
  d <- nofix_data()
  args <- list(data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
               random_var_cols = c("w1", "w2"), rc_mean = TRUE, S = 40L, seed = 9L)
  fit0 <- suppressMessages(do.call("run_mxlogit", c(args, list(covariate_cols = NULL,
                                                               polish = "newton"))))
  expect_true(fit0$optimizer$polish$converged)
  expect_identical(fit0$param_map$beta, integer(0))
  expect_identical(names(coef(fit0)),
                   c("Mu_w1", "Mu_w2", "L_11", "L_22", "ASC_2", "ASC_3"))
  # character(0) is the same call
  fit_e <- suppressMessages(do.call("run_mxlogit", c(args, list(covariate_cols = character(0),
                                                                polish = "newton"))))
  expect_equal(coef(fit_e), coef(fit0), tolerance = 1e-8)

  held <- suppressMessages(do.call("run_mxlogit", c(args, list(
    covariate_cols = "x", theta_init = c(0, unname(coef(fit0))),
    lower = c(x = 0), upper = c(x = 0), polish = "newton"))))
  expect_equal(unname(coef(held))[-1L], unname(coef(fit0)), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(held)), as.numeric(logLik(fit0)), tolerance = 1e-9)

  # The kernels at one point, one thread, draws in generate mode as the fits'
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  i0 <- prepare_mxl_data(d, "id", "alt", "choice", NULL, c("w1", "w2"))
  ix <- prepare_mxl_data(d, "id", "alt", "choice", "x", c("w1", "w2"))
  eta <- get_halton_normals(40L, length(i0$M), 2L)
  th <- unname(coef(fit0)) + 0.05
  common <- function(fn, input, theta, ...) {
    fn(theta = theta, X = input$X, W = input$W, alt_idx = input$alt_idx,
       choice_idx = input$choice_idx, M = input$M, weights = input$weights,
       eta_draws = eta, rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = TRUE,
       use_asc = TRUE, include_outside_option = FALSE, ...)
  }
  g0 <- common(choicer:::mxl_loglik_gradient_parallel, i0, th)
  gx <- common(choicer:::mxl_loglik_gradient_parallel, ix, c(0, th))
  expect_equal(g0$objective, gx$objective, tolerance = 1e-12)
  expect_equal(as.numeric(g0$gradient), as.numeric(gx$gradient)[-1L], tolerance = 1e-10)
  h0 <- common(choicer:::mxl_hessian_parallel, i0, th)
  hx <- common(choicer:::mxl_hessian_parallel, ix, c(0, th))
  expect_equal(unname(h0), unname(hx)[-1L, -1L], tolerance = 1e-10)
  b0 <- common(choicer:::mxl_bhhh_parallel, i0, th)
  bx <- common(choicer:::mxl_bhhh_parallel, ix, c(0, th))
  expect_equal(unname(b0), unname(bx)[-1L, -1L], tolerance = 1e-10)

  # Post-estimation: predictions and the elasticities of a random coefficient's
  # variable, all read at the same point
  expect_equal(predict(fit0, type = "shares"), predict(held, type = "shares"), tolerance = 1e-6)
  expect_equal(elasticities(fit0, elast_var = "w1", is_random_coef = TRUE),
               elasticities(held, elast_var = "w1", is_random_coef = TRUE), tolerance = 1e-6)
  expect_equal(diversion_ratios(fit0, wrt_var = "w1", is_random_coef = TRUE),
               diversion_ratios(held, wrt_var = "w1", is_random_coef = TRUE), tolerance = 1e-6)
  # the variance routes exist, and a fixed variable is not there to ask for
  for (type in c("hessian", "bhhh", "robust")) {
    expect_identical(dim(vcov(fit0, type = type)), c(6L, 6L))
  }
  expect_error(elasticities(fit0, elast_var = "w1", is_random_coef = FALSE), "not found")
})

test_that("an all-random mixed logit fits in every layout", {
  d <- nofix_data()
  fit <- function(...) {
    suppressMessages(run_mxlogit(d, "id", "alt", "choice", NULL, c("w1", "w2"), S = 30L,
                                 seed = 4L, ...))
  }
  ok <- function(f) {
    expect_true(is.finite(f$loglik))
    expect_true(all(is.finite(f$se)))
  }
  ok(fit())
  ok(fit(rc_mean = TRUE, person_col = "person"))
  ok(fit(rc_mean = TRUE, rc_correlation = TRUE))
  ok(fit(rc_mean = TRUE, rc_dist = c(1L, 0L), use_asc = FALSE,
         theta_init = c(-0.5, 0, -1, -1)))
  ok(suppressWarnings(fit(rc_mean = TRUE, draws = "store")))
  ok(fit(rc_mean = TRUE, se_method = "cluster", cluster_col = "person"))
  ok(fit(rc_mean = TRUE, concentrate_out_asc = "profile_likelihood"))
  ok(fit(rc_mean = TRUE, concentrate_out_asc = "blp"))
})
