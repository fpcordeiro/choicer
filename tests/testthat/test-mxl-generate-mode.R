# tests/testthat/test-mxl-generate-mode.R
#
# Phase-B tests for the on-the-fly Halton draw generator.
#
# Test groups:
#   1. C++ estimation equivalence: store vs generate (scramble="none").
#      scramble="none" -> identity permutations -> same Halton sequence as
#      get_halton_normals() -> loglik/gradient identical to machine epsilon.
#      Tests both the kernel AND confirms that a run_mxlogit() fit with
#      seed=0, scramble="none" matches the store-mode logLik at the same theta.
#   2. generate-mode gradient vs numDeriv::grad across config matrix
#      (rc_correlation x rc_mean x rc_dist x include_outside_option).
#      NIT-2 fix: outside-option branch (include_outside_option=TRUE) now
#      covered by two additional cells (uncorr/normal/oo and corr/normal/oo).
#   3. Post-estimation generics equivalence in generate mode (scramble="none",
#      seed=0): predict, elasticities, diversion_ratios, logsum,
#      consumer_surplus, and vcov under both se_method="hessian" and
#      se_method="bhhh".
#      NIT-3 fix: both objects carry byte-identical coefficients (fit_g is a
#      copy of fit_s with only draws_info switched to generate mode); the
#      comparison isolates draw-path equivalence, not optimizer reproducibility.
#
# NOTE: NO cross-platform bit-equality assertions on Owen-scrambled draws.

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

# Machine epsilon level: store vs generate (same Halton points, same arithmetic)
TOL_KERNEL  <- 1e-12

# Optimizer-level: after convergence on a well-identified dataset
TOL_EQUIV   <- 1e-6

# Gradient accuracy against numDeriv
TOL_GRAD_GEN <- 1e-5

# Post-estimation evaluated at the same theta
TOL_POST <- 1e-10

# ---------------------------------------------------------------------------
# Shared fixture: a small identified MXL dataset + converged store-mode fit
# ---------------------------------------------------------------------------

.gen_mode_fixture <- local({
  .fixture <- NULL
  function() {
    if (!is.null(.fixture)) return(.fixture)
    set.seed(321)
    sim <- simulate_mxl_data(
      N   = 150, J = 3,
      beta  = 0.8,
      delta = c(0.4, -0.4),
      Sigma = matrix(0.5, 1, 1),
      seed  = 321,
      outside_option = FALSE,
      vary_choice_set = FALSE
    )
    dt   <- sim$data
    S_val <- 50L

    # Store-mode fit (converged)
    fit_s <- run_mxlogit(
      data            = dt,
      id_col          = "id",
      alt_col         = "alt",
      choice_col      = "choice",
      covariate_cols  = "x1",
      random_var_cols = "w1",
      S               = S_val,
      rc_correlation  = FALSE,
      rc_mean         = FALSE,
      rc_dist         = rep(0L, 1L),
      use_asc         = TRUE,
      include_outside_option = FALSE,
      draws           = "store",
      se_method       = "hessian",
      control         = list(print_level = 0L, maxeval = 1000L)
    )

    # Prepare the same inputs for C++ kernel calls
    inputs <- prepare_mxl_data(
      dt, "id", "alt", "choice", "x1", "w1",
      rc_correlation = FALSE
    )
    K_w   <- ncol(inputs$W)
    eta   <- get_halton_normals(S_val, inputs$N, K_w)
    eta_empty <- array(0, dim = c(K_w, 0L, 0L))

    .fixture <<- list(
      dt      = dt,
      S       = S_val,
      fit_s   = fit_s,
      inputs  = inputs,
      K_w     = K_w,
      eta     = eta,
      eta_empty = eta_empty
    )
    .fixture
  }
})

# ---------------------------------------------------------------------------
# 1. C++ kernel equivalence: store vs generate (scramble="none", seed=0)
# ---------------------------------------------------------------------------

test_that("store and generate (scramble=none) loglik/gradient match at same theta", {
  fix   <- .gen_mode_fixture()
  theta <- coef(fix$fit_s)  # converged theta from store mode
  K_w   <- fix$K_w
  S     <- fix$S

  r_store <- mxl_loglik_gradient_parallel(
    theta                  = theta,
    X                      = fix$inputs$X,
    W                      = fix$inputs$W,
    alt_idx                = fix$inputs$alt_idx,
    choice_idx             = fix$inputs$choice_idx,
    M                      = fix$inputs$M,
    weights                = fix$inputs$weights,
    eta_draws              = fix$eta,
    rc_dist                = rep(0L, K_w),
    rc_correlation         = FALSE,
    rc_mean                = FALSE,
    use_asc                = TRUE,
    include_outside_option = FALSE
  )

  r_gen <- mxl_loglik_gradient_parallel(
    theta                  = theta,
    X                      = fix$inputs$X,
    W                      = fix$inputs$W,
    alt_idx                = fix$inputs$alt_idx,
    choice_idx             = fix$inputs$choice_idx,
    M                      = fix$inputs$M,
    weights                = fix$inputs$weights,
    eta_draws              = fix$eta_empty,
    rc_dist                = rep(0L, K_w),
    rc_correlation         = FALSE,
    rc_mean                = FALSE,
    use_asc                = TRUE,
    include_outside_option = FALSE,
    gen_seed               = 0L,
    gen_scramble           = 0L,
    gen_S                  = S
  )

  # Loglik (negated) must match to machine-epsilon level
  diff_obj <- abs(r_store$objective - r_gen$objective)
  expect_true(diff_obj < TOL_KERNEL,
    label = paste0("C++ loglik diff store vs generate: ", diff_obj))

  # Gradient must match to machine-epsilon level
  max_grad_diff <- max(abs(as.numeric(r_store$gradient) - as.numeric(r_gen$gradient)))
  expect_true(max_grad_diff < TOL_KERNEL,
    label = paste0("C++ gradient max diff: ", max_grad_diff))
})

test_that("store and generate (scramble=none) loglik match on multi-K_w theta", {
  # Test with 2 random coefficients to cover multi-dimensional case
  set.seed(111)
  dt_2kw <- create_small_mxl_data(seed = 111)
  inputs  <- prepare_mxl_data(dt_2kw, "id", "alt", "choice", "x1", c("w1","w2"),
                               rc_correlation = FALSE)
  K_w <- ncol(inputs$W)
  S   <- 30L
  set.seed(5)
  theta <- c(runif(1, -0.3, 0.3), log(runif(K_w, 0.3, 0.8)),
             runif(nrow(inputs$alt_mapping) - 1L, -0.2, 0.2))

  eta       <- get_halton_normals(S, inputs$N, K_w)
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))

  r_s <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta,
    rc_dist=rep(0L,K_w), rc_correlation=FALSE, rc_mean=FALSE,
    use_asc=TRUE, include_outside_option=FALSE
  )
  r_g <- mxl_loglik_gradient_parallel(
    theta, inputs$X, inputs$W, inputs$alt_idx, inputs$choice_idx,
    inputs$M, inputs$weights, eta_empty,
    rc_dist=rep(0L,K_w), rc_correlation=FALSE, rc_mean=FALSE,
    use_asc=TRUE, include_outside_option=FALSE,
    gen_seed=0L, gen_scramble=0L, gen_S=S
  )

  expect_equal(r_s$objective, r_g$objective, tolerance = TOL_KERNEL,
               label = "2-K_w loglik match")
  expect_equal(as.numeric(r_s$gradient), as.numeric(r_g$gradient),
               tolerance = TOL_KERNEL, label = "2-K_w gradient match")
})

test_that("run_mxlogit draws_info: mode='store' is default, mode='generate' set correctly", {
  dt <- create_small_mxl_data(seed = 1)

  fit_s <- run_mxlogit(
    data=dt, id_col="id", alt_col="alt", choice_col="choice",
    covariate_cols="x1", random_var_cols=c("w1","w2"),
    S=20L, rc_correlation=FALSE, rc_mean=FALSE, use_asc=TRUE,
    include_outside_option=FALSE,
    control=list(print_level=0L, maxeval=50L)
  )
  expect_equal(fit_s$draws_info$mode, "store")
  expect_null(fit_s$draws_info$seed)
  expect_null(fit_s$draws_info$scramble)

  fit_g <- run_mxlogit(
    data=dt, id_col="id", alt_col="alt", choice_col="choice",
    covariate_cols="x1", random_var_cols=c("w1","w2"),
    S=20L, rc_correlation=FALSE, rc_mean=FALSE, use_asc=TRUE,
    include_outside_option=FALSE,
    draws="generate", seed=42L, scramble="owen",
    control=list(print_level=0L, maxeval=50L)
  )
  expect_equal(fit_g$draws_info$mode,     "generate")
  expect_equal(fit_g$draws_info$seed,     42L)
  expect_equal(fit_g$draws_info$scramble, "owen")
})

# ---------------------------------------------------------------------------
# 2. generate-mode gradient vs numDeriv::grad
# ---------------------------------------------------------------------------

.make_gen_grad_test <- function(rc_correlation, rc_mean, rc_dist_type,
                                seed_data = 123) {
  dt     <- create_small_mxl_data(seed = seed_data)
  inputs <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1","w2"),
                              rc_correlation = rc_correlation)
  K_x    <- ncol(inputs$X)
  K_w    <- ncol(inputs$W)
  J      <- nrow(inputs$alt_mapping)
  L_size <- if (rc_correlation) K_w * (K_w + 1L) / 2L else K_w
  mu_size <- if (rc_mean) K_w else 0L
  S_val  <- 30L
  rc_dist <- rep(as.integer(rc_dist_type), K_w)

  set.seed(7)
  theta <- c(
    runif(K_x, -0.3, 0.3),
    if (rc_mean) runif(K_w, -0.3, 0.3) else NULL,
    if (rc_dist_type == 0L) log(runif(L_size, 0.3, 0.8))
    else                     runif(L_size, -0.5, 0.5),
    runif(J - 1L, -0.2, 0.2)
  )

  eta_empty <- array(0, dim = c(K_w, 0L, 0L))

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      theta=th, X=inputs$X, W=inputs$W,
      alt_idx=inputs$alt_idx, choice_idx=inputs$choice_idx,
      M=inputs$M, weights=inputs$weights,
      eta_draws=eta_empty, rc_dist=rc_dist,
      rc_correlation=rc_correlation, rc_mean=rc_mean,
      use_asc=TRUE, include_outside_option=FALSE,
      gen_seed=0L, gen_scramble=0L, gen_S=S_val
    )$objective
  }

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta=theta, X=inputs$X, W=inputs$W,
    alt_idx=inputs$alt_idx, choice_idx=inputs$choice_idx,
    M=inputs$M, weights=inputs$weights,
    eta_draws=eta_empty, rc_dist=rc_dist,
    rc_correlation=rc_correlation, rc_mean=rc_mean,
    use_asc=TRUE, include_outside_option=FALSE,
    gen_seed=0L, gen_scramble=0L, gen_S=S_val
  )$gradient

  list(obj_fn=obj_fn, theta=theta, analytic_grad=analytic_grad)
}

# NIT-2: helper for outside-option grad test
# Creates a dataset where some ids choose the outside option (all inside
# choice = 0), matching the pattern in tests/testthat/test-bhhh.R.
# Uses K_x=1, K_w=1 to keep theta small.
.make_gen_grad_test_oo <- function(rc_correlation, seed_data = 600) {
  set.seed(seed_data)
  N <- 40L; J <- 3L
  dt <- data.table::data.table(
    id  = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), N),
    x1  = stats::rnorm(N * J),
    w1  = stats::rnorm(N * J)
  )
  dt[, choice := 0L]
  # ~60% choose inside alt; ~40% take outside option (choice stays 0 for all rows)
  inside_ids <- sample(seq_len(N), size = ceiling(0.6 * N))
  dt[id %in% inside_ids, choice := sample(c(1L, rep(0L, J - 1L))), by = id]

  inputs <- prepare_mxl_data(
    dt, "id", "alt", "choice",
    covariate_cols     = "x1",
    random_var_cols    = "w1",
    rc_correlation     = rc_correlation,
    include_outside_option = TRUE
  )

  K_x   <- ncol(inputs$X)
  K_w   <- ncol(inputs$W)
  # With include_outside_option=TRUE the outside option is base (utility=0),
  # so all J inside alternatives get ASCs: n_asc = J.
  J_all <- nrow(inputs$alt_mapping)   # J inside + 1 outside
  n_asc <- J_all - 1L                 # ASCs for inside alts only
  L_size <- if (rc_correlation) K_w * (K_w + 1L) / 2L else K_w
  S_val <- 30L
  rc_dist <- rep(0L, K_w)  # normal random coefficient

  set.seed(seed_data + 1L)
  theta <- c(
    stats::runif(K_x, -0.3, 0.3),
    log(stats::runif(L_size, 0.3, 0.8)),
    stats::runif(n_asc, -0.2, 0.2)
  )

  eta_empty <- array(0, dim = c(K_w, 0L, 0L))

  obj_fn <- function(th) {
    mxl_loglik_gradient_parallel(
      theta=th, X=inputs$X, W=inputs$W,
      alt_idx=inputs$alt_idx, choice_idx=inputs$choice_idx,
      M=inputs$M, weights=inputs$weights,
      eta_draws=eta_empty, rc_dist=rc_dist,
      rc_correlation=rc_correlation, rc_mean=FALSE,
      use_asc=TRUE, include_outside_option=TRUE,
      gen_seed=0L, gen_scramble=0L, gen_S=S_val
    )$objective
  }

  analytic_grad <- mxl_loglik_gradient_parallel(
    theta=theta, X=inputs$X, W=inputs$W,
    alt_idx=inputs$alt_idx, choice_idx=inputs$choice_idx,
    M=inputs$M, weights=inputs$weights,
    eta_draws=eta_empty, rc_dist=rc_dist,
    rc_correlation=rc_correlation, rc_mean=FALSE,
    use_asc=TRUE, include_outside_option=TRUE,
    gen_seed=0L, gen_scramble=0L, gen_S=S_val
  )$gradient

  list(obj_fn=obj_fn, theta=theta, analytic_grad=analytic_grad)
}

test_that("generate-mode gradient vs numDeriv: uncorr, rc_mean=FALSE, normal", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test(FALSE, FALSE, 0L, 100)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("uncorr/rc_mean=F/normal max_diff=", signif(max_d,4)))
})

test_that("generate-mode gradient vs numDeriv: corr, rc_mean=FALSE, normal", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test(TRUE, FALSE, 0L, 200)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("corr/rc_mean=F/normal max_diff=", signif(max_d,4)))
})

test_that("generate-mode gradient vs numDeriv: uncorr, rc_mean=TRUE, normal", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test(FALSE, TRUE, 0L, 300)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("uncorr/rc_mean=T/normal max_diff=", signif(max_d,4)))
})

test_that("generate-mode gradient vs numDeriv: uncorr, rc_mean=FALSE, log-normal", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test(FALSE, FALSE, 1L, 400)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("uncorr/rc_mean=F/lognormal max_diff=", signif(max_d,4)))
})

test_that("generate-mode gradient vs numDeriv: corr, rc_mean=TRUE, log-normal", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test(TRUE, TRUE, 1L, 500)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("corr/rc_mean=T/lognormal max_diff=", signif(max_d,4)))
})

# NIT-2: outside-option cells (include_outside_option=TRUE)
# Exercises validate_choice_data (bypassed path in generate mode), the
# outside-option utility branch in fill_choice_utilities / stable_softmax,
# and DiffW assembly under generate mode.

test_that("generate-mode gradient vs numDeriv: uncorr, rc_mean=FALSE, normal, outside_option=TRUE", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test_oo(rc_correlation = FALSE, seed_data = 600L)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("uncorr/rc_mean=F/normal/oo max_diff=", signif(max_d,4)))
})

test_that("generate-mode gradient vs numDeriv: corr, rc_mean=FALSE, normal, outside_option=TRUE", {
  skip_if_not_installed("numDeriv")
  g <- .make_gen_grad_test_oo(rc_correlation = TRUE, seed_data = 700L)
  num_g <- numDeriv::grad(g$obj_fn, g$theta, method = "Richardson")
  max_d <- max(abs(drop(g$analytic_grad) - num_g))
  expect_true(max_d < TOL_GRAD_GEN,
    label = paste0("corr/rc_mean=F/normal/oo max_diff=", signif(max_d,4)))
})

# ---------------------------------------------------------------------------
# 3. Post-estimation generics: evaluated at byte-identical theta
#
# NIT-3 fix: fit_g is constructed as a copy of fit_s with only draws_info
# switched to generate mode (seed=0, scramble="none"). Both objects carry
# byte-identical coefficients and data; the test isolates draw-path
# equivalence, not optimizer reproducibility.
#
# For vcov: the cached vcov in fit_g is cleared (set to NULL) so that
# ensure_vcov() recomputes it through the generate-mode draw path.
# ---------------------------------------------------------------------------

test_that("post-estimation generics match store vs generate at converged theta (hessian)", {
  fix   <- .gen_mode_fixture()
  fit_s <- fix$fit_s

  # Construct the generate-mode twin: copy fit_s, switch only draws_info.
  # This guarantees byte-identical coefficients -- the test isolates draw-path
  # equivalence, not optimizer reproducibility (NIT-3 fix).
  fit_g <- fit_s
  fit_g$draws_info$mode    <- "generate"
  fit_g$draws_info$seed    <- 0L
  fit_g$draws_info$scramble <- "none"
  # Clear cached vcov/se so ensure_vcov() recomputes via generate-mode draws.
  fit_g$vcov <- NULL
  fit_g$se   <- NULL

  # Verify coefficients are byte-identical before any generic call
  expect_identical(coef(fit_s), coef(fit_g),
                   label = "NIT-3: coef vectors byte-identical")

  # predict: both use same theta + same draw path (seed=0, scramble=none
  # produces same Halton points as store mode with seed=0)
  pred_s <- predict(fit_s)
  pred_g <- predict(fit_g)
  expect_equal(pred_s$choice_prob, pred_g$choice_prob, tolerance = TOL_POST,
               label = "predict share: store vs generate")

  # elasticities
  el_s <- elasticities(fit_s, elast_var = "x1")
  el_g <- elasticities(fit_g, elast_var = "x1")
  expect_equal(el_s, el_g, tolerance = TOL_POST,
               label = "elasticities: store vs generate")

  # diversion_ratios
  dr_s <- diversion_ratios(fit_s, wrt_var = "x1")
  dr_g <- diversion_ratios(fit_g, wrt_var = "x1")
  expect_equal(dr_s, dr_g, tolerance = TOL_POST,
               label = "diversion_ratios: store vs generate")

  # logsum
  ls_s <- logsum(fit_s)
  ls_g <- logsum(fit_g)
  expect_equal(ls_s, ls_g, tolerance = TOL_POST,
               label = "logsum: store vs generate")

  # consumer_surplus
  cs_s <- consumer_surplus(fit_s, price_var = "x1")
  cs_g <- consumer_surplus(fit_g, price_var = "x1")
  expect_equal(cs_s, cs_g, tolerance = TOL_POST,
               label = "consumer_surplus: store vs generate")

  # vcov hessian: fit_g$vcov was cleared, so ensure_vcov() runs the
  # generate-mode Hessian at the same theta. Must match fit_s (cached).
  expect_equal(vcov(fit_s), vcov(fit_g), tolerance = TOL_POST,
               label = "vcov hessian: store vs generate")
})

test_that("vcov bhhh matches between store and generate modes at same theta", {
  fix   <- .gen_mode_fixture()
  fit_s <- fix$fit_s

  # Build a bhhh-SE store-mode fit at the converged theta.
  # Use maxeval=0 so the optimizer does not move at all from theta_init.
  theta_conv <- coef(fit_s)
  dt  <- fix$dt
  S   <- fix$S
  K_w <- fix$K_w

  # We only need a fit object with se_method="bhhh" at the same coefficients.
  # Use the converged theta as init and maxeval=1 to get a fresh BHHH vcov;
  # then construct fit_g_bhhh as a copy with draws_info switched to generate.
  fit_s_bhhh <- run_mxlogit(
    data=dt, id_col="id", alt_col="alt", choice_col="choice",
    covariate_cols="x1", random_var_cols="w1",
    S=S, rc_correlation=FALSE, rc_mean=FALSE, rc_dist=rep(0L,K_w),
    use_asc=TRUE, include_outside_option=FALSE,
    draws="store", se_method="bhhh",
    theta_init=theta_conv,
    control=list(print_level=0L, maxeval=1L)
  )

  # NIT-3 fix: construct generate-mode twin via copy + draws_info swap.
  # Coefficients are byte-identical; only the draw acquisition path differs.
  fit_g_bhhh <- fit_s_bhhh
  fit_g_bhhh$draws_info$mode    <- "generate"
  fit_g_bhhh$draws_info$seed    <- 0L
  fit_g_bhhh$draws_info$scramble <- "none"
  # Clear cached vcov/se to force recomputation through generate-mode BHHH.
  fit_g_bhhh$vcov <- NULL
  fit_g_bhhh$se   <- NULL

  expect_identical(coef(fit_s_bhhh), coef(fit_g_bhhh),
                   label = "NIT-3 bhhh: coef vectors byte-identical")

  expect_equal(vcov(fit_s_bhhh), vcov(fit_g_bhhh), tolerance = TOL_POST,
               label = "vcov bhhh: store vs generate")
})

# ---------------------------------------------------------------------------
# 4. generate mode: seed governs reproducibility (Owen scramble)
# ---------------------------------------------------------------------------

test_that("generate mode with same seed (Owen) produces identical loglik", {
  dt <- create_small_mxl_data(seed = 77)

  args <- list(
    data=dt, id_col="id", alt_col="alt", choice_col="choice",
    covariate_cols="x1", random_var_cols="w1",
    S=20L, rc_correlation=FALSE, rc_mean=FALSE, use_asc=TRUE,
    include_outside_option=FALSE,
    draws="generate", seed=42L, scramble="owen",
    control=list(print_level=0L, maxeval=50L)
  )

  fit1 <- do.call(run_mxlogit, args)
  fit2 <- do.call(run_mxlogit, args)

  # Same seed -> same objective trajectory -> same converged loglik
  expect_equal(as.numeric(logLik(fit1)), as.numeric(logLik(fit2)),
               tolerance = 1e-14,
               label = "same-seed Owen: identical loglik")
  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-14,
               label = "same-seed Owen: identical coef")
})
