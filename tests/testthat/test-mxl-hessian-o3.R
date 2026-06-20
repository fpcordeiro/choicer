# tests/testthat/test-mxl-hessian-o3.R
#
# Durable regression tests for mxl_hessian_parallel (O3 BLAS-3 restructure).
# This file is the permanent test-suite guard for the O3 Hessian path; it
# was authored by the tester (independently of implementation) from
# runs/2026-06-20-mxl-hessian-blas3-o3/test-spec.md.
#
# Test groups
# -----------
# A. numDeriv oracle — analytical vs finite-difference Hessian of the MXL
#    log-likelihood, across a representative config slice:
#      A1: normal RC, uncorrelated, no mean (baseline)
#      A2: normal RC, correlated, with mean
#      A3: log-normal RC, uncorrelated, no mean
#      A4: outside option (include_outside_option = TRUE)
#      A5: use_asc = FALSE (Jd == 0 guard)
#    Tolerance: 1e-4 (TOL_HESS from setup.R — the standard package FD tolerance).
#    Justification: Richardson FD of the negated simulated log-likelihood can
#    have O(h^2) errors; for Hessian elements of magnitude ~10-100, absolute
#    error reaches ~1e-4 to 1e-5.  1e-4 is safely above FD noise without being
#    loose enough to hide algorithmic errors (the prior ad-hoc validation showed
#    worst-case analytical-vs-oracle deviation of 4.6e-13, 9 orders of magnitude
#    smaller than this gate).
#    Skip: skip_on_cran() + skip_if_not_installed("numDeriv")
#
# B. Property invariants — always-on, tight tolerances, no CRAN skip:
#      B1: symmetry (tol 1e-10) across 6 representative configs
#      B2: all-finite across 5 configs
#      B3: negation-sign (H[1,1] > 0) for normal and log-normal rc_dist
#      B4: weight-scaling linearity — H(2w) == 2*H(w) to within 1e-10
#
# C. wesml_vcov() coverage — closes reviewer Note 1:
#    Explicit end-to-end check that wesml_vcov() sandwich SEs are finite,
#    symmetric PSD, and consistent with fit$se on a WESML-weighted MXL fit.
#    Gated skip_on_cran() (requires fitting).
#
# OpenMP threads are capped by setup.R (2 threads); no new omp calls here.

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

# FD oracle tolerance — matches package-wide TOL_HESS in setup.R; see §A above
TOL_O3_FD    <- 1e-4

# Property-invariant tolerances
TOL_O3_SYMM  <- 1e-10  # symmetry max|H - H'|
TOL_O3_SCALE <- 1e-10  # weight-scaling linearity max|H(2w) - 2*H(w)|

# ---------------------------------------------------------------------------
# Shared fixture helper
# ---------------------------------------------------------------------------

# Compute the number of Cholesky parameters for K_w random coefficients.
# When rc_correlation=TRUE: K_w*(K_w+1)/2 (lower-triangular matrix).
# When rc_correlation=FALSE: K_w (diagonal only).
# Note: uses explicit parentheses to avoid R's operator-precedence ambiguity
# between * and %/% (both have the same precedence; explicit parens are safer).
.o3_K_L <- function(K_w, rc_correlation) {
  if (rc_correlation) as.integer((K_w * (K_w + 1L)) / 2L) else as.integer(K_w)
}

# Build a small MXL dataset and all C++-level inputs for mxl_hessian_parallel.
# Returns a list with inputs, eta_draws, theta, and configuration flags.
.o3_mxl_fixture <- function(
    N, J, K_w,
    rc_dist, rc_correlation, rc_mean,
    use_asc, include_outside_option,
    S = 30L,
    weights = NULL,
    seed_data = 42L,
    seed_theta = 7L
) {
  set.seed(seed_data)

  w_cols <- paste0("w", seq_len(K_w))

  if (include_outside_option) {
    J_inside <- J
    n_cols   <- N * (J_inside + 1L)
    dt <- data.table::data.table(
      id  = rep(seq_len(N), each = J_inside + 1L),
      alt = rep(0:J_inside, N),
      x1  = c(rep(0, N), rnorm(N * J_inside)),
      w1  = c(rep(0, N), rnorm(N * J_inside)),
      w2  = c(rep(0, N), runif(N * J_inside, -1, 1))
    )
    dt[, choice := 0L]
    dt[, choice := {
      ch <- sample(0:J_inside, 1L)
      as.integer(alt == ch)
    }, by = id]
    dt[, n_ch := sum(choice), by = id]
    dt[n_ch == 0L, choice := as.integer(alt == 0L)]
    dt[, n_ch := NULL]
    inputs <- prepare_mxl_data(
      dt, "id", "alt", "choice",
      covariate_cols         = "x1",
      random_var_cols        = w_cols,
      rc_correlation         = rc_correlation,
      include_outside_option = TRUE
    )
  } else {
    n_rows <- N * J
    dt <- data.table::data.table(
      id  = rep(seq_len(N), each = J),
      alt = rep(seq_len(J), N),
      x1  = rnorm(n_rows),
      w1  = rnorm(n_rows),
      w2  = runif(n_rows, -1, 1)
    )
    dt[, choice := 0L]
    dt[, choice := sample(c(1L, rep(0L, J - 1L))), by = id]
    inputs <- prepare_mxl_data(
      dt, "id", "alt", "choice",
      covariate_cols         = "x1",
      random_var_cols        = w_cols,
      rc_correlation         = rc_correlation,
      include_outside_option = FALSE
    )
  }

  if (!is.null(weights)) inputs$weights <- weights

  eta_draws <- get_halton_normals(S, N, K_w)

  # Build theta: beta(K_x) | [mu(K_w) if rc_mean] | L(K_L) | [asc(J-1 or J)]
  set.seed(seed_theta)
  K_x_val <- ncol(inputs$X)
  K_w_val <- ncol(inputs$W)
  K_L_val <- .o3_K_L(K_w_val, rc_correlation)
  n_asc   <- if (!use_asc) 0L else (nrow(inputs$alt_mapping) - 1L)
  n_params <- K_x_val +
              (if (rc_mean) K_w_val else 0L) +
              K_L_val +
              n_asc
  theta <- runif(n_params, -0.3, 0.3)

  list(
    inputs                 = inputs,
    eta_draws              = eta_draws,
    theta                  = theta,
    rc_dist                = rc_dist,
    rc_correlation         = rc_correlation,
    rc_mean                = rc_mean,
    use_asc                = use_asc,
    include_outside_option = include_outside_option
  )
}

# Call mxl_hessian_parallel from a fixture.
.o3_call_H <- function(fx, weights_override = NULL) {
  wts <- if (!is.null(weights_override)) weights_override else fx$inputs$weights
  mxl_hessian_parallel(
    theta                  = fx$theta,
    X                      = fx$inputs$X,
    W                      = fx$inputs$W,
    alt_idx                = fx$inputs$alt_idx,
    choice_idx             = fx$inputs$choice_idx,
    M                      = fx$inputs$M,
    weights                = wts,
    eta_draws              = fx$eta_draws,
    rc_dist                = fx$rc_dist,
    rc_correlation         = fx$rc_correlation,
    rc_mean                = fx$rc_mean,
    use_asc                = fx$use_asc,
    include_outside_option = fx$include_outside_option
  )
}

# Negated log-likelihood objective for numDeriv.
.o3_obj_fn <- function(theta, fx) {
  mxl_loglik_gradient_parallel(
    theta                  = theta,
    X                      = fx$inputs$X,
    W                      = fx$inputs$W,
    alt_idx                = fx$inputs$alt_idx,
    choice_idx             = fx$inputs$choice_idx,
    M                      = fx$inputs$M,
    weights                = fx$inputs$weights,
    eta_draws              = fx$eta_draws,
    rc_dist                = fx$rc_dist,
    rc_correlation         = fx$rc_correlation,
    rc_mean                = fx$rc_mean,
    use_asc                = fx$use_asc,
    include_outside_option = fx$include_outside_option
  )$objective
}

# ===========================================================================
# Section A — numDeriv oracle (skip_on_cran + skip_if_not_installed)
# Tolerance: TOL_O3_FD = 1e-4, matching the package-wide FD tolerance TOL_HESS.
# ===========================================================================

test_that("O3 A1: normal-uncorr-no-mean — analytical Hessian matches numDeriv (tol=1e-4)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  fx <- .o3_mxl_fixture(
    N = 60L, J = 4L, K_w = 2L,
    rc_dist        = c(0L, 0L),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 30L, seed_data = 42L, seed_theta = 7L
  )

  H_anal <- .o3_call_H(fx)
  H_num  <- numDeriv::hessian(function(th) .o3_obj_fn(th, fx), fx$theta,
                               method = "Richardson")
  diff <- max(abs(H_anal - H_num))

  expect_true(all(is.finite(H_anal)), label = "A1: Hessian all-finite")
  expect_lt(diff, TOL_O3_FD,
    label = paste0("A1 max|H_anal-H_num| = ", format(diff, scientific = TRUE),
                   " (tol=", TOL_O3_FD, ")"))
})

test_that("O3 A2: normal-corr-with-mean — analytical Hessian matches numDeriv (tol=1e-4)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  fx <- .o3_mxl_fixture(
    N = 50L, J = 3L, K_w = 2L,
    rc_dist        = c(0L, 0L),
    rc_correlation = TRUE, rc_mean = TRUE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 30L, seed_data = 43L, seed_theta = 8L
  )

  H_anal <- .o3_call_H(fx)
  H_num  <- numDeriv::hessian(function(th) .o3_obj_fn(th, fx), fx$theta,
                               method = "Richardson")
  diff <- max(abs(H_anal - H_num))

  expect_true(all(is.finite(H_anal)), label = "A2: Hessian all-finite")
  expect_lt(diff, TOL_O3_FD,
    label = paste0("A2 max|H_anal-H_num| = ", format(diff, scientific = TRUE)))
})

test_that("O3 A3: log-normal-uncorr-no-mean — analytical Hessian matches numDeriv (tol=1e-4)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  fx <- .o3_mxl_fixture(
    N = 50L, J = 3L, K_w = 1L,
    rc_dist        = c(1L),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 30L, seed_data = 44L, seed_theta = 9L
  )

  H_anal <- .o3_call_H(fx)
  H_num  <- numDeriv::hessian(function(th) .o3_obj_fn(th, fx), fx$theta,
                               method = "Richardson")
  diff <- max(abs(H_anal - H_num))

  expect_true(all(is.finite(H_anal)), label = "A3: Hessian all-finite")
  expect_lt(diff, TOL_O3_FD,
    label = paste0("A3 max|H_anal-H_num| = ", format(diff, scientific = TRUE)))
})

test_that("O3 A4: outside-option — analytical Hessian matches numDeriv (tol=1e-4)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  fx <- .o3_mxl_fixture(
    N = 50L, J = 3L, K_w = 2L,
    rc_dist        = c(0L, 0L),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = TRUE,
    S = 30L, seed_data = 45L, seed_theta = 10L
  )

  H_anal <- .o3_call_H(fx)
  H_num  <- numDeriv::hessian(function(th) .o3_obj_fn(th, fx), fx$theta,
                               method = "Richardson")
  diff <- max(abs(H_anal - H_num))

  expect_true(all(is.finite(H_anal)), label = "A4: Hessian all-finite")
  expect_lt(diff, TOL_O3_FD,
    label = paste0("A4 max|H_anal-H_num| = ", format(diff, scientific = TRUE)))
})

test_that("O3 A5: use_asc=FALSE (Jd=0) — analytical Hessian matches numDeriv (tol=1e-4)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  fx <- .o3_mxl_fixture(
    N = 50L, J = 4L, K_w = 2L,
    rc_dist        = c(0L, 0L),
    rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = FALSE, include_outside_option = FALSE,
    S = 30L, seed_data = 46L, seed_theta = 11L
  )

  H_anal <- .o3_call_H(fx)
  H_num  <- numDeriv::hessian(function(th) .o3_obj_fn(th, fx), fx$theta,
                               method = "Richardson")
  diff <- max(abs(H_anal - H_num))

  expect_true(all(is.finite(H_anal)), label = "A5: Hessian all-finite")
  expect_lt(diff, TOL_O3_FD,
    label = paste0("A5 max|H_anal-H_num| = ", format(diff, scientific = TRUE)))
})

# ===========================================================================
# Section B — Property invariants (always-on, tight tolerances, no CRAN skip)
# ===========================================================================

# -- B1: Symmetry (tol 1e-10) across 6 configs --------------------------------

test_that("O3 B1a: symmetry — normal-uncorr-no-mean", {
  fx <- .o3_mxl_fixture(
    N = 40L, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 25L, seed_data = 101L, seed_theta = 21L
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1a: max|H-H'| for normal-uncorr-no-mean")
})

test_that("O3 B1b: symmetry — normal-corr-with-mean", {
  fx <- .o3_mxl_fixture(
    N = 40L, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = TRUE, rc_mean = TRUE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 25L, seed_data = 102L, seed_theta = 22L
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1b: max|H-H'| for normal-corr-with-mean")
})

test_that("O3 B1c: symmetry — log-normal-uncorr-no-mean", {
  fx <- .o3_mxl_fixture(
    N = 40L, J = 3L, K_w = 1L,
    rc_dist = c(1L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 25L, seed_data = 103L, seed_theta = 23L
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1c: max|H-H'| for log-normal-uncorr-no-mean")
})

test_that("O3 B1d: symmetry — outside-option", {
  fx <- .o3_mxl_fixture(
    N = 40L, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = TRUE,
    S = 25L, seed_data = 104L, seed_theta = 24L
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1d: max|H-H'| for outside-option")
})

test_that("O3 B1e: symmetry — use_asc=FALSE (Jd=0)", {
  fx <- .o3_mxl_fixture(
    N = 40L, J = 4L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = FALSE, include_outside_option = FALSE,
    S = 25L, seed_data = 105L, seed_theta = 25L
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1e: max|H-H'| for use_asc=FALSE")
})

test_that("O3 B1f: symmetry — WESML non-uniform weights", {
  N <- 40L
  set.seed(106L)
  w_i <- sample(c(0.5, 1.0, 2.0), N, replace = TRUE)
  fx <- .o3_mxl_fixture(
    N = N, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 25L, seed_data = 106L, seed_theta = 26L,
    weights = w_i
  )
  H <- .o3_call_H(fx)
  expect_lt(max(abs(H - t(H))), TOL_O3_SYMM,
    label = "B1f: max|H-H'| for WESML weights")
})

# -- B2: All-finite across 5 configs -------------------------------------------

test_that("O3 B2: all-finite — 5 representative configs", {
  configs <- list(
    list(K_w=2L, rc_dist=c(0L,0L), rc_corr=FALSE, rc_mean=FALSE,
         use_asc=TRUE,  oo=FALSE, J=4L, sd=201L, st=31L, label="normal-uncorr-no-mean"),
    list(K_w=2L, rc_dist=c(0L,0L), rc_corr=TRUE,  rc_mean=TRUE,
         use_asc=TRUE,  oo=FALSE, J=3L, sd=202L, st=32L, label="normal-corr-with-mean"),
    list(K_w=1L, rc_dist=c(1L),   rc_corr=FALSE, rc_mean=FALSE,
         use_asc=TRUE,  oo=FALSE, J=3L, sd=203L, st=33L, label="log-normal-uncorr"),
    list(K_w=2L, rc_dist=c(0L,0L), rc_corr=FALSE, rc_mean=FALSE,
         use_asc=TRUE,  oo=TRUE,  J=3L, sd=204L, st=34L, label="outside-option"),
    list(K_w=2L, rc_dist=c(0L,0L), rc_corr=FALSE, rc_mean=FALSE,
         use_asc=FALSE, oo=FALSE, J=4L, sd=205L, st=35L, label="use_asc=FALSE")
  )
  for (cfg in configs) {
    fx <- .o3_mxl_fixture(
      N = 40L, J = cfg$J, K_w = cfg$K_w,
      rc_dist        = cfg$rc_dist,
      rc_correlation = cfg$rc_corr,
      rc_mean        = cfg$rc_mean,
      use_asc        = cfg$use_asc,
      include_outside_option = cfg$oo,
      S = 25L, seed_data = cfg$sd, seed_theta = cfg$st
    )
    H <- .o3_call_H(fx)
    expect_true(all(is.finite(H)),
      label = paste0("B2 all-finite for config '", cfg$label, "'"))
  }
})

# -- B3: Negation-sign — H[1,1] > 0 at generic theta --------------------------
# mxl_hessian_parallel returns Hess of the NEGATED log-likelihood;
# H[1,1] > 0 at any theta confirms the sign convention is correct.

test_that("O3 B3: negation-sign — H[1,1] > 0 for normal and log-normal rc_dist", {
  fx_norm <- .o3_mxl_fixture(
    N = 50L, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 30L, seed_data = 301L, seed_theta = 41L
  )
  H_norm <- .o3_call_H(fx_norm)
  expect_gt(H_norm[1L, 1L], 0,
    label = "B3: H[1,1] > 0 for normal RC (negated log-likelihood convention)")

  fx_ln <- .o3_mxl_fixture(
    N = 50L, J = 3L, K_w = 1L,
    rc_dist = c(1L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 30L, seed_data = 302L, seed_theta = 42L
  )
  H_ln <- .o3_call_H(fx_ln)
  expect_gt(H_ln[1L, 1L], 0,
    label = "B3: H[1,1] > 0 for log-normal RC (negated log-likelihood convention)")
})

# -- B4: Weight-scaling linearity — H(2w) = 2*H(w) to within 1e-10 -----------
# H = sum_i w_i H_i is linear in weights; H_i does not depend on w_i.
# Doubling all weights must therefore exactly double the Hessian.

test_that("O3 B4: weight-scaling linearity — H(2w) = 2*H(w) within 1e-10", {
  N <- 40L
  set.seed(401L)
  w_base <- runif(N, 0.5, 2.0)
  fx <- .o3_mxl_fixture(
    N = N, J = 3L, K_w = 2L,
    rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE,
    use_asc = TRUE, include_outside_option = FALSE,
    S = 25L, seed_data = 401L, seed_theta = 51L,
    weights = w_base
  )
  H_1w <- .o3_call_H(fx, weights_override = w_base)
  H_2w <- .o3_call_H(fx, weights_override = 2.0 * w_base)
  dev  <- max(abs(H_2w - 2.0 * H_1w))
  expect_lt(dev, TOL_O3_SCALE,
    label = paste0("B4: max|H(2w)-2*H(w)| = ", format(dev, scientific = TRUE),
                   " (tol=", TOL_O3_SCALE, ")"))
})

# ===========================================================================
# Section C — wesml_vcov() coverage (closes reviewer Note 1)
# Explicit end-to-end check of the sandwich/WESML SE path on a weighted fit.
# ===========================================================================

test_that("O3 C1: wesml_vcov() returns finite, symmetric PSD sandwich vcov on WESML-weighted fit", {
  skip_on_cran()

  set.seed(501L)
  N <- 150L; J <- 3L
  sim <- simulate_mxl_data(
    N   = N, J = J,
    beta  = 0.8,
    delta = c(0.3, -0.3),
    Sigma = matrix(0.4, 1L, 1L),
    seed  = 501L,
    outside_option  = FALSE,
    vary_choice_set = FALSE
  )
  dt <- sim$data

  # WESML weights: population shares (0.40, 0.35, 0.25) vs empirical shares
  pop_shares  <- c(0.40, 0.35, 0.25)
  chosen_alts <- dt[choice == 1L, alt]
  emp_tab     <- table(factor(chosen_alts, levels = seq_len(J)))
  emp_shares  <- as.numeric(emp_tab) / sum(emp_tab)
  w_i         <- pop_shares[chosen_alts] / emp_shares[chosen_alts]

  fit_sw <- suppressWarnings(suppressMessages(
    run_mxlogit(dt, "id", "alt", "choice", "x1", "w1",
                S = 40L, weights = w_i,
                se_method = "sandwich",
                keep_data = TRUE)
  ))

  # (C1a) sandwich SEs are all-finite and positive
  expect_true(all(is.finite(fit_sw$se)),
    label = "C1a: sandwich SEs all-finite on WESML-weighted fit")
  expect_true(all(fit_sw$se > 0),
    label = "C1a: sandwich SEs all-positive on WESML-weighted fit")

  # (C1b) wesml_vcov() vcov is finite and symmetric PSD
  vcov_sw <- suppressWarnings(wesml_vcov(fit_sw, "vcov"))
  expect_true(all(is.finite(vcov_sw)),
    label = "C1b: wesml_vcov() vcov all-finite")
  expect_lt(max(abs(vcov_sw - t(vcov_sw))), 1e-10,
    label = "C1b: wesml_vcov() vcov symmetric to 1e-10")
  ev <- eigen(vcov_sw, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-6,
    label = paste0("C1b: wesml_vcov() PSD, min eigenvalue = ",
                   format(min(ev), scientific = TRUE)))

  # (C1c) wesml_vcov("se") is consistent with sqrt(diag(fit$vcov)) to 1e-6 rel tol
  se_from_vcov <- wesml_vcov(fit_sw, "se")
  reldiff <- max(abs(se_from_vcov - fit_sw$se) /
                   pmax(abs(fit_sw$se), 1e-10))
  expect_lt(reldiff, 1e-6,
    label = paste0("C1c: wesml_vcov('se') vs fit$se max rel diff = ",
                   format(reldiff, scientific = TRUE)))
})

