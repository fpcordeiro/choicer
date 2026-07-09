# Phase-1 HB diagnostics / parallel-chains regression tests: ess()/mcse()/
# rhat(rank=) correctness (incl. the fixed Geyer lag-0 bug), chain retention &
# back-compat, Lever G determinism, the consolidated summary() table,
# traceplot(), and the corrected keep_beta_i="draws" memory guard.
#
# CRAN-fast: small fixtures (N<=60, R<=300). Kernel-level bitwise-invariance
# and end-to-end recovery live in test-hmnl-gibbs.R / test-hmnp-gibbs.R /
# test-hmnl-recovery.R / test-hmnp-recovery.R; this file is scoped to the
# NEW Phase-1 surface (ess, mcse, rhat(rank=), object$chains, mcmc$cores,
# the consolidated diagnostic table, traceplot(), and the memory guard).

# --- Fixtures ----------------------------------------------------------------

.ar1_series <- function(phi, n = 5000, seed = 1) {
  set.seed(seed)
  eps <- stats::rnorm(n)
  x <- numeric(n)
  x[1] <- eps[1]
  for (i in 2:n) x[i] <- phi * x[i - 1] + eps[i]
  matrix(x, ncol = 1, dimnames = list(NULL, "x"))
}

.hmnl_fx_phase1 <- function(N = 60, T = 2, J = 4, seed = 1) {
  sim <- simulate_hmnl_data(N = N, T = T, J = J, seed = seed,
                            theta = c(0.5, -0.4), sigma_d = 0.5)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  list(sim = sim, d = d)
}

.hmnp_fx_phase1 <- function(N = 60, T = 2, J = 4, seed = 1) {
  sim <- simulate_hmnp_data(N = N, T = T, J = J, seed = seed,
                            theta = c(0.5, -0.4), sigma_d = 0.5)
  d <- prepare_hmnp_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1")
  list(sim = sim, d = d)
}

# --- 1. ess(): closed-form AR(1) validation (no coda dependency) ------------

test_that("ess() bulk matches closed-form AR(1) ESS within tolerance", {
  n <- 5000
  for (phi in c(0, 0.4, 0.8)) {
    m <- .ar1_series(phi, n = n, seed = 100 + round(phi * 10))
    e <- ess(m)
    closed_form <- n * (1 - phi) / (1 + phi)
    expect_equal(unname(e[1, "bulk"]), closed_form, tolerance = 0.20)
  }
})

test_that("ess() recovers ~n for iid draws and is far below n for phi=0.95", {
  n <- 5000
  e_iid <- ess(.ar1_series(0, n = n, seed = 11))
  expect_equal(unname(e_iid[1, "bulk"]), n, tolerance = 0.10)

  e_hi <- ess(.ar1_series(0.95, n = n, seed = 12))
  expect_lt(unname(e_hi[1, "bulk"]), n * 0.10)
})

test_that("ess() regression pin: phi=0.4 bulk ESS is clearly below n (Geyer lag-0 fix)", {
  # The pre-fix bug omitted rho_0 = 1 from the Geyer pairing (pairing
  # rho_1+rho_2, rho_3+rho_4, ... instead of rho_0+rho_1, rho_2+rho_3, ...),
  # which silently dropped 1 from tau_hat and (a) returned ess ~= n with NO
  # autocorrelation penalty at low/moderate phi, and (b) roughly doubled
  # ess_bulk at higher phi. Pin that both symptoms are gone.
  n <- 5000
  e <- ess(.ar1_series(0.4, n = n, seed = 21))
  expect_lt(unname(e[1, "bulk"]), n * 0.85)   # old bug: ~= n, no penalty
  expect_gt(unname(e[1, "bulk"]), n * 0.30)   # sanity: not degenerate either
})

test_that("ess() return shape and NA-for-constant-column contract", {
  set.seed(5)
  m <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("a", "b")))
  e <- ess(m)
  expect_true(is.matrix(e))
  expect_equal(colnames(e), c("bulk", "tail"))
  expect_equal(rownames(e), c("a", "b"))
  expect_true(all(is.finite(e)))

  const_col <- cbind(a = rep(1, 100), b = rnorm(100))
  e_const <- ess(const_col)
  expect_true(is.na(e_const["a", "bulk"]))
  expect_true(is.na(e_const["a", "tail"]))
  expect_true(is.finite(e_const["b", "bulk"]))
})

test_that("ess() accepts a multi-chain list with identical dims", {
  set.seed(6)
  ch <- list(matrix(rnorm(400), ncol = 1, dimnames = list(NULL, "a")),
             matrix(rnorm(400), ncol = 1, dimnames = list(NULL, "a")))
  e <- ess(ch)
  expect_true(is.finite(e["a", "bulk"]))
  expect_error(ess(list(matrix(1:8, ncol = 2), matrix(1:6, ncol = 2))),
               "identical dimensions")
  expect_error(ess(matrix(1:3, ncol = 1)), "at least 4 draws")
})

# --- 2. mcse() ----------------------------------------------------------------

test_that("mcse(kind = 'mean') matches sd_pooled / sqrt(ess_bulk)", {
  set.seed(30)
  m <- matrix(rnorm(2000), ncol = 2, dimnames = list(NULL, c("a", "b")))
  e <- ess(m)
  mc <- mcse(m, kind = "mean")
  sd_pooled <- apply(m, 2, stats::sd)
  expect_equal(unname(mc), unname(sd_pooled / sqrt(e[, "bulk"])),
               tolerance = 1e-10)
  expect_true(all(is.finite(mc)) && all(mc > 0))
})

test_that("mcse() NA propagates for a constant column", {
  set.seed(31)
  const_col <- cbind(a = rep(2, 100), b = rnorm(100))
  mc <- mcse(const_col)
  expect_true(is.na(mc["a"]))
  expect_true(is.finite(mc["b"]))

  mc_med <- mcse(const_col, kind = "median")
  expect_true(is.na(mc_med["a"]))
  expect_true(is.finite(mc_med["b"]))
})

# --- 3. rhat(rank =) ----------------------------------------------------------

test_that("rhat() ~1 for white noise, both rank = FALSE and rank = TRUE", {
  set.seed(42)
  draws <- matrix(rnorm(4000), ncol = 1, dimnames = list(NULL, "a"))
  expect_lt(abs(rhat(draws)[["a"]] - 1), 0.05)
  expect_lt(abs(rhat(draws, rank = TRUE)[["a"]] - 1), 0.05)
})

test_that("rhat() >> 1.2 for a random-walk (non-stationary) fixture", {
  set.seed(7)
  drifting <- matrix(cumsum(rnorm(2000)), ncol = 1, dimnames = list(NULL, "a"))
  expect_gt(rhat(drifting)[["a"]], 1.2)
  expect_gt(rhat(drifting, rank = TRUE)[["a"]], 1.2)
})

test_that("rhat(draws) is identical to rhat(draws, rank = FALSE) [legacy pin]", {
  set.seed(9)
  draws <- matrix(rnorm(2000), ncol = 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(rhat(draws), rhat(draws, rank = FALSE))

  # multi-chain form too
  ch <- list(draws, matrix(rnorm(2000), ncol = 2, dimnames = list(NULL, c("a", "b"))))
  expect_identical(rhat(ch), rhat(ch, rank = FALSE))
})

# --- 4. Chain retention / back-compat ----------------------------------------

test_that("chains = 2 HMNL fit exposes object$chains with documented per-chain fields", {
  fx <- .hmnl_fx_phase1()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2)
  )
  expect_length(fit$chains, 2)
  expect_equal(fit$mcmc$chains, 2L)
  for (ch in fit$chains) {
    expect_setequal(names(ch),
                    c("b", "w_vech", "delta", "theta", "sigma_d2", "loglik"))
    expect_equal(nrow(ch$b), fit$mcmc$R_keep)
    expect_equal(colnames(ch$b), colnames(fit$draws$b))
    expect_equal(colnames(ch$delta), colnames(fit$draws$delta))
    expect_equal(colnames(ch$theta), colnames(fit$draws$theta))
    expect_true(all(is.finite(ch$b)))
  }
  # object$beta_i is NOT replicated per chain (chain-1-only, unchanged field)
  expect_true(is.null(fit$beta_i) || is.list(fit$beta_i))
})

test_that("chains = 2 HMNP fit exposes object$chains with documented per-chain fields", {
  fx <- .hmnp_fx_phase1()
  fit <- suppressWarnings(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2)
  )
  expect_length(fit$chains, 2)
  expect_equal(fit$mcmc$chains, 2L)
  for (ch in fit$chains) {
    expect_setequal(names(ch),
                    c("b", "w_vech", "delta", "theta", "sigma_d2", "sigma2"))
    expect_equal(nrow(ch$b), fit$mcmc$R_keep)
    expect_true(all(is.finite(ch$b)))
    expect_true(all(ch$sigma2 > 0))
  }
})

test_that("chains = 1 fixed-seed HMNL: draws are identical to chains[[1]] (bitwise back-compat)", {
  fx <- .hmnl_fx_phase1()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 300, burn = 100, seed = 5))
  )
  expect_length(fit$chains, 1)
  expect_identical(fit$draws$b, fit$chains[[1]]$b)
  expect_identical(fit$draws$theta, fit$chains[[1]]$theta)
  expect_identical(fit$draws$delta, fit$chains[[1]]$delta)
  expect_identical(fit$draws$sigma_d2, fit$chains[[1]]$sigma_d2)
  expect_identical(fit$draws$w_vech, fit$chains[[1]]$w_vech)
  expect_identical(fit$draws$loglik, fit$chains[[1]]$loglik)
})

test_that("chains = 1 fixed-seed HMNP: draws are identical to chains[[1]] (bitwise back-compat)", {
  fx <- .hmnp_fx_phase1()
  fit <- suppressWarnings(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = 300, burn = 100, seed = 5))
  )
  expect_length(fit$chains, 1)
  expect_identical(fit$draws$b, fit$chains[[1]]$b)
  expect_identical(fit$draws$theta, fit$chains[[1]]$theta)
  expect_identical(fit$draws$delta, fit$chains[[1]]$delta)
  expect_identical(fit$draws$sigma_d2, fit$chains[[1]]$sigma_d2)
  expect_identical(fit$draws$w_vech, fit$chains[[1]]$w_vech)
})

# --- 5. Lever G (parallel chains) was DEFERRED after a libgomp+fork hazard --
# was found during Phase-1 validation (a later fit in the same R session
# after a cores > 1 fit could segfault). There is no `mcmc$cores` argument
# and no parallel dispatch -- chains run sequentially (tested in section 4
# above); no tests target it here by design.

# --- 6. summary() consolidated diagnostic table ------------------------------

test_that("summary() prints the consolidated diagnostics table (HMNL)", {
  fx <- .hmnl_fx_phase1()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2)
  )
  s <- summary(fit)
  expect_output(print(s), "Convergence diagnostics")
  expect_output(print(s), "ESS_bulk")
  expect_output(print(s), "R-hat")
  expect_output(print(s), "Acceptance:")
})

test_that("summary() prints the conjugate line and never warns about sigma^2 (raw) (HMNP)", {
  fx <- .hmnp_fx_phase1()

  warnings_seen <- character(0)
  fit <- withCallingHandlers(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("sigma\\^2 \\(raw\\)|sigma2", warnings_seen)))

  s <- summary(fit)
  expect_output(print(s), "Convergence diagnostics")
  expect_output(print(s), "ESS_bulk")
  expect_output(print(s), "Acceptance: conjugate")
  expect_output(print(s), "sigma\\^2 \\(raw\\)")
  expect_output(print(s), "non-identified")
})

# --- 7. traceplot() -----------------------------------------------------------

test_that("traceplot() runs onto a pdf device, returns invisibly, dispatches on choicer_hb", {
  fx <- .hmnl_fx_phase1()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2)
  )
  expect_true(inherits(fit, "choicer_hb"))

  pdf_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_path)
  on.exit({ grDevices::dev.off(); unlink(pdf_path) }, add = TRUE)

  ret <- traceplot(fit)
  expect_identical(ret, fit)

  ret2 <- traceplot(fit, block = c("b", "delta"))
  expect_identical(ret2, fit)

  expect_error(traceplot(fit, block = "bogus"), "must be one or more of")
})

test_that("traceplot() handles HMNP's sigma2 block and an explicit delta `which`", {
  fx <- .hmnp_fx_phase1(J = 5)
  fit <- suppressWarnings(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = 300, burn = 100), chains = 2)
  )

  pdf_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_path)
  on.exit({ grDevices::dev.off(); unlink(pdf_path) }, add = TRUE)

  ret <- traceplot(fit, block = c("sigma_d2", "sigma2"))
  expect_identical(ret, fit)

  lbl <- colnames(fit$draws$delta)[1:2]
  ret2 <- traceplot(fit, block = "delta", which = lbl)
  expect_identical(ret2, fit)
})

test_that("traceplot() falls back gracefully when object$chains is unavailable", {
  fx <- .hmnl_fx_phase1()
  fit <- suppressWarnings(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = 300, burn = 100))
  )
  fit$chains <- NULL

  pdf_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_path)
  on.exit({ grDevices::dev.off(); unlink(pdf_path) }, add = TRUE)

  expect_message(traceplot(fit), "object\\$chains.*unavailable")
})

# --- 8. keep_beta_i = "draws" memory guard (corrected 1.9x x chains formula) --

test_that("memory guard stop() fires above 4 GB and reports the 1.9x-chains formula", {
  fx <- .hmnl_fx_phase1(N = 10, T = 1, J = 3)
  K <- fx$d$K_struct
  N <- fx$d$N_persons
  chains <- 2L
  R <- 7000001L
  burn <- 1L
  thin <- 1L
  R_keep_est <- (R - burn + thin - 1L) %/% thin
  bytes_per_chain <- 8 * as.numeric(K) * N * R_keep_est * .HB_BETA_I_WRAP_FACTOR
  bytes_total <- bytes_per_chain * chains
  expect_gt(bytes_total, 4e9)

  expect_error(
    run_hmnlogit(input_data = fx$d, mcmc = list(R = R, burn = burn, thin = thin),
                 chains = chains, keep_beta_i = "draws"),
    regexp = sprintf("would allocate an estimated %.1f GB", bytes_total / 1e9),
    fixed = TRUE
  )
})

test_that("memory guard warning() fires between 1 GB and 4 GB with the corrected formula", {
  # The > 1 GB branch only warns (it does not stop()), so the wrapper would
  # normally continue on to actually run the sampler for R = 4,000,001
  # iterations -- far too slow (and memory-hungry) for a unit test. We only
  # need to confirm the warning ITSELF carries the correct message; a
  # calling handler that converts that specific warning into an error
  # aborts the wrapper immediately after it fires, before the (expensive)
  # Gibbs sampler is ever invoked.
  fx <- .hmnl_fx_phase1(N = 10, T = 1, J = 3)
  K <- fx$d$K_struct
  N <- fx$d$N_persons
  chains <- 1L
  R <- 4000001L
  burn <- 1L
  thin <- 1L
  R_keep_est <- (R - burn + thin - 1L) %/% thin
  bytes_per_chain <- 8 * as.numeric(K) * N * R_keep_est * .HB_BETA_I_WRAP_FACTOR
  bytes_total <- bytes_per_chain * chains
  expect_gt(bytes_total, 1e9)
  expect_lt(bytes_total, 4e9)

  expected_msg <- sprintf("allocates an estimated %.1f GB", bytes_total / 1e9)
  abort_after_guard_warning <- function(w) {
    if (grepl(expected_msg, conditionMessage(w), fixed = TRUE)) {
      stop("test-harness-abort-after-guard-warning")
    }
  }
  expect_error(
    withCallingHandlers(
      run_hmnlogit(input_data = fx$d, mcmc = list(R = R, burn = burn, thin = thin),
                   chains = chains, keep_beta_i = "draws"),
      warning = abort_after_guard_warning
    ),
    "test-harness-abort-after-guard-warning"
  )
})

test_that("memory guard also fires for run_hmnprobit with the same corrected formula", {
  fx <- .hmnp_fx_phase1(N = 10, T = 1, J = 3)
  K <- fx$d$K_struct
  N <- fx$d$N_persons
  chains <- 2L
  R <- 7000001L
  burn <- 1L
  thin <- 1L
  R_keep_est <- (R - burn + thin - 1L) %/% thin
  bytes_per_chain <- 8 * as.numeric(K) * N * R_keep_est * .HB_BETA_I_WRAP_FACTOR
  bytes_total <- bytes_per_chain * chains
  expect_gt(bytes_total, 4e9)

  expect_error(
    run_hmnprobit(input_data = fx$d, mcmc = list(R = R, burn = burn, thin = thin),
                  chains = chains, keep_beta_i = "draws"),
    regexp = sprintf("would allocate an estimated %.1f GB", bytes_total / 1e9),
    fixed = TRUE
  )
})
