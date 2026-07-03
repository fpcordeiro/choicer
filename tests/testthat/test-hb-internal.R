# Tests for the shared C++ HB infrastructure (src/hb_internal.h) via the
# hb_test_* shims in src/hb_test_exports.cpp:
#   choicer:::hb_test_chol(A)
#   choicer:::hb_test_trisolve(L, b, transpose)
#   choicer:::hb_test_spd_solve(A, b)
#   choicer:::hb_test_logsumexp(v, include_outside)
#   choicer:::hb_test_sigma_d2_gibbs(xi, n_iter, seed, half_cauchy, s_d, c0, d0)
# The hand-rolled Cholesky/trisolves are the worker-thread-safe replacements
# for LAPACK inside the OpenMP region, so they are validated against the
# LAPACK-backed base R references here.

test_that("hand-rolled Cholesky matches base chol() on random SPD matrices", {
  set.seed(42)
  for (K in c(1L, 2L, 4L, 7L)) {
    A <- crossprod(matrix(rnorm(K * K), K, K)) + diag(K)
    res <- choicer:::hb_test_chol(A)
    expect_true(res$ok)
    expect_equal(res$L, t(chol(A)), tolerance = 1e-10)
    # exact reconstruction A = L L'
    expect_equal(res$L %*% t(res$L), A, tolerance = 1e-10)
  }
})

test_that("hand-rolled Cholesky reports failure on non-PD input", {
  A_indef <- matrix(c(1, 2, 2, 1), 2, 2)     # eigenvalues 3 and -1
  expect_false(choicer:::hb_test_chol(A_indef)$ok)

  A_singular <- matrix(1, 3, 3)              # rank 1
  expect_false(choicer:::hb_test_chol(A_singular)$ok)

  A_nan <- diag(2); A_nan[2, 2] <- NaN       # non-finite pivot
  expect_false(choicer:::hb_test_chol(A_nan)$ok)
})

test_that("hand-rolled triangular solves match forwardsolve/backsolve", {
  set.seed(7)
  K <- 6
  A <- crossprod(matrix(rnorm(K * K), K, K)) + K * diag(K)
  L <- t(chol(A))
  b <- rnorm(K)
  expect_equal(as.vector(choicer:::hb_test_trisolve(L, b, FALSE)),
               forwardsolve(L, b), tolerance = 1e-10)
  expect_equal(as.vector(choicer:::hb_test_trisolve(L, b, TRUE)),
               backsolve(t(L), b), tolerance = 1e-10)
})

test_that("SPD solve via two trisolves matches base solve()", {
  set.seed(11)
  for (K in c(1L, 3L, 5L)) {
    A <- crossprod(matrix(rnorm(K * K), K, K)) + K * diag(K)
    b <- rnorm(K)
    res <- choicer:::hb_test_spd_solve(A, b)
    expect_true(res$ok)
    expect_equal(as.vector(res$x), solve(A, b), tolerance = 1e-10)
  }
  # failure path: not SPD -> ok = FALSE, never an exception
  res_bad <- choicer:::hb_test_spd_solve(matrix(c(1, 2, 2, 1), 2, 2), c(1, 1))
  expect_false(res_bad$ok)
})

test_that("fixed-order LSE matches log(sum(exp())) incl. the implicit outside term", {
  set.seed(3)
  v <- rnorm(9, sd = 2)
  expect_equal(choicer:::hb_test_logsumexp(v, FALSE), log(sum(exp(v))),
               tolerance = 1e-12)
  expect_equal(choicer:::hb_test_logsumexp(v, TRUE), log(1 + sum(exp(v))),
               tolerance = 1e-12)

  # overflow-safe at large utilities
  big <- c(1000, 1000.5)
  expect_equal(choicer:::hb_test_logsumexp(big, FALSE),
               1000.5 + log(1 + exp(-0.5)), tolerance = 1e-12)
  expect_equal(choicer:::hb_test_logsumexp(big, TRUE),
               1000.5 + log(exp(-1000.5) + 1 + exp(-0.5)), tolerance = 1e-12)

  # empty inside range: outside-only task has denominator exp(0)
  expect_equal(choicer:::hb_test_logsumexp(numeric(0), TRUE), 0)
  expect_identical(choicer:::hb_test_logsumexp(numeric(0), FALSE), -Inf)

  # NaN propagates, +Inf dominates
  expect_true(is.nan(choicer:::hb_test_logsumexp(c(0, NaN), TRUE)))
  expect_identical(choicer:::hb_test_logsumexp(c(0, Inf), FALSE), Inf)
})

test_that("sigma_d2 scale-mixture Gibbs is positive and in the right ballpark", {
  set.seed(99)
  J <- 200
  sigma_d_true <- 0.7
  xi <- rnorm(J, 0, sigma_d_true)
  target <- mean(xi^2)      # the likelihood concentrates sigma_d2 near this

  # half-Cauchy(0, 1) via the Makalic-Schmidt two-block mixture
  draws <- choicer:::hb_test_sigma_d2_gibbs(
    xi, n_iter = 2000L, seed = 42, half_cauchy = TRUE,
    s_d = 1, c0 = 3, d0 = 3
  )
  expect_length(draws, 2000L)
  expect_true(all(is.finite(draws)))
  expect_true(all(draws > 0))
  post <- draws[-(1:200)]
  expect_gt(mean(post), 0.5 * target)
  expect_lt(mean(post), 1.5 * target)

  # IG(c0, d0) fallback branch
  draws_ig <- choicer:::hb_test_sigma_d2_gibbs(
    xi, n_iter = 2000L, seed = 42, half_cauchy = FALSE,
    s_d = 1, c0 = 3, d0 = 3
  )
  expect_true(all(draws_ig > 0))
  post_ig <- draws_ig[-(1:200)]
  expect_gt(mean(post_ig), 0.5 * target)
  expect_lt(mean(post_ig), 1.5 * target)

  # deterministic given the master seed (own Xoshiro streams, not R's RNG)
  again <- choicer:::hb_test_sigma_d2_gibbs(
    xi, n_iter = 2000L, seed = 42, half_cauchy = TRUE,
    s_d = 1, c0 = 3, d0 = 3
  )
  expect_identical(draws, again)
})
