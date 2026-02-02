# Tests for Cholesky decomposition and covariance matrix utilities:
# - build_L_mat()
# - build_var_mat()
# - jacobian_vech_Sigma()

test_that("build_L_mat creates lower-triangular matrix (correlated)", {
  K_w <- 3
  # L_params for correlated: diagonal log-params at positions 1, 3, 6
  # Off-diagonal at positions 2, 4, 5
  L_params <- c(log(1.0), 0.5, log(0.8), 0.3, -0.2, log(0.6))

  L <- build_L_mat(L_params, K_w, rc_correlation = TRUE)

  # Should be K_w x K_w

  expect_equal(dim(L), c(K_w, K_w))

  # Should be lower triangular (upper triangle is zero)
  expect_true(all(L[upper.tri(L)] == 0))

  # Diagonal elements should be exp(log-params) = positive
  expect_equal(L[1, 1], exp(log(1.0)))
  expect_equal(L[2, 2], exp(log(0.8)))
  expect_equal(L[3, 3], exp(log(0.6)))
  expect_true(all(diag(L) > 0))

  # Off-diagonal elements should match params directly
  expect_equal(L[2, 1], 0.5)
  expect_equal(L[3, 1], 0.3)
  expect_equal(L[3, 2], -0.2)
})

test_that("build_L_mat creates diagonal matrix (uncorrelated)", {
  K_w <- 3
  L_params <- c(log(1.0), log(0.5), log(0.8))

  L <- build_L_mat(L_params, K_w, rc_correlation = FALSE)

  # Should be K_w x K_w
  expect_equal(dim(L), c(K_w, K_w))

  # Should be diagonal (off-diagonal all zero)
  expect_true(all(L[row(L) != col(L)] == 0))

  # Diagonal should be exp(params)
  expect_equal(diag(L), exp(L_params))
})

test_that("build_var_mat produces symmetric positive definite matrix", {
  K_w <- 2
  L_params <- c(log(1.0), 0.3, log(0.5))

  Sigma <- build_var_mat(L_params, K_w, rc_correlation = TRUE)

  # Should be K_w x K_w
  expect_equal(dim(Sigma), c(K_w, K_w))

  # Should be symmetric
  expect_equal(Sigma, t(Sigma), tolerance = 1e-12)

  # Should be positive definite (all eigenvalues > 0)
  eig <- eigen(Sigma, only.values = TRUE)$values
  expect_true(all(eig > 0))
})

test_that("build_var_mat equals L %*% t(L)", {
  K_w <- 3
  L_params <- c(log(1.0), 0.5, log(0.8), 0.3, -0.2, log(0.6))

  L <- build_L_mat(L_params, K_w, rc_correlation = TRUE)
  Sigma <- build_var_mat(L_params, K_w, rc_correlation = TRUE)

  # Sigma should equal L %*% t(L)
  expected_Sigma <- L %*% t(L)
  expect_equal(Sigma, expected_Sigma, tolerance = 1e-12)
})

test_that("build_var_mat diagonal case produces diagonal variances", {
  K_w <- 3
  L_params <- c(log(0.5), log(1.0), log(0.8))

  Sigma <- build_var_mat(L_params, K_w, rc_correlation = FALSE)

  # Should be diagonal
  expect_true(all(Sigma[row(Sigma) != col(Sigma)] == 0))

  # Diagonal should be exp(L_params)^2 = exp(2*L_params)
  expect_equal(diag(Sigma), exp(L_params)^2, tolerance = 1e-12)
})

test_that("jacobian_vech_Sigma is numerically accurate", {
  skip_if_not_installed("numDeriv")

  K_w <- 2
  L_params <- c(log(0.8), 0.2, log(0.6))

  # Analytical Jacobian from C++
  J_anal <- jacobian_vech_Sigma(L_params, K_w, rc_correlation = TRUE)

  # Numerical Jacobian using numDeriv
  vech_sigma_fn <- function(lp) {
    S <- build_var_mat(lp, K_w, rc_correlation = TRUE)
    # Extract lower triangle including diagonal (column-major order)
    S[lower.tri(S, diag = TRUE)]
  }

  J_num <- numDeriv::jacobian(vech_sigma_fn, L_params)

  expect_equal(J_anal, J_num, tolerance = TOL_GRAD)
})

test_that("jacobian_vech_Sigma works for larger K_w", {
  skip_if_not_installed("numDeriv")

  K_w <- 3
  L_size <- K_w * (K_w + 1) / 2  # 6

  set.seed(42)
  L_params <- runif(L_size, -0.5, 0.5)

  J_anal <- jacobian_vech_Sigma(L_params, K_w, rc_correlation = TRUE)

  vech_sigma_fn <- function(lp) {
    S <- build_var_mat(lp, K_w, rc_correlation = TRUE)
    S[lower.tri(S, diag = TRUE)]
  }

  J_num <- numDeriv::jacobian(vech_sigma_fn, L_params)

  expect_equal(J_anal, J_num, tolerance = TOL_GRAD)
})

test_that("jacobian_vech_Sigma works for diagonal case", {
  skip_if_not_installed("numDeriv")

  K_w <- 3
  L_params <- c(log(0.5), log(1.0), log(0.8))

  J_anal <- jacobian_vech_Sigma(L_params, K_w, rc_correlation = FALSE)

  vech_sigma_fn <- function(lp) {
    S <- build_var_mat(lp, K_w, rc_correlation = FALSE)
    S[lower.tri(S, diag = TRUE)]
  }

  J_num <- numDeriv::jacobian(vech_sigma_fn, L_params)

  expect_equal(J_anal, J_num, tolerance = TOL_GRAD)
})

test_that("build_L_mat handles K_w = 1", {
  K_w <- 1
  L_params <- log(0.7)

  L <- build_L_mat(L_params, K_w, rc_correlation = TRUE)

  expect_equal(dim(L), c(1, 1))
  expect_equal(L[1, 1], exp(L_params))

  L_uncorr <- build_L_mat(L_params, K_w, rc_correlation = FALSE)
  expect_equal(L, L_uncorr)
})

test_that("build_var_mat handles K_w = 1", {
  K_w <- 1
  L_params <- log(0.7)

  Sigma <- build_var_mat(L_params, K_w, rc_correlation = TRUE)

  expect_equal(dim(Sigma), c(1, 1))
  expect_equal(Sigma[1, 1], exp(L_params)^2, tolerance = 1e-12)
})

test_that("build_L_mat parameter ordering is correct", {
  # For K_w = 3 with correlation:
  # L = [ exp(p1)  0       0     ]
  #     [ p2       exp(p3) 0     ]
  #     [ p4       p5      exp(p6)]

  K_w <- 3
  p1 <- log(1); p2 <- 0.1; p3 <- log(2)
  p4 <- 0.2; p5 <- 0.3; p6 <- log(3)
  L_params <- c(p1, p2, p3, p4, p5, p6)

  L <- build_L_mat(L_params, K_w, rc_correlation = TRUE)

  expect_equal(L[1, 1], exp(p1))
  expect_equal(L[2, 1], p2)
  expect_equal(L[2, 2], exp(p3))
  expect_equal(L[3, 1], p4)
  expect_equal(L[3, 2], p5)
  expect_equal(L[3, 3], exp(p6))
})
