# Regression tests for the row-major Cholesky packing convention shared by
# build_L_mat() / jacobian_vech_Sigma() (C++) and the R-side simulator and
# label generator. Guards against the K_w >= 3 bug where simulator and
# optimizer disagreed on which theta-slot held which Cholesky entry.

test_that("simulator's L_params round-trip through build_var_mat for K_w in 2..4", {
  for (K_w in 2:4) {
    set.seed(K_w)
    A <- matrix(stats::rnorm(K_w * K_w), K_w, K_w)
    Sigma <- crossprod(A) + diag(0.1, K_w)
    sim <- simulate_mxl_data(
      N = 100L, J = 3L,
      beta = c(0.1, -0.1),
      Sigma = Sigma, seed = 1L,
      outside_option = TRUE, vary_choice_set = FALSE
    )
    expect_equal(
      build_var_mat(sim$true_params$L_params, K_w, TRUE),
      Sigma, tolerance = 1e-10, ignore_attr = TRUE
    )
  }
})

test_that("vech_row matches jacobian_vech_Sigma row ordering at K_w=3", {
  skip_if_not_installed("numDeriv")
  K_w <- 3L
  L_params <- c(log(1.0), 0.4, log(1.158), 0.2, 0.190, log(0.851))
  J_ana <- jacobian_vech_Sigma(L_params, K_w, TRUE)
  J_num <- numDeriv::jacobian(
    function(p) choicer:::vech_row(build_var_mat(p, K_w, TRUE)),
    L_params
  )
  expect_equal(J_ana, J_num, tolerance = 1e-6, ignore_attr = TRUE)
})
