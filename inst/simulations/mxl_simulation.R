# Mixed Logit - Solver Hardening Workflow
#
# Demonstrates the recommended hardening for production specifications:
#   1. Warm-up MNL fit on the same X block -> structured theta_init for MXL.
#   2. Cholesky-diagonal bounds via the named-vector form of `lower`/`upper`.
#   3. `scale_vars = "sd"` for Hessian conditioning across blocks.
#
# Run from package root: Rscript inst/simulations/mxl_simulation.R

library(choicer)

# 1) DGP ======================================================================
sim <- simulate_mxl_data(N = 2e4, J = 8, seed = 123)
print(sim)

# 2) Data preparation =========================================================
mxl_inputs <- prepare_mxl_data(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  random_var_cols        = c("w1", "w2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE,
  rc_correlation         = TRUE
)

# S > sqrt(N) for MSL consistency
S <- 50 * ceiling((1.5 * sqrt(sim$settings$N)) / 50)
eta_draws <- get_halton_normals(S, mxl_inputs$N, sim$settings$K_w)
cat("Simulation draws: S =", S, "\n\n")

# 3) Warm-up MNL fit ==========================================================
# Fit a fixed-coefficient MNL on the same X block. The MNL coefficients give
# a near-optimum starting point for beta and the ASCs, which keeps the first
# L-BFGS step in the MXL small. Cold-starting from zeros on a correlated-RC
# spec can land the first step somewhere the simulated likelihood evaluates
# to NaN, and L-BFGS's line-search needs a finite reference value to backtrack
# - the warm-up sidesteps that failure mode.
cat("--- Warm-up MNL fit ---\n")
mnl_fit <- run_mnlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "alt",
  choice_col             = "choice",
  covariate_cols         = c("x1", "x2"),
  outside_opt_label      = 0L,
  include_outside_option = FALSE
)
mnl_coef <- coef(mnl_fit)

# 4) Build structured theta_init for MXL ======================================
# Splice MNL coefficients into the MXL parameter vector block-by-block. The
# MXL ordering is (beta, [mu], L, ASC); the MNL ordering is (beta, ASC).
# Default `log(0.5)` on the Cholesky diagonal (L_pp = 0.5, RC variance 0.25)
# is a moderate prior - much safer than `0` (L_pp = 1, unit variance), which
# can let the first L-BFGS step push ell_pp toward -Inf.
K_x    <- ncol(mxl_inputs$X)
K_w    <- ncol(mxl_inputs$W)
J      <- nrow(mxl_inputs$alt_mapping)
L_size <- K_w * (K_w + 1) / 2          # rc_correlation = TRUE
n_asc  <- J - 1
n_params <- K_x + L_size + n_asc

theta_init <- numeric(n_params)

# beta block
theta_init[seq_len(K_x)] <- mnl_coef[seq_len(K_x)]

# Cholesky block: log(0.5) on diagonals, 0 on off-diagonals
sigma_idx          <- K_x + seq_len(L_size)
diag_pos_in_sigma  <- cumsum(seq_len(K_w))   # row-major diag positions
theta_init[sigma_idx[diag_pos_in_sigma]] <- log(0.5)

# ASC block: MNL ASCs (last n_asc entries of mnl_coef)
asc_idx <- K_x + L_size + seq_len(n_asc)
theta_init[asc_idx] <- mnl_coef[K_x + seq_len(n_asc)]

# 5) Cholesky bounds (named partial form) =====================================
# Clip ell_pp in [-5, 5] so each L_pp in [exp(-5), exp(5)] = [0.007, 148].
# Keeps the optimizer from diverging the diagonal early without constraining
# the eventual MLE. Bounds are passed in natural units; the package forward-
# transforms them into scaled space internally.
L_diag_names <- paste0("L_", seq_len(K_w), seq_len(K_w))
lower <- setNames(rep(-5, K_w), L_diag_names)
upper <- setNames(rep( 5, K_w), L_diag_names)

# 6) Final MXL fit ============================================================
cat("\n--- Final MXL fit ---\n")
fit <- run_mxlogit(
  input_data = mxl_inputs,
  eta_draws  = eta_draws,
  use_asc    = TRUE,
  theta_init = theta_init,
  lower      = lower,
  upper      = upper,
  scale_vars = "sd",
  se_method  = "bhhh",
  control    = list(print_level = 1L)
)

cat("\n")
summary(fit)

# 7) Parameter Recovery =======================================================
cat("\n--- Parameter Recovery ---\n")
print(recovery_table(fit, sim$true_params))
