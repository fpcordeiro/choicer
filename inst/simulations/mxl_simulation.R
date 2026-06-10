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

# 8) Post-estimation ==========================================================
# x2 plays the role of price: it has a *fixed* coefficient (true beta_x2 =
# -0.6 < 0), as wtp() and consumer_surplus() require - a random price
# coefficient would make the WTP ratio ill-defined (no finite moments).

# Goodness of fit (equal-shares null; this DGP draws varying choice sets,
# so the closed-form market-shares null does not apply)
cat("\n--- Goodness of fit ---\n")
print(gof(fit))

# Predicted market shares (simulated over the Halton draws)
cat("\n--- Predicted shares ---\n")
print(predict(fit, type = "shares"))

# Elasticities and diversion ratios. Unlike MNL, MXL diversion is not
# share-proportional and depends on the perturbed covariate (the random
# coefficients break IIA), so wrt_var is required.
cat("\n--- Elasticities (x2) ---\n")
print(elasticities(fit, elast_var = "x2"))
cat("\n--- Diversion ratios (w.r.t. x2) ---\n")
print(diversion_ratios(fit, wrt_var = "x2"))

# Willingness to pay for x1 in units of x2 with delta-method SEs.
# True WTP_x1 = -0.8 / (-0.6) = 1.33. The random coefficients (w1, w2) were
# fitted without estimated means (rc_mean = FALSE, normal), so their mean WTP
# is 0 by construction and they are excluded with a note.
cat("\n--- Willingness to pay ---\n")
print(wtp(fit, price_var = "x2"))

# Expected consumer surplus from the simulated logsum: the per-draw
# log-sum-exp is averaged across draws (averaging utilities first would
# understate the logsum by Jensen's inequality). Point estimate for MXL.
cs0 <- consumer_surplus(fit, price_var = "x2")
cat("\n--- Consumer surplus (baseline) ---\n")
print(cs0)

# Counterfactual: raise x2 ("price") by 0.5 for alternative 2 only ============
dt_cf <- data.table::copy(sim$data)[alt == 2, x2 := x2 + 0.5]

cat("\n--- Counterfactual shares (x2 + 0.5 for alt 2) ---\n")
shares_cf  <- predict(fit, type = "shares", newdata = dt_cf)
shares_cmp <- rbind(baseline       = drop(predict(fit, type = "shares")),
                    counterfactual = drop(shares_cf))
colnames(shares_cmp) <- as.character(fit$alt_mapping[[2]])
print(shares_cmp)

# Welfare change: negative, since the price of one alternative rose
cs1 <- consumer_surplus(fit, price_var = "x2", newdata = dt_cf)
cat("\nDelta mean consumer surplus:", round(cs1$mean_cs - cs0$mean_cs, 4),
    "(negative: the price increase lowers expected surplus)\n")
