# Setup environment ------------------------------------------------------------
rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(nloptr)
library(logitr)

devtools::load_all()

# Simulation settings ----------------------------------------------------------
N <- 1e+4              # Number of individuals
J_global <- 6          # Number of total alternatives (excluding outside good)
S <- 50 * ceiling((1.5 * sqrt(N)) / 50)     # Number of simulation draws per individual; S > sqrt(N) for consistency

# Define seed for reproducibility
set.seed(123)

# Simulate DGP -----------------------------------------------------------------
delta_true <- rep(c(0.5, -0.5), J_global/2)
beta_true <- c(0.8, -0.6)
mu_true <- c(0.3, 0.7)
Sigma_true <- matrix(
  c(1.5, 0.5,
    0.5, 1.0),
  nrow = 2
)

# Sigma_true <- matrix(
#   c(1.0, 0.0,
#     0.0, 1.0),
#   nrow = K_w
# )

# First random coefficient is log-normal, second is normal
rc_dist <- c(1L, 0L)

K_x <- length(beta_true)                # Number of variables with "fixed" coefficients
K_w <- ncol(Sigma_true)                 # Number of random-coefficient variables

if (is.null(Sigma_true[2,1]) || is.na(Sigma_true[2,1]) || Sigma_true[2,1] == 0) {
  rc_correlation <- FALSE
} else {
  rc_correlation <- TRUE
}

# Cholesky decomposition of Sigma_true
L_true <- t(chol(Sigma_true))
L_params_true <- c(log(L_true[1,1]),L_true[2,1],log(L_true[2,2]))
L_size <- if (rc_correlation) {
  K_w * (K_w + 1) / 2
} else {
  K_w
}

# Simulate gamma_i for each individual
gamma_i <- L_true  %*% matrix(rnorm(N * K_w), nrow = K_w, ncol = N)

# transform to log-normal distribution for selected dimensions
for (r in seq_len(nrow(gamma_i)) ){
  if (rc_dist[r] == 1) {
    # Ensure it's negative
    gamma_i[r, ] <- exp(gamma_i[r, ])
  } else if (rc_dist[r] != 0L) {
    warning("rc_dist must take values either 0 or 1.")
  }
}

# Add mean shifter to random coefficients
gamma_i <- kronecker(matrix(mu_true, ncol=1), matrix(1, nrow=1, ncol=N)) + gamma_i
gamma_i <- t(gamma_i)

# Draw alternatives
choice_set_sizes <- sample(2:J_global, N, replace = TRUE)
J_sum <- sum(choice_set_sizes)
alt_indices <- vector("list", length = N)
for (i in 1:N) {
  alt_indices[[i]] <- sample(J_global, choice_set_sizes[i]) |> sort()
}

# Create data.table
dt <- data.table(id = rep(seq_len(N), times = choice_set_sizes))

dt[, `:=`(
  w1 = -runif(.N, min = 0.1, max = 3), # NEGATIVE PRICE
  w2 = runif(.N, min = -1, max = 1),
  x1 = runif(.N, min = 0, max = 3),
  x2 = runif(.N, min = -1, max = 1)
)]

dt[, `:=`(
  gamma1 = gamma_i[id, 1],
  gamma2 = gamma_i[id, 2]
)]

# Add alternatives and deltas
dt[, `:=`(
  alt = alt_indices[[id]],
  delta = delta_true[alt_indices[[id]]]
),
by = id
]

# Add outside option
dt <- dt[
  , .(alt = 0L),
  keyby = id
] |>
  list(dt) |>
  rbindlist(use.names = TRUE, fill = TRUE)

# Set variables to zero for outside option
cols <- setdiff(names(dt), c("id","alt","choice"))
dt[alt == 0L, (cols) := 0]

# Set key for data.table
setkey(dt, id, alt)

# Logit error shock
dt[, epsilon := -log(-log(runif(.N)))]

# Compute utility per individual and alternative
dt[, utility := delta + x1 * beta_true[1] + x2 * beta_true[2] + w1 * (gamma1) + w2 * (gamma2) + epsilon]

# Find choice
dt[, choice := fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]

# logitr -----------------------------------------------------------------------

# library(logitr)

# # Optimization settings
# nloptr_opts <- list(
#   "algorithm" = "NLOPT_LD_LBFGS",
#   "xtol_rel" = 1.0e-8,
#   "maxeval" = 1e+3,
#   "print_level" = 1L,
#   "check_derivatives" = TRUE,
#   "check_derivatives_print" = "none"
# )

# dt[, alt_factor := factor(alt)]

# halton_seq <- randtoolbox::halton(n = S, dim = J_global + K_x + K_w, normal = TRUE)

# logitr_test <- logitr(
#   data = dt,
#   outcome = "choice",
#   obsID = "id",
#   pars = c("x1","x2","w1","w2","alt_factor"),
#   randPars = c("w1" = "ln", "w2" = "n"),
#   correlation = TRUE,
#   numDraws = S,
#   standardDraws = halton_seq,
#   options = nloptr_opts
# )

# logitr_test |> summary() |> print()

# choicer Version --------------------------------------------------------------

# Generate eta_draws with dimensions K_w x S x N
eta_draws <- get_halton_normals(S, N, K_w)

mxl_inputs <- prepare_mxl_data(
  data = dt,
  id_col = "id",
  alt_col = "alt",
  choice_col = "choice",
  covariate_cols = c("x1", "x2"),
  random_var_cols = c("w1", "w2"),
  outside_opt_label = 0L,
  include_outside_option = FALSE,
  rc_correlation = rc_correlation
)

nloptr_opts <- list(
  "algorithm" = "NLOPT_LD_LBFGS",
  "xtol_rel" = 1.0e-8,
  "maxeval" = 1e+3,
  "print_level" = 1L,
  "check_derivatives" = TRUE,
  "check_derivatives_print" = "none"
)

result_mnl <- nloptr::nloptr(
  x0 = rep(0, ncol(mxl_inputs$X) + nrow(mxl_inputs$alt_mapping) - 1),
  eval_f = mnl_loglik_gradient_parallel,
  opts = nloptr_opts,
  X = mxl_inputs$X,
  alt_idx = mxl_inputs$alt_idx,
  choice_idx = mxl_inputs$choice_idx,
  M = mxl_inputs$M,
  weights = mxl_inputs$weights,
  use_asc = TRUE,
  include_outside_option = mxl_inputs$include_outside_option
)

# Debug: Verify dimensions before optimization
cat("=== Dimension Check ===\n")
cat("N from simulation:", N, "\n")
cat("N from data:", mxl_inputs$N, "\n")
cat("N for eta_draws:", dim(eta_draws)[3], "\n")
cat("K_x (from X):", ncol(mxl_inputs$X), "\n")
cat("K_w (from W):", ncol(mxl_inputs$W), "\n")
cat("J (from alt_mapping):", nrow(mxl_inputs$alt_mapping), "\n")
expected_len <- ncol(mxl_inputs$X) + K_w + L_size + nrow(mxl_inputs$alt_mapping) - 1
cat("Expected theta length:", expected_len, "\n")
cat("=======================\n\n")

# Run model estimation (let run_mxlogit compute theta_init for safety)
result <- run_mxlogit(
  input_data = mxl_inputs,
  eta_draws = eta_draws,
  rc_dist = rc_dist,
  rc_mean = TRUE,
  use_asc = TRUE,
  # theta_init = runif(ncol(mxl_inputs$X) + K_w + L_size + nrow(mxl_inputs$alt_mapping) - 1, -1, 1),
  theta_init = c(
    result_mnl$solution[1:ncol(mxl_inputs$X)],
    rep(0, K_w + L_size),
    result_mnl$solution[(ncol(mxl_inputs$X)+1):length(result_mnl$solution)]
  ),
  nloptr_opts = nloptr_opts
)

# Compare likelihoods
llk_true <- mxl_loglik_gradient_parallel(
  theta = c(beta_true, mu_true, L_params_true, delta_true),
  X = mxl_inputs$X,
  W = mxl_inputs$W,
  alt_idx = mxl_inputs$alt_idx,
  choice_idx = mxl_inputs$choice_idx,
  M = mxl_inputs$M,
  weights = mxl_inputs$weights,
  eta_draws = eta_draws,
  rc_dist = rc_dist,
  rc_correlation = rc_correlation,
  rc_mean = TRUE,
  use_asc = TRUE,
  include_outside_option = mxl_inputs$include_outside_option
)

cat("Log Likelihood at true parameters:", -llk_true$objective, "\n")
cat("Log Likelihood at estimated parameters:", -result$objective, "\n\n")

# Extract the estimated parameters
theta_est <- result$solution

# Separate the estimated parameters
beta_est <- theta_est[1:K_x]
mu_est <- theta_est[(K_x+1):(K_x+K_w)]
L_params_est <- theta_est[(K_x+K_w + 1):(K_x+K_w + L_size)]
delta_est <- theta_est[(K_x+K_w + L_size + 1):length(theta_est)]

# Compute estimated Sigma matrix
Sigma_est <- build_var_mat(L_params_est, K_w, rc_correlation)

# Print true and estimated beta
cat("True beta:", beta_true, "\n")
cat("Estimated beta:", beta_est, "\n\n")

# Print true and estimated mu
cat("True mu:", mu_true, "\n")
cat("Estimated mu:", c(exp(mu_est[1]),mu_est[2]), "\n\n")

# Print true and estimated Sigma
cat("True Sigma:\n")
print(Sigma_true)
cat("Estimated Sigma:\n")
print(Sigma_est)

# Print true and estimated delta
cat("\nTrue delta:", delta_true, "\n")
cat("Estimated delta:", delta_est, "\n\n")

# Get result table
get_mxl_result(
  result,
  X = mxl_inputs$X,
  W = mxl_inputs$W,
  alt_idx = mxl_inputs$alt_idx,
  choice_idx = mxl_inputs$choice_idx,
  M = mxl_inputs$M,
  weights = mxl_inputs$weights,
  eta_draws = eta_draws,
  rc_dist = rc_dist,
  rc_correlation = rc_correlation,
  rc_mean = TRUE,
  use_asc = TRUE,
  include_outside_option = mxl_inputs$include_outside_option,
  omit_asc_output = FALSE,
  param_names = NULL
)
