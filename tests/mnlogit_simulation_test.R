# Setup environment ------------------------------------------------------------
rm(list = ls(all.names = TRUE))

library(data.table)
library(nloptr)
library(choicer)

# Simulation settings ----------------------------------------------------------
N <- 1e+4                 # Number of choice situations
J_global <- 10            # Number of total alternatives (excluding outside good)

# For reproducibility
set.seed(123)

# Simulate DGP -----------------------------------------------------------------

# True parameter values
delta_true <- rep(c(0.5, -0.5), J_global/2)
beta_true <- c(0.8, -0.6)
K_x <- length(beta_true)        # Number of covariate variables in W

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
  x1 = runif(.N, min = -1, max = 1),
  x2 = runif(.N, min = -1, max = 1)
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
dt[, utility := delta + x1 * beta_true[1] + x2 * beta_true[2] + epsilon]

# Find choice
dt[, choice_id:= fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]

# choicer estimation ----------------------------------------------------------

cat("Starting choicer.\n\n")

# Optimization settings
nloptr_opts <- list(
  "algorithm" = "NLOPT_LD_LBFGS",
  "xtol_rel" = 1.0e-8,
  "maxeval" = 1e+3,
  "print_level" = 0L,
  "check_derivatives" = TRUE,
  "check_derivatives_print" = "none"
)

input_list <- prepare_mnl_data(
  data = dt,
  id_col = "id",
  alt_col = "alt",
  choice_col = "choice_id",
  covariate_cols = c("x1", "x2"),
  outside_opt_label = 0L,
  include_outside_option = FALSE
)

gc()

# Initial parameter vector theta_init
theta_init <- runif(J_global + K_x, -1, 1)

logit_test <- mnl_loglik_gradient_parallel(
  theta = theta_init,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  choice_idx = input_list$choice_idx,
  M = input_list$M,
  weights = input_list$weights,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option
)

# Run likelihood maximization
result <- nloptr::nloptr(
  x0 = theta_init,
  eval_f = mnl_loglik_gradient_parallel,
  opts = nloptr_opts,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  choice_idx = input_list$choice_idx,
  weights = input_list$weights,
  M = input_list$M,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option
)

# Extract the estimated parameters
theta_est <- result$solution

# Separate the estimated parameters
beta_est <- theta_est[1:K_x]
delta_est <- theta_est[(K_x+1):length(theta_est)]

# Log-likelihood
cat("Log-likelihood:", -result$objective, "\n\n")

# Print true and estimated delta
cat("True delta:", delta_true, "\n")
cat("Estimated delta:", delta_est, "\n\n")

# Print true and estimated beta
cat("True beta:", beta_true, "\n")
cat("Estimated beta:", beta_est, "\n\n")

param_names <- c("x1", "x2", as.character(input_list$alt_mapping[alt != 0L, alt]))

get_mnl_result(
  opt_result = result,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  choice_idx = input_list$choice_idx,
  M = input_list$M,
  weights = input_list$weights,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option,
  param_names = param_names,
  omit_asc_output = FALSE
)

# Post-Estimation --------------------------------------------------------------

# Predict individual-level choice probabilities
model_individual_predict <- mnl_predict(
  theta = result$solution,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  M = input_list$M,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option
)

# Predict market shares
model_shares_predict <- mnl_predict_shares(
  theta = result$solution,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  M = input_list$M,
  weights = input_list$weights,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option
)

dt[, `:=`(
  model_logit_prob = model_individual_predict$choice_prob,
  model_v = model_individual_predict$utility
)]

dt[, .(
  model_v = mean(model_v),
  data_v = mean(utility - epsilon),
  model_mkt_share = mean(model_logit_prob),
  data_mkt_share = mean(choice_id)
),
keyby = alt
]

alt_mapping <- input_list$alt_mapping[] |> copy()

alt_mapping[, MODEL_MKT_SHARE := model_shares_predict]

# BLP contraction --------------------------------------------------------------

delta0 <- rep(0, nrow(alt_mapping))

delta_contraction <- blp_contraction(
  delta = delta0,
  target_shares = alt_mapping$MKT_SHARE,
  X = input_list$X,
  beta = result$solution[1:K_x],
  alt_idx = input_list$alt_idx,
  M = input_list$M,
  weights = input_list$weights,
  include_outside_option = input_list$include_outside_option,
  tol = 1e-10
)

alt_mapping[, `:=`(
  delta_mnl_est = c(0, delta_est),
  delta_contr = delta_contraction
)]

if (exists('logitr_test')) {
  coefs <- logitr_test$coefficients
  coefs_len <- length(coefs)
  alt_mapping[, delta_logitr := c(0, coefs[(K_x + 1):coefs_len])]
  alt_mapping[, delta_diff := delta_mnl_est - delta_logitr]
}
