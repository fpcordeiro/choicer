# Setup environment ------------------------------------------------------------
rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(mlogit)
library(choicer)

set.seed(42)

# Simulation Parameters --------------------------------------------------------

# Number of decision makers
N <- 1e+5

# True parameters for systematic utility
beta_X <- 1.5
beta_W <- -0.8

# True alternative-specific constants (deltas)
# delta_1 is normalized to 0
deltas <- c(
  "1" = 0.0,
  "2" = 0.3,
  "3" = -0.2,
  "4" = -0.5,
  "5" = 0.4
)

# Nest structure and dissimilarity parameters (lambdas)
# 0 < lambda <= 1. A lambda of 1 implies no correlation (standard logit).
lambda_1 <- 0.6 # Nest 1: {1, 2}
lambda_2 <- 0.3 # Nest 2: {3, 4, 5}

## Create Data Structure & Systematic Utility ----------------------------------

# Create the main data.table
# We have 6 alternatives: j=0 (outside) and j=1...5 (inside)
DT <- CJ(id = 1:N, j = 1:5)

# Define the nests
# Nest 1: {1, 2}
# Nest 2: {3, 4, 5}
DT[, nest := fcase(
  j %in% c(1, 2), 1,
  j %in% c(3, 4, 5), 2
)]

# Assign lambdas
lambda_vec_all <- c(lambda_1, lambda_2)
DT[, lambda := lambda_vec_all[nest]]

# Add deltas (ASCs)
delta_vec_all <- c(deltas)
DT[, delta := delta_vec_all[j]]

# Simulate covariates
# Covariates are 0 for the outside option
DT[, X := 0][j > 0, X := rnorm(N * 5, mean = 1, sd = 1)]
DT[, W := 0][j > 0, W := runif(N * 5, min = -2, max = 2)]

# Calculate systematic utility (V_ij)
DT[, V := delta + beta_X * X + beta_W * W]

# Simulate GEV Errors (epsilon_ij) ---------------------------------------------

# Special case for lambda=0.5
DT_nests <- DT[, .(id, nest)] |> unique()
DT_nests[, nest_shock := -log(sqrt(2) * qnorm(1 - runif(.N)*0.5))]
DT[
  DT_nests,
  nest_shock := i.nest_shock,
  on = .(id, nest)
]
DT[, epsilon := nest_shock - 0.5 * log(-log(runif(.N)))]

# Calculate Final Utility and Simulate Choices ---------------------------------

# Final utility is the sum of systematic and random components
DT[, U := V + epsilon]

# Find the chosen alternative for each individual
# This is the alternative 'j' that provides the maximum utility
DT[, choice:= fifelse(seq_len(.N) == which.max(U), 1L, 0L), by = id]

# Display the simulated choice shares
cat("\n--- Simulated Choice Shares ---\n")
choice_shares <- DT[, .(Share = sum(choice) / N), by = j][order(j)]
print(choice_shares)

# choicer ----------------------------------------------------------------------

dt_inputs <- prepare_mnl_data(
  data = DT,
  id_col = "id",
  alt_col = "j",
  choice_col = "choice",
  covariate_cols = c("X", "W")
)

alt_mapping <- dt_inputs$alt_mapping |> copy()

alt_mapping[, nest := fcase(
  j %in% c(1,2), 1L,
  j %in% c(3,4,5), 2L
)]

dt_inputs$nest_idx <- alt_mapping$nest

nloptr_opts <- list(
  "algorithm" = "NLOPT_LD_LBFGS",
  "xtol_rel" = 1.0e-8,
  "maxeval" = 1e+3,
  "print_level" = 0L,
  "check_derivatives" = TRUE,
  "check_derivatives_print" = 'all'
)

param_names <- c("X","W","Lambda_1","Lambda_2", paste0("ASC_", alt_mapping[alt_int > 1]$j))

nestlogit_result <- run_nestlogit(
  input_data = dt_inputs,
  use_asc = TRUE,
  param_names = param_names,
  nloptr_opts = nloptr_opts
)


# mlogit------------------------------------------------------------------------

DT[, alt_factor := factor(j)]
DT_idx <- dfidx(as.data.frame(DT), idx = c("id", "alt_factor"))

nests_spec <- list(
  # nest0 = "0",
  nest1 = c("1","2"),
  nest2 = c("3","4","5")
)

nl <- mlogit(
  choice ~ X + W,
  data       = DT_idx,
  nests      = nests_spec,
  un.nest.el = FALSE
)

summary(nl) |> print()
