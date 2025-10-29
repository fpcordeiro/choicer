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
  "1" = 0.5,
  "2" = 0.3,
  "3" = -0.2,
  "4" = -0.5,
  "5" = 0.4
)

# Nest structure and dissimilarity parameters (lambdas)
lambda_0 <- 1.0 # singleton nest for outside option
lambda_1 <- 0.8 # Nest 1: {1, 2}
lambda_2 <- 0.2 # Nest 2: {3, 4, 5}

# Create Data Structure & Systematic Utility -----------------------------------

# Create the main data.table
# We have 6 alternatives: j=0 (outside) and j=1...5 (inside)
DT <- CJ(id = 1:N, j = 0:5)

# Define the nests
# Nest 1: {1, 2}
# Nest 2: {3, 4, 5}
DT[, nest := fcase(
  j == 0, 0,
  j %in% c(1, 2), 1,
  j %in% c(3, 4, 5), 2
)]

# Assign lambdas
lambda_vec_all <- c(lambda_0, lambda_1, lambda_2)
DT[, lambda := lambda_vec_all[nest + 1]]

# Add deltas (ASCs)
delta_vec_all <- c(0, deltas)
DT[, delta := delta_vec_all[j + 1]]

# Simulate covariates
# Covariates are 0 for the outside option
DT[, X := 0][j > 0, X := rnorm(N * 5, mean = 1, sd = 1)]
DT[, W := 0][j > 0, W := runif(N * 5, min = -2, max = 2)]

# Calculate systematic utility (V_ij)
DT[, V := delta + beta_X * X + beta_W * W]

# Calculate choice probabilities -----------------------------------------------

lse <- function(x, na.rm = FALSE) {
  # Empty input: log( sum(exp(numeric(0))) ) = log(0) = -Inf
  if (length(x) == 0L) return(-Inf)

  # Handle NA / NaN per na.rm
  if (!na.rm && anyNA(x)) return(NA_real_)
  if (na.rm) x <- x[!is.na(x)]

  a <- max(x)

  # If any +Inf is present, result is +Inf regardless of NA elsewhere
  if (a == Inf) return(Inf)

  # If only -Inf are include, return -Inf
  if (a == -Inf) return(-Inf)

  a + log(sum(exp(x - a)))
}

# Calculate conditional probability: P(j|k)
DT[, V_over_lambda := V / lambda]
DT_iv <- DT[, .(IV = lse(V_over_lambda)), by = .(id, nest, lambda)]
DT_iv <- unique(DT_iv, by = c("id","nest"))

# Calculate marginal probability: P(k)
DT_iv[, log_denom := lse(lambda * IV), by = id]
DT_iv[, nest_prob := exp(lambda * IV - log_denom)]

# Join probabilities back
DT[
  DT_iv[, .(id, nest, nest_prob, IV)],
  on = c("id", "nest"),
  `:=`(nest_prob = i.nest_prob,
       IV = i.IV)
]

# Calculate full probability: P(j) = P(j|k) * P(k)
DT[, cond_prob := exp(V_over_lambda - IV)]
DT[, P_j := cond_prob * nest_prob]

# For j=0 (nest 0), P(j=0|k=0) = 1, so P(j=0) = P(k=0)
DT[j == 0, P_j := nest_prob]

# Calculate Final Utility and Simulate Choices ---------------------------------

# Find the chosen alternative for each individual
DT[, choice:= fifelse(j == sample(j, size=1, prob=P_j), 1L, 0L), by = id]

# Calculate the average probability for each alternative
choice_shares <- DT[, .(Share = sum(choice) / N), by = j][order(j)]
avg_probs <- DT[, .(Mean_Prob = mean(P_j)), by = j][order(j)]
cat("\n--- Verification ---\n")
setnames(choice_shares, "j", "Alternative")
verification <- merge(avg_probs, choice_shares, by.x = "j", by.y = "Alternative")
setnames(verification, "j", "Alternative")
print(verification, digits = 4)

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
  j == 0, 1L,
  j %in% c(1,2), 2L,
  j %in% c(3,4,5), 3L
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

# Derivatives -----------------------------------------------------------------

library(numDeriv)

J <- nrow(dt_inputs$alt_mapping)
K_x <- ncol(dt_inputs$X)
K_l <- sum(table(dt_inputs$nest_idx) > 1)
x0 <- c(rep(0, K_x), rep(0.5, K_l),  rep(0, J - 1))

out_choicer <- nl_loglik_gradient_parallel(
  theta = x0,
  X = dt_inputs$X,
  alt_idx = dt_inputs$alt_idx,
  choice_idx = dt_inputs$choice_idx,
  nest_idx = dt_inputs$nest_idx,
  M = dt_inputs$M,
  weights = dt_inputs$weights,
  use_asc = TRUE,
  include_outside_option = dt_inputs$include_outside_option
)

wrapper_fun <- function(x) {
  output <- nl_loglik_gradient_parallel(
    theta = x,
    X = dt_inputs$X,
    alt_idx = dt_inputs$alt_idx,
    choice_idx = dt_inputs$choice_idx,
    nest_idx = dt_inputs$nest_idx,
    M = dt_inputs$M,
    weights = dt_inputs$weights,
    use_asc = TRUE,
    include_outside_option = dt_inputs$include_outside_option
  )
  return(output$objective)
}
grad_num <- numDeriv::grad(
  func = wrapper_fun,
  x = x0,
  method = "Richardson"
)

grad_num
drop(out_choicer$gradient)
abs(drop(out_choicer$gradient) - grad_num)
