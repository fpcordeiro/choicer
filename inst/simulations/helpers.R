# helpers.R - Shared DGP functions for simulation scripts
# Source this file before running any simulation script.

library(data.table)
library(choicer)

# Utilities ====================================================================

#' Log-sum-exp (numerically stable)
lse <- function(x) {
  a <- max(x)
  if (a == Inf) return(Inf)
  if (a == -Inf) return(-Inf)
  a + log(sum(exp(x - a)))
}

#' Print parameter comparison table
compare_params <- function(true, estimated, names = NULL) {
  if (is.null(names)) names <- paste0("param_", seq_along(true))
  tbl <- data.table(
    Parameter  = names,
    True       = round(true, 4),
    Estimated  = round(estimated, 4),
    Difference = round(estimated - true, 4)
  )
  print(tbl)
  cat("Max |difference|:", round(max(abs(estimated - true)), 4), "\n")
  invisible(tbl)
}

# MNL DGP =====================================================================

#' Simulate multinomial logit data with varying choice sets and outside option
#'
#' @param N Number of choice situations.
#' @param J Number of inside alternatives.
#' @param beta True fixed coefficients (length K_x).
#' @param delta True ASCs for inside alternatives (length J).
#' @param seed Random seed.
#' @return List with `data`, `true_params`, and `settings`.
simulate_mnl_data <- function(
    N     = 5000,
    J     = 5,
    beta  = c(0.8, -0.6),
    delta = NULL,
    seed  = 123
) {
  set.seed(seed)
  K_x <- length(beta)
  if (is.null(delta)) delta <- rep(c(0.5, -0.5), length.out = J)

  # Varying choice set size (2 to J inside alternatives per individual)
  choice_set_sizes <- sample(2:J, N, replace = TRUE)
  alt_indices <- vector("list", N)
  for (i in seq_len(N)) {
    alt_indices[[i]] <- sort(sample(J, choice_set_sizes[i]))
  }

  # Inside alternatives
  dt <- data.table(id = rep(seq_len(N), times = choice_set_sizes))
  dt[, `:=`(x1 = runif(.N, -1, 1), x2 = runif(.N, -1, 1))]
  dt[, `:=`(
    alt       = alt_indices[[id]],
    delta_val = delta[alt_indices[[id]]]
  ), by = id]

  # Outside option (alt = 0, all covariates and delta = 0)
  dt_outside <- dt[, .(alt = 0L), keyby = id]
  dt <- rbindlist(list(dt_outside, dt), use.names = TRUE, fill = TRUE)
  fill_cols <- setdiff(names(dt), c("id", "alt"))
  dt[alt == 0L, (fill_cols) := 0]
  setkey(dt, id, alt)

  # Gumbel errors and argmax choice
  dt[, epsilon := -log(-log(runif(.N)))]
  dt[, utility := delta_val + x1 * beta[1] + x2 * beta[2] + epsilon]
  dt[, choice := fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]
  dt[, c("delta_val", "epsilon", "utility") := NULL]

  list(
    data = dt,
    true_params = list(beta = beta, delta = delta),
    settings = list(N = N, J = J, K_x = K_x)
  )
}

# MXL DGP =====================================================================

#' Simulate mixed logit data with correlated random coefficients
#'
#' @param N Number of choice situations.
#' @param J Number of inside alternatives.
#' @param beta True fixed coefficients for X variables.
#' @param delta True ASCs for inside alternatives.
#' @param Sigma True variance-covariance matrix of random coefficients.
#' @param seed Random seed.
#' @return List with `data`, `true_params`, and `settings`.
simulate_mxl_data <- function(
    N     = 5000,
    J     = 4,
    beta  = c(0.8, -0.6),
    delta = NULL,
    Sigma = matrix(c(1.0, 0.5, 0.5, 1.5), nrow = 2),
    seed  = 123
) {
  set.seed(seed)
  K_x <- length(beta)
  K_w <- ncol(Sigma)
  if (is.null(delta)) delta <- rep(c(0.5, -0.5), length.out = J)

  # Cholesky decomposition and individual random coefficients
  L_true <- t(chol(Sigma))
  gamma_i <- t(L_true %*% matrix(rnorm(N * K_w), nrow = K_w, ncol = N))

  # Cholesky parameterization (log-diagonal, free off-diagonal)
  L_params <- c(log(L_true[1, 1]), L_true[2, 1], log(L_true[2, 2]))

  # Varying choice sets
  choice_set_sizes <- sample(2:J, N, replace = TRUE)
  alt_indices <- vector("list", N)
  for (i in seq_len(N)) {
    alt_indices[[i]] <- sort(sample(J, choice_set_sizes[i]))
  }

  # Inside alternatives
  dt <- data.table(id = rep(seq_len(N), times = choice_set_sizes))
  dt[, `:=`(
    w1 = runif(.N, -1, 1), w2 = runif(.N, -1, 1),
    x1 = runif(.N, -1, 1), x2 = runif(.N, -1, 1)
  )]
  dt[, `:=`(gamma1 = gamma_i[id, 1], gamma2 = gamma_i[id, 2])]
  dt[, `:=`(
    alt       = alt_indices[[id]],
    delta_val = delta[alt_indices[[id]]]
  ), by = id]

  # Outside option
  dt_outside <- dt[, .(alt = 0L), keyby = id]
  dt <- rbindlist(list(dt_outside, dt), use.names = TRUE, fill = TRUE)
  fill_cols <- setdiff(names(dt), c("id", "alt"))
  dt[alt == 0L, (fill_cols) := 0]
  setkey(dt, id, alt)

  # Gumbel errors and argmax choice
  dt[, epsilon := -log(-log(runif(.N)))]
  dt[, utility := delta_val + x1 * beta[1] + x2 * beta[2] +
       w1 * gamma1 + w2 * gamma2 + epsilon]
  dt[, choice := fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]
  dt[, c("delta_val", "gamma1", "gamma2", "epsilon", "utility") := NULL]

  list(
    data = dt,
    true_params = list(beta = beta, delta = delta, Sigma = Sigma, L_params = L_params),
    settings = list(N = N, J = J, K_x = K_x, K_w = K_w)
  )
}

# NL DGP ======================================================================

#' Simulate nested logit data with outside option
#'
#' Computes nested logit probabilities analytically and samples choices.
#' The outside option (j = 0) is in a singleton nest with lambda = 1.
#'
#' @param N Number of choice situations.
#' @param beta True coefficients (length 2, for X and W).
#' @param delta True ASCs for inside alternatives (named vector).
#' @param nests List of integer vectors defining nest membership for inside alts.
#' @param lambdas Dissimilarity parameters for each nest (length = length(nests)).
#' @param seed Random seed.
#' @return List with `data`, `true_params`, `nest_structure`, and `settings`.
simulate_nl_data <- function(
    N       = 10000,
    beta    = c(1.5, -0.8),
    delta   = c("1" = 0.5, "2" = 0.3, "3" = -0.2, "4" = -0.5, "5" = 0.4),
    nests   = list(c(1, 2), c(3, 4, 5)),
    lambdas = c(0.8, 0.2),
    seed    = 123
) {
  set.seed(seed)
  J_inside <- length(delta)

  # Full balanced panel: j = 0 (outside) and j = 1..J_inside
  dt <- CJ(id = 1:N, j = 0:J_inside)

  # Nest assignments (outside option = nest 0, lambda = 1)
  dt[, nest := 0L]
  for (g in seq_along(nests)) {
    dt[j %in% nests[[g]], nest := as.integer(g)]
  }
  lambda_all <- c(1.0, lambdas)
  dt[, lambda := lambda_all[nest + 1]]

  # ASCs (outside option = 0)
  delta_all <- c(0, delta)
  dt[, delta_val := delta_all[j + 1]]

  # Covariates (zero for outside option)
  dt[, X := 0.0][j > 0, X := rnorm(.N, mean = 1, sd = 1)]
  dt[, W := 0.0][j > 0, W := runif(.N, -2, 2)]

  # Systematic utility
  dt[, V := delta_val + beta[1] * X + beta[2] * W]

  # --- Nested logit probabilities ---
  # Inclusive values per nest
  dt[, V_over_lambda := V / lambda]
  dt_iv <- dt[, .(IV = lse(V_over_lambda)), by = .(id, nest, lambda)]
  dt_iv <- unique(dt_iv, by = c("id", "nest"))

  # Nest selection probabilities
  dt_iv[, log_denom := lse(lambda * IV), by = id]
  dt_iv[, nest_prob := exp(lambda * IV - log_denom)]

  # Join back and compute full probabilities
  dt[dt_iv[, .(id, nest, nest_prob, IV)],
     on = c("id", "nest"),
     `:=`(nest_prob = i.nest_prob, IV = i.IV)]
  dt[, cond_prob := exp(V_over_lambda - IV)]
  dt[, P_j := cond_prob * nest_prob]
  dt[j == 0, P_j := nest_prob]  # singleton nest: P(j|k) = 1

  # Sample choices from probabilities
  dt[, choice := fifelse(j == sample(j, size = 1, prob = P_j), 1L, 0L), by = id]

  # Clean up
  dt[, c("nest", "lambda", "delta_val", "V", "V_over_lambda",
         "IV", "nest_prob", "cond_prob", "P_j") := NULL]

  list(
    data = dt,
    true_params = list(beta = beta, delta = delta, lambdas = lambdas),
    nest_structure = nests,
    settings = list(N = N, J_inside = J_inside)
  )
}
