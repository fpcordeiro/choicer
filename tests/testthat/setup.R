# Setup file for testthat - loaded before all tests
# Contains fixtures, tolerances, and helper data generators

library(data.table)

# Helper to create nested logit inputs from prepare_nl_data() (used in multiple test files)
# Creates 6 explicit alternatives in 3 nests:
#   Nest "A": j=0 (singleton, zero covariates — acts as outside option)
#   Nest "B": j=1,2
#   Nest "C": j=3,4,5
create_nl_inputs <- function(seed = 123) {
  set.seed(seed)
  N <- 30
  J <- 6

  dt <- data.table(
    id = rep(1:N, each = J),
    j = rep(0:(J-1), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  # "Outside option" (j=0) has zero covariates
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  # Add nest column
  dt[, nest := fifelse(j == 0, "A", fifelse(j <= 2, "B", "C"))]

  prepare_nl_data(
    dt, "id", "j", "choice", c("x1", "x2"),
    nest_col = "nest"
  )
}

# Numerical tolerances for comparisons
TOL_GRAD <- 1e-5       # Gradient vs numerical derivative
TOL_HESS <- 1e-4       # Hessian vs numerical Hessian
TOL_LOGLIK <- 1e-10    # Likelihood value precision
TOL_PROB <- 1e-8       # Probability sum-to-one

# Create a small deterministic MNL dataset (N=20, J=3, K=2)
create_small_mnl_data <- function(seed = 12345) {
  set.seed(seed)
  N <- 20
  J <- 3
  dt <- data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[, choice := 0L]
  dt[, choice := {
    probs <- rep(1/J, J)
    sample(c(1L, rep(0L, J - 1)))
  }, by = id]
  dt[]
}

# Create a small deterministic MXL dataset (N=30, J=3, K_x=1, K_w=2)
create_small_mxl_data <- function(seed = 42) {
  set.seed(seed)
  N <- 30
  J <- 3
  dt <- data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1 = rnorm(N * J),
    w1 = rnorm(N * J),
    w2 = runif(N * J, -1, 1)
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
  dt[]
}

# Create a small nested logit dataset (N=30, J=6, 2 nests)
create_small_nl_data <- function(seed = 123) {
  set.seed(seed)
  N <- 30
  J <- 6
  dt <- data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    nest = rep(c(1L, 1L, 1L, 2L, 2L, 2L), N),
    x1 = rnorm(N * J),
    x2 = runif(N * J, -1, 1)
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
  dt[]
}

# Create an identified MXL dataset with real signal (N=200, J=3, K_x=1, K_w=1).
# Unlike `create_small_mxl_data()` (uniform-random choices, flat likelihood),
# this fixture places the MLE at an interior point so standard-error estimators
# can be meaningfully compared. Delegates to `simulate_mxl_data()` so the
# package-level DGP is the single source of truth.
create_identified_mxl_data <- function(seed = 101, N = 200, J = 3) {
  sim <- simulate_mxl_data(
    N              = N,
    J              = J,
    beta           = 1.0,
    delta          = rep(0, J),
    Sigma          = matrix(0.25, 1, 1),
    seed           = seed,
    outside_option = FALSE,
    vary_choice_set = FALSE
  )
  sim$data
}

# Create dataset with known probabilities (equal utilities -> equal probs)
create_equal_prob_data <- function(N = 100, J = 3, seed = 999) {
  set.seed(seed)
  dt <- data.table(
    id = rep(1:N, each = J),
    alt = rep(1:J, N),
    x1 = 0  # Zero covariate -> no effect
  )
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  dt[]
}

# Create an identified MNP dataset with real signal: probit DGP with
# correlated differenced errors. `beta` is on the identified (sigma_11 = 1)
# scale when Sigma[1, 1] = 1. Returns the long-format data plus the true
# identified parameters.
create_mnp_sim_data <- function(N = 800, J = 3,
                                beta = c(1.0, -0.5),
                                Sigma = matrix(c(1.0, 0.5, 0.5, 1.5), 2, 2),
                                seed = 99) {
  stopifnot(J >= 2, nrow(Sigma) == J - 1)
  set.seed(seed)
  p <- J - 1
  K <- length(beta)
  cov_cols <- paste0("x", seq_len(K))

  dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
  for (k in seq_len(K)) dt[, (cov_cols[k]) := rnorm(.N)]

  V <- matrix(as.matrix(dt[, ..cov_cols]) %*% beta, nrow = J)   # J x N
  L <- t(chol(Sigma))
  eps <- L %*% matrix(rnorm(p * N), p, N)
  W <- V[2:J, , drop = FALSE] - rep(1, p) %o% V[1, ] + eps      # p x N

  chosen <- apply(W, 2, function(w) if (max(w) < 0) 1L else which.max(w) + 1L)
  dt[, choice := as.integer(alt == chosen[id])]

  list(
    data = dt[],
    beta_id = beta / sqrt(Sigma[1, 1]),
    Sigma_id = Sigma / Sigma[1, 1],
    covariate_cols = cov_cols
  )
}
