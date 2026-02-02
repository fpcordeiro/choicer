# Setup file for testthat - loaded before all tests
# Contains fixtures, tolerances, and helper data generators

library(data.table)

# Helper to create nested logit inputs from MNL data (used in multiple test files)
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
  # Outside option (j=0) has zero covariates
  dt[j == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]

  # Prepare using MNL data prep
  inputs <- prepare_mnl_data(
    dt, "id", "j", "choice", c("x1", "x2"),
    outside_opt_label = 0L,
    include_outside_option = TRUE
  )

  # Add nest structure:
  # Nest 1: j=0 (singleton - outside option)
  # Nest 2: j=1,2
  # Nest 3: j=3,4,5
  inputs$nest_idx <- c(1L, 2L, 2L, 3L, 3L, 3L)

  inputs
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
