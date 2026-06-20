# Analytical Hessian of the log-likelihood v2

Computes the Hessian of the log-likelihood for the Mixed Logit model
using OpenMP for parallelization. Mirrors the parameters of
mxl_loglik_gradient_parallel.

## Usage

``` r
mxl_hessian_parallel(
  theta,
  X,
  W,
  alt_idx,
  choice_idx,
  M,
  weights,
  eta_draws,
  rc_dist,
  rc_correlation = TRUE,
  rc_mean = FALSE,
  use_asc = TRUE,
  include_outside_option = FALSE,
  gen_seed = -1L,
  gen_scramble = 1L,
  gen_S = 0L
)
```

## Arguments

- theta:

  vector collecting model parameters (beta, mu, L, delta (ASCs))

- X:

  design matrix for covariates with fixed coefficients; sum(M_i) x K_x

- W:

  design matrix for covariates with random coefficients; sum(M_i) x K_w
  or J x K_w

- alt_idx:

  sum(M) x 1 vector with indices of alternatives within each choice set;
  1-based indexing

- choice_idx:

  N x 1 vector with indices of chosen alternatives; 1-based indexing
  relative to X; 0 is used if include_outside_option=True

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- eta_draws:

  Array with choice situation draws; K_w x S x N

- rc_dist:

  K_w x 1 integer vector indicating distribution of random coefficients:
  0 = normal, 1 = log-normal

- rc_correlation:

  whether random coefficients should be correlated

- rc_mean:

  whether to estimate means for random coefficients.

- use_asc:

  whether to use alternative-specific constants.

- include_outside_option:

  whether to include outside option normalized to 0 (if so, the outside
  option is not included in the data)

- gen_seed:

  Integer master seed for the on-the-fly Halton generator. `< 0`
  (default) uses the materialized `eta_draws` cube; `>= 0` generates
  draws on the fly from this seed.

- gen_scramble:

  Integer scramble mode for on-the-fly generation: `0` = identity
  permutations (plain Halton, compat), `1` = Owen (2017) digit
  scrambling.

- gen_S:

  Integer number of draws per individual, used only when
  `gen_seed >= 0`.

## Value

Hessian evaluated at input arguments

## Note

For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE, the
distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k \*
eta), where exp(mu_k) shifts the location and exp(L_k \* eta) ~
LogNormal(0, sigma_k^2). This differs from the textbook parameterization
exp(mu_k + L_k \* eta).

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#>         id   alt         x1          w1
#>      <int> <int>      <num>       <num>
#>   1:     1     1  1.3709584 -0.04069848
#>   2:     1     2 -0.5646982 -1.55154482
#>   3:     1     3  0.3631284  1.16716955
#>   4:     2     1  0.6328626 -0.27364570
#>   5:     2     2  0.4042683 -0.46784532
#>  ---                                   
#> 146:    49     2  1.1133860 -0.47733551
#> 147:    49     3 -0.4809928 -0.16626149
#> 148:    50     1 -0.4331690  0.86256338
#> 149:    50     2  0.6968626  0.09734049
#> 150:    50     3 -1.0563684 -1.62561674
dt[, choice := 0L]
#>         id   alt         x1          w1 choice
#>      <int> <int>      <num>       <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848      0
#>   2:     1     2 -0.5646982 -1.55154482      0
#>   3:     1     3  0.3631284  1.16716955      0
#>   4:     2     1  0.6328626 -0.27364570      0
#>   5:     2     2  0.4042683 -0.46784532      0
#>  ---                                          
#> 146:    49     2  1.1133860 -0.47733551      0
#> 147:    49     3 -0.4809928 -0.16626149      0
#> 148:    50     1 -0.4331690  0.86256338      0
#> 149:    50     2  0.6968626  0.09734049      0
#> 150:    50     3 -1.0563684 -1.62561674      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1          w1 choice
#>      <int> <int>      <num>       <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848      0
#>   2:     1     2 -0.5646982 -1.55154482      0
#>   3:     1     3  0.3631284  1.16716955      1
#>   4:     2     1  0.6328626 -0.27364570      0
#>   5:     2     2  0.4042683 -0.46784532      0
#>  ---                                          
#> 146:    49     2  1.1133860 -0.47733551      0
#> 147:    49     3 -0.4809928 -0.16626149      0
#> 148:    50     1 -0.4331690  0.86256338      1
#> 149:    50     2  0.6968626  0.09734049      0
#> 150:    50     3 -1.0563684 -1.62561674      0
d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
eta <- get_halton_normals(50, d$N, ncol(d$W))
theta <- rep(0, ncol(d$X) + ncol(d$W) + nrow(d$alt_mapping) - 1)
H <- mxl_hessian_parallel(theta, d$X, d$W, d$alt_idx, d$choice_idx,
  d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
  rc_correlation = FALSE, rc_mean = FALSE)
dim(H)
#> [1] 4 4
# }
```
