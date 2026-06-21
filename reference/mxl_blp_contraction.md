# BLP contraction mapping for mixed logit

Finds the ASC (delta) parameters such that predicted market shares match
target shares, using the contraction mapping of Berry, Levinsohn, and
Pakes (1995).

## Usage

``` r
mxl_blp_contraction(
  delta,
  target_shares,
  X,
  W,
  beta,
  mu,
  L_params,
  alt_idx,
  M,
  weights,
  eta_draws,
  rc_dist,
  rc_correlation = TRUE,
  rc_mean = FALSE,
  include_outside_option = FALSE,
  tol = 1e-08,
  max_iter = 1000L
)
```

## Arguments

- delta:

  J-1 or J vector with initial guess for deltas (ASCs)

- target_shares:

  J vector with target market shares

- X:

  design matrix for fixed coefficients; sum(M_i) x K_x

- W:

  design matrix for random coefficients; sum(M_i) x K_w or J x K_w

- beta:

  K_x vector with fixed coefficients

- mu:

  K_w vector with mean parameters (raw, will be transformed if
  log-normal)

- L_params:

  Cholesky parameters vector

- alt_idx:

  sum(M) x 1 vector with indices of alternatives; 1-based indexing

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- eta_draws:

  Array with draws; K_w x S x N

- rc_dist:

  K_w vector indicating distribution (0=normal, 1=log-normal)

- rc_correlation:

  whether random coefficients are correlated

- rc_mean:

  whether mu parameters represent means (TRUE) or are zero (FALSE)

- include_outside_option:

  whether outside option is included

- tol:

  convergence tolerance (default 1e-8)

- max_iter:

  maximum iterations (default 1000)

## Value

vector with converged delta (ASC) values

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
fit <- run_mxlogit(input_data = d, eta_draws = eta)
#> Optimization run time 0h:0m:0.01s
pm <- fit$param_map
delta <- mxl_blp_contraction(rep(0, J), rep(1/J, J), d$X, d$W,
  coef(fit)[pm$beta], rep(0, ncol(d$W)), coef(fit)[pm$sigma],
  d$alt_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
  rc_correlation = FALSE, rc_mean = FALSE)
delta
#>             [,1]
#> [1,]  0.00000000
#> [2,] -0.01837654
#> [3,] -0.04093111
# }
```
