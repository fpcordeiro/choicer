# BLP95 contraction mapping to find delta given target shares

BLP95 contraction mapping to find delta given target shares

## Usage

``` r
blp_contraction(
  delta,
  target_shares,
  X,
  beta,
  alt_idx,
  M,
  weights,
  include_outside_option = FALSE,
  tol = 1e-08,
  max_iter = 1000L
)
```

## Arguments

- delta:

  J x 1 vector with initial guess for deltas (ASCs)

- target_shares:

  J x 1 vector with target shares for each alternative

- X:

  sum(M) x K design matrix with covariates. M\[i\] x K matrix for
  individual i

- beta:

  K x 1 vector with model parameters

- alt_idx:

  sum(M) x 1 vector with indices of alternatives within each choice set;
  1-based indexing

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- include_outside_option:

  whether to include outside option normalized to 0 (if so, the outside
  option is not included in the data)

- tol:

  convergence tolerance

- max_iter:

  maximum number of iterations

## Value

vector with contraction's delta (ASCs) output

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1          x2
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
#>         id   alt         x1          x2 choice
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
#>         id   alt         x1          x2 choice
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
fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s
beta <- coef(fit)[fit$param_map$beta]
delta <- blp_contraction(rep(0, J), rep(1/J, J), fit$data$X,
  beta, fit$data$alt_idx, fit$data$M, fit$data$weights)
delta
#>             [,1]
#> [1,]  0.00000000
#> [2,] -0.03013133
#> [3,]  0.01204404
# }
```
