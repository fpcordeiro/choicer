# BHHH/OPG information matrix for multinomial logit model

Computes the weighted outer product of per-individual scores \\\sum_i
w_i\\ s_i s_i^\top\\ for the Multinomial Logit model. The per-individual
score \\s_i\\ is the (positive) gradient of individual \\i\\'s
log-likelihood contribution and is weight-free; the supplied `weights`
enter only as the leading multiplier. Passing `weights = w` yields the
ordinary weighted BHHH/OPG information; passing `weights = w^2` yields
the sandwich *meat* \\B = \sum_i w_i^2 s_i s_i^\top\\ used for robust
(WESML) inference.

## Usage

``` r
mnl_bhhh_parallel(
  theta,
  X,
  alt_idx,
  choice_idx,
  M,
  weights,
  use_asc = TRUE,
  include_outside_option = FALSE
)
```

## Arguments

- theta:

  K + J - 1 or K + J vector with model parameters

- X:

  sum(M) x K design matrix with covariates. Stacks M\[i\] x K matrices
  for individual i.

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

- use_asc:

  whether to use alternative-specific constants

- include_outside_option:

  whether to include outside option normalized to 0 (if so, the outside
  option is not included in the data)

## Value

A symmetric positive-semidefinite information matrix \\\sum_i w_i\\ s_i
s_i^\top\\ (same sign convention as the negated Hessian).

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
B <- mnl_bhhh_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
  fit$data$choice_idx, fit$data$M, fit$data$weights)
dim(B)
#> [1] 4 4
# }
```
