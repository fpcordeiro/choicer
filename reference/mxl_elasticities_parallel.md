# Compute aggregate elasticities for mixed logit model

Computes the aggregate elasticity matrix (weighted average of individual
elasticities) for the Mixed Logit model. The elasticity E(i,j)
represents the percentage change in the probability of choosing
alternative i when the attribute of alternative j changes by 1%.

## Usage

``` r
mxl_elasticities_parallel(
  theta,
  X,
  W,
  alt_idx,
  choice_idx,
  M,
  weights,
  eta_draws,
  rc_dist,
  elast_var_idx,
  is_random_coef,
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

  parameter vector (beta, \[mu\], L, delta)

- X:

  design matrix for fixed coefficients; sum(M_i) x K_x

- W:

  design matrix for random coefficients; sum(M_i) x K_w or J x K_w

- alt_idx:

  sum(M) x 1 vector with indices of alternatives; 1-based indexing

- choice_idx:

  N x 1 vector (kept for API consistency, not used)

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- eta_draws:

  Array with draws; K_w x S x N

- rc_dist:

  K_w vector indicating distribution (0=normal, 1=log-normal)

- elast_var_idx:

  1-based index of the variable for elasticity computation

- is_random_coef:

  TRUE if variable is in W (random coef), FALSE if in X (fixed coef)

- rc_correlation:

  whether random coefficients are correlated

- rc_mean:

  whether mu parameters are estimated

- use_asc:

  whether ASCs are included

- include_outside_option:

  whether outside option is included

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

J x J matrix of aggregate elasticities

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
elas <- mxl_elasticities_parallel(coef(fit), d$X, d$W, d$alt_idx,
  d$choice_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
  elast_var_idx = 1L, is_random_coef = FALSE,
  rc_correlation = FALSE, rc_mean = FALSE)
elas
#>              [,1]         [,2]         [,3]
#> [1,]  0.004514083 -0.002819572  0.001085349
#> [2,] -0.004301106 -0.002986786  0.000274385
#> [3,] -0.003179751 -0.003194349 -0.007922687
# }
```
