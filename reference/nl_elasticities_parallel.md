# Compute aggregate elasticities for the Nested Logit model

Computes the aggregate (weighted-average) elasticity matrix for the
Nested Logit model. Reduces to the MNL elasticities when all lambda = 1.

## Usage

``` r
nl_elasticities_parallel(
  theta,
  X,
  alt_idx,
  choice_idx,
  nest_idx,
  M,
  weights,
  elast_var_idx,
  use_asc = TRUE,
  include_outside_option = FALSE
)
```

## Arguments

- theta:

  (K + n_non_singleton_nests + n_delta) vector with model parameters.
  Order: `[beta (K), lambda (non-singleton), delta]`.

- X:

  sum(M) x K design matrix with covariates.

- alt_idx:

  sum(M) x 1 vector with indices of alternatives; 1-based indexing.

- choice_idx:

  N x 1 vector (kept for API consistency, not used).

- nest_idx:

  J x 1 vector with nest indices for each alternative; 1-based indexing.

- M:

  N x 1 vector with number of alternatives for each individual.

- weights:

  N x 1 vector with weights for each observation.

- elast_var_idx:

  1-based index of the column in X for which to compute the elasticity.

- use_asc:

  whether to use alternative-specific constants.

- include_outside_option:

  whether to include outside option normalized to V=0, lambda=1.

## Value

J x J matrix of aggregate elasticities (row = responding alt, col =
perturbed alt).

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 4
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1         x2
#>      <int> <int>      <num>      <num>
#>   1:     1     1  1.3709584 -2.0009292
#>   2:     1     2 -0.5646982  0.3337772
#>   3:     1     3  0.3631284  1.1713251
#>   4:     1     4  0.6328626  2.0595392
#>   5:     2     1  0.4042683 -1.3768616
#>  ---                                  
#> 196:    49     4  1.0857749  1.0965134
#> 197:    50     1  0.4037749  0.4420131
#> 198:    50     2  0.5864875  0.2410163
#> 199:    50     3  1.8152284 -0.2556077
#> 200:    50     4  0.1288214  0.9310329
dt[, nest := ifelse(alt <= 2, "A", "B")]
#>         id   alt         x1         x2   nest
#>      <int> <int>      <num>      <num> <char>
#>   1:     1     1  1.3709584 -2.0009292      A
#>   2:     1     2 -0.5646982  0.3337772      A
#>   3:     1     3  0.3631284  1.1713251      B
#>   4:     1     4  0.6328626  2.0595392      B
#>   5:     2     1  0.4042683 -1.3768616      A
#>  ---                                         
#> 196:    49     4  1.0857749  1.0965134      B
#> 197:    50     1  0.4037749  0.4420131      A
#> 198:    50     2  0.5864875  0.2410163      A
#> 199:    50     3  1.8152284 -0.2556077      B
#> 200:    50     4  0.1288214  0.9310329      B
dt[, choice := 0L]
#>         id   alt         x1         x2   nest choice
#>      <int> <int>      <num>      <num> <char>  <int>
#>   1:     1     1  1.3709584 -2.0009292      A      0
#>   2:     1     2 -0.5646982  0.3337772      A      0
#>   3:     1     3  0.3631284  1.1713251      B      0
#>   4:     1     4  0.6328626  2.0595392      B      0
#>   5:     2     1  0.4042683 -1.3768616      A      0
#>  ---                                                
#> 196:    49     4  1.0857749  1.0965134      B      0
#> 197:    50     1  0.4037749  0.4420131      A      0
#> 198:    50     2  0.5864875  0.2410163      A      0
#> 199:    50     3  1.8152284 -0.2556077      B      0
#> 200:    50     4  0.1288214  0.9310329      B      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1         x2   nest choice
#>      <int> <int>      <num>      <num> <char>  <int>
#>   1:     1     1  1.3709584 -2.0009292      A      0
#>   2:     1     2 -0.5646982  0.3337772      A      0
#>   3:     1     3  0.3631284  1.1713251      B      0
#>   4:     1     4  0.6328626  2.0595392      B      1
#>   5:     2     1  0.4042683 -1.3768616      A      0
#>  ---                                                
#> 196:    49     4  1.0857749  1.0965134      B      1
#> 197:    50     1  0.4037749  0.4420131      A      0
#> 198:    50     2  0.5864875  0.2410163      A      0
#> 199:    50     3  1.8152284 -0.2556077      B      0
#> 200:    50     4  0.1288214  0.9310329      B      1
fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#> Optimization run time 0h:0m:0.01s
elas <- nl_elasticities_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
  fit$data$choice_idx, fit$data$nest_idx, fit$data$M, fit$data$weights,
  elast_var_idx = 1L)
elas
#>              [,1]         [,2]         [,3]         [,4]
#> [1,]  0.001361410 -0.002540128  0.001811961 -0.013019968
#> [2,]  0.001254277 -0.002616362  0.001811961 -0.013019968
#> [3,] -0.006952266  0.003469046 -0.043995957 -0.001679743
#> [4,] -0.006952266  0.003469046 -0.015185111  0.005494552
# }
```
