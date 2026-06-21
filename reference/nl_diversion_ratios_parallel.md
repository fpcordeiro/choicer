# Compute Nested Logit diversion ratios (parallelized over individuals)

Computes the diversion ratio matrix DR(j-\>k) for the Nested Logit
model. Entry (k, j) = fraction of demand lost by alternative j captured
by k. Reduces to the MNL diversion ratios when all lambda = 1.

## Usage

``` r
nl_diversion_ratios_parallel(
  theta,
  X,
  alt_idx,
  nest_idx,
  M,
  weights,
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

- nest_idx:

  J x 1 vector with nest indices for each alternative; 1-based indexing.

- M:

  N x 1 vector with number of alternatives for each individual.

- weights:

  N x 1 vector with weights for each observation.

- use_asc:

  whether to use alternative-specific constants.

- include_outside_option:

  whether to include outside option normalized to V=0, lambda=1.

## Value

J x J matrix where entry (k, j) = DR(j-\>k). Diagonal is 0.

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
dr <- nl_diversion_ratios_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
  fit$data$nest_idx, fit$data$M, fit$data$weights)
dr
#>            [,1]       [,2]       [,3]      [,4]
#> [1,]  0.0000000 -0.8981740 0.45017345 0.4192658
#> [2,] -1.0621318  0.0000000 0.49007051 0.4562837
#> [3,]  1.4249183  1.3117489 0.00000000 0.1244504
#> [4,]  0.6372135  0.5864251 0.05975604 0.0000000
# }
```
