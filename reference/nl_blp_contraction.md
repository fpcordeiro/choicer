# BLP95 contraction mapping for the Nested Logit model

Damped iterative fixed point recovering delta given target shares, using
the NL probability structure. `damping = 1` reproduces the plain BLP
update.

## Usage

``` r
nl_blp_contraction(
  delta,
  target_shares,
  X,
  beta,
  lambda,
  alt_idx,
  nest_idx,
  M,
  weights,
  include_outside_option = FALSE,
  damping = 1,
  tol = 1e-08,
  max_iter = 1000L
)
```

## Arguments

- delta:

  J x 1 vector with initial guess for deltas (ASCs).

- target_shares:

  vector with target shares (outside-option share first when present).

- X:

  sum(M) x K design matrix with covariates.

- beta:

  K x 1 vector with fixed coefficients.

- lambda:

  full nest dissimilarity vector of length n_nests (singletons = 1).

- alt_idx:

  sum(M) x 1 vector with indices of alternatives; 1-based indexing.

- nest_idx:

  J x 1 vector with nest indices for each alternative; 1-based indexing.

- M:

  N x 1 vector with number of alternatives for each individual.

- weights:

  N x 1 vector with weights for each observation.

- include_outside_option:

  whether to include outside option normalized to V=0, lambda=1.

- damping:

  damping factor for the update (default 1.0 = plain BLP).

- tol:

  convergence tolerance.

- max_iter:

  maximum number of iterations.

## Value

vector with contraction's delta (ASCs) output.

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
beta <- coef(fit)[fit$param_map$beta]
lambda <- rep(1, length(unique(fit$data$nest_idx)))
lambda[as.integer(names(which(table(fit$data$nest_idx) > 1)))] <-
  coef(fit)[fit$param_map$lambda]
delta <- nl_blp_contraction(rep(0, J), rep(1/J, J), fit$data$X, beta, lambda,
  fit$data$alt_idx, fit$data$nest_idx, fit$data$M, fit$data$weights)
#> Warning: Maximum iterations reached without convergence.
delta
#>             [,1]
#> [1,]   0.0000000
#> [2,]  -0.1201658
#> [3,] 108.3182046
#> [4,] 108.3827250
# }
```
