# Analytical Hessian of the negated log-likelihood for the Nested Logit model

Computes the exact (analytical) Hessian of the negated log-likelihood
for the Nested Logit model using OpenMP parallelisation with
thread-local accumulators. Covers all parameter blocks: beta-beta,
beta-lambda, beta-delta, lambda-lambda, lambda-delta, and delta-delta.
Singleton nests (lambda fixed to 1, not estimated) contribute no rows or
columns to the lambda blocks.

## Usage

``` r
nl_loglik_hessian_parallel(
  theta,
  X,
  alt_idx,
  choice_idx,
  nest_idx,
  M,
  weights,
  use_asc = TRUE,
  include_outside_option = FALSE
)
```

## Arguments

- theta:

  (K + n_non_singleton_nests + n_delta) parameter vector. Order:
  `[beta (K), lambda (n_non_singleton_nests), delta (n_delta)]`. Same
  layout as `nl_loglik_gradient_parallel`.

- X:

  sum(M) x K design matrix of covariates.

- alt_idx:

  sum(M)-length integer vector of 1-based alternative indices.

- choice_idx:

  N-length integer vector of 1-based chosen alternative indices; 0
  indicates the outside option was chosen.

- nest_idx:

  J-length integer vector of 1-based nest indices for each inside
  alternative.

- M:

  N-length integer vector of alternative-set sizes.

- weights:

  N-length numeric vector of individual weights.

- use_asc:

  Logical; whether alternative-specific constants are included.

- include_outside_option:

  Logical; whether an outside option (V=0) is present.

## Value

A symmetric (P x P) matrix: the Hessian of the negated log-likelihood
evaluated at `theta`. Structurally identical to the output of
`nl_loglik_numeric_hessian`; suitable for `invert_hessian()`.

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
d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
K_x <- ncol(d$X)
K_l <- sum(table(d$nest_idx) > 1)   # number of non-singleton nests (= 2)
theta <- c(rep(0, K_x), rep(0.8, K_l), rep(0, J - 1))
H <- nl_loglik_hessian_parallel(theta, d$X, d$alt_idx, d$choice_idx,
  d$nest_idx, d$M, d$weights)
dim(H)
#> [1] 7 7
# }
```
