# Gibbs sampler for the Bayesian multinomial probit model

Runs the McCulloch-Rossi (1994) Gibbs sampler with Albert-Chib data
augmentation for the multinomial probit model in utility differences
against a base alternative. The chain operates on the non-identified
parameterization (unrestricted `Sigma`); identified quantities are
obtained by normalizing each draw by `sigma_11` (handled by
[`run_mnprobit`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md)).

## Usage

``` r
mnp_gibbs(X, y, p, beta_bar, A, nu, V, R, burn, thin, seed, trace = 0L)
```

## Arguments

- X:

  (N\*p) x K stacked design matrix of utility differences. Rows are
  grouped by choice situation, with the p = J - 1 difference rows of
  situation i ordered by alternative.

- y:

  N vector of choices: 0 for the base alternative, j in 1..p for the
  j-th non-base alternative.

- p:

  Number of utility differences (J - 1).

- beta_bar:

  K vector, prior mean of beta.

- A:

  K x K prior precision matrix of beta.

- nu:

  Inverse-Wishart prior degrees of freedom (\>= p).

- V:

  p x p inverse-Wishart prior scale matrix.

- R:

  Total number of Gibbs iterations.

- burn:

  Number of initial iterations to discard (0 \<= burn \< R).

- thin:

  Keep every thin-th post-burn-in draw.

- seed:

  Master RNG seed (non-negative; all streams derive from it).

- trace:

  Print progress every `trace` iterations (0 = silent).

## Value

List with `betadraw` (R_keep x K), `sigmadraw` (R_keep x p(p+1)/2, lower
triangle of Sigma in row-major order), and `R_keep`.

## Details

The latent-utility sweep is parallelized with OpenMP across choice
situations (they are conditionally independent given `beta` and
`Sigma`). Each (iteration, observation) pair uses its own RNG stream, so
draws are reproducible independent of the number of threads (see
[`set_num_threads()`](https://fpcordeiro.github.io/choicer/reference/set_num_threads.md)).

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 100; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt          x1           x2
#>      <int> <int>       <num>        <num>
#>   1:     1     1  1.37095845 -0.004620768
#>   2:     1     2 -0.56469817  0.760242168
#>   3:     1     3  0.36312841  0.038990913
#>   4:     2     1  0.63286260  0.735072142
#>   5:     2     2  0.40426832 -0.146472627
#>  ---                                     
#> 296:    99     2 -0.47733551  0.160327395
#> 297:    99     3 -0.16626149 -0.433641942
#> 298:   100     1  0.86256338  1.537412419
#> 299:   100     2  0.09734049 -2.170246577
#> 300:   100     3 -1.62561674  1.027004619
dt[, choice := 0L]
#>         id   alt          x1           x2 choice
#>      <int> <int>       <num>        <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768      0
#>   2:     1     2 -0.56469817  0.760242168      0
#>   3:     1     3  0.36312841  0.038990913      0
#>   4:     2     1  0.63286260  0.735072142      0
#>   5:     2     2  0.40426832 -0.146472627      0
#>  ---                                            
#> 296:    99     2 -0.47733551  0.160327395      0
#> 297:    99     3 -0.16626149 -0.433641942      0
#> 298:   100     1  0.86256338  1.537412419      0
#> 299:   100     2  0.09734049 -2.170246577      0
#> 300:   100     3 -1.62561674  1.027004619      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt          x1           x2 choice
#>      <int> <int>       <num>        <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768      0
#>   2:     1     2 -0.56469817  0.760242168      0
#>   3:     1     3  0.36312841  0.038990913      1
#>   4:     2     1  0.63286260  0.735072142      0
#>   5:     2     2  0.40426832 -0.146472627      1
#>  ---                                            
#> 296:    99     2 -0.47733551  0.160327395      0
#> 297:    99     3 -0.16626149 -0.433641942      0
#> 298:   100     1  0.86256338  1.537412419      0
#> 299:   100     2  0.09734049 -2.170246577      0
#> 300:   100     3 -1.62561674  1.027004619      1
d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
out <- mnp_gibbs(d$X, d$y, d$p,
  beta_bar = rep(0, d$K), A = 0.01 * diag(d$K),
  nu = d$p + 3, V = (d$p + 3) * diag(d$p),
  R = 500, burn = 100, thin = 1, seed = 42)
colMeans(out$betadraw)
#> [1]  0.104035441  0.006645003 -0.661673087 -0.708334851
# }
```
