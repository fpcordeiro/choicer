# Split-\\\widehat{R}\\ convergence diagnostic

Computes the split-\\\widehat{R}\\ (potential scale reduction factor) of
Gelman et al. for each column of a matrix of posterior draws. Every
chain is split in half, so the diagnostic detects non-stationarity
within a single chain as well as disagreement across chains; values near
1 indicate convergence, and values above roughly 1.05 warrant a longer
run.

## Usage

``` r
rhat(draws)
```

## Arguments

- draws:

  A matrix of posterior draws (rows = iterations, columns = parameters)
  for a single chain, or a list of such matrices (one per chain,
  identical dimensions).

## Value

Named numeric vector with one \\\widehat{R}\\ per parameter (`NA` for
parameters with zero variance).

## Examples

``` r
set.seed(42)
draws <- matrix(rnorm(2000), ncol = 2,
                dimnames = list(NULL, c("a", "b")))
rhat(draws)          # ~1: white noise is stationary
#>         a         b 
#> 0.9990172 1.0003443 
drifting <- cbind(a = cumsum(rnorm(1000)))
rhat(drifting)       # >> 1: a random walk is not
#>        a 
#> 1.817129 
```
