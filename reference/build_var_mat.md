# Reconstruct variance matrix L from L_params

Reconstruct variance matrix L from L_params

## Usage

``` r
build_var_mat(L_params, K_w, rc_correlation)
```

## Arguments

- L_params:

  flattened choleski decomposition version of the random coefficient
  parameters matrix

- K_w:

  dimension of the random coefficient parameter (symmetric) matrix

- rc_correlation:

  whether random coefficients are correlated

## Value

matrix equal to LL', where L is the choleski decomposition of random
coefficient matrix

## Examples

``` r
L_params <- c(log(1.0), 0.3, log(0.5))
Sigma <- build_var_mat(L_params, K_w = 2, rc_correlation = TRUE)
Sigma  # 2x2 covariance matrix
#>      [,1] [,2]
#> [1,]  1.0 0.30
#> [2,]  0.3 0.34
```
