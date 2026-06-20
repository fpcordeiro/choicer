# Utility to compute analytical Jacobian of random coefficient matrix transformed by vech (dVech(Sigma) / dTheta)

Utility to compute analytical Jacobian of random coefficient matrix
transformed by vech (dVech(Sigma) / dTheta)

## Usage

``` r
jacobian_vech_Sigma(L_params, K_w, rc_correlation = TRUE)
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

Jacobian (dVech(Sigma) / dTheta)

## Examples

``` r
L_params <- c(log(0.8), 0.2, log(0.6))
J_mat <- jacobian_vech_Sigma(L_params, K_w = 2, rc_correlation = TRUE)
dim(J_mat)  # 3 x 3 for K_w=2 correlated
#> [1] 3 3
```
