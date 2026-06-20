# Predicted aggregate market shares for Mixed Logit

Exported wrapper around the internal `mxl_predict_shares_internal`.
Parses `theta` using the standard parameter ordering and returns the
simulated weighted-average market shares.

## Usage

``` r
mxl_predict_shares(
  theta,
  X,
  W,
  alt_idx,
  M,
  weights,
  eta_draws,
  rc_dist,
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

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- eta_draws:

  Array with draws; K_w x S x N

- rc_dist:

  K_w vector indicating distribution (0=normal, 1=log-normal)

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

Vector of length J (or J+1 with outside option) of predicted shares.
