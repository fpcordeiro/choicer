# Diversion ratios for Mixed Logit (simulated, derivative-based)

Computes the matrix of attribute-based diversion ratios for a fitted
Mixed Logit model. DR(k, j) is the fraction of demand lost by
alternative `j` that is captured by alternative `k` when a marginal
change in alternative j's `elast_var` attribute reduces s_j.

## Usage

``` r
mxl_diversion_ratios_parallel(
  theta,
  X,
  W,
  alt_idx,
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

- M:

  N x 1 vector with number of alternatives for each individual

- weights:

  N x 1 vector with weights for each observation

- eta_draws:

  Array with draws; K_w x S x N

- rc_dist:

  K_w vector indicating distribution (0=normal, 1=log-normal)

- elast_var_idx:

  1-based index of the perturbed variable

- is_random_coef:

  TRUE if the variable is in W (random coef), FALSE if in X (fixed)

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

J x J (or (J+1) x (J+1)) matrix of diversion ratios with zero diagonal.

## Details

In MNL the per-draw realized coefficient is a constant, so it cancels in
the ratio and the result is independent of the variable chosen. In MXL,
the realized coefficient \\\beta\_{ik}^s\\ varies across individuals and
draws, so the diversion ratio depends on which attribute is perturbed.
For a variable with a fixed coefficient the dependence again vanishes
(the constant cancels); for a random-coefficient variable it does not.
