# Per-observation simulated choice probabilities for Mixed Logit

Returns the simulated choice probability for each (individual,
alternative) row of `X`, averaged over the supplied Halton draws.
Mirrors `mnl_predict`.

## Usage

``` r
mxl_predict(
  theta,
  X,
  W,
  alt_idx,
  M,
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

  whether the outside option is present

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

List with `choice_prob` (length sum(M)), `utility` (length sum(M),
simulated mean of the deterministic + W\*gamma component), and, when
`include_outside_option = TRUE`, `choice_prob_outside` (length N).
