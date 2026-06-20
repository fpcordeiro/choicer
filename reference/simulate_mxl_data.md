# Simulate mixed logit data

Generates synthetic choice data with random coefficients drawn from a
multivariate normal (optionally log-normal per dimension) and an
additional mean shifter `mu`. Random coefficients are parameterized via
the lower Cholesky factor of `Sigma`. Covariates are Uniform(-1, 1) by
default; columns named in `price_cols` are drawn as `-Uniform(0.1, 3)`
to mimic strictly-negative price variables.

## Usage

``` r
simulate_mxl_data(
  N = 5000,
  J = 4,
  beta = c(0.8, -0.6),
  delta = NULL,
  mu = NULL,
  Sigma = matrix(c(1, 0.5, 0.5, 1.5), nrow = 2),
  rc_dist = NULL,
  rc_correlation = NULL,
  price_cols = NULL,
  seed = 123,
  outside_option = TRUE,
  vary_choice_set = TRUE
)
```

## Arguments

- N:

  Number of choice situations.

- J:

  Number of inside alternatives.

- beta:

  Fixed coefficients for `x1..x{K_x}` (length `K_x = length(beta)`).

- delta:

  ASCs for inside alternatives (length `J`). Defaults to an alternating
  pattern of `c(0.5, -0.5)`.

- mu:

  Mean shifter for random coefficients (length `K_w = ncol(Sigma)`).
  Defaults to a zero vector.

- Sigma:

  Covariance matrix of random coefficients (square, `K_w x K_w`).

- rc_dist:

  Integer vector (length `K_w`): `0L` for normal, `1L` for log-normal.
  Default `NULL` is treated as all-normal.

- rc_correlation:

  Logical; if `NULL` (default) it is auto-detected from the off-diagonal
  entries of `Sigma`.

- price_cols:

  Character vector of `w*` column names to draw as `-Uniform(0.1, 3)`
  instead of `Uniform(-1, 1)`. Default `NULL`.

- seed:

  Random seed (`NULL` skips
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

- outside_option:

  Logical; include outside option with `alt = 0`.

- vary_choice_set:

  Logical; if `TRUE` (default) choice set size is sampled uniformly from
  `2:J`.

## Value

A `choicer_sim` object. `true_params` includes `beta`, `delta`, `Sigma`,
`L_params` (packed Cholesky parameters), `mu`, `rc_dist`,
`rc_correlation`.

## Details

Random coefficients are constructed to match the estimator's
parameterization in `src/mxlogit.cpp`. For every dimension the raw draw
is `L %*% eta` where `eta ~ N(0, I)`. A normal random coefficient
(`rc_dist = 0`) is then `gamma_k = mu_k + (L %*% eta)_k`. A log-normal
random coefficient (`rc_dist = 1`) follows the shifted log-normal
`beta_k = exp(mu_k) + exp((L %*% eta)_k)` – not the textbook
`exp(mu_k + sigma_k * eta)` – so `mu_k` in `true_params$mu` is on the
same scale the estimator recovers and
[`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md)
can compare like-for-like.

## Examples

``` r
# \donttest{
sim <- simulate_mxl_data(N = 1000, J = 4, seed = 123)
print(sim)
#> <choicer_sim: mxl>
#>   settings:
#>     N = 1000
#>     J = 4
#>     K_x = 2
#>     K_w = 2
#>     outside_option = TRUE
#>     vary_choice_set = TRUE
#>   rows in $data: 3989
#>   true_params: beta, delta, Sigma, L_params, mu, rc_dist, rc_correlation
# }
```
