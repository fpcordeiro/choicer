# Simulate hierarchical multinomial probit data

Generates synthetic panel choice data from the hierarchical probit DGP
with iid normal utility shocks: \$\$U\_{ijt} = x\_{ijt}'\beta_i +
\delta_j + \epsilon\_{ijt}, \qquad U\_{iot} = \epsilon\_{iot}, \qquad
\epsilon \sim N(0, \sigma^2),\$\$ choice by argmax within the task. The
outside option is stochastic — it carries its own \\N(0, \sigma^2)\\
shock on top of systematic utility 0, exactly as in the estimator.
\\\beta_i \sim N(\beta, W)\\ (normal only) and \\\delta_j = z_j'\theta +
\xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\, as in
[`simulate_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_hmnl_data.md).

## Usage

``` r
simulate_hmnp_data(
  N = 500,
  T = 10,
  J = 4,
  beta = c(0.8, -0.6),
  W = NULL,
  theta = c(0.5, -0.4),
  sigma_d = 0.5,
  Z = NULL,
  include_outside = TRUE,
  seed = 123,
  vary_choice_set = FALSE,
  sigma = 1
)
```

## Arguments

- N:

  Number of respondents.

- T:

  Number of choice situations per respondent.

- J:

  Number of inside alternatives.

- beta:

  Population means of the structural random coefficients (length
  `K_x = length(beta)`); all coordinates are normal.

- W:

  Covariance of the random coefficients (`K_x` x `K_x`). Defaults to
  `diag(0.5, K_x)`.

- theta:

  Mean-function coefficients for \\\delta_j = z_j'\theta + \xi_j\\; the
  first entry is the intercept, entries `2..P` load on the
  alternative-level covariates.

- sigma_d:

  Standard deviation of the alternative-level effects \\\xi_j\\.

- Z:

  Optional `J x (length(theta) - 1)` matrix of alternative-level
  covariates (excluding the intercept). Default `NULL` draws them
  Uniform(-1, 1).

- include_outside:

  Logical; if `TRUE` (default) every choice set also contains an outside
  option with systematic utility 0 and its own normal shock. The outside
  good is *implicit* in the returned data (matching the
  [`prepare_hmnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnp_data.md)
  convention): no physical row is emitted, and a choice situation the
  outside option wins has an all-zeros `choice` column.

- seed:

  Random seed (`NULL` skips
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

- vary_choice_set:

  Logical; if `TRUE` choice-set size is sampled uniformly from `2:J` per
  task. Default `FALSE`.

- sigma:

  Standard deviation of the iid utility shocks (DGP scale).

## Value

A `choicer_sim` object. `true_params` contains `beta`, `W`, `theta`,
`sigma_d`, the realized `delta` and `xi`, and the full mean-function
design `Z` — all on the identified scale (see Details).

## Details

The iid-probit likelihood identifies parameters only up to the common
scale \\\sigma\\, so `true_params` is reported on the *identified*
scale: `beta` \\= \beta/\sigma\\, `W` \\= W/\sigma^2\\, `theta` \\=
\theta/\sigma\\, `sigma_d` \\= \sigma_d/\sigma\\, `delta` \\=
\delta/\sigma\\, `xi` \\= \xi/\sigma\\. With the default `sigma = 1` the
DGP scale and the identified scale coincide.

## Examples

``` r
# \donttest{
sim <- simulate_hmnp_data(N = 100, T = 4, J = 4, seed = 123)
print(sim)
#> <choicer_sim: hmnp>
#>   settings:
#>     N = 100
#>     T = 4
#>     J = 4
#>     K_x = 2
#>     P = 2
#>     include_outside = TRUE
#>     vary_choice_set = FALSE
#>     sigma = 1
#>   rows in $data: 1600
#>   true_params: beta, W, theta, sigma_d, delta, xi, Z
sim$true_params$delta
#> [1] 1.4492921 0.3046101 0.6374623 1.0511186
# }
```
