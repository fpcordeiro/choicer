# Simulate hierarchical multinomial logit data

Generates synthetic panel choice data from the hierarchical (random
coefficients + alternative-level random effects) logit DGP: respondents
`i = 1..N` face `T` choice situations each, with utilities \$\$U\_{ijt}
= x\_{ijt}'\gamma_i + \delta_j + \epsilon\_{ijt}, \qquad U\_{iot} =
\epsilon\_{iot},\$\$ i.i.d. Gumbel shocks (including a shock on the
outside option, whose systematic utility is 0), \\\beta_i \sim N(\beta,
W)\\ with \\\gamma\_{ik} = \beta\_{ik}\\ or \\\exp(\beta\_{ik})\\ per
`rc_dist`, and \\\delta_j = z_j'\theta + \xi_j\\ with \\\xi_j \sim N(0,
\sigma_d^2)\\. Covariates are Uniform(-1, 1); the alternative-level
covariates `z*` are constant within each alternative.

## Usage

``` r
simulate_hmnl_data(
  N = 500,
  T = 10,
  J = 4,
  beta = c(0.8, -0.6),
  W = NULL,
  theta = c(0.5, -0.4),
  sigma_d = 0.5,
  Z = NULL,
  rc_dist = NULL,
  include_outside = TRUE,
  seed = 123,
  vary_choice_set = FALSE
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
  `K_x = length(beta)`; chain scale for log-normal coordinates).

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

- rc_dist:

  Integer vector (length `K_x`): `0L` for normal, `1L` for log-normal
  coordinates. Default `NULL` is all-normal.

- include_outside:

  Logical; if `TRUE` (default) every choice set also contains an outside
  option with systematic utility 0 and its own Gumbel shock. The outside
  good is *implicit* in the returned data (matching the
  [`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md)
  convention): no physical row is emitted, and a choice situation the
  outside option wins has an all-zeros `choice` column.

- seed:

  Random seed (`NULL` skips
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

- vary_choice_set:

  Logical; if `TRUE` choice-set size is sampled uniformly from `2:J` per
  task. Default `FALSE`.

## Value

A `choicer_sim` object. `true_params` contains `beta`, `W`, `theta`,
`sigma_d`, the realized `delta` and `xi` vectors, the full mean-function
design `Z` (intercept first), and `rc_dist`.

## Details

Log-normal coordinates are reported on the *chain* (log) scale in
`true_params$beta` — the scale on which the estimator's hierarchy
operates — while entering utility as `exp(beta_ik)`.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 100, T = 4, J = 4, seed = 123)
print(sim)
#> <choicer_sim: hmnl>
#>   settings:
#>     N = 100
#>     T = 4
#>     J = 4
#>     K_x = 2
#>     P = 2
#>     include_outside = TRUE
#>     vary_choice_set = FALSE
#>   rows in $data: 1600
#>   true_params: beta, W, theta, sigma_d, delta, xi, Z, rc_dist
sim$true_params$delta
#> [1] 1.4492921 0.3046101 0.6374623 1.0511186
# }
```
