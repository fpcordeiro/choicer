# Simulate nested logit data

Generates synthetic choice data with nested logit probabilities computed
analytically (log-sum-exp over inclusive values), then samples choices
from the implied multinomial. The outside option (`j = 0`) sits in a
singleton nest with `lambda = 1`.

## Usage

``` r
simulate_nl_data(
  N = 10000,
  beta = c(1.5, -0.8),
  delta = c(`1` = 0.5, `2` = 0.3, `3` = -0.2, `4` = -0.5, `5` = 0.4),
  nests = list(c(1, 2), c(3, 4, 5)),
  lambdas = c(0.8, 0.2),
  seed = 123
)
```

## Arguments

- N:

  Number of choice situations.

- beta:

  Fixed coefficients for covariates `X, W` (length 2 by default).

- delta:

  Named numeric vector of ASCs for inside alternatives.

- nests:

  List of integer vectors defining nest membership for inside
  alternatives.

- lambdas:

  Numeric vector of dissimilarity parameters, one per nest.

- seed:

  Random seed (`NULL` skips
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

## Value

A `choicer_sim` object. `true_params` includes `beta`, `delta`,
`lambdas`; `settings` includes the `nest_structure`. The returned `data`
retains a `nest` column (integer, with `0L` for the outside option) for
convenient use with
[`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md).

## Note

Unlike
[`simulate_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnl_data.md)
and
[`simulate_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mxl_data.md),
this function does not expose `outside_option` or `vary_choice_set`
flags. The outside option (`j = 0`) is always present as a singleton
nest with `lambda = 1`, and every individual faces the full set of
inside alternatives. Add these flags if downstream use cases need them.

## Examples

``` r
# \donttest{
sim <- simulate_nl_data(N = 2000, seed = 123)
print(sim)
#> <choicer_sim: nl>
#>   settings:
#>     N = 2000
#>     J_inside = 5
#>     nest_structure = (1,2), (3,4,5)
#>   rows in $data: 12000
#>   true_params: beta, delta, lambdas
# }
```
