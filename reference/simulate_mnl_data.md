# Simulate multinomial logit data

Generates synthetic choice data with i.i.d. Gumbel errors, optionally
with varying choice-set sizes and an outside option (alt = 0). Choices
are determined by argmax of utility; covariates are drawn as Uniform(-1,
1).

## Usage

``` r
simulate_mnl_data(
  N = 5000,
  J = 5,
  beta = c(0.8, -0.6),
  delta = NULL,
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

  Alternative-specific constants for inside alternatives (length `J`).
  Defaults to an alternating pattern of `c(0.5, -0.5)`.

- seed:

  Random seed. Pass `NULL` to skip
  [`set.seed()`](https://rdrr.io/r/base/Random.html) (useful inside
  [`monte_carlo()`](https://fpcordeiro.github.io/choicer/reference/monte_carlo.md)
  where the caller manages RNG).

- outside_option:

  Logical; if `TRUE` (default) an outside option with `alt = 0` and zero
  covariates is added to every choice set.

- vary_choice_set:

  Logical; if `TRUE` (default) choice set size is sampled uniformly from
  `2:J`; if `FALSE` every individual faces all `J` inside alternatives.

## Value

A `choicer_sim` object.

## Examples

``` r
# \donttest{
sim <- simulate_mnl_data(N = 1000, J = 5, seed = 123)
print(sim)
#> <choicer_sim: mnl>
#>   settings:
#>     N = 1000
#>     J = 5
#>     K_x = 2
#>     outside_option = TRUE
#>     vary_choice_set = TRUE
#>   rows in $data: 4500
#>   true_params: beta, delta
# }
```
