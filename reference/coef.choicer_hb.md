# Extract posterior means from a hierarchical Bayes fit

Extract posterior means from a hierarchical Bayes fit

## Usage

``` r
# S3 method for class 'choicer_hb'
coef(object, component = c("beta", "theta", "delta", "xi"), ...)
```

## Arguments

- object:

  A `choicer_hmnl` or `choicer_hmnp` object.

- component:

  Which block to return: `"beta"` (population means b, default),
  `"theta"` (delta mean function), `"delta"` (alternative effects), or
  `"xi"` (unobserved quality, delta - z'theta).

- ...:

  Additional arguments (ignored).

## Value

Named numeric vector of posterior means.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 50, T = 2, J = 3, seed = 42)
fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                    person_col = "pid",
                    mcmc = list(R = 300, burn = 100))
#> MCMC run time 0h:0m:0.01s
#> Warning: split-R-hat exceeds 1.05 for: (Intercept). Consider a longer run.
coef(fit)
#>         x1         x2 
#>  0.5880772 -0.5520164 
coef(fit, component = "delta")
#>           1           2           3 
#>  0.26708115 -0.02969955  0.22955716 
# }
```
