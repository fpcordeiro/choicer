# Print a hierarchical Bayes fit

Print a hierarchical Bayes fit

## Usage

``` r
# S3 method for class 'choicer_hb'
print(x, ...)
```

## Arguments

- x:

  A `choicer_hmnl` or `choicer_hmnp` object.

- ...:

  Additional arguments (ignored).

## Value

The object invisibly.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 50, T = 2, J = 3, seed = 42)
fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                    person_col = "pid",
                    mcmc = list(R = 300, burn = 100))
#> MCMC run time 0h:0m:0.01s
#> Warning: split-R-hat exceeds 1.05 for: (Intercept). Consider a longer run.
print(fit)
#> Hierarchical Bayesian Multinomial Logit (HMNL) model
#> Posterior-mean population coefficients (b):
#>      x1      x2 
#>  0.5881 -0.5520 
#> 
#> Alternatives: 3  Respondents: 50  Choice situations: 100 
#> Draws kept: 200  (R = 300 , burn = 100 , thin = 1 )
# }
```
