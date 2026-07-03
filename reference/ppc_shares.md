# Posterior-predictive share check for hierarchical Bayes fits

Compares each alternative's observed take rate (share of choice
situations in which it was chosen, including the outside option) with
its posterior-predictive share from
[`predict.choicer_hb()`](https://fpcordeiro.github.io/choicer/reference/predict.choicer_hb.md).
Large systematic gaps indicate model misfit — e.g. a missing covariate
or an outside-option share the delta level cannot rationalize.

## Usage

``` r
ppc_shares(object, n_draws = 200L)
```

## Arguments

- object:

  A `choicer_hmnl` or `choicer_hmnp` fit (with `keep_data = TRUE`).

- n_draws:

  Posterior draws to integrate over (default 200).

## Value

A `data.table` with columns `alternative`, `observed`, `predicted`,
`lower`, `upper` (95% posterior-predictive interval), and `covered` (is
the observed share inside the interval).

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 100, T = 3, J = 4, seed = 42)
fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                    person_col = "pid",
                    mcmc = list(R = 500, burn = 200))
#> MCMC run time 0h:0m:0.04s
#> Warning: split-R-hat exceeds 1.05 for: (Intercept). Consider a longer run.
ppc_shares(fit)
#>    alternative  observed predicted      lower     upper covered
#>         <char>     <num>     <num>      <num>     <num>  <lgcl>
#> 1:           1 0.1900000 0.1993650 0.16013925 0.2340801    TRUE
#> 2:           2 0.2133333 0.2080083 0.16242845 0.2597030    TRUE
#> 3:           3 0.2766667 0.2596282 0.21971533 0.3081084    TRUE
#> 4:           4 0.2000000 0.2000520 0.15744978 0.2346050    TRUE
#> 5:   (outside) 0.1200000 0.1329464 0.09376908 0.1955407    TRUE
# }
```
