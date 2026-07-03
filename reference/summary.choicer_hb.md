# Summarize a hierarchical Bayes fit

Posterior summaries (mean, SD, equal-tailed credible interval) for the
population coefficients \\b\\, the mean-function coefficients
\\\theta\\, the alternative-effect variance \\\sigma_d^2\\ (and, for the
HMNP, the raw shock variance trace), plus the \\\delta_j\\ / \\\xi_j\\
quality ladder, acceptance diagnostics, and the split-R-hat table when
multiple chains were run.

## Usage

``` r
# S3 method for class 'choicer_hb'
summary(object, prob = 0.95, ...)
```

## Arguments

- object:

  A `choicer_hmnl` or `choicer_hmnp` object.

- prob:

  Probability mass of the equal-tailed credible interval (default 0.95).

- ...:

  Additional arguments (ignored).

## Value

A `summary.choicer_hb` object.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 50, T = 2, J = 3, seed = 42)
fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                    person_col = "pid",
                    mcmc = list(R = 300, burn = 100))
#> MCMC run time 0h:0m:0.01s
#> Warning: split-R-hat exceeds 1.05 for: (Intercept). Consider a longer run.
summary(fit)
#> Hierarchical Bayesian Multinomial Logit (HMNL) model
#> 
#> Population coefficients b (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> x1           0.588077   0.301117  -0.002989   0.583537   1.150561
#> x2          -0.552016   0.307059  -1.189243  -0.567522   0.043198
#> 
#> Delta mean function theta (posterior):
#> Parameter          Mean         SD       2.5%     Median      97.5%
#> (Intercept)    0.158710   0.418980  -0.776411   0.218686   0.696385
#> 
#> Alternative-effect variance (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> sigma_d^2    0.260585   0.860165   0.000090   0.064074   1.381520
#> 
#> Quality ladder (delta = mean utility vs the outside option; xi = delta - z'theta):
#>  alternative delta_mean delta_sd xi_mean  xi_sd
#>            1     0.2671   0.2954  0.1084 0.2976
#>            2    -0.0297   0.4123 -0.1884 0.3071
#>            3     0.2296   0.3069  0.0708 0.2729
#> 
#> Mean acceptance: beta 0.24, delta 0.43
#> 
#> split-R-hat:
#>          x1          x2 (Intercept)    sigma_d2 
#>       1.037       1.048       1.114       1.016 
#> 
#> Respondents: 50  Choice situations: 100  Alternatives: 3 
#> Draws kept: 200  Chains: 1 
#> MCMC run time 0h:0m:0.01s 
# }
```
