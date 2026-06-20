# Mixed logit and preference heterogeneity

The mixed (random-coefficients) logit lets tastes vary across people.
Instead of a single coefficient on each random attribute, choicer
estimates a *distribution* of coefficients, which relaxes the IIA
property and allows realistic substitution patterns. Estimation is by
simulated maximum likelihood using Halton draws, with the likelihood,
gradient and Hessian evaluated in parallel C++.

``` r

library(choicer)
set_num_threads(2)
```

## Simulate correlated random coefficients

[`simulate_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mxl_data.md)
draws choices in which the coefficients on `w1` and `w2` are themselves
random and *correlated* across the population.

``` r

sim <- simulate_mxl_data(N = 2000, J = 4, seed = 1)
sim
#> <choicer_sim: mxl>
#>   settings:
#>     N = 2000
#>     J = 4
#>     K_x = 2
#>     K_w = 2
#>     outside_option = TRUE
#>     vary_choice_set = TRUE
#>   rows in $data: 7982
#>   true_params: beta, delta, Sigma, L_params, mu, rc_dist, rc_correlation
```

## Fit

A robust recipe for mixed logit: warm-start from a plain MNL, scale the
variables so the Hessian is well conditioned, and use enough Halton
draws. Here we estimate a full (correlated) covariance of the random
coefficients.

``` r

fit <- run_mxlogit(
  data            = sim$data,
  id_col          = "id",
  alt_col         = "alt",
  choice_col      = "choice",
  covariate_cols  = c("x1", "x2"),  # fixed coefficients
  random_var_cols = c("w1", "w2"),  # random coefficients
  rc_correlation  = TRUE,           # estimate their full covariance
  S               = 100L,           # Halton draws per person
  draws           = "generate",     # generate draws on the fly (low memory)
  seed            = 7L,
  scale_vars      = "sd",           # condition the Hessian across blocks
  se_method       = "bhhh"
)
#> Optimization run time 0h:0m:1s
summary(fit)
#> Mixed Logit (MXL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           0.776408   0.070138  11.0697  0.00e+00  ***
#> x2          -0.675643   0.068318  -9.8897  0.00e+00  ***
#> Sigma_11     0.905183   0.451785   2.0036  4.51e-02  *
#> Sigma_21     0.894700   0.310679   2.8798  3.98e-03  **
#> Sigma_22     1.884033   0.575835   3.2718  1.07e-03  **
#> ASC_1        0.536159   0.078972   6.7893  1.13e-11  ***
#> ASC_2       -0.477738   0.105075  -4.5466  5.45e-06  ***
#> ASC_3        0.591089   0.077545   7.6225  2.49e-14  ***
#> ASC_4       -0.513880   0.105365  -4.8771  1.08e-06  ***
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Random coefficient covariance (Sigma):
#>        w1     w2
#> w1 0.9052 0.8947
#> w2 0.8947 1.8840
#> 
#> Std. Errors: BHHH (OPG) 
#> Log-likelihood: -2435.84 
#> AIC: 4889.67  | BIC: 4940.08 
#> McFadden R2: 0.106 (adj: 0.102) | Hit rate: 0.452 
#> N: 2000  | Parameters: 9 
#> Optimization time: 1.43 s
#> Convergence: 3 ( NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. )
```

> **Tip.** For real applications increase the number of draws (`S`)
> until your estimates are stable, and keep `scale_vars = "sd"`. If the
> solver struggles, pass an explicit `theta_init` (for example the MNL
> coefficients) and bounds on the Cholesky diagonal. See
> `inst/simulations/mxl_simulation.R` for a fully hardened example.

## Parameter recovery

``` r

recovery_table(fit, sim$true_params)
#> <choicer_recovery> model=choicer_mxl level=0.95
#>    parameter  group    true estimate     se    bias rel_bias_pct z_vs_true
#>       <char> <char>   <num>    <num>  <num>   <num>        <num>     <num>
#> 1:        x1   beta  0.8000   0.7764 0.0701 -0.0236       -2.949   -0.3364
#> 2:        x2   beta -0.6000  -0.6756 0.0683 -0.0756       12.607   -1.1072
#> 3:      L_11  sigma  0.0000  -0.0498 0.2496 -0.0498           NA   -0.1996
#> 4:      L_21  sigma  0.5000   0.9404 0.3057  0.4404       88.078    1.4406
#> 5:      L_22  sigma  0.1116  -0.0002 0.3566 -0.1117     -100.136   -0.3133
#> 6:     ASC_1    asc  0.5000   0.5362 0.0790  0.0362        7.232    0.4579
#> 7:     ASC_2    asc -0.5000  -0.4777 0.1051  0.0223       -4.452    0.2119
#> 8:     ASC_3    asc  0.5000   0.5911 0.0775  0.0911       18.218    1.1747
#> 9:     ASC_4    asc -0.5000  -0.5139 0.1054 -0.0139        2.776   -0.1317
#>    lower_ci upper_ci covers
#>       <num>    <num> <lgcl>
#> 1:   0.6389   0.9139   TRUE
#> 2:  -0.8095  -0.5417   TRUE
#> 3:  -0.5389   0.4393   TRUE
#> 4:   0.3412   1.5396   TRUE
#> 5:  -0.6991   0.6988   TRUE
#> 6:   0.3814   0.6909   TRUE
#> 7:  -0.6837  -0.2718   TRUE
#> 8:   0.4391   0.7431   TRUE
#> 9:  -0.7204  -0.3074   TRUE
```

The `beta` rows are the fixed coefficients, the `sigma` rows describe
the covariance of the random coefficients (its Cholesky elements), and
the `asc` rows are the alternative-specific constants.

## Substitution is no longer share-proportional

Because tastes are heterogeneous, mixed logit breaks IIA: people who
like one alternative tend to like its close substitutes, so demand
diverts to similar alternatives rather than in proportion to shares.
Diversion therefore depends on *which* attribute is changing, so
[`diversion_ratios()`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.md)
takes a `wrt_var`:

``` r

elasticities(fit, elast_var = "x2")
#>   0         1         2         3         4
#> 0 0 -0.015656 -0.009793 -0.025914 -0.009414
#> 1 0 -0.023929 -0.005506 -0.017254 -0.005718
#> 2 0 -0.010723 -0.016649 -0.017014 -0.007051
#> 3 0 -0.008498 -0.006243 -0.006971 -0.006645
#> 4 0 -0.011956 -0.006111 -0.021578 -0.020221
diversion_ratios(fit, wrt_var = "x2")
#>        0      1      2      3      4
#> 0 0.0000 0.3435 0.3183 0.3533 0.3177
#> 1 0.3112 0.0000 0.2616 0.3195 0.2629
#> 2 0.1789 0.1623 0.0000 0.1658 0.1434
#> 3 0.3342 0.3336 0.2790 0.0000 0.2760
#> 4 0.1757 0.1605 0.1412 0.1614 0.0000
```

The rest of the toolkit —
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md),
[`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md),
[`blp()`](https://fpcordeiro.github.io/choicer/reference/blp.md) —
behaves exactly as in the [getting-started
vignette](https://fpcordeiro.github.io/choicer/articles/choicer.md),
integrating over the distribution of tastes automatically.
