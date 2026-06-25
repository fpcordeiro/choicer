# Nested logit and grouped substitution

Nested logit groups alternatives into *nests* of close substitutes.
Within a nest, alternatives share an unobserved component, so
substitution is stronger inside a nest than across nests. Each nest has
a positive dissimilarity parameter λ. The usual random-utility
interpretation of nested logit focuses on `0 < λ <= 1`: λ = 1 gives the
MNL limit, and smaller λ means tighter substitution inside the nest.
choicer imposes λ \> 0 rather than an upper bound. Estimates above one
are mathematically computable but imply negative within-nest correlation
and should be treated as evidence against the proposed nesting
structure, not as the usual “close substitutes” interpretation.

Think of nested logit as the **parsimonious middle ground** between the
multinomial and mixed logits. It introduces genuine within-nest
correlation in unobserved utility at the cost of just one extra
parameter per nest, and it stays globally well-behaved and cheap to
estimate (a closed-form GEV likelihood, no simulation). The price you
pay is a strong prior: *you* must specify the nesting tree in advance,
and the model only permits correlation within the nests you draw. When
the right grouping is obvious from the application (travel modes,
product categories), that is a defensible restriction; when it is not,
the result can be sensitive to how you nest. A good nested-logit
application therefore treats the tree as an economic assumption, not as
a tuning parameter chosen after looking at fit statistics. See [Choosing
among choice
models](https://fpcordeiro.github.io/choicer/articles/choicer.html#choosing-among-choice-models)
for how this tradeoff compares with the alternatives.

``` r

library(choicer)
set_num_threads(2)
```

## Simulate a nested process

[`simulate_nl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_nl_data.md)
builds two nests of inside goods plus an outside option, with known
dissimilarity parameters.

``` r

sim <- simulate_nl_data(N = 4000, seed = 1)
sim
#> <choicer_sim: nl>
#>   settings:
#>     N = 4000
#>     J_inside = 5
#>     nest_structure = (1,2), (3,4,5)
#>   rows in $data: 24000
#>   true_params: beta, delta, lambdas
```

## Fit

Supply the nest membership column; choicer estimates the coefficients,
the alternative-specific constants and the nest dissimilarity parameters
jointly, using an analytical gradient and Hessian.

``` r

fit <- run_nestlogit(
  data                   = sim$data,
  id_col                 = "id",
  alt_col                = "j",
  choice_col             = "choice",
  covariate_cols         = c("X", "W"),
  nest_col               = "nest",
  use_asc                = TRUE,
  include_outside_option = TRUE,
  outside_opt_label      = 0L
)
#> Optimization run time 0h:0m:0.12s
summary(fit)
#> Nested Logit (NL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> X            1.502372   0.047212  31.8220  0.00e+00  ***
#> W           -0.806760   0.027163 -29.7012  0.00e+00  ***
#> Lambda_1     0.794271   0.043175  18.3964  0.00e+00  ***
#> Lambda_2     0.185597   0.013846  13.4041  0.00e+00  ***
#> ASC_1        0.547211   0.123706   4.4235  9.71e-06  ***
#> ASC_2        0.316589   0.128397   2.4657  1.37e-02  *
#> ASC_3       -0.137415   0.133417  -1.0300  3.03e-01  
#> ASC_4       -0.449435   0.140065  -3.2088  1.33e-03  **
#> ASC_5        0.426713   0.122610   3.4802  5.01e-04  ***
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -3063.62 
#> AIC: 6145.24  | BIC: 6201.89 
#> McFadden R2: 0.573 (adj: 0.571) | Hit rate: 0.694 
#> N: 4000  | Parameters: 9 
#> Optimization time: 0.12 s
#> Convergence: 3 ( NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. )
```

## Parameter recovery

The `lambda` rows are the nest dissimilarity parameters — the part that
is unique to nested logit.

``` r

recovery_table(fit, sim$true_params)
#> <choicer_recovery> model=choicer_nl level=0.95
#>    parameter  group  true estimate     se    bias rel_bias_pct z_vs_true
#>       <char> <char> <num>    <num>  <num>   <num>        <num>     <num>
#> 1:         X   beta   1.5   1.5024 0.0472  0.0024       0.1581    0.0502
#> 2:         W   beta  -0.8  -0.8068 0.0272 -0.0068       0.8450   -0.2489
#> 3:  Lambda_1 lambda   0.8   0.7943 0.0432 -0.0057      -0.7161   -0.1327
#> 4:  Lambda_2 lambda   0.2   0.1856 0.0138 -0.0144      -7.2013   -1.0402
#> 5:     ASC_1    asc   0.5   0.5472 0.1237  0.0472       9.4422    0.3816
#> 6:     ASC_2    asc   0.3   0.3166 0.1284  0.0166       5.5296    0.1292
#> 7:     ASC_3    asc  -0.2  -0.1374 0.1334  0.0626     -31.2926    0.4691
#> 8:     ASC_4    asc  -0.5  -0.4494 0.1401  0.0506     -10.1130    0.3610
#> 9:     ASC_5    asc   0.4   0.4267 0.1226  0.0267       6.6782    0.2179
#>    lower_ci upper_ci covers
#>       <num>    <num> <lgcl>
#> 1:   1.4098   1.5949   TRUE
#> 2:  -0.8600  -0.7535   TRUE
#> 3:   0.7096   0.8789   TRUE
#> 4:   0.1585   0.2127   TRUE
#> 5:   0.3048   0.7897   TRUE
#> 6:   0.0649   0.5682   TRUE
#> 7:  -0.3989   0.1241   TRUE
#> 8:  -0.7240  -0.1749   TRUE
#> 9:   0.1864   0.6670   TRUE
```

## Nest-consistent elasticities

This is the payoff of nesting: cross-elasticities are larger for
alternatives in the *same* nest than for alternatives in different
nests. choicer’s
[`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md)
respects the nest structure automatically.

``` r

elasticities(fit, elast_var = "W")
#>   0       1       2        3        4       5
#> 0 0 -0.1118 -0.1097 -0.07679 -0.07381 -0.1138
#> 1 0 -0.1529 -0.1444 -0.07679 -0.07381 -0.1138
#> 2 0 -0.1423 -0.1271 -0.07679 -0.07381 -0.1138
#> 3 0 -0.1118 -0.1097 -0.79069 -0.65893 -0.7997
#> 4 0 -0.1118 -0.1097 -0.63239 -0.66527 -0.7997
#> 5 0 -0.1118 -0.1097 -0.63239 -0.65893 -0.8606
diversion_ratios(fit)
#>        0      1       2       3       4       5
#> 0 0.0000 0.0553 0.05027 0.03971 0.04039 0.04317
#> 1 0.2725 0.0000 0.36377 0.23303 0.22071 0.27854
#> 2 0.2300 0.3379 0.00000 0.20647 0.21384 0.24121
#> 3 0.1533 0.1826 0.17416 0.00000 0.24129 0.24856
#> 4 0.1319 0.1463 0.15263 0.20418 0.00000 0.18852
#> 5 0.2123 0.2780 0.25917 0.31661 0.28378 0.00000
```

Those substitution patterns are credible only to the extent that the
nesting tree is credible. If plausible alternative trees imply
materially different diversion or welfare conclusions, that sensitivity
is part of the empirical result rather than a nuisance to hide.

## Share inversion with BLP

[`blp()`](https://fpcordeiro.github.io/choicer/reference/blp.md) runs
the Berry-Levinsohn-Pakes contraction to recover the mean utilities that
reproduce a set of target shares — useful for calibration and demand
estimation from aggregate data. For strongly-nested models a damping
factor below 1 stabilizes the iteration.

``` r

target_shares <- predict(fit, type = "shares")
head(blp(fit, target_shares, damping = 0.5))
#>         [,1]
#> [1,]  0.5472
#> [2,]  0.3166
#> [3,] -0.1374
#> [4,] -0.4494
#> [5,]  0.4267
```

As always, [`predict()`](https://rdrr.io/r/stats/predict.html),
[`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md) and
[`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md)
are available on the fitted object with the same syntax used throughout
choicer.
