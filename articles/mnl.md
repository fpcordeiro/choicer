# Multinomial logit

The multinomial logit (MNL) is the workhorse of discrete choice. choicer
fits it by maximum likelihood with a C++ core, analytical gradients and
an analytical Hessian, so estimation and standard errors are fast even
with many alternative-specific constants.

This vignette does two things: it shows the full MNL workflow, and —
because the data are simulated from a known process — it checks that
choicer **recovers the true parameters**.

``` r

library(choicer)
set_num_threads(2)
```

## Simulate a known data-generating process

[`simulate_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnl_data.md)
draws choices from a logit model with i.i.d. Gumbel errors. The returned
object carries both the data and the true parameters.

``` r

sim <- simulate_mnl_data(N = 2000, J = 4, seed = 1)
sim
#> <choicer_sim: mnl>
#>   settings:
#>     N = 2000
#>     J = 4
#>     K_x = 2
#>     outside_option = TRUE
#>     vary_choice_set = TRUE
#>   rows in $data: 7967
#>   true_params: beta, delta
```

## Fit

``` r

fit <- run_mnlogit(
  data           = sim$data,
  id_col         = "id",
  alt_col        = "alt",
  choice_col     = "choice",
  covariate_cols = c("x1", "x2")
)
#> Optimization run time 0h:0m:0.05s
summary(fit)
#> Multinomial Logit (MNL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           0.790385   0.057681  13.7026  0.00e+00  ***
#> x2          -0.568065   0.056300 -10.0900  0.00e+00  ***
#> ASC_1        0.459936   0.067736   6.7902  1.12e-11  ***
#> ASC_2       -0.570242   0.085511  -6.6686  2.58e-11  ***
#> ASC_3        0.488752   0.067826   7.2060  5.76e-13  ***
#> ASC_4       -0.467326   0.081197  -5.7555  8.64e-09  ***
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -2437.1 
#> AIC: 4886.2  | BIC: 4919.8 
#> McFadden R2: 0.104 (adj: 0.102) | Hit rate: 0.447 
#> N: 2000  | Parameters: 6 
#> Optimization time: 0.05 s
#> Convergence: 1 ( NLOPT_SUCCESS: Generic success return value. )
```

## Did we recover the truth?

[`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md)
lines up each estimate against the value that generated the data, with
the bias, a z-score, and whether the 95% confidence interval covers the
truth.

``` r

recovery_table(fit, sim$true_params)
#> <choicer_recovery> model=choicer_mnl level=0.95
#>    parameter  group  true estimate     se    bias rel_bias_pct z_vs_true
#>       <char> <char> <num>    <num>  <num>   <num>        <num>     <num>
#> 1:        x1   beta   0.8   0.7904 0.0577 -0.0096       -1.202   -0.1667
#> 2:        x2   beta  -0.6  -0.5681 0.0563  0.0319       -5.322    0.5672
#> 3:     ASC_1    asc   0.5   0.4599 0.0677 -0.0401       -8.013   -0.5915
#> 4:     ASC_2    asc  -0.5  -0.5702 0.0855 -0.0702       14.048   -0.8214
#> 5:     ASC_3    asc   0.5   0.4888 0.0678 -0.0112       -2.249   -0.1658
#> 6:     ASC_4    asc  -0.5  -0.4673 0.0812  0.0327       -6.535    0.4024
#>    lower_ci upper_ci covers
#>       <num>    <num> <lgcl>
#> 1:   0.6773   0.9034   TRUE
#> 2:  -0.6784  -0.4577   TRUE
#> 3:   0.3272   0.5927   TRUE
#> 4:  -0.7378  -0.4026   TRUE
#> 5:   0.3558   0.6217   TRUE
#> 6:  -0.6265  -0.3082   TRUE
```

In this run the intervals cover their true values, which is indicative
of correct behavior. What matters formally, though, is coverage over
repeated simulations; a full Monte Carlo exercise is left outside this
vignette for brevity.

## Post-estimation

The full demand and welfare toolkit is available on the fitted object.
Treating `x2` as price:

``` r

predict(fit, type = "shares")        # predicted market shares
#>        [,1]
#> [1,] 0.2405
#> [2,] 0.2610
#> [3,] 0.1080
#> [4,] 0.2645
#> [5,] 0.1260
elasticities(fit, elast_var = "x2")  # own- and cross-price elasticities
#>   0        1         2         3         4
#> 0 0 -0.01624 -0.010341 -0.013943 -0.010855
#> 1 0 -0.01633 -0.006815 -0.008984 -0.006413
#> 2 0 -0.01336 -0.002613 -0.010000 -0.008505
#> 3 0 -0.01081 -0.006343 -0.020585 -0.006275
#> 4 0 -0.01114 -0.007366 -0.010306 -0.010266
diversion_ratios(fit)                # where demand goes
#>        0      1      2      3      4
#> 0 0.0000 0.3627 0.3252 0.3652 0.3403
#> 1 0.3227 0.0000 0.2684 0.3245 0.2587
#> 2 0.1620 0.1503 0.0000 0.1417 0.1364
#> 3 0.3233 0.3229 0.2519 0.0000 0.2646
#> 4 0.1920 0.1641 0.1545 0.1686 0.0000
wtp(fit, price_var = "x2")           # willingness to pay, with delta-method SEs
#> Willingness to pay (WTP), price variable: 'x2' (95% CI)
#>    Estimate Std_Error z_value CI_lower CI_upper
#> x1    1.391    0.1643   8.467    1.069    1.713
gof(fit)                             # McFadden R2 and hit rate
#> Goodness of fit
#>   Log-likelihood:       -2437.1 
#>   Null log-likelihood: -2721.21 (equal shares)
#>   McFadden R2:          0.1044 (adj: 0.1022)
#>   Hit rate:             0.4465
#>   N: 2000  | Parameters: 6
```

## A note on IIA

It is worth being precise about what IIA is and where it bites, because
the usual one-liner (“logit assumes IIA, so it’s wrong”) obscures more
than it explains.

**In micro MNL, IIA is imposed on individual conditional choice
probabilities.** IIA can also be defined for aggregate shares or other
choice-probability systems, but the logit restriction starts at the
individual level. For a single decision-maker, the ratio of any two
choice probabilities, $`P_{ij}/P_{ik} = \exp(V_{ij} - V_{ik})`$, depends
only on the attributes of $`j`$ and $`k`$ — never on what other
alternatives exist or what they look like. That is the individual-level
independence-of-irrelevant-alternatives restriction, and the MNL
satisfies it by construction.

**At the aggregate level the picture is different.** The diversion
ratios and cross-elasticities reported above are sample averages over
individuals. The familiar share-proportional formula
$`DR(j\to k) = s_k / (1 - s_j)`$ holds *only* when every individual
faces the same covariates, so that $`P_{ij}`$ does not vary across
$`i`$. As soon as covariates vary across choice situations — different
attributes, demographics, or choice sets — the aggregate substitution
patterns are **no longer share-proportional** and need not satisfy a
market-level IIA property, even though each individual conditional
probability still obeys IIA. (The diversion derivation in the [math
companion](https://fpcordeiro.github.io/choicer/articles/multinomial_logit_math.md)
makes this explicit.)

So the MNL is more flexible in aggregate than “IIA” suggests — but only
along dimensions you have *observed*. Its genuine limitations are two.
First, all of its substitution flexibility must be carried by observed
covariates; it cannot represent substitution driven by *unobserved*
correlation in tastes. Second, for a single decision-maker — and
therefore for any individual-level counterfactual — substitution remains
rigidly proportional. That second point matters most for counterfactuals
(a new product, a removed alternative), where there is no data to
discipline the substitution pattern and the model’s structure does all
the work.

That is why the MNL is best treated as a disciplined baseline rather
than a straw man. If the variation in the data is rich enough to explain
the substitution margin through observed attributes, the MNL can be the
most transparent model. If the economics turns on unobserved closeness —
products that share latent tastes, schools with unmeasured peer groups,
hospitals with unobserved networks — a more structured model is doing
substantive work, and that work needs to be defended.

When those limitations bite — close substitutes, rich unobserved
heterogeneity, counterfactuals that hinge on *who* substitutes to *what*
— reach for the [nested
logit](https://fpcordeiro.github.io/choicer/articles/nl.md) or [mixed
logit](https://fpcordeiro.github.io/choicer/articles/mxl.md). They keep
the same [`predict()`](https://rdrr.io/r/stats/predict.html) /
[`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md)
/ [`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md)
interface, and the [getting-started
vignette](https://fpcordeiro.github.io/choicer/articles/choicer.html#choosing-among-logit-models)
lays out what each relaxation costs.
