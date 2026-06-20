# Discrete choice from data to policy, in a dozen lines

**choicer** estimates discrete-choice models for applied economics. The
likelihoods, analytical gradients and Hessians are written in C++ with
OpenMP, so estimation is fast and scales to specifications with hundreds
of alternative-specific constants. Just as important, *one consistent
post-estimation interface* takes you from a fitted model to the
quantities you actually report: market shares, elasticities, diversion
ratios, willingness-to-pay with standard errors, and the welfare effects
of a policy counterfactual.

This vignette runs the whole workflow on a real data set in a few lines.

``` r

library(choicer)
set_num_threads(2) # be a good CRAN citizen; raise this on your own machine
```

## The data

`mode_choice` is the classic Greene & Hensher intercity travel-mode
data: 210 travellers, each choosing among **air, train, bus and car**.
It ships with the package in choicer’s long layout — one row per
traveller and alternative.

``` r

data(mode_choice)
head(mode_choice, 4)
#>   id  mode choice wait travel vcost gcost income size
#> 1  1   air      0   69    100    59    70     35    1
#> 2  1 train      0   34    372    31    71     35    1
#> 3  1   bus      0   35    417    25    70     35    1
#> 4  1   car      1    0    180    10    30     35    1
```

`wait` (terminal waiting time) and `travel` (in-vehicle time) are in
minutes; `vcost` is the travel cost. These vary across modes within a
traveller, which is exactly what a choice model needs.

## Fit a multinomial logit

Point at the identifier, alternative and choice columns, list the
covariates, and fit. Alternative-specific constants are included by
default.

``` r

fit <- run_mnlogit(
  data           = mode_choice,
  id_col         = "id",
  alt_col        = "mode",
  choice_col     = "choice",
  covariate_cols = c("wait", "travel", "vcost")
)
#> Optimization run time 0h:0m:0.01s
summary(fit)
#> Multinomial Logit (MNL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> wait        -0.096887   0.010342  -9.3683  0.00e+00  ***
#> travel      -0.003995   0.000849  -4.7043  2.55e-06  ***
#> vcost       -0.013912   0.006651  -2.0916  3.65e-02  *
#> ASC_train   -0.786669   0.602607  -1.3054  1.92e-01  
#> ASC_bus     -1.433640   0.680713  -2.1061  3.52e-02  *
#> ASC_car     -4.739865   0.867532  -5.4636  4.67e-08  ***
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Log-likelihood: -192.889 
#> AIC: 397.777  | BIC: 417.86 
#> McFadden R2: 0.337 (adj: 0.317) | Hit rate: 0.738 
#> N: 210  | Parameters: 6 
#> Optimization time: 0.01 s
#> Convergence: 3 ( NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. )
```

Estimation takes a fraction of a second. Travellers dislike waiting
time, in-vehicle time and cost — all three coefficients are negative —
and the [`summary()`](https://rdrr.io/r/base/summary.html) footer
reports McFadden’s R² and the in-sample hit rate alongside the usual fit
statistics.

## Predicted shares

``` r

shares <- drop(predict(fit, type = "shares"))
names(shares) <- as.character(fit$alt_mapping[[2]])
round(shares, 3)
#>   air train   bus   car 
#> 0.276 0.300 0.143 0.281

barplot(shares, col = "steelblue", ylab = "Predicted share",
        main = "Predicted mode shares")
```

![](choicer_files/figure-html/shares-1.png)

## Elasticities and diversion

How does demand respond to cost, and where does it go? The same fitted
object answers both, with no extra bookkeeping.

``` r

elasticities(fit, elast_var = "vcost")
#>           air   train      bus      car
#> air   -0.8309  0.1679  0.06729  0.06912
#> train  0.3551 -0.5463  0.06729  0.06912
#> bus    0.3551  0.1679 -0.39815  0.06912
#> car    0.3551  0.1679  0.06729 -0.22296
diversion_ratios(fit)
#>          air  train    bus    car
#> air   0.0000 0.2776 0.2513 0.3860
#> train 0.3080 0.0000 0.3556 0.4144
#> bus   0.1719 0.2192 0.0000 0.1996
#> car   0.5201 0.5032 0.3931 0.0000
```

The own-cost elasticities sit on the diagonal; off-diagonal entries are
the cross-elasticities. Because this is a plain logit, substitution
follows the independence-of-irrelevant-alternatives (IIA) property —
diversion is share-proportional. (Relaxing IIA is exactly what the
[nested logit](https://fpcordeiro.github.io/choicer/articles/nl.md) and
[mixed logit](https://fpcordeiro.github.io/choicer/articles/mxl.md)
vignettes are about.)

## Willingness to pay

With `vcost` playing the role of price,
[`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md) returns
the marginal value of each attribute in money units, with delta-method
standard errors.

``` r

wtp(fit, price_var = "vcost")
#> Willingness to pay (WTP), price variable: 'vcost' (95% CI)
#>        Estimate Std_Error z_value CI_lower  CI_upper
#> wait    -6.9645    3.4085  -2.043 -13.6450 -0.283902
#> travel  -0.2871    0.1436  -2.000  -0.5685 -0.005757
```

These are the travellers’ implied **value of time**: how much they would
pay to shave a minute of in-vehicle or terminal time.

## A policy counterfactual

Suppose a subsidy cuts the cost of **train** travel by 25%. Perturb the
data and predict — no refitting required.

``` r

mc_cf <- mode_choice
mc_cf$vcost[mc_cf$mode == "train"] <- mc_cf$vcost[mc_cf$mode == "train"] * 0.75

shares_cf <- drop(predict(fit, type = "shares", newdata = mc_cf))
names(shares_cf) <- names(shares)

rbind(baseline = shares, counterfactual = shares_cf) |> round(3)
#>                  air train   bus   car
#> baseline       0.276 0.300 0.143 0.281
#> counterfactual 0.269 0.324 0.138 0.270
```

Train’s share rises, drawing travellers away from the other modes. We
can put a money value on the change in welfare with the expected
consumer surplus:

``` r

cs0 <- consumer_surplus(fit, price_var = "vcost")
cs1 <- consumer_surplus(fit, price_var = "vcost", newdata = mc_cf)

cs1$mean_cs - cs0$mean_cs # change in mean consumer surplus
#> [1] 3.2
```

The cheaper train fare raises expected consumer surplus, as it should.

## One interface, every model

Everything above — [`predict()`](https://rdrr.io/r/stats/predict.html),
[`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md),
[`diversion_ratios()`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.md),
[`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md),
[`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md),
plus [`blp()`](https://fpcordeiro.github.io/choicer/reference/blp.md)
for share inversion — works identically on **mixed logit**
([`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md))
and **nested logit**
([`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md))
fits, and choicer also offers a Bayesian **multinomial probit** sampler
([`run_mnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md)).
Learn each model in its own vignette:

- [Multinomial
  logit](https://fpcordeiro.github.io/choicer/articles/mnl.md) — the
  workhorse, and a parameter-recovery check
- [Mixed logit](https://fpcordeiro.github.io/choicer/articles/mxl.md) —
  random coefficients and preference heterogeneity
- [Nested logit](https://fpcordeiro.github.io/choicer/articles/nl.md) —
  correlated alternatives and relaxed IIA
- [Bayesian multinomial
  probit](https://fpcordeiro.github.io/choicer/articles/mnp.md) —
  flexible error covariance via MCMC

## How choicer compares

Several excellent R packages estimate discrete-choice models — among
them [mlogit](https://CRAN.R-project.org/package=mlogit),
[logitr](https://CRAN.R-project.org/package=logitr),
[gmnl](https://CRAN.R-project.org/package=gmnl),
[apollo](https://CRAN.R-project.org/package=apollo) and
[mixl](https://CRAN.R-project.org/package=mixl). choicer’s focus is a
fast C++ core with analytical gradients and Hessians, and a single,
consistent post-estimation toolkit aimed at the demand and welfare
quantities applied economists report.
