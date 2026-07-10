
# choicer: fast discrete-choice models with a focus on economic applications

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/choicer)](https://CRAN.R-project.org/package=choicer)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`choicer` provides fast estimation of discrete-choice models for applied economics. Likelihoods, analytical gradients and Hessians are implemented in C++ with OpenMP parallelism, scaling efficiently to specifications with many alternative-specific constants. Post-estimation routines return predicted shares, own- and cross-price elasticities, diversion ratios, willingness-to-pay with delta-method standard errors, goodness-of-fit measures, counterfactual predictions, and the BLP contraction. Supports multinomial logit (MNL), mixed logit (MXL), and nested logit (NL), plus Bayesian multinomial probit (MNP) and hierarchical Bayesian multinomial logit and probit (HMNL, HMNP) via Gibbs sampling; more models will be added.

## Installation

Install the released version from CRAN:

``` r
install.packages("choicer")
```

Or install the development version from GitHub:

``` r
pak::pkg_install("fpcordeiro/choicer")
```

## Example

Estimate a multinomial logit model on the shipped `mode_choice` data (the
classic Greene & Hensher intercity travel-mode study: 210 travellers choosing
among air, train, bus, and car) and compute elasticities and diversion ratios:

``` r
library(choicer)
data(mode_choice)

# Estimate
fit <- run_mnlogit(
    data           = mode_choice,
    id_col         = "id",
    alt_col        = "mode",
    choice_col     = "choice",
    covariate_cols = c("wait", "travel", "vcost")
)

summary(fit)

# Post-estimation
predict(fit, type = "shares")           # predicted market shares
elasticities(fit, elast_var = "vcost")  # own- and cross-price elasticities
diversion_ratios(fit)                   # diversion ratio matrix
wtp(fit, price_var = "vcost")           # willingness-to-pay with delta-method SEs
gof(fit)                                # McFadden R2 and hit rate (also in summary())

# Counterfactual / policy prediction: cut train fares 25% and predict
mc_cf <- mode_choice
mc_cf$vcost[mc_cf$mode == "train"] <- 0.75 * mc_cf$vcost[mc_cf$mode == "train"]
predict(fit, type = "shares", newdata = mc_cf)

# Welfare: change in expected consumer surplus from the counterfactual
cs0 <- consumer_surplus(fit, price_var = "vcost")
cs1 <- consumer_surplus(fit, price_var = "vcost", newdata = mc_cf)
cs1$mean_cs - cs0$mean_cs
```

The same post-estimation toolkit is available for nested logit. Elasticities
respect the nest structure — within-nest and cross-nest cross-elasticities
differ, so IIA does not hold across nests. The `blp()` contraction for NL
accepts a `damping` argument (default 1) which can be reduced for models with
strong nesting:

``` r
# Nested logit — simulate data and fit
sim <- simulate_nl_data(N = 5e4, seed = 123)

fit_nl <- run_nestlogit(
    data                   = sim$data,
    id_col                 = "id",
    alt_col                = "j",
    choice_col             = "choice",
    covariate_cols         = c("X", "W"),
    nest_col               = "nest",
    use_asc                = TRUE,
    include_outside_option = TRUE,
    outside_opt_label      = 0L,
    keep_data              = TRUE   # required for post-estimation
)

# Post-estimation
predict(fit_nl, type = "shares")        # predicted market shares
elasticities(fit_nl, elast_var = "X")   # J×J elasticity matrix (nest-consistent)
diversion_ratios(fit_nl)                # J×J diversion matrix
# BLP share inversion: recover mean utilities matching target shares
target_shares <- predict(fit_nl, type = "shares")
blp(fit_nl, target_shares, damping = 0.5)  # use damping < 1 for strongly-nested models
```

`choicer` also fits hierarchical Bayesian models with a BLP-style random-effect
alternative intercept `delta_j = z_j'theta + xi_j` on top of respondent-level
random coefficients `beta_i ~ N(b, W)`, estimated by Gibbs sampling. Both
models share an implicit outside option (no base alternative) and the
`choicer_hb` post-estimation suite:

``` r
# Hierarchical Bayesian MNL — simulate panel data and fit
sim_hb <- simulate_hmnl_data(N = 300, T = 6, J = 4, seed = 123)

set.seed(42)
fit_hmnl <- run_hmnlogit(
    data               = sim_hb$data,
    id_col             = "task",
    alt_col            = "alt",
    choice_col         = "choice",
    covariate_cols     = c("x1", "x2"),
    person_col         = "pid",
    alt_covariate_cols = "z1",
    mcmc               = list(R = 4000, burn = 1000, thin = 2)
)

summary(fit_hmnl)

# Posterior-aware post-estimation
wtp(fit_hmnl, price_var = "x1")               # posterior median + quantile interval
consumer_surplus(fit_hmnl, price_var = "x1")  # expected consumer surplus

# MCMC diagnostics
rhat(fit_hmnl$draws$b)       # split R-hat for the population-mean draws
ppc_shares(fit_hmnl)         # observed vs. posterior-predictive choice shares
```

The Bayesian multinomial probit `run_mnprobit()` (non-hierarchical, class
`choicer_mnp`) is fit the same way, via `prior=`/`mcmc=` instead of an
optimizer; see `?run_mnprobit`.

## Supported models

| Model | Function | Post-estimation |
|-------|----------|-----------------|
| Multinomial Logit | `run_mnlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()`, `wtp()`, `gof()`, `logsum()`, `consumer_surplus()` |
| Mixed Logit | `run_mxlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()`, `wtp()`, `gof()`, `logsum()`, `consumer_surplus()` |
| Nested Logit | `run_nestlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()`, `wtp()`, `gof()`, `logsum()`, `consumer_surplus()` |
| Bayesian MNP | `run_mnprobit()` | `summary()`, `coef()`, `vcov()`, `recovery_table()` |
| Hierarchical Bayesian MNL | `run_hmnlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `wtp()`, `logsum()`, `consumer_surplus()`, `recovery_table()` |
| Hierarchical Bayesian MNP | `run_hmnprobit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `wtp()`, `logsum()`, `consumer_surplus()`, `recovery_table()` |

All fitted models support `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `nobs()`. `summary()` reports McFadden R2 and the hit rate alongside the usual fit statistics, and `predict()` accepts `newdata` (a long data.frame or a modified design list) for counterfactual and policy prediction, even on fits with `keep_data = FALSE`. The Bayesian models (`choicer_mnp`, and `choicer_hmnl`/`choicer_hmnp` under the shared `choicer_hb` class) are posterior-draws objects rather than `choicer_fit` models: they support `summary()`, `coef()`, `vcov()`, and `nobs()`, but not `logLik()`, `AIC()`, or `BIC()`.

## Alternative packages

There are multiple R packages that offer similar functionalities:

- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [gmnl](https://CRAN.R-project.org/package=gmnl)
- [apollo](https://CRAN.R-project.org/package=apollo)
- [mixl](https://cran.r-project.org/package=mixl)
