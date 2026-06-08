
# choicer: fast discrete-choice models with a focus on economic applications

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/choicer)](https://CRAN.R-project.org/package=choicer)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`choicer` provides fast estimation of discrete-choice models for applied economics. Likelihoods, analytical gradients and Hessians are implemented in C++ with OpenMP parallelism, scaling efficiently to specifications with many alternative-specific constants. Post-estimation routines return predicted shares, own- and cross-price elasticities, diversion ratios, and the BLP contraction. Supports multinomial logit (MNL), mixed logit (MXL), and nested logit (NL); more models will be added.

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

Estimate a multinomial logit model and compute elasticities and diversion ratios:

``` r
library(choicer)
library(data.table)

# Estimate
fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_col = c("x1", "x2")
)

summary(fit)

# Post-estimation
predict(fit, type = "shares")        # predicted market shares
elasticities(fit, elast_var = "x1")  # own- and cross-price elasticities
diversion_ratios(fit)                # diversion ratio matrix
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

## Supported models

| Model | Function | Post-estimation |
|-------|----------|-----------------|
| Multinomial Logit | `run_mnlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()` |
| Mixed Logit | `run_mxlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()` |
| Nested Logit | `run_nestlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()` |

All fitted models support `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `nobs()`.

## Alternative packages

There are multiple R packages that offer similar functionalities:

- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [gmnl](https://CRAN.R-project.org/package=gmnl)
- [apollo](https://CRAN.R-project.org/package=apollo)
- [mixl](https://cran.r-project.org/package=mixl)
