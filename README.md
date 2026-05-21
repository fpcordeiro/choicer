
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

## Supported models

| Model | Function | Post-estimation |
|-------|----------|-----------------|
| Multinomial Logit | `run_mnlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()` |
| Mixed Logit | `run_mxlogit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `blp()` |
| Nested Logit | `run_nestlogit()` | — |

All fitted models support `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `nobs()`.

## Alternative packages

There are multiple R packages that offer similar functionalities:

- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [gmnl](https://CRAN.R-project.org/package=gmnl)
- [apollo](https://CRAN.R-project.org/package=apollo)
- [mixl](https://cran.r-project.org/package=mixl)
