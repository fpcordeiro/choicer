
# choicer: fast discrete-choice models with a focus on economic applications

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`choicer` provides implementations of discrete-choice models with a focus on economic applications. Computationally intensive likelihoods are written in C++ and exposed for use with generic optimizers. Special care is taken to handle high-dimensional alternative-specific constants efficiently. Post-estimation routines compute elasticities, diversion ratios, and BLP contraction. Currently supports multinomial logit (MNL), mixed logit (MXL), and nested logit (NL); more models will be added.

## Installation

You can install the development version of `choicer` with:

``` r
# Using remotes
remotes::install_github("fpcordeiro/choicer")

# Or using pak
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
| Mixed Logit | `run_mxlogit()` | `elasticities()`, `blp()` |
| Nested Logit | `run_nestlogit()` | — |

All fitted models support `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `nobs()`.

## Alternative packages

There are multiple R packages that offer similar functionalities:

- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [gmnl](https://CRAN.R-project.org/package=gmnl)
- [apollo](https://CRAN.R-project.org/package=apollo)
