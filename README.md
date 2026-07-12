
# choicer: fast discrete-choice models with a focus on economic applications

<!-- badges: start -->
[![CRAN version](https://www.r-pkg.org/badges/version/choicer)](https://CRAN.R-project.org/package=choicer)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/last-month/choicer?color=blue)](https://CRAN.R-project.org/package=choicer)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fpcordeiro/choicer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`choicer` provides fast discrete-choice estimation for applied economics. Multinomial Logit (MNL), Mixed Logit (MXL), and Nested Logit (NL) use C++ likelihoods, analytical gradients and
Hessians, and OpenMP parallelism. Bayesian Multinomial Probit (MNP), Hierarchical Multinomial Logit (HMNL), and Hierarchical Multinomial Probit (HMNP) use dedicated C++ MCMC
kernels and posterior diagnostics. A common research interface covers fitted
shares, elasticities, diversion, willingness to pay, welfare, counterfactuals,
and BLP inversion wherever each model supports the economic object. Frequentist
inference includes analytical-Hessian, BHHH, robust/WESML, and cluster-robust
covariances, recomputable through `vcov()` without refitting.

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
```

`mode_choice` is a **choice-based sample**: minority modes were over-sampled and
car was under-sampled. With the full set of ASCs used here, the MNL slope
coefficients—and therefore the WTP ratios—remain consistent, but the unweighted
constants, fitted shares, elasticity levels, and surplus levels reflect the
sampling design rather than the population. The following outputs demonstrate
the workflow; do not report them as population mode shares or
population-average welfare without
external population shares and WESML weights (see `vignette("wesml")`).

``` r

# Post-estimation
predict(fit, type = "shares")           # fitted sample shares, not population shares
elasticities(fit, elast_var = "vcost")  # own- and cross-price elasticities
diversion_ratios(fit)                   # diversion ratio matrix
wtp(fit, price_var = "vcost")           # willingness-to-pay with delta-method SEs
gof(fit)                                # McFadden R2 and hit rate (also in summary())

# Counterfactual / policy prediction: cut train fares 25% and predict
mc_cf <- mode_choice
mc_cf$vcost[mc_cf$mode == "train"] <- 0.75 * mc_cf$vcost[mc_cf$mode == "train"]
predict(fit, type = "shares", newdata = mc_cf)

# Illustrative sample-average welfare contrast (not a population estimate)
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
predict(fit_nl, type = "shares")        # aggregate fitted shares in this simulation
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
`choicer_hb` post-estimation suite — including entry counterfactuals, where a
never-observed alternative receives a posterior-predictive `delta` from its
characteristics (see `vignette("hb", package = "choicer")`):

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
    chains             = 2,
    mcmc               = list(R = 4000, burn = 1000, thin = 2)
)

summary(fit_hmnl)

# Posterior-aware post-estimation (x2 is the price-like attribute in this DGP)
wtp(fit_hmnl, price_var = "x2")               # posterior median + quantile interval
consumer_surplus(fit_hmnl, price_var = "x2")  # posterior-median CS by task

# Multi-chain diagnostics and posterior-predictive fit
b_chains <- lapply(fit_hmnl$chains, function(ch) ch$b)
rhat(b_chains, rank = TRUE)
ess(b_chains)
mcse(b_chains)
ppc_shares(fit_hmnl)
```

For serious posterior work, use multiple chains, inspect rank-normalized R-hat,
bulk/tail ESS, MCSE and trace plots for all reported blocks, and increase the run
length until the slowest block is stable. choicer's chains use distinct seeds but
the same data-driven initialization, so they are not deliberately overdispersed;
a one-chain split diagnostic is weaker still because it cannot reveal disagreement
between separately simulated runs.

For HMNL/HMNP, v0.2.0 retains every chain in `fit$chains` for diagnostics but
uses chain 1 for the top-level posterior summaries and post-estimation methods;
it does not pool chains automatically.

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
| Hierarchical Bayesian MNP | `run_hmnprobit()` | `predict()`, `elasticities()`, `diversion_ratios()`, `wtp()`, `recovery_table()` |

The frequentist `choicer_fit` models (MNL, MXL, NL) support `summary()`,
`coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, goodness-of-fit,
and counterfactual `predict(newdata = )`. Bayesian `choicer_mnp` and
`choicer_hb` objects are posterior-draws objects: they support `summary()`,
`coef()`, `vcov()`, and `nobs()`, but not likelihood-based AIC/BIC. MNP does
not yet have prediction or economic post-estimation; HMNL and HMNP have their
own posterior prediction and substitution methods.

### Capability map

| Capability | MNL | MXL | NL | MNP | HMNL | HMNP |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| Panel likelihood / persistent person tastes | — | — | — | — | yes | yes |
| Flexible unobserved substitution | iid EV1 benchmark | random tastes | nested shocks | full differenced-normal covariance | panel tastes | panel tastes; iid normal shocks |
| WESML / cluster-robust frequentist inference | yes | yes | yes | — | — | — |
| Counterfactual prediction | yes | yes | yes | — | yes | yes |
| Elasticities and diversion | yes | yes | yes | — | yes | yes |
| Logsum / consumer-surplus welfare | yes | yes | yes | — | yes | — |
| Entry prediction for unseen alternatives | — | — | — | — | yes | yes |
| Native multi-chain posterior diagnostics | — | — | — | — | yes | yes |

## Alternative packages

Several excellent R packages estimate discrete-choice models:

- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [gmnl](https://CRAN.R-project.org/package=gmnl)
- [apollo](https://CRAN.R-project.org/package=apollo)
- [mixl](https://cran.r-project.org/package=mixl)
- [bayesm](https://cran.r-project.org/package=bayesm)

choicer's distinguishing focus is a fast C++ core with analytical gradients and
Hessians, and one consistent post-estimation toolkit aimed at the demand and
welfare quantities applied economists report: elasticities, diversion,
willingness to pay, consumer surplus, counterfactual aggregate shares, and share
inversion, with the matching robust, clustered, and choice-based-sampling
inference.
