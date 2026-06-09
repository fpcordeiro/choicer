# choicer (development version)

## Post-estimation

- `wtp()` — willingness-to-pay for MNL, MXL, and NL with analytic delta-method standard errors. For MXL, log-normal random coefficients report the *median* WTP under the package's shifted log-normal parameterization; random price coefficients are rejected
- `gof()` — goodness of fit: McFadden pseudo R-squared (plain and adjusted) and in-sample hit rate, with `"equal_shares"` (default) and `"market_shares"` null models; now also shown in the `summary()` footer
- `predict(..., newdata = )` — counterfactual prediction for all three models, from either a long data.frame in the fit-time format or a modified-design list (`X`, `alt_idx`, `M`, ...); works even with `keep_data = FALSE`. NL fits now store `nest_idx` top-level to support this
- `logsum()` — expected maximum utility (inclusive value) per choice situation for MNL (closed form), MXL (simulated with a dedicated per-draw kernel, avoiding the Jensen bias of averaging utilities first), and NL (nested inclusive-value formula)
- `consumer_surplus()` — expected consumer surplus `logsum / (-alpha)` (Train 2009, Ch. 3) with a delta-method standard error of the mean CS for MNL; supports `newdata` for policy ΔCS analysis

# choicer 0.1.0

Initial CRAN release.

## Supported models

- **Multinomial Logit** (`run_mnlogit()`) — estimation, prediction, elasticities, diversion ratios, BLP contraction
- **Mixed Logit** (`run_mxlogit()`) — normal and log-normal random coefficients, correlated random coefficients via Cholesky, Halton draws, elasticities, BLP contraction
- **Nested Logit** (`run_nestlogit()`) — estimation with nest-specific dissimilarity (lambda) parameters

## S3 class system

- Parent class `choicer_fit` with subclasses `choicer_mnl`, `choicer_mxl`, `choicer_nl`
- Standard methods: `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `predict()`
- Classed data objects: `choicer_data_mnl`, `choicer_data_mxl`, `choicer_data_nl` from `prepare_*_data()`

## Post-estimation generics

- `elasticities()` — methods for MNL and MXL
- `diversion_ratios()` — method for MNL
- `blp()` — BLP contraction for MNL and MXL

## API

- Dual workflow for all `run_*logit()` functions: convenience (pass `data` + column names) or advanced (pass pre-prepared `input_data`)
- Pluggable optimizer via `optimizer = "nloptr" | "optim" | <function>`
- `prepare_nl_data()` added for nested logit data preparation

## Computation

- C++ likelihoods, gradients, and analytical Hessians via Rcpp/RcppArmadillo
- OpenMP parallelization over individuals
- Log-sum-exp trick throughout for numerical stability
