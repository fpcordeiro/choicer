# choicer (development version)

## Mixed logit — on-the-fly randomized Halton draws

- `run_mxlogit()` gains three new arguments: `draws`, `seed`, and `scramble`.
  - `draws = "store"` (default) keeps the existing behavior: a full
    K_w × S × N Halton cube is pre-materialized and stored in memory.
  - `draws = "generate"` activates the new mode: each individual's S draws are
    computed on the fly in C++ from a compact seed, eliminating the O(N)
    `eta_draws` cube. This is the recommended choice when N is large or memory
    is constrained.
  - `scramble = "owen"` (default when `draws = "generate"`) applies Owen (2017)
    digit scrambling to the Halton sequence for better high-dimensional
    uniformity; recommended for K_w > 5 (Bhat 2003). `scramble = "none"`
    reproduces the randtoolbox sequence exactly and is intended for testing.
  - `seed` sets the integer master seed for the on-the-fly generator; if
    `NULL` (default), a seed is drawn from R's RNG so `set.seed()` governs
    reproducibility. Ignored when `draws = "store"`.
- The `draws_info` field on `choicer_mxl` objects gains three new elements:
  `mode` (`"store"` or `"generate"`), `seed`, and `scramble`. Existing code
  that accesses `draws_info$S`, `draws_info$N`, or `draws_info$K_w` is
  unaffected; the new fields are `NULL` for objects fitted before this release.
- All post-estimation generics (`predict()`, `elasticities()`,
  `diversion_ratios()`, `logsum()`, `consumer_surplus()`, `vcov()`) propagate
  the fitted draw mode automatically — no user action required.
- Default behavior (`draws = "store"`) is bitwise unchanged.


## Bayesian models

- `run_mnprobit()` — Bayesian multinomial probit via Gibbs sampling with data augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). Runs the non-identified chain with conjugate priors and reports identified quantities normalized per draw by `sigma_11`. New `choicer_mnp` posterior object with `summary()` (posterior mean, SD, credible intervals), `coef()`, `vcov()`, `nobs()`; math note in `docs/bayesian_multinomial_probit_math.md`
- C++ MCMC infrastructure built from scratch (`src/rng.h`, `src/bayes_samplers.h`): xoshiro256++/splitmix64 RNG with one stream per (iteration, observation) — draws are bitwise reproducible independent of the OpenMP thread count — plus exact truncated-normal, multivariate-normal, Wishart (Bartlett), and inverse-Wishart samplers. The truncated normal picks the cheapest exact method per region (naive normal rejection in high-mass regions, Robert 1995 exponential rejection in the tail, inverse CDF for narrow intervals)
- The Gibbs chain runs inside a single persistent OpenMP region: the latent-utility sweep and mean refresh are work-shared across choice situations, the conjugate beta/Sigma conditionals run on the master thread between lightweight barriers with hand-rolled fixed-order linear algebra (no BLAS inside the region), and the truncated-normal conditional moments are hoisted per Sigma draw. On the `_benchmarks/` MNP preset this samples ~2x faster than `MNP::mnp()` and ~3.5x faster than `bayesm::rmnpGibbs()` on one thread, and scales with threads on top
- `simulate_mnp_data()` — probit DGP returning a `choicer_sim` with truth on the identified scale; `recovery_table()` gains a `choicer_mnp` method (posterior mean/SD, normal-approximation credible intervals, `sigma` block); parameter-recovery walkthrough in `inst/simulations/mnp_simulation.R`

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
