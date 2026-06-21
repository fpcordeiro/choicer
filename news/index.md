# Changelog

## choicer (development version)

### WESML sandwich inference for MNL and nested logit

- [`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md)
  and
  [`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md)
  now compute the robust (Huber–White / WESML) sandwich variance
  `V = A^{-1} B A^{-1}`, at full feature parity with
  [`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md).
  Both gain `se_method = "sandwich"` (bread = weighted negated Hessian,
  meat = weight-squared OPG) and `se_method = "bhhh"` (ordinary
  outer-product-of-gradients).
  [`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md)’s
  `se_method` therefore now accepts `"hessian"` (default), `"bhhh"`, or
  `"sandwich"`;
  [`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md)’s
  accepts `"hessian"` (default), `"numeric"`, `"bhhh"`, or `"sandwich"`.
- New C++ kernels
  [`mnl_bhhh_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_bhhh_parallel.md)
  and
  [`nl_bhhh_parallel()`](https://fpcordeiro.github.io/choicer/reference/nl_bhhh_parallel.md)
  accumulate the weighted outer product of per-individual scores. Their
  per-individual score is weight-free, so passing `weights = w` yields
  the BHHH information and `weights = w^2` yields the sandwich meat. The
  NL kernel includes the full beta/lambda/delta score blocks
  (singleton-nest lambdas fixed to 1 contribute no score).
- [`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md),
  [`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md),
  [`prepare_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
  and
  [`prepare_nl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md)
  gain a `weights_col` argument: a row-level weight column is collapsed
  to one weight per choice situation (validated constant within `id`). A
  choice-based-sampling provenance guard auto-adopts the recorded WESML
  weight column and errors rather than silently fitting unweighted under
  a WESML label; non-uniform weights under a non-sandwich `se_method`
  emit a warning.
- Behavior change: weighted fits now emit a steering warning
  recommending `se_method = "sandwich"` when non-uniform weights are
  supplied with the default (`"hessian"`) or `"bhhh"` method; the
  `"bhhh"` case gets a sharper message explaining that BHHH/OPG is not a
  valid WESML correction (its meat is `w^1`, not `w^2`). Point estimates
  and standard errors are unchanged.
- [`wesml_vcov()`](https://fpcordeiro.github.io/choicer/reference/wesml_vcov.md)
  now dispatches on `choicer_mnl` and `choicer_nl` (in addition to
  `choicer_mxl`), returning the post-hoc sandwich variance from a fit
  stored with `keep_data = TRUE`.
- [`summary()`](https://rdrr.io/r/base/summary.html) for MNL and NL fits
  now reports the standard-error method and any WESML weighting in the
  printed footer.
- Added “Choice-Based Sampling and WESML Weighting” sections to the
  multinomial logit and nested logit derivation vignettes.

### Weighting safety hardening

- Weights are now validated to be finite and strictly positive in
  [`prepare_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md),
  [`prepare_nl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md),
  and
  [`prepare_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mxl_data.md)
  (covering both the `weights=` and `weights_col=` paths). Zero,
  negative, or non-finite weights previously could slip through and
  silently invalidate weighted and WESML sandwich inference (weight `w`
  enters the bread, `w^2` the meat); they now error with an actionable
  message.
- Advanced-mode fits (`input_data=` passed directly to
  [`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md),
  [`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md),
  or
  [`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md))
  that carry WESML `choice_sampling` provenance but resolve to uniform
  weights now error instead of warning. The message explains how to
  proceed: bake the non-uniform WESML weights into `input_data` via
  `prepare_*_data(weights=/weights_col=)`, or strip the provenance with
  `attr(input_data, "choice_sampling") <- NULL` for a deliberate
  unweighted fit. Convenience-mode behavior is unchanged.

### Documentation

- Added five vignettes: a getting-started tour (“Discrete choice from
  data to policy, in a dozen lines”) plus one per model (multinomial
  logit, mixed logit, nested logit, Bayesian multinomial probit).
- Added the `mode_choice` data set: the classic Greene & Hensher
  intercity travel-mode choice data (210 travellers x 4 modes), in
  choicer’s long layout, used by the getting-started vignette.
- Added a pkgdown website configuration, including the model derivation
  notes as “The math behind choicer” articles.

### Mixed logit — on-the-fly randomized Halton draws

- [`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)
  gains three new arguments: `draws`, `seed`, and `scramble`.
  - `draws = "store"` (default) keeps the existing behavior: a full K_w
    × S × N Halton cube is pre-materialized and stored in memory.
  - `draws = "generate"` activates the new mode: each individual’s S
    draws are computed on the fly in C++ from a compact seed,
    eliminating the O(N) `eta_draws` cube. This is the recommended
    choice when N is large or memory is constrained.
  - `scramble = "owen"` (default when `draws = "generate"`) applies
    Owen (2017) digit scrambling to the Halton sequence for better
    high-dimensional uniformity; recommended for K_w \> 5 (Bhat 2003).
    `scramble = "none"` reproduces the randtoolbox sequence exactly and
    is intended for testing.
  - `seed` sets the integer master seed for the on-the-fly generator; if
    `NULL` (default), a seed is drawn from R’s RNG so
    [`set.seed()`](https://rdrr.io/r/base/Random.html) governs
    reproducibility. Ignored when `draws = "store"`.
- The `draws_info` field on `choicer_mxl` objects gains three new
  elements: `mode` (`"store"` or `"generate"`), `seed`, and `scramble`.
  Existing code that accesses `draws_info$S`, `draws_info$N`, or
  `draws_info$K_w` is unaffected; the new fields are `NULL` for objects
  fitted before this release.
- All post-estimation generics
  ([`predict()`](https://rdrr.io/r/stats/predict.html),
  [`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md),
  [`diversion_ratios()`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.md),
  [`logsum()`](https://fpcordeiro.github.io/choicer/reference/logsum.md),
  [`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)) propagate the fitted
  draw mode automatically — no user action required.
- Default behavior (`draws = "store"`) is bitwise unchanged.

### Bayesian models

- [`run_mnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md)
  — Bayesian multinomial probit via Gibbs sampling with data
  augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). Runs the
  non-identified chain with conjugate priors and reports identified
  quantities normalized per draw by `sigma_11`. New `choicer_mnp`
  posterior object with
  [`summary()`](https://rdrr.io/r/base/summary.html) (posterior mean,
  SD, credible intervals),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html); math note in
  `docs/bayesian_multinomial_probit_math.md`
- C++ MCMC infrastructure built from scratch (`src/rng.h`,
  `src/bayes_samplers.h`): xoshiro256++/splitmix64 RNG with one stream
  per (iteration, observation) — draws are bitwise reproducible
  independent of the OpenMP thread count — plus exact truncated-normal,
  multivariate-normal, Wishart (Bartlett), and inverse-Wishart samplers.
  The truncated normal picks the cheapest exact method per region (naive
  normal rejection in high-mass regions, Robert 1995 exponential
  rejection in the tail, inverse CDF for narrow intervals)
- The Gibbs chain runs inside a single persistent OpenMP region: the
  latent-utility sweep and mean refresh are work-shared across choice
  situations, the conjugate beta/Sigma conditionals run on the master
  thread between lightweight barriers with hand-rolled fixed-order
  linear algebra (no BLAS inside the region), and the truncated-normal
  conditional moments are hoisted per Sigma draw. On the `_benchmarks/`
  MNP preset this samples ~2x faster than `MNP::mnp()` and ~3.5x faster
  than `bayesm::rmnpGibbs()` on one thread, and scales with threads on
  top
- [`simulate_mnp_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnp_data.md)
  — probit DGP returning a `choicer_sim` with truth on the identified
  scale;
  [`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md)
  gains a `choicer_mnp` method (posterior mean/SD, normal-approximation
  credible intervals, `sigma` block); parameter-recovery walkthrough in
  `inst/simulations/mnp_simulation.R`

### Post-estimation

- [`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md) —
  willingness-to-pay for MNL, MXL, and NL with analytic delta-method
  standard errors. For MXL, log-normal random coefficients report the
  *median* WTP under the package’s shifted log-normal parameterization;
  random price coefficients are rejected
- [`gof()`](https://fpcordeiro.github.io/choicer/reference/gof.md) —
  goodness of fit: McFadden pseudo R-squared (plain and adjusted) and
  in-sample hit rate, with `"equal_shares"` (default) and
  `"market_shares"` null models; now also shown in the
  [`summary()`](https://rdrr.io/r/base/summary.html) footer
- `predict(..., newdata = )` — counterfactual prediction for all three
  models, from either a long data.frame in the fit-time format or a
  modified-design list (`X`, `alt_idx`, `M`, …); works even with
  `keep_data = FALSE`. NL fits now store `nest_idx` top-level to support
  this
- [`logsum()`](https://fpcordeiro.github.io/choicer/reference/logsum.md)
  — expected maximum utility (inclusive value) per choice situation for
  MNL (closed form), MXL (simulated with a dedicated per-draw kernel,
  avoiding the Jensen bias of averaging utilities first), and NL (nested
  inclusive-value formula)
- [`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md)
  — expected consumer surplus `logsum / (-alpha)` (Train 2009, Ch. 3)
  with a delta-method standard error of the mean CS for MNL; supports
  `newdata` for policy ΔCS analysis

## choicer 0.1.0

CRAN release: 2026-05-20

Initial CRAN release.

### Supported models

- **Multinomial Logit**
  ([`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md))
  — estimation, prediction, elasticities, diversion ratios, BLP
  contraction
- **Mixed Logit**
  ([`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md))
  — normal and log-normal random coefficients, correlated random
  coefficients via Cholesky, Halton draws, elasticities, BLP contraction
- **Nested Logit**
  ([`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md))
  — estimation with nest-specific dissimilarity (lambda) parameters

### S3 class system

- Parent class `choicer_fit` with subclasses `choicer_mnl`,
  `choicer_mxl`, `choicer_nl`
- Standard methods: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html)
- Classed data objects: `choicer_data_mnl`, `choicer_data_mxl`,
  `choicer_data_nl` from `prepare_*_data()`

### Post-estimation generics

- [`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md)
  — methods for MNL and MXL
- [`diversion_ratios()`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.md)
  — method for MNL
- [`blp()`](https://fpcordeiro.github.io/choicer/reference/blp.md) — BLP
  contraction for MNL and MXL

### API

- Dual workflow for all `run_*logit()` functions: convenience (pass
  `data` + column names) or advanced (pass pre-prepared `input_data`)
- Pluggable optimizer via `optimizer = "nloptr" | "optim" | <function>`
- [`prepare_nl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md)
  added for nested logit data preparation

### Computation

- C++ likelihoods, gradients, and analytical Hessians via
  Rcpp/RcppArmadillo
- OpenMP parallelization over individuals
- Log-sum-exp trick throughout for numerical stability
