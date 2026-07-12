# choicer 0.2.0

## Public API cleanup (breaking)

- Low-level C++ likelihood, gradient, Hessian, score, prediction,
  post-estimation, and Gibbs-sampler wrappers are now internal implementation
  details rather than exported package functions. High-level fitting,
  preparation, S3 post-estimation, simulation, and recovery APIs are unchanged.
  The expert-facing raw BLP contractions, `get_halton_normals()`,
  `set_num_threads()`, and `thread_info()` remain public. With no downstream
  CRAN dependencies, v0.2.0 is the least disruptive point to narrow this
  surface before applications depend on unstable kernel signatures.

## Robust and clustered standard errors (MNL / MXL / NL)

- `vcov()` on a fitted MNL, MXL, or NL model gains `type =` and `cluster =`
  arguments for post-hoc variance recomputation without refitting (requires
  `keep_data = TRUE`): `"hessian"` (inverse analytical Hessian), `"bhhh"`
  (OPG), `"robust"` (Huber-White / WESML sandwich), and `"cluster"`
  (cluster-robust sandwich over within-cluster sums of weighted scores — use
  it when the same decision maker contributes several choice situations).
  With no arguments, `vcov()` returns the as-fitted variance, unchanged.
- New `cluster_col=` argument on `run_mnlogit()` / `run_mxlogit()` /
  `run_nestlogit()` (and the corresponding `prepare_*_data()` functions):
  supplies per-situation cluster labels at fit time and selects the new
  `se_method = "cluster"`.
- All score-based variances (BHHH, robust, cluster) are now assembled in R
  from one per-situation score matrix, computed by new internal C++ kernels
  that reuse the BHHH loop bodies. The existing `se_method = "sandwich"` /
  `wesml_vcov()` results are unchanged (the robust meat is the `w^2` special
  case of the shared path).
- Scope note (MXL): clustering repairs the inference, not the estimand. The
  MXL simulated likelihood treats each choice situation as an independent
  draw from the mixing distribution (cross-sectional MSL, not the panel
  product form), so clustered standard errors on panel data are robust to
  within-person dependence but do not turn the fit into a panel mixed logit.
  For panel random coefficients use `run_hmnlogit()` (`person_col=`).

## Corrections

- `blp()` for an MXL fit using `draws = "generate"` now carries the fitted
  draw count, seed, and digit-permutation setting through every contraction
  evaluation. Previously that path could evaluate counterfactual shares with
  fallback draw metadata rather than the simulation design used for the fit.
- On-the-fly MXL draws now accept all 128 prime-base dimensions implemented by
  the generator (the prior guard rejected the 128th dimension); 129 or more
  random coefficients still fail early with an explicit limit.
- HMNL prediction and welfare now use a taskwise max-shifted softmax/logsum, so
  extreme counterfactual utilities remain finite. Counterfactual compensating
  variation accepts explicit non-negative task weights and matches baseline
  and policy choice situations by (`person_col`, `id_col`) before subtraction;
  it now rejects added, dropped, or substituted tasks instead of pairing them
  silently by sorted position.
- The documentation previously claimed that MCMC draws from the Gibbs
  samplers (`run_mnprobit()`, `run_hmnlogit()`, `run_hmnprobit()`) are
  bitwise reproducible regardless of the OpenMP thread count. Independent
  validation disconfirmed this: draws are reproducible given the seed and a
  fixed thread count, and across different thread counts are invariant only
  up to floating-point reduction-order round-off (~1e-15), not bitwise. All
  man pages, vignettes, and NEWS entries now state the correct guarantee.
- Documentation now distinguishes the frequentist models' analytical
  derivatives from the Bayesian C++ samplers, labels `run_mxlogit()` as a
  cross-sectional simulated likelihood, and no longer directs users to
  WTP-space estimation or bounded/censored mixing distributions as if choicer
  implemented them. The README also labels `mode_choice` shares and welfare as
  sample-design quantities unless external population shares and WESML weights
  are supplied.
- The MNP tutorial no longer recommends unavailable user-specified starts and
  now states its current post-estimation boundary. HMNP documentation now makes
  explicit that it uses iid utility-level normal shocks rather than the full
  differenced-error covariance estimated by `run_mnprobit()`.

## Hierarchical Bayes (HMNL/HMNP) convergence diagnostics and multi-chain performance

- New `ess(draws)` (rank-normalized bulk and tail effective sample size) and
  `mcse(draws, kind = c("mean", "median"))` (Monte Carlo standard error of
  the mean or median), following Vehtari, Gelman, Simpson, Carpenter &
  Bürkner (2021, *Bayesian Analysis*). `rhat()` gains a `rank = TRUE` option
  for the same paper's rank-normalized, folded R-hat; the default
  `rank = FALSE` reproduces the original split R-hat exactly.
- New `traceplot()` generic, with a `choicer_hb` method, for visual
  inspection of the `b`/`theta`/`sigma_d2` chains (and user-selected `delta`
  columns).
- `summary()` on a `choicer_hmnl`/`choicer_hmnp` fit now prints one
  consolidated multi-chain diagnostics table (R-hat, ESS bulk, ESS tail, and
  MCSE per parameter block, plus a single worst-case summary row spanning
  all `J` `delta_j` alternative effects) in place of the previous bare
  split-R-hat printout. HMNP now always prints an explicit
  "conjugate — no acceptance step" line where HMNL prints its beta/delta
  acceptance rates. The fit-time convergence warning now checks rank-R-hat
  and ESS bulk across every tracked parameter, including all `J` `delta_j`
  columns (previously only `b`/`theta`/`sigma_d2` were checked).
- `run_hmnlogit()` / `run_hmnprobit()` fits now retain **all** requested
  chains' hierarchical draws in a new `object$chains` field (a list, one
  element per chain: `b`, `w_vech`, `delta`, `theta`, `sigma_d2`, plus
  `loglik` for HMNL / `sigma2` for HMNP), which the multi-chain diagnostics
  above consume. `object$draws` (chain 1) and `object$beta_i`
  (chain-1-only respondent taste summaries/draws) are unchanged.
- Chains requested via `chains=` are sampled sequentially, and all of their
  hierarchical draws are retained (`object$chains`) for the multi-chain
  diagnostics above. Parallel multi-chain execution is planned for a future
  release.
- The `keep_beta_i = "draws"` memory guard on both models is now based on a
  measured ~1.9x `Rcpp::wrap()`/R-list overhead factor (the previous formula
  under-estimated true memory use by roughly 2x) and accounts for all
  requested chains; it still fails fast before the chain runs rather than
  after an out-of-memory blow-up.

## Hierarchical Bayesian MNL and MNP (BLP-style random-effects ASCs)

- New `run_hmnlogit()` (class `choicer_hmnl`) and `run_hmnprobit()` (class
  `choicer_hmnp`): hierarchical Bayes discrete choice with two decoupled
  random-effect levels — respondent-level structural tastes
  `beta_i ~ N(b, W)` (`W` is `K x K`, no ASC columns) and a global BLP-style
  alternative effect `delta_j = z_j' theta + xi_j`, `xi_j ~ N(0, sigma_d^2)`
  (one scalar variance regardless of `J`; partial pooling toward the
  characteristics-based mean; posterior-predictive `delta` for alternatives
  outside the estimation sample). Both models carry a first-class implicit
  outside option (systematic utility 0 plus its own shock) that anchors the
  location of `delta` — no base alternative or sum-to-zero constraint.
  Panel and cross-sectional (`person_col = NULL`, `T_i = 1`) modes share
  one code path.
- HMNL: adaptive RW-Metropolis-within-Gibbs on the mnprobit engine pattern
  (single persistent OpenMP region, per-(iteration, unit) RNG streams;
  draws are reproducible given the seed and a fixed thread count, and
  invariant across thread counts only up to floating-point reduction-order
  round-off, ~1e-15), with a strictly serial `delta` sweep
  (the conditionals are coupled through the softmax denominators) at O(1)
  incremental cost per affected task; supports log-normal coordinates via
  `rc_dist`. HMNP: fully conjugate Albert-Chib augmentation in
  un-differenced utility space with iid `N(0, sigma^2)` shocks, a
  non-identified `sigma^2` chain (parameter expansion), and per-draw scale
  normalization of every identified quantity.
- Priors: `b ~ N(b_bar, A^-1)`, `W ~ IW(nu, V)`,
  `theta ~ N(theta_bar, A_theta^-1)`, and `sigma_d ~ half-Cauchy(0, s_d)`
  via the Makalic-Schmidt conjugate scale mixture (IG fallback available).
- New preps `prepare_hmnl_data()` / `prepare_hmnp_data()` (two-level
  person/task indexing, implicit-outside convention shared with
  `prepare_mnl_data()`, alternative-level design `Z`, optional
  control-function residual column for price endogeneity per Petrin &
  Train 2010) and DGPs `simulate_hmnl_data()` / `simulate_hmnp_data()`.
- Post-estimation on the shared `choicer_hb` class, all posterior-draw
  based: `predict()` (population/individual, entry counterfactuals via the
  posterior-predictive `delta`; HMNP probabilities by deterministic 20-node
  1-D Gauss-Hermite approximation), `wtp()` (posterior median + quantile
  intervals), `logsum()` / `consumer_surplus()` (HMNL-only; probit Emax is
  roadmapped), `elasticities()` / `diversion_ratios()` (common-random-path
  perturbation engine including the outside option), `recovery_table()`,
  plus `rhat()` (split R-hat) and `ppc_shares()` diagnostics.
- Two new math articles: `vignettes/articles/hierarchical_mnl_math.Rmd`
  and `hierarchical_mnp_math.Rmd`; recovery demos in
  `inst/simulations/hmnl_simulation.R` / `hmnp_simulation.R`.

## WESML sandwich inference for MNL and nested logit

- `run_mnlogit()` and `run_nestlogit()` now compute the robust (Huber–White /
  WESML) sandwich variance `V = A^{-1} B A^{-1}`, at full feature parity with
  `run_mxlogit()`. Both gain `se_method = "sandwich"` (bread = weighted negated
  Hessian, meat = weight-squared OPG) and `se_method = "bhhh"` (ordinary
  outer-product-of-gradients). `run_mnlogit()`'s `se_method` therefore now
  accepts `"hessian"` (default), `"bhhh"`, or `"sandwich"`; `run_nestlogit()`'s
  accepts `"hessian"` (default), `"numeric"`, `"bhhh"`, or `"sandwich"`.
- New C++ kernels `mnl_bhhh_parallel()` and `nl_bhhh_parallel()` accumulate the
  weighted outer product of per-individual scores. Their per-individual score is
  weight-free, so passing `weights = w` yields the BHHH information and
  `weights = w^2` yields the sandwich meat. The NL kernel includes the full
  beta/lambda/delta score blocks (singleton-nest lambdas fixed to 1 contribute
  no score).
- `run_mnlogit()`, `run_nestlogit()`, `prepare_mnl_data()` and
  `prepare_nl_data()` gain a `weights_col` argument: a row-level weight column
  is collapsed to one weight per choice situation (validated constant within
  `id`). A choice-based-sampling provenance guard auto-adopts the recorded
  WESML weight column and errors rather than silently fitting unweighted under a
  WESML label; non-uniform weights under a non-sandwich `se_method` emit a
  warning.
- Behavior change: weighted fits now emit a steering warning recommending
  `se_method = "sandwich"` when non-uniform weights are supplied with the default
  (`"hessian"`) or `"bhhh"` method; the `"bhhh"` case gets a sharper message
  explaining that BHHH/OPG is not a valid WESML correction (its meat is `w^1`,
  not `w^2`). Point estimates and standard errors are unchanged.
- `wesml_vcov()` now dispatches on `choicer_mnl` and `choicer_nl` (in addition
  to `choicer_mxl`), returning the post-hoc sandwich variance from a fit stored
  with `keep_data = TRUE`.
- `summary()` for MNL and NL fits now reports the standard-error method and any
  WESML weighting in the printed footer.
- Added "Choice-Based Sampling and WESML Weighting" sections to the multinomial
  logit and nested logit derivation vignettes.

## Weighting safety hardening

- Weights are now validated to be finite and strictly positive in
  `prepare_mnl_data()`, `prepare_nl_data()`, and `prepare_mxl_data()` (covering
  both the `weights=` and `weights_col=` paths). Zero, negative, or non-finite
  weights previously could slip through and silently invalidate weighted and
  WESML sandwich inference (weight `w` enters the bread, `w^2` the meat); they
  now error with an actionable message.
- Advanced-mode fits (`input_data=` passed directly to `run_mnlogit()`,
  `run_nestlogit()`, or `run_mxlogit()`) that carry WESML `choice_sampling`
  provenance but resolve to uniform weights now error instead of warning. The
  message explains how to proceed: bake the non-uniform WESML weights into
  `input_data` via `prepare_*_data(weights=/weights_col=)`, or strip the
  provenance with `attr(input_data, "choice_sampling") <- NULL` for a deliberate
  unweighted fit. Convenience-mode behavior is unchanged.

## Documentation

- Added eight vignettes: a getting-started tour ("Discrete choice from data to
  policy, in a dozen lines"), one per model (multinomial logit, mixed logit,
  nested logit, Bayesian multinomial probit, hierarchical Bayes), one on
  choice-based sampling and WESML weights, and one on standard errors
  (Hessian / BHHH / robust / cluster-robust, and when to use which).
- Added the `mode_choice` data set: the classic Greene & Hensher intercity
  travel-mode choice data (210 travellers x 4 modes), in choicer's long layout,
  used by the getting-started vignette. `?mode_choice` and the vignettes now
  document that the sample is choice-based (car under-sampled): slopes and WTP
  ratios are unaffected in the ASC-saturated logit, while constants, shares,
  and surplus levels inherit the design — see the WESML vignette for the
  correction.
- Added a pkgdown website configuration, including the model derivation notes as
  "The math behind choicer" articles.
- Expanded the identification discussions across the documentation: price
  endogeneity and the control-function route in the getting-started
  model-choice section (Petrin & Train 2010; `blp()` as the Berry 1994
  inversion), what identifies the nested-logit dissimilarity parameters, the
  Keane (1992) covariance-identification caveat for the multinomial probit,
  and the classical MSL asymptotics (fixed-`S` bias, `S` growth conditions) in
  the mixed logit math note.

## Mixed logit — on-the-fly digit-permuted Halton draws

- `run_mxlogit()` gains three new arguments: `draws`, `seed`, and `scramble`.
  - `draws = "store"` (default) keeps the existing behavior: a full
    K_w × S × N Halton cube is pre-materialized and stored in memory.
  - `draws = "generate"` activates the new mode: each individual's S draws are
    computed on the fly in C++ from a compact seed, eliminating the O(N)
    `eta_draws` cube. This is the recommended choice when N is large or memory
    is constrained.
  - `scramble = "permuted"` (default when `draws = "generate"`) applies one
    seeded permutation per dimension and base-digit position, shared across
    sequence indices. This is a deterministic position-wise digit permutation,
    not Owen's nested-uniform scramble, and no standard randomized-QMC
    unbiasedness or replicate-error guarantee is claimed. The historical value
    `"owen"` remains as a deprecated alias for compatibility.
    `scramble = "none"` reproduces the randtoolbox sequence exactly.
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

- `run_mnprobit()` — Bayesian multinomial probit via Gibbs sampling with data augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). Runs the non-identified chain with conjugate priors and reports identified quantities normalized per draw by `sigma_11`. New `choicer_mnp` posterior object with `summary()` (posterior mean, SD, credible intervals), `coef()`, `vcov()`, `nobs()`; math note in `vignettes/articles/bayesian_multinomial_probit_math.Rmd`
- C++ MCMC infrastructure built from scratch (`src/rng.h`, `src/bayes_samplers.h`): xoshiro256++/splitmix64 RNG with one stream per (iteration, observation) — draws are reproducible given the seed and a fixed OpenMP thread count (invariant across thread counts only up to floating-point reduction-order round-off, ~1e-15) — plus exact truncated-normal, multivariate-normal, Wishart (Bartlett), and inverse-Wishart samplers. The truncated normal picks the cheapest exact method per region (naive normal rejection in high-mass regions, Robert 1995 exponential rejection in the tail, inverse CDF for narrow intervals)
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
