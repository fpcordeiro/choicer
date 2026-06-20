# Runs mixed logit estimation

Estimates a mixed logit model via simulated maximum likelihood.

## Usage

``` r
run_mxlogit(
  data = NULL,
  id_col = NULL,
  alt_col = NULL,
  choice_col = NULL,
  covariate_cols = NULL,
  random_var_cols = NULL,
  input_data = NULL,
  eta_draws = NULL,
  S = 100L,
  rc_dist = NULL,
  rc_mean = FALSE,
  rc_correlation = FALSE,
  use_asc = TRUE,
  theta_init = NULL,
  lower = NULL,
  upper = NULL,
  optimizer = NULL,
  control = list(),
  se_method = c("hessian", "bhhh", "sandwich"),
  scale_vars = c("none", "sd", "mad", "iqr"),
  weights = NULL,
  outside_opt_label = NULL,
  include_outside_option = FALSE,
  draws = c("store", "generate"),
  seed = NULL,
  scramble = c("owen", "none"),
  keep_data = TRUE,
  nloptr_opts = NULL,
  weights_col = NULL
)
```

## Arguments

- data:

  Data frame containing choice data (convenience workflow). Mutually
  exclusive with `input_data`.

- id_col:

  Name of the column identifying choice situations.

- alt_col:

  Name of the column identifying alternatives.

- choice_col:

  Name of the column indicating chosen alternative (1/0).

- covariate_cols:

  Vector of column names for fixed covariates.

- random_var_cols:

  Vector of column names for random coefficients.

- input_data:

  List output from
  [`prepare_mxl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mxl_data.md)
  (advanced workflow). Mutually exclusive with `data`.

- eta_draws:

  Array of shape K_w x S x N with standard normal draws. Required for
  the advanced workflow; auto-generated from `S` in the convenience
  workflow.

- S:

  Integer number of Halton draws per individual (convenience workflow
  only). Default 100.

- rc_dist:

  Integer vector indicating distribution of random coefficients (0 =
  normal, 1 = log-normal). Default: all normal.

- rc_mean:

  Logical indicating whether to estimate means for random coefficients.

- rc_correlation:

  Logical indicating whether random coefficients are correlated
  (convenience workflow). Ignored when `input_data` is used (taken from
  the prepared data).

- use_asc:

  Logical indicating whether to include alternative-specific constants.

- theta_init:

  Initial parameter vector in natural-scale units. If `NULL`, defaults
  to zeros for the \\\beta\\, \\\mu\\, and ASC blocks, and `log(0.5)` on
  the Cholesky diagonal (so each diagonal factor \\L\_{pp} = 0.5\\, i.e.
  a moderate random-coefficient variance of `0.25`). The
  zero-on-diagonal alternative corresponds to \\L\_{pp} = 1\\ (unit RC
  variance), which often lets the first L-BFGS step overshoot.

- lower, upper:

  Optional parameter bounds for the optimizer, in natural-scale units
  (forward-transformed internally to scaled space when
  `scale_vars != "none"`). Each accepts three forms:

  `NULL`

  :   (default) Unbounded (`-Inf`/`Inf`).

  Unnamed numeric vector of length `n_params`

  :   Full-length vector ordered exactly like `theta_init` (the
      nloptr-native form).

  Named numeric vector

  :   Names must be a subset of the parameter names (\\\beta\\ block:
      column names of `X`; \\\mu\\ block: `Mu_<col>` (if
      `rc_mean = TRUE`); Cholesky block: `L_<i><j>` for \\i \ge j\\; ASC
      block: `ASC_<level>`). Unlisted parameters default to
      \\\pm\infty\\. This is the recommended form for typical use, e.g.
      `lower = c(L_11 = -5, L_22 = -5)` to clip Cholesky diagonals.

- optimizer:

  Optimizer to use: `"nloptr"` (default), `"optim"`, or a custom
  function. See
  [`run_mnlogit`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md)
  for details.

- control:

  List of optimizer-specific control parameters.

- se_method:

  Method for computing standard errors. One of `"hessian"` (default) for
  the analytical Hessian of the simulated log-likelihood, `"bhhh"` for
  the BHHH/outer-product-of-gradients (OPG) estimator, or `"sandwich"`
  for the robust (Huber-White) variance \\V = A^{-1} B A^{-1}\\ (bread
  \\A\\ = weighted negated Hessian, meat \\B\\ = weight-squared OPG).
  Use `"sandwich"` for valid inference under choice-based / WESML
  weighting, where the inverse-Hessian and ordinary BHHH are invalid; it
  reduces to the usual robust variance under uniform weights. BHHH
  scales better to large problems (many alternatives or simulation
  draws) but may underestimate standard errors in finite samples or away
  from the optimum.

- scale_vars:

  Pre-estimation column scaling for design matrices. One of `"none"`
  (default), `"sd"` (sample standard deviation), `"mad"`
  ([`stats::mad`](https://rdrr.io/r/stats/mad.html), i.e. 1.4826
  \\\times\\ median absolute deviation; SD-equivalent under normality),
  or `"iqr"` (`stats::IQR(x) / 1.349`; also SD-equivalent under
  normality). When not `"none"`, every column of `X` and `W` is divided
  by the chosen scale before optimization to improve Hessian
  conditioning. Robust scales (`"mad"`/`"iqr"`) better capture the bulk
  for heavy-tailed columns where SD is dominated by outliers, but
  [`stats::mad`](https://rdrr.io/r/stats/mad.html) can return zero when
  more than half of a column's entries are identical (e.g., a sparse 0/1
  dummy) and will then trigger the same near-constant-column error as
  `"sd"`. Coefficients and standard errors are back-transformed to the
  user's natural units via the delta method, so reported quantities are
  invariant to this choice. Columns of `W` associated with log-normal
  random coefficients (`rc_dist == 1`) are passed through unchanged,
  since the shifted log-normal parameterization does not admit a
  closed-form back-transform under multiplicative scaling.

- weights:

  Optional weight vector (convenience workflow). If `NULL`, equal
  weights are used.

- outside_opt_label:

  Label for the outside option (convenience workflow).

- include_outside_option:

  Logical whether to include an outside option (convenience workflow).

- draws:

  Draw storage mode. One of `"store"` (default) or `"generate"`.
  `"store"` pre-materializes the full \\K_w \times S \times N\\ Halton
  cube (existing behavior, exact reproducibility). `"generate"` computes
  each individual's draws on-the-fly in C++ from a stored seed,
  eliminating the O(N) cube; recommended for memory-constrained or
  large-N settings. Uses randomized digit-scrambled (Owen 2017) Halton
  sequences when `scramble = "owen"`. Only supported in the convenience
  workflow.

- seed:

  Integer master seed for the on-the-fly generator. Used only when
  `draws = "generate"`. If `NULL` (default), a seed is drawn from R's
  RNG at call time (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) governs
  reproducibility). Ignored when `draws = "store"`.

- scramble:

  Scrambling mode for on-the-fly Halton draws. One of `"owen"` (default)
  for Owen (2017) digit scrambling, or `"none"` for plain unscrambled
  Halton (identity permutations). `"none"` reproduces the randtoolbox
  sequence exactly and is intended for testing; `"owen"` is strongly
  recommended for estimation with \\K_w \> 5\\. Used only when
  `draws = "generate"`.

- keep_data:

  Logical. If `TRUE` (default), stores prepared data in the returned
  object for post-estimation functions.

- nloptr_opts:

  Deprecated. Use `optimizer` and `control` instead.

- weights_col:

  Optional name of a column in `data` holding a per-row weight (constant
  within each choice situation). Mutually exclusive with `weights`; the
  recommended way to pass WESML weights from
  [`sample_by_choice`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
  /
  [`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md),
  since alignment is by id rather than by position. Convenience workflow
  only. If `data` carries choice-based-sampling provenance (a
  `"choice_sampling"` attribute, as attached by
  [`sample_by_choice`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
  /
  [`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md))
  and neither `weights` nor `weights_col` is supplied, the recorded
  weight column is auto-detected and applied (with a message); if that
  column is absent the call errors rather than silently fitting an
  unweighted model under a WESML label.

## Value

A `choicer_mxl` object (inherits from `choicer_fit`). Standard S3
methods available: [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html).

## Details

Two workflows are supported:

- Convenience:

  Supply `data` and column names. Data preparation
  ([`prepare_mxl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mxl_data.md))
  and Halton draw generation
  ([`get_halton_normals`](https://fpcordeiro.github.io/choicer/reference/get_halton_normals.md))
  are handled automatically.

- Advanced:

  Call
  [`prepare_mxl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mxl_data.md)
  and
  [`get_halton_normals`](https://fpcordeiro.github.io/choicer/reference/get_halton_normals.md)
  yourself, then pass the results via `input_data` and `eta_draws`.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 100; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#>         id   alt          x1           w1         w2
#>      <int> <int>       <num>        <num>      <num>
#>   1:     1     1  1.37095845 -0.004620768 -0.2484829
#>   2:     1     2 -0.56469817  0.760242168  0.4223204
#>   3:     1     3  0.36312841  0.038990913  0.9876533
#>   4:     2     1  0.63286260  0.735072142  0.8355682
#>   5:     2     2  0.40426832 -0.146472627 -0.6605219
#>  ---                                                
#> 296:    99     2 -0.47733551  0.160327395  0.1704735
#> 297:    99     3 -0.16626149 -0.433641942  1.2006682
#> 298:   100     1  0.86256338  1.537412419 -0.1634059
#> 299:   100     2  0.09734049 -2.170246577  1.2824759
#> 300:   100     3 -1.62561674  1.027004619  2.7271964
dt[, choice := 0L]
#>         id   alt          x1           w1         w2 choice
#>      <int> <int>       <num>        <num>      <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768 -0.2484829      0
#>   2:     1     2 -0.56469817  0.760242168  0.4223204      0
#>   3:     1     3  0.36312841  0.038990913  0.9876533      0
#>   4:     2     1  0.63286260  0.735072142  0.8355682      0
#>   5:     2     2  0.40426832 -0.146472627 -0.6605219      0
#>  ---                                                       
#> 296:    99     2 -0.47733551  0.160327395  0.1704735      0
#> 297:    99     3 -0.16626149 -0.433641942  1.2006682      0
#> 298:   100     1  0.86256338  1.537412419 -0.1634059      0
#> 299:   100     2  0.09734049 -2.170246577  1.2824759      0
#> 300:   100     3 -1.62561674  1.027004619  2.7271964      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt          x1           w1         w2 choice
#>      <int> <int>       <num>        <num>      <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768 -0.2484829      0
#>   2:     1     2 -0.56469817  0.760242168  0.4223204      1
#>   3:     1     3  0.36312841  0.038990913  0.9876533      0
#>   4:     2     1  0.63286260  0.735072142  0.8355682      0
#>   5:     2     2  0.40426832 -0.146472627 -0.6605219      0
#>  ---                                                       
#> 296:    99     2 -0.47733551  0.160327395  0.1704735      1
#> 297:    99     3 -0.16626149 -0.433641942  1.2006682      0
#> 298:   100     1  0.86256338  1.537412419 -0.1634059      0
#> 299:   100     2  0.09734049 -2.170246577  1.2824759      0
#> 300:   100     3 -1.62561674  1.027004619  2.7271964      1

fit <- run_mxlogit(
  data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 50L
)
#> Optimization run time 0h:0m:0.02s
summary(fit)
#> Mixed Logit (MXL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           0.005704   0.121525   0.0469  9.63e-01  
#> Sigma_11     0.000000   0.000000   0.0000  1.00e+00  
#> Sigma_22     0.000000   0.000000   0.0000  1.00e+00  
#> ASC_2        0.031083   0.248189   0.1252  9.00e-01  
#> ASC_3        0.089237   0.244715   0.3647  7.15e-01  
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Random coefficient covariance (Sigma):
#>              w1           w2
#> w1 3.520212e-18 0.000000e+00
#> w2 0.000000e+00 1.385186e-15
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -109.79 
#> AIC: 229.581  | BIC: 242.607 
#> McFadden R2: 0.001 (adj: -0.045) | Hit rate: 0.350 
#> N: 100  | Parameters: 5 
#> Optimization time: 0.02 s
#> Convergence: 1 ( NLOPT_SUCCESS: Generic success return value. )
# }
```
