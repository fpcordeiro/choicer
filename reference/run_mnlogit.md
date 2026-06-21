# Runs multinomial logit estimation

Estimates a multinomial logit model via maximum likelihood.

## Usage

``` r
run_mnlogit(
  data = NULL,
  id_col = NULL,
  alt_col = NULL,
  choice_col = NULL,
  covariate_cols = NULL,
  input_data = NULL,
  optimizer = NULL,
  control = list(),
  weights = NULL,
  weights_col = NULL,
  outside_opt_label = NULL,
  include_outside_option = FALSE,
  use_asc = TRUE,
  keep_data = TRUE,
  scale_vars = c("none", "sd", "mad", "iqr"),
  se_method = c("hessian", "bhhh", "sandwich"),
  nloptr_opts = NULL
)
```

## Arguments

- data:

  Data frame containing choice data (convenience workflow). Mutually
  exclusive with `input_data`.

- id_col:

  Name of the column identifying choice situations (individuals).

- alt_col:

  Name of the column identifying alternatives.

- choice_col:

  Name of the column indicating chosen alternative (1 = chosen, 0 = not
  chosen).

- covariate_cols:

  Vector of names of columns to be used as covariates.

- input_data:

  List output from
  [`prepare_mnl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
  (advanced workflow). Mutually exclusive with `data`.

- optimizer:

  Optimizer to use: `"nloptr"` (default), `"optim"`, or a custom
  function with signature `f(theta_init, eval_f, lower, upper, control)`
  where `eval_f(theta)` returns `list(objective, gradient)`. Must return
  a list with `par`/`value` (or `solution`/`objective`). If the custom
  function accepts `control` or `...`, the `control` argument is
  forwarded; otherwise it is silently ignored.

- control:

  List of optimizer-specific control parameters passed to the chosen
  optimizer (e.g., `list(maxeval = 2000)` for nloptr).

- weights:

  Optional vector of weights for each choice situation. If `NULL`, equal
  weights are used. All weights must be finite and strictly positive.

- weights_col:

  Optional name of a column in `data` holding per-row weights
  (convenience workflow only). The column must be constant within each
  `id_col` (one weight per choice situation) and is collapsed
  accordingly. Mutually exclusive with `weights`. All weights must be
  finite and strictly positive. Used for choice-based / WESML weighting;
  pair with `se_method = "sandwich"` for valid inference.

- outside_opt_label:

  Label for the outside option (if any). If `NULL`, no outside option is
  assumed.

- include_outside_option:

  Logical indicating whether to include an outside option in the model.

- use_asc:

  Logical indicating whether to include alternative-specific constants
  (ASCs) in the model.

- keep_data:

  Logical. If `TRUE` (default), stores prepared data in the returned
  object for [`predict()`](https://rdrr.io/r/stats/predict.html) and
  post-estimation functions.

- scale_vars:

  Pre-estimation column scaling for the design matrix. One of `"none"`
  (default), `"sd"` (sample standard deviation), `"mad"`
  ([`stats::mad`](https://rdrr.io/r/stats/mad.html)), or `"iqr"`
  (`stats::IQR(x) / 1.349`). When not `"none"`, every column of `X` is
  divided by the chosen scale before optimization to improve Hessian
  conditioning. Coefficients and standard errors are back-transformed to
  the user's natural units via the delta method, so reported quantities
  are invariant to this choice.

- se_method:

  Method for computing standard errors: `"hessian"` (default, analytical
  Hessian), `"bhhh"` (outer product of gradients), or `"sandwich"`
  (robust Huber–White / WESML variance \\A^{-1} B A^{-1}\\). Use
  `"sandwich"` under choice-based / WESML weighting.

- nloptr_opts:

  Deprecated. Use `optimizer` and `control` instead.

## Value

A `choicer_mnl` object (inherits from `choicer_fit`). Standard S3
methods available: [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`predict()`](https://rdrr.io/r/stats/predict.html).

## Details

Two workflows are supported:

- Convenience (default):

  Supply `data` and column names. Data preparation
  ([`prepare_mnl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md))
  is handled automatically.

- Advanced:

  Call
  [`prepare_mnl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
  yourself and pass the result via `input_data`.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 100; J <- 3; beta_true <- c(1.0, -0.5)
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt          x1           x2
#>      <int> <int>       <num>        <num>
#>   1:     1     1  1.37095845 -0.004620768
#>   2:     1     2 -0.56469817  0.760242168
#>   3:     1     3  0.36312841  0.038990913
#>   4:     2     1  0.63286260  0.735072142
#>   5:     2     2  0.40426832 -0.146472627
#>  ---                                     
#> 296:    99     2 -0.47733551  0.160327395
#> 297:    99     3 -0.16626149 -0.433641942
#> 298:   100     1  0.86256338  1.537412419
#> 299:   100     2  0.09734049 -2.170246577
#> 300:   100     3 -1.62561674  1.027004619
dt[, V := drop(as.matrix(.SD) %*% beta_true), .SDcols = c("x1","x2")]
#>         id   alt          x1           x2           V
#>      <int> <int>       <num>        <num>       <num>
#>   1:     1     1  1.37095845 -0.004620768  1.37326883
#>   2:     1     2 -0.56469817  0.760242168 -0.94481926
#>   3:     1     3  0.36312841  0.038990913  0.34363295
#>   4:     2     1  0.63286260  0.735072142  0.26532653
#>   5:     2     2  0.40426832 -0.146472627  0.47750464
#>  ---                                                 
#> 296:    99     2 -0.47733551  0.160327395 -0.55749920
#> 297:    99     3 -0.16626149 -0.433641942  0.05055948
#> 298:   100     1  0.86256338  1.537412419  0.09385717
#> 299:   100     2  0.09734049 -2.170246577  1.18246377
#> 300:   100     3 -1.62561674  1.027004619 -2.13911905
dt[, prob := exp(V) / sum(exp(V)), by = id]
#>         id   alt          x1           x2           V       prob
#>      <int> <int>       <num>        <num>       <num>      <num>
#>   1:     1     1  1.37095845 -0.004620768  1.37326883 0.68700257
#>   2:     1     2 -0.56469817  0.760242168 -0.94481926 0.06764341
#>   3:     1     3  0.36312841  0.038990913  0.34363295 0.24535402
#>   4:     2     1  0.63286260  0.735072142  0.26532653 0.33940231
#>   5:     2     2  0.40426832 -0.146472627  0.47750464 0.41962617
#>  ---                                                            
#> 296:    99     2 -0.47733551  0.160327395 -0.55749920 0.21814327
#> 297:    99     3 -0.16626149 -0.433641942  0.05055948 0.40069908
#> 298:   100     1  0.86256338  1.537412419  0.09385717 0.24525785
#> 299:   100     2  0.09734049 -2.170246577  1.18246377 0.72844833
#> 300:   100     3 -1.62561674  1.027004619 -2.13911905 0.02629382
dt[, choice := as.integer(alt == sample(alt, 1, prob = prob)), by = id]
#>         id   alt          x1           x2           V       prob choice
#>      <int> <int>       <num>        <num>       <num>      <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768  1.37326883 0.68700257      1
#>   2:     1     2 -0.56469817  0.760242168 -0.94481926 0.06764341      0
#>   3:     1     3  0.36312841  0.038990913  0.34363295 0.24535402      0
#>   4:     2     1  0.63286260  0.735072142  0.26532653 0.33940231      1
#>   5:     2     2  0.40426832 -0.146472627  0.47750464 0.41962617      0
#>  ---                                                                   
#> 296:    99     2 -0.47733551  0.160327395 -0.55749920 0.21814327      0
#> 297:    99     3 -0.16626149 -0.433641942  0.05055948 0.40069908      0
#> 298:   100     1  0.86256338  1.537412419  0.09385717 0.24525785      0
#> 299:   100     2  0.09734049 -2.170246577  1.18246377 0.72844833      1
#> 300:   100     3 -1.62561674  1.027004619 -2.13911905 0.02629382      0

fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s
summary(fit)
#> Multinomial Logit (MNL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           1.091325   0.192975   5.6553  1.56e-08  ***
#> x2          -0.511225   0.152893  -3.3437  8.27e-04  ***
#> ASC_2       -0.061376   0.288667  -0.2126  8.32e-01  
#> ASC_3       -0.074843   0.289651  -0.2584  7.96e-01  
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -81.2408 
#> AIC: 170.482  | BIC: 180.902 
#> McFadden R2: 0.261 (adj: 0.224) | Hit rate: 0.590 
#> N: 100  | Parameters: 4 
#> Optimization time: 0 s
#> Convergence: 1 ( NLOPT_SUCCESS: Generic success return value. )
coef(fit)
#>          x1          x2       ASC_2       ASC_3 
#>  1.09132487 -0.51122457 -0.06137617 -0.07484335 
AIC(fit)
#> [1] 170.4816
predict(fit, type = "shares")
#>      [,1]
#> [1,] 0.34
#> [2,] 0.32
#> [3,] 0.34
# }
```
