# Runs nested logit estimation

Estimates a nested logit model via maximum likelihood.

## Usage

``` r
run_nestlogit(
  data = NULL,
  id_col = NULL,
  alt_col = NULL,
  choice_col = NULL,
  covariate_cols = NULL,
  nest_col = NULL,
  input_data = NULL,
  use_asc = TRUE,
  theta_init = NULL,
  param_names = NULL,
  optimizer = NULL,
  control = list(),
  weights = NULL,
  weights_col = NULL,
  outside_opt_label = NULL,
  include_outside_option = FALSE,
  keep_data = TRUE,
  se_method = c("hessian", "numeric", "bhhh", "sandwich"),
  nloptr_opts = NULL
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

  Vector of column names for covariates.

- nest_col:

  Name of the column mapping each alternative to its nest (convenience
  workflow).

- input_data:

  List containing prepared input data for estimation (advanced
  workflow). Mutually exclusive with `data`.

- use_asc:

  Logical indicating whether to include alternative specific constants
  (ASCs).

- theta_init:

  Optional initial parameter vector. If `NULL`, a default vector is
  used.

- param_names:

  Optional vector of parameter names. If `NULL`, default names are
  generated.

- optimizer:

  Optimizer to use: `"nloptr"` (default), `"optim"`, or a custom
  function. See
  [`run_mnlogit`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md)
  for details.

- control:

  List of optimizer-specific control parameters.

- weights:

  Optional weight vector (convenience workflow). If `NULL`, equal
  weights are used. All weights must be finite and strictly positive.

- weights_col:

  Optional name of a column in `data` holding per-row weights
  (convenience workflow only). The column must be constant within each
  `id_col` (one weight per choice situation) and is collapsed
  accordingly. Mutually exclusive with `weights`. All weights must be
  finite and strictly positive. Used for choice-based / WESML weighting;
  pair with `se_method = "sandwich"` for valid inference.

- outside_opt_label:

  Label for the outside option (convenience workflow).

- include_outside_option:

  Logical whether to include an outside option (convenience workflow).

- keep_data:

  Logical. If `TRUE` (default), stores prepared data in the returned
  object for post-estimation functions.

- se_method:

  Method for computing standard errors: `"hessian"` (default, analytical
  Hessian via `nl_loglik_hessian_parallel`), `"numeric"`
  (finite-difference oracle via `nl_loglik_numeric_hessian`), `"bhhh"`
  (outer product of gradients via `nl_bhhh_parallel`), or `"sandwich"`
  (robust Huber–White / WESML variance \\A^{-1} B A^{-1}\\). Use
  `"sandwich"` under choice-based / WESML weighting.

- nloptr_opts:

  Deprecated. Use `optimizer` and `control` instead.

## Value

A `choicer_nl` object (inherits from `choicer_fit`). Standard S3 methods
available: [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html).

## Details

Two workflows are supported:

- Convenience:

  Supply `data` and column names (including `nest_col`). Data
  preparation
  ([`prepare_nl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md))
  is handled automatically.

- Advanced:

  Call
  [`prepare_nl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md)
  (or build the input list manually) and pass it via `input_data`.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 100; J <- 4
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1          x2
#>      <int> <int>      <num>       <num>
#>   1:     1     1  1.3709584  1.33491259
#>   2:     1     2 -0.5646982 -0.86927176
#>   3:     1     3  0.3631284  0.05548695
#>   4:     1     4  0.6328626  0.04906691
#>   5:     2     1  0.4042683 -0.57835573
#>  ---                                   
#> 396:    99     4  1.0965134 -1.07287540
#> 397:   100     1  0.4420131 -2.29297143
#> 398:   100     2  0.2410163 -1.20720685
#> 399:   100     3 -0.2556077  0.11410943
#> 400:   100     4  0.9310329 -1.03329708
dt[, nest := ifelse(alt <= 2, "A", "B")]
#>         id   alt         x1          x2   nest
#>      <int> <int>      <num>       <num> <char>
#>   1:     1     1  1.3709584  1.33491259      A
#>   2:     1     2 -0.5646982 -0.86927176      A
#>   3:     1     3  0.3631284  0.05548695      B
#>   4:     1     4  0.6328626  0.04906691      B
#>   5:     2     1  0.4042683 -0.57835573      A
#>  ---                                          
#> 396:    99     4  1.0965134 -1.07287540      B
#> 397:   100     1  0.4420131 -2.29297143      A
#> 398:   100     2  0.2410163 -1.20720685      A
#> 399:   100     3 -0.2556077  0.11410943      B
#> 400:   100     4  0.9310329 -1.03329708      B
dt[, choice := 0L]
#>         id   alt         x1          x2   nest choice
#>      <int> <int>      <num>       <num> <char>  <int>
#>   1:     1     1  1.3709584  1.33491259      A      0
#>   2:     1     2 -0.5646982 -0.86927176      A      0
#>   3:     1     3  0.3631284  0.05548695      B      0
#>   4:     1     4  0.6328626  0.04906691      B      0
#>   5:     2     1  0.4042683 -0.57835573      A      0
#>  ---                                                 
#> 396:    99     4  1.0965134 -1.07287540      B      0
#> 397:   100     1  0.4420131 -2.29297143      A      0
#> 398:   100     2  0.2410163 -1.20720685      A      0
#> 399:   100     3 -0.2556077  0.11410943      B      0
#> 400:   100     4  0.9310329 -1.03329708      B      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1          x2   nest choice
#>      <int> <int>      <num>       <num> <char>  <int>
#>   1:     1     1  1.3709584  1.33491259      A      1
#>   2:     1     2 -0.5646982 -0.86927176      A      0
#>   3:     1     3  0.3631284  0.05548695      B      0
#>   4:     1     4  0.6328626  0.04906691      B      0
#>   5:     2     1  0.4042683 -0.57835573      A      0
#>  ---                                                 
#> 396:    99     4  1.0965134 -1.07287540      B      0
#> 397:   100     1  0.4420131 -2.29297143      A      0
#> 398:   100     2  0.2410163 -1.20720685      A      0
#> 399:   100     3 -0.2556077  0.11410943      B      1
#> 400:   100     4  0.9310329 -1.03329708      B      0

fit <- run_nestlogit(
  data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = c("x1", "x2"), nest_col = "nest"
)
#> Optimization run time 0h:0m:0.01s
summary(fit)
#> Nested Logit (NL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1          -0.135699   0.188113  -0.7214  4.71e-01  
#> x2           0.295655   0.222668   1.3278  1.84e-01  
#> Lambda_1     6.296863  23.550214   0.2674  7.89e-01  
#> Lambda_2     1.640148   1.813470   0.9044  3.66e-01  
#> ASC_2       -2.778324  10.488956  -0.2649  7.91e-01  
#> ASC_3        2.033850  11.737696   0.1733  8.62e-01  
#> ASC_4        1.520384  11.782281   0.1290  8.97e-01  
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -134.492 
#> AIC: 282.983  | BIC: 301.219 
#> McFadden R2: 0.030 (adj: -0.021) | Hit rate: 0.380 
#> N: 100  | Parameters: 7 
#> Optimization time: 0.01 s
#> Convergence: 1 ( NLOPT_SUCCESS: Generic success return value. )
# }
```
