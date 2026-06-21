# Parameter recovery table

Compares fitted coefficients to a set of true parameter values on the
same scale as the estimator's internal parameterization. Returns one row
per estimated parameter with true value, estimate, standard error, bias,
relative bias (%), z-score against the truth, Wald CI, and a coverage
indicator.

## Usage

``` r
recovery_table(object, truth = NULL, level = 0.95, ...)

# S3 method for class 'choicer_fit'
recovery_table(object, truth = NULL, level = 0.95, ...)

# S3 method for class 'choicer_mnp'
recovery_table(object, truth = NULL, level = 0.95, ...)

# S3 method for class 'choicer_mc'
recovery_table(object, truth = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `choicer_fit` object (MNL, MXL, or NL) or a `choicer_mc` result.

- truth:

  Either a `choicer_sim` object (whose `$true_params` will be used) or a
  named list of true parameter values.

- level:

  Confidence level for the Wald CI and coverage indicator. Default
  `0.95`.

- ...:

  Unused.

## Value

See class-specific methods.

## Details

For MXL fits the `sigma` block compares the raw Cholesky parameters
(`L_params`), not the reconstructed covariance matrix. For log-normal
random-coefficient means the raw `mu` estimate is compared directly;
callers who want recovery on the DGP scale (`exp(mu)`) should transform
both sides before calling.

When the estimator has normalized the first inside alternative's ASC to
zero (which happens for MNL/MXL with `include_outside_option = FALSE`
and no outside option baked into the fit), the first entry of
`truth$delta` is dropped before the comparison so lengths match.

## Methods (by class)

- `recovery_table(choicer_fit)`: Returns a `choicer_recovery` object (a
  `data.table`) with columns `parameter`, `group`, `true`, `estimate`,
  `se`, `bias`, `rel_bias_pct`, `z_vs_true`, `lower_ci`, `upper_ci`,
  `covers`.

- `recovery_table(choicer_mnp)`: Method for Bayesian MNP fits
  (`choicer_mnp`). The `estimate` column holds posterior means of the
  identified draws and `se` holds their posterior standard deviations,
  so `lower_ci` / `upper_ci` are normal-approximation credible
  intervals. In addition to the `beta` and `asc` blocks, a `sigma` block
  compares the identified covariance of the utility differences (lower
  triangle in the estimator's row-major `Sigma_ij` order); its first row
  is the `sigma_11 = 1` normalization and is exact by construction.
  `truth` must be on the identified scale, as returned by
  [`simulate_mnp_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnp_data.md).

- `recovery_table(choicer_mc)`: For a `choicer_mc` object, delegates to
  `summary(object, level)` and returns a `choicer_mc_summary`. Inspect
  `object$replications` directly for per-rep detail.

## Examples

``` r
# \donttest{
sim <- simulate_mnl_data(N = 2000, J = 4, seed = 123)
fit <- run_mnlogit(
  data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = c("x1", "x2"),
  outside_opt_label = 0L, include_outside_option = FALSE, use_asc = TRUE
)
#> Optimization run time 0h:0m:0.14s
recovery_table(fit, sim)
#> <choicer_recovery> model=choicer_mnl level=0.95
#>    parameter  group  true estimate     se    bias rel_bias_pct z_vs_true
#>       <char> <char> <num>    <num>  <num>   <num>        <num>     <num>
#> 1:        x1   beta   0.8   0.7963 0.0561 -0.0037      -0.4612   -0.0657
#> 2:        x2   beta  -0.6  -0.5433 0.0558  0.0567      -9.4436    1.0160
#> 3:     ASC_1    asc   0.5   0.5680 0.0687  0.0680      13.5941    0.9895
#> 4:     ASC_2    asc  -0.5  -0.2788 0.0821  0.2212     -44.2448    2.6959
#> 5:     ASC_3    asc   0.5   0.5590 0.0686  0.0590      11.8065    0.8603
#> 6:     ASC_4    asc  -0.5  -0.4733 0.0859  0.0267      -5.3362    0.3106
#>    lower_ci upper_ci covers
#>       <num>    <num> <lgcl>
#> 1:   0.6863   0.9064   TRUE
#> 2:  -0.6526  -0.4340   TRUE
#> 3:   0.4333   0.7026   TRUE
#> 4:  -0.4396  -0.1179  FALSE
#> 5:   0.4245   0.6935   TRUE
#> 6:  -0.6417  -0.3049   TRUE
# }
```
