# Monte Carlo parameter recovery

Replicates a `(DGP -> fit)` cycle `R` times with independent seeds and
collects per-parameter estimates, standard errors, bias, and coverage.
Returns a `choicer_mc` object; call
[`summary()`](https://rdrr.io/r/base/summary.html) for aggregated
statistics (mean estimate, bias, RMSE, coverage rate, convergence rate).

## Usage

``` r
monte_carlo(
  sim_fun,
  fit_fun,
  R = 100,
  seed = 1L,
  parallel = FALSE,
  progress = TRUE,
  ...
)
```

## Arguments

- sim_fun:

  Function of `seed` returning a `choicer_sim`.

- fit_fun:

  Function of a `choicer_sim` returning a `choicer_fit`.

- R:

  Number of replications.

- seed:

  Base integer seed. Replication `r` uses `seed + r - 1L`.

- parallel:

  Logical; if `TRUE` and `future.apply` is available, run replications
  in parallel using the user's active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html).

- progress:

  Logical; print a one-line progress update per iteration in serial
  mode. Ignored when `parallel = TRUE`.

- ...:

  Unused.

## Value

A `choicer_mc` object: a list with elements `replications` (a long
`data.table` with one row per estimated parameter per replication) and
`meta` (run metadata).

## Details

Each iteration calls `sim_fun(seed = seed + r - 1L)`, then
`fit_fun(sim)`. Write `sim_fun` as a closure that captures `N`, `J`, and
other DGP settings and forwards `seed`. Write `fit_fun` as a closure
that takes a `choicer_sim` and returns a fitted `choicer_fit` object,
wrapping any data-preparation, draws, or optimizer-control setup.

## Examples

``` r
# \donttest{
sim_fun <- function(seed) simulate_mnl_data(N = 1000, J = 4, seed = seed)
fit_fun <- function(sim) run_mnlogit(
  data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = c("x1", "x2"), outside_opt_label = 0L,
  include_outside_option = FALSE, use_asc = TRUE,
  control = list(print_level = 0L)
)
mc <- monte_carlo(sim_fun, fit_fun, R = 5, seed = 1L, progress = FALSE)
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
summary(mc)
#> <choicer_mc_summary> R=5 convergence_rate=1 coverage_level=0.95
#>    parameter  group  true R_success mean_est median_est sd_est mean_se    bias
#>       <char> <char> <num>     <int>    <num>      <num>  <num>   <num>   <num>
#> 1:        x1   beta   0.8         5   0.8162     0.8058 0.1233  0.0818  0.0162
#> 2:        x2   beta  -0.6         5  -0.6625    -0.6725 0.0574  0.0801 -0.0625
#> 3:     ASC_1    asc   0.5         5   0.5444     0.5149 0.0996  0.0975  0.0444
#> 4:     ASC_2    asc  -0.5         5  -0.3858    -0.3487 0.1162  0.1177  0.1142
#> 5:     ASC_3    asc   0.5         5   0.5362     0.5979 0.1316  0.0978  0.0362
#> 6:     ASC_4    asc  -0.5         5  -0.4600    -0.5735 0.2051  0.1208  0.0400
#>      rmse coverage
#>     <num>    <num>
#> 1: 0.1115      0.8
#> 2: 0.0809      1.0
#> 3: 0.0995      1.0
#> 4: 0.1544      1.0
#> 5: 0.1232      1.0
#> 6: 0.1878      0.8
# }
```
