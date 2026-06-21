# Asymptotic diagnostics for a Monte Carlo study

Consumes a `choicer_mc` object and returns per-parameter asymptotic
diagnostics: Monte Carlo bias (with MC standard error), empirical SD of
the estimates, mean of the reported standard errors, SE-to-SD ratio
(information-matrix-equality check), Wald coverage at nominal 90 / 95 /
99 percent with Wilson confidence bands, moments of the studentized
statistic `z = (theta_hat - theta_0) / se`, and four normality tests on
`z` (Shapiro-Wilk, Anderson-Darling via
[`goftest::ad.test`](https://rdrr.io/pkg/goftest/man/ad.test.html), a
hand-coded Jarque-Bera statistic, and a one-sample Kolmogorov-Smirnov
test against `N(0, 1)`).

## Usage

``` r
mc_asymptotics(
  mc,
  level = 0.95,
  se_col = "se",
  conv_threshold = 0.99,
  se_ratio_threshold_floor = 0.1
)
```

## Arguments

- mc:

  A `choicer_mc` object returned by
  [`monte_carlo()`](https://fpcordeiro.github.io/choicer/reference/monte_carlo.md).

- level:

  Confidence level for the Wilson bands on coverage rates. Defaults to
  `0.95`.

- se_col:

  Name of the column in `mc$replications` to use as the standard-error
  source. Defaults to `"se"` (the Hessian-based SE stored by
  [`monte_carlo()`](https://fpcordeiro.github.io/choicer/reference/monte_carlo.md)).
  Callers that augment replications with an alternative SE flavor (e.g.,
  `"se_bhhh"` for a BHHH/OPG comparison) can pass that column name to
  recompute every SE-dependent diagnostic (`mean_se`, `se_ratio`,
  `mean_se_w`, `cov90/95/99`, z-moments, normality tests, pass flags)
  against that flavor. Useful for the information-matrix-equality check
  in Claim 4 of the MXL validation suite.

- conv_threshold:

  Numeric in `[0, 1]`. Minimum fraction of replications that must
  converge for the per-parameter `pass_convergence` flag to be `TRUE`.
  The flag compares `R_used / R_total` (per parameter) against this
  threshold. Defaults to `0.99`.

- se_ratio_threshold_floor:

  Numeric scalar. Minimum half-width for the `pass_se_ratio` band. The
  actual band used is
  `max(se_ratio_threshold_floor, 3 * 1.4 / sqrt(R_used))`, where the
  `1.4 / sqrt(R)` term approximates the large-sample SD of
  `mean_se / sd_emp`. The floor guarantees the band is never tighter
  than the historical hard cutoff. Defaults to `0.10`.

## Value

An object of class `choicer_mc_asymptotics` — a `data.table` with one
row per unique parameter and columns documented above — with `meta`
attached as an attribute (`attr(x, "meta")`).

## Details

Six logical pass / fail flags are attached to every parameter row:
`pass_bias` requires `|bias_mc_se| < 3`; `pass_se_ratio` requires
`|se_ratio - 1|` to lie within
`max(se_ratio_threshold_floor, 3 * 1.4 / sqrt(R_used))` (a noise-aware
band that widens at small `R_used` and tightens to the floor at large
`R_used`); `pass_cov95` requires the nominal 95 percent level to lie in
the Wilson band for empirical coverage; `pass_skew` requires
`|skew_z| < 0.3`; `pass_kurt` requires excess kurtosis of `z` in
`[-0.5, 1.0]`; `pass_convergence` requires the per-parameter convergence
rate (`R_used / R_total`) to meet `conv_threshold`.

Non-converged replications are excluded per parameter (reported in
`R_excluded`). Winsorized (5 percent / 95 percent) versions of `bias`,
`sd_emp`, and `mean_se` are reported in parallel columns (`bias_w`,
`sd_emp_w`, `mean_se_w`) so silent outlier exclusion is transparent to
the reader. Two robust SE-to-SD ratios accompany the Hessian-mean-based
`se_ratio`: `se_ratio_med` (median SE divided by the empirical SD) and
`se_ratio_w` (winsorized mean SE divided by the winsorized empirical
SD); both stay near `1` when 1-2 replications produce near-singular
Hessians that inflate `mean_se`. The companion `se_med` column reports
the median per-replication SE used by `se_ratio_med`. Neither robust
ratio drives a `pass_*` flag — they are purely informational.

Winsorized z-moment counterparts (`mean_z_w`, `sd_z_w`, `skew_z_w`,
`kurt_excess_z_w`) are reported alongside the raw z-moments and feed an
additional `pass_z_w` flag (Winsorized skew within the same band as
`pass_skew` AND Winsorized excess kurtosis within the same band as
`pass_kurt`). A companion `pass_cov95_w` flag is `TRUE` when either
`pass_cov95` is `TRUE` OR the per-rep Winsorized z-CI (the empirical 2.5
/ 97.5 percentiles of the Winsorized z) covers truth-zero. These two
flags are designed for boundary scenarios (e.g., near-zero variance
components) where a small number of reps with vanishing SE inflate the
raw z-moments without indicating an estimator defect.

## Examples

``` r
# \donttest{
sim_fun <- function(seed) simulate_mnl_data(N = 1000, J = 3, seed = seed)
fit_fun <- function(sim) run_mnlogit(
  data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = c("x1", "x2"), outside_opt_label = 0L,
  include_outside_option = FALSE, use_asc = TRUE,
  control = list(print_level = 0L)
)
mc <- monte_carlo(sim_fun, fit_fun, R = 50L, seed = 1L, progress = FALSE)
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0.01s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
#> Optimization run time 0h:0m:0s
mc_asymptotics(mc)
#> <choicer_mc_asymptotics> R_total=50 wilson_level=0.95
#>   pass_bias     : 5 / 5
#>   pass_se_ratio : 5 / 5
#>   pass_cov95    : 5 / 5
#>   pass_skew     : 2 / 5
#>   pass_kurt     : 4 / 5
#>   pass_convergence: 5 / 5
#>   all_pass      : 1 / 5
#>    parameter  group  true R_used    bias bias_mc_se se_ratio cov95  skew_z
#>       <char> <char> <num>  <int>   <num>      <num>    <num> <num>   <num>
#> 1:        x1   beta   0.8     50 -0.0050    -0.3853   0.9090  0.92  0.3072
#> 2:        x2   beta  -0.6     50 -0.0180    -1.5758   1.0296  0.92  0.3884
#> 3:     ASC_1    asc   0.5     50 -0.0259    -2.1113   1.0486  0.96  0.2840
#> 4:     ASC_2    asc  -0.5     50 -0.0066    -0.3587   0.8646  0.90  0.5034
#> 5:     ASC_3    asc   0.5     50  0.0009     0.0683   0.9421  0.94 -0.2596
#>    kurt_excess_z shapiro_p
#>            <num>     <num>
#> 1:        0.3865    0.8896
#> 2:        0.5349    0.3681
#> 3:       -0.4024    0.7106
#> 4:        0.3318    0.3830
#> 5:       -0.6189    0.5001
# }
```
