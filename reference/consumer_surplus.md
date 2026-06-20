# Expected consumer surplus

Computes the expected consumer surplus per choice situation (Train 2009,
Ch. 3): \$\$E\[CS_i\] = \frac{logsum_i}{-\alpha},\$\$ where \\logsum_i\\
is the expected maximum utility (see
[`logsum`](https://fpcordeiro.github.io/choicer/reference/logsum.md))
and \\\alpha\\ is the (fixed) price coefficient, so that \\-\alpha\\ is
the marginal utility of income. The formula assumes *no income effects*:
utility is linear in price, and the marginal utility of income is
constant across the price changes considered.

## Usage

``` r
consumer_surplus(
  object,
  price_var,
  newdata = NULL,
  level = 0.95,
  weights = NULL,
  ...
)

# S3 method for class 'choicer_mnl'
consumer_surplus(
  object,
  price_var,
  newdata = NULL,
  level = 0.95,
  weights = NULL,
  ...
)

# S3 method for class 'choicer_mxl'
consumer_surplus(
  object,
  price_var,
  newdata = NULL,
  level = 0.95,
  weights = NULL,
  ...
)

# S3 method for class 'choicer_nl'
consumer_surplus(
  object,
  price_var,
  newdata = NULL,
  level = 0.95,
  weights = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object (`choicer_mnl`, `choicer_mxl`, or `choicer_nl`).

- price_var:

  Name of the price variable. Must be a fixed-coefficient variable (a
  column of the design matrix `X`).

- newdata:

  Optional counterfactual data (data.frame or list), as in
  [`logsum`](https://fpcordeiro.github.io/choicer/reference/logsum.md)
  and [`predict()`](https://rdrr.io/r/stats/predict.html). When `NULL`
  (default), the data stored at fit time is used (requires
  `keep_data = TRUE`).

- level:

  Confidence level for the normal-approximation interval around the mean
  CS (MNL only). Default 0.95.

- weights:

  Optional numeric vector with one weight per choice situation, used for
  the mean CS (and its SE), as in
  [`predict()`](https://rdrr.io/r/stats/predict.html): for a data.frame
  `newdata`, one weight per id in order of first appearance. Defaults to
  equal weights. Ignored when `newdata` is `NULL` (the stored fit
  weights apply).

- ...:

  Additional arguments passed to methods.

## Value

A `choicer_cs` object: a list with `cs` (per-choice- situation surplus,
length N), `mean_cs` (weighted mean), `se_mean_cs` (delta-method SE; NA
for MXL/NL or when the variance-covariance matrix is unavailable), `ci`
(confidence interval for the mean), `price_var`, `level`, and `n`.

## Details

Consumer surplus *levels* inherit the additive utility normalization (in
particular the ASC normalization), so the level is only defined up to a
constant; *differences* in CS between scenarios — e.g.
`consumer_surplus(fit, "price", newdata = scenario)` minus the baseline
— are the economically meaningful quantity.

For MNL fits, a delta-method standard error of the weighted mean CS is
reported (weights are the stored fit weights, or the resolved `newdata`
weights). For MXL and NL fits only point estimates are returned
(`se_mean_cs = NA`): the delta method for the simulated MXL logsum and
the nested logsum is deferred; simulation-based intervals (Krinsky-Robb:
resample coefficients from their asymptotic distribution and recompute
the mean CS) are a practical alternative.

The price variable must have a *fixed* coefficient. For mixed logit a
random price coefficient is rejected (as in
[`wtp`](https://fpcordeiro.github.io/choicer/reference/wtp.md)): with a
random denominator \\1/(-\alpha)\\ generally has no finite moments.

## References

Train, K. (2009). *Discrete Choice Methods with Simulation*, 2nd ed.,
Ch. 3. Cambridge University Press.

## See also

[`logsum`](https://fpcordeiro.github.io/choicer/reference/logsum.md),
[`wtp`](https://fpcordeiro.github.io/choicer/reference/wtp.md)

## Examples

``` r
# \donttest{
library(data.table)
sim <- simulate_mnl_data(N = 1000, J = 3, beta = c(0.8, -0.6), seed = 123,
                         outside_option = FALSE, vary_choice_set = FALSE)
fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s

# treat x2 as the price variable
cs0 <- consumer_surplus(fit, price_var = "x2")
cs0
#> Consumer surplus, price variable: 'x2'
#>   Mean CS: 1.692 
#>   SE (delta method): 0.2281 
#>   95% CI: [1.245, 2.139]
#>   N: 1000 
#> Note: CS levels depend on the utility (ASC) normalization; differences between scenarios are the meaningful quantity.

# Change in consumer surplus from a price increase on alternative 2:
# levels depend on the ASC normalization, differences do not.
dt_cf <- copy(sim$data)[alt == 2, x2 := x2 + 0.5]
cs1 <- consumer_surplus(fit, price_var = "x2", newdata = dt_cf)
delta_cs <- cs1$mean_cs - cs0$mean_cs
delta_cs  # negative: the price increase lowers expected surplus
#> [1] -0.07643608
# }
```
