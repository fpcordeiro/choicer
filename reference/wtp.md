# Compute willingness to pay

Computes willingness-to-pay (WTP) ratios with delta-method standard
errors from a fitted choice model. For an attribute coefficient
\\\theta_k\\ and a price coefficient \\\theta_p\\, the WTP is \$\$WTP_k
= -\theta_k / \theta_p,\$\$ the marginal rate of substitution between
the attribute and price. Standard errors use the delta method with
analytic gradients \\\partial g/\partial \theta_k = -1/\theta_p\\ and
\\\partial g/\partial \theta_p = \theta_k/\theta_p^2\\, applied to the
corresponding 2x2 block of `vcov(object)`.

## Usage

``` r
wtp(object, price_var, attr_vars = NULL, level = 0.95, ...)

# S3 method for class 'choicer_fit'
wtp(object, price_var, attr_vars = NULL, level = 0.95, ...)

# S3 method for class 'choicer_mxl'
wtp(object, price_var, attr_vars = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model object (`choicer_mnl`, `choicer_mxl`, or `choicer_nl`).

- price_var:

  Name of the price variable. Must be a fixed-coefficient variable (a
  column of the design matrix `X`).

- attr_vars:

  Character vector of attributes to report. Defaults to all
  fixed-coefficient variables other than `price_var` (plus, for mixed
  logit with `rc_mean = TRUE`, all random coefficients). ASC names (e.g.
  `"ASC_2"`) may also be supplied; the WTP of an ASC is \\-ASC_j /
  \theta_p\\.

- level:

  Confidence level for the normal-approximation interval \\Estimate \pm
  z\_{1-(1-level)/2} \times SE\\. Default 0.95.

- ...:

  Additional arguments passed to methods.

## Value

A data.frame of class `choicer_wtp` with one row per attribute and
columns `Estimate`, `Std_Error`, `z_value`, `CI_lower`, `CI_upper`.
Attributes `price_var` and `level` record the inputs; `median_rows`
lists rows that are median (rather than mean) WTP. Standard errors are
NA when the variance-covariance matrix is unavailable.

## Details

For mixed logit models, random coefficients are included via their
estimated location parameters. The package's log-normal random
coefficient is the *shifted* log-normal \\\beta_k = \exp(\mu_k) +
\exp((L\eta)\_k)\\ (see
[`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)),
so:

- Normal random coefficient \\k\\ (`rc_mean = TRUE`): mean WTP \\-\mu_k
  / \theta_p\\, labeled `Mu_x`.

- Log-normal random coefficient \\k\\ (`rc_mean = TRUE`): **median** WTP
  \\-(\exp(\mu_k) + 1) / \theta_p\\, since the median of
  \\\exp((L\eta)\_k)\\ is 1. (The mean, \\\exp(\mu_k) +
  \exp(\sigma_k^2/2)\\, is highly sensitive to the estimated variance;
  the median is the more robust summary.) These rows are labeled by the
  attribute name and flagged as medians when printed.

- Log-normal random coefficient with `rc_mean = FALSE`: \\\beta_k =
  \exp((L\eta)\_k)\\ has median 1, so the median WTP is \\-1/\theta_p\\
  with uncertainty driven solely by \\\theta_p\\.

Normal random coefficients with `rc_mean = FALSE` have mean 0 by
construction and are excluded from the table.

The price variable must have a *fixed* coefficient. A random price
coefficient is rejected: the ratio of two random coefficients generally
has no finite moments (the denominator has positive density at 0), so
mean or median WTP computed from location parameters would be
meaningless. Use a fixed price coefficient, or estimate the model in WTP
space.

## Examples

``` r
# \donttest{
library(data.table)
sim <- simulate_mnl_data(N = 1000, J = 4, beta = c(0.8, -0.6), seed = 123,
                         outside_option = FALSE, vary_choice_set = FALSE)
fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s
# treat x2 as the price variable
wtp(fit, price_var = "x2")
#> Willingness to pay (WTP), price variable: 'x2' (95% CI)
#>    Estimate Std_Error z_value CI_lower CI_upper
#> x1     1.59    0.2336   6.809    1.132    2.048
wtp(fit, price_var = "x2", attr_vars = c("x1", "ASC_2"), level = 0.90)
#> Willingness to pay (WTP), price variable: 'x2' (90% CI)
#>       Estimate Std_Error z_value CI_lower CI_upper
#> x1       1.590    0.2336   6.809    1.206    1.974
#> ASC_2   -2.027    0.3157  -6.420   -2.546   -1.507
# }
```
