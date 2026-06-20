# Goodness of fit for a fitted choice model

Computes McFadden's pseudo R-squared (plain and adjusted) and the
in-sample hit rate for a fitted model.

## Usage

``` r
gof(object, null = c("equal_shares", "market_shares"), ...)

# S3 method for class 'choicer_fit'
gof(object, null = c("equal_shares", "market_shares"), ...)
```

## Arguments

- object:

  A fitted model object (`choicer_mnl`, `choicer_mxl`, or `choicer_nl`).

- null:

  Null model for the pseudo R-squared: `"equal_shares"` (default) or
  `"market_shares"`.

- ...:

  Additional arguments passed to methods.

## Value

A `choicer_gof` object: a list with `loglik`, `loglik_null`, `null`,
`mcfadden_r2`, `mcfadden_r2_adj`, `hit_rate`, `nobs`, and `n_params`.

## Details

Two null models are available for the pseudo R-squared \\R^2 = 1 - LL /
LL_0\\ (adjusted: \\R^2\_{adj} = 1 - (LL - K) / LL_0\\ with \\K\\ the
number of estimated parameters):

- `"equal_shares"` (default): every alternative in individual \\i\\'s
  choice set is equally likely, so \\LL_0 = -\sum_i w_i \log(M_i +
  1\_{outside})\\. This is exact for unbalanced choice sets and
  arbitrary weights.

- `"market_shares"`: the maximized log-likelihood of an ASC-only model,
  \\LL_0 = \sum_j N_j \log(s_j)\\ with \\N_j\\ the choice counts and
  \\s_j\\ the observed market shares (including the outside option when
  present). This closed form is valid only for balanced choice sets and
  uniform weights; otherwise an error suggests refitting an ASC-only
  model.

The hit rate is the weighted share of individuals whose observed choice
has the highest predicted probability. When the model includes an
outside option, the outside good competes for the predicted maximum (its
probability is \\1 - \sum_j p\_{ij}\\), and an individual predicted to
choose the outside good is a hit when they actually did.

Both the null log-likelihood and the hit rate require the stored
estimation data; models fitted with `keep_data = FALSE` return NA fields
with a message.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1          x2
#>      <int> <int>      <num>       <num>
#>   1:     1     1  1.3709584 -0.04069848
#>   2:     1     2 -0.5646982 -1.55154482
#>   3:     1     3  0.3631284  1.16716955
#>   4:     2     1  0.6328626 -0.27364570
#>   5:     2     2  0.4042683 -0.46784532
#>  ---                                   
#> 146:    49     2  1.1133860 -0.47733551
#> 147:    49     3 -0.4809928 -0.16626149
#> 148:    50     1 -0.4331690  0.86256338
#> 149:    50     2  0.6968626  0.09734049
#> 150:    50     3 -1.0563684 -1.62561674
dt[, choice := 0L]
#>         id   alt         x1          x2 choice
#>      <int> <int>      <num>       <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848      0
#>   2:     1     2 -0.5646982 -1.55154482      0
#>   3:     1     3  0.3631284  1.16716955      0
#>   4:     2     1  0.6328626 -0.27364570      0
#>   5:     2     2  0.4042683 -0.46784532      0
#>  ---                                          
#> 146:    49     2  1.1133860 -0.47733551      0
#> 147:    49     3 -0.4809928 -0.16626149      0
#> 148:    50     1 -0.4331690  0.86256338      0
#> 149:    50     2  0.6968626  0.09734049      0
#> 150:    50     3 -1.0563684 -1.62561674      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1          x2 choice
#>      <int> <int>      <num>       <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848      0
#>   2:     1     2 -0.5646982 -1.55154482      0
#>   3:     1     3  0.3631284  1.16716955      1
#>   4:     2     1  0.6328626 -0.27364570      0
#>   5:     2     2  0.4042683 -0.46784532      0
#>  ---                                          
#> 146:    49     2  1.1133860 -0.47733551      0
#> 147:    49     3 -0.4809928 -0.16626149      0
#> 148:    50     1 -0.4331690  0.86256338      1
#> 149:    50     2  0.6968626  0.09734049      0
#> 150:    50     3 -1.0563684 -1.62561674      0
fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s
gof(fit)
#> Goodness of fit
#>   Log-likelihood:       -52.4423 
#>   Null log-likelihood: -54.9306 (equal shares)
#>   McFadden R2:          0.0453 (adj: -0.0275)
#>   Hit rate:             0.4600
#>   N: 50  | Parameters: 4 
gof(fit, null = "market_shares")
#> Goodness of fit
#>   Log-likelihood:       -52.4423 
#>   Null log-likelihood: -53.6827 (market shares)
#>   McFadden R2:          0.0231 (adj: -0.0514)
#>   Hit rate:             0.4600
#>   N: 50  | Parameters: 4 
# }
```
