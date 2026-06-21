# Summary for multinomial logit model

Computes and returns a coefficient summary table with standard errors,
z-values, p-values, and significance codes. Triggers lazy Hessian
computation if standard errors have not been computed yet.

## Usage

``` r
# S3 method for class 'choicer_mnl'
summary(object, gof = TRUE, ...)
```

## Arguments

- object:

  A choicer_mnl object.

- gof:

  Logical; compute goodness-of-fit measures (McFadden R-squared, hit
  rate) for the summary footer. Involves an in-sample prediction pass
  (for mixed logit, a full simulation over draws); set to FALSE to skip.

- ...:

  Additional arguments (ignored).

## Value

A summary.choicer_mnl object (list with coefficients table and metadata,
including a `gof` element with goodness-of-fit measures from
[`gof`](https://fpcordeiro.github.io/choicer/reference/gof.md); its
fields are NA when the model was fitted with `keep_data = FALSE`).

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
summary(fit)
#> Multinomial Logit (MNL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1          -0.105496   0.182728  -0.5773  5.64e-01  
#> x2          -0.229495   0.175867  -1.3049  1.92e-01  
#> ASC_2        0.188575   0.331919   0.5681  5.70e-01  
#> ASC_3       -0.344221   0.383461  -0.8977  3.69e-01  
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -52.4423 
#> AIC: 112.885  | BIC: 120.533 
#> McFadden R2: 0.045 (adj: -0.028) | Hit rate: 0.460 
#> N: 50  | Parameters: 4 
#> Optimization time: 0 s
#> Convergence: 1 ( NLOPT_SUCCESS: Generic success return value. )
# }
```
