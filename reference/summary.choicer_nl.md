# Summary for nested logit model

Triggers lazy Hessian computation if standard errors have not been
computed yet.

## Usage

``` r
# S3 method for class 'choicer_nl'
summary(object, gof = TRUE, ...)
```

## Arguments

- object:

  A choicer_nl object.

- gof:

  Logical; compute goodness-of-fit measures (McFadden R-squared, hit
  rate) for the summary footer. Involves an in-sample prediction pass
  (for mixed logit, a full simulation over draws); set to FALSE to skip.

- ...:

  Additional arguments (ignored).

## Value

A summary.choicer_nl object (includes a `gof` element with
goodness-of-fit measures from
[`gof`](https://fpcordeiro.github.io/choicer/reference/gof.md); its
fields are NA when the model was fitted with `keep_data = FALSE`).

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 4
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1         x2
#>      <int> <int>      <num>      <num>
#>   1:     1     1  1.3709584 -2.0009292
#>   2:     1     2 -0.5646982  0.3337772
#>   3:     1     3  0.3631284  1.1713251
#>   4:     1     4  0.6328626  2.0595392
#>   5:     2     1  0.4042683 -1.3768616
#>  ---                                  
#> 196:    49     4  1.0857749  1.0965134
#> 197:    50     1  0.4037749  0.4420131
#> 198:    50     2  0.5864875  0.2410163
#> 199:    50     3  1.8152284 -0.2556077
#> 200:    50     4  0.1288214  0.9310329
dt[, nest := ifelse(alt <= 2, "A", "B")]
#>         id   alt         x1         x2   nest
#>      <int> <int>      <num>      <num> <char>
#>   1:     1     1  1.3709584 -2.0009292      A
#>   2:     1     2 -0.5646982  0.3337772      A
#>   3:     1     3  0.3631284  1.1713251      B
#>   4:     1     4  0.6328626  2.0595392      B
#>   5:     2     1  0.4042683 -1.3768616      A
#>  ---                                         
#> 196:    49     4  1.0857749  1.0965134      B
#> 197:    50     1  0.4037749  0.4420131      A
#> 198:    50     2  0.5864875  0.2410163      A
#> 199:    50     3  1.8152284 -0.2556077      B
#> 200:    50     4  0.1288214  0.9310329      B
dt[, choice := 0L]
#>         id   alt         x1         x2   nest choice
#>      <int> <int>      <num>      <num> <char>  <int>
#>   1:     1     1  1.3709584 -2.0009292      A      0
#>   2:     1     2 -0.5646982  0.3337772      A      0
#>   3:     1     3  0.3631284  1.1713251      B      0
#>   4:     1     4  0.6328626  2.0595392      B      0
#>   5:     2     1  0.4042683 -1.3768616      A      0
#>  ---                                                
#> 196:    49     4  1.0857749  1.0965134      B      0
#> 197:    50     1  0.4037749  0.4420131      A      0
#> 198:    50     2  0.5864875  0.2410163      A      0
#> 199:    50     3  1.8152284 -0.2556077      B      0
#> 200:    50     4  0.1288214  0.9310329      B      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1         x2   nest choice
#>      <int> <int>      <num>      <num> <char>  <int>
#>   1:     1     1  1.3709584 -2.0009292      A      0
#>   2:     1     2 -0.5646982  0.3337772      A      0
#>   3:     1     3  0.3631284  1.1713251      B      0
#>   4:     1     4  0.6328626  2.0595392      B      1
#>   5:     2     1  0.4042683 -1.3768616      A      0
#>  ---                                                
#> 196:    49     4  1.0857749  1.0965134      B      1
#> 197:    50     1  0.4037749  0.4420131      A      0
#> 198:    50     2  0.5864875  0.2410163      A      0
#> 199:    50     3  1.8152284 -0.2556077      B      0
#> 200:    50     4  0.1288214  0.9310329      B      1
fit <- run_nestlogit(
  data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = c("x1", "x2"), nest_col = "nest"
)
#> Optimization run time 0h:0m:0.01s
summary(fit)
#> Nested Logit (NL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           0.368065   0.337997   1.0890  2.76e-01  
#> x2           0.658149   0.371606   1.7711  7.65e-02  
#> Lambda_1   158.730748 9738.696580   0.0163  9.87e-01  
#> Lambda_2     2.096862   1.990404   1.0535  2.92e-01  
#> ASC_2       13.362918 829.500434   0.0161  9.87e-01  
#> ASC_3      116.246676 7172.663321   0.0162  9.87e-01  
#> ASC_4      114.565371 7172.763864   0.0160  9.87e-01  
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Std. Errors: Analytical Hessian 
#> Log-likelihood: -63.8964 
#> AIC: 141.793  | BIC: 155.177 
#> McFadden R2: 0.078 (adj: -0.023) | Hit rate: 0.360 
#> N: 50  | Parameters: 7 
#> Optimization time: 0.01 s
#> Convergence: 3 ( NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. )
# }
```
