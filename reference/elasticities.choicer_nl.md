# Elasticities for nested logit model

Elasticities for nested logit model

## Usage

``` r
# S3 method for class 'choicer_nl'
elasticities(object, elast_var, ...)
```

## Arguments

- object:

  A `choicer_nl` object fitted with `keep_data = TRUE`.

- elast_var:

  Variable for elasticity computation: a column name (character) or
  1-based index into the design matrix X.

- ...:

  Additional arguments (ignored).

## Value

A J x J elasticity matrix with alternative labels.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 4
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, nest := rep(c(1L, 1L, 2L, 2L), N)]
#>         id   alt  nest
#>      <int> <int> <int>
#>   1:     1     1     1
#>   2:     1     2     1
#>   3:     1     3     2
#>   4:     1     4     2
#>   5:     2     1     1
#>  ---                  
#> 196:    49     4     2
#> 197:    50     1     1
#> 198:    50     2     1
#> 199:    50     3     2
#> 200:    50     4     2
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt  nest         x1         x2
#>      <int> <int> <int>      <num>      <num>
#>   1:     1     1     1  1.3709584 -2.0009292
#>   2:     1     2     1 -0.5646982  0.3337772
#>   3:     1     3     2  0.3631284  1.1713251
#>   4:     1     4     2  0.6328626  2.0595392
#>   5:     2     1     1  0.4042683 -1.3768616
#>  ---                                        
#> 196:    49     4     2  1.0857749  1.0965134
#> 197:    50     1     1  0.4037749  0.4420131
#> 198:    50     2     1  0.5864875  0.2410163
#> 199:    50     3     2  1.8152284 -0.2556077
#> 200:    50     4     2  0.1288214  0.9310329
dt[, choice := 0L]
#>         id   alt  nest         x1         x2 choice
#>      <int> <int> <int>      <num>      <num>  <int>
#>   1:     1     1     1  1.3709584 -2.0009292      0
#>   2:     1     2     1 -0.5646982  0.3337772      0
#>   3:     1     3     2  0.3631284  1.1713251      0
#>   4:     1     4     2  0.6328626  2.0595392      0
#>   5:     2     1     1  0.4042683 -1.3768616      0
#>  ---                                               
#> 196:    49     4     2  1.0857749  1.0965134      0
#> 197:    50     1     1  0.4037749  0.4420131      0
#> 198:    50     2     1  0.5864875  0.2410163      0
#> 199:    50     3     2  1.8152284 -0.2556077      0
#> 200:    50     4     2  0.1288214  0.9310329      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt  nest         x1         x2 choice
#>      <int> <int> <int>      <num>      <num>  <int>
#>   1:     1     1     1  1.3709584 -2.0009292      0
#>   2:     1     2     1 -0.5646982  0.3337772      0
#>   3:     1     3     2  0.3631284  1.1713251      0
#>   4:     1     4     2  0.6328626  2.0595392      1
#>   5:     2     1     1  0.4042683 -1.3768616      0
#>  ---                                               
#> 196:    49     4     2  1.0857749  1.0965134      1
#> 197:    50     1     1  0.4037749  0.4420131      0
#> 198:    50     2     1  0.5864875  0.2410163      0
#> 199:    50     3     2  1.8152284 -0.2556077      0
#> 200:    50     4     2  0.1288214  0.9310329      1
fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#> Optimization run time 0h:0m:0.01s
elasticities(fit, "x1")
#>              1            2            3            4
#> 1  0.001361410 -0.002540128  0.001811961 -0.013019967
#> 2  0.001254278 -0.002616362  0.001811961 -0.013019967
#> 3 -0.006952266  0.003469046 -0.043995956 -0.001679743
#> 4 -0.006952266  0.003469046 -0.015185111  0.005494552
# }
```
