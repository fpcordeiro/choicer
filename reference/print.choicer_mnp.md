# Print a choicer_mnp object

Prints a brief summary of the fitted Bayesian multinomial probit model.

## Usage

``` r
# S3 method for class 'choicer_mnp'
print(x, ...)
```

## Arguments

- x:

  A choicer_mnp object.

- ...:

  Additional arguments (ignored).

## Value

The object invisibly.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 100; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt          x1           x2
#>      <int> <int>       <num>        <num>
#>   1:     1     1  1.37095845 -0.004620768
#>   2:     1     2 -0.56469817  0.760242168
#>   3:     1     3  0.36312841  0.038990913
#>   4:     2     1  0.63286260  0.735072142
#>   5:     2     2  0.40426832 -0.146472627
#>  ---                                     
#> 296:    99     2 -0.47733551  0.160327395
#> 297:    99     3 -0.16626149 -0.433641942
#> 298:   100     1  0.86256338  1.537412419
#> 299:   100     2  0.09734049 -2.170246577
#> 300:   100     3 -1.62561674  1.027004619
dt[, choice := 0L]
#>         id   alt          x1           x2 choice
#>      <int> <int>       <num>        <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768      0
#>   2:     1     2 -0.56469817  0.760242168      0
#>   3:     1     3  0.36312841  0.038990913      0
#>   4:     2     1  0.63286260  0.735072142      0
#>   5:     2     2  0.40426832 -0.146472627      0
#>  ---                                            
#> 296:    99     2 -0.47733551  0.160327395      0
#> 297:    99     3 -0.16626149 -0.433641942      0
#> 298:   100     1  0.86256338  1.537412419      0
#> 299:   100     2  0.09734049 -2.170246577      0
#> 300:   100     3 -1.62561674  1.027004619      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt          x1           x2 choice
#>      <int> <int>       <num>        <num>  <int>
#>   1:     1     1  1.37095845 -0.004620768      0
#>   2:     1     2 -0.56469817  0.760242168      0
#>   3:     1     3  0.36312841  0.038990913      1
#>   4:     2     1  0.63286260  0.735072142      0
#>   5:     2     2  0.40426832 -0.146472627      1
#>  ---                                            
#> 296:    99     2 -0.47733551  0.160327395      0
#> 297:    99     3 -0.16626149 -0.433641942      0
#> 298:   100     1  0.86256338  1.537412419      0
#> 299:   100     2  0.09734049 -2.170246577      0
#> 300:   100     3 -1.62561674  1.027004619      1
fit <- run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                    mcmc = list(R = 300, burn = 100))
#> MCMC run time 0h:0m:0.01s
print(fit)
#> Bayesian Multinomial Probit (MNP) model
#>   N obs: 100  | Parameters: 4 
#>   Posterior draws kept: 200 (R = 300, burn = 100, thin = 1)
#>   Base alternative: 1 
#>   Estimates are posterior means of identified parameters (beta / sqrt(sigma_11)).
# }
```
