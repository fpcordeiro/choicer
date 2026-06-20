# Robust (sandwich) variance for a weighted / choice-based mixed logit fit

Recomputes the robust Huber-White sandwich variance \\V = A^{-1} B
A^{-1}\\ for a fitted mixed logit, where the bread \\A = \sum_i w_i
(-H_i)\\ is the weighted negated Hessian and the meat \\B = \sum_i w_i^2
s_i s_i'\\ is the weight-squared outer product of the per-individual
scores. This is the appropriate variance under choice-based (endogenous
stratified) / WESML weighting, where the inverse-Hessian and the
ordinary BHHH variance are invalid. It can be called on any fitted model
(e.g. one estimated with `se_method = "hessian"`) to obtain robust
standard errors post hoc, without refitting.

## Usage

``` r
wesml_vcov(object, ...)

# S3 method for class 'choicer_mxl'
wesml_vcov(object, type = c("vcov", "se"), ...)
```

## Arguments

- object:

  A fitted `choicer_mxl` object (requires `keep_data = TRUE`).

- ...:

  Unused.

- type:

  Either `"vcov"` (default) to return the variance-covariance matrix or
  `"se"` to return the standard-error vector.

## Value

A variance-covariance matrix (`type = "vcov"`) or a named numeric vector
of standard errors (`type = "se"`), in the raw parameter space.

## Details

If the stored weights are uniform (all equal), a warning is emitted: the
returned variance is then the ordinary robust (Huber-White) variance,
not a WESML-weighted variance. Refit with WESML weights for a
choice-based-sampling correction.

## See also

[`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md),
[`sample_by_choice`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md),
[`run_mxlogit`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)

## Examples

``` r
# \donttest{
library(data.table)
set.seed(1)
N <- 200L; J <- 3L
dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#>         id   alt         x1          w1
#>      <int> <int>      <num>       <num>
#>   1:     1     1 -0.6264538 -0.34106698
#>   2:     1     2  0.1836433  1.50242453
#>   3:     1     3 -0.8356286  0.52830771
#>   4:     2     1  1.5952808  0.54219136
#>   5:     2     2  0.3295078 -0.13667336
#>  ---                                   
#> 596:   199     2  0.7682782 -0.30824994
#> 597:   199     3 -0.8161606  0.01551524
#> 598:   200     1 -0.4361069 -0.44231772
#> 599:   200     2  0.9047050 -1.63800773
#> 600:   200     3 -0.7630863 -0.64140116
dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#>         id   alt         x1          w1 choice
#>      <int> <int>      <num>       <num>  <int>
#>   1:     1     1 -0.6264538 -0.34106698      0
#>   2:     1     2  0.1836433  1.50242453      0
#>   3:     1     3 -0.8356286  0.52830771      1
#>   4:     2     1  1.5952808  0.54219136      0
#>   5:     2     2  0.3295078 -0.13667336      0
#>  ---                                          
#> 596:   199     2  0.7682782 -0.30824994      1
#> 597:   199     3 -0.8161606  0.01551524      0
#> 598:   200     1 -0.4361069 -0.44231772      1
#> 599:   200     2  0.9047050 -1.63800773      0
#> 600:   200     3 -0.7630863 -0.64140116      0
fit <- run_mxlogit(dt, "id", "alt", "choice", "x1", "w1", S = 50L)
#> Optimization run time 0h:0m:0.02s
wesml_vcov(fit, "se")
#> Warning: Stored weights are uniform; wesml_vcov() returns the ordinary robust (Huber-White) variance, not a WESML-weighted variance. Refit with WESML weights for a choice-based-sampling correction.
#>         x1       L_11      ASC_2      ASC_3 
#> 0.09022091 1.50408165 0.17217370 0.17296024 
# }
```
