# Diversion ratios for mixed logit model

Computes the attribute-based diversion ratio matrix. Entry (k, j) is the
fraction of demand lost by alternative j that is captured by alternative
k when a marginal change in alternative j's `wrt_var` attribute reduces
s_j.

## Usage

``` r
# S3 method for class 'choicer_mxl'
diversion_ratios(object, wrt_var, is_random_coef = FALSE, ...)
```

## Arguments

- object:

  A `choicer_mxl` object fitted with `keep_data = TRUE`.

- wrt_var:

  Variable used to perturb alternative j's utility: a column name
  (character) or 1-based index. Indexes into X columns for fixed
  coefficients, or W columns for random coefficients (when
  `is_random_coef = TRUE`).

- is_random_coef:

  Logical. `TRUE` if the variable has a random coefficient (is in W),
  `FALSE` if fixed (in X). Default `FALSE`.

- ...:

  Additional arguments (ignored).

## Value

A J x J diversion ratio matrix with alternative labels. Cross-products
are averaged across simulation draws inside the integration to avoid
Jensen-style bias.

## Details

Unlike MNL, the MXL diversion ratio depends on which variable is
perturbed: the realised coefficient \\\beta\_{ik}^s\\ varies across
individuals and draws and does not cancel in the ratio. For a variable
with a fixed coefficient the result is independent of the variable
(\\\beta\\ cancels); for a random-coefficient variable it is not.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 50; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#>         id   alt         x1          w1
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
#>         id   alt         x1          w1 choice
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
#>         id   alt         x1          w1 choice
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
fit <- run_mxlogit(
  data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = "x1", random_var_cols = "w1", S = 50L
)
#> Optimization run time 0h:0m:0.01s
diversion_ratios(fit, "x1")
#>           1         2         3
#> 1 0.0000000 0.6361155 0.4677707
#> 2 0.6654416 0.0000000 0.5322293
#> 3 0.3345584 0.3638845 0.0000000
diversion_ratios(fit, "w1", is_random_coef = TRUE)
#>           1         2         3
#> 1 0.0000000 0.5102239 0.7924648
#> 2 0.2143425 0.0000000 0.2075352
#> 3 0.7856575 0.4897761 0.0000000
# }
```
