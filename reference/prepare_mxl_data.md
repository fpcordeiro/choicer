# Prepare inputs for `mxl_loglik_gradient_parallel()`

Prepares and validates inputs for mixed logit estimation routine.

## Usage

``` r
prepare_mxl_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  random_var_cols,
  weights = NULL,
  outside_opt_label = NULL,
  include_outside_option = FALSE,
  rc_correlation = FALSE,
  weights_col = NULL
)
```

## Arguments

- data:

  Data frame containing choice data

- id_col:

  Name of the column identifying choice situations (individuals)

- alt_col:

  Name of the column identifying alternatives

- choice_col:

  Name of the column indicating chosen alternative (1 = chosen, 0 = not
  chosen)

- covariate_cols:

  Vector of names of columns to be used as covariates

- random_var_cols:

  Vector of names of columns to be used as random variables

- weights:

  Optional vector of weights for each choice situation. If NULL, equal
  weights are used. All weights must be finite and strictly positive.

- outside_opt_label:

  Label for the outside option (if any). If NULL, no outside option is
  assumed.

- include_outside_option:

  Logical indicating whether to include an outside option in the model.

- rc_correlation:

  Logical indicating whether random coefficients are correlated. Default
  is FALSE.

- weights_col:

  Optional name of a column in `data` holding a per-row weight (constant
  within each choice situation, finite and strictly positive). Mutually
  exclusive with `weights`.

## Value

A `choicer_data_mxl` object (list) containing:

- `X`: Fixed-coefficient design matrix (sum(M) x K_x).

- `W`: Random-coefficient design matrix (sum(M) x K_w).

- `alt_idx`: Integer vector of alternative indices.

- `choice_idx`: Integer vector of chosen alternative indices.

- `M`: Integer vector with number of alternatives per choice situation.

- `N`: Number of choice situations.

- `weights`: Vector of weights.

- `include_outside_option`: Logical flag.

- `rc_correlation`: Logical flag.

- `alt_mapping`: data.table mapping alternatives to summary statistics.

- `dropped_cols`: Names of columns dropped due to collinearity, if any.

- `data_spec`: List with column-name metadata.

## Examples

``` r
library(data.table)
set.seed(42)
N <- 50; J <- 3
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#>         id   alt         x1          w1           w2
#>      <int> <int>      <num>       <num>        <num>
#>   1:     1     1  1.3709584 -0.04069848 -0.004620768
#>   2:     1     2 -0.5646982 -1.55154482  0.760242168
#>   3:     1     3  0.3631284  1.16716955  0.038990913
#>   4:     2     1  0.6328626 -0.27364570  0.735072142
#>   5:     2     2  0.4042683 -0.46784532 -0.146472627
#>  ---                                                
#> 146:    49     2  1.1133860 -0.47733551 -0.585011509
#> 147:    49     3 -0.4809928 -0.16626149  0.320957523
#> 148:    50     1 -0.4331690  0.86256338 -0.299396017
#> 149:    50     2  0.6968626  0.09734049 -0.278543083
#> 150:    50     3 -1.0563684 -1.62561674  0.546115158
dt[, choice := 0L]
#>         id   alt         x1          w1           w2 choice
#>      <int> <int>      <num>       <num>        <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848 -0.004620768      0
#>   2:     1     2 -0.5646982 -1.55154482  0.760242168      0
#>   3:     1     3  0.3631284  1.16716955  0.038990913      0
#>   4:     2     1  0.6328626 -0.27364570  0.735072142      0
#>   5:     2     2  0.4042683 -0.46784532 -0.146472627      0
#>  ---                                                       
#> 146:    49     2  1.1133860 -0.47733551 -0.585011509      0
#> 147:    49     3 -0.4809928 -0.16626149  0.320957523      0
#> 148:    50     1 -0.4331690  0.86256338 -0.299396017      0
#> 149:    50     2  0.6968626  0.09734049 -0.278543083      0
#> 150:    50     3 -1.0563684 -1.62561674  0.546115158      0
dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#>         id   alt         x1          w1           w2 choice
#>      <int> <int>      <num>       <num>        <num>  <int>
#>   1:     1     1  1.3709584 -0.04069848 -0.004620768      0
#>   2:     1     2 -0.5646982 -1.55154482  0.760242168      1
#>   3:     1     3  0.3631284  1.16716955  0.038990913      0
#>   4:     2     1  0.6328626 -0.27364570  0.735072142      1
#>   5:     2     2  0.4042683 -0.46784532 -0.146472627      0
#>  ---                                                       
#> 146:    49     2  1.1133860 -0.47733551 -0.585011509      0
#> 147:    49     3 -0.4809928 -0.16626149  0.320957523      0
#> 148:    50     1 -0.4331690  0.86256338 -0.299396017      1
#> 149:    50     2  0.6968626  0.09734049 -0.278543083      0
#> 150:    50     3 -1.0563684 -1.62561674  0.546115158      0
input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
str(input$X)
#>  num [1:150, 1] 1.371 -0.565 0.363 0.633 0.404 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr "x1"
str(input$W)
#>  num [1:150, 1:2] -0.0407 -1.5515 1.1672 -0.2736 -0.4678 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "w1" "w2"
```
