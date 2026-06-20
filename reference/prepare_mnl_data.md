# Prepare inputs for `mnl_loglik_gradient_parallel()`

Prepares and validates inputs for multinomial logit estimation routine.

## Usage

``` r
prepare_mnl_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  weights = NULL,
  outside_opt_label = NULL,
  include_outside_option = FALSE
)
```

## Arguments

- data:

  Data frame containing choice data.

- id_col:

  Name of the column identifying choice situations (individuals).

- alt_col:

  Name of the column identifying alternatives.

- choice_col:

  Name of the column indicating chosen alternative (1 = chosen, 0 = not
  chosen).

- covariate_cols:

  Vector of names of columns to be used as covariates.

- weights:

  Optional vector of weights for each choice situation. If `NULL`, equal
  weights are used.

- outside_opt_label:

  Label for the outside option (if any). If `NULL`, no outside option is
  assumed.

- include_outside_option:

  Logical indicating whether to include an outside option in the model.

## Value

A list containing:

- `X`: Design matrix (sum(M) x K).

- `alt_idx`: Integer vector of alternative indices.

- `choice_idx`: Integer vector of chosen alternative indices.

- `M`: Integer vector with number of alternatives per choice situation.

- `N`: Number of choice situations.

- `weights`: Vector of weights.

- `include_outside_option`: Logical flag.

- `alt_mapping`: Data.table mapping alternatives to summary statistics.

- `dropped_cols`: Names of columns dropped due to collinearity, if any.

## Examples

``` r
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
input <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
str(input$X)
#>  num [1:150, 1:2] 1.371 -0.565 0.363 0.633 0.404 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "x1" "x2"
input$alt_mapping
#> Key: <alt_int, alt>
#>    alt_int   alt N_OBS N_CHOICES TAKE_RATE MKT_SHARE
#>      <int> <int> <int>     <int>     <num>     <num>
#> 1:       1     1    50        17      0.34      0.34
#> 2:       2     2    50        21      0.42      0.42
#> 3:       3     3    50        12      0.24      0.24
```
