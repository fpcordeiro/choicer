# Prepare inputs for `mnp_gibbs()`

Prepares and validates inputs for Bayesian multinomial probit
estimation. Covariates are differenced against the base alternative, so
the design matrix has one row per (choice situation, non-base
alternative) pair. Balanced choice sets are required: every choice
situation must contain the same \\J\\ alternatives.

## Usage

``` r
prepare_mnp_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  base_alt = NULL,
  use_asc = TRUE
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

- base_alt:

  Label of the base (reference) alternative used for utility
  differencing. If `NULL` (default), the first alternative in sort order
  is used.

- use_asc:

  Logical indicating whether to include alternative-specific constants
  (one intercept per non-base alternative).

## Value

A list containing:

- `X`: Stacked differenced design matrix ((N \* p) x K), covariate
  columns first, then ASC columns when `use_asc = TRUE`.

- `y`: Integer vector of choices (0 = base alternative, j in 1..p for
  the j-th non-base alternative), one per choice situation.

- `p`: Number of utility differences (J - 1).

- `J`: Number of alternatives.

- `N`: Number of choice situations.

- `K`: Number of columns of `X`.

- `alt_mapping`: Data.table mapping alternatives to summary statistics
  (the base alternative is `alt_int = 1`).

- `base_alt`: Resolved label of the base alternative.

- `param_map`: Named list of integer index vectors (beta, asc).

- `use_asc`: Logical flag.

- `dropped_cols`: Names of columns dropped due to collinearity, if any.

- `data_spec`: List with column name metadata.

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
input <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
str(input$X)
#>  num [1:100, 1:4] -1.936 -1.008 -0.229 -0.739 -1.606 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:4] "x1" "x2" "ASC_2" "ASC_3"
input$alt_mapping
#> Key: <alt_int, alt>
#>    alt_int   alt N_OBS N_CHOICES TAKE_RATE MKT_SHARE
#>      <int> <int> <int>     <int>     <num>     <num>
#> 1:       1     1    50        17      0.34      0.34
#> 2:       2     2    50        21      0.42      0.42
#> 3:       3     3    50        12      0.24      0.24
```
