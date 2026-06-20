# Prepare inputs for nested logit estimation

Validates inputs, builds design matrices, and constructs nest structure
for nested logit estimation. Calls
[`prepare_mnl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
internally for base data preparation, then adds nest-specific fields.

## Usage

``` r
prepare_nl_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  nest_col,
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

- nest_col:

  Name of the column mapping each alternative to its nest. Every
  alternative must belong to exactly one nest.

- weights:

  Optional vector of weights for each choice situation. If `NULL`, equal
  weights are used.

- outside_opt_label:

  Label for the outside option (if any). If `NULL`, no outside option is
  assumed.

- include_outside_option:

  Logical indicating whether to include an outside option in the model.

## Value

A `choicer_data_nl` object (list) containing:

- All fields from
  [`prepare_mnl_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
  (`X`, `alt_idx`, `choice_idx`, `M`, `N`, `weights`,
  `include_outside_option`, `alt_mapping`, `dropped_cols`).

- `nest_idx`: Integer vector of length J mapping each alternative (in
  `alt_mapping` row order) to its nest.

- `data_spec`: List with column name metadata including `nest_col`.

## Examples

``` r
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
input <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
input$nest_idx
#> [1] 1 1 2 2
input$alt_mapping
#> Key: <alt_int, alt>
#>    alt_int   alt N_OBS N_CHOICES TAKE_RATE MKT_SHARE
#>      <int> <int> <int>     <int>     <num>     <num>
#> 1:       1     1    50        10      0.20      0.20
#> 2:       2     2    50        11      0.22      0.22
#> 3:       3     3    50        20      0.40      0.40
#> 4:       4     4    50         9      0.18      0.18
```
