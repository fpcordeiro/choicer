# Prepare inputs for hierarchical multinomial probit estimation

Prepares and validates panel (or cross-sectional) choice data for the
hierarchical Bayesian multinomial probit with iid \\N(0, \sigma^2)\\
utility shocks. The model shares its two-level random-effect structure
with
[`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md):
respondent-level structural tastes \\\beta_i \sim N(b, W)\\ over the
`covariate_cols` (normal only — the probit keeps full conjugacy), and a
global alternative-level effect \\\delta_j = z_j'\theta + \xi_j\\,
\\\xi_j \sim N(0, \sigma_d^2)\\.

## Usage

``` r
prepare_hmnp_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  person_col = NULL,
  alt_covariate_cols = NULL,
  outside_opt_label = NULL,
  cf_residual_col = NULL,
  include_outside_option = TRUE
)
```

## Arguments

- data:

  Data frame containing choice data.

- id_col:

  Name of the column identifying choice situations (tasks). Task ids
  only need to be unique within a respondent.

- alt_col:

  Name of the column identifying alternatives.

- choice_col:

  Name of the column indicating the chosen alternative (1 = chosen, 0 =
  not chosen).

- covariate_cols:

  Vector of names of structural covariate columns (the
  random-coefficient dimensions).

- person_col:

  Name of the respondent column grouping choice situations. `NULL`
  (default) makes each choice situation its own respondent.

- alt_covariate_cols:

  Names of alternative-level covariate columns (constant within each
  alternative) forming the \\\delta\\ mean function. `NULL` (default)
  gives an intercept-only design (P = 1).

- outside_opt_label:

  Label of physical outside-option rows, removed when
  `include_outside_option = TRUE` (the outside good is implicit).

- cf_residual_col:

  Name of a first-stage residual column (control function for an
  endogenous covariate), appended to `X`. Default `NULL`.

- include_outside_option:

  Logical; if `TRUE` (default) an implicit outside option with
  systematic utility 0 is part of every choice set.

## Value

A list of class `c("choicer_data_hmnp", "list")` with the same
components as
[`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md)
(minus `rc_dist`).

## Details

The returned structure is identical to
[`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md)
(both preps share one internal engine), except there is no `rc_dist`
field. Unlike
[`prepare_mnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnp_data.md),
utilities are NOT differenced against a base alternative: the iid-shock
model works in un-differenced utility space, so unbalanced choice sets
are supported and the outside option is implicit (its latent utility is
a stochastic \\N(0, \sigma^2)\\ draw in the kernel, systematic utility
0).

## See also

[`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md)
for the component-by-component description.

## Examples

``` r
library(data.table)
set.seed(42)
N <- 20; T <- 3; J <- 4
dt <- data.table(
  pid  = rep(1:N, each = T * J),
  task = rep(seq_len(N * T), each = J),
  alt  = rep(1:J, N * T)
)
dt[, `:=`(x1 = rnorm(.N), x2 = runif(.N, -1, 1))]
#>        pid  task   alt         x1         x2
#>      <int> <int> <int>      <num>      <num>
#>   1:     1     1     1  1.3709584 -0.5341312
#>   2:     1     1     2 -0.5646982  0.1540964
#>   3:     1     1     3  0.3631284  0.6817541
#>   4:     1     1     4  0.6328626 -0.7355924
#>   5:     1     2     1  0.4042683  0.7917824
#>  ---                                        
#> 236:    20    59     4 -0.3909654  0.4483229
#> 237:    20    60     1  1.3487070 -0.2291228
#> 238:    20    60     2 -0.0227647 -0.4036954
#> 239:    20    60     3  0.2442259 -0.8850577
#> 240:    20    60     4 -0.9423717  0.3473859
dt[, choice := 0L]
#>        pid  task   alt         x1         x2 choice
#>      <int> <int> <int>      <num>      <num>  <int>
#>   1:     1     1     1  1.3709584 -0.5341312      0
#>   2:     1     1     2 -0.5646982  0.1540964      0
#>   3:     1     1     3  0.3631284  0.6817541      0
#>   4:     1     1     4  0.6328626 -0.7355924      0
#>   5:     1     2     1  0.4042683  0.7917824      0
#>  ---                                               
#> 236:    20    59     4 -0.3909654  0.4483229      0
#> 237:    20    60     1  1.3487070 -0.2291228      0
#> 238:    20    60     2 -0.0227647 -0.4036954      0
#> 239:    20    60     3  0.2442259 -0.8850577      0
#> 240:    20    60     4 -0.9423717  0.3473859      0
dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
   by = task]
#>        pid  task   alt         x1         x2 choice
#>      <int> <int> <int>      <num>      <num>  <int>
#>   1:     1     1     1  1.3709584 -0.5341312      0
#>   2:     1     1     2 -0.5646982  0.1540964      0
#>   3:     1     1     3  0.3631284  0.6817541      1
#>   4:     1     1     4  0.6328626 -0.7355924      0
#>   5:     1     2     1  0.4042683  0.7917824      1
#>  ---                                               
#> 236:    20    59     4 -0.3909654  0.4483229      1
#> 237:    20    60     1  1.3487070 -0.2291228      0
#> 238:    20    60     2 -0.0227647 -0.4036954      1
#> 239:    20    60     3  0.2442259 -0.8850577      0
#> 240:    20    60     4 -0.9423717  0.3473859      0
input <- prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "x2"),
                           person_col = "pid")
input$Ti[1:5]
#> [1] 3 3 3 3 3
input$alt_mapping
#>    alt_int   alt N_OBS N_CHOICES TAKE_RATE MKT_SHARE
#>      <int> <int> <int>     <int>     <num>     <num>
#> 1:       0    NA    60        14 0.2333333 0.2333333
#> 2:       1     1    60        15 0.2500000 0.2500000
#> 3:       2     2    60         8 0.1333333 0.1333333
#> 4:       3     3    60        10 0.1666667 0.1666667
#> 5:       4     4    60        13 0.2166667 0.2166667
```
