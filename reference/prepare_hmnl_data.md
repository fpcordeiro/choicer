# Prepare inputs for hierarchical multinomial logit estimation

Prepares and validates panel (or cross-sectional) choice data for the
hierarchical Bayesian multinomial logit. The model has two random-effect
levels: respondent-level structural tastes \\\beta_i \sim N(b, W)\\ over
the `covariate_cols`, and a global alternative-level effect \\\delta_j =
z_j'\theta + \xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\, with
mean-function design \\z_j\\ built from `alt_covariate_cols`.

## Usage

``` r
prepare_hmnl_data(
  data,
  id_col,
  alt_col,
  choice_col,
  covariate_cols,
  person_col = NULL,
  alt_covariate_cols = NULL,
  outside_opt_label = NULL,
  cf_residual_col = NULL,
  include_outside_option = TRUE,
  rc_dist = NULL
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

- rc_dist:

  Integer vector, one entry per column of `covariate_cols`: `0` for a
  normal random coefficient, `1` for log-normal (the coefficient enters
  utility as `exp(beta_ik)`; hierarchy normal on the log scale). Default
  `NULL` is all-normal. Automatically aligned through dropped columns; a
  `cf_residual_col` coordinate is always normal.

## Value

A list of class `c("choicer_data_hmnl", "list")` containing:

- `X`: Structural design matrix (total_rows x K_struct), no ASC columns;
  `cf_residual_col` last when supplied.

- `alt_of_row`: Integer alternative code per row (`1..J`).

- `alt_idx`: Alias of `alt_of_row` for the pooled-MLE init.

- `Z`: Alternative-level design (J x P), intercept first.

- `M`: Inside alternatives per choice situation.

- `choice_pos`: 1-based within-task position of the chosen row; `0` =
  outside option chosen.

- `Ti`: Choice situations per respondent.

- `person_ids`, `N_persons`, `n_tasks`, `J`, `K_struct`, `P`.

- `include_outside_option`: Logical flag.

- `alt_mapping`: Data.table mapping alternatives to summary statistics
  (outside option is `alt_int = 0`).

- `param_map`: Named list of index vectors (`beta`, `theta`), robust to
  collinearity drops.

- `rc_dist`: Integer vector aligned with the columns of `X`.

- `dropped_cols`, `dropped_z_cols`: Dropped column names, if any.

- `data_spec`: Column-name metadata (incl. `person_col`,
  `outside_opt_label`, `cf_residual_col`, `alt_covariate_cols`).

## Details

**Structure.** The design matrix `X` carries structural covariates only
— no alternative-specific-constant dummies. The alternative effect
\\\delta_j\\ is indexed by `alt_of_row` (integer codes `1..J`), so
memory and compute scale with the number of rows, not with `J` extra
design columns.

**Outside option.** With `include_outside_option = TRUE` (the default)
the outside good is modelled implicitly, following the
[`prepare_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
convention: physical outside rows (identified by `outside_opt_label`)
are removed, the estimation kernels add the outside term (systematic
utility 0), and a choice situation whose inside rows are all `0` in
`choice_col` is coded as "outside chosen" (`choice_pos = 0`). The
outside option anchors the location of \\\delta\\ (mean utility relative
to the outside good).

**Cross-section vs panel.** `person_col` groups choice situations into
respondents sharing one \\\beta_i\\. With `person_col = NULL` (default)
every choice situation is its own respondent (`Ti` all 1) — the
cross-sectional random-coefficients mode.

**Control function.** `cf_residual_col` (a user-supplied first-stage
residual, Petrin & Train 2010) is appended to `X` as an ordinary
covariate; its provenance is recorded in `data_spec`. The first stage is
NOT run here — supplying a valid residual is the user's responsibility.

## See also

[`prepare_hmnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnp_data.md)
for the hierarchical probit counterpart.

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
dt[, quality := 0.1 * alt]                # alternative-level covariate
#>        pid  task   alt         x1         x2 quality
#>      <int> <int> <int>      <num>      <num>   <num>
#>   1:     1     1     1  1.3709584 -0.5341312     0.1
#>   2:     1     1     2 -0.5646982  0.1540964     0.2
#>   3:     1     1     3  0.3631284  0.6817541     0.3
#>   4:     1     1     4  0.6328626 -0.7355924     0.4
#>   5:     1     2     1  0.4042683  0.7917824     0.1
#>  ---                                                
#> 236:    20    59     4 -0.3909654  0.4483229     0.4
#> 237:    20    60     1  1.3487070 -0.2291228     0.1
#> 238:    20    60     2 -0.0227647 -0.4036954     0.2
#> 239:    20    60     3  0.2442259 -0.8850577     0.3
#> 240:    20    60     4 -0.9423717  0.3473859     0.4
dt[, choice := 0L]
#>        pid  task   alt         x1         x2 quality choice
#>      <int> <int> <int>      <num>      <num>   <num>  <int>
#>   1:     1     1     1  1.3709584 -0.5341312     0.1      0
#>   2:     1     1     2 -0.5646982  0.1540964     0.2      0
#>   3:     1     1     3  0.3631284  0.6817541     0.3      0
#>   4:     1     1     4  0.6328626 -0.7355924     0.4      0
#>   5:     1     2     1  0.4042683  0.7917824     0.1      0
#>  ---                                                       
#> 236:    20    59     4 -0.3909654  0.4483229     0.4      0
#> 237:    20    60     1  1.3487070 -0.2291228     0.1      0
#> 238:    20    60     2 -0.0227647 -0.4036954     0.2      0
#> 239:    20    60     3  0.2442259 -0.8850577     0.3      0
#> 240:    20    60     4 -0.9423717  0.3473859     0.4      0
# leave some tasks all-zero: outside option chosen
dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
   by = task]
#>        pid  task   alt         x1         x2 quality choice
#>      <int> <int> <int>      <num>      <num>   <num>  <int>
#>   1:     1     1     1  1.3709584 -0.5341312     0.1      0
#>   2:     1     1     2 -0.5646982  0.1540964     0.2      0
#>   3:     1     1     3  0.3631284  0.6817541     0.3      1
#>   4:     1     1     4  0.6328626 -0.7355924     0.4      0
#>   5:     1     2     1  0.4042683  0.7917824     0.1      1
#>  ---                                                       
#> 236:    20    59     4 -0.3909654  0.4483229     0.4      1
#> 237:    20    60     1  1.3487070 -0.2291228     0.1      0
#> 238:    20    60     2 -0.0227647 -0.4036954     0.2      1
#> 239:    20    60     3  0.2442259 -0.8850577     0.3      0
#> 240:    20    60     4 -0.9423717  0.3473859     0.4      0
input <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                           person_col = "pid",
                           alt_covariate_cols = "quality")
str(input$Z)
#>  num [1:4, 1:2] 1 1 1 1 0.1 0.2 0.3 0.4
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "(Intercept)" "quality"
input$alt_mapping
#>    alt_int   alt N_OBS N_CHOICES TAKE_RATE MKT_SHARE
#>      <int> <int> <int>     <int>     <num>     <num>
#> 1:       0    NA    60        14 0.2333333 0.2333333
#> 2:       1     1    60        15 0.2500000 0.2500000
#> 3:       2     2    60         8 0.1333333 0.1333333
#> 4:       3     3    60        10 0.1666667 0.1666667
#> 5:       4     4    60        13 0.2166667 0.2166667
```
