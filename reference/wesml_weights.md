# WESML weights for choice-based (endogenous stratified) samples

Computes Manski-Lerman (1977) Weighted Exogenous Sample Maximum
Likelihood (WESML) weights for a choice-based sample. The weight for a
choice situation whose chosen alternative is \\j\\ is \\w = Q(j) /
H(j)\\, where \\Q(j)\\ is the population share of alternative \\j\\ and
\\H(j)\\ its sample share among choosers. Using these weights in
[`run_mxlogit`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)
restores consistency under choice-based sampling; pair them with
`se_method = "sandwich"` for valid (robust) standard errors (the plain
inverse-Hessian is invalid under weighting).

## Usage

``` r
wesml_weights(
  data,
  id_col,
  alt_col,
  choice_col,
  Q,
  H = NULL,
  normalize = TRUE,
  attach = FALSE,
  weight_name = ".wesml_weight",
  outside_opt_label = NULL,
  include_outside_option = FALSE
)
```

## Arguments

- data:

  A long-format choice data set (data.frame or data.table), one row per
  alternative per choice situation.

- id_col, alt_col, choice_col:

  Column names identifying the choice situation, the alternative, and
  the 0/1 chosen indicator.

- Q:

  Named numeric vector of population shares, one entry per chosen
  stratum (names matched to `as.character(alt)`), each strictly
  positive. Renormalized to sum 1 if needed.

- H:

  Optional named numeric vector of sample shares. If `NULL` (default) it
  is computed from `data` as the fraction of choice situations choosing
  each alternative.

- normalize:

  If `TRUE` (default) the returned weights are scaled to mean 1. This
  does not affect the point estimates or the sandwich variance.

- attach:

  If `TRUE`, return `data` with a row-level weight column appended (the
  per-situation weight repeated across all rows of a situation), ready
  to pass to `run_mxlogit(weights_col = ...)`. If `FALSE` (default)
  return an id-keyed table of weights.

- weight_name:

  Name of the weight column (default `".wesml_weight"`).

- outside_opt_label, include_outside_option:

  Set `include_outside_option = TRUE` and supply `outside_opt_label`
  when the outside good is implicit (choice situations with no `1` in
  `choice_col` are treated as having chosen the outside good).

## Value

Either an id-keyed `data.table` with columns `id_col` and `weight_name`
(default), or, when `attach = TRUE`, a copy of `data` with the weight
column appended. The result carries `"Q"`, `"H"`, and
`"choice_sampling"` attributes recording provenance.

## Details

Strata are defined by the chosen alternative and keyed by
`as.character(alt)` so numeric and character alternative codes match
supplied share names unambiguously.

## References

Manski, C. F. and Lerman, S. R. (1977). The Estimation of Choice
Probabilities from Choice Based Samples. *Econometrica* 45(8),
1977-1988. Train, K. E. (2009). *Discrete Choice Methods with
Simulation*, Section 3.7. Cambridge University Press.

## See also

[`sample_by_choice`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md),
[`run_mxlogit`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md),
[`wesml_vcov`](https://fpcordeiro.github.io/choicer/reference/wesml_vcov.md)

## Examples

``` r
library(data.table)
set.seed(1)
N <- 300L; J <- 3L
pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
pop[, x1 := rnorm(.N)]
#>         id   alt         x1
#>      <int> <int>      <num>
#>   1:     1     1 -0.6264538
#>   2:     1     2  0.1836433
#>   3:     1     3 -0.8356286
#>   4:     2     1  1.5952808
#>   5:     2     2  0.3295078
#>  ---                       
#> 896:   299     2 -0.0132557
#> 897:   299     3 -0.9314803
#> 898:   300     1  1.2146891
#> 899:   300     2 -2.0883404
#> 900:   300     3 -0.5261525
pop[, w1 := rnorm(.N)]
#>         id   alt         x1         w1
#>      <int> <int>      <num>      <num>
#>   1:     1     1 -0.6264538 -1.5414026
#>   2:     1     2  0.1836433  0.1943211
#>   3:     1     3 -0.8356286  0.2644225
#>   4:     2     1  1.5952808 -1.1187352
#>   5:     2     2  0.3295078  0.6509530
#>  ---                                  
#> 896:   299     2 -0.0132557  1.9363973
#> 897:   299     3 -0.9314803 -1.4558838
#> 898:   300     1  1.2146891  1.4819057
#> 899:   300     2 -2.0883404  1.0761196
#> 900:   300     3 -0.5261525 -0.7574884
pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#>         id   alt         x1         w1 choice
#>      <int> <int>      <num>      <num>  <int>
#>   1:     1     1 -0.6264538 -1.5414026      0
#>   2:     1     2  0.1836433  0.1943211      1
#>   3:     1     3 -0.8356286  0.2644225      0
#>   4:     2     1  1.5952808 -1.1187352      0
#>   5:     2     2  0.3295078  0.6509530      1
#>  ---                                         
#> 896:   299     2 -0.0132557  1.9363973      0
#> 897:   299     3 -0.9314803 -1.4558838      0
#> 898:   300     1  1.2146891  1.4819057      0
#> 899:   300     2 -2.0883404  1.0761196      1
#> 900:   300     3 -0.5261525 -0.7574884      0

# Population shares of the chosen alternative
Q <- prop.table(table(pop[choice == 1, alt]))
wt <- wesml_weights(pop, "id", "alt", "choice", Q = Q)
head(wt)
#> Key: <id>
#>       id .wesml_weight
#>    <int>         <num>
#> 1:     1             1
#> 2:     2             1
#> 3:     3             1
#> 4:     4             1
#> 5:     5             1
#> 6:     6             1
```
