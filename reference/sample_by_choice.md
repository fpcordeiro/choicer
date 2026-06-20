# Draw a choice-based sample stratified by the chosen alternative

Subsamples whole choice situations from a population data set according
to fixed per-stratum quotas, where strata are defined by the *chosen*
alternative. The input data set is treated as the population, so the
population shares \\Q(j)\\ are known exactly; the returned sample
carries a ready-to-use WESML weight column (see
[`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)).

## Usage

``` r
sample_by_choice(
  data,
  id_col,
  alt_col,
  choice_col,
  n_per_alt = NULL,
  frac_per_alt = NULL,
  seed = NULL,
  weight_name = ".wesml_weight",
  outside_opt_label = NULL,
  include_outside_option = FALSE
)
```

## Arguments

- data, id_col, alt_col, choice_col:

  As in
  [`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md).

- n_per_alt:

  Either a single integer applied to every stratum, or a named integer
  vector of per-stratum counts (names matched to `as.character(alt)`,
  covering all strata). Mutually exclusive with `frac_per_alt`.

- frac_per_alt:

  Either a single fraction in `[0, 1]` applied to every stratum, or a
  named numeric vector of per-stratum fractions. Mutually exclusive with
  `n_per_alt`.

- seed:

  Optional integer seed for reproducible sampling.

- weight_name:

  Name of the attached weight column (default `".wesml_weight"`).

- outside_opt_label, include_outside_option:

  As in
  [`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)
  (for an implicit outside good).

## Value

A `data.table` subsample with the weight column appended and `"Q"`,
`"H"`, and `"choice_sampling"` attributes (the last records the scheme,
shares, quotas, and `meat = "robust"`).

## Details

Sampling is by choice situation (id), never by row: all alternative-rows
of a sampled situation are kept together. Sampling is **without
replacement**.

## References

Manski, C. F. and Lerman, S. R. (1977). *Econometrica* 45(8), 1977-1988.

## See also

[`wesml_weights`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md),
[`run_mxlogit`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)

## Examples

``` r
library(data.table)
set.seed(1)
N <- 600L; J <- 3L
pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
pop[, x1 := rnorm(.N)]
#>          id   alt         x1
#>       <int> <int>      <num>
#>    1:     1     1 -0.6264538
#>    2:     1     2  0.1836433
#>    3:     1     3 -0.8356286
#>    4:     2     1  1.5952808
#>    5:     2     2  0.3295078
#>   ---                       
#> 1796:   599     2  1.9363973
#> 1797:   599     3 -1.4558838
#> 1798:   600     1  1.4819057
#> 1799:   600     2  1.0761196
#> 1800:   600     3 -0.7574884
pop[, w1 := rnorm(.N)]
#>          id   alt         x1         w1
#>       <int> <int>      <num>      <num>
#>    1:     1     1 -0.6264538  0.7140855
#>    2:     1     2  0.1836433  0.5813846
#>    3:     1     3 -0.8356286 -0.1467239
#>    4:     2     1  1.5952808  1.5069818
#>    5:     2     2  0.3295078 -0.2795326
#>   ---                                  
#> 1796:   599     2  1.9363973 -0.5294766
#> 1797:   599     3 -1.4558838  0.7047356
#> 1798:   600     1  1.4819057 -0.9388858
#> 1799:   600     2  1.0761196  0.8752661
#> 1800:   600     3 -0.7574884 -0.7443670
pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#>          id   alt         x1         w1 choice
#>       <int> <int>      <num>      <num>  <int>
#>    1:     1     1 -0.6264538  0.7140855      0
#>    2:     1     2  0.1836433  0.5813846      0
#>    3:     1     3 -0.8356286 -0.1467239      1
#>    4:     2     1  1.5952808  1.5069818      1
#>    5:     2     2  0.3295078 -0.2795326      0
#>   ---                                         
#> 1796:   599     2  1.9363973 -0.5294766      0
#> 1797:   599     3 -1.4558838  0.7047356      0
#> 1798:   600     1  1.4819057 -0.9388858      1
#> 1799:   600     2  1.0761196  0.8752661      0
#> 1800:   600     3 -0.7574884 -0.7443670      0

s <- sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 50L, seed = 1L)
attr(s, "choice_sampling")$H   # realized sample shares
#>         3         1         2 
#> 0.3333333 0.3333333 0.3333333 
head(s[[".wesml_weight"]])
#> [1] 0.970 0.970 0.970 0.985 0.985 0.985
```
