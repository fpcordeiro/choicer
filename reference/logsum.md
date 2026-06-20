# Expected logsum (inclusive value) per choice situation

Computes the expected maximum utility ("logsum" or inclusive value) for
each choice situation, up to the additive constant of the extreme-value
error:

- MNL: \\\log \sum_j \exp(V\_{ij})\\.

- MXL: \\E\_\beta\[\log \sum_j \exp(V\_{ij}(\beta))\]\\, simulated by
  averaging the log-sum-exp *across* the deterministic Halton draws
  (regenerated from `object$draws_info`). Taking the log-sum-exp of
  draw-averaged utilities would understate the expectation (Jensen's
  inequality), so a dedicated kernel
  ([`mxl_logsum`](https://fpcordeiro.github.io/choicer/reference/mxl_logsum.md))
  is used.

- NL: \\\log \sum_b \exp(\lambda_b I\_{ib})\\ with the nest inclusive
  value \\I\_{ib} = \log \sum\_{j \in b} \exp(V\_{ij}/\lambda_b)\\
  (singleton nests have \\\lambda_b = 1\\).

When the model includes an outside option, its normalized utility \\V =
0\\ contributes an \\\exp(0)\\ term to the sum.

## Usage

``` r
logsum(object, newdata = NULL, ...)

# S3 method for class 'choicer_mnl'
logsum(object, newdata = NULL, ...)

# S3 method for class 'choicer_mxl'
logsum(object, newdata = NULL, ...)

# S3 method for class 'choicer_nl'
logsum(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted model object (`choicer_mnl`, `choicer_mxl`, or `choicer_nl`).

- newdata:

  Optional counterfactual data: a data.frame in the fit-time long format
  or a list with `X`, `alt_idx`, `M` (plus `W` for MXL), as in
  [`predict()`](https://rdrr.io/r/stats/predict.html). When `NULL`
  (default), the data stored at fit time is used (requires
  `keep_data = TRUE`).

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector with one logsum per choice situation. With a data.frame
`newdata`, choice situations are ordered by id (as in
[`predict()`](https://rdrr.io/r/stats/predict.html)).

## Details

Logsum *levels* depend on the ASC normalization (and, more generally, on
any additive utility normalization), so only logsum *differences*
between scenarios (e.g. via `newdata`) are meaningful.

## See also

[`consumer_surplus`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md)

## Examples

``` r
# \donttest{
library(data.table)
sim <- simulate_mnl_data(N = 500, J = 3, beta = c(0.8, -0.6), seed = 1,
                         outside_option = FALSE, vary_choice_set = FALSE)
fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#> Optimization run time 0h:0m:0s
head(logsum(fit))
#> [1] 1.1893076 1.5710443 1.4462352 0.2674369 0.6492900 1.4766210
# }
```
