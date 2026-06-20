# Simulate multinomial probit data

Generates synthetic choice data from the MNP data-generating process
estimated by
[`run_mnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md):
latent utility differences against the base alternative (alternative 1),
\$\$w_i = X_i \beta + \delta + \varepsilon_i, \qquad \varepsilon_i \sim
N\_{J-1}(0, \Sigma),\$\$ with alternative \\j \> 1\\ chosen iff
\\w\_{ij} \> \max(0, \max\_{k \neq j} w\_{ik})\\ and the base chosen iff
all \\w\_{ij} \< 0\\. Covariates are Uniform(-1, 1). Choice sets are
balanced (every individual faces all `J` alternatives), as the MNP
estimator requires; there is no outside-option flag — model an outside
good as a zero-covariate base alternative instead.

## Usage

``` r
simulate_mnp_data(
  N = 5000,
  J = 3,
  beta = c(0.8, -0.6),
  delta = NULL,
  Sigma = matrix(c(1, 0.5, 0.5, 1.5), nrow = 2),
  seed = 123
)
```

## Arguments

- N:

  Number of choice situations.

- J:

  Number of alternatives (alternative 1 is the base).

- beta:

  Fixed coefficients for `x1..x{K_x}` (length `K_x = length(beta)`).

- delta:

  ASCs of the differenced utilities, one per non-base alternative
  (length `J - 1`). Defaults to an alternating pattern of
  `c(0.5, -0.5)`.

- Sigma:

  Covariance matrix of the differenced errors (`(J-1) x (J-1)`).

- seed:

  Random seed (`NULL` skips
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

## Value

A `choicer_sim` object. `true_params` contains `beta`, `delta`, and
`Sigma` on the identified scale (see Details).

## Details

The MNP likelihood only identifies parameters up to scale, so
`true_params` is reported on the *identified* scale (normalized by
\\\sigma\_{11}\\): `beta` \\= \beta / \sqrt{\sigma\_{11}}\\, `delta` \\=
\delta / \sqrt{\sigma\_{11}}\\, and `Sigma` \\= \Sigma / \sigma\_{11}\\
— the scale on which
[`run_mnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md)
reports its posterior. With the default `Sigma` (\\\sigma\_{11} = 1\\)
the DGP scale and the identified scale coincide.

## Examples

``` r
# \donttest{
sim <- simulate_mnp_data(N = 1000, J = 3, seed = 123)
print(sim)
#> <choicer_sim: mnp>
#>   settings:
#>     N = 1000
#>     J = 3
#>     K_x = 2
#>     base_alt = 1
#>   rows in $data: 3000
#>   true_params: beta, delta, Sigma
# }
```
