# Runs Bayesian multinomial probit estimation

Estimates a multinomial probit model by Gibbs sampling with data
augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). The model is
specified in utility differences against a base alternative: for choice
situation \\i\\ with \\J\\ alternatives, \\w_i = X_i \beta +
\epsilon_i\\ with \\\epsilon_i \sim N\_{J-1}(0, \Sigma)\\.

## Usage

``` r
run_mnprobit(
  data = NULL,
  id_col = NULL,
  alt_col = NULL,
  choice_col = NULL,
  covariate_cols = NULL,
  input_data = NULL,
  base_alt = NULL,
  use_asc = TRUE,
  prior = list(),
  mcmc = list(),
  keep_data = TRUE
)
```

## Arguments

- data:

  Data frame containing choice data (convenience workflow). Mutually
  exclusive with `input_data`.

- id_col:

  Name of the column identifying choice situations (individuals).

- alt_col:

  Name of the column identifying alternatives.

- choice_col:

  Name of the column indicating chosen alternative (1 = chosen, 0 = not
  chosen).

- covariate_cols:

  Vector of names of columns to be used as covariates.

- input_data:

  List output from
  [`prepare_mnp_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnp_data.md)
  (advanced workflow). Mutually exclusive with `data`.

- base_alt:

  Label of the base (reference) alternative used for utility
  differencing. If `NULL` (default), the first alternative in sort order
  is used.

- use_asc:

  Logical indicating whether to include alternative-specific constants
  (one intercept per non-base alternative in the differenced utilities).

- prior:

  Named list of prior settings, merged over defaults:

  `beta_bar`

  :   Prior mean of \\\beta\\ (default `rep(0, K)`).

  `A`

  :   Prior precision of \\\beta\\ (default `0.01 * diag(K)`).

  `nu`

  :   Inverse-Wishart degrees of freedom (default `p + 3`).

  `V`

  :   Inverse-Wishart scale matrix (default `nu * diag(p)`).

- mcmc:

  Named list of MCMC settings, merged over defaults:

  `R`

  :   Total Gibbs iterations (default `10000`).

  `burn`

  :   Burn-in iterations discarded (default `floor(R / 5)`).

  `thin`

  :   Keep every `thin`-th post-burn-in draw (default `1`).

  `seed`

  :   Master RNG seed (default: drawn from R's RNG).

  `trace`

  :   Print progress every `trace` iterations (default `0`, silent).

- keep_data:

  Logical. If `TRUE` (default), stores prepared data in the returned
  object.

## Value

A `choicer_mnp` object. S3 methods available:
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html) (posterior means of
identified coefficients), [`vcov()`](https://rdrr.io/r/stats/vcov.html)
(posterior covariance of identified coefficient draws),
[`nobs()`](https://rdrr.io/r/stats/nobs.html). Posterior draws are
stored in `$draws` (`beta` / `sigma` on the identified scale, `beta_raw`
/ `sigma_raw` unnormalized).

## Details

Two workflows are supported:

- Convenience (default):

  Supply `data` and column names. Data preparation
  ([`prepare_mnp_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnp_data.md))
  is handled automatically.

- Advanced:

  Call
  [`prepare_mnp_data`](https://fpcordeiro.github.io/choicer/reference/prepare_mnp_data.md)
  yourself and pass the result via `input_data`.

**Identification.** The multinomial probit likelihood is invariant to a
common rescaling \\(\beta, \Sigma) \to (c\beta, c^2\Sigma)\\. The
sampler runs on the non-identified parameterization (unrestricted
\\\Sigma\\ with an inverse-Wishart prior) and identified quantities are
computed by normalizing each kept draw by \\\sigma\_{11}\\: \\\beta /
\sqrt{\sigma\_{11}}\\ and \\\Sigma / \sigma\_{11}\\. This is the
McCulloch & Rossi (1994) default, which keeps all Gibbs conditionals
conjugate and mixes better than the fully identified sampler of
McCulloch, Polson & Rossi (2000). Reported coefficients, standard
deviations, and credible intervals are posterior summaries of the
identified draws.

**Reproducibility.** The sampler uses its own thread-safe RNG with one
stream per (iteration, observation), so results are reproducible
independent of the number of OpenMP threads (see
[`set_num_threads()`](https://fpcordeiro.github.io/choicer/reference/set_num_threads.md)).
When `mcmc$seed` is not supplied, a master seed is drawn from R's RNG,
so [`set.seed()`](https://rdrr.io/r/base/Random.html) controls the run.

**Scope.** Balanced choice sets are required: every choice situation
must contain the same \\J\\ alternatives. To model an outside option,
include it as explicit rows with zero covariates and set `base_alt` to
its label.

## References

Albert, J. H., & Chib, S. (1993). Bayesian Analysis of Binary and
Polychotomous Response Data. *Journal of the American Statistical
Association*, 88(422), 669-679.

McCulloch, R., & Rossi, P. E. (1994). An exact likelihood analysis of
the multinomial probit model. *Journal of Econometrics*, 64(1-2),
207-240.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(42)
N <- 200; J <- 3; beta_true <- c(1.0, -0.5)
dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#>         id   alt         x1         x2
#>      <int> <int>      <num>      <num>
#>   1:     1     1  1.3709584 -0.2484829
#>   2:     1     2 -0.5646982  0.4223204
#>   3:     1     3  0.3631284  0.9876533
#>   4:     2     1  0.6328626  0.8355682
#>   5:     2     2  0.4042683 -0.6605219
#>  ---                                  
#> 596:   199     2  0.1603274  1.4683500
#> 597:   199     3 -0.4336419 -0.8764555
#> 598:   200     1  1.5374124 -1.2266047
#> 599:   200     2 -2.1702466  0.3378379
#> 600:   200     3  1.0270046  0.4408241
dt[, U := drop(as.matrix(.SD) %*% beta_true) + rnorm(.N), .SDcols = c("x1", "x2")]
#>         id   alt         x1         x2           U
#>      <int> <int>      <num>      <num>       <num>
#>   1:     1     1  1.3709584 -0.2484829  0.74868346
#>   2:     1     2 -0.5646982  0.4223204 -0.73925224
#>   3:     1     3  0.3631284  0.9876533  0.19261139
#>   4:     2     1  0.6328626  0.8355682  0.59475455
#>   5:     2     2  0.4042683 -0.6605219  1.61108575
#>  ---                                              
#> 596:   199     2  0.1603274  1.4683500 -0.54232500
#> 597:   199     3 -0.4336419 -0.8764555 -1.90026895
#> 598:   200     1  1.5374124 -1.2266047  2.02695959
#> 599:   200     2 -2.1702466  0.3378379 -2.33883716
#> 600:   200     3  1.0270046  0.4408241 -0.03718577
dt[, choice := as.integer(U == max(U)), by = id]
#>         id   alt         x1         x2           U choice
#>      <int> <int>      <num>      <num>       <num>  <int>
#>   1:     1     1  1.3709584 -0.2484829  0.74868346      1
#>   2:     1     2 -0.5646982  0.4223204 -0.73925224      0
#>   3:     1     3  0.3631284  0.9876533  0.19261139      0
#>   4:     2     1  0.6328626  0.8355682  0.59475455      0
#>   5:     2     2  0.4042683 -0.6605219  1.61108575      1
#>  ---                                                     
#> 596:   199     2  0.1603274  1.4683500 -0.54232500      0
#> 597:   199     3 -0.4336419 -0.8764555 -1.90026895      0
#> 598:   200     1  1.5374124 -1.2266047  2.02695959      1
#> 599:   200     2 -2.1702466  0.3378379 -2.33883716      0
#> 600:   200     3  1.0270046  0.4408241 -0.03718577      0

fit <- run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
                    mcmc = list(R = 500, burn = 100))
#> MCMC run time 0h:0m:0.01s
summary(fit)
#> Bayesian Multinomial Probit (MNP) model
#> 
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> x1           0.772101   0.123511   0.566925   0.752453   1.031363
#> x2          -0.541413   0.114047  -0.778688  -0.536530  -0.356168
#> ASC_2       -0.265011   0.189260  -0.625681  -0.261778   0.089678
#> ASC_3       -0.224619   0.203268  -0.686996  -0.193081   0.095245
#> 
#> Covariance of utility differences (Sigma, identified scale):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> Sigma_11     1.000000   0.000000   1.000000   1.000000   1.000000
#> Sigma_21    -0.348633   0.446461  -1.302204  -0.293338   0.458029
#> Sigma_22     1.732825   0.809985   0.727891   1.553463   3.737842
#> 
#> Posterior mean Sigma:
#>           w_2       w_3
#> w_2  1.000000 -0.348633
#> w_3 -0.348633  1.732825
#> 
#> Base alternative: 1 
#> Draws kept: 400 (R = 500, burn = 100, thin = 1, seed = 476683043)
#> N: 200  | Parameters: 4 
#> Sampling time: 0.01 s
#> Identification: per-draw normalization by sigma_11 (McCulloch-Rossi 1994).
coef(fit)
#>         x1         x2      ASC_2      ASC_3 
#>  0.7721010 -0.5414133 -0.2650109 -0.2246186 
# }
```
