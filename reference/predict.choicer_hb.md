# Posterior choice probabilities and shares for hierarchical Bayes fits

Computes counterfactual choice probabilities, integrating over the
posterior draws and (at the population level) over the
random-coefficient distribution: for each kept draw \\(b_r, W_r,
\delta_r, \ldots)\\ one \\\beta \sim N(b_r, W_r)\\ is drawn and the
model probabilities are averaged. Alternatives in `newdata` that were
not in the estimation sample receive a posterior-predictive
\\\delta\_{new} \sim N(z\_{new}'\theta_r, \sigma\_{d,r}^2)\\ — the entry
counterfactual unlocked by the random-effects \\\delta\\. Price or
subsidy counterfactuals are just modified covariate columns in
`newdata`.

## Usage

``` r
# S3 method for class 'choicer_hb'
predict(
  object,
  newdata = NULL,
  level = c("population", "individual"),
  n_draws = 200L,
  aggregate = TRUE,
  ...
)
```

## Arguments

- object:

  A `choicer_hmnl` or `choicer_hmnp` fit.

- newdata:

  Data frame with the estimation columns (choice column not required).
  `NULL` (default) predicts on the estimation data.

- level:

  `"population"` (default) integrates over \\N(b, W)\\; `"individual"`
  uses the respondent-level \\\beta_i\\ (requires the prediction rows to
  belong to estimation respondents; fully posterior-integrated when the
  fit kept `keep_beta_i = "draws"`).

- n_draws:

  Number of posterior draws to integrate over (thinned evenly from the
  kept draws; default 200).

- aggregate:

  If `TRUE` (default) return a per-alternative posterior share table
  (including the outside option); if `FALSE` return the posterior-mean
  probability per row of the prediction data.

- ...:

  Ignored.

## Value

With `aggregate = TRUE`, a `data.table` with columns `alternative`,
`share` (posterior mean), `sd`, `lower`, `upper` (95% equal-tailed
interval); the posterior share draws are attached as `attr(, "draws")`.
With `aggregate = FALSE`, a numeric vector of posterior-mean choice
probabilities, one per prediction row.

## Details

HMNL probabilities are closed-form logit; HMNP probabilities use the 1-D
Gauss-Hermite representation of the iid-probit integral \\P(j) = \int
\phi(u) \prod\_{k \ne j} \Phi(V_j - V_k + u) du\\.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 100, T = 3, J = 4, seed = 42)
fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                    person_col = "pid", alt_covariate_cols = "z1",
                    mcmc = list(R = 500, burn = 200))
#> MCMC run time 0h:0m:0.04s
#> Warning: Few alternatives relative to the delta mean function (J = 4, P = 2): theta and sigma_d^2 lean on the prior.
#> Warning: split-R-hat exceeds 1.05 for: z1. Consider a longer run.
predict(fit)                       # posterior shares, estimation data
#>    alternative     share         sd     lower     upper
#>         <char>     <num>      <num>     <num>     <num>
#> 1:           1 0.1861200 0.01664956 0.1560522 0.2199947
#> 2:           2 0.1895338 0.02042708 0.1565832 0.2316484
#> 3:           3 0.2841487 0.02256055 0.2359523 0.3257916
#> 4:           4 0.1968078 0.01553763 0.1642162 0.2280392
#> 5:   (outside) 0.1433896 0.02301090 0.1042867 0.1899193
cf <- sim$data
cf$x1 <- cf$x1 + 0.5               # a counterfactual attribute change
predict(fit, newdata = cf)
#>    alternative     share         sd      lower     upper
#>         <char>     <num>      <num>      <num>     <num>
#> 1:           1 0.1935651 0.02206788 0.15326894 0.2375603
#> 2:           2 0.1947846 0.02428569 0.15318034 0.2462673
#> 3:           3 0.2949920 0.02519159 0.23801711 0.3437924
#> 4:           4 0.2049632 0.01821599 0.17148835 0.2399649
#> 5:   (outside) 0.1116951 0.04918297 0.03281954 0.2079229
# }
```
