# Fit a hierarchical Bayesian multinomial logit (HMNL)

Runs the adaptive RW-Metropolis-within-Gibbs sampler for the
hierarchical (random-coefficients, panel or cross-sectional) multinomial
logit with a BLP-style alternative-level random effect: \$\$U\_{ijt} =
x\_{ijt}'\gamma_i + \delta_j + \epsilon\_{ijt}, \qquad U\_{iot} =
\epsilon\_{iot},\$\$ with i.i.d. Gumbel shocks (including on the
implicit outside option, whose systematic utility is 0), \\\beta_i \sim
N(b, W)\\ over the structural covariates (\\\gamma\_{ik} = \beta\_{ik}\\
or \\\exp(\beta\_{ik})\\ per `rc_dist`), and \\\delta_j = z_j'\theta +
\xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\. Partial pooling shrinks each
\\\delta_j\\ toward its characteristics-based mean \\z_j'\theta\\; the
outside option anchors the level of \\\delta\\ (mean utility relative to
the outside good), so no base alternative or sum-to-zero constraint is
needed.

## Usage

``` r
run_hmnlogit(
  data = NULL,
  id_col = NULL,
  alt_col = NULL,
  choice_col = NULL,
  covariate_cols = NULL,
  person_col = NULL,
  alt_covariate_cols = NULL,
  outside_opt_label = NULL,
  cf_residual_col = NULL,
  input_data = NULL,
  include_outside_option = TRUE,
  rc_dist = NULL,
  prior = list(),
  mcmc = list(),
  chains = 1,
  keep_beta_i = c("means", "draws", "none"),
  keep_data = TRUE
)
```

## Arguments

- data:

  Data frame (convenience pathway). Supply either `data` (with the
  column names) or `input_data`, not both.

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

- input_data:

  A `choicer_data_hmnl` object from
  [`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md)
  (advanced pathway).

- include_outside_option:

  Logical; if `TRUE` (default) an implicit outside option with
  systematic utility 0 is part of every choice set.

- rc_dist:

  Integer vector, one entry per column of `covariate_cols`: `0` for a
  normal random coefficient, `1` for log-normal (the coefficient enters
  utility as `exp(beta_ik)`; hierarchy normal on the log scale). Default
  `NULL` is all-normal. Automatically aligned through dropped columns; a
  `cf_residual_col` coordinate is always normal.

- prior:

  Named list overriding prior defaults: `b_bar` (0), `A` (0.01 I), `nu`
  (K + 3), `V` (nu I), `theta_bar` (0), `A_theta` (0.01 I), `sd_prior`
  (list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3)).

- mcmc:

  Named list overriding MCMC defaults: `R` (10000), `burn` (R %/% 5),
  `thin` (1), `seed` (drawn via
  [`sample.int()`](https://rdrr.io/r/base/sample.html) so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) governs), `trace`
  (0), `s_init` (2.38 / sqrt(K)), `accept_target` (0.234).

- chains:

  Number of independent chains (seeds offset by 1). Chain 1 provides the
  reported draws; all chains feed the split-R-hat table.

- keep_beta_i:

  `"means"` (default) stores posterior means/SDs of the individual-level
  \\\beta_i\\; `"draws"` additionally stores the full (K, N, R_keep)
  draw cube (memory-guarded); `"none"` stores neither.

- keep_data:

  Logical; keep the prepared data on the fit (default `TRUE`, needed by
  post-estimation methods).

## Value

A `choicer_hmnl` object (classed `c("choicer_hmnl", "choicer_hb")`) with
posterior summaries (`coefficients`, `se`, `vcov` for \\b\\;
`theta_summary`; `sigma_d2_summary`; `W_mean`; `delta` and `xi`
quality-ladder tables; `beta_i`), the raw thinned `draws`, acceptance
diagnostics in `accept`, the split-R-hat table in `rhat` (when
`chains > 1`), and sampler metadata.

## Details

**Initialization.** \\\beta_i\\ start at the pooled MNL maximum
likelihood estimate over the structural covariates (log-normal
coordinates transformed to the chain scale with a warn-and-clamp at
0.05); \\\delta\\ starts at shrunk log choice-share contrasts against
the outside option; \\\theta\\ at the OLS regression of the initial
\\\delta\\ on `Z`.

**Priors.** \\b \sim N(b\\bar, A^{-1})\\, \\W \sim IW(\nu, V)\\,
\\\theta \sim N(\theta\\bar, A\_\theta^{-1})\\, and \\\sigma_d \sim\\
half-Cauchy\\(0, s_d)\\ via the Makalic-Schmidt scale mixture (set
`sd_prior$half_cauchy = FALSE` for a plain \\IG(c_0, d_0)\\ on
\\\sigma_d^2\\).

**Endogeneity.** If a price-like covariate is endogenous (correlated
with \\\xi_j\\), supply a first-stage residual via `cf_residual_col`
(Petrin & Train 2010); posterior uncertainty does NOT propagate
first-stage estimation error.

## See also

[`prepare_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnl_data.md),
[`simulate_hmnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_hmnl_data.md),
[`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md),
[`rhat()`](https://fpcordeiro.github.io/choicer/reference/rhat.md)

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
summary(fit)
#> Hierarchical Bayesian Multinomial Logit (HMNL) model
#> 
#> Population coefficients b (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> x1           0.753790   0.161606   0.451668   0.734653   1.088809
#> x2          -0.597166   0.213346  -0.980699  -0.598832  -0.192214
#> 
#> Delta mean function theta (posterior):
#> Parameter          Mean         SD       2.5%     Median      97.5%
#> (Intercept)    0.414833   0.184127   0.070274   0.441049   0.711342
#> z1            -0.408093   0.257467  -0.969343  -0.429151   0.130759
#> 
#> Alternative-effect variance (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> sigma_d^2    0.071588   0.297311   0.000085   0.009473   0.599862
#> 
#> Quality ladder (delta = mean utility vs the outside option; xi = delta - z'theta):
#>  alternative delta_mean delta_sd xi_mean  xi_sd
#>            1     0.0568   0.1506 -0.0195 0.1403
#>            2     0.0780   0.1652  0.0199 0.1702
#>            3     0.5801   0.1510 -0.0092 0.2189
#>            4     0.1382   0.1493 -0.0069 0.1319
#> 
#> Mean acceptance: beta 0.24, delta 0.51
#> 
#> split-R-hat:
#>          x1          x2 (Intercept)          z1    sigma_d2 
#>       1.023       1.009       1.003       1.103       1.016 
#> 
#> Respondents: 100  Choice situations: 300  Alternatives: 4 
#> Draws kept: 300  Chains: 1 
#> MCMC run time 0h:0m:0.04s 
coef(fit, component = "delta")
#>          1          2          3          4 
#> 0.05678103 0.07797479 0.58013663 0.13823079 
# }
```
