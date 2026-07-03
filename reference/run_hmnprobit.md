# Fit a hierarchical Bayesian multinomial probit (HMNP)

Runs the fully conjugate Albert-Chib Gibbs sampler for the hierarchical
multinomial probit with iid normal utility shocks in un-differenced
utility space: \$\$U\_{ijt} = x\_{ijt}'\beta_i + \delta_j +
\epsilon\_{ijt}, \qquad U\_{iot} = \epsilon\_{iot}, \qquad \epsilon \sim
N(0, \sigma^2),\$\$ choice by argmax within the task including the
stochastic implicit outside option, \\\beta_i \sim N(b, W)\\ (normal
coordinates only — log-normal would break conjugacy), and \\\delta_j =
z_j'\theta + \xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\.

## Usage

``` r
run_hmnprobit(
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

  A `choicer_data_hmnp` object from
  [`prepare_hmnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnp_data.md).

- include_outside_option:

  Logical; if `TRUE` (default) an implicit outside option with
  systematic utility 0 is part of every choice set.

- prior:

  As in
  [`run_hmnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md),
  plus `a0` (3) and `s0` (3), the inverse-gamma shape/scale on the
  non-identified \\\sigma^2\\.

- mcmc:

  Named list overriding MCMC defaults: `R` (10000), `burn` (R %/% 5),
  `thin` (1), `seed` (drawn via
  [`sample.int()`](https://rdrr.io/r/base/sample.html) so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) governs), `trace`
  (0). No proposal-scale settings: the sampler is fully conjugate, with
  no Metropolis steps.

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

A `choicer_hmnp` object (classed `c("choicer_hmnp", "choicer_hb")`); the
same layout as
[`run_hmnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md)'s
return, with all reported summaries on the identified scale, raw chains
in `draws$*_raw`, and the non-identified `draws$sigma2` trace.

## Details

**Identification.** The probit likelihood is invariant to a common
rescaling of utilities and \\\sigma\\, so the chain runs on the
non-identified parameterization (free \\\sigma^2\\, better mixing via
parameter expansion) and every kept draw is normalized by the matching
power of the CURRENT \\\sigma\\: reported \\b/\sigma\\, \\W/\sigma^2\\,
\\\delta/\sigma\\, \\\theta/\sigma\\, \\\sigma_d^2/\sigma^2\\. Raw
chains are kept in `draws$*_raw`. The outside option anchors the
location of \\\delta\\ exactly as in
[`run_hmnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md).

## See also

[`prepare_hmnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_hmnp_data.md),
[`simulate_hmnp_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_hmnp_data.md),
[`run_hmnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md)

## Examples

``` r
# \donttest{
sim <- simulate_hmnp_data(N = 100, T = 3, J = 4, seed = 42)
fit <- run_hmnprobit(sim$data, "task", "alt", "choice", c("x1", "x2"),
                     person_col = "pid", alt_covariate_cols = "z1",
                     mcmc = list(R = 500, burn = 200))
#> MCMC run time 0h:0m:0.02s
#> Warning: Few alternatives relative to the delta mean function (J = 4, P = 2): theta and sigma_d^2 lean on the prior.
summary(fit)
#> Hierarchical Bayesian Multinomial Probit (HMNP) model
#> 
#> Population coefficients b (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> x1           0.680063   0.134663   0.453139   0.673162   0.969598
#> x2          -0.559368   0.140942  -0.824908  -0.551598  -0.280261
#> 
#> Delta mean function theta (posterior):
#> Parameter          Mean         SD       2.5%     Median      97.5%
#> (Intercept)    0.574660   0.403289  -0.118313   0.607332   1.059542
#> z1            -0.369009   0.476782  -1.099452  -0.431757   0.515994
#> 
#> Alternative-effect variance (posterior):
#> Parameter        Mean         SD       2.5%     Median      97.5%
#> sigma_d^2    0.399914   2.645053   0.000433   0.081124   1.375950
#> 
#> Raw shock variance (non-identified chain):
#> Parameter            Mean         SD       2.5%     Median      97.5%
#> sigma^2 (raw)    2.771550   1.112982   1.140880   2.933388   4.968923
#> 
#> Quality ladder (delta = mean utility vs the outside option; xi = delta - z'theta):
#>  alternative delta_mean delta_sd xi_mean  xi_sd
#>            1     0.1041   0.1806 -0.1645 0.4163
#>            2     0.4782   0.1288  0.2261 0.4299
#>            3     0.7919   0.1613  0.0594 0.4974
#>            4     0.2462   0.1660 -0.0846 0.3790
#> 
#> split-R-hat:
#>          x1          x2 (Intercept)          z1    sigma_d2 
#>       1.014       0.997       1.003       1.014       1.001 
#> 
#> Respondents: 100  Choice situations: 300  Alternatives: 4 
#> Draws kept: 300  Chains: 1 
#> MCMC run time 0h:0m:0.02s 
# }
```
