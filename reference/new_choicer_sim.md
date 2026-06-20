# Construct a `choicer_sim` object

Wraps simulated data, true parameter values, and DGP settings into a
classed list. Returned by
[`simulate_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnl_data.md),
[`simulate_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mxl_data.md),
and
[`simulate_nl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_nl_data.md),
and consumed by
[`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md).

## Usage

``` r
new_choicer_sim(data, true_params, settings, model)
```

## Arguments

- data:

  A `data.table` of simulated choice observations.

- true_params:

  Named list of true DGP parameters (e.g. `beta`, `delta`, `Sigma`,
  `mu`, `lambdas`).

- settings:

  Named list of DGP settings (e.g. `N`, `J`, `K_x`).

- model:

  Character scalar: `"mnl"`, `"mxl"`, `"nl"`, or `"mnp"`.

## Value

A list of class `choicer_sim`.
