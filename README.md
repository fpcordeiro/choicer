
# choicer: fast discrete-choice models with a focus on economic applications

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`choicer` provides implementations of discrete-choice models with a focus on economic applications. Computationally intensive likelihoods are written in C++ and exposed for use with generic optimizers. Special care is taken to handle high-dimensional alternative-specific constants efficiently. Currently supports multinomial logit (MNL) and mixed logit (MXL); more models will be added.

*The current version is a work in progress and it is not ready for production*

## Installation

You can install the development version of `choicer` with:

``` r
# Using `remotes`
remotes::install_github("fpcordeiro/choicer")

# Or using `pak`
pak::pkg_install("fpcordeiro/choicer")
```

## Example

View `tests/mnlogit_simulation_test.R` for details
``` r
library(choicer)

# Optimization settings
nloptr_opts <- list(
  "algorithm" = "NLOPT_LD_LBFGS",
  "xtol_rel" = 1.0e-8,
  "maxeval" = 1e+3,
  "print_level" = 0L,
  "check_derivatives" = TRUE,
  "check_derivatives_print" = "none"
)

# Prepare and validate data
input_list <- prepare_mnl_data(
  data = dt,
  id_col = "id",
  alt_col = "alt",
  choice_col = "choice_id",
  covariate_cols = c("x1", "x2"),
  outside_opt_label = 0L,
  include_outside_option = FALSE
)

# Initial parameter vector theta_init
theta_init <- runif(J_global + K_x, -1, 1)

# Run the optimization
result <- nloptr::nloptr(
  x0 = theta_init,
  eval_f = mnl_loglik_gradient_parallel,
  opts = nloptr_opts,
  X = input_list$X,
  alt_idx = input_list$alt_idx,
  choice_idx = input_list$choice_idx,
  weights = input_list$weights,
  M = input_list$M,
  use_asc = TRUE,
  include_outside_option = input_list$include_outside_option
)
```

## Alternative packages
There are multiple R packages that offer similar functionalities. You should definitely check them out:
- [mlogit](https://CRAN.R-project.org/package=mlogit)
- [logitr](https://CRAN.R-project.org/package=logitr)
- [apollo](https://CRAN.R-project.org/package=apollo)

# Whishlist:

* summary, predict
* Robust standard errors
* Nested Logit and Generalized Extreme-Value models
* BLP contraction mapping
  - Allow to be used within mixed logit
* Goodness of fit stats
