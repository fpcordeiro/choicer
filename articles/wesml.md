# Choice-based sampling and WESML weights

Choice data are often sampled by outcome. A transport researcher running
an on-site survey interviews travellers at the terminal of the mode they
actually chose; a hospital-choice study may oversample patients of rare
hospitals; a marketing team may recruit equal numbers of buyers of each
brand. In each case the unit is drawn *conditional on the alternative it
chose*, so the sample choice shares are not the population choice
shares. Treating such a sample as random changes the likelihood target
and, in general, biases the estimates.

WESML fixes that sampling problem; it does not fix every econometric
problem. The weighted likelihood still relies on the maintained utility
specification and on whatever exogeneity assumptions justify
interpreting the covariates, especially prices, as demand shifters
rather than equilibrium outcomes.

Manski and Lerman’s (1977) weighted exogenous sample maximum likelihood
(WESML) correction weights each choice situation by

``` math
w_i = \frac{Q_{j(i)}}{H_{j(i)}},
```

where $`j(i)`$ is the alternative chosen by situation $`i`$, $`Q_j`$ is
the population share choosing alternative $`j`$, and $`H_j`$ is the
corresponding sample share. Maximizing the weighted log-likelihood
$`\sum_i w_i \log P_i`$ recovers the population parameters. choicer
provides two helpers:

- [`sample_by_choice()`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
  draws a choice-based sample from a population frame and attaches WESML
  weights.
- [`wesml_weights()`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)
  computes the same weights when you already have a sample and know the
  population shares `Q`.

Both helpers normalize the weights to mean 1 by default. Normalization —
and indeed any rescaling of the weights by a common factor — leaves the
point estimates and the robust (sandwich) variance unchanged, so the
attached `.wesml_weight` need not equal $`Q/H`$ literally; only the
*relative* weights across strata matter.

``` r

library(choicer)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
set_num_threads(2)
```

## Build a population

For exposition, start from a simulated population in which tastes are
heterogeneous (a random coefficient on `w1` and `w2`), so a mixed logit
is the natural estimator. We turn off the outside option and fix the
choice set so that every situation has exactly one chosen alternative
and the strata are clean. In empirical work the population shares `Q`
usually come from administrative totals, market shares, or survey
weights external to the choice-based estimation sample.

``` r

sim <- simulate_mxl_data(
  N               = 3000,
  J               = 4,
  Sigma           = diag(c(1.0, 1.5)),  # two uncorrelated random coefficients
  seed            = 11,
  outside_option  = FALSE,
  vary_choice_set = FALSE
)

pop <- as.data.table(sim$data)
Q <- prop.table(table(pop[choice == 1, alt]))
round(Q, 3)
#> 
#>     1     2     3     4 
#> 0.355 0.150 0.330 0.165
```

## Draw a choice-based sample

Now sample the same number of choice situations from each chosen
alternative. This keeps whole choice situations together: if an id is
sampled, all of its alternative rows are retained.

``` r

cb <- sample_by_choice(
  pop,
  id_col     = "id",
  alt_col    = "alt",
  choice_col = "choice",
  n_per_alt  = 300L,
  seed       = 12L
)

strata <- sort(names(attr(cb, "Q")))
rbind(
  population = attr(cb, "Q")[strata],
  sample     = attr(cb, "H")[strata]
) |> round(3)
#>                1    2    3     4
#> population 0.355 0.15 0.33 0.165
#> sample     0.250 0.25 0.25 0.250

cb[choice == 1, .(id, chosen_alt = alt, .wesml_weight)][1:8]
#>       id chosen_alt .wesml_weight
#>    <int>      <int>         <num>
#> 1:     1          2        0.6000
#> 2:     4          1        1.4200
#> 3:     7          1        1.4200
#> 4:     8          1        1.4200
#> 5:    10          3        1.3187
#> 6:    11          2        0.6000
#> 7:    12          4        0.6613
#> 8:    13          3        1.3187
```

The sample choice shares are deliberately equalized, but the attached
weights restore the population shares in the weighted likelihood. The
weight is constant within an id and repeated across that id’s
alternative rows, which is exactly the row-level layout
[`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)
expects through `weights_col`.

## Weighted estimation and inference

We fit two mixed logits on the choice-based sample: an ordinary
(unweighted) fit that ignores the sampling design, and a WESML fit that
passes the weight column and requests the robust sandwich covariance.
Passing `weights_col` by name keeps the estimation target visible in the
script, which is the recommended style even when the data already carry
a `choice_sampling` attribute from
[`sample_by_choice()`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md).

``` r

common <- list(
  data            = cb,
  id_col          = "id",
  alt_col         = "alt",
  choice_col      = "choice",
  covariate_cols  = c("x1", "x2"),  # fixed coefficients
  random_var_cols = c("w1", "w2"),  # random coefficients
  S               = 100L,
  draws           = "generate",
  seed            = 7L,
  scale_vars      = "sd"
)

fit_unweighted <- do.call(run_mxlogit, c(common, list(se_method = "bhhh")))
#> Detected WESML choice-based-sampling provenance; applying attached weights from column '.wesml_weight'.
#> Optimization run time 0h:0m:0.67s
#> Warning: Non-uniform weights detected with se_method = 'bhhh': BHHH/OPG
#> standard errors use the w^1 meat (sum w_i s_i s_i')^{-1}, which is NOT a valid
#> choice-based-sampling (WESML) correction; the correct sandwich meat is w^2. Use
#> se_method = 'sandwich' for valid WESML inference.

fit_wesml <- do.call(run_mxlogit, c(common, list(
  weights_col = ".wesml_weight",
  se_method   = "sandwich"
)))
#> Optimization run time 0h:0m:0.65s
```

> **Tip.** As in the [mixed logit
> vignette](https://fpcordeiro.github.io/choicer/articles/mxl.md), raise
> the number of draws `S` until the estimates are stable and warm-start
> a stubborn solver with `theta_init`. `S = 100` here keeps the package
> build quick.

The unweighted estimator treats the equalized sample shares as if they
were the population shares; WESML reweights the sampled situations back
to the population. With alternative-specific constants in the model the
correction is most visible in the constants and, through them, in the
fitted shares:

``` r

round(cbind(
  unweighted = coef(fit_unweighted),
  wesml      = coef(fit_wesml)
), 3)
#>       unweighted  wesml
#> x1         0.810  0.810
#> x2        -0.561 -0.561
#> L_11       0.169  0.169
#> L_22       0.127  0.127
#> ASC_2     -1.107 -1.107
#> ASC_3     -0.080 -0.080
#> ASC_4     -0.978 -0.978
```

``` r

share_compare <- rbind(
  population = as.numeric(Q),
  wesml      = drop(predict(fit_wesml, type = "shares")),
  unweighted = drop(predict(fit_unweighted, type = "shares"))
)
colnames(share_compare) <- names(Q)
round(share_compare, 3)
#>                1     2    3     4
#> population 0.355 0.150 0.33 0.165
#> wesml      0.358 0.149 0.33 0.163
#> unweighted 0.358 0.149 0.33 0.163
```

The WESML fit reproduces the population shares `Q`, while the unweighted
fit reproduces the equalized *sample* shares — a direct picture of the
bias the correction removes. In a single finite sample the WESML
estimates need not be closer to the truth parameter by parameter, but
they target the population likelihood under the choice-based sampling
design.

For inference, the point of `se_method = "sandwich"` is that under
non-uniform weights the inverse weighted Hessian and the ordinary BHHH
variance are *not* valid covariance estimators. The sandwich uses the
weighted Hessian as bread, $`A = \sum_i w_i(-H_i)`$, and the
weight-squared outer product of the per-situation scores as meat,
$`B = \sum_i w_i^2 s_i s_i'`$, giving $`V = A^{-1} B A^{-1}`$. Because
$`A`$ scales linearly and $`B`$ quadratically in the weights, $`V`$ is
invariant to any common rescaling of them — consistent with the mean-1
normalization above.

``` r

summary(fit_wesml)
#> Mixed Logit (MXL) model
#> 
#> Parameter    Estimate  Std.Error  z-value  Pr(>|z|)  
#> x1           0.809770   0.080489  10.0607  0.00e+00  ***
#> x2          -0.561483   0.075121  -7.4743  7.75e-14  ***
#> Sigma_11     1.402991   0.539452   2.6008  9.30e-03  **
#> Sigma_22     1.289004   0.533570   2.4158  1.57e-02  *
#> ASC_2       -1.107296   0.107727 -10.2787  0.00e+00  ***
#> ASC_3       -0.079835   0.101007  -0.7904  4.29e-01  
#> ASC_4       -0.978022   0.107984  -9.0571  0.00e+00  ***
#> ---
#> Signif. codes:  '***' 0.001 '**' 0.01 '*' 0.05
#> 
#> Random coefficient covariance (Sigma):
#>       w1    w2
#> w1 1.403 0.000
#> w2 0.000 1.289
#> 
#> Std. Errors: Sandwich (robust) 
#> Weighting: WESML choice-based 
#> Log-likelihood: -1471.18 
#> AIC: 2956.36  | BIC: 2991.99 
#> McFadden R2: 0.116 (adj: 0.111) | Hit rate: 0.436 
#> N: 1200  | Parameters: 7 
#> Optimization time: 0.65 s
#> Convergence: 3 ( NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. )
```

The same robust variance is available post hoc on any fitted mixed logit
via
[`wesml_vcov()`](https://fpcordeiro.github.io/choicer/reference/wesml_vcov.md),
so you can obtain choice-based-sampling standard errors even from a fit
estimated with `se_method = "hessian"` without refitting.

> **A note on the multinomial logit.** choicer implements WESML
> weighting and the robust sandwich for the *mixed* logit (which nests
> the plain logit as the degenerate, zero-variance case). For the plain
> multinomial logit there is a classical and convenient result (Manski
> and Lerman, 1977): when the model includes a full set of
> alternative-specific constants, choice-based sampling leaves the slope
> coefficients consistently estimated *even without weighting* — only
> the ASCs are inconsistent. Each constant is shifted by
> $`\ln\!\big(H_j / Q_j\big)`$ and can be corrected by subtracting that
> term. So for an MNL with ASCs the substantive marginal-utility
> parameters are unaffected by the sampling scheme; only the constants
> (and the predicted shares they drive) need correcting.

## Starting from an existing sample

When the choice-based sample already exists, provide the population
shares `Q` directly:

``` r

cb2 <- copy(cb)
cb2[, .wesml_weight := NULL]

cb2 <- wesml_weights(
  cb2,
  id_col     = "id",
  alt_col    = "alt",
  choice_col = "choice",
  Q          = attr(cb, "Q"),
  attach     = TRUE
)

attr(cb2, "choice_sampling")
#> $scheme
#> [1] "wesml"
#> 
#> $Q
#>      2      3      1      4 
#> 0.1500 0.3297 0.3550 0.1653 
#> 
#> $H
#>    1    2    3    4 
#> 0.25 0.25 0.25 0.25 
#> 
#> $meat
#> [1] "robust"
#> 
#> $source
#> [1] "wesml_weights"
#> 
#> $weight_name
#> [1] ".wesml_weight"
```

The names of `Q` must match the chosen-alternative strata exactly after
coercion to character. This strict matching is intentional: silently
dropping a realized stratum would change the target population.

## References

Manski, C. F. and Lerman, S. R. (1977). The estimation of choice
probabilities from choice based samples. *Econometrica*, 45(8),
1977-1988.

Train, K. E. (2009). *Discrete Choice Methods with Simulation* (2nd
ed.). Cambridge University Press, Section 3.7.
