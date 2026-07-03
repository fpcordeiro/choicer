# The math behind choicer: hierarchical Bayesian multinomial probit

This document describes the hierarchical Bayesian multinomial probit
with i.i.d. shocks (HMNP) as implemented in `src/hmnprobit.cpp`, with
shared infrastructure in `src/hb_internal.h`. The R-side entry point is
[`run_hmnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_hmnprobit.md).
The hierarchical structure ($`\beta_i \sim N(b, W)`$;
$`\delta_j = z_j'\theta + \xi_j`$) is identical to the HMNL (see the
companion article); this note focuses on what differs: the shock family,
the data augmentation, identification, and the parallelization
asymmetry.

## Table of Contents

1.  [Notation](#notation)
2.  [Model Definition](#id_1-model-definition)
3.  [Identification](#id_2-identification)
4.  [Priors and Joint Posterior](#id_3-priors-and-joint-posterior)
5.  [The Gibbs Sampler](#id_4-the-gibbs-sampler)
6.  [Sampling Primitives](#id_5-sampling-primitives)
7.  [Post-Processing and Reported
    Quantities](#id_6-post-processing-and-reported-quantities)
8.  [Implementation Details](#id_7-implementation-details)
9.  [References](#references)

------------------------------------------------------------------------

## Notation

As in the HMNL article, plus: $`U_{ijt}`$ the latent utilities
(augmented data), $`U_{iot}`$ the per-task outside latent, $`\sigma^2`$
the (non-identified) common shock variance, $`G_i = X_i'X_i`$ the
respondent Gram block.

## 1. Model Definition

``` math
U_{ijt} = x_{ijt}'\beta_i + \delta_j + \varepsilon_{ijt}, \qquad U_{iot} = \varepsilon_{iot}, \qquad \varepsilon \sim \text{i.i.d. } N(0, \sigma^2),
```

with choice by argmax within the task **including the outside latent**.
The outside option is *stochastic* — its own shock on top of systematic
utility 0, not a deterministic bound. This keeps the model structurally
parallel to the HMNL (whose outside good necessarily carries an EV1
shock): utilities differenced against the outside are equicorrelated
($`\operatorname{Var} = 2\sigma^2`$, $`\operatorname{Cov} = \sigma^2`$),
the exact probit analog of logit-with-outside-good, and it is what makes
the $`\pi/\sqrt6`$ scale comparison between the two models meaningful.

The sampler works in **un-differenced utility space**: with i.i.d.
shocks the latent conditionals are univariate normals whose only
coupling is the argmax truncation, unbalanced choice sets are trivial,
and there is no error covariance matrix to sample — one scalar
$`\sigma^2`$. Respondent coefficients are normal only
($`\beta_i \sim N(b, W)`$; a log-normal coordinate would break
conjugacy), and the alternative level is exactly the HMNL’s:
$`\delta_j = z_j'\theta + \xi_j`$, $`\xi_j \sim N(0, \sigma_d^2)`$.

## 2. Identification

- **Scale.** The likelihood is invariant to a common rescaling of
  $`(U, \sigma)`$. The chain runs on the *non-identified*
  parameterization (free $`\sigma^2`$ — parameter expansion, which mixes
  better than pinning $`\sigma^2 = 1`$), and every kept draw is
  normalized by the matching power of the **current** $`\sigma`$:
  reported $`b/\sigma`$, $`W/\sigma^2`$, $`\delta/\sigma`$,
  $`\theta/\sigma`$, $`\sigma_d^2/\sigma^2`$, $`\beta_i/\sigma`$. This
  is the McCulloch–Rossi $`\sigma_{11}`$ discipline: identified
  summaries are computed on per-draw-normalized chains, never by
  normalizing posterior means ($`E[b/\sigma] \ne E[b]/E[\sigma]`$). Raw
  chains are kept in `draws$*_raw`.
- **Location.** As in the HMNL, the outside option anchors the level of
  $`\delta`$; version 1 requires it.
- **Degenerate cross-check.** With one inside alternative ($`J = 1`$)
  the model is a binary probit whose composite error
  $`\varepsilon_1 - \varepsilon_o`$ has variance $`2\sigma^2`$, so the
  identified coefficients equal $`\sqrt2`$ times the
  `glm(..., binomial(link = "probit"))` coefficients — a closed-form
  validation used in the test suite.

## 3. Priors and Joint Posterior

As in the HMNL for $`(b, W, \theta, \sigma_d^2)`$, plus
$`\sigma^2 \sim \mathcal{IG}(a_0, s_0)`$ (defaults $`a_0 = s_0 = 3`$) on
the non-identified variance. All conditionals below are conjugate — the
sampler has no Metropolis steps.

## 4. The Gibbs Sampler

One systematic-scan iteration:

**(a) Fused respondent pass** (work-shared over respondents): for each
task, a truncated-normal sweep over the inside rows *and the outside
latent* in fixed order — the chosen latent is drawn on
$`(\max \text{others}, \infty)`$, every other latent on
$`(-\infty, U_{\text{chosen}})`$, each conditional a plain
$`N(\text{mean}, \sigma^2)`$ with mean $`x'\beta_i + \delta_j`$ (or 0
for the outside latent, which participates exactly like a row and enters
the bounds at its current value). Then the conjugate $`\beta_i`$
regression draw on the $`\delta`$-residualized utilities,

``` math
\beta_i \mid \cdot \sim N\big(Q_i^{-1} r_i,\; Q_i^{-1}\big), \qquad Q_i = W^{-1} + G_i/\sigma^2, \qquad r_i = W^{-1} b + X_i'(U_i - \delta)/\sigma^2,
```

with $`G_i = X_i'X_i`$ precomputed once; then a refresh of the row cache
$`\mu_x(r) = x_r'\beta_i`$.

**(b) $`\delta_j`$ — work-shared conjugate draws.** Given the augmented
$`U`$, the $`\delta_j`$ are **conditionally independent across $`j`$**
(the softmax coupling of the HMNL does not exist here), so the draws
parallelize validly:

``` math
\delta_j \mid \cdot \sim N\left(\frac{\sum_{r \in j}(U_r - \mu_x(r))/\sigma^2 + z_j'\theta/\sigma_d^2}{n_j/\sigma^2 + 1/\sigma_d^2},\; \big(n_j/\sigma^2 + 1/\sigma_d^2\big)^{-1}\right),
```

each $`j`$’s sufficient statistic summed over its CSR row list in fixed
row order by a single thread (bitwise thread-count invariance).

**(c) RSS pass** (work-shared, fixed order) with the *new* $`\beta`$ and
$`\delta`$ — a separate pass so $`\sigma^2`$ conditions on current
values, not the stale $`\delta`$ from (a). Outside latents contribute
$`U_{iot}^2`$ (their residual is the latent itself) and count in
$`n_{\text{latents}}`$.

**(d) Hierarchy** (master): $`b`$, $`W`$, $`\theta`$,
$`(\sigma_d^2, a_d)`$ exactly as in the HMNL, then

``` math
\sigma^2 \mid U, \beta, \delta \sim \mathcal{IG}\big(a_0 + \tfrac12 n_{\text{latents}},\; s_0 + \tfrac12 \text{RSS}\big).
```

**Recording.** Raw chains for
$`(b, W, \delta, \theta, \sigma_d^2, \sigma^2)`$; the
$`\beta_i / \delta / \xi`$ summaries are accumulated on the **identified
scale** (divided by the current $`\sigma`$) at recording time.

## 5. Sampling Primitives

The truncated-normal draw is the three-regime sampler in
`src/bayes_samplers.h` (naive rejection / Robert exponential rejection /
inverse-CDF), drawing exclusively from per-(iteration, unit)
Xoshiro256++ streams. The worker-side $`\beta_i`$ draw uses the
hand-rolled Cholesky machinery of `src/hb_internal.h` with per-thread
scratch — never BLAS/LAPACK off the master thread.

## 6. Post-Processing and Reported Quantities

As for the HMNL, with two differences. Choice probabilities have no
closed form; with i.i.d. shocks they reduce to the 1-D integral

``` math
P(j) = \int \phi(u) \prod_{k \ne j} \Phi(V_j - V_k + u)\, du \qquad (V_o = 0,\ \text{identified scale } \sigma = 1),
```

evaluated by fixed-node Gauss–Hermite quadrature in
[`predict.choicer_hb()`](https://fpcordeiro.github.io/choicer/reference/predict.choicer_hb.md)
(deterministic; verified against brute-force Monte Carlo argmax). And
[`logsum()`](https://fpcordeiro.github.io/choicer/reference/logsum.md) /
[`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md)
refuse with an informative error: the logsum formula is EV1-specific,
probit $`E[\max U]`$ has no closed form, and the simulated-Emax variant
is roadmapped — the HMNL is the welfare workhorse.

## 7. Implementation Details

Same engine contract, persistent-region structure, abort protocol, and
thread-invariance discipline as the HMNL (see its §7). The RNG partition
adds tag $`N + i`$ for respondent $`i`$’s $`\beta_i`$ regression draw;
the master block remains based at $`2N + J`$ ($`\sigma_d^2`$ and
$`\sigma^2`$ share the last tag’s stream in fixed order). The fused pass
writes only respondent-owned slots ($`U`$, $`U_{\text{out}}`$,
$`\mu_x`$, $`\beta_i`$ column, RSS slot), the $`\delta`$ pass only
$`\delta_j`$, and the RSS reduction happens on the master thread in
fixed respondent order.

## References

- Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
  polychotomous response data. *JASA*, 88(422), 669–679.
- Liu, J. S., & Wu, Y. N. (1999). Parameter expansion for data
  augmentation. *JASA*, 94(448), 1264–1274.
- McCulloch, R., & Rossi, P. E. (1994). An exact likelihood analysis of
  the multinomial probit model. *Journal of Econometrics*, 64(1–2),
  207–240.
- Rossi, P. E., Allenby, G. M., & McCulloch, R. (2005). *Bayesian
  Statistics and Marketing*. Wiley.
- Train, K. (2009). *Discrete Choice Methods with Simulation* (2nd ed.).
  Cambridge University Press.
