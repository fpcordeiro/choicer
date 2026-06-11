# Mathematical Documentation: Bayesian Multinomial Probit Model

This document provides a detailed mathematical description of the Bayesian multinomial probit (MNP) model as implemented in `src/mnprobit.cpp`, with sampling primitives in `src/bayes_samplers.h` and the RNG core in `src/rng.h`.

## Table of Contents

1. [Notation](#notation)
2. [Model Definition](#1-model-definition)
3. [Identification](#2-identification)
4. [Priors and Joint Posterior](#3-priors-and-joint-posterior)
5. [The Gibbs Sampler](#4-the-gibbs-sampler)
6. [Sampling Primitives](#5-sampling-primitives)
7. [Post-Processing and Reported Quantities](#6-post-processing-and-reported-quantities)
8. [Implementation Details](#7-implementation-details)
9. [References](#references)

---

## Notation

| Symbol | Description |
|--------|-------------|
| $i = 1, \ldots, N$ | Index for individuals (choice situations) |
| $j = 1, \ldots, J$ | Index for alternatives; alternative 1 is the **base** alternative |
| $p = J - 1$ | Dimension of the utility differences |
| $y_i \in \{0, 1, \ldots, p\}$ | Choice of individual $i$: $0$ for the base alternative, $j$ for the $j$-th non-base alternative |
| $w_i$ | Latent utility differences of individual $i$ ($p \times 1$) |
| $X_i$ | Differenced design matrix of individual $i$ ($p \times K$) |
| $\beta$ | Coefficient vector ($K \times 1$) |
| $\Sigma$ | Covariance matrix of the differenced errors ($p \times p$) |
| $\Omega = \Sigma^{-1}$ | Precision matrix |
| $\bar{\beta}, A$ | Prior mean and prior precision of $\beta$ |
| $\nu, V$ | Inverse-Wishart prior degrees of freedom and scale |
| $R$, `burn`, `thin` | Total Gibbs iterations, burn-in, thinning interval |
| `seed` | Master RNG seed (all streams derive from it) |
| $\sigma_{11}$ | The $(1,1)$ element of $\Sigma$, used for identification |

---

## 1. Model Definition

### 1.1 Utility Specification and Differencing

Each individual $i$ derives latent utility $U_{ij} = Z_{ij}\beta + e_{ij}$ from alternative $j$, with jointly normal errors, and chooses the alternative with the highest utility. Only utility *differences* matter for the choice, so the model is estimated in differences against the base alternative (alternative 1 in the package's integer coding):

$$
w_{ij} = U_{i,j+1} - U_{i,1}, \qquad j = 1, \ldots, p,
$$

which gives the $p$-dimensional latent regression

$$
w_i = X_i \beta + \varepsilon_i, \qquad \varepsilon_i \sim N_p(0, \Sigma),
$$

where row $j$ of $X_i$ contains the differenced covariates $Z_{i,j+1} - Z_{i,1}$.

### 1.2 Choice Rule

The observed choice is determined by the latent differences:

$$
y_i = \begin{cases}
j & \text{if } w_{ij} > \max\!\left(0, \max_{k \neq j} w_{ik}\right) \\
0 \;(\text{base}) & \text{if } \max_j w_{ij} < 0.
\end{cases}
$$

The "0" threshold is the base alternative's (differenced) utility.

*Code reference: [mnprobit.cpp:233-243](../src/mnprobit.cpp#L233-L243)*

### 1.3 Alternative-Specific Constants and the Base Alternative

ASCs enter the differenced model as $p$ intercepts, one per non-base alternative: the ASC column for alternative $j$ is an indicator that equals 1 in row $j$ of every $X_i$. This mirrors the package's reference-alternative convention for the frequentist models (first sorted alternative is the reference); the `base_alt` argument moves a chosen label to the front of the factor levels, exactly as `outside_opt_label` does for `prepare_mnl_data()`.

Balanced choice sets are required: every choice situation must contain the same $J$ alternatives. A user who wants an outside good includes it as explicit rows with zero covariates and sets `base_alt` to its label.

*Code reference: [mnprobit_utils.R (prepare_mnp_data)](../R/mnprobit_utils.R)*

---

## 2. Identification

### 2.1 Scale Invariance

The choice rule is invariant to a common rescaling of the latent utilities: for any $c > 0$,

$$
(\beta, \Sigma) \quad \text{and} \quad (c\,\beta, \; c^2\,\Sigma)
$$

imply identical choice probabilities. One scale restriction is therefore required.

### 2.2 The Two Standard Strategies

1. **Fully identified sampler** (McCulloch, Polson & Rossi 2000): restrict $\sigma_{11} = 1$ inside the sampler, placing a prior directly on the identified parameters. This requires a non-standard prior on a decomposition of $\Sigma$ and a non-conjugate conditional, and the resulting chain is known to mix more slowly.

2. **Non-identified chain with post-processing** (McCulloch & Rossi 1994; the `bayesm::rmnpGibbs` approach): run the Gibbs sampler on the unrestricted $(\beta, \Sigma)$ with a conjugate inverse-Wishart prior, and report the identified quantities

$$
\tilde{\beta}^{(r)} = \frac{\beta^{(r)}}{\sqrt{\sigma_{11}^{(r)}}},
\qquad
\tilde{\Sigma}^{(r)} = \frac{\Sigma^{(r)}}{\sigma_{11}^{(r)}},
$$

computed **per draw** $r$.

### 2.3 The Package Default

`choicer` implements strategy 2. Rationale: every Gibbs conditional remains standard (truncated normal, multivariate normal, inverse-Wishart), the chain navigates the posterior more freely, and the reported quantities $\tilde{\beta}, \tilde{\Sigma}$ are fully identified. The caveat, stated in McCulloch, Polson & Rossi (2000), is that the prior on the *identified* parameters is induced by the prior on the non-identified ones rather than specified directly; with the diffuse defaults below this is innocuous in practice.

*Code reference: [mnprobit_utils.R (run_mnprobit post-processing)](../R/mnprobit_utils.R)*

---

## 3. Priors and Joint Posterior

### 3.1 Priors

$$
\beta \sim N(\bar{\beta}, A^{-1}), \qquad \Sigma \sim \mathrm{IW}(\nu, V),
$$

where the inverse-Wishart density is parameterized (matching `bayesm`) as

$$
\pi(\Sigma) \propto |\Sigma|^{-(\nu + p + 1)/2} \exp\!\left(-\tfrac{1}{2}\operatorname{tr}(V \Sigma^{-1})\right),
\qquad
E[\Sigma] = \frac{V}{\nu - p - 1}.
$$

Defaults: $\bar{\beta} = 0$, $A = 0.01\, I_K$, $\nu = p + 3$, $V = \nu I_p$ — proper but diffuse, with prior mean of $\Sigma$ close to $I_p$.

### 3.2 Augmented Joint Posterior

With the latent $w = (w_1, \ldots, w_N)$ treated as additional unknowns (Albert & Chib 1993), the joint posterior is

$$
\pi(\beta, \Sigma, w \mid y) \;\propto\;
\prod_{i=1}^{N} \Big[ \phi_p(w_i \mid X_i \beta, \Sigma)\; \mathbf{1}\{w_i \in B(y_i)\} \Big]\;
\pi(\beta)\, \pi(\Sigma),
$$

where $\phi_p$ is the $p$-variate normal density and $B(y_i)$ is the cone implied by the choice rule of §1.2. All three blocks of full conditionals are standard, which is exactly what data augmentation buys: no MNP choice probabilities (multivariate normal rectangle probabilities) are ever evaluated.

---

## 4. The Gibbs Sampler

One iteration cycles through three blocks.

### 4.1 Latent Utilities: $w_{ij} \mid w_{i,-j}, \beta, \Sigma, y_i$

Given $(\beta, \Sigma)$, the $w_i$ are independent across $i$, and each component $w_{ij}$ given the others is univariate normal, truncated to the region implied by $y_i$. Writing $\mu_i = X_i\beta$ and using the precision matrix $\Omega = \Sigma^{-1}$, the untruncated conditional moments are

$$
\tau_j^2 = \frac{1}{\Omega_{jj}},
\qquad
m_{ij} = \mu_{ij} - \tau_j^2 \sum_{k \neq j} \Omega_{jk}\,(w_{ik} - \mu_{ik}).
$$

The truncation bounds, using the **current** values of the other components, are:

| Case | Lower bound $a$ | Upper bound $b$ |
|------|-----------------|-----------------|
| $y_i = j$ (this component chosen) | $\max\!\left(0, \max_{k \neq j} w_{ik}\right)$ | $+\infty$ |
| $y_i = 0$ (base chosen) | $-\infty$ | $0$ |
| $y_i = c$, $c \neq j$ (another non-base chosen) | $-\infty$ | $w_{ic}$ |

Each $w_{ij}$ is drawn from $\mathrm{TN}(m_{ij}, \tau_j^2; a, b)$ and immediately replaces the old value (a within-$i$ Gibbs sweep over $j = 1, \ldots, p$).

Because the $w_i$ are conditionally independent across $i$, the sweep over choice situations is parallelized with OpenMP — this is an exact Gibbs step, not an approximation. The within-$i$ sweep over $j$ stays sequential.

*Code reference: [mnprobit.cpp:204-252](../src/mnprobit.cpp#L204-L252)*

### 4.2 Coefficients: $\beta \mid w, \Sigma$

Stacking the latent regressions and applying standard conjugate GLS algebra:

$$
\beta \mid w, \Sigma \;\sim\; N\!\left(\tilde{\beta},\, Q^{-1}\right),
\qquad
Q = A + \sum_{i=1}^{N} X_i' \Omega X_i,
\qquad
\tilde{\beta} = Q^{-1}\!\left(A\bar{\beta} + \sum_{i=1}^{N} X_i' \Omega\, w_i\right).
$$

Two computational identities keep this cheap:

1. **Precomputed Gram blocks.** Let $X_{(j)}$ be the $N \times K$ matrix collecting the component-$j$ rows of all $X_i$, and $G_{jk} = X_{(j)}' X_{(k)}$. Then

$$
\sum_{i} X_i' \Omega X_i = \sum_{j=1}^{p}\sum_{k=1}^{p} \Omega_{jk}\, G_{jk}.
$$

The $G_{jk}$ do not depend on the chain and are computed once at startup; assembling $Q$ then costs $O(p^2 K^2)$ per iteration, **independent of $N$**.

2. **No inverse.** With the Cholesky factorization $Q = U'U$ ($U$ upper triangular), the mean solves two triangular systems and the draw is

$$
\beta = \tilde{\beta} + U^{-1} z, \qquad z \sim N(0, I_K),
$$

so $\operatorname{Var}(U^{-1}z) = Q^{-1}$ without ever forming $Q^{-1}$.

*Code reference: [mnprobit.cpp:153-176](../src/mnprobit.cpp#L153-L176) (Gram blocks), [mnprobit.cpp:254-278](../src/mnprobit.cpp#L254-L278) (draw)*

### 4.3 Covariance: $\Sigma \mid w, \beta$

With $\varepsilon_i = w_i - X_i\beta$ and $S = \sum_i \varepsilon_i \varepsilon_i'$ (one rank-$N$ BLAS product of the $p \times N$ residual matrix), conjugacy of the inverse-Wishart gives

$$
\Sigma \mid w, \beta \;\sim\; \mathrm{IW}\!\left(\nu + N,\; V + S\right).
$$

$V + S$ is symmetrized before factorization to remove floating-point asymmetry.

*Code reference: [mnprobit.cpp:280-292](../src/mnprobit.cpp#L280-L292)*

---

## 5. Sampling Primitives

All primitives are implemented from scratch (header-only) and draw from the package's own RNG, never from R's RNG.

### 5.1 Truncated Univariate Normal

To draw from $\mathrm{TN}(\mu, \sigma^2; a, b)$, standardize to bounds $(\alpha, \beta_u) = ((a-\mu)/\sigma, (b-\mu)/\sigma)$ and mirror the interval into the upper half-line when $\beta_u < 0$ (by symmetry of the normal). Two regimes:

- **Central case** ($\alpha < 4$): inverse CDF on the *upper-tail* scale,

$$
z = \bar{\Phi}^{-1}\!\big(u\big), \qquad u \sim \mathrm{Uniform}\big(\bar{\Phi}(\beta_u),\, \bar{\Phi}(\alpha)\big),
$$

where $\bar{\Phi}(x) = 1 - \Phi(x)$. Working with upper-tail probabilities keeps full floating-point precision when both bounds are far in the tail, where $\Phi(x)$ rounds to 1. (The naive $\Phi^{-1}(\mathrm{Uniform}(\Phi(a), \Phi(b)))$ returns $\infty$ for $a \gtrsim 8$.)

- **Tail case** ($\alpha \geq 4$): Robert (1995) exponential rejection. With rate

$$
\lambda = \frac{\alpha + \sqrt{\alpha^2 + 4}}{2},
$$

propose $z = \alpha + E/\lambda$ with $E \sim \mathrm{Exp}(1)$, reject $z > \beta_u$, and accept with probability $\exp\!\left(-(z - \lambda)^2 / 2\right)$. This is exact and has acceptance probability approaching 1 as $\alpha \to \infty$.

$\bar{\Phi}$ and $\bar{\Phi}^{-1}$ are evaluated with Rmath's `pnorm`/`qnorm` (pure functions, thread-safe — only R's *RNG* is unusable in threads).

*Code reference: [bayes_samplers.h:36-72](../src/bayes_samplers.h#L36-L72)*

### 5.2 Standard Normal, Gamma, Chi-Squared

- **Normal**: Marsaglia polar method — exact, no tables; generates pairs and caches the spare. *Code reference: [rng.h:69-85](../src/rng.h#L69-L85)*
- **Gamma**$(a, 1)$: Marsaglia & Tsang (2000) squeeze method — exact for any $a \geq 1$; for $a < 1$ the boosting identity $X = \mathrm{Gamma}(a+1)\, U^{1/a}$ is used. *Code reference: [rng.h:93-113](../src/rng.h#L93-L113)*
- **Chi-squared**$(k)$ $= 2\,\mathrm{Gamma}(k/2, 1)$, valid for non-integer $k$ (needed by Bartlett). *Code reference: [rng.h:115-117](../src/rng.h#L115-L117)*

### 5.3 Multivariate Normal

$x = \mu + Lz$ with $L$ the lower Cholesky factor of $\Sigma$ and $z \sim N(0, I)$.

*Code reference: [bayes_samplers.h:75-86](../src/bayes_samplers.h#L75-L86)*

### 5.4 Wishart and Inverse-Wishart

**Wishart**$(\mathrm{df}, S)$ via the Bartlett decomposition: with $L = \mathrm{chol}(S)$ (lower) and $T$ lower triangular where

$$
T_{ii} = \sqrt{\chi^2_{\mathrm{df} - i + 1}} \;\; (1\text{-based } i), \qquad T_{ij} \sim N(0, 1) \;\; (j < i),
$$

the draw is $W = (LT)(LT)'$. **Inverse-Wishart**$(\mathrm{df}, V)$: draw $W \sim \mathrm{Wishart}(\mathrm{df}, V^{-1})$ and return $W^{-1}$ (symmetrized).

*Code reference: [bayes_samplers.h:92-130](../src/bayes_samplers.h#L92-L130)*

---

## 6. Post-Processing and Reported Quantities

The kept draws $(\beta^{(r)}, \Sigma^{(r)})$ live on the non-identified scale. The R wrapper computes, for every kept draw,

$$
\tilde{\beta}^{(r)} = \frac{\beta^{(r)}}{\sqrt{\sigma_{11}^{(r)}}},
\qquad
\tilde{\Sigma}^{(r)} = \frac{\Sigma^{(r)}}{\sigma_{11}^{(r)}},
$$

and all reported summaries — posterior means (`coef()`), posterior standard deviations, posterior covariance (`vcov()`), and equal-tailed credible intervals (`summary()`) — are computed on these identified draws.

The normalization must be per draw: the identified quantity is a nonlinear function of the chain state, so

$$
E\!\left[\frac{\beta}{\sqrt{\sigma_{11}}}\right] \;\neq\; \frac{E[\beta]}{\sqrt{E[\sigma_{11}]}} ,
$$

and normalizing posterior means after averaging would be biased. By construction $\tilde{\sigma}_{11}^{(r)} = 1$ for every draw.

*Code reference: [mnprobit_utils.R (run_mnprobit post-processing)](../R/mnprobit_utils.R)*

---

## 7. Implementation Details

### 7.1 Data Layout and Parameter Vector

`prepare_mnp_data()` differences the covariates against the base alternative and stacks the result into a single $(N \cdot p) \times K$ matrix `X`, rows grouped by choice situation with components in alternative order. The coefficient vector is $\theta = [\beta_{\text{covariates}}, \mathrm{ASC}_2, \ldots, \mathrm{ASC}_J]$ ($K = K_x + p$ columns when `use_asc = TRUE`); the C++ engine sees only the generic latent regression and has no ASC logic. Alternative-invariant covariates difference into the span of the ASC columns and are removed by the collinearity check.

`y` is an integer $N$-vector with 0 for the base and $j \in \{1, \ldots, p\}$ for the $j$-th non-base alternative.

### 7.2 In-Place Chain State

The latent matrix $W$ ($p \times N$) is allocated once and swept in place for the entire run — iteration $r+1$ starts from iteration $r$'s values, which is precisely the Markov property. The systematic-utility matrix `Mu`, the residual matrix, and all $\beta$-draw workspaces are likewise pre-allocated; the hot loop performs no heap allocation. Where a stacked-vector view of a $p \times N$ matrix is needed (e.g. $X'\mathrm{vec}(\Omega W)$), Armadillo's non-copying alias constructor is used (column-major $p \times N$ memory *is* the stacked $N p$-vector).

*Code reference: [mnprobit.cpp:178-201](../src/mnprobit.cpp#L178-L201)*

### 7.3 RNG Streams and Reproducibility

The RNG core is xoshiro256++ seeded via splitmix64 (the generator authors' recommended initialization). Every consumer derives an independent stream from the master seed:

- the latent update of observation $i$ in iteration $r$ uses stream $(\texttt{seed}, r, i)$;
- the $\beta$ and $\Sigma$ draws of iteration $r$ use the tagged streams $(\texttt{seed}, r, N)$ and $(\texttt{seed}, r, N+1)$.

Streams are deterministic functions of their key, so the draws are **bitwise reproducible regardless of the number of OpenMP threads or the loop schedule** — verified by a thread-invariance test (1 vs 4 threads). When the user does not pass `mcmc$seed`, the master seed is drawn from R's RNG, so `set.seed()` governs the run end-to-end, consistent with the rest of the package.

*Code reference: [rng.h:119-131](../src/rng.h#L119-L131), [mnprobit.cpp:217](../src/mnprobit.cpp#L217)*

### 7.4 Parallelization

The latent sweep runs under `#pragma omp parallel for schedule(static)` over choice situations; threads write disjoint columns of $W$ and share read-only $\Omega$ and `Mu`, so no synchronization is needed. Thread count follows the package-wide convention (`set_num_threads()` / `get_num_threads()`). The $\beta$ and $\Sigma$ draws, `Rcpp::checkUserInterrupt()`, and all error paths stay on the master thread. Unlike the frequentist likelihood engines, the Gibbs chain itself is inherently sequential across iterations; only the within-iteration sweep parallelizes.

*Code reference: [mnprobit.cpp:208-216](../src/mnprobit.cpp#L208-L216)*

### 7.5 Draw Storage

`betadraw` is $R_{\text{keep}} \times K$ with $R_{\text{keep}} = \lceil (R - \texttt{burn}) / \texttt{thin} \rceil$. `sigmadraw` is $R_{\text{keep}} \times p(p+1)/2$, storing the lower triangle of $\Sigma$ in row-major order $(1,1), (2,1), (2,2), (3,1), \ldots$ — the same `vech_row()` convention used by the mixed logit code, with display names `Sigma_ij`.

*Code reference: [mnprobit.cpp:294-305](../src/mnprobit.cpp#L294-L305)*

### 7.6 Per-Iteration Complexity

| Step | Cost | Notes |
|------|------|-------|
| Latent sweep | $O(N p^2)$ | parallelized over $i$ |
| $\Omega W$, residual product | $O(N p^2)$ | single BLAS calls |
| $X'\mathrm{vec}(\cdot)$ | $O(N p K)$ | single BLAS gemv |
| Assemble $Q$ | $O(p^2 K^2)$ | independent of $N$ (Gram blocks) |
| Cholesky / solves for $\beta$ | $O(K^3)$ | |
| $\Sigma$ draw | $O(p^3)$ | |

Total cost is linear in $N$ with no $N K^2$ term, and the latent sweep — the dominant term for large $N$ — scales near-linearly with threads.

---

## References

- Albert, J. H., & Chib, S. (1993). Bayesian Analysis of Binary and Polychotomous Response Data. *Journal of the American Statistical Association*, 88(422), 669–679.
- Allenby, G. M., & Rossi, P. E. (1999). Marketing models of consumer heterogeneity. *Journal of Econometrics*, 89(1–2), 57–78.
- Blackman, D., & Vigna, S. (2021). Scrambled Linear Pseudorandom Number Generators. *ACM Transactions on Mathematical Software*, 47(4), 1–32.
- Marsaglia, G., & Tsang, W. W. (2000). A Simple Method for Generating Gamma Variables. *ACM Transactions on Mathematical Software*, 26(3), 363–372.
- McCulloch, R., & Rossi, P. E. (1994). An exact likelihood analysis of the multinomial probit model. *Journal of Econometrics*, 64(1–2), 207–240.
- McCulloch, R. E., Polson, N. G., & Rossi, P. E. (2000). A Bayesian analysis of the multinomial probit model with fully identified parameters. *Journal of Econometrics*, 99(1), 173–193.
- Robert, C. P. (1995). Simulation of truncated normal variables. *Statistics and Computing*, 5(2), 121–125.
- Rossi, P. E., Allenby, G. M., & McCulloch, R. (2005). *Bayesian Statistics and Marketing*. John Wiley & Sons.
- Train, K. E. (2009). *Discrete Choice Methods with Simulation* (2nd ed.). Cambridge University Press. Chapter 5.
