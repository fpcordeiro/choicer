# The math behind choicer: multinomial logit

This document provides a detailed mathematical description of the
multinomial logit (MNL) model as implemented in `src/mnlogit.cpp`.

## Table of Contents

1.  [Notation](#notation)
2.  [Model Definition](#id_1-model-definition)
3.  [Log-Likelihood Function](#id_2-log-likelihood-function)
4.  [Gradient Computation](#id_3-gradient-computation)
5.  [Hessian Computation](#id_4-hessian-computation)
6.  [Elasticity Computation](#id_5-elasticity-computation)
7.  [Diversion Ratio Computation](#id_6-diversion-ratio-computation)
8.  [BLP Contraction Mapping](#id_7-blp-contraction-mapping)
9.  [Willingness to Pay](#id_8-willingness-to-pay)
10. [Consumer Surplus and the
    Logsum](#id_9-consumer-surplus-and-the-logsum)
11. [Goodness of Fit](#id_10-goodness-of-fit)
12. [Choice-Based Sampling and WESML
    Weighting](#id_11-choice-based-sampling-and-wesml-weighting)
13. [Implementation Details](#id_12-implementation-details)

------------------------------------------------------------------------

## Notation

| Symbol | Description |
|----|----|
| $`i = 1, \ldots, N`$ | Index for individuals (choice situations) |
| $`j = 1, \ldots, J_i`$ | Index for alternatives available to individual $`i`$ |
| $`j_i`$ | The alternative chosen by individual $`i`$ |
| $`w_i`$ | Weight for individual $`i`$ |
| $`X_{ij}`$ | Row vector of covariates for individual $`i`$ and alternative $`j`$ ($`1 \times K`$) |
| $`\beta`$ | Coefficient vector ($`K \times 1`$) |
| $`\delta_j`$ | Alternative-specific constant (ASC) for alternative $`j`$ |
| $`V_{ij}`$ | Systematic (deterministic) utility |
| $`P_{ij}`$ | Probability that individual $`i`$ chooses alternative $`j`$ |

------------------------------------------------------------------------

## 1. Model Definition

### 1.1 Utility Specification

The utility that individual $`i`$ derives from alternative $`j`$ is:

``` math
U_{ij} = V_{ij} + \varepsilon_{ij}
```

where the systematic utility is:

``` math
V_{ij} = X_{ij}\beta + \delta_j
```

and $`\varepsilon_{ij}`$ is an i.i.d. Type I Extreme Value (Gumbel)
error term with location 0 and scale 1.

### 1.2 Alternative-Specific Constants (ASCs)

The ASC parameters $`\delta_j`$ capture the average effect of unobserved
factors for each alternative. For identification, one ASC must be
normalized:

- **Without outside option** (`include_outside_option = FALSE`): The
  first inside alternative’s ASC is fixed to zero ($`\delta_1 = 0`$).
  The parameter vector contains $`J-1`$ free ASC parameters.

- **With outside option** (`include_outside_option = TRUE`): The outside
  option (alternative 0) has utility $`V_{i0} = 0`$, serving as the
  reference. All $`J`$ inside alternatives have free ASC parameters.

*Code reference:
[mnlogit.cpp:50-70](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L50-L70)*

------------------------------------------------------------------------

## 2. Log-Likelihood Function

### 2.1 Choice Probability

Under the Type I Extreme Value distribution assumption, the probability
that individual $`i`$ chooses alternative $`j`$ follows the multinomial
logit (softmax) formula:

``` math
P_{ij} = \frac{\exp(V_{ij})}{\sum_{k=1}^{J_i} \exp(V_{ik})}
```

If an outside option is included, the denominator includes
$`\exp(V_{i0}) = \exp(0) = 1`$.

*Code reference:
[mnlogit.cpp:119-122](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L119-L122)*

### 2.2 Log-Likelihood

The log-likelihood is the weighted sum of the log-probabilities of the
observed choices:

``` math
\ell(\theta) = \sum_{i=1}^{N} w_i \log P_{ij_i}
```

where $`\theta = (\beta, \delta)`$ is the full parameter vector and
$`j_i`$ is the alternative chosen by individual $`i`$.

Expanding the log-probability:

``` math
\log P_{ij_i} = V_{ij_i} - \log\left(\sum_{k=1}^{J_i} \exp(V_{ik})\right)
```

*Code reference:
[mnlogit.cpp:133-142](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L133-L142)*

### 2.3 Log-Sum-Exp Trick for Numerical Stability

To prevent numerical overflow when computing the denominator, the
implementation subtracts the maximum utility before exponentiation:

``` math
\sum_k \exp(V_{ik}) = \exp(V_{\max}) \sum_k \exp(V_{ik} - V_{\max})
```

where $`V_{\max} = \max_k V_{ik}`$.

This gives:

``` math
\log\left(\sum_k \exp(V_{ik})\right) = V_{\max} + \log\left(\sum_k \exp(V_{ik} - V_{\max})\right)
```

The log-probability becomes:

``` math
\log P_{ij_i} = (V_{ij_i} - V_{\max}) - \log\left(\sum_k \exp(V_{ik} - V_{\max})\right)
```

*Code reference:
[mnlogit.cpp:119-122](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L119-L122)*

------------------------------------------------------------------------

## 3. Gradient Computation

### 3.1 General Formula

The gradient of the log-probability of the chosen alternative with
respect to utility is:

``` math
\frac{\partial \log P_{ij_i}}{\partial V_{ia}} = \mathbf{1}_{a = j_i} - P_{ia}
```

where $`\mathbf{1}_{a = j_i}`$ is the indicator function (1 if $`a`$ is
the chosen alternative, 0 otherwise).

### 3.2 Gradient with Respect to $`\beta`$

Since $`\frac{\partial V_{ij}}{\partial \beta} = X_{ij}^T`$, the
contribution to the gradient from individual $`i`$ is:

``` math
\frac{\partial \ell_i}{\partial \beta} = w_i \sum_{j=1}^{J_i} \left(\mathbf{1}_{j = j_i} - P_{ij}\right) X_{ij}^T
```

This simplifies to:

``` math
\frac{\partial \ell_i}{\partial \beta} = w_i \left( X_{ij_i}^T - \sum_{j=1}^{J_i} P_{ij} X_{ij}^T \right)
```

The full gradient is:

``` math
\nabla_\beta \ell = \sum_{i=1}^{N} w_i \left( X_{ij_i}^T - \sum_{j=1}^{J_i} P_{ij} X_{ij}^T \right)
```

**Implementation note.** Define the difference vector
$`d_i \in \mathbb{R}^{J_i}`$ with entries
$`d_{ij} = \mathbf{1}_{j = j_i} - P_{ij}`$. Then the sum above is
$`X_i^T d_i`$, a single BLAS matrix-vector product
(`X_i.t() * diff_vec`). When an outside option is present, $`X_i`$
covers only the $`M_i`$ inside alternatives and $`d_i`$ is sliced
accordingly (`diff_vec.subvec(1, m_i)`).

*Code reference:
[mnlogit.cpp:146-156](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L146-L156)*

### 3.3 Gradient with Respect to $`\delta`$

Since
$`\frac{\partial V_{ij}}{\partial \delta_a} = \mathbf{1}_{j = a}`$, the
contribution to the gradient from individual $`i`$ for alternative $`a`$
is:

``` math
\frac{\partial \ell_i}{\partial \delta_a} = w_i \left(\mathbf{1}_{j_i = a} - P_{ia}\right)
```

The gradient is computed only for the free ASC parameters: - Without
outside option: $`\delta_1 = 0`$ (fixed), so we compute for
$`a = 2, \ldots, J`$ - With outside option: we compute for all inside
alternatives $`a = 1, \ldots, J`$

The delta gradient reuses the same difference vector $`d_i`$ computed
for the beta block; its entries are scattered into the appropriate
positions of the gradient vector via an irregular index mapping.

*Code reference:
[mnlogit.cpp:158-172](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L158-L172)*

------------------------------------------------------------------------

## 4. Hessian Computation

### 4.1 Analytical Hessian — Full-Vector Form

The Hessian of the log-likelihood for individual $`i`$ is:

``` math
H_i = \frac{\partial^2 \ell_i}{\partial \theta \,\partial \theta'} = -w_i \left[ \sum_{a \in A_i} P_{ia} Z_a Z_a^T - \left(\sum_{a \in A_i} P_{ia} Z_a\right)\left(\sum_{a \in A_i} P_{ia} Z_a\right)^T \right]
```

where $`A_i`$ is individual $`i`$’s choice set (including the outside
option when present) and $`Z_a \in \mathbb{R}^{K + J_\delta}`$ is the
“score vector” for alternative $`a`$:

``` math
Z_a = \begin{pmatrix} x_a \\ e_{d(a)} \end{pmatrix}
```

Here $`x_a \in \mathbb{R}^K`$ is the covariate vector for alternative
$`a`$, $`e_{d(a)}`$ is the standard basis vector in
$`\mathbb{R}^{J_\delta}`$ indicating the free ASC position $`d(a)`$ for
alternative $`a`$, and $`J_\delta`$ is the number of free ASC parameters
($`J-1`$ without outside option, $`J`$ with outside option). For the
outside option itself, $`Z_0 = 0`$ (no covariates, no free ASC).

This can be written compactly as:

``` math
H_i = -w_i \operatorname{Var}_{P_i}(Z)
```

where the variance is taken under the discrete distribution
$`P_i = (P_{i0}, P_{i1}, \ldots, P_{iJ_i})`$. The full Hessian is
$`\nabla^2_\theta \ell = \sum_{i=1}^N H_i`$ and the function returns
$`-\nabla^2_\theta \ell`$ (Hessian of the negative log-likelihood).

### 4.2 Block Decomposition

The parameter vector is $`\theta = (\beta, \delta)`$ with
$`\beta \in \mathbb{R}^K`$ (covariate coefficients) and
$`\delta \in \mathbb{R}^{J_\delta}`$ (free ASC parameters). Because
$`Z_a`$ has at most $`K+1`$ nonzero entries — $`K`$ in the $`\beta`$
block and at most one in the $`\delta`$ block — the dense
$`n_\mathrm{params} \times n_\mathrm{params}`$ outer product
$`Z_a Z_a^T`$ has only three nontrivial submatrices. Summing over
alternatives, the individual Hessian decomposes into three blocks:

``` math
-H_i / w_i = \begin{pmatrix} \mathbf{BB}_i & \mathbf{BD}_i \\ \mathbf{BD}_i^T & \mathbf{DD}_i \end{pmatrix}
```

**Beta-beta block** ($`K \times K`$, symmetric):

``` math
\mathbf{BB}_i = \sum_{a=1}^{m_i} p_{ia} \, x_a x_a^T - \mu_\beta \mu_\beta^T, \qquad \mu_\beta = \sum_{a=1}^{m_i} p_{ia} \, x_a
```

where the sum runs over the $`m_i`$ inside alternatives only (since
$`Z_0 = 0`$). This is the probability-weighted covariance of the
covariates:

``` math
\mathbf{BB}_i = \mathbb{E}_{P_i}[x x^T] - \mathbb{E}_{P_i}[x]\, \mathbb{E}_{P_i}[x]^T = \operatorname{Cov}_{P_i}(x)
```

**Delta-delta block** ($`J_\delta \times J_\delta`$, symmetric):

``` math
\mathbf{DD}_i = \operatorname{diag}(\mu_\delta) - \mu_\delta \mu_\delta^T, \qquad \mu_\delta[r] = \sum_{\{a : d(a)=r\}} p_{ia}
```

Because at most one alternative maps to each free ASC index,
$`\mu_\delta`$ is the vector of choice probabilities for ASC-bearing
alternatives. The diagonal of $`\mathbf{DD}_i`$ is
$`\mu_\delta \odot (1 - \mu_\delta)`$ and the off-diagonal $`(r,s)`$
entry is $`-\mu_\delta[r]\,\mu_\delta[s]`$.

**Beta-delta block** ($`K \times J_\delta`$, generally not symmetric):

``` math
\mathbf{BD}_i[:,r] = \sum_{\{a : d(a)=r\}} p_{ia} \, x_a - \mu_\beta \, \mu_\delta[r]
```

This block measures the probability-weighted covariance between the
covariates and the ASC indicators.

**Why the decomposition equals the full-vector form.** Partition
$`Z_a = (x_a^T, e_{d(a)}^T)^T`$. Then:

``` math
\sum_a p_{ia} Z_a Z_a^T = \begin{pmatrix} \sum_a p_{ia} x_a x_a^T & \sum_a p_{ia} x_a e_{d(a)}^T \\ \sum_a p_{ia} e_{d(a)} x_a^T & \sum_a p_{ia} e_{d(a)} e_{d(a)}^T \end{pmatrix}
```

``` math
\left(\sum_a p_{ia} Z_a\right)\left(\sum_a p_{ia} Z_a\right)^T = \begin{pmatrix} \mu_\beta \mu_\beta^T & \mu_\beta \mu_\delta^T \\ \mu_\delta \mu_\beta^T & \mu_\delta \mu_\delta^T \end{pmatrix}
```

Subtracting block by block yields exactly $`\mathbf{BB}_i`$,
$`\mathbf{BD}_i`$, $`\mathbf{DD}_i`$ above.

### 4.3 Symmetry Exploitation

$`\mathbf{BB}_i`$ and $`\mathbf{DD}_i`$ are symmetric. The
implementation accumulates only the upper triangles of these blocks
during the inner loop over alternatives and reflects them once after the
loop, halving the write operations for those blocks:

- **BB**: the double loop `for r in 0..K-1, for c in r..K-1` accumulates
  $`p_{ia} x_{ar} x_{ac}`$ into the upper triangle; the lower triangle
  is filled by mirroring.
- **DD**: $`\mu_\delta`$ is a length-$`J_\delta`$ vector; the matrix
  $`\operatorname{diag}(\mu_\delta) - \mu_\delta \mu_\delta^T`$ is
  assembled from the vector after all alternatives are processed, with
  upper triangle only and then mirrored.
- **BD**: $`K \times J_\delta`$ rectangular — the full block is computed
  (no symmetry available across the two index spaces).

The existing final symmetrization
`global_hess = 0.5*(global_hess + global_hess^T)` remains as a numerical
tidy-up.

**Complexity.** The block assembly costs $`O(J \cdot K^2 / 2)`$ for BB,
$`O(J \cdot K)`$ for BD scatter, and $`O(J_\delta^2 / 2)`$ for DD —
versus $`O(J \cdot (K + J_\delta)^2)`$ for the naive dense outer
product. The asymptotic gain scales as $`(1 + J_\delta/K)^2`$, which is
material when the number of free ASC parameters is comparable to or
larger than $`K`$.

### 4.4 ASC Normalization and Indexing

The free-delta position $`d(a)`$ for inside alternative $`a`$ (0-based
local index within individual $`i`$’s choice set, with
$`a \in \{0,\ldots,m_i-1\}`$) is:

``` math
d(a) = \begin{cases} \text{alt\_idx0\_i}[a] & \text{if } \texttt{include\_outside\_option} = \text{TRUE} \\ \text{alt\_idx0\_i}[a] - 1 & \text{if } \texttt{include\_outside\_option} = \text{FALSE and alt\_idx0\_i}[a] > 0 \\ \text{(skip — no ASC contribution)} & \text{if } \texttt{include\_outside\_option} = \text{FALSE and alt\_idx0\_i}[a] = 0 \end{cases}
```

where `alt_idx0_i[a]` is the 0-based global alternative ID from
`alt_idx0`. When there is no outside option, the first inside
alternative ($`\text{alt\_idx0\_i}[a] = 0`$) has its ASC fixed to zero
by the identification convention, so it contributes no free ASC entry
and is skipped in the BD and DD accumulation.

**Implementation note.** The optimization is internal and numerically
equivalent (within $`\approx 10^{-9}`$ in absolute terms, well within
the $`10^{-6}`$ tolerance verified at test time) to the prior dense
outer-product implementation. The public results — vcov, se, logLik, and
all post-estimation quantities — are unchanged.

*Code reference: `src/mnlogit.cpp`, function
`mnl_loglik_hessian_parallel`.*

------------------------------------------------------------------------

## 5. Elasticity Computation

### 5.1 Elasticity Definitions

The elasticity measures the percentage change in choice probability with
respect to a percentage change in an attribute. For attribute $`k`$ with
coefficient $`\beta_k`$:

**Own-Elasticity** (elasticity of $`P_{ij}`$ with respect to
$`x_{ijk}`$):

``` math
E_{jj}^k = \frac{\partial P_{ij}}{\partial x_{ijk}} \cdot \frac{x_{ijk}}{P_{ij}} = \beta_k \cdot x_{ijk} \cdot (1 - P_{ij})
```

**Cross-Elasticity** (elasticity of $`P_{ij}`$ with respect to
$`x_{imk}`$ where $`m \neq j`$):

``` math
E_{jm}^k = \frac{\partial P_{ij}}{\partial x_{imk}} \cdot \frac{x_{imk}}{P_{ij}} = -\beta_k \cdot x_{imk} \cdot P_{im}
```

### 5.2 Derivation

Starting from the choice probability:

``` math
P_{ij} = \frac{\exp(V_{ij})}{\sum_k \exp(V_{ik})}
```

Taking the derivative with respect to $`x_{imk}`$:

``` math
\frac{\partial P_{ij}}{\partial x_{imk}} = P_{ij} \left( \mathbf{1}_{j=m} - P_{im} \right) \beta_k
```

For own-elasticity ($`j = m`$):

``` math
\frac{\partial P_{ij}}{\partial x_{ijk}} = P_{ij} (1 - P_{ij}) \beta_k
```

Therefore:

``` math
E_{jj}^k = \beta_k \cdot x_{ijk} \cdot (1 - P_{ij})
```

For cross-elasticity ($`j \neq m`$):

``` math
\frac{\partial P_{ij}}{\partial x_{imk}} = -P_{ij} P_{im} \beta_k
```

Therefore:

``` math
E_{jm}^k = -\beta_k \cdot x_{imk} \cdot P_{im}
```

### 5.3 Aggregate Elasticities

The implementation computes weighted average elasticities across all
individuals:

``` math
\bar{E}_{jm}^k = \frac{\sum_{i=1}^{N} w_i \cdot E_{ijm}^k}{\sum_{i=1}^{N} w_i}
```

where $`E_{ijm}^k`$ is the elasticity for individual $`i`$.

*Code reference:
[mnlogit.cpp:906-930](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L906-L930)*

------------------------------------------------------------------------

## 6. Diversion Ratio Computation

### 6.1 Definition

The diversion ratio from alternative $`j`$ to alternative $`k`$ measures
the fraction of demand lost by $`j`$ that is captured by $`k`$ when
$`j`$ becomes less attractive (e.g., due to a price increase). It is a
key metric in antitrust analysis and merger simulation.

``` math
DR(j \to k) = -\frac{\partial Q_k / \partial p_j}{\partial Q_j / \partial p_j}
```

where $`Q_j = \sum_i w_i P_{ij}`$ is the (weighted) aggregate demand for
alternative $`j`$ and $`p_j`$ is the price of $`j`$.

### 6.2 Derivation for MNL

From the MNL choice probability derivatives (Section 5.2), the demand
derivatives with respect to price $`p_j`$ (with coefficient $`\beta_p`$)
are:

``` math
\frac{\partial Q_k}{\partial p_j} = -\beta_p \sum_{i=1}^{N} w_i \, P_{ij} \, P_{ik} \quad (k \neq j)
```

``` math
\frac{\partial Q_j}{\partial p_j} = \beta_p \sum_{i=1}^{N} w_i \, P_{ij} \, (1 - P_{ij})
```

Substituting into the diversion ratio formula:

``` math
DR(j \to k) = -\frac{-\beta_p \sum_i w_i \, P_{ij} \, P_{ik}}{\beta_p \sum_i w_i \, P_{ij} \, (1 - P_{ij})} = \frac{\sum_i w_i \, P_{ij} \, P_{ik}}{\sum_i w_i \, P_{ij} \, (1 - P_{ij})}
```

The price coefficient $`\beta_p`$ cancels, so the diversion ratio does
not depend on which covariate we differentiate with respect to. This is
a consequence of the IIA property.

### 6.3 Aggregate Diversion Ratio Matrix

The implementation computes the full $`J \times J`$ diversion ratio
matrix $`D`$ where:

``` math
D_{kj} = DR(j \to k) = \frac{\sum_{i=1}^{N} w_i \, P_{ij} \, P_{ik}}{\sum_{i=1}^{N} w_i \, P_{ij} \, (1 - P_{ij})} \quad (k \neq j)
```

``` math
D_{jj} = 0
```

### 6.4 Properties

**Column sums equal one.** For a given alternative $`j`$, the
off-diagonal entries in column $`j`$ sum to 1, meaning all diverted
demand is accounted for:

``` math
\sum_{k \neq j} DR(j \to k) = \frac{\sum_i w_i \, P_{ij} \sum_{k \neq j} P_{ik}}{\sum_i w_i \, P_{ij} (1 - P_{ij})} = \frac{\sum_i w_i \, P_{ij} (1 - P_{ij})}{\sum_i w_i \, P_{ij} (1 - P_{ij})} = 1
```

where we used $`\sum_{k \neq j} P_{ik} = 1 - P_{ij}`$.

**IIA proportionality.** When all individuals face the same choice set
and the same covariates (so $`P_{ij}`$ does not vary across $`i`$), the
diversion ratio simplifies to:

``` math
DR(j \to k) = \frac{P_j \, P_k}{P_j (1 - P_j)} = \frac{P_k}{1 - P_j} = \frac{s_k}{1 - s_j}
```

where $`s_j`$ is the market share of alternative $`j`$. This means
diversion is proportional to the receiving alternative’s market share,
independent of any characteristics of the losing alternative $`j`$ other
than its own share. This is a well-known limitation of MNL due to IIA.

### 6.5 Implementation

The function accumulates two quantities in parallel across individuals:

1.  **Numerator matrix**: $`N_{kj} = \sum_i w_i \, P_{ij} \, P_{ik}`$
    for each pair $`(k, j)`$
2.  **Denominator vector**:
    $`d_j = \sum_i w_i \, P_{ij} \, (1 - P_{ij})`$ for each $`j`$

The final matrix is computed as $`D_{kj} = N_{kj} / d_j`$ for
$`k \neq j`$, with $`D_{jj} = 0`$.

*Code reference:
[mnlogit.cpp:1085-1100](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L1085-L1100)*

------------------------------------------------------------------------

## 7. BLP Contraction Mapping

### 7.1 Problem Statement

Given observed market shares $`s_j`$, find the ASC parameters
$`\delta_j`$ such that the model-predicted shares match the observed
shares:

``` math
\hat{s}_j(\delta) = s_j \quad \forall j
```

where $`\hat{s}_j(\delta)`$ is the predicted market share for
alternative $`j`$.

### 7.2 Contraction Mapping Algorithm

Berry, Levinsohn, and Pakes (1995) show that the following iteration
converges to the solution:

``` math
\delta_j^{(t+1)} = \delta_j^{(t)} + \log(s_j) - \log(\hat{s}_j^{(t)})
```

This is equivalent to:

``` math
\delta^{(t+1)} = \delta^{(t)} + \log(s) - \log(\hat{s}^{(t)})
```

### 7.3 Implementation Details

1.  **Initialization**: Start with initial guess $`\delta^{(0)}`$
2.  **Prediction**: Compute predicted shares $`\hat{s}(\delta^{(t)})`$
    using `mnl_predict_shares_internal`
3.  **Update**: Apply the contraction:
    $`\delta^{(t+1)} = \delta^{(t)} + \log(s) - \log(\hat{s}^{(t)})`$
4.  **Convergence**: Check if
    $`\max_j |\delta_j^{(t+1)} - \delta_j^{(t)}| < \text{tol}`$
5.  **Normalization**: Subtract $`\delta_1`$ from all ASCs to maintain
    identification

*Code reference:
[mnlogit.cpp:547-567](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L547-L567)*

------------------------------------------------------------------------

## 8. Willingness to Pay

### 8.1 Definition

When one covariate is a price $`p`$ with coefficient
$`\alpha = \beta_p`$, the willingness to pay (WTP) for attribute $`k`$
is the marginal rate of substitution between the attribute and price —
the price change that leaves utility unchanged after a unit change in
the attribute:

``` math
\mathrm{WTP}_k = -\frac{\partial V / \partial x_k}{\partial V / \partial p} = -\frac{\beta_k}{\alpha}
```

Since utility is linear, the same ratio applies to ASCs: the WTP for
alternative $`j`$’s unobserved quality is $`-\delta_j / \alpha`$.

### 8.2 Delta-Method Standard Errors

WTP is a nonlinear function $`g(\theta) = -\theta_k / \theta_p`$ of the
estimated coefficients. By the delta method, with $`\hat{V}`$ the
estimated coefficient covariance matrix,

``` math
\widehat{\mathrm{Var}}(g(\hat\theta)) = \nabla g^T \, \hat{V}_{(k,p)} \, \nabla g,
\qquad
\nabla g = \begin{pmatrix} \partial g / \partial \theta_k \\ \partial g / \partial \theta_p \end{pmatrix}
= \begin{pmatrix} -1/\theta_p \\ \theta_k/\theta_p^2 \end{pmatrix}
```

where $`\hat{V}_{(k,p)}`$ is the $`2 \times 2`$ block of $`\hat{V}`$ for
$`(\theta_k, \theta_p)`$. The gradient is analytic (exact), so no
numerical differentiation is involved. Confidence intervals use the
normal approximation:

``` math
\mathrm{WTP}_k \pm z_{1-(1-\text{level})/2} \cdot \widehat{\mathrm{SE}}
```

Because the estimates and covariance are stored in natural (unscaled)
units, no scaling adjustment is needed even when the model was estimated
with `scale_vars`.

**Caveat.** The ratio of two asymptotically normal estimators has heavy
tails when the denominator is imprecisely estimated; the delta-method
interval is a first-order approximation that deteriorates as
$`|\alpha| / \mathrm{SE}(\hat\alpha)`$ falls. For a weakly identified
price coefficient, simulation methods (Krinsky–Robb) or Fieller
intervals are more reliable.

*Code reference:
[R/wtp.R](https://fpcordeiro.github.io/choicer/R/wtp.R)*

------------------------------------------------------------------------

## 9. Consumer Surplus and the Logsum

### 9.1 The Logsum (Expected Maximum Utility)

Under Type I Extreme Value errors, the expected maximum utility over
individual $`i`$’s choice set is, up to an additive constant (Euler’s
constant $`\gamma`$):

``` math
\mathbb{E}\left[\max_j U_{ij}\right] = \log\left(\sum_{j \in C_i} \exp(V_{ij})\right) + \gamma \equiv \mathrm{logsum}_i + \gamma
```

When the model includes an outside option with normalized utility
$`V_{i0} = 0`$, the sum includes its $`\exp(0) = 1`$ term. The
implementation uses the same max-subtraction trick as the log-likelihood
(Section 2.3).

### 9.2 Expected Consumer Surplus

With utility linear in price (no income effects), the marginal utility
of income is $`-\alpha`$ (positive for a negative price coefficient),
and expected consumer surplus in money units is (Train 2009, Ch. 3):

``` math
\mathbb{E}[CS_i] = \frac{\mathrm{logsum}_i}{-\alpha}
```

**Identification caveat.** Like utility itself, the logsum is only
defined up to an additive normalization (in particular the ASC
normalization), so CS *levels* are not interpretable on their own. The
economically meaningful quantity is the *difference* between scenarios,

``` math
\Delta \mathbb{E}[CS_i] = \frac{\mathrm{logsum}_i^{(1)} - \mathrm{logsum}_i^{(0)}}{-\alpha},
```

computed by evaluating the logsum on counterfactual data (`newdata`) and
on the baseline. The normalization constant cancels in the difference.

### 9.3 Delta-Method SE for the Mean Consumer Surplus

For the weighted mean $`m(\theta) = \sum_i w_i CS_i / \sum_i w_i`$, the
implementation reports a delta-method standard error. The building block
is the derivative of the logsum, which by the envelope-style identity is
a probability-weighted average of utility derivatives:

``` math
\frac{\partial\, \mathrm{logsum}_i}{\partial \theta_r}
= \sum_{j \in C_i} P_{ij} \, \frac{\partial V_{ij}}{\partial \theta_r}
```

with $`\partial V_{ij}/\partial \beta_r = x_{ijr}`$,
$`\partial V_{ij}/\partial \delta_a = \mathbf{1}_{j = a}`$, and a zero
contribution from the outside option row ($`V_{i0} \equiv 0`$). Then:

- For non-price parameters:
  $`\dfrac{\partial CS_i}{\partial \theta_r} = \dfrac{1}{-\alpha} \dfrac{\partial\, \mathrm{logsum}_i}{\partial \theta_r}`$.
- For the price coefficient, which enters both the logsum and the
  $`1/(-\alpha)`$ factor:

``` math
\frac{\partial CS_i}{\partial \alpha}
= \frac{\mathrm{logsum}_i}{\alpha^2}
+ \frac{1}{-\alpha} \sum_{j} P_{ij}\, x_{ij,p}
```

The SE is $`\sqrt{G^T \hat{V} G}`$, where $`G`$ is the weighted average
of the per-individual gradient vectors and $`\hat{V}`$ the full
coefficient covariance. All ingredients ($`P_{ij}`$, $`V_{ij}`$, $`X`$)
come from a single prediction pass.

*Code reference:
[R/surplus.R](https://fpcordeiro.github.io/choicer/R/surplus.R)*

------------------------------------------------------------------------

## 10. Goodness of Fit

### 10.1 McFadden Pseudo R-Squared

``` math
R^2 = 1 - \frac{\ell(\hat\theta)}{\ell_0},
\qquad
R^2_{\text{adj}} = 1 - \frac{\ell(\hat\theta) - K}{\ell_0}
```

where $`K`$ is the number of estimated parameters and $`\ell_0`$ is a
null-model log-likelihood. Two nulls are available:

**Equal shares** (default). Every alternative in individual $`i`$’s
choice set is equally likely:

``` math
\ell_0 = -\sum_{i=1}^{N} w_i \log\left(M_i + \mathbf{1}_{\text{outside}}\right)
```

where $`M_i`$ is the number of inside alternatives. This closed form is
exact for unbalanced choice sets and arbitrary weights.

**Market shares.** The maximized log-likelihood of a constants-only
(ASC-only) model. When every alternative is available in every choice
situation, the ASC-only model’s fitted probabilities equal the observed
market shares $`s_j`$, giving the closed form

``` math
\ell_0 = \sum_{j} N_j \log(s_j)
```

with $`N_j`$ the choice counts (the outside option enters as its own
category when present; never-chosen alternatives contribute
$`0 \cdot \log 0 = 0`$). This closed form requires identical choice-set
*composition* across individuals — equal set sizes are not sufficient —
and uniform weights; otherwise an ASC-only model must be refit
explicitly.

### 10.2 Hit Rate

The hit rate is the weighted share of choice situations in which the
observed choice has the highest predicted probability:

``` math
\mathrm{HR} = \frac{\sum_i w_i \, \mathbf{1}\{\hat{j}_i = j_i\}}{\sum_i w_i},
\qquad
\hat{j}_i = \arg\max_{j \in C_i} P_{ij}
```

With an outside option, the outside good competes for the maximum with
probability $`P_{i0} = 1 - \sum_j P_{ij}`$, and a predicted outside
choice is a hit when the outside good was in fact chosen ($`j_i = 0`$).

*Code reference:
[R/gof.R](https://fpcordeiro.github.io/choicer/R/gof.R)*

------------------------------------------------------------------------

## 11. Choice-Based Sampling and WESML Weighting

### 11.1 Endogenous Stratified (Choice-Based) Sampling

Under *random* (exogenous) sampling each choice situation is drawn
independently of its outcome. Under **choice-based** (endogenous
stratified) sampling the strata are defined by the **chosen
alternative** $`j_i`$, and situations are drawn with stratum-specific
frequencies — for example, deliberately oversampling individuals who
chose a rare alternative. Let $`Q(j)`$ be the population share of
alternative $`j`$ and $`H(j)`$ the sample share of choosers of $`j`$.

For the multinomial logit there is a well-known special case: with a
**full set of alternative-specific constants**, choice-based sampling
biases **only the constants**, and the slope coefficients $`\beta`$
remain consistent (Manski & McFadden 1981). Outside that case — most
importantly when the ASCs are not saturated, or when the constants
themselves are of interest (welfare, counterfactual entry) — naive
(unweighted) maximum likelihood is **inconsistent for the population
parameters**, and the weighting below is required.

### 11.2 The WESML Weighted Log-Likelihood

Manski and Lerman (1977) restore consistency by weighting each
situation’s contribution by the ratio of population to sample share of
its chosen alternative,

``` math
w_i = \frac{Q(j_i)}{H(j_i)},
\qquad
\ell^{W}(\theta) = \sum_{i=1}^{N} w_i \log P_{i j_i}(\theta).
```

The maximizer $`\hat\theta`$ is invariant to multiplying every $`w_i`$
by a common positive constant, so the weights may be normalized to mean
1 without changing the estimates.

### 11.3 Robust (Sandwich) Variance: the $`w^2`$ Meat

WESML is a weighted M-estimator solving
$`\sum_i w_i\, s_i(\hat\theta) = 0`$, where
$`s_i = \partial \log P_{i j_i}/\partial\theta`$ is the per-situation
score (\$$`3) and`$H_i = ^2 P\_{i j_i}/,^\$ its Hessian (\$\$4). Its
robust (Huber–White) asymptotic variance is the **sandwich**

``` math
V = A^{-1} B A^{-1}, \qquad
A = \sum_{i} w_i\,(-H_i), \qquad
B = \sum_{i} w_i^{2}\, s_i s_i^{\top}.
```

The weight enters the **bread** $`A`$ linearly, but the **meat** $`B`$
carries the weight **squared**: the contribution of situation $`i`$ to
the estimating equation is $`w_i s_i`$, so the variance of that
contribution involves $`(w_i s_i)(w_i s_i)^{\top}
= w_i^{2}\, s_i s_i^{\top}`$.

This is exactly why the two naive variances are wrong under weighting:

- the inverse-Hessian $`A^{-1}`$ assumes the information-matrix equality
  $`A = B`$, which **fails** once the $`w_i`$ are non-degenerate;
- the ordinary BHHH/OPG $`\big(\sum_i w_i\, s_i s_i^{\top}\big)^{-1}`$
  uses the weight to the *first* power, not the second.

A consistency check: rescaling all weights by a constant $`c`$ sends
$`A \to cA`$ and $`B \to c^{2} B`$, so
$`V \to (cA)^{-1}(c^{2}B)(cA)^{-1} = A^{-1}BA^{-1}`$ — the sandwich is
**invariant** to the weight scale, whereas $`A^{-1} \to c^{-1} A^{-1}`$
is not.

### 11.4 Scope: Robust vs. Design-Based Variance

The estimator implemented here is the **robust weighted-M-estimator
(Huber–White) variance**, whose meat
$`B = \sum_i w_i^{2}\, s_i s_i^{\top}`$ is uncentered. This is the
asymptotic variance under **variable-probability sampling** (each
population unit sampled independently with a stratum-dependent
probability). Under a **fixed-quota** design the design-based variance
centers the scores within strata and may be smaller; that
stratum-centered estimator is **out of scope** here. See Manski &
McFadden (1981) for the design-based treatment and Cosslett (1981) for
the asymptotically efficient conditional-maximum-likelihood alternative.

### 11.5 Implementation

The per-individual score accumulated inside `mnl_bhhh_parallel` is
**weight-free**, so the two summed matrices are obtained from the
existing exports evaluated at $`\hat\theta`$:

| Matrix | Call |
|----|----|
| Bread $`A = \sum_i w_i(-H_i)`$ | `mnl_loglik_hessian_parallel(..., weights = w)` |
| Meat $`B = \sum_i w_i^{2} s_i s_i^{\top}`$ | `mnl_bhhh_parallel(..., weights = w^2)` |

and combined in R as $`V = A^{-1} B A^{-1}`$. User-facing entry points:
[`wesml_weights()`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)
computes the Manski–Lerman weights from the population shares $`Q`$,
while
[`sample_by_choice()`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
draws a choice-based sample and attaches those weights;
`run_mnlogit(..., se_method = "sandwich")` estimates with the robust
variance;
[`wesml_vcov()`](https://fpcordeiro.github.io/choicer/reference/wesml_vcov.md)
returns it post hoc. Passing `run_mnlogit(..., weights_col = )`
collapses a row-level weight column to one weight per choice situation
(validated constant within `id`), and a provenance guard prevents a
WESML-labeled sample from being silently fit unweighted.

*Code reference:
[`mnl_bhhh_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_bhhh_parallel.md)
in
[src/mnlogit.cpp](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp);
`compute_sandwich_vcov()` and `.sandwich_combine()` in
[R/classes.R](https://fpcordeiro.github.io/choicer/R/classes.R);
[`wesml_weights()`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)
/
[`sample_by_choice()`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
in [R/sampling.R](https://fpcordeiro.github.io/choicer/R/sampling.R).*

------------------------------------------------------------------------

## 12. Implementation Details

### 11.1 Parameter Vector Structure

The full parameter vector $`\theta`$ is organized as:

| Block | Indices | Length | Description |
|----|----|----|----|
| $`\beta`$ | $`[0, K)`$ | $`K`$ | Coefficients for design matrix $`X`$ |
| $`\delta`$ | $`[K, K + J_{asc})`$ | $`J-1`$ or $`J`$ | Alternative-specific constants |

where $`J_{asc} = J - 1`$ without outside option (first ASC normalized
to 0), or $`J_{asc} = J`$ with outside option.

*Code reference:
[mnlogit.cpp:50-70](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L50-L70)*

### 11.2 Data Organization

- **Design matrix $`X`$**: Stacked matrix of dimension
  $`(\sum_i M_i) \times K`$, where $`M_i`$ is the number of (inside)
  alternatives for individual $`i`$
- **Alternative indices**: 1-based indexing in R, converted to 0-based
  in C++
- **Choice indices**: 1-based for inside alternatives, 0 for outside
  option when included
- **Prefix sums $`S`$**: Used for efficient indexing into the stacked
  data: $`S_i = \sum_{j < i} M_j`$

### 11.3 OpenMP Parallelization

The implementation parallelizes over individuals using OpenMP: - Each
thread maintains local accumulators for log-likelihood, gradient, and
Hessian - Thread results are combined using `#pragma omp critical`
sections - Dynamic scheduling is used for load balancing:
`#pragma omp for schedule(dynamic)`

*Code reference:
[mnlogit.cpp:85-183](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L85-L183)*

### 11.4 Negated Objectives

The C++ functions return: - $`-\ell(\theta)`$ (negated log-likelihood) -
$`-\nabla\ell(\theta)`$ (negated gradient) - $`-\nabla^2\ell(\theta)`$
(negated Hessian)

This is for compatibility with minimization routines (e.g., `nloptr`)
that expect a loss function to minimize rather than a likelihood to
maximize.

*Code reference:
[mnlogit.cpp:185-189](https://fpcordeiro.github.io/choicer/src/mnlogit.cpp#L185-L189)*

------------------------------------------------------------------------

## References

- McFadden, D. (1974). Conditional logit analysis of qualitative choice
  behavior. In P. Zarembka (Ed.), *Frontiers in Econometrics*
  (pp. 105-142). Academic Press.
- Train, K. E. (2009). *Discrete Choice Methods with Simulation* (2nd
  ed.). Cambridge University Press.
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in
  market equilibrium. *Econometrica*, 63(4), 841-890.
- Manski, C. F., & Lerman, S. R. (1977). The estimation of choice
  probabilities from choice based samples. *Econometrica*, 45(8),
  1977-1988.
- Manski, C. F., & McFadden, D. (1981). Alternative estimators and
  sample designs for discrete choice analysis. In C. F. Manski & D.
  McFadden (Eds.), *Structural Analysis of Discrete Data with
  Econometric Applications* (pp. 2-50). MIT Press.
- Cosslett, S. R. (1981). Maximum likelihood estimator for choice-based
  samples. *Econometrica*, 49(5), 1289-1316.
