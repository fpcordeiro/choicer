# Mathematical Documentation: Nested Logit Model

This document provides a detailed mathematical description of the nested logit (NL) model as implemented in `src/nestlogit.cpp`.

## Table of Contents

1. [Notation](#notation)
2. [Model Definition](#1-model-definition)
3. [Log-Likelihood Function](#2-log-likelihood-function)
4. [Gradient Computation](#3-gradient-computation)
5. [Hessian Computation](#4-hessian-computation)
6. [Implementation Details](#5-implementation-details)
7. [Post-Estimation Quantities](#6-post-estimation-quantities)

---

## Notation

| Symbol | Description |
|--------|-------------|
| $i = 1, \ldots, N$ | Index for individuals (choice situations) |
| $j = 1, \ldots, J_i$ | Index for inside alternatives available to individual $i$ |
| $j_i$ | The alternative chosen by individual $i$ |
| $w_i$ | Weight for individual $i$ |
| $X_{ij}$ | Row vector of covariates for individual $i$ and alternative $j$ ($1 \times K$) |
| $\beta$ | Coefficient vector ($K \times 1$) |
| $\delta_j$ | Alternative-specific constant (ASC) for alternative $j$ |
| $V_{ij}$ | Systematic (deterministic) utility |
| $B_k$ | Set of alternatives belonging to nest $k$ |
| $k(j)$ | The nest containing alternative $j$ |
| $\lambda_k$ | Inclusive value (dissimilarity) coefficient for nest $k$ |
| $I_k$ | Inclusive value for nest $k$ |
| $P(j \mid k)$ | Probability of choosing $j$ conditional on choosing nest $k$ |
| $P_k$ | Marginal probability of choosing nest $k$ |
| $P_{ij}$ | Probability that individual $i$ chooses alternative $j$ |

---

## 1. Model Definition

### 1.1 Utility Specification

The utility that individual $i$ derives from alternative $j$ is:

$$
U_{ij} = V_{ij} + \varepsilon_{ij}
$$

where the systematic utility is:

$$
V_{ij} = X_{ij}\beta + \delta_j
$$

The error term $\varepsilon_{ij}$ follows a Generalized Extreme Value (GEV) distribution that induces within-nest correlation.

### 1.2 Nest Structure

The $J$ inside alternatives are partitioned into $K$ mutually exclusive and exhaustive nests $B_1, \ldots, B_K$. An outside option (alternative 0) with $V_{i0} = 0$ may optionally be included, and it is always assigned its own singleton nest with $\lambda_0 = 1$.

### 1.3 Inclusive Value Coefficients

Each nest $k$ has a dissimilarity parameter $\lambda_k > 0$ that governs within-nest substitution:

- $\lambda_k = 1$ corresponds to independent (MNL-like) alternatives within the nest.
- $0 < \lambda_k < 1$ indicates positive within-nest correlation (the typical case).
- $\lambda_k > 1$ is theoretically permissible but implies negative correlation.

**Singleton nests** (nests with exactly one alternative) have $\lambda_k$ fixed to 1 and are not estimated. Only non-singleton nests contribute $\lambda$ parameters to $\theta$.

*Code reference: [nestlogit.cpp:61-104](../src/nestlogit.cpp#L61-L104)*

### 1.4 Alternative-Specific Constants (ASCs)

The ASC identification follows the same convention as MNL:

- **Without outside option** (`include_outside_option = FALSE`): The first inside alternative's ASC is fixed to zero ($\delta_1 = 0$). The parameter vector contains $J-1$ free ASC parameters.
- **With outside option** (`include_outside_option = TRUE`): The outside option serves as the reference with $V_{i0} = 0$. All $J$ inside alternatives have free ASC parameters.

*Code reference: [nestlogit.cpp:107-135](../src/nestlogit.cpp#L107-L135)*

---

## 2. Log-Likelihood Function

### 2.1 Probability Decomposition

The nested logit choice probability decomposes as:

$$
P_{ij} = P(j \mid k(j)) \cdot P_{k(j)}
$$

The two components are derived from the GEV generating function.

### 2.2 Inclusive Value

The inclusive value for nest $k$ is the log-sum of scaled utilities within the nest:

$$
I_k = \sum_{j \in B_k} \exp\!\left(\frac{V_{ij}}{\lambda_k}\right)
$$

$$
\log I_k = \log\!\left(\sum_{j \in B_k} \exp\!\left(\frac{V_{ij}}{\lambda_k}\right)\right)
$$

### 2.3 Conditional Choice Probability

The probability of choosing alternative $j$ given that nest $k(j)$ is chosen:

$$
P(j \mid k) = \frac{\exp(V_{ij}/\lambda_k)}{I_k}
$$

In log form:

$$
\log P(j \mid k) = \frac{V_{ij}}{\lambda_k} - \log I_k
$$

### 2.4 Nest (Marginal) Probability

The probability of choosing nest $k$:

$$
P_k = \frac{\exp(\lambda_k \log I_k)}{\sum_{l} \exp(\lambda_l \log I_l)}
$$

If an outside option is included, its contribution enters the denominator as $\exp(\lambda_0 \log I_0) = \exp(1 \cdot \log(\exp(0/1))) = \exp(0) = 1$.

In log form:

$$
\log P_k = \lambda_k \log I_k - \log\!\left(\sum_l \exp(\lambda_l \log I_l)\right)
$$

### 2.5 Joint Choice Probability

Taking logs of the joint probability:

$$
\log P_{ij} = \log P(j \mid k(j)) + \log P_{k(j)}
$$

$$
\log P_{ij} = \frac{V_{ij}}{\lambda_{k(j)}} - \log I_{k(j)} + \lambda_{k(j)} \log I_{k(j)} - \log D_i
$$

where $D_i = \sum_l \exp(\lambda_l \log I_l)$ (plus 1 if an outside option is present).

*Code reference: [nestlogit.cpp:235-260](../src/nestlogit.cpp#L235-L260)*

### 2.6 Log-Likelihood

The log-likelihood is the weighted sum of log-probabilities of the observed choices:

$$
\ell(\theta) = \sum_{i=1}^{N} w_i \log P_{ij_i}
$$

where $\theta = (\beta, \lambda, \delta)$ and $j_i$ is the alternative chosen by individual $i$.

*Code reference: [nestlogit.cpp:263-266](../src/nestlogit.cpp#L263-L266)*

### 2.7 Log-Sum-Exp Trick for Numerical Stability

Two log-sum-exp stabilizations are applied.

**Within-nest inclusive value.** To compute $\log I_k$ without overflow, the implementation subtracts the maximum scaled utility in each nest:

$$
\log I_k = V_{\max,k} + \log\!\left(\sum_{j \in B_k} \exp\!\left(\frac{V_{ij}}{\lambda_k} - V_{\max,k}\right)\right)
$$

where $V_{\max,k} = \max_{j \in B_k} V_{ij}/\lambda_k$.

**Across-nest denominator.** The denominator $D_i = \sum_l \exp(\lambda_l \log I_l)$ is similarly stabilized by subtracting $\max_l \lambda_l \log I_l$ before exponentiation.

*Code reference: [nestlogit.cpp:190-230](../src/nestlogit.cpp#L190-L230)*

---

## 3. Gradient Computation

### 3.1 Gradient with Respect to $V_{ia}$ (Utility Scores)

To derive the gradients with respect to $\beta$ and $\delta$, we first compute the derivative of the log-probability of the chosen alternative with respect to the systematic utility $V_{ia}$ of each alternative $a$.

Let $k_i$ denote the chosen nest (i.e., $k(j_i)$). The log-probability of the chosen alternative is:

$$
\log P_{ij_i} = \log P(j_i \mid k_i) + \log P_{k_i}
$$

**Derivative of the conditional term.** Since $\log P(j \mid k) = V_{ij}/\lambda_k - \log I_k$ and $\frac{\partial \log I_k}{\partial V_{ia}} = \frac{1}{\lambda_k} \mathbf{1}_{a \in B_k} P(a \mid k)$:

$$
\frac{\partial \log P(j_i \mid k_i)}{\partial V_{ia}} = \frac{1}{\lambda_{k_i}}\left(\mathbf{1}_{a = j_i} - \mathbf{1}_{a \in B_{k_i}} P(a \mid k_i)\right)
$$

**Derivative of the nest term.** Using $\frac{\partial \log D_i}{\partial V_{ia}} = P_{k(a)} P(a \mid k(a)) = P_{ia}$:

$$
\frac{\partial \log P_{k_i}}{\partial V_{ia}} = \mathbf{1}_{a \in B_{k_i}} P(a \mid k_i) - P_{ia}
$$

**Combined.** Adding both terms:

$$
\frac{\partial \log P_{ij_i}}{\partial V_{ia}} = \frac{\mathbf{1}_{a = j_i}}{\lambda_{k_i}} + \left(1 - \frac{1}{\lambda_{k_i}}\right)\mathbf{1}_{a \in B_{k_i}} P(a \mid k_i) - P_{ia}
$$

This is implemented as the vector `grad_vec` of length $J_i$:

$$
\text{grad\_vec}[a] = -P_{ia} + \left(1 - \frac{1}{\lambda_{k_i}}\right) \mathbf{1}_{a \in B_{k_i}} P(a \mid k_i) + \frac{1}{\lambda_{k_i}} \mathbf{1}_{a = j_i}
$$

*Code reference: [nestlogit.cpp:284-295](../src/nestlogit.cpp#L284-L295)*

### 3.2 Gradient with Respect to $\beta$

Since $\frac{\partial V_{ij}}{\partial \beta} = X_{ij}^T$, the contribution from individual $i$ is:

$$
\frac{\partial \ell_i}{\partial \beta} = w_i \sum_{j=1}^{J_i} \text{grad\_vec}[j] \cdot X_{ij}^T = w_i \, X_i^T \, \text{grad\_vec}
$$

The full gradient is accumulated as a single BLAS matrix-vector product, identically to MNL.

*Code reference: [nestlogit.cpp:297-298](../src/nestlogit.cpp#L297-L298)*

### 3.3 Gradient with Respect to $\delta$

Since $\frac{\partial V_{ij}}{\partial \delta_a} = \mathbf{1}_{j = a}$, the gradient entries are the corresponding elements of `grad_vec`, scattered to the appropriate positions in $\theta$:

$$
\frac{\partial \ell_i}{\partial \delta_a} = w_i \cdot \text{grad\_vec}[a]
$$

The same free-ASC convention as MNL applies (first inside alternative's ASC excluded when no outside option).

*Code reference: [nestlogit.cpp:300-312](../src/nestlogit.cpp#L300-L312)*

### 3.4 Gradient with Respect to $\lambda_k$

Define the within-nest probability-weighted mean utility:

$$
\bar{V}_k = \sum_{j \in B_k} P(j \mid k) \cdot V_{ij}
$$

**For the chosen nest ($k = k_i$).** The log-probability decomposes as:

$$
\log P_{ij_i} = \underbrace{\frac{V_{ij_i}}{\lambda_{k_i}} - \log I_{k_i}}_{\log P(j_i \mid k_i)} + \underbrace{\lambda_{k_i} \log I_{k_i} - \log D_i}_{\log P_{k_i}}
$$

Using $\frac{\partial \log I_k}{\partial \lambda_k} = -\bar{V}_k / \lambda_k^2$:

$$
\frac{\partial \log P(j_i \mid k_i)}{\partial \lambda_{k_i}} = \frac{\bar{V}_{k_i} - V_{ij_i}}{\lambda_{k_i}^2}
$$

$$
\frac{\partial \log P_{k_i}}{\partial \lambda_{k_i}} = (1 - P_{k_i})\!\left(\log I_{k_i} - \frac{\bar{V}_{k_i}}{\lambda_{k_i}}\right)
$$

Combining:

$$
\frac{\partial \ell_i}{\partial \lambda_{k_i}} = w_i \left[(1 - P_{k_i})\!\left(\log I_{k_i} - \frac{\bar{V}_{k_i}}{\lambda_{k_i}}\right) + \frac{\bar{V}_{k_i} - V_{ij_i}}{\lambda_{k_i}^2}\right]
$$

**For a non-chosen nest ($k \neq k_i$).** The conditional probability $P(j_i \mid k_i)$ is independent of $\lambda_k$, so only the nest probability contributes:

$$
\frac{\partial \log P_{k_i}}{\partial \lambda_k} = -P_k \!\left(\log I_k - \frac{\bar{V}_k}{\lambda_k}\right)
$$

$$
\frac{\partial \ell_i}{\partial \lambda_k} = -w_i \, P_k \!\left(\log I_k - \frac{\bar{V}_k}{\lambda_k}\right)
$$

The common factor $\log I_k - \bar{V}_k / \lambda_k$ is labeled `term_in_brackets` in the code.

*Code reference: [nestlogit.cpp:314-342](../src/nestlogit.cpp#L314-L342)*

---

## 4. Hessian Computation

The nested logit Hessian is computed numerically via central finite differences on the gradient:

$$
\frac{\partial^2 \ell}{\partial \theta_i \, \partial \theta_j} \approx \frac{(\nabla_\theta \ell)_j\big|_{\theta + \epsilon_i e_i} - (\nabla_\theta \ell)_j\big|_{\theta - \epsilon_i e_i}}{2\epsilon_i}
$$

where $e_i$ is the $i$-th standard basis vector and the step size is scaled to the parameter magnitude:

$$
\epsilon_i = \epsilon \cdot \max(|\theta_i|, 1)
$$

with default $\epsilon = 10^{-6}$.

This approach avoids the considerably more complex analytical derivation required for a closed-form nested logit Hessian. The resulting matrix is used to compute the variance-covariance matrix of the estimates.

*Code reference: [nestlogit.cpp:396-458](../src/nestlogit.cpp#L396-L458)*

---

## 5. Implementation Details

### 5.1 Parameter Vector Structure

The full parameter vector $\theta$ is organized as:

| Block | Indices | Length | Description |
|-------|---------|--------|-------------|
| $\beta$ | $[0, K)$ | $K$ | Coefficients for design matrix $X$ |
| $\lambda$ | $[K, K + K_\lambda)$ | $K_\lambda$ | Non-singleton nest dissimilarity parameters |
| $\delta$ | $[K + K_\lambda, K + K_\lambda + J_{\text{asc}})$ | $J-1$ or $J$ | Alternative-specific constants |

where $K_\lambda$ is the number of non-singleton nests and $J_{\text{asc}} = J - 1$ without outside option or $J$ with outside option.

*Code reference: [nestlogit.cpp:57-135](../src/nestlogit.cpp#L57-L135)*

### 5.2 Singleton Nest Handling

At startup, the implementation identifies which nests are singletons (contain exactly one alternative). These nests:
- Have their $\lambda$ fixed to 1 (no contribution to the upper-level logit denominator beyond their single utility term).
- Are excluded from the estimated parameter vector.
- Are tracked via `nest_k_to_theta_idx`, an integer vector of length $K$ where entry $k$ holds the corresponding index in $\theta$ for non-singleton nests and $-1$ for singletons.

This ensures that the NL model degenerates gracefully to MNL on a per-nest basis when nests are singletons.

*Code reference: [nestlogit.cpp:61-104](../src/nestlogit.cpp#L61-L104)*

### 5.3 Data Organization

- **Design matrix $X$**: Stacked matrix of dimension $(\sum_i M_i) \times K$, where $M_i$ is the number of inside alternatives for individual $i$.
- **`alt_idx`**: 1-based indices mapping rows of $X$ to global alternative IDs.
- **`nest_idx`**: 1-based nest assignment of length $J$ (number of unique inside alternatives), not length $\sum M_i$. Accessed via `nest_idx.elem(alt_idx0_i)` to handle varying choice sets.
- **`choice_idx`**: 1-based index of the chosen inside alternative; 0 for outside option.
- **Prefix sums $S$**: Used for efficient slicing into the stacked data: $S_i = \sum_{j < i} M_j$.

*Code reference: [nestlogit.cpp:137-146](../src/nestlogit.cpp#L137-L146)*

### 5.4 Base Utility Pre-computation

Prior to the parallel loop, all systematic utilities $V_{ij} = X_{ij}\beta + \delta_j$ are computed in a single pass:

```
base_util = X * beta                 // single BLAS dgemm
if use_asc: base_util += delta[alt_idx0]
```

This avoids redundant matrix-vector products inside the per-individual loop.

*Code reference: [nestlogit.cpp:144-146](../src/nestlogit.cpp#L144-L146)*

### 5.5 OpenMP Parallelization

The implementation parallelizes over individuals using OpenMP:
- Each thread maintains local accumulators for log-likelihood and gradient.
- `nest_k_to_theta_idx` is copied thread-privately to avoid false sharing.
- Thread results are combined using `#pragma omp critical` sections.
- Dynamic scheduling balances load across individuals with varying choice set sizes.

*Code reference: [nestlogit.cpp:152-353](../src/nestlogit.cpp#L152-L353)*

### 5.6 Negated Objectives

The C++ function returns:
- $-\ell(\theta)$ (negated log-likelihood)
- $-\nabla\ell(\theta)$ (negated gradient)

for compatibility with minimization routines (e.g., `nloptr`) that minimize a loss function rather than maximize a likelihood.

*Code reference: [nestlogit.cpp:356-359](../src/nestlogit.cpp#L356-L359)*

---

## 6. Post-Estimation Quantities

### 6.1 Choice Probabilities (predict)

The joint probability is the decomposition from §2.1:

$$
P_{ij} = P(j \mid k(j)) \cdot P_{k(j)}
$$

with $P(j \mid k)$ and $P_k$ given in §§2.3–2.4. `nl_predict` returns individual-level probabilities $P_{ij}$ and systematic utilities $V_{ij}$ for every observation row. `nl_predict_shares` aggregates over individuals to obtain weighted market shares:

$$
\hat{s}_j = \frac{\sum_{i=1}^{N} w_i \, P_{ij}}{\sum_{i=1}^{N} w_i}
$$

*Code reference: [nestlogit.cpp:521-708](../src/nestlogit.cpp#L521-L708)*

---

### 6.2 Elasticities

#### 6.2.1 Semi-elasticity with respect to utility

Let $r = k(j)$ be the nest of the responding alternative $j$ and $s = k(a)$ be the nest of the perturbed alternative $a$. From §3.1, the derivative of $\log P_{ij}$ with respect to the systematic utility of any alternative $a$ is:

$$
\frac{\partial \log P_{ij}}{\partial V_{ia}}
= \frac{\mathbf{1}_{a=j}}{\lambda_r}
+ \left(1 - \frac{1}{\lambda_r}\right)\mathbf{1}_{s=r} \, P(a \mid r)
- P_{ia}
$$

This is the general expression; the gradient formula in §3.1 applies it only to the *chosen* alternative, whereas here we evaluate it for all pairs $(j, a)$ to construct the full elasticity matrix.

#### 6.2.2 Elasticity of $P_{ij}$ with respect to covariate $k$ of alternative $a$

Since $\partial V_{ia} / \partial x_{ia,k} = \beta_k$:

$$
E_{ij,ia}^{k}
= \frac{\partial P_{ij}}{\partial x_{ia,k}} \cdot \frac{x_{ia,k}}{P_{ij}}
= \beta_k \, x_{ia,k} \cdot \frac{\partial \log P_{ij}}{\partial V_{ia}}
$$

Substituting the three cases:

| Case | Condition | Semi-elasticity $\partial \log P_{ij}/\partial V_{ia}$ |
|------|-----------|--------------------------------------------------------|
| Own | $a = j$ | $\dfrac{1}{\lambda_r} + \left(1 - \dfrac{1}{\lambda_r}\right) P(j \mid r) - P_{ij}$ |
| Cross, same nest | $s = r,\; a \neq j$ | $\left(1 - \dfrac{1}{\lambda_r}\right) P(a \mid r) - P_{ia}$ |
| Cross, different nest | $s \neq r$ | $-P_{ia}$ |

**Reduction to MNL.** When $\lambda_k = 1$ for all nests the within-nest term $(1 - 1/\lambda_r)$ vanishes and the conditional probability drops out, giving $\partial \log P_{ij}/\partial V_{ia} = \mathbf{1}_{a=j} - P_{ia}$, which is exactly the MNL semi-elasticity (§3.1 of the MNL document).

**Outside option.** When an outside option is present it occupies its own singleton nest with $\lambda_0 = 1$ and $x_{i0,k} = 0$. Its row (responding = outside option) has semi-elasticity $-P_{ia}$ for every inside alternative $a$ (different nest, §3.1); its column (perturbed = outside option) is identically zero because $\beta_k \cdot 0 = 0$.

The implementation accumulates $w_i \cdot E_{ij,ia}^k$ across individuals and divides by $\sum_i w_i$ to produce weighted average elasticities.

*Code reference: [nestlogit.cpp:749-892](../src/nestlogit.cpp#L749-L892)*

---

### 6.3 Diversion Ratios

#### 6.3.1 Definition

The diversion ratio from alternative $k$ to alternative $j$ is:

$$
DR(k \to j) = -\frac{\partial Q_k / \partial V_{ij}}{\partial Q_j / \partial V_{ij}}
$$

where $Q_m = \sum_i w_i P_{im}$ is aggregate demand and the perturbation is through the systematic utility of the *losing* alternative $j$.

#### 6.3.2 Partial derivatives of demand

Using $\partial P_{im}/\partial V_{ij} = P_{im} \cdot (\partial \log P_{im}/\partial V_{ij})$ and the semi-elasticity from §6.2.1, the three cases with $r = k(j)$:

$$
\frac{\partial P_{im}}{\partial V_{ij}} =
\begin{cases}
P_{ij}\!\left[\dfrac{1}{\lambda_r} + \left(1 - \dfrac{1}{\lambda_r}\right) P(j \mid r) - P_{ij}\right] & m = j \\[6pt]
P_{im}\!\left[\left(1 - \dfrac{1}{\lambda_r}\right) P(j \mid r) - P_{ij}\right] & k(m) = r,\; m \neq j \\[6pt]
-P_{im} \, P_{ij} & k(m) \neq r
\end{cases}
$$

#### 6.3.3 Covariate-independence of diversion ratios

For a concrete attribute $k$ the perturbation is $\partial V_{ij} = \beta_k \, \partial x_{ij,k}$. Both numerator and denominator of $DR(k \to j)$ are proportional to the *same* factor $\beta_k \, x_{ij,k}$ (the same attribute of the same losing alternative $j$), so it cancels:

$$
DR(k \to j)
= -\frac{\beta_k \, x_{ij,k} \sum_i w_i \,(\partial P_{ik}/\partial V_{ij})}
         {\beta_k \, x_{ij,k} \sum_i w_i \,(\partial P_{ij}/\partial V_{ij})}
= -\frac{\sum_i w_i \,(\partial P_{ik}/\partial V_{ij})}
         {\sum_i w_i \,(\partial P_{ij}/\partial V_{ij})}
$$

The NL diversion ratio is therefore independent of which covariate is perturbed — a property it shares with MNL (see §6.2 of the MNL document) but which does **not** hold for mixed logit, where $\beta_k$ is individual-specific.

#### 6.3.4 Population aggregation and column-sum property

Define:

$$
\text{num}(k, j) = \sum_{i=1}^{N} w_i \left(-\frac{\partial P_{ik}}{\partial V_{ij}}\right), \qquad
\text{den}(j) = \sum_{i=1}^{N} w_i \frac{\partial P_{ij}}{\partial V_{ij}}
$$

so that $DR(k \to j) = \text{num}(k,j) \,/\, \text{den}(j)$ for $k \neq j$.

**Columns sum to 1.** The fundamental identity $\sum_{m} \partial P_{im}/\partial V_{ij} = 0$ (probabilities sum to 1 for every individual and every perturbation direction, including the outside option when present) gives:

$$
\sum_{k \neq j} \frac{\partial P_{ik}}{\partial V_{ij}} = -\frac{\partial P_{ij}}{\partial V_{ij}}
\implies
\sum_{k \neq j} \text{num}(k, j) = \text{den}(j)
\implies
\sum_{k \neq j} DR(k \to j) = 1
$$

When an outside option is present it acts as both a source (it can lose demand) and a destination (it can receive diverted demand), so it appears in both the row index $k$ and the column index $j$ of the $J_{\text{total}} \times J_{\text{total}}$ diversion matrix, and the column sums still equal 1 over all alternatives including the outside option.

#### 6.3.5 Reduction to MNL at $\lambda = 1$

When all $\lambda_k = 1$ the cross-same-nest case collapses to the cross-different-nest case ($\partial P_{im}/\partial V_{ij} = -P_{im} P_{ij}$ for $m \neq j$) and the own case gives $\partial P_{ij}/\partial V_{ij} = P_{ij}(1 - P_{ij})$. Substituting:

$$
DR(k \to j)\big|_{\lambda=1}
= \frac{\sum_i w_i \, P_{ij} \, P_{ik}}{\sum_i w_i \, P_{ij} (1 - P_{ij})}
$$

which is the MNL diversion ratio (§6.3 of the MNL document).

*Code reference: [nestlogit.cpp:927-1073](../src/nestlogit.cpp#L927-L1073)*

---

### 6.4 BLP Share Inversion

#### 6.4.1 Problem statement

Given observed market shares $s_j^*$, find the ASC vector $\delta$ such that the NL model reproduces those shares:

$$
\hat{s}_j(\delta) = s_j^* \quad \forall j
$$

#### 6.4.2 Damped contraction mapping

Berry, Levinsohn, and Pakes (1995) establish that the fixed-point iteration

$$
\delta^{(t+1)} = \delta^{(t)} + \kappa \bigl(\log s^* - \log \hat{s}(\delta^{(t)})\bigr)
$$

converges to the unique solution. Setting $\kappa = 1$ reproduces the original BLP95 update. Predicted shares $\hat{s}(\delta^{(t)})$ are computed via the NL formula: for each alternative $j$,

$$
\hat{s}_j(\delta^{(t)}) = \frac{\sum_{i=1}^{N} w_i \, P_{ij}(\delta^{(t)})}{\sum_{i=1}^{N} w_i}
$$

where $P_{ij}(\delta^{(t)}) = P(j \mid k(j)) \cdot P_{k(j)}$ is evaluated at the current $\delta^{(t)}$ with $\beta$ and $\lambda$ held fixed.

#### 6.4.3 Damping for nested models

For strongly nested data (small $\lambda_k$) the plain BLP update ($\kappa = 1$) can overshoot. Setting $\kappa < 1$ improves contraction. A practical heuristic is:

$$
\kappa \approx \min_k \lambda_k
$$

since the curvature of the NL log-share map is sharpest in nests with small dissimilarity parameters. The function argument is called `damping`.

#### 6.4.4 Baseline re-normalization

After convergence the returned $\delta$ is re-centered so that the baseline is identified: $\delta \leftarrow \delta - \delta_0$ (subtract the outside-option entry, or equivalently $\delta_1$ when there is no outside option). This matches the ASC identification convention from §1.4.

**Note on the `lambda` argument.** `nl_blp_contraction` expects a *full* length-$n_{\text{nests}}$ vector where singleton nests have their entry set to 1. This differs from the estimated parameter block $\lambda$ in $\theta$, which contains only non-singleton entries. The caller is responsible for expanding the estimated lambdas into the full vector before invoking the contraction.

*Code reference: [nestlogit.cpp:1115-1186](../src/nestlogit.cpp#L1115-L1186)*

### 6.5 Willingness to Pay

Utility is linear in covariates, so willingness to pay is the same marginal rate of substitution as in the MNL model: for attribute $k$ and price coefficient $\alpha = \beta_p$,

$$
\mathrm{WTP}_k = -\frac{\beta_k}{\alpha},
\qquad
\widehat{\mathrm{SE}} = \sqrt{\nabla g^T \, \hat{V}_{(k,p)} \, \nabla g},
\quad
\nabla g = \begin{pmatrix} -1/\alpha \\ \beta_k/\alpha^2 \end{pmatrix}
$$

with $\hat{V}_{(k,p)}$ the $2 \times 2$ block of the coefficient covariance. The nest parameters $\lambda$ rescale within-nest substitution but do not enter the ratio of marginal utilities, so they affect WTP only through the estimated covariance of $(\beta_k, \alpha)$. ASCs may also be used as the numerator ($-\delta_j/\alpha$). See §8 of the MNL document for the full derivation and caveats.

*Code reference: [R/wtp.R](../R/wtp.R)*

### 6.6 Consumer Surplus and the Nested Logsum

For the nested logit, the expected maximum utility (logsum) over individual $i$'s choice set replaces the flat MNL log-sum with the two-level structure of §2:

$$
\mathrm{logsum}_i = \log D_i = \log\left( \sum_{l} \exp\big(\lambda_l \log I_l\big) \; [+\, 1] \right),
\qquad
\log I_l = \log \sum_{j \in B_l \cap C_i} \exp\!\left(\frac{V_{ij}}{\lambda_l}\right)
$$

where the $+1 = \exp(0)$ term appears when an outside option is present and singleton nests have $\lambda_l = 1$ (§5.2). Both levels are computed with the max-subtraction trick (§2.7).

Two useful identities:

- **Reduction to MNL.** With all $\lambda_l = 1$, $\lambda_l \log I_l = \log \sum_{j \in B_l} \exp(V_{ij})$ and the nested logsum collapses to the MNL logsum $\log \sum_j \exp(V_{ij})$.
- **Outside-option link.** Since $P_{i0} = 1 / D_i$ (§2.5), the logsum equals $-\log P_{i0}$ — a direct consistency check against the probability kernel.

Expected consumer surplus follows Train (2009, Ch. 3) exactly as in the MNL case:

$$
\mathbb{E}[CS_i] = \frac{\mathrm{logsum}_i}{-\alpha}
$$

with the same caveats: no income effects, and CS levels inherit the ASC normalization, so only differences across scenarios (counterfactual `newdata` vs. baseline) are meaningful. The implementation reports point estimates for NL (the delta-method SE of the mean CS is currently provided for MNL only; Krinsky–Robb resampling of the coefficients is the practical interval method here).

*Code reference: [R/surplus.R](../R/surplus.R)*

### 6.7 Goodness of Fit

The measures are model-agnostic and identical to §10 of the MNL document, computed from the NL log-likelihood and predicted probabilities:

- **McFadden pseudo R²**: $R^2 = 1 - \ell(\hat\theta)/\ell_0$ and $R^2_{\text{adj}} = 1 - (\ell(\hat\theta) - K)/\ell_0$, with $K$ counting all estimated parameters including the $\lambda$ block.
- **Equal-shares null** (default, exact for unbalanced sets and weights): $\ell_0 = -\sum_i w_i \log(M_i + \mathbf{1}_{\text{outside}})$.
- **Market-shares null** (constants-only closed form): $\ell_0 = \sum_j N_j \log(s_j)$; requires identical choice-set composition across individuals and uniform weights. The constants-only model is itself an MNL — at $\beta = 0$ and free ASCs the nest structure is undetermined — so the same closed form applies.
- **Hit rate**: weighted share of choice situations whose observed choice has the highest predicted NL probability; with an outside option the outside good competes with $P_{i0} = 1 - \sum_j P_{ij}$.

*Code reference: [R/gof.R](../R/gof.R)*

---

## References

- McFadden, D. (1978). Modelling the choice of residential location. In A. Karlqvist et al. (Eds.), *Spatial Interaction Theory and Planning Models* (pp. 75-96). North-Holland.
- Ben-Akiva, M., & Lerman, S. R. (1985). *Discrete Choice Analysis: Theory and Application to Travel Demand*. MIT Press.
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica*, 63(4), 841–890.
- Train, K. E. (2009). *Discrete Choice Methods with Simulation* (2nd ed.). Cambridge University Press. Chapter 4.
