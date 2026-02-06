# Mathematical Documentation: Multinomial Logit Model

This document provides a detailed mathematical description of the multinomial logit (MNL) model as implemented in `src/mnlogit.cpp`.

## Table of Contents

1. [Notation](#notation)
2. [Model Definition](#1-model-definition)
3. [Log-Likelihood Function](#2-log-likelihood-function)
4. [Gradient Computation](#3-gradient-computation)
5. [Hessian Computation](#4-hessian-computation)
6. [Elasticity Computation](#5-elasticity-computation)
7. [Diversion Ratio Computation](#6-diversion-ratio-computation)
8. [BLP Contraction Mapping](#7-blp-contraction-mapping)
9. [Implementation Details](#8-implementation-details)

---

## Notation

| Symbol | Description |
|--------|-------------|
| $i = 1, \ldots, N$ | Index for individuals (choice situations) |
| $j = 1, \ldots, J_i$ | Index for alternatives available to individual $i$ |
| $j_i$ | The alternative chosen by individual $i$ |
| $w_i$ | Weight for individual $i$ |
| $X_{ij}$ | Row vector of covariates for individual $i$ and alternative $j$ ($1 \times K$) |
| $\beta$ | Coefficient vector ($K \times 1$) |
| $\delta_j$ | Alternative-specific constant (ASC) for alternative $j$ |
| $V_{ij}$ | Systematic (deterministic) utility |
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

and $\varepsilon_{ij}$ is an i.i.d. Type I Extreme Value (Gumbel) error term with location 0 and scale 1.

### 1.2 Alternative-Specific Constants (ASCs)

The ASC parameters $\delta_j$ capture the average effect of unobserved factors for each alternative. For identification, one ASC must be normalized:

- **Without outside option** (`include_outside_option = FALSE`): The first inside alternative's ASC is fixed to zero ($\delta_1 = 0$). The parameter vector contains $J-1$ free ASC parameters.

- **With outside option** (`include_outside_option = TRUE`): The outside option (alternative 0) has utility $V_{i0} = 0$, serving as the reference. All $J$ inside alternatives have free ASC parameters.

*Code reference: [mnlogit.cpp:39-55](../src/mnlogit.cpp#L39-L55)*

---

## 2. Log-Likelihood Function

### 2.1 Choice Probability

Under the Type I Extreme Value distribution assumption, the probability that individual $i$ chooses alternative $j$ follows the multinomial logit (softmax) formula:

$$
P_{ij} = \frac{\exp(V_{ij})}{\sum_{k=1}^{J_i} \exp(V_{ik})}
$$

If an outside option is included, the denominator includes $\exp(V_{i0}) = \exp(0) = 1$.

*Code reference: [mnlogit.cpp:102-105](../src/mnlogit.cpp#L102-L105)*

### 2.2 Log-Likelihood

The log-likelihood is the weighted sum of the log-probabilities of the observed choices:

$$
\ell(\theta) = \sum_{i=1}^{N} w_i \log P_{ij_i}
$$

where $\theta = (\beta, \delta)$ is the full parameter vector and $j_i$ is the alternative chosen by individual $i$.

Expanding the log-probability:

$$
\log P_{ij_i} = V_{ij_i} - \log\left(\sum_{k=1}^{J_i} \exp(V_{ik})\right)
$$

*Code reference: [mnlogit.cpp:117-125](../src/mnlogit.cpp#L117-L125)*

### 2.3 Log-Sum-Exp Trick for Numerical Stability

To prevent numerical overflow when computing the denominator, the implementation subtracts the maximum utility before exponentiation:

$$
\sum_k \exp(V_{ik}) = \exp(V_{\max}) \sum_k \exp(V_{ik} - V_{\max})
$$

where $V_{\max} = \max_k V_{ik}$.

This gives:

$$
\log\left(\sum_k \exp(V_{ik})\right) = V_{\max} + \log\left(\sum_k \exp(V_{ik} - V_{\max})\right)
$$

The log-probability becomes:

$$
\log P_{ij_i} = (V_{ij_i} - V_{\max}) - \log\left(\sum_k \exp(V_{ik} - V_{\max})\right)
$$

*Code reference: [mnlogit.cpp:103-104](../src/mnlogit.cpp#L103-L104)*

---

## 3. Gradient Computation

### 3.1 General Formula

The gradient of the log-probability of the chosen alternative with respect to utility is:

$$
\frac{\partial \log P_{ij_i}}{\partial V_{ia}} = \mathbf{1}_{a = j_i} - P_{ia}
$$

where $\mathbf{1}_{a = j_i}$ is the indicator function (1 if $a$ is the chosen alternative, 0 otherwise).

### 3.2 Gradient with Respect to $\beta$

Since $\frac{\partial V_{ij}}{\partial \beta} = X_{ij}^T$, the contribution to the gradient from individual $i$ is:

$$
\frac{\partial \ell_i}{\partial \beta} = w_i \sum_{j=1}^{J_i} \left(\mathbf{1}_{j = j_i} - P_{ij}\right) X_{ij}^T
$$

This simplifies to:

$$
\frac{\partial \ell_i}{\partial \beta} = w_i \left( X_{ij_i}^T - \sum_{j=1}^{J_i} P_{ij} X_{ij}^T \right)
$$

The full gradient is:

$$
\nabla_\beta \ell = \sum_{i=1}^{N} w_i \left( X_{ij_i}^T - \sum_{j=1}^{J_i} P_{ij} X_{ij}^T \right)
$$

*Code reference: [mnlogit.cpp:132-139](../src/mnlogit.cpp#L132-L139)*

### 3.3 Gradient with Respect to $\delta$

Since $\frac{\partial V_{ij}}{\partial \delta_a} = \mathbf{1}_{j = a}$, the contribution to the gradient from individual $i$ for alternative $a$ is:

$$
\frac{\partial \ell_i}{\partial \delta_a} = w_i \left(\mathbf{1}_{j_i = a} - P_{ia}\right)
$$

The gradient is computed only for the free ASC parameters:
- Without outside option: $\delta_1 = 0$ (fixed), so we compute for $a = 2, \ldots, J$
- With outside option: we compute for all inside alternatives $a = 1, \ldots, J$

*Code reference: [mnlogit.cpp:141-152](../src/mnlogit.cpp#L141-L152)*

---

## 4. Hessian Computation

### 4.1 Analytical Hessian

The implementation provides an analytical Hessian for efficient computation of standard errors.

Define the "score vector" for alternative $a$ as $Z_a$, which concatenates:
- $X_{ia}^T$ for the $\beta$ components
- Indicator $\mathbf{1}_{a = j}$ for the $\delta_j$ components

The Hessian of the log-likelihood for individual $i$ is:

$$
H_i = \frac{\partial^2 \ell_i}{\partial \theta \partial \theta'} = -w_i \left[ \sum_{j=1}^{J_i} P_{ij} Z_j Z_j^T - \left(\sum_{j=1}^{J_i} P_{ij} Z_j\right)\left(\sum_{j=1}^{J_i} P_{ij} Z_j\right)^T \right]
$$

This can be written more compactly as:

$$
H_i = -w_i \left[ \mathbb{E}_P[Z Z^T] - \mathbb{E}_P[Z] \mathbb{E}_P[Z]^T \right] = -w_i \cdot \text{Var}_P(Z)
$$

where the expectation is taken with respect to the choice probabilities $P_{ij}$.

The full Hessian is:

$$
\nabla^2_\theta \ell = \sum_{i=1}^{N} H_i
$$

In the code:
- `sum_P_Z` = $\sum_j P_{ij} Z_j$
- `sum_P_Z_Zt` = $\sum_j P_{ij} Z_j Z_j^T$
- `H_i` = $w_i \cdot (\text{sum\_P\_Z} \cdot \text{sum\_P\_Z}^T - \text{sum\_P\_Z\_Zt})$

*Code reference: [mnlogit.cpp:695-736](../src/mnlogit.cpp#L695-L736)*

### 4.2 Numerical Hessian

The implementation also provides a numerical Hessian via central finite differences for verification:

$$
H_{jk} \approx \frac{g_k(\theta + \epsilon_j e_j) - g_k(\theta - \epsilon_j e_j)}{2\epsilon_j}
$$

where:
- $e_j$ is the unit vector in direction $j$
- $\epsilon_j = \epsilon \cdot \max(|\theta_j|, 1)$ provides adaptive step size
- $g$ is the gradient function

*Code reference: [mnlogit.cpp:203-246](../src/mnlogit.cpp#L203-L246)*

---

## 5. Elasticity Computation

### 5.1 Elasticity Definitions

The elasticity measures the percentage change in choice probability with respect to a percentage change in an attribute. For attribute $k$ with coefficient $\beta_k$:

**Own-Elasticity** (elasticity of $P_{ij}$ with respect to $x_{ijk}$):

$$
E_{jj}^k = \frac{\partial P_{ij}}{\partial x_{ijk}} \cdot \frac{x_{ijk}}{P_{ij}} = \beta_k \cdot x_{ijk} \cdot (1 - P_{ij})
$$

**Cross-Elasticity** (elasticity of $P_{ij}$ with respect to $x_{imk}$ where $m \neq j$):

$$
E_{jm}^k = \frac{\partial P_{ij}}{\partial x_{imk}} \cdot \frac{x_{imk}}{P_{ij}} = -\beta_k \cdot x_{imk} \cdot P_{im}
$$

### 5.2 Derivation

Starting from the choice probability:

$$
P_{ij} = \frac{\exp(V_{ij})}{\sum_k \exp(V_{ik})}
$$

Taking the derivative with respect to $x_{imk}$:

$$
\frac{\partial P_{ij}}{\partial x_{imk}} = P_{ij} \left( \mathbf{1}_{j=m} - P_{im} \right) \beta_k
$$

For own-elasticity ($j = m$):

$$
\frac{\partial P_{ij}}{\partial x_{ijk}} = P_{ij} (1 - P_{ij}) \beta_k
$$

Therefore:

$$
E_{jj}^k = \beta_k \cdot x_{ijk} \cdot (1 - P_{ij})
$$

For cross-elasticity ($j \neq m$):

$$
\frac{\partial P_{ij}}{\partial x_{imk}} = -P_{ij} P_{im} \beta_k
$$

Therefore:

$$
E_{jm}^k = -\beta_k \cdot x_{imk} \cdot P_{im}
$$

### 5.3 Aggregate Elasticities

The implementation computes weighted average elasticities across all individuals:

$$
\bar{E}_{jm}^k = \frac{\sum_{i=1}^{N} w_i \cdot E_{ijm}^k}{\sum_{i=1}^{N} w_i}
$$

where $E_{ijm}^k$ is the elasticity for individual $i$.

*Code reference: [mnlogit.cpp:892-912](../src/mnlogit.cpp#L892-L912)*

---

## 6. Diversion Ratio Computation

### 6.1 Definition

The diversion ratio from alternative $j$ to alternative $k$ measures the fraction of demand lost by $j$ that is captured by $k$ when $j$ becomes less attractive (e.g., due to a price increase). It is a key metric in antitrust analysis and merger simulation.

$$
DR(j \to k) = -\frac{\partial Q_k / \partial p_j}{\partial Q_j / \partial p_j}
$$

where $Q_j = \sum_i w_i P_{ij}$ is the (weighted) aggregate demand for alternative $j$ and $p_j$ is the price of $j$.

### 6.2 Derivation for MNL

From the MNL choice probability derivatives (Section 5.2), the demand derivatives with respect to price $p_j$ (with coefficient $\beta_p$) are:

$$
\frac{\partial Q_k}{\partial p_j} = -\beta_p \sum_{i=1}^{N} w_i \, P_{ij} \, P_{ik} \quad (k \neq j)
$$

$$
\frac{\partial Q_j}{\partial p_j} = \beta_p \sum_{i=1}^{N} w_i \, P_{ij} \, (1 - P_{ij})
$$

Substituting into the diversion ratio formula:

$$
DR(j \to k) = -\frac{-\beta_p \sum_i w_i \, P_{ij} \, P_{ik}}{\beta_p \sum_i w_i \, P_{ij} \, (1 - P_{ij})} = \frac{\sum_i w_i \, P_{ij} \, P_{ik}}{\sum_i w_i \, P_{ij} \, (1 - P_{ij})}
$$

The price coefficient $\beta_p$ cancels, so the diversion ratio does not depend on which covariate we differentiate with respect to. This is a consequence of the IIA property.

### 6.3 Aggregate Diversion Ratio Matrix

The implementation computes the full $J \times J$ diversion ratio matrix $D$ where:

$$
D_{kj} = DR(j \to k) = \frac{\sum_{i=1}^{N} w_i \, P_{ij} \, P_{ik}}{\sum_{i=1}^{N} w_i \, P_{ij} \, (1 - P_{ij})} \quad (k \neq j)
$$

$$
D_{jj} = 0
$$

### 6.4 Properties

**Column sums equal one.** For a given alternative $j$, the off-diagonal entries in column $j$ sum to 1, meaning all diverted demand is accounted for:

$$
\sum_{k \neq j} DR(j \to k) = \frac{\sum_i w_i \, P_{ij} \sum_{k \neq j} P_{ik}}{\sum_i w_i \, P_{ij} (1 - P_{ij})} = \frac{\sum_i w_i \, P_{ij} (1 - P_{ij})}{\sum_i w_i \, P_{ij} (1 - P_{ij})} = 1
$$

where we used $\sum_{k \neq j} P_{ik} = 1 - P_{ij}$.

**IIA proportionality.** When all individuals face the same choice set and the same covariates (so $P_{ij}$ does not vary across $i$), the diversion ratio simplifies to:

$$
DR(j \to k) = \frac{P_j \, P_k}{P_j (1 - P_j)} = \frac{P_k}{1 - P_j} = \frac{s_k}{1 - s_j}
$$

where $s_j$ is the market share of alternative $j$. This means diversion is proportional to the receiving alternative's market share, independent of any characteristics of the losing alternative $j$ other than its own share. This is a well-known limitation of MNL due to IIA.

### 6.5 Implementation

The function accumulates two quantities in parallel across individuals:

1. **Numerator matrix**: $N_{kj} = \sum_i w_i \, P_{ij} \, P_{ik}$ for each pair $(k, j)$
2. **Denominator vector**: $d_j = \sum_i w_i \, P_{ij} \, (1 - P_{ij})$ for each $j$

The final matrix is computed as $D_{kj} = N_{kj} / d_j$ for $k \neq j$, with $D_{jj} = 0$.

*Code reference: [mnlogit.cpp:938-1069](../src/mnlogit.cpp#L938-L1069)*

---

## 7. BLP Contraction Mapping

### 7.1 Problem Statement

Given observed market shares $s_j$, find the ASC parameters $\delta_j$ such that the model-predicted shares match the observed shares:

$$
\hat{s}_j(\delta) = s_j \quad \forall j
$$

where $\hat{s}_j(\delta)$ is the predicted market share for alternative $j$.

### 7.2 Contraction Mapping Algorithm

Berry, Levinsohn, and Pakes (1995) show that the following iteration converges to the solution:

$$
\delta_j^{(t+1)} = \delta_j^{(t)} + \log(s_j) - \log(\hat{s}_j^{(t)})
$$

This is equivalent to:

$$
\delta^{(t+1)} = \delta^{(t)} + \log(s) - \log(\hat{s}^{(t)})
$$

### 7.3 Implementation Details

1. **Initialization**: Start with initial guess $\delta^{(0)}$
2. **Prediction**: Compute predicted shares $\hat{s}(\delta^{(t)})$ using `mnl_predict_shares_internal`
3. **Update**: Apply the contraction: $\delta^{(t+1)} = \delta^{(t)} + \log(s) - \log(\hat{s}^{(t)})$
4. **Convergence**: Check if $\max_j |\delta_j^{(t+1)} - \delta_j^{(t)}| < \text{tol}$
5. **Normalization**: Subtract $\delta_1$ from all ASCs to maintain identification

*Code reference: [mnlogit.cpp:566-586](../src/mnlogit.cpp#L566-L586)*

---

## 8. Implementation Details

### 8.1 Parameter Vector Structure

The full parameter vector $\theta$ is organized as:

| Block | Indices | Length | Description |
|-------|---------|--------|-------------|
| $\beta$ | $[0, K)$ | $K$ | Coefficients for design matrix $X$ |
| $\delta$ | $[K, K + J_{asc})$ | $J-1$ or $J$ | Alternative-specific constants |

where $J_{asc} = J - 1$ without outside option (first ASC normalized to 0), or $J_{asc} = J$ with outside option.

*Code reference: [mnlogit.cpp:30-55](../src/mnlogit.cpp#L30-L55)*

### 8.2 Data Organization

- **Design matrix $X$**: Stacked matrix of dimension $(\sum_i M_i) \times K$, where $M_i$ is the number of (inside) alternatives for individual $i$
- **Alternative indices**: 1-based indexing in R, converted to 0-based in C++
- **Choice indices**: 1-based for inside alternatives, 0 for outside option when included
- **Prefix sums $S$**: Used for efficient indexing into the stacked data: $S_i = \sum_{j < i} M_j$

### 8.3 OpenMP Parallelization

The implementation parallelizes over individuals using OpenMP:
- Each thread maintains local accumulators for log-likelihood, gradient, and Hessian
- Thread results are combined using `#pragma omp critical` sections
- Dynamic scheduling is used for load balancing: `#pragma omp for schedule(dynamic)`

*Code reference: [mnlogit.cpp:66-164](../src/mnlogit.cpp#L66-L164)*

### 8.4 Negated Objectives

The C++ functions return:
- $-\ell(\theta)$ (negated log-likelihood)
- $-\nabla\ell(\theta)$ (negated gradient)
- $-\nabla^2\ell(\theta)$ (negated Hessian)

This is for compatibility with minimization routines (e.g., `nloptr`) that expect a loss function to minimize rather than a likelihood to maximize.

*Code reference: [mnlogit.cpp:167-170](../src/mnlogit.cpp#L167-L170)*

---

## References

- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), *Frontiers in Econometrics* (pp. 105-142). Academic Press.
- Train, K. E. (2009). *Discrete Choice Methods with Simulation* (2nd ed.). Cambridge University Press.
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica*, 63(4), 841-890.
