# Mathematical Documentation: Mixed Logit Model

This document provides a detailed mathematical description of the mixed logit (random coefficients logit) model as implemented in `src/mxlogit.cpp`.

## Table of Contents

1. [Notation](#notation)
2. [Model Definition](#1-model-definition)
3. [Log-Likelihood Function](#2-log-likelihood-function)
4. [Gradient Computation](#3-gradient-computation)
5. [Hessian Computation](#4-hessian-computation)
6. [Implementation Details](#5-implementation-details)

---

## Notation

| Symbol | Description |
|--------|-------------|
| $i = 1, \ldots, N$ | Index for individuals (choice situations) |
| $j = 1, \ldots, J_i$ | Index for alternatives available to individual $i$ |
| $s = 1, \ldots, S$ | Index for simulation draws |
| $j_i$ | The alternative chosen by individual $i$ |
| $w_i$ | Weight for individual $i$ |
| $X_{ij}$ | Row vector of covariates with fixed coefficients ($1 \times K_x$) |
| $W_{ij}$ | Row vector of covariates with random coefficients ($1 \times K_w$) |
| $\beta$ | Fixed coefficient vector ($K_x \times 1$) |
| $\mu$ | Mean parameters for random coefficients ($K_w \times 1$) |
| $L$ | Lower-triangular Cholesky factor ($K_w \times K_w$) |
| $\Sigma = LL'$ | Covariance matrix of random coefficients |
| $\delta_j$ | Alternative-specific constant (ASC) for alternative $j$ |
| $\eta_i$ | Standard normal draw vector for individual $i$ ($K_w \times 1$) |
| $\gamma_i = L\eta_i$ | Random coefficient deviations for individual $i$ |

---

## 1. Model Definition

### 1.1 Utility Specification

The utility that individual $i$ derives from alternative $j$ under simulation draw $s$ is:

$$
V_{ij}^s = X_{ij}\beta + W_{ij}\left(\mu^* + \gamma_i^{s*}\right) + \delta_j
$$

where:
- $X_{ij}\beta$ captures the systematic utility from fixed coefficients
- $W_{ij}\mu^*$ captures the mean contribution of random coefficients
- $W_{ij}\gamma_i^{s*}$ captures the individual-specific random deviation
- $\delta_j$ is the alternative-specific constant (with one normalized to zero for identification)

### 1.2 Random Coefficient Structure

Random coefficients are generated through a Cholesky decomposition structure:

$$
\gamma_i^s = L \eta_i^s, \quad \text{where } \eta_i^s \sim N(0, I_{K_w})
$$

This implies:

$$
\gamma_i^s \sim N(0, \Sigma), \quad \text{where } \Sigma = LL'
$$

The Cholesky parameterization ensures $\Sigma$ is always positive semi-definite.

### 1.3 Distribution Transformations

The implementation supports two distributions for each random coefficient $k$:

**Normal Distribution** (`rc_dist[k] = 0`):
$$
\mu_k^* = \mu_k, \quad \gamma_{ik}^{s*} = \gamma_{ik}^s
$$

**Log-Normal Distribution** (`rc_dist[k] = 1`):
$$
\mu_k^* = \exp(\mu_k), \quad \gamma_{ik}^{s*} = \exp(\gamma_{ik}^s)
$$

> **Note:** This is a *shifted* log-normal parameterization: $\beta_k = \exp(\mu_k) + \exp(L_k \eta)$, which differs from the textbook parameterization $\beta_k = \exp(\mu_k + L_k \eta)$.

### 1.4 Cholesky Parameterization

For a full covariance matrix (`rc_correlation = TRUE`), the lower-triangular matrix $L$ is parameterized as:

$$
L = \begin{pmatrix}
\exp(\ell_{11}) & 0 & 0 & \cdots \\
\ell_{21} & \exp(\ell_{22}) & 0 & \cdots \\
\ell_{31} & \ell_{32} & \exp(\ell_{33}) & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}
$$

where:
- **Diagonal elements**: $L_{pp} = \exp(\ell_{pp})$ (ensures positivity)
- **Off-diagonal elements**: $L_{pq} = \ell_{pq}$ for $p > q$ (unconstrained)

For diagonal-only covariance (`rc_correlation = FALSE`):

$$
L = \text{diag}\left(\exp(\ell_{11}), \exp(\ell_{22}), \ldots, \exp(\ell_{K_w K_w})\right)
$$

The parameter vector for $L$ has length:
- Full: $K_w(K_w + 1)/2$
- Diagonal: $K_w$

---

## 2. Log-Likelihood Function

### 2.1 Choice Probability

For a given draw $s$, the probability that individual $i$ chooses alternative $j$ follows the logit formula:

$$
P_{ij}^s = \frac{\exp(V_{ij}^s)}{\sum_{k=1}^{J_i} \exp(V_{ik}^s)}
$$

If an outside option is included, it is normalized to $V_{i0}^s = 0$.

### 2.2 Simulated Probability

Since the random coefficients $\gamma_i$ are unobserved, we integrate over their distribution using simulation:

$$
\bar{P}_{ij_i} = \frac{1}{S} \sum_{s=1}^{S} P_{ij_i}^s
$$

where $j_i$ is the alternative chosen by individual $i$.

### 2.3 Log-Likelihood

The simulated log-likelihood is:

$$
\ell(\theta) = \sum_{i=1}^{N} w_i \log\left(\bar{P}_{ij_i}\right) = \sum_{i=1}^{N} w_i \log\left(\frac{1}{S} \sum_{s=1}^{S} P_{ij_i}^s\right)
$$

where $\theta = (\beta, \mu, \text{vec}(L), \delta)$ is the full parameter vector.

### 2.4 Log-Sum-Exp Trick for Numerical Stability

To avoid numerical overflow/underflow, the implementation uses:

1. **Probability computation**: Before computing $P_{ij}^s$, subtract $\max_k V_{ik}^s$:
   $$
   P_{ij}^s = \frac{\exp(V_{ij}^s - \max_k V_{ik}^s)}{\sum_k \exp(V_{ik}^s - \max_k V_{ik}^s)}
   $$

2. **Log-probability averaging**: Use log-sum-exp to compute $\log(\sum_s P_{ij_i}^s)$:
   $$
   \text{logSumExp}(a, b) = \max(a,b) + \log\left(\exp(a - \max(a,b)) + \exp(b - \max(a,b))\right)
   $$

   Applied iteratively across draws $s$.

---

## 3. Gradient Computation

### 3.1 General Formula

Define for draw $s$ the log-probability of the chosen alternative:

$$
\log P_{ij_i}^s = V_{ij_i}^s - \log\left(\sum_k \exp(V_{ik}^s)\right)
$$

The gradient of this with respect to utility is:

$$
\frac{\partial \log P_{ij_i}^s}{\partial V_{ia}^s} = \mathbf{1}_{a = j_i} - P_{ia}^s
$$

### 3.2 Per-Draw Gradient

Define $z_{ij}^s = \frac{\partial V_{ij}^s}{\partial \theta}$ as the gradient of utility with respect to parameters. The per-draw gradient is:

$$
g_i^s = \frac{\partial \log P_{ij_i}^s}{\partial \theta} = \sum_{j=1}^{J_i} \left(\mathbf{1}_{j = j_i} - P_{ij}^s\right) z_{ij}^s
$$

### 3.3 Gradient of Simulated Log-Likelihood

Using the identity for the derivative of a logarithm of a sum:

$$
\frac{\partial}{\partial \theta} \log\left(\frac{1}{S}\sum_s P_{ij_i}^s\right) = \frac{\sum_s P_{ij_i}^s \cdot g_i^s}{\sum_s P_{ij_i}^s}
$$

The full gradient is:

$$
\nabla_\theta \ell = \sum_{i=1}^{N} w_i \cdot \frac{\sum_s P_{ij_i}^s \cdot g_i^s}{\sum_s P_{ij_i}^s}
$$

### 3.4 Utility Derivatives by Parameter Block

#### 3.4.1 Fixed Coefficients ($\beta$)

$$
\frac{\partial V_{ij}^s}{\partial \beta} = X_{ij}^T
$$

*Code reference: [mxlogit.cpp:303-310](../src/mxlogit.cpp#L303-L310)*

#### 3.4.2 Random Coefficient Means ($\mu$)

For Normal distribution:
$$
\frac{\partial V_{ij}^s}{\partial \mu_p} = W_{ijp}
$$

For Log-Normal distribution:
$$
\frac{\partial V_{ij}^s}{\partial \mu_p} = W_{ijp} \cdot \exp(\mu_p)
$$

*Code reference: [mxlogit.cpp:333-339](../src/mxlogit.cpp#L333-L339)*

#### 3.4.3 Cholesky Parameters ($L$)

Let $\ell_{pq}$ denote the parameter for $L_{pq}$ (where $p \geq q$).

**Step 1**: Derivative of $L_{pq}$ with respect to $\ell_{pq}$:
$$
\frac{\partial L_{pq}}{\partial \ell_{pq}} = \begin{cases}
L_{pp} = \exp(\ell_{pp}) & \text{if } p = q \text{ (diagonal)} \\
1 & \text{if } p > q \text{ (off-diagonal)}
\end{cases}
$$

**Step 2**: Derivative of $\gamma_p^s$ with respect to $\ell_{pq}$:
$$
\frac{\partial \gamma_p^s}{\partial \ell_{pq}} = \frac{\partial L_{pq}}{\partial \ell_{pq}} \cdot \eta_q^s
$$

Note: $\gamma_p^s = \sum_{r \leq p} L_{pr} \eta_r^s$, so only the term with $r = q$ contributes.

**Step 3**: Derivative of transformed $\gamma_p^{s*}$:

For Normal distribution:
$$
\frac{\partial \gamma_p^{s*}}{\partial \gamma_p^s} = 1
$$

For Log-Normal distribution:
$$
\frac{\partial \gamma_p^{s*}}{\partial \gamma_p^s} = \exp(\gamma_p^s) = \gamma_p^{s*}
$$

**Step 4**: Chain rule for utility derivative:
$$
\frac{\partial V_{ij}^s}{\partial \ell_{pq}} = W_{ijp} \cdot \frac{\partial \gamma_p^{s*}}{\partial \gamma_p^s} \cdot \frac{\partial \gamma_p^s}{\partial \ell_{pq}}
$$

*Code reference: [mxlogit.cpp:342-371](../src/mxlogit.cpp#L342-L371)*

#### 3.4.4 Alternative-Specific Constants ($\delta$)

$$
\frac{\partial V_{ij}^s}{\partial \delta_a} = \mathbf{1}_{j = a}
$$

where $\delta_1 = 0$ (normalized) if no outside option, or all $\delta_j$ for $j \geq 1$ are free if outside option is included.

*Code reference: [mxlogit.cpp:312-323](../src/mxlogit.cpp#L312-L323)*

---

## 4. Hessian Computation

### 4.1 General Formula

Define:
- $g_i = \frac{\partial \log \bar{P}_i}{\partial \theta}$ — the gradient for individual $i$
- $g_i^s = \frac{\partial \log P_{ij_i}^s}{\partial \theta}$ — the per-draw gradient
- $H_i^s = \frac{\partial^2 \log P_{ij_i}^s}{\partial \theta \partial \theta'}$ — the per-draw Hessian

The Hessian of the simulated log-probability for individual $i$ is:

$$
H_i = \frac{\partial^2 \log \bar{P}_i}{\partial \theta \partial \theta'} = \frac{\sum_s P_{ij_i}^s \left(g_i^s (g_i^s)' + H_i^s\right)}{\sum_s P_{ij_i}^s} - g_i g_i'
$$

The full Hessian is:

$$
\nabla^2_\theta \ell = \sum_{i=1}^{N} w_i \cdot H_i
$$

*Code reference: [mxlogit.cpp:829-840](../src/mxlogit.cpp#L829-L840)*

### 4.2 Per-Draw Hessian

The per-draw Hessian $H_i^s$ has the form:

$$
H_i^s = -\sum_j P_{ij}^s z_{ij}^s (z_{ij}^s)' + \left(\sum_j P_{ij}^s z_{ij}^s\right)\left(\sum_j P_{ij}^s z_{ij}^s\right)' + \sum_j \left(\mathbf{1}_{j=j_i} - P_{ij}^s\right) H_{V,ij}^s
$$

where $H_{V,ij}^s = \frac{\partial^2 V_{ij}^s}{\partial \theta \partial \theta'}$ is the Hessian of utility.

In the code, this is computed as:
- `sum_Pzz` $= \sum_j P_{ij}^s z_{ij}^s (z_{ij}^s)'$
- `sum_Pz` $= \sum_j P_{ij}^s z_{ij}^s$
- `sum_diff_H_V` $= \sum_j (\mathbf{1}_{j=j_i} - P_{ij}^s) H_{V,ij}^s$

Then: $H_i^s = -\text{sum\_Pzz} + \text{sum\_Pz} \cdot \text{sum\_Pz}' + \text{sum\_diff\_H\_V}$

*Code reference: [mxlogit.cpp:826](../src/mxlogit.cpp#L826)*

### 4.3 Second Derivatives of Utility

Most second derivatives are zero since utility is linear in most parameters. Non-zero second derivatives arise from:

#### 4.3.1 Mean Parameters ($\mu$) — Log-Normal Only

For log-normal distribution:
$$
\frac{\partial^2 V_{ij}^s}{\partial \mu_p^2} = W_{ijp} \cdot \exp(\mu_p)
$$

*Code reference: [mxlogit.cpp:750-752](../src/mxlogit.cpp#L750-L752)*

#### 4.3.2 Cholesky Diagonal Parameters ($\ell_{pp}$)

Since $L_{pp} = \exp(\ell_{pp})$:
$$
\frac{\partial^2 L_{pp}}{\partial \ell_{pp}^2} = \exp(\ell_{pp}) = L_{pp}
$$

For normal distribution:
$$
\frac{\partial^2 V_{ij}^s}{\partial \ell_{pp}^2} = W_{ijp} \cdot L_{pp} \cdot \eta_p^s
$$

For log-normal distribution (combining chain rule):
$$
\frac{\partial^2 V_{ij}^s}{\partial \ell_{pp}^2} = W_{ijp} \left[ \gamma_p^{s*} \left(\frac{\partial \gamma_p^s}{\partial \ell_{pp}}\right)^2 + \frac{\partial \gamma_p^{s*}}{\partial \gamma_p^s} \cdot \frac{\partial^2 \gamma_p^s}{\partial \ell_{pp}^2} \right]
$$

where $\frac{\partial^2 \gamma_p^s}{\partial \ell_{pp}^2} = L_{pp} \cdot \eta_p^s$.

*Code reference: [mxlogit.cpp:774-778](../src/mxlogit.cpp#L774-L778)*

#### 4.3.3 Cholesky Cross-Derivatives ($\ell_{pq}$, $\ell_{pr}$)

For off-diagonal elements with the same row $p$ but different columns $q \neq r$:

This contributes through the log-normal transformation (when applicable):
$$
\frac{\partial^2 V_{ij}^s}{\partial \ell_{pq} \partial \ell_{pr}} = W_{ijp} \cdot \gamma_p^{s*} \cdot \frac{\partial \gamma_p^s}{\partial \ell_{pq}} \cdot \frac{\partial \gamma_p^s}{\partial \ell_{pr}}
$$

*Code reference: [mxlogit.cpp:782-790](../src/mxlogit.cpp#L782-L790)*

---

## 5. Implementation Details

### 5.1 Parameter Vector Structure

The full parameter vector $\theta$ is organized as:

| Block | Indices | Length | Description |
|-------|---------|--------|-------------|
| $\beta$ | $[0, K_x)$ | $K_x$ | Fixed coefficients |
| $\mu$ | $[K_x, K_x + K_w)$ | $K_w$ (if `rc_mean=TRUE`) | Random coefficient means |
| $L$ | next block | $K_w(K_w+1)/2$ or $K_w$ | Cholesky parameters |
| $\delta$ | remaining | $J-1$ or $J$ | ASCs |

### 5.2 Numerical Hessian (Alternative)

The implementation also provides a numerical Hessian via central finite differences:

$$
H_{jk} \approx \frac{g_k(\theta + \epsilon_j e_j) - g_k(\theta - \epsilon_j e_j)}{2\epsilon_j}
$$

where $\epsilon_j = \epsilon \cdot \max(|\theta_j|, 1)$ for adaptive step size.

*Code reference: [mxlogit.cpp:419-482](../src/mxlogit.cpp#L419-L482)*

### 5.3 Delta Method for Variance Parameters

When reporting results, the Cholesky parameters $\ell$ are transformed to variance matrix elements $\Sigma = LL'$.

The Jacobian $J = \frac{\partial \text{vech}(\Sigma)}{\partial \ell}$ is computed analytically in `jacobian_vech_Sigma()`.

For element $\Sigma_{ij} = \sum_k L_{ik} L_{jk}$:
$$
\frac{\partial \Sigma_{ij}}{\partial \ell_{pq}} = \frac{\partial L_{pq}}{\partial \ell_{pq}} \left( L_{jq} \mathbf{1}_{i=p} + L_{iq} \mathbf{1}_{j=p} \right)
$$

The variance-covariance matrix of $\text{vech}(\Sigma)$ is then:
$$
V_\Sigma = J \cdot V_\ell \cdot J'
$$

*Code reference: [mxlogit.cpp:503-541](../src/mxlogit.cpp#L503-L541)*

### 5.4 OpenMP Parallelization

The implementation parallelizes over individuals using OpenMP:
- Each thread maintains local accumulators for log-likelihood and gradient/Hessian
- Thread results are combined using `#pragma omp critical` sections
- Dynamic scheduling is used for load balancing: `#pragma omp for schedule(dynamic)`

### 5.5 Practical Notes

1. **Negated Objectives**: The C++ functions return $-\ell(\theta)$ and $-\nabla\ell(\theta)$ (and $-\nabla^2\ell(\theta)$ for the Hessian) for compatibility with minimization routines that expect a loss function.

2. **ASC Identification**: When `include_outside_option = FALSE`, the first inside alternative's ASC ($\delta_1$) is fixed to zero for identification. When `include_outside_option = TRUE`, all inside ASCs are free parameters (the outside option with $V=0$ serves as the reference).

---

## References

- Train, K. E. (2009). *Discrete Choice Methods with Simulation*. Cambridge University Press.
- McFadden, D., & Train, K. (2000). Mixed MNL models for discrete response. *Journal of Applied Econometrics*, 15(5), 447-470.
