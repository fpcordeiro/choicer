# Gibbs sampler for the hierarchical Bayesian multinomial logit model

Runs the adaptive random-walk Metropolis-within-Gibbs sampler for the
hierarchical (random-coefficients, panel) multinomial logit with a
BLP-style alternative-level random effect: inside utilities \\U\_{ijt} =
x\_{ijt}'\gamma_i + \delta_j + EV1\\ against an implicit outside option
with systematic utility 0, \\\beta_i \sim N(b, W)\\ (\\\gamma\_{ik} =
\beta\_{ik}\\ or \\\exp(\beta\_{ik})\\ per `rc_dist`), and \\\delta_j =
z_j'\theta + \xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\.

## Usage

``` r
hmnl_gibbs(
  X,
  Z,
  M,
  choice_pos,
  include_outside_option,
  alt_of_row,
  Ti,
  rc_dist,
  beta_pooled,
  delta_init,
  theta_init,
  b_bar,
  A,
  nu,
  V,
  theta_bar,
  A_theta,
  sd_prior,
  R,
  burn,
  thin,
  seed,
  keep_beta_i,
  s_init,
  accept_target,
  trace = 0L
)
```

## Arguments

- X:

  total_rows x K_struct structural design matrix (inside rows only, no
  ASC columns), rows sorted by (person, task, alternative).

- Z:

  J x P alternative-level mean-function design (intercept first).

- M:

  Integer vector: inside alternatives per choice situation.

- choice_pos:

  Integer vector: 1-based within-task position of the chosen row; 0 =
  outside option chosen.

- include_outside_option:

  Must be `TRUE` (the implicit outside good anchors the location of
  delta; a no-outside mode is roadmapped).

- alt_of_row:

  Integer vector: 1-based alternative code per row of X.

- Ti:

  Integer vector: choice situations per respondent.

- rc_dist:

  Integer vector (length K_struct): 0 = normal coordinate, 1 =
  log-normal (enters utility as `exp(beta_ik)`).

- beta_pooled:

  Pooled MNL MLE on the chain scale (log scale for log-normal
  coordinates); centers the H_i proposal information.

- delta_init:

  Initial delta (length J).

- theta_init:

  Initial theta (length P).

- b_bar:

  K vector, prior mean of b.

- A:

  K x K prior precision matrix of b.

- nu:

  Inverse-Wishart prior degrees of freedom for W (\>= K).

- V:

  K x K inverse-Wishart prior scale matrix for W.

- theta_bar:

  P vector, prior mean of theta.

- A_theta:

  P x P prior precision matrix of theta.

- sd_prior:

  List with elements `half_cauchy` (logical), `s_d` (half-Cauchy scale),
  `c0`, `d0` (IG fallback).

- R:

  Total number of Gibbs iterations.

- burn:

  Number of initial iterations discarded (0 \<= burn \< R);
  proposal-scale adaptation happens during burn-in only.

- thin:

  Keep every thin-th post-burn-in draw.

- seed:

  Master RNG seed (non-negative; all streams derive from it).

- keep_beta_i:

  0 = no beta_i output, 1 = online means/SDs, 2 = means/SDs plus the
  full (K, N, R_keep) draw cube.

- s_init:

  Initial per-respondent proposal scale.

- accept_target:

  Robbins-Monro acceptance target for the beta_i updates (the delta_j
  target is fixed at 0.44).

- trace:

  Print progress every `trace` iterations (0 = silent).

## Value

List with `bdraw` (R_keep x K), `wdraw` (R_keep x K(K+1)/2, lower
triangle of W in row-major order), `deltadraw` (R_keep x J), `thetadraw`
(R_keep x P), `sigma_d2draw`, `loglik_trace`, acceptance rates and final
proposal scales (`accept_rate_beta`, `accept_rate_delta`, `s_final`,
`s_delta_final`), posterior summaries `beta_i_mean` / `beta_i_sd` (K x
N, `NULL` when `keep_beta_i = 0`), `beta_i_draws` (K x N x R_keep cube
when `keep_beta_i = 2`), `delta_mean` / `delta_sd` / `xi_mean` / `xi_sd`
(J x 1), and `R_keep`.

## Details

The per-respondent \\\beta_i\\ updates are parallelized with OpenMP; the
\\\delta_j\\ updates run as a strictly serial sweep (their conditionals
are coupled through the shared softmax denominators). Each (iteration,
unit) pair uses its own RNG stream, so draws are bitwise reproducible
independent of the number of threads (see
[`set_num_threads()`](https://fpcordeiro.github.io/choicer/reference/set_num_threads.md)).
This is the low-level engine behind
[`run_hmnlogit`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md),
which handles initialization and post-processing.

## Examples

``` r
# \donttest{
sim <- simulate_hmnl_data(N = 20, T = 2, J = 3, seed = 42)
d <- prepare_hmnl_data(sim$data, "task", "alt", "choice",
                       c("x1", "x2"), person_col = "pid")
out <- hmnl_gibbs(d$X, d$Z, d$M, d$choice_pos, TRUE, d$alt_of_row, d$Ti,
  rc_dist = d$rc_dist, beta_pooled = rep(0, d$K_struct),
  delta_init = rep(0, d$J), theta_init = rep(0, d$P),
  b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
  nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
  theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
  sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
  R = 300, burn = 100, thin = 1, seed = 7, keep_beta_i = 1,
  s_init = 2.38 / sqrt(d$K_struct), accept_target = 0.234)
colMeans(out$bdraw)
#> [1]  0.7338366 -1.4694685
# }
```
