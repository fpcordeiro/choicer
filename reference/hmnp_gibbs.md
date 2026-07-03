# Gibbs sampler for the hierarchical Bayesian multinomial probit model

Runs the fully conjugate Albert-Chib Gibbs sampler for the hierarchical
multinomial probit with iid \\N(0, \sigma^2)\\ utility shocks in
un-differenced utility space: inside utilities \\U\_{ijt} =
x\_{ijt}'\beta_i + \delta_j + \epsilon\\ against a stochastic implicit
outside option \\U\_{iot} = \epsilon\\, \\\beta_i \sim N(b, W)\\, and
\\\delta_j = z_j'\theta + \xi_j\\, \\\xi_j \sim N(0, \sigma_d^2)\\. The
chain runs on the non-identified parameterization (free \\\sigma^2\\);
identified quantities are obtained by normalizing each draw by the
matching power of \\\sigma\\ (handled by
[`run_hmnprobit`](https://fpcordeiro.github.io/choicer/reference/run_hmnprobit.md)).

## Usage

``` r
hmnp_gibbs(
  X,
  Z,
  M,
  choice_pos,
  include_outside_option,
  alt_of_row,
  Ti,
  delta_init,
  theta_init,
  b_bar,
  A,
  nu,
  V,
  theta_bar,
  A_theta,
  sd_prior,
  a0,
  s0,
  R,
  burn,
  thin,
  seed,
  keep_beta_i,
  trace = 0L
)
```

## Arguments

- X:

  total_rows x K_struct structural design matrix (inside rows only),
  rows sorted by (person, task, alternative).

- Z:

  J x P alternative-level mean-function design (intercept first).

- M:

  Integer vector: inside alternatives per choice situation.

- choice_pos:

  Integer vector: 1-based within-task position of the chosen row; 0 =
  outside option chosen.

- include_outside_option:

  Must be `TRUE` (the outside good anchors the location of delta; a
  no-outside mode is roadmapped).

- alt_of_row:

  Integer vector: 1-based alternative code per row of X.

- Ti:

  Integer vector: choice situations per respondent.

- delta_init:

  Initial delta (length J), raw scale.

- theta_init:

  Initial theta (length P), raw scale.

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

  List with elements `half_cauchy` (logical), `s_d`, `c0`, `d0` — the
  sigma_d prior (see
  [`run_hmnlogit`](https://fpcordeiro.github.io/choicer/reference/run_hmnlogit.md)).

- a0, s0:

  Inverse-gamma prior shape/scale for the (non-identified) shock
  variance \\\sigma^2\\.

- R:

  Total number of Gibbs iterations.

- burn:

  Number of initial iterations discarded (0 \<= burn \< R).

- thin:

  Keep every thin-th post-burn-in draw.

- seed:

  Master RNG seed (non-negative; all streams derive from it).

- keep_beta_i:

  0 = no beta_i output, 1 = online means/SDs (identified scale), 2 =
  means/SDs plus the full (K, N, R_keep) cube of per-draw-normalized
  \\\beta_i / \sigma\\ draws.

- trace:

  Print progress every `trace` iterations (0 = silent).

## Value

List with RAW draw matrices `bdraw`, `wdraw` (lower triangle,
row-major), `deltadraw`, `thetadraw`, `sigma_d2draw`, `sigma2draw`,
identified-scale summaries `beta_i_mean` / `beta_i_sd` / `beta_i_draws`
/ `delta_mean` / `delta_sd` / `xi_mean` / `xi_sd`, and `R_keep`.

## Details

The latent sweep and the \\\beta_i\\ draws are parallelized with OpenMP
across respondents; the \\\delta_j\\ draws are parallelized across
alternatives (conditionally independent given the augmented utilities —
unlike the HMNL, whose delta sweep must be serial). Each (iteration,
unit) pair uses its own RNG stream, so draws are bitwise reproducible
independent of the number of threads.

## Examples

``` r
# \donttest{
sim <- simulate_hmnp_data(N = 30, T = 2, J = 3, seed = 42)
d <- prepare_hmnp_data(sim$data, "task", "alt", "choice",
                       c("x1", "x2"), person_col = "pid")
out <- hmnp_gibbs(d$X, d$Z, d$M, d$choice_pos, TRUE, d$alt_of_row, d$Ti,
  delta_init = rep(0, d$J), theta_init = rep(0, d$P),
  b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
  nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
  theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
  sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
  a0 = 3, s0 = 3, R = 300, burn = 100, thin = 1, seed = 7,
  keep_beta_i = 1)
colMeans(out$bdraw / sqrt(as.numeric(out$sigma2draw)))
#> [1]  1.2266535 -0.7392549
# }
```
