// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "hb_internal.h"

// ============================================================================
// Hierarchical Bayesian multinomial logit (HMNL) via adaptive RW-Metropolis-
// within-Gibbs (see vignettes/articles/hierarchical_mnl_math.Rmd for the
// model definition and sampler derivation).
//
// Model (inside alternatives j, implicit outside option o):
//   U_ijt = x_ijt' gamma_i + delta_j + EV1,     U_iot = EV1
//   gamma_ik = beta_ik (rc_dist = 0) or exp(beta_ik) (rc_dist = 1)
//   beta_i ~ N(b, W)                       (respondent level, chain scale)
//   delta_j = z_j' theta + xi_j,  xi_j ~ N(0, sigma_d^2)   (alternative level)
// Priors: b ~ N(b_bar, A^{-1}), W ~ IW(nu, V), theta ~ N(theta_bar,
// A_theta^{-1}), sigma_d ~ half-Cauchy(0, s_d) via the Makalic-Schmidt
// scale mixture (or IG(c0, d0) fallback) -- see hb_internal.h.
//
// Gibbs sweep per iteration (systematic scan):
//   (0) fixed-order cache rebuild of the linear predictors v, level-form
//       task denominators D_t = 1 + sum_inside exp(v) (the +1 is the
//       implicit outside good), per-task logliks, per-respondent logliks --
//       needed every iteration because the beta phase changes them, and it
//       bounds incremental floating-point drift to a single sweep;
//   (a) per-respondent RW-MH on beta_i with proposal covariance
//       s_i^2 (H_i + W^{-1})^{-1}, H_i the respondent information at the
//       pooled MNL MLE (precomputed once; Jacobian-adjusted to the chain
//       scale for log-normal coordinates); ONE fresh likelihood pass per
//       candidate, cached loglik for the current state;
//   (b) SERIAL sweep of RW-MH updates on delta_j. The delta_j conditionals
//       are coupled through the shared per-task softmax denominators, so a
//       work-shared update across j would NOT leave the posterior invariant:
//       each delta_j update must see the latest delta and denominator
//       caches (contrast with the HMNP conjugate delta draw -- the
//       asymmetry is documented in hb_internal.h). Cost per alternative is
//       O(#tasks containing j) via the cached D_t and exp(v) terms.
//   (c) master draws: b | beta, W; W | beta, b; theta | delta, sigma_d^2;
//       sigma_d^2 | delta, theta (hb_internal.h conditionals).
//
// Parallelism / reproducibility contract (mirrors src/mnprobit.cpp:20-35):
// ONE persistent OpenMP parallel region for the whole chain. Work-shared
// schedule(static) loops over respondents write only respondent-owned slots
// (each task belongs to exactly one respondent). Worker code paths call only
// the hand-rolled hb_internal.h helpers -- NEVER arma::chol / arma::solve /
// BLAS / LAPACK / the R API off the master thread. The serial delta sweep,
// the hierarchy draws (which MAY use armadillo -- master thread), recording,
// trace, and the interrupt poll run in a master block between barriers. All
// cross-respondent reductions (the loglik trace, the sufficient statistics
// inside draw_b/W_conditional) are accumulated on the master thread in fixed
// respondent order, so draws are reproducible given the seed and a fixed
// thread count; across thread counts they are invariant only up to
// floating-point reduction-order round-off (~1e-15), not bitwise.
//
// RNG partition per iteration r via make_stream(seed, r, tag) (documented in
// hb_internal.h; N respondents, J inside alternatives):
//   tag i           (i = 0..N-1)  respondent i's RW-MH work
//   tag N + i                     UNUSED here (reserved: HMNP beta_i draw)
//   tag 2N + j      (j = 0..J-1)  delta_j RW-MH
//   tag 2N + J + 0..3             b, W, theta, sigma_d^2 (master)
// The master block is based at 2N + J so the partition is collision-free
// for ANY J, including J > N.
//
// Adaptation: log s_i (target `accept_target`, default 0.234) and log s_dj
// (scalar target 0.44) are Robbins-Monro-adapted during burn-in only, with
// vanishing gain (r + 1)^{-0.6}, clamped to [log 0.01, log 10], NaN
// acceptance treated as 0, and FROZEN at r == burn.
//
// Numerical guards: a non-finite candidate loglik (exp overflow of a
// log-normal coordinate or an extreme delta) auto-rejects; D_t >= 1 always
// because of the implicit outside term, so log(D_t) is safe; the startup
// loglik is checked finite (abort otherwise).
// ============================================================================

static void validate_hmnl_inputs(const arma::mat& X, const arma::mat& Z,
                                 const Rcpp::IntegerVector& M,
                                 const Rcpp::IntegerVector& choice_pos,
                                 const bool include_outside_option,
                                 const Rcpp::IntegerVector& alt_of_row,
                                 const Rcpp::IntegerVector& Ti,
                                 const Rcpp::IntegerVector& rc_dist,
                                 const arma::vec& beta_pooled,
                                 const arma::vec& delta_init,
                                 const arma::vec& theta_init,
                                 const arma::vec& b_bar, const arma::mat& A,
                                 const double nu, const arma::mat& V,
                                 const arma::vec& theta_bar,
                                 const arma::mat& A_theta,
                                 const int R, const int burn, const int thin,
                                 const double seed, const int keep_beta_i,
                                 const double s_init,
                                 const double accept_target) {
  const int n_tasks = M.size();
  const int N = Ti.size();
  const int K = X.n_cols;
  const int J = Z.n_rows;
  const int P = Z.n_cols;

  if (!include_outside_option) {
    Rcpp::stop("hmnl_gibbs: include_outside_option = FALSE is not supported: "
               "the implicit outside option is the location anchor that "
               "identifies the level of delta.");
  }
  if (n_tasks < 1) {
    Rcpp::stop("M must contain at least one choice situation.");
  }
  if (N < 1) {
    Rcpp::stop("Ti must contain at least one respondent.");
  }
  if (K < 1) {
    Rcpp::stop("X must have at least one column.");
  }
  if (J < 1 || P < 1) {
    Rcpp::stop("Z must have at least one row (inside alternative) and one "
               "column.");
  }

  long long total_rows = 0;
  long long ti_sum = 0;
  for (int t = 0; t < n_tasks; ++t) {
    if (M[t] < 1) {
      Rcpp::stop("M must be positive for every choice situation (M[%d] = %d).",
                 t + 1, M[t]);
    }
    total_rows += M[t];
  }
  for (int i = 0; i < N; ++i) {
    if (Ti[i] < 1) {
      Rcpp::stop("Ti must be positive for every respondent (Ti[%d] = %d).",
                 i + 1, Ti[i]);
    }
    ti_sum += Ti[i];
  }
  if (total_rows != static_cast<long long>(X.n_rows)) {
    Rcpp::stop("X has %d rows but sum(M) is %d.",
               static_cast<int>(X.n_rows), static_cast<int>(total_rows));
  }
  if (ti_sum != n_tasks) {
    Rcpp::stop("sum(Ti) (%d) does not match the number of choice situations "
               "(%d).", static_cast<int>(ti_sum), n_tasks);
  }
  if (alt_of_row.size() != total_rows) {
    Rcpp::stop("alt_of_row length (%d) does not match the number of rows of "
               "X (%d).", static_cast<int>(alt_of_row.size()),
               static_cast<int>(total_rows));
  }
  for (int rr = 0; rr < alt_of_row.size(); ++rr) {
    if (alt_of_row[rr] < 1 || alt_of_row[rr] > J) {
      Rcpp::stop("alt_of_row[%d] = %d is outside {1, ..., %d}.",
                 rr + 1, alt_of_row[rr], J);
    }
  }
  if (choice_pos.size() != n_tasks) {
    Rcpp::stop("choice_pos length (%d) does not match the number of choice "
               "situations (%d).", static_cast<int>(choice_pos.size()),
               n_tasks);
  }
  for (int t = 0; t < n_tasks; ++t) {
    if (choice_pos[t] < 0 || choice_pos[t] > M[t]) {
      Rcpp::stop("choice_pos[%d] = %d is outside {0, ..., M[%d] = %d} "
                 "(0 = outside option chosen).",
                 t + 1, choice_pos[t], t + 1, M[t]);
    }
  }
  if (!X.is_finite()) {
    Rcpp::stop("X contains non-finite values.");
  }
  if (!Z.is_finite()) {
    Rcpp::stop("Z contains non-finite values.");
  }
  if (rc_dist.size() != K) {
    Rcpp::stop("rc_dist must have length %d (one entry per column of X).", K);
  }
  for (int k = 0; k < K; ++k) {
    if (rc_dist[k] != 0 && rc_dist[k] != 1) {
      Rcpp::stop("rc_dist entries must be 0 (normal) or 1 (log-normal).");
    }
  }
  if (static_cast<int>(beta_pooled.n_elem) != K || !beta_pooled.is_finite()) {
    Rcpp::stop("beta_pooled must be a finite vector of length %d.", K);
  }
  if (static_cast<int>(delta_init.n_elem) != J || !delta_init.is_finite()) {
    Rcpp::stop("delta_init must be a finite vector of length %d.", J);
  }
  if (static_cast<int>(theta_init.n_elem) != P || !theta_init.is_finite()) {
    Rcpp::stop("theta_init must be a finite vector of length %d.", P);
  }
  if (static_cast<int>(b_bar.n_elem) != K) {
    Rcpp::stop("b_bar length (%d) does not match the number of columns of X "
               "(%d).", static_cast<int>(b_bar.n_elem), K);
  }
  if (static_cast<int>(A.n_rows) != K || static_cast<int>(A.n_cols) != K) {
    Rcpp::stop("A must be a %d x %d matrix.", K, K);
  }
  if (static_cast<int>(V.n_rows) != K || static_cast<int>(V.n_cols) != K) {
    Rcpp::stop("V must be a %d x %d matrix.", K, K);
  }
  if (nu < K) {
    Rcpp::stop("nu (%g) must be >= K (%d) for a proper inverse-Wishart prior.",
               nu, K);
  }
  if (static_cast<int>(theta_bar.n_elem) != P) {
    Rcpp::stop("theta_bar length (%d) does not match the number of columns "
               "of Z (%d).", static_cast<int>(theta_bar.n_elem), P);
  }
  if (static_cast<int>(A_theta.n_rows) != P ||
      static_cast<int>(A_theta.n_cols) != P) {
    Rcpp::stop("A_theta must be a %d x %d matrix.", P, P);
  }
  if (R < 1) {
    Rcpp::stop("R must be >= 1.");
  }
  if (burn < 0 || burn >= R) {
    Rcpp::stop("burn must satisfy 0 <= burn < R.");
  }
  if (thin < 1) {
    Rcpp::stop("thin must be >= 1.");
  }
  if (!std::isfinite(seed) || seed < 0) {
    Rcpp::stop("seed must be a finite non-negative number.");
  }
  if (keep_beta_i < 0 || keep_beta_i > 2) {
    Rcpp::stop("keep_beta_i must be 0 (none), 1 (means), or 2 (draws).");
  }
  if (!std::isfinite(s_init) || s_init <= 0) {
    Rcpp::stop("s_init must be a positive number.");
  }
  if (!std::isfinite(accept_target) || accept_target <= 0 ||
      accept_target >= 1) {
    Rcpp::stop("accept_target must lie in (0, 1).");
  }
}

//' Gibbs sampler for the hierarchical Bayesian multinomial logit model
//'
//' Runs the adaptive random-walk Metropolis-within-Gibbs sampler for the
//' hierarchical (random-coefficients, panel) multinomial logit with a
//' BLP-style alternative-level random effect: inside utilities
//' \eqn{U_{ijt} = x_{ijt}'\gamma_i + \delta_j + EV1} against an implicit
//' outside option with systematic utility 0, \eqn{\beta_i \sim N(b, W)}
//' (\eqn{\gamma_{ik} = \beta_{ik}} or \eqn{\exp(\beta_{ik})} per
//' \code{rc_dist}), and \eqn{\delta_j = z_j'\theta + \xi_j},
//' \eqn{\xi_j \sim N(0, \sigma_d^2)}.
//'
//' The per-respondent \eqn{\beta_i} updates are parallelized with OpenMP;
//' the \eqn{\delta_j} updates run as a strictly serial sweep (their
//' conditionals are coupled through the shared softmax denominators). Each
//' (iteration, unit) pair uses its own RNG stream, so draws are
//' reproducible given the seed and a fixed thread count; across different
//' thread counts they are invariant only up to floating-point
//' reduction-order round-off (~1e-15), not bitwise (see
//' \code{set_num_threads()}). This is the low-level engine behind
//' \code{\link{run_hmnlogit}}, which handles initialization and
//' post-processing.
//'
//' @param X total_rows x K_struct structural design matrix (inside rows
//'   only, no ASC columns), rows sorted by (person, task, alternative).
//' @param Z J x P alternative-level mean-function design (intercept first).
//' @param M Integer vector: inside alternatives per choice situation.
//' @param choice_pos Integer vector: 1-based within-task position of the
//'   chosen row; 0 = outside option chosen.
//' @param include_outside_option Must be \code{TRUE} (the implicit outside
//'   good anchors the location of delta; a no-outside mode is roadmapped).
//' @param alt_of_row Integer vector: 1-based alternative code per row of X.
//' @param Ti Integer vector: choice situations per respondent.
//' @param rc_dist Integer vector (length K_struct): 0 = normal coordinate,
//'   1 = log-normal (enters utility as \code{exp(beta_ik)}).
//' @param beta_pooled Pooled MNL MLE on the chain scale (log scale for
//'   log-normal coordinates); centers the H_i proposal information.
//' @param delta_init Initial delta (length J).
//' @param theta_init Initial theta (length P).
//' @param b_bar K vector, prior mean of b.
//' @param A K x K prior precision matrix of b.
//' @param nu Inverse-Wishart prior degrees of freedom for W (>= K).
//' @param V K x K inverse-Wishart prior scale matrix for W.
//' @param theta_bar P vector, prior mean of theta.
//' @param A_theta P x P prior precision matrix of theta.
//' @param sd_prior List with elements \code{half_cauchy} (logical),
//'   \code{s_d} (half-Cauchy scale), \code{c0}, \code{d0} (IG fallback).
//' @param R Total number of Gibbs iterations.
//' @param burn Number of initial iterations discarded (0 <= burn < R);
//'   proposal-scale adaptation happens during burn-in only.
//' @param thin Keep every thin-th post-burn-in draw.
//' @param seed Master RNG seed (non-negative; all streams derive from it).
//' @param keep_beta_i 0 = no beta_i output, 1 = online means/SDs,
//'   2 = means/SDs plus the full (K, N, R_keep) draw cube.
//' @param s_init Initial per-respondent proposal scale.
//' @param accept_target Robbins-Monro acceptance target for the beta_i
//'   updates (the delta_j target is fixed at 0.44).
//' @param trace Print progress every \code{trace} iterations (0 = silent).
//' @returns List with \code{bdraw} (R_keep x K), \code{wdraw} (R_keep x
//'   K(K+1)/2, lower triangle of W in row-major order), \code{deltadraw}
//'   (R_keep x J), \code{thetadraw} (R_keep x P), \code{sigma_d2draw},
//'   \code{loglik_trace}, acceptance rates and final proposal scales
//'   (\code{accept_rate_beta}, \code{accept_rate_delta}, \code{s_final},
//'   \code{s_delta_final}), posterior summaries \code{beta_i_mean} /
//'   \code{beta_i_sd} (K x N, \code{NULL} when \code{keep_beta_i = 0}),
//'   \code{beta_i_draws} (K x N x R_keep cube when \code{keep_beta_i = 2}),
//'   \code{delta_mean} / \code{delta_sd} / \code{xi_mean} / \code{xi_sd}
//'   (J x 1), and \code{R_keep}.
//' @examples
//' \donttest{
//' sim <- simulate_hmnl_data(N = 20, T = 2, J = 3, seed = 42)
//' d <- prepare_hmnl_data(sim$data, "task", "alt", "choice",
//'                        c("x1", "x2"), person_col = "pid")
//' out <- choicer:::hmnl_gibbs(d$X, d$Z, d$M, d$choice_pos, TRUE, d$alt_of_row, d$Ti,
//'   rc_dist = d$rc_dist, beta_pooled = rep(0, d$K_struct),
//'   delta_init = rep(0, d$J), theta_init = rep(0, d$P),
//'   b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
//'   nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
//'   theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
//'   sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
//'   R = 300, burn = 100, thin = 1, seed = 7, keep_beta_i = 1,
//'   s_init = 2.38 / sqrt(d$K_struct), accept_target = 0.234)
//' colMeans(out$bdraw)
//' }
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::List hmnl_gibbs(const arma::mat& X,
                      const arma::mat& Z,
                      const Rcpp::IntegerVector& M,
                      const Rcpp::IntegerVector& choice_pos,
                      const bool include_outside_option,
                      const Rcpp::IntegerVector& alt_of_row,
                      const Rcpp::IntegerVector& Ti,
                      const Rcpp::IntegerVector& rc_dist,
                      const arma::vec& beta_pooled,
                      const arma::vec& delta_init,
                      const arma::vec& theta_init,
                      const arma::vec& b_bar,
                      const arma::mat& A,
                      const double nu,
                      const arma::mat& V,
                      const arma::vec& theta_bar,
                      const arma::mat& A_theta,
                      const Rcpp::List& sd_prior,
                      const int R,
                      const int burn,
                      const int thin,
                      const double seed,
                      const int keep_beta_i,
                      const double s_init,
                      const double accept_target,
                      const int trace = 0) {
  validate_hmnl_inputs(X, Z, M, choice_pos, include_outside_option,
                       alt_of_row, Ti, rc_dist, beta_pooled, delta_init,
                       theta_init, b_bar, A, nu, V, theta_bar, A_theta,
                       R, burn, thin, seed, keep_beta_i, s_init,
                       accept_target);

  const int n_tasks = M.size();
  const int N = Ti.size();
  const int K = X.n_cols;
  const int J = Z.n_rows;
  const int P = Z.n_cols;
  const uint64_t useed = static_cast<uint64_t>(seed);

  // sigma_d prior settings (validated here, not in the worker paths).
  HbSigmaDPrior sdp;
  if (sd_prior.containsElementNamed("half_cauchy")) {
    sdp.half_cauchy = Rcpp::as<bool>(sd_prior["half_cauchy"]);
  }
  if (sd_prior.containsElementNamed("s_d")) {
    sdp.s_d = Rcpp::as<double>(sd_prior["s_d"]);
  }
  if (sd_prior.containsElementNamed("c0")) {
    sdp.c0 = Rcpp::as<double>(sd_prior["c0"]);
  }
  if (sd_prior.containsElementNamed("d0")) {
    sdp.d0 = Rcpp::as<double>(sd_prior["d0"]);
  }
  if (!std::isfinite(sdp.s_d) || sdp.s_d <= 0 ||
      !std::isfinite(sdp.c0) || sdp.c0 <= 0 ||
      !std::isfinite(sdp.d0) || sdp.d0 <= 0) {
    Rcpp::stop("sd_prior: s_d, c0, and d0 must be positive numbers.");
  }

  // Panel indexing (master thread: build() allocates R vectors, then the
  // struct is plain-C++ and read-only inside the region).
  HbPanel panel;
  if (!panel.build(M, Ti, alt_of_row, choice_pos, J,
                   include_outside_option)) {
    Rcpp::stop("hmnl_gibbs: inconsistent panel indexing "
               "(M / Ti / alt_of_row / choice_pos).");
  }
  const int total_rows = panel.total_rows;

  // task -> owning respondent (for the delta sweep's cached-loglik updates).
  std::vector<int> person_of_task(n_tasks);
  for (int i = 0; i < N; ++i) {
    for (int t = panel.task_offsets[i]; t < panel.task_offsets[i + 1]; ++t) {
      person_of_task[t] = i;
    }
  }

  // Plain-int copy of rc_dist so workers never touch an Rcpp vector.
  std::vector<int> rc(K);
  for (int k = 0; k < K; ++k) rc[k] = rc_dist[k];

  // --- Chain state (allocated once, updated in place) -----------------------
  arma::mat beta_i(K, N);
  beta_i.each_col() = beta_pooled;
  arma::vec b = beta_pooled;
  arma::mat W = arma::eye(K, K);
  arma::mat Winv = arma::eye(K, K);
  arma::vec delta = delta_init;
  arma::vec theta = theta_init;
  double sigma_d2 = 1.0;
  double a_d = 1.0;

  // Adaptive proposal scales (log scale, clamped, frozen at r == burn).
  const double LOG_S_MIN = std::log(0.01);
  const double LOG_S_MAX = std::log(10.0);
  const double DELTA_TARGET = 0.44;   // scalar-update optimum
  arma::vec s_log(N);
  s_log.fill(std::log(s_init));
  arma::vec s_delta_log(J, arma::fill::zeros);   // s_delta = 1 at start

  // Per-iteration caches (level form; disjoint respondent-owned slots).
  arma::vec v(total_rows), e(total_rows);
  arma::vec Dv(n_tasks), ll_task(n_tasks), ll_pers(N);
  // Candidate buffers for the beta phase (same disjoint-slot layout).
  arma::vec v_cand(total_rows), e_cand(total_rows);
  arma::vec D_cand(n_tasks), ll_task_cand(n_tasks);

  // Respondent information at the pooled init (filled in the region).
  arma::cube H(K, K, N);

  // Post-burn acceptance counters (respondent- / alternative-owned slots).
  arma::vec acc_beta(N, arma::fill::zeros);
  arma::vec acc_delta(J, arma::fill::zeros);

  // Master scratch for the delta sweep / hierarchy draws.
  arma::vec mu_delta(J), xi_tmp(J);

  // --- Recording -------------------------------------------------------------
  const int R_keep = (R - burn + thin - 1) / thin;
  arma::mat bdraw(R_keep, K);
  arma::mat wdraw(R_keep, K * (K + 1) / 2);
  arma::mat deltadraw(R_keep, J);
  arma::mat thetadraw(R_keep, P);
  arma::vec sigma_d2draw(R_keep);
  arma::vec loglik_trace(R_keep);
  HbWelford beta_wf, delta_wf, xi_wf;
  if (keep_beta_i >= 1) beta_wf.reset(K, N);
  delta_wf.reset(J, 1);
  xi_wf.reset(J, 1);
  arma::cube beta_draws;
  if (keep_beta_i == 2) beta_draws.set_size(K, N, R_keep);

  // Abort protocol (mnprobit pattern): failures inside the region set a code
  // (workers via a critical section), every thread re-reads it after a
  // barrier and leaves together; Rcpp::stop happens after the region.
  enum { ABORT_NONE = 0, ABORT_CHOL, ABORT_INTERRUPT, ABORT_B, ABORT_W,
         ABORT_WINV, ABORT_THETA, ABORT_INIT };
  int abort_code = ABORT_NONE;
  int abort_iter = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local workspaces (never shared across respondents).
    arma::vec gamma_loc(K), z_loc(K), t_loc(K), d_loc(K), xrow(K), xbar(K);
    arma::vec beta_star(K);
    arma::mat P_loc(K, K), L_loc(K, K);

    // --- Startup: respondent information H_i at the pooled init -------------
    // I_t = sum_j p_jt x_jt x_jt' - xbar_t xbar_t', xbar_t = sum_j p_jt x_jt,
    // with the implicit outside contributing x = 0 and its softmax mass;
    // H_i = sum_{t in i} I_t, then Jacobian-adjusted to the chain scale for
    // log-normal coordinates: multiply row & column k by exp(beta_pooled_k).
    // Deterministic (no RNG), work-shared, hand-rolled fixed-order loops.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < N; ++i) {
      for (int k = 0; k < K; ++k) {
        gamma_loc(k) = (rc[k] == 1) ? std::exp(beta_pooled(k))
                                    : beta_pooled(k);
      }
      arma::mat& Hi = H.slice(i);
      Hi.zeros();
      for (int t = panel.task_offsets[i]; t < panel.task_offsets[i + 1];
           ++t) {
        const int r0 = panel.row_offsets[t];
        const int r1 = panel.row_offsets[t + 1];
        double vmax = 0.0;                       // outside term v = 0
        for (int rr = r0; rr < r1; ++rr) {
          double s = delta_init(panel.alt_of_row[rr]);
          for (int k = 0; k < K; ++k) s += X(rr, k) * gamma_loc(k);
          v(rr) = s;                             // scratch use; rebuilt at r=0
          if (s > vmax) vmax = s;
        }
        double denom = std::exp(0.0 - vmax);
        for (int rr = r0; rr < r1; ++rr) denom += std::exp(v(rr) - vmax);
        xbar.zeros();
        for (int rr = r0; rr < r1; ++rr) {
          const double p = std::exp(v(rr) - vmax) / denom;
          for (int k = 0; k < K; ++k) xrow(k) = X(rr, k);
          hb_rank1_accum(Hi, xrow, p);
          for (int k = 0; k < K; ++k) xbar(k) += p * xrow(k);
        }
        hb_rank1_accum(Hi, xbar, -1.0);
      }
      for (int k = 0; k < K; ++k) {
        if (rc[k] == 1) {
          const double g = std::exp(beta_pooled(k));
          for (int l = 0; l < K; ++l) {
            Hi(k, l) *= g;
            Hi(l, k) *= g;
          }
        }
      }
    }   // implied barrier: H complete before the chain starts

    for (int r = 0; r < R; ++r) {
      // --- (0) Cache rebuild: v, e, D_t, ll_t, cached respondent logliks ----
      // Fixed sub-order per respondent; slots are respondent-owned, so this
      // pass does not depend on the thread schedule.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
          gamma_loc(k) = (rc[k] == 1) ? std::exp(beta_i(k, i))
                                      : beta_i(k, i);
        }
        double lli = 0.0;
        for (int t = panel.task_offsets[i]; t < panel.task_offsets[i + 1];
             ++t) {
          const int r0 = panel.row_offsets[t];
          const int r1 = panel.row_offsets[t + 1];
          double Dt = 1.0;                       // implicit outside exp(0)
          for (int rr = r0; rr < r1; ++rr) {
            double s = delta(panel.alt_of_row[rr]);
            for (int k = 0; k < K; ++k) s += X(rr, k) * gamma_loc(k);
            v(rr) = s;
            e(rr) = std::exp(s);
            Dt += e(rr);
          }
          Dv(t) = Dt;
          const int cr = panel.chosen_row[t];
          const double llt = ((cr < 0) ? 0.0 : v(cr)) - std::log(Dt);
          ll_task(t) = llt;
          lli += llt;
        }
        ll_pers(i) = lli;
      }   // implied barrier

      // Startup guard: the initial state must have a finite loglik (a bad
      // beta_pooled / delta_init would otherwise poison the MH ratios).
      if (r == 0) {
        CHOICER_OMP_MASKED
        {
          for (int i = 0; i < N; ++i) {
            if (!std::isfinite(ll_pers(i))) {
              abort_code = ABORT_INIT;
              abort_iter = 0;
              break;
            }
          }
        }
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (abort_code != ABORT_NONE) break;
      }

      // --- (a) beta phase: per-respondent RW-MH (work-shared) ---------------
      // Proposal beta* = beta_i + s_i L^{-T} z with L = chol(H_i + W^{-1})
      // (hand-rolled -- worker thread). One fixed-order candidate loglik
      // pass over the respondent's own rows; non-finite => auto-reject.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        Xoshiro256pp rng = make_stream(useed, static_cast<uint64_t>(r),
                                       static_cast<uint64_t>(i));
        const arma::mat& Hi = H.slice(i);
        for (int cj = 0; cj < K; ++cj) {
          for (int ci = 0; ci < K; ++ci) {
            P_loc(ci, cj) = Hi(ci, cj) + Winv(ci, cj);
          }
        }
        if (!hb_chol_lower(P_loc, L_loc)) {
#ifdef _OPENMP
#pragma omp critical
#endif
          {
            if (abort_code == ABORT_NONE) {
              abort_code = ABORT_CHOL;
              abort_iter = r;
            }
          }
          continue;
        }
        const double s_i = std::exp(s_log(i));
        for (int k = 0; k < K; ++k) z_loc(k) = rng.rnorm();
        hb_back_solve(L_loc, z_loc, t_loc);
        for (int k = 0; k < K; ++k) {
          beta_star(k) = beta_i(k, i) + s_i * t_loc(k);
        }
        for (int k = 0; k < K; ++k) {
          gamma_loc(k) = (rc[k] == 1) ? std::exp(beta_star(k))
                                      : beta_star(k);
        }

        double cand_ll = 0.0;
        for (int t = panel.task_offsets[i]; t < panel.task_offsets[i + 1];
             ++t) {
          const int r0 = panel.row_offsets[t];
          const int r1 = panel.row_offsets[t + 1];
          double Dt = 1.0;
          for (int rr = r0; rr < r1; ++rr) {
            double s = delta(panel.alt_of_row[rr]);
            for (int k = 0; k < K; ++k) s += X(rr, k) * gamma_loc(k);
            v_cand(rr) = s;
            const double es = std::exp(s);
            e_cand(rr) = es;
            Dt += es;
          }
          D_cand(t) = Dt;
          const int cr = panel.chosen_row[t];
          const double llt = ((cr < 0) ? 0.0 : v_cand(cr)) - std::log(Dt);
          ll_task_cand(t) = llt;
          cand_ll += llt;
        }
        if (!std::isfinite(cand_ll)) cand_ll = -arma::datum::inf;

        for (int k = 0; k < K; ++k) d_loc(k) = beta_i(k, i) - b(k);
        const double q_cur = hb_quad_form(Winv, d_loc);
        for (int k = 0; k < K; ++k) d_loc(k) = beta_star(k) - b(k);
        const double q_star = hb_quad_form(Winv, d_loc);
        const double log_alpha =
          (cand_ll - ll_pers(i)) + 0.5 * (q_cur - q_star);

        const double u = rng.runif();
        if (std::log(u) < log_alpha) {           // NaN log_alpha => reject
          for (int k = 0; k < K; ++k) beta_i(k, i) = beta_star(k);
          for (int t = panel.task_offsets[i];
               t < panel.task_offsets[i + 1]; ++t) {
            const int r0 = panel.row_offsets[t];
            const int r1 = panel.row_offsets[t + 1];
            for (int rr = r0; rr < r1; ++rr) {
              v(rr) = v_cand(rr);
              e(rr) = e_cand(rr);
            }
            Dv(t) = D_cand(t);
            ll_task(t) = ll_task_cand(t);
          }
          ll_pers(i) = cand_ll;
          if (r >= burn) acc_beta(i) += 1.0;
        }

        if (r < burn) {                          // adapt during burn-in ONLY
          double ap = std::exp(std::min(0.0, log_alpha));
          if (std::isnan(ap)) ap = 0.0;
          const double gain = std::pow(static_cast<double>(r + 1), -0.6);
          double ls = s_log(i) + gain * (ap - accept_target);
          if (ls < LOG_S_MIN) ls = LOG_S_MIN;
          if (ls > LOG_S_MAX) ls = LOG_S_MAX;
          s_log(i) = ls;
        }
      }   // implied barrier: beta phase complete
      if (abort_code != ABORT_NONE) break;

      // --- (b) delta phase (SERIAL, master) + (c) hierarchy + recording -----
      CHOICER_OMP_MASKED
      {
        // Prior means z_j' theta, fixed alternative order (theta is fixed
        // for the whole sweep -- it is redrawn only after delta).
        for (int j = 0; j < J; ++j) {
          double m = 0.0;
          for (int p = 0; p < P; ++p) m += Z(j, p) * theta(p);
          mu_delta(j) = m;
        }

        // Serial sweep: delta_{j+1} sees the latest delta and the updated
        // caches. NOT parallelizable (softmax coupling across j).
        for (int j = 0; j < J; ++j) {
          Xoshiro256pp rng_d = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + j));
          const double sdj = std::exp(s_delta_log(j));
          const double dstar = delta(j) + sdj * rng_d.rnorm();
          const double dd = dstar - delta(j);

          // Delta-loglik over the tasks containing alternative j only:
          // each task contributes dd * 1{chosen row is alt j}
          // - log(D_t - e_jt + e*_jt) + log(D_t), all from the caches.
          double dll = 0.0;
          bool ok = true;
          for (int idx = panel.alt_row_offsets[j];
               idx < panel.alt_row_offsets[j + 1]; ++idx) {
            const int rr = panel.alt_rows[idx];
            const int t = panel.task_of_row[rr];
            const double e_new = std::exp(v(rr) + dd);
            if (!std::isfinite(e_new)) {         // overflow => auto-reject
              ok = false;
              break;
            }
            const double D_new = Dv(t) - e(rr) + e_new;
            dll += std::log(Dv(t)) - std::log(D_new);
            if (panel.chosen_row[t] == rr) dll += dd;
          }

          double log_alpha = -arma::datum::inf;
          if (ok) {
            const double rc_cur = delta(j) - mu_delta(j);
            const double rc_star = dstar - mu_delta(j);
            log_alpha = dll +
              0.5 * (rc_cur * rc_cur - rc_star * rc_star) / sigma_d2;
          }

          const double u = rng_d.runif();
          if (std::log(u) < log_alpha) {         // NaN log_alpha => reject
            for (int idx = panel.alt_row_offsets[j];
                 idx < panel.alt_row_offsets[j + 1]; ++idx) {
              const int rr = panel.alt_rows[idx];
              const int t = panel.task_of_row[rr];
              const double e_old = e(rr);
              v(rr) += dd;
              const double e_new = std::exp(v(rr));   // == exp(v_old + dd)
              e(rr) = e_new;
              Dv(t) = Dv(t) - e_old + e_new;
              const int cr = panel.chosen_row[t];
              const double ll_new =
                ((cr < 0) ? 0.0 : v(cr)) - std::log(Dv(t));
              ll_pers(person_of_task[t]) += ll_new - ll_task(t);
              ll_task(t) = ll_new;
            }
            delta(j) = dstar;
            if (r >= burn) acc_delta(j) += 1.0;
          }

          if (r < burn) {                        // adapt during burn-in ONLY
            double ap = ok ? std::exp(std::min(0.0, log_alpha)) : 0.0;
            if (std::isnan(ap)) ap = 0.0;
            const double gain = std::pow(static_cast<double>(r + 1), -0.6);
            double ls = s_delta_log(j) + gain * (ap - DELTA_TARGET);
            if (ls < LOG_S_MIN) ls = LOG_S_MIN;
            if (ls > LOG_S_MAX) ls = LOG_S_MAX;
            s_delta_log(j) = ls;
          }
        }

        // --- (c) hierarchy draws (master; conditionals from hb_internal.h,
        // fixed-order sufficient statistics; arma::inv_sympd is fine HERE).
        {
          Xoshiro256pp rng_b = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 0));
          if (!draw_b_conditional(rng_b, beta_i, b_bar, A, Winv, b)) {
            abort_code = ABORT_B;
            abort_iter = r;
          }
        }
        if (abort_code == ABORT_NONE) {
          Xoshiro256pp rng_w = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 1));
          if (!draw_W_conditional(rng_w, beta_i, b, nu, V, W)) {
            abort_code = ABORT_W;
            abort_iter = r;
          } else if (!arma::inv_sympd(Winv, W)) {
            abort_code = ABORT_WINV;
            abort_iter = r;
          }
        }
        if (abort_code == ABORT_NONE) {
          Xoshiro256pp rng_t = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 2));
          if (!draw_theta_conditional(rng_t, Z, delta, sigma_d2, theta_bar,
                                      A_theta, theta)) {
            abort_code = ABORT_THETA;
            abort_iter = r;
          }
        }
        if (abort_code == ABORT_NONE) {
          Xoshiro256pp rng_s = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 3));
          for (int j = 0; j < J; ++j) {          // xi = delta - Z theta
            double m = 0.0;
            for (int p = 0; p < P; ++p) m += Z(j, p) * theta(p);
            xi_tmp(j) = delta(j) - m;
          }
          draw_sigma_d2_conditional(rng_s, xi_tmp, sdp, sigma_d2, a_d);
        }

        // --- Recording and housekeeping (master, inside the region) --------
        if (abort_code == ABORT_NONE && r >= burn &&
            (r - burn) % thin == 0) {
          const int row = (r - burn) / thin;
          for (int k = 0; k < K; ++k) bdraw(row, k) = b(k);
          int idx = 0;
          for (int a1 = 0; a1 < K; ++a1) {
            for (int a2 = 0; a2 <= a1; ++a2) {
              wdraw(row, idx++) = W(a1, a2);
            }
          }
          for (int j = 0; j < J; ++j) deltadraw(row, j) = delta(j);
          for (int p = 0; p < P; ++p) thetadraw(row, p) = theta(p);
          sigma_d2draw(row) = sigma_d2;
          double llsum = 0.0;                    // fixed respondent order
          for (int i = 0; i < N; ++i) llsum += ll_pers(i);
          loglik_trace(row) = llsum;
          if (keep_beta_i >= 1) beta_wf.update(beta_i);
          if (keep_beta_i == 2) beta_draws.slice(row) = beta_i;
          delta_wf.update(delta);
          xi_wf.update(xi_tmp);
        }
        if (abort_code == ABORT_NONE) {
          if (trace > 0 && (r + 1) % trace == 0) {
            Rprintf("hmnl_gibbs: iteration %d / %d\n", r + 1, R);
          }
          if ((r + 1) % 100 == 0 && hb_pending_interrupt()) {
            abort_code = ABORT_INTERRUPT;
            abort_iter = r;
          }
        }
      }   // end master block
#ifdef _OPENMP
#pragma omp barrier
#endif
      if (abort_code != ABORT_NONE) break;
    }   // end iteration loop
  }   // end parallel region

  switch (abort_code) {
    case ABORT_CHOL:
      Rcpp::stop("hmnl_gibbs: proposal precision H_i + W^{-1} is not "
                 "positive definite at iteration %d.", abort_iter + 1);
    case ABORT_B:
      Rcpp::stop("hmnl_gibbs: posterior precision of b is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_W:
      Rcpp::stop("hmnl_gibbs: W scale matrix is not positive definite at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_WINV:
      Rcpp::stop("hmnl_gibbs: W draw is numerically singular at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_THETA:
      Rcpp::stop("hmnl_gibbs: posterior precision of theta is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_INIT:
      Rcpp::stop("hmnl_gibbs: initial log-likelihood is not finite; check "
                 "beta_pooled and delta_init.");
    case ABORT_INTERRUPT:
      Rcpp::stop("hmnl_gibbs: interrupted by user at iteration %d.",
                 abort_iter + 1);
    default:
      break;
  }

  const double denom_post = static_cast<double>(R - burn);
  arma::vec accept_rate_beta = acc_beta / denom_post;
  arma::vec accept_rate_delta = acc_delta / denom_post;
  arma::vec s_final = arma::exp(s_log);
  arma::vec s_delta_final = arma::exp(s_delta_log);

  SEXP beta_i_mean = R_NilValue, beta_i_sd = R_NilValue,
       beta_i_draws = R_NilValue;
  if (keep_beta_i >= 1) {
    beta_i_mean = Rcpp::wrap(beta_wf.mean);
    beta_i_sd = Rcpp::wrap(beta_wf.sd());
  }
  if (keep_beta_i == 2) {
    beta_i_draws = Rcpp::wrap(beta_draws);
  }

  return Rcpp::List::create(
    Rcpp::Named("bdraw") = bdraw,
    Rcpp::Named("wdraw") = wdraw,
    Rcpp::Named("deltadraw") = deltadraw,
    Rcpp::Named("thetadraw") = thetadraw,
    Rcpp::Named("sigma_d2draw") = sigma_d2draw,
    Rcpp::Named("loglik_trace") = loglik_trace,
    Rcpp::Named("accept_rate_beta") = accept_rate_beta,
    Rcpp::Named("accept_rate_delta") = accept_rate_delta,
    Rcpp::Named("s_final") = s_final,
    Rcpp::Named("s_delta_final") = s_delta_final,
    Rcpp::Named("beta_i_mean") = beta_i_mean,
    Rcpp::Named("beta_i_sd") = beta_i_sd,
    Rcpp::Named("beta_i_draws") = beta_i_draws,
    Rcpp::Named("delta_mean") = delta_wf.mean,
    Rcpp::Named("delta_sd") = delta_wf.sd(),
    Rcpp::Named("xi_mean") = xi_wf.mean,
    Rcpp::Named("xi_sd") = xi_wf.sd(),
    Rcpp::Named("R_keep") = R_keep
  );
}
