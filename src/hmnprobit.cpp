// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "hb_internal.h"

// ============================================================================
// Hierarchical Bayesian multinomial probit with iid utility shocks (HMNP)
// via Albert-Chib data augmentation in UN-differenced utility space (see
// _plans/hierarchical_bayes_plan.md, "Models"/"C++ kernel design").
//
// Model (inside alternatives j, implicit outside option o):
//   U_ijt = x_ijt' beta_i + delta_j + eps,   U_iot = eps,   eps ~ iid N(0, s2)
//   beta_i ~ N(b, W)                          (respondent level, normal only)
//   delta_j = z_j' theta + xi_j,  xi_j ~ N(0, sigma_d^2)  (alternative level)
// Choice = argmax within the task, INCLUDING the outside latent: the outside
// good is stochastic (its own shock on top of systematic utility 0), NOT a
// deterministic bound — this keeps the model parallel to the HMNL (whose
// outside carries an EV1 shock) and gives the equicorrelated differenced-
// error structure that is the probit analog of logit-with-outside.
//
// Because the shocks are iid, the latent conditionals are univariate normals
// whose only coupling is through the truncation bounds (chosen > all others,
// all others < chosen), so the within-task sweep is a plain scalar TN scan.
// Un-differenced space keeps unbalanced choice sets trivial.
//
// The chain runs on the NON-identified parameterization (free s2 = sigma^2);
// the probit likelihood is invariant to a common rescaling of (U, sigma), so
// run_hmnprobit divides every kept draw by the matching power of the current
// sigma (the McCulloch-Rossi sigma_11 pattern, doubling as PX-DA for better
// mixing). Raw chains are returned; the per-draw-normalized beta_i / delta /
// xi summaries are accumulated here at recording time (identified scale).
//
// Gibbs sweep per iteration (systematic scan, all conjugate):
//   (a) fused work-shared pass over respondents: TN sweep over the
//       respondent's rows plus each task's implicit outside latent (fixed
//       sub-order; the outside latent participates exactly like a row with
//       mean 0), then the conjugate beta_i regression draw on the
//       delta-residualized utilities (precision W^{-1} + G_i / s2, G_i =
//       X_i'X_i precomputed), then a refresh of the row cache
//       mu_x(r) = x_r' beta_i with the new beta_i;
//   (b) work-shared conjugate delta_j draws — VALID in parallel here: given
//       the augmented U the delta_j are conditionally independent across j
//       (contrast the HMNL's coupled softmax conditionals; the asymmetry is
//       documented in hb_internal.h). Each j's sufficient statistic is
//       summed over its CSR row list in fixed row order by a single thread;
//   (c) work-shared RSS pass with the NEW beta and delta (a separate pass so
//       sigma^2 conditions on the current values — reusing stats from (a)
//       would condition on the stale delta), reduced on master in fixed
//       respondent order;
//   (d) master draws: b, W, theta, sigma_d^2 (hb_internal.h conditionals)
//       and s2 ~ IG(a0 + n_latents/2, s0 + RSS/2), where n_latents counts
//       the inside rows AND the outside latents (their residual is the
//       latent itself).
//
// Parallelism / reproducibility contract (mirrors src/mnprobit.cpp:20-35):
// ONE persistent OpenMP region; workers touch only hand-rolled hb_internal.h
// helpers (no BLAS/LAPACK/R API off master); every cross-unit reduction is
// fixed-order; per-(iteration, unit) RNG streams make draws bitwise
// independent of the thread count.
//
// RNG partition per iteration r (make_stream(seed, r, tag), documented in
// hb_internal.h):
//   tag i           (i = 0..N-1)  respondent i's TN sweep
//   tag N + i       (i = 0..N-1)  respondent i's beta_i regression draw
//   tag 2N + j      (j = 0..J-1)  delta_j draw
//   tag 2N + J + 0..3             b, W, theta, (sigma_d^2 then s2) (master)
// ============================================================================

static void validate_hmnp_inputs(const arma::mat& X, const arma::mat& Z,
                                 const Rcpp::IntegerVector& M,
                                 const Rcpp::IntegerVector& choice_pos,
                                 const bool include_outside_option,
                                 const Rcpp::IntegerVector& alt_of_row,
                                 const Rcpp::IntegerVector& Ti,
                                 const arma::vec& delta_init,
                                 const arma::vec& theta_init,
                                 const arma::vec& b_bar, const arma::mat& A,
                                 const double nu, const arma::mat& V,
                                 const arma::vec& theta_bar,
                                 const arma::mat& A_theta,
                                 const double a0, const double s0,
                                 const int R, const int burn, const int thin,
                                 const double seed, const int keep_beta_i) {
  const int n_tasks = M.size();
  const int N = Ti.size();
  const int K = X.n_cols;
  const int J = Z.n_rows;
  const int P = Z.n_cols;

  if (!include_outside_option) {
    Rcpp::stop("hmnp_gibbs: include_outside_option = FALSE is not supported: "
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
  if (!std::isfinite(a0) || a0 <= 0 || !std::isfinite(s0) || s0 <= 0) {
    Rcpp::stop("a0 and s0 (the sigma^2 prior) must be positive numbers.");
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
}

//' Gibbs sampler for the hierarchical Bayesian multinomial probit model
//'
//' Runs the fully conjugate Albert-Chib Gibbs sampler for the hierarchical
//' multinomial probit with iid \eqn{N(0, \sigma^2)} utility shocks in
//' un-differenced utility space: inside utilities
//' \eqn{U_{ijt} = x_{ijt}'\beta_i + \delta_j + \epsilon} against a
//' stochastic implicit outside option \eqn{U_{iot} = \epsilon},
//' \eqn{\beta_i \sim N(b, W)}, and \eqn{\delta_j = z_j'\theta + \xi_j},
//' \eqn{\xi_j \sim N(0, \sigma_d^2)}. The chain runs on the non-identified
//' parameterization (free \eqn{\sigma^2}); identified quantities are
//' obtained by normalizing each draw by the matching power of \eqn{\sigma}
//' (handled by \code{\link{run_hmnprobit}}).
//'
//' The latent sweep and the \eqn{\beta_i} draws are parallelized with
//' OpenMP across respondents; the \eqn{\delta_j} draws are parallelized
//' across alternatives (conditionally independent given the augmented
//' utilities — unlike the HMNL, whose delta sweep must be serial). Each
//' (iteration, unit) pair uses its own RNG stream, so draws are bitwise
//' reproducible independent of the number of threads.
//'
//' @param X total_rows x K_struct structural design matrix (inside rows
//'   only), rows sorted by (person, task, alternative).
//' @param Z J x P alternative-level mean-function design (intercept first).
//' @param M Integer vector: inside alternatives per choice situation.
//' @param choice_pos Integer vector: 1-based within-task position of the
//'   chosen row; 0 = outside option chosen.
//' @param include_outside_option Must be \code{TRUE} (the outside good
//'   anchors the location of delta; a no-outside mode is roadmapped).
//' @param alt_of_row Integer vector: 1-based alternative code per row of X.
//' @param Ti Integer vector: choice situations per respondent.
//' @param delta_init Initial delta (length J), raw scale.
//' @param theta_init Initial theta (length P), raw scale.
//' @param b_bar K vector, prior mean of b.
//' @param A K x K prior precision matrix of b.
//' @param nu Inverse-Wishart prior degrees of freedom for W (>= K).
//' @param V K x K inverse-Wishart prior scale matrix for W.
//' @param theta_bar P vector, prior mean of theta.
//' @param A_theta P x P prior precision matrix of theta.
//' @param sd_prior List with elements \code{half_cauchy} (logical),
//'   \code{s_d}, \code{c0}, \code{d0} — the sigma_d prior (see
//'   \code{\link{run_hmnlogit}}).
//' @param a0,s0 Inverse-gamma prior shape/scale for the (non-identified)
//'   shock variance \eqn{\sigma^2}.
//' @param R Total number of Gibbs iterations.
//' @param burn Number of initial iterations discarded (0 <= burn < R).
//' @param thin Keep every thin-th post-burn-in draw.
//' @param seed Master RNG seed (non-negative; all streams derive from it).
//' @param keep_beta_i 0 = no beta_i output, 1 = online means/SDs
//'   (identified scale), 2 = means/SDs plus the full (K, N, R_keep) cube of
//'   per-draw-normalized \eqn{\beta_i / \sigma} draws.
//' @param trace Print progress every \code{trace} iterations (0 = silent).
//' @returns List with RAW draw matrices \code{bdraw}, \code{wdraw} (lower
//'   triangle, row-major), \code{deltadraw}, \code{thetadraw},
//'   \code{sigma_d2draw}, \code{sigma2draw}, identified-scale summaries
//'   \code{beta_i_mean} / \code{beta_i_sd} / \code{beta_i_draws} /
//'   \code{delta_mean} / \code{delta_sd} / \code{xi_mean} / \code{xi_sd},
//'   and \code{R_keep}.
//' @examples
//' \donttest{
//' sim <- simulate_hmnp_data(N = 30, T = 2, J = 3, seed = 42)
//' d <- prepare_hmnp_data(sim$data, "task", "alt", "choice",
//'                        c("x1", "x2"), person_col = "pid")
//' out <- hmnp_gibbs(d$X, d$Z, d$M, d$choice_pos, TRUE, d$alt_of_row, d$Ti,
//'   delta_init = rep(0, d$J), theta_init = rep(0, d$P),
//'   b_bar = rep(0, d$K_struct), A = 0.01 * diag(d$K_struct),
//'   nu = d$K_struct + 3, V = (d$K_struct + 3) * diag(d$K_struct),
//'   theta_bar = rep(0, d$P), A_theta = 0.01 * diag(d$P),
//'   sd_prior = list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
//'   a0 = 3, s0 = 3, R = 300, burn = 100, thin = 1, seed = 7,
//'   keep_beta_i = 1)
//' colMeans(out$bdraw / sqrt(as.numeric(out$sigma2draw)))
//' }
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hmnp_gibbs(const arma::mat& X,
                      const arma::mat& Z,
                      const Rcpp::IntegerVector& M,
                      const Rcpp::IntegerVector& choice_pos,
                      const bool include_outside_option,
                      const Rcpp::IntegerVector& alt_of_row,
                      const Rcpp::IntegerVector& Ti,
                      const arma::vec& delta_init,
                      const arma::vec& theta_init,
                      const arma::vec& b_bar,
                      const arma::mat& A,
                      const double nu,
                      const arma::mat& V,
                      const arma::vec& theta_bar,
                      const arma::mat& A_theta,
                      const Rcpp::List& sd_prior,
                      const double a0,
                      const double s0,
                      const int R,
                      const int burn,
                      const int thin,
                      const double seed,
                      const int keep_beta_i,
                      const int trace = 0) {
  validate_hmnp_inputs(X, Z, M, choice_pos, include_outside_option,
                       alt_of_row, Ti, delta_init, theta_init, b_bar, A, nu,
                       V, theta_bar, A_theta, a0, s0, R, burn, thin, seed,
                       keep_beta_i);

  const int n_tasks = M.size();
  const int N = Ti.size();
  const int K = X.n_cols;
  const int J = Z.n_rows;
  const int P = Z.n_cols;
  const uint64_t useed = static_cast<uint64_t>(seed);

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

  HbPanel panel;
  if (!panel.build(M, Ti, alt_of_row, choice_pos, J,
                   include_outside_option)) {
    Rcpp::stop("hmnp_gibbs: inconsistent panel indexing "
               "(M / Ti / alt_of_row / choice_pos).");
  }
  const int total_rows = panel.total_rows;
  // n_latents counts the outside latents too: their residual is the latent
  // itself and they enter the sigma^2 conditional.
  const double n_latents =
    static_cast<double>(total_rows) + static_cast<double>(n_tasks);

  // --- Chain state (allocated once, updated in place) -----------------------
  arma::mat beta_i(K, N, arma::fill::zeros);
  arma::vec b(K, arma::fill::zeros);
  arma::mat W = arma::eye(K, K);
  arma::mat Winv = arma::eye(K, K);
  arma::vec delta = delta_init;
  arma::vec theta = theta_init;
  double sigma_d2 = 1.0;
  double a_d = 1.0;
  double s2 = 1.0;

  // Latent utilities: inside rows (U) and the per-task outside latent
  // (U_out). Feasible start: chosen latent at +1, everything else at -1
  // (the mnprobit convention).
  arma::vec U(total_rows), U_out(n_tasks);
  for (int t = 0; t < n_tasks; ++t) {
    const int cr = panel.chosen_row[t];
    for (int rr = panel.row_offsets[t]; rr < panel.row_offsets[t + 1]; ++rr) {
      U(rr) = (rr == cr) ? 1.0 : -1.0;
    }
    U_out(t) = (cr < 0) ? 1.0 : -1.0;
  }

  // Row cache mu_x(r) = x_r' beta_{person(r)}; refreshed inside the fused
  // pass right after each beta_i draw, and consumed by the delta and RSS
  // passes (which then only add the current delta).
  arma::vec mu_x(total_rows, arma::fill::zeros);

  // Per-respondent RSS slots, reduced on master in fixed order.
  arma::vec rss_slot(N, arma::fill::zeros);

  // Precomputed Gram blocks G_i = X_i' X_i (structural K only).
  arma::cube G(K, K, N);

  // --- Recording -------------------------------------------------------------
  const int R_keep = (R - burn + thin - 1) / thin;
  arma::mat bdraw(R_keep, K);
  arma::mat wdraw(R_keep, K * (K + 1) / 2);
  arma::mat deltadraw(R_keep, J);
  arma::mat thetadraw(R_keep, P);
  arma::vec sigma_d2draw(R_keep);
  arma::vec sigma2draw(R_keep);
  HbWelford beta_wf, delta_wf, xi_wf;
  if (keep_beta_i >= 1) beta_wf.reset(K, N);
  delta_wf.reset(J, 1);
  xi_wf.reset(J, 1);
  arma::cube beta_draws;
  if (keep_beta_i == 2) beta_draws.set_size(K, N, R_keep);
  arma::mat beta_id_scratch(K, N);
  arma::vec delta_id_scratch(J), xi_id_scratch(J);

  enum { ABORT_NONE = 0, ABORT_CHOL, ABORT_INTERRUPT, ABORT_B, ABORT_W,
         ABORT_WINV, ABORT_THETA };
  int abort_code = ABORT_NONE;
  int abort_iter = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local workspaces.
    arma::vec rhs_loc(K), mean_loc(K), z_loc(K), t_loc(K), bnew(K);
    arma::mat Q_loc(K, K), L_loc(K, K);

    // --- Startup: Gram blocks G_i (work-shared, deterministic) --------------
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < N; ++i) {
      arma::mat& Gi = G.slice(i);
      Gi.zeros();
      const int t0 = panel.task_offsets[i];
      const int t1 = panel.task_offsets[i + 1];
      for (int t = t0; t < t1; ++t) {
        for (int rr = panel.row_offsets[t]; rr < panel.row_offsets[t + 1];
             ++rr) {
          for (int cj = 0; cj < K; ++cj) {         // fixed-order rank-1
            const double xj = X(rr, cj);
            for (int ci = 0; ci < K; ++ci) {
              Gi(ci, cj) += X(rr, ci) * xj;
            }
          }
        }
      }
    }   // implied barrier

    for (int r = 0; r < R; ++r) {
      const double sd_eps = std::sqrt(s2);
      const double inv_s2 = 1.0 / s2;

      // --- (a) Fused pass: TN sweep + beta_i draw + mu_x refresh ------------
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        Xoshiro256pp rng_u = make_stream(useed, static_cast<uint64_t>(r),
                                         static_cast<uint64_t>(i));
        const int t0 = panel.task_offsets[i];
        const int t1 = panel.task_offsets[i + 1];

        // TN sweep, fixed sub-order: inside rows then the outside latent,
        // task by task. The iid structure makes each conditional a plain
        // N(mean, s2) truncated by the argmax constraint given the CURRENT
        // values of the other latents in the task.
        for (int t = t0; t < t1; ++t) {
          const int r0 = panel.row_offsets[t];
          const int r1 = panel.row_offsets[t + 1];
          const int cr = panel.chosen_row[t];   // -1 = outside chosen

          for (int rr = r0; rr < r1; ++rr) {
            const double m = mu_x(rr) + delta(panel.alt_of_row[rr]);
            double lo = -arma::datum::inf;
            double hi = arma::datum::inf;
            if (rr == cr) {
              lo = U_out(t);                    // beat the outside latent
              for (int kk = r0; kk < r1; ++kk) {
                if (kk != rr && U(kk) > lo) lo = U(kk);
              }
            } else if (cr < 0) {
              hi = U_out(t);                    // outside chosen: below it
            } else {
              hi = U(cr);                       // below the chosen row
            }
            U(rr) = rtruncnorm(rng_u, m, sd_eps, lo, hi);
          }

          // Outside latent: mean 0, full participant in the augmentation.
          {
            double lo = -arma::datum::inf;
            double hi = arma::datum::inf;
            if (cr < 0) {
              lo = -arma::datum::inf;
              for (int kk = r0; kk < r1; ++kk) {
                if (U(kk) > lo) lo = U(kk);
              }
              if (r1 == r0) lo = -arma::datum::inf;  // defensive (M >= 1)
            } else {
              hi = U(cr);
            }
            U_out(t) = rtruncnorm(rng_u, 0.0, sd_eps, lo, hi);
          }
        }

        // beta_i | U_i, b, W, s2, delta: conjugate regression on the
        // delta-residualized utilities. Precision W^{-1} + G_i / s2 (worker
        // thread: hand-rolled Cholesky only).
        Xoshiro256pp rng_b = make_stream(useed, static_cast<uint64_t>(r),
                                         static_cast<uint64_t>(N + i));
        const arma::mat& Gi = G.slice(i);
        for (int cj = 0; cj < K; ++cj) {
          for (int ci = 0; ci < K; ++ci) {
            Q_loc(ci, cj) = Winv(ci, cj) + Gi(ci, cj) * inv_s2;
          }
        }
        for (int kk = 0; kk < K; ++kk) {         // rhs = W^{-1} b
          double s = 0.0;
          for (int ll = 0; ll < K; ++ll) s += Winv(kk, ll) * b(ll);
          rhs_loc(kk) = s;
        }
        for (int t = t0; t < t1; ++t) {          // + X_i'(U - delta)/s2
          for (int rr = panel.row_offsets[t]; rr < panel.row_offsets[t + 1];
               ++rr) {
            const double resid =
              (U(rr) - delta(panel.alt_of_row[rr])) * inv_s2;
            for (int kk = 0; kk < K; ++kk) rhs_loc(kk) += X(rr, kk) * resid;
          }
        }
        if (!hb_chol_lower(Q_loc, L_loc)) {
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
        hb_chol_solve(L_loc, rhs_loc, mean_loc, t_loc);
        hb_mvn_precision_draw(rng_b, mean_loc, L_loc, bnew, z_loc, t_loc);
        for (int kk = 0; kk < K; ++kk) beta_i(kk, i) = bnew(kk);

        // Refresh the row cache with the new beta_i (consumed by the delta
        // and RSS passes).
        for (int t = t0; t < t1; ++t) {
          for (int rr = panel.row_offsets[t]; rr < panel.row_offsets[t + 1];
               ++rr) {
            double s = 0.0;
            for (int kk = 0; kk < K; ++kk) s += X(rr, kk) * bnew(kk);
            mu_x(rr) = s;
          }
        }
      }   // implied barrier: U, beta_i, mu_x complete
      if (abort_code != ABORT_NONE) break;

      // --- (b) delta_j | U, beta, theta, sigma_d2, s2 (work-shared) ---------
      // Conditionally independent across j given U; each j's sufficient
      // statistic is summed over its CSR row list in fixed row order by one
      // thread, so the draw is bitwise thread-count invariant.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int j = 0; j < J; ++j) {
        Xoshiro256pp rng_d = make_stream(
          useed, static_cast<uint64_t>(r),
          static_cast<uint64_t>(2LL * N + j));
        double resid_sum = 0.0;
        const int i0 = panel.alt_row_offsets[j];
        const int i1 = panel.alt_row_offsets[j + 1];
        for (int idx = i0; idx < i1; ++idx) {
          const int rr = panel.alt_rows[idx];
          resid_sum += U(rr) - mu_x(rr);
        }
        const double n_j = static_cast<double>(i1 - i0);
        double zt = 0.0;                          // prior mean z_j' theta
        for (int p = 0; p < P; ++p) zt += Z(j, p) * theta(p);
        const double prec = n_j * inv_s2 + 1.0 / sigma_d2;
        const double mean = (resid_sum * inv_s2 + zt / sigma_d2) / prec;
        delta(j) = mean + rng_d.rnorm() / std::sqrt(prec);
      }   // implied barrier: delta complete

      // --- (c) RSS with the NEW beta and delta (work-shared, fixed order) ---
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        double rss = 0.0;
        for (int t = panel.task_offsets[i]; t < panel.task_offsets[i + 1];
             ++t) {
          for (int rr = panel.row_offsets[t]; rr < panel.row_offsets[t + 1];
               ++rr) {
            const double e = U(rr) - mu_x(rr) - delta(panel.alt_of_row[rr]);
            rss += e * e;
          }
          const double eo = U_out(t);             // outside residual = latent
          rss += eo * eo;
        }
        rss_slot(i) = rss;
      }   // implied barrier

      // --- (d) Hierarchy + variance draws + recording (master) --------------
#ifdef _OPENMP
#pragma omp master
#endif
      {
        {
          Xoshiro256pp rng_bb = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 0));
          if (!draw_b_conditional(rng_bb, beta_i, b_bar, A, Winv, b)) {
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
          // sigma_d^2 then s2 share the last tag (documented partition):
          // one stream, two draws in fixed order.
          Xoshiro256pp rng_s = make_stream(
            useed, static_cast<uint64_t>(r),
            static_cast<uint64_t>(2LL * N + J + 3));
          for (int j = 0; j < J; ++j) {           // xi = delta - Z theta
            double m = 0.0;
            for (int p = 0; p < P; ++p) m += Z(j, p) * theta(p);
            xi_id_scratch(j) = delta(j) - m;      // raw-scale xi (scratch)
          }
          draw_sigma_d2_conditional(rng_s, xi_id_scratch, sdp, sigma_d2, a_d);

          double rss_total = 0.0;                 // fixed respondent order
          for (int i = 0; i < N; ++i) rss_total += rss_slot(i);
          s2 = hb_rinvgamma(rng_s, a0 + 0.5 * n_latents,
                            s0 + 0.5 * rss_total);
        }

        // --- Recording (raw chains; identified-scale summaries) ------------
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
          sigma2draw(row) = s2;

          // Identified-scale accumulation: divide by the CURRENT sigma so
          // summaries are E[beta_i / sigma], never E[beta_i] / E[sigma].
          const double inv_sig = 1.0 / std::sqrt(s2);
          if (keep_beta_i >= 1) {
            for (int i = 0; i < N; ++i) {
              for (int k = 0; k < K; ++k) {
                beta_id_scratch(k, i) = beta_i(k, i) * inv_sig;
              }
            }
            beta_wf.update(beta_id_scratch);
            if (keep_beta_i == 2) beta_draws.slice(row) = beta_id_scratch;
          }
          for (int j = 0; j < J; ++j) {
            delta_id_scratch(j) = delta(j) * inv_sig;
            xi_id_scratch(j) *= inv_sig;          // xi computed above
          }
          delta_wf.update(delta_id_scratch);
          xi_wf.update(xi_id_scratch);
        }
        if (abort_code == ABORT_NONE) {
          if (trace > 0 && (r + 1) % trace == 0) {
            Rprintf("hmnp_gibbs: iteration %d / %d\n", r + 1, R);
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
      Rcpp::stop("hmnp_gibbs: posterior precision of beta_i is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_B:
      Rcpp::stop("hmnp_gibbs: posterior precision of b is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_W:
      Rcpp::stop("hmnp_gibbs: W scale matrix is not positive definite at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_WINV:
      Rcpp::stop("hmnp_gibbs: W draw is numerically singular at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_THETA:
      Rcpp::stop("hmnp_gibbs: posterior precision of theta is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_INTERRUPT:
      Rcpp::stop("hmnp_gibbs: interrupted by user at iteration %d.",
                 abort_iter + 1);
    default:
      break;
  }

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
    Rcpp::Named("sigma2draw") = sigma2draw,
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
