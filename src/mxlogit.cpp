// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "halton.h"

// Reconstruct lower-triangular choleski factor L from L_params
// [[Rcpp::export]]
arma::mat build_L_mat(const arma::vec &L_params, const int K_w,
                      const bool rc_correlation) {
  arma::mat L = arma::zeros(K_w, K_w);
  int idx = 0;
  if (rc_correlation) { // full (lower-triangular) factor
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx) {
        double val = L_params(idx);
        if (i == j) { // diagonal - exp()
          L(i, j) = std::exp(val);
        } else {
          L(i, j) = val; // off-diagonal stays unconstrained
        }
      }
    }
  } else { // diagonal-only (Sigma is diagonal)
    for (int k = 0; k < K_w; ++k) {
      L(k, k) = std::exp(L_params(k));
    }
  }
  return L;
}

//' Reconstruct variance matrix L from L_params
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @returns matrix equal to LL', where L is the choleski decomposition of random coefficient matrix
//' @examples
//' L_params <- c(log(1.0), 0.3, log(0.5))
//' Sigma <- choicer:::build_var_mat(L_params, K_w = 2, rc_correlation = TRUE)
//' Sigma  # 2x2 covariance matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat build_var_mat(const arma::vec &L_params, const int K_w,
                        const bool rc_correlation) {
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  // Return the variance matrix
  return L * L.t();
}

// ============================================================================
// Likelihood units: the panel structure shared by the estimation kernels
//
// A likelihood unit is one decision maker. Unit u owns the contiguous choice
// situations t = off[u], ..., off[u+1]-1 (situations are sorted by decision
// maker) and ONE K_w x S block of draws eta_u, shared by all of them. With
// Ti = NULL every situation is its own unit (off = 0, 1, ..., N): the
// cross-sectional likelihood, with the same draw block and weight per
// situation as before. With beta_us = mu_final + Gamma_us the tastes at draw
// s and j_t the alternative chosen in situation t:
//
//   lambda_us = sum_t log P_ts(j_t)               log-probability of u's choices
//   ell_u     = log (1/S) sum_s exp(lambda_us)    simulated log-likelihood
//   omega_us  = exp(lambda_us - LSE_s lambda_u.)  posterior draw weights (sum 1)
//   s_u       = sum_s omega_us sum_t g_uts        score, where
//   g_uts     = sum_a (1{a = j_t} - P_ts(a)) z_tsa,  z_tsa = dV_tsa / dtheta
//
// Draw weights stay in log space (lambda sums, max-shifted LSE), so products
// of many choice probabilities never underflow.
// ============================================================================

// Unit offsets (CSR, half-open) from Ti; NULL gives the cross-section.
// Master thread only (reads an R vector); workers see plain ints. The Ti
// checks mirror hmnl_gibbs (src/hmnlogit.cpp): at least one respondent,
// every Ti positive, and sum(Ti) equal to the number of choice situations.
inline std::vector<int> mxl_unit_offsets(
    const Rcpp::Nullable<Rcpp::IntegerVector>& Ti, const int N) {
  if (Ti.isNull()) { // unit u = choice situation u
    std::vector<int> off(N + 1);
    std::iota(off.begin(), off.end(), 0);
    return off;
  }
  const Rcpp::IntegerVector T_u(Ti.get());
  if (T_u.size() == 0) {
    Rcpp::stop("Ti must contain at least one respondent.");
  }
  long long total = 0;
  for (int u = 0; u < T_u.size(); ++u) {
    if (T_u[u] == NA_INTEGER || T_u[u] < 1) {
      Rcpp::stop("Ti must be positive for every respondent (Ti[%d] = %s).",
                 u + 1, T_u[u] == NA_INTEGER ? "NA" : std::to_string(T_u[u]));
    }
    total += T_u[u];
  }
  if (total != N) {
    Rcpp::stop("sum(Ti) (%d) does not match the number of choice situations "
               "(%d).", total, N);
  }
  std::vector<int> off(T_u.size() + 1, 0);
  std::partial_sum(T_u.begin(), T_u.end(), off.begin() + 1);
  return off;
}

// Weights are decision-maker weights (objective sum_u w_u ell_u): constant
// within each unit, whose weight is weights[off[u]]. O(N), master thread; a
// no-op in the cross-section. NaN-aware: two NaN weights count as equal, so
// only genuinely varying weights are reported (finiteness is checked in R).
inline void check_unit_weights(const arma::vec& weights,
                               const std::vector<int>& off) {
  for (std::size_t u = 0; u + 1 < off.size(); ++u) {
    const double w_u = weights[off[u]];
    for (int t = off[u] + 1; t < off[u + 1]; ++t) {
      const double w_t = weights[t];
      if (w_t != w_u && !(std::isnan(w_t) && std::isnan(w_u))) {
        Rcpp::stop("weights must be constant within each decision maker "
                   "(unit %d).", u + 1);
      }
    }
  }
}

// Per-thread scratch grows with the stacked rows R_u of the largest unit:
// WGamma and DiffW are R_u x S doubles each. An allocation failure inside
// the parallel region would terminate R, so refuse on the primary thread
// when the two would exceed 2 GiB, i.e. max_u R_u x S > 2^27 (134 million);
// in the cross-section R_u = M_i, far below that for realistic data.
inline void check_unit_scratch(const std::vector<int>& row_off,
                               const std::vector<int>& off, const int S) {
  int max_R = 0;
  for (std::size_t u = 0; u + 1 < off.size(); ++u) {
    max_R = std::max(max_R, row_off[off[u + 1]] - row_off[off[u]]);
  }
  const double bytes = 2.0 * sizeof(double) * max_R * static_cast<double>(S);
  if (bytes > 2.0 * 1024 * 1024 * 1024) {
    Rcpp::stop("The largest decision maker stacks %d alternative rows across "
               "%d draws, which needs %.1f GB of scratch memory per thread. "
               "Reduce S, or check person_col: it should identify decision "
               "makers, not markets.", max_R, S, bytes / 1e9);
  }
}

// Read-only inputs of the per-unit routines, shared by all threads: the
// stacked data, the model at theta, the draw source and the unit layout.
// Built on the master thread; holds no SEXP.
struct MxlUnitData {
  const arma::mat& X;
  const arma::mat& W;              // row-aligned with X, or J x K_w
  const arma::uvec& alt_idx0;      // 0-based alternative of each stacked row
  const arma::uvec& choice_idx;    // per situation, as passed to the kernel
  const arma::vec& base_util;      // X beta + W mu_final + delta, per row
  const std::vector<int>& row_off; // situation t: rows [row_off[t], row_off[t+1])
  const std::vector<int>& off;     // unit u: situations [off[u], off[u+1])
  const arma::cube& eta_draws;     // store mode: slice u
  const HaltonGen& gen;            // generate mode: Halton block u + 1
  const MxlParams& par;
  const arma::uvec& rc_dist;
  const int n_params;
  const bool use_generate, rc_correlation, rc_mean, use_asc,
      include_outside_option;
};

// Thread-private state of the unit in hand; declare it inside the parallel
// region. Members are resized per unit and keep their memory across units.
struct MxlUnitScratch {
  int t0 = 0, t1 = 0;   // situations [t0, t1)
  int r0 = 0, R = 0;    // stacked rows [r0, r0 + R)
  arma::mat eta;        // K_w x S draws of the unit
  arma::mat W_u;        // R x K_w
  arma::mat Gamma;      // K_w x S random-coefficient draws (Gamma_final)
  arma::mat Dgamma1;    // K_w x S first derivative of the RC transform
  arma::mat Dgamma2;    // K_w x S second derivative (Hessian only)
  arma::mat WGamma;     // R x S: W_u * Gamma
  arma::vec util, V, P; // one situation at one draw
  arma::vec lambda;     // S: sum_t log P_ts(j_t)
  arma::vec omega;      // S: posterior draw weights
  arma::mat DiffW;      // R x S unweighted residuals 1{r = j_t} - P_ts(r)
  arma::vec d_bar;      // R: DiffW * omega
  arma::mat BW, A;      // K_w x S and K_w x K_w: Cholesky-block collapse
};

// Load unit u: W_u, the draws eta_u (cube slice u or on-the-fly Halton block
// u + 1), Gamma_u = batch_gamma_draws(L, eta_u) with its derivatives, and
// WGamma_u = W_u Gamma_u in a single dgemm.
inline void mxl_unit_load(const MxlUnitData& ud, const int u,
                          MxlUnitScratch& sc, const bool with_Dgamma2 = false) {
  sc.t0 = ud.off[u];
  sc.t1 = ud.off[u + 1];
  sc.r0 = ud.row_off[sc.t0];
  sc.R = ud.row_off[sc.t1] - sc.r0;
  const int r1 = sc.r0 + sc.R - 1;
  sc.W_u = make_W_i(ud.W, ud.X.n_rows, sc.r0, r1,
                    ud.alt_idx0.subvec(sc.r0, r1));
  if (ud.use_generate) {
    ud.gen.fill_eta_i(sc.eta, u + 1);
  } else {
    sc.eta = ud.eta_draws.slice(u);
  }
  sc.Gamma = batch_gamma_draws(ud.par.L, sc.eta, ud.rc_dist, &sc.Dgamma1,
                               with_Dgamma2 ? &sc.Dgamma2 : nullptr);
  sc.WGamma = sc.W_u * sc.Gamma;
}

// Logit probabilities of situation t at draw s, into sc.P: inside utilities
// base_util + WGamma(rows of t, s), with the outside option (V = 0) in slot 0
// when present.
inline void mxl_situation_probs(const MxlUnitData& ud, MxlUnitScratch& sc,
                                const int t, const int s) {
  const int first = ud.row_off[t], last = ud.row_off[t + 1] - 1;
  const int num_choices =
      ud.include_outside_option ? last - first + 2 : last - first + 1;
  sc.util = ud.base_util.subvec(first, last) +
            sc.WGamma.col(s).subvec(first - sc.r0, last - sc.r0);
  sc.V.set_size(num_choices);
  fill_choice_utilities(sc.V, sc.util, num_choices, ud.include_outside_option);
  stable_softmax(sc.V, sc.P);
}

// Draw loop of the loaded unit: lambda_s = sum_t log P_ts(j_t), the same
// log(P_choice) primitive as the cross-section, and, when `residuals`, the
// UNWEIGHTED residuals DiffW (R x S): over the rows of situation t, column s
// holds 1{r = j_t} - P_ts(r) for the inside alternatives (the outside option
// carries no parameter, so its slot is dropped). Fills
// omega_s = exp(lambda_s - lse) and returns lse = log sum_s exp(lambda_s).
inline double mxl_unit_simulate(const MxlUnitData& ud, MxlUnitScratch& sc,
                                const bool residuals) {
  const int S = sc.Gamma.n_cols;
  const int o = ud.include_outside_option ? 1 : 0; // slot of 1st inside alt
  sc.lambda.zeros(S);
  if (residuals) sc.DiffW.set_size(sc.R, S);
  for (int t = sc.t0; t < sc.t1; ++t) {
    const int r = ud.row_off[t] - sc.r0;             // first row of t in unit
    const int m = ud.row_off[t + 1] - ud.row_off[t]; // inside alternatives
    int chosen_alt = static_cast<int>(ud.choice_idx[t]); // validated serially
    if (!ud.include_outside_option) chosen_alt -= 1;      // slot in P
    for (int s = 0; s < S; ++s) {
      mxl_situation_probs(ud, sc, t, s);
      sc.lambda(s) += std::log(sc.P(chosen_alt));
      if (residuals) {
        sc.DiffW.col(s).subvec(r, r + m - 1) = -sc.P.tail(m);
        if (chosen_alt >= o) sc.DiffW(r + chosen_alt - o, s) += 1.0;
      }
    }
  }
  const double lse = logSumExp(sc.lambda);
  sc.omega = arma::exp(sc.lambda - lse);
  return lse;
}

// Score of the loaded unit, s_u = sum_s omega_s sum_t g_uts, by the BLAS-3
// collapse of the residuals. Utilities are linear in beta, mu and the ASC
// dummies, so those blocks need only d_bar = DiffW omega; the Cholesky block
// couples the residuals of draw s with eta_s:
//   beta : X_u' d_bar
//   mu   : (W_u' d_bar) % dmu_final_dmu                          (rc_mean)
//   L    : A = ((W_u' DiffW % Dgamma1) diag(omega)) eta_u',  s[L_pq] = dL_pq A(p, q)
//   delta: d_bar scattered by alternative over the stacked rows
inline void mxl_unit_score(const MxlUnitData& ud, MxlUnitScratch& sc,
                           arma::vec& score) {
  const MxlParams& par = ud.par;
  const int K_w = ud.W.n_cols;
  const int r1 = sc.r0 + sc.R - 1;
  score.zeros(ud.n_params);
  sc.d_bar = sc.DiffW * sc.omega; // R x 1, one dgemv

  // Beta block
  score.subvec(par.idx_beta_start, par.idx_mu_start - 1) =
      ud.X.rows(sc.r0, r1).t() * sc.d_bar;

  if (K_w > 0) {
    // Mu block (only if rc_mean)
    if (ud.rc_mean) {
      score.subvec(par.idx_mu_start, par.idx_L_start - 1) =
          (sc.W_u.t() * sc.d_bar) % par.dmu_final_dmu;
    }

    // L block: the only block with per-draw eta coupling
    sc.BW = sc.W_u.t() * sc.DiffW;     // K_w x S, one dgemm
    sc.BW %= sc.Dgamma1;               // Dgamma1 = 1 for normal rows
    sc.BW.each_row() %= sc.omega.t();  // draw weights
    sc.A = sc.BW * sc.eta.t();         // K_w x K_w, one dgemm
    if (ud.rc_correlation) {
      int lp = 0;
      for (int p = 0; p < K_w; ++p) {
        for (int q = 0; q <= p; ++q, ++lp) {
          const double dLpq = (p == q) ? par.L(p, p) : 1.0;
          score[par.idx_L_start + lp] = dLpq * sc.A(p, q);
        }
      }
    } else {
      // Diagonal L: s[L_p] = L(p,p) * A(p,p)
      for (int p = 0; p < K_w; ++p) {
        score[par.idx_L_start + p] = par.L(p, p) * sc.A(p, p);
      }
    }
  }

  // Delta block (scatter -- irregular alt-index mapping)
  if (ud.use_asc) {
    scatter_delta_grad(score, par.idx_delta_start, sc.d_bar,
                       ud.alt_idx0.subvec(sc.r0, r1), sc.R,
                       ud.include_outside_option, 1.0);
  }
}

//' Log-likelihood and gradient for Mixed Logit
//'
//' Computes the log-likelihood and its gradient for the Mixed Logit model using
//' OpenMP for parallelization. Allows for inclusion of alternative-specific
//' constants, outside option, observation weights, correlated random
//' coefficients, and panel data (one draw block per decision maker, see Ti).
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients. If so, mean parameters (mu) should be included in theta after beta parameters.
//' @param use_asc whether to use alternative-specific constants. If so, parameters should be included in theta after beta and L (and mu, if applicable).
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @returns List with loglikelihood and gradient evaluated at input arguments
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' K_x <- ncol(d$X); K_w <- ncol(d$W); J <- nrow(d$alt_mapping)
//' theta <- rep(0, K_x + K_w + J - 1)
//' result <- choicer:::mxl_loglik_gradient_parallel(theta, d$X, d$W, d$alt_idx,
//'   d$choice_idx, d$M, d$weights, eta, rc_dist = rep(0L, K_w),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' result$objective
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mxl_loglik_gradient_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue) {

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  const int n_params = theta.n_elem;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL)
  const std::vector<int> off = mxl_unit_offsets(Ti, N);
  const int n_units = static_cast<int>(off.size()) - 1;
  // In generate mode bypass the cube check; otherwise validate normally.
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights, &choice_idx, Ti.isNull() ? -1 : n_units);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights, &choice_idx);
    check_rc_dist_length(rc_dist, K_w);
  }
  check_unit_weights(weights, off);

  // Convenience objects shared by all threads
  arma::uvec alt_idx0 = alt_idx - 1; // 0-based
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);
  const std::vector<int> row_off(S_prefix.begin(), S_prefix.end());
  check_unit_scratch(row_off, off, Sdraw); // primary thread, before any scratch

  // Pre-compute base utility for all individuals (single BLAS call)
  // base_util = X*beta + W*mu_final + delta, computed once for all rows
  arma::vec base_util = compute_base_util_mxl(X, W, par.beta, par.mu_final,
                                              alt_idx0, use_asc, par.delta);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  // Rcpp::stop() is only safe outside parallel regions.
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d (mxl_loglik_gradient_parallel)", i);
    }
  }

  // Construct on-the-fly Halton generator (outside parallel region; no shared mutable state).
  // In cube mode (gen_seed < 0) this object is never used.
  const bool use_generate = (gen_seed >= 0);
  HaltonGen halton_gen;
  if (use_generate) {
    halton_gen = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const MxlUnitData ud{X, W, alt_idx0, choice_idx, base_util, row_off, off,
                       eta_draws, halton_gen, par, rc_dist, n_params,
                       use_generate, rc_correlation, rc_mean, use_asc,
                       include_outside_option};
  const double log_S = std::log(static_cast<double>(Sdraw));

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators and unit scratch
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int u = 0; u < n_units; ++u) {
      mxl_unit_load(ud, u, sc);
      const double lse = mxl_unit_simulate(ud, sc, true); // log sum_s exp(lambda_s)
      mxl_unit_score(ud, sc, s_u);

      // ell_u = lse - log(S); unit weight = weight of its first situation
      const double w_u = weights[off[u]];
      local_loglik += w_u * (lse - log_S);
      local_grad += w_u * s_u;
    } // end unit loop

// Combine thread-local results
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
    }
  } // end parallel region

  // Sanitize NaN/Inf so the optimizer's line-search can backtrack instead
  // of stalling at an undefined objective. The minimizer expects a finite
  // reference value; a "very bad" sentinel lets it shrink the step.
  double obj = -global_loglik;
  arma::vec grad = -global_grad;
  if (!std::isfinite(obj)) {
    obj = 1e10;
    grad.zeros();
  } else {
    grad.elem(arma::find_nonfinite(grad)).zeros();
  }
  return Rcpp::List::create(Rcpp::Named("objective") = obj,
                            Rcpp::Named("gradient") = grad);
}

// vech(): lower-triangular vectorisation (including the diagonal), row-major
// Matches the row-major packing of L_params used by build_L_mat() above and
// by the outer parameter loop in jacobian_vech_Sigma() below.
inline arma::vec vech(const arma::mat &M) {
  arma::uword K = M.n_rows;
  arma::vec out(K * (K + 1) / 2);
  arma::uword idx = 0;
  for (arma::uword i = 0; i < K; ++i)
    for (arma::uword j = 0; j <= i; ++j)
      out(idx++) = M(i, j);
  return out;
}

//' Utility to compute analytical Jacobian of random coefficient matrix transformed by vech (dVech(Sigma) / dTheta)
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @returns Jacobian (dVech(Sigma) / dTheta)
//' @examples
//' L_params <- c(log(0.8), 0.2, log(0.6))
//' J_mat <- choicer:::jacobian_vech_Sigma(L_params, K_w = 2, rc_correlation = TRUE)
//' dim(J_mat)  # 3 x 3 for K_w=2 correlated
//' @keywords internal
// [[Rcpp::export]]
arma::mat jacobian_vech_Sigma(const arma::vec &L_params, const int K_w,
                              const bool rc_correlation = true) {
  // dimensions
  const int L_size = rc_correlation ? K_w * (K_w + 1) / 2 : K_w;

  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  arma::mat J(L_size, L_size, arma::fill::zeros);

  // loop over parameters
  arma::mat E(K_w, K_w, arma::fill::zeros); // holds dL / dtheta_m
  std::size_t idx_param = 0;

  if (rc_correlation) {
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx_param) {
        // reset E
        E.zeros();
        if (i == j) {        // diagonal: L_ii = exp(z_i)
          E(i, j) = L(i, i); // dL_ii / dz_i = exp(z_i)
        } else {             // off-diagonal parameter
          E(i, j) = 1.0;
        }
        arma::mat dSigma = E * L.t() + L * E.t(); // product rule
        J.col(idx_param) = vech(dSigma);
      }
    }

  } else {
    // diagonal Sigma only (no correlations)
    // Jacobian of Sigma wrt L_params  (diagonal-only case)
    for (int k = 0; k < K_w; ++k) {
      // Sigma_kk = L_kk^2  ,  L_kk = exp(z_k)
      // dSigma_kk/dz_k = 2 * exp(2 z_k) = 2 * L_kk^2
      double deriv = 2.0 * L(k, k) * L(k, k);
      J(k, k) = deriv;
    }
  }
  return J;
}

//' Analytical Hessian of the log-likelihood v2
//'
//' Computes the Hessian of the log-likelihood for the Mixed Logit model using
//' OpenMP for parallelization. Mirrors the parameters of mxl_loglik_gradient_parallel.
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @returns Hessian evaluated at input arguments
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' theta <- rep(0, ncol(d$X) + ncol(d$W) + nrow(d$alt_mapping) - 1)
//' H <- choicer:::mxl_hessian_parallel(theta, d$X, d$W, d$alt_idx, d$choice_idx,
//'   d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_hessian_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue) {
  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  const int n_params = theta.n_elem;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  // Reuse the parser's layout when assembling derivatives.
  const int idx_beta_start = par.idx_beta_start;
  const int idx_mu_start = par.idx_mu_start;
  const int idx_L_start = par.idx_L_start;
  const int idx_delta_start = par.idx_delta_start;
  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL)
  const std::vector<int> off = mxl_unit_offsets(Ti, N);
  const int n_units = static_cast<int>(off.size()) - 1;
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights, &choice_idx, Ti.isNull() ? -1 : n_units);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights, &choice_idx);
    check_rc_dist_length(rc_dist, K_w);
  }
  check_unit_weights(weights, off);
  const arma::mat& L = par.L;
  const arma::vec& dmu_final_dmu = par.dmu_final_dmu;
  const arma::vec& dmu2_final_dmu2 = par.dmu2_final_dmu2;

  arma::uvec alt_idx0 = alt_idx - 1;
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);
  const std::vector<int> row_off(S_prefix.begin(), S_prefix.end());
  check_unit_scratch(row_off, off, Sdraw); // primary thread, before any scratch

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_h = compute_base_util_mxl(X, W, par.beta, par.mu_final,
                                                alt_idx0, use_asc, par.delta);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d (mxl_hessian_parallel)", i);
    }
  }

  // Block layout: continuous block c = [beta | mu | L], size Kc = idx_delta_start
  //               delta block d = [delta ASCs],         size Jd = n_params - idx_delta_start
  // H_V (second derivative of utility w.r.t. theta) is nonzero only in the
  // (mu, L) x (mu, L) sub-block, which lives entirely within the continuous block.
  const int Kc = idx_delta_start;           // size of continuous block
  const int Jd = n_params - idx_delta_start; // size of delta block (0 when !use_asc)

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_h = (gen_seed >= 0);
  HaltonGen halton_gen_h;
  if (use_generate_h) {
    halton_gen_h = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const MxlUnitData ud{X, W, alt_idx0, choice_idx, base_util_h, row_off, off,
                       eta_draws, halton_gen_h, par, rc_dist, n_params,
                       use_generate_h, rc_correlation, rc_mean, use_asc,
                       include_outside_option};

  // Per unit u (Louis 1982), with omega_s the posterior draw weights and
  // g_s = sum_t g_ts, H_s = sum_t H_ts the per-draw score and Hessian of
  // lambda_s = sum_t log P_ts(j_t):
  //   H_u = sum_s omega_s (H_s + g_s g_s') - g_bar g_bar',  g_bar = sum_s omega_s g_s,
  //   H_ts = -sum_Pzz_ts + sum_Pz_ts sum_Pz_ts' + sum_diff_H_V_ts.
  // Pass 1 (the shared draw loop) yields omega; pass 2 accumulates the O3
  // buffers with omega_s where the cross-section used P_choice_s / P_i_hat.

  // Global accumulator
  arma::mat global_hess = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_hess = arma::zeros(n_params, n_params);

    // Thread-private scratch — sized once, reset inside loops as needed.
    MxlUnitScratch sc; // unit slices, draws, probabilities P_ts, weights omega

    // Per-(t, s) block accumulators (reset at the start of each situation-draw).
    // cc: Kc x Kc dense outer-product sum
    arma::mat sum_Pzz_cc(Kc, Kc);
    // cd: Kc x Jd; each col j accumulates P_a * zc_a for the alt whose delta is j
    // Note: the "Jd > 0 ? Jd : 1" sentinel below (and for all Jd-sized buffers) avoids
    // zero-size allocation when there are no ASC parameters; these buffers are never read
    // when Jd == 0 because every access is guarded by "if (Jd > 0)" / "if (j_delta >= 0)".
    arma::mat sum_Pzz_cd(Kc, Jd > 0 ? Jd : 1);
    // dd: diagonal only — stored as length-Jd vector
    arma::vec sum_Pzz_dd(Jd > 0 ? Jd : 1);
    // P*z sums
    arma::vec sum_Pz_c(Kc);
    arma::vec sum_Pz_d(Jd > 0 ? Jd : 1);
    // gradient components
    arma::vec g_c(Kc);
    arma::vec g_d(Jd > 0 ? Jd : 1);
    // H_V in the continuous block (mu,L sub-block; beta and delta rows/cols are zero)
    arma::mat sum_diff_H_V_cc(Kc, Kc);

    // Per-alt scratch for the continuous-block z vector
    arma::vec zc_a(Kc);

    // O3: Per-unit block buffers — accumulate the pieces linear in omega_s
    // across the unit's situations and draws. These replace the per-draw H_is
    // allocation and a running sum_s omega_s (g_s g_sT + H_s) accumulation.
    arma::mat buf_Pzz_cc(Kc, Kc);               // Identity C: Σ_ts ω_s sum_Pzz_cc_ts
    arma::mat buf_Pzz_cd(Kc, Jd > 0 ? Jd : 1); // Identity C: Σ_ts ω_s sum_Pzz_cd_ts
    arma::vec buf_Pzz_dd(Jd > 0 ? Jd : 1);      // Identity C: Σ_ts ω_s sum_Pzz_dd_ts
    arma::mat buf_diff_HV_cc(Kc, Kc);           // Identity C: Σ_ts ω_s sum_diff_H_V_cc_ts
    // O3: Column stashes for BLAS-3 batching.
    // G_stash(:,s) = sqrt(ω_s) * g_s, g_s = Σ_t g_ts → G Gᵀ = Σ_s ω_s g_s g_sᵀ (Identity A)
    // F_stash(:,s) = sqrt(ω_s) * [sum_Pz_c; sum_Pz_d]_ts → F Fᵀ = Σ_s ω_s sum_Pz sum_PzT
    //                                                     for situation t (Identity B)
    // These are allocated once per thread at max size (n_params x Sdraw).
    // F column s is fully written in draw s of every situation; G is zeroed
    // per unit and accumulates over the unit's situations.
    arma::mat G_stash(n_params, Sdraw);
    arma::mat F_stash(n_params, Sdraw);
    arma::mat opg_pz(n_params, n_params); // G Gᵀ + Σ_t F_t F_tᵀ
    arma::vec g_bar(n_params);            // Σ_s ω_s g_s: the unit's score

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int u = 0; u < n_units; ++u) {
      const double w_u = weights[off[u]]; // unit weight

      // --- Pass 1: posterior draw weights omega_s of the unit's choices ---
      mxl_unit_load(ud, u, sc, true);
      const double lse = mxl_unit_simulate(ud, sc, false);
      if (!std::isfinite(lse)) continue; // zero simulated likelihood: skip unit

      // O3: Initialize per-unit block buffers.
      buf_Pzz_cc.zeros();
      buf_diff_HV_cc.zeros();
      if (Jd > 0) {
        buf_Pzz_cd.zeros();
        buf_Pzz_dd.zeros();
      }
      G_stash.zeros();
      opg_pz.zeros();
      g_bar.zeros();

      // --- Pass 2: situations t (outer) x draws s (inner) ---
      for (int t = sc.t0; t < sc.t1; ++t) {
        const int m_t = row_off[t + 1] - row_off[t];
        const int num_choices = include_outside_option ? m_t + 1 : m_t;

        // chosen alternative index (validated serially above)
        int chosen_alt = choice_idx[t];
        if (!include_outside_option)
          chosen_alt -= 1;

        // Loop over simulations (draws) s
        for (int s = 0; s < Sdraw; ++s) {
          // Column views into the batched matrices (zero-copy)
          const auto eta_i_s                = sc.eta.col(s);
          const auto dgamma_final_dgamma    = sc.Dgamma1.col(s);
          const auto d2gamma_final_dgamma2  = sc.Dgamma2.col(s);
          const double omega_s = sc.omega(s);

          // === 3. Utility and probabilities P_ts (as in pass 1) ===
          mxl_situation_probs(ud, sc, t, s);
          const arma::vec& P_s = sc.P;

          // === 4. Per-(t, s) Gradient (g_ts) and Hessian (H_ts) — block form ===
          // Reset per-draw block accumulators.
          sum_Pzz_cc.zeros();
          sum_Pz_c.zeros();
          g_c.zeros();
          sum_diff_H_V_cc.zeros();
          if (Jd > 0) {
            sum_Pzz_cd.zeros();
            sum_Pzz_dd.zeros();
            sum_Pz_d.zeros();
            g_d.zeros();
          }

          for (int a = 0; a < num_choices; ++a) {
            // Determine if this alt contributes a nonzero z (outside option at a==0 skips)
            if (include_outside_option && a == 0) {
              // Outside option: z_as == 0 and H_V_as == 0, so this alternative
              // contributes nothing to the Hessian accumulators — skip entirely.
              continue;
            }

            const int current_a_idx = include_outside_option ? a - 1 : a;
            const int row = row_off[t] + current_a_idx; // stacked row of alt a
            const arma::rowvec w_ap_row = sc.W_u.row(row - sc.r0);

            // --- Build continuous-block vector zc_a ---
            zc_a.zeros();

            // beta sub-block
            zc_a.subvec(idx_beta_start, idx_mu_start - 1) = X.row(row).t();

            // mu sub-block  (nonzero only when rc_mean)
            if (rc_mean) {
              for (int p = 0; p < K_w; ++p) {
                zc_a[idx_mu_start + p] = w_ap_row(p) * dmu_final_dmu(p);
              }
            }

            // L sub-block
            if (rc_correlation) {
              int lp_idx_r = 0;
              for (int p = 0; p < K_w; ++p) {
                const double f_p_prime = dgamma_final_dgamma(p);
                for (int q = 0; q <= p; ++q, ++lp_idx_r) {
                  const double dLpq_dparam_r = (p == q) ? L(p, p) : 1.0;
                  zc_a[idx_L_start + lp_idx_r] =
                      w_ap_row(p) * f_p_prime * dLpq_dparam_r * eta_i_s(q);
                }
              }
            } else { // Diagonal L
              for (int p = 0; p < K_w; ++p) {
                zc_a[idx_L_start + p] =
                    w_ap_row(p) * dgamma_final_dgamma(p) * L(p, p) * eta_i_s(p);
              }
            }

            // --- Delta index for this alt (may be -1 meaning no delta entry) ---
            int j_delta = -1; // index into delta block; -1 = no entry
            if (use_asc && Jd > 0) {
              const int id = alt_idx0[row];
              if (include_outside_option) {
                j_delta = id;           // always valid (id >= 0)
              } else if (id > 0) {
                j_delta = id - 1;       // first inside alt (id==0) has no ASC
              }
              // id==0 and !include_outside_option => j_delta stays -1
            }

            // --- Accumulate block sums ---
            const double P_a = P_s(a);
            const double diff = (a == chosen_alt ? 1.0 : 0.0) - P_a;

            // cc block: P_a * zc_a * zc_a^T  (rank-1 update)
            sum_Pzz_cc += P_a * (zc_a * zc_a.t());
            sum_Pz_c   += P_a * zc_a;
            g_c        += diff * zc_a;

            // cd and dd blocks (delta scatter)
            if (j_delta >= 0) {
              sum_Pzz_cd.col(j_delta) += P_a * zc_a;
              sum_Pzz_dd(j_delta)     += P_a;
              sum_Pz_d(j_delta)       += P_a;
              g_d(j_delta)            += diff;
            }

            // H_V accumulation — nonzero only in (mu,L)x(mu,L) within cc block
            // Reuse same logic as before but write directly into sum_diff_H_V_cc
            if (rc_mean) {
              for (int p = 0; p < K_w; ++p) {
                sum_diff_H_V_cc(idx_mu_start + p, idx_mu_start + p) +=
                    diff * w_ap_row(p) * dmu2_final_dmu2(p);
              }
            }
            if (rc_correlation) {
              int lp_idx_r = 0;
              for (int p = 0; p < K_w; ++p) {
                const double f_p_prime        = dgamma_final_dgamma(p);
                const double f_p_double_prime = d2gamma_final_dgamma2(p);
                for (int q = 0; q <= p; ++q, ++lp_idx_r) {
                  const int r = idx_L_start + lp_idx_r;
                  const double dLpq_dparam_r = (p == q) ? L(p, p) : 1.0;
                  const double d2Lpq_dparam2_r = (p == q) ? L(p, p) : 0.0;
                  const double dgamma_p_dparam_r = dLpq_dparam_r * eta_i_s(q);
                  const double d2gamma_p_dparam2_r = d2Lpq_dparam2_r * eta_i_s(q);

                  // Diagonal H_V entry
                  sum_diff_H_V_cc(r, r) +=
                      diff * w_ap_row(p) *
                      (f_p_double_prime * std::pow(dgamma_p_dparam_r, 2) +
                       f_p_prime * d2gamma_p_dparam2_r);

                  // Off-diagonal L-L entries (same row p, column < q)
                  for (int q_s = 0; q_s < q; ++q_s) {
                    const int param_s = idx_L_start + lp_idx_r - (q - q_s);
                    const double dLpq_s_dparam = (p == q_s) ? L(p, p) : 1.0;
                    const double dgamma_p_dparam_s = dLpq_s_dparam * eta_i_s(q_s);
                    const double d2V = diff * w_ap_row(p) * f_p_double_prime *
                                       dgamma_p_dparam_r * dgamma_p_dparam_s;
                    sum_diff_H_V_cc(r, param_s) += d2V;
                    sum_diff_H_V_cc(param_s, r) += d2V;
                  }
                }
              }
            } else { // Diagonal L
              for (int p = 0; p < K_w; ++p) {
                const int r = idx_L_start + p;
                const double f_p_prime        = dgamma_final_dgamma(p);
                const double f_p_double_prime = d2gamma_final_dgamma2(p);
                const double dgamma_p_dparam  = L(p, p) * eta_i_s(p);
                const double d2gamma_p_dparam = L(p, p) * eta_i_s(p);
                sum_diff_H_V_cc(r, r) +=
                    diff * w_ap_row(p) *
                    (f_p_double_prime * std::pow(dgamma_p_dparam, 2) +
                     f_p_prime * d2gamma_p_dparam);
              }
            }
          } // end alt loop

          // === 5. O3: Accumulate per-unit buffers (Identity C) and fill
          //          column stashes G_stash/F_stash (Identities A + B). ===
          //
          // Instead of assembling H_ts (n_params x n_params) and accumulating
          //   hess_term1 += omega_s * (g_s g_sT + sum_t H_ts),
          // we split the linear-in-omega_s and the outer-product pieces:
          //
          //   Linear (Identity C): buf_Pzz_cc     += omega_s * sum_Pzz_cc
          //                        buf_diff_HV_cc += omega_s * sum_diff_H_V_cc
          //                        (and cd/dd variants)
          //   Outer-product (Identity A): G_stash(:,s) += g_ts, scaled by
          //     sqrt(omega_s) after the unit, so G GT = Σ_s omega_s g_s g_sT.
          //   Outer-product (Identity B): F_stash(:,s) = sqrt(omega_s) * [sum_Pz_c; sum_Pz_d]
          //     so that F FT = Σ_s omega_s sum_Pz sum_PzT after situation t.

          const double sqrt_ws = std::sqrt(omega_s);

          // Identity C — scalar-times-matrix accumulation into per-unit buffers.
          buf_Pzz_cc    += omega_s * sum_Pzz_cc;
          buf_diff_HV_cc += omega_s * sum_diff_H_V_cc;
          if (Jd > 0) {
            buf_Pzz_cd += omega_s * sum_Pzz_cd;
            buf_Pzz_dd += omega_s * sum_Pzz_dd;
          }

          // Identity A — accumulate g_s = Σ_t g_ts in column s (g_ts = [g_c; g_d]).
          G_stash.col(s).head(Kc) += g_c;
          if (Jd > 0) {
            G_stash.col(s).tail(Jd) += g_d;
          }

          // Identity B — fill column s of F_stash with sqrt(omega_s) * [sum_Pz_c; sum_Pz_d].
          F_stash.col(s).head(Kc) = sqrt_ws * sum_Pz_c;
          if (Jd > 0) {
            F_stash.col(s).tail(Jd) = sqrt_ws * sum_Pz_d;
          }

          // Accumulate the unit score g_bar = Σ_s omega_s g_s.
          g_bar.head(Kc) += omega_s * g_c;
          if (Jd > 0) {
            g_bar.tail(Jd) += omega_s * g_d;
          }
        } // end S loop

        // Identity B for situation t (one BLAS-3 product per situation: the
        // outer products do not combine across situations before the product).
        opg_pz += F_stash * F_stash.t();
      } // end situation loop

      // Identity A: G Gᵀ = Σ_s omega_s g_s g_sᵀ.
      G_stash.each_row() %= arma::sqrt(sc.omega).t();
      opg_pz += G_stash * G_stash.t();

      // === 6. O3: Per-unit finalization — assemble Hessian once from buffers ===
      // 6a. Assemble hess_term1 block-by-block.
      //     hess_term1 = (G Gᵀ + Σ_t F_t F_tᵀ) + (-buf_Pzz) + buf_diff_HV_cc (cc block only)
      //     This is the batched equivalent of Σ_s omega_s (g_s g_sT + H_s).
      arma::mat hess_t1(n_params, n_params, arma::fill::zeros);

      // cc block: contributions from OPG, sum-Pz outer product, -Pzz, and H_V.
      hess_t1.submat(0, 0, Kc - 1, Kc - 1) =
          opg_pz.submat(0, 0, Kc - 1, Kc - 1)
          - buf_Pzz_cc
          + buf_diff_HV_cc;

      if (Jd > 0) {
        // cd block: -buf_Pzz_cd + opg_pz cd block (NOT symmetric — full rectangular).
        arma::mat cd = opg_pz.submat(0, Kc, Kc - 1, n_params - 1) - buf_Pzz_cd;
        hess_t1.submat(0, Kc, Kc - 1, n_params - 1)     = cd;
        hess_t1.submat(Kc, 0, n_params - 1, Kc - 1)     = cd.t();  // dc = (cd)ᵀ

        // dd block: -diag(buf_Pzz_dd) + opg_pz dd block (sum_Pz_d outer products).
        hess_t1.submat(Kc, Kc, n_params - 1, n_params - 1) =
            opg_pz.submat(Kc, Kc, n_params - 1, n_params - 1)
            - arma::diagmat(buf_Pzz_dd);
      }

      // 6b. Louis identity: H_u = hess_term1 - g_bar g_barᵀ. The draw weights
      //     are normalized, so there is no division by the simulated P_u.
      local_hess += w_u * (hess_t1 - g_bar * g_bar.t());
    } // end unit loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_hess += local_hess;
    }
  } // end parallel region

  return -global_hess;
}

//' BHHH (outer product of gradients) information matrix for Mixed Logit
//'
//' Computes the BHHH approximation to the observed information matrix for the
//' Mixed Logit model: \eqn{H_{BHHH} = \sum_i w_i \cdot s_i s_i^\top}, where
//' \eqn{s_i} is the score of likelihood unit i (gradient of its simulated
//' log-likelihood \eqn{\log \bar{P}_i}; a decision maker when Ti is supplied,
//' a choice situation otherwise).
//' This outer product of gradients (OPG) estimator provides an alternative to
//' the analytical Hessian for standard error computation that scales to large
//' problems where the analytical Hessian is infeasible (e.g., many alternatives
//' or simulation draws).
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @returns n_params x n_params PSD matrix representing the observed information
//'   matrix estimated by the outer product of gradients (same sign convention
//'   as the negated Hessian returned by \code{mxl_hessian_parallel}, so it can
//'   be inverted directly to obtain vcov).
//' @note The BHHH/OPG estimator is only asymptotically equivalent to the
//'   Hessian-based information matrix at the true MLE. In finite samples it can
//'   underestimate standard errors, particularly when the model is mis-specified
//'   or away from the optimum.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' theta <- rep(0, ncol(d$X) + ncol(d$W) + nrow(d$alt_mapping) - 1)
//' H <- choicer:::mxl_bhhh_parallel(theta, d$X, d$W, d$alt_idx, d$choice_idx,
//'   d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_bhhh_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue) {

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  const int n_params = theta.n_elem;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL)
  const std::vector<int> off = mxl_unit_offsets(Ti, N);
  const int n_units = static_cast<int>(off.size()) - 1;
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights, &choice_idx, Ti.isNull() ? -1 : n_units);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights, &choice_idx);
    check_rc_dist_length(rc_dist, K_w);
  }
  check_unit_weights(weights, off);

  // Convenience objects shared by all threads
  arma::uvec alt_idx0 = alt_idx - 1; // 0-based
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);
  const std::vector<int> row_off(S_prefix.begin(), S_prefix.end());
  check_unit_scratch(row_off, off, Sdraw); // primary thread, before any scratch

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, par.beta, par.mu_final,
                                              alt_idx0, use_asc, par.delta);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d (mxl_bhhh_parallel)", i);
    }
  }

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_b = (gen_seed >= 0);
  HaltonGen halton_gen_b;
  if (use_generate_b) {
    halton_gen_b = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const MxlUnitData ud{X, W, alt_idx0, choice_idx, base_util, row_off, off,
                       eta_draws, halton_gen_b, par, rc_dist, n_params,
                       use_generate_b, rc_correlation, rc_mean, use_asc,
                       include_outside_option};

  // Global BHHH accumulator
  arma::mat global_bhhh = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulator and unit scratch
    arma::mat local_bhhh = arma::zeros(n_params, n_params);
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int u = 0; u < n_units; ++u) {
      mxl_unit_load(ud, u, sc);
      mxl_unit_simulate(ud, sc, true);
      mxl_unit_score(ud, sc, s_u);
      local_bhhh += weights[off[u]] * s_u * s_u.t();
    } // end unit loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_bhhh += local_bhhh;
    }
  } // end parallel region

  // Return PSD information matrix (same sign convention as negated Hessian).
  return global_bhhh;
}

// Per-unit score matrix for the mixed logit model (internal).
//
// Returns the U x n_params matrix whose row u is the weight-free score
// s_u = d ell_u / d theta of likelihood unit u (a decision maker when Ti is
// supplied, a choice situation otherwise) over the beta, mu, L and delta/ASC
// blocks. Same per-unit routine as mxl_bhhh_parallel with the outer-product
// accumulator replaced by a row write; weights are applied on the R side (see
// .assemble_score_vcov in R/classes.R).
// [[Rcpp::export]]
arma::mat mxl_scores_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue) {

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  const int n_params = theta.n_elem;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL)
  const std::vector<int> off = mxl_unit_offsets(Ti, N);
  const int n_units = static_cast<int>(off.size()) - 1;
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        nullptr, &choice_idx, Ti.isNull() ? -1 : n_units);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, nullptr, &choice_idx);
    check_rc_dist_length(rc_dist, K_w);
  }

  // Convenience objects shared by all threads
  arma::uvec alt_idx0 = alt_idx - 1; // 0-based
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);
  const std::vector<int> row_off(S_prefix.begin(), S_prefix.end());
  check_unit_scratch(row_off, off, Sdraw); // primary thread, before any scratch

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, par.beta, par.mu_final,
                                              alt_idx0, use_asc, par.delta);

  // --- Serial pre-loop validation of chosen-alternative indices ---
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d (mxl_scores_parallel)", i);
    }
  }

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_s = (gen_seed >= 0);
  HaltonGen halton_gen_s;
  if (use_generate_s) {
    halton_gen_s = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const MxlUnitData ud{X, W, alt_idx0, choice_idx, base_util, row_off, off,
                       eta_draws, halton_gen_s, par, rc_dist, n_params,
                       use_generate_s, rc_correlation, rc_mean, use_asc,
                       include_outside_option};

  // Output: one row per likelihood unit (each written by exactly one
  // iteration, so no accumulator or critical section is needed).
  arma::mat scores(n_units, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int u = 0; u < n_units; ++u) {
      mxl_unit_load(ud, u, sc);
      mxl_unit_simulate(ud, sc, true);
      mxl_unit_score(ud, sc, s_u);
      scores.row(u) = s_u.t();
    } // end unit loop
  } // end parallel region

  return scores;
}

// Conditional tastes of each likelihood unit (internal).
//
// Simulated mean and standard deviation of the random coefficients given the
// unit's observed choices (Revelt & Train 2000; Train 2009, ch. 11), on the
// scale the coefficients enter utility, beta_us = mu_final + Gamma_us (the
// log-normal rows of Gamma exponentiated). With the posterior draw weights
// omega_us of the shared draw loop:
//   mean_u = mu_final + Gamma_u omega_u
//   sd_u   = sqrt(sum_s omega_us (Gamma_us - Gamma_u omega_u)^2)   (elementwise)
// Takes the inputs of mxl_scores_parallel and returns list(mean, sd), each
// K_w x U (U = decision makers when Ti is supplied, choice situations
// otherwise). Units whose choices have zero simulated probability get NA.
// [[Rcpp::export]]
Rcpp::List mxl_conditional_tastes_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue) {

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  const int n_params = theta.n_elem;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL)
  const std::vector<int> off = mxl_unit_offsets(Ti, N);
  const int n_units = static_cast<int>(off.size()) - 1;
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        nullptr, &choice_idx, Ti.isNull() ? -1 : n_units);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, nullptr, &choice_idx);
    check_rc_dist_length(rc_dist, K_w);
  }

  // Convenience objects shared by all threads
  arma::uvec alt_idx0 = alt_idx - 1; // 0-based
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);
  const std::vector<int> row_off(S_prefix.begin(), S_prefix.end());
  check_unit_scratch(row_off, off, Sdraw); // primary thread, before any scratch

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, par.beta, par.mu_final,
                                              alt_idx0, use_asc, par.delta);

  // --- Serial pre-loop validation of chosen-alternative indices ---
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d (mxl_conditional_tastes_parallel)", i);
    }
  }

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_c = (gen_seed >= 0);
  HaltonGen halton_gen_c;
  if (use_generate_c) {
    halton_gen_c = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const MxlUnitData ud{X, W, alt_idx0, choice_idx, base_util, row_off, off,
                       eta_draws, halton_gen_c, par, rc_dist, n_params,
                       use_generate_c, rc_correlation, rc_mean, use_asc,
                       include_outside_option};
  const double na = NA_REAL; // read on the master thread

  // Output: one column per likelihood unit (disjoint writes, no reduction).
  arma::mat taste_mean(K_w, n_units);
  arma::mat taste_sd(K_w, n_units);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    MxlUnitScratch sc;
    arma::vec gamma_bar; // Gamma_u omega_u: conditional mean of Gamma

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int u = 0; u < n_units; ++u) {
      mxl_unit_load(ud, u, sc);
      const double lse = mxl_unit_simulate(ud, sc, false);
      if (!std::isfinite(lse)) {
        taste_mean.col(u).fill(na);
        taste_sd.col(u).fill(na);
        continue;
      }
      gamma_bar = sc.Gamma * sc.omega;
      taste_mean.col(u) = par.mu_final + gamma_bar;
      taste_sd.col(u) =
          arma::sqrt(arma::square(sc.Gamma.each_col() - gamma_bar) * sc.omega);
    } // end unit loop
  } // end parallel region

  return Rcpp::List::create(Rcpp::Named("mean") = taste_mean,
                            Rcpp::Named("sd") = taste_sd);
}

// ============================================================================
// Mixed Logit: Share Prediction and BLP Contraction
// ============================================================================

// Internal function for computing simulated market shares
arma::vec mxl_predict_shares_internal(
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu_final,         // Transformed mu (exp(mu) for log-normal)
    const arma::mat& L,                // Cholesky factor
    const arma::uvec& alt_idx0,        // 0-based indexing
    const Rcpp::IntegerVector& M,
    const Rcpp::IntegerVector& S_prefix,
    const arma::vec& weights,
    const arma::vec& delta,            // Full J-element delta
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int num_alts,
    const bool use_asc,
    const bool include_outside_option,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  const int N = M.size();
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  // C++ runtime guards for generate mode
  if (gen_seed >= 0 && gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
  if (gen_seed >= 0 && K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");

  const bool use_generate_s = (gen_seed >= 0);
  HaltonGen halton_gen_s;
  if (use_generate_s) {
    halton_gen_s = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const double weight_sum = arma::accu(weights);

  if (weight_sum <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_s = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Initialize global accumulator for predicted shares
  arma::vec global_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulator
    arma::vec local_shares = arma::zeros(num_alts);
    arma::mat eta_i_buf_s;    // generate mode scratch
    arma::mat eta_i_store_s;  // cube mode materialised slice

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_s.subvec(start_idx, end_idx);

      // Accumulate probabilities over draws
      arma::vec P_bar_i = arma::zeros(num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_s;
      if (use_generate_s) {
        halton_gen_s.fill_eta_i(eta_i_buf_s, i + 1);
        eta_i_ptr_s = &eta_i_buf_s;
      } else {
        eta_i_store_s = eta_draws.slice(i);
        eta_i_ptr_s = &eta_i_store_s;
      }
      const arma::mat& eta_i_s_ref = *eta_i_ptr_s;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_s_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      for (int s = 0; s < Sdraw; ++s) {
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Compute probabilities with numerical stability
        stable_softmax(V_s, P_s);

        P_bar_i += P_s;
      }  // end s loop

      P_bar_i /= static_cast<double>(Sdraw);

      // Accumulate shares by alternative
      if (include_outside_option) {
        local_shares(0) += w_i * P_bar_i(0);
      }
      for (int a = 0; a < m_i; ++a) {
        if (include_outside_option) {
          local_shares(alt_idx0_i(a) + 1) += w_i * P_bar_i(a + 1);
        } else {
          local_shares(alt_idx0_i(a)) += w_i * P_bar_i(a);
        }
      }
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_shares += local_shares;
    }
  }  // end parallel region

  return global_shares / weight_sum;
}

//' Per-observation simulated choice probabilities for Mixed Logit
//'
//' Returns the simulated choice probability for each (individual, alternative)
//' row of `X`, averaged over the supplied Halton draws. Mirrors `mnl_predict`.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether the outside option is present
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns List with `choice_prob` (length sum(M)), `utility` (length sum(M),
//'   simulated mean of the deterministic + W*gamma component), and, when
//'   `include_outside_option = TRUE`, `choice_prob_outside` (length N).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mxl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based alt indices and prefix sums
  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_p = (gen_seed >= 0);
  HaltonGen halton_gen_p;
  if (use_generate_p) {
    halton_gen_p = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Output accumulators (each individual writes to a disjoint subvec)
  arma::vec choice_prob = arma::zeros(X.n_rows);
  arma::vec utility = arma::zeros(X.n_rows);
  arma::vec choice_prob_outside;
  if (include_outside_option) {
    choice_prob_outside = arma::zeros(N);
  }

  // Thread-private buffers (declared outside the parallel loop for OpenMP)
  arma::mat eta_i_buf_p;
  arma::mat eta_i_store_p;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) firstprivate(eta_i_buf_p, eta_i_store_p)
#endif
  for (int i = 0; i < N; ++i) {
    const int m_i = M[i];
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const int start_idx = S_prefix[i];
    const int end_idx = start_idx + m_i - 1;
    const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

    arma::mat W_i =
        make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

    // Pre-computed base utility for this individual
    const arma::vec base_util_i = base_util.subvec(start_idx, end_idx);

    // Per-individual accumulators (averaged over draws)
    arma::vec P_inside_avg = arma::zeros(m_i);
    arma::vec util_inside_avg = arma::zeros(m_i);
    double P_outside_avg = 0.0;

    // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
    const arma::mat* eta_i_ptr_p;
    if (use_generate_p) {
      halton_gen_p.fill_eta_i(eta_i_buf_p, i + 1);
      eta_i_ptr_p = &eta_i_buf_p;
    } else {
      eta_i_store_p = eta_draws.slice(i);
      eta_i_ptr_p = &eta_i_store_p;
    }
    const arma::mat& eta_i_p_ref = *eta_i_ptr_p;
    arma::mat Gamma_final = batch_gamma_draws(L, eta_i_p_ref, rc_dist);

    // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
    const arma::mat WGamma = W_i * Gamma_final;

    arma::vec inside_utils(m_i);
    arma::vec V_s(num_choices);
    arma::vec P_s;

    for (int s = 0; s < Sdraw; ++s) {
      inside_utils = base_util_i + WGamma.col(s);

      fill_choice_utilities(V_s, inside_utils, num_choices,
                            include_outside_option);

      // Stable softmax
      stable_softmax(V_s, P_s);

      // Accumulate inside probabilities and utilities
      if (include_outside_option) {
        P_outside_avg += P_s(0);
        P_inside_avg += P_s.subvec(1, num_choices - 1);
      } else {
        P_inside_avg += P_s;
      }
      util_inside_avg += inside_utils;
    }

    // Average over draws
    const double S_d = static_cast<double>(Sdraw);
    P_inside_avg /= S_d;
    util_inside_avg /= S_d;
    if (include_outside_option) {
      P_outside_avg /= S_d;
      choice_prob_outside(i) = P_outside_avg;
    }

    // Disjoint writes by individual — no race
    choice_prob.subvec(start_idx, end_idx) = P_inside_avg;
    utility.subvec(start_idx, end_idx) = util_inside_avg;
  }

  Rcpp::List out;
  out["choice_prob"] = choice_prob;
  out["utility"] = utility;
  if (include_outside_option) {
    out["choice_prob_outside"] = choice_prob_outside;
  }
  return out;
}

//' Simulated expected logsum (inclusive value) for Mixed Logit
//'
//' Computes the simulated expected logsum (expected maximum utility, up to an
//' additive constant) for each choice situation:
//' \deqn{logsum_i = (1/S) \sum_s \log \sum_j \exp(V_{ij}^s),}
//' where the inner sum runs over individual i's alternatives and includes the
//' outside option's \eqn{\exp(0)} term when `include_outside_option = TRUE`.
//' The log-sum-exp must be averaged *across draws*: applying log-sum-exp to
//' the draw-averaged utilities returned by `mxl_predict` understates the
//' expectation because log-sum-exp is convex (Jensen's inequality).
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether the outside option is present
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns Vector of length N with the simulated expected logsum per choice
//'   situation.
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' ls <- choicer:::mxl_logsum(coef(fit), d$X, d$W, d$alt_idx, d$M, eta,
//'   rc_dist = rep(0L, ncol(d$W)), rc_correlation = FALSE, rc_mean = FALSE)
//' head(ls)
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::vec mxl_logsum(const arma::vec &theta, const arma::mat &X, const arma::mat &W,
                     const arma::uvec &alt_idx, const Rcpp::IntegerVector &M,
                     const arma::cube &eta_draws, const arma::uvec &rc_dist,
                     const bool rc_correlation = true, const bool rc_mean = false,
                     const bool use_asc = true, const bool include_outside_option = false,
                     const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0) {
  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based alt indices and prefix sums
  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_ls = (gen_seed >= 0);
  HaltonGen halton_gen_ls;
  if (use_generate_ls) {
    halton_gen_ls = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Thread-private buffers (declared here; copied per-thread via firstprivate)
  arma::mat eta_i_buf_ls;
  arma::mat eta_i_store_ls;

  // Output accumulator (each individual writes only its own slot)
  arma::vec logsum = arma::zeros(N);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) firstprivate(eta_i_buf_ls, eta_i_store_ls)
#endif
  for (int i = 0; i < N; ++i) {
    const int m_i = M[i];
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const int start_idx = S_prefix[i];
    const int end_idx = start_idx + m_i - 1;
    const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

    arma::mat W_i =
        make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

    // Pre-computed base utility for this individual
    const arma::vec base_util_i = base_util.subvec(start_idx, end_idx);

    // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
    const arma::mat* eta_i_ptr_ls;
    if (use_generate_ls) {
      halton_gen_ls.fill_eta_i(eta_i_buf_ls, i + 1);
      eta_i_ptr_ls = &eta_i_buf_ls;
    } else {
      eta_i_store_ls = eta_draws.slice(i);
      eta_i_ptr_ls = &eta_i_store_ls;
    }
    const arma::mat& eta_i_ls_ref = *eta_i_ptr_ls;
    arma::mat Gamma_final = batch_gamma_draws(L, eta_i_ls_ref, rc_dist);

    // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
    const arma::mat WGamma = W_i * Gamma_final;

    // Accumulate the per-draw log-sum-exp (NOT the logsum of averaged
    // utilities; see the Jensen note in the docs above).
    double logsum_acc = 0.0;

    // CHANGE #5: hoist per-draw temporaries outside the s loop
    arma::vec inside_utils(m_i);
    arma::vec V_s(num_choices);

    for (int s = 0; s < Sdraw; ++s) {
      inside_utils = base_util_i + WGamma.col(s);

      // Build full V_s with the outside option's V = 0 slot when present
      if (include_outside_option) {
        V_s(0) = 0.0; // outside option fixed at 0
        V_s.subvec(1, num_choices - 1) = inside_utils;
      } else {
        V_s = inside_utils;
      }

      // Stable log-sum-exp (max-subtraction)
      const double V_max = V_s.max();
      logsum_acc += V_max + std::log(arma::accu(arma::exp(V_s - V_max)));
    }

    // Average over draws; disjoint write by individual — no race
    logsum(i) = logsum_acc / static_cast<double>(Sdraw);
  }

  return logsum;
}

//' Predicted aggregate market shares for Mixed Logit
//'
//' Exported wrapper around the internal `mxl_predict_shares_internal`. Parses
//' `theta` using the standard parameter ordering and returns the simulated
//' weighted-average market shares.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns Vector of length J (or J+1 with outside option) of predicted shares.
//' @keywords internal
// [[Rcpp::export]]
arma::vec mxl_predict_shares(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Basic dimensions
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based indexing and prefix sums
  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Number of alternatives for the output
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int num_alts = compute_J_total(J_inside, include_outside_option);

  return mxl_predict_shares_internal(
    X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
    delta, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
    gen_seed, gen_scramble, gen_S
  );
}

//' Diversion ratios for Mixed Logit (simulated, derivative-based)
//'
//' Computes the matrix of attribute-based diversion ratios for a fitted
//' Mixed Logit model. DR(k, j) is the fraction of demand lost by alternative
//' `j` that is captured by alternative `k` when a marginal change in
//' alternative j's `elast_var` attribute reduces s_j.
//'
//' In MNL the per-draw realized coefficient is a constant, so it cancels in
//' the ratio and the result is independent of the variable chosen. In MXL,
//' the realized coefficient \eqn{\beta_{ik}^s} varies across individuals
//' and draws, so the diversion ratio depends on which attribute is perturbed.
//' For a variable with a fixed coefficient the dependence again vanishes
//' (the constant cancels); for a random-coefficient variable it does not.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param elast_var_idx 1-based index of the perturbed variable
//' @param is_random_coef TRUE if the variable is in W (random coef), FALSE if in X (fixed)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns J x J (or (J+1) x (J+1)) matrix of diversion ratios with zero diagonal.
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int elast_var_idx,
    const bool is_random_coef,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Attribute-based diversion ratio (simulated):
  //   DR(k, j) = E_i[ w_i * (1/S) sum_s beta_{ik}^s * P_ij(s) * P_ik(s) ]
  //           /  E_i[ w_i * (1/S) sum_s beta_{ik}^s * P_ij(s) * (1 - P_ij(s)) ]
  // where beta_{ik}^s is the realized coefficient on the perturbed variable
  // for individual i and draw s. For a fixed coefficient beta_{ik}^s is a
  // constant scalar and cancels (MNL property). For a random coefficient it
  // varies and does not cancel. Cross-products P_ij(s)*P_ik(s) must be
  // accumulated INSIDE the draw loop; averaging across draws first and
  // multiplying later is biased. See docs/mixed_logit_math.md.

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

  // Validate the perturbed variable index. Catch the empty-block cases
  // first (K_x=0 with is_random_coef=FALSE, or K_w=0 with is_random_coef=TRUE)
  // with an actionable message before the index-range check.
  if (!is_random_coef && K_x == 0) {
    Rcpp::stop("Cannot compute diversion ratios w.r.t. a fixed-coefficient "
               "variable: the model has no fixed coefficients (K_x = 0). "
               "Did you mean is_random_coef = TRUE?");
  }
  if (is_random_coef && K_w == 0) {
    Rcpp::stop("Cannot compute diversion ratios w.r.t. a random-coefficient "
               "variable: the model has no random coefficients (K_w = 0).");
  }
  const int var_idx = elast_var_idx - 1;
  if (is_random_coef) {
    if (var_idx < 0 || var_idx >= K_w) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for W matrix (K_w=%d).",
                 elast_var_idx, K_w);
    }
  } else {
    if (var_idx < 0 || var_idx >= K_x) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for X matrix (K_x=%d).",
                 elast_var_idx, K_x);
    }
  }

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const double beta_k = is_random_coef ? 0.0 : beta(var_idx);
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based alt indices
  arma::uvec alt_idx0 = alt_idx - 1;

  // Total alternatives for output matrix
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Prefix sums
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_d = (gen_seed >= 0);
  HaltonGen halton_gen_d;
  if (use_generate_d) {
    halton_gen_d = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Global accumulators
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_numerator = arma::zeros(J_total, J_total);
    arma::vec local_denominator = arma::zeros(J_total);
    arma::mat eta_i_buf_d;
    arma::mat eta_i_store_d;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util.subvec(start_idx, end_idx);

      // Map local indices to global alternative indices
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // Per-individual accumulators (sum across draws, divided by S below)
      arma::mat ind_num = arma::zeros(num_choices, num_choices);
      arma::vec ind_den = arma::zeros(num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_d;
      if (use_generate_d) {
        halton_gen_d.fill_eta_i(eta_i_buf_d, i + 1);
        eta_i_ptr_d = &eta_i_buf_d;
      } else {
        eta_i_store_d = eta_draws.slice(i);
        eta_i_ptr_d = &eta_i_store_d;
      }
      const arma::mat& eta_i_d_ref = *eta_i_ptr_d;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_d_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      // Loop over draws — cross-products MUST be accumulated INSIDE this loop
      for (int s = 0; s < Sdraw; ++s) {
        const auto gamma_i_s_final = Gamma_final.col(s); // still needed for random coef value

        // Realized coefficient on the perturbed variable for this (i, s).
        // For a fixed coef this is the constant beta_k; for a random coef
        // it is mu_final(var_idx) + gamma_i_s_final(var_idx), already
        // transformed (exp(.) applied upstream when rc_dist == 1).
        const double beta_k_eff = is_random_coef
            ? (mu_final(var_idx) + gamma_i_s_final(var_idx))
            : beta_k;

        // CHANGE #2: use pre-computed WGamma column instead of W_i * gamma_i_s_final
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Stable softmax
        stable_softmax(V_s, P_s);

        // Accumulate cross-products weighted by beta_k_eff inside the draw loop
        for (int j_local = 0; j_local < num_choices; ++j_local) {
          const double P_j = P_s(j_local);
          ind_den(j_local) += beta_k_eff * P_j * (1.0 - P_j);
          for (int k_local = 0; k_local < num_choices; ++k_local) {
            if (k_local == j_local) continue;
            ind_num(k_local, j_local) += beta_k_eff * P_j * P_s(k_local);
          }
        }
      }  // end s loop

      // Average over draws
      const double S_d = static_cast<double>(Sdraw);
      ind_num /= S_d;
      ind_den /= S_d;

      // Scatter individual contribution into thread-local globals
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map(j_local);
        local_denominator(global_j) += w_i * ind_den(j_local);
        for (int k_local = 0; k_local < num_choices; ++k_local) {
          if (k_local == j_local) continue;
          const int global_k = global_j_map(k_local);
          local_numerator(global_k, global_j) += w_i * ind_num(k_local, j_local);
        }
      }
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_numerator += local_numerator;
      global_denominator += local_denominator;
    }
  }  // end parallel region

  // Final ratios with numerical guard (denominator can be negative when
  // beta_k_eff is negative, e.g. price; check magnitude, not sign)
  arma::mat DR = arma::zeros(J_total, J_total);
  for (int j = 0; j < J_total; ++j) {
    if (std::abs(global_denominator(j)) > 1e-15) {
      for (int k = 0; k < J_total; ++k) {
        if (k != j) {
          DR(k, j) = global_numerator(k, j) / global_denominator(j);
        }
      }
    }
  }

  return DR;
}

//' BLP contraction mapping for mixed logit
//'
//' Finds the ASC (delta) parameters such that predicted market shares
//' match target shares, using the contraction mapping of Berry, Levinsohn,
//' and Pakes (1995).
//'
//' @param delta J-1 or J vector with initial guess for deltas (ASCs)
//' @param target_shares J vector with target market shares
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param beta K_x vector with fixed coefficients
//' @param mu K_w vector with mean parameters (raw, will be transformed if log-normal)
//' @param L_params Cholesky parameters vector
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters represent means (TRUE) or are zero (FALSE)
//' @param include_outside_option whether outside option is included
//' @param tol convergence tolerance (default 1e-8)
//' @param max_iter maximum iterations (default 1000)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations, \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns vector with converged delta (ASC) values
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' pm <- fit$param_map
//' delta <- mxl_blp_contraction(rep(0, J), rep(1/J, J), d$X, d$W,
//'   coef(fit)[pm$beta], rep(0, ncol(d$W)), coef(fit)[pm$sigma],
//'   d$alt_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' delta
//' }
//' @export
// [[Rcpp::export]]
arma::vec mxl_blp_contraction(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu,
    const arma::vec& L_params,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool include_outside_option = false,
    const double tol = 1e-8,
    const int max_iter = 1000,
    const int gen_seed = -1,
    const int gen_scramble = 1,
    const int gen_S = 0
) {
  const int K_w = W.n_cols;
  const bool use_asc = true;

  // delta is harmonized to cover every referenced alternative below, so the
  // ASC-coverage check is skipped here (use_asc = false, empty delta).
  check_rc_dist_length(rc_dist, K_w);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws,
                        /*use_asc=*/false, arma::vec(), &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, /*use_asc=*/false, arma::vec(),
                         &weights);
    check_rc_dist_length(rc_dist, K_w);
  }

  // Build L matrix
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  // Transform mu for log-normal coefficients
  arma::vec mu_final = mu;
  if (rc_mean) {
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) {  // log-normal
        mu_final(k) = std::exp(mu(k));
      }
    }
  }

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Deduce number of alternatives from data
  const int J_inside = static_cast<int>(arma::max(alt_idx0)) + 1;  // inside options
  const int num_alts = include_outside_option ? (J_inside + 1) : J_inside;  // total options incl. outside

  // Validate target shares length
  if (static_cast<int>(target_shares.n_elem) != num_alts) {
    Rcpp::stop("Error: target_shares must have length %d (total alternatives, incl. outside if present).", num_alts);
  }

  // Compute prefix sums
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // ---------------------------------------------------------------------------
  // Harmonize delta input:
  // - include_outside_option = TRUE: delta must have length J_inside (inside
  //   alternatives only). The outside option has no ASC (utility normalized to 0).
  // - include_outside_option = FALSE: allow J-1 (free) or J (full) length. Pad
  //   a leading zero if only free deltas are provided; keep baseline anchored.
  // ---------------------------------------------------------------------------
  arma::vec delta_current;  // length J_inside
  if (include_outside_option) {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else {
      Rcpp::stop("Error: delta must have length %d (inside alternatives only).", J_inside);
    }
  } else {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else if (static_cast<int>(delta.n_elem) == J_inside - 1) {
      delta_current = arma::zeros(J_inside);
      delta_current.subvec(1, J_inside - 1) = delta;  // pad baseline with 0
    } else {
      Rcpp::stop("Error: delta must have length %d (full) or %d (free, with baseline omitted).",
                 J_inside, J_inside - 1);
    }
    // Anchor baseline at zero for identification
    delta_current -= delta_current(0);
  }

  // Compute initial predicted shares
  arma::vec shares_pred = mxl_predict_shares_internal(
    X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
    delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
    gen_seed, gen_scramble, gen_S
  );

  // Work with inside shares only for the contraction step
  arma::vec shares_pred_inside = include_outside_option
                                 ? shares_pred.subvec(1, num_alts - 1)
                                 : shares_pred;
  arma::vec target_shares_inside = include_outside_option
                                   ? target_shares.subvec(1, num_alts - 1)
                                   : target_shares;

  // Guard against zeros before taking logs
  auto validate_positive = [](const arma::vec& v, const char* name) {
    if (v.min() <= 0.0) {
      Rcpp::stop("%s contains non-positive entries; cannot take logarithm.", name);
    }
  };
  validate_positive(shares_pred_inside, "Predicted shares");
  validate_positive(target_shares_inside, "Target shares");

  arma::vec log_shares_old = arma::log(shares_pred_inside);
  arma::vec log_shares_target = arma::log(target_shares_inside);

  // Iteration
  int iter = 0;
  double residual = 10.0;

  while (iter < max_iter) {
    Rcpp::checkUserInterrupt(); // H4: allow user to interrupt long-running contraction
    arma::vec delta_new = delta_current + (log_shares_target - log_shares_old);

    // Re-anchor baseline each iteration when there is no outside option
    if (!include_outside_option) {
      delta_new -= delta_new(0);
    }

    residual = arma::max(arma::abs(delta_new - delta_current));

    if (residual < tol) {
      break;
    }

    delta_current = delta_new;
    shares_pred = mxl_predict_shares_internal(
      X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
      delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
      gen_seed, gen_scramble, gen_S
    );

    shares_pred_inside = include_outside_option
                         ? shares_pred.subvec(1, num_alts - 1)
                         : shares_pred;
    validate_positive(shares_pred_inside, "Predicted shares");
    log_shares_old = arma::log(shares_pred_inside);
    ++iter;
  }

  if (iter >= max_iter) {
    Rcpp::Rcout << "Warning: Maximum iterations reached without convergence. Residual: "
                << residual << std::endl;
  }

  return delta_current;
}

//' Compute aggregate elasticities for mixed logit model
//'
//' Computes the aggregate elasticity matrix (weighted average of individual
//' elasticities) for the Mixed Logit model. The elasticity E(i,j) represents
//' the percentage change in the probability of choosing alternative i when
//' the attribute of alternative j changes by 1%.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param choice_idx N x 1 vector (kept for API consistency, not used)
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param elast_var_idx 1-based index of the variable for elasticity computation
//' @param is_random_coef TRUE if variable is in W (random coef), FALSE if in X (fixed coef)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns J x J matrix of aggregate elasticities
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' elas <- choicer:::mxl_elasticities_parallel(coef(fit), d$X, d$W, d$alt_idx,
//'   d$choice_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   elast_var_idx = 1L, is_random_coef = FALSE,
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' elas
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int elast_var_idx,
    const bool is_random_coef,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  (void)choice_idx;  // unused, kept for API consistency

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

  // Convert 1-based R index to 0-based C++ index
  const int var_idx = elast_var_idx - 1;

  // Validate variable index
  if (is_random_coef) {
    if (var_idx < 0 || var_idx >= K_w) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for W matrix (K_w=%d).",
                 elast_var_idx, K_w);
    }
  } else {
    if (var_idx < 0 || var_idx >= K_x) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for X matrix (K_x=%d).",
                 elast_var_idx, K_x);
    }
  }

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const double beta_k = is_random_coef ? 0.0 : beta(var_idx);
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Determine total number of alternatives
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Compute prefix sums
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_e = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_e = (gen_seed >= 0);
  HaltonGen halton_gen_e;
  if (use_generate_e) {
    halton_gen_e = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Global accumulators
  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;
    arma::mat eta_i_buf_e;
    arma::mat eta_i_store_e;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const auto X_i = X.rows(start_idx, end_idx);
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_e.subvec(start_idx, end_idx);

      // Map local indices to global alternative indices
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // Get attribute values for the elasticity variable
      arma::vec x_k_i = arma::zeros(num_choices);
      if (is_random_coef) {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = W_i.col(var_idx);
        } else {
          x_k_i = W_i.col(var_idx);
        }
      } else {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = X_i.col(var_idx);
        } else {
          x_k_i = X_i.col(var_idx);
        }
      }

      // Compute P_bar (average probabilities) and accumulate elasticity terms
      arma::vec P_bar_i = arma::zeros(num_choices);
      arma::mat elas_accum = arma::zeros(num_choices, num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_e;
      if (use_generate_e) {
        halton_gen_e.fill_eta_i(eta_i_buf_e, i + 1);
        eta_i_ptr_e = &eta_i_buf_e;
      } else {
        eta_i_store_e = eta_draws.slice(i);
        eta_i_ptr_e = &eta_i_store_e;
      }
      const arma::mat& eta_i_e_ref = *eta_i_ptr_e;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_e_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      for (int s = 0; s < Sdraw; ++s) {
        const auto gamma_i_s_final = Gamma_final.col(s); // still needed for random coef value

        // Get effective coefficient for this draw
        double beta_k_eff;
        if (is_random_coef) {
          beta_k_eff = mu_final(var_idx) + gamma_i_s_final(var_idx);
        } else {
          beta_k_eff = beta_k;
        }

        // CHANGE #2: use pre-computed WGamma column instead of W_i * gamma_i_s_final
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Compute probabilities
        stable_softmax(V_s, P_s);

        P_bar_i += P_s;

        // Accumulate elasticity terms for this draw
        for (int j_local = 0; j_local < num_choices; ++j_local) {
          const double P_j = P_s(j_local);

          for (int m_local = 0; m_local < num_choices; ++m_local) {
            const double P_m = P_s(m_local);
            const double x_km = x_k_i(m_local);

            double elas_term;
            if (j_local == m_local) {
              // Own-elasticity: beta_k * x_k * P_j * (1 - P_j)
              elas_term = beta_k_eff * x_km * P_j * (1.0 - P_j);
            } else {
              // Cross-elasticity: -beta_k * x_km * P_j * P_m
              elas_term = -beta_k_eff * x_km * P_j * P_m;
            }

            elas_accum(j_local, m_local) += elas_term;
          }
        }
      }  // end s loop

      P_bar_i /= static_cast<double>(Sdraw);
      elas_accum /= static_cast<double>(Sdraw);

      // Compute final elasticities: E = elas_accum / P_bar
      // and map to global indices
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map(j_local);
        const double P_bar_j = P_bar_i(j_local);

        if (P_bar_j > 1e-12) {  // Avoid division by zero
          for (int m_local = 0; m_local < num_choices; ++m_local) {
            const int global_m = global_j_map(m_local);
            double elasticity = elas_accum(j_local, m_local) / P_bar_j;
            local_elas_matrix(global_j, global_m) += w_i * elasticity;
          }
        }
      }

      local_total_weight += w_i;
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_elas_matrix += local_elas_matrix;
      global_total_weight += local_total_weight;
    }
  }  // end parallel region

  // Compute weighted average
  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}
