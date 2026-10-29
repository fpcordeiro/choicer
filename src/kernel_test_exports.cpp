// kernel_test_exports.cpp — Thin Rcpp wrappers exposing choicer_internal.h
// and choicer_derivs.h helpers for unit testing. These functions are NOT
// user-facing API; they are @noRd and only exported so that tests can pin
// the raw-array softmax, log-sum-exp and nested logit probabilities to their
// Armadillo counterparts bit for bit, check how Armadillo is configured (its
// index word), check how the derivative kernels choose and size their
// accumulators, read the physical memory the estimation kernels' memory
// check compares with, and run the fixed-point solver on a test map.
//
// DO NOT add any of these to the public documentation or NAMESPACE.

// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "choicer_derivs.h"
#include "blp_engine.h"
#include <algorithm>

//' stable_softmax_n(), log_sum_exp_n() and max_shifted_lse_n() next to
//' stable_softmax(), logSumExp() and the Armadillo expression of the mixed
//' logit log-sum, for a bitwise comparison in tests
//'
//' @param v Numeric vector of utilities (length >= 1).
//' @return List with the shifted utilities, probabilities and log-denominator
//'   from both softmax versions, and both versions of each log-sum-exp of v.
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_softmax_n(const arma::vec& v) {
  arma::vec v_arma = v, p_arma;
  const double ld_arma = stable_softmax(v_arma, p_arma);
  arma::vec v_raw = v, p_raw(v.n_elem);
  const int n = static_cast<int>(v.n_elem);
  const double ld_raw = stable_softmax_n(v_raw.memptr(), p_raw.memptr(), n);
  return Rcpp::List::create(
      Rcpp::Named("v_arma") = v_arma, Rcpp::Named("v_raw") = v_raw,
      Rcpp::Named("p_arma") = p_arma, Rcpp::Named("p_raw") = p_raw,
      Rcpp::Named("log_denom_arma") = ld_arma,
      Rcpp::Named("log_denom_raw") = ld_raw,
      Rcpp::Named("lse_arma") = logSumExp(v),
      Rcpp::Named("lse_raw") = log_sum_exp_n(v.memptr(), n),
      Rcpp::Named("lse_shift_arma") =
          v.max() + std::log(arma::accu(arma::exp(v - v.max()))),
      Rcpp::Named("lse_shift_raw") = max_shifted_lse_n(v.memptr(), n));
}

// nl_individual_probs() as it was written on Armadillo expressions (src/choicer.h
// before the raw-array rewrite), kept verbatim as the reference the rewrite
// must reproduce bit for bit.
static void nl_individual_probs_ref(
    const arma::vec& V_inside,
    const arma::uvec& nest_idx0_i,
    const arma::vec& lambda,
    const int n_nests,
    const bool include_outside_option,
    arma::vec& P_i,
    arma::vec& P_j_given_k,
    arma::vec& P_k,
    arma::vec& log_I_k,
    arma::vec& log_P_i,
    double& log_P_outside
) {
  const int m_i = V_inside.n_elem;

  // V_ij / lambda_k  (lambda_k = 1 for singletons)
  arma::vec V_over_lambda = V_inside / lambda.elem(nest_idx0_i);

  // --- log_I_k (inclusive value) via log-sum-exp within each nest ---
  arma::vec max_V_k = arma::vec(n_nests).fill(-arma::datum::inf);
  for (int j = 0; j < m_i; ++j) {
    const int k = nest_idx0_i[j];
    if (V_over_lambda[j] > max_V_k[k]) {
      max_V_k[k] = V_over_lambda[j];
    }
  }

  arma::vec I_k_unscaled = arma::zeros(n_nests);
  for (int j = 0; j < m_i; ++j) {
    const int k = nest_idx0_i[j];
    if (std::isfinite(max_V_k[k])) {
      I_k_unscaled[k] += std::exp(V_over_lambda[j] - max_V_k[k]);
    }
  }

  log_I_k = arma::vec(n_nests).fill(-arma::datum::inf);
  for (int k = 0; k < n_nests; ++k) {
    if (I_k_unscaled[k] > 0) {
      log_I_k[k] = max_V_k[k] + std::log(I_k_unscaled[k]);
    }
  }

  // --- log(P_k) (nest probability) ---
  arma::vec nest_terms = lambda % log_I_k;

  double max_nest_term = nest_terms.max();
  if (!std::isfinite(max_nest_term)) {
    max_nest_term = 0;
  }

  double sum_exp_nest_terms = arma::accu(arma::exp(nest_terms - max_nest_term));
  if (include_outside_option) {
    // Outside option: V=0, lambda=1 -> term = 0
    sum_exp_nest_terms += std::exp(0.0 - max_nest_term);
  }

  const double log_denom_P_nest = max_nest_term + std::log(sum_exp_nest_terms);

  arma::vec log_P_k = nest_terms - log_denom_P_nest;
  log_P_outside = include_outside_option ? (0.0 - log_denom_P_nest)
                                         : -arma::datum::inf;

  // --- log(P_j|k) and log(P_ij) ---
  arma::vec log_P_j_given_k = V_over_lambda - log_I_k.elem(nest_idx0_i);
  log_P_i = log_P_j_given_k + log_P_k.elem(nest_idx0_i);

  P_i = arma::exp(log_P_i);
  P_j_given_k = arma::exp(log_P_j_given_k);
  P_k = arma::exp(log_P_k);
}

//' nl_individual_probs() next to its Armadillo-expression reference, for a
//' bitwise comparison in tests
//'
//' @param V Inside utilities of one individual.
//' @param nest0 0-based nest of each inside alternative (same length as V).
//' @param lambda Full lambda vector, one entry per nest.
//' @param include_outside_option Whether the outside option is in the set.
//' @return List with each output of both versions (`*_ref` the reference).
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_nl_individual_probs(const arma::vec& V,
                                    const Rcpp::IntegerVector& nest0,
                                    const arma::vec& lambda,
                                    const bool include_outside_option) {
  const int m = V.n_elem;
  const int n_nests = lambda.n_elem;
  if (nest0.size() != m) Rcpp::stop("nest0 must have one entry per utility.");
  for (int j = 0; j < m; ++j) {
    if (nest0[j] < 0 || nest0[j] >= n_nests) Rcpp::stop("nest0 out of range.");
  }
  arma::uvec nest_u(m);
  for (int j = 0; j < m; ++j) nest_u[j] = static_cast<arma::uword>(nest0[j]);
  arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
  double log_P_outside;
  nl_individual_probs_ref(V, nest_u, lambda, n_nests, include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i,
                          log_P_outside);
  // Buffers larger than this individual needs and filled with R's NA (a NaN
  // with a payload no arithmetic produces), as a thread's are after earlier,
  // larger individuals: an element read or left before this call writes it
  // would show.
  NlProbs pr(m + 3, n_nests);
  for (std::vector<double>* b :
       {&pr.P_i, &pr.P_j_given_k, &pr.log_P_i, &pr.P_k, &pr.log_I_k,
        &pr.V_over_lambda, &pr.log_P_j_given_k, &pr.max_V_k,
        &pr.I_k_unscaled, &pr.nest_terms, &pr.log_P_k}) {
    std::fill(b->begin(), b->end(), NA_REAL);
  }
  pr.log_P_outside = NA_REAL;
  nl_individual_probs(V.memptr(), nest0.begin(), m, lambda, n_nests,
                      include_outside_option, pr);
  auto num = [](const std::vector<double>& x, const int n) {
    return Rcpp::NumericVector(x.begin(), x.begin() + n);
  };
  return Rcpp::List::create(
      Rcpp::Named("P_i_ref") = Rcpp::NumericVector(P_i.begin(), P_i.end()),
      Rcpp::Named("P_i") = num(pr.P_i, m),
      Rcpp::Named("P_j_given_k_ref") =
          Rcpp::NumericVector(P_j_given_k.begin(), P_j_given_k.end()),
      Rcpp::Named("P_j_given_k") = num(pr.P_j_given_k, m),
      Rcpp::Named("P_k_ref") = Rcpp::NumericVector(P_k.begin(), P_k.end()),
      Rcpp::Named("P_k") = num(pr.P_k, n_nests),
      Rcpp::Named("log_I_k_ref") =
          Rcpp::NumericVector(log_I_k.begin(), log_I_k.end()),
      Rcpp::Named("log_I_k") = num(pr.log_I_k, n_nests),
      Rcpp::Named("log_P_i_ref") =
          Rcpp::NumericVector(log_P_i.begin(), log_P_i.end()),
      Rcpp::Named("log_P_i") = num(pr.log_P_i, m),
      Rcpp::Named("log_P_outside_ref") = log_P_outside,
      Rcpp::Named("log_P_outside") = pr.log_P_outside);
}

//' nl_individual_probs_offered() over the nests of one individual's rows, for
//' a bitwise comparison with nl_individual_probs() in tests
//'
//' @inheritParams test_nl_individual_probs
//' @return List with the outputs, the per-nest ones over every nest: R's NA
//'   (the buffers' fill) where a nest is not offered, which the routine must
//'   not write; and `nests`, the offered nests (0-based, increasing).
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_nl_individual_probs_offered(const arma::vec& V,
                                            const Rcpp::IntegerVector& nest0,
                                            const arma::vec& lambda,
                                            const bool include_outside_option) {
  const int m = V.n_elem;
  const int n_nests = lambda.n_elem;
  if (nest0.size() != m) Rcpp::stop("nest0 must have one entry per utility.");
  for (int j = 0; j < m; ++j) {
    if (nest0[j] < 0 || nest0[j] >= n_nests) Rcpp::stop("nest0 out of range.");
  }
  if (!nl_lambda_positive(lambda)) Rcpp::stop("every lambda must be positive.");
  std::vector<int> nests(nest0.begin(), nest0.end());
  std::sort(nests.begin(), nests.end());
  nests.erase(std::unique(nests.begin(), nests.end()), nests.end());
  const int n_o = static_cast<int>(nests.size());
  // Buffers filled with R's NA, as in test_nl_individual_probs()
  NlProbs pr(m + 3, n_nests);
  for (std::vector<double>* b :
       {&pr.P_i, &pr.P_j_given_k, &pr.log_P_i, &pr.P_k, &pr.log_I_k,
        &pr.V_over_lambda, &pr.log_P_j_given_k, &pr.max_V_k,
        &pr.I_k_unscaled, &pr.nest_terms, &pr.log_P_k}) {
    std::fill(b->begin(), b->end(), NA_REAL);
  }
  pr.log_P_outside = NA_REAL;
  nl_individual_probs_offered(V.memptr(), nest0.begin(), m, lambda, nests.data(),
                              n_o, include_outside_option, pr);
  auto num = [](const std::vector<double>& x, const int n) {
    return Rcpp::NumericVector(x.begin(), x.begin() + n);
  };
  return Rcpp::List::create(
      Rcpp::Named("P_i") = num(pr.P_i, m),
      Rcpp::Named("P_j_given_k") = num(pr.P_j_given_k, m),
      Rcpp::Named("P_k") = num(pr.P_k, n_nests),
      Rcpp::Named("log_I_k") = num(pr.log_I_k, n_nests),
      Rcpp::Named("log_P_i") = num(pr.log_P_i, m),
      Rcpp::Named("log_P_outside") = pr.log_P_outside,
      Rcpp::Named("nests") = Rcpp::IntegerVector(nests.begin(), nests.end()));
}

//' Bytes in Armadillo's index type, arma::uword: 8 when choicer is built with
//' ARMA_64BIT_WORD (src/Makevars) on a 64-bit platform
//'
//' @return An integer.
//' @noRd
// [[Rcpp::export(rng = false)]]
int test_arma_word_bytes() {
  return static_cast<int>(sizeof(arma::uword));
}

//' Element count of a no-copy Armadillo view of a (2^31 - 1) x 3 matrix, as
//' RcppArmadillo views R's matrices: 6442450941 under a 64-bit word, wrapped
//' to 2147483645 under a 32-bit one. The view is never read.
//'
//' @return A double.
//' @noRd
// [[Rcpp::export(rng = false)]]
double test_arma_view_n_elem() {
  double cell = 0.0;
  const arma::mat view(&cell, 2147483647u, 3u, false, true);
  return static_cast<double>(view.n_elem);
}

//' The derivative kernels' choice of accumulator, as SymAcc makes it:
//' the rows s = sym_acc_rows() a thread keeps per column (n, or K_c when the
//' threads share the result) and its buffer's column offsets
//' sym_acc_col_off(j, s), j = 0..n
//'
//' @param n,K_c Numbers of parameters and of continuous parameters.
//' @param T Threads.
//' @param acc_bytes The kernels' budget for the private triangles.
//' @return List with \code{shared} (s < n) and the n + 1 offsets (doubles),
//'   the last one the buffer's length.
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List test_sym_acc(const double n, const int K_c, const int T,
                        const double acc_bytes) {
  const choicer_off nn = static_cast<choicer_off>(n);
  const choicer_off s = sym_acc_rows(nn, K_c, T, acc_bytes);
  const bool shared = s < nn;
  Rcpp::NumericVector off(static_cast<R_xlen_t>(nn) + 1);
  for (choicer_off j = 0; j <= nn; ++j) {
    off[static_cast<R_xlen_t>(j)] = static_cast<double>(sym_acc_col_off(j, s));
  }
  return Rcpp::List::create(Rcpp::Named("shared") = shared,
                            Rcpp::Named("col_off") = off);
}

//' The machine's physical memory in bytes as the mixed logit kernels' memory
//' check reads it (sysmem.c), 0 when it cannot be read
//'
//' @return A double.
//' @noRd
// [[Rcpp::export(rng = false)]]
double test_physical_memory() {
  return choicer_physical_memory();
}

// A map for the fixed-point solver's tests (fixed_point.h): the step
// f(x) = b - A x, its residual max |f(x)| or, for the first evaluations, a
// scripted one (NA: the computed one), evaluations that fail on a schedule
// (the k-th, 1-based, fails where `fail` lists k), and a fallback step f(x) / 2
// while the solver asks for it (where it has one). Every evaluated point is
// recorded, with its residual and whether it was the fallback's.
struct FpTestMap {
  const arma::mat& A;
  const arma::vec& b;
  const std::vector<int>& fail;
  const std::vector<double>& script;
  const bool fallback;
  bool fb = false;
  int count = 0;
  std::vector<double> trace_x, trace_res;
  std::vector<int> trace_fb;

  FpTestMap(const arma::mat& A_, const arma::vec& b_, const std::vector<int>& fail_,
            const std::vector<double>& script_, const bool fallback_)
      : A(A_), b(b_), fail(fail_), script(script_), fallback(fallback_) {}

  arma::uword n() const { return b.n_elem; }
  bool has_fallback() const { return fallback; }
  void use_fallback(const bool on) { fb = on; }
  bool eval(const arma::vec& x, double& res, arma::vec& step) {
    ++count;
    trace_x.insert(trace_x.end(), x.begin(), x.end());
    trace_fb.push_back(fb ? 1 : 0);
    if (std::find(fail.begin(), fail.end(), count) != fail.end()) {
      res = arma::datum::inf;
      trace_res.push_back(res);
      return false;
    }
    step = b - A * x;
    res = arma::abs(step).max();
    const std::size_t k = static_cast<std::size_t>(count - 1);
    if (k < script.size() && !std::isnan(script[k])) res = script[k];
    if (fb) step *= 0.5;
    trace_res.push_back(res);
    return true;
  }
};

//' The fixed-point solver on a test map, for tests
//'
//' @param A,b The linear step f(x) = b - A x (A square, b its length).
//' @param x0 The start.
//' @param control List from .blp_control() (its method, budget, tol,
//'   damping, patience, window and the methods' settings).
//' @param fail 1-based indices of the evaluations that fail.
//' @param script Residuals of the first evaluations (NA: computed).
//' @param fallback Whether the map has a fallback step, f(x) / 2.
//' @return List: x (the point returned), the solver's record (code,
//'   iterations, evaluations, backtracks, stalls, race, resumes, restarts,
//'   fallback, residual) and the evaluations' trace (trace_x, an n x
//'   evaluations matrix; trace_res, the residuals; trace_fb, whether each was
//'   the fallback's).
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_fp_solve(const arma::mat& A, const arma::vec& b,
                         const arma::vec& x0, const Rcpp::List& control,
                         const Rcpp::IntegerVector& fail,
                         const Rcpp::NumericVector& script, const bool fallback) {
  if (A.n_rows != b.n_elem || A.n_cols != b.n_elem || x0.n_elem != b.n_elem) {
    Rcpp::stop("A must be n x n and b and x0 of length n.");
  }
  const FpControl fc = blp_fp_control(blp_control(control));
  const std::vector<int> fail_v(fail.begin(), fail.end());
  const std::vector<double> script_v(script.begin(), script.end());
  FpTestMap map(A, b, fail_v, script_v, fallback);
  arma::vec x = x0;
  const FpResult r = fp_solve(map, x, fc);
  arma::mat tx(map.trace_x.data(), b.n_elem, static_cast<arma::uword>(map.count));
  return Rcpp::List::create(
      Rcpp::Named("x") = Rcpp::NumericVector(x.begin(), x.end()),
      Rcpp::Named("code") = r.code, Rcpp::Named("iterations") = r.iterations,
      Rcpp::Named("evaluations") = r.evaluations,
      Rcpp::Named("backtracks") = r.backtracks, Rcpp::Named("stalls") = r.stalls,
      Rcpp::Named("race") = r.race, Rcpp::Named("resumes") = r.resumes,
      Rcpp::Named("restarts") = r.restarts, Rcpp::Named("fallback") = r.fallback,
      Rcpp::Named("residual") = r.residual, Rcpp::Named("trace_x") = tx,
      Rcpp::Named("trace_res") = Rcpp::wrap(map.trace_res),
      Rcpp::Named("trace_fb") = Rcpp::wrap(map.trace_fb));
}
