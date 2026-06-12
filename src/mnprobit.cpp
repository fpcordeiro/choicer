// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "bayes_samplers.h"
#include <R_ext/Utils.h>   // R_CheckUserInterrupt

// ============================================================================
// Bayesian multinomial probit (MNP) via Gibbs sampling with data augmentation
// (Albert & Chib 1993; McCulloch & Rossi 1994).
//
// Model, in utility differences against a base alternative (p = J - 1):
//   w_i = X_i beta + eps_i,   eps_i ~ N_p(0, Sigma)
//   y_i = j  (j in 1..p)  iff  w_ij > max(0, max_{k != j} w_ik)
//   y_i = 0  (base)       iff  max_j w_ij < 0
//
// Priors: beta ~ N(beta_bar, A^{-1}), Sigma ~ IW(nu, V). The chain runs on
// the non-identified parameterization (unrestricted Sigma); identified
// quantities beta / sqrt(sigma_11) and Sigma / sigma_11 are computed per
// draw in the R wrapper (run_mnprobit).
//
// Parallelism: given (beta, Sigma) the latent w_i are conditionally
// independent across choice situations, so the augmentation sweep is an
// exact Gibbs step when parallelized over i. Each (iteration, observation)
// task draws from its own RNG stream (rng.h), so results are bitwise
// reproducible regardless of the number of OpenMP threads.
//
// The whole chain runs inside ONE OpenMP parallel region: the latent sweep
// and the Mu refresh are work-shared loops, while the beta / Sigma
// conditionals run on the master thread between barriers. Per-iteration
// region setup/teardown is avoided, and no BLAS is called inside the region
// (the master-path linear algebra is hand-rolled in a fixed summation
// order), so worker threads only ever wait at lightweight barriers and the
// draws stay independent of the thread count. Errors and user interrupts
// inside the region are reported through a shared abort flag and raised
// after the region ends (exceptions must not cross an OpenMP boundary;
// R_ToplevelExec lets the master poll for interrupts without a longjmp).
// ============================================================================

namespace {

void chk_int_fn(void*) { R_CheckUserInterrupt(); }

// True if a user interrupt is pending. Safe inside the parallel region as
// long as only the master (R's main) thread calls it.
bool pending_interrupt() {
  return R_ToplevelExec(chk_int_fn, nullptr) == FALSE;
}

}  // namespace

static void validate_mnp_inputs(const arma::mat& X, const Rcpp::IntegerVector& y,
                                const int p, const arma::vec& beta_bar,
                                const arma::mat& A, const double nu,
                                const arma::mat& V, const int R,
                                const int burn, const int thin,
                                const double seed) {
  const int N = y.size();
  const int K = X.n_cols;
  if (p < 1) {
    Rcpp::stop("p must be >= 1.");
  }
  if (N < 1) {
    Rcpp::stop("y must contain at least one choice situation.");
  }
  if (K < 1) {
    Rcpp::stop("X must have at least one column.");
  }
  if (X.n_rows != static_cast<arma::uword>(N) * p) {
    Rcpp::stop("X has %d rows but N * p is %d.",
               static_cast<int>(X.n_rows), N * p);
  }
  if (!X.is_finite()) {
    Rcpp::stop("X contains non-finite values.");
  }
  if (static_cast<int>(beta_bar.n_elem) != K) {
    Rcpp::stop("beta_bar length (%d) does not match the number of columns of X (%d).",
               static_cast<int>(beta_bar.n_elem), K);
  }
  if (static_cast<int>(A.n_rows) != K || static_cast<int>(A.n_cols) != K) {
    Rcpp::stop("A must be a %d x %d matrix.", K, K);
  }
  if (static_cast<int>(V.n_rows) != p || static_cast<int>(V.n_cols) != p) {
    Rcpp::stop("V must be a %d x %d matrix.", p, p);
  }
  if (nu < p) {
    Rcpp::stop("nu (%g) must be >= p (%d) for a proper inverse-Wishart prior.",
               nu, p);
  }
  for (int i = 0; i < N; ++i) {
    if (y[i] < 0 || y[i] > p) {
      Rcpp::stop("y[%d] = %d is outside {0, ..., p} (0 = base alternative).",
                 i + 1, y[i]);
    }
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
}

//' Gibbs sampler for the Bayesian multinomial probit model
//'
//' Runs the McCulloch-Rossi (1994) Gibbs sampler with Albert-Chib data
//' augmentation for the multinomial probit model in utility differences
//' against a base alternative. The chain operates on the non-identified
//' parameterization (unrestricted \code{Sigma}); identified quantities are
//' obtained by normalizing each draw by \code{sigma_11} (handled by
//' \code{\link{run_mnprobit}}).
//'
//' The latent-utility sweep is parallelized with OpenMP across choice
//' situations (they are conditionally independent given \code{beta} and
//' \code{Sigma}). Each (iteration, observation) pair uses its own RNG
//' stream, so draws are reproducible independent of the number of threads
//' (see \code{set_num_threads()}).
//'
//' @param X (N*p) x K stacked design matrix of utility differences. Rows are
//'   grouped by choice situation, with the p = J - 1 difference rows of
//'   situation i ordered by alternative.
//' @param y N vector of choices: 0 for the base alternative, j in 1..p for
//'   the j-th non-base alternative.
//' @param p Number of utility differences (J - 1).
//' @param beta_bar K vector, prior mean of beta.
//' @param A K x K prior precision matrix of beta.
//' @param nu Inverse-Wishart prior degrees of freedom (>= p).
//' @param V p x p inverse-Wishart prior scale matrix.
//' @param R Total number of Gibbs iterations.
//' @param burn Number of initial iterations to discard (0 <= burn < R).
//' @param thin Keep every thin-th post-burn-in draw.
//' @param seed Master RNG seed (non-negative; all streams derive from it).
//' @param trace Print progress every \code{trace} iterations (0 = silent).
//' @returns List with \code{betadraw} (R_keep x K), \code{sigmadraw}
//'   (R_keep x p(p+1)/2, lower triangle of Sigma in row-major order), and
//'   \code{R_keep}.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 100; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
//' out <- mnp_gibbs(d$X, d$y, d$p,
//'   beta_bar = rep(0, d$K), A = 0.01 * diag(d$K),
//'   nu = d$p + 3, V = (d$p + 3) * diag(d$p),
//'   R = 500, burn = 100, thin = 1, seed = 42)
//' colMeans(out$betadraw)
//' }
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List mnp_gibbs(const arma::mat& X,
                     const Rcpp::IntegerVector& y,
                     const int p,
                     const arma::vec& beta_bar,
                     const arma::mat& A,
                     const double nu,
                     const arma::mat& V,
                     const int R,
                     const int burn,
                     const int thin,
                     const double seed,
                     const int trace = 0) {
  validate_mnp_inputs(X, y, p, beta_bar, A, nu, V, R, burn, thin, seed);

  const int N = y.size();
  const int K = X.n_cols;
  const uint64_t useed = static_cast<uint64_t>(seed);

  // One-time cross-product blocks G[j, k] = X_j' X_k, where X_j collects the
  // component-j difference rows of all situations. These let the posterior
  // precision for beta be assembled in O(p^2 K^2) per iteration, independent
  // of N: sum_i X_i' Omega X_i = sum_{j,k} Omega(j, k) G[j, k].
  std::vector<arma::mat> G(static_cast<size_t>(p) * p);
  {
    std::vector<arma::mat> Xj(p);
    arma::uvec idx(N);
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < N; ++i) {
        idx(i) = static_cast<arma::uword>(i) * p + j;
      }
      Xj[j] = X.rows(idx);
    }
    for (int j = 0; j < p; ++j) {
      for (int k = j; k < p; ++k) {
        G[static_cast<size_t>(j) * p + k] = Xj[j].t() * Xj[k];
        if (k != j) {
          G[static_cast<size_t>(k) * p + j] =
            G[static_cast<size_t>(j) * p + k].t();
        }
      }
    }
  }
  const arma::vec Ab = A * beta_bar;

  // Chain state, allocated once and updated in place for the whole run.
  arma::vec beta = beta_bar;
  arma::mat Sigma = arma::eye(p, p);
  arma::mat Omega = arma::eye(p, p);
  arma::mat W(p, N);                       // latent utility differences
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < p; ++j) {
      W(j, i) = (y[i] == j + 1) ? 1.0 : -1.0;   // feasible start
    }
  }

  // Pre-allocated workspaces (all shared; only the master thread writes the
  // small ones, work-shared loops write disjoint columns of W / Mu).
  arma::mat Mu(p, N), Vm(p, N);
  arma::mat Q(K, K), U(K, K), S(p, p), Vpost(p, p);
  arma::vec bvec(K), btilde(K), z(K);

  // Hoisted per-iteration invariants of the latent sweep: the univariate
  // conditionals of N_p(mu_i, Sigma) depend on i only through mu_i and the
  // residuals, so tau_j = sqrt(1 / Omega(j, j)) and Ratio(k, j) =
  // Omega(k, j) / Omega(j, j) are computed once per Sigma draw instead of
  // once per (i, j).
  arma::mat Ratio(p, p);
  arma::vec tau(p);
  for (int j = 0; j < p; ++j) {
    const double ojj = Omega(j, j);
    tau(j) = std::sqrt(1.0 / ojj);
    for (int k = 0; k < p; ++k) Ratio(k, j) = Omega(k, j) / ojj;
  }

  const double* Xv = X.memptr();
  const arma::uword NP = static_cast<arma::uword>(N) * p;

  // Initial Mu = X beta, same per-observation summation order as the
  // in-loop refresh.
  for (int i = 0; i < N; ++i) {
    double* mucol = Mu.colptr(i);
    const arma::uword row0 = static_cast<arma::uword>(i) * p;
    for (int j = 0; j < p; ++j) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += Xv[k * NP + row0 + j] * beta(k);
      mucol[j] = s;
    }
  }

  const int R_keep = (R - burn + thin - 1) / thin;
  arma::mat betadraw(R_keep, K);
  arma::mat sigmadraw(R_keep, p * (p + 1) / 2);

  // Abort protocol: the master thread records a reason and iteration, every
  // thread re-reads the flag after a barrier and leaves the loop together;
  // the corresponding Rcpp::stop happens after the parallel region.
  enum { ABORT_NONE = 0, ABORT_CHOL, ABORT_INTERRUPT, ABORT_IWISHART,
         ABORT_OMEGA };
  int abort_code = ABORT_NONE;
  int abort_iter = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::vec d(p);   // thread-local residuals w_i - mu_i
    arma::vec e(p);   // master-only scratch for the Sigma step

    for (int r = 0; r < R; ++r) {
      // --- (a) Latent sweep: w_ij | w_i,-j, beta, Sigma, y_i -----------------
      // Univariate truncated-normal Gibbs steps using the precision-based
      // conditionals of N_p(mu_i, Sigma): m_ij = mu_ij -
      // sum_{k != j} Ratio(k, j) (w_ik - mu_ik), sd tau_j.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        Xoshiro256pp rng = make_stream(useed, r, static_cast<uint64_t>(i));
        const int yi = y[i];
        double* wcol = W.colptr(i);
        const double* mucol = Mu.colptr(i);
        for (int j = 0; j < p; ++j) {
          d(j) = wcol[j] - mucol[j];
        }
        for (int j = 0; j < p; ++j) {
          const double* rcol = Ratio.colptr(j);
          double rdot = 0.0;
          for (int k = 0; k < p; ++k) {
            rdot += rcol[k] * d(k);
          }
          const double m = mucol[j] - (rdot - d(j));   // Ratio(j, j) == 1

          // Truncation region implied by y_i, given the other components.
          double lo = -arma::datum::inf;
          double hi = arma::datum::inf;
          if (yi == j + 1) {
            lo = 0.0;                              // chosen: beat base and all others
            for (int k = 0; k < p; ++k) {
              if (k != j && wcol[k] > lo) lo = wcol[k];
            }
          } else if (yi == 0) {
            hi = 0.0;                              // base chosen: below 0
          } else {
            hi = wcol[yi - 1];                     // below the chosen component
          }

          wcol[j] = rtruncnorm(rng, m, tau(j), lo, hi);
          d(j) = wcol[j] - mucol[j];
        }
      }   // implied barrier: sweep complete before the beta step

      // --- (b) beta | w, Sigma: N(btilde, Q^{-1}) ----------------------------
      // Q = A + sum_i X_i' Omega X_i, btilde = Q^{-1} (A beta_bar +
      // sum_i X_i' Omega w_i); drawn via Cholesky solves, no inverse formed.
      // Master-only; the X' (Omega w) products are hand-rolled in a fixed
      // order so results do not depend on the BLAS or the thread count.
#ifdef _OPENMP
#pragma omp master
#endif
      {
        Xoshiro256pp rng_beta = make_stream(useed, r, static_cast<uint64_t>(N));
        Q = A;
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < p; ++k) {
            Q += Omega(j, k) * G[static_cast<size_t>(j) * p + k];
          }
        }
        for (int i = 0; i < N; ++i) {              // Vm = Omega W
          const double* wcol = W.colptr(i);
          double* vcol = Vm.colptr(i);
          for (int j = 0; j < p; ++j) {
            const double* ocol = Omega.colptr(j);  // = row j (Omega symmetric)
            double s = 0.0;
            for (int k = 0; k < p; ++k) s += ocol[k] * wcol[k];
            vcol[j] = s;
          }
        }
        const double* Vmv = Vm.memptr();           // bvec = Ab + X' vec(Vm)
        for (int k = 0; k < K; ++k) {
          const double* xc = X.colptr(k);
          double s = 0.0;
          for (arma::uword t = 0; t < NP; ++t) s += xc[t] * Vmv[t];
          bvec(k) = Ab(k) + s;
        }
        if (!arma::chol(U, Q)) {
          abort_code = ABORT_CHOL;
          abort_iter = r;
        } else {
          btilde = arma::solve(arma::trimatl(U.t()), bvec);
          btilde = arma::solve(arma::trimatu(U), btilde);
          for (int k = 0; k < K; ++k) {
            z(k) = rng_beta.rnorm();
          }
          beta = btilde + arma::solve(arma::trimatu(U), z);
        }
      }
#ifdef _OPENMP
#pragma omp barrier
#endif
      if (abort_code != ABORT_NONE) break;

      // --- Mu refresh: mu_i = X_i beta with the new beta ---------------------
      // Per-observation summation order is fixed, so the result is bitwise
      // independent of the number of threads.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < N; ++i) {
        double* mucol = Mu.colptr(i);
        const arma::uword row0 = static_cast<arma::uword>(i) * p;
        for (int j = 0; j < p; ++j) {
          double s = 0.0;
          for (int k = 0; k < K; ++k) s += Xv[k * NP + row0 + j] * beta(k);
          mucol[j] = s;
        }
      }   // implied barrier: Mu complete before the Sigma step

      // --- (c) Sigma | w, beta: IW(nu + N, V + S), S = sum_i eps_i eps_i' ----
#ifdef _OPENMP
#pragma omp master
#endif
      {
        Xoshiro256pp rng_sigma = make_stream(useed, r, static_cast<uint64_t>(N) + 1);
        S.zeros();
        for (int i = 0; i < N; ++i) {
          const double* wcol = W.colptr(i);
          const double* mucol = Mu.colptr(i);
          for (int j = 0; j < p; ++j) e(j) = wcol[j] - mucol[j];
          for (int j = 0; j < p; ++j) {
            for (int k = 0; k <= j; ++k) S(j, k) += e(j) * e(k);
          }
        }
        for (int j = 0; j < p; ++j) {
          for (int k = j + 1; k < p; ++k) S(j, k) = S(k, j);
        }
        Vpost = V + S;
        Vpost = 0.5 * (Vpost + Vpost.t());
        if (!riwishart_nothrow(rng_sigma, nu + N, Vpost, Sigma)) {
          abort_code = ABORT_IWISHART;
          abort_iter = r;
        } else if (!arma::inv_sympd(Omega, Sigma)) {
          abort_code = ABORT_OMEGA;
          abort_iter = r;
        } else {
          for (int j = 0; j < p; ++j) {
            const double ojj = Omega(j, j);
            tau(j) = std::sqrt(1.0 / ojj);
            for (int k = 0; k < p; ++k) Ratio(k, j) = Omega(k, j) / ojj;
          }

          // --- Record and housekeeping (master, inside the region) ----------
          if (r >= burn && (r - burn) % thin == 0) {
            const int row = (r - burn) / thin;
            betadraw.row(row) = beta.t();
            int idx = 0;
            for (int a = 0; a < p; ++a) {
              for (int b = 0; b <= a; ++b) {
                sigmadraw(row, idx++) = Sigma(a, b);
              }
            }
          }
          if (trace > 0 && (r + 1) % trace == 0) {
            Rprintf("mnp_gibbs: iteration %d / %d\n", r + 1, R);
          }
          if ((r + 1) % 100 == 0 && pending_interrupt()) {
            abort_code = ABORT_INTERRUPT;
            abort_iter = r;
          }
        }
      }
#ifdef _OPENMP
#pragma omp barrier
#endif
      if (abort_code != ABORT_NONE) break;
    }
  }

  switch (abort_code) {
    case ABORT_CHOL:
      Rcpp::stop("mnp_gibbs: posterior precision of beta is not positive "
                 "definite at iteration %d.", abort_iter + 1);
    case ABORT_IWISHART:
      Rcpp::stop("mnp_gibbs: Sigma scale matrix is not positive definite at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_OMEGA:
      Rcpp::stop("mnp_gibbs: Sigma draw is numerically singular at "
                 "iteration %d.", abort_iter + 1);
    case ABORT_INTERRUPT:
      Rcpp::stop("mnp_gibbs: interrupted by user at iteration %d.",
                 abort_iter + 1);
    default:
      break;
  }

  return Rcpp::List::create(
    Rcpp::Named("betadraw") = betadraw,
    Rcpp::Named("sigmadraw") = sigmadraw,
    Rcpp::Named("R_keep") = R_keep
  );
}

// ============================================================================
// Internal wrappers exposing the sampling primitives to R for testing
// (reached via choicer:::; not part of the package API).
// ============================================================================

// [[Rcpp::export(rng = false)]]
arma::vec rnorm_cpp(const int n, const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  arma::vec out(n);
  for (int i = 0; i < n; ++i) out(i) = rng.rnorm();
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::vec rgamma_cpp(const int n, const double a, const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  arma::vec out(n);
  for (int i = 0; i < n; ++i) out(i) = rng.rgamma(a);
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::vec rtruncnorm_cpp(const int n, const double mu, const double sigma,
                         const double a, const double b, const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  arma::vec out(n);
  for (int i = 0; i < n; ++i) out(i) = rtruncnorm(rng, mu, sigma, a, b);
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::mat rmvnorm_cpp(const int n, const arma::vec& mu, const arma::mat& Sigma,
                      const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  arma::mat out(n, mu.n_elem);
  for (int i = 0; i < n; ++i) {
    out.row(i) = rmvnorm(rng, mu, Sigma).t();
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::mat rwishart_cpp(const double df, const arma::mat& S, const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  return rwishart(rng, df, S);
}

// [[Rcpp::export(rng = false)]]
arma::mat riwishart_cpp(const double df, const arma::mat& V, const double seed) {
  Xoshiro256pp rng(static_cast<uint64_t>(seed));
  return riwishart(rng, df, V);
}
