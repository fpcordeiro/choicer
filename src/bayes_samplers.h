#ifndef CHOICER_BAYES_SAMPLERS_HPP
#define CHOICER_BAYES_SAMPLERS_HPP

#include "choicer.h"
#include "rng.h"

// ============================================================================
// Header-only sampling primitives for the Bayesian (MCMC) engines.
//
// Every function draws exclusively from the Xoshiro256pp stream passed in,
// never from R's RNG, so rtruncnorm() is safe inside OpenMP regions provided
// each task uses its own stream (see make_stream() in rng.h). R::pnorm and
// R::qnorm below are pure Rmath distribution functions — thread-safe; only
// R's *RNG* is off-limits in threaded code.
//
// rmvnorm / rwishart / riwishart may call Rcpp::stop and must therefore only
// be called from the master thread (they are used for the serial beta and
// Sigma conditionals).
//
// Inverse-Wishart convention (matches bayesm): Sigma ~ IW(df, V) has density
// proportional to |Sigma|^{-(df+p+1)/2} exp(-tr(V Sigma^{-1})/2) and prior
// mean V / (df - p - 1).
// ============================================================================

// ----------------------------------------------------------------------------
// Truncated univariate normal on (a, b); either bound may be +/-Inf.
//
// After standardising, the interval is mirrored when it lies below zero so
// the hard case is always a lower bound far in the *upper* tail. Central
// case: inverse CDF on the upper-tail scale (R::pnorm/R::qnorm with
// lower_tail = 0), which keeps full precision for bounds deep in the tail
// where lower-tail probabilities round to 1. Tail case (lo >= 4): Robert
// (1995) exponential rejection with rate lambda = (lo + sqrt(lo^2 + 4)) / 2,
// acceptance probability exp(-(z - lambda)^2 / 2).
// ----------------------------------------------------------------------------
inline double rtruncnorm(Xoshiro256pp& rng, const double mu, const double sigma,
                         const double a, const double b) {
  double lo = (a - mu) / sigma;
  double hi = (b - mu) / sigma;

  // Mirror intervals lying below zero into the upper tail.
  bool flipped = false;
  if (hi < 0.0) {
    const double tmp = lo;
    lo = -hi;
    hi = -tmp;
    flipped = true;
  }

  double z;
  if (lo < 4.0) {
    // Central case: inverse CDF on upper-tail probabilities.
    // P(Z > lo) > P(Z > hi); u is uniform on (P(Z > hi), P(Z > lo)).
    const double pu_lo = R::pnorm(lo, 0.0, 1.0, 0, 0);
    const double pu_hi = R::pnorm(hi, 0.0, 1.0, 0, 0);
    const double u = pu_hi + (pu_lo - pu_hi) * rng.runif();
    z = R::qnorm(u, 0.0, 1.0, 0, 0);
  } else {
    // Tail case: Robert (1995) exponential rejection on (lo, hi).
    const double lambda = 0.5 * (lo + std::sqrt(lo * lo + 4.0));
    for (;;) {
      z = lo + rng.rexp() / lambda;
      if (z > hi) continue;
      const double diff = z - lambda;
      if (std::log(rng.runif()) <= -0.5 * diff * diff) break;
    }
  }

  return mu + sigma * (flipped ? -z : z);
}

// ----------------------------------------------------------------------------
// Multivariate normal N(mu, Sigma) via the lower Cholesky factor.
// ----------------------------------------------------------------------------
inline arma::vec rmvnorm(Xoshiro256pp& rng, const arma::vec& mu,
                         const arma::mat& Sigma) {
  arma::mat L;
  if (!arma::chol(L, Sigma, "lower")) {
    Rcpp::stop("rmvnorm: covariance matrix is not positive definite.");
  }
  arma::vec z(mu.n_elem);
  for (arma::uword i = 0; i < z.n_elem; ++i) z(i) = rng.rnorm();
  return mu + L * z;
}

// ----------------------------------------------------------------------------
// Wishart(df, S) via the Bartlett decomposition: W = (L T)(L T)' where
// L = chol(S, "lower") and T is lower triangular with T(i, i) =
// sqrt(chisq(df - i)) (0-based i) and iid N(0, 1) below the diagonal.
// Requires df >= p; non-integer df allowed.
// ----------------------------------------------------------------------------
inline arma::mat rwishart(Xoshiro256pp& rng, const double df,
                          const arma::mat& S) {
  const int p = S.n_rows;
  if (df < p) {
    Rcpp::stop("rwishart: degrees of freedom (%g) must be >= dimension (%d).",
               df, p);
  }
  arma::mat L;
  if (!arma::chol(L, S, "lower")) {
    Rcpp::stop("rwishart: scale matrix is not positive definite.");
  }
  arma::mat T(p, p, arma::fill::zeros);
  for (int i = 0; i < p; ++i) {
    T(i, i) = std::sqrt(rng.rchisq(df - i));
    for (int j = 0; j < i; ++j) T(i, j) = rng.rnorm();
  }
  const arma::mat LT = L * T;
  return LT * LT.t();
}

// ----------------------------------------------------------------------------
// Inverse-Wishart(df, V): draw W ~ Wishart(df, V^{-1}) and invert.
// The result is symmetrised to remove floating-point asymmetry.
// ----------------------------------------------------------------------------
inline arma::mat riwishart(Xoshiro256pp& rng, const double df,
                           const arma::mat& V) {
  arma::mat Vinv;
  if (!arma::inv_sympd(Vinv, V)) {
    Rcpp::stop("riwishart: scale matrix is not positive definite.");
  }
  const arma::mat W = rwishart(rng, df, Vinv);
  arma::mat out;
  if (!arma::inv_sympd(out, W)) {
    Rcpp::stop("riwishart: Wishart draw is numerically singular.");
  }
  return 0.5 * (out + out.t());
}

#endif
