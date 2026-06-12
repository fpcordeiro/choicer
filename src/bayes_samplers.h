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
// the hard case is always a lower bound in the *upper* tail. Three regimes,
// chosen to keep the expected cost low without ever leaving exact sampling:
//
//   * lo < 0.45 and wide (hi - lo > 0.5): naive rejection — draw N(0, 1)
//     until it lands in (lo, hi). Acceptance is at least 1 - Phi(0.45) ~ 0.33
//     one-sided and Phi(0.95) - Phi(0.45) ~ 0.16 two-sided, and the polar
//     normal draw is ~5x cheaper than a pnorm/qnorm round trip. This is the
//     regime the MNP latent sweep hits almost always.
//   * lo >= 0.45 and wide (or ultra-deep, lo > 30): Robert (1995)
//     exponential rejection with rate lambda = (lo + sqrt(lo^2 + 4)) / 2,
//     acceptance probability exp(-(z - lambda)^2 / 2); exact arbitrarily far
//     into the tail.
//   * narrow intervals (hi - lo <= 0.5, lo <= 30): inverse CDF on the
//     upper-tail scale (R::pnorm/R::qnorm with lower_tail = 0), where
//     rejection methods would stall; upper-tail probabilities at lo <= 30.5
//     stay far above double underflow, so u and qnorm remain accurate.
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

  const bool wide = (hi - lo) > 0.5;   // true whenever either bound is Inf
  double z;
  if (wide && lo < 0.45) {
    // Naive rejection: cheapest exact method when the region holds mass.
    do {
      z = rng.rnorm();
    } while (z <= lo || z >= hi);
  } else if (wide || lo > 30.0) {
    // Tail case: Robert (1995) exponential rejection on (lo, hi).
    const double lambda = 0.5 * (lo + std::sqrt(lo * lo + 4.0));
    for (;;) {
      z = lo + rng.rexp() / lambda;
      if (z > hi) continue;
      const double diff = z - lambda;
      if (std::log(rng.runif()) <= -0.5 * diff * diff) break;
    }
  } else {
    // Narrow interval: inverse CDF on upper-tail probabilities.
    // P(Z > lo) > P(Z > hi); u is uniform on (P(Z > hi), P(Z > lo)).
    const double pu_lo = R::pnorm(lo, 0.0, 1.0, 0, 0);
    const double pu_hi = R::pnorm(hi, 0.0, 1.0, 0, 0);
    const double u = pu_hi + (pu_lo - pu_hi) * rng.runif();
    z = R::qnorm(u, 0.0, 1.0, 0, 0);
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
//
// The *_nothrow variants report failure through the return value instead of
// Rcpp::stop, so they are callable from inside an OpenMP parallel region
// (where a C++ exception would cross the region boundary and terminate).
// ----------------------------------------------------------------------------
inline bool rwishart_nothrow(Xoshiro256pp& rng, const double df,
                             const arma::mat& S, arma::mat& out) {
  const int p = S.n_rows;
  if (df < p) return false;
  arma::mat L;
  if (!arma::chol(L, S, "lower")) return false;
  arma::mat T(p, p, arma::fill::zeros);
  for (int i = 0; i < p; ++i) {
    T(i, i) = std::sqrt(rng.rchisq(df - i));
    for (int j = 0; j < i; ++j) T(i, j) = rng.rnorm();
  }
  const arma::mat LT = L * T;
  out = LT * LT.t();
  return true;
}

inline arma::mat rwishart(Xoshiro256pp& rng, const double df,
                          const arma::mat& S) {
  const int p = S.n_rows;
  if (df < p) {
    Rcpp::stop("rwishart: degrees of freedom (%g) must be >= dimension (%d).",
               df, p);
  }
  arma::mat out;
  if (!rwishart_nothrow(rng, df, S, out)) {
    Rcpp::stop("rwishart: scale matrix is not positive definite.");
  }
  return out;
}

// ----------------------------------------------------------------------------
// Inverse-Wishart(df, V): draw W ~ Wishart(df, V^{-1}) and invert.
// The result is symmetrised to remove floating-point asymmetry.
// ----------------------------------------------------------------------------
inline bool riwishart_nothrow(Xoshiro256pp& rng, const double df,
                              const arma::mat& V, arma::mat& out) {
  arma::mat Vinv;
  if (!arma::inv_sympd(Vinv, V)) return false;
  arma::mat W;
  if (!rwishart_nothrow(rng, df, Vinv, W)) return false;
  if (!arma::inv_sympd(out, W)) return false;
  out = 0.5 * (out + out.t());
  return true;
}

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
