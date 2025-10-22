#ifndef CHOICER_HPP
#define CHOICER_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

// Function Declarations ------------------------------------------------------
Rcpp::IntegerVector compute_prefix_sum(const Rcpp::IntegerVector& M);

// Inline Function Definitions ------------------------------------------------
inline double logSumExp(const arma::vec& x) {
  double a = x.max();
  if (!arma::is_finite(a)) {
    return arma::datum::nan;
  }
  return a + std::log(arma::sum(arma::exp(x - a)));
}

#endif