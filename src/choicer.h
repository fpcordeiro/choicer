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
  if (x.n_elem == 0) {
    return -arma::datum::inf;                   // log(sum(empty)) = log(0)
  }

  const double a = x.max();

  if (std::isnan(a))            return arma::datum::nan;  // propagate NaN
  if (a ==  arma::datum::inf)   return  arma::datum::inf; // any +inf -> +inf
  if (a == -arma::datum::inf)   return -arma::datum::inf; // all -inf -> -inf

  const double s = arma::accu(arma::exp(x - a)); // s >= 1 because at least one term is exp(0)=1
  return a + std::log(s);
}

#endif