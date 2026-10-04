#ifndef CHOICER_HPP
#define CHOICER_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

// ---------------------------------------------------------------------------
// Primary-thread block.
//
// OpenMP 5.1 deprecated `master` in favour of `masked`; GCC 16 warns on the
// old spelling. `masked` with no `filter` clause is semantically identical --
// the block runs on the primary thread (filter(0)) with no implied barrier --
// so we simply emit whichever spelling the compiler in hand expects. Without
// OpenMP it expands to nothing: the sole thread is already the primary one.
#ifndef _OPENMP
#  define CHOICER_OMP_MASKED
#elif _OPENMP >= 202011
#  define CHOICER_OMP_MASKED _Pragma("omp masked")
#else
#  define CHOICER_OMP_MASKED _Pragma("omp master")
#endif

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