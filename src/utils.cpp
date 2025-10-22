#include "choicer.h"

// Helper to safely get environment variables
// Inspired by data.table
const char* mygetenv(const char* name, const char* default_val) {
  const char* value = std::getenv(name);
  return (value == nullptr) ? default_val : value;
}

// Inspired by data.table
// [[Rcpp::export]]
void get_num_threads() {  
  #ifndef _OPENMP
    Rprintf("This installation has not been compiled with OpenMP support.\n");
  #else
  #pragma omp parallel
  {
    #pragma omp master
    {
      Rprintf("  OpenMP version (_OPENMP)       %d\n", _OPENMP);
      Rprintf("  omp_get_num_threads()          %d\n", omp_get_num_threads());
      Rprintf("  omp_get_num_procs()            %d\n", omp_get_num_procs());
      Rprintf("  omp_get_thread_limit()         %d\n", omp_get_thread_limit());
      Rprintf("  omp_get_max_threads()          %d\n", omp_get_max_threads());
      Rprintf("  OMP_THREAD_LIMIT               %s\n", mygetenv("OMP_THREAD_LIMIT", "unset"));
      Rprintf("  OMP_NUM_THREADS                %s\n", mygetenv("OMP_NUM_THREADS", "unset"));
    }
  }
  #endif
}

// [[Rcpp::export]]
void set_num_threads(int n_threads) {
  #ifndef _OPENMP
    Rprintf("This installation has not been compiled with OpenMP support.\n");
  #else
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  #endif
}

Rcpp::IntegerVector compute_prefix_sum(const Rcpp::IntegerVector& M) {
  int N = M.size();
  Rcpp::IntegerVector S(N + 1);
  S[0] = 0;
  std::partial_sum(M.begin(), M.end(), S.begin() + 1);
  return S;
}
