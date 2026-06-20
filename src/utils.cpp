#include "choicer.h"

// Helper to safely get environment variables
// Inspired by data.table
const char* mygetenv(const char* name, const char* default_val) {
  const char* value = std::getenv(name);
  return (value == nullptr) ? default_val : value;
}

Rcpp::String env_or_na(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') return Rcpp::String(NA_STRING);
  return Rcpp::String(value);
}

Rcpp::List collect_thread_info() {
#ifdef _OPENMP
  int active_threads = 1;
#pragma omp parallel
  {
#pragma omp master
    {
      active_threads = omp_get_num_threads();
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("openmp_enabled") = true,
    Rcpp::Named("_OPENMP") = _OPENMP,
    Rcpp::Named("omp_get_num_threads") = active_threads,
    Rcpp::Named("omp_get_max_threads") = omp_get_max_threads(),
    Rcpp::Named("omp_get_num_procs") = omp_get_num_procs(),
    Rcpp::Named("omp_get_thread_limit") = omp_get_thread_limit(),
    Rcpp::Named("OMP_THREAD_LIMIT") = env_or_na("OMP_THREAD_LIMIT"),
    Rcpp::Named("OMP_NUM_THREADS") = env_or_na("OMP_NUM_THREADS")
  );
#else
  return Rcpp::List::create(
    Rcpp::Named("openmp_enabled") = false,
    Rcpp::Named("_OPENMP") = NA_INTEGER,
    Rcpp::Named("omp_get_num_threads") = NA_INTEGER,
    Rcpp::Named("omp_get_max_threads") = 1,
    Rcpp::Named("omp_get_num_procs") = NA_INTEGER,
    Rcpp::Named("omp_get_thread_limit") = NA_INTEGER,
    Rcpp::Named("OMP_THREAD_LIMIT") = env_or_na("OMP_THREAD_LIMIT"),
    Rcpp::Named("OMP_NUM_THREADS") = env_or_na("OMP_NUM_THREADS")
  );
#endif
}

//' Query choicer OpenMP thread settings
//'
//' @return A list with OpenMP availability, active/max thread settings, CPU
//'   thread capacity reported by OpenMP, thread limits, and relevant
//'   environment variables.
//' @export
// [[Rcpp::export]]
Rcpp::List thread_info() {
  return collect_thread_info();
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

//' Set the number of OpenMP threads used by choicer
//'
//' @param n_threads Positive integer number of threads.
//' @return Invisibly returns `NULL`.
//' @export
// [[Rcpp::export]]
void set_num_threads(int n_threads) {
  if (n_threads < 1) Rcpp::stop("`n_threads` must be a positive integer.");
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#else
  Rcpp::warning("This installation was not compiled with OpenMP support.");
#endif
}

Rcpp::IntegerVector compute_prefix_sum(const Rcpp::IntegerVector& M) {
  int N = M.size();
  Rcpp::IntegerVector S(N + 1);
  S[0] = 0;
  std::partial_sum(M.begin(), M.end(), S.begin() + 1);
  return S;
}
