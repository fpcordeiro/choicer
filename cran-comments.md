## choicer 0.2.0

This is the first update release after choicer 0.1.0 was accepted to CRAN.
It adds hierarchical Bayesian models (HMNL/HMNP) with convergence
diagnostics, robust and cluster-robust standard errors, and corrects a
documentation claim about thread-count reproducibility of MCMC draws.

## Test environments

- local macOS, R 4.6.0
- win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

Tests cap OpenMP at two threads (`OMP_THREAD_LIMIT=2`, `OMP_NUM_THREADS=2`
in `tests/testthat.R`) in line with CRAN policy.

## Downstream dependencies

None.
