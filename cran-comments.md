## choicer 0.2.0

This is the first update release after choicer 0.1.0 was accepted to CRAN.
It adds hierarchical Bayesian models (HMNL/HMNP) with convergence
diagnostics, robust and cluster-robust standard errors, and corrects a
documentation claim about thread-count reproducibility of MCMC draws.

## Test environments

- local macOS, R 4.6.0

## R CMD check results

Current local release-candidate results:

- `R CMD check`: 0 errors | 0 warnings | 0 notes
- `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes

Both checks were run on 2026-07-12 against the version-stamped
`choicer_0.2.0.tar.gz`. The `--as-cran` run included CRAN incoming-feasibility,
HTML-manual/math-rendering, examples including `--run-donttest`, tests, and
vignette rebuilding.

Tests cap OpenMP at two threads (`OMP_THREAD_LIMIT=2`, `OMP_NUM_THREADS=2`
in `tests/testthat.R`) in line with CRAN policy.

## Downstream dependencies

None.
