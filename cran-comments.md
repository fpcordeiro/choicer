## choicer 0.2.0

This is the first update release after choicer 0.1.0 was accepted to CRAN.
It adds hierarchical Bayesian models (HMNL/HMNP) with convergence
diagnostics, robust and cluster-robust standard errors, and corrects a
documentation claim about thread-count reproducibility of MCMC draws.

## Test environments

- local macOS Sequoia 15.7.7 (aarch64), R 4.6.0
- win-builder, Windows Server 2022 x64, R-devel (2026-07-11 r90235)

## R CMD check results

- local macOS: `R CMD check --as-cran choicer_0.2.0.tar.gz`
  - 0 errors | 0 warnings | 0 notes
- win-builder R-devel:
  - 0 errors | 0 warnings | 1 note
  - `counterfactuals` in `DESCRIPTION` is a correctly spelled, standard econometric term.

Both checks used the version-stamped `choicer_0.2.0.tar.gz`. The local
`--as-cran` run included CRAN incoming-feasibility, HTML-manual/math-rendering,
examples including `--run-donttest`, tests, and vignette rebuilding.

Tests cap OpenMP at two threads (`OMP_THREAD_LIMIT=2`, `OMP_NUM_THREADS=2`
in `tests/testthat.R`) in line with CRAN policy.

## Downstream dependencies

None.
