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
- previous win-builder R-devel (before the spelling dictionary):
  - 0 errors | 0 warnings | 1 note
  - `counterfactuals` in `DESCRIPTION` is a correctly spelled, standard econometric term.

The current local `--as-cran` run used a freshly rebuilt
`choicer_0.2.0.tar.gz` with the spelling and thread-cap fixes. It included CRAN
incoming-feasibility, HTML-manual/math-rendering, examples including
`--run-donttest`, tests, and vignette rebuilding. A fresh win-builder R-devel
check should be run before resubmission; `.aspell/defaults.R` registers the
legitimate term `counterfactuals` for the incoming DESCRIPTION spell check.

Tests cap OpenMP and data.table at two threads in line with CRAN policy:
environment variables are set before loading the package in `tests/testthat.R`,
and authoritative runtime caps are set in `tests/testthat/setup.R`.

## Downstream dependencies

None.
