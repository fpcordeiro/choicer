## Resubmission

This is a resubmission of choicer 0.1.0. The previous submission produced
two NOTEs on the CRAN incoming pre-tests:

1. *Possibly misspelled words in DESCRIPTION* (Hessians, elasticities).
   These are technical terms, not misspellings:
   - "Hessians" is the plural of the Hessian matrix, a standard term
     in numerical optimization.
   - "elasticities" is standard econometric terminology (plural of
     elasticity).
   ("BLP" was removed from DESCRIPTION in this resubmission; it is
   still referenced in the documentation for the `blp()` method,
   which now cites Berry, Levinsohn, and Pakes (1995)
   <doi:10.2307/2171802>.)

2. *CPU time 14.5 times elapsed time when running tests* on
   r-devel-linux-x86_64-debian-gcc. The package uses OpenMP for
   parallelization. `tests/testthat.R` now sets `OMP_THREAD_LIMIT=2`
   and `OMP_NUM_THREADS=2` before loading the package, so tests use
   at most two cores in line with CRAN policy.

## Test environments

- local macOS, R 4.6.0
- win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission of a new package (first-time submission).

## Downstream dependencies

None — this is a new package.
