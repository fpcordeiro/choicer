## Resubmission

This is a resubmission of choicer 0.1.0. The previous submission produced
two NOTEs on the CRAN incoming pre-tests:

1. *Possibly misspelled words in DESCRIPTION* (BLP, Hessians, elasticities).
   These are technical terms, not misspellings:
   - "BLP" is a proper-noun acronym for the Berry-Levinsohn-Pakes
     contraction (Berry, Levinsohn, and Pakes, 1995, Econometrica).
   - "Hessians" is the plural of the Hessian matrix, a standard term
     in numerical optimization.
   - "elasticities" is standard econometric terminology (plural of
     elasticity).

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
