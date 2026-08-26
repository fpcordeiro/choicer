## choicer 0.2.1

This patch release fixes the WARNING reported for choicer 0.2.0 on the
`r-devel-linux-x86_64-fedora-gcc` and `r-devel-linux-x86_64-debian-gcc` check
flavors:

```
hmnlogit.cpp:554:19: warning: 'master' construct deprecated since OpenMP 5.1, use 'masked'
```

OpenMP 5.1 deprecated the `master` construct in favor of `masked`, and GCC 16
warns on the old spelling. The seven affected directives (in `hmnlogit.cpp`,
`hmnprobit.cpp`, `mnprobit.cpp`, and `utils.cpp`) now go through a
`CHOICER_OMP_MASKED` macro in `src/choicer.h` that emits `masked` when
`_OPENMP >= 202011` and `master` otherwise, so the package compiles without
warnings on newer toolchains while remaining correct on OpenMP 4.5 compilers
and on builds without OpenMP.

`masked` with no `filter` clause is semantically identical to `master` — the
block executes on the primary thread with no implied barrier — so the Gibbs
samplers' barrier structure, fixed-order reductions, and posterior draws are
unchanged. There are no user-visible behavior, API, or numerical changes in
this release.

## Test environments

- local macOS Sequoia 15.7.7 (aarch64), R 4.6.0, Apple clang
- TODO BEFORE SUBMITTING: win-builder, Windows Server 2022 x64, R-devel
- TODO BEFORE SUBMITTING: R-hub `gcc16` container
  (`ghcr.io/r-hub/containers/gcc16:latest`) — the toolchain that reported the
  WARNING; this is the check that confirms the fix on the affected flavor

## R CMD check results

- local macOS: `R CMD check --as-cran choicer_0.2.1.tar.gz`
  - 0 errors | 0 warnings | 0 notes

The run included CRAN incoming-feasibility, HTML-manual/math-rendering,
examples including `--run-donttest`, tests, and vignette rebuilding.
`.aspell/defaults.R` registers the legitimate econometric term
`counterfactuals` for the incoming DESCRIPTION spell check.

Apple clang defaults to `-fopenmp-version=51`, so `_OPENMP` is 202011 on the
development machine and the local check above compiles the **new `masked`
branch** of the guard, not the old one. That build passes the full suite
(`FAIL 0 | WARN 0 | SKIP 80 | PASS 1718` under `R CMD check`; 2045 assertions
with no skips when run outside the check environment), including the MNP,
HMNL, and HMNP sampler reproducibility and parameter-recovery tests.

Clang does not emit the deprecation diagnostic, so the reported warning itself
could not be reproduced locally. The other branch of the guard was exercised
separately: the affected translation units compile warning-free under `-Wall`
with clang at `-fopenmp-version=45` (`_OPENMP` 201511, selecting `master`),
under GCC 15, and with OpenMP disabled entirely, and in every configuration
the guarded block executes on thread 0. GCC sets `_OPENMP` to 201511 through
the 15 series and to 202111 in GCC 16, so GCC 16 takes the `masked` branch and
no longer warns, while GCC 15 and earlier — which do not warn — continue to
compile `master`.

Tests cap OpenMP and data.table at two threads in line with CRAN policy:
environment variables are set before loading the package in `tests/testthat.R`,
and authoritative runtime caps are set in `tests/testthat/setup.R`.

## Downstream dependencies

None.
