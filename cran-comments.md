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
- win-builder, x86_64-w64-mingw32, R Under development (2026-08-31 r90457 ucrt)
- R-hub `gcc16` container (`ghcr.io/r-hub/containers/gcc16:latest`),
  R Under development (2026-09-01 r90464), GCC 16.2.1 — the toolchain that
  reported the WARNING

## R CMD check results

- local macOS: `R CMD check --as-cran choicer_0.2.1.tar.gz`
  - 0 errors | 0 warnings | 0 notes
- win-builder R-devel: Status: OK
  - 0 errors | 0 warnings | 0 notes
- R-hub gcc16 (`--as-cran`, with `rcmdcheck(error_on = "warning")` so that a
  compiler warning fails the run rather than passing as a check WARNING):
    Status: OK
  - 0 errors | 0 warnings | 0 notes
  - `checking whether package 'choicer' can be installed ... OK`, with no
    "Found the following significant warnings" section; the reported
    `master`-construct diagnostic does not appear anywhere in the build log

The local run included CRAN incoming-feasibility, HTML-manual/math-rendering,
examples including `--run-donttest`, tests, and vignette rebuilding.
`.aspell/defaults.R` registers the legitimate econometric term
`counterfactuals` for the incoming DESCRIPTION spell check.

GCC sets `_OPENMP` to 201511 through the 15 series and to 202111 in GCC 16, so
GCC 16 takes the `masked` branch while GCC 15 and earlier — which do not warn
— continue to compile `master`. Both branches were exercised before
submission: the affected translation units compile warning-free under `-Wall`
with clang at `-fopenmp-version=45` (selecting `master`) and at
`-fopenmp-version=51`/`52` (selecting `masked`), under GCC 15, and with OpenMP
disabled entirely, and in every configuration the guarded block executes on
the primary thread. The full test suite passes (`FAIL 0 | WARN 0 | SKIP 80 |
PASS 1718` under `R CMD check`), including the MNP, HMNL, and HMNP sampler
reproducibility and parameter-recovery tests.

Tests cap OpenMP and data.table at two threads in line with CRAN policy:
environment variables are set before loading the package in `tests/testthat.R`,
and authoritative runtime caps are set in `tests/testthat/setup.R`.

## Downstream dependencies

None.
