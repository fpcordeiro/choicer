# HB sampler profiling harness (run 20260706-hb-profiling)

Evidence-gathering harness for `hmnl_gibbs` (src/hmnlogit.cpp) and
`hmnp_gibbs` (src/hmnprobit.cpp) at healthcare-provider-choice scale
(N ~ 30,000, J ~ 200, K = P = 10, R = 10,000 / burn = 5,000 / 2 sequential
chains). **No package source is modified** -- see `patches/hb_profile.patch`
and `00_build_instrumented.R`.

Phase 1 (current): harness built, smoke-tested, bitwise-gated. The full grid
in `configs/proposed_grid.csv` has **NOT** been run -- it requires user
approval (see the run's `harness-summary.md`).

## Files

| File | Purpose |
|---|---|
| `patches/hb_profile.patch` | Unified diff adding `#ifdef HB_PROFILE` timer blocks to `src/hb_internal.h`, `src/hmnlogit.cpp`, `src/hmnprobit.cpp`, and `-DHB_PROFILE` to `src/Makevars(.win)`. Applied ONLY to an rsync'd scratch copy, never to this repo. |
| `00_build_instrumented.R` | rsyncs the repo source to a scratch dir, applies the patch, `R CMD INSTALL`s into a throwaway library. |
| `01_bitwise_gate.R` | Bitwise sanity gate: runs a tiny fixed-seed config through a non-instrumented reference build and the instrumented build (separate subprocesses) and asserts identical draws; also a timer-overhead check. |
| `10_run_config.R` | Per-config WORKER: simulates data, prepares it, calls `hmnl_gibbs`/`hmnp_gibbs` directly (bypassing `run_hmnlogit`/`run_hmnprobit` to reach the `profile` list element), writes one JSON result. Meant to be launched as a fresh subprocess under `/usr/bin/time -l`. |
| `20_harness.R` | Config-driven orchestrator: reads a grid CSV, launches `10_run_config.R` per row wrapped in `/usr/bin/time -l`, parses peak RSS, appends to `results/hb_profiling_results.csv`. |
| `30_analyze.R` | Fits log-log scaling slopes per block over the N/J sweeps and extrapolates to the target scale + full feasibility budget. NOT yet run (grid not executed). |
| `configs/smoke_grid.csv` | The ONE-config-per-model smoke test allowed in Phase 1 (N=500, J=20, K=10, P=10, R=100). |
| `configs/build_proposed_grid.R`, `configs/proposed_grid*.csv` | The proposed full grid (N sweep, J sweep, thread points, target-scale anchor, memory probe, run_mnprobit contrast) -- generated, NOT executed. |
| `results/hb_profiling_results.csv` | Accumulated per-config results (currently: smoke configs only). |
| `results/timer_overhead.csv` | HB_PROFILE timer overhead measurement. |
| `results/bitwise_gate_result.rds` | Saved gate comparison + overhead table. |

## Reproducing the instrumented build

```r
Rscript _benchmarks/hb_profiling/00_build_instrumented.R
# installs into .../tmp/hb_prof_build/templib by default (outside the repo)
```

## Running the (approved) smoke config

```r
Rscript _benchmarks/hb_profiling/20_harness.R \
  --grid _benchmarks/hb_profiling/configs/smoke_grid.csv \
  --lib-dir /path/to/templib
```

## Running the bitwise gate

```r
Rscript _benchmarks/hb_profiling/01_bitwise_gate.R \
  --lib-dir /path/to/templib --ref-lib-dir /path/to/templib_ref
```

`--ref-lib-dir` is a plain (non-instrumented) `R CMD INSTALL` of the SAME
source tree into a separate throwaway library -- used instead of "the
currently installed choicer" because the system-installed copy at this
machine's default library path predates the HMNL/HMNP feature.

## Safety notes

- Every `#ifdef HB_PROFILE` timer reads the wall clock on the thread that
  would have executed that code anyway; it never touches an RNG stream or
  reorders computation (verified by the bitwise gate).
- `results/tmp/` (worker JSON/RDS/`time` logs) is gitignored -- ephemeral,
  reproducible by re-running the scripts above.
- The scratch build directories (`choicer_src`, `choicer_src_ref`,
  `templib`, `templib_ref`) live entirely outside this repo, under the run's
  workspace `tmp/hb_prof_build/`.
