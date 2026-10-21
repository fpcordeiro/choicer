# Population-scale benchmarks

Development only, excluded from the package build. These scripts measure the
mixed logit's compiled kernels, its fits and the data preparations on
synthetic data at the scale of health-claims and school-census records
(10^6 to 10^8 stacked rows, hundreds to thousands of alternatives), compare two
builds of choicer, and run two measurements meant for a many-core server.

They produced the figures of the scale audit's report and of the NEWS section
"Population-scale mixed logit".

## Data

`gen_data.cpp` builds data directly in the prepared (stacked) layout that the
kernels take, with no call to `prepare_*_data()`, deterministically for each
unit (seed 1). The designs (`gen_regime()` in `scale_common.R`) have 10 fixed
coefficients and 3 correlated normal random coefficients with means:

| Regime | Decision makers | Alternatives | Choice set | Situations each | Parameters |
|---|---|---|---|---|---|
| `claims` | patients | 500 hospitals, 50 regions | 20-200 of the region's, 90% kept | Pareto(1.1), up to 500 | 518 |
| `census` | students | 3,000 schools, 150 regions (about 2,600 appear) | 10-100 | 1 (85%), else 2 or 3 | about 2,620 |
| `census20k` | students | 20,000 schools | 10-100 | as `census` | about 17,000 |
| `small` | individuals | 4 | all 4 | 1 | 22 |

Sizes in use: claims with 1,500 / 15,000 / 150,000 patients (1.1·10^6 /
8.7·10^6 / 8.5·10^7 rows); census with 20,000 / 200,000 / 2,000,000 students
(1.0·10^6 / 9.8·10^6 / 9.8·10^7 rows); small with 250,000 / 2,500,000
individuals (10^6 / 10^7 rows). `fit_study.R` also drops the alternatives that
are never chosen, and the situations that offer one chosen wherever it is
offered, so that every constant has a finite estimate.

## Protocol

- **Builds in isolation.** Install each build from a copy of its tree into a
  library of its own (`R CMD INSTALL --library=<dir> <copy>`), never in the
  working tree, and point `R_LIBS` at it. Every script asserts that choicer was
  loaded from that library: an installed release in the site library would
  otherwise shadow the build.
- **One configuration per R process.**
- **Time.** Alternate the two builds (AB, BA, AB: three rounds or more), take
  the median of each run's timed calls, compare the builds by the median over
  rounds of their ratio, and report the load beside every figure. Run
  configurations of minutes per call once per build, and flag them. Time runs
  go without the shim (it costs every `free()` a lookup).
- **Memory.** On macOS, the malloc-counting shim (`memcount.c`, through the
  `Rmc` front end that `make_rmc.sh` writes) counts the heap blocks and bytes a
  call moves and its peak live heap, and samples the peak resident memory.
  These do not depend on the machine's load.
- **Load.** The machine's other processes distort wall times, most at many
  threads. The grid waits for a quiet machine and logs the load, CPU, available
  memory, swap and top processes of every run; a run during which the machine
  swapped is set aside. Under load, compare fits by their kernels' calls
  (`fit_study.R`) and model their time as calls times per-pass costs measured
  in alternated rounds.
- **Decompose a cost by removal**, timing a variant without the part in
  question, not by a sampling profile: one profile credited the inverse normal
  with 33.5% of a gradient evaluation where removing it saved 8.2%.
- **Identity.** At more than one thread, sums over decision makers vary in
  their last bits with the order in which threads take them; compare results
  at one thread where bit identity matters. The kernels' check values in the
  logs (the objective, or a sum of the result) agree across builds at the
  rounding of the work's summation order.

## Scripts

| Script | What it runs |
|---|---|
| `scale_common.R` | Shared helpers: the build's namespace check, `gen_regime()`, the shim's counters |
| `gen_data.cpp` | The synthetic data (Rcpp, compiled on first use into `output/cpp-cache`) |
| `bench_kernels.R` | One build, one configuration: the kernels behind the log-likelihood and gradient, Hessian, BHHH, scores, conditional tastes, prediction and the meats (`mode=time` or `mode=mem`) |
| `run_grid.sh` | The grid of the report: two builds alternated over the gradient grid, heavy panels, the other kernels at scale, data preparation and the per-call layout (written for macOS: `vm_stat`, `sysctl`, the shim) |
| `prep_mem.R` | One `prepare_mxl_data()`, `prepare_mnl_data()` or `prepare_nl_data()` call from long data |
| `fit_study.R` | `run_mxlogit()` fits (default, `polish`, `concentrate_out_asc`) with every kernel's calls counted; per-pass costs; a table against a Newton-certified reference |
| `budget_sweep.R` | Server run: the draw-batch budget (below) |
| `layout_cost.R` | Server run: the per-call layout at 10^9 rows (below) |
| `make_rmc.sh`, `memcount.c` | macOS: the malloc shim and the `Rmc` front end that preloads it |

```sh
# Two builds, A and B, each from a copy of its tree:
R CMD INSTALL --no-docs --library=/tmp/lib-A /tmp/choicer-A
R CMD INSTALL --no-docs --library=/tmp/lib-B /tmp/choicer-B
sh _benchmarks/scale/make_rmc.sh                    # macOS, for memory runs
MC_SHIM=$PWD/_benchmarks/scale/output/memcount.dylib \
  sh _benchmarks/scale/run_grid.sh A=/tmp/lib-A B=/tmp/lib-B
# One configuration:
R_LIBS=/tmp/lib-B Rscript _benchmarks/scale/bench_kernels.R regime=census U=200000 S=100 \
  threads=1,11 reps=3 what=grad,tastes
# Fits on one build:
R_LIBS=/tmp/lib-B Rscript _benchmarks/scale/fit_study.R fit regime=census U=200000 cell=lbfgs
```

Outputs go to `output/` (ignored by git), or to `CHOICER_SCALE_OUT`.

## Runs for a many-core server

### The draw-batch budget

The estimation kernels form a decision maker's rows-by-draws matrices a batch
of draws at a time, within a budget of 4 MiB per thread (`MXL_BATCH_BYTES`,
`src/mxlogit.cpp`): a unit of R rows is split into batches when
n_mats · 8 · R · S bytes exceed it (n_mats = 2 for the gradient and BHHH
kernels, 1 for the Hessian and the conditional tastes): at 4 MiB, a unit of
more than 2^18 / S rows for the gradient. Long panels split (a patient with
thousands of rows; at S = 1,000, any unit of more than 262 rows, so census sets
only their largest units, 0.4%); census sets at S = 100 and small choice sets
do not. A batch that stays in a core's cache is
read several times per draw; on an 11-core laptop (Apple M3 Pro) budgets of
1-16 MiB ran fastest on the heavy claims panel at S = 1,000, and 64-256 MiB
slower. Where many cores share a cache, the best budget may be smaller.

The budget is fixed at compile time, so `budget_sweep.R` builds one library per
budget from a source tree (with the machine's own Makevars and
`-DCHOICER_MXL_BATCH_BYTES=<bytes>`, from 1 KiB to 16 GiB), checks that each
carries its budget, and times the same kernels on the same data with each (the
package's tests assume the default budget, so don't run them against these
libraries):

```sh
Rscript _benchmarks/scale/budget_sweep.R --src=<choicer source tree> --quick  # setup check, minutes
Rscript _benchmarks/scale/budget_sweep.R --src=<choicer source tree>          # the sweep
```

The defaults: budgets of 0.25, 1, 2, 4, 8, 16 and 64 MiB; half and all of the
physical cores the process may use (on Linux counted from `lscpu` within its
CPU affinity, as `parallel::detectCores()` counts logical CPUs there); the
cases `heavy` (claims, 1,500 patients, S = 1,000), `claims`, `claims_s500` and
`claims_s1000` (15,000 patients, S = 100, 500 and 1,000), `census_s1000`
(200,000 students, S = 1,000), and the controls `census` (S = 100) and
`small`; the gradient, conditional tastes and BHHH kernels; three rounds of
three timed calls. A budget that batches every unit exactly as another does
runs once per case, and is reported with that budget's times (`same_as`). The
plan it prints gives each case's largest unit's share of the rows: where it
exceeds one thread's share at the largest thread count, as for `heavy`, a few
units set the wall time. Expect a few hours on 64 cores. A run resumes in the
same `--out` with `--no-build`, on the libraries already there, if its
configuration is unchanged; it refuses to rebuild into an `--out` that holds
results. Send back `machine.txt`, `compile.txt`, `summary.csv`, `raw.csv`,
the `install-*.log` files and `logs/` from the output directory. On a
multi-socket server, consider `numactl --interleave=all`, since the synthetic
data are written by one thread, and on SMT machines `OMP_PLACES=cores
OMP_PROC_BIND=spread`, so that threads counted in physical cores do not share
one. choicer's dependencies must be in `R_LIBS`, `R_LIBS_USER` or the site
library: the processes run without `~/.Rprofile`. A `~/.Renviron` that sets
`R_LIBS` overrides the libraries the script passes to its processes, which
its checks then refuse. `vm_add` (Linux) is a kernel's rise in resident memory
over the process's size before it ran; memory kept from earlier calls in the
same process is reused without raising it, so compare it across budgets only
within one kernel and thread count.

**Decision rule.** Keep 4 MiB unless another budget, at the largest thread
count, is more than 5% faster in every round on each long panel
(`claims_s500`, `claims_s1000`, `census_s1000`) where it batches at least a
fifth of the rows differently from 4 MiB (there must be one), and at most 2%
slower wherever else it batches at least a twentieth of the rows differently
(`heavy` included). Where it batches fewer rows differently, any effect is
within the timings' noise, and where it batches every row as 4 MiB does, its
times are 4 MiB's by construction; both are left out. The results must be complete, agree across budgets to 1e-10, and
come from thread teams of the size asked. A budget that depends on the cache
size would need that evidence on more than one machine.

### The per-call layout at 10^9 rows

Every estimation kernel call first builds the layout of the stacked design
(`MxlLayout`, `src/mxlogit.cpp`). It fills N + 1 eight-byte row offsets (N
choice situations; a panel adds U + 1 unit offsets), checks the alternatives'
codes in one pass over the rows (in parallel above 10^6 rows), and makes serial
passes over the situations and units: the offsets, the choices and the largest
unit, plus the unit data's check that a panel's weights are constant within
each decision maker. None of it depends on the parameters, so a fit could
build it once; a cache would pay only above about 5% of an evaluation. The
prediction kernels build the same offsets, and five of them also form a
base-utility vector of 8 bytes a row per call (`mxl_predict()` reuses its
output).

Calibration here (Apple M3 Pro, 5 performance + 6 efficiency cores; K_x =
K_w = 1, so the gradient is at its cheapest; the layout timed as the call that
`options(choicer.max_memory = 1)` stops right after the memory check, medians
of five; October 2026):

| Design | Situations | Rows | Layout, 1 / 11 threads | Gradient at 11 threads, S = 25 / 100 | Layout's share |
|---|---|---|---|---|---|
| 4 alternatives, cross-section | 2.5·10^6 | 10^7 | 12.9 / 12.4 ms | 0.80 / 2.42 s | 1.6% / 0.5% |
| census-style panel (2·10^6 students) | 2.45·10^6 | 9.8·10^7 | 20.5 / 17.3 ms | 2.22 / 7.60 s | 0.8% / 0.2% |

The layout takes about the same time at one thread and at eleven: its serial
passes dominate, at 5 to 8 ns per situation. At 10^9 rows they grow with the
situations, to about 1.2 s for 2.5·10^8 situations of four alternatives and
0.2 s for census-like sets (2.5·10^7 situations). The gradient grows with the
rows and the draws: at 10^9 rows it would take 100 and 10 times the times
above on this machine. On a 64-core server whose cores matched these and
scaled as they do here (6.4 to 6.5 times at eleven threads), an evaluation at
10^9 rows would take about 14 s (four alternatives) and 3.8 s (census-like) at
S = 25, of which the layout would be about 9% and 5%; at S = 100, about 3% and
1.5%. With perfect scaling the shares would be about 1.7 times as large. The
offsets take 2 GB at 2.5·10^8 situations, and the prediction kernels' base
vector 8 GB at 10^9 rows.

So at 10^9 rows on many threads and few draws, a layout cache, or a parallel
prefix sum and reductions in place of the serial passes, may pay. On a machine
with the memory for it, measure:

```sh
R_LIBS=<lib> Rscript _benchmarks/scale/layout_cost.R --units=2.5e8            # 10^9 rows; peaks near 40 GB
R_LIBS=<lib> Rscript _benchmarks/scale/layout_cost.R --units=2.5e8 --panel    # as a panel of single choices
```

It reports the layout's time against gradients at S = 25 and 100 at each
thread count given (`--threads=`, default all cores).

## Results of record

The report "Mixed Logit at Scale" (September and October 2026) holds the grid's
results for choicer 0.2.1.9000 before and after the scale audit, with each
run's load.
