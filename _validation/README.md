# Validation suite

Developer infrastructure for `choicer`, not shipped with the package (`_validation/` is in
`.Rbuildignore`). The suite asks two things of the frequentist and Bayesian estimators: do they
estimate what theory says they estimate, and is the uncertainty they report honest? It asks on
designs that resemble the package's target applications: stated-preference panels,
revealed-preference mode choice, hospital choice from claims, school choice from a census.
`DESIGN.md` holds the derivations behind the claims, the metrics' formulas, the designs' calibration
and what the build measured.

## What it can and cannot say

A pass says that, on these data-generating processes and at this number of replications, an
estimator's bias, the calibration of its standard errors and the size of its tests do not deviate
from what theory predicts by more than a tolerance fixed in advance. Agreement with truth on
realistic designs is strong evidence that the implementation is right and that its inference can be
trusted there.

It cannot say the following.

- That a model suits an application. Except in the misspecification designs (F8's unweighted fits
  on the NL population, F11), the fitted model is the true one. The suite validates the estimators, not any modeling
  choice.
- That a result extends beyond its design. J, N, S and the share of the outside option are those of
  the tables below.
- What a user reaches from choicer's default start. Every mixed logit fit runs from the default
  start and from the truth and keeps the higher certified likelihood (S1 in `DESIGN.md`); F7 takes
  Newton steps from the truth (S2). So the verdicts concern the maximum of the likelihood, not the
  optimizer's path. How often the default start lost is recorded (`diag$start`), and it is a finding
  about the package.
- That simulation error is inside the standard errors. The variance routes estimate the sampling
  variance given the draws; family F5 measures the simulation error beside them.
- Anything where a claim is measured, not proved. The incidental-constants prediction for the MXL
  and the NL (C6) is derived by analogy to panel models; no published theorem covers it, so it is a
  finding, never a gate.
- That derivatives are exact, or that malformed inputs are handled. Those checks live in
  `tests/testthat/` and are not repeated here.

**The long tier has not been run yet.** Numbers quoted in this directory come from the build's
probes and from the quick tier's calibration, each with the R it was measured at. A server run comes
later.

## The three gates

| Gate | Question | Compared with | Cost |
|---|---|---|---|
| `golden/` identity harness | Did this change move any number between two builds? | another build, bit for bit | a minute or two per build |
| `quick.R` | Are the estimates and the inference still right? | the truth, exact references, other packages | estimated 9.5 minutes on 8 workers (budget 15) |
| `mc.R` | Do the asymptotic claims hold on realistic designs? | the truth, over replications | about 790 CPU-hours |

The golden harness is not the quick run. It catches regressions but never compares anything with the
truth, so it cannot say whether either build is right. The other two can, and neither can prove that
a change moved nothing.

## Running it

**The build under test.** Install it from a copy of its tree into a library of its own, and name
that library first on `R_LIBS`:

```sh
mkdir -p /tmp/lib-dev && R CMD INSTALL --library=/tmp/lib-dev /tmp/copy-of-the-tree
```

Every runner stops unless `R_LIBS` is set, its first entry holds choicer, the namespace was loaded
from there and the version is at least 0.2.1.9000. The site library's 0.2.1 lacks `person_col`,
`polish` and `concentrate_out_asc`, and it would otherwise shadow the development build. Commands
below write `<lib>` for that library.

**Dependencies.** Besides choicer's: `survival`, `mlogit` (with `dfidx`) and `bayesm` for the quick
tier's references, `data.table`, `Rcpp` and a C++ compiler for the generator, `numDeriv` in the
self-checks.

**The generator's build cache.** Families F5, F6 and F7, and Q0 and Q2 of the quick tier, simulate
through `lib/geo_gen.cpp`, compiled once by `Rcpp::sourceCpp()`. `VS_GEO_CACHE` names its cache
(default `_validation/output/geo_gen_cache`). Processes sharing a cache compile once, under a lock;
give a machine or a checkout one cache.

**Outputs** go to `_validation/output/` (ignored by git) or to `--out=DIR`. Inside a package tree
the suite accepts only a directory git ignores.

### The quick tier

```sh
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/quick.R
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/quick.R --only=Q0,Q1.3,Q2:F1,Q3:F2 --workers=4
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/quick.R --only=Q3:F1,Q4 --reps=2
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/report.R --run=<out>/quick_<stamp>_reps --tier=quick
```

`--only` selects tests, and `Q2:<entry>` and `Q3:<entry>` select one of their entries. `--workers`
(default 8) sets the PSOCK workers, each at one thread. `--reps=N` runs every cell of Q3 and Q4 at N
reps instead of 100: that exercises the plumbing, since the gates hold at 100. `--perturb=x` adds x
to the truth on choicer's side of Q0, a failure injection for `checks/wp9_quick.R`. `--fixtures`
adds the plumbing fixture T0 to Q3. `--ref=NAME` is recorded in the metadata, and `--out=DIR`
replaces the default `_validation/output/quick/`.

It writes `quick_<stamp>.csv` (every check row, with the test's seconds), `quick_<stamp>_info.csv`
(rows that carry information and no verdict), `quick_<stamp>_Q3_<entry>.csv` and
`quick_<stamp>_Q4.csv` (the metrics), and the reps of Q3 and Q4 under `quick_<stamp>_reps/`. It
prints a summary, the gate tables and the report's command. It exits 1 if any check fails. A missing
reference, a family that does not load, a test that returns no row and a gated row without a value
all count as failures.

| Test | What it compares | Claim | Pass rule |
|---|---|---|---|
| Q0 | Every DGP against choicer at the truth: the probabilities its choices are drawn from against `predict()` (MNL, NL), and the suite's simulated log-likelihood on a shared cube of draws, and on the DGP's own tastes as a one-draw cube, against a fit with every parameter fixed (MXL) | the conventions | max abs diff below 1e-12; relative difference below 1e-12 |
| Q1.1 | MNL on `mode_choice` against `survival::clogit` and `mlogit` | same estimator | \|Δθ\|/SE below 1e-4, relative SE difference below 1e-4, \|ΔlogL\| below 1e-6 |
| Q1.2 | NL (ground nest) against `mlogit` with `un.nest.el = TRUE` | same estimator | \|Δθ\|/SE below 1e-3, \|ΔlogL\| below 1e-5 |
| Q1.3 | MXL at S = 4,000 against the suite's exact Gauss-Hermite maximum likelihood (D4-40, N = 2,000) | C4 | \|Δθ\|/SE below 0.05 |
| Q1.4 | MNP against `bayesm::rmnpGibbs`, same prior, base alternative aligned | same posterior | means within 4 combined MCSEs, SDs within 10% |
| Q2 | One large fit per family: F1 (N = 10^6), F2 (200,000), F3 (16k), F4 (3.2k persons), F6 MXL (20k), F6 MNL (80k) | C1 | polish certified; max \|z\| below 4; BHHH and robust SEs within 10% of the Hessian's (for a constant with n choices, within the larger of 10% and 2/√n, 4 SD of the routes' difference) |
| Q3 | Mini Monte Carlos, R = 100: the bottom rungs of F1, F2, F3 (S = 200), F4 (S = 200), and F8's unweighted MNL under quotas | C1, C2 | the gate rule at the quick tolerances, on each family's `gates_quick`, F8's including the exact constant shift; at most 2% failed reps per cell and method; information rows say how often the default start failed or lost |
| Q4 | Mini SBC, R = 100: MNP, HMNL, HMNP on small data and short chains | C7 | ECDF test, Holm p above 0.001; at most 2% failed reps; each quantity's ranks cover at least 95% of the reps |

Each family's `gates_quick` holds only the checks theory supports at its bottom rung, confirmed at R
= 300 (R = 3,000 for the cheap families) on the family's own seeds; a comment beside it lists what
was dropped and why, with the numbers.

### The long tier

```sh
# One family, or part of one, in one process:
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/mc.R --family=F1 --run=try
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/mc.R --family=F3 --cells=16k --methods=ml --reps=1:100 --run=try
# Shard 2 of 8 of a run, by hand:
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/mc.R --family=F2 --shard=2/8 --run=f2
# Background shards, after a warm pass that compiles the generator once:
R_LIBS=<lib> VS_GEO_CACHE=<dir> _validation/run_shards.sh --family=F1,F2,F8,F11 --shards=32 --run=cheap
R_LIBS=<lib> VS_GEO_CACHE=<dir> _validation/run_shards.sh --family=F7 --shards=4 --threads=16 --run=f7
# More replications of F7's mixed logit, into the same run:
R_LIBS=<lib> VS_GEO_CACHE=<dir> _validation/run_shards.sh --family=F7 --cells=claims,census --methods=mxl --reps=101:1000 --shards=4 --threads=16 --run=f7
# Failed reps again; seed ranges checked:
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/mc.R --family=F3 --run=f3 --retry-failed
R_LIBS=<lib> Rscript _validation/mc.R --check-seeds
# The report:
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/report.R --run=f3
```

`--family` takes an id, a list or `all`. A unit (family, cell, method, rep) writes one file,
`<out>/<run>/<family>/<cell>/<method>/rep_NNNNN.rds`, atomically. Each rep is simulated once and its
data go to every method of the rep, so comparisons between methods are paired. Running again with
the same `--run` skips the units with a readable file, so a killed or preempted run resumes where it
stopped and extending `--reps` adds to it. Shard k of K takes the reps with `rep %% K == k - 1`. A
failed rep is a record (error, warnings, optimizer code) that counts in the failure rate;
`--retry-failed` reruns those. The seed is a function of family, cell and rep and not of the method,
and a rep is reproducible at a fixed thread count, bit for bit at one thread.

`run_shards.sh` passes other arguments to `mc.R`, starts each shard with `nohup` from a generated
script under `<run>/logs/`, and prints the commands to watch, stop, resume and report the run. It
starts nothing while a shard of the same run is alive. Its `--threads` sets each shard's OpenMP
threads.

A sensible order is the cheap families with hard gates first and F7 last.

### Calibration and self-checks

`lib/calibrate.R` prints a DGP's moments against their targets, and its truth:

```sh
R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/lib/calibrate.R --dgp=D4-40 --variant=nl --N=200000
```

The suite's own checks exit 1 on a failure, and each file's header lists its sections:

```sh
R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 Rscript _validation/checks/wp0_common.R
```

| File | What it checks |
|---|---|
| `checks/wp0_common.R` | the build guard, arguments, seeds, rep files, the output guard, the family contract (under 2 minutes) |
| `checks/wp0_simulate.R` | the suite's probability engines against choicer, and choice frequencies over 10^6 situations |
| `checks/wp1_metrics_null.R` | the metrics on fake records with known distributions: MCSE calibration, gate size and power, the exact SBC test |
| `checks/wp2_plumbing.R` | the runners and report on a fixture: shards, resume after a kill, failures, `--retry-failed` (under 5 minutes) |
| `checks/wp3_mode.R` | D2, F1, F8, F11, Q1.1, Q1.2, the pseudo-true value and the exact constant shift |
| `checks/wp4_sp.R` | D1, D1-L, F4, the start-policy study |
| `checks/wp5_geo.R` | D4's designs, F2, F3, Q1.3, identification strength |
| `checks/wp6_geo.R` | the C++ generator against its R reference, D3, D3-claims, D4-census, F6, F7 |
| `checks/wp7_sim.R` | F5 and its rate uncertainty on synthetic records, F10, and F10 against the script it replaced (read from git history, commit ae38e12) |
| `checks/wp8_bayes.R` | the priors and data of the Bayesian samplers, the SBC loop, F9, Q1.4, Q4 |
| `checks/wp9_quick.R` | `quick.R` through its command line, with injected failures; the report's verdicts against the runner's on fixture records |

## The claims

`DESIGN.md` states each claim's conditions and derivation, with references.

| Claim | What theory predicts | Tested in |
|---|---|---|
| C1 First order | Under correct specification, √N(θ̂ − θ0) is asymptotically N(0, H⁻¹) and H = B, so Wald intervals and the joint Wald test reach nominal level, LR is χ²_q for interior restrictions (λ = 1 is interior in choicer's NL), and the delta method holds | F1, F2, F3, F4 panel, F7's MNL; Q1, Q2, Q3 |
| C2 Sampling designs | Choice-based samples: unweighted MNL slopes are consistent and the constants shift exactly by log(H_j/Q_j) − log(H_ref/Q_ref); WESML is consistent, its plain robust variance exact under random stratum sizes and conservative under quotas. A panel fitted as a cross-section is consistent and inefficient, with cluster SEs valid and Hessian SEs too small | F8; F4 composite |
| C3 Misspecification | θ̂ converges to the pseudo-true θ*, and only the sandwich estimates its variance | F11 |
| C4 Simulation | MSL's bias is about b/S, so bias/SD is O(√N/S), and MSL is equivalent to ML only if √N/S goes to 0. choicer's permuted Halton has no randomized-QMC guarantee, so the error is measured. Panel integrands raise it | F5; the gates' condition √N/S ≤ 0.25; Q1.3 |
| C5 Higher order | N·bias is stable; one-sided tail errors are O(N^-½) and two-sided coverage errors O(N⁻¹); E[LR] = q(1 + c/N); ratios such as WTP are skewed when the denominator is weakly determined | ladders in F1 to F4 (findings) |
| C6 Many alternatives | With J constants the MNL's slopes carry no leading-order incidental bias, but its non-linear functionals carry O(J/N) and its sandwich is biased down. For the MXL and NL the prediction is bias/SD of O(J/√N), derived and measured, not a theorem | F6, F7 |
| C7 Bayesian samplers | SBC ranks are uniform at any N. Bernstein-von Mises gives coverage and posterior SD over ESE of 1 for quantities whose information grows with N | F9; Q1.4, Q4 |

## Metrics and the gate rule

Every metric comes with a Monte Carlo SE (Morris, White & Crowther 2019). A failed rep is recorded
and counted in the failure rate, while the other metrics read the reps that survive, so read the two
together. Formulas as coded are in `DESIGN.md`.

| Metric | Definition | Target | mc tolerance | quick tolerance |
|---|---|---|---|---|
| `bias_std` | (mean θ̂ − θ0) / ESE | 0 | 0.10 | 0.30 |
| `ese` | sd(θ̂) | stable as √N·ESE | finding | finding |
| `se_ratio` | √mean(ŝe²) / ESE, per variance route | 1 | 0.05 | 0.15 |
| `tail_lower`, `tail_upper` | P(θ0 below the 95% interval), P(θ0 above it), per route | 0.025 each | 0.75 pp | 2.5 pp |
| `coverage` | 1 minus the tails | 0.95 | 1.5 pp | 5 pp |
| `joint_wald_size`, `lr_size` | rate of p below 0.05 | 0.05 | 1.5 pp | 2.5 pp |
| `lr_mean_q` | mean(LR)/q | 1 | finding | finding |
| `route_failure_rate`, `joint_wald_failure_rate` | reps where a variance route or the Wald statistic failed | 0 | 0.5 pp | 2 pp |
| `failure_rate` | failed reps | 0 in the quick tier | finding | 2 pp per cell and method |
| `global_max_rate` | multi-start refits that found a higher maximum at another point | | finding | finding |
| `sbc_ecdf` | the exact ECDF test of SBC ranks | uniform | Holm p below 0.01 alone | Holm p below 0.001 (Q4), ranks in at least 95% of the reps |
| `pooled_*` | coverage pooled over the constants within each rep, by expected-count tercile | 0.95 | finding | finding |

**The gate rule.** A gated row fails when its Holm-adjusted p-value, within the family, is below
0.01 and its deviation from target exceeds the tolerance. The p-value is the exact binomial test for
proportions, the t test for `bias_std`, and the normal reference on deviation over MCSE otherwise.
Gates apply only where the theory's conditions hold at the design: the top rung with √N/S ≤ 0.25,
SBC, the exact constant shift, cluster SEs in F4, robust SEs in F11. Everything else is a finding:
tabulated with its MCSE, to be compared with the predicted sign or rate, never a verdict.

**The quick tier's verdicts.** `quick.R` and `report.R --tier=quick` share one definition of them,
`gate_quick()` in `lib/metrics.R`, so the report's gated rows are the runner's checks. It is the
gate rule at the family's `gates_quick$alpha` (0.01 unless the family sets one; F9's is 0.001, as
Q4 asks) and two more rules, because the metrics read only the reps that survive:

- a cell and method fails when more than 2% of its reps failed: its `failure_rate` row is gated at
  target 0, outside the Holm adjustment (the T0 fixture, whose `fail` method fails on purpose, is
  exempt);
- an SBC quantity fails when its ranks cover under 95% of the cell and method's reps, whatever its
  p-value says.

A findings-only family has no gated rows in either tier. `checks/wp9_quick.R` compares the two sides
row by row on fixture records. The long tier's verdicts are the gate rule alone.

## The designs

The DGPs never call choicer's kernels. Their probabilities come from `lib/simulate.R`, in choicer's
parameterization (Q0 ties the two), and θ0 is frozen in the registry. No DGP writes physical outside
rows. Alternatives never chosen in a replication are dropped, with the situations offering one
chosen wherever offered: this is ML's limit, and the counts are recorded.

| DGP | Real-data analogue | Structure | Used by |
|---|---|---|---|
| D1 | Stated-preference conjoint (Revelt & Train 1998) | N persons × 8 tasks, 3 unlabeled alternatives and an opt-out, level-balanced. Fixed price, an irrelevant attribute, a shifted log-normal on −time and a normal on quality, correlated | F4, F5 |
| D1-L | Labeled products | D1's variant `labeled`: J = 6 products with one characteristic, an implicit outside option, normal tastes; logit or probit shocks | F9 |
| D2 | Intercity mode choice (`mode_choice`) | The 210 travellers' covariates resampled to Greene's population shares, J = 4. Variants: MNL, NL (ground nest), MXL (normal travel coefficient with σ = \|μ\|) | F1, F8, F11, F9's MNP |
| D3 | Hospital choice from claims | 200 hospitals over 200 × 200 km in 10 towns, the hospitals within 30 km (3 to 30, mean 21.5), no outside option, the largest hospital the reference, a cross-section | F6, F5 |
| D3-claims | The same at population scale | 2,000 hospitals, sets of 3 to 60 (mean 39), a panel: 625,000 patients, about 10^6 admissions, 3.9 × 10^7 rows | F7, F5 |
| D4 | School choice from a census | "D4-small": 4 markets of 40 schools, 160 constants. No family uses it now; it is the design of S1's measurement of C6 | checks only |
| D4-40 | The same, few constants | One market of 40 schools, each student's 15 nearest and an outside option; sectors public, charter, private; NL nests by sector | F2, F3, F5, Q1.3 |
| D4-census | The same at population scale | 150 markets of 20 schools (J = 3,000), 10^6 students, 1.5 × 10^7 rows | F7, F5 |
| DB | The samplers' data | MNP on D2's design (base `car`), HMNL and HMNP on D1-L, truth from the prior for SBC | F9, Q1.4, Q4 |

F10 keeps the DGP of the Monte Carlo behind §5.9 of the mixed logit article, inline in
`lib/families/F10.R`.

## The families

R is the long tier's replications per cell. Cells ending in `.q` or starting with `q_` belong to the
quick tier.

| Family | Claims | DGP | Cells | R | Methods | Gated |
|---|---|---|---|---|---|---|
| F1 MNL | C1, C5 | D2 | N = 500, 2k, 8k | 5,000 | `ml` | at 8k: bias, Hessian SE ratio and tails of every coefficient, joint Wald on the slopes, two-sided VTTS coverage, LR test of IIA |
| F2 NL | C1, C5 | D4-40 NL | N = 1k, 4k, 16k | 5,000 | `ml` | at 16k: bias, Hessian SE ratio and tails of β, λ and the constants with expected count at least 50; joint Wald on β and λ; WTP tails |
| F3 MXL | C1, C4, C5 | D4-40 | N = 1k, 4k, 16k; S = 500, 500, 520 | 1,000; `blp` 500 at 16k | `ml`; `blp` at 16k | at 16k, `ml`: bias, Hessian SE ratio and tails of β and L, joint Wald, Σ's tails |
| F4 MXL panel | C1, C2, C4 | D1 | persons = 200, 800, 3.2k; panel S = 500; composite S = 800 at 200 and 800 | 1,000 | `panel`; `xsec` (composite, cluster) | `panel` at 3.2k: bias, Hessian SE ratio and tails of every parameter, joint Wald, LR exclusion, Σ's tails; `xsec` at 800: bias, cluster SE ratio and tails |
| F5 simulation | C4 | regimes on F3 16k, F4 3.2k, F6 20k, F7 claims, F7 census | S = 50 to 800 × 2 seeds, store at 100, 200, against S = 3,200; population regimes S = 50 to 400 against 800 | 30 datasets; population 1 | `msl` | findings only |
| F6 many alternatives | C6 | D3 | N = 5k, 20k, 80k, J = 200 (N/J = 25, 100, 400); S = 200; MNL data in cells `.mnl` | 1,000, 1,000, 500; `blp` 500 at 20k | `mnl`, `mxl`, `blp` at 20k | `mnl` at every rung: slopes' bias and Hessian SE ratio |
| F7 population | C1, C6 | D3-claims, D4-census | `claims.mnl`, `claims`, `census.mnl`, `census`; S = 100 | MNL 1,000; MXL 100 (extendable); `blp` 25 | `mnl`, `mxl`, `blp` on the census | `mnl` on both: slopes' bias, Hessian SE ratio and tails, joint Wald |
| F8 choice-based | C2 | D2 | `quota`, `random`, `random_nl`; N = 2,000 | 2,000 | `mnl_unweighted`, `mnl_wesml`, `nl_wesml`, `nl_unweighted`; `mnl_misspec` on `random_nl` | at `random` and `random_nl`: unweighted MNL's bias against the exact shift, Hessian SE ratio and tails, joint Wald; WESML robust tails |
| F9 Bayesian | C7 | DB | SBC: MNP N = 300, HMNL and HMNP 100 persons × 8. BvM: MNP N = 500, 2k, HMNL 200, 800 persons | SBC 1,000; BvM 500 | `mnp`, `hmnl`, `hmnp` | SBC cells: the ECDF test on every scalar and the log-likelihood; BvM findings |
| F10 share matching | C1 | inline | `cross`, `panel`, `base`, `markets`; N = 2,000, J = 10, S = 100 | 1,000 | `ml`, `blp_sample`, `blp_exact`, `blp_n10`, `blp_n1` | findings only |
| F11 misspecification | C3 | D2 MXL fitted by MNL | N = 500, 8k | 2,000 | `ml` | at 8k: bias against θ*, robust SE ratio and tails, robust joint Wald |

F5's rates of simulation error take their SEs from complete units, the spread of each dataset's own
rate over the datasets (over the seeds in the population regimes, which have one dataset). A
ladder's S points share datasets, reference fits and draw seeds, and the S-level MCSEs miss that
covariance.

## Reading a report

`report.R --run=<dir or tag> [--tier=mc|quick]` writes into the run directory:

- `REPORT.md`: a summary per family (gated rows passing and failing, reps missing), then for each
  family its claims and gates, its units against the design (absent units and missing reps flagged),
  its failing rows, a section per cell and method (gated rows first, findings below, blocks of more
  than ten constants summarized and pooled, failures and warnings by class), its own findings tables
  where it has them (F5, F10), and plots; last the run's metadata (choicer's path and version,
  `--ref`, git SHA and uncommitted paths, threads, BLAS, CPU);
- `tables/`: `<family>_metrics.csv` with every row, `<family>_units.csv`, `summary.csv`, `meta.csv`,
  `<family>_findings_<name>.csv`, and failures, warnings and messages;
- `plots/`: QQ plots of the Hessian z at the gate cells, ladders against each varying cell column
  (ESE·√N, N·bias, tail imbalance, with fitted slopes), SBC ECDF differences with simultaneous
  bands.

`--tier=quick` applies the quick tolerances to the units quick.R ran, with quick.R's own verdicts
(`gate_quick()`, under "The quick tier's verdicts" above): its Holm level, the bound on failed reps
and the rank coverage, row for row.
The report's exit status says whether it rendered (1 on a section's error or a run without rep
files), not whether a gate passed: read the summary.

A row's `verdict` is one of:

- **`pass`**: gated, and either the deviation is within tolerance or it is not significant after
  Holm.
- **`fail`**: gated, Holm p below 0.01, and the deviation beyond tolerance (SBC: Holm p alone; in
  the quick tier, `gate_quick()`'s rules above). Look at the estimand, its MCSE and whether the
  theory's conditions hold before blaming the estimator.
- **`na`**: gated, but the value, its MCSE or its tolerance is missing: a gate entry that matched no
  row, or every rep failed. It is never a pass, and the quick tier fails on it.
- **`finding`**: not gated. Read it against the predicted sign or rate.

Gating an SE metric also gates that route's failure rate, and gating a joint Wald size its failure
rate, so that reps dropped by a failed route count against the gate.

The **`band`** column, on gated rows only, is another rule's verdict, shown for comparison: it tests
whether the true value lies inside the tolerance band, so a true value anywhere in the band fails a
family with probability at most 0.01, at any R. The plan's rule tests the value against its target
and then asks the estimate to exceed the tolerance. `fail` with `band` `pass` therefore says that
the estimate is significantly off target and past the tolerance, but not significantly past it. S1
kept the plan's rule for the verdict.

## Budget and server

CPU-hours at one thread, estimated from measured per-fit times, not from a run. F7's figure is an
extrapolation, linear in rows and draws: no full-size fit at S = 100 has been run, so the first
long-tier reps' `diag$newton_steps` replace the extrapolated step counts.

| F1, F2, F8, F11 | F3 | F4 | F5 | F6 | F7 | F9 | F10 | Total |
|---|---|---|---|---|---|---|---|---|
| about 4 | about 78 | about 59 | about 35 | about 234 | about 340 | about 28 | about 12 | about 790 |

The long tier needs one machine of about 64 cores and at least 64 GB of RAM, so as to hold four
concurrent population fits of at most 15 GB at 16 threads. It will be set up in a later session;
`run_shards.sh` and the rep files are what it relies on (resume after preemption, one process per
shard).

## Excluded by design

The suite tests well-identified designs and leaves out three things on purpose:

- weak identification near σ = 0, a boundary: log σ goes to −∞ in choicer's parameterization, and
  the LR test of σ = 0 is a mixture of chi-squares (Self & Liang 1987);
- Fieller and Krinsky-Robb alternatives to the delta method;
- cross-package MXL comparisons under other draw schemes, which measure the draws. Q1.3's exact
  quadrature takes their place.

The MXL-only Monte Carlo this suite replaced, its helpers, and the one-off script F10 ports are in
git history (commit ae38e12 holds all three).

## Layout

| Path | Purpose |
|---|---|
| `quick.R`, `mc.R`, `run_shards.sh`, `report.R` | the runners and the report |
| `lib/common.R`, `lib/families.R`, `lib/fit.R` | the build guard, seeds and rep files; the DGP registry and family contract; the rep record, variance routes, two starts, multi-start |
| `lib/simulate.R`, `lib/dgp_mode.R`, `lib/dgp_sp.R`, `lib/dgp_geo.R`, `lib/geo_gen.cpp`, `lib/dgp_bayes.R` | the probability engines and the DGPs, each header giving its design and calibration |
| `lib/calibrate.R` | a DGP's moments against its targets |
| `lib/metrics.R`, `lib/sbc.R` | the metrics, the gate rule, the SBC test and loop |
| `lib/references.R`, `lib/references_quad.R`, `lib/quick_checks.R` | Q1.1 and Q1.2, Q1.3, Q0 and Q2 (Q1.4 and Q4's check are in `lib/dgp_bayes.R`) |
| `lib/families/` | one file per family, its header stating its design, methods and gates; `T0.R` is a plumbing fixture |
| `checks/` | the suite's self-checks |
| `golden/` | the identity harness |
| `DESIGN.md` | derivations, metrics' formulas, calibration, measurements, status, decisions, package findings |

**Cholesky packing.** All mixed logit functions pack `L` row-major lower-triangular, `[L_11, L_21,
L_22, ...]`, with the diagonal on the log scale. The suite's `sigma_vech()` and `L_from_theta()`
follow it.

## Golden identity harness (`golden/`)

`golden/refactor_identity_check.R` records the outputs of one build of choicer
and checks another build's against them, leaf by leaf: every compiled entry
point of the multinomial, mixed and nested logit kernels over grids of flags,
draw modes, weights, panels, edge inputs and error cases, plus the data
preparations, the Bayesian samplers and fits through the `run_*()` wrappers
(about 35,000 leaves). Leaves must be bit-identical, except a tolerance group
(the mixed logit estimation kernels' outputs whose summation order a change may
alter), which passes when max|new - base| / max(1, max|base|) is at most 1e-10.
It is the identity gate of the scale audit: a change that should not move
results passes `check` against a baseline of the build before it.

```sh
# Each build installed from a copy of its tree into its own library:
R_LIBS=/tmp/lib-main Rscript _validation/golden/refactor_identity_check.R baseline --ref=main
R_LIBS=/tmp/lib-dev  Rscript _validation/golden/refactor_identity_check.R check --ref=main
# Everything bit-identical, the tolerance group included:
R_LIBS=/tmp/lib-dev  Rscript _validation/golden/refactor_identity_check.R check --ref=main --tol-slots=none
```

The build under test is the first entry of `R_LIBS`; the harness stops if
choicer was loaded from anywhere else, and runs at one thread. Baselines,
per-leaf reports (`refactor_check_report_<ref>.csv`) and call logs go to
`_validation/output/golden/` (ignored by git), or to `--out=<dir>`; the harness
refuses an output directory inside a package tree that git does not ignore.
The batteries are the `*_batteries.R` files beside it (the header lists them);
a new one defines `run_<name>()` and is added to `results` there. `check` exits
with status 0 on a pass and 1 on a failure, and prints a summary per group, the
failing leaves and the worst tolerance-group deviations.

The share inversion (`blp()` and the three exported contractions) is checked by
a solver tier instead of identity, because its results may move within its
tolerance (BLP engine plan, user decision U1). Every BLP call of the batteries
goes through `golden_blp()` (`golden/blp_batteries.R`, sourced first), which
records one leaf per call: the value, warnings, messages, error and printed
output; the residual at the returned delta recomputed through the untouched
`*_predict_shares()` kernel (the largest absolute log-share residual,
re-anchored on each component's reference without an outside option, and each
component's outside share's with one); a reference solution with its residual
floor and conditioning kappa (the delta the targets were computed at, where the
call site knows it, or else Newton's on the same kernel); and the engine's
input checks that main did not make, read off the call's arguments by the
harness itself (`gb_expect_error()`). A leaf passes if it is identical to the
baseline's, or if: it converged within `tol` with every constant within
2 kappa (tol + floor) of the reference and each component's reference held at
its start; or it did not converge, warns with class `choicer_blp_unconverged`
and reports the residual it stopped at after `max_iter` evaluations, where the
baseline did not meet its `tol` either (or `max_iter` binds by design); or it
stops with the baseline's error; or it stops exactly where its inputs fail one
of those checks, with that check's message (or one its call site expects).
`check` prints the tier's report: counts, the worst residual and bound ratios,
and the entries less accurate than the baseline (reported, not gated: the
baseline's contraction returned one unverified step past its last check).
`golden/blp_batteries.R`'s own battery, `blp`, makes the inversion iterate:
counterfactual, symmetrized, uniform and tiny targets, regional choice sets,
two components with and without an outside option, repeated codes, nested
logits at lambda 0.3 to 3 with damping, mixed logits with stored, generated and
on-the-fly draws, `blp()` on fits, `max_iter` 0 and 3, and the inputs the
engine's checks stop.
