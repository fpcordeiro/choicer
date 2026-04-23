# MXL Monte Carlo Validation Suite

Developer infrastructure for the `choicer` package. This suite exercises every
supported feature of the Mixed Logit estimator against the classical
asymptotic claims from MLE / MSL theory, at enough replications to separate
signal from MC noise.

**Not shipped with the package.** `_validation/` is listed in
`.Rbuildignore` and only exists in the development tree.

## What this validates

Seven claims operationalized as pass / fail checks on every scenario:

1. Consistency: `|bias / MC-SE| < 3` at the largest N in Scenario A.
2. Asymptotic normality: moments of `z = (theta_hat - theta_0) / se`
   (skew, excess kurt) plus four normality tests (Shapiro, AD, JB, KS).
3. Wald coverage: Wilson CI on empirical coverage contains the nominal
   level at 90 / 95 / 99 pct.
4. Information-matrix equality: `mean(SE) / sd_emp` near 1 (both Hessian
   and BHHH).
5. Simulation bias: `|bias|` approximately linear in `1 / S` (Scenario D).
6. LR statistic distribution: empirical CDF matches chi-squared via KS;
   size of nominal-5 pct test near 0.05 (Scenario A).
7. Analytical Hessian correctness: max relative error vs `numDeriv` under
   1e-4 (one rep per scenario).

## Scenario grid

| id | purpose | N | S | K_w | rc_dist | rc_corr | rc_mean | R |
|----|---------|---|---|-----|---------|---------|---------|---|
| A  | consistency across N | 1k / 2.5k / 5k / 10k | ceil(1.5 sqrt N) rounded up to 50 | 2 | normal, normal | FALSE | FALSE | 1000 (2000 at 10k) |
| B  | correlated Sigma, K_w=2 | 5000 | 100 | 2 | normal, normal | TRUE | FALSE | 1000 |
| B2 | correlated Sigma, K_w=3 | 5000 | 100 | 3 | normal x 3 | TRUE | FALSE | 1000 |
| C  | log-normal RC with mu | 5000 | 100 | 2 | log-n, normal | FALSE | TRUE | 1000 |
| D  | simulation bias O(1/S) | 5000 | 25 / 50 / 100 / 250 / 500 | 2 | normal, normal | FALSE | FALSE | 300 |
| F  | weak identification | 5000 | 100 | 2 | normal, normal | FALSE | FALSE | 1000 |

All scenarios use J = 8 inside alternatives, `outside_option = TRUE`,
`vary_choice_set = TRUE`, `use_asc = TRUE`. See
`.claude/plans/MXL_MC_VALIDATION_PLAN.md` (in the main worktree) for the
theoretical motivation and the justification of the fixed-J, `outside_option = TRUE` conventions.

## How to run

### Pre-flight (QUICK mode, < 90 minutes)

```bash
OMP_NUM_THREADS=2 QUICK=TRUE Rscript _validation/mxl_monte_carlo.R
```

QUICK mode reduces every R to 50 and trims Scenario A's N-grid to {1k, 5k}
and Scenario D's S-grid. The output artifacts are written but the study's
headline claims cannot be fully adjudicated at R = 50 alone; inspect the
per-rep RDS files for debugging only.

### Full overnight run (~75 hours on an 8-worker, 2-thread setup)

```bash
OMP_NUM_THREADS=2 Rscript _validation/mxl_monte_carlo.R
```

The driver wires `future::plan(multisession, workers = floor(cores/2))` so
parallel R-level replications combine with 2-thread OpenMP per fit without
oversubscribing the machine. Set `OMP_NUM_THREADS=1` if you plan to use
more parallel workers; set higher if your cores are hyper-threaded.

## Output artifacts

Under `_validation/output/`:

- `mc_<scenario>_raw.rds`, `mc_<scenario>_natural.rds`: long per-rep
  `choicer_mc` objects on the raw (Cholesky / log-mu) and natural
  (vech Sigma / exp mu) scales. For Scenario A / D, a named list keyed by
  N or S.
- `aux_<scenario>.rds`: auxiliary tables (LR statistics for A; first-rep
  Hessian-agreement records for all scenarios).
- `asymptotics_<scenario>_raw.csv`, `asymptotics_<scenario>_natural.csv`:
  `mc_asymptotics()` output.
- `plot_qq_<scenario>.png`: QQ of studentized z vs N(0, 1), faceted by
  parameter.
- `plot_consistency.png`: log |bias| and log RMSE vs log N (Scenario A).
- `plot_coverage.png`: empirical 95 pct coverage with Wilson bands vs N
  (Scenario A).
- `plot_sim_bias.png`: |bias| vs 1 / S (Scenario D).
- `plot_se_ratios.png`: `mean_se / sd_emp` ratios across scenarios
  (Hessian only in v1; BHHH slot reserved).
- `plot_lr_cdf.png`: empirical LR CDF vs chi-squared 1 (Scenario A).
- `REPORT.md`: per-scenario pass / fail matrix and one-paragraph narratives.

## What this suite does NOT validate

**Cross-package validation**: comparing `choicer` point estimates and SEs
against `mlogit`, `gmnl`, and Apollo on a canonical dataset is a
complementary artifact and a high-leverage sanity check, but it is
orthogonal to Monte Carlo parameter recovery. Add a separate
`_validation/cross_package.R` when desired.

**Incidental ASCs (J -> infinity) asymptotics**: Scenario A fixes J = 8
across the N-grid on purpose — standard MSL theory is a fixed-J, N -> inf
regime. A companion scenario that scales J with N would turn each ASC
`delta_j` into an incidental parameter (Neyman-Scott) whose information
scales with per-alternative frequency rather than total N, so the
sqrt(N)-rate claim would not apply. That study is out of scope here.

**No-outside-option identification**: every scenario uses
`outside_option = TRUE`. A separate no-OO scenario would exercise the
first-ASC-normalization branch of `run_mxlogit()` but is deliberately
excluded here to keep dimensionality down.
