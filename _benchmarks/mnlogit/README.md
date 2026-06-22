# Multinomial Logit Benchmark

`run.R` compares `choicer` against reference multinomial/conditional logit
implementations on a common simulated balanced long-format dataset.

The benchmark uses two default sweeps:

- N-scaling: varies the number of choice situations at a fixed J.
- J-scaling: varies the number of alternatives at a fixed N.

Runtime output separates deterministic package setup from the estimator call:

- `setup_time_sec`: reshaping, package-specific data objects, compiled utility
  setup where relevant.
- `fit_time_sec`: model estimation only; this is the plotted headline metric.
- `total_time_sec`: setup plus fit.

Each package/spec/run attempt is executed in its own sequential child R process.
The runner writes per-attempt partial CSVs as soon as each child exits, then
reconciles those partials into the existing aggregate outputs:

- `run_status.csv`: progress ledger for pending/running/ok/error/timeout/not_installed attempts.
- `partials/`: one raw-result CSV and, when estimates exist, one coefficient CSV
  per attempt.
- `logs/`: child stdout/stderr for diagnosing package failures.

Use `--max-run-sec=<seconds>` to cap each child process. Timed-out attempts are
recorded with `status = "timeout"` and are excluded from timing summaries in the
same way as other non-`ok` attempts. The default is no timeout.

## Example

```sh
Rscript _benchmarks/mnlogit/run.R \
  --n-runs=1 \
  --n-vec=50 \
  --j-vec=3 \
  --fixed-n=50 \
  --fixed-j=3 \
  --packages=choicer,mlogit,logitr \
  --max-run-sec=60 \
  --tag=smoke
```

Packages that are not installed are recorded as `status = "not_installed"`
when `skip_missing = TRUE`.

Coefficient checks compare each package to `choicer` on the same simulated
sample after canonicalizing coefficient names to `x1`, `x2`, and `ASC_2`,
`ASC_3`, ...
