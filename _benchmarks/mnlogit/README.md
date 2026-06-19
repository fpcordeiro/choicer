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

## Example

```sh
Rscript _benchmarks/mnlogit/run.R \
  --n-runs=1 \
  --n-vec=50 \
  --j-vec=3 \
  --fixed-n=50 \
  --fixed-j=3 \
  --packages=choicer,mlogit,logitr \
  --tag=smoke
```

Packages that are not installed are recorded as `status = "not_installed"`
when `skip_missing = TRUE`.

Coefficient checks compare each package to `choicer` on the same simulated
sample after canonicalizing coefficient names to `x1`, `x2`, and `ASC_2`,
`ASC_3`, ...
