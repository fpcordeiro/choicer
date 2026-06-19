# choicer Benchmarks

This directory contains development-only benchmark scripts. It is excluded from
package builds via `.Rbuildignore`, so reference packages used here are not
listed in `DESCRIPTION`.

## Structured Benchmarks

The revamped benchmark infrastructure starts with multinomial logit:

```sh
Rscript _benchmarks/mnlogit/run.R \
  --n-runs=3 \
  --n-vec=1000,5000,10000 \
  --j-vec=5,10,25 \
  --fixed-n=10000 \
  --fixed-j=10 \
  --packages=choicer,mlogit,logitr \
  --tag=local
```

Each run writes timestamped outputs under
`_benchmarks/mnlogit/output/YYYYmmdd_HHMMSS[_tag]/`.

Default outputs:

- `metadata.json`, `metadata.txt`, `session_info.txt`
- `raw_results.csv`
- `coef_results.csv`
- `summary_results.csv`
- `summary_table.md`
- `runtime_scaling.png`, `runtime_scaling.pdf`

## Legacy Benchmarks

The older ad hoc scripts remain available through their original top-level
paths. The scripts themselves now live in model-specific folders and the
top-level files are shims.

```sh
Rscript _benchmarks/mxlogit_simulation_benchmark.R
Rscript _benchmarks/nestlogit_benchmark.R
Rscript _benchmarks/mnprobit_simulation_benchmark.R
Rscript _benchmarks/halton_benchmark.R
```
