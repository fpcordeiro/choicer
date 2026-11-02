# Validation suite: design and status

The derivations behind the suite's claims (plan §1: C1–C7), the metrics and gate rule, and the
specifications are written here by WP10. Until then this file holds the build's status table.

## Status

| WP | Group | Owns | State | Commit | Notes |
|---|---|---|---|---|---|
| WP0a | A | `lib/common.R`, `lib/families.R`, `lib/fit.R`, `lib/calibrate.R`, `lib/families/T0.R`, `checks/wp0_common.R` | done | (this commit) | build guard, seeds (FNV-1a + rep, method-free), atomic rep files, metadata, contracts as code; gates per method and over several cells; T0 is a plumbing fixture; 38/38 checks |
| WP0b | A | `lib/simulate.R`, `checks/wp0_simulate.R` | done | (this commit) | the suite's probability engines in choicer's parameterizations, equal to `predict()` to ≤ 7e-16 and to the fixed-parameter MXL likelihood to ≤ 3e-16 rel.; 117/117 checks; shared file, lead-serialized afterwards |
| WP1 | B | `lib/metrics.R`, `lib/sbc.R` (test half); `checks/wp1_metrics_null.R` (test-writer) | done | (this commit) | null study 85 PASS / 0 FAIL / 11 INFO (four bugs it found fixed: an MCSE copied across rows, missing na rows, SBC mirror ties, infinite W counted as a rejection); reviewer: SBC p floor → exact null by DP, route-failure rates, global-max criterion; `gate_family(rule = "plan" | "band")`, plan by default, for the user at S1 |
| WP2 | B | `mc.R`, `report.R`, `run_shards.sh`, `quick.R` (skeleton), `checks/wp2_plumbing.R` | done | (this commit) | plumbing 16/16; `run_rep()` and `vs_boot()` moved into the lib; `R_by_method`; reviewer: absent units and missing reps flagged, `--tier`, launcher checks; plumbing 20/20 |
| WP3 | C | `lib/dgp_mode.R`, `lib/families/F1.R`, `F8.R`, `F11.R`, `lib/references.R` (clogit/mlogit) | in progress | | |
| WP4 | C | `lib/dgp_sp.R`, `lib/families/F4.R` | in progress | | |
| WP5 | C | `lib/dgp_geo.R`, `lib/families/F2.R`, `F3.R`, `lib/references_quad.R` | done, S1 decisions open | | the quadrature reference gets its own file so WP3 and WP5 own disjoint files |
| S1 | — | | pending | | designs' realism, θ0, identification, fit times |
| WP6 | D | `lib/geo_gen.cpp`, `lib/families/F6.R`, `F7.R` | pending | | |
| WP7 | D | `lib/families/F5.R`, `F10.R` | pending | | |
| WP8 | D | `lib/dgp_bayes.R`, `lib/sbc.R` (loop half), `lib/families/F9.R` | pending | | |
| S2 | — | | pending | | budget from R = 2 smoke; F7 reps, server plan |
| WP9 | E | `quick.R` (assembly), `lib/fit.R` additions via the lead | pending | | |
| WP10 | E | `README.md`, `DESIGN.md`, AGENTS.md rows, article note, deletions | pending | | |
| S3 | — | | pending | | quick-tier run, then a PR |

WP0 runs as two implementers (WP0a, WP0b) on disjoint files: the probability engines in
`lib/simulate.R` are shared by every DGP of group C, so they are fixed before group C starts
rather than serialized through it.
