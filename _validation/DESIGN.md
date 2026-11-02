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
| WP3 | C | `lib/dgp_mode.R`, `lib/families/F1.R`, `F8.R`, `F11.R`, `lib/references.R`, `checks/wp3_mode.R` | done | (this commit) | 25/25 checks; Q1.1 7e-10 SE, Q1.2 6e-8 SE; θ* exact (210 traveller types), checked against an N = 10^7 fit; `require_converged()` |
| WP4 | C | `lib/dgp_sp.R`, `lib/families/F4.R`, `checks/wp4_sp.R` | done | (this commit) | 29/29 checks; D1 and D1-L (HB, logit and probit); `sigma_vech()`; start-policy study |
| WP5 | C | `lib/dgp_geo.R`, `lib/families/F2.R`, `F3.R`, `lib/references_quad.R`, `checks/wp5_geo.R` | done | (this commit) | the quadrature reference gets its own file so WP3 and WP5 own disjoint files |
| S1 | — | | decided 2026-11-02 | | see "Decisions" below; applied in the next commits |
| WP6 | D | `lib/geo_gen.cpp`, `lib/families/F6.R`, `F7.R`, `lib/dgp_geo.R` (from WP5), `checks/wp6_geo.R` | done | (this commit) | C++ port equals the R reference (39/39 and 300 lattice cases); 187/187 checks; contraction self-test and compile lock; F7 ≈ 1,060 CPU-h vs plan 650 (S2) |
| WP7 | D | `lib/families/F5.R`, `F10.R`, `checks/wp7_sim.R` | done | (this commit) | F10 reproduces `blp_concentration_mc.R` rep by rep (store cube kept); F5 warm-started Newton ladders on F3, F4, F6, F7; `findings` tables; two-start tie rule (tol 1e-6); `"blp"` fits certified by their decrement; F4's gate S = 500 has simulation bias < 0.1 SE (12-dataset probe) |
| WP8 | D | `lib/dgp_bayes.R`, `lib/sbc.R` (loop half), `lib/families/F9.R`, `checks/wp8_bayes.R` | done | (this commit) | 55/55 checks; Q1.4 vs bayesm 2.2 MCSE / SD 2.8%; Q4 passes; chains 5k + 30k thinned by 3 (10k left ESS below 99); MNP base `car`; BvM tail ESS recorded, chains lengthened (F9 ≈ 28 CPU-h) |
| S2 | — | | in progress | | R = 2 smoke of every family: 142/142 ok, report renders, 0 gated failures |
| WP9 | E | `quick.R` (assembly), `lib/fit.R` additions via the lead | pending | | |
| WP10 | E | `README.md`, `DESIGN.md`, AGENTS.md rows, article note, deletions | pending | | |
| S3 | — | | pending | | quick-tier run, then a PR |

WP0 runs as two implementers (WP0a, WP0b) on disjoint files: the probability engines in
`lib/simulate.R` are shared by every DGP of group C, so they are fixed before group C starts
rather than serialized through it.

## Decisions

**S1 (2026-11-02), the user's answers:**
- Every mixed logit fit starts twice, from choicer's default and from θ0, and keeps the higher
  likelihood; how often the default lost is recorded (a package finding). From θ0 no fit stopped
  on a plateau; from the default, F4's composite fit at 200 persons ended below the θ0 fit in 17
  of 30 reps, and both 3.2k composite fits run were 13 log-likelihood units below.
- F2 and F3 move to a 40-school design: with 160 constants the distance coefficient's bias/SD was
  −0.74 at 1k and −0.41 at 4k, scaling as J/√N (C6) and vanishing with the constants known, so a
  gate at 16k (predicted −0.20) would test C6, not C1. The many-constants regime stays in F6/F7.
- The plan's gate rule stays (Holm p < 0.01 and the estimate beyond tolerance); the band rule's
  verdicts are reported beside it.
- The quick tier gates, at its small rungs, only checks theory supports there (`gates_quick` per
  family), confirmed once at R = 300.
- Applied by the lead's recommendation without objection: F4's composite fit is gated at the
  800-person rung with S = 800 (√U/S = 0.1); F3's top-rung S = 520 (√N/S ≤ 0.25); F1 gates the value
  of travel time's two-sided coverage (C5 skews its tails); D2-MXL's σ = |μ|; F8's random cell by
  Bernoulli sampling, plus an NL-population cell; F2 gates constants with expected count ≥ 50;
  D4's `rc_mean = FALSE` with a free-location shifted log-normal on distance.

## Package findings (for one follow-up session after this branch)

1. nloptr's NLOPT_FAILURE (code −1) returns the start values without a warning, in the MNL and NL
   (and once at the optimum in 150 MNL fits).
2. `run_nestlogit()`'s default `scale_vars = "none"` failed at the start in about half of D2's WESML
   NL fits (code −1); `"bhhh"` failed in 0 of 120.
3. The mixed logit's Newton polish certifies points on a zero-variance plateau (an L diagonal near
   −19, SE ≈ 10^4), and the default start reaches them often in composite fits of panel data.
4. `prepare_mxl_data(covariate_cols = NULL)` fails ("length(alt_idx) == nrow(X) is not TRUE").
5. `sample_by_choice(frac_per_alt =)` keeps round(f · available) per stratum, i.e. quotas, not
   Bernoulli sampling: worth a sentence in its documentation.
6. `ess()` and `mcse()` return NA for chains longer than about 65,000 draws: `.fft_autocov()`
   (`R/hb_diagnostics.R:64`) divides by `L * n`, an integer product (`nextn()` returns an integer)
   that overflows; `as.numeric(L) * n` fixes it.
7. From the plan: `summary.choicer_mc(level =)` is ignored; `monte_carlo()` on MNP/HB fits calls
   `logLik` outside its `tryCatch`.
