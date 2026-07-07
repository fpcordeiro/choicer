# HB Gibbs sampler profiling report

**Run:** `20260706-hb-profiling` (Phase 2 — full grid). **Target application:**
healthcare-provider-choice, cross-sectional (`person_col = NULL`, `T_i = 1`),
N ≈ 30,000, J ≈ 200 inside alternatives, K = 10 random-coefficient
covariates, P = 10 provider-level covariates. **Feasibility budget:**
R = 10,000, burn = 5,000, **2 sequential chains** (confirmed:
`R/hmnlogit_utils.R`'s and `hmnprobit_utils.R`'s `run_chain()` extra chains
are drawn via a plain `lapply`, not run in parallel). **Hardware:** Apple M3
Pro, 11 cores (5 performance), 36 GB RAM, R 4.6.0, macOS. **Models:**
`hmnl_gibbs` (`src/hmnlogit.cpp`), `hmnp_gibbs` (`src/hmnprobit.cpp`),
`run_mnprobit` contrast.

No package source was modified to produce this report — see the Methods
appendix (§7) for the instrumentation and reproduction procedure.

---

## 1. Hotspot breakdown (C++-block granularity)

All tables below are measured wall-clock `profile` block totals (seconds),
divided by `R`, from the **target-scale anchor** run (N=30000, J=200, K=10,
P=10, R=30, 5 threads) — the only run in this study conducted at the actual
target N and J. Percentages are each block's share of the anchor's
per-iteration total.

### HMNL (`hmnl_gibbs`) — per-iteration cost at N=30000, J=200, 5 threads

| block | sec/iter | share |
|---|---|---|
| (0) cache rebuild | 0.01427 | 9.1% |
| (a) beta RW-MH (per-respondent) | 0.02093 | 13.3% |
| (b) delta serial sweep | **0.12125** | **77.0%** |
| (c1) hierarchy draws (b/W/theta/sigma_d2) | 0.00092 | 0.6% |
| (c2) recording | 0.00001 | <0.1% |
| **total** | **0.15738** | 100% |

**Beta phase (a) sub-split** (thread-summed CPU time, at N=4000/J=200 — the
largest J-sweep point measured; not separately re-measured at the anchor,
see §2 for why the sub-split's *shape* is what matters, not its absolute
magnitude at N=4000):

| sub-block | measured behavior |
|---|---|
| proposal construction (Cholesky of `H_i+W⁻¹` + back-solve, O(K³)/respondent) | **flat in J** (J-slope 0.062) — confirms this term is J-invariant, as expected from fixed K=10 |
| candidate loglik pass (O(rows_i) = O(J) per respondent) | scales ~linearly with J (J-slope 0.901) |
| accept/reject + adaptation bookkeeping | nearly flat in J (J-slope 0.194), small magnitude |

This is the direct, measured confirmation of the leader's structural
hypothesis: `beta_mh`'s *total* J-slope (0.500, see §2) is a **blend** of a
J-invariant O(K³) term and a J-linear O(J) term — not a single scaling
regime — exactly the ambiguity flagged in Phase 1 as needing the real
J-sweep to resolve.

### HMNP (`hmnp_gibbs`) — per-iteration cost at N=30000, J=200, 5 threads

| block | sec/iter | share |
|---|---|---|
| (a) fused pass (TN sweep + beta_i draw + mu_x refresh, wall) | 0.04478 | 82.0% |
| (b) delta_j draws (parallel) | 0.00312 | 5.7% |
| (c) RSS pass (parallel) | 0.00220 | 4.0% |
| (d1) hierarchy draws | 0.00195 | 3.6% |
| (d2) recording | 0.00007 | 0.1% |
| startup Gram blocks (one-time, not per-iter) | n/a (once) | — |
| **total** | **0.05465** | 100% |

**Fused-pass (a) sub-split** (thread-summed CPU time; from the J-sweep at
N=4000, the largest measured J):

| sub-block | measured behavior |
|---|---|
| TN latent sweep (O(N·J)) | near-linear in J (J-slope 0.856) |
| conjugate beta_i regression draw (O(N·K²), K fixed) | markedly sub-linear in J (J-slope 0.392) — same K³-dominance signature as HMNL's Cholesky term |
| mu_x cache refresh (O(N·J)) | near-linear in J (J-slope 0.939) |

The instrumentation added in Phase 2 (mirroring HMNP's existing sub-split
pattern onto HMNL's beta phase) worked exactly as intended: both models show
the same qualitative signature — a K³-dominated, J-invariant sub-block
sitting alongside a J-linear sub-block — confirming this is a structural
feature of per-respondent RW-MH/regression updates with a fixed K, not an
artifact of one kernel's implementation.

---

## 2. N-scaling and J-scaling curves (fitted empirical exponents)

**Measured points:** N-sweep at J=50, N ∈ {1000, 2000, 4000, 8000}; J-sweep
at N=4000, J ∈ {25, 50, 100, 200}, all at R=200/burn=100/threads=5.
**Extrapolated:** the target-scale numbers in §1 and §5 (N=30000, J=200) are
the one MEASURED anchor point at that scale, not an extrapolation; the
*sub-split* percentages quoted in §1 for the beta/fused phases are read off
the N=4000/J=200 sweep point (the instrumentation for those specific
sub-splits was not separately re-run at the N=30000 anchor scale — their
*qualitative* pattern, not exact magnitude, is what's being asserted).
Fitted exponents come from `log(per-iteration block cost) ~ log(sweep var)`,
n=4 points per fit (full table: `results/analysis_scaling_slopes.csv`).

| model | sweep (other var fixed) | block | exponent | interpretation |
|---|---|---|---|---|
| HMNL | N (J=50) | total | 0.975 | linear in N — no N² surprises |
| HMNL | N (J=50) | delta_serial | 0.994 | linear in N, as expected (O(N·J)) |
| HMNL | J (N=4000) | cache_rebuild | 0.907 | linear in J, matches O(N·J) |
| HMNL | J (N=4000) | **delta_serial** | **1.085** | **linear-or-slightly-super-linear in J** — the serial bottleneck does not merely stay proportional, its *share* actively grows (see §3) |
| HMNL | J (N=4000) | beta_mh (total) | 0.500 | mixed regime (see sub-split above) |
| HMNP | N (J=50) | total | 0.956 | linear in N |
| HMNP | J (N=4000) | fused_pass (total) | 0.621 | sub-linear, driven by the K³-dominated beta_draw sub-component |
| HMNP | J (N=4000) | delta_parallel | 0.697 | sub-linear (thread-count fixed at 5 across the sweep; some load-imbalance at low J) |
| HMNP | J (N=4000) | rss | 0.815 | close to linear, some noise |

**Caveat on the J-sweep sub-unity exponents** (`hierarchy`, `recording`,
`rss`, `delta_parallel` for both models): these fits use only 4 points at
absolute magnitudes of a few to a few tens of milliseconds, where OS/
scheduler jitter is a non-trivial fraction of the signal. Treat any single
exponent outside {0, 1} with caution unless — like `beta_mh`/`fused_pass` —
it is corroborated by a sub-split that explains *why* (K³ vs J-linear
components). The N-sweep exponents are far more reliable (all magnitudes
≥ 1ms even at N=1000, and all land within [0.75, 1.10] of linear).

---

## 3. The serial delta sweep: CONFIRMED as the large-J bottleneck

**Verdict: CONFIRMED**, with numbers, both as a *share* and via a formal
Amdahl decomposition.

### Share of per-iteration wall time vs J (N=4000, 5 threads)

| J | cache_rebuild | beta_mh | **delta_serial** |
|---|---|---|---|
| 25 | 9.0% | 30.0% | **56.2%** |
| 50 | 8.8% | 21.8% | **66.6%** |
| 100 | 8.9% | 15.9% | **73.8%** |
| 200 | 8.6% | 12.4% | **78.1%** |

The serial delta sweep's share of total per-iteration wall time grows
monotonically from 56% at J=25 to **78% at J=200** — and reaches **77.0% at
the actual target scale** (N=30000, J=200; §1). This is not merely "the
biggest block" — its *dominance* is actively increasing with J, exactly the
architectural risk flagged in the leader's pre-analysis.

### Amdahl decomposition at N=4000, J=100 (threads 1 vs 5, measured)

| model | serial block | sec/iter @1thr | sec/iter @5thr | measured speedup | serial fraction @1thr | Amdahl ceiling (∞ threads) |
|---|---|---|---|---|---|---|
| **HMNL** | delta_serial | 0.01437 | 0.00794 | **1.81×** | 40.9% | **2.44×** |
| HMNP | hierarchy (master-only) | 0.01232 | 0.00299 | 4.12× | 0.85% | 118×|

HMNL is already close to Amdahl-saturated at 5 threads (1.81× of a 2.44×
ceiling — 74% of the achievable headroom already captured); **adding more
cores buys almost nothing further** for HMNL specifically because of this
one serial block. HMNP, by contrast, has essentially no serial bottleneck
(0.85% serial fraction, near-ideal 4.12× speedup on 5 threads, ceiling
118×) — it would keep scaling well with more cores, HMNL would not.

**Bottom line:** the serial delta sweep is confirmed, with growing severity
as J increases, as the single dominant architectural constraint on HMNL's
per-iteration cost at the target J=200 scale. It is *not* a fixable-by-more-
cores problem (see §6, lever #1, for the estimated ceiling of actually
addressing it).

---

## 4. Peak-memory findings

### Measured RSS across the grid

Peak RSS was measured via `/usr/bin/time -l` (macOS "maximum resident set
size", reported in bytes) around the **entire worker process**, so it
includes R startup, package load, `simulate_h*_data()`/`prepare_h*_data()`
data generation/preparation, and the Gibbs run itself — not just the C++
kernel's internal working set.

| config | N | J | R | peak RSS |
|---|---|---|---|---|
| nsweep (J=50) | 1000→8000 | 50 | 200 | 0.19 → 0.75-0.78 GB (both models, near-identical) |
| jsweep (N=4000) | 4000 | 25→200 | 200 | 0.28 → 1.06-1.18 GB |
| **target anchor** | **30000** | **200** | 30 | **HMNL 5.92 GB · HMNP 6.54 GB** |

The measured target-scale peak RSS (~6-6.5 GB, `keep_beta_i="means"`) is
**far above** a pure C++-working-set estimate (X ≈ 480 MB + 3-4 candidate/
cache vectors of ≈48 MB each + a K×K×N cube ≈ 24 MB ≈ under 1 GB total).
This gap is dominated by the **R-side DGP and data-prep overhead**
(`simulate_hmnl_data()`/`prepare_hmnl_data()` building and rbinding a
~6-million-row `data.table`, plus transient copies) — a real number to
budget for (`prepare_h*_data()` runs regardless of DGP), but specific to
how data reaches the sampler, not to the sampler's own steady-state memory.

### `keep_beta_i = "draws"` memory probe vs the analytic formula

Probe: N=8000, J=50, R=400/burn=200 (R_keep=200), compared against the
matched `keep_beta_i="means"` config at the same N, J (isolates the
incremental cost of the K×N×R_keep draw cube from the base N,J-dependent
footprint):

| model | measured Δ (cube-attributable) | analytic K·N·R_keep·8B | ratio |
|---|---|---|---|
| HMNL | 244.8 MB | 128.0 MB | **1.91×** |
| HMNP | 220.0 MB | 128.0 MB | **1.72×** |

**The naive analytic formula underestimates actual memory by ~1.7-1.9×** —
almost certainly the `Rcpp::wrap()` copy of the C++ `arma::cube` into a new
R array plus R's list-construction overhead, which the analytic formula
does not account for. **Recommendation: budget ~2× the naive
K·N·R_keep·8B figure for `keep_beta_i="draws"`.**

### Projected RSS at target scale, `keep_beta_i = "draws"`

At the target budget (K=10, N=30000, R_keep = (10000-5000)/1 = 5000):

- Naive analytic cube: 8 × 10 × 30000 × 5000 = **12.0 GB**
- With the measured ×1.7-1.9 correction: **≈ 20.4 - 22.8 GB**, ON TOP of the
  ≈ 6-6.5 GB base footprint → **≈ 26.4 - 29.3 GB for a single chain**.

**Critical compounding risk for 2 sequential chains:** `run_hmnlogit()`/
`run_hmnprobit()` run the second chain via `lapply()` in the *same* R
process, and the first chain's complete result object (`out`, including its
full `beta_i_draws` cube) remains referenced in the calling frame — R's
garbage collector cannot reclaim it — while the second chain runs and grows
its own cube. **Peak RSS during chain 2 could therefore be ≈ 2 × (26.4-29.3
GB) ≈ 53-59 GB — exceeding the target machine's 36 GB RAM.**

### Memory or time — which binds first?

**Depends entirely on `keep_beta_i`:**
- **Default (`keep_beta_i = "means"`):** peak RSS ≈ 6-6.5 GB per chain,
  ≈ 12-13 GB even if (hypothetically) both chains' results were
  simultaneously retained — comfortably within 36 GB. **Time binds first**
  (§5: ~30-55 min full budget).
- **`keep_beta_i = "draws"` with 2 chains:** projected peak RSS (~53-59 GB)
  **exceeds available RAM** — **memory binds first**, and would need to be
  addressed (streaming to disk, or running each chain in a fresh
  subprocess) before a 2-chain `keep_beta_i="draws"` run at this scale is
  safe to attempt.

---

## 5. Per-model feasibility at target scale + full budget

Extrapolation method: the **measured** target-scale anchor (N=30000, J=200,
R=30) per-iteration cost, multiplied by the full R=10,000 budget (linear in
R by construction — R only appears as a loop count, not as a scale
parameter of any block) and by 2 sequential chains. This is a measured
base rate scaled linearly, not a cross-(N,J) power-law guess.

| model | sec/iter @ target | per chain | both chains (sequential) |
|---|---|---|---|
| **HMNL** | 0.1574 | 1573.8 s = **26.2 min** | 3147.6 s = **52.5 min** |
| **HMNP** | 0.0546 | 546.5 s ≈ **9.1 min** | ~1093 s ≈ **18.2 min** |

*(Phase-1's pre-grid analytical prediction was 35.5-50.5 min for HMNL and
~17.3 min for HMNP — the measured anchor lands close to, and for HMNL just
above, that range: a good validation of the Phase-1 methodology, and
confirmation that HMNL's `beta_mh` block leans toward the O(N·J) end of its
ambiguous range at J=200, consistent with its measured J-slope of 0.5.)*

**Both models comfortably fit inside a single overnight run, or well within
an hour on this hardware, at the default `keep_beta_i="means"` setting.**
HMNL is ~3× slower than HMNP per iteration at target scale — entirely
attributable to the serial delta sweep (§3) — but neither model is anywhere
near an infeasible wall-clock regime.

### HMNP Albert-Chib tractability at N·J = 6,000,000 latents/sweep

**Tractable, with clear numbers.** The fused pass (TN latent sweep +
conjugate beta_i draw + mu_x refresh) over all 6 million (N×J) latent
utilities plus their per-task outside latents completes in **0.0448
sec/iteration** at 5 threads (measured at the anchor) — i.e., roughly **134
million latent-adjacent operations per second** in aggregate across the
fused pass's three sub-phases. The Albert-Chib augmentation itself is not a
scaling risk at this N×J; HMNP's total per-iteration cost is dominated by
this fully-parallel pass (82% share, §1), which is exactly the well-behaved
(near-linear-in-N-and-J, near-ideal thread speedup) regime confirmed in
§§2-3.

---

## 6. `run_mnprobit` contrast: not viable at J=200

**Measured** (N=2000 fixed, J ∈ {4, 8, 12}, R=200, wall-time only, not
separately instrumented per spec):

| J | p = J-1 | sec/iter | peak RSS |
|---|---|---|---|
| 4 | 3 | 0.000405 | 107 MB |
| 8 | 7 | 0.000653 | 117 MB |
| 12 | 11 | 0.000970 | 139 MB |

A naive power-law fit to just these three points gives a shallow exponent
(≈0.66 in p) — but **this badly underestimates the true asymptotic cost**,
because at p ≤ 11 the model is still solidly in a regime where the O(N·p)
work of the fully-differenced latent sweep dominates a comparatively
negligible O(p³) inverse-Wishart/Cholesky term (a p³ FLOP count at p=11 is
~1,300 flops — trivial next to N·p ≈ 22,000 truncated-normal draws).

**Reading `src/mnprobit.cpp`'s latent-sweep loop directly** (not guessing)
shows each observation's TN-conditional draw requires a **full p×p
precision-matrix reduction** (`Ratio(k,j)` summed over all `k = 0..p-1` for
every `j`), because `run_mnprobit` differences ALL utilities against a base
alternative with a *general* (J-1)×(J-1) covariance Σ — unlike HMNP's iid
per-shock model, which needs no such cross-dimension conditioning. This
makes the latent sweep **O(N·p²) per iteration**, not O(N·p).

Fitting `sec_per_iter = c·p²` (grounded in that code-level structure, not a
blind exponent search) to the 3 measured points gives c ≈ 8.94e-6,
predicting at the target-adjacent (N=30000, p=199 i.e. J=200):

**≈ 5.3 sec/iteration — about 100× slower per iteration than HMNP at the
same (N, J)** (0.0546 sec/iter, §5). A single 10,000-iteration chain would
take **~885 minutes (≈ 14.75 hours)** — before even accounting for the
`(J-1)×(J-1)` Σ's own O(p³) inverse-Wishart draw or its `O(p²K²)`
beta-precision assembly (both noted in the kernel's own comments as
additional, non-negligible costs at large p).

**Verdict: `run_mnprobit` is not a viable engine at J≈200 for this
application.** This is a structural conclusion (fully-differenced,
general-covariance MNP inherently costs O(N·p²) per iteration in its latent
sweep), not merely a matter of insufficient hardware — it is exactly why
the hierarchical models (HMNL/HMNP) with their scalar or block-structured
alternative effects were the right architectural choice for J~200. (Caveat:
the O(N·p²) coefficient was fit to only 3 low-p points with 1 residual
degree of freedom; a moderate-J, e.g. J=30-50, `run_mnprobit` run would
sharpen this considerably if a precise number is ever needed — but the
qualitative conclusion, "1-2 orders of magnitude slower than HMNP at
J=200", is robust to that uncertainty.)

---

## 7. Predicted-vs-actual (Phase 1 → Phase 2 validation)

All 27 grid configs' actual `profile$t_total` were compared against the
Phase-1 analytical pre-grid prediction (`configs/predicted_grid_times.csv`,
computed *before* any full-grid execution). **Zero configs exceeded a 3×
deviation** in either direction (full table:
`results/analysis_predicted_vs_actual.csv`); actual/predicted-high ratios
ranged from 0.49× to 1.48× across all 27 configs — the scaling model used
for the pre-grid prediction held up well, with the target anchor landing
almost exactly at its predicted-high bound (ratio 1.04×) as expected given
the beta_mh ambiguity resolved toward its O(N·J) end (§2).

---

## 8. Ranked shortlist of optimization levers (no implementation)

Each lever below is backed by a measured number from this study; ceiling
estimates are explicitly flagged as rough (back-of-envelope from the
measured Amdahl/share data), not engineering commitments.

1. **Parallelize (or otherwise restructure) HMNL's serial delta sweep — highest-value lever.**
   Measured: 77.0% of per-iteration wall time at the target scale (§1);
   share grows monotonically with J (56%→78% from J=25→200, §3); Amdahl
   ceiling at infinite threads is only 2.44× vs the 1.81× already captured
   at 5 threads (§3) — **more cores alone cannot help further**. **Rough
   ceiling if this block matched the ~4.3× speedup already achieved by
   HMNL's OTHER (already-parallel) blocks:** total per-iteration cost could
   drop from the measured 0.1574 s/iter toward roughly 0.097 s/iter (≈38%
   reduction), moving the full 2-chain budget from ~52.5 min toward
   ~32-35 min. This requires an algorithmic change (the softmax-coupling
   that currently forces seriality is a genuine mathematical constraint,
   documented in `hb_internal.h`/`hmnlogit.cpp` — not an oversight), e.g. a
   different conditional structure, batching, or an approximation scheme.

2. **No action needed for HMNP's core loop.** Amdahl ceiling 118× at
   infinite threads, 0.85% serial fraction, near-ideal 4.12× measured
   speedup at 5 threads (§3) — already well-parallelized; more cores would
   continue to help roughly linearly.

3. **Avoid `keep_beta_i = "draws"` with 2 sequential chains at target scale
   without a code change.** Measured: the naive K·N·R_keep·8B formula
   underestimates actual RSS by 1.7-1.9× (§4); projected single-chain peak
   ≈ 26-29 GB, and because chain 1's result stays referenced in R while
   chain 2 runs, 2-chain peak RSS could reach ≈ 53-59 GB — **exceeding the
   target machine's 36 GB RAM**. **Ceiling:** avoiding this failure mode
   entirely (streaming per-chain draws to disk, or running each chain as a
   separate subprocess so chain 1's memory is released before chain 2
   starts) removes the risk at effectively zero runtime cost — this is a
   correctness/safety lever, not a speed lever.

4. **Do not use `run_mnprobit` for J≈200 problems.** Measured+structural:
   ~100× slower per iteration than HMNP at the same (N,J) (§6), rooted in
   an O(N·p²) latent-sweep cost inherent to its fully-differenced,
   general-covariance formulation. **Ceiling: not fixable by tuning** —
   this is the right reason to prefer HMNP's architecture for this
   application, already the apparent choice; no further action, just
   confirmation.

5. **(Cleared, not a lever) HMNL's beta-phase proposal construction
   (Cholesky, O(K³)) is J-invariant** (J-slope 0.062, §1) — there is
   nothing to gain here as J grows; correctly not a target for
   optimization. Included for completeness: a lever explicitly checked and
   found not to matter, rather than left unexamined.

---

## 9. Methods appendix

### Instrumentation

- All timers are `#ifdef HB_PROFILE`-guarded blocks added ONLY to a
  scratch copy of the package source (never the repo), captured as
  `_benchmarks/hb_profiling/patches/hb_profile.patch` (a unified diff).
  Phase 1 added master-thread wall-clock timers at existing barrier/phase
  boundaries in `hmnl_gibbs`/`hmnp_gibbs`, plus HMNP's fused-pass
  thread-local CPU-time sub-split. Phase 2 added the analogous
  thread-local sub-split to HMNL's beta phase (proposal construction /
  candidate loglik / accept-adapt), mirroring HMNP's existing pattern.
- No timer touches an RNG stream (`Xoshiro256pp` draws untouched) or
  reorders computation.

### Bitwise-gate evidence

- Tiny fixed-seed config (N=200, J=10, K=4, P=3, R=200, burn=100, seed=42)
  run through a non-instrumented reference build and the instrumented
  build (separate subprocesses, identical prepared inputs). **Result: all
  of `bdraw`/`wdraw`/`deltadraw`/`thetadraw`/`sigma_d2draw`
  (+`sigma2draw` for HMNP) bitwise IDENTICAL (`max_abs_diff == 0`), for
  both kernels, both BEFORE and AFTER the Phase-2 beta-phase sub-split
  extension.** Re-gated after the Phase-2 patch extension, per the
  coordinator's requirement, before the full grid was run.
- Timer overhead: negligible relative to per-iteration work at production
  scale (a handful of `omp_get_wtime()` calls per iteration vs O(N·J) work
  per iteration); measured overhead at a tiny smoke scale (R=2000,
  N=200/J=10) was ±7-15%, dominated by run-to-run OS noise at that tiny
  absolute magnitude (sub-200ms total), not a systematic timer cost.

### Hardware, seeds, reproduction

- Apple M3 Pro, 11 cores (5 performance), 36 GB RAM, R 4.6.0, macOS.
  `choicer::set_num_threads()` used to fix thread counts per config.
- All configs use fixed seeds (`configs/*.csv`'s `seed` column; kernel
  seeds offset deterministically from it) for reproducibility.
- **To reproduce:**
  1. `Rscript _benchmarks/hb_profiling/00_build_instrumented.R` — builds
     the instrumented library from the current patch.
  2. `Rscript _benchmarks/hb_profiling/01_bitwise_gate.R --lib-dir <templib> --ref-lib-dir <templib_ref>`
     — must show `max_abs_diff == 0` for both kernels before proceeding.
  3. `Rscript _benchmarks/hb_profiling/20_harness.R --grid _benchmarks/hb_profiling/configs/full_grid_execution_order.csv --lib-dir <templib>`
     — runs all 27 configs (target-scale anchor last).
  4. `Rscript _benchmarks/hb_profiling/30_analyze.R` — produces the
     `results/analysis_*.csv` tables this report is built from.
- Raw results: `_benchmarks/hb_profiling/results/hb_profiling_results.csv`
  (27 rows). Analysis tables:
  `results/analysis_scaling_slopes.csv`, `results/analysis_amdahl.csv`,
  `results/analysis_memory_probe.csv`,
  `results/analysis_predicted_vs_actual.csv`,
  `results/analysis_target_extrapolation.csv`.

### Scope note

This report is evidence-gathering only, per the run's request: no source
changes, no optimization implementation. All findings above are backed by
measured numbers from this repo's actual kernels at the stated hardware and
configuration; extrapolations are explicitly labeled as such, with the
target-scale numbers in §§1, 4, 5 being *measured directly* (the anchor
config), not extrapolated.
