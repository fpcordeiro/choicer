# HB Bayesian engine — two plans: runtime/memory optimization & convergence diagnostics

**Scope.** Planning only. No source changes until approved. Both plans target
the hierarchical-Bayes engine: `hmnl_gibbs` (`src/hmnlogit.cpp`), `hmnp_gibbs`
(`src/hmnprobit.cpp`), their wrappers `run_hmnlogit()` / `run_hmnprobit()`
(`R/hmnlogit_utils.R`, `R/hmnprobit_utils.R`), the shared `choicer_hb`
post-estimation (`R/hb_postest.R`, `R/methods.R`), the existing diagnostics
(`R/hb_diagnostics.R`), and the `run_mnprobit` contrast (`src/mnprobit.cpp`).

**Evidence base.** Plan 1 is grounded strictly in
`_benchmarks/hb_profiling_report.md` (run `20260706-hb-profiling`, Apple M3 Pro,
5 performance cores, 36 GB, R 4.6.0). Every lever below cites a measured number
from that report; nothing is a speculative optimization. Where I add a lever the
report's §8 shortlist does not name explicitly (Lever G), it is derived
arithmetically from the report's own §1/§3 Amdahl numbers, and I show the
derivation.

**Reproducibility contract (applies to both plans).** The engine's guarantees —
bitwise draw-equality independent of thread count, via one persistent OpenMP
region, fixed-order cross-unit reductions, and a per-(iteration, unit) RNG
partition `make_stream(seed, r, tag)` — are a hard constraint. Every Plan-1 lever
carries an explicit gate stating whether it (a) leaves the posterior *bitwise*
identical, (b) leaves the posterior *law* identical but changes the realized
sample path (requires statistical re-validation, not bitwise), or (c) changes the
target (rejected unless the change is a provably-correct MCMC operator with
quantified error).

---

# PLAN 1 — Runtime & memory optimization

## 1.1 Success target

Target application (from the report header): healthcare-provider choice,
cross-sectional (`person_col = NULL`, `T_i = 1`), **N ≈ 30,000**, **J ≈ 200**
inside alternatives, **K = 10** random-coefficient covariates, **P = 10**
provider-level covariates, **R = 10,000** (burn = 5,000), **2 chains**, on the
**5 performance cores** of an M3 Pro with **36 GB** RAM.

Measured baseline (report §5, §4), `keep_beta_i = "means"`:

| model | 2-chain wall-clock (sequential) | peak RSS / chain |
|---|---|---|
| HMNL | **52.5 min** | 5.92 GB |
| HMNP | **18.2 min** | 6.54 GB |

**Success target (this plan):**

> Fit N ≈ 30,000, J ≈ 200, K = 10, P = 10, R = 10,000 × **2 chains** on **5
> cores** in **under 30 min for HMNL** and **under 15 min for HMNP**, at **peak
> RSS under 8 GB** in the default `keep_beta_i = "means"` mode; and make a
> 2-chain `keep_beta_i = "draws"` run at this scale **either feasible under 36 GB
> or fail fast with an accurate memory budget** (never a silent 53–59 GB blow-up).

Stretch (if Lever A also lands): HMNL 2-chain **under 22 min**.

Rationale for the numbers: HMNL's 52.5 min is dominated by a serial block that
leaves ~66% of the 5-core capacity idle (derived in §1.3, Lever G); reclaiming
that idle capacity is the ~2× that gets HMNL under 30 min at essentially zero
correctness risk. HMNP is already near-ideal; "under 15 min" is a modest tighten
from adding cores, not a code change. 8 GB is chosen above the measured 6–6.5 GB
means-mode footprint with headroom for two simultaneous chains.

## 1.2 Ranking methodology

The brief asks to rank by **measured hotspot share × scaling exponent × expected
payoff × (1 / implementation risk)**. I score each lever on those four factors
(1–5 each; risk enters inverted, so "low risk" scores high), multiply, and sort.
Hotspot share and scaling exponent come directly from report §1–§3; payoff is the
report's own ceiling estimate; risk reflects whether the posterior stays bitwise
identical (low) vs. requires a new MCMC operator and full re-validation (high).

| # | Lever | Hotspot share | Scaling exp. | Payoff | Risk (inv.) | Score | Correctness class |
|---|---|---|---|---|---|---|---|
| **G** | Run the 2 chains in parallel processes (HMNL) | 4 (exploits the 77% serial block's idle cores) | 4 | 4 (~2× on 2-chain budget) | 5 (kernel untouched) | **320** | (a) bitwise |
| **A** | Restructure HMNL's serial `delta` sweep | 5 (77%, growing) | 5 (J-exp 1.085) | 3 (block ceiling 2.44× → ~38% total) | 2 (new operator) | **150** | (b) law-preserving / (c) if approx |
| **B** | `keep_beta_i="draws"` memory safety + streaming | 2 | 3 | 3 (removes 53–59 GB failure mode) | 4 | **72** | (a) bitwise |
| **C** | Cut R-side data-prep peak RSS | 2 | 3 | 2 (memory only; means-mode already fits) | 4 | **48** | (a) bitwise |
| **D** | HMNP: add cores, no code change | 5 (82% but already parallel) | 1 (nothing to fix) | 2 | 5 | **100**→N/A | confirm only |
| **E** | Do **not** use `run_mnprobit` at J≈200 | — | — (O(N·p²) structural) | — | — | N/A | architectural, no action |
| **F** | HMNL beta-phase Cholesky O(K³) | 0 (J-invariant, slope 0.062) | 0 | 0 | — | **0** | cleared, no action |

Ranked order for action: **G → A → B → C**, with **D** as a
one-line configuration note, **E/F** as confirmations. G outranks A on
payoff-per-risk even though A targets the larger hotspot: G *hides* the 77% block
behind idle cores at zero correctness risk, A *shrinks* it but needs a new
sampler operator and statistical re-validation. They compound — do G first, then A
raises the ceiling further.

## 1.3 Lever detail

### Lever G — parallelize the two chains (HMNL) — **rank 1**

**What it changes.** R-orchestration only. Today (confirmed at
`R/hmnlogit_utils.R:413-420` and `R/hmnprobit_utils.R:246-253`) the extra chains
run via a plain `lapply` in the *same* R process, strictly sequentially. Change:
resolve the per-chain seeds in the parent (`seed`, `seed+1`, …), then dispatch
each chain to a separate process (e.g. `parallel::mclapply` / a `future` plan),
collect results, assemble as today.

**Why the evidence supports it.** From report §1, HMNL per-iteration at target =
0.1574 s, of which the serial `delta` block is 0.121 s (77%) on **one** thread
while the other four cores idle; the parallel blocks (cache 0.014 + beta 0.021 +
hierarchy 0.001 ≈ 0.036 s) use all five. Core-seconds available per iter = 5 ×
0.1574 = 0.787; useful work ≈ 0.121 (serial, 1 core) + ~0.15 (parallel) ≈ 0.27 →
**~34% core utilization**. The report's own Amdahl decomposition (§3) says HMNL
captures only **1.81× of a 2.44× ceiling** at 5 threads — i.e. more cores inside
one chain buy almost nothing. Two chains time-sharing the machine fill the idle
~66%: their serial `delta` blocks run concurrently on two cores while the parallel
blocks fill the rest, so the 2-chain wall-clock collapses toward a *single*
chain's time.

**Estimated payoff.** 2-chain HMNL from 52.5 min toward **~27–30 min** (~2×). No
gain for HMNP (report §3: 4.12× of a 118× ceiling — cores already saturated
within a chain; two HMNP chains would merely split the same busy cores).

**Effort.** Low–medium (R only; no kernel edit).

**Correctness gate — class (a), bitwise.** The kernel is untouched and is
`[[Rcpp::export(rng = false)]]`: its draws depend only on the explicit `seed`
argument via `make_stream`, never on R's RNG state. As long as the parent
resolves each chain's seed *before* dispatch and passes it explicitly, every
chain's draws are **bitwise identical** to the current sequential run, and the
per-(iteration,unit) RNG partition and thread-count invariance are preserved
untouched. **Pitfall to gate in review:** do not let the parallel backend reseed
or fork-inherit an R RNG that feeds the default-seed `sample.int()` — resolve
seeds in the parent first. **Validation:** run 1 vs 2-process execution on a fixed
small config; assert `max_abs_diff == 0` on all draw matrices (the exact
bitwise-gate protocol from report §9). Memory: two means-mode chains = ~2 × 6 GB ≈
12 GB < 36 GB (fits); and separate processes *release* chain-1's memory at join,
which also defuses Lever B's compounding risk for free.

### Lever A — restructure HMNL's serial `delta` sweep — **rank 2**

**What it changes.** The `delta_j` update loop (`src/hmnlogit.cpp:665-745`),
currently a strictly serial RW-MH sweep. It is serial by necessity: the `delta_j`
conditionals are coupled through the shared per-task softmax denominators `D_t`
(documented at `hmnlogit.cpp:31-38` and in `hb_internal.h`), so a naive
work-shared sweep would not leave the posterior invariant.

**Why the evidence supports it.** Report §3: this block is 77.0% of per-iteration
wall time at target, its share **grows** monotonically with J (56%→78%, J 25→200),
its J-exponent is 1.085 (linear-to-super-linear), and it is the sole reason HMNL
is ~3× slower than HMNP and Amdahl-capped at 2.44×. It is *the* structural
constraint.

**Three sub-options, in increasing effort / decreasing risk-of-bias:**

- **A1 — batched parallel proposals with deferred denominator refresh
  (approximate).** Propose all `delta_j` against stale caches, accept/reject in
  parallel, then rebuild `D_t` once. **Class (c): changes the target** (the
  accept ratios use stale denominators). *Flagged: fails the "posterior provably
  unchanged" gate* unless recast as a valid Metropolis-within-Gibbs with a
  correct joint accept. Recommend **reject** unless bias is formally bounded and
  measured negligible.
- **A2 — blocked / joint `delta` Metropolis (exact).** Propose a block of
  `delta_j` from a Gaussian approximation to their joint conditional; one joint
  accept step. Likelihood evaluation parallelizes; the operator still targets the
  exact posterior. **Class (b): law-preserving, sample path changes.** Medium–high
  effort.
- **A3 — Pólya-Gamma augmentation for the logit `delta` block (exact, best
  asymptotics).** PG data augmentation renders the logit `delta` conditional
  Gaussian → a conjugate, *work-shareable* draw exactly like HMNP's parallel
  `delta` pass, collapsing the serial block. Already on the roadmap
  (`ROADMAP.md` / `BAYESIAN_PLAN.md` note "Pólya-Gamma sampler"). **Class (b).**
  High effort (new kernel machinery), highest structural payoff.

**Estimated payoff.** Report §8 lever #1: if the block reached the ~4.3× speedup
of HMNL's already-parallel blocks, per-iter drops 0.157→~0.097 s (~38%), 2-chain
budget 52.5→~32–35 min (before Lever G). A3 could beat this (removes the serial
regime entirely) but is a rewrite.

**Effort.** A1 low but rejected; A2 medium–high; A3 high.

**Correctness gate — class (b).** Any A-variant changes the realized draws — it is
**not** bitwise-comparable to the current serial sweep. Gate: (1) the new operator
must provably target the same conditional (A2/A3 do; A1 does not); (2) re-establish
thread-count invariance for the new parallel structure (fixed-order reductions,
per-(iteration,unit) `make_stream` tags — reuse the HMNP `delta` pass's proven
pattern at `hmnprobit.cpp:547-571`); (3) validate the *law* statistically, not
bitwise: `recovery_table()` on `simulate_hmnl_data()` truth must still cover, and
posterior mean/SD/quantiles of `b, theta, delta, sigma_d2` must agree with the
current serial implementation within MCMC error on a fixed dataset (pre-vs-post
posterior comparison). A1 additionally requires an explicit bias measurement.

### Lever B — `keep_beta_i="draws"` memory safety + streaming — **rank 3**

**What it changes.** The `beta_i` draw-cube path. Two concrete fixes: (i) make the
R-side guard (`R/hmnlogit_utils.R:341-354`, `R/hmnprobit_utils.R:214-227`)
budget-accurate; (ii) offer a streaming-to-disk (or per-chain-subprocess) path so
draws are obtainable at scale without a blow-up.

**Why the evidence supports it.** Report §4: the naive `K·N·R_keep·8B` formula
underestimates true RSS by **1.7–1.9×** (the `Rcpp::wrap` cube copy + R list
overhead); projected single-chain peak ≈ 26–29 GB, and because chain 1's result
stays referenced while chain 2 runs in the same process, 2-chain peak could reach
**53–59 GB > 36 GB**.

**Important nuance the report frames as live but the code already blocks.** The
current guard hard-stops when `bytes > 4e9`: at target (K=10, N=30,000, R_keep=5,000)
`bytes = 8·10·30000·5000 = 12 GB > 4 GB` → `stop()` **before the chain runs**. So
the catastrophic 53–59 GB path is *already unreachable* at target scale under
present code. Lever B is therefore: (a) correct the guard's *message* to reflect
the measured ~2× wrap overhead and the 2-chain retention (so users get an honest
number), and (b) provide a genuine large-scale draws path (stream each kept slice
to disk, or run each chain in a fresh subprocess — which Lever G already gives) for
users who truly need per-respondent draws at N=30k.

**Estimated payoff / effort.** Memory-safety, not speed. Low–medium effort.

**Correctness gate — class (a), bitwise.** Memory layout / IO only; the draw
*values* are unchanged. Streaming must byte-preserve each slice. Validation: reload
streamed draws and assert equality with an in-memory reference at a small scale.

### Lever C — cut R-side data-prep peak RSS — **rank 4**

**What it changes.** `simulate_h*_data()` / `prepare_h*_data()`
(`R/hb_data.R`, `.prepare_hb_panel`), which build and `rbind` a ~6-million-row
`data.table` plus transient copies.

**Why the evidence supports it.** Report §4: measured target-scale peak (6–6.5 GB)
is *far above* the C++ working set (<1 GB) and is dominated by this R-side DGP /
prep overhead — "a real number to budget for, but specific to how data reaches the
sampler, not the sampler's steady state."

**Estimated payoff / effort.** Memory only. **Low priority:** at means-mode the
6.5 GB already fits 36 GB comfortably and *time* binds first (report §4). Worth
doing only to enable draws-mode headroom or very-large-N. Medium effort
(in-place `data.table` construction, avoid `rbind` copies).

**Correctness gate — class (a), bitwise.** Prep must emit identical
`X / Z / M / choice_pos / alt_of_row / Ti`. Validation: assert prepared-object
equality (all.equal on every field) against the current path on a fixed seed →
posterior is then bitwise identical downstream.

### Lever D — HMNP core loop: no code change

Report §2–§3: 0.85% serial fraction, near-ideal 4.12× on 5 threads, Amdahl ceiling
118×. Action: **none** beyond documenting that HMNP scales ~linearly with added
cores (`set_num_threads()`), unlike HMNL. Confirmation lever.

### Lever E — `run_mnprobit` is not a viable engine at J≈200

Report §6: ~100× slower per iteration than HMNP at matched (N, J), rooted in an
**O(N·p²)** latent sweep inherent to its fully-differenced, general-covariance
formulation (`src/mnprobit.cpp` `Ratio(k,j)` reduction). A single 10k chain ≈ 14.75 h.
**Not fixable by tuning** — architectural. Action: none; this confirms HMNP's
scalar/iid structure was the right choice for J≈200.

### Lever F — HMNL beta-phase Cholesky: cleared

Report §1/§8 #5: the O(K³) proposal construction is J-invariant (J-slope 0.062).
Nothing to gain as J grows. Explicitly checked and excluded.

## 1.4 Per-model feasibility at target scale

- **HMNL — feasible, time-bound.** 52.5 min (2 chains) today; **Lever G alone
  reaches the <30 min target** at zero correctness risk; Lever A2/A3 pushes toward
  <22 min. Memory comfortable in means-mode (6 GB/chain). The serial `delta` sweep
  is the only structural constraint and is *not* a more-cores problem (Amdahl 2.44×).
- **HMNP — feasible, already well inside budget.** 18.2 min (2 chains), near-ideal
  parallel scaling; the Albert–Chib augmentation over 6M latents/sweep is tractable
  (report §5: 0.0448 s/iter fused pass, ~134M latent-ops/s). Adds cores linearly.
  No structural work needed.
- **`keep_beta_i="draws"` at target — blocked today, needs Lever B** to become
  either feasible (streaming / per-chain subprocess) or fail-fast with an honest
  budget. Means-mode is the supported production path at this scale.

---

# PLAN 2 — Convergence diagnostics

## 2.1 Inventory of what already exists

| Component | Location | Covers |
|---|---|---|
| `rhat()` (split-R̂, exported) | `R/hb_diagnostics.R:26` | Any draw matrix / list of matrices; split-half between/within variance. **Not** rank-normalized. |
| `.hb_rhat_table()` (internal) | `R/hb_diagnostics.R:76` | Builds split-R̂ over **b, theta, sigma_d2 only** — *not* `delta`, *not* `W`, *not* `sigma2`. |
| `ppc_shares()` (exported) | `R/hb_diagnostics.R:115` | Posterior-predictive choice-share check (fit-quality, not convergence). |
| Acceptance rates | `hmnl_gibbs` → `object$accept` (`R/hmnlogit_utils.R:526`) | HMNL only: per-`beta_i`, per-`delta_j`, means, final proposal scales. HMNP has none (fully conjugate). |
| Fit-time warning | wrappers (`hmnlogit_utils.R:502`, `hmnprobit_utils.R:345`) | Warns if any tracked R̂ > 1.05, or acceptance outside (0.10, 0.60). |
| `summary.choicer_hb` display | `R/methods.R:1950`, `:1989` | Prints b/theta/sigma_d2(/sigma2) tables, quality ladder, mean acceptance (HMNL), split-R̂ (if present). |

**Two structural gaps found while reading the code:**

1. **Only chain 1's draws are retained.** In both wrappers, `extra_chains` are
   consumed to build the R̂ table (`hmnlogit_utils.R:442`, `hmnprobit_utils.R:292`)
   and then **discarded**; `object$draws` holds chain 1 only. Any *multi-chain*
   diagnostic (ESS bulk/tail, multi-chain traceplots, rank-normalized R̂) requires
   retaining all chains' hierarchical draws.
2. **The hardest-mixing block is undiagnosed.** `delta_j` is exactly the parameter
   most at risk (RW-MH, serial sweep, J≈200 of them) yet it is absent from
   `.hb_rhat_table()`. Diagnostics must reach it.

## 2.2 The basic set (what a practicing econometrician wants — and no more)

Add exactly four things; reuse everything above.

1. **Effective sample size — bulk and tail.** Rank-normalized ESS
   (Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021, *Bayesian Analysis*):
   `ess_bulk` (rank-normalized draws) and `ess_tail` (min of the 5% and 95%
   quantile-indicator efficiencies). Multi-chain aware (pools within/between like
   R̂). This is the single most useful "is 10,000 iterations actually enough?"
   number, and it is what makes MCSE meaningful.
2. **Monte Carlo standard error (MCSE) on posterior summaries.** MCSE(mean) =
   posterior SD / √ess_bulk; report alongside each posterior mean so a user can
   see whether a reported coefficient is resolved to the precision they quote.
   Optionally MCSE(median) via the tail ESS.
3. **Traceplots.** A `traceplot()` generic + `traceplot.choicer_hb` method: chains
   overlaid, user-selectable parameter blocks (`b`, `theta`, `sigma_d2`, and a
   representative subset of `delta`). **Base `graphics`/`grDevices` only — no new
   dependency** (this would be the package's first plot method; it must stay
   CRAN-lean, consistent with HMNL/HMNP adding no new hard dependency).
4. **Acceptance rates.** Already computed for HMNL — surface them in the
   consolidated table. **HMNP: explicitly reported as "n/a — fully conjugate, no
   Metropolis step"** (not blank, so the user knows it is by design).

**Also (small, high-value):** upgrade the R̂ in the diagnostic table to the same
rank-normalized definition (pairs correctly with ESS) and **extend coverage to the
`delta` block** — but summarized (worst-case max R̂ / min ESS across the J
alternatives plus the offending index), never J≈200 raw rows.

## 2.3 Consolidated `summary()` diagnostic table (multi-chain aware)

One table, appended to `print.summary.choicer_hb`, driven by all retained chains:

```
Convergence diagnostics (C chains, R_keep draws each)
Block            R-hat   ESS_bulk  ESS_tail   MCSE(mean)
b[x1]            1.002      4180      3920      0.0011
b[...]            ...        ...       ...        ...
theta[...]        ...        ...       ...        ...
sigma_d^2        1.004      3010      2650      0.0008
delta (J=200)    1.031*     612*      540*      —          *worst: delta[hosp_44]
Acceptance: beta 0.28, delta 0.41          (HMNP: conjugate — no acceptance step)
```

- Per-parameter rows for the small blocks (`b`, `theta`, `sigma_d2`); a single
  **worst-case summary row** for `delta` (and optionally `W`) so J≈200 never
  floods the console, with the offending parameter named.
- Multi-chain when `chains > 1`; single-chain falls back to the within-chain split
  (as `rhat()` already does).
- Keep the fit-time warning but broaden it to the Vehtari thresholds
  (warn if any R̂ > 1.01 **or** min ESS_bulk < 400), including the `delta` block.

## 2.4 Files that change and new exports

| File | Change |
|---|---|
| `R/hb_diagnostics.R` | Add `ess()` (bulk+tail, rank-normalized; autocovariance via base `stats::fft`), `mcse()`, an internal rank-normalized R̂ (or upgrade `rhat()` with a `rank = TRUE` arg, default preserving current behavior), `.hb_diagnostic_table()` builder, and `traceplot()` generic + `traceplot.choicer_hb`. |
| `R/hmnlogit_utils.R`, `R/hmnprobit_utils.R` | **Retain all chains' hierarchical draws** (`b, w_vech, delta, theta, sigma_d2, (sigma2), loglik`) in the fit object — cheap: `delta` is R_keep×J×8B ≈ 8 MB/chain at target, the rest are tiny. Keep `beta_i` summaries for chain 1 only (do **not** replicate the big cube per chain). Feed `.hb_diagnostic_table()`; broaden the fit-time warning. |
| `R/methods.R` | `summary.choicer_hb` computes/stores the diagnostic table; `print.summary.choicer_hb` renders it and the acceptance/conjugate line. |
| `R/classes.R` | `new_choicer_hmnl` / `new_choicer_hmnp` carry the retained per-chain draws (additive field, backward compatible). |
| `NAMESPACE` (via `devtools::document()`) | New exports: **`ess`**, **`mcse`**, **`traceplot`** (+ registered `traceplot.choicer_hb`). Keep `rhat`, `ppc_shares`. |
| `tests/testthat/` | New test file: ESS/MCSE against known-autocorrelation series (AR(1) closed form), R̂ rank-normalized sanity, multi-chain shape, `traceplot` returns invisibly / draws without error, HMNP acceptance = "conjugate". |
| Docs (scriber) | `?ess`, `?mcse`, `?traceplot`; a short "Convergence diagnostics" section in the HMNL/HMNP vignettes. |

**No new package dependency.** ESS/MCSE/rank-normalized-R̂ are ~50 lines each on
base R (`stats::fft`, `stats::qnorm`); traceplots use base graphics. This matches
the package's deliberate CRAN-lean, dependency-minimal posture.

## 2.5 Explicit exclusions (do not build)

- **All HMC/NUTS-only machinery: divergences, energy / E-BFMI plots, tree-depth,
  step-size adaptation diagnostics.** Irrelevant to a Gibbs / RW-Metropolis
  sampler — there is no Hamiltonian trajectory to diagnose. Stated plainly so a
  reviewer expecting `bayesplot`-style output understands the omission is
  deliberate, not an oversight.
- No autocorrelation-function *plots*, no posterior density/pairs panels, no
  `posterior` / `bayesplot` / `coda` dependency. ESS already summarizes
  autocorrelation numerically; the four items in §2.2 are the whole set.
- No re-derivation of `ppc_shares` (fit quality, already present and orthogonal to
  convergence).

## 2.6 Correctness / reproducibility note

Diagnostics are read-only post-processing over stored draws — they do **not**
touch the sampler, so no reproducibility-contract concern arises. The only sampler-
adjacent change is *retaining* the extra chains that are presently discarded; this
is additive to the return object and leaves every existing draw and summary
bitwise identical. Validation: on a fixed multi-chain config, assert that chain 1's
stored draws and all existing summaries are unchanged vs. the current code, and
that the new `ess`/`mcse` match closed-form values on synthetic AR(1) input.
