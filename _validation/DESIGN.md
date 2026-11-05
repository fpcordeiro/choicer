# Validation suite: design and status

`README.md` says what the suite is for and how to run it. This file holds what a reader needs to
check it: each claim with its condition, derivation and references, the metrics and gate rule as
coded, the designs' calibration and identification, what the build measured, then the build's
status, the user's decisions and the findings about the package. The approved plan is
`~/.claude/plans/plan-a-revamp-of-ancient-platypus.md`; the code is authoritative where the two
differ.

**The long tier has not run.** Every number below comes from the build's probes or the quick tier's
calibration, with the R it was measured at, and is evidence for the design, not a result of the
suite.

## The claims

Notation. Likelihood units u = 1, ..., N (choice situations; decision makers in a panel), weights
w_u, ℓ(θ) = Σ w_u log L_u(θ), score s_u = ∇ log L_u, H = −E ∇² log L_u, B = E s_u s_u', both at θ0
unless stated. J alternatives, n_j ≈ N s_j choices of alternative j, S draws per unit.

### C1: first order, correct specification

*Claim.* √N(θ̂ − θ0) converges to N(0, H⁻¹), and H = B.

*Condition.* Independent units, a correct model, θ0 identified and interior, and a number of
parameters fixed as N grows (Newey & McFadden 1994, Thm 3.3; White 1982). The last condition is what
C4 and C6 qualify: simulation adds N-dependent terms, and J constants are not a fixed number once J
grows with N.

*Derivation.* A mean-value expansion of the first-order condition gives 0 = N^-½ Σ s_u(θ0) − H_N(θ̃)
√N(θ̂ − θ0). The central limit theorem applies to the scores (variance B), a uniform law of large
numbers to H_N (limit H). Differentiating ∫ L_u dy = 1 twice gives E ∇² log L_u + E s_u s_u' = 0, so
H = B.

*Consequences.* The Hessian, BHHH and sandwich variances estimate the same matrix; Wald intervals
reach their level; W = (θ̂ − θ0)' V̂⁻¹ (θ̂ − θ0) is asymptotically χ²_p. For q interior restrictions
LR is asymptotically χ²_q (Wilks 1938). In choicer's NL, λ lives on (0, ∞) with no upper bound, so λ
= 1 is interior and the LR test of IIA is regular with df = 1. For smooth g, √N(g(θ̂) − g(θ0))
converges to N(0, ∇g' H⁻¹ ∇g), which covers WTP and Σ = LL'.

*Tested.* Gated at the top rung of F1, F2, F3, F4's panel, and in F7's MNL, where C4 and C6 leave C1
standing: F2 and F3 at 16k have J/√N = 0.32, F3 has √N/S = 0.24, F4 at 3.2k has √U/S = 0.11. Q2
checks one large fit per family, Q1.1 and Q1.2 the estimator against other packages. The gated Wald,
WTP and F1's LR rows are C1's consequences at the same rungs.

### C2: sampling designs

*Choice-based sampling, MNL with a full set of constants, unweighted.* Let stratum j be sampled at
rate f_j ∝ H_j/Q_j, with H the design's shares and Q the population's. Then the sampled probability
given the covariates is P_s(j | x) = P(j | x) f_j / Σ_k P(k | x) f_k. With P(j | x) ∝ exp(x_j'β +
δ_j) this is an MNL with the same β and constants δ_j + log(H_j/Q_j), so the constants shift exactly
by log(H_j/Q_j) − log(H_ref/Q_ref) and the slopes are untouched (Manski & Lerman 1977; Prentice &
Pyke 1979; Hsieh, Manski & McFadden 1985). With random stratum sizes the Hessian variance is the
sampling variance for every coefficient. Under fixed quotas it remains so for the slopes only
(Prentice & Pyke 1979).

*WESML.* With weights w_j = Q_j/H_j, E_s[w log P(y | x; θ)] = E[log P(y | x; θ)], maximized at θ0,
so the estimator is consistent in any correctly specified model (Manski & Lerman 1977). Its bread is the weighted
Hessian. The meat depends on the sampling. Under Bernoulli sampling of a frame the kept units with
their weights are i.i.d., and the plain robust meat Σ w_u² s_u s_u' is exact. Under fixed quotas the
score's variance is Σ_j H_j w_j² Var(s | y = j), which subtracts Σ_j H_j w_j² m_j m_j' (m_j = E[s |
y = j]) from the plain meat (Manski & McFadden 1981). choicer's meat (`R/sampling.R`) is the plain
one, so under quotas it is conservative: an SE ratio above 1.

*A panel fitted as a cross-section.* Each task's marginal probability is the mixed logit, so Σ_t log
P_ut(θ) has a score with mean 0 at θ0: a composite likelihood (Varin, Reid & Firth 2011; MXL article
§1.6). The estimator is consistent, with variance H_c⁻¹ B_c H_c⁻¹, where B_c holds the within-person
covariance of the tasks' scores, so H_c ≠ B_c when tastes are shared. Cluster-robust SEs are valid
and Hessian SEs too small, and it is less efficient than the panel likelihood.

*Tested.* F8: unweighted MNL under random sizes (gated: bias against the exact shift, Hessian SE
ratio, tails, joint Wald), WESML robust tails under random sizes for the MNL and NL populations
(gated), the SE ratio above 1 under quotas (a finding; the population predicts 1.008 to 1.018 for
the constants, 1.000 for the slopes), and the unweighted NL on the NL population, which misses the
shift (a finding; its pseudo-true value is exact). F4's `xsec`: bias, cluster SE ratio and tails
(gated at 800 persons), the Hessian SE ratio below 1 and the efficiency loss (findings).

### C3: misspecification

*Claim.* θ̂ converges to θ* = argmax E log L_u(θ), and √N(θ̂ − θ*) converges to N(0, H*⁻¹ B* H*⁻¹)
with H* ≠ B* in general (White 1982). Only the sandwich estimates this variance, so this is the one
design that separates the robust route from the Hessian.

*Condition.* θ* unique and interior; independent units.

*Tested.* F11 fits an MNL to D2's mixed logit population. θ* is exact, not simulated: the population
is 210 traveller types, so E log L_u is a finite sum (`D2_THETA_STAR`, checked in
`checks/wp3_mode.R` against a fit of 10^7 travellers). Gated at 8k: bias against θ*, robust SE ratio
and tails, the robust joint Wald. A finding: the Hessian SE ratio, whose population values at θ* are
0.79 (travel), 0.84 (income × air), 0.85 to 0.88 (constants), 0.97 (wait), 0.99 (vcost).

### C4: simulation

*Claim.* With P̄_S the simulated probability, a Taylor expansion gives E log P̄_S = log P −
Var(P̄_S)/(2P²) + O(S⁻²). With i.i.d. draws Var(P̄_S) = σ²/S, so the expected simulated
log-likelihood is displaced by O(N/S), its score by the same, and θ̂_S − θ̂_ML has a bias b/S.
Against the sampling SD, O(N^-½), that is bias/SD = O(√N/S). MSL is equivalent to ML only if √N/S
goes to 0, and with S ∝ √N a bias remains (Lee 1992, 1995; Gouriéroux & Monfort 1996; Train 2009
§10.5). The simulation noise in the score is independent across units with i.i.d. draws, so the
simulation SD over the sampling SD is O(S^-½).

*What is not proved.* These are i.i.d. rates. choicer's permuted Halton draws have no randomized-QMC
guarantee (MXL article §2.5, §5.4), so the error is measured, and common choices of S are often too
low (Czajkowski & Budziński 2019). A panel's integrand, a product over a person's tasks, has a
larger σ²/P², which raises the constant.

*What it means for the variance routes.* The routes estimate the sampling variance given the draws.
If r is the simulation SD over the sampling SE, then ESE is inflated by √(1 + r²), and the SE ratio
reads 1/√(1 + r²), below 1, with nothing wrong in the standard errors.

*Tested.* The gates run where √N/S ≤ 0.25. F5 asks whether that is the right condition: on the gated
regimes' own data, it measures bias/SE and SD/SE against S, warm-started from a reference at S =
3,200 (800 at population scale), with their rate in S, the S at which |bias| reaches 0.1 SE, and
that S against 4√U. Findings only. Q1.3 compares MSL at S = 4,000 with exact Gauss-Hermite ML.

### C5: higher order, regular models

*Claim.* E θ̂ − θ0 = b/N + O(N⁻²) (Cox & Snell 1968; for the MNL, Kosmidis & Firth 2011), so N·bias
is stable net of simulation bias, and bias/SD = O(N^-½). The Edgeworth expansion of the t statistic
has an even leading term (Hall 1992 §3.5), so a one-sided tail errs by O(N^-½) and the symmetric
two-sided interval by O(N⁻¹): the tails' imbalance falls as N^-½. E[LR] = q(1 + c/N) (Bartlett 1937;
Lawley 1956). Ratios such as WTP are skewed when |α|/SE(α) is small (Daly, Hess & de Jong 2012).

*Consequences for the design.* F1 gates the value of travel time's two-sided coverage, not its
tails: |b_vcost|/SE is 2.0, 4.0, 7.9 at 500, 2k, 8k, and under the normal approximation of the
coefficients the tails at 8k are 5.1% and 0.1% (coverage 94.8%), at 500 10.9% and 0.0%.

*Tested.* Findings on the rungs of F1 to F4: N·bias, √N·ESE, tail imbalance, mean LR/q, with fitted
slopes in the report. N·bias needs F5's simulation bias subtracted in the mixed logit. The quick
tier's dropped rows (comments beside each `gates_quick`) are measurements of this claim at the
bottom rungs.

### C6: many alternatives

*The constants.* Var(δ̂_j) and the bias of δ̂_j are O(1/n_j).

*The common parameters: a derived prediction.* Let θ be the common parameters and δ̂(θ) the
constants given them. Expanding the score for θ in δ̂ − δ0 to second order, the expected score at
the truth is Σ_j [ℓ_θδj E(δ̂_j − δ_j) + ½ ℓ_θδjδj Var(δ̂_j)] plus smaller terms. Each derivative is
of order n_j, and the moments are O(1/n_j), so each alternative contributes O(1), and the sum O(J).
The information is O(N). So bias(θ̂) = O(J/N), SD(θ̂) = O(N^-½) and bias/SD = O(J/√N). This is the
Hahn & Newey (2004) expansion with alternatives as the units and n_j as the observations per unit.
Centered intervals need J²/N to go to 0, and the expansion assumes min_j n_j goes to infinity. The
plan's review found no published theorem for the mixed logit with constants, which is why this is
measured. The claims and census regimes have J/√N of 2.0 (admissions; 2.5 on patients) and 3.0.

*The MNL's slopes.* The MNL is a Poisson model with situation and alternative effects (Baker 1994).
In the two-way Poisson model the O(1) terms above vanish: the slopes have no leading-order
incidental bias (Fernández-Val & Weidner 2016; Weidner & Zylkin 2021).

*What stays biased in the MNL.* For g(δ̂) non-linear, E g(δ̂) − g(δ) ≈ ½ Σ_j g_jj Var(δ̂_j), of
order J/N in functionals that sum over alternatives (closure diversion, the change in consumer
surplus). Sandwich and cluster variances are biased down in finite samples (Pfaffermayr 2019),
through the constants' leverage.

*Tested.* F6 on J = 200 with N/J = 25, 100, 400 (J/√N = 2.83, 1.41, 0.71). Gated: the MNL's slopes,
bias and Hessian SE ratio, at every rung. Findings: the MNL's robust SE ratio below 1 at small N/J;
the MXL's N·bias/J and bias/ESE against J/√N; whether `"blp"` reduces it; pooled constant coverage
by expected-count tercile; the closure's bias. F7 repeats this at population scale, with the MNL
gated as F1.

### C7: Bayesian samplers

*Exact at any N.* Let θ̃ be drawn from π, y from p(y | θ̃), and θ̃_1, ..., θ̃_L from p(θ | y). By
Bayes' rule θ̃ given y has the posterior distribution, so θ̃ and the draws are exchangeable given y,
and the rank of f(θ̃) among the f(θ̃_l) is uniform on {0, ..., L} for any f, normalized MNP and HMNP
quantities included (Cook, Gelman & Rubin 2006; Talts et al. 2018). The conditions: the prior that
simulates is the prior that fits; the draws are nearly independent (F9's chains were lengthened until the
measured bulk ESS was above L = 99, and a rep below it is flagged, not dropped); the chain's law is the
posterior. HMNL starts from a pooled MNL fit,
but its adaptation freezes at burn-in, so this holds.

*Bernstein-von Mises* (van der Vaart 1998, Thm 10.1). The posterior approaches N(θ̂, H⁻¹/N) in total
variation, for quantities whose information grows with N: MNP's β and Σ; HMNL's and HMNP's b, W and
δ at fixed J. Credible intervals then attain frequentist coverage and posterior SD/ESE goes to 1. θ
and σ_d² gain information only with J, so for them SBC alone applies.

*Tested.* F9: SBC on every scalar quantity and the log-likelihood (Modrák et al. 2023), gated by the
ECDF test (Säilynoja, Bürkner & Vehtari 2022) with the exact null; BvM coverage and SD/ESE as
findings, with chains lengthened until the smallest tail ESS is about 400, where a credible
interval's tail reads 2.60% in place of 2.5%. Rates of the fit-time R-hat and ESS warning are
reported. Q1.4 compares MNP with `bayesm`, Q4 runs a mini SBC.

### References

Baker (1994), *The Statistician* 43. Bartlett (1937), *Proc. R. Soc. A* 160. Cook, Gelman & Rubin
(2006), *JCGS* 15. Cox & Snell (1968), *JRSS B* 30. Czajkowski & Budziński (2019), *J. Choice
Modelling* 31. Daly, Hess & de Jong (2012), *Transp. Res. B* 46. Fernández-Val & Weidner (2016), *J.
Econometrics* 192. Gouriéroux & Monfort (1996), *Simulation-Based Econometric Methods*. Hahn & Newey
(2004), *Econometrica* 72. Hall (1992), *The Bootstrap and Edgeworth Expansion*. Hsieh, Manski &
McFadden (1985), *JASA* 80. Kosmidis & Firth (2011), *Biometrika* 98. Lawley (1956), *Biometrika*
43. Lee (1992), *Econometric Theory* 8. Lee (1995), *Econometric Theory* 11. Manski & Lerman (1977),
*Econometrica* 45. Manski & McFadden (1981), in *Structural Analysis of Discrete Data*. Modrák et
al. (2023), *Bayesian Analysis*. Morris, White & Crowther (2019), *Statistics in Medicine* 38. Newey
& McFadden (1994), *Handbook of Econometrics* 4. Pfaffermayr (2019), *Applied Economics Letters* 26.
Prentice & Pyke (1979), *Biometrika* 66. Revelt & Train (1998), *REStat* 80. Säilynoja, Bürkner &
Vehtari (2022), *Statistics and Computing* 32. Self & Liang (1987), *JASA* 82. Talts et al. (2018),
arXiv:1804.06788. Train (2009), *Discrete Choice Methods with Simulation*, 2nd ed. van der Vaart
(1998), *Asymptotic Statistics*. Varin, Reid & Firth (2011), *Statistica Sinica* 21. Weidner &
Zylkin (2021), *J. International Economics* 132. White (1982), *Econometrica* 50. Wilks (1938),
*Ann. Math. Statist.* 9.

## Metrics and the gate rule, as coded

`lib/metrics.R`. A cell's metrics read its OK reps; a failed rep enters `failure_rate` only. R below
is the count of reps entering the row: those with a finite estimate and, for a route's rows, a
finite SE. d_r = θ̂_r − θ0_r (each rep's own truth), ESE = sd(d), ŝe_r the route's SE, z =
Φ⁻¹(0.975).

| Metric | Value | MCSE | p-value for the gate |
|---|---|---|---|
| `bias_std` | mean(d)/ESE | 1/√R | exact under normal errors: √R·`bias_std` against t with R − 1 df |
| `ese` | sd(d) | ESE/√(2(R − 1)) | never gated |
| `se_ratio` | M/ESE with M² = mean(ŝe²) | (M/ESE) √(Var(ŝe²)/(4 R M⁴) + 1/(2(R − 1))) | normal: 2Φ(−\|value − 1\|/MCSE) |
| `tail_lower`, `tail_upper` | mean(d > z ŝe), mean(d < −z ŝe) | smoothed proportion | exact binomial |
| `coverage` | 1 − the tails | smoothed proportion | exact binomial |
| `joint_wald_size`, `lr_size` | rate of p < 0.05 over finite statistics | smoothed proportion | exact binomial |
| `joint_wald_failure_rate`, `route_failure_rate` | share of reps with a statistic or SE not finite | smoothed proportion | exact binomial, target 0 |
| `lr_mean_q` | mean(LR/q) | sd(LR/q)/√R | never gated |
| `failure_rate` | failed reps over all reps | smoothed proportion | never gated |
| `global_max_rate` | share of multi-start reps with ΔlogL > 0.01 and max \|Δθ\|/SE > 0.1 | smoothed proportion | never gated |
| `sbc_ecdf` | the smallest two-sided binomial tail over the evaluation points | none | exact (below) |

The *smoothed proportion*: for x of R, p̃ = (x + 1)/(R + 2) and MCSE = √(p̃(1 − p̃)/R), which is
Morris, White & Crowther's (2019, Table 6) value at large R and never 0 for a count of 0 or R. The
gate does not read it. The *exact binomial* p is 2 min(P(X ≤ x), P(X ≥ x)) capped at 1, X ~
Binomial(R, target). The exact references matter at small R: with the normal one 34% of null
families fail at R = 2, and about 1e-3 per tail row at R = 100.

The *joint Wald* statistic is computed on the Jacobi-equilibrated V̂, on the common parameters.
*Pooled coverage* (`pooled_tails()`) averages each rep's share of constants outside the interval, by
expected-count tercile, with MCSE the between-rep SD over √R, since a rep's constants share its
data. *Ladders* (report): ESE·√N, N·bias (MCSE N·ESE/√R) and the tail imbalance upper minus lower
(its MCSE from the multinomial). `rate_fit()` is weighted least squares of log|value| on log x with
weights (|value|/MCSE)², its slope's SE the known-variance one.

*SBC.* Ranks among L = 99 equally spaced draws, ties split at random. At z_i = i/(L + 1), #{rank <
i} is Binomial(R, z_i) under uniformity. The statistic is the smallest two-sided binomial tail over
i, and p = P(min_i tail_i ≤ stat) under the exact discrete null: the rank counts are multinomial, so
the cumulative count is a Poisson walk conditioned on its sum, carried over each point's acceptance
interval (`lib/sbc.R`), with no Monte Carlo floor. The band in the plots is this test's: the ECDF
leaves it exactly when p ≤ α.

**The gate.** `gate_family(tab, alpha = 0.01, rule = "plan")`: p_holm is the Holm adjustment over a
family's gated rows with a p-value. A gated row is:

- `na` where its value, MCSE, target or tolerance is missing (SBC: its p);
- `fail` where p_holm < 0.01 and |value − target| > tol (1 + 1e-9), the factor keeping rounding from
  deciding a verdict (SBC: p_holm < 0.01 alone);
- `pass` otherwise.

A row not gated is `finding`. Tolerances (`tolerances()`):

| Set | `bias_std` | `se_ratio` | tail | coverage | size | failure rate |
|---|---|---|---|---|---|---|
| mc | 0.10 | 0.05 | 0.0075 | 0.015 | 0.015 | 0.005 |
| quick | 0.30 | 0.15 | 0.025 | 0.05 | 0.025 | 0.02 |

Significance alone would fail harmless deviations at R = 5,000; a tolerance alone would fail noise
at R = 100. Under this rule a true value inside the band but off target fails at a rate that grows
with R. `rule = "band"` tests instead that the true value lies inside the band, one-sided against
the nearer edge, with the tests above at that edge (a central t approximation for `bias_std`), and
fails where p_holm < 0.01, so a true value anywhere in the band fails a family with probability at
most 0.01. The report shows it in the column `band`; S1 kept the plan's rule for the verdict.
`checks/wp1_metrics_null.R` measures both: an exact null family of 33 gated rows fails in at most 2%
of runs, and the gate's power at the tolerances equals its closed form.

A gate entry is `metric[:route][@block or estimand]`, and a typo in one stops the run. An entry that
matches no row adds a row with value NA, so that a gated check that cannot be computed is `na`,
never `pass`. Gating an SE metric also gates its route's failure rate, and gating a joint Wald size
its failure rate, so that reps dropped by a failed route count against the gate. Gates sit only
where the theory's conditions hold at the design (each family's header says why); the rest are
findings.

## The designs: calibration and identification

The DGPs never call choicer's kernels: probabilities come from `lib/simulate.R`, an independent
implementation in choicer's parameterizations (MXL log-normal β = exp(μ) + exp((Lη)); NL RU2 with
singleton nests at λ = 1), and Q0 checks the two against each other at θ0, within 1e-12 for
probabilities and 1e-12 relative for simulated log-likelihoods. θ0 is frozen as literals in the
registry. D2's, and the pseudo-true values of F8 and F11, came from fits through choicer, once, and
are checked in `checks/wp3_mode.R`; as literals they do not depend on the build under test. D3 and D4 share one generator: `lib/geo_gen.cpp` gives data identical to the R
reference (`checks/wp6_geo.R`). No DGP writes outside rows. `lib/calibrate.R` prints each DGP's
moments against its targets, which are stated ranges, not fits to the cited papers.

| DGP | θ0 | Targets, and the moments at θ0 (N in the last column) | N |
|---|---|---|---|
| D1 | literals; the constants' level calibrated to the opt-out share | opt-out 0.10-0.20: 0.150; SD/mean of the time coefficient 0.5-1: 0.571; of quality: 0.667; WTP for quality 1.4-1.6: 1.50 | 10^5 persons × 8 |
| D2 | the WESML MNL fit of `mode_choice` at Q = (air .14, bus .09, car .64, train .13); NL likewise; MXL below | population shares within 1e-8 of Q (a full set of constants makes this exact at the WESML first-order conditions) | the 210 types, exact |
| D3 | literals; the distance location calibrated to the nearest-hospital share | nearest 0.30-0.55: 0.398; distance elasticity −3 to −1: −1.97; J Σ s_j² 2-6: 2.21 | 2 × 10^5 |
| D3-claims | as D3 | 0.398; −1.99; 2.31; admitted once 0.65-0.85: 0.750; five times or more 0.01-0.06: 0.040; mean admissions 1.3-1.8: 1.60 | 2 × 10^5 patients |
| D4-40 | literals; sector constants calibrated per variant to outside .15, public .65, charter .09, private .11 | MXL: outside 0.10-0.20: 0.147; public 0.55-0.75: 0.652; charter 0.05-0.15: 0.091; private 0.07-0.15: 0.110; quality gradient 0.15-0.40: 0.304 (NL: 0.150; 0.650; 0.090; 0.110; 0.381) | 2 × 10^5 |
| D4-census | as D4-40 | 0.148; 0.654; 0.091; 0.108; 0.280 | 2 × 10^5 |

**D1.** The attributes are level-balanced within person and position. The tastes' parameters are
literals: b_p = −1, μ_q = 1.5 (mean WTP for quality 1.5 price units), sd(g_2) = 1, exp(μ_t) = 0.5,
sd(g_1) = 0.7, ρ = 0.4. Price is fixed because `wtp()` refuses a random price. Measured σ/SE(σ) at
200 persons, median of five panel fits at S = 500: 5.96 (sd g_1), 10.97 (sd g_2), 3.37 (the
correlation of g), and at 800 persons, one fit: 14.3, 20.0, 3.4. Position constants are equal in
truth and estimated separately.

**D2.** The population resamples traveller i with probability Q_{y_i}/n_{y_i}, undoing the sample's
choice-based design. The free WESML MXL fit puts σ at 0 (σ/SE = 0.05), so that variant holds σ at
|μ_travel|, a chosen coefficient of variation of 1, with the constants re-solved so that its shares
equal Q. It serves F11, where it must make the MNL wrong, not as an estimate of the data. Its
integral uses Gauss-Hermite with 2,000 nodes (1e-12; 1,000 give 2e-9).

**D3.** Sets cross every boundary, connecting the constants (`checks/wp6_geo.R`: one component).
Both random columns vary by patient within a hospital, so `rc_mean = TRUE` is identified beside the
constants. Distance is in units of 10 km, because choicer's shifted log-normal has a median of at
least 1 per unit of its column. The targets are plausible ranges for US inpatient choice, not
estimates.

**D4.** With a constant for every school, a mean on quality, a school attribute, is collinear with
the constants. So `rc_mean = FALSE`: quality's normal term is mean-zero, its mean effect is in the
constants; distance's coefficient is the fixed b_d minus a log-normal, a shifted log-normal with
free location. σ_q is identified by how substitution among schools of different quality varies with
the composition of the nearest-K sets, which is why D4-40 is one market of 40 schools with K = 15,
not four of 10 (with 10, K = 15 would offer everyone the whole market). The plan's target
σ_q/SE(σ_q) ≥ 5 at N = 1,000 is met after raising s_d, s_q and ρ_q: 6.7 and 6.6 on D4-40 (a 16k
fit's SE times 4, and the median of five 1k fits), 5.5 and 5.1 on D4-small (`checks/wp5_geo.R
--only=ident`).

**Population scale.** Expected counts per constant, minimum, 10th percentile and median: claims 1.3,
41, 308; census 11.8, 65, 219. So min_j n_j is small in the claims even though the median is not,
and about 0.9 hospitals are expected never to be chosen there.

**Alternatives never chosen.** ML's limit as δ_j goes to −∞ is the fit without j, so
`drop_unbounded_alts()` drops such alternatives, and situations offering one chosen wherever
offered, until stable; `"blp"` requires it. Counts go to `diag$drop` (F6 expects about 10 of 200 at
5k, 1.5 at 20k, 0.1 at 80k).

**The Bayesian designs.** Priors are proper and moderately informative, passed explicitly
(`lib/dgp_bayes.R`, whose header holds the prior-predictive checks behind them). The MNP's base is
`car`, not choicer's default `air`: the normalization would otherwise fall on bus − air, bus being
the least chosen mode, and the identified draws mixed 2 to 10 times slower. D1-L's J = 6 is 3P with
one characteristic, the fewest that avoid `run_hmnlogit()`'s few-alternatives warning. F9's chains
are 5,000 + 30,000 kept every third for SBC, since 10,000 straight draws left the smallest bulk ESS
at 73 to 173.

## What the build measured

**C6 on D4-small (J = 160).** The distance coefficient's bias/SD in the mixed logit at R = 100, the
same from choicer's default start and from θ0: −0.74 ± 0.10 at N = 1,000 (J/√N = 5.06) and −0.41 ±
0.10 at 4,000 (2.53), about −0.16 J/√N, the rate C6 predicts. With the constants held at their
truth, at 4,000: −0.02 ± 0.10. So the bias is the constants', and a gate at 16k (predicted −0.20)
would test C6, not C1: S1 moved F2 and F3 to J = 40. On D4-40 at J/√N = 1.26, the rule predicts
−0.20, and the mixed logit's measure is −0.175 ± 0.058 (R = 300); the nested logit's, −0.111 ± 0.018
(R = 3,000), is smaller. `"blp"` with the sample's shares on D4-small at 4,000 (R = 60, paired with
ML) moved the coefficient by 0.103 ± 0.004 ML SE toward the truth. These are probes; F6 and F7 are
the measurement.

**The default start's plateaus.** With an L diagonal going to −∞, the score in log σ vanishes, since
it is σ times the score in σ. The polish then certifies a point on a plateau: SE about 10^4, on D1
up to 14 log-likelihood units below the maximum (package finding 3). From choicer's default start,
fits per cell, and those from θ0 in the last column:

| Fit (D1; panel S = 500, composite S = 200) | Reps | Not certified | On a plateau | Certified, below θ0's by over 0.01 | Their median and max shortfall | From θ0: not certified, plateau |
|---|---|---|---|---|---|---|
| panel, 200 persons | 30 | 0 | 2 | 2 | 3.84, 5.08 | 0, 0 |
| panel, 800 | 16 | 0 | 0 | 0 | | 0, 0 |
| panel, 3.2k | 4 | 0 | 0 | 0 | | 0, 0 |
| composite, 200 | 30 | 7 | 15 | 17 | 1.23, 4.52 | 0, 0 |
| composite, 800 | 16 | 11 | 6 | 5 | 1.57, 4.03 | 0, 0 |
| composite, 3.2k | 4 | 2 | 2 | 2 | 12.95, 13.63 | 0, 0 |

In one panel fit and one composite fit at 200 persons the default start ended above θ0's by over
0.01: θ0 is not always the better start, which is why both run. On D4-small at N = 1,000, S = 500, 2
of 100 default-start fits ended on the plateau (L_11 near −18 and −20), 7.03 and 0.36 log-likelihood
units below θ0's, and none at 16k (4 fits). Two starts do not remove the plateau where the sample's
maximum sits there: on F4's quick cell (200 persons, S = 200, R = 300) the composite's default start
lost in 239 fits and 20 kept fits still ended on a plateau; on F3's, 1 of 300, whose SE alone takes
L_11's SE ratio to 716. That is one reason the quick tier does not gate L_11 at those rungs.

**F4's simulation variance.** On the panel at 3.2k persons, 12 datasets, against a reference at S =
3,200 (each dataset's Hessian SE the unit): the bias at S = 500 is at most about 0.1 SE: L_11's
−0.100 ± 0.091 in the ladder's fit, L_21's −0.109 ± 0.076 in F4's own, every other parameter below
0.07. The SD is not negligible:

| Parameter | SD/SE at S = 500 | at 1,000 | at 2,000 | S where \|bias\| reaches 0.1 SE (rate fit) |
|---|---|---|---|---|
| price | 0.029 | 0.028 | 0.020 | 146 |
| irr | 0.012 | 0.011 | 0.005 | 1 |
| μ, −time | 0.023 | 0.022 | 0.018 | 25 |
| μ, quality | 0.066 | 0.058 | 0.032 | 188 |
| L_11 | 0.31 ± 0.07 | 0.23 | 0.23 | 114 |
| L_21 | 0.16 ± 0.04 | 0.17 | 0.13 | 463 ± 260 |
| L_22 | 0.074 | 0.051 | 0.046 | 20 |

These SDs are of the difference from the reference, which carries the reference's own error, so they
fall more slowly than S^-½ as S nears 3,200. By C4, L_11's r = 0.31 ± 0.07 would read as an SE ratio
of 0.955 at S = 500, inside the tolerance of 0.05 but not by much: if F4's gate on it fails at the
top rung, look here before suspecting the standard errors. The pair of draw seeds at S = 500 gives 0.28 ±
0.06 and 0.23 ± 0.05 for L_11 and L_21. And the rule of thumb 4√U = 226 would have allowed a lower S
than the rate fit's largest value, 463 for L_21, which is what C4 says of panels; S = 500 was kept.
F5's 30 datasets replace this probe.

## Status

| WP | Group | Owns | State | Commit | Notes |
|---|---|---|---|---|---|
| WP0a | A | `lib/common.R`, `lib/families.R`, `lib/fit.R`, `lib/calibrate.R`, `lib/families/T0.R`, `checks/wp0_common.R` | done | 34db843 | build guard, seeds (FNV-1a + rep, method-free), atomic rep files, metadata, contracts as code; gates per method and over several cells; T0 is a plumbing fixture; 38/38 checks |
| WP0b | A | `lib/simulate.R`, `checks/wp0_simulate.R` | done | 34db843 | the suite's probability engines in choicer's parameterizations, equal to `predict()` to ≤ 7e-16 and to the fixed-parameter MXL likelihood to ≤ 3e-16 rel.; 117/117 checks; shared file, lead-serialized afterwards |
| WP1 | B | `lib/metrics.R`, `lib/sbc.R` (test half); `checks/wp1_metrics_null.R` (test-writer) | done | 9d71868 | null study 85 PASS / 0 FAIL / 11 INFO (four bugs it found fixed: an MCSE copied across rows, missing na rows, SBC mirror ties, infinite W counted as a rejection); reviewer: SBC p floor → exact null by DP, route-failure rates, global-max criterion; `gate_family(rule = "plan" | "band")`, plan by default, for the user at S1 |
| WP2 | B | `mc.R`, `report.R`, `run_shards.sh`, `quick.R` (skeleton), `checks/wp2_plumbing.R` | done | 9d71868 | plumbing 16/16; `run_rep()` and `vs_boot()` moved into the lib; `R_by_method`; reviewer: absent units and missing reps flagged, `--tier`, launcher checks; plumbing 20/20 |
| WP3 | C | `lib/dgp_mode.R`, `lib/families/F1.R`, `F8.R`, `F11.R`, `lib/references.R`, `checks/wp3_mode.R` | done | 8ff4d74, 83b6f74 | 25/25 checks; Q1.1 7e-10 SE, Q1.2 6e-8 SE; θ* exact (210 traveller types), checked against an N = 10^7 fit; `require_converged()` |
| WP4 | C | `lib/dgp_sp.R`, `lib/families/F4.R`, `checks/wp4_sp.R` | done | 8ff4d74, 83b6f74 | 29/29 checks; D1 and D1-L (HB, logit and probit); `sigma_vech()`; start-policy study |
| WP5 | C | `lib/dgp_geo.R`, `lib/families/F2.R`, `F3.R`, `lib/references_quad.R`, `checks/wp5_geo.R` | done | 8ff4d74, 83b6f74 | the quadrature reference gets its own file so WP3 and WP5 own disjoint files |
| S1 | — | | decided 2026-11-02 | | see "Decisions" below; applied in the next commits |
| WP6 | D | `lib/geo_gen.cpp`, `lib/families/F6.R`, `F7.R`, `lib/dgp_geo.R` (from WP5), `checks/wp6_geo.R` | done | c60451e | C++ port equals the R reference (39/39 and 300 lattice cases); 187/187 checks; contraction self-test and compile lock; F7 ≈ 340 CPU-h after S2's Newton route, plan 650 (WP6b) |
| WP7 | D | `lib/families/F5.R`, `F10.R`, `checks/wp7_sim.R` | done | c60451e | F10 reproduces the old share-matching script rep by rep (store cube kept); F5 warm-started Newton ladders on F3, F4, F6, F7; `findings` tables; two-start tie rule (tol 1e-6); `"blp"` fits certified by their decrement; F4's gate S = 500 has simulation bias < 0.1 SE (12-dataset probe) |
| WP8 | D | `lib/dgp_bayes.R`, `lib/sbc.R` (loop half), `lib/families/F9.R`, `checks/wp8_bayes.R` | done | c60451e | 55/55 checks; Q1.4 vs bayesm 2.2 MCSE / SD 2.8%; Q4 passes; chains 5k + 30k thinned by 3 (10k left ESS below 99); MNP base `car`; BvM tail ESS recorded, chains lengthened (F9 ≈ 28 CPU-h) |
| S2 | — | | decided 2026-11-05 | | R = 2 smoke of every family: 142/142 ok, report renders, 0 gated failures; see "Decisions" |
| WP9 | E | `quick.R` (assembly), `lib/quick_checks.R`, `gates_quick` in F1–F4, F8, `checks/wp9_quick.R` | done | (this commit) | every gated row passes at R = 300 and every piece passes on its own; 6.5 to 9.5 min estimated on 8 workers against 15 |
| WP6b | E | `fit_newton()`, F7's route, F5's F7 regimes | done | (this commit) | the S2 decision, measured at N/25 before freezing |
| WP10 | E | `README.md`, `DESIGN.md`, article and NEWS lines, deletions, `.claude/pending-validation/` (AGENTS.md patch and the cloud files' removal, run by the user) | done | (this commit) | the bundle was rehearsed on copies only |
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

**S2 (2026-11-05), the user's answers:**
- F7's mixed logit fits take Newton steps from θ0 (`fit_newton()`). Measured at N/25, on 10 fits:
  all certified with no fallback, the same maximum within 9e-5 SE, 8.3× (claims) and 3.6× (census)
  faster than two starts; at full size with small S, 5–11 Newton steps at 7.1 GB and 3.4 GB. F7 falls
  from about 1,060 to about 340 CPU-h (extrapolated linearly in rows and draws; no full-size fit at
  S = 100 was run, so the first mc reps' `diag$newton_steps` replace the extrapolated step counts).
  The S1 two-start rule stays everywhere else; at population scale the default start runs on every
  tenth rep only, so that its behaviour there stays measured.
- F7's MXL cells run 100 reps first, extendable with `--reps=101:1000`.
- F9's 800-person HMNL chains stay at 90k iterations thinned by 9 (one quantity's tails may run
  about 0.2 pp wide; BvM is a finding).
- The mc tier runs later on one server of about 64 cores and at least 64 GB of RAM (four concurrent
  population fits of at most 15 GB at 16 threads), sized and set up in a later session.

## Known gaps

- **The quick tier's coverage is uneven by design.** Its gate lists hold what theory supports at the
  small rungs; measured at R = 100, one gated slope moved 0.5 ESE is caught 85–95% of the time, one
  dropped slope 28–52%, SEs × 0.8 on a whole coefficient block 85–99.5% (F4's panel 50–66%). F2's λ, F3's
  L and F4's L block are not gated there. Q0's likelihood equality and the polish certificate fix the
  estimators, Q1.3 and Q2 cover most of the rest; a defect shared by the Hessian and BHHH routes of the
  panel's L block has no quick check (the mc tier's F4 has one).
- **Q3 spends α = 0.01 per family**, the plan's rule, so about 5% of quick runs would fail somewhere under
  new draws with a correct build. The seeds are frozen, so a correct build meets this only when a change
  moves the draws.
- **F4's gate carries simulation variance.** At S = 500 and 3.2k persons L_11's simulation SD is about
  0.3 SE, which the Hessian SE leaves out: its SE ratio should read about 0.955 with nothing wrong,
  inside the 0.05 tolerance, near its edge.
- **F7's global-maximum evidence is thin.** With the default start on every tenth rep, R = 100 gives 10
  comparisons (0 of 10 still allows a missed-maximum rate up to about 26%), and the census `"blp"` fits
  start once, from the truth.

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
8. From the quick tier's confirmation: at F4's small rung the composite fit's default start lost to θ0 in
   239 of 300 reps, and 20 kept fits ended on a plateau (the same family as item 3).
