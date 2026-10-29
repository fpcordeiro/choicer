# choicer (development version)

## Share inversion

`blp()` and the exported `blp_contraction()`, `nl_blp_contraction()` and
`mxl_blp_contraction()` solve the share inversion, the constants at which a
model's market shares equal target shares (Berry 1994), with one engine in
place of three copies of Berry's contraction. The solution is the same; the
iteration that reaches it, its stopping rule and its handling of inputs that
have no solution are new, so results move within the tolerance. The
contraction remains available, under the engine's rules, as
`blp(method = "original")`, and the new step accelerated as
`method = "anderson"` (Anderson acceleration) and `method = "squarem"`
(SQUAREM).

- **The step.** The iteration solves for the free constants at which each
  log-share residual against a base, r_j - r_b (r = log S - log s), is zero;
  the base is the outside option (the correction of Fukasawa 2024) or,
  without one, the reference of the alternative's group. Each iteration
  evaluates the shares once and takes a step in two parts. It first shifts
  the group's free constants together by tau, chosen to match the base's
  share to first order: residual j responds to such a shift at the rate
  -h_j, with h_j = C_jb / s_j + sum_k C_kb / s_b exactly and C_jb the
  co-choice of j with the base. It then divides the rest of each residual
  by that residual's exact derivative in its own constant, A_jj / s_j +
  C_jb / s_b, where A_jj is the share Jacobian's diagonal. Without the
  shift, that diagonal misjudges a common move of the constants: where the
  base's share is small and each alternative dominates choice situations of
  its own, it multiplies the shift's error by less than -1 and diverges
  where Berry's contraction converges. Berry's contraction moves each
  constant by its own log-share residual, r_j (r_j - r_b with the reference
  held), and remains available as `blp(method = "original")`.

  On 121 simulated designs (multinomial, nested and mixed logits; 2 to
  5,000 alternatives; full, regional and chained choice sets; outside
  shares from 0.02 to 0.5; far starts), solved to 1e-12, the step converged
  on every one and needed a median 0.36 of the contraction's evaluations of
  the shares (10th to 90th percentiles 0.04 and 0.81; medians 32 against
  161 for the multinomial logit, 15 against 184 for the nested and 174
  against 593 for the mixed logit, over the designs where both converged).
  It never needed more than 1.05 times the contraction's evaluations: the
  most was on binary choice sets chained over ten alternatives without an
  outside option, 1,357 against 1,299, both past the default `max_iter`.
  The contraction did not converge within 3,000 evaluations on 8 nested
  logits with lambda of 0.4 or less, where its step overshoots by about 1 /
  lambda; the step needed 23 to 119. On 80 designs in which each
  alternative dominates choice situations of its own and the base is rarely
  chosen, solved to 1e-10 for each model, the step did not converge within
  3,000 evaluations on 1, 2 and 1 (multinomial, nested and mixed logits;
  each was still approaching its solution, and reached it within 26,000,
  the racing contraction's turns included; the contraction alone failed on
  them too), and the contraction on 22, 27 and 20; where both converged,
  the step needed a median 0.16 of the contraction's evaluations, at most
  0.52.
- **Safeguards.** A step is at most 5 log units long, or as long as its
  residual against the base, |r_j - r_b|, where that is longer. An
  evaluation whose shares are not finite halves the step. An iteration whose
  best residual has not improved for 30 iterations (a stall) halves its later
  steps. Where its residual at the stall is more than twice the best (it
  diverges), it first returns to its best point, once per best point;
  otherwise it goes on from where it is, and each such stall doubles the
  iterations the next one waits for, so that constants climbing back from
  where their shares barely respond are not pulled back. Two more safeguards
  apply in nests with lambda above 1, whose nest-mates are complements where
  the nest's choice probability is below 1 - 1 / lambda. An alternative whose
  co-choice with the reference is negative divides by at least 1. The
  step's divisor is at least (1 + lambda) / (2 lambda) times the one the
  same probabilities give with lambda = 1, where the exact one would
  overshoot a common shift of the nest's constants and, for lambda above 2,
  diverge. On 80 nested logits with lambda from 1.05 to 5 the step converged
  on all, in a median 0.32 of the contraction's evaluations and never more
  (with an outside option the contraction needed 129 to 170 evaluations, the
  step 9 to 70). `damping`, for the nested logit, now multiplies the step and
  must lie in (0, 1]; with the default step, small lambda no longer needs it
  (`method = "original"` still may).
- **Berry's contraction as a fallback.** The step can amplify the small
  residual of an alternative that nearly owns its choice situations into a
  step that wrecks its neighbours' shares there, and it can carry the
  constants of alternatives that make up the whole choice set of some choice
  situations so far down that the residual barely responds to their common
  level; the contraction does neither. So once 100 iterations fail to halve
  the best residual, or a trial point's shares stay non-finite after every
  halving (the step then stops, and the contraction goes on alone), the
  contraction races the step from the start: the two take steps in turn,
  each with its halvings, and the first to meet `tol` is returned. The
  contraction stops if 100 of its own evaluations fail to halve its best
  residual while the step's best is lower; the step then has 200 iterations
  at a time to halve its own best residual, and whenever it fails to, the
  contraction goes on from where it stopped. Where the contraction converges
  in E evaluations, the default converges within T + 4E, and within T + 2E
  where the contraction halves its best residual within every 100 of its
  evaluations, T being the evaluations before the contraction joined; the
  bound leaves out the step's halvings after the contraction joined and any
  200 iterations in which the step alone halved its best residual. No design
  above stalled at the budgets reported; the contraction joined on 6 of the
  240 in which each alternative dominates situations of its own, at
  evaluation 201: four that neither solved within 3,000 evaluations, and two
  that it made 1.08 times as long. On 2,995 small designs drawn to be hard
  (3 to 6 alternatives in 2 to 8 choice situations of random sizes,
  single-alternative ones included; integer covariates and constants in
  [-10, 10] and [-3, 3], the mixed logit's random-coefficient covariate in
  [-3, 3]; half with an outside option; the three models, lambda from 0.2 to
  2), at the default tolerance and 1,000 evaluations, the default failed on
  296 and the contraction on 1,742, the default on only one where the
  contraction converged (a nested logit that the contraction needed 770
  evaluations for); where both converged, the default needed about a tenth
  of the contraction's evaluations. With 200,000 evaluations the default
  failed on 73, 15 of them designs that the contraction solved, in 97,115 to
  177,800 evaluations, where T + 4E exceeds the budget. The contraction's
  turns cost evaluations where the step is slow but gets there: 1.07 to 1.20
  times the step's own from 10 far starts in a segmented market, and 1.06 to
  1.34 times on the four of the 240 that took more than 3,000 (3,607 to
  25,225 evaluations; both solved to 1e-10); more where the contraction
  keeps halving its steps (1.45 times on a nested logit with lambda 0.06
  whose contraction never improves on its start).
- **Anderson acceleration.** `blp(method = "anderson")`, and the `method`
  argument the exported functions gain, accelerate the step by Anderson
  acceleration (type II; Walker and Ni 2011). Each iteration combines the last
  10 changes in the step (`control = list(memory = )`) to cancel the current
  step in least squares, dropping the oldest while that problem's condition
  estimate exceeds 1e8 (`max_condition`), and so extrapolates along the slow
  directions the step crawls: the level shifts between groups of alternatives
  that share few choice situations. It restarts, clearing its history, when
  five iterations bring no residual below 0.9 of the best so far or the shares
  at an extrapolated point are not finite. It restarts from where it is, and
  goes back to its best point only where its residual has grown past twice the
  best, since the way to the solution can climb through larger residuals
  first; restarting from its best point every time failed on four to five
  times as many of the small designs above. It races the contraction as the
  default does, so the bound T + 4E holds (its trial points whose shares are
  not finite left out, as the step's halvings are). On the 121 simulated
  designs above it converged on every one, in medians of 13, 12 and 16.5
  evaluations of the shares (multinomial, nested and mixed logits; the default
  32, 18.5 and 174), never more than 1.25 times the default's; on the 80
  nested logits with lambda above 1 in a median 0.67 of the default's
  evaluations, never more. From 10 far starts in a market of four regions it
  needed 27 to 36 evaluations, the default 3,143 to 3,655. On the 5,985 small
  hard designs (the 2,995 above and 2,990 more drawn alike) it failed on 163,
  85 and 62 at 1,000, 20,000 and 200,000 evaluations, the default on 590, 236
  and 122, and where both converged it needed a median 0.45 of the default's
  evaluations (90th percentile 1.14). Where the targets lie within working
  precision of what the data can attain, the constants are not identified to
  that precision (conditioning above 10^7), and the extrapolation can stop
  farther from the constants the targets came from than the default does, its
  shares still within `tol`. `"simple"` remains the default.
- **SQUAREM.** `method = "squarem"` extrapolates the step by SQUAREM (Varadhan
  and Roland 2008), the default of PyBLP (Conlon and Gortmaker 2020), with its
  settings in `control` (`scheme` 3; `step_min` and `step_max` 1 at first,
  `step_max` multiplied by `step_factor`, 4, each time an accepted step reaches
  it), and races the contraction too. Each cycle takes two steps and
  extrapolates along them, three evaluations of the shares in all. It is there
  for comparison with PyBLP and the literature rather than as a default: on the
  121 simulated designs it converged on every one, in medians of 17, 14.5 and
  29 evaluations (multinomial, nested and mixed logits), never more than the
  default's; from the 10 far starts it needed 1,311 to 1,466 evaluations; and
  on the 5,985 small hard designs it failed on 630, 295 and 199 at 1,000,
  20,000 and 200,000 evaluations, more often than the default (590, 236 and
  122), above all on nested logits. With these bounds it never takes a step
  shorter than two plain ones, so it cannot damp a step that oscillates within
  a nest, as Anderson's least squares can.
- **Stopping rule.** The iteration stops at the first evaluation whose
  largest absolute log-share residual is at most `tol` (relative to the
  reference of each group without an outside option; with one, each group's
  outside share included), and returns that point. The contraction stopped
  when the largest change in the constants fell below `tol`, which left the
  outside share unchecked, and returned one more step, never evaluated.
  `max_iter` counts evaluations of the shares (0 returns the start).
- **Not converging warns.** A call that stops without meeting `tol` returns
  the evaluated point of smallest residual and warns with class
  `choicer_blp_unconverged`; the condition's `convergence` element holds the
  status (`"max_iter"`, `"non_finite"` or `"non_finite_start"`), the
  iterations, the evaluations and the residual. The contraction printed
  "Warning: Maximum iterations reached without convergence." to the console
  and returned its last step. A mixed logit start at which a predicted share
  is zero, an error before, is now such a stop.
- **Groups of alternatives.** Without an outside option, alternatives that
  share no choice situation with another group's (regional markets, say)
  have constants identified only up to the group's level, and the shares of
  each group add up to the weight of its situations. Each group's reference,
  its smallest code, is held at its start (after the start's first constant
  is subtracted, as before), where the contraction moved every level with
  the iteration.
- **Targets that cannot be attained are errors, before any evaluation.** An
  alternative offered in no choice situation of positive weight; targets
  that do not sum to 1, the outside share included, or, without an outside
  option, a group's targets that do not add up to its situations' weight,
  within a relative 1.5e-8 (within it they are rescaled); targets that leave
  no share to the outside option in a group; and a target not below the
  weight of the situations that offer its alternative, or, without an
  outside option, not above the weight of those that offer it alone. These
  conditions are necessary, not sufficient. Also errors now: a start or
  weights that are not finite, negative weights, a nested logit's lambda that
  is not positive and finite, a negative or missing `tol`, a `max_iter` that
  is not a whole number, and a mixed logit's target share that is not
  positive (the outside option's included), which stopped it only after an
  evaluation of the shares. Such inputs gave NaN, constants that drifted
  without converging, or, with lambda of 0 or below, shares without
  meaning.
- **`blp()`** checks its targets in R before the model's solver runs (and so
  before a mixed logit forms its draws): one share per alternative, the
  outside option's first, finite, positive and summing to 1. It checks
  `delta_init` the same way, and gains `method` and `control`.
  `method = "simple"`, the default, is the step above; `"anderson"` and
  `"squarem"` its accelerations, whose settings go in `control`; `"original"`
  is the contraction of Berry, Levinsohn and Pakes (1995, eq. 6.8),
  delta + log S - log s(delta), normalized as the step is and, for the
  nested logit, multiplied by `damping`. All share the checks, the stopping
  rule, the evaluation budget, the warning, the point returned and the
  halving of a step whose shares are not finite (the accelerations first
  fall back where an extrapolated point's are not); the step's stalls belong
  to the default, its cap to the step that the accelerations extrapolate
  too, and the race with the contraction to all three. Its result keeps its
  n x 1 shape and gains the inside alternatives' labels as row names and a
  `"convergence"` attribute: the
  warning's record (converged, status, iterations, evaluations, residual,
  `tol` and method), there whether or not it converged.
- **Cost per evaluation.** Each evaluation of the shares also forms the
  step's derivatives. On simulated claims- and census-style designs (8.7 and
  9.8 million rows; the mixed logit's 1.4 and 1.5 million, with 20 draws), an
  evaluation took about 4% longer than an iteration of the contraction at one
  thread and 20-28% less at eleven for the multinomial logit, about 30% and
  20% longer for the nested logit, and 3-7% longer for the mixed logit
  (minima of three runs alternated with the previous code). Setting a call
  up (the groups of alternatives, the checks) takes about half an evaluation
  at one thread and two to three at eleven (0.3 s at 10^8 rows). Fewer
  evaluations more than make up for it.
- **The exported functions** keep their arguments, defaults, layouts of
  `delta` and n x 1 results; they are now R functions over the engine, and
  gain `method` and `control` as their last arguments.
  `mxl_blp_contraction()` converts an integer or logical cube of draws to
  double, where it read such a cube from a temporary copy that R could free
  while it was read.

## Population-scale mixed logit

This release reworks the mixed logit for the data it is most needed on:
health-claims and school-census records with 10^8 stacked rows, thousands of
alternative-specific constants and many threads. The estimator is unchanged.
Two defaults change the optimizer's path (the coordinates it works in, and its
start): on the designs we measured, a fit that uses them moved within the
optimizer's tolerance or reached a maximum it used to stop short of, with the
exceptions stated below; and every fit's information is now equilibrated
before it is inverted, which moves standard errors by rounding, more where the
information is ill-conditioned ("Estimation and inference in every model",
below). The rest leaves results as they were, up to the rounding stated, or
adds an option: `polish`, a Newton phase that finishes and certifies the
optimizer's maximum; `concentrate_out_asc`, the profile likelihood of the
constants; `blp(keep_draws_bytes = )`; and `options(choicer.max_memory = )`.
`vignette("mxl")` ("Many alternatives: constants, convergence and standard
errors") puts them to work, and `?set_num_threads` covers threads and memory.
Unless stated otherwise, the figures below come from simulated designs:
claims-style panels of patients choosing among 20 to 200 of 500 hospitals, up
to 500 visits each (518 parameters), and census-style data of students
choosing among 10 to 100 of 2,500 to 2,600 schools (2,530 or 2,620 parameters,
depending on the study), with three correlated random coefficients, timed at
eleven threads, mostly as medians of runs alternated with the previous code,
some as single runs; the multinomial and nested logits' figures come from
designs of their own, described where they appear.

### Defaults that move estimates

- **`scale_vars = "bhhh"` is the default** of `run_mxlogit()`, as of
  `run_mnlogit()` (below). The optimizer works on rescaled parameters, theta =
  m * theta_t + c, each scaled by the inverse square root of the diagonal of
  the BHHH (outer product of gradients) matrix at the start values, rounded to
  a power of two, from one extra gradient pass. The likelihood is evaluated at
  the parameters, on the design as supplied; start values, bounds, estimates
  and standard errors are in natural units. L-BFGS starts every iteration from
  a multiple of the identity as its curvature, and with thousands of constants
  the curvatures of large and small alternatives' constants and of the
  Cholesky entries span orders of magnitude that no column scale reaches. With
  the default optimizer, 50 draws, a budget of 2,000 likelihood evaluations
  and the 0.2.x start values, `"bhhh"` reached the maximum (log-likelihood
  within 1e-6 of a Newton-polished reference) in 42 and 49 evaluations on
  census-style designs of 20,000 and 200,000 students, where `"none"` and
  `"sd"` had not arrived after 2,000 (up to 90 log-likelihood units short; one
  run reported convergence), and in 198 on the claims-style panel, against 988
  for `"none"`; from the new start (below), in 56, 44 and 157. A fit that uses
  the default reaches the same maximum through other coordinates: where the
  likelihood is well identified its estimates move within the optimizer's
  tolerance, and on the designs above fits that stopped short of the maximum,
  or reported a convergence they had not reached, reached it. Standard errors
  move with the estimates, and where the information is near singular (a
  variance at its boundary) they can differ widely; where the data barely
  identify a direction, estimates can move along it; and on small, degenerate
  likelihoods the two paths can stop at different points (among the test
  suite's fits `"bhhh"` stopped short in two, by 0.09 and 0.02 log-likelihood
  units: an alternative never chosen, a random coefficient's variance going to
  zero). `scale_vars = "none"` with the 0.2.x start (below) restores the 0.2.x
  path. With any other choice, settings in parameter units in `control`
  (nloptr's `xtol_abs`; optim's `parscale` and `pgtol`) and custom optimizers
  act on the scaled coordinates. `"sd"`, `"mad"` and `"iqr"` keep their
  meaning and take the same route, so the data are never divided and a scaled
  fit no longer holds a scaled copy of `X` and `W` (0.9 to 1.0 GB at 10^7
  rows). Where the diagonal is unusable at the start values, the parameters
  are left unscaled, with a message. Fits record the map in `param_scale` and
  `param_shift`; a scaled fit's `theta_init` must have one finite value per
  parameter (a vector of another length used to be recycled).
- **The default start is in the columns' units.** With `theta_init = NULL`,
  each normal random coefficient starts at L_pp = 1 / h_p on the Cholesky
  diagonal, h_p the typical step of its column: the distance between its two
  values for a column with two (1 for a dummy), the mean distance of the other
  rows from a value that fills more than half of the column, and its standard
  deviation otherwise. It started at 0.5 whatever the column's units. A step
  of h_p now moves the utility by a random amount with standard deviation 1 at
  the start, whether a price is in cents or in thousands, and rescaling a
  column moves the start with it. The coordinates above cannot repair a start
  in the wrong units: with a census-style design's columns rescaled by 10^3
  and 10^-3 the 0.2.x start stopped 233 log-likelihood units short of the
  maximum, where the new one reaches it in 54 evaluations (56 in the original
  units). Over 16 designs and rescalings the new start took at most 1.35 times
  the evaluations of the best of five starts we compared, while the 0.2.x
  start missed four of them. A log-normal coefficient keeps L_pp = 0.5, the
  standard deviation of the logarithm of its random part, which has no units;
  such a covariate still belongs in units where its coefficient is of order
  one. Over the test suite's default-start fits, those that reached the same
  maximum took 6.2% fewer iterations; 14 reached a higher log-likelihood (by
  up to 0.81), and five, all with `scale_vars = "none"` on deliberately badly
  scaled designs, a lower one (by 0.18 to 0.27), where the unscaled optimizer
  can drive a variance toward zero. On a weakly identified design whose
  supremum lies on the boundary (a degenerate covariance), `"none"` and `"sd"`
  now stop at the interior maximum that `"bhhh"` reaches, 0.021 below the
  supremum, where the second random coefficient's variance given the first
  (L_22^2) is about 1.4 instead of nearly zero. A default start outside
  `lower`/`upper` is moved to the nearest bound, with a message (nloptr used
  to refuse it), and bounds that hold no start (`lower` above `upper`,
  `lower = Inf`, `upper = -Inf`) stop the call before the fit. A `theta_init`
  of zeros with log(0.5) on the Cholesky diagonal restores the 0.2.x start.

### Convergence with thousands of constants

- **`polish`: a Newton phase after the optimizer.** `run_mxlogit()` gains
  `polish` and `polish_control`. With `polish = "newton"` the optimizer is
  followed by Newton steps on the analytic Hessian of the simulated
  log-likelihood, and with `"bhhh"` on the BHHH matrix, with the fit's own
  draws and weights. The estimator does not change: the polish finishes the
  climb that L-BFGS makes slowly where curvatures differ by orders of
  magnitude, and certifies the maximum. It converges when the Newton decrement
  nu = g'C^-1 g / 2 (g the gradient, C the curvature), over the mean weight,
  falls below `polish_control$tol_decrement` (1e-8); on a positive definite
  Hessian, every linear combination of the estimates then lies within sqrt(2
  nu) Hessian standard errors (1.4e-4) of the maximum of the log-likelihood's
  quadratic model, whatever the sample size. BHHH is positive semi-definite
  everywhere and cheaper per step, but converges linearly, and its decrement
  on its own shows only a stationary point; `"bhhh"` needs at least as many
  likelihood units as the parameters it moves. Where the curvature is not
  positive definite the steps are damped (Levenberg-Marquardt, on the
  equilibrated information); their length is capped in the optimizer's
  coordinates, which under `scale_vars = "none"` (the parameters' own units,
  in which a badly scaled column would stall them) come from one BHHH pass at
  the optimizer's solution; they respect `lower` and `upper` (projected
  Newton, Bertsekas 1982), and near the maximum of a likelihood summed over
  millions of units a line search that reads the gradient (Hager and Zhang
  2005) accepts steps whose gain is below the objective's rounding. A
  converged polish sets `convergence` to the optimizer's success code and
  leads `message` with its verdict; `fit$optimizer$polish` records the phase
  and the optimizer's own code, message and iterations. A polish that stops
  short keeps the optimizer's code and says why; one that fails (the Hessian's
  memory check, say) keeps its last accepted point and warns (class
  `choicer_polish_failed`). Its last curvature serves as the variance's own
  matrix where `se_method` needs it (after `"newton"`, the Hessian of
  `"hessian"` and the bread of `"sandwich"` and `"cluster"`; after `"bhhh"`,
  the BHHH matrix of `"bhhh"`), so where it serves the variance, a polish that
  certifies the optimizer's solution at once adds one evaluation of the
  log-likelihood and its gradient and one factorization. The constants of
  inside alternatives that have no finite estimate (see the warning under
  "Estimation and inference in every model") are held at the optimizer's
  values; where the outside option or the base alternative is the one without,
  the polish does not report convergence. A variance whose maximum is at zero
  drives its log-Cholesky diagonal toward -Inf: Newton may certify a point
  that `tol_decrement` sets, and BHHH holds the diagonal and stops unconverged
  (`?run_mxlogit`; §5.7 of the mixed logit's math article). In a validation
  study (census-style designs of 20,000 and 200,000 students, claims-style
  panels of 1,500 and 15,000 patients, alternatives never chosen dropped and
  the situations that offer one chosen wherever offered, 50 stored draws,
  eleven threads), the default fit stopped within 1.1e-4 standard errors of
  the maximum on all four designs. A Newton polish after it certified that
  maximum at once on three of them (0 to 3% of the fit, its Hessian then
  serving the variance) and in one step on the largest claims-style design
  (+20%). Handing over earlier rarely paid: of ten Newton hand-overs after 25
  to 200 evaluations one was faster (by a fifth) and one as fast; after 25
  evaluations Newton took 1.5 to 10 times as long as the default fit and, on
  the smaller claims-style design, certified another local maximum, 13.4
  log-likelihood units lower. (Times are modeled from the kernels' calls and
  per-pass costs, the machine being loaded by other work.) So
  `polish = "newton"` after the optimizer's own stop is the recommended check
  for a reportable fit with many constants; the default stays `"none"`, which
  leaves every fit as it was, bit for bit.
- **`concentrate_out_asc`: the constants concentrated out.** With
  `concentrate_out_asc = "profile_likelihood"` the optimizer moves only the
  other parameters, on the profile log-likelihood l*(theta_c) = max_delta
  l(theta_c, delta): at every point it evaluates, the constants solve their
  first-order condition, which sets each alternative's weighted chosen count
  equal to its predicted count, each decision maker's probabilities averaged
  over the draws with the weight of the draw given their choices. That is not
  BLP's share matching, which equates unconditional shares and defines another
  estimator (`"blp"` is reserved for it, and stops). The profile's maximizer
  is the joint one, and the standard errors come from the full information at
  the estimates, so the estimation of the constants enters the other
  parameters' standard errors; treating the constants as known would
  understate them (at the true parameters of a census-style design, to 45% of
  the correct standard error for a covariate that varies mostly between
  alternatives, and to 91% for an alternative x decision-maker interaction).
  Each solve takes damped Newton steps on the constants' block of the
  information, which an internal kernel forms with the objective and gradient
  in about three gradient passes (the full Hessian costs 30 to 40), from the
  constants of the best point the optimizer has reached, until the decrement
  over the mean weight falls below 1e-12 or every constant's score is below
  1e-10 of its alternative's weighted chosen count plus the mean weight
  (`concentrate_control`; at that decrement the other parameters' solution
  moves by at most 1.4e-5 standard errors, even where estimating the other
  parameters inflates the constants' variance a hundredfold, §5.8 of the math
  article). Where the constants cannot be solved at the start values the fit
  stops (class `choicer_profile_unconverged`); a solve that fails later
  reaches the optimizer as an overflow, so that it backtracks, is counted in
  `fit$optimizer$concentration`, and warns (same class) when it comes after
  the optimizer's best point or at its solution. Every constant needs a finite
  estimate: an alternative (or the outside option, or the base alternative)
  never chosen, or chosen wherever it is offered (counted by choice
  situation), is an error before the draws (class
  `choicer_unbounded_constants`), and so are bounds on the constants; up to
  `concentrate_control$max_dense = 10000` constants are factored densely (more
  is an error, before the draws; a solve holds two J x J matrices, 0.8 GB each
  at J = 10,000). It is an expert option: in the validation study above the
  profile reached the joint maximum, ending closer to it, in 21 to 49
  evaluations of the optimizer, but took 2.0 to 3.8 times the joint fit's time
  (modeled; 1.8 to 4.9 times on the clock), L-BFGS in the default coordinates
  needing 46 to 224 evaluations of the joint problem. The default, `"none"`,
  leaves every fit as it was.

### Standard errors with thousands of constants

- **The Hessian and the BHHH matrix work, per decision maker, on the
  parameters that decision maker's likelihood involves**: the fixed
  coefficients, the random coefficients' means and Cholesky factor, and the
  constants of the alternatives that decision maker faced (59 of 2,620
  parameters per student on average). They used to carry every constant
  through every decision maker. The Hessian behind the default standard errors
  took 0.35 s instead of 24.8 s for 2,000 students (1.9 s instead of 101 s at
  one thread) and 3.2 s instead of 228 s for 20,000; the BHHH matrix 0.15 s
  instead of 25.1 s for 20,000 and 1.4 s instead of 246 s for 200,000. On the
  claims-style panel, where a patient's visits involve about 60 of the
  hospitals, the Hessian took about half the time (26 s instead of 55 s for
  15,000 patients).
- **Robust and clustered variances no longer form the matrix of every decision
  maker's score** (42 GB for 2 million students with 2,620 parameters, held
  twice while it was weighted). The robust meat of `vcov(type = "robust")` and
  `wesml_vcov()` is now the BHHH matrix with squared weights, as that of
  `se_method = "sandwich"` already was, and the clustered meat of
  `vcov(type = "cluster")` and `se_method = "cluster"` sums the weighted
  scores within clusters in C++: for 200,000 students, the robust meat took
  1.4 s and 0.4 GB of working memory instead of 12.2 s and 8.4 GB. With many
  thousands of constants the Hessian, the BHHH matrix and the clustered meat
  keep one copy of the result per thread only while the copies fit in 2 GiB
  (about 7,000 parameters at eleven threads); above that the threads add into
  one result: at 17,156 parameters, the Hessian for 20,000 students took 2.9 s
  and 2.4 GB at eleven threads, where per-thread copies would need 13 GB. The
  internal score kernel refuses to form a matrix above 2 GiB.
- The compact blocks and the meats computed in C++ alter standard errors by
  floating-point rounding only: on our reference battery the Hessian by at
  most 3.5e-15 relative, the meats by 7e-16, and the BHHH matrix not at all.
- A mixed logit's variance at fit time is computed by the route
  `vcov(fit, type = )` takes post hoc (`"robust"` for
  `se_method = "sandwich"`), with the fit's own draws, so the two agree, bit
  for bit on one thread (draws supplied through `eta_draws` are not
  regenerated post hoc). A scaled fit used to invert the information of its
  scaled design and transform the result back, while `vcov(fit, type = )`
  inverted the natural information as it was; both lost digits where the
  columns span orders of magnitude, which the equilibrated inversion
  ("Estimation and inference in every model") keeps (on one test's robust
  variance, whose information had condition number 1.2e12, standard errors
  computed raw were off by 2.1e-7 relative against an exact computation,
  equilibrated by 9e-16).

### Memory and speed

- **The estimation kernels** keep each thread's working arrays across decision
  makers, read the stacked design's integer indices in place, and compute each
  decision maker's base utilities inside the parallel loop, so a likelihood
  evaluation allocates nothing the length of the design (what remains is 8
  bytes of row offsets per choice situation, and 8 per decision maker in a
  panel). Each decision maker's draws are processed in batches within a 4 MiB
  budget per thread, the score accumulated across batches by a streaming
  log-sum-exp, so a thread's memory no longer grows with a long panel times
  `S`. On a claims panel with up to 500 visits per patient, peak working
  memory at `S = 1000` fell from 6.4 GB to 0.12 GB and an evaluation ran 40%
  faster; on 10^8-row panels peak working memory fell from about 3 GB to under
  0.2 GB; the analytic Hessian's heap traffic on a 10^6-row panel fell from
  106 million allocations per call to about three thousand. These kernel
  changes by themselves alter results only by floating-point rounding (at most
  3.4e-15 relative on our reference battery, and within 1e-10 for decision
  makers split into draw batches).
- **Generate-mode draws** (`draws = "generate"`) are formed in two passes over
  each block, all the Halton uniforms first, then the inverse normal, with the
  first dimension stepped by an exact integer odometer and the digits of the
  next ten taken by division by a compile-time constant. The draws are
  bit-identical on 64-bit platforms; gradient evaluations ran 32-36% faster on
  a four-alternative cross-section, where regenerating the draws took about
  half the time, and 10-14% faster on census-style data. On the
  four-alternative data a generate-mode evaluation now takes about 1.35 times
  as long as a store-mode one, down from 2.
- **Before allocating, the estimation routines check their memory.** The
  routines behind the log-likelihood and gradient, the Hessian (and the
  profile's block of the constants), the BHHH matrix, the score matrix, the
  meats and the conditional tastes compare an upper estimate of their working
  memory, at the threads set, with the machine's physical memory, and stop
  with an error that says how much they need and how many threads fit, where
  their threads' memory used to grow until an allocation failed and terminated
  R. They never lower the thread count on their own. Memory already in use is
  not counted, nor is the inversion in R of the resulting
  parameters-by-parameters matrices. `options(choicer.max_memory = <bytes>)`
  sets the limit instead (a container's or a batch job's), and `Inf` lifts the
  check (`?set_num_threads`).
- **Prediction and post-estimation at scale.** The kernels behind `predict()`,
  `logsum()`, `consumer_surplus()`, `elasticities()`, `diversion_ratios()`,
  `blp()`, `gof()` and the fit statistics of `summary()` read the alternative
  codes in place, simulate each choice situation in buffers allocated once per
  thread (where they allocated about ten objects per situation), and allocate
  their working arrays before the parallel loop, so that running out of memory
  is an R error. On census-style data of 9.8 million rows, choice
  probabilities took 16-20% less time and moved 0.16 GB through the heap
  instead of 11.8 GB; log-sums 26-27% less, 0.08 GB instead of 18.8 GB. With
  stored draws (`draws = "store"`, the default) these functions now form each
  choice situation's block of Halton draws on the fly, in the thread that
  simulates it, instead of building a `K_w x S x N` cube on every call:
  `predict()` on census-style data of 98 million rows (a 5.9 GB cube) took
  9.0 s instead of 76 s and its working memory peaked at 1.6 GB, the results
  it returns, instead of 14.4 GB; and these functions no longer stop past
  store mode's 2^31 - 1 points. The draws, and so the results, are the cube's
  bit for bit where randtoolbox and choicer are compiled with the same
  floating-point contraction, as on CRAN's platforms; otherwise some of them
  move by a few parts in 10^9, in the far tails. With more than 128 random
  coefficients they build the cube, as before. `blp()` gains
  `keep_draws_bytes` (default 2^32 bytes, 4 GiB): when the draws of all choice
  situations fit, its contraction forms them once and reuses them at every
  iteration, leaving the deltas as they were (on the claims-style panel at 85
  million rows, whose 3.6 GB of draws fit the default, each evaluation of the
  shares took 5.2 s instead of 6.4 s with stored draws and 5.6 s instead of
  7.0 s with generated ones, in single runs, two alternated ones for the
  stored draws); the memory is taken per call, and a budget above the memory
  actually free can end the R session. The exported `mxl_blp_contraction()`
  accepts `gen_scramble = 2` for these draws, and now rejects codes other than
  0, 1 and 2, which it treated as 0. Given the same draws, the prediction
  kernels' results are unchanged by this rework, bit for bit, with OpenBLAS or
  the reference BLAS: each choice situation's probabilities, mean utilities
  and log-sums at any number of threads, and sums over situations at one
  thread; with other BLAS libraries they may differ in their last bits.
  `elasticities()` and `diversion_ratios()` keep a J x J matrix per thread
  only while the threads' copies fit in 2 GiB, and check that memory first
  ("Estimation and inference in every model", below).
- **Stored draws for estimation.** `get_halton_normals()`, which builds the
  estimation cube of `draws = "store"` (and again for `vcov(type = )`,
  `wesml_vcov()`, `conditional_tastes()` and standard errors computed later),
  fills it a block of units at a time: 21 s instead of 27 s, and 1.4 instead
  of 3.0 times the cube's size at peak, for a 4.8 GB cube of 2 million units.
  It checks its inputs and stops before allocating more than 2^31 - 1 points
  (`S * N`), the largest start `randtoolbox::halton()` accepts; generate mode
  has no such limit, and choicer now requires randtoolbox 1.31.0 or later.
  `run_mxlogit(draws = "store")` (the convenience workflow) warns before it
  builds draws that would take more than 1 GiB (`8 * K_w * S * U` bytes): the
  cube is held while the model is fitted and rebuilt by the functions above.
  Under `options(warn = 2)` such fits stop before estimating.

## Estimation and inference in every model

- **`scale_vars` in the multinomial and nested logits.** `run_mnlogit()` now
  defaults to `scale_vars = "bhhh"`, as `run_mxlogit()` does (above), and
  `run_nestlogit()` gains `scale_vars` (its last argument) and keeps `"none"`.
  With the default optimizer, the multinomial logit reached the maximum in 26
  and 39 evaluations on the claims- and census-style designs (509 and 2,521
  parameters), against 787 and 1,358 with `"none"`; on that census-style
  design with its columns rescaled by 10^3 and 10^-3, `"none"` reported
  convergence 1,920 log-likelihood units short, where `"bhhh"` took 38
  evaluations. The nested logit keeps `"none"`: its `"bhhh"` scales are fixed
  at the start values, and nest parameters that end far from their start of
  0.5 can mislead the optimizer (14% more iterations over the test suite's
  nested logit fits, one of them 540 against 29). As for the mixed logit,
  settings in parameter units in `control` and custom optimizers act on the
  scaled coordinates, and `scale_vars = "none"` restores the 0.2.x path.
  `run_mnlogit(scale_vars = )` no longer holds a scaled copy of `X` (on 8.7
  million rows and 10 columns, a scaled fit's peak heap fell from 1.39 to
  0.10 GB, an unscaled fit's being 0.08 GB). With `optimizer = "optim"`,
  BFGS's `reltol` (multinomial logit, about 1.5e-8 of the log-likelihood) and
  L-BFGS-B's `factr` (mixed and nested logits, about 2.2e-9) can stop before
  full accuracy, scaled or not: on the claims-style multinomial logit optim
  reported convergence 0.44 log-likelihood units short with `"bhhh"` and 0.52
  short with `"none"`; tighten them through `control`, where `maxit` (100 by
  default) caps the iterations.
- **`polish` in the multinomial and nested logits.** `run_mnlogit()` and
  `run_nestlogit()` gain `polish` and `polish_control` (their last arguments),
  the Newton phase of `run_mxlogit()` (above): from wherever the optimizer
  stopped, `polish = "newton"` takes Newton steps on the analytic Hessian and
  `"bhhh"` on the BHHH matrix, with the fit's weights, until the Newton
  decrement over the mean weight falls below `tol_decrement` (1e-8), which
  certifies a maximum (a local one in the nested logit, whose log-likelihood
  is not concave; with `"bhhh"`, on its own, a stationary point). The
  estimator does not change, and the default `"none"` leaves every fit as it
  was, bit for bit. The record (`fit$optimizer$polish`), the convergence code
  and message, the variance's reuse of the last curvature (never for the
  nested logit's `se_method = "numeric"`, whose Hessian is numerical) and the
  constants held at the optimizer's values are the mixed logit's. In the
  nested logit a never-chosen alternative's constant is held only where its
  nest is a singleton or its dissimilarity is at most 1 at the optimizer's
  solution: where a nest-mate is chosen its score is P(j|k)(1 - 1/lambda_k -
  P_k), which can be positive above 1, so the constant can then have a finite
  estimate, and a polish whose held constant's nest ends above 1, or whose
  moved one's ends at or below 1, does not report convergence. Under
  `scale_vars = "none"`, the nested logit's default, the steps' scales come
  from one BHHH pass at the optimizer's solution, as for the mixed logit
  (with half of the covariates rescaled by 10^3 and the other half by 10^-3,
  20 steps at unit scales left a decrement of 1,054 to 2,060, where 8 or 9 at
  these certified the maximum). A trial point fails where the likelihood
  kernel reports that its objective or gradient is not the log-likelihood's
  (a log probability clamped to -1e10, a total that is not finite, a gradient
  entry set to zero), which the optimizers do not see. In a study on census-
  and claims-style designs (2,504 and 2,602 alternatives, in 50 nests for the
  nested logit, and 24,560 and 245,031 choice situations; 499 and 500
  alternatives, and 17,339 and 153,721 situations; constants without a finite
  estimate removed; single runs at eleven threads), the multinomial logit's
  default fit stopped at the maximum and a Newton polish certified it at
  once: with its Hessian then serving the variance, it added one gradient
  pass and one Cholesky factorization, 3 to 13% of the fit and its variance;
  with `scale_vars = "none"` the optimizer stopped up to 9.2 log-likelihood
  units short, a gap the polish closed in 1 to 6 steps (up to a tenth of the
  fit). The nested logit's default fit reached its 1,000-evaluation limit on
  three of the four designs, 0.35 and 15.8 units short on the census-style
  ones (up to 6.7 standard errors); a Newton polish certified the maximum
  there in 3 and 5 steps, for about a tenth of the fit's time, and in one
  step on the claims-style ones. Stopping the optimizer sooner
  (`control = list(maxeval = 100)`) and polishing took about a fifth of the
  time of the default fit and its polish, and with `scale_vars = "bhhh"` the
  optimizer stopped at the maximum in 50 to 66 evaluations on three of the
  four designs. (Times are modeled from the kernels' calls and per-pass
  costs, the machine being loaded by other work.) So `polish = "newton"`
  after the optimizer's own stop is the recommended check for a reportable
  fit with many constants (`?run_mnlogit`, `?run_nestlogit`; §12.7 of the
  multinomial logit's math article, §5.9 of the nested logit's).
- **The multinomial and nested logit kernels** read the alternative codes, the
  choices and the nests in place, where they copied them into index vectors on
  every call, and the nested logit computes each situation's probabilities in
  buffers allocated once per thread. On 17.5 million rows a multinomial logit
  gradient evaluation took 23% less time at one thread and 39% less at eleven,
  and peaked at 0.14 instead of 0.57 GB of working memory. On 7 million rows a
  nested logit gradient evaluation took 37% and 46% less. This change by
  itself leaves results unchanged, bit for bit. Malformed index vectors, which only a hand-built `input_data` or a
  direct call can supply, are now reported in a defined way (a missing or
  negative code as such, an index beyond the integer range as `NA`).
- **The multinomial and nested logits' variances at population scale.** The
  Hessian and BHHH kernels of `run_mnlogit()` and `run_nestlogit()`, and their
  new clustered-meat kernels, form each choice situation's terms only over the
  parameters its likelihood involves (the coefficients, the constants of the
  alternatives it offers and, in the nested logit, the dissimilarities of the
  nests it offers) and add them into upper-triangular accumulators, as the
  mixed logit's kernels do (above). The work per situation no longer grows
  with the number of alternatives: on a census-style design of 24,560 choice
  situations and 2,611 parameters (2,661 in the nested logit, with its 50
  nests), the multinomial logit's Hessian took 0.13 s instead of 281 s at one
  thread and its BHHH matrix 0.07 s instead of 122 s, and the nested logit's
  BHHH matrix 0.14 s instead of 90 s; at five threads the multinomial logit's
  Hessian moved 0.22 GB through the heap instead of 1,336 GB. With 17,147
  parameters, the multinomial logit's Hessian, which the dense kernel took
  358 s to form on 500 of the situations at five threads, takes 0.30 s on all
  24,560. On designs with a handful of alternatives, where a situation touches
  every parameter, the BHHH matrix takes less time than with the dense kernel:
  0.35 instead of 0.65 s in the multinomial logit and 0.67 instead of 0.92 s
  in the nested logit on 2.5 million situations at one thread (0.68 and
  1.10 s in 0.2.1; see the next item). The robust
  meat of `vcov(type = "robust")` and `wesml_vcov()` is now the BHHH kernel
  with squared weights, as that of `se_method = "sandwich"` already was, and
  the clustered meat comes from a kernel, so no variance route forms the
  matrix of every situation's scores (on the census-style design the clustered
  meat took 0.08 s instead of 2.8 s at one thread, and the post-hoc robust
  meat peaked at 0.20 instead of 1.08 GB of heap at five threads); the
  internal score kernels, which no route calls, stop above 2 GiB. A fit's
  variance is the one `vcov(fit, type = )` returns post hoc (`"robust"` for
  `se_method = "sandwich"`, `"hessian"` for the nested logit's `"numeric"`),
  bit for bit at one thread, for every `se_method`: before, the fit's BHHH and
  sandwich variances differed from the post-hoc ones in their last bits. Every
  multinomial and nested logit kernel behind the log-likelihood, its gradient,
  Hessian, BHHH matrix, scores and clustered meat now compares an estimate of
  its memory with `options(choicer.max_memory)` or the machine's physical
  memory before allocating, and stops with an error that gives the number of
  threads that fits, or says that even one does not (see `?set_num_threads`);
  their accumulators are allocated before the parallel loops. Estimates are
  unchanged, bit for bit at one thread. The variances move by rounding: on our
  reference battery, relative to each matrix's largest entry, the multinomial
  logit's Hessian by at most 4.3e-18 and only in designs with a single
  coefficient or a single constant, the nested logit's by 3e-17 (1.4e-16 where
  a utility is infinite, a covariate missing or a weight not finite), the BHHH
  kernels' matrices not at all, and every variance by at most 1.5e-12 of its
  entries' scale; among the test suite's fits, where an ill-conditioned
  information amplifies such moves, by up to 1.7e-5 (a nested logit whose
  information has condition number 2e11). Where a variance is zero in exact
  arithmetic (a clustered variance with fewer clusters than parameters, say),
  rounding decides its sign, so its standard error can turn from `NA` to a
  tiny number, or back.
- **The multinomial and nested logits' work per choice situation.** Their
  compiled routines (the log-likelihood and gradient, the Hessian, BHHH
  matrix, scores and clustered meat, predictions and shares, elasticities and
  diversion ratios) no longer allocate memory
  for each choice situation: every thread's working arrays are allocated once
  per call, before the parallel loop (on 250,000 four-alternative situations
  a multinomial logit BHHH call made 250,023 heap allocations and moved 90 MB
  through the heap, against 26 allocations and 10 MB now; the nested logit's
  predictions, shares and diversion ratios already allocated nothing per
  situation). A thread keeps a situation's bookkeeping while the
  next situations offer the same alternatives, and the loops hand a thread
  about 256 stacked rows at a time instead of a single situation. The nested
  logit's loops over nests cover only the nests a situation offers; so do its
  probabilities, where every dissimilarity is positive and the nests
  outnumber a situation's alternatives on average. At one thread, on 2.5
  million four-alternative situations and on census- and claims-style designs
  (24,560 situations of 40 alternatives on average, and 153,721 of 56.5), the
  multinomial logit's gradient took 44%, 46% and 34% less time than before
  this change, its BHHH matrix 53%, 23% and 20% less, and predictions 34%,
  26% and 26% less; the nested
  logit's gradient 22%, 18% and 3% less and its BHHH matrix 39%, 25% and 23%
  less (two nests on the four-alternative design, 50 on the others). At five
  threads the gains are smaller: on the census- and claims-style designs the
  multinomial logit's gradient took 27% and 20% less time and its BHHH matrix
  12% and 16% less. On the census-style design with 500 nests, the nested
  logit's gradient, predictions, shares and contraction iterations (the
  contraction since replaced: "Share inversion", above) took 41-48% of the
  time they took with the probabilities computed over every nest.
  Results are unchanged, bit for bit at one thread.
- **Elasticities and diversion ratios with many alternatives.**
  `elasticities()` and `diversion_ratios()` of multinomial, nested and mixed
  logit fits add each choice situation's terms into one J x J matrix per
  thread only while those copies fit in 2 GiB, and beyond that into one
  shared matrix; with one thread they add into the result itself, which
  halves their peak memory (62 instead of 117 MB with J = 2,602
  alternatives). They now check that memory before allocating, as the
  estimation routines do, and stop with the same error, which names the
  number of threads that fits (`?set_num_threads`): the multinomial and nested
  logits used to allocate every thread's matrix inside the parallel loop,
  where running out of memory ended the R session, and the mixed logit
  allocated its copies without a check. With 17,138 alternatives and the
  outside option (a 2.35 GB matrix), the multinomial logit's elasticities took
  0.35 s at five threads, where five per-thread copies and the result would
  have taken about 14 GB. With more than one thread the last bits of the
  result vary between calls, as before: with the situations each thread takes
  and, beyond the 2 GiB, with the order of the shared additions; at one
  thread results are unchanged. For a hand-built design whose largest
  alternative code is 2^31 - 1, with the outside option, the multinomial and
  nested logits' elasticity and diversion-ratio routines now stop with the
  mixed logit's message.
- **Cluster labels.** `vcov(type = "cluster")` and `cluster_col` code plain
  factor, integer, logical and character labels without converting them to
  strings: coding 10^7 integer labels in 2 million clusters took 0.41 s
  instead of 4.9-5.6 s. Other labels (doubles, dates, integer64, factor
  subclasses, factors with repeated levels or levels marked `"bytes"`) go
  through `as.character()` as before, and the clusters are the same either
  way.
- **Data preparation at population scale.** `prepare_mnl_data()`,
  `prepare_mxl_data()` and `prepare_nl_data()` copy only the index columns
  into their working table, read the covariates in place from the caller's
  data, scan them for missing values a column at a time, and gather `X` and
  `W` in C++ once the rows are filtered and sorted; `X` and `W` are now always
  double precision (an all-integer design used to be an integer matrix), which
  by itself leaves estimates unchanged, bit for bit. The collinearity check
  factors a design of a million or more values by row chunks and applies
  `qr()`'s own rank rule to the resulting factor: it drops the same columns
  (the two can disagree only on a column whose dependence lies within rounding
  of the tolerance), makes none of `qr()`'s copies of the design, and no
  longer stops at its limit of 2^31 - 1 values ("too large a matrix for
  LINPACK"). The checks that a weight, cluster or decision-maker column is
  constant within each choice situation no longer sort every situation's
  values (87 of the 97 GB `prepare_mxl_data()` moved through the heap). On
  claims- and census-style panels of 8.7 and 9.8 million rows,
  `prepare_mxl_data()` took 2.9-3.7 s instead of 5.7-7.6 s, allocated
  7.5-8.8 GB instead of 97-149 GB, and peaked at 1.9 times its input on the
  heap instead of 4.0-4.5 times (2.0-2.1 instead of 4.4-5.7 for the
  multinomial and nested logits); per-situation counts now computed in one
  pass took another 9-15% off these three. The preparations of
  `run_mnprobit()`, `run_hmnlogit()` and `run_hmnprobit()` were reworked the
  same way (on 10^7 rows, `prepare_mnp_data()` took 5.0 s instead of 13.8 s
  and `prepare_hmnl_data()` 14.5 s instead of 31.7 s, at peak heap 3.0 to 3.1
  times the input instead of 4.9; counting per situation in one pass then
  halved the first's time and cut the second's by two fifths), with three
  changes: the choice-situation id used as a covariate is recognized as
  constant within situations (kept with a message, or dropped with a warning
  without an outside option); integer covariates whose differences overflow
  32-bit integers no longer fail; and inputs whose column names clashed with
  the preparations' own are covered under Corrections. Prepared objects are
  otherwise unchanged on every call the test suite makes and on a battery of
  edge cases.
- **`optimizer = "optim"` reuses the likelihood pass for the gradient.**
  `stats::optim()` asks for the objective and its gradient through separate
  functions, and each call ran the likelihood pass, which computes both; the
  gradient now comes from the pass that computed the objective at the same
  point. L-BFGS-B, which fits the mixed and nested logits, takes about half
  the evaluations (14,104 to 7,087 for the mixed logit over the test suite's
  fits refitted with optim); BFGS, which fits the multinomial logit, one fewer
  per gradient. This change by itself leaves estimates, standard errors,
  iteration counts and convergence codes unchanged, bit for bit on one thread.
  nloptr, the default optimizer, is unchanged.
- **Every fit's information is equilibrated before it is inverted.** The
  standard errors of every fit of `run_mnlogit()`, `run_mxlogit()` and
  `run_nestlogit()`, at fit time and post hoc (`vcov(fit, type = )`,
  `wesml_vcov()`), come from the information matrix scaled on both sides by
  the powers of two nearest the inverse square roots of its diagonal's
  magnitudes (which adds no rounding), inverted, and scaled back; the sandwich
  is solved twice against the bread's factors without forming its inverse. In
  natural units an information matrix's curvatures can span many orders of
  magnitude (columns in different units, the constants of rarely and often
  chosen alternatives); inverted as it was, such a matrix could fail as
  computationally singular and leave the standard errors `NA`. Against an
  exact computation on 963 inversions from our reference battery and the test
  suite, the equilibrated route recovered 23 of the 30 the old one refused
  (standard errors within about 1e-15 relative of exact), lost none, and where
  an inversion changed its error stayed within 2.7 times the old one. (The
  math articles count against a raw `solve()`, which refused 52 of the 963;
  the old route already equilibrated 22 of those.) Where the information is
  ill-conditioned for reasons no rescaling removes, standard errors keep only
  the digits its conditioning allows, and the last of them move (among the
  test suite's fits by up to 1e-4 relative, on nested logits); a variance that
  is zero in exact arithmetic can come out `NA` with a message, or the
  reverse. An inverse holds two parameters-by-parameters matrices besides the
  information, where `solve()` held three (four or more when the information
  was scaled), and the sandwich two besides the bread and meat, where it held
  four or five: 4.6 GB against up to 11.6 GB at 17,000 parameters.
- **A variance that fails no longer discards the estimates.** An error in the
  variance step (`se_method`), such as a kernel's memory check or an
  allocation that fails in R, used to stop the call and lose the estimates.
  The fit now returns them without a variance (`vcov` is `NULL`, the standard
  errors `NA`, and `variance_error` holds the message) and warns (class
  `choicer_vcov_failed`) with the error and how to compute the variance later,
  `vcov(fit, type = )`; nothing the size of the variance is allocated after
  the failure. `summary()` and the post-estimation methods work as on any fit,
  with `NA` wherever a standard error enters, and do not rerun the failed
  step. An allocation that fails while the information is inverted, which was
  reported as a singular information (and recomputed at every `summary()`),
  takes the same route at fit time; after the fit it is an error, from
  `vcov(fit, type = )` or from the methods that compute a missing variance.
  The variance's inputs are now checked before the likelihood is first
  evaluated: the guard on choice-based-sampling weights (whose warnings are
  now printed as they occur) and the fit-time cluster labels. Fits whose
  variance can be computed are unchanged by this, bit for bit. An allocation
  that fails inside the multinomial and nested logits' Hessian, BHHH and
  cluster-meat loops, which now allocate only small per-situation buffers
  there, still terminates R, as one in their elasticity and diversion-ratio
  kernels (a J x J matrix per thread, unchecked) can. Other inputs that can
  fail are checked before the expensive work too: `run_mxlogit()` computes the
  column scales of `"sd"`, `"mad"` and `"iqr"` (and stops on a column without
  spread) and checks a supplied start where the optimizer works in scaled
  coordinates before it builds stored draws, and the three wrappers check that
  every parameter has a name before the optimizer runs, where a hand-built
  `input_data` whose `X` lacks column names, or a `param_names` of another
  length, used to fail only when the estimates were named, after it (an
  unnamed one-column `W` with `rc_mean = TRUE` gave a mean named `Mu_`; it now
  needs a name too).
- **A warning for constants without a finite maximum-likelihood estimate.**
  `run_mnlogit()`, `run_mxlogit()` and `run_nestlogit()` warn (class
  `choicer_unbounded_constants`), before the optimizer runs, when the counts
  in `alt_mapping` show an alternative never chosen, or chosen in every choice
  situation that offers it: its constant's score has one sign at every value
  of the parameters, so the optimizer stops it wherever its tolerance does and
  its standard error means nothing, while the other parameters approach the
  maximum without that alternative (or without the situations that offer it).
  Where the outside option or the base alternative is the one, the other
  constants diverge together. The warning names the alternatives and how to
  reach finite constants: drop a never-chosen alternative or pool it with a
  chosen one; drop the situations that offer an always-chosen one or pool it;
  fit without an outside option that is never chosen, as when the data hold
  the outside option's own rows and `include_outside_option = TRUE` adds
  another. In the mixed logit a finite bound in the direction a constant
  diverges holds it, and that alternative is not named; in the nested logit
  the never-chosen case needs the nest's dissimilarity parameter at most 1. It
  is printed as it occurs, before a fit that can take hours. The counts are of
  rows, so an alternative listed twice in a situation can hide being chosen
  wherever offered, and other patterns (a group of alternatives chosen only in
  situations that offer no alternative outside it) are not detected. Estimates
  are unchanged.

## Mixed logit — panel likelihood (`person_col`) and conditional tastes

- `run_mxlogit()` / `prepare_mxl_data()` gain `person_col`: the frequentist
  panel mixed logit of Revelt and Train (1998). When supplied, all choice
  situations of a decision maker share one draw of the random coefficients,
  and the likelihood integrates the *joint* probability of that person's
  choices over the taste distribution, instead of integrating every choice
  situation separately as before. `id_col` must still identify choice
  situations uniquely across the whole data set (ids do not restart within
  person): an id under two decision makers is an error, found before
  situations with missing or infinite values are dropped (they are dropped
  by id, which would remove such an id from every decision maker), over the
  rows whose id and decision maker are known. A choice situation whose
  decision maker is missing is dropped like one with any other missing
  value, so a shared id one of whose copies has no decision maker in any
  row is dropped whole rather than reported. `NULL` (the default) keeps
  every choice situation its own decision maker, reproducing the existing
  cross-sectional likelihood.
  - `simulate_mxl_data()` gains `T` (choice situations per decision maker;
    default `1`, a cross-section) and records the realized random coefficients
    in `true_params$gamma_i` (`K_w x N`), the realized tastes that
    `conditional_tastes()`'s conditional means can be compared with (they
    track them with shrinkage). `T = 1` reproduces earlier versions' data
    exactly.
  - New exported generic `conditional_tastes()` (class `choicer_tastes`,
    method for `choicer_mxl`): the mean and SD of a decision maker's random
    coefficients conditional on the choices they made (Revelt and Train 2000;
    Train 2009, ch. 11) — the frequentist analogue of the hierarchical Bayes
    `beta_i` summaries of `run_hmnlogit()`. Its `print()` reports, per
    coefficient, the share of population taste variance revealed by choices
    and the residual share (law of total variance), which should sum to
    about one — a numerical check on convergence and the draw count `S`, not
    a specification test (both identities are first-order conditions of the
    likelihood, so they hold at a converged fit whatever the true taste
    distribution).
  - Inference: with `person_col`, scores are per decision maker, so
    `se_method = "sandwich"` / `vcov(type = "robust")` is already robust to
    within-person dependence; `cluster_col` / `cluster` labels must then nest
    decision makers (be constant within each). Weights become decision-maker
    weights, constant within person (the objective is `sum_n w_n log L_n`).
    WESML (choice-based) weighting is not supported with `person_col`: a
    choice-based weight varies with the alternative a person chose, which can
    differ across their situations.
  - Fitted-object fields: a panel fit gains `n_persons` (the number of
    decision makers; printed as "Respondents"), `data$Ti` (choice situations
    per decision maker, in prepared order), and `data$person_ids`. `nobs()`
    is unchanged — it still counts choice situations, as for HMNL.
  - New recovery script `inst/simulations/mxl_panel_simulation.R`, in the
    style of `mxl_simulation.R`: simulates a panel, fits it with
    `person_col`, reports `recovery_table()` and `conditional_tastes()`, and
    contrasts standard errors against a cluster-robust cross-sectional fit of
    the same data.
  - The panel likelihood leaves the cross-sectional likelihood
    (`person_col = NULL`) as it was, up to floating-point reassociation
    (<= 1e-10 relative in the estimation kernels on our reference battery of
    test configurations). The other numerical changes of this release (the
    new defaults, the equilibrated variance, the corrections) are listed under
    "Population-scale mixed logit", "Estimation and inference in every model"
    and "Corrections".
  - New input guards: a non-`NULL` empty `Ti` is an error, and a
    primary-thread memory check now stops when a single decision maker
    stacks so many alternative rows (tens of millions) that its design rows
    would need more than 2 GiB of scratch in every thread — a signal that
    `person_col` identifies markets or some other high-cardinality grouping
    rather than decision makers.
  - Alignment guards for panel fits, whose prepared situations are ordered by
    decision maker rather than by id: a positional `weights` vector and an
    unnamed post-hoc `cluster` vector are rejected (with a pointer to
    `weights_col` / id-named labels) whenever that order differs from id
    order, and `predict()` / `logsum()` / `consumer_surplus()` note when
    `newdata` lacks the person column (rows are then ordered by id rather than
    in the in-sample, decision-maker-first order).
  - `conditional_tastes()` on other fits errors informatively, pointing
    hierarchical Bayes users to their `beta_i` summaries.
  - Supersedes the v0.2.0 "Scope note (MXL)" below: clustering a
    cross-sectional fit still only repairs the inference, not the efficiency
    lost by ignoring within-person correlation, but the likelihood itself is
    now a choice: pass `person_col` to use the within-person variation the
    panel provides.

## Corrections

The following were found during this development cycle, while implementing
the panel likelihood and the population-scale work above or reviewing it,
and independently verified. Each affected released versions unless it says
otherwise.

- `mxl_hessian_parallel()` silently dropped any choice situation whose
  simulated choice probability, summed over draws, was `<= 1e-12`, while the
  objective and gradient kept it. A situation that was merely very poorly fit
  by the current parameter vector — not one with genuinely zero simulated
  probability — could therefore be silently excluded from the analytical
  Hessian, which could then be substantially wrong, and so could the
  inverse-Hessian and sandwich standard errors built from it. Such situations
  are now kept; only units whose log-likelihood is not finite, which now
  takes utilities that overflow (see the log-probability item below), are
  skipped.
- When a choice situation's simulated probability was subnormal (utility gaps
  of several hundred), the gradient came back identically zero alongside a
  finite objective value — a false stationary point that could stop the
  optimizer early — and the score-based variances (BHHH, robust, cluster)
  were non-finite. The posterior draw weights are now normalized in log space
  (matching the panel likelihood's `omega_ns = exp(lambda_ns - LSE_s
  lambda_n.)`), so these quantities are finite whenever the simulated
  probability is positive (and, with the next item, accurate).
- The mixed logit kernels took each choice situation's log probability as
  `log(P)`. A choice probability that underflowed to zero, at a utility gap
  beyond about 745, made its decision maker's log-likelihood `-Inf` and with
  it the whole objective the optimizer's `1e10` sentinel with a zero
  gradient; subnormal probabilities, at gaps beyond about 708, lost
  precision. The chance of hitting this grows with the number of choice
  situations: on a synthetic claims-style panel of 8.7 million stacked rows
  with a single choice made against a covariate value of 2000, the fit
  returned the sentinel at 342 of 540 evaluations and stopped where that
  choice's probability underflows (a coefficient of -0.37, truth -0.5). The
  log probability is now the chosen utility minus the log-sum-exp of the
  choice set, as in the multinomial logit kernels, so the objective,
  gradient, Hessian, scores and conditional tastes are exact whenever the
  utilities are finite, however improbable the observed choices; the same
  fit now converges (-0.48). The sentinel, the Hessian's skip of a unit and
  `NA` conditional tastes now occur only when the utilities themselves
  overflow (for instance an exploding Cholesky factor during a line search),
  and because finite objectives are no longer bounded below the sentinel,
  `run_mxlogit()` reports it as ten times the largest objective seen along
  the optimizer's path. An explicit kernel overflow flag distinguishes the
  sentinel from a valid objective equal to `1e10`. Fits away from such
  regions agree with earlier versions to optimizer tolerance rather than
  bit for bit, since the optimizer's path shifts slightly.
- `mxl_hessian_parallel()` now centers the draw scores in the Louis
  identity, `sum_s omega_s (g_s - g_bar)(g_s - g_bar)'` in place of
  `sum_s omega_s g_s g_s' - g_bar g_bar'`. The two agree in exact arithmetic,
  but the former draw weights summed to one only to the precision of their
  log-scale arguments, and the uncentered form multiplied that error by
  `|g|^2`: for a choice far below its competitors, an entry of about 1 came
  out 0.15 off.
  Posterior weights now divide max-shifted exponentials by their sum,
  including in the batched score, and Hessian centering first removes a
  highest-weight draw's score. This also prevents spurious curvature when a
  finite utility gap is so large that subtracting the log-sum-exp loses
  normalization accuracy.
  The weights described here were this development version's log-space
  weights; released versions normalized linear weights, which sum to one to
  rounding precision, so there the uncentered form lost accuracy only at far
  larger scores.
- `run_mxlogit()`'s advanced workflow (`input_data` + `eta_draws`) recorded
  the `S` argument (default 100) instead of the number of draws in
  `eta_draws`, so post-hoc recomputation (`vcov(fit, type = )`, lazily
  computed standard errors, prediction) could regenerate a different draw
  set from the one used in estimation. `draws_info$S` now records
  `dim(eta_draws)[2]`. Post-hoc methods reproduce the estimation draws only
  when `eta_draws` was built with `get_halton_normals(S, U, K_w)`.
- `run_mxlogit()`'s advanced workflow (`input_data` with `eta_draws`) passed
  the cube of draws to its compiled routines through RcppArmadillo's
  conversion, which reads a cube stored as integers, rather than doubles, from
  a temporary copy that R may free while the routines read it: such a fit
  could be computed from freed memory, or crash. `get_halton_normals()`
  returns doubles, so only a cube built otherwise (an integer array) was
  exposed. `run_mxlogit()` now converts such a cube to double first, and so
  does `mxl_blp_contraction()` ("Share inversion", above).
- `run_mnlogit()`, `run_mxlogit()` and `run_nestlogit()` did not check the
  values of an `input_data` (the advanced workflow): weights, `X` or `W`
  changed in a prepared object after `prepare_*_data()`, and the values of
  `eta_draws`, reached the compiled routines as they were. In earlier
  results the usual symptom is a fit whose estimates equal its start
  values, reported as converged, with a log-likelihood of `-1e10` (near 0
  with integer64 weights). With the default optimizer (nloptr):
  - A weight that was `NA`, `NaN` or infinite left each model at its start
    values: the optimizer stopped after one iteration reporting success,
    with `NA` or `NaN` standard errors. Since 0.2.0 the only warning, about
    non-uniform weights, pointed elsewhere, and with
    `se_method = "sandwich"` or `"cluster"` there was none. In 0.1.0 the
    multinomial and nested logits stopped instead, with nloptr's error that
    the objective is `NA` at the start values (an infinite weight ended in
    `NLOPT_FAILURE`).
  - Weights of 0 or below were used as given; the preparations reject
    them, as they do non-finite weights, since 0.2.0.
  - integer64 weights were read as their raw bit patterns (a weight of 1 as
    4.9e-324), so the log-likelihood came out near 0, again "success" at
    the start values.
  - A non-finite value in `X` left the multinomial logit at its start values
    in the same way, with warnings that log probabilities had been clamped,
    and ended the nested logit's fit in `NLOPT_FAILURE` with some standard
    errors finite; in 0.1.0 both stopped with nloptr's error that the
    gradient is `NA`. In the mixed logit, a non-finite value in `X`, `W` or
    `eta_draws` gave "success" at the start values with finite standard
    errors and no warning.

  The three functions now check `input_data` before they first evaluate the
  likelihood. Its weights must be finite and strictly positive, with the
  preparations' messages (integer64 weights are converted to double first,
  as the preparations convert them), and its `X` (and `W`) and `eta_draws`
  must be finite, with an error that names the one at fault (for instance
  "`input_data$X` must contain only finite values."). The check reads the
  matrices in place, without a copy. Objects as the preparations return
  them fit as before, bit for bit.
- `prepare_mnl_data()`, `prepare_nl_data()`, `prepare_mxl_data()` and
  `prepare_mnp_data()`, and the `run_*()` functions that call them, stopped
  on an infinite covariate value (`Inf` or `-Inf`) with an internal error
  that did not name the column: "NA/NaN/Inf in foreign function call
  (arg 1)", from the collinearity check's `qr()`, or, for a design of one
  column, "all(is.finite(X)) is not TRUE" (or `W`). They did so in every
  release (`prepare_mnp_data()` since it was added in 0.2.0), while a `NaN`
  in the same place dropped its choice situation as a missing value, and
  `prepare_hmnl_data()` and `prepare_hmnp_data()`, also added in 0.2.0,
  dropped such situations with a warning. The four now drop them too, with
  a warning that names the covariates holding such values ("Removed 2
  choice situations containing non-finite covariate values (in x1, w1).";
  the hierarchical preparations' warning names them now as well), and stop
  when none is left; the error for an `outside_opt_label` that no remaining
  situation offers counts these situations too. A new section of the
  documentation of all six preparations and their fitting functions,
  "Missing and non-finite values", says so, and notes that dropping such
  situations selects the sample on the covariate (an infinite value usually
  comes from a transformation such as `log(0)`). The scan reads each
  covariate column in place first, without its missing values, so data
  without infinite values cost no copy: on 10^7 rows and 10 covariates it
  takes 0.10 s and 40 MB (one logical flag per row), with or without
  missing values, where the hierarchical preparations' former scan took
  0.15 s and 440 MB, and 1.1 GB with missing values. Data without infinite
  covariate values are prepared as before, bit for bit.
- The data preparations (`prepare_mnl_data()`, `prepare_mxl_data()`,
  `prepare_nl_data()`, `prepare_mnp_data()`, `prepare_hmnl_data()`,
  `prepare_hmnp_data()`, and the `run_*()` functions that call them) could
  corrupt or misread an input column whose name matched one of their own.
  They added working columns with fixed names (`alt_int`, `idx_in_group`,
  `HAS_NA`, and in the hierarchical preparations `HB_PERSON`, `task_idx`,
  `TASK_HAS_NA`, `HAS_BAD`, `TASK_HAS_BAD`) to a copy of the data, read
  per-situation results back under names (`N`, `V1`, `chosen`, `pos`) that
  an id column of the same name shadowed, and looked up variables of their
  own (`levels`, `J`, `outside_opt_label`, `ids_to_drop`, `id_col`,
  `choice_col`, and in the hierarchical preparations `person_col` and
  `task_by`) where a column of the same name took their place. `predict()`,
  `logsum()` and `consumer_surplus()` with a data.frame `newdata`,
  `wesml_weights()`, `sample_by_choice()`, and the hierarchical fits and
  their post-estimation methods did the same: they looked up variables of
  their own (`am`, `pos`, `spec`, `d`, `choice_col`, `id_col`, `keep_ids`)
  inside the data or the fit's `alt_mapping`, and read counts and strata
  back under names (`V1`, `.strat`) that an id column of the same name
  shadowed. These clashes changed the results without an error:
  - in `prepare_hmnl_data()` / `prepare_hmnp_data()`, a covariate or
    control-function residual named `alt_int`, `idx_in_group` or
    `task_idx`, or one named `HB_PERSON` when the respondent ids were
    numeric, was replaced by the alternative codes, the within-situation
    positions, the task indices or the respondent ids (values -1.2, -0.69,
    2.29 entered the design matrix as 1, 2, 3), as was an alternative-level
    covariate named `idx_in_group` when every choice situation offered the
    same alternatives; `prepare_mnp_data()` differenced the codes of a
    covariate named `alt_int`;
  - in the same two functions, a column named `task_by`, other than the
    task column of a fit without `person_col`, regrouped the rows by its own
    values: a covariate or control-function residual whose values differ
    from row to row made every row a choice situation of its own, without
    an error whenever the outside option was modelled, as it is by default
    (repeated values usually stopped the preparation);
  - in `prepare_mnl_data()`, `prepare_mxl_data()` and `prepare_nl_data()`,
    a choice column also listed in `covariate_cols` left the working table
    before the choices were counted: the preparation stopped with an
    unrelated error or, when an object of the column's name was visible to
    the call (a vector `choice` in the workspace, say), counted that object
    instead, so that `alt_mapping` credited every alternative with every
    choice and, unless the rows were already in prepared order, the chosen
    alternatives were wrong;
  - a choice column named `alt_int` or `idx_in_group` recorded the first
    alternative of every choice situation as chosen (`idx_in_group` was
    harmless in `prepare_mnp_data()`), and in the hierarchical preparations
    one named `task_idx` recorded the outside good in every task but the
    first;
  - an id column named `alt_int` or `idx_in_group` with the outside option,
    and in the hierarchical preparations a task column named `alt_int` or
    `task_idx`, regrouped the rows into spurious choice situations;
  - with the outside option, any id, choice, covariate, weight, cluster or
    decision-maker column named `outside_opt_label` kept the outside-option
    rows as an extra inside alternative;
  - an alternative column named `idx_in_group` with unequal choice sets
    added unidentified constants, leaving every standard error `NA`;
  - in `predict()`, `logsum()` and `consumer_surplus()` with a data.frame
    `newdata`, a covariate named `pos` whose values were valid positions in
    the fit's `alt_mapping` (without an outside option, any permutation of 1
    to J within each choice situation, such as a rank) gave each row the
    alternative at that position, so the alternative-specific constants and
    the shares went to the wrong alternatives;
  - a column named `keep_ids` made `sample_by_choice()` return the wrong
    choice situations (the whole population when it was the id column).

  These changed labels and summaries: an alternative column named
  `alt_int`, or with equal choice sets `idx_in_group`, had its labels
  replaced by the codes; one named `TAKE_RATE` or `MKT_SHARE` had its
  labels, and so the ASC names, overwritten by the take rates or shares;
  one named `N_OBS` or `N_CHOICES` with numeric labels misstated the take
  rates or shares, and made `gof(null = "market_shares")` stop (`N_OBS`)
  or return a wrong null log-likelihood (`N_CHOICES`); and one named `J`
  mislabelled the constants of `run_mnlogit()`, `run_mxlogit()`,
  `run_nestlogit()` and `prepare_mnp_data()` (with three alternatives
  labelled 1 to 3 or stored as a factor, the third alternative's constant
  took the base alternative's label). Whatever matched alternatives by
  label inherited the error: `prepare_nl_data()` assigns nests by label, so
  it stopped or, with text labels such as "1" to "10", assigned wrong
  nests; `predict()`, `logsum()` and `consumer_surplus()` with `newdata`
  stopped or attached another alternative's constant (labels 0 to 2 had
  become 1 to 3); and `predict()` on the hierarchical fits, which takes an
  unknown label for a new alternative, predicted every alternative from
  posterior-predictive constants without an error.

  An id column named `N` gave the ids as the choice-set sizes `M`, in the
  preparations and in `predict()`, `logsum()` and `consumer_surplus()` with
  `newdata`. The kernels' checks usually stopped the fit or prediction, but
  a prediction whose ids summed to the number of rows (ids 1 to 11 with 6
  alternatives each, for instance) ran on wrong choice sets without an
  error, and more rarely so could a fit. A covariate named `levels` in
  `prepare_mnl_data()`, `prepare_mxl_data()` and `prepare_nl_data()` left
  every alternative code `NA` (or, when its values repeated, stopped the
  preparation), and with the outside option an id column named `pos` gave
  the ids as the chosen positions; fitting then failed.
  Other clashes stopped with unrelated errors: an id column named `V1` (the
  name `data.table::fread(header = FALSE)` and `as.data.frame()` of an
  unnamed matrix give the first column) in `prepare_mnl_data()`,
  `prepare_mxl_data()` and `prepare_nl_data()` whenever a weight or cluster
  column was given, and always in the hierarchical preparations; a
  covariate named `alt_int` or `idx_in_group` in `prepare_mnl_data()`,
  `prepare_mxl_data()` and `prepare_nl_data()`; a column named
  `choice_col` that a model used (other than the nested logit's nest
  column), as did one named `id_col` outside the hierarchical preparations
  and, in those, one named `person_col` when `person_col` was given or
  `id_col` when it was not; an id or alternative
  column whose name contains a comma, which data.table's grouping splits
  into several names; and an id, weight or cluster column named like a
  working column, among others. So did, outside the preparations, a
  covariate named `am` in a data.frame `newdata`, most covariates there
  named `pos`, and with an outside option one named `spec`; an id column
  named `V1` or `.strat`, or a choice column named `choice_col`, in
  `wesml_weights()` and `sample_by_choice()`; and an alternative column
  named `am`, which stopped `run_hmnlogit()` and `run_hmnprobit()`, or
  `d`, which stopped their `predict()` without `newdata`, `elasticities()`,
  `diversion_ratios()` and `ppc_shares()`, and for the HMNL `logsum()`
  without `newdata` and `consumer_surplus()`. Fits whose data used a name
  in the first group, and fits with an id column named `N` that ran, should
  be re-estimated, and the predictions and samples in that group
  recomputed; for the second group, re-estimate nested logit fits and
  recompute predictions made with `newdata`.

  Covariates are now read from the data rather than from the working table;
  the remaining working columns, per-situation results, counts and strata
  carry the reserved prefix `.choicer_` or are read by position; and the
  preparations, post-estimation and the WESML helpers compute their own
  values, and read the columns they are given by name, before indexing the
  data or `alt_mapping`, so the names listed above are no longer looked up
  among the columns. Inputs without a clash give the same results as
  before. A choice column stored as text or as a factor, which the
  preparations, `wesml_weights()` and `sample_by_choice()` never accepted,
  and in `prepare_mnl_data()`, `prepare_mxl_data()` and `prepare_nl_data()`
  a covariate that is also the id, alternative, choice, weight, cluster or
  decision-maker column (see above for the choice column), now stop them
  with an error that says why. Two inputs are now errors: a column whose
  name starts with `.choicer_`, used as an id, alternative, choice,
  covariate, weight, cluster, nest or decision-maker column; and an
  alternative column named `alt_int`, `N_OBS`, `N_CHOICES`, `TAKE_RATE` or
  `MKT_SHARE`, the fixed columns of the returned `alt_mapping`, which must
  now be renamed. Released versions handled `.choicer_` columns correctly
  (they used no such names) and, of the alternative-column names, only an
  `alt_int` column already coded 1 to J, without the outside option.
- A covariate named like a parameter the fit generates gave two parameters
  one name. `prepare_mnp_data()` (and so `run_mnprobit()`) names the
  constant of each non-base alternative `ASC_<label>` and builds `param_map`
  by column name: since 0.2.0, a covariate such as `ASC_b` next to
  alternative `b`'s constant left that constant in neither the `beta` nor
  the `asc` block, so `recovery_table()` could not align the ASCs, and
  `summary()` stopped with a duplicate row-name error; when the covariate
  was itself alternative `b`'s dummy, one of the two columns was dropped as
  collinear without the usual message. `run_mnlogit()`, `run_mxlogit()` and
  `run_nestlogit()` generate `ASC_<label>`, `Lambda_<k>` (nested logit),
  `Mu_<variable>` and `L_<i><j>` (mixed logit), and the mixed logit's
  `summary()` prints `Sigma_<i><j>` for `L_<i><j>` and `exp(Mu_<variable>)`
  for the mean of a log-normal coefficient; a covariate could take any of
  these names. Estimation addresses their parameters by position, so
  without named bounds (below) the estimates, standard errors and
  predictions were right. But `summary()` stopped with a duplicate
  row-name error, unless it relabelled the generated parameter
  (`L_<i><j>`, or a log-normal `Mu_<variable>`),
  and indexing `coef()` or `vcov()` by name, or `wtp(attr_vars = )`, found
  the covariate: `wtp(attr_vars = "ASC_b")` returned the covariate's WTP,
  and the constant's could not be reached by name. A named bound in
  `run_mxlogit()` goes to the first parameter of that name, and the fixed
  coefficients come first, so it bounded the covariate's coefficient and
  left the generated parameter free. With a covariate named `ASC_b`,
  `lower = c(ASC_b = 0)` held that covariate's coefficient at 0 whenever its
  unconstrained estimate was negative, and so moved the other estimates,
  while alternative `b`'s constant stayed unbounded; a covariate named
  `L_11` would likewise have taken the documented `lower = c(L_11 = -5)`
  meant for the first Cholesky diagonal. Such names are now an error that
  lists them: `prepare_mnp_data()` stops before building `X` (rename the
  covariates, or set `use_asc = FALSE` for ASC dummies built by hand), and
  the three other functions stop before optimizing, `run_mxlogit()` before
  building its draws. The same check stops an `input_data` whose `X`
  repeats a column name, and repeated `run_nestlogit(param_names = )`, which
  used to fit with the same ambiguities. Preparations and fits whose
  parameter names and `summary()` labels are all distinct are unchanged,
  `summary()` included.
- Columns of class `integer64` (bit64), which `data.table::fread()` returns
  for integers beyond 2^31 - 1 and several database drivers return for
  `BIGINT`, passed the check that covariates are numeric but entered the
  design matrices as their raw 64-bit patterns read as doubles, in every
  model of released versions. Positive values became tiny numbers (1 read as
  4.9e-324, 2^40 as 5.4e-312), so fits, standard errors and predictions were
  silently wrong. Negative values above -2^52 became `NaN`, which stopped
  the preparation with an unrelated error (`NA/NaN/Inf in foreign function
  call`); with bit64 loaded, `prepare_mnp_data()` differenced these columns
  as integers and so usually stopped the same way. When bit64 was not
  loaded, as for data read from an `.rds` file, the check for missing
  values misread these columns too: it took negative values for missing
  ones, dropping their choice situations with the warning about missing
  values, and let a missing value into the design as 0. The design
  builders now read integer64 columns as their values, as
  `as.double()` converts them: `prepare_mnl_data()`, `prepare_mxl_data()`,
  `prepare_nl_data()`, `prepare_mnp_data()`, `prepare_hmnl_data()` and
  `prepare_hmnp_data()` (including `alt_covariate_cols` and
  `cf_residual_col`), and the `newdata` of `predict()`, `logsum()` and
  `consumer_surplus()`. The conversion is exact below 2^53 in magnitude;
  from 2^53 up, values are rounded to the nearest double, with a warning
  naming the column. An integer64 column among those a model uses now loads
  bit64 first, so its missing values are found and handled as in other
  columns; if bit64 is not installed, such a column is an error. Weights
  (`weights_col`, a `weights` vector) and integer64 `X` and `W` matrices in
  the list form of `newdata` had the same defect and are now converted the
  same way, as are prediction `weights`. bit64 is now a suggested package.
- `sample_by_choice()` misread `integer64` ids, as `data.table::fread()`
  returns ids beyond 2^31 - 1, ever since it was added in 0.2.0: it matched
  the sampled ids back to the rows after `unlist()` had dropped their
  class. With bit64 loaded, the sample came back empty without an error;
  without it, each id drawn from among ids of large magnitude also selected
  its neighbours (about six choice situations for each one requested, for
  ids near 10^18). With
  the outside option, `wesml_weights()` returned its ids as raw bits
  (doubles such as 7.8e-242) whenever a choice situation chose the outside
  good and bit64 was not loaded, and with `attach = TRUE` stopped with an
  unrelated error. Both functions now load bit64 when the id, alternative
  or choice column is `integer64`, as the preparations do, and keep the
  ids' class, so they sample and weight such ids by value; if bit64 is not
  installed, such a column is an error. Samples of other ids are drawn as
  before.
- `nl_blp_contraction()` takes `lambda` and `nest_idx` as given, not from a
  fitted model, and never checked them against each other: a nest code
  below 1 or missing, or a `lambda` with other than one entry per nest,
  failed an Armadillo index or size check inside the contraction's parallel
  loop, which aborted R in builds with OpenMP and was an Armadillo error in
  builds without it. Both are now errors before the contraction starts. As
  in the other nested logit kernels, this includes a nest code below 1 on
  an alternative that no choice situation offers, which used to be
  ignored.
- `run_mnprobit()`, `run_hmnlogit()` and `run_hmnprobit()` miscounted the
  draws to keep when `R - burn + thin - 1` exceeded `.Machine$integer.max`,
  as for a `thin` near it: the count overflowed an integer. Where it came
  out 0 this aborted R in builds with OpenMP; otherwise the run stopped
  with an unrelated Armadillo error (with `keep_beta_i = "draws"`, the
  memory estimate failed first, with another unrelated error). The count,
  `ceiling((R - burn) / thin)`, is now formed without overflow.
- Mixed logit fits with stored draws (`draws = "store"`, the default) read
  some units' draws from the wrong place, without an error, once the draw
  cube (`K_w * S * N` values) was large enough that a unit's draws began
  past its first 2^32 values. The cube was built for any `S * N` up to
  2^31 - 1 points (past that, building it failed). With three random
  coefficients, the default 100 draws and 15 million choice situations, for
  instance, it holds 36 GB, and the last 683,442 situations read other
  units' draws in place of their own. Those draws entered the estimates and
  their standard errors, and every result that rebuilds the cube: standard
  errors recomputed with `vcov(type = )` or `wesml_vcov()`, predictions,
  elasticities, diversion ratios, BLP inversions, log-sums, consumer surplus
  and fit statistics. Such results should be recomputed. The compiled code
  used RcppArmadillo's default 32-bit word, under which the kernels'
  in-place view of a matrix or cube of more than 2^32 - 1 values wrapped its
  element count and element offsets: a unit's draws, or a choice situation's
  covariate values, that began past the first 2^32 values were read from the
  wrong place. A hand-built `eta_draws` that large was misread the same way.
  A design matrix that large (for instance 10^9 rows and 5 covariates)
  reached the kernels only through a hand-built `input_data`, the `newdata`
  of a prediction or welfare function, or a direct call to an exported
  contraction function: the preparations stopped it earlier, at `qr()`'s
  limit of 2^31 - 1 values. The compiled code now uses Armadillo's 64-bit
  word (`ARMA_64BIT_WORD`), so the kernels read matrices and cubes of any
  size R can hold, and outputs of more than 2^32 - 1 values, such as the
  output of `conditional_tastes()`, are allocated where Armadillo stopped
  with "requested size is too large". Below 2^32 values, this correction by
  itself leaves results unchanged, bit for bit.
- `summary()` of a mixed logit fit without a variance, one whose information
  matrix is singular, printed the Cholesky factor's parameters (the diagonal
  on its log scale) under the labels of the covariance matrix,
  `Sigma_<i><j>`, and the means of log-normal coefficients untransformed
  under `exp(Mu_<variable>)`: it transformed the estimates only when it had a
  variance for their standard errors. It now reports the covariance matrix and
  `exp(Mu)` whether or not the variance is available, with `NA` standard
  errors where it is not.
- `run_mxlogit(use_asc = FALSE)`, the mixed logit without alternative-specific
  constants, counted the `J - 1` constants of the model with them whatever
  `use_asc` said. Since 0.2.0 the first evaluation of the likelihood, at start
  values that many entries too long, stopped with "Theta vector too long: ...
  Did you mean use_asc = TRUE?". 0.1.0, whose kernels did not check the
  length, returned the constants as coefficients held at their start of 0, on
  which the likelihood does not depend: its other estimates were those of the
  model without constants, but `n_params`, AIC and BIC counted the constants,
  and every standard error was `NA` (the information is zero in their rows and
  columns, so singular). It now fits the model without constants: its
  parameters are the fixed coefficients, the random coefficients' means (with
  `rc_mean = TRUE`) and their Cholesky factor, `param_map` has no `asc` block,
  as in `run_mnlogit()` and `run_nestlogit()`, and the default start,
  `theta_init` and `lower`/`upper` have one entry per parameter, with or
  without the outside option. Prediction, elasticities, diversion ratios, BLP
  inversions, log-sums, consumer surplus, WTP, conditional tastes and every
  `vcov(type = )` work on such fits. Fits with constants (`use_asc = TRUE`,
  the default) are unchanged, bit for bit.
- `blp()` on a model fitted without constants (`use_asc = FALSE`) and with
  the outside option (`include_outside_option = TRUE`) started its
  contraction from one zero too many, counting the outside option, which
  has no constant, among the alternatives. The multinomial and nested logit
  methods stopped with an error that blamed `target_shares` ("must have the
  same length as the total number of alternatives") ever since `blp()`
  gained them (the multinomial logit in 0.1.0, the nested logit in 0.2.0);
  the mixed logit's had the same start, which no fit of a released version
  reached. The contraction now starts from one zero per inside alternative;
  other fits are unchanged.
- `predict(type = "shares")` with `newdata` that does not offer the
  alternatives that come last in `object$alt_mapping`, as when an
  alternative exits in a counterfactual, returned one share fewer per such
  alternative for multinomial and nested logit fits, with or without
  constants, since `newdata` was added in 0.2.0; mixed logit fits without
  constants, possible since the `use_asc = FALSE` correction above, did the
  same. The shares were right, but the result stopped at the last
  alternative that `newdata` offers, so it no longer lined up with the
  baseline's shares; a mixed logit fit with constants, and any model when an
  earlier alternative is missing, gave each alternative that `newdata` does
  not offer a 0. Multinomial, nested and mixed logit fits now return one
  share per alternative of `object$alt_mapping`, in its order (the outside
  option's first, when it is modelled), 0 for an alternative that `newdata`
  does not offer; other predictions are unchanged, bit for bit. The
  documentation called the result a named numeric vector; it is, as it
  always was, a one-column matrix without names (the labels are in
  `object$alt_mapping[[2]]`).
- An `outside_opt_label` that names none of the alternatives, without the
  outside option (`include_outside_option = FALSE`), became a phantom base:
  the preparations code that label first, so it took code 1 with no rows.
  Fits with constants then stopped with an error about the constants' block
  (with three alternatives, "Theta's delta (ASC) block implies 3
  alternatives but alt_idx references alternative 4"; an Armadillo index
  error in 0.1.0), nested logit fits without them with an error about the
  nests (in 0.1.0 R aborted), and multinomial logit fits without constants,
  like mixed logit fits without them since the `use_asc = FALSE` correction
  above, carried an alternative that is never offered: a predicted share of
  0, and elasticity and diversion-ratio matrices with an extra row and
  column and no labels. An `NA` label silently made the first alternative
  the base. A label of another class than the alternative column (a
  character label with a column of dates, a plain number with an
  `integer64` column, a factor label with a character column) could leave
  alternatives without a code, which the preparations returned: with dates
  and a character label, checked on Apple silicon, fits since 0.2.0 stopped
  with an error about the alternative codes, and 0.1.0 fitted a model
  without constants on them. `prepare_mnl_data()`, `prepare_nl_data()` and
  `prepare_mxl_data()`, and the `run_*()` functions that call them, now stop
  in these cases with an error that says why. Labels that coded the
  alternatives before still do: the check reads the codes the preparations
  build.

# choicer 0.2.1

Patch release addressing a compilation warning reported by CRAN's GCC check
flavors. No user-visible behavior, API, or numerical results change.

## Compilation

- Replaced the seven `#pragma omp master` directives in the Gibbs samplers and
  thread-introspection helpers with a `CHOICER_OMP_MASKED` macro that emits
  `masked` on OpenMP 5.1 and later, `master` on earlier versions, and nothing
  when the package is built without OpenMP. OpenMP 5.1 deprecated `master` in
  favor of `masked`, and GCC 16 — now the compiler on the
  `r-devel-linux-x86_64-fedora-gcc` and `r-devel-linux-x86_64-debian-gcc`
  check flavors — warns on the old spelling, which `R CMD check` reports as a
  significant install warning, escalating the check result to WARNING.
  `masked` with no `filter` clause is semantically identical to `master`: the
  block runs on the primary thread with no implied barrier, so the samplers'
  fixed-order reductions, barrier structure, and
  no-R-API-off-the-primary-thread contract are unchanged, as are posterior
  draws for a given seed and thread count.

# choicer 0.2.0

## Public API cleanup (breaking)

- Low-level C++ likelihood, gradient, Hessian, score, prediction,
  post-estimation, and Gibbs-sampler wrappers are now internal implementation
  details rather than exported package functions. High-level fitting,
  preparation, S3 post-estimation, simulation, and recovery APIs are unchanged.
  The expert-facing raw BLP contractions, `get_halton_normals()`,
  `set_num_threads()`, and `thread_info()` remain public. With no downstream
  CRAN dependencies, v0.2.0 is the least disruptive point to narrow this
  surface before applications depend on unstable kernel signatures.

## Robust and clustered standard errors (MNL / MXL / NL)

- `vcov()` on a fitted MNL, MXL, or NL model gains `type =` and `cluster =`
  arguments for post-hoc variance recomputation without refitting (requires
  `keep_data = TRUE`): `"hessian"` (inverse analytical Hessian), `"bhhh"`
  (OPG), `"robust"` (Huber-White / WESML sandwich), and `"cluster"`
  (cluster-robust sandwich over within-cluster sums of weighted scores — use
  it when the same decision maker contributes several choice situations).
  With no arguments, `vcov()` returns the as-fitted variance, unchanged.
- New `cluster_col=` argument on `run_mnlogit()` / `run_mxlogit()` /
  `run_nestlogit()` (and the corresponding `prepare_*_data()` functions):
  supplies per-situation cluster labels at fit time and selects the new
  `se_method = "cluster"`.
- All score-based variances (BHHH, robust, cluster) are now assembled in R
  from one per-situation score matrix, computed by new internal C++ kernels
  that reuse the BHHH loop bodies. The existing `se_method = "sandwich"` /
  `wesml_vcov()` results are unchanged (the robust meat is the `w^2` special
  case of the shared path).
- Scope note (MXL): clustering repairs the inference, not the estimand. The
  MXL simulated likelihood treats each choice situation as an independent
  draw from the mixing distribution (cross-sectional MSL, not the panel
  product form), so clustered standard errors on panel data are robust to
  within-person dependence but do not turn the fit into a panel mixed logit.
  For panel random coefficients use `run_hmnlogit()` (`person_col=`).

## Corrections

- `blp()` for an MXL fit using `draws = "generate"` now carries the fitted
  draw count, seed, and digit-permutation setting through every contraction
  evaluation. Previously that path could evaluate counterfactual shares with
  fallback draw metadata rather than the simulation design used for the fit.
- On-the-fly MXL draws now accept all 128 prime-base dimensions implemented by
  the generator (the prior guard rejected the 128th dimension); 129 or more
  random coefficients still fail early with an explicit limit.
- HMNL prediction and welfare now use a taskwise max-shifted softmax/logsum, so
  extreme counterfactual utilities remain finite. Counterfactual compensating
  variation accepts explicit non-negative task weights and matches baseline
  and policy choice situations by (`person_col`, `id_col`) before subtraction;
  it now rejects added, dropped, or substituted tasks instead of pairing them
  silently by sorted position.
- The documentation previously claimed that MCMC draws from the Gibbs
  samplers (`run_mnprobit()`, `run_hmnlogit()`, `run_hmnprobit()`) are
  bitwise reproducible regardless of the OpenMP thread count. Independent
  validation disconfirmed this: draws are reproducible given the seed and a
  fixed thread count, and across different thread counts are invariant only
  up to floating-point reduction-order round-off (~1e-15), not bitwise. All
  man pages, vignettes, and NEWS entries now state the correct guarantee.
- Documentation now distinguishes the frequentist models' analytical
  derivatives from the Bayesian C++ samplers, labels `run_mxlogit()` as a
  cross-sectional simulated likelihood, and no longer directs users to
  WTP-space estimation or bounded/censored mixing distributions as if choicer
  implemented them. The README also labels `mode_choice` shares and welfare as
  sample-design quantities unless external population shares and WESML weights
  are supplied.
- The MNP tutorial no longer recommends unavailable user-specified starts and
  now states its current post-estimation boundary. HMNP documentation now makes
  explicit that it uses iid utility-level normal shocks rather than the full
  differenced-error covariance estimated by `run_mnprobit()`.

## Hierarchical Bayes (HMNL/HMNP) convergence diagnostics and multi-chain performance

- New `ess(draws)` (rank-normalized bulk and tail effective sample size) and
  `mcse(draws, kind = c("mean", "median"))` (Monte Carlo standard error of
  the mean or median), following Vehtari, Gelman, Simpson, Carpenter &
  Bürkner (2021, *Bayesian Analysis*). `rhat()` gains a `rank = TRUE` option
  for the same paper's rank-normalized, folded R-hat; the default
  `rank = FALSE` reproduces the original split R-hat exactly.
- New `traceplot()` generic, with a `choicer_hb` method, for visual
  inspection of the `b`/`theta`/`sigma_d2` chains (and user-selected `delta`
  columns).
- `summary()` on a `choicer_hmnl`/`choicer_hmnp` fit now prints one
  consolidated multi-chain diagnostics table (R-hat, ESS bulk, ESS tail, and
  MCSE per parameter block, plus a single worst-case summary row spanning
  all `J` `delta_j` alternative effects) in place of the previous bare
  split-R-hat printout. HMNP now always prints an explicit
  "conjugate — no acceptance step" line where HMNL prints its beta/delta
  acceptance rates. The fit-time convergence warning now checks rank-R-hat
  and ESS bulk across every tracked parameter, including all `J` `delta_j`
  columns (previously only `b`/`theta`/`sigma_d2` were checked).
- `run_hmnlogit()` / `run_hmnprobit()` fits now retain **all** requested
  chains' hierarchical draws in a new `object$chains` field (a list, one
  element per chain: `b`, `w_vech`, `delta`, `theta`, `sigma_d2`, plus
  `loglik` for HMNL / `sigma2` for HMNP), which the multi-chain diagnostics
  above consume. `object$draws` (chain 1) and `object$beta_i`
  (chain-1-only respondent taste summaries/draws) are unchanged.
- Chains requested via `chains=` are sampled sequentially, and all of their
  hierarchical draws are retained (`object$chains`) for the multi-chain
  diagnostics above. Parallel multi-chain execution is planned for a future
  release.
- The `keep_beta_i = "draws"` memory guard on both models is now based on a
  measured ~1.9x `Rcpp::wrap()`/R-list overhead factor (the previous formula
  under-estimated true memory use by roughly 2x) and accounts for all
  requested chains; it still fails fast before the chain runs rather than
  after an out-of-memory blow-up.

## Hierarchical Bayesian MNL and MNP (BLP-style random-effects ASCs)

- New `run_hmnlogit()` (class `choicer_hmnl`) and `run_hmnprobit()` (class
  `choicer_hmnp`): hierarchical Bayes discrete choice with two decoupled
  random-effect levels — respondent-level structural tastes
  `beta_i ~ N(b, W)` (`W` is `K x K`, no ASC columns) and a global BLP-style
  alternative effect `delta_j = z_j' theta + xi_j`, `xi_j ~ N(0, sigma_d^2)`
  (one scalar variance regardless of `J`; partial pooling toward the
  characteristics-based mean; posterior-predictive `delta` for alternatives
  outside the estimation sample). Both models carry a first-class implicit
  outside option (systematic utility 0 plus its own shock) that anchors the
  location of `delta` — no base alternative or sum-to-zero constraint.
  Panel and cross-sectional (`person_col = NULL`, `T_i = 1`) modes share
  one code path.
- HMNL: adaptive RW-Metropolis-within-Gibbs on the mnprobit engine pattern
  (single persistent OpenMP region, per-(iteration, unit) RNG streams;
  draws are reproducible given the seed and a fixed thread count, and
  invariant across thread counts only up to floating-point reduction-order
  round-off, ~1e-15), with a strictly serial `delta` sweep
  (the conditionals are coupled through the softmax denominators) at O(1)
  incremental cost per affected task; supports log-normal coordinates via
  `rc_dist`. HMNP: fully conjugate Albert-Chib augmentation in
  un-differenced utility space with iid `N(0, sigma^2)` shocks, a
  non-identified `sigma^2` chain (parameter expansion), and per-draw scale
  normalization of every identified quantity.
- Priors: `b ~ N(b_bar, A^-1)`, `W ~ IW(nu, V)`,
  `theta ~ N(theta_bar, A_theta^-1)`, and `sigma_d ~ half-Cauchy(0, s_d)`
  via the Makalic-Schmidt conjugate scale mixture (IG fallback available).
- New preps `prepare_hmnl_data()` / `prepare_hmnp_data()` (two-level
  person/task indexing, implicit-outside convention shared with
  `prepare_mnl_data()`, alternative-level design `Z`, optional
  control-function residual column for price endogeneity per Petrin &
  Train 2010) and DGPs `simulate_hmnl_data()` / `simulate_hmnp_data()`.
- Post-estimation on the shared `choicer_hb` class, all posterior-draw
  based: `predict()` (population/individual, entry counterfactuals via the
  posterior-predictive `delta`; HMNP probabilities by deterministic 20-node
  1-D Gauss-Hermite approximation), `wtp()` (posterior median + quantile
  intervals), `logsum()` / `consumer_surplus()` (HMNL-only; probit Emax is
  roadmapped), `elasticities()` / `diversion_ratios()` (common-random-path
  perturbation engine including the outside option), `recovery_table()`,
  plus `rhat()` (split R-hat) and `ppc_shares()` diagnostics.
- Two new math articles: `vignettes/articles/hierarchical_mnl_math.Rmd`
  and `hierarchical_mnp_math.Rmd`; recovery demos in
  `inst/simulations/hmnl_simulation.R` / `hmnp_simulation.R`.

## WESML sandwich inference for MNL and nested logit

- `run_mnlogit()` and `run_nestlogit()` now compute the robust (Huber–White /
  WESML) sandwich variance `V = A^{-1} B A^{-1}`, at full feature parity with
  `run_mxlogit()`. Both gain `se_method = "sandwich"` (bread = weighted negated
  Hessian, meat = weight-squared OPG) and `se_method = "bhhh"` (ordinary
  outer-product-of-gradients). `run_mnlogit()`'s `se_method` therefore now
  accepts `"hessian"` (default), `"bhhh"`, or `"sandwich"`; `run_nestlogit()`'s
  accepts `"hessian"` (default), `"numeric"`, `"bhhh"`, or `"sandwich"`.
- New C++ kernels `mnl_bhhh_parallel()` and `nl_bhhh_parallel()` accumulate the
  weighted outer product of per-individual scores. Their per-individual score is
  weight-free, so passing `weights = w` yields the BHHH information and
  `weights = w^2` yields the sandwich meat. The NL kernel includes the full
  beta/lambda/delta score blocks (singleton-nest lambdas fixed to 1 contribute
  no score).
- `run_mnlogit()`, `run_nestlogit()`, `prepare_mnl_data()` and
  `prepare_nl_data()` gain a `weights_col` argument: a row-level weight column
  is collapsed to one weight per choice situation (validated constant within
  `id`). A choice-based-sampling provenance guard auto-adopts the recorded
  WESML weight column and errors rather than silently fitting unweighted under a
  WESML label; non-uniform weights under a non-sandwich `se_method` emit a
  warning.
- Behavior change: weighted fits now emit a steering warning recommending
  `se_method = "sandwich"` when non-uniform weights are supplied with the default
  (`"hessian"`) or `"bhhh"` method; the `"bhhh"` case gets a sharper message
  explaining that BHHH/OPG is not a valid WESML correction (its meat is `w^1`,
  not `w^2`). Point estimates and standard errors are unchanged.
- `wesml_vcov()` now dispatches on `choicer_mnl` and `choicer_nl` (in addition
  to `choicer_mxl`), returning the post-hoc sandwich variance from a fit stored
  with `keep_data = TRUE`.
- `summary()` for MNL and NL fits now reports the standard-error method and any
  WESML weighting in the printed footer.
- Added "Choice-Based Sampling and WESML Weighting" sections to the multinomial
  logit and nested logit derivation vignettes.

## Weighting safety hardening

- Weights are now validated to be finite and strictly positive in
  `prepare_mnl_data()`, `prepare_nl_data()`, and `prepare_mxl_data()` (covering
  both the `weights=` and `weights_col=` paths). Zero, negative, or non-finite
  weights previously could slip through and silently invalidate weighted and
  WESML sandwich inference (weight `w` enters the bread, `w^2` the meat); they
  now error with an actionable message.
- Advanced-mode fits (`input_data=` passed directly to `run_mnlogit()`,
  `run_nestlogit()`, or `run_mxlogit()`) that carry WESML `choice_sampling`
  provenance but resolve to uniform weights now error instead of warning. The
  message explains how to proceed: bake the non-uniform WESML weights into
  `input_data` via `prepare_*_data(weights=/weights_col=)`, or strip the
  provenance with `attr(input_data, "choice_sampling") <- NULL` for a deliberate
  unweighted fit. Convenience-mode behavior is unchanged.

## Documentation

- Added eight vignettes: a getting-started tour ("Discrete choice from data to
  policy, in a dozen lines"), one per model (multinomial logit, mixed logit,
  nested logit, Bayesian multinomial probit, hierarchical Bayes), one on
  choice-based sampling and WESML weights, and one on standard errors
  (Hessian / BHHH / robust / cluster-robust, and when to use which).
- Added the `mode_choice` data set: the classic Greene & Hensher intercity
  travel-mode choice data (210 travellers x 4 modes), in choicer's long layout,
  used by the getting-started vignette. `?mode_choice` and the vignettes now
  document that the sample is choice-based (car under-sampled): slopes and WTP
  ratios are unaffected in the ASC-saturated logit, while constants, shares,
  and surplus levels inherit the design — see the WESML vignette for the
  correction.
- Added a pkgdown website configuration, including the model derivation notes as
  "The math behind choicer" articles.
- Expanded the identification discussions across the documentation: price
  endogeneity and the control-function route in the getting-started
  model-choice section (Petrin & Train 2010; `blp()` as the Berry 1994
  inversion), what identifies the nested-logit dissimilarity parameters, the
  Keane (1992) covariance-identification caveat for the multinomial probit,
  and the classical MSL asymptotics (fixed-`S` bias, `S` growth conditions) in
  the mixed logit math note.

## Mixed logit — on-the-fly digit-permuted Halton draws

- `run_mxlogit()` gains three new arguments: `draws`, `seed`, and `scramble`.
  - `draws = "store"` (default) keeps the existing behavior: a full
    K_w × S × N Halton cube is pre-materialized and stored in memory.
  - `draws = "generate"` activates the new mode: each individual's S draws are
    computed on the fly in C++ from a compact seed, eliminating the O(N)
    `eta_draws` cube. This is the recommended choice when N is large or memory
    is constrained.
  - `scramble = "permuted"` (default when `draws = "generate"`) applies one
    seeded permutation per dimension and base-digit position, shared across
    sequence indices. This is a deterministic position-wise digit permutation,
    not Owen's nested-uniform scramble, and no standard randomized-QMC
    unbiasedness or replicate-error guarantee is claimed. The historical value
    `"owen"` remains as a deprecated alias for compatibility.
    `scramble = "none"` reproduces the randtoolbox sequence exactly.
  - `seed` sets the integer master seed for the on-the-fly generator; if
    `NULL` (default), a seed is drawn from R's RNG so `set.seed()` governs
    reproducibility. Ignored when `draws = "store"`.
- The `draws_info` field on `choicer_mxl` objects gains three new elements:
  `mode` (`"store"` or `"generate"`), `seed`, and `scramble`. Existing code
  that accesses `draws_info$S`, `draws_info$N`, or `draws_info$K_w` is
  unaffected; the new fields are `NULL` for objects fitted before this release.
- All post-estimation generics (`predict()`, `elasticities()`,
  `diversion_ratios()`, `logsum()`, `consumer_surplus()`, `vcov()`) propagate
  the fitted draw mode automatically — no user action required.
- Default behavior (`draws = "store"`) is bitwise unchanged.


## Bayesian models

- `run_mnprobit()` — Bayesian multinomial probit via Gibbs sampling with data augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). Runs the non-identified chain with conjugate priors and reports identified quantities normalized per draw by `sigma_11`. New `choicer_mnp` posterior object with `summary()` (posterior mean, SD, credible intervals), `coef()`, `vcov()`, `nobs()`; math note in `vignettes/articles/bayesian_multinomial_probit_math.Rmd`
- C++ MCMC infrastructure built from scratch (`src/rng.h`, `src/bayes_samplers.h`): xoshiro256++/splitmix64 RNG with one stream per (iteration, observation) — draws are reproducible given the seed and a fixed OpenMP thread count (invariant across thread counts only up to floating-point reduction-order round-off, ~1e-15) — plus exact truncated-normal, multivariate-normal, Wishart (Bartlett), and inverse-Wishart samplers. The truncated normal picks the cheapest exact method per region (naive normal rejection in high-mass regions, Robert 1995 exponential rejection in the tail, inverse CDF for narrow intervals)
- The Gibbs chain runs inside a single persistent OpenMP region: the latent-utility sweep and mean refresh are work-shared across choice situations, the conjugate beta/Sigma conditionals run on the master thread between lightweight barriers with hand-rolled fixed-order linear algebra (no BLAS inside the region), and the truncated-normal conditional moments are hoisted per Sigma draw. On the `_benchmarks/` MNP preset this samples ~2x faster than `MNP::mnp()` and ~3.5x faster than `bayesm::rmnpGibbs()` on one thread, and scales with threads on top
- `simulate_mnp_data()` — probit DGP returning a `choicer_sim` with truth on the identified scale; `recovery_table()` gains a `choicer_mnp` method (posterior mean/SD, normal-approximation credible intervals, `sigma` block); parameter-recovery walkthrough in `inst/simulations/mnp_simulation.R`

## Post-estimation

- `wtp()` — willingness-to-pay for MNL, MXL, and NL with analytic delta-method standard errors. For MXL, log-normal random coefficients report the *median* WTP under the package's shifted log-normal parameterization; random price coefficients are rejected
- `gof()` — goodness of fit: McFadden pseudo R-squared (plain and adjusted) and in-sample hit rate, with `"equal_shares"` (default) and `"market_shares"` null models; now also shown in the `summary()` footer
- `predict(..., newdata = )` — counterfactual prediction for all three models, from either a long data.frame in the fit-time format or a modified-design list (`X`, `alt_idx`, `M`, ...); works even with `keep_data = FALSE`. NL fits now store `nest_idx` top-level to support this
- `logsum()` — expected maximum utility (inclusive value) per choice situation for MNL (closed form), MXL (simulated with a dedicated per-draw kernel, avoiding the Jensen bias of averaging utilities first), and NL (nested inclusive-value formula)
- `consumer_surplus()` — expected consumer surplus `logsum / (-alpha)` (Train 2009, Ch. 3) with a delta-method standard error of the mean CS for MNL; supports `newdata` for policy ΔCS analysis

# choicer 0.1.0

Initial CRAN release.

## Supported models

- **Multinomial Logit** (`run_mnlogit()`) — estimation, prediction, elasticities, diversion ratios, BLP contraction
- **Mixed Logit** (`run_mxlogit()`) — normal and log-normal random coefficients, correlated random coefficients via Cholesky, Halton draws, elasticities, BLP contraction
- **Nested Logit** (`run_nestlogit()`) — estimation with nest-specific dissimilarity (lambda) parameters

## S3 class system

- Parent class `choicer_fit` with subclasses `choicer_mnl`, `choicer_mxl`, `choicer_nl`
- Standard methods: `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `predict()`
- Classed data objects: `choicer_data_mnl`, `choicer_data_mxl`, `choicer_data_nl` from `prepare_*_data()`

## Post-estimation generics

- `elasticities()` — methods for MNL and MXL
- `diversion_ratios()` — method for MNL
- `blp()` — BLP contraction for MNL and MXL

## API

- Dual workflow for all `run_*logit()` functions: convenience (pass `data` + column names) or advanced (pass pre-prepared `input_data`)
- Pluggable optimizer via `optimizer = "nloptr" | "optim" | <function>`
- `prepare_nl_data()` added for nested logit data preparation

## Computation

- C++ likelihoods, gradients, and analytical Hessians via Rcpp/RcppArmadillo
- OpenMP parallelization over individuals
- Log-sum-exp trick throughout for numerical stability
