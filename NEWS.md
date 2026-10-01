# choicer (development version)

## Mixed logit — panel likelihood (`person_col`) and conditional tastes

- `run_mxlogit()` / `prepare_mxl_data()` gain `person_col`: the frequentist
  panel mixed logit of Revelt and Train (1998). When supplied, all choice
  situations of a decision maker share one draw of the random coefficients,
  and the likelihood integrates the *joint* probability of that person's
  choices over the taste distribution, instead of integrating every choice
  situation separately as before. `id_col` must still identify choice
  situations uniquely across the whole data set (ids do not restart within
  person). `NULL` (the default) keeps every choice situation its own decision
  maker, reproducing the existing cross-sectional likelihood.
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
  - Apart from the corrections below, default (cross-sectional) behavior is
    unchanged up to floating-point reassociation (<= 1e-10 relative in the
    estimation kernels, measured on the reference battery of test
    configurations); prediction kernels are bit-identical.
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

## Mixed logit — kernels for population-scale data

- The estimation kernels keep each thread's working arrays across decision
  makers instead of allocating them per unit, and the analytical Hessian no
  longer allocates a temporary for every alternative-and-draw outer product
  or per-unit block: the draws are read in place from the store-mode cube,
  and a unit's design rows and random-coefficient draws go into reused
  buffers. On a synthetic 10^6-row panel the Hessian's heap traffic fell
  from 106 million allocations (374 GB) per call to about three thousand, a
  gradient evaluation on 10^7 rows allocates a quarter to a third as often,
  and the kernels ran 3-24% faster in our benchmarks.
- Each decision maker's draws are now processed in batches sized to a small
  per-thread memory budget, with the score accumulated across batches by a
  streaming log-sum-exp, so a thread's working memory no longer grows with
  the product of a decision maker's number of choice situations and the
  number of draws `S`. The draw loop also sweeps a decision maker's rows in
  memory order and no longer allocates per decision maker. On a synthetic
  claims panel with heavy-tailed histories (up to 500 choice situations per
  decision maker), peak working memory at `S = 1000` fell from 6.4 GB to
  0.12 GB and gradient evaluations ran 40% faster; on a 10^7-row panel they
  ran 12-25% faster at `S = 100`-`500`. The analytical Hessian no longer
  holds a decision maker's rows-by-draws matrix either: its peak working
  memory on the claims panel at `S = 500` fell from 1.7 GB to 0.2 GB. A
  decision maker is split into batches only when its stacked alternative
  rows times `S` exceed 2^18 (a long panel, a choice set of hundreds of
  alternatives, or a very large `S`).
- The estimation kernels read the stacked design's integer indices in place
  instead of copying them on every call, and compute each decision maker's
  base utilities (`X beta + W mu + delta`) inside the parallel loop instead
  of for all rows up front, so a likelihood evaluation no longer allocates
  anything the length of the stacked design (what remains is 8 bytes of row
  offsets per choice situation, and in a panel 8 per decision maker). On
  synthetic 10^8-row panels, peak working memory fell from about 3 GB to
  under 0.2 GB and an evaluation ran 19% faster at `S = 25` (7-9% at
  `S = 100`); fits at `S = 50` on 10^7 rows took 9% less time per
  evaluation. The kernels' offsets into the stacked design are now 64-bit.
- Generate-mode draws (`draws = "generate"`) are now formed in two passes
  over each block: all the Halton uniforms first, then the inverse normal
  CDF. The first (base-2) dimension is stepped by an exact integer odometer,
  and the digits of the next ten (bases 3 to 31) are taken by division by a
  compile-time constant. The draws are bit-identical to before on 64-bit
  platforms, so this change by itself alters no estimate, standard error,
  prediction or other post-estimation quantity; it lowers the cost of
  regenerating the draws at each likelihood evaluation. With three random
  coefficients and `S = 100`, gradient evaluations ran 32-36% faster on a
  synthetic four-alternative cross-section (10^6 rows at 1 and 11 threads,
  10^7 rows at 11), where regenerating the draws had taken about half the
  time; 10-14% faster on a school-census-like design (thousands of schools,
  10-100 per choice set); and showed no measurable change (under 3%, within
  run-to-run noise) on claims-like panels with 20-200 hospitals per choice
  set, whose cost lies elsewhere. On the four-alternative data a
  generate-mode evaluation now takes about 1.35 times as long as a
  store-mode one, down from about 2 times. The generator's
  digit-permutation tables are also smaller: 1.3 MB instead of 11 MB at the
  128-dimension maximum.
- Together these kernel changes alter results only by floating-point
  rounding: at most 3.4e-15 relative on our reference battery of kernel
  configurations, and within the 1e-10 our tests allow for decision makers
  split into draw batches. The numerical changes of substance are the
  corrections below.

## Data preparation at population scale

- `get_halton_normals()` builds the draw cube of `draws = "store"`, at fit
  time and again for post-estimation. It now fills the cube a block of units
  at a time; before, it generated the whole Halton sequence in one call and
  copied it into the cube unit by unit in an R loop. The draws are
  bit-identical wherever the old code was well defined. For 2 million draw
  units at `S = 100` with three random coefficients (a 4.8 GB cube),
  building the cube took 21 s instead of 27 s, peak memory fell from 3.0 to
  1.4 times the size of the cube (the rest is R's garbage-collection slack),
  and 1.3 thousand allocations replaced 8 million. The single call also
  computed positions in its output with 32-bit integers, which overflow past
  2^31 - 1 values (a 17 GB cube); each block now stays far below that. The
  function now checks that `S`, `N` and `K_w` are positive whole numbers,
  and stops before allocating a cube that store mode cannot handle: more
  than 2^31 - 1 points (`S * N`), the largest starting index
  `randtoolbox::halton()` accepts, where the old code failed with an
  unrelated error; or more than 2^32 - 1 values (`K_w * S * N`), the most
  the kernels can address with their 32-bit indices.
  `draws = "generate"` has neither limit. choicer now requires randtoolbox
  1.31.0 or later, the release that added `halton(start = )`.
- `prepare_mnl_data()`, `prepare_mxl_data()` and `prepare_nl_data()` (and so
  `run_mnlogit()`, `run_mxlogit()` and `run_nestlogit()`) copy only the
  columns the model uses from a data.frame, tibble or data.table, where they
  copied the whole data set; they scan those columns for missing values one
  at a time instead of through a rows-by-columns `is.na()` matrix, and check
  that the covariates are numeric without copying them. `prepare_nl_data()`
  no longer makes a second full copy to read the alternative-to-nest map. On
  a synthetic claims-style panel of 8.7 million rows, peak heap use fell from
  5.7 to 4.3 times the size of the input for `prepare_nl_data()` and from
  4.5 to 4.0 times for `prepare_mxl_data()`; the saving grows with the
  number of columns the model does not use. Apart from the storage change
  below, the prepared objects are unchanged on our reference battery and on
  every call the test suite makes.
- The same three functions now copy only the index columns (ids,
  alternatives, choices, and any weight, cluster or decision-maker column)
  into their working table, and gather the design matrices `X` and `W` in
  C++ straight from the caller's covariate columns once the rows are
  filtered and sorted. `X` and `W` are now always double precision: a
  design matrix whose columns were all integer used to be an integer
  matrix, which the estimation kernels converted to double on every call.
  Estimates, standard errors and predictions are unchanged, bit for bit.
  The functions now stop with an error when `X` or `W` would hold more than
  2^32 - 1 values (for instance 10^9 rows and 5 covariates), the most the
  estimation kernels can address with their 32-bit indices.
- The collinearity check that drops dependent covariates (the columns
  `qr(X, tol = 1e-7)` moves past its rank) now factors a design of a
  million or more values by row chunks and applies `qr()`'s own rank rule
  to the resulting factor, which has the same cross-products. It drops the
  same columns on our reference battery, on every call the test suite makes
  and on a sweep of near-collinear designs; the two can disagree only on a
  column whose dependence lies within rounding of the tolerance, where
  `qr()`'s own answer already changes with the order of the rows. The check
  no longer makes `qr()`'s copies of the design (about three times its
  size), nor stops at its limit of 2^31 - 1 values, which ended the
  preparation of any larger design (for instance 2.2 * 10^8 rows and 10
  covariates) with "too large a matrix for LINPACK". Designs up to the
  kernels' 2^32 - 1 values can now be prepared. `prepare_mnp_data()`,
  `prepare_hmnl_data()` and `prepare_hmnp_data()`, for which `qr()`'s limit
  was the only stop, now also stop with an error before building an `X` of
  more than 2^32 - 1 values.
- The checks that a weight, cluster or decision-maker column is constant
  within each choice situation now count distinct (situation, value) pairs
  in one grouping, instead of sorting every situation's values separately,
  which alone allocated 87 of the 97 GB that `prepare_mxl_data()` moved
  through the heap on the 8.7-million-row claims-style panel; the
  per-situation value is read from each situation's first row. With all the
  changes above, on that panel and on a census-style panel of 9.8 million
  rows, `prepare_mxl_data()` now takes 2.9-3.7 s instead of 5.7-7.6 s,
  allocates 7.5-8.8 GB instead of 97-149 GB, and its peak heap use is 1.9
  times the size of its input instead of 4.0-4.5 times; for
  `prepare_mnl_data()` and `prepare_nl_data()` it is 2.0-2.1 times instead
  of 4.4-5.7 times. Peak resident memory, which also counts memory R has freed
  but not returned, is 2.3-2.4 times the input for all three instead of
  5.3-8.1 times (medians of three runs).
- `prepare_mnp_data()`, `prepare_hmnl_data()` and `prepare_hmnp_data()` (and
  so `run_mnprobit()`, `run_hmnlogit()` and `run_hmnprobit()`) now also copy
  only the columns the model uses, scan them for missing values one at a
  time, and check that the covariates are numeric without copying them. The
  hierarchical preparations also scan for non-finite covariate values a
  column at a time, and look up the choice situations to drop only when a
  row is flagged, where they grouped every row by choice situation twice.
  Apart from the changes below, the prepared objects are unchanged on every
  call the test suite makes and on a battery of edge cases; a used column
  that shared its name with a working column of the old code (`HAS_NA`,
  `TASK_HAS_NA`, `HAS_BAD`, `TASK_HAS_BAD`) stopped the preparation with an
  error, and now works like any other.
- The same three functions now copy only the index columns into their
  working table and gather `X` in C++ from the caller's covariate columns,
  as the MNL, MXL and NL preparations do; `prepare_mnp_data()` writes each
  choice situation's differences from its base alternative in the same
  pass. `X` is now always double precision there too (posterior draws are
  unchanged, bit for bit). The hierarchical preparations also read the
  alternative-level covariates by source row, call a covariate constant
  within choice situations when every row equals its situation's first row,
  and look for an alternative listed twice in one pass over the keys; both
  checks used to evaluate R code for every choice situation. A few inputs
  now prepare differently. The choice-situation id used as a covariate (a
  task-order term, say) is now recognized as constant within every
  situation: it is kept with a message, or dropped with a warning when there
  is no outside option, where it used to be kept silently. Inputs whose
  column names clashed with the preparations' own names are covered under
  Corrections. Integer covariates whose differences or within-situation
  ranges overflow 32-bit integers no longer fail. On synthetic inputs of
  10 million rows in which the model uses every column, `prepare_mnp_data()`
  (2.5 million choice situations of four alternatives) now takes 5.0 s instead
  of 13.8 s, allocates 6.4 GB instead of 9.0 GB, and its peak heap use is
  3.0 times the size of its input instead of 4.9 times; `prepare_hmnl_data()`
  (500,000 respondents with five choice situations of four alternatives
  each) takes 14.5 s instead of 31.7 s, allocates 12.2 GB instead of
  15.5 GB, and peaks at 3.1 times its input instead of 4.9 times. Peak
  resident memory is 4.5 times the input instead of 8.0 times for
  `prepare_mnp_data()` and 4.7 times instead of 7.6 times for
  `prepare_hmnl_data()` (medians of three runs).
- All six preparations now count each choice situation's alternatives and
  chosen rows on private tables of the columns involved, which data.table
  sums in one pass, and number the rows within each situation in one pass
  too (`data.table::rowidv()`; `prepare_mnp_data()` has none to number);
  they used to evaluate R code for every choice situation. On the inputs of
  the previous item this halves the time `prepare_mnp_data()` takes and
  cuts that of `prepare_hmnl_data()` by two fifths; on the claims- and
  census-style panels above, `prepare_mnl_data()`, `prepare_mxl_data()` and
  `prepare_nl_data()` take 9-15% less time. Memory traffic grows by up to 3%
  (data.table's working arrays for the sums), and peak heap use is unchanged
  or lower (medians of three alternated runs). The prepared objects are
  unchanged on every call the test suite makes and on a battery of edge
  cases.

## Corrections

The following were found while implementing the panel likelihood above and
the kernel and data-preparation work for population-scale data that
followed, and independently verified. Each affected released versions
unless it says otherwise.

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
  with an error that says why. Two
  inputs are now errors: a column whose name starts with `.choicer_`, used
  as an id, alternative, choice, covariate, weight, cluster or
  decision-maker column; and an alternative column named `alt_int`,
  `N_OBS`, `N_CHOICES`, `TAKE_RATE` or `MKT_SHARE`, the fixed columns of
  the returned `alt_mapping`, which must now be renamed. Released versions
  handled `.choicer_` columns correctly (they used no
  such names) and, of the alternative-column names, only an `alt_int`
  column already coded 1 to J, without the outside option.
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
