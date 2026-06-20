# Package index

## Model fitting

Estimate a discrete-choice model.

- [`run_mnlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mnlogit.md)
  : Runs multinomial logit estimation
- [`run_mxlogit()`](https://fpcordeiro.github.io/choicer/reference/run_mxlogit.md)
  : Runs mixed logit estimation
- [`run_nestlogit()`](https://fpcordeiro.github.io/choicer/reference/run_nestlogit.md)
  : Runs nested logit estimation
- [`run_mnprobit()`](https://fpcordeiro.github.io/choicer/reference/run_mnprobit.md)
  : Runs Bayesian multinomial probit estimation

## Data preparation

Build design matrices and inputs for the estimators.

- [`prepare_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnl_data.md)
  :

  Prepare inputs for
  [`mnl_loglik_gradient_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_loglik_gradient_parallel.md)

- [`prepare_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mxl_data.md)
  :

  Prepare inputs for
  [`mxl_loglik_gradient_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_loglik_gradient_parallel.md)

- [`prepare_nl_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_nl_data.md)
  : Prepare inputs for nested logit estimation

- [`prepare_mnp_data()`](https://fpcordeiro.github.io/choicer/reference/prepare_mnp_data.md)
  :

  Prepare inputs for
  [`mnp_gibbs()`](https://fpcordeiro.github.io/choicer/reference/mnp_gibbs.md)

## Demand and substitution

Predicted shares, elasticities, diversion and share inversion.

- [`predict(`*`<choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/predict.choicer_mnl.md)
  : Predict from a multinomial logit model
- [`predict(`*`<choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/predict.choicer_mxl.md)
  : Predict from a mixed logit model
- [`predict(`*`<choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/predict.choicer_nl.md)
  : Predict from a nested logit model
- [`elasticities()`](https://fpcordeiro.github.io/choicer/reference/elasticities.md)
  : Compute aggregate elasticities
- [`elasticities(`*`<choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/elasticities.choicer_mnl.md)
  : Elasticities for multinomial logit model
- [`elasticities(`*`<choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/elasticities.choicer_mxl.md)
  : Elasticities for mixed logit model
- [`elasticities(`*`<choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/elasticities.choicer_nl.md)
  : Elasticities for nested logit model
- [`diversion_ratios()`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.md)
  : Compute aggregate diversion ratios
- [`diversion_ratios(`*`<choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.choicer_mnl.md)
  : Diversion ratios for multinomial logit model
- [`diversion_ratios(`*`<choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.choicer_mxl.md)
  : Diversion ratios for mixed logit model
- [`diversion_ratios(`*`<choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/diversion_ratios.choicer_nl.md)
  : Diversion ratios for nested logit model
- [`blp()`](https://fpcordeiro.github.io/choicer/reference/blp.md) : BLP
  contraction mapping
- [`blp(`*`<choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/blp.choicer_mnl.md)
  : BLP contraction mapping for multinomial logit model
- [`blp(`*`<choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/blp.choicer_mxl.md)
  : BLP contraction mapping for mixed logit model
- [`blp(`*`<choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/blp.choicer_nl.md)
  : BLP contraction mapping for nested logit model
- [`blp_contraction()`](https://fpcordeiro.github.io/choicer/reference/blp_contraction.md)
  : BLP95 contraction mapping to find delta given target shares

## Welfare

Willingness-to-pay and consumer surplus.

- [`wtp()`](https://fpcordeiro.github.io/choicer/reference/wtp.md) :
  Compute willingness to pay
- [`consumer_surplus()`](https://fpcordeiro.github.io/choicer/reference/consumer_surplus.md)
  : Expected consumer surplus
- [`logsum()`](https://fpcordeiro.github.io/choicer/reference/logsum.md)
  : Expected logsum (inclusive value) per choice situation

## Goodness of fit and methods

- [`gof()`](https://fpcordeiro.github.io/choicer/reference/gof.md) :
  Goodness of fit for a fitted choice model
- [`summary(`*`<choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/summary.choicer_mnl.md)
  : Summary for multinomial logit model
- [`summary(`*`<choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/summary.choicer_mnp.md)
  : Summary for Bayesian multinomial probit model
- [`summary(`*`<choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/summary.choicer_mxl.md)
  : Summary for mixed logit model
- [`summary(`*`<choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/summary.choicer_nl.md)
  : Summary for nested logit model
- [`coef(`*`<choicer_fit>`*`)`](https://fpcordeiro.github.io/choicer/reference/coef.choicer_fit.md)
  : Extract coefficients from a choicer_fit object
- [`coef(`*`<choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/coef.choicer_mnp.md)
  : Extract coefficients from a choicer_mnp object
- [`vcov(`*`<choicer_fit>`*`)`](https://fpcordeiro.github.io/choicer/reference/vcov.choicer_fit.md)
  : Extract variance-covariance matrix from a choicer_fit object
- [`vcov(`*`<choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/vcov.choicer_mnp.md)
  : Extract variance-covariance matrix from a choicer_mnp object
- [`logLik(`*`<choicer_fit>`*`)`](https://fpcordeiro.github.io/choicer/reference/logLik.choicer_fit.md)
  : Extract log-likelihood from a choicer_fit object
- [`nobs(`*`<choicer_fit>`*`)`](https://fpcordeiro.github.io/choicer/reference/nobs.choicer_fit.md)
  : Extract number of observations from a choicer_fit object
- [`nobs(`*`<choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/nobs.choicer_mnp.md)
  : Extract number of observations from a choicer_mnp object
- [`print(`*`<choicer_cs>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.choicer_cs.md)
  : Print a consumer surplus summary
- [`print(`*`<choicer_fit>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.choicer_fit.md)
  : Print a choicer_fit object
- [`print(`*`<choicer_gof>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.choicer_gof.md)
  : Print goodness-of-fit measures
- [`print(`*`<choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.choicer_mnp.md)
  : Print a choicer_mnp object
- [`print(`*`<choicer_wtp>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.choicer_wtp.md)
  : Print a WTP table
- [`print(`*`<summary.choicer_mnl>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.summary.choicer_mnl.md)
  : Print summary for multinomial logit model
- [`print(`*`<summary.choicer_mnp>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.summary.choicer_mnp.md)
  : Print summary for Bayesian multinomial probit model
- [`print(`*`<summary.choicer_mxl>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.summary.choicer_mxl.md)
  : Print summary for mixed logit model
- [`print(`*`<summary.choicer_nl>`*`)`](https://fpcordeiro.github.io/choicer/reference/print.summary.choicer_nl.md)
  : Print summary for nested logit model

## Simulation and recovery

Data-generating processes and parameter-recovery diagnostics.

- [`simulate_mnl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnl_data.md)
  : Simulate multinomial logit data

- [`simulate_mxl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mxl_data.md)
  : Simulate mixed logit data

- [`simulate_nl_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_nl_data.md)
  : Simulate nested logit data

- [`simulate_mnp_data()`](https://fpcordeiro.github.io/choicer/reference/simulate_mnp_data.md)
  : Simulate multinomial probit data

- [`recovery_table()`](https://fpcordeiro.github.io/choicer/reference/recovery_table.md)
  : Parameter recovery table

- [`monte_carlo()`](https://fpcordeiro.github.io/choicer/reference/monte_carlo.md)
  : Monte Carlo parameter recovery

- [`mc_asymptotics()`](https://fpcordeiro.github.io/choicer/reference/mc_asymptotics.md)
  : Asymptotic diagnostics for a Monte Carlo study

- [`new_choicer_sim()`](https://fpcordeiro.github.io/choicer/reference/new_choicer_sim.md)
  :

  Construct a `choicer_sim` object

## Choice-based sampling

- [`wesml_weights()`](https://fpcordeiro.github.io/choicer/reference/wesml_weights.md)
  : WESML weights for choice-based (endogenous stratified) samples
- [`sample_by_choice()`](https://fpcordeiro.github.io/choicer/reference/sample_by_choice.md)
  : Draw a choice-based sample stratified by the chosen alternative
- [`wesml_vcov()`](https://fpcordeiro.github.io/choicer/reference/wesml_vcov.md)
  : Robust (sandwich) variance for a weighted / choice-based mixed logit
  fit

## Configuration and helpers

- [`set_num_threads()`](https://fpcordeiro.github.io/choicer/reference/set_num_threads.md)
  : Set the number of OpenMP threads used by choicer
- [`thread_info()`](https://fpcordeiro.github.io/choicer/reference/thread_info.md)
  : Query choicer OpenMP thread settings
- [`get_halton_normals()`](https://fpcordeiro.github.io/choicer/reference/get_halton_normals.md)
  : Halton draws for mixed logit
- [`build_var_mat()`](https://fpcordeiro.github.io/choicer/reference/build_var_mat.md)
  : Reconstruct variance matrix L from L_params
- [`jacobian_vech_Sigma()`](https://fpcordeiro.github.io/choicer/reference/jacobian_vech_Sigma.md)
  : Utility to compute analytical Jacobian of random coefficient matrix
  transformed by vech (dVech(Sigma) / dTheta)

## Data

- [`mode_choice`](https://fpcordeiro.github.io/choicer/reference/mode_choice.md)
  : Intercity travel mode choice

## Low-level C++ kernels

Direct access to the C++ likelihood, gradient, Hessian and
post-estimation kernels. Most users do not need these.

- [`mnl_diversion_ratios_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_diversion_ratios_parallel.md)
  : Compute MNL diversion ratios (parallelized over individuals)
- [`mnl_elasticities_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_elasticities_parallel.md)
  : Compute aggregate elasticities for MNL model
- [`mnl_loglik_gradient_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_loglik_gradient_parallel.md)
  : Log-likelihood and gradient for multinomial logit model
- [`mnl_loglik_hessian_parallel()`](https://fpcordeiro.github.io/choicer/reference/mnl_loglik_hessian_parallel.md)
  : Hessian matrix for multinomial logit model
- [`mnl_predict()`](https://fpcordeiro.github.io/choicer/reference/mnl_predict.md)
  : Prediction of choice probabilities and utilities based on fitted
  model
- [`mnl_predict_shares()`](https://fpcordeiro.github.io/choicer/reference/mnl_predict_shares.md)
  : Prediction of market shares based on fitted model
- [`mxl_bhhh_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_bhhh_parallel.md)
  : BHHH (outer product of gradients) information matrix for Mixed Logit
- [`mxl_blp_contraction()`](https://fpcordeiro.github.io/choicer/reference/mxl_blp_contraction.md)
  : BLP contraction mapping for mixed logit
- [`mxl_diversion_ratios_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_diversion_ratios_parallel.md)
  : Diversion ratios for Mixed Logit (simulated, derivative-based)
- [`mxl_elasticities_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_elasticities_parallel.md)
  : Compute aggregate elasticities for mixed logit model
- [`mxl_hessian_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_hessian_parallel.md)
  : Analytical Hessian of the log-likelihood v2
- [`mxl_loglik_gradient_parallel()`](https://fpcordeiro.github.io/choicer/reference/mxl_loglik_gradient_parallel.md)
  : Log-likelihood and gradient for Mixed Logit
- [`mxl_logsum()`](https://fpcordeiro.github.io/choicer/reference/mxl_logsum.md)
  : Simulated expected logsum (inclusive value) for Mixed Logit
- [`mxl_predict()`](https://fpcordeiro.github.io/choicer/reference/mxl_predict.md)
  : Per-observation simulated choice probabilities for Mixed Logit
- [`mxl_predict_shares()`](https://fpcordeiro.github.io/choicer/reference/mxl_predict_shares.md)
  : Predicted aggregate market shares for Mixed Logit
- [`nl_blp_contraction()`](https://fpcordeiro.github.io/choicer/reference/nl_blp_contraction.md)
  : BLP95 contraction mapping for the Nested Logit model
- [`nl_diversion_ratios_parallel()`](https://fpcordeiro.github.io/choicer/reference/nl_diversion_ratios_parallel.md)
  : Compute Nested Logit diversion ratios (parallelized over
  individuals)
- [`nl_elasticities_parallel()`](https://fpcordeiro.github.io/choicer/reference/nl_elasticities_parallel.md)
  : Compute aggregate elasticities for the Nested Logit model
- [`nl_loglik_gradient_parallel()`](https://fpcordeiro.github.io/choicer/reference/nl_loglik_gradient_parallel.md)
  : Log-likelihood and gradient for Nested Logit model
- [`nl_loglik_hessian_parallel()`](https://fpcordeiro.github.io/choicer/reference/nl_loglik_hessian_parallel.md)
  : Analytical Hessian of the negated log-likelihood for the Nested
  Logit model
- [`nl_loglik_numeric_hessian()`](https://fpcordeiro.github.io/choicer/reference/nl_loglik_numeric_hessian.md)
  : Numerical Hessian of the log-likelihood via finite differences
- [`nl_predict()`](https://fpcordeiro.github.io/choicer/reference/nl_predict.md)
  : Prediction of choice probabilities and utilities for the Nested
  Logit model
- [`nl_predict_shares()`](https://fpcordeiro.github.io/choicer/reference/nl_predict_shares.md)
  : Prediction of market shares for the Nested Logit model
- [`mnp_gibbs()`](https://fpcordeiro.github.io/choicer/reference/mnp_gibbs.md)
  : Gibbs sampler for the Bayesian multinomial probit model
