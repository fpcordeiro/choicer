#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table NSE variables
utils::globalVariables(c(
  ".", "..covariate_cols", "..random_var_cols",
  "HAS_NA", "N_CHOICES", "N_OBS", "alt_int", "idx_in_group",
  # DGP / simulation NSE variables (R/simulation.R)
  "alt", "choice", "delta_val", "epsilon", "utility", "id",
  "w1", "w2", "x1", "x2", "X", "W", "V", "V_over_lambda",
  "IV", "i.IV", "nest", "nest_prob", "i.nest_prob", "log_denom",
  "cond_prob", "P_j", "j", "lambda",
  # Recovery / Monte Carlo NSE variables (R/recovery.R)
  "parameter", "converged", "estimate", "se", "true", "covers",
  "group", "rep_id", "R_success", "mean_est", "median_est",
  "sd_est", "mean_se", "bias", "rmse", "coverage"
))

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table data.table
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef logLik nobs predict vcov pnorm
#' @useDynLib choicer, .registration = TRUE
## usethis namespace: end
NULL
