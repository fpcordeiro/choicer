#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table NSE variables
utils::globalVariables(c(
  ".", "..covariate_cols", "..random_var_cols",
  "HAS_NA", "N_CHOICES", "N_OBS", "alt_int", "idx_in_group"
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
