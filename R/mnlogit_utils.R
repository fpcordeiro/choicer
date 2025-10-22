
#' Runs multinomial logit estimation
#'
#' Runs data preparation and validation, likelihood maximization, and result summary
#'
#' @export
run_mnlogit <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    nloptr_opts = NULL,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    use_asc = TRUE,
    path_output = NULL
) {

  # Validate inputs and prepare data
  input_list <- prepare_mnl_data(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    weights = weights,
    outside_opt_label,
    include_outside_option
  )

  # Initial parameter vector theta_init
  J <- nrow(input_list$alt_mapping)
  theta_init <- rep(0, ncol(input_list$X) + J - 1)

  if (is.null(nloptr_opts)) {
    nloptr_opts <- list(
      "algorithm" = "NLOPT_LD_LBFGS",
      "xtol_rel" = 1.0e-8,
      "maxeval" = 1e+3,
      "print_level" = 0L
    )
  }

  time <- system.time({
    # Run the optimization
    result <- nloptr::nloptr(
      x0 = theta_init,
      eval_f = mnl_loglik_gradient_parallel,
      opts = nloptr_opts,
      X = input_list$X,
      alt_idx = input_list$alt_idx,
      choice_idx = input_list$choice_idx,
      M = input_list$M,
      weights = input_list$weights,
      use_asc = use_asc,
      include_outside_option = input_list$include_outside_option
    )
  })

  cat("Optimization run time", convertTime(time), "\n\n")

  asc_names <- paste0("ASC_", input_list$alt_mapping[2:J][[alt_col]])
  param_names <- c(covariate_cols, asc_names)

  get_mnl_result(
    opt_result = result,
    X = input_list$X,
    alt_idx = input_list$alt_idx,
    choice_idx = input_list$choice_idx,
    M = input_list$M,
    weights = input_list$weights,
    use_asc = use_asc,
    include_outside_option = input_list$include_outside_option,
    param_names = param_names,
    eps = 1e-8,
    file_name = path_output
  )

  result$alt_mapping <- input_list$alt_mapping

  return(result)
}

#' Prepare inputs for `mnl_loglik_gradient_parallel()`
#'
#' Prepares and validates inputs for mixed logit estimation routine.
#'
#' @export
prepare_mnl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE
) {
  ## Preliminary housekeeping --------------------------------------------------
  dt <- as.data.table(data)[]

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols)
  if (!all(needed %in% names(dt)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(dt)), collapse = ", "))

  # Drop non-relevant variables
  vars_to_drop <- setdiff(names(dt), needed)
  if (length(vars_to_drop) > 0) {
    dt[, (vars_to_drop) := NULL]
  }

  ## Drop ids with missing observations ----------------------------------------
  dt[, HAS_NA := rowSums(is.na(.SD)) > 0]
  ids_to_drop <- dt[HAS_NA==TRUE, get(id_col)] |> unique()
  if (length(ids_to_drop) > 0) {
    dt <- dt[!(get(id_col) %in% ids_to_drop)]
    warning("Removed ", length(ids_to_drop),
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }
  dt[, HAS_NA := NULL]

  ## Sanity checks -------------------------------------------------------------

  ## Covariates must be numeric
  if (!all(vapply(dt[, ..covariate_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")

  ## choice column must be 0/1 and exactly one '1' per choice situation
  bad_choice <- dt[[choice_col]] %in% c(0, 1) == FALSE
  if (any(bad_choice))
    stop("`", choice_col, "` must contain only 0 and 1.")

  by_id <- dt[, .(chosen = sum(get(choice_col))), by = id_col]
  if (include_outside_option == FALSE && any(by_id$chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(by_id$chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). An id with no explicit choice is assumed to be outside option.")
  }

  ## Create integer alternative codes ------------------------------------------

  if (!is.null(outside_opt_label) && include_outside_option==FALSE) {
    levels <- c(outside_opt_label, sort(setdiff(unique(dt[[alt_col]]), outside_opt_label)))
  } else {
    levels <- sort(unique(dt[[alt_col]]))
  }

  dt[, alt_int := as.integer(factor(get(alt_col), levels = levels))]

  ## Order rows ----------------------------------------------------------------
  ##   within each id: ascending alternative id
  ##   between ids   : ascending id
  setorderv(dt, c(id_col, "alt_int"))

  ## index of each row within its choice set
  dt[, idx_in_group := seq_len(.N), by = id_col]

  ## Build objects -------------------------------------------------------------
  ## design matrix
  X <- as.matrix(dt[, ..covariate_cols])                        # sum(M) × K
  dt[, (covariate_cols) := NULL]
  X <- check_collinearity(X)

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$alt_int)                             # length == sum(M)

  ## M[i] – # alternatives per choice situation
  M <- dt[, .N, by = id_col][["N"]]                             # length N

  ## N – number of individuals / choice situations
  ids <- dt[, get(id_col)][!duplicated(dt[[id_col]])]  # vector of ids in *current* order
  N   <- length(ids)

  ## choice_idx[i] – 1-based index *within* the choice set data
  ## 0 == outside option (only if chosen = 0 for all inside options & include_outside_option == TRUE)
  if (include_outside_option) {
    # start with all-zero (everyone assumed to pick the outside good)
    choice_idx <- integer(N)
    chosen_dt <- dt[get(choice_col) == 1, .(pos = idx_in_group), by = id_col]

    # match chosen ids back to the master index vector
    setkeyv(chosen_dt, id_col)
    choice_idx[match(chosen_dt[[id_col]], ids)] <- chosen_dt$pos
  } else {
    # exactly one explicit choice per id
    choice_idx <- dt[get(choice_col) == 1, idx_in_group]
  }

  if (is.null(weights)) weights <- rep(1, length(M))

  ## Alternatives summary ------------------------------------------------------

  if (include_outside_option) {
    inside_alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
    outside_alt_mapping <- data.table(alt_int=0L, N_OBS = N, N_CHOICES = sum(choice_idx == 0L))
    outside_alt_mapping[[alt_col]] <- outside_opt_label
    alt_mapping <- list(outside_alt_mapping, inside_alt_mapping) |>
      rbindlist(use.names = TRUE, fill = TRUE)
    setcolorder(alt_mapping, c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  } else {
    alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
  }

  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    length(alt_idx)    == nrow(X),
    length(choice_idx) == N,
    length(M)          == N,
    length(weights)    == N,
    all(is.finite(X))
  )

  ## return output -------------------------------------------------------------
  list(
    X           = X,
    alt_idx     = alt_idx,
    choice_idx  = as.integer(choice_idx),
    M           = M,
    N           = N,
    weights     = weights,
    include_outside_option = include_outside_option,
    alt_mapping = alt_mapping
  )
}

# Converts elapsed time to a nicely formatted string
convertTime <- function(time) {
  et <- time["elapsed"]
  if (et < 1) {
    s <- round(et, 2)
  } else {
    s <- round(et, 0)
  }
  h <- s %/% 3600
  s <- s - 3600 * h
  m <- s %/% 60
  s <- s - 60 * m
  return(paste(h, "h:", m, "m:", s, "s", sep = ""))
}

#' Coefficient summary table for multinomial logit model
#'
#' Prints and saves coefficient summary table for multinomial logit model
#'
#' @export
get_mnl_result <- function(
    opt_result,
    X, alt_idx, choice_idx, M, weights,
    use_asc = TRUE,
    include_outside_option = TRUE,
    omit_asc_output = FALSE,
    eps = 1e-6,
    param_names = NULL,
    file_name = NULL
) {
  # Extract the estimated parameter vector
  if (is.null(opt_result$solution)) {
    stop("opt_result must contain '$solution' with the parameter estimates.")
  }
  est_theta <- opt_result$solution
  p_len <- length(est_theta)

  # Compute Hessian at est_theta
  hess <- mnl_loglik_hessian_parallel(
    theta = est_theta,
    X = X,
    alt_idx = alt_idx,
    choice_idx = choice_idx,
    M = M,
    weights = weights,
    use_asc = use_asc,
    include_outside_option = include_outside_option
  )

  # Try to invert Hessian for standard errors
  vcov_mat <- NULL
  se <- rep(NA, p_len)  # default in case of errors

  singular_flag <- FALSE
  tryCatch({
    vcov_mat <- solve(hess)
  }, error = function(e) {
    singular_flag <<- TRUE
    message("Error computing vcov_mat (likely singular Hessian): ", e$message)
  })

  if (!singular_flag && !is.null(vcov_mat)) {
    # successfully computed Hessian inverse
    se <- sqrt(diag(vcov_mat))
  } else {
    # either singular Hessian or some inversion error
    message("Standard errors set to NA due to Hessian inversion failure.")
  }

  # Compute z-values and p-values
  zval <- est_theta / se
  # two-sided p-value under normal approximation
  pval <- 2 * (1 - pnorm(abs(zval)))

  # significance codes
  significance_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("")
  }

  # Determine parameter names if not provided
  if (is.null(param_names)) {
    # We assume the first K_x = ncol(X) are X_1,...,X_K_w
    K_x <- ncol(X)
    param_names <- character(p_len)

    # Beta parameters
    for (i in seq_len(K_x)) {
      param_names[i] <- paste0("X_", i)
    }

    # ASC parameters (if used)
    if (use_asc && !omit_asc_output) {
      delta_length <- p_len - K_x
      if (delta_length > 0) {
        for (d in seq_len(delta_length)) {
          # If outside option is included, name them delta_1,... else skip 'delta_1'
          if (include_outside_option) {
            param_names[K_x + d] <- paste0("ASC_", d)
          } else {
            param_names[K_x + d] <- paste0("ASC_", d + 1)
          }
        }
      }
    }
  } else {
    # If param_names is provided, ensure it has the correct length
    if (length(param_names) != p_len) {
      stop("param_names must have the same length as the number of parameters.")
    }
  }

  # Build a data frame for results
  # Add an "Index" column that goes 1, 2, ..., p_len
  final_len <- if (omit_asc_output) ncol(X) else p_len
  res_df <- data.frame(
    Index     = seq_len(final_len),
    Parameter = param_names[1:final_len],
    Estimate  = est_theta[1:final_len],
    Std_Error = se[1:final_len],
    z_value   = zval[1:final_len],
    Pr_z      = pval[1:final_len],
    Signif    = sapply(pval[1:final_len], significance_code),
    stringsAsFactors = FALSE
  )

  # Write CSV
  if (!is.null(file_name)) write.csv(res_df, file_name, row.names = FALSE)

  # rint a formatted table to screen

  # compute dynamic column widths:
  # - index_colwidth: enough for the largest index
  index_colwidth <- max(nchar(as.character(p_len)), nchar("Index"))

  # - param_colwidth: depends on the longest parameter name
  param_colwidth <- max(nchar(res_df$Parameter), nchar("Parameter"))

  # Helper to print each row
  print_row <- function(x) {
    cat(sprintf(
      # Format: index, param, estimate, std.error, z-value, Pr(>|z|), Signif
      paste0(
        "%", index_colwidth, "d  ",    # Index
        "%-", param_colwidth, "s  ",  # Parameter name (left-justified)
        "%10.6f %10.6f %8.4f %9.2e  %s\n"
      ),
      as.integer(x["Index"]),
      x["Parameter"],
      as.numeric(x["Estimate"]),
      as.numeric(x["Std_Error"]),
      as.numeric(x["z_value"]),
      as.numeric(x["Pr_z"]),
      x["Signif"]
    ))
  }

  # Header row
  cat("Model Coefficients:\n")
  cat(sprintf(
    paste0(
      "%", index_colwidth, "s  ",
      "%-", param_colwidth, "s  ",
      "%10s %10s %8s %9s  %s\n"
    ),
    "Index", "Parameter", "Estimate", "Std.Error", "z-value", "Pr(>|z|)", ""
  ))

  # Print each row
  for (i in seq_len(nrow(res_df))) {
    print_row(res_df[i, ])
  }

  invisible(res_df)
}

remove_nullspace_cols <- function(mat, tol = 1e-7) {
  if (is.null(mat)) return(mat)
  if (ncol(mat)==1) return(mat)
  qrdecomp <- qr(mat, tol = tol)
  rank <- qrdecomp$rank
  if (rank == ncol(mat)) return(mat)
  bad_cols_idx  <- qrdecomp$pivot[(rank + 1):ncol(mat)]
  mat <- mat[, setdiff(1:ncol(mat), bad_cols_idx), drop=FALSE]
  return(mat)
}

check_collinearity <- function(X) {
  colnames_before <- colnames(X)
  X <- remove_nullspace_cols(X)
  colnames_after <- colnames(X)
  colnames_diff <- setdiff(colnames_before, colnames_after)
  if (length(colnames_diff) > 0) {
    cat("The following variables were dropped due to collinearity:\n")
    cat(colnames_diff, "\n")
    vars_drop <- colnames_diff
  }
  return(X)
}

