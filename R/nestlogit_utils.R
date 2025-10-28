
#' Runs nested logit estimation
#'
#' Wrapper for nested logit model estimation using nloptr optimizer
#'
#' @param input_data List containing prepared input data for estimation
#' @param use_asc Logical indicating whether to include alternative specific constants (ASCs)
#' @param theta_init Optional initial parameter vector for optimization. If NULL, a default vector is used.
#' @param param_names Optional vector of parameter names for result summary. If NULL, default names are generated.
#' @param nloptr_opts Optional list of options for the nloptr optimizer. If NULL, default options are used.
#' @param path_output Optional file path to save the coefficient summary table as a CSV file. If NULL, no file is saved
#' @return List containing optimization result and additional information
#' @importFrom nloptr nloptr
#' @export
run_nestlogit <- function(
    input_data,
    use_asc = TRUE,
    theta_init = NULL,
    param_names = NULL,
    nloptr_opts = NULL,
    path_output = NULL
) {

  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_l <- sum(table(input_data$nest_idx) > 1)

  # Initial parameter vector theta_init
  if (is.null(theta_init) && use_asc) {
    theta_init <- c(rep(0, K_x), rep(0.5, K_l),  rep(0, J - 1))
  } else if (is.null(theta_init) && !use_asc) {
    theta_init <- c(rep(0, K_x), rep(0.5, K_l))
  }

  if (use_asc) {
    theta_lb <- c(rep(-Inf, K_x), rep(1e-16, K_l), rep(-Inf, J - 1))
  } else {
    theta_lb <- c(rep(-Inf, K_x), rep(1e-16, K_l))
  }

  if (is.null(nloptr_opts)) {
    nloptr_opts <- list(
      "algorithm" = "NLOPT_LD_LBFGS",
      "xtol_rel" = 1.0e-8,
      "maxeval" = 1e+3,
      "print_level" = 1L
    )
  }

  time <- system.time({
    # Run the optimization
    result <- nloptr::nloptr(
      x0 = theta_init,
      eval_f = nl_loglik_gradient_parallel,
      lb = theta_lb,
      opts = nloptr_opts,
      X = input_data$X,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      nest_idx = input_data$nest_idx,
      M = input_data$M,
      weights = input_data$weights,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
  })

  cat("Optimization run time", convertTime(time), "\n\n")

  get_nl_result(
    opt_result = result,
    X = input_data$X,
    alt_idx = input_data$alt_idx,
    choice_idx = input_data$choice_idx,
    nest_idx = input_data$nest_idx,
    M = input_data$M,
    weights = input_data$weights,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option,
    omit_asc_output = FALSE,
    param_names = param_names,
    file_name = path_output
  )

  result$alt_mapping <- input_data$alt_mapping

  return(result)
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

#' Coefficient summary table for nested logit model
#'
#' Prints and saves coefficient summary table for nested logit model
#'
#' @param opt_result Result object from nloptr optimization containing at least 'solution' element
#' @param X Design matrix used in estimation
#' @param alt_idx Alternative indices for each observation
#' @param choice_idx Chosen alternative indices for each choice situation
#' @param nest_idx Nest indices for each alternative (same length as number of alternatives)
#' @param M Vector of number of alternatives per choice situation
#' @param weights Vector of weights for each choice situation
#' @param use_asc Logical indicating whether ASCs were included in the model
#' @param include_outside_option Logical indicating whether an outside option was included in the model
#' @param omit_asc_output Logical indicating whether to omit ASC parameters from the output table
#' @param param_names Optional vector of parameter names. If NULL, default names are generated.
#' @param file_name Optional file path to save the coefficient summary table as a CSV file. If NULL, no file is saved.
#' @importFrom stats pnorm
#' @export
get_nl_result <- function(
    opt_result,
    X, alt_idx, choice_idx, nest_idx, M, weights,
    use_asc = TRUE,
    include_outside_option = TRUE,
    omit_asc_output = FALSE,
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
  hess <- nl_loglik_numeric_hessian(
    theta = est_theta,
    X = X,
    alt_idx = alt_idx,
    choice_idx = choice_idx,
    nest_idx = nest_idx,
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
    K_l <- length(unique(nest_idx))
    param_names <- character(p_len)

    # Beta parameters
    for (i in seq_len(K_x)) {
      param_names[i] <- paste0("X_", i)
    }

    # Lambda parameters (inclusive value)
    for (i in seq_len(K_l)) {
      param_names[K_x + i] <- paste0("Lambda_", i)
    }

    # ASC parameters (if used)
    if (use_asc && !omit_asc_output) {
      delta_length <- p_len - K_x -K_l
      if (delta_length > 0) {
        for (d in seq_len(delta_length)) {
          # If outside option is included, name them delta_1,... else skip 'delta_1'
          if (include_outside_option) {
            param_names[K_x + K_l + d] <- paste0("ASC_", d)
          } else {
            param_names[K_x + K_l + d] <- paste0("ASC_", d + 1)
          }
        }
      }
    } else {
      delta_length <- 0
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
  if (!is.null(file_name)) utils::write.csv(res_df, file_name, row.names = FALSE)

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

