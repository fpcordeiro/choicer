
#' Runs multinomial logit estimation
#'
#' Estimates a multinomial logit model via maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience (default)}{Supply \code{data} and column names. Data
#'     preparation (\code{\link{prepare_mnl_data}}) is handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mnl_data}} yourself and pass the
#'     result via \code{input_data}.}
#' }
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations (individuals).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of columns to be used as covariates.
#' @param input_data List output from \code{\link{prepare_mnl_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"}, or
#'   a custom function with signature \code{f(theta_init, eval_f, lower, upper, control)}
#'   where \code{eval_f(theta)} returns \code{list(objective, gradient)}.
#'   Must return a list with \code{par}/\code{value} (or \code{solution}/\code{objective}).
#'   If the custom function accepts \code{control} or \code{...}, the \code{control}
#'   argument is forwarded; otherwise it is silently ignored.
#' @param control List of optimizer-specific control parameters passed to the
#'   chosen optimizer (e.g., \code{list(maxeval = 2000)} for nloptr).
#' @param weights Optional vector of weights for each choice situation. If \code{NULL}, equal weights are used.
#' @param outside_opt_label Label for the outside option (if any). If \code{NULL}, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param use_asc Logical indicating whether to include alternative-specific constants (ASCs) in the model.
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in the
#'   returned object for \code{predict()} and post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control} instead.
#' @returns A \code{choicer_mnl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()}, \code{vcov()},
#'   \code{logLik()}, \code{AIC()}, \code{BIC()}, \code{nobs()}, \code{predict()}.
#' @importFrom nloptr nloptr
#' @export
run_mnlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    input_data = NULL,
    optimizer = NULL,
    control = list(),
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    use_asc = TRUE,
    keep_data = TRUE,
    nloptr_opts = NULL
) {
  cl <- match.call()

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "and covariate_cols.")
    }
    input_list <- prepare_mnl_data(
      data, id_col, alt_col, choice_col, covariate_cols,
      weights = weights,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option
    )
  } else {
    # Advanced workflow: use input_data directly
    input_list <- input_data
  }

  # Resolve alt_col for parameter naming (needed by both workflows)
  if (is.null(alt_col)) alt_col <- names(input_list$alt_mapping)[2]

  # Parameter dimensions
  K_x <- ncol(input_list$X)
  J <- nrow(input_list$alt_mapping)
  n_asc <- if (use_asc) J - 1 else 0
  n_params <- K_x + n_asc
  theta_init <- rep(0, n_params)

  # Build eval_f closure (captures data in environment)
  eval_f <- function(theta) {
    mnl_loglik_gradient_parallel(
      theta = theta,
      X = input_list$X,
      alt_idx = input_list$alt_idx,
      choice_idx = input_list$choice_idx,
      M = input_list$M,
      weights = input_list$weights,
      use_asc = use_asc,
      include_outside_option = input_list$include_outside_option
    )
  }

  # Run optimizer
  elapsed <- system.time({
    opt <- run_optimizer(
      optimizer = optimizer,
      theta_init = theta_init,
      eval_f = eval_f,
      control = control
    )
  })

  message("Optimization run time ", convertTime(elapsed))

  # Parameter names and index map
  theta_hat <- opt$par
  asc_names <- if (use_asc) {
    paste0("ASC_", input_list$alt_mapping[2:J][[alt_col]])
  } else {
    character(0)
  }
  param_names <- c(colnames(input_list$X), asc_names)
  names(theta_hat) <- param_names

  param_map <- list(beta = seq_len(K_x))
  if (n_asc > 0) param_map$asc <- K_x + seq_len(n_asc)

  # Compute vcov eagerly
  hess <- mnl_loglik_hessian_parallel(
    theta = theta_hat,
    X = input_list$X,
    alt_idx = input_list$alt_idx,
    choice_idx = input_list$choice_idx,
    M = input_list$M,
    weights = input_list$weights,
    use_asc = use_asc,
    include_outside_option = input_list$include_outside_option
  )
  vcov_result <- invert_hessian(hess)
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # Build S3 object
  new_choicer_mnl(
    call = cl,
    coefficients = theta_hat,
    loglik = -opt$value,
    nobs = input_list$N,
    n_params = n_params,
    convergence = opt$convergence,
    message = opt$message,
    data_spec = input_list$data_spec %||% list(
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      outside_opt_label = outside_opt_label
    ),
    alt_mapping = input_list$alt_mapping,
    param_map = param_map,
    use_asc = use_asc,
    include_outside_option = input_list$include_outside_option,
    optimizer = list(
      name = if (is.function(optimizer)) "custom" else (optimizer %||% "nloptr"),
      control = control,
      elapsed_time = elapsed[["elapsed"]],
      iterations = opt$iterations
    ),
    vcov = vcov_result$vcov,
    se = vcov_result$se,
    data = if (keep_data) {
      list(
        X = input_list$X,
        alt_idx = input_list$alt_idx,
        choice_idx = input_list$choice_idx,
        M = input_list$M,
        weights = input_list$weights
      )
    }
  )
}

#' Prepare inputs for `mnl_loglik_gradient_parallel()`
#'
#' Prepares and validates inputs for multinomial logit estimation routine.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (individuals).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of columns to be used as covariates.
#' @param weights Optional vector of weights for each choice situation. If `NULL`, equal weights are used.
#' @param outside_opt_label Label for the outside option (if any). If `NULL`, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @returns A list containing:
#'   \itemize{
#'     \item `X`: Design matrix (sum(M) x K).
#'     \item `alt_idx`: Integer vector of alternative indices.
#'     \item `choice_idx`: Integer vector of chosen alternative indices.
#'     \item `M`: Integer vector with number of alternatives per choice situation.
#'     \item `N`: Number of choice situations.
#'     \item `weights`: Vector of weights.
#'     \item `include_outside_option`: Logical flag.
#'     \item `alt_mapping`: Data.table mapping alternatives to summary statistics.
#'     \item `dropped_cols`: Names of columns dropped due to collinearity, if any.
#'   }
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

  ## Remove outside-option rows when modelling it implicitly ------------------
  if (include_outside_option && !is.null(outside_opt_label)) {
    dt <- dt[get(alt_col) != outside_opt_label]
    if (nrow(dt) == 0) {
      stop("No inside alternatives remain after removing outside option rows.")
    }
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
  X <- as.matrix(dt[, ..covariate_cols])                        # sum(M) x K
  dt[, (covariate_cols) := NULL]
  X_res <- check_collinearity(X)
  X <- X_res$mat
  if (!is.null(X_res$dropped)) dropped_vars <- X_res$dropped

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$alt_int)                             # length == sum(M)

  ## M[i] - # alternatives per choice situation
  M <- dt[, .N, by = id_col][["N"]]                             # length N

  ## N - number of individuals / choice situations
  ids <- dt[, get(id_col)][!duplicated(dt[[id_col]])]  # vector of ids in *current* order
  N   <- length(ids)

  ## choice_idx[i] - 1-based index *within* the choice set data
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
  structure(
    list(
      X           = X,
      alt_idx     = alt_idx,
      choice_idx  = as.integer(choice_idx),
      M           = M,
      N           = N,
      weights     = weights,
      include_outside_option = include_outside_option,
      alt_mapping = alt_mapping,
      dropped_cols = if(exists("dropped_vars")) dropped_vars else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        outside_opt_label = outside_opt_label
      )
    ),
    class = "choicer_data_mnl"
  )
}

