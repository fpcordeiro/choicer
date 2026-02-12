
#' Runs mixed logit estimation
#'
#' Estimates a mixed logit model via simulated maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience}{Supply \code{data} and column names. Data preparation
#'     (\code{\link{prepare_mxl_data}}) and Halton draw generation
#'     (\code{\link{get_halton_normals}}) are handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mxl_data}} and
#'     \code{\link{get_halton_normals}} yourself, then pass the results via
#'     \code{input_data} and \code{eta_draws}.}
#' }
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1/0).
#' @param covariate_cols Vector of column names for fixed covariates.
#' @param random_var_cols Vector of column names for random coefficients.
#' @param input_data List output from \code{\link{prepare_mxl_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param eta_draws Array of shape K_w x S x N with standard normal draws.
#'   Required for the advanced workflow; auto-generated from \code{S} in the
#'   convenience workflow.
#' @param S Integer number of Halton draws per individual (convenience workflow
#'   only). Default 100.
#' @param rc_dist Integer vector indicating distribution of random coefficients
#'   (0 = normal, 1 = log-normal). Default: all normal.
#' @param rc_mean Logical indicating whether to estimate means for random
#'   coefficients.
#' @param rc_correlation Logical indicating whether random coefficients are
#'   correlated (convenience workflow). Ignored when \code{input_data} is used
#'   (taken from the prepared data).
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants.
#' @param theta_init Initial parameter vector. If \code{NULL}, zeros are used.
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters.
#' @param weights Optional weight vector (convenience workflow). If \code{NULL},
#'   equal weights are used.
#' @param outside_opt_label Label for the outside option (convenience workflow).
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @returns A \code{choicer_mxl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#'
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 50L
#' )
#' summary(fit)
#' }
#' @importFrom nloptr nloptr
#' @export
run_mxlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    random_var_cols = NULL,
    input_data = NULL,
    eta_draws = NULL,
    S = 100L,
    rc_dist = NULL,
    rc_mean = FALSE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    theta_init = NULL,
    optimizer = NULL,
    control = list(),
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
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
        is.null(covariate_cols) || is.null(random_var_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and random_var_cols.")
    }
    input_data <- prepare_mxl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      random_var_cols = random_var_cols,
      weights = weights,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation
    )
    K_w <- ncol(input_data$W)
    eta_draws <- get_halton_normals(S, input_data$N, K_w)
  } else {
    # Advanced workflow
    if (is.null(eta_draws)) {
      stop("'eta_draws' is required when using 'input_data' (advanced workflow).")
    }
  }

  # Parameter dimensions
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_w <- ncol(input_data$W)
  rc_correlation <- input_data$rc_correlation
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- J - 1
  n_params <- K_x + mu_size + L_size + n_asc

  if (is.null(theta_init)) theta_init <- rep(0, n_params)
  if (is.null(rc_dist)) rc_dist <- rep(0L, K_w)

  # Build eval_f closure
  eval_f <- function(theta) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      eta_draws = eta_draws,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
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

  # Build parameter names
  beta_names <- colnames(input_data$X)
  mu_names <- if (rc_mean) paste0("Mu_", colnames(input_data$W)) else character(0)

  # Cholesky parameter names (lower triangular order)
  if (rc_correlation) {
    sigma_names <- character(L_size)
    idx <- 1
    for (j in seq_len(K_w)) {
      for (i in j:K_w) {
        sigma_names[idx] <- sprintf("L_%d%d", i, j)
        idx <- idx + 1
      }
    }
  } else {
    sigma_names <- paste0("L_", seq_len(K_w), seq_len(K_w))
  }

  alt_col <- names(input_data$alt_mapping)[2]
  asc_names <- paste0("ASC_", input_data$alt_mapping[2:J][[alt_col]])

  param_names <- c(beta_names, mu_names, sigma_names, asc_names)
  names(theta_hat) <- param_names

  # Parameter index map
  pos <- 0
  param_map <- list(beta = seq_len(K_x))
  pos <- K_x
  if (mu_size > 0) {
    param_map$mu <- pos + seq_len(mu_size)
    pos <- pos + mu_size
  }
  param_map$sigma <- pos + seq_len(L_size)
  pos <- pos + L_size
  param_map$asc <- pos + seq_len(n_asc)

  # Compute vcov eagerly
  hess <- mxl_hessian_parallel(
    theta = theta_hat,
    X = input_data$X,
    W = input_data$W,
    alt_idx = input_data$alt_idx,
    choice_idx = input_data$choice_idx,
    M = input_data$M,
    weights = input_data$weights,
    eta_draws = eta_draws,
    rc_dist = rc_dist,
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option
  )
  vcov_result <- invert_hessian(hess)
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # Reconstruct Sigma for display
  L_params <- theta_hat[param_map$sigma]
  sigma_mat <- build_var_mat(L_params, K_w, rc_correlation)
  w_names <- colnames(input_data$W)
  if (!is.null(w_names)) {
    rownames(sigma_mat) <- w_names
    colnames(sigma_mat) <- w_names
  }

  # Draws info (metadata only, not the full array)
  draws_info <- list(
    S = dim(eta_draws)[2],
    N = dim(eta_draws)[3],
    K_w = K_w
  )

  # Build S3 object
  new_choicer_mxl(
    call = cl,
    coefficients = theta_hat,
    loglik = -opt$value,
    nobs = input_data$N,
    n_params = n_params,
    convergence = opt$convergence,
    message = opt$message,
    data_spec = input_data$data_spec,
    alt_mapping = input_data$alt_mapping,
    param_map = param_map,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option,
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
        X = input_data$X,
        W = input_data$W,
        alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx,
        M = input_data$M,
        weights = input_data$weights
      )
    },
    draws_info = draws_info,
    rc_dist = rc_dist,
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    sigma = sigma_mat
  )
}


#' Prepare inputs for `mxl_loglik_gradient_parallel()`
#'
#' Prepares and validates inputs for mixed logit estimation routine.
#'
#' @param data Data frame containing choice data
#' @param id_col Name of the column identifying choice situations (individuals)
#' @param alt_col Name of the column identifying alternatives
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen)
#' @param covariate_cols Vector of names of columns to be used as covariates
#' @param random_var_cols Vector of names of columns to be used as random variables
#' @param weights Optional vector of weights for each choice situation. If NULL, equal weights are used.
#' @param outside_opt_label Label for the outside option (if any). If NULL, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param rc_correlation Logical indicating whether random coefficients are correlated. Default is FALSE.
#' @return A list containing prepared inputs for mixed logit estimation, including rc_correlation.
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
#' str(input$X)
#' str(input$W)
#' @import data.table
#' @export
prepare_mxl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    random_var_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    rc_correlation = FALSE
) {

  ## Preliminary housekeeping --------------------------------------------------
  dt <- as.data.table(data)[]

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols, random_var_cols)
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

  ## Sanity checks ---------------------------------------------------------

  ## covariates must be numeric
  if (!all(vapply(dt[, ..covariate_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")
  if (!all(vapply(dt[, ..random_var_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")

  ## choice column must be 0 or 1
  bad_choice <- dt[[choice_col]] %in% c(0, 1) == FALSE
  if (any(bad_choice))
    stop("`", choice_col, "` must contain only 0 and 1.")

  ## Exactly one '1' per choice situation
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
  X <- as.matrix(dt[, ..covariate_cols])
  X_res <- check_collinearity(X)
  X <- X_res$mat
  if (!is.null(X_res$dropped)) dropped_vars <- X_res$dropped # accumulate dropped vars if we had multiple checks


  W <- as.matrix(dt[, ..random_var_cols])
  W_res <- check_collinearity(W)
  W <- W_res$mat
  if (!is.null(W_res$dropped)) {
     if(exists("dropped_vars")) dropped_vars <- c(dropped_vars, W_res$dropped)
     else dropped_vars <- W_res$dropped
  }

  cols_to_drop <- union(covariate_cols, random_var_cols)
  dt[, (cols_to_drop) := NULL]

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$alt_int)                             # length == sum(M)

  ## M[i] - # alternatives per choice situation
  M <- dt[, .N, by = id_col][["N"]]                             # length N

  ## N: number of individuals / choice situations
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

  # Weights default = 1
  if (is.null(weights)) weights <- rep(1, N)

  ## Alternative summary -------------------------------------------------------

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
    all(is.finite(X)),
    all(is.finite(W))
  )

  ## Return output -------------------------------------------------------------
  structure(
    list(
      X           = X,
      W           = W,
      alt_idx     = alt_idx,
      choice_idx  = as.integer(choice_idx),
      M           = M,
      N           = N,
      weights     = weights,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      alt_mapping = alt_mapping[],
      dropped_cols = if(exists("dropped_vars")) dropped_vars else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        random_var_cols = random_var_cols,
        outside_opt_label = outside_opt_label
      )
    ),
    class = "choicer_data_mxl"
  )
}

#' Halton draws for mixed logit
#'
#' Create halton normal draws in appropriate format for mixed logit estimation
#'
#' @param S Number of draws for each choice situation
#' @param N number of choice situations
#' @param K_w dimension of random coefficients (number of columns in W matrix)
#' @return K_w x S x N array with halton standard normal draws
#' @examples
#' draws <- get_halton_normals(S = 50, N = 10, K_w = 2)
#' dim(draws)  # 2 x 50 x 10
#' @importFrom randtoolbox halton
#' @export
get_halton_normals <- function(S, N, K_w) {
  # Generate all needed Halton draws at once
  # We need S * N draws for each of K_w dimensions
  total_draws <- S * N

  # Generate Halton sequence
  # randtoolbox::halton returns a matrix of size total_draws x K_w
  # (but drops to vector when K_w = 1, so ensure matrix)
  halton_seq <- randtoolbox::halton(n = total_draws, dim = K_w, normal = TRUE)
  if (!is.matrix(halton_seq)) halton_seq <- matrix(halton_seq, ncol = 1)

  # Initialize the eta_draws array
  eta_draws <- array(0, dim = c(K_w, S, N))
  
  # Fill the array
  # The original code used: start_index = (i - 1) * S + 1 for each individual
  # This corresponds to taking chunks of S rows from the halton sequence
  for (i in 1:N) {
     start_row <- (i - 1) * S + 1
     end_row   <- i * S
     
     # halton_seq[start:end, ] is S x K_x
     # we want K_x x S for eta_draws[, , i]
     eta_draws[, , i] <- t(halton_seq[start_row:end_row, , drop=FALSE])
  }

  return(eta_draws)
}


#' Compute aggregate elasticities for mixed logit model
#'
#' Computes the aggregate elasticity matrix (weighted average of individual
#' elasticities) for the Mixed Logit model. The elasticity E(i,j) represents
#' the percentage change in the probability of choosing alternative i when
#' the attribute of alternative j changes by 1%.
#'
#' @param opt_result Result object from `run_mxlogit()` containing `$solution`.
#' @param input_data List output from `prepare_mxl_data()`.
#' @param eta_draws Array of shape K_w x S x N with standard normal draws.
#' @param rc_dist Integer vector indicating distribution (0 = normal, 1 = log-normal).
#' @param elast_var_idx 1-based index of variable for elasticity computation.
#'   If `is_random_coef = FALSE`, this indexes into X columns. If `is_random_coef = TRUE`,
#'   this indexes into W columns.
#' @param is_random_coef Logical. `TRUE` if the variable has a random coefficient (is in W),
#'   `FALSE` if it has a fixed coefficient (is in X). Default is `FALSE`.
#' @param rc_mean Logical indicating whether means were estimated for random coefficients.
#' @param use_asc Logical indicating whether ASCs were included in the model.
#' @returns A J x J elasticity matrix where entry (i, j) is the elasticity of P_i
#'   with respect to the attribute of alternative j. Row and column names are set
#'   from the alternative labels if available.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
#' eta <- get_halton_normals(50, input$N, ncol(input$W))
#' fit <- run_mxlogit(input_data = input, eta_draws = eta)
#' elas <- mxl_elasticities(
#'   list(solution = coef(fit)), input, eta,
#'   rc_dist = rep(0L, ncol(input$W)), elast_var_idx = 1L
#' )
#' }
#' @export
mxl_elasticities <- function(
    opt_result,
    input_data,
    eta_draws,
    rc_dist,
    elast_var_idx,
    is_random_coef = FALSE,
    rc_mean = FALSE,
    use_asc = TRUE
) {
  if (is.null(opt_result$solution)) {
    stop("opt_result must contain '$solution' with the parameter estimates.")
  }

  elas_mat <- mxl_elasticities_parallel(
    theta = opt_result$solution,
    X = input_data$X,
    W = input_data$W,
    alt_idx = input_data$alt_idx,
    choice_idx = input_data$choice_idx,
    M = input_data$M,
    weights = input_data$weights,
    eta_draws = eta_draws,
    rc_dist = rc_dist,
    elast_var_idx = elast_var_idx,
    is_random_coef = is_random_coef,
    rc_correlation = input_data$rc_correlation,
    rc_mean = rc_mean,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option
  )

  # Add row and column names from alt_mapping if available
  if (!is.null(input_data$alt_mapping)) {
    alt_labels <- input_data$alt_mapping[[2]]  # Second column has alternative labels
    if (length(alt_labels) == nrow(elas_mat)) {
      rownames(elas_mat) <- alt_labels
      colnames(elas_mat) <- alt_labels
    }
  }

  return(elas_mat)
}

#' BLP contraction mapping for mixed logit
#'
#' Finds the ASC (delta) parameters such that predicted market shares
#' match target shares, using the contraction mapping of Berry, Levinsohn,
#' and Pakes (1995). This is useful for demand estimation and counterfactual
#' simulations.
#'
#' @param delta_init Initial guess for delta (ASC) values. Should be J-1 elements
#'   if `include_outside_option = FALSE`, or J elements if `include_outside_option = TRUE`.
#' @param target_shares Target market shares. Must have J elements (including
#'   outside option if applicable).
#' @param input_data List output from `prepare_mxl_data()`.
#' @param beta Fixed coefficients vector (K_x elements).
#' @param mu Mean parameters for random coefficients (K_w elements). Use zeros
#'   if `rc_mean = FALSE`.
#' @param L_params Cholesky parameters vector. Length is K_w * (K_w + 1) / 2 if
#'   `rc_correlation = TRUE`, or K_w if `rc_correlation = FALSE`.
#' @param eta_draws Array of shape K_w x S x N with standard normal draws.
#' @param rc_dist Integer vector indicating distribution (0 = normal, 1 = log-normal).
#' @param rc_mean Logical indicating whether means are estimated for random coefficients.
#' @param tol Convergence tolerance. Default is 1e-8.
#' @param max_iter Maximum number of iterations. Default is 1000.
#' @returns Converged delta (ASC) vector. Length is J-1 if `include_outside_option = FALSE`,
#'   or J if `include_outside_option = TRUE`.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
#' eta <- get_halton_normals(50, input$N, ncol(input$W))
#' fit <- run_mxlogit(input_data = input, eta_draws = eta)
#' pm <- fit$param_map
#' delta <- mxl_blp(
#'   delta_init = rep(0, J), target_shares = rep(1/J, J),
#'   input, beta = coef(fit)[pm$beta], mu = rep(0, ncol(input$W)),
#'   L_params = coef(fit)[pm$sigma], eta, rc_dist = rep(0L, ncol(input$W))
#' )
#' }
#' @export
mxl_blp <- function(
    delta_init,
    target_shares,
    input_data,
    beta,
    mu,
    L_params,
    eta_draws,
    rc_dist,
    rc_mean = FALSE,
    tol = 1e-8,
    max_iter = 1000
) {
  delta_out <- mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = input_data$X,
    W = input_data$W,
    beta = beta,
    mu = mu,
    L_params = L_params,
    alt_idx = input_data$alt_idx,
    M = input_data$M,
    weights = input_data$weights,
    eta_draws = eta_draws,
    rc_dist = rc_dist,
    rc_correlation = input_data$rc_correlation,
    rc_mean = rc_mean,
    include_outside_option = input_data$include_outside_option,
    tol = tol,
    max_iter = max_iter
  )

  return(delta_out)
}
