
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
#' @param scale_vars Pre-estimation column scaling for the design matrix. One of
#'   \code{"none"} (default), \code{"sd"} (sample standard deviation),
#'   \code{"mad"} (\code{stats::mad}), or \code{"iqr"}
#'   (\code{stats::IQR(x) / 1.349}). When not \code{"none"}, every column of
#'   \code{X} is divided by the chosen scale before optimization to improve
#'   Hessian conditioning. Coefficients and standard errors are back-transformed
#'   to the user's natural units via the delta method, so reported quantities
#'   are invariant to this choice.
#' @param weights Optional vector of weights for each choice situation. If \code{NULL}, equal weights are used. All weights must be finite and strictly positive.
#' @param weights_col Optional name of a column in \code{data} holding per-row
#'   weights (convenience workflow only). The column must be constant within each
#'   \code{id_col} (one weight per choice situation) and is collapsed accordingly.
#'   Mutually exclusive with \code{weights}. All weights must be finite and strictly
#'   positive. Used for choice-based / WESML
#'   weighting; pair with \code{se_method = "sandwich"} for valid inference.
#' @param outside_opt_label Label for the outside option (if any). If \code{NULL}, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param use_asc Logical indicating whether to include alternative-specific constants (ASCs) in the model.
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in the
#'   returned object for \code{predict()} and post-estimation functions.
#' @param se_method Method for computing standard errors: \code{"hessian"}
#'   (default, analytical Hessian), \code{"bhhh"} (outer product of gradients),
#'   or \code{"sandwich"} (robust Huber--White / WESML variance
#'   \eqn{A^{-1} B A^{-1}}). Use \code{"sandwich"} under choice-based / WESML
#'   weighting.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control} instead.
#' @returns A \code{choicer_mnl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()}, \code{vcov()},
#'   \code{logLik()}, \code{AIC()}, \code{BIC()}, \code{nobs()}, \code{predict()}.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3; beta_true <- c(1.0, -0.5)
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, V := drop(as.matrix(.SD) %*% beta_true), .SDcols = c("x1","x2")]
#' dt[, prob := exp(V) / sum(exp(V)), by = id]
#' dt[, choice := as.integer(alt == sample(alt, 1, prob = prob)), by = id]
#'
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' summary(fit)
#' coef(fit)
#' AIC(fit)
#' predict(fit, type = "shares")
#' }
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
    weights_col = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    use_asc = TRUE,
    keep_data = TRUE,
    scale_vars = c("none", "sd", "mad", "iqr"),
    se_method = c("hessian", "bhhh", "sandwich"),
    nloptr_opts = NULL
) {
  cl <- match.call()

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  scale_vars <- match.arg(scale_vars)
  se_method <- match.arg(se_method)

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  cs_meta <- if (has_data) attr(data, "choice_sampling") else attr(input_data, "choice_sampling")

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_input && !is.null(weights_col)) {
    stop("`weights_col` is only supported in the convenience (data) workflow. ",
         "Bake weights into `input_data` via prepare_mnl_data(weights_col = ) ",
         "or supply `weights` to prepare_mnl_data().")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "and covariate_cols.")
    }
    # WESML provenance present but no weights supplied: auto-adopt the recorded
    # weight column, or error -- never silently fit unweighted under a WESML label.
    if (!is.null(cs_meta) && is.null(weights) && is.null(weights_col)) {
      wn <- cs_meta$weight_name
      if (!is.null(wn) && wn %in% names(data)) {
        weights_col <- wn
        message("Detected WESML choice-based-sampling provenance; applying attached ",
                "weights from column '", wn, "'.")
      } else {
        stop("Data carries WESML choice-based-sampling provenance but no weights were ",
             "supplied, and the recorded weight column (",
             if (is.null(wn)) "unknown" else paste0("'", wn, "'"),
             ") is not present in `data`. Pass `weights_col=` or `weights=` explicitly.",
             call. = FALSE)
      }
    }
    input_list <- prepare_mnl_data(
      data, id_col, alt_col, choice_col, covariate_cols,
      weights = weights,
      weights_col = weights_col,
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

  # Parameter names and index map (built early so the scaling layer can address
  # blocks and so theta_hat/vcov/se can be named downstream).
  asc_names <- if (use_asc) {
    paste0("ASC_", input_list$alt_mapping[2:J][[alt_col]])
  } else {
    character(0)
  }
  param_names <- c(colnames(input_list$X), asc_names)
  param_map <- list(beta = seq_len(K_x))
  if (n_asc > 0) param_map$asc <- K_x + seq_len(n_asc)

  # --- Variable scaling (optional) --------------------------------------------
  # Scale columns of X by their sample SD (or robust SD-equivalent) to improve
  # Hessian conditioning. Keep the natural-scale matrix for storage; theta_hat
  # and vcov are back-transformed after optimization. The back-transform is
  # purely multiplicative (1/sX on the beta block, identity on ASCs).
  natural_X <- input_list$X
  sX <- rep(1, K_x); names(sX) <- colnames(input_list$X)
  bt_mult <- rep(1, n_params)
  bt_shift <- rep(0, n_params)
  if (scale_vars != "none" && K_x > 0) {
    sX <- .column_scales(input_list$X, scale_vars)
    .assert_scales_ok(sX, scale_vars, "fixed-coefficient")
    input_list$X <- sweep(input_list$X, 2, sX, "/")
    bt_mult[param_map$beta] <- 1 / sX
  }

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

  theta_hat <- opt$par
  names(theta_hat) <- param_names

  # Choice-based-sampling provenance and a guardrail for weighted inference.
  weights_nonuniform <- length(unique(input_list$weights)) > 1
  if (weights_nonuniform && se_method == "bhhh") {
    warning("Non-uniform weights detected with se_method = 'bhhh': BHHH/OPG ",
            "standard errors use the w^1 meat (sum w_i s_i s_i')^{-1}, which is ",
            "NOT a valid choice-based-sampling (WESML) correction; the correct ",
            "sandwich meat is w^2. Use se_method = 'sandwich' for valid WESML ",
            "inference.",
            call. = FALSE)
  } else if (weights_nonuniform && se_method != "sandwich") {
    warning("Non-uniform weights detected. If these are sampling/WESML ",
            "weights, use se_method = 'sandwich' for valid inference.",
            call. = FALSE)
  }
  choice_sampling <- if (!is.null(cs_meta)) {
    utils::modifyList(as.list(cs_meta),
                      list(se_method = se_method, weights_applied = weights_nonuniform))
  } else if (weights_nonuniform) {
    list(scheme = "user", se_method = se_method, weights_applied = TRUE)
  } else {
    NULL
  }
  if (!is.null(cs_meta) && !weights_nonuniform) {
    if (has_input) {
      stop("`input_data` is flagged as a WESML choice-based sample (it carries ",
           "`choice_sampling` provenance), but the resolved weights are uniform. ",
           "Fitting would produce an invalid unweighted estimator mislabeled as ",
           "WESML. To proceed, either bake the non-uniform WESML weights into ",
           "`input_data` via prepare_mnl_data(weights = ) / prepare_mnl_data(weights_col = ), ",
           "or, if you deliberately want an unweighted fit, strip the provenance with ",
           "`attr(input_data, \"choice_sampling\") <- NULL`.",
           call. = FALSE)
    }
    warning("WESML provenance is present but the applied weights are uniform; the fit ",
            "is effectively unweighted and is NOT a WESML-corrected estimator.",
            call. = FALSE)
  }

  # Compute vcov eagerly using the selected SE method. For "sandwich"
  # (robust / WESML) errors, form V = A^{-1} B A^{-1} with bread A = weighted
  # negated Hessian and meat B = weight-squared OPG (pass weights^2 to the
  # weight-free BHHH routine). Computed in scaled space; back-transform applies.
  if (se_method == "sandwich") {
    A_bread <- mnl_loglik_hessian_parallel(
      theta = theta_hat, X = input_list$X, alt_idx = input_list$alt_idx,
      choice_idx = input_list$choice_idx, M = input_list$M,
      weights = input_list$weights, use_asc = use_asc,
      include_outside_option = input_list$include_outside_option
    )
    B_meat <- mnl_bhhh_parallel(
      theta = theta_hat, X = input_list$X, alt_idx = input_list$alt_idx,
      choice_idx = input_list$choice_idx, M = input_list$M,
      weights = input_list$weights^2, use_asc = use_asc,
      include_outside_option = input_list$include_outside_option
    )
    vcov_result <- .sandwich_combine(A_bread, B_meat)
  } else {
    hess <- switch(
      se_method,
      hessian = mnl_loglik_hessian_parallel(
        theta = theta_hat, X = input_list$X, alt_idx = input_list$alt_idx,
        choice_idx = input_list$choice_idx, M = input_list$M,
        weights = input_list$weights, use_asc = use_asc,
        include_outside_option = input_list$include_outside_option
      ),
      bhhh = mnl_bhhh_parallel(
        theta = theta_hat, X = input_list$X, alt_idx = input_list$alt_idx,
        choice_idx = input_list$choice_idx, M = input_list$M,
        weights = input_list$weights, use_asc = use_asc,
        include_outside_option = input_list$include_outside_option
      )
    )
    vcov_result <- invert_hessian(hess)
  }
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # --- Back-transform to natural scale ----------------------------------------
  if (scale_vars != "none") {
    bt <- .backtransform_estimates(theta_hat, vcov_result, bt_mult, bt_shift, param_names)
    theta_hat <- bt$theta
    vcov_result <- bt$vcov_result
    input_list$X <- natural_X
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
    },
    scale_vars = scale_vars,
    sX = sX,
    se_method = se_method,
    choice_sampling = choice_sampling
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
#' @param weights Optional vector of weights for each choice situation. If `NULL`, equal weights are used. All weights must be finite and strictly positive.
#' @param weights_col Optional name of a column in `data` holding per-row
#'   weights. The column must be constant within each `id_col` (one weight per
#'   choice situation) and is collapsed accordingly. Mutually exclusive with
#'   `weights`. All weights must be finite and strictly positive.
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
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
#' str(input$X)
#' input$alt_mapping
#' @export
prepare_mnl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    weights_col = NULL
) {
  ## Preliminary housekeeping --------------------------------------------------
  # Capture any choice-based-sampling provenance before column drops / coercion,
  # so it can be carried onto the returned object for the advanced pathway.
  cs_provenance <- attr(data, "choice_sampling")
  dt <- data.table::as.data.table(data)[]

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols)
  if (!is.null(weights) && !is.null(weights_col)) {
    stop("Supply only one of `weights` or `weights_col`.")
  }
  if (!is.null(weights_col)) needed <- c(needed, weights_col)
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
  data.table::setorderv(dt, c(id_col, "alt_int"))

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

  ## Collapse a row-level weight column to one weight per choice situation.
  ## Done AFTER ordering/filtering so alignment is by id, never by position.
  if (!is.null(weights_col)) {
    if (!is.numeric(dt[[weights_col]])) {
      stop("`", weights_col, "` must be numeric.")
    }
    nuniq <- dt[, data.table::uniqueN(get(weights_col)), by = id_col][["V1"]]
    if (any(nuniq != 1L)) {
      stop("`", weights_col, "` must be constant within each '", id_col,
           "' (one weight per choice situation).")
    }
    wmap <- dt[, get(weights_col)[1L], by = id_col]
    weights <- wmap[["V1"]][match(ids, wmap[[id_col]])]
    if (any(!is.finite(weights))) {
      stop("`", weights_col, "` produced non-finite weights.")
    }
  }

  ## choice_idx[i] - 1-based index *within* the choice set data
  ## 0 == outside option (only if chosen = 0 for all inside options & include_outside_option == TRUE)
  if (include_outside_option) {
    # start with all-zero (everyone assumed to pick the outside good)
    choice_idx <- integer(N)
    chosen_dt <- dt[get(choice_col) == 1, .(pos = idx_in_group), by = id_col]

    # match chosen ids back to the master index vector
    data.table::setkeyv(chosen_dt, id_col)
    choice_idx[match(chosen_dt[[id_col]], ids)] <- chosen_dt$pos
  } else {
    # exactly one explicit choice per id
    choice_idx <- dt[get(choice_col) == 1, idx_in_group]
  }

  if (is.null(weights)) weights <- rep(1, length(M))

  ## Weights must be finite and strictly positive. Zero/negative weights would
  ## silently invalidate weighted and WESML sandwich inference (w in the bread,
  ## w^2 in the meat). Validated here so every resolution path (weights=,
  ## weights_col=, and the uniform default) is covered.
  if (any(!is.finite(weights))) {
    stop("Weights must be finite, but non-finite values (NA/NaN/Inf) were found.",
         call. = FALSE)
  }
  if (any(weights <= 0)) {
    stop("Weights must be strictly positive, but values <= 0 were found.",
         call. = FALSE)
  }

  ## Alternatives summary ------------------------------------------------------

  if (include_outside_option) {
    inside_alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
    outside_alt_mapping <- data.table::data.table(alt_int=0L, N_OBS = N, N_CHOICES = sum(choice_idx == 0L))
    outside_alt_mapping[[alt_col]] <- outside_opt_label
    alt_mapping <- list(outside_alt_mapping, inside_alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping, c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
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
  out <- structure(
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
  if (!is.null(cs_provenance)) {
    attr(out, "choice_sampling") <- cs_provenance
  }
  out
}

