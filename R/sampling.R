# Choice-based (endogenous stratified) sampling and WESML weighting
#
# Implements Manski-Lerman (1977) Weighted Exogenous Sample Maximum Likelihood
# (WESML) weights for choice-based samples, plus a helper to draw such samples
# from a population data set. Strata are defined by the *chosen* alternative.

# --- internal helpers --------------------------------------------------------

#' Validate the id / alt / choice columns of a long choice data set
#' @noRd
.validate_choice_columns <- function(dt, id_col, alt_col, choice_col,
                                     include_outside_option) {
  needed <- c(id_col, alt_col, choice_col)
  miss <- setdiff(needed, names(dt))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  if (!all(dt[[choice_col]] %in% c(0, 1))) {
    stop("`", choice_col, "` must contain only 0 and 1.")
  }
  counts <- dt[, sum(get(choice_col)), by = id_col][["V1"]]
  if (!include_outside_option && any(counts != 1)) {
    stop("Each '", id_col, "' must have exactly one chosen alternative ",
         "(one 1 in '", choice_col, "').")
  }
  if (include_outside_option && any(counts > 1)) {
    stop("Each '", id_col, "' must have at most one chosen alternative ",
         "(one 1 in '", choice_col, "').")
  }
  invisible(TRUE)
}

#' Chosen stratum (chosen alternative) for each choice situation
#'
#' Returns a data.table with one row per choice situation: the id column and a
#' character stratum key `.strat` equal to the chosen alternative coerced to
#' character. Situations with no explicit choice are mapped to the outside
#' stratum when `include_outside_option = TRUE`.
#' @noRd
.choicer_chosen_strata <- function(dt, id_col, alt_col, choice_col,
                                   include_outside_option = FALSE,
                                   outside_opt_label = NULL) {
  chosen <- dt[get(choice_col) == 1]
  id_vec    <- chosen[[id_col]]
  strat_vec <- as.character(chosen[[alt_col]])

  all_ids   <- unique(dt[[id_col]])
  no_choice <- setdiff(all_ids, id_vec)
  if (length(no_choice) > 0) {
    if (!include_outside_option) {
      stop("Some choice situations have no chosen alternative.")
    }
    if (is.null(outside_opt_label)) {
      stop("include_outside_option = TRUE but outside_opt_label is NULL; ",
           "provide the label used for the (implicit) outside good.")
    }
    id_vec    <- c(id_vec, no_choice)
    strat_vec <- c(strat_vec,
                   rep(as.character(outside_opt_label), length(no_choice)))
  }

  per_id <- data.table::data.table(.id = id_vec, .strat = strat_vec)
  data.table::setnames(per_id, ".id", id_col)
  data.table::setkeyv(per_id, id_col)
  per_id
}

#' Validate a named vector of shares against the realized strata
#'
#' Coerces names to character, checks finiteness/non-negativity, requires the
#' names to match `realized` exactly, and (optionally) renormalizes to sum 1.
#' @noRd
.check_shares <- function(x, realized, name, normalize = TRUE) {
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("`", name, "` must be a named numeric vector ",
         "(names = alternatives / strata).")
  }
  xn <- as.numeric(x)
  names(xn) <- as.character(names(x))
  if (any(!is.finite(xn)) || any(xn <= 0)) {
    stop("`", name, "` must be finite and strictly positive ",
         "(a realized stratum cannot have a zero share).")
  }
  .match_strata(names(xn), realized, name)
  if (normalize) {
    s <- sum(xn)
    if (!is.finite(s) || s <= 0) {
      stop("`", name, "` must have a positive sum.")
    }
    if (abs(s - 1) > 1e-8) {
      message("`", name, "` did not sum to 1 (sum = ",
              format(s, digits = 6), "); renormalizing.")
      xn <- xn / s
    }
  }
  xn
}

#' Require two character sets to match exactly, with a clear diff on failure
#' @noRd
.match_strata <- function(nm, realized, argname) {
  nm <- as.character(nm)
  if (anyDuplicated(nm)) {
    dups <- unique(nm[duplicated(nm)])
    stop("`", argname, "` has duplicate names: ", paste(dups, collapse = ", "))
  }
  missing_s <- setdiff(realized, nm)
  extra_s   <- setdiff(nm, realized)
  if (length(missing_s) > 0 || length(extra_s) > 0) {
    stop("`", argname, "` names must match the chosen strata exactly.\n",
         if (length(missing_s) > 0) {
           paste0("  Missing strata: ", paste(missing_s, collapse = ", "), "\n")
         },
         if (length(extra_s) > 0) {
           paste0("  Unexpected names: ", paste(extra_s, collapse = ", "), "\n")
         },
         "  Realized strata: ", paste(realized, collapse = ", "))
  }
  invisible(TRUE)
}

# --- exported functions ------------------------------------------------------

#' WESML weights for choice-based (endogenous stratified) samples
#'
#' Computes Manski-Lerman (1977) Weighted Exogenous Sample Maximum Likelihood
#' (WESML) weights for a choice-based sample. The weight for a choice situation
#' whose chosen alternative is \eqn{j} is \eqn{w = Q(j) / H(j)}, where
#' \eqn{Q(j)} is the population share of alternative \eqn{j} and \eqn{H(j)} its
#' sample share among choosers. Using these weights in
#' \code{\link{run_mxlogit}} restores consistency under choice-based sampling;
#' pair them with \code{se_method = "sandwich"} for valid (robust) standard
#' errors (the plain inverse-Hessian is invalid under weighting).
#'
#' Strata are defined by the chosen alternative and keyed by
#' \code{as.character(alt)} so numeric and character alternative codes match
#' supplied share names unambiguously.
#'
#' @param data A long-format choice data set (data.frame or data.table), one
#'   row per alternative per choice situation.
#' @param id_col,alt_col,choice_col Column names identifying the choice
#'   situation, the alternative, and the 0/1 chosen indicator.
#' @param Q Named numeric vector of population shares, one entry per chosen
#'   stratum (names matched to \code{as.character(alt)}), each strictly
#'   positive. Renormalized to sum 1 if needed.
#' @param H Optional named numeric vector of sample shares. If \code{NULL}
#'   (default) it is computed from \code{data} as the fraction of choice
#'   situations choosing each alternative.
#' @param normalize If \code{TRUE} (default) the returned weights are scaled to
#'   mean 1. This does not affect the point estimates or the sandwich variance.
#' @param attach If \code{TRUE}, return \code{data} with a row-level weight
#'   column appended (the per-situation weight repeated across all rows of a
#'   situation), ready to pass to \code{run_mxlogit(weights_col = ...)}. If
#'   \code{FALSE} (default) return an id-keyed table of weights.
#' @param weight_name Name of the weight column (default \code{".wesml_weight"}).
#' @param outside_opt_label,include_outside_option Set
#'   \code{include_outside_option = TRUE} and supply \code{outside_opt_label}
#'   when the outside good is implicit (choice situations with no \code{1} in
#'   \code{choice_col} are treated as having chosen the outside good).
#'
#' @returns Either an id-keyed \code{data.table} with columns \code{id_col} and
#'   \code{weight_name} (default), or, when \code{attach = TRUE}, a copy of
#'   \code{data} with the weight column appended. The result carries \code{"Q"},
#'   \code{"H"}, and \code{"choice_sampling"} attributes recording provenance.
#'
#' @references Manski, C. F. and Lerman, S. R. (1977). The Estimation of Choice
#'   Probabilities from Choice Based Samples. \emph{Econometrica} 45(8),
#'   1977-1988. Train, K. E. (2009). \emph{Discrete Choice Methods with
#'   Simulation}, Section 3.7. Cambridge University Press.
#' @seealso \code{\link{sample_by_choice}}, \code{\link{run_mxlogit}},
#'   \code{\link{wesml_vcov}}
#' @examples
#' library(data.table)
#' set.seed(1)
#' N <- 300L; J <- 3L
#' pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
#' pop[, x1 := rnorm(.N)]
#' pop[, w1 := rnorm(.N)]
#' pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#'
#' # Population shares of the chosen alternative
#' Q <- prop.table(table(pop[choice == 1, alt]))
#' wt <- wesml_weights(pop, "id", "alt", "choice", Q = Q)
#' head(wt)
#' @export
wesml_weights <- function(data, id_col, alt_col, choice_col, Q,
                          H = NULL, normalize = TRUE, attach = FALSE,
                          weight_name = ".wesml_weight",
                          outside_opt_label = NULL,
                          include_outside_option = FALSE) {
  dt <- data.table::as.data.table(data)
  .validate_choice_columns(dt, id_col, alt_col, choice_col,
                           include_outside_option)

  per_id   <- .choicer_chosen_strata(dt, id_col, alt_col, choice_col,
                                     include_outside_option, outside_opt_label)
  strata   <- per_id[[".strat"]]
  N        <- nrow(per_id)
  realized <- sort(unique(strata))

  Qn <- .check_shares(Q, realized, "Q", normalize = TRUE)

  if (is.null(H)) {
    tab <- per_id[, .N, by = ".strat"]
    Hn  <- stats::setNames(tab[["N"]] / N, tab[[".strat"]])
  } else {
    Hn <- .check_shares(H, realized, "H", normalize = FALSE)
  }
  if (any(Hn[realized] <= 0)) {
    bad <- realized[Hn[realized] <= 0]
    stop("Sample share H is zero for stratum/strata: ",
         paste(bad, collapse = ", "), " (no choosers in the sample).")
  }

  w <- as.numeric(Qn[strata] / Hn[strata])
  if (normalize) w <- w / mean(w)

  cs <- list(scheme = "wesml", Q = Qn, H = Hn[realized],
             meat = "robust", source = "wesml_weights")

  wt <- data.table::data.table(.id = per_id[[id_col]], .w = w)
  data.table::setnames(wt, c(".id", ".w"), c(id_col, weight_name))

  if (attach) {
    if (weight_name %in% names(dt)) {
      stop("Column '", weight_name, "' already exists in `data`; ",
           "choose a different `weight_name`.")
    }
    out <- merge(data.table::copy(dt), wt, by = id_col,
                 all.x = TRUE, sort = FALSE)
  } else {
    out <- wt
    data.table::setkeyv(out, id_col)
  }
  data.table::setattr(out, "Q", Qn)
  data.table::setattr(out, "H", Hn[realized])
  data.table::setattr(out, "choice_sampling", cs)
  out[]
}

#' Draw a choice-based sample stratified by the chosen alternative
#'
#' Subsamples whole choice situations from a population data set according to
#' fixed per-stratum quotas, where strata are defined by the \emph{chosen}
#' alternative. The input data set is treated as the population, so the
#' population shares \eqn{Q(j)} are known exactly; the returned sample carries a
#' ready-to-use WESML weight column (see \code{\link{wesml_weights}}).
#'
#' Sampling is by choice situation (id), never by row: all alternative-rows of a
#' sampled situation are kept together. Sampling is \strong{without
#' replacement}.
#'
#' @param data,id_col,alt_col,choice_col As in \code{\link{wesml_weights}}.
#' @param n_per_alt Either a single integer applied to every stratum, or a named
#'   integer vector of per-stratum counts (names matched to
#'   \code{as.character(alt)}, covering all strata). Mutually exclusive with
#'   \code{frac_per_alt}.
#' @param frac_per_alt Either a single fraction in \code{[0, 1]} applied to every
#'   stratum, or a named numeric vector of per-stratum fractions. Mutually
#'   exclusive with \code{n_per_alt}.
#' @param seed Optional integer seed for reproducible sampling.
#' @param weight_name Name of the attached weight column (default
#'   \code{".wesml_weight"}).
#' @param outside_opt_label,include_outside_option As in
#'   \code{\link{wesml_weights}} (for an implicit outside good).
#'
#' @returns A \code{data.table} subsample with the weight column appended and
#'   \code{"Q"}, \code{"H"}, and \code{"choice_sampling"} attributes (the last
#'   records the scheme, shares, quotas, and \code{meat = "robust"}).
#'
#' @references Manski, C. F. and Lerman, S. R. (1977). \emph{Econometrica}
#'   45(8), 1977-1988.
#' @seealso \code{\link{wesml_weights}}, \code{\link{run_mxlogit}}
#' @examples
#' library(data.table)
#' set.seed(1)
#' N <- 600L; J <- 3L
#' pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
#' pop[, x1 := rnorm(.N)]
#' pop[, w1 := rnorm(.N)]
#' pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#'
#' s <- sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 50L, seed = 1L)
#' attr(s, "choice_sampling")$H   # realized sample shares
#' head(s[[".wesml_weight"]])
#' @export
sample_by_choice <- function(data, id_col, alt_col, choice_col,
                             n_per_alt = NULL, frac_per_alt = NULL,
                             seed = NULL, weight_name = ".wesml_weight",
                             outside_opt_label = NULL,
                             include_outside_option = FALSE) {
  if (is.null(n_per_alt) == is.null(frac_per_alt)) {
    stop("Supply exactly one of `n_per_alt` or `frac_per_alt`.")
  }
  dt <- data.table::as.data.table(data)
  if (weight_name %in% names(dt)) {
    stop("Column '", weight_name, "' already exists in `data`; ",
         "choose a different `weight_name`.")
  }
  .validate_choice_columns(dt, id_col, alt_col, choice_col,
                           include_outside_option)

  per_id   <- .choicer_chosen_strata(dt, id_col, alt_col, choice_col,
                                     include_outside_option, outside_opt_label)
  strata   <- per_id[[".strat"]]
  realized <- sort(unique(strata))
  ids      <- per_id[[id_col]]
  avail    <- stats::setNames(
    as.integer(table(factor(strata, levels = realized))), realized)

  target <- .resolve_targets(n_per_alt, frac_per_alt, realized, avail)

  zero_strata <- realized[target[realized] == 0L]
  if (length(zero_strata) > 0) {
    stop("Zero sampling target for population stratum/strata: ",
         paste(zero_strata, collapse = ", "),
         ". Every population stratum must be sampled for the WESML weights to ",
         "target the full population; increase the quota (with frac_per_alt, ",
         "round(frac * available) must be >= 1) or use `n_per_alt`.")
  }

  if (!is.null(seed)) set.seed(seed)

  keep_ids <- list()
  for (s in realized) {
    n_s <- target[[s]]
    if (n_s <= 0L) next
    ids_s <- ids[strata == s]
    if (n_s > length(ids_s)) {
      stop("Stratum '", s, "' has only ", length(ids_s),
           " choice situations but ", n_s,
           " were requested (sampling without replacement).")
    }
    keep_ids[[s]] <- if (n_s == length(ids_s)) ids_s else sample(ids_s, n_s)
  }
  keep_ids <- unlist(keep_ids, use.names = FALSE)
  if (length(keep_ids) == 0L) {
    stop("No choice situations selected; check `n_per_alt`/`frac_per_alt`.")
  }

  sub <- dt[get(id_col) %in% keep_ids]

  # Population shares Q(j) come from the full data (input == population). Every
  # realized stratum is retained (zero targets are rejected above), so Q is used
  # as-is -- never renormalized over a subset, which would silently change the
  # WESML target population.
  sub_per_id <- .choicer_chosen_strata(sub, id_col, alt_col, choice_col,
                                       include_outside_option,
                                       outside_opt_label)
  sub_strata <- sub_per_id[[".strat"]]
  N_sub      <- nrow(sub_per_id)

  tabQ  <- per_id[, .N, by = ".strat"]
  Q_pop <- stats::setNames(tabQ[["N"]] / nrow(per_id), tabQ[[".strat"]])

  tabH  <- sub_per_id[, .N, by = ".strat"]
  H_sub <- stats::setNames(tabH[["N"]] / N_sub, tabH[[".strat"]])

  w <- as.numeric(Q_pop[sub_strata] / H_sub[sub_strata])
  w <- w / mean(w)

  wt <- data.table::data.table(.id = sub_per_id[[id_col]], .w = w)
  data.table::setnames(wt, c(".id", ".w"), c(id_col, weight_name))
  out <- merge(sub, wt, by = id_col, all.x = TRUE, sort = FALSE)

  cs <- list(scheme = "wesml", Q = Q_pop, H = H_sub, meat = "robust",
             source = "sample_by_choice", quotas = target)
  data.table::setattr(out, "Q", Q_pop)
  data.table::setattr(out, "H", H_sub)
  data.table::setattr(out, "choice_sampling", cs)
  out[]
}

#' Resolve per-stratum sampling targets from n_per_alt / frac_per_alt
#' @noRd
.resolve_targets <- function(n_per_alt, frac_per_alt, realized, avail) {
  if (!is.null(n_per_alt)) {
    if (is.null(names(n_per_alt))) {
      if (length(n_per_alt) != 1L) {
        stop("`n_per_alt` must be a single number (applied to every stratum) ",
             "or a named vector covering all strata.")
      }
      if (!is.finite(n_per_alt) || n_per_alt < 0) {
        stop("`n_per_alt` must be a finite, non-negative number.")
      }
      if (abs(n_per_alt - round(n_per_alt)) > 1e-8) {
        stop("`n_per_alt` must be a whole number (a fixed count); ",
             "use `frac_per_alt` for fractional sampling.")
      }
      tgt <- stats::setNames(rep(as.integer(round(n_per_alt)), length(realized)),
                             realized)
    } else {
      .match_strata(names(n_per_alt), realized, "n_per_alt")
      if (any(!is.finite(n_per_alt)) || any(n_per_alt < 0)) {
        stop("`n_per_alt` values must be finite and non-negative.")
      }
      if (any(abs(n_per_alt - round(n_per_alt)) > 1e-8)) {
        stop("`n_per_alt` values must be whole numbers (fixed counts); ",
             "use `frac_per_alt` for fractional sampling.")
      }
      nm  <- as.character(names(n_per_alt))
      tgt <- stats::setNames(as.integer(round(n_per_alt))[match(realized, nm)],
                             realized)
    }
  } else {
    if (is.null(names(frac_per_alt))) {
      if (length(frac_per_alt) != 1L) {
        stop("`frac_per_alt` must be a single number in [0, 1] ",
             "or a named vector covering all strata.")
      }
      fr <- stats::setNames(rep(as.numeric(frac_per_alt), length(realized)),
                            realized)
    } else {
      .match_strata(names(frac_per_alt), realized, "frac_per_alt")
      nm <- as.character(names(frac_per_alt))
      fr <- stats::setNames(as.numeric(frac_per_alt)[match(realized, nm)],
                            realized)
    }
    if (any(!is.finite(fr)) || any(fr < 0 | fr > 1)) {
      stop("`frac_per_alt` values must be finite and in [0, 1].")
    }
    tgt <- stats::setNames(as.integer(round(fr * avail[realized])), realized)
  }
  if (any(tgt < 0L)) stop("Sampling targets must be non-negative.")
  tgt
}
