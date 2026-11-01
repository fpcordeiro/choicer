# ASC display selection. Named row indices keep alternative labels independent
# of parameter names (which callers may customize).
.fit_asc_rows <- function(object) {
  if (inherits(object, "choicer_hb")) {
    return(stats::setNames(seq_len(nrow(object$delta)),
                           as.character(object$delta$alternative)))
  }
  asc <- object$param_map$asc %||% integer()
  if (!length(asc)) return(asc)
  stats::setNames(asc, as.character(object$alt_mapping[[2L]][-1L]))
}

.select_asc_rows <- function(n, asc, print_asc) {
  invalid <- function() {
    stop("`print_asc` must be TRUE, FALSE, a non-negative integer count, ",
         "a character vector of alternative labels, or an integer vector ",
         "of ASC row positions in the full table.", call. = FALSE)
  }
  if (is.object(print_asc) || !is.null(dim(print_asc)) || anyNA(print_asc)) {
    invalid()
  }
  if (is.logical(print_asc)) {
    if (length(print_asc) != 1L) invalid()
    selected <- if (print_asc) asc else integer()
  } else if (is.character(print_asc)) {
    unknown <- setdiff(print_asc, names(asc))
    if (length(unknown)) {
      stop("`print_asc` contains labels without an ASC in this table: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    selected <- asc[names(asc) %in% print_asc]
  } else if (is.numeric(print_asc) && !is.complex(print_asc)) {
    if (any(!is.finite(print_asc) | print_asc < 0 |
            print_asc != trunc(print_asc))) invalid()
    if (length(print_asc) == 1L) {
      selected <- asc[seq_len(min(print_asc, length(asc)))]
    } else {
      if (any(!print_asc %in% asc)) {
        stop("`print_asc` row positions must identify ASC rows in the full table.",
             call. = FALSE)
      }
      selected <- asc[asc %in% print_asc]
    }
  } else {
    invalid()
  }
  rows <- seq_len(n)
  list(rows = rows[!rows %in% asc | rows %in% selected],
       shown = length(selected), hidden = length(asc) - length(selected))
}

.print_asc_note <- function(selection) {
  if (selection$hidden > 0L) {
    cat(sprintf("ASCs: %d shown, %d hidden (print output truncated).\n",
                selection$shown, selection$hidden))
  }
}

# Recovery tables may omit truth blocks or be subset/reordered by the caller.
# Match labels by parameter name rather than retaining stale row positions.
.recovery_asc_rows <- function(x) {
  asc <- which(x$group %in% c("asc", "delta"))
  labels <- attr(x, "asc_labels", exact = TRUE)
  stats::setNames(asc, unname(labels[x$parameter[asc]]))
}

.recovery_asc_print <- function(out, object, print_asc) {
  asc <- .fit_asc_rows(object)
  parameters <- if (inherits(object, "choicer_hb")) {
    colnames(object$draws$delta)
  } else {
    names(object$coefficients)[asc]
  }
  data.table::setattr(out, "asc_labels",
                     stats::setNames(names(asc) %||% character(), parameters))
  .select_asc_rows(nrow(out), .recovery_asc_rows(out), print_asc)
  data.table::setattr(out, "print_asc", print_asc)
  out
}
