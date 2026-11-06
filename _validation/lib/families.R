# ============================================================================
# The library loader, the registry of data-generating processes (DGPs) and
# the families (lib/families/<ID>.R): loading and validating the contracts
# that let several families, runners and reports be written independently.
# Definitions only. A script loads the library with
#   source(file.path(<_validation>, "lib", "families.R"), local = (boot <- new.env()))
#   L <- boot$vs_source_lib(file.path(<_validation>, "lib"))
# ============================================================================

# Sources the library into one environment (parent: the global environment,
# so attached packages are visible) and returns it: common.R, families.R
# (whose registry the DGP files fill), simulate.R, fit.R, metrics.R, sbc.R,
# references.R, then dgp_*.R alphabetically, then any other lib/*.R
# alphabetically; calibrate.R (a script) never, and files that do not exist
# yet are skipped. Sets `.vs_lib` and `.vs_root` (the `_validation/`
# directory) there.
vs_source_lib <- function(lib_dir) {
  lib_dir <- normalizePath(lib_dir, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$.vs_lib <- lib_dir
  env$.vs_root <- dirname(lib_dir)
  files <- setdiff(list.files(lib_dir, pattern = "\\.R$"), "calibrate.R")
  first <- c("common.R", "families.R", "simulate.R", "fit.R", "metrics.R", "sbc.R",
             "references.R")
  dgp <- sort(grep("^dgp_", files, value = TRUE), method = "radix")
  rest <- sort(setdiff(files, c(first, dgp)), method = "radix")
  for (f in c(intersect(first, files), dgp, rest)) {
    sys.source(file.path(lib_dir, f), envir = env, keep.source = FALSE)
  }
  env
}

# --- DGP registry -----------------------------------------------------------

.vs_dgps <- new.env(parent = emptyenv())

# Registers a DGP: list(id, name, theta0, simulate, moments, targets,
# variants). `theta0` is a named numeric, or a list of them named by variant;
# `simulate(N, seed, variant = "default", ...)` returns list(data, truth,
# meta); `moments(sim)` a named numeric; `targets` a data.frame(moment, lo,
# hi); `variants` the variant names, "default" among them.
dgp_register <- function(entry, replace = FALSE) {
  where <- if (is.list(entry) && is.character(entry$id)) sprintf("DGP '%s'", entry$id[1]) else "DGP"
  bad <- function(field, what) stop(where, ": `", field, "` must be ", what, call. = FALSE)
  if (!is.list(entry)) bad("entry", "a list")
  need <- c("id", "name", "theta0", "simulate", "moments", "targets", "variants")
  miss <- setdiff(need, names(entry))
  if (length(miss)) stop(where, ": missing field(s) ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.character(entry$id) || length(entry$id) != 1L || !grepl("^[A-Za-z0-9._-]+$", entry$id)) {
    bad("id", "a single name of letters, digits, '.', '_' or '-'")
  }
  if (!is.character(entry$name) || length(entry$name) != 1L) bad("name", "a single string")
  if (!is.character(entry$variants) || !length(entry$variants) || anyNA(entry$variants) ||
      anyDuplicated(entry$variants) || !"default" %in% entry$variants) {
    bad("variants", "unique variant names including \"default\"")
  }
  th <- if (is.list(entry$theta0)) entry$theta0 else list(default = entry$theta0)
  if (is.null(names(th)) || !all(names(th) %in% entry$variants) || !"default" %in% names(th)) {
    bad("theta0", "a named numeric, or a list of them named by variant (\"default\" included)")
  }
  for (v in th) {
    if (!is.numeric(v) || !length(v) || is.null(names(v)) || any(!nzchar(names(v))) ||
        anyDuplicated(names(v)) || any(!is.finite(v))) {
      bad("theta0", "finite numerics with unique parameter names")
    }
  }
  if (!is.function(entry$simulate) ||
      !all(c("N", "seed", "variant") %in% names(formals(entry$simulate)))) {
    bad("simulate", "function(N, seed, variant = \"default\", ...)")
  }
  if (!is.function(entry$moments)) bad("moments", "function(sim)")
  tg <- entry$targets
  if (!is.data.frame(tg) || !all(c("moment", "lo", "hi") %in% names(tg)) ||
      !is.character(tg$moment) || anyDuplicated(tg$moment) || !is.numeric(tg$lo) ||
      !is.numeric(tg$hi) || any(tg$lo > tg$hi, na.rm = TRUE)) {
    bad("targets", "a data.frame(moment, lo, hi) with unique moments and lo <= hi")
  }
  if (!replace && exists(entry$id, envir = .vs_dgps, inherits = FALSE)) {
    stop(where, " is already registered", call. = FALSE)
  }
  assign(entry$id, entry, envir = .vs_dgps)
  invisible(entry)
}

list_dgps <- function() sort(ls(.vs_dgps), method = "radix")

get_dgp <- function(id) {
  if (!is.character(id) || length(id) != 1L || !exists(id, envir = .vs_dgps, inherits = FALSE)) {
    reg <- list_dgps()
    stop("no DGP '", paste(id, collapse = ","), "' registered; registered: ",
         if (length(reg)) paste(reg, collapse = ", ") else "none", call. = FALSE)
  }
  get(id, envir = .vs_dgps, inherits = FALSE)
}

# A DGP's truth for one variant (theta0 itself when it is not per variant).
dgp_theta0 <- function(dgp, variant = "default") {
  if (!is.list(dgp$theta0)) return(dgp$theta0)
  dgp$theta0[[variant]] %||% stop("DGP '", dgp$id, "' has no theta0 for variant '", variant, "'")
}

# --- Families ---------------------------------------------------------------

# Family ids: the lib/families/*.R files; T* files are plumbing fixtures,
# listed only with fixtures = TRUE.
list_families <- function(lib_dir = get0(".vs_lib", envir = environment(list_families)),
                          fixtures = FALSE) {
  ids <- sub("\\.R$", "", list.files(file.path(lib_dir, "families"), pattern = "\\.R$"))
  if (!fixtures) ids <- ids[!startsWith(ids, "T")]
  sort(ids, method = "radix")
}

# Sources <dir>/<id>.R into a child of `lib_env` (so its functions see the
# library), takes the `family` it assigns and validates it against the
# family-entry contract (validate_family()).
load_family <- function(id, lib_env, dir = file.path(lib_env$.vs_lib, "families")) {
  path <- file.path(dir, paste0(id, ".R"))
  if (!file.exists(path)) {
    stop("no family file ", path, "; families there: ",
         paste(sub("\\.R$", "", list.files(dir, pattern = "\\.R$")), collapse = ", "),
         call. = FALSE)
  }
  env <- new.env(parent = lib_env)
  sys.source(path, envir = env, keep.source = FALSE)
  fam <- get0("family", envir = env, inherits = FALSE)
  if (is.null(fam)) stop(path, " does not assign `family`", call. = FALSE)
  validate_family(fam, id = id, where = path)
}

# Checks a family entry and returns it:
#   list(id, claims, dgp,
#        cells = data.frame(cell, N, ...),             unique cell names
#        methods = list(<name> = function(sim, cell, seed)),
#        methods_by_cell = NULL | list(<cell> = <method names>),
#        R = list(quick = c(<cell> = n), mc = c(<cell> = n)),
#        R_by_method = NULL | list(<cell> = c(<method> = n)),   optional
#        simulate = function(cell, seed) list(data, truth, meta),
#        estimands = function(fit, sim) data.frame,
#        gates = list(cell = <cells, one or more> | list(<method> = <cells>),
#                     metrics = <character> | list(<method> = <character>))
#                | NULL (a findings-only family, with `findings`),
#        gates_quick = NULL | <as gates, cells among those of R$quick,
#                               and optionally alpha = <its Holm level>>,   optional
#        findings = NULL | function(records) <named list of data.frames>,  optional
#        pooled_groups = NULL | function(records) <named vector constant -> group>)  optional
# A method's and simulate()'s `cell` is the cell's row as a list (cell, N
# and the other columns). R_by_method (optional, mc tier) runs some methods
# of a cell at fewer reps than R$mc[cell], on the same paired reps 1..n:
# its cells exist, its methods run in the cell, and n is a whole number in
# 1..R$mc[cell]. The gates' cells are the cells where every method is
# gated, or a list naming each gated method's cells (cells where it runs;
# a method the list does not name is gated nowhere); gates_quick (optional)
# replaces gates in the quick tier, where gate_quick() (lib/metrics.R) reads its
# alpha (default 0.01; the mc tier's level is fixed). findings (optional) turns
# the family's records (every cell and method) into its own tables, which report.R
# renders in place of the generic ladders and QQ plots; a family that has
# it may have no gates. pooled_groups (optional) maps a (cell, method)'s
# records to groups of the constants (e.g. expected-count terciles), so
# that report.R also pools their coverage by group (pooled_tails(groups =)).
# The gates' metrics apply to every method when
# they are a character vector, else to the methods the list names; an entry
# is a metric name, optionally with a route suffix and a block filter
# ("se_ratio:hessian", "bias_std@beta"): only the shape (non-empty strings)
# is checked here, not the vocabulary. Errors name the offending field.
validate_family <- function(fam, id = NULL, where = "family") {
  bad <- function(field, what) {
    stop(where, ": family field `", field, "` must be ", what, call. = FALSE)
  }
  name_ok <- function(x) is.character(x) && !anyNA(x) && all(grepl("^[A-Za-z0-9._+-]+$", x))
  if (!is.list(fam)) bad("family", "a list")
  need <- c("id", "claims", "dgp", "cells", "methods", "methods_by_cell", "R", "simulate",
            "estimands", "gates")
  miss <- setdiff(need, names(fam))
  if (length(miss)) {
    stop(where, ": family is missing field(s) ", paste0("`", miss, "`", collapse = ", "),
         call. = FALSE)
  }
  if (!name_ok(fam$id) || length(fam$id) != 1L) bad("id", "a single name")
  if (!is.null(id) && fam$id != id) bad("id", sprintf("\"%s\", the file's name", id))
  if (!is.character(fam$claims) || !length(fam$claims)) bad("claims", "a character vector")
  if (!is.character(fam$dgp) || length(fam$dgp) != 1L) bad("dgp", "a single DGP id")
  cells <- fam$cells
  if (!is.data.frame(cells) || !nrow(cells)) bad("cells", "a data.frame with a row per cell")
  if (!name_ok(cells$cell)) bad("cells$cell", "a character column of names")
  if (anyDuplicated(cells$cell)) bad("cells$cell", "unique")
  if (!is.numeric(cells$N) || anyNA(cells$N) || any(cells$N < 1) ||
      any(cells$N != floor(cells$N))) {
    bad("cells$N", "a column of positive whole numbers")
  }
  m <- fam$methods
  if (!is.list(m) || !length(m) || !name_ok(names(m)) || anyDuplicated(names(m))) {
    bad("methods", "a list of functions with unique names")
  }
  for (k in names(m)) {
    if (!is.function(m[[k]]) || !identical(names(formals(m[[k]])), c("sim", "cell", "seed"))) {
      bad(paste0("methods$", k), "function(sim, cell, seed)")
    }
  }
  mbc <- fam$methods_by_cell
  if (!is.null(mbc)) {
    if (!is.list(mbc) || is.null(names(mbc)) || !all(names(mbc) %in% cells$cell)) {
      bad("methods_by_cell", "NULL or a list named by cells")
    }
    for (k in names(mbc)) {
      if (!is.character(mbc[[k]]) || !all(mbc[[k]] %in% names(m))) {
        bad(paste0("methods_by_cell$", k), "names of `methods`")
      }
    }
  }
  R <- fam$R
  if (!is.list(R) || !length(R) || is.null(names(R)) || !all(names(R) %in% c("quick", "mc"))) {
    bad("R", "list(quick = c(<cell> = n), mc = c(<cell> = n))")
  }
  for (k in names(R)) {
    r <- R[[k]]
    if (is.null(r)) next
    if (!is.numeric(r) || is.null(names(r)) || anyNA(r) || any(r < 0) || any(r != floor(r))) {
      bad(paste0("R$", k), "whole numbers named by cell")
    }
    if (!all(names(r) %in% cells$cell)) {
      bad(paste0("R$", k), paste0("named by cells (unknown: ",
                                  paste(setdiff(names(r), cells$cell), collapse = ", "), ")"))
    }
  }
  rbm <- fam$R_by_method
  if (!is.null(rbm)) {
    if (!is.list(rbm) || is.null(names(rbm)) || anyDuplicated(names(rbm)) ||
        !all(names(rbm) %in% cells$cell)) {
      bad("R_by_method", "NULL or list(<cell> = c(<method> = n)) named by cells")
    }
    for (k in names(rbm)) {
      r <- rbm[[k]]
      run <- mbc[[k]] %||% names(m)
      if (!is.numeric(r) || !length(r) || is.null(names(r)) || anyDuplicated(names(r)) ||
          !all(names(r) %in% run)) {
        bad(paste0("R_by_method$", k), paste0("numbers named by the cell's methods (",
                                              paste(run, collapse = ", "), ")"))
      }
      cap <- unname(R$mc[k])
      if (!length(cap) || is.na(cap) || anyNA(r) || any(r < 1) || any(r != floor(r)) ||
          any(r > cap)) {
        bad(paste0("R_by_method$", k),
            sprintf("whole numbers in 1..R$mc[\"%s\"] (%s)", k,
                    if (length(cap) && !is.na(cap)) cap else "not set"))
      }
    }
  }
  if (!is.function(fam$simulate) ||
      !identical(names(formals(fam$simulate)), c("cell", "seed"))) {
    bad("simulate", "function(cell, seed)")
  }
  if (!is.function(fam$estimands) ||
      !identical(names(formals(fam$estimands)), c("fit", "sim"))) {
    bad("estimands", "function(fit, sim)")
  }
  # gates, and gates_quick (optional, the quick tier's, its cells among those
  # of R$quick): cell is cells where every method is gated, or a list naming
  # each gated method's cells, cells where the method runs; alpha, if given, a
  # level strictly between 0 and 1.
  check_gates <- function(g, field, allowed) {
    if (!is.list(g) || !all(c("cell", "metrics") %in% names(g))) {
      bad(field, paste("list(cell = <cells> | list(<method> = <cells>),",
                       "metrics = <character> | list(<method> = <character>))"))
    }
    cells_ok <- function(x, ok) {
      is.character(x) && length(x) && !anyNA(x) && !anyDuplicated(x) && all(x %in% ok)
    }
    gc <- g$cell
    if (is.list(gc)) {
      if (!length(gc) || is.null(names(gc)) || anyDuplicated(names(gc)) ||
          !all(names(gc) %in% names(m))) {
        bad(paste0(field, "$cell"), paste0("cells, or a list of them named by methods (methods: ",
                                           paste(names(m), collapse = ", "), ")"))
      }
      for (k in names(gc)) {
        runs <- Filter(function(cl) k %in% (mbc[[cl]] %||% names(m)), allowed)
        if (!cells_ok(gc[[k]], runs)) {
          bad(paste0(field, "$cell$", k), paste0("one or more distinct cells where method ", k,
                                                 " runs (", paste(runs, collapse = ", "), ")"))
        }
      }
    } else if (!cells_ok(gc, allowed)) {
      bad(paste0(field, "$cell"), paste0("one or more distinct cells (cells: ",
                                         paste(allowed, collapse = ", "), ")"))
    }
    if (!is.null(g$alpha) && !(is.numeric(g$alpha) && length(g$alpha) == 1L && !is.na(g$alpha) &&
                                 g$alpha > 0 && g$alpha < 1)) {
      bad(paste0(field, "$alpha"), "a single number between 0 and 1")
    }
    entries_ok <- function(x) is.character(x) && !anyNA(x) && all(nzchar(trimws(x)))
    gm <- g$metrics
    if (is.list(gm)) {
      if (length(gm) && (is.null(names(gm)) || anyDuplicated(names(gm)) ||
                         !all(names(gm) %in% names(m)))) {
        bad(paste0(field, "$metrics"), paste0("a list named by methods (methods: ",
                                              paste(names(m), collapse = ", "), ")"))
      }
      for (k in names(gm)) {
        if (!entries_ok(gm[[k]])) bad(paste0(field, "$metrics$", k), "non-empty strings")
      }
    } else if (!entries_ok(gm)) {
      bad(paste0(field, "$metrics"), "non-empty strings, or a list of them named by methods")
    }
  }
  if (!is.null(fam$findings) &&
      (!is.function(fam$findings) || !identical(names(formals(fam$findings)), "records"))) {
    bad("findings", "NULL or function(records)")
  }
  if (!is.null(fam$pooled_groups) &&
      (!is.function(fam$pooled_groups) || length(formals(fam$pooled_groups)) != 1L)) {
    bad("pooled_groups", "NULL or a function of one argument (the records)")
  }
  if (is.null(fam$gates)) {
    if (is.null(fam$findings)) {
      bad("gates", "a list(cell, metrics), or NULL in a findings-only family (one with `findings`)")
    }
  } else {
    check_gates(fam$gates, "gates", cells$cell)
  }
  if (!is.null(fam$gates_quick)) {
    check_gates(fam$gates_quick, "gates_quick", names(R$quick)[R$quick > 0])
  }
  fam
}

# A cell's row as a list (what simulate() and the methods receive).
family_cell <- function(fam, cell) {
  i <- match(cell, fam$cells$cell)
  if (is.na(i)) {
    stop("family ", fam$id, " has no cell '", cell, "'; cells: ",
         paste(fam$cells$cell, collapse = ", "), call. = FALSE)
  }
  as.list(fam$cells[i, , drop = FALSE])
}

# The methods run in a cell: methods_by_cell's list, else every method.
cell_methods <- function(fam, cell) {
  fam$methods_by_cell[[cell]] %||% names(fam$methods)
}
