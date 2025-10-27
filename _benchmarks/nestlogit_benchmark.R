rm(list = ls(all.names = TRUE))
gc()

library(mlogit)
library(data.table)
devtools::load_all()

data("TravelMode", package = "AER")

# Hensher and Greene (2002), table 1 p.8-9 model 5
TravelMode$incomeother <- with(TravelMode, ifelse(mode %in% c('air', 'car'), income, 0))

nl <- mlogit(
  choice ~ gcost + wait + incomeother,
  TravelMode,
  nests = list(public = c('train', 'bus'), other = c('car','air'))
  )

summary(nl) |> print()

# Choicer
dt_travel <- as.data.table(TravelMode)

dt_travel[, choice := fifelse(choice=="yes", 1L, 0L) |> as.integer()]

dt_inputs <- prepare_mnl_data(
  data = dt_travel,
  id_col = "individual",
  alt_col = "mode",
  choice_col = "choice",
  covariate_cols = c("gcost", "wait", "incomeother")
)

alt_mapping <- dt_inputs$alt_mapping |> copy()

alt_mapping[, nest := fcase(
  mode %chin% c("train","bus"), 1L,
  mode %chin% c("car","air"), 2L
)]

dt_inputs$nest_idx <- alt_mapping$nest

nloptr_opts <- list(
  "algorithm" = "NLOPT_LD_LBFGS",
  "xtol_rel" = 1.0e-8,
  "maxeval" = 1e+3,
  "print_level" = 0L,
  "check_derivatives" = TRUE,
  "check_derivatives_print" = 'all'
)

param_names <- c("gcost", "wait", "incomeother","Lambda_public","Lambda_other", paste0("ASC_", alt_mapping[alt_int > 1]$mode))

nestlogit_result <- run_nestlogit(
    input_data = dt_inputs,
    use_asc = TRUE,
    param_names = param_names,
    nloptr_opts = nloptr_opts
)
