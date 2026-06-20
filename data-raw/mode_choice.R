# Prepare the `mode_choice` dataset shipped with choicer.
#
# Source: the classic intercity travel-mode choice data of Greene & Hensher
# (1997), distributed in the AER package as `TravelMode`. 210 travellers each
# face the same four modes (air, train, bus, car) and choose one. We reshape it
# into choicer's long layout (one row per decider x alternative) and give the
# columns self-describing names. AER is used only to source the data here; it is
# NOT a dependency of choicer.
#
# Re-run with: Rscript data-raw/mode_choice.R

stopifnot(requireNamespace("AER", quietly = TRUE))
data("TravelMode", package = "AER")
dat <- TravelMode

mode_choice <- data.frame(
  id     = as.integer(as.character(dat$individual)),
  mode   = factor(as.character(dat$mode),
                  levels = c("air", "train", "bus", "car")),
  choice = as.integer(dat$choice == "yes"),
  wait   = as.numeric(dat$wait),    # terminal waiting time (minutes; 0 for car)
  travel = as.numeric(dat$travel),  # in-vehicle travel time (minutes)
  vcost  = as.numeric(dat$vcost),   # in-vehicle cost component (currency units)
  gcost  = as.numeric(dat$gcost),   # generalized cost measure (currency units)
  income = as.numeric(dat$income),  # household income (decider-level, 1000s)
  size   = as.integer(dat$size),    # travelling party size (decider-level)
  stringsAsFactors = FALSE
)

# Order by decider then by the canonical mode ordering.
mode_choice <- mode_choice[order(mode_choice$id, as.integer(mode_choice$mode)), ]
rownames(mode_choice) <- NULL

# Sanity checks: exactly one chosen mode per decider, four modes each.
stopifnot(
  all(tapply(mode_choice$choice, mode_choice$id, sum) == 1L),
  all(table(mode_choice$id) == 4L)
)

save(mode_choice, file = "data/mode_choice.rda",
     compress = "bzip2", version = 2)
message("Wrote data/mode_choice.rda (", nrow(mode_choice), " rows)")
