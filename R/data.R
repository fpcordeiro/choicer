#' Intercity travel mode choice
#'
#' Stated choices of intercity travel mode for 210 travellers, each choosing
#' among the same four modes: air, train, bus and car. This is the classic
#' Greene & Hensher (1997) data set, reshaped into choicer's long layout (one
#' row per traveller-by-alternative). It is a convenient, recognizable example
#' for multinomial and nested logit models and the demand/welfare toolkit
#' (elasticities, diversion ratios, willingness-to-pay, counterfactuals).
#'
#' @format A data frame with 840 rows (210 travellers x 4 modes) and 9 columns:
#' \describe{
#'   \item{id}{Integer traveller (choice situation) identifier, 1-210.}
#'   \item{mode}{Factor giving the travel mode: \code{"air"}, \code{"train"},
#'     \code{"bus"} or \code{"car"}. Use as the alternative column.}
#'   \item{choice}{Integer indicator, 1 for the chosen mode and 0 otherwise.
#'     Exactly one mode is chosen per traveller.}
#'   \item{wait}{Terminal waiting time in minutes (0 for car).}
#'   \item{travel}{In-vehicle travel time in minutes.}
#'   \item{vcost}{In-vehicle cost component, in currency units.}
#'   \item{gcost}{Generalized cost measure, in currency units.}
#'   \item{income}{Household income (traveller level, in thousands).}
#'   \item{size}{Size of the travelling party (traveller level).}
#' }
#'
#' @details
#' \code{wait}, \code{travel}, \code{vcost} and \code{gcost} vary across modes
#' within a traveller, while \code{income} and \code{size} are traveller-level
#' attributes that are constant across modes. A standard specification regresses
#' the choice on \code{wait}, \code{travel} and \code{vcost} with
#' alternative-specific constants; \code{vcost} then plays the role of price for
#' willingness-to-pay and consumer-surplus calculations.
#'
#' @source
#' Greene, W. H. and Hensher, D. A. (1997). Reshaped from the \code{TravelMode}
#' data distributed with the \pkg{AER} package
#' (\url{https://CRAN.R-project.org/package=AER}). The same data appear in
#' Greene's \emph{Econometric Analysis} and in several other choice-modelling
#' packages.
#'
#' @examples
#' data(mode_choice)
#' head(mode_choice)
#' table(mode_choice$mode[mode_choice$choice == 1L])
"mode_choice"
