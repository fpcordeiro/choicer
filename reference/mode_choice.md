# Intercity travel mode choice

Stated choices of intercity travel mode for 210 travellers, each
choosing among the same four modes: air, train, bus and car. This is the
classic Greene & Hensher (1997) data set, reshaped into choicer's long
layout (one row per traveller-by-alternative). It is a convenient,
recognizable example for multinomial and nested logit models and the
demand/welfare toolkit (elasticities, diversion ratios,
willingness-to-pay, counterfactuals).

## Usage

``` r
mode_choice
```

## Format

A data frame with 840 rows (210 travellers x 4 modes) and 9 columns:

- id:

  Integer traveller (choice situation) identifier, 1-210.

- mode:

  Factor giving the travel mode: `"air"`, `"train"`, `"bus"` or `"car"`.
  Use as the alternative column.

- choice:

  Integer indicator, 1 for the chosen mode and 0 otherwise. Exactly one
  mode is chosen per traveller.

- wait:

  Terminal waiting time in minutes (0 for car).

- travel:

  In-vehicle travel time in minutes.

- vcost:

  In-vehicle cost component, in currency units.

- gcost:

  Generalized cost measure, in currency units.

- income:

  Household income (traveller level, in thousands).

- size:

  Size of the travelling party (traveller level).

## Source

Greene, W. H. and Hensher, D. A. (1997). Reshaped from the `TravelMode`
data distributed with the AER package
(<https://CRAN.R-project.org/package=AER>). The same data appear in
Greene's *Econometric Analysis* and in several other choice-modelling
packages.

## Details

`wait`, `travel`, `vcost` and `gcost` vary across modes within a
traveller, while `income` and `size` are traveller-level attributes that
are constant across modes. A standard specification regresses the choice
on `wait`, `travel` and `vcost` with alternative-specific constants;
`vcost` then plays the role of price for willingness-to-pay and
consumer-surplus calculations.

## Examples

``` r
data(mode_choice)
head(mode_choice)
#>   id  mode choice wait travel vcost gcost income size
#> 1  1   air      0   69    100    59    70     35    1
#> 2  1 train      0   34    372    31    71     35    1
#> 3  1   bus      0   35    417    25    70     35    1
#> 4  1   car      1    0    180    10    30     35    1
#> 5  2   air      0   64     68    58    68     30    2
#> 6  2 train      0   44    354    31    84     30    2
table(mode_choice$mode[mode_choice$choice == 1L])
#> 
#>   air train   bus   car 
#>    58    63    30    59 
```
